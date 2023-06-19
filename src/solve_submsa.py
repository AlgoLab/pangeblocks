#!/usr/bin/env python3
"""
Compute a Variation Graph from an MSA
"""
import cProfile
import time
import json
import argparse
from blocks import Block
from ilp.input import InputBlockSet
from ilp.optimization import Optimization
# from ilp.variaton_graph_parser import asGFA
# from ilp.postprocessing import postprocessing
from maximal_blocks import compute_maximal_blocks
from pathlib import Path
from dataclasses import astuple
import logging

from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from collections import namedtuple        

# def solve_msa(args):
def solve_submsa(path_msa, start_column, end_column, path_save_ilp, path_opt_solution, solve_ilp, obj_function, penalization, min_len, min_coverage, time_limit, **kwargs):
    logging.info(f"Working on: {Path(path_msa).stem} | columns [{start_column},{end_column}]")

    # # Load set of decomposed blocks
    # path_msa= args.path_msa
    # start_column = args.start_column
    # end_column = args.end_column
    # path_save_ilp = args.path_save_ilp
    # path_opt_solution = args.path_opt_solution
    # solve_ilp = args.solve_ilp
    kwargs_opt = dict(
        obj_function=obj_function,
        penalization=penalization,
        min_len=min_len,
        min_coverage=min_coverage,
        time_limit=time_limit
    )

    # 1. compute maximal blocks
    logging.info(f"Computing maximal blocks")
    maximal_blocks = compute_maximal_blocks(filename=path_msa, start_column=start_column, end_column=end_column, only_vertical=False)
    
    # 2. compute input set of blocks for the ILP (decomposition is included)
    logging.info("Generating input set")
    inputset_gen = InputBlockSet()
    inputset = inputset_gen(path_msa, maximal_blocks, start_column, end_column)

    # 3. solve the ILP / output ILP model
    # find optimal coverage of the MSA by blocks
    logging.info("Starting optimization")
    opt = Optimization(blocks=inputset, path_msa=path_msa, start_column=start_column, end_column=end_column,
                       log_level=args.log_level, path_save_ilp=path_save_ilp, **kwargs_opt)
    opt_coverage = opt(solve_ilp=solve_ilp)

    for b in opt_coverage:
        logging.info(f"Optimal Coverage block {b.str()}")

    if path_opt_solution:
        Path(path_opt_solution).parent.mkdir(exist_ok=True, parents=True)
        with open(path_opt_solution, "w") as fp:
            blocks = [astuple(block) for block in opt_coverage]
            blocks = [[ [int(s) for s in b[0]],int(b[1]), int(b[2]),b[3]] for b in blocks] 
            json.dump(blocks, fp)


if __name__=="__main__":
    ## Command line options
    parser = argparse.ArgumentParser()
    # Inputs
    parser.add_argument("--path-msa", help="path to MSA in .fa format", dest="path_msa")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")
    # ILP
    parser.add_argument("--obj-function", help="objective function (nodes/strings/weighted/depth)", dest="obj_function")
    parser.add_argument("--penalization", help="penalization for shorter blocks when using 'weighted', and under-covered blocks when using 'depth' as obj_function", dest="penalization", type=int)
    parser.add_argument("--min-len", help="minimum length of shorter blocks when using 'weighted' as obj_function to be penalized", dest="min_len", type=int)
    parser.add_argument("--min-coverage", help="minimum percentage of sequence when using 'depth' as obj_function to be penalized", dest="min_coverage", type=float)
    parser.add_argument("--time-limit", help="time limit in minutes to run the ILP, after this the best solution so far will be returned", dest="time_limit", type=int, default=180)
    # outputs
    parser.add_argument("--solve-ilp", help="decide wether to solve the ILP or just generate the ILP model", type=bool, default=True, dest="solve_ilp")
    parser.add_argument("--path-save-ilp", help="path to save the ILP formulation", default=None, dest="path_save_ilp")
    parser.add_argument("--path-opt-solution", help="file to save optimal solution (Blocks)", dest="path_opt_solution")
    parser.add_argument("--log-level", help='set log level ERROR/WARNING/INFO/DEBUG', default='ERROR', dest='log_level')
    
    parser.add_argument("--submsa-index", help="file with start-end positions of vertical blocks in the MSA", dest="submsa_index")
    parser.add_argument("--workers", help="Workers for ThreadPoolExecutor to solve subMSAs", dest="workers", type=int, default=16)
    args = parser.parse_args()


    logging.basicConfig(level=args.log_level,
                    format='[Solve subMSA] %(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

    # If index with start-end pairs in between vertical blocks is given, run in parallel all subMSAs 
    if args.submsa_index:
        OptArgs=namedtuple("OptArgs",["obj_function", "penalization", "min_len", "min_coverage", "time_limit"])
        ArgsPool=namedtuple("Args",["start_column", "end_column", "path_save_ilp", "path_opt_solution"])

        submsa_index = []
        with open(args.submsa_index, "r") as fp:
            for line in fp.readlines():
                start, end = line.replace("\n","").split("\t")
                submsa_index.append((int(start),int(end)))

        opt_args=OptArgs(args.obj_function, args.penalization, args.min_len, args.min_coverage, args.time_limit)
        submsa = partial(solve_submsa, path_msa=args.path_msa, solve_ilp=args.solve_ilp, 
                         obj_function=args.obj_function, penalization=args.penalization,
                         min_len=args.min_len, min_coverage=args.min_coverage, time_limit=args.time_limit
                         )
        def run(argspool):
            "Function to run with ThreadPoolExecutor"
            submsa(start_column=argspool.start_column, end_column=argspool.end_column, path_save_ilp=argspool.path_save_ilp, path_opt_solution=argspool.path_opt_solution)
        
        # for start, end in tqdm(submsa_index, desc="Solving MSA"):
        #     path_save_ilp = args.path_save_ilp + "_{start}-{end}.mps"
        #     path_opt_solution = args.path_opt_solution + "_{start}-{end}.json"
        #     args_submsa = ArgsPool(start, end, path_save_ilp, path_opt_solution)
        #     run(args_submsa)
        
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            with tqdm(total=len(submsa_index)) as progress:
                progress.set_description("Solving subMSA")
                futures=[]
                for start,end in submsa_index:
                    path_save_ilp = args.path_save_ilp + f"_{start}-{end}.mps"
                    path_opt_solution = args.path_opt_solution + f"_{start}-{end}.json"
                    args_submsa = ArgsPool(start, end, path_save_ilp, path_opt_solution)
                    future = pool.submit(run, args_submsa)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)
                
                for future in futures:
                    future.result()
    
    else:
        solve_submsa(**vars(args))

    