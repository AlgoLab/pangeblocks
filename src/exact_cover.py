#!/usr/bin/env python3
"""
Compute a Variation Graph from an MSA
"""
import sys
import gc
import os
import cProfile
import time
import json
import argparse
from Bio import AlignIO
from blocks import Block
from ilp.input import InputBlockSet
from ilp.optimization import Optimization
from maximal_blocks import compute_maximal_blocks as maximal_blocks_suffixtree # FIXME: did not touch this
from blocks.maximal_blocks.wild_pbwt import compute_maximal_blocks as maximal_blocks_pbwt
from pathlib import Path
from dataclasses import astuple
import logging

from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from collections import namedtuple, defaultdict
from typing import Optional 

logging.basicConfig(level=logging.INFO,
                    format='[Solve SubMSA] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

def label_from_block(block, msa):
    K, start, end = block
    start, end = int(start), int(end)
    return str(msa[int(K[0])].seq[start:end+1])

def solve_submsa(path_msa, start_column, end_column, 
                 solve_ilp, path_save_ilp, path_opt_solution, 
                 obj_function, penalization, min_len, min_coverage, 
                 time_limit, threads_ilp=8,
                 use_wildpbwt: bool = True, bin_wildpbwt: Optional[str] = None, **kwargs ):
    logging.info(f"Working on: {Path(path_msa).stem} | columns [{start_column},{end_column}]")
    
    # # Load set of decomposed blocks
    kwargs_opt = dict(
        obj_function=obj_function,
        penalization=penalization,
        min_len=min_len,
        min_coverage=min_coverage,
        time_limit=time_limit,
        threads_ilp=threads_ilp
    )

    # if only one column is involved, return the blocks with one character
    if start_column==end_column:
        logging.info("start and end columns are the same, skipping ILP")
        msa = load_submsa(path_msa, start_column, end_column)

        blocks_one_char = []
        n_cols=msa.get_alignment_length()
        n_seqs=len(msa)

        for col in range(n_cols):
            seq_by_char = defaultdict(list)
            for row in range(n_seqs):
                seq_by_char[msa[row,col]].append(row)

            for c, K in seq_by_char.items():
                # ommit vertical blocks, they will be part of a maximal one
                if len(K) < n_seqs:
                    blocks_one_char.append(
                            # Block(K=K, start=col+start_column, end=col+start_column, label=c) # FIXME: Block
                            Block(K=K, start=col+start_column, end=col+start_column)
                    )
        opt_coverage = blocks_one_char

    else:
        # 1. compute maximal blocks
        logging.info(f"Computing maximal blocks")
        # Return positions w.r.t. the full MSA
        if use_wildpbwt:
            maximal_blocks = maximal_blocks_pbwt(
                        filename=path_msa,
                        start_column=start_column, end_column=end_column,
                        alphabet_to_ascii = {"-":0,"A":1,"C":2,"G":3,"T":4,"N":5},
                        bin_wildpbwt = bin_wildpbwt,
                        label_blocks = True
                    )
        else:
            maximal_blocks = maximal_blocks_suffixtree(
                        filename=path_msa, 
                        start_column=start_column, end_column=end_column, only_vertical=False
                        )
        
        # 2. compute input set of blocks for the ILP (decomposition is included)
        logging.info(f"Generating input set ({start_column},{end_column})")
        logging.info(f"Number of maximal blocks {len(maximal_blocks)} ({start_column},{end_column})")
        logging.info(f"Size [bytes] of maximal blocks {sys.getsizeof(maximal_blocks)} ({start_column},{end_column})")
        inputset_gen = InputBlockSet()
        inputset = inputset_gen(path_msa, maximal_blocks, start_column, end_column)
        logging.info(f"Generated input set ({start_column},{end_column})")
        del maximal_blocks
        
        # 3. solve the ILP / output ILP model
        # find optimal coverage of the MSA by blocks
        logging.info("Starting optimization")
        opt = Optimization(blocks=inputset, path_msa=path_msa, start_column=start_column, end_column=end_column,
                        log_level=args.log_level, path_save_ilp=path_save_ilp, **kwargs_opt)
        opt_coverage = opt(solve_ilp=solve_ilp)

    for b in opt_coverage:
        logging.info(f"Optimal Coverage block {b.str()}")

    if path_opt_solution:
        full_msa = load_submsa(path_msa)
        Path(path_opt_solution).parent.mkdir(exist_ok=True, parents=True)
        with open(path_opt_solution, "w") as fp:
            blocks = [astuple(block) for block in opt_coverage]
            # blocks = [[ [int(s) for s in b[0]],int(b[1]), int(b[2]),b[3]] for b in blocks]
            blocks = [[ [int(s) for s in b[0]],int(b[1]), int(b[2]), label_from_block(b, full_msa)] for b in blocks] 
            json.dump(blocks, fp)

    del blocks, opt_coverage
    gc.collect()

def load_submsa(filename, start_column=0, end_column=-1):
    "Return 0-indexed MSA from start_column to end_column (both included)"
    # load MSA
    msa = AlignIO.read(filename, "fasta")
    
    # filter sub-MSA if start/end columns are given
    if start_column>0 or end_column!=-1:
        # get last column
        n_cols = msa.get_alignment_length()
        assert start_column < n_cols and end_column < n_cols, f"start_column={start_column}, end_column={end_column}. Must be < {n_cols} (number of columns in the MSA)"
        if end_column == -1:
            msa = msa[:, start_column:] # end_column included
        else:
            msa = msa[:, start_column:end_column+1] # end_column included
    return msa

if __name__=="__main__":
    ## Command line options
    parser = argparse.ArgumentParser()
    # Inputs
    parser.add_argument("--path-msa", help="path to MSA in .fa format", dest="path_msa")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")
    parser.add_argument("-pbwt","--use-wildpbwt", help="Compute maximal blocks with WildPBWT, otherwise use Suffix Tree. Default True", default=True, type=bool, dest="use_wildpbwt")
    parser.add_argument("--bin-wildpbwt", help="path to bin/wild-pbwt", dest="bin_wildpbwt", default="Wild-pBWT/bin/wild-pbwt")
    # ILP
    parser.add_argument("--obj-function", help="objective function", dest="obj_function", choices=["nodes","strings","weighted","depth"])
    parser.add_argument("--penalization", help="penalization for shorter blocks when using 'weighted', and under-covered blocks when using 'depth' as obj_function", dest="penalization", type=int)
    parser.add_argument("--min-len", help="minimum length of shorter blocks when using 'weighted' as obj_function to be penalized", dest="min_len", type=int)
    parser.add_argument("--min-coverage", help="minimum percentage of sequence when using 'depth' as obj_function to be penalized", dest="min_coverage", type=float)
    parser.add_argument("--time-limit", help="time limit in minutes to run the ILP, after this the best solution so far will be returned", dest="time_limit", type=int, default=180)
    parser.add_argument("--threads-ilp", help="threads used by gurobi to solve an ILPi", dest="threads_ilp", type=int, default=4)
    # outputs
    parser.add_argument("--prefix-output", help="prefix to save optimization results. Parent folder will be created if it does not exists.", dest="prefix_output")
    parser.add_argument("--solve-ilp", help="decide wether to solve the ILP or just generate the ILP model", type=bool, default=True, dest="solve_ilp")
    # parser.add_argument("--path-save-ilp", help="path to save the ILP formulation", default=None, dest="path_save_ilp")
    # parser.add_argument("--path-opt-solution", help="file to save optimal solution (Blocks)", dest="path_opt_solution")
    parser.add_argument("--log-level", help='set log level ERROR/WARNING/INFO/DEBUG', default='ERROR', dest='log_level')
    
    parser.add_argument("--submsa-index", help="file with start-end positions of vertical blocks in the MSA", dest="submsa_index")
    parser.add_argument("--workers", help="Workers for ThreadPoolExecutor to solve subMSAs", dest="workers", type=int, default=16)
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                    format='[solve_submsa] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

    OptArgs=namedtuple("OptArgs",["obj_function", "penalization", "min_len", "min_coverage", "time_limit"])
    ArgsPool=namedtuple("Args",["start_column", "end_column", "path_save_ilp", "path_opt_solution"])

    opt_args=OptArgs(args.obj_function, args.penalization, args.min_len, args.min_coverage, args.time_limit)
    submsa = partial(solve_submsa, path_msa=args.path_msa, solve_ilp=args.solve_ilp, 
                         obj_function=args.obj_function, penalization=args.penalization,
                         min_len=args.min_len, min_coverage=args.min_coverage, time_limit=args.time_limit,
                         use_wildpbwt=args.use_wildpbwt, bin_wildpbwt=args.bin_wildpbwt,
                         threads_ilp=args.threads_ilp
                         )
    
    # If index with start-end pairs in between vertical blocks is given, run in parallel all subMSAs 
    if args.submsa_index:

        submsa_index = []
        with open(args.submsa_index, "r") as fp:
            for line in fp.readlines():
                start, end = line.replace("\n","").split("\t")
                start, end = int(start),int(end)
                
                if args.start_column <= start and end <= args.end_column:
                    submsa_index.append((start,end))
        
        if args.workers > 1:
            logging.info("Running with ThreadPoolExecutor: {args.workers} workers")
            def run(argspool):
                "Function to run with ThreadPoolExecutor"
                submsa(start_column=argspool.start_column, end_column=argspool.end_column, path_save_ilp=argspool.path_save_ilp, path_opt_solution=argspool.path_opt_solution)
            
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                with tqdm(total=len(submsa_index), leave=True, ncols=100, bar_format='{l_bar}{bar}| [{elapsed}{postfix}]') as pbar:
                    
                    futures=[]
                    for start,end in submsa_index:
                        path_save_ilp = args.prefix_output + f"_{start}-{end}.mps"
                        path_opt_solution = args.prefix_output + f"_{start}-{end}.json"
                        args_submsa = ArgsPool(start, end, path_save_ilp, path_opt_solution,)
                        future = pool.submit(run, args_submsa)
                        future.add_done_callback(lambda p: pbar.update())
                        future.add_done_callback(lambda p: pbar.set_description(f"Solved subMSA: [{start},{end}]"))
                        futures.append(future)
                    
                    for future in futures:
                        future.result()
        else:
            logging.info("Running with sequentially")
            for start_column, end_column in tqdm(submsa_index, desc="Solving SubMSAs"):
                
                solve_submsa(
                    path_msa=args.path_msa,
                    start_column=start_column,
                    end_column=end_column,
                    solve_ilp=args.solve_ilp,
                    path_save_ilp=args.prefix_output + f"_{start_column}-{end_column}.mps", 
                    path_opt_solution=args.prefix_output + f"_{start_column}-{end_column}.json",
                    obj_function=args.obj_function,
                    penalization=args.penalization,
                    min_len=args.min_len,
                    min_coverage=args.min_coverage,
                    time_limit=args.time_limit,
                    use_wildpbwt=args.use_wildpbwt,
                    bin_wildpbwt=args.bin_wildpbwt,
                    threads_ilp=args.threads_ilp,
                )

    else:

        solve_submsa(
            path_msa=args.path_msa,
            start_column=args.start_column,
            end_column=args.end_column,
            solve_ilp=args.solve_ilp,
            path_save_ilp=args.prefix_output + f"_{args.start_column}-{args.end_column}.mps", 
            path_opt_solution=args.prefix_output + f"_{args.start_column}-{args.end_column}.json",
            obj_function=args.obj_function,
            penalization=args.penalization,
            min_len=args.min_len,
            min_coverage=args.min_coverage,
            time_limit=args.time_limit,
            use_wildpbwt=args.use_wildpbwt,
            bin_wildpbwt=args.bin_wildpbwt,
            threads_ilp=args.threads_ilp,
        )
        # solve_submsa(**vars(args))