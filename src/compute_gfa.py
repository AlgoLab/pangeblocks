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
from ilp.variaton_graph_parser import asGFA
from ilp.postprocessing import postprocessing
from pathlib import Path
from dataclasses import astuple
import logging

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("--path_blocks", help="json file with blocks covering the MSA")
parser.add_argument("--path_msa", help="path to MSA in .fa format")
parser.add_argument("--path_gfa", help="path to save the output in GFA format")
parser.add_argument("--path_oc", help="path to save the optimal coverage in json format")
parser.add_argument("--obj_function", help="objective function (nodes/strings/weighted)", dest="obj_function")
parser.add_argument("--penalization", help="penalization for shorter blocks when using 'weighted' as obj_function", dest="penalization", type=int)
parser.add_argument("--min_len", help="minimum length of shorter blocks when using 'weighted' as obj_function to be penalized", dest="min_len", type=int)
parser.add_argument("--time_limit", help="time limit in minutes to run the ILP, after this the best solution so far will be returned", dest="time_limit", type=int)

parser.add_argument(
    "--path_ilp", help="path to save the ILP formulation", default=None)
parser.add_argument("--log_level", default='ERROR', dest='log_level',
                    help='set log level (ERROR/WARNING/INFO/DEBUG')
args = parser.parse_args()
logging.basicConfig(level=args.log_level, )

kwargs_opt = dict(
    obj_function=args.obj_function,
    penalization=args.penalization,
    min_len=args.min_len,
    time_limit=args.time_limit
)

def main():
    # Load set of decomposed blocks
    path_blocks = args.path_blocks
    path_msa= args.path_msa
    path_gfa = args.path_gfa
    path_oc = args.path_oc

    logging.info(f"Computing Pangeblocks for: {Path(path_msa).stem}")
    logging.info("Loading blocks")
    with open(path_blocks) as fp:
        blocks = [Block(*block) for block in json.load(fp)]

    ti = time.time()
    logging.info("Generating input set")
    inputset_gen = InputBlockSet()
    inputset = inputset_gen(path_msa, blocks)
    tf = time.time()
    t_input = tf - ti
    logging.info("Generating input set Done!")
    print(f"time inputblockset: {t_input:0.2}")

    # find optimal coverage of the MSA by blocks
    ti = time.time()
    logging.info("Starting optimization")
    opt = Optimization(blocks=inputset, path_msa=path_msa,
                       log_level=args.log_level, path_save_ilp=args.path_ilp, **kwargs_opt
                       )
    opt_coverage, times = opt(return_times=True)
    tf = time.time()
    t_ilp = tf - ti
    print(f"time optimization: {t_ilp:0.2}")
    
    name_msa = Path(path_msa).stem

    for b in opt_coverage:
        logging.info(f"Optimal Coverage block {b.str()}")

    # save optimal coverage
    if path_oc: 
        Path(path_oc).parent.mkdir(parents=True, exist_ok=True)
        blocks = [astuple(block) for block in opt_coverage]
        # int64 are not json serializable
        blocks = [[ [int(s) for s in b[0]],int(b[1]), int(b[2]),b[3]] for b in blocks] 
        with open(path_oc, "w") as fp:
            json.dump(blocks, fp)

    # parse optimal coverage as GFA
    ti = time.time()
    logging.info("Parsing graph as GFA")
    parser=asGFA()
    parser(opt_coverage, path_gfa, path_msa)
    tf = time.time()
    t_gfa = tf - ti
    print(f"time GFA: {t_gfa:0.2}")

    path_time = Path(path_gfa).parent
    path_time.mkdir(exist_ok=True, parents=True)
    logging.info("saving times")
    with open(path_time.joinpath(name_msa + ".txt"), "w") as fp:
        fp.write(f"time_input\t{t_input}\n")
        fp.write(f"time_ilp\t{t_ilp}\n")
        for t_name, t in times.items():
            fp.write(f"time_ilp-{t_name}\t{t}\n")
        fp.write(f"time_gfa\t{t_gfa}\n")
    logging.info(f"Done with {Path(path_msa).stem}!")

if __name__=="__main__":
    main()