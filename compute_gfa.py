"""
Compute a Variation Graph from an MSA
"""
import cProfile
import time
import json
import argparse
from src.blocks import Block
from src.ilp.input import InputBlockSet
from src.ilp.optimization import Optimization
from src.ilp.variaton_graph_parser import asGFA

from pathlib import Path

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("--path_blocks", help="json file with blocks covering the MSA")
parser.add_argument("--path_msa", help="path to MSA in .fa format")
parser.add_argument("--path_gfa", help="path to save the output in GFA format")
args = parser.parse_args()

def main():
    # Load set of decomposed blocks
    path_blocks = args.path_blocks #f"../experiment/block_decomposition/{NAME_MSA}.json"
    path_msa= args.path_msa#f"../msas/{NAME_MSA}.fa"
    path_gfa = args.path_gfa

    with open(path_blocks) as fp:
        blocks = [Block(*block) for block in json.load(fp)] 

    ti = time.time()
    inputset_gen = InputBlockSet()
    inputset = inputset_gen(path_msa, blocks)
    tf = time.time()
    t_input = tf - ti
    print(f"time inputblockset: {t_input:0.2}")

    ti = time.time()
    opt = Optimization(blocks=inputset, path_msa=path_msa)
    opt_coverage, times = opt(return_times=True)
    tf = time.time()
    t_ilp = tf - ti
    print(f"time optimization: {t_ilp:0.2}")

    ti = time.time()
    parser=asGFA()
    parser(opt_coverage,path_gfa)
    tf = time.time()
    t_gfa = tf - ti
    print(f"time GFA: {t_gfa:0.2}")


    path_time = Path(path_gfa).parent
    path_time.mkdir(exist_ok=True, parents=True)
    name_msa = Path(path_msa).stem
    # import pdb; pdb.set_trace()
    # times = [round(t,3) for t in times]
    with open(path_time.joinpath(name_msa + ".txt"), "w") as fp:
        fp.write(f"time_input: {t_input}\n")
        fp.write(f"time_ilp: {t_ilp}\n")
        for t_name, t in times.items():
            fp.write(f"time_ilp - {t_name}: {t}\n")
        fp.write(f"time_gfa: {t_gfa}\n")

if __name__=="__main__":
    # import cProfile, pstats
    # profiler = cProfile.Profile()
    # profiler.enable()
    main()
    # profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats('ncalls')
    # stats.print_stats()