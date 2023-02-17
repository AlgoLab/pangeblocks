import cProfile
from time import time
import json
import argparse
from blocks import Block
from ilp.input import InputBlockSet
from ilp.optimization import Optimization
from ilp.variaton_graph_parser import asGFA
from plot.coverage_msa import CoverageMSA

from pathlib import Path

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("--path_blocks", help="json file with blocks covering the MSA")
parser.add_argument("--path_msa", help="path to MSA in .fa format")
parser.add_argument("--path_gray", help="path save grayscale coverage of the MSA")
parser.add_argument("--path_color", help="path save RGB image with coverage of the MSA by blocks")
args = parser.parse_args()

def main():
    with open(args.path_blocks, "r") as fp:
        blocks=json.load(fp)
    blocks=[Block(*args) for args in blocks]

    cov_msa = CoverageMSA()
    _, n_seqs, n_cols = cov_msa.load_msa(args.path_msa)
        
    if args.path_gray:
        cov_panel = cov_msa.get_coverage_panel(n_seqs, n_cols, blocks)
        cov_msa.save_img(cov_panel, path_save=args.path_gray)

    if args.path_color: 
        cov_color_panel = cov_msa.get_coverage_color_panel(n_seqs, n_cols, blocks)
        cov_msa.save_colored_img(cov_color_panel, path_save=args.path_color)

if __name__=="__main__":
    main()