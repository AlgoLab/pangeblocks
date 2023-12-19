#!/usr/bin/env python3
"""
Compute a Variation Graph from an MSA
"""
import cProfile
import time
import json
import argparse
from blocks import Block
from ilp.variaton_graph_parser import asGFA
from pathlib import Path
from dataclasses import astuple
import logging


def main(args):
    # Load set of decomposed blocks
    path_msa = args.path_msa
    path_gfa = args.path_gfa
    path_vert_blocks = args.path_vert_blocks
    dir_subsols = args.dir_subsolutions

    logging.info(f"Computing Pangeblocks for: {Path(path_msa).stem}")
    logging.info("Loading blocks")
    with open(path_vert_blocks) as fp:
        opt_coverage = [Block(*b[:3]) for b in json.load(fp)]

    name_msa = Path(path_msa).stem
    for path_subsol in Path(dir_subsols).rglob(f"*.json"):
        with open(path_subsol) as fp:
            opt_coverage.extend([Block(*b[:3]) for b in json.load(fp)])

    # parse optimal coverage as GFA
    ti = time.time()
    logging.info("Parsing graph as GFA")
    parser = asGFA()
    parser(opt_coverage, path_gfa, path_msa)
    tf = time.time()
    t_gfa = tf - ti
    print(f"time GFA: {t_gfa:0.2}")


if __name__ == "__main__":
    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-msa", help="path to MSA in .fa format")
    parser.add_argument("--path-gfa", help="path to save the output in GFA format")
    parser.add_argument(
        "--dir-subsolutions",
        help="directory where suboptimal coverage are saved",
        dest="dir_subsolutions",
    )
    parser.add_argument(
        "--path-vert-blocks",
        help="json file with vertical blocks",
        dest="path_vert_blocks",
        type=str,
    )
    parser.add_argument(
        "--log-level",
        default="ERROR",
        dest="log_level",
        help="set log level (ERROR/WARNING/INFO/DEBUG)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="[Solve subMSA] %(asctime)s.%(msecs)03d | %(message)s",
        datefmt="%Y-%m-%d@%H:%M:%S",
    )

    logging.basicConfig(
        level=args.log_level,
    )

    main(args)
