#!/usr/bin/env python3
import argparse
import logging
from Bio import AlignIO

# pangeblocks
from msa.fix_blocks import modify_msa

if __name__=="__main__":
    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-msa","--path-msa", help="MSA to edit", dest="path_msa")
    parser.add_argument("-b","--path-blocks", help="A txt file with maximal blocks to fix in the MSA", dest="path_blocks")
    parser.add_argument("-o","--output", help="json file to save vertical blocks", dest="output")
    parser.add_argument("--log-level", default='ERROR', help="set log level (ERROR/WARNING/INFO/DEBUG)", dest="log_level")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                    format='[modify_msa] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
    logging.info(f"filename MSA: '{args.path_msa}'") 
    msa = modify_msa(msa=args.path_msa, blocks_to_fix=args.path_blocks)
    
    with open(args.output,"w") as fp:
        AlignIO.write(msa, fp, "fasta")