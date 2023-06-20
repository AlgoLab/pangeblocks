#!/usr/bin/env python3
import argparse
import json
import time
from Bio import AlignIO
from pathlib import Path
import numpy as np
from collections import defaultdict, namedtuple
from typing import Tuple, List, Union, Optional
import logging

def timer(func):
    "returns output and execution time of 'func'"
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        result = func(*args, **kwargs)
        t2 = time.time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result, round(t2-t1, 4)
    return wrap_func

def load_submsa(filename, start_column=0, end_column=-1):
    "Return MSA from start_column to end_column (both included)"
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

@timer
def compute_vertical_blocks(filename: Union[str,Path], output: Optional[Union[str,Path]] = None, 
                           start_column: int = 0, end_column: int = -1, threshold_vertical_blocks: int = 1):
    "Compute maximal blocks in a submsa"

    logging.info(f"threshold vertical blocks: {threshold_vertical_blocks}")
    # load subMSA
    msa=load_submsa(filename, start_column, end_column)
    n_cols=msa.get_alignment_length()
    n_seqs=len(msa)

    vertical_blocks = []
    
    chars_block = []
    start_col = 0
    for col in range(n_cols):
        chars_col = msa[:,col]
        chars_col = list(set(chars_col))

        if len(chars_col) == 1: 
            chars_block.append(chars_col[0])
        else:
            
            if len(chars_block)>= args.threshold_vertical_blocks:
                end_col = col-1 # end column is included
                vertical_blocks.append(
                    (list(range(n_seqs)), start_col, end_col, "".join(chars_block))
                )
        
            start_col = col+1 # update starting column for the next iteration
            # update initial values for a new block
            chars_block = []
            

    # Save maximal blocks
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w") as fp:
            json.dump(vertical_blocks, fp,)

    return vertical_blocks

if __name__=="__main__":
    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="MSA filename")
    parser.add_argument("-o","--output", help="json file to save vertical blocks", dest="output")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")
    # parser.add_argument("-vb","--only-vertical-blocks", help="Output only vertical blocks: those using all sequences", type=bool, default=False, dest="only_vertical")
    parser.add_argument("--threshold-vertical-blocks", help="vertical blocks with length at least the threshold will be considered to split the MSA", type=int, dest="threshold_vertical_blocks", default=1)
    parser.add_argument("--log-level", default='ERROR', help="set log level (ERROR/WARNING/INFO/DEBUG)", dest="log_level")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                    format='[Maximal Blocks] %(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
    logging.info(f"filename MSA: '{args.filename}'") 
    logging.info(f"subMSA columns [start, end] = {[args.start_column, args.end_column]}")

    vertical_blocks, times = compute_vertical_blocks(
        filename=args.filename, output=args.output, 
        start_column=args.start_column, end_column=args.end_column, 
        threshold_vertical_blocks=args.threshold_vertical_blocks
        )
    
    