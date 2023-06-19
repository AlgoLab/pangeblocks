#!/usr/bin/env python3
import argparse
import json
import time
from Bio import AlignIO
from pathlib import Path
import numpy as np
from collections import defaultdict, namedtuple
from suffix_tree import Tree
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
        # n_seqs = len(align)
        n_cols = msa.get_alignment_length()
        assert start_column < n_cols and end_column < n_cols, f"start_column={start_column}, end_column={end_column}. Must be < {n_cols} (number of columns in the MSA)"
        msa = msa[:, start_column:end_column+1] # end_column included
    # print(msa)
    return msa

@timer
def compute_pos_strings(seqs):
    "positional strings from maximal repeats in a suffix tree"
    tree = Tree({num: enumerate(seq) for num, seq in enumerate(seqs)})
    blocks = [path for (c, path) in tree.maximal_repeats()]
    decoded_blocks = [
        (b[0][0],
         "".join([c[1] for c in b if type(c) == tuple][:len(b)])
         ) for b in blocks
    ]
    
    return decoded_blocks

def set_K_from_pos_string(pos_string, seqs):
    """
    recover set of sequences in the MSA given a positional string (decoded_block)
    seqs: list of sequences in order like in the MSA
    """
    start, label = pos_string
    end = start + len(label) - 1
    K = [row for row, seq in enumerate(seqs) if seq[start:end+1] == label]
    return K, (start, end, label)

@timer
def maximal_blocks_from_pos_strings(pos_strings: list, seqs: list, only_vertical: bool, start_column: int):
    "generate the set of sequences K that each positional string match and return it as a block"
    max_blocks = []
    n_seqs = len(seqs)
    
    # iterate over positional strings and sequences to obtain blocks
    for ps in pos_strings:
        K, extended_ps = set_K_from_pos_string(ps, seqs)
        size_K = len(K)
        
        # if start column is > 0 then the coordinates in the MSA must be updated
        if start_column>0:
            extended_ps = (extended_ps[0] + start_column, extended_ps[1] + start_column, extended_ps[-1])
        
        # check if only vertical blocks are required
        if only_vertical:
            # add block if it is a vertical block
            if size_K == n_seqs:
                max_blocks.append(
                    (K, *extended_ps)
                )
        else:
            max_blocks.append(
                (K, *extended_ps)
            )

    return max_blocks

def compute_maximal_blocks(filename: Union[str,Path], output: Optional[Union[str,Path]] = None, 
                           start_column: int = 0, end_column: int = -1, 
                           only_vertical: bool = False):# multi: bool = True):
    "Compute maximal blocks in a submsa"

    # load subMSA
    msa=load_submsa(filename, start_column, end_column)
    n_cols=msa.get_alignment_length()
    n_seqs=len(msa)

    # # identify unique sequences to create the Tree and compute maximal blocks
    seq_by_string = defaultdict(list)
    all_seqs = []
    for seq, record in enumerate(msa):
        str_seq = str(record.seq)
        seq_by_string[str_seq].append(seq)
        all_seqs.append(str_seq)
    # n_unique_seqs = len(seq_by_string)
    # seqs = seq_by_string.keys()

    ## compute max blocks by first get maximal repeats
    pos_strings, t_pos_strings = compute_pos_strings(all_seqs)

    # to create the blocks from positional strings, we match the string in the positional string
    # to all the rows in the original (without removing duplicates) MSA
    max_blocks, t_max_blocks = maximal_blocks_from_pos_strings(pos_strings=pos_strings, seqs=all_seqs, only_vertical=only_vertical, start_column=start_column)

    # Save maximal blocks
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w") as fp:
            json.dump(max_blocks, fp,)

        # save times
        path_time = Path(output).parent / (Path(output).stem + ".txt")
        with open(path_time, "w") as fp:
            fp.write(f"t_pos_string\t{t_pos_strings}\n")
            fp.write(f"t_max_blocks\t{t_max_blocks}\n")
    
    return max_blocks

if __name__=="__main__":
    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="MSA filename")
    parser.add_argument("-o","--output", help="output file", dest="output")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")
    parser.add_argument("-vb","--only-vertical-blocks", help="Output only vertical blocks: those using all sequences", type=bool, default=False, dest="only_vertical")
    parser.add_argument("--log-level", default='ERROR', help="set log level (ERROR/WARNING/INFO/DEBUG)", dest="log_level")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                    format='[Maximal Blocks] %(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
    logging.info(f"filename MSA: '{args.filename}'") 
    logging.info(f"subMSA columns [start, end] = {[args.start_column, args.end_column]}")
    logging.info(f"Output only vertical blocks = {args.only_vertical}")

    compute_maximal_blocks(
        filename=args.filename, output=args.output, 
        start_column=args.start_column, end_column=args.end_column, 
        only_vertical=args.only_vertical
        )