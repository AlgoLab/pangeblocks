#!/usr/bin/env python3

import argparse
import json
import time
from Bio import AlignIO
from pathlib import Path
import numpy as np
from collections import defaultdict, namedtuple
from suffix_tree import Tree
from typing import Tuple, List, Union

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="MSA_filename")
parser.add_argument("--output", help="output file")
parser.add_argument("--multi", help="split alignment into regions", action="store_true")
args = parser.parse_args()


align = AlignIO.read(args.filename, "fasta")
n_cols = align.get_alignment_length()
n_seqs = len(align) 
all_seqs = [str(record.seq) for record in align]
seqs = list(set(all_seqs)) # removing duplicates
n_unique_seqs = len(seqs)

def timer(func):
    "returns output and execution time of 'func'"
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        result = func(*args, **kwargs)
        t2 = time.time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result, round(t2-t1,4)
    return wrap_func

@timer
def compute_pos_strings(seqs):
    "positional strings from maximal repeats in a suffix tree"
    tree = Tree({num: enumerate(seq) for num, seq in enumerate(seqs)})
    blocks = [path for (c, path) in tree.maximal_repeats()]
    decoded_blocks = [
        (b[0][0], 
         b[-1][0], 
        "".join([c[1] for c in b if type(c) == tuple])
        ) for b in blocks
    ]
    return decoded_blocks

def set_K_from_pos_string(pos_string, seqs):
    """
    recover set of sequences in the MSA given a positional string (decoded_block)
    seqs: list of sequences in order like in the MSA 
    """
    start, end, label=pos_string
    K = [row for row,seq in enumerate(seqs) if seq[start:end+1]==label]
    return K

@timer
def maximal_blocks_from_pos_strings(pos_strings: list, seqs: list):
    "generate the set of sequences K that each positional string match and return it as a block"
    max_blocks = []
    for ps in pos_strings:
        K=set_K_from_pos_string(ps,seqs)
        max_blocks.append(
            (K,*ps)
            )
    return max_blocks

# block_end is a boolean list that is true iff the corresponding column consists
# of a single character (excluding indels)
# Moreover, block_end has an additional final element that is True
# Therefore each true position of block_end is also the final end of a block
if args.multi:
    block_end = [len(set(col) - set("-")) <= 1 for col in zip(*seqs)]
else:
    block_end = [False] * n_cols
block_end.append(True)

# convert the block_end list into a list of pairs, called block_boundaries.
# Each pair is the [begin, end) extremes of the range of columns characterizing
# a block
ends = [pos + 1 for pos, val in enumerate(block_end) if val]
block_boundaries = zip([pos - 1 for pos in [1] + ends[:-1]], ends)

# compute max blocks and count them
# TODO: run in parallel
pos_strings = []
for (start, end) in block_boundaries:
    seqs_sub_msa = list(set([seq[start:end] for seq in seqs]))
    n_cols_sub_msa = end - start
    pos_strings_sub_msa, t_pos_strings = compute_pos_strings(seqs)    
    pos_strings.extend(pos_strings_sub_msa)

# to create the blocks from positional strings, we match the string in 
# the positional string to all the rows in the original (without removing duplicates)
# MSA 
max_blocks, t_max_blocks = maximal_blocks_from_pos_strings(pos_strings, all_seqs)

# Save maximal blocks
n_max_blocks = len(max_blocks)
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
with open(args.output, "w") as fp:
    json.dump(max_blocks, fp)

# save times
path_time = Path(args.output).parent / (Path(args.output).stem + ".txt")
with open(path_time,"w") as fp:
    fp.write(f"t_pos_string\t{t_pos_strings}\n")
    fp.write(f"t_max_blocks\t{t_max_blocks}\n")