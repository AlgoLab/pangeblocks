#!/usr/bin/env python3

import argparse
import json
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
seqs = list(set([str(record.seq) for record in align]))
n_unique_seqs = len(seqs)


def compute_max_blocks(seqs):
    tree = Tree({num: enumerate(seq) for num, seq in enumerate(seqs)})
    blocks = [path for (c, path) in tree.maximal_repeats()]
    decoded_blocks = [
        (b[0][0], 
         b[-1][0], 
        "".join([c[1] for c in b if type(c) == tuple])
        ) for b in blocks
    ]
    return decoded_blocks


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
max_blocks = []
for (start, end) in block_boundaries:
    seqs_sub_msa = list(set([seq[start:end] for seq in seqs]))
    n_cols_sub_msa = end - start
    max_blocks_sub_msa = compute_max_blocks(seqs)
    n_max_blocks_sub_msa = len(max_blocks_sub_msa)

    max_blocks.extend(max_blocks_sub_msa)

n_max_blocks = len(max_blocks)
with open(args.output, "w") as fp:
    json.dump(max_blocks, fp)
