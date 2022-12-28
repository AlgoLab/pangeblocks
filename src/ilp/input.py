import numpy as np
from collections import defaultdict
from src.blocks import Block
from typing import Union
from pathlib import Path
from Bio import AlignIO
from itertools import cycle
from PIL import Image

COLORS = [
    (255,0,0), # red
    (0,255,0), # green 
    (0,0,255), # blue
    (255,255,0), # yellow
    (0,255,255), # cyan
    (255,0,255), # magenta
    (255,125,0), # orange
    (125,0,255), # blue-magenta
    (255,0,125), # red-magenta
]

class InputBlockSet:

    def __call__(self, path_msa: Union[str,Path], blocks: list[Block]) -> list[Block]:
        
        msa, n_seqs, n_cols = self.load_msa(path_msa)
        coverage_panel = self.get_coverage_panel(n_seqs, n_cols, blocks)
        missing_blocks = self.get_missing_blocks(coverage_panel, msa)
        blocks_one_char = self.get_blocks_one_char(msa, n_seqs, n_cols)
        # set B: input blocks (maximal blocks, the decompositions under intersection by pairs and blocks of one position in the MSA)
        set_B = blocks + missing_blocks# + blocks_one_char  #[block for block in missing_blocks if block.j-block.i+1 > 1]

        return set_B

    @staticmethod
    def get_coverage_panel(n_seqs, n_cols, blocks):
        """returns a matrix of size equal to msa (n_seq x n_cols) with 
        the number of blocks in the list_blocks that covers each position"""

        # coverage_by_pos = defaultdict(int)
        coverage_panel = np.zeros((n_seqs, n_cols))
        for block in blocks:
            for r in block.K:
                for c in range(block.i,block.j+1):
                    coverage_panel[r,c] += 1
        return coverage_panel

    def get_missing_blocks(self, coverage_panel, msa): 
        """return the missing blocks to cover the MSA
        all consecutives one character not covered positions are
        clustered in one block
        """
        rows, cols=np.where(coverage_panel == 0)
        missing_blocks = [(r,c) for r,c in zip(rows,cols)]
        # missing_blocks = get_missing_blocks(coverage_panel)
        missing_blocks = sorted(missing_blocks, key= lambda d: (d[0],d[1]))
        idx_missing_blocks_by_seq = defaultdict(list)
        for pos in missing_blocks: 
            idx_missing_blocks_by_seq[pos[0]].append(pos[1])

        # now split each row into separate sublist of (row,col)
        missing_blocks=[]
        for seq, cols in idx_missing_blocks_by_seq.items():
            consecutive_pos = self.get_list_consecutive_pos(cols)
            for pos in consecutive_pos: 
                label = str(msa[int(seq)].seq[pos[0]:pos[-1]+1])
                missing_blocks.append(
                    Block(K=(seq,), i=pos[0],j=pos[-1], label=label)
                )

        return missing_blocks

    @staticmethod
    def get_list_consecutive_pos(positions: list[int]):
        "given a list with positions, split it in sublists of consecutive numbers"
        sublists = []

        # Set up current list with first element of input
        curr = [positions[0]]

        # For each remaining element:
        for x in positions[1:]:
            # If the next element is not 1 greater than the last seen element
            if x - 1 != curr[-1]:
                # Append the list to the return variable and start a new list
                sublists.append(curr)
                curr = [x]
            # Otherwise, append the element to the current list.
            else:
                curr.append(x)
        sublists.append(curr)
        return sublists

    def get_blocks_one_char(self, msa, n_seqs, n_cols):
        "generate trivial blocks, one seq and one col"
        blocks_one_char = []
        for col in range(n_cols):
            seq_by_char = defaultdict(list)
            for row in range(n_seqs):
                seq_by_char[msa[row,col]].append(row)

            for c, K in seq_by_char.items():
                blocks_one_char.append(
                        Block(K=K, i=col, j=col, label=c)
                )

        return blocks_one_char

    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols