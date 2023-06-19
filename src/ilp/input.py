import numpy as np
from collections import defaultdict
from typing import Union
from pathlib import Path
from Bio import AlignIO
import sys
PATH=Path(__file__).parent.parent.parent
sys.path.append(str(PATH))

from src.blocks import Block 
from src.blocks.block_decomposer import Decomposer

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='[Input Set] %(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class InputBlockSet:
    """Given a set of maximal blocks and the subMSA they come from,
    generate the Input set of block to be used in the ILP
    """    
    def __call__(self, path_msa: Union[str,Path], maximal_blocks: list[Block],
                 start_column: int, end_column: int) -> list[Block]:
        
        # parse maximal blocks as Block objects
        maximal_blocks = [Block(*b) for b in maximal_blocks]
        # load MSA
        msa = self.load_submsa(path_msa, start_column, end_column)
        n_seqs = len(msa) 
        n_cols = msa.get_alignment_length()
        logging.info("subMSA loaded, nrows:(%s), ncols:(%s)" % (n_seqs, n_cols))
        
        # blocks not covered by maximal blocks
        coverage_panel = self.get_coverage_panel(n_seqs, n_cols, maximal_blocks, start_column)
        missing_blocks = self.get_missing_blocks(coverage_panel, msa, start_column)
        
        # block decomposition
        decomposer = Decomposer(row_maximal_decomposition=True)
        decomposed_blocks = decomposer(list_blocks=maximal_blocks)

        # glue missing blocks of one character in the same column
        missing_blocks = self.glue_vertical_blocks(missing_blocks)
        
        # # check missing blocks are correctly labeled
        # for block in missing_blocks: 
        #     # logging.debug("missing block %s", block)
        #     for r in block.K:
        #         r=int(r)
        #         start=int(block.start)
        #         end=int(block.end)
        #         label_msa = msa[r].seq[start:end+1].upper()

        #         for pos_block, c in enumerate(range(start, end+1)):
        #             char_msa = str(msa[r,start:end+1].seq)[pos_block].upper()
        #             char_block = block.label[pos_block].upper() 
        #             if char_msa != char_block:
        #                 logging.debug("incorrect missing block covering position (%s,%s): %s" % (r,c,block.str()))


        blocks_one_char = self.get_blocks_one_char(msa, start_column)
        # check blocks one char are correctly labeled
        for block in blocks_one_char:
            # logging.debug("block one char %s", block)
            for r in block.K:
                r=int(r)
                start=int(block.start)
                end=int(block.end)
                for pos_block, c in enumerate(range(start, end+1)):
                    char_msa = msa[r].seq[c-start_column].upper()
                    char_block = block.label[pos_block].upper() 
                    if char_msa != char_block:
                        logging.debug("incorrect block one char covering position (%s,%s): %s" % (r,c,block.str()))

        # Input set: input blocks:decomposition of blocks  of one position in the MSA)
        input_set_ilp = decomposed_blocks + missing_blocks + blocks_one_char

        return input_set_ilp
    
    def glue_vertical_blocks(self,list_blocks,):
        "Glue blocks of length 1 that shares column"
        new_blocks = []
        blocks_by_start = defaultdict(list)
        for block in list_blocks:
            if len(block.label) == 1: # only one character blocks
                blocks_by_start[(block.start,block.label)].extend(list(block.K))
            else: 
                new_blocks.append(block)
        
        for i_label,K in blocks_by_start.items():
            i, label = i_label
            new_blocks.append(Block(K,i,i,label))

        return new_blocks

    @staticmethod
    def get_coverage_panel(n_seqs, n_cols, blocks, start_column): #FIXME: memory consuming?
        """returns a matrix of size equal to msa (n_seq x n_cols) with 
        the number of blocks in the list_blocks that covers each position"""

        coverage_panel = np.zeros((n_seqs, n_cols))
        for block in blocks:
            for r in block.K:
                for c in range(block.start-start_column,block.end+1-start_column):
                    coverage_panel[r,c] += 1
        return coverage_panel

    def get_missing_blocks(self, coverage_panel, msa, start_column): 
        """return the missing blocks to cover the MSA
        all consecutives one character not covered positions are
        clustered in one block
        """
        rows, cols=np.where(coverage_panel == 0)
        missing_blocks = [(r,c+start_column) for r,c in zip(rows,cols)]

        missing_blocks = sorted(missing_blocks, key= lambda d: (d[0],d[1]))
        idx_missing_blocks_by_seq = defaultdict(list)
        for row,col in missing_blocks: 
            idx_missing_blocks_by_seq[row].append(col)

        # now split each row into separate sublist of (row,col)
        missing_blocks=[]
        for row, cols in idx_missing_blocks_by_seq.items():
            consecutive_pos = self.get_list_consecutive_pos(cols)
            for pos in consecutive_pos: 
                start_submsa = pos[0] - start_column
                end_submsa = pos[-1] - start_column
                label = str(msa[int(row)].seq[start_submsa:end_submsa+1])
                missing_blocks.append(
                    Block(K=(row,), start=pos[0],end=pos[-1], label=label)
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

    def get_blocks_one_char(self, msa, start_column):
        """
        generate trivial blocks, one seq and one col
        'start_column' is used to generate the blocks w.r.t. original MSA
        """
        blocks_one_char = []
        n_cols=msa.get_alignment_length()
        n_seqs=len(msa)

        for col in range(n_cols):
            seq_by_char = defaultdict(list)
            for row in range(n_seqs):
                seq_by_char[msa[row,col]].append(row)

            for c, K in seq_by_char.items():
                # ommit vertical blocks, they will be part of a maximal one
                if len(K) < n_seqs:
                    blocks_one_char.append(
                            Block(K=K, start=col+start_column, end=col+start_column, label=c)
                    )

        return blocks_one_char

    @staticmethod
    def load_submsa(filename, start_column=0, end_column=-1):
        "Return 0-indexed MSA from start_column to end_column (both included)"
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