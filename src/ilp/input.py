import numpy as np
from collections import defaultdict
from typing import Union
from pathlib import Path
from Bio import AlignIO

# ------
# FIXME:  better way to import this?
import sys
PATH=Path(__file__).parent.parent.parent
sys.path.append(str(PATH))

from src.blocks import Block # FIXME: Block
from src.blocks.block_decomposer import Decomposer
# ------

import logging
logging.basicConfig(level=logging.INFO,
                    format='[Input] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class InputBlockSet:
    """Given a set of maximal blocks and the subMSA they come from,
    generate the Input set of block to be used in the ILP
    """    

    def __init__(self, standard_decomposition,
                 min_nrows_to_fix_block=0, min_ncols_to_fix_block=0):
        self.standard_decomposition=standard_decomposition 
        logging.info(f">>>> InputBlockSet standard_decomposition={standard_decomposition}")
        
        # criteria to fix some maximal blocks
        self.min_nrows_to_fix_block=min_nrows_to_fix_block
        self.min_ncols_to_fix_block=min_ncols_to_fix_block

    def __call__(self, path_msa: Union[str,Path], maximal_blocks: list[Block],
                 start_column: int, end_column: int) -> list[Block]:
        
        self.start_column = start_column
        self.end_column = end_column

        # parse maximal blocks as Block objects
        maximal_blocks = [Block(*b[:3]) for b in maximal_blocks]
        # load MSA
        self.msa = self.load_submsa(path_msa, start_column, end_column)
        n_seqs = len(self.msa) 
        n_cols = self.msa.get_alignment_length()
        logging.info("subMSA loaded, nrows:(%s), ncols:(%s)" % (n_seqs, n_cols))
        
        # blocks not covered by maximal blocks
        logging.info(f"get coverage panel ({self.start_column},{self.end_column})")
        coverage_panel = self.get_coverage_panel(n_seqs, n_cols, maximal_blocks, start_column)
        logging.info(f"end coverage panel ({self.start_column},{self.end_column})")

        logging.info(f"get missing blocks ({self.start_column},{self.end_column})")
        missing_blocks = self.get_missing_blocks(coverage_panel, self.msa, start_column)
        logging.info(f"end missing blocks ({self.start_column},{self.end_column})")
        
        # block decomposition
        logging.info(f"decomposer ({self.start_column},{self.end_column})")
        decomposer = Decomposer(
                                standard_decomposition=self.standard_decomposition,        #  if False, row maximal decomposition is used
                                min_nrows_to_fix_block=self.min_nrows_to_fix_block,        #  minimum number of rows to fix a maximal block
                                min_ncols_to_fix_block=self.min_ncols_to_fix_block         #  minimum number of columns to fix a maximal block
                                ) 
        decomposed_blocks, fixed_blocks = decomposer(
                                                    list_blocks=maximal_blocks,            # NOTE: decomposed_blocks contains the maximal blocks, check Decomposer.decomposition_from_inter_blocks()
                                                    start=start_column, end=end_column     # to track with logging
                                                    )                

        logging.info(f"end decomposer ({self.start_column},{self.end_column})")
        logging.info(f"Number of decomposed blocks {len(decomposed_blocks)} ({self.start_column},{self.end_column})")        
        logging.info(f"Number of fixed blocks {len(fixed_blocks)} ({self.start_column},{self.end_column})")
        # glue missing blocks of one character in the same column 
        logging.info(f"glue missing blocks ({self.start_column},{self.end_column})")
        missing_blocks = self.glue_vertical_blocks(missing_blocks)
        logging.info(f"Number of missing blocks {len(missing_blocks)} ({self.start_column},{self.end_column})")        

        # missing blocks and fixed blocks will be joined to ommit them from the ILP
        missing_blocks.extend(fixed_blocks)

        # row maximal decomposition
        if not self.standard_decomposition:
            logging.info(f"get blocks one char ({self.start_column},{self.end_column})")
            blocks_one_char = self.get_blocks_one_char(self.msa, start_column, missing_blocks)
            logging.info(f"Number of blocks one char {len(blocks_one_char)} ({self.start_column},{self.end_column})")

            # Input set: input blocks:decomposition of blocks  of one position in the MSA)
            # to avoid having duplicated one_char blocks
            from dataclasses import astuple
            decomposed_blocks=set(astuple(b) for b in decomposed_blocks)
            for block in blocks_one_char:
                decomposed_blocks.add(astuple(block))
            decomposed_blocks = [Block(*b) for b in decomposed_blocks]
            input_set_ilp = decomposed_blocks + blocks_one_char 
            
        else:
            # standard decomposition 
            from dataclasses import astuple
            logging.info(f"# decomposed blocks {len(decomposed_blocks)}")
            decomposed_blocks=set(astuple(b) for b in decomposed_blocks)
            decomposed_blocks = [Block(*b) for b in decomposed_blocks]
            input_set_ilp = decomposed_blocks
        
        return input_set_ilp , missing_blocks
    
    def glue_vertical_blocks(self,list_blocks,):
        "Glue blocks of length 1 that shares column"
        new_blocks = []
        blocks_by_start = defaultdict(list)
        for block in list_blocks:

            if block.len() == 1: # only one character blocks
                # print(block, block.K[0], block.start)
                # print(self.msa)
                label = self.msa[int(block.K[0])].seq[int(block.start) - self.start_column]
                blocks_by_start[(block.start, label)].extend(list(block.K))
            else: 
                new_blocks.append(block)
        
        for i_label,K in blocks_by_start.items():
            i, label = i_label
            # new_blocks.append(Block(K,i,i,label)) # FIXME: Block
            new_blocks.append(Block(K,i,i))
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
                # start_submsa = pos[0] - start_column
                # end_submsa = pos[-1] - start_column
                # # label = str(msa[int(row)].seq[start_submsa:end_submsa+1])
                missing_blocks.append(
                    Block(K=(row,), start=pos[0],end=pos[-1])
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

    def get_blocks_one_char(self, msa, start_column, ommit_blocks):
        """
        generate trivial blocks: one character in each column
        'start_column' is used to generate the blocks w.r.t. original MSA
        'ommit_blocks' is used to ommit positions (r,c) that are covered by any of these blocks
        which are meant to be removed from the ILP
        """
        pos_ommit_blocks = []
        for block in ommit_blocks:
            pos_ommit_blocks.extend(
                (row,col) for row in block.K for col in range(block.start, block.end+1)  
            )

        blocks_one_char = []
        n_cols=msa.get_alignment_length()
        n_seqs=len(msa)

        for col in range(n_cols):
            seq_by_char = defaultdict(list)
            for row in range(n_seqs):
                if (row,col) not in pos_ommit_blocks: 
                    seq_by_char[msa[row,col]].append(row)

            for c, K in seq_by_char.items():
                # ommit vertical blocks, they will be part of a maximal one
                if len(K) < n_seqs:
                    blocks_one_char.append(
                            Block(K=K, start=col+start_column, end=col+start_column) 
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