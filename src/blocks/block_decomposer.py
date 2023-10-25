import json
from pathlib import Path
from typing import Optional, Union
from . import Block
from .analyzer import BlockAnalyzer
# from .block_decomposition import block_decomposition
from .decompositions import (
    block_decomposition_row_maximal,
    block_decomposition_standard
)

from dataclasses import astuple

import sys # sys.getsizeof()
import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block Decomposer] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class Decomposer:
    """
    Given a list of maximal blocks, apply a decomposition by pairs and return a list with all new blocks + the input set of blocks

    Optionally, a list of maximal blocks can be fixed based on a criteria of minimum number of rows and/or cols (a subset of maximal blocks) 
    which are intended to be fixed in the solution of the ILP, so any block resulting from the decomposition that intersects, 
    or is contained in any of these blocks will be removed from the output list. 
    """

    def __init__(self, return_positional_strings: bool=False, standard_decomposition: bool = False,
                 min_nrows_to_fix_block: int = 0, min_ncols_to_fix_block: int = 0
                 ):
        self.return_positional_strings=return_positional_strings
        self.standard_decomposition=standard_decomposition
        self.min_nrows_to_fix_block = min_nrows_to_fix_block
        self.min_ncols_to_fix_block = min_ncols_to_fix_block

        logging.info(f">>>> Decomposer standard_decomposition={standard_decomposition}")
        if standard_decomposition is True:
            logging.info("Using standard decomposition")
            self.block_decomposition=block_decomposition_standard
        else:
            logging.info("Using row maximal decomposition")
            self.block_decomposition=block_decomposition_row_maximal

        if min_nrows_to_fix_block>0 or min_ncols_to_fix_block>0:
            logging.info(f"Fix blocks: blocks with at least {min_nrows_to_fix_block} rows and {min_ncols_to_fix_block} columns will be fixed")
        else:
            logging.info("Fix blocks: no blocks will be fixed")

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> dict:
        start,end = kwargs.get("start",0), kwargs.get("end",0)
        self.start, self.end = start, end # only for logging 
        if path_blocks:
            list_blocks = self._load_list_blocks(path_blocks)
            logging.info(f"Blocks loaded from {str(path_blocks)}")

        # compute list with (idx1,idx2) of intersected blocks
        logging.info(f"Computing pairs of overlapping blocks ({start},{end})")
        inter_blocks, list_blocks=self._list_inter_blocks(list_blocks, return_sorted_list=True) 
        logging.info(f"Computed pairs of overlapping blocks ({start},{end})")
        logging.info(f"Number of pairs of overlapping blocks {len(inter_blocks)} ({start},{end})")
        logging.info(f"Size [bytes] list of pairs indexes of overlapping blocks {sys.getsizeof(inter_blocks)} ({start},{end})")
        logging.info(f"Size [bytes] list of blocks {sys.getsizeof(list_blocks)} ({start},{end})")
        
        # find blocks to fix: indexes of the 'list_blocks'
        if self.min_ncols_to_fix_block>0 or self.min_nrows_to_fix_block>0:
            pos_blocks_to_fix = self.find_blocks_to_fix(list_blocks)
        else: 
            pos_blocks_to_fix = []
        for pos in pos_blocks_to_fix:
            logging.info(f"block fixed {list_blocks[pos]} ({start},{end})")

        # decompose blocks: do not contain fixed blocks (given in 'pos_blocks_to_fix' or any block intersecting one of those ones) 
        logging.info(f"Computing decomposition from intersections of blocks ({start},{end})")
        decomposed_blocks=self.decomposition_from_inter_blocks(list_blocks, inter_blocks, pos_blocks_to_fix) 
        logging.info(f"Computed decomposition from intersections of blocks ({start},{end})")
        logging.info(f"Size [bytes] list of decomposed blocks {sys.getsizeof(list_blocks)} ({start},{end})")
        
        if self.return_positional_strings is True:
            return [block.to_positional_string() for block in decomposed_blocks]        
        
        fixed_blocks = [list_blocks[pos] for pos in pos_blocks_to_fix]
        logging.info(f"intersection between decomposed_blocks and fixed_blocks {len(set(decomposed_blocks).intersection(fixed_blocks))}")
        return decomposed_blocks, fixed_blocks

    def decomposition_from_inter_blocks(self, list_blocks, inter_blocks, pos_blocks_to_fix):
        """Return a new list of blocks that arise from the decomposition of the intersections. This list:
        - DO CONTAIN maximal blocks
        - DO NOT CONTAIN the blocks that are fixed (pos_blocks_to_fix), or any block intersecting one of the fixed blocks
        """
        # decomposed_blocks=set(astuple(block) for pos,block in enumerate(list_blocks) if pos not in pos_blocks_to_fix)
        decomposed_blocks=set()
        pos_discard_block = set() # positions of maximal blocks that intersects with one of the fixed blocks
        # decompose pairs of blocks with non-empty intersection
        for pos1, pos2 in inter_blocks:
            block1 = list_blocks[pos1]
            block2 = list_blocks[pos2]
            list_blocks_decomposition = self.block_decomposition(block1, block2)

            # remove any block intersecting one of the fixed blocks
            if pos1 in pos_blocks_to_fix:
                pos_discard_block.add(pos2)
                fixed_block = block1
                other_block = block2

                for block in list_blocks_decomposition:
                    if self._blocks_intersect(block, fixed_block):
                        list_blocks_decomposition.remove(block)
                if fixed_block in list_blocks_decomposition:
                    list_blocks_decomposition.remove(fixed_block)
                if other_block in list_blocks_decomposition:
                    list_blocks_decomposition.remove(other_block)
                # decomposed_blocks.discard(other_block)

            if pos2 in pos_blocks_to_fix:
                pos_discard_block.add(pos1)
                fixed_block = block2
                other_block = block1

                for block in list_blocks_decomposition:
                    if self._blocks_intersect(block, fixed_block):
                        list_blocks_decomposition.remove(block)
                if fixed_block in list_blocks_decomposition:
                    list_blocks_decomposition.remove(fixed_block)
                if other_block in list_blocks_decomposition:
                    list_blocks_decomposition.remove(other_block)
                # decomposed_blocks.discard(other_block)
            
            for block in list_blocks_decomposition:
                if not any([self._blocks_intersect(block,fix_block) for fix_block in [list_blocks[pos] for pos in pos_blocks_to_fix]]):
                    logging.info(f"> > block from decomposition in input set {block}")
                    decomposed_blocks.add(astuple(block))
        
        # add maximal blocks not intersecting fixed blocks
        pos_discard_block.update(pos_blocks_to_fix)

        for pos,block in enumerate(list_blocks):
            if pos not in pos_discard_block:
                logging.info(f"> > maximal block in input set {block}")
                decomposed_blocks.add(astuple(block))

        return [Block(*b) for b in decomposed_blocks]
    
    def _list_inter_blocks(self, list_blocks: list[Block], return_sorted_list: bool = False) -> list[tuple]:
        # FIXME: remove return_sorted_list param, does not make sense to not return the sorted list, since 'intersections' correspond to indexes in the sorted list
        "list of indexes (in a sorted list by i) of pairs of blocks with non-empty intersection"
        blocks = sorted(list_blocks, key=lambda block: (block.start, len(block.K)))
        
        # save pairs of indexes for the sorted blocks that intersect
        intersections = [] 

        # Lemma1: given two overlapping maximal blocks (K1,b1,e1) and (K2,b2,e2) 
        # [b1,e1] subset [b2,e2] iff K2 subset K1
        # primary intersection: the one that fits the above equivalence, secondary intersection: negation of the equivalence

        n_primary_intersections = 0
        n_secondary_intersections = 0            

        for pos1, block1 in enumerate(blocks[:-1]):
            # compare against the next blocks in the sorted list 
            for rel_pos, block2 in enumerate(blocks[pos1+1:]):
                pos2 = rel_pos + pos1 + 1
                block2 = blocks[pos2]

                # TODO: change this for the method does_block_intersect(block1, block2)
                # check for not empty intersection
                common_rows = list(set(block1.K).intersection(set(block2.K))) # intersection set K
                common_cols = list(set(range(block1.start,block1.end+1)).intersection(set(range(block2.start,block2.end+1)))) # intersection columns [i,j]

                if (common_rows and common_cols):
                    intersections.append((pos1,pos2))

                    # count primary and second intersections
                    if set(block1.K).issubset(block2.K) or set(block2.K).issubset(block1.K):
                        n_primary_intersections +=1                        
                    else:
                        n_secondary_intersections +=1

        logging.info(f"Number of primary intersections {n_primary_intersections} ({self.start},{self.end})")
        logging.info(f"Number of secondary intersections {n_secondary_intersections} ({self.start},{self.end})")

        if return_sorted_list is True:
            return intersections, blocks
        return intersections

    @staticmethod
    def _load_list_blocks(path_list_blocks: Union[str,Path]) -> list[Block]:
        "load blocks in a list from a json file"
        with open(path_list_blocks,"r") as fp:
            list_blocks=[Block(*args) for args in json.load(fp)]
        return list_blocks
    
    def find_blocks_to_fix(self, list_blocks):
        "Return positions in the list_blocks of the blocks that will be fixed"
        # filter (position in the list of) blocks that satisfy the criteria of minimum number of rows and columns
        potential_blocks_to_fix = [pos for pos,block in enumerate(list_blocks) if block.nrows() >= self.min_nrows_to_fix_block and block.ncols() >= self.min_ncols_to_fix_block]
        blocks_to_discard = [] # list of (positions in the list of) blocks that will be discarded
        
        # idx<n>: index in the current list | pos<n>: index in the list of blocks
        for idx1, pos1 in enumerate(potential_blocks_to_fix[:-1]):
            block1 = list_blocks[pos1]
            
            for rel_pos, pos2 in enumerate(potential_blocks_to_fix[idx1+1:]):
                idx2 = rel_pos + pos1 + 1
                pos2 = potential_blocks_to_fix[idx2]
                block2 = list_blocks[pos2]

                # resolve overlapping blocks: the one with the more number of positions (nrows x ncols) is preferred
                if self._blocks_intersect(block1,block2):
                    if block1.ncells() < block2.ncells():
                        blocks_to_discard.append(pos1)
                    else:
                        blocks_to_discard.append(pos2)

        return [pos for pos in potential_blocks_to_fix if pos not in blocks_to_discard]

    @staticmethod
    def _blocks_intersect(block1, block2):
        "Return True if block1 and block2 intersect"
        # check for not empty intersection
        common_rows = list(set(block1.K).intersection(set(block2.K))) # intersection set K
        common_cols = list(set(range(block1.start,block1.end+1)).intersection(set(range(block2.start,block2.end+1)))) # intersection columns [i,j]

        return common_rows and common_cols