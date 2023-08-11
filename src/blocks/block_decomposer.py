import numpy as np
from pathlib import Path
from typing import Optional, Union
from . import Block
from .analyzer import BlockAnalyzer
# from .block_decomposition import block_decomposition
from .decompositions import (
    block_decomposition_row_maximal,
    block_decomposition_not_row_maximal,
)

from dataclasses import astuple

import sys # sys.getsizeof()
import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block Decomposer] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class Decomposer(BlockAnalyzer):

    def __init__(self, return_positional_strings: bool=False, row_maximal_decomposition: bool = True):
        self.return_positional_strings=return_positional_strings
        self.row_maximal_decomposition=row_maximal_decomposition

        if row_maximal_decomposition:
            logging.info("Using row maximal decomposition")
            self.block_decomposition=block_decomposition_row_maximal
        else:
            logging.info("Using not row maximal decomposition")
            self.block_decomposition=block_decomposition_not_row_maximal

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> dict:
        
        if path_blocks:
            list_blocks = self._load_list_blocks(path_blocks)
            logging.info(f"Blocks loaded from {str(path_blocks)}")

        # compute list with (idx1,idx2) of intersected blocks
        logging.info(f"Computing pairs of overlapping blocks")
        inter_blocks, list_blocks=self._list_inter_blocks(list_blocks, return_sorted_list=True) 
        logging.info(f"Computed pairs of overlapping blocks")
        logging.info(f"Number of pairs of overlapping blocks {len(inter_blocks)}")
        logging.info(f"Size [bytes] list of pairs indexes of overlapping blocks {sys.getsizeof(inter_blocks)}")
        logging.info(f"Size [bytes] list of blocks {sys.getsizeof(list_blocks)}")
        
        # decompose blocks
        logging.info(f"Computing decomposition from intersections of blocks")
        decomposed_blocks=self.decomposition_from_inter_blocks(list_blocks, inter_blocks) 
        logging.info(f"Computed decomposition from intersections of blocks")
        logging.info(f"Size [bytes] list of decomposed blocks {sys.getsizeof(list_blocks)}")
        
        if self.return_positional_strings is True:
            return [block.to_positional_string() for block in decomposed_blocks]        
        
        return decomposed_blocks

    def decomposition_from_inter_blocks(self, list_blocks, inter_blocks):
        "Return a new list of blocks that arise from the decomposition of the intersections"
        decomposed_blocks=set(astuple(b) for b in list_blocks)

        # decompose pairs of blocks with non-empty intersection
        for pos1, pos2 in inter_blocks:
            block1 = list_blocks[pos1]
            block2 = list_blocks[pos2]
            # get new blocks (not maximal)
            blocks_from_inter = [block for block in self.block_decomposition(block1, block2) if block not in [block1,block2]]
            for b in blocks_from_inter:
                decomposed_blocks.add(astuple(b))
            # decomposed_blocks.extend(
            #     blocks_from_inter       
            # )

        
        # # join decomposed blocks and input blocks
        # decomposed_blocks.extend(list_blocks)

        # decomposed_blocks = set([astuple(b) for b in decomposed_blocks])
        
        # return non-duplicated blocks
        # return list([Block(*b) for b in set(decomposed_blocks)])
        return [Block(*b) for b in decomposed_blocks]