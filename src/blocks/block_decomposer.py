import numpy as np
from pathlib import Path
from typing import Optional, Union
from . import Block
from .analyzer import BlockAnalyzer
# from .block_decomposition import block_decomposition
from .decompositions import (
    block_decomposition_row_maximal,
    block_decomposition_not_row_maximal,
    block_decomposition_standard
)

from dataclasses import astuple

import sys # sys.getsizeof()
import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block Decomposer] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class Decomposer(BlockAnalyzer):

    def __init__(self, return_positional_strings: bool=False, standard_decomposition: bool = False):
        self.return_positional_strings=return_positional_strings
        self.standard_decomposition=standard_decomposition
        if standard_decomposition:
            logging.info("Using standard decomposition")
            self.block_decomposition=block_decomposition_standard
        else:
            logging.info("Using row maximal decomposition")
            self.block_decomposition=block_decomposition_row_maximal

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> dict:
        start,end = kwargs.get("start",0), kwargs.get("end",0)

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
        
        # decompose blocks
        logging.info(f"Computing decomposition from intersections of blocks ({start},{end})")
        decomposed_blocks=self.decomposition_from_inter_blocks(list_blocks, inter_blocks) 
        logging.info(f"Computed decomposition from intersections of blocks ({start},{end})")
        logging.info(f"Size [bytes] list of decomposed blocks {sys.getsizeof(list_blocks)} ({start},{end})")
        
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

        return [Block(*b) for b in decomposed_blocks]