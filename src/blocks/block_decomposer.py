import numpy as np
from pathlib import Path
from typing import Optional, Union
from .block import Block
from .analyzer import BlockAnalyzer
# from .block_decomposition import block_decomposition
from .decompositions import (
    block_decomposition_row_maximal,
    block_decomposition_not_row_maximal,
)

class Decomposer(BlockAnalyzer):

    def __init__(self, return_positional_strings: bool=False, row_maximal_decomposition: bool = True):
        self.return_positional_strings=return_positional_strings
        self.row_maximal_decomposition=row_maximal_decomposition

        if row_maximal_decomposition:
            self.block_decomposition=block_decomposition_row_maximal
        else:
            self.block_decomposition=block_decomposition_not_row_maximal

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> dict:
        
        if path_blocks:
            list_blocks = self._load_list_blocks(path_blocks)

        # compute list with (idx1,idx2) of intersected blocks
        inter_blocks, sorted_blocks =self._list_inter_blocks(list_blocks, return_sorted_list=True) 

        # decompose blocks
        decomposed_blocks=self.decomposition_from_inter_blocks(sorted_blocks, inter_blocks) 
        
        if self.return_positional_strings is True:
            return [block.to_positional_string() for block in decomposed_blocks]        
        
        return decomposed_blocks

    def decomposition_from_inter_blocks(self, list_blocks, inter_blocks):
        "Return a new list of blocks that arise from the decomposition of the intersections"
        decomposed_blocks=[]

        # decompose pairs of blocks with non-empty intersection
        for pos1, pos2 in inter_blocks:
            block1 = list_blocks[pos1]
            block2 = list_blocks[pos2]
            # get new blocks (not maximal)
            blocks_from_inter = [block for block in self.block_decomposition(block1, block2) if block not in [block1,block2]]
            decomposed_blocks.extend(
                blocks_from_inter       
            )
        
        # join decomposed blocks and input blocks
        decomposed_blocks.extend(list_blocks)

        # return non-duplicated blocks
        return list(set(decomposed_blocks))