import numpy as np
from pathlib import Path
from typing import Optional, Union
from .block import Block
from .analyzer import BlockAnalyzer
from .block_decomposition import block_decomposition

class Decomposer(BlockAnalyzer):

    def __init__(self, return_positional_strings: bool=False):
        self.return_positional_strings=return_positional_strings

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> dict:
        
        if path_blocks is not None:
            list_blocks = self._load_list_blocks(path_blocks)

        # compute matrix with intersections
        inter_blocks=self._matrix_inter_blocks(list_blocks)

        # decompose blocks
        decomposed_blocks=self.decomposition_from_inter_blocks(list_blocks, inter_blocks)
        
        if self.return_positional_strings is True:
            return [block.to_positional_string() for block in decomposed_blocks]        
        
        return decomposed_blocks

    @staticmethod
    def decomposition_from_inter_blocks(list_blocks, inter_blocks):
        "Return a new list of blocks that arise from the decomposition of the intersections"
        decomposed_blocks=[]

        rows, cols = np.where(inter_blocks>0)
        for row, col in zip(rows, cols):
            block1 = list_blocks[row]
            block2 = list_blocks[col]
            # get new blocks (not maximal)
            blocks_from_inter = [block for block in block_decomposition(block1, block2) if block not in [block1,block2]]
            decomposed_blocks.extend(
                blocks_from_inter       
            )
        
        # join decomposed blocks and input blocks
        decomposed_blocks.extend(list_blocks)

        # return non-duplicated blocks
        return list(set(decomposed_blocks))