import numpy as np
from pathlib import Path
from typing import Optional, Union
from .block import Block
from .analyzer import BlockAnalyzer
from .block_decomposition import block_decomposition
from ..positional_strings.utils import positional_string_from_block

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
            return [positional_string_from_block(block) for block in decomposed_blocks]        
        
        return decomposed_blocks

    @staticmethod
    def decomposition_from_inter_blocks(list_blocks, inter_blocks):
        # TODO: decidir si retornar todos los bloques, o solo aquellos distintos del input
        # ademas, decidir como retornar aquellos bloques que no intersectan, si hacerlo aqui o en una funcion aparte
        # (yo creo que aqui, aquellos con inter_blocks==0)
        decomposed_blocks=[]

        rows, cols = np.where(inter_blocks>0)
        for row, col in zip(rows, cols):
            block1 = list_blocks[row]
            block2 = list_blocks[col]
            
            decomposed_blocks.extend(
                block_decomposition(block1, block2)    
            )    
        return decomposed_blocks