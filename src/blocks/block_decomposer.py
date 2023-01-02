import numpy as np
from pathlib import Path
from typing import Optional, Union
from .block import Block
from .analyzer import BlockAnalyzer
from .block_decomposition import block_decomposition
from collections import defaultdict
import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')



class Decomposer(BlockAnalyzer):

    def __init__(self, return_positional_strings: bool = False):
        self.return_positional_strings = return_positional_strings

    def __call__(self, list_blocks: Optional[list[Block]] = None, path_blocks: Optional[Union[str, Path]] = None, **kwargs) -> set:

        if path_blocks is not None:
            list_blocks = self._load_list_blocks(path_blocks)

        # decompose blocks

        decomposed_blocks = set(list_blocks)
        for (idx1, idx2) in self._list_inter_blocks(list_blocks):
            for block in block_decomposition(list_blocks[idx1], list_blocks[idx2]):
                decomposed_blocks.add(block)

        return decomposed_blocks
