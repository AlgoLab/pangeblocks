from importlib.resources import path
import json
import numpy as np
import pandas as pd
from typing import Optional, Union
from pathlib import Path
from .block import Block


class BlockStats:
    """Compute some stats for a list of blocks"""

    def __call__(self, list_blocks: Optional[list[Block]]=None, path_blocks: Optional[Union[str, Path]]=None, **kwargs) -> dict:
        """Compute statistics for a set of blocks
        """
        
        if path_blocks is not None:
            list_blocks = self._load_list_blocks(path_blocks)

        # compute matrix with intersections
        inter_blocks=self._matrix_inter_blocks(list_blocks)

        # stats
        return dict(
            number_of_blocks=len(list_blocks),
            blocks_with_overlap=self._blocks_with_overlap(inter_blocks),
            inter_between_blocks=self._inter_between_blocks(inter_blocks),
        )
        
    @staticmethod
    def _blocks_with_overlap(inter_blocks) -> int:
        "number of blocks that has at least one overlap"
        blocks_with_overlap=(inter_blocks.max(axis=1)>0).sum()
        return int(blocks_with_overlap)

    @staticmethod
    def _inter_between_blocks(inter_blocks) -> int:
        "number of intersections between pairs of blocks"
        inter_between_blocks=(inter_blocks).sum().sum()
        return int(inter_between_blocks)
    
    @staticmethod
    def _matrix_inter_blocks(list_blocks) -> np.ndarray:
        "matrix to count overlaps by pairs of blocks"
        n_blocks=len(list_blocks)
        inter_blocks=np.zeros((n_blocks, n_blocks))

        # overlap between blocks
        for l,bl in enumerate(list_blocks[:-1]):
            for m, bm in enumerate(list_blocks[1:]):
            
                if set(bl.K).intersection(set(bm.K)) and \
                    set(range(bl.i,bl.j+1)).intersection(set(range(bm.i,bm.j+1))) and \
                    bl != bm:
                    inter_blocks[l, m+1] += 1

        return inter_blocks

    @staticmethod
    def _load_list_blocks(path_list_blocks: Union[str,Path])->list[Block]:
        "load blocks in a list from a json file"
        with open(path_list_blocks,"r") as fp:
            list_blocks=[Block(*args) for args in json.load(fp)]
        return list_blocks