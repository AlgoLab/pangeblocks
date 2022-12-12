import json
import numpy as np
from typing import Optional, Union
from pathlib import Path
from .block import Block

# TODO: sort blocks and avoid using the matrix fro all vs all comparison change it for a list
class BlockAnalyzer:
    """Compute some stats for a list of blocks"""

    def __call__(self, list_blocks: Optional[list[Block]]=None, path_blocks: Optional[Union[str, Path]]=None, **kwargs) -> dict:
        """Compute statistics for a list of blocks.
            If a list of blocks is not provided

        Args:
            list_blocks (Optional[list[Block]], optional): a collection of elements of the type Block. Defaults to None.
            path_blocks (Optional[Union[str, Path]], optional): a path to a json file containing blocks (will be parsed as Block when loaded). 
                                                                Defaults to None.

        Returns:
            dict: number of blocks, number of blocks with overlap and number of intersections between blocks
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
    def _blocks_with_overlap(inter_blocks: np.ndarray) -> int:
        "number of blocks that has at least one overlap"
        blocks_rows, blocks_cols = np.where(inter_blocks>0)
        blocks_with_overlap = set(blocks_rows).union(set(blocks_cols))
        n_blocks_with_overlap = len(blocks_with_overlap)
        return n_blocks_with_overlap

    @staticmethod
    def _inter_between_blocks(inter_blocks: np.ndarray) -> int:
        "number of intersections between pairs of blocks"
        inter_between_blocks=(inter_blocks).sum().sum()
        return int(inter_between_blocks)
    
    @staticmethod
    def _matrix_inter_blocks(list_blocks: list[Block]) -> np.ndarray:
        "matrix to count overlaps by pairs of blocks"
        n_blocks=len(list_blocks)
        inter_blocks=np.zeros((n_blocks, n_blocks))

        # overlap between blocks
        for l,bl in enumerate(list_blocks[:-1]):
            for m, bm in enumerate(list_blocks[l+1:]):
                
                if set(bl.K).intersection(set(bm.K)) and \
                    set(range(bl.i,bl.j+1)).intersection(set(range(bm.i,bm.j+1))) and \
                    bl != bm:
                    inter_blocks[l, m+l+1] += 1

        return inter_blocks

    @staticmethod
    def _load_list_blocks(path_list_blocks: Union[str,Path]) -> list[Block]:
        "load blocks in a list from a json file"
        with open(path_list_blocks,"r") as fp:
            list_blocks=[Block(*args) for args in json.load(fp)]
        return list_blocks