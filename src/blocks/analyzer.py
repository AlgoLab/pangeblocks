import json
import numpy as np
from typing import Optional, Union
from pathlib import Path
from .block import Block
from collections import defaultdict
import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

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
        inter_blocks=self._list_inter_blocks(list_blocks)

        # stats
        return dict(
            number_of_blocks=len(list_blocks),
            blocks_with_overlap=self._blocks_with_overlap(inter_blocks),
            inter_between_blocks=self._inter_between_blocks(inter_blocks),
        )

    @staticmethod
    def _blocks_with_overlap(inter_blocks: np.ndarray) -> int:
        "number of blocks that has at least one overlap"
        overlapping_blocks = []
        for idx_blocks in inter_blocks:
            overlapping_blocks.extend(idx_blocks)
        blocks_with_overlap = set(overlapping_blocks)

        return len(blocks_with_overlap)

    @staticmethod
    def _inter_between_blocks(inter_blocks: np.ndarray) -> int:
        "number of intersections between pairs of blocks"
        return len(inter_blocks)

    @staticmethod
    def _load_list_blocks(path_list_blocks: Union[str,Path]) -> list[Block]:
        "load blocks in a list from a json file"
        with open(path_list_blocks,"r") as fp:
            list_blocks=[Block(*args) for args in json.load(fp)]
        return list_blocks