# TODO: move to src folder
import argparse
import json
import time
from collections import defaultdict
from dataclasses import astuple
from pathlib import Path
from src.blocks import BlockAnalyzer, Block, block_decomposition, Decomposer
from src.utils import MonitorValuesPlus
import logging

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="max blocks filename")
parser.add_argument("--output", help="output file")
parser.add_argument("--output-stats", help="output file stats")
parser.add_argument("--log_level", default='ERROR', dest='log_level',
                    help='set log level (ERROR/WARNING/INFO/DEBUG')
args = parser.parse_args()
logging.basicConfig(level=args.log_level,
                    format='%(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

mv_decomposed_blocks = MonitorValuesPlus(
    list_vars=[
        "filename", "total_blocks"
    ],
    out_file=args.output_stats,
    overwrite=False
)


def list_intersecting_blocks(list_blocks: list[Block]) -> list[tuple]:
    """returns a the tuples(idx_block1, idx_block2) for each pair
    of blocks that intersects.
    It is implemented as a generator, since the list might grow too much."""

    # We organize the blocks in two dictionaries:
    # blocks_ending: for each column, the blocks that end in that column
    # block_by_start: for each column, the blocks that start in that column
    blocks_ending = defaultdict(set)
    block_by_start = defaultdict(list)
    block_by_end = defaultdict(list)
    for idx, block in enumerate(list_blocks):
        block_by_start[block.i].append(idx)
        block_by_end[block.j].append(idx)
    for column in block_by_end.keys():
        blocks_ending[column] = set(block_by_end[column])
    n_cols = max(block_by_end.keys()) + 1
    # We will iterate over the columns of the MSA.
    # We maintain a set of blocks that span the current column
    # (current_blocks) and a set of blocks that have been visited previously
    # (visited_blocks) and that span the current column.
    # For each column, we will
    # (1) compute the intersections between blocks starting in that column
    # and the blocks in visited_blocks
    # (2) remove from visited_blocks the blocks that end in the current
    # column (so that the next iteration will not consider them)
    visited_blocks = set()
    for column in range(n_cols):
        logging.info(f"current column: {column}")
        # notice that in block_by_start[column] we have the blocks that
        # start in the current column.
        # We are going to save them in the current_blocks set
        current_blocks = set([(idx, list_blocks[idx])
                              for idx in block_by_start[column]])
        logging.debug(f"current blocks: {current_blocks}")

        # (1) compute the intersections between blocks in block_by_start[column] and
        # blocks in visited_blocks
        # Notice that all visited blocks span the current column, so
        # there is no need to check if they intersect in the current column
        for idx1, block1 in current_blocks:
            k1 = set(block1.K)
            logging.debug(f"current block: {idx1}, {block1.str()}")
            for idx2, block2 in visited_blocks:
                logging.debug(f"first block: {idx1}, {block1.str()}")
                logging.debug(f"second block: {idx2}, {block2.str()}")
                logging.debug(f"sets: {k1}, {set(block2.K)}")
                if not k1.isdisjoint(set(block2.K)):
                    yield (idx1, idx2)
                    logging.debug(
                        f"intersect: {idx1}, {block1.str()}, {idx2}, {block2.str()}")
        # Add all blocks in current_blocks to the set of visited blocks
        visited_blocks |= current_blocks
        #  Now that we have processed all the blocks that start in the
        #  current column, we can remove the blocks that end in the current
        #  column from the list of visited blocks
        visited_blocks -= blocks_ending[column]


block_analyzer = BlockAnalyzer()
list_blocks = block_analyzer._load_list_blocks(args.filename)
logging.debug(f"len(list_blocks)={len(list_blocks)}, {list_blocks}")
# decompose blocks

decomposed_blocks = set(list_blocks)
for (idx1, idx2) in list_intersecting_blocks(list_blocks):
    for block in block_decomposition(list_blocks[idx1], list_blocks[idx2]):
        decomposed_blocks.add(block)

logging.debug(
    f"len(decomposed_blocks)={len(decomposed_blocks)}, {decomposed_blocks}")


Path(args.output).parent.mkdir(parents=True, exist_ok=True)
print(args.output)
with open(args.output, "w") as fp:
    json.dump([astuple(block) for block in decomposed_blocks], fp)

# stats for decomposed blocks
if len(decomposed_blocks) > 0:
    filename = args.output  # filename of decomposed blocks
    total_blocks = len(decomposed_blocks)
    mv_decomposed_blocks()
