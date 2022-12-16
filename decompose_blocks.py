# TODO: move to src folder
import argparse
import json
import time
from dataclasses import astuple
from pathlib import Path
from src.blocks import BlockAnalyzer, Block, block_decomposition, Decomposer
from src.utils import MonitorValuesPlus

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="max blocks filename")
parser.add_argument("--output", help="output file")
parser.add_argument("--output-stats", help="output file stats")
args = parser.parse_args()

mv_decomposed_blocks = MonitorValuesPlus(
            list_vars=[
                "filename", "total_blocks"
            ],
            out_file=args.output_stats,
            overwrite=False
            )

t1 = time.time()

block_analyzer=BlockAnalyzer()
list_blocks = block_analyzer._load_list_blocks(args.filename)

# decompose blocks
decomposer=Decomposer()
decomposed_blocks = decomposer(list_blocks)
t2 = time.time()
t_decomp_blocks = t2-t1

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
print(args.output)
with open(args.output, "w") as fp:
    json.dump([astuple(block) for block in decomposed_blocks], fp)

# save times
path_time = Path(args.output).parent / (Path(args.output).stem + ".txt")
with open(path_time,"w") as fp:
    fp.write(f"t_decomp_blocks\t{t_decomp_blocks}\n")

# stats for decomposed blocks
if len(decomposed_blocks)>0:
    filename=args.output # filename of decomposed blocks
    total_blocks = len(decomposed_blocks)
    mv_decomposed_blocks()