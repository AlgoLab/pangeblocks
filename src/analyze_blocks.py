import argparse
import json 
from pathlib import Path
from blocks import Block
from blocks.analyzer import BlockAnalyzer
from utils.timer_decorator import timer
from utils import MonitorValuesPlus

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="max blocks filename")
parser.add_argument("--output", help="output file")
args = parser.parse_args()

filename = args.filename
name_msa = Path(filename).stem
mv_blocks = MonitorValuesPlus(
            list_vars=[
                "filename","total_blocks","blocks_with_overlap", "inter_between_blocks"
            ],
            out_file=args.output,
            overwrite=False
            )

with open(filename) as fp:
    max_blocks = json.load(fp)

# analyze list of blocks
if len(max_blocks)>0:
    block_analyzer = BlockAnalyzer()
    dict_analysis=block_analyzer([Block(*block) for block in max_blocks])
    total_blocks = len(max_blocks)
    blocks_with_overlap=dict_analysis["blocks_with_overlap"]
    inter_between_blocks=dict_analysis["inter_between_blocks"]

    mv_blocks()