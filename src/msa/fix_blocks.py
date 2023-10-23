""" 
Given a list of maximal blocks and an MSA, modify the MSA with new characters in the position covered by the maximal blocks
"""
import json 

from collections import defaultdict
from copy import deepcopy
from typing import Union
from pathlib import Path

from Bio.Align import MultipleSeqAlignment as MSA
from Bio import AlignIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

# pangeblocks
from ..blocks import Block

CHAR="X"

def load_blocks(filename: Union[str,Path]):
    filename = Path(filename)
    extension = filename.suffix
    if extension == ".json":
        with open(filename, "r") as fp:
            blocks = json.load(fp)

    elif extension == ".txt":
        blocks = []
        with open(filename, "r") as fp:
            for line in fp.readlines:
                block = eval(line.replace("\n"))
                blocks.append(block)
                
    else:
        raise Exception("Valid format for blocks are .json and .txt")

    # parse blocks as Block instances
    return [Block(*b) for b in blocks]

def modify_msa(msa: Union[str, Path, MSA], blocks_to_fix: Union[str, Path, list]):

    if type(msa) in (str,Path):
        msa = AlignIO.read(msa, "fasta")
    n_seqs = len(msa)

    if type(blocks_to_fix) in (str, Path):
        blocks_to_fix = load_blocks(blocks_to_fix)
    
    # collect intervals of columns for each row to modify sequences of the MSA
    interval_by_row = defaultdict(list)
    for block in blocks_to_fix:
        for r in block.K:
            interval_by_row[r].append([block.start,block.end])

    # create new MSA
    records = []
 
    for row in range(n_seqs):
        if row in interval_by_row:
            seq = MutableSeq(msa[row].seq)
            for cols in interval_by_row[row]:
                seq[cols[0]:cols[1]+1] = "X"*(cols[1]-cols[0]+1)
            records.append(
                SeqRecord(
                    seq=Seq(seq), id=msa[row].id, name=msa[row].name, description=msa[row].description 
                    )
            )
    
    return MSA(records)