"""
Compute maximal blocks
"""
import time 
import pandas as pd
from pathlib import Path
from Bio import AlignIO #load MSA
# from src.utils.monitor_values import MonitorValues
from src.utils.monitor_values_plus import MonitorValuesPlus
from src.blocks.blocks_suffix_tree import compute_max_blocks
from rich.progress import track

mv_blocks = MonitorValuesPlus(
            list_vars=["path_msa","n_seqs","n_unique_seqs","n_cols","n_max_blocks","t_max_blocks"],
            out_file="out/stats_max_blocks.tsv",
            overwrite=False
            )

df = pd.read_csv("out/stats_msas.tsv", sep="\t", usecols=["path_msa","n_unique_seqs"])
list_msa = df[df["n_unique_seqs"]>1]["path_msa"][10000:10050]

t_max_blocks = 0
for path_msa in track(list_msa, description="Working on MSAs"):
    if t_max_blocks > 60*3:
        break
    ti = time.time()
    # load MSA, count seqs and columns
    align=AlignIO.read(path_msa, "fasta")
    n_cols = align.get_alignment_length()
    n_seqs = len(align)
    seqs = list(set([str(record.seq) for record in align]))
    n_unique_seqs = len(seqs)

    # compute max blocks and count them
    max_blocks = compute_max_blocks(seqs, n_cols)
    n_max_blocks = len(max_blocks)
    # TODO: save max blocks in a file linked to the MSA's name 
    t_max_blocks = time.time()-ti
    mv_blocks()