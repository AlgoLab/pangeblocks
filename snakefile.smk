

import json
import time
import random; random.seed(0)
import pandas as pd
from pathlib import Path

from Bio import AlignIO #load MSA
from src.utils.monitor_values_plus import MonitorValuesPlus
from src.blocks.blocks_suffix_tree import compute_max_blocks


df = pd.read_csv("out/stats_msas.tsv", sep="\t", usecols=["path_msa","n_unique_seqs"])
list_msa = df[df["n_unique_seqs"]>1]["path_msa"][:30]

# list_msa=list(Path("/home/disco/Data/pandora-msas/msas/").rglob("*fa"))
# random.shuffle(list_msa)

NAMES=[Path(path).stem for path in list_msa]

path_msa_from_name = lambda name_msa: f"/home/disco/Data/pandora-msas/msas/{name_msa}.fa"

# Monitor Values

mv_blocks = MonitorValuesPlus(
            list_vars=["path_msa","n_seqs","n_unique_seqs","n_cols","n_max_blocks","t_max_blocks"],
            out_file="out/smk_max_blocks.tsv",
            overwrite=False
            )

rule all:
    input:
        expand("out/blocks/{msa_name}.json",msa_name=NAMES)

rule compute_blocks:
    input: 
        "/home/disco/Data/pandora-msas/msas/{name_msa}.fa"
    output:
        "out/blocks/{name_msa}.json"
    run:   
        ti = time.time()
        # load MSA, count seqs and columns
        path_msa=path_msa_from_name(f"{wildcards.name_msa}")
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)
        seqs = list(set([str(record.seq) for record in align]))
        n_unique_seqs = len(seqs)

        # compute max blocks and count them
        max_blocks = compute_max_blocks(seqs, n_cols)
        n_max_blocks = len(max_blocks)

        t_max_blocks = time.time()-ti        
        mv_blocks()
        # TODO: save max blocks in a file linked to the MSA's name 
        with open(f"out/blocks/{wildcards.name_msa}.json","w") as fp:
            json.dump(max_blocks,fp)
