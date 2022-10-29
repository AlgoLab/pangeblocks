configfile: "params.yaml"

import json
import time
import random; random.seed(0)
import pandas as pd
from pathlib import Path
from rich.progress import track
from Bio import AlignIO #load MSA
from src.utils.monitor_values_plus import MonitorValuesPlus
from src.blocks.blocks_suffix_tree import compute_max_blocks

PATH_STATS="out/stats_msas.tsv"
PATH_MSAS=config["PATH_MSAS"]
MAX_MSAS=config["MAX_MSAS"]
RANDOMIZE=config["RANDOMIZE"]

def read_stats_msa(path_stats):

    df = pd.read_csv(path_stats, sep="\t", usecols=["path_msa","n_unique_seqs"])
    list_msa = df[df["n_unique_seqs"]>1]["path_msa"].tolist()
    if RANDOMIZE is True:
        random.shuffle(list_msa)
    if MAX_MSAS:
        return list_msa[:MAX_MSAS]
    return list_msa

def path_msa_from_name(name_msa):
    return Path(PATH_MSAS).joinpath(f"{name_msa}.fa")


list_msa = read_stats_msa(PATH_STATS)
NAMES=[Path(path).stem for path in list_msa]


rule all:
    input:
        expand("out/blocks/{name_msa}.json",name_msa=NAMES)

# --- Compute blocks
mv_blocks = MonitorValuesPlus(
            list_vars=["path_msa","n_seqs","n_unique_seqs","n_cols","n_max_blocks","t"],
            out_file="out/smk_max_blocks.tsv",
            overwrite=False
            )

rule compute_blocks:
    input: 
        "out/stats_msas.tsv"
    output:
        "out/blocks/{name_msa}.json",
    run:   
        ti = time.time()
        # load MSA, count seqs and columns
        path_msa=path_msa_from_name(f"{wildcards.name_msa}")
        try:
            align=AlignIO.read(path_msa, "fasta")
            n_cols = align.get_alignment_length()
            n_seqs = len(align)
            seqs = list(set([str(record.seq) for record in align]))
            n_unique_seqs = len(seqs)
            # compute max blocks and count them
            max_blocks = compute_max_blocks(seqs, n_cols)
            n_max_blocks = len(max_blocks)

        except:
            n_cols, n_seqs, n_unique_seqs, n_max_blocks = None, None, None, None

        t = time.time()-ti        
        mv_blocks()

        # save max blocks in a file linked to the MSA's name 
        with open(f"out/blocks/{wildcards.name_msa}.json","w") as fp:
            json.dump(max_blocks,fp)

# # --- EDA MSA
# rule eda_msa: 
#     input: 
#         PATH_MSAS
#     output:
#         "out/stats_msas.tsv"
#     shell:
#         "python3 eda_msas.py {path_msas}".format(path_msas=PATH_MSAS)