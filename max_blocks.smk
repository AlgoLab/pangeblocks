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
from src.blocks import BlockAnalyzer, Block, block_decomposition, Decomposer
from dataclasses import astuple
PATH_STATS="out/analysis-msa/stats_msas.tsv"
PATH_MSAS=config["PATH_MSAS"]
MAX_MSAS=config["MAX_MSAS"]
RANDOMIZE=config["RANDOMIZE"]

def read_stats_msa(path_stats):

    df = pd.read_csv(path_stats, sep="\t",index_col=False)# usecols=["path_msa","n_unique_seqs"])
    df.sort_values("n_cols", inplace=True)
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
        # expand("out/blocks/{name_msa}.json",name_msa=NAMES)
        expand("out/decomposed-blocks/{name_msa}.json",name_msa=NAMES)

# --- Compute blocks
mv_blocks = MonitorValuesPlus(
            list_vars=[
                        "path_msa","path_blocks","n_seqs",
                        "n_unique_seqs","n_cols","n_max_blocks","t", 
                        "blocks_with_overlap", "inter_between_blocks"
                        ],
            out_file="out/smk_max_blocks.tsv",
            overwrite=False
            )

rule compute_blocks:
    input: 
        "out/analysis-msa/stats_msas.tsv"
    output:
        "out/blocks/{name_msa}.json",
    run:   
        
        # load MSA, count seqs and columns
        path_msa=path_msa_from_name(f"{wildcards.name_msa}")
        # try:
        ti = time.time()
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)
        seqs = list(set([str(record.seq) for record in align]))
        n_unique_seqs = len(seqs)
        
        # compute max blocks and count them
        max_blocks = compute_max_blocks(seqs, n_cols)
        n_max_blocks = len(max_blocks)

        # save max blocks in a file linked to the MSA's name
        path_blocks = f"out/blocks/{wildcards.name_msa}.json" 
        with open(path_blocks,"w") as fp:
            json.dump(max_blocks,fp)
        
        # analyze list of blocks
        dict_analysis=block_analyzer(max_blocks)
        blocks_with_overlap=dict_analysis["blocks_with_overlap"]
        inter_between_blocks=dict_analysis["inter_between_blocks"]
        t = time.time()-ti  

        # log values
        mv_blocks()
        # except:
        #     pass
            # n_cols, n_seqs, n_unique_seqs, n_max_blocks, blocks_with_overlap, inter_between_blocks = None, None, None, None, None, None
       
        

# --- Decompose blocks
mv_decompose_blocks = MonitorValuesPlus(
            list_vars=["path_blocks","path_decomposed_blocks","n_max_blocks","n_decomposed_blocks","t"],
            out_file="out/smk_decomposed_blocks.tsv",
            overwrite=False
            )
block_analyzer=BlockAnalyzer()

rule decompose_blocks:
    input:
        "out/blocks/{name_msa}.json"
    output: 
        "out/decomposed-blocks/{name_msa}.json"
    run:
        ti = time.time() 
        # load blocks
        path_blocks=f"out/blocks/{wildcards.name_msa}.json"
        list_blocks = BlockAnalyzer()._load_list_blocks(path_blocks)        
        n_max_blocks=len(list_blocks)

        # decompose blocks
        decomposer=Decomposer()
        new_blocks = decomposer(list_blocks)
        n_decomposed_blocks = len(new_blocks)
        t = time.time()-ti

        # save decomposed blocks in a file linked to the MSA's name 
        path_decomposed_blocks=f"out/decomposed-blocks/{wildcards.name_msa}.json"
        with open(path_decomposed_blocks,"w") as fp:
            json.dump([astuple(block) for block in new_blocks],fp)

        # log values
        mv_decompose_blocks()