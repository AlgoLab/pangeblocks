"""
Compute maximal blocks in subalignments of an MSA by considering 
only consecutive columns that contains more than one character in {A,C,G,T}
"""
# params
configfile: "params.yaml"

PATH_STATS="out/analysis-msa/stats_msas.tsv"
PATH_MSAS=config["PATH_MSAS"]
MAX_MSAS=config["MAX_MSAS"]
RANDOMIZE=config["RANDOMIZE"]

# libs
import json
import time
import random; random.seed(0)

from dataclasses import astuple
from Bio import AlignIO
from pathlib import Path

from src.utils.monitor_values_plus import MonitorValuesPlus
from src.blocks.blocks_suffix_tree import compute_max_blocks
from src.blocks import (
    BlockAnalyzer, 
    Block, 
    Decomposer
)
from src.msa import AnalyzerMSA

from src.utils_smk import (
    read_stats_msa, 
    split_vec_by_consecutive_values
)

def path_msa_from_name(name_msa):
    return Path(PATH_MSAS).joinpath(f"{name_msa}.fa")

# read msas to process
list_msa = read_stats_msa(path_stats=PATH_STATS, randomize=RANDOMIZE, max_msas=MAX_MSAS)
print(list_msa[0])
NAMES=[Path(path).stem for path in list_msa] # msas names

rule all: 
    input:
        expand("out-sub/blocks/{name_msa}.json", name_msa=NAMES)

# Compute blocks
mv_blocks = MonitorValuesPlus(
            list_vars=[
                "path_msa","path_blocks","n_seqs",
                "n_unique_seqs","n_cols","n_max_blocks","t", 
                "blocks_with_overlap", "inter_between_blocks"
            ],
            out_file="out-sub/smk_max_blocks.tsv",
            overwrite=False
            )

block_analyzer=BlockAnalyzer()

rule compute_blocks:
    input:
        "out/analysis-msa/stats_msas.tsv"
    output:
        "out-sub/blocks/{name_msa}.json"
    run:
        ti = time.time() # time until max blocks are computed
        # load MSA, count seqs and columns
        path_msa = path_msa_from_name(f"{wildcards.name_msa}")
        align  = AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)
        seqs   = list(set([str(record.seq) for record in align]))
        n_unique_seqs = len(seqs)

        # identify sub-alignments
        analyzer = AnalyzerMSA()
        bool_id_cols = analyzer.check_identical_columns(seqs, n_cols)

        # get position of consecutive values
        pos_splits = split_vec_by_consecutive_values(bool_id_cols)
        
        pos_sub_msas = []
        for start, end in pos_splits:
            # False (0) means cols are not identical
            # and we will compute max block in the sub-msa
            if sum(bool_id_cols[start:end+1])==0: 
                pos_sub_msas.append((start,end))
        
        # compute max blocks and count them
        # TODO: run in parallel
        max_blocks = []
        if len(pos_sub_msas)>0:
            for start_sub_msa, end_sub_msa in pos_sub_msas:
                seqs_sub_msa = [seq[start:end+1] for seq in seqs]
                n_cols_sub_msa = end-start+1
                max_blocks_sub_msa=compute_max_blocks(seqs,n_cols_sub_msa)
                n_max_blocks_sub_msa=len(max_blocks_sub_msa)

                max_blocks.extend(max_blocks_sub_msa)
        t = time.time()-ti
        n_max_blocks = len(max_blocks)
        path_blocks=f"out-sub/blocks/{wildcards.name_msa}.json"
        with open(path_blocks,"w") as fp:
            json.dump(max_blocks, fp)

        # analyze list of blocks
        if len(max_blocks)>0:
            dict_analysis=block_analyzer(max_blocks)
            blocks_with_overlap=dict_analysis["blocks_with_overlap"]
            inter_between_blocks=dict_analysis["inter_between_blocks"]
         

        mv_blocks()