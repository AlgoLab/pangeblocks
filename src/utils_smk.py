"""
Utilities for snakemake pipelines
"""
from typing import Optional
import random 
import pandas as pd 

def read_stats_msa(path_stats: str, randomize: bool, max_msas: Optional[int]=None):

    df = pd.read_csv(path_stats, sep="\t",index_col=False)
    df.sort_values("n_cols", inplace=True)
    list_msa = df[df["n_unique_seqs"]>1]["path_msa"].tolist()
    if randomize is True:
        random.shuffle(list_msa)
    if max_msas:
        return list_msa[:max_msas]
    return list_msa

def split_vec_by_consecutive_values(vec):
    """split a vector (start:end) by consecutive values"""
    splits=[]
    curr_pos = 0
    start = 0
    end   = 0 

    while curr_pos < len(vec)-1:
        
        if vec[curr_pos] == vec[curr_pos+1]:
            end = curr_pos + 1 
        else:
            splits.append((start,end))
            start = end + 1 
            end = start

        # move one position 
        curr_pos +=1

    # append last consecutive (positions) of values
    splits.append((start, end))

    return splits