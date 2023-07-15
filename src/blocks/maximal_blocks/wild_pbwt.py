"Compute Maximal Blocks from an (sub)MSA using wild-pBWT"
import os
import tempfile
import subprocess
import json
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional
from src.blocks.maximal_blocks.utils import (
    load_submsa,
) 

import logging

def compute_maximal_blocks(filename: Union[str,Path], output: Optional[Union[str,Path]] = None, 
                           start_column: int = 0, end_column: int = -1, 
                           only_vertical: bool = False,
                           alphabet_to_ascii: dict = {"-":0,"A":1,"C":2,"G":3,"T":4},
                           bin_wildpbwt: str = "Wild-pBWT/bin/wild-pbwt",
                           label_blocks: bool = True
                           ):
    "Compute maximal blocks in a submsa"
    
    PATH_WILD_PBWT=bin_wildpbwt
    logging.info(f"WILD-PBWT: {PATH_WILD_PBWT}")
    if PATH_WILD_PBWT:
        assert os.path.isfile(PATH_WILD_PBWT), f"it seems that {PATH_WILD_PBWT} is not the correct path to /bin/wild-pbwt"
    SIZE_ALPHABET=len(alphabet_to_ascii)
    
    # load subMSA
    msa=load_submsa(filename, start_column, end_column)
    n_cols=msa.get_alignment_length()
    n_seqs=len(msa)

    # create file to store input matrix for wild-pBWT
    fd, path = tempfile.mkstemp()
    logging.info(f"tmp file [{start_column},{end_column}] {path}")
    try:
        with os.fdopen(fd, 'w') as fileTemp:
            for seq, record in enumerate(msa):
                row_panel=str(record.seq)
                for a, c in alphabet_to_ascii.items():
                    row_panel = row_panel.replace(a,str(c))
                # print(row_panel)
                # write rows in ASCII alphabetto the tmp file
                fileTemp.write(row_panel + "\n")

        # compute maximal blocks with wild-pBWT
        logging.info("")
        COMMANDS = f"{PATH_WILD_PBWT} -a {SIZE_ALPHABET} -f {path} -o y".split(" ")
        # print(COMMANDS)
        result = subprocess.run(COMMANDS, capture_output=True, text=True)
        # print(result.stderr)
        # print(result.stdout)
        max_blocks = [eval(f"[{r}]") for r in result.stdout.split("\n") if len(r)>1]

    finally:
        # print(os.path.exists(path))
        os.remove(path)

    if label_blocks: 
        max_blocks_strings = [
            [b[0], start_column + b[1], start_column + b[2], str(msa[b[0][0], b[1]:b[2]+1].seq) ] for b in max_blocks
        ]
        return max_blocks_strings #, msa
    
    return max_blocks #, msa