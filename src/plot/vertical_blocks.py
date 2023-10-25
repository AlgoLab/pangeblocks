"Plot vertical blocks in an MSA"
import json
import numpy as np
from Bio import AlignIO
from PIL import Image
from pathlib import Path
from typing import Optional, Union
import logging


_Path = Union[str, Path]

def vertical_blocks(
                    path_msa: _Path, 
                    path_blocks: _Path, 
                    path_submsas: Optional[_Path] = None, 
                    path_img: Optional[_Path] = None, 
                    alpha: int=1, 
                    rows_plot=100,
                    gray=True,
                ):
    """Plot a matrix of the size of the MSA showing vertical blocks in white

    Args:
        path_msa (_type_): fasta file with the MSA
        path_blocks (_type_): path to json file with set of blocks to be considered. Only vertical blocks will be used.
        path_img (_type_): path where to save the image
        alpha (int, optional): minimum length (columns) threshold to consider a vertical block. Defaults to 1.
    """    
    if gray:
        COLOR_VERT_BLOCK = 0     # black (one column subMSA, no ILP, greedy solution)
        COLOR_SUBMSA     = 125   # gray  (ILP acts)
        COLOR_ONE_COL    = 255   # white (fixed blocks of length at least alpha)
    else:
        COLOR_VERT_BLOCK = (255,255,255)    # white (fixed blocks of length at least alpha)
        COLOR_SUBMSA     = (128,128,128)    # gray  (ILP acts)
        COLOR_ONE_COL    = (204,0,0)      # red (one column subMSA, no ILP, greedy solution)


    logging.basicConfig(level= "DEBUG",
                    format=f'[plot-alpha{alpha}] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S'
                    )
    msa=AlignIO.read(path_msa, "fasta")
    n_cols = msa.get_alignment_length()
    n_seqs = len(msa)

    with open(path_blocks) as fp:
        blocks = json.load(fp) 

    # create image. 0 (black) subMSA, 255 (white) vertical block
    if gray:
        numpy_image = np.full((1, n_cols), COLOR_SUBMSA)
    else:
        channel = np.full((1, n_cols), COLOR_SUBMSA[0])
        numpy_image = np.dstack((channel.copy(), channel.copy(), channel.copy())) 

    # Filter vertical blocks of len at least alpha
    forbidden_cols = []
    for b in blocks:
        try:
            K, start, end, label = b
        except:
            K, start, end = b
        
        if len(K) == n_seqs and (end-start+1) >= alpha:
            if gray: 
                numpy_image[:,start:end+1] = COLOR_VERT_BLOCK
            else:
                for col in range(start, end+1):
                    numpy_image[:,col,:] = COLOR_VERT_BLOCK
            forbidden_cols.extend(range(start,end+1,1))

    # ONE COLUMN SUBMSAS
    if path_submsas:
        with open(path_submsas) as fp:
            for line in fp.readlines():
                start, end = line.replace("\n","").split("\t")
                start, end=int(start), int(end)
                if start == end and start not in forbidden_cols:
                    logging.info(f"column between vertical blocks {start} ")
                    if gray:
                        numpy_image[:, start] = COLOR_ONE_COL
                    else:
                        numpy_image[:,start,:] = COLOR_ONE_COL


    # generate array for image
    img_array = np.tile(numpy_image, rows_plot).reshape(rows_plot,-1)
    img = Image.fromarray(img_array.astype('uint8')) #, 'RGB')
    if path_img: 
        img.save(path_img)
    
    return img