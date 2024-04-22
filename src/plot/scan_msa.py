"""
Plot vertical blocks:
- in black those fixed by alpha
- in white those chosen by the ILP
"""

from Bio import AlignIO
from PIL import Image

import numpy as np
import matplotlib.pyplot as plt

BLACK = 0
GRAY  = 125
WHITE = 230

def vertical_blocks(blocks_alpha, blocks_optimum , path_msa, rowsfig=100, figsize=(8,3), return_pil=False): 
    """blocks alpha will be colored with Black, blocks optimum by White color,
    and the remaining part of the MSA will be colored with gray
    Assumption: 
    - blocks_alpha are vertical blocks
    - blocks_optimum can contain not only vertical blocks, but they are not considered
    """
    msa = AlignIO.read(path_msa ,"fasta")
    cols = msa.get_alignment_length()
    rows = len(msa)
    panel = np.full((1, cols), GRAY)

    # fill with black color those columns covered by vertical blocks chosen by alpha
    for b in blocks_alpha:
        K, start, end, label = b
        panel[:, start:end+1] = BLACK

    # and with white color those columns covered by vertical blocks chosen by the ILP
    for b in blocks_optimum:
        K, start, end, label = b 
        if len(K) == rows:
            panel[:, start:end+1] = WHITE


    img_array = np.tile(panel,(rowsfig,1))

    if return_pil:
        # to generate a subplot with many of these
        return Image.fromarray(img_array.astype('uint8'))

    f, ax = plt.subplots(figsize=figsize)
    ax.imshow(img_array, cmap='gray', vmin=0, vmax=255)
    ax.set_yticks([])
    return ax
