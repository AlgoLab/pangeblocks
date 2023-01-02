from PIL import Image
from Bio import AlignIO
import numpy as np 
from pathlib import Path
from itertools import cycle

COLORS = [
    (255,0,0), # red
    (0,255,0), # green 
    (0,0,255), # blue
    (255,255,0), # yellow
    (0,255,255), # cyan
    (255,0,255), # magenta
    (255,125,0), # orange
    (125,0,255), # blue-magenta
    (255,0,125), # red-magenta
]

class CoverageMSA:

    # def __call__(self, path_msa, blocks, path_save grayscale: bool=True, color: bool=False):
    #     msa, n_seqs, n_cols = self.load_msa(path_msa)
    #     coverage_gray = self.get_coverage_panel(n_seqs, n_cols, blocks)
    #     save_img = 

    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols

    @staticmethod
    def array2img(array):
        "Array to PIL image"
        m, M = array.min(), array.max()
        # rescale to [0,1]
        img_rescaled = (array - m) / (M-m)
        
        # invert colors black->white
        img_array = np.ceil(255 - img_rescaled*255)
        dtype = eval(f"np.int8")
        img_array = np.array(img_array, dtype=dtype)
        
        # convert to Image 
        img_pil = Image.fromarray(img_array,'L')
        return img_pil

    def save_img(self, array, path_save):
        "Save image in grayscale for the FCGR provided as input"
        Path(path_save).parent.mkdir(exist_ok=True, parents=True)
        img_pil = self.array2img(array)
        img_pil.save(path_save)


    def save_colored_img(self, array, path_save):
        img_pil = Image.fromarray(array.astype("uint8"),"RGB")
        img_pil.save(path_save)

    @staticmethod
    def get_coverage_panel(n_seqs, n_cols, blocks):
        """returns a matrix of size equal to msa (n_seq x n_cols) with 
        the number of blocks in the list_blocks that covers each position"""

        # coverage_by_pos = defaultdict(int)
        coverage_panel = np.zeros((n_seqs, n_cols))
        for block in blocks:
            for r in block.K:
                for c in range(block.i,block.j+1):
                    coverage_panel[r,c] += 1
        return coverage_panel

    @staticmethod
    def get_coverage_color_panel( n_seqs, n_cols, blocks, colors=COLORS):
        """returns a matrix of size equal to msa (n_seq x n_cols) with 
        the number of blocks in the list_blocks that covers each position"""
        colors = iter(cycle(colors))
        # coverage_by_pos = defaultdict(int)
        coverage_panel = 255*np.ones((n_seqs, n_cols,3))
        for block in blocks:
            color = next(colors)
            for r in block.K:
                for c in range(block.i,block.j+1):
                    coverage_panel[r,c,:] = color
        return coverage_panel