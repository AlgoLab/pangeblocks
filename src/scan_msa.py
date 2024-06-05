from Bio import AlignIO
from collections import namedtuple
from pathlib import Path

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from plot.scan_msa import vertical_blocks as plot_vertical_blocks



def load_msa(path_msa):
    with open(path_msa) as fp:
        msa = AlignIO.read(fp, "fasta")
    return msa

def compute_vertical_blocks(msa, threshold_vertical_blocks=1):
    "Compute maximal blocks in a submsa"

    # info subMSA
    n_cols=msa.get_alignment_length()
    n_seqs=len(msa)

    vertical_blocks = []
    
    chars_block = []
    start_col = 0
    for col in range(n_cols):
        chars_col = msa[:,col]
        chars_col = list(set(chars_col))

        if len(chars_col) == 1: 
            chars_block.append(chars_col[0])

        else:
            
            if len(chars_block)>= threshold_vertical_blocks:
                end_col = col-1 # end column is included
                vertical_blocks.append(
                    (list(range(n_seqs)), start_col, end_col, "".join(chars_block))
                )
        
            start_col = col+1 # update starting column for the next iteration
            # update initial values for a new block
            chars_block = []
        
    # special case for block ending in the last column
    if len(chars_block)>= threshold_vertical_blocks:
        end_col = col # end column is included
        vertical_blocks.append(
            (list(range(n_seqs)), start_col, end_col, "".join(chars_block))
        )

    return vertical_blocks

def compute_submsas(vertical_blocks, msa):
    SM = namedtuple("SubMSA",["start","end","nrows","ncols"])
    ncols=msa.get_alignment_length()
    nseqs=len(msa)

    if len(vertical_blocks) == 0:
        start = 0
        end = ncols

        return [SM(start, end, nseqs, ncols)]

    vertical_blocks = sorted(vertical_blocks, key=lambda b: (b[1],b[2]))

    # Include auxiliar first and last block if needed
    first_vb = vertical_blocks[0]
    if first_vb[1] > 0:
        vertical_blocks.insert(0, [[],-1,-1,"N"])
    
    last_vb = vertical_blocks[-1]
    if last_vb[2] < ncols-1: 
        vertical_blocks.append([[],ncols,ncols,"N"])
    
    pairs_blocks = zip(vertical_blocks[:-1], vertical_blocks[1:])

    # save starting and ending positions for each subMSA
    submsas = []
    
    for j, pair in enumerate(pairs_blocks):
        left_block, right_block = pair

        start = left_block[2] + 1 # start submsa = end of left block + 1  
        end   = right_block[1] - 1 # end submsa = start of right block - 1
        ncols = end - start + 1
        submsas.append(SM(start, end, nseqs, ncols))
    return submsas

def collect_info_breakpoints(msa, vertical_blocks):

    ncols=msa.get_alignment_length()
    nrows=len(msa)

    submsas = [1,2]
    alpha = 1

    info = []
    longest_submsa = 0
    alpha_breakpoints = []
    while len(submsas) > 1:
        
        blocks_alpha = [b for b in vertical_blocks if b[2]-b[1]+1 >= alpha]
        if len(blocks_alpha) == 0: 
            info.append(
                    dict(alpha=alpha, longest_submsa=ncols, blocks_longest_submsa=[[range(nrows),0,ncols-1,""]], n_submsas=1)
                )
            alpha_breakpoints.append(alpha)
            break
        submsas = compute_submsas(vertical_blocks=blocks_alpha, msa=msa)

        # 4. get the longest subMSA
        longest_submsa_alpha = max([submsa.ncols for submsa in submsas])
        blocks_longest_submsa = []
        for submsa in submsas:
            if submsa.ncols == longest_submsa_alpha:
                blocks_longest_submsa.append([blocks_alpha[0][0], submsa.start, submsa.end,""])

        if longest_submsa_alpha > longest_submsa:
            longest_submsa = longest_submsa_alpha
            alpha_breakpoints.append(alpha)
        

        info.append(
            dict(alpha=alpha, longest_submsa=longest_submsa_alpha, blocks_longest_submsa=blocks_longest_submsa, n_submsas=len(submsas))
        )
        alpha+=1

    return info, alpha_breakpoints

def plot_longest_submsa_by_alpha(info, dirsave, name):
    PATH_SAVE = Path(dirsave)
    sns.set_style("darkgrid")
    fig, ax = plt.subplots(1,1)
    plot_longest_submsa=pd.DataFrame(info).plot(
        x="alpha",
        y="longest_submsa",
        ax=ax
        )

    plot_longest_submsa.text(x=0.5, y=1.1, s=f"Longest subMSA by {chr(945)} | MSA {name}", fontsize=14, weight='bold', ha='center', va='bottom', transform=plot_longest_submsa.transAxes)
    ax.set_xlabel(f"{chr(945)}")
    ax.set_ylabel("number of columns")
    ax.set_ylim(bottom=0, top = 1.05*(ncols))
    fig.savefig(PATH_SAVE.joinpath("longest_submsa_by_alpha"), dpi=300)


def barplot_longest_submsa_breakpoint(info, alpha_breakpoints, dirsave, name):
    sns.set_style("darkgrid")
    PATH_SAVE=Path(dirsave)
    fig, ax = plt.subplots(1,1)
    data_plot = pd.DataFrame(info).query(f"alpha in {alpha_breakpoints}")

    breakpoints_barplot = data_plot.plot(
        x="alpha",y="longest_submsa", kind="bar", 
        ax=ax
        )

    breakpoints_barplot.text(x=0.5, y=1.1, s=f"Longest subMSA by {chr(945)} | MSA {name}", fontsize=14, weight='bold', ha='center', va='bottom', transform=breakpoints_barplot.transAxes)
    breakpoints_barplot.text(x=0.5, y=1.05, s=f"breakpoints correspond to values of {chr(945)} where the longest subMSA change", fontsize=10, alpha=0.75, ha='center', va='bottom', transform=breakpoints_barplot.transAxes)

    ax.set_xlabel(f"{chr(945)} breakpoints")
    ax.set_ylabel("number of columns longest subMSA")

    plt.xticks(rotation=0)
    fig.savefig(PATH_SAVE.joinpath("barplot_longest_submsa_by_alpha_breakpoints"), dpi=300)


def msa_blocks_submsas_by_alpha(path_msa, alpha_breakpoints, name):

    fig = plt.figure(1,(20,16))
    fig.subplots_adjust(hspace=.4)

    # alpha_breakpoints
    blocks_longest_submsas_by_alpha = dict()
    longest_submsa_by_alpha = dict()
    number_of_submsas_by_alpha = dict()
    for d in list(filter(lambda d: d["alpha"] in alpha_breakpoints, info)):
        longest_submsa_by_alpha[d["alpha"]] = d["longest_submsa"]
        blocks_longest_submsas_by_alpha[d["alpha"]] = d["blocks_longest_submsa"]
        number_of_submsas_by_alpha[d["alpha"]] = d["n_submsas"]
        
    nplots = len(alpha_breakpoints)
    for j, alpha in enumerate(alpha_breakpoints):
        blocks_alpha = [b for b in vertical_blocks if b[2]-b[1]+1 >= alpha] 
        n_blocks_alpha = len(blocks_alpha)
        blocks_longest_submsa = blocks_longest_submsas_by_alpha[alpha]
        n_submsas = number_of_submsas_by_alpha[alpha]
        img = plot_vertical_blocks(blocks_alpha, blocks_optimum=blocks_longest_submsa, path_msa=path_msa, rowsfig=1000, return_pil=True)
        plt.subplot(nplots, 1, j+1)
        plt.title(f"{chr(945)}={alpha} | largest subMSA ( {nrows} x {longest_submsa_by_alpha[alpha]} ) | Fixed Vertical Blocks = {n_blocks_alpha} | subMSAs = {n_submsas}")
        plt.yticks([])
        plt.imshow(np.asarray(img), cmap="gray", vmin=0, vmax=255, interpolation="none")

    plt.suptitle(f"""Vertical Blocks and subMSAs | {name} | MSA ({nrows},{ncols})""", y=1.5, fontsize=18, weight='bold', ha='center', va='bottom', transform=fig.axes[0].transAxes)
    plt.tight_layout()

    fig.savefig(PATH_SAVE.joinpath("vertical_blocks_by_alpha.png"),dpi=300)


if __name__ == "__main__":
    import argparse
    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("path_msa", help="MSA filename",)
    parser.add_argument("-o","--output-dir", help="output directory", dest="outdir")
    args = parser.parse_args()

    PATH_MSA = Path(args.path_msa)
    name = PATH_MSA.stem
    PATH_SAVE = Path(args.outdir).joinpath(name)
    PATH_SAVE.mkdir(exist_ok=True, parents=True)
    
    msa = load_msa(PATH_MSA)
    nrows = len(msa)
    ncols = msa.get_alignment_length()   

    vertical_blocks = compute_vertical_blocks(msa)
    info, alpha_breakpoints = collect_info_breakpoints(msa, vertical_blocks)

    plot_longest_submsa_by_alpha(info, PATH_SAVE, name)

    barplot_longest_submsa_breakpoint(info, alpha_breakpoints, PATH_SAVE, name)

    msa_blocks_submsas_by_alpha(PATH_MSA, alpha_breakpoints, name)