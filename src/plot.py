import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt 
from plot.vertical_blocks import vertical_blocks

def main(args):
    "plot vertical blocks of length at least alpha in the MSA"
    path_msa = args.path_msa
    path_blocks = args.path_blocks
    dir_output = args.dir_output
    # path_img = args.path_img
    path_submsas = args.path_submsas
    alphas = [int(a) for a in  args.alpha]
    rows_plot = args.rows_plot
    one_plot = args.one_plot



    Path(dir_output).mkdir(exist_ok=True, parents=True)
    if one_plot:

        imgs = []
        n_plots = len(alphas)
        # fig, axs = plt.subplots(n_plots,1)
        fig = plt.figure(1,(30,25))
        for j,alpha in enumerate(alphas):
            print("type of alpha", type(alpha))
            name_msa = Path(path_msa).stem 
            img = vertical_blocks(path_msa=path_msa, path_blocks=path_blocks, path_submsas=path_submsas, alpha=alpha, rows_plot=rows_plot)
            # axs[j] = plt.imshow(np.asarray(img), cmap="gray", interpolation="none")
            plt.subplot(n_plots, 1, j+1)
            plt.imshow(np.asarray(img), cmap="gray", interpolation="none")
        
        path_img = Path(dir_output).joinpath(f"{name_msa}-alpha{'_'.join(args.alpha)}.png")
        plt.savefig(path_img, dpi=300)
        
    else:
    
        for alpha in alphas:
            print("type of alpha", type(alpha))
            name_msa = Path(path_msa).stem 

            path_img = Path(dir_output).joinpath(f"{name_msa}-alpha{alpha}.png")
            vertical_blocks(path_msa, path_blocks, path_submsas, path_img, alpha, rows_plot)

if __name__ == "__main__":

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-msa", help="path to GFA", type=str, dest="path_msa")
    parser.add_argument("--path-blocks", help="path to save postprocess GFA", type=str, dest="path_blocks")
    parser.add_argument("--dir-output", help="prefix output. Will be used like 'dir_output/<name_msa>-<alpha>.png'", type=str, dest="dir_output", default=None)
    parser.add_argument("--path-submsas", help="path to submsa index", type=str, dest="path_submsas", default=None)
    parser.add_argument("-a", "--alpha", nargs="+",help="minimum threshold (columns) to consider a vertical block", dest="alpha", default=1)
    parser.add_argument("--one-plot", help="output one figure with all alpha values",action="store_true", dest="one_plot")
    parser.add_argument("--rows-plot", help="number of rows to generate the plot", type=int, dest="rows_plot", default=1000)
    args = parser.parse_args()
    print(args)
    main(args)  