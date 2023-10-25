"""
Given a set of blocks, compute the value of the objective functions under this set
"""
import argparse
import json
from Bio import AlignIO

from pathlib import Path

def main(args):

    PATH_MSA=args.filename
    PATH_BLOCKS=args.path_blocks

    MIN_LEN=args.min_len           # minimum length of blocks (after removing the '-' from the label of the block)
    MIN_COVERAGE=args.min_coverage # minimum number of out of the total total number of sequences (N_SEQS is needed to compute this)
    PENALIZATION=args.penalization # for [weighted] blocks shorter than MIN_LEN and 


    with open(PATH_BLOCKS, "r") as fp:
        blocks = json.load(fp)

    results = dict()

    # [Nodes]: minimize number of blocks
    if args.nodes:
        nodes = len(blocks)
        results["nodes"] = nodes

    # [Strings]: minimize total length of the graph (sum of labels without '-')
    if args.strings: 
        strings = sum([ 
                        len(b[-1].replace("-","")) 
                        for b in blocks
                    ])
        results["strings"] = strings

    # [Weighted]: same as strings but penalizing nodes with label shorter than MIN_LEN
    if args.weighted: 
        weighted = sum([ 
                        PENALIZATION if len(b[-1].replace("-","")) < MIN_LEN else 1
                        for b in blocks 
                    ])
        results["weighted"] = weighted

    # [Depth]
    if args.depth:
        
        msa = AlignIO.read(PATH_MSA, "fasta")
        N_SEQS=len(msa)

        depth = sum([
                        1 if len(b[0])/N_SEQS > MIN_COVERAGE else PENALIZATION
                        for b in blocks
                    ])
        results["depth"] = depth
    
    return results

if __name__== "__main__":

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="MSA filename")
    parser.add_argument("--path-blocks", help="json file with blocks", dest="path_blocks")
    parser.add_argument("-o","--output", help="json file to save results", dest="output")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")
    
    # parameteres of objective functions
    parser.add_argument("--min-len", default=50, dest="min_len")
    parser.add_argument("--min-coverage", default=0.1, dest="min_coverage")
    parser.add_argument("--penalization", default=128, dest="penalization")
    parser.add_argument("--nodes", action="store_true", dest="nodes")
    parser.add_argument("--strings", action="store_true", dest="strings")
    parser.add_argument("--weighted", action="store_true", dest="weighted")
    parser.add_argument("--depth", action="store_true", dest="depth")
    args = parser.parse_args()

    results = main(args)
    
    if args.output:
        Path(args.output).parent.mkdir(exist_ok=True, parents=True)
        with open(args.output, "w") as fp:
            json.dump(results, fp, indent=1)
    
    print(results)