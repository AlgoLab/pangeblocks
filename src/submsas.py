#!/usr/bin/env python3
import json
import logging
import argparse 
from pathlib import Path
from Bio import AlignIO


def info_msa(filename):
    "Return MSA from start_column to end_column (both included)"
    # load MSA
    msa = AlignIO.read(filename, "fasta")
    ncols = msa.get_alignment_length()
    nrows = len(msa)
    return nrows, ncols

def split_submsa(max_positions_msa, nrows, start_column, end_column):
    """Divide a subMSA in smaller subMSAs based on the number of positions it can have given by
    max_positions_msa. Returns a list with new subMSAs
    """
    cols_submsa = max_positions_msa // nrows
    end_submsa = 0 
    start_end = []

    start_submsa = start_column
    while end_submsa < end_column:  
        end_submsa = start_submsa + cols_submsa - 1 
        if end_submsa > end_column:
            end_submsa = end_column
        print(start_submsa, end_submsa)
        start_end.append(
            [start_submsa, end_submsa]
        )
        # update next start of the submsa
        start_submsa = end_submsa + 1 

    return start_end
    
if __name__=="__main__":
    """
    Return plain text with starting and ending column for each subMSA to be computed.
    This is obtained from positions in between vertical blocks
    """

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-vertical-blocks", help="path to vertical_blocks.json", dest="path_vertical_blocks")
    parser.add_argument("--path-msa", help="path to msa  in fasta format", dest="path_msa")
    parser.add_argument("-o","--output", help="output file .txt", dest="output")
    parser.add_argument("--threshold-vertical-blocks", help="vertical blocks with length greather or equal than the threshold \
                        will be considered to split the MSA", type=int, dest="threshold_vertical_blocks", default=1)
    parser.add_argument("--max-positions-msa", help="subMSAs with more than this number of positions will be divided into \
                        smaller subMSAs", type=int, dest="max_positions_msa", default=100000)
    parser.add_argument("--log-level", default='ERROR', help="set log level (ERROR/WARNING/INFO/DEBUG). Default ERROR", dest="log_level")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                    format='[submsas] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
    logging.info(f"filename MSA: '{args.path_msa}'")
    logging.info(f"filename vertical blocks: '{args.path_vertical_blocks}'") 
    # logging.info(f"subMSA columns [start, end] = {[args.start_column, args.end_column]}")
    # logging.info(f"Output only vertical blocks = {args.only_vertical}")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.path_vertical_blocks, "r") as fp:
        vertical_blocks = [vb for vb in json.load(fp) if vb[2]-vb[1]+1 >= args.threshold_vertical_blocks] 

    nrows, ncols = info_msa(args.path_msa)

    if len(vertical_blocks) == 0:
        
        with open(args.output, "w") as fp: 
            
            start = 0 # start MSA  
            end   = -1 # end MSA 
            logging.info("No Vertical Blocks", start, end)
            fp.writelines(f"{start}\t{end}\n")

    else:   
        # sort vertical blocks by starting position
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
        with open(args.output, "w") as fp: 
            for j, pair in enumerate(pairs_blocks):
                left_block, right_block = pair

                start = left_block[2] + 1 # start submsa = end of left block + 1  
                end   = right_block[1] - 1 # end submsa = start of right block - 1

                submsas = split_submsa(max_positions_msa=args.max_positions_msa, nrows=nrows, start_column=start, end_column=end)
                for submsa in submsas:
                    _start, _end = submsa
                    logging.info("Vertical Block",left_block, right_block, _start, _end)
        
                    fp.writelines(f"{_start}\t{_end}\n")