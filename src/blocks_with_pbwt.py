import argparse
from blocks.maximal_blocks.wild_pbwt import compute_maximal_blocks
    
if __name__== "__main__":

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="MSA filename")
    parser.add_argument("-o","--output", help="json or txt file to save maximal blocks", dest="output")
    parser.add_argument("-sc","--start-column", help="First column in the MSA to consider. Default=0", type=int, default=0, dest="start_column")
    parser.add_argument("-ec","--end-column", help="Last column in the MSA to consider. Default=-1", type=int, default=-1, dest ="end_column")


    args = parser.parse_args()
    blocks = compute_maximal_blocks(
        msa=args.filename, 
        output=args.output, 
        start_column=args.start_column, 
        end_column=args.end_column, 
        label_blocks=False,
        )