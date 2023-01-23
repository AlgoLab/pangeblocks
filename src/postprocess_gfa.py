import argparse
# from ilp.postprocessing import postprocessing
from ilp.gfanx import postprocessing

def main():

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_gfa", help="path to GFA")
    # parser.add_argument("--path_save", help="path to save postprocess GFA")
    args = parser.parse_args()

    postprocessing(args.path_gfa)

if __name__=="__main__":
    main()