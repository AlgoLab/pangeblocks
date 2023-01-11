import argparse
from ilp.postprocessing import postprocessing

def main():

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_gfa", help="path to GFA")
    parser.add_argument("--path_save", help="path to save postprocess GFA")
    args = parser.parse_args()

    postprocessing(args.path_gfa, args.path_save)

if __name__=="__main__":
    main()