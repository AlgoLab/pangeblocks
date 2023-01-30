from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from src.utils import MonitorValuesPlus
from src.msa import AnalyzerMSA
from tqdm import tqdm
import argparse

def main(path_msas: str, path_output: str):

    mv = MonitorValuesPlus(list_vars=["path_msa","n_seqs","n_unique_seqs","n_identical_cols","n_cols", "perc_identical_cols"],
                            out_file=f"{path_output}/analysis-msa/stats_msas.tsv",
                            overwrite=True)

    mv_problems = MonitorValuesPlus(list_vars=["path_msa"],
                                    out_file=f"{path_output}/analysis-msa/problematic_files.tsv",
                                    overwrite=True)
    
    analyzer = AnalyzerMSA()
    def run(path_msa):
        "Function to run with ThreadPoolExecutor"
        try:
            n_cols, n_seqs, n_unique_seqs, n_identical_cols = analyzer(path_msa)
            perc_identical_cols = round(100*n_identical_cols/n_cols, 2)
            mv()
        except:     
            mv_problems()

    # list of (path to) msas
    list_paths = list(Path(path_msas).rglob("*.[fa]*"))

    with ThreadPoolExecutor(max_workers=16) as pool:
        with tqdm(total=len(list_paths)) as progress:
            progress.set_description("Analyzing MSA")
            futures=[]
            for path_msa in list_paths:

                future = pool.submit(run, path_msa)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            
            for future in futures:
                future.result()


if __name__ == "__main__":

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("path_msas", help="directory to MSAs")
    parser.add_argument("path_output", help="output directory")
    args = parser.parse_args()

    main(args.path_msas, args.path_output)    