import sys
import time 

from pathlib import Path
from Bio import AlignIO #load MSA
from utils.monitor_values_plus import MonitorValuesPlus as MonitorValues
from rich.progress import track

def run(path_msas: str, prefix_output: str):
    dir_output = Path(prefix_output).parent
    dir_output.mkdir(exist_ok=True, parents=True)

    out_file = Path(prefix_output + "-stats_msa.tsv")
    # out_file.parent.mkdir(exist_ok=True, parents=True)
    problematic_out_file = Path(prefix_output + "-problematic_files.tsv")

    # log values
    mv_msas = MonitorValues(["path_msa","n_seqs","n_unique_seqs","n_cols"], out_file)
    mv_msas_problems = MonitorValues(["path_msa"], problematic_out_file)

    list_msa=list(Path(path_msas).glob("*.[fa]*"))
    for path_msa in track(list_msa, description="EDA MSAs"):

        try:
            # load MSA, count seqs and columns
            align=AlignIO.read(path_msa, "fasta")
            n_cols = align.get_alignment_length()
            n_seqs = len(align)
            n_unique_seqs = len(set([str(record.seq) for record in align]))

            mv_msas()
        except:
            mv_msas_problems()

    # df = mv_msas.get_values_asdf()
    # df.sort_values(by=["n_cols","n_seqs"], inplace=True)
    # df.to_csv(out_file,sep="\t")
    # mv_msas_problems.get_values_asdf().to_csv(problematic_out_file, sep="\t")

if __name__ == "__main__":

    # TODO: use argparse
    path_msas = sys.argv[-1]
    run(path_msas)