import sys
import time 

from pathlib import Path
from Bio import AlignIO #load MSA
from utils.monitor_values import MonitorValues
from rich.progress import track

def run(path_msas: str):
    
    out_file = Path("out/stats_msas.tsv")
    out_file.parent.mkdir(exist_ok=True, parents=True)
    problematic_out_file = Path("out/problematic_files.tsv")

    # log values
    mv_msas = MonitorValues(["path_msa","n_seqs","n_unique_seqs","n_cols"])
    mv_msas_problems = MonitorValues(["path_msa"])

    list_msa=list(Path(path_msas).rglob("*fa"))
    for path_msa in track(list_msa, description="Working on MSAs"):

        try:
            # load MSA, count seqs and columns
            align=AlignIO.read(path_msa, "fasta")
            n_cols = align.get_alignment_length()
            n_seqs = len(align)
            n_unique_seqs = len(set([str(record.seq) for record in align]))

            mv_msas()
        except:
            mv_msas_problems()

    df = mv_msas.get_values_asdf()
    df.sort_values(by=["n_cols","n_seqs"], inplace=True)
    df.to_csv(out_file,sep="\t")
    mv_msas_problems.get_values_asdf().to_csv(problematic_out_file, sep="\t")

if __name__ == "__main__":

    # TODO: use argparse
    path_msas = sys.argv[-1]
    run(path_msas)