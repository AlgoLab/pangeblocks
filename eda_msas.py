import time 
from pathlib import Path
from Bio import AlignIO #load MSA
from src.utils.monitor_values import MonitorValues
from rich.progress import track

# log values
mv_msas = MonitorValues(["path_msa","n_seqs","n_unique_seqs","n_cols"])#,"n_max_blocks","t_max_blocks"])
mv_msas_problems = MonitorValues(["path_msa"])
mv_blocks = MonitorValues(["path_msa","n_seqs","n_cols","n_max_blocks","t_max_blocks"])

list_msa=list(Path("/home/disco/Data/pandora-msas/msas/").rglob("*fa"))#GC00000173_5.fa"
# list_msa = ["/home/disco/Data/pandora-msas/msas/GC00000708_59.fa"]

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
df.to_csv("stats_msas.tsv",sep="\t")
mv_msas_problems.get_values_asdf().to_csv("problem_msas.tsv", sep="\t")