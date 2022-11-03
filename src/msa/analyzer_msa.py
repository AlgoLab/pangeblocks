from Bio import AlignIO #load MSA
from typing import Union
from pathlib import Path

# TODO: add what is computed in src/eda_msas.py in this class
class AnalyzerMSA:

    # def __init__(self,):
    #     self.n_cols = None # number of columns
    #     self.n_seqs = None # number of sequences

    def __call__(self, path_msa: Union[str, Path]) -> tuple:
        "Compute metrics/analysis over an MSA"
        # load MSA
        align, n_seqs, n_cols = self.load_msa(path_msa)
        
        # sequences and number of columns in the MSA
        seqs = self.get_seqs(align)

        # unique sequences
        n_unique_seqs = len(set([str(seq) for seq in seqs]))

        # identical columns
        identical_cols = self.check_identical_columns(seqs, n_cols)
        n_identical_cols = sum(identical_cols)
        return (
                n_cols, 
                n_seqs, 
                n_unique_seqs, 
                n_identical_cols
                )

    def get_column(self, idx: int, seqs: list) -> list:
        return [seq[idx] for seq in seqs]

    def get_seqs(self,align):
        "get sequences from an alignment"
        # extract sequences
        seqs = []
        for record in align:
            seqs.append(record.seq)
        return seqs

    def is_one_character(self, column) -> bool:
        "Return True if a column (list) contain only one character ('-' are ommited)"
        chars = list(set(column)-set("-"))
        if len(chars)==1:
            return True
        return False

    def check_identical_columns(self, seqs, n_cols) -> list[bool]:
        """evaluate each column and verify if it contains only one character 
        (indels are not taken into consideration)"""
        results = []
        for idx in range(n_cols):
            column = self.get_column(idx, seqs)
            results.append(
                self.is_one_character(column)
            )
        return results

    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols