from Bio import AlignIO #load MSA

# TODO: add what is computed in src/eda_msas.py in this class
class AnalyzerMSA:

    def __init__(self, path_msa: str):
        self.path_msa = path_msa
        self.cols_msa = None # will be updated when the MSA is loaded
        self.seqs = self.get_seqs(path_msa)# sequences in the MSA

    def get_column(self, idx: int) -> list:
        return [seq[idx] for seq in self.seqs]

    def get_seqs(self,path_msa):
        # load MSA
        align=AlignIO.read(path_msa, "fasta")
        self.cols_msa = align.get_alignment_length()
        
        # extract sequences
        seqs = []
        for record in align:
            seqs.append(record.seq)
        return seqs

    def is_one_character(self, column) -> bool:
        chars = list(set(column)-set("-"))
        if len(chars)==1:
            return True
        return False

    def analyze_msa(self,) -> list[bool]:
        """evaluate each column and verify if it contains only one character 
        (indels are not taken into consideration)"""
        results = []
        for idx in range(self.cols_msa):
            column = self.get_column(idx)
            results.append(
                self.is_one_character(column)
            )
        return results