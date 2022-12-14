# Pangeblocks (work in progress)
Pangenome graph construction from maximal blocks in an MSAs 

Set the parameters in `params.yaml`:
```yaml
PATH_MSAS: "msas"          # folder containing MSAs in fasta format
PATH_OUTPUT: "experiment"  # folder where to save the results
```
___

### Create a virtual environment and install dependencies
```bash
python3 -m venv pbenv
source pbenv/bin/activate
pip install -r requirements.txt
```

### Run pipelines
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
```
The above smk pipeline will analyze the MSAs and output two files in `PATH_OUTPUT/analysis-msas`:
1. `stats_msas.tsv` with basic information about the MSAS: path, number of columns and rows (sequences), number of identical columns, and number of unique sequences
2. `problematic_msas.tsv`: contains a list of MSAs that has no information

After the previous pipeline has run, the computation of maximal blocks will be done only in the MSAs in `stats_msas.tsv` 
```bash
snakemake -s compute_blocks.smk -c16 # compute max blocks for MSAs
```