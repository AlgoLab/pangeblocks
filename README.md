# Pangeblocks (work in progress)
Pangenome graph construction from maximal blocks in an MSAs 

Set the parameters in `params-grid-exp.yaml`:
```yaml
PATH_MSAS: "msas" # folder containing MSAs in .fasta/.fa format
PATH_OUTPUT: "output-pangeblocks" # folder where to save the results
OPTIMIZATION:
  OBJECTIVE_FUNCTION: "strings" # one of: nodes, strings, weighted
  # used only with "weighted"
  PENALIZATION: # for blocks shorter than MIN_LEN: cost=PENALIZATION in the objective function, otherwise cost=1 
    - 3 
  MIN_LEN: 
    - 2 
  TIME_LIMIT: 30 # time limit to run the ILP (minutes)
LOG_LEVEL: "INFO"
```
___

### Create a virtual environment and install dependencies
```bash
python3 -m venv .pbenv
source .pbenv/bin/activate
pip install -r requirements.txt
```

### Run pipelines
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
```
The above smk pipeline will analyze the MSAs and output two files in `PATH_OUTPUT/analysis-msas`:
1. `stats_msas.tsv` with basic information about the MSAS: path, number of columns and rows (sequences), number of identical columns, and number of unique sequences
2. `problematic_msas.tsv`: contains a list of MSAs that has no information

After the previous pipeline has run, the computation of pangeblock graphs will be done only in the MSAs in `stats_msas.tsv` 
```bash
snakemake -s pangeblock-grid-exp.smk -c16 # variation graph as GFA
```