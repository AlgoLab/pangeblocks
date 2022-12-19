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

After the previous pipeline has run, the computation of pangeblock graphs will be done only in the MSAs in `stats_msas.tsv` 
```bash
snakemake -s pangeblock.smk -c16 # variation graph as GFA
```

___
### Independent pipelines
`make_prg.smk` will download version as indicated here https://github.com/iqbal-lab-org/make_prg#download (time accession: 19/12/2022)
, run make_prg and extract the gfa outputs for the MSAs that are in `PATH_MSAS`
