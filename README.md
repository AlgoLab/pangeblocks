# Pangeblocks (work in progress)
Pangenome graph construction from maximal blocks in an MSA 

Set the parameters in `params.yml`:
```yaml
PATH_MSAS: "msas" # folder containing MSAs in .fasta/.fa format
PATH_OUTPUT: "output" # folder where to save the results
OPTIMIZATION:
  OBJECTIVE_FUNCTION:
    - "nodes"    # minimize number of nodes 
    - "strings"  # minimize length of the graph
    - "weighted" # penalize by PENALIZATION blocks shorter that MIN_LEN (other blocks cost=1)
    - "depth"    # penalize by PENALIZATION blocks used by less than MIN_COVERAGE (other blocks cost=1)
  PENALIZATION: # used only with "weighted" and "depth"
    - 3
    - 5
    - 7 
    - 128
  MIN_LEN: # used only with "weighted"
    - 3
    - 5
    - 10   
  MIN_COVERAGE: # used only with "depth"
    - 0.1
    - 0.2
    - 0.3
  THRESHOLD_VERTICAL_BLOCKS: # minimum length for vertical blocks to be fixed in the optimal solution
    - 1
    - 2
    - 8
    - 16
  TIME_LIMIT: 30 # time limit to run each ILP (minutes)
LOG_LEVEL: "INFO"
THREADS: 
  TOTAL: 32
  SUBMSAS: 16
  ILP: 8
```
___

### Create a virtual environment and install dependencies
```bash
python3 -m venv .pbenv
source .pbenv/bin/activate
pip install -r requirements.txt
```

or
```bash
mamba create env -n pangeblocks -f envs/pangeblocks.yml
mamba activate pangeblocks
```

### Create Variation Graphs from MSAs

[PANGEBLOCKS]
To construct variation graphs from MSAs, run `PangeBlocks``:
```bash
snakemake -s pangeblock.smk -c32 --use-conda # variation graph as GFA
```
___
[OPTIONAL]
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
```
The above smk pipeline will analyze the MSAs and output two files in `PATH_OUTPUT/analysis-msas`:
1. `stats_msas.tsv` with basic information about the MSAS: path, number of columns and rows (sequences), number of identical columns, and number of unique sequences
2. `problematic_msas.tsv`: contains a list of MSAs that has no information

### Running under docker

To run `pangeblocks` on a small example (MSA with 10 rows and 200 columns), run
the following command, replacing `/tmp/pgb` with the directory that will contain
the results.

```
docker run -it --user $(id -u):$(id -g) \ 
-v ./test/sars-cov-2-subMSA/:/data \
--mount type=bind,source=/tmp/pgb,target=/results \
algolab/pangeblocks:latest
```

If you want to run  `pangeblocks` on your data, you also have to provide the
directory containing the MSA, replacing `./test/sars-cov-2-subMSA/` with your
directory and adding the correct  `pangeblocks` call, specifying the arguments.
For example:

```
docker run -it --user $(id -u):$(id -g) \ 
-v ./DATADIR/:/data \
--mount type=bind,source=/tmp/pgb,target=/results \
algolab/pangeblocks:latest
/app/pangeblocks --path-msa /data/my.msa
```
