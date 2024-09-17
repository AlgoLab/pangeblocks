<img src="img/logo-pangeblocks-no-background.png" width="300" height="200">

# Customized construction of pangenome graphs via maximal blocks

## How to run

### Prerequisites
Clone the repo and setup the conda environment:
```bash
git clone https://github.com/AlgoLab/pangeblocks.git
cd pangeblocks
mamba env create -n pangeblocks -f envs/snakemake.yml
conda activate pangeblocks
```

#### Gurobi
`pangeblocks` requires gurobi ([check license](https://www.gurobi.com/)).

#### Example 
To run `pangeblocks` on the example data we provide (MSA with 10 rows and 200 columns), run:
```bash
./pangeblocks --dir-msa test/sars-cov-2-subMSA --dir-output output-sars-cov-2 --obj-function weighted --penalization 100 --min-len 20
```
You will find the final graph in `output-sars-cov-2/gfa-unchop`.

Alternatively, you can run `pangeblocks` on the small example using docker:
```
mkdir /tmp/pgb-out
docker run -it --user $(id -u):$(id -g) -v ./test/sars-cov-2-subMSA/:/data \
    --mount type=bind,source=/tmp/pgb-out,target=/results algolab/pangeblocks:latest
ls /tmp/pgb-out/sars-cov-2.gfa # <- this is the graph
```

### Tweak the execution

`pangeblocks` is a snakemake pipeline. We provide a CLI that parses command line options and exectues the snakemake pipeline to compute graphs from a directory with <msa>.fa files.

```bash
usage: pangeblocks [-h] [--dir-msa DIR_MSA] [--dir-output DIR_OUTPUT] [--log-level LOG_LEVEL] [--obj-function {nodes,strings,weighted,depth,depth_and_len}] [--penalization PENALIZATION]
                   [--min-len MIN_LEN] [--time-limit TIME_LIMIT] [--threshold-vertical-blocks ALPHA] [--min-coverage MIN_COVERAGE] [--larger-decomposition] [--consistent]
                   [--cores THREADS] [--submsa_threads SUBMSA_THREADS] [--ilp-threads ILP_THREADS] [--max-memory MAX_MEMORY] [--min-rows-block MIN_ROWS_BLOCK]
                   [--max-rows-block MAX_ROWS_BLOCK] [--max-msa-size MAX_MSA_SIZE]

options:
  -h, --help            show this help message and exit
  --dir-msa DIR_MSA     directory to MSAs in .fa format
  --dir-output DIR_OUTPUT
                        path to save the outputs in GFA format
  --log-level LOG_LEVEL
                        set log level (ERROR/WARNING/INFO/DEBUG)
  --obj-function {nodes,strings,weighted,depth,depth_and_len}
                        the objective function to optimize
  --penalization PENALIZATION
                        used only with the 'weighted', 'depth', and 'depth_and_len' obj function
  --min-len MIN_LEN     used only with the 'weighted' obj function
  --time-limit TIME_LIMIT
                        Timeout (in minutes) to stop ILP
  --threshold-vertical-blocks ALPHA
                        minimum width of a vertical block
  --min-coverage MIN_COVERAGE
                        used only with 'depth' obj function
  --larger-decomposition
                        if True, use complete-decomposition of blocks, otherwise use row-maximal decomposition
  --consistent          use an alpha-consistent strategy
  --cores THREADS       Number of cores to be used
  --submsa_threads SUBMSA_THREADS
  --ilp-threads ILP_THREADS
  --max-memory MAX_MEMORY
                        Maximum RAM used (in MBytes)
  --min-rows-block MIN_ROWS_BLOCK
  --max-rows-block MAX_ROWS_BLOCK
  --max-msa-size MAX_MSA_SIZE
```

### Run a grid experiment
___

To construct variation graphs from MSAs, **set parameters** in `params.yml` and then run `pangeblocks` 
```bash
snakemake -s pangeblocks.smk -c16 --use-conda # variation graph as GFA
```
for each file in the directory defined at `PATH_INPUT` a .gfa file will be created in `PATH_OUTPUT/gfa-unchop`

### Running under docker
___

If you want to run  `pangeblocks` on your data using docker, you have to provide the
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

___

## Considerations

- Maximal blocks are computed with the [Wild-PBWT](https://github.com/AlgoLab/Wild-pBWT)
- **Troubleshooting** Wild-PBWT requires SDSL, [check this to install it](https://github.com/msgr0/Wild-pBWT?tab=readme-ov-file#prerequisites)
- ILPs are solved using [Gurobi](https://www.gurobi.com/), you might need a license 
- Each MSA must be in the **ALPHABET** $\{A,C,G,T,-,N\}$ **Not** case sensitive. We recommend to map all characters not in the alphabet to N.

___

## How does it work?
___
`pangeblocks` creates a variation graph from an MSA by selecting a set of blocks. 
It creates a search space of blocks from **maximal blocks**, and then an Integer Linear Programming model selects the best subset of blocks to cover all cells of the MSA.

<img src="img/matrix-cover-style.svg" width="600" height="300">

- For each block we create a node with its label.
- Consecutive blocks are connected by an arc.
- Each input sequence in the MSA is spelled by a path in the graph 

<img src="img/variation-graph.svg" width="600" height="300">

Finally, indels are removed from the graph (could be the remotion of an entire node, or the removal of indels in the label of a node), and non-branching paths are collapsed. 

<img src="img/variation-graph-postprocessed.svg" width="600" height="400">

## Supplementary tools
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
```
The above smk pipeline will analyze the MSAs and output two files in `PATH_OUTPUT/analysis-msas`:
1. `stats_msas.tsv` with basic information about the MSAS: path, number of columns and rows (sequences), number of identical columns, and number of unique sequences
2. `problematic_msas.tsv`: contains a list of MSAs that has no information
