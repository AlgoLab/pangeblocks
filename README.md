<img src="img/logo-pangeblocks.png" width="300" height="200">

**Customized construction of pangenome graphs via maximal blocks**

## Pangenome graphs
___
A pangenome graph is a data structure for a set of genomic sequences, where nodes are labeled by substrings in the alphabet $\{A,C,G,T\}$ in the case of DNA (and $\{A,C,G,U\}$ in the case of RNA). Input sequences are represented by paths in the graph, where each input sequence is spell by the concatenation of labels of the nodes in its path. 

## What is a block?
___
In the context of an MSA, a **block** is defined as set of rows
that shares the same string in an interval of columns.

Formally, a block composed of the set of rows $K$, in the interval $[b,e]$ 
is defined as the triplet $(K, b, e)$. 

In the image below, each block is highlighted by a different color. For example, 

- The yellow one in the middle: $(\{\square, \lozenge, \triangle \}, 3, 5)$  spelling **A - -**
- The blue one spanning all rows:  $(\{ \circ ,\square, \lozenge, \triangle, \triangledown \}, 7, 8)$  spelling **G A**
- The green one in column  5: $(\{\circ,\triangledown \}, 5, 5)$ spelling spelling the single character **C**


## How it works?
___

`pangeblocks` creates a variation graph from an MSA by selecting a set of blocks. 
It creates a search space of blocks from **maximal blocks**, and then an Integer Linear Programming model selects the best subset of blocks to cover all cells of the MSA.

<img src="img/matrix-cover-style.svg" width="400" height="200">

- For each block we create a node with its label.
- Consecutive blocks are connected by an arc.
- Each input sequence in the MSA is spelled by a path in the graph 

<img src="img/variation-graph.svg" width="400" height="200">

Finally, indels are removed from the graph (could be the remotion of an entire node, or the removal of indels in the label of a node), and non-branching paths are collapsed. 

<img src="img/variation-graph-postprocessed.svg" width="400" height="200">


### Create a virtual environment
```bash
mamba create env -n pangeblocks -f envs/snakemake.yml
mamba activate pangeblocks
```

### Create Variation Graphs from MSAs

To construct variation graphs from MSAs, run `pangeblocks`:
```bash
snakemake -s pangeblock.smk -c32 --use-conda # variation graph as GFA
```

### Parameters

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


___

## Considerations


- **ALPHABET** $\{A,C,G,T,-,N\}$ **Not** case sensitive. We recommend to map all characters not in the alphabet to N.


___
[OPTIONAL]
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
```
The above smk pipeline will analyze the MSAs and output two files in `PATH_OUTPUT/analysis-msas`:
1. `stats_msas.tsv` with basic information about the MSAS: path, number of columns and rows (sequences), number of identical columns, and number of unique sequences
2. `problematic_msas.tsv`: contains a list of MSAs that has no information
