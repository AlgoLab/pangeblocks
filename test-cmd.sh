#!/bin/sh

#/usr/bin/time --verbose python src/compute_gfa.py  --path_blocks test/slow/block_decomposition/DMA-3108.json --path_msa test/slow/msas-HLA-zoo/DMA-3108.fa --path_gfa test/slow/penalization7-min_len5/DMA-3108.gfa --path_oc test/slow/penalization7-min_len5/DMA-3108.json --obj_function weighted --penalization 7 --min_len 5 --time_limit 10 --log_level INFO  > out.log 2> err.log
#scp pangeblocks:/data/pangeblocks-experiments/output-HLA-zoo-mafft.op3-ep0/block_decomposition/DMA-3108.json test/slow/block_decomposition/slow.json
#scp pangeblocks:/data/msas-pangeblocks/msas-HLA-zoo/mafft.op3-ep0/DMA-3108.fa test/slow/msas-HLA-zoo/slow.fa
#scp pangeblocks:/data/pangeblocks-experiments/output-HLA-zoo-mafft.op3-ep0/opt-coverage/weighted/penalization3-min_len3/DMA-3108.json test/slow/penalization7-min_len5/slow.json

# /usr/bin/time --verbose python src/compute_gfa.py  --path_blocks test/slow/block_decomposition/slow.json --path_msa test/slow/msas-HLA-zoo/slow.fa --path_gfa test/slow/penalization7-min_len5/slow.gfa --path_oc test/slow/penalization7-min_len5/slow.json --obj_function weighted --penalization 7 --min_len 5 --time_limit 10 --log_level INFO  > out.log 2> err.log

# op=$1
# ep=0
# path_msas="/data/msas-pangeblocks/msas-tests"
# path_output="output"
# echo $path_msas
# echo "waiting 7200 secs"
# sleep 7200
# snakemake -s eda.smk -c8 
# snakemake -s pangeblock-grid-exp.smk -c16 --use-conda
# --config PATH_MSAS=$path_msas --config PATH_OUTPUT=$path_output

from src.ilp.input import InputBlockSet
from src.maximal_blocks import compute_maximal_blocks

ibs = InputBlockSet()
path_msa = "/home/avila/pangeblocks/test/test3.fa"
sc,ec=4,6
maximal_blocks = compute_maximal_blocks(path_msa, start_column=sc,end_column=ec)
decomp, missing, one_char = ibs(path_msa, maximal_blocks, sc, ec)

src/solve_submsa.py --path-msa test/test3.fa -sc 4  -ec 6 --obj-function nodes --path-save-ilp output-test/ilp/test3_4-6.mps --path-opt-solution output-test/ilp/test3_4-6.json

src/solve_submsa.py --path-msa test/test3.fa --submsa-index output-test/submsas/test3.txt --obj-function nodes --path-save-ilp output-test/ilp/test3 --path-opt-solution output-test/ilp/test3