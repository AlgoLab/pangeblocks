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

# from src.ilp.input import InputBlockSet
# from src.maximal_blocks import compute_maximal_blocks

# ibs = InputBlockSet()
# path_msa = "/home/avila/pangeblocks/test/test3.fa"
# sc,ec=4,6
# maximal_blocks = compute_maximal_blocks(path_msa, start_column=sc,end_column=ec)
# decomp, missing, one_char = ibs(path_msa, maximal_blocks, sc, ec)

# src/solve_submsa.py --path-msa test/test3.fa -sc 4  -ec 6 --obj-function nodes --path-save-ilp output-test/ilp/test3_4-6.mps --path-opt-solution output-test/ilp/test3_4-6.json

# src/solve_submsa.py --path-msa test/test3.fa --submsa-index output-test/submsas/test3.txt --obj-function nodes --path-save-ilp output-test/ilp/test3 --path-opt-solution output-test/ilp/test3

# /home/avila/Wild-pBWT/bin/wild-pbwt -a 5 -f test-gfa/haplotypes.txt -o y
# > haplotypes_maxblocks.txt
# git clone https://github.com/illoxian/Wild-pBWT.git


# fil-profile --no-browser run submsas.py \
# --path-vertical-blocks output-test-grid/maximal-blocks/test/vertical_blocks_alpha1.json\
# --path-msa test/test.fa \
# --output out.txt


# ENVIRONMENT
# mamba activate pangeblocks

# fil-profile --no-browser run src/solve_submsa.py --path-msa /data/alessia/covid/big-msa/500-SARS-CoV2-MSA.fasta \
# -sc 100 -ec 200 \
# --obj-function nodes \
# --path-save-ilp /data/alessia/covid/output-big-msa/ilp/500-SARS-CoV2-MSA/alpha1/500-SARS-CoV2-MSA \
# --path-opt-solution /data/alessia/covid/output-big-msa/ilp/500-SARS-CoV2-MSA/alpha1/500-SARS-CoV2-MSA \
# --penalization 0 \
# --min-len 0 \
# --min-coverage 0 \
# --time-limit 30 \
# --solve-ilp true \
# --use-wildpbwt True \
# --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
# --workers 16 2> profile-100-200-nobug.log

# --submsa-index /data/alessia/covid/output-big-msa/submsas/500-SARS-CoV2-MSA_alpha1.txt \

python src/solve_submsa.py --path-msa test/test.fa -sc 5 -ec 8 --path-save-ilp test-opt --path-opt-solution test-opt 2> test.err.log > test.stdout.log 