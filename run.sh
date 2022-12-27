# pggb -i msas/GC00002971_r1_r1_1.fa.gz \  #input file in fasta format
# -o subsample-experiment/output-pggb \  # output directory  
# -n 103 \  # number of haplotypes          
# -s 10 \  # segment length for scaffolding the graph
# -t 16 \  # number of threads
# -p 10 \  # minimum average nucleotide identity for segments       
# pggb -i msas-pggb/GC00002971_r1_r1_1.fa.gz -p 70 -s 500 -n 10 -t 16 -o subsample-experiment/output-pggb
# NAME_MSA="GC00004575"
# python -i  compute_gfa.py --path_blocks /home/avila/pangeblocks/subsample-experiment/block_decomposition/$NAME_MSA.json \
#  --path_msa ilp-msas/$NAME_MSA.fa \
#  --path_gfa ilp-experiments/gfa/$NAME_MSA-2.gfa

# python compute_gfa.py --path_blocks ilp-experiment/block_decomposition/GC00002971_r1_r1_1.json \
# --path_msa msa-exp/GC00002971_r1_r1_1.fa \
# --path_gfa experiment/gfa/GC00002971_r1_r1_1.gfa

python src/check_sequence_in_msa.py --source 90 --sink 100 --path-msa msas-didelot/toyexample.fa --path-gfa experiment-didelot/gfa/toyexample.gfa
python src/check_sequence_in_msa.py --source 0 --sink 975 --path-msa msas-didelot/coli27-86.fa --path-gfa output-didelot/output-pandora/coli27-86.gfa