# python src/check_sequence_in_msa.py --source 90 --sink 100 --path-msa msas-didelot/toyexample.fa --path-gfa experiment-didelot/gfa/toyexample.gfa
# python src/check_sequence_in_msa.py --source 0 --sink 975 --path-msa msas-didelot/coli27-86.fa --path-gfa output-didelot/output-pandora/coli27-86.gfa

# python compute_gfa.py --log_level=DEBUG --path_blocks output-msas/block_decomposition/Cluster_12313.json \
# --path_msa msas/Cluster_12313.fa --path_gfa output-msas/gfa/Cluster_12313.gfa --path_ilp output-msas/tmp/ilp.lp 2> output-msas/tmp/log.txt

# python compute_gfa.py --log_level=DEBUG --path_blocks output-msas/block_decomposition/GC00004574.json \
# --path_msa msas/GC00004574.fa --path_gfa output-msas/gfa/GC00004574.gfa --path_ilp output-msas/tmp/ilp.lp 2> output-msas/tmp/log.txt

python compute_gfa.py --log_level=DEBUG --path_blocks output-test/block_decomposition/test.json \
--path_msa test/test.fa --path_gfa output-test/gfa/test.gfa --path_ilp output-test/tmp/ilp.lp 2> output-test/tmp/log.txt
