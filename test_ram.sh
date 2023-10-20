#!/bin/bash

# start=0
# standard_decomposition=True
# if [ $standard_decomposition = "True" ]; then
#     decomposition="standard"
# else    
#     decomposition="row-maximal"
# fi

# # Define a list of file paths
# file_paths=(
#     "/data/msas-pangeblocks/sars-cov-2/20-SARS-CoV2-MSA.fa"
#     "/data/msas-pangeblocks/sars-cov-2/50-SARS-CoV2-MSA.fa"
#     # "/data/msas-pangeblocks/sars-cov-2/100-SARS-CoV2-MSA.fa"
# )

# #  10 20 50 100 200 300 500 700 1000 1500 2000
# for path_msa in "${file_paths[@]}"
# do
#     name_ext="${path_msa##*/}"
#     name="${name_ext%.*}"
#     echo $name    

#     dirout="/data/analysis-paper/experiments/ram-usage/$name-$decomposition-decomp-not-alpha-consistent"
#     mkdir -p $dirout/logs

#     for end in 50 100 200 300 500 700 1000 1500 2000 2500
#     do 
#         echo $start,$end,$decomposition

#         /usr/bin/time --verbose src/exact_cover.py --path-msa $path_msa --obj-function nodes \
#             --prefix-output $dirout/$name \
#             --penalization 1000 --min-len 0 --min-coverage 0 \
#             --start-column $start --end-column $end --time-limit 10800 --solve-ilp True \
#             --use-wildpbwt True --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
#             --standard-decomposition $standard_decomposition --alpha-consistent False --threads-ilp 8 \
#             --workers 1 > $dirout/logs/$name-$decomposition-$start-$end.std.log 2> $dirout/logs/$name-$decomposition-$start-$end.err.log
#     done
# done

/usr/bin/time -v src/exact_cover.py --path-msa /data/msas-pangeblocks/sars-cov-2-clean/10-SARS-CoV2-MSA.fa \
    --obj-function nodes \
    --prefix-output test-ram \
    --penalization 0 --min-len 0 --min-coverage 0 \
    --start-column 0 --end-column 100 \
    --time-limit 180 --solve-ilp True  \
    --use-wildpbwt True --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
    --standard-decomposition False \
    --alpha-consistent False \
    --threads-ilp 8 --workers 1 2> test_ram.log
    