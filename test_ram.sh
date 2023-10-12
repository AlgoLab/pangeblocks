#!/bin/bash
# path_submsas_index="output-test/submsas/test3_alpha50.txt" # "/data/pangeblocks-experiments/sars-cov-2/submsas/100-SARS-CoV2-MSA_alpha50.txt"  
path_msa="/data/msas-pangeblocks/sars-cov-2/100-SARS-CoV2-MSA.fasta" # "test/test3.fa"
name="100-SARS-CoV2-MSA"

# cat $path_submsas_index |\
# while read start end; do 
#     echo "start: $start, end: $end";
#     # run exact cover in the current submsa

#     /usr/bin/time --verbose src/exact_cover.py --path-msa $path_msa --obj-function nodes \
#         --prefix-output {params.dir_subsols}/{wildcards.name_msa} \
#         --penalization 1000 --min-len 0 --min-coverage 0 \
#         --start-column $start --end-column $end --time-limit {params.time_limit} --solve-ilp True \
#         --use-wildpbwt True --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
#         --standard-decomposition False --threads-ilp {params.threads_ilp} \
#         --workers {params.workers} > {output.auxfile} 2> {log.stderr}

# done

# TODO: identify other subMSAs that are long enough, or where the previous implementation crushed and try them

# start=18957
# end=21681
start=0
# end=400

dirout="output-ram-standard_decomp"
mkdir -p $dirout/logs

#  10 20 
for end in 50 100 200 300 500 700 1000 1500 2000 2500
do 
    echo $start,$end
    # /usr/bin/time --verbose src/exact_cover.py --path-msa $path_msa --obj-function nodes \
    #     --prefix-output $dirout/opt-$name \
    #     --penalization 1000 --min-len 0 --min-coverage 0 \
    #     --start-column $start --end-column $end --time-limit 1800 --solve-ilp True \
    #     --use-wildpbwt True --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
    #     --standard-decomposition False --threads-ilp 8 \
    #     --workers 1 > $dirout/logs/$name-$start-$end.std.log 2> $dirout/logs/$name-$start-$end.err.log

    /usr/bin/time --verbose src/exact_cover_notU.py --path-msa $path_msa --obj-function nodes \
        --prefix-output $dirout/opt-$name-notU \
        --penalization 1000 --min-len 0 --min-coverage 0 \
        --start-column $start --end-column $end --time-limit 1800 --solve-ilp True \
        --use-wildpbwt True --bin-wildpbwt Wild-pBWT/bin/wild-pbwt \
        --standard-decomposition True --threads-ilp 8 \
        --workers 1 > $dirout/logs/$name-notU-$start-$end.std.log 2> $dirout/logs/$name-notU-$start-$end.err.log
done