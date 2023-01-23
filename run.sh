# cat info-all-msas/analysis-msa/filtered_msas.txt | shuf --random-source=<(yes 42) -n 100 > 100msas.txt
# cat 100msas.txt | tr '\n' '\0' | xargs -0 -L1 -I '$' echo '$' > 100msas.txt
mkdir -p 100msas
cat 100msas.txt | xargs -I {} cp {} 100msas/