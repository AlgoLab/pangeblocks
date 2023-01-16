ls ../data/msas | shuf --random-source=<(yes 42) -n 100 > ../data/100msas.txt
cat ../data/100msas.txt | tr '\n' '\0' | xargs -0 -L1 -I '$' echo '../data/msas/$' > 100msas.txt
mkdir -p 100msas
cat 100msas.txt | xargs -I {} cp {} 100msas/