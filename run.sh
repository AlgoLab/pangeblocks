# cat info-all-msas/analysis-msa/filtered_msas.txt | shuf --random-source=<(yes 42) -n 100 > 100msas.txt
# cat 100msas.txt | tr '\n' '\0' | xargs -0 -L1 -I '$' echo '$' > 100msas.txt
# mkdir -p 100msas
# cat 100msas.txt | xargs -I {} cp {} 100msas/

path_gfa="output-msa-difficile-weighted/gfa/slpa-basis.mafft.gfa" 

# Chop all nodes to 1 base, - are replaced by N
../vg mod -X 1 $path_gfa > slpa-basis.mafft.chop1.gfa

# Load to networkx to rm nodes with N, and let it do the hard work of nodes removal
python src/ilp/gfanx.py slpa-basis.mafft.chop1.gfa > slpa-basis.mafft.chop1.noN.gfa

# Unchop nodes
../vg mod -u slpa-basis.mafft.chop1.noN.gfa > slpa-basis.mafft.unchop.noN.gfa

# Generate bandaga labels
python src/graph/bandage_labels_from_gfa.py --path_gfa slpa-basis.mafft.unchop.noN.gfa --path_save slpa-basis.mafft.unchop.noN.csv

python src/graph/bandage_labels_from_gfa.py --path_gfa slpa-basis.mafft.chop1.gfa --path_save slpa-basis.mafft.chop1.csv

python src/graph/bandage_labels_from_gfa.py --path_gfa slpa-basis.mafft.chop1.noN.gfa --path_save slpa-basis.mafft.chop1.noN.csv
