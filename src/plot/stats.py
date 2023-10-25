# import argparse
from collections import Counter

def dist_len(path_submsas):

    lens = []
    with open(path_submsas) as fp:
        for line in fp.readlines():
            start, end = line.replace("\n","").split("\t")
            start, end = int(start), int(end)
            length = end - start + 1
            lens.append(length)
        
    return Counter(lens)

c = dist_len("/data/alessia/covid/output-big-msa-2/submsas/1000-SARS-CoV2-MSA_alpha1.txt")
