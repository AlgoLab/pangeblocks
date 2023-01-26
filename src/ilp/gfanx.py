"""
Postprocessing of the GFA using Networkx
Remove nodes with indels 
"""
# import sys
import networkx as nx
# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument("--path_gfa", help="path to save the output in GFA format", dest="path_gfa")
# args = parser.parse_args()

# path_gfa = args.path_gfa
def postprocessing(path_gfa):
    G = nx.DiGraph()

    torm = []

    for line in open(path_gfa):
        line = line.strip("\n").split("\t")
        if line[0] == "S":
            _, idx, seq = line
            seq = seq.replace("-","")
            G.add_node(idx, seq=seq)
            if seq == '':
                torm.append(idx)

    for line in open(path_gfa):
        line = line.strip("\n").split("\t")
        if line[0] == "L":
            _, idx1, _, idx2, _, _ = line
            G.add_edge(idx1, idx2)


    G.remove_nodes_from(torm)
    nodeattrs = nx.get_node_attributes(G, 'seq')
    for n in G.nodes:
        print('S', n, nodeattrs[n], sep='\t')


    new_paths = []
    for line in open(path_gfa):
        if line.startswith('P'):
            line = line.strip()
            d = line.split('\t')
            p = d[2][:-1].split('+,')
            p = [x for x in p if x not in torm]
            print(*d[:2], '+,'.join(p)+'+', '*', sep='\t')
            new_paths.append(p)

    # add missing edges
    for path in new_paths:
        for u,v in zip(path[:-1],path[1:]):
            if (u,v) not in G.edges:
                G.add_edge(u,v)

    for e in G.edges:
        print('L', e[0], '+', e[1], '+', '*', sep='\t')