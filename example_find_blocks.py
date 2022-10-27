from tkinter import FIRST
import numpy as np
import networkx as nx

from collections import defaultdict
from graphviz import Digraph
from pprint import pprint 
from networkx.algorithms.flow import maximum_flow

from pangeblocks.blocks_suffix_tree import (
    Block,
    map_coverage, 
    maximal_blocks_from, 
    is_maximal_block, 
)

from pangeblocks.graph_max_flow import (
    nodes_edges_from_blocks, 
    is_consecutive,
    missing_edge,
    Node, 
    Edge, 
)

## Inputs
ALPHABET = {"A","C","G","T"}
FIRST_CHAR = chr(1000)

seqs = [
        "ACCGTCGTAAAATT--ACG",
        "ACCGTCGTCCCCTT--ACG",
        "ACCGTCGTCCCCTTACACT"  
    ]   

num_seqs = len(seqs)
id_seqs = list(range(num_seqs))
# f = lambda n : "a" if n=="0" else "b"
# seqs = [ "".join([f(n) for n in seq]) for seq in seqs]

block2node = lambda block: Node(block.K, block.i, block.j, block.seq)
node2str = lambda node: f"{'-'.join(str(_) for _ in node.K)},{node.i},{node.j},{node.seq}"

lens = [len(seq) for seq in seqs]
size_msa = max(lens)
assert all([len_seq==size_msa for len_seq in lens]), "all seqs must have same length"

## Experiment
# 0. source and sink blocks
source_block = Block(K=id_seqs, i=-2,j=-2,seq='source')
post_blocks = [Block(K=[id_seq], i=-1,j=-1, seq=f'post-source-{id_seq}') for id_seq,_ in enumerate(seqs)]

prev_blocks = [Block(K=[id_seq], i=size_msa,j=size_msa, seq=f'prev-sink-{id_seq}') for id_seq,_ in enumerate(seqs)]
sink_block = Block(K=id_seqs, i=size_msa+1,j=size_msa+1,seq='sink')

# 1. Find maximal blocks
# following "Linear‑time method I: based on suffix trees" from https://doi.org/10.1186/s13015-020-0163-6
max_blocks = []

# iterate over all starting positions
for i in range(size_msa-1):
    blocks = maximal_blocks_from(i, seqs = seqs)
    # add new maximal blocks    
    for block in blocks:
        if is_maximal_block(block, max_blocks):
            max_blocks.append(block)

# --- 
# 2. Missing non-maximal blocks needed to cover the MSA
# what is covered with the max blocks
coverage = map_coverage(max_blocks, size_msa = size_msa , n_seqs=len(seqs))

# add blocks of |K| = 1 (coverage 0)
no_cov = np.where(coverage==0)
list_seq, list_col = no_cov
list_seq = np.append(list_seq, -1)
list_col  = np.append(list_col, -1)

# find positions of blocks
pos_extra_blocks = defaultdict(list)
cols_block = [list_col[0]]
for idx1, idx2 in zip(range(len(list_col)-1), range(1, len(list_col))):
    pair1 = list_seq[idx1], list_col[idx1]
    pair2 = list_seq[idx2], list_col[idx2]
    
    if is_consecutive(pair1, pair2):
        cols_block.append(pair2[1])
    else:
        pos_extra_blocks[pair1[0]].append(cols_block)
        cols_block = [pair2[1]]
    # print(pair1 + pair2 + is_consecutive(pair1,pair2))    
    pair1 = pair2

extra_blocks = []
# Create extra blocks
for id_seq, list_cols_blocks in pos_extra_blocks.items():
    for cols in list_cols_blocks:  
        i = min(cols)
        j = max(cols)
        seq = seqs[id_seq][i:j+1]
        block = Block([id_seq], i, j, seq)
        print(block)
        extra_blocks.append(block)

# sort max_blocks by starting position (to intersect and create nodes)
all_blocks = []
all_blocks.extend([source_block])
all_blocks.extend(post_blocks)
all_blocks.extend(prev_blocks)
all_blocks.extend([sink_block])
all_blocks.extend(max_blocks)
all_blocks.extend(extra_blocks)
all_blocks = sorted(all_blocks, key=lambda block: block.i)

# ---   
# Graph
list_nodes = []
list_edges = []

for j, block1 in enumerate(all_blocks):    

    for block2 in all_blocks:
        if block1.seq!="source" and block2.seq != "sink":
            nodes, edges = nodes_edges_from_blocks(block1, block2)

            list_nodes.extend(nodes)
            list_edges.extend(edges)

# agregar las aristas desde el nodo fuente
source_node = block2node(source_block)
list_nodes.append(source_node)
for block in post_blocks:
    node = block2node(block)
    capacity = 1
    list_nodes.append(node)
    list_edges.append(Edge(source_block, node, capacity))

# connect each post-block-node to nodes that contains the id_seq in .K
for block1 in post_blocks:
    for block2 in all_blocks:
        if block2.i==0 and set(block1.K).intersection(block2.K):
            node1 = block2node(block1)
            node2 = block2node(block2)
            capacity = 1 #len(set(source_block.K).intersection(set(block.K)))
            list_nodes.extend([node1,node2])
            list_edges.append(Edge(node1, node2, capacity))

# agregar las aristas hasta el nodo sink
sink_node = block2node(sink_block)
list_nodes.append(sink_node)
for block in prev_blocks:
    node = block2node(block)
    capacity = 1
    list_nodes.append(node)
    list_edges.append(Edge(node, sink_node, capacity))

# Add missing edges between consecutives block-nodes
# consecutive blocks : connect them

list_missing_edges = []
for node1 in list_nodes:
    for node2 in list_nodes:
        if node1 == block2node(source_block):
            print(node2)
        edge = missing_edge(node1, node2)
        list_missing_edges.extend(edge)
        list_edges.extend(edge)

# _____________________________________________________________       
# Create graph and save as dot
from graphviz import Digraph
pg = Digraph(format="png")

hash_node = dict()
id_node = 0
for node in list_nodes:
    hn = hash_node.get(node2str(node), False)
    if hn is False:
        hash_node[node2str(node)] = (str(id_node), node)
        id_node +=1

for str_node, info_node in hash_node.items():
    id_node, node = info_node
    pg.node(id_node, str_node) #f"{node.i},{node.j},{node.seq}")

# list_edges = list(set(list_edges))
edges_added = []
capacities = dict()
for edge in list_edges:
    hash_edge = (
        hash_node[node2str(edge.node1)][0] , 
        hash_node[node2str(edge.node2)][0]
        )
    if hash_edge not in edges_added:
        pg.edge(
            hash_edge[0], hash_edge[1], label=str(edge.capacity)
            )
        edges_added.append(hash_edge)
        capacities[hash_edge] = edge.capacity

print(pg.source)
pg.render(directory='output/graph', view=True).replace('\\', '/')

# ________________________
# Build pangeblock_
# edges in networkx format
idnode2node = {idnode: info_node for idnode, info_node in hash_node.values()}
edges_nx = []
for edge in edges_added:
    id_node1, id_node2 = edge
    capacity = capacities[edge]
    node1 = idnode2node[id_node1]
    node2 = idnode2node[id_node2]
    edges_nx.append((id_node1, id_node2, {"capacity": capacity, "weight": -len(node2.K)*(len(node1.seq)+len(node2.seq))})) # round((len(node1.K)+len(node2.K)/2)) }))#1}))


# Auxiliar nodes
# one source node
id_source = node2str(block2node(source_block))
edges_nx.append(("s", hash_node[id_source][0], {"capacity":num_seqs,"weight":1}))

# one prev-sink node for each sequence with in-out-capacity-edges=1
last_nodes_by_seq = defaultdict(list)
for node in list_nodes:
    if node.j == size_msa:
        for k in node.K:
            last_nodes_by_seq[k].append(node.K)

# one sink node connected to all the prev-sink
id_sink = node2str(block2node(sink_block))
edges_nx.append((hash_node[id_sink][0], "t", {"capacity":num_seqs,"weight":1}))

#  --------------       Max Flow min cost         ----------------

G = nx.DiGraph()
G.add_edges_from(
    edges_nx
)
mincostFlow = nx.max_flow_min_cost(G, "s", "t")
mincost = nx.cost_of_flow(G, mincostFlow)
maxFlow = maximum_flow(G,"s","t")[1]
nx.cost_of_flow(G, maxFlow) >= mincost

# Pangeblock from max-flow-min-cost
pangeblock = defaultdict(set)
flow_by_edge = dict()
for k,v in maxFlow.items():
    id_node1 = k 
    for id_node2, flow in v.items():

        if flow>0 and id_node1 not in ["s","t"] and id_node2 not in ["s","t"]:
            pangeblock["nodes"] = pangeblock["nodes"].union(set([id_node1, id_node2]))
            pangeblock["edges"] = pangeblock["edges"].union(set([(id_node1, id_node2)]))
            flow_by_edge[(id_node1, id_node2)] = flow


pb = Digraph(format="png")

for id_node in pangeblock["nodes"]:
    node = idnode2node[id_node]
    str_node = node2str(node)
    pb.node(id_node, str_node) 

for edge in pangeblock["edges"]:
    pb.edge(edge[0],edge[1], label = str(flow_by_edge[edge]))

pb.render(directory='output/pb', view=True).replace('\\', '/')