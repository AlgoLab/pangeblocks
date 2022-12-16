"""Given a set of non overlapping and consecutives blocks by pairs,
create edges and nodes for DAG"""

from collections import namedtuple
from ..blocks import Block

Node = namedtuple("Node",["K","i","j","label"]) # is a block
Edge = namedtuple("Edge",["node1","node2","seqs"])

def nodes_edges_from_blocks(block1: Block, block2: Block):
    b1, b2 = block1, block2
    nodes = []
    edges = []
    # not empty intersection
    common_rows = set(b1.K).intersection(set(b2.K))
    K = len(common_rows)# number of seqs that can traverse the edge
    K = list(common_rows)
    K.sort()

    if b1.j == b2.i-1:
        print("Condicion- consecutive blocks")
        node1 = Node(b1.K, b1.i, b1.j, b1.label)
        node2 = Node(b2.K, b2.i, b2.j, b2.label)
        nodes.extend([node1, node2])
        edges.append(Edge(node1, node2, K))
    else: 
        print("Not consecutive blocks")
    return nodes, edges