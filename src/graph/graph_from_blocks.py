# not used
"""Given a set of non overlapping and consecutives blocks by pairs,
create edges and nodes for DAG"""

from collections import namedtuple
# from ...blocks import Block

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='[Graph from blocks] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

Node = namedtuple("Node",["K","i","j","label"]) # is a block
Edge = namedtuple("Edge",["node1","node2","seqs"])

# def label_from_msa()

def nodes_edges_from_blocks(block1, block2):
    b1, b2 = block1, block2
    nodes = []
    edges = []
    # not empty intersection
    common_rows = set(b1.K).intersection(set(b2.K))
    # K = len(common_rows)# number of seqs that can traverse the edge
    K = list(common_rows)
    K.sort()

    if b1.end == b2.start-1 and len(K)>0:
        # print("Condicion- consecutive blocks")
        node1 = Node(b1.K, b1.start, b1.end, b1.label)
        node2 = Node(b2.K, b2.start, b2.end, b2.label)
        nodes.extend([node1, node2])
        edges.append(Edge(node1, node2, K))
    
    if b2.end == b1.start-1 and len(K)>0:
        logging.debug("block %s and block %s not connected" % (b1.str(), b2.str()))
    # else: 
    #     # print("Not consecutive blocks")
    return nodes, edges