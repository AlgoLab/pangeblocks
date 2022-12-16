"""
Create nodes and edges from blocks
"""
from collections import namedtuple
from ..blocks import Block

Node = namedtuple("Node",["K","i","j","label"]) # is a block
Edge = namedtuple("Edge",["node1","node2","seqs"])

def nodes_edges_from_blocks(block1, block2):
    "Generate nodes and edges from 2 blocks based on their location and common shared sequences" 
    b1, b2 = block1, block2
    nodes = []
    edges = []
    # not empty intersection
    common_rows = set(b1.K).intersection(set(b2.K))
    K = len(common_rows)# number of seqs that can traverse the edge

    if common_rows: 
        K = list(common_rows)
        K.sort()

        # Condicion1
        if b1.i == b2.i and b1.j < b2.j:
            print("Condicion1")
            # nodes
            node1 = Node(b1.K, b1.i, b1.j, b1.label)
            node2 = Node(b2.K, b1.j+1, b2.j, b2.label[b1.j-b1.i+1:])
            nodes.extend([node1, node2])
            
            # edges
            edges.append(Edge(node1, node2, K))

        # Condicion2
        elif b1.i < b2.i and b1.j == b2.j:
            print("Condicion2")
            node1 = Node(b1.K, b1.i, b2.i-1, b1.label[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.label)
            nodes.extend([node1, node2])

            # edges
            edges.append(Edge(node1,node2,K))
        
        # Condicion3
        elif b1.i < b2.i and b2.j < b1.j:
            print("Condicion3")
            node1 = Node(b1.K, b1.i, b2.i-1, b1.label[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.label)
            node3 = Node(b1.K, b2.j+1, b1.i, b1.label[:b1.i-(b2.j+1)+1])
            
            nodes.extend([node1, node2, node3])
            # edges
            edges.append(Edge(node1, node2, K))
            edges.append(Edge(node2, node3, K))     
        
        # Condicion4
        elif b1.i < b2.i and b2.i < b1.j and b1.j < b2.j:
            print("Condicion4")
            # option 1 
            node1 = Node(b1.K, b1.i, b2.i-1, b1.label[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.label)
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2,K))

            # option 2
            node1 = Node(b1.K, b1.i, b1.j, b1.label)
            node2 = Node(b2.K, b1.j+1, b2.j, b2.label[b1.j+1-b2.i:])
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2, K))

        # consecutive blocks -> connect them
        # TODO: verificar si debe moverse al inicio, puede tener problemas con la condicion 4
        if b1.j == b2.i-1:
            print("Condicion- consecutive blocks")
            node1 = Node(b1.K, b1.i, b1.j, b1.label)
            node2 = Node(b2.K, b2.i, b2.j, b2.label)
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2, K))

    return nodes, edges

def is_consecutive(pair1, pair2):
    "consecutive cols in the coverage, pair=(id_seq, num_col)"
    seq1, col1 = pair1
    seq2, col2 = pair2 
    if seq1 == seq2 and col1 == col2-1:
        return True
    return False

def missing_edge(node1, node2):
    "missing edges between consecutives block-nodes"
    b1, b2 = node1, node2
    common_rows = set(b1.K).intersection(set(b2.K))
    capacity = len(common_rows)
    if b1.j == b2.i-1 and capacity>0:
        # if node1 == Node(source_block.K, source_block.i, source_block.j, source_block.label):
        #     print(node1,node2)
        return [Edge(node1, node2, capacity)]
    return []