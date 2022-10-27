"""
Build graph for maximum flow
"""
from collections import namedtuple

Node = namedtuple("Node",["K","i","j","seq"]) # is a block
Edge = namedtuple("Edge",["node1","node2","capacity"])
Block = namedtuple("Block",["K","i","j","seq"])

def nodes_edges_from_blocks(block1, block2):
    "Generate nodes and edges from 2 blocks based on their location and common shared sequences" 
    b1, b2 = block1, block2
    nodes = []
    edges = []
    # not empty intersection
    common_rows = set(b1.K).intersection(set(b2.K))
    capacity = len(common_rows)# number of seqs that uses the nodes that form he edge

    if common_rows: 
        K = list(common_rows)
        K.sort()

        # Condicion1
        if b1.i == b2.i and b1.j < b2.j:
            print("Condicion1")
            # nodes
            node1 = Node(b1.K, b1.i, b1.j, b1.seq)
            node2 = Node(b2.K, b1.j+1, b2.j, b2.seq[b1.j-b1.i+1:])
            nodes.extend([node1, node2])
            
            # edges
            edges.append(Edge(node1, node2, capacity))

        # Condicion2
        elif b1.i < b2.i and b1.j == b2.j:
            print("Condicion2")
            node1 = Node(b1.K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.seq)
            nodes.extend([node1, node2])

            # edges
            edges.append(Edge(node1,node2,capacity))
        
        # Condicion3
        elif b1.i < b2.i and b2.j < b1.j:
            print("Condicion3")
            node1 = Node(b1.K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.seq)
            node3 = Node(b1.K, b2.j+1, b1.i, b1.seq[:b1.i-(b2.j+1)+1])
            
            nodes.extend([node1, node2, node3])
            # edges
            edges.append(Edge(node1, node2, capacity))
            edges.append(Edge(node2, node3, capacity))     
        
        # Condicion4
        elif b1.i < b2.i and b2.i < b1.j and b1.j < b2.j:
            print("Condicion4")
            # option 1 
            node1 = Node(b1.K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            node2 = Node(b2.K, b2.i, b2.j, b2.seq)
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2,capacity))

            # option 2
            node1 = Node(b1.K, b1.i, b1.j, b1.seq)
            node2 = Node(b2.K, b1.j+1, b2.j, b2.seq[b1.j+1-b2.i:])
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2, capacity))

        # consecutive blocks -> connect them
        # TODO: verificar si debe moverse al inicio, puede tener problemas con la condicion 4
        if b1.j == b2.i-1:
            print("Condicion- consecutive blocks")
            node1 = Node(b1.K, b1.i, b1.j, b1.seq)
            node2 = Node(b2.K, b2.i, b2.j, b2.seq)
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2, capacity))

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
        # if node1 == Node(source_block.K, source_block.i, source_block.j, source_block.seq):
        #     print(node1,node2)
        return [Edge(node1, node2, capacity)]
    return []