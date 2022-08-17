"""
Build graph for maximum flow
"""
from collections import namedtuple

Node = namedtuple("Node",["K","i","j","seq"]) # is a block
Edge = namedtuple("Edge",["node1","node2","capacity"])
Block = namedtuple("Block",["K","i","j","seq"])

def nodes_edges_from_blocks(block1, block2): 
    b1, b2 = block1, block2
    nodes = []
    edges = []
    # not empty intersection
    common_rows = set(b1.K).intersection(set(b2.K))
    capacity = len(common_rows)# number of seqs that uses the nodes that form he edge
    # list_blocks = [block1, block2]
    # list_blocks.sort(key = lambda: )
    if common_rows: 
        K = list(common_rows)
        K.sort()

        # Condicion1
        if b1.i == b2.i and b1.j < b2.j:
            print("Condicion1")
            # nodes
            node1 = Node(b1.K, b1.i, b1.j, b1.seq)
            node2 = Node(b2.K, b1.j+1, b2.j, b2.seq[b1.j-b1.i+1:])#b2.seq[b2.j-b1.j+1:])
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

def disect_blocks(block1, block2): 
    b1, b2 = block1, block2
    blocks = []
    
    # not empty intersection
    K = list(set(b1.K).intersection(set(b2.K)))
    K.sort()
    if len(K)>0: 
        
        if b1.i == b2.i and b1.j < b2.j:
            block1 = Block(b1.K, b1.i, b1.j, b1.seq)
            block2 = Block(K, b1.j+1, b2.j, b2.seq[:b2.j-(b1.j+1)+1])
            blocks.extend([block1, block2])
            
        elif b1.i < b2.i and b1.j == b2.j:
            block1 = Block(K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            block2 = Block(b2.K, b2.i, b2.j, b2.seq)
            blocks.extend([block1, block2])
        
        elif b1.i < b2.i and b2.j < b1.j:
            block1 = Block(K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            block2 = Block(b2.K, b2.i, b2.j, b2.seq)
            block3 = Block(K, b2.j+1, b1.i, b2.seq[:b1.i-(b2.j+1)+1])
            blocks.extend([block1, block2, block3])
        
        elif b1.i < b2.i and b2.i < b1.j and b1.j < b2.j:            
            # option 1 
            block1 = Block(K, b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i+1])
            block2 = Block(b2.K, b2.i, b2.j, b2.seq)
            blocks.extend([block1, block2])

            # option 2
            block1 = Block(b1.K, b1.i, b1.j, b1.seq)
            block2 = Block(K, b1.j+1, b2.j, b2.seq[:b2.j-(b1.j+1)+1])
            blocks.extend([block1, block2])

        # consecutive blocks -> connect them
        elif b1.j == b2.i -1:
            block1 = Block(K, b1.i, b1.j, b1.seq)
            block2 = Block(K, b2.i, b2.j, b2.seq)
            blocks.extend([block1, block2])

    return blocks

# class GraphBuilder:

#     def __init__(self, max_blocks):
#         self.max_blocks = max_blocks

#     def nodes_from_blocks(self, block1, block2):
#         "Create nodes based on blocks intersection"
#         b1, b2 = block1, block2
#         nodes = []
#         edges = []
#         # not empty intersection
#         common_rows = set(b1.K).intersection(set(b2.K))
#         capacity = len(common_rows)# number of seqs that uses the nodes that form he edge
#         if common_rows: 
            
#             if b1.i == b2.i and b1.j < b2.j:
#                 # nodes
#                 node1 = Node(b1.i, b1.j, b1.seq)
#                 node2 = Node(b1.j+1, b2.j, b2.seq[:b2.j-(b1.j+1)])
#                 nodes.extend([node1, node2])
                
#                 # edges
#                 edges.append(Edge(node1, node2, capacity))

#             elif b1.i < b2.i and b1.j == b2.j:
#                 node1 = Node(b1.i, b2.i-1, b1.seq[b2.i-1-b1.i])
#                 node2 = Node(b2.i, b2.j, b2.seq)
#                 nodes.extend([node1, node2])

#                 # edges
#                 edges.append(Edge(node1,node2,capacity))
            
#             elif b1.i < b2.i and b2.j < b1.j:
#                 node1 = Node(b1.i, b2.i-1, b1.seq[:b2.i-1-b1.i])
#                 node2 = Node(b2.i, b2.j, b2.seq)
#                 node3 = Node(b2.j+1, b1.i, b2.seq[:b1.i-(b2.j+1)])

#                 # edges
#                 edges.append(Edge(node1, node2, capacity))
#                 edges.append(Edge(node2, node3, capacity))     
            
#             elif b1.i < b2.i and  b2.i < b1.j and b1.j < b2.j:
#                 # option 1 
#                 node1 = Node(b1.i, b2.i-1, b1.seq[b2.i-1-b1.i])
#                 node2 = Node(b2.i, b2.j, b2.seq)
#                 edges.append(Edge(node1, node2,capacity))

#                 # option 2
#                 node1 = Node(b1.i, b1.j, b1.seq)
#                 node2 = Node(b1.i+1, b2.j, b2.seq[:b2.j-(b1.j)])
#                 edges.append(Edge(node1, node2, capacity))
