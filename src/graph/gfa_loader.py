from collections import defaultdict

def load_gfa(path_gfa):
    nodes=dict()
    edges=list()
    paths=dict()
    in_nodes = defaultdict(list)
    out_nodes = defaultdict(list)

    with open(path_gfa,"r") as fp: 
        for line in fp.readlines():
            # print(line)
            line = line.replace("\n","")
            line = line.split("\t")
            # node
            if line[0]=="S":
                id_node = int(line[1])
                label_node = line[2].upper()
                nodes[id_node]=dict(label=label_node)
            # edge
            elif line[0]=="L":
                start_node = int(line[1])
                end_node = int(line[3])
                in_nodes[end_node].append(start_node)
                out_nodes[start_node].append(end_node)

                edges.append((start_node,end_node))

    return nodes, edges, paths, in_nodes, out_nodes
