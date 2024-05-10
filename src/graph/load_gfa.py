from collections import defaultdict

NUC_COMPLEMENT = {n:c for n,c in zip ("ACGT","TGCA")}
def reverse_complement(seq):
    seq = seq[::-1]
    return "".join(NUC_COMPLEMENT[n] for n in seq)


def load_gfa(path_gfa, return_nodes_path=True):
    nodes=dict()
    edges=[]
    paths=dict()
    with open(path_gfa, "r") as fp:
        for line in fp.readlines():    
            # nodes
            if line.startswith("S"):
                try:
                    _, nodeid, label = line.replace("\n","").split("\t")
                except:
                    _, nodeid, label, _ = line.replace("\n","").split("\t")
                nodes[nodeid] = {"label": label, "len": len(label)}

            # edges: L	4	+	86	+	0M
            if line.startswith("L"):
                _, nodeid1, _, nodeid2, _, _ = line.replace("\n","").split("\t") 
                edges.append((nodeid1, nodeid2))

    # once nodes are loaded, check the paths
    with open(path_gfa, "r") as fp:
        for line in fp.readlines():    

            # paths
            nodes_path=[]
            if line.startswith("P"):
                _, seq_id, path, *_ = line.replace("\n","").split("\t")

                labels_path=[]
                for node in path.split(","):
                    if "-" in node: # reverse label '-' next to node
                        nodeid = node.replace("-","")
                        label  = reverse_complement(nodes[nodeid]["label"])
                        labels_path.append(label.upper())
                    else: # forward label '+' next to node
                        nodeid = node.replace("+","")
                        label  = nodes[nodeid]["label"]
                        labels_path.append(label.upper())
                    nodes_path.append(nodeid)
                
                if return_nodes_path is True:
                    paths[seq_id] = nodes_path
                else:
                    seq_path   = "".join(labels_path)
                    paths[seq_id] = seq_path

    node_depth = defaultdict(int)
    for seq_id, _nodes in paths.items():
        for node in _nodes: 
            node_depth[node] += 1

    n_paths = len(paths)
    node_depth = {k: v / n_paths for k,v in node_depth.items()} # return node depth in [0,1]
    return nodes, edges, paths, node_depth