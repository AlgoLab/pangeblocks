"""Postprocessing of GFA
modification of 
https://github.com/AlgoLab/RecGraph-exps/blob/main/scripts/clean_gfa_from_ast.py
"""
from collections import defaultdict

def postprocessing(path_gfa): 
    """
    remove nodes with '-'
    remove '-' in the labels of nodes
    """
    to_remove = []
    predecessors = defaultdict(list)
    successors = defaultdict(list)

    # segments (nodes)
    for line in open(path_gfa):
        if not line.startswith("S"):
            continue
        
        _, idx, seq = line.split("\t")
        seq = seq.replace("\n","")
        
        if set(seq) == set("-"):
            # nodes labeled with "-" will be removed
            to_remove.append(idx)
            predecessors[idx] = []
            successors[idx] = []
        elif "-" in seq:
            # remove "-" from the label
            seq = seq.replace("-","")
            print("S", idx, seq, sep="\t") 
        else:
            # otherwise, keep the node
            print(line, end="")

    # links (edges)
    for line in open(path_gfa):
        if not line.startswith("L"):
            continue
        
        _, idx1, _, idx2, _, _ = line.split("\t")
        if idx1 in to_remove:
            successors[idx1].append(idx2)
        elif idx2 in to_remove:
            predecessors[idx2].append(idx1)
        else:
            print(line, end="")


    for idx in successors:
        while any([x in to_remove for x in successors[idx]]):
            new_successors = []
            for sidx in successors[idx]:
                if sidx in to_remove:
                    new_successors.extend(successors[sidx])
                else:
                    new_successors.append(sidx)
            successors[idx] = new_successors

    for idx in predecessors:
        while any([x in to_remove for x in predecessors[idx]]):
            new_predecessors = []
            for sidx in predecessors[idx]:
                if sidx in to_remove:
                    new_predecessors.extend(predecessors[sidx])
                else:
                    new_predecessors.append(sidx)
            predecessors[idx] = new_predecessors

    # adds new edges
    for idx in to_remove:
        if len(predecessors[idx]) == 0 or len(successors[idx]) == 0:
            continue
        for p in predecessors[idx]:
            for s in successors[idx]:
                print("L", p, "+", s, "+", "0M", sep="\t")

    # modify paths
    for line in open(path_gfa):
        if not line.startswith("P"):
            continue

        _, seq_id, path = line.split("\t")
        
        path     = path.replace("\n","").replace("+","")
        idx_path = path.split(",")
        
        new_path = []
        for idx in idx_path:
            if idx in to_remove:
                # find predecessor and successor in the path
                # there should be only one predecessor and one successor
                predecessor = [pred for pred in predecessors[idx] if pred in idx_path][0] 
                successor   = [suc  for suc  in successors[idx] if suc in idx_path][0]
                new_path.extend([predecessor, successor])
            else: 
                new_path.append(idx)
        # from pdb import set_trace
        # set_trace()
        new_path = ",".join([idx+"+" for idx in new_path])
        print("P", seq_id, new_path, sep="\t")