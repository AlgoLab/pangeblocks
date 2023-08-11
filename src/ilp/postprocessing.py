"""Postprocessing of GFA
modification of 
https://github.com/AlgoLab/RecGraph-exps/blob/main/scripts/clean_gfa_from_ast.py
"""
from pathlib import Path
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO,
                    format='[Postprocessing] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
                    

def postprocessing(path_gfa, path_save): 
    """
    remove nodes with '-'
    remove '-' in node labels
    """
    HEADER = ""
    lines_new_gfa = []

    to_remove = []
    predecessors = defaultdict(list) # new predecessors of nodes
    successors = defaultdict(list)   # new successors of nodes

    # segments (nodes)
    for line in open(path_gfa):
        # if header exists, we save it to include it in the new GFA
        if line.startswith("H"):
            HEADER = line
        if not line.startswith("S"):
            continue
        
        _, idx, seq = line.split("\t")
        seq = seq.replace("\n","")
        
        if set(seq) == set("-"):
            # nodes labeled with "-" will be removed
            to_remove.append(idx)
            # predecessors[idx] = []
            # successors[idx] = []
        elif "-" in seq:
            # remove "-" from the label
            seq = seq.replace("-","")
            # print("S", idx, seq, sep="\t") 
            lines_new_gfa.append( 
                "\t".join(["S", idx, seq]) + "\n"
            )
        else:
            # otherwise, keep the node
            # print(line, end="")
            lines_new_gfa.append( 
                line
            )

    for idx in to_remove:
        logging.info(
            f"to remove {idx}"
        )

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
            # print(line, end="")
            lines_new_gfa.append( 
                line
            )

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
                # print("L", p, "+", s, "+", "0M", sep="\t")
                lines_new_gfa.append( 
                    "\t".join(["L", p, "+", s, "+", "0M"]) + "\n"
                )

    # modify paths
    for line in open(path_gfa):
        if not line.startswith("P"):
            continue

        _, seq_id, path, _ = line.split("\t")
        
        path     = path.replace("\n","").replace("+","")
        idx_path = path.split(",")
        
        new_path = []
        for idx in idx_path:
            if idx not in to_remove:
                new_path.append(idx)

        # FIXME: remove duplicates (should not be necessary)
        clean_path = []
        for idx in new_path:
            if idx not in clean_path:
                clean_path.append(idx)

        clean_path = ",".join([idx+"+" for idx in clean_path])
        clean_path = clean_path.strip()
        # print("P", seq_id, clean_path, sep="\t")
        lines_new_gfa.append( 
                "\t".join(["P", seq_id, clean_path, "*"]) + "\n"
            )

        # save new gfa
        Path(path_save).parent.mkdir(exist_ok=True, parents=True)
        with open(path_save, "w") as fp:
            if HEADER:
                fp.write(HEADER)
            fp.writelines(lines_new_gfa)     