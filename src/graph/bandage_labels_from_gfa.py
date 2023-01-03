import argparse
from pathlib import Path
from collections import namedtuple, defaultdict

def bandage_labels(path_gfa: str, path_save_labels= str):
    """
    Writes a CSV file with labels for Bandage
    """
    assert Path(path_gfa).suffix == ".gfa", "Input must be a gfa file"

    BandageLabels = namedtuple("BandageLables",["node_name","first_base","last_base", "seq"])
    info_labels = []
    seqs_by_node = defaultdict(list)

    with open(path_gfa) as fp:
        for line in fp:
            line=line.replace("\n","")
            if line.startswith("S"):
                split_line  = line.split("\t")
                node_name = split_line[1]
                seq = split_line[2]
                first_base=seq[0]
                last_base=seq[-1]
                info_labels.append(BandageLabels(node_name, first_base, last_base, seq))
            
            if line.startswith("P"):
                split_line = line.split("\t")
                seq_id = split_line[1]
                path   = split_line[-1].replace("+","")
                nodes_path = [int(node_id) for node_id in path.split(",")]
                for node_id in nodes_path:
                    seqs_by_node[node_id].append(seq_id)


    # Write info in csv file
    colnames = ",".join(BandageLabels._fields)+ ",list_seqs" + "\n"

    path_save = Path(path_save_labels)
    path_save.parent.mkdir(exist_ok=True, parents=True)
    with open(path_save,"w") as fp:
        fp.write(colnames)
        for label in info_labels:
            list_seqs = seqs_by_node[int(label.node_name)]
            list_seqs.sort()
            list_seqs = "-".join(list_seqs)
            fp.write(f"{label.node_name},{label.first_base},{label.last_base},{label.seq},{list_seqs}\n")

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--path_gfa")
    parser.add_argument("--path_save")
    args = parser.parse_args()
    bandage_labels(args.path_gfa, args.path_save)
