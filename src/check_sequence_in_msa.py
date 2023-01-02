
import argparse 
from collections import defaultdict
from Bio import AlignIO
from graph.gfa_loader import load_gfa

# Command line options
parser = argparse.ArgumentParser()
parser.add_argument("--source","-s", help="id source node")
parser.add_argument("--sink","-S", help="id sink node")
parser.add_argument("--path-gfa", help="gfa file")
parser.add_argument("--path-msa", help="msa file in fasta format")
args = parser.parse_args()

nodes, edges, paths, in_nodes, out_nodes = load_gfa(args.path_gfa)
msa=AlignIO.read(args.path_msa,"fasta")

def check_sequence(seq, source, sink):
    source_node=source
    sink_node=sink
    # booleans
    _visited = defaultdict(bool)
    _current_pos_seq = defaultdict(int) # save path and position until the seq spelt by the path matches the input seq 

    _visited[source_node] = True
    _current_pos_seq[source_node] = 0

    nodes_to_visit = [source_node]

    node_id = -1 # reset node_id

    while nodes_to_visit and node_id != sink_node:
        
        current_node = nodes_to_visit.pop()
        current_pos_seq = _current_pos_seq[current_node]
        
        for node_id in out_nodes[current_node]:    
            
            if node_id == sink_node:
                break # path exists

            label = nodes[node_id]["label"].upper() # string node
            subseq = str(seq[current_pos_seq:current_pos_seq+len(label)]).upper() # string sequence

            if label == "*": # case make_prg
                # _current_pos_seq[node_id] = current_pos_seq
                nodes_to_visit.append(node_id)
                _current_pos_seq[node_id] = current_pos_seq
            
            elif subseq == label: # there is a match
                # print(subseq, label)
                _current_pos_seq[node_id] = current_pos_seq + len(label)
                nodes_to_visit.append(node_id)
            
            else: # no match
                _visited[node_id] = True

    return node_id == sink_node

for num,record in enumerate(msa):
    seq = record.seq
    print(num, check_sequence(seq, int(args.source), int(args.sink)))