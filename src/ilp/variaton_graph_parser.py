"""Parse solution of the ILP formulation as a variation graph"""

from dataclasses import astuple 
from Bio import AlignIO
from collections import defaultdict
from ..graph import nodes_edges_from_blocks
from pathlib import Path
from ..blocks.block import Block
class asGFA:

    def __call__(self, optimal_coverage, path_gfa, path_msa, header="pangeblock"):
        
        list_nodes, list_edges = self.create_graph(optimal_coverage, path_msa)
        self.parse(list_nodes, list_edges, path_gfa, header)

    def create_graph(self, optimal_coverage, path_msa):
        list_nodes = []
        list_edges = []

        _, n_seqs, n_cols= self.load_msa(path_msa)
        optimal_coverage.append(Block((r for r in range(n_seqs)), -1,-1,"s")) # source node
        optimal_coverage.append(Block((r for r in range(n_seqs)), n_cols,n_cols,"S")) # sink node
        
        sorted_solution = sorted(optimal_coverage, key=lambda block: block.i )
        for pos1, block1 in enumerate(sorted_solution[:-1]):

            for rel_pos, block2 in enumerate(sorted_solution[pos1+1:]):
                pos2 = rel_pos + pos1 + 1

                nodes, edges = nodes_edges_from_blocks(block1, block2)
                list_nodes.extend(nodes)
                list_edges.extend(edges)

        list_nodes = set([(node.K,node.i,node.j,node.label) for node in list_nodes])
        
        
        _list_edges = []
        for edge in list_edges:
            node1=edge.node1
            node2=edge.node2

            _list_edges.append(
                ((node1.K,node1.i,node1.j,node1.label),(node2.K,node2.i,node2.j,node2.label), tuple(edge.seqs))
            )

        list_edges = list(set(_list_edges))

        return list_nodes, list_edges

    def parse(self, list_nodes, list_edges, path_gfa,header="pangeblock",):
        
        # graph in GAF format
        lines_gfa = []

        HEADER = f"H\t{header}"

        # segments (nodes)
        lines_segments = [] 
        node2id={} # will be useful to add the lines/edges
        for id_node, node in enumerate(list_nodes):
            lines_segments.append(
                f"S\t{id_node}\t{node[-1]}"
            )

            node2id[node]=id_node

        # links (edges)
        lines_links=[]
        data_paths = defaultdict(list)
        for edge in list_edges:
            node1, node2, seqs = edge[0], edge[1], edge[2]
            lines_links.append(
                f"L\t{node2id[node1]}\t+\t{node2id[node2]}\t+\t0M"
            )

            # paths
            for seq in seqs:
                data_paths[seq].extend([node1,node2])

        lines_paths=[]
        paths = dict()
        for seq, nodes_seq in data_paths.items():
            nodes_seq = list(set(nodes_seq))
            id_nodes = [node2id[node] for node in sorted(nodes_seq, key=lambda node: node[1])]
            paths[seq] = ",".join([str(id_node)+"+" for id_node in id_nodes])
            lines_paths.append(
                f"P\tseq{seq}\t\t{paths[seq]}"
            )
        lines_gfa.append(HEADER)
        lines_gfa.extend(list(set(lines_segments)))
        lines_gfa.extend(list(set(lines_links)))
        lines_gfa.extend(lines_paths)

        path_gfa = Path(path_gfa)
        path_gfa.parent.mkdir(exist_ok=True, parents=True)
        with open(path_gfa,"w") as fp:
            for line in lines_gfa:
                fp.write(line + "\n")

    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align = AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols