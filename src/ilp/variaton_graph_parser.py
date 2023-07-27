"""Parse solution of the ILP formulation as a variation graph"""

from dataclasses import astuple 
from Bio import AlignIO
from collections import defaultdict, namedtuple
from pathlib import Path

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

class asGFA:

    def __call__(self, optimal_coverage, path_gfa, path_msa, header="VN:Z:1.0"):
        
        list_nodes, list_edges = self.create_graph(optimal_coverage, path_msa)
        list_nodes = sorted(list_nodes, key=lambda node: node[1]) # sort nodes by starting position
        self.parse(list_nodes, list_edges, path_gfa, header)

    def create_graph(self, optimal_coverage, path_msa):
        list_nodes = []
        list_edges = []

        self.msa, n_seqs, n_cols= self.load_msa(path_msa)
        
        self.idx2seqid = dict()
        for idx, record in enumerate(self.msa): 
            self.idx2seqid[idx] = record.id 

        sorted_solution = sorted(optimal_coverage, key=lambda block: block.start )
        for pos1, block1 in enumerate(sorted_solution[:-1]):

            for rel_pos, block2 in enumerate(sorted_solution[pos1+1:]):
                pos2 = rel_pos + pos1 + 1

                if block2.start > block1.end + 1:
                    logging.debug("breaking for loop block1.end=%s << block2.start=%s" % (block1.end, block2.start) )
                    break

                nodes, edges = self.nodes_edges_from_blocks(block1, block2)
                list_nodes.extend(nodes)
                list_edges.extend(edges)

                if len(nodes)>0:
                    logging.debug("adding nodes: %s", nodes)
                else:
                    logging.debug("blocks not connected: %s and %s " % (block1, block2))

        logging.debug("Number of nodes: %s", len(list_nodes))
        list_nodes = list(set([(node.K,node.i,node.j,
                                str(self.msa[int(node.K[0]),int(node.i):int(node.j)+1].seq) ) for node in list_nodes]))
        logging.debug("Number of unique nodes %s", len(list_nodes))
        
        _list_edges = []
        for edge in list_edges:
            node1=edge.node1
            node2=edge.node2

            _list_edges.append(
                ((node1.K,node1.i,node1.j,node1.label),(node2.K,node2.i,node2.j,node2.label), tuple(edge.seqs))
            )
        
        logging.debug("Number of edges %s", len(_list_edges))
        list_edges = list(set(_list_edges))
        logging.debug("Number of unique edges %s", len(list_edges))

        logging.debug("Number of nodes == Number of blocks in opt coverage? %s", len(list_nodes)==len(optimal_coverage))
        return list_nodes, list_edges

    def parse(self, list_nodes, list_edges, path_gfa,header="VN:Z:1.0",):
        
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
            node1, node2, idxseqs = edge[0], edge[1], edge[2]
            lines_links.append(
                f"L\t{node2id[node1]}\t+\t{node2id[node2]}\t+\t0M"
            )

            # paths
            for idx in idxseqs:
                data_paths[self.idx2seqid[idx]].extend([node1,node2])

        lines_paths=[]
        paths = dict()
        for seqid, nodes_seq in data_paths.items():
            nodes_seq = list(set(nodes_seq))
            nodes_seq = sorted(nodes_seq, key=lambda node: node[1]) # sort by starting position
            id_nodes = [node2id[node] for node in nodes_seq]
            paths[seqid] = ",".join([str(id_node)+"+" for id_node in id_nodes])
            lines_paths.append(
                f"P\t{seqid}\t{paths[seqid]}\t*"
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

    def nodes_edges_from_blocks(self, block1, block2):
        Node = namedtuple("Node",["K","i","j","label"]) # is a block
        Edge = namedtuple("Edge",["node1","node2","seqs"])

        b1, b2 = block1, block2
        nodes = []
        edges = []
        # not empty intersection
        common_rows = set(b1.K).intersection(set(b2.K))
        # K = len(common_rows)# number of seqs that can traverse the edge
        K = list(common_rows)
        K.sort()

        if b1.end == b2.start-1 and len(K)>0:
            b1_label = self.msa[b1.K[0]].seq[b1.start:b1.end+1]
            b2_label = self.msa[b2.K[0]].seq[b2.start:b2.end+1]

            # print("Condicion- consecutive blocks")
            node1 = Node(b1.K, b1.start, b1.end, b1_label)
            node2 = Node(b2.K, b2.start, b2.end, b2_label)
            nodes.extend([node1, node2])
            edges.append(Edge(node1, node2, K))
        
        if b2.end == b1.start-1 and len(K)>0:
            logging.debug("block %s and block %s not connected" % (b1.str(), b2.str()))
        # else: 
        #     # print("Not consecutive blocks")
        return nodes, edges