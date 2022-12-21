"Read a graph in gfa format and compute some stats"
from typing import Union
from pathlib import Path
from collections import Counter

class GFAStats:

    def __call__(self, path_gfa: Union[str,Path]):
        nodes_by_id = self.get_nodes(path_gfa)
        lens = [len(label) for label in nodes_by_id.values()]
        return lens

    def get_nodes(self, path_gfa):
        nodes = dict()
        with open(path_gfa, "r") as fp:
            for line in fp.readlines():
                line = line.replace("\n","")
                line_split = line.split("\t")
                
                if line_split[0] == "S": 
                    id_node = line_split[1]
                    label   = line_split[2]
                    nodes[id_node] = label

        return nodes