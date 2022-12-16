# https://graphviz.readthedocs.io/en/stable/manual.html
from graphviz import Digraph
from pathlib import Path

class PlotGraph:
    "Plot graph from a set of nodes and a list of edges"

    def __call__(self, nodes, edges, path_save: str):
        # Create an empty graph
        self.digraph  = Digraph("Pangeblock")

        nodes_plot = [self.node2str(node) for node in nodes]
        edges_plot = [(self.node2str(edge[0]),self.node2str(edge[1])) for edge in edges]
        self.digraph = self.create_dot_graph(self.digraph,nodes_plot, edges_plot)
        
        path_save = Path(path_save)
        if path_save.suffix == ".dot":
            path_save.parent.mkdir(parents=True, exist_ok=True)
            self.save_as_dot(path_save)
        else: 
            raise("Path not valid, must be  a '.dot' file")
        

    def create_dot_graph(self, graph, nodes, edges):
        "Create graph in 'dot' language to be plotted"
        # Add nodes to the graph
        for node in nodes: 
            graph.node(node)
        # Add edges to the graph
        for edge in edges: 
            u, v = edge
            graph.edge(u,v)

        return graph

    def show(self, graph):
        "Print graph"
        # self.create_dot_graph(self.nodes, self.edges)
        return graph #self.digraph

    def save_as_dot(self, path_save: str = "graph.dot"):
        " Save tree as .dot file to generate plot with grahviz"
        self.digraph.render(filename=path_save)

    def node2str(self, node):
        return f"K=({','.join(str(_) for _ in node.K)}),i={node.i},j={node.j},{node.label}"