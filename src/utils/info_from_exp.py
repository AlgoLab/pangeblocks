import json
import pandas as pd
from pathlib import Path

from src.graph.info_from_gfa import load_info_graph
from .params_from_path import params_from_path

def get_info_from_exp(path_exp, name):
    NAME=name
    info = []
    # config info from the experiment
    with open(Path(path_exp).joinpath("config.json"), "r") as fp:
        config = json.load(fp)
    decomposition = "complete" if config["DECOMPOSITION"]["STANDARD"] is True else "row-maximal"
    alpha_consistent = True if config["DECOMPOSITION"]["ALPHA_CONSISTENT"] is True else False

    # load info of graphs for each experiment
    path_graphs = Path(path_exp).joinpath("gfa-unchop").rglob(f"*{NAME}.gfa")
    for path in path_graphs:
        params = params_from_path(path, output_dict=True)
        params["alpha_consistent"] = alpha_consistent
        params["decomposition"] = decomposition
        params["path"] = path

        info_pb = load_info_graph(path)        
        params.update(info_pb)

        info.append(params)

    return pd.DataFrame(info)

def get_info_from_ilp(path_exp, name):
    NAME=name
    info = []
    # config info from the experiment
    with open(Path(path_exp).joinpath("config.json"), "r") as fp:
        config = json.load(fp)
    decomposition = "complete" if config["DECOMPOSITION"]["STANDARD"] is True else "row-maximal"
    alpha_consistent = True if config["DECOMPOSITION"]["ALPHA_CONSISTENT"] is True else False

    # load info of graphs for each experiment
    path_ilp = Path(path_exp).joinpath("ilp").rglob(f"{NAME}*.gfa")
    for path in path_ilp:
        params = params_from_path(path, output_dict=True)
        params["alpha_consistent"] = alpha_consistent
        params["decomposition"] = decomposition
        params["path"] = path

        info_pb = load_info_graph(path)        
        params.update(info_pb)

        info.append(params)

    return pd.DataFrame(info)