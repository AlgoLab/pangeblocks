configfile: "params.yaml"

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from src.utils import MonitorValuesPlus
from src.msa import AnalyzerMSA
from tqdm import tqdm
from rich.progress import (
    BarColumn,
    # DownloadColumn,
    Progress,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

progress = Progress(
    TextColumn("[bold blue]{task.fields[filename]}", justify="left"),
    BarColumn(bar_width=None),
    "[progress.percentage]{task.percentage:>3.1f}%",
    "•",
    TimeRemainingColumn(),
)

# --- EDA MSA
rule eda_msa: 
    input: 
        path_msas=config["PATH_MSAS"]
    output:
        "out/analysis-msa/stats_msas.tsv",
        "out/analysis-msa/problematic_files.tsv"
    run:
        
        mv = MonitorValuesPlus(list_vars=["path_msa","n_identical_cols","n_cols", "proportion"],
                                out_file="out/analysis-msa/stats_msas.tsv",
                                overwrite=True)

        mv_problems = MonitorValuesPlus(list_vars=["path_msa"],
                                        out_file="out/analysis-msa/problematic_files.tsv",
                                        overwrite=True)
        
        analyzer = AnalyzerMSA()
        def run(path_msa):
            "Function to run with ThreadPoolExecutor"
            try:
                n_cols, n_seqs, n_unique_seqs, n_identical_cols = analyzer(path_msa)
                proportion = round(100*n_identical_cols/n_cols, 2)
                mv()
            except:     
                mv_problems()

        # list of (path to) msas
        list_paths = list(Path(config["PATH_MSAS"]).rglob("*.fa"))
        
        # with progress:
        
        with ThreadPoolExecutor(max_workers=16) as pool:
            with tqdm(total=len(list_paths)) as progress:
                futures=[]
                for path_msa in list_paths:

                    future = pool.submit(run, path_msa)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)
                
                for future in futures:
                    future.result()

                
        