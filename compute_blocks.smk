configfile: "params.yaml"
from pathlib import Path
import pandas as pd

# PATH_MSAS = Path(config["PATH_MSAS"]).rglob("*fa")
# NAMES=[Path(path).stem for path in PATH_MSAS["path_msa"]] # msas names
NAMES = ["GC00001810"]
PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]

rule all:
    input:
        expand("{path_output}/max_blocks/stats/{name_msa}.tsv", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/block_decomposition/{name_msa}.json", name_msa=NAMES, path_output=PATH_OUTPUT)

rule compute_blocks:
    input:
        expand("{path_msas}/{{name_msa}}.fa", path_msas=PATH_MSAS)
    output: 
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    shell:
        "python src/compute_blocks.py {input} --output {output}"

rule analyze_blocks:
    input:
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    output:
        expand("{path_output}/max_blocks/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
    shell: 
        "python analyze_blocks.py {input} --output {output}"

rule decompose_blocks:
    input:
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    output:
        expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
        expand("{path_output}/block_decomposition/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
    shell:
        "python decompose_blocks.py {input} --output {output[0]} --output-stats {output[1]}"
