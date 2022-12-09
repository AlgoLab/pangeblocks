configfile: "params.yaml"
from pathlib import Path

PATH_MSAS = Path(config["PATH_MSAS"]).rglob("*fa")
NAMES=[Path(path).stem for path in PATH_MSAS] # msas names

rule all:
    input:
        expand("output/max_blocks/{name_msa}.json", name_msa=NAMES)

rule compute_blocks:
    input:
        "/home/avila/data/{name_msa}.fa"
    output: 
        "output/max_blocks/{name_msa}.json"
    shell:
        "python src/compute_blocks.py {input} --output {output}"