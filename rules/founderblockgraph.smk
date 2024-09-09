configfile: "rules/params/founderblockgraph.yaml"
from pathlib import Path

BIN=config["BIN"]
PATH_MSAS=Path(config["FBG_INPUT"])
PATH_OUTPUT=Path(config["FBG_OUTPUT"])
Path(PATH_OUTPUT).mkdir(exist_ok=True, parents=True)

list_fasta = list(Path(PATH_MSAS).glob("*.fa")) #+ list(Path(PATH_MSAS).glob("*.fasta"))
NAMES = [path.stem for path in list_fasta]
print(NAMES)



rule all:
    input:
        expand( PATH_OUTPUT.joinpath("{name_msa}.xgfa"), name_msa=NAMES)

rule create_graph:
    input:
        path_msa=PATH_MSAS.joinpath("{name_msa}.fa")
    output:
        path_graph=PATH_OUTPUT.joinpath("{name_msa}.xgfa")
    log:
        PATH_OUTPUT.joinpath("logs/{name_msa}.log")
    threads:
        4
    params:
        founderblockgraph=BIN
    shell:
        "/usr/bin/time -v {params.founderblockgraph} --input={input.path_msa} --output={output.path_graph} --elastic --gfa --output-paths --threads={threads} 2> {log}"