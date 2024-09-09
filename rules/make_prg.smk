configfile: "rules/params/make_prg.yaml"
from pathlib import Path

BIN=config["BIN"]
PATH_MSAS=Path(config["MAKEPRG_INPUT"])
PATH_OUTPUT=Path(config["MAKEPRG_OUTPUT"])
Path(PATH_OUTPUT).mkdir(exist_ok=True, parents=True)

list_fasta = list(Path(PATH_MSAS).glob("*.fa")) #+ list(Path(PATH_MSAS).glob("*.fasta"))
NAMES = [path.stem for path in list_fasta]
print(NAMES)



rule all:
    input:
        expand( PATH_OUTPUT.joinpath("{name_msa}.prg.gfa"), name_msa=NAMES)

rule create_graph:
    input:
        path_msa=PATH_MSAS.joinpath("{name_msa}.fa")
    output:
        path_graph=PATH_OUTPUT.joinpath("{name_msa}.prg.gfa")
    log:
        PATH_OUTPUT.joinpath("logs/{name_msa}.log")
    threads:
        4
    params:
        make_prg=BIN,
        output_prefix=lambda w: PATH_OUTPUT.joinpath(f"{w.name_msa}")
    shell:
        "/usr/bin/time -v {params.make_prg} from_msa --input {input.path_msa} --output-prefix {params.output_prefix} --output-type g --threads {threads} 2> {log}"