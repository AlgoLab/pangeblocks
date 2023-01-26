configfile: "params.yaml"
from pathlib import Path

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
Path(PATH_OUTPUT).joinpath("output-pandora").mkdir(exist_ok=True, parents=True)

NAMES = [path.stem for path in Path(PATH_MSAS).rglob("*.fa")]

rule all:
    input:
        expand("{path_output}/output-pandora/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT)

rule download_make_prg:
    output:
        "make_prg_0.4.0"
    shell:
        """
        wget https://github.com/iqbal-lab-org/make_prg/releases/download/0.4.0/make_prg_0.4.0
        chmod +x make_prg_0.4.0
        """

rule generate_prg:
    input:
        "make_prg_0.4.0"
    output:
        f"{PATH_OUTPUT}/output-pandora/.prg.gfa.zip"
    benchmark:
        f"{PATH_OUTPUT}/output-pandora/benchmarks/make_prg.benchmark.txt"
    log:
        f"{PATH_OUTPUT}/output-pandora/make_prg.log"
    shell:
        f" /usr/bin/time --verbose ./{{input}} from_msa -i {PATH_MSAS} -o {PATH_OUTPUT}/output-pandora/ --output-type g 2> {{log}}"

rule extract_gfa: 
    input:
        path_zip=f"{PATH_OUTPUT}/output-pandora/.prg.gfa.zip",
        path_output=PATH_OUTPUT
    output: 
        expand("{path_output}/output-pandora/{name_msa}.gfa", name_msa=NAMES,path_output=PATH_OUTPUT)
    shell:
        "unzip {input.path_zip} -d {input.path_output}/output-pandora"