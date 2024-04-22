configfile: "rules/params/vg.yaml"
from pathlib import Path
from os.path import join as pjoin

PATH_MSAS=config["PATH_MSAS"]

PATH_OUTPUT=config["VG_OUTPUT"]
Path(PATH_OUTPUT).mkdir(exist_ok=True, parents=True)

list_fasta = list(Path(PATH_MSAS).glob("*.fa")) #+ list(Path(PATH_MSAS).glob("*.fasta"))
NAMES = [path.stem for path in list_fasta]
print(NAMES)

rule all:
    input:
        expand( pjoin(PATH_OUTPUT,"{name_msa}.gfa"), name_msa=NAMES)

rule vg_graph: 
    input: 
        path_msa=pjoin(PATH_MSAS, "{name_msa}.fa")
    output: 
        path_vg=pjoin(PATH_OUTPUT, "{name_msa}.vg")
    log:
        # out=pjoin(PATH_OUTPUT,"logs","{name_msa}-rule-vg_graph.out.log"),
        err=pjoin(PATH_OUTPUT,"logs","{name_msa}-rule-vg_graph.err.log"),
    conda:
        "../envs/pggb.yaml"
    shell:
        "/usr/bin/time -v vg construct --msa {input.path_msa} 2> {log.err} > {output.path_vg}"

rule vg_to_gfa:
    input:
        path_vg=pjoin(PATH_OUTPUT, "{name_msa}.vg")
    output: 
        path_gfa=temp(pjoin(PATH_OUTPUT, "{name_msa}-aux.gfa"))
    log:
        # out=pjoin(PATH_OUTPUT,"logs","{name_msa}-rule-vg_to_gfa.out.log"),
        err=pjoin(PATH_OUTPUT,"logs","{name_msa}-rule-vg_to_gfa.err.log"),
    conda:
        "../envs/pggb.yaml"
    shell:
        "/usr/bin/time -v vg view --gfa {input.path_vg} 2> {log.err} > {output.path_gfa}"

rule unchop:
    input:
        path_gfa=pjoin(PATH_OUTPUT, "{name_msa}-aux.gfa")
    output:
        path_unchop_gfa=pjoin(PATH_OUTPUT, "{name_msa}.gfa")
    conda:
        "../envs/pggb.yaml"
    shell:
        "vg mod -u {input.path_gfa} > {output.path_unchop_gfa}"