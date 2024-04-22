configfile: "rules/params/mafft.yaml"
from pathlib import Path
from os.path import join as pjoin

PATH_SEQS=config["PATH_SEQS"]
OP=config["OP"]
EP=config["EP"]

PATH_OUTPUT=config["MAFFT_OUTPUT"]
Path(PATH_OUTPUT).mkdir(exist_ok=True, parents=True)

list_fasta = list(Path(PATH_SEQS).glob("*.fa"))# + list(Path(PATH_SEQS).glob("*.fa"))
NAMES = [path.stem for path in list_fasta]
print(NAMES)

rule all:
    input:
        expand(pjoin(PATH_OUTPUT, "mafft.op{op}-ep{ep}" ,"{name_msa}.fa"), op=OP, ep=EP, name_msa=NAMES)

rule create_msa:
    input: 
        pjoin(PATH_SEQS, "{name_msa}.fa")
    output: 
        pjoin(PATH_OUTPUT, "mafft.op{op}-ep{ep}", "{name_msa}.fa")
    params:
        op=OP,
        ep=EP
    log:
        err=pjoin(PATH_OUTPUT, "mafft.op{op}-ep{ep}", "logs", "{name_msa}.err.log")
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        /usr/bin/time -v mafft --op {params.op} --ep {params.ep} {input} > {output} 2> {log.err}
        """