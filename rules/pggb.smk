configfile: "rules/params/pggb.yaml"
from pathlib import Path
from os.path import join as pjoin

PATH_INPUT=config["PGGB_INPUT"]
PATH_OUTPUT=config["PGGB_OUTPUT"]
Path(PATH_OUTPUT).mkdir(exist_ok=True, parents=True)

list_fasta = list(Path(PATH_INPUT).glob("*.fa"))# + list(Path(PATH_INPUT).glob("*.fasta")) 
NAMES = [path.stem for path in list_fasta]
print(NAMES)

rule all:
    input:
        expand( pjoin(PATH_OUTPUT,"{name_input}.gfa"), name_input=NAMES)

rule index: 
    input: 
        pjoin(PATH_INPUT, "{name_input}.fa")
    output:
        pjoin(PATH_INPUT, "{name_input}.fa.gz"),
        pjoin(PATH_INPUT, "{name_input}.fa.gz.fai")
    # resources:
    #     mem_mb=200000
    shell:
        """
        f={input}
        /usr/bin/time bgzip -@ 16 -c $f > $f.gz && samtools faidx $f.gz
        """

rule pggb: 
    input: 
        pjoin(PATH_INPUT, "{name_input}.fa"),
        pjoin(PATH_INPUT, "{name_input}.fa.gz")
    output: 
        path_gfa=pjoin(PATH_OUTPUT, "{name_input}.gfa")
    log:
        out=pjoin(PATH_OUTPUT,"logs","{name_input}-rule-pggb.out.log"),
        err=pjoin(PATH_OUTPUT,"logs","{name_input}-rule-pggb.err.log"),
    params:
        path_output=PATH_OUTPUT
    resources:
        # mem_mb=200000
        limit_space=1,
    conda:
        "../envs/pggb.yml"
    shell:
        """
        f={input[0]}
        /usr/bin/time -v pggb -i $f.gz -t 16 -n $(cat $f.gz.fai | wc -l) -o {params.path_output}/$(basename $f .fa.gz) 2> {log.err} 
        """


