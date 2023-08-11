configfile: "test/params.yaml"
from pathlib import Path
import pandas as pd
from os.path import join as pjoin

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
LOG_LEVEL = config["LOG_LEVEL"]

# path msas
MSAS = list(Path(PATH_MSAS).glob("*.[fa]*"))
NAMES = [path.stem for path in MSAS]
EXT_MSA = MSAS[0].suffix

# --- 

rule all: 
    input:
        expand(pjoin(PATH_OUTPUT, "gfa-unchop","{name_msa}.gfa"), name_msa=NAMES)

rule compute_vertical_blocks:
    input: 
        pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA)
    output: 
        pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json")
    params:
        log_level=LOG_LEVEL,
        threshold_vertical_blocks=config["THRESHOLD_VERTICAL_BLOCKS"]
    log:
        err=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.log")
    shell:
        """/usr/bin/time --verbose src/greedy_vertical_blocks.py {input} --output {output} \
        --threshold-vertical-blocks {params.threshold_vertical_blocks} --log-level {params.log_level} > {log.err} 2>&1"""

rule submsa_index:
    input: 
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vertical_blocks=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json"),
    output:
        path_submsa_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}.txt")
    params:
        threshold_vertical_blocks=config["THRESHOLD_VERTICAL_BLOCKS"]
    log:
        out=pjoin(PATH_OUTPUT, "submsas", "{name_msa}.log")
    shell:
        """/usr/bin/time --verbose src/submsas.py --path-msa {input.path_msa} \
        --path-vertical-blocks {input.path_vertical_blocks} --threshold-vertical-blocks {params} --output {output} > {log.out} 2>&1
        """

rule ilp:
    input: 
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_submsas_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}.txt")
    output:
        auxfile=pjoin(PATH_OUTPUT, "ilp","{name_msa}","{name_msa}.log")
    params:
        path_save=pjoin(PATH_OUTPUT, "ilp","{name_msa}"),
        obj_function=config["OPTIMIZATION"]["OBJECTIVE_FUNCTION"],
        penalization=config["OPTIMIZATION"]["PENALIZATION"],
        min_len=config["OPTIMIZATION"]["MIN_LEN"],
        min_coverage=config["OPTIMIZATION"]["MIN_COVERAGE"],
        log_level=config["LOG_LEVEL"],
        time_limit=config["OPTIMIZATION"]["TIME_LIMIT"]
    threads:
        config["THREADS"]
    log:
        err=pjoin(PATH_OUTPUT, "ilp","{name_msa}","{name_msa}.err.log")
    shell:
        """/usr/bin/time --verbose src/solve_submsa.py --path-msa {input.path_msa} --obj-function {params.obj_function} \
        --path-save-ilp {params.path_save}/{wildcards.name_msa} --path-opt-solution {params.path_save}/{wildcards.name_msa} \
        --penalization {params.penalization} --min-len {params.min_len} --min-coverage {params.min_coverage} \
        --submsa-index {input.path_submsas_index} --time-limit {params.time_limit} --solve-ilp true \
        --workers {threads} > {output.auxfile} 2> {log.err}"""

rule coverage_to_graph:
    input:
        auxfile=pjoin(PATH_OUTPUT, "ilp","{name_msa}","{name_msa}.log"),
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vb=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks.json")
    output:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{name_msa}.gfa"),
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp","{name_msa}")
    log:
        auxfile=pjoin(PATH_OUTPUT, "gfa","{name_msa}.log"),
    shell:
        """/usr/bin/time --verbose src/compute_gfa.py --path-msa {input.path_msa} \
        --dir-subsolutions {params.dir_subsols} --path-vert-blocks {input.path_vb} \
        --path-gfa {output} > {log} 2>&1"""

rule postprocessing_gfa:
    input:
        path_gfa=pjoin(PATH_OUTPUT, "gfa", "{name_msa}.gfa")
    output:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post", "{name_msa}.gfa") 
    shell:
        "/usr/bin/time --verbose python src/postprocess_gfa.py --path_gfa {input.path_gfa} > {output.path_post_gfa}" 

rule unchop_gfa:
    input:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post", "{name_msa}.gfa")
    output:
        path_unchop_gfa=pjoin(PATH_OUTPUT, "gfa-unchop", "{name_msa}.gfa"),
        path_labels=pjoin(PATH_OUTPUT, "gfa-unchop", "{name_msa}.csv")
    conda:
        "envs/vg.yaml"
    shell:
        """
        /usr/bin/time --verbose vg mod -u {input} > {output.path_unchop_gfa} 
        python src/graph/bandage_labels_from_gfa.py --path_gfa {output.path_unchop_gfa} --path_save {output.path_labels}
        """