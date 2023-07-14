configfile: "params-grid-exp.yaml"
from pathlib import Path
import pandas as pd
from os.path import join as pjoin

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
LOG_LEVEL = config["LOG_LEVEL"]

ALPHA=config["OPTIMIZATION"]["THRESHOLD_VERTICAL_BLOCKS"]

# 'weighted' and 'depth' loss functions
PENALIZATION=config["OPTIMIZATION"]["PENALIZATION"] 
MIN_LEN=config["OPTIMIZATION"]["MIN_LEN"]
MIN_COVERAGE=config["OPTIMIZATION"]["MIN_COVERAGE"]

# path msas
SUBSET_HLA = ["A-3105", "B-3106", "C-3107", "DQA1-3117", "DQB1-3119", "DRB1-3123"]
MSAS = list(Path(PATH_MSAS).glob("*.[fa]*"))
NAMES = [path.stem for path in MSAS]
# NAMES = [path.stem for path in MSAS if path.stem in SUBSET_HLA]
EXT_MSA = MSAS[0].suffix
print(NAMES)
rule all:
    input:
        expand(pjoin(PATH_OUTPUT, "gfa-unchop", "{obj_func}", "penalization0-min_len0-min_coverage0-alpha{alpha}", "{name_msa}.gfa"), 
        obj_func=["strings","nodes"], name_msa=NAMES, alpha=ALPHA),
        expand(pjoin(PATH_OUTPUT, "gfa-unchop", "weighted", "penalization{penalization}-min_len{min_len}-min_coverage0-alpha{alpha}" ,"{name_msa}.gfa"), 
        penalization=PENALIZATION, min_len=MIN_LEN, min_coverage=MIN_COVERAGE, alpha=ALPHA, name_msa=NAMES),
        expand(pjoin(PATH_OUTPUT, "gfa-unchop", "depth", "penalization{penalization}-min_len0-min_coverage{min_coverage}-alpha{alpha}" ,"{name_msa}.gfa"),
        penalization=PENALIZATION, min_coverage=MIN_COVERAGE, alpha=ALPHA, name_msa=NAMES)

rule compute_vertical_blocks:
    input:
        pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA)
    output: 
        pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks_alpha{alpha}.json")
    params:
        log_level=LOG_LEVEL
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{name_msa}-rule-compute_blocks_alpha{alpha}.err.log"),
    shell:
        """/usr/bin/time --verbose src/greedy_vertical_blocks.py {input} --output {output} \
        --threshold-vertical-blocks {wildcards.alpha} --log-level {params.log_level} > {log.stderr} 2>&1"""

rule submsa_index:
    input:
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vertical_blocks=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks_alpha{alpha}.json"),
    output:
        path_submsa_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}_alpha{alpha}.txt")
    log: 
        stdout=pjoin(PATH_OUTPUT, "logs", "{name_msa}-alpha{alpha}-rule-submsa_index.out.log"),
    shell:
        """/usr/bin/time --verbose src/submsas.py --path-msa {input.path_msa} \
        --path-vertical-blocks {input.path_vertical_blocks} --threshold-vertical-blocks {wildcards.alpha} --output {output} > {log.stdout} 2>&1
        """

rule ilp:
    input:
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_submsas_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}_alpha{alpha}.txt")
    output:
        auxfile=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log")
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "alpha{alpha}"),
        log_level=config["LOG_LEVEL"],
        time_limit=config["OPTIMIZATION"]["TIME_LIMIT"],
        threads_ilp=config["THREADS"]["ILP"],
        use_wildpbwt=config["USE_WILDPBWT"],
        bin_wildpbwt=config["BIN_WILDPBWT"],

    threads:
        config["THREADS"]["TOTAL"]
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log"),
    shell:
        """/usr/bin/time --verbose src/solve_submsa.py --path-msa {input.path_msa} --obj-function {wildcards.obj_func} \
        --path-save-ilp {params.dir_subsols}/{wildcards.name_msa} --path-opt-solution {params.dir_subsols}/{wildcards.name_msa} \
        --penalization {wildcards.penalization} --min-len {wildcards.min_len} --min-coverage {wildcards.min_coverage} \
        --submsa-index {input.path_submsas_index} --time-limit {params.time_limit} --solve-ilp true \
        --use-wildpbwt {params.use_wildpbwt} --bin-wildpbwt {params.bin_wildpbwt}\
        --workers {threads} > {output.auxfile} 2> {log.stderr}"""

rule coverage_to_graph:
    input:
        auxfile=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log"),
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vb=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks_alpha{alpha}.json")
    output:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "alpha{alpha}"),
    log:
        stdout=pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-coverage_to_graph.log")
    shell:
        """/usr/bin/time --verbose src/compute_gfa.py --path-msa {input.path_msa} \
        --dir-subsolutions {params.dir_subsols} --path-vert-blocks {input.path_vb} \
        --path-gfa {output} > {log} 2>&1"""

rule postprocessing_gfa:
    input:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    output:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    shell:
        "/usr/bin/time --verbose python src/postprocess_gfa.py --path_gfa {input.path_gfa} > {output.path_post_gfa}" 

rule unchop_gfa:
    input:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    output:
        path_unchop_gfa=pjoin(PATH_OUTPUT, "gfa-unchop","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa"),
        path_labels=pjoin(PATH_OUTPUT, "gfa-unchop","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.csv")
    conda:
        "envs/vg.yaml"
    shell:
        """
        /usr/bin/time --verbose vg mod -u {input} > {output.path_unchop_gfa} 
        python src/graph/bandage_labels_from_gfa.py --path_gfa {output.path_unchop_gfa} --path_save {output.path_labels}
        """