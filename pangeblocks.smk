configfile: "params.yml"
from pathlib import Path
import json
from os.path import join as pjoin

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
LOG_LEVEL = config["LOG_LEVEL"]

ALPHA=config["OPTIMIZATION"]["THRESHOLD_VERTICAL_BLOCKS"]

OBJ_FUNCTIONS=config["OPTIMIZATION"]["OBJECTIVE_FUNCTION"]
print(OBJ_FUNCTIONS)

# 'weighted' and 'depth' loss functions
PENALIZATION=config["OPTIMIZATION"]["PENALIZATION"] 
MIN_LEN=config["OPTIMIZATION"]["MIN_LEN"]
MIN_COVERAGE=config["OPTIMIZATION"]["MIN_COVERAGE"]

# path msas
MSAS = list(Path(PATH_MSAS).glob("*.fa"))
NAMES = [path.stem for path in MSAS]
# SUBSET_HLA = ["A-3105", "B-3106", "C-3107", "DQA1-3117", "DQB1-3119", "DRB1-3123"]
# NAMES = [path.stem for path in MSAS if path.stem in SUBSET_HLA]
print(NAMES)
EXT_MSA = MSAS[0].suffix

Path(PATH_OUTPUT).mkdir(parents=True, exist_ok=True)
config["NAMES"] = NAMES
with open(Path(PATH_OUTPUT).joinpath("config.json"), "w") as fp:
    json.dump(config, fp, indent=1)

def get_graphs(wildcards):
    "Return a list of graphs to be generated based on parameters provided in the config file"
    graphs = []
    if "nodes" in OBJ_FUNCTIONS:
        graphs.extend(
            pjoin(PATH_OUTPUT, "gfa-unchop", "nodes", f"penalization0-min_len0-min_coverage0-alpha{alpha}", f"{name_msa}.gfa")
            for alpha in ALPHA for name_msa in NAMES
        )
            
    if "strings" in OBJ_FUNCTIONS:
        graphs.extend(
            pjoin(PATH_OUTPUT, "gfa-unchop", "strings", f"penalization0-min_len0-min_coverage0-alpha{alpha}", f"{name_msa}.gfa")
            for alpha in ALPHA for name_msa in NAMES
        )

    if "weighted" in OBJ_FUNCTIONS:
        graphs.extend(
            pjoin(PATH_OUTPUT, "gfa-unchop", "weighted", f"penalization{penalization}-min_len{min_len}-min_coverage0-alpha{alpha}", f"{name_msa}.gfa")
            for alpha in ALPHA for penalization in PENALIZATION for min_len in MIN_LEN for name_msa in NAMES
        )

    if "depth" in OBJ_FUNCTIONS:
        graphs.extend(
            pjoin(PATH_OUTPUT, "gfa-unchop", "depth", f"penalization{penalization}-min_len0-min_coverage{min_coverage}-alpha{alpha}", f"{name_msa}.gfa")
            for alpha in ALPHA for penalization in PENALIZATION for min_coverage in MIN_COVERAGE for name_msa in NAMES
        )
    if "depth_and_len" in OBJ_FUNCTIONS:
        graphs.extend(
            pjoin(PATH_OUTPUT, "gfa-unchop", "depth_and_len", f"penalization0-min_len0-min_coverage0-alpha{alpha}", f"{name_msa}.gfa")
            for alpha in ALPHA for name_msa in NAMES
        )
    
    return graphs

rule all:
    input:
        get_graphs,
        "Wild-pBWT/bin/wild-pbwt"


rule install_wild_pbwt:
    params: 
        github = "https://github.com/illoxian/Wild-pBWT.git"
    output:
        "Wild-pBWT/bin/wild-pbwt"
    log:
        pjoin(PATH_OUTPUT, "logs", "rule-install_wild_pbwt.err.log"),
    conda:
        "envs/wild-pbwt.yml"
    shell:
        """
        if ! [ -f "Wild-pBWT/bin/wild-pbwt"]; then
            rm -rf Wild-pBWT/
            git clone {params.github} && cd Wild-pBWT
            make wild-pbwt
        else
            echo "wild-pbwt already installed"
        fi
        """

rule compute_vertical_blocks:
    input:
        msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
    output: 
        pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks_alpha{alpha}.json")
    params:
        log_level=LOG_LEVEL
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{name_msa}-rule-compute_blocks_alpha{alpha}.err.log"),
    conda: 
        "envs/pangeblocks.yml"
    shell:
        """/usr/bin/time --verbose src/greedy_vertical_blocks.py {input.msa} --output {output} \
        --threshold-vertical-blocks {wildcards.alpha} --log-level {params.log_level} > {log.stderr} 2>&1"""

rule submsa_index:
    input:
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vertical_blocks=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks_alpha{alpha}.json"),
    output:
        path_submsa_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}_alpha{alpha}.txt")
    params:
        max_positions=config["MAX_POSITIONS_SUBMSAS"]
    log: 
        stdout=pjoin(PATH_OUTPUT, "logs", "{name_msa}-alpha{alpha}-rule-submsa_index.out.log"),
    conda: 
        "envs/pangeblocks.yml"
    shell:
        """/usr/bin/time --verbose src/submsas.py --path-msa {input.path_msa} \
        --path-vertical-blocks {input.path_vertical_blocks} --threshold-vertical-blocks {wildcards.alpha} \
        --max-positions-submsas {params.max_positions} --output {output} > {log.stdout} 2>&1
        """

rule ilp:
    input:
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_submsas_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}_alpha{alpha}.txt"),
        bin_wildpbwt="Wild-pBWT/bin/wild-pbwt"
    output:
        auxfile=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log")
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{obj_func}","penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}"),
        log_level=config["LOG_LEVEL"],
        time_limit=config["OPTIMIZATION"]["TIME_LIMIT"],
        threads_ilp=config["THREADS"]["ILP"],
        workers=config["THREADS"]["SUBMSAS"],
        use_wildpbwt=config["USE_WILDPBWT"],
        standard_decomposition=config["DECOMPOSITION"]["STANDARD"],
        alpha_consistent=config["DECOMPOSITION"]["ALPHA_CONSISTENT"],
        min_nrows_fix_block=config["MIN_ROWS_FIX_BLOCK"],
        min_ncols_fix_block=config["MIN_COLS_FIX_BLOCK"]
    threads:
        config["THREADS"]["TOTAL"]
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log"),
    conda: 
        "envs/pangeblocks.yml"
    resources: 
        mem_mb=config["RESOURCES"]["MEM_MB"]
    shell:
        """
        /usr/bin/time --verbose src/exact_cover.py --path-msa {input.path_msa} --obj-function {wildcards.obj_func} \
        --prefix-output {params.dir_subsols}/{wildcards.name_msa} \
        --penalization {wildcards.penalization} --min-len {wildcards.min_len} --min-coverage {wildcards.min_coverage} \
        --submsa-index {input.path_submsas_index} --time-limit {params.time_limit} --solve-ilp True \
        --use-wildpbwt {params.use_wildpbwt} --bin-wildpbwt {input.bin_wildpbwt} \
        --standard-decomposition {params.standard_decomposition} --threads-ilp {params.threads_ilp} \
        --workers {params.workers} --alpha-consistent {params.alpha_consistent} \
        --min-rows-fixblock {params.min_nrows_fix_block} --min-columns-fixblock {params.min_ncols_fix_block} > {output.auxfile} 2> {log.stderr}
        """

rule coverage_to_graph:
    input:
        auxfile=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-ilp.log"),
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vb=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks_alpha{alpha}.json")
    output:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{obj_func}","penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}"),
    log:
        stdout=pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-coverage_to_graph.log")
    conda: 
        "envs/pangeblocks.yml"
    shell:
        """/usr/bin/time --verbose src/compute_gfa.py --path-msa {input.path_msa} \
        --dir-subsolutions {params.dir_subsols} --path-vert-blocks {input.path_vb} \
        --path-gfa {output} > {log} 2>&1"""

rule postprocessing_gfa:
    input:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    output:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    conda: 
        "envs/pangeblocks.yml"
    shell:
        "/usr/bin/time --verbose python src/postprocess_gfa.py --path_gfa {input.path_gfa} > {output.path_post_gfa}" 

rule unchop_gfa:
    input:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa")
    output:
        path_unchop_gfa=pjoin(PATH_OUTPUT, "gfa-unchop","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.gfa"),
        path_labels=pjoin(PATH_OUTPUT, "gfa-unchop","{obj_func}", "penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}", "{name_msa}.csv")
    log:
        pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-alpha{alpha}-rule-unchop_gfa.log")
    conda:
        "envs/vg.yml"
    shell:
        """
        mkdir -p "$(dirname "{output.path_unchop_gfa}")" 
        /usr/bin/time --verbose vg mod -u {input} > {output.path_unchop_gfa} 2> {log}
        src/graph/bandage_labels_from_gfa.py --path-gfa {output.path_unchop_gfa} --path-save {output.path_labels}
        """