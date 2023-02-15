configfile: "params-grid-exp.yaml"
from pathlib import Path
import pandas as pd
from os.path import join as pjoin

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
LOG_LEVEL = config["LOG_LEVEL"]

# weighted
PENALIZATION=config["OPTIMIZATION"]["PENALIZATION"]
MIN_LEN=config["OPTIMIZATION"]["MIN_LEN"]


# load names of MSAs
STATS_MSAS = pd.read_csv(
                Path(PATH_OUTPUT).joinpath("analysis-msa/stats_msas.tsv"), 
                index_col=False, sep="\t"
                )

NAMES = STATS_MSAS[["path_msa","n_seqs"]].query("n_seqs>1")["path_msa"].apply(lambda path: Path(path).stem)
EXT_MSA = Path(STATS_MSAS["path_msa"][0]).suffix 
# --- 

rule all:
    input:
        expand(pjoin(PATH_OUTPUT, "gfa-unchop", "{obj_func}", "penalization0-min_len0", "{name_msa}.gfa"), obj_func=["strings","nodes"], name_msa=NAMES),
        expand(pjoin(PATH_OUTPUT, "gfa-unchop", "weighted", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa"), penalization=PENALIZATION, min_len=MIN_LEN, name_msa=NAMES)

rule compute_blocks:
    input:
        pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA)
    output: 
        pjoin(PATH_OUTPUT, "max_blocks", "{name_msa}.json")
    log: 
        stderr=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-compute_blocks.err.log"),
        stdout=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-compute_blocks.out.log")
    shell:
        "/usr/bin/time --verbose python src/compute_blocks.py {input} --output {output} 2> {log.stderr} > {log.stdout}"

rule analyze_blocks:
    input:
        pjoin(PATH_OUTPUT, "max_blocks", "{name_msa}.json")
    output:
        pjoin(PATH_OUTPUT, "max_blocks", "stats", "{name_msa}.tsv")
    log: 
        stderr=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-analyze_blocks.err.log"),
        stdout=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-analyze_blocks.out.log")
    shell: 
        "/usr/bin/time --verbose python analyze_blocks.py {input} --output {output}" # 2> {log.stderr} > {log.stdout}"

rule decompose_blocks:
    input:
        path_max_blocks=pjoin(PATH_OUTPUT, "max_blocks", "{name_msa}.json")
    output:
        path_blocks=pjoin(PATH_OUTPUT, "block_decomposition", "{name_msa}.json"),
        path_blocks_stats=pjoin(PATH_OUTPUT, "block_decomposition", "stats", "{name_msa}.tsv")
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-decompose_blocks.err.log"),
        stdout=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-decompose_blocks.out.log")
    shell:
        "/usr/bin/time --verbose python decompose_blocks.py {input.path_max_blocks} --output {output.path_blocks} --output-stats {output.path_blocks_stats} 2> {log.stderr} > {log.stdout}"

rule pangeblock:
    input:
        path_blocks=pjoin(PATH_OUTPUT, "block_decomposition", "{name_msa}.json"),
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA)
    output: 
        path_gfa=pjoin(PATH_OUTPUT, "gfa", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa"),
        path_oc=pjoin(PATH_OUTPUT, "opt-coverage","{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}.json"),
        path_labels=pjoin(PATH_OUTPUT, "gfa", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.csv")
    log: 
        stderr=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-pangeblock.err.log"),
        stdout=pjoin(PATH_OUTPUT, "logs", "{obj_func}","penalization{penalization}-min_len{min_len}","{name_msa}-rule-pangeblock.out.log")
    params: 
        obj_function= "{obj_func}", 
        penalization= "{penalization}",
        min_len= "{min_len}",
        log_level=config["LOG_LEVEL"],
        time_limit=config["OPTIMIZATION"]["TIME_LIMIT"]
    shell: 
        """
        /usr/bin/time --verbose python compute_gfa.py --path_blocks {input.path_blocks} \
        --path_msa {input.path_msa} --path_gfa {output.path_gfa} --path_oc {output.path_oc} \
        --obj_function {params.obj_function} --penalization {params.penalization} --min_len {params.min_len} \
        --time_limit {params.time_limit} \
        --log_level {params.log_level}  > {log.stdout} 2> {log.stderr}
        python src/graph/bandage_labels_from_gfa.py --path_gfa {input.path_gfa} --path_save {output.path_labels}
        """

rule postprocessing_gfa:
    input:
        path_gfa=pjoin(PATH_OUTPUT, "gfa", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa")
    output:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa") 
    shell:
        "/usr/bin/time --verbose python src/postprocess_gfa.py --path_gfa {input.path_gfa} > {output.path_post_gfa}" 

rule unchop_gfa:
    input:
        path_post_gfa=pjoin(PATH_OUTPUT, "gfa-post", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa")
    output:
        path_unchop_gfa=pjoin(PATH_OUTPUT, "gfa-unchop", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.gfa"),
        path_labels=pjoin(PATH_OUTPUT, "gfa-unchop", "{obj_func}", "penalization{penalization}-min_len{min_len}" ,"{name_msa}.csv")
    shell:
        """
        ../vg mod -u {input} > {output.path_unchop_gfa}
        python src/graph/bandage_labels_from_gfa.py --path_gfa {output.path_unchop_gfa} --path_save {output.path_labels}
        """