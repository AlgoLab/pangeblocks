from pathlib import Path
import json
from os.path import join as pjoin
from os import getcwd
import pprint

config['gfa'] = Path(config["PATH_OUTPUT"])
config['PATH_OUTPUT'] = config['gfa'].parent
config['msa'] = Path(config["PATH_MSA"])
config['alpha']=config["THRESHOLD_VERTICAL_BLOCKS"]

# path msas
config['tmpdir']=pjoin(config["PATH_OUTPUT"], "tmp")
config['logdir']=pjoin(config["PATH_OUTPUT"], "logs")
Path(config['tmpdir']).mkdir(parents=True, exist_ok=True)
Path(config['logdir']).mkdir(parents=True, exist_ok=True)
pprint.pprint(config)


rule all:
    input:
        gfa=config['gfa'],

rule compute_vertical_blocks:
    input:
        msa=config['msa'],
    output: 
        pjoin(config['tmpdir'], "maximal-blocks", "vertical_blocks.json")
    params:
        root_dir=config['root_dir'],
        alpha=config['alpha'],
        log_level=config['LOG_LEVEL']
    log:
        stderr=pjoin(config['logdir'], "rule-compute_vertical_blocks.log"),
    shell:
        """/usr/bin/time --verbose {params.root_dir}/src/greedy_vertical_blocks.py {input.msa} --output {output} \
        --threshold-vertical-blocks {params.alpha} --log-level {params.log_level} > {log.stderr} 2>&1"""

rule submsa_index:
    input:
        path_msa=config['msa'],
        path_vertical_blocks=pjoin(config['tmpdir'], "maximal-blocks", "vertical_blocks.json"),
    output:
        path_submsa_index=pjoin(config['tmpdir'], "submsas.txt")
    params:
        root_dir=config['root_dir'],
        alpha=config['alpha'],
        max_positions=config["MAX_POSITIONS_SUBMSAS"]
    log: 
        stdout=pjoin(config['logdir'], "rule-submsa_index.log"),
    shell:
        """/usr/bin/time --verbose {params.root_dir}/src/submsas.py --path-msa {input.path_msa} \
        --path-vertical-blocks {input.path_vertical_blocks} --threshold-vertical-blocks {params.alpha} \
        --max-positions-submsas {params.max_positions} --output {output} > {log.stdout} 2>&1
        """

rule ilp:
    input:
        path_msa=config['msa'],
        path_submsas_index=pjoin(config['tmpdir'], "submsas.txt")
    output:
        auxfile=pjoin(config['tmpdir'], "ilp", "rule-ilp.log")
    params:
        root_dir=config['root_dir'],
        dir_subsols=pjoin(config['tmpdir'], "ilp", "subsol"),
        log_level=config["LOG_LEVEL"],
        time_limit=config["TIME_LIMIT"],
        threads_ilp=config["ILP"],
        workers=config["SUBMSAS"],
        standard_decomposition=config["STANDARD"],
        alpha_consistent=config["ALPHA_CONSISTENT"],
        min_nrows_fix_block=config["MIN_ROWS_FIX_BLOCK"],
        obj_function=config['OBJECTIVE_FUNCTION'],
        penalization=config['PENALIZATION'],
        min_len=config['MIN_LEN'],
        min_coverage=config['MIN_COVERAGE'],
        min_ncols_fix_block=config["MIN_COLS_FIX_BLOCK"]
    log:
        stderr=pjoin(config['logdir'], "rule-ilp.log"),
    resources:
        mem_mb=config["MEM_MB"]
    shell:
        """
        /usr/bin/time --verbose {params.root_dir}/src/exact_cover.py --path-msa {input.path_msa} --obj-function {params.obj_function} \
        --prefix-output {params.dir_subsols}/ilp-output \
        --penalization {params.penalization} --min-len {params.min_len} --min-coverage {params.min_coverage} \
        --submsa-index {input.path_submsas_index} --time-limit {params.time_limit} --solve-ilp True \
        --use-wildpbwt True --bin-wildpbwt "{params.root_dir}/lib/Wild-pBWT/bin/wild-pbwt" \
        --standard-decomposition {params.standard_decomposition} --threads-ilp {params.threads_ilp} \
        --workers {params.workers} --alpha-consistent {params.alpha_consistent} \
        --min-rows-fixblock {params.min_nrows_fix_block} --min-columns-fixblock {params.min_ncols_fix_block} > {output.auxfile} 2> {log.stderr}
        """

rule coverage_to_graph:
    input:
        auxfile=pjoin(config['tmpdir'], "ilp", "rule-ilp.log"),
        path_vb=pjoin(config['tmpdir'], "maximal-blocks", "vertical_blocks.json"),
        msa=config['msa'],
    output:
        path_gfa=pjoin(config['tmpdir'], "gfa", "raw.gfa")
    params:
        root_dir=config['root_dir'],
        dir_subsols=pjoin(config['tmpdir'], "ilp"),
    log:
        stdout=pjoin(config['logdir'], "coverage_to_graph.log")
    shell:
        """/usr/bin/time --verbose {params.root_dir}/src/compute_gfa.py --path-msa {input.msa} \
        --dir-subsolutions {params.dir_subsols} --path-vert-blocks {input.path_vb} \
        --path-gfa {output} > {log} 2>&1"""

rule postprocessing_gfa:
    input:
        path_gfa=pjoin(config['tmpdir'], "gfa", "raw.gfa"),
    params:
            root_dir=config['root_dir'],
    output:
        path_post_gfa=pjoin(config['tmpdir'], "gfa-post", "chopped.gfa")
    shell:
        "/usr/bin/time --verbose python {params.root_dir}/src/postprocess_gfa.py --path_gfa {input.path_gfa} > {output.path_post_gfa}"

rule unchop_gfa:
    input:
        path_post_gfa=pjoin(config['tmpdir'], "gfa-post", "chopped.gfa"),
    output:
        path_unchop_gfa=config['gfa'],
        path_labels=f"{config['gfa']}.csv"
    params:
        root_dir=config['root_dir'],
    log:
        pjoin(config['logdir'], "unchop_gfa.log")
    shell:
        """
        /usr/bin/time --verbose vg mod -u {input} > {output.path_unchop_gfa} 2> {log}
        {params.root_dir}/src/graph/bandage_labels_from_gfa.py --path-gfa {output.path_unchop_gfa} --path-save {output.path_labels}
        """
