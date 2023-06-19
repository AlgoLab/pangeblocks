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
        # expand(pjoin(PATH_OUTPUT, "submsas", "{name_msa}.txt"), name_msa=NAMES)
        # expand(pjoin(PATH_OUTPUT, "ilp","{name_msa}","{name_msa}.log"), name_msa=NAMES),
        # expand(pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json"), name_msa=NAMES)
        expand(pjoin(PATH_OUTPUT, "gfa","{name_msa}.gfa"), name_msa=NAMES)

rule compute_vertical_blocks:
    input: 
        pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA)
    output: 
        pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json")
    params:
        log_level=LOG_LEVEL
    shell:
        """/usr/bin/time --verbose src/maximal_blocks.py {input} --output {output} \
        --only-vertical-blocks true --log-level {params.log_level}""" 

rule submsa_index:
    input: 
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vertical_blocks=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json"),
    output: 
        submsa_index=pjoin(PATH_OUTPUT, "submsas", "{name_msa}.txt")
    shell:
        """/usr/bin/time --verbose src/submsas.py --path-msa {input.path_msa} \
        --path-vertical-blocks {input.path_vertical_blocks} --output {output}
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
    shell:
        """/usr/bin/time --verbose src/solve_submsa.py --path-msa {input.path_msa} --obj-function {params.obj_function} \
        --path-save-ilp {params.path_save}/{wildcards.name_msa} --path-opt-solution {params.path_save}/{wildcards.name_msa} \
        --penalization {params.penalization} --min-len {params.min_len} --min-coverage {params.min_coverage} \
        --submsa-index {input.path_submsas_index} --time-limit {params.time_limit} --solve-ilp true \
        --workers {threads} > {output.auxfile}"""

# TODO
rule coverage_to_graph:
    input:
        auxfile=pjoin(PATH_OUTPUT, "ilp","{name_msa}","{name_msa}.log"),
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        path_vb=pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}","vertical_blocks.json")
    output:
        path_gfa=pjoin(PATH_OUTPUT, "gfa","{name_msa}.gfa"),
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp","{name_msa}")
    shell:
        """/usr/bin/time --verbose src/compute_gfa.py --path-msa {input.path_msa} \
        --dir-subsolutions {params.dir_subsols} --path-vert-blocks {input.path_vb} --path-gfa {output}"""

# rule compute_blocks:
#     input:
#         pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
#         pjoin(PATH_OUTPUT, "submsas", "{name_msa}.txt")
#         # pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "vertical_blocks.json")
#     output: 
#         pjoin(PATH_OUTPUT, "maximal-blocks", "{name_msa}", "maximal_blocks_submsa_{start_column}-{end_column}.json")
#     params:
#         start_column=0,
#         end_column=-1,
#         log_level=LOG_LEVEL
#     shell:
#         """/usr/bin/time --verbose src/maximal_blocks.py {input} --output {output} 
#         --start-column {params.start_column} --end-column {params.end_column} 
#         --only-vertical-blocks false --log-level {params.log_level}"""

# rule all:
#     input:
#         expand("{path_output}/max_blocks/stats/{name_msa}.tsv", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/block_decomposition/{name_msa}.json", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/gfa/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/coverage/{name_msa}-gray.jpg", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/coverage/{name_msa}-color.jpg", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/gfa/{name_msa}.csv", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/gfa-post/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT),
#         expand("{path_output}/gfa-unchop/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT)

# rule compute_blocks:
#     input:
#         expand("{path_msas}/{{name_msa}}{ext_msa}", path_msas=PATH_MSAS, ext_msa=EXT_MSA)
#     output: 
#         expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
#     log: 
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-compute_blocks.out.log', path_output=PATH_OUTPUT),
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-compute_blocks.err.log', path_output=PATH_OUTPUT)
#     shell:
#         "/usr/bin/time --verbose python src/compute_blocks.py {input} --output {output} 2> {log.stderr} > {log.stdout}"

# rule analyze_blocks:
#     input:
#         expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
#     output:
#         expand("{path_output}/max_blocks/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
#     log: 
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-analyze_blocks.out.log', path_output=PATH_OUTPUT),
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-analyze_blocks.err.log', path_output=PATH_OUTPUT)
#     shell: 
#         "/usr/bin/time --verbose python analyze_blocks.py {input} --output {output} 2> {log.stderr} > {log.stdout}"

# rule decompose_blocks:
#     input:
#         expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
#     output:
#         expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
#         expand("{path_output}/block_decomposition/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
#     log: 
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-decompose_blocks.out.log', path_output=PATH_OUTPUT),
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-decompose_blocks.err.log', path_output=PATH_OUTPUT)
#     shell:
#         "/usr/bin/time --verbose python decompose_blocks.py {input} --output {output[0]} --output-stats {output[1]} 2> {log.stderr} > {log.stdout}"

# rule pangeblock:
#     input:
#         path_blocks=expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
#         path_msa=expand("{path_msas}/{{name_msa}}{ext_msa}", path_msas=PATH_MSAS, ext_msa=EXT_MSA)
#     output: 
#         path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT),
#         path_oc=expand("{path_output}/opt-coverage/{{name_msa}}.json", path_output=PATH_OUTPUT)
#     log: 
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-pangeblock.err.log', path_output=PATH_OUTPUT),
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-pangeblock.out.log', path_output=PATH_OUTPUT)        
#     params: 
#         obj_function=config["OPTIMIZATION"]["OBJECTIVE_FUNCTION"],
#         penalization=config["OPTIMIZATION"]["PENALIZATION"],
#         min_len=config["OPTIMIZATION"]["MIN_LEN"],
#         log_level=config["LOG_LEVEL"],
#         time_limit=config["OPTIMIZATION"]["TIME_LIMIT"]
#     shell: 
#         """/usr/bin/time --verbose python compute_gfa.py --path_blocks {input.path_blocks} \
#         --path_msa {input.path_msa} --path_gfa {output.path_gfa} --path_oc {output.path_oc} \
#         --obj_function {params.obj_function} --penalization {params.penalization} --min_len {params.min_len} \
#         --time_limit {params.time_limit} \
#         --log_level {params.log_level} > {log.stdout} 2> {log.stderr}
#         """

# rule bandage_labels:
#     input: 
#         path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
#     output: 
#         path_labels=expand("{path_output}/gfa/{{name_msa}}.csv", path_output=PATH_OUTPUT)
#     shell:
#         "python src/graph/bandage_labels_from_gfa.py --path_gfa {input} --path_save {output}"

# rule coverage:
#     input:
#         path_blocks=expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
#         path_msa=expand("{path_msa}/{{name_msa}}{ext_msa}", path_msa=PATH_MSAS, ext_msa=EXT_MSA)
#     output:
#         path_gray=expand("{path_output}/coverage/{{name_msa}}-gray.jpg", path_output=PATH_OUTPUT),
#         path_color=expand("{path_output}/coverage/{{name_msa}}-color.jpg", path_output=PATH_OUTPUT)
#     log: 
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-coverage.out.log', path_output=PATH_OUTPUT),
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-coverage.err.log', path_output=PATH_OUTPUT)
#     shell:
#         "/usr/bin/time --verbose python coverage.py --path_blocks {input[0]} --path_msa {input[1]} --path_gray {output[0]} --path_color {output[1]} 2> {log.stderr} > {log.stdout}"

# rule postprocessing_gfa:
#     input:
#         path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
#     output:
#         path_post_gfa=expand("{path_output}/gfa-post/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
#     log:
#         stderr=expand('{path_output}/logs/{{name_msa}}-rule-postprocessing_gfa.err.log', path_output=PATH_OUTPUT),
#         stdout=expand('{path_output}/logs/{{name_msa}}-rule-postprocessing_gfa.out.log', path_output=PATH_OUTPUT)    
#     shell:
#         "/usr/bin/time --verbose python src/postprocess_gfa.py --path_gfa {input.path_gfa} 2> {log.stderr} > {output.path_post_gfa}"

# rule unchop_gfa:
#     input:
#         path_post_gfa=expand("{path_output}/gfa-post/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
#     output:
#         path_unchop_gfa=expand("{path_output}/gfa-unchop/{{name_msa}}.gfa", path_output=PATH_OUTPUT),
#         path_labels=expand("{path_output}/gfa-unchop/{{name_msa}}.csv", path_output=PATH_OUTPUT)
#     # log:
#     #     stdout=expand('{path_output}/logs/{{name_msa}}-rule-unchop_gfa.out.log', path_output=PATH_OUTPUT),
#     #     stderr=expand('{path_output}/logs/{{name_msa}}-rule-unchop_gfa.err.log', path_output=PATH_OUTPUT)
#     shell:
#         """
#         ../vg mod -u {input} > {output.path_unchop_gfa}
#         python src/graph/bandage_labels_from_gfa.py --path_gfa {output.path_unchop_gfa} --path_save {output.path_labels}
#         """