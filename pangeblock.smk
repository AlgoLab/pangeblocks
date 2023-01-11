configfile: "params.yaml"
from pathlib import Path
import pandas as pd

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]

# load names of MSAs
STATS_MSAS = pd.read_csv(
                Path(PATH_OUTPUT).joinpath("analysis-msa/stats_msas.tsv"), 
                index_col=False, sep="\t"
                )
NAMES = STATS_MSAS["path_msa"].apply(lambda path: Path(path).stem) # txt file with names of the MSAs

# --- 

rule all:
    input:
        expand("{path_output}/max_blocks/stats/{name_msa}.tsv", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/block_decomposition/{name_msa}.json", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/gfa/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/coverage/{name_msa}-gray.jpg", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/coverage/{name_msa}-color.jpg", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/gfa/{name_msa}.csv", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/gfa-post/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT),
        expand("{path_output}/gfa-unchop/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT)

rule compute_blocks:
    input:
        expand("{path_msas}/{{name_msa}}.fa", path_msas=PATH_MSAS)
    output: 
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    shell:
        "python src/compute_blocks.py {input} --output {output}"

rule analyze_blocks:
    input:
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    output:
        expand("{path_output}/max_blocks/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
    shell: 
        "python analyze_blocks.py {input} --output {output}"

rule decompose_blocks:
    input:
        expand("{path_output}/max_blocks/{{name_msa}}.json", path_output=PATH_OUTPUT)
    output:
        expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
        expand("{path_output}/block_decomposition/stats/{{name_msa}}.tsv", path_output=PATH_OUTPUT)
    shell:
        "python decompose_blocks.py {input} --output {output[0]} --output-stats {output[1]}"

rule pangeblock:
    input:
        path_blocks=expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
        path_msa=expand("{path_msas}/{{name_msa}}.fa", path_msas=PATH_MSAS)
    output: 
        path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    shell: 
        "python compute_gfa.py --path_blocks {input[0]} --path_msa {input[1]} --path_gfa {output[0]}"

rule bandage_labels:
    input: 
        path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    output: 
        path_labels=expand("{path_output}/gfa/{{name_msa}}.csv", path_output=PATH_OUTPUT)
    shell:
        "python src/graph/bandage_labels_from_gfa.py --path_gfa {input} --path_save {output}"

rule coverage:
    input:
        path_blocks=expand("{path_output}/block_decomposition/{{name_msa}}.json", path_output=PATH_OUTPUT),
        path_msa=expand("{path_msa}/{{name_msa}}.fa", path_msa=PATH_MSAS)
    output:
        path_gray=expand("{path_output}/coverage/{{name_msa}}-gray.jpg", path_output=PATH_OUTPUT),
        path_color=expand("{path_output}/coverage/{{name_msa}}-color.jpg", path_output=PATH_OUTPUT)
    shell:
        "python coverage.py --path_blocks {input[0]} --path_msa {input[1]} --path_gray {output[0]} --path_color {output[1]}"

rule postprocessing_gfa:
    input:
        path_gfa=expand("{path_output}/gfa/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    output:
        path_post_gfa=expand("{path_output}/gfa-post/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    shell:
        "python src/postprocess_gfa.py --path_gfa {input} --path_save {output}"

rule unchop_gfa:
    input:
        path_post_gfa=expand("{path_output}/gfa-post/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    output:
        path_unchop_gfa=expand("{path_output}/gfa-unchop/{{name_msa}}.gfa", path_output=PATH_OUTPUT),
        path_labels=expand("{path_output}/gfa-unchop/{{name_msa}}.csv", path_output=PATH_OUTPUT)
    shell:
        """
        ../vg mod -u {input} > {output.path_unchop_gfa}
        python src/graph/bandage_labels_from_gfa.py --path_gfa {output.path_unchop_gfa} --path_save {output.path_labels}
        """