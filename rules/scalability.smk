configfile: "rules/params/scalability.yaml"
from pathlib import Path
import pandas as pd
from os.path import join as pjoin

DECOMPOSITION = "complete" if config["DECOMPOSITION"]["STANDARD"] else "row-maximal"
print(DECOMPOSITION)

PATH_OUTPUT = config["PATH_OUTPUT"]
PATH_MSAS   = config["PATH_MSAS"]
LOG_LEVEL = config["LOG_LEVEL"]

OBJ_FUNCTIONS=config["OPTIMIZATION"]["OBJECTIVE_FUNCTION"]
print(OBJ_FUNCTIONS)

# 'weighted' and 'depth' loss functions
PENALIZATION=config["OPTIMIZATION"]["PENALIZATION"] 
MIN_LEN=config["OPTIMIZATION"]["MIN_LEN"]
MIN_COVERAGE=config["OPTIMIZATION"]["MIN_COVERAGE"]

# path msas
MSAS = list(Path(PATH_MSAS).glob("*.fa"))
NAMES = [path.stem for path in MSAS] # if "100" not in path.stem]
# NAMES = [path.stem for path in MSAS if "83" in path.stem]

print(NAMES)
EXT_MSA = MSAS[0].suffix

# subMSAS (start,end) columns
SUBMSAS = config["SUBMSAS"]
for submsa in SUBMSAS:
    print(submsa)

Path(PATH_OUTPUT).mkdir(parents=True, exist_ok=True)
config["NAMES"] = NAMES
with open(Path(PATH_OUTPUT).joinpath("config.json"), "w") as fp:
    json.dump(config, fp, indent=1)


def get_mem_mb(wildcards, attempt):
    return attempt*10000


def get_files(wildcards):
    "Return a list of files to be generated based on parameters provided in the config file"
    files = []
    if "nodes" in OBJ_FUNCTIONS:
        files.extend(
            pjoin(PATH_OUTPUT, "ilp", f"{name_msa}", f"{name_msa}-nodes-penalization0-min_len0-min_coverage0-start_{submsa[0]}-end_{submsa[1]}-rule-ilp.log")
            for name_msa in NAMES for submsa in SUBMSAS
        )
            
    if "strings" in OBJ_FUNCTIONS:
        files.extend(
            pjoin(PATH_OUTPUT, "ilp", f"{name_msa}", f"{name_msa}-strings-penalization0-min_len0-min_coverage0-start_{submsa[0]}-end_{submsa[1]}-rule-ilp.log")
            for name_msa in NAMES
        )

    if "weighted" in OBJ_FUNCTIONS:
        files.extend(
            pjoin(PATH_OUTPUT, "ilp", f"{name_msa}", f"{name_msa}-weighted-penalization{penalization}-min_len{min_len}-min_coverage0-start_{submsa[0]}-end_{submsa[1]}-rule-ilp.log")
            for penalization in PENALIZATION for min_len in MIN_LEN for name_msa in NAMES
        )

    if "depth" in OBJ_FUNCTIONS:
        files.extend(
            pjoin(PATH_OUTPUT, "ilp", f"{name_msa}", f"{name_msa}-depth-penalization{penalization}-min_len0-min_coverage{min_coverage}-start_{submsa[0]}-end_{submsa[1]}-rule-ilp.log")
            for penalization in PENALIZATION for min_coverage in MIN_COVERAGE for name_msa in NAMES
        )
    
    return files

rule all:
    input:
        get_files,
        "Wild-pBWT/bin/wild-pbwt"

rule install_wild_pbwt:
    params: 
        github = "https://github.com/illoxian/Wild-pBWT.git"
    output:
        "Wild-pBWT/bin/wild-pbwt"
    log:
        pjoin(PATH_OUTPUT, "logs", "rule-install_wild_pbwt.err.log"),
    conda:
        "../envs/wild-pbwt.yml"
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

rule ilp:
    input:
        path_msa=pjoin(PATH_MSAS, "{name_msa}" + EXT_MSA),
        bin_wildpbwt="Wild-pBWT/bin/wild-pbwt"
    output:
        auxfile=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-start_{start}-end_{end}-rule-ilp.log")
    params:
        dir_subsols=pjoin(PATH_OUTPUT, "ilp", "{name_msa}", "{obj_func}","penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}"),
        log_level=config["LOG_LEVEL"],
        time_limit=config["OPTIMIZATION"]["TIME_LIMIT"],
        threads_ilp=config["THREADS"]["ILP"],
        workers=config["THREADS"]["SUBMSAS"],
        use_wildpbwt=config["USE_WILDPBWT"],
        standard_decomposition=config["DECOMPOSITION"]["STANDARD"],
        alpha_consistent=config["DECOMPOSITION"]["ALPHA_CONSISTENT"]
    threads:
        config["THREADS"]["ILP"]
    log:
        stderr=pjoin(PATH_OUTPUT, "logs", "{name_msa}-{obj_func}-penalization{penalization}-min_len{min_len}-min_coverage{min_coverage}-start_{start}-end_{end}-rule-ilp.log"),
    conda: 
        "../envs/pangeblocks.yml"
    resources: 
        mem_mb=get_mem_mb
    shell:
        """
        /usr/bin/time --verbose src/exact_cover.py --path-msa {input.path_msa} --obj-function {wildcards.obj_func} \
        --prefix-output {params.dir_subsols}/{wildcards.name_msa} \
        --penalization {wildcards.penalization} --min-len {wildcards.min_len} --min-coverage {wildcards.min_coverage} \
        --start-column {wildcards.start} --end-column {wildcards.end} \
        --time-limit {params.time_limit} --solve-ilp True \
        --use-wildpbwt {params.use_wildpbwt} --bin-wildpbwt {input.bin_wildpbwt} \
        --standard-decomposition {params.standard_decomposition} --threads-ilp {params.threads_ilp} \
        --solve-ilp True \
        --workers {params.workers} --alpha-consistent {params.alpha_consistent} > {output.auxfile} 2> {log.stderr}
        """
