configfile: "params.yaml"
from pathlib import Path
PATH_MSAS=config["PATH_MSAS"]
PATH_OUTPUT=config["PATH_OUTPUT"]

NAMES = [path.stem for path in Path(PATH_MSAS).rglob("*.fa")]
NAMES.extend([Path(path.stem).with_suffix("") for path in Path(PATH_MSAS).rglob("*.fa.gz")])
NAMES = [str(name_msa) for name_msa in NAMES]
print(NAMES)

rule all:
    input:
        expand("{path_output}/output-pggb/{name_msa}.gfa", name_msa=NAMES, path_output=PATH_OUTPUT)
        # expand("{path_msas}/{name_msa}.fa.gz.fai", name_msa=NAMES, path_msas=PATH_MSAS),
        # expand("{path_msas}/{name_msa}.fa.gz.gzi", name_msa=NAMES, path_msas=PATH_MSAS)

rule install_pggb:
    input: 
        remove_cache = "/opt/mambaforge/pkgs/cache/"
    output: 
        path_env= f"{str(Path().resolve())}/pggb-env"
    shell:
        """
        rm -rf input.remove_cache
        mamba create --prefix {output.path_env} -c conda-forge -c bioconda pggb
        """

rule compress_fasta:
    input:
        expand("{path_msas}/{{name_msa}}.fa", path_msas=PATH_MSAS)
    output:
        expand("{path_msas}/{{name_msa}}.fa.gz", path_msas=PATH_MSAS)
    shell:
        "bgzip -@ 16 {input}"

rule index_fasta: 
    input: 
        expand("{path_msas}/{{name_msa}}.fa.gz", path_msas=PATH_MSAS)
    output:
        expand("{path_msas}/{{name_msa}}.fa.gz.fai", path_msas=PATH_MSAS),
        expand("{path_msas}/{{name_msa}}.fa.gz.gzi", path_msas=PATH_MSAS)
    shell:
        "samtools faidx {input}"

rule run_pggb:
    input: 
        expand("{path_msas}/{{name_msa}}.fa.gz", path_msas=PATH_MSAS)
    output:
        gfa_output = expand("{path_output}/output-pggb/{{name_msa}}.gfa", path_output=PATH_OUTPUT)
    shell:
        """
        pggb -i {input} -p 70 -s 500 -n 10 -t 16 -o {PATH_OUTPUT}/output-pggb/
        mv {PATH_OUTPUT}/output-pggb/{wildcards.name_msa}*.gfa {PATH_OUTPUT}/output-pggb/{wildcards.name_msa}.gfa 
        """