configfile: "params.yaml"

# --- EDA MSA
rule eda_msa: 
    input: 
        path_msas=config["PATH_MSAS"]
    output:
        "out/stats_msas.tsv"
    shell:
        "python3 src/eda_msas.py {input.path_msas}"