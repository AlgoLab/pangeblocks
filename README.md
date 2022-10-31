```bash
python3 -m venv pbenv
source pbenv/bin/activate
pip install -r requirements.txt
```

Fix parameters in `params.yaml`
```bash
snakemake -s eda.smk -c1         # compute stats for each MSA
snakemake -s max_blocks.smk -c10 # compute max blocks for MSAs
```