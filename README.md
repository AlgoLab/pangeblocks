```bash
python3 -m venv pbenv
source pbenv/bin/activate
pip install -r requirements.txt
```

Fix parameters in `params.yaml`
```bash
snakemake -s eda.smk -c16         # compute stats for each MSA
snakemake -s max_blocks.smk -c16 # compute max blocks for MSAs
```