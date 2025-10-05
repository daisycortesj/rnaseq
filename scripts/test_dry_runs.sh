
#!/usr/bin/env bash
set -euo pipefail

echo "[1/3] Env check"
bash scripts/inspect_env.sh || true

echo "[2/3] DRY RUN — STAR (toy reference)"
python -m rna_pipeline.cli   --fasta data/reference/toy.fa   --gtf   data/reference/toy.gtf   --outdir results/star_index_toy   --threads 1   --readlen 50   --dry

echo "[3/3] DRY RUN — Trinity (toy reads)"
python -m rna_pipeline.cli   --reads-left  data/reads/toy_R1.fq.gz   --reads-right data/reads/toy_R2.fq.gz   --outdir results/trinity_toy   --threads 1   --mem-gb 4   --dry

echo "✓ Dry runs completed. To execute the tiny STAR index for real, run:"
echo "python -m rna_pipeline.cli --fasta data/reference/toy.fa --gtf data/reference/toy.gtf --outdir results/star_index_toy --threads 1 --readlen 50"
