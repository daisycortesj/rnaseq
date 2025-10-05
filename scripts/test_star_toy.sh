
#!/usr/bin/env bash
set -euo pipefail

python -m rna_pipeline.cli   --fasta data/reference/toy.fa   --gtf   data/reference/toy.gtf   --outdir results/star_index_toy   --threads 1   --readlen 50
