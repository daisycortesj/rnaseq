
#!/bin/bash
#SBATCH --job-name=star-index
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

# Activate your conda env (adjust path/env name)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Optional: sanity check
bash scripts/inspect_env.sh || true

# Example: reference present → STAR index
python -m rna_pipeline.cli   --fasta data/reference/genome.fa   --gtf   data/reference/genes.gtf   --outdir results/star_index   --threads ${SLURM_CPUS_PER_TASK}   --readlen 150

# Example Trinity command (comment STAR above and uncomment below if no reference):
# python -m rna_pipeline.cli #   --reads-left  data/reads/sample_R1.fq.gz #   --reads-right data/reads/sample_R2.fq.gz #   --outdir results/trinity #   --threads ${SLURM_CPUS_PER_TASK} #   --mem-gb 64
