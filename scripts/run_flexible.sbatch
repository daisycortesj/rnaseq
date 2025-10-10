#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Ultra-flexible SLURM script for RNA-seq pipeline
# USAGE: sbatch scripts/run_flexible.sbatch [SLURM_OPTIONS] -- [PIPELINE_OPTIONS]
# 
# EXAMPLES:
#   # Basic STAR indexing
#   sbatch scripts/run_flexible.sbatch \
#     --fasta /path/to/genome.fa \
#     --gtf /path/to/genes.gtf \
#     --outdir star_index
#
#   # STAR with custom SLURM resources
#   sbatch scripts/run_flexible.sbatch \
#     --job-name=carrot_star \
#     --cpus-per-task=32 \
#     --mem=256G \
#     --time=24:00:00 \
#     --chdir=/path/to/working/dir \
#     -- \
#     --fasta /path/to/carrot_genome.fa \
#     --gtf /path/to/carrot_genes.gtf \
#     --outdir carrot_star_index \
#     --threads 32
#
#   # Trinity de novo assembly
#   sbatch scripts/run_flexible.sbatch \
#     --job-name=trinity_assembly \
#     --mem=512G \
#     --time=48:00:00 \
#     -- \
#     --reads-left /path/to/reads_R1.fq.gz \
#     --reads-right /path/to/reads_R2.fq.gz \
#     --outdir trinity_results \
#     --threads 16 \
#     --mem-gb 128
# ---------------------------------------------------------------------

# ===== Default SLURM directives =====
#SBATCH --job-name=rna_pipeline
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH -o rna_pipeline_%j.out
#SBATCH -e rna_pipeline_%j.err

set -euo pipefail

echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}  MEM: ${SLURM_MEM_PER_NODE:-128G}  TIME: ${SLURM_JOB_TIMELIMIT:-12:00:00}"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"

# 1) Activate conda environment
echo "Activating conda environment..."
# Try multiple conda initialization methods
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
elif [ -f ~/anaconda3/etc/profile.d/conda.sh ]; then
    source ~/anaconda3/etc/profile.d/conda.sh
elif [ -f /opt/miniconda3/etc/profile.d/conda.sh ]; then
    source /opt/miniconda3/etc/profile.d/conda.sh
else
    echo "Warning: Could not find conda initialization script"
    echo "Trying to source ~/.bashrc..."
    set +u  # Temporarily disable unbound variable checking
    source ~/.bashrc || echo "Warning: Could not source ~/.bashrc"
    set -u  # Re-enable unbound variable checking
fi

conda activate rnaseq

# 2) Set up Python path
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH}"
echo "Python path set to: ${PYTHONPATH}"

# 3) Check that required tools are available
echo "Checking tool versions:"
STAR --version || echo "STAR not found."
Trinity --version || echo "Trinity not found."

# 4) Run the RNA-seq pipeline with all provided arguments
echo "============================================================"
echo "▶ Arguments received: $*"
echo "▶ Running: python -m rna_pipeline.cli $*"
echo "============================================================"

srun python -m rna_pipeline.cli $*

echo "============================================================"
echo "[INFO] Pipeline completed successfully."
echo "Job finished at: $(date)"
echo "============================================================"
