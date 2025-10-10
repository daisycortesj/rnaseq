#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Requests cluster resources to activate RNA seq environment and runs rnaseq pipeline
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=index_star
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/04_reference
#SBATCH -o bldstridx_output.txt
#SBATCH -e bldstridx_error.txt

set -euo pipefail

echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}  MEM: ${SLURM_MEM_PER_NODE:-128G}  TIME: ${SLURM_JOB_TIMELIMIT:-unset}"

# Check and fix .bashrc if needed
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    echo "Loading .bashrc successfully"
    source ~/.bashrc
    conda activate rnaseq
fi

# Verify conda activation worked
conda info --envs | grep rnaseq || echo "ERROR: conda environment not activated"

# Go to pipeline directory (from reference directory)
cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code

# Check tools
echo "Checking tool versions:"
STAR --version || echo "STAR not found."
Trinity --version || echo "Trinity not found."

# Run Python pipeline with all arguments
echo "▶ Running: python -m rna_pipeline.cli $@"
python -m rna_pipeline.cli "$@"
echo "------------------------------------------------------------"
echo "[INFO] Pipeline completed successfully."
echo "------------------------------------------------------------"
