#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Requests cluster resources to activate RNA seq enviroment and runs rnaseq pipeline. It forwards all agruments to rna_pipeline.cli to execute STAR or Trinity
# ---------------------------------------------------------------------

# ===== SLURM directives (ignored when run with "bash ...", used by sbatch) =====
#SBATCH --job-name=index_star
#SBATCH --account=tholl_lab_1      # <-- change to your ARC allocation
#SBATCH --partition=normal_q
## (optional) QoS (≤24h): uncomment next line if desired
## #SBATCH --qos=tc_normal_short
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
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}  MEM: ${SLURM_MEM_PER_NODE:-64G}  TIME: ${SLURM_JOB_TIMELIMIT:-unset}"


# 1) Activate your conda env on the compute node
source ~/.bashrc
conda activate rnaseq   # <-- change if your env name is different


# 2) Go to repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"


# 3) Check that STAR and Trinity are avaliable
echo "Checking tool versions:"
STAR --version || echo "STAR not found."
Trinity --version || echo "Trinity not found."


# 4) Forward ALL CLI args to the Python pipeline
echo "▶ Running: python -m rna_pipeline.cli $@"
srun python -m rna_pipeline.cli "$@"
echo "------------------------------------------------------------"
echo "[INFO] Pipeline completed successfully."
echo "------------------------------------------------------------"

