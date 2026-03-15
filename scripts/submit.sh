#!/bin/bash
# =====================================================================
# submit.sh — Submit sbatch jobs with settings from config.sh
# =====================================================================
#
# This wrapper automatically passes --account and --mail-user
# from config.sh so those #SBATCH lines don't need to be hardcoded.
# Each sbatch script controls its own --chdir via #SBATCH directives.
#
# Usage:
#   ./scripts/submit.sh scripts/01_qc/run_fastp.sbatch
#   ./scripts/submit.sh scripts/05_pydeseq2/run_step1_analysis.sbatch DC
#   ./scripts/submit.sh scripts/01_qc/run_fastp.sbatch /other/rawdata
#
# =====================================================================

# Load config
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# First argument = the .sbatch script, rest = passed to the script
SBATCH_SCRIPT="$1"
shift

if [ -z "${SBATCH_SCRIPT}" ]; then
    echo "Usage: ./scripts/submit.sh <script.sbatch> [arguments...]"
    echo ""
    echo "Examples:"
    echo "  ./scripts/submit.sh scripts/01_qc/run_fastqc.sbatch"
    echo "  ./scripts/submit.sh scripts/01_qc/run_fastp.sbatch"
    echo "  ./scripts/submit.sh scripts/05_pydeseq2/run_step1_analysis.sbatch DC"
    exit 1
fi

if [ ! -f "${SBATCH_SCRIPT}" ]; then
    echo "ERROR: Script not found: ${SBATCH_SCRIPT}"
    exit 1
fi

echo "Submitting: ${SBATCH_SCRIPT}"
echo "  Account:  ${SLURM_ACCOUNT}"
echo "  Email:    ${MAIL_USER}"
echo "  Base dir: ${BASE_DIR}"
echo ""

sbatch \
    "${SBATCH_SCRIPT}" \
    "$@"
