#!/bin/bash
# =============================================================================
# SCRIPT:  run_hmmer_genefinder.sbatch
#
# WHAT THIS DOES:
#   SLURM wrapper for trinity_hmmer_cyp_omt.sh (CYP/OMT pipeline Script 02).
#   Runs targeted HMMER hmmsearch on TransDecoder proteins from a Trinity
#   de novo assembly to find CYP450 and OMT candidate transcripts.
#
# INPUTS REQUIRED (must exist before you submit):
#   1. TransDecoder protein FASTA (.pep) — from run_transdecoder.sbatch
#      MF example:
#        01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
#
#   2. HMMER in rnaseq conda env (hmmsearch, hmmpress)
#        conda activate rnaseq && hmmsearch -h
#
#   3. Pfam HMM profiles — downloaded automatically on first run to:
#        07_NRdatabase/hmmerdb/pfam_profiles/PF00067.hmm  (CYP450)
#        07_NRdatabase/hmmerdb/pfam_profiles/PF00891.hmm  (OMT)
#        07_NRdatabase/hmmerdb/pfam_profiles/PF08100.hmm  (OMT)
#        07_NRdatabase/hmmerdb/pfam_profiles/PF01596.hmm  (CCoA-OMT)
#
# UPSTREAM PIPELINE (run these first for MF nutmeg):
#   trinity_pooled.sbatch → run_cdhit.sbatch → run_rsem.sbatch MF
#   → run_transdecoder.sbatch MF
#
# HOW TO SUBMIT:
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF
#
# TO RE-RUN:
#   sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force
#
# OUTPUT:
#   06_analysis/hmmer_genefinder_MF/
#     cyp450_hmmer_ids.txt
#     omt_hmmer_ids.txt
#     *.domtbl
# =============================================================================

#SBATCH --job-name=hmmer_genefinder
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o hmmer_genefinder_%j.out
#SBATCH -e hmmer_genefinder_%j.err

set -euo pipefail


# =============================================================================
# SECTION 1: FIND config.sh
# =============================================================================

CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not find scripts/config.sh"
    exit 1
fi

source "${CONFIG_FILE}"


# =============================================================================
# SECTION 2: ACTIVATE CONDA (HMMER lives in rnaseq env)
# =============================================================================

echo "Activating conda environment: ${CONDA_ENV}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"


# =============================================================================
# SECTION 3: RUN THE HMMER SCRIPT
# =============================================================================

SPECIES_CODE="${1:-}"
shift || true

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: Species code required."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force"
    echo ""
    echo "Required input (MF example):"
    echo "  ${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep"
    echo ""
    echo "Run TransDecoder first if missing:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF"
    exit 1
fi

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}/scripts/07_domains"
HMMER_SCRIPT="${SCRIPT_DIR}/trinity_hmmer_cyp_omt.sh"

if [[ ! -f "${HMMER_SCRIPT}" ]]; then
    echo "ERROR: Script not found: ${HMMER_SCRIPT}"
    exit 1
fi

chmod +x "${HMMER_SCRIPT}"

echo ""
echo "Running: bash ${HMMER_SCRIPT} ${SPECIES_CODE} $*"
echo ""

bash "${HMMER_SCRIPT}" "${SPECIES_CODE}" "$@"
