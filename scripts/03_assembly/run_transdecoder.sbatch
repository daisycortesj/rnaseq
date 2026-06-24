#!/bin/bash
# =============================================================================
# SCRIPT:  run_transdecoder.sbatch
#
# WHAT THIS DOES:
#   Trinity gives you DNA sequences (nucleotides: A, T, G, C).
#   BLAST and HMMER need PROTEIN sequences (amino acids).
#   This script runs TransDecoder on the HPC cluster: it finds the coding
#   region (ORF) in each transcript and writes a protein FASTA file.
#
# WHEN TO RUN (Step 7 in docs/README_trinity.md):
#   AFTER Step 5 (RSEM) and Step 6 (count matrix).
#
# HOW TO SUBMIT (MF nutmeg example):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/03_assembly/run_transdecoder.sbatch MF
#
# TO RE-RUN FROM SCRATCH:
#   sbatch scripts/03_assembly/run_transdecoder.sbatch MF --force
#
# OUTPUT (MF example):
#   01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
#
# NOTE ON FILE NAMING:
#   Every pipeline step in this repo uses ONE file named run_<tool>.sbatch
#   (e.g. run_cdhit.sbatch, run_rsem.sbatch). You only ever submit the
#   .sbatch file — there is no separate .sh helper file.
#
# PATHS come from scripts/config.sh (BASE_DIR, PROCESSED_DIR, etc.)
# =============================================================================


# ---- SLURM: tell the cluster how much time/memory this job needs ----
# TransDecoder can run many hours on a large assembly — we request 3 days.
#SBATCH --job-name=transdecoder
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit
#SBATCH -o transdecoder_%j.out
#SBATCH -e transdecoder_%j.err

# set -e  = stop if any command fails
# set -u  = stop if we use a variable that was never set
# set -o pipefail = treat a failed command in a pipe as a failure
set -euo pipefail


# =============================================================================
# SECTION 1: LOAD config.sh
# =============================================================================
# config.sh defines BASE_DIR, PROCESSED_DIR, CONDA_ENV, and get_sample_info()

CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../../config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not find scripts/config.sh."
    echo "Submit from the repo folder:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# SECTION 2: READ WHAT THE USER TYPED AFTER sbatch
# =============================================================================
# Example: sbatch run_transdecoder.sbatch MF cdhit --force
#   $1 = MF       (species code — required)
#   $2 = cdhit    (optional — which assembly to use; default is cdhit)
#   --force       (optional — re-run even if output exists)

SPECIES_CODE="${1:-}"       # e.g. MF, DC, DG
ASSEMBLY_TYPE="${2:-cdhit}" # cdhit (default) or pooled
FORCE=0                     # 0 = skip if output exists; 1 = re-run

# Look through ALL arguments for --force
for arg in "$@"; do
    if [[ "${arg}" == "--force" ]]; then
        FORCE=1
    fi
done

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: You must give a species code (e.g. MF for nutmeg)."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch <CODE> [cdhit|pooled] [--force]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF cdhit --force"
    exit 1
fi

# get_sample_info is defined in config.sh — sets SPECIES_NAME, GENOME_TYPE, etc.
get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# SECTION 3: CHECK THAT RSEM FINISHED (Step 5 must be done first)
# =============================================================================
# TransDecoder reads the CD-HIT FASTA, not RSEM files — but we check RSEM
# finished so you don't skip ahead in the pipeline by accident.

RSEM_MATRIX="${PROCESSED_DIR}/00_7_RSEM/RSEM.gene.counts.matrix"

echo ""
echo "Checking that RSEM finished (Step 5)..."
if [[ ! -f "${RSEM_MATRIX}" ]]; then
    echo "ERROR: RSEM count matrix not found:"
    echo "  ${RSEM_MATRIX}"
    echo ""
    echo "Run these first:"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch ${SPECIES_CODE_UPPER}"
    echo "  sbatch scripts/04_counting/run_count_matrix.sbatch ${SPECIES_CODE_UPPER} --type rsem"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch ${SPECIES_CODE_UPPER}"
    exit 1
fi
echo "  OK — RSEM matrix found: ${RSEM_MATRIX}"


# =============================================================================
# SECTION 4: PICK THE INPUT FASTA FOR THIS SPECIES
# =============================================================================
# Same path logic as run_rsem.sbatch so transcript IDs match DESeq2 counts.

case "${SPECIES_CODE_UPPER}" in
    MF)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_MF_trinity"
        ;;
    DC)
        TRINITY_DIR="${PROCESSED_DIR}/00_1_DC_trinity"
        ;;
    DG)
        TRINITY_DIR="${PROCESSED_DIR}/00_2_DG_trinity"
        ;;
    SK)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_Sukman_trinity"
        ;;
    DCDG)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_DC_DG_trinity"
        ;;
    *)
        echo "ERROR: Species '${SPECIES_CODE_UPPER}' is not set up in this script."
        echo "Copy a case block from run_rsem.sbatch and add it here."
        exit 1
        ;;
esac

case "${ASSEMBLY_TYPE}" in
    cdhit|CDHIT)
        TRINITY_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"
        ASSEMBLY_LABEL="CD-HIT (95% — same IDs as RSEM)"
        ;;
    pooled|POOL)
        TRINITY_FASTA="${TRINITY_DIR}/${SPECIES_CODE_UPPER}_trinity_pooled.Trinity.fasta"
        ASSEMBLY_LABEL="pooled Trinity (no CD-HIT — not recommended)"
        ;;
    --force)
        # User typed: sbatch ... MF --force  (skipped assembly type)
        TRINITY_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"
        ASSEMBLY_LABEL="CD-HIT (95% — same IDs as RSEM)"
        ;;
    *)
        echo "ERROR: Unknown assembly type '${ASSEMBLY_TYPE}'."
        echo "Use cdhit (default) or pooled."
        exit 1
        ;;
esac

# TransDecoder writes output files next to the input FASTA
WORK_DIR="$(dirname "${TRINITY_FASTA}")"
FASTA_BASENAME="$(basename "${TRINITY_FASTA}")"
PEP_OUTPUT="${TRINITY_FASTA}.transdecoder.pep"

JOB_START=$(date +%s)


# =============================================================================
# SECTION 5: ACTIVATE CONDA (TransDecoder lives in the rnaseq environment)
# =============================================================================
echo ""
echo "Activating conda environment: ${CONDA_ENV}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"


# =============================================================================
# SECTION 6: CHECK THAT TransDecoder IS INSTALLED
# =============================================================================
LONGORFS="TransDecoder.LongOrfs"
PREDICT="TransDecoder.Predict"

if ! command -v "${LONGORFS}" &> /dev/null; then
    echo ""
    echo "ERROR: ${LONGORFS} not found."
    echo "Install once on the login node:"
    echo "  conda activate rnaseq"
    echo "  conda install -y -c bioconda transdecoder"
    exit 1
fi

if ! command -v "${PREDICT}" &> /dev/null; then
    echo ""
    echo "ERROR: ${PREDICT} not found."
    echo "Install once on the login node:"
    echo "  conda activate rnaseq"
    echo "  conda install -y -c bioconda transdecoder"
    exit 1
fi


# =============================================================================
# SECTION 7: CHECK INPUT + SKIP IF ALREADY DONE
# =============================================================================
echo "========================================================"
echo "TransDecoder — turn nucleotide transcripts into proteins"
echo "Species:     ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Assembly:    ${ASSEMBLY_LABEL}"
echo "Input:       ${TRINITY_FASTA}"
echo "Output:      ${PEP_OUTPUT}"
echo "Force:       ${FORCE}  (0=skip if exists, 1=re-run)"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"
echo "Started:     $(date)"
echo "========================================================"

if [[ ! -f "${TRINITY_FASTA}" ]]; then
    echo "ERROR: Input FASTA not found:"
    echo "  ${TRINITY_FASTA}"
    echo "Run CD-HIT first: sbatch scripts/04_busco/run_cdhit.sbatch ${SPECIES_CODE_UPPER}"
    exit 1
fi

TRANSCRIPT_COUNT=$(grep -c "^>" "${TRINITY_FASTA}" || true)
echo "Transcripts in input: ${TRANSCRIPT_COUNT}"

# Skip if protein file already exists (TransDecoder is slow — don't redo by accident)
if [[ -f "${PEP_OUTPUT}" && "${FORCE}" -eq 0 ]]; then
    PEP_COUNT=$(grep -c "^>" "${PEP_OUTPUT}" || true)
    echo ""
    echo "SKIP — protein file already exists (${PEP_COUNT} proteins):"
    echo "  ${PEP_OUTPUT}"
    echo "Nothing was re-run. To start over, add --force."
    exit 0
fi

if [[ -f "${PEP_OUTPUT}" && "${FORCE}" -eq 1 ]]; then
    echo ""
    echo "NOTE: --force set. Deleting old TransDecoder files..."
    rm -f \
        "${TRINITY_FASTA}.transdecoder.pep" \
        "${TRINITY_FASTA}.transdecoder.cds" \
        "${TRINITY_FASTA}.transdecoder.gff3" \
        "${TRINITY_FASTA}.transdecoder.bed" \
        "${TRINITY_FASTA}.transdecoder.dir" \
        2>/dev/null || true
fi


# =============================================================================
# SECTION 8: RUN TransDecoder (two commands)
# =============================================================================
# TransDecoder must run in the folder where the FASTA lives
cd "${WORK_DIR}"

echo ""
echo "Step 1 of 2: TransDecoder.LongOrfs"
echo "  Scans every transcript for start/stop codons (slow step)."
echo "  Command: TransDecoder.LongOrfs -t ${FASTA_BASENAME}"
echo "  Started: $(date)"

STEP1_START=$(date +%s)
TransDecoder.LongOrfs -t "${FASTA_BASENAME}"
STEP1_END=$(date +%s)
STEP1_MINS=$(( (STEP1_END - STEP1_START) / 60 ))
echo "Step 1 finished in ${STEP1_MINS} minutes."

echo ""
echo "Step 2 of 2: TransDecoder.Predict"
echo "  Keeps one best protein per transcript."
echo "  Command: TransDecoder.Predict -t ${FASTA_BASENAME} --single_best_only"
echo "  Started: $(date)"

STEP2_START=$(date +%s)
TransDecoder.Predict -t "${FASTA_BASENAME}" --single_best_only
STEP2_END=$(date +%s)
STEP2_MINS=$(( (STEP2_END - STEP2_START) / 60 ))
echo "Step 2 finished in ${STEP2_MINS} minutes."


# =============================================================================
# SECTION 9: PRINT SUMMARY
# =============================================================================
if [[ ! -f "${PEP_OUTPUT}" ]]; then
    echo "ERROR: Expected output not created: ${PEP_OUTPUT}"
    exit 1
fi

PEP_COUNT=$(grep -c "^>" "${PEP_OUTPUT}" || true)
JOB_END=$(date +%s)
TOTAL_MINS=$(( (JOB_END - JOB_START) / 60 ))

echo ""
echo "========================================================"
echo "DONE — TransDecoder finished successfully"
echo "========================================================"
echo "Input transcripts:  ${TRANSCRIPT_COUNT}"
echo "Output proteins:    ${PEP_COUNT}"
echo "Main file:          ${PEP_OUTPUT}"
echo "Step 1 runtime:     ${STEP1_MINS} min"
echo "Step 2 runtime:     ${STEP2_MINS} min"
echo "Total runtime:      ${TOTAL_MINS} min"
echo "Finished:           $(date)"
echo "========================================================"

exit 0
