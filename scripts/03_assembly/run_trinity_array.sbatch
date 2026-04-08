#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity de novo assembly as a SLURM job array
#          Each sample runs as an independent job with its own time limit
#
# USAGE (two steps):
#   1. Generate the sample list for one species (auto-detects clean vs raw):
#        bash scripts/03_assembly/run_trinity_array.sbatch --generate MF
#      This creates 01_processed/sample_list.txt
#
#   2. Submit the array (replace N with the sample count printed above):
#        sbatch --array=1-N scripts/03_assembly/run_trinity_array.sbatch
#
# RESUME:  If jobs time out, simply resubmit - incomplete assemblies
#          will automatically resume from checkpoint.
#
# ADVANTAGES over run_trinity_all.sbatch:
#   - Each sample has its own 72-hour time limit
#   - Samples can run in parallel (faster overall)
#   - If one sample fails/times out, others are unaffected
#
# PATHS: All directory paths come from scripts/config.sh
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_arr
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_array_%A_%a.out
#SBATCH -e trinity_array_%A_%a.err

set -euo pipefail

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# ===== Sample list generator mode =====
# Run with: bash run_trinity_array.sbatch --generate MF
# This creates the sample_list.txt with the correct files (clean or raw)
if [[ "${1:-}" == "--generate" ]]; then
    SPECIES_CODE="${2:-}"
    if [[ -z "${SPECIES_CODE}" ]]; then
        echo "Usage: bash run_trinity_array.sbatch --generate <SPECIES_CODE>"
        echo "Example: bash run_trinity_array.sbatch --generate MF"
        exit 1
    fi

    get_sample_info "${SPECIES_CODE}"
    SAMPLE_LIST="${PROCESSED_DIR}/sample_list.txt"

    # Auto-detect: clean reads if available, otherwise raw
    CLEAN_DIR="${PROCESSED_DIR}/${SPECIES_DIR}"
    CLEAN_COUNT=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" 2>/dev/null | wc -l)

    if [[ "${CLEAN_COUNT}" -gt 0 ]]; then
        find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" | sort > "${SAMPLE_LIST}"
        echo "Generated sample list using FASTP-CLEANED reads"
    else
        find "${RAWDATA_DIR}/${SPECIES_DIR}" -maxdepth 1 -name "*_1.fq.gz" | sort > "${SAMPLE_LIST}"
        echo "Generated sample list using RAW reads"
    fi

    SAMPLE_COUNT=$(wc -l < "${SAMPLE_LIST}")
    echo "Samples: ${SAMPLE_COUNT}"
    echo "File: ${SAMPLE_LIST}"
    cat "${SAMPLE_LIST}"
    echo ""
    echo "Now submit with:"
    echo "  sbatch --array=1-${SAMPLE_COUNT} scripts/03_assembly/run_trinity_array.sbatch"
    exit 0
fi

# ===== Configuration =====
SAMPLE_LIST="${PROCESSED_DIR}/sample_list.txt"
TRINITY_OUTDIR="${PROCESSED_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
MEM_GB="128"

echo "============================================================"
echo "Trinity De Novo Assembly - Job Array"
echo "============================================================"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Validate sample list =====
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    echo ""
    echo "Generate it first with:"
    echo "  bash scripts/03_assembly/run_trinity_array.sbatch --generate MF"
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Total samples in list: ${TOTAL_SAMPLES}"

if [[ ${SLURM_ARRAY_TASK_ID:-1} -gt ${TOTAL_SAMPLES} ]]; then
    echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} exceeds sample count ${TOTAL_SAMPLES}"
    exit 1
fi

# ===== Get sample for this task =====
LEFT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" "${SAMPLE_LIST}")
sample_dir=$(dirname "${LEFT_FILE}")
filename=$(basename "${LEFT_FILE}")

# Handle both clean and raw filenames:
#   MFF1_1_clean.fq.gz → base=MFF1, suffix=_clean
#   MFF1_1.fq.gz       → base=MFF1, suffix=
if [[ "${filename}" == *"_clean.fq.gz" ]]; then
    base="${filename%_1_clean.fq.gz}"
    FILE_SUFFIX="_clean"
else
    base="${filename%_1.fq.gz}"
    FILE_SUFFIX=""
fi
RIGHT_FILE="${sample_dir}/${base}_2${FILE_SUFFIX}.fq.gz"

echo ""
echo "Processing sample ${SLURM_ARRAY_TASK_ID:-1} of ${TOTAL_SAMPLES}: ${base}"
echo "  Left:  ${LEFT_FILE}"
echo "  Right: ${RIGHT_FILE}"

# ===== Validate files =====
if [[ ! -r "${LEFT_FILE}" ]]; then
    echo "ERROR: Cannot read left file: ${LEFT_FILE}"
    exit 1
fi

if [[ ! -r "${RIGHT_FILE}" ]]; then
    echo "ERROR: Cannot read right file: ${RIGHT_FILE}"
    exit 1
fi

echo "  ✓ Both read files found and readable"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    source ~/.bashrc
    conda activate rnaseq
fi

# Make sure Trinity is available
if ! command -v Trinity &> /dev/null; then
    echo "ERROR: Trinity not found. Install with: conda install -c bioconda trinity"
    exit 1
fi

echo "Trinity version: $(Trinity --version 2>&1 | head -1)"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Set output directory =====
OUTPUT_DIR="${TRINITY_OUTDIR}/${base}_trinity_assembly"
echo "  Output: ${OUTPUT_DIR}"

# ===== Check for existing/incomplete assembly =====
RESUME_FLAG=""
if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
        echo "  ✓ Assembly already complete (Trinity.fasta exists)"
        echo "  Exiting successfully - no work needed"
        exit 0
    else
        echo "  ⚠ Incomplete assembly detected - will RESUME from checkpoint"
        RESUME_FLAG="--resume"
    fi
fi

# ===== Run Trinity =====
echo ""
echo "Running Trinity assembly..."
CMD="python -m rna_pipeline.cli \
  --mode trinity \
  --reads-left \"${LEFT_FILE}\" \
  --reads-right \"${RIGHT_FILE}\" \
  --outdir \"${OUTPUT_DIR}\" \
  --threads ${THREADS} \
  --mem-gb ${MEM_GB} \
  ${RESUME_FLAG}"

echo "Command: ${CMD}"
echo "------------------------------------------------------------"

if eval ${CMD} 2>&1; then
    echo ""
    echo "✓ ${base} completed successfully"
    
    # Check for Trinity output
    if [[ -f "${OUTPUT_DIR}/Trinity.fasta" ]]; then
        CONTIG_COUNT=$(grep -c ">" "${OUTPUT_DIR}/Trinity.fasta" || echo "0")
        echo "✓ Trinity.fasta: ${CONTIG_COUNT} contigs"
    else
        echo "⚠ Warning: Trinity.fasta not found"
    fi
else
    EXIT_CODE=$?
    echo ""
    echo "✗ ${base} failed (exit code: ${EXIT_CODE})"
    exit ${EXIT_CODE}
fi

echo ""
echo "============================================================"
echo "End time: $(date)"
echo "============================================================"

exit 0
