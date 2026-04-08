#!/bin/bash
# =====================================================================
# SCRIPT:  run_fastp.sbatch
# PURPOSE: Trim adapters and filter low-quality reads for ONE species
#
# WHAT IS fastp?
#   fastp takes raw FASTQ files and cleans them up by:
#   - Removing adapter sequences (leftover sequencing artifacts)
#   - Trimming low-quality bases from the ends of reads
#   - Removing reads that are too short after trimming
#   This gives you higher-quality data for alignment and counting.
#
# USAGE:
#   sbatch scripts/01_qc/run_fastp.sbatch <SPECIES_CODE>
#
# EXAMPLES:
#   sbatch scripts/01_qc/run_fastp.sbatch DC    # D. carota
#   sbatch scripts/01_qc/run_fastp.sbatch DG    # D. glaber
#   sbatch scripts/01_qc/run_fastp.sbatch MF    # M. fragrans
#
# INPUTS:
#   Raw .fq.gz files in: 00_rawdata/<SPECIES_DIR>/
#
# OUTPUTS (saved to 01_processed/<SPECIES_DIR>/):
#   - *_1_clean.fq.gz, *_2_clean.fq.gz  ← cleaned reads
#   - *_fastp.html                        ← per-sample QC report
#   - *_fastp.json                        ← machine-readable stats
# =====================================================================


# ===== SLURM directives =====
#SBATCH --job-name=fastp
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 01_processed/fastp_%j.out
#SBATCH -e 01_processed/fastp_%j.err

set -euo pipefail


# =====================================================================
# LOAD CONFIG AND PARSE SPECIES CODE
# =====================================================================

source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# Read species code from command line (e.g., MF, DC, DG)
SPECIES_CODE="${1:-}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:  sbatch run_fastp.sbatch <SPECIES_CODE>"
    echo ""
    echo "Valid codes: DC, DG, MF, SK, DCDG"
    echo "To add a new code, edit SECTION 4 in scripts/config.sh"
    exit 1
fi

# Look up species info (sets SPECIES_DIR, SPECIES_NAME, GENOME_TYPE)
get_sample_info "${SPECIES_CODE}"


# =====================================================================
# BUILD PATHS AND SET PARAMETERS
# =====================================================================

# Where the raw FASTQ files are
SPECIES_FASTQ_DIR="${RAWDATA_DIR}/${SPECIES_DIR}"

# Where to save cleaned reads — same folder pattern as rawdata
# e.g., 01_processed/00_3_MF/
FASTP_OUTPUT="${PROCESSED_DIR}/${SPECIES_DIR}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- Trimming parameters (change these if needed) ---
MIN_QUALITY=20    # Drop bases below this quality (Q20 = 99% accuracy)
MIN_LENGTH=50     # Drop reads shorter than this after trimming


# =====================================================================
# SAFETY CHECKS
# =====================================================================

if [[ ! -d "${SPECIES_FASTQ_DIR}" ]]; then
    echo "ERROR: Input directory not found: ${SPECIES_FASTQ_DIR}"
    echo "Run: ls ${RAWDATA_DIR}/"
    exit 1
fi

FASTQ_COUNT=$(find "${SPECIES_FASTQ_DIR}" -name "*_1.fq.gz" -o -name "*_R1.fq.gz" | wc -l)
if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
    echo "ERROR: No paired FASTQ files found in ${SPECIES_FASTQ_DIR}"
    exit 1
fi


# =====================================================================
# PRINT JOB INFO
# =====================================================================

echo "============================================================"
echo "fastp — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "============================================================"
echo "Host: $(hostname)  CPUs: ${THREADS}  Start: $(date)"
echo "Input:  ${SPECIES_FASTQ_DIR}"
echo "Output: ${FASTP_OUTPUT}"
echo "Samples found: ${FASTQ_COUNT}"
echo "Quality cutoff: Q${MIN_QUALITY}  Min length: ${MIN_LENGTH}bp"
echo "============================================================"


# =====================================================================
# ACTIVATE CONDA
# =====================================================================

mkdir -p "${FASTP_OUTPUT}"

if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    source ~/.bashrc
    conda activate rnaseq
fi

# Make sure fastp is available
if ! command -v fastp &> /dev/null; then
    echo "ERROR: fastp not found. Install with: conda install -c bioconda fastp"
    exit 1
fi

echo "fastp version: $(fastp --version 2>&1)"
echo "------------------------------------------------------------"


# =====================================================================
# RUN fastp ON EACH SAMPLE
# =====================================================================
# fastp processes one sample (pair of R1+R2 files) at a time.
# We loop through all R1 files and find the matching R2 file.

SAMPLE_COUNT=0
CLEANED_COUNT=0
FAILED_COUNT=0

for R1 in $(find "${SPECIES_FASTQ_DIR}" -name "*_1.fq.gz" -o -name "*_R1.fq.gz" | sort); do

    # --- Find the matching R2 file ---
    # If R1 is MFF1_1.fq.gz, R2 is MFF1_2.fq.gz
    R2="${R1/_1.fq.gz/_2.fq.gz}"
    R2="${R2/_R1.fq.gz/_R2.fq.gz}"

    # Extract sample name from filename (e.g., MFF1 from MFF1_1.fq.gz)
    SAMPLE=$(basename "${R1}" | sed 's/_1\.fq\.gz$//' | sed 's/_R1\.fq\.gz$//')

    # Check that the R2 file exists
    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: R2 not found for ${SAMPLE}. Expected: ${R2}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi

    # --- Set output filenames ---
    OUT_R1="${FASTP_OUTPUT}/${SAMPLE}_1_clean.fq.gz"
    OUT_R2="${FASTP_OUTPUT}/${SAMPLE}_2_clean.fq.gz"
    JSON_REPORT="${FASTP_OUTPUT}/${SAMPLE}_fastp.json"
    HTML_REPORT="${FASTP_OUTPUT}/${SAMPLE}_fastp.html"

    # --- Skip if already processed ---
    if [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]] && [[ -f "${HTML_REPORT}" ]]; then
        echo "Skipping ${SAMPLE} (already processed)"
        SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
        continue
    fi

    echo "------------------------------------------------------------"
    echo "Processing: ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"

    # --- Run fastp ---
    # -i / -I            = input R1 / R2
    # -o / -O            = output cleaned R1 / R2
    # --detect_adapter_for_pe = auto-detect adapters for paired-end reads
    # --qualified_quality_phred = minimum base quality to keep
    # --length_required  = minimum read length after trimming
    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT_R1}" \
        -O "${OUT_R2}" \
        --json "${JSON_REPORT}" \
        --html "${HTML_REPORT}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred ${MIN_QUALITY} \
        --length_required ${MIN_LENGTH} \
        --thread ${THREADS}

    # Check if it worked
    if [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]]; then
        echo "Done: ${SAMPLE}"
        CLEANED_COUNT=$((CLEANED_COUNT + 1))
    else
        echo "FAILED: ${SAMPLE}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done


# =====================================================================
# SUMMARY
# =====================================================================

echo "============================================================"
echo "fastp completed — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "End time: $(date)"
echo ""
echo "  Samples processed: ${SAMPLE_COUNT}"
echo "  Successful:        ${CLEANED_COUNT}"
echo "  Failed:            ${FAILED_COUNT}"
echo ""
echo "RESULTS in: ${FASTP_OUTPUT}/"
echo "  Cleaned reads: *_clean.fq.gz"
echo "  Reports:       *_fastp.html"
echo ""
echo "NEXT STEPS:"
echo "1. Re-run FastQC on cleaned data to compare before/after:"
echo "   sbatch scripts/01_qc/run_fastqc.sbatch ${SPECIES_CODE}"
echo "2. Proceed to alignment (STAR) with the cleaned reads"
echo "============================================================"

# --- Copy log files into the species folder ---
cp "01_processed/fastp_${SLURM_JOB_ID}.out" "${FASTP_OUTPUT}/" 2>/dev/null || true
cp "01_processed/fastp_${SLURM_JOB_ID}.err" "${FASTP_OUTPUT}/" 2>/dev/null || true
