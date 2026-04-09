#!/bin/bash
# =====================================================================
# SCRIPT:  run_ribodetector.sbatch
# PURPOSE: Remove residual ribosomal RNA (rRNA) from fastp-cleaned reads
#
# WHAT IS RiboDetector?
#   RiboDetector is a deep-learning tool that identifies and removes
#   ribosomal RNA sequences from your FASTQ files. Even after standard
#   library prep, some rRNA can sneak through. Shows up as high
#   rRNA mapping rates in STAR. Running RiboDetector gives you cleaner
#   reads so your downstream analysis (alignment, counting, DESeq2)
#   focuses on real mRNA transcripts, not rRNA contamination.
# USAGE:
#   sbatch scripts/01_qc/run_ribodetector.sbatch <SPECIES_CODE>
#
# EXAMPLES:
#   sbatch scripts/01_qc/run_ribodetector.sbatch MF    # M. fragrans
#   sbatch scripts/01_qc/run_ribodetector.sbatch DC    # D. carota
#   sbatch scripts/01_qc/run_ribodetector.sbatch DG    # D. glaber
#
# INPUTS:
#   Cleaned reads from fastp in: 01_processed/<SPECIES_DIR>/
#     Files matching: *_1_clean.fq.gz and *_2_clean.fq.gz
#
# OUTPUTS (saved to 01_processed/<SPECIES_DIR>/ribofree/):
#   - *_1_ribofree.fq.gz, *_2_ribofree.fq.gz  ← reads with rRNA removed
#   - *_1_rrna.fq.gz, *_2_rrna.fq.gz          ← the rRNA reads (for QC)
#   - ribodetector_summary.txt                  ← summary of rRNA removal
#
# INSTALL (run once before first use):
#   conda activate rnaseq
#   pip install ribodetector
# =====================================================================


# ===== SLURM directives =====
#SBATCH --job-name=ribodetect
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
#SBATCH -o 01_processed/ribodetector_%j.out
#SBATCH -e 01_processed/ribodetector_%j.err

# Exit immediately if any command fails, if an unset variable is used,
# or if any command in a pipeline fails
set -euo pipefail


# =====================================================================
# LOAD CONFIG AND PARSE SPECIES CODE
# =====================================================================

# Source the shared config file (gives us BASE_DIR, PROCESSED_DIR, etc.)
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# Read species code from the command line (e.g., MF, DC, DG)
SPECIES_CODE="${1:-}"

# If no species code was provided, show usage instructions and exit
if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:  sbatch scripts/01_qc/run_ribodetector.sbatch <SPECIES_CODE>"
    echo ""
    echo "Valid codes: DC, DG, MF, SK, DCDG"
    echo "To add a new code, edit SECTION 4 in scripts/config.sh"
    exit 1
fi

# Look up species info from config.sh
# This sets: SPECIES_DIR, SPECIES_NAME, GENOME_TYPE
get_sample_info "${SPECIES_CODE}"


# =====================================================================
# BUILD PATHS AND SET PARAMETERS
# =====================================================================

# Where fastp put the cleaned reads (our INPUT)
# e.g., /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/
CLEAN_DIR="${PROCESSED_DIR}/${SPECIES_DIR}"

# Where to save rRNA-free reads (our OUTPUT)
# We make a subfolder called "ribofree" inside the processed directory
RIBO_OUTPUT="${PROCESSED_DIR}/${SPECIES_DIR}/ribofree"

# How many CPU threads to use (matches SLURM --cpus-per-task)
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- RiboDetector parameters ---
# READ_LENGTH: The average length of your reads AFTER fastp trimming.
#   Standard Illumina paired-end reads are usually 150bp.
#   If your reads are a different length, change this number.
READ_LENGTH=150

# ENSURE_MODE: How to handle paired-end reads where R1 and R2 disagree.
#   "rrna"   = only remove a pair if BOTH reads look like rRNA
#              (conservative — keeps more data, may leave trace rRNA)
#   "both"   = discard pairs where R1 and R2 disagree
#              (strict — removes more data but cleaner result)
#   "norrna" = only keep a pair if BOTH reads look like non-rRNA
#   "none"   = use the average probability of the pair (default)
ENSURE_MODE="rrna"


# =====================================================================
# SAFETY CHECKS
# =====================================================================

# Make sure the fastp output directory exists
if [[ ! -d "${CLEAN_DIR}" ]]; then
    echo "ERROR: Fastp output directory not found: ${CLEAN_DIR}"
    echo "Did you run fastp first?"
    echo "  sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
    exit 1
fi

# Count how many cleaned R1 files we have
SAMPLE_COUNT=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" | wc -l)

# If no cleaned files found, tell the user
if [[ "${SAMPLE_COUNT}" -eq 0 ]]; then
    echo "ERROR: No fastp-cleaned files found in ${CLEAN_DIR}"
    echo "Expected files like: *_1_clean.fq.gz and *_2_clean.fq.gz"
    echo ""
    echo "Did you run fastp first?"
    echo "  sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
    exit 1
fi


# =====================================================================
# PRINT JOB INFO
# =====================================================================

echo "============================================================"
echo "RiboDetector — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "============================================================"
echo "Host: $(hostname)  CPUs: ${THREADS}  Start: $(date)"
echo "Input:       ${CLEAN_DIR}"
echo "Output:      ${RIBO_OUTPUT}"
echo "Samples:     ${SAMPLE_COUNT}"
echo "Read length: ${READ_LENGTH}bp"
echo "Ensure mode: ${ENSURE_MODE}"
echo "============================================================"


# =====================================================================
# ACTIVATE CONDA AND CHECK RIBODETECTOR
# =====================================================================

# Create the output directory
mkdir -p "${RIBO_OUTPUT}"

# Activate conda environment
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate "${CONDA_ENV}"
else
    source ~/.bashrc
    conda activate "${CONDA_ENV}"
fi

# Make sure ribodetector is installed
if ! command -v ribodetector_cpu &> /dev/null; then
    echo "ERROR: ribodetector not found!"
    echo "Install it with:  pip install ribodetector"
    echo "  (make sure your conda env '${CONDA_ENV}' is active first)"
    exit 1
fi

echo "RiboDetector version: $(pip show ribodetector 2>/dev/null | grep Version || echo 'unknown')"
echo "------------------------------------------------------------"


# =====================================================================
# RUN RIBODETECTOR ON EACH SAMPLE
# =====================================================================
#
# How this loop works:
#   1. We find every *_1_clean.fq.gz file (the R1 read from fastp)
#   2. For each R1, we figure out the matching R2 file name
#   3. We run ribodetector_cpu to separate rRNA from non-rRNA reads
#   4. Non-rRNA reads go into *_ribofree.fq.gz (what you want to keep)
#   5. rRNA reads go into *_rrna.fq.gz (saved for QC, you can check them)
#   6. AFTER RiboDetector finishes, we count reads from the two output
#      files (ribofree + rrna) and add them to get the input total.
#      We do it this way instead of counting the big input file BEFORE
#      RiboDetector, because zcat on a large .fq.gz can take 5-10 min
#      of silent waiting with no log output — making it look stuck.
#
# We use ribodetector_cpu (not ribodetector) because HPC compute nodes
# usually don't have GPUs. If your cluster has GPU nodes, you could
# switch to the GPU version for faster processing.

PROCESSED_COUNT=0
FAILED_COUNT=0

# Start the summary file
SUMMARY_FILE="${RIBO_OUTPUT}/ribodetector_summary.txt"
echo "RiboDetector Summary — ${SPECIES_NAME} (${SPECIES_CODE})" > "${SUMMARY_FILE}"
echo "Date: $(date)" >> "${SUMMARY_FILE}"
echo "Read length: ${READ_LENGTH}bp | Ensure mode: ${ENSURE_MODE}" >> "${SUMMARY_FILE}"
echo "------------------------------------------------------------" >> "${SUMMARY_FILE}"
printf "%-20s %12s %12s %12s %8s\n" "Sample" "Input_Reads" "NonrRNA" "rRNA" "%rRNA" >> "${SUMMARY_FILE}"
echo "------------------------------------------------------------" >> "${SUMMARY_FILE}"

for R1_CLEAN in $(find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" | sort); do

    # --- Figure out the matching R2 file ---
    # If R1 is MFF1_1_clean.fq.gz, R2 is MFF1_2_clean.fq.gz
    R2_CLEAN="${R1_CLEAN/_1_clean.fq.gz/_2_clean.fq.gz}"

    # Extract the sample name (e.g., "MFF1" from "MFF1_1_clean.fq.gz")
    SAMPLE=$(basename "${R1_CLEAN}" "_1_clean.fq.gz")

    # Check that the R2 file exists
    if [[ ! -f "${R2_CLEAN}" ]]; then
        echo "WARNING: R2 not found for ${SAMPLE}. Expected: ${R2_CLEAN}"
        echo "  Skipping this sample."
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi

    # --- Set output file names ---
    # Non-rRNA reads (the GOOD reads you want to keep)
    OUT_R1="${RIBO_OUTPUT}/${SAMPLE}_1_ribofree.fq.gz"
    OUT_R2="${RIBO_OUTPUT}/${SAMPLE}_2_ribofree.fq.gz"

    # rRNA reads (saved so you can check how much rRNA was removed)
    RRNA_R1="${RIBO_OUTPUT}/${SAMPLE}_1_rrna.fq.gz"
    RRNA_R2="${RIBO_OUTPUT}/${SAMPLE}_2_rrna.fq.gz"

    # --- Skip if already processed ---
    if [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]]; then
        echo "Skipping ${SAMPLE} (already processed)"
        PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
        continue
    fi

    echo "------------------------------------------------------------"
    echo "Processing: ${SAMPLE}"
    echo "  R1 input: ${R1_CLEAN}"
    echo "  R2 input: ${R2_CLEAN}"

    # --- Run RiboDetector ---
    # ribodetector_cpu  = CPU version (use 'ribodetector' if you have GPU)
    # -t THREADS        = number of CPU threads to use
    # -l READ_LENGTH    = average read length (needed for the neural network)
    # -i R1 R2          = input paired-end files
    # -e ENSURE_MODE    = how to handle pairs where R1/R2 disagree
    # -o OUT_R1 OUT_R2  = output non-rRNA reads (what you keep)
    # -r RRNA_R1 RRNA_R2 = output rRNA reads (for QC checking)
    ribodetector_cpu \
        -t "${THREADS}" \
        -l "${READ_LENGTH}" \
        -i "${R1_CLEAN}" "${R2_CLEAN}" \
        -e "${ENSURE_MODE}" \
        -o "${OUT_R1}" "${OUT_R2}" \
        -r "${RRNA_R1}" "${RRNA_R2}"

    # --- Check if it worked ---
    if [[ -f "${OUT_R1}" ]] && [[ -f "${OUT_R2}" ]]; then
        echo "  Done: ${SAMPLE}"
        PROCESSED_COUNT=$((PROCESSED_COUNT + 1))

        # Count reads from the SMALLER output files (much faster than input)
        echo "  Counting output reads..."
        NONRRNA_READS=$(zcat "${OUT_R1}" | wc -l)
        NONRRNA_READS=$((NONRRNA_READS / 4))

        RRNA_READS=$(zcat "${RRNA_R1}" | wc -l)
        RRNA_READS=$((RRNA_READS / 4))

        # Total input = non-rRNA kept + rRNA removed
        INPUT_READS=$((NONRRNA_READS + RRNA_READS))

        # Calculate percentage of rRNA (using awk for decimal math)
        if [[ "${INPUT_READS}" -gt 0 ]]; then
            RRNA_PCT=$(awk "BEGIN {printf \"%.2f\", (${RRNA_READS}/${INPUT_READS})*100}")
        else
            RRNA_PCT="0.00"
        fi

        echo "  Input reads:   ${INPUT_READS}"
        echo "  Non-rRNA kept: ${NONRRNA_READS}"
        echo "  rRNA removed:  ${RRNA_READS} (${RRNA_PCT}%)"

        # Write to summary file
        printf "%-20s %12d %12d %12d %7s%%\n" \
            "${SAMPLE}" "${INPUT_READS}" "${NONRRNA_READS}" "${RRNA_READS}" "${RRNA_PCT}" \
            >> "${SUMMARY_FILE}"
    else
        echo "  FAILED: ${SAMPLE}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi

done


# =====================================================================
# SUMMARY
# =====================================================================

echo "------------------------------------------------------------" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Total processed: ${PROCESSED_COUNT}" >> "${SUMMARY_FILE}"
echo "Total failed:    ${FAILED_COUNT}" >> "${SUMMARY_FILE}"

echo ""
echo "============================================================"
echo "RiboDetector completed — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "End time: $(date)"
echo ""
echo "  Samples processed: ${PROCESSED_COUNT}"
echo "  Failed:            ${FAILED_COUNT}"
echo ""
echo "RESULTS in: ${RIBO_OUTPUT}/"
echo "  rRNA-free reads: *_ribofree.fq.gz  (use these for STAR)"
echo "  rRNA reads:      *_rrna.fq.gz      (check these for QC)"
echo "  Summary:         ribodetector_summary.txt"
echo ""
echo "NEXT STEPS:"
echo "  1. Check the summary to see how much rRNA was removed:"
echo "     cat ${RIBO_OUTPUT}/ribodetector_summary.txt"
echo ""
echo "  2. Re-run STAR alignment using the rRNA-free reads:"
echo "     (you will need to update STAR input to point to the ribofree/ folder)"
echo "============================================================"

# --- Copy log files into the species folder ---
cp "01_processed/ribodetector_${SLURM_JOB_ID}.out" "${RIBO_OUTPUT}/" 2>/dev/null || true
cp "01_processed/ribodetector_${SLURM_JOB_ID}.err" "${RIBO_OUTPUT}/" 2>/dev/null || true
