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
#
# WHERE DOES THIS FIT IN THE PIPELINE?
#   Raw reads → fastp (adapter/quality trim) → RiboDetector (rRNA removal) → STAR alignment
#   Step 1        Step 2 (already done)            Step 3 (THIS SCRIPT)         Step 4
#
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
#   - *_ribodetector.log                        ← per-sample progress log
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
#SBATCH --time=72:00:00
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

# How many CPU threads to use
# Set to 4 — a safe value that avoids memory issues on shared nodes
THREADS=4

# SAMPLES: The list of sample names to process.
# Add or remove names here if your set of samples changes.
# These names must match the beginning of your fastp output files.
# e.g., MFF1 → MFF1_1_clean.fq.gz and MFF1_2_clean.fq.gz
SAMPLES="MFF1 MFF2 MFF3 MFL1 MFL2 MFL3"

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

# CHUNK_SIZE: Number of reads loaded into memory at a time.
#   256 is the value recommended in the RiboDetector docs.
#   Without this flag RiboDetector may silently hang on large files.
#   Lower this number (e.g. 128) if you get out-of-memory errors.
CHUNK_SIZE=256


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

# Count samples from the SAMPLES list (just count the words)
SAMPLE_COUNT=$(echo "${SAMPLES}" | wc -w)


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
echo "Chunk size:  ${CHUNK_SIZE}"
echo "============================================================"


# =====================================================================
# ACTIVATE CONDA AND CHECK RIBODETECTOR
# =====================================================================

# Create the output directory
mkdir -p "${RIBO_OUTPUT}"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
    source /etc/profile.d/conda.sh 2>/dev/null
conda activate "${CONDA_ENV}"

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
#   1. We loop over each sample name in the SAMPLES list above
#   2. We build the R1 and R2 input paths using the sample name
#   3. We run ribodetector_cpu to separate rRNA from non-rRNA reads
#   4. Non-rRNA reads go into *_ribofree.fq.gz (what you want to keep)
#   5. rRNA reads go into *_rrna.fq.gz (saved for QC, you can check them)
#
# We use ribodetector_cpu (not ribodetector) because HPC compute nodes
# usually don't have GPUs. If your cluster has GPU nodes, you could
# switch to the GPU version for faster processing.

PROCESSED_COUNT=0
FAILED_COUNT=0

# NOTE: MFF1 was already processed manually — starting from MFF2
# -t 4           = 4 CPU threads (safe value for shared nodes)
# -l 150         = read length 150bp (standard Illumina paired-end)
# -e rrna        = only remove pair if BOTH reads look like rRNA (conservative)
# --chunk_size   = 256 reads loaded at a time — REQUIRED to prevent silent hang
# -i             = input R1 and R2 cleaned reads from fastp
# -o             = output rRNA-free reads (use these for STAR)
# -r             = output rRNA reads (saved for QC)
for SAMPLE in MFF1 MFL1 MFL2 MFL3; do
    echo "========================================"
    echo "Starting: ${SAMPLE}  $(date)"
    echo "========================================"
    ribodetector_cpu \
        -t 4 \
        -l 150 \
        -i /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/${SAMPLE}_1_clean.fq.gz \
           /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/${SAMPLE}_2_clean.fq.gz \
        -e rrna \
        --chunk_size 256 \
        -o /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/ribofree/${SAMPLE}_1_ribofree.fq.gz \
           /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/ribofree/${SAMPLE}_2_ribofree.fq.gz \
        -r /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/ribofree/${SAMPLE}_1_rrna.fq.gz \
           /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/ribofree/${SAMPLE}_2_rrna.fq.gz
    echo "Done: ${SAMPLE}  $(date)"
    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
done


# =====================================================================
# SUMMARY
# =====================================================================

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
echo "  Per-sample logs: *_ribodetector.log"
echo ""
echo "NEXT STEPS:"
echo "  1. Check per-sample logs:"
echo "     cat ${RIBO_OUTPUT}/MFF1_ribodetector.log"
echo ""
echo "  2. Re-run STAR alignment using the rRNA-free reads:"
echo "     (update STAR input to point to the ribofree/ folder)"
echo "============================================================"

# --- Copy log files into the species folder ---
cp "01_processed/ribodetector_${SLURM_JOB_ID}.out" "${RIBO_OUTPUT}/" 2>/dev/null || true
cp "01_processed/ribodetector_${SLURM_JOB_ID}.err" "${RIBO_OUTPUT}/" 2>/dev/null || true
