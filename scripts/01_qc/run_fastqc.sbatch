#!/bin/bash
# =====================================================================
# SCRIPT:  run_fastqc.sbatch
# PURPOSE: Run FastQC + MultiQC on FASTQ files for ONE species.
#          Works on either raw reads or fastp-cleaned reads.
#
# USAGE:
#   sbatch scripts/01_qc/run_fastqc.sbatch <SPECIES_CODE> [raw|clean|kraken|kraken-classified]
#
# EXAMPLES:
#   sbatch scripts/01_qc/run_fastqc.sbatch MF                      # raw (default)
#   sbatch scripts/01_qc/run_fastqc.sbatch MF raw                  # raw reads
#   sbatch scripts/01_qc/run_fastqc.sbatch MF clean                # after fastp
#   sbatch scripts/01_qc/run_fastqc.sbatch MF kraken               # kraken2-filtered (nutmeg reads)
#   sbatch scripts/01_qc/run_fastqc.sbatch MF kraken-classified    # kraken2-classified (contamination)
#
# INPUTS:
#   raw               → reads from 00_rawdata/<SPECIES_DIR>/
#   clean             → reads from 01_processed/<SPECIES_DIR>/*_clean.fq.gz
#   kraken            → reads from 01_processed/00_00_kraken/<SPECIES_CODE>*_filtered.fq.gz
#                       (unclassified = your nutmeg reads — feed these to Trinity)
#   kraken-classified → reads from 01_processed/00_00_kraken/<SPECIES_CODE>*_classified*.fq.gz
#                       (classified = contamination — useful for comparing vs filtered QC)
#
#   NOTE: both kraken stages check that files exist before running.
#         If they are missing, the job will stop with a helpful message.
#
# OUTPUTS (saved to 01_processed/<SPECIES_DIR>/):
#   fastqc_raw/  or  fastqc_clean/  or  fastqc_kraken/  or  fastqc_kraken_classified/
#   ├── fastqc_results/
#   └── multiqc_report.html
# =====================================================================
# SECTION 1: SLURM DIRECTIVES
# =====================================================================

# ===== SLURM directives =====
#SBATCH --job-name=fastqc
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 01_processed/fastqc_%j.out
#SBATCH -e 01_processed/fastqc_%j.err

# =====================================================================
# SECTION 2: SAFETY SETTINGS
# =====================================================================

set -euo pipefail
# ↑ This is a bash safety net. It means:
#   -e  = stop the script immediately if ANY command fails
#   -u  = treat unset variables as errors (catches typos)
#   -o pipefail = catch errors in piped commands (cmd1 | cmd2)
#   

# =====================================================================
# SECTION 3: LOAD CONFIG AND PARSE SPECIES CODE
# =====================================================================

# --- Load the shared config file ---
# This file defines all our directory paths (RAWDATA_DIR, PROCESSED_DIR, etc.)
# and the get_sample_info() function that maps species codes to folder names.
# ${SLURM_SUBMIT_DIR:-.} means "use SLURM_SUBMIT_DIR if set, otherwise use
# the current directory (.)"
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# --- Read arguments from the command line ---
# $1 = species code (required), e.g., MF, DC, DG
# $2 = data stage (optional), "raw" or "clean" — defaults to "raw"
SPECIES_CODE="${1:-}"
DATA_STAGE="${2:-raw}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:  sbatch run_fastqc.sbatch <SPECIES_CODE> [raw|clean|kraken|kraken-classified]"
    echo ""
    echo "  raw               = check raw reads (before fastp)        ← default"
    echo "  clean             = check cleaned reads (after fastp)"
    echo "  kraken            = check kraken2-filtered reads (your nutmeg reads)"
    echo "  kraken-classified = check kraken2-classified reads (contamination)"
    echo ""
    echo "Valid codes: DC, DG, MF, SK, DCDG"
    exit 1
fi

# Make sure DATA_STAGE is one of the four valid options
if [[ "${DATA_STAGE}" != "raw" ]] && \
   [[ "${DATA_STAGE}" != "clean" ]] && \
   [[ "${DATA_STAGE}" != "kraken" ]] && \
   [[ "${DATA_STAGE}" != "kraken-classified" ]]; then
    echo "ERROR: Second argument must be 'raw', 'clean', 'kraken', or 'kraken-classified'."
    echo "       Got '${DATA_STAGE}'"
    echo "Usage:  sbatch run_fastqc.sbatch ${SPECIES_CODE} [raw|clean|kraken|kraken-classified]"
    exit 1
fi

# --- Look up species info from config.sh ---
# This function reads the species code (e.g., "MF") and sets three variables:
#   SPECIES_DIR   = "00_3_MF"            (folder name under 00_rawdata/)
#   SPECIES_NAME  = "Myristica fragrans" (display name for logs)
#   GENOME_TYPE   = "nutmeg"             (which genome to use later)
get_sample_info "${SPECIES_CODE}"


# =====================================================================
# SECTION 4: BUILD FILE PATHS
# =====================================================================
# Choose input/output based on whether we're checking raw or clean data.
#
#   raw:   input = 00_rawdata/00_3_MF/         output = 01_processed/00_3_MF/fastqc_raw/
#   clean: input = 01_processed/00_3_MF/       output = 01_processed/00_3_MF/fastqc_clean/

if [[ "${DATA_STAGE}" == "raw" ]]; then
    # Raw reads live under 00_rawdata/
    SPECIES_FASTQ_DIR="${RAWDATA_DIR}/${SPECIES_DIR}"
elif [[ "${DATA_STAGE}" == "clean" ]]; then
    # Cleaned reads live under 01_processed/ (from fastp)
    SPECIES_FASTQ_DIR="${PROCESSED_DIR}/${SPECIES_DIR}"
elif [[ "${DATA_STAGE}" == "kraken" ]]; then
    # Kraken2-filtered reads (unclassified = nutmeg) live in the shared kraken folder.
    # We point here only so we can check files exist in SECTION 5.
    # After the check passes, SECTION 5 resets SPECIES_FASTQ_DIR to a temp
    # folder containing ONLY the *_filtered.fq.gz files.
    SPECIES_FASTQ_DIR="${PROCESSED_DIR}/00_00_kraken"
else
    # kraken-classified: the contamination reads (bacteria, fungi, viruses, etc.)
    # Same kraken directory — SECTION 5 will symlink only *_classified*.fq.gz.
    SPECIES_FASTQ_DIR="${PROCESSED_DIR}/00_00_kraken"
fi

# QC results go into a subfolder so no stage overwrites another.
# Note: "kraken-classified" becomes "fastqc_kraken-classified/" on disk.
QC_OUTPUT="${PROCESSED_DIR}/${SPECIES_DIR}/fastqc_${DATA_STAGE}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"


# =====================================================================
# SECTION 5: SAFETY CHECKS
# =====================================================================
# These checks catch common mistakes BEFORE the job wastes time running.

# --- Check: Does the input directory exist? ---
if [[ ! -d "${SPECIES_FASTQ_DIR}" ]]; then
    echo "ERROR: Input directory not found: ${SPECIES_FASTQ_DIR}"
    if [[ "${DATA_STAGE}" == "clean" ]]; then
        echo "Have you run fastp first?"
        echo "  sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
    elif [[ "${DATA_STAGE}" == "kraken" ]]; then
        echo "Have you run kraken2 first?"
        echo "  Expected kraken output directory: ${SPECIES_FASTQ_DIR}"
    fi
    exit 1
fi

# --- Check: Are there FASTQ files to process? ---
# For kraken / kraken-classified: we look only for the specific file type so
#   the Python pipeline doesn't accidentally grab the wrong set of reads.
#   We build a temp folder with symlinks to ONLY the files we want.
# For raw / clean: count all .fq.gz / .fastq.gz files in the directory.

if [[ "${DATA_STAGE}" == "kraken" ]]; then
    # Count only the *_filtered.fq.gz files for this species (nutmeg reads)
    FASTQ_COUNT=$(find "${SPECIES_FASTQ_DIR}" -maxdepth 1 -name "${SPECIES_CODE}*_filtered.fq.gz" | wc -l)

    if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
        echo "ERROR: No kraken2-filtered files found for species '${SPECIES_CODE}'"
        echo ""
        echo "  Looked in:  ${SPECIES_FASTQ_DIR}"
        echo "  Pattern:    ${SPECIES_CODE}*_filtered.fq.gz"
        echo ""
        echo "Expected files like:"
        echo "  ${SPECIES_FASTQ_DIR}/${SPECIES_CODE}F1_1_filtered.fq.gz"
        echo "  ${SPECIES_FASTQ_DIR}/${SPECIES_CODE}L1_1_filtered.fq.gz"
        echo ""
        echo "Have you run kraken2 for this species?"
        exit 1
    fi

    # Build a temp folder with symlinks to ONLY the filtered files.
    # A symlink is a shortcut — the file appears in two places without
    # copying any data. The $$ makes the folder name unique per job.
    KRAKEN_FILTERED_DIR="${SLURM_TMPDIR:-/tmp}/kraken_filtered_$$"
    mkdir -p "${KRAKEN_FILTERED_DIR}"
    echo "Creating temp folder with only filtered reads: ${KRAKEN_FILTERED_DIR}"

    for filtered_file in "${SPECIES_FASTQ_DIR}/${SPECIES_CODE}"*_filtered.fq.gz; do
        ln -sf "${filtered_file}" "${KRAKEN_FILTERED_DIR}/"
    done

    echo "Symlinked ${FASTQ_COUNT} filtered files into temp folder."
    SPECIES_FASTQ_DIR="${KRAKEN_FILTERED_DIR}"

elif [[ "${DATA_STAGE}" == "kraken-classified" ]]; then
    # Count only the *_classified*.fq.gz files for this species (contamination reads)
    # Note the naming pattern from kraken2: MFF1_classified_1.fq.gz, MFF1_classified_2.fq.gz
    FASTQ_COUNT=$(find "${SPECIES_FASTQ_DIR}" -maxdepth 1 -name "${SPECIES_CODE}*_classified*.fq.gz" | wc -l)

    if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
        echo "ERROR: No kraken2-classified files found for species '${SPECIES_CODE}'"
        echo ""
        echo "  Looked in:  ${SPECIES_FASTQ_DIR}"
        echo "  Pattern:    ${SPECIES_CODE}*_classified*.fq.gz"
        echo ""
        echo "Expected files like:"
        echo "  ${SPECIES_FASTQ_DIR}/${SPECIES_CODE}F1_classified_1.fq.gz"
        echo "  ${SPECIES_FASTQ_DIR}/${SPECIES_CODE}L1_classified_1.fq.gz"
        echo ""
        echo "Have you run kraken2 for this species?"
        exit 1
    fi

    # Build a temp folder with symlinks to ONLY the classified files.
    KRAKEN_CLASSIFIED_DIR="${SLURM_TMPDIR:-/tmp}/kraken_classified_$$"
    mkdir -p "${KRAKEN_CLASSIFIED_DIR}"
    echo "Creating temp folder with only classified (contamination) reads: ${KRAKEN_CLASSIFIED_DIR}"

    for classified_file in "${SPECIES_FASTQ_DIR}/${SPECIES_CODE}"*_classified*.fq.gz; do
        ln -sf "${classified_file}" "${KRAKEN_CLASSIFIED_DIR}/"
    done

    echo "Symlinked ${FASTQ_COUNT} classified files into temp folder."
    SPECIES_FASTQ_DIR="${KRAKEN_CLASSIFIED_DIR}"

else
    # For raw and clean stages: count all standard FASTQ files
    FASTQ_COUNT=$(find "${SPECIES_FASTQ_DIR}" -name "*.fq.gz" -o -name "*.fastq.gz" | wc -l)

    if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in ${SPECIES_FASTQ_DIR}"
        if [[ "${DATA_STAGE}" == "clean" ]]; then
            echo "Have you run fastp first?"
            echo "  sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
        fi
        exit 1
    fi
fi


# =====================================================================
# SECTION 6: PRINT JOB INFO
# =====================================================================
# This prints a summary so you can see what the job is doing when you
# check the log file (01_processed/qc_<jobid>.out).

echo "============================================================"
echo "FastQC — ${SPECIES_NAME} (${SPECIES_CODE}) [${DATA_STAGE} reads]"
echo "============================================================"
echo "Host: $(hostname)  CPUs: ${THREADS}  Start: $(date)"
echo "Data stage: ${DATA_STAGE}"
echo "Input:  ${SPECIES_FASTQ_DIR}"
echo "Output: ${QC_OUTPUT}"
echo "FASTQ files found: ${FASTQ_COUNT}"
echo "============================================================"


# =====================================================================
# SECTION 7: CREATE OUTPUT DIRECTORY
# =====================================================================

# mkdir = make directory
# -p    = create parent directories too if they don't exist, and
#         don't error if the directory already exists
mkdir -p "${QC_OUTPUT}"


# =====================================================================
# SECTION 8: ACTIVATE CONDA ENVIRONMENT
# =====================================================================
# Conda manages our Python packages and bioinformatics tools.
# We need to activate the "rnaseq" environment which has FastQC,
# MultiQC, and our pipeline code installed.

if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    # If .bashrc has issues, activate conda directly
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    # Normal path: load .bashrc (which sets up conda), then activate
    echo "Loading .bashrc successfully"
    source ~/.bashrc
    conda activate rnaseq
fi

# Double-check that conda activated correctly
conda info --envs | grep rnaseq || echo "ERROR: conda environment not activated"

echo "------------------------------------------------------------"


# =====================================================================
# SECTION 9: RUN FastQC + MultiQC
# =====================================================================
# This is the main part! We call our Python pipeline which:
#   1. Finds all .fq.gz files in the species directory
#   2. Runs FastQC on each file (checks quality of each FASTQ file)
#   3. Runs MultiQC to combine all FastQC reports into one summary

# First, move to the directory where our pipeline code lives
cd "${CODE_DIR}"

# Build the command as a string so we can print it before running it
# (this helps with debugging — you can copy-paste and run it manually)
#
# Flags explained:
#   --mode qc        = tell the pipeline to run quality control mode
#   --fastq-dir ...  = where to find the raw FASTQ files
#   --outdir ...     = where to save the QC results
#   --threads ...    = how many CPU cores to use (faster with more)
CMD="python -m rna_pipeline.cli \
  --mode qc \
  --fastq-dir \"${SPECIES_FASTQ_DIR}\" \
  --outdir \"${QC_OUTPUT}\" \
  --threads ${THREADS}"

echo "Running QC pipeline:"
echo "${CMD}"
echo "------------------------------------------------------------"

# eval runs the command stored in the CMD string
eval ${CMD}

echo "------------------------------------------------------------"


# =====================================================================
# SECTION 10: PRINT RESULTS SUMMARY
# =====================================================================

echo "============================================================"
echo "FastQC completed — ${SPECIES_NAME} (${SPECIES_CODE}) [${DATA_STAGE}]"
echo "End time: $(date)"
echo ""
echo "RESULTS:"
echo "  FastQC reports:   ${QC_OUTPUT}/fastqc_results/"
echo "  MultiQC summary:  ${QC_OUTPUT}/multiqc_report.html"
echo ""
if [[ "${DATA_STAGE}" == "raw" ]]; then
    echo "NEXT STEPS:"
    echo "1. Download and review: ${QC_OUTPUT}/multiqc_report.html"
    echo "2. Run fastp to clean the reads:"
    echo "   sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
    echo "3. Then re-run FastQC on cleaned reads to compare:"
    echo "   sbatch scripts/01_qc/run_fastqc.sbatch ${SPECIES_CODE} clean"
elif [[ "${DATA_STAGE}" == "clean" ]]; then
    echo "NEXT STEPS:"
    echo "1. Compare raw vs clean reports:"
    echo "   Raw:   ${PROCESSED_DIR}/${SPECIES_DIR}/fastqc_raw/multiqc_report.html"
    echo "   Clean: ${QC_OUTPUT}/multiqc_report.html"
    echo "2. If quality looks good, proceed to alignment (STAR)"
elif [[ "${DATA_STAGE}" == "kraken" ]]; then
    echo "NEXT STEPS:"
    echo "1. Download and review: ${QC_OUTPUT}/multiqc_report.html"
    echo "2. Run FastQC on the classified (contamination) reads to compare:"
    echo "   sbatch scripts/01_qc/run_fastqc.sbatch ${SPECIES_CODE} kraken-classified"
    echo "3. Then compare the two reports:"
    echo "   Nutmeg reads:       ${QC_OUTPUT}/multiqc_report.html"
    echo "   Contamination reads: ${PROCESSED_DIR}/${SPECIES_DIR}/fastqc_kraken-classified/multiqc_report.html"
    echo "4. If quality looks good, proceed to Trinity assembly"
else
    # kraken-classified stage
    echo "NEXT STEPS:"
    echo "1. Download both reports and compare them side by side:"
    echo ""
    echo "   Nutmeg reads (filtered):"
    echo "     ${PROCESSED_DIR}/${SPECIES_DIR}/fastqc_kraken/multiqc_report.html"
    echo ""
    echo "   Contamination reads (classified):"
    echo "     ${QC_OUTPUT}/multiqc_report.html"
    echo ""
    echo "   Things to look for in the comparison:"
    echo "     GC content:   contamination often shows a different GC% peak than nutmeg"
    echo "     Sequence length: contaminants may be shorter/different distribution"
    echo "     Duplication:  microbial reads can have very high duplication rates"
    echo "     Adapter content: should be similar if both came from the same fastp run"
fi
echo "============================================================"
