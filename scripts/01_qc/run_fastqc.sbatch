#!/bin/bash
# =====================================================================
# SCRIPT:  run_fastqc.sbatch
# PURPOSE: Run FastQC + MultiQC on the raw FASTQ files for ONE species
#
# WHAT IS FastQC?
#   FastQC reads each .fq.gz file and checks the quality of the
#   sequencing reads (base quality, adapter contamination, GC content,
#   etc.).  It produces one HTML report per file.
#
# WHAT IS MultiQC?
#   MultiQC takes ALL the individual FastQC reports and combines them
#   into ONE summary HTML page so you can compare all samples at once.
#
# USAGE:
#   sbatch scripts/01_qc/run_fastqc.sbatch <SPECIES_CODE>
#
# EXAMPLES:
#   sbatch scripts/01_qc/run_fastqc.sbatch DC    # D. carota
#   sbatch scripts/01_qc/run_fastqc.sbatch DG    # D. glaber
#   sbatch scripts/01_qc/run_fastqc.sbatch MF    # M. fragrans
#   sbatch scripts/01_qc/run_fastqc.sbatch SK    # Sukman samples
#
# INPUTS:
#   Raw .fq.gz files in: 00_rawdata/<SPECIES_DIR>/
#   (The folder name is looked up automatically from config.sh)
#
# OUTPUTS:
#   01_processed/qc_<CODE>/
#   ├── fastqc_results/        ← one HTML + ZIP per FASTQ file
#   └── multiqc_report.html    ← the combined report (open this one!)
#
# HOW TO MODIFY:
#   - To change resources (memory/time), edit the #SBATCH lines below
#   - To add a new species, edit SECTION 4 in scripts/config.sh
#   - To change the output folder, edit the QC_OUTPUT variable
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

# --- Read the species code from the command line ---
# When you run: sbatch run_fastqc.sbatch MF
# $1 = "MF" (the first argument after the script name)
# ${1:-} means "use $1 if provided, otherwise use empty string"
# This prevents the -u safety setting from crashing on missing args.
SPECIES_CODE="${1:-}"

# --- Check if the user actually provided a species code ---
# -z tests if a string is empty (zero length)
if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:  sbatch run_fastqc.sbatch <SPECIES_CODE>"
    echo ""
    echo "Valid codes: DC, DG, MF, SK, DCDG"
    echo "  DC   → D. carota ssp. maximus  (00_1_DC)"
    echo "  DG   → Daucus glaber            (00_2_DG)"
    echo "  MF   → Myristica fragrans       (00_3_MF)"
    echo "  SK   → Sukman samples            (00_4_Sukman)"
    echo "  DCDG → DC + DG combined          (00_4_DC_DG)"
    echo ""
    echo "To add a new code, edit SECTION 4 in scripts/config.sh"
    exit 1
    # ↑ exit 1 means "stop the script with an error status"
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

# --- Where to find the raw FASTQ files ---
# Combines the base rawdata path with the species-specific folder
# Example: /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_3_MF
SPECIES_FASTQ_DIR="${RAWDATA_DIR}/${SPECIES_DIR}"

# --- Where to save the QC output ---
# Saves into the same numbered folder as rawdata, but under 01_processed/
# This keeps input and output folder names matching:
#   Input:  00_rawdata/00_3_MF/       ← raw FASTQ files
#   Output: 01_processed/00_3_MF/     ← QC reports go here
QC_OUTPUT="${PROCESSED_DIR}/${SPECIES_DIR}"

# --- How many CPU threads to use ---
# Uses the number of CPUs requested in the SLURM directives above.
# If running outside SLURM (for testing), defaults to 8.
THREADS="${SLURM_CPUS_PER_TASK:-8}"


# =====================================================================
# SECTION 5: SAFETY CHECKS
# =====================================================================
# These checks catch common mistakes BEFORE the job wastes time running.

# --- Check 1: Does the input directory exist? ---
# -d tests if a directory exists
if [[ ! -d "${SPECIES_FASTQ_DIR}" ]]; then
    echo "ERROR: Input directory not found: ${SPECIES_FASTQ_DIR}"
    echo ""
    echo "This means the folder '${SPECIES_DIR}' doesn't exist under"
    echo "${RAWDATA_DIR}/"
    echo ""
    echo "To fix: check that your raw data folder name matches what's"
    echo "in config.sh (SECTION 4, get_sample_info function)."
    echo ""
    echo "Run this to see what folders exist:"
    echo "  ls ${RAWDATA_DIR}/"
    exit 1
fi

# --- Check 2: Are there actually FASTQ files in that directory? ---
# find looks for files ending in .fq.gz or .fastq.gz
# wc -l counts how many files were found
FASTQ_COUNT=$(find "${SPECIES_FASTQ_DIR}" -name "*.fq.gz" -o -name "*.fastq.gz" | wc -l)

if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in ${SPECIES_FASTQ_DIR}"
    echo ""
    echo "Expected files like: MFF1_1.fq.gz, MFF1_2.fq.gz, etc."
    echo ""
    echo "Run this to check what's in the folder:"
    echo "  ls ${SPECIES_FASTQ_DIR}/"
    exit 1
fi


# =====================================================================
# SECTION 6: PRINT JOB INFO
# =====================================================================
# This prints a summary so you can see what the job is doing when you
# check the log file (01_processed/qc_<jobid>.out).

echo "============================================================"
echo "FastQC Quality Control — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-16G}"
echo "Start time: $(date)"
echo ""
echo "Species:    ${SPECIES_NAME}"
echo "Code:       ${SPECIES_CODE}"
echo "Input dir:  ${SPECIES_FASTQ_DIR}"
echo "Output dir: ${QC_OUTPUT}"
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
echo "Quality control completed successfully!"
echo "Species: ${SPECIES_NAME} (${SPECIES_CODE})"
echo "End time: $(date)"
echo ""
echo "RESULTS:"
echo "  Individual FastQC reports: ${QC_OUTPUT}/fastqc_results/"
echo "  MultiQC summary report:    ${QC_OUTPUT}/multiqc_report.html"
echo "  Job logs:                  01_processed/qc_${SLURM_JOB_ID}.out"
echo ""
echo "NEXT STEPS:"
echo "1. Download the MultiQC report to your computer:"
echo "   scp daisycortesj@owl1:${QC_OUTPUT}/multiqc_report.html ."
echo "2. Open multiqc_report.html in your web browser"
echo "3. Look for:"
echo "   - Per base quality scores (should be >28, ideally >30)"
echo "   - Adapter contamination (should be <5%)"
echo "   - Per sequence GC content (should match ${SPECIES_NAME})"
echo "4. If quality is good → proceed to fastp trimming:"
echo "   sbatch scripts/01_qc/run_fastp.sbatch ${SPECIES_CODE}"
echo "5. If quality issues found → consider adjusting fastp parameters"
echo "============================================================"
