#!/bin/bash
# =============================================================================
# SCRIPT:  run_busco.sbatch
# WHAT IT DOES:
#   Runs BUSCO (Benchmarking Universal Single-Copy Orthologs) on your
#   pooled Trinity transcriptome to measure how complete the assembly is.
#
# WHAT BUSCO MEASURES:
#   BUSCO compares your assembly against a set of genes that are expected
#   to be present in nearly every plant genome (embryophyta_odb10 = 1,614 genes).
#   It reports what percentage of those genes it found in your assembly:
#
#     Complete (C)    — found once, full length  → good
#     Duplicated (D)  — found more than once     → expected for transcriptomes
#                       (Trinity makes multiple isoforms per gene)
#     Fragmented (F)  — found but incomplete     → some assembly issues
#     Missing (M)     — not found at all         → potential gaps
#
#   A good plant transcriptome typically shows:
#     C > 70%  (ideally > 85%)
#     High D is NORMAL for transcriptomes (not a problem like it is for genomes)
#     M < 20%
#
# INPUT:
#   Pooled Trinity assembly:
#     01_processed/00_4_MF_trinity/trinity_out.Trinity.fasta
#
# OUTPUTS:
#   06_analysis/busco_MF/                  ← main output folder
#   06_analysis/busco_MF/short_summary*    ← the key result file (read this first)
#   06_analysis/busco_MF/full_table.tsv    ← gene-by-gene details
#
# HOW TO SUBMIT:
#   sbatch scripts/04_busco/run_busco.sbatch
#
# BEFORE SUBMITTING — install BUSCO if not already done:
#   See the install instructions at the bottom of this file header.
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


# ---- SLURM cluster resource requests ----
# BUSCO is faster than Trinity but still needs decent resources.
# 16 CPUs and 32 GB RAM is enough for a transcriptome this size.
#SBATCH --job-name=busco_MF
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o busco_MF_%j.out
#SBATCH -e busco_MF_%j.err

# Stop the script immediately if any command fails
set -euo pipefail

#showing a new change

# Load shared folder paths (gives us BASE_DIR, ANALYSIS_DIR, PROCESSED_DIR, etc.)
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"


# =============================================================================
# PATHS — everything in one place, easy to update
# =============================================================================

# Input: your pooled Trinity FASTA (produced by trinity_pooled.sbatch)
ASSEMBLY="${PROCESSED_DIR}/00_4_MF_trinity/trinity_out.Trinity.fasta"

# Output folder for BUSCO results
BUSCO_OUT_DIR="${ANALYSIS_DIR}/busco_MF"

# Where BUSCO stores the lineage database after downloading it.
# This is saved to 04_reference/ so it can be reused for future runs.
# BUSCO auto-downloads embryophyta_odb10 (~700 MB) the first time.
# On subsequent runs it reads from this folder (much faster).
BUSCO_LINEAGE_DIR="${BASE_DIR}/04_reference/busco_lineages"

# The lineage database to use — embryophyta_odb10 covers land plants (1,614 genes)
LINEAGE="embryophyta_odb10"

# Number of CPUs — must match --cpus-per-task above
THREADS="${SLURM_CPUS_PER_TASK:-16}"


# =============================================================================
# JOB HEADER
# =============================================================================
echo "========================================================"
echo "BUSCO — Transcriptome Completeness Assessment"
echo "Species:  Myristica fragrans (MF)"
echo "Lineage:  ${LINEAGE}"
echo "Job ID:   ${SLURM_JOB_ID:-N/A}"
echo "Host:     $(hostname)"
echo "Started:  $(date)"
echo "========================================================"


# =============================================================================
# ACTIVATE CONDA
# BUSCO is installed in your rnaseq conda environment.
# =============================================================================
source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

# Confirm BUSCO is available
if ! command -v busco &> /dev/null; then
    echo ""
    echo "ERROR: busco not found in your conda environment."
    echo ""
    echo "Install it with:"
    echo "  conda activate rnaseq"
    echo "  conda install -c conda-forge -c bioconda busco"
    echo ""
    echo "Then resubmit this script."
    exit 1
fi

echo "BUSCO version: $(busco --version 2>&1)"
echo "------------------------------------------------------------"


# =============================================================================
# CHECK BUSCO'S UNDERLYING DEPENDENCIES
#
# BUSCO does not run on its own — it calls several other tools internally.
# If any of these are missing, BUSCO will either crash mid-run or give wrong
# results. We check for all of them here and give clear install instructions
# if anything is missing.
#
# Tool         | What BUSCO uses it for
# -------------|------------------------------------------------------
# Python 3     | BUSCO itself is written in Python
# HMMER        | Searches for BUSCO gene profiles using Hidden Markov Models
# AUGUSTUS     | Predicts gene structures (used in genome mode)
# BLAST+       | Finds genomic regions where BUSCO genes might live (tBLASTn)
# SEPP         | Used by the auto-lineage selection pipeline
# =============================================================================

echo ""
echo "Checking BUSCO dependencies..."

# A helper flag — if any tool is missing we flip this to 1 and exit at the end.
# This way you see ALL missing tools at once instead of one at a time.
MISSING_DEPS=0

# ---- Python 3 ----
if command -v python3 &> /dev/null; then
    echo "  ✓ Python 3 found: $(python3 --version 2>&1)"
else
    echo "  ✗ Python 3 NOT found."
    echo "    Install with:  conda install -c conda-forge python=3"
    MISSING_DEPS=1
fi

# ---- HMMER (the 'hmmsearch' command is the one BUSCO actually calls) ----
if command -v hmmsearch &> /dev/null; then
    echo "  ✓ HMMER found: $(hmmsearch -h 2>&1 | head -2 | tail -1 | xargs)"
else
    echo "  ✗ HMMER NOT found."
    echo "    Install with:  conda install -c bioconda hmmer"
    MISSING_DEPS=1
fi

# ---- AUGUSTUS (used in genome mode; also needed for some transcriptome runs) ----
if command -v augustus &> /dev/null; then
    echo "  ✓ AUGUSTUS found: $(augustus --version 2>&1 | head -1)"
else
    echo "  ✗ AUGUSTUS NOT found."
    echo "    Install with:  conda install -c bioconda augustus"
    echo "    Note: AUGUSTUS is optional for transcriptome mode but"
    echo "          required for genome mode."
    MISSING_DEPS=1
fi

# ---- BLAST+ (the 'tblastn' command is the one BUSCO calls) ----
if command -v tblastn &> /dev/null; then
    echo "  ✓ BLAST+ found: $(tblastn -version 2>&1 | head -1)"
else
    echo "  ✗ BLAST+ NOT found."
    echo "    Install with:  conda install -c bioconda blast"
    MISSING_DEPS=1
fi

# ---- SEPP (used only when --auto-lineage is enabled) ----
# SEPP is a placement tool for automatic lineage detection.
# It is only strictly required if you use --auto-lineage or --auto-lineage-euk.
# We just warn here instead of failing, since we use a fixed lineage above.
if command -v run_sepp.py &> /dev/null || python3 -c "import sepp" &> /dev/null 2>&1; then
    echo "  ✓ SEPP found"
else
    echo "  ~ SEPP not found (optional — only needed for --auto-lineage mode)."
    echo "    We are using a fixed lineage (${LINEAGE}), so this is OK."
    echo "    Install with:  conda install -c bioconda sepp    (if you need auto-lineage)"
fi

# ---- Did any required dependency fail? ----
if [[ "${MISSING_DEPS}" -eq 1 ]]; then
    echo ""
    echo "ERROR: One or more required BUSCO dependencies are missing."
    echo "Install the missing tools (commands shown above), then resubmit."
    echo ""
    echo "Easiest fix — install everything at once inside your conda environment:"
    echo "  conda activate rnaseq"
    echo "  conda install -c conda-forge -c bioconda busco hmmer blast augustus sepp"
    exit 1
fi

echo "All required dependencies found."
echo "------------------------------------------------------------"


# =============================================================================
# CHECK THAT THE TRINITY ASSEMBLY EXISTS
# =============================================================================
echo ""
echo "Checking input file..."
if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Trinity assembly not found at:"
    echo "  ${ASSEMBLY}"
    echo ""
    echo "Make sure trinity_pooled.sbatch has finished successfully before"
    echo "running this script."
    exit 1
fi

CONTIG_COUNT=$(grep -c ">" "${ASSEMBLY}" || echo "0")
echo "✓ Assembly found: ${CONTIG_COUNT} transcripts"
echo "  ${ASSEMBLY}"
echo "------------------------------------------------------------"


# =============================================================================
# CREATE DIRECTORIES
# =============================================================================
mkdir -p "${BUSCO_OUT_DIR}"
mkdir -p "${BUSCO_LINEAGE_DIR}"


# =============================================================================
# RUN BUSCO
# =============================================================================
echo ""
echo "========================================================"
echo "Running BUSCO..."
echo "========================================================"

# The BUSCO command:
#   -i  → input FASTA file (your Trinity assembly)
#   -l  → lineage database (embryophyta_odb10 = land plants)
#   -o  → output folder name
#   --out_path → where to put that output folder
#   -m  → mode: "transcriptome" (not genome or protein)
#   -c  → number of CPUs
#   --download_path → where to save/find the lineage database
#
# NOTE: The first time you run this, BUSCO downloads embryophyta_odb10
# (~700 MB) to BUSCO_LINEAGE_DIR. This can take 5–10 minutes.
# Subsequent runs skip the download and start immediately.

BUSCO_CMD="busco \
    -i ${ASSEMBLY} \
    -l ${LINEAGE} \
    -o busco_MF_results \
    --out_path ${BUSCO_OUT_DIR} \
    -m transcriptome \
    -c ${THREADS} \
    --download_path ${BUSCO_LINEAGE_DIR} \
    --force"
# --force overwrites a previous incomplete run in the output folder.
# Remove --force if you want BUSCO to error instead of overwriting.

echo "Command:"
echo "${BUSCO_CMD}"
echo "------------------------------------------------------------"

if eval ${BUSCO_CMD} 2>&1; then
    echo ""
    echo "✓ BUSCO finished successfully"
else
    BUSCO_EXIT=$?
    echo ""
    echo "✗ BUSCO failed (exit code: ${BUSCO_EXIT})"
    echo "  Check: ${BUSCO_OUT_DIR}/busco_MF_results/logs/busco.log"
    exit ${BUSCO_EXIT}
fi


# =============================================================================
# PRINT THE SUMMARY RESULT
# The short_summary file is the key output — it shows C/D/F/M percentages.
# =============================================================================
echo ""
echo "========================================================"
echo "BUSCO Summary"
echo "========================================================"

SUMMARY_FILE=$(find "${BUSCO_OUT_DIR}/busco_MF_results" \
    -name "short_summary*.txt" 2>/dev/null | head -1)

if [[ -f "${SUMMARY_FILE}" ]]; then
    cat "${SUMMARY_FILE}"
else
    echo "Summary file not found — check the output folder:"
    echo "  ${BUSCO_OUT_DIR}/busco_MF_results/"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ COMPLETE"
echo "========================================================"
echo "Results folder:  ${BUSCO_OUT_DIR}/busco_MF_results/"
echo "Summary:         ${SUMMARY_FILE:-see results folder}"
echo "Full table:      ${BUSCO_OUT_DIR}/busco_MF_results/full_table.tsv"
echo ""
echo "What to look for:"
echo "  C > 70% (ideally 85%+) = good assembly"
echo "  High D% is NORMAL for transcriptomes (isoforms inflate duplicates)"
echo "  M > 30% = possible assembly gaps worth investigating"
echo ""
echo "Next step: Kraken2 decontamination on Trinity.fasta"
echo "Finished: $(date)"
echo "========================================================"

exit 0
