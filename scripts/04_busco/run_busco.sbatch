#!/bin/bash
# =============================================================================
# SCRIPT:  run_busco.sbatch
# WHAT IT DOES:
#   Runs BUSCO (Benchmarking Universal Single-Copy Orthologs) on your
#   pooled Trinity transcriptome to measure how complete the assembly is.
#
# WHAT BUSCO MEASURES:
#   BUSCO compares your assembly against a set of genes expected to be present
#   in nearly every member of a given taxonomic group, then reports what
#   percentage it found:
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
# WHICH LINEAGE DATABASE GETS USED:
#   The lineage is chosen automatically from the species code you pass on the
#   command line (which maps to GENOME_TYPE in scripts/config.sh):
#
#     Species code   GENOME_TYPE   Lineage             # of genes
#     ------------   -----------   ------------------  ----------
#     MF             nutmeg        embryophyta_odb10   1,614
#     DC, DG, SK,    carrot        eudicots_odb10      ~2,326
#       DCDG
#
#   Why different lineages?
#     - Nutmeg (Myristica fragrans) is a basal angiosperm — not a eudicot —
#       so we use the broader "land plants" set (embryophyta).
#     - Carrot (Daucus carota) IS a eudicot, so we can use the more specific
#       eudicots set (more genes = stricter quality check).
#
# INPUT:
#   Pooled Trinity assembly (produced by trinity_pooled.sbatch):
#     MF  →  01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
#     DC  →  01_processed/00_1_DC_trinity/DC_trinity_pooled.Trinity.fasta
#     DG  →  01_processed/00_2_DG_trinity/DG_trinity_pooled.Trinity.fasta
#     etc.
#
# OUTPUTS:
#   01_processed/00_5_BUSCO/busco_<CODE>/                  ← main output folder
#   01_processed/00_5_BUSCO/busco_<CODE>/short_summary*    ← the key result file
#   01_processed/00_5_BUSCO/busco_<CODE>/full_table.tsv    ← gene-by-gene details
#
#   SLURM logs (.out / .err) also go in 01_processed/00_5_BUSCO/.
#
# HOW TO SUBMIT:
#   sbatch scripts/04_busco/run_busco.sbatch MF      # nutmeg
#   sbatch scripts/04_busco/run_busco.sbatch DC      # carrot (D. carota)
#   sbatch scripts/04_busco/run_busco.sbatch DG      # carrot (D. glaber)
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


# ---- SLURM cluster resource requests ----
# BUSCO is faster than Trinity but still needs decent resources.
# 16 CPUs and 32 GB RAM is enough for a transcriptome this size.
#SBATCH --job-name=busco
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
# All BUSCO outputs and SLURM logs go to 01_processed/00_5_BUSCO/.
# NOTE: SBATCH directives are parsed by SLURM BEFORE bash runs, so we can't use
# variables here — the path must be hardcoded. If you change BASE_DIR in
# config.sh, also update this path.
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_5_BUSCO
#SBATCH -o busco_%j.out
#SBATCH -e busco_%j.err

# Stop the script immediately if any command fails
set -euo pipefail

# Load shared folder paths and the get_sample_info() function.
# Use this script's location so config.sh loads even if you sbatch from another folder.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"


# =============================================================================
# READ THE SPECIES CODE FROM THE COMMAND LINE
# =============================================================================
# $1 = species code (required). Valid codes are defined in config.sh:
#      MF, DC, DG, SK, DCDG.

SPECIES_CODE="${1:-}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/04_busco/run_busco.sbatch <SPECIES_CODE>"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_busco/run_busco.sbatch MF      # nutmeg"
    echo "  sbatch scripts/04_busco/run_busco.sbatch DC      # carrot (D. carota)"
    echo "  sbatch scripts/04_busco/run_busco.sbatch DG      # carrot (D. glaber)"
    exit 1
fi

# get_sample_info populates SPECIES_DIR, SPECIES_NAME, GENOME_TYPE.
# It exits non-zero if the species code isn't recognized — set -e will catch that.
get_sample_info "${SPECIES_CODE}"

# Normalize the species code to uppercase so paths/names are consistent.
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# PICK LINEAGE + TRINITY INPUT PATH BASED ON SPECIES CODE
# =============================================================================
# Why a case statement?
#   It keeps every species-specific setting (lineage + Trinity FASTA path) in
#   ONE place. To add a new species, just add a new case block here.
#
# Why not auto-derive the Trinity path from SPECIES_DIR in config.sh?
#   Because SPECIES_DIR (e.g., 00_3_MF) refers to the RAW DATA folder, while
#   Trinity outputs live in a different folder under 01_processed/ that doesn't
#   share the same number prefix (e.g., 00_4_MF_trinity). So we set the
#   Trinity path explicitly per species.

case "${SPECIES_CODE_UPPER}" in
    MF)
        # Nutmeg is a basal angiosperm (not a eudicot, not a monocot),
        # so we stay broad and use the land-plant set.
        LINEAGE="embryophyta_odb10"
        LINEAGE_NOTE="land plants (1,614 genes) — used for basal angiosperms like nutmeg"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_MF_trinity"
        ASSEMBLY="${TRINITY_DIR}/MF_trinity_pooled.Trinity.fasta"
        TRINITY_WORK_DIR="${TRINITY_DIR}/MF_trinity_pooled"
        ;;
    DC)
        # Daucus carota (carrot) is a eudicot — use the tighter eudicot set.
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes) — used for carrots and other eudicot plants"
        TRINITY_DIR="${PROCESSED_DIR}/00_1_DC_trinity"
        ASSEMBLY="${TRINITY_DIR}/DC_trinity_pooled.Trinity.fasta"
        TRINITY_WORK_DIR="${TRINITY_DIR}/DC_trinity_pooled"
        ;;
    DG)
        # Daucus glaber (carrot relative) is also a eudicot.
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes) — used for carrots and other eudicot plants"
        TRINITY_DIR="${PROCESSED_DIR}/00_2_DG_trinity"
        ASSEMBLY="${TRINITY_DIR}/DG_trinity_pooled.Trinity.fasta"
        TRINITY_WORK_DIR="${TRINITY_DIR}/DG_trinity_pooled"
        ;;
    SK)
        # Sukman samples — carrot, eudicot.
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes) — used for carrots and other eudicot plants"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_Sukman_trinity"
        ASSEMBLY="${TRINITY_DIR}/SK_trinity_pooled.Trinity.fasta"
        TRINITY_WORK_DIR="${TRINITY_DIR}/SK_trinity_pooled"
        ;;
    DCDG)
        # Combined DC + DG carrot assembly — also eudicot.
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes) — used for carrots and other eudicot plants"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_DC_DG_trinity"
        ASSEMBLY="${TRINITY_DIR}/DCDG_trinity_pooled.Trinity.fasta"
        TRINITY_WORK_DIR="${TRINITY_DIR}/DCDG_trinity_pooled"
        ;;
    *)
        # Safe fallback if a new code is added to config.sh but not here.
        # Pick lineage from GENOME_TYPE, but we can't guess the Trinity path,
        # so error out clearly.
        echo "ERROR: No Trinity input path configured for species '${SPECIES_CODE_UPPER}'."
        echo ""
        echo "Add a new case block to run_busco.sbatch for this species,"
        echo "specifying:"
        echo "  - LINEAGE             (which BUSCO lineage to use)"
        echo "  - TRINITY_DIR         (folder containing the assembly)"
        echo "  - ASSEMBLY            (full path to the .Trinity.fasta file)"
        echo "  - TRINITY_WORK_DIR    (Trinity's working folder, same prefix)"
        exit 1
        ;;
esac


# =============================================================================
# REMAINING PATHS — same for all species
# =============================================================================

# Top-level folder where ALL BUSCO results for ALL species live.
# Sits next to the Trinity output folders under 01_processed/.
BUSCO_BASE_DIR="${PROCESSED_DIR}/00_5_BUSCO"

# Output folder for BUSCO results — one subfolder per species so runs don't
# overwrite each other (e.g., busco_MF/, busco_DC/, busco_DG/).
BUSCO_OUT_DIR="${BUSCO_BASE_DIR}/busco_${SPECIES_CODE_UPPER}"

# Name of the inner results folder BUSCO creates
BUSCO_RESULTS_NAME="busco_${SPECIES_CODE_UPPER}_results"

# Where BUSCO stores the lineage database after downloading it.
# This is saved to 04_reference/ so it can be reused for future runs.
# BUSCO auto-downloads the lineage (~700 MB) the first time.
# On subsequent runs it reads from this folder (much faster).
BUSCO_LINEAGE_DIR="${BASE_DIR}/04_reference/busco_lineages"

# Number of CPUs — must match --cpus-per-task above
THREADS="${SLURM_CPUS_PER_TASK:-16}"


# =============================================================================
# JOB HEADER
# =============================================================================
echo "========================================================"
echo "BUSCO — Transcriptome Completeness Assessment"
echo "Species code:  ${SPECIES_CODE_UPPER}"
echo "Species:       ${SPECIES_NAME}"
echo "Genome type:   ${GENOME_TYPE}"
echo "Lineage:       ${LINEAGE}"
echo "               (${LINEAGE_NOTE})"
echo "Job ID:        ${SLURM_JOB_ID:-N/A}"
echo "Host:          $(hostname)"
echo "Started:       $(date)"
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
    if [[ -d "${TRINITY_WORK_DIR}" ]]; then
        echo "Found Trinity working directory (intermediate files only):"
        echo "  ${TRINITY_WORK_DIR}/"
        echo ""
        echo "The final FASTA should be one level UP (same folder as the work dir):"
        echo "  ${ASSEMBLY}"
        echo ""
        echo "If that file is missing, trinity_pooled.sbatch may not have finished yet."
    else
        echo "Expected folder: ${TRINITY_DIR}/"
        echo "Run trinity_pooled.sbatch for ${SPECIES_CODE_UPPER} first, then resubmit BUSCO."
    fi
    exit 1
fi

CONTIG_COUNT=$(grep -c ">" "${ASSEMBLY}" || echo "0")
echo "✓ Assembly found: ${CONTIG_COUNT} transcripts"
echo "  ${ASSEMBLY}"
echo "------------------------------------------------------------"


# =============================================================================
# CREATE DIRECTORIES
# =============================================================================
# mkdir -p creates parent folders too, so this also creates 00_5_BUSCO/
# the first time you run BUSCO on any species.
mkdir -p "${BUSCO_BASE_DIR}"
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
#   -l  → lineage database (chosen above based on species)
#   -o  → output folder name (named after the species code)
#   --out_path → where to put that output folder
#   -m  → mode: "transcriptome" (not genome or protein)
#   -c  → number of CPUs
#   --download_path → where to save/find the lineage database
#
# NOTE: The first time you run this for a given lineage, BUSCO downloads
# that lineage (~700 MB) to BUSCO_LINEAGE_DIR. This can take 5–10 minutes.
# Subsequent runs skip the download and start immediately.

BUSCO_CMD="busco \
    -i ${ASSEMBLY} \
    -l ${LINEAGE} \
    -o ${BUSCO_RESULTS_NAME} \
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
    echo "  Check: ${BUSCO_OUT_DIR}/${BUSCO_RESULTS_NAME}/logs/busco.log"
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

SUMMARY_FILE=$(find "${BUSCO_OUT_DIR}/${BUSCO_RESULTS_NAME}" \
    -name "short_summary*.txt" 2>/dev/null | head -1)

if [[ -f "${SUMMARY_FILE}" ]]; then
    cat "${SUMMARY_FILE}"
else
    echo "Summary file not found — check the output folder:"
    echo "  ${BUSCO_OUT_DIR}/${BUSCO_RESULTS_NAME}/"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ COMPLETE"
echo "========================================================"
echo "Species:         ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Lineage used:    ${LINEAGE}"
echo "Results folder:  ${BUSCO_OUT_DIR}/${BUSCO_RESULTS_NAME}/"
echo "Summary:         ${SUMMARY_FILE:-see results folder}"
echo "Full table:      ${BUSCO_OUT_DIR}/${BUSCO_RESULTS_NAME}/full_table.tsv"
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
