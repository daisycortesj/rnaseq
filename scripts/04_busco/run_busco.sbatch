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
# INPUT (two options — pick one per job):
#
#   1. Pooled Trinity (default):
#     MF  →  01_processed/00_4_MF_trinity/MF_trinity_pooled.Trinity.fasta
#     DC  →  01_processed/00_1_DC_trinity/DC_trinity_pooled.Trinity.fasta
#     DG  →  01_processed/00_2_DG_trinity/DG_trinity_pooled.Trinity.fasta
#     etc.
#
#   2. CD-HIT deduplicated (optional second argument "cdhit"):
#     MF  →  01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
#     DC  →  01_processed/00_6_cdhit/DC_trinity_cdhit95.fasta
#     etc.
#     Run CD-HIT first: sbatch scripts/04_busco/run_cdhit.sbatch MF
#
#   3. Longest isoform per gene (one script — builds FASTA + runs BUSCO):
#     sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF
#     Output FASTA: 01_processed/00_6_cdhit/<CODE>_longest_isoform.fasta
#
# OUTPUTS (flat — one folder per species, results live directly in it):
#   01_processed/00_5_BUSCO/busco_<CODE>/                       ← pooled Trinity
#   01_processed/00_5_BUSCO/busco_<CODE>_cdhit/                 ← CD-HIT input
#   01_processed/00_5_BUSCO/busco_<CODE>_longest_isoform/        ← longest isoform
#   01_processed/00_5_BUSCO/busco_<CODE>/short_summary*.txt     ← read this first
#   01_processed/00_5_BUSCO/busco_<CODE>/run_<lineage>/         ← per-lineage details
#       └── full_table.tsv                                       (gene-by-gene)
#       └── missing_busco_list.tsv                               (genes not found)
#   01_processed/00_5_BUSCO/busco_<CODE>/logs/busco.log         ← full debug log
#
#   SLURM logs (.out / .err) also go in 01_processed/00_5_BUSCO/.
#
# HOW TO SUBMIT:
#   sbatch scripts/04_busco/run_busco.sbatch MF              # pooled Trinity (default)
#   sbatch scripts/04_busco/run_busco.sbatch MF cdhit        # CD-HIT assembly
#   sbatch scripts/04_busco/run_cdhit.sbatch MF              # run CD-HIT first
#   sbatch scripts/04_busco/run_busco.sbatch DC              # carrot (D. carota)
#   sbatch scripts/04_busco/run_busco.sbatch DG              # carrot (D. glaber)
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
#

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
    echo "ERROR: Could not locate scripts/config.sh."
    echo "Submit this job from inside the rnaseq repo, e.g.:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/04_busco/run_busco.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# READ THE SPECIES CODE FROM THE COMMAND LINE
# =============================================================================
# $1 = species code (required). Valid codes are defined in config.sh:
#      MF, DC, DG, SK, DCDG.
# $2 = input type (optional): "cdhit" for CD-HIT assembly; omit for pooled Trinity.

SPECIES_CODE="${1:-}"
INPUT_TYPE="${2:-}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/04_busco/run_busco.sbatch <SPECIES_CODE> [cdhit]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_busco/run_busco.sbatch MF          # pooled Trinity"
    echo "  sbatch scripts/04_busco/run_busco.sbatch MF cdhit      # CD-HIT assembly"
    echo "  sbatch scripts/04_busco/run_busco.sbatch DC          # carrot (D. carota)"
    echo "  sbatch scripts/04_busco/run_busco.sbatch DG          # carrot (D. glaber)"
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

# Default input = pooled Trinity (ASSEMBLY set in case block above).
# Pass "cdhit" as $2 (or USE_CDHIT=1) to use the CD-HIT deduplicated FASTA.
USE_CDHIT="${USE_CDHIT:-0}"

CDHIT_DIR="${PROCESSED_DIR}/00_6_cdhit"
CDHIT_ASSEMBLY="${CDHIT_DIR}/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"

INPUT_TYPE_LOWER=$(echo "${INPUT_TYPE}" | tr '[:upper:]' '[:lower:]')

if [[ "${INPUT_TYPE_LOWER}" == "cdhit" ]] || [[ "${USE_CDHIT}" == "1" ]]; then
    ASSEMBLY_SOURCE="cdhit"
    ASSEMBLY="${CDHIT_ASSEMBLY}"
    BUSCO_RUN_SUFFIX="_cdhit"
else
    ASSEMBLY_SOURCE="pooled_trinity"
    BUSCO_RUN_SUFFIX=""
fi


# =============================================================================
# REMAINING PATHS — same for all species
# =============================================================================

# Top-level folder where ALL BUSCO results for ALL species live.
# Sits next to the Trinity output folders under 01_processed/.
BUSCO_BASE_DIR="${PROCESSED_DIR}/00_5_BUSCO"

# Per-species folder name (just the name, no path — BUSCO will create it).
# Examples: busco_MF (pooled), busco_MF_cdhit (CD-HIT input).
BUSCO_RUN_NAME="busco_${SPECIES_CODE_UPPER}${BUSCO_RUN_SUFFIX}"

# Full path to the per-species output folder. BUSCO writes all of its outputs
# (logs/, run_embryophyta_odb10/, short_summary*.txt, etc.) directly into this
# folder — no extra "busco_MF_results/" nesting.
BUSCO_OUT_DIR="${BUSCO_BASE_DIR}/${BUSCO_RUN_NAME}"

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
echo "Input source:  ${ASSEMBLY_SOURCE}"
echo "Job ID:        ${SLURM_JOB_ID:-N/A}"
echo "Host:          $(hostname)"
echo "Started:       $(date)"
echo "========================================================"


# =============================================================================
# ACTIVATE CONDA
#
# BUSCO uses its own conda env (busco_env) because the main rnaseq env uses
# Python 3.14, but BUSCO/AUGUSTUS/SEPP don't support Python 3.14 yet.
#
# Create busco_env ONCE on a login node BEFORE running this script:
#   bash scripts/04_busco/setup_busco_env.sh
#
# That helper script handles everything: checks for internet, picks
# conda/mamba, creates the env, and verifies all tools work.
#
# If you want to use a different env name, set BUSCO_CONDA_ENV before sbatch:
#   BUSCO_CONDA_ENV=my_other_env sbatch scripts/04_busco/run_busco.sbatch MF
# =============================================================================
# Default to busco_env, but allow override via the BUSCO_CONDA_ENV variable.
BUSCO_CONDA_ENV="${BUSCO_CONDA_ENV:-busco_env}"

source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${BUSCO_CONDA_ENV}"

# Confirm BUSCO is available in the activated env
if ! command -v busco &> /dev/null; then
    echo ""
    echo "ERROR: busco not found in conda env '${BUSCO_CONDA_ENV}'."
    echo ""
    echo "Run the one-time setup helper on a login node, then resubmit:"
    echo "  bash scripts/04_busco/setup_busco_env.sh"
    echo ""
    echo "Why a separate env? The main rnaseq env uses Python 3.14, but BUSCO"
    echo "and its dependencies (SEPP, AUGUSTUS) don't support Python 3.14 yet."
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
    echo "ERROR: One or more required BUSCO dependencies are missing from"
    echo "       conda env '${BUSCO_CONDA_ENV}'."
    echo ""
    echo "Easiest fix — run the one-time setup helper on a login node:"
    echo "  bash scripts/04_busco/setup_busco_env.sh"
    echo ""
    echo "(That helper rebuilds the env from scratch with all 5 tools."
    echo " Then resubmit this BUSCO job.)"
    exit 1
fi

echo "All required dependencies found."
echo "------------------------------------------------------------"


# =============================================================================
# CHECK THAT THE INPUT ASSEMBLY EXISTS
# =============================================================================
echo ""
echo "Checking input file..."
if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found at:"
    echo "  ${ASSEMBLY}"
    echo ""
    if [[ "${ASSEMBLY_SOURCE}" == "cdhit" ]]; then
        echo "Run CD-HIT on the pooled Trinity assembly first:"
        echo "  sbatch scripts/04_busco/run_cdhit.sbatch ${SPECIES_CODE_UPPER}"
        echo ""
        echo "Or run BUSCO on pooled Trinity instead (omit the cdhit argument):"
        echo "  sbatch scripts/04_busco/run_busco.sbatch ${SPECIES_CODE_UPPER}"
    elif [[ -d "${TRINITY_WORK_DIR}" ]]; then
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
# mkdir -p creates parent folders too, and is silent if the folder exists.
# We create BUSCO_BASE_DIR (the parent), but NOT BUSCO_OUT_DIR — BUSCO will
# create that itself, and BUSCO refuses to start if its output folder already
# exists (unless --force is used, which we do use, but it's cleaner to let
# BUSCO own the folder).
mkdir -p "${BUSCO_BASE_DIR}"
mkdir -p "${BUSCO_LINEAGE_DIR}"


# =============================================================================
# RUN BUSCO
# =============================================================================
echo ""
echo "========================================================"
echo "Running BUSCO..."
echo "========================================================"

# The BUSCO command — how the path arguments combine:
#   --out_path  → parent folder where BUSCO will create its output folder
#   -o          → name of the output folder itself
#
# So with --out_path 00_5_BUSCO and -o busco_MF, BUSCO writes to:
#   00_5_BUSCO/busco_MF/
#       ├── logs/
#       ├── run_<lineage>/
#       └── short_summary*.txt
#
# Other flags:
#   -i              → input FASTA (your Trinity assembly)
#   -l              → lineage database (chosen above based on species)
#   -m              → mode: "transcriptome" (not genome or protein)
#   -c              → number of CPUs
#   --download_path → where to save/find the lineage database
#
# NOTE: The first time you run this for a given lineage, BUSCO downloads
# that lineage (~700 MB) to BUSCO_LINEAGE_DIR. This can take 5–10 minutes.
# Subsequent runs skip the download and start immediately.

BUSCO_CMD="busco \
    -i ${ASSEMBLY} \
    -l ${LINEAGE} \
    -o ${BUSCO_RUN_NAME} \
    --out_path ${BUSCO_BASE_DIR} \
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
    echo "  Check: ${BUSCO_OUT_DIR}/logs/busco.log"
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

SUMMARY_FILE=$(find "${BUSCO_OUT_DIR}" \
    -name "short_summary*.txt" 2>/dev/null | head -1)

if [[ -f "${SUMMARY_FILE}" ]]; then
    cat "${SUMMARY_FILE}"
else
    echo "Summary file not found — check the output folder:"
    echo "  ${BUSCO_OUT_DIR}/"
fi


# =============================================================================
# DONE
# =============================================================================
# Find the per-lineage results subfolder (e.g., run_embryophyta_odb10/) so we
# can show the full_table.tsv path. BUSCO names this folder after the lineage,
# which can change between odb10/odb12/etc., so we glob for it.
LINEAGE_DIR=$(find "${BUSCO_OUT_DIR}" -maxdepth 1 -type d -name "run_*" 2>/dev/null | head -1)

echo ""
echo "========================================================"
echo "✓ COMPLETE"
echo "========================================================"
echo "Species:         ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Lineage used:    ${LINEAGE}"
echo "Results folder:  ${BUSCO_OUT_DIR}/"
echo "Summary:         ${SUMMARY_FILE:-see results folder}"
echo "Full table:      ${LINEAGE_DIR:-${BUSCO_OUT_DIR}/run_<lineage>}/full_table.tsv"
echo ""
echo "What to look for:"
echo "  C > 70% (ideally 85%+) = good assembly"
if [[ "${ASSEMBLY_SOURCE}" == "cdhit" ]]; then
    echo "  D% should be LOWER than pooled Trinity (CD-HIT collapsed isoforms)"
else
    echo "  High D% is NORMAL for pooled Trinity (isoforms inflate duplicates)"
fi
echo "  M > 30% = possible assembly gaps worth investigating"
echo ""
echo "Next step: Kraken2 decontamination on Trinity.fasta"
echo "Finished: $(date)"
echo "========================================================"

exit 0
