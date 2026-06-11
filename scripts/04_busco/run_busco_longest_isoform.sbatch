#!/bin/bash
# =============================================================================
# SCRIPT:  run_busco_longest_isoform.sbatch
# PURPOSE: Two-step job in one submission:
#          Step 1 — build longest-isoform-per-gene FASTA from pooled Trinity
#          Step 2 — run BUSCO on that assembly (third D% benchmark)
#
# WHY THIS BUSCO RUN MATTERS:
#   Compare THREE assembly versions:
#     1. Raw pooled Trinity     → high D (isoforms + gene families)
#     2. CD-HIT 95%             → lower D (near-duplicates removed)
#     3. Longest isoform        → lowest D from isoforms (gene families remain)
#
# INPUT (Step 1):
#   Pooled Trinity FASTA + gene_trans_map
#
# OUTPUT:
#   Step 1: 01_processed/00_6_cdhit/<CODE>_longest_isoform.fasta
#   Step 2: 01_processed/00_5_BUSCO/busco_<CODE>_longest_isoform/
#
# HOW TO SUBMIT:
#   sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF
#
# NOTE: Independent of RSEM — both can run in parallel.
#       Use CD-HIT assembly for RSEM, not longest isoform.
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


#SBATCH --job-name=busco_longest
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_5_BUSCO
#SBATCH -o busco_longest_isoform_%j.out
#SBATCH -e busco_longest_isoform_%j.err

set -euo pipefail


# =============================================================================
# LOAD CONFIG
# =============================================================================
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
    echo "Submit from inside the rnaseq repo, e.g.:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# READ SPECIES CODE
# =============================================================================
SPECIES_CODE="${1:-}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/04_busco/run_busco_longest_isoform.sbatch <SPECIES_CODE>"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_busco/run_busco_longest_isoform.sbatch MF"
    exit 1
fi

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# PICK LINEAGE + TRINITY PATHS (same mapping as run_busco.sbatch)
# =============================================================================
case "${SPECIES_CODE_UPPER}" in
    MF)
        LINEAGE="embryophyta_odb10"
        LINEAGE_NOTE="land plants (1,614 genes) — basal angiosperm nutmeg"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_MF_trinity"
        POOLED_FASTA="${TRINITY_DIR}/MF_trinity_pooled.Trinity.fasta"
        GENE_MAP="${TRINITY_DIR}/MF_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DC)
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes)"
        TRINITY_DIR="${PROCESSED_DIR}/00_1_DC_trinity"
        POOLED_FASTA="${TRINITY_DIR}/DC_trinity_pooled.Trinity.fasta"
        GENE_MAP="${TRINITY_DIR}/DC_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DG)
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes)"
        TRINITY_DIR="${PROCESSED_DIR}/00_2_DG_trinity"
        POOLED_FASTA="${TRINITY_DIR}/DG_trinity_pooled.Trinity.fasta"
        GENE_MAP="${TRINITY_DIR}/DG_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    SK)
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes)"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_Sukman_trinity"
        POOLED_FASTA="${TRINITY_DIR}/SK_trinity_pooled.Trinity.fasta"
        GENE_MAP="${TRINITY_DIR}/SK_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DCDG)
        LINEAGE="eudicots_odb10"
        LINEAGE_NOTE="eudicots (~2,326 genes)"
        TRINITY_DIR="${PROCESSED_DIR}/00_4_DC_DG_trinity"
        POOLED_FASTA="${TRINITY_DIR}/DCDG_trinity_pooled.Trinity.fasta"
        GENE_MAP="${TRINITY_DIR}/DCDG_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    *)
        echo "ERROR: No paths configured for '${SPECIES_CODE_UPPER}'."
        echo "Add a case block to run_busco_longest_isoform.sbatch."
        exit 1
        ;;
esac


# =============================================================================
# PATHS
# =============================================================================
CDHIT_DIR="${PROCESSED_DIR}/00_6_cdhit"
INPUT="${CDHIT_DIR}/${SPECIES_CODE_UPPER}_longest_isoform.fasta"
BUSCO_BASE_DIR="${PROCESSED_DIR}/00_5_BUSCO"
BUSCO_RUN_NAME="busco_${SPECIES_CODE_UPPER}_longest_isoform"
BUSCO_OUT_DIR="${BUSCO_BASE_DIR}/${BUSCO_RUN_NAME}"
BUSCO_LINEAGE_DIR="${BASE_DIR}/04_reference/busco_lineages"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
FORCE_REBUILD="${FORCE_REBUILD:-0}"

mkdir -p "${CDHIT_DIR}"
mkdir -p "${BUSCO_BASE_DIR}"
mkdir -p "${BUSCO_LINEAGE_DIR}"


# =============================================================================
# FIND PYTHON HELPER
# =============================================================================
SCRIPT_DIR=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/04_busco" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/04_busco"
do
    if [[ -f "${candidate}/extract_longest_isoform.py" ]]; then
        SCRIPT_DIR="${candidate}"
        break
    fi
done

if [[ -z "${SCRIPT_DIR}" ]]; then
    echo "ERROR: extract_longest_isoform.py not found."
    exit 1
fi


# =============================================================================
# JOB HEADER
# =============================================================================
echo "========================================================"
echo "Longest Isoform + BUSCO — Two-Step Job"
echo "Species code:  ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Lineage:       ${LINEAGE}"
echo "               (${LINEAGE_NOTE})"
echo "Job ID:        ${SLURM_JOB_ID:-N/A}"
echo "Host:          $(hostname)"
echo "Started:       $(date)"
echo "========================================================"


# =============================================================================
# STEP 1: BUILD LONGEST-ISOFORM FASTA
# =============================================================================
echo ""
echo "========================================================"
echo "STEP 1: Extract longest isoform per gene"
echo "========================================================"

if [[ ! -f "${POOLED_FASTA}" ]]; then
    echo "ERROR: Pooled Trinity FASTA not found: ${POOLED_FASTA}"
    echo "Run trinity_pooled.sbatch for ${SPECIES_CODE_UPPER} first."
    exit 1
fi

if [[ ! -f "${GENE_MAP}" ]]; then
    echo "ERROR: Gene map not found: ${GENE_MAP}"
    exit 1
fi

POOLED_COUNT=$(grep -c ">" "${POOLED_FASTA}" || echo "0")
echo "Pooled Trinity: ${POOLED_COUNT} transcripts"
echo "  ${POOLED_FASTA}"
echo "Gene map: $(wc -l < "${GENE_MAP}") lines"
echo "  ${GENE_MAP}"
echo "------------------------------------------------------------"

if [[ -f "${INPUT}" ]] && [[ "${FORCE_REBUILD}" != "1" ]]; then
    INPUT_COUNT=$(grep -c ">" "${INPUT}" || echo "0")
    echo "✓ Longest isoform FASTA already exists — skipping extraction"
    echo "  ${INPUT} (${INPUT_COUNT} transcripts)"
    echo "  Set FORCE_REBUILD=1 to regenerate."
else
    echo "Running extract_longest_isoform.py..."
    python3 "${SCRIPT_DIR}/extract_longest_isoform.py" \
        --fasta "${POOLED_FASTA}" \
        --gene-map "${GENE_MAP}" \
        --output "${INPUT}"

    INPUT_COUNT=$(grep -c ">" "${INPUT}" || echo "0")
    echo ""
    echo "✓ Step 1 complete: ${POOLED_COUNT} → ${INPUT_COUNT} transcripts"
    echo "  Output: ${INPUT}"
fi
echo "------------------------------------------------------------"


# =============================================================================
# THREE-WAY COMPARISON TABLE (before BUSCO)
# =============================================================================
echo ""
echo "Three-way BUSCO comparison:"
echo "─────────────────────────────────────────────────────────"
echo "1. Pooled Trinity         | busco_${SPECIES_CODE_UPPER}/"
echo "2. CD-HIT 95%             | busco_${SPECIES_CODE_UPPER}_cdhit/"
echo "3. Longest isoform        | busco_${SPECIES_CODE_UPPER}_longest_isoform/  ← THIS RUN"
echo "─────────────────────────────────────────────────────────"

if [[ "${SPECIES_CODE_UPPER}" == "MF" ]]; then
    echo ""
    echo "MF reference numbers from prior runs:"
    echo "  Raw Trinity:     547,893 transcripts | C:95.5% | D:78.1%"
    echo "  CD-HIT 95%:      467,314 transcripts | C:95.4% | D:65.5%"
    echo "  Longest isoform: ${INPUT_COUNT} transcripts | C:? | D:?"
fi
echo "------------------------------------------------------------"


# =============================================================================
# STEP 2: ACTIVATE BUSCO ENV + RUN BUSCO
# =============================================================================
echo ""
echo "========================================================"
echo "STEP 2: Run BUSCO on longest isoform assembly"
echo "========================================================"

BUSCO_CONDA_ENV="${BUSCO_CONDA_ENV:-busco_env}"

echo "Activating BUSCO conda environment..."
source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${BUSCO_CONDA_ENV}"

if ! command -v busco &>/dev/null; then
    echo "ERROR: busco not found in conda env '${BUSCO_CONDA_ENV}'"
    echo "Set up with: bash scripts/04_busco/setup_busco_env.sh"
    exit 1
fi

echo "✓ BUSCO version: $(busco --version 2>&1)"
echo ""
echo "Input:    ${INPUT} (${INPUT_COUNT} transcripts)"
echo "Lineage:  ${LINEAGE}"
echo "Output:   ${BUSCO_OUT_DIR}/"
echo "Mode:     transcriptome"
echo "CPUs:     ${THREADS}"
echo "------------------------------------------------------------"

busco \
    -i "${INPUT}" \
    -l "${LINEAGE}" \
    -o "${BUSCO_RUN_NAME}" \
    --out_path "${BUSCO_BASE_DIR}" \
    -m transcriptome \
    -c "${THREADS}" \
    --download_path "${BUSCO_LINEAGE_DIR}" \
    --force

echo ""
echo "✓ Step 2 complete — BUSCO finished successfully"


# =============================================================================
# DISPLAY RESULTS + THREE-WAY COMPARISON
# =============================================================================
echo ""
echo "========================================================"
echo "BUSCO Results"
echo "========================================================"

SUMMARY_FILE=$(find "${BUSCO_OUT_DIR}" \
    -name "short_summary*.txt" 2>/dev/null | head -1)

if [[ -f "${SUMMARY_FILE}" ]]; then
    cat "${SUMMARY_FILE}"
else
    echo "Summary file not found — check: ${BUSCO_OUT_DIR}/"
fi

echo ""
echo "========================================================"
echo "Three-Way BUSCO Comparison"
echo "========================================================"

POOL_SUMMARY=$(find "${BUSCO_BASE_DIR}/busco_${SPECIES_CODE_UPPER}" \
    -maxdepth 1 -name "short_summary*.txt" 2>/dev/null | head -1)
CDHIT_SUMMARY=$(find "${BUSCO_BASE_DIR}/busco_${SPECIES_CODE_UPPER}_cdhit" \
    -maxdepth 1 -name "short_summary*.txt" 2>/dev/null | head -1)

for LABEL in "Pooled Trinity:${POOL_SUMMARY}" "CD-HIT:${CDHIT_SUMMARY}" "Longest isoform:${SUMMARY_FILE}"; do
    NAME="${LABEL%%:*}"
    FILE="${LABEL#*:}"
    if [[ -f "${FILE}" ]]; then
        echo "${NAME}:"
        grep "^C:" "${FILE}" | head -1
    else
        echo "${NAME}: summary not found (run that BUSCO job first)"
    fi
    echo ""
done

echo "HOW TO INTERPRET:"
echo "  D dropped sharply vs CD-HIT  → isoforms drove high D"
echo "  D dropped moderately         → gene family expansion is real biology"
echo "  C below ~88%                 → check extraction log; use CD-HIT for RSEM"


# =============================================================================
# PHYLOGENOMICS — single-copy BUSCO sequences
# =============================================================================
LINEAGE_DIR=$(find "${BUSCO_OUT_DIR}" \
    -maxdepth 1 -type d -name "run_*" 2>/dev/null | head -1)
SINGLE_COPY_DIR="${LINEAGE_DIR}/busco_sequences/single_copy_busco_sequences"

echo ""
echo "========================================================"
echo "Phylogenomics Resource"
echo "========================================================"

if [[ -d "${SINGLE_COPY_DIR}" ]]; then
    FAA_COUNT=$(find "${SINGLE_COPY_DIR}" -name "*.faa" 2>/dev/null | wc -l)
    FNA_COUNT=$(find "${SINGLE_COPY_DIR}" -name "*.fna" 2>/dev/null | wc -l)
    echo "Single-copy BUSCO sequences:"
    echo "  Protein (.faa):    ${FAA_COUNT} genes"
    echo "  Nucleotide (.fna): ${FNA_COUNT} genes"
    echo "  Location: ${SINGLE_COPY_DIR}/"
else
    echo "Single-copy sequences folder not found."
    echo "Check: ${SINGLE_COPY_DIR}/"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ COMPLETE (Step 1 + Step 2)"
echo "========================================================"
echo "Species:         ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Longest isoform: ${INPUT} (${INPUT_COUNT} transcripts)"
echo "BUSCO results:   ${BUSCO_OUT_DIR}/"
echo "Summary:         ${SUMMARY_FILE:-see results folder}"
echo ""
echo "Next: compare D% vs busco_${SPECIES_CODE_UPPER}_cdhit/"
echo "RSEM uses CD-HIT assembly:"
echo "  sbatch scripts/04_rsem/run_rsem.sbatch ${SPECIES_CODE_UPPER}"
echo "Finished: $(date)"
echo "========================================================"

exit 0
