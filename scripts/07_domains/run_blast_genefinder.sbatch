#!/bin/bash
# =============================================================================
# SCRIPT:  run_blast_genefinder.sbatch  (CYP/OMT pipeline Script 03)
#
# WHAT THIS DOES:
#   Takes a few known CYP450 and OMT PROTEIN sequences and searches them
#   against your Trinity assembly with tblastn (protein query vs translated
#   nucleotide database).
#
#   This is a second opinion / safety net alongside HMMER (Script 02):
#   the final candidate list should be the UNION of HMMER + BLAST hits so
#   you don't miss family members.
#
#   Unlike HMMER, this step does NOT need TransDecoder — only the assembly
#   FASTA. You can run it as soon as CD-HIT finishes (parallel with RSEM or
#   TransDecoder).
#
# WHEN TO RUN (Step 10 in docs/README_trinity.md):
#   AFTER CD-HIT (Step 3) — does not wait for TransDecoder:
#     sbatch scripts/04_busco/run_cdhit.sbatch MF
#
# HOW TO SUBMIT (MF nutmeg example):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/07_domains/run_blast_genefinder.sbatch MF
#
# TO RE-RUN FROM SCRATCH:
#   sbatch scripts/07_domains/run_blast_genefinder.sbatch MF --force
#
# OPTIONAL FILTERS (off by default — keep all hits at e-value 1e-5):
#   MIN_PIDENT=30 MIN_LENGTH=100 sbatch scripts/07_domains/run_blast_genefinder.sbatch MF
#
# INPUTS REQUIRED:
#   1. CD-HIT Trinity assembly (same transcript IDs as RSEM / PyDESeq2)
#      MF: 01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
#   2. Known query proteins (YOU must populate these from UniProt/NCBI):
#      reference/queries/known_cyp450.fasta
#      reference/queries/known_omt.fasta
#   3. BLAST+ in rnaseq conda env (makeblastdb, tblastn)
#
# OUTPUT (MF example):
#   06_analysis/blast_genefinder_MF/cyp450_blast_hits.tsv
#   06_analysis/blast_genefinder_MF/omt_blast_hits.tsv
#   06_analysis/blast_genefinder_MF/cyp450_blast_ids.txt
#   06_analysis/blast_genefinder_MF/omt_blast_ids.txt
#   06_analysis/blast_genefinder_MF/trinity_db.*  (local BLAST database)
#
# PATHS come from scripts/config.sh (BASE_DIR, PROCESSED_DIR, etc.)
# =============================================================================


# ---- SLURM: tell the cluster how much time/memory this job needs ----
# tblastn with a handful of queries is FAST (minutes). We request 1 hour to be safe.
#SBATCH --job-name=blast_genefinder
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o blast_genefinder_%j.out
#SBATCH -e blast_genefinder_%j.err

# set -e  = stop if any command fails
# set -u  = stop if we use a variable that was never set
# set -o pipefail = treat a failed command in a pipe as a failure
set -euo pipefail


# =============================================================================
# SECTION 1: LOAD config.sh
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
    echo "ERROR: Could not find scripts/config.sh."
    echo "Submit from the repo folder:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/07_domains/run_blast_genefinder.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# SECTION 2: READ WHAT THE USER TYPED AFTER sbatch
# =============================================================================
# Example: sbatch run_blast_genefinder.sbatch MF --force
#   $1 = MF       (species code — required)
#   --force       (optional — re-run even if output files already exist)

SPECIES_CODE="${1:-}"
FORCE=0

for arg in "$@"; do
    if [[ "${arg}" == "--force" ]]; then
        FORCE=1
    fi
done

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: You must give a species code (e.g. MF for nutmeg)."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/07_domains/run_blast_genefinder.sbatch <CODE> [--force]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/07_domains/run_blast_genefinder.sbatch MF"
    echo "  sbatch scripts/07_domains/run_blast_genefinder.sbatch MF --force"
    echo ""
    echo "Optional filters (default OFF):"
    echo "  MIN_PIDENT=30 MIN_LENGTH=100 sbatch scripts/07_domains/run_blast_genefinder.sbatch MF"
    echo ""
    echo "Required inputs (MF example):"
    echo "  ${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta"
    echo "  ${CODE_DIR}/reference/queries/known_cyp450.fasta"
    echo "  ${CODE_DIR}/reference/queries/known_omt.fasta"
    echo ""
    echo "Run CD-HIT first if the assembly is missing:"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch MF"
    exit 1
fi

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# SECTION 3: ACTIVATE CONDA (BLAST+ lives in rnaseq env)
# =============================================================================

echo ""
echo "Activating conda environment: ${CONDA_ENV}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"


# =============================================================================
# SECTION 4: SET FILE PATHS
# =============================================================================
# Use the CD-HIT assembly — same transcript IDs as RSEM counts and HMMER output.

ASSEMBLY_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"
QUERY_DIR="${CODE_DIR}/reference/queries"
CYP_QUERY="${QUERY_DIR}/known_cyp450.fasta"
OMT_QUERY="${QUERY_DIR}/known_omt.fasta"

OUTPUT_DIR="${ANALYSIS_DIR}/blast_genefinder_${SPECIES_CODE_UPPER}"
BLAST_DB="${OUTPUT_DIR}/trinity_db"

CYP_HITS="${OUTPUT_DIR}/cyp450_blast_hits.tsv"
OMT_HITS="${OUTPUT_DIR}/omt_blast_hits.tsv"
CYP_OUT="${OUTPUT_DIR}/cyp450_blast_ids.txt"
OMT_OUT="${OUTPUT_DIR}/omt_blast_ids.txt"

# Optional filters — empty means OFF (keep all hits at e-value cutoff)
MIN_PIDENT="${MIN_PIDENT:-}"
MIN_LENGTH="${MIN_LENGTH:-}"

# tblastn settings from the brief
EVALUE="1e-5"
MAX_TARGET_SEQS=50


# =============================================================================
# SECTION 5: CHECK THAT BLAST+ IS INSTALLED
# =============================================================================

echo ""
echo "=========================================="
echo "BLAST CYP/OMT gene finder (Script 03)"
echo "=========================================="
echo "Species:     ${SPECIES_NAME} (${SPECIES_CODE_UPPER})"
echo "Assembly:    ${ASSEMBLY_FASTA}"
echo "CYP queries: ${CYP_QUERY}"
echo "OMT queries: ${OMT_QUERY}"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Force:       ${FORCE}"
echo "MIN_PIDENT:  ${MIN_PIDENT:-OFF}"
echo "MIN_LENGTH:  ${MIN_LENGTH:-OFF}"
echo "------------------------------------------"

if ! command -v makeblastdb &> /dev/null; then
    echo "ERROR: makeblastdb not found on PATH."
    echo "Install BLAST+ in your conda env:"
    echo "  conda activate ${CONDA_ENV}"
    echo "  conda install -c bioconda blast"
    exit 1
fi

if ! command -v tblastn &> /dev/null; then
    echo "ERROR: tblastn not found on PATH."
    echo "Install BLAST+ in your conda env:"
    echo "  conda activate ${CONDA_ENV}"
    echo "  conda install -c bioconda blast"
    exit 1
fi

# tblastn uses SINGLE-dash -version (not --version)
BLAST_VERSION=$(tblastn -version 2>&1 | head -1)
echo "BLAST+:      ${BLAST_VERSION}"


# =============================================================================
# SECTION 6: CHECK THAT INPUT FILES EXIST
# =============================================================================

if [[ ! -f "${ASSEMBLY_FASTA}" ]]; then
    echo ""
    echo "ERROR: Trinity assembly not found:"
    echo "  ${ASSEMBLY_FASTA}"
    echo ""
    echo "Run CD-HIT first (Step 3 in docs/README_trinity.md):"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch ${SPECIES_CODE_UPPER}"
    exit 1
fi

if [[ ! -f "${CYP_QUERY}" ]]; then
    echo ""
    echo "ERROR: CYP450 query file not found:"
    echo "  ${CYP_QUERY}"
    echo ""
    echo "Add 3-5 characterized plant CYP450 proteins. See:"
    echo "  ${QUERY_DIR}/README.md"
    exit 1
fi

if [[ ! -f "${OMT_QUERY}" ]]; then
    echo ""
    echo "ERROR: OMT query file not found:"
    echo "  ${OMT_QUERY}"
    echo ""
    echo "Add 3-5 characterized plant OMT proteins. See:"
    echo "  ${QUERY_DIR}/README.md"
    exit 1
fi

# Warn if query files still contain PLACEHOLDER headers
if grep -q '^>PLACEHOLDER' "${CYP_QUERY}"; then
    echo ""
    echo "WARNING: ${CYP_QUERY} still has PLACEHOLDER sequences."
    echo "Replace them with real CYP450 proteins from UniProt/NCBI before trusting results."
    echo "See: ${QUERY_DIR}/README.md"
fi

if grep -q '^>PLACEHOLDER' "${OMT_QUERY}"; then
    echo ""
    echo "WARNING: ${OMT_QUERY} still has PLACEHOLDER sequences."
    echo "Replace them with real OMT proteins from UniProt/NCBI before trusting results."
    echo "See: ${QUERY_DIR}/README.md"
fi

N_TRANSCRIPTS=$(grep -c '^>' "${ASSEMBLY_FASTA}" || true)
N_CYP_QUERIES=$(grep -c '^>' "${CYP_QUERY}" || true)
N_OMT_QUERIES=$(grep -c '^>' "${OMT_QUERY}" || true)

echo "Assembly:    ${N_TRANSCRIPTS} transcript sequences"
echo "CYP queries: ${N_CYP_QUERIES} protein sequences"
echo "OMT queries: ${N_OMT_QUERIES} protein sequences"

mkdir -p "${OUTPUT_DIR}"


# =============================================================================
# SECTION 7: SKIP IF OUTPUTS ALREADY EXIST (--force to overwrite)
# =============================================================================

if [[ "${FORCE}" -eq 0 ]] && [[ -f "${CYP_OUT}" ]] && [[ -f "${OMT_OUT}" ]]; then
    echo ""
    echo "Outputs already exist — skipping (use --force to re-run):"
    echo "  ${CYP_OUT}"
    echo "  ${OMT_OUT}"
    echo ""
    echo "Candidate counts:"
    echo "  CYP450: $(wc -l < "${CYP_OUT}" | tr -d ' ') transcripts"
    echo "  OMT:    $(wc -l < "${OMT_OUT}" | tr -d ' ') transcripts"
    exit 0
fi

if [[ "${FORCE}" -eq 1 ]]; then
    echo ""
    echo "--force set: rebuilding BLAST database and overwriting outputs."
fi


# =============================================================================
# SECTION 8: BUILD LOCAL BLAST DATABASE FROM ASSEMBLY
# =============================================================================
# makeblastdb indexes the nucleotide FASTA so tblastn can search it quickly.
# -dbtype nucl = nucleotide database (Trinity transcripts are DNA/RNA)

if [[ "${FORCE}" -eq 1 ]] || [[ ! -f "${BLAST_DB}.nin" ]]; then
    echo ""
    echo "Building BLAST database (makeblastdb)..."
    makeblastdb \
        -in "${ASSEMBLY_FASTA}" \
        -dbtype nucl \
        -out "${BLAST_DB}"
else
    echo ""
    echo "BLAST database already exists — skipping makeblastdb:"
    echo "  ${BLAST_DB}"
    echo "  (use --force to rebuild)"
fi


# =============================================================================
# SECTION 9: RUN tblastn FOR EACH FAMILY
# =============================================================================
# tblastn = protein query searched against a translated nucleotide database.
# outfmt 6 = tab-separated table (easy to parse).
# Column 2 (sseqid) = matched Trinity transcript ID.

run_tblastn() {
    # $1 = label (cyp450 or omt)
    # $2 = query FASTA path
    # $3 = output TSV path
    local label="$1"
    local query="$2"
    local outfile="$3"

    echo ""
    echo "  [${label}] tblastn"
    echo "    Query:  ${query}"
    echo "    DB:     ${BLAST_DB}"
    echo "    Output: ${outfile}"

    tblastn \
        -query "${query}" \
        -db "${BLAST_DB}" \
        -evalue "${EVALUE}" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        -max_target_seqs "${MAX_TARGET_SEQS}" \
        -num_threads "${SLURM_CPUS_PER_TASK:-4}" \
        > "${outfile}"

    echo "    Done. $(wc -l < "${outfile}" | tr -d ' ') alignment rows written."
}

echo ""
echo "Running tblastn (usually a few minutes)..."

run_tblastn "cyp450" "${CYP_QUERY}" "${CYP_HITS}"
run_tblastn "omt"    "${OMT_QUERY}" "${OMT_HITS}"


# =============================================================================
# SECTION 10: PARSE TSV FILES → UNIQUE TRANSCRIPT ID LISTS
# =============================================================================
# Column 2 = sseqid (Trinity transcript ID).
# Optional MIN_PIDENT / MIN_LENGTH filters apply here if you set them.

parse_blast_to_ids() {
    # $1 = input TSV
    # $2 = output ID list (one per line)
    local infile="$1"
    local outfile="$2"

    if [[ ! -s "${infile}" ]]; then
        : > "${outfile}"
        return 0
    fi

    if [[ -n "${MIN_PIDENT}" ]] && [[ -n "${MIN_LENGTH}" ]]; then
        awk -v min_pident="${MIN_PIDENT}" -v min_len="${MIN_LENGTH}" \
            '$3 >= min_pident && $4 >= min_len {print $2}' "${infile}" \
            | sort -u > "${outfile}"
    elif [[ -n "${MIN_PIDENT}" ]]; then
        awk -v min_pident="${MIN_PIDENT}" \
            '$3 >= min_pident {print $2}' "${infile}" \
            | sort -u > "${outfile}"
    elif [[ -n "${MIN_LENGTH}" ]]; then
        awk -v min_len="${MIN_LENGTH}" \
            '$4 >= min_len {print $2}' "${infile}" \
            | sort -u > "${outfile}"
    else
        awk '{print $2}' "${infile}" | sort -u > "${outfile}"
    fi
}

echo ""
echo "Parsing BLAST tables → transcript ID lists..."

parse_blast_to_ids "${CYP_HITS}" "${CYP_OUT}"
parse_blast_to_ids "${OMT_HITS}" "${OMT_OUT}"


# =============================================================================
# SECTION 11: PRINT SUMMARY
# =============================================================================

N_CYP=$(wc -l < "${CYP_OUT}" | tr -d ' ')
N_OMT=$(wc -l < "${OMT_OUT}" | tr -d ' ')

echo ""
echo "=========================================="
echo "BLAST gene finder complete"
echo "=========================================="
echo "CYP450 candidates: ${N_CYP} transcript IDs"
echo "  → ${CYP_OUT}"
echo ""
echo "OMT candidates:    ${N_OMT} transcript IDs"
echo "  → ${OMT_OUT}"
echo ""
echo "Raw BLAST tables (keep for e-value filtering in Script 04):"
echo "  ${CYP_HITS}"
echo "  ${OMT_HITS}"
echo ""
echo "Next: combine with HMMER hits (union = final candidate set):"
echo "  sort -u ${ANALYSIS_DIR}/hmmer_genefinder_${SPECIES_CODE_UPPER}/cyp450_hmmer_ids.txt ${CYP_OUT} > combined_cyp450_ids.txt"
echo "  sort -u ${ANALYSIS_DIR}/hmmer_genefinder_${SPECIES_CODE_UPPER}/omt_hmmer_ids.txt ${OMT_OUT} > combined_omt_ids.txt"
echo "=========================================="
