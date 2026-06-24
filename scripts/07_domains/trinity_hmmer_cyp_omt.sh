#!/bin/bash
# =============================================================================
# SCRIPT:  trinity_hmmer_cyp_omt.sh  (CYP/OMT pipeline — Script 02 of 4)
#
# WHAT THIS DOES:
#   Searches TransDecoder protein sequences for conserved Pfam domains that
#   DEFINE cytochrome P450 (CYP) and O-methyltransferase (OMT) families.
#
#   HMMER profiles are trained on hundreds of known family members, so they
#   catch even divergent nutmeg genes better than keyword search alone.
#   This is the main, most sensitive first-pass net for finding candidates.
#
# INPUTS (from Script 01 — TransDecoder):
#   MF_trinity_cdhit95.fasta.transdecoder.pep
#
# PFAM PROFILES (downloaded automatically if missing):
#   PF00067  — cytochrome P450 domain (CYP450)
#   PF00891  — Methyltransf_2 (OMT)
#   PF08100  — dimerisation domain (OMT)
#   PF01596  — CCoA-type OMT
#
#   NOTE: PF00067 is well verified. Confirm OMT accessions on InterPro for
#   your specific OMT subtype before fully trusting them:
#   https://www.ebi.ac.uk/interpro/
#
# OUTPUTS (in 06_analysis/hmmer_genefinder_<SPECIES>/):
#   cyp450_hmmer_ids.txt   — one Trinity transcript ID per line
#   omt_hmmer_ids.txt      — one Trinity transcript ID per line (3 OMT domains merged)
#   *.domtbl               — raw HMMER domain tables (e-values for Script 04)
#
# HOW TO RUN (login node or inside sbatch):
#   bash scripts/07_domains/trinity_hmmer_cyp_omt.sh MF
#   bash scripts/07_domains/trinity_hmmer_cyp_omt.sh MF --force
#
# ON HPC (recommended):
#   sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF
#
# INPUTS REQUIRED:
#   - TransDecoder .pep file (see run_transdecoder.sbatch)
#   - HMMER in rnaseq conda env (hmmsearch, hmmpress)
#   - Pfam HMM profiles (auto-downloaded to 07_NRdatabase/hmmerdb/pfam_profiles/)
# =============================================================================

set -euo pipefail


# =============================================================================
# SECTION 1: PARSE ARGUMENTS
# =============================================================================

SPECIES_CODE=""
FORCE=0

for arg in "$@"; do
    case "${arg}" in
        --force)
            FORCE=1
            ;;
        -h|--help)
            echo "Usage: trinity_hmmer_cyp_omt.sh <SPECIES_CODE> [--force]"
            echo ""
            echo "Examples:"
            echo "  bash scripts/07_domains/trinity_hmmer_cyp_omt.sh MF"
            echo "  bash scripts/07_domains/trinity_hmmer_cyp_omt.sh MF --force"
            exit 0
            ;;
        *)
            if [[ -z "${SPECIES_CODE}" ]]; then
                SPECIES_CODE="${arg}"
            else
                echo "ERROR: Unknown argument '${arg}'"
                exit 1
            fi
            ;;
    esac
done

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: Species code required (e.g. MF for nutmeg)."
    echo "Usage: trinity_hmmer_cyp_omt.sh <SPECIES_CODE> [--force]"
    exit 1
fi


# =============================================================================
# SECTION 2: LOAD config.sh
# =============================================================================

CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "$(dirname "$0")/../config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not find scripts/config.sh"
    exit 1
fi

source "${CONFIG_FILE}"
get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# SECTION 3: SET PATHS
# =============================================================================

# TransDecoder writes .pep next to the CD-HIT FASTA (same IDs as RSEM / PyDESeq2)
PEP_FILE="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta.transdecoder.pep"

OUTPUT_DIR="${ANALYSIS_DIR}/hmmer_genefinder_${SPECIES_CODE_UPPER}"
HMM_DIR="${DB_DIR}/hmmerdb/pfam_profiles"

CYP_OUT="${OUTPUT_DIR}/cyp450_hmmer_ids.txt"
OMT_OUT="${OUTPUT_DIR}/omt_hmmer_ids.txt"

# Pfam profiles used in this step
declare -A HMM_FILES=(
    [cyp450]="PF00067.hmm"
    [omt_a]="PF00891.hmm"
    [omt_b]="PF08100.hmm"
    [omt_ccoa]="PF01596.hmm"
)


# =============================================================================
# SECTION 4: CHECK DEPENDENCIES
# =============================================================================

echo ""
echo "=========================================="
echo "Trinity CYP/OMT HMMER scan (Script 02)"
echo "=========================================="
echo "Species:     ${SPECIES_NAME} (${SPECIES_CODE_UPPER})"
echo "Proteins:    ${PEP_FILE}"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Force:       ${FORCE}"
echo "------------------------------------------"

if ! command -v hmmsearch &> /dev/null; then
    echo "ERROR: hmmsearch not found on PATH."
    echo "Install HMMER:"
    echo "  conda activate ${CONDA_ENV}"
    echo "  conda install -c bioconda hmmer"
    exit 1
fi

HMMER_VERSION=$(hmmsearch -h 2>&1 | grep -m1 "^# HMMER" | sed 's/^# //')
echo "HMMER:       ${HMMER_VERSION}"

if ! command -v hmmpress &> /dev/null; then
    echo "ERROR: hmmpress not found on PATH (needed to index HMM files)."
    exit 1
fi


# =============================================================================
# SECTION 5: CHECK INPUTS
# =============================================================================

if [[ ! -f "${PEP_FILE}" ]]; then
    echo ""
    echo "ERROR: TransDecoder protein file not found:"
    echo "  ${PEP_FILE}"
    echo ""
    echo "Run TransDecoder first (Script 01):"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch ${SPECIES_CODE_UPPER}"
    exit 1
fi

N_PROTEINS=$(grep -c '^>' "${PEP_FILE}" || true)
echo "Proteins:    ${N_PROTEINS} sequences in .pep file"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${HMM_DIR}"


# =============================================================================
# SECTION 6: SKIP IF OUTPUTS ALREADY EXIST (--force to overwrite)
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
    echo "--force set: re-running hmmsearch and overwriting outputs."
fi


# =============================================================================
# SECTION 7: DOWNLOAD PFAM HMM PROFILES IF MISSING
# =============================================================================

download_hmm() {
    # $1 = Pfam accession (e.g. PF00067)
    local pfam_acc="$1"
    local hmm_path="${HMM_DIR}/${pfam_acc}.hmm"

    if [[ -f "${hmm_path}" ]] && [[ "$(wc -c < "${hmm_path}" | tr -d ' ')" -gt 100 ]]; then
        echo "  Using cached: ${hmm_path}"
        return 0
    fi

    echo "  Downloading ${pfam_acc} from InterPro/Pfam..."

    local urls=(
        "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pfam_acc}?annotation=hmm"
        "https://pfam.xfam.org/family/${pfam_acc}/hmm"
    )

    local url
    for url in "${urls[@]}"; do
        echo "    Trying: ${url}"
        if command -v wget &> /dev/null; then
            wget -q -O "${hmm_path}" "${url}" || true
        elif command -v curl &> /dev/null; then
            curl -sL -o "${hmm_path}" "${url}" || true
        else
            echo "ERROR: Need wget or curl to download HMM profiles."
            exit 1
        fi

        if [[ -f "${hmm_path}" ]] && [[ "$(wc -c < "${hmm_path}" | tr -d ' ')" -gt 100 ]]; then
            # InterPro API may return gzip
            if file "${hmm_path}" | grep -q gzip; then
                gunzip -c "${hmm_path}" > "${hmm_path}.tmp"
                mv "${hmm_path}.tmp" "${hmm_path}"
            fi
            if head -1 "${hmm_path}" | grep -q '^HMMER'; then
                echo "    Downloaded: ${hmm_path}"
                return 0
            fi
        fi
        rm -f "${hmm_path}"
    done

    echo ""
    echo "ERROR: Could not download ${pfam_acc}.hmm automatically."
    echo "Download manually from InterPro and save to:"
    echo "  ${hmm_path}"
    echo "  https://www.ebi.ac.uk/interpro/entry/pfam/${pfam_acc}/"
    exit 1
}

echo ""
echo "Checking Pfam HMM profiles..."

for key in "${!HMM_FILES[@]}"; do
    pfam_acc="${HMM_FILES[$key]%.hmm}"
    download_hmm "${pfam_acc}"
done


# =============================================================================
# SECTION 8: PREPARE HMM FILES (hmmpress)
# =============================================================================

echo ""
echo "Preparing HMM indexes (hmmpress)..."

for key in "${!HMM_FILES[@]}"; do
    hmm_file="${HMM_DIR}/${HMM_FILES[$key]}"
    if [[ ! -f "${hmm_file}.h3m" ]]; then
        echo "  hmmpress ${hmm_file}"
        hmmpress -f "${hmm_file}"
    else
        echo "  Already indexed: ${hmm_file}"
    fi
done


# =============================================================================
# SECTION 9: RUN hmmsearch (--cut_ga = Pfam curated gathering threshold)
# =============================================================================
# --cut_ga is stricter and more biologically meaningful than a raw E-value.
# It uses the threshold Pfam curators set for each family.

echo ""
echo "Running hmmsearch (this usually takes a few minutes)..."

run_hmmsearch() {
    # $1 = label (e.g. cyp450)
    # $2 = HMM filename (e.g. PF00067.hmm)
    local label="$1"
    local hmm_file="${HMM_DIR}/$2"
    local domtbl="${OUTPUT_DIR}/${label}.domtbl"

    echo ""
    echo "  [${label}] hmmsearch --cut_ga"
    echo "    HMM:      ${hmm_file}"
    echo "    Proteins: ${PEP_FILE}"
    echo "    Output:   ${domtbl}"

    hmmsearch \
        --cut_ga \
        --domtblout "${domtbl}" \
        --cpu "${SLURM_CPUS_PER_TASK:-4}" \
        --noali \
        "${hmm_file}" \
        "${PEP_FILE}" \
        > "${OUTPUT_DIR}/${label}.out"

    echo "    Done."
}

run_hmmsearch "cyp450"   "PF00067.hmm"
run_hmmsearch "omt_a"    "PF00891.hmm"
run_hmmsearch "omt_b"    "PF08100.hmm"
run_hmmsearch "omt_ccoa" "PF01596.hmm"


# =============================================================================
# SECTION 10: PARSE .domtbl → TRANSCRIPT ID LISTS
# =============================================================================
# domtbl column 1 = protein ID (e.g. TRINITY_DN123_c0_g1_i1.p1)
# Strip trailing .p1 / .p2 to recover the Trinity transcript ID used in RSEM.

parse_domtbl_to_ids() {
    # $1 = input .domtbl, $2 = output .txt
    grep -v '^#' "$1" \
        | awk '{print $1}' \
        | sed -E 's/\.p[0-9]+$//' \
        | sort -u \
        > "$2"
}

echo ""
echo "Parsing domain tables → transcript ID lists..."

parse_domtbl_to_ids "${OUTPUT_DIR}/cyp450.domtbl" "${CYP_OUT}"

# Merge all three OMT domain hits into ONE candidate set
OMT_TMP="${OUTPUT_DIR}/omt_all.tmp"
cat "${OUTPUT_DIR}/omt_a.domtbl" \
    "${OUTPUT_DIR}/omt_b.domtbl" \
    "${OUTPUT_DIR}/omt_ccoa.domtbl" \
    | grep -v '^#' \
    | awk '{print $1}' \
    | sed -E 's/\.p[0-9]+$//' \
    | sort -u \
    > "${OMT_TMP}"
mv "${OMT_TMP}" "${OMT_OUT}"

N_CYP=$(wc -l < "${CYP_OUT}" | tr -d ' ')
N_OMT=$(wc -l < "${OMT_OUT}" | tr -d ' ')

echo ""
echo "=========================================="
echo "HMMER scan complete"
echo "=========================================="
echo "CYP450 candidates: ${N_CYP} transcript IDs → ${CYP_OUT}"
echo "OMT candidates:    ${N_OMT} transcript IDs → ${OMT_OUT}"
echo ""
echo "Raw domain tables (keep for Script 04 e-value filtering):"
echo "  ${OUTPUT_DIR}/cyp450.domtbl"
echo "  ${OUTPUT_DIR}/omt_a.domtbl"
echo "  ${OUTPUT_DIR}/omt_b.domtbl"
echo "  ${OUTPUT_DIR}/omt_ccoa.domtbl"
echo "=========================================="
