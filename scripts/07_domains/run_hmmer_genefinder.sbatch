#!/bin/bash
# =============================================================================
# SCRIPT:  run_hmmer_genefinder.sbatch
#
# WHAT THIS DOES:
#   Searches TransDecoder protein sequences for conserved Pfam domains that
#   DEFINE cytochrome P450 (CYP) and O-methyltransferase (OMT) families.
#
#   HMMER profiles are trained on hundreds of known family members, so they
#   catch even divergent nutmeg genes better than keyword search alone.
#   This is your main, most sensitive first-pass filter for CYP/OMT candidates.
#
# WHEN TO RUN (Step 9 in docs/README_trinity.md):
#   AFTER TransDecoder (Step 7):
#     sbatch scripts/03_assembly/run_transdecoder.sbatch MF
#
# HOW TO SUBMIT (MF nutmeg example):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF
#
# TO RE-RUN FROM SCRATCH:
#   sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force
#
# INPUTS REQUIRED:
#   1. TransDecoder protein FASTA (.pep)
#      MF: 01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
#   2. HMMER in rnaseq conda env (hmmsearch, hmmpress)
#   3. Pfam HMM profiles (auto-downloaded on first run):
#        PF00067  — cytochrome P450 (CYP)
#        PF00891  — Methyltransf_2 (OMT)
#        PF08100  — OMT dimerisation domain
#        PF01596  — CCoA-type OMT
#
# OUTPUT (MF example):
#   06_analysis/hmmer_genefinder_MF/cyp450_hmmer_ids.txt
#   06_analysis/hmmer_genefinder_MF/omt_hmmer_ids.txt
#   06_analysis/hmmer_genefinder_MF/*.domtbl
#
# NOTE ON FILE NAMING:
#   Like run_transdecoder.sbatch and run_rsem.sbatch, this is ONE self-contained
#   file. You only ever submit the .sbatch file — there is no separate .sh helper.
#
# PATHS come from scripts/config.sh (BASE_DIR, PROCESSED_DIR, etc.)
# =============================================================================


# ---- SLURM: tell the cluster how much time/memory this job needs ----
# HMMER is fast (minutes) — we request 2 hours to be safe.
#SBATCH --job-name=hmmer_genefinder
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o hmmer_genefinder_%j.out
#SBATCH -e hmmer_genefinder_%j.err

# set -e  = stop if any command fails
# set -u  = stop if we use a variable that was never set
# set -o pipefail = treat a failed command in a pipe as a failure
set -euo pipefail


# =============================================================================
# SECTION 1: LOAD config.sh
# =============================================================================
# config.sh defines BASE_DIR, PROCESSED_DIR, CONDA_ENV, and get_sample_info()

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
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# SECTION 2: READ WHAT THE USER TYPED AFTER sbatch
# =============================================================================
# Example: sbatch run_hmmer_genefinder.sbatch MF --force
#   $1 = MF       (species code — required)
#   --force       (optional — re-run even if output files already exist)

SPECIES_CODE="${1:-}"   # e.g. MF for nutmeg
FORCE=0                 # 0 = skip if outputs exist; 1 = re-run

# Look through ALL arguments for --force
for arg in "$@"; do
    if [[ "${arg}" == "--force" ]]; then
        FORCE=1
    fi
done

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: You must give a species code (e.g. MF for nutmeg)."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch <CODE> [--force]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force"
    echo ""
    echo "Required input (MF example):"
    echo "  ${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep"
    echo ""
    echo "Run TransDecoder first if that file is missing:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch MF"
    exit 1
fi

# get_sample_info is defined in config.sh — sets SPECIES_NAME, GENOME_TYPE, etc.
get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# SECTION 3: ACTIVATE CONDA (HMMER lives in rnaseq env)
# =============================================================================

echo ""
echo "Activating conda environment: ${CONDA_ENV}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"


# =============================================================================
# SECTION 4: SET FILE PATHS
# =============================================================================
# TransDecoder writes .pep next to the CD-HIT FASTA.
# Those transcript IDs match RSEM counts and PyDESeq2 — do not use a different FASTA.

PEP_FILE="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta.transdecoder.pep"
OUTPUT_DIR="${ANALYSIS_DIR}/hmmer_genefinder_${SPECIES_CODE_UPPER}"
HMM_DIR="${DB_DIR}/hmmerdb/pfam_profiles"

CYP_OUT="${OUTPUT_DIR}/cyp450_hmmer_ids.txt"
OMT_OUT="${OUTPUT_DIR}/omt_hmmer_ids.txt"

# Pfam profile files — one HMM per gene family domain
# Keys (left) are output file prefixes; values (right) are Pfam filenames.
declare -A HMM_FILES=(
    [cyp450]="PF00067.hmm"
    [omt_a]="PF00891.hmm"
    [omt_b]="PF08100.hmm"
    [omt_ccoa]="PF01596.hmm"
)


# =============================================================================
# SECTION 5: CHECK THAT HMMER IS INSTALLED
# =============================================================================

echo ""
echo "=========================================="
echo "HMMER CYP/OMT gene finder (Script 02)"
echo "=========================================="
echo "Species:     ${SPECIES_NAME} (${SPECIES_CODE_UPPER})"
echo "Proteins:    ${PEP_FILE}"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Force:       ${FORCE}"
echo "------------------------------------------"

if ! command -v hmmsearch &> /dev/null; then
    echo "ERROR: hmmsearch not found on PATH."
    echo "Install HMMER in your conda env:"
    echo "  conda activate ${CONDA_ENV}"
    echo "  conda install -c bioconda hmmer"
    exit 1
fi

# hmmsearch -h prints version on line 2 (there is no --version flag)
HMMER_VERSION=$(hmmsearch -h 2>&1 | grep -m1 "^# HMMER" | sed 's/^# //')
echo "HMMER:       ${HMMER_VERSION}"

if ! command -v hmmpress &> /dev/null; then
    echo "ERROR: hmmpress not found on PATH."
    echo "hmmpress indexes HMM files before searching — install with:"
    echo "  conda install -c bioconda hmmer"
    exit 1
fi


# =============================================================================
# SECTION 6: CHECK THAT THE .pep INPUT FILE EXISTS
# =============================================================================

if [[ ! -f "${PEP_FILE}" ]]; then
    echo ""
    echo "ERROR: TransDecoder protein file not found:"
    echo "  ${PEP_FILE}"
    echo ""
    echo "Run TransDecoder first (Step 7 in docs/README_trinity.md):"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch ${SPECIES_CODE_UPPER}"
    exit 1
fi

# Count how many protein sequences are in the file (lines starting with >)
N_PROTEINS=$(grep -c '^>' "${PEP_FILE}" || true)
echo "Proteins:    ${N_PROTEINS} sequences in .pep file"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${HMM_DIR}"


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
    echo "--force set: re-running hmmsearch and overwriting outputs."
fi


# =============================================================================
# SECTION 8: DOWNLOAD PFAM HMM PROFILES IF MISSING
# =============================================================================
# Each Pfam family has a "profile HMM" — a statistical model of the domain.
# We download them once and cache them in 07_NRdatabase/hmmerdb/pfam_profiles/.

download_hmm() {
    # Argument $1 = Pfam accession, e.g. PF00067
    local pfam_acc="$1"
    local hmm_path="${HMM_DIR}/${pfam_acc}.hmm"

    # If we already downloaded this file and it looks valid, reuse it
    if [[ -f "${hmm_path}" ]] && [[ "$(wc -c < "${hmm_path}" | tr -d ' ')" -gt 100 ]]; then
        echo "  Using cached: ${hmm_path}"
        return 0
    fi

    echo "  Downloading ${pfam_acc} from InterPro/Pfam..."

    # Try two URLs — Pfam mirrors change sometimes
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
            # InterPro sometimes returns a gzip file — decompress if needed
            if file "${hmm_path}" | grep -q gzip; then
                gunzip -c "${hmm_path}" > "${hmm_path}.tmp"
                mv "${hmm_path}.tmp" "${hmm_path}"
            fi
            # A valid HMM file always starts with "HMMER"
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
    # Strip ".hmm" to get accession, e.g. PF00067
    pfam_acc="${HMM_FILES[$key]%.hmm}"
    download_hmm "${pfam_acc}"
done


# =============================================================================
# SECTION 9: PREPARE HMM FILES WITH hmmpress
# =============================================================================
# hmmpress builds index files (.h3m, .h3i, etc.) so hmmsearch runs faster.
# You only need to run this once per HMM file.

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
# SECTION 10: RUN hmmsearch
# =============================================================================
# hmmsearch = search ONE HMM profile against MANY protein sequences.
#
# --cut_ga  = use Pfam's curated "gathering threshold" (recommended)
# --domtblout = write a domain hit table we parse later
# --cpu     = number of threads (SLURM gives us 8)
# --noali   = skip printing alignments (saves disk space)

run_hmmsearch() {
    # $1 = short label for output files (e.g. cyp450)
    # $2 = Pfam HMM filename (e.g. PF00067.hmm)
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

echo ""
echo "Running hmmsearch (this usually takes a few minutes)..."

run_hmmsearch "cyp450"   "PF00067.hmm"
run_hmmsearch "omt_a"    "PF00891.hmm"
run_hmmsearch "omt_b"    "PF08100.hmm"
run_hmmsearch "omt_ccoa" "PF01596.hmm"


# =============================================================================
# SECTION 11: PARSE .domtbl FILES → TRANSCRIPT ID LISTS
# =============================================================================
# TransDecoder names proteins like:  TRINITY_DN123_c0_g1_i1.p1
# RSEM / PyDESeq2 use transcript IDs: TRINITY_DN123_c0_g1_i1
# We strip the trailing .p1 / .p2 so the IDs match your count matrix.

parse_domtbl_to_ids() {
    # $1 = input .domtbl file
    # $2 = output .txt file (one transcript ID per line)
    grep -v '^#' "$1" \
        | awk '{print $1}' \
        | sed -E 's/\.p[0-9]+$//' \
        | sort -u \
        > "$2"
}

echo ""
echo "Parsing domain tables → transcript ID lists..."

# CYP450: one Pfam profile → one ID list
parse_domtbl_to_ids "${OUTPUT_DIR}/cyp450.domtbl" "${CYP_OUT}"

# OMT: three Pfam profiles → merge into ONE combined ID list
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


# =============================================================================
# SECTION 12: PRINT SUMMARY
# =============================================================================

N_CYP=$(wc -l < "${CYP_OUT}" | tr -d ' ')
N_OMT=$(wc -l < "${OMT_OUT}" | tr -d ' ')

echo ""
echo "=========================================="
echo "HMMER gene finder complete"
echo "=========================================="
echo "CYP450 candidates: ${N_CYP} transcript IDs"
echo "  → ${CYP_OUT}"
echo ""
echo "OMT candidates:    ${N_OMT} transcript IDs"
echo "  → ${OMT_OUT}"
echo ""
echo "Raw domain tables (keep for later e-value filtering):"
echo "  ${OUTPUT_DIR}/cyp450.domtbl"
echo "  ${OUTPUT_DIR}/omt_a.domtbl"
echo "  ${OUTPUT_DIR}/omt_b.domtbl"
echo "  ${OUTPUT_DIR}/omt_ccoa.domtbl"
echo "=========================================="
