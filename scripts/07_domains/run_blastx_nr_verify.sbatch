#!/bin/bash
# =============================================================================
# SCRIPT:  run_blastx_nr_verify.sbatch   (CYP/OMT pipeline Script 05 - NAMING step)
#
# WHAT THIS DOES:
#   Puts a real NAME + ORGANISM on your candidate transcripts by searching them
#   against NCBI nr with BLASTX (nucleotide transcript -> protein database).
#   It turns "TRINITY_DN10886_c0_g1_i2" into something readable like
#   "cytochrome P450 CYP736A-like, Cinnamomum, 68% identity".
#
#   This does NOT find new candidates. Your wide net (HMMER + InterProScan)
#   already found the family members. This step just names them.
#
# TWO MODES (same script, different folders so results never collide):
#   (default)          name only your candidate transcripts   -> FAST
#   --whole-assembly   name EVERY transcript in the assembly   -> SLOW (hours)
#
#   Output folders:
#     candidates    -> ${ANALYSIS_DIR}/blastx_nr_MF_candidates/
#     whole assembly-> ${ANALYSIS_DIR}/blastx_nr_MF_wholeassembly/
#
# HOW TO SUBMIT:
#   FAST (run this now):
#     sbatch scripts/07_domains/run_blastx_nr_verify.sbatch MF
#   WHOLE ASSEMBLY (save for later):
#     sbatch scripts/07_domains/run_blastx_nr_verify.sbatch MF --whole-assembly
#
# INPUTS:
#   1. CD-HIT assembly (same IDs as your candidate lists):
#        ${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta
#   2. Candidate transcript ID lists (TRANSCRIPT level, with the _i suffix).
#      Defaults point at your combined HMMER+BLAST lists. Override with env vars
#      if your paths differ:
#        CYP_IDS=/path/to/combined_cyp450_ids.txt
#        OMT_IDS=/path/to/combined_omt_ids.txt
#      NOTE: use the *_ids.txt (transcripts), NOT the *_genes.txt files.
#      You blast actual sequences, and only transcripts have sequences.
#   3. A formatted nr BLAST database. Set BLAST_NR_DB, or rely on $BLASTDB and
#      leave it as "nr". Find it first (see the echo $BLASTDB check below).
#   4. blastx, blastdbcmd, seqkit in the rnaseq conda env.
# =============================================================================

#SBATCH --job-name=blastx_nr
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00        # candidate mode finishes fast; whole-assembly may need more
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis/00_3_MF
#SBATCH -o blastx_nr_%j.out
#SBATCH -e blastx_nr_%j.err

set -euo pipefail


# ---- SECTION 1: load config.sh --------------------------------------------
CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../../config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    [[ -f "${candidate}" ]] && { CONFIG_FILE="${candidate}"; break; }
done
[[ -z "${CONFIG_FILE}" ]] && { echo "ERROR: config.sh not found. Submit from the repo folder."; exit 1; }
echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# ---- SECTION 2: args -------------------------------------------------------
SPECIES_CODE="${1:-}"
FORCE=0
WHOLE_ASSEMBLY=0
for arg in "$@"; do
    [[ "${arg}" == "--force" ]] && FORCE=1
    [[ "${arg}" == "--whole-assembly" ]] && WHOLE_ASSEMBLY=1
done
[[ -z "${SPECIES_CODE}" ]] && { echo "ERROR: give a species code, e.g.  sbatch ... MF"; exit 1; }

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# ---- SECTION 3: conda ------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"


# ---- SECTION 4: paths ------------------------------------------------------
ASSEMBLY_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"

# nr database. If your cluster sets $BLASTDB to the folder holding nr.*, the
# plain name "nr" works. Otherwise set BLAST_NR_DB to the full path.
BLAST_NR_DB="${BLAST_NR_DB:-nr}"

# Candidate TRANSCRIPT lists (with _i suffix). Defaults = your combined lists.
CYP_IDS="${CYP_IDS:-${BASE_DIR}/combined_cyp450_ids.txt}"
OMT_IDS="${OMT_IDS:-${BASE_DIR}/combined_omt_ids.txt}"

# Mode-tagged output folder so the two runs never overwrite each other.
if [[ ${WHOLE_ASSEMBLY} -eq 1 ]]; then
    MODE_TAG="wholeassembly"
else
    MODE_TAG="candidates"
fi
OUTPUT_DIR="${ANALYSIS_DIR}/blastx_nr_${SPECIES_CODE_UPPER}_${MODE_TAG}"
mkdir -p "${OUTPUT_DIR}"

EVALUE="1e-5"
MAX_TARGET_SEQS=25          # set generously; we sort for the true best hit ourselves
THREADS="${SLURM_CPUS_PER_TASK:-16}"
OUTFMT="6 qseqid sacc pident length evalue bitscore qcovs staxids stitle"


# ---- SECTION 5: checks -----------------------------------------------------
echo ""
echo "=========================================="
echo "BLASTX vs nr  -  NAMING step  (${MODE_TAG})"
echo "=========================================="
echo "Species:  ${SPECIES_NAME} (${SPECIES_CODE_UPPER})"
echo "Assembly: ${ASSEMBLY_FASTA}"
echo "nr DB:    ${BLAST_NR_DB}"
echo "Output:   ${OUTPUT_DIR}"
echo "------------------------------------------"

command -v blastx >/dev/null || { echo "ERROR: blastx not found. conda install -c bioconda blast"; exit 1; }
command -v seqkit >/dev/null || { echo "ERROR: seqkit not found. conda install -c bioconda seqkit"; exit 1; }
echo "blastx: $(blastx -version 2>&1 | head -1)"
echo "seqkit: $(seqkit version 2>&1 | head -1)"

[[ -f "${ASSEMBLY_FASTA}" ]] || { echo "ERROR: assembly not found: ${ASSEMBLY_FASTA}"; exit 1; }

if ! blastdbcmd -db "${BLAST_NR_DB}" -info >/dev/null 2>&1; then
    echo "ERROR: nr database '${BLAST_NR_DB}' not reachable."
    echo "  Find it:   echo \$BLASTDB   (and:  ls \$BLASTDB | grep nr )"
    echo "  Then set:  BLAST_NR_DB=/full/path/to/nr sbatch ... ${SPECIES_CODE_UPPER}"
    echo "  Or ask ARC helpdesk where the shared nr lives before downloading."
    exit 1
fi

if [[ ${WHOLE_ASSEMBLY} -eq 0 ]]; then
    [[ -f "${CYP_IDS}" ]] || { echo "ERROR: CYP list not found: ${CYP_IDS} (set CYP_IDS=...)"; exit 1; }
    [[ -f "${OMT_IDS}" ]] || { echo "ERROR: OMT list not found: ${OMT_IDS} (set OMT_IDS=...)"; exit 1; }
    echo "CYP candidates: $(wc -l < "${CYP_IDS}" | tr -d ' ') transcripts"
    echo "OMT candidates: $(wc -l < "${OMT_IDS}" | tr -d ' ') transcripts"
fi


# ---- SECTION 6: extract candidate transcripts by ID ------------------------
extract_candidates() {
    local ids="$1" out_fa="$2" label="$3"
    if [[ ${FORCE} -eq 0 && -s "${out_fa}" ]]; then
        echo "  [${label}] candidate FASTA exists, skipping extract"
        return 0
    fi
    # seqkit grep -f: keep FASTA records whose ID is in the ID-list file.
    seqkit grep -f "${ids}" "${ASSEMBLY_FASTA}" > "${out_fa}"
    echo "  [${label}] pulled $(grep -c '^>' "${out_fa}" || true) transcripts"
}


# ---- SECTION 7: blastx + best-hit-per-transcript ---------------------------
run_blastx() {
    local query_fa="$1" out_tsv="$2" best_tsv="$3" label="$4"

    if [[ ${FORCE} -eq 0 && -s "${out_tsv}" ]]; then
        echo "  [${label}] blastx output exists, skipping (--force to redo)"
    else
        echo "  [${label}] blastx vs nr ... (slow part)"
        blastx \
            -query "${query_fa}" \
            -db "${BLAST_NR_DB}" \
            -out "${out_tsv}" \
            -evalue "${EVALUE}" \
            -max_target_seqs "${MAX_TARGET_SEQS}" \
            -num_threads "${THREADS}" \
            -outfmt "${OUTFMT}"
    fi

    # Best hit per transcript = highest bitscore (col 6). This is the "name".
    if [[ -s "${out_tsv}" ]]; then
        sort -t $'\t' -k1,1 -k6,6nr "${out_tsv}" | awk -F'\t' '!seen[$1]++' > "${best_tsv}"
        echo "  [${label}] named $(cut -f1 "${best_tsv}" | sort -u | wc -l | tr -d ' ') transcripts -> ${best_tsv}"
    else
        : > "${best_tsv}"
        echo "  [${label}] no hits"
    fi
}

process_family() {
    local label="$1" ids="$2"
    local fa="${OUTPUT_DIR}/${label}_query.fasta"
    local tsv="${OUTPUT_DIR}/${label}_vs_nr.tsv"
    local best="${OUTPUT_DIR}/${label}_vs_nr.besthit.tsv"
    echo ""
    echo "### ${label}"
    if [[ ${WHOLE_ASSEMBLY} -eq 1 ]]; then
        run_blastx "${ASSEMBLY_FASTA}" "${tsv}" "${best}" "${label}"
    else
        extract_candidates "${ids}" "${fa}" "${label}"
        run_blastx "${fa}" "${tsv}" "${best}" "${label}"
    fi
}

# In whole-assembly mode both families search the same input, so run once.
if [[ ${WHOLE_ASSEMBLY} -eq 1 ]]; then
    process_family "allhits" ""
else
    process_family "cyp450" "${CYP_IDS}"
    process_family "omt"    "${OMT_IDS}"
fi


# ---- SECTION 8: summary ----------------------------------------------------
echo ""
echo "=========================================="
echo "Naming complete  (${MODE_TAG})"
echo "Look at the *.besthit.tsv files. Columns:"
echo "  qseqid  sacc  pident  length  evalue  bitscore  qcovs  staxids  stitle"
echo "The stitle column (last) is the readable protein name."
echo "Output folder: ${OUTPUT_DIR}"
echo "=========================================="
