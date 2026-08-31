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
#   (default)          name only your candidate transcripts   -> hours–days vs nr
#   --whole-assembly   name EVERY transcript in the assembly  -> days vs nr
#
#   Output folders:
#     candidates    -> ${ANALYSIS_DIR}/blastx_nr_MF_candidates/
#     whole assembly-> ${ANALYSIS_DIR}/blastx_nr_MF_wholeassembly/
#
# HOW TO SUBMIT:
#   FAST (candidates only — start here):
#     sbatch scripts/07_domains/run_blastx_nr_verify.sbatch MF
#   WHOLE ASSEMBLY (save for later):
#     sbatch scripts/07_domains/run_blastx_nr_verify.sbatch MF --whole-assembly
#
# INPUTS:
#   1. CD-HIT assembly (same IDs as your candidate lists):
#        ${PROCESSED_DIR}/00_6_cdhit/MF_trinity_cdhit95.fasta
#   2. Candidate ID lists. Defaults sit next to your SLURM logs, in the
#      species folder (MF = 00_3_MF):
#        ${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450_ids.txt
#        ${ANALYSIS_DIR}/${SPECIES_DIR}/combined_omt_ids.txt
#      which for MF is:
#        /projects/tholl_lab_1/daisy_analysis/06_analysis/00_3_MF/combined_cyp450_ids.txt
#        /projects/tholl_lab_1/daisy_analysis/06_analysis/00_3_MF/combined_omt_ids.txt
#      Use combined_*_ids.txt (transcripts, with _i). Do NOT use *_genes.txt.
#      Override only with a REAL path if your files live somewhere else:
#        CYP_IDS=/projects/tholl_lab_1/daisy_analysis/06_analysis/00_3_MF/combined_cyp450_ids.txt \
#        OMT_IDS=/projects/tholl_lab_1/daisy_analysis/06_analysis/00_3_MF/combined_omt_ids.txt \
#        sbatch scripts/07_domains/run_blastx_nr_verify.sbatch MF
#   3. A formatted nr BLAST database. This script sets BLASTDB to
#        ${DB_DIR}/blastdb
#      which is /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/blastdb
#      (same folder your other BLAST jobs use). Override with:
#        BLAST_NR_DB=/full/path/to/nr  sbatch ... MF
#   4. blastx, blastdbcmd, seqkit in the rnaseq conda env.
# =============================================================================

#SBATCH --job-name=blastx_nr
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G             # nr is huge; more RAM = less disk thrashing
#SBATCH --time=6-00:00:00      # nr vs 757 transcripts can take days
#                              # Override if needed:
#                              #   sbatch --time=6-00:00:00 ... MF
#                              # Whole-assembly vs nr often needs 6 days.
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

# Tell BLAST where the downloaded databases live (same as run_blastx_discovery).
# BLAST then looks for files named nr.* inside this folder.
# DB_DIR comes from config.sh = /projects/tholl_lab_1/daisy_analysis/07_NRdatabase
export BLASTDB="${DB_DIR}/blastdb"

# Database NAME (not a folder). Override with BLAST_NR_DB=/full/path/to/nr if needed.
BLAST_NR_DB="${BLAST_NR_DB:-nr}"

# Candidate lists. Default = the files you already have in 00_3_MF/
# (SPECIES_DIR comes from get_sample_info: MF -> 00_3_MF).
# If a leftover CYP_IDS=/full/path/... is still in your shell, unset it:
#   unset CYP_IDS OMT_IDS
CYP_IDS="${CYP_IDS:-${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450_ids.txt}"
OMT_IDS="${OMT_IDS:-${ANALYSIS_DIR}/${SPECIES_DIR}/combined_omt_ids.txt}"

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
echo "BLASTDB:  ${BLASTDB}"
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
    echo "  BLAST looked in: ${BLASTDB}"
    echo "  Check that nr is there:"
    echo "    ls ${BLASTDB}/nr.pal"
    echo "  If nr lives somewhere else, resubmit with the full path:"
    echo "    BLAST_NR_DB=/full/path/to/nr sbatch scripts/07_domains/run_blastx_nr_verify.sbatch ${SPECIES_CODE_UPPER}"
    echo "  If nr was never downloaded, run:"
    echo "    sbatch scripts/setup/download_blastdb.sbatch"
    exit 1
fi

if [[ ${WHOLE_ASSEMBLY} -eq 0 ]]; then
    if [[ ! -f "${CYP_IDS}" ]]; then
        echo "ERROR: CYP list not found: ${CYP_IDS}"
        echo "  That path is not a real file. Common causes:"
        echo "    1. You copied a placeholder like CYP_IDS=/full/path/...  (unset it:  unset CYP_IDS OMT_IDS)"
        echo "    2. The combined list lives somewhere else. Find it with:"
        echo "         ls ${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450_ids.txt"
        echo "         ls ${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450*.txt"
        echo "  Then resubmit with the REAL path, for example:"
        echo "    CYP_IDS=${ANALYSIS_DIR}/${SPECIES_DIR}/combined_cyp450_ids.txt \\"
        echo "    OMT_IDS=${ANALYSIS_DIR}/${SPECIES_DIR}/combined_omt_ids.txt \\"
        echo "    sbatch scripts/07_domains/run_blastx_nr_verify.sbatch ${SPECIES_CODE_UPPER}"
        exit 1
    fi
    if [[ ! -f "${OMT_IDS}" ]]; then
        echo "ERROR: OMT list not found: ${OMT_IDS}"
        echo "  Unset a leftover placeholder with:  unset OMT_IDS"
        echo "  Or pass the real file:  OMT_IDS=${ANALYSIS_DIR}/${SPECIES_DIR}/combined_omt_ids.txt sbatch ... ${SPECIES_CODE_UPPER}"
        exit 1
    fi
    echo "CYP candidates: $(wc -l < "${CYP_IDS}" | tr -d ' ') lines in ${CYP_IDS}"
    echo "OMT candidates: $(wc -l < "${OMT_IDS}" | tr -d ' ') lines in ${OMT_IDS}"
    echo "First CYP ID:   $(head -1 "${CYP_IDS}")"
    echo "First OMT ID:   $(head -1 "${OMT_IDS}")"
fi


# ---- SECTION 6: extract candidate transcripts by ID ------------------------
extract_candidates() {
    local ids="$1" out_fa="$2" label="$3"
    local n_ids n_fa
    # Count non-empty ID lines (the list we SHOULD search).
    n_ids=$(grep -cve '^[[:space:]]*$' "${ids}" || true)
    n_fa=0
    if [[ -s "${out_fa}" ]]; then
        n_fa=$(grep -c '^>' "${out_fa}" || true)
    fi

    # Skip only when the FASTA already has the SAME number of sequences as IDs.
    # Old runs extracted the BLAST-genefinder list (79 CYP) and later jobs
    # skipped forever even after combined_cyp450_ids.txt grew to 757.
    if [[ ${FORCE} -eq 0 && -s "${out_fa}" && "${n_fa}" -eq "${n_ids}" ]]; then
        echo "  [${label}] candidate FASTA exists (${n_fa} seqs = ${n_ids} IDs), skipping extract"
        return 0
    fi
    echo "  [${label}] extracting ${n_ids} IDs from assembly (old FASTA had ${n_fa} seqs)"
    # seqkit grep -f: keep FASTA records whose ID is in the ID-list file.
    seqkit grep -f "${ids}" "${ASSEMBLY_FASTA}" > "${out_fa}"
    echo "  [${label}] pulled $(grep -c '^>' "${out_fa}" || true) transcripts"
}


# ---- SECTION 7: blastx + best-hit-per-transcript ---------------------------
run_blastx() {
    local query_fa="$1" out_tsv="$2" best_tsv="$3" label="$4"

    # Re-run BLAST if the query FASTA is newer than the TSV (extract just grew).
    local rerun=0
    if [[ ${FORCE} -eq 1 ]]; then
        rerun=1
    elif [[ ! -s "${out_tsv}" ]]; then
        rerun=1
    elif [[ -f "${query_fa}" && "${query_fa}" -nt "${out_tsv}" ]]; then
        echo "  [${label}] query FASTA is newer than blast TSV — re-running blastx"
        rerun=1
    fi

    if [[ ${rerun} -eq 0 ]]; then
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
