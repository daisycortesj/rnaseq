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
#   2. HMMER in rnaseq conda env (hmmsearch)
#   3. Pfam HMM profiles (auto-downloaded on first run):
#        PF00067  — cytochrome P450 (CYP)          [catalytic]
#        PF00891  — Methyltransf_2 (Class II OMT)  [catalytic]
#        PF01596  — Methyltransf_3 (CCoAOMT-like)  [catalytic]
#        PF08100  — Dimerisation domain            [STRUCTURAL, not catalytic]
#
# OUTPUT (MF example) — TWO LEVELS, see CHANGE 3 below:
#   06_analysis/hmmer_genefinder_MF/cyp450_hmmer_ids.txt        (isoform-level)
#   06_analysis/hmmer_genefinder_MF/cyp450_hmmer_gene_ids.txt   (GENE-level) <-- use this
#   06_analysis/hmmer_genefinder_MF/omt_hmmer_ids.txt           (isoform-level)
#   06_analysis/hmmer_genefinder_MF/omt_hmmer_gene_ids.txt      (GENE-level) <-- use this
#   06_analysis/hmmer_genefinder_MF/omt_dimer_only_gene_ids.txt (for inspection)
#   06_analysis/hmmer_genefinder_MF/*.domtbl
#
# PATHS come from scripts/config.sh (BASE_DIR, PROCESSED_DIR, etc.)
#
# =============================================================================
# ##### REVISION HISTORY #####
# 2026-07-17  Fixed 5 issues. Each is marked inline with ">>> CHANGE N".
#
#   CHANGE 1 (Section 4)  Added GENE-level output paths.
#   CHANGE 2 (Section 7)  Skip-check now also requires the gene-level files,
#                         otherwise a re-run would skip and never build them.
#   CHANGE 3 (Section 11) *** THE REAL BUG ***  The old comment claimed
#                         "RSEM / PyDESeq2 use transcript IDs". FALSE. RSEM was
#                         run with --transcript-to-gene-map, so the count matrix
#                         is GENE-level (TRINITY_DN0_c0_g1). The script stripped
#                         ".p1" but left "_i1", producing isoform IDs that could
#                         NEVER match. Overlap with the DE table was literally 0.
#                         Now writes both isoform- and gene-level lists.
#   CHANGE 4 (Section 11) OMT definition tightened. PF08100 is the DIMERISATION
#                         domain -- structural, not catalytic. Including it in a
#                         plain union let non-OMTs into the list. Core OMT is now
#                         PF00891 OR PF01596; PF08100-only hits are written to a
#                         separate file for you to inspect, not silently merged.
#   CHANGE 5 (Sections 5,8,9)  Removed dead code: hmmpress (only needed for
#                         hmmscan, we run hmmsearch) and the retired
#                         pfam.xfam.org fallback URL.
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

# >>> CHANGE 5a: guard against `sbatch script.sbatch --force` with no species.
#     Without this, SPECIES_CODE would literally become the string "--force"
#     and you'd get a confusing error much further down.
if [[ "${SPECIES_CODE}" == "--force" ]]; then
    SPECIES_CODE=""
fi

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: You must give a species code (e.g. MF for nutmeg)."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch <CODE> [--force]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF"
    echo "  sbatch scripts/07_domains/run_hmmer_genefinder.sbatch MF --force"
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

PEP_FILE="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta.transdecoder.pep"
OUTPUT_DIR="${ANALYSIS_DIR}/hmmer_genefinder_${SPECIES_CODE_UPPER}"
HMM_DIR="${DB_DIR}/hmmerdb/pfam_profiles"

# ---------------------------------------------------------------------------
# >>> CHANGE 1: two levels of output instead of one.
#
# Trinity IDs nest like this:
#
#     TRINITY_DN123_c0_g1_i1.p1
#                        |   +--- ORF number    (added by TransDecoder)
#                        +------- isoform       (added by Trinity)
#     TRINITY_DN123_c0_g1  <---- GENE  (this is where RSEM stops)
#
# RSEM was run with --transcript-to-gene-map, so it summed all isoforms per
# gene. Your count matrix and PyDESeq2 results are therefore GENE-level.
#
# The .domtbl records which ISOFORM carried the domain -- that is real
# information and we do not want to throw it away. So we write both:
#   *_hmmer_ids.txt       isoform-level  (raw record; keep for the record)
#   *_hmmer_gene_ids.txt  gene-level     (JOINS TO THE COUNT MATRIX; use this)
# ---------------------------------------------------------------------------

# ISOFORM-level (raw record: which isoform actually carried the domain)
CYP_OUT="${OUTPUT_DIR}/cyp450_hmmer_ids.txt"
OMT_OUT="${OUTPUT_DIR}/omt_hmmer_ids.txt"

# GENE-level (joins to the RSEM/PyDESeq2 count matrix — USE THESE DOWNSTREAM)
CYP_GENE_OUT="${OUTPUT_DIR}/cyp450_hmmer_gene_ids.txt"
OMT_GENE_OUT="${OUTPUT_DIR}/omt_hmmer_gene_ids.txt"

# >>> CHANGE 4a: PF08100-only hits get their own file (see CHANGE 4 in Sec 11)
OMT_DIMER_ONLY_OUT="${OUTPUT_DIR}/omt_dimer_only_gene_ids.txt"

# Pfam profile files — one HMM per gene family domain
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
echo "HMMER CYP/OMT gene finder"
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

# >>> CHANGE 5b: REMOVED the hmmpress check that used to live here.
#     We do not use hmmpress at all -- see CHANGE 5c in Section 9.


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

N_PROTEINS=$(grep -c '^>' "${PEP_FILE}" || true)
echo "Proteins:    ${N_PROTEINS} sequences in .pep file"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${HMM_DIR}"


# =============================================================================
# SECTION 7: SKIP IF OUTPUTS ALREADY EXIST (--force to overwrite)
# =============================================================================
# >>> CHANGE 2: the skip test now ALSO requires the gene-level files.
#     Why this matters: without it, the old isoform-level files (which already
#     exist on disk from your June run) would satisfy the test, the script
#     would exit here, and the new gene-level files would never be created.
#     You would see "skipping" and think it worked.

if [[ "${FORCE}" -eq 0 ]] \
   && [[ -f "${CYP_OUT}" ]] && [[ -f "${OMT_OUT}" ]] \
   && [[ -f "${CYP_GENE_OUT}" ]] && [[ -f "${OMT_GENE_OUT}" ]]; then
    echo ""
    echo "Outputs already exist — skipping (use --force to re-run)."
    echo ""
    echo "Candidate counts:"
    echo "  CYP450: $(wc -l < "${CYP_OUT}"      | tr -d ' ') isoforms"
    echo "  CYP450: $(wc -l < "${CYP_GENE_OUT}" | tr -d ' ') genes    <-- use this"
    echo "  OMT:    $(wc -l < "${OMT_OUT}"      | tr -d ' ') isoforms"
    echo "  OMT:    $(wc -l < "${OMT_GENE_OUT}" | tr -d ' ') genes    <-- use this"
    exit 0
fi

if [[ "${FORCE}" -eq 1 ]]; then
    echo ""
    echo "--force set: re-running hmmsearch and overwriting outputs."
fi


# =============================================================================
# SECTION 8: DOWNLOAD PFAM HMM PROFILES IF MISSING
# =============================================================================
# Each Pfam family has a "profile HMM" — a statistical model of the domain,
# built from an alignment of hundreds of known family members. That is why it
# finds divergent nutmeg genes that a plain BLAST search would miss.

download_hmm() {
    # Argument $1 = Pfam accession, e.g. PF00067
    local pfam_acc="$1"
    local hmm_path="${HMM_DIR}/${pfam_acc}.hmm"

    # If we already downloaded this file and it looks valid, reuse it
    if [[ -f "${hmm_path}" ]] && [[ "$(wc -c < "${hmm_path}" | tr -d ' ')" -gt 100 ]]; then
        echo "  Using cached: ${hmm_path}"
        return 0
    fi

    echo "  Downloading ${pfam_acc} from InterPro..."

    # >>> CHANGE 5c: removed the second URL (pfam.xfam.org). That site was
    #     retired when Pfam moved into InterPro, so it could never succeed.
    #     Dead code in a retry loop is worse than no retry loop: it makes you
    #     think you have a fallback when you do not.
    local url="https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pfam_acc}?annotation=hmm"
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
    pfam_acc="${HMM_FILES[$key]%.hmm}"   # strip ".hmm" -> PF00067
    download_hmm "${pfam_acc}"
done


# =============================================================================
# SECTION 9: (REMOVED — used to run hmmpress)
# =============================================================================
# >>> CHANGE 5d: the old Section 9 ran hmmpress on each .hmm file, with a
#     comment saying it made hmmsearch faster. That is not true.
#
#       hmmpress  builds a compressed, indexed database of MANY profiles.
#       hmmscan   needs that index  (one sequence  vs  many profiles).
#       hmmsearch does NOT          (one profile   vs  many sequences).
#
#     We run hmmsearch, which reads the plain .hmm file directly. The hmmpress
#     step did nothing except create .h3m/.h3i/.h3f/.h3p files nobody read.
#     Harmless, but confusing -- so it is gone.
#
#     (Old .h3* files may still be lying around in ${HMM_DIR}. Safe to delete.)


# =============================================================================
# SECTION 10: RUN hmmsearch
# =============================================================================
# hmmsearch = search ONE HMM profile against MANY protein sequences.
#
# --cut_ga    = use Pfam's curated "gathering threshold" -- the cutoff Pfam
#               curators chose by hand for THIS family. Much better than
#               guessing an E-value yourself.
# --domtblout = write a per-domain hit table we parse later
# --cpu       = threads (SLURM gives us 8)
# --noali     = skip printing alignments (saves disk)

run_hmmsearch() {
    # $1 = short label for output files (e.g. cyp450)
    # $2 = Pfam HMM filename (e.g. PF00067.hmm)
    local label="$1"
    local hmm_file="${HMM_DIR}/$2"
    local domtbl="${OUTPUT_DIR}/${label}.domtbl"

    echo ""
    echo "  [${label}] hmmsearch --cut_ga  (${2})"
    echo "    Output: ${domtbl}"

    hmmsearch \
        --cut_ga \
        --domtblout "${domtbl}" \
        --cpu "${SLURM_CPUS_PER_TASK:-4}" \
        --noali \
        "${hmm_file}" \
        "${PEP_FILE}" \
        > "${OUTPUT_DIR}/${label}.out"

    echo "    Hits: $(grep -vc '^#' "${domtbl}" || true) domain rows"
}

echo ""
echo "Running hmmsearch (this usually takes a few minutes)..."

run_hmmsearch "cyp450"   "PF00067.hmm"
run_hmmsearch "omt_a"    "PF00891.hmm"
run_hmmsearch "omt_b"    "PF08100.hmm"
run_hmmsearch "omt_ccoa" "PF01596.hmm"


# =============================================================================
# SECTION 11: PARSE .domtbl FILES -> ID LISTS (isoform-level AND gene-level)
# =============================================================================
# >>> CHANGE 3: *** THIS SECTION CONTAINED THE BUG ***
#
# The old comment here said:
#     "RSEM / PyDESeq2 use transcript IDs: TRINITY_DN123_c0_g1_i1"
# That is FALSE. RSEM ran with --transcript-to-gene-map, so the count matrix
# is GENE-level:  TRINITY_DN0_c0_g1  (no _i1).
#
# The code faithfully implemented the wrong belief: it stripped ".p1" and
# stopped. The resulting IDs could never match the DE table -- overlap was
# exactly 0. Nothing errored. The files looked fine. That is the dangerous
# kind of bug.
#
# THE SCIENTIFIC DECISION ENCODED BELOW:
#   A gene is called CYP/OMT if ANY of its isoforms carries the domain.
#   This is the standard convention for de novo transcriptomes. Three isoforms
#   of one gene hitting PF00067 is ONE gene, not three -- so collapsing them
#   is not cosmetic, it is the correct count.

# --- helper: .domtbl -> ISOFORM ids (strip only TransDecoder's .pN) ---------
parse_domtbl_to_isoform_ids() {
    # $1 = input .domtbl   $2 = output .txt
    grep -v '^#' "$1" \
        | awk '{print $1}' \
        | sed -E 's/\.p[0-9]+$//' \
        | sort -u \
        > "$2"
}

# --- helper: ISOFORM ids -> GENE ids (strip Trinity's _iN, then dedupe) -----
isoform_ids_to_gene_ids() {
    # $1 = input isoform id list   $2 = output gene id list
    # The '$' anchors the match to END of line, so an ID that happens to
    # contain "_i" somewhere in the middle is never mangled.
    sed -E 's/_i[0-9]+$//' "$1" \
        | grep -v '^[[:space:]]*$' \
        | sort -u \
        > "$2"
}

echo ""
echo "Parsing domain tables -> ID lists..."

# --------------------------- CYP450 ----------------------------------------
# One Pfam profile (PF00067), so this is straightforward.
parse_domtbl_to_isoform_ids "${OUTPUT_DIR}/cyp450.domtbl" "${CYP_OUT}"
isoform_ids_to_gene_ids     "${CYP_OUT}" "${CYP_GENE_OUT}"

# ----------------------------- OMT -----------------------------------------
# >>> CHANGE 4: the OMT definition used to be a plain union of all three
#     domains:  PF00891 OR PF08100 OR PF01596.
#
#     Problem: PF08100 is the DIMERISATION domain. It is STRUCTURAL, not
#     catalytic. In plant Class II OMTs (COMT / EOMT / CVOMT) it sits at the
#     N-terminus alongside the PF00891 catalytic domain -- the two together
#     are diagnostic. But PF08100 on its OWN does not make something an OMT.
#     A plain union therefore let non-OMTs into your candidate list.
#
#     New definition:
#         CORE OMT      = PF00891  OR  PF01596     (the catalytic domains)
#         PF08100-only  -> written to a SEPARATE file for you to inspect
#
#     Why separate rather than deleted: a PF08100-only hit can be a REAL OMT
#     whose ORF was truncated by TransDecoder so only the N-terminal half was
#     predicted. Those are worth a look, not an automatic discard. Deciding
#     what to do with them is your call, not the script's.

# Core OMT: catalytic domains only (PF00891 = omt_a, PF01596 = omt_ccoa)
OMT_TMP="${OUTPUT_DIR}/omt_core.tmp"
cat "${OUTPUT_DIR}/omt_a.domtbl" \
    "${OUTPUT_DIR}/omt_ccoa.domtbl" \
    | grep -v '^#' \
    | awk '{print $1}' \
    | sed -E 's/\.p[0-9]+$//' \
    | sort -u \
    > "${OMT_TMP}"
mv "${OMT_TMP}" "${OMT_OUT}"
isoform_ids_to_gene_ids "${OMT_OUT}" "${OMT_GENE_OUT}"

# Dimerisation domain (PF08100 = omt_b) -> gene ids, for inspection
DIMER_ISO_TMP="${OUTPUT_DIR}/omt_dimer_iso.tmp"
DIMER_GENE_TMP="${OUTPUT_DIR}/omt_dimer_gene.tmp"
parse_domtbl_to_isoform_ids "${OUTPUT_DIR}/omt_b.domtbl" "${DIMER_ISO_TMP}"
isoform_ids_to_gene_ids     "${DIMER_ISO_TMP}" "${DIMER_GENE_TMP}"

# Genes with PF08100 but WITHOUT any catalytic OMT domain = the suspicious set
comm -23 "${DIMER_GENE_TMP}" "${OMT_GENE_OUT}" > "${OMT_DIMER_ONLY_OUT}"

N_DIMER_TOTAL=$(wc -l < "${DIMER_GENE_TMP}" | tr -d ' ')
N_DIMER_ONLY=$(wc -l < "${OMT_DIMER_ONLY_OUT}" | tr -d ' ')
rm -f "${DIMER_ISO_TMP}" "${DIMER_GENE_TMP}"


# =============================================================================
# SECTION 12: PRINT SUMMARY
# =============================================================================
# >>> CHANGE 4b/1b: report BOTH levels, and report the PF08100 diagnostic.

N_CYP=$(wc -l < "${CYP_OUT}" | tr -d ' ')
N_OMT=$(wc -l < "${OMT_OUT}" | tr -d ' ')
N_CYP_GENE=$(wc -l < "${CYP_GENE_OUT}" | tr -d ' ')
N_OMT_GENE=$(wc -l < "${OMT_GENE_OUT}" | tr -d ' ')

echo ""
echo "================================================================"
echo "HMMER gene finder complete"
echo "================================================================"
echo ""
echo "ISOFORM-level (raw record of which isoform carried the domain):"
echo "  CYP450: ${N_CYP} isoforms  -> ${CYP_OUT}"
echo "  OMT:    ${N_OMT} isoforms  -> ${OMT_OUT}"
echo ""
echo "GENE-level  *** USE THESE for PyDESeq2 joins ***"
echo "  CYP450: ${N_CYP_GENE} genes  (collapsed from ${N_CYP} isoforms)"
echo "    -> ${CYP_GENE_OUT}"
echo "  OMT:    ${N_OMT_GENE} genes  (collapsed from ${N_OMT} isoforms)"
echo "    -> ${OMT_GENE_OUT}"
echo ""
echo "OMT domain definition used:"
echo "  CORE    = PF00891 (Methyltransf_2) OR PF01596 (Methyltransf_3)"
echo "  PF08100 (Dimerisation) is STRUCTURAL and is NOT counted on its own."
echo "    genes with PF08100 total          : ${N_DIMER_TOTAL}"
echo "    genes with PF08100 but no catalytic domain : ${N_DIMER_ONLY}"
echo "    -> ${OMT_DIMER_ONLY_OUT}"
echo "  (If that last number is large, inspect it -- some may be real OMTs"
echo "   whose ORF was truncated. If it is near zero, your old union was fine.)"
echo ""
echo "Raw domain tables (keep for later e-value filtering):"
echo "  ${OUTPUT_DIR}/cyp450.domtbl     (PF00067)"
echo "  ${OUTPUT_DIR}/omt_a.domtbl      (PF00891)"
echo "  ${OUTPUT_DIR}/omt_b.domtbl      (PF08100)"
echo "  ${OUTPUT_DIR}/omt_ccoa.domtbl   (PF01596)"
echo ""
echo "NEXT STEP -- check the join actually works:"
echo "  tail -n +2 ${ANALYSIS_DIR}/pydeseq2_${SPECIES_CODE_UPPER}_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv \\"
echo "    | cut -f1 | sort -u > /tmp/de_ids.txt"
echo "  comm -12 /tmp/de_ids.txt ${CYP_GENE_OUT} | wc -l"
echo "================================================================"
