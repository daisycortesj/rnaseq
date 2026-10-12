#!/bin/bash
# =====================================================================
# SCRIPT:  run_check_mdb_site.sbatch
# PURPOSE: Ask whether P450 proteins in your nutmeg assembly have the
#          methylenedioxy-bridge (MDB) site reported by the Ibdah group.
#
# WHAT THE SITE IS:
#   Sesamin synthase (SS) builds the methylenedioxy bridge. The residue
#   that matters is Ala308, inside the I-helix motif G-T-D-[A]-S.
#   The nutmeg candidates Mf30795, Mf31663, and Mf41180 have Ser there
#   instead of Ala. Both A and S are called MDB-like.
#   A Thr (T) at that spot is the ordinary P450 residue. In SS, changing
#   Ala308 to Thr reduced MDB activity.
#
# WHAT GETS CHECKED:
#   Every protein in the TransDecoder .pep file that looks like a P450
#   (it has the heme-binding motif) and is at least 200 amino acids.
#   Each one is aligned to SS. The table reports the residue that lines
#   up with Ala308.
#
# HOW TO SUBMIT (on the cluster, from the rnaseq repo):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   conda activate rnaseq
#   conda install -c conda-forge biopython     # only the first time
#   sbatch scripts/11_verify/run_check_mdb_site.sbatch MF
#
# OPTIONAL:
#   sbatch scripts/11_verify/run_check_mdb_site.sbatch MF nofungi
#   IDS=/path/to/ids.txt sbatch scripts/11_verify/run_check_mdb_site.sbatch MF
#   The ids file is one protein ID per line. Only those proteins are checked.
#
# INPUT:
#   01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta.transdecoder.pep
#
# OUTPUT:
#   06_analysis/mdb_site_MF/mdb_site.tsv
#   06_analysis/mdb_site_MF/mdb_site_SUMMARY.txt
#
# HOW TO READ mdb_site.tsv:
#   class "MDB-like (A/S)"  -> this protein has the Ibdah MDB residue
#   class "typical P450 (T)" -> ordinary P450, not the MDB site
#   class "region missing"   -> the protein is a partial ORF
#   class "other"            -> look at the I-helix window yourself
# =====================================================================


# ===== SLURM directives =====
#SBATCH --job-name=mdb_site
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/mdb_site_%j.out
#SBATCH -e 06_analysis/mdb_site_%j.err

set -euo pipefail


# =====================================================================
# LOAD CONFIG
# =====================================================================
# SLURM copies this .sbatch file into a spool folder, so the script's
# own path is not the repo. Search the folder you submitted from, then
# the known repo path.

CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not find scripts/config.sh"
    echo "Submit from: /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    exit 1
fi

source "${CONFIG_FILE}"
REPO_ROOT="$(cd "$(dirname "${CONFIG_FILE}")/.." && pwd)"
CHECK_PY="${REPO_ROOT}/scripts/11_verify/check_mdb_site.py"


# =====================================================================
# SPECIES AND WHICH ASSEMBLY
# =====================================================================
# First argument is the species code. Second argument "nofungi" reads
# the fungal-filtered transcripts instead of the original CD-HIT file.

SPECIES_CODE="${1:-MF}"
ASSEMBLY_CHOICE="${2:-}"

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')
apply_assembly_adjustment "${ASSEMBLY_CHOICE}"

# TransDecoder writes the protein file next to the transcript FASTA.
PEP_FILE="${ACTIVE_ASSEMBLY}.transdecoder.pep"
OUT_DIR="${ANALYSIS_DIR}/mdb_site_${SPECIES_CODE}${ASSEMBLY_TAG}"
OUT_TSV="${OUT_DIR}/mdb_site.tsv"
SUMMARY="${OUT_DIR}/mdb_site_SUMMARY.txt"


# =====================================================================
# SAFETY CHECKS
# =====================================================================

if [[ ! -f "${CHECK_PY}" ]]; then
    echo "ERROR: Checker not found: ${CHECK_PY}"
    echo "Sync the repo so scripts/11_verify/check_mdb_site.py is there."
    exit 1
fi

if [[ ! -f "${PEP_FILE}" ]]; then
    echo "ERROR: TransDecoder proteins not found: ${PEP_FILE}"
    echo "Run this first:"
    echo "  sbatch scripts/03_assembly/run_transdecoder.sbatch ${SPECIES_CODE} ${ASSEMBLY_CHOICE}"
    exit 1
fi

if [[ -n "${IDS:-}" && ! -f "${IDS}" ]]; then
    echo "ERROR: ID list not found: ${IDS}"
    exit 1
fi

mkdir -p "${OUT_DIR}"

# set +u lets conda start. set -u turns the strict check back on.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq
set -u

if ! python3 -c "import Bio" 2>/dev/null; then
    echo "ERROR: Biopython is not installed in the rnaseq environment."
    echo "On a login node, run:"
    echo "  conda activate rnaseq"
    echo "  conda install -c conda-forge biopython"
    exit 1
fi

echo "============================================================"
echo "MDB site check — ${SPECIES_NAME} (${SPECIES_CODE})"
echo "============================================================"
echo "Host: $(hostname)  Start: $(date)"
echo "Proteins: ${PEP_FILE}"
echo "Output:   ${OUT_TSV}"
if [[ -n "${IDS:-}" ]]; then
    echo "ID list:  ${IDS}"
fi
echo "============================================================"


# =====================================================================
# RUN THE CHECK
# =====================================================================

ID_ARGS=()
if [[ -n "${IDS:-}" ]]; then
    ID_ARGS=(--ids "${IDS}")
fi

python3 "${CHECK_PY}" \
    --pep "${PEP_FILE}" \
    --out "${OUT_TSV}" \
    "${ID_ARGS[@]}"


# =====================================================================
# SHORT SUMMARY
# =====================================================================
# Column 5 is the class. This counts how many proteins fell in each class
# and copies the MDB-like rows into a short file you can read first.

{
    echo "MDB site check"
    echo "Proteins: ${PEP_FILE}"
    echo "Full table: ${OUT_TSV}"
    echo ""
    echo "Counts by class"
    awk -F '\t' 'NR > 1 { count[$5] = count[$5] + 1 } END { for (name in count) print count[name] "\t" name }' "${OUT_TSV}"
    echo ""
    echo "Proteins with an MDB-like residue (A or S) at the SS Ala308 position"
    awk -F '\t' 'NR == 1 || $5 ~ /MDB-like/' "${OUT_TSV}"
} | tee "${SUMMARY}"

echo "Finished: $(date)"
echo "Read this first: ${SUMMARY}"
