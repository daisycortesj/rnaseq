#!/bin/bash
################################################################################
# BUILD GENE-FAMILY DATABASE — GFF keyword search + HMMER domain scan + Combine
################################################################################
#
# CHANGED: 2026-03-31 — Renamed from run_cyp450_database.sbatch → run_gene_database.sbatch.
#   Now accepts FAMILY arg (cyp or omt) instead of being CYP-only.
#   Renamed Python script calls: cyp_gff_keyword_search.py → gene_gff_keyword_search.py,
#   cyp_hmmer_scan.py → gene_hmmer_scan.py, cyp_combine_results.py → gene_combine_results.py.
#
# Supports: CYP (cytochrome P450) and OMT (O-methyltransferase)
#
# What this does:
#   Step 1: Searches the GTF annotation for family-specific keywords
#   Step 2: Runs HMMER domain scan against all carrot proteins
#           CYP uses PF00067 (P450 domain), OMT uses PF00891 (O-methyltransferase)
#   Step 3: Combines both lists into a master gene-family database
#
# How to run (FAMILY is required: cyp or omt):
#   sbatch scripts/08_gene_families/run_gene_database.sbatch cyp
#   sbatch scripts/08_gene_families/run_gene_database.sbatch omt
#
# What you get (in 07_NRdatabase/<family>_database/):
#   <family>_gtf_candidates.csv   — genes found by keyword search
#   <family>_hmmer_confirmed.csv  — genes confirmed by domain scan
#   <family>_master_list.csv      — combined master list with evidence flags
#
# Next steps:
#   sbatch scripts/08_gene_families/run_gene_extract.sbatch cyp   (or omt)
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

#SBATCH --job-name=gene_db
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/gene_db_%j.out
#SBATCH -e 07_NRdatabase/gene_db_%j.err


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

# ===== Load shared config =====
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

set -e

source ~/.bashrc
conda activate rnaseq


# ══════════════════════════════════════════════════════════════════════════════
# FAMILY ARGUMENT (required: cyp or omt)
# ══════════════════════════════════════════════════════════════════════════════

FAMILY="${1:-}"

if [ -z "${FAMILY}" ]; then
    echo "ERROR: You need to specify a gene family!"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/08_gene_families/run_gene_database.sbatch cyp"
    echo "  sbatch scripts/08_gene_families/run_gene_database.sbatch omt"
    exit 1
fi

# Normalize to lowercase for paths, uppercase for labels
FAMILY=$(echo "${FAMILY}" | tr '[:upper:]' '[:lower:]')
FAMILY_UPPER=$(echo "${FAMILY}" | tr '[:lower:]' '[:upper:]')

if [ "${FAMILY}" != "cyp" ] && [ "${FAMILY}" != "omt" ]; then
    echo "ERROR: Unknown gene family: ${FAMILY}"
    echo "  Supported families: cyp, omt"
    exit 1
fi

# ══════════════════════════════════════════════════════════════════════════════
# FAMILY-SPECIFIC SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

if [ "${FAMILY}" = "cyp" ]; then
    OUTDIR="07_NRdatabase/cyp450_database"
    PFAM_ACC="PF00067"
    FAMILY_LABEL="CYP450 (cytochrome P450)"
elif [ "${FAMILY}" = "omt" ]; then
    OUTDIR="07_NRdatabase/omt_database"
    PFAM_ACC="PF00891"
    FAMILY_LABEL="OMT (O-methyltransferase)"
fi

GTF_CANDIDATES="${BASE_DIR}/${OUTDIR}/${FAMILY}_gtf_candidates.csv"
HMMER_CONFIRMED="${BASE_DIR}/${OUTDIR}/${FAMILY}_hmmer_confirmed.csv"
MASTER_LIST="${BASE_DIR}/${OUTDIR}/${FAMILY}_master_list.csv"
HMM_CACHE="${BASE_DIR}/${OUTDIR}/hmm_profiles"


# ══════════════════════════════════════════════════════════════════════════════
# SETUP
# ══════════════════════════════════════════════════════════════════════════════

cd "${BASE_DIR}"
mkdir -p "${OUTDIR}"

echo "========================================================================"
echo "  ${FAMILY_UPPER} DATABASE BUILD"
echo "  Family: ${FAMILY_LABEL}"
echo "  Pfam:   ${PFAM_ACC}"
echo "  Output: ${OUTDIR}/"
echo "  Started: $(date)"
echo "========================================================================"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: GFF KEYWORD SEARCH
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 1: GFF keyword search (${FAMILY_UPPER})..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/gene_gff_keyword_search.py" \
    --gff "${BASE_DIR}/04_reference/dc_genomic.gtf" \
    --output "${GTF_CANDIDATES}" \
    --family "${FAMILY}"

echo ""


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: HMMER DOMAIN SCAN
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 2: HMMER domain scan (${PFAM_ACC})..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/gene_hmmer_scan.py" \
    --proteins "${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa" \
    --output "${HMMER_CONFIRMED}" \
    --pfam "${PFAM_ACC}" \
    --family-name "${FAMILY_UPPER}" \
    --hmm-cache "${HMM_CACHE}" \
    --evalue 1e-5 \
    --cpus 4

echo ""


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: COMBINE RESULTS
# ══════════════════════════════════════════════════════════════════════════════

echo "STEP 3: Combining results..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/08_gene_families/gene_combine_results.py" \
    --gtf-csv "${GTF_CANDIDATES}" \
    --hmmer-csv "${HMMER_CONFIRMED}" \
    --gff "${BASE_DIR}/04_reference/dc_genomic.gtf" \
    --output "${MASTER_LIST}" \
    --family-name "${FAMILY_UPPER}"

echo ""
echo "========================================================================"
echo "  DONE! $(date)"
echo "  Family: ${FAMILY_LABEL}"
echo "  Results in: ${BASE_DIR}/${OUTDIR}/"
echo ""
echo "  NEXT STEPS:"
echo "    # Path B (short — straight to plots):"
echo "    sbatch scripts/08_gene_families/run_gene_extract.sbatch ${FAMILY}"
echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC ${FAMILY}_expressed"
echo ""
echo "    # Path C (full — with BLAST):"
echo "    sbatch scripts/08_gene_families/run_gene_extract.sbatch ${FAMILY}"
echo "    sbatch scripts/06_blast/run_blastp_discovery.sbatch DC swissprot \\"
echo "        ${BASE_DIR}/${OUTDIR}/${FAMILY}_proteins.fasta"
echo "    sbatch scripts/06_blast/run_combine_blast_deseq.sbatch DC swissprot discovery standard ${FAMILY}"
echo "    sbatch scripts/05_pydeseq2/run_step3_plots.sbatch DC ${FAMILY}_discovery"
echo "========================================================================"
