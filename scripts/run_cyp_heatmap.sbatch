#!/bin/bash
################################################################################
# CYP450 HEATMAP PIPELINE — Intersect + Extract + BLAST + Combine
################################################################################
#
# Prerequisite:
#   Run run_cyp450_database.sbatch FIRST to produce cyp_master_list.csv
#
# What this does (4 steps in one job):
#   1. Intersect CYP master list with PyDESeq2 results → cyp_expressed_list.tsv
#   2. Extract CYP protein sequences from the full FASTA → cyp_proteins.fasta
#   3. Run blastp on the CYP subset (~396 proteins → fast) → cyp_blast_results.tsv
#   4. Combine BLAST + expression (reuses combine_blast_deseq.py) → cyp_blast_annotated.tsv
#
# After this job finishes, run step 5a and 5b manually:
#   5a. Filter for DE genes:
#       sbatch scripts/run_pydeseq2_step2_filter.sbatch DC \
#           /projects/tholl_lab_1/daisy_analysis/07_NRdatabase/cyp450_database/cyp_blast_annotated.tsv
#
#   5b. Generate plots (heatmap, volcano, MA):
#       sbatch scripts/run_pydeseq2_step3_plots.sbatch DC \
#           /projects/tholl_lab_1/daisy_analysis/06_analysis/pydeseq2_DC_step2_filtered/cyp_blast_annotated_FILTERED.tsv
#
# How to run:
#   sbatch scripts/run_cyp_heatmap.sbatch
#
# Output (in 07_NRdatabase/cyp450_database/):
#   cyp_expressed_list.tsv    — CYPs in count matrix with expression stats
#   cyp_proteins.fasta        — subset protein sequences (gene_id headers)
#   cyp_blast_results.tsv     — raw BLAST hits (outfmt 6)
#   cyp_blast_annotated.tsv   — combined: expression + BLAST annotation
#
################################################################################

#SBATCH --job-name=cyp_heatmap
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 07_NRdatabase/cyp_heatmap_%j.out
#SBATCH -e 07_NRdatabase/cyp_heatmap_%j.err

set -euo pipefail

source ~/.bashrc
conda activate rnaseq

# ============================================================================
# PATHS — edit these if your layout differs
# ============================================================================
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"
OUTDIR="${BASE_DIR}/07_NRdatabase/cyp450_database"

CYP_LIST="${OUTDIR}/cyp_master_list.csv"
GFF="${BASE_DIR}/04_reference/dc_genomic.gtf"
PROTEIN_FASTA="${BASE_DIR}/04_reference/GCF_001625215.2_DH1_v3.0_protein.faa"
DESEQ="${BASE_DIR}/06_analysis/pydeseq2_DC_step1_unfiltered/pydeseq2_results_UNFILTERED.tsv"

EXPRESSED="${OUTDIR}/cyp_expressed_list.tsv"
CYP_FASTA="${OUTDIR}/cyp_proteins.fasta"
BLAST_OUT="${OUTDIR}/cyp_blast_results.tsv"
ANNOTATED="${OUTDIR}/cyp_blast_annotated.tsv"

BLAST_DB="swissprot"
BLAST_EVALUE="1e-4"
BLAST_THREADS="${SLURM_CPUS_PER_TASK:-8}"

export BLASTDB="${BASE_DIR}/07_NRdatabase/blastdb"

# ============================================================================
# VALIDATION
# ============================================================================
echo "========================================================================"
echo "  CYP450 HEATMAP PIPELINE"
echo "  Started: $(date)"
echo "========================================================================"

for f in "$CYP_LIST" "$PROTEIN_FASTA" "$DESEQ" "$GFF"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

mkdir -p "$OUTDIR"
cd "$BASE_DIR"

# ============================================================================
# STEP 1: Intersect CYP master list with PyDESeq2
# ============================================================================
echo ""
echo "STEP 1/4: Intersecting CYP master list with PyDESeq2..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/cyp_intersect_pydeseq2.py" \
    --cyp-list "$CYP_LIST" \
    --deseq "$DESEQ" \
    --gff "$GFF" \
    --output "$EXPRESSED"

echo ""

# ============================================================================
# STEP 2: Extract CYP protein sequences
# ============================================================================
echo "STEP 2/4: Extracting CYP protein sequences..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/cyp_extract_proteins.py" \
    --expressed "$EXPRESSED" \
    --fasta "$PROTEIN_FASTA" \
    --output "$CYP_FASTA"

N_PROTS=$(grep -c "^>" "$CYP_FASTA" || true)
echo "  Proteins in subset FASTA: $N_PROTS"
echo ""

# ============================================================================
# STEP 3: Run BLASTp on CYP proteins (discovery mode)
# ============================================================================
echo "STEP 3/4: Running BLASTp (${BLAST_DB}, ~${N_PROTS} proteins)..."
echo "------------------------------------------------------------------------"

which blastp
blastp -version

blastp \
    -query "$CYP_FASTA" \
    -db "$BLAST_DB" \
    -evalue "$BLAST_EVALUE" \
    -max_target_seqs 25 \
    -qcov_hsp_perc 40 \
    -max_hsps 1 \
    -num_threads "$BLAST_THREADS" \
    -outfmt "6 qseqid sseqid stitle sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp" \
    -out "$BLAST_OUT"

BLAST_LINES=$(wc -l < "$BLAST_OUT" || echo 0)
BLAST_QUERIES=$(cut -f1 "$BLAST_OUT" | sort -u | wc -l || echo 0)
echo "  BLAST hits: $BLAST_LINES lines, $BLAST_QUERIES unique queries"
echo ""

# ============================================================================
# STEP 4: Combine BLAST + expression (reuse existing script)
# ============================================================================
echo "STEP 4/4: Combining BLAST annotations with expression data..."
echo "------------------------------------------------------------------------"

python "${CODE_DIR}/scripts/combine_blast_deseq.py" \
    --deseq "$EXPRESSED" \
    --blast "$BLAST_OUT" \
    --output "$ANNOTATED"

echo ""
echo "========================================================================"
echo "  STEPS 1-4 DONE! $(date)"
echo "  Results in: ${OUTDIR}/"
echo ""
echo "  Key output files:"
echo "    ${EXPRESSED}     — CYP genes with expression stats"
echo "    ${CYP_FASTA}     — subset protein FASTA"
echo "    ${BLAST_OUT}     — raw BLAST hits"
echo "    ${ANNOTATED}     — expression + BLAST combined"
echo ""
echo "  NEXT STEPS (run these after this job finishes):"
echo ""
echo "  Step 5a — Filter for DE genes:"
echo "    sbatch scripts/run_pydeseq2_step2_filter.sbatch DC ${ANNOTATED}"
echo ""
echo "  Step 5b — Generate plots (heatmap, volcano, MA):"
echo "    sbatch scripts/run_pydeseq2_step3_plots.sbatch DC <path to FILTERED output from 5a>"
echo "========================================================================"
