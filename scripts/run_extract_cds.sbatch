#!/bin/bash
################################################################################
# EXTRACT CDS SEQUENCES FROM GTF + GENOME
################################################################################
#
# PURPOSE:
#   Extract coding sequences (CDS) for all genes from GTF + genome FASTA.
#   Produces a nucleotide FASTA file for BLAST (blastx vs NR).
#
# PREREQUISITES (run these first):
#   1. featurecounts.sbatch DC
#      → gene_count_matrix.tsv + sample_metadata.tsv
#
#   2. run_pydeseq2_step1_analysis.sbatch DC
#      → pydeseq2_results_UNFILTERED.tsv + all_gene_ids.txt
#
# HOW TO RUN:
#   # Default (DC species):
#   sbatch scripts/run_extract_cds.sbatch DC
#
#   # Other species:
#   sbatch scripts/run_extract_cds.sbatch DG
#   sbatch scripts/run_extract_cds.sbatch MF
#
#   # Chain after DESeq2:
#   JOB1=$(sbatch --parsable scripts/run_pydeseq2_step1_analysis.sbatch DC)
#   sbatch --dependency=afterok:$JOB1 scripts/run_extract_cds.sbatch DC
#
# INPUT FILES (where they come from):
#
#   File                          Created by                                HPC path
#   ───────────────────────────── ───────────────────────────────────────── ─────────────────────────────────────
#   dc_genomic.gtf                NCBI reference download                  04_reference/
#   GCF_..._genomic.fna           NCBI reference download                  04_reference/
#   all_gene_ids.txt              run_pydeseq2_step1_analysis.sbatch DC    06_analysis/pydeseq2_DC_step1_unfiltered/
#
# OUTPUT:
#   06_analysis/blast_input_DC/all_genes_cds.fasta   — nucleotide CDS for BLAST
#
# NEXT STEPS:
#   # Translate to protein (for blastp):
#   sbatch scripts/run_translate_cds.sbatch DC
#
#   # Or run blastx directly (nucleotide vs protein DB):
#   sbatch scripts/blastx_discoveryfilter.sbatch DC swissprot
#
################################################################################

#SBATCH --job-name=extract_cds
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/extract_cds_%j.out
#SBATCH -e 06_analysis/extract_cds_%j.err

set -euo pipefail

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

SPECIES="${1:-DC}"
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# ── Activate environment ─────────────────────────────────────────────────────
source ~/.bashrc
conda activate rnaseq

echo "============================================================"
echo "Extract CDS Sequences — ${SPECIES}"
echo "============================================================"
echo ""

# ── Run with BASE_DIR set so paths resolve correctly ─────────────────────────
BASE_DIR="${BASE_DIR}" python3 "${CODE_DIR}/scripts/extract_cds_from_gtf.py" "${SPECIES}"

echo ""
echo "============================================================"
echo "COMPLETE — Extract CDS"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  # Translate to protein:"
echo "  sbatch 05_rnaseq-code/scripts/run_translate_cds.sbatch ${SPECIES}"
echo ""
echo "  # Or run blastx directly:"
echo "  sbatch 05_rnaseq-code/scripts/blastx_discoveryfilter.sbatch ${SPECIES} swissprot"
echo ""
echo "End time: $(date)"
