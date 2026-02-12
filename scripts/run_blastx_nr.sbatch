#!/bin/bash
################################################################################
# Run BLASTX against NCBI nr (nucleotide query from .fasta CDS file)
################################################################################
#
# Use this for NUCLEOTIDE FASTA (.fasta) — e.g. all_genes_cds.fasta from
# extract_cds_from_gtf.py. BLAST translates your DNA in silico and compares to nr.
#
# For PROTEIN FASTA (.faa) use run_blastp_nr.sbatch instead.
#
# Options for nr:
#   - Local: set NR_DB to path to your nr database.
#   - Remote: leave NR_DB unset to use NCBI's nr via -remote.
#
# How to run:
#   sbatch scripts/run_blastx_nr.sbatch
#   sbatch scripts/run_blastx_nr.sbatch path/to/all_genes_cds.fasta
#   QUERY=06_analysis/blast_input_DC/all_genes_cds.fasta sbatch scripts/run_blastx_nr.sbatch
#
################################################################################

#SBATCH --job-name=blastx_nr
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/blastx_nr_%j.out
#SBATCH -e 06_analysis/blastx_nr_%j.err

set -e

BASE_DIR="${BASE_DIR:-/projects/tholl_lab_1/daisy_analysis}"
cd "${BASE_DIR}"

# Query: nucleotide .fasta (CDS). Default carrot CDS; override with first arg or QUERY env.
QUERY="${1:-${QUERY:-06_analysis/blast_input_DC/all_genes_cds.fasta}}"
if [ ! -f "$QUERY" ]; then
    echo "ERROR: Query file not found: $QUERY"
    echo "Usage: sbatch run_blastx_nr.sbatch [path/to/cds.fasta]"
    exit 1
fi

OUT_DIR=$(dirname "$QUERY")
OUT_BASE=$(basename "$QUERY" .fasta)
OUT_BASE=$(basename "$OUT_BASE" .fna)
BLAST_OUT="${OUT_DIR}/${OUT_BASE}_blast_nr.txt"

source ~/.bashrc
conda activate rnaseq
[ -z "${CONDA_DEFAULT_ENV}" ] && { echo "ERROR: conda env rnaseq not activated"; exit 1; }

echo "============================================================"
echo "BLASTX vs nr (nucleotide query → translated vs nr)"
echo "============================================================"
echo "Query:  $QUERY"
echo "Output: $BLAST_OUT"
echo "Started: $(date)"
echo ""

if [ -n "${NR_DB}" ] && [ -d "$(dirname "${NR_DB}")" ]; then
    echo "Using local nr database: $NR_DB"
    DB_ARG="-db ${NR_DB}"
else
    echo "Using NCBI nr via -remote."
    DB_ARG="-db nr -remote"
fi

blastx \
    -query "$QUERY" \
    $DB_ARG \
    -out "$BLAST_OUT" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    -evalue 1e-5 \
    -max_target_seqs 20 \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}"

echo ""
echo "Finished: $(date)"
echo "Output: $BLAST_OUT"
echo "Next: python scripts/combine_blast_deseq.py --blast $BLAST_OUT ..."
echo "============================================================"
