#!/bin/bash
################################################################################
# FEATURECOUNTS — Count mapped reads per gene + build count matrix
################################################################################
#
# PURPOSE:
#   This script does TWO things in one job:
#     1. Runs featureCounts (Subread package) on samtools-sorted BAMs
#        to count how many read pairs overlap each gene
#     2. Calls build_count_matrix.py to convert the featureCounts output
#        into PyDESeq2-ready files (gene_count_matrix.tsv + sample_metadata.tsv)
#
#
# WHAT featureCounts DOES:
#   featureCounts reads each BAM file and a GTF annotation file. For every
#   read pair in the BAM, it checks which gene's exons the read overlaps.
#   If a read overlaps multiple genes, --largestOverlap assigns it to the
#   gene with the biggest overlap. The result is an integer count per gene
#   per sample — exactly what differential expression tools need.
#
# WHAT build_count_matrix.py DOES (the reformatting):
#   featureCounts output has extra columns (Chr, Start, End, Strand, Length)
#   and uses full BAM file paths as column headers. build_count_matrix.py:
#     - Strips the annotation columns
#     - Converts BAM paths to clean sample names (DC1L1, DC1R1, etc.)
#     - Extracts biological metadata (species, condition, replicate) from
#       the sample names using a regex
#     - Writes gene_count_matrix.tsv + sample_metadata.tsv + count_summary.txt
#
# PIPELINE ORDER:
#   Step 1: STAR alignment                → *_Aligned.sortedByCoord.out.bam
#   Step 2: samtools_sort.sbatch          → *_Aligned.sortedBySamtools.bam
#   Step 3: featureCounts.sbatch (THIS)   → featurecounts.txt
#                                          → gene_count_matrix.tsv  (PyDESeq2 input)
#                                          → sample_metadata.tsv    (PyDESeq2 design)
#
# HOW TO RUN:
#   sbatch scripts/featureCounts.sbatch DC    (for D. carota)
#   sbatch scripts/featureCounts.sbatch DG    (for D. glaber)
#
# NOTE:
#   featureCounts requires a GTF annotation file. MF (nutmeg) does not
#   have a GTF, so this script only works for DC and DG (carrot genome).
#   For MF, use the STAR GeneCounts path instead:
#     sbatch scripts/build_count_matrix.sbatch MF --type star
#
# CUSTOM GTF PATH:
#   GTF=/path/to/custom.gtf sbatch scripts/featureCounts.sbatch DC
#
# INPUT:
#   03_count_tables/{species}/*_Aligned.sortedBySamtools.bam  (from samtools_sort.sbatch)
#   04_reference/dc_genomic.gtf                         (gene annotation)
#
# OUTPUT:
#   03_count_tables/{species}/featurecounts.txt          (raw featureCounts output)
#   03_count_tables/{species}/featurecounts.txt.summary  (assignment stats)
#   03_count_tables/{species}/gene_count_matrix.tsv      (clean matrix for PyDESeq2)
#   03_count_tables/{species}/sample_metadata.tsv        (sample info for PyDESeq2)
#   03_count_tables/{species}/count_summary.txt          (human-readable stats)
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════
# --nodes=1         : featurecounts is single-node (multithreaded via -T)
# --cpus-per-task=16 : featureCounts uses threads (-T) to parallelize counting
# --mem=16G         : enough for loading all BAMs + GTF into memory

#SBATCH --job-name=featurecounts
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

# Working directory and log files
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/03_count_tables
#SBATCH -o featurecounts_%j.out
#SBATCH -e featurecounts_%j.err

set -euo pipefail


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════
# BASE_DIR      : root of the project on the HPC
# BAM_DIR       : where sorted BAMs live (one subdir per species)
# COUNT_DIR     : where count tables and matrices are saved
# REFERENCE_DIR : where genome files and GTF annotations live
# THREADS       : CPU threads for featureCounts; inherited from SLURM

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
BAM_DIR="${BASE_DIR}/02_mapped"
COUNT_DIR="${BASE_DIR}/03_count_tables"
REFERENCE_DIR="${BASE_DIR}/04_reference"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════
# $1 is the first argument: sbatch scripts/featureCounts.sbatch DC  ← $1 = "DC"

SPECIES="${1:-}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: No species code provided"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/featureCounts.sbatch DC"
    echo "  sbatch scripts/featureCounts.sbatch DG"
    echo ""
    echo "NOTE: MF (nutmeg) has no GTF — use the STAR GeneCounts path instead."
    exit 1
fi

# Convert to uppercase so 'dc' and 'DC' both work
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Map species code → directory name, display name, and GTF path
# DC and DG both use the carrot genome (dc_genomic.gtf)
# MF has no GTF, so featureCounts cannot be used
case "${SPECIES}" in
    DC)
        SPECIES_DIR="00_1_DC"
        SPECIES_NAME="D. carota ssp. maximus"
        DEFAULT_GTF="${REFERENCE_DIR}/dc_genomic.gtf"
        ;;
    DG)
        SPECIES_DIR="00_2_DG"
        SPECIES_NAME="Daucus glaber"
        DEFAULT_GTF="${REFERENCE_DIR}/dc_genomic.gtf"   # DG aligns to carrot genome too
        ;;
    MF)
        echo "ERROR: MF (nutmeg) has no GTF annotation file."
        echo "featureCounts requires a GTF. Use the STAR GeneCounts path instead:"
        echo "  sbatch scripts/build_count_matrix.sbatch MF"
        exit 1
        ;;
    *)
        echo "ERROR: Unknown species code '${SPECIES}'"
        echo "Valid options: DC, DG"
        exit 1
        ;;
esac

# Allow overriding the GTF path via environment variable:
#   GTF=/path/to/other.gtf sbatch scripts/featureCounts.sbatch DC
GTF="${GTF:-${DEFAULT_GTF}}"

# Input: sorted BAMs from samtools_sort.sbatch
INPUT_DIR="${BAM_DIR}/${SPECIES_DIR}"

# Output: count tables go here (same dir where STAR ReadsPerGene files live)
OUTPUT_DIR="${COUNT_DIR}/${SPECIES_DIR}"
OUTPUT_FILE="${OUTPUT_DIR}/featurecounts.txt"


# ══════════════════════════════════════════════════════════════════════════════
# ACTIVATE CONDA ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════
# featureCounts (from the Subread package) and Python are in the rnaseq env

source ~/.bashrc
conda activate rnaseq

# Install Subread (contains featureCounts) if not already available
if ! command -v featureCounts &> /dev/null; then
    echo "featureCounts not found — installing Subread into rnaseq env..."
    conda install -y -c bioconda subread
fi

echo "============================================================"
echo "featureCounts — ${SPECIES_NAME} (${SPECIES})"
echo "============================================================"
echo "Start time:  $(date)"
echo "BAM dir:     ${INPUT_DIR}"
echo "GTF:         ${GTF}"
echo "Output:      ${OUTPUT_FILE}"
echo "Threads:     ${THREADS}"
echo ""

# Verify featureCounts is available (via rna_pipeline wrapper)
python3 -c "from rna_pipeline.tools.featurecounts import check_featurecounts; assert check_featurecounts(), 'featureCounts not found — run: conda install -c bioconda subread'"
echo ""


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATE INPUT FILES
# ══════════════════════════════════════════════════════════════════════════════
# Check that the GTF and sorted BAMs exist before running featureCounts.

# 1. Check GTF exists
if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF file not found: ${GTF}"
    echo ""
    echo "Expected location: ${DEFAULT_GTF}"
    echo "Override with: GTF=/path/to/your.gtf sbatch scripts/featureCounts.sbatch ${SPECIES}"
    exit 1
fi

# 2. Check BAM directory exists
if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: BAM directory not found: ${INPUT_DIR}"
    echo "Did you run the alignment and samtools sort steps first?"
    exit 1
fi

# 3. Collect all samtools-sorted BAMs into an array
BAM_FILES=( "${INPUT_DIR}"/*_Aligned.sortedBySamtools.bam )

# 4. Check that the glob matched real files
if [ ${#BAM_FILES[@]} -eq 0 ] || [ ! -f "${BAM_FILES[0]}" ]; then
    echo "ERROR: No *_Aligned.sortedBySamtools.bam files found in ${INPUT_DIR}"
    echo ""
    echo "Run samtools sort first:"
    echo "  sbatch scripts/samtools_sort.sbatch ${SPECIES}"
    exit 1
fi

# Print which BAMs will be counted
echo "Found ${#BAM_FILES[@]} sorted BAM file(s):"
for b in "${BAM_FILES[@]}"; do
    echo "  - $(basename "$b")"
done
echo ""

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: RUN FEATURECOUNTS
# ══════════════════════════════════════════════════════════════════════════════
#
# featureCounts reads the GTF annotation and all BAMs at once, producing
# a single combined output file with one column per BAM.
#
# FLAG REFERENCE (matching previous student's approach):
#   -a <GTF>          : annotation file — tells featureCounts where genes/exons are
#   -o <file>         : output file path
#   --largestOverlap  : when a read overlaps multiple genes, assign it to the
#                        gene with the largest overlap (default: skip ambiguous)
#   -t exon           : count reads overlapping "exon" features in the GTF
#                        (standard for RNA-seq — exons are the expressed parts)
#   -g gene_id        : group exon-level counts by the "gene_id" attribute,
#                        so the output is gene-level (not exon-level)
#   -p                : input data is paired-end
#   --countReadPairs  : count read PAIRS (fragments) instead of individual reads
#                        (recommended for paired-end — avoids double-counting)
#   -T <threads>      : number of CPU threads for parallel processing
#
# The last argument "${BAM_FILES[@]}" expands to all BAM file paths.
# featureCounts processes them all in one run and produces one output file
# with a column for each BAM.

echo "Running featureCounts..."
echo "------------------------------------------------------------"

featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_FILE}" \
    --largestOverlap \
    -t exon \
    -g gene_id \
    -p --countReadPairs \
    -T "${THREADS}" \
    "${BAM_FILES[@]}"


# ══════════════════════════════════════════════════════════════════════════════
# FEATURECOUNTS OUTPUT CHECK
# ══════════════════════════════════════════════════════════════════════════════
# featureCounts produces two files:
#   featurecounts.txt          — the main count table
#   featurecounts.txt.summary  — how many reads were assigned/unassigned and why

echo ""
echo "============================================================"
echo "DONE — featureCounts"
echo "============================================================"

if [ -f "${OUTPUT_FILE}" ]; then
    # Line count minus 2 header lines = number of genes
    N_GENES=$(tail -n +3 "${OUTPUT_FILE}" | wc -l)
    N_SAMPLES=${#BAM_FILES[@]}
    echo "Genes counted:   ${N_GENES}"
    echo "Samples:         ${N_SAMPLES}"
    echo "Output:          ${OUTPUT_FILE}"
    echo "Summary:         ${OUTPUT_FILE}.summary"
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_FILE}" "${OUTPUT_FILE}.summary" 2>/dev/null
else
    echo "ERROR: featureCounts output file not created"
    exit 1
fi


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: BUILD COUNT MATRIX + SAMPLE METADATA
# ══════════════════════════════════════════════════════════════════════════════
#
# featureCounts output (featurecounts.txt) is NOT directly usable by PyDESeq2:
#   - It has extra annotation columns: Geneid, Chr, Start, End, Strand, Length
#   - Column headers are full BAM file paths, not clean sample names
#   - There's no sample_metadata.tsv telling PyDESeq2 what to compare
#
# build_count_matrix.py (--type featurecounts) fixes all of this:
#   1. Reads featurecounts.txt
#   2. Drops the annotation columns (Chr, Start, End, Strand, Length)
#   3. Renames BAM paths → sample names:
#        /path/DC1L1_Aligned.sortedBySamtools.bam → DC1L1
#   4. Parses each sample name with a regex to extract metadata:
#        DC1L1 → group=DC, group_number=1, condition=L, replicate=1
#   5. Writes three files:
#        gene_count_matrix.tsv  — clean genes × samples integer counts
#        sample_metadata.tsv    — one row per sample with biological info
#        count_summary.txt      — quick stats (total reads, zero-count genes, etc.)

echo ""
echo "Building count matrix + sample metadata from featureCounts output..."
echo "------------------------------------------------------------"

python3 "${CODE_DIR}/scripts/04_counting/build_count_matrix.py" \
    "${OUTPUT_DIR}" \
    -o "${OUTPUT_DIR}" \
    --type featurecounts

if [ $? -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "COMPLETE — featureCounts + count matrix + metadata"
    echo "============================================================"
    echo ""
    echo "Output files in ${OUTPUT_DIR}/:"
    echo "  featurecounts.txt        — raw featureCounts output"
    echo "  gene_count_matrix.tsv    — clean count matrix (PyDESeq2 input)"
    echo "  sample_metadata.tsv      — sample info (PyDESeq2 design)"
    echo "  count_summary.txt        — stats summary"
    echo ""
    echo "Samples:"
    head -1 "${OUTPUT_DIR}/gene_count_matrix.tsv" | cut -f2- | tr '\t' '\n' | sed 's/^/  - /'
    echo ""
    echo "Next step — run PyDESeq2:"
    echo "  sbatch scripts/run_pydeseq2_step1_analysis.sbatch ${SPECIES}"
else
    echo "WARNING: featureCounts succeeded but count matrix creation failed"
    echo "You can retry with:"
    echo "  sbatch scripts/build_count_matrix.sbatch ${SPECIES} --type featurecounts"
fi

echo ""
echo "End time: $(date)"
