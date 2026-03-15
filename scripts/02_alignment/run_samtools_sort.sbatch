#!/bin/bash
################################################################################
# SAMTOOLS SORT — Re-sort STAR BAMs with samtools
################################################################################
#
# PURPOSE:
#   Takes STAR-produced BAMs (already sorted by coordinate during alignment)
#   and re-sorts them with samtools.
#
# WHY THIS STEP EXISTS:
#   STAR sorts BAMs during alignment (--outSAMtype BAM SortedByCoordinate),
#   but its sort header can differ slightly from the samtools standard.
#   Re-sorting with samtools ensures maximum compatibility with downstream
#   tools like featureCounts, GATK, or IGV. This also creates a .bai index
#   file which is needed for visualization in genome browsers.
#
# PIPELINE ORDER:
#   Step 1: STAR alignment               → *_Aligned.sortedByCoord.out.bam
#   Step 2: samtools sort (THIS SCRIPT)   → *_Aligned.sortedBySamtools.bam + .bai
#   Step 3: featurecounts.sbatch          → featurecounts.txt → gene_count_matrix.tsv
#
# HOW TO RUN:
#   sbatch scripts/samtools_sort.sbatch DC    (for D. carota)
#   sbatch scripts/samtools_sort.sbatch DG    (for D. glaber)
#   sbatch scripts/samtools_sort.sbatch MF    (for M. fragrans)
#
# INPUT:
#   02_mapped/{species}/*_Aligned.sortedByCoord.out.bam
#   (produced by STAR alignment in run_alignment_all.sbatch)
#
# OUTPUT:
#   02_mapped/{species}/*_Aligned.sortedBySamtools.bam      (re-sorted BAM)
#   02_mapped/{species}/*_Aligned.sortedBySamtools.bam.bai  (BAM index)
#
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# SLURM SETTINGS
# ══════════════════════════════════════════════════════════════════════════════
# --nodes=1        : samtools sort is single-node (multithreaded, not MPI)
# --cpus-per-task=8: samtools uses threads (-@) for parallel sorting
# --mem=16G        : sorting large BAMs needs memory for buffering chunks

#SBATCH --job-name=samtools_sort
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu

# Working directory and log files
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/03_count_tables
#SBATCH -o samtools_sort_%j.out
#SBATCH -e samtools_sort_%j.err

# Exit immediately if any command fails (-e), treat unset variables as
# errors (-u), and fail pipelines if any stage fails (-o pipefail)
set -euo pipefail


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════
# BASE_DIR : root of the project on the HPC
# BAM_DIR  : where STAR wrote the aligned BAM files (one subdir per species)
# THREADS  : how many CPU threads samtools can use; inherited from SLURM

# ===== Load shared config =====
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

BAM_DIR="${MAPPED_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


# ══════════════════════════════════════════════════════════════════════════════
# GET SPECIES FROM COMMAND LINE
# ══════════════════════════════════════════════════════════════════════════════
# $1 is the first argument you pass after the script name, e.g.:
#   sbatch scripts/samtools_sort.sbatch DC   ← $1 = "DC"
# ${1:-} means "use $1, or empty string if nothing was provided"

SPECIES="${1:-}"

if [ -z "${SPECIES}" ]; then
    echo "ERROR: No species code provided"
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/samtools_sort.sbatch DC"
    echo "  sbatch scripts/samtools_sort.sbatch DG"
    echo "  sbatch scripts/samtools_sort.sbatch MF"
    exit 1
fi

# Convert to uppercase so 'dc' and 'DC' both work
SPECIES=$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')

# Map species code → directory name and display name (uses shared config)
get_sample_info "${SPECIES}" || exit 1

# Full path to the directory containing this species' BAM files
INPUT_DIR="${BAM_DIR}/${SPECIES_DIR}"


# ══════════════════════════════════════════════════════════════════════════════
# ACTIVATE CONDA ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════
# We need samtools from the 'rnaseq' conda environment.
# source ~/.bashrc loads conda into the shell, then we activate.

source ~/.bashrc
conda activate rnaseq

echo "============================================================"
echo "Samtools Sort — ${SPECIES_NAME} (${SPECIES})"
echo "============================================================"
echo "Start time:  $(date)"
echo "Input dir:   ${INPUT_DIR}"
echo "Threads:     ${THREADS}"
echo ""

# Print samtools version for reproducibility
samtools --version | head -1


# ══════════════════════════════════════════════════════════════════════════════
# FIND BAM FILES TO SORT
# ══════════════════════════════════════════════════════════════════════════════
# STAR produces files named: {sample}_Aligned.sortedByCoord.out.bam
# We glob for all of them in the species directory.

if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: BAM directory not found: ${INPUT_DIR}"
    echo "Did you run the alignment step first?"
    echo "  sbatch scripts/run_alignment_all.sbatch"
    exit 1
fi

# Bash array: collect all STAR BAMs into BAM_FILES
BAM_FILES=( "${INPUT_DIR}"/*_Aligned.sortedByCoord.out.bam )

# Check that the glob actually matched real files
if [ ${#BAM_FILES[@]} -eq 0 ] || [ ! -f "${BAM_FILES[0]}" ]; then
    echo "ERROR: No *_Aligned.sortedByCoord.out.bam files found in ${INPUT_DIR}"
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM file(s) to sort"
echo "------------------------------------------------------------"

# Counters for the summary at the end
SORTED_COUNT=0
FAILED_COUNT=0


# ══════════════════════════════════════════════════════════════════════════════
# SORT EACH BAM FILE
# ══════════════════════════════════════════════════════════════════════════════
# For each STAR BAM:
#   1. Extract sample name from filename (strip suffix)
#   2. Define output path with _sortedBySamtools suffix
#   3. Skip if output already exists (idempotent — safe to re-run)
#   4. samtools sort: re-sort by coordinate, using $THREADS parallel threads
#   5. samtools index: create .bai index for the sorted BAM

for bam in "${BAM_FILES[@]}"; do

    # Strip the suffix to get just the sample name
    # e.g. DC1L1_Aligned.sortedByCoord.out.bam → DC1L1
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)

    # Output BAM goes in the same directory with a different suffix
    outbam="${INPUT_DIR}/${sample}_Aligned.sortedBySamtools.bam"

    echo ""
    echo "▶ ${sample}"
    echo "  Input:  $(basename "${bam}")"
    echo "  Output: $(basename "${outbam}")"

    # Skip if already sorted (delete the output file to force re-sort)
    if [ -f "${outbam}" ]; then
        echo "  ⏭ Already exists, skipping (delete to re-sort)"
        SORTED_COUNT=$((SORTED_COUNT + 1))
        continue
    fi

    # samtools sort flags:
    #   -@  : number of sorting/compression threads
    #   -o  : output file path
    # The input BAM is the last positional argument
    if samtools sort -@ "${THREADS}" -o "${outbam}" "${bam}"; then
        echo "  ✓ Sorted successfully"

        # samtools index creates a .bai file alongside the BAM
        # This is needed for IGV, GATK, and quick random access
        samtools index -@ "${THREADS}" "${outbam}"
        echo "  ✓ Indexed (.bai created)"

        SORTED_COUNT=$((SORTED_COUNT + 1))
    else
        echo "  ✗ FAILED"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
done


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "============================================================"
echo "DONE — samtools sort"
echo "============================================================"
echo "Sorted:  ${SORTED_COUNT}"
echo "Failed:  ${FAILED_COUNT}"
echo "End time: $(date)"
echo ""
echo "Sorted BAMs:"
ls -lh "${INPUT_DIR}"/*_sortedBySamtools.bam 2>/dev/null || echo "  (none)"
echo ""
echo "Next step — run featureCounts:"
echo "  sbatch scripts/featurecounts.sbatch ${SPECIES}"
