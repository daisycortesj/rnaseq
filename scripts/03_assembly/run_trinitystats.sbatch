#!/bin/bash
# =============================================================================
# SCRIPT:  run_trinitystats.sbatch
# WHAT IT DOES:
#   Automatically finds every completed Trinity assembly in 01_processed/
#   and runs TrinityStats.pl on each one to report assembly quality metrics.
#
# HOW TO SUBMIT:
#   All species:    sbatch scripts/03_assembly/run_trinitystats.sbatch
#   One species:    sbatch scripts/03_assembly/run_trinitystats.sbatch MF
#                   sbatch scripts/03_assembly/run_trinitystats.sbatch DG
#
# FASTA FILE LOCATION:
#   The FASTA is inside the assembly folder, named after the sample:
#     01_processed/MFF1_trinity_assembly/MFF1_trinity_assembly.Trinity.fasta
#
# OUTPUTS:
#   Per-sample:  01_processed/MFF1_trinity_assembly/Trinity_assembly_stats.txt
#   Summary:     01_processed/trinitystats_summary.txt
# =============================================================================


# ---- SLURM cluster resource requests ----
#SBATCH --job-name=trinitystats
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinitystats_%j.out
#SBATCH -e trinitystats_%j.err

# Stop the script immediately if any command fails
set -euo pipefail

# Load shared folder paths (gives us PROCESSED_DIR and CONDA_ENV)
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# Optional species filter — pass on the command line, no editing needed
# Example: sbatch run_trinitystats.sbatch MF   → only process MF* samples
#          sbatch run_trinitystats.sbatch       → process everything
SPECIES_FILTER="${1:-}"


# ---- Activate conda so TrinityStats.pl is available ----
source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

# Check that TrinityStats.pl is actually installed before starting
if ! command -v TrinityStats.pl &> /dev/null; then
    echo "ERROR: TrinityStats.pl not found in your conda environment."
    echo "Fix: conda install -c bioconda trinity"
    exit 1
fi


# ---- Print a job header to the log ----
echo "========================================================"
echo "TrinityStats — Assembly Quality Metrics"
echo "Started:  $(date)"
echo "Species:  ${SPECIES_FILTER:-ALL}"
echo "========================================================"


# ---- Create the summary table file ----
# This collects one row per sample so you can compare everything at once.
SUMMARY_FILE="${PROCESSED_DIR}/trinitystats_summary.txt"
printf "%-20s %14s %14s %10s %10s\n" \
    "Sample" "Transcripts" "Genes" "N50" "GC%" > "${SUMMARY_FILE}"
printf "%-20s %14s %14s %10s %10s\n" \
    "--------------------" "--------------" "--------------" "----------" "----------" \
    >> "${SUMMARY_FILE}"


# ---- Counters ----
count_total=0
count_done=0
count_skipped=0
count_failed=0


# ============================================================================
# MAIN LOOP
# Searches 01_processed/ for every folder ending in "_trinity_assembly"
# and runs TrinityStats on each one.
# ============================================================================
for assembly_folder in "${PROCESSED_DIR}"/*_trinity_assembly; do

    # If the glob found nothing, the variable contains a literal string — skip it
    if [[ ! -d "${assembly_folder}" ]]; then
        continue
    fi

    # Get the sample name (strip the folder path and the "_trinity_assembly" suffix)
    # Example: /path/01_processed/MFF1_trinity_assembly  →  MFF1
    sample=$(basename "${assembly_folder}" "_trinity_assembly")

    # If a species filter was given, skip samples that don't contain that string.
    # This is set by the command-line argument — no editing needed.
    # Example: SPECIES_FILTER=MF skips DGL1, DGL2, etc.
    if [[ -n "${SPECIES_FILTER}" && "${sample}" != *"${SPECIES_FILTER}"* ]]; then
        continue
    fi

    count_total=$((count_total + 1))

    # The FASTA is inside the assembly folder, named after the sample.
    # Example: MFF1_trinity_assembly/MFF1_trinity_assembly.Trinity.fasta
    fasta="${assembly_folder}/${sample}_trinity_assembly.Trinity.fasta"
    stats="${assembly_folder}/Trinity_assembly_stats.txt"

    echo ""
    echo "▶ ${sample}"
    echo "  FASTA: ${fasta}"

    # If the FASTA file is missing, the assembly hasn't finished yet — skip
    if [[ ! -f "${fasta}" ]]; then
        echo "  SKIPPED — FASTA not found (assembly may still be running)"
        echo "  Expected: ${fasta}"
        count_skipped=$((count_skipped + 1))
        continue
    fi

    # If stats already exist from a previous run, reuse them to save time
    if [[ -f "${stats}" ]]; then
        echo "  Reusing existing stats file"
    else
        # Run TrinityStats.pl — saves output to the stats text file
        echo "  Running TrinityStats.pl ..."
        if ! TrinityStats.pl "${fasta}" > "${stats}" 2>&1; then
            echo "  FAILED — TrinityStats.pl returned an error"
            count_failed=$((count_failed + 1))
            continue
        fi
        echo "  Saved: ${stats}"
    fi

    # Print the full stats to the job log
    cat "${stats}"

    # ---- Extract key numbers from the stats file ----
    # grep finds the line containing the label
    # awk '{print $NF}' grabs the last word on that line (the number)
    # The   || echo "N/A"   handles the case where the line isn't found

    transcripts=$(grep -i "Total trinity transcripts" "${stats}" | awk '{print $NF}' || echo "N/A")
    genes=$(grep -i "Total trinity 'genes'"        "${stats}" | awk '{print $NF}' || echo "N/A")
    n50=$(grep  -i "^N50:"                          "${stats}" | awk '{print $NF}' || echo "N/A")
    gc=$(grep   -i "Percent GC"                     "${stats}" | awk '{print $NF}' || echo "N/A")

    # Add a row for this sample to the summary table
    printf "%-20s %14s %14s %10s %10s\n" \
        "${sample}" "${transcripts}" "${genes}" "${n50}" "${gc}" \
        >> "${SUMMARY_FILE}"

    count_done=$((count_done + 1))
    echo "  ✓ Added to summary table"

done


# ---- Final report ----
echo ""
echo "========================================================"
echo "DONE"
echo "========================================================"
echo "Total assemblies found : ${count_total}"
echo "Successfully processed : ${count_done}"
echo "Skipped (no FASTA yet) : ${count_skipped}"
echo "Failed                 : ${count_failed}"
echo ""
echo "Summary table → ${SUMMARY_FILE}"
echo ""
cat "${SUMMARY_FILE}"
echo ""
echo "Finished: $(date)"
echo "========================================================"

exit 0
