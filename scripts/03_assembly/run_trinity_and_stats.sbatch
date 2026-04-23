#!/bin/bash
# =============================================================================
# SCRIPT:  run_trinity_and_stats.sbatch
# WHAT IT DOES:
#   Runs Trinity de novo assembly AND TrinityStats for each sample.
#   Uses a SLURM job array so all samples run at the same time in parallel.
#   Based on run_trinity_array.sbatch — do not modify that script, use this one.
#
# HOW TO USE (two steps):
#
#   Step 1 — Generate the sample list (run once per species, no sbatch):
#     bash scripts/03_assembly/run_trinity_and_stats.sbatch --generate MF
#     bash scripts/03_assembly/run_trinity_and_stats.sbatch --generate DG
#
#     This scans your reads and writes 01_processed/sample_list.txt
#     It also prints the number of samples — you need that for Step 2.
#
#   Step 2 — Submit the array (replace N with the number from Step 1):
#     sbatch --array=1-N scripts/03_assembly/run_trinity_and_stats.sbatch
#     Example: sbatch --array=1-6 scripts/03_assembly/run_trinity_and_stats.sbatch
#
# RESUME:
#   If a job times out, resubmit the same command.
#   If Trinity.fasta already exists it skips Trinity and goes straight to stats.
#
# AFTER ALL JOBS FINISH — get one combined summary table:
#   sbatch scripts/03_assembly/run_trinitystats.sbatch MF
#
# OUTPUTS (per sample):
#   01_processed/MFF1_trinity_assembly/                                       ← folder
#   01_processed/MFF1_trinity_assembly/MFF1_trinity_assembly.Trinity.fasta    ← sequences
#   01_processed/MFF1_trinity_assembly/Trinity_assembly_stats.txt             ← quality metrics
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


# ---- SLURM cluster resource requests ----
# Trinity needs a lot of memory and time.
# TrinityStats runs at the end of the same job — it finishes in seconds.
#SBATCH --job-name=trinity_stats
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_stats_%A_%a.out
#SBATCH -e trinity_stats_%A_%a.err

# Stop immediately if any command fails
set -euo pipefail

# Load shared folder paths (PROCESSED_DIR, RAWDATA_DIR, CODE_DIR, CONDA_ENV)
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"


# =============================================================================
# STEP 1 MODE: --generate
# Run as: bash run_trinity_and_stats.sbatch --generate MF
# Writes sample_list.txt and tells you how many samples were found.
# =============================================================================
if [[ "${1:-}" == "--generate" ]]; then

    SPECIES_CODE="${2:-}"
    if [[ -z "${SPECIES_CODE}" ]]; then
        echo "Usage: bash run_trinity_and_stats.sbatch --generate <SPECIES_CODE>"
        echo "Example: bash run_trinity_and_stats.sbatch --generate MF"
        exit 1
    fi

    # Look up folder name and species info from config.sh
    get_sample_info "${SPECIES_CODE}"
    SAMPLE_LIST="${PROCESSED_DIR}/sample_list.txt"

    # Auto-detect which reads to use — best available wins:
    #   Priority 1: rRNA-removed reads (ribofree/) — cleanest
    #   Priority 2: fastp-cleaned reads (*_clean.fq.gz)
    #   Priority 3: raw reads — fallback
    CLEAN_DIR="${PROCESSED_DIR}/${SPECIES_DIR}"
    RIBOFREE_DIR="${PROCESSED_DIR}/${SPECIES_DIR}/ribofree"
    RIBOFREE_COUNT=$(find "${RIBOFREE_DIR}" -maxdepth 1 -name "*_1_ribofree.fq.gz" 2>/dev/null | wc -l)
    CLEAN_COUNT=$(find "${CLEAN_DIR}"    -maxdepth 1 -name "*_1_clean.fq.gz"   2>/dev/null | wc -l)

    if [[ "${RIBOFREE_COUNT}" -gt 0 ]]; then
        find "${RIBOFREE_DIR}" -maxdepth 1 -name "*_1_ribofree.fq.gz" | sort > "${SAMPLE_LIST}"
        echo "Using RIBOFREE reads (rRNA removed) — best quality"
    elif [[ "${CLEAN_COUNT}" -gt 0 ]]; then
        find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" | sort > "${SAMPLE_LIST}"
        echo "Using FASTP-CLEANED reads"
    else
        find "${RAWDATA_DIR}/${SPECIES_DIR}" -maxdepth 1 -name "*_1.fq.gz" | sort > "${SAMPLE_LIST}"
        echo "Using RAW reads (no cleaned reads found)"
    fi

    SAMPLE_COUNT=$(wc -l < "${SAMPLE_LIST}")
    echo "Species: ${SPECIES_CODE} (${SPECIES_NAME})"
    echo "Samples: ${SAMPLE_COUNT}"
    echo "List:    ${SAMPLE_LIST}"
    echo ""
    cat "${SAMPLE_LIST}"
    echo ""
    echo "Now submit with:"
    echo "  sbatch --array=1-${SAMPLE_COUNT} scripts/03_assembly/run_trinity_and_stats.sbatch"
    exit 0
fi


# =============================================================================
# STEP 2 MODE: assembly + stats (runs when submitted with sbatch --array)
# =============================================================================

SAMPLE_LIST="${PROCESSED_DIR}/sample_list.txt"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
MEM_GB="128"

echo "========================================================"
echo "Trinity Assembly + TrinityStats"
echo "========================================================"
echo "Array job:  ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "This task:  ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Host:       $(hostname)"
echo "CPUs: ${THREADS}   MEM: ${SLURM_MEM_PER_NODE:-128G}"
echo "Started:    $(date)"
echo "========================================================"


# ---- Validate sample list ----
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    echo "Run Step 1 first:"
    echo "  bash scripts/03_assembly/run_trinity_and_stats.sbatch --generate MF"
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Total samples in list: ${TOTAL_SAMPLES}"

if [[ ${SLURM_ARRAY_TASK_ID:-1} -gt ${TOTAL_SAMPLES} ]]; then
    echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} exceeds sample count ${TOTAL_SAMPLES}"
    exit 1
fi


# ---- Get this task's sample from the list ----
# Each task picks a different line from sample_list.txt.
# Task 1 gets line 1 (MFF1), task 2 gets line 2 (MFF2), etc.
LEFT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" "${SAMPLE_LIST}")
sample_dir=$(dirname "${LEFT_FILE}")
filename=$(basename "${LEFT_FILE}")

# Figure out the sample name and read suffix from the filename
# MFF1_1_ribofree.fq.gz → base=MFF1, suffix=_ribofree
# MFF1_1_clean.fq.gz    → base=MFF1, suffix=_clean
# MFF1_1.fq.gz          → base=MFF1, suffix=
if [[ "${filename}" == *"_1_ribofree.fq.gz" ]]; then
    base="${filename%_1_ribofree.fq.gz}"
    FILE_SUFFIX="_ribofree"
elif [[ "${filename}" == *"_1_clean.fq.gz" ]]; then
    base="${filename%_1_clean.fq.gz}"
    FILE_SUFFIX="_clean"
else
    base="${filename%_1.fq.gz}"
    FILE_SUFFIX=""
fi

# Build the R2 (right/reverse) file path from the R1 path
RIGHT_FILE="${sample_dir}/${base}_2${FILE_SUFFIX}.fq.gz"

echo ""
echo "Sample: ${base}  (task ${SLURM_ARRAY_TASK_ID:-1} of ${TOTAL_SAMPLES})"
echo "  R1: ${LEFT_FILE}"
echo "  R2: ${RIGHT_FILE}"

# Confirm both read files exist before starting
if [[ ! -r "${LEFT_FILE}" ]]; then
    echo "ERROR: Cannot read R1 file: ${LEFT_FILE}"
    exit 1
fi
if [[ ! -r "${RIGHT_FILE}" ]]; then
    echo "ERROR: Cannot read R2 file: ${RIGHT_FILE}"
    exit 1
fi
echo "  ✓ Both read files found"


# ---- Activate conda ----
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    source ~/miniconda3/etc/profile.d/conda.sh
fi
source ~/.bashrc
conda activate "${CONDA_ENV}"

# Confirm Trinity is available
if ! command -v Trinity &> /dev/null; then
    echo "ERROR: Trinity not found. Install: conda install -c bioconda trinity"
    exit 1
fi
echo "Trinity: $(Trinity --version 2>&1 | head -1)"

# Go to the code directory (required by the Python CLI wrapper)
cd "${CODE_DIR}"

# Paths for this sample's output
OUTPUT_DIR="${PROCESSED_DIR}/${base}_trinity_assembly"
FASTA="${OUTPUT_DIR}/${base}_trinity_assembly.Trinity.fasta"
STATS_FILE="${OUTPUT_DIR}/Trinity_assembly_stats.txt"


# =============================================================================
# PART 1 — TRINITY ASSEMBLY
# =============================================================================
echo ""
echo "------------------------------------------------------------"
echo "Part 1: Trinity Assembly"
echo "------------------------------------------------------------"
echo "Output folder: ${OUTPUT_DIR}"

# Check if assembly already finished — skip Trinity if so
RESUME_FLAG=""
if [[ -f "${FASTA}" ]]; then
    echo "✓ FASTA already exists — skipping Trinity, going straight to TrinityStats"
elif [[ -d "${OUTPUT_DIR}" ]]; then
    echo "⚠ Partial assembly found — resuming from checkpoint"
    RESUME_FLAG="--resume"
fi

# Run Trinity only if the FASTA doesn't exist yet
if [[ ! -f "${FASTA}" ]]; then
    CMD="python -m rna_pipeline.cli \
      --mode trinity \
      --reads-left \"${LEFT_FILE}\" \
      --reads-right \"${RIGHT_FILE}\" \
      --outdir \"${OUTPUT_DIR}\" \
      --threads ${THREADS} \
      --mem-gb ${MEM_GB} \
      ${RESUME_FLAG}"

    echo "Running: ${CMD}"
    echo "------------------------------------------------------------"

    if eval ${CMD} 2>&1; then
        echo "✓ Trinity finished for ${base}"
    else
        EXIT_CODE=$?
        echo "✗ Trinity failed for ${base} (exit code: ${EXIT_CODE})"
        exit ${EXIT_CODE}
    fi
fi

# Confirm the FASTA was created
if [[ ! -f "${FASTA}" ]]; then
    echo "✗ ERROR: Trinity finished but FASTA not found at:"
    echo "  ${FASTA}"
    exit 1
fi

CONTIG_COUNT=$(grep -c ">" "${FASTA}" || echo "0")
echo "✓ FASTA: ${CONTIG_COUNT} contigs"
echo "  ${FASTA}"


# =============================================================================
# PART 2 — TRINITYSTATS
# =============================================================================
echo ""
echo "------------------------------------------------------------"
echo "Part 2: TrinityStats"
echo "------------------------------------------------------------"

if ! command -v TrinityStats.pl &> /dev/null; then
    echo "⚠ WARNING: TrinityStats.pl not found — skipping stats"
    echo "  Install: conda install -c bioconda trinity"
else
    # If stats already exist reuse them, otherwise calculate fresh
    if [[ -f "${STATS_FILE}" ]]; then
        echo "Stats file already exists — reusing it"
    else
        echo "Running TrinityStats.pl ..."
        if TrinityStats.pl "${FASTA}" > "${STATS_FILE}" 2>&1; then
            echo "✓ Stats saved: ${STATS_FILE}"
        else
            echo "⚠ WARNING: TrinityStats.pl failed (assembly is still valid)"
            echo "  Re-run later: sbatch scripts/03_assembly/run_trinitystats.sbatch"
        fi
    fi

    echo ""
    echo "---- Stats for ${base} ----"
    cat "${STATS_FILE}"
    echo "---------------------------"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ ${base} complete"
echo "  Assembly: ${FASTA}"
echo "  Stats:    ${STATS_FILE}"
echo ""
echo "When ALL samples are done, run this for a combined summary:"
echo "  sbatch scripts/03_assembly/run_trinitystats.sbatch MF"
echo "Finished: $(date)"
echo "========================================================"

exit 0
