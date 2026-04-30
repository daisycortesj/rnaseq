#!/bin/bash
# =============================================================================
# SCRIPT:  trinity_pooled.sbatch
# PURPOSE: Run ONE pooled Trinity de novo transcriptome assembly using all
#          6 Myristica fragrans samples (MFF1-3 + MFL1-3) together.
#
# WHY POOLED (not per-sample):
#   Differential expression analysis requires all samples to share the same
#   transcript IDs. Running Trinity separately per sample produces 6 different
#   assemblies with incompatible IDs — you cannot compare counts across them.
#   Pooling produces one reference transcriptome that all 6 samples are then
#   quantified against (e.g. with RSEM), enabling DESeq2 comparison.
#
# OUTPUT LOCATION:
#   All outputs go to 01_processed/00_4_MF_trinity/
#
#   Trinity appends ".Trinity.fasta" to whatever you give --output, so:
#     --output .../trinity_out
#   produces:
#     .../trinity_out.Trinity.fasta      ← your final assembly
#     .../trinity_out.Trinity.fasta.gene_trans_map
#
# HOW TO SUBMIT:
#   sbatch scripts/03_assembly/trinity_pooled.sbatch
#
# !! DO NOT submit this job twice for the same output directory !!
#    Two Trinity jobs writing to the same folder will corrupt each other.
#    Check squeue first: squeue -u daisycortesj
#
# HOW TO MONITOR PROGRESS:
#   squeue -u daisycortesj                      # is the job running?
#   tail -f 02_trinity_pooled/trinity_pooled_<JOBID>.out   # live log
#   jobload <JOBID>                             # CPU/memory usage on the node
#
# AFTER THIS JOB FINISHES, next steps are:
#   1. BUSCO  — assembly completeness (embryophyta_odb10)
#   2. Kraken2 — decontamination of Trinity.fasta
#   3. TransDecoder — ORF prediction
#   4. RSEM per sample — quantify against the pooled assembly
#   5. DESeq2 — differential expression MFF vs MFL
# =============================================================================


# =============================================================================
# SLURM RESOURCE REQUESTS
# =============================================================================
# QoS note:
#   tc_long  → Tinkercliffs, allows up to 14 days wall time on normal_q
#              This is what you need — the default QoS (tc_normal_base) caps
#              at a much shorter wall time and was likely why your previous
#              job hit the 72-hour limit and did not finish.
#
#   If you have access to OWL (768 GB RAM, 96 cores per node):
#     Change --partition to OWL's partition name
#     Change --qos to owl_long (or equivalent — check ARC docs)
#     Raise --mem to 500G and --cpus-per-task to 48 or 96
#     OWL's 8 GB/core ratio is much better suited for Trinity than Tinkercliffs
#
#   To check available partitions/QoS on ARC:
#     sinfo
#     sacctmgr show qos
# =============================================================================
#SBATCH --job-name=trinity_pooled
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1                   # Trinity must run on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # used by Trinity butterfly parallelization
#SBATCH --mem=240G                  # ~1 GB per million read pairs + overhead
#SBATCH --time=7-00:00:00           # 7 days — generous buffer for ~190M read pairs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity
#SBATCH -o trinity_pooled_%j.out    # stdout log (%j = SLURM job ID)
#SBATCH -e trinity_pooled_%j.err    # stderr log

# Stop the script immediately if any command fails
set -euo pipefail


# =============================================================================
# ABSOLUTE PATHS
# All paths are absolute so this script works regardless of where it is called.
# =============================================================================

# Where your cleaned input reads live
#INPUT_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF"
INPUT_MFF="/projects/tholl_lab_1/daisy_analysis/01_processed/00_00_kraken"
INPUT_MFL="/projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF"

# Where all Trinity output will be saved (separate from per-sample runs)
OUT_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/00_4_MF_trinity"

# The --output argument passed to Trinity
# Trinity will create: ${TRINITY_OUT_PREFIX}.Trinity.fasta
# Output goes directly into OUT_DIR — no extra subfolder created.
TRINITY_OUT_PREFIX="${OUT_DIR}/MF_trinity_pooled"

# R1 (left/forward) reads — comma-separated, NO spaces
# !! ORDER MATTERS: the Nth R1 file must be paired with the Nth R2 file below !!
# If you add or reorder files here, update --right in the exact same order.
#LEFT_READS="\
#${INPUT_DIR}/MFF1_1_clean.fq.gz,\
#${INPUT_DIR}/MFF2_1_clean.fq.gz,\
#${INPUT_DIR}/MFF3_1_clean.fq.gz,\
#${INPUT_DIR}/MFL1_1_clean.fq.gz,\
#${INPUT_DIR}/MFL2_1_clean.fq.gz,\
#${INPUT_DIR}/MFL3_1_clean.fq.gz"

LEFT_READS="\
${INPUT_MFF}/MFF1_1_filtered.fq.gz,\
${INPUT_MFF}/MFF2_1_filtered.fq.gz,\
${INPUT_MFF}/MFF3_1_filtered.fq.gz,\
${INPUT_MFL}/MFL1_1_clean.fq.gz,\
${INPUT_MFL}/MFL2_1_clean.fq.gz,\
${INPUT_MFL}/MFL3_1_clean.fq.gz"

# R2 (right/reverse) reads — must be in the SAME order as LEFT_READS above
# Trinity pairs file 1 of --left with file 1 of --right, file 2 with file 2, etc.
# Mismatched order = wrong reads paired = garbage assembly with no error message.
#RIGHT_READS="\
#${INPUT_DIR}/MFF1_2_clean.fq.gz,\
#${INPUT_DIR}/MFF2_2_clean.fq.gz,\
#${INPUT_DIR}/MFF3_2_clean.fq.gz,\
#${INPUT_DIR}/MFL1_2_clean.fq.gz,\
#${INPUT_DIR}/MFL2_2_clean.fq.gz,\
#${INPUT_DIR}/MFL3_2_clean.fq.gz"

RIGHT_READS="\
${INPUT_MFF}/MFF1_2_filtered.fq.gz,\
${INPUT_MFF}/MFF2_2_filtered.fq.gz,\
${INPUT_MFF}/MFF3_2_filtered.fq.gz,\
${INPUT_MFL}/MFL1_2_clean.fq.gz,\
${INPUT_MFL}/MFL2_2_clean.fq.gz,\
${INPUT_MFL}/MFL3_2_clean.fq.gz"

# NOTE on strand-specificity (--SS_lib_type):
#   If your library is strand-specific (e.g. Illumina TruSeq Stranded),
#   add "--SS_lib_type RF" to the Trinity command below.
#   If you are unsure, contact your sequencing facility and confirm before running.
#   Using the wrong setting (or omitting it when it should be set) reduces assembly
#   quality. When in doubt, omit it — unstranded is a safe default.
#   Uncomment the line below once confirmed:
# SS_LIB_ARG="--SS_lib_type RF"
SS_LIB_ARG=""   # leave empty if unstranded or unconfirmed


# =============================================================================
# INODE SAFETY CHECK
# Inodes are "file slots" on the filesystem — separate from disk space.
# Trinity creates 200,000+ small files during the Chrysalis phase.
# If the inode quota is too low, Trinity crashes mid-run with:
#   "Error, cannot write to file ... exit code 122"
# This check warns you before wasting compute time.
# =============================================================================
echo "Checking inode availability before starting..."

# Get the number of free inodes on the /projects filesystem
INODE_FREE=$(df -i /projects/tholl_lab_1/ | awk 'NR==2 {print $4}')
INODE_PCT=$(df -i /projects/tholl_lab_1/ | awk 'NR==2 {gsub(/%/,"",$5); print $5}')

echo "  Free inodes:    ${INODE_FREE}"
echo "  Inode usage:    ${INODE_PCT}%"

# Warn if inode usage is above 85% — Trinity needs room to create ~300K new files
if [[ "${INODE_PCT}" -ge 85 ]]; then
    echo ""
    echo "WARNING: Inode usage is ${INODE_PCT}% — Trinity may fail mid-run!"
    echo "  Trinity needs ~300,000 free file slots to complete Chrysalis."
    echo "  Clean up old Trinity intermediate folders before submitting."
    echo "  Folders to check:"
    echo "    find /projects/tholl_lab_1/daisy_analysis/ -name 'read_partitions' -type d"
    echo ""
    echo "Proceeding anyway — but expect a possible exit code 122 failure."
fi

# Hard stop if over 95% — almost certain to fail
if [[ "${INODE_PCT}" -ge 95 ]]; then
    echo ""
    echo "ERROR: Inode usage is ${INODE_PCT}% — stopping now."
    echo "  Trinity will definitely fail. Clean up intermediate files first."
    exit 1
fi
echo "Inode check passed."
echo "------------------------------------------------------------"


# =============================================================================
# JOB HEADER — printed to log so you know when/where it ran
# =============================================================================
echo "========================================================"
echo "Trinity Pooled Assembly — Myristica fragrans (MF)"
echo "========================================================"
echo "Job ID:      ${SLURM_JOB_ID:-N/A}"
echo "Host:        $(hostname)"
echo "Working dir: $(pwd)"
echo "Started:     $(date)"
echo "CPUs:        ${SLURM_CPUS_PER_TASK:-32}"
echo "Memory:      ${SLURM_MEM_PER_NODE:-240G}"
echo "========================================================"


# =============================================================================
# MODULE LOADING
# ARC uses the module system — NOT conda — for Trinity.
# Run this first to find the current available version:
#   module spider Trinity
# Then update the version string below to match what is available.
# =============================================================================
echo ""
echo "Loading modules..."
module reset
module load Trinity/2.15.2-foss-2023a   # update version if needed (module spider Trinity)

echo "Trinity version: $(Trinity --version 2>&1 | head -1)"
echo "TrinityStats.pl: $(which TrinityStats.pl)"
echo "------------------------------------------------------------"


# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================
echo ""
echo "Creating output directory: ${OUT_DIR}"
mkdir -p "${OUT_DIR}"


# =============================================================================
# VERIFY ALL INPUT FILES EXIST BEFORE STARTING
# Fails fast with a clear message if any file is missing.
# This prevents Trinity from running for hours and then crashing on a typo.
# =============================================================================
echo ""
echo "Checking input files..."

# Build a combined list of all 12 files (6 R1 + 6 R2) to check
#ALL_FILES=(
#    "${INPUT_DIR}/MFF1_1_clean.fq.gz"
#    "${INPUT_DIR}/MFF2_1_clean.fq.gz"
#    "${INPUT_DIR}/MFF3_1_clean.fq.gz"
#    "${INPUT_DIR}/MFL1_1_clean.fq.gz"
#    "${INPUT_DIR}/MFL2_1_clean.fq.gz"
#    "${INPUT_DIR}/MFL3_1_clean.fq.gz"
#    "${INPUT_DIR}/MFF1_2_clean.fq.gz"
#    "${INPUT_DIR}/MFF2_2_clean.fq.gz"
#    "${INPUT_DIR}/MFF3_2_clean.fq.gz"
#    "${INPUT_DIR}/MFL1_2_clean.fq.gz"
#    "${INPUT_DIR}/MFL2_2_clean.fq.gz"
#    "${INPUT_DIR}/MFL3_2_clean.fq.gz"
#)
ALL_FILES=(
    "${INPUT_MFF}/MFF1_1_filtered.fq.gz"
    "${INPUT_MFF}/MFF2_1_filtered.fq.gz"
    "${INPUT_MFF}/MFF3_1_filtered.fq.gz"
    "${INPUT_MFL}/MFL1_1_clean.fq.gz"
    "${INPUT_MFL}/MFL2_1_clean.fq.gz"
    "${INPUT_MFL}/MFL3_1_clean.fq.gz"
    "${INPUT_MFF}/MFF1_2_filtered.fq.gz"
    "${INPUT_MFF}/MFF2_2_filtered.fq.gz"
    "${INPUT_MFF}/MFF3_2_filtered.fq.gz"
    "${INPUT_MFL}/MFL1_2_clean.fq.gz"
    "${INPUT_MFL}/MFL2_2_clean.fq.gz"
    "${INPUT_MFL}/MFL3_2_clean.fq.gz"
)

# NOTE: these _clean files are adapter-trimmed but NOT decontaminated.
# Decontamination will be done post-assembly with Kraken2 + BLAST on Trinity.fasta.
# This is intentional — do not add a contamination filter step here.

ALL_OK=true
for f in "${ALL_FILES[@]}"; do
    if [[ ! -f "${f}" ]]; then
        echo "  ✗ MISSING: ${f}"
        ALL_OK=false
    else
        SIZE=$(du -sh "${f}" | cut -f1)
        echo "  ✓ ${f}  (${SIZE})"
    fi
done

if [[ "${ALL_OK}" == "false" ]]; then
    echo ""
    echo "ERROR: One or more input files are missing. Fix the paths above and resubmit."
    exit 1
fi
echo "All 12 input files confirmed."
echo "------------------------------------------------------------"


# =============================================================================
# CHECK FOR EXISTING ASSEMBLY — skip or resume
# =============================================================================
FASTA="${TRINITY_OUT_PREFIX}.Trinity.fasta"
RESUME_FLAG=""

if [[ -f "${FASTA}" ]]; then
    echo ""
    echo "✓ Trinity.fasta already exists — skipping assembly, running TrinityStats only"
    echo "  ${FASTA}"
elif [[ -d "${TRINITY_OUT_PREFIX}" ]]; then
    echo ""
    echo "⚠ Partial assembly found — will resume from checkpoint"
    RESUME_FLAG=""
    # Note: --resume was removed in Trinity 2.15.x
    # Trinity automatically detects and resumes from existing checkpoint files
    # in the output directory — no flag needed.
fi


# =============================================================================
# RUN TRINITY (only if FASTA doesn't already exist)
# =============================================================================
if [[ ! -f "${FASTA}" ]]; then

    echo ""
    echo "========================================================"
    echo "Running Trinity..."
    echo "========================================================"

    # Build the Trinity command
    # Each flag is on its own line for readability — the backslash \ continues the command
    TRINITY_CMD="Trinity \
        --seqType fq \
        --left  ${LEFT_READS} \
        --right ${RIGHT_READS} \
        --max_memory 230G \
        --CPU ${SLURM_CPUS_PER_TASK:-32} \
        --output ${TRINITY_OUT_PREFIX} \
        --min_contig_length 200 \
        ${SS_LIB_ARG} \
        ${RESUME_FLAG}"
    # Note: --no_distributed_trinity_exec was REMOVED.
    #   That flag caused Trinity to stop after Chrysalis (Phase 2) and never
    #   run Butterfly (Phase 3), so Trinity.fasta was never created.
    #   It does NOT run Butterfly locally — it skips Butterfly entirely.
    #
    #   When Trinity runs inside a SLURM batch job (which this script is),
    #   it automatically runs all three phases locally on the compute node
    #   without submitting extra jobs to the scheduler. No flag needed.
    #
    # Note: --no_normalize_reads is intentionally OMITTED.
    #   Trinity's in-silico read normalization is helpful for datasets this size.
    #   Skipping it (--no_normalize_reads) would increase memory and runtime
    #   without improving assembly quality for ~190M read pairs.
    #
    # Note: --trimmomatic is intentionally OMITTED.
    #   Your reads already have the _clean suffix (adapter-trimmed by fastp).
    #   Running Trimmomatic again would be redundant.

    # Print the full command to the log before running (useful for debugging)
    echo "Command:"
    echo "${TRINITY_CMD}"
    echo "------------------------------------------------------------"

    # Run Trinity
    if eval ${TRINITY_CMD} 2>&1; then
        echo ""
        echo "✓ Trinity finished successfully"
    else
        TRINITY_EXIT=$?
        echo ""
        echo "✗ Trinity failed (exit code: ${TRINITY_EXIT})"
        echo "  Check the error log: ${OUT_DIR}/trinity_pooled_${SLURM_JOB_ID:-?}.err"
        exit ${TRINITY_EXIT}
    fi

fi


# =============================================================================
# CONFIRM TRINITY OUTPUT EXISTS
# =============================================================================
if [[ ! -f "${FASTA}" ]]; then
    echo ""
    echo "✗ ERROR: Trinity reported success but the FASTA file was not found."
    echo "  Expected: ${FASTA}"
    echo "  Check the error log for details."
    exit 1
fi

CONTIG_COUNT=$(grep -c ">" "${FASTA}" || echo "0")
echo ""
echo "✓ Assembly confirmed: ${CONTIG_COUNT} contigs"
echo "  ${FASTA}"


# =============================================================================
# RUN TRINITYSTATS
# =============================================================================
echo ""
echo "========================================================"
echo "Running TrinityStats..."
echo "========================================================"

STATS_FILE="${OUT_DIR}/Trinity_assembly_stats.txt"

if TrinityStats.pl "${FASTA}" > "${STATS_FILE}" 2>&1; then
    echo "✓ Stats saved: ${STATS_FILE}"
    echo ""
    echo "---- Assembly quality summary ----"
    cat "${STATS_FILE}"
    echo "----------------------------------"
else
    echo "⚠ WARNING: TrinityStats.pl failed (your assembly is still valid)"
    echo "  Re-run manually: TrinityStats.pl ${FASTA}"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ COMPLETE"
echo "========================================================"
echo "Assembly FASTA:  ${FASTA}"
echo "Gene-trans map:  ${FASTA}.gene_trans_map"
echo "Stats:           ${STATS_FILE}"
echo ""
echo "IMPORTANT — Free up inodes now that Trinity is done:"
echo "  The intermediate files in MF_trinity_pooled/ are no longer needed."
echo "  Run this to delete them and free ~200,000 file slots:"
echo ""
echo "    rm -rf ${TRINITY_OUT_PREFIX}/read_partitions"
echo "    rm -rf ${TRINITY_OUT_PREFIX}/insilico_read_normalization"
echo "    rm -rf ${TRINITY_OUT_PREFIX}/chrysalis"
echo "    rm -rf ${TRINITY_OUT_PREFIX}/jellyfish.kmers.fa"
echo ""
echo "Next steps:"
echo "  1. BUSCO:        busco -i ${FASTA} -l embryophyta_odb10 -m transcriptome"
echo "  2. Kraken2:      kraken2 --db <DB> ${FASTA}   (decontamination)"
echo "  3. TransDecoder: TransDecoder.LongOrfs -t ${FASTA}"
echo "  4. RSEM per sample against this pooled assembly"
echo "  5. DESeq2: MFF vs MFL differential expression"
echo ""
echo "Finished: $(date)"
echo "========================================================"

exit 0
