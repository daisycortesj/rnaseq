#!/bin/bash
# =============================================================================
# SCRIPT:  kraken2_prefilter.sbatch
# PURPOSE: Run Kraken2 on 6 paired-end RNA-seq samples to remove microbial
#          contamination BEFORE Trinity de novo assembly.
#
# WHY WE DO THIS BEFORE TRINITY:
#   MultiQC flagged "warn" on per-sequence GC content for all 6 samples.
#   MFF (fruit) samples showed ~50% GC vs MFL (leaf) at ~41% — a 9-point gap
#   that strongly suggests microbial contamination in the fruit samples (likely
#   endophytic fungi or bacteria). Plants typically run 40–45% GC.
#
#   If you assemble without decontamination, Trinity builds transcripts that
#   mix nutmeg sequence with microbial sequence ("chimeric transcripts"). These
#   chimeras produce false gene counts, corrupt DESeq2 results, and make
#   downstream annotation (BLAST, HMMER) unreliable.
#
#   Kraken2 classifies every read by comparing k-mers against a reference
#   database. Reads that DON'T match anything (unclassified) are almost
#   certainly your plant reads → those go to Trinity.
#   Reads that DO match (classified) are contamination → saved separately.
#
# WHERE THIS FITS IN YOUR PIPELINE:
#   Raw reads → fastp (trim) → RiboDetector (rRNA removal) → Kraken2 (THIS) → Trinity
#
# INPUTS:
#   /projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF/
#     MFF1_1_clean.fq.gz    MFF1_2_clean.fq.gz
#     MFF2_1_clean.fq.gz    MFF2_2_clean.fq.gz
#     MFF3_1_clean.fq.gz    MFF3_2_clean.fq.gz
#     MFL1_1_clean.fq.gz    MFL1_2_clean.fq.gz
#     MFL2_1_clean.fq.gz    MFL2_2_clean.fq.gz
#     MFL3_1_clean.fq.gz    MFL3_2_clean.fq.gz
#   (adapter-trimmed and rRNA-removed; NOT yet decontaminated)
#
# OUTPUTS (all in /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/):
#   {SAMPLE}_1_filtered.fq.gz       ← unclassified R1 — FEED THESE TO TRINITY
#   {SAMPLE}_2_filtered.fq.gz       ← unclassified R2 — FEED THESE TO TRINITY
#   {SAMPLE}_classified_1.fq.gz     ← contamination R1 (saved for your records)
#   {SAMPLE}_classified_2.fq.gz     ← contamination R2 (saved for your records)
#   {SAMPLE}_kraken_output.txt.gz   ← per-read classification (large; compressed)
#   {SAMPLE}_kraken_report.txt      ← per-taxon summary (small; human-readable)
#   contamination_summary.tsv       ← comparison table across all 6 samples
#
# HOW TO SUBMIT:
#   sbatch scripts/01_qc/kraken2_prefilter.sbatch
#
# HOW TO MONITOR:
#   squeue -u daisycortesj
#   tail -f /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/kraken2_prefilter_<JOBID>.out
#   jobload <JOBID>
#
# AFTER THIS FINISHES:
#   Update trinity_pooled.sbatch to point INPUT_DIR to:
#     /projects/tholl_lab_1/daisy_analysis/02_kraken_filtered/
#   and change file suffix from _clean.fq.gz to _filtered.fq.gz
#
# REFERENCE:
#   Wood, Lu, Langmead 2019. Genome Biology 20:257. DOI 10.1186/s13059-019-1891-0
# =============================================================================


# =============================================================================
# DATABASE DOWNLOAD INSTRUCTIONS
# Run these ONCE on a login node BEFORE submitting this script.
# The PlusPF database includes: archaea, bacteria, viruses, plasmids, human,
# UniVec_Core, protozoa, and FUNGI. Fungi is the most important addition for
# plant samples — endophytic fungi are the most common nutmeg contaminant.
#
# OPTION A — Download PlusPF yourself (~77 GB download, ~100 GB on disk):
#   mkdir -p /projects/tholl_lab_1/daisy_analysis/databases/kraken2_pluspf
#   cd /projects/tholl_lab_1/daisy_analysis/databases/kraken2_pluspf
#   wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz
#   tar -xzvf k2_pluspf_20240904.tar.gz
#   rm k2_pluspf_20240904.tar.gz
#   Check https://benlangmead.github.io/aws-indexes/k2 for newer release dates
#   and update the URL/filename accordingly.
#
# OPTION B — Check if ARC already has a shared Kraken2 database (saves time!):
#   ls /projects/arcsingularity/kraken2/ 2>/dev/null
#   ls /apps/ 2>/dev/null | grep -i kraken
#   Or email arc@vt.edu and ask: "Do you have a shared Kraken2 PlusPF database?"
#
# After downloading, verify these 3 required files exist:
#   ls /projects/tholl_lab_1/daisy_analysis/databases/kraken2_pluspf/hash.k2d
#   ls /projects/tholl_lab_1/daisy_analysis/databases/kraken2_pluspf/opts.k2d
#   ls /projects/tholl_lab_1/daisy_analysis/databases/kraken2_pluspf/taxo.k2d
# =============================================================================


# =============================================================================
# SLURM RESOURCE REQUESTS
#
# Memory note: PlusPF requires ~100 GB RAM to load into memory at runtime.
#   128G gives a comfortable buffer. If the job is killed with OOM (out of
#   memory), raise this to 150G or 160G.
#
# Threads note: Kraken2 parallelizes well. 16 threads is a good balance
#   between speed and not monopolizing the node for a single job.
#
# Time note: Each sample takes ~5–15 minutes with 16 threads + DB in RAM.
#   6 samples + compression + summary ≈ 1–2 hours. The 12-hour limit is a
#   very generous safety buffer.
#
# To check available QoS options on OWL:
#   sacctmgr show qos
#   sinfo
# =============================================================================
#SBATCH --job-name=kraken2_prefilter
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --qos=owl_normal_base
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/02_kraken_filtered
#SBATCH -o kraken2_prefilter_%j.out
#SBATCH -e kraken2_prefilter_%j.err

# Exit immediately if any command fails, if an unset variable is used,
# or if any command in a pipeline fails
set -euo pipefail


# =============================================================================
# CONDA ENVIRONMENT ACTIVATION
#
# Kraken2 was installed via conda into the "rnaseq" environment:
#   conda activate rnaseq
#   conda install -c bioconda kraken2
#
# HOW conda activation works inside a SLURM job:
#   SLURM jobs do NOT automatically run your ~/.bashrc (where conda normally
#   gets set up when you log in). Without the two "source" and "activate"
#   lines below, the job would fail with "conda: command not found".
#   We manually load the conda shell functions and then activate the environment.
# =============================================================================
echo "Activating conda environment..."

# Find where conda is installed on this system
# Works for both miniconda3 and anaconda3 installations
CONDA_BASE=$(conda info --base)

# Load conda shell functions into this bash session
# This is the equivalent of what your ~/.bashrc does when you log in normally
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate the rnaseq environment where kraken2 is installed
conda activate rnaseq

echo "Kraken2 version: $(kraken2 --version | head -1)"
echo "Kraken2 path:    $(which kraken2)"
echo "Conda env:       ${CONDA_DEFAULT_ENV}"
echo "------------------------------------------------------------"


# =============================================================================
# VARIABLE DEFINITIONS
# All paths and settings are defined here so you only need to edit ONE place
# if anything changes (different database, different input folder, etc.)
# =============================================================================

# Where your cleaned input reads live (output from fastp + RiboDetector)
INPUT_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/00_3_MF"

# Where all Kraken2 output will be saved
OUTPUT_DIR="/projects/tholl_lab_1/daisy_analysis/02_kraken_filtered"

# Path to the Kraken2 PlusPF database (see download instructions above)
# PlusPF includes: bacteria, archaea, viruses, plasmids, human, fungi, protozoa
# Database size: ~77 GB on disk; needs ~100 GB RAM at runtime
KRAKEN_DB="/projects/tholl_lab_1/daisy_analysis/07_NRdatabase/00_1_kraken2_pluspf"

# The 6 sample name prefixes — used to build all input/output file paths
SAMPLES=("MFF1" "MFF2" "MFF3" "MFL1" "MFL2" "MFL3")

# Number of threads to use (matches --cpus-per-task above)
THREADS=16


# =============================================================================
# JOB HEADER — printed to your .out log so you know when/where it ran
# =============================================================================
echo "========================================================"
echo "Kraken2 Pre-Assembly Decontamination — Myristica fragrans (MF)"
echo "========================================================"
echo "Job ID:       ${SLURM_JOB_ID:-N/A}"
echo "Job name:     ${SLURM_JOB_NAME:-N/A}"
echo "Host:         $(hostname)"
echo "Working dir:  $(pwd)"
echo "Started:      $(date)"
echo "CPUs:         ${THREADS}"
echo "Memory:       128G"
echo "Database:     ${KRAKEN_DB}"
echo "Input dir:    ${INPUT_DIR}"
echo "Output dir:   ${OUTPUT_DIR}"
echo "Samples:      ${SAMPLES[*]}"
echo "========================================================"


# =============================================================================
# CREATE OUTPUT DIRECTORY
# mkdir -p creates the directory (and any missing parent directories).
# If it already exists, this does nothing — no error.
# =============================================================================
echo ""
echo "Creating output directory: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"


# =============================================================================
# PRE-FLIGHT CHECKS — verify everything is in place BEFORE processing starts
#
# WHY: Without these checks, the script might run for 30 minutes and then
# crash because of a missing file or wrong database path. Checking upfront
# saves you hours of wasted cluster time.
# =============================================================================
echo ""
echo "Running pre-flight checks..."

# ---- Check 1: Kraken2 database directory and required files ----
# Kraken2 needs exactly these 3 files to work:
#   hash.k2d  — the k-mer index (~68 GB, the main database)
#   opts.k2d  — database build settings
#   taxo.k2d  — NCBI taxonomy tree
#
# If those files already exist, we skip the download entirely.
# If they are missing, we download and extract automatically.
# This means it is safe to re-submit this script — it will never
# re-download a database that is already there.
echo "  Checking Kraken2 database at: ${KRAKEN_DB}"

DB_TAR_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz"
DB_TAR_FILE="${KRAKEN_DB}/k2_pluspf_20240904.tar.gz"

# Check if all 3 required files already exist
if [ -f "${KRAKEN_DB}/hash.k2d" ] && \
   [ -f "${KRAKEN_DB}/opts.k2d" ] && \
   [ -f "${KRAKEN_DB}/taxo.k2d" ]; then

    # Database is already there — skip download completely
    echo "  Database already exists and is complete. Skipping download."

else
    # Database is missing or incomplete — download it now
    echo "  Database not found or incomplete. Downloading PlusPF database..."
    echo "  This is a 66 GB download and may take 1-2 hours."
    echo "  Download started at: $(date)"

    # Create the database directory if it does not exist yet
    mkdir -p "${KRAKEN_DB}"

    # Only download the .tar.gz if it is not already sitting there
    # (handles the case where a previous download finished but extraction failed)
    if [ ! -f "${DB_TAR_FILE}" ]; then
        wget -O "${DB_TAR_FILE}" "${DB_TAR_URL}"
    else
        echo "  Archive already downloaded, skipping wget. Proceeding to extract."
    fi

    # Extract the database
    echo "  Extracting database... $(date)"
    tar -xzvf "${DB_TAR_FILE}" -C "${KRAKEN_DB}"

    # Delete the compressed archive to free up ~66 GB of disk space
    echo "  Removing compressed archive to free disk space..."
    rm "${DB_TAR_FILE}"

    echo "  Database download and extraction complete at: $(date)"

    # Verify the 3 required files are now present after extraction
    for db_file in "hash.k2d" "opts.k2d" "taxo.k2d"; do
        if [ ! -f "${KRAKEN_DB}/${db_file}" ]; then
            echo "ERROR: Extraction finished but ${db_file} is still missing."
            echo "       The download may have been corrupted. Delete the folder"
            echo "       and re-submit: rm -rf ${KRAKEN_DB}"
            exit 1
        fi
    done
fi

echo "  Database OK (all 3 required .k2d files found)"

# ---- Check 2: All 12 input FASTQ files exist ----
# Checks both R1 (_1_clean.fq.gz) and R2 (_2_clean.fq.gz) for every sample.
echo "  Checking input FASTQ files..."
ALL_INPUTS_OK=true

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${INPUT_DIR}/${SAMPLE}_1_clean.fq.gz"
    R2="${INPUT_DIR}/${SAMPLE}_2_clean.fq.gz"

    if [ ! -f "${R1}" ]; then
        echo "  ERROR: Missing R1 for ${SAMPLE}: ${R1}"
        ALL_INPUTS_OK=false
    fi
    if [ ! -f "${R2}" ]; then
        echo "  ERROR: Missing R2 for ${SAMPLE}: ${R2}"
        ALL_INPUTS_OK=false
    fi
done

if [ "${ALL_INPUTS_OK}" = false ]; then
    echo "ERROR: One or more input files are missing. Fix the paths and resubmit."
    exit 1
fi
echo "  Input files OK (12/12 found)"

# ---- Check 3: Output directory is writable ----
if [ ! -w "${OUTPUT_DIR}" ]; then
    echo "ERROR: Output directory is not writable: ${OUTPUT_DIR}"
    echo "       Check permissions with: ls -ld ${OUTPUT_DIR}"
    exit 1
fi
echo "  Output directory is writable"

echo "All pre-flight checks passed. Starting Kraken2..."
echo ""


# =============================================================================
# MAIN LOOP — run Kraken2 on each sample
#
# For each sample, Kraken2:
#   1. Loads the database into RAM (~100 GB — this is why we need 128G)
#   2. Reads both R1 and R2 together as a paired unit
#   3. Compares k-mers (short sequences) from each read against the database
#   4. Classifies each pair as either:
#        UNCLASSIFIED — k-mers don't match anything → likely your nutmeg reads
#        CLASSIFIED   — k-mers match something in the DB → likely contamination
#
# NOTE on the `#` in --unclassified-out and --classified-out:
#   The `#` is a Kraken2 placeholder — Kraken2 replaces it with `1` for R1
#   and `2` for R2. This is NOT a shell wildcard. Do not quote it away or
#   change it. The shell sees the literal `#` and passes it to Kraken2 as-is.
#
# NOTE on pairing order:
#   R1 and R2 MUST be given in the correct order. Kraken2 classifies each
#   R1/R2 pair together. If you accidentally pass R2 where R1 should go,
#   Kraken2 will not error — it will silently produce incorrect results.
# =============================================================================

# Record when the main processing loop starts (used for runtime calculation)
LOOP_START=$(date +%s)

for SAMPLE in "${SAMPLES[@]}"; do

    # Build full paths to input files for this sample
    R1="${INPUT_DIR}/${SAMPLE}_1_clean.fq.gz"
    R2="${INPUT_DIR}/${SAMPLE}_2_clean.fq.gz"

    echo "=========================================="
    echo "Processing sample: ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo "  Start time: $(date)"
    echo "=========================================="

    # ------------------------------------------------------------------
    # Run Kraken2
    # ------------------------------------------------------------------
    # Flag explanations (important for understanding the output):
    #
    #   --db             path to the PlusPF database directory
    #   --threads        number of CPU threads to use (matches SLURM request)
    #   --paired         treat R1 and R2 as paired reads (classify together)
    #   --gzip-compressed  input files are .fq.gz; Kraken2 decompresses on the fly
    #   --unclassified-out reads that did NOT match the DB (your nutmeg reads)
    #                    The `#` is replaced by 1 and 2 for R1 and R2 files
    #                    THESE are the reads you feed to Trinity
    #   --classified-out   reads that DID match the DB (contamination)
    #                    Saved for your records — useful for identifying what
    #                    was living on your nutmeg (fungi, bacteria, etc.)
    #                    !! Do NOT feed these to Trinity !!
    #   --output         per-read classification (one line per read pair)
    #                    Large file — compressed immediately after
    #   --report         per-taxon summary (how many reads matched each organism)
    #                    Small, human-readable — NOT compressed (you'll want to
    #                    read this with cat or a text editor)
    #   --use-names      include scientific names in the output (easier to read
    #                    than taxon IDs alone)
    # ------------------------------------------------------------------
    kraken2 \
        --db "${KRAKEN_DB}" \
        --threads "${THREADS}" \
        --paired \
        --gzip-compressed \
        --unclassified-out "${OUTPUT_DIR}/${SAMPLE}_#_filtered.fq" \
        --classified-out   "${OUTPUT_DIR}/${SAMPLE}_classified_#.fq" \
        --output           "${OUTPUT_DIR}/${SAMPLE}_kraken_output.txt" \
        --report           "${OUTPUT_DIR}/${SAMPLE}_kraken_report.txt" \
        --use-names \
        "${R1}" "${R2}"

    # Check that Kraken2 finished successfully before continuing
    if [ $? -ne 0 ]; then
        echo "ERROR: Kraken2 failed on sample ${SAMPLE}. Check the .err log."
        exit 1
    fi

    # ------------------------------------------------------------------
    # Compress outputs to save disk space
    # ------------------------------------------------------------------
    # Uncompressed Kraken2 outputs can be many GB per sample. We compress
    # immediately after each sample so we don't run out of disk space
    # before all 6 samples finish.
    #
    # NOTE on file names after compression:
    #   The `#` in --unclassified-out produced these actual files:
    #     ${SAMPLE}_1_filtered.fq   (R1 unclassified)
    #     ${SAMPLE}_2_filtered.fq   (R2 unclassified)
    #   So we compress those specific names here.
    # ------------------------------------------------------------------
    echo "  Compressing outputs for ${SAMPLE}..."

    gzip "${OUTPUT_DIR}/${SAMPLE}_1_filtered.fq"        # unclassified R1 → .fq.gz
    gzip "${OUTPUT_DIR}/${SAMPLE}_2_filtered.fq"        # unclassified R2 → .fq.gz
    gzip "${OUTPUT_DIR}/${SAMPLE}_classified_1.fq"      # classified R1   → .fq.gz
    gzip "${OUTPUT_DIR}/${SAMPLE}_classified_2.fq"      # classified R2   → .fq.gz
    gzip "${OUTPUT_DIR}/${SAMPLE}_kraken_output.txt"    # per-read table  → .txt.gz

    echo "  Finished ${SAMPLE} at $(date)"
    echo ""

done


# =============================================================================
# SUMMARY TABLE — contamination_summary.tsv
#
# After all 6 samples are processed, parse each sample's _kraken_report.txt
# and extract the key contamination metrics into one easy-to-read table.
#
# Kraken2 report format (tab-separated, 6 columns):
#   Col 1: % of reads covered by this clade (cumulative)
#   Col 2: number of reads covered by this clade (cumulative)
#   Col 3: number of reads directly assigned to this taxon
#   Col 4: rank code (U=unclassified, R=root, D=domain, K=kingdom, P=phylum...)
#   Col 5: NCBI taxon ID
#   Col 6: scientific name (indented by spaces to show taxonomy depth)
#
# Rank codes we care about:
#   U  = Unclassified (didn't match anything → likely your nutmeg reads)
#   D  = Domain level entries (Bacteria, Archaea, Eukaryota, Viruses)
#   K  = Kingdom level entries (Fungi lives here under Eukaryota > Fungi)
#
# To look up something specific in a report later:
#   grep "Fungi" MFF1_kraken_report.txt
#   grep "Bacteria" MFF1_kraken_report.txt
# =============================================================================
echo "=========================================="
echo "Building contamination_summary.tsv..."
echo "=========================================="

SUMMARY_FILE="${OUTPUT_DIR}/contamination_summary.tsv"

# Write the header row (tab-separated)
printf "Sample\tTotal_Reads\tUnclassified_pct\tClassified_pct\tBacterial_pct\tFungal_pct\tViral_pct\n" \
    > "${SUMMARY_FILE}"

for SAMPLE in "${SAMPLES[@]}"; do

    REPORT="${OUTPUT_DIR}/${SAMPLE}_kraken_report.txt"

    if [ ! -f "${REPORT}" ]; then
        echo "  WARNING: Report not found for ${SAMPLE}: ${REPORT}. Skipping."
        continue
    fi

    # Extract % unclassified — rank code "U" in column 4
    # awk: split each line by tab, if column 4 is exactly "U", print column 1
    UNCLASSIFIED_PCT=$(awk -F'\t' '$4 == "U" {printf "%.2f", $1}' "${REPORT}")

    # Extract % classified = 100 - unclassified
    # bc is a command-line calculator that handles decimal math
    CLASSIFIED_PCT=$(echo "scale=2; 100 - ${UNCLASSIFIED_PCT}" | bc)

    # Extract % bacterial — look for domain "Bacteria" (rank code D)
    # $6 may have leading spaces (indentation in the report); gsub removes them
    BACTERIAL_PCT=$(awk -F'\t' '$4 == "D" && gsub(/^ +/, "", $6) && $6 == "Bacteria" {printf "%.2f", $1}' "${REPORT}")
    BACTERIAL_PCT="${BACTERIAL_PCT:-0.00}"   # default to 0.00 if not found

    # Extract % fungal — Fungi is at Kingdom rank (rank code K)
    FUNGAL_PCT=$(awk -F'\t' '$4 == "K" && gsub(/^ +/, "", $6) && $6 == "Fungi" {printf "%.2f", $1}' "${REPORT}")
    FUNGAL_PCT="${FUNGAL_PCT:-0.00}"

    # Extract % viral — look for domain "Viruses" (rank code D)
    VIRAL_PCT=$(awk -F'\t' '$4 == "D" && gsub(/^ +/, "", $6) && $6 == "Viruses" {printf "%.2f", $1}' "${REPORT}")
    VIRAL_PCT="${VIRAL_PCT:-0.00}"

    # Get total reads from the report — sum of reads in all root-level entries
    # The first "U" line contains the total unclassified count; the "R" (root)
    # line contains total classified. We add them to get overall total reads.
    UNCLASSIFIED_COUNT=$(awk -F'\t' '$4 == "U" {print $2}' "${REPORT}")
    CLASSIFIED_COUNT=$(awk -F'\t' '$4 == "R" {print $2}' "${REPORT}")
    TOTAL_READS=$(echo "${UNCLASSIFIED_COUNT:-0} + ${CLASSIFIED_COUNT:-0}" | bc)

    # Append this sample's row to the summary file
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE}" \
        "${TOTAL_READS}" \
        "${UNCLASSIFIED_PCT}" \
        "${CLASSIFIED_PCT}" \
        "${BACTERIAL_PCT}" \
        "${FUNGAL_PCT}" \
        "${VIRAL_PCT}" \
        >> "${SUMMARY_FILE}"

    echo "  ${SAMPLE}: ${UNCLASSIFIED_PCT}% unclassified (kept for Trinity), ${CLASSIFIED_PCT}% classified (removed)"

done

echo ""
echo "Summary table written to: ${SUMMARY_FILE}"
echo ""


# =============================================================================
# COMPLETION MESSAGE
# =============================================================================

# Calculate total runtime
LOOP_END=$(date +%s)
ELAPSED=$(( LOOP_END - LOOP_START ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo "========================================================"
echo "Kraken2 pre-assembly decontamination COMPLETE"
echo "========================================================"
echo "Finished:      $(date)"
echo "Total runtime: ${ELAPSED_MIN} minutes ${ELAPSED_SEC} seconds"
echo "Processed:     ${#SAMPLES[@]} samples (${SAMPLES[*]})"
echo "Output dir:    ${OUTPUT_DIR}"
echo ""
echo "Files created per sample:"
echo "  {SAMPLE}_1_filtered.fq.gz       ← UNCLASSIFIED R1 — use these for Trinity"
echo "  {SAMPLE}_2_filtered.fq.gz       ← UNCLASSIFIED R2 — use these for Trinity"
echo "  {SAMPLE}_classified_1.fq.gz     ← classified R1 (contamination, saved)"
echo "  {SAMPLE}_classified_2.fq.gz     ← classified R2 (contamination, saved)"
echo "  {SAMPLE}_kraken_output.txt.gz   ← per-read classification (compressed)"
echo "  {SAMPLE}_kraken_report.txt      ← per-taxon summary (human-readable)"
echo ""
echo "REVIEW CONTAMINATION:"
echo "  cat ${OUTPUT_DIR}/contamination_summary.tsv"
echo "  Compare MFF (fruit) vs MFL (leaf) — expect MFF to show higher"
echo "  fungal/bacterial contamination based on the GC% discrepancy."
echo ""
echo "IMPORTANT REMINDERS:"
echo "  - Do NOT delete the original _clean.fq.gz files. Keep them as backup."
echo "  - 'Classified' does not always mean contamination. If a read matched"
echo "    a plant genome in PlusPF (rare — nutmeg is not in the DB), it would"
echo "    appear as classified even though it is a real nutmeg read. This"
echo "    effect is minimal with PlusPF but worth knowing."
echo ""
echo "NEXT STEP — update trinity_pooled.sbatch:"
echo "  Change INPUT_DIR to: ${OUTPUT_DIR}"
echo "  Change file suffix from _clean.fq.gz to _filtered.fq.gz"
echo "  Then submit: sbatch scripts/03_assembly/trinity_pooled.sbatch"
echo "========================================================"
