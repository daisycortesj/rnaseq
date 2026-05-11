#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Align RNA-seq reads using STAR (one species or all)
# USAGE:   sbatch run_star_align_all.sbatch [SPECIES_CODE]
#
# EXAMPLES:
#   sbatch run_star_align_all.sbatch MF    # only nutmeg
#   sbatch run_star_align_all.sbatch DC    # only D. carota
#   sbatch run_star_align_all.sbatch       # all species
#
# AUTO-DETECTION (no manual settings needed):
#   For each species, the script checks if fastp-cleaned reads exist
#   in 01_processed/. If they do, it uses those. If not, it falls
#   back to raw reads from 00_rawdata/.
#
#   Example:
#     00_1_DC → no *_clean.fq.gz in 01_processed/ → uses RAW reads
#     00_3_MF → has *_clean.fq.gz in 01_processed/ → uses CLEAN reads
#
# GENOME MAPPING:
#   - 00_1_DC → carrot genome (with GTF → gene counts)
#   - 00_2_DG → carrot genome (with GTF → gene counts)
#   - 00_3_MF → nutmeg genome (no GTF  → BAM only)
#
# PATHS: All directory paths come from scripts/config.sh
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=star_align_all
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/02_mapped
#SBATCH -o alignment_sample_%j.out  
#SBATCH -e alignment_sample_%j.err       

set -euo pipefail

# ===== Load shared config =====
# config.sh defines all our directory paths:
#   RAWDATA_DIR   = 00_rawdata/    (raw FASTQ files)
#   PROCESSED_DIR = 01_processed/  (fastp-cleaned reads)
#   MAPPED_DIR    = 02_mapped/     (STAR output — BAMs go here)
#   REFERENCE_DIR = 04_reference/  (genome indexes)
#   CODE_DIR      = 05_rnaseq-code/ (this repo / Python pipeline)
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# --- Optional: filter to a single species ---
# If a species code is given (e.g., MF), only process that species.
# If no code is given, process all species directories.
SPECIES_FILTER="${1:-}"
if [[ -n "${SPECIES_FILTER}" ]]; then
    get_sample_info "${SPECIES_FILTER}"
    FILTER_DIR="${SPECIES_DIR}"
fi

# Directory paths
ALIGNMENT_DIR="${MAPPED_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Kraken2-filtered reads directory (used for nutmeg FRUIT samples only)
# MFF* samples had elevated GC content (~50%) suggesting microbial contamination,
# so they were run through Kraken2 to remove non-plant reads before alignment.
# MFL* (leaf) samples stay on the normal ribofree → clean → raw path.
KRAKEN_DIR="${PROCESSED_DIR}/00_00_kraken"

echo "============================================================"
if [[ -n "${SPECIES_FILTER}" ]]; then
    echo "STAR Read Alignment — ${SPECIES_NAME} (${SPECIES_FILTER})"
else
    echo "STAR Read Alignment Pipeline - Process All Samples"
fi
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo ""
echo "Data source: AUTO-DETECT per species"
echo "  (uses clean reads if available, otherwise raw reads)"
echo ""
echo "Genome mapping:"
echo "  - 00_1_DC → carrot (with GTF)"
echo "  - 00_2_DG → carrot (with GTF)"
echo "  - 00_3_MF → nutmeg (no GTF)"
echo "Start time: $(date)"
echo "============================================================"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    echo "Loading .bashrc successfully"
    source ~/.bashrc
    conda activate rnaseq
fi

# Verify conda activation
conda info --envs | grep rnaseq || echo "ERROR: conda environment not activated"

# Check tools
echo "Checking tool versions:"
STAR --version || echo "STAR not found."
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Process All Sample Directories =====
echo "Scanning species directories in ${RAWDATA_DIR}/ ..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

# ====================================================================
# MAIN PROCESSING LOOP: Process all sample directories
# ====================================================================
# We always discover species from 00_rawdata/ (the canonical list).
# For each species, we then auto-detect whether fastp-cleaned reads
# exist in 01_processed/. If they do, we use those. Otherwise we
# fall back to the raw reads.
#
# Per-species logic:
#   1. Determine genome (carrot or nutmeg) and quant mode
#   2. Auto-detect: clean reads in 01_processed/ or raw in 00_rawdata/
#   3. Find all R1 files, pair with R2
#   4. Run STAR alignment for each sample pair
# ====================================================================
for sample_dir_path in "${RAWDATA_DIR}"/00_*/ ; do
    # Extract directory name (e.g., 00_1_DC, 00_3_MF)
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    # If a species filter was given, skip directories that don't match
    if [[ -n "${SPECIES_FILTER}" ]] && [[ "${SAMPLE_DIR}" != "${FILTER_DIR}" ]]; then
        continue
    fi
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
    # --- Determine genome based on directory name ---
    case "${SAMPLE_DIR}" in
        00_1_DC|00_2_DG)
            GENOME_TYPE="carrot"
            GENOME_INDEX="${REFERENCE_DIR}/carrot_star_index"
            QUANT_MODE="--quant-mode"
            echo "  Genome: CARROT (with GTF → gene quantification enabled)"
            ;;
        00_3_MF)
            GENOME_TYPE="nutmeg"
            GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
            QUANT_MODE="--no-quant-mode"
            echo "  Genome: NUTMEG (no GTF → alignment only, no gene counts)"
            ;;
        *)
            echo "  ⚠ WARNING: Unknown directory '${SAMPLE_DIR}', skipping..."
            continue
            ;;
    esac
    
    # --- Validate genome index exists ---
    if [[ ! -d "${GENOME_INDEX}" ]]; then
        echo "  ✗ ERROR: Genome index not found: ${GENOME_INDEX}"
        echo "  Please build the index first using run_genome_index.sbatch"
        continue
    fi
    echo "  ✓ Genome index: ${GENOME_INDEX}"
    
    # --- Auto-detect input reads: ribofree > clean > raw ---
    # Priority 1: RiboDetector-filtered reads (ribofree/) — use if available
    # Priority 2: Fastp-cleaned reads (*_clean.fq.gz) — use if no ribofree
    # Priority 3: Raw reads — fallback if nothing else exists
    CLEAN_DIR="${PROCESSED_DIR}/${SAMPLE_DIR}"
    RIBOFREE_DIR="${PROCESSED_DIR}/${SAMPLE_DIR}/ribofree"
    RIBOFREE_COUNT=$(find "${RIBOFREE_DIR}" -maxdepth 1 -name "*_1_ribofree.fq.gz" 2>/dev/null | wc -l)
    CLEAN_COUNT=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*_1_clean.fq.gz" 2>/dev/null | wc -l)

    if [[ "${RIBOFREE_COUNT}" -gt 0 ]]; then
        # RiboDetector-filtered reads found → use these (rRNA removed)
        SPECIES_INPUT_DIR="${RIBOFREE_DIR}"
        FILE_PATTERN="*_1_ribofree.fq.gz"
        FILE_SUFFIX="_ribofree"
        INPUT_SOURCE="RIBOFREE (rRNA removed)"
    elif [[ "${CLEAN_COUNT}" -gt 0 ]]; then
        # Fastp-cleaned reads found → use them
        SPECIES_INPUT_DIR="${CLEAN_DIR}"
        FILE_PATTERN="*_1_clean.fq.gz"
        FILE_SUFFIX="_clean"
        INPUT_SOURCE="FASTP-CLEANED"
    else
        # No clean reads → fall back to raw reads
        SPECIES_INPUT_DIR="${RAWDATA_DIR}/${SAMPLE_DIR}"
        FILE_PATTERN="*_1.fq.gz"
        FILE_SUFFIX=""
        INPUT_SOURCE="RAW"
    fi
    
    echo "  Data source: ${INPUT_SOURCE}"
    echo "  Input dir:   ${SPECIES_INPUT_DIR}"
    
    # --- Set output directory ---
    # All alignment output goes to 02_mapped/<SPECIES_DIR>/
    #   e.g., 02_mapped/00_1_DC/  or  02_mapped/00_3_MF/
    OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
    echo "  Output dir:  ${OUTDIR}"
    
    mkdir -p "${OUTDIR}"
    
    # --- Count samples ---
    SAMPLE_COUNT=$(find "${SPECIES_INPUT_DIR}" -maxdepth 1 -name "${FILE_PATTERN}" 2>/dev/null | wc -l)
    
    if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
        echo "  ⚠ No samples found in ${SAMPLE_DIR}, skipping..."
        continue
    fi
    
    echo "  Samples found: ${SAMPLE_COUNT}"
    echo "------------------------------------------------------------"
    
    # Loop through all R1 files matching the pattern in this species' input dir
    for file in "${SPECIES_INPUT_DIR}"/${FILE_PATTERN} ; do
        
        # Extract sample base name from the file path
        # Example for raw:   MFF1_1.fq.gz       → base = MFF1
        # Example for clean: MFF1_1_clean.fq.gz  → base = MFF1
        filename="${file##*/}"
        base="${filename%_1${FILE_SUFFIX}.fq.gz}"
        
        # Construct full paths for paired-end reads (R1 + R2)
        # R1 = forward reads, R2 = reverse reads (from the same biological sample)
        READS_R1="${SPECIES_INPUT_DIR}/${base}_1${FILE_SUFFIX}.fq.gz"
        READS_R2="${SPECIES_INPUT_DIR}/${base}_2${FILE_SUFFIX}.fq.gz"

        # ---------------------------------------------------------------
        # NUTMEG FRUIT OVERRIDE: MFF* samples use Kraken2-filtered reads
        # ---------------------------------------------------------------
        # MFF (fruit) samples showed ~50% GC content vs MFL (leaf) at ~41%,
        # which indicated microbial contamination. Those samples were run
        # through Kraken2 to remove non-plant reads. We use those cleaner
        # filtered reads here instead of the ribofree reads.
        # MFL (leaf) samples are fine — they keep the auto-detected path above.
        if [[ "${GENOME_TYPE}" == "nutmeg" ]] && [[ "${base}" == MFF* ]]; then
            READS_R1="${KRAKEN_DIR}/${base}_1_filtered.fq.gz"
            READS_R2="${KRAKEN_DIR}/${base}_2_filtered.fq.gz"
            echo "  [FRUIT] Using Kraken2-filtered reads: ${KRAKEN_DIR}"
        fi

        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))  # Increment sample counter
        
        # STEP 5: Display sample information
        # ------------------------------------------------------------
        echo ""
        echo "▶ Sample ${TOTAL_SAMPLES}: ${base}"      # Shows sample base name (e.g., DC1L1)
        echo "  R1: ${READS_R1}"                        # Shows forward reads path
        echo "  R2: ${READS_R2}"                        # Shows reverse reads path
        
        # STEP 6: Validate that both paired files exist on disk
        # ------------------------------------------------------------
        # Check R1 file (forward reads)
        if [[ ! -f "${READS_R1}" ]]; then
            echo "  ✗ ERROR: R1 file not found"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue  # Skip to next sample if R1 missing
        fi
        
        # Check R2 file (reverse reads)
        if [[ ! -f "${READS_R2}" ]]; then
            echo "  ✗ ERROR: R2 file not found"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue  # Skip to next sample if R2 missing
        fi
        
        echo "  ✓ Both read files found"  # Both files exist, proceed with alignment
        
        # STEP 7: Build the alignment command for this sample
        # ------------------------------------------------------------
        # This constructs the Python command that will:
        #   1. Call the rna_pipeline CLI
        #   2. Pass the genome index path
        #   3. Pass both R1 and R2 files (--reads-left and --reads-right)
        #   4. Specify output directory and number of threads
        #   5. Enable/disable gene quantification based on genome type
        CMD="python -m rna_pipeline.cli \
          --mode align \
          --genome-index ${GENOME_INDEX} \
          --reads-left ${READS_R1} \
          --reads-right ${READS_R2} \
          --outdir ${OUTDIR} \
          --threads ${THREADS} \
          ${QUANT_MODE} \
          --sample-name ${base}"
        
        # Example command for DC1L1 using carrot genome:
        # python -m rna_pipeline.cli --mode align \
        #   --genome-index /path/to/carrot_star_index \
        #   --reads-left /path/to/DC1L1_1.fq.gz \
        #   --reads-right /path/to/DC1L1_2.fq.gz \
        #   --outdir /path/to/output \
        #   --threads 16 \
        #   --quant-mode \
        #   --sample-name DC1L1
        
        echo "  Running alignment..."
        
        # STEP 8: Execute the alignment command
        # ------------------------------------------------------------
        # This runs the command and tracks success/failure
        # The command will:
        #   → Call Python CLI (rna_pipeline/cli.py)
        #   → Which calls main.py (validates files)
        #   → Which calls star.py (builds STAR command with --readFilesIn R1 R2)
        #   → Which calls local.py (executes STAR via subprocess)
        #   → STAR reads both FASTQ files and aligns them as paired-end reads
        if eval ${CMD}; then
            echo "  ✓ ${base} completed successfully"
            SUCCESSFUL_SAMPLES=$((SUCCESSFUL_SAMPLES + 1))
            
            # OUTPUT FILES CREATED (for sample DC1L1 example):
            #   ${OUTDIR}/DC1L1_Aligned.sortedByCoord.out.bam  - Main BAM alignment file
            #   ${OUTDIR}/DC1L1_ReadsPerGene.out.tab           - Gene count matrix (if quant-mode)
            #   ${OUTDIR}/DC1L1_Log.final.out                  - Alignment statistics
            #   ${OUTDIR}/DC1L1_Log.out                        - Detailed log
            #   ${OUTDIR}/DC1L1_SJ.out.tab                     - Splice junctions
        else
            echo "  ✗ ${base} failed (exit code: $?)"
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            echo "  Continuing with next sample..."
        fi
        
        echo "  ------------------------------------------------------------"
    done
    
    echo ""
    echo "Completed directory: ${SAMPLE_DIR}"
    echo "============================================================"
done

# ===== Final Summary =====
echo ""
echo "╔════════════════════════════════════════════════════════════"
echo "║ FINAL SUMMARY"
echo "╚════════════════════════════════════════════════════════════"
echo "Total samples processed: ${TOTAL_SAMPLES}"
echo "Successful: ${SUCCESSFUL_SAMPLES}"
echo "Failed: ${FAILED_SAMPLES}"
echo ""
echo "Data source: auto-detected per species (clean if available, else raw)"
echo ""
echo "Genome mapping used:"
echo "  - 00_1_DC → carrot (with GTF)"
echo "  - 00_2_DG → carrot (with GTF)"
echo "  - 00_3_MF → nutmeg (no GTF)"
echo "End time: $(date)"
echo "============================================================"

# List all output directories
echo ""
echo "Output directories created:"
ls -ld "${ALIGNMENT_DIR}"/00_*/ 2>/dev/null || echo "(none)"

exit 0
