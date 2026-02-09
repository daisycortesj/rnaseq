#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Align RNA-seq reads for ALL sample directories using STAR
# USAGE:   sbatch run_alignment_all.sbatch
# NOTE:    This processes ALL sample directories with auto-mapped genomes:
#          - 00_1_DC → carrot genome
#          - 00_2_DG → carrot genome
#          - 00_3_MF → nutmeg genome
#          For a single directory, use run_alignment.sbatch
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
#SBATCH -o alignment_all_%j.out  
#SBATCH -e alignment_all_%j.err       

set -euo pipefail

# ===== Configuration =====
# ┌─────────────────────────────────────────────────────────────┐
# │ IMPORTANT: Choose your input data source                   │
# │ Set INPUT_DATA to:                                         │
# │   "raw"   → Use original FASTQ files (if QC is good)      │
# │   "fastp" → Use cleaned FASTQ files (if QC shows problems)│
# └─────────────────────────────────────────────────────────────┘

INPUT_DATA="raw"           # Change to "fastp" after cleaning

# Directory paths
RAWDATA_DIR="/projects/tholl_lab_1/daisy_analysis/00_rawdata"
FASTP_DIR="/projects/tholl_lab_1/daisy_analysis/01_processed/fastp"
REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
ALIGNMENT_DIR="/projects/tholl_lab_1/daisy_analysis/02_mapped"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Determine input directory based on choice
if [[ "${INPUT_DATA}" == "fastp" ]]; then
    INPUT_DIR="${FASTP_DIR}"
    INPUT_SOURCE="FASTP-CLEANED"
    FILE_SUFFIX="_clean"
else
    INPUT_DIR="${RAWDATA_DIR}"
    INPUT_SOURCE="RAW"
    FILE_SUFFIX=""
fi

echo "============================================================"
echo "STAR Read Alignment Pipeline - Process All Samples"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo ""
echo "┌─────────────────────────────────────────────────────────┐"
echo "│ DATA SOURCE: ${INPUT_SOURCE}"
echo "└─────────────────────────────────────────────────────────┘"
echo "  Input: ${INPUT_DIR}"
echo ""
echo "Genome mapping:"
echo "  - 00_1_DC → carrot"
echo "  - 00_2_DG → carrot"
echo "  - 00_3_MF → nutmeg"
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
echo "Searching for sample directories in ${INPUT_DIR}..."
echo "============================================================"

TOTAL_SAMPLES=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

# ====================================================================
# MAIN PROCESSING LOOP: Process all sample directories
# ====================================================================
# This loop processes each sample directory (00_1_DC, 00_2_DG, 00_3_MF)
# For each directory, it:
#   1. Determines which genome to use (carrot or nutmeg)
#   2. Finds all R1 files in the directory
#   3. Pairs each R1 with its matching R2 file
#   4. Runs STAR alignment for each sample pair
# ====================================================================
for sample_dir_path in "${INPUT_DIR}"/00_*/ ; do
    # Extract directory name (e.g., 00_1_DC)
    SAMPLE_DIR=$(basename "${sample_dir_path}")
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════"
    echo "║ Processing directory: ${SAMPLE_DIR}"
    echo "╚════════════════════════════════════════════════════════════"
    
    # Determine genome based on directory name
    case "${SAMPLE_DIR}" in
        00_1_DC|00_2_DG)
            GENOME_TYPE="carrot"
            GENOME_INDEX="${REFERENCE_DIR}/carrot_star_index"
            QUANT_MODE="--quant-mode"
            echo "Using CARROT genome (with gene quantification)"
            ;;
        00_3_MF)
            GENOME_TYPE="nutmeg"
            GENOME_INDEX="${REFERENCE_DIR}/nutmeg_star_index"
            QUANT_MODE="--no-quant-mode"
            echo "Using NUTMEG genome (MYU assembly)"
            ;;
        *)
            echo "⚠ WARNING: Unknown directory '${SAMPLE_DIR}', skipping..."
            continue
            ;;
    esac
    
    # Validate genome index exists
    if [[ ! -d "${GENOME_INDEX}" ]]; then
        echo "✗ ERROR: Genome index not found: ${GENOME_INDEX}"
        echo "Please build the index first using run_genome_index.sbatch"
        continue
    fi
    echo "✓ Genome index: ${GENOME_INDEX}"
    
    # Set output directory for this sample directory
    OUTDIR="${ALIGNMENT_DIR}/${SAMPLE_DIR}"
    echo "Output directory: ${OUTDIR}"
    
    # Create output directory
    mkdir -p "${OUTDIR}"
    
    # Count samples in this directory (adjust pattern based on data source)
    if [[ "${INPUT_DATA}" == "fastp" ]]; then
        SAMPLE_COUNT=$(find "${sample_dir_path}" -name "*_1_clean.fq.gz" 2>/dev/null | wc -l)
    else
        SAMPLE_COUNT=$(find "${sample_dir_path}" -name "*_1.fq.gz" 2>/dev/null | wc -l)
    fi
    
    if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
        echo "⚠ No samples found in ${SAMPLE_DIR}, skipping..."
        continue
    fi
    
    echo "Found ${SAMPLE_COUNT} sample(s) in ${SAMPLE_DIR}"
    echo "------------------------------------------------------------"
    
    # ================================================================
    # SAMPLE EXTRACTION AND PAIRING SECTION
    # ================================================================
    
    # STEP 1: Determine which file pattern to search for
    # This depends on whether we're using raw or cleaned data
    if [[ "${INPUT_DATA}" == "fastp" ]]; then
        FILE_PATTERN="*_1_clean.fq.gz"    # Pattern for cleaned files
    else
        FILE_PATTERN="*_1.fq.gz"          # Pattern for raw R1 files (forward reads)
    fi
    
    # STEP 2: Loop through all R1 files matching the pattern
    # Example: finds DC1L1_1.fq.gz, DC1L2_1.fq.gz, DC1R1_1.fq.gz, etc.
    for file in "${sample_dir_path}"${FILE_PATTERN} ; do
        
        # STEP 3: Extract sample base name from the file path
        # ------------------------------------------------------------
        withpath="${file}"                          # Full path: /projects/.../00_1_DC/DC1L1_1.fq.gz
        filename="${withpath##*/}"                  # Strips path, keeps filename: DC1L1_1.fq.gz
        base="${filename%_1${FILE_SUFFIX}.fq.gz}"  # Removes suffix, extracts base: DC1L1
        
        # Example walkthrough for DC1L1:
        #   Input file:  /path/to/00_1_DC/DC1L1_1.fq.gz
        #   withpath  =  /path/to/00_1_DC/DC1L1_1.fq.gz
        #   filename  =  DC1L1_1.fq.gz        (removed directory path)
        #   base      =  DC1L1                (removed _1.fq.gz suffix)
        
        # STEP 4: Construct full paths for paired-end reads (R1 + R2)
        # ------------------------------------------------------------
        #LOOKHERE : This is where we pair R1 and R2 files for the same sample
        READS_R1="${sample_dir_path}${base}_1${FILE_SUFFIX}.fq.gz"  # Forward reads: DC1L1_1.fq.gz
        READS_R2="${sample_dir_path}${base}_2${FILE_SUFFIX}.fq.gz"  # Reverse reads: DC1L1_2.fq.gz
        
        # IMPORTANT: R1 and R2 are PAIRED sequencing reads from the SAME biological sample
        #   - R1 (_1.fq.gz) = Forward/left reads (first read in pair)
        #   - R2 (_2.fq.gz) = Reverse/right reads (mate pair)
        #   - They are NOT biological replicates - they're technical pairs from sequencing
        
        # Example for sample DC1L1:
        #   READS_R1 = /path/to/00_1_DC/DC1L1_1.fq.gz
        #   READS_R2 = /path/to/00_1_DC/DC1L1_2.fq.gz
        
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
echo "DATA USED: ${INPUT_SOURCE}"
if [[ "${INPUT_SOURCE}" == "RAW" ]]; then
    echo "  (If alignment rates are low, consider using fastp-cleaned data)"
fi
echo ""
echo "Genome mapping used:"
echo "  - 00_1_DC → carrot"
echo "  - 00_2_DG → carrot"
echo "  - 00_3_MF → nutmeg"
echo "End time: $(date)"
echo "============================================================"

# List all output directories
echo ""
echo "Output directories created:"
ls -ld "${ALIGNMENT_DIR}"/00_*/ 2>/dev/null || echo "(none)"

exit 0
