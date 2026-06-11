#!/bin/bash
# =============================================================================
# SCRIPT:  run_rsem.sbatch
# PURPOSE: Quantify expression for each sample using RSEM + Bowtie2 against
#          a pooled Trinity assembly (CD-HIT filtered by default).
#
# WHAT RSEM DOES:
#   1. Aligns trimmed reads back to the Trinity assembly (Bowtie2 internally)
#   2. Uses the EM algorithm to handle multi-mapping reads probabilistically
#   3. Produces expected read counts per transcript AND per gene
#   4. Gene-level counts go into DESeq2 / PyDESeq2
#
# INPUT (MF example):
#   Reference: MF_trinity_cdhit95.fasta (CD-HIT output — default)
#   Gene map:  MF_trinity_pooled.Trinity.fasta.gene_trans_map
#   Reads:     6 samples × R1 + R2 (trimmed / Kraken-filtered)
#
# OUTPUT per sample (in 01_processed/00_7_RSEM/<SAMPLE>/):
#   <SAMPLE>.genes.results      ← USE THIS for DESeq2
#   <SAMPLE>.isoforms.results   ← isoform level (not for DESeq2)
#   <SAMPLE>.transcript.bam     ← alignment file
#
# HOW TO SUBMIT:
#   sbatch scripts/04_rsem/run_rsem.sbatch MF
#   sbatch scripts/04_rsem/run_rsem.sbatch MF cdhit     # same as default
#   sbatch scripts/04_rsem/run_rsem.sbatch MF pooled    # use raw pooled Trinity
#
# AFTER THIS JOB:
#   RSEM.gene.counts.matrix is written in 01_processed/00_7_RSEM/
#   → feed that matrix into scripts/05_pydeseq2/
#
# PATHS: All folder paths come from scripts/config.sh
# =============================================================================


# ---- SLURM cluster resource requests ----
# RSEM + Bowtie2 alignment is CPU-heavy; 16 cores and 64 GB RAM is a safe default.
#SBATCH --job-name=rsem
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
# NOTE: SBATCH directives are parsed before bash runs — path must be hardcoded.
# If you change BASE_DIR in config.sh, also update --chdir below.
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed/00_7_RSEM
#SBATCH -o rsem_%j.out
#SBATCH -e rsem_%j.err

set -euo pipefail


# =============================================================================
# LOAD CONFIG
# =============================================================================
CONFIG_FILE=""
for candidate in \
    "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../config.sh" \
    "${SLURM_SUBMIT_DIR:-.}/../../config.sh" \
    "/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code/scripts/config.sh"
do
    if [[ -f "${candidate}" ]]; then
        CONFIG_FILE="${candidate}"
        break
    fi
done

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "ERROR: Could not locate scripts/config.sh."
    echo "Submit from inside the rnaseq repo, e.g.:"
    echo "  cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF"
    exit 1
fi

echo "Loading config from: ${CONFIG_FILE}"
source "${CONFIG_FILE}"


# =============================================================================
# READ SPECIES CODE + ASSEMBLY TYPE FROM COMMAND LINE
# =============================================================================
SPECIES_CODE="${1:-}"
ASSEMBLY_TYPE="${2:-cdhit}"

if [[ -z "${SPECIES_CODE}" ]]; then
    echo "ERROR: No species code provided."
    echo ""
    echo "Usage:"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch <SPECIES_CODE> [cdhit|pooled]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF cdhit"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF pooled"
    exit 1
fi

get_sample_info "${SPECIES_CODE}"
SPECIES_CODE_UPPER=$(echo "${SPECIES_CODE}" | tr '[:lower:]' '[:upper:]')


# =============================================================================
# PICK ASSEMBLY + GENE MAP PATHS (same Trinity folders as run_busco.sbatch)
# =============================================================================
case "${SPECIES_CODE_UPPER}" in
    MF)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_MF_trinity"
        GENE_MAP="${TRINITY_DIR}/MF_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DC)
        TRINITY_DIR="${PROCESSED_DIR}/00_1_DC_trinity"
        GENE_MAP="${TRINITY_DIR}/DC_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DG)
        TRINITY_DIR="${PROCESSED_DIR}/00_2_DG_trinity"
        GENE_MAP="${TRINITY_DIR}/DG_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    SK)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_Sukman_trinity"
        GENE_MAP="${TRINITY_DIR}/SK_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    DCDG)
        TRINITY_DIR="${PROCESSED_DIR}/00_4_DC_DG_trinity"
        GENE_MAP="${TRINITY_DIR}/DCDG_trinity_pooled.Trinity.fasta.gene_trans_map"
        ;;
    *)
        echo "ERROR: No Trinity paths configured for species '${SPECIES_CODE_UPPER}'."
        echo "Add a new case block to run_rsem.sbatch (same as run_busco.sbatch)."
        exit 1
        ;;
esac

case "${ASSEMBLY_TYPE}" in
    cdhit|CDHIT)
        ASSEMBLY="${PROCESSED_DIR}/00_6_cdhit/${SPECIES_CODE_UPPER}_trinity_cdhit95.fasta"
        ASSEMBLY_LABEL="CD-HIT (95% identity)"
        ;;
    pooled|POOL)
        ASSEMBLY="${TRINITY_DIR}/${SPECIES_CODE_UPPER}_trinity_pooled.Trinity.fasta"
        ASSEMBLY_LABEL="pooled Trinity (no CD-HIT)"
        ;;
    *)
        echo "ERROR: Unknown assembly type '${ASSEMBLY_TYPE}'."
        echo "Use 'cdhit' (default) or 'pooled'."
        exit 1
        ;;
esac


# =============================================================================
# OUTPUT PATHS
# =============================================================================
OUT_DIR="${PROCESSED_DIR}/00_7_RSEM"
RSEM_REF_DIR="${OUT_DIR}/${SPECIES_CODE_UPPER}_rsem_ref"
RSEM_REF="${RSEM_REF_DIR}/${SPECIES_CODE_UPPER}_rsem"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${RSEM_REF_DIR}"
mkdir -p "${OUT_DIR}"


# =============================================================================
# ACTIVATE CONDA ENVIRONMENT
#
# RSEM is NOT available as an ARC module on owl1 (module spider RSEM returns
# nothing). Use the rnaseq conda env instead — environment.yml lists rsem,
# bowtie2, and trinity together so Perl/GCC versions stay compatible.
#
# Do NOT use "module load RSEM" + "module load Bowtie2" — that caused:
#   Cwd.c: loadable library and perl binaries are mismatched
# =============================================================================
echo "Activating conda environment: ${CONDA_ENV}..."
source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "ERROR: conda env '${CONDA_ENV}' did not activate."
    echo "Create it on a login node: conda env create -f environment.yml"
    exit 1
fi

# Install RSEM + Bowtie2 into rnaseq env if missing (one-time, needs internet)
if ! command -v rsem-calculate-expression &>/dev/null; then
    echo "RSEM not found — installing into ${CONDA_ENV}..."
    conda install -y -c bioconda -c conda-forge rsem bowtie2
fi

if ! command -v rsem-calculate-expression &>/dev/null; then
    echo "ERROR: rsem-calculate-expression still not found."
    echo "On a login node run:"
    echo "  conda activate ${CONDA_ENV}"
    echo "  conda install -y -c bioconda -c conda-forge rsem bowtie2"
    exit 1
fi

if ! command -v bowtie2 &>/dev/null; then
    echo "ERROR: bowtie2 not found in ${CONDA_ENV}."
    echo "  conda install -y -c bioconda bowtie2"
    exit 1
fi

echo "Conda env:       ${CONDA_DEFAULT_ENV}"
echo "RSEM version:    $(rsem-calculate-expression --version 2>&1)"
echo "Bowtie2 version: $(bowtie2 --version | head -1)"
echo "------------------------------------------------------------"


# =============================================================================
# JOB HEADER
# =============================================================================
echo "========================================================"
echo "RSEM expression quantification"
echo "Species code:   ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Assembly type:  ${ASSEMBLY_LABEL}"
echo "Job ID:         ${SLURM_JOB_ID:-N/A}"
echo "Host:           $(hostname)"
echo "Started:        $(date)"
echo "CPUs:           ${THREADS}"
echo "========================================================"


# =============================================================================
# VERIFY CRITICAL INPUT FILES EXIST
# =============================================================================
echo ""
echo "Checking input files..."

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found: ${ASSEMBLY}"
    if [[ "${ASSEMBLY_TYPE}" == "cdhit" ]] || [[ "${ASSEMBLY_TYPE}" == "CDHIT" ]]; then
        echo "Run CD-HIT first:"
        echo "  sbatch scripts/04_busco/run_cdhit.sbatch ${SPECIES_CODE_UPPER}"
    else
        echo "Run pooled Trinity first:"
        echo "  sbatch scripts/03_assembly/trinity_pooled.sbatch"
    fi
    exit 1
fi
echo "✓ Assembly: $(grep -c '>' "${ASSEMBLY}") transcripts"
echo "  ${ASSEMBLY}"

if [[ ! -f "${GENE_MAP}" ]]; then
    echo "ERROR: Gene-transcript map not found: ${GENE_MAP}"
    echo "  This file is created automatically by Trinity."
    echo "  It must come from the SAME pooled Trinity run as your assembly."
    exit 1
fi
echo "✓ Gene map: $(wc -l < "${GENE_MAP}") entries"
echo "  ${GENE_MAP}"
echo "------------------------------------------------------------"


# =============================================================================
# STEP 1: BUILD RSEM REFERENCE (one-time per assembly)
#
# --transcript-to-gene-map: tells RSEM how to aggregate isoform counts
#   into gene-level counts (required for DESeq2).
# --bowtie2: use Bowtie2 as the aligner (recommended for transcriptomes).
# =============================================================================
if [[ -f "${RSEM_REF}.seq" ]]; then
    echo ""
    echo "RSEM reference already exists — skipping build step"
    echo "  ${RSEM_REF}.seq"
else
    echo ""
    echo "========================================================"
    echo "STEP 1: Building RSEM reference index (Bowtie2)..."
    echo "========================================================"

    rsem-prepare-reference \
        --transcript-to-gene-map "${GENE_MAP}" \
        --bowtie2 \
        "${ASSEMBLY}" \
        "${RSEM_REF}"

    echo "✓ Reference built: ${RSEM_REF}"
fi
echo "------------------------------------------------------------"


# =============================================================================
# STEP 2: DEFINE SAMPLES AND READ PATHS
#
# CRITICAL: R1 and R2 arrays must be in IDENTICAL ORDER.
# The Nth R1 file is paired with the Nth R2 file.
#
# MF is fully configured below (matches trinity_pooled.sbatch):
#   MFF = fruit → Kraken2-filtered reads in 00_00_kraken/
#   MFL = leaf  → fastp-cleaned reads in 00_3_MF/
#
# For other species, add a new case block with your sample names and read paths.
# =============================================================================
SAMPLE_NAMES=()
R1_READS=()
R2_READS=()

case "${SPECIES_CODE_UPPER}" in
    MF)
        READS_MFF="${PROCESSED_DIR}/00_00_kraken"
        READS_MFL="${PROCESSED_DIR}/00_3_MF"

        SAMPLE_NAMES=(MFF1 MFF2 MFF3 MFL1 MFL2 MFL3)

        R1_READS=(
            "${READS_MFF}/MFF1_1_filtered.fq.gz"
            "${READS_MFF}/MFF2_1_filtered.fq.gz"
            "${READS_MFF}/MFF3_1_filtered.fq.gz"
            "${READS_MFL}/MFL1_1_clean.fq.gz"
            "${READS_MFL}/MFL2_1_clean.fq.gz"
            "${READS_MFL}/MFL3_1_clean.fq.gz"
        )

        R2_READS=(
            "${READS_MFF}/MFF1_2_filtered.fq.gz"
            "${READS_MFF}/MFF2_2_filtered.fq.gz"
            "${READS_MFF}/MFF3_2_filtered.fq.gz"
            "${READS_MFL}/MFL1_2_clean.fq.gz"
            "${READS_MFL}/MFL2_2_clean.fq.gz"
            "${READS_MFL}/MFL3_2_clean.fq.gz"
        )
        ;;
    *)
        echo "ERROR: Sample/read paths not configured for '${SPECIES_CODE_UPPER}'."
        echo "Add a SAMPLE_NAMES + R1/R2 block in run_rsem.sbatch for this species."
        exit 1
        ;;
esac

NUM_SAMPLES=${#SAMPLE_NAMES[@]}


# =============================================================================
# STEP 3: QUANTIFY EACH SAMPLE
#
# RSEM flags:
#   --paired-end          reads are paired (R1 + R2)
#   --bowtie2             Bowtie2 alignment (handles mismatches better than Bowtie1)
#   --strandedness none   non-directional library (Novogene MF libraries)
#   --estimate-rspd       correct read-start position bias (recommended)
#   --append-names        add gene/transcript names to output files
#   --output-genome-bam   save alignment BAM for IGV
# =============================================================================
echo ""
echo "========================================================"
echo "STEP 2: Quantifying expression for ${NUM_SAMPLES} samples..."
echo "========================================================"

for i in "${!SAMPLE_NAMES[@]}"; do

    SAMPLE="${SAMPLE_NAMES[$i]}"
    R1="${R1_READS[$i]}"
    R2="${R2_READS[$i]}"
    SAMPLE_OUT="${OUT_DIR}/${SAMPLE}"

    mkdir -p "${SAMPLE_OUT}"

    echo ""
    echo "----------------------------------------"
    echo "Processing sample: ${SAMPLE} ($((i + 1))/${NUM_SAMPLES})"
    echo "  R1: $(basename "${R1}")"
    echo "  R2: $(basename "${R2}")"
    echo "  Output: ${SAMPLE_OUT}/"
    echo "----------------------------------------"

    if [[ ! -f "${R1}" ]]; then
        echo "ERROR: R1 not found: ${R1}"
        exit 1
    fi
    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: R2 not found: ${R2}"
        exit 1
    fi

    GENES_RESULT="${SAMPLE_OUT}/${SAMPLE}.genes.results"
    if [[ -f "${GENES_RESULT}" ]]; then
        echo "✓ Already complete — skipping ${SAMPLE}"
        echo "  Found: ${GENES_RESULT}"
        continue
    fi

    rsem-calculate-expression \
        --paired-end \
        --bowtie2 \
        --strandedness none \
        --estimate-rspd \
        --append-names \
        --num-threads "${THREADS}" \
        --output-genome-bam \
        "${R1}" \
        "${R2}" \
        "${RSEM_REF}" \
        "${SAMPLE_OUT}/${SAMPLE}"

    echo "✓ Completed: ${SAMPLE}"
    echo "  Gene results:    ${SAMPLE_OUT}/${SAMPLE}.genes.results"
    echo "  Isoform results: ${SAMPLE_OUT}/${SAMPLE}.isoforms.results"

done

echo ""
echo "========================================================"
echo "All ${NUM_SAMPLES} samples quantified"
echo "========================================================"


# =============================================================================
# STEP 4: ALIGNMENT STATISTICS SUMMARY
# Good alignment rate: >70%    Concerning: <60%
# =============================================================================
echo ""
echo "========================================================"
echo "Alignment Statistics Summary"
echo "========================================================"

for SAMPLE in "${SAMPLE_NAMES[@]}"; do
    STAT_FILE="${OUT_DIR}/${SAMPLE}/${SAMPLE}.stat/${SAMPLE}.cnt"
    if [[ -f "${STAT_FILE}" ]]; then
        echo "--- ${SAMPLE} ---"
        awk 'NR==1 {
            total=$3
            aligned=$2
            pct = (aligned/total)*100
            printf "  Total reads:    %d\n", total
            printf "  Aligned reads:  %d\n", aligned
            printf "  Alignment rate: %.1f%%\n", pct
        }' "${STAT_FILE}"
    fi
done


# =============================================================================
# STEP 5: COMBINE INTO COUNT MATRIX (Trinity utility script)
# =============================================================================
echo ""
echo "========================================================"
echo "STEP 3: Building count matrix for DESeq2..."
echo "========================================================"

MATRIX_DIR="${OUT_DIR}"
cd "${MATRIX_DIR}"

GENES_RESULT_ARGS=()
for SAMPLE in "${SAMPLE_NAMES[@]}"; do
    GENES_RESULT_ARGS+=("${OUT_DIR}/${SAMPLE}/${SAMPLE}.genes.results")
done

# Trinity ships abundance_estimates_to_matrix.pl (included in rnaseq conda env)
if ! command -v abundance_estimates_to_matrix.pl &>/dev/null; then
    echo "Trinity utility not found — installing trinity into ${CONDA_ENV}..."
    conda install -y -c bioconda -c conda-forge trinity
fi

if ! command -v abundance_estimates_to_matrix.pl &>/dev/null; then
    echo "WARNING: abundance_estimates_to_matrix.pl not in PATH"
    echo "Install Trinity in conda: conda install -y -c bioconda trinity"
    echo "Then run manually from ${MATRIX_DIR}:"
    echo ""
    echo "abundance_estimates_to_matrix.pl \\"
    echo "    --est_method RSEM \\"
    echo "    --gene_trans_map ${GENE_MAP} \\"
    echo "    --name_sample_by_last_column_name \\"
    echo "    --cross_sample_norm TMM \\"
    for SAMPLE in "${SAMPLE_NAMES[@]}"; do
        echo "    ${OUT_DIR}/${SAMPLE}/${SAMPLE}.genes.results \\"
    done
    echo ""
    echo "Skipping matrix creation — run manually after loading Trinity"
else
    abundance_estimates_to_matrix.pl \
        --est_method RSEM \
        --gene_trans_map "${GENE_MAP}" \
        --name_sample_by_last_column_name \
        --cross_sample_norm TMM \
        "${GENES_RESULT_ARGS[@]}"

    echo "✓ Count matrix created in ${MATRIX_DIR}"
    ls -lh RSEM.gene.counts.matrix 2>/dev/null || \
        echo "  Check ${MATRIX_DIR} for RSEM.gene.counts.matrix"
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ RSEM COMPLETE"
echo "========================================================"
echo "Species:        ${SPECIES_CODE_UPPER} (${SPECIES_NAME})"
echo "Assembly:       ${ASSEMBLY}"
echo "Samples:        ${SAMPLE_NAMES[*]}"
echo "Results in:     ${OUT_DIR}/"
echo ""
echo "For DESeq2 / PyDESeq2, use:"
echo "  ${MATRIX_DIR}/RSEM.gene.counts.matrix"
echo ""
echo "Next steps:"
echo "  1. Check alignment rates above (target >70% per sample)"
echo "  2. Run PyDESeq2: scripts/05_pydeseq2/"
echo "  3. Run annotation (TransDecoder / BLAST)"
echo ""
echo "Finished: $(date)"
echo "========================================================"

exit 0
