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
#   sbatch scripts/04_rsem/run_rsem.sbatch MF nofungi   # cleaned FASTA, new folder
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
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch <SPECIES_CODE> [cdhit|pooled|nofungi]"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF cdhit"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF pooled"
    echo "  sbatch scripts/04_rsem/run_rsem.sbatch MF nofungi"
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

apply_assembly_adjustment "${ASSEMBLY_TYPE}"

case "${ASSEMBLY_TYPE}" in
    cdhit|CDHIT)
        ASSEMBLY="${CDHIT_FASTA}"
        ASSEMBLY_LABEL="CD-HIT (95% identity)"
        ;;
    pooled|POOL)
        ASSEMBLY="${TRINITY_DIR}/${SPECIES_CODE_UPPER}_trinity_pooled.Trinity.fasta"
        ASSEMBLY_LABEL="pooled Trinity (no CD-HIT)"
        ;;
    nofungi|NOFUNGI)
        # Cleaned FASTA. Counts go to 00_7_RSEM_nofungi, not 00_7_RSEM.
        ASSEMBLY="${NOFUNGI_FASTA}"
        ASSEMBLY_LABEL="CD-HIT after fungal filter"
        ;;
    *)
        echo "ERROR: Unknown assembly type '${ASSEMBLY_TYPE}'."
        echo "Use 'cdhit' (default), 'pooled', or 'nofungi'."
        exit 1
        ;;
esac


# =============================================================================
# OUTPUT PATHS
# =============================================================================
# nofungi uses RSEM_DIR from config.sh (00_7_RSEM_nofungi).
# cdhit and pooled keep writing to the original 00_7_RSEM folder.
if [[ "${ASSEMBLY_TYPE}" == "nofungi" || "${ASSEMBLY_TYPE}" == "NOFUNGI" ]]; then
    OUT_DIR="${RSEM_DIR}"
else
    OUT_DIR="${PROCESSED_DIR}/00_7_RSEM"
fi
RSEM_REF_DIR="${OUT_DIR}/${SPECIES_CODE_UPPER}_rsem_ref"
RSEM_REF="${RSEM_REF_DIR}/${SPECIES_CODE_UPPER}_rsem"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${RSEM_REF_DIR}"
mkdir -p "${OUT_DIR}"


# =============================================================================
# SET UP RSEM + BOWTIE2 (rnaseq conda env — same as STAR / featureCounts)
#
# On ARC you MUST module purge BEFORE conda activate. Otherwise HPC modules
# stay on PATH and you get mismatched tools, e.g.:
#   rsem    → /apps/arch/software/RSEM/1.3.3-foss-2023a/...
#   bowtie2 → /apps/arch/software/Bowtie2/2.5.4-GCC-13.3.0/...  ← Perl mismatch
#
# Correct conda paths look like:
#   ~/miniconda3/envs/rnaseq/bin/rsem-calculate-expression
#   ~/miniconda3/envs/rnaseq/bin/bowtie2
#
# One-time install on login node (if missing):
#   module purge && conda activate rnaseq
#   conda install -y -c bioconda -c conda-forge rsem bowtie2
# =============================================================================
echo "Clearing inherited HPC modules..."
module purge 2>/dev/null || true

echo "Activating conda environment: ${CONDA_ENV}..."
if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: conda.sh not found under ~/miniconda3 or ~/anaconda3"
    exit 1
fi

conda activate "${CONDA_ENV}"

if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "ERROR: conda env '${CONDA_ENV}' did not activate."
    exit 1
fi

export PATH="${CONDA_PREFIX}/bin:${PATH}"
hash -r

RSEM_BIN="${CONDA_PREFIX}/bin/rsem-calculate-expression"
BOWTIE2_BIN="${CONDA_PREFIX}/bin/bowtie2"

if [[ ! -x "${RSEM_BIN}" ]] || [[ ! -x "${BOWTIE2_BIN}" ]]; then
    echo "ERROR: RSEM or Bowtie2 not found in ${CONDA_PREFIX}/bin/"
    echo "On a login node run:"
    echo "  module purge && conda activate ${CONDA_ENV}"
    echo "  conda install -y -c bioconda -c conda-forge rsem bowtie2"
    exit 1
fi

echo "Conda env:       ${CONDA_DEFAULT_ENV}"
echo "RSEM:            ${RSEM_BIN}"
echo "Bowtie2:         ${BOWTIE2_BIN}"

if ! "${RSEM_BIN}" --version &>/dev/null; then
    echo "ERROR: RSEM failed — check for module contamination (use module purge)."
    exit 1
fi

echo "RSEM version:    $("${RSEM_BIN}" --version 2>&1)"
echo "Bowtie2 version: $("${BOWTIE2_BIN}" --version 2>&1 | head -1)"
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
    elif [[ "${ASSEMBLY_TYPE}" == "nofungi" || "${ASSEMBLY_TYPE}" == "NOFUNGI" ]]; then
        echo "The cleaned FASTA is missing. It is written by the fungal filter:"
        echo "  sbatch scripts/01_qc/run_minimap2_fungi.sbatch"
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

# The Trinity gene map still lists fungal transcripts that were removed.
# RSEM will error if the map names a transcript that is not in the FASTA.
# Write a smaller map in the new folder. The original Trinity map stays put.
if [[ "${ASSEMBLY_TYPE}" == "nofungi" || "${ASSEMBLY_TYPE}" == "NOFUNGI" ]]; then
    FILTERED_GENE_MAP="${OUT_DIR}/${SPECIES_CODE_UPPER}_nofungi.gene_trans_map"
    echo ""
    echo "Keeping gene-map rows whose transcript is still in the cleaned FASTA..."
    python3 - "${ASSEMBLY}" "${GENE_MAP}" "${FILTERED_GENE_MAP}" << 'PY'
import sys

fasta_path = sys.argv[1]
map_path = sys.argv[2]
out_path = sys.argv[3]

# Transcript names still present after the fungal filter.
keep_names = set()
fasta_handle = open(fasta_path, "r")
for line in fasta_handle:
    if line.startswith(">"):
        name = line[1:].strip().split()[0]
        keep_names.add(name)
fasta_handle.close()

kept = 0
map_handle = open(map_path, "r")
out_handle = open(out_path, "w")
for line in map_handle:
    parts = line.rstrip("\n").split("\t")
    # Trinity map columns: gene id, then transcript id.
    if len(parts) < 2:
        continue
    transcript_id = parts[1]
    if transcript_id in keep_names:
        out_handle.write(line)
        kept = kept + 1
map_handle.close()
out_handle.close()
print("Gene-map rows kept: " + str(kept))
PY
    if [[ ! -s "${FILTERED_GENE_MAP}" ]]; then
        echo "ERROR: Filtered gene map is empty: ${FILTERED_GENE_MAP}"
        exit 1
    fi
    GENE_MAP="${FILTERED_GENE_MAP}"
    echo "Using filtered gene map: ${GENE_MAP}"
fi
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

    "${CONDA_PREFIX}/bin/rsem-prepare-reference" \
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
#   (no --output-genome-bam — Trinity is transcriptome-only, no genome FASTA)
#   RSEM still writes .transcript.bam by default for IGV
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

    "${RSEM_BIN}" \
        --paired-end \
        --bowtie2 \
        --strandedness none \
        --estimate-rspd \
        --append-names \
        --num-threads "${THREADS}" \
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
        # RSEM .cnt line 1 (skip # comments): N0 N1 N2 N_tot
        #   N0=unalignable  N1=alignable  N2=filtered  N_tot=total
        awk '
            $1 ~ /^#/ { next }
            NF >= 4 && total == 0 {
                unalignable = $1 + 0
                aligned     = $2 + 0
                filtered    = $3 + 0
                total       = $4 + 0
                if (total > 0) {
                    pct = (aligned / total) * 100
                    printf "  Total reads:     %d\n", total
                    printf "  Alignable (N1):  %d\n", aligned
                    printf "  Unalignable:     %d\n", unalignable
                    printf "  Filtered (N2):   %d\n", filtered
                    printf "  Alignment rate:  %.1f%%\n", pct
                } else {
                    print "  (could not parse .cnt — check file manually)"
                }
                exit
            }
        ' total=0 "${STAT_FILE}"
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

# Build the gene matrix from each sample's .genes.results file.
# RSEM already summed isoforms into genes in those files.
# Do not pass them through abundance_estimates_to_matrix.pl with the
# gene map: that tool looks up transcript names (..._i1). The gene file's
# first column is the gene name (no _i1), so every row is skipped and
# RSEM.gene.counts.matrix is left with only the sample header.
python3 - "${OUT_DIR}" "${SAMPLE_NAMES[@]}" << 'PY'
import sys

out_dir = sys.argv[1]
sample_names = sys.argv[2:]

# gene_id -> list of expected counts, one per sample, in sample_names order
counts_by_gene = {}
sample_lengths = []

for sample_name in sample_names:
    path = out_dir + "/" + sample_name + "/" + sample_name + ".genes.results"
    handle = open(path, "r")
    header = handle.readline().rstrip("\n").split("\t")
    gene_col = header.index("gene_id")
    count_col = header.index("expected_count")
    n_genes = 0
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= count_col:
            continue
        gene_id = parts[gene_col]
        expected_count = parts[count_col]
        if gene_id not in counts_by_gene:
            counts_by_gene[gene_id] = {}
        counts_by_gene[gene_id][sample_name] = expected_count
        n_genes = n_genes + 1
    handle.close()
    sample_lengths.append(n_genes)
    print(sample_name + " genes: " + str(n_genes))

if len(counts_by_gene) == 0:
    sys.stderr.write("ERROR: No gene rows were read from the .genes.results files.\n")
    sys.exit(1)

out_path = out_dir + "/RSEM.gene.counts.matrix"
out_handle = open(out_path, "w")
# First cell is blank so column 1 is the gene id, matching the old matrix.
header_line = "\t" + "\t".join(sample_names) + "\n"
out_handle.write(header_line)
for gene_id in sorted(counts_by_gene):
    values = []
    for sample_name in sample_names:
        values.append(counts_by_gene[gene_id].get(sample_name, "0"))
    out_handle.write(gene_id + "\t" + "\t".join(values) + "\n")
out_handle.close()
print("Wrote " + str(len(counts_by_gene)) + " genes to " + out_path)
PY

echo "✓ Count matrix created in ${MATRIX_DIR}"
ls -lh RSEM.gene.counts.matrix
# A header-only file is about 31 bytes. A real matrix is many megabytes.
MATRIX_BYTES=$(wc -c < RSEM.gene.counts.matrix | tr -d ' ')
if [[ "${MATRIX_BYTES}" -lt 1000 ]]; then
    echo "ERROR: RSEM.gene.counts.matrix is only ${MATRIX_BYTES} bytes."
    echo "The per-sample .genes.results files did not contain gene rows."
    exit 1
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
if [[ "${ASSEMBLY_TYPE}" == "nofungi" || "${ASSEMBLY_TYPE}" == "NOFUNGI" ]]; then
    NEXT_TAG=" nofungi"
else
    NEXT_TAG=""
fi
echo "Next steps:"
echo "  1. Check alignment rates above (target >70% per sample)"
echo "  2. Build PyDESeq2 inputs:"
echo "     sbatch scripts/04_counting/run_count_matrix.sbatch ${SPECIES_CODE_UPPER} --type rsem${NEXT_TAG}"
echo "  3. Run PyDESeq2 (MF example — Fruit vs Leaf):"
echo "     CONTRAST_A=F CONTRAST_B=L sbatch scripts/05_pydeseq2/run_step1_analysis.sbatch ${SPECIES_CODE_UPPER}${NEXT_TAG}"
echo ""
echo "Finished: $(date)"
echo "========================================================"

exit 0
