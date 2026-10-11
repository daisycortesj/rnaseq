#!/bin/bash
# =====================================================================
# SCRIPT:  search_mf30795_reads.sbatch
# PURPOSE: Search two read sets, one after the other, without copying
#          either set to a laptop.
#
#   PART 1 — your reads (fruit + leaf FASTQs that went into Trinity)
#   PART 2 — the Ibdah public run (PRJDB17003, SRA run DRR513061)
#
#   Both parts search the SAME four proteins, then summarize with the
#   SAME rules (scripts/11_verify/summarize_read_hits.py):
#     Mf41180, Mf31663, Mf30795, and SS (sesamin synthase, the control)
#
#   The count is unique fragments (a read pair counts once), kept only
#   when the alignment is at least 30 amino acids, and assigned to the
#   single best protein. That is why the two SUMMARY_fragments.txt
#   files can be compared.
#
# WHY THIS STEP:
#   A missing Trinity hit only says the protein is absent from your
#   assembly. These two searches ask whether each candidate is in your
#   reads, in the Ibdah reads, or in both.
#
# WHICH READS (default):
#   The same 12 files used by scripts/03_assembly/trinity_pooled.sbatch
#     fruit (MFF1-3): 01_processed/00_00_kraken/*_filtered.fq.gz
#     leaf  (MFL1-3): 01_processed/00_3_MF/*_clean.fq.gz
#   These are the reads Trinity actually saw. They are the right files
#   for "did assembly miss a gene that was in its input?"
#
# OPTIONAL:
#   sbatch scripts/11_verify/search_mf30795_reads.sbatch raw
#   searches 00_rawdata/00_3_MF instead (the files from the sequencer,
#   before fastp and Kraken). Use that only if the Trinity-input search
#   finds nothing and you want to know whether QC removed the gene.
#
# HOW TO SUBMIT (on the cluster, from the rnaseq repo):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   conda activate rnaseq
#   conda install -c bioconda diamond sra-tools   # only the first time
#   sbatch scripts/11_verify/search_mf30795_reads.sbatch
#
# OUTPUTS (one pair of files per dataset; the columns match):
#   Part 1: 06_analysis/mf_candidates_read_search_<trinity|raw>/SUMMARY_fragments.txt
#   Part 2: 06_analysis/prjdb17003_DRR513061_search/SUMMARY_fragments.txt
#
# HOW TO READ THE SHARED TABLE:
#   f_95_100   the gene itself
#   f_70_95    a close relative
#   f_35_70    a more distant P450 relative (SS usually lands here)
#   f_below35  a weak match
#   Compare the same protein's f_95_100 column between the two files.
# =====================================================================


# ===== SLURM directives =====
# Do not run diamond on the login node. sbatch sends this to a compute node.
#SBATCH --job-name=mf_candidates
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/mf_candidates_%j.out
#SBATCH -e 06_analysis/mf_candidates_%j.err

# Stop immediately if a command fails, and treat unset variables as errors.
set -euo pipefail


# =====================================================================
# LOAD CONFIG
# =====================================================================
# config.sh holds BASE_DIR, RAWDATA_DIR, PROCESSED_DIR, and ANALYSIS_DIR.
# This finds it whether you submitted from the repo folder or from
# daisy_analysis.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [[ -f "${REPO_ROOT}/scripts/config.sh" ]]; then
    source "${REPO_ROOT}/scripts/config.sh"
elif [[ -f "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh" ]]; then
    source "${SLURM_SUBMIT_DIR}/scripts/config.sh"
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    echo "ERROR: Could not find scripts/config.sh"
    echo "Submit from: /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
    exit 1
fi


# =====================================================================
# CHOOSE WHICH READS TO SEARCH
# =====================================================================
# First argument wins. If you omit it, the job searches the Trinity inputs.
#   trinity  (default)  reads that went into the pooled assembly
#   raw                 original FASTQ files in 00_rawdata/00_3_MF

READ_SET="${1:-trinity}"
READ_SET="$(echo "${READ_SET}" | tr '[:upper:]' '[:lower:]')"

if [[ "${READ_SET}" != "trinity" && "${READ_SET}" != "raw" ]]; then
    echo "ERROR: Unknown read set '${READ_SET}'"
    echo "Use:  sbatch scripts/11_verify/search_mf30795_reads.sbatch"
    echo "Or:   sbatch scripts/11_verify/search_mf30795_reads.sbatch raw"
    exit 1
fi

# One protein file for both searches. If the databases differed, a read
# could be assigned to a different "best" protein in each dataset.
QUERY_FAA="${REPO_ROOT}/reference/queries/ibdah_mf_candidates.faa"
OUT_DIR="${ANALYSIS_DIR}/mf_candidates_read_search_${READ_SET}"
SUMMARIZE_PY="${SCRIPT_DIR}/summarize_read_hits.py"

# Part 2: the Ibdah SRA run. Hits are saved under 06_analysis, not in rawdata.
IBDAH_QUERY="${QUERY_FAA}"
IBDAH_DIR="${RAWDATA_DIR}/00_6_Ibdah"
SRA_FILE="${IBDAH_DIR}/DRR513061/DRR513061.sra"
IBDAH_OUT="${ANALYSIS_DIR}/prjdb17003_DRR513061_search"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Same folders as trinity_pooled.sbatch
INPUT_MFF="${PROCESSED_DIR}/00_00_kraken"
INPUT_MFL="${PROCESSED_DIR}/00_3_MF"
RAW_MF_DIR="${RAWDATA_DIR}/00_3_MF"


# =====================================================================
# BUILD THE LIST OF FASTQ FILES
# =====================================================================

READ_FILES=()

if [[ "${READ_SET}" == "trinity" ]]; then
    # Order matches trinity_pooled.sbatch: fruit first, then leaf.
    # Both mates are searched. A read pair can therefore contribute 2 hits.
    READ_FILES=(
        "${INPUT_MFF}/MFF1_1_filtered.fq.gz"
        "${INPUT_MFF}/MFF1_2_filtered.fq.gz"
        "${INPUT_MFF}/MFF2_1_filtered.fq.gz"
        "${INPUT_MFF}/MFF2_2_filtered.fq.gz"
        "${INPUT_MFF}/MFF3_1_filtered.fq.gz"
        "${INPUT_MFF}/MFF3_2_filtered.fq.gz"
        "${INPUT_MFL}/MFL1_1_clean.fq.gz"
        "${INPUT_MFL}/MFL1_2_clean.fq.gz"
        "${INPUT_MFL}/MFL2_1_clean.fq.gz"
        "${INPUT_MFL}/MFL2_2_clean.fq.gz"
        "${INPUT_MFL}/MFL3_1_clean.fq.gz"
        "${INPUT_MFL}/MFL3_2_clean.fq.gz"
    )
else
    # Raw sequencer files. Names in this project look like MFF1_1.fq.gz.
    while IFS= read -r fastq_path; do
        READ_FILES+=("${fastq_path}")
    done < <(find "${RAW_MF_DIR}" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort)
fi


# =====================================================================
# SAFETY CHECKS
# =====================================================================

if [[ ! -f "${QUERY_FAA}" ]]; then
    echo "ERROR: Query proteins not found: ${QUERY_FAA}"
    echo "Sync the repo so reference/queries/ibdah_mf_candidates.faa is there."
    exit 1
fi

if [[ ! -f "${SUMMARIZE_PY}" ]]; then
    echo "ERROR: Summary script not found: ${SUMMARIZE_PY}"
    echo "Sync the repo so scripts/11_verify/summarize_read_hits.py is there."
    exit 1
fi

# Part 2 needs either the packed SRA file or FASTQ files already unpacked from it.
# Checking now means the job stops before the long Part 1 search if this path is wrong.
IBDAH_READY=0
if [[ -f "${SRA_FILE}" ]]; then
    IBDAH_READY=1
fi
if [[ -f "${IBDAH_DIR}/DRR513061_1.fastq" || -f "${IBDAH_DIR}/DRR513061_1.fastq.gz" ]]; then
    IBDAH_READY=1
fi
if [[ -f "${IBDAH_DIR}/DRR513061/DRR513061_1.fastq" || -f "${IBDAH_DIR}/DRR513061/DRR513061_1.fastq.gz" ]]; then
    IBDAH_READY=1
fi
if [[ "${IBDAH_READY}" -eq 0 ]]; then
    echo "ERROR: Could not find the Ibdah reads."
    echo "Expected the SRA file: ${SRA_FILE}"
    echo "Or unpacked FASTQ files named DRR513061_1.fastq in ${IBDAH_DIR}"
    exit 1
fi

if [[ "${#READ_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found for read set '${READ_SET}'."
    if [[ "${READ_SET}" == "raw" ]]; then
        echo "Looked in: ${RAW_MF_DIR}"
    fi
    exit 1
fi

MISSING=0
for fastq_path in "${READ_FILES[@]}"; do
    if [[ ! -f "${fastq_path}" ]]; then
        echo "ERROR: Missing FASTQ: ${fastq_path}"
        MISSING=$((MISSING + 1))
    fi
done
if [[ "${MISSING}" -gt 0 ]]; then
    echo "Fix the paths above before resubmitting. Trinity inputs are required for the default run."
    exit 1
fi


# =====================================================================
# ACTIVATE CONDA AND CHECK DIAMOND
# =====================================================================

mkdir -p "${OUT_DIR}/tmp"

# set +u lets conda's startup script use variables this job has not set.
# set -u turns that strict check back on after the environment is active.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq
set -u

if ! command -v diamond &> /dev/null; then
    echo "ERROR: diamond is not installed in the rnaseq environment."
    echo "On a login node (not in this job), run:"
    echo "  conda activate rnaseq"
    echo "  conda install -c bioconda diamond"
    echo "Then resubmit this script."
    exit 1
fi

echo "============================================================"
echo "DIAMOND blastx — same four proteins in your reads, then in the Ibdah reads"
echo "============================================================"
echo "Host: $(hostname)  CPUs: ${THREADS}  Start: $(date)"
echo "Read set: ${READ_SET}"
echo "Query:    ${QUERY_FAA}"
echo "Output:   ${OUT_DIR}"
echo "FASTQ files: ${#READ_FILES[@]}"
# diamond version can return a pipe error. Ignore that so the search still starts.
diamond version > "${OUT_DIR}/diamond_version.txt" 2>&1 || true
echo "diamond: $(head -n 1 "${OUT_DIR}/diamond_version.txt")"
echo "============================================================"


# =====================================================================
# BUILD A TINY PROTEIN DATABASE
# =====================================================================
# DIAMOND searches reads (DNA) against proteins, so the database is the
# two proteins, not the reads. Two sequences finish in seconds.
# --dbtype is implied by a protein FASTA for makedb.

diamond makedb \
    --in "${QUERY_FAA}" \
    --db "${OUT_DIR}/queries" \
    --threads "${THREADS}"

echo "Database written: ${OUT_DIR}/queries.dmnd"


# =====================================================================
# SEARCH EACH FASTQ
# =====================================================================
# blastx translates the read in 6 frames and aligns it to the proteins.
# --evalue 1e-5          keep hits that are unlikely to be chance
# --max-target-seqs 1    keep only the best protein hit for each read
# --sensitive            needed so the ~40% SS control is not missed
# --outfmt 6             one tab-separated hit per line:
#                        read_id  protein  %identity  length  evalue  bitscore

for fastq_path in "${READ_FILES[@]}"; do
    sample_name="$(basename "${fastq_path}")"
    sample_name="${sample_name%.gz}"
    sample_name="${sample_name%.fq}"
    sample_name="${sample_name%.fastq}"

    hit_file="${OUT_DIR}/${sample_name}_hits.tsv"

    echo "------------------------------------------------------------"
    echo "Searching: ${fastq_path}"
    echo "Hits file: ${hit_file}"

    diamond blastx \
        --query "${fastq_path}" \
        --db "${OUT_DIR}/queries" \
        --out "${hit_file}" \
        --threads "${THREADS}" \
        --evalue 1e-5 \
        --max-target-seqs 1 \
        --sensitive \
        --tmpdir "${OUT_DIR}/tmp" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore

    # An empty file is a real result (zero hits), not a crash.
    if [[ ! -f "${hit_file}" ]]; then
        : > "${hit_file}"
    fi
    hit_count="$(wc -l < "${hit_file}" | tr -d ' ')"
    echo "Hits in this file: ${hit_count}"
done


# =====================================================================
# SUMMARIZE YOUR READS
# =====================================================================
# Same script and same rules as the Ibdah summary below.
# --tissue auto splits MFF (fruit) from MFL (leaf).

python3 "${SUMMARIZE_PY}" \
    --hits-dir "${OUT_DIR}" \
    --queries "${QUERY_FAA}" \
    --label "Our reads (${READ_SET})"

echo "Part 1 finished: $(date)"
echo "Your-read summary: ${OUT_DIR}/SUMMARY_fragments.txt"


# =====================================================================
# PART 2: SEARCH THE IBDAH PUBLIC RUN (PRJDB17003 / DRR513061)
# =====================================================================
# This is the other lab's nutmeg RNA-seq run. fasterq-dump turns the
# packed .sra file into normal FASTQ files. If those FASTQ files are
# already on disk, the dump is skipped.
#
# Proteins are the same four used in Part 1: Mf41180, Mf31663,
# Mf30795, and SS. The summary script is also the same.

echo "============================================================"
echo "PART 2 — Ibdah run DRR513061"
echo "============================================================"

mkdir -p "${IBDAH_OUT}/tmp"

IBDAH_READS=()
for fastq_path in \
    "${IBDAH_DIR}/DRR513061_1.fastq" \
    "${IBDAH_DIR}/DRR513061_2.fastq" \
    "${IBDAH_DIR}/DRR513061_1.fastq.gz" \
    "${IBDAH_DIR}/DRR513061_2.fastq.gz" \
    "${IBDAH_DIR}/DRR513061/DRR513061_1.fastq" \
    "${IBDAH_DIR}/DRR513061/DRR513061_2.fastq" \
    "${IBDAH_DIR}/DRR513061/DRR513061_1.fastq.gz" \
    "${IBDAH_DIR}/DRR513061/DRR513061_2.fastq.gz"
do
    if [[ -f "${fastq_path}" ]]; then
        IBDAH_READS+=("${fastq_path}")
    fi
done

if [[ "${#IBDAH_READS[@]}" -eq 0 ]]; then
    if ! command -v fasterq-dump &> /dev/null; then
        echo "ERROR: fasterq-dump is not installed, and the FASTQ files are not unpacked yet."
        echo "On a login node, run:"
        echo "  conda activate rnaseq"
        echo "  conda install -c bioconda sra-tools"
        exit 1
    fi
    echo "Unpacking ${SRA_FILE}"
    echo "This can take a while. Temporary files go in ${IBDAH_DIR}/tmp"
    mkdir -p "${IBDAH_DIR}/tmp"
    # --split-files writes read 1 and read 2 as separate FASTQ files.
    # -e is the number of threads. -O is the folder that receives the FASTQ files.
    fasterq-dump "${SRA_FILE}" \
        --split-files \
        -e "${THREADS}" \
        -t "${IBDAH_DIR}/tmp" \
        -O "${IBDAH_DIR}"
    IBDAH_READS=(
        "${IBDAH_DIR}/DRR513061_1.fastq"
        "${IBDAH_DIR}/DRR513061_2.fastq"
    )
fi

for fastq_path in "${IBDAH_READS[@]}"; do
    if [[ ! -f "${fastq_path}" ]]; then
        echo "ERROR: Expected Ibdah FASTQ was not created: ${fastq_path}"
        exit 1
    fi
done

echo "Ibdah FASTQ files:"
for fastq_path in "${IBDAH_READS[@]}"; do
    echo "  ${fastq_path}"
done

diamond makedb \
    --in "${IBDAH_QUERY}" \
    --db "${IBDAH_OUT}/queries" \
    --threads "${THREADS}"

for fastq_path in "${IBDAH_READS[@]}"; do
    sample_name="$(basename "${fastq_path}")"
    sample_name="${sample_name%.gz}"
    sample_name="${sample_name%.fastq}"
    sample_name="${sample_name%.fq}"
    hit_file="${IBDAH_OUT}/${sample_name}_hits.tsv"

    echo "------------------------------------------------------------"
    echo "Searching: ${fastq_path}"
    echo "Hits file: ${hit_file}"

    diamond blastx \
        --query "${fastq_path}" \
        --db "${IBDAH_OUT}/queries" \
        --out "${hit_file}" \
        --threads "${THREADS}" \
        --evalue 1e-5 \
        --max-target-seqs 1 \
        --sensitive \
        --tmpdir "${IBDAH_OUT}/tmp" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore

    if [[ ! -f "${hit_file}" ]]; then
        : > "${hit_file}"
    fi
    hit_count="$(wc -l < "${hit_file}" | tr -d ' ')"
    echo "Hits in this file: ${hit_count}"
done

# Same summary script as Part 1. --tissue leaf labels this whole run as
# leaf, because DRR513061 is not split into MFF and MFL sample names.
python3 "${SUMMARIZE_PY}" \
    --hits-dir "${IBDAH_OUT}" \
    --queries "${IBDAH_QUERY}" \
    --label "Ibdah DRR513061" \
    --tissue leaf

echo "Part 2 finished: $(date)"
echo "Your-read summary:  ${OUT_DIR}/SUMMARY_fragments.txt"
echo "Ibdah-read summary: ${IBDAH_OUT}/SUMMARY_fragments.txt"
