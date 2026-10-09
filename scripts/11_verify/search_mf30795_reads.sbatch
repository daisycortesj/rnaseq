#!/bin/bash
# =====================================================================
# SCRIPT:  search_mf30795_reads.sbatch
# PURPOSE: Ask whether the Mf30795 protein is present in the nutmeg
#          RNA-seq reads, without copying those reads to a laptop.
#
# WHY THIS STEP:
#   tblastn against the pooled Trinity assembly did not find Mf30795.
#   That only says the gene is missing from the assembled transcripts.
#   This job looks one step earlier: the FASTQ reads Trinity was given.
#   DIAMOND blastx translates each read and compares it to two proteins.
#
# THE TWO PROTEINS (reference/queries/mf30795_ss.faa):
#   Mf30795  — the protein you are looking for
#   SS       — santalene synthase, a positive control. An SS-like gene
#              (DN662) is already in this assembly, so some reads should
#              hit SS at about 40-55% identity. If SS gets zero hits,
#              the search setup failed.
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
#   conda install -c bioconda diamond          # only the first time
#   sbatch scripts/11_verify/search_mf30795_reads.sbatch
#
# OUTPUTS (under 06_analysis/mf30795_read_search_<trinity|raw>/):
#   *_hits.tsv     one file per FASTQ
#   SUMMARY.txt    how many reads hit each protein, plus a conclusion
#
# HOW TO READ SUMMARY.txt:
#   SS has hits around 40-55% identity
#       The search worked.
#   Mf30795 has many reads at 95-100% identity
#       The gene is in the reads. Trinity (or a later step) missed it.
#   Mf30795 has only a few reads at 95-100% identity
#       It is present, but too rare to assemble reliably.
#   Mf30795 hits are all around 40-50% identity
#       Those reads are relatives. Mf30795 itself is not in this data.
# =====================================================================


# ===== SLURM directives =====
# Do not run diamond on the login node. sbatch sends this to a compute node.
#SBATCH --job-name=mf30795_reads
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/mf30795_reads_%j.out
#SBATCH -e 06_analysis/mf30795_reads_%j.err

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

QUERY_FAA="${REPO_ROOT}/reference/queries/mf30795_ss.faa"
OUT_DIR="${ANALYSIS_DIR}/mf30795_read_search_${READ_SET}"
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
    echo "Sync the repo to the cluster so reference/queries/mf30795_ss.faa is there."
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

if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    source ~/.bashrc
    conda activate rnaseq
fi

if ! command -v diamond &> /dev/null; then
    echo "ERROR: diamond is not installed in the rnaseq environment."
    echo "On a login node (not in this job), run:"
    echo "  conda activate rnaseq"
    echo "  conda install -c bioconda diamond"
    echo "Then resubmit this script."
    exit 1
fi

echo "============================================================"
echo "DIAMOND blastx — is Mf30795 in the nutmeg reads?"
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
# --max-target-seqs 1    keep only the better of Mf30795 vs SS
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
# SUMMARIZE
# =====================================================================
# Python counts hits per protein and writes a plain-language conclusion.
# The job prints the same text to the log and saves SUMMARY.txt.

export SEARCH_OUT_DIR="${OUT_DIR}"
export SEARCH_READ_SET="${READ_SET}"

python3 << 'PY'
import glob
import os

out_dir = os.environ["SEARCH_OUT_DIR"]
read_set = os.environ["SEARCH_READ_SET"]

records = []
hit_paths = sorted(glob.glob(os.path.join(out_dir, "*_hits.tsv")))

for path in hit_paths:
    sample = os.path.basename(path).replace("_hits.tsv", "")
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line == "":
                continue
            columns = line.split("\t")
            protein = columns[1]
            identity = float(columns[2])
            records.append((sample, protein, identity))

def count_where(protein_name, minimum, maximum, tissue):
    """Count hits for one protein, identity window, and tissue."""
    total = 0
    for sample, protein, identity in records:
        if protein != protein_name:
            continue
        if identity < minimum or identity > maximum:
            continue
        if tissue == "fruit" and "MFF" not in sample:
            continue
        if tissue == "leaf" and "MFL" not in sample:
            continue
        total = total + 1
    return total

proteins = ["Mf30795", "SS"]
lines = []
lines.append("Mf30795 read search")
lines.append("Read set: " + read_set)
lines.append("Hit files: " + str(len(hit_paths)))
lines.append("Total hit rows: " + str(len(records)))
lines.append("")
lines.append("Counts are reads, not pairs. R1 and R2 were both searched,")
lines.append("so one fragment can add 2 to the count.")
lines.append("")
lines.append("protein   tissue   reads_95_to_100   reads_35_to_70   all_hits   max_identity")

for protein_name in proteins:
    for tissue in ["fruit", "leaf", "all"]:
        high = count_where(protein_name, 95, 100, tissue)
        mid = count_where(protein_name, 35, 70, tissue)
        everything = count_where(protein_name, 0, 100, tissue)
        max_identity = 0
        for sample, protein, identity in records:
            if protein != protein_name:
                continue
            if tissue == "fruit" and "MFF" not in sample:
                continue
            if tissue == "leaf" and "MFL" not in sample:
                continue
            if identity > max_identity:
                max_identity = identity
        row = protein_name.ljust(10)
        row = row + tissue.ljust(9)
        row = row + str(high).rjust(16)
        row = row + str(mid).rjust(17)
        row = row + str(everything).rjust(11)
        row = row + "   " + format(max_identity, ".1f")
        lines.append(row)

mf_high = count_where("Mf30795", 95, 100, "all")
ss_mid = count_where("SS", 35, 70, "all")
ss_all = count_where("SS", 0, 100, "all")
mf_max = 0
mf_all = count_where("Mf30795", 0, 100, "all")
for sample, protein, identity in records:
    if protein == "Mf30795" and identity > mf_max:
        mf_max = identity

lines.append("")
lines.append("Conclusion")
if ss_all == 0:
    lines.append("SS has zero hits. The search did not work. Check that diamond ran and that the FASTQ paths are the nutmeg reads.")
elif mf_high >= 20:
    lines.append("Mf30795 has many near-perfect reads (" + str(mf_high) + " at >=95% identity). The gene is in these reads. Trinity, or a later step such as CD-HIT or TransDecoder, failed to keep it.")
elif mf_high >= 1:
    lines.append("Mf30795 has only " + str(mf_high) + " near-perfect read(s). It is present but rare, which often means expression was too low to assemble.")
elif mf_all > 0 and mf_max < 70:
    lines.append("Mf30795 hits top out at " + format(mf_max, ".1f") + "% identity. Those reads are relatives, not Mf30795. This protein is not in this read set.")
elif mf_all == 0 and ss_mid > 0:
    lines.append("SS was found (the control worked) and Mf30795 had no hits. Mf30795 is not in this read set.")
elif mf_all == 0:
    lines.append("Mf30795 had no hits. SS also had no hits in the 35-70% window, so confirm the control before treating this as a true absence.")
else:
    lines.append("Look at the table. Mf30795 max identity is " + format(mf_max, ".1f") + "%. Near-perfect (>=95) count is " + str(mf_high) + ".")

if ss_mid > 0:
    lines.append("Positive control: SS has " + str(ss_mid) + " hits between 35% and 70% identity, which is the expected range for the SS-like gene already in this project.")
elif ss_all > 0:
    lines.append("SS has hits, but none fell in the 35-70% window. Open the per-file tables before deciding the control failed.")

text = "\n".join(lines) + "\n"
summary_path = os.path.join(out_dir, "SUMMARY.txt")
with open(summary_path, "w") as handle:
    handle.write(text)
print(text)
print("Wrote " + summary_path)
PY

echo "Finished: $(date)"
echo "Read this file: ${OUT_DIR}/SUMMARY.txt"
