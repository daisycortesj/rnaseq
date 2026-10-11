#!/bin/bash
# =====================================================================
# SCRIPT:  search_mf30795_reads.sbatch
# PURPOSE: Search two read sets, one after the other, without copying
#          either set to a laptop.
#
#   PART 1 — your reads (fruit + leaf FASTQs that went into Trinity)
#            proteins: Mf41180, Mf31663, and SS
#            Mf30795 is skipped here. You already searched it.
#
#   PART 2 — the Ibdah public run (PRJDB17003, SRA run DRR513061)
#            proteins: Mf41180, Mf31663, Mf30795, and SS
#            This asks whether those proteins are in their dataset.
#
# WHY THIS STEP:
#   A missing Trinity hit only says the protein is absent from your
#   assembly. Part 1 asks if Mf41180 or Mf31663 is in your reads.
#   Part 2 asks if any of the three candidates is in the Ibdah reads.
#
# THE THREE PROTEINS (reference/queries/mf41180_mf31663_ss.faa):
#   Mf41180  — candidate not searched yet
#   Mf31663  — candidate not searched yet
#   SS       — santalene synthase, a positive control. An SS-like gene
#              (DN662) is already in this assembly, so some reads should
#              hit SS at about 40-55% identity. If SS gets zero hits,
#              the search setup failed. SS stays in the file even though
#              it was searched before, because it tells you the job worked.
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
# OUTPUTS:
#   Part 1: 06_analysis/mf41180_mf31663_read_search_<trinity|raw>/SUMMARY.txt
#   Part 2: 06_analysis/prjdb17003_DRR513061_search/SUMMARY.txt
#
# HOW TO READ SUMMARY.txt:
#   SS has hits around 40-55% identity
#       The search worked.
#   A candidate has many reads at 95-100% identity
#       That gene is in the reads. Trinity (or a later step) missed it.
#   A candidate has only a few reads at 95-100% identity
#       It is present, but too rare to assemble reliably.
#   A candidate's hits are all around 40-50% identity
#       Those reads are relatives. That protein itself is not in this data.
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

QUERY_FAA="${REPO_ROOT}/reference/queries/mf41180_mf31663_ss.faa"
OUT_DIR="${ANALYSIS_DIR}/mf41180_mf31663_read_search_${READ_SET}"

# Part 2: the Ibdah SRA run. Hits are saved under 06_analysis, not in rawdata.
IBDAH_QUERY="${REPO_ROOT}/reference/queries/ibdah_mf_candidates.faa"
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
    echo "Sync the repo to the cluster so reference/queries/mf41180_mf31663_ss.faa is there."
    exit 1
fi

if [[ ! -f "${IBDAH_QUERY}" ]]; then
    echo "ERROR: Ibdah query proteins not found: ${IBDAH_QUERY}"
    echo "Sync the repo so reference/queries/ibdah_mf_candidates.faa is there."
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
echo "DIAMOND blastx — are Mf41180 and Mf31663 in the nutmeg reads?"
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
# SUMMARIZE
# =====================================================================
# Python counts hits per protein and writes a plain-language conclusion.
# The job prints the same text to the log and saves SUMMARY.txt.

export SEARCH_OUT_DIR="${OUT_DIR}"
export SEARCH_READ_SET="${READ_SET}"
export SEARCH_QUERY_FAA="${QUERY_FAA}"

python3 << 'PY'
import glob
import os

out_dir = os.environ["SEARCH_OUT_DIR"]
read_set = os.environ["SEARCH_READ_SET"]
query_faa = os.environ["SEARCH_QUERY_FAA"]

# Read protein names from the FASTA headers, in file order.
# The first word after ">" is the name DIAMOND writes in the hit table.
protein_names = []
with open(query_faa) as handle:
    for line in handle:
        if line.startswith(">"):
            protein_names.append(line[1:].split()[0])

candidates = []
for name in protein_names:
    if name != "SS":
        candidates.append(name)

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

proteins = protein_names
lines = []
lines.append("Candidate read search: " + ", ".join(candidates))
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
        row = protein_name.ljust(12)
        row = row + tissue.ljust(9)
        row = row + str(high).rjust(16)
        row = row + str(mid).rjust(17)
        row = row + str(everything).rjust(11)
        row = row + "   " + format(max_identity, ".1f")
        lines.append(row)

ss_mid = count_where("SS", 35, 70, "all")
ss_all = count_where("SS", 0, 100, "all")

lines.append("")
lines.append("Conclusion")
if ss_all == 0:
    lines.append("SS has zero hits. The search did not work. Check that diamond ran and that the FASTQ paths are the nutmeg reads.")
else:
    for candidate in candidates:
        high = count_where(candidate, 95, 100, "all")
        everything = count_where(candidate, 0, 100, "all")
        max_identity = 0
        for sample, protein, identity in records:
            if protein == candidate and identity > max_identity:
                max_identity = identity
        if high >= 20:
            lines.append(candidate + " has many near-perfect reads (" + str(high) + " at >=95% identity). The gene is in these reads. Trinity, or a later step such as CD-HIT or TransDecoder, failed to keep it.")
        elif high >= 1:
            lines.append(candidate + " has only " + str(high) + " near-perfect read(s). It is present but rare, which often means expression was too low to assemble.")
        elif everything > 0 and max_identity < 70:
            lines.append(candidate + " hits top out at " + format(max_identity, ".1f") + "% identity. Those reads are relatives, not " + candidate + ". This protein is not in this read set.")
        elif everything == 0 and ss_mid > 0:
            lines.append("SS was found (the control worked) and " + candidate + " had no hits. " + candidate + " is not in this read set.")
        elif everything == 0:
            lines.append(candidate + " had no hits. SS also had no hits in the 35-70% window, so confirm the control before treating this as a true absence.")
        else:
            lines.append("Look at the table. " + candidate + " max identity is " + format(max_identity, ".1f") + "%. Near-perfect (>=95) count is " + str(high) + ".")

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

echo "Part 1 finished: $(date)"
echo "Your-read summary: ${OUT_DIR}/SUMMARY.txt"


# =====================================================================
# PART 2: SEARCH THE IBDAH PUBLIC RUN (PRJDB17003 / DRR513061)
# =====================================================================
# This is the other lab's nutmeg RNA-seq run. fasterq-dump turns the
# packed .sra file into normal FASTQ files. If those FASTQ files are
# already on disk, the dump is skipped.
#
# Proteins in this search: Mf41180, Mf31663, Mf30795, and SS.
# Mf30795 is included only here, because it was already searched in
# your own reads.

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

export SEARCH_OUT_DIR="${IBDAH_OUT}"
export SEARCH_QUERY_FAA="${IBDAH_QUERY}"

python3 << 'PY'
import glob
import os

out_dir = os.environ["SEARCH_OUT_DIR"]
query_faa = os.environ["SEARCH_QUERY_FAA"]

protein_names = []
with open(query_faa) as handle:
    for line in handle:
        if line.startswith(">"):
            protein_names.append(line[1:].split()[0])

candidates = []
for name in protein_names:
    if name != "SS":
        candidates.append(name)

records = []
hit_paths = sorted(glob.glob(os.path.join(out_dir, "*_hits.tsv")))
for path in hit_paths:
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line == "":
                continue
            columns = line.split("\t")
            # columns: read, protein, identity, length, evalue, bitscore
            protein = columns[1]
            identity = float(columns[2])
            length = float(columns[3])
            evalue = columns[4]
            bitscore = columns[5]
            read_name = columns[0]
            records.append((read_name, protein, identity, length, evalue, bitscore))

lines = []
lines.append("PRJDB17003 run DRR513061 (Ibdah) vs " + ", ".join(protein_names))
lines.append("Hit files: " + str(len(hit_paths)))
lines.append("Total hit rows: " + str(len(records)))
lines.append("")
lines.append("Counts are reads. R1 and R2 were both searched,")
lines.append("so one fragment can add 2.")
lines.append("")
lines.append("protein    all_hits   reads_95pct_and_30aa   max_identity")

for protein_name in protein_names:
    all_hits = 0
    strong_hits = 0
    max_identity = 0
    for read_name, protein, identity, length, evalue, bitscore in records:
        if protein != protein_name:
            continue
        all_hits = all_hits + 1
        if identity >= 95 and length >= 30:
            strong_hits = strong_hits + 1
        if identity > max_identity:
            max_identity = identity
    row = protein_name.ljust(11)
    row = row + str(all_hits).rjust(10)
    row = row + str(strong_hits).rjust(24)
    row = row + "   " + format(max_identity, ".1f")
    lines.append(row)

lines.append("")
lines.append("A strong hit is >=95% identity and at least 30 amino acids long.")
lines.append("Many strong hits: that protein is in the Ibdah reads.")
lines.append("Only ~40-50% identity: those reads are relatives, not that protein.")
lines.append("SS hits around 40-55%: the search worked.")
lines.append("")

for candidate in candidates:
    strong_rows = []
    for read_name, protein, identity, length, evalue, bitscore in records:
        if protein == candidate and length >= 30:
            strong_rows.append((identity, read_name, protein, length, evalue, bitscore))
    strong_rows.sort(reverse=True)
    lines.append("Top long " + candidate + " matches (up to 10):")
    if len(strong_rows) == 0:
        lines.append("  none with alignment length >= 30")
    else:
        shown = 0
        for identity, read_name, protein, length, evalue, bitscore in strong_rows:
            if shown == 10:
                break
            lines.append(
                "  "
                + read_name
                + "\t"
                + protein
                + "\t"
                + format(identity, ".1f")
                + "\t"
                + str(int(length))
                + "\t"
                + evalue
                + "\t"
                + bitscore
            )
            shown = shown + 1
    lines.append("")

text = "\n".join(lines) + "\n"
summary_path = os.path.join(out_dir, "SUMMARY.txt")
with open(summary_path, "w") as handle:
    handle.write(text)
print(text)
print("Wrote " + summary_path)
PY

echo "Part 2 finished: $(date)"
echo "Your-read summary:  ${OUT_DIR}/SUMMARY.txt"
echo "Ibdah-read summary: ${IBDAH_OUT}/SUMMARY.txt"
