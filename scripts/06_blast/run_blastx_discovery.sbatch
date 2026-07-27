#!/bin/bash
# =============================================================================
# SCRIPT:  run_blastx_discovery.sbatch
#
# WHAT IT DOES (plain English):
#   Take EVERY transcript in your CD-HIT Trinity FASTA (DNA/RNA letters) and
#   search it against a big protein database (like NCBI nr) using BLASTX.
#   BLASTX translates your DNA in all reading frames and looks for similar
#   proteins. The result is an annotation table: "this transcript looks like
#   protein X from species Y".
#
# WHY BATCHES?
#   Searching hundreds of thousands of transcripts vs nr can take ~2 weeks.
#   ARC jobs often max out at 7 days. So we:
#     1. Split the big FASTA into small batches (e.g. 2000 sequences each)
#     2. Run BLASTX on one batch
#     3. SAVE progress (checkpoint)
#     4. Start the next batch
#   When time runs out, resubmit the SAME command. Finished batches are
#   SKIPPED. Only unfinished sequences are run. Nothing is repeated.
#
# -----------------------------------------------------------------------------
# INPUTS (what you need before running)
# -----------------------------------------------------------------------------
#   1. CD-HIT Trinity FASTA (nucleotide), for nutmeg MF usually:
#        01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta
#
#      Example lines:
#        >TRINITY_DN100_c0_g1_i1 len=1200
#        ATGCGTACGTACGTACGTACGT...
#        >TRINITY_DN101_c0_g1_i1 len=900
#        GCTAGCTAGCTAGCTAGCTA...
#
#   2. A BLAST protein database already downloaded on ARC, e.g. nr:
#        07_NRdatabase/blastdb/   (BLASTDB env points here)
#
#   3. conda env "rnaseq" with BLAST+ installed:
#        conda activate rnaseq
#        conda install -c bioconda blast
#
# -----------------------------------------------------------------------------
# OUTPUTS (files created in 06_analysis/blastx_MF_cdhit/ for MF + nr)
# -----------------------------------------------------------------------------
#   blastx_MF_nr_discovery.tsv
#       Main results table (tab-separated). One row = one BLAST hit.
#       Columns include: query ID, subject ID, title, % identity, e-value, etc.
#
#   blastx_MF_nr_discovery.done_ids.txt
#       Checkpoint list: every transcript ID already searched (hit OR no-hit).
#       Used so the next job does NOT repeat finished sequences.
#
#   blastx_MF_nr_discovery.progress.txt
#       Short status file you can open anytime to see how far you are.
#
#   blastx_discoveryfilter_<JOBID>.out / .err
#       SLURM log files (stdout / stderr) in 06_analysis/
#
# -----------------------------------------------------------------------------
# HOW TO RUN (on ARC)
# -----------------------------------------------------------------------------
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#
#   # First time OR resume after timeout (same command!):
#   sbatch scripts/06_blast/run_blastx_discovery.sbatch MF nr
#
#   # Check progress later:
#   cat 06_analysis/blastx_MF_cdhit/blastx_MF_nr_discovery.progress.txt
#
#   # Check if job is running:
#   squeue -u $USER
#
# WHERE THIS RUNS: HPC / ARC only (needs SLURM + large BLAST DB)
# =============================================================================


# -----------------------------------------------------------------------------
# SLURM RESOURCE REQUESTS
# Edit these if your allocation / needs change.
# -----------------------------------------------------------------------------
# --job-name     = name shown in squeue
# --account      = your lab allocation
# --partition    = queue name on ARC
# --cpus-per-task= CPU cores for BLAST threads
# --mem          = RAM
# --time         = max wall clock (7 days is ARC default max without long QoS)
# --mail-*       = email when job ends or fails
# --chdir        = working directory for log files
# -o / -e        = stdout and stderr log file names (%j = job ID)
#SBATCH --job-name=blastx_discoveryfilter
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/06_analysis
#SBATCH -o blastx_discoveryfilter_%j.out
#SBATCH -e blastx_discoveryfilter_%j.err


# Stop on undefined variables; keep going on some command failures so resume works
set -uo pipefail


# =============================================================================
# SECTION A — EASY SETTINGS TO EDIT
# Change these numbers/names when you want different behavior.
# You can also override BATCH_SIZE on the command line without editing:
#   BATCH_SIZE=1000 sbatch scripts/06_blast/run_blastx_discovery.sbatch MF nr
# =============================================================================

# Default species if you forget to type one after sbatch
#   MF = nutmeg (Trinity / CD-HIT)
#   DC / DG = carrot species (reference CDS path)
DEFAULT_SPECIES="MF"

# Default database if you forget to type one
#   nr         = huge NCBI non-redundant proteins (slow, most complete)
#   swissprot  = smaller curated set (faster, good first pass)
#   refseq_protein / DCdb = other options already on your system
DEFAULT_DATABASE="nr"

# How many transcripts to BLAST in ONE batch before saving a checkpoint.
#   Smaller number  = safer (lose less work if job dies mid-batch), slightly slower
#   Larger number   = fewer BLAST startups, but more work lost if a batch is killed
#   Beginner tip: leave at 2000 unless a single batch takes many hours
BATCH_SIZE="${BATCH_SIZE:-2000}"

# If less than this many minutes of job time remain, do NOT start a new batch.
#   Why: we only want to start a batch we can finish before SLURM kills the job.
MIN_REMAINING_MINUTES="${MIN_REMAINING_MINUTES:-90}"

# BLAST discovery filters (more permissive = more hits)
EVALUE="1e-4"            # max e-value (smaller = stricter)
QCOV_HSP_PERC="40"       # min % of query covered by the hit
MAX_TARGET_SEQS="25"     # max hits kept per query sequence


# =============================================================================
# SECTION B — ACTIVATE SOFTWARE + LOAD PATHS
# =============================================================================

# Turn on conda, then activate the env that has blastx
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# config.sh sets BASE_DIR, PROCESSED_DIR, ANALYSIS_DIR, DB_DIR, etc.
source "${SLURM_SUBMIT_DIR:-.}/scripts/config.sh"

# Tell BLAST where the downloaded databases live
export BLASTDB="${DB_DIR}/blastdb"


# =============================================================================
# SECTION C — READ WHAT YOU TYPED AFTER sbatch
# Example: sbatch ... run_blastx_discovery.sbatch MF nr
#   $1 = species code
#   $2 = database name
#   $3 = optional custom FASTA path
# =============================================================================

SPECIES="${1:-${DEFAULT_SPECIES}}"
DBNAME="${2:-${DEFAULT_DATABASE}}"
CUSTOM_QUERY="${3:-}"

# Check species is one we know
if [[ ! "$SPECIES" =~ ^(DC|DG|MF)$ ]]; then
    echo "ERROR: Invalid species '${SPECIES}'. Use DC, DG, or MF."
    echo "Usage:"
    echo "  sbatch scripts/06_blast/run_blastx_discovery.sbatch MF nr"
    exit 1
fi

# Check database name is one we know
if [[ ! "$DBNAME" =~ ^(nr|swissprot|refseq_protein|DCdb)$ ]]; then
    echo "ERROR: Invalid database '${DBNAME}'."
    echo "Use: nr, swissprot, refseq_protein, or DCdb"
    exit 1
fi


# =============================================================================
# SECTION D — CHOOSE THE INPUT FASTA (ALL transcripts / genes)
# Priority:
#   1) path you passed as argument 3
#   2) for MF: CD-HIT Trinity assembly
#   3) for DC/DG: extracted CDS FASTA from the reference genome
# =============================================================================

CDHIT_FASTA="${PROCESSED_DIR}/00_6_cdhit/${SPECIES}_trinity_cdhit95.fasta"
CDS_FASTA="${ANALYSIS_DIR}/blast_input_${SPECIES}/all_genes_cds.fasta"

if [ -n "${CUSTOM_QUERY}" ]; then
    # You gave an exact file path
    QUERY_FASTA="${CUSTOM_QUERY}"
    QUERY_BASENAME=$(basename "${QUERY_FASTA}")
    QUERY_BASENAME="${QUERY_BASENAME%.fasta}"
    QUERY_BASENAME="${QUERY_BASENAME%.fa}"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_${QUERY_BASENAME}"

elif [[ "${SPECIES}" == "MF" ]] && [ -f "${CDHIT_FASTA}" ]; then
    # Nutmeg default: ALL CD-HIT-collapsed Trinity transcripts
    QUERY_FASTA="${CDHIT_FASTA}"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_cdhit"

elif [ -f "${CDS_FASTA}" ]; then
    # Carrot / reference-genome CDS
    QUERY_FASTA="${CDS_FASTA}"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}"

elif [ -f "${CDHIT_FASTA}" ]; then
    # Fallback if a CD-HIT file exists for this species
    QUERY_FASTA="${CDHIT_FASTA}"
    OUTDIR="${ANALYSIS_DIR}/blastx_${SPECIES}_cdhit"

else
    echo "ERROR: Could not find an input FASTA for species ${SPECIES}."
    echo "Looked for:"
    echo "  ${CDHIT_FASTA}"
    echo "  ${CDS_FASTA}"
    echo ""
    echo "For nutmeg, run CD-HIT first:"
    echo "  sbatch scripts/04_busco/run_cdhit.sbatch MF"
    exit 1
fi

if [ ! -f "${QUERY_FASTA}" ]; then
    echo "ERROR: Query FASTA not found: ${QUERY_FASTA}"
    exit 1
fi


# =============================================================================
# SECTION E — PREPARE OUTPUT FOLDER + PRINT A JOB HEADER
# =============================================================================

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "=========================================="
echo "BLASTX  —  all transcripts  —  batch+resume"
echo "=========================================="
echo "Species code:         ${SPECIES}"
echo "Input FASTA:          ${QUERY_FASTA}"
echo "Protein database:     ${DBNAME}"
echo "Output folder:        ${OUTDIR}"
echo "Batch size:           ${BATCH_SIZE} sequences per batch"
echo "Safe-stop threshold:  ${MIN_REMAINING_MINUTES} minutes left"
echo "E-value cutoff:       ${EVALUE}"
echo "------------------------------------------"
echo "Started:   $(date)"
echo "Host:      $(hostname)"
echo "BLASTDB:   ${BLASTDB}"
echo "Job ID:    ${SLURM_JOB_ID:-not running under SLURM}"
echo "------------------------------------------"

# Quick check that blastx exists in this conda env
python3 -c "from rna_pipeline.tools.blast import check_blast; assert check_blast(), 'blastx not found — run: conda install -c bioconda blast'"
echo "=========================================="


# Pass Bash settings into Python through environment variables
export QUERY_FASTA
export DBNAME
export OUTDIR
export EVALUE
export QCOV_HSP_PERC
export MAX_TARGET_SEQS
export SPECIES
export BATCH_SIZE
export MIN_REMAINING_MINUTES

# If SLURM is about to kill the job (time limit), print a friendly message.
# Finished batches are already saved on disk.
trap 'echo ">>> SLURM sent SIGTERM (usually = time limit reached)."; \
      echo ">>> Finished batches are saved. Resubmit the SAME command to continue."; \
      exit 0' SIGTERM


# =============================================================================
# SECTION F — PYTHON: BATCH BLASTX + CHECKPOINT RESUME
# =============================================================================

python - <<'PY'
# ---------------------------------------------------------------------------
# This Python block does the real work:
#   1. Read the big FASTA
#   2. Skip IDs already in the "done" checkpoint file
#   3. Split remaining sequences into batches
#   4. For each batch: run blastx, save hits, mark IDs done
#   5. Stop early if not enough walltime left for another batch
# ---------------------------------------------------------------------------

import os
import subprocess
import sys
import time
from pathlib import Path


# ---- Read settings the Bash half exported ----
query_fasta = os.environ["QUERY_FASTA"]
database_name = os.environ["DBNAME"]
output_directory = os.environ["OUTDIR"]
evalue_cutoff = os.environ["EVALUE"]
query_coverage_percent = os.environ["QCOV_HSP_PERC"]
max_hits_per_query = os.environ["MAX_TARGET_SEQS"]
species_code = os.environ["SPECIES"]
number_of_cpu_threads = os.environ.get("SLURM_CPUS_PER_TASK", "1")
sequences_per_batch = int(os.environ.get("BATCH_SIZE", "2000"))
minimum_minutes_to_start_batch = int(os.environ.get("MIN_REMAINING_MINUTES", "90"))

# ---- Names of the three main output files ----
# (all live inside output_directory)
results_table_path = os.path.join(
    output_directory,
    "blastx_" + species_code + "_" + database_name + "_discovery.tsv",
)
done_ids_path = os.path.join(
    output_directory,
    "blastx_" + species_code + "_" + database_name + "_discovery.done_ids.txt",
)
progress_path = os.path.join(
    output_directory,
    "blastx_" + species_code + "_" + database_name + "_discovery.progress.txt",
)


# =============================================================================
# HELPER FUNCTIONS (small tools used below)
# =============================================================================

def parse_slurm_time_to_seconds(time_text):
    """
    SLURM prints remaining time in weird formats, for example:
      2-12:30:00   = 2 days + 12 hours + 30 minutes
      12:30:00     = 12 hours + 30 minutes
      45:00        = 45 minutes
    This function turns that text into a number of seconds.
    Returns None if we cannot understand the text.
    """
    time_text = time_text.strip()

    if time_text in ("", "NOT_SET", "UNLIMITED", "Invalid"):
        return None

    number_of_days = 0
    if "-" in time_text:
        day_part, time_text = time_text.split("-", 1)
        number_of_days = int(day_part)

    time_parts = time_text.split(":")

    if len(time_parts) == 3:
        hours = int(time_parts[0])
        minutes = int(time_parts[1])
        seconds = int(time_parts[2])
    elif len(time_parts) == 2:
        hours = 0
        minutes = int(time_parts[0])
        seconds = int(time_parts[1])
    elif len(time_parts) == 1:
        hours = 0
        minutes = 0
        seconds = int(time_parts[0])
    else:
        return None

    total_seconds = (
        number_of_days * 86400
        + hours * 3600
        + minutes * 60
        + seconds
    )
    return total_seconds


def get_remaining_walltime_seconds():
    """
    Ask SLURM: "How much time does THIS job still have?"
    Returns seconds as a number, or None if we are not under SLURM.
    """
    job_id = os.environ.get("SLURM_JOB_ID")
    if job_id is None:
        return None

    try:
        command = ["squeue", "-j", job_id, "-h", "-o", "%L"]
        result_text = subprocess.check_output(
            command,
            text=True,
            stderr=subprocess.DEVNULL,
        )
        result_text = result_text.strip()
        return parse_slurm_time_to_seconds(result_text)
    except Exception as error:
        print("WARNING: Could not read remaining walltime (" + str(error) + ")")
        return None


def load_done_sequence_ids(file_path):
    """
    Read the checkpoint file into a Python set of IDs.
    Example line in the file: TRINITY_DN100_c0_g1_i1
    """
    finished_ids = set()

    if not os.path.exists(file_path):
        return finished_ids

    checkpoint_file = open(file_path)
    for line in checkpoint_file:
        sequence_id = line.strip()
        if sequence_id != "":
            finished_ids.add(sequence_id)
    checkpoint_file.close()

    return finished_ids


def append_done_sequence_ids(file_path, sequence_id_list):
    """
    Add newly finished IDs to the checkpoint file.
    We APPEND (do not overwrite) so old progress is kept.
    """
    checkpoint_file = open(file_path, "a")
    for sequence_id in sequence_id_list:
        checkpoint_file.write(sequence_id + "\n")
    checkpoint_file.close()


def write_progress_file(file_path, total_count, done_count, batch_number, total_batches):
    """
    Write a tiny status file you can open with 'cat' anytime.
    """
    remaining_count = total_count - done_count
    progress_file = open(file_path, "w")
    progress_file.write("total_sequences\t" + str(total_count) + "\n")
    progress_file.write("done_sequences\t" + str(done_count) + "\n")
    progress_file.write("remaining_sequences\t" + str(remaining_count) + "\n")
    progress_file.write(
        "last_finished_batch\t" + str(batch_number) + "/" + str(total_batches) + "\n"
    )
    progress_file.write(
        "updated_at\t" + time.strftime("%Y-%m-%d %H:%M:%S") + "\n"
    )
    progress_file.close()


def read_fasta_one_record_at_a_time(fasta_path):
    """
    Walk through a FASTA file and yield one record at a time.

    For each record we return:
      sequence_id   = text after '>' up to the first space
      record_text   = the full FASTA entry (header + sequence lines)

    Why yield? So we do not need to hold the whole huge file in memory twice.
    """
    current_id = None
    current_lines = []

    fasta_file = open(fasta_path)
    for line in fasta_file:
        if line.startswith(">"):
            # We hit a new header — first save the previous record (if any)
            if current_id is not None:
                full_record = "".join(current_lines)
                yield current_id, full_record

            # New record starts here
            # Example header: >TRINITY_DN100_c0_g1_i1 len=1200
            # sequence_id becomes: TRINITY_DN100_c0_g1_i1
            header_without_arrow = line[1:]
            current_id = header_without_arrow.split()[0]
            current_lines = [line]
        else:
            # Sequence line (A/T/G/C...) belongs to the current record
            if current_id is not None:
                current_lines.append(line)

    # Do not forget the last record at the end of the file
    if current_id is not None:
        full_record = "".join(current_lines)
        yield current_id, full_record

    fasta_file.close()


def run_blastx_on_one_batch(batch_fasta_path, batch_results_path):
    """
    Call the blastx program on one small FASTA file.
    Writes a tab-separated hit table to batch_results_path.
    """
    blastx_command = [
        "blastx",
        "-query", batch_fasta_path,
        "-db", database_name,
        "-evalue", evalue_cutoff,
        "-max_target_seqs", max_hits_per_query,
        "-qcov_hsp_perc", query_coverage_percent,
        "-max_hsps", "1",
        "-num_threads", number_of_cpu_threads,
        # outfmt 6 = easy-to-parse table. Column names listed after '6'.
        "-outfmt",
        "6 qseqid sseqid stitle sscinames pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen qcovhsp",
        "-out", batch_results_path,
    ]

    print("Running BLASTX command:")
    print("  " + " ".join(blastx_command))
    subprocess.run(blastx_command, check=True)


# =============================================================================
# STEP 1 — LOAD CHECKPOINT + BUILD LIST OF SEQUENCES STILL NEEDED
# =============================================================================

print("")
print("STEP 1: Loading checkpoint and scanning the input FASTA...")

finished_sequence_ids = load_done_sequence_ids(done_ids_path)

# Backward compatible: if an older run only wrote the TSV (hits), treat those
# query IDs as done too. (No-hit sequences from old runs may be searched once.)
if os.path.exists(results_table_path):
    if os.path.getsize(results_table_path) > 0:
        results_file = open(results_table_path)
        for line in results_file:
            line = line.strip()
            if line != "":
                query_id = line.split("\t")[0]
                finished_sequence_ids.add(query_id)
        results_file.close()

total_sequence_count = 0
pending_records = []  # list of [sequence_id, fasta_text] still needing BLAST

for sequence_id, record_text in read_fasta_one_record_at_a_time(query_fasta):
    total_sequence_count = total_sequence_count + 1

    # Only keep sequences we have NOT finished yet
    if sequence_id not in finished_sequence_ids:
        pending_records.append([sequence_id, record_text])

already_done_count = total_sequence_count - len(pending_records)

print("  Total sequences in FASTA:     " + str(total_sequence_count))
print("  Already finished (skip):      " + str(already_done_count))
print("  Still need BLAST:             " + str(len(pending_records)))

if len(pending_records) == 0:
    print("")
    print("All sequences are already done. Nothing left to BLAST.")
    write_progress_file(
        progress_path,
        total_sequence_count,
        total_sequence_count,
        0,
        0,
    )
    sys.exit(0)


# =============================================================================
# STEP 2 — SPLIT REMAINING SEQUENCES INTO BATCHES
# =============================================================================

print("")
print("STEP 2: Splitting remaining sequences into batches...")

batch_list = []
start_index = 0
while start_index < len(pending_records):
    end_index = start_index + sequences_per_batch
    one_batch = pending_records[start_index:end_index]
    batch_list.append(one_batch)
    start_index = end_index

number_of_batches = len(batch_list)

print("  Sequences per batch:          " + str(sequences_per_batch))
print("  Batches left in this job:     " + str(number_of_batches))
print(
    "  Will safe-stop if walltime < "
    + str(minimum_minutes_to_start_batch)
    + " minutes"
)


# =============================================================================
# STEP 3 — RUN EACH BATCH (save after every one)
# =============================================================================

print("")
print("STEP 3: Running BLASTX batch by batch...")

# Make sure the main results file exists so we can always append to it
Path(results_table_path).touch(exist_ok=True)

batch_number = 0
for one_batch in batch_list:
    batch_number = batch_number + 1

    # ---- Safe-stop check: enough time left to start another batch? ----
    remaining_seconds = get_remaining_walltime_seconds()
    if remaining_seconds is not None:
        remaining_minutes = remaining_seconds / 60.0
        print("")
        print(
            "  Walltime left before batch "
            + str(batch_number)
            + "/"
            + str(number_of_batches)
            + ": "
            + str(round(remaining_minutes, 1))
            + " minutes"
        )

        if remaining_minutes < minimum_minutes_to_start_batch:
            done_so_far = min(len(finished_sequence_ids), total_sequence_count)
            print("")
            print(">>> SAFE STOP")
            print(
                ">>> Only "
                + str(round(remaining_minutes, 1))
                + " minutes left."
            )
            print(
                ">>> Need at least "
                + str(minimum_minutes_to_start_batch)
                + " minutes to start a new batch."
            )
            print(
                ">>> Progress saved: "
                + str(done_so_far)
                + "/"
                + str(total_sequence_count)
                + " sequences"
            )
            print(">>> Resubmit the SAME command to continue:")
            print(
                ">>>   sbatch scripts/06_blast/run_blastx_discovery.sbatch "
                + species_code
                + " "
                + database_name
            )
            write_progress_file(
                progress_path,
                total_sequence_count,
                done_so_far,
                batch_number - 1,
                number_of_batches,
            )
            sys.exit(0)

    # ---- Collect IDs and write a temporary batch FASTA ----
    batch_sequence_ids = []
    for record in one_batch:
        batch_sequence_ids.append(record[0])

    batch_fasta_path = os.path.join(
        output_directory,
        "_batch_"
        + species_code
        + "_"
        + database_name
        + "_"
        + str(batch_number).zfill(5)
        + ".fasta",
    )
    batch_tsv_path = os.path.join(
        output_directory,
        "_batch_"
        + species_code
        + "_"
        + database_name
        + "_"
        + str(batch_number).zfill(5)
        + ".tsv",
    )

    print("")
    print(
        "========== BATCH "
        + str(batch_number)
        + "/"
        + str(number_of_batches)
        + " =========="
    )
    print("  Sequences in this batch: " + str(len(batch_sequence_ids)))
    print("  First ID: " + batch_sequence_ids[0])
    print("  Last ID:  " + batch_sequence_ids[-1])

    batch_fasta_file = open(batch_fasta_path, "w")
    for record in one_batch:
        record_text = record[1]
        batch_fasta_file.write(record_text)
    batch_fasta_file.close()

    # ---- Run BLASTX on this batch ----
    batch_start_time = time.time()
    try:
        run_blastx_on_one_batch(batch_fasta_path, batch_tsv_path)
    except subprocess.CalledProcessError as error:
        print(
            "ERROR: BLASTX failed on batch "
            + str(batch_number)
            + " with exit code "
            + str(error.returncode),
            file=sys.stderr,
        )
        print(
            "These IDs were NOT marked done, so they will be retried next time.",
            file=sys.stderr,
        )
        if os.path.exists(batch_fasta_path):
            os.remove(batch_fasta_path)
        if os.path.exists(batch_tsv_path):
            os.remove(batch_tsv_path)
        sys.exit(1)

    batch_elapsed_seconds = time.time() - batch_start_time
    batch_elapsed_minutes = batch_elapsed_seconds / 60.0

    # ---- Append hit rows into the MAIN results table (permanent save) ----
    hit_row_count = 0
    if os.path.exists(batch_tsv_path):
        if os.path.getsize(batch_tsv_path) > 0:
            batch_results_file = open(batch_tsv_path)
            main_results_file = open(results_table_path, "a")
            for line in batch_results_file:
                main_results_file.write(line)
                hit_row_count = hit_row_count + 1
            batch_results_file.close()
            main_results_file.close()

    # ---- Mark EVERY ID in this batch as done (hit OR no-hit) ----
    # This is why we do not repeat sequences on the next job.
    append_done_sequence_ids(done_ids_path, batch_sequence_ids)
    for sequence_id in batch_sequence_ids:
        finished_sequence_ids.add(sequence_id)

    # ---- Delete temporary batch files to save disk space ----
    if os.path.exists(batch_fasta_path):
        os.remove(batch_fasta_path)
    if os.path.exists(batch_tsv_path):
        os.remove(batch_tsv_path)

    done_display = min(len(finished_sequence_ids), total_sequence_count)
    write_progress_file(
        progress_path,
        total_sequence_count,
        done_display,
        batch_number,
        number_of_batches,
    )

    print("  Finished in " + str(round(batch_elapsed_minutes, 1)) + " minutes")
    print("  Hit rows saved:     " + str(hit_row_count))
    print("  IDs checkpointed:   " + str(len(batch_sequence_ids)))
    print(
        "  Overall progress:   "
        + str(done_display)
        + "/"
        + str(total_sequence_count)
    )


# =============================================================================
# STEP 4 — FINAL SUMMARY
# =============================================================================

print("")
print("========================================")
print("BLASTX finished for all remaining batches in THIS job")
print("========================================")
print("Results table:   " + results_table_path)
print("Done-ID file:    " + done_ids_path)
print("Progress file:   " + progress_path)
print("")

if os.path.exists(results_table_path):
    total_lines = 0
    unique_query_ids = set()
    results_file = open(results_table_path)
    for line in results_file:
        line = line.strip()
        if line != "":
            total_lines = total_lines + 1
            query_id = line.split("\t")[0]
            unique_query_ids.add(query_id)
    results_file.close()

    print("Total hit rows in table:              " + str(total_lines))
    print("Unique queries that had at least 1 hit: " + str(len(unique_query_ids)))
    print("")
    print("First 5 lines of the results table:")
    print("------------------------------------------")
    results_file = open(results_table_path)
    line_number = 0
    for line in results_file:
        print(line.rstrip())
        line_number = line_number + 1
        if line_number >= 5:
            break
    results_file.close()
    print("------------------------------------------")

print("")
print("If progress is not 100% yet, resubmit the SAME sbatch command.")
print("Finished sequences will be skipped automatically.")
PY


# =============================================================================
# SECTION G — JOB FOOTER (always runs after Python exits)
# =============================================================================

echo "=========================================="
echo "BLASTX SLURM job finished"
echo "Time: $(date)"
echo ""
echo "Check progress anytime with:"
echo "  cat ${OUTDIR}/blastx_${SPECIES}_${DBNAME}_discovery.progress.txt"
echo ""
echo "Main results table:"
echo "  ${OUTDIR}/blastx_${SPECIES}_${DBNAME}_discovery.tsv"
echo ""
echo "To resume after a 7-day timeout (same command):"
echo "  sbatch scripts/06_blast/run_blastx_discovery.sbatch ${SPECIES} ${DBNAME}"
echo "=========================================="
