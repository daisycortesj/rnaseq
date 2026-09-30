#!/bin/bash
# =============================================================================
# SCRIPT:  run_minimap2_fungi.sbatch
# PURPOSE: Align nutmeg transcripts to your fungal genomes, then keep the
#          transcripts that did not match.
#
# Professor's command, with your real file names:
#   minimap2 -ax splice combined_fungi.fasta Trinity.fasta > transcripts_vs_fungi.sam
#
#   combined_fungi.fasta  = fungal_contaminants.fna   (the 12 GB file)
#   Trinity.fasta         = MF_trinity_cdhit95.fasta
#   transcripts_vs_fungi.sam = the alignment file (already created)
#
# Mapped transcripts match a fungus. Unmapped transcripts did not.
# Those unmapped transcripts are saved as MF_unmapped_transcripts.fasta.
#
# The SAM file already exists, so this script SKIPS the minimap2 command.
# It starts at the "keep unmapped" step.
#
# HOW TO RUN (on ARC, from the rnaseq repo):
#   cd /projects/tholl_lab_1/daisy_analysis/05_rnaseq-code
#   sbatch scripts/01_qc/run_minimap2_fungi.sbatch
#
# TO CHANGE INPUTS OR THE OUTPUT NAME:
#   Edit the "CHANGE THESE" section below. Then submit the job again.
# =============================================================================


# =============================================================================
# SLURM — how much of the computer this job is allowed to use
#
# --cpus-per-task   number of CPUs. If you change this, also change THREADS
#                   in the CHANGE THESE section so the two numbers match.
# --mem             memory. Indexing the 12 GB fungal FASTA was killed at 64G
#                   (out of memory). 240G matches the Trinity job on this queue.
# --time            longest the job may run (hours:minutes:seconds)
# --chdir           folder the job starts in
# -o / -e           log files. %j is replaced with the job number.
# =============================================================================
#SBATCH --job-name=minimap2_fungi
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 00_rawdata/00_5_Fungus/minimap2_fungi_%j.out
#SBATCH -e 00_rawdata/00_5_Fungus/minimap2_fungi_%j.err

# Stop the script if a command fails, if a variable is empty, or if a
# command in a pipeline fails.
set -euo pipefail


# =============================================================================
# CHANGE THESE — edit the paths here if your files move or you want a new name
# =============================================================================

# combined_fungi.fasta
# This is the fungal genomes. Do NOT point this at fungal_contaminants.mmi.
# That .mmi file is 0 bytes (empty), so minimap2 would align to nothing.
FUNGI_FASTA="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/fungal_contaminants.fna"

# Trinity.fasta
# CD-HIT transcripts (the same FASTA used for RSEM).
TRANSCRIPTS="/projects/tholl_lab_1/daisy_analysis/01_processed/00_6_cdhit/MF_trinity_cdhit95.fasta"

# transcripts_vs_fungi.sam
# Already written by the minimap2 job. This script reads it and does not
# overwrite it.
OUTPUT_SAM="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/transcripts_vs_fungi.sam"

# How many CPUs minimap2 should use. Keep this equal to --cpus-per-task above.
# Unused while the minimap2 command stays skipped.
THREADS=16

# Names that matched a fungus.
MAPPED_IDS="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/MF_mapped_ids.txt"

# Every transcript name in the CD-HIT file.
ALL_IDS="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/MF_all_ids.txt"

# Names that did not match a fungus.
UNMAPPED_IDS="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/MF_unmapped_ids.txt"

# FASTA of those unmapped transcripts. This is the file to keep.
UNMAPPED_FASTA="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/MF_unmapped_transcripts.fasta"


# =============================================================================
# START THE ENVIRONMENT
# =============================================================================
echo "Activating conda environment rnaseq..."
source ~/.bashrc
conda activate rnaseq

# command -v checks that the program is installed and prints its path
if ! command -v minimap2 > /dev/null 2>&1; then
    echo "ERROR: minimap2 was not found."
    echo "Install it once on a login node:"
    echo "  conda activate rnaseq"
    echo "  conda install -c bioconda minimap2"
    exit 1
fi

echo "minimap2: $(command -v minimap2)"
echo "Started:  $(date)"
echo ""


# =============================================================================
# CHECK THE INPUT FILES BEFORE ALIGNING
#
# -s means "this path exists and the file is not empty."
# An empty file would make minimap2 finish in a few seconds with 0 sequences.
# =============================================================================
if [[ ! -s "${FUNGI_FASTA}" ]]; then
    echo "ERROR: Fungal FASTA is missing or empty:"
    echo "  ${FUNGI_FASTA}"
    exit 1
fi

if [[ ! -s "${TRANSCRIPTS}" ]]; then
    echo "ERROR: Transcript FASTA is missing or empty:"
    echo "  ${TRANSCRIPTS}"
    exit 1
fi

if [[ ! -s "${OUTPUT_SAM}" ]]; then
    echo "ERROR: SAM file is missing or empty:"
    echo "  ${OUTPUT_SAM}"
    echo "The minimap2 step is skipped, so this file must already exist."
    exit 1
fi

echo "Fungal genomes:  ${FUNGI_FASTA}"
echo "                 $(ls -lh "${FUNGI_FASTA}" | awk '{print $5}')"
echo "Transcripts:     ${TRANSCRIPTS}"
echo "Existing SAM:    ${OUTPUT_SAM}"
echo "                 $(ls -lh "${OUTPUT_SAM}" | awk '{print $5}')"
echo ""


# =============================================================================
# ALIGN — SKIPPED
#
# The SAM file was already created. Leave the command here so it can be
# turned back on later. Delete the # symbols on the two minimap2 lines
# if you need to build transcripts_vs_fungi.sam again.
#
# This is the professor's command with your files filled in:
#   minimap2 -ax splice combined_fungi.fasta Trinity.fasta > transcripts_vs_fungi.sam
#
# Flags:
#   -a   write a SAM alignment file (text)
#   -x splice   transcripts are spliced mRNA, the fungal file is genomic DNA
#   -t   number of CPUs
#
# The first FASTA is the reference (fungi).
# The second FASTA is the queries (your transcripts).
# > sends the alignment into OUTPUT_SAM instead of onto the screen.
# =============================================================================
echo "Skipping minimap2. The SAM file is already made."
echo "  ${OUTPUT_SAM}"

# minimap2 -ax splice -t "${THREADS}" "${FUNGI_FASTA}" "${TRANSCRIPTS}" > "${OUTPUT_SAM}"


# =============================================================================
# STEP 4 — list transcript names
#
# A SAM flag stores yes/no facts about an alignment.
# Dividing the flag by 4 and looking at the remainder tells us one fact:
#   remainder 0  -> this row MATCHED a fungus
#   remainder 1  -> this row did not match
#
# sort -u keeps each name once. The SAM lists some names twice.
# =============================================================================
echo ""
echo "Finding transcript names that matched a fungus..."

grep -v '^@' "${OUTPUT_SAM}" \
  | awk '{ if (int($2/4) % 2 == 0) print $1 }' \
  | sort -u > "${MAPPED_IDS}"

echo "Finding every transcript name in the CD-HIT FASTA..."

grep '^>' "${TRANSCRIPTS}" \
  | sed 's/^>//' \
  | awk '{print $1}' \
  | sort -u > "${ALL_IDS}"

# comm -23 keeps names that are in the first file and absent from the second.
# Those names never matched a fungus.
echo "Keeping names that did not match..."

comm -23 "${ALL_IDS}" "${MAPPED_IDS}" > "${UNMAPPED_IDS}"

echo "Matched a fungus:  $(wc -l < "${MAPPED_IDS}" | tr -d ' ')"
echo "Did not match:     $(wc -l < "${UNMAPPED_IDS}" | tr -d ' ')"
echo ""


# =============================================================================
# STEP 5 — copy the unmapped sequences into MF_unmapped_transcripts.fasta
#
# Python reads the CD-HIT file one line at a time.
# When a header name is on the unmapped list, it writes that transcript.
# The original CD-HIT file is left as it is.
# =============================================================================
echo "Writing ${UNMAPPED_FASTA}"

python3 - "${TRANSCRIPTS}" "${UNMAPPED_IDS}" "${UNMAPPED_FASTA}" << 'PY'
import sys

transcript_fasta = sys.argv[1]
unmapped_id_file = sys.argv[2]
output_fasta = sys.argv[3]

# A set stores each name once, so looking up a name is simple.
unmapped_ids = set()
id_handle = open(unmapped_id_file, "r")
for line in id_handle:
    name = line.strip()
    if name != "":
        unmapped_ids.add(name)
id_handle.close()

keep_this_transcript = False
written = 0

fasta_handle = open(transcript_fasta, "r")
out_handle = open(output_fasta, "w")

for line in fasta_handle:
    # A header line starts with > and holds the transcript name.
    if line.startswith(">"):
        header = line[1:].strip()
        name = header.split()[0]
        if name in unmapped_ids:
            keep_this_transcript = True
            written = written + 1
            out_handle.write(line)
        else:
            keep_this_transcript = False
    else:
        # Sequence lines belong to the header above them.
        if keep_this_transcript:
            out_handle.write(line)

fasta_handle.close()
out_handle.close()
print("Wrote " + str(written) + " transcripts")
PY

if [[ ! -s "${UNMAPPED_FASTA}" ]]; then
    echo "ERROR: Unmapped FASTA was not written:"
    echo "  ${UNMAPPED_FASTA}"
    exit 1
fi

echo ""
echo "Finished: $(date)"
ls -lh "${OUTPUT_SAM}" "${UNMAPPED_IDS}" "${UNMAPPED_FASTA}"
echo "Plant transcripts saved as: ${UNMAPPED_FASTA}"
