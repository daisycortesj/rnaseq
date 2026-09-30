#!/bin/bash
# =============================================================================
# SCRIPT:  run_minimap2_fungi.sbatch
# PURPOSE: Align nutmeg transcripts to your fungal genomes.
#
# Professor's command, with your real file names:
#   minimap2 -ax splice combined_fungi.fasta Trinity.fasta > transcripts_vs_fungi.sam
#
#   combined_fungi.fasta  = fungal_contaminants.fna   (the 12 GB file)
#   Trinity.fasta         = MF_trinity_cdhit95.fasta
#   transcripts_vs_fungi.sam = the alignment file this job writes
#
# Mapped transcripts match a fungus. Unmapped transcripts did not.
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
# The alignment is written here.
OUTPUT_SAM="/projects/tholl_lab_1/daisy_analysis/00_rawdata/00_5_Fungus/transcripts_vs_fungi.sam"

# How many CPUs minimap2 should use. Keep this equal to --cpus-per-task above.
THREADS=16


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

echo "Fungal genomes:  ${FUNGI_FASTA}"
echo "                 $(ls -lh "${FUNGI_FASTA}" | awk '{print $5}')"
echo "Transcripts:     ${TRANSCRIPTS}"
echo "Output SAM:      ${OUTPUT_SAM}"
echo ""


# =============================================================================
# ALIGN
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
echo "Running minimap2..."
echo "  minimap2 -ax splice -t ${THREADS} ${FUNGI_FASTA} ${TRANSCRIPTS} > ${OUTPUT_SAM}"

minimap2 -ax splice -t "${THREADS}" "${FUNGI_FASTA}" "${TRANSCRIPTS}" > "${OUTPUT_SAM}"

if [[ ! -s "${OUTPUT_SAM}" ]]; then
    echo "ERROR: The SAM file is missing or empty:"
    echo "  ${OUTPUT_SAM}"
    exit 1
fi

echo ""
echo "Finished: $(date)"
ls -lh "${OUTPUT_SAM}"
echo "Alignment file is ready: ${OUTPUT_SAM}"
