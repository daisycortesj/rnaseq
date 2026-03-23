#!/bin/bash
# =====================================================================
# config.sh — Shared configuration for all pipeline scripts
# =====================================================================
#
# HOW TO USE:
#   All .sbatch scripts source this file automatically.
#   You don't need to run this file directly.
#
# FOR A NEW PROJECT:
#   Edit Section 1 (3 lines) and everything else updates automatically.
#
# =====================================================================


# ─────────────────────────────────────────────────────────────────────
# SECTION 1: EDIT THESE FOR YOUR PROJECT (only 3 lines to change)
# ─────────────────────────────────────────────────────────────────────

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"   # Root of your HPC data
SLURM_ACCOUNT="tholl_lab_1"                        # Your SLURM allocation
MAIL_USER="daisycortesj@vt.edu"                    # Email for job notifications


# ─────────────────────────────────────────────────────────────────────
# SECTION 2: DIRECTORIES (auto-derived from BASE_DIR — don't edit)
# ─────────────────────────────────────────────────────────────────────
#
# These match the HPC data layout:
#   00_rawdata/      → raw FASTQ files
#   01_processed/    → QC reports, fastp output, Trinity assemblies
#   02_mapped/       → STAR BAM files
#   03_count_tables/ → gene count matrices per species
#   04_reference/    → genome FASTA, GTF, STAR indices
#   05_rnaseq-code/  → this repository
#   06_analysis/     → DESeq2, BLAST, HMMER, plots
#   07_NRdatabase/   → BLAST databases, CYP450 database

RAWDATA_DIR="${BASE_DIR}/00_rawdata"
PROCESSED_DIR="${BASE_DIR}/01_processed"
MAPPED_DIR="${BASE_DIR}/02_mapped"
COUNT_DIR="${BASE_DIR}/03_count_tables"
REFERENCE_DIR="${BASE_DIR}/04_reference"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"
ANALYSIS_DIR="${BASE_DIR}/06_analysis"
DB_DIR="${BASE_DIR}/07_NRdatabase"


# ─────────────────────────────────────────────────────────────────────
# SECTION 3: CONDA ENVIRONMENT (change if your env has a different name)
# ─────────────────────────────────────────────────────────────────────

CONDA_ENV="rnaseq"

# ─────────────────────────────────────────────────────────────────────
# SECTION 4: SAMPLE GROUPS (add new samples here)
# ─────────────────────────────────────────────────────────────────────
#
# This is the ONE place to register sample groups.
# When you add new samples (e.g., 00_5_NewSpecies), add a new block
# to the function below and all scripts will recognize it.
#
# Each sample group needs:
#   CODE             — short code you type on the command line (DC, DG, etc.)
#   SPECIES_DIR      — folder name under 00_rawdata/ (e.g., 00_1_DC)
#   SPECIES_NAME     — full species name (for display in output)
#   GENOME_TYPE      — which genome to align to (carrot, nutmeg, etc.)
#   SAMPLE_CONDITIONS — (optional) how to extract condition from sample names
#
# SAMPLE_CONDITIONS tells build_count_matrix.py how to figure out which
# samples are Leaf vs Root (or any other condition). The format is:
#
#     SAMPLE_CONDITIONS="substring1=COND1 substring2=COND2"
#
# For each sample name, the script checks if it contains the substring,
# and assigns that condition. Examples:
#
#   Names like DC1L1, DC1R1      → leave SAMPLE_CONDITIONS="" (auto-detected)
#   Names like T_L_R1, T_R_R1   → SAMPLE_CONDITIONS="_L_=L _R_=R"
#   Names like CtrlA1, TreatA1  → SAMPLE_CONDITIONS="Ctrl=C Treat=T"

get_sample_info() {
    local code
    code=$(echo "$1" | tr '[:lower:]' '[:upper:]')

    # Reset so each call starts clean
    SAMPLE_CONDITIONS=""

    case "${code}" in
        DC)
            SPECIES_DIR="00_1_DC"
            SPECIES_NAME="D. carota ssp. maximus"
            GENOME_TYPE="carrot"
            ;;
        DG)
            SPECIES_DIR="00_2_DG"
            SPECIES_NAME="Daucus glaber"
            GENOME_TYPE="carrot"
            ;;
        MF)
            SPECIES_DIR="00_3_MF"
            SPECIES_NAME="Myristica fragrans"
            GENOME_TYPE="nutmeg"
            ;;
        SK)
            SPECIES_DIR="00_4_Sukman"
            SPECIES_NAME="Sukman samples"
            GENOME_TYPE="carrot"
            SAMPLE_CONDITIONS="_L_=L _R_=R"
            ;;
        DCDG)
            SPECIES_DIR="00_4_DC_DG"
            SPECIES_NAME="D. carota + D. glaber (combined)"
            GENOME_TYPE="carrot"
            ;;
        # ── ADD NEW SAMPLE GROUPS BELOW ──
        # XX)
        #     SPECIES_DIR="00_5_NewName"
        #     SPECIES_NAME="Species name"
        #     GENOME_TYPE="carrot"
        #     SAMPLE_CONDITIONS="_L_=L _R_=R"   # if sample names are non-standard
        #     ;;
        *)
            echo "ERROR: Unknown species code '${code}'"
            echo "Valid codes: DC, DG, MF, SK, DCDG"
            echo "To add a new code, edit SECTION 4 in scripts/config.sh"
            return 1
            ;;
    esac
}
