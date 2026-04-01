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
# SECTION 4: SAMPLE NAMING CONVENTION + SAMPLE GROUPS
# ─────────────────────────────────────────────────────────────────────
#
# ▸ NAMING CONVENTION (follow this when adding new FASTQ files)
#
#   Name your paired-end FASTQ files like this:
#
#     {VARIETY}_{CONDITION}_{REPLICATE}_{READ}.fq.gz
#
#     VARIETY   = species code + variety number  (DC1, DC2, DG, MF, …)
#     CONDITION = single letter for tissue type  (L=leaf, R=root, F=fruit, …)
#     REPLICATE = replicate number               (1, 2, 3)
#     READ      = paired-end read number         (1 or 2)
#
#   Examples:
#     DC1_L_1_1.fq.gz   DC1_L_1_2.fq.gz    (DC variety 1, Leaf, replicate 1)
#     DC2_R_3_1.fq.gz   DC2_R_3_2.fq.gz    (DC variety 2, Root, replicate 3)
#     DG_L_1_1.fq.gz    DG_L_1_2.fq.gz     (DG, Leaf, replicate 1)
#     MF_F_2_1.fq.gz    MF_F_2_2.fq.gz     (MF, Fruit, replicate 2)
#
#   If you only have ONE variety, drop the number:
#     DG_L_1_1.fq.gz  (not DG1_L_1_1.fq.gz)
#
#   If you have MULTIPLE varieties, include the number:
#     DC1_L_1_1.fq.gz  DC2_L_1_1.fq.gz
#
#   The pipeline auto-detects variety, condition, and replicate from
#   these names. No SAMPLE_CONDITIONS needed if you follow this format.
#
# ▸ LEGACY NAMES (already processed — handled by SAMPLE_CONDITIONS)
#
#   Older samples may not follow the convention above. The pipeline
#   handles them via regex or SAMPLE_CONDITIONS:
#
#     DC1L1_1.fq.gz     → auto-detected (variety=DC1, condition=L, rep=1)
#     DGL1_1.fq.gz      → auto-detected (variety=DG, condition=L, rep=1)
#     MFF1_1.fq.gz      → needs SAMPLE_CONDITIONS="MFF=F MFL=L"
#     T_L_R1_1.fq.gz    → needs SAMPLE_CONDITIONS="_L_=L _R_=R"
#
# ▸ SAMPLE GROUPS
#
#   This is the ONE place to register sample groups.
#   When you add new samples, add a new block to the function below.
#
#   Each sample group needs:
#     CODE             — short code you type on the command line (DC, DG, etc.)
#     SPECIES_DIR      — folder name under 00_rawdata/ (e.g., 00_1_DC)
#     SPECIES_NAME     — full species name (for display in output)
#     GENOME_TYPE      — which genome to align to (carrot, nutmeg, etc.)
#     SAMPLE_CONDITIONS — (optional) only needed for legacy/non-standard names

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
            SAMPLE_CONDITIONS="MFF=F MFL=L"
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
        #
        # If your FASTQ files follow the naming convention above
        # (VARIETY_CONDITION_REPLICATE_READ.fq.gz), no SAMPLE_CONDITIONS needed:
        #
        # XX)
        #     SPECIES_DIR="00_5_NewName"
        #     SPECIES_NAME="Species name"
        #     GENOME_TYPE="carrot"
        #     ;;
        #
        # If your files use non-standard names, add SAMPLE_CONDITIONS:
        #
        # XX)
        #     SPECIES_DIR="00_5_NewName"
        #     SPECIES_NAME="Species name"
        #     GENOME_TYPE="carrot"
        #     SAMPLE_CONDITIONS="_L_=L _R_=R"
        #     ;;
        *)
            echo "ERROR: Unknown species code '${code}'"
            echo "Valid codes: DC, DG, MF, SK, DCDG"
            echo "To add a new code, edit SECTION 4 in scripts/config.sh"
            return 1
            ;;
    esac
}
