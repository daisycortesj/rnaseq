#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Build STAR index for different genomes
# USAGE:   sbatch run_genome_index.sbatch <genome_name>
#          where genome_name is: nutmeg, carrot, or custom
# GENOMES:
#          nutmeg = MYU_GWHGDHP00000000.1_genomic.fna (583M, no GTF)
#          carrot = GCF_001625215.2_DH1_v3.0_genomic.fna (426M, with GTF)
# EXAMPLES:
#          sbatch run_genome_index.sbatch nutmeg
#          sbatch run_genome_index.sbatch carrot
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=star_index
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/04_reference
#SBATCH -o star_index_%x_%j.out
#SBATCH -e star_index_%x_%j.err

set -euo pipefail

# ===== Configuration =====
GENOME_TYPE="${1:-nutmeg}"  # Default to nutmeg if no argument provided
REFERENCE_DIR="/projects/tholl_lab_1/daisy_analysis/04_reference"
CODE_DIR="/projects/tholl_lab_1/daisy_analysis/05_rnaseq-code"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "============================================================"
echo "STAR Index Building Pipeline"
echo "============================================================"
echo "SLURM job running on host: $(hostname)"
echo "CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-128G}  TIME: ${SLURM_JOB_TIMELIMIT:-unset}"
echo "Genome type: ${GENOME_TYPE}"
echo "Start time: $(date)"
echo "============================================================"

# ===== Activate Conda Environment =====
if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Warning: .bashrc has issues, using alternative conda activation"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate rnaseq
else
    echo "Loading .bashrc successfully"
    source ~/.bashrc
    conda activate rnaseq
fi

# Verify conda activation
conda info --envs | grep rnaseq || echo "ERROR: conda environment not activated"

# Check tools
echo "Checking tool versions:"
STAR --version || echo "STAR not found."
Trinity --version || echo "Trinity not found."
echo "------------------------------------------------------------"

# ===== Set Genome-Specific Parameters =====
case "${GENOME_TYPE}" in
    nutmeg)
        echo "Building STAR index for NUTMEG genome (MYU assembly, 583M, no GTF)"
        FASTA="${REFERENCE_DIR}/MYU_GWHGDHP00000000.1_genomic.fna"
        GTF=""  # No GTF for nutmeg
        OUTDIR="${REFERENCE_DIR}/nutmeg_star_index"
        READLEN=150
        ;;
    
    carrot)
        echo "Building STAR index for CARROT genome (GCF_001625215.2, 426M, with GTF)"
        FASTA="${REFERENCE_DIR}/GCF_001625215.2_DH1_v3.0_genomic.fna"
        GTF="${REFERENCE_DIR}/dc_genomic.gtf"
        OUTDIR="${REFERENCE_DIR}/carrot_star_index"
        READLEN=150
        ;;
    
    custom)
        echo "Using custom genome parameters from environment variables"
        # Set these before running: export CUSTOM_FASTA=..., CUSTOM_GTF=..., CUSTOM_OUTDIR=...
        FASTA="${CUSTOM_FASTA:-}"
        GTF="${CUSTOM_GTF:-}"
        OUTDIR="${CUSTOM_OUTDIR:-results/custom_star_index}"
        READLEN="${CUSTOM_READLEN:-150}"
        ;;
    
    *)
        echo "ERROR: Unknown genome type '${GENOME_TYPE}'"
        echo "Valid options: nutmeg, carrot, custom"
        exit 1
        ;;
esac

# ===== Validate Input Files =====
echo "Validating input files..."
if [[ -z "${FASTA}" ]]; then
    echo "ERROR: FASTA file not specified"
    exit 1
fi

if [[ ! -f "${FASTA}" ]]; then
    echo "ERROR: FASTA file not found: ${FASTA}"
    exit 1
fi

echo "✓ FASTA: ${FASTA}"

if [[ -n "${GTF}" ]]; then
    if [[ ! -f "${GTF}" ]]; then
        echo "ERROR: GTF file not found: ${GTF}"
        exit 1
    fi
    echo "✓ GTF: ${GTF}"
    GTF_ARG="--gtf ${GTF}"
else
    echo "✓ GTF: None (building index without annotations)"
    GTF_ARG=""
fi

echo "✓ Output directory: ${OUTDIR}"
echo "✓ Read length: ${READLEN}"
echo "------------------------------------------------------------"

# ===== Go to Pipeline Directory =====
cd "${CODE_DIR}"

# ===== Build STAR Command =====
CMD="python -m rna_pipeline.cli --fasta ${FASTA} ${GTF_ARG} --outdir ${OUTDIR} --threads ${THREADS} --readlen ${READLEN}"

echo "Running command:"
echo "${CMD}"
echo "------------------------------------------------------------"

# ===== Run Pipeline =====
eval ${CMD}

# ===== Summary =====
echo "============================================================"
echo "Pipeline completed successfully!"
echo "End time: $(date)"
echo "Output directory: ${OUTDIR}"
echo "============================================================"

# List output files
if [[ -d "${OUTDIR}" ]]; then
    echo "Generated files:"
    ls -lh "${OUTDIR}"
fi

