#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run Trinity for missing samples only
# USAGE:   sbatch run_trinity_missing.sbatch
# NOTES:   Only processes samples that don't have assemblies yet
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_missing
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis/01_processed
#SBATCH -o trinity_missing_%j.out
#SBATCH -e trinity_missing_%j.err

set -euo pipefail

echo "============================================================"
echo "Trinity assembly for missing samples"
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

# Verify tools are available
which Trinity || { echo "ERROR: Trinity not found"; exit 1; }
which java || { echo "ERROR: Java not found"; exit 1; }

# ===== Missing DG samples =====
echo "Processing missing DG samples..."
for sample in DGL1 DGL3 DGR1 DGR2 DGR3; do
    if [[ ! -d "${sample}_trinity_assembly" ]]; then
        echo "Running Trinity for $sample"
        Trinity \
            --seqType fq \
            --left /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_2_DG/${sample}_1.fq.gz \
            --right /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_2_DG/${sample}_2.fq.gz \
            --CPU 32 \
            --max_memory 120G \
            --min_kmer_cov 2 \
            --normalize_reads \
            --output ${sample}_trinity_assembly
        echo "Completed: $sample"
    else
        echo "Skipping $sample (assembly already exists)"
    fi
done

# ===== Missing MF samples =====
echo "Processing missing MF samples..."
for sample in MFF2 MFF3 MFL1 MFL2 MFL3; do
    if [[ ! -d "${sample}_trinity_assembly" ]]; then
        echo "Running Trinity for $sample"
        Trinity \
            --seqType fq \
            --left /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_3_MF/${sample}_1.fq.gz \
            --right /projects/tholl_lab_1/daisy_analysis/00_rawdata/00_3_MF/${sample}_2.fq.gz \
            --CPU 32 \
            --max_memory 120G \
            --min_kmer_cov 2 \
            --normalize_reads \
            --output ${sample}_trinity_assembly
        echo "Completed: $sample"
    else
        echo "Skipping $sample (assembly already exists)"
    fi
done

echo "============================================================"
echo "Trinity assembly completed"
echo "End time: $(date)"
echo "============================================================"
