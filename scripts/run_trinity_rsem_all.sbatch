#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Quantify expression with RSEM using Trinity-assembled transcripts
# USAGE:   sbatch run_trinity_rsem_all.sbatch
# INPUTS:  - Trinity assemblies in 01_processed/*_trinity_assembly/Trinity.fasta
#          - Raw FASTQs in 00_rawdata/00_*/<sample>_{1,2}.fq.gz
# OUTPUTS: 03_count_tables/<group>_trinity_rsem/<sample>.genes.results
# NOTES:   Uses Trinity's helper 'align_and_estimate_abundance.pl' (RSEM+bowtie2)
#          Set SS_lib_type to match library: FR (forward) or RF (reverse); omit for none
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=trinity_rsem_all
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 03_count_tables/trinity_rsem_all_%j.out
#SBATCH -e 03_count_tables/trinity_rsem_all_%j.err

set -euo pipefail

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
RAW_DIR="${BASE_DIR}/00_rawdata"
PROC_DIR="${BASE_DIR}/01_processed"
OUT_BASE="${BASE_DIR}/03_count_tables"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# SS_lib_type: FR (forward), RF (reverse), or empty for unstranded
SS_LIB_TYPE="${SS_LIB_TYPE:-}"  # e.g., RF for Illumina TruSeq stranded

echo "============================================================"
echo "Trinity-guided RSEM quantification for all samples"
echo "Host: $(hostname)  CPUs: ${THREADS}  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo "SS_lib_type: ${SS_LIB_TYPE:-unset}"
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

which Trinity || { echo "ERROR: Trinity not found"; exit 1; }
which bowtie2 || { echo "ERROR: bowtie2 not found"; exit 1; }
which rsem-calculate-expression || { echo "ERROR: RSEM not found"; exit 1; }

mkdir -p "${OUT_BASE}"

TOTAL=0
OK=0
FAIL=0

# Helper: run Trinity's abundance estimation per sample
run_trinity_rsem() {
    local transcripts="$1"   # Trinity.fasta
    local left_fq="$2"
    local right_fq="$3"
    local out_dir="$4"

    mkdir -p "${out_dir}"

    local ss_arg=()
    if [[ -n "${SS_LIB_TYPE}" ]]; then
        ss_arg=("--SS_lib_type" "${SS_LIB_TYPE}")
    fi

    align_and_estimate_abundance.pl \
        --transcripts "${transcripts}" \
        --seqType fq \
        --left "${left_fq}" \
        --right "${right_fq}" \
        --est_method RSEM \
        --aln_method bowtie2 \
        --thread_count "${THREADS}" \
        --prep_reference \
        "${ss_arg[@]}" \
        --output_dir "${out_dir}"
}

# Iterate over Trinity assemblies in 01_processed/*_trinity_assembly
shopt -s nullglob
for asm_dir in "${PROC_DIR}"/*_trinity_assembly; do
    [[ -d "${asm_dir}" ]] || continue
    transcripts="${asm_dir}/Trinity.fasta"
    [[ -f "${transcripts}" ]] || { echo "No Trinity.fasta in ${asm_dir}, skipping"; continue; }

    # Try to infer sample name and group from assembly directory name
    asm_name=$(basename "${asm_dir}")    # e.g., DC1L1_trinity_assembly
    sample_base=${asm_name%_trinity_assembly}

    # Group determination based on raw data layout (00_1_DC, 00_2_DG, 00_3_MF)
    # Search for matching FASTQs in any group
    found=0
    for group in 00_1_DC 00_2_DG 00_3_MF; do
        left_fq="${RAW_DIR}/${group}/${sample_base}_1.fq.gz"
        right_fq="${RAW_DIR}/${group}/${sample_base}_2.fq.gz"
        if [[ -f "${left_fq}" && -f "${right_fq}" ]]; then
            out_dir="${OUT_BASE}/${group}_trinity_rsem/${sample_base}"
            echo "Running Trinity RSEM: ${group}/${sample_base}"
            TOTAL=$((TOTAL + 1))
            set +e
            run_trinity_rsem "${transcripts}" "${left_fq}" "${right_fq}" "${out_dir}"
            exit_code=$?
            set -e
            if [[ ${exit_code} -eq 0 ]]; then
                echo "  ✓ ${sample_base}"
                OK=$((OK + 1))
            else
                echo "  ✗ ${sample_base} (exit ${exit_code})"
                FAIL=$((FAIL + 1))
            fi
            found=1
            break
        fi
    done

    if [[ ${found} -eq 0 ]]; then
        echo "Could not locate FASTQs for ${sample_base} in 00_rawdata/*/, skipping"
    fi
done
shopt -u nullglob

echo ""
echo "╔════════════════════════════════════════════════════════════"
echo "║ FINAL SUMMARY"
echo "╚════════════════════════════════════════════════════════════"
echo "Total samples attempted: ${TOTAL}"
echo "Successful: ${OK}   Failed: ${FAIL}"
echo "Outputs in: ${OUT_BASE}/*_trinity_rsem/*/*.genes.results"
echo "End time: $(date)"
echo "============================================================"

exit 0


