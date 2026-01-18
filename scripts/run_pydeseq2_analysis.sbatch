#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run PyDESeq2 differential expression analysis with CYP heatmap
# USAGE:   sbatch run_pydeseq2_analysis.sbatch
#          
#          # For Figure 6A-style CYP heatmap:
#          CYP_FAMILY_MAP=/path/to/cyp_families.tsv ROOT_UP_ONLY=true sbatch run_pydeseq2_analysis.sbatch
#
# INPUTS:  - Count matrix: 06_analysis/count_matrices/gene_count_matrix.tsv
#          - Metadata: 06_analysis/count_matrices/sample_metadata.tsv
#          - CYP family map (optional): TSV with gene_id, cyp_family columns
#
# OUTPUTS: 06_analysis/pydeseq2_results/
#          - pydeseq2_results.tsv (full DE results)
#          - deg_filtered.tsv (DEGs passing filters)
#          - cyp_deg_filtered.tsv (CYP DEGs only, if CYP_FAMILY_MAP provided)
#          - cyp_heatmap.pdf/png (Figure 6A-style heatmap)
#          - cyp_heatmap_matrix.tsv (heatmap values)
#          - *.pdf (QC and analysis plots)
#
# ENVIRONMENT VARIABLES (all optional):
#   DESIGN          - Design formula (default: auto-detect from metadata)
#   CONTRAST_FACTOR - Metadata column for contrast (default: condition)
#   CONTRAST_A      - Numerator condition, positive log2FC (default: root)
#   CONTRAST_B      - Denominator condition (default: leaf)
#   PADJ            - Adjusted p-value cutoff (default: 0.05)
#   LFC             - Log2 fold change cutoff (default: 2.0, set 0 to disable)
#   ROOT_UP_ONLY    - Keep only genes upregulated in CONTRAST_A (default: false)
#   CYP_FAMILY_MAP  - TSV file mapping gene_id to cyp_family (enables CYP heatmap)
#   CYP_GENES       - Text file with CYP gene IDs, one per line (alternative)
#   SAMPLE_ORDER    - Text file with sample names in desired column order
#   SCALE           - Scaling method: center or zscore (default: center)
#   ROW_CLUSTER     - Cluster rows: true/false (default: true)
#   COL_CLUSTER     - Cluster columns: true/false (default: false)
# ---------------------------------------------------------------------

# ===== SLURM directives =====
#SBATCH --job-name=pydeseq2
#SBATCH --account=tholl_lab_1
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daisycortesj@vt.edu
#SBATCH --chdir=/projects/tholl_lab_1/daisy_analysis
#SBATCH -o 06_analysis/pydeseq2_analysis_%j.out
#SBATCH -e 06_analysis/pydeseq2_analysis_%j.err

set -euo pipefail

# ===== Configuration =====
BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
COUNT_MATRIX="${BASE_DIR}/06_analysis/count_matrices/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/06_analysis/count_matrices/sample_metadata.tsv"
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_results"

# Design and contrast (optional - will auto-detect if not set)
DESIGN="${DESIGN:-}"
CONTRAST_FACTOR="${CONTRAST_FACTOR:-condition}"
CONTRAST_A="${CONTRAST_A:-root}"
CONTRAST_B="${CONTRAST_B:-leaf}"

# DEG filtering thresholds
PADJ="${PADJ:-0.05}"
LFC="${LFC:-2.0}"
ROOT_UP_ONLY="${ROOT_UP_ONLY:-false}"

# CYP gene filtering (provide one or the other for CYP heatmap)
CYP_FAMILY_MAP="${CYP_FAMILY_MAP:-}"
CYP_GENES="${CYP_GENES:-}"

# Heatmap options
SAMPLE_ORDER="${SAMPLE_ORDER:-}"
SCALE="${SCALE:-center}"
ROW_CLUSTER="${ROW_CLUSTER:-true}"
COL_CLUSTER="${COL_CLUSTER:-false}"

echo "============================================================"
echo "PyDESeq2 Differential Expression Analysis"
echo "============================================================"
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-4}  MEM: ${SLURM_MEM_PER_NODE:-32G}"
echo ""
echo "Input files:"
echo "  Count matrix: ${COUNT_MATRIX}"
echo "  Metadata: ${METADATA}"
echo ""
echo "Contrast settings:"
echo "  Factor: ${CONTRAST_FACTOR}"
echo "  Comparison: ${CONTRAST_A} vs ${CONTRAST_B}"
echo "  Design formula: ${DESIGN:-auto-detect}"
echo ""
echo "DEG filtering:"
echo "  padj cutoff: ${PADJ}"
echo "  log2FC cutoff: ${LFC}"
echo "  Root-up only: ${ROOT_UP_ONLY}"
echo ""
echo "CYP heatmap:"
if [ -n "${CYP_FAMILY_MAP}" ]; then
    echo "  CYP family map: ${CYP_FAMILY_MAP}"
elif [ -n "${CYP_GENES}" ]; then
    echo "  CYP gene list: ${CYP_GENES}"
else
    echo "  (No CYP mapping provided - CYP heatmap will be skipped)"
fi
echo "  Scale method: ${SCALE}"
echo "  Row clustering: ${ROW_CLUSTER}"
echo "  Column clustering: ${COL_CLUSTER}"
if [ -n "${SAMPLE_ORDER}" ]; then
    echo "  Sample order: ${SAMPLE_ORDER}"
fi
echo ""
echo "Output directory: ${OUTPUT_DIR}"
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

# Verify conda environment is activated
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "WARNING: Conda environment may not be activated properly"
    echo "Attempting to activate rnaseq environment..."
    conda activate rnaseq || {
        echo "ERROR: Failed to activate rnaseq conda environment"
        exit 1
    }
fi

echo "Active conda environment: ${CONDA_DEFAULT_ENV:-none}"

# ===== Check for required files =====
if [ ! -f "${COUNT_MATRIX}" ]; then
    echo "ERROR: Count matrix not found: ${COUNT_MATRIX}"
    echo "Please run build_count_matrix.py first or update COUNT_MATRIX path"
    exit 1
fi

if [ ! -f "${METADATA}" ]; then
    echo "ERROR: Metadata file not found: ${METADATA}"
    echo "Please run build_count_matrix.py first or update METADATA path"
    exit 1
fi

# Check for optional CYP files if specified
if [ -n "${CYP_FAMILY_MAP}" ] && [ ! -f "${CYP_FAMILY_MAP}" ]; then
    echo "ERROR: CYP family map file not found: ${CYP_FAMILY_MAP}"
    echo "Expected format (TSV):"
    echo "  gene_id    cyp_family    cyp_clan"
    echo "  LOC1234    CYP71         CYP71_clan"
    exit 1
fi

if [ -n "${CYP_GENES}" ] && [ ! -f "${CYP_GENES}" ]; then
    echo "ERROR: CYP genes file not found: ${CYP_GENES}"
    exit 1
fi

if [ -n "${SAMPLE_ORDER}" ] && [ ! -f "${SAMPLE_ORDER}" ]; then
    echo "ERROR: Sample order file not found: ${SAMPLE_ORDER}"
    exit 1
fi

# ===== Locate analysis script =====
SCRIPT_NAME="pydeseq2_analysis.py"

if [ -f "${BASE_DIR}/05_rnaseq-code/${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="${BASE_DIR}/05_rnaseq-code/${SCRIPT_NAME}"
elif [ -f "${BASE_DIR}/${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="${BASE_DIR}/${SCRIPT_NAME}"
elif [ -f "$(dirname "$0")/../${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="$(dirname "$0")/../${SCRIPT_NAME}"
else
    echo "ERROR: ${SCRIPT_NAME} not found"
    echo "Searched in:"
    echo "  - ${BASE_DIR}/05_rnaseq-code/${SCRIPT_NAME}"
    echo "  - ${BASE_DIR}/${SCRIPT_NAME}"
    echo "  - $(dirname "$0")/../${SCRIPT_NAME}"
    exit 1
fi
echo "Using analysis script: ${ANALYSIS_SCRIPT}"

# ===== Check and install dependencies =====
echo ""
echo "Checking for PyDESeq2 and dependencies..."
echo "Python version: $(python3 --version 2>&1)"
echo "Python path: $(which python3 2>&1)"

# Quick check for required packages
PACKAGES_MISSING=false

for pkg in pydeseq2 matplotlib seaborn pandas numpy scipy; do
    if pip list 2>/dev/null | grep -qi "^${pkg} "; then
        echo "  ✓ ${pkg} found"
    else
        echo "  ✗ ${pkg} not found"
        PACKAGES_MISSING=true
    fi
done

if [ "${PACKAGES_MISSING}" = true ]; then
    echo ""
    echo "Installing missing packages..."
    pip install pydeseq2 matplotlib seaborn pandas numpy scipy 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to install required packages"
        exit 1
    fi
    echo "✓ Packages installed successfully"
fi

# Verify imports work
echo ""
echo "Testing imports..."
timeout 120 python3 -c "
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
print('✓ All required imports successful')
" 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Package imports failed"
    exit 1
fi

# ===== Create output directory =====
mkdir -p "${OUTPUT_DIR}"

# ===== Build command line arguments =====
echo ""
echo "============================================================"
echo "Running PyDESeq2 analysis..."
echo "============================================================"

# Start building the command
CMD="python3 \"${ANALYSIS_SCRIPT}\" \"${COUNT_MATRIX}\" \"${METADATA}\" -o \"${OUTPUT_DIR}\""

# Add design formula if specified
if [ -n "${DESIGN}" ]; then
    CMD="${CMD} --design \"${DESIGN}\""
fi

# Add contrast settings
CMD="${CMD} --contrast-factor \"${CONTRAST_FACTOR}\""
CMD="${CMD} --contrast-A \"${CONTRAST_A}\""
CMD="${CMD} --contrast-B \"${CONTRAST_B}\""

# Add DEG filtering
CMD="${CMD} --padj ${PADJ}"
CMD="${CMD} --lfc ${LFC}"

if [ "${ROOT_UP_ONLY}" = "true" ] || [ "${ROOT_UP_ONLY}" = "True" ] || [ "${ROOT_UP_ONLY}" = "1" ]; then
    CMD="${CMD} --root-up-only"
fi

# Add CYP gene mapping
if [ -n "${CYP_FAMILY_MAP}" ]; then
    CMD="${CMD} --cyp-family-map \"${CYP_FAMILY_MAP}\""
elif [ -n "${CYP_GENES}" ]; then
    CMD="${CMD} --cyp-genes \"${CYP_GENES}\""
fi

# Add sample ordering
if [ -n "${SAMPLE_ORDER}" ]; then
    CMD="${CMD} --sample-order \"${SAMPLE_ORDER}\""
fi

# Add heatmap options
CMD="${CMD} --scale ${SCALE}"

if [ "${ROW_CLUSTER}" = "false" ] || [ "${ROW_CLUSTER}" = "False" ] || [ "${ROW_CLUSTER}" = "0" ]; then
    CMD="${CMD} --no-row-cluster"
fi

if [ "${COL_CLUSTER}" = "true" ] || [ "${COL_CLUSTER}" = "True" ] || [ "${COL_CLUSTER}" = "1" ]; then
    CMD="${CMD} --col-cluster"
fi

# Print and execute command
echo ""
echo "Command:"
echo "${CMD}"
echo ""

eval ${CMD}

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: PyDESeq2 analysis failed with exit code ${EXIT_CODE}"
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════"
echo "║ ANALYSIS COMPLETE"
echo "╚════════════════════════════════════════════════════════════"
echo ""
echo "Results directory: ${OUTPUT_DIR}"
echo ""
echo "Output files:"
echo "  Core DE results:"
echo "    - pydeseq2_results.tsv (full differential expression results)"
echo "    - deg_filtered.tsv (DEGs passing padj/lfc filters)"
echo "    - analysis_summary.txt (summary report)"
echo ""
echo "  Plots:"
echo "    - pydeseq2_ma_plot.pdf (MA plot)"
echo "    - pydeseq2_volcano_plot.pdf (Volcano plot)"
echo "    - qc_total_counts.pdf (QC: total counts per sample)"
echo "    - qc_gene_counts.pdf (QC: gene count distributions)"

if [ -n "${CYP_FAMILY_MAP}" ] || [ -n "${CYP_GENES}" ]; then
    echo ""
    echo "  CYP heatmap (Figure 6A style):"
    echo "    - cyp_deg_filtered.tsv (CYP DEGs with family annotations)"
    echo "    - cyp_heatmap_matrix.tsv (transformed values for plotting)"
    echo "    - cyp_heatmap.pdf (high-quality heatmap)"
    echo "    - cyp_heatmap.png (for quick viewing)"
fi

echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0
