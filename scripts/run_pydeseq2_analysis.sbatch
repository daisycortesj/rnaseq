#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run PyDESeq2 differential expression analysis with CYP heatmap
#
# USAGE:
#   # Basic DE analysis (no CYP heatmap)
#   sbatch run_pydeseq2_analysis.sbatch
#
#   # CYP heatmap (auto-extract CYP genes from GTF)
#   GTF_FILE=/path/to/genomic.gtf ROOT_UP_ONLY=true sbatch run_pydeseq2_analysis.sbatch
#
#   # Or provide pre-made CYP family map
#   CYP_FAMILY_MAP=/path/to/cyp_families.tsv ROOT_UP_ONLY=true sbatch run_pydeseq2_analysis.sbatch
#
# INPUTS:
#   Required:
#   - Count matrix: 06_analysis/count_matrices/gene_count_matrix.tsv
#   - Metadata: 06_analysis/count_matrices/sample_metadata.tsv
#
#   For CYP heatmap (one of):
#   - GTF_FILE: Path to GTF annotation (CYP families extracted automatically)
#   - CYP_FAMILY_MAP: Pre-made TSV with gene_id, cyp_family columns
#
# OUTPUTS: 06_analysis/pydeseq2_results/
#   Core DE:
#   - pydeseq2_results.tsv (full DE results)
#   - deg_filtered.tsv (DEGs passing filters)
#   - analysis_summary.txt
#   - MA plot, volcano plot, QC plots
#
#   CYP heatmap (if GTF_FILE or CYP_FAMILY_MAP provided):
#   - cyp_families.tsv (auto-generated from GTF)
#   - cyp_deg_filtered.tsv (CYP DEGs)
#   - cyp_heatmap.pdf/png (CYP gene expression heatmap)
#   - cyp_heatmap_matrix.tsv
#
# ENVIRONMENT VARIABLES:
#   GTF_FILE        - GTF annotation file (auto-extracts CYP families)
#   CYP_FAMILY_MAP  - Pre-made CYP family mapping (skips auto-extraction)
#   CONTRAST_FACTOR - Metadata column for contrast (default: condition)
#   CONTRAST_A      - Numerator condition (default: root)
#   CONTRAST_B      - Denominator condition (default: leaf)
#   PADJ            - Adjusted p-value cutoff (default: 0.05)
#   LFC             - Log2 fold change cutoff (default: 2.0)
#   ROOT_UP_ONLY    - Only root-upregulated genes (default: false)
#   SCALE           - Scaling: center or zscore (default: center)
#   ROW_CLUSTER     - Cluster rows (default: true)
#   COL_CLUSTER     - Cluster columns (default: false)
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

# ╔════════════════════════════════════════════════════════════════════╗
# ║                         CONFIGURATION                               ║
# ╚════════════════════════════════════════════════════════════════════╝

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
CODE_DIR="${BASE_DIR}/05_rnaseq-code"

# Input files
COUNT_MATRIX="${BASE_DIR}/06_analysis/count_matrices/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/06_analysis/count_matrices/sample_metadata.tsv"
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_results"

# GTF file for automatic CYP extraction (set this OR CYP_FAMILY_MAP)
GTF_FILE="${GTF_FILE:-${BASE_DIR}/04_reference/dc_genomic.gtf}"

# Contrast settings
DESIGN="${DESIGN:-}"
CONTRAST_FACTOR="${CONTRAST_FACTOR:-condition}"
CONTRAST_A="${CONTRAST_A:-root}"
CONTRAST_B="${CONTRAST_B:-leaf}"

# DEG filtering thresholds
PADJ="${PADJ:-0.05}"
LFC="${LFC:-2.0}"
ROOT_UP_ONLY="${ROOT_UP_ONLY:-false}"

# CYP gene filtering
# If CYP_FAMILY_MAP is set, use it directly
# If not set but GTF_FILE exists, auto-extract
CYP_FAMILY_MAP="${CYP_FAMILY_MAP:-}"
CYP_GENES="${CYP_GENES:-}"

# Heatmap options
SAMPLE_ORDER="${SAMPLE_ORDER:-}"
SCALE="${SCALE:-center}"
ROW_CLUSTER="${ROW_CLUSTER:-true}"
COL_CLUSTER="${COL_CLUSTER:-false}"

# ╔════════════════════════════════════════════════════════════════════╗
# ║                         PRINT HEADER                                ║
# ╚════════════════════════════════════════════════════════════════════╝

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║        PyDESeq2 Differential Expression Analysis                   ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-4}  MEM: ${SLURM_MEM_PER_NODE:-32G}"
echo "Start time: $(date)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "INPUT FILES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Count matrix: ${COUNT_MATRIX}"
echo "  Metadata:     ${METADATA}"
echo "  GTF file:     ${GTF_FILE}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "CONTRAST SETTINGS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Factor:       ${CONTRAST_FACTOR}"
echo "  Comparison:   ${CONTRAST_A} vs ${CONTRAST_B}"
echo "  Design:       ${DESIGN:-auto-detect}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DEG FILTERING"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  padj cutoff:  ${PADJ}"
echo "  log2FC cutoff: ${LFC}"
echo "  Root-up only: ${ROOT_UP_ONLY}"
echo ""

# ╔════════════════════════════════════════════════════════════════════╗
# ║                    ACTIVATE CONDA ENVIRONMENT                       ║
# ╚════════════════════════════════════════════════════════════════════╝

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ACTIVATING ENVIRONMENT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if ! bash -c "source ~/.bashrc" 2>/dev/null; then
    echo "Using alternative conda activation..."
    source ~/miniconda3/etc/profile.d/conda.sh
else
    source ~/.bashrc
fi
conda activate rnaseq

if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "ERROR: Failed to activate rnaseq conda environment"
    exit 1
fi
echo "✓ Conda environment: ${CONDA_DEFAULT_ENV}"

# ╔════════════════════════════════════════════════════════════════════╗
# ║                      VALIDATE INPUT FILES                           ║
# ╚════════════════════════════════════════════════════════════════════╝

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "VALIDATING FILES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -f "${COUNT_MATRIX}" ]; then
    echo "ERROR: Count matrix not found: ${COUNT_MATRIX}"
    echo "Run build_count_matrix.py first."
    exit 1
fi
echo "✓ Count matrix found"

if [ ! -f "${METADATA}" ]; then
    echo "ERROR: Metadata not found: ${METADATA}"
    exit 1
fi
echo "✓ Metadata found"

# ╔════════════════════════════════════════════════════════════════════╗
# ║              AUTO-EXTRACT CYP FAMILIES FROM GTF                     ║
# ╚════════════════════════════════════════════════════════════════════╝

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "CYP HEATMAP SETUP"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Determine CYP family map source
if [ -n "${CYP_FAMILY_MAP}" ] && [ -f "${CYP_FAMILY_MAP}" ]; then
    # User provided a pre-made CYP family map
    echo "Using provided CYP family map: ${CYP_FAMILY_MAP}"
    
elif [ -f "${GTF_FILE}" ]; then
    # Auto-extract from GTF
    echo "GTF file found: ${GTF_FILE}"
    echo "Auto-extracting CYP gene families..."
    
    # Set output path for extracted CYP families
    CYP_FAMILY_MAP="${OUTPUT_DIR}/cyp_families.tsv"
    mkdir -p "${OUTPUT_DIR}"
    
    # Find the extraction script
    EXTRACT_SCRIPT=""
    if [ -f "${CODE_DIR}/extract_cyp_families.py" ]; then
        EXTRACT_SCRIPT="${CODE_DIR}/extract_cyp_families.py"
    elif [ -f "$(dirname "$0")/../extract_cyp_families.py" ]; then
        EXTRACT_SCRIPT="$(dirname "$0")/../extract_cyp_families.py"
    fi
    
    if [ -n "${EXTRACT_SCRIPT}" ]; then
        echo "Running: python3 ${EXTRACT_SCRIPT}"
        python3 "${EXTRACT_SCRIPT}" "${GTF_FILE}" -o "${CYP_FAMILY_MAP}"
        
        if [ $? -eq 0 ] && [ -f "${CYP_FAMILY_MAP}" ]; then
            CYP_COUNT=$(wc -l < "${CYP_FAMILY_MAP}")
            echo "✓ Extracted $((CYP_COUNT - 1)) CYP genes to: ${CYP_FAMILY_MAP}"
        else
            echo "WARNING: CYP extraction failed. Continuing without CYP heatmap."
            CYP_FAMILY_MAP=""
        fi
    else
        echo "WARNING: extract_cyp_families.py not found. Skipping CYP extraction."
        CYP_FAMILY_MAP=""
    fi
else
    echo "No GTF_FILE or CYP_FAMILY_MAP provided."
    echo "CYP heatmap will be skipped."
    echo ""
    echo "To enable CYP heatmap, set one of:"
    echo "  GTF_FILE=/path/to/genomic.gtf"
    echo "  CYP_FAMILY_MAP=/path/to/cyp_families.tsv"
fi

# Show heatmap settings if CYP map is available
if [ -n "${CYP_FAMILY_MAP}" ]; then
    echo ""
    echo "CYP heatmap settings:"
    echo "  Scale method:     ${SCALE}"
    echo "  Row clustering:   ${ROW_CLUSTER}"
    echo "  Column clustering: ${COL_CLUSTER}"
    if [ -n "${SAMPLE_ORDER}" ]; then
        echo "  Sample order:     ${SAMPLE_ORDER}"
    fi
fi

# ╔════════════════════════════════════════════════════════════════════╗
# ║                      LOCATE ANALYSIS SCRIPT                         ║
# ╚════════════════════════════════════════════════════════════════════╝

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "LOCATING SCRIPTS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

SCRIPT_NAME="pydeseq2_analysis.py"

if [ -f "${CODE_DIR}/${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="${CODE_DIR}/${SCRIPT_NAME}"
elif [ -f "${BASE_DIR}/${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="${BASE_DIR}/${SCRIPT_NAME}"
elif [ -f "$(dirname "$0")/../${SCRIPT_NAME}" ]; then
    ANALYSIS_SCRIPT="$(dirname "$0")/../${SCRIPT_NAME}"
else
    echo "ERROR: ${SCRIPT_NAME} not found"
    exit 1
fi
echo "✓ Analysis script: ${ANALYSIS_SCRIPT}"

# ╔════════════════════════════════════════════════════════════════════╗
# ║                    CHECK PYTHON DEPENDENCIES                        ║
# ╚════════════════════════════════════════════════════════════════════╝

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "CHECKING DEPENDENCIES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

PACKAGES_MISSING=false
for pkg in pydeseq2 matplotlib seaborn pandas numpy scipy; do
    if pip list 2>/dev/null | grep -qi "^${pkg} "; then
        echo "  ✓ ${pkg}"
    else
        echo "  ✗ ${pkg} (missing)"
        PACKAGES_MISSING=true
    fi
done

if [ "${PACKAGES_MISSING}" = true ]; then
    echo ""
    echo "Installing missing packages..."
    pip install pydeseq2 matplotlib seaborn pandas numpy scipy 2>&1
fi

# Quick import test
echo ""
echo "Testing imports..."
timeout 120 python3 -c "
import pandas, numpy, matplotlib, seaborn, scipy
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
print('✓ All imports successful')
" 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Package imports failed"
    exit 1
fi

# ╔════════════════════════════════════════════════════════════════════╗
# ║                       RUN PYDESEQ2 ANALYSIS                         ║
# ╚════════════════════════════════════════════════════════════════════╝

mkdir -p "${OUTPUT_DIR}"

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    RUNNING PYDESEQ2 ANALYSIS                       ║"
echo "╚════════════════════════════════════════════════════════════════════╝"

# Build command
CMD="python3 \"${ANALYSIS_SCRIPT}\" \"${COUNT_MATRIX}\" \"${METADATA}\" -o \"${OUTPUT_DIR}\""

# Design formula
if [ -n "${DESIGN}" ]; then
    CMD="${CMD} --design \"${DESIGN}\""
fi

# Contrast settings
CMD="${CMD} --contrast-factor \"${CONTRAST_FACTOR}\""
CMD="${CMD} --contrast-A \"${CONTRAST_A}\""
CMD="${CMD} --contrast-B \"${CONTRAST_B}\""

# DEG filtering
CMD="${CMD} --padj ${PADJ}"
CMD="${CMD} --lfc ${LFC}"

if [ "${ROOT_UP_ONLY}" = "true" ] || [ "${ROOT_UP_ONLY}" = "True" ] || [ "${ROOT_UP_ONLY}" = "1" ]; then
    CMD="${CMD} --root-up-only"
fi

# CYP gene mapping
if [ -n "${CYP_FAMILY_MAP}" ] && [ -f "${CYP_FAMILY_MAP}" ]; then
    CMD="${CMD} --cyp-family-map \"${CYP_FAMILY_MAP}\""
elif [ -n "${CYP_GENES}" ] && [ -f "${CYP_GENES}" ]; then
    CMD="${CMD} --cyp-genes \"${CYP_GENES}\""
fi

# Sample ordering
if [ -n "${SAMPLE_ORDER}" ] && [ -f "${SAMPLE_ORDER}" ]; then
    CMD="${CMD} --sample-order \"${SAMPLE_ORDER}\""
fi

# Heatmap options
CMD="${CMD} --scale ${SCALE}"

if [ "${ROW_CLUSTER}" = "false" ] || [ "${ROW_CLUSTER}" = "False" ] || [ "${ROW_CLUSTER}" = "0" ]; then
    CMD="${CMD} --no-row-cluster"
fi

if [ "${COL_CLUSTER}" = "true" ] || [ "${COL_CLUSTER}" = "True" ] || [ "${COL_CLUSTER}" = "1" ]; then
    CMD="${CMD} --col-cluster"
fi

# Execute
echo ""
echo "Command:"
echo "${CMD}"
echo ""

eval ${CMD}
EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo ""
    echo "ERROR: PyDESeq2 analysis failed with exit code ${EXIT_CODE}"
    exit 1
fi

# ╔════════════════════════════════════════════════════════════════════╗
# ║                         SUMMARY                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                      ANALYSIS COMPLETE                             ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Results: ${OUTPUT_DIR}/"
echo ""
echo "Core DE outputs:"
echo "  • pydeseq2_results.tsv      - Full differential expression results"
echo "  • deg_filtered.tsv          - DEGs passing padj/lfc filters"
echo "  • pydeseq2_ma_plot.pdf      - MA plot"
echo "  • pydeseq2_volcano_plot.pdf - Volcano plot"
echo "  • qc_*.pdf                  - Quality control plots"
echo "  • analysis_summary.txt      - Summary report"

if [ -n "${CYP_FAMILY_MAP}" ] && [ -f "${CYP_FAMILY_MAP}" ]; then
    echo ""
    echo "CYP heatmap outputs:"
    echo "  • cyp_families.tsv          - CYP gene family mapping (auto-generated)"
    echo "  • cyp_deg_filtered.tsv      - CYP DEGs with annotations"
    echo "  • cyp_heatmap.pdf           - High-quality heatmap"
    echo "  • cyp_heatmap.png           - Quick-view heatmap"
    echo "  • cyp_heatmap_matrix.tsv    - Transformed values"
fi

echo ""
echo "End time: $(date)"
echo "════════════════════════════════════════════════════════════════════"

exit 0
