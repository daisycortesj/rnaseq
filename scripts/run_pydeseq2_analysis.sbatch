#!/bin/bash
# ---------------------------------------------------------------------
# PURPOSE: Run PyDESeq2 differential expression analysis
# USAGE:   sbatch run_pydeseq2_analysis.sbatch
# INPUTS:  - Count matrix: 06_analysis/count_matrices/gene_count_matrix.tsv
#          - Metadata: 06_analysis/count_matrices/sample_metadata.tsv
# OUTPUTS: 06_analysis/pydeseq2_results/pydeseq2_results.tsv
#         06_analysis/pydeseq2_results/*.pdf (plots)
# NOTES:   Performs differential expression analysis using PyDESeq2
#          Set DESIGN environment variable to override auto-detection
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

BASE_DIR="/projects/tholl_lab_1/daisy_analysis"
COUNT_MATRIX="${BASE_DIR}/06_analysis/count_matrices/gene_count_matrix.tsv"
METADATA="${BASE_DIR}/06_analysis/count_matrices/sample_metadata.tsv"
OUTPUT_DIR="${BASE_DIR}/06_analysis/pydeseq2_results"

# Design formula (optional - will auto-detect if not set)
# Examples: "treatment", "group", "condition", "treatment+condition"
DESIGN="${DESIGN:-}"

echo "============================================================"
echo "PyDESeq2 Differential Expression Analysis"
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-4}  MEM: ${SLURM_MEM_PER_NODE:-32G}"
echo "Count matrix: ${COUNT_MATRIX}"
echo "Metadata: ${METADATA}"
echo "Output directory: ${OUTPUT_DIR}"
if [ -n "${DESIGN}" ]; then
    echo "Design formula: ${DESIGN}"
else
    echo "Design formula: auto-detect"
fi
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
if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "WARNING: Conda environment may not be activated properly"
    echo "Attempting to activate rnaseq environment..."
    conda activate rnaseq || {
        echo "ERROR: Failed to activate rnaseq conda environment"
        exit 1
    }
fi

echo "Active conda environment: ${CONDA_DEFAULT_ENV:-none}"

# Check for required files
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

# Check for pydeseq2_analysis.py script
# Try multiple possible locations (scripts are in 05_rnaseq-code/)
if [ -f "${BASE_DIR}/05_rnaseq-code/pydeseq2_analysis.py" ]; then
    ANALYSIS_SCRIPT="${BASE_DIR}/05_rnaseq-code/pydeseq2_analysis.py"
elif [ -f "${BASE_DIR}/pydeseq2_analysis.py" ]; then
    ANALYSIS_SCRIPT="${BASE_DIR}/pydeseq2_analysis.py"
elif [ -f "${SCRIPT_DIR}/../pydeseq2_analysis.py" ]; then
    ANALYSIS_SCRIPT="${SCRIPT_DIR}/../pydeseq2_analysis.py"
else
    echo "ERROR: pydeseq2_analysis.py not found"
    echo "Searched in:"
    echo "  - ${BASE_DIR}/05_rnaseq-code/pydeseq2_analysis.py"
    echo "  - ${BASE_DIR}/pydeseq2_analysis.py"
    echo "  - ${SCRIPT_DIR}/../pydeseq2_analysis.py"
    exit 1
fi
echo "Using analysis script: ${ANALYSIS_SCRIPT}"

# Check and install PyDESeq2 if needed
echo "Checking for PyDESeq2 and dependencies..."
echo "Python version: $(python3 --version 2>&1)"
echo "Python path: $(which python3 2>&1)"
echo "Pip path: $(which pip 2>&1)"
echo ""

# Quick check: see if packages are installed (faster than import)
echo "Checking if packages are installed..."
PYDESEQ2_INSTALLED=1
MATPLOTLIB_INSTALLED=1
SEABORN_INSTALLED=1

# Check using pip list (faster than import)
if pip list 2>/dev/null | grep -q "^pydeseq2 "; then
    PYDESEQ2_INSTALLED=0
    echo "  ✓ PyDESeq2 package found"
else
    echo "  ✗ PyDESeq2 package not found"
fi

if pip list 2>/dev/null | grep -q "^matplotlib "; then
    MATPLOTLIB_INSTALLED=0
    echo "  ✓ matplotlib package found"
else
    echo "  ✗ matplotlib package not found"
fi

if pip list 2>/dev/null | grep -q "^seaborn "; then
    SEABORN_INSTALLED=0
    echo "  ✓ seaborn package found"
else
    echo "  ✗ seaborn package not found"
fi

# If packages are installed, do a quick import test (with timeout)
if [ ${PYDESEQ2_INSTALLED} -eq 0 ] && [ ${MATPLOTLIB_INSTALLED} -eq 0 ] && [ ${SEABORN_INSTALLED} -eq 0 ]; then
    echo ""
    echo "Testing imports (this may take a moment on first import)..."
    timeout 60 python3 -c "import pydeseq2; import matplotlib; import seaborn; print('All imports successful')" 2>&1
    if [ $? -ne 0 ]; then
        echo "WARNING: Packages installed but import failed. Will try to reinstall."
        PYDESEQ2_INSTALLED=1
        MATPLOTLIB_INSTALLED=1
        SEABORN_INSTALLED=1
    else
        echo "✓ All packages import successfully"
    fi
fi

if [ ${PYDESEQ2_INSTALLED} -ne 0 ] || [ ${MATPLOTLIB_INSTALLED} -ne 0 ] || [ ${SEABORN_INSTALLED} -ne 0 ]; then
    echo ""
    echo "Missing required packages. Installing..."
    echo "This may take a few minutes (PyDESeq2 is a large package)..."
    
    # Install missing packages
    MISSING_PACKAGES=""
    [ ${PYDESEQ2_INSTALLED} -ne 0 ] && MISSING_PACKAGES="${MISSING_PACKAGES} pydeseq2"
    [ ${MATPLOTLIB_INSTALLED} -ne 0 ] && MISSING_PACKAGES="${MISSING_PACKAGES} matplotlib"
    [ ${SEABORN_INSTALLED} -ne 0 ] && MISSING_PACKAGES="${MISSING_PACKAGES} seaborn"
    
    echo "Installing: ${MISSING_PACKAGES}"
    echo "Note: PyDESeq2 installation can take 5-10 minutes..."
    
    # Install with progress output
    pip install ${MISSING_PACKAGES} 2>&1 | tee /tmp/pip_install.log
    
    INSTALL_EXIT=$?
    if [ ${INSTALL_EXIT} -ne 0 ]; then
        echo ""
        echo "ERROR: Failed to install required packages: ${MISSING_PACKAGES}"
        echo "Exit code: ${INSTALL_EXIT}"
        echo "Last 20 lines of pip output:"
        tail -20 /tmp/pip_install.log
        echo ""
        echo "Please install manually: pip install pydeseq2 matplotlib seaborn"
        exit 1
    fi
    
    echo ""
    echo "Verifying installation (this may take a moment)..."
    # Verify installation with longer timeout (PyDESeq2 first import can be slow)
    timeout 120 python3 -c "import pydeseq2, matplotlib, seaborn; print('✓ All packages imported successfully')" 2>&1
    VERIFY_EXIT=$?
    if [ ${VERIFY_EXIT} -ne 0 ]; then
        if [ ${VERIFY_EXIT} -eq 124 ]; then
            echo "WARNING: Import test timed out (PyDESeq2 may be slow on first import)"
            echo "Continuing anyway - packages should work during analysis"
        else
            echo "ERROR: Installation completed but packages still not importable"
            echo "Try importing manually to see the error:"
            echo "  python3 -c 'import pydeseq2'"
            exit 1
        fi
    fi
    
    echo "✓ Successfully installed: ${MISSING_PACKAGES}"
else
    echo ""
    echo "✓ All required packages are installed (PyDESeq2, matplotlib, seaborn)"
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# ===== Run PyDESeq2 Analysis =====
echo ""
echo "Running PyDESeq2 analysis..."

if [ -n "${DESIGN}" ]; then
    python3 "${ANALYSIS_SCRIPT}" \
        "${COUNT_MATRIX}" \
        "${METADATA}" \
        -o "${OUTPUT_DIR}" \
        --design "${DESIGN}"
else
    python3 "${ANALYSIS_SCRIPT}" \
        "${COUNT_MATRIX}" \
        "${METADATA}" \
        -o "${OUTPUT_DIR}"
fi

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: PyDESeq2 analysis failed with exit code ${EXIT_CODE}"
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════"
echo "║ FINAL SUMMARY"
echo "╚════════════════════════════════════════════════════════════"
echo "Results directory: ${OUTPUT_DIR}"
echo ""
echo "Output files:"
echo "  - pydeseq2_results.tsv (differential expression results)"
echo "  - pydeseq2_ma_plot.pdf (MA plot)"
echo "  - pydeseq2_volcano_plot.pdf (Volcano plot)"
echo "  - qc_total_counts.pdf (Quality control plots)"
echo "  - qc_gene_counts.pdf (Gene count distributions)"
echo "  - heatmap_top_variable_genes.pdf (Heatmap)"
echo "  - analysis_summary.txt (Summary report)"
echo ""
echo "End time: $(date)"
echo "============================================================"

exit 0

