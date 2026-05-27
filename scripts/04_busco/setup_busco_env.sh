#!/bin/bash
# =============================================================================
# SCRIPT:  setup_busco_env.sh
#
# WHAT IT DOES:
#   One-time setup helper that creates the dedicated `busco_env` conda
#   environment containing BUSCO and all its dependencies.
#
# WHY THIS IS A SEPARATE SCRIPT (NOT INSIDE run_busco.sbatch):
#   On HPC clusters, compute nodes (where SLURM jobs run) usually have NO
#   internet access. Conda needs internet to download packages, so the
#   environment MUST be created on a login node BEFORE you sbatch the job.
#
# HOW TO USE:
#   Run this ONCE on a login node (e.g., owl1):
#       bash scripts/04_busco/setup_busco_env.sh
#
#   After it finishes successfully, you can submit BUSCO jobs:
#       sbatch scripts/04_busco/run_busco.sbatch MF
#
#   This script is IDEMPOTENT — safe to re-run. If busco_env already exists
#   and is healthy, it just verifies and exits.
#
# WHAT IT INSTALLS:
#   - Python 3.11      (BUSCO doesn't support Python 3.14 yet)
#   - busco            (the main tool)
#   - augustus         (gene prediction)
#   - sepp             (auto-lineage placement)
#   - hmmer            (profile HMM searches)
#   - blast            (sequence searches)
#
# WHERE BUSCO_ENV LIVES:
#   ~/miniconda3/envs/busco_env/   (or wherever your conda is installed)
# =============================================================================

# Stop the script immediately if any command fails
set -euo pipefail


# =============================================================================
# CONFIG — name of the env and the tools to install
# =============================================================================
# You can override the env name from the command line:
#   ENV_NAME=my_busco_env bash setup_busco_env.sh
ENV_NAME="${ENV_NAME:-busco_env}"

# Python version — 3.11 has the broadest bioconda compatibility right now.
# Don't use 3.14 — BUSCO/SEPP haven't been packaged for it.
PYTHON_VERSION="3.11"

# Packages to install. Listed one per line for readability.
PACKAGES=(
    "python=${PYTHON_VERSION}"
    "busco"
    "augustus"
    "sepp"
    "hmmer"
    "blast"
)


# =============================================================================
# HEADER
# =============================================================================
echo "========================================================"
echo "BUSCO Environment Setup"
echo "========================================================"
echo "Env name:       ${ENV_NAME}"
echo "Python version: ${PYTHON_VERSION}"
echo "Packages:       ${PACKAGES[*]}"
echo "Host:           $(hostname)"
echo "Started:        $(date)"
echo "========================================================"


# =============================================================================
# STEP 1 — VERIFY WE'RE ON A LOGIN NODE WITH INTERNET ACCESS
# =============================================================================
echo ""
echo "[1/5] Checking internet access..."

# Quick test — try to reach the conda channel.
# We use --max-time so this fails fast if there's no internet.
if curl -s --max-time 10 https://conda.anaconda.org/ > /dev/null; then
    echo "  ✓ Internet looks OK."
else
    echo "  ✗ No internet access detected."
    echo ""
    echo "  This script must run on a LOGIN NODE (e.g., owl1), NOT a compute node."
    echo "  Compute nodes don't have internet, so conda can't download packages."
    echo "  SSH to a login node, then re-run this script."
    exit 1
fi


# =============================================================================
# STEP 2 — MAKE SURE CONDA IS AVAILABLE
# =============================================================================
echo ""
echo "[2/5] Locating conda..."

# Source conda's shell functions so 'conda activate' works inside this script.
source ~/.bashrc 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh

if ! command -v conda &> /dev/null; then
    echo "  ✗ conda command not found."
    echo ""
    echo "  Install miniconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo "  ✓ conda found: $(conda --version)"

# Prefer mamba if available — it's 5-10x faster at solving env conflicts.
# This is optional; conda alone will work, just slower.
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
    echo "  ✓ Using mamba (faster solver): $(mamba --version | head -1)"
else
    INSTALLER="conda"
    echo "  ~ mamba not found — using conda (slower). Install mamba to speed this up:"
    echo "      conda install -n base -c conda-forge mamba -y"
fi


# =============================================================================
# STEP 3 — CHECK IF busco_env ALREADY EXISTS
# =============================================================================
echo ""
echo "[3/5] Checking if env '${ENV_NAME}' already exists..."

# conda env list prints each env on its own line. We grep for the exact env name
# at the start of a line (^${ENV_NAME}) to avoid false matches like "busco_env_old".
if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    ENV_EXISTS=1
    echo "  ✓ Env '${ENV_NAME}' already exists. Skipping creation."
    echo "    (To rebuild from scratch: conda env remove -n ${ENV_NAME} -y)"
else
    ENV_EXISTS=0
    echo "  ~ Env '${ENV_NAME}' not found — will create it."
fi


# =============================================================================
# STEP 4 — CREATE THE ENV (if it doesn't exist yet)
# =============================================================================
if [[ "${ENV_EXISTS}" -eq 0 ]]; then
    echo ""
    echo "[4/5] Creating env '${ENV_NAME}'..."
    echo "      This takes 5-10 minutes. Conda will download ~1 GB of packages."
    echo ""

    # The actual create command. The "${PACKAGES[@]}" expansion passes each
    # package as a separate argument (proper bash array expansion).
    "${INSTALLER}" create -n "${ENV_NAME}" \
        -c conda-forge \
        -c bioconda \
        "${PACKAGES[@]}" \
        -y

    echo ""
    echo "  ✓ Env '${ENV_NAME}' created."
else
    echo ""
    echo "[4/5] Skipping creation (env already exists)."
fi


# =============================================================================
# STEP 5 — VERIFY THE ENV WORKS
# =============================================================================
echo ""
echo "[5/5] Verifying env contents..."

conda activate "${ENV_NAME}"

# Track any tool that's missing so we can report all problems at once.
MISSING=0

# Why hardcode each tool's version check?
#   Each bioinformatics tool uses a different "print my version" command:
#     python3   → python3 --version           (line 1 has version)
#     busco     → busco --version             (line 1 has version)
#     hmmsearch → hmmsearch -h    (NO --version!  Version is on line 2)
#     tblastn   → tblastn -version  (SINGLE dash, not double!)
#     augustus  → augustus --version          (line 1 has version)
#   Hardcoding each one keeps the output clean and easy to debug if a
#   tool ever changes its version output format.

# ---- Python ----
if command -v python3 &> /dev/null; then
    echo "  ✓ Python: $(python3 --version 2>&1)"
else
    echo "  ✗ Python NOT found in ${ENV_NAME}."
    MISSING=1
fi

# ---- BUSCO ----
if command -v busco &> /dev/null; then
    echo "  ✓ BUSCO: $(busco --version 2>&1 | head -1)"
else
    echo "  ✗ BUSCO NOT found in ${ENV_NAME}."
    MISSING=1
fi

# ---- HMMER ----
# hmmsearch prints its version in the help banner on line 2:
#   line 1: # hmmsearch :: search profile(s) against a sequence database
#   line 2: # HMMER 3.4 (Aug 2023); http://hmmer.org/   ← this one
# We strip the leading "# " for a cleaner display.
if command -v hmmsearch &> /dev/null; then
    HMMER_VERSION=$(hmmsearch -h 2>&1 | grep -m1 "^# HMMER" | sed 's/^# //')
    echo "  ✓ HMMER: ${HMMER_VERSION:-installed (version line not found)}"
else
    echo "  ✗ HMMER NOT found in ${ENV_NAME}."
    MISSING=1
fi

# ---- BLAST+ ----
# tblastn uses SINGLE-dash -version. Line 1 looks like: "tblastn: 2.17.0+"
if command -v tblastn &> /dev/null; then
    echo "  ✓ BLAST+: $(tblastn -version 2>&1 | head -1)"
else
    echo "  ✗ BLAST+ NOT found in ${ENV_NAME}."
    MISSING=1
fi

# ---- AUGUSTUS ----
if command -v augustus &> /dev/null; then
    echo "  ✓ AUGUSTUS: $(augustus --version 2>&1 | head -1)"
else
    echo "  ✗ AUGUSTUS NOT found in ${ENV_NAME}."
    MISSING=1
fi

# SEPP is special — it's a Python module, not a command. Test by importing.
if python3 -c "import sepp" &> /dev/null 2>&1 || command -v run_sepp.py &> /dev/null; then
    echo "  ✓ SEPP found"
else
    echo "  ✗ SEPP NOT found in ${ENV_NAME}."
    MISSING=1
fi

if [[ "${MISSING}" -eq 1 ]]; then
    echo ""
    echo "ERROR: Some tools are missing from '${ENV_NAME}'."
    echo "Try rebuilding the env from scratch:"
    echo "  conda env remove -n ${ENV_NAME} -y"
    echo "  bash scripts/04_busco/setup_busco_env.sh"
    exit 1
fi


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "========================================================"
echo "✓ COMPLETE — '${ENV_NAME}' is ready"
echo "========================================================"
echo "You can now submit BUSCO jobs:"
echo "  sbatch scripts/04_busco/run_busco.sbatch MF      # nutmeg"
echo "  sbatch scripts/04_busco/run_busco.sbatch DC      # carrot"
echo ""
echo "Finished: $(date)"
echo "========================================================"
