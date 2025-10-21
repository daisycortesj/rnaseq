#!/bin/bash
# Fix Trinity path issues after conda installation
# This script creates the necessary symlinks for Trinity to work properly

set -e

echo "Fixing Trinity installation paths..."

# Get the conda environment path
CONDA_PREFIX=${CONDA_PREFIX:-$CONDA_ENV_PATH}

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: CONDA_PREFIX not set. Please activate the conda environment first."
    exit 1
fi

echo "Conda environment: $CONDA_PREFIX"

# Check if Trinity is installed
if [ ! -f "$CONDA_PREFIX/bin/Trinity" ]; then
    echo "Error: Trinity not found in $CONDA_PREFIX/bin/Trinity"
    exit 1
fi

# Create the expected directory structure
echo "Creating directory structure..."
mkdir -p "$CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts"

# Create symlinks to fix path issues
echo "Creating symlinks..."

# Main Trinity executable
if [ ! -L "$CONDA_PREFIX/opt/trinity-2.9.1/Trinity" ]; then
    ln -sf "$CONDA_PREFIX/bin/Trinity" "$CONDA_PREFIX/opt/trinity-2.9.1/Trinity"
    echo "Created symlink: $CONDA_PREFIX/opt/trinity-2.9.1/Trinity"
else
    echo "Symlink already exists: $CONDA_PREFIX/opt/trinity-2.9.1/Trinity"
fi

# Trinity executable for support scripts
if [ ! -L "$CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts/../../Trinity" ]; then
    ln -sf "$CONDA_PREFIX/bin/Trinity" "$CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts/../../Trinity"
    echo "Created symlink: $CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts/../../Trinity"
else
    echo "Symlink already exists: $CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts/../../Trinity"
fi

# Test Trinity
echo "Testing Trinity installation..."
if "$CONDA_PREFIX/bin/Trinity" --version > /dev/null 2>&1; then
    echo "✓ Trinity is working correctly"
    "$CONDA_PREFIX/bin/Trinity" --version
else
    echo "✗ Trinity test failed"
    exit 1
fi

echo ""
echo "Trinity path fix completed successfully!"
echo "You can now run Trinity commands without path errors."