#!/bin/bash
# ---------------------------------------------------------------------
# Automated conda environment setup for RNA-seq pipeline
# Source this script in your sbatch files to ensure all dependencies
# are properly installed
# ---------------------------------------------------------------------

setup_rnaseq_env() {
    echo "Setting up RNA-seq conda environment..."
    
    # Activate conda
    if ! bash -c "source ~/.bashrc" 2>/dev/null; then
        echo "Using alternative conda activation"
        source ~/miniconda3/etc/profile.d/conda.sh
    else
        source ~/.bashrc
    fi
    
    # Activate environment
    conda activate rnaseq || {
        echo "ERROR: Failed to activate rnaseq environment"
        echo "Create it with: conda env create -f environment.yml"
        return 1
    }
    
    echo "✓ Conda environment activated: rnaseq"
    
    # Check critical dependencies and fix if needed
    local needs_fix=0
    
    # Check bowtie2-build-s (needed for RSEM)
    if ! command -v bowtie2-build-s &> /dev/null; then
        echo "⚠ bowtie2-build-s not found - will reinstall bowtie2"
        needs_fix=1
    fi
    
    # Check samtools version
    if command -v samtools &> /dev/null; then
        samtools_version=$(samtools --version 2>&1 | head -1 | awk '{print $2}')
        samtools_major=$(echo "$samtools_version" | cut -d. -f1)
        samtools_minor=$(echo "$samtools_version" | cut -d. -f2)
        if [[ $samtools_major -lt 1 ]] || [[ $samtools_major -eq 1 && $samtools_minor -lt 3 ]]; then
            echo "⚠ samtools version $samtools_version is too old (need >= 1.3)"
            needs_fix=1
        fi
    else
        echo "⚠ samtools not found"
        needs_fix=1
    fi
    
    # Auto-fix if needed
    if [[ $needs_fix -eq 1 ]]; then
        echo ""
        echo "Auto-fixing dependencies (this may take a few minutes)..."
        echo "This only needs to happen once."
        
        # Reinstall bowtie2 with all components
        conda install -y -c bioconda bowtie2 2>&1 | tail -5
        
        # Ensure samtools is recent
        conda install -y -c bioconda 'samtools>=1.3' 2>&1 | tail -5
        
        # Ensure RSEM is installed
        conda install -y -c bioconda rsem 2>&1 | tail -5
        
        echo "✓ Dependencies fixed!"
    fi
    
    # Final verification
    echo ""
    echo "Environment check:"
    echo "  Trinity:         $(command -v Trinity &> /dev/null && echo '✓' || echo '✗ MISSING')"
    echo "  bowtie2:         $(command -v bowtie2 &> /dev/null && echo '✓' || echo '✗ MISSING')"
    echo "  bowtie2-build-s: $(command -v bowtie2-build-s &> /dev/null && echo '✓' || echo '✗ MISSING')"
    echo "  samtools:        $(command -v samtools &> /dev/null && samtools --version | head -1 || echo '✗ MISSING')"
    echo "  RSEM:            $(command -v rsem-calculate-expression &> /dev/null && echo '✓' || echo '✗ MISSING')"
    echo ""
    
    return 0
}

# Run setup if script is executed (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    setup_rnaseq_env
fi
