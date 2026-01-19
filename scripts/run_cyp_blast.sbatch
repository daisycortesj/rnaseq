#!/bin/bash
#SBATCH --job-name=cyp_blast
#SBATCH --output=cyp_blast_%j.out
#SBATCH --error=cyp_blast_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=normal_q
#SBATCH --account=tholl_lab_1

################################################################################
# CYP Gene Identification via BLAST
#
# This script identifies Cytochrome P450 (CYP) genes in your genome by:
# 1. Extracting protein sequences from your annotation
# 2. BLASTing against a plant P450 reference database
# 3. Creating a comprehensive CYP family mapping file
#
# Usage:
#   # Basic (uses defaults)
#   sbatch scripts/run_cyp_blast.sbatch
#
#   # With custom paths
#   PROTEIN_FASTA=/path/to/proteins.faa \
#   GTF_FILE=/path/to/annotation.gtf \
#   OUTPUT_DIR=/path/to/output \
#   sbatch scripts/run_cyp_blast.sbatch
#
# Output:
#   - cyp_families_blast.tsv  (gene_id, cyp_family, blast_evalue, identity)
#   - cyp_blast_results.txt   (raw BLAST output)
#   - cyp_genes_list.txt      (simple list of CYP gene IDs)
#
################################################################################

set -euo pipefail

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║              CYP Gene Identification via BLAST                     ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-8}  MEM: ${SLURM_MEM_PER_NODE:-16G}"
echo "Start time: $(date)"
echo ""

# =============================================================================
# Configuration
# =============================================================================

# Input files
PROTEIN_FASTA="${PROTEIN_FASTA:-/projects/tholl_lab_1/daisy_analysis/04_reference/ncbi_dataset/data/GCF_001625215.2/protein.faa}"
GTF_FILE="${GTF_FILE:-/projects/tholl_lab_1/daisy_analysis/04_reference/dc_genomic.gtf}"
GENOME_FASTA="${GENOME_FASTA:-/projects/tholl_lab_1/daisy_analysis/04_reference/GCF_001625215.2_DH1_v3.0_genomic.fna}"

# Output directory
OUTPUT_DIR="${OUTPUT_DIR:-/projects/tholl_lab_1/daisy_analysis/07_cyp_analysis}"

# BLAST parameters
EVALUE="${EVALUE:-1e-20}"
IDENTITY_CUTOFF="${IDENTITY_CUTOFF:-40}"  # 40% identity = same CYP family
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Reference database (will be downloaded if not present)
P450_REF_DIR="${OUTPUT_DIR}/p450_reference"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "CONFIGURATION"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Protein FASTA: ${PROTEIN_FASTA}"
echo "  GTF file:      ${GTF_FILE}"
echo "  Output dir:    ${OUTPUT_DIR}"
echo "  E-value:       ${EVALUE}"
echo "  Identity:      ${IDENTITY_CUTOFF}%"
echo "  Threads:       ${THREADS}"
echo ""

# =============================================================================
# Setup
# =============================================================================

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${P450_REF_DIR}"
cd "${OUTPUT_DIR}"

# Activate conda environment
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ACTIVATING ENVIRONMENT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate rnaseq 2>/dev/null || conda activate base
    echo "✓ Conda environment activated"
fi

# Check for BLAST
if ! command -v blastp &> /dev/null; then
    echo "Loading BLAST module..."
    module load BLAST/2.13.0-gompi-2022a 2>/dev/null || \
    module load blast 2>/dev/null || \
    module load ncbi-blast 2>/dev/null || \
    { echo "ERROR: BLAST not found. Please install or load BLAST module."; exit 1; }
fi
echo "✓ BLAST available: $(which blastp)"

# =============================================================================
# Step 1: Validate input files
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Validating input files"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check for protein FASTA
if [[ ! -f "${PROTEIN_FASTA}" ]]; then
    echo "WARNING: Protein FASTA not found at ${PROTEIN_FASTA}"
    echo "Searching for protein files..."
    
    # Try to find protein file
    FOUND_PROTEIN=$(find /projects/tholl_lab_1/daisy_analysis/04_reference -name "*.faa" -o -name "*protein*.fa*" 2>/dev/null | head -1)
    
    if [[ -n "${FOUND_PROTEIN}" ]]; then
        PROTEIN_FASTA="${FOUND_PROTEIN}"
        echo "Found: ${PROTEIN_FASTA}"
    else
        echo "ERROR: No protein FASTA found."
        echo ""
        echo "You need to download protein sequences from NCBI:"
        echo "  1. Go to: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001625215.2/"
        echo "  2. Download 'Protein FASTA (.faa)'"
        echo "  3. Set PROTEIN_FASTA=/path/to/protein.faa"
        exit 1
    fi
fi
echo "✓ Protein FASTA: ${PROTEIN_FASTA}"
echo "  $(grep -c "^>" "${PROTEIN_FASTA}") protein sequences"

# Check GTF
if [[ -f "${GTF_FILE}" ]]; then
    echo "✓ GTF file: ${GTF_FILE}"
else
    echo "WARNING: GTF not found (optional for BLAST)"
fi

# =============================================================================
# Step 2: Download/prepare P450 reference database
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2: Preparing P450 reference database"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

P450_REF="${P450_REF_DIR}/plant_p450_reference.fasta"

if [[ ! -f "${P450_REF}" ]]; then
    echo "Downloading plant P450 reference sequences..."
    
    # Create a curated plant P450 reference from UniProt
    # This includes well-annotated plant P450s
    cat > "${P450_REF_DIR}/download_p450.py" << 'PYEOF'
#!/usr/bin/env python3
"""Download plant P450 reference sequences from UniProt"""
import urllib.request
import sys

# UniProt query for reviewed plant P450s
# taxonomy:3193 = Embryophyta (land plants)
# family:cytochrome P450
query = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A3193%29+AND+%28family%3A%22cytochrome+p450%22%29+AND+%28reviewed%3Atrue%29"

print("Downloading plant P450 sequences from UniProt...")
try:
    with urllib.request.urlopen(query, timeout=300) as response:
        data = response.read().decode('utf-8')
    
    output_file = sys.argv[1] if len(sys.argv) > 1 else "plant_p450_reference.fasta"
    with open(output_file, 'w') as f:
        f.write(data)
    
    # Count sequences
    n_seqs = data.count('>')
    print(f"Downloaded {n_seqs} P450 sequences to {output_file}")
    
except Exception as e:
    print(f"Error downloading: {e}")
    print("Creating minimal reference from known P450 families...")
    
    # Fallback: create a minimal reference with key P450 motifs
    # This is the conserved heme-binding region
    minimal_ref = """>CYP71_consensus cytochrome P450 71 family consensus
FGAGRRICAG
>CYP72_consensus cytochrome P450 72 family consensus  
FGAGRRICAG
>CYP82_consensus cytochrome P450 82 family consensus
FGAGRRICAG
>CYP86_consensus cytochrome P450 86 family consensus
FGAGRRICAG
"""
    with open(sys.argv[1] if len(sys.argv) > 1 else "plant_p450_reference.fasta", 'w') as f:
        f.write(minimal_ref)
PYEOF
    
    python3 "${P450_REF_DIR}/download_p450.py" "${P450_REF}"
    
    if [[ ! -s "${P450_REF}" ]]; then
        echo "WARNING: Could not download P450 reference."
        echo "Using GTF-based extraction as fallback."
    fi
else
    echo "✓ P450 reference exists: ${P450_REF}"
fi

N_REF=$(grep -c "^>" "${P450_REF}" 2>/dev/null || echo "0")
echo "  Reference contains ${N_REF} P450 sequences"

# =============================================================================
# Step 3: Create BLAST database from carrot proteins
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 3: Creating BLAST database"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

BLAST_DB="${OUTPUT_DIR}/carrot_proteins"

if [[ ! -f "${BLAST_DB}.phr" ]]; then
    echo "Building BLAST database from carrot proteins..."
    makeblastdb -in "${PROTEIN_FASTA}" \
                -dbtype prot \
                -out "${BLAST_DB}" \
                -title "Daucus_carota_proteins"
    echo "✓ BLAST database created"
else
    echo "✓ BLAST database exists"
fi

# =============================================================================
# Step 4: Run BLAST search
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 4: Running BLAST search"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

BLAST_OUTPUT="${OUTPUT_DIR}/cyp_blast_results.txt"

if [[ -s "${P450_REF}" ]] && [[ $(grep -c "^>" "${P450_REF}") -gt 10 ]]; then
    echo "Running BLAST: P450 references vs carrot proteins..."
    echo "  This may take 10-30 minutes..."
    
    blastp -query "${P450_REF}" \
           -db "${BLAST_DB}" \
           -out "${BLAST_OUTPUT}" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
           -evalue "${EVALUE}" \
           -num_threads "${THREADS}" \
           -max_target_seqs 5
    
    echo "✓ BLAST completed"
    echo "  $(wc -l < "${BLAST_OUTPUT}") hits found"
else
    echo "Skipping BLAST (no reference database)"
    touch "${BLAST_OUTPUT}"
fi

# =============================================================================
# Step 5: Also search GTF for CYP annotations (complementary approach)
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 5: Extracting CYP genes from GTF annotations"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

GTF_CYP="${OUTPUT_DIR}/cyp_from_gtf.txt"

if [[ -f "${GTF_FILE}" ]]; then
    # Broader search pattern to catch more CYP genes
    grep -i "cytochrome P450\|CYP[0-9]\|monooxygenase.*CYP\|hydroxylase.*CYP" "${GTF_FILE}" | \
        grep -w "gene" | \
        grep -v -i "cyclophilin\|peptidyl-prolyl" | \
        grep -oP 'gene_id "\K[^"]+' | \
        sort -u > "${GTF_CYP}"
    
    echo "✓ Found $(wc -l < "${GTF_CYP}") CYP genes from GTF annotations"
else
    touch "${GTF_CYP}"
    echo "  GTF not available, skipping"
fi

# =============================================================================
# Step 6: Parse results and create CYP family mapping
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 6: Creating CYP family mapping"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Python script to parse BLAST results and create family mapping
cat > "${OUTPUT_DIR}/parse_blast_results.py" << 'PYEOF'
#!/usr/bin/env python3
"""
Parse BLAST results and GTF annotations to create comprehensive CYP family mapping.
"""
import re
import sys
from collections import defaultdict

def extract_cyp_family(description):
    """Extract CYP family number from a description string."""
    # Look for CYP followed by numbers
    match = re.search(r'CYP(\d+)', description, re.IGNORECASE)
    if match:
        return f"CYP{match.group(1)}"
    return None

def parse_blast_results(blast_file, identity_cutoff=40):
    """Parse BLAST output and extract CYP gene hits."""
    cyp_genes = {}  # gene_id -> (family, evalue, identity)
    
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                
                query_id = parts[0]  # P450 reference
                subject_id = parts[1]  # Carrot protein
                identity = float(parts[2])
                evalue = float(parts[10])
                
                if identity < identity_cutoff:
                    continue
                
                # Extract CYP family from query description
                family = extract_cyp_family(query_id)
                if not family:
                    family = "CYP_unknown"
                
                # Extract gene ID from subject (protein ID -> gene ID)
                # NCBI protein IDs are like XP_017229978.1
                gene_id = subject_id.split('.')[0]
                
                # Keep best hit per gene
                if gene_id not in cyp_genes or evalue < cyp_genes[gene_id][1]:
                    cyp_genes[gene_id] = (family, evalue, identity)
    except FileNotFoundError:
        pass
    
    return cyp_genes

def parse_gtf_cyp(gtf_cyp_file, gtf_file=None):
    """Parse GTF-extracted CYP genes and get their families."""
    cyp_genes = {}
    
    # Read gene IDs from GTF extraction
    gene_ids = set()
    try:
        with open(gtf_cyp_file, 'r') as f:
            for line in f:
                gene_ids.add(line.strip())
    except FileNotFoundError:
        return cyp_genes
    
    # If GTF provided, extract family info
    if gtf_file:
        try:
            with open(gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    # Check if this line contains a gene we care about
                    gene_match = re.search(r'gene_id "([^"]+)"', line)
                    if not gene_match:
                        continue
                    
                    gene_id = gene_match.group(1)
                    if gene_id not in gene_ids:
                        continue
                    
                    # Extract description
                    desc_match = re.search(r'description "([^"]+)"', line)
                    if desc_match:
                        description = desc_match.group(1)
                        family = extract_cyp_family(description)
                        if family and gene_id not in cyp_genes:
                            cyp_genes[gene_id] = (family, 0, 100)  # Annotation = perfect match
        except FileNotFoundError:
            pass
    
    # For genes without family info, assign "CYP_unknown"
    for gene_id in gene_ids:
        if gene_id not in cyp_genes:
            cyp_genes[gene_id] = ("CYP_unknown", 0, 100)
    
    return cyp_genes

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', default='cyp_blast_results.txt')
    parser.add_argument('--gtf-cyp', default='cyp_from_gtf.txt')
    parser.add_argument('--gtf', default=None)
    parser.add_argument('--output', default='cyp_families_blast.tsv')
    parser.add_argument('--identity', type=float, default=40)
    args = parser.parse_args()
    
    # Combine results from both methods
    all_cyp = {}
    
    # BLAST results (higher confidence)
    blast_cyp = parse_blast_results(args.blast, args.identity)
    all_cyp.update(blast_cyp)
    print(f"CYP genes from BLAST: {len(blast_cyp)}")
    
    # GTF annotations (fills gaps)
    gtf_cyp = parse_gtf_cyp(args.gtf_cyp, args.gtf)
    for gene_id, info in gtf_cyp.items():
        if gene_id not in all_cyp:
            all_cyp[gene_id] = info
    print(f"CYP genes from GTF: {len(gtf_cyp)}")
    print(f"Total unique CYP genes: {len(all_cyp)}")
    
    # Write output
    with open(args.output, 'w') as f:
        f.write("gene_id\tcyp_family\tevalue\tidentity\n")
        for gene_id in sorted(all_cyp.keys()):
            family, evalue, identity = all_cyp[gene_id]
            f.write(f"{gene_id}\t{family}\t{evalue}\t{identity}\n")
    
    print(f"Output written to: {args.output}")
    
    # Summary by family
    family_counts = defaultdict(int)
    for gene_id, (family, _, _) in all_cyp.items():
        family_counts[family] += 1
    
    print("\nCYP families found:")
    for family in sorted(family_counts.keys(), key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0):
        print(f"  {family}: {family_counts[family]} genes")

if __name__ == '__main__':
    main()
PYEOF

python3 "${OUTPUT_DIR}/parse_blast_results.py" \
    --blast "${BLAST_OUTPUT}" \
    --gtf-cyp "${GTF_CYP}" \
    --gtf "${GTF_FILE}" \
    --output "${OUTPUT_DIR}/cyp_families_blast.tsv" \
    --identity "${IDENTITY_CUTOFF}"

# =============================================================================
# Step 7: Create simple gene list for pydeseq2_analysis.py
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 7: Creating output files"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Simple CYP gene list
tail -n +2 "${OUTPUT_DIR}/cyp_families_blast.tsv" | cut -f1 > "${OUTPUT_DIR}/cyp_genes_list.txt"

# Create compatible format for pydeseq2_analysis.py (just gene_id and cyp_family)
echo -e "gene_id\tcyp_family" > "${OUTPUT_DIR}/cyp_families_comprehensive.tsv"
tail -n +2 "${OUTPUT_DIR}/cyp_families_blast.tsv" | cut -f1,2 >> "${OUTPUT_DIR}/cyp_families_comprehensive.tsv"

echo "✓ Created output files:"
echo "  - cyp_families_blast.tsv      (full results with e-value/identity)"
echo "  - cyp_families_comprehensive.tsv (for pydeseq2_analysis.py)"
echo "  - cyp_genes_list.txt          (simple gene list)"

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                         SUMMARY                                    ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Total CYP genes identified: $(tail -n +2 "${OUTPUT_DIR}/cyp_families_comprehensive.tsv" | wc -l)"
echo ""
echo "Output files:"
echo "  ${OUTPUT_DIR}/cyp_families_comprehensive.tsv"
echo ""
echo "To use with PyDESeq2 heatmap:"
echo ""
echo "  CYP_FAMILY_MAP=${OUTPUT_DIR}/cyp_families_comprehensive.tsv \\"
echo "  CONTRAST_A=R \\"
echo "  CONTRAST_B=L \\"
echo "  ROOT_UP_ONLY=true \\"
echo "  sbatch scripts/run_pydeseq2_analysis.sbatch"
echo ""
echo "End time: $(date)"
echo "Done!"
