#!/usr/bin/env python3
"""
Extract CYP (Cytochrome P450) Gene Families from GTF

This script automatically extracts CYP gene information from a GTF annotation file
and creates the cyp_families.tsv file needed for the CYP heatmap in pydeseq2_analysis.py.

Usage:
    python extract_cyp_families.py <gtf_file> -o <output_file>

Example:
    python extract_cyp_families.py reference/genomic.gtf -o cyp_families.tsv

The script:
1. Reads the GTF file
2. Finds all genes annotated as "cytochrome P450" or containing CYP family names
3. Extracts the CYP family (e.g., CYP71, CYP72, CYP86)
4. Excludes cyclophilins (peptidyl-prolyl isomerases that have "CYP" in name but are NOT P450s)
5. Outputs a TSV with gene_id and cyp_family columns

Output format:
    gene_id     cyp_family
    LOC108193079    CYP71
    LOC108196352    CYP86
    LOC108201942    CYP707
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def extract_cyp_from_gtf(gtf_file, output_file, verbose=True):
    """
    Extract CYP gene families from a GTF file.
    
    Parameters:
        gtf_file: Path to GTF annotation file
        output_file: Path to output TSV file
        verbose: Print progress messages
        
    Returns:
        Dictionary mapping gene_id to cyp_family
    """
    if verbose:
        print(f"Reading GTF file: {gtf_file}")
    
    # Patterns to identify CYP genes
    # Match "cytochrome P450" in description OR CYP followed by numbers (family)
    cyp_pattern = re.compile(r'CYP(\d+)', re.IGNORECASE)
    
    # Words that indicate this is NOT a cytochrome P450 (cyclophilins have CYP in name)
    exclude_patterns = [
        'peptidyl-prolyl',
        'isomerase',
        'cyclophilin',
        'PPIase'
    ]
    
    cyp_genes = {}  # gene_id -> (cyp_family, full_description)
    genes_seen = set()
    lines_processed = 0
    
    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            lines_processed += 1
            
            # Only process gene entries (not transcripts, exons, etc.)
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            attributes = fields[8]
            
            # Extract gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)
            
            # Skip if we've already processed this gene
            if gene_id in genes_seen:
                continue
            genes_seen.add(gene_id)
            
            # Extract description
            desc_match = re.search(r'description "([^"]+)"', attributes)
            if not desc_match:
                continue
            description = desc_match.group(1)
            
            # Check if this is a cytochrome P450 (not a cyclophilin)
            description_lower = description.lower()
            
            # Skip if it matches exclusion patterns (cyclophilins)
            if any(excl in description_lower for excl in exclude_patterns):
                continue
            
            # Must contain "cytochrome P450" or "cytochrome p450" to be a real P450
            if 'cytochrome p450' not in description_lower:
                # Also accept if it just has CYP### pattern without "cytochrome P450"
                # but only if it doesn't match exclusions
                if not cyp_pattern.search(description):
                    continue
            
            # Extract CYP family number
            cyp_match = cyp_pattern.search(description)
            if cyp_match:
                cyp_family = f"CYP{cyp_match.group(1)}"
                cyp_genes[gene_id] = (cyp_family, description)
    
    if verbose:
        print(f"  Processed {lines_processed:,} lines")
        print(f"  Found {len(genes_seen):,} unique genes")
        print(f"  Found {len(cyp_genes):,} CYP (cytochrome P450) genes")
    
    # Count genes per family
    family_counts = defaultdict(int)
    for gene_id, (family, desc) in cyp_genes.items():
        family_counts[family] += 1
    
    if verbose and cyp_genes:
        print(f"\n  CYP families found:")
        for family in sorted(family_counts.keys(), key=lambda x: int(re.search(r'\d+', x).group())):
            print(f"    {family}: {family_counts[family]} genes")
    
    # Write output file
    if verbose:
        print(f"\nWriting output: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("gene_id\tcyp_family\n")
        for gene_id in sorted(cyp_genes.keys()):
            cyp_family, _ = cyp_genes[gene_id]
            f.write(f"{gene_id}\t{cyp_family}\n")
    
    if verbose:
        print(f"  Wrote {len(cyp_genes)} CYP gene entries")
    
    return cyp_genes


def main():
    parser = argparse.ArgumentParser(
        description="Extract CYP gene families from GTF annotation file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python extract_cyp_families.py genomic.gtf -o cyp_families.tsv
    
    # With full paths
    python extract_cyp_families.py /path/to/reference/genomic.gtf \\
        -o /path/to/output/cyp_families.tsv

Output format (TSV):
    gene_id         cyp_family
    LOC108193079    CYP71
    LOC108196352    CYP86
    LOC108201942    CYP707

Notes:
    - Only extracts true cytochrome P450 genes
    - Excludes cyclophilins (peptidyl-prolyl isomerases) which have CYP in name
    - CYP family is extracted as the base number (e.g., CYP71D312 -> CYP71)
        """
    )
    
    parser.add_argument(
        "gtf_file",
        help="Path to GTF annotation file"
    )
    parser.add_argument(
        "-o", "--output",
        default="cyp_families.tsv",
        help="Output TSV file (default: cyp_families.tsv)"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress messages"
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not Path(args.gtf_file).exists():
        print(f"ERROR: GTF file not found: {args.gtf_file}")
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Run extraction
    try:
        cyp_genes = extract_cyp_from_gtf(
            args.gtf_file, 
            args.output, 
            verbose=not args.quiet
        )
        
        if len(cyp_genes) == 0:
            print("\nWARNING: No CYP genes found!")
            print("Check that your GTF file contains 'cytochrome P450' or 'CYP' in gene descriptions.")
            sys.exit(1)
        
        if not args.quiet:
            print("\n✓ CYP family extraction complete!")
            print(f"  Output: {args.output}")
            print(f"  Total CYP genes: {len(cyp_genes)}")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
