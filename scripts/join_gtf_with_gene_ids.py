#!/usr/bin/env python3
"""
Join parsed GTF data with gene IDs list (inner join).

This script takes the output from parse_refseq_gtf.py and filters it
to keep only genes that appear in your all_gene_ids.txt file.

USAGE:
    python join_gtf_with_gene_ids.py --parsed genes_parsed.tsv --gene-ids all_gene_ids.txt --out filtered_genes.tsv

Author: Daisy Cortes
"""

import argparse
import sys


def load_gene_ids(gene_ids_file):
    """
    Load gene IDs from a text file (one gene ID per line).
    
    Returns a set of gene IDs for fast lookup.
    """
    gene_ids = set()
    
    print(f"Loading gene IDs from: {gene_ids_file}")
    
    with open(gene_ids_file, 'r') as f:
        for line in f:
            # Strip whitespace and skip empty lines
            gene_id = line.strip()
            if gene_id:
                gene_ids.add(gene_id)
    
    print(f"  Loaded {len(gene_ids)} gene IDs")
    return gene_ids


def inner_join_with_gene_ids(parsed_file, gene_ids_file, output_file):
    """
    Perform inner join between parsed GTF and gene IDs list.
    
    Only keeps rows where gene_id appears in the gene_ids_file.
    """
    # Load gene IDs into a set for fast lookup
    target_gene_ids = load_gene_ids(gene_ids_file)
    
    # Process the parsed GTF file
    print(f"\nProcessing parsed GTF: {parsed_file}")
    
    matched_rows = []
    total_rows = 0
    
    with open(parsed_file, 'r') as f:
        # Read and save header
        header = f.readline().strip()
        
        # Process each data row
        for line in f:
            total_rows += 1
            fields = line.strip().split('\t')
            
            # Assuming format: local_number, gene_id, transcript_id, description
            if len(fields) >= 2:
                gene_id = fields[1]  # Second column is gene_id
                
                # Check if this gene_id is in our target list
                if gene_id in target_gene_ids:
                    matched_rows.append(line.strip())
    
    print(f"  Total rows in parsed file: {total_rows}")
    print(f"  Matched rows (in gene IDs list): {len(matched_rows)}")
    print(f"  Filtered out: {total_rows - len(matched_rows)}")
    
    # Write output
    print(f"\nWriting output to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write(header + '\n')
        
        # Write matched rows
        for row in matched_rows:
            f.write(row + '\n')
    
    print(f"✓ Inner join complete!")
    print(f"  Output file has {len(matched_rows)} genes")
    
    # Check for missing genes
    found_gene_ids = set()
    for row in matched_rows:
        fields = row.split('\t')
        if len(fields) >= 2:
            found_gene_ids.add(fields[1])
    
    missing_genes = target_gene_ids - found_gene_ids
    if missing_genes:
        print(f"\n⚠️  Warning: {len(missing_genes)} genes from gene IDs list were NOT found in GTF")
        if len(missing_genes) <= 10:
            print(f"  Missing genes: {', '.join(sorted(list(missing_genes)))}")
        else:
            print(f"  First 10 missing genes: {', '.join(sorted(list(missing_genes))[:10])}")


def main():
    parser = argparse.ArgumentParser(
        description='Inner join parsed GTF with gene IDs list',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python join_gtf_with_gene_ids.py --parsed genes_parsed.tsv --gene-ids all_gene_ids.txt --out filtered_genes.tsv
  
  # With full paths
  python join_gtf_with_gene_ids.py \\
    --parsed 06_analysis/genes_parsed.tsv \\
    --gene-ids 06_analysis/pydeseq2_DC_step1_unfiltered/all_gene_ids.txt \\
    --out 06_analysis/filtered_genes_with_annotations.tsv
        """
    )
    
    parser.add_argument(
        '--parsed',
        required=True,
        help='Parsed GTF file (output from parse_refseq_gtf.py)'
    )
    
    parser.add_argument(
        '--gene-ids',
        required=True,
        help='Gene IDs file (one gene ID per line, e.g., all_gene_ids.txt)'
    )
    
    parser.add_argument(
        '--out',
        required=True,
        help='Output file (filtered genes with annotations)'
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    try:
        with open(args.parsed, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Parsed GTF file not found: {args.parsed}", file=sys.stderr)
        sys.exit(1)
    
    try:
        with open(args.gene_ids, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Gene IDs file not found: {args.gene_ids}", file=sys.stderr)
        sys.exit(1)
    
    # Perform inner join
    try:
        inner_join_with_gene_ids(args.parsed, args.gene_ids, args.out)
    except Exception as e:
        print(f"Error during inner join: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
