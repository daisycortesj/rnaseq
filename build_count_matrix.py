#!/usr/bin/env python3
"""
Build gene count matrix from STAR ReadsPerGene.out.tab files or Trinity RSEM .genes.results files
and prepare for statistical analysis with DESeq2/edgeR/PyDESeq2.
"""

import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re

def parse_sample_name(filename, file_type='star'):
    """
    Parse sample name from filename.
    
    Args:
        filename: Path to count file
        file_type: 'star' for STAR files, 'rsem' for RSEM files
    
    Returns:
        Sample name string
    """
    basename = os.path.basename(filename)
    if file_type == 'star':
        sample_name = basename.replace('_ReadsPerGene.out.tab', '')
    elif file_type == 'rsem':
        # RSEM files: sample.genes.results or sample/rsem.genes.results
        sample_name = basename.replace('.genes.results', '')
        # If it's a directory path, extract from parent directory
        if '/' in str(filename):
            parent_dir = os.path.basename(os.path.dirname(filename))
            # Try to use parent directory name if it looks like a sample name
            if not sample_name or sample_name == 'rsem':
                sample_name = parent_dir
    else:
        sample_name = basename
    return sample_name

def extract_sample_info(sample_name):
    """
    Extract sample information from sample name.
    Expected format: DC1L1, DC1R1, DGL1, DGR1, etc.
    """
    # Pattern: (Group)(Number)(Condition)(Replicate)
    pattern = r'([A-Z]+)(\d*)([LR])(\d+)'
    match = re.match(pattern, sample_name)
    
    if match:
        group = match.group(1)  # DC, DG, etc.
        group_num = match.group(2) if match.group(2) else ""  # 1, 2, etc.
        condition = match.group(3)  # L, R
        replicate = match.group(4)  # 1, 2, 3
        
        return {
            'sample': sample_name,
            'group': group,
            'group_number': group_num,
            'condition': condition,
            'replicate': replicate,
            'treatment': f"{group}{group_num}",
            'full_condition': f"{group}{group_num}_{condition}"
        }
    else:
        # Fallback for different naming conventions
        return {
            'sample': sample_name,
            'group': sample_name[:2] if len(sample_name) >= 2 else sample_name,
            'group_number': "",
            'condition': "",
            'replicate': "",
            'treatment': sample_name,
            'full_condition': sample_name
        }

def read_star_counts(filepath):
    """
    Read STAR ReadsPerGene.out.tab file and return counts.
    
    STAR output format:
    Column 1: gene_id
    Column 2: unstranded counts
    Column 3: first strand counts  
    Column 4: second strand counts
    
    We'll use unstranded counts (column 2) by default.
    """
    try:
        df = pd.read_csv(filepath, sep='\t', header=None, 
                        names=['gene_id', 'unstranded', 'first_strand', 'second_strand'])
        
        # Remove the first 4 rows which contain summary statistics
        df = df.iloc[4:].reset_index(drop=True)
        
        # Convert counts to integer
        df['unstranded'] = pd.to_numeric(df['unstranded'], errors='coerce').fillna(0).astype(int)
        
        return df[['gene_id', 'unstranded']]
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def read_rsem_counts(filepath):
    """
    Read RSEM .genes.results file and return expected counts.
    
    RSEM output format:
    - gene_id
    - transcript_id(s)
    - length
    - effective_length
    - expected_count (this is what we want for count matrices)
    - TPM
    - FPKM
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        
        # Check if required columns exist
        if 'expected_count' not in df.columns:
            raise ValueError(f"RSEM file {filepath} missing 'expected_count' column")
        if 'gene_id' not in df.columns:
            # Sometimes it's just the first column without a header
            if len(df.columns) >= 5:
                df.columns = ['gene_id', 'transcript_id', 'length', 'effective_length', 
                             'expected_count', 'TPM', 'FPKM'][:len(df.columns)]
            else:
                raise ValueError(f"RSEM file {filepath} has unexpected format")
        
        # Convert expected_count to integer (round first)
        df['expected_count'] = pd.to_numeric(df['expected_count'], errors='coerce').fillna(0)
        df['expected_count'] = df['expected_count'].round().astype(int)
        
        return df[['gene_id', 'expected_count']].rename(columns={'expected_count': 'counts'})
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def build_count_matrix(count_dir, output_dir="count_matrices", count_type="auto"):
    """
    Build gene count matrix from STAR or Trinity/RSEM count files.
    
    Args:
        count_dir: Directory containing count files
        output_dir: Output directory for count matrix and metadata
        count_type: 'star', 'rsem', or 'auto' (auto-detect)
    """
    count_dir = Path(count_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Check if output files already exist - skip if they do
    count_matrix_file = output_dir / "gene_count_matrix.tsv"
    metadata_file = output_dir / "sample_metadata.tsv"
    
    if count_matrix_file.exists() and metadata_file.exists():
        print(f"\nOutput files already exist in {output_dir}/")
        print(f"  - {count_matrix_file.name}")
        print(f"  - {metadata_file.name}")
        print("Skipping matrix creation. Files will not be overwritten.")
        print("\nTo regenerate, delete the existing files or use a different output directory.")
        # Load and return existing files
        existing_matrix = pd.read_csv(count_matrix_file, sep='\t', index_col=0)
        existing_metadata = pd.read_csv(metadata_file, sep='\t')
        return existing_matrix, existing_metadata
    
    # Auto-detect count type if not specified
    if count_type == "auto":
        star_files = list(count_dir.glob("*ReadsPerGene.out.tab"))
        rsem_files = list(count_dir.glob("**/*.genes.results"))
        
        if star_files and not rsem_files:
            count_type = "star"
        elif rsem_files and not star_files:
            count_type = "rsem"
        elif star_files and rsem_files:
            # Both found, prefer STAR
            print("Warning: Both STAR and RSEM files found. Using STAR files.")
            count_type = "star"
        else:
            raise ValueError(f"No count files found in {count_dir}")
    
    # Find count files based on type
    if count_type == "star":
        count_files = list(count_dir.glob("*ReadsPerGene.out.tab"))
        if not count_files:
            raise ValueError(f"No ReadsPerGene.out.tab files found in {count_dir}")
    elif count_type == "rsem":
        # RSEM files (from Trinity workflow) can be in subdirectories
        count_files = list(count_dir.glob("**/*.genes.results"))
        if not count_files:
            raise ValueError(f"No .genes.results files found in {count_dir}")
    else:
        raise ValueError(f"Unknown count_type: {count_type}")
    
    print(f"Found {len(count_files)} {count_type.upper()} count files")
    
    # Read all count files
    count_data = {}
    sample_info = []
    
    for filepath in count_files:
        sample_name = parse_sample_name(str(filepath), file_type=count_type)
        print(f"Processing {sample_name}...")
        
        if count_type == "star":
            counts = read_star_counts(filepath)
            if counts is not None:
                counts = counts.rename(columns={'unstranded': 'counts'})
        else:  # rsem
            counts = read_rsem_counts(filepath)
        
        if counts is not None:
            count_data[sample_name] = counts
            sample_info.append(extract_sample_info(sample_name))
    
    if not count_data:
        raise ValueError("No valid count files found")
    
    # Create count matrix
    print("Building count matrix...")
    
    # Get all gene IDs (should be the same across samples)
    all_genes = set()
    for counts in count_data.values():
        all_genes.update(counts['gene_id'].tolist())
    
    all_genes = sorted(list(all_genes))
    print(f"Found {len(all_genes)} genes")
    
    # Create matrix
    count_matrix = pd.DataFrame(index=all_genes)
    
    for sample_name, counts in count_data.items():
        # Set gene_id as index for easier merging
        counts_indexed = counts.set_index('gene_id')
        count_matrix[sample_name] = counts_indexed['counts']
    
    # Fill missing values with 0
    count_matrix = count_matrix.fillna(0).astype(int)
    
    # Create sample metadata
    metadata = pd.DataFrame(sample_info)
    
    # Save outputs
    summary_file = output_dir / "count_summary.txt"
    
    # Save outputs
    print("Saving outputs...")
    
    # Count matrix
    count_matrix.to_csv(count_matrix_file, sep='\t')
    print(f"Count matrix saved: {count_matrix_file}")
    
    # Sample metadata
    metadata.to_csv(metadata_file, sep='\t', index=False)
    print(f"Sample metadata saved: {metadata_file}")
    
    # Summary statistics
    with open(summary_file, 'w') as f:
        f.write("Gene Count Matrix Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total samples: {len(count_data)}\n")
        f.write(f"Total genes: {len(all_genes)}\n")
        f.write(f"Matrix dimensions: {count_matrix.shape}\n\n")
        
        f.write("Sample information:\n")
        f.write("-" * 30 + "\n")
        for info in sample_info:
            f.write(f"{info['sample']}: {info['treatment']} - {info['condition']} (rep {info['replicate']})\n")
        
        f.write(f"\nCount statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total reads: {count_matrix.sum().sum():,}\n")
        f.write(f"Mean reads per sample: {count_matrix.sum().mean():,.0f}\n")
        f.write(f"Median reads per sample: {count_matrix.sum().median():,.0f}\n")
        f.write(f"Min reads per sample: {count_matrix.sum().min():,}\n")
        f.write(f"Max reads per sample: {count_matrix.sum().max():,}\n")
        
        f.write(f"\nGene statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Genes with zero counts: {(count_matrix.sum(axis=1) == 0).sum()}\n")
        f.write(f"Genes with < 10 total counts: {(count_matrix.sum(axis=1) < 10).sum()}\n")
        f.write(f"Genes with < 100 total counts: {(count_matrix.sum(axis=1) < 100).sum()}\n")
    
    print(f"Summary saved: {summary_file}")
    
    return count_matrix, metadata

def main():
    parser = argparse.ArgumentParser(
        description="Build gene count matrix from STAR or RSEM output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # STAR counts (default)
  python build_count_matrix.py star_counts/ -o count_matrices
  
  # RSEM counts
  python build_count_matrix.py rsem_counts/ -o count_matrices --type rsem
  
  # Auto-detect (prefers STAR if both found)
  python build_count_matrix.py counts/ -o count_matrices --type auto
        """
    )
    parser.add_argument("count_dir", help="Directory containing count files")
    parser.add_argument("-o", "--output", default="count_matrices", 
                       help="Output directory (default: count_matrices)")
    parser.add_argument("--type", choices=["star", "rsem", "auto"], default="auto",
                       help="Count file type: 'star' (ReadsPerGene.out.tab), "
                            "'rsem' (Trinity/RSEM .genes.results), or 'auto' (default: auto-detect)")
    
    args = parser.parse_args()
    
    try:
        count_matrix, metadata = build_count_matrix(args.count_dir, args.output, args.type)
        print("\n" + "="*60)
        print("SUCCESS: Gene count matrix created!")
        print("="*60)
        print(f"Matrix shape: {count_matrix.shape}")
        print(f"Samples: {list(count_matrix.columns)}")
        print(f"\nNext steps:")
        print(f"1. Review the count matrix: {args.output}/gene_count_matrix.tsv")
        print(f"2. Review sample metadata: {args.output}/sample_metadata.tsv")
        print(f"3. Run DESeq2/edgeR analysis using the provided R script")
        print(f"   OR run PyDESeq2 analysis:")
        print(f"   python pydeseq2_run_analysis.py {args.output}/gene_count_matrix.tsv {args.output}/sample_metadata.tsv -o pydeseq2_results")
        
    except Exception as e:
        print(f"ERROR: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
