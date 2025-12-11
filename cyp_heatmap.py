#!/usr/bin/env python3
"""
CYP Gene Heatmap Generator

This script creates a heatmap showing expression levels of CYP 
(Cytochrome P450) genes across your samples.

WHAT ARE CYP GENES?
CYP genes encode cytochrome P450 enzymes, which are important for:
- Drug metabolism
- Hormone synthesis
- Detoxification of foreign compounds

HOW TO USE:
    python cyp_heatmap.py <count_matrix> <metadata> [output_dir]

EXAMPLE:
    python cyp_heatmap.py gene_count_matrix.tsv sample_metadata.tsv results/

INPUTS:
    count_matrix: A TSV file with genes as rows and samples as columns
    metadata: A TSV file with sample information (must have 'sample' column)
    output_dir: Where to save the heatmap (default: cyp_heatmap_results)
"""

# =============================================================================
# STEP 1: IMPORT REQUIRED LIBRARIES
# =============================================================================
# These are the tools we need to run this script

import os                    # For working with file paths
import sys                   # For command line arguments and exit codes
import argparse              # For parsing command line arguments
import pandas as pd          # For working with data tables (DataFrames)
import numpy as np           # For numerical calculations
from pathlib import Path     # For easier file path handling

# Suppress warnings to keep output clean
import warnings
warnings.filterwarnings('ignore')

# Try to import plotting libraries
try:
    import matplotlib.pyplot as plt  # For creating plots
    import seaborn as sns            # For pretty heatmaps
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install: pip install matplotlib seaborn")
    sys.exit(1)


# =============================================================================
# STEP 2: DEFINE HELPER FUNCTIONS
# =============================================================================

def read_count_matrix(count_file):
    """
    Read the gene count matrix from a TSV file.
    
    The count matrix has:
    - Rows = genes (gene names are in the first column)
    - Columns = samples (sample names are in the header)
    - Values = read counts (how many times each gene was detected)
    
    Parameters:
        count_file: Path to the TSV file containing gene counts
        
    Returns:
        A pandas DataFrame with genes as rows and samples as columns
    """
    print(f"Reading count matrix: {count_file}")
    
    # Read the file (sep='\t' means it's tab-separated)
    # index_col=0 means the first column contains row names (gene names)
    count_matrix = pd.read_csv(count_file, sep='\t', index_col=0)
    
    print(f"  Found {len(count_matrix)} genes")
    print(f"  Found {len(count_matrix.columns)} samples")
    
    return count_matrix


def read_metadata(metadata_file):
    """
    Read the sample metadata from a TSV file.
    
    The metadata contains information about each sample, like:
    - sample: The sample name (must match column names in count matrix)
    - treatment: What treatment was applied (e.g., 'control' vs 'treated')
    - group: Sample grouping information
    
    Parameters:
        metadata_file: Path to the TSV file containing sample information
        
    Returns:
        A pandas DataFrame with sample information
    """
    print(f"Reading metadata: {metadata_file}")
    
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    print(f"  Found {len(metadata)} samples in metadata")
    print(f"  Columns: {list(metadata.columns)}")
    
    return metadata


def find_cyp_genes(count_matrix):
    """
    Find all genes that contain 'CYP' in their name.
    
    CYP genes are named with 'CYP' followed by numbers and letters,
    like CYP1A1, CYP2D6, CYP71D55, etc.
    
    Parameters:
        count_matrix: DataFrame with gene names as row index
        
    Returns:
        A list of gene names that contain 'CYP'
    """
    print("\nSearching for CYP genes...")
    
    # Get all gene names from the row index
    all_genes = count_matrix.index.tolist()
    
    # Find genes that contain 'CYP' (case-insensitive search)
    # str.upper() converts to uppercase so we catch CYP, Cyp, cyp, etc.
    cyp_genes = [gene for gene in all_genes if 'CYP' in str(gene).upper()]
    
    print(f"  Found {len(cyp_genes)} CYP genes out of {len(all_genes)} total genes")
    
    # Show first few CYP genes found
    if len(cyp_genes) > 0:
        print(f"  Examples: {cyp_genes[:5]}")
    
    return cyp_genes


def normalize_counts(count_matrix):
    """
    Normalize the count data to make samples comparable.
    
    WHY NORMALIZE?
    Different samples may have different total reads (sequencing depth).
    Normalization accounts for this so we can compare expression levels.
    
    This uses a simple method: divide by total counts, multiply by 1 million
    (this gives "Counts Per Million" or CPM)
    
    Parameters:
        count_matrix: DataFrame with raw counts
        
    Returns:
        DataFrame with normalized counts
    """
    print("\nNormalizing counts...")
    
    # Calculate total counts for each sample (sum down each column)
    total_counts = count_matrix.sum(axis=0)
    
    # Divide each value by the sample's total counts, multiply by 1 million
    # This gives us "Counts Per Million" (CPM)
    normalized = (count_matrix / total_counts) * 1_000_000
    
    print(f"  Sample total counts range: {total_counts.min():.0f} to {total_counts.max():.0f}")
    
    return normalized


def create_cyp_heatmap(count_matrix, metadata, cyp_genes, output_dir):
    """
    Create a heatmap showing CYP gene expression across samples.
    
    WHAT IS A HEATMAP?
    A heatmap shows data as colors in a grid:
    - Each row is a gene
    - Each column is a sample
    - Color intensity shows expression level
    - Red = higher expression, Blue = lower expression
    
    Parameters:
        count_matrix: DataFrame with normalized counts
        metadata: DataFrame with sample information
        cyp_genes: List of CYP gene names to include
        output_dir: Where to save the heatmap
    """
    print("\nCreating CYP gene heatmap...")
    
    # ----- STEP A: Prepare the data -----
    
    # Select only the CYP genes from our count matrix
    cyp_data = count_matrix.loc[cyp_genes].copy()
    print(f"  Selected {len(cyp_data)} CYP genes")
    
    # Log-transform the data (log2)
    # WHY? Gene expression data is often skewed, log helps visualize it better
    # We add 1 to avoid log(0) which is undefined
    cyp_data_log = np.log2(cyp_data + 1)
    
    # Center the data: subtract the mean of each gene
    # WHY? This shows relative expression (above/below average) rather than absolute levels
    # axis=1 means calculate mean across columns (samples) for each row (gene)
    gene_means = cyp_data_log.mean(axis=1)
    cyp_data_centered = cyp_data_log.subtract(gene_means, axis=0)
    
    print(f"  Data transformed: log2(counts + 1), then centered")
    
    # ----- STEP B: Set up sample colors based on treatment/group -----
    
    # Find which column to use for coloring samples
    color_column = None
    if 'treatment' in metadata.columns:
        color_column = 'treatment'
    elif 'group' in metadata.columns:
        color_column = 'group'
    elif 'condition' in metadata.columns:
        color_column = 'condition'
    
    # Create color mapping for samples if we found a grouping column
    col_colors = None
    if color_column is not None:
        print(f"  Using '{color_column}' to color samples")
        
        # Make sure metadata sample names match our data columns
        metadata_indexed = metadata.set_index('sample')
        metadata_indexed = metadata_indexed.reindex(cyp_data_centered.columns)
        
        # Create a color for each unique group/treatment
        unique_groups = metadata_indexed[color_column].unique()
        colors = sns.color_palette("Set2", len(unique_groups))
        color_dict = dict(zip(unique_groups, colors))
        
        # Map samples to their colors
        col_colors = metadata_indexed[color_column].map(color_dict)
        
        print(f"  Groups found: {list(unique_groups)}")
    
    # ----- STEP C: Create the heatmap -----
    
    # Set up the figure size based on number of genes and samples
    # More genes = taller figure, more samples = wider figure
    fig_height = max(8, len(cyp_genes) * 0.3)  # At least 8 inches tall
    fig_width = max(10, len(cyp_data_centered.columns) * 0.8)  # At least 10 inches wide
    
    print(f"  Figure size: {fig_width:.1f} x {fig_height:.1f} inches")
    
    # Create a clustered heatmap using seaborn
    # Clustering groups similar genes and samples together
    g = sns.clustermap(
        cyp_data_centered,              # The data to plot
        cmap='RdBu_r',                  # Color scheme: Red-Blue reversed (red=high)
        center=0,                        # Center the colormap at 0
        col_colors=col_colors,          # Color bar showing sample groups
        figsize=(fig_width, fig_height),
        method='ward',                   # Clustering method
        metric='euclidean',              # Distance metric for clustering
        yticklabels=True,               # Show gene names on y-axis
        xticklabels=True,               # Show sample names on x-axis
        cbar_kws={'label': 'Log2 Expression (centered)'},  # Colorbar label
        dendrogram_ratio=(0.15, 0.1),   # Size of clustering trees
        linewidths=0.5,                 # Lines between cells
        linecolor='lightgray'
    )
    
    # ----- STEP D: Make it look nice -----
    
    # Add a title
    g.fig.suptitle('CYP Gene Expression Heatmap', fontsize=14, fontweight='bold', y=1.02)
    
    # Rotate sample labels for readability
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), 
        rotation=45, 
        ha='right',
        fontsize=9
    )
    
    # Format gene labels
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), 
        rotation=0, 
        fontsize=8
    )
    
    # Add a legend if we have sample colors
    if col_colors is not None and color_column is not None:
        from matplotlib.patches import Patch
        
        legend_elements = [
            Patch(facecolor=color_dict[group], label=group)
            for group in unique_groups
        ]
        
        g.fig.legend(
            handles=legend_elements,
            title=color_column.capitalize(),
            loc='upper right',
            bbox_to_anchor=(1.15, 0.9),
            fontsize=9
        )
    
    # ----- STEP E: Save the heatmap -----
    
    output_file = output_dir / "cyp_heatmap.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  Saved heatmap: {output_file}")
    
    # Also save as PNG for easier viewing
    output_file_png = output_dir / "cyp_heatmap.png"
    
    # Recreate for PNG (need to do this because we closed the previous figure)
    g = sns.clustermap(
        cyp_data_centered,
        cmap='RdBu_r',
        center=0,
        col_colors=col_colors,
        figsize=(fig_width, fig_height),
        method='ward',
        metric='euclidean',
        yticklabels=True,
        xticklabels=True,
        cbar_kws={'label': 'Log2 Expression (centered)'},
        dendrogram_ratio=(0.15, 0.1),
        linewidths=0.5,
        linecolor='lightgray'
    )
    g.fig.suptitle('CYP Gene Expression Heatmap', fontsize=14, fontweight='bold', y=1.02)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    
    plt.savefig(output_file_png, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  Saved PNG: {output_file_png}")
    
    return cyp_data_centered


def save_cyp_gene_list(cyp_genes, cyp_data, output_dir):
    """
    Save the list of CYP genes found to a text file.
    
    This is useful for:
    - Checking which genes were included
    - Using the list in other analyses
    """
    # Save just the gene names
    gene_list_file = output_dir / "cyp_genes_found.txt"
    with open(gene_list_file, 'w') as f:
        f.write("# CYP genes found in your data\n")
        f.write(f"# Total: {len(cyp_genes)} genes\n\n")
        for gene in sorted(cyp_genes):
            f.write(f"{gene}\n")
    
    print(f"\nSaved gene list: {gene_list_file}")
    
    # Save the expression data as a table
    expression_file = output_dir / "cyp_expression_data.tsv"
    cyp_data.to_csv(expression_file, sep='\t')
    
    print(f"Saved expression data: {expression_file}")


# =============================================================================
# STEP 3: MAIN FUNCTION - TIES EVERYTHING TOGETHER
# =============================================================================

def main():
    """
    Main function that runs when you execute this script.
    
    This coordinates all the steps:
    1. Parse command line arguments
    2. Read the input files
    3. Find CYP genes
    4. Normalize and create heatmap
    5. Save results
    """
    
    # ----- Set up command line argument parser -----
    parser = argparse.ArgumentParser(
        description="Create a heatmap of CYP gene expression",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLE USAGE:
    python cyp_heatmap.py gene_counts.tsv metadata.tsv -o results/
    
This will create:
    results/cyp_heatmap.pdf       - The heatmap as a PDF
    results/cyp_heatmap.png       - The heatmap as a PNG
    results/cyp_genes_found.txt   - List of CYP genes found
    results/cyp_expression_data.tsv - Expression data for CYP genes
        """
    )
    
    # Required arguments
    parser.add_argument(
        "count_matrix", 
        help="Path to gene count matrix TSV file (genes as rows, samples as columns)"
    )
    parser.add_argument(
        "metadata", 
        help="Path to sample metadata TSV file (must have 'sample' column)"
    )
    
    # Optional arguments
    parser.add_argument(
        "-o", "--output", 
        default="cyp_heatmap_results",
        help="Output directory (default: cyp_heatmap_results)"
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # ----- Validate input files exist -----
    if not os.path.exists(args.count_matrix):
        print(f"ERROR: Count matrix file not found: {args.count_matrix}")
        return 1
    
    if not os.path.exists(args.metadata):
        print(f"ERROR: Metadata file not found: {args.metadata}")
        return 1
    
    # ----- Create output directory -----
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # ----- Print header -----
    print("=" * 60)
    print("CYP Gene Heatmap Generator")
    print("=" * 60)
    print(f"Count matrix: {args.count_matrix}")
    print(f"Metadata: {args.metadata}")
    print(f"Output directory: {output_dir}")
    print("=" * 60)
    
    # ----- Run the analysis -----
    try:
        # Step 1: Read the data files
        count_matrix = read_count_matrix(args.count_matrix)
        metadata = read_metadata(args.metadata)
        
        # Step 2: Find CYP genes
        cyp_genes = find_cyp_genes(count_matrix)
        
        # Check if we found any CYP genes
        if len(cyp_genes) == 0:
            print("\nERROR: No CYP genes found in your data!")
            print("Make sure your gene names contain 'CYP' (case-insensitive)")
            print("Example gene names: CYP1A1, CYP2D6, CYP71D55")
            return 1
        
        # Step 3: Normalize the counts
        normalized_counts = normalize_counts(count_matrix)
        
        # Step 4: Create the heatmap
        cyp_data = create_cyp_heatmap(
            normalized_counts, 
            metadata, 
            cyp_genes, 
            output_dir
        )
        
        # Step 5: Save the gene list and expression data
        save_cyp_gene_list(cyp_genes, cyp_data, output_dir)
        
        # ----- Print completion message -----
        print("\n" + "=" * 60)
        print("SUCCESS! Analysis complete.")
        print("=" * 60)
        print(f"\nResults saved to: {output_dir}/")
        print("\nOutput files:")
        print("  - cyp_heatmap.pdf       : High-quality heatmap (for publications)")
        print("  - cyp_heatmap.png       : Heatmap image (for quick viewing)")
        print("  - cyp_genes_found.txt   : List of CYP genes in your data")
        print("  - cyp_expression_data.tsv : Expression values for CYP genes")
        print("=" * 60)
        
        return 0  # Success!
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1  # Error occurred


# =============================================================================
# STEP 4: RUN THE SCRIPT
# =============================================================================

# This block runs when you execute the script directly (not when importing it)
if __name__ == "__main__":
    # Call main() and exit with its return code
    sys.exit(main())