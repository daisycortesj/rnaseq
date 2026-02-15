#!/usr/bin/env python3
"""
===============================================================================
parse_refseq_gtf.py 
===============================================================================

This script parses RefSeq/Gnomon GTF annotation files and extracts structured
data into an easy-to-read TSV (tab-separated values) file.

WHAT IS A GTF FILE?
A GTF (Gene Transfer Format) file contains genome annotations. Each line
describes a genomic feature (gene, transcript, CDS, etc.) with 9 columns.

OUTPUT FORMAT:
The script produces a simplified 4-column TSV file:
    1. local_number   - Unique index (1, 2, 3, ...) for each row (KEY COLUMN)
    2. gene_id        - Gene identifier (e.g., LOC135151205)
    3. transcript_id  - Transcript XM number (e.g., XM_064088976.1)
    4. description    - Gene/product description

USAGE EXAMPLES:
    python parse_refseq_gtf.py --gtf annotation.gtf --feature gene --out genes_parsed.tsv
    python parse_refseq_gtf.py --gtf annotation.gtf --feature transcript --out transcripts_parsed.tsv
    python parse_refseq_gtf.py --gtf annotation.gtf --feature CDS --out cds_parsed.tsv

Author: Daisy Cortes
===============================================================================
"""

# ============================================================================
# IMPORTS - These are Python modules we need
# ============================================================================
# argparse: helps us handle command-line arguments (--gtf, --out, --feature)
# sys: provides system functions like exit() and stderr for error messages
import argparse
import sys


# ============================================================================
# FUNCTION 1: Parse GTF Attributes (Column 9)
# ============================================================================
def parse_gtf_attributes(attr_string):
    """
    This function takes the messy attribute string from column 9 of a GTF
    and turns it into a clean Python dictionary.
    
    EXAMPLE INPUT:
        gene_id "LOC135151205"; transcript_id "XM_064088976.1"; gene "LOC135151205";
    
    EXAMPLE OUTPUT:
        {
            'gene_id': 'LOC135151205',
            'transcript_id': 'XM_064088976.1',
            'gene': 'LOC135151205'
        }
    
    HOW IT WORKS:
    1. Split the string by semicolons (;) to get individual attributes
    2. For each attribute, split by the first space to separate key and value
    3. Remove the quotes around values
    4. Store in a dictionary for easy access
    """
    
    # Step 1: Create an empty dictionary to store our results
    # A dictionary in Python is like a lookup table: {key: value, key: value}
    attributes = {}
    
    # Step 2: Split the string by semicolons
    # .split(';') breaks a string into a list at each semicolon
    # Example: "a; b; c" becomes ["a", " b", " c"]
    parts = attr_string.split(';')
    
    # Step 3: Loop through each part
    for part in parts:
        # .strip() removes leading/trailing whitespace
        # Example: "  hello  " becomes "hello"
        part = part.strip()
        
        # Skip empty parts (happens when there's a trailing semicolon)
        if not part:
            continue
        
        # Step 4: Split on the FIRST space only to separate key from value
        # Example: 'gene_id "LOC135151205"' splits into ['gene_id', '"LOC135151205"']
        # The '1' means "split only once" so descriptions with spaces stay intact
        tokens = part.split(' ', 1)
        
        # Make sure we got both a key and a value
        if len(tokens) == 2:
            key = tokens[0].strip()      # The attribute name (e.g., "gene_id")
            value = tokens[1].strip()    # The attribute value (e.g., '"LOC135151205"')
            
            # Step 5: Remove the surrounding quotes from the value
            # GTF format puts values in quotes: "value"
            # We check if it starts AND ends with quotes, then remove them
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]  # [1:-1] means "from position 1 to second-to-last"
            
            # Step 6: Store in our dictionary
            attributes[key] = value
    
    return attributes


# ============================================================================
# FUNCTION 2: Parse One Line of GTF File
# ============================================================================
def parse_gtf_line(line):
    """
    Parse a single line from a GTF file into a Python dictionary.
    
    GTF FORMAT (9 columns separated by tabs):
    1. seqname     - chromosome/contig name (e.g., NC_030381.2)
    2. source      - who annotated this (e.g., Gnomon, RefSeq)
    3. feature     - what type of feature (e.g., gene, transcript, CDS)
    4. start       - start position (1-based)
    5. end         - end position (inclusive)
    6. score       - usually a dot (.)
    7. strand      - plus (+) or minus (-)
    8. frame       - reading frame (0, 1, 2, or .)
    9. attributes  - semicolon-separated key-value pairs
    
    RETURNS:
        A dictionary with all the parsed information, or None if invalid
    """
    
    # Step 1: Remove newline and split by tabs
    # .strip() removes the newline character at the end
    # .split('\t') splits on tab characters
    fields = line.strip().split('\t')
    
    # Step 2: Check if we have exactly 9 columns (GTF requirement)
    if len(fields) != 9:
        return None  # Invalid line, skip it
    
    # Step 3: Create a dictionary with the first 8 columns
    # We use index numbers: fields[0] is first column, fields[1] is second, etc.
    record = {
        'seqname': fields[0],   # Column 1: chromosome/scaffold
        'source': fields[1],    # Column 2: source (Gnomon, RefSeq, etc.)
        'feature': fields[2],   # Column 3: feature type (gene, CDS, etc.)
        'start': fields[3],     # Column 4: start position
        'end': fields[4],       # Column 5: end position
        'score': fields[5],     # Column 6: score (usually ".")
        'strand': fields[6],    # Column 7: strand (+ or -)
        'frame': fields[7],     # Column 8: reading frame
    }
    
    # Step 4: Parse the complex attributes column (column 9)
    attributes = parse_gtf_attributes(fields[8])
    
    # Step 5: Add the attributes to our record dictionary
    # .update() merges two dictionaries together
    record.update(attributes)
    
    return record


# ============================================================================
# FUNCTION 3: Extract Output Fields in Correct Order
# ============================================================================
def extract_output_fields(record, local_number):
    """
    Pull out specific fields from our parsed record in the correct order
    for our output TSV file.
    
    We want these columns (in this exact order):
        local_number, gene_id, transcript_id, description
    
    IMPORTANT: 
    - local_number is a unique index (1, 2, 3, ...) for each row
    - gene_id is the gene identifier (e.g., "LOC135151205")
    - transcript_id is the XM number (e.g., "XM_064088976.1")
    - description is the gene/product description
    
    For gene rows: transcript_id is usually empty ""
    For transcript/CDS rows: transcript_id has the XM number
    """
    
    # Define the exact order of columns we want (simplified output)
    output_fields = [
        'gene_id',       # Gene identifier (e.g., LOC135151205)
        'transcript_id', # Transcript XM number (may be empty for gene rows)
        'description',   # Gene/product description
    ]
    
    # Start with the local number (unique index)
    values = [str(local_number)]
    
    # Extract each field from the record
    # .get(field, '') means "get this field, or use empty string if missing"
    for field in output_fields:
        value = record.get(field, '')  # Get value or empty string
        values.append(value)
    
    return values


# ============================================================================
# FUNCTION 4: Parse Entire GTF File
# ============================================================================
def parse_gtf_file(gtf_path, feature_filter='gene'):
    """
    Read through the entire GTF file and extract rows matching our feature type.
    
    INPUTS:
        gtf_path: Path to the GTF file (e.g., "annotation.gtf")
        feature_filter: What type of rows to keep (e.g., "gene", "transcript", "CDS")
    
    RETURNS:
        - parsed_records: List of dictionaries (one per matching row)
        - total_lines: How many data lines we read (excluding comments)
        - filtered_lines: How many lines matched our feature filter
    """
    
    # Create empty list to store our results
    parsed_records = []
    
    # Initialize counters
    total_lines = 0      # Total non-comment lines
    filtered_lines = 0   # Lines that match our feature filter
    
    # Print status message
    print(f"Reading GTF file: {gtf_path}")
    print(f"Filtering for feature type: {feature_filter}")
    
    # Open the file and read it line by line
    # 'with' automatically closes the file when done
    with open(gtf_path, 'r') as f:
        for line in f:
            # Skip comment lines that start with #
            if line.startswith('#'):
                continue
            
            # Skip empty lines
            if not line.strip():
                continue
            
            # Count this data line
            total_lines += 1
            
            # Parse this line
            record = parse_gtf_line(line)
            
            # If parsing failed, print a warning and skip
            if record is None:
                print(f"Warning: Could not parse line {total_lines}", file=sys.stderr)
                continue
            
            # Check if this row matches our feature filter
            # Example: if we want "gene" rows, only keep rows where feature == "gene"
            if record['feature'] == feature_filter:
                filtered_lines += 1
                parsed_records.append(record)
    
    return parsed_records, total_lines, filtered_lines


# ============================================================================
# FUNCTION 5: Write TSV Output File
# ============================================================================
def write_tsv_output(records, output_path):
    """
    Write our parsed records to a TSV (tab-separated values) file.
    TSV is like CSV but uses tabs instead of commas.
    
    The output file will have:
    - A header row with column names: local_number, gene_id, transcript_id, description
    - One row per record
    - Tabs separating columns
    - local_number is a unique index (1, 2, 3, ...) serving as the key column
    """
    
    # Define our column headers (simplified output)
    headers = [
        'local_number',   # Unique index (key column)
        'gene_id',        # Gene identifier
        'transcript_id',  # Transcript XM number
        'description'     # Gene/product description
    ]
    
    # Open output file for writing
    with open(output_path, 'w') as f:
        # Write header row
        # '\t'.join(headers) puts tabs between column names
        # '\n' adds a newline at the end
        f.write('\t'.join(headers) + '\n')
        
        # Write each record as a row with unique index
        # enumerate() gives us both the index and the item
        # We start counting from 1 (not 0) by using start=1
        for index, record in enumerate(records, start=1):
            # Get values in correct order, passing the unique local_number
            values = extract_output_fields(record, index)
            
            # Join with tabs and write
            f.write('\t'.join(values) + '\n')


# ============================================================================
# FUNCTION 6: Count Unique Genes
# ============================================================================
def count_unique_genes(records):
    """
    Count how many unique gene_id values we have in our data.
    
    Uses a Python 'set' which automatically keeps only unique values.
    For example: {1, 2, 2, 3} becomes {1, 2, 3}
    """
    
    # Create an empty set (like a list but only unique values)
    gene_ids = set()
    
    # Loop through all records
    for record in records:
        # Get the gene_id (or empty string if missing)
        gene_id = record.get('gene_id', '')
        
        # Only add non-empty gene_ids
        if gene_id:
            gene_ids.add(gene_id)  # Sets ignore duplicates automatically
    
    # Return the count
    return len(gene_ids)


# ============================================================================
# MAIN FUNCTION: This runs when you execute the script
# ============================================================================
def main():
    """
    This is the main function that:
    1. Parses command-line arguments
    2. Reads the GTF file
    3. Filters for the feature type we want
    4. Writes the output TSV file
    5. Prints a summary
    """
    
    # ========================================================================
    # STEP 1: Set up command-line argument parser
    # ========================================================================
    # argparse makes it easy to handle command-line arguments like --gtf, --out, etc.
    parser = argparse.ArgumentParser(
        description='Parse RefSeq/Gnomon GTF files and extract structured data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract gene-level annotations
  python parse_refseq_gtf.py --gtf annotation.gtf --feature gene --out genes_parsed.tsv
  
  # Extract transcript-level annotations
  python parse_refseq_gtf.py --gtf annotation.gtf --feature transcript --out transcripts_parsed.tsv
  
  # Extract CDS annotations
  python parse_refseq_gtf.py --gtf annotation.gtf --feature CDS --out cds_parsed.tsv
        """
    )
    
    # Define the --gtf argument (REQUIRED)
    parser.add_argument(
        '--gtf',
        required=True,
        help='Path to input GTF file'
    )
    
    # Define the --out argument (REQUIRED)
    parser.add_argument(
        '--out',
        required=True,
        help='Path to output TSV file'
    )
    
    # Define the --feature argument (OPTIONAL, defaults to "gene")
    parser.add_argument(
        '--feature',
        default='gene',
        choices=['gene', 'transcript', 'CDS', 'exon', 'start_codon', 'stop_codon'],
        help='Feature type to extract (default: gene)'
    )
    
    # ========================================================================
    # STEP 2: Parse the command-line arguments
    # ========================================================================
    args = parser.parse_args()
    
    # ========================================================================
    # STEP 3: Parse the GTF file
    # ========================================================================
    # Use try/except to handle errors gracefully
    try:
        records, total_lines, filtered_lines = parse_gtf_file(args.gtf, args.feature)
    except FileNotFoundError:
        # If file doesn't exist, print error and exit
        print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # If any other error occurs, print it and exit
        print(f"Error reading GTF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # ========================================================================
    # STEP 4: Count unique genes
    # ========================================================================
    unique_genes = count_unique_genes(records)
    
    # ========================================================================
    # STEP 5: Write output file
    # ========================================================================
    try:
        write_tsv_output(records, args.out)
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # ========================================================================
    # STEP 6: Print summary
    # ========================================================================
    print("\n" + "="*60)
    print("PARSING COMPLETE")
    print("="*60)
    print(f"Total non-comment rows read:     {total_lines}")
    print(f"Rows kept after feature filter:  {filtered_lines}")
    print(f"Unique gene_id values:           {unique_genes}")
    print(f"\nOutput columns:")
    print(f"  1. local_number (unique index)")
    print(f"  2. gene_id")
    print(f"  3. transcript_id (XM number)")
    print(f"  4. description")
    print(f"\nOutput written to:               {args.out}")
    print("="*60)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
# This special if statement means: "only run main() if this script is executed directly"
# (not if it's imported as a module in another script)
if __name__ == '__main__':
    main()
