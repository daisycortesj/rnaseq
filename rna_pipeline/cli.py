import argparse
from .main import run_workflow

def build_parser():
    ap = argparse.ArgumentParser(
        description="RNA-seq pipeline: STAR indexing for reference-based analysis or Trinity for de novo assembly"
    )
    # Input files
    ap.add_argument("--fasta", help="Reference genome FASTA file")
    ap.add_argument("--gtf",   help="Gene annotation GTF/GFF file")
    ap.add_argument("--reads-left",  help="Left FASTQ reads for Trinity de novo assembly")
    ap.add_argument("--reads-right", help="Right FASTQ reads for Trinity de novo assembly")
    
    # Output settings
    ap.add_argument("--outdir", default="results", help="Output directory for STAR index or Trinity assembly")
    
    # Resource settings (optimized for larger genomes)
    ap.add_argument("--threads", type=int, default=16, help="Number of CPU threads (default: 16)")
    ap.add_argument("--mem-gb",  type=int, default=128, help="RAM in GB for Trinity (default: 128GB)")
    ap.add_argument("--readlen", type=int, default=150, help="Read length for STAR overhang calculation (default: 150bp)")
    
    # Execution options
    ap.add_argument("--dry", action="store_true", help="Print commands without running (for testing)")
    ap.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    return ap

def main():
    args = build_parser().parse_args()
    run_workflow(args)

if __name__ == "__main__":
    main()
