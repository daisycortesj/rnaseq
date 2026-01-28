
import argparse
from .main import run_workflow

def build_parser():
    ap = argparse.ArgumentParser(
        description="RNA-seq pipeline: STAR index building, STAR alignment, or Trinity de novo assembly"
    )
    
    # Mode selection
    ap.add_argument("--mode", choices=["index", "align", "trinity", "qc", "auto"], default="auto",
                    help="Pipeline mode: 'index' (build STAR index), 'align' (align reads), 'trinity' (de novo assembly), 'qc' (quality control), 'auto' (detect from inputs)")
    
    # Config file
    ap.add_argument("--config", help="Path to config.yaml file (default: auto-detect)")
    
    # Genome selection (config-driven approach)
    ap.add_argument("--genome", help="Genome type from config (e.g., 'carrot', 'nutmeg')")
    
    # Index building arguments
    ap.add_argument("--fasta", help="Reference genome FASTA (for index building)")
    ap.add_argument("--gtf",   help="Annotation GTF/GFF (optional for index building)")
    ap.add_argument("--readlen", type=int, default=150, help="Read length for STAR splice junction overhang (index building)")
    
    # Alignment arguments
    ap.add_argument("--genome-index", help="Path to STAR genome index directory (for alignment)")
    ap.add_argument("--sample-name", help="Sample name for output prefix (for alignment)")
    ap.add_argument("--quant-mode", action="store_true", default=None, 
                    help="Enable gene quantification in alignment (requires GTF in index)")
    ap.add_argument("--no-quant-mode", action="store_false", dest="quant_mode",
                    help="Disable gene quantification (use for genomes without GTF)")
    
    # Read files (used for alignment, Trinity, or QC)
    ap.add_argument("--reads-left",  help="R1/forward reads FASTQ (for alignment, Trinity, or QC)")
    ap.add_argument("--reads-right", help="R2/reverse reads FASTQ (for paired-end alignment, Trinity, or QC)")
    
    # QC-specific arguments
    ap.add_argument("--fastq-dir", help="Directory to search for all FASTQ files (for QC mode)")
    
    # General arguments
    ap.add_argument("--outdir", default="results", help="Output directory")
    ap.add_argument("--threads", type=int, default=8, help="Number of threads")
    ap.add_argument("--mem-gb",  type=int, default=64, help="RAM for Trinity (GB)")
    ap.add_argument("--resume", action="store_true", help="Resume incomplete Trinity assembly from checkpoint")
    ap.add_argument("--dry", action="store_true", help="Print commands without running")
    
    return ap

def main():
    args = build_parser().parse_args()
    run_workflow(args)

if __name__ == "__main__":
    main()
