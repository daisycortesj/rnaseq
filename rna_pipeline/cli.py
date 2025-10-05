
import argparse
from .main import run_workflow

def build_parser():
    ap = argparse.ArgumentParser(
        description="Choose STAR index (if reference present) else Trinity de novo"
    )
    ap.add_argument("--fasta", help="Reference genome FASTA (if available)")
    ap.add_argument("--gtf",   help="Annotation GTF/GFF (if available)")
    ap.add_argument("--reads-left",  help="Left FASTQ for Trinity (if no reference)")
    ap.add_argument("--reads-right", help="Right FASTQ for Trinity (if no reference)")
    ap.add_argument("--outdir", default="results/star_index", help="Output directory (STAR or Trinity)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--mem-gb",  type=int, default=64, help="RAM for Trinity")
    ap.add_argument("--readlen", type=int, default=150, help="Read length (STAR overhang)")
    ap.add_argument("--dry", action="store_true", help="Print commands without running")
    return ap

def main():
    args = build_parser().parse_args()
    run_workflow(args)

if __name__ == "__main__":
    main()
