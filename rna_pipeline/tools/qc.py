"""
Quality control tools for RNA-seq data using FastQC and MultiQC.
"""
import os
import subprocess
from pathlib import Path
from typing import List, Tuple


def check_qc_tools() -> Tuple[bool, List[str]]:
    """
    Check if FastQC and MultiQC are available.
    
    Returns:
        Tuple of (all_available, list_of_missing_tools)
    """
    missing = []
    
    try:
        subprocess.run(["fastqc", "--version"], 
                      capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append("fastqc")
    
    try:
        subprocess.run(["multiqc", "--version"], 
                      capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append("multiqc")
    
    return (len(missing) == 0, missing)


def install_qc_tools():
    """
    Install FastQC and MultiQC using conda.
    
    Returns:
        Return code from conda install command
    """
    print("[QC] Installing FastQC and MultiQC via conda...")
    cmd = [
        "conda", "install", 
        "-c", "bioconda",
        "fastqc", "multiqc",
        "-y"  # Auto-approve installation
    ]
    
    result = subprocess.run(cmd, capture_output=False)
    return result.returncode


def find_fastq_files(search_dir: str) -> List[str]:
    """
    Recursively find all FASTQ files in a directory.
    
    Args:
        search_dir: Directory to search for FASTQ files
        
    Returns:
        List of paths to FASTQ files (sorted)
    """
    search_path = Path(search_dir)
    
    # Find all fastq.gz and fq.gz files
    fastq_files = []
    for pattern in ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]:
        fastq_files.extend(search_path.rglob(pattern))
    
    # Convert to strings and sort
    return sorted([str(f) for f in fastq_files])


def build_fastqc_cmd(fastq_files: List[str], outdir: str, threads: int) -> List[str]:
    """
    Build FastQC command for quality control.
    
    Args:
        fastq_files: List of paths to FASTQ files
        outdir: Output directory for FastQC reports
        threads: Number of threads to use
        
    Returns:
        Command as list of strings
    """
    cmd = [
        "fastqc",
        "--outdir", outdir,
        "--threads", str(threads),
        "--quiet"  # Less verbose output
    ]
    
    # Add all FASTQ files
    cmd.extend(fastq_files)
    
    return cmd


def build_multiqc_cmd(fastqc_dir: str, outdir: str, title: str = "RNA-seq QC Report") -> List[str]:
    """
    Build MultiQC command to aggregate FastQC reports.
    
    Args:
        fastqc_dir: Directory containing FastQC reports
        outdir: Output directory for MultiQC report
        title: Title for the MultiQC report
        
    Returns:
        Command as list of strings
    """
    cmd = [
        "multiqc",
        fastqc_dir,
        "--outdir", outdir,
        "--filename", "multiqc_report",
        "--title", title,
        "--comment", "Quality control metrics for RNA-seq samples",
        "--force"  # Overwrite existing reports
    ]
    
    return cmd


def count_fastq_reads(fastq_file: str) -> int:
    """
    Count number of reads in a FASTQ file (works with gzipped files).
    
    Args:
        fastq_file: Path to FASTQ file (.fastq, .fq, .fastq.gz, .fq.gz)
        
    Returns:
        Number of reads (lines / 4)
    """
    import gzip
    
    try:
        # Determine if file is gzipped
        if fastq_file.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        with opener(fastq_file, mode) as f:
            line_count = sum(1 for _ in f)
        
        # FASTQ has 4 lines per read
        return line_count // 4
    
    except Exception as e:
        print(f"Warning: Could not count reads in {fastq_file}: {e}")
        return 0


def summarize_fastq_files(fastq_files: List[str]) -> dict:
    """
    Generate summary statistics for FASTQ files.
    
    Args:
        fastq_files: List of FASTQ file paths
        
    Returns:
        Dictionary with summary stats
    """
    summary = {
        'total_files': len(fastq_files),
        'total_reads': 0,
        'files': []
    }
    
    print(f"[QC] Analyzing {len(fastq_files)} FASTQ files...")
    
    for fq in fastq_files:
        # Get file size
        size_mb = os.path.getsize(fq) / (1024 * 1024)
        
        # Count reads (can be slow for large files, so make it optional)
        # reads = count_fastq_reads(fq)
        # summary['total_reads'] += reads
        
        file_info = {
            'path': fq,
            'name': os.path.basename(fq),
            'size_mb': round(size_mb, 2),
            # 'reads': reads
        }
        
        summary['files'].append(file_info)
    
    return summary
