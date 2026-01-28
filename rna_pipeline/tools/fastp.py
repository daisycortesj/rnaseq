"""
Read trimming and filtering using fastp.
"""
import os
from typing import List


def build_fastp_cmd(reads_left: str, reads_right: str, 
                   out_left: str, out_right: str,
                   threads: int = 8,
                   min_quality: int = 20,
                   min_length: int = 50) -> List[str]:
    """
    Build fastp command for paired-end read trimming/filtering.
    
    Args:
        reads_left: Path to R1/forward reads
        reads_right: Path to R2/reverse reads
        out_left: Output path for trimmed R1 reads
        out_right: Output path for trimmed R2 reads
        threads: Number of threads (default: 8)
        min_quality: Minimum base quality score (default: 20)
        min_length: Minimum read length after trimming (default: 50)
        
    Returns:
        Command as list of strings
        
    Note:
        fastp automatically detects and trims adapters for paired-end data!
    """
    # Generate JSON report filename from output prefix
    out_dir = os.path.dirname(out_left)
    sample_name = os.path.basename(out_left).replace('_R1_trimmed.fq.gz', '').replace('_1_trimmed.fq.gz', '')
    json_report = os.path.join(out_dir, f"{sample_name}_fastp.json")
    html_report = os.path.join(out_dir, f"{sample_name}_fastp.html")
    
    cmd = [
        "fastp",
        # Input files
        "-i", reads_left,
        "-I", reads_right,
        # Output files
        "-o", out_left,
        "-O", out_right,
        # Reports
        "--json", json_report,
        "--html", html_report,
        # Adapter detection (automatic for paired-end!)
        "--detect_adapter_for_pe",
        # Quality filtering
        "--qualified_quality_phred", str(min_quality),
        # Length filtering
        "--length_required", str(min_length),
        # Performance
        "--thread", str(threads),
        # Disable unnecessary features for speed
        "--disable_quality_filtering",  # We set quality via qualified_quality_phred
        "--disable_length_filtering",   # We set length via length_required
    ]
    
    return cmd


def build_fastp_single_end_cmd(reads: str, out_reads: str,
                               threads: int = 8,
                               min_quality: int = 20,
                               min_length: int = 50) -> List[str]:
    """
    Build fastp command for single-end read trimming/filtering.
    
    Args:
        reads: Path to input reads
        out_reads: Output path for trimmed reads
        threads: Number of threads (default: 8)
        min_quality: Minimum base quality score (default: 20)
        min_length: Minimum read length after trimming (default: 50)
        
    Returns:
        Command as list of strings
    """
    out_dir = os.path.dirname(out_reads)
    sample_name = os.path.basename(out_reads).replace('_trimmed.fq.gz', '')
    json_report = os.path.join(out_dir, f"{sample_name}_fastp.json")
    html_report = os.path.join(out_dir, f"{sample_name}_fastp.html")
    
    cmd = [
        "fastp",
        # Input file
        "-i", reads,
        # Output file
        "-o", out_reads,
        # Reports
        "--json", json_report,
        "--html", html_report,
        # Quality filtering
        "--qualified_quality_phred", str(min_quality),
        # Length filtering
        "--length_required", str(min_length),
        # Performance
        "--thread", str(threads),
    ]
    
    return cmd


def get_trimmed_filename(original_path: str, suffix: str = "_trimmed") -> str:
    """
    Generate output filename for trimmed reads.
    
    Args:
        original_path: Path to original FASTQ file
        suffix: Suffix to add before extension (default: "_trimmed")
        
    Returns:
        Path to trimmed file
        
    Example:
        >>> get_trimmed_filename("/path/sample_R1.fq.gz")
        "/path/sample_R1_trimmed.fq.gz"
    """
    base, ext = os.path.splitext(original_path)
    
    # Handle .fq.gz and .fastq.gz
    if ext == '.gz':
        base2, ext2 = os.path.splitext(base)
        if ext2 in ['.fq', '.fastq']:
            return f"{base2}{suffix}{ext2}{ext}"
    
    # Fallback
    return f"{base}{suffix}{ext}"


def summarize_fastp_report(json_path: str) -> dict:
    """
    Parse fastp JSON report and extract key statistics.
    
    Args:
        json_path: Path to fastp JSON report
        
    Returns:
        Dictionary with summary statistics
    """
    import json
    
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        summary = {
            'total_reads_before': data['summary']['before_filtering']['total_reads'],
            'total_reads_after': data['summary']['after_filtering']['total_reads'],
            'reads_passed': data['filtering_result']['passed_filter_reads'],
            'reads_failed': data['filtering_result']['failed_filter_reads'],
            'q30_rate_before': data['summary']['before_filtering']['q30_rate'],
            'q30_rate_after': data['summary']['after_filtering']['q30_rate'],
            'adapter_trimmed_reads': data['adapter_cutting']['adapter_trimmed_reads'],
        }
        
        # Calculate percentage
        if summary['total_reads_before'] > 0:
            summary['percent_passed'] = (summary['reads_passed'] / summary['total_reads_before']) * 100
        else:
            summary['percent_passed'] = 0
        
        return summary
    
    except Exception as e:
        print(f"Warning: Could not parse fastp report {json_path}: {e}")
        return {}
