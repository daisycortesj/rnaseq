
from .logging_setup import setup_logging
from .utils.io_utils import file_nonempty, ensure_dir
from .utils.sys_utils import which
from .tools.star import build_star_index_cmd, build_star_align_cmd
from .tools.trinity import build_trinity_cmd
from .tools.qc import (check_qc_tools, install_qc_tools, find_fastq_files,
                       build_fastqc_cmd, build_multiqc_cmd, summarize_fastq_files)
from .tools.featurecounts import check_featurecounts
from .tools.blast import check_blast
from .tools.hmmer import check_hmmer
from .tools.prosite import check_prosite
from .tools.pydeseq2 import check_pydeseq2, get_pydeseq2_version
from .runners.local import run as run_local
import os
import subprocess


def verify_tool(name, check_fn, log, version_cmd=None):
    """
    Check if a tool is available. Log result and optionally print its version.
    Returns True if found, False if missing.
    """
    ok = check_fn()
    if ok:
        version = ""
        if version_cmd:
            try:
                result = subprocess.run(
                    version_cmd, capture_output=True, text=True, check=False
                )
                raw = (result.stdout + result.stderr).strip().split("\n")[0]
                version = f"  ({raw})"
            except Exception:
                pass
        log.info(f"  OK   {name}{version}")
    else:
        log.warning(f"  MISS {name}  <-- not found")
    return ok


def run_check_tools(log):
    """Verify every tool the pipeline can use and report status."""
    log.info("=" * 55)
    log.info("Pipeline Tool Check")
    log.info("=" * 55)

    results = {}

    log.info("")
    log.info("-- Upstream (QC / Alignment / Assembly) --")
    results["FastQC/MultiQC"] = check_qc_tools()[0]
    if results["FastQC/MultiQC"]:
        log.info("  OK   FastQC / MultiQC")
    else:
        log.warning("  MISS FastQC / MultiQC")
    results["STAR"]     = verify_tool("STAR",     lambda: which("STAR"),    log, ["STAR", "--version"])
    results["samtools"] = verify_tool("samtools",  lambda: which("samtools"), log, ["samtools", "--version"])
    results["Trinity"]  = verify_tool("Trinity",   lambda: which("Trinity"), log)

    log.info("")
    log.info("-- Counting & Differential Expression --")
    results["featureCounts"] = verify_tool("featureCounts", check_featurecounts, log)
    results["PyDESeq2"]      = verify_tool(
        f"PyDESeq2 ({get_pydeseq2_version()})", check_pydeseq2, log
    )

    log.info("")
    log.info("-- Annotation & Domains --")
    results["BLAST"]   = verify_tool("blastp / blastx", check_blast,  log, ["blastp", "-version"])
    results["HMMER"]   = verify_tool("hmmscan (HMMER)", check_hmmer,  log)
    results["PROSITE"] = verify_tool("patmatmotifs (EMBOSS)", check_prosite, log)

    log.info("")
    log.info("=" * 55)
    n_ok   = sum(1 for v in results.values() if v)
    n_total = len(results)
    if n_ok == n_total:
        log.info(f"All {n_total} tools found. Pipeline is ready.")
    else:
        n_miss = n_total - n_ok
        log.warning(f"{n_ok}/{n_total} tools found. {n_miss} missing.")
        log.warning("Install missing tools: conda env update -f environment.yml")
    log.info("=" * 55)
    return n_ok == n_total


def run_workflow(args):
    log = setup_logging()

    if args.mode == "check-tools":
        ok = run_check_tools(log)
        raise SystemExit(0 if ok else 1)

    # Load config if --genome is specified
    config = None
    if args.genome:
        try:
            from .config import PipelineConfig
            config = PipelineConfig(args.config)
            log.info(f"Loaded configuration for genome: {args.genome}")
            
            # Auto-populate genome-specific parameters
            genome_cfg = config.get_genome_config(args.genome)
            
            # Set genome index if not explicitly provided
            if not args.genome_index and args.mode in ["align", "auto"]:
                args.genome_index = config.get_genome_index(args.genome)
                log.info(f"Using genome index: {args.genome_index}")
            
            # Set quant_mode based on genome config if not explicitly set
            if args.quant_mode is None:
                args.quant_mode = genome_cfg.get('quant_mode', True)
                log.info(f"Gene quantification: {'enabled' if args.quant_mode else 'disabled'}")
        
        except ImportError:
            log.warning("PyYAML not installed. Install with: pip install pyyaml")
            log.warning("Continuing without config support...")
        except Exception as e:
            log.error(f"Error loading config: {e}")
            raise SystemExit(1)

    # Determine mode
    mode = determine_mode(args, log)
    log.info(f"Pipeline mode: {mode}")

    # Check required tools based on mode
    check_tools(mode, log)

    # Route to appropriate workflow
    if mode == "index":
        run_index_workflow(args, log)
    elif mode == "align":
        run_align_workflow(args, log)
    elif mode == "trinity":
        run_trinity_workflow(args, log)
    elif mode == "qc":
        run_qc_workflow(args, log)


def determine_mode(args, log):
    """Determine which mode to run based on arguments."""
    if args.mode != "auto":
        return args.mode
    
    # Auto-detect mode from arguments
    have_fasta = file_nonempty(args.fasta) if args.fasta else False
    have_genome_index = args.genome_index and os.path.isdir(args.genome_index)
    have_reads = (args.reads_left and args.reads_right)
    
    if have_genome_index and have_reads:
        return "align"
    elif have_fasta:
        return "index"
    elif have_reads:
        return "trinity"
    else:
        raise SystemExit("Cannot determine mode. Provide: --fasta (for index), --genome-index + --reads (for align), or just --reads (for Trinity)")


def check_tools(mode, log):
    """Check if required tools are available."""
    missing = []
    
    if mode in ["index", "align"]:
        if not which("STAR"):
            missing.append("STAR")
    
    if mode == "trinity":
        if not which("Trinity"):
            missing.append("Trinity")
    
    if mode == "qc":
        # QC tools should be in environment.yml (primary installation method)
        # But we provide auto-install as a fallback/safety net for robustness
        tools_available, missing_qc = check_qc_tools()
        if not tools_available:
            log.warning(f"QC tools not found: {', '.join(missing_qc)}")
            log.info("Note: These tools should be in your conda environment.")
            log.info("Recommended: conda env update -f environment.yml")
            log.info("Attempting automatic installation as fallback...")
            rc = install_qc_tools()
            if rc != 0:
                log.error("Failed to install QC tools. Please install manually:")
                log.error("  conda env update -f environment.yml")
                log.error("  OR: conda install -c bioconda fastqc multiqc")
                raise SystemExit(1)
            log.info("✓ QC tools installed successfully!")
            log.info("✓ For future runs, update your environment to avoid this delay.")
    
    if missing:
        raise SystemExit("Missing tools: " + ", ".join(missing))
    
    if mode != "qc":
        log.info(f"Tools available: {', '.join(['STAR'] if mode in ['index', 'align'] else ['Trinity'])}")


def run_index_workflow(args, log):
    """Build STAR genome index."""
    have_fasta = file_nonempty(args.fasta) if args.fasta else False
    have_gtf   = file_nonempty(args.gtf)   if args.gtf   else False
    
    if not have_fasta:
        raise SystemExit("Index mode requires --fasta")
    
    if have_gtf:
        log.info("Building STAR index with gene annotations")
    else:
        log.info("Building STAR index WITHOUT gene annotations (no GTF)")
    
    ensure_dir(args.outdir)
    cmd = build_star_index_cmd(args.fasta, args.gtf, args.outdir, args.threads, args.readlen)
    rc = run_local(cmd, dry=args.dry)
    
    if rc != 0:
        raise SystemExit(rc)
    
    log.info("Done: STAR index at %s", args.outdir)


def run_align_workflow(args, log):
    """Align reads using STAR."""
    if not args.genome_index:
        raise SystemExit("Align mode requires --genome-index")
    
    if not os.path.isdir(args.genome_index):
        raise SystemExit(f"Genome index not found: {args.genome_index}")
    
    if not args.reads_left:
        raise SystemExit("Align mode requires --reads-left")
    
    if not file_nonempty(args.reads_left):
        raise SystemExit(f"Read file not found or empty: {args.reads_left}")
    
    # Check if reads_right exists (for paired-end)
    reads_right = ""
    if args.reads_right:
        if file_nonempty(args.reads_right):
            reads_right = args.reads_right
        else:
            log.warning(f"Read file not found: {args.reads_right}, treating as single-end")
    
    # Determine output prefix
    if args.sample_name:
        sample_name = args.sample_name
    else:
        # Extract sample name from reads_left filename
        sample_name = os.path.basename(args.reads_left).split('.')[0].replace('_R1', '').replace('_1', '')
        log.info(f"Auto-detected sample name: {sample_name}")
    
    ensure_dir(args.outdir)
    outprefix = os.path.join(args.outdir, sample_name + "_")
    
    # Set default quant_mode if still None
    quant_mode = args.quant_mode if args.quant_mode is not None else True
    
    log.info(f"Aligning reads from: {args.reads_left}" + (f" + {reads_right}" if reads_right else " (single-end)"))
    log.info(f"Output prefix: {outprefix}")
    
    cmd = build_star_align_cmd(
        args.genome_index, 
        args.reads_left, 
        reads_right, 
        outprefix, 
        args.threads,
        quant_mode
    )
    
    rc = run_local(cmd, dry=args.dry)
    
    if rc != 0:
        raise SystemExit(rc)
    
    log.info("Done: Alignment at %s", args.outdir)
    log.info("Main output: %sAligned.sortedByCoord.out.bam", outprefix)
    if quant_mode:
        log.info("Gene counts: %sReadsPerGene.out.tab", outprefix)


def run_trinity_workflow(args, log):
    """Run Trinity de novo assembly."""
    have_reads = (args.reads_left and args.reads_right)
    
    if not have_reads:
        raise SystemExit("Trinity mode requires --reads-left and --reads-right")
    
    # Check for resume mode
    resume = getattr(args, 'resume', False)
    
    log.info("Running Trinity de novo assembly")
    log.info(f"Left reads: {args.reads_left}")
    log.info(f"Right reads: {args.reads_right}")
    log.info(f"Output directory: {args.outdir}")
    log.info(f"Threads: {args.threads}, Memory: {args.mem_gb}G")
    if resume:
        log.info("RESUME MODE: Trinity will auto-detect checkpoints and continue")
    
    # Check if output directory exists
    if os.path.exists(args.outdir):
        trinity_marker = os.path.join(args.outdir, "trinity_out_dir")
        if os.path.exists(trinity_marker):
            if resume:
                log.info(f"Found previous Trinity run - resuming from checkpoint")
            else:
                log.warning(f"Output directory contains a previous Trinity run: {args.outdir}")
                log.warning("Use --resume to continue from checkpoint, or remove the directory: rm -rf %s", args.outdir)
    
    ensure_dir(args.outdir)
    
    cmd = build_trinity_cmd(args.reads_left, args.reads_right, args.outdir, args.threads, args.mem_gb, resume=resume)
    rc = run_local(cmd, dry=args.dry)
    
    if rc != 0:
        log.error("Trinity failed with exit code %d", rc)
        log.error("Common issues:")
        log.error("  1. Output directory already exists from previous run")
        log.error("  2. Input files not found or inaccessible")
        log.error("  3. Insufficient memory or disk space")
        log.error("  4. File format issues (ensure files are gzipped FASTQ)")
        raise SystemExit(rc)
    
    log.info("Done: Trinity assembly at %s", args.outdir)


def run_qc_workflow(args, log):
    """Run FastQC and MultiQC quality control."""
    # Determine input: either specific FASTQ files or search directory
    fastq_files = []
    
    if args.fastq_dir:
        # Search directory for all FASTQ files
        log.info(f"Searching for FASTQ files in: {args.fastq_dir}")
        fastq_files = find_fastq_files(args.fastq_dir)
        
        if not fastq_files:
            raise SystemExit(f"No FASTQ files found in: {args.fastq_dir}")
        
        log.info(f"Found {len(fastq_files)} FASTQ files")
    
    elif args.reads_left or args.reads_right:
        # Use specific files provided
        if args.reads_left and file_nonempty(args.reads_left):
            fastq_files.append(args.reads_left)
        if args.reads_right and file_nonempty(args.reads_right):
            fastq_files.append(args.reads_right)
        
        if not fastq_files:
            raise SystemExit("QC mode requires either --fastq-dir or --reads-left/--reads-right")
        
        log.info(f"Running QC on {len(fastq_files)} specified file(s)")
    
    else:
        raise SystemExit("QC mode requires either --fastq-dir or --reads-left/--reads-right")
    
    # Show summary of files
    summary = summarize_fastq_files(fastq_files)
    log.info(f"Total files: {summary['total_files']}")
    for file_info in summary['files'][:5]:  # Show first 5
        log.info(f"  - {file_info['name']} ({file_info['size_mb']} MB)")
    if summary['total_files'] > 5:
        log.info(f"  ... and {summary['total_files'] - 5} more files")
    
    # Create output directories
    ensure_dir(args.outdir)
    fastqc_dir = os.path.join(args.outdir, "fastqc_results")
    ensure_dir(fastqc_dir)
    
    # Run FastQC
    log.info("=" * 60)
    log.info("Step 1: Running FastQC...")
    log.info("=" * 60)
    
    fastqc_cmd = build_fastqc_cmd(fastq_files, fastqc_dir, args.threads)
    rc = run_local(fastqc_cmd, dry=args.dry)
    
    if rc != 0:
        log.error("FastQC failed with exit code %d", rc)
        raise SystemExit(rc)
    
    log.info("✓ FastQC completed successfully")
    log.info(f"  Reports: {fastqc_dir}/")
    
    # Run MultiQC
    log.info("=" * 60)
    log.info("Step 2: Running MultiQC to aggregate reports...")
    log.info("=" * 60)
    
    multiqc_cmd = build_multiqc_cmd(fastqc_dir, args.outdir, 
                                    title="RNA-seq Quality Control Report")
    rc = run_local(multiqc_cmd, dry=args.dry)
    
    if rc != 0:
        log.error("MultiQC failed with exit code %d", rc)
        raise SystemExit(rc)
    
    log.info("✓ MultiQC completed successfully")
    
    # Final summary
    log.info("=" * 60)
    log.info("Quality Control Complete!")
    log.info("=" * 60)
    log.info(f"Individual FastQC reports: {fastqc_dir}/")
    log.info(f"MultiQC summary report:    {os.path.join(args.outdir, 'multiqc_report.html')}")
    log.info("")
    log.info("NEXT STEPS:")
    log.info("1. Download and open: multiqc_report.html")
    log.info("2. Check 'General Statistics' and 'Sequence Quality' sections")
    log.info("3. Look for quality scores >28 (green zones)")
    log.info("4. Check for adapter contamination (<5% is good)")
    log.info("5. If quality is good → proceed with alignment/assembly")
    log.info("6. If quality issues found → consider trimming with fastp")
    log.info("=" * 60)
