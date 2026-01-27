
from .logging_setup import setup_logging
from .utils.io_utils import file_nonempty, ensure_dir
from .utils.sys_utils import which
from .tools.star import build_star_index_cmd, build_star_align_cmd
from .tools.trinity import build_trinity_cmd
from .runners.local import run as run_local
import os

def run_workflow(args):
    log = setup_logging()
    
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
    
    if missing:
        raise SystemExit("Missing tools: " + ", ".join(missing))
    
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
        log.info("RESUME MODE: Will continue from checkpoint if available")
    
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
