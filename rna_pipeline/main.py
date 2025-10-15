
from .logging_setup import setup_logging
from .utils.io_utils import file_nonempty, ensure_dir
from .utils.sys_utils import which
from .tools.star import build_star_index_cmd
from .tools.trinity import build_trinity_cmd
from .runners.local import run as run_local

def run_workflow(args):
    log = setup_logging()

    # Check for both tools since we may run both
    missing = []
    if not which("STAR"): missing.append("STAR")
    if not which("Trinity"): missing.append("Trinity")
    if missing:
        raise SystemExit("Missing tools: " + ", ".join(missing))
    log.info("Tools available: STAR + Trinity")

    have_fasta = file_nonempty(args.fasta) if args.fasta else False
    have_gtf   = file_nonempty(args.gtf)   if args.gtf   else False
    have_reads = (args.reads_left and args.reads_right)

    # Run STAR if genome reference is provided (GTF is optional)
    if have_fasta:
        if have_gtf:
            log.info("Step 1: Building STAR index with gene annotations")
        else:
            log.info("Step 1: Building STAR index WITHOUT gene annotations (no GTF)")
        ensure_dir(args.outdir)
        cmd = build_star_index_cmd(args.fasta, args.gtf, args.outdir, args.threads, args.readlen)
        rc = run_local(cmd, dry=args.dry)
        if rc != 0:
            raise SystemExit(rc)
        log.info("Done: STAR index at %s", args.outdir)
    else:
        log.info("Skipping STAR: No fasta file provided")

    # Run Trinity if reads are provided
    if have_reads:
        log.info("Step 2: Running Trinity de novo assembly")
        trinity_outdir = args.outdir + "_trinity" if (have_fasta and have_gtf) else args.outdir
        ensure_dir(trinity_outdir)
        cmd = build_trinity_cmd(args.reads_left, args.reads_right, trinity_outdir, args.threads, args.mem_gb)
        rc = run_local(cmd, dry=args.dry)
        if rc != 0:
            raise SystemExit(rc)
        log.info("Done: Trinity assembly at %s", trinity_outdir)
    else:
        log.info("Skipping Trinity: No reads (--reads-left/--reads-right) provided")

    if not have_fasta and not have_reads:
        raise SystemExit("No inputs provided! Need either fasta (for STAR) or reads (for Trinity) or both.")
