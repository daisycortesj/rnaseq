
from .logging_setup import setup_logging
from .utils.io_utils import file_nonempty, ensure_dir
from .utils.sys_utils import which
from .tools.star import build_star_index_cmd
from .tools.trinity import build_trinity_cmd
from .runners.local import run as run_local

def run_workflow(args):
    log = setup_logging()

    require_trinity = not (args.fasta and args.gtf)
    missing = []
    if not which("STAR"): missing.append("STAR")
    if require_trinity and not which("Trinity"): missing.append("Trinity")
    if missing:
        raise SystemExit("Missing tools: " + ", ".join(missing))
    log.info("Tools available: %s", "STAR" + (" + Trinity" if require_trinity else ""))

    have_fasta = file_nonempty(args.fasta) if args.fasta else False
    have_gtf   = file_nonempty(args.gtf)   if args.gtf   else False

    if have_fasta and have_gtf:
        log.info("Mode: Reference detected → STAR index")
        ensure_dir(args.outdir)
        cmd = build_star_index_cmd(args.fasta, args.gtf, args.outdir, args.threads, args.readlen)
        rc = run_local(cmd, dry=args.dry)
        if rc != 0:
            raise SystemExit(rc)
        log.info("Done: STAR index at %s", args.outdir)
    else:
        log.info("Mode: No complete reference detected → Trinity de novo")
        if not (args.reads_left and args.reads_right):
            raise SystemExit("Missing --reads-left/--reads-right for Trinity mode.")
        ensure_dir(args.outdir)
        cmd = build_trinity_cmd(args.reads_left, args.reads_right, args.outdir, args.threads, args.mem_gb)
        rc = run_local(cmd, dry=args.dry)
        if rc != 0:
            raise SystemExit(rc)
        log.info("Done: Trinity assembly at %s", args.outdir)
