
def build_trinity_cmd(left: str, right: str, outdir: str, threads: int, mem_gb: int):
    return [
        "Trinity",
        "--seqType", "fq",
        "--left", left,
        "--right", right,
        "--CPU", str(threads),
        "--max_memory", f"{mem_gb}G",
        "--output", outdir,
    ]
