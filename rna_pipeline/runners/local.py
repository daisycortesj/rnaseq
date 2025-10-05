
import subprocess, sys

def run(cmd: list[str], dry: bool = False) -> int:
    print("[run] " + " ".join(cmd), flush=True)
    if dry:
        return 0
    try:
        subprocess.check_call(cmd)
        return 0
    except subprocess.CalledProcessError as e:
        print(f"[error] command failed with code {e.returncode}", file=sys.stderr, flush=True)
        return e.returncode
