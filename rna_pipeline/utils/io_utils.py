
import os

def file_nonempty(p: str) -> bool:
    return bool(p) and os.path.isfile(p) and os.path.getsize(p) > 0

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)
