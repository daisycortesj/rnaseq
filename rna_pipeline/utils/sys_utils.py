
import shutil

def which(tool: str) -> bool:
    """Return True if `tool` is on PATH."""
    return shutil.which(tool) is not None
