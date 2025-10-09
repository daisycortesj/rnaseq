__version__ = "0.1.0"
__author__ = "Daisy Cortes"
__email__ = "daisycortesj@vt.edu"

# Expose a programmatic API only (
from .main import run_workflow as run  # Import the actual function name

__all__ = ["run", "__version__"]