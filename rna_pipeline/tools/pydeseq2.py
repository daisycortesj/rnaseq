"""
PyDESeq2 wrapper — differential expression analysis.

Unlike the other tools in this package, PyDESeq2 is a Python library
(not a command-line tool). This module provides helper functions to
run the 3-step DESeq2 workflow programmatically.
"""

from typing import Optional
import os


def check_pydeseq2() -> bool:
    """Return True if PyDESeq2 is importable."""
    try:
        import pydeseq2
        return True
    except ImportError:
        return False


def get_pydeseq2_version() -> str:
    """Return the installed PyDESeq2 version string, or 'not installed'."""
    try:
        import pydeseq2
        return getattr(pydeseq2, "__version__", "unknown")
    except ImportError:
        return "not installed"


def run_deseq2(
    count_matrix_path: str,
    metadata_path: str,
    output_dir: str,
    design_factor: str = "condition",
    ref_level: Optional[str] = None,
    min_counts: int = 10,
):
    """
    Run PyDESeq2 differential expression on a count matrix.

    Parameters
    ----------
    count_matrix_path : str
        Path to gene_count_matrix.tsv (genes x samples).
    metadata_path : str
        Path to sample_metadata.tsv with a design factor column.
    output_dir : str
        Directory for results output.
    design_factor : str
        Column name in metadata to test (default "condition").
    ref_level : str or None
        Reference level for the design factor. If None, PyDESeq2
        chooses alphabetically.
    min_counts : int
        Minimum total counts across samples to keep a gene.

    Returns
    -------
    DeseqStats
        The fitted stats object (contains results_df with log2FC, padj, etc.).
    """
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = pd.read_csv(count_matrix_path, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)

    counts = counts.loc[:, metadata.index]

    keep = counts.sum(axis=1) >= min_counts
    counts = counts.loc[keep]
    print(f"[PyDESeq2] {keep.sum()} genes pass min_counts={min_counts} filter "
          f"(dropped {(~keep).sum()})")

    dds = DeseqDataSet(
        counts=counts.T,
        metadata=metadata,
        design_factors=design_factor,
        ref_level=[design_factor, ref_level] if ref_level else None,
    )

    dds.deseq2()

    stat = DeseqStats(dds)
    stat.summary()

    os.makedirs(output_dir, exist_ok=True)
    results = stat.results_df
    out_path = os.path.join(output_dir, "pydeseq2_results_UNFILTERED.tsv")
    results.to_csv(out_path, sep="\t")
    print(f"[PyDESeq2] Results written to {out_path}")

    return stat


def filter_results(
    results_path: str,
    output_path: str,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
):
    """
    Filter a DESeq2 results TSV by significance and fold-change.

    Parameters
    ----------
    results_path : str
        Path to the unfiltered results TSV.
    output_path : str
        Path for the filtered output TSV.
    padj_cutoff : float
        Maximum adjusted p-value (FDR).
    lfc_cutoff : float
        Minimum absolute log2 fold change.
    direction : str
        "up" for upregulated only, "down" for downregulated, "both" for either.

    Returns
    -------
    pandas.DataFrame
        The filtered results.
    """
    import pandas as pd

    df = pd.read_csv(results_path, sep="\t", index_col=0)

    mask = df["padj"] < padj_cutoff

    if direction == "up":
        mask &= df["log2FoldChange"] > lfc_cutoff
    elif direction == "down":
        mask &= df["log2FoldChange"] < -lfc_cutoff
    else:
        mask &= df["log2FoldChange"].abs() > lfc_cutoff

    filtered = df.loc[mask].sort_values("padj")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    filtered.to_csv(output_path, sep="\t")
    print(f"[PyDESeq2] {len(filtered)} genes pass filters "
          f"(padj<{padj_cutoff}, |log2FC|>{lfc_cutoff}, direction={direction})")

    return filtered
