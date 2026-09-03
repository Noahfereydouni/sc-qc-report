"""
qc_metrics.py
=============

Everything that *reads data and computes numbers* lives here.
No plotting, no HTML - just the biology-facing measurements.

The six metrics this module produces:

  1. UMI counts per cell          -> adata.obs["total_counts"]
  2. Genes detected per cell      -> adata.obs["n_genes_by_counts"]
  3. Percent mitochondrial reads  -> adata.obs["pct_counts_mt"]
  4. Percent ribosomal reads      -> adata.obs["pct_counts_ribo"]
  5. Doublet score (Scrublet)     -> adata.obs["doublet_score"]
  6. Cell / gene counts before and after filtering -> QCResult.counts
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


# --------------------------------------------------------------------------
# Filtering thresholds
# --------------------------------------------------------------------------
# These are ordinary 10x scRNA-seq starting points, not universal truths.
# Change them from the command line (see run_qc.py --help) if your tissue
# behaves differently - for example heart and kidney tolerate more
# mitochondrial signal than PBMCs do.

@dataclass
class Thresholds:
    min_genes_per_cell: int = 200     # a cell must detect at least this many genes
    max_pct_mito: float = 15.0        # above this, the cell is likely dying/lysed
    min_cells_per_gene: int = 3       # a gene must appear in at least this many cells
    drop_predicted_doublets: bool = True   # remove barcodes Scrublet calls doublets


# --------------------------------------------------------------------------
# The object we hand to the report writer
# --------------------------------------------------------------------------

@dataclass
class QCResult:
    """Everything the report needs, and nothing it doesn't."""

    sample_name: str
    matrix_dir: str
    thresholds: Thresholds

    # Per-cell metrics for every barcode in the input matrix (metrics 1-5).
    # One row per cell; this is what the violins/histograms are drawn from.
    per_cell: pd.DataFrame

    # Metric 6: the before/after table.
    counts: dict = field(default_factory=dict)

    # True if Scrublet ran successfully. If it didn't (too few cells, for
    # instance) the doublet column is all-NaN and the report says so.
    scrublet_ran: bool = True
    scrublet_note: str = ""


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_matrix(matrix_dir: str | Path) -> sc.AnnData:
    """
    Read a 10x 'filtered_feature_bc_matrix' folder.

    That folder is expected to hold the three files CellRanger writes:
    barcodes.tsv.gz, features.tsv.gz (or genes.tsv.gz) and matrix.mtx.gz.
    """
    matrix_dir = Path(matrix_dir)
    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"Not a folder: {matrix_dir}")

    expected = {"matrix.mtx", "matrix.mtx.gz"}
    if not any((matrix_dir / name).exists() for name in expected):
        raise FileNotFoundError(
            f"{matrix_dir} does not look like a 10x matrix folder "
            "(no matrix.mtx or matrix.mtx.gz inside)."
        )

    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    # Two cells can share a gene symbol in some references; make names unique
    # so downstream indexing is unambiguous.
    adata.var_names_make_unique()
    return adata


# --------------------------------------------------------------------------
# Metrics 1-4: the standard per-cell QC numbers
# --------------------------------------------------------------------------

def annotate_gene_families(adata: sc.AnnData) -> None:
    """
    Flag mitochondrial and ribosomal genes by name.

    Human symbols are MT-* / RPS* / RPL*; mouse are mt-* / Rps* / Rpl*.
    Matching case-insensitively covers both without asking the user which
    species they ran.
    """
    names = adata.var_names.str.upper()
    adata.var["mt"] = names.str.startswith("MT-")
    adata.var["ribo"] = names.str.startswith(("RPS", "RPL"))


def compute_per_cell_metrics(adata: sc.AnnData) -> None:
    """Fill in total_counts, n_genes_by_counts, pct_counts_mt, pct_counts_ribo."""
    annotate_gene_families(adata)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )


# --------------------------------------------------------------------------
# Metric 5: doublet score
# --------------------------------------------------------------------------

def compute_doublet_scores(adata: sc.AnnData) -> tuple[bool, str]:
    """
    Run Scrublet and write 'doublet_score' + 'predicted_doublet' into .obs.

    We use scanpy's built-in sc.pp.scrublet, which is the Scrublet algorithm
    (Wolock et al. 2019) maintained inside scanpy - same method, one less
    package to install.

    Scrublet needs a reasonable number of cells to build its simulated-doublet
    neighbourhood. On a very small matrix it can fail outright, so we catch
    that and let the report say "not available" instead of crashing.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.scrublet(adata, verbose=False)
        return True, ""
    except Exception as exc:  # noqa: BLE001 - any failure means "no score"
        adata.obs["doublet_score"] = np.nan
        adata.obs["predicted_doublet"] = False
        return False, f"Scrublet could not run on this matrix ({exc})."


# --------------------------------------------------------------------------
# Metric 6: filtering, and the before/after tally
# --------------------------------------------------------------------------

def apply_filters(adata: sc.AnnData, thr: Thresholds) -> tuple[sc.AnnData, dict]:
    """
    Apply the standard first-pass filters and count what survived.

    Order matters and follows the usual scanpy workflow:
      1. drop low-complexity cells (too few genes)
      2. drop high-mitochondrial cells
      3. drop predicted doublets
      4. drop genes seen in almost no cells

    Returns the filtered object plus the before/after numbers (metric 6).
    """
    n_cells_before, n_genes_before = adata.n_obs, adata.n_vars

    keep = adata.obs["n_genes_by_counts"] >= thr.min_genes_per_cell
    dropped_low_genes = int((~keep).sum())

    keep_mito = adata.obs["pct_counts_mt"] <= thr.max_pct_mito
    dropped_high_mito = int((keep & ~keep_mito).sum())
    keep = keep & keep_mito

    dropped_doublets = 0
    if thr.drop_predicted_doublets and "predicted_doublet" in adata.obs:
        keep_singlet = ~adata.obs["predicted_doublet"].astype(bool)
        dropped_doublets = int((keep & ~keep_singlet).sum())
        keep = keep & keep_singlet

    filtered = adata[keep.values].copy()

    # Gene-level filter is applied only after cells are chosen, so "genes
    # after" means "genes still detected in the cells we kept".
    if filtered.n_obs > 0:
        sc.pp.filter_genes(filtered, min_cells=thr.min_cells_per_gene)

    counts = {
        "cells_before": int(n_cells_before),
        "cells_after": int(filtered.n_obs),
        "genes_before": int(n_genes_before),
        "genes_after": int(filtered.n_vars),
        "dropped_low_genes": dropped_low_genes,
        "dropped_high_mito": dropped_high_mito,
        "dropped_doublets": dropped_doublets,
    }
    counts["cells_removed"] = counts["cells_before"] - counts["cells_after"]
    counts["pct_cells_kept"] = (
        100.0 * counts["cells_after"] / counts["cells_before"]
        if counts["cells_before"]
        else 0.0
    )
    return filtered, counts


# --------------------------------------------------------------------------
# The one function run_qc.py actually calls
# --------------------------------------------------------------------------

def run_qc(
    matrix_dir: str | Path,
    sample_name: str | None = None,
    thresholds: Thresholds | None = None,
) -> QCResult:
    """Load a 10x folder, compute all six metrics, return them in one object."""
    thr = thresholds or Thresholds()
    matrix_dir = Path(matrix_dir)
    sample_name = sample_name or matrix_dir.name

    adata = load_matrix(matrix_dir)
    compute_per_cell_metrics(adata)
    scrublet_ran, scrublet_note = compute_doublet_scores(adata)
    _, counts = apply_filters(adata, thr)

    per_cell = pd.DataFrame(
        {
            "umi_counts": adata.obs["total_counts"].to_numpy(dtype=float),
            "genes_detected": adata.obs["n_genes_by_counts"].to_numpy(dtype=float),
            "pct_mito": adata.obs["pct_counts_mt"].to_numpy(dtype=float),
            "pct_ribo": adata.obs["pct_counts_ribo"].to_numpy(dtype=float),
            "doublet_score": adata.obs["doublet_score"].to_numpy(dtype=float),
        },
        index=adata.obs_names,
    )
    # Fraction of barcodes Scrublet flagged - used by the verdict.
    counts["pct_predicted_doublets"] = (
        100.0 * float(np.asarray(adata.obs["predicted_doublet"]).astype(bool).mean())
        if scrublet_ran
        else float("nan")
    )

    return QCResult(
        sample_name=sample_name,
        matrix_dir=str(matrix_dir),
        thresholds=thr,
        per_cell=per_cell,
        counts=counts,
        scrublet_ran=scrublet_ran,
        scrublet_note=scrublet_note,
    )
