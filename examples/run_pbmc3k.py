#!/usr/bin/env python3
"""
examples/run_pbmc3k.py
======================

The worked example, on real data.

Downloads the 3k PBMC dataset that 10x Genomics published in 2016, runs the
same QC pipeline the command line tool runs, and writes the report next to
this file.

    python examples/run_pbmc3k.py

scanpy fetches the data the first time you run this, about 6 MB, and caches
it in ./data so later runs are offline. The default thresholds were chosen
for exactly this kind of sample, so no tuning is needed here.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite

# qc_metrics and qc_report live in the repo root, one level up from examples/.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qc_metrics import Thresholds, run_qc              # noqa: E402
from qc_report import STATUS, build_checks, overall_verdict, write_report  # noqa: E402

HERE = Path(__file__).resolve().parent


def write_10x_folder(adata, out_dir: Path) -> None:
    """
    Write an in-memory AnnData back out as a 10x matrix folder.

    scanpy hands us pbmc3k as an AnnData, but run_qc reads a CellRanger
    folder, which is what a real user actually has. Rather than add a second
    entry point to the library for this one example, we write the object back
    out to the format the tool already reads. The example then exercises the
    exact code path the command line uses, loader included.

    The folder is temporary and thrown away as soon as the report is written.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(adata.obs_names) + "\n")

    # Keep the real Ensembl IDs where the dataset carries them, so the folder
    # is a faithful copy rather than something that merely parses.
    ids = (adata.var["gene_ids"].astype(str)
           if "gene_ids" in adata.var else adata.var_names.astype(str))
    with gzip.open(out_dir / "features.tsv.gz", "wt") as fh:
        for gene_id, symbol in zip(ids, adata.var_names):
            fh.write(f"{gene_id}\t{symbol}\tGene Expression\n")

    # 10x stores genes as rows and cells as columns, the transpose of AnnData.
    # The values are UMI counts, so they are integers even though scanpy may
    # hand them over as floats.
    counts = adata.X
    counts = counts.T.tocoo() if sparse.issparse(counts) else sparse.coo_matrix(counts.T)
    counts = counts.astype(np.int32)
    with gzip.open(out_dir / "matrix.mtx.gz", "wb") as fh:
        mmwrite(fh, counts, field="integer")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Run the QC report on the public 10x pbmc3k dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--out", default=str(HERE / "pbmc3k_report.html"),
                   help="Where to write the HTML report.")
    args = p.parse_args(argv)

    print("Fetching pbmc3k from scanpy (cached in ./data after the first run) ...",
          flush=True)
    adata = sc.datasets.pbmc3k()
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        matrix_dir = Path(tmp) / "pbmc3k_filtered_gene_bc_matrices"
        write_10x_folder(adata, matrix_dir)
        result = run_qc(
            matrix_dir,
            sample_name="pbmc3k (10x Genomics, healthy donor PBMCs)",
            thresholds=Thresholds(),
        )
        # The temporary folder is an implementation detail, so do not print it
        # into the report as though it were where the data came from.
        result.matrix_dir = "scanpy.datasets.pbmc3k()"
        out = Path(args.out).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        write_report(result, str(out))

    checks = build_checks(result)
    status, sentence = overall_verdict(checks)
    print()
    print(f"  {STATUS[status]['icon']} {STATUS[status]['label'].upper()}")
    print(f"  {sentence}")
    print()
    for c in checks:
        print(f"    [{STATUS[c.status]['icon']}] {c.metric}: {c.headline}")
    print()
    print(f"Report written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
