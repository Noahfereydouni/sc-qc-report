#!/usr/bin/env python3
"""
run_qc.py
=========

The command you actually type.

    python run_qc.py path/to/filtered_feature_bc_matrix -o qc_report.html

Everything else is optional. The defaults are ordinary 10x starting points;
override them when your tissue calls for it, e.g. heart or kidney routinely
sit above a 15% mitochondrial cutoff:

    python run_qc.py my_matrix/ --max-pct-mito 25 --min-genes 300
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from qc_metrics import Thresholds, run_qc
from qc_report import build_checks, overall_verdict, write_report, STATUS


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Make a single self-contained HTML QC report from a 10x "
                    "filtered matrix folder.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("matrix_dir",
                   help="Folder holding barcodes.tsv.gz, features.tsv.gz and "
                        "matrix.mtx.gz (CellRanger's filtered_feature_bc_matrix).")
    p.add_argument("-o", "--out", default="qc_report.html",
                   help="Where to write the HTML report.")
    p.add_argument("-n", "--sample-name", default=None,
                   help="Name to show at the top of the report. "
                        "Defaults to the matrix folder's name.")
    p.add_argument("--min-genes", type=int, default=200,
                   help="Drop cells detecting fewer genes than this.")
    p.add_argument("--max-pct-mito", type=float, default=15.0,
                   help="Drop cells above this percent mitochondrial reads.")
    p.add_argument("--min-cells", type=int, default=3,
                   help="Drop genes detected in fewer cells than this.")
    p.add_argument("--keep-doublets", action="store_true",
                   help="Score doublets but do not filter them out.")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    thresholds = Thresholds(
        min_genes_per_cell=args.min_genes,
        max_pct_mito=args.max_pct_mito,
        min_cells_per_gene=args.min_cells,
        drop_predicted_doublets=not args.keep_doublets,
    )

    print(f"Reading {args.matrix_dir} ...", flush=True)
    try:
        result = run_qc(args.matrix_dir, args.sample_name, thresholds)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"  {result.counts['cells_before']:,} cells x "
          f"{result.counts['genes_before']:,} genes", flush=True)
    if not result.scrublet_ran:
        print(f"  note: {result.scrublet_note}", flush=True)

    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    write_report(result, str(out))

    # Echo the verdict to the terminal too, so you know whether the report is
    # worth opening before you open it.
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
