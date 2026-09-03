#!/usr/bin/env python3
"""
make_test_matrix.py
===================

Writes a fake 10x 'filtered_feature_bc_matrix' folder so you can try the
report without downloading a real dataset.

    python make_test_matrix.py test_data/good_sample
    python make_test_matrix.py test_data/poor_sample --quality poor

The three files it writes (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
are exactly the CellRanger v3 layout, so scanpy reads the folder the same way
it reads a real one.

This is a plausible-looking simulation, not a biologically faithful one - it
exists to exercise the report, not to benchmark anything.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

# The 13 protein-coding mitochondrial genes, named as they appear in the
# 10x human reference - the report finds them by the "MT-" prefix.
MT_GENES = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
            "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6",
            "MT-CYB"]

# Two quality settings, so you can see the report reach both verdicts.
PROFILES = {
    "good": dict(median_umis=4200, mito_mean=0.05, doublet_rate=0.05,
                 debris_frac=0.05),
    "poor": dict(median_umis=750,  mito_mean=0.24, doublet_rate=0.18,
                 debris_frac=0.22),
}


def gene_names(n_other: int, n_ribo: int = 60) -> list[str]:
    """Mitochondrial + ribosomal + filler genes, in that recognisable style."""
    ribo = ([f"RPS{i}" for i in range(1, n_ribo // 2 + 1)] +
            [f"RPL{i}" for i in range(1, n_ribo - n_ribo // 2 + 1)])
    other = [f"GENE{i:05d}" for i in range(n_other)]
    return MT_GENES + ribo + other


def simulate(n_cells: int, n_genes: int, profile: dict,
             seed: int = 0) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """
    Build a genes x cells count matrix with the structure QC is meant to find.

    Deliberately baked in:
      * three cell types, each over-expressing its own marker genes
      * a debris population - low depth, high mitochondrial fraction
      * doublets, made by literally adding two cells' counts together
    """
    rng = np.random.default_rng(seed)

    names = gene_names(n_other=n_genes, n_ribo=60)
    n_mt, n_ribo = len(MT_GENES), 60
    n_nuclear = len(names) - n_mt          # ribosomal + filler
    total_genes = len(names)

    # --- expression profile per cell type ---------------------------------
    base = rng.lognormal(mean=0.0, sigma=1.4, size=n_nuclear)
    base[:n_ribo] *= 25          # ribosomal genes are genuinely abundant

    # Real references carry thousands of genes that are barely expressed in a
    # given tissue. Without them the min-cells gene filter has nothing to do,
    # which makes the before/after gene panel look broken.
    barely_expressed = rng.random(n_nuclear) < 0.20
    barely_expressed[:n_ribo] = False
    base[barely_expressed] *= 1e-4

    n_types = 3
    profiles = []
    for t in range(n_types):
        p = base.copy()
        markers = rng.choice(n_nuclear - n_ribo, size=40, replace=False) + n_ribo
        p[markers] *= 12
        profiles.append(p / p.sum())
    cell_type = rng.integers(0, n_types, size=n_cells)

    # --- sequencing depth and mitochondrial fraction per cell -------------
    depth = rng.lognormal(mean=np.log(profile["median_umis"]), sigma=0.45,
                          size=n_cells)

    # Debris/dying cells: shallow and mitochondria-heavy.
    is_debris = rng.random(n_cells) < profile["debris_frac"]
    depth[is_debris] /= 7

    # Beta gives a fraction in (0, 1) with a realistic right-hand tail.
    conc = 12.0
    m = profile["mito_mean"]
    mito_frac = rng.beta(m * conc, (1 - m) * conc, size=n_cells)
    mito_frac[is_debris] = np.clip(mito_frac[is_debris] * 3.0, 0, 0.9)

    mt_p = rng.dirichlet(np.ones(n_mt) * 2.0)

    # --- draw the counts ---------------------------------------------------
    counts = np.zeros((total_genes, n_cells), dtype=np.int32)
    for c in range(n_cells):
        n_nuc = max(int(depth[c] * (1 - mito_frac[c])), 1)
        n_mito = int(depth[c] * mito_frac[c])
        counts[n_mt:, c] = rng.multinomial(n_nuc, profiles[cell_type[c]])
        if n_mito > 0:
            counts[:n_mt, c] = rng.multinomial(n_mito, mt_p)

    # --- doublets: two barcodes' worth of RNA under one barcode -----------
    n_doublets = int(n_cells * profile["doublet_rate"])
    doublet_idx = rng.choice(n_cells, size=n_doublets, replace=False)
    for c in doublet_idx:
        partner = rng.integers(0, n_cells)
        counts[:, c] += counts[:, partner]

    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n_cells)]
    return sparse.csr_matrix(counts), names, barcodes


def write_10x(out_dir: Path, matrix: sparse.csr_matrix,
              genes: list[str], barcodes: list[str]) -> None:
    """Write the three gzipped files CellRanger produces."""
    out_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(barcodes) + "\n")

    # features.tsv.gz is three columns: ID, symbol, feature type.
    with gzip.open(out_dir / "features.tsv.gz", "wt") as fh:
        for i, symbol in enumerate(genes):
            fh.write(f"ENSGSYNTH{i:08d}\t{symbol}\tGene Expression\n")

    # 10x stores the matrix genes-by-cells; scanpy transposes it on read.
    with gzip.open(out_dir / "matrix.mtx.gz", "wb") as fh:
        mmwrite(fh, matrix.tocoo(), field="integer")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Write a synthetic 10x filtered matrix folder for testing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("out_dir", help="Folder to create.")
    p.add_argument("--cells", type=int, default=2000, help="Number of barcodes.")
    p.add_argument("--genes", type=int, default=1500,
                   help="Number of filler genes (mito and ribosomal are added "
                        "on top).")
    p.add_argument("--quality", choices=sorted(PROFILES), default="good",
                   help="'good' makes a healthy-looking sample; 'poor' makes one "
                        "the report should complain about.")
    p.add_argument("--seed", type=int, default=0, help="Random seed.")
    args = p.parse_args(argv)

    matrix, genes, barcodes = simulate(
        args.cells, args.genes, PROFILES[args.quality], seed=args.seed)
    out = Path(args.out_dir)
    write_10x(out, matrix, genes, barcodes)

    print(f"Wrote a '{args.quality}' synthetic sample to {out.resolve()}")
    print(f"  {len(barcodes):,} barcodes x {len(genes):,} genes")
    print(f"  files: {', '.join(sorted(f.name for f in out.iterdir()))}")
    print()
    print("Next:")
    print(f"  python run_qc.py {out} -o qc_report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
