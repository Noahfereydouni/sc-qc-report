# sc-qc-report

Point this at a 10x filtered matrix folder and it tells you, in one HTML file, whether the sample is worth analysing.

![Example QC report](examples/example_report.png)

The screenshot is a deliberately poor quality sample, so you can see what a failing verdict looks like. The file itself is [examples/example_report.html](examples/example_report.html). Everything is embedded in it, so you can email it or drop it in Slack and it opens anywhere.

## The six metrics

| Metric | Why it matters |
| --- | --- |
| UMI counts per cell | A low median means the chemistry underperformed or the cells were fragile. Small populations sink below the noise, so you cannot trust their absence. |
| Genes detected per cell | A low median, especially with a second peak near zero, means many of your barcodes are empty droplets carrying ambient RNA rather than real cells. |
| Percent mitochondrial reads | A high value means cells were dying before they reached the instrument, usually because dissociation was too harsh or too slow. Filtering afterwards does not bring those cells back. |
| Percent ribosomal reads | High ribosomal content alongside low gene counts is the signature of stressed cells whose transcriptome has collapsed onto housekeeping genes. On its own it varies widely between tissues and means little. |
| Doublet score | A high flagged fraction means too many cells were loaded. The cost turns up later as clusters that look like novel intermediate cell types but are two cells stuck together. |
| Cells and genes before and after filtering | A large drop says the problem happened before sequencing, not in the analysis. Which filter removed them says which problem: low complexity, dying cells, or overloading. |

Each metric gets a violin and a histogram. A plain language verdict at the top says whether the sample looks healthy and, if not, which metric is the concern.

## Getting started

```bash
pip install -r requirements.txt
python examples/run_pbmc3k.py
```

That runs the pipeline on the 3k PBMC dataset 10x Genomics published in 2016, which scanpy downloads and caches for you, and writes `examples/pbmc3k_report.html`.

On your own data:

```bash
python run_qc.py path/to/filtered_feature_bc_matrix --out qc_report.html
```

Tested on Python 3.11.

## A note on the thresholds

The cutoffs are the usual 10x rules of thumb, not validated numbers. 200 genes per cell and 15% mitochondrial reads are sensible for PBMCs and wrong for plenty of other tissues. Heart, kidney and tumour samples often sit well above 15% mitochondrial as normal biology. Change them with `--min-genes` and `--max-pct-mito`, which move the filtering and the verdict together. Treat the verdict as triage, not as a decision.

## Testing

`make_test_matrix.py` writes a synthetic 10x folder, so the report can be exercised with no download at all. `python make_test_matrix.py test_data/demo --quality poor` produces a sample the report should complain about, which is where the screenshot above came from.

MIT licensed. See [LICENSE](LICENSE).

Built with Claude Code. The QC logic and thresholds are mine.
