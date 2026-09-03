# sc-qc-report

Point this at a 10x filtered matrix folder and it tells you, in one HTML file, whether the sample is worth analysing.

![Example QC report](examples/example_report.png)

The report above is the synthetic poor quality sample bundled with this repo. The file itself is [examples/example_report.html](examples/example_report.html). Everything is embedded in it, so you can email it or drop it in Slack and it opens anywhere.

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

## Running it

```bash
pip install -r requirements.txt
python run_qc.py path/to/filtered_feature_bc_matrix --out qc_report.html
```

No data to hand? `python make_test_matrix.py test_data/demo --quality poor` writes a synthetic 10x folder you can point at. Tested on Python 3.11.

## A note on the thresholds

The cutoffs are the usual 10x rules of thumb, not validated numbers. 200 genes per cell and 15% mitochondrial reads are sensible for PBMCs and wrong for plenty of other tissues. Heart, kidney and tumour samples often sit well above 15% mitochondrial as normal biology. Change them with `--min-genes` and `--max-pct-mito`, which move the filtering and the verdict together. Treat the verdict as triage, not as a decision.

MIT licensed. See [LICENSE](LICENSE).

Built with Claude Code. The QC logic and thresholds are mine.
