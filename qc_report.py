"""
qc_report.py
============

Everything that *turns numbers into something a human reads* lives here:
the plain-language verdict, the plots, and the HTML file itself.

The output is one self-contained .html - every image is embedded as base64,
the CSS is inline, nothing is fetched from the internet. You can email it,
drop it in a Slack thread, or open it on a machine with no Python at all.
"""

from __future__ import annotations

import base64
import datetime as _dt
import html
import io
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")          # render to files, never try to open a window
import matplotlib.pyplot as plt
import numpy as np

from qc_metrics import QCResult


# --------------------------------------------------------------------------
# Colour palette
# --------------------------------------------------------------------------
# One data colour (every plot shows a single distribution, so there is nothing
# to tell apart by hue) plus a fixed status palette for the verdict. Status
# colours are always paired with an icon and a word, never used alone.

SERIES = "#2a78d6"          # the one data blue
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_FAINT = "#84837d"
GRID = "#e8e7e3"
RULE = "#d6d5d1"

STATUS = {
    "good":    {"color": "#0ca30c", "icon": "✓", "label": "Looks healthy"},
    "watch":   {"color": "#fab219", "icon": "!",      "label": "Worth a look"},
    "concern": {"color": "#d03b3b", "icon": "✕", "label": "Concern"},
}
_RANK = {"good": 0, "watch": 1, "concern": 2}


# --------------------------------------------------------------------------
# The verdict
# --------------------------------------------------------------------------

@dataclass
class Check:
    """One metric's judgement, in words a biologist would actually say."""
    metric: str
    status: str          # "good" | "watch" | "concern"
    headline: str        # the number, stated plainly
    comment: str         # what it means / what to do


def _median(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def build_checks(res: QCResult) -> list[Check]:
    """
    Score each metric against ordinary 10x rules of thumb.

    These are deliberately conservative and deliberately few. They are a
    triage aid, not a substitute for looking at the plots - a "concern" on an
    unusual tissue may be perfectly normal biology.
    """
    pc, counts, thr = res.per_cell, res.counts, res.thresholds
    checks: list[Check] = []

    # 1. UMI counts per cell -------------------------------------------------
    med_umi = _median(pc["umi_counts"].to_numpy())
    if med_umi < 500:
        s, c = "concern", ("Very low. This usually means poor RNA capture - "
                           "check cell viability and library prep before analysing.")
    elif med_umi < 1000:
        s, c = "watch", ("On the low side. Usable, but rare cell types may be "
                         "hard to call.")
    else:
        s, c = "good", "A normal sequencing depth per cell."
    checks.append(Check("UMI counts per cell", s, f"median {med_umi:,.0f} UMIs", c))

    # 2. Genes detected per cell --------------------------------------------
    med_genes = _median(pc["genes_detected"].to_numpy())
    if med_genes < 250:
        s, c = "concern", ("Very low complexity. Often ambient RNA or debris "
                           "rather than intact cells.")
    elif med_genes < 500:
        s, c = "watch", ("Modest complexity. Fine for abundant cell types, thin "
                         "for subtle ones.")
    else:
        s, c = "good", "A healthy number of genes per cell."
    checks.append(Check("Genes detected per cell", s, f"median {med_genes:,.0f} genes", c))

    # 3. Percent mitochondrial ----------------------------------------------
    mito = pc["pct_mito"].to_numpy()
    med_mito = _median(mito)
    frac_over = 100.0 * float(np.mean(mito > thr.max_pct_mito)) if mito.size else 0.0
    if med_mito > thr.max_pct_mito:
        s, c = "concern", ("More than half the cells are above your cutoff - the "
                           "sample was probably stressed or dying at capture.")
    elif med_mito > 10:
        s, c = "watch", "Elevated. Some cell stress, but the bulk of cells are usable."
    else:
        s, c = "good", "Low mitochondrial signal - cells were intact."
    checks.append(Check(
        "Percent mitochondrial reads", s,
        f"median {med_mito:.1f}%  ({frac_over:.1f}% of cells above the "
        f"{thr.max_pct_mito:g}% cutoff)", c))

    # 4. Percent ribosomal ---------------------------------------------------
    # Ribosomal content varies enormously between tissues, so this never rises
    # above "worth a look" on its own.
    med_ribo = _median(pc["pct_ribo"].to_numpy())
    if med_ribo > 50:
        s, c = "watch", ("Unusually ribosome-dominated. Often goes with low "
                         "complexity - read this next to the genes-per-cell plot.")
    elif med_ribo < 5:
        s, c = "watch", ("Unusually low. Check that the reference actually "
                         "contains RPS/RPL gene symbols.")
    else:
        s, c = "good", "A typical ribosomal fraction."
    checks.append(Check("Percent ribosomal reads", s, f"median {med_ribo:.1f}%", c))

    # 5. Doublet score -------------------------------------------------------
    if not res.scrublet_ran:
        checks.append(Check("Doublet score (Scrublet)", "watch",
                            "not available", res.scrublet_note))
    else:
        pct_dbl = counts["pct_predicted_doublets"]
        med_score = _median(pc["doublet_score"].to_numpy())
        if pct_dbl > 15:
            s, c = "concern", ("A high doublet rate. Consider loading fewer cells "
                               "next time, and treat mixed-identity clusters with "
                               "suspicion.")
        elif pct_dbl > 8:
            s, c = "watch", "Somewhat above the usual 10x range, but manageable."
        else:
            s, c = "good", "A normal doublet rate for a 10x run."
        checks.append(Check(
            "Doublet score (Scrublet)", s,
            f"{pct_dbl:.1f}% of barcodes flagged (median score {med_score:.3f})", c))

    # 6. Before / after filtering -------------------------------------------
    kept = counts["pct_cells_kept"]
    if kept < 50:
        s, c = "concern", ("Filtering removed most of the barcodes. Something "
                           "upstream went wrong, or the thresholds are too strict "
                           "for this tissue.")
    elif kept < 70:
        s, c = "watch", "A sizeable loss. Check which filter is doing the removing."
    else:
        s, c = "good", "Most cells survived filtering."
    checks.append(Check(
        "Cells retained after filtering", s,
        f"{counts['cells_after']:,} of {counts['cells_before']:,} cells kept "
        f"({kept:.1f}%)", c))

    return checks


def overall_verdict(checks: list[Check]) -> tuple[str, str]:
    """Reduce the per-metric checks to one status plus one sentence."""
    worst = max(checks, key=lambda c: _RANK[c.status])
    status = worst.status

    if status == "good":
        return "good", ("All six QC metrics are within normal ranges for a 10x "
                        "run. This sample looks fine to take forward.")

    offenders = [c.metric for c in checks if c.status == status]
    if len(offenders) == 1:
        who = offenders[0].lower()
    else:
        who = ", ".join(o.lower() for o in offenders[:-1]) + f" and {offenders[-1].lower()}"

    panels = "that panel" if len(offenders) == 1 else "those panels"
    if status == "concern":
        return status, (f"The sample has a real problem and the concern is "
                        f"{who}. Read {panels} below before using this data.")
    return status, (f"The sample is broadly usable, but keep an eye on {who}.")


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight",
                facecolor=SURFACE)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _style_axis(ax, value_axis: str = "y") -> None:
    """Recessive grid and axes, so the data is the loudest thing on screen."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(RULE)
        ax.spines[side].set_linewidth(0.8)
    ax.grid(True, axis=value_axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=8, length=3, width=0.8)


def _decade_ticks(ax, lo: float, hi: float, axis: str = "y") -> None:
    """
    Label a log10-transformed axis with the real numbers (10, 100, 1000...).

    We transform the data rather than using set_yscale('log') because a violin's
    density is computed in whatever space it is given; transforming first keeps
    the shape honest.
    """
    lo_e, hi_e = int(np.floor(lo)), int(np.ceil(hi))
    ticks = list(range(lo_e, hi_e + 1))
    labels = [f"{10**t:,.0f}" if t >= 0 else f"{10**t:g}" for t in ticks]
    # Setting ticks outside the data range makes matplotlib widen the view, so
    # remember the real limits and put them back afterwards.
    if axis == "y":
        keep = ax.get_ylim()
        ax.set_yticks(ticks); ax.set_yticklabels(labels)
        ax.set_ylim(keep)
    else:
        keep = ax.get_xlim()
        ax.set_xticks(ticks); ax.set_xticklabels(labels)
        ax.set_xlim(keep)


def metric_figure(values: np.ndarray, title: str, unit: str,
                  log: bool = False, cutoff: float | None = None,
                  cutoff_label: str = "") -> str:
    """
    One metric, drawn twice: a violin (shape of the distribution) beside a
    histogram (where the cells actually pile up). Same data, two readings.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        v = np.zeros(1)

    plot_v = np.log10(np.clip(v, 1e-9, None)) if log else v
    plot_cut = (np.log10(cutoff) if (log and cutoff) else cutoff)

    fig, (ax_v, ax_h) = plt.subplots(
        1, 2, figsize=(7.4, 2.7), gridspec_kw={"width_ratios": [1, 2.2]})
    fig.patch.set_facecolor(SURFACE)

    # --- violin -----------------------------------------------------------
    parts = ax_v.violinplot([plot_v], showextrema=False, showmedians=False,
                            widths=0.75)
    for body in parts["bodies"]:
        body.set_facecolor(SERIES)
        body.set_alpha(0.32)
        body.set_edgecolor(SERIES)
        body.set_linewidth(1.2)

    q1, med, q3 = np.percentile(plot_v, [25, 50, 75])
    ax_v.vlines(1, q1, q3, color=SERIES, linewidth=3, zorder=3)   # IQR
    ax_v.plot(1, med, "o", markersize=7, color=SERIES,            # median,
              markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)  # ringed
    _style_axis(ax_v, "y")
    ax_v.set_xticks([])
    ax_v.spines["bottom"].set_visible(False)
    ax_v.set_ylabel(unit, color=INK_SOFT, fontsize=8)
    if log:
        _decade_ticks(ax_v, plot_v.min(), plot_v.max(), "y")
    if plot_cut is not None:
        ax_v.axhline(plot_cut, color=STATUS["concern"]["color"],
                     linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)

    # --- histogram --------------------------------------------------------
    ax_h.hist(plot_v, bins=60, color=SERIES, alpha=0.85, edgecolor="none")
    _style_axis(ax_h, "y")
    ax_h.set_xlabel(unit, color=INK_SOFT, fontsize=8)
    ax_h.set_ylabel("cells", color=INK_SOFT, fontsize=8)
    if log:
        _decade_ticks(ax_h, plot_v.min(), plot_v.max(), "x")
    if plot_cut is not None:
        ax_h.axvline(plot_cut, color=STATUS["concern"]["color"],
                     linewidth=1.2, linestyle=(0, (4, 3)), zorder=3)
        ax_h.annotate(cutoff_label, xy=(plot_cut, 1), xycoords=("data", "axes fraction"),
                      xytext=(4, -10), textcoords="offset points",
                      fontsize=7.5, color=STATUS["concern"]["color"], ha="left")

    fig.suptitle(title, x=0.005, ha="left", fontsize=10.5, color=INK,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _fig_to_base64(fig)


def filtering_figure(counts: dict) -> str:
    """
    Metric 6, as two small panels rather than one chart.

    Cells and genes are different quantities on different scales. Putting them
    on a shared axis - or worse, two axes - would invite a comparison that
    doesn't mean anything, so they get a panel each.
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.7))
    fig.patch.set_facecolor(SURFACE)

    panels = [
        ("Cells", counts["cells_before"], counts["cells_after"]),
        ("Genes", counts["genes_before"], counts["genes_after"]),
    ]
    for ax, (name, before, after) in zip(axes, panels):
        bars = ax.bar(["Before\nfiltering", "After\nfiltering"], [before, after],
                      width=0.55, color=SERIES, edgecolor=SURFACE, linewidth=2)
        _style_axis(ax, "y")
        ax.set_title(name, fontsize=9.5, color=INK, loc="left", pad=8)
        ax.set_ylim(0, max(before, 1) * 1.18)
        for rect, value in zip(bars, [before, after]):
            ax.annotate(f"{value:,}",
                        xy=(rect.get_x() + rect.get_width() / 2, value),
                        xytext=(0, 4), textcoords="offset points",
                        ha="center", fontsize=9, color=INK, fontweight="bold")

    fig.suptitle("Cell and gene counts, before and after filtering",
                 x=0.005, ha="left", fontsize=10.5, color=INK, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return _fig_to_base64(fig)


# --------------------------------------------------------------------------
# Summary table
# --------------------------------------------------------------------------

def summary_rows(res: QCResult) -> list[dict]:
    """
    The same distributions as numbers.

    A picture of a distribution is easy to read and impossible to quote; this
    table is what goes into a methods section, and it is also the fallback for
    anyone reading the report with a screen reader.
    """
    spec = [
        ("UMI counts per cell", "umi_counts", "{:,.0f}"),
        ("Genes detected per cell", "genes_detected", "{:,.0f}"),
        ("Percent mitochondrial reads", "pct_mito", "{:.1f}"),
        ("Percent ribosomal reads", "pct_ribo", "{:.1f}"),
        ("Doublet score (Scrublet)", "doublet_score", "{:.3f}"),
    ]
    rows = []
    for label, col, fmt in spec:
        v = res.per_cell[col].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            rows.append({"metric": label, "median": "n/a", "iqr": "n/a",
                         "range": "n/a"})
            continue
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        rows.append({
            "metric": label,
            "median": fmt.format(med),
            "iqr": f"{fmt.format(q1)} – {fmt.format(q3)}",
            "range": f"{fmt.format(v.min())} – {fmt.format(v.max())}",
        })
    return rows


# --------------------------------------------------------------------------
# HTML assembly
# --------------------------------------------------------------------------

_CSS = f"""
:root {{ color-scheme: light; }}
* {{ box-sizing: border-box; }}
body {{
  margin: 0; background: #f2f1ee; color: {INK};
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
        Helvetica, Arial, sans-serif;
}}
.wrap {{ max-width: 880px; margin: 0 auto; padding: 40px 20px 72px; }}
h1 {{ font-size: 22px; margin: 0 0 4px; letter-spacing: -0.01em; }}
.sub {{ color: {INK_SOFT}; font-size: 13px; margin: 0 0 28px; }}
.sub code {{ background: #e6e5e1; padding: 1px 5px; border-radius: 4px;
            font-size: 12px; }}
.card {{ background: {SURFACE}; border: 1px solid {RULE}; border-radius: 10px;
        padding: 22px 24px; margin-bottom: 18px; }}
.verdict {{ border-left: 5px solid var(--accent); }}
.verdict .flag {{ display: inline-flex; align-items: center; gap: 8px;
                 font-weight: 700; color: var(--accent); font-size: 13px;
                 text-transform: uppercase; letter-spacing: .06em; }}
.verdict .flag .icon {{
  width: 20px; height: 20px; border-radius: 50%; background: var(--accent);
  color: #fff; display: inline-flex; align-items: center;
  justify-content: center; font-size: 12px; line-height: 1; }}
.verdict p {{ margin: 10px 0 0; font-size: 17px; line-height: 1.5; }}
h2 {{ font-size: 13px; text-transform: uppercase; letter-spacing: .07em;
     color: {INK_FAINT}; margin: 34px 0 12px; font-weight: 700; }}
ul.checks {{ list-style: none; margin: 0; padding: 0; }}
ul.checks li {{ display: flex; gap: 12px; padding: 11px 0;
               border-top: 1px solid {GRID}; }}
ul.checks li:first-child {{ border-top: none; }}
.dot {{ flex: 0 0 auto; width: 18px; height: 18px; border-radius: 50%;
       margin-top: 3px; color: #fff; font-size: 11px; line-height: 18px;
       text-align: center; }}
.checktext strong {{ display: block; font-size: 14px; }}
.checktext .num {{ color: {INK_SOFT}; font-size: 13.5px; }}
.checktext .say {{ color: {INK_SOFT}; font-size: 13.5px; }}
.plot {{ width: 100%; height: auto; display: block; }}
.figcap {{ color: {INK_FAINT}; font-size: 12.5px; margin: 6px 2px 0; }}
.tablewrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13.5px;
        min-width: 460px; }}
th, td {{ text-align: right; padding: 8px 10px;
         border-bottom: 1px solid {GRID}; }}
th:first-child, td:first-child {{ text-align: left; }}
th {{ color: {INK_FAINT}; font-weight: 600; font-size: 12px;
     text-transform: uppercase; letter-spacing: .05em; }}
tbody tr:last-child td {{ border-bottom: none; }}
.foot {{ color: {INK_FAINT}; font-size: 12px; margin-top: 26px;
        line-height: 1.7; }}
"""


def _esc(s) -> str:
    return html.escape(str(s))


def _check_list_html(checks: list[Check]) -> str:
    items = []
    for c in checks:
        st = STATUS[c.status]
        items.append(
            f'<li><span class="dot" style="background:{st["color"]}" '
            f'title="{_esc(st["label"])}">{st["icon"]}</span>'
            f'<span class="checktext"><strong>{_esc(c.metric)} '
            f'&mdash; {_esc(st["label"])}</strong>'
            f'<span class="num">{_esc(c.headline)}</span><br>'
            f'<span class="say">{_esc(c.comment)}</span></span></li>'
        )
    return '<ul class="checks">' + "".join(items) + "</ul>"


def _figure_html(b64: str, alt: str, caption: str) -> str:
    return (f'<div class="card"><img class="plot" alt="{_esc(alt)}" '
            f'src="data:image/png;base64,{b64}">'
            f'<p class="figcap">{caption}</p></div>')


def build_html(res: QCResult) -> str:
    """Assemble the whole report as one string of HTML."""
    checks = build_checks(res)
    status, sentence = overall_verdict(checks)
    st = STATUS[status]
    thr, counts = res.thresholds, res.counts
    pc = res.per_cell

    figs = [
        _figure_html(
            metric_figure(pc["umi_counts"].to_numpy(),
                          "1. UMI counts per cell", "UMIs per cell (log scale)",
                          log=True),
            "Violin and histogram of UMI counts per cell",
            "Total transcripts captured per barcode. A single tight peak is what "
            "you want; a second peak down at the low end is usually ambient RNA."),
        _figure_html(
            metric_figure(pc["genes_detected"].to_numpy(),
                          "2. Genes detected per cell",
                          "genes per cell (log scale)", log=True,
                          cutoff=float(thr.min_genes_per_cell),
                          cutoff_label=f"cutoff {thr.min_genes_per_cell}"),
            "Violin and histogram of genes detected per cell",
            "How many distinct genes each barcode saw. The dashed line is the "
            "minimum-genes filter; cells to its left were removed."),
        _figure_html(
            metric_figure(pc["pct_mito"].to_numpy(),
                          "3. Percent mitochondrial reads", "% of UMIs from MT- genes",
                          cutoff=thr.max_pct_mito,
                          cutoff_label=f"cutoff {thr.max_pct_mito:g}%"),
            "Violin and histogram of percent mitochondrial reads",
            "High mitochondrial fraction means the cell was stressed or its "
            "membrane had broken before capture. The dashed line is your cutoff."),
        _figure_html(
            metric_figure(pc["pct_ribo"].to_numpy(),
                          "4. Percent ribosomal reads", "% of UMIs from RPS/RPL genes"),
            "Violin and histogram of percent ribosomal reads",
            "Ribosomal fraction is strongly tissue-dependent, so there is no "
            "cutoff here - it is shown for context, next to the panels above."),
    ]

    if res.scrublet_ran:
        figs.append(_figure_html(
            metric_figure(pc["doublet_score"].to_numpy(),
                          "5. Doublet score (Scrublet)", "Scrublet doublet score"),
            "Violin and histogram of Scrublet doublet scores",
            "Scrublet compares each barcode against simulated doublets. A long "
            "right-hand tail separated from the main body is the signature of "
            "real doublets."))
    else:
        figs.append(
            f'<div class="card"><strong>5. Doublet score (Scrublet)</strong>'
            f'<p class="figcap">{_esc(res.scrublet_note)}</p></div>')

    figs.append(_figure_html(
        filtering_figure(counts),
        "Bar charts of cell and gene counts before and after filtering",
        f"Filters applied: at least {thr.min_genes_per_cell} genes per cell, "
        f"at most {thr.max_pct_mito:g}% mitochondrial reads, genes kept if seen "
        f"in at least {thr.min_cells_per_gene} cells"
        + (", predicted doublets removed." if thr.drop_predicted_doublets else ".")))

    rows = "".join(
        f"<tr><td>{_esc(r['metric'])}</td><td>{_esc(r['median'])}</td>"
        f"<td>{_esc(r['iqr'])}</td><td>{_esc(r['range'])}</td></tr>"
        for r in summary_rows(res)
    )

    drop_rows = "".join(
        f"<tr><td>{_esc(label)}</td><td>{value:,}</td></tr>"
        for label, value in [
            ("Barcodes in the input matrix", counts["cells_before"]),
            (f"Removed: fewer than {thr.min_genes_per_cell} genes",
             counts["dropped_low_genes"]),
            (f"Removed: above {thr.max_pct_mito:g}% mitochondrial",
             counts["dropped_high_mito"]),
            ("Removed: predicted doublets", counts["dropped_doublets"]),
            ("Cells remaining", counts["cells_after"]),
        ]
    )

    stamp = _dt.datetime.now().strftime("%d %b %Y, %H:%M")

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>scRNA-seq QC &mdash; {_esc(res.sample_name)}</title>
<style>{_CSS}</style>
</head><body><div class="wrap">

<h1>Single-cell QC report &mdash; {_esc(res.sample_name)}</h1>
<p class="sub">Matrix: <code>{_esc(res.matrix_dir)}</code> &middot; generated {stamp}</p>

<div class="card verdict" style="--accent:{st['color']}">
  <span class="flag"><span class="icon">{st['icon']}</span>{_esc(st['label'])}</span>
  <p>{_esc(sentence)}</p>
</div>

<h2>Metric by metric</h2>
<div class="card">{_check_list_html(checks)}</div>

<h2>Distributions</h2>
{''.join(figs)}

<h2>The same numbers, as a table</h2>
<div class="card"><div class="tablewrap"><table>
<thead><tr><th>Metric</th><th>Median</th><th>Middle 50%</th><th>Full range</th></tr></thead>
<tbody>{rows}</tbody></table></div></div>

<h2>What filtering removed</h2>
<div class="card"><div class="tablewrap"><table>
<thead><tr><th>Step</th><th>Cells</th></tr></thead>
<tbody>{drop_rows}</tbody></table></div>
<p class="figcap">Each cell is counted once, against the first filter it failed.</p></div>

<p class="foot">
Six metrics, no more: UMI counts, genes detected, percent mitochondrial,
percent ribosomal, Scrublet doublet score, and cell/gene counts around
filtering. The verdict is a triage aid based on ordinary 10x rules of thumb -
unusual tissues legitimately break those rules, so read the plots before
acting on it.<br>
Built with scanpy. Doublet scores from Scrublet (Wolock, Lopez &amp; Klein 2019)
as implemented in <code>scanpy.pp.scrublet</code>.
</p>

</div></body></html>
"""


def write_report(res: QCResult, out_path: str) -> str:
    """Render the report and write it to disk. Returns the path written."""
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(build_html(res))
    return out_path
