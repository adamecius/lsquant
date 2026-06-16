#!/usr/bin/env python3
"""lsquant_report - present a unified lsquant results.json nicely.

Reads one (or several) ``<system>.results.json`` produced by ``lsquant run`` and:

  1. prints a summary table (system, fingerprint, observables, timing, peak RSS,
     any errors);
  2. plots the final spectra it carries (DOS, conductivities);
  3. writes a self-contained HTML dashboard ``<system>.report.html`` with the
     plots embedded (base64 PNG -- no server, no JS dependency, opens anywhere).

The style mirrors examples/lsqstyle.py so a report looks like the tutorials.

Usage:
    lsquant_report.py graphene.results.json                # summary + HTML
    lsquant_report.py *.results.json -o dashboard.html     # several systems
    lsquant_report.py run.results.json --no-html           # summary only
"""

import argparse
import base64
import html
import io
import json
import os
import sys

# --- style (mirrors examples/lsqstyle.py) ------------------------------------
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]


def _apply_style():
    import matplotlib as mpl
    from cycler import cycler
    mpl.rcParams.update({
        "savefig.dpi": 130, "savefig.bbox": "tight",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.prop_cycle": cycler(color=PALETTE),
        "axes.grid": True, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": "#DDDDDD", "grid.linewidth": 0.6,
        "lines.linewidth": 1.7, "legend.frameon": False, "legend.fontsize": 9,
    })


def _human_bytes(n):
    try:
        n = float(n)
    except (TypeError, ValueError):
        return "n/a"
    if n < 0:
        return "n/a"
    for u in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024 or u == "TiB":
            return f"{n:.1f} {u}"
        n /= 1024


# --- plotting ----------------------------------------------------------------
def _plot_observable(name, block):
    """Return a PNG (bytes) for a known observable block, or None."""
    import matplotlib.pyplot as plt
    fig = None
    if "energy_eV" in block and "dos" in block:
        fig, ax = plt.subplots(figsize=(7.0, 3.4))
        ax.plot(block["energy_eV"], block["dos"], color=PALETTE[0])
        ax.set_xlabel("energy  (eV)"); ax.set_ylabel("DOS  (1/eV)")
        ax.set_title("Density of states")
    elif "energy_eV" in block and "sigma" in block:
        fig, ax = plt.subplots(figsize=(7.0, 3.4))
        ax.plot(block["energy_eV"], block["sigma"], color=PALETTE[1])
        ax.axhline(0.0, color="#BBBBBB", lw=0.8)
        ax.set_xlabel("energy  (eV)")
        ax.set_ylabel(r"$\sigma$  (%s)" % block.get("units", "e^2/h"))
        ax.set_title("Conductivity %s (%s)" % (name.replace("sigma_", ""),
                                               block.get("formula", "")))
    if fig is None:
        return None
    buf = io.BytesIO()
    fig.tight_layout(); fig.savefig(buf, format="png", dpi=130)
    plt.close(fig)
    return buf.getvalue()


# --- text summary ------------------------------------------------------------
def summarize(doc, path):
    sys_ = doc.get("system", {})
    run = doc.get("run", {})
    print(f"\n=== {path} ===")
    print(f"  system      : {sys_.get('label', '?')}")
    fp = sys_.get("fingerprint", "")
    print(f"  fingerprint : {fp[:16]}{'...' if fp else '(none)'}")
    print(f"  size        : {sys_.get('size', '?')}  "
          f"bandwidth {sys_.get('bandwidth', '?')}  bandcenter {sys_.get('bandcenter', '?')}")
    results = doc.get("results", {})
    print(f"  observables : {len(results)}")
    for nm, blk in results.items():
        t = blk.get("timing_s")
        tt = f"{t:.2f} s" if isinstance(t, (int, float)) else "?"
        print(f"      - {nm:<16} M={blk.get('num_moments', '?'):<6} time={tt}")
    print(f"  run         : wall {run.get('wall_s', '?')} s  "
          f"peak RSS {_human_bytes(run.get('peak_rss_bytes', -1))}")
    errs = run.get("errors", [])
    if errs:
        print(f"  errors      : {len(errs)}")
        for e in errs:
            print(f"      ! {e.get('observable', '?')}: {e.get('error', '')}")


# --- HTML dashboard ----------------------------------------------------------
_CSS = """
body{font-family:DejaVu Sans,Arial,sans-serif;margin:2rem;color:#222;background:#fafafa}
h1{font-size:1.4rem}h2{font-size:1.1rem;margin-top:1.6rem}
.card{background:#fff;border:1px solid #e3e3e3;border-radius:8px;padding:1rem 1.2rem;margin:1rem 0}
table{border-collapse:collapse}td,th{padding:.25rem .8rem;text-align:left}
th{color:#666;font-weight:600}code{background:#f0f0f0;padding:.1rem .3rem;border-radius:4px}
img{max-width:100%;height:auto}.err{color:#b00}
.small{color:#888;font-size:.85rem}
"""


def _img_tag(png):
    b64 = base64.b64encode(png).decode("ascii")
    return f'<img src="data:image/png;base64,{b64}"/>'


def build_html(docs):
    parts = ["<!doctype html><meta charset='utf-8'>",
             "<title>lsquant report</title>",
             f"<style>{_CSS}</style>",
             "<h1>lsquant results</h1>"]
    for path, doc in docs:
        sys_ = doc.get("system", {})
        run = doc.get("run", {})
        results = doc.get("results", {})
        fp = sys_.get("fingerprint", "")
        parts.append("<div class='card'>")
        parts.append(f"<h2>{html.escape(str(sys_.get('label', '?')))}</h2>")
        parts.append(f"<p class='small'>{html.escape(path)}</p>")
        parts.append("<table>")
        parts.append(f"<tr><th>fingerprint</th><td><code>{fp[:24]}</code></td></tr>")
        parts.append(f"<tr><th>size</th><td>{sys_.get('size', '?')}</td></tr>")
        parts.append(f"<tr><th>bandwidth / center</th><td>{sys_.get('bandwidth', '?')} / "
                     f"{sys_.get('bandcenter', '?')} eV</td></tr>")
        parts.append(f"<tr><th>wall / peak RSS</th><td>{run.get('wall_s', '?')} s / "
                     f"{_human_bytes(run.get('peak_rss_bytes', -1))}</td></tr>")
        parts.append("</table>")

        for nm, blk in results.items():
            t = blk.get("timing_s")
            tt = f"{t:.2f} s" if isinstance(t, (int, float)) else "?"
            parts.append(f"<h3>{html.escape(nm)} <span class='small'>"
                         f"M={blk.get('num_moments', '?')}, {tt}</span></h3>")
            try:
                png = _plot_observable(nm, blk)
            except Exception as exc:  # plotting must never crash the report
                png = None
                parts.append(f"<p class='err'>plot failed: {html.escape(str(exc))}</p>")
            if png:
                parts.append(_img_tag(png))
            elif "note" in blk:
                parts.append(f"<p class='small'>{html.escape(str(blk['note']))}</p>")

        errs = run.get("errors", [])
        for e in errs:
            parts.append(f"<p class='err'>! {html.escape(str(e.get('observable', '?')))}: "
                         f"{html.escape(str(e.get('error', '')))}</p>")
        parts.append("</div>")
    return "\n".join(parts)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Present lsquant results.json")
    ap.add_argument("results", nargs="+", help="one or more *.results.json")
    ap.add_argument("-o", "--output", help="HTML output path "
                    "(default: <first-system>.report.html)")
    ap.add_argument("--no-html", action="store_true", help="summary only")
    args = ap.parse_args(argv)

    docs = []
    for p in args.results:
        with open(p) as f:
            docs.append((p, json.load(f)))
        summarize(docs[-1][1], p)

    if args.no_html:
        return 0
    try:
        _apply_style()
        html_text = build_html(docs)
    except ImportError:
        print("\n(matplotlib not available: skipping HTML; install numpy+matplotlib)",
              file=sys.stderr)
        return 0
    out = args.output
    if not out:
        label = docs[0][1].get("system", {}).get("label", "lsquant")
        out = f"{label}.report.html"
    with open(out, "w") as f:
        f.write(html_text)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
