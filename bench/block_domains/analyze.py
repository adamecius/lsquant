"""Summarise results.jsonl of run_scaling.sh as Markdown tables."""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

here = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
rows = [json.loads(line) for line in (here / "results.jsonl").read_text().splitlines() if line.strip()]
runs = defaultdict(dict)
for r in rows:
    runs[(r["coupled"], r["block"], r["mode"], r["method"])][r["threads"]] = r

threads = [1, 2, 4, 8, 16, 32, 64]
out = []
for coupled in (False, True):
    out.append(f"\n### {'With' if coupled else 'Without'} neighbour-block coupling\n")
    out.append("| block | observable | method | nnz(H) | " + " | ".join(f"{t} c (s)" for t in threads)
               + " | speedup 64 | efficiency 64 |")
    out.append("|" + "---|" * (6 + len(threads)))
    for w in (5, 10, 30):
        for mode in ("dos", "cond"):
            for method in ("global", "local"):
                d = runs.get((coupled, w, mode, method), {})
                if not d:
                    continue
                cells = [f"{d[t]['total_s']:.2f}" if t in d else "-" for t in threads]
                sp = d[1]["total_s"] / d[64]["total_s"] if 1 in d and 64 in d else float("nan")
                nnz = next(iter(d.values()))["h_nonzeros"]
                out.append(f"| {w} | {mode} | {method} | {nnz:,} | " + " | ".join(cells)
                           + f" | {sp:.1f} | {sp / 64:.0%} |")

out.append("\n### Local / global time ratio without coupling (below 1: local faster)\n")
out.append("| block | observable | " + " | ".join(f"{t} c" for t in threads) + " |")
out.append("|" + "---|" * (2 + len(threads)))
for w in (5, 10, 30):
    for mode in ("dos", "cond"):
        g = runs.get((False, w, mode, "global"), {})
        loc = runs.get((False, w, mode, "local"), {})
        cells = [f"{loc[t]['total_s'] / g[t]['total_s']:.2f}" if t in g and t in loc else "-" for t in threads]
        out.append(f"| {w} | {mode} | " + " | ".join(cells) + " |")

out.append("\n### Cost of the coupling for the global method (coupled / uncoupled time)\n")
out.append("| block | observable | " + " | ".join(f"{t} c" for t in threads) + " |")
out.append("|" + "---|" * (2 + len(threads)))
for w in (5, 10, 30):
    for mode in ("dos", "cond"):
        a = runs.get((False, w, mode, "global"), {})
        b = runs.get((True, w, mode, "global"), {})
        cells = [f"{b[t]['total_s'] / a[t]['total_s']:.2f}" if t in a and t in b else "-" for t in threads]
        out.append(f"| {w} | {mode} | " + " | ".join(cells) + " |")


def diff(p, q):
    a = np.fromfile(p, dtype=np.complex128)
    b = np.fromfile(q, dtype=np.complex128)
    return np.max(np.abs(a - b)) / np.max(np.abs(a))


out.append("\n### Moment agreement (max |difference| / max |moment|)\n")
out.append("| block | observable | local vs global, no coupling, 64 c | global 1 c vs 64 c | local ignoring coupling vs coupled global |")
out.append("|---|---|---|---|---|")
for w in (5, 10, 30):
    for mode in ("dos", "cond"):
        o = here / "out"
        f = lambda m, c, t: o / f"{mode}-{m}-w{w}-c{c}-t{t}.bin"  # noqa: E731
        vals = []
        for p, q in ((f("global", 0, 64), f("local", 0, 64)), (f("global", 0, 1), f("global", 0, 64)),
                     (f("global", 1, 64), f("local", 1, 64))):
            vals.append(f"{diff(p, q):.1e}" if p.exists() and q.exists() else "-")
        out.append(f"| {w} | {mode} | " + " | ".join(vals) + " |")

print("\n".join(out))
