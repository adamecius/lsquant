#!/usr/bin/env python3
"""Thin wrapper: wannier2sparse tools/hr_exactdiag.py on the Gamma-centred k mesh k = n/N.

hr_exactdiag.py evaluates `dos`, `sigma`, `ahc` and `shc` on the shifted mesh
k = (n + 1/2)/N, written inline as `(np.arange(nk) + .5) / nk` (in `_spectral`) and
`(np.arange(a.nk) + .5) / a.nk` (in the commands). The stage-5 gate (G9, handoff 0005
codex.md section 2.2) needs the mesh of the 48 x 48 supercell, k = n/N. This wrapper reads
the tool's source, replaces exactly those mesh expressions by `np.arange(nk) / nk` and
`np.arange(a.nk) / a.nk` (it aborts if the expected occurrences are not found, so a change
of the tool cannot pass silently), and runs the tool's own `main` with the remaining
arguments. Nothing else of the tool changes: same O(k), eigensolver, Gaussian, prefactors
and output.

usage: scripts/hr_exactdiag_gamma.py <path to hr_exactdiag.py> <hr_exactdiag arguments...>
"""
import sys

src_path = sys.argv[1]
src = open(src_path).read()
subs = [("(np.arange(nk) + .5) / nk", "np.arange(nk) / nk", 1),
        ("(np.arange(a.nk) + .5) / a.nk", "np.arange(a.nk) / a.nk", 3)]
for old, new, count in subs:
    found = src.count(old)
    if found != count:
        sys.exit("hr_exactdiag_gamma: expected %d occurrence(s) of %r in %s, found %d" % (count, old, src_path, found))
    src = src.replace(old, new)
if "+ .5" in src:
    sys.exit("hr_exactdiag_gamma: a shifted-mesh expression is left in %s" % src_path)
ns = {"__name__": "hr_exactdiag_gamma_module", "__file__": src_path}
exec(compile(src, src_path, "exec"), ns)
ns["main"](sys.argv[2:])
