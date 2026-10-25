#!/usr/bin/env python3
"""Numeric localisation of a byte mismatch between a produced file and its golden.

Used by scripts/check_goldens.sh when `cmp` reports a difference. It never writes
anything; it prints one line:

    first-diff line L, element <index>, new=<a> golden=<b>; max |abs| diff D at <index>,
    max rel diff R at <index> (<class>)

File layouts understood (all legacy LSQUANT text formats):

  .chebmom1D   line 1 header "N BandWidth BandCenter", line 2 "M", then one moment
               "re im" per line; element k -> m = k.
  .chebmom2D   line 1 header, line 2 "M0 M1", then the flat moment vector "re im" per
               line; element k -> (m, n) = divmod(k, M1)   (2D row-major index m*M1 + n).
  .chebmomTD   line 1 header "N BW BC Nt dt", line 2 "M Nt", then "re im" per line,
               m-major; element k -> (m, t) = divmod(k, Nt).
  anything else (.dat, .CSR, BOUNDS): whitespace-separated tokens; element = (row, column).

Relative difference of two numbers a, b is |a - b| / max(|a|, |b|) (0 when both are 0).
`nan` equals `nan`; a nan/inf facing a finite number counts as an infinite difference.
Class: "identical values (formatting only)" when every token parses to the same number,
"rounding level (< 1e-12)" when the max relative difference is below 1e-12 (a pure
floating-point reordering, AGENTS.md section 8: needs Fable's explicit acceptance), and
"DEFECT level (>= 1e-12)" otherwise.

Usage: golden_diff.py <new_file> <golden_file>
"""
import math
import os
import sys


def _num(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def _rel(a, b):
    if math.isnan(a) and math.isnan(b):
        return 0.0, 0.0
    if not (math.isfinite(a) and math.isfinite(b)):
        return (0.0, 0.0) if a == b else (math.inf, math.inf)
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d, (d / s if s > 0 else 0.0)


def _layout(path, lines):
    ext = os.path.splitext(path)[1]
    if ext in (".chebmom1D", ".chebmom2D", ".chebmomTD") and len(lines) >= 2:
        second = lines[1].split()
        if ext == ".chebmom1D":
            return 2, lambda k: "m=%d" % k
        if ext == ".chebmom2D":
            m1 = int(second[1])
            return 2, lambda k: "(m,n)=(%d,%d)" % divmod(k, m1)
        nt = int(second[1])
        return 2, lambda k: "(m,t)=(%d,%d)" % divmod(k, nt)
    return 0, None


def _where(header, mapper, line_idx, col):
    """Human-readable element index of token `col` on 0-based line `line_idx`."""
    if mapper is not None and line_idx >= header:
        part = "re" if col == 0 else "im"
        return "line %d %s %s" % (line_idx + 1, mapper(line_idx - header), part)
    return "line %d col %d" % (line_idx + 1, col + 1)


def main():
    if len(sys.argv) != 3:
        print(__doc__.strip().splitlines()[-1])
        return 2
    new_path, gold_path = sys.argv[1], sys.argv[2]
    with open(new_path) as f:
        new = f.read().splitlines()
    with open(gold_path) as f:
        gold = f.read().splitlines()
    header, mapper = _layout(gold_path, gold)

    first = None
    for i in range(max(len(new), len(gold))):
        a = new[i] if i < len(new) else "<EOF>"
        b = gold[i] if i < len(gold) else "<EOF>"
        if a != b:
            ta, tb = a.split(), b.split()
            col = next((j for j in range(max(len(ta), len(tb)))
                        if (ta[j] if j < len(ta) else None) != (tb[j] if j < len(tb) else None)), 0)
            va = ta[col] if col < len(ta) else "<none>"
            vb = tb[col] if col < len(tb) else "<none>"
            first = "first-diff %s, new=%s golden=%s" % (_where(header, mapper, i, col), va, vb)
            break
    if first is None:
        print("no line differs (difference is in line endings or trailing bytes)")
        return 0

    if len(new) != len(gold):
        print("%s; STRUCTURE differs: %d lines new vs %d golden" % (first, len(new), len(gold)))
        return 1

    max_abs = (0.0, None)
    max_rel = (0.0, None)
    for i, (la, lb) in enumerate(zip(new, gold)):
        ta, tb = la.split(), lb.split()
        if len(ta) != len(tb):
            print("%s; STRUCTURE differs at line %d (%d vs %d tokens)" % (first, i + 1, len(ta), len(tb)))
            return 1
        for j, (x, y) in enumerate(zip(ta, tb)):
            if x == y:
                continue
            a, b = _num(x), _num(y)
            if a is None or b is None:
                print("%s; NON-NUMERIC token differs at line %d: %r vs %r" % (first, i + 1, x, y))
                return 1
            d, r = _rel(a, b)
            if d > max_abs[0]:
                max_abs = (d, _where(header, mapper, i, j))
            if r > max_rel[0]:
                max_rel = (r, _where(header, mapper, i, j))
    if max_rel[0] == 0.0:
        cls = "identical values (formatting only)"
    elif max_rel[0] < 1e-12:
        cls = "rounding level (< 1e-12)"
    else:
        cls = "DEFECT level (>= 1e-12)"
    print("%s; max abs diff %.3e at %s; max rel diff %.3e at %s; %s"
          % (first, max_abs[0], max_abs[1], max_rel[0], max_rel[1], cls))
    return 1


if __name__ == "__main__":
    sys.exit(main())
