#!/usr/bin/env bash
# Golden-master gate (AGENTS.md section 8). Runs the stage's programs on the committed
# golden inputs in a temporary directory and compares every produced file byte for byte
# with goldens/stage0/. It never writes into goldens/ and never modifies a golden.
#
# usage: scripts/check_goldens.sh <stage> [build_dir]
#
#   stage 0  legacy build (build_dir holds the inline_* drivers): reruns the G1, G2
#            (chain1d, exact trace, M = 64) and G5 (Haldane, exact trace, M = 128)
#            pipelines at OMP_NUM_THREADS=16, then G1 and G2 again at 1 and 64 threads.
#   stage 1  new core (build_dir/lsquant_core/stage1_moments, handoff 0001): G1 and G2
#            at OMP_NUM_THREADS = 1, 8, 64; G2 again with --batch 16 at 8 threads (the
#            batched result must equal the unbatched G2); random-phase determinism
#            (two runs, cmp);
#            grep for "basis|std::cout|getenv" under lsquant/core (must be empty).
#   stage 2  Python binding (build_dir/lsquant_core/_core*.so, handoff 0002): stage 1 in
#            full, then python3 -m lsquant.apps.stage2_dos in its two forms at
#            OMP_NUM_THREADS = 1 and 16 (PYTHONPATH=build_dir/lsquant_core:python):
#            (a) --moments G1 -> JACKSON.dat vs G3; (b) --operators -> .chebmom1D vs G1
#            and JACKSON.dat vs G3; grep for "delta|DOS|Kubo|basis|std::cout|getenv"
#            under lsquant/core and lsquant/io on code lines only (lines whose first
#            non-blank characters are //, * or /* are comments, excluded and listed;
#            the code-line result must be empty).
#   stage 3  rank-2 reconstruction (handoff 0003): stages 1 and 2 in full, then
#            python3 -m lsquant.apps.stage3_kubo at OMP_NUM_THREADS = 16 and 1:
#            G5 moments (Haldane VX-VY and VX-VX, window [-4.5, 4.5], exact trace,
#            512 states, M = 128) vs goldens/stage0/haldane; G5 curves bastin, sea,
#            surf (VX-VY) and greenwood (VX-VX) at 10 meV from the golden moment
#            files; G4 moments (graphene VX-VX, --window gershgorin, random phase
#            seed 12345, 1 state, M = 64); G4 curves bastin at 64 meV, sea and surf
#            at 10 meV (sea_surf/), greenwood_fft and bastin_fft at 10 meV
#            (fft_grid/); G4c (graphene VX-VY, window from goldens/stage0/graphene/
#            offcentre_enclosing/BOUNDS = -8.5 7.5, exact trace, 256 states, M = 128)
#            moments (legacy), and its Bastin curve at 10 meV from the golden and from
#            the new moments vs the curve blessed from the new code (legacy defect D2) (a DIFFERS curve also gets an informational value/energy column
#            note); grep -rn "compat/" lsquant/core lsquant/io must be empty.
#   stage 4  time evolution (handoff 0004): stages 1 to 3 in full, then
#            python3 -m lsquant.apps.stage4_time at OMP_NUM_THREADS = 16 and 1:
#            G6 moments (chain1d_ballistic, --word msd --velocity VX, window [-2, 2],
#            exact trace, 256 states, M = 256, 64 steps, tmax 20) vs the golden
#            .chebmomTD; the four G6 curves (10 meV, E_F = 0, 0.5, 1.0, 1.5) from the
#            golden moment file and again from the produced one; G7 moments
#            (chain1d_mag, --word projected, SX, SY, SZ with PX, --window gershgorin,
#            exact trace, 128 states, M = 32, 96 steps, tmax 95); grep -rnwE
#            "time|velocity|spin|projector|bessel|gsl" (whole lowercase words) on code
#            lines of lsquant/core and lsquant/io must be empty.
#   stage 5  Bloch and eigen bases (handoff 0005, codex.md revised after run 1): stages 1 to
#            4 in full, then python3 -m lsquant.apps.stage5_bases on graphene and Haldane
#            (goldens/stage0/bases/<model>/model, unchanged wannier2sparse gen_models.py files),
#            supercell 48 48 1, M = 256, --window eigen, 40 meV, real basis = the wannier2sparse
#            CSR (--csr-dir), at OMP_NUM_THREADS = 16 and 1: W = 0 real and Bloch (random phase,
#            32 states, seed 12345), Bloch and eigen (exact trace), DOS and kuboxx; W = 1 eV
#            (--disorder 1 --disorder-seed 7) real and Bloch; FFTW convention and the
#            wannier2sparse CSR rows (counted); real-path files of T = 1 and 16 byte-identical
#            (Bloch and eigen reported); figures/stage_5_table.py: PASS/FAIL rows counted
#            (real vs Bloch moments 1e-9, Bloch vs eigen curves 1e-8, EigenBasis exact vs G9
#            1e-8), INFO rows reported (KPM vs G9, finite size); grep -rni fftw
#            lsquant/core/chebyshev lsquant/io must be empty. Needs W2S_BIN and W2S_TOOLS
#            (default ~/wannier2sparse). S5_KEEP=<dir> keeps the outputs for figures/stage_5.py.
#            Stage "5only" runs the stage-5 part alone.
#
#   stage 6  non-orthogonal metric (handoff 0006): stages 1 to 5 in full, then
#            build_dir/lsquant_core/stage6_nonorth on the rotated fixtures
#            goldens/stage0/nonorth/{chain1d,haldane}/{real,complex}: routes lowdin and
#            generalised (naive as evidence only), appliers exact, ischo, chebyshev,
#            moments vs G1, G2 (chain1d dos, kubo) and G5 (Haldane kuboxy = VX-VY,
#            kubo = VX-VX) at the brief's tolerances (not bytes; see the stage-6 block),
#            the two routes against each other, and the Haldane Hall plateau C from the
#            stage-3 Bastin reconstruction of every VX-VY file (|C - 1| <= 0.05).
#            The chain1d runs are repeated at OMP_NUM_THREADS=1 and cmp'd with the main
#            runs (S6_T1_CHAIN1D=0 skips this). S6_KEEP=<dir> keeps the runs;
#            S6_THREADS (default 16) sets OMP_NUM_THREADS of the main runs.
#            Stage "6only" runs the stage-6 part alone.
#
#   stage 7  words, run files, CLI and MCP (handoff 0007): for every run file under
#            runs/goldens/, `python3 -m lsquant run <file> --out <tmp>` at OMP_NUM_THREADS =
#            16 and 1 (S7_THREADS), then cmp of every produced file against its golden: G1, G2,
#            G3, G4, G4c, G5, G6, G7 (goldens/stage0), tutorials 01 to 06 (every non-PNG file
#            of goldens/stage0/tutorials/0N, same relative path; captured/index.tsv by sha1);
#            G9 run files against stage5_bases run with the run file's arguments, including its
#            csr_dir (byte), the csr_dir CSR files against a fresh wannier2sparse run, plus
#            figures/stage_5_table.py on the run-file outputs; G10 run files against
#            stage6_nonorth with the run file's resolved metric (byte; computed once, at the first
#            thread count, into S7_REF6) plus max relative difference vs G1/G2/G5 at the stage-6
#            brief's tolerance (exact 1e-10, ischo 1e-9, chebyshev 1e-8). Run files listed in
#            S7_T1_SKIP (default: the two dense exact Lowdin chain1d files) run at the first thread
#            count only. Then the MCP server (g1.json, g5.json through the `run` and `status`
#            tools, outputs cmp'd; a model dir outside --root must be refused), the CLI refusals
#            (non-zero exit, stderr, no output dir), the greps of stages 1 to 5 and the stage-7
#            greps, and codex.md section 2.5 (stage7_deleted: no src/ include/ test/, the
#            legacy-final tag, a worktree of it in S7_LEGACY_WT builds the legacy drivers and
#            reproduces the stage-0 G1, G2, G5 pipelines). Totals per golden and thread count.
#            S7_KEEP=<dir> keeps the run-file outputs (figures/stage_7.py). Stages "7aux" (MCP and
#            refusals), "7deleted" (section 2.5) and "7tail" (MCP, refusals, greps, section 2.5)
#            run one part. S7_SET=list runs the run files named in S7_LIST. S7_SET=quick runs the
#            closure subset of handoff 0007 (g1 g2 g3 g4 g4c g5 g7_s{x,y,z}
#            g9_graphene_{real,eigen} g10_haldane_real_generalised tutorial01 tutorial03
#            tutorial04) at 16 threads; S7_SET=all (default) runs every run file.
#
#            Stage "7runs" runs the run-file part alone (setup, the run files of S7_SET at
#            S7_THREADS, the per-golden comparison, totals; no MCP, refusals, greps, section 2.5):
#            independent "7runs" invocations with disjoint S7_LIST plus one "7tail" make up "7".
#
# Tolerance mode (GOLDEN_TOL=1; final validation on a machine that is not byte-identical to the
# fat node where the goldens were blessed, e.g. tcn-aronnode: gcc 9.4, conda Eigen 3.4, glibc 2.31).
# Every compare() that is not byte-identical is classified by a helper (tolcmp.py, written into
# $TMP) and reported as TOLPASS (counted as passing) or DIFFERS (counted as failing), with the
# scripts/golden_diff.py numbers either way:
#   moment files (*.chebmom*, 15 significant digits): the two header lines identical as text,
#     same layout (lines, tokens per line), and max_k |new_k - gold_k| <= 1e-12 * max_k |gold_k|
#     over the complex moments (the stage-6 measure);
#   text curves (*.dat, 6 significant digits): same layout, the energy column (first column)
#     identical as text, every other value within one unit of the last printed significant digit
#     (unit = 10^(floor(log10|x|) - P + 1), P = the largest number of significant digits printed
#     in the two files, the larger unit of the two values; absolute, never relative to the value);
#   every other file (JSON, TSV, CSR, index files): non-numeric text identical, every number within
#     one unit of the last printed significant digit as above.
#   captured/index.tsv (sha1 column): an entry whose sha1 does not match passes if the produced
#     capture passes the rule above against the golden capture of the same name.
# The stage-6 thread rows (OMP_NUM_THREADS=1 vs S6_THREADS on the same machine) stay byte
# comparisons. Totals report byte-identical and tolerance-PASS counts separately.
#
# Output: one line per compared file, IDENTICAL or DIFFERS; for DIFFERS the first
# differing bytes (cmp -l | head -5) and the numeric localisation of
# scripts/golden_diff.py (first differing line, element index, both values, max
# relative difference over the file). Exit code 0 only if every file is IDENTICAL
# (and, for stage >= 1, the grep is empty).
#
# KPM_SEED=12345 always. Seeds, tolerances and goldens are never changed here.
set -uo pipefail

usage() { echo "usage: scripts/check_goldens.sh <stage> [build_dir]"; exit 2; }
[ $# -ge 1 ] || usage
STAGE="$1"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${2:-$REPO/build}" 2>/dev/null && pwd)" || { echo "ERROR: build dir '${2:-build}' not found"; exit 2; }
GOLD="$REPO/goldens/stage0"
export KPM_SEED=12345
[ -d "$GOLD" ] || { echo "ERROR: $GOLD missing"; exit 2; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/check_goldens.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
NFAIL=0
NCMP=0
GOLDEN_TOL="${GOLDEN_TOL:-0}"
NBYTE=0; NTOLP=0; CMP_LAST=""; FORCE_BYTE=0
TOLCMP="$TMP/tolcmp.py"
cat >"$TOLCMP" <<'TOLEOF'
# tolerance classifier of check_goldens.sh (GOLDEN_TOL=1); prints "TOLPASS ..." or "DIFFERS ...", exit 0/1
import math, os, re, sys
new_p, gold_p = sys.argv[1], sys.argv[2]
NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|[-+]?(?:nan|inf)', re.I)
def toks(line):
    out, pos = [], 0
    for m in NUM.finditer(line):
        if m.start() > pos: out.append(("t", line[pos:m.start()]))
        out.append(("n", m.group(0))); pos = m.end()
    if pos < len(line): out.append(("t", line[pos:]))
    return out
def sig(tok):
    m = re.match(r'[-+]?(\d*)\.?(\d*)', tok)
    d = (m.group(1) + m.group(2)).lstrip("0") if m else ""
    return len(d)
def done(ok, msg):
    print(("TOLPASS " if ok else "DIFFERS ") + msg); sys.exit(0 if ok else 1)
new = open(new_p, errors="replace").read().split("\n")
gold = open(gold_p, errors="replace").read().split("\n")
if len(new) != len(gold):
    done(False, "layout: %d lines new vs %d golden" % (len(new), len(gold)))
ext = os.path.splitext(gold_p)[1]
if ext.startswith(".chebmom"):
    if new[:2] != gold[:2]:
        done(False, "moment header differs: new %r golden %r" % (new[:2], gold[:2]))
    a, b = [], []
    for i, (x, y) in enumerate(zip(new[2:], gold[2:])):
        tx, ty = x.split(), y.split()
        if len(tx) != len(ty):
            done(False, "layout: line %d has %d vs %d tokens" % (i + 3, len(tx), len(ty)))
        if not tx: continue
        a.append(complex(float(tx[0]), float(tx[1]) if len(tx) > 1 else 0.0))
        b.append(complex(float(ty[0]), float(ty[1]) if len(ty) > 1 else 0.0))
    g = max((abs(v) for v in b), default=0.0)
    d = [abs(x - y) for x, y in zip(a, b)]
    k = max(range(len(d)), key=lambda i: d[i]) if d else 0
    r = (d[k] / g) if g > 0 else (0.0 if not d or d[k] == 0 else math.inf)
    ok = all(math.isfinite(abs(v)) for v in a) and r <= 1e-12
    done(ok, "moments: max|d|/max|gold| = %.3e (<= 1e-12 required), max|d| = %.3e at moment %d, max|gold| = %.3e, header identical"
         % (r, d[k] if d else 0.0, k, g))
P = 1
for y in gold + new:
    for kind, t in toks(y):
        if kind == "n": P = max(P, sig(t))
def unit(x):
    if x == 0 or not math.isfinite(x): return 0.0
    return 10.0 ** (math.floor(math.log10(abs(x))) - P + 1)
curve = ext == ".dat"
worst, nd, nbad, where, maxabs = 0.0, 0, 0, "-", 0.0
for i, (x, y) in enumerate(zip(new, gold)):
    if x == y: continue
    tx, ty = toks(x), toks(y)
    if len(tx) != len(ty) or any(p[0] != q[0] for p, q in zip(tx, ty)):
        done(False, "layout: line %d tokenises differently: %r vs %r" % (i + 1, x[:120], y[:120]))
    col = -1
    for (kx, sx), (ky, sy) in zip(tx, ty):
        if kx == "n": col += 1
        if sx == sy: continue
        if kx == "t":
            done(False, "text differs at line %d: %r vs %r" % (i + 1, sx[:80], sy[:80]))
        if curve and col == 0:
            done(False, "energy column differs at line %d: %s vs %s" % (i + 1, sx, sy))
        a, b = float(sx), float(sy)
        if not (math.isfinite(a) and math.isfinite(b)):
            done(False, "non-finite value differs at line %d col %d: %s vs %s" % (i + 1, col + 1, sx, sy))
        u = max(unit(a), unit(b))
        dd = abs(a - b); nd += 1; maxabs = max(maxabs, dd)
        q = dd / u if u > 0 else math.inf
        nbad += q > 1.0 + 1e-9
        if q > worst: worst, where = q, "line %d col %d new=%s golden=%s" % (i + 1, col + 1, sx, sy)
ok = worst <= 1.0 + 1e-9
gmax = 0.0
for y in gold:
    col = -1
    for kind, t in toks(y):
        if kind != "n": continue
        col += 1
        if curve and col == 0: continue
        try: v = abs(float(t))
        except ValueError: continue
        if math.isfinite(v): gmax = max(gmax, v)
done(ok, "%s: %d values differ, %d beyond one unit, max %.3g units of the last printed digit (P = %d significant digits; <= 1 required) at %s, max|d| = %.3e, max|d|/max|gold| = %.3e (informative)"
     % ("curve" if curve else "text", nd, nbad, worst, P, where, maxabs, maxabs / gmax if gmax > 0 else 0.0))
TOLEOF

# compare <produced file> <golden file> <tag>
compare() {
  local new="$1" gold="$2" tag="$3"
  NCMP=$((NCMP + 1))
  if [ ! -f "$new" ]; then
    echo "DIFFERS   $tag  (not produced: $(basename "$new"))"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS; return
  fi
  if [ ! -f "$gold" ]; then
    echo "DIFFERS   $tag  (no golden: ${gold#$REPO/})"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS; return
  fi
  if cmp -s "$new" "$gold"; then
    echo "IDENTICAL $tag"; NBYTE=$((NBYTE + 1)); CMP_LAST=IDENTICAL
  elif [ "$GOLDEN_TOL" = 1 ] && [ "$FORCE_BYTE" = 0 ]; then
    local t; t="$(python3 "$TOLCMP" "$new" "$gold" 2>&1)"
    if [ "${t%% *}" = TOLPASS ]; then
      echo "TOLPASS   $tag"; NTOLP=$((NTOLP + 1)); CMP_LAST=TOLPASS
    else
      echo "DIFFERS   $tag"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS
    fi
    echo "          tolerance: ${t#* }"
    echo "          golden_diff.py: $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  else
    NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS
    echo "DIFFERS   $tag"
    echo "          cmp: $(cmp "$new" "$gold" 2>&1 | head -1)"
    cmp -l "$new" "$gold" 2>/dev/null | head -5 | sed 's/^/          cmp -l: /'
    echo "          $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  fi
}

# exact-trace state file over the full basis of a CSR (legacy "local" state set)
write_exact() { local n; n=$(head -1 "$1" | awk '{print $1}'); { echo local; echo "$n"; seq 0 $((n - 1)); } > "$2"; }

# ---------------------------------------------------------------- stage 0 (legacy)
# G1 + G2: the chain1d part of scripts/regenerate_goldens.sh, on the golden operators.
legacy_chain1d() {   # <threads> <workdir>
  local T="$1" W="$2"
  mkdir -p "$W/operators"
  cp "$GOLD/chain1d/operators/"*.CSR "$W/operators/"
  cp "$GOLD/chain1d/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/chain1d.HAM.CSR" "$W/exact"
  ( cd "$W"
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-nonEqOp"   chain1d VX VX 64 exact >nonEq.log 2>&1
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-spectralOp" chain1d 1 64 exact   >spectral.log 2>&1 )
  compare "$W/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "[T=$T] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)"
  compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "[T=$T] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)"
}

# G5: Haldane, exact trace, M = 128, bounds [-4.5, 4.5] (test/hall_haldane.sh), 10 meV.
legacy_haldane() {   # <threads> <workdir>
  local T="$1" W="$2" f
  mkdir -p "$W/operators"
  cp "$GOLD/haldane/operators/"haldane.{HAM,VX,VY}.CSR "$W/operators/"
  cp "$GOLD/haldane/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/haldane.HAM.CSR" "$W/exact"
  ( cd "$W"
    export OMP_NUM_THREADS=$T
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VY 128 exact >xy.log 2>&1
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VX 128 exact >xx.log 2>&1
    XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
    XX=NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
    "$BUILD/inline_kuboBastinFromChebmom"     "$XY" 10 >bastin.log 2>&1
    "$BUILD/inline_kuboBastinSeaFromChebmom"  "$XY" 10 >sea.log 2>&1
    "$BUILD/inline_kuboBastinSurfFromChebmom" "$XY" 10 >surf.log 2>&1
    "$BUILD/inline_kuboGreenwoodFromChebmom"  "$XX" 10 >greenwood.log 2>&1 )
  for f in NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
           NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D; do
    compare "$W/$f" "$GOLD/haldane/$f" "[T=$T] haldane/$f (G5)"
  done
  for f in "$GOLD"/haldane/*.dat; do
    compare "$W/$(basename "$f")" "$f" "[T=$T] haldane/$(basename "$f") (G5)"
  done
}

stage0() {
  for b in inline_compute-kpm-nonEqOp inline_compute-kpm-spectralOp inline_kuboBastinFromChebmom \
           inline_kuboBastinSeaFromChebmom inline_kuboBastinSurfFromChebmom inline_kuboGreenwoodFromChebmom; do
    [ -x "$BUILD/$b" ] || { echo "ERROR: legacy driver $BUILD/$b missing (build with the legacy tree)"; exit 2; }
  done
  echo "== stage 0: legacy pipelines vs goldens/stage0 (KPM_SEED=$KPM_SEED, build $BUILD)"
  legacy_chain1d 16 "$TMP/g12_T16"
  legacy_haldane 16 "$TMP/g5_T16"
  legacy_chain1d 1  "$TMP/g12_T1"
  legacy_chain1d 64 "$TMP/g12_T64"
}

# ---------------------------------------------------------------- stage 1 (new core)
stage1() {
  local EXE="$BUILD/lsquant_core/stage1_moments" T W w R
  [ -x "$EXE" ] || EXE="$BUILD/lsquant/stage1_moments"
  [ -x "$EXE" ] || { echo "ERROR: stage1_moments not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 1: $EXE vs goldens/stage0/chain1d (KPM_SEED=$KPM_SEED)"
  for T in 1 8 64; do
    W="$TMP/s1_T$T"; mkdir -p "$W"
    for w in dos kubo; do
      OMP_NUM_THREADS=$T "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
        --trace exact --states 2048 --moments 64 --word $w --out "$W" >"$W/$w.log" 2>&1 \
        || echo "          [T=$T] stage1_moments --word $w exited non-zero: $(tail -1 "$W/$w.log")"
    done
    compare "$W/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
            "$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
            "[T=$T] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)"
    compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
            "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
            "[T=$T] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)"
  done
  # G2 with batching (batchSize 16 < M = 64), T = 8: must equal the unbatched G2
  W="$TMP/s1_batch16"; mkdir -p "$W"
  OMP_NUM_THREADS=8 "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
    --trace exact --seed 12345 --states 2048 --moments 64 --word kubo --batch 16 --out "$W" >"$W/kubo.log" 2>&1 \
    || echo "          [T=8 batch 16] stage1_moments --word kubo exited non-zero: $(tail -1 "$W/kubo.log")"
  compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "[T=8, batch 16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2, batched)"
  # random-phase determinism: no golden; two identical runs are the criterion
  for R in 1 2; do
    W="$TMP/s1_random$R"; mkdir -p "$W"
    for w in dos kubo; do
      OMP_NUM_THREADS=8 "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
        --trace random --seed 12345 --states 4 --moments 64 --word $w --out "$W" >"$W/$w.log" 2>&1 \
        || echo "          [random run $R] stage1_moments --word $w exited non-zero: $(tail -1 "$W/$w.log")"
    done
  done
  local n=0 f
  for f in "$TMP/s1_random1"/*.chebmom*; do
    [ -f "$f" ] || continue; n=$((n + 1))
    compare "$TMP/s1_random2/$(basename "$f")" "$f" "[random, run2 vs run1] $(basename "$f")"
  done
  [ "$n" -gt 0 ] || { echo "DIFFERS   [random] no .chebmom output produced"; NFAIL=$((NFAIL + 1)); }
}

core_grep() {
  local hits
  if [ ! -d "$REPO/lsquant/core" ]; then
    echo "GREP      lsquant/core missing"; NFAIL=$((NFAIL + 1)); return
  fi
  hits="$(cd "$REPO" && grep -rn "basis\|std::cout\|getenv" lsquant/core)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"basis\\|std::cout\\|getenv\" lsquant/core: empty"
  else
    echo "GREP      grep -rn \"basis\\|std::cout\\|getenv\" lsquant/core: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

# ---------------------------------------------------------------- stage 2 (Python binding)
G1F=SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
G3F=meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat

stage2() {
  local PYDIR="$BUILD/lsquant_core" T W rc
  ls "$PYDIR"/_core*.so >/dev/null 2>&1 || PYDIR="$BUILD/lsquant"
  ls "$PYDIR"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 2: python3 -m lsquant.apps.stage2_dos, PYTHONPATH=$PYDIR:python, vs goldens/stage0/chain1d (KPM_SEED=$KPM_SEED)"
  for T in 1 16; do
    # form (a): reconstruct the G3 curve from the committed G1 moment file
    W="$TMP/s2_T${T}_a"; mkdir -p "$W"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR:$REPO/python" python3 -m lsquant.apps.stage2_dos \
        --moments "$GOLD/chain1d/$G1F" --broadening-meV 120 --out "$W" ) >"$W/run.log" 2>&1; rc=$?
    [ $rc -eq 0 ] || echo "          [T=$T] stage2_dos form (a) exited $rc: $(tail -1 "$W/run.log")"
    sed "s|$W/||; s|^|          [T=$T] (a) printed: |" "$W/run.log"
    compare "$W/$G3F" "$GOLD/chain1d/$G3F" "[T=$T] (a) --moments G1 -> chain1d/$G3F (G3)"
    # form (b): moments through MomentEngine from the golden operators, then reconstruct
    W="$TMP/s2_T${T}_b"; mkdir -p "$W"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR:$REPO/python" python3 -m lsquant.apps.stage2_dos \
        --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 --trace exact --states 2048 \
        --moments 64 --broadening-meV 120 --out "$W" ) >"$W/run.log" 2>&1; rc=$?
    [ $rc -eq 0 ] || echo "          [T=$T] stage2_dos form (b) exited $rc: $(tail -1 "$W/run.log")"
    sed "s|$W/||; s|^|          [T=$T] (b) printed: |" "$W/run.log"
    compare "$W/$G1F" "$GOLD/chain1d/$G1F" "[T=$T] (b) --operators -> chain1d/$G1F (G1)"
    compare "$W/$G3F" "$GOLD/chain1d/$G3F" "[T=$T] (b) --operators -> chain1d/$G3F (G3)"
  done
}

# physics names and side channels in the engine, on code lines only: a line whose first
# non-blank characters are //, * or /* is a comment (e.g. the notation delta_{m0}) and is
# excluded from the verdict, but listed below it.
stage2_grep() {
  local all hits coms pat='delta\|DOS\|Kubo\|basis\|std::cout\|getenv'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  if [ ! -d "$REPO/lsquant/core" ] || [ ! -d "$REPO/lsquant/io" ]; then
    echo "GREP      lsquant/core or lsquant/io missing"; NFAIL=$((NFAIL + 1)); return
  fi
  all="$(cd "$REPO" && grep -rn "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else
    echo "GREP      grep -rn \"$pat\" lsquant/core lsquant/io, code lines: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
}

# ---------------------------------------------------------------- stage 3 (rank-2 reconstruction)
# One stage3_kubo invocation in its own directory; prints the app's output lines.
run_s3() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc; shift 3
  mkdir -p "$W"
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR3:$REPO/python" python3 -m lsquant.apps.stage3_kubo \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  [ $rc -eq 0 ] || echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"
  sed "s|$W/||; s|^|          [T=$T] $tag printed: |" "$W/run.log"
}

G5XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact
G5XX=NonEqOpVX-VXhaldaneKPM_M128x128_stateexact
G4XX=NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault
G4CXY=NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact

# Informational, not a verdict: for a two-column curve that DIFFERS, say whether the value
# column is byte-identical as text and what the energy-column offset (golden - new) is.
columns_note() {   # <new> <golden>
  [ -f "$1" ] && [ -f "$2" ] || return 0
  cmp -s "$1" "$2" && return 0
  local v; cmp -s <(awk '{print $2}' "$1") <(awk '{print $2}' "$2") && v="identical" || v="DIFFERS"
  echo "          note: value column $v as text; energy column golden - new: $(paste "$1" "$2" | awk '{d=$3-$1; s+=d; if(NR==1||d<mn)mn=d; if(NR==1||d>mx)mx=d} END{printf "min %.6g max %.6g mean %.6g over %d lines", mn, mx, s/NR, NR}')"
}

stage3() {
  local T W G
  PYDIR3="$BUILD/lsquant_core"
  ls "$PYDIR3"/_core*.so >/dev/null 2>&1 || PYDIR3="$BUILD/lsquant"
  ls "$PYDIR3"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 3: python3 -m lsquant.apps.stage3_kubo, PYTHONPATH=$PYDIR3:python, vs goldens/stage0/{haldane,graphene}"
  for T in 16 1; do
    # 1. G5 moments: Haldane, exact trace over the 512 sites, M = 128, window [-4.5, 4.5]
    W="$TMP/s3_T${T}_g5m"
    run_s3 $T "$W" "G5 VX-VY moments" --operators "$GOLD/haldane/operators" --label haldane --left VX --right VY \
      --window -4.5 4.5 --trace exact --states 512 --moments 128
    run_s3 $T "$W" "G5 VX-VX moments" --operators "$GOLD/haldane/operators" --label haldane --left VX --right VX \
      --window -4.5 4.5 --trace exact --states 512 --moments 128
    compare "$W/$G5XY.chebmom2D" "$GOLD/haldane/$G5XY.chebmom2D" "[T=$T] haldane/$G5XY.chebmom2D (G5 moments)"
    compare "$W/$G5XX.chebmom2D" "$GOLD/haldane/$G5XX.chebmom2D" "[T=$T] haldane/$G5XX.chebmom2D (G5 moments)"
    # 2. G5 curves from the committed golden moment files, 10 meV
    W="$TMP/s3_T${T}_g5c"
    run_s3 $T "$W" "G5 bastin"    --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula bastin    --broadening-meV 10
    run_s3 $T "$W" "G5 sea"       --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula sea       --broadening-meV 10
    run_s3 $T "$W" "G5 surf"      --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula surf      --broadening-meV 10
    run_s3 $T "$W" "G5 greenwood" --moments "$GOLD/haldane/$G5XX.chebmom2D" --formula greenwood --broadening-meV 10
    for G in KuboBastin_$G5XY KuboBastinSea_$G5XY KuboBastinSurf_$G5XY KuboGreenWood_$G5XX; do
      compare "$W/${G}JACKSON.dat" "$GOLD/haldane/${G}JACKSON.dat" "[T=$T] haldane/${G}JACKSON.dat (G5 curve)"
    done
    # 3. G4 moments: graphene, random phase seed 12345, 1 state, M = 64, Gershgorin window (pad 0.10)
    W="$TMP/s3_T${T}_g4m"
    run_s3 $T "$W" "G4 VX-VX moments" --operators "$GOLD/graphene/operators" --label graphene_golden --left VX --right VX \
      --window gershgorin --trace random --seed 12345 --states 1 --moments 64
    compare "$W/$G4XX.chebmom2D" "$GOLD/graphene/$G4XX.chebmom2D" "[T=$T] graphene/$G4XX.chebmom2D (G4 moments)"
    # 4. G4 curves from the committed golden moment file
    W="$TMP/s3_T${T}_g4c"
    run_s3 $T "$W" "G4 bastin 64 meV"     --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula bastin        --broadening-meV 64
    run_s3 $T "$W" "G4 sea 10 meV"        --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula sea           --broadening-meV 10
    run_s3 $T "$W" "G4 surf 10 meV"       --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula surf          --broadening-meV 10
    run_s3 $T "$W" "G4 greenwood_fft 10 meV" --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula greenwood_fft --broadening-meV 10
    run_s3 $T "$W" "G4 bastin_fft 10 meV" --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula bastin_fft    --broadening-meV 10
    compare "$W/KuboBastin_${G4XX}JACKSON.dat" "$GOLD/graphene/KuboBastin_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/KuboBastin_${G4XX}JACKSON.dat (G4 curve, 64 meV)"
    compare "$W/KuboBastinSea_${G4XX}JACKSON.dat" "$GOLD/graphene/sea_surf/KuboBastinSea_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/sea_surf/KuboBastinSea_${G4XX}JACKSON.dat (G4 curve, 10 meV)"
    compare "$W/KuboBastinSurf_${G4XX}JACKSON.dat" "$GOLD/graphene/sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat (G4 curve, 10 meV)"
    compare "$W/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat" "$GOLD/graphene/fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat" \
            "[T=$T] graphene/fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat (G4 curve, 10 meV)"
    compare "$W/KuboBastin_${G4XX}JACKSON_FFTgrid.dat" "$GOLD/graphene/fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat" \
            "[T=$T] graphene/fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat (G4 curve, 10 meV)"
    # 5. G4c: graphene VX-VY, off-centre window from the BOUNDS file (-8.5 7.5), exact trace
    #    256 states, M = 128: moments vs the legacy moments; Bastin at 10 meV, from the golden
    #    and from the produced moments, vs the curve blessed from the new code (legacy defect
    #    D2, goldens/stage0/legacy_defects/offcentre_axis/).
    for G in offcentre_enclosing; do
      W="$TMP/s3_T${T}_g4c_$G"
      run_s3 $T "$W" "G4c $G VX-VY moments" --operators "$GOLD/graphene/operators" --label graphene_golden --left VX --right VY \
        --window $(cat "$GOLD/graphene/$G/BOUNDS") --trace exact --states 256 --moments 128
      compare "$W/$G4CXY.chebmom2D" "$GOLD/graphene/$G/$G4CXY.chebmom2D" "[T=$T] graphene/$G/$G4CXY.chebmom2D (G4c moments)"
      run_s3 $T "$W/cg" "G4c $G bastin (golden moments)" --moments "$GOLD/graphene/$G/$G4CXY.chebmom2D" --formula bastin --broadening-meV 10
      run_s3 $T "$W/cn" "G4c $G bastin (new moments)" --moments "$W/$G4CXY.chebmom2D" --formula bastin --broadening-meV 10
      for c in cg cn; do
        compare "$W/$c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/$G/KuboBastin_${G4CXY}JACKSON.dat" \
                "[T=$T] graphene/$G/KuboBastin_${G4CXY}JACKSON.dat (G4c curve, 10 meV, from $([ $c = cg ] && echo golden || echo new) moments)"
        columns_note "$W/$c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/$G/KuboBastin_${G4CXY}JACKSON.dat"
      done
    done
  done
}

# lsquant/compat holds legacy formulas for the gate only; nothing in core/ or io/ may include it.
stage3_grep() {
  local hits
  hits="$(cd "$REPO" && grep -rn "compat/" lsquant/core lsquant/io)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"compat/\" lsquant/core lsquant/io: empty"
  else
    echo "GREP      grep -rn \"compat/\" lsquant/core lsquant/io: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

# ---------------------------------------------------------------- stage 4 (time evolution)
# One stage4_time invocation in its own directory; prints the app's output lines.
run_s4() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc; shift 3
  mkdir -p "$W"
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR4:$REPO/python" python3 -m lsquant.apps.stage4_time \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  [ $rc -eq 0 ] || echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"
  sed "s|$W/||; s|^|          [T=$T] $tag printed: |" "$W/run.log"
  : >"$W/run.log"
}

G6M=CorrelationVXchain1dKPM_M256_state20exact
G6EF="0.000000 0.500000 1.000000 1.500000"

stage4() {
  local T W ef OP s e
  PYDIR4="$BUILD/lsquant_core"
  ls "$PYDIR4"/_core*.so >/dev/null 2>&1 || PYDIR4="$BUILD/lsquant"
  ls "$PYDIR4"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 4: python3 -m lsquant.apps.stage4_time, PYTHONPATH=$PYDIR4:python, vs goldens/stage0/{chain1d_ballistic,chain1d_mag}"
  for T in 16 1; do
    # 1. G6 moments: chain1d ballistic MSD, VX, window [-2, 2], exact trace 256 states, M = 256, 64 steps, tmax 20
    W="$TMP/s4_T${T}_g6m"; s=$(date +%s)
    run_s4 $T "$W" "G6 msd moments" --operators "$GOLD/chain1d_ballistic/operators" --label chain1d --velocity VX \
      --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd
    e=$(date +%s); echo "          [T=$T] G6 msd moments wall $((e - s)) s"
    compare "$W/$G6M.chebmomTD" "$GOLD/chain1d_ballistic/$G6M.chebmomTD" "[T=$T] chain1d_ballistic/$G6M.chebmomTD (G6 moments)"
    # 1b. G6 curves (10 meV, E_F = 0, 0.5, 1.0, 1.5) from the committed golden moment file
    W="$TMP/s4_T${T}_g6c_gold"
    run_s4 $T "$W" "G6 curves (golden moments)" --moments "$GOLD/chain1d_ballistic/$G6M.chebmomTD" \
      --broadening-meV 10 --fermi-energy 0 0.5 1.0 1.5
    for ef in $G6EF; do
      compare "$W/mean${G6M}EF${ef}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat" \
              "[T=$T] chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat (G6 curve from golden moments)"
    done
    # 1c. the same four curves from the moment file just produced in step 1
    W="$TMP/s4_T${T}_g6c_new"
    run_s4 $T "$W" "G6 curves (new moments)" --moments "$TMP/s4_T${T}_g6m/$G6M.chebmomTD" \
      --broadening-meV 10 --fermi-energy 0 0.5 1.0 1.5
    for ef in $G6EF; do
      compare "$W/mean${G6M}EF${ef}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat" \
              "[T=$T] chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat (G6 curve from new moments)"
    done
    # 2. G7 moments: chain1d_mag, S_a with projector PX, Gershgorin window (legacy SpectralBounds
    #    without BOUNDS file: [-2.31, 2.31]), exact trace 128 states, M = 32, 96 steps, tmax 95
    W="$TMP/s4_T${T}_g7m"
    for OP in SX SY SZ; do
      s=$(date +%s)
      run_s4 $T "$W" "G7 $OP-PX moments" --operators "$GOLD/chain1d_mag/operators" --label chain1d_mag --operator $OP \
        --projector PX --window gershgorin --trace exact --states 128 --moments 32 --steps 96 --tmax 95 --word projected
      e=$(date +%s); echo "          [T=$T] G7 $OP-PX moments wall $((e - s)) s"
    done
    for OP in SX SY SZ; do
      compare "$W/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD" \
              "$GOLD/chain1d_mag/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD" \
              "[T=$T] chain1d_mag/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)"
    done
  done
}

# physics names in the engine (stage 4), on code lines only, whole lowercase words; the
# identifiers TimeAxis, EvolveRegister, EvolveDerived, ApplyDerivedSeries are allowed (they
# are not lowercase whole words, so -w with a case-sensitive pattern never matches them).
stage4_grep() {
  local all hits coms pat='time|velocity|spin|projector|bessel|gsl'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  all="$(cd "$REPO" && grep -rnwE "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rnwE \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else
    echo "GREP      grep -rnwE \"$pat\" lsquant/core lsquant/io, code lines: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
}

# ---------------------------------------------------------------- stage 5 (Bloch and eigen bases)
# G9: graphene and Haldane (wannier2sparse examples/gen_models.py, committed under
# goldens/stage0/bases/<model>/model/), supercell 48 x 48 x 1, M = 256, window
# `--window eigen` (the app's default), broadening 40 meV (= G9's Gaussian eta).
S5N=48
S5_W2S="${W2S_BIN:-$HOME/wannier2sparse/build/wannier2sparse}"
S5_BMEV=40

# Python helper of the stage-5 gate, written into $TMP (it never touches the repository).
s5_helpers() {
  mkdir -p "$TMP/s5py"
  # FFTW convention of BlochTransform against the dense B (6 x 4 x 1, non-square on purpose),
  # and the real basis against the wannier2sparse supercell CSR at 6 x 4 x 1 and 48 x 48 x 1
  # (brief section 3: the real basis IS the wannier2sparse CSR, B VX_w2s B^+ = v_x to 1e-10).
  # Prints PASS/FAIL rows (counted by the gate) and INFO rows. The helper name starts with s5_
  # so that it never shadows a standard-library module.
  cat >"$TMP/s5py/s5_conv.py" <<'PY'
import sys, numpy as np
from lsquant import _core
from lsquant.models.wannier import load_wannier_model
from lsquant.bases import make_basis
mdir, seed, w2s_dir, w2s_small = sys.argv[1:5]
m = load_wannier_model(mdir, seed=seed)
uc, tau, W = m.uc, m.tau, m.num_wann
def gate(ok, text): print(("PASS  " if ok else "FAIL  ") + text)
rng = np.random.default_rng(1)
def rv(n): return rng.standard_normal(n) + 1j * rng.standard_normal(n)
def rel(a, b): return np.max(np.abs(a - b)) / np.max(np.abs(b))
# dense B on 6 x 4 x 1, real index (R, alpha) -> iR * W + alpha with iR = n1 + N1 n2 + N1 N2 n3
N = (6, 4, 1); Nk = 24
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
n = np.array([(n1, n2, n3) for n3 in range(N[2]) for n2 in range(N[1]) for n1 in range(N[0])], float)
Rc = n @ uc
for sgn in (-1, +1):
    Bd = np.zeros((Nk * W, Nk * W), complex)
    for al in range(W):
        Bd[al::W, al::W] = np.exp(sgn * 1j * (k @ (Rc + tau[al]).T)) / np.sqrt(Nk)
    xs = [rv(Nk * W) for _ in range(4)]
    e_f = max(rel(Bf.multiply(x), Bd @ x) for x in xs)
    e_b = max(rel(Bb.multiply(x), Bd.conj().T @ x) for x in xs)
    line = ("CONV %s 6x4x1 dense B = exp(%si k.(R+tau))/sqrt(Nk), iR = n1 + N1 n2: RealToBloch vs B x %.3e, BlochToReal vs B^+ x %.3e"
            % (seed, '+' if sgn > 0 else '-', e_f, e_b))
    if sgn < 0:
        gate(max(e_f, e_b) <= 1e-12, line + " (brief.md section 3, tol 1e-12)")
    else:
        print("INFO  " + line + " (opposite sign, reported)")
def csr_checks(N, d, Bf=None, Bb=None):
    dim = N[0] * N[1] * N[2] * W
    Hr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.HAM.CSR" % seed))
    Vr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.VX.CSR" % seed))
    s = "%dx%dx%d" % N
    xs = [rv(dim) for _ in range(4)]
    # the real basis as stage5_bases builds it: with the wannier2sparse CSRs (--csr-dir) and without
    for lab, cd in (("loaded from --csr-dir", d), ("built by lsquant", None)):
        rb = make_basis("real", m, N, csr_dir=cd) if cd else make_basis("real", m, N)   # model dir holds no CSR
        Hp = rb.hamiltonian(); Vp = rb.operator("VX")
        eh = max(rel(Hp.multiply(x), Hr.multiply(x)) for x in xs); ev = max(rel(Vp.multiply(x), Vr.multiply(x)) for x in xs)
        gate(max(eh, ev) <= 1e-12, "CONV %s %s RealBasis (%s) vs wannier2sparse CSR: H %.3e, VX %.3e (tol 1e-12)" % (seed, s, lab, eh, ev))
    Hp = _core.SparseMatrixOperator(m.hamiltonian_csr(N)); Vp = _core.SparseMatrixOperator(m.velocity_csr(N, "VX"))
    eh = max(rel(Hp.multiply(x), Hr.multiply(x)) for x in xs); ev = max(rel(Vp.multiply(x), Vr.multiply(x)) for x in xs)
    gate(max(eh, ev) <= 1e-12, "CONV %s %s lsquant supercell builder (WannierModel.hamiltonian_csr, velocity_csr) vs wannier2sparse CSR: H %.3e, VX %.3e (tol 1e-12)" % (seed, s, eh, ev))
    if Bf is not None:
        Hk = _core.BlockDiagonalOperator(m.Hk(m.k_mesh(N))); Vk = _core.BlockDiagonalOperator(m.dHk(m.k_mesh(N), "VX"))
        e1 = max(rel(Bf.multiply(Hr.multiply(Bb.multiply(x))), Hk.multiply(x)) for x in xs)
        e2 = max(rel(Bf.multiply(Vr.multiply(Bb.multiply(x))), Vk.multiply(x)) for x in xs)
        gate(e1 <= 1e-10, "CONV %s %s B H_w2s B^+ x vs H_k x: %.3e (tol 1e-10)" % (seed, s, e1))
        gate(e2 <= 1e-10, "CONV %s %s B VX_w2s B^+ x vs v_x(k) x: %.3e (tol 1e-10, brief.md section 3)" % (seed, s, e2))
csr_checks((6, 4, 1), w2s_small)
N = (48, 48, 1); Nk = 48 * 48
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
xs = [rv(Nk * W) for _ in range(4)]
e = max(rel(Bb.multiply(Bf.multiply(x)), x) for x in xs)
gate(e <= 1e-12, "CONV %s 48x48x1 max|B^+ B x - x|/max|x| = %.3e (tol 1e-12)" % (seed, e))
csr_checks(N, w2s_dir, Bf, Bb)
PY
}

# One stage5_bases run in its own directory (wall time printed)
run_s5() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc s e; shift 3
  mkdir -p "$W"; s=$(date +%s.%N)
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR5:$REPO/python" python3 -m lsquant.apps.stage5_bases \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  e=$(date +%s.%N)
  [ $rc -eq 0 ] || { echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"; NFAIL=$((NFAIL + 1)); }
  printf "          [T=%s] %-34s wall %7.2f s  -> %s\n" "$T" "$tag" "$(echo "$e - $s" | bc)" "$(ls "$W" | grep -v '\.log$' | tr '\n' ' ')"
}

stage5() {
  local T m w W C MD CSR rel f g
  PYDIR5="$BUILD/lsquant_core"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || PYDIR5="$BUILD/lsquant"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  [ -x "$S5_W2S" ] || { echo "ERROR: wannier2sparse binary $S5_W2S missing (set W2S_BIN)"; exit 2; }
  echo "== stage 5: python3 -m lsquant.apps.stage5_bases, PYTHONPATH=$PYDIR5:python, vs goldens/stage0/bases (G9)"
  s5_helpers
  for m in graphene haldane; do
    MD="$GOLD/bases/$m/model"; CSR="$TMP/s5w2s/$m"
    mkdir -p "$CSR" "$TMP/s5w2s/${m}_6x4"
    ( cd "$MD" && "$S5_W2S" $m $S5N $S5N 1 VX -o "$CSR" --quiet --no-log-file \
      && "$S5_W2S" $m 6 4 1 VX -o "$TMP/s5w2s/${m}_6x4" --quiet --no-log-file ) >"$TMP/s5w2s/$m.log" 2>&1 \
      || echo "          wannier2sparse $m VX failed: $(tail -1 "$TMP/s5w2s/$m.log")"
    echo "   $m: model $MD (unchanged gen_models.py files), real basis = wannier2sparse CSR ($m $S5N $S5N 1 VX) via --csr-dir;"
    echo "      --window eigen (default), --broadening-meV $S5_BMEV, M=256, 32 random-phase states seed 12345, W=1: --disorder 1 --disorder-seed 7"
    ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR5:$REPO/python" python3 "$TMP/s5py/s5_conv.py" "$MD" $m \
        "$CSR" "$TMP/s5w2s/${m}_6x4" ) >"$TMP/s5_conv_$m.txt" 2>&1
    sed 's/^/          /' "$TMP/s5_conv_$m.txt"
    NCMP=$((NCMP + $(grep -c '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt")))
    NFAIL=$((NFAIL + $(grep -c '^FAIL' "$TMP/s5_conv_$m.txt")))
    grep -q '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt" || { echo "DIFFERS   $m: convention checks not produced"; NFAIL=$((NFAIL + 1)); }
    for T in 16 1; do
      for w in dos kuboxx; do
        C="--model $MD --seed $m --seed 12345 --supercell $S5N $S5N 1 --word $w --moments 256 --broadening-meV $S5_BMEV --csr-dir $CSR"
        W="$TMP/s5/T$T/$m/W0"
        run_s5 $T "$W/real/$w"   "$m W=0 real  $w (random 32)" $C --basis real  --trace random --states 32
        run_s5 $T "$W/bloch/$w"  "$m W=0 bloch $w (random 32)" $C --basis bloch --trace random --states 32
        run_s5 $T "$W/blochx/$w" "$m W=0 bloch $w (exact)"     $C --basis bloch --trace exact
        run_s5 $T "$W/eigen/$w"  "$m W=0 eigen $w (exact)"     $C --basis eigen --trace exact
        W="$TMP/s5/T$T/$m/W1"
        run_s5 $T "$W/real/$w"   "$m W=1 real  $w (random 32)" $C --basis real  --trace random --states 32 --disorder 1 --disorder-seed 7
        run_s5 $T "$W/bloch/$w"  "$m W=1 bloch $w (random 32)" $C --basis bloch --trace random --states 32 --disorder 1 --disorder-seed 7
      done
    done
    # thread counts: the real path must be byte-identical between T = 1 and 16 (counted);
    # the Bloch and eigen paths are reported
    for W in "$TMP/s5/T16/$m"/W*/*/*; do
      rel="${W#$TMP/s5/T16/}"
      for f in "$W"/*.chebmom "$W"/*.chebmom1D "$W"/*.chebmom2D "$W"/*.dat "$W"/*.basis.json; do
        [ -f "$f" ] || continue
        g="$TMP/s5/T1/$rel/$(basename "$f")"
        case "$rel" in
          */real/*) compare "$g" "$f" "[T=1 vs T=16] $rel/$(basename "$f")" ;;
          *) if cmp -s "$g" "$f"; then echo "INFO      [T=1 vs T=16] $rel/$(basename "$f") identical"
             else echo "INFO      [T=1 vs T=16] $rel/$(basename "$f") differs: $(python3 "$REPO/scripts/golden_diff.py" "$g" "$f" 2>&1)"; fi ;;
        esac
      done
    done
  done
  # the tolerance table of codex.md section 2, at T = 16 and T = 1
  for T in 16 1; do
    echo "   tolerance table, T=$T"
    ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR5:$REPO/python" python3 "$REPO/figures/stage_5_table.py" "$TMP/s5/T$T" "$GOLD/bases" ) >"$TMP/s5_table_T$T.txt" 2>&1
    sed 's/^/          /' "$TMP/s5_table_T$T.txt"
    NCMP=$((NCMP + $(grep -c '^\(PASS\|FAIL\)' "$TMP/s5_table_T$T.txt")))
    NFAIL=$((NFAIL + $(grep -c '^FAIL' "$TMP/s5_table_T$T.txt")))
    grep -q '^\(PASS\|FAIL\)' "$TMP/s5_table_T$T.txt" || { echo "DIFFERS   [T=$T] stage-5 table not produced"; NFAIL=$((NFAIL + 1)); }
  done
  [ -z "${S5_KEEP:-}" ] || { rm -rf "$S5_KEEP"; cp -r "$TMP/s5" "$S5_KEEP"; echo "   stage-5 outputs kept in $S5_KEEP"; }
}

# FFTW belongs to the algebra layer (BlochTransform) only; the Chebyshev engine and io never see it.
stage5_grep() {
  local hits
  hits="$(cd "$REPO" && grep -rni "fftw" lsquant/core/chebyshev lsquant/io)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rni \"fftw\" lsquant/core/chebyshev lsquant/io: empty"
  else
    echo "GREP      grep -rni \"fftw\" lsquant/core/chebyshev lsquant/io: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

# ---------------------------------------------------------------- stage 6 (non-orthogonal metric)
# Rotated models (goldens/stage0/nonorth/<model>/<variant>, fixtures of handoff 0006) must
# reproduce the orthogonal goldens G1, G2 (chain1d) and the G5 moments (Haldane VX-VY,
# VX-VX). These runs cannot be byte-identical (the applier is a factorisation, an
# iteration or a truncated series), so the criterion is the brief's tolerance (section 6):
#   exact 1e-10; ischo (--tol 1e-12) 1e-9; chebyshev (ctor truncation coeffTol 1e-12) 1e-8.
# Measure (tester's decision, see tests.md):
#   rel = max_k |mu_new[k] - mu_gold[k]| / max_k |mu_gold[k]|  over all complex moments;
# the per-element relative difference of scripts/golden_diff.py is printed as well, but it
# is meaningless on moments that vanish by symmetry (e.g. the odd DOS moments).
# The two header lines (N BandWidth BandCenter / M) must be numerically identical.
# Rows: PASS / FAIL (counted); route naive = the legacy symmetric estimator, EVIDENCE (not
# counted); lowdin + ischo = N/A (ISCHO provides S^-1 only; stage6_nonorth rejects it).
# Then lowdin vs generalised with the same applier (generalised/ischo vs lowdin/exact) at
# the same levels, and the Hall plateau C = (mean sigma_xy over |E| < 0.5) / A * 2 pi,
# A = (dim/2) 3 sqrt(3)/2 (test/hall_haldane.sh, tutorial03.chern_from_plateau), from the
# stage-3 Bastin reconstruction at 10 meV of every produced VX-VY file: PASS if |C-1| <= 0.05.
# With S6_KEEP=<dir> the run directories, gate_table.tsv and chern.tsv are kept there
# (figures/stage_6.py reads them). S6_THREADS sets OMP_NUM_THREADS (default 16).
S6CMP=""
S6_LAST=""
s6_tools() {
  S6CMP="$TMP/s6cmp.py"
  cat > "$S6CMP" <<'PYEOF'
import sys
def load(p):
    L = open(p).read().split("\n")
    h = [L[0].split(), L[1].split()]
    v = []
    for l in L[2:]:
        t = l.split()
        if t:
            v.append(complex(float(t[0]), float(t[1]) if len(t) > 1 else 0.0))
    return h, v
a, b = load(sys.argv[1]), load(sys.argv[2])
hdr = [[float(x) for x in r] for r in a[0]] == [[float(x) for x in r] for r in b[0]]
if len(a[1]) != len(b[1]):
    print("nan header_ok=%s length %d vs %d" % (hdr, len(a[1]), len(b[1])))
    sys.exit()
g = max(abs(x) for x in b[1])
d = [abs(x - y) for x, y in zip(a[1], b[1])]
k = max(range(len(d)), key=lambda i: d[i])
print("%.3e header_ok=%s max_abs=%.3e at_element=%d max|gold|=%.3e" % (d[k] / g, hdr, d[k], k, g))
PYEOF
  S6CHERN="$TMP/s6chern.py"
  cat > "$S6CHERN" <<'PYEOF'
import sys, math
import numpy as np
d = np.loadtxt(sys.argv[1])
dim = int(sys.argv[2])
if not np.all(np.isfinite(d)):
    print("nan nonfinite")
    sys.exit()
E, s = d[:, 0], d[:, 1]
m = np.abs(E) < 0.5
A = (dim / 2.0) * 3.0 * math.sqrt(3.0) / 2.0
p = float(np.mean(s[m]))
print("%.6f plateau=%.6f points=%d" % (p / A * 2 * math.pi, p, int(m.sum())))
PYEOF
}

# s6_judge <new> <gold> <tol> <tag> <gate|evidence>; sets S6_LAST="<rel> <verdict>"
s6_judge() {
  local new="$1" gold="$2" tol="$3" tag="$4" kind="$5" out rel hdr verdict
  if [ ! -f "$new" ] || [ ! -f "$gold" ]; then
    verdict=FAIL; out="nan (missing: $([ -f "$new" ] || basename "$new") $([ -f "$gold" ] || echo "${gold#$TMP/}"))"
  else
    out="$(python3 "$S6CMP" "$new" "$gold" 2>&1)"
    rel="${out%% *}"; hdr="$(echo "$out" | grep -o 'header_ok=[A-Za-z]*')"
    if python3 -c "import sys; r=float('$rel'); sys.exit(0 if (r == r and r <= $tol and '$hdr' == 'header_ok=True') else 1)" 2>/dev/null
    then verdict=PASS; else verdict=FAIL; fi
  fi
  rel="${out%% *}"
  if [ "$kind" = evidence ]; then
    echo "EVIDENCE  $tag  rel=$out  (would be $verdict at tol $tol; not counted)"
  else
    NCMP=$((NCMP + 1))
    [ "$verdict" = PASS ] || NFAIL=$((NFAIL + 1))
    printf '%-9s %s  rel=%s  tol=%s\n' "$verdict" "$tag" "$out" "$tol"
  fi
  if [ -f "$new" ] && [ -f "$gold" ]; then
    echo "          golden_diff.py: $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  fi
  S6_LAST="$rel $verdict"
}

s6_field() { echo "$1" | grep -o " $2=[^ ]*" | head -1 | cut -d= -f2; }

stage6() {
  local EXE="$BUILD/lsquant_core/stage6_nonorth" PYDIR6 model variant word route applier W tol
  local args lo hi M gold tag line rc res it eo f b P C kind
  local XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact
  [ -x "$EXE" ] || EXE="$BUILD/lsquant/stage6_nonorth"
  [ -x "$EXE" ] || { echo "ERROR: stage6_nonorth not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  PYDIR6="$(dirname "$EXE")"
  local RUNEXE="$EXE"
  s6_tools
  local KEEP="${S6_KEEP:-}" TSV="$TMP/s6_gate_table.tsv" CTSV="$TMP/s6_chern.tsv" T6="${S6_THREADS:-16}"
  printf 'model\tvariant\tword\troute\tapplier\tmax_rel\ttol\tverdict\tresidual\titerations\teffective_order\n' > "$TSV"
  echo "== stage 6: $EXE on goldens/stage0/nonorth vs G1, G2, G5 moments (OMP_NUM_THREADS=$T6)"
  for model in chain1d haldane; do
    if [ $model = chain1d ]; then lo=-2; hi=2; M=64; else lo=-4.5; hi=4.5; M=128; fi
    for variant in real complex; do
      for word in $([ $model = chain1d ] && echo "dos kubo" || echo "kuboxy kubo"); do
        case "$model/$word" in
          chain1d/dos)    gold="$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" ;;
          chain1d/kubo)   gold="$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" ;;
          haldane/kuboxy) gold="$GOLD/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D" ;;
          haldane/kubo)   gold="$GOLD/haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D" ;;
        esac
        b="$(basename "$gold")"; P="$TMP/s6/${model}_${variant}_${word}"
        for route in lowdin generalised naive; do
          for applier in exact ischo chebyshev; do
            [ $route = naive ] && [ $applier != exact ] && continue
            tag="$model/$variant $word $route/$applier"
            if [ $route = lowdin ] && [ $applier = ischo ]; then
              echo "N/A       $tag  (ISCHO provides S^-1 only; the Lowdin route needs S^-1/2)"
              printf '%s\t%s\t%s\t%s\t%s\tN/A\t-\tN/A\t-\t-\t-\n' $model $variant $word $route $applier >> "$TSV"
              continue
            fi
            case $applier in
              exact) tol=1e-10; args="--max-dim 4096" ;;
              ischo) tol=1e-9; args="--tol 1e-12 --max-iter 200" ;;
              chebyshev) tol=1e-8; args="--order 60 --coeff-tol 1e-12" ;;
            esac
            W="${P}_${route}_${applier}"; mkdir -p "$W"
            echo "$args" > "$W/args"
            ( cd "$REPO" && OMP_NUM_THREADS=$T6 "$RUNEXE" --fixture "$GOLD/nonorth/$model/$variant" \
                --label $model --window $lo $hi --route $route --applier $applier $args --trace exact \
                --moments $M --word $word --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
            [ $rc -eq 0 ] || echo "          [$tag] stage6_nonorth exited $rc: $(tail -1 "$W/err.log")"
            line="$(grep '^stage6_nonorth:' "$W/run.log" | tail -1)"
            echo "          stdout: ${line//$W\//}"
            kind=gate; [ $route = naive ] && kind=evidence
            s6_judge "$W/$b" "$gold" $tol "$tag" $kind
            res="$(s6_field "$line" residual)"; it="$(s6_field "$line" iterations)"; eo="$(s6_field "$line" effective_order)"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $model $variant $word $route $applier \
              "${S6_LAST% *}" $tol "$([ $kind = evidence ] && echo "EVIDENCE(${S6_LAST#* })" || echo "${S6_LAST#* }")" \
              "${res:--}" "${it:--}" "${eo:--}" >> "$TSV"
          done
        done
        # the two routes against each other
        s6_judge "${P}_generalised_exact/$b" "${P}_lowdin_exact/$b" 1e-10 \
          "$model/$variant $word generalised/exact vs lowdin/exact" gate
        printf '%s\t%s\t%s\tgeneralised_vs_lowdin\texact\t%s\t1e-10\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
        s6_judge "${P}_generalised_chebyshev/$b" "${P}_lowdin_chebyshev/$b" 1e-8 \
          "$model/$variant $word generalised/chebyshev vs lowdin/chebyshev" gate
        printf '%s\t%s\t%s\tgeneralised_vs_lowdin\tchebyshev\t%s\t1e-8\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
        s6_judge "${P}_generalised_ischo/$b" "${P}_lowdin_exact/$b" 1e-9 \
          "$model/$variant $word generalised/ischo vs lowdin/exact" gate
        printf '%s\t%s\t%s\tgeneralised_ischo_vs_lowdin_exact\tischo\t%s\t1e-9\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
      done
    done
  done
  # Thread-count repeat, chain1d rows only (Haldane skipped to save time): the same runs at
  # OMP_NUM_THREADS=1 must be byte-identical to the runs above (cmp).
  if [ "${S6_T1_CHAIN1D:-1}" = 1 ]; then
    echo "== stage 6: chain1d rows repeated at OMP_NUM_THREADS=1, cmp vs OMP_NUM_THREADS=$T6"
    for W in "$TMP"/s6/chain1d_*; do
      [ -f "$W/args" ] || continue
      tag="$(basename "$W")"; f="$(ls "$W"/*.chebmom1D "$W"/*.chebmom2D 2>/dev/null | head -1)"
      set -- $(echo "$tag" | tr '_' ' ')   # chain1d variant word route applier
      mkdir -p "$W/T1"
      ( cd "$REPO" && OMP_NUM_THREADS=1 "$RUNEXE" --fixture "$GOLD/nonorth/chain1d/$2" \
          --label chain1d --window -2 2 --route $4 --applier $5 $(cat "$W/args") --trace exact \
          --moments 64 --word $3 --out "$W/T1" ) >"$W/T1/run.log" 2>"$W/T1/err.log" \
        || echo "          [T=1 $tag] stage6_nonorth exited non-zero: $(tail -1 "$W/T1/err.log")"
      if [ -z "$f" ]; then
        NCMP=$((NCMP + 1)); NFAIL=$((NFAIL + 1)); echo "DIFFERS   [T=1 vs T=$T6] $tag (no T=$T6 output)"
      else
        FORCE_BYTE=1; compare "$W/T1/$(basename "$f")" "$f" "[T=1 vs T=$T6] $tag/$(basename "$f")"; FORCE_BYTE=0
      fi
    done
  fi
  # Hall plateau from every produced Haldane VX-VY file (stage-3 Bastin, 10 meV)
  echo "== stage 6: Hall plateau, python3 -m lsquant.apps.stage3_kubo --formula bastin --broadening-meV 10"
  printf 'source\tC\tverdict\n' > "$CTSV"
  C="$(python3 "$S6CHERN" "$GOLD/haldane/KuboBastin_${XY}JACKSON.dat" 512 2>&1)"
  echo "REFERENCE Hall golden G5 KuboBastin (legacy)  C=$C"
  printf 'golden_G5\t%s\tREF\n' "${C%% *}" >> "$CTSV"
  for W in "$TMP"/s6/haldane_*_kuboxy_*; do
    tag="$(basename "$W")"
    if [ -f "$W/$XY.chebmom2D" ]; then
      ( cd "$REPO" && OMP_NUM_THREADS=$T6 PYTHONPATH="$PYDIR6:$REPO/python" python3 -m lsquant.apps.stage3_kubo \
          --moments "$W/$XY.chebmom2D" --formula bastin --broadening-meV 10 --out "$W" ) >"$W/bastin.log" 2>&1 \
        || echo "          [$tag] stage3_kubo exited non-zero: $(tail -1 "$W/bastin.log")"
      C="$(python3 "$S6CHERN" "$W/KuboBastin_${XY}JACKSON.dat" 512 2>&1 | tail -1)"
    else
      C="nan (no VX-VY moments produced)"
    fi
    if python3 -c "import sys; c=float('${C%% *}'); sys.exit(0 if abs(c - 1) <= 0.05 else 1)" 2>/dev/null
    then rc=PASS; else rc=FAIL; fi
    case "$tag" in
      *_naive_*) echo "EVIDENCE  Hall $tag  C=$C  (would be $rc; not counted)"; rc="EVIDENCE($rc)" ;;
      *) NCMP=$((NCMP + 1)); [ $rc = PASS ] || NFAIL=$((NFAIL + 1))
         printf '%-9s Hall %s  C=%s  (|C - 1| <= 0.05)\n' $rc "$tag" "$C" ;;
    esac
    printf '%s\t%s\t%s\n' "$tag" "${C%% *}" "$rc" >> "$CTSV"
  done
  if [ -n "$KEEP" ]; then
    mkdir -p "$KEEP"; rm -rf "$KEEP/runs"; cp -r "$TMP/s6" "$KEEP/runs"
    cp "$TSV" "$KEEP/gate_table.tsv"; cp "$CTSV" "$KEEP/chern.tsv"
    echo "          kept: $KEEP/{runs,gate_table.tsv,chern.tsv}"
  fi
}


# ---------------------------------------------------------------- stage 7 (words, run files, CLI, MCP)
# Every golden from words: `python -m lsquant run runs/goldens/<name>.json --out <tmp>` and a
# cmp of each produced file with its golden. The old per-stage sections (stage1..stage5) are
# not run at stage 7; the mapping from their compared files to the run files is:
#   stage 1 G1, G2 (stage1_moments)          -> g1.json, g2.json
#   stage 2 G3 (+ G1 via --operators)        -> g3.json (G3 curve, DOS_reconstructed.dat, G1)
#   stage 3 G5 moments + 4 curves            -> g5.json;  G4 moments + 5 curves -> g4.json;
#           G4c moments + Bastin curve       -> g4c.json
#   stage 4 G6 moments + 4 curves            -> g6.json;  G7 3 moments -> g7.json (and g7_s{x,y,z}.json)
#   stage 5 G9 (stage5_bases, no byte golden)-> g9_<model>_<cfg>.json, cmp'd against stage5_bases run
#           with the run file's own arguments, and figures/stage_5_table.py on the run-file outputs
#   stage 6 G10 (stage6_nonorth, exact)      -> g10_<model>_<variant>_<route>.json, cmp'd against
#           stage6_nonorth with the run file's arguments, and vs G1/G2/G5 at 1e-10 (stage-6 brief s. 6)
#   tutorials 01..05 (stage-0 captures)      -> tutorial0N.json, every non-PNG file under
#           goldens/stage0/tutorials/0N at the same relative path
# Not compared from the old sections: stage-1 batch-16 and random-phase determinism, stage-2
# form (a), the stage-3/4 "curves from the golden moment file" rows (no run file reconstructs
# from a golden file; `lsquant reconstruct` is not part of the gate), the stage-5 convention
# rows (s5_conv.py). Totals are printed per golden and per thread count.
S7_THREADS="${S7_THREADS:-16 1}"
# G10 run files whose dense exact Lowdin route takes ~45 min per word at 16 threads are run at the
# first thread count only (their 1-thread repeat would take many hours on a shared node).
S7_T1_SKIP="${S7_T1_SKIP-g10_chain1d_real_lowdin g10_chain1d_complex_lowdin}"
S7_REF6="${S7_REF6:-}"
s7skip() { [ "$1" = "${S7_THREADS%% *}" ] && return 1; case " $S7_T1_SKIP " in *" $2 "*) return 0 ;; esac; return 1; }
S7_RUNS="$REPO/runs/goldens"
G1F2=NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
declare -A S7I S7D S7ORDER S7B S7T
# S7_SET=all (default): every run file. S7_SET=quick: the closure subset of handoff 0007 (16 threads
# by default; g6, g7, the other g9 and g10 files, tutorial02 and tutorial05 deferred to the final run).
# S7_SET=list: the run files named in S7_LIST (space separated), at S7_THREADS.
S7_SET="${S7_SET:-all}"
S7_QUICK="g1 g2 g3 g4 g4c g5 g7_sx g7_sy g7_sz g9_graphene_real g9_graphene_eigen g10_haldane_real_generalised tutorial01 tutorial03 tutorial04"
[ "$S7_SET" = quick ] && S7_THREADS="${S7_THREADS_QUICK:-16}"
s7want() {
  [ "$S7_SET" = all ] && return 0
  case " $([ "$S7_SET" = list ] && echo "$S7_LIST" || echo "$S7_QUICK") " in *" $1 "*) return 0 ;; esac; return 1
}

# captured/index.tsv: every line "<NNNN_name> <sha1> <source>" must name a produced capture with that sha1
s7idx() {   # <group> <produced captured dir> <golden index.tsv> <tag>
  local g="$1" d="$2" idx="$3" tag="$4" bad=0 n=0 tp=0 name sha t
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  NCMP=$((NCMP + 1))
  while IFS=$'\t' read -r name sha _; do
    n=$((n + 1))
    [ -f "$d/$name" ] && [ "$(sha1sum "$d/$name" | cut -d' ' -f1)" = "$sha" ] && continue
    if [ "$GOLDEN_TOL" = 1 ] && [ -f "$d/$name" ] && [ -f "$(dirname "$idx")/$name" ] \
       && t="$(python3 "$TOLCMP" "$d/$name" "$(dirname "$idx")/$name" 2>&1)"; then
      tp=$((tp + 1)); echo "          index.tsv line $n: $name sha1 differs, capture within tolerance: ${t#* }"; continue
    fi
    bad=$((bad + 1)); echo "          index.tsv line $n: $name sha1 $sha not matched by $d/$name"
  done <"$idx"
  if [ $bad -eq 0 ] && [ $tp -eq 0 ] && [ $n -gt 0 ]; then echo "IDENTICAL $tag ($n entries, every sha1 matches)"; S7I[$g]=$(( ${S7I[$g]:-0} + 1 )); S7B[$g]=$(( ${S7B[$g]:-0} + 1 )); NBYTE=$((NBYTE + 1))
  elif [ $bad -eq 0 ] && [ $n -gt 0 ]; then echo "TOLPASS   $tag ($n entries: $((n - tp)) sha1 match, $tp captures within tolerance)"; S7I[$g]=$(( ${S7I[$g]:-0} + 1 )); S7T[$g]=$(( ${S7T[$g]:-0} + 1 )); NTOLP=$((NTOLP + 1))
  else echo "DIFFERS   $tag ($bad of $n entries do not match)"; NFAIL=$((NFAIL + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); fi
}
S7GROUPS=()

s7cmp() {   # <group> <new> <golden> <tag>   (compare(), plus a per-group tally)
  local g="$1" before=$NFAIL
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  compare "$2" "$3" "$4"
  if [ "$NFAIL" -eq "$before" ]; then S7I[$g]=$(( ${S7I[$g]:-0} + 1 )); else S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); fi
  case "$CMP_LAST" in IDENTICAL) S7B[$g]=$(( ${S7B[$g]:-0} + 1 )) ;; TOLPASS) S7T[$g]=$(( ${S7T[$g]:-0} + 1 )) ;; esac
}

# max relative difference as the stage-6 gate defines it (handoff 0006 tests.md 2.1):
# rel = max_k |new_k - gold_k| / max_k |gold_k| over the data (moment files: lines after the
# two header lines; curves: the value column), 0 when cmp-identical. The pointwise relative
# difference of scripts/golden_diff.py is ~1 on moments that vanish by symmetry, so it is not used.
s7relval() {
  cmp -s "$1" "$2" && { echo 0; return; }
  python3 - "$1" "$2" <<'PY'
import sys
def data(p):
    L = open(p).read().split("\n")
    if p.endswith((".chebmom1D", ".chebmom2D", ".chebmomTD")):
        return [float(t) for l in L[2:] for t in l.split()]
    return [float(l.split()[1]) for l in L if len(l.split()) >= 2]
a, b = data(sys.argv[1]), data(sys.argv[2])
if len(a) != len(b): print("inf"); sys.exit()
m = max(abs(x) for x in b) or 1.0
print("%.3e" % (max(abs(x - y) for x, y in zip(a, b)) / m))
PY
}

s7rel() {   # <group> <new> <golden> <tol|INFO> <tag>
  local g="$1" new="$2" gold="$3" tol="$4" tag="$5" r
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  if [ ! -f "$new" ]; then
    [ "$tol" = INFO ] && { echo "INFO      $tag  (not produced)"; return; }
    echo "FAIL      $tag  (not produced)"; NFAIL=$((NFAIL + 1)); NCMP=$((NCMP + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); return
  fi
  r="$(s7relval "$new" "$gold")"; [ -n "$r" ] || r=inf
  if [ "$tol" = INFO ]; then echo "INFO      $tag  max rel diff $r"; return; fi
  NCMP=$((NCMP + 1))
  if awk -v r="$r" -v t="$tol" 'BEGIN{exit !(r <= t)}'; then
    echo "PASS      $tag  max rel diff $r <= $tol"; S7I[$g]=$(( ${S7I[$g]:-0} + 1 ))
  else
    echo "FAIL      $tag  max rel diff $r > $tol"; NFAIL=$((NFAIL + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 ))
  fi
}

run_s7() {   # <threads> <run file name> <outdir>
  local T="$1" n="$2" W="$3" rc s e
  rm -rf "$W"; mkdir -p "$(dirname "$W")"; s=$(date +%s.%N)
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant run "runs/goldens/$n.json" --out "$W" ) \
    >"$W.out" 2>"$W.err"; rc=$?
  e=$(date +%s.%N)
  [ $rc -eq 0 ] || { echo "          [T=$T] lsquant run runs/goldens/$n.json exited $rc: $(tail -1 "$W.err")"; NFAIL=$((NFAIL + 1)); }
  printf "          [T=%s] python -m lsquant run runs/goldens/%-36s wall %8.2f s, %4d files written\n" "$T" "$n.json" "$(echo "$e - $s" | bc)" "$(find "$W" -type f 2>/dev/null | wc -l)"
}

# The G9 real-basis run files name a csr_dir (wannier2sparse output) that lsquant refuses to
# generate itself: the gate generates it, explicitly, for every selected g9 run file whose csr_dir
# lacks <seed>.HAM.CSR, with the command the refusal names (wannier2sparse <seed> N1 N2 N3 VX
# -p <model dir> -o <csr_dir>), and prints what it did.
s7_prep_csr() {
  local f n C m sc
  for f in "$S7_RUNS"/g9_*.json; do
    n="$(basename "$f" .json)"; s7want $n || continue
    read -r C m sc <<<"$(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["model"].get("csr_dir") or "-", d["model"]["seed"], " ".join(map(str, d["basis"].get("supercell", [1,1,1]))))' "$f")"
    [ "$C" = - ] && continue
    case "$C" in /*) ;; *) C="$REPO/$C" ;; esac
    if [ -f "$C/$m.HAM.CSR" ]; then echo "          csr_dir ${C#$REPO/} for $n: present"; continue; fi
    mkdir -p "$C"
    "$S5_W2S" $m $sc VX -p "$GOLD/bases/$m/model" -o "$C" --quiet --no-log-file >"$TMP/s7prep_$m.log" 2>&1 \
      && echo "          csr_dir ${C#$REPO/} for $n: generated by the gate: wannier2sparse $m $sc VX -p goldens/stage0/bases/$m/model -o ${C#$REPO/}" \
      || { echo "DIFFERS   csr_dir ${C#$REPO/}: wannier2sparse failed: $(tail -1 "$TMP/s7prep_$m.log")"; NFAIL=$((NFAIL + 1)); }
  done
}

s7_ref5() {   # <threads> <model> <cfg> <outdir> [csr dir]: stage5_bases with the run file's arguments
  local T="$1" m="$2" c="$3" O="$4" C="${5:-}" w a
  case "$c" in
    real) a="--basis real --trace random --states 32" ;;
    bloch) a="--basis bloch --trace random --states 32" ;;
    bloch_exact) a="--basis bloch --trace exact" ;;
    eigen) a="--basis eigen --trace exact" ;;
    real_w1) a="--basis real --trace random --states 32 --disorder 1 --disorder-seed 7" ;;
    bloch_w1) a="--basis bloch --trace random --states 32 --disorder 1 --disorder-seed 7" ;;
  esac
  [ -z "$C" ] || a="$a --csr-dir $C"
  for w in dos kuboxx; do
    mkdir -p "$O/$w"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant.apps.stage5_bases \
        --model "$GOLD/bases/$m/model" --seed $m --seed 12345 --supercell 48 48 1 --word $w --moments 256 \
        --broadening-meV 40 $a --out "$O/$w" ) >"$O/$w.log" 2>&1 || echo "          [T=$T] stage5_bases $m $c $w failed: $(tail -1 "$O/$w.log")"
  done
}

# stage6_nonorth with the run file's resolved metric (run.resolved.json of the run): route, applier,
# order, coeff_tol, tol, max_iter, max_dim; words dos+kubo (chain1d) or kuboxy+kubo (haldane)
s7_ref6() {   # <threads> <model> <variant> <resolved run file> <outdir>
  local T="$1" m="$2" v="$3" J="$4" O="$5" w win st M ma
  if [ "$m" = chain1d ]; then win="-2 2"; st=2048; M=64; set -- dos kubo; else win="-4.5 4.5"; st=512; M=128; set -- kuboxy kubo; fi
  ma="$(python3 -c 'import json,sys; x=json.load(open(sys.argv[1]))["metric"]; print("--route %s --applier %s --order %d --coeff-tol %r --tol %r --max-iter %d --max-dim %d" % (x["route"], x["applier"], x["order"], x["coeff_tol"], x["tol"], x["max_iter"], x["max_dim"]))' "$J")"
  mkdir -p "$O"; echo "stage6_nonorth $ma" >"$O.args"
  for w in "$@"; do
    OMP_NUM_THREADS=$T "$EXE6" --fixture "$GOLD/nonorth/$m/$v" --label $m --window $win $ma \
      --trace exact --states $st --moments $M --word $w --out "$O" >>"$O.log" 2>&1 \
      || echo "          [T=$T] stage6_nonorth $m $v $w failed: $(tail -1 "$O.log")"
  done
}

s7_goldens() {   # <threads>: every run file, every golden
  local T="$1" B="$TMP/s7/T$1" n f g m c v r W R C tol a
  [ "$T" = "${S7_THREADS%% *}" ] && s7_prep_csr
  for f in "$S7_RUNS"/*.json; do
    n="$(basename "$f" .json)"; s7want $n || continue
    if s7skip $T $n; then echo "          [T=$T] $n.json not run at this thread count (S7_T1_SKIP)"; continue; fi
    run_s7 $T $n "$B/$n"
  done
  echo "   [T=$T] per-golden comparison"
  s7want g1 && s7cmp "T$T G1" "$B/g1/$G1F" "$GOLD/chain1d/$G1F" "[T=$T] g1.json -> chain1d/$G1F (G1)"
  s7want g2 && s7cmp "T$T G2" "$B/g2/$G1F2" "$GOLD/chain1d/$G1F2" "[T=$T] g2.json -> chain1d/$G1F2 (G2)"
  for g in $G3F DOS_reconstructed.dat $G1F; do
    s7want g3 && s7cmp "T$T G3" "$B/g3/$g" "$GOLD/chain1d/$g" "[T=$T] g3.json -> chain1d/$g (G3)"
  done
  if s7want g4; then
    s7cmp "T$T G4" "$B/g4/$G4XX.chebmom2D" "$GOLD/graphene/$G4XX.chebmom2D" "[T=$T] g4.json -> graphene/$G4XX.chebmom2D (G4)"
    s7cmp "T$T G4" "$B/g4/KuboBastin_${G4XX}JACKSON.dat" "$GOLD/graphene/KuboBastin_${G4XX}JACKSON.dat" "[T=$T] g4.json -> graphene/KuboBastin_${G4XX}JACKSON.dat (G4, 64 meV)"
    for g in sea_surf/KuboBastinSea_${G4XX}JACKSON.dat sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat \
             fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat; do
      s7cmp "T$T G4" "$B/g4/$g" "$GOLD/graphene/$g" "[T=$T] g4.json -> graphene/$g (G4, 10 meV)"
    done
  fi
  if s7want g4c; then
    s7cmp "T$T G4c" "$B/g4c/$G4CXY.chebmom2D" "$GOLD/graphene/offcentre_enclosing/$G4CXY.chebmom2D" "[T=$T] g4c.json -> graphene/offcentre_enclosing/$G4CXY.chebmom2D (G4c)"
    s7cmp "T$T G4c" "$B/g4c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/offcentre_enclosing/KuboBastin_${G4CXY}JACKSON.dat" "[T=$T] g4c.json -> graphene/offcentre_enclosing/KuboBastin_${G4CXY}JACKSON.dat (G4c)"
  fi
  for g in $G5XY.chebmom2D $G5XX.chebmom2D KuboBastin_${G5XY}JACKSON.dat KuboBastinSea_${G5XY}JACKSON.dat \
           KuboBastinSurf_${G5XY}JACKSON.dat KuboGreenWood_${G5XX}JACKSON.dat; do
    s7want g5 && s7cmp "T$T G5" "$B/g5/$g" "$GOLD/haldane/$g" "[T=$T] g5.json -> haldane/$g (G5)"
  done
  if s7want g6; then
    s7cmp "T$T G6" "$B/g6/$G6M.chebmomTD" "$GOLD/chain1d_ballistic/$G6M.chebmomTD" "[T=$T] g6.json -> chain1d_ballistic/$G6M.chebmomTD (G6)"
    for g in $G6EF; do
      s7cmp "T$T G6" "$B/g6/mean${G6M}EF${g}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${g}JACKSON.dat" "[T=$T] g6.json -> chain1d_ballistic/mean${G6M}EF${g}JACKSON.dat (G6)"
    done
  fi
  for g in SX SY SZ; do
    f=TimeProj$g-PXchain1d_magKPM_M32_stateexact.chebmomTD
    s7want g7 && s7cmp "T$T G7" "$B/g7/$f" "$GOLD/chain1d_mag/$f" "[T=$T] g7.json -> chain1d_mag/$f (G7)"
    n=g7_$(echo $g | tr 'A-Z' 'a-z')
    s7want $n && s7cmp "T$T G7" "$B/$n/$f" "$GOLD/chain1d_mag/$f" "[T=$T] $n.json -> chain1d_mag/$f (G7)"
  done
  # tutorials: every non-PNG file of the stage-0 tutorial golden, at the same relative path; the
  # capture index captured/index.tsv of the stage-0 shim (name, sha1, source path in the legacy
  # pipeline's temporary directory) is metadata, checked by sha1 against the produced captures
  for n in 01 02 03 04 05 06; do
    s7want tutorial$n || continue
    W="$B/tutorial$n"
    while IFS= read -r f; do
      case "$f" in
        captured/index.tsv) s7idx "T$T Tut$n" "$W/captured" "$GOLD/tutorials/$n/$f" "[T=$T] tutorial$n.json -> tutorials/$n/$f (Tut $n, sha1 column vs the produced captures)" ;;
        *) s7cmp "T$T Tut$n" "$W/$f" "$GOLD/tutorials/$n/$f" "[T=$T] tutorial$n.json -> tutorials/$n/$f (Tut $n)" ;;
      esac
    done < <(cd "$GOLD/tutorials/$n" && find . -type f ! -name '*.png' | sed 's|^\./||' | sort)
    while IFS= read -r f; do
      [ -f "$GOLD/tutorials/$n/$f" ] || echo "INFO      [T=$T] tutorial$n.json wrote $f (no golden at that path)"
    done < <(cd "$W" 2>/dev/null && find . -type f | sed 's|^\./||' | sort)
  done
  # G9: run files vs stage5_bases with the same arguments (byte), then the stage-5 table on the run files.
  # A run file's model.csr_dir (resolved against the repository) is passed to stage5_bases as --csr-dir;
  # at the first thread count the CSR files there are cmp'd with a fresh wannier2sparse run in $TMP.
  for m in graphene haldane; do
    for c in real bloch bloch_exact eigen real_w1 bloch_w1; do
      s7want g9_${m}_$c || continue
      C="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["model"].get("csr_dir") or "")' "$S7_RUNS/g9_${m}_$c.json")"
      [ -z "$C" ] || case "$C" in /*) ;; *) C="$REPO/$C" ;; esac
      if [ -n "$C" ] && [ ! -e "$TMP/s7w2s/$m.done" ]; then
        mkdir -p "$TMP/s7w2s/$m"; : >"$TMP/s7w2s/$m.done"
        ( cd "$GOLD/bases/$m/model" && "$S5_W2S" $m 48 48 1 VX -o "$TMP/s7w2s/$m" --quiet --no-log-file ) >"$TMP/s7w2s/$m.log" 2>&1 \
          || echo "          wannier2sparse $m failed: $(tail -1 "$TMP/s7w2s/$m.log")"
        for f in "$TMP/s7w2s/$m"/*.CSR; do
          s7cmp "T$T G9 csr_dir" "$C/$(basename "$f")" "$f" "[T=$T] g9 run-file csr_dir ${C#$REPO/}/$(basename "$f") vs fresh wannier2sparse $m 48 48 1 VX"
        done
      fi
      R="$TMP/s7ref5/T$T/$m/$c"; s7_ref5 $T $m $c "$R" "$C"
      for f in "$R"/dos/* "$R"/kuboxx/*; do
        [ -f "$f" ] || continue
        s7cmp "T$T G9" "$B/g9_${m}_$c/$(basename "$f")" "$f" "[T=$T] g9_${m}_$c.json -> stage5_bases $c $(basename "$(dirname "$f")")/$(basename "$f") (G9)"
      done
    done
    for c in real:W0/real bloch:W0/bloch bloch_exact:W0/blochx eigen:W0/eigen real_w1:W1/real bloch_w1:W1/bloch; do
      R="$TMP/s7tab/T$T/$m/${c#*:}"; mkdir -p "$R/dos" "$R/kuboxx"
      for f in "$B/g9_${m}_${c%%:*}"/*; do
        case "$(basename "$f")" in
          *SpectralOp1*) ln -sf "$f" "$R/dos/" ;;
          *NonEqOp*) ln -sf "$f" "$R/kuboxx/" ;;
        esac
      done
    done
  done
  n=0; for m in graphene haldane; do for c in real bloch bloch_exact eigen real_w1 bloch_w1; do s7want g9_${m}_$c && n=$((n + 1)); done; done
  if [ $n -eq 12 ]; then
  echo "   [T=$T] figures/stage_5_table.py on the run-file outputs"
  ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR7:$REPO/python" python3 "$REPO/figures/stage_5_table.py" "$TMP/s7tab/T$T" "$GOLD/bases" ) >"$TMP/s7_table_T$T.txt" 2>&1
  sed 's/^/          /' "$TMP/s7_table_T$T.txt"
  g="T$T G9 table"; [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  S7I[$g]=$(grep -c '^PASS' "$TMP/s7_table_T$T.txt"); S7D[$g]=$(grep -c '^FAIL' "$TMP/s7_table_T$T.txt")
  NCMP=$((NCMP + S7I[$g] + S7D[$g])); NFAIL=$((NFAIL + S7D[$g]))
  [ $((S7I[$g] + S7D[$g])) -gt 0 ] || { echo "DIFFERS   [T=$T] stage-5 table not produced"; NFAIL=$((NFAIL + 1)); }
  else
    echo "INFO      [T=$T] figures/stage_5_table.py skipped (S7_SET=$S7_SET: the table needs all six G9 configurations of both models)"
  fi
  # G10: run files vs stage6_nonorth with the run file's resolved metric (byte), then vs the orthogonal
  # goldens at the stage-6 brief's tolerance (exact 1e-10, ischo 1e-9, chebyshev 1e-8). stage6_nonorth
  # (dense exact S^-1/2 on 2048 sites: ~45 min per word on a loaded TCN) runs once per run file, at the
  # first thread count, into S7_REF6 (default $TMP/s7ref6; an existing <dir>/<run file> is reused).
  for f in "$S7_RUNS"/g10_*.json; do
    n="$(basename "$f" .json)"; s7want $n || continue
    s7skip $T $n && continue
    m="$(echo $n | cut -d_ -f2)"; v="$(echo $n | cut -d_ -f3)"
    R="$S7_REF6/$n"
    if ! ls "$B/$n"/*.chebmom* >/dev/null 2>&1; then
      g="T$T G10"; [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
      echo "DIFFERS   [T=$T] $n.json produced no moment file ($(tail -1 "$B/$n.err")); stage6_nonorth reference not run"
      NFAIL=$((NFAIL + 1)); NCMP=$((NCMP + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 ))
      continue
    fi
    # a reference pre-computed outside the gate (S7_REF6 set by the caller) is used once its log says DONE
    if [ -n "${S7_REF6_WAIT_LOG:-}" ]; then until grep -q "^DONE $n " "$S7_REF6_WAIT_LOG"; do sleep 60; done; fi
    [ -d "$R" ] || s7_ref6 $T $m $v "$B/$n/run.resolved.json" "$R"
    a="stage6_nonorth $(python3 -c 'import json,sys; x=json.load(open(sys.argv[1]))["metric"]; print("--route %s --applier %s --order %d --coeff-tol %r --tol %r --max-iter %d --max-dim %d" % (x["route"], x["applier"], x["order"], x["coeff_tol"], x["tol"], x["max_iter"], x["max_dim"]))' "$B/$n/run.resolved.json")"
    [ "$a" = "$(cat "$R.args" 2>/dev/null)" ]; s7chk $? "[T=$T] $n.json: reference arguments equal the run's resolved metric ($a)"
    for r in "$R"/*; do
      [ -f "$r" ] || continue
      s7cmp "T$T G10" "$B/$n/$(basename "$r")" "$r" "[T=$T] $n.json -> stage6_nonorth $(basename "$r") (G10)"
    done
    case "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["metric"]["applier"])' "$B/$n/run.resolved.json")" in
      exact) tol=1e-10 ;; ischo) tol=1e-9 ;; *) tol=1e-8 ;;
    esac
    if [ $m = chain1d ]; then
      s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G1F" "$GOLD/chain1d/$G1F" $tol "[T=$T] $n.json $G1F vs G1"
      s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G1F2" "$GOLD/chain1d/$G1F2" $tol "[T=$T] $n.json $G1F2 vs G2"
      s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G3F" "$GOLD/chain1d/$G3F" INFO "[T=$T] $n.json $G3F vs G3 (curve, 6 digits)"
    else
      s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G5XY.chebmom2D" "$GOLD/haldane/$G5XY.chebmom2D" $tol "[T=$T] $n.json $G5XY.chebmom2D vs G5"
      s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G5XX.chebmom2D" "$GOLD/haldane/$G5XX.chebmom2D" $tol "[T=$T] $n.json $G5XX.chebmom2D vs G5"
      for g in KuboBastin_${G5XY} KuboBastinSea_${G5XY} KuboBastinSurf_${G5XY} KuboGreenWood_${G5XX}; do
        s7rel "T$T G10 vs G1/G2/G5" "$B/$n/${g}JACKSON.dat" "$GOLD/haldane/${g}JACKSON.dat" INFO "[T=$T] $n.json ${g}JACKSON.dat vs G5 (curve, 6 digits)"
      done
    fi
  done
}

# the MCP server: run g1.json and g5.json through `run`, poll `status`, cmp the outputs
s7_mcp() {
  local R="$TMP/s7mcp" g
  mkdir -p "$R/goldens/stage0/chain1d" "$R/goldens/stage0/haldane"
  cp -r "$GOLD/chain1d/operators" "$R/goldens/stage0/chain1d/"
  cp -r "$GOLD/haldane/operators" "$R/goldens/stage0/haldane/"
  cat >"$TMP/s7_mcp_client.py" <<'PY'
import asyncio, json, os, sys, time
def show(prefix, obj):
    print(prefix + " " + json.dumps(obj, separators=(", ", ": "))); sys.stdout.flush()
async def main(root, runfiles):
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client
    params = StdioServerParameters(command=sys.executable, args=["-m", "lsquant.mcp_server", "--root", root], env=dict(os.environ))
    async with stdio_client(params) as (rd, wr):
        async with ClientSession(rd, wr) as s:
            await s.initialize()
            show("MCP <- list_tools:", [t.name for t in (await s.list_tools()).tools])
            async def call(name, args):
                show("MCP -> call_tool:", {"tool": name, "arguments": args})
                t0 = time.time(); res = await s.call_tool(name, args); dt = time.time() - t0
                text = res.content[0].text if res.content else ""
                try: body = json.loads(text)
                except ValueError: body = text
                err = bool(getattr(res, "isError", False) or getattr(res, "is_error", False))
                show("MCP <- %s (isError=%s, returned in %.3f s):" % (name, err, dt), body)
                return body, err
            for rfp in runfiles:
                name = os.path.splitext(os.path.basename(rfp))[0]
                raw = json.load(open(rfp, encoding="utf-8"))
                raw.setdefault("output", {})["dir"] = "out/" + name
                body, err = await call("run", {"runfile_json": json.dumps(raw)})
                if err or not isinstance(body, dict):
                    print("MCP RESULT %s: run refused" % name); continue
                st, _ = await call("status", {"run_id": body["run_id"]}); n = 1
                while st["status"] == "running":
                    await asyncio.sleep(2.0)
                    st = json.loads((await s.call_tool("status", {"run_id": body["run_id"]})).content[0].text); n += 1
                show("MCP <- status poll %d (final):" % n, st)
                print("MCP RESULT %s: %s" % (name, st["status"]))
            raw = json.load(open(runfiles[0], encoding="utf-8"))
            raw["model"]["dir"] = "/etc"; raw["output"]["dir"] = "out/outside"
            body, err = await call("run", {"runfile_json": json.dumps(raw)})
            refused = err or (isinstance(body, str) and "error" in body.lower())
            print("MCP RESULT confinement (model.dir=/etc): %s (isError=%s)" % ("refused" if refused else "NOT refused", err))
            print("MCP RESULT out dir created: %s" % os.path.exists(os.path.join(root, "out", "outside")))
asyncio.run(main(os.path.realpath(sys.argv[1]), sys.argv[2:]))
PY
  echo "== stage 7 MCP: python -m lsquant.mcp_server --root $R (copies of the chain1d and haldane golden operators), OMP_NUM_THREADS=16"
  ( cd "$R" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR7:$REPO/python" python3 "$TMP/s7_mcp_client.py" "$R" "$S7_RUNS/g1.json" "$S7_RUNS/g5.json" ) >"$TMP/s7_mcp.txt" 2>"$TMP/s7_mcp.err"
  echo "          (server stderr: $(wc -l <"$TMP/s7_mcp.err") lines; first error line: $(grep -m1 -E "Error|ValueError" "$TMP/s7_mcp.err" | sed "s/  */ /g"))"
  sed "s|$TMP|\$TMP|g; s/^/          /" "$TMP/s7_mcp.txt"
  { grep -q "MCP RESULT confinement (model.dir=/etc): refused" "$TMP/s7_mcp.txt" && grep -q "MCP RESULT out dir created: False" "$TMP/s7_mcp.txt"; } || { echo "FAIL      MCP confinement"; NFAIL=$((NFAIL + 1)); }
  NCMP=$((NCMP + 1))
  s7cmp "MCP" "$R/out/g1/$G1F" "$GOLD/chain1d/$G1F" "[MCP] g1.json -> chain1d/$G1F (G1)"
  for g in $G5XY.chebmom2D $G5XX.chebmom2D KuboBastin_${G5XY}JACKSON.dat KuboBastinSea_${G5XY}JACKSON.dat \
           KuboBastinSurf_${G5XY}JACKSON.dat KuboGreenWood_${G5XX}JACKSON.dat; do
    s7cmp "MCP" "$R/out/g5/$g" "$GOLD/haldane/$g" "[MCP] g5.json -> haldane/$g (G5)"
  done
}

# refusals through the CLI: non-zero exit, message on stderr, nothing written
s7_refuse() {   # <tag> <args...>
  local tag="$1" rc O="$TMP/s7ref/out_$1"; shift
  ( cd "$REPO" && PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant "$@" ) >"$TMP/s7ref/$tag.out" 2>"$TMP/s7ref/$tag.err"; rc=$?
  NCMP=$((NCMP + 1))
  if [ $rc -ne 0 ] && [ -s "$TMP/s7ref/$tag.err" ] && [ ! -e "$O" ]; then
    echo "REFUSED   [$tag] exit $rc, stderr: $(cat "$TMP/s7ref/$tag.err")"
  else
    echo "FAIL      [$tag] exit $rc, stdout: $(head -c 300 "$TMP/s7ref/$tag.out" | tr '\n' ' ') stderr: $(cat "$TMP/s7ref/$tag.err"); output dir exists: $([ -e "$O" ] && echo yes || echo no)"
    NFAIL=$((NFAIL + 1))
  fi
}

s7_refusals() {
  local R="$TMP/s7ref" f
  mkdir -p "$R"
  echo "== stage 7 refusals through the CLI (python -m lsquant run <bad>.json --out <dir>; the dir must not appear)"
  # a run file built from g1.json / g2.json / g6.json with one change; jq-free, by python
  mk() {   # <name> <base> <python expression on d>
    python3 - "$S7_RUNS/$2.json" "$R/$1.json" "$3" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
for stmt in sys.argv[3].split(";"):
    k, v = stmt.split("=", 1); keys = k.strip().split("."); t = d
    for kk in keys[:-1]: t = t[kk]
    if v.strip() == "DEL": t.pop(keys[-1], None)
    else: t[keys[-1]] = json.loads(v)
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
  }
  mk unknown_operator g1 'word="Tr[FOO delta(E-H)]"'
  mk rank_3 g2 'model.operators=["VX","VY","VZ"]; word="Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]"'
  mk mixed_basis g1 'model.operators=["VX"]; word="Tr[VX@bloch delta(E-H)]"'
  mk ut_without_time g6 'time=DEL; reconstruction=DEL; word="Tr[VX delta(E-H) U(t)^+ VX U(t)]"'
  mk unknown_key_top g1 'unknown_extra_key=1; model.dir="does/not/exist"'
  mk unknown_key_nested g1 'expansion.windw=[-2,2]; model.dir="does/not/exist"'
  mk g_plus g5 'words=DEL; reconstruction=DEL; word="Tr[VX delta(E-H) VY G+(E)]"'
  mk g_minus g5 'words=DEL; reconstruction=DEL; word="Tr[VX delta(E-H) VY G-(E)]"'
  mk dg_plus g5 'words=DEL; reconstruction=DEL; word="Tr[VX delta(E-H) VY dG+(E)/dE]"'
  mk window_eigen_csr g1 'expansion.window="eigen"'
  mk formula_unmatched g1 'reconstruction={"broadening_meV": 120, "formulas": ["bastin"]}'
  # csr_dir refusals (review 3): an empty csr_dir and one holding a 6x4x1 wannier2sparse output
  mkdir -p "$R/csr_empty"
  ( cd "$GOLD/bases/graphene/model" && "$S5_W2S" graphene 6 4 1 VX -o "$R/csr_6x4" --quiet --no-log-file ) >/dev/null 2>&1
  ( cd "$R/csr_6x4" && sha1sum * ) >"$R/csr_6x4.sha1" 2>/dev/null
  mk csr_empty g9_graphene_real "model.csr_dir=\"$R/csr_empty\""
  mk csr_mismatch g9_graphene_real "model.csr_dir=\"$R/csr_6x4\""
  for f in unknown_operator rank_3 mixed_basis ut_without_time unknown_key_top unknown_key_nested \
           g_plus g_minus dg_plus window_eigen_csr formula_unmatched csr_empty csr_mismatch; do
    s7_refuse $f run "$R/$f.json" --out "$R/out_$f"
  done
  # the csr_dir refusals must name the directory, the problem and the wannier2sparse command, and write nothing
  local e
  e="$(cat "$R/csr_empty.err")"
  case "$e" in *"$R/csr_empty"*wannier2sparse\ graphene\ 48\ 48\ 1*) [ -z "$(ls -A "$R/csr_empty")" ] ;; *) false ;; esac
  s7chk $? "csr_empty: message names the directory and the command 'wannier2sparse graphene 48 48 1 ...'; csr_empty still empty ($(ls -A "$R/csr_empty" | wc -l) files)"
  e="$(cat "$R/csr_mismatch.err")"
  case "$e" in *"$R/csr_6x4"*) case "$e" in *upercell*wannier2sparse\ graphene\ 48\ 48\ 1*|*wannier2sparse\ graphene\ 48\ 48\ 1*upercell*) ( cd "$R/csr_6x4" && sha1sum -c --quiet "$R/csr_6x4.sha1" ) >/dev/null 2>&1 && [ "$(ls "$R/csr_6x4" | wc -l)" -eq "$(wc -l <"$R/csr_6x4.sha1")" ] ;; *) false ;; esac ;; *) false ;; esac
  s7chk $? "csr_mismatch: message names the directory, the supercell mismatch and the command; csr_6x4 unchanged"
  # the same words through `words describe` (no run file): reported
  for f in "Tr[FOO delta(E-H)]" "Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]" "Tr[VX@bloch delta(E-H)]" \
           "Tr[VX delta(E-H) U(t)^+ VX U(t)]" "Tr[VX delta(E-H) VY G+(E)]" "Tr[VX delta(E-H) VY dG+(E)/dE]" \
           "Tr[VX delta(E-H) VY G-(E)]" "Tr[VX G+(E) VY delta(E-H)]"; do
    ( cd "$REPO" && PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant words describe "$f" ) >"$R/d.out" 2>"$R/d.err"
    echo "INFO      words describe \"$f\": exit $?, $( [ -s "$R/d.err" ] && cat "$R/d.err" || python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d.get("catalogue",{}); print("ACCEPTED as", c.get("shape_id"), "formulas", c.get("default_formulas"), "engine word", repr(c.get("engine_word")))' "$R/d.out")"
  done
}

stage7_setup() {
  PYDIR7="$BUILD/lsquant_core"
  ls "$PYDIR7"/_core*.so >/dev/null 2>&1 || PYDIR7="$BUILD/lsquant"
  ls "$PYDIR7"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  EXE6="$BUILD/lsquant_core/stage6_nonorth"
  [ -x "$EXE6" ] || { echo "ERROR: $EXE6 missing"; exit 2; }
  [ -n "$S7_REF6" ] || S7_REF6="$TMP/s7ref6"
  [ -x "$S5_W2S" ] || { echo "ERROR: wannier2sparse binary $S5_W2S missing (set W2S_BIN)"; exit 2; }
}

s7_totals() {
  local g
  echo "== stage 7 totals per golden (IDENTICAL or PASS / DIFFERS or FAIL)"
  for g in "${S7GROUPS[@]}"; do
    if [ "$GOLDEN_TOL" = 1 ]; then
      printf "   %-26s %4d pass (%4d byte-identical, %4d tolerance-PASS)  %4d differ\n" "$g" "${S7I[$g]:-0}" "${S7B[$g]:-0}" "${S7T[$g]:-0}" "${S7D[$g]:-0}"
    else printf "   %-26s %4d identical  %4d differ\n" "$g" "${S7I[$g]:-0}" "${S7D[$g]:-0}"; fi
  done
}

stage7() {
  local T
  stage7_setup
  echo "== stage 7: python3 -m lsquant run runs/goldens/*.json ($(ls "$S7_RUNS"/*.json | wc -l) run files, S7_SET=$S7_SET), PYTHONPATH=$PYDIR7:python, threads: $S7_THREADS"
  [ "$S7_SET" = all ] || echo "   deferred (not run in this set): $(for f in "$S7_RUNS"/*.json; do n=$(basename "$f" .json); s7want $n || printf "%s " $n; done)"
  for T in $S7_THREADS; do s7_goldens $T; done
  s7_mcp
  s7_refusals
  s7_totals
  [ -z "${S7_KEEP:-}" ] || { rm -rf "$S7_KEEP"; cp -r "$TMP/s7" "$S7_KEEP"; echo "   stage-7 outputs kept in $S7_KEEP"; }
}


# codex.md section 2.5: the deleted legacy tree. The build under test comes from this tree without
# src/, include/, test/; the legacy-final tag still builds the legacy drivers in a separate worktree
# (S7_LEGACY_WT, default /tmp/legacy-opus7), and that legacy build reproduces G1, G2, G5 (stage 0).
S7_LEGACY_WT="${S7_LEGACY_WT:-/tmp/legacy-opus7}"
s7chk() {   # <ok 0|1> <text>
  NCMP=$((NCMP + 1))
  if [ "$1" -eq 0 ]; then echo "PASS      $2"; else echo "FAIL      $2"; NFAIL=$((NFAIL + 1)); fi
}
stage7_deleted() {
  local d n t b home rc
  echo "== stage 7 section 2.5: deleted legacy tree (HEAD $(git -C "$REPO" rev-parse --short HEAD))"
  for d in src include test; do [ ! -e "$REPO/$d" ]; s7chk $? "$d/ absent from the working tree"; done
  n="$(git -C "$REPO" ls-files src include test | wc -l)"; [ "$n" -eq 0 ]; s7chk $? "git ls-files src include test: $n files"
  home="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$BUILD/CMakeCache.txt")"
  [ "$(cd "$home" && pwd -P)" = "$(cd "$REPO" && pwd -P)" ]; s7chk $? "build $BUILD configured from this tree ($home)"
  n="$(grep -cE 'inline_|kuboBastinFromChebmom|add_subdirectory\((src|test)' "$REPO/CMakeLists.txt")"; [ "$n" -eq 0 ]; s7chk $? "CMakeLists.txt: $n legacy target lines"
  [ ! -e "$REPO/scripts/regenerate_goldens.sh" ]; s7chk $? "scripts/regenerate_goldens.sh removed"
  n="$(git -C "$REPO" ls-files 'examples/*.sh' 'examples/**/*.sh' | wc -l)"; [ "$n" -eq 0 ]; s7chk $? "legacy shell scripts under examples/: $n"
  n="$(git -C "$REPO" ls-files legacy | wc -l)"; [ "$n" -eq 0 ]; s7chk $? "git ls-files legacy: $n files"
  # inline_* (legacy driver names): allowed only under goldens/stage0/**/legacy_defects and in
  # goldens/stage0/PROVENANCE.md. Reported but not counted: handoff/ (the project's record),
  # figures/out/**/*.log (recorded build and run logs) and the stage-0 harness, which must name the
  # legacy drivers to rebuild the goldens from legacy-final (scripts/check_goldens.sh stage 0,
  # figures/stage_0.py). Everything else (code, documentation) fails the check.
  local all ex bad
  all="$(git -C "$REPO" grep -l 'inline_[A-Za-z]' -- . ':!goldens/stage0/PROVENANCE.md' | grep -v '^goldens/stage0/\(.*/\)\?legacy_defects/')"
  ex="$(printf '%s\n' "$all" | grep -E '^(handoff/|figures/out/.*\.log$|scripts/check_goldens\.sh$|figures/stage_0\.py$)')"
  bad="$(printf '%s\n' "$all" | grep -vE '^(handoff/|figures/out/.*\.log$|scripts/check_goldens\.sh$|figures/stage_0\.py$)' | grep -v '^$')"
  [ -z "$bad" ]; s7chk $? "no inline_* mention outside goldens/stage0/**/legacy_defects and PROVENANCE.md (code and documentation): $(printf '%s\n' "$bad" | grep -c .) files$( [ -n "$bad" ] && echo ": $(echo $bad)")"
  [ -n "$bad" ] && for f in $bad; do git -C "$REPO" grep -n 'inline_[A-Za-z]' -- "$f" | head -3 | sed 's/^/          /'; done
  echo "INFO      inline_* in the record and the stage-0 harness (not counted): $(printf '%s\n' "$ex" | grep -c .) files: $(echo $ex | cut -c1-600)"
  t="$(git -C "$REPO" rev-parse --verify -q legacy-final^{commit})"; [ -n "$t" ]; s7chk $? "tag legacy-final exists ($t)"
  [ -n "$t" ] || return
  n="$(git -C "$REPO" ls-tree -r --name-only legacy-final -- src include test | wc -l)"; [ "$n" -gt 0 ]; s7chk $? "legacy-final holds the legacy tree ($n files under src include test)"
  git -C "$REPO" merge-base --is-ancestor legacy-final HEAD; s7chk $? "legacy-final is an ancestor of HEAD"
  rm -rf "$S7_LEGACY_WT"; git -C "$REPO" worktree prune
  git -C "$REPO" worktree add --detach "$S7_LEGACY_WT" legacy-final >"$TMP/s7legacy_wt.log" 2>&1; s7chk $? "git worktree add $S7_LEGACY_WT legacy-final"
  b="$S7_LEGACY_WT/build"
  ( cmake -S "$S7_LEGACY_WT" -B "$b" -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON && cmake --build "$b" -j 32 ) >"$TMP/s7legacy_build.log" 2>&1; rc=$?
  s7chk $rc "legacy-final: cmake -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON && cmake --build (exit $rc; $(tail -1 "$TMP/s7legacy_build.log"))"
  for d in inline_compute-kpm-nonEqOp inline_compute-kpm-spectralOp inline_kuboBastinFromChebmom \
           inline_kuboBastinSeaFromChebmom inline_kuboBastinSurfFromChebmom inline_kuboGreenwoodFromChebmom; do
    [ -x "$b/$d" ]; s7chk $? "legacy binary $d built"
  done
  echo "   legacy provenance: the stage-0 pipelines with the legacy-final build (G1, G2 at 16, 1, 64 threads; G5 at 16)"
  BUILD="$b" stage0
  git -C "$REPO" worktree remove --force "$S7_LEGACY_WT" && echo "INFO      worktree $S7_LEGACY_WT removed"
}
# words and run files live in python/lsquant; the engine knows none of it, and nothing evaluates text
stage7_grep() {
  local all hits coms pat='catalogue|runfile|grammar|Tr\[|mcp'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  all="$(cd "$REPO" && grep -rnE "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then echo "GREP      grep -rnE \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else echo "GREP      grep -rnE \"$pat\" lsquant/core lsquant/io, code lines: FOUND"; echo "$hits" | sed 's/^/          /'; NFAIL=$((NFAIL + 1)); fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
  hits="$(cd "$REPO" && grep -rnE --include='*.py' "\b(eval|exec)\(|literal_eval|exec_module" python/lsquant/words python/lsquant/runfile.py python/lsquant/cli.py python/lsquant/mcp_server.py)"
  if [ -z "$hits" ]; then echo "GREP      grep -rnE \"eval(|exec(|literal_eval|exec_module\" python/lsquant/{words,runfile.py,cli.py,mcp_server.py}: empty"
  else echo "GREP      eval/exec in the word language or run file: FOUND"; echo "$hits" | sed 's/^/          /'; NFAIL=$((NFAIL + 1)); fi
  hits="$(cd "$REPO" && grep -rnE --include='*.py' "\bexec_module\b" python/lsquant)"
  [ -z "$hits" ] || { echo "INFO      exec_module elsewhere in python/lsquant (tutorial generators imported by path):"; echo "$hits" | sed 's/^/          /'; }
}
case "$STAGE" in
  0) stage0 ;;
  1) stage1; core_grep ;;
  2) stage1; core_grep; stage2; stage2_grep ;;
  3) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep ;;
  4) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep ;;
  5) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep; stage5; stage5_grep ;;
  6) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep; stage5; stage5_grep; stage6 ;;
  6only) stage6 ;;
  5only) stage5; stage5_grep ;;
  7) stage7; core_grep; stage2_grep; stage3_grep; stage4_grep; stage5_grep; stage7_grep; stage7_deleted ;;
  7deleted) stage7_deleted ;;
  7tail) stage7_setup; s7_mcp; s7_refusals; s7_totals; core_grep; stage2_grep; stage3_grep; stage4_grep; stage5_grep; stage7_grep; stage7_deleted ;;
  7aux) stage7_setup; s7_mcp; s7_refusals; s7_totals ;;
  7runs) stage7_setup; echo "== stage 7 (run files only): S7_SET=$S7_SET${S7_LIST:+ S7_LIST=$S7_LIST}, threads: $S7_THREADS"
         for T in $S7_THREADS; do s7_goldens $T; done; s7_totals
         [ -z "${S7_KEEP:-}" ] || { rm -rf "$S7_KEEP"; cp -r "$TMP/s7" "$S7_KEEP"; echo "   stage-7 outputs kept in $S7_KEEP"; } ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$GOLDEN_TOL" = 1 ] && echo "== tolerance mode (GOLDEN_TOL=1): $NBYTE byte-identical, $NTOLP tolerance-PASS (compare() and index.tsv rows)"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
