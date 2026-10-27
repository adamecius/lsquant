#!/usr/bin/env python3
"""Stage-0 figure: the six tutorial figures regenerated from the legacy build.

Stage 0 (AGENTS.md section 9) proves the reference pipeline: the legacy drivers, run
exactly as the six tutorial READMEs say, reproduce the tutorial figures, and the data
those figures are drawn from match the committed goldens under goldens/stage0/tutorials.

For each tutorial this script
  1. copies the tutorial's inputs (*.py, *.json, README.md, requirements.txt, never the
     committed figures or .dat files) into a fresh work tree that mirrors the repository
     layout (work/examples/<tutorial>, work/examples/prl_style.py, work/test -> test),
     all taken from the legacy-final tag (git archive into <work>/legacy-final), never
     from the current tree,
     so the scripts' relative paths (../../test/golden, ../prl_style.py) resolve and
     examples/ itself is never written;
  2. runs every command block of the README in order (shell blocks through bash, Python
     scripts through a capture shim), with LSQUANT_BIN=<build>, PATH=<build>:$PATH,
     KPM_SEED=12345, OMP_NUM_THREADS=16 (64 where the README says so), MPLBACKEND=Agg;
  3. the capture shim wraps numpy.loadtxt and builtins.open (read mode) so that every
     .dat / .chebmom1D / .chebmom2D / .chebmomTD file a tutorial script reads for its
     figure is copied, in read order, to <work>/<tutorial>/captured/NNNN_<name>; this
     is needed because the scripts reconstruct in temporary directories and delete the
     curves after reading them. The tutorial scripts themselves are not modified;
  4. writes the PNGs to <out>/0N/ with a caption file next to each
     (<name>.caption.txt, beginning "FIG. 0.N"), giving model, size, M, seed, the
     data comparison against goldens/stage0/tutorials/0N and the section-4 stability
     result (scripts/check_goldens.sh 0).

Usage:
  python3 figures/stage_0.py [--build build] [--out figures/out/stage_0]
                             [--work DIR] [--only 1,3] [--stability-log FILE]

--stability-log reuses a saved `scripts/check_goldens.sh 0` output instead of running
it. The script never writes into goldens/.
"""
import argparse
import builtins
import filecmp
import glob
import hashlib
import os
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAPTURE_SUFFIXES = (".dat", ".chebmom1D", ".chebmom2D", ".chebmomTD")

# ---------------------------------------------------------------------------------
# README command blocks, in order. ("sh", command) runs through bash in the tutorial
# directory; ("py", argv) runs `python3 argv...` through the capture shim. Commands are
# copied from the READMEs; where a README step cannot be run as written the deviation
# is stated in "notes" and goes into the caption.
# ---------------------------------------------------------------------------------
TUTORIALS = {
    1: {
        "dir": "01_dos_1d_chain",
        "model": "1D nearest-neighbour chain, t = -1 eV, periodic",
        "size": "N = 256, 1024, 4096",
        "moments": "M = 512 (exact trace); reconstruction Jackson 10-120 meV",
        "seed": "KPM_SEED=12345 (exact trace, no random vectors)",
        "steps": [
            ("py", ["make_chain.py", "256"]),
            ("py", ["make_chain.py", "1024"]),
            ("py", ["make_chain.py", "4096"]),
            ("sh", "for N in 256 1024 4096; do "
                   "printf '{ \"mode\": \"spectral\", \"label\": \"chain1d_N%s\", \"operator\": \"1\",\\n"
                   "  \"num_moments\": 512, \"state\": \"exact_chain1d_N%s\" }\\n' $N $N > run_dos.json; "
                   "lsquant compute --config run_dos.json || exit 1; done"),
            ("py", ["-c", "import lsqplot; lsqplot.figure_moments([256, 1024, 4096], M=512)"]),
            ("sh", "lsquant reconstruct SpectralOp1chain1d_N256KPM_M512_state*.chebmom1D dos 120"),
            ("py", ["lsqplot.py"]),
        ],
        "notes": "README step 2 'Repeat for chain1d_N1024 and chain1d_N4096' run as a loop "
                 "over the three labels.",
    },
    2: {
        "dir": "02_localization_1d_chain",
        "model": "1D Anderson chain, box disorder W = 1 eV",
        "size": "N = 512 (README step 2); N = 256, 1024, 4096 x 16 realisations (lsqloc.py)",
        "moments": "M = 256, 128 times, tmax = 600, Jackson 10 meV",
        "seed": "KPM_SEED=12345 for the README CLI step; lsqloc.py sets KPM_SEED = disorder seed 1..16",
        "steps": [
            ("py", ["make_disordered_chain.py", "512", "1.0", "1"]),
            ("sh", "printf '{ \"mode\": \"msd\", \"label\": \"chain1d_dis_N512_W1\", \"operator\": \"VX\",\\n"
                   "  \"num_moments\": 256, \"num_times\": 128, \"tmax\": 600 }\\n' > run_msd.json && "
                   "lsquant compute --config run_msd.json && "
                   "inline_timeCorrelationsFromChebmom Correlation*chain1d_dis_N512_W1*chebmomTD 10 0.0"),
            ("py", ["lsqloc.py"]),
        ],
        "notes": "",
    },
    3: {
        "dir": "03_haldane_hall",
        "model": "Haldane model t1 = -1, t2 = 0.15, phi = pi/2 (test/golden/haldane operators)",
        "size": "dim = 512 (256 cells), bounds [-4.5, 4.5] eV from the .desc",
        "moments": "M = 128 x 128 (exact trace), Jackson 10 meV",
        "seed": "KPM_SEED=12345 (exact trace)",
        "steps": [
            ("py", ["make_haldane.py"]),
            ("sh", "lsquant inspect operators/haldane.HAM.desc"),
            ("sh", "printf '{ \"label\": \"haldane\", \"operator_right\": \"VX\", \"operator_left\": \"VY\",\\n"
                   "  \"num_moments\": 128, \"kernel\": \"jackson\", \"state\": \"exact\" }\\n' > run_hall.json && "
                   "printf '{ \"label\": \"haldane\", \"operator_right\": \"VX\", \"operator_left\": \"VX\",\\n"
                   "  \"num_moments\": 128, \"kernel\": \"jackson\", \"state\": \"exact\" }\\n' > run_long.json && "
                   "lsquant compute --config run_hall.json && lsquant compute --config run_long.json"),
            ("sh", "lsquant reconstruct NonEqOpVX-VYhaldane*.chebmom2D bastin 10 && "
                   "inline_kuboBastinSeaFromChebmom  NonEqOpVX-VYhaldane*.chebmom2D 10 && "
                   "inline_kuboBastinSurfFromChebmom NonEqOpVX-VYhaldane*.chebmom2D 10 && "
                   "lsquant reconstruct NonEqOpVX-VXhaldane*.chebmom2D greenwood 10"),
            ("py", ["lsqhall.py"]),
            ("py", ["haldane_kspace.py"]),
        ],
        "notes": "",
    },
    4: {
        "dir": "04_spin_precession_magnetic_chain",
        "model": "magnetic chain, on-site exchange J_ex = 0.1 eV (test/models/chain1d_mag)",
        "size": "chain1d_mag operators (HAM, SX, SY, SZ, PX)",
        "moments": "M = 32, 96 times, tmax = 95 fs (exact trace)",
        "seed": "KPM_SEED=12345 (exact trace)",
        "steps": [
            ("stage_mag", None),
            ("sh", "ls operators"),
            ("sh", "for OP in SX SY SZ; do "
                   "inline_compute-kpm-TimeEvProjetedOp chain1d_mag $OP PX 32 96 95 exact || exit 1; done"),
            ("py", ["lsqspin.py"]),
        ],
        "notes": "The README says the operators are 'shipped here, staged under operators/', but "
                 "examples/04 has no operators/ directory and no 'exact' state file; they are "
                 "staged from test/models/chain1d_mag/operators and an exact local-state file is "
                 "written, as test/spin_precession.sh does.",
    },
    5: {
        "dir": "05_charge_transport_square",
        "model": "2D square Anderson lattice, W = 0.5, 1, 2, 4 (README step: W = 1)",
        "size": "L = 200 (40000 sites)",
        "moments": "MSD M = 256, 160 times, tmax = 150; DOS M = 1024; Jackson 80 meV at E = 2",
        "seed": "KPM_SEED=12345 for the README CLI step; lsqtransport.py sets KPM_SEED=1",
        "steps": [
            ("py", ["make_square.py", "200", "1", "1"]),
            ("sh", "printf '{ \"mode\": \"msd\", \"label\": \"square_L200_W1\", \"operator\": \"VX\",\\n"
                   "  \"num_moments\": 256, \"num_times\": 160, \"tmax\": 150 }\\n' > run_msd.json && "
                   "lsquant compute --config run_msd.json && "
                   "inline_timeCorrelationsFromChebmom Correlation*square_L200_W1*chebmomTD 80 2"),
            ("py", ["lsqtransport.py"]),
        ],
        "notes": "README step 3 ('lsquant compute --config dos.json; lsquant reconstruct "
                 "SpectralOp1*.chebmom1D dos 80') references a dos.json the README never "
                 "defines; it is not run standalone. lsqtransport.py writes and runs that "
                 "dos.json itself.",
    },
    6: {
        "dir": "06_spin_relaxation_square",
        "model": "2D square lattice: Elliott-Yafet random spin-flip field; Rashba + Anderson (DP)",
        "size": "L = 20 (EY, DP); rate law L = 50 (N = 5000), exact cross-check L = 8",
        "moments": "EY M = 128, 80 times, tmax 260; DP M = 160, 120 times, tmax 45; Jackson 140 meV at E = 2",
        "seed": "model seed 1; rp_spin.py sets KPM_SEED = 1..20 for the random-phase trace",
        "steps": [
            ("py", ["make_ey.py", "20", "0.1", "1"]),
            ("py", ["lsqrelax.py"]),
            ("py", ["make_rashba.py", "20", "0.2", "3", "1"]),
            ("py64", ["rp_validate.py"]),
            ("py", ["plot_ratelaw.py"]),
        ],
        "notes": "rp_validate.py run with OMP_NUM_THREADS=64 as the README says.",
    },
}


# ---------------------------------------------------------------------------------
# capture shim: runs inside the tutorial's Python process
# ---------------------------------------------------------------------------------
def run_shim(capdir, argv):
    """Execute a tutorial script (or `-c code`) with read capture enabled."""
    import runpy
    import numpy as np

    os.makedirs(capdir, exist_ok=True)
    seen = set()
    index = os.path.join(capdir, "index.tsv")
    real_open = builtins.open
    real_loadtxt = np.loadtxt

    def capture(path):
        try:
            p = os.fspath(path)
        except TypeError:
            return
        if not isinstance(p, str) or not p.endswith(CAPTURE_SUFFIXES) or not os.path.isfile(p):
            return
        with real_open(p, "rb") as fh:
            digest = hashlib.sha1(fh.read()).hexdigest()
        key = (os.path.abspath(p), digest)
        if key in seen:
            return
        seen.add(key)
        n = len([f for f in os.listdir(capdir) if f != "index.tsv"])
        dst = os.path.join(capdir, "%04d_%s" % (n, os.path.basename(p)))
        shutil.copyfile(p, dst)
        with real_open(index, "a") as fh:
            ap = os.path.abspath(p)
            src = (os.path.relpath(ap) if ap.startswith(os.getcwd() + os.sep)
                   else "<tmpdir>/" + os.path.basename(p))
            fh.write("%s\t%s\t%s\n" % (os.path.basename(dst), digest, src))

    def open_capture(file, mode="r", *a, **k):
        if isinstance(mode, str) and "r" in mode and "+" not in mode:
            capture(file)
        return real_open(file, mode, *a, **k)

    def loadtxt_capture(fname, *a, **k):
        capture(fname)
        return real_loadtxt(fname, *a, **k)

    builtins.open = open_capture
    np.loadtxt = loadtxt_capture
    sys.path.insert(0, os.getcwd())
    if argv[0] == "-c":
        sys.argv = ["-c"] + argv[2:]
        exec(compile(argv[1], "<string>", "exec"), {"__name__": "__main__"})
    else:
        sys.argv = list(argv)
        runpy.run_path(argv[0], run_name="__main__")


# ---------------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------------
LEGACY_TAG = "legacy-final"
LEGACY = None  # extracted legacy-final tree (examples/, test/, utilities/); set in main()


def extract_legacy(work):
    """Extract examples/, test/ and utilities/ of the legacy-final tag into <work>/legacy-final.

    Stage 0 is legacy provenance: the tutorial scripts, prl_style.py and the test/ inputs
    they read are taken from the tag, never from the current tree (which no longer holds
    the legacy drivers). `git archive` leaves the repository and its worktrees untouched."""
    dst = os.path.join(work, LEGACY_TAG)
    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    arch = subprocess.run(["git", "-C", REPO, "archive", LEGACY_TAG, "examples", "test", "utilities"],
                          stdout=subprocess.PIPE, check=True).stdout
    subprocess.run(["tar", "-x", "-C", dst], input=arch, check=True)
    return dst


def prepare(work, tut):
    """Fresh work tree mirroring the legacy-final layout around one tutorial."""
    src = os.path.join(LEGACY, "examples", tut)
    dst = os.path.join(work, "examples", tut)
    if os.path.exists(os.path.join(work, "examples")):
        shutil.rmtree(os.path.join(work, "examples", tut), ignore_errors=True)
    os.makedirs(dst)
    for f in os.listdir(src):
        if f.endswith((".py", ".json", ".md", ".txt")) and os.path.isfile(os.path.join(src, f)):
            shutil.copy(os.path.join(src, f), dst)
    shutil.copy(os.path.join(LEGACY, "examples", "prl_style.py"), os.path.join(work, "examples"))
    for link in ("test", "utilities"):
        lp = os.path.join(work, link)
        if os.path.lexists(lp):
            os.remove(lp)
        os.symlink(os.path.join(LEGACY, link), lp)
    return dst


def stage_mag(cwd):
    opd = os.path.join(LEGACY, "test", "models", "chain1d_mag", "operators")
    os.makedirs(os.path.join(cwd, "operators"), exist_ok=True)
    for op in ("HAM", "SX", "SY", "SZ", "PX"):
        shutil.copy(os.path.join(opd, "chain1d_mag.%s.CSR" % op), os.path.join(cwd, "operators"))
    n = int(open(os.path.join(opd, "chain1d_mag.HAM.CSR")).readline().split()[0])
    with open(os.path.join(cwd, "exact"), "w") as f:
        f.write("local\n%d\n%s\n" % (n, "\n".join(map(str, range(n)))))


def run_tutorial(num, spec, build, work, log):
    tut = spec["dir"]
    cwd = prepare(work, tut)
    capdir = os.path.join(cwd, "captured")
    env = dict(os.environ, LSQUANT_BIN=build, KPM_SEED="12345", OMP_NUM_THREADS="16",
               MPLBACKEND="Agg", PATH=build + os.pathsep + os.environ.get("PATH", ""))
    failures = []
    t0 = time.time()
    for kind, arg in spec["steps"]:
        if kind == "stage_mag":
            stage_mag(cwd)
            log.write("[%s] staged test/models/chain1d_mag/operators + exact\n" % tut)
            continue
        e = dict(env)
        if kind == "sh":
            argv, shown = ["bash", "-c", arg], arg
        else:
            if kind == "py64":
                e["OMP_NUM_THREADS"] = "64"
            argv = [sys.executable, os.path.abspath(__file__), "--shim", capdir] + arg
            shown = "python3 " + " ".join(arg)
        ts = time.time()
        log.write("[%s] $ %s\n" % (tut, shown)); log.flush()
        r = subprocess.run(argv, cwd=cwd, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log.write(r.stdout)
        log.write("[%s] exit %d (%.1f s)\n" % (tut, r.returncode, time.time() - ts)); log.flush()
        if r.returncode != 0:
            failures.append("`%s` exited %d: %s" % (shown, r.returncode,
                                                     (r.stdout.strip().splitlines() or ["(no output)"])[-1]))
    return cwd, failures, time.time() - t0


def collect(cwd):
    """Figure data produced in the tutorial dir: captured reads + README-step .dat files."""
    data = {}
    capdir = os.path.join(cwd, "captured")
    if os.path.isdir(capdir):
        for f in sorted(os.listdir(capdir)):
            data["captured/" + f] = os.path.join(capdir, f)
    for f in sorted(glob.glob(os.path.join(cwd, "*.dat"))):
        data[os.path.basename(f)] = f
    pngs = sorted(glob.glob(os.path.join(cwd, "fig_*.png")))
    return data, pngs


def compare_with_goldens(num, data, pngs):
    gdir = os.path.join(REPO, "goldens", "stage0", "tutorials", "%02d" % num)
    if not os.path.isdir(gdir):
        return "no goldens/stage0/tutorials/%02d to compare with" % num, []
    golden = set()
    for root, _, files in os.walk(gdir):
        for f in files:
            golden.add(os.path.relpath(os.path.join(root, f), gdir))
    lines, nid, ndiff = [], 0, 0
    for name, path in sorted(data.items()):
        g = os.path.join(gdir, name)
        if name not in golden:
            lines.append("DIFFERS   %s (not in golden)" % name); ndiff += 1
        elif filecmp.cmp(path, g, shallow=False):
            nid += 1
        else:
            out = subprocess.run([sys.executable, os.path.join(REPO, "scripts", "golden_diff.py"), path, g],
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout.strip()
            lines.append("DIFFERS   %s: %s" % (name, out)); ndiff += 1
    produced = set(data) | {os.path.basename(p) for p in pngs}
    for name in sorted(golden - produced):
        if not name.endswith(".png"):
            lines.append("DIFFERS   %s (golden not produced)" % name); ndiff += 1
    for p in pngs:
        g = os.path.join(gdir, os.path.basename(p))
        tag = "IDENTICAL" if os.path.isfile(g) and filecmp.cmp(p, g, shallow=False) else "differs"
        lines.append("png %s %s" % (os.path.basename(p), tag))
    summary = "figure data vs goldens/stage0/tutorials/%02d: %d IDENTICAL, %d DIFFERS" % (num, nid, ndiff)
    return summary, lines


def stability(build, logfile):
    if logfile:
        text = open(logfile).read()
    else:
        r = subprocess.run(["bash", os.path.join(REPO, "scripts", "check_goldens.sh"), "0", build],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=REPO)
        text = r.stdout
    lines = [l for l in text.splitlines() if l.startswith(("IDENTICAL", "DIFFERS", "GATE", "=="))]
    nid = sum(l.startswith("IDENTICAL") for l in lines)
    nd = sum(l.startswith("DIFFERS") for l in lines)
    verdict = next((l for l in lines if l.startswith("GATE")), "GATE result missing")
    return text, "scripts/check_goldens.sh 0: %d IDENTICAL, %d DIFFERS (G1, G2, G5 at 16 threads; " \
                 "G1, G2 at 1 and 64 threads) -> %s" % (nid, nd, verdict)


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--shim":
        run_shim(sys.argv[2], sys.argv[3:])
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_0"))
    ap.add_argument("--work", default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "lsquant_stage_0_work"))
    ap.add_argument("--only", default="1,2,3,4,5,6")
    ap.add_argument("--stability-log", default=None)
    a = ap.parse_args()
    build, out, work = (os.path.abspath(p) for p in (a.build, a.out, a.work))
    os.makedirs(out, exist_ok=True)
    os.makedirs(work, exist_ok=True)
    global LEGACY
    LEGACY = extract_legacy(work)

    stab_text, stab = stability(build, a.stability_log)
    with open(os.path.join(out, "check_goldens_stage0.log"), "w") as f:
        f.write(stab_text)
    print(stab)

    head = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"],
                          stdout=subprocess.PIPE, text=True).stdout.strip()
    status = 0
    with open(os.path.join(out, "run.log"), "w") as log:
        for num in [int(x) for x in a.only.split(",")]:
            spec = TUTORIALS[num]
            cwd, failures, wall = run_tutorial(num, spec, build, work, log)
            data, pngs = collect(cwd)
            summary, detail = compare_with_goldens(num, data, pngs)
            odir = os.path.join(out, "%02d" % num)
            shutil.rmtree(odir, ignore_errors=True)
            os.makedirs(odir)
            for i, p in enumerate(pngs):
                shutil.copy(p, odir)
                name = os.path.basename(p)
                cap = ["FIG. 0.%d%s (%s). Tutorial %02d (%s), regenerated by the legacy build "
                       "(HEAD %s) through the README pipeline." % (num, "abcdefgh"[i] if len(pngs) > 1 else "",
                                                                    name, num, spec["dir"], head),
                       "Model: %s." % spec["model"],
                       "Size: %s." % spec["size"],
                       "Moments: %s." % spec["moments"],
                       "Seed: %s." % spec["seed"],
                       "Data: %s." % summary,
                       "Stability (section 4): %s." % stab]
                if spec["notes"]:
                    cap.append("Note: %s" % spec["notes"])
                if failures:
                    cap.append("FAILED steps: " + "; ".join(failures))
                with open(os.path.join(odir, name[:-4] + ".caption.txt"), "w") as f:
                    f.write("\n".join(cap) + "\n")
            with open(os.path.join(odir, "compare.txt"), "w") as f:
                f.write(summary + "\n" + "\n".join(detail) + "\n")
            print("tutorial %02d: %d figures, %d data files, %.0f s; %s%s"
                  % (num, len(pngs), len(data), wall, summary,
                     "" if not failures else "; FAILED: " + "; ".join(failures)))
            if failures:
                status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
