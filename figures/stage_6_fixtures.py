#!/usr/bin/env python3
"""Tester's independent fixture verification for handoff 0006 (codex.md section 1).
Own CSR reader (legacy text layout), dense numpy/scipy only."""
import hashlib, json, os, sys
import numpy as np
import scipy.linalg as sla

REPO = os.path.expanduser("~/lsquant-unified")
G = os.path.join(REPO, "goldens/stage0")

def read_csr(path):
    with open(path) as f:
        lines = f.read().split("\n")
    n, nnz = (int(x) for x in lines[0].split()[:2])
    vals = np.array(lines[1].split(), dtype=float)
    assert vals.size == 2 * nnz, (path, vals.size, nnz)
    cols = np.array(lines[2].split(), dtype=np.int64)
    ptr = np.array(lines[3].split(), dtype=np.int64)
    assert cols.size == nnz and ptr.size == n + 1 and ptr[-1] == nnz, path
    v = vals[0::2] + 1j * vals[1::2]
    A = np.zeros((n, n), dtype=complex)
    for i in range(n):
        A[i, cols[ptr[i]:ptr[i + 1]]] += v[ptr[i]:ptr[i + 1]]
    return A

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()

def herm(A):
    return float(np.max(np.abs(A - A.conj().T)))

ok_all = True
rows = []
for model, ops in (("chain1d", ["VX"]), ("haldane", ["VX", "VY"])):
    orig = {k: read_csr(os.path.join(G, model, "operators", f"{model}.{k}.CSR")) for k in ["HAM"] + ops}
    eH = np.linalg.eigvalsh(orig["HAM"])
    for variant in ("real", "complex"):
        d = os.path.join(G, "nonorth", model, variant)
        man = json.load(open(os.path.join(d, "manifest.json")))
        r = {"fixture": f"{model}/{variant}", "N": orig["HAM"].shape[0]}
        fails = []
        # sha256
        shaok = True
        for fn, ref in man["sha256"].items():
            if sha(os.path.join(d, fn)) != ref:
                shaok = False; fails.append(f"sha256 {fn}")
        files = sorted(f for f in os.listdir(d) if f.endswith(".CSR"))
        if sorted(man["sha256"]) != files:
            shaok = False; fails.append(f"manifest sha list {sorted(man['sha256'])} != files {files}")
        r["sha256"] = "all match" if shaok else "MISMATCH"
        S = read_csr(os.path.join(d, f"{model}.S.CSR"))
        Hp = read_csr(os.path.join(d, f"{model}.HAM.CSR"))
        Vp = {k: read_csr(os.path.join(d, f"{model}.{k}.CSR")) for k in ops}
        r["herm_S"] = herm(S)
        r["herm_H'"] = herm(Hp)
        r["herm_V'"] = max(herm(Vp[k]) for k in ops)
        if r["herm_S"] > 1e-14: fails.append("S not Hermitian")
        if r["herm_H'"] > 1e-14: fails.append("H' not Hermitian to 1e-14")
        if r["herm_V'"] > 1e-14: fails.append("V' not Hermitian to 1e-14")
        S = 0.5 * (S + S.conj().T)
        eS = np.linalg.eigvalsh(S)
        r["lmin"], r["lmax"] = float(eS[0]), float(eS[-1])
        r["cond"] = float(eS[-1] / eS[0])
        r["cond_manifest"] = man["cond(S')"]
        r["cond_reldiff"] = abs(r["cond"] - man["cond(S')"]) / man["cond(S')"]
        if not eS[0] > 0: fails.append("S not positive definite")
        if r["cond_reldiff"] > 1e-6: fails.append("cond(S) != manifest to 1e-6")
        if not (5 <= r["cond"] <= 10): fails.append("cond(S) outside [5,10]")
        r["lmin_manifest_diff"] = abs(eS[0] - man["lambda_min(S')"])
        r["lmax_manifest_diff"] = abs(eS[-1] - man["lambda_max(S')"])
        Hh = 0.5 * (Hp + Hp.conj().T)
        eHp = sla.eigh(Hh, S, eigvals_only=True)
        r["spec_H"] = float(np.max(np.abs(np.sort(eHp) - eH)))
        if r["spec_H"] > 1e-10: fails.append("eigvalsh(H',S') != eigvalsh(H) to 1e-10")
        # extra (tester): same T for every operator: spectrum of a fixed combination
        cH, cV = 1.0, {"VX": 0.37, "VY": 0.21}
        Cp = cH * Hh + sum(cV[k] * 0.5 * (Vp[k] + Vp[k].conj().T) for k in ops)
        C0 = cH * orig["HAM"] + sum(cV[k] * orig[k] for k in ops)
        r["spec_comb"] = float(np.max(np.abs(np.sort(sla.eigh(Cp, S, eigvals_only=True)) - np.linalg.eigvalsh(0.5 * (C0 + C0.conj().T)))))
        if r["spec_comb"] > 1e-10: fails.append("combined H+V spectrum (extra) differs")
        r["window"] = man.get("window")
        r["verdict"] = "ACCEPTED" if not fails else "REJECTED: " + "; ".join(fails)
        ok_all &= not fails
        rows.append(r)
        print(json.dumps(r))
        sys.stdout.flush()
print("ALL_ACCEPTED" if ok_all else "SOME_REJECTED")
