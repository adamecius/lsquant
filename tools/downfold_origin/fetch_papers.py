#!/usr/bin/env python
"""Look up the downfolding / Wannier references in OpenAlex and download open-access PDFs."""
import csv, os, re, sys, time
import requests
from pyalex import Works

OUT = os.path.expanduser("~/papers/downfolding")
REFS = [  # (short key, DOI or None, title for search)
    ("Qiao2023_PDWF", "10.1038/s41524-023-01146-w", "Projectability disentanglement for accurate and automated electronic-structure Hamiltonians"),
    ("Jiang2025_PDWF_SOC", None, "Robust Wannierization including magnetization and spin-orbit coupling via projectability disentanglement"),
    ("Damle2015_SCDM", "10.1021/ct500985f", "Compressed representation of Kohn-Sham orbitals via selected columns of the density matrix"),
    ("Damle2018_SCDM_entangled", "10.1137/17M1129696", "Disentanglement via entanglement: a unified method for Wannier localization"),
    ("Vitale2020_autoSCDM", "10.1038/s41524-020-0312-y", "Automated high-throughput Wannierisation"),
    ("Thygesen2005_POWF_PRL", "10.1103/PhysRevLett.94.026405", "Partly occupied Wannier functions"),
    ("Thygesen2005_POWF_PRB", "10.1103/PhysRevB.72.125119", "Partly occupied Wannier functions: Construction and applications"),
    ("Souza2001_disentanglement", "10.1103/PhysRevB.65.035109", "Maximally localized Wannier functions for entangled energy bands"),
    ("Agapito2013_PAO_projection", "10.1103/PhysRevB.88.165127", "Effective and accurate representation of extended Bloch states on finite Hilbert spaces"),
    ("Agapito2016_PAO_2D", "10.1103/PhysRevB.93.125137", "Accurate tight-binding Hamiltonians for two-dimensional and layered materials"),
    ("BuongiornoNardelli2018_PAOFLOW", "10.1016/j.commatsci.2017.11.034", "PAOFLOW: A utility to construct and operate on ab initio Hamiltonians from the projections of electronic wavefunctions on atomic orbital bases"),
    ("Eschrig2009_FPLO_projectiveWF", "10.1103/PhysRevB.79.104503", "Tight-binding models for the iron-based superconductors"),
    ("Korytar2010_SIESTA_Wannier", "10.1088/0953-8984/22/38/385601", "Band selection and disentanglement using maximally localized Wannier functions: the cases of Co impurities in bulk copper and the Cu(111) surface"),
    ("Rayson2009_filtration", "10.1103/PhysRevB.80.205104", "Highly efficient method for Kohn-Sham density functional calculations of 500-10000 atom systems"),
    ("Nakata2015_CONQUEST_MSSF", "10.1039/c5cp00934k", "Optimized multi-site local orbitals in the large-scale DFT program CONQUEST"),
    ("Lu2004_QUAMBO", "10.1103/PhysRevB.70.041101", "Molecule intrinsic minimal basis sets"),
    ("Qian2008_QUAMBO", "10.1103/PhysRevB.78.245112", "Quasiatomic orbitals for ab initio tight-binding analysis"),
    ("Andersen2000_NMTO", "10.1103/PhysRevB.62.R16219", "Muffin-tin orbitals of arbitrary order"),
    ("Fan2021_linear_scaling_transport", "10.1016/j.physrep.2020.12.001", "Linear scaling quantum transport methodologies"),
]

HEAD = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) papers-fetch/1.0"}


def norm(t):
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def lookup(doi, title):
    if doi:
        try:
            return Works()[f"https://doi.org/{doi}"], "doi"
        except Exception:
            pass
    res = Works().search(title).get(per_page=5)
    for w in res:
        if norm(w.get("title"))[:60] == norm(title)[:60]:
            return w, "título"
    return (res[0], "título (aprox.)") if res else (None, "no encontrado")


def pdf_candidates(w):
    urls = []
    for loc in [w.get("best_oa_location")] + (w.get("locations") or []):
        if loc and loc.get("pdf_url"):
            urls.append(loc["pdf_url"])
    for loc in w.get("locations") or []:            # arXiv landing -> pdf
        lp = (loc or {}).get("landing_page_url") or ""
        m = re.search(r"arxiv\.org/abs/([\w.\-/]+)", lp)
        if m:
            urls.append(f"https://arxiv.org/pdf/{m.group(1)}")
    seen = []
    return [u for u in urls if not (u in seen or seen.append(u))]


def download(url, fn):
    r = requests.get(url, headers=HEAD, timeout=60, allow_redirects=True)
    if r.status_code == 200 and r.content[:5] == b"%PDF-":
        open(fn, "wb").write(r.content)
        return len(r.content)
    return None


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for key, doi, title in REFS:
        w, how = lookup(doi, title)
        if w is None:
            rows.append(dict(key=key, status="no encontrado en OpenAlex", doi=doi or "", title=title))
            print(f"[--] {key}: no encontrado"); continue
        oa = w.get("open_access") or {}
        info = dict(key=key, title=w.get("title"), year=w.get("publication_year"),
                    doi=(w.get("doi") or "").replace("https://doi.org/", ""), openalex=w.get("id"),
                    oa_status=oa.get("oa_status"), cited_by=w.get("cited_by_count"), found_by=how)
        fn = os.path.join(OUT, f"{key}.pdf")
        got = None
        for u in pdf_candidates(w):
            try:
                n = download(u, fn)
            except Exception:
                n = None
            if n:
                got = u; info["status"] = f"descargado ({n/1e6:.1f} MB)"; info["pdf_url"] = u
                break
            time.sleep(1)
        if not got:
            info["status"] = "sin PDF abierto"
            info["pdf_url"] = ""
        rows.append(info)
        print(f"[{'OK' if got else '--'}] {key}: {info['status']}  (OA={info['oa_status']}, {how})")
        time.sleep(0.5)
    with open(os.path.join(OUT, "index.csv"), "w", newline="") as fh:
        keys = ["key", "status", "year", "title", "doi", "oa_status", "cited_by", "pdf_url", "openalex", "found_by"]
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
