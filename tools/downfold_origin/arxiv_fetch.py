import re, time, requests, urllib.parse, os, sys
from pypdf import PdfReader
OUT = os.path.expanduser("~/papers/downfolding")
want = {
 "Eschrig2009_FPLO_projectiveWF": "Tight-binding models for the iron-based superconductors",
 "Rayson2009_filtration": "Highly efficient method for Kohn-Sham density functional calculations of 500-10000 atom systems",
 "BuongiornoNardelli2018_PAOFLOW": "PAOFLOW: A utility to construct and operate on ab initio Hamiltonians",
 "Qian2008_QUAMBO": "Quasiatomic orbitals for ab initio tight-binding analysis",
 "Lu2004_QUAMBO": "Molecule intrinsic minimal basis sets",
}
norm = lambda t: re.sub(r"[^a-z0-9]", "", t.lower())
for key, title in want.items():
    words = [w for w in re.findall(r"[A-Za-z0-9]+", title) if len(w) > 3][:7]
    q = "+AND+".join("ti:" + urllib.parse.quote(w) for w in words)
    r = requests.get(f"https://export.arxiv.org/api/query?search_query={q}&max_results=5", timeout=60)
    entries = re.findall(r"<entry>(.*?)</entry>", r.text, re.S)
    hit = None
    for e in entries:
        t = re.sub(r"\s+", " ", re.search(r"<title>(.*?)</title>", e, re.S).group(1))
        if norm(t)[:40] == norm(title)[:40]:
            hit = (t, re.search(r"<id>http://arxiv.org/abs/(.*?)</id>", e).group(1)); break
    if not hit:
        print(f"[--] {key}: no está en arXiv (o título distinto)"); time.sleep(3); continue
    t, aid = hit
    pdf = requests.get(f"https://arxiv.org/pdf/{aid}", timeout=60, headers={"User-Agent": "Mozilla/5.0"})
    fn = os.path.join(OUT, key + ".pdf")
    if pdf.content[:5] == b"%PDF-":
        open(fn, "wb").write(pdf.content)
        first = re.sub(r"\s+", " ", PdfReader(fn).pages[0].extract_text() or "")[:90]
        print(f"[OK] {key}: arXiv:{aid}  «{t}»  | p1: {first}")
    else:
        print(f"[--] {key}: arXiv:{aid} encontrado pero sin PDF")
    time.sleep(3)
