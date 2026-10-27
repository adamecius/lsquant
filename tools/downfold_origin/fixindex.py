import csv, os
fn = os.path.expanduser("~/papers/downfolding/index.csv")
rows = list(csv.DictReader(open(fn)))
for r in rows:
    if r["key"] == "Eschrig2009_FPLO_projectiveWF":
        r["status"] = "descargado (arXiv:0905.4844; el enlace de OpenAlex apuntaba a otro paper)"
        r["pdf_url"] = "https://arxiv.org/pdf/0905.4844"
    if r["key"] == "BuongiornoNardelli2018_PAOFLOW":
        r["status"] = "sin PDF abierto (ver PAOFLOW 2.0 y 2026 abajo)"
rows += [
 dict(key="Cerasoli2021_PAOFLOW2", status="descargado (arXiv:2107.13026)", year="2021",
      title="Advanced modeling of materials with PAOFLOW 2.0: New features and software design", pdf_url="https://arxiv.org/pdf/2107.13026"),
 dict(key="PAOFLOW2026_automated_suite", status="descargado (arXiv:2609.21984)", year="2026",
      title="PAOFLOW: an automated suite for ab initio electronic, transport, and topological properties of materials", pdf_url="https://arxiv.org/pdf/2609.21984"),
]
w = csv.DictWriter(open(fn, "w", newline=""), fieldnames=list(rows[0].keys()), extrasaction="ignore")
w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['key']:36s} {r['status']}")
