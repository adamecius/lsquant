import json, glob, os
rows=[]
for f in sorted(glob.glob("final/*/*summary.json")):
    tool = "SCDM" if "_scdm" in f else "capas"
    for r in json.load(open(f)):
        rows.append((r["tag"].replace("_scdm",""), tool, r))
print(f"{'sistema':16s} {'método':6s} {'tipo':18s} {'J':>3s} {'σ_min':>6s} {'r_cut':>6s} {'err meV':>8s} {'ancho SIESTA→modelo (eV)':>26s}")
for tag, tool, r in sorted(rows, key=lambda x:(x[0],x[1])):
    print(f"{tag:16s} {tool:6s} {r['kind']+' '+r['block'][:20]:38.38s} {r['J']:3d} {r['smin']:6.3f} {r['rcut']:6.2f} {r['err_meV']:8.1f} {r['width_siesta']:10.1f} → {r['width_model']:6.2f}")
