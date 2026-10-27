import numpy as np, sisl, spglib, warnings; warnings.filterwarnings("ignore")
for f in ["geom_relaxed.fdf", "relax/WTe2.XV"]:
    g = sisl.get_sile(f).read_geometry()
    cell = (g.cell, g.fxyz, [a.Z for a in g.atoms])
    for tol in [1e-4, 1e-3, 1e-2, 5e-2]:
        ds = spglib.get_symmetry_dataset(cell, symprec=tol)
        print(f"{f:22s} symprec={tol:<6} grupo: {ds.international:8s} (#{ds.number}) ops={len(ds.rotations)}")
    # best inversion centre: for each candidate centre (midpoints of same-species pairs) measure max mismatch
    xyz = g.xyz; best = (9, None)
    for i in range(g.na):
        for j in range(g.na):
            if g.atoms[i].Z != g.atoms[j].Z: continue
            c = 0.5 * (xyz[i] + xyz[j]); err = 0
            for k in range(g.na):
                r = 2 * c - xyz[k]
                d = [np.min(np.linalg.norm((r - xyz[l]) - np.round((r - xyz[l]) @ np.linalg.inv(g.cell)) @ g.cell * np.array([1, 1, 0]), axis=-1))
                     for l in range(g.na) if g.atoms[l].Z == g.atoms[k].Z]
                err = max(err, min(d))
            if err < best[0]: best = (err, c)
    print(f"   desviación máxima respecto a inversión exacta: {best[0]*1e3:.1f} mÅ")
