import numpy as np, sisl, spglib, warnings; warnings.filterwarnings("ignore")
g = sisl.get_sile("relax/WTe2.XV").read_geometry()
cell = (g.cell, g.fxyz, [a.Z for a in g.atoms])
ds = spglib.get_symmetry_dataset(cell, symprec=1e-2)
# symmetrize positions: average each atom over the images produced by all operations
R, T = ds.rotations, ds.translations
f = g.fxyz.copy(); fs = np.zeros_like(f)
for i in range(g.na):
    acc = []
    for r, t in zip(R, T):
        # find atom j such that r f_j + t = f_i (mod 1), same species
        for j in range(g.na):
            if g.atoms[j].Z != g.atoms[i].Z: continue
            d = r @ f[j] + t - f[i]; d -= np.round(d)
            if np.linalg.norm(d @ g.cell) < 2e-2:
                acc.append(f[i] + d); break
    fs[i] = np.mean(acc, axis=0)
g2 = g.copy(); g2.xyz = fs @ g.cell
for tol in [1e-5, 1e-4, 1e-3]:
    d2 = spglib.get_symmetry_dataset((g2.cell, g2.fxyz, [a.Z for a in g2.atoms]), symprec=tol)
    print(f"simetrizada, symprec={tol}: {d2.international} ops={len(d2.rotations)}")
print("desplazamiento máximo aplicado: %.2f mÅ" % (1e3 * np.abs(g2.xyz - g.xyz).max()))
sp = {74: 1, 52: 2}
geo = "LatticeConstant 1.0 Ang\n%block LatticeVectors\n" + "".join(" %14.8f %14.8f %14.8f\n" % tuple(v) for v in g2.cell) \
    + "%endblock LatticeVectors\nAtomicCoordinatesFormat Ang\n%block AtomicCoordinatesAndAtomicSpecies\n" \
    + "".join(" %14.8f %14.8f %14.8f %d\n" % (*x, sp[a.Z]) for x, a in zip(g2.xyz, g2.atoms)) \
    + "%endblock AtomicCoordinatesAndAtomicSpecies\n"
open("geom_sym.fdf", "w").write(geo)
