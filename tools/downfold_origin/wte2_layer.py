import numpy as np
from ase.io import read, write
bulk = read("WTe2_COD2310355.cif")
print("bulk:", bulk.get_chemical_formula(), bulk.cell.cellpar().round(3))
f = bulk.get_scaled_positions(); f[:, 2] = (f[:, 2] + 0.25) % 1.0 - 0.25    # center layer at z=0
sel = np.abs(f[:, 2]) < 0.2
layer = bulk[sel]
pos = f[sel] @ bulk.cell.array; pos[:, 2] -= pos[:, 2].mean()
a, b = bulk.cell.lengths()[:2]
vac = 20.0
c = (pos[:, 2].max() - pos[:, 2].min()) + vac
layer.set_cell([a, b, c]); pos[:, 2] += c / 2; layer.set_positions(pos)
layer.wrap()
print("layer:", layer.get_chemical_formula(), "thickness %.3f A" % (pos[:,2].max()-pos[:,2].min()))
for s, p in zip(layer.get_chemical_symbols(), layer.get_positions()):
    print(f"{s:2s} {p[0]:9.4f} {p[1]:9.4f} {p[2]:9.4f}")
# distances W-Te check
from ase.neighborlist import neighbor_list
i, j, d = neighbor_list("ijd", layer, 3.0)
print("W-Te bonds <3A:", sorted(set(np.round(d[[layer[x].symbol=='W' and layer[y].symbol=='Te' for x,y in zip(i,j)]],3))))
write("WTe2_1Tp_layer.xyz", layer)
sym = {"W": 1, "Te": 2}
with open("geom.fdf", "w") as fh:
    fh.write("LatticeConstant 1.0 Ang\n%block LatticeVectors\n")
    for v in layer.cell: fh.write(" %12.6f %12.6f %12.6f\n" % tuple(v))
    fh.write("%endblock LatticeVectors\nAtomicCoordinatesFormat Ang\n%block AtomicCoordinatesAndAtomicSpecies\n")
    for s, p in zip(layer.get_chemical_symbols(), layer.get_positions()):
        fh.write(" %12.6f %12.6f %12.6f %d\n" % (*p, sym[s]))
    fh.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
