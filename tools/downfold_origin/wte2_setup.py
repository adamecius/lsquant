import sisl, numpy as np, os
g = sisl.get_sile("relax/WTe2.XV").read_geometry()
EF = -3.636140   # eV (relaxed run)
sp = {"W": 1, "Te": 2}
geo = "LatticeConstant 1.0 Ang\n%block LatticeVectors\n" + "".join(" %12.6f %12.6f %12.6f\n" % tuple(v) for v in g.cell) \
    + "%endblock LatticeVectors\nAtomicCoordinatesFormat Ang\n%block AtomicCoordinatesAndAtomicSpecies\n" \
    + "".join(" %12.6f %12.6f %12.6f %d\n" % (*xyz, sp[a.symbol]) for xyz, a in zip(g.xyz, g.atoms)) \
    + "%endblock AtomicCoordinatesAndAtomicSpecies\n"
open("geom_relaxed.fdf", "w").write(geo)
bands = """BandLinesScale ReciprocalLatticeVectors
%block BandLines
 1   0.0 0.0 0.0  \\Gamma
 50  0.5 0.0 0.0  X
 50  0.5 0.5 0.0  S
 50  0.0 0.5 0.0  Y
 50  0.0 0.0 0.0  \\Gamma
%endblock BandLines
SaveHS true
"""
trial = "trial-orbitals [10 11 12 13 14]\n  trial-orbitals [29 30 31 32 33]\n  trial-orbitals [41 42 43]\n  trial-orbitals [59 60 61]\n  trial-orbitals [77 78 79]\n  trial-orbitals [95 96 97]"
lo, hi = EF - 6.4, EF + 12.5
flo, fhi = EF - 6.4, EF + 1.6
w90 = "%block Wannier.Manifolds\n  mlwf\n  proj\n%endblock Wannier.Manifolds\n"
for name, nitt in [("mlwf", 200), ("proj", 0)]:
    w90 += f"""%block Wannier.Manifold.{name}
  bands 33 70
  {trial}
  spreading.nitt {nitt}
  write_hr true
  write_tb true
  window  {lo:.3f} {hi:.3f} eV
  window.frozen  {flo:.3f} {fhi:.3f} eV
%endblock Wannier.Manifold.{name}
"""
w90 += "Wannier.k [8 14 1]\nWannier.Manifolds.Unk .false.\n"
cases = {
 "A_nopol": ("", w90),
 "B_spinpol": ("Spin polarized\n%block DM.InitSpin\n 1 +1.0\n 2 +1.0\n%endblock DM.InitSpin\n", w90),
 "C_soc": ("Spin spin-orbit\n", ""),
}
for c, (spin, wan) in cases.items():
    os.makedirs(c, exist_ok=True)
    open(f"{c}/WTe2.fdf", "w").write("SystemLabel WTe2\n%include ../common.fdf\n%include ../geom_relaxed.fdf\n" + spin + bands + wan)
    os.system(f"cp pseudo/*.psml {c}/")
print(open("A_nopol/WTe2.fdf").read()[:1500])
