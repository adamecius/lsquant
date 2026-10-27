import sisl, numpy as np, scipy.linalg as sl, warnings; warnings.filterwarnings("ignore")
for case, nocc in [("A_nopol", 46), ("C_soc", 92)]:
    H = sisl.get_sile(f"{case}/WTe2.HSX").read_hamiltonian()
    n1, n2 = 40, 70
    E = []
    for i in range(n1):
        for j in range(n2):
            k = [i/n1, j/n2, 0]
            Hk = H.Hk(k, format="array"); Sk = H.Sk(k, format="array")
            if Sk.shape[0] != Hk.shape[0]: Sk = np.kron(Sk, np.eye(2))
            E.append(sl.eigh(Hk, Sk, eigvals_only=True)[nocc-2:nocc+2])
    E = np.array(E)
    vb, cb = E[:, 1], E[:, 2]
    print(f"{case}: VBM={vb.max():+.3f}  CBM={cb.min():+.3f}  gap indirecto={cb.min()-vb.max():+.3f} eV  "
          f"gap directo mínimo={np.min(cb-vb):.3f} eV")
