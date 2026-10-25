import sisl, numpy as np, warnings; warnings.filterwarnings("ignore")
H=sisl.get_sile("runs/bismuthene/bismuthene.HSX").read_hamiltonian()
print(H.spin, H.no, H.spin.is_spinorbit, H.spin.is_noncolinear)
Hk=H.Hk([0.1,0.2,0],format="array"); Sk=H.Sk([0.1,0.2,0],format="array")
print("Hk",Hk.shape,"Sk",Sk.shape, "herm", np.abs(Hk-Hk.conj().T).max())
import scipy.linalg as sl
S2=np.kron(H.Sk([0.1,0.2,0],format="array")[:H.no,:H.no] if Sk.shape[0]==H.no else Sk, np.eye(1))
if Sk.shape[0]==H.no:
    Sa=np.kron(Sk,np.eye(2)); Sb=np.kron(np.eye(2),Sk)
    for n,S in [("interleaved",Sa),("blocked",Sb)]:
        try: e=sl.eigh(Hk,S,eigvals_only=True); print(n, np.round(e[:6],3))
        except Exception as x: print(n,"fail",x)
e=sl.eigh(Hk,np.kron(Sk,np.eye(2)) if Sk.shape[0]==H.no else Sk,eigvals_only=True)
Ef=float(open("runs/bismuthene/bismuthene.bands").readline().split()[0]) if False else 0
for k,n in [([0,0,0],"G"),([2/3,1/3,0],"K"),([.5,0,0],"M")]:
    hk=H.Hk(k,format="array"); sk=H.Sk(k,format="array"); sk=np.kron(sk,np.eye(2)) if sk.shape[0]==H.no else sk
    print(n, np.round(sl.eigh(hk,sk,eigvals_only=True)[10:22],3))
