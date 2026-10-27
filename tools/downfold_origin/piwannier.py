"""Downfold SIESTA DZP graphene -> orthogonal 1-orbital/atom pi model (projected Wannier functions)."""
import sisl, numpy as np, warnings, sys; warnings.filterwarnings("ignore")
d = sys.argv[1] if len(sys.argv)>1 else "DZP"
N = int(sys.argv[2]) if len(sys.argv)>2 else 36
sile=sisl.get_sile(f"{d}/graphene.HSX"); H=sile.read_hamiltonian()
orbs=H.geometry.atoms.atom[0].orbitals; no=len(orbs)
odd=[ia*no+j for ia in range(2) for j,o in enumerate(orbs) if (o.l==1 and o.m==0) or (o.l==2 and abs(o.m)==1)]
trial=[ia*no+j for ia in range(2) for j,o in enumerate(orbs) if o.l==1 and o.m==0 and o.zeta==1]
odd=np.array(odd); tr=np.array([list(odd).index(t) for t in trial])
ks=[(i/N,j/N) for i in range(N) for j in range(N)]
HW={}; gapmin=1e9; detmin=1e9; leak=0
for k1,k2 in ks:
    k=[k1,k2,0]; Hk=H.Hk(k).toarray(); Sk=H.Sk(k).toarray()
    leak=max(leak,np.abs(Sk[np.ix_(odd,np.setdiff1d(range(2*no),odd))]).max(),np.abs(Hk[np.ix_(odd,np.setdiff1d(range(2*no),odd))]).max())
    h=Hk[np.ix_(odd,odd)]; s=Sk[np.ix_(odd,odd)]
    import scipy.linalg as sl
    e,v=sl.eigh(h,s); gapmin=min(gapmin,e[2]-e[1])
    A=v[:,:2].conj().T@s[:,tr]          # <psi_n|S|g_j>
    u,sv,wh=np.linalg.svd(A); U=u@wh; detmin=min(detmin,sv.min())
    HW[(k1,k2)]=U.conj().T@np.diag(e[:2])@U
# real space hoppings, R in [-N/2,N/2)
Rs=[(a,b) for a in range(-N//2,N//2) for b in range(-N//2,N//2)]
T={}
for R in Rs:
    T[R]=sum(np.exp(-2j*np.pi*(k1*R[0]+k2*R[1]))*HW[(k1,k2)] for k1,k2 in ks)/len(ks)
print(f"[{d}] fuga pi<->sigma max |H|,|S| = {leak:.1e} ; gap minimo pi*/siguiente = {gapmin:.2f} eV ; min sing(A) = {detmin:.3f}")
geo=H.geometry
xyz=geo.xyz; cell=geo.cell
rows=[]
for R,t in T.items():
    for i in range(2):
        for j in range(2):
            dist=np.linalg.norm(xyz[j]+R[0]*cell[0]+R[1]*cell[1]-xyz[i])
            rows.append((dist,t[i,j].real,abs(t[i,j].imag)))
rows.sort()
print("  hoppings (dist Ang : t eV)  [max imag %.1e]"%max(r[2] for r in rows))
shells={}
for dd,t,_ in rows:
    key=round(dd,2); shells.setdefault(key,[]).append(t)
for i,(dd,ts) in enumerate(sorted(shells.items())[:10]):
    print(f"   {dd:6.2f} Å : {np.mean(ts):+8.4f}  (n={len(ts)}, spread {np.ptp(ts):.1e})")
np.save(f"{d}/pi_hoppings.npy",np.array([(R[0],R[1],i,j,T[R][i,j].real) for R in Rs for i in range(2) for j in range(2)]))
# test: rebuild bands with truncation radius and compare to DZP pi bands
def bands_from(rcut, k):
    Hk=np.zeros((2,2),complex)
    for R,t in T.items():
        for i in range(2):
            for j in range(2):
                dist=np.linalg.norm(xyz[j]+R[0]*cell[0]+R[1]*cell[1]-xyz[i])
                if dist<=rcut+1e-6: Hk[i,j]+=t[i,j]*np.exp(2j*np.pi*(k[0]*R[0]+k[1]*R[1]))
    return np.linalg.eigvalsh(Hk)
path=sisl.BandStructure(H,[[0,0,0],[.5,0,0],[2/3,1/3,0],[0,0,0]],90)
import scipy.linalg as sl
ref=[]
for k in path.k:
    Hk=H.Hk(k).toarray()[np.ix_(odd,odd)]; Sk=H.Sk(k).toarray()[np.ix_(odd,odd)]; ref.append(sl.eigh(Hk,Sk,eigvals_only=True)[:2])
ref=np.array(ref)
print("  rcut(Å)  n_hop/átomo   err_max(eV)  err_RMS  ancho_espectro(eV)  E(K)")
for rc in [1.43,2.47,2.85,3.77,4.27,4.93,5.69,6.2,8.0,10.0,12.0]:
    E=np.array([bands_from(rc,k) for k in path.k]); err=E-ref
    nh=sum(1 for dd,_,_ in rows if 1e-3<dd<=rc+1e-6)//2
    eK=bands_from(rc,[2/3,1/3,0])
    print(f"  {rc:6.2f}   {nh:5d}   {np.abs(err).max():10.4f}  {np.sqrt((err**2).mean()):8.4f}   {E.max()-E.min():8.2f}   {np.round(eK,4)}")
