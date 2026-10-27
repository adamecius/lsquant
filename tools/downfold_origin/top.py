import sisl, numpy as np, glob, os, warnings; warnings.filterwarnings("ignore")
np.set_printoptions(precision=2,suppress=True,linewidth=200)
rows=[]
for d in sorted(glob.glob("study/*")):
    sile=sisl.get_sile(f"{d}/graphene.HSX"); H=sile.read_hamiltonian(); S=sile.read_overlap()
    orbs=H.geometry.atoms.atom[0].orbitals; no=len(orbs)
    pz=np.array([ia*no+j for ia in range(2) for j,o in enumerate(orbs) if o.l==1 and o.m==0])
    pzs=sorted(set(j for j,o in enumerate(orbs) if o.l==1 and o.m==0))
    n=12; Emax=-1e9; lmin=1e9; kmax=None
    for i in range(n):
        for j in range(n):
            k=[i/n,j/n,0]; e=H.eigh(k); w=np.linalg.eigvalsh(S.Sk(k).toarray())
            if e.max()>Emax: Emax,kmax=e.max(),k
            lmin=min(lmin,w.min())
    e,v=H.eigh([0,0,0],eigvals_only=False); Sk=H.Sk([0,0,0]).toarray(); c=np.real(np.conj(v)*(Sk@v))[pz].sum(0)
    # top state composition at Gamma
    ct=np.real(np.conj(v[:,-1])*(Sk@v[:,-1])); lab=[f"{'spd'[o.l]}{o.m}z{o.zeta}" for o in orbs]
    comp={}
    for ia in range(2):
        for jj,l in enumerate(lab): comp[l]=comp.get(l,0)+ct[ia*no+jj]
    top=sorted(comp.items(),key=lambda x:-abs(x[1]))[:3]
    rows.append((lmin,os.path.basename(d),Emax,e[c>.5],top))
for lmin,nm,Emax,pzG,top in sorted(rows,key=lambda r:-r[0]):
    print(f"{nm:>12} λmin={lmin:8.2e}  E_max(ZB)={Emax:7.2f} eV   pz@Γ: {np.array2string(pzG,precision=2)}   estado más alto en Γ: {', '.join(f'{a}:{b:+.2f}' for a,b in top)}")
