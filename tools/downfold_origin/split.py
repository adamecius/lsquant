import sisl, numpy as np, glob, warnings, re; warnings.filterwarnings("ignore")
def onsite_lmin(S):
    no=S.geometry.atoms.atom[0].no; M=S.tocsr(0).toarray()[:no,:no]  # R=0 block, atom 1
    return np.linalg.eigvalsh(M).min()
print(f"{'caso':>16} {'lmin atomo':>10} {'lmin(G)':>9} {'lmin(K)':>9} {'min BZ':>9} {'pi_occ':>7} {'gapM':>6} {'sigma(G)':>8} {'E_tot':>9}")
for d in ["DZP"]+sorted(glob.glob("split/*")):
    S=sisl.get_sile(f"{d}/graphene.HSX").read_overlap(); H=sisl.get_sile(f"{d}/graphene.HSX").read_hamiltonian()
    n=24; lm=min(np.linalg.eigvalsh(S.Sk([i/n,j/n,0]).toarray()).min() for i in range(n) for j in range(n))
    lK=np.linalg.eigvalsh(S.Sk([2/3,1/3,0]).toarray()).min(); lG=np.linalg.eigvalsh(S.Sk([0,0,0]).toarray()).min()
    orbs=H.geometry.atoms.atom[0].orbitals; no=len(orbs); pz=np.array([ia*no+j for ia in range(2) for j,o in enumerate(orbs) if o.l==1 and o.m==0])
    def ev(k):
        e,v=H.eigh(k,eigvals_only=False); Sk=H.Sk(k).toarray(); w=np.real(np.conj(v)*(Sk@v)); return e,w[pz].sum(0)
    eG,cG=ev([0,0,0]); eM,cM=ev([.5,0,0]); piG=eG[cG>.5]; piM=eM[cM>.5]
    sig=eG[(cG<.5)&(eG<0)].max()
    Et=float(re.findall(r"siesta:\s+Total =\s+(-?[\d.]+)",open(f"{d}/graphene.out").read())[-1])
    print(f"{d:>16} {onsite_lmin(S):10.4f} {lG:9.2e} {lK:9.2e} {lm:9.2e} {-piG[0]:7.2f} {piM[1]-piM[0]:6.2f} {sig:8.2f} {Et:9.3f}")
