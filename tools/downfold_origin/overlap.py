import sisl, numpy as np, glob, warnings, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
def load(d): return sisl.get_sile(f"{d}/graphene.HSX").read_overlap()
mp = sisl.MonkhorstPack(sisl.Lattice(np.eye(3)), [1,1,1])  # placeholder
def bz_stats(S, n=24):
    ks=[[i/n,j/n,0] for i in range(n) for j in range(n)]
    lmin=[];lmax=[]
    for k in ks:
        w=np.linalg.eigvalsh(S.Sk(k).toarray()); lmin.append(w.min()); lmax.append(w.max())
    lmin=np.array(lmin); lmax=np.array(lmax); i=lmin.argmin()
    return lmin.min(), lmax.max(), ks[i]
runs=[("DZP","DZP"),("SZ","SZ")]+[(f"SZ es={d.split('es')[1]}",d) for d in sorted(glob.glob("scan/SZ_es*"),key=lambda s:float(s.split("es")[1]))]
print(f"{'base':>14} {'norb':>4} {'diag S(G)':>10} {'lmin(G)':>8} {'lmax(G)':>8} {'lmin BZ':>8} {'lmax BZ':>8} {'cond BZ':>8}  k(lmin)")
for name,d in runs:
    S=load(d); SG=S.Sk([0,0,0]).toarray(); w=np.linalg.eigvalsh(SG)
    a,b,k=bz_stats(S)
    print(f"{name:>14} {S.no:4d} {np.diag(SG).mean():10.3f} {w.min():8.4f} {w.max():8.3f} {a:8.4f} {b:8.3f} {b/a:8.1f}  {np.round(k[:2],3)}")
# figures for DZP and SZ
fig,ax=plt.subplots(2,3,figsize=(13,8))
for r,(name,d) in enumerate(runs[:2]):
    S=load(d); orbs=S.geometry.atoms.atom[0].orbitals
    lab=[f"{ia+1}:{'spd'[o.l]}{'' if o.l==0 else ['y','z','x'][o.m+1] if o.l==1 else o.m}{o.zeta}" for ia in range(2) for o in orbs]
    SG=S.Sk([0,0,0]).toarray(); SK=S.Sk([2/3,1/3,0]).toarray()
    im=ax[r,0].imshow(SG,cmap="RdBu_r",vmin=-2,vmax=2); ax[r,0].set_title(f"{name}: S(Γ)")
    ax[r,0].set_xticks(range(len(lab))); ax[r,0].set_xticklabels(lab,rotation=90,fontsize=6); ax[r,0].set_yticks(range(len(lab))); ax[r,0].set_yticklabels(lab,fontsize=6); plt.colorbar(im,ax=ax[r,0])
    for kk,kn,c in [([0,0,0],"Γ","C0"),([.5,0,0],"M","C1"),([2/3,1/3,0],"K","C2")]:
        w=np.linalg.eigvalsh(S.Sk(kk).toarray()); ax[r,1].semilogy(range(1,len(w)+1),w,"o-",ms=3,c=c,label=kn)
    ax[r,1].axhline(1,c="0.5",lw=.6,ls=":"); ax[r,1].set_title(f"{name}: espectro de S(k)"); ax[r,1].set_xlabel("índice"); ax[r,1].legend()
    n=48; L=np.zeros((n,n))
    for i in range(n):
        for j in range(n): L[i,j]=np.linalg.eigvalsh(S.Sk([i/n,j/n,0]).toarray()).min()
    im=ax[r,2].imshow(L.T,origin="lower",extent=[0,1,0,1],cmap="viridis"); ax[r,2].set_title(f"{name}: λ_min(k) en la ZB (coord. reducidas)"); plt.colorbar(im,ax=ax[r,2])
fig.tight_layout(); fig.savefig("overlap_analysis.png",dpi=120)
