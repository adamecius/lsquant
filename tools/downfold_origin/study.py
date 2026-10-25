import sisl, numpy as np, glob, re, os, warnings, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
G,M,K=[0,0,0],[.5,0,0],[2/3,1/3,0]
path=[G,M,K,G]
def analyze(d):
    sile=sisl.get_sile(f"{d}/graphene.HSX"); H=sile.read_hamiltonian(); S=sile.read_overlap()
    orbs=H.geometry.atoms.atom[0].orbitals; no=len(orbs)
    pz=np.array([ia*no+j for ia in range(2) for j,o in enumerate(orbs) if o.l==1 and o.m==0])
    n=24; lm=[np.linalg.eigvalsh(S.Sk([i/n,j/n,0]).toarray()).min() for i in range(n) for j in range(n)]
    lmax=max(np.linalg.eigvalsh(S.Sk(k).toarray()).max() for k in [G,M,K])
    on=np.linalg.eigvalsh(S.tocsr(0).toarray()[:no,:no]).min()
    def ev(k):
        e,v=H.eigh(k,eigvals_only=False); Sk=H.Sk(k).toarray(); w=np.real(np.conj(v)*(Sk@v)); return e,w[pz].sum(0)
    eG,cG=ev(G); eM,cM=ev(M); piG=eG[cG>.5]; piM=eM[cM>.5]
    Kp=np.array(K); dd=1e-3; vF=np.sort(np.abs(H.eigh(Kp*(1-dd))))[:2].mean()/(np.linalg.norm(Kp@H.geometry.rcell)*dd)
    bs=sisl.BandStructure(H,path,300,["Γ","M","K","Γ"]); E=bs.apply.array.eigh()
    out=open(f"{d}/graphene.out").read(); Et=float(re.findall(r"siesta:\s+Total =\s+(-?[\d.]+)",out)[-1])
    rcs=[float(x)*0.529177 for x in re.findall(r"rc =\s+([\d.]+)",out)]
    return dict(name=os.path.basename(d),no=no,on=on,lmin=min(lm),cond=lmax/min(lm),val=-E[:,:4].min(),pi=-piG[0],
                pis=piG[1],gapM=piM[1]-piM[0],sig=eG[(cG<.5)&(eG<0)].max(),vF=vF/6.582e-16*1e-16,Et=Et,rmax=max(rcs),E=E,
                x=bs.lineark(),t=bs.lineark(True)[1])
R=[analyze(d) for d in sorted(glob.glob("study/*"))]
ref=[r for r in R if r["name"]=="TZDP_es1e-3"][0]
for r in R:
    r["rms_occ"]=np.sqrt(np.mean((r["E"][:,:4]-ref["E"][:,:4])**2))
    r["rms_pi"]=np.sqrt(np.mean((r["E"][:,4]-ref["E"][:,4])**2))   # pi* band
R.sort(key=lambda r:-r["lmin"])
hdr=f"{'base':>12} {'no':>3} {'rmax':>5} {'λ_at':>7} {'λ_min BZ':>9} {'cond':>8} {'val':>6} {'π occ':>6} {'π*(Γ)':>6} {'gapM':>5} {'σ(Γ)':>6} {'vF':>6} {'RMSocc':>7} {'RMSπ*':>6} {'Wπ':>6} {'E_tot':>9}"
print(hdr)
for r in R: print(f"{r['name']:>12} {r['no']:3d} {r['rmax']:5.2f} {r['on']:7.4f} {r['lmin']:9.2e} {r['cond']:8.1e} {r['val']:6.2f} {r['pi']:6.2f} {r['pis']:6.2f} {r['gapM']:5.2f} {r['sig']:6.2f} {r['vF']:6.3f} {r['rms_occ']:7.3f} {r['rms_pi']:6.3f} {r['pi']+r['pis']:6.2f} {r['Et']:9.3f}")
# figures
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
fam=lambda n:"SZ" if n.startswith("SZ") else ("DZ" if n.startswith("DZ") else "TZ")
col={"SZ":"#c0504d","DZ":"#1f4e79","TZ":"#2e7d32"}
for r in R:
    c=col[fam(r["name"])]
    ax[0].semilogx(r["lmin"],r["val"],"o",c=c); ax[0].annotate(r["name"],(r["lmin"],r["val"]),fontsize=7,xytext=(3,3),textcoords="offset points")
    ax[1].semilogx(r["lmin"],r["pi"],"o",c=c); ax[1].annotate(r["name"],(r["lmin"],r["pi"]),fontsize=7,xytext=(3,3),textcoords="offset points")
    ax[2].loglog(r["lmin"],max(r["rms_occ"],1e-3),"o",c=c); ax[2].annotate(r["name"],(r["lmin"],max(r["rms_occ"],1e-3)),fontsize=7,xytext=(3,3),textcoords="offset points")
ax[0].set(xlabel="λ_min de S(k) en la ZB",ylabel="ancho de valencia (eV)")
ax[1].set(xlabel="λ_min de S(k) en la ZB",ylabel="ancho π ocupada (eV)")
ax[2].set(xlabel="λ_min de S(k) en la ZB",ylabel="RMS bandas ocupadas vs TZDP_es1e-3 (eV)")
for a in ax: a.grid(alpha=.3)
fig.tight_layout(); fig.savefig("study_lmin_vs_bandwidth.png",dpi=120)
fig,ax=plt.subplots(figsize=(7,5.5)); sel=["SZ_es0.5","SZ","DZP","DZP_es1e-4","TZDP"]; cs=["#e08a00","#c0504d","#1f4e79","#7b3fa0","#2e7d32"]
for nm,c in zip(sel,cs):
    r=[q for q in R if q["name"]==nm][0]
    for i in range(r["E"].shape[1]): ax.plot(r["x"],r["E"][:,i],c=c,lw=1,label=f"{nm} (λmin={r['lmin']:.1e})" if i==0 else None)
ax.set_xticks(ref["t"]); ax.set_xticklabels(["Γ","M","K","Γ"]); [ax.axvline(t,c="0.7",lw=.5) for t in ref["t"]]
ax.axhline(0,c="0.4",ls=":",lw=.6); ax.set_ylim(-22,15); ax.set_xlim(ref["x"][0],ref["x"][-1]); ax.set_ylabel("E − E$_F$ (eV)"); ax.legend(fontsize=7,loc="lower right",frameon=False)
fig.tight_layout(); fig.savefig("study_bands.png",dpi=120)
