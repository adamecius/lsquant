import sisl, numpy as np, warnings; warnings.filterwarnings("ignore")
S=sisl.get_sile("DZP/graphene.HSX").read_overlap(); orbs=S.geometry.atoms.atom[0].orbitals
lab=[f"{ia+1}:{'spd'[o.l]}{['y','z','x'][o.m+1] if o.l==1 else ('' if o.l==0 else o.m)}z{o.zeta}" for ia in range(2) for o in orbs]
print("orbitales:", [(f"{'spd'[o.l]}z{o.zeta}", round(o.R,2)) for o in orbs if o.m==0])
for kn,k in [("G",[0,0,0]),("M",[.5,0,0]),("K",[2/3,1/3,0])]:
    w,v=np.linalg.eigh(S.Sk(k).toarray()); print(f"\n{kn}: 5 menores autovalores", np.array2string(w[:5],precision=6))
    c=np.abs(v[:,0])**2; idx=np.argsort(-c)[:8]
    print("   vector de λ_min:", ", ".join(f"{lab[i]}={c[i]:.2f}" for i in idx))
# check: near-zero only from d? remove d orbitals
keep=[i for i,l in enumerate(lab) if ":d" not in l]
for kn,k in [("G",[0,0,0]),("K",[2/3,1/3,0])]:
    Sk=S.Sk(k).toarray()[np.ix_(keep,keep)]; print(f"DZ sin d, {kn}: λ_min={np.linalg.eigvalsh(Sk).min():.5f}")
# raw overlap between d orbitals on same atom (should be ~ delta) and onsite S block
S0=S.Sk([0,0,0]).toarray()
