import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(6,5))
for b,c,ls in [("DZP","#1f4e79","-"),("SZ","#c0504d","--")]:
    d=np.load(f"{b}/bands.npz"); k,E,tk=d["k"],d["E"],d["tk"]
    for i in range(E.shape[1]): ax.plot(k,E[:,i],c=c,ls=ls,lw=1.2,label=b if i==0 else None)
ax.set_xticks(k[tk]); ax.set_xticklabels(["Γ","M","K","Γ"]); [ax.axvline(k[t],c="0.7",lw=.6) for t in tk]
ax.axhline(0,c="0.4",lw=.6,ls=":"); ax.set_xlim(k[0],k[-1]); ax.set_ylim(-21,12)
ax.set_ylabel("E − E$_F$ (eV)"); ax.legend(frameon=False,loc="lower right"); ax.set_title("Grafeno PBE, SIESTA 5.4.2")
fig.tight_layout(); fig.savefig("bands_DZP_vs_SZ.png",dpi=130)
