import numpy as np, sys
def readbands(fn):
    L=open(fn).read().split()
    Ef=float(L[0]); i=1+2+2
    nb,ns,nk=int(L[i]),int(L[i+1]),int(L[i+2]); i+=3
    k=[];E=[]
    for _ in range(nk):
        k.append(float(L[i])); E.append([float(x) for x in L[i+1:i+1+nb*ns]]); i+=1+nb*ns
    nt=int(L[i]); i+=1; ticks=[]
    for _ in range(nt): ticks.append((float(L[i]),L[i+1])); i+=2
    return Ef,np.array(k),np.array(E)-Ef,ticks
for b in sys.argv[1:]:
    Ef,k,E,t=readbands(f"{b}/graphene.bands")
    tk=[np.argmin(abs(k-x)) for x,_ in t]; iG,iM,iK=tk[0],tk[1],tk[2]
    np.set_printoptions(precision=2,suppress=True,linewidth=150)
    print(f"\n### {b}  EF={Ef:.3f} eV  nbands={E.shape[1]}")
    for n,i in [("G",iG),("M",iM),("K",iK)]: print(f"  {n}:",E[i,:10])
    # pi band = band touching EF at K among occupied: the highest occupied band (index 3)
    pi=E[:,3]; pis=E[:,4]
    print(f"  Fondo valencia (G)       : {E[:,0].min():.2f} eV")
    print(f"  Banda 4 (pi) en G        : {pi[iG]:.2f} eV ; en M {pi[iM]:.2f} ; en K {pi[iK]:.3f}")
    print(f"  Banda 5 (pi*) en K {pis[iK]:.3f} ; en M {pis[iM]:.2f}")
    np.savez(f"{b}/bands.npz",k=k,E=E,tk=tk)
