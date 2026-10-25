import sys, numpy as np
sys.argv=["downfold.py","runs/Fe","Fe","--out","dbg"]
import downfold as D
# monkeypatch window_error to report where errors happen
orig=D.window_error
def we(model,kval,Eref,W,rcut):
    return orig(model,kval,Eref,W,rcut)
H,orb,ns=D.load("runs/Fe","Fe")
import argparse
p=argparse.Namespace(n2d=24,n3d=10,w=1.0,gap=0.05,smin=0.15,tol=0.01,tol_ent=0.05,jmax=40,emax_above=25.0,emax_below=30.0,nval=40,rsel=10.0,tau_grid=[1e9],mu_grid=[1],sigma_grid=[1])
# grid-point check: use kval = grid points
geom=H.geometry; N,ks=D.kgrid(geom,24,10)
spec=D.Spectrum(H,0,[("full",np.eye(H.no))],ks)
E=spec.E["full"]
W=(-1.01,1.0)
nF=((E>=W[0])&(E<=W[1])).sum(1)
print("frozen count per k: min",nF.min(),"max",nF.max())
print("states < 25 eV per k: min",(E<25).sum(1).min())
lowE=np.sort(E[:, :])[:, :14]
print("band energies at Gamma:",np.round(E[0][:16],2))
