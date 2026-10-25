import sisl, numpy as np, warnings
warnings.filterwarnings("ignore")
for b in ["DZP","SZ"]:
    H = sisl.get_sile(f"{b}/graphene.HSX").read_hamiltonian()
    Ef = float(open(f"{b}/graphene.bands").readline().split()[0])
    orbs = H.geometry.atoms.atom[0].orbitals; no=len(orbs)
    pz = np.array([ia*no+j for ia in range(2) for j,o in enumerate(orbs) if o.l==1 and o.m==0])
    def ev(k):
        e,v = H.eigh(k, eigvals_only=False); S=H.Sk(k).toarray()
        w = np.real(np.conj(v)*(S@v)); return e, w[pz].sum(0)
    G,M,K = [0,0,0],[.5,0,0],[2/3,1/3,0]
    eG,cG=ev(G); eM,cM=ev(M); eK,cK=ev(K)
    bs = sisl.BandStructure(H,[G,M,K,G],600,["G","M","K","G"])
    E = bs.apply.array.eigh()
    piG=eG[cG>0.5]; piM=eM[cM>0.5]
    rc=H.geometry.rcell; d=1e-3; kk=np.array(K)*(1-d)
    e1=np.abs(np.sort(np.abs(H.eigh(kk)))[:2]).mean(); dk=np.linalg.norm(np.array(K)@rc)*d
    vF=e1/dk
    print(f"\n### {b}: {no} orb/atom, radios (Ang): " + ", ".join(f"{'spd'[o.l]}{o.zeta}:{o.R:.2f}" for o in orbs if o.m==0))
    print(f" E_F = {Ef:.3f} eV | pi en K (E-EF): {np.round(eK[cK>0.5][:2],3)}")
    print(f" pi en Gamma: {np.round(piG,2)} ; pi en M: {np.round(piM,2)}")
    print(f" Ancho de valencia (EF - fondo en Gamma) : {-E[:, :4].min():.2f} eV")
    print(f" Ancho pi ocupada (EF - pi(Gamma))       : {-piG[0]:.2f} eV")
    print(f" pi* en Gamma                             : {piG[1]:.2f} eV")
    print(f" Ancho total pi+pi* en Gamma              : {piG[1]-piG[0]:.2f} eV")
    print(f" Gap pi-pi* en M                          : {piM[1]-piM[0]:.2f} eV")
    print(f" hbar*vF = {vF:.2f} eV*Ang -> vF = {vF/6.582e-16*1e-10/1e6:.3f}e6 m/s, t = {2*vF/(3*1.42):.2f} eV")
    x=bs.lineark(); t=bs.lineark(True)[1]
    np.savez(f"{b}/bands.npz",E=E,x=x,ticks=t)
