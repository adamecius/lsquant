import sisl, numpy as np, warnings; warnings.filterwarnings("ignore")
H = sisl.get_sile("relax/WTe2.HSX").read_hamiltonian()
g = H.geometry
print("cell", g.cell.round(3)); print("no", H.no)
# orbital index table (1-based) from ORB_INDX
rows=[l.split() for l in open("relax/WTe2.ORB_INDX").readlines()[3:] if l.strip() and l.split()[0].isdigit() and len(l.split())>10]
print("ORB_INDX cols example:", rows[0][:10])
first={}
for r in rows[:H.no]:
    io, ia, isp, sp, iao, n, l, m, z, p = int(r[0]), int(r[1]), int(r[2]), r[3], int(r[4]), int(r[5]), int(r[6]), int(r[7]), int(r[8]), r[9]
    if z==1 and ((sp=="W" and l==2 and n==5) or (sp=="Te" and l==1 and n==5)):
        first.setdefault(ia, []).append(io)
print("trial orbitals (1-based), W 5d z1 / Te 5p z1:", first)
n=8; E=np.array([H.eigh([i/n,j/n,0]) for i in range(n) for j in range(n)])
print("nbands", E.shape[1])
for b in range(26, 72):
    print(b+1, round(E[:,b].min(),2), round(E[:,b].max(),2), end=" | ")
print()
