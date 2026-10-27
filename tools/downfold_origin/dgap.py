import sisl, numpy as np, warnings; warnings.filterwarnings("ignore")
H=sisl.get_sile("runs/diamond/diamond.HSX").read_hamiltonian()
n=24; best=(9,None)
for i in range(n):
  for j in range(n):
    for l in range(n):
      k=[i/n,j/n,l/n]; e=H.eigh(k); g=e[8]-e[7]
      if g<best[0]: best=(g,k,e[6:10])
print("min gap band8->9 on 24^3:",best[0],"at",best[1],best[2])
from scipy.optimize import minimize
f=lambda k: (lambda e: e[8]-e[7])(H.eigh(k))
r=minimize(f,best[1],method="Nelder-Mead",options=dict(xatol=1e-6,fatol=1e-7))
print("refined min gap:",r.fun,"at",np.round(r.x,4))
