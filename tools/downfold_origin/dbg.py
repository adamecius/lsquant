import sys; sys.argv=["x"]; sys.path.insert(0,".")
import downfold as D, numpy as np
H,orb,ns=D.load("/home/ICN2/jgarcia/graphene_basis/study/DZP","graphene")
SH=D.shells(orb); ZR=D.zeta_rows(orb)
for k in SH: print(k, len(SH[k]), {am:len(v) for am,v in list(ZR[k].items())[:3]}, len(ZR[k]))
