import numpy as np
from draws import draws, SEED
from het import bHet, bSF, bStair, rhoHet
from rev import bRev, aHet
D = draws(20000, seed=SEED)
print("draw seed", SEED, " kept", len(D))
rows = []
for th, c, w2 in D:
    rows.append((rhoHet(th**2, c, w2), bHet(c, w2), bSF(c, w2), bStair(c, w2),
                 bRev(c, w2), float(np.sum(th**4/c))))
A = np.array(rows)
rho, bH, sf, st, rv, thr = A.T
mn = np.minimum(st, rv)
for nm, v in [("bSF", sf), ("stair", st), ("rev(ideal)", rv), ("min(stair,rev)", mn)]:
    print("%-16s coverage %5.1f%%  mean gap %.4f  beats stair on %5.1f%% of draws"
          % (nm, 100*np.mean(rho > v), np.mean((v-bH)/bH),
             100*np.mean(v < st - 1e-12)))
m = thr < 1.5
print("low-signal slice n=%d: stair %.1f%%  min(stair,rev) %.1f%%"
      % (m.sum(), 100*np.mean(rho[m] > st[m]), 100*np.mean(rho[m] > mn[m])))
