"""Same protocol, 20000 candidate draws, fast bounds only (closed forms)."""
import numpy as np, json
from het import bHet, bSF, rhoHet, bPart, bStair, optW2
SEED = 20260831
rng = np.random.default_rng(SEED)
R = []
for _ in range(20000):
    M = int(rng.integers(2, 5))
    th = rng.uniform(0.3, 2.0, M); c = rng.uniform(0.3, 3.0, M)
    thr = float(np.sum(th**4/c))
    if thr <= 1.0:
        continue
    w2 = optW2(th**2, c)
    R.append(dict(M=M, thr=thr, rho=float(rhoHet(th**2, c, w2)),
                  bHet=float(bHet(c, w2)), bSF=float(bSF(c, w2)),
                  bPart=float(bPart(c, w2)[0]), stair=float(bStair(c, w2))))
json.dump(R, open("bigsweep.json", "w"))
print("seed", SEED, "kept", len(R), "of 20000")
