"""Monte Carlo: lamMax((1/d) sum_i w_i^2 E_i^T E_i) against bHet and the bounds."""
import numpy as np
from scipy.linalg import eigh
from het import bHet, bSF, bPart
from split import bSplit, envelope

SEED = 20260902
d = 2000
reps = 5
inst = [(np.array([0.5, 2.0]),      np.array([1.0, 0.25])),
        (np.array([0.3, 3.0]),      np.array([1.0, 0.05])),
        (np.array([0.5, 1.0, 2.0]), np.array([1.0, 0.4, 0.05]))]
print("SEED", SEED, " d =", d, " reps =", reps)
rng = np.random.default_rng(SEED)
for c, w2 in inst:
    n = np.round(c*d).astype(int)
    lam = []
    for r in range(reps):
        G = np.zeros((d, d))
        for i in range(len(c)):
            E = rng.standard_normal((n[i], d))
            G += (w2[i]/d)*(E.T @ E)
        lam.append(eigh(G, eigvals_only=True, subset_by_index=[d-1, d-1])[0])
    lam = np.array(lam)
    K4 = bSplit(c, w2, 4, seed=20260830, nrestart=25)[0]
    env = envelope(c, w2, seed=20260830)[1]
    print("c=%s w2=%s n=%s" % (np.round(c,3), np.round(w2,3), n))
    print("   lamMax mean %.4f max %.4f | bHet %.4f | bSF %.4f bPart %.4f K4 %.4f env %.4f"
          % (lam.mean(), lam.max(), bHet(c, w2), bSF(c, w2), bPart(c, w2)[0], K4, env))
