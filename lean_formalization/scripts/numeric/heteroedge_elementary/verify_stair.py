import itertools, numpy as np
from scipy.optimize import linprog
from het import *
from split import fvec, total
SEED = 20260903; rng = np.random.default_rng(SEED); print("SEED", SEED)

def lpEnv(c, w2):
    M = len(c); A = []; b = []
    for r in range(1, M+1):
        for T in itertools.combinations(range(M), r):
            row = np.zeros(M); row[list(T)] = 1.0; A.append(row)
            b.append((1+np.sqrt(c[list(T)].sum()))**2)
    r = linprog(-w2, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)]*M,
                method='highs')
    return float(-r.fun)

bad = 0
for t in range(400):
    M = int(rng.integers(1, 6)); c = rng.uniform(0.05, 5, M); w2 = rng.uniform(0, 2, M)
    s = bStair(c, w2); l = lpEnv(c, w2); tot = total(c, stairSplit(c, w2))
    if abs(s-l) > 1e-9*max(1, s) or abs(s-tot) > 1e-9*max(1, s):
        bad += 1; print("MISMATCH", M, c, w2, s, l, tot)
print("greedy == LP == staircase-split value, 400 instances, mismatches:", bad)

bad2 = 0
for t in range(200):
    M = int(rng.integers(2, 6)); c = rng.uniform(0.05, 5, M); y = rng.normal(0, 2, M)
    G = rng.dirichlet(np.full(M, 0.5), 20000)
    best = np.max(G@y - fvec(c, G))
    vs = []
    for r in range(1, M+1):
        for T in itertools.combinations(range(M), r):
            bv = np.zeros(M); bv[list(T)] = 1.0/r
            vs.append(bv@y - fvec(c, bv[None, :])[0])
    if best > max(vs) + 1e-9:
        bad2 += 1
print("uniform-on-subset vertices dominate a 20000-point grid, failures:", bad2)

bad3 = 0; ratios = []
for t in range(300):
    M = int(rng.integers(1, 5)); c = rng.uniform(0.05, 4, M); w2 = rng.uniform(0.01, 2, M)
    s = bStair(c, w2); h = bHet(c, w2); sf = bSF(c, w2); bp = bPart(c, w2)[0]
    if s < h-1e-9 or s > sf+1e-9 or s > bp+1e-9:
        bad3 += 1; print("ORDER FAIL", c, w2, h, s, sf, bp)
    ratios.append((s-h)/h)
print("bHet <= bStair <= min(bSF,bPart) failures:", bad3,
      " mean (bStair-bHet)/bHet %.4f max %.4f" % (np.mean(ratios), np.max(ratios)))
