import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np
SEED = 20260901
print("seed =", SEED)
rng = np.random.default_rng(SEED)

M, d = 4, 60

# --- check 1: beta^T A_beta^{-1} beta = S/(S+1) ---
worst = 0.0
for _ in range(200):
    b = rng.uniform(0, 0.95, M)
    A = np.outer(b, b) + np.diag(1 - b**2)
    lhs = b @ np.linalg.solve(A, b)
    S = np.sum(b**2 / (1 - b**2))
    worst = max(worst, abs(lhs - S/(S+1)))
    assert np.linalg.det(A) > 0
print("check1  max |b^T A^-1 b - S/(S+1)| =", worst)

# --- check 2: single table case  S/(S+1) = beta_k^2 ---
b = np.zeros(M); b[2] = 0.73
S = np.sum(b**2/(1-b**2))
print("check2  S/(S+1) =", S/(S+1), " beta_k^2 =", b[2]**2)

# --- check 3: row-space bound dominates every weighted perf, and equals g'G^-1 g ---
V = rng.normal(size=(M, d))
V = V / np.linalg.norm(V, axis=1, keepdims=True)      # unit rows  vhat_i
v = rng.normal(size=d); v = v / np.linalg.norm(v)      # the spike
g = V @ v
G = V @ V.T
rb_closed = g @ np.linalg.solve(G, g)
Q = V.T @ np.linalg.solve(G, V)                        # projector on row space
rb_proj = np.linalg.norm(Q @ v)**2
print("check3a rowBound projector =", rb_proj, " closed =", rb_closed,
      " diff =", abs(rb_proj - rb_closed))
viol = 0; best = 0.0
for _ in range(4000):
    w = rng.normal(size=M) * rng.choice([1.0, 0.0, 3.0], M)
    if not np.any(w != 0):
        continue
    B = np.diag(w) @ V
    ev, evec = np.linalg.eigh(B.T @ B)
    lam = ev[-1]
    # top eigenspace, with a tolerance for multiplicity
    sel = ev >= lam - 1e-9
    P = evec[:, sel] @ evec[:, sel].T
    perf = np.linalg.norm(P @ v)**2
    best = max(best, perf)
    if perf > rb_proj + 1e-8:
        viol += 1
print("check3b violations of perfW <= rowBound over 4000 random w:", viol,
      " best perfW =", best, " rowBound =", rb_proj)

# --- check 4: w = 0 breaks the bound (so the hypothesis is needed) ---
B = np.zeros((M, d))
ev, evec = np.linalg.eigh(B.T @ B)          # all zero; top eigenspace is everything
print("check4  perf at w = 0 is ||v||^2 =", float(v @ v), " > rowBound", rb_proj)
