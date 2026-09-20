# Monte Carlo: are the columns of the component-weighted stacksvd estimate asymptotically
# orthogonal (main_paper.tex:926)? Model: M=2 tables, r=2 shared components, R_i = I,
# X_i = sqrt(d)^{-1}... we use the paper's scaling X_i = U_i Theta_i V^T + Z_i/sqrt(d), n_i = c_i d.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np, sys
seed = 20260905
rng = np.random.default_rng(seed)
print("seed", seed)
def run(d, theta, c, W, reps):
    # theta: M x r, W: M x r weights (w_ij for stack j), returns mean <vhat_1, vhat_2>^2 and alignments
    M, r = theta.shape
    n = [int(ci * d) for ci in c]
    out = []; al = []
    for _ in range(reps):
        V, _ = np.linalg.qr(rng.standard_normal((d, r)))
        X = []
        for i in range(M):
            U, _ = np.linalg.qr(rng.standard_normal((n[i], r)))
            Z = rng.standard_normal((n[i], d)) / np.sqrt(d)
            X.append(U @ np.diag(theta[i]) @ V.T + Z)
        vh = []
        for j in range(r):
            S = np.vstack([W[i, j] * X[i] for i in range(M)])
            G = S.T @ S
            evals, evecs = np.linalg.eigh(G)
            vh.append(evecs[:, -1 - j])   # j-th largest (ordered case: ell_j = j)
        out.append((vh[0] @ vh[1]) ** 2)
        al.append([(vh[j] @ V[:, j]) ** 2 for j in range(r)])
    return np.mean(out), np.std(out) / np.sqrt(reps), np.mean(al, axis=0)
theta = np.array([[3.0, 1.5], [2.0, 1.2]]); c = (1.0, 1.0)
for d in (200, 400, 800):
    same = np.ones((2, 2))
    m, se, al = run(d, theta, c, same, 20)
    print(f"d={d} equal weights: <v1,v2>^2 = {m:.2e} (exact 0 expected), align {al}")
    W = np.array([[1.0, 0.3], [0.4, 1.0]])   # different weights per component
    m, se, al = run(d, theta, c, W, 20)
    print(f"d={d} weights {W.tolist()}: <v1,v2>^2 = {m:.4f} +- {se:.4f}, align {al}")
