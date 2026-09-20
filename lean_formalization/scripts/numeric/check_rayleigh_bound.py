"""Numeric check of the Rayleigh-quotient (trace) upper bound for weighted SVDstack.

Claims tested (see notes/archive/external_audit_v3_review_verdicts_2026-09-01.md, notes/paper_edits.md E1, E2):

  C1 (rank 1). For every weight vector w (including data-dependent and singular w), the
      finite-N performance ||v_hat(w)^T v||^2 is at most g^T G^{-1} g, where
      g = V_tilde v and G = V_tilde V_tilde^T. That bound does not depend on w and converges
      to S/(S+1) with S = sum beta_i^2 / (1 - beta_i^2).
  C2 (rank r, projector form). For every weight matrix W (block diagonal, one block per
      table), ||V_hat(W)^T V||_F^2 <= tr(g^T G^{-1} g) with g = V_tilde V, and that bound
      converges to L_opt = tr(B_R^T A^{-1} B_R), the paper's optimal value.
  C3 (audit's counterexample). In the audit's example (R_1 = e_1 with a subcritical theta,
      R_2 = R_3 = e_2 with supercritical theta, so rank B_R = 1 < r = 2), the paper's
      removal remark is a proof gap and not a false value: the estimator at W_opt still
      attains L_opt, with no rank condition on B_R.

Run:  python3 scripts/check_rayleigh_bound.py
"""
import sys

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np

SEED = 20260901
rng = np.random.default_rng(SEED)
print(f"seed {SEED}")


def spiked_table(theta_vec, R, n, d, rng):
    """X = sum_k theta_k u_k v_k^T + E / sqrt(d) with v_k = e_k (aligned columns), E iid N(0,1).

    theta_vec: length r, theta_k = 0 for components absent from R.
    Returns X (n x d) with the paper's normalization X = U Theta V^T + E, E_ij ~ N(0, 1/d).
    """
    r = len(theta_vec)
    U = np.linalg.qr(rng.standard_normal((n, r)))[0]
    V = np.eye(d)[:, :r]
    E = rng.standard_normal((n, d)) / np.sqrt(d)
    return U @ np.diag(theta_vec) @ V.T + E


def beta_sq(theta, c):
    """Squared right-singular-vector correlation, paper's prop:single_table."""
    if theta ** 4 <= c:
        return 0.0
    return (theta ** 4 - c) / (theta ** 4 + theta ** 2)


def top_right_vectors(X, k):
    _, s, Vt = np.linalg.svd(X, full_matrices=False)
    return Vt[:k].T  # d x k


# ---------------------------------------------------------------- C1: rank 1
print("\n== C1: rank 1, bound uniform in w ==")
d = 1500
c = np.array([0.5, 1.0, 2.0])
theta = np.array([1.6, 1.4, 1.5])
M = len(c)
n = (c * d).astype(int)
beta2 = np.array([beta_sq(t, ci) for t, ci in zip(theta, c)])
S = np.sum(beta2 / (1 - beta2))
print(f"beta^2 = {np.round(beta2, 4)}, S/(S+1) = {S / (S + 1):.4f}")

v = np.zeros(d)
v[0] = 1.0
Vt_rows = []
for i in range(M):
    X = spiked_table(np.array([theta[i]]), None, n[i], d, rng)
    vh = top_right_vectors(X, 1)[:, 0]
    if vh @ v < 0:
        vh = -vh
    Vt_rows.append(vh)
Vtil = np.array(Vt_rows)  # M x d
g = Vtil @ v
G = Vtil @ Vtil.T
bound = g @ np.linalg.solve(G, g)
print(f"g^T G^-1 g = {bound:.4f}  (limit S/(S+1) = {S / (S + 1):.4f})")

worst = -1.0
for trial in range(400):
    kind = trial % 4
    if kind == 0:
        w = rng.standard_normal(M)
    elif kind == 1:
        w = np.abs(rng.standard_normal(M))
    elif kind == 2:  # singular weight
        w = rng.standard_normal(M)
        w[rng.integers(M)] = 0.0
    else:  # data-dependent weight: proportional to the Rayleigh optimizer, perturbed
        w = np.linalg.solve(G, g) + 0.05 * rng.standard_normal(M)
    Wt = np.diag(w) @ Vtil
    vhat = top_right_vectors(Wt, 1)[:, 0]
    perf = (vhat @ v) ** 2
    worst = max(worst, perf - bound)
print(f"max over 400 w of perf(w) - g^T G^-1 g = {worst:.3e}  (must be <= ~1e-12)")
w_opt = 1 / np.sqrt(1 - beta2)
vhat = top_right_vectors(np.diag(w_opt) @ Vtil, 1)[:, 0]
print(f"perf(w_opt) = {(vhat @ v) ** 2:.4f}")
ok1 = worst <= 1e-9

# --------------------------------------------- C2, C3: rank r, audit's example
print("\n== C2/C3: rank r = 2, audit's example (rank B_R = 1 < r) ==")
d = 1500
r = 2
# table 1 sees component 1 with subcritical theta; tables 2 and 3 see component 2, supercritical
c = np.array([1.0, 1.0, 1.0])
thetas = [np.array([0.8, 0.0]), np.array([0.0, 1.8]), np.array([0.0, 1.6])]
M = 3
n = (c * d).astype(int)
V = np.eye(d)[:, :r]
blocks = []
for i in range(M):
    X = spiked_table(thetas[i], None, n[i], d, rng)
    blocks.append(top_right_vectors(X, 1))  # r_i = 1 each: d x 1
Vtil = np.hstack(blocks).T  # r_tilde x d, r_tilde = 3
g = Vtil @ V  # 3 x 2
G = Vtil @ Vtil.T
trace_bound = np.trace(g.T @ np.linalg.solve(G, g))

# paper's limit objects
beta = np.array([beta_sq(np.max(np.abs(t)), ci) for t, ci in zip(thetas, c)])  # per table
B_R = np.zeros((3, r))
B_R[0, 0] = np.sqrt(beta[0])  # 0: subcritical
B_R[1, 1] = np.sqrt(beta[1])
B_R[2, 1] = np.sqrt(beta[2])
A = B_R @ B_R.T + np.diag(1 - np.sum(B_R ** 2, axis=1))
L_opt = np.trace(B_R.T @ np.linalg.solve(A, B_R))
S2 = beta[1] / (1 - beta[1]) + beta[2] / (1 - beta[2])
print(f"beta^2 per table = {np.round(beta, 4)}  (table 1 subcritical)")
print(f"rank B_R = {np.linalg.matrix_rank(B_R)} < r = {r}")
print(f"L_opt = tr(B_R^T A^-1 B_R) = {L_opt:.4f};  S_2/(S_2+1) = {S2 / (S2 + 1):.4f}")
print(f"finite-N trace bound tr(g^T G^-1 g) = {trace_bound:.4f}")

worst = -1.0
for trial in range(300):
    if trial % 3 == 0:
        w = rng.standard_normal(3)
    elif trial % 3 == 1:
        w = np.abs(rng.standard_normal(3)) + 0.1
    else:
        w = rng.standard_normal(3)
        w[rng.integers(3)] = 0.0
    Vh = top_right_vectors(np.diag(w) @ Vtil, r)
    perf = np.sum((Vh.T @ V) ** 2)
    worst = max(worst, perf - trace_bound)
print(f"max over 300 W of perf(W) - trace bound = {worst:.3e}  (must be <= ~1e-12)")
ok2 = worst <= 1e-9

D = np.diag(1 - np.sum(B_R ** 2, axis=1))
W_opt = np.diag(1 / np.sqrt(np.diag(D)))
Vh = top_right_vectors(W_opt @ Vtil, r)
perf_opt = np.sum((Vh.T @ V) ** 2)
print(f"perf(W_opt) = {perf_opt:.4f}  vs L_opt = {L_opt:.4f}  (audit's example, no rank condition)")
ok3 = abs(perf_opt - L_opt) < 0.03

# also: with table 1 subcritical, dropping it (paper's removal remark) gives the same value
Vh_drop = top_right_vectors(W_opt[1:, 1:] @ Vtil[1:], r)
perf_drop = np.sum((Vh_drop.T @ V) ** 2)
print(f"perf after dropping the subcritical table = {perf_drop:.4f}")

# paper's closed form L_opt = r - sum of the r smallest eigenvalues of A^{-1/2} D A^{-1/2}
evals, evecs = np.linalg.eigh(A)
A_mhalf = evecs @ np.diag(evals ** -0.5) @ evecs.T
mu = np.sort(np.linalg.eigvalsh(A_mhalf @ D @ A_mhalf))
L_paper = r - np.sum(mu[:r])
print(f"paper's L_opt formula = {L_paper:.6f}  vs trace form = {L_opt:.6f}")
ok3 = ok3 and abs(L_paper - L_opt) < 1e-9

print("\nRESULT:", "PASS" if (ok1 and ok2 and ok3) else "FAIL")
sys.exit(0 if (ok1 and ok2 and ok3) else 1)
