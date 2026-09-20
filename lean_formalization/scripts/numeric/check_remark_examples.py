"""Closed-form checks of the four examples in remark:stack_outperform_svd and
remark:svd_outperform_stack (main_paper.tex lines 566 to 620). Deterministic; no seed.

Formulas (paper): beta_i^2 = (theta^4 - c)/(theta^4 + theta^2) if theta^4 > c else 0;
stackSVD limit on a set S: (T^2 - C)/(T^2 + T) with T = sum theta_i^2, C = sum c_i, if T^2 > C
else 0; unweighted svdstack limit: (beta . x)^2 / lam_max with x the top unit eigenvector of
A_beta = beta beta^T + diag(1 - beta^2); weighted svdstack optimum: S/(S+1),
S = sum beta_i^2/(1 - beta_i^2).

Finding: in example (iii) the paper takes theta_3 = c_3^{1/4}, so theta_3^4 = c_3 and cor.2's
binary set {i : theta_i^4 > c_i} excludes table 3. Binary stackSVD is then the stack of
tables 1 and 2 with limit 62/72 = 0.8611, above the svdstack value 6/7 = 0.8571. The value
the paper computes is the unweighted stack of all three tables. The fix theta_3 = (c_3+1)^{1/4}
puts table 3 inside the binary set and keeps the stated conclusion (notes/paper_edits.md E6).

Run:  python3 scripts/check_remark_examples.py
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np


def beta_sq(theta, c, theta4=None):
    """theta4 overrides theta**4 so that an example at the exact threshold is not decided by
    floating-point rounding (c3**0.25 to the fourth power is c3 + 1e-14)."""
    t4 = theta ** 4 if theta4 is None else theta4
    return (t4 - c) / (t4 + theta ** 2) if t4 > c else 0.0


def stack_limit(theta, c):
    T, C = np.sum(np.square(theta)), np.sum(c)
    return (T ** 2 - C) / (T ** 2 + T) if T ** 2 > C else 0.0


def binary_set(theta, c, theta4):
    return [i for i in range(len(theta)) if theta4[i] > c[i]]


def svdstack_limit(beta):
    A = np.outer(beta, beta) + np.diag(1 - beta ** 2)
    w, V = np.linalg.eigh(A)
    return (beta @ V[:, -1]) ** 2 / w[-1]


def svdstack_opt(beta):
    S = np.sum(beta ** 2 / (1 - beta ** 2))
    return S / (S + 1)


ok = True
print("== (i) theta_i = c_i = 1 ==")
for M in (2, 3, 10):
    th, c = np.ones(M), np.ones(M)
    b = np.array([beta_sq(t, ci) for t, ci in zip(th, c)])
    s = stack_limit(th, c)
    print(f"M={M}: stackSVD {s:.4f} vs 1-2/(M+1) = {1 - 2 / (M + 1):.4f}; beta = {b}; svdstack opt {svdstack_opt(b):.4f}")
    ok &= abs(s - (1 - 2 / (M + 1))) < 1e-12 and np.all(b == 0)

print("\n== (ii) c_i = c0, theta_1 = theta_2 = theta0 > c0^{1/4}, rest 0 ==")
th0, c0, M = 1.5, 1.0, 25
th = np.zeros(M); th[:2] = th0; c = c0 * np.ones(M)
b = np.array([beta_sq(t, ci) for t, ci in zip(th, c)])
b0sq = b[0]
print(f"M={M} > 4 theta0^4/c0 = {4 * th0 ** 4 / c0:.2f}: stackSVD {stack_limit(th, c):.4f}; "
      f"svdstack {svdstack_limit(np.sqrt(b)):.4f} vs 2b0^2/(1+b0^2) = {2 * b0sq / (1 + b0sq):.4f}")
ok &= stack_limit(th, c) == 0 and abs(svdstack_limit(np.sqrt(b)) - 2 * b0sq / (1 + b0sq)) < 1e-9

print("\n== (iii) c = (1, 1, c3), theta = (2, 2, c3^{1/4}) ==")
for c3 in (100.0, 1e4, 1e6):
    th = np.array([2.0, 2.0, c3 ** 0.25]); c = np.array([1.0, 1.0, c3])
    th4 = np.array([16.0, 16.0, c3])  # exact fourth powers
    b = np.array([beta_sq(t, ci, t4) for t, ci, t4 in zip(th, c, th4)])
    bset = binary_set(th, c, th4)
    s_full = stack_limit(th, c)
    s_bin = stack_limit(th[bset], c[bset])
    paper = (16 * np.sqrt(c3) + 62) / (c3 + 17 * np.sqrt(c3) + 72)
    sv = svdstack_limit(np.sqrt(b))
    print(f"c3={c3:g}: beta^2 = {np.round(b, 4)}, binary set {bset}, unweighted stack {s_full:.4f} "
          f"(paper formula {paper:.4f}), binary stack {s_bin:.4f}, svdstack {sv:.4f} (6/7 = {6 / 7:.4f})")
    ok &= abs(s_full - paper) < 1e-9 and abs(sv - 6 / 7) < 1e-9 and bset == [0, 1]
    ok &= s_bin > sv  # binary stackSVD beats svdstack in the paper's own example
print("  fix theta_3 = (c3+1)^{1/4}:")
for c3 in (100.0, 1e4, 1e6):
    th = np.array([2.0, 2.0, (c3 + 1) ** 0.25]); c = np.array([1.0, 1.0, c3])
    th4 = np.array([16.0, 16.0, c3 + 1])
    b = np.array([beta_sq(t, ci, t4) for t, ci, t4 in zip(th, c, th4)])
    bset = binary_set(th, c, th4)
    s_bin = stack_limit(th[bset], c[bset])
    sv = svdstack_limit(np.sqrt(b))
    print(f"  c3={c3:g}: beta_3^2 = {b[2]:.2e}, binary set {bset}, binary stack {s_bin:.4f}, svdstack {sv:.4f}")
    ok &= bset == [0, 1, 2] and sv > s_bin

print("\n== (iv) theta = (sqrt5, 4), c = (1, 38.4) ==")
th = np.array([np.sqrt(5), 4.0]); c = np.array([1.0, 38.4])
b = np.array([beta_sq(t, ci) for t, ci in zip(th, c)])
subsets = {"both": [0, 1], "table 1": [0], "table 2": [1]}
vals = {k: stack_limit(th[s], c[s]) for k, s in subsets.items()}
print(f"beta^2 = {np.round(b, 4)}; binary stack {vals}; weighted svdstack {svdstack_opt(np.sqrt(b)):.4f}; "
      f"unweighted svdstack {svdstack_limit(np.sqrt(b)):.4f}")
ok &= abs(vals['both'] - 0.8693) < 1e-3 and abs(vals['table 1'] - 0.8) < 1e-9 and abs(svdstack_opt(np.sqrt(b)) - 8 / 9) < 1e-9
ok &= svdstack_opt(np.sqrt(b)) > max(vals.values())
print("\nRESULT:", "PASS" if ok else "FAIL")
