#!/usr/bin/env python3
"""Numeric check of the three regimes of `swLimitEx` (item F18b, notes/paper_edits.md E12).

The instance of `prop:singleweight_suboptimality` (`main_paper.tex:2171` to `:2194`):
`R_1 = e_1`, `R_2 = e_2`, `Theta_1 = Theta_2 = theta_0`, `c_1 = c_2 = c_0`, on the Lean witness
`theta_0 = 8/5`, `c_0 = 1`, `n_i = d`. Single-weight stackSVD with weights `w = (1, ratio)`.

Sections:
  1. the closed forms `gammaEx`, `DetectableEx`, `swTermEx`, `swLimitEx` of
     `RankR/SingleWeight/Example.lean`, and two analytic checks (the tie equals the
     unweighted value at the doubled aspect ratio; the one-detectable term is below beta_0^2);
  2. the boundary ratio at which the smaller root stops being detectable;
  3. Monte Carlo at d = 600 and d = 2000 in each regime, with the per-index overlaps.

Run: python3 scripts/check_singleweight_regimes.py   (about 6 minutes with 2 threads)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import time
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np

SEED = 20260907
TH = 8 / 5
C = 1.0


# 1. closed forms (Lean names in RankR/SingleWeight/Example.lean and Defs.lean)
def betaSq(theta, c):
    return (theta ** 4 - c) / (theta ** 4 + theta ** 2) if theta ** 4 > c else 0.0


def gammaEx(w, l):
    return w[l] ** 2 * (1 + TH ** 2)


def thresh_sum(w, l):
    g = gammaEx(w, l)
    return sum(C * w[i] ** 4 / (g - w[i] ** 2) ** 2 for i in range(2))


def detectableEx(w, l):
    g = gammaEx(w, l)
    return (max(w[0] ** 2, w[1] ** 2) < g) and (thresh_sum(w, l) < 1)


def swTermEx(w, l):
    j = 1 - l
    g = gammaEx(w, l)
    return (TH ** 4 - C - C * w[j] ** 4 * TH ** 4 / (g - w[j] ** 2) ** 2) / (TH ** 2 * (1 + TH ** 2))


def swLimitEx(w):
    return sum(swTermEx(w, l) for l in range(2) if detectableEx(w, l))


def regime(w):
    d0, d1 = detectableEx(w, 0), detectableEx(w, 1)
    if abs(w[0] - w[1]) < 1e-12:
        return "tie"
    return {2: "both detectable", 1: "one detectable", 0: "none detectable"}[int(d0) + int(d1)]


# 3. Monte Carlo
def sim(w, d, reps, rng):
    perf, ov = [], []
    V = np.eye(d)[:, :2]                      # V = (e_1, e_2): R_1 = e_1, R_2 = e_2
    for _ in range(reps):
        n = d                                  # n_i = d, so c_i = 1
        Xs = []
        for i in range(2):
            u = rng.standard_normal(n)
            u /= np.linalg.norm(u)
            Z = rng.standard_normal((n, d)) / np.sqrt(d)
            Xs.append(w[i] * (TH * np.outer(u, V[:, i]) + Z))
        X = np.vstack(Xs)
        _, S, Vt = np.linalg.svd(X, full_matrices=False)
        assert S[0] >= S[1] >= S[2], "singular values are not sorted"
        Vh = Vt[:2].T
        perf.append(np.linalg.norm(Vh.T @ V, "fro") ** 2)
        ov.append([(Vt[0] @ V[:, 0]) ** 2, (Vt[0] @ V[:, 1]) ** 2,
                   (Vt[1] @ V[:, 0]) ** 2, (Vt[1] @ V[:, 1]) ** 2])
    return np.mean(perf), np.std(perf, ddof=1) / np.sqrt(reps), np.mean(ov, axis=0)


def main():
    t0 = time.time()
    print(f"seed {SEED}   theta_0 = {TH}   c_0 = {C}   beta_0^2 = {betaSq(TH, C):.6f}   "
          f"2 beta(theta_0, 2 c_0)^2 = {2 * betaSq(TH, 2 * C):.6f}")

    print("\n1. analytic checks")
    for a in (0.5, 1.0, 3.0):
        w = np.array([a, a])
        v = swLimitEx(w)
        print(f"   tie w = ({a}, {a}): swLimitEx {v:.9f}   2 beta(theta_0, 2c_0)^2 "
              f"{2 * betaSq(TH, 2 * C):.9f}   |diff| {abs(v - 2 * betaSq(TH, 2 * C)):.1e}")
        assert abs(v - 2 * betaSq(TH, 2 * C)) < 1e-12
    for r in (1.5, 2.0, 5.0):
        w = np.array([1.0, r])
        t = swTermEx(w, 1)
        print(f"   ratio {r}: term of the larger weight {t:.6f} < beta_0^2 {betaSq(TH, C):.6f}: "
              f"{t < betaSq(TH, C)}")
        assert t < betaSq(TH, C)
        assert detectableEx(w, 1), "the larger root must be detectable (2 c_0 < theta_0^4)"

    print("\n2. boundary ratio (threshold sum of the smaller root equals 1)")
    lo, hi = 1.0, 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if thresh_sum(np.array([1.0, mid]), 0) < 1:
            lo = mid
        else:
            hi = mid
    print(f"   ratio* = {lo:.6f}; smaller root detectable at 1.30 -> "
          f"{detectableEx(np.array([1.0, 1.30]), 0)}, at 1.31 -> "
          f"{detectableEx(np.array([1.0, 1.31]), 0)}")

    reps = {600: 12, 2000: 8}
    print("\n3. Monte Carlo (12 draws at d = 600, 8 at d = 2000), columns of the overlap "
          "block: [v0.idx0, v1.idx0, v0.idx1, v1.idx1]")
    rng = np.random.default_rng(SEED)
    print("   ratio  regime            d     limit     MC       +-se     overlaps")
    for r in (1.0, 1.2, 1.5, 2.0):
        w = np.array([1.0, r])
        for d in (600, 2000):
            m, s, ov = sim(w, d, reps[d], rng)
            print(f"   {r:4.1f}   {regime(w):16s} {d:5d}   {swLimitEx(w):.4f}   {m:.4f}  "
                  f"{s:.4f}   {np.round(ov, 4)}")
    print(f"\nelapsed {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
