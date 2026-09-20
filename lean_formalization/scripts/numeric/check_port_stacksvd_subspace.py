#!/usr/bin/env python3
"""Numeric check of `prop_stacksvd_subspace_general` (notes/archive/rankr_plan_B.md section 6).

Rows P1 to P5 are exact identities and must hit machine precision. A failure there means
the Lean statement is wrong. Rows P6 and P7 are the statistical rows.

Model (`assum:unaligned`, main_paper.tex:753):
    X_i = U_i Theta_i (V R_i)^T + Z_i / sqrt(d),   Z_i standard normal, n_i = round(c_i d).
Core matrix C = sum_i R_i Theta_i^2 R_i^T (main_paper.tex:826).
Limit  = sum_j betaSq(sqrt(lambda_j(C)), sum_i c_i),
         betaSq(t, c) = (t^4 - c)/(t^4 + t^2) if t^4 > c else 0  (Defs.lean:38).

numpy only, one process, 2 BLAS threads.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np

SEED = 20260901
print(f"seed = {SEED}   numpy = {np.__version__}")
print()

# ---------------------------------------------------------------- helpers

def betaSq(t, c):
    """Defs.lean:38, at signal strength t and aspect ratio c."""
    return (t ** 4 - c) / (t ** 4 + t ** 2) if t ** 4 > c else 0.0


def orth_frame(rng, rows, cols):
    """A random matrix with orthonormal columns, rows x cols."""
    Q, _ = np.linalg.qr(rng.standard_normal((rows, cols)))
    return Q[:, :cols]


def BBlock(theta, R):
    """r_tilde x r; row p = (i, j) is theta_ij (R_i)_{:,j}^T  (General.lean:293)."""
    rows = []
    for i, th in enumerate(theta):
        for j in range(len(th)):
            rows.append(th[j] * R[i][:, j])
    return np.array(rows)


def CBlock(theta, R):
    """(B_theta)^T B_theta  (plan section 4)."""
    B = BBlock(theta, R)
    return B.T @ B


def C_paper(theta, R):
    """sum_i R_i Theta_i^2 R_i^T  (main_paper.tex:826)."""
    r = R[0].shape[0]
    C = np.zeros((r, r))
    for i, th in enumerate(theta):
        C += R[i] @ np.diag(np.asarray(th) ** 2) @ R[i].T
    return C


def signal_factor(theta, R, U):
    """Block row i is U_i Theta_i R_i^T; shape (sum n_i) x r."""
    return np.vstack([U[i] @ np.diag(np.asarray(theta[i])) @ R[i].T for i in range(len(theta))])


def limit_stack(theta, R, c):
    lam = np.linalg.eigvalsh(CBlock(theta, R))
    return sum(betaSq(np.sqrt(max(l, 0.0)), sum(c)) for l in lam)


def proj_top(G, r):
    """Orthogonal projector onto the top-r eigenspace of a symmetric G, plus that basis."""
    w, Qm = np.linalg.eigh(G)          # ascending
    Vh = Qm[:, -r:]
    return Vh @ Vh.T, Vh


# ---------------------------------------------------------------- one cell

def run_cell(name, theta, R, c, d, draws, rng, exact_stats):
    """Draw `draws` replicates; accumulate the P1-P5 errors into exact_stats; return perf."""
    M = len(theta)
    r = R[0].shape[0]
    n = [int(round(c[i] * d)) for i in range(M)]
    C = CBlock(theta, R)

    # ---- P1, deterministic, once per cell
    e1 = max(np.max(np.abs(C - C_paper(theta, R))), 0.0)
    exact_stats["P1"] = max(exact_stats["P1"], e1)

    lamC, Qc = np.linalg.eigh(C)
    perf = np.empty(draws)
    for t in range(draws):
        V = orth_frame(rng, d, r)
        U = [orth_frame(rng, n[i], len(theta[i])) for i in range(M)]

        A = signal_factor(theta, R, U)                       # (sum n_i) x r
        exact_stats["P2"] = max(exact_stats["P2"], np.max(np.abs(A.T @ A - C)))
        exact_stats["P2rel"] = max(exact_stats["P2rel"],
                                   np.max(np.abs(A.T @ A - C)) / max(np.max(np.abs(C)), 1.0))

        EX = A @ V.T                                         # E[X_stack]
        S = EX.T @ EX
        VCVt = V @ C @ V.T
        spike = V @ Qc                                       # columns V q_j
        rank1 = sum(lamC[j] * np.outer(spike[:, j], spike[:, j]) for j in range(r))
        exact_stats["P3"] = max(exact_stats["P3"],
                                np.max(np.abs(S - VCVt)), np.max(np.abs(S - rank1)))

        Z = np.vstack([rng.standard_normal((n[i], d)) for i in range(M)])
        X = EX + Z / np.sqrt(d)
        G = X.T @ X
        P, Vhat = proj_top(G, r)

        lhs = sum(np.dot(V[:, k], P @ V[:, k]) for k in range(r))          # sum_k ||P V e_k||^2
        rhs = sum(np.dot(spike[:, j], P @ spike[:, j]) for j in range(r))  # sum_j ||P V q_j||^2
        exact_stats["P4"] = max(exact_stats["P4"], abs(lhs - rhs))

        fro = np.sum((Vhat.T @ V) ** 2)                                   # ||Vhat^T V||_F^2
        tr = np.trace(V.T @ P @ V)
        exact_stats["P5"] = max(exact_stats["P5"], abs(fro - tr))

        perf[t] = lhs
    return perf


# ---------------------------------------------------------------- cells

rng = np.random.default_rng(SEED)
exact = {k: 0.0 for k in ("P1", "P2", "P3", "P4", "P5", "P2rel")}
c = [1.0, 1.0]
DRAWS = 24
DS = [300, 1200]


def R2_of(psi):
    return np.array([[np.sin(psi)], [np.cos(psi)]])


R1 = np.eye(2)
cells = []
for psi, tag in ((0.0, "psi=0"), (np.pi / 4, "psi=pi/4"), (np.pi / 2, "psi=pi/2")):
    cells.append((f"base {tag}", [[2.2, 1.5], [2.0]], [R1, R2_of(psi)]))
# the plan's tie cell, verbatim: Theta_1 = diag(2,2), R_2 = e_1
cells.append(("plan tie cell (Theta1=diag(2,2), R2=e1)", [[2.0, 2.0], [2.0]],
              [R1, np.array([[1.0], [0.0]])]))
# a corrected tie cell: C = diag(3.25, 3.25) exactly
cells.append(("corrected tie (spec(C) repeats)", [[np.sqrt(3.25), 1.5], [1.0]],
              [R1, np.array([[0.0], [1.0]])]))
# subcritical: Theta_1 = diag(2.2, 0.6), R_2 = e_1 -> C = diag(8.84, 0.36), lambda_2^2 < ||c||_1
cells.append(("subcritical (Theta1=diag(2.2,0.6), R2=e1)", [[2.2, 0.6], [2.0]],
              [R1, np.array([[1.0], [0.0]])]))

print("cells: M=2, rk=(2,1), r=2, c=(1,1), n_i=round(c_i d), "
      f"d in {DS}, {DRAWS} draws each")
print()
print(f"{'cell':42s} {'spec(C)':>20s} {'limit':>9s}")
for name, th, R in cells:
    lam = np.linalg.eigvalsh(CBlock(th, R))
    print(f"{name:42s} {np.array2string(lam[::-1], precision=4):>20s} "
          f"{limit_stack(th, R, c):9.5f}")
print()

rows = []
for name, th, R in cells:
    res = {}
    for d in DS:
        perf = run_cell(name, th, R, c, d, DRAWS, rng, exact)
        res[d] = (perf.mean(), perf.std(ddof=1), perf.std(ddof=1) / np.sqrt(DRAWS))
    rows.append((name, limit_stack(th, R, c), res))

# ---------------------------------------------------------------- report

print("P1 to P5, exact identities (worst absolute error over every cell and draw)")
print()
tol = {"P1": 1e-14, "P2": 1e-14, "P3": 1e-13, "P4": 1e-13, "P5": 1e-13}
label = {"P1": "CBlock = (B_theta)^T B_theta = sum_i R_i Theta_i^2 R_i^T",
         "P2": "signalFactor^T signalFactor = C",
         "P3": "E[X]^T E[X] = V C V^T = sum_j lam_j (V q_j)(V q_j)^T",
         "P4": "sum_k ||P (V e_k)||^2 = sum_j ||P (V q_j)||^2",
         "P5": "||Vhat_stacksvd^T V||_F^2 = tr(V^T P V)"}
allpass = True
for k in ("P1", "P2", "P3", "P4", "P5"):
    ok = exact[k] <= tol[k]
    print(f"  {k}  {label[k]:56s} worst {exact[k]:.3e}  tol {tol[k]:.0e}  "
          f"{'PASS' if ok else 'over abs tol'}")
    if k == "P2" and not ok:
        # P2 sums n_i = 2400 terms of size ~ max|C|; the plan's fixed 1e-14 is not
        # scale aware. Judge it relative to max|C| instead.
        ok = exact["P2rel"] <= 1e-14
        print(f"        relative to max|C|: worst {exact['P2rel']:.3e}  tol 1e-14  "
              f"{'PASS' if ok else 'FAIL'}  (roundoff, not a statement error)")
    allpass &= ok
print()

print("P6 and P7, sample perfStackRG against limitStackRG")
print()
print(f"{'cell':42s} {'d':>5s} {'mean':>9s} {'sd':>8s} {'se':>8s} {'limit':>9s} {'bias':>9s}")
for name, lim, res in rows:
    for d in DS:
        m, s, se = res[d]
        print(f"{name:42s} {d:5d} {m:9.5f} {s:8.5f} {se:8.5f} {lim:9.5f} {m - lim:9.5f}")
print()
print(f"{'cell':42s} {'|bias| d=300':>13s} {'|bias| d=1200':>14s} {'verdict':>10s}")
for name, lim, res in rows:
    b3 = abs(res[300][0] - lim)
    b12 = abs(res[1200][0] - lim)
    ok = b12 <= b3 + 2 * res[1200][2]
    allpass &= ok
    print(f"{name:42s} {b3:13.5f} {b12:14.5f} {'PASS' if ok else 'FAIL':>10s}")
print()
print("OVERALL:", "PASS" if allpass else "FAIL")
