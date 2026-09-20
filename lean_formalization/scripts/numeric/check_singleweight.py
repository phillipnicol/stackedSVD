#!/usr/bin/env python3
"""Numeric audit of the two single-weight appendix results of "Stacked SVD or SVD stacked?".

Targets (main_paper.tex, read-only snapshot in the repo root):

  * assum:gen_rank_stacksvd_eig_sep   (line 2104)
  * prop:gen_rank_stacksvd_singleweight (line 2112) and its proof (2124 to 2167)
  * prop:singleweight_suboptimality   (line 915), proof at 2169 to 2194

Model (assum:unaligned, line 753), weighted per table:

    X_i = U_i Theta_i (V R_i)^T + E_i,  E_i entries iid N(0, 1/d),  c_i = n_i / d
    X_stack(w) = [w_1 X_1; ...; w_M X_M]
    Vhat_stacksvd(w) = top r eigenvectors of sum_i w_i^2 X_i^T X_i
    performance = || Vhat(w)^T V ||_F^2

Deterministic layer of the proposition:

    A     = [w_1 U_1 Theta_1 R_1^T; ...; w_M U_M Theta_M R_M^T]   (n x r)
    Sigma = diag(w_i^2 I_{n_i}),   G = A A^T + Sigma
    Mmat(g)  = sum_i w_i^2/(g - w_i^2)  R_i Theta_i^2 R_i^T       (r x r)
    Kmat(g)  = sum_i w_i^2/(g - w_i^2)^2 R_i Theta_i^2 R_i^T      (r x r)
    gamma_l  : roots of det(I_r - Mmat(g)) = 0 above max_i w_i^2
    z_l      : unit eigenvector of Mmat(gamma_l) with eigenvalue 1
    eta_l    = 1 - sum_i c_i w_i^4/(gamma_l - w_i^2)^2
    perf     = sum_l eta_l / (gamma_l * z_l^T Kmat(gamma_l) z_l)

Runs on 2 cores in under 5 minutes. Prints its seed.
"""

import os

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["NUMEXPR_NUM_THREADS"] = "2"
os.environ["VECLIB_MAXIMUM_THREADS"] = "2"

import time  # noqa: E402

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np  # noqa: E402

SEED = 20260905
BIG = 1.0e12


# ----------------------------------------------------------------------------
# 1. Deterministic scalar layer
# ----------------------------------------------------------------------------


def signal_blocks(R, Theta):
    """R_i Theta_i^2 R_i^T for each table."""
    return [Ri @ np.diag(np.asarray(Ti, float) ** 2) @ Ri.T for Ri, Ti in zip(R, Theta)]


def Mmat(g, w, S):
    r = S[0].shape[0]
    out = np.zeros((r, r))
    for wi, Si in zip(w, S):
        if wi != 0.0:
            out += (wi**2 / (g - wi**2)) * Si
    return out


def Kmat(g, w, S):
    r = S[0].shape[0]
    out = np.zeros((r, r))
    for wi, Si in zip(w, S):
        if wi != 0.0:
            out += (wi**2 / (g - wi**2) ** 2) * Si
    return out


def factor_blocks(R, Theta):
    """`C_i = R_i Theta_i`, so that `S_i = C_i C_i^T`."""
    return [Ri @ np.diag(np.asarray(Ti, float)) for Ri, Ti in zip(R, Theta)]


def secular_roots(w, R, Theta, r):
    """All roots of `det(I_r - Mmat(g)) = 0` with `g > max_i w_i^2`, largest first.

    Exact linear algebra, no bisection.  With `C_i = R_i Theta_i` (`r x rk_i`), stack
    `B = [w_1 C_1^T; ...; w_M C_M^T]` (`rtot x r`) and set `W = diag(w_i^2 I_{rk_i})`.  Then
    `B^T (g I - W)^{-1} B = sum_i w_i^2/(g - w_i^2) C_i C_i^T = Mmat(g)`, so the matrix
    determinant lemma gives

        det(W + B B^T - g I) = det(W - g I) det(I_r - Mmat(g)).

    The roots above `max_i w_i^2` are therefore exactly the eigenvalues of the symmetric
    `rtot x rtot` matrix `H = W + B B^T` above that threshold, with multiplicity.  This is the
    paper's own reduction (`main_paper.tex:2124`) compressed from `n x n` to `rtot x rtot`.

    An earlier bisection implementation was wrong twice: it invented a root from `eigvalsh`
    roundoff on `Mmat` evaluated ~1e-13 from the pole, and it divided by zero once
    `wmax2 + step == wmax2` in floating point.  Both are gone with the pole.
    """
    C = factor_blocks(R, Theta)
    rows, diag = [], []
    for wi, Ci in zip(w, C):
        rows.append(wi * Ci.T)                       # rk_i x r
        diag.append(np.full(Ci.shape[1], wi ** 2))   # rk_i entries
    B = np.vstack(rows)
    H = np.diag(np.concatenate(diag)) + B @ B.T
    ev = np.linalg.eigvalsh(H)[::-1]
    wmax2 = max(wi ** 2 for wi in w)
    thr = wmax2 * (1.0 + 1e-12) + 1e-14
    roots = [float(x) for x in ev if x > thr]
    return roots[:r]


def unit_eigvecs_one(Mg, mult=1, Kg=None):
    """`mult` eigenvectors of `Mg` for the eigenvalues closest to 1, K-orthogonal at a tie.

    `mult > 1` is the tied-root case: the eigenvalue 1 of `Mmat(gamma)` then has multiplicity
    `mult` and no single `z` is defined.  Euclidean orthonormality is NOT enough there.  The
    proof sets `xi_l = (gamma - Sigma)^-1 A z_l / ||.||`, and

        <xi_l, xi_m> proportional to z_l^T A^T (gamma - Sigma)^-2 A z_m = z_l^T Kmat(gamma) z_m,

    so the `xi_l` form an orthonormal eigenbasis of `G` only when the `z_l` are **K-orthogonal**.
    The performance sums `1/(z_l^T K z_l)`, a sum of reciprocals, which is not invariant under a
    rotation inside the tied plane.  So the basis must be the eigenbasis of `Kg` restricted to
    the eigenspace: take an orthonormal basis `Q` of the eigenspace, diagonalise `Q^T Kg Q`, and
    rotate.  That basis is Euclidean-orthonormal and K-orthogonal at once.

    Measured on a tie (`M=2, r=2, R_1=e_1, R_2=e_2, theta^2=(3,15), w=(1,.5), c=(.05,.05)`,
    where `Mmat(4) = I` and `Kmat(4) = diag(1/3, 4/15)`): the K-eigenbasis gives 1.677750 and a
    45 degree rotation of it gives 1.657037, while the Monte Carlo gives 1.680 +/- 0.002 at
    `d = 1600`.  The rotated basis is wrong by 13 standard errors.
    """
    vals, vecs = np.linalg.eigh(Mg)
    order = np.argsort(np.abs(vals - 1.0))[:mult]
    Q = np.column_stack([vecs[:, k] / np.linalg.norm(vecs[:, k]) for k in order])
    if mult > 1 and Kg is not None:
        _, W = np.linalg.eigh(Q.T @ Kg @ Q)
        Q = Q @ W
    Z = [Q[:, j] for j in range(Q.shape[1])]
    return Z, max(abs(vals[k] - 1.0) for k in order)


def group_roots(gam, rtol=1e-9):
    """Group equal roots (a tie is a repeated eigenvalue of `H`), largest first."""
    groups = []
    for g in gam:
        if groups and abs(g - groups[-1][0]) <= rtol * max(1.0, abs(g)):
            groups[-1][1] += 1
        else:
            groups.append([g, 1])
    return [(g, m) for g, m in groups]


def eta_of(g, w, c):
    return 1.0 - sum(ci * wi**4 / (g - wi**2) ** 2 for wi, ci in zip(w, c))


def perf_formula(w, c, R, Theta):
    """The right-hand side of prop:gen_rank_stacksvd_singleweight, plus diagnostics."""
    S = signal_blocks(R, Theta)
    r = S[0].shape[0]
    gam = secular_roots(w, R, Theta, r)
    terms, zres = [], []
    for g, mult in group_roots(gam):
        Kg = Kmat(g, w, S)
        Z, res = unit_eigvecs_one(Mmat(g, w, S), mult, Kg)
        zres.append(res)
        for z in Z:
            den = g * float(z @ Kg @ z)
            terms.append(eta_of(g, w, c) / den)
    a4 = [1.0 - eta_of(g, w, c) for g in gam]  # sum_i c_i w_i^4/(gamma-w_i^2)^2
    return {
        "gamma": gam,
        "terms": terms,
        "perf": float(sum(terms)),
        "assum4": a4,
        "n_roots": len(gam),
        "z_residual": max(zres) if zres else float("nan"),
    }


def build_G(w, R, Theta, nsmall):
    """G = A A^T + Sigma on small blocks, plus A and Sigma, for the exact identity checks."""
    rng = np.random.default_rng(SEED + 7)
    r = R[0].shape[0]
    blocks, sig = [], []
    for i, (Ri, Ti, ni) in enumerate(zip(R, Theta, nsmall)):
        ri = Ri.shape[1]
        Q, _ = np.linalg.qr(rng.standard_normal((ni, ri)))
        blocks.append(w[i] * Q @ np.diag(np.asarray(Ti, float)) @ Ri.T)
        sig.append(np.full(ni, w[i] ** 2))
    A = np.vstack(blocks)
    sigma = np.concatenate(sig)
    G = A @ A.T + np.diag(sigma)
    return G, A, sigma


def roots_via_G(w, R, Theta, r, nsmall, rng):
    """The roots found a second way: eigenvalues of an explicitly built `G` above `max w_i^2`.

    `secular_roots` compresses the problem to `H = W + B B^T` (`rtot x rtot`); this builds the
    full `G = A A^T + Sigma` (`n x n`) and calls LAPACK on it.  The two agree only if the
    matrix determinant lemma of the paper's proof holds and `secular_roots` misses no root and
    invents none.
    """
    blocks, sig = [], []
    for i, (Ri, Ti, ni) in enumerate(zip(R, Theta, nsmall)):
        ri = Ri.shape[1]
        Q, _ = np.linalg.qr(rng.standard_normal((ni, ri)))
        blocks.append(w[i] * Q @ np.diag(np.asarray(Ti, float)) @ Ri.T)
        sig.append(np.full(ni, w[i] ** 2))
    A = np.vstack(blocks)
    G = A @ A.T + np.diag(np.concatenate(sig))
    ev = np.linalg.eigvalsh(G)[::-1]
    wmax2 = max(wi ** 2 for wi in w)
    return [x for x in ev if x > wmax2 * (1 + 1e-9)]


def root_cross_check(rng, trials=300):
    """`secular_roots` against `roots_via_G` on random instances. Returns (mismatches, worst)."""
    bad, worst = 0, 0.0
    for _ in range(trials):
        M = int(rng.integers(1, 5))
        r = int(rng.integers(1, 5))
        rk = [int(rng.integers(1, r + 1)) for _ in range(M)]
        R = [np.linalg.qr(rng.standard_normal((r, k)))[0][:, :k] for k in rk]
        TH = [np.abs(rng.standard_normal(k)) * 2.0 + 0.2 for k in rk]
        w = list(np.abs(rng.standard_normal(M)) * 1.0 + 0.05)
        nsmall = [int(rng.integers(k + 2, k + 14)) for k in rk]
        got = secular_roots(w, R, TH, r)
        ref = roots_via_G(w, R, TH, r, nsmall, rng)
        if len(got) != len(ref):
            bad += 1
            continue
        for a, b in zip(got, ref):
            worst = max(worst, abs(a - b) / max(1.0, abs(b)))
    return bad, worst


# ----------------------------------------------------------------------------
# 2. Rank-1 reference closed forms (Lean: StackSVDWeighted.lean, Scalars.lean)
# ----------------------------------------------------------------------------


def betaSq(theta, c):
    """Defs.lean betaSq: (th^4-c)/(th^4+th^2) above threshold, 0 below."""
    return (theta**4 - c) / (theta**4 + theta**2) if theta**4 > c else 0.0


def gammaTop_rank1(theta, w):
    """Root of sum_j th_j^2 w_j^2/(g - w_j^2) = 1 above max w_j^2 (Secular.lean)."""
    wmax2 = max(wi**2 for wi in w)

    def f(g):
        return sum(t**2 * wi**2 / (g - wi**2) for t, wi in zip(theta, w)) - 1.0

    lo = wmax2 * (1 + 1e-14) + 1e-14
    step = max(1.0, wmax2) * 1e-3
    while f(lo) <= 0 and step > 1e-16:
        lo = wmax2 + step
        step *= 0.5
    if f(lo) <= 0:
        return None
    hi = max(2 * lo, lo + 1.0)
    while f(hi) > 0 and hi < BIG:
        hi *= 2.0
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if f(m) > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def Lw_rank1(theta, c, w):
    """StackSVDWeighted.lean Lw: eta1 / LwDen under Assumption4, else 0."""
    g = gammaTop_rank1(theta, w)
    if g is None:
        return 0.0
    a4 = sum(ci * wi**4 / (g - wi**2) ** 2 for wi, ci in zip(w, c))
    if a4 >= 1.0:
        return 0.0
    den = g * sum(t**2 * wi**2 / (g - wi**2) ** 2 for t, wi in zip(theta, w))
    return (1.0 - a4) / den


def stackSVDLimitW(theta, c):
    """Scalars.lean stackSVDLimitW: root x in (0,1) of sum_i th_i^4(1-x)/(c_i + x th_i^2) = 1."""

    def gWf(x):
        return sum(t**4 * (1 - x) / (ci + x * t**2) for t, ci in zip(theta, c)) - 1.0

    if gWf(0.0) <= 0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if gWf(m) > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def optWstack(theta, c):
    return [t / np.sqrt(t**2 + ci) for t, ci in zip(theta, c)]


# ----------------------------------------------------------------------------
# 3. Monte Carlo
# ----------------------------------------------------------------------------


def draw_tables(d, c, R, Theta, rng):
    """One draw of the M tables with V = first r columns of I_d (rotation invariance)."""
    r = R[0].shape[0]
    Xs = []
    for Ri, Ti, ci in zip(R, Theta, c):
        ri = Ri.shape[1]
        ni = int(round(ci * d))
        U, _ = np.linalg.qr(rng.standard_normal((ni, ri)))
        VR = np.zeros((d, ri))
        VR[:r, :] = Ri  # V R_i with V = I_d[:, :r]
        X = U @ np.diag(np.asarray(Ti, float)) @ VR.T
        X += rng.standard_normal((ni, d)) / np.sqrt(d)
        Xs.append(X)
    return Xs


def perf_stacksvd(Xs, w, r):
    d = Xs[0].shape[1]
    S = np.zeros((d, d))
    for X, wi in zip(Xs, w):
        if wi != 0.0:
            S += (wi**2) * (X.T @ X)
    vals, vecs = np.linalg.eigh(S)
    Vhat = vecs[:, -r:]
    return float(np.sum(Vhat[:r, :] ** 2))  # ||Vhat^T V||_F^2, V = I_d[:, :r]


def perf_svdstack_unweighted(Xs, R, r):
    """Top r right singular vectors of Vtilde, the stack of the per-table top r_i ones."""
    d = Xs[0].shape[1]
    rows = []
    for X, Ri in zip(Xs, R):
        ri = Ri.shape[1]
        vals, vecs = np.linalg.eigh(X.T @ X)
        rows.append(vecs[:, -ri:][:, ::-1].T)  # r_i x d, strongest first
    Vt = np.vstack(rows)
    _, _, Vh = np.linalg.svd(Vt, full_matrices=False)
    Vhat = Vh[:r, :].T
    return float(np.sum(Vhat[:r, :] ** 2))


def mc(d, c, w, R, Theta, reps, rng, also_svdstack=False):
    r = R[0].shape[0]
    a, b = [], []
    for _ in range(reps):
        Xs = draw_tables(d, c, R, Theta, rng)
        a.append(perf_stacksvd(Xs, w, r))
        if also_svdstack:
            b.append(perf_svdstack_unweighted(Xs, R, r))
    a = np.asarray(a)
    out = {"mean": a.mean(), "se": a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0}
    if also_svdstack:
        b = np.asarray(b)
        out["sv_mean"] = b.mean()
        out["sv_se"] = b.std(ddof=1) / np.sqrt(len(b)) if len(b) > 1 else 0.0
    return out


# ----------------------------------------------------------------------------
# 4. Instances
# ----------------------------------------------------------------------------


def rot(psi):
    return np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]])


# Instance A: general rotations, r = 2, r_1 = r_2 = 2, non-commuting signal blocks.
A_R = [np.eye(2), rot(0.6)]
A_TH = [np.array([1.9, 1.35]), np.array([1.7, 1.15])]
A_C = [1.0, 1.0]
A_W = [1.0, 0.75]

# Instance A2: three tables, r = 3, r_i = 2, orthonormal columns from a fixed QR.
_rngA2 = np.random.default_rng(SEED + 11)
A2_R = [np.linalg.qr(_rngA2.standard_normal((3, 2)))[0] for _ in range(3)]
A2_TH = [np.array([2.6, 2.0]), np.array([2.4, 1.9]), np.array([2.7, 1.8])]
A2_C = [0.5, 0.75, 1.0]
A2_W = [1.0, 0.85, 0.7]

# Instance B: the suboptimality instance of the paper (R_1 = e_1, R_2 = e_2).
B_TH0 = 1.6
B_C0 = 1.0
B_R = [np.array([[1.0], [0.0]]), np.array([[0.0], [1.0]])]
B_TH = [np.array([B_TH0]), np.array([B_TH0])]
B_C = [B_C0, B_C0]


def sep(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def main():
    t0 = time.time()
    print(f"check_singleweight.py   seed = {SEED}   numpy {np.__version__}")
    print(f"OMP_NUM_THREADS = {os.environ['OMP_NUM_THREADS']}")

    # ---------------------------------------------------------------- 0. sanity
    sep("0. SANITY: rank-1 reduction (r = 1, r_i = 1, R_i = [1])")
    print("Compared with: StackSVDWeighted.lean Lw (gammaTop / eta1 / LwDen),")
    print("               Defs.lean betaSq, Scalars.lean stackSVDLimitW + optWstack.")
    theta = [1.4, 0.9, 1.7]
    cc = [1.0, 0.6, 1.5]
    R1 = [np.array([[1.0]]) for _ in theta]
    TH1 = [np.array([t]) for t in theta]
    worst = 0.0
    for w in ([1.0, 1.0, 1.0], [1.0, 0.7, 1.3], optWstack(theta, cc)):
        gen = perf_formula(w, cc, R1, TH1)
        ref = Lw_rank1(theta, cc, w)
        g1 = gammaTop_rank1(theta, w)
        dg = abs(gen["gamma"][0] - g1)
        dp = abs(gen["perf"] - ref)
        worst = max(worst, dg, dp)
        print(f"  w = {np.round(w, 6).tolist()}")
        print(
            f"    gamma  general {gen['gamma'][0]:.12f}  rank-1 {g1:.12f}  |diff| {dg:.3e}"
        )
        print(f"    perf   general {gen['perf']:.12f}  Lw     {ref:.12f}  |diff| {dp:.3e}")
    # single table, w = 1: must be betaSq
    g1 = perf_formula([1.0], [1.0], [np.array([[1.0]])], [np.array([1.4])])
    print(f"  M = 1, theta = 1.4, c = 1: general {g1['perf']:.12f}  betaSq "
          f"{betaSq(1.4, 1.0):.12f}  |diff| {abs(g1['perf'] - betaSq(1.4, 1.0)):.3e}")
    worst = max(worst, abs(g1["perf"] - betaSq(1.4, 1.0)))
    # optimal weights: Lw(optWstack) must equal stackSVDLimitW
    wo = optWstack(theta, cc)
    print(f"  Lw(optWstack) {Lw_rank1(theta, cc, wo):.12f}  stackSVDLimitW "
          f"{stackSVDLimitW(theta, cc):.12f}  |diff| "
          f"{abs(Lw_rank1(theta, cc, wo) - stackSVDLimitW(theta, cc)):.3e}")
    worst = max(worst, abs(Lw_rank1(theta, cc, wo) - stackSVDLimitW(theta, cc)))
    print(f"  SANITY worst deviation: {worst:.3e}")

    # ------------------------------------------------- 1. exact secular identity
    sep("1. EXACT: det(G - g I) = det(Sigma - g I) det(I_r - Mmat(g))  (proof, line 2124)")
    print("Top eigenvalues of G above max_i w_i^2 vs the roots of det(I_r - Mmat) = 0;")
    print("eigenvector identity A^T xi_l = z_l / ||(g - Sigma)^-1 A z_l|| (line 2143).")
    for name, R, TH, w, c, nsm in (
        ("A  (M=2, r=2, r_i=2)", A_R, A_TH, A_W, A_C, [17, 13]),
        ("A2 (M=3, r=3, r_i=2)", A2_R, A2_TH, A2_W, A2_C, [15, 14, 12]),
        ("B  (M=2, r=2, r_i=1)", B_R, B_TH, [1.0, 0.85], B_C, [11, 9]),
    ):
        S = signal_blocks(R, TH)
        r = S[0].shape[0]
        G, Amat, sigma = build_G(w, R, TH, nsm)
        ev = np.linalg.eigvalsh(G)[::-1]
        wmax2 = max(wi**2 for wi in w)
        ev_out = [x for x in ev if x > wmax2 * (1 + 1e-9)]
        roots = secular_roots(w, R, TH, r)
        k = min(len(ev_out), len(roots))
        dev = max([abs(ev_out[i] - roots[i]) for i in range(k)], default=float("nan"))
        # eigenvector identity
        derr = 0.0
        for g, mult in group_roots(roots):
            for z in unit_eigvecs_one(Mmat(g, w, S), mult, Kmat(g, w, S))[0]:
                y = Amat @ z / (g - sigma)
                xi = y / np.linalg.norm(y)
                lhs = Amat.T @ xi
                rhs = z / np.linalg.norm(y)
                derr = max(derr, float(np.linalg.norm(lhs - rhs)))
                derr = max(derr, float(np.linalg.norm(G @ xi - g * xi)))
                derr = max(derr,
                           abs(np.dot(lhs, lhs) - 1.0 / float(z @ Kmat(g, w, S) @ z)))
        print(f"  {name}:  eigs(G) above wmax^2: {len(ev_out)}   secular roots: {len(roots)}")
        print(f"      roots {['%.10f' % x for x in roots]}")
        print(f"      max |eig(G) - root| = {dev:.3e}   eigenvector identities max err "
              f"= {derr:.3e}")

    sep("1b. CROSS-CHECK of secular_roots on random instances")
    print("secular_roots uses H = W + B B^T (rtot x rtot); roots_via_G builds the full")
    print("G = A A^T + Sigma (n x n) and calls LAPACK. Agreement checks that no root is")
    print("missed and none is invented, at every rank and every alignment.")
    rngx = np.random.default_rng(SEED + 3131)
    bad, worst = root_cross_check(rngx, trials=300)
    print(f"  300 random instances (M, r in 1..4, r_i in 1..r): count mismatches = {bad}, "
          f"worst relative root error = {worst:.3e}")
    print("  targeted edge cases:")
    edge = {
        "rank-deficient sum R R^T (r=2, both e1)":
            ([1.0, 0.8], [np.array([[1.], [0.]]), np.array([[1.], [0.]])],
             [np.array([1.6]), np.array([1.6])], 2, [11, 9]),
        "one weight zero":
            ([1.0, 0.0], [np.eye(2), np.array([[0.], [1.]])],
             [np.array([1.9, 1.35]), np.array([1.7])], 2, [12, 9]),
        "very disparate weights (1 vs 1e-3)":
            ([1.0, 1e-3], [np.eye(2), np.eye(2)],
             [np.array([1.9, 1.35]), np.array([1.7, 1.15])], 2, [12, 11]),
        "root barely above wmax^2":
            ([1.0, 1.0], [np.array([[1.], [0.]]), np.array([[0.], [1.]])],
             [np.array([3.0]), np.array([0.02])], 2, [11, 9]),
        "exactly tied roots (paper instance, w = (1,1))":
            ([1.0, 1.0], [np.array([[1.], [0.]]), np.array([[0.], [1.]])],
             [np.array([1.6]), np.array([1.6])], 2, [11, 9]),
    }
    rnge = np.random.default_rng(SEED + 4242)
    for name, (w, R, TH, r, nsm) in edge.items():
        got = secular_roots(w, R, TH, r)
        ref = roots_via_G(w, R, TH, r, nsm, rnge)
        ok = len(got) == len(ref) and all(
            abs(a - b) <= 1e-9 * max(1, abs(b)) for a, b in zip(got, ref))
        print(f"    [{'OK  ' if ok else 'DIFF'}] {name}: {len(got)} root(s) "
              f"{['%.8f' % x for x in got]}")

    # ---------------------------------------------- 2. Monte Carlo, proposition
    sep("2. MONTE CARLO: prop:gen_rank_stacksvd_singleweight  (line 2112)")
    reps = {400: 30, 800: 16, 1600: 10}
    rng = np.random.default_rng(SEED)
    rows = []
    for name, R, TH, w, c in (
        ("A  w=(1,.75)", A_R, A_TH, A_W, A_C),
        ("A2 w=(1,.85,.7)", A2_R, A2_TH, A2_W, A2_C),
        ("B  w=(1,.85)", B_R, B_TH, [1.0, 0.85], B_C),
        ("B  w=(1,1)", B_R, B_TH, [1.0, 1.0], B_C),
    ):
        f = perf_formula(w, c, R, TH)
        print(f"\n  {name}")
        print(f"    Theta            {[np.round(t, 4).tolist() for t in TH]}")
        print(f"    c                {c}    w  {np.round(w, 4).tolist()}")
        print(f"    gamma            {['%.6f' % x for x in f['gamma']]}")
        print(f"    per-root sum_i c_i w_i^4/(g-w_i^2)^2 (assum4 < 1) "
              f"{['%.4f' % x for x in f['assum4']]}")
        print(f"    ||Mmat(g)z - z|| max {f['z_residual']:.2e}")
        print(f"    formula          {f['perf']:.6f}   terms "
              f"{['%.6f' % x for x in f['terms']]}")
        for d in (400, 800, 1600):
            m = mc(d, c, w, R, TH, reps[d], rng)
            print(f"    d = {d:5d}  MC {m['mean']:.6f} +/- {m['se']:.6f}   "
                  f"dev {m['mean'] - f['perf']:+.6f}")
            rows.append((name, d, f["perf"], m["mean"], m["se"]))

    # ----------------------------------- 3. the suboptimality instance, 3 methods
    sep("3. prop:singleweight_suboptimality  (line 915, proof 2169 to 2194)")
    b0 = betaSq(B_TH0, B_C0)
    print(f"  theta_0 = {B_TH0}, c_0 = {B_C0}, beta_0^2 = {b0:.6f}, 2 beta_0^2 = {2*b0:.6f}")
    print(f"  detectability theta_0^4 = {B_TH0**4:.4f} > c_0 = {B_C0}: "
          f"{B_TH0**4 > B_C0}")
    print("  scale invariance: w_1 = 1 fixed, t = w_2^2 varies; the assumption needs")
    print(f"  t > 1/(1+theta_0^2) = {1/(1+B_TH0**2):.6f}")

    def paper_stacksvd(t, th0=B_TH0, c0=B_C0):
        """The paper's display at line 2190 with w_1 = 1, w_2^2 = t."""
        w1s, w2s = 1.0, t
        d1 = (w1s * (1 + th0**2) - w2s) ** 2
        d2 = (w2s * (1 + th0**2) - w1s) ** 2
        term1 = (th0**4 - c0 - c0 * w2s**2 * th0**4 / d1) / (th0**2 * (1 + th0**2))
        term2 = (th0**4 - c0 - c0 * w1s**2 * th0**4 / d2) / (th0**2 * (1 + th0**2))
        return term1 + term2

    def assum4_B(t, th0=B_TH0, c0=B_C0):
        """The two clause-2 sums at gamma_1 = 1+th0^2, gamma_2 = t(1+th0^2), w_1 = 1."""
        g = [1.0 * (1 + th0**2), t * (1 + th0**2)]
        return [c0 * 1.0 / (gi - 1.0) ** 2 + c0 * t**2 / (gi - t) ** 2 for gi in g]

    ts = np.linspace(1.0 / (1 + B_TH0**2) + 1e-4, 1.0, 4001)
    vals = np.array([paper_stacksvd(t) for t in ts])
    ok = np.array([max(assum4_B(t)) < 1.0 for t in ts])
    tbest = ts[ok][int(np.argmax(vals[ok]))]
    print(f"  sup over admissible t of the paper display: {vals[ok].max():.6f} at "
          f"t = {tbest:.4f} (w_2 = {np.sqrt(tbest):.4f})")
    print(f"  strict?  sup < 2 beta_0^2 : {vals[ok].max() < 2*b0}  gap "
          f"{2*b0 - vals[ok].max():.6f}")
    print(f"  equal weights t = 1 give 2*betaSq(theta_0, 2c_0) = "
          f"{2*betaSq(B_TH0, 2*B_C0):.6f}; the sup over w sits at that boundary,")
    print("  which is the case the paper's display excludes and its last sentence covers.")
    print("  general formula vs the paper display, three weights:")
    for t in (0.75, 0.9, 1.0):
        f = perf_formula([1.0, np.sqrt(t)], B_C, B_R, B_TH)
        p = paper_stacksvd(t)
        print(f"    t = {t:.4f}: general {f['perf']:.9f}  paper {p:.9f}  |diff| "
              f"{abs(f['perf'] - p):.3e}  n_roots {f['n_roots']}")
    print("  term identity (the route recommended for the Lean proof):")
    print("    term_l = beta_0^2 - c_0 w_j^4 theta_0^2 / ((1+theta_0^2)(w_l^2(1+theta_0^2)-w_j^2)^2)")
    worst_t = 0.0
    for t in (0.35, 0.5, 0.75, 0.9, 1.0):
        w = [1.0, np.sqrt(t)]
        for l, j in ((0, 1), (1, 0)):
            num = (
                B_TH0**4 - B_C0
                - B_C0 * w[j] ** 4 * B_TH0**4
                / (w[l] ** 2 * (1 + B_TH0**2) - w[j] ** 2) ** 2
            )
            direct = num / (B_TH0**2 * (1 + B_TH0**2))
            shifted = b0 - B_C0 * w[j] ** 4 * B_TH0**2 / (
                (1 + B_TH0**2) * (w[l] ** 2 * (1 + B_TH0**2) - w[j] ** 2) ** 2
            )
            worst_t = max(worst_t, abs(direct - shifted))
    print(f"    max |paper term - (beta_0^2 - penalty)| over 10 cases: {worst_t:.3e}")
    print("    every penalty is > 0 when the other weight is nonzero, so each term < beta_0^2")

    rng = np.random.default_rng(SEED + 1)
    print("  Monte Carlo of the three methods (svdstack is unweighted = optimally weighted):")
    for d in (400, 800, 1600):
        m_opt = mc(d, B_C, [1.0, np.sqrt(tbest)], B_R, B_TH, reps[d], rng,
                   also_svdstack=True)
        m_eq = mc(d, B_C, [1.0, 1.0], B_R, B_TH, reps[d], rng)
        print(f"    d = {d:5d}  stacksvd(w*) {m_opt['mean']:.4f}+/-{m_opt['se']:.4f} "
              f"(pred {paper_stacksvd(tbest):.4f})   "
              f"stacksvd(1,1) {m_eq['mean']:.4f}+/-{m_eq['se']:.4f} "
              f"(pred {paper_stacksvd(1.0):.4f})   "
              f"svdstack {m_opt['sv_mean']:.4f}+/-{m_opt['sv_se']:.4f} "
              f"(pred {2*b0:.4f})")

    # ------------------------------------------- 4. hypothesis necessity scan
    sep("4. HYPOTHESIS NECESSITY SCAN")

    print("\n  (a) distinctness of the gamma_l  (drop the first clause of assum 2104)")
    f = perf_formula([1.0, 1.0], B_C, B_R, B_TH)
    unw = 2 * betaSq(B_TH0, 2 * B_C0)
    print(f"      instance B at w = (1,1): gamma = {['%.6f' % x for x in f['gamma']]} "
          f"(tied)")
    print(f"      formula {f['perf']:.6f}   unweighted rank-r law 2*betaSq(th,2c) "
          f"{unw:.6f}   |diff| {abs(f['perf'] - unw):.3e}")
    print("      -> the value survives an exact tie HERE only because K(gamma) is a multiple")
    print("         of the identity on the tied plane (theta and c equal across the tables).")

    print("\n      A tie where it does NOT survive: M=2, r=2, R_1=e_1, R_2=e_2,")
    print("      theta^2 = (3, 15), w = (1, 0.5), c = (0.05, 0.05).")
    Rt2 = [np.array([[1.0], [0.0]]), np.array([[0.0], [1.0]])]
    THt2 = [np.array([np.sqrt(3.0)]), np.array([np.sqrt(15.0)])]
    wt2, ct2 = [1.0, 0.5], [0.05, 0.05]
    St2 = signal_blocks(Rt2, THt2)
    gt = 4.0
    Kt = Kmat(gt, wt2, St2)
    print(f"      Mmat(4) = I (double root); Kmat(4) = diag({Kt[0, 0]:.6f}, {Kt[1, 1]:.6f}), "
          f"not a multiple of I")
    for lab, th in (("K-eigenbasis {e1,e2}", 0.0), ("rotated 30 deg", np.pi / 6),
                    ("rotated 45 deg", np.pi / 4)):
        z1 = np.array([np.cos(th), np.sin(th)])
        z2 = np.array([-np.sin(th), np.cos(th)])
        v = sum(eta_of(gt, wt2, ct2) / (gt * float(z @ Kt @ z)) for z in (z1, z2))
        print(f"        {lab:22s}: formula {v:.6f}   z_1^T K z_2 = {float(z1 @ Kt @ z2):+.6f}")
    print(f"        perf_formula (uses the K-eigenbasis): "
          f"{perf_formula(wt2, ct2, Rt2, THt2)['perf']:.6f}")
    rngt = np.random.default_rng(SEED + 7)
    for d in (400, 800, 1600):
        mt = mc(d, ct2, wt2, Rt2, THt2, reps[d], rngt)
        print(f"        MC d = {d:5d}: {mt['mean']:.6f} +/- {mt['se']:.6f}")
    print("      -> the sum is over RECIPROCALS of z^T K z, so it is not rotation invariant.")
    print("         xi_l = (gamma - Sigma)^-1 A z_l are orthogonal iff z_l^T K z_m = 0, so the")
    print("         z_l must be K-orthogonal, not merely orthonormal. Distinctness of the")
    print("         gamma_l is needed for the VALUE, not only for the paper's proof route.")

    print("\n  (b) sup_l sum_i c_i w_i^4/(gamma_l - w_i^2)^2 < 1  (second clause)")
    for c0 in (1.0, 1.9, 2.2):
        THb = [np.array([B_TH0]), np.array([B_TH0])]
        Cb = [c0, c0]
        t = 0.62
        f = perf_formula([1.0, np.sqrt(t)], Cb, B_R, THb)
        print(f"      c_0 = {c0}: gamma {['%.4f' % x for x in f['gamma']]}  assum4 "
              f"{['%.4f' % x for x in f['assum4']]}  formula {f['perf']:.6f}  terms "
              f"{['%.6f' % x for x in f['terms']]}")
    rng = np.random.default_rng(SEED + 2)
    Cb = [2.2, 2.2]
    THb = [np.array([B_TH0]), np.array([B_TH0])]
    fbad = perf_formula([1.0, np.sqrt(0.62)], Cb, B_R, THb)
    for d in (400, 800, 1600):
        m = mc(d, Cb, [1.0, np.sqrt(0.62)], B_R, THb, reps[d], rng)
        print(f"      c_0 = 2.2, t = 0.62, d = {d:5d}: MC {m['mean']:.4f}+/-{m['se']:.4f}"
              f"   formula {fbad['perf']:.4f}   first term only "
              f"{fbad['terms'][0]:.4f}")
    print("      -> when a root violates the clause, that component is in the bulk and the")
    print("         formula's term for it is wrong (it can even be negative).")

    print("\n  (c) gamma_r > sup_i w_i^2  (third clause)")
    for t in (0.40, 0.32, 0.25):
        f = perf_formula([1.0, np.sqrt(t)], B_C, B_R, B_TH)
        print(f"      t = {t:.2f}: roots above wmax^2 = {f['n_roots']} (need r = 2), "
              f"gamma {['%.4f' % x for x in f['gamma']]}, formula sum over the roots "
              f"found {f['perf']:.6f}")
    print(f"      threshold t = 1/(1+theta_0^2) = {1/(1+B_TH0**2):.4f}")
    rng = np.random.default_rng(SEED + 3)
    for d in (400, 800, 1600):
        m = mc(d, B_C, [1.0, np.sqrt(0.25)], B_R, B_TH, reps[d], rng)
        print(f"      t = 0.25, d = {d:5d}: MC {m['mean']:.4f}+/-{m['se']:.4f}   "
              f"one-root formula {perf_formula([1.0,0.5], B_C, B_R, B_TH)['perf']:.4f}   "
              f"beta_0^2 = {b0:.4f}")

    print("\n  (d) Rank(sum_i R_i R_i^T) = r  (assum:unaligned, line 758)")
    Rdeg = [np.array([[1.0], [0.0]]), np.array([[1.0], [0.0]])]
    f = perf_formula([1.0, 0.8], B_C, Rdeg, B_TH)
    print(f"      R_1 = R_2 = e_1, r = 2: rank(sum R_i R_i^T) = "
          f"{np.linalg.matrix_rank(sum(Ri @ Ri.T for Ri in Rdeg))}")
    print(f"      roots above wmax^2 = {f['n_roots']} (need r = 2); Mmat is singular, so "
          f"det(I - Mmat) has\n      only {f['n_roots']} root and the sum over l = 1..r is "
          f"undefined.")

    print("\n  (e) ties inside Theta_i  (assum:unaligned line 770 asks for distinct entries)")
    Rt = [np.eye(2), rot(0.6)]
    THt = [np.array([1.6, 1.6]), np.array([1.7, 1.15])]
    f = perf_formula([1.0, 0.75], A_C, Rt, THt)
    print(f"      Theta_1 = diag(1.6,1.6): gamma {['%.6f' % x for x in f['gamma']]}  "
          f"formula {f['perf']:.6f}  z residual {f['z_residual']:.2e}")
    rng = np.random.default_rng(SEED + 4)
    for d in (400, 800, 1600):
        m = mc(d, A_C, [1.0, 0.75], Rt, THt, reps[d], rng)
        print(f"      d = {d:5d}: MC {m['mean']:.4f}+/-{m['se']:.4f}   dev "
              f"{m['mean'] - f['perf']:+.4f}")

    print("\n  (f) w_i = 0")
    f0 = perf_formula([1.0, 0.0], A_C, A_R, A_TH)
    f1 = perf_formula([1.0], [A_C[0]], [A_R[0]], [A_TH[0]])
    print(f"      instance A with w = (1,0): roots {f0['n_roots']}, formula "
          f"{f0['perf']:.9f}")
    print(f"      table 1 alone           : roots {f1['n_roots']}, formula "
          f"{f1['perf']:.9f}   |diff| {abs(f0['perf']-f1['perf']):.3e}")
    fB0 = perf_formula([1.0, 0.0], B_C, B_R, B_TH)
    print(f"      instance B with w = (1,0): roots {fB0['n_roots']} (need r = 2) -> "
          f"dropping a table can break the rank condition.")
    rng = np.random.default_rng(SEED + 5)
    for d in (400, 800, 1600):
        m = mc(d, A_C, [1.0, 0.0], A_R, A_TH, reps[d], rng)
        print(f"      A, w = (1,0), d = {d:5d}: MC {m['mean']:.4f}+/-{m['se']:.4f}   "
              f"formula {f0['perf']:.4f}")

    sep("5. SUMMARY at d = 1600 (Monte Carlo minus formula)")
    print("    instance                 formula      MC(d=1600)     dev      se")
    for name, d, pred, mean, se in rows:
        if d == 1600:
            print(f"    {name:22s}  {pred:9.6f}   {mean:9.6f}  {mean-pred:+8.6f}  "
                  f"{se:.6f}")
    worst_row = max((abs(m - p_), n) for n, d, p_, m, s_ in rows if d == 1600)
    print(f"    largest |MC - formula| at d = 1600: {worst_row[0]:.6f} ({worst_row[1]})")

    print(f"\nelapsed {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
