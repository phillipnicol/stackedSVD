#!/usr/bin/env python3
"""Numeric check for task U7 of `notes/archive/rankr_plan_A.md` (section 4).

Written 2026-09-01 for the adversarial audit `notes/archive/audit_rankr_plan_A_2026-09-01.md`.
numpy only. Seed 2026090101, printed. Every claim is an assertion, not an eyeball.

Model (plan section 4, and `notes/archive/plan_subspacelaw.md` section 1.1):
  X = U D Vᵀ + E,  E is n x d with iid N(0, 1/d),  D = diag(sqrt(lam)),
  U is n x r orthonormal, V is d x r orthonormal, c = n / d.
  Eperp = (I - U Uᵀ) E,  W0 = Eperpᵀ Eperp,  Q = V D + Eᵀ U.
  Then Xᵀ X = W0 + Q Qᵀ  (checked, cell C0).
  Spike j is supercritical iff lam_j ** 2 > c;  s = number of supercritical spikes.
  bulkEdge b = (1 + sqrt(c)) ** 2,  betaSq(sqrt(lam), c) = (lam^2 - c)/(lam^2 + lam) if
  lam^2 > c else 0,  rhoSq(sqrt(lam), c) = lam + 1 + c + c/lam if lam^2 > c else b.
  The two closed forms are `StackedSVD/Defs.lean:38` and `:46` at theta^2 = lam.

Cells
  C0  the split identity, without which nothing below means anything.
  C1  the count claim inside part 3 of U7: for every index k < r of S = Xᵀ X,
      #{a : mu_a(W0) >= lam_k(S)} <= k + 1 <= r.  Plan calls this the riskiest step.
  C2  gap G1: #{i : lam_i(S) > b + eps} == s, and lamMax(W0 + Qsub Qsubᵀ) <= b + eps.
  C3  U5(b): max_{a<r} ||Qᵀ u_a(W0)||^2 against (sum lam + r)/d.
  C4  ||P_edge v_j||^2 at a subcritical j, and its decay exponent in d.
  C5  the target of the whole track: ||P_top v_j||^2 -> betaSq(sqrt(lam_j), c).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import math
import resource
import time

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np

SEED = 2026090101
C = 1.5
R = 3
B = (1.0 + math.sqrt(C)) ** 2
EPS_LIST = (0.05, 0.10, 0.20)
EPS_MAIN = 0.10

DRAWS = {400: 200, 800: 100, 1600: 40}

SPECTRA = {
    # name -> (core eigenvalues lam_j, deterministic stream offset)
    "mixed_tied": ((4.0, 4.0, 0.5), 1),
    "rank_deficient": ((4.0, 0.0, 0.0), 2),
    "critical": ((math.sqrt(C), 3.0, 0.5), 3),
    "all_supercritical": ((4.0, 3.0, 2.0), 4),
}


def beta_sq(lam, c):
    """betaSq (sqrt lam) c, `StackedSVD/Defs.lean:38` at theta^2 = lam."""
    t4 = lam * lam
    if t4 > c:
        return (t4 - c) / (t4 + lam)
    return 0.0


def rho_sq(lam, c):
    """rhoSq (sqrt lam) c, `StackedSVD/Defs.lean:46` at theta^2 = lam."""
    if lam * lam > c:
        return lam + 1.0 + c + c / lam
    return (1.0 + math.sqrt(c)) ** 2


def n_super(lams, c):
    return sum(1 for x in lams if x * x > c)


def one_draw(rng, d, n):
    """Shared randomness for every spectrum at this draw."""
    E = rng.standard_normal((n, d)) / math.sqrt(d)
    U = np.linalg.qr(rng.standard_normal((n, R)))[0]
    V = np.linalg.qr(rng.standard_normal((d, R)))[0]
    Ep = E - U @ (U.T @ E)
    W0 = Ep.T @ Ep
    mu, uvec = np.linalg.eigh(W0)  # ascending
    mu = mu[::-1].copy()
    uvec = uvec[:, ::-1].copy()  # columns sorted descending
    EtU = E.T @ U
    return E, U, V, W0, mu, uvec, EtU


def run_cell(name, lams, off, d, ndraw):
    n = int(round(C * d))
    lams = np.asarray(lams, dtype=float)
    s = n_super(lams, C)
    sup = [j for j in range(R) if lams[j] ** 2 > C]
    sub = [j for j in range(R) if lams[j] ** 2 <= C]
    Dg = np.sqrt(lams)
    rng = np.random.default_rng(SEED + 1000 * d + off)
    print(f"      [stream] cell={name} d={d} rng seed = {SEED + 1000 * d + off}")

    max_count_all = 0
    max_count_edge = 0
    viol_c1 = 0
    cnt_ok = {e: 0 for e in EPS_LIST}
    w1_ok = {e: 0 for e in EPS_LIST}
    kappa_sum = 0.0
    pedge_sum = np.zeros(R)
    ptop_sum = np.zeros(R)
    lam1_sum = 0.0
    c0_err = None

    for t in range(ndraw):
        E, U, V, W0, mu, uvec, EtU = one_draw(rng, d, n)
        Q = V * Dg[None, :] + EtU
        S = W0 + Q @ Q.T

        if t == 0:
            X = U @ (Dg[:, None] * V.T) + E
            c0_err = float(np.max(np.abs(X.T @ X - S)))

        lamS, vecS = np.linalg.eigh(S)
        lamS = lamS[::-1].copy()
        vecS = vecS[:, ::-1].copy()

        # C1: the count claim of U7 part 3.
        for k in range(R):
            cnt = int(np.count_nonzero(mu >= lamS[k]))
            max_count_all = max(max_count_all, cnt)
            if cnt > k + 1:
                viol_c1 += 1
            if lamS[k] <= B + 2 * EPS_MAIN:
                max_count_edge = max(max_count_edge, cnt)

        # C2: gap G1.
        if sub:
            Qsub = Q[:, sub]
            w1max = float(np.max(np.linalg.eigvalsh(W0 + Qsub @ Qsub.T)))
        else:
            w1max = float(mu[0])
        for e in EPS_LIST:
            if int(np.count_nonzero(lamS > B + e)) == s:
                cnt_ok[e] += 1
            if w1max <= B + e:
                w1_ok[e] += 1

        # C3: U5(b), the delocalization constant kappa.
        QtU = Q.T @ uvec[:, :R]  # r x r, column a is Qᵀ u_a
        kappa_sum += float(np.max(np.sum(QtU ** 2, axis=0)))

        # C4 / C5: the projectors applied to the spike directions.
        top = vecS[:, :R]
        coef = top.T @ V  # r x r, entry (k, j) is <u_k(S), v_j>
        edge_mask = lamS[:R] <= B + EPS_MAIN
        ptop_sum += np.sum(coef ** 2, axis=0)
        pedge_sum += np.sum(coef[edge_mask, :] ** 2, axis=0)
        lam1_sum += float(lamS[0])

    out = {
        "name": name,
        "d": d,
        "n": n,
        "ndraw": ndraw,
        "s": s,
        "sup": sup,
        "sub": sub,
        "c0_err": c0_err,
        "max_count_all": max_count_all,
        "max_count_edge": max_count_edge,
        "viol_c1": viol_c1,
        "cnt_ok": {e: cnt_ok[e] / ndraw for e in EPS_LIST},
        "w1_ok": {e: w1_ok[e] / ndraw for e in EPS_LIST},
        "kappa": kappa_sum / ndraw,
        "kappa_target": (float(np.sum(lams)) + R) / d,
        "pedge": pedge_sum / ndraw,
        "ptop": ptop_sum / ndraw,
        "lam1": lam1_sum / ndraw,
        "rho1": rho_sq(float(lams[0]), C) if lams[0] ** 2 > C else max(
            rho_sq(float(x), C) for x in lams),
    }
    return out


def main():
    t0 = time.time()
    print("=" * 78)
    print("check_edge_count.py   task U7 numeric gate, notes/archive/rankr_plan_A.md section 4")
    print(f"seed = {SEED}   numpy {np.__version__}   c = {C}   r = {R}")
    print(f"bulkEdge b = (1+sqrt c)^2 = {B:.6f}   OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']}")
    print(f"draws per d: {DRAWS}")
    print("=" * 78)

    fails = []
    results = {}

    for name, (lams, off) in SPECTRA.items():
        s = n_super(lams, C)
        print()
        print(f"--- cell {name}: lam = {tuple(round(x, 4) for x in lams)}, "
              f"s = {s} supercritical, betaSq = "
              f"{tuple(round(beta_sq(x, C), 4) for x in lams)}")
        for d in sorted(DRAWS):
            res = run_cell(name, lams, off, d, DRAWS[d])
            results[(name, d)] = res

            # C0 split identity
            ok0 = res["c0_err"] < 1e-10
            print(f"  d={d:5d} n={res['n']:5d} draws={res['ndraw']:4d}  "
                  f"C0 split identity max|XtX - W0 - QQt| = {res['c0_err']:.3e}  "
                  f"[{'PASS' if ok0 else 'FAIL'}]")
            if not ok0:
                fails.append(f"C0 {name} d={d} err={res['c0_err']:.3e}")

            # C1 count
            ok1 = res["viol_c1"] == 0
            print(f"          C1 count  max #{{a: mu_a >= lam_k}} over k<r = "
                  f"{res['max_count_all']} (edge-only {res['max_count_edge']}), "
                  f"violations of (<= k+1) = {res['viol_c1']} / {res['ndraw'] * R}  "
                  f"[{'PASS' if ok1 else 'FAIL'}]")
            if not ok1:
                fails.append(f"C1 {name} d={d} violations={res['viol_c1']}")

            # C2 G1
            for e in EPS_LIST:
                fc = res["cnt_ok"][e]
                fw = res["w1_ok"][e]
                bar = "PASS" if fc == 1.0 else "FAIL"
                print(f"          C2 G1 eps={e:.2f}  P[#(lam_i > b+eps) == s] = {fc:.3f}  "
                      f"P[lamMax(W0+QsubQsubt) <= b+eps] = {fw:.3f}  "
                      f"[plan bar {bar}]")
                if fc != 1.0:
                    fails.append(f"C2 {name} d={d} eps={e} count-rate={fc:.3f}")
                if fw != 1.0:
                    fails.append(f"C2w {name} d={d} eps={e} w1-rate={fw:.3f}")

            # C3 U5b
            ratio = res["kappa"] / res["kappa_target"]
            ok3 = ratio < 5.0
            print(f"          C3 U5b  mean max_a ||Qt u_a||^2 = {res['kappa']:.5f}  "
                  f"target (sum lam + r)/d = {res['kappa_target']:.5f}  "
                  f"ratio {ratio:.2f}  [{'PASS' if ok3 else 'FAIL'}]")
            if not ok3:
                fails.append(f"C3 {name} d={d} ratio={ratio:.2f}")

            # C4 / C5
            pe = ", ".join(f"{x:.5f}" for x in res["pedge"])
            pt = ", ".join(f"{x:.5f}" for x in res["ptop"])
            tgt = ", ".join(f"{beta_sq(x, C):.5f}" for x in lams)
            print(f"          C4 mean ||P_edge v_j||^2 = [{pe}]")
            print(f"          C5 mean ||P_top  v_j||^2 = [{pt}]  target betaSq [{tgt}]")

    # C4 decay exponent and C5 convergence, across d
    print()
    print("=" * 78)
    print("cross-d trends")
    ds = sorted(DRAWS)
    for name, (lams, off) in SPECTRA.items():
        sub = [j for j in range(R) if lams[j] ** 2 <= C]
        if sub:
            j = sub[0]
            v_lo = results[(name, ds[0])]["pedge"][j]
            v_hi = results[(name, ds[-1])]["pedge"][j]
            expo = math.log(v_lo / v_hi) / math.log(ds[-1] / ds[0]) if v_hi > 0 else float("inf")
            okd = v_hi < v_lo
            print(f"  {name:18s} C4 ||P_edge v_{j}||^2: {v_lo:.5f} (d={ds[0]}) -> "
                  f"{v_hi:.5f} (d={ds[-1]})   fitted exponent {expo:.2f}  "
                  f"[{'decaying' if okd else 'NOT DECAYING'}]")
            if not okd:
                fails.append(f"C4 {name} not decaying")
        errs = []
        for d in ds:
            e = float(np.max(np.abs(results[(name, d)]["ptop"]
                                    - np.array([beta_sq(x, C) for x in lams]))))
            errs.append(e)
        ok5 = errs[-1] < errs[0]
        print(f"  {name:18s} C5 max_j |mean ||P_top v_j||^2 - betaSq_j| = "
              f"{[round(x, 5) for x in errs]}  [{'converging' if ok5 else 'NOT CONVERGING'}]")
        if not ok5:
            fails.append(f"C5 {name} not converging")

    print()
    print("=" * 78)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0
    print(f"wall {time.time() - t0:.1f} s   maxrss {rss:.2f} GB   seed {SEED}")
    if fails:
        print(f"FAILURES ({len(fails)}):")
        for f in fails:
            print(f"  {f}")
    else:
        print("all assertions PASS")
    print("=" * 78)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
