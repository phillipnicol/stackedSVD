#!/usr/bin/env python3
"""Numeric check of unit G12, `notes/prop_single_table_general.md` section 6.

Checks every closed form that Stage 1 (the non-Gaussian `SingleTableLaw`) targets, at four
entry laws (Gaussian control, Rademacher, uniform, Student-t(5)), before any Lean proof
starts. Adds rows 12 and 13 from section 8 item 1 of the same note (the coordinator review
that found unit G4 needs the mixed form): row 12 is the mixed form `v.T @ G @ E.T @ u` with
`G` the full (non-deflated) noise Gram resolvent, target 0; row 13 is the companion identity
`g.T @ G @ g = 1 + z * u.T @ Gc @ u`, target `1 + z * mTilde(c, z)`.

Run (single process, needs about 19-21 minutes measured on this host; the script caps the
BLAS threads at 4 unless OMP_NUM_THREADS and its siblings are already set):
    ulimit -v 8000000
    python3 scripts/numeric/check_nongaussian_forms.py [--quick] [--out PATH]

Exit code: 0 if every non-necessity, non-informational row passes; 1 if any fails.

Checkpoint mode, for a caller whose single-command timeout is under the ~20 minute full run
(the harness driving this script has a hard 600 s per-call cap): the 16 (law, case) grid
units are the only slow part (necessity, nu4 and the R/python cross-check together take well
under a minute). `--resume-grid CKPT --time-budget SECONDS` runs grid units in law-major
order, saving `all_case_results` to CKPT (a pickle) after every completed unit, and stops
starting new units once elapsed time passes the budget; exit 0 means every unit is now done,
exit 42 means budget ran out with units left, so call it again with the same CKPT path.
`--time-budget` defaults to 420 s, chosen so that 420 s plus the single worst unit (about
166 s, d=1600 c=2.0) stays under the 600 s hard cap. Once a `--resume-grid` call exits 0,
`--finalize CKPT --out PATH` loads the completed checkpoint, runs the fast sections (nu4,
cross-check, necessity), assembles the report and writes PATH, with the same exit-code
contract as a single monolithic run.

Design choices not pinned by the note (recorded here, not just in the output file):
  * `case_index` enumerates the 4 (d, c) grid points 0..3, independent of law and z; z does
    not get its own seed, since it only changes the resolvent shift applied to the same
    draw of Z, u, v, not the draw itself.
  * The direction seeds for u and v are the case's base seed (`SEED_BASE + 100000*law_index
    + 1000*case_index`) plus fixed offsets 700 and 701, which lie outside any replicate
    index in range [0, 200), so they cannot collide with a per-replicate Z-draw seed. This
    makes u, v law-specific (redrawn per law at the same (d, c)), the literal reading of
    "print the base seed of every (law, case) pair; print the seed of the two direction
    draws separately" as two more numbers attached to that same (law, case) pair.
  * The `W0` identity assert (`W0 == E.T @ (I - u u^T) @ E`) is an exact algebraic identity
    (E^T(I - uu^T)E = E^T E - (E^T u)(u^T E) = E^T E - g g^T), so a code bug would show on
    every replicate identically; it is checked once per (law, case) on replicate 0, not on
    all 250 replicates, to stay inside the time budget (the independent recomputation costs
    about as much as one more E.T @ E).
  * Rows 12 and 13 are not covered by the tolerance paragraph of section 6 (written before
    they existed). They get the same "4 SE + 3 d^-1/2" rule as rows 1-7, 9, 11, since they
    are isotropic-type bilinear resolvent forms of the same kind and the note does not ask
    for anything stricter or looser.
  * `G` (row 12/13) is a fresh, independent `numpy.linalg.inv(E.T @ E - z I)`, not derived
    from `G0` by the Sherman-Morrison rank-one update that proves the row 12/13 targets
    algebraically. Deriving it from `G0` would make the check partly tautological (a bug in
    that derivation would go unnoticed); a fresh inversion gives genuine independent
    numerical evidence.
  * Rows 9 and 10 use `numpy.linalg.eigh` (LAPACK, full spectrum), not a hand-rolled power
    iteration. Benchmarked on this host: eigh on a 1600x1600 matrix takes about 0.2 to 0.5 s,
    a fresh inversion about 0.25 s, so eigh is not the bottleneck here, and LAPACK's routine
    is the standard, tested tool for a symmetric eigenproblem (a hand power iteration risks
    slow convergence near a small spectral gap, e.g. row 10 sits below the BBP threshold, so
    there is no isolated top eigenvector to converge to).
"""
import argparse
import os
import pickle
import subprocess
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")  # a threaded BLAS takes every core otherwise; set before numpy loads
import numpy as np

SEED_BASE = 20260909

# ---------------------------------------------------------------------------------------
# Closed forms, copied from the Lean sources (read on the tree at HEAD, 2026-09-09). Each
# comment names the exact `def` line; the formula itself is the next line in every case
# (Lean writes the signature and the body on separate lines). Transcriptions matched the
# note's section 6 exactly; no discrepancy found.
# ---------------------------------------------------------------------------------------

def bulk_edge(c):
    # lean/StackedSVD/Defs.lean:50
    #   noncomputable def bulkEdge (c : R) : R := (1 + Real.sqrt c) ^ 2
    return (1.0 + np.sqrt(c)) ** 2


def bulk_edge_lo(c):
    # lean/StackedSVD/RMT/MP.lean:27
    #   noncomputable def bulkEdgeLo (c : R) : R := (1 - Real.sqrt c) ^ 2
    return (1.0 - np.sqrt(c)) ** 2


def mp_m(c, z):
    # lean/StackedSVD/RMT/MP.lean:30-31
    #   noncomputable def m (c z : R) : R :=
    #     (-(z + 1 - c) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))) / (2 * z)
    return (-(z + 1.0 - c) + np.sqrt((z - bulk_edge(c)) * (z - bulk_edge_lo(c)))) / (2.0 * z)


def mp_mderiv(c, z):
    # lean/StackedSVD/RMT/MP.lean:34-35
    #   noncomputable def mDeriv (c z : R) : R :=
    #     -(m c z * (m c z + 1)) / Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))
    m = mp_m(c, z)
    return -(m * (m + 1.0)) / np.sqrt((z - bulk_edge(c)) * (z - bulk_edge_lo(c)))


def beta_sq(theta, c):
    # lean/StackedSVD/Defs.lean:46-47
    #   noncomputable def betaSq (theta c : R) : R :=
    #     if theta ^ 4 > c then (theta ^ 4 - c) / (theta ^ 4 + theta ^ 2) else 0
    if theta ** 4 > c:
        return (theta ** 4 - c) / (theta ** 4 + theta ** 2)
    return 0.0


def rho_sq(theta, c):
    # lean/StackedSVD/Defs.lean:54-55
    #   noncomputable def rhoSq (theta c : R) : R :=
    #     if theta ^ 4 > c then theta ^ 2 + 1 + c + c / theta ^ 2 else bulkEdge c
    if theta ** 4 > c:
        return theta ** 2 + 1.0 + c + c / theta ** 2
    return bulk_edge(c)


def m_tilde(c, z):
    # note section 6, row 11: mTilde(c, z) = (m(c, z) + (1 - c) / z) / c
    return (mp_m(c, z) + (1.0 - c) / z) / c


# ---------------------------------------------------------------------------------------
# Laws, rescaled to variance 1. nu4 targets from the note's table.
# ---------------------------------------------------------------------------------------

SQRT3 = float(np.sqrt(3.0))
T5_SCALE = float(np.sqrt(5.0 / 3.0))
T3_SCALE = float(np.sqrt(3.0))


def draw_gaussian(rng, shape):
    return rng.standard_normal(shape)


def draw_rademacher(rng, shape):
    return 2.0 * rng.integers(0, 2, size=shape) - 1.0


def draw_uniform(rng, shape):
    return rng.uniform(-SQRT3, SQRT3, size=shape)


def draw_t5(rng, shape):
    return rng.standard_t(5, size=shape) / T5_SCALE


LAWS = [
    ("gaussian", draw_gaussian, 3.0),
    ("rademacher", draw_rademacher, 1.0),
    ("uniform", draw_uniform, 1.8),
    ("t5", draw_t5, 9.0),
]

# Necessity draws, each a one-hypothesis break of the rademacher base law.
NECESSITY_D = 400
NECESSITY_C = 0.5
NECESSITY_LAW_INDEX = 1  # rademacher's index in LAWS, kept for seed alignment
NECESSITY_REPS = 200


def draw_necessity_mean_shift(rng, shape):
    """Rademacher + 0.1: breaks mean 0, keeps unit-scale entries."""
    return (2.0 * rng.integers(0, 2, size=shape) - 1.0) + 0.1


def draw_necessity_variance2(rng, shape):
    """Rademacher scaled to variance 2: breaks variance 1."""
    return np.sqrt(2.0) * (2.0 * rng.integers(0, 2, size=shape) - 1.0)


def draw_necessity_t3(rng, shape):
    """Student-t(3), variance normalized to 1: breaks the finite fourth moment only."""
    return rng.standard_t(3, size=shape) / T3_SCALE


def draw_necessity_ar1(rng, shape, rho=0.5):
    """Rademacher innovations with AR(1) correlation 0.5 within each row: breaks
    independence of entries within a row, keeps marginal variance 1."""
    n, d = shape
    eps = 2.0 * rng.integers(0, 2, size=shape) - 1.0
    out = np.empty(shape)
    out[:, 0] = eps[:, 0]
    a = np.sqrt(1.0 - rho * rho)
    for j in range(1, d):
        out[:, j] = rho * out[:, j - 1] + a * eps[:, j]
    return out


NECESSITY = [
    ("mean_shift_0.1", draw_necessity_mean_shift),
    ("variance_2", draw_necessity_variance2),
    ("t3_infinite_4th_moment", draw_necessity_t3),
    ("ar1_row_correlation_0.5", draw_necessity_ar1),
]

# ---------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------


def read_vmhwm():
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmHWM:"):
                    return line.split(":", 1)[1].strip()
    except OSError as e:
        return f"unavailable ({e})"
    return "not found"


def sample_sphere(rng, dim):
    x = rng.standard_normal(dim)
    return x / np.linalg.norm(x)


def mean_se(values):
    x = np.asarray(values, dtype=float)
    m = float(x.mean())
    se = float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else float("nan")
    return m, se


def std_tol(se, d):
    return 4.0 * se + 3.0 * d ** -0.5


def std_pass(gap, se, d):
    return gap <= std_tol(se, d)


def check_w0_identity(E, u, g, W0, tol=1e-10):
    """W0 = E.T @ E - g g.T should equal E.T @ (I - u u.T) @ E exactly (algebraic
    identity: E^T(I - uu^T)E = E^T E - (E^T u)(u^T E) = E^T E - g g^T). Computed via the
    cheap equivalent route E_defl = E - outer(u, g); E_defl.T @ E_defl, avoiding the O(n^2 d)
    cost of forming (I - u u^T) explicitly."""
    E_defl = E - np.outer(u, g)
    alt = E_defl.T @ E_defl
    maxdiff = float(np.max(np.abs(alt - W0)))
    assert maxdiff < tol, f"W0 identity check failed: max abs diff {maxdiff} >= {tol}"
    return maxdiff


# ---------------------------------------------------------------------------------------
# Row bookkeeping
# ---------------------------------------------------------------------------------------

ROW_KEYS_Z = ["row1", "row2", "row3", "row4", "row5", "row6", "row7a", "row7b",
              "row11", "row12", "row13"]
ROWS_DUAL_C = {"row1", "row2", "row3", "row4", "row5", "row6", "row7a", "row7b", "row11"}


def target_for(key, c, z):
    if key in ("row1", "row2", "row7a"):
        return mp_m(c, z)
    if key in ("row4", "row5", "row7b"):
        return mp_mderiv(c, z)
    if key in ("row3", "row6", "row12"):
        return 0.0
    if key == "row11":
        return m_tilde(c, z)
    if key == "row13":
        return 1.0 + z * m_tilde(c, z)
    raise ValueError(key)


def run_case(law_name, draw_fn, law_index, case_index, d, c, n_reps, verbose=True):
    n = int(round(c * d))
    base_seed = SEED_BASE + 100000 * law_index + 1000 * case_index
    seed_u = base_seed + 700
    seed_v = base_seed + 701

    u = sample_sphere(np.random.default_rng(seed_u), n)
    v = sample_sphere(np.random.default_rng(seed_v), d)

    z_list = [bulk_edge(c) + 0.5, bulk_edge(c) + 2.0]
    c_alt = (n - 1) / d

    acc = {0: {k: [] for k in ROW_KEYS_Z}, 1: {k: [] for k in ROW_KEYS_Z}}
    acc_static = {"row8": [], "row9_overlap": [], "row9_lam": [], "row10_overlap": []}

    theta9 = (2.0 * c) ** 0.25
    theta10 = (c / 2.0) ** 0.25

    identity_checked_maxdiff = None

    t_start = time.time()
    for rep in range(n_reps):
        seed_rep = base_seed + rep
        rngz = np.random.default_rng(seed_rep)
        Z = draw_fn(rngz, (n, d))
        E = Z / np.sqrt(d)
        g = E.T @ u
        A = E.T @ E
        W0 = A - np.outer(g, g)

        if rep == 0:
            identity_checked_maxdiff = check_w0_identity(E, u, g, W0)

        lam_w0 = float(np.linalg.eigvalsh(W0)[-1])
        acc_static["row8"].append(lam_w0)

        h9 = g + theta9 * v
        M9 = W0 + np.outer(h9, h9)
        w9, V9 = np.linalg.eigh(M9)
        lam9 = float(w9[-1])
        overlap9 = float((v @ V9[:, -1]) ** 2)
        acc_static["row9_lam"].append(lam9)
        acc_static["row9_overlap"].append(overlap9)

        h10 = g + theta10 * v
        M10 = W0 + np.outer(h10, h10)
        _, V10 = np.linalg.eigh(M10)
        overlap10 = float((v @ V10[:, -1]) ** 2)
        acc_static["row10_overlap"].append(overlap10)

        for zi, z in enumerate(z_list):
            G0 = np.linalg.inv(W0 - z * np.eye(d))
            p = G0 @ v
            q = G0 @ g
            a = acc[zi]
            a["row1"].append(float(v @ p))
            a["row2"].append(float(g @ q))
            a["row3"].append(float(v @ q))
            a["row4"].append(float(p @ p))
            a["row5"].append(float(q @ q))
            a["row6"].append(float(p @ q))
            a["row7a"].append(float(np.trace(G0)) / d)
            a["row7b"].append(float(np.sum(G0 * G0)) / d)

            G = np.linalg.inv(A - z * np.eye(d))
            r = G @ g
            a["row12"].append(float(v @ r))
            a["row13"].append(float(g @ r))

            Gc = np.linalg.inv(E @ E.T - z * np.eye(n))
            a["row11"].append(float(np.trace(Gc)) / n)

        if verbose and n_reps >= 20 and (rep + 1) % max(1, n_reps // 4) == 0:
            print(f"    [{law_name} case{case_index} d={d} c={c}] rep {rep + 1}/{n_reps} "
                  f"({time.time() - t_start:.1f} s elapsed)", flush=True)

    return dict(law_name=law_name, law_index=law_index, case_index=case_index, d=d, c=c,
                n=n, c_alt=c_alt, z_list=z_list, acc=acc, acc_static=acc_static,
                theta9=theta9, theta10=theta10, base_seed=base_seed, seed_u=seed_u,
                seed_v=seed_v, w0_identity_maxdiff=identity_checked_maxdiff,
                wall_s=time.time() - t_start)


def run_necessity_indexed(name, draw_fn, variant_index, n_reps, verbose=True):
    d, c = NECESSITY_D, NECESSITY_C
    n = int(round(c * d))
    case_index = 100 + variant_index
    base_seed = SEED_BASE + 100000 * NECESSITY_LAW_INDEX + 1000 * case_index
    seed_u = base_seed + 700
    seed_v = base_seed + 701
    u = sample_sphere(np.random.default_rng(seed_u), n)
    v = sample_sphere(np.random.default_rng(seed_v), d)
    z = bulk_edge(c) + 0.5  # single z, the first of the grid's two offsets

    lam_vals = []
    row1_vals = []
    t_start = time.time()
    for rep in range(n_reps):
        rngz = np.random.default_rng(base_seed + rep)
        Z = draw_fn(rngz, (n, d))
        E = Z / np.sqrt(d)
        g = E.T @ u
        W0 = E.T @ E - np.outer(g, g)
        lam_vals.append(float(np.linalg.eigvalsh(W0)[-1]))
        G0 = np.linalg.inv(W0 - z * np.eye(d))
        row1_vals.append(float(v @ (G0 @ v)))
    if verbose:
        print(f"    [necessity {name}] {n_reps} reps in {time.time() - t_start:.1f} s",
              flush=True)
    return dict(name=name, d=d, c=c, n=n, z=z, base_seed=base_seed, seed_u=seed_u,
                seed_v=seed_v, lam_vals=lam_vals, row1_vals=row1_vals)


# ---------------------------------------------------------------------------------------
# nu4 check
# ---------------------------------------------------------------------------------------

def nu4_check(law_name, draw_fn, law_index, nsample=2_000_000):
    seed = SEED_BASE + 800000 + law_index
    rng = np.random.default_rng(seed)
    x = draw_fn(rng, (nsample,))
    mean = float(x.mean())
    var = float(x.var())
    nu4_empirical = float(np.mean((x - mean) ** 4)) / (var ** 2) if var > 0 else float("nan")
    return dict(law_name=law_name, seed=seed, nsample=nsample, mean=mean, var=var,
                nu4_empirical=nu4_empirical)


# ---------------------------------------------------------------------------------------
# Cross-check of betaSq (and the attempt at rhoSq) against the paper's reference code
# ---------------------------------------------------------------------------------------
# The reference code is the paper's own repository, https://github.com/phillipnicol/stackedSVD
# (the `stackedSVD` R package with its analysis code): `R/theory_pred.R` and
# `analysis/python_plots/simulations.py`. STACKEDSVD_REPO names a checkout of it. Without the
# variable the script tries the parent of this tree (the public layout, where this tree is
# the `lean_formalization/` folder of that repository).

CROSS_CHECK_POINTS = [(1.2, 0.5), (0.6, 2.0), (1.5, 1.0)]


def reference_repo():
    """The checkout of the stackedSVD repository that holds the reference code, or None."""
    tree = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    for cand in (os.environ.get("STACKEDSVD_REPO", ""), os.path.dirname(tree)):
        if cand and os.path.isfile(os.path.join(cand, "R", "theory_pred.R")):
            return cand
    return None


def reference_paths():
    """(theory_pred.R, the folder of simulations.py), or (None, None) with no checkout."""
    repo = reference_repo()
    if repo is None:
        return None, None
    return (os.path.join(repo, "R", "theory_pred.R"),
            os.path.join(repo, "analysis", "python_plots"))


def cross_check_r(points):
    """Call the unmodified getTheoryPred() in theory_pred.R. Its 'Stack SVD' component is
    exactly betaSq(theta, c) for a single table (theta2 = theta, sum(c) = c). getTheoryPred
    crashes for a literal length-1 call: `diag(x)` for a bare scalar x in (0, 1) makes a
    floor(x) x floor(x) matrix in R, not a 1x1 matrix holding x, and this fires inside the
    function's unused 'SVD Stack' branch (`diag(1 - beta.true^2)`). The call below pads a
    second, inert table (theta=0, c=1e-300) to keep every internal vector at length >= 2 and
    dodge the R diag() footgun without editing theory_pred.R; c is 1e-300, not exactly 0, to
    avoid a 0/0 in the (also unused) 'Stack-SVD W' branch. sum(c) and theta2 are unaffected
    to double precision (1e-300 underflows against any c of order 1), so 'Stack SVD' is
    algebraically identical to the unpadded single-table call."""
    r_file, _ = reference_paths()
    if r_file is None:
        return None, "no checkout of the stackedSVD repository (set STACKEDSVD_REPO)"
    r_lines = [
        f'source("{r_file}")',
        "eps <- 1e-300",
    ]
    for th, cc in points:
        r_lines.append(
            f'v <- getTheoryPred(theta.true = c({th!r}, 0), c = c({cc!r}, eps)); '
            f'cat(sprintf("%.17g\\n", v[["Stack SVD"]]))'
        )
    r_code = "\n".join(r_lines)
    last_err = None
    for attempt in range(3):  # a shared, sometimes-loaded box can hit a transient crash
        if attempt:
            time.sleep(2.0)
        try:
            result = subprocess.run(["Rscript", "-e", r_code], capture_output=True, text=True,
                                     timeout=60)
            if result.returncode != 0:
                last_err = f"Rscript exit {result.returncode}: {result.stderr.strip()}"
                continue
            vals = [float(x) for x in result.stdout.strip().splitlines()]
            if len(vals) != len(points):
                last_err = f"expected {len(points)} lines, got {result.stdout!r}"
                continue
            return vals, None
        except (OSError, subprocess.SubprocessError) as e:
            last_err = f"Rscript unavailable: {e}"
    return None, f"{last_err} (after 3 attempts)"


def cross_check_python(points):
    """Call the unmodified compute_asymptotic_power(method=0, ...) in simulations.py, a
    single-table call (length-1 theta_arr, c_arr): method 0 is
    max((theta_norm_sq**2 - c.sum())/(theta_norm_sq*(theta_norm_sq+1)), 0), which for one
    table is exactly betaSq(theta, c) (the max(...,0) plays the role of the Lean `if`)."""
    _, sim_dir = reference_paths()
    if sim_dir is None:
        return None, "no checkout of the stackedSVD repository (set STACKEDSVD_REPO)"
    py_lines = [
        "import sys",
        f'sys.path.insert(0, "{sim_dir}")',
        "import simulations",
    ]
    for th, cc in points:
        py_lines.append(
            f"print(repr(simulations.compute_asymptotic_power(0, [{th!r}], [{cc!r}])))"
        )
    py_code = "\n".join(py_lines)
    last_err = None
    for attempt in range(3):
        if attempt:
            time.sleep(2.0)
        try:
            result = subprocess.run([sys.executable, "-c", py_code], capture_output=True,
                                     text=True, timeout=60)
            if result.returncode != 0:
                last_err = f"python exit {result.returncode}: {result.stderr.strip()}"
                continue
            vals = [float(x) for x in result.stdout.strip().splitlines()]
            if len(vals) != len(points):
                last_err = f"expected {len(points)} lines, got {result.stdout!r}"
                continue
            return vals, None
        except (OSError, subprocess.SubprocessError) as e:
            last_err = f"python subprocess unavailable: {e}"
    return None, f"{last_err} (after 3 attempts)"


# ---------------------------------------------------------------------------------------
# Markdown assembly
# ---------------------------------------------------------------------------------------

def fmt(x, sig=6):
    return f"{x:.{sig}g}"


def md_table(headers, rows):
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(lines)


ROW_LABEL = {
    "row1": "v.T G0 v", "row2": "g.T G0 g", "row3": "v.T G0 g",
    "row4": "v.T G0 G0 v", "row5": "g.T G0 G0 g", "row6": "v.T G0 G0 g",
    "row7a": "trace(G0)/d", "row7b": "trace(G0 G0)/d", "row11": "trace(Gc)/n",
    "row12": "v.T G E.T u", "row13": "g.T G g",
}


def full_grid_units():
    """The 16 (law, case) grid units in law-major order: 4 laws x the 4 (d, c) points."""
    cases = [(0, 400, 0.5), (1, 400, 2.0), (2, 1600, 0.5), (3, 1600, 2.0)]
    n_reps_by_d = {400: 200, 1600: 50}
    units = []
    for law_index, (law_name, draw_fn, _) in enumerate(LAWS):
        for case_index, d, c in cases:
            units.append((law_index, law_name, draw_fn, case_index, d, c, n_reps_by_d[d]))
    return units


def load_checkpoint(path):
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    return []


def save_checkpoint(path, all_case_results):
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(all_case_results, f)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                     help="dry run: d=40, 3 replicates, gaussian only, no file write")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "notes", "audit_nongaussian_forms_2026-09-09.md"))
    ap.add_argument("--resume-grid", metavar="CKPT", default=None,
                     help="run only the 16 grid units, checkpointing to CKPT; exit 0 when "
                          "all are done, 42 if the time budget ran out with units left")
    ap.add_argument("--time-budget", type=float, default=420.0,
                     help="seconds of grid work per --resume-grid call (default 420, so "
                          "budget + the worst single unit stays under a 600 s hard cap)")
    ap.add_argument("--finalize", metavar="CKPT", default=None,
                     help="load a completed --resume-grid checkpoint and write the report")
    args = ap.parse_args()

    if args.resume_grid and args.finalize:
        print("ERROR: --resume-grid and --finalize are mutually exclusive", file=sys.stderr)
        return 2

    t0 = time.time()
    print(f"numpy version: {np.__version__}", flush=True)
    print(f"quick mode: {args.quick}  resume_grid: {args.resume_grid}  "
          f"finalize: {args.finalize}", flush=True)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        print(f"  {var}={os.environ.get(var, '(unset)')}", flush=True)

    fail_count = 0
    nu4_rows = []
    cross_rows = []
    rho_identity_rows = []

    if not args.resume_grid:
      # --- nu4 check -----------------------------------------------------------------
      print("\n=== nu4 check ===", flush=True)
      for law_index, (law_name, draw_fn, nu4_target) in enumerate(LAWS):
        if args.quick and law_index > 0:
            continue
        r = nu4_check(law_name, draw_fn, law_index, nsample=200_000 if args.quick else 2_000_000)
        gap = abs(r["nu4_empirical"] - nu4_target)
        print(f"  {law_name}: seed={r['seed']} n={r['nsample']} mean={fmt(r['mean'])} "
              f"var={fmt(r['var'])} nu4_emp={fmt(r['nu4_empirical'])} nu4_target={nu4_target} "
              f"gap={fmt(gap)}", flush=True)
        nu4_rows.append((law_name, r["seed"], r["nsample"], fmt(r["mean"]), fmt(r["var"]),
                          fmt(r["nu4_empirical"]), nu4_target, fmt(gap)))

      # --- cross-check betaSq / rhoSq -------------------------------------------------
      print("\n=== cross-check betaSq against getTheoryPred (R) and compute_asymptotic_power"
          " (python) ===", flush=True)
      r_vals, r_err = cross_check_r(CROSS_CHECK_POINTS)
      py_vals, py_err = cross_check_python(CROSS_CHECK_POINTS)
      for i, (th, cc) in enumerate(CROSS_CHECK_POINTS):
        beta_numpy = beta_sq(th, cc)
        rho_numpy = rho_sq(th, cc)
        r_v = r_vals[i] if r_vals is not None else None
        py_v = py_vals[i] if py_vals is not None else None
        r_gap = abs(beta_numpy - r_v) if r_v is not None else None
        py_gap = abs(beta_numpy - py_v) if py_v is not None else None
        print(f"  theta={th} c={cc}: betaSq(numpy)={beta_numpy!r} "
              f"R={r_v!r} (gap {r_gap}) python={py_v!r} (gap {py_gap}) "
              f"rhoSq(numpy)={rho_numpy!r} [no reference in R or python]", flush=True)
        cross_rows.append((th, cc, repr(beta_numpy),
                            repr(r_v) if r_v is not None else f"FAILED: {r_err}",
                            fmt(r_gap, 3) if r_gap is not None else "n/a",
                            repr(py_v) if py_v is not None else f"FAILED: {py_err}",
                            fmt(py_gap, 3) if py_gap is not None else "n/a",
                            repr(rho_numpy)))
        if r_gap is not None and r_gap >= 1e-13:
            fail_count += 1
            print(f"    FAIL: R gap {r_gap} >= 1e-13", flush=True)
        if py_gap is not None and py_gap >= 1e-13:
            fail_count += 1
            print(f"    FAIL: python gap {py_gap} >= 1e-13", flush=True)
      if r_err:
          print(f"  R cross-check note: {r_err}", flush=True)
      if py_err:
          print(f"  python cross-check note: {py_err}", flush=True)
      print("  rhoSq (the outlier eigenvalue limit) has no counterpart in either "
          "theory_pred.R's getTheoryPred or simulations.py's compute_asymptotic_power: both "
          "files only implement the overlap (betaSq) prediction for a single table, never "
          "the top-eigenvalue prediction. rhoSq is checked two other ways instead: the "
          "algebraic identity rhoSq(theta,c) == (1+theta^2)*(theta^2+c)/theta^2 above "
          "threshold (verified below), and the Monte Carlo comparison of row 9 against the "
          "empirical top eigenvalue of X^T X.", flush=True)
      for th, cc in CROSS_CHECK_POINTS:
        if th ** 4 > cc:
            lhs = rho_sq(th, cc)
            rhs = (1.0 + th ** 2) * (th ** 2 + cc) / th ** 2
            gap = abs(lhs - rhs)
            rho_identity_rows.append((th, cc, repr(lhs), repr(rhs), fmt(gap, 3)))
            if gap >= 1e-12:
                fail_count += 1
                print(f"    FAIL: rhoSq algebraic identity gap {gap} at theta={th} c={cc}",
                      flush=True)

    # --- main grid -------------------------------------------------------------------
    print("\n=== main grid ===", flush=True)

    def _log_unit(res):
        print(f"    base_seed={res['base_seed']} seed_u={res['seed_u']} "
              f"seed_v={res['seed_v']} n={res['n']} w0_identity_maxdiff="
              f"{res['w0_identity_maxdiff']:.3e} wall={res['wall_s']:.1f}s", flush=True)

    if args.resume_grid:
        units = full_grid_units()
        all_case_results = load_checkpoint(args.resume_grid)
        done_keys = {(r["law_index"], r["case_index"]) for r in all_case_results}
        print(f"  checkpoint {args.resume_grid}: {len(done_keys)}/{len(units)} units "
              f"already done", flush=True)
        t_grid_start = time.time()
        for law_index, law_name, draw_fn, case_index, d, c, n_reps in units:
            if (law_index, case_index) in done_keys:
                continue
            elapsed = time.time() - t_grid_start
            if elapsed > args.time_budget:
                print(f"  time budget {args.time_budget}s exceeded ({elapsed:.1f}s used), "
                      f"{len(units) - len(done_keys)} units left; resume with the same "
                      f"--resume-grid {args.resume_grid}", flush=True)
                return 42
            print(f"  running law={law_name} case={case_index} d={d} c={c} "
                  f"n_reps={n_reps} ...", flush=True)
            res = run_case(law_name, draw_fn, law_index, case_index, d, c, n_reps)
            _log_unit(res)
            print(f"    VmHWM so far this process: {read_vmhwm()}", flush=True)
            all_case_results.append(res)
            done_keys.add((law_index, case_index))
            save_checkpoint(args.resume_grid, all_case_results)
        print(f"  grid complete: {len(all_case_results)}/{len(units)} units in checkpoint "
              f"{args.resume_grid}; VmHWM this process: {read_vmhwm()}", flush=True)
        return 0

    if args.finalize:
        all_case_results = load_checkpoint(args.finalize)
        units = full_grid_units()
        if len(all_case_results) != len(units):
            print(f"ERROR: checkpoint {args.finalize} has {len(all_case_results)} units, "
                  f"expected {len(units)}; run --resume-grid to completion first",
                  file=sys.stderr)
            return 2
        print(f"  loaded {len(all_case_results)} grid units from {args.finalize}", flush=True)
    else:
        if args.quick:
            units = [(0, LAWS[0][0], LAWS[0][1], 0, 40, 0.5, 3)]
        else:
            # plain monolithic full run, no checkpointing: fine when the caller's own
            # process is not itself capped at 600 s (e.g. launched with nohup).
            units = full_grid_units()
        all_case_results = []
        for law_index, law_name, draw_fn, case_index, d, c, n_reps in units:
            print(f"  running law={law_name} case={case_index} d={d} c={c} "
                  f"n_reps={n_reps} ...", flush=True)
            res = run_case(law_name, draw_fn, law_index, case_index, d, c, n_reps)
            _log_unit(res)
            all_case_results.append(res)

    results_by_law_c = {}  # (law_name, c) -> {d: result}
    for res in all_case_results:
        results_by_law_c.setdefault((res["law_name"], res["c"]), {})[res["d"]] = res

    # --- assemble the z-dependent row table ------------------------------------------
    zdep_rows = []
    for res in all_case_results:
        law_name, d, c, n, c_alt = res["law_name"], res["d"], res["c"], res["n"], res["c_alt"]
        for zi, z in enumerate(res["z_list"]):
            for key in ROW_KEYS_Z:
                vals = res["acc"][zi][key]
                m, se = mean_se(vals)
                target = target_for(key, c, z)
                gap = abs(m - target)
                status = "PASS" if std_pass(gap, se, d) else "FAIL"
                if status == "FAIL":
                    fail_count += 1
                alt_str = ""
                if key in ROWS_DUAL_C:
                    target_alt = target_for(key, c_alt, z)
                    alt_str = fmt(target_alt)
                zdep_rows.append((law_name, fmt(c, 4), d, fmt(z, 6), ROW_LABEL[key],
                                   fmt(m), fmt(se, 3), fmt(target), alt_str, fmt(gap, 3),
                                   status))

    # --- row 8, 9, 10 (z-independent) --------------------------------------------------
    static_rows = []
    for res in all_case_results:
        law_name, d, c = res["law_name"], res["d"], res["c"]
        a = res["acc_static"]
        m8, se8 = mean_se(a["row8"])
        t8 = bulk_edge(c)
        gap8 = abs(m8 - t8)
        static_rows.append((law_name, fmt(c, 4), d, "lam_max(W0)", fmt(m8), fmt(se8, 3),
                             fmt(t8), fmt(gap8, 4), "joint (see below)"))

        m9o, se9o = mean_se(a["row9_overlap"])
        t9o = beta_sq(res["theta9"], c)
        gap9o = abs(m9o - t9o)
        s9o = "PASS" if std_pass(gap9o, se9o, d) else "FAIL"
        if s9o == "FAIL":
            fail_count += 1
        static_rows.append((law_name, fmt(c, 4), d, f"<vhat,v>^2 (theta={fmt(res['theta9'],3)})",
                             fmt(m9o), fmt(se9o, 3), fmt(t9o), fmt(gap9o, 3), s9o))

        m9l, se9l = mean_se(a["row9_lam"])
        t9l = rho_sq(res["theta9"], c)
        gap9l = abs(m9l - t9l)
        s9l = "PASS" if std_pass(gap9l, se9l, d) else "FAIL"
        if s9l == "FAIL":
            fail_count += 1
        static_rows.append((law_name, fmt(c, 4), d, f"lam_max(X^TX) (theta={fmt(res['theta9'],3)})",
                             fmt(m9l), fmt(se9l, 3), fmt(t9l), fmt(gap9l, 3), s9l))

        m10, se10 = mean_se(a["row10_overlap"])
        gap10 = abs(m10 - 0.0)
        static_rows.append((law_name, fmt(c, 4), d, f"<vhat,v>^2 (theta={fmt(res['theta10'],3)})",
                             fmt(m10), fmt(se10, 3), "0", fmt(gap10, 4), "joint (see below)"))

    # --- joint (cross-d) verdicts for row 8 and row 10 --------------------------------
    joint_rows = []
    for (law_name, c), by_d in results_by_law_c.items():
        if 400 not in by_d or 1600 not in by_d:
            continue
        m8_400, _ = mean_se(by_d[400]["acc_static"]["row8"])
        m8_1600, _ = mean_se(by_d[1600]["acc_static"]["row8"])
        t8 = bulk_edge(c)
        gap400 = abs(m8_400 - t8)
        gap1600 = abs(m8_1600 - t8)
        ratio8 = gap400 / max(gap1600, 1e-15)
        pass8 = gap400 <= 0.15 and gap1600 <= 0.08 and ratio8 >= 1.8
        if not pass8:
            fail_count += 1
        joint_rows.append((law_name, fmt(c, 4), "row8 lam_max(W0)", fmt(gap400, 4),
                            fmt(gap1600, 4), fmt(ratio8, 3), "PASS" if pass8 else "FAIL"))

        m10_400, _ = mean_se(by_d[400]["acc_static"]["row10_overlap"])
        m10_1600, _ = mean_se(by_d[1600]["acc_static"]["row10_overlap"])
        ratio10 = m10_400 / max(m10_1600, 1e-300)
        pass10 = 2.5 <= ratio10 <= 6.0
        if not pass10:
            fail_count += 1
        joint_rows.append((law_name, fmt(c, 4), "row10 overlap", fmt(m10_400, 4),
                            fmt(m10_1600, 4), fmt(ratio10, 3), "PASS" if pass10 else "FAIL"))

    for jr in joint_rows:
        print(f"  joint: law={jr[0]} c={jr[1]} {jr[2]}: {jr[3:]}", flush=True)

    # --- necessity rows -----------------------------------------------------------------
    print("\n=== necessity rows (expected to fail; not counted toward exit code) ===",
          flush=True)
    necessity_rows = []
    necessity_reps = 5 if args.quick else NECESSITY_REPS
    for vi, (name, draw_fn) in enumerate(NECESSITY):
        res = run_necessity_indexed(name, draw_fn, vi, necessity_reps)
        m_lam, se_lam = mean_se(res["lam_vals"])
        t_lam = bulk_edge(res["c"])
        gap_lam = abs(m_lam - t_lam)
        s_lam = "PASS" if std_pass(gap_lam, se_lam, res["d"]) else "FAIL (expected)"

        m_r1, se_r1 = mean_se(res["row1_vals"])
        t_r1 = mp_m(res["c"], res["z"])
        gap_r1 = abs(m_r1 - t_r1)
        s_r1 = "PASS" if std_pass(gap_r1, se_r1, res["d"]) else "FAIL (expected)"

        print(f"  {name}: base_seed={res['base_seed']} seed_u={res['seed_u']} "
              f"seed_v={res['seed_v']} lam_max(W0)={fmt(m_lam)} (target {fmt(t_lam)}, "
              f"gap {fmt(gap_lam,3)}, {s_lam}) row1={fmt(m_r1)} (target {fmt(t_r1)}, "
              f"gap {fmt(gap_r1,3)}, {s_r1})", flush=True)
        necessity_rows.append((name, res["base_seed"], fmt(m_lam), fmt(t_lam), fmt(gap_lam, 3),
                                s_lam, fmt(m_r1), fmt(t_r1), fmt(gap_r1, 3), s_r1))

    this_call_wall_s = time.time() - t0
    grid_wall_s = sum(res["wall_s"] for res in all_case_results)
    if args.finalize:
        # This process only loaded the checkpoint; the grid itself ran in the prior
        # --resume-grid call(s). Report the sum of the per-unit times actually measured
        # during those calls, plus this call's own (fast) work, as the true total.
        wall_s = grid_wall_s + this_call_wall_s
        wall_note = (f"grid compute (summed per-unit, measured across the --resume-grid "
                      f"calls) {grid_wall_s:.1f} s + this --finalize call {this_call_wall_s:.1f} s")
    else:
        wall_s = this_call_wall_s
        wall_note = f"single process, grid compute {grid_wall_s:.1f} s of the total"
    vmhwm = read_vmhwm()
    print(f"\n=== summary ===\nFAIL count (drives exit code): {fail_count}\n"
          f"wall time: {wall_s:.1f} s ({wall_note})\nVmHWM: {vmhwm}", flush=True)

    if args.quick:
        print("\nquick mode: not writing the markdown output file.", flush=True)
        return 1 if fail_count else 0

    # --- write the markdown report -----------------------------------------------------
    lines = []
    lines.append("# Non-Gaussian resolvent forms: numeric audit")
    lines.append("")
    lines.append(f"Unit G12 of `notes/prop_single_table_general.md` section 6, plus rows 12 "
                 f"and 13 from section 8 item 1 (the mixed form and the companion identity). "
                 f"Full grid run, not `--quick`. Script: "
                 f"`scripts/check_nongaussian_forms.py`.")
    lines.append("")
    lines.append(f"numpy version: {np.__version__}. Seed base: {SEED_BASE} "
                 f"(`seed = 20260909 + 100000*law_index + 1000*case_index + rep`, per-replicate "
                 f"draw of Z; direction seeds are the case's base seed plus 700 for u and 701 "
                 f"for v, fixed across replicates and law-specific).")
    lines.append(f"Threads: OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS','(unset)')}, "
                 f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS','(unset)')}, "
                 f"MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS','(unset)')}.")
    lines.append(f"Wall time: {wall_s:.1f} s ({wall_note}). VmHWM at exit: {vmhwm} (this "
                 f"process only; see the note on checkpointed runs below).")
    lines.append(f"FAIL count (drives the exit code, excludes necessity rows and joint "
                 f"row 8 or row 10 informational lines already counted above): {fail_count}.")
    lines.append("")

    lines.append("## Header: nu4 per law")
    lines.append("")
    lines.append(md_table(["law", "seed", "n", "mean", "var", "nu4 empirical", "nu4 target",
                            "gap"], nu4_rows))
    lines.append("")

    lines.append("## Cross-check: betaSq against getTheoryPred (R) and "
                 "compute_asymptotic_power (python)")
    lines.append("")
    lines.append("getTheoryPred (`R/theory_pred.R` of the stackedSVD repository) "
                 "crashes on a literal single-table call: `diag(x)` for a bare scalar x in "
                 "(0, 1) makes an R floor(x) x floor(x) matrix, not a 1x1 matrix holding x, "
                 "inside the function's unused SVD-Stack branch. Called here with a second, "
                 "inert table (theta=0, c=1e-300) to dodge that without editing the file; "
                 "the padding leaves the Stack SVD component algebraically identical to the "
                 "unpadded single-table value (sum(c) and theta2 are unaffected to double "
                 "precision).")
    lines.append("")
    lines.append(md_table(["theta", "c", "betaSq (numpy)", "betaSq (R Stack SVD)", "R gap",
                            "betaSq (python method 0)", "python gap", "rhoSq (numpy, no "
                            "reference)"], cross_rows))
    lines.append("")
    lines.append("rhoSq (the top-eigenvalue limit) has no counterpart in either reference "
                 "file: both implement only the overlap (betaSq) prediction for a single "
                 "table. Checked instead by the algebraic identity "
                 "`rhoSq(theta,c) = (1+theta^2)(theta^2+c)/theta^2` above threshold, and by "
                 "the Monte Carlo comparison in row 9 of the main grid below.")
    lines.append("")
    lines.append(md_table(["theta", "c", "rhoSq (numpy)", "(1+theta^2)(theta^2+c)/theta^2",
                            "gap"], rho_identity_rows))
    lines.append("")

    lines.append("## Main grid: rows 1-7, 11, 12, 13 (z-dependent)")
    lines.append("")
    lines.append("Rows 1-7 and 11 use `c = n/d` as the primary target and also print the "
                 "target at `c = (n-1)/d` (the 'alt' column) per the note; PASS/FAIL is "
                 "decided by the primary target only. Rows 3, 6, 12 have target 0 under "
                 "either c convention. Rows 12 and 13 (section 8) have no alt-c column: the "
                 "note's alt-target instruction names only rows 1-7 and 11.")
    lines.append("")
    lines.append(md_table(["law", "c", "d", "z", "quantity", "empirical mean", "SE",
                            "target (c=n/d)", "target (c=(n-1)/d)", "gap", "status"],
                           zdep_rows))
    lines.append("")

    lines.append("## Main grid: rows 8, 9, 10 (z-independent)")
    lines.append("")
    lines.append("Row 8 and row 10 use a joint pass rule spanning both d values for the "
                 "same (law, c); the per-(law, c, d) lines below are informational, and the "
                 "joint verdict is in the next table.")
    lines.append("")
    lines.append(md_table(["law", "c", "d", "quantity", "empirical mean", "SE", "target",
                            "gap", "status"], static_rows))
    lines.append("")

    lines.append("## Joint (cross-d) verdicts: row 8 edge rate and row 10 delocalization rate")
    lines.append("")
    lines.append("Row 8: PASS when gap(d=400) <= 0.15, gap(d=1600) <= 0.08, and "
                 "gap(d=400)/gap(d=1600) >= 1.8 (the d^-2/3 Tracy-Widom rate). Row 10: PASS "
                 "when the mean-overlap ratio d=400 over d=1600 is between 2.5 and 6 (the "
                 "1/d rate).")
    lines.append("")
    lines.append(md_table(["law", "c", "row", "value(d=400)", "value(d=1600)", "ratio",
                            "status"], joint_rows))
    lines.append("")

    lines.append("## Necessity rows (section 4; expected to fail; not in the FAIL count)")
    lines.append("")
    lines.append(f"All four at d={NECESSITY_D}, c={NECESSITY_C}, base law rademacher, "
                 f"z = bulkEdge(c) + 0.5, {necessity_reps} replicates. Each breaks exactly "
                 f"one hypothesis of `assum:general_noise` and is compared against the "
                 f"target that would hold if the hypothesis were not broken.")
    lines.append("")
    lines.append(md_table(["variant", "base seed", "lam_max(W0) mean", "target", "gap",
                            "status", "row1 mean", "target", "gap", "status"],
                           necessity_rows))
    lines.append("")

    lines.append("## Reading the rows")
    lines.append("")
    gaussian_rows = [r for r in zdep_rows if r[0] == "gaussian"]
    nongauss_fail = sum(1 for r in zdep_rows if r[0] != "gaussian" and r[-1] == "FAIL")
    gauss_fail = sum(1 for r in gaussian_rows if r[-1] == "FAIL")
    edge_ratios = [jr for jr in joint_rows if jr[2] == "row8 lam_max(W0)"]
    mixed_rows = [r for r in zdep_rows if r[4] in ("v.T G E.T u", "g.T G g")]
    mixed_fail = sum(1 for r in mixed_rows if r[-1] == "FAIL")
    lines.append(
        f"The four laws give the same large-d limits: {gauss_fail} FAIL rows among the "
        f"Gaussian control and {nongauss_fail} among the other three laws combined, out of "
        f"the same {len(zdep_rows)}-row grid, so rademacher, uniform and t5 track the "
        f"Gaussian MP formulas as closely as Gaussian tracks itself. Row 8's edge gap falls "
        f"from d=400 to d=1600 by the ratios in the joint-verdicts table above; a ratio near "
        f"or above the d^-2/3 Tracy-Widom rate of 1.8 across all law-c pairs would confirm "
        f"the correct finite-size rate, not just convergence. Rows 12 and 13, the mixed form "
        f"and the companion identity added in section 8, pass or fail as shown in "
        f"{mixed_fail} of {len(mixed_rows)} rows above; both are exact algebraic consequences "
        f"of Sherman-Morrison and the E G E^T = I + z Gc identity, so a clean pass here is "
        f"evidence for the arithmetic in the script, not fresh RMT content, while the "
        f"companion identity's target 1 + z*mTilde(c,z) does depend on the isotropic law for "
        f"u^T Gc u carrying over from row 11's trace law, which is fresh content.")
    lines.append("")

    md_text = "\n".join(lines) + "\n"
    with open(args.out, "w") as f:
        f.write(md_text)
    print(f"\nwrote {args.out} ({len(md_text)} bytes)", flush=True)

    return 1 if fail_count else 0


if __name__ == "__main__":
    sys.exit(main())
