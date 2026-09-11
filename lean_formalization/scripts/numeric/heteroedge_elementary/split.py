"""bSplit(K): minimize sum_k f(a_k) over a_{k,i} >= 0 with sum_k a_{k,i} = w2_i,
plus a cutting-plane bracket for the convex envelope B(w2) = lim_K bSplit(K)."""
import numpy as np
from scipy.optimize import minimize, linprog
from het import fSF, bPart, set_partitions

EPS = 1e-300

def fvec(c, B):
    """f on each row of B (n x M), vectorized."""
    return (np.sqrt(np.maximum(B.max(axis=1), 0.0)) + np.sqrt(np.maximum(B @ c, 0.0)))**2

def total(c, A):
    return float(fvec(c, A).sum())

def _softmax(L):
    L = L - L.max(axis=0, keepdims=True)
    P = np.exp(L)
    return P / P.sum(axis=0, keepdims=True)

def _objgrad(x, c, w2, K, M):
    P = _softmax(x.reshape(K, M))
    A = P * w2[None, :]
    u = A.max(axis=1); v = A @ c
    su = np.sqrt(np.maximum(u, 0.0)); sv = np.sqrt(np.maximum(v, 0.0))
    val = float(((su + sv)**2).sum())
    ru = np.where(su > 1e-150, sv/np.maximum(su, 1e-150), 0.0)   # sqrt(v/u)
    rv = np.where(sv > 1e-150, su/np.maximum(sv, 1e-150), 0.0)   # sqrt(u/v)
    G = (1.0 + rv)[:, None] * c[None, :]
    am = A.argmax(axis=1)
    G[np.arange(K), am] += 1.0 + ru
    gL = w2[None, :] * P * (G - (P*G).sum(axis=0, keepdims=True))
    return val, gL.ravel()

def bSplit(c, w2, K, seed=0, nrestart=25, maxiter=400):
    M = len(w2)
    rng = np.random.default_rng(seed)
    best, bestA = np.inf, None
    starts = []
    for P in set_partitions(M):
        if len(P) > K:
            continue
        L = np.full((K, M), -5.0)
        for k, g in enumerate(P):
            for i in g:
                L[k, i] = 5.0
        starts.append(L.ravel())
    for _ in range(nrestart):
        starts.append(rng.normal(0, 2.0, K*M))
    for x0 in starts:
        r = minimize(_objgrad, x0, args=(c, w2, K, M), jac=True, method='L-BFGS-B',
                     options=dict(maxiter=maxiter, ftol=1e-13, gtol=1e-11))
        if r.fun < best:
            best = float(r.fun); bestA = _softmax(r.x.reshape(K, M))*w2[None, :]
    return best, bestA

def envelope(c, w2, iters=80, ngrid=6000, seed=0):
    """SUPERSEDED, do not trust: cutting planes against a random grid. The grid
    can miss the violated direction, so the value can sit BELOW the true envelope.
    Use het.bStair, which is exact (see the report, section 4). Kept only because
    it was the route that first suggested the closed form."""
    M = len(w2)
    rng = np.random.default_rng(seed)
    grid = np.vstack([np.eye(M), rng.dirichlet(np.full(M, 0.4), ngrid),
                      rng.dirichlet(np.ones(M), ngrid)])
    fg = fvec(c, grid)
    B = [np.eye(M)[i] for i in range(M)]
    y = None
    for _ in range(iters):
        A_ub = np.array(B); b_ub = fvec(c, A_ub)
        r = linprog(-w2, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)]*M, method='highs')
        if not r.success:
            return -np.inf, np.inf
        y = r.x
        viol = grid @ y - fg
        j = int(np.argmax(viol))
        if viol[j] <= 1e-11:
            break
        B.append(grid[j])
    up = float(y @ w2)
    ratio = max(1.0, float(np.max((grid @ y)/fg)))
    return up/ratio, up

def _merge_starts(c, w2, K):
    """splits with at most K pieces obtained by merging the M staircase layers."""
    from het import stairSplit, set_partitions
    A = stairSplit(c, w2); M = len(w2)
    out = []
    for P in set_partitions(M):
        if len(P) > K:
            continue
        B = np.zeros((K, M))
        for k, g in enumerate(P):
            B[k] = A[list(g)].sum(axis=0)
        out.append(B)
    return out

def bSplitK(c, w2, K, seed=0, nrestart=25, maxiter=500):
    """best split into at most K pieces. Starts: index partitions, staircase-layer
    merges, random. K >= M returns the exact staircase optimum."""
    from het import bStair, stairSplit
    M = len(w2)
    if K >= M:
        return bStair(c, w2), stairSplit(c, w2)
    best, bestA = bSplit(c, w2, K, seed=seed, nrestart=nrestart, maxiter=maxiter)
    for B in _merge_starts(c, w2, K):
        v = total(c, B)
        if v < best:
            best, bestA = v, B
        L = np.log(np.maximum(B/np.maximum(w2[None, :], 1e-300), 1e-8))
        r = minimize(_objgrad, L.ravel(), args=(c, w2, K, M), jac=True,
                     method='L-BFGS-B', options=dict(maxiter=maxiter, ftol=1e-13))
        if r.fun < best:
            best = float(r.fun); bestA = _softmax(r.x.reshape(K, M))*w2[None, :]
    return best, bestA
