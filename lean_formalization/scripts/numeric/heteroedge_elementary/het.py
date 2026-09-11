"""Scalars transcribed from StackedSVD/RMT/Het/MPhet.lean and StackSVDWeighted.lean.

zfun(s)      = -1/s + sum_i c_i w_i^2 / (1 + w_i^2 s)
phi(s)       = 1 - sum_i c_i w_i^4 s^2 / (1 + w_i^2 s)^2      (= s^2 zfun'(s))
sLo          = -1 / max_i w_i^2
sStar        = unique zero of phi in (sLo, 0)
bHet         = zfun(sStar)
gammaTop     = unique root g > max_i w_i^2 of 1 + sum_j th_j^2 w_j^2/(w_j^2 - g)
rhoHet       = gammaTop * (1 + sum_i c_i w_i^2/(gammaTop - w_i^2))
bSF          = (max_i |w_i| + sqrt(sum_i c_i w_i^2))^2
optWstack_i  = th_i / sqrt(th_i^2 + c_i)
"""
import numpy as np
from scipy.optimize import brentq

def zfun(c, w2, s):
    return -1.0/s + np.sum(c*w2/(1.0 + w2*s))

def phi(c, w2, s):
    return 1.0 - np.sum(c*w2**2*s**2/(1.0 + w2*s)**2)

def sStar(c, w2):
    w2m = np.max(w2)
    lo = -1.0/w2m
    # phi -> 1 as s -> 0^-, phi -> -inf as s -> lo^+ (block(s) with max weight blows up)
    a, b = lo*(1-1e-14), -1e-14
    # bracket: walk in from the left
    fa = phi(c, w2, a)
    assert fa < 0, ("phi(lo+) not negative", fa)
    fb = phi(c, w2, b)
    assert fb > 0, ("phi(0-) not positive", fb)
    return brentq(lambda x: phi(c, w2, x), a, b, xtol=1e-17, rtol=8.9e-16, maxiter=300)

def bHet(c, w2):
    return zfun(c, w2, sStar(c, w2))

def bSF(c, w2):
    return (np.sqrt(np.max(w2)) + np.sqrt(np.sum(c*w2)))**2

def secular(th2, w2, g):
    return 1.0 + np.sum(th2*w2/(w2 - g))

def gammaTop(th2, w2):
    w2m = np.max(w2)
    a = w2m*(1+1e-13) + 1e-15
    fa = secular(th2, w2, a)
    if not (fa < 0):
        return None                     # no outlier root
    b = a + 1.0
    while secular(th2, w2, b) < 0:
        b *= 2.0
        if b > 1e12:
            return None
    return brentq(lambda x: secular(th2, w2, x), a, b, xtol=1e-15, rtol=8.9e-16, maxiter=300)

def rhoHet(th2, c, w2):
    g = gammaTop(th2, w2)
    if g is None:
        return None
    return g*(1.0 + np.sum(c*w2/(g - w2)))

def optW2(th2, c):
    """optWstack squared: w_i^2 = th_i^2/(th_i^2 + c_i)."""
    return th2/(th2 + c)

# ---- the elementary bounds ----------------------------------------------
def fSF(c, a):
    """f(a) = (max_i sqrt(a_i) + sqrt(sum_i c_i a_i))^2, the SF edge at weights sqrt(a)."""
    if np.all(a <= 0):
        return 0.0
    return (np.sqrt(np.max(a)) + np.sqrt(np.sum(c*a)))**2

def bPart(c, w2):
    """min over set partitions of sum_g f(a restricted to g). M <= 8."""
    M = len(w2)
    best = np.inf
    bestP = None
    for P in set_partitions(M):
        tot = 0.0
        for g in P:
            a = np.zeros(M); a[list(g)] = w2[list(g)]
            tot += fSF(c, a)
        if tot < best:
            best, bestP = tot, P
    return best, bestP

def set_partitions(M):
    """all set partitions of {0..M-1} (restricted growth strings)."""
    if M == 0:
        yield []
        return
    a = [0]*M; m = [0]*M
    while True:
        P = {}
        for i, b in enumerate(a):
            P.setdefault(b, []).append(i)
        yield [tuple(v) for v in P.values()]
        i = M - 1
        while i > 0 and a[i] > m[i-1]:
            i -= 1
        if i == 0:
            return
        a[i] += 1
        mi = max(m[i-1], a[i])
        for j in range(i+1, M):
            a[j] = 0; m[j] = mi
        m[i] = mi

# ---- the staircase (layer-cake) bound: closed form ----------------------
def bStair(c, w2):
    """B = sum_k (w2_(k) - w2_(k+1)) * (1 + sqrt(c_(1)+...+c_(k)))^2, indices
    sorted by decreasing w2_i, w2_(M+1) = 0.  Greedy optimum of the polymatroid
    {y >= 0 : y(T) <= (1+sqrt(c(T)))^2 for every nonempty T}."""
    o = np.argsort(-np.asarray(w2, float))
    ws = np.asarray(w2, float)[o]; cs = np.asarray(c, float)[o]
    C = np.cumsum(cs)
    lev = np.append(ws, 0.0)
    return float(np.sum((lev[:-1] - lev[1:])*(1.0 + np.sqrt(C))**2))

def stairSplit(c, w2):
    """the M-piece split attaining bStair: row k is the top-k indices at height
    w2_(k) - w2_(k+1)."""
    M = len(w2)
    o = np.argsort(-np.asarray(w2, float))
    ws = np.asarray(w2, float)[o]; lev = np.append(ws, 0.0)
    A = np.zeros((M, M))
    for k in range(M):
        A[k, o[:k+1]] = lev[k] - lev[k+1]
    return A
