"""The 'reverse Weyl' bound (NOT in the formalized toolkit: it needs a lamMin
lower edge, i.e. Gordon's minimax, which we do not have). Measured only to size
the prize.  lamMax(N) <= lamMax(N+B) - lamMin(B) with B = (1/d) sum (t-w_i^2) E_i^T E_i
and t = max_i w_i^2, so N+B = t (1/d) E^T E and lamMax(N+B) -> t (1+sqrt(sum c))^2."""
import numpy as np, json
from scipy.optimize import brentq
from het import bHet, phi, zfun, bStair, bSF

def aHet(c, w2):
    """lower edge of the generalized MP law: zfun at the critical point s > 0.
    0 when sum c_i <= 1 (hard edge at 0) or every weight is 0."""
    if np.max(w2) <= 0 or np.sum(c) <= 1.0:
        return 0.0
    hi = 1.0
    while phi(c, w2, hi) > 0:
        hi *= 2.0
        if hi > 1e12:
            return 0.0
    s = brentq(lambda x: phi(c, w2, x), 1e-14, hi, xtol=1e-16, rtol=8.9e-16)
    return float(zfun(c, w2, s))

def bRev(c, w2):
    t = float(np.max(w2))
    return t*(1+np.sqrt(np.sum(c)))**2 - aHet(c, t - w2)

R = json.load(open("bigsweep.json"))
