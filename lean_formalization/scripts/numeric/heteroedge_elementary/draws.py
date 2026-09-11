"""Reused protocol (notes/plan_heterolaw_A.md section 3.6, seed 20260831):
M in {2,3,4}, theta_i ~ U(0.3,2), c_i ~ U(0.3,3), w = optWstack, keep the
supercritical draws (sum theta_i^4/c_i > 1)."""
import numpy as np
from het import optW2

SEED = 20260831

def draws(ntry=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(ntry):
        M = int(rng.integers(2, 5))
        th = rng.uniform(0.3, 2.0, M)
        c = rng.uniform(0.3, 3.0, M)
        if np.sum(th**4/c) > 1.0:
            out.append((th, c, optW2(th**2, c)))
    return out
