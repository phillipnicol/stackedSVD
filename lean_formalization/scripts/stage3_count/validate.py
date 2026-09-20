"""Validation of the shape enumeration of cells.py, and the smallest exponent p that makes
the cell-refined claim L2 true.  Exact integers; no randomness, so no seed.

Run as `python3 scripts/stage3_count/validate.py` from the repository root; it imports
`cells.py` from its own folder, not the current working directory."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cells import shapes, falling
from math import comb

# 1. tree shapes by row count are the Narayana numbers; the totals are the Catalan numbers
print("k | tree shapes by s | Narayana | total | Catalan | all NoSingle shapes")
for k in range(1, 8):
    Sh, Tr = shapes(k)
    tr = [Tr.get((s, k + 1 - s), 0) for s in range(1, k + 1)]
    nar = [comb(k, s) * comb(k, s - 1) // k for s in range(1, k + 1)]
    cat = comb(2 * k, k) // (k + 1)
    print(f"{k} | {tr} | {nar} | {sum(tr)} | {cat} | {sum(Sh.values())}"
          f" | match={tr == nar and sum(tr) == cat}")

# 2. smallest exponent p with  cellCard(s,t) m^D <= ((2k)^p)^D cellCard(s+D,t)
print()
print("smallest exponent p that survives every tested cell:")
for k in range(1, 8):
    Sh, Tr = shapes(k)
    worst = 0.0
    for (n, d) in [(2 * k, 2 * k), (2 * k, 2 * k + 1), (10 ** 3, 2 * k), (2 * k, 10 ** 3),
                   (10 ** 12, 10 ** 9), (10 ** 9, 10 ** 12), (3 * k, 7 * k), (10 ** 6, 10 ** 3)]:
        m = min(n, d)
        if 2 * k > m:
            continue
        for (s, t), c in sorted(Sh.items()):
            v = s + t
            if v > k:
                continue
            D = k + 1 - v
            lhs = c * falling(n, s) * falling(d, t) * m ** D
            rhs = Tr.get((s + D, t), 0) * falling(n, s + D) * falling(d, t)
            if rhs == 0:
                print(f"  EMPTY TARGET k={k} s={s} t={t} D={D}")
                continue
            worst = max(worst, (lhs / rhs) ** (1.0 / D) / (2 * k))
    # worst is the needed (2k)^p per unit of D, as a power of 2k
    import math
    p = math.log(worst) / math.log(2 * k) + 1 if worst > 0 else 0
    print(f"k={k}: needed factor per unit D = (2k)^{p:.3f}")
