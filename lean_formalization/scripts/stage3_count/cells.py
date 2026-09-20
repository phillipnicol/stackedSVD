"""Unit Xcount: numeric check of the cell-refined form of the Furedi-Komlos count.

No randomness, so no seed.  Exact integer arithmetic throughout.

Shapes: a walk (i, j) : [k] -> [n], [k] -> [d] is put in canonical form by relabeling the
rows in order of first appearance and the columns in order of first appearance.  A canonical
shape with s rows and t columns has falling(n, s) * falling(d, t) labeled copies, so

    cellCard(n, d, k, s, t) = Sh(k, s, t) * falling(n, s) * falling(d, t).

Claim L2 (all D extra vertices on the row side), for s + t = v <= k, D = k + 1 - v:

    cellCard(s, t) * min(n, d)^D  <=  ((2k)^6)^D * cellCard(s + D, t).

Claim L1: clsCard(v) = sum over s of cellCard(s, v - s).  (Trivially true by construction.)
Claim TR: Tr(k, s, t) >= 1 for every s, t >= 1 with s + t = k + 1.
"""
from itertools import product
from math import comb

def rgs(k):
    """All restricted growth strings of length k (canonical set partitions of [k])."""
    if k == 0:
        return [()]
    out = []
    def rec(pref, mx):
        if len(pref) == k:
            out.append(tuple(pref)); return
        for a in range(mx + 2):
            rec(pref + [a], max(mx, a))
    rec([0], 0)
    return out

def falling(n, s):
    p = 1
    for u in range(s):
        p *= (n - u)
    return p if p > 0 else 0

def shapes(k):
    """Sh[(s,t)] and Tr[(s,t)]: canonical NoSingle shapes, and the tree ones (v = k+1)."""
    Sh, Tr = {}, {}
    R = list(rgs(k))
    for i in R:
        s = max(i) + 1
        for j in R:
            t = max(j) + 1
            mult = {}
            for u in range(k):
                for e in ((i[u], j[u]), (i[u], j[(u + 1) % k])):
                    mult[e] = mult.get(e, 0) + 1
            if any(m == 1 for m in mult.values()):
                continue
            Sh[(s, t)] = Sh.get((s, t), 0) + 1
            if s + t == k + 1:
                Tr[(s, t)] = Tr.get((s, t), 0) + 1
    return Sh, Tr

def main():
    bad_l2, bad_tr, cells = 0, 0, 0
    for k in range(1, 8):
        Sh, Tr = shapes(k)
        # TR: every cell of the tree class is nonempty
        for s in range(1, k + 1):
            t = k + 1 - s
            if Tr.get((s, t), 0) < 1:
                print(f"  TR FAIL k={k} s={s} t={t}")
                bad_tr += 1
        for (n, d) in [(2 * k, 2 * k), (2 * k, 2 * k + 1), (2 * k + 1, 2 * k), (10 ** 3, 2 * k),
                       (2 * k, 10 ** 3), (10 ** 6, 10 ** 3), (10 ** 3, 10 ** 6),
                       (10 ** 12, 10 ** 9), (10 ** 9, 10 ** 12), (3 * k, 7 * k)]:
            m = min(n, d)
            if 2 * k > m:
                continue
            F = (2 * k) ** 6
            for (s, t), c in sorted(Sh.items()):
                v = s + t
                if v > k:            # D >= 1 only; D = 0 is the trivial case
                    continue
                D = k + 1 - v
                lhs = c * falling(n, s) * falling(d, t) * m ** D
                rhs = (F ** D) * Tr.get((s + D, t), 0) * falling(n, s + D) * falling(d, t)
                cells += 1
                if lhs > rhs:
                    bad_l2 += 1
                    if bad_l2 <= 8:
                        print(f"  L2 FAIL k={k} n={n} d={d} s={s} t={t} D={D} "
                              f"lhs/rhs={lhs / rhs:.4g}")
        # worst ratio at this k over the tested (n, d)
        worst = 0.0
        for (n, d) in [(2 * k, 2 * k), (10 ** 12, 10 ** 9), (10 ** 3, 10 ** 6)]:
            m = min(n, d); F = (2 * k) ** 6
            for (s, t), c in sorted(Sh.items()):
                v = s + t
                if v > k:
                    continue
                D = k + 1 - v
                lhs = c * falling(n, s) * falling(d, t) * m ** D
                rhs = (F ** D) * Tr.get((s + D, t), 0) * falling(n, s + D) * falling(d, t)
                if rhs:
                    worst = max(worst, lhs / rhs)
        print(f"k={k}: shapes cells={len(Sh)} tree cells={len(Tr)} worst L2 ratio={worst:.4g}")
    print(f"TOTAL cells tested={cells}  L2 violations={bad_l2}  TR violations={bad_tr}")

main()
