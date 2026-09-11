#!/usr/bin/env python3
"""Xcount2: the cell-wise shape bound that the step encoding has to deliver.

Exact enumeration, no randomness, so there is no seed.

Write S(k, s, t) for the number of NoSingle shapes of a closed walk of length
2k with s rows and t columns, and N(k, s) for the number of tree shapes with s
rows (s + t = k + 1; this is the Narayana number C(k,s) C(k,s-1) / k).  Put
D = k + 1 - s - t.

The route of this unit needs, for every cell with D >= 1,

    (1)   S(k, s, t)  <=  ((2k)^6 / 2)^D  *  N(k, s)              [cell wise]

because the labeled counts then follow from
    n^(s) d^(t) * min(n,d)  <=  2 * n^(s) d^(t+1)      (2k <= min(n,d)),
which moves one unit of deficiency from the count to the column labels and
keeps the row count s fixed.

Also printed, to show that the cell grading cannot be dropped:
    (2)   the total shape count sum_{s+t=v} S(k,s,t) against max_s N(k,s),
which is what a cell-blind bound would need; it fails.

This file is renamed from cells.py of the Xcount2 stage-3 session (2026-09-10)
to avoid a name clash with the unrelated cells.py of Xcount, copied alongside it
in this folder. walk_mult is inlined from rules.py and rgs from shapes_odd.py
(the same session); vertex_seq was imported there but never called, so it is
dropped here. Standard library only, no sibling imports.
"""
from collections import Counter
from fractions import Fraction
from math import comb


def walk_mult(i, j, k):
    cnt = Counter()
    for t in range(k):
        cnt[(i[t], j[t])] += 1
        cnt[(i[t], j[(t + 1) % k])] += 1
    return cnt


def rgs(k):
    """Every restricted growth string of length k (canonical labelings)."""
    out = []

    def rec(pref, mx):
        if len(pref) == k:
            out.append(tuple(pref))
            return
        for a in range(mx + 2):
            pref.append(a)
            rec(pref, max(mx, a))
            pref.pop()

    if k == 0:
        return [()]
    rec([], -1)
    return out


def narayana(k, s):
    if s < 1 or s > k:
        return 0
    return comb(k, s) * comb(k, s - 1) // k


def main():
    print("=== Xcount2: cell-wise shape bound (exact, no seed) ===")
    for k in range(2, 8):
        strings = rgs(k)
        S = Counter()
        for i in strings:
            for j in strings:
                cnt = walk_mult(i, j, k)
                if any(c == 1 for c in cnt.values()):
                    continue
                S[(len(set(i)), len(set(j)))] += 1
        # validation: the D = 0 row must be the Narayana numbers
        ok_tree = all(S[(s, k + 1 - s)] == narayana(k, s) for s in range(1, k + 1))
        base = Fraction((2 * k) ** 6, 2)
        worst = Fraction(0)
        worst_cell = None
        viol = 0
        for (s, t), cnt in sorted(S.items()):
            D = k + 1 - s - t
            if D <= 0:
                continue
            rhs = base ** D * narayana(k, s)
            if rhs == 0:
                viol += 1
                print("   EMPTY TREE CELL k=%d s=%d t=%d D=%d S=%d" % (k, s, t, D, cnt))
                continue
            r = Fraction(cnt, 1) / rhs
            if r > worst:
                worst = r
                worst_cell = (s, t, D, cnt, narayana(k, s))
            if r > 1:
                viol += 1
                print("   VIOLATION k=%d s=%d t=%d D=%d S=%d N=%d"
                      % (k, s, t, D, cnt, narayana(k, s)))
        # the cell-blind version, for comparison
        tot = Counter()
        for (s, t), cnt in S.items():
            tot[s + t] += cnt
        maxN = max(narayana(k, s) for s in range(1, k + 1))
        blind_worst = Fraction(0)
        for v, cnt in tot.items():
            D = k + 1 - v
            if D <= 0:
                continue
            r = Fraction(cnt, 1) / (base ** D * maxN)
            blind_worst = max(blind_worst, r)
        # a cell-blind bound must beat N(k,1) = 1, not max_s N(k,s)
        blind1 = Fraction(0)
        for v, cnt in tot.items():
            D = k + 1 - v
            if D <= 0:
                continue
            blind1 = max(blind1, Fraction(cnt, 1) / (base ** D * narayana(k, 1)))
        print("k=%d  tree row = Narayana: %s   cells with D>=1: %d   violations: %d"
              % (k, ok_tree, sum(1 for c in S if k + 1 - c[0] - c[1] >= 1), viol))
        print("    worst cell ratio S / ((2k)^6/2)^D N(k,s) = %.3e   at (s,t,D)=%s S=%d N=%d"
              % (float(worst), worst_cell[:3] if worst_cell else None,
                 worst_cell[3] if worst_cell else 0, worst_cell[4] if worst_cell else 0))
        print("    cell-blind against N(k,1)=1: worst ratio %.3e  (must be <= 1)"
              % float(blind1))


if __name__ == '__main__':
    main()
