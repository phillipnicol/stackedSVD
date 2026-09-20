#!/usr/bin/env python3
"""Xcount2: is the sorried one-step lemma true?  Exact, no randomness, no seed.

The lemma left as a sorry in Count2.lean is

    clsCard n d k v * min(n,d)  <=  (2k)^6 * clsCard n d k (v+1),     2 <= v <= k,
    under  1 <= k  and  2k <= min(n,d).

clsCard is rebuilt from the shape table as
    clsCard(n,d,k,v) = sum_{s+t=v} S(k,s,t) * n^(s) * d^(t)     (falling factorials),
which is first cross-checked against a direct enumeration of every walk.

Scanned: k = 2..7; min(n,d) from 2k (the smallest the hypothesis allows) up to
10^9, both n >= d and d >= n, and a few extreme ratios.  Printed: the number of
violations and the worst ratio  lhs / rhs.

walk_mult is inlined from rules.py and rgs from shapes_odd.py (both scratch
modules of the Xcount2 stage-3 session, 2026-09-10) so this script is standard
library only and needs no sibling file.
"""
import itertools
from collections import Counter
from fractions import Fraction


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


def shape_table(k):
    S = Counter()
    for i in rgs(k):
        for j in rgs(k):
            cnt = walk_mult(i, j, k)
            if any(c == 1 for c in cnt.values()):
                continue
            S[(len(set(i)), len(set(j)))] += 1
    return S


def falling(n, s):
    r = 1
    for a in range(s):
        r *= (n - a)
        if r <= 0:
            return 0
    return r


def cls_from_shapes(S, n, d, v):
    tot = 0
    for (s, t), c in S.items():
        if s + t == v:
            tot += c * falling(n, s) * falling(d, t)
    return tot


def cls_direct(n, d, k, v):
    tot = 0
    for i in itertools.product(range(n), repeat=k):
        for j in itertools.product(range(d), repeat=k):
            cnt = walk_mult(i, j, k)
            if any(c == 1 for c in cnt.values()):
                continue
            if len(set(i)) + len(set(j)) == v:
                tot += 1
    return tot


def main():
    print("=== Xcount2: the sorried one-step lemma, exact, no seed ===")
    # 1. cross-check the shape reconstruction against a direct enumeration
    print("-- cross-check clsCard from shapes against direct enumeration --")
    bad = 0
    for k in (2, 3, 4):
        S = shape_table(k)
        for n in (1, 2, 3, 4):
            for d in (1, 2, 3, 4):
                if (n ** k) * (d ** k) > 400000:
                    continue
                for v in range(0, k + 3):
                    a = cls_from_shapes(S, n, d, v)
                    b = cls_direct(n, d, k, v)
                    if a != b:
                        bad += 1
                        print("   MISMATCH k=%d n=%d d=%d v=%d shapes=%d direct=%d"
                              % (k, n, d, v, a, b))
    print("   mismatches:", bad)

    # 2. the one-step lemma itself
    print("-- one-step lemma over the admissible sizes --")
    sizes = []
    for m in (1, 2, 4, 16, 1000, 10 ** 6, 10 ** 9):
        for ratio in (1, 2, 10, 1000, 10 ** 6):
            sizes.append((m, ratio))
    for k in range(2, 8):
        S = shape_table(k)
        viol = 0
        worst = Fraction(0)
        worst_at = None
        undef = 0
        for (mm, ratio) in sizes:
            m = 2 * k * mm                  # min(n,d) >= 2k, the hypothesis
            for (n, d) in ((m, m * ratio), (m * ratio, m)):
                for v in range(2, k + 1):
                    lhs = cls_from_shapes(S, n, d, v) * min(n, d)
                    rhs = (2 * k) ** 6 * cls_from_shapes(S, n, d, v + 1)
                    if lhs > rhs:
                        viol += 1
                        if viol <= 3:
                            print("   VIOLATION k=%d n=%d d=%d v=%d lhs=%d rhs=%d"
                                  % (k, n, d, v, lhs, rhs))
                    if rhs == 0:
                        if lhs > 0:
                            undef += 1
                        continue
                    r = Fraction(lhs, rhs)
                    if r > worst:
                        worst = r
                        worst_at = (n, d, v)
        print("k=%d  violations: %d   empty-right-side with nonzero left: %d"
              "   worst lhs/rhs = %.4f at (n,d,v)=%s"
              % (k, viol, undef, float(worst), worst_at))

    # 3. the same at the smallest admissible size, min(n,d) = 2k exactly
    print("-- smallest admissible size, min(n,d) = 2k --")
    for k in range(2, 8):
        S = shape_table(k)
        worst = Fraction(0)
        worst_at = None
        viol = 0
        for extra in range(0, 6):
            n = 2 * k + extra
            d = 2 * k
            for (a, b) in ((n, d), (d, n)):
                for v in range(2, k + 1):
                    lhs = cls_from_shapes(S, a, b, v) * min(a, b)
                    rhs = (2 * k) ** 6 * cls_from_shapes(S, a, b, v + 1)
                    if lhs > rhs:
                        viol += 1
                    if rhs > 0:
                        r = Fraction(lhs, rhs)
                        if r > worst:
                            worst = r
                            worst_at = (a, b, v)
        print("k=%d  violations: %d  worst lhs/rhs = %.4f at (n,d,v)=%s"
              % (k, viol, float(worst), worst_at))


if __name__ == '__main__':
    main()
