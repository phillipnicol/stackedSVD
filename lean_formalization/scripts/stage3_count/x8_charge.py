#!/usr/bin/env python3
"""X8: the charging decomposition behind #bad <= 3 D.  Exact, no randomness, no seed.

Write E for the distinct entries of the walk, e = |E|, v = s + t, and
    X = k - e = (1/2) sum_e (mult(e) - 2)      (excess multiplicity units)
    c = e - v + 1                              (cycle rank of the edge graph)
The identity D = X + c is checked at every shape.

Then the largest #bad is tabulated against (X, c), and the three case counts of the
charging argument are printed:
    case1 : bad step whose source has odd degree at least 3 in the odd set
    case2 : bad step whose source has odd degree exactly 1 (the walk leaves the open
            edge behind)
    case3 : bad step at a source with odd degree 0 (the walk is at the root, closed)
"""
import sys, time, os
from collections import Counter
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("x8_enum", os.path.join(HERE, "x8_enum.py"))
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)


def analyze(useq, k):
    n2 = 2 * k
    seq = list(useq) + [useq[0]]
    mult = Counter()
    for s in range(n2):
        x, y = seq[s], seq[s + 1]
        mult[(y[1], x[1]) if s % 2 == 0 else (x[1], y[1])] += 1
    e = len(mult)
    visited = {seq[0]}
    odd = {}
    nbad = 0
    cases = [0, 0, 0]
    for s in range(n2):
        x, y = seq[s], seq[s + 1]
        ox = odd.get(x, set())
        if y not in visited:
            visited.add(y)
        else:
            forced = False
            if len(ox) == 1:
                f = next(iter(ox))
                other = ('r', f[0]) if x[0] == 'c' else ('c', f[1])
                forced = (other == y)
            if not forced:
                nbad += 1
                if len(ox) >= 3:
                    cases[0] += 1
                elif len(ox) == 1:
                    cases[1] += 1
                else:
                    cases[2] += 1
        f = (y[1], x[1]) if s % 2 == 0 else (x[1], y[1])
        for z in (x, y):
            sz = odd.setdefault(z, set())
            if f in sz:
                sz.discard(f)
            else:
                sz.add(f)
    return e, nbad, cases


def main(kmax):
    print("=== X8: the charging decomposition (exact, no seed) ===")
    for k in range(2, kmax + 1):
        t0 = time.time()
        best = {}
        idok = True
        c3_at_root = True
        worst_lin = (0, None)
        for (useq, nr, nc) in E.enumerate_shapes(k):
            v = nr + nc
            D = k + 1 - v
            if D < 1:
                continue
            e, nbad, cases = analyze(useq, k)
            X = k - e
            c = e - v + 1
            if X + c != D:
                idok = False
            key = (X, c)
            if nbad > best.get(key, (0, None))[0]:
                best[key] = (nbad, useq, cases)
            # test the linear form #bad <= 2 X + 3 c
            slack = 2 * X + 3 * c - nbad
            if worst_lin[1] is None or slack < worst_lin[0]:
                worst_lin = (slack, (useq, nbad, X, c))
        print("k=%d  D = X + c at every shape: %s   [%.1f s]"
              % (k, idok, time.time() - t0))
        print("    max #bad by (X, c): %s"
              % sorted((key, val[0]) for key, val in best.items()))
        print("    worst slack of #bad <= 2X + 3c : %d  at %s (#bad=%d X=%d c=%d)"
              % (worst_lin[0], E.fmt(worst_lin[1][0]), worst_lin[1][1],
                 worst_lin[1][2], worst_lin[1][3]))
        cc = Counter()
        for key, val in best.items():
            cc[key] = val[2]
        print("    cases (deg>=3, deg=1, deg=0) at the worst shape of each (X,c): %s"
              % sorted(cc.items()))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
