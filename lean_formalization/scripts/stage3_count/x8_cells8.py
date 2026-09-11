#!/usr/bin/env python3
"""X8 numerics, part 2: the exponent p the absolute count needs, cell by cell.

Exact integer and Fraction arithmetic over the complete shape enumeration of
x8_enum.py, no randomness, so there is no seed.

For a cell (s, t) with D = k + 1 - s - t >= 1 the absolute route pays

  tight : sum_b T_b (2k)^b / N(k, s+D)          T_b = number of DISTINCT type
          words of shapes of the cell with b bad steps (what the encoding
          really costs if the type words were counted exactly)
  proof : C(k,s) C(k,t-1) sum_{b<=3D} C(2k,b) (k+1)^b / N(k, s+D)
          (the bound Lean can prove: the innovative positions split by parity,
          the bad positions a subset, the endpoint one of the k+1 visited
          vertices)
  crude : C(k,s) C(k,t-1) (3D+1) (2k)^(6D) / N(k, s+D)   (the closed form)

The inequality the proof needs, after the falling factorial step
falling(n,s) min(n,d)^D <= 2^D falling(n,s+D), is  ratio * 2^D <= (2k)^(p D).
Printed per k: the largest ratio over the cells and the smallest integer p that
works at every cell, for each accounting.  Also Q = k C(k,s) / C(k,s+D-1), the
binomial-to-Narayana ratio, against the claim Q <= k^D.
"""
import sys, time, os
from fractions import Fraction
from math import comb
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("x8_enum", os.path.join(HERE, "x8_enum.py"))
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)


def narayana(k, s):
    if s < 1 or s > k:
        return 0
    return comb(k, s) * comb(k, s - 1) // k


def smallest_p(ratio, D, k):
    """least integer p >= 0 with ratio * 2^D <= (2k)^(p D)."""
    lhs = ratio * 2 ** D
    p = 0
    while Fraction((2 * k) ** (p * D)) < lhs:
        p += 1
        if p > 40:
            return None
    return p


def main(kmax):
    print("=== X8: the exponent p per cell (exact, no seed) ===")
    for k in range(2, kmax + 1):
        t0 = time.time()
        shapes = E.enumerate_shapes(k)
        cells = {}
        for (useq, nr, nc) in shapes:
            D = k + 1 - nr - nc
            if D < 1:
                continue
            typ, nb, nd = E.classify(useq, k)
            cells.setdefault((nr, nc), {}).setdefault(nb, set()).add(typ)
        pmax = {'tight': 0, 'proof': 0, 'crude': 0}
        rmax = {'tight': Fraction(0), 'proof': Fraction(0), 'crude': Fraction(0)}
        argm = {'tight': None, 'proof': None, 'crude': None}
        qmax = Fraction(0)
        qok = True
        for (s, t), byb in sorted(cells.items()):
            D = k + 1 - s - t
            N = narayana(k, s + D)
            assert N > 0, (k, s, t, D)
            r = {}
            r['tight'] = Fraction(sum(len(v) * (2 * k) ** b for b, v in byb.items()), N)
            r['proof'] = Fraction(comb(k, s) * comb(k, t - 1)
                                  * sum(comb(2 * k, b) * (k + 1) ** b
                                        for b in range(0, 3 * D + 1)), N)
            r['crude'] = Fraction(comb(k, s) * comb(k, t - 1) * (3 * D + 1)
                                  * (2 * k) ** (6 * D), N)
            Q = Fraction(k * comb(k, s), comb(k, s + D - 1))
            qmax = max(qmax, Q / k ** D)
            if Q > k ** D:
                qok = False
            for name in r:
                p = smallest_p(r[name], D, k)
                if p > pmax[name]:
                    pmax[name] = p
                if r[name] > rmax[name]:
                    rmax[name] = r[name]
                    argm[name] = (s, t, D)
        print("k=%2d cells(D>=1)=%3d   p_tight=%d  p_proof=%d  p_crude=%d"
              "   max Q/k^D=%s (Q<=k^D: %s)   [%.1f s]"
              % (k, len(cells), pmax['tight'], pmax['proof'], pmax['crude'],
                 qmax, qok, time.time() - t0))
        for name in ('tight', 'proof', 'crude'):
            print("     %-5s worst ratio %.4g at (s,t,D)=%s" %
                  (name, float(rmax[name]), argm[name]))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
