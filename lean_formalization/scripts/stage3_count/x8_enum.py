#!/usr/bin/env python3
"""X8 numerics: every NoSingle walk shape of length 2k, the parity step rule.

Exact enumeration by depth first search, no randomness, so there is no seed.

A shape is a closed walk u_0 u_1 ... u_{2k-1} u_0 on the complete bipartite graph,
u_0 = column 0, even index = column, odd index = row, with the rows labeled in the
order of first visit and the columns likewise.  That is the canonical form of a walk
under a separate injective relabeling of rows and of columns, so shape counts here
match the rgs(k) x rgs(k) enumeration of the earlier units.

Step s goes from u_s to u_{s+1} and traverses the entry (row, column) of the pair.
Classification (Eulerian parity rule):
  I  u_{s+1} is new;
  F  u_{s+1} is old, u_s has exactly one incident entry with an odd number of
     traversals among the steps 0..s-1, and that entry leads to u_{s+1};
  B  otherwise.
A deviation is an F step whose source u_s still has an innovative step of its own
later in the walk (an unvisited child in the plane tree of innovative steps).

Prints, per k: the shape count, max #bad/D, max #dev/D, max (2 #bad + #dev)/D.
"""
import sys, time
from collections import Counter
from fractions import Fraction
from math import comb


def enumerate_shapes(k):
    """Yield (types, nbad, ndev, s, t) for every NoSingle shape of length 2k."""
    n2 = 2 * k
    out = []
    # state: u list, mult dict edge->count, nr, nc
    u = [('c', 0)]
    mult = Counter()
    ones = [0]          # number of edges with multiplicity exactly 1

    def edge_of(s, x, y):
        # step s from x to y; s even: x column, y row; s odd: x row, y column
        return (y[1], x[1]) if s % 2 == 0 else (x[1], y[1])

    def rec(s, nr, nc):
        if s == n2:
            if ones[0] == 0:
                out.append((tuple(u), nr, nc))
            return
        x = u[s]
        if s % 2 == 0:
            cands = [('r', a) for a in range(nr)] + [('r', nr)]
        else:
            if s == n2 - 1:
                cands = [('c', 0)]
            else:
                cands = [('c', a) for a in range(nc)] + [('c', nc)]
        rem = n2 - 1 - s
        for y in cands:
            e = edge_of(s, x, y)
            c0 = mult[e]
            # prune: at most k distinct edges (every multiplicity is at least 2)
            if c0 == 0 and len(mult) >= k:
                continue
            mult[e] = c0 + 1
            if c0 == 0:
                ones[0] += 1
            elif c0 == 1:
                ones[0] -= 1
            if ones[0] <= rem:
                u.append(y)
                nr2 = nr + 1 if (s % 2 == 0 and y[1] == nr) else nr
                nc2 = nc + 1 if (s % 2 == 1 and y[1] == nc) else nc
                rec(s + 1, nr2, nc2)
                u.pop()
            if c0 == 0:
                ones[0] -= 1
                del mult[e]
            else:
                mult[e] = c0
                if c0 == 1:
                    ones[0] += 1

    rec(0, 0, 1)   # the root column c0 is already visited
    return out


def classify(useq, k):
    """types string, #bad, #dev for one shape given as the vertex sequence."""
    n2 = 2 * k
    seq = list(useq) + [useq[0]]
    visited = {seq[0]}
    odd = {}                       # vertex -> set of incident edges of odd count
    types = []
    innov_from = Counter()         # source vertex -> number of innovative steps
    # first pass needs the innovative positions, so do two passes
    vis = {seq[0]}
    innov_pos = []
    for s in range(n2):
        y = seq[s + 1]
        if y not in vis:
            vis.add(y)
            innov_pos.append(s)
    later_innov = [0] * (n2 + 1)   # later_innov[s] handled below per vertex
    # per vertex, the list of innovative step indices leaving it
    innov_at = {}
    for s in innov_pos:
        innov_at.setdefault(seq[s], []).append(s)
    nbad = ndev = 0
    for s in range(n2):
        x, y = seq[s], seq[s + 1]
        ox = odd.get(x, set())
        if y not in visited:
            visited.add(y)
            types.append('I')
        else:
            forced = False
            if len(ox) == 1:
                e = next(iter(ox))
                other = ('r', e[0]) if x[0] == 'c' else ('c', e[1])
                if other == y:
                    forced = True
            if forced:
                types.append('F')
                if any(p > s for p in innov_at.get(x, ())):
                    ndev += 1
            else:
                types.append('B')
                nbad += 1
        e = (y[1], x[1]) if s % 2 == 0 else (x[1], y[1])
        for z in (x, y):
            sz = odd.setdefault(z, set())
            if e in sz:
                sz.discard(e)
            else:
                sz.add(e)
    return ''.join(types), nbad, ndev


def narayana(k, s):
    if s < 1 or s > k:
        return 0
    return comb(k, s) * comb(k, s - 1) // k


def main(kmax):
    print("=== X8: shapes, parity rule, bad steps and deviations (exact, no seed) ===")
    for k in range(2, kmax + 1):
        t0 = time.time()
        shapes = enumerate_shapes(k)
        wb = wd = wm = Fraction(0)
        wbw = wdw = wmw = None
        cells = {}                      # (s,t) -> {b: set of type words}
        tree_count = Counter()
        for (useq, nr, nc) in shapes:
            s, t = nr, nc
            D = k + 1 - s - t
            typ, nbad, ndev = classify(useq, k)
            if D == 0:
                tree_count[s] += 1
                continue
            for val, cur, name in ((nbad, wb, 'b'), (ndev, wd, 'd'),
                                   (2 * nbad + ndev, wm, 'm')):
                pass
            r = Fraction(nbad, D)
            if r > wb:
                wb, wbw = r, (useq, nbad, D)
            r = Fraction(ndev, D)
            if r > wd:
                wd, wdw = r, (useq, ndev, D)
            r = Fraction(2 * nbad + ndev, D)
            if r > wm:
                wm, wmw = r, (useq, nbad, ndev, D)
            cells.setdefault((s, t), {}).setdefault(nbad, set()).add(typ)
        ok_tree = all(tree_count[s] == narayana(k, s) for s in range(1, k + 1))
        print("k=%d shapes=%d  treecells=Narayana:%s  max #bad/D=%s  max #dev/D=%s"
              "  max (2#bad+#dev)/D=%s   [%.1f s]"
              % (k, len(shapes), ok_tree, wb, wd, wm, time.time() - t0))
        if wbw:
            print("    bad witness: %s  #bad=%d D=%d" % (fmt(wbw[0]), wbw[1], wbw[2]))
        if wmw:
            print("    mixed witness: %s #bad=%d #dev=%d D=%d"
                  % (fmt(wmw[0]), wmw[1], wmw[2], wmw[3]))
        yield k, cells, shapes


def fmt(useq):
    return ''.join('%s%d' % (a, b) for (a, b) in useq)


if __name__ == '__main__':
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for _ in main(kmax):
        pass
