#!/usr/bin/env python3
"""Xcount2: the doubled 2L-cycle refutes "#bad <= 3 D" for the first-visit rules.

Exact, no randomness, so there is no seed.

The walk: rows r_0..r_{L-1}, columns c_0..c_{L-1}, the 2L-cycle
c_0 r_0 c_1 r_1 ... c_{L-1} r_{L-1} c_0 traversed twice.  So k = 2L,
i = (r_0,...,r_{L-1},r_0,...,r_{L-1}), j = (c_0,...,c_{L-1},c_0,...,c_{L-1}).
Every entry has multiplicity 2 (NoSingle holds), v = 2L = k, D = k+1-v = 1.

Printed per L: #bad under the parent rule, the stack rule and the odd (Eulerian
parity) rule, against D = 1.

walk_mult, vertex_seq, edge_of, bad_parent, bad_stack and bad_odd are inlined here
from the original scratch module rules.py (Xcount2 stage-3 session, 2026-09-10) so
this script is standard library only and needs no sibling file.
"""
from collections import Counter


def walk_mult(i, j, k):
    cnt = Counter()
    for t in range(k):
        cnt[(i[t], j[t])] += 1
        cnt[(i[t], j[(t + 1) % k])] += 1
    return cnt


def vertex_seq(i, j, k):
    u = []
    for t in range(k):
        u.append(('c', j[t]))
        u.append(('r', i[t]))
    return u


def edge_of(x, y):
    return (x, y) if x[0] == 'r' else (y, x)


def bad_parent(u, k):
    n2 = 2 * k
    visited = {u[0]}
    parent = {}
    nb = 0
    for s in range(n2):
        x = u[s]
        y = u[(s + 1) % n2]
        if y not in visited:
            visited.add(y)
            parent[y] = x
        elif parent.get(x, None) == y:
            pass
        else:
            nb += 1
    return nb


def bad_stack(u, k):
    n2 = 2 * k
    visited = {u[0]}
    stack = []
    nb = 0
    for s in range(n2):
        x = u[s]
        y = u[(s + 1) % n2]
        if y not in visited:
            visited.add(y)
            stack.append(x)
        elif stack and stack[-1] == y:
            stack.pop()
        else:
            nb += 1
    return nb


def bad_odd(u, k):
    n2 = 2 * k
    visited = {u[0]}
    odd = {}                                # vertex -> set of odd incident edges
    nb = 0
    for s in range(n2):
        x = u[s]
        y = u[(s + 1) % n2]
        ox = odd.get(x, set())
        if y not in visited:
            visited.add(y)
        elif len(ox) == 1:
            e = next(iter(ox))
            other = e[0] if e[1] == x else e[1]
            if other != y:
                nb += 1
        else:
            nb += 1
        e = edge_of(x, y)
        for z in (x, y):
            sz = odd.setdefault(z, set())
            if e in sz:
                sz.discard(e)
            else:
                sz.add(e)
    return nb


def main():
    print("=== doubled 2L-cycle: k = 2L, v = k, D = 1 (exact, no seed) ===")
    print(" L    k   v   D  | #bad parent  #bad stack  #bad odd")
    for L in range(2, 13):
        k = 2 * L
        i = tuple(list(range(L)) + list(range(L)))
        j = tuple(list(range(L)) + list(range(L)))
        cnt = walk_mult(i, j, k)
        assert all(c == 2 for c in cnt.values()), sorted(cnt.values())
        v = len(set(i)) + len(set(j))
        D = k + 1 - v
        u = vertex_seq(i, j, k)
        print("%2d  %3d %3d %3d  |    %3d         %3d        %3d"
              % (L, k, v, D, bad_parent(u, k), bad_stack(u, k), bad_odd(u, k)))


if __name__ == '__main__':
    main()
