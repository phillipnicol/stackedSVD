"""Unit Xcount: numeric check of every statement proved in Count.lean.

Exhaustive enumeration, exact integers, no randomness, so there is no seed to print.
walkMult is recomputed by hand exactly as
lean/StackedSVD/RMT/General/Edge/Defs.lean writes it (not parsed from the file):
the 2k factors are Y (i t) (j t) and Y (i t) (j (cycSucc t)) with
cycSucc t = (t + 1) % k. Standard library only, no sibling imports.
"""
from itertools import product

def mult(i, j, k):
    m = {}
    for t in range(k):
        for e in ((i[t], j[t]), (i[t], j[(t + 1) % k])):
            m[e] = m.get(e, 0) + 1
    return m

def catRow(s, x):
    return x + 1 if x + 2 <= s else 0

def catCol(s, x):
    return 0 if x + 1 <= s else x + 1 - s

fails = {n: 0 for n in ["sum_walkMult", "mem_walkEdges", "card_walkEdges_le",
                        "walkVerts_le_edges_succ", "walkVerts_le_of_excess",
                        "isPaired_of_walkVerts_eq", "two_le_walkVerts",
                        "noSingle_of_families", "clsCard_eq_sum_cellCard",
                        "one_le_cellCard"]}
walks = 0

CASES = [(n, d, k) for n in (1, 2, 3) for d in (1, 2, 3) for k in (1, 2, 3, 4)]
CASES += [(4, 3, 4), (3, 4, 4), (4, 4, 4), (2, 3, 5), (3, 2, 5), (3, 3, 5)]
for (n, d, k) in CASES:
    cls, cell = {}, {}
    for i in product(range(n), repeat=k):
        for j in product(range(d), repeat=k):
            walks += 1
            m = mult(i, j, k)
            E = [e for e, c in m.items() if c > 0]
            v = len(set(i)) + len(set(j))
            nosingle = all(c != 1 for c in m.values())
            excess = nosingle and any(c >= 3 for c in m.values())
            # sum_walkMult
            if sum(m.values()) != 2 * k:
                fails["sum_walkMult"] += 1
            # mem_walkEdges_fst / _snd
            for t in range(k):
                if (i[t], j[t]) not in E or (i[t], j[(t + 1) % k]) not in E:
                    fails["mem_walkEdges"] += 1
            # walkVerts_le_edges_succ
            if v > len(E) + 1:
                fails["walkVerts_le_edges_succ"] += 1
            # two_le_walkVerts (k >= 1 always here)
            if v < 2:
                fails["two_le_walkVerts"] += 1
            if nosingle:
                # card_walkEdges_le
                if len(E) > k:
                    fails["card_walkEdges_le"] += 1
                cls[v] = cls.get(v, 0) + 1
                key = (len(set(i)), len(set(j)))
                cell[key] = cell.get(key, 0) + 1
                # isPaired_of_walkVerts_eq
                if v == k + 1 and any(c not in (0, 2) for c in m.values()):
                    fails["isPaired_of_walkVerts_eq"] += 1
            # walkVerts_le_of_excess
            if excess and v > k:
                fails["walkVerts_le_of_excess"] += 1
            # noSingle_of_families: the hypothesis implies the conclusion
            F1 = {(i[t], j[t]) for t in range(k)}
            F2 = {(i[t], j[(t + 1) % k]) for t in range(k)}
            if F1 == F2 and not nosingle:
                fails["noSingle_of_families"] += 1
    # clsCard_eq_sum_cellCard
    for v, c in cls.items():
        tot = sum(cell.get((s, v - s), 0) for s in range(v + 1))
        if tot != c:
            fails["clsCard_eq_sum_cellCard"] += 1

# one_le_cellCard: the caterpillar witness, over a wide range of (s, t, n, d)
cats = 0
for s in range(1, 9):
    for t in range(1, 9):
        k = s + t - 1
        for (n, d) in [(s, t), (s + 3, t + 1), (2 * k, 2 * k), (2 * k + 5, 2 * k + 2)]:
            i = tuple(catRow(s, x) for x in range(k))
            j = tuple(catCol(s, x) for x in range(k))
            cats += 1
            ok = (max(i) < n and max(j) < d and len(set(i)) == s and len(set(j)) == t
                  and all(c != 1 for c in mult(i, j, k).values())
                  and len(set(i)) + len(set(j)) == k + 1)
            if not ok:
                fails["one_le_cellCard"] += 1
                if fails["one_le_cellCard"] <= 5:
                    print(f"  caterpillar FAIL s={s} t={t} k={k} n={n} d={d} i={i} j={j} "
                          f"mult={sorted(mult(i, j, k).values())}")

print(f"walks enumerated: {walks};  caterpillar cases: {cats}")
for name, c in fails.items():
    print(f"  {name}: {c} violations")
