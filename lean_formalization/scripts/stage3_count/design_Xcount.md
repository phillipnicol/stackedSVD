# Unit Xcount: the design of the Furedi-Komlos count in `Count.lean`

Answer first: the full step encoding does not fit the session. The file is restructured so that
`clsCard_mul_le` is a proved consequence of ONE named `sorry`, `exists_fkEncode`, which is the
encoding as an explicit injection. Everything else in the file is proved.

## The shape of the reduction

`clsCard n d k v = (walkSet n d k v).card` by `rfl`, with
`walkSet n d k v = univ.filter (fun p => NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = v)`.

Write `D = k + 1 - v` and `m = min n d`. Then

```
card (walkSet n d k v ×ˢ univ(Fin D → Fin m))          = clsCard n d k v * m ^ D
card (univ(Fin D → Fin ((2k)^6)) ×ˢ walkSet n d k (k+1)) = ((2k)^6) ^ D * clsCard n d k (k+1)
```

so `Finset.card_le_card_of_injOn` turns the theorem into: an injection

```
F : (walk with v vertices, NoSingle) × (Fin D → Fin m)
      → (Fin D → Fin ((2k)^6)) × (walk with k+1 vertices, NoSingle)
```

that lands in the tree class. `D = 0` is proved separately (both sides are
`clsCard n d k (k+1)`), so `exists_fkEncode` may assume `v ≤ k`, that is `1 ≤ D`.

## What the injection must do (the content left as `sorry`)

Furedi and Komlos 1981 page 237; Anderson, Guionnet and Zeitouni Lemma 2.1.23; Tao section 2.3.
Walk the closed walk `j 0 → i 0 → j 1 → ... → j (k-1) → i (k-1) → j 0` once and classify each of
its `2k` steps:

1. innovative: the endpoint is visited for the first time. There are exactly `v - 1` of these and
   they build a spanning tree of the walk graph, with the same row and column counts as the walk.
2. matched: the step retraces the innovative step of the same edge that is still open, in the
   bracket order of a depth first search.
3. bad: everything else. With no entry of multiplicity 1 there are at most `3 D` bad steps.

The tree skeleton has `v - 1 = k - D` edges, so the image walk needs `D` more edges: attach `D`
pendant edges at fresh labels taken from the `(Fin D → Fin m)` input, alternating sides so that
the target cell `(s + a, t + b)`, `a + b = D`, stays nonempty (the tree class is nonempty in every
cell `s + t = k + 1` with `s, t >= 1`, and `2k <= m` gives at least `k + 1` free labels per side).
The code records, per unit of deficiency, the position of a bad step (at most `2k` values), its
endpoint among the visited vertices (at most `k + 1 <= 2k` values) and the collision repair of a
supplied label that is already used (at most `k + 1 <= 2k` values); three bad steps per unit of
deficiency at `(2k)^2` each is `(2k)^6`, the exponent of `fkFactor`.

## Why the cheap route is not used

The one-step merge count is refuted in section 5 of the unit plan: the walk `i = (0,1,0,1)`,
`j = (0,1,0,1)` at `k = 4` (the 4-cycle walked twice in the same direction, every multiplicity 2,
`v = 4 = k`) has no vertex split that keeps every multiplicity away from 1, so the fiber of the
merge map over it is empty. The encoding above is not a local surgery and does not have this hole.

## Numeric checks reused

`verify2.py` of the unit plan: the `D`-step statement has 0 violations for `k <= 7` at `(n, d)`
from `(k+1, k+1)` to `(10^12, 10^9)`, and it still holds with `(2k)^2` in place of `(2k)^6`. No
randomness, so no seed. This session added `count_checks.py` beside this file for the four new
lemmas (`sum_walkMult`, `card_walkEdges_le`, `two_le_walkVerts`, `walkVerts_le_edges_succ`).

## What this session actually landed (2026-09-10, 13:20 EDT)

`Count.lean` is 604 lines and compiles with ONE `sorry`, in `cellCard_mul_le_pos`.

The reduction that is now proved:

1. `clsCard n d k v = (walkSet n d k v).card` (`rfl`), and `walkSet` splits by row count into
   the cells `cellSet n d k s t`: `clsCard_eq_sum_cellCard`.
2. `cellCard_mul_le` sends every extra vertex to the ROW side: cell `(s, t)` is compared with
   cell `(s + D, t)`, `D = k + 1 - (s + t)`.  The index map `s -> s + D` is injective and lands
   in `range (k + 2)`, and `(s + D) + t = k + 1`, so the terms of the two decompositions match
   one for one and no `4 ^ k` is lost.
3. The degenerate cells (`s = 0` or `t = 0`) are proved empty, so the `sorry` carries
   `2 <= k`, `1 <= s`, `1 <= t`.
4. `one_le_cellCard`: the target cell is nonempty, by the explicit caterpillar walk.

What is left in `cellCard_mul_le_pos`, in the order a next agent should take it:

- a. the labels: `falling n s * min n d ^ D <= 2 ^ D * falling n (s + D)`, since
  `min n d <= n` and `n - s >= n - k >= n / 2` at `2k <= min n d`.  Pure arithmetic.
- b. the factorization `cellCard n d k s t = Sh k s t * falling n s * falling d t`, where
  `Sh` counts the walks up to a separate relabeling of the rows and of the columns.  This is
  an orbit count; it is the part with no mathematical risk and a real Lean cost.
- c. the shape bound `Sh k s t <= (2k)^(4 D)` for `D >= 1`, which is Furedi-Komlos proper.
  `one_le_cellCard` then supplies the `Tr(s + D, t) >= 1` on the other side.

Numerics for the split: `cells.py` (0 violations in 560 cells) and `validate.py` (the tree
shapes reproduce the Narayana numbers and the Catalan totals for `k = 1..7`, and the total
28549 NoSingle shapes at `k = 7` agrees with the unit plan's independent enumeration).  The
factor the cell form actually needs is at most `(2k)^1.37` at `k = 7`, against `(2k)^6`.
