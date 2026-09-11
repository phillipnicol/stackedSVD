# Unit Xcount2: the design of the second, independent attempt at X8

Status: the arithmetic reduction is proved and compiles. The combinatorial core is one
tracked `sorry` (`clsCard_step_core`). Every number below is exact integer or `Fraction`
arithmetic over a complete enumeration, so there is no seed to print. Scripts and outputs
sit beside this file (`steps.py`, `rules.py`, `cycle.py`, `shapes_odd.py`, `cells.py`,
`shapes_odd.out`).

## 0. Answers first

1. **The bad-step bound in the `Count.lean` docstring is false as written.** The docstring
   says "a walk with no single entry has at most `3 D` bad steps" for the first-visit tree
   classification (innovative, matched, bad). The doubled `2L`-cycle has `D = 1` and
   `k + 1 = 2L + 1` bad steps under that rule, so no constant works. Section 2.
2. **One rule change repairs it.** Call a step forced when the current vertex has exactly one
   incident entry traversed an odd number of times so far, and that entry leads to the
   destination (Eulerian parity, not the tree parent). Over every NoSingle shape with
   `k ≤ 7`, the worst `#bad / D` is `8/3 = 2.67`, and `#bad = 2` at every `D = 1` shape.
   Section 3.
3. **The cell grading cannot be dropped, and the cell-wise bound the route needs is true**
   with a margin of `1.8e5` at `k = 7`. Section 4.
4. **The route is not finishable in one session.** The reduction of the `D`-step form to `D`
   copies of a one-step form is proved (`Count2.lean`, 105 lines, exit 0). The one-step form
   itself is the single `sorry`. Section 5 gives the decomposition it still needs.

## 1. What the file proves

`Count2.lean` states `clsCard_mul_le'`, the exact statement of `clsCard_mul_le`, and proves
it from a one-step form by induction on the deficiency `D = k + 1 - v`:

```
clsCard n d k v * min n d  ≤  (2k)^6 * clsCard n d k (v+1)          (clsCard_step)
```

The chain loses nothing: `D` copies of a factor `(2k)^6 / min(n,d)` give exactly
`((2k)^6)^D / (min n d)^D`. Below `v = 2` the class is empty (a walk of positive length uses
at least one row and at least one column), so the `sorry` carries the hypothesis `2 ≤ v ≤ k`.

The one-step form is not the form the unit X plan recommends. The plan warns that the
one-step form "invites the merge-count proof, which is refuted". That warning is about a
*proof*, not about the *statement*: the one-step statement is true (0 violations, smallest
working exponent 2, unit X plan section 2(b) and `steps.py` here), and the merge-count proof
is refuted whichever form is stated. Splitting the induction out separates 105 lines of
arithmetic from the combinatorics and makes the remaining gap one lemma about two adjacent
classes.

## 2. The refutation: `#bad ≤ 3 D` fails for the first-visit rules

Take the `2L`-cycle `c_0 r_0 c_1 r_1 ... c_{L-1} r_{L-1} c_0` traversed twice:

```
k = 2L,   i = (r_0,...,r_{L-1}, r_0,...,r_{L-1}),   j = (c_0,...,c_{L-1}, c_0,...,c_{L-1})
```

Every entry has multiplicity 2, so `NoSingle` holds; `v = 2L = k` and `D = k + 1 - v = 1`.
`cycle.py` prints the bad-step count of three rules on this family:

| L | k | D | parent rule | stack rule | parity rule |
|---|---|---|---|---|---|
| 2 | 4 | 1 | 5 | 4 | 2 |
| 4 | 8 | 1 | 9 | 8 | 2 |
| 8 | 16 | 1 | 17 | 16 | 2 |
| 12 | 24 | 1 | 25 | 24 | 2 |

"parent" is the rule of the `Count.lean` docstring: a step is free when the destination is
new (innovative) or is the parent of the current vertex in the tree of the innovative steps.
"stack" is the bracket-matching variant: free when the destination is the source on top of
the stack of open innovative steps. Both give `#bad = k + 1` and `#bad = k` at `D = 1`, so
`#bad ≤ c D` is false for every constant `c`.

The reason is structural, not an artifact of the family. On the second traversal the walk
moves *away* from the root along tree edges it has already used, and both first-visit rules
only ever free a move *toward* the root.

Exhaustive confirmation, `shapes_odd.py`, every NoSingle shape (a walk up to a separate
injective relabeling of rows and of columns; the shape count reproduces the unit X plan
exactly, 28549 at `k = 7`):

| k | shapes | worst `#bad/D` parent | stack | parity | max `#bad` at `D=1`, parity |
|---|---|---|---|---|---|
| 4 | 62 | 5 | 4 | 2 | 2 |
| 5 | 398 | 5 | 5 | 2 | 2 |
| 6 | 3108 | 7 | 6 | 7/3 | 2 |
| 7 | 28549 | 7 | 7 | 8/3 | 2 |

## 3. The rule that works: Eulerian parity

Index the `2k` steps of the closed walk by `s = 0 .. 2k-1`, with the cyclic vertex sequence
`u_{2t} = column (j t)`, `u_{2t+1} = row (i t)`. Step `2t` traverses the entry `(i t, j t)`
and step `2t+1` traverses `(i t, j (cycSucc t))`, so the `2k` steps are exactly the `2k`
factors that `walkMult` counts.

Classify step `s`, from `u_s` to `u_{s+1}`:

- **I** (innovative): `u_{s+1}` is not in `{u_0, ..., u_s}`.
- **F** (forced): `u_{s+1}` is old, the vertex `u_s` has exactly one incident entry traversed
  an odd number of times among the steps `0 .. s-1`, and that entry leads to `u_{s+1}`.
- **B** (bad): otherwise. The code records `u_{s+1}`.

Why the parity rule frees the second traversal of a cycle: after `s` steps, the set `O_s` of
entries traversed an odd number of times has odd degree exactly at `u_0` and at `u_s`. So at
`u_s ≠ u_0` the degree in `O_s` is odd, hence at least 1, and the rule fires whenever it is
exactly 1. On the second traversal of a doubled cycle every vertex has one incident odd entry
(the edge ahead), so every such step is forced.

The code is `(type sequence, the new vertices in first-visit order, the destinations of the
bad steps)`. The decoder replays the walk: at an **I** step it takes the next new label, at an
**F** step it recomputes `O_s` from the history and follows the unique odd entry, at a **B**
step it reads the recorded destination. So the encoding is injective by construction; the
same reconstruction was run and checked for the stack rule on every walk with `k ≤ 5`,
`n, d ≤ 5` (`steps.py`: 0 reconstruction failures, 0 code collisions, and `#I = v - 1` with 0
failures over 1.4 million NoSingle walks).

## 4. The accounting, and why the cell grading is not optional

Write `S(k, s, t)` for the number of NoSingle shapes with `s` rows and `t` columns,
`N(k, s) = C(k,s) C(k,s-1) / k` for the Narayana number (the tree shapes with `s` rows, so
`t = k + 1 - s`), and `D = k + 1 - s - t`. The labeled count splits as
`clsCard n d k v = Σ_{s+t=v} S(k,s,t) * n^(s) * d^(t)` with falling factorials, and the tree
class is `clsCard n d k (k+1) = Σ_s N(k,s) * n^(s) * d^(k+1-s)`.

The route needs the **cell-wise** shape bound

```
(1)     S(k, s, t)  ≤  ((2k)^6 / 2)^D * N(k, s).
```

Given (1), the one-step form follows from one falling-factorial step that keeps `s` fixed and
spends the deficiency on the columns: `d^(t+1) = d^(t) (d - t)` and `t ≤ k ≤ min(n,d)/2 ≤ d/2`
give `min(n,d) ≤ d ≤ 2 (d - t)`, hence `n^(s) d^(t) min(n,d) ≤ 2 n^(s) d^(t+1)`.

`cells.py` checks (1) on every cell with `k ≤ 7`. **0 violations**; the worst ratio is
`5.6e-06` at `k = 7`, cell `(s,t) = (1,6)`, `D = 1`, so the margin is about `1.8e5`. The
`D = 0` row of the table reproduces the Narayana numbers exactly, which validates the
enumerator.

**A cell-blind bound fails, and the failure is not far away.** Replacing (1) by
`Σ_{s+t=v} S(k,s,t) ≤ ((2k)^6/2)^D * N(k,1)` is what a bound that forgets `(s,t)` has to
survive, because at `n ≫ d` the tree class is dominated by the cell `s = k`, `t = 1`, where
`N(k,k) = 1`. That ratio is 2.5e-04, 3.4e-04, 5.6e-04, 1.0e-03 at `k = 4,5,6,7`: it grows by
about 1.85 per unit of `k` while the denominator grows polynomially, so it crosses 1 near
`k ≈ 18`. Stage 3 runs at `k = momOrder C d = ⌈C log d⌉`, which passes 18 quickly, so the
cell grading has to be carried through the whole encoding, exactly as unit X plan section 6
note 2 says.

## 5. What is left, as a decomposition

| # | Statement | Estimate |
|---|---|---|
| Y1 | `cellCard n d k s t`, and `clsCard n d k v = Σ_{s+t=v} cellCard n d k s t` | 150 |
| Y2 | `cellCard n d k s t = S k s t * n^(s) * d^(t)` (the shape / label split) | 500 |
| Y3 | `clsCard n d k (k+1) = Σ_s N k s * n^(s) * d^(k+1-s)` (the tree class, exactly) | 400 |
| Y4 | the parity step classification, and `#I = v - 1` | 600 |
| Y5 | `#B ≤ 3 D` under the parity rule (the one real theorem) | 900 |
| Y6 | the reconstruction: a shape is determined by its code | 700 |
| Y7 | (1) assembled, then the falling-factorial step and `clsCard_step` | 400 |

Y5 is the item to attack first, because it is the one that the first-visit rules get wrong
and the one no numeric check can replace. Nothing above needs `1 ≤ K`, `hn` or `hd`.

## 6. What was ruled out

- The one-vertex merge with a crude split fiber (unit X plan section 5): refuted there, and
  reconfirmed here on the doubled 4-cycle `i = j = (0,1,0,1)`, `k = 4`.
- Relabeling two column slots at once to a fresh label: on the doubled 4-cycle the only
  admissible pair of slots is the whole preimage of one column, so the old column disappears
  and `v` does not rise.
- A crude pattern count (`s^k` surjections for the rows, `t^k` for the columns): it costs
  `k^{2k}` where the budget is `(2k)^{6D}`, so it dies at `D = 1`.
- A bound through the total shape count instead of the cell (section 4).
