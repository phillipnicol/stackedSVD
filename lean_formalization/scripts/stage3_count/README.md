# Stage 3 numeric evidence: the Furedi-Komlos count

`cellCard_mul_le_pos` is proved (2026-09-10 evening), in
`lean/StackedSVD/RMT/General/Edge/CountBound.lean` (moved from `Count.lean`). It proves one
step of the sharp upper edge of Bai-Yin's law at four moments: the Furedi-Komlos bound
on walk counts in the moment method. A closed walk of length `2k` on rows times columns,
with no entry visited exactly once, is grouped by its row count `s` and column count `t`
(`s + t <= k`). Then `D = k + 1 - s - t` extra vertices are still needed to reach a tree
walk. The claim: the count of walks in cell `(s, t)`, times `min(n, d)^D`, is at most
`((2k)^12)^D` times the count of walks in the tree cell `(s + D, t)`; the exponent was 6
until 2026-09-10 evening, when the X8 plan (`notes/x8_plan.md`, a six-part code of a walk
and a Narayana lower bound on the tree cell) replaced it with 12. The six scripts below
still check the earlier exponent 6; each finds a wide margin, so both values hold.

## Scripts

All six are exact integer (or `Fraction`) enumeration; none use randomness.

- `cells.py`: claim L2 (the cell bound) and claim TR (every tree cell nonempty), `k = 1..7`, 4 s. 560 cells tested, 0 violations, worst ratio 4.9e-06.
- `validate.py`: the shape counts against the Narayana and Catalan numbers, and the smallest exponent that survives, `k = 1..7`, 12 s. Exact match at every `k`; needed exponent only `(2k)^1.37` at `k = 7`, against the `(2k)^6` the proof used before 2026-09-10 evening (now `(2k)^12`, a wider margin still).
- `count_checks.py`: all 10 lemmas and definitions of `Count.lean` by direct enumeration, `k` up to 5 (up to 15 for the caterpillar witness), 2 s. 192,741 walks, 256 caterpillar cases, 0 violations.
- `cycle.py`: the doubled `2L`-cycle family at `D = 1`, `k = 4..24`, under 1 s. Bad-step count `k + 1` (parent rule), `k` (stack rule), 2 (Eulerian parity rule).
- `verify_step.py`: the one-step form of the independent `Count2.lean` route, `clsCard * min(n,d) <= (2k)^6 * clsCard(v+1)`, `k = 2..7`, 17 s. 0 mismatches, 0 violations, worst ratio 2.0e-04 at `k = 2`.
- `cells_onestep.py`: the cell-wise bound that same route needs, `S(k,s,t) <= ((2k)^6/2)^D * N(k,s)`, `k = 2..7`, 6 s. 0 violations, worst ratio 5.6e-06 at `k = 7`.

Runtimes are wall-clock, `ulimit -v 16000000`, 2 cores. No script needed a cap; the
slowest (`verify_step.py`) finishes in 17 seconds. Outputs sit beside each script as
`output_<name>.txt`. Design notes: `design_Xcount.md` (the merged `Count.lean` route)
and `design_Xcount2.md` (the independent, not-merged one-step route, `Count2.lean`).

## X8: the proof plan evidence

Three more scripts back `notes/x8_plan.md`, the proof plan that proved
`cellCard_mul_le_pos` (the per-cell step count) with the exponent `p = 12` in place of `6`
(status: proved, 2026-09-10 evening). Exact integer or `Fraction` enumeration, no
randomness. Run with
no argument for `k = 2..7`; pass an integer argument for a larger `k` (the plan's own run
went to `k = 8`, `notes/x8_plan.md` section 0).

- `x8_enum.py`: every walk shape of length `2k` under the Eulerian parity step rule
  (innovative, forced, bad). Reports the shape count, the Narayana check on the tree
  cells, and the worst `#bad/D`, `#dev/D` and `(2#bad+#dev)/D` ratios, with a witness
  walk for each. `k = 2..7`, 1 s, exit 0. Output: `output_x8_enum.txt`.
- `x8_charge.py`: the charging decomposition behind the bad-step bound. Imports
  `x8_enum.py`. Checks the identity `D = X + c` (`X` the excess multiplicity units, `c`
  the cycle rank) at every shape, and tabulates the worst `#bad` per `(X, c)` and the
  three charging cases (source degree at least 3, exactly 1, or 0, in the odd set).
  `k = 2..7`, 1 s, exit 0. Output: `output_x8_charge.txt`.
- `x8_cells8.py`: the exponent `p` each cell needs, three accountings (tight: distinct
  type words counted exactly; proof: the binomial bound Lean can prove; crude: the
  closed form the Lean proof uses) against the Narayana count of the target cell.
  Imports `x8_enum.py`. `k = 2..7`, 2 s, exit 0. Output: `output_x8_cells8.txt`.

Result (`notes/x8_plan.md` section 0, full `k <= 8` sweep): max `#bad/D = 3` at `k = 8`
under the Eulerian parity rule, with `D = 4` and 12 bad steps; max `#dev/D = 4/3`; max
`(2 #bad + #dev)/D = 6` for `k <= 8`. The per-cell exponent needed is 3 by the tight
accounting, 6 by the exact binomial, 8 by the closed form; the plan takes `p = 12`.

`X8_skeleton.lean` is the compiled Lean skeleton of the plan (five `sorry` warnings, the
assembly theorem proved from them). It is not part of the build; its header comment says
how to check it.

## Caution

`design_Xcount.md` proposes a bracket-order matched-step rule for the encoding; the proof
of `cellCard_mul_le_pos` (the per-cell step count, proved 2026-09-10 evening) does not use
it. `cycle.py` shows that rule is false: its bad-step
count grows with `k`, not with a constant times `D = 1`. The stack (bracket-matching)
variant fails the same way. Only the Eulerian parity rule keeps the bad-step count
bounded (2, at every `k` tested), and the proof (`Code.lean`, the six-part step code) uses
that rule, not bracket matching.

## x8_scratch/ - the X8 plan's scratch Lean files

`x8_scratch/` holds the scratch Lean files and agent reports behind `notes/x8_plan.md`, the
proof plan that proved `cellCard_mul_le_pos` (the per-cell step count). None of it is part of the Lake build; each
`.lean` file was checked on its own with `lean-local.sh` against the Mathlib snapshot on
2026-09-10. See `x8_scratch/README.md` for the file list and the status of each plan unit.
