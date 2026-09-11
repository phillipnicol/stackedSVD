/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

/- This file sits outside the build.

To check it: copy it to `lean/StackedSVD/RMT/General/Edge/X8.lean`, then run
`LEAN_J=2 "$LEAN_TOOLS/lean-local.sh" <that absolute path>` (the development server).

It compiled on 2026-09-10 13:10 EDT with five `sorry` warnings. `docs/SORRIES.md` does
not track it, because it sits outside `lean/StackedSVD/`. -/
import StackedSVD.RMT.General.Edge.Count

/-! # X8, the skeleton of the Furedi-Komlos cell count

The plan is `notes/x8_plan.md`. The target is `cellCard_mul_le_pos` of
`Count.lean`, restated here as `cellCard_mul_le_pos'` with the exponent `12` in place of
`6`. The route is the absolute count of one cell:

1. `cellCard_le_code`: the step encoding (innovative, forced by the Eulerian parity rule,
   bad) bounds the cell above by the number of codes; at most `4 D` steps are bad.
2. `descFac_mul_min_pow_le`: one unit of deficiency buys one more row label at the price 2.
3. `choose_mul_le_pow_choose`: the two binomial factors of the code are at most `k ^ D`
   times the Narayana product of the target cell.
4. `tree_cell_lower`: the target cell holds at least `C(k,s) C(k,s-1) / k ^ 3` labeled tree
   walks (the Narayana lower bound, with a `k ^ 2` slack).
5. `code_pow_arith`: the polynomial factors fit in `(2k) ^ (4D)`.

The exponent is `p = 2 c + 4` where `c` is the constant of the bad-step bound `#bad <= c D`.
Here `c = 4`, so `p = 12`. The sharp `c = 3` holds for every shape with `k <= 8` and gives
`p = 10`; then the numerals `4 * D + 1`, `(2 * k) ^ 8` and `(2 * k) ^ 12` below read
`3 * D + 1`, `(2 * k) ^ 6` and `(2 * k) ^ 10`, and nothing else changes.

Every lemma above is a `sorry` here. The final theorem is proved from them, so the exponent
`12` and the shape of every hypothesis are settled. -/

open Finset

namespace StackedSVD
namespace Edge

/-- L2a, the step encoding. A walk of the cell `(s, t)` is determined by the positions of
its innovative steps (the row ones among the `k` even steps, the column ones among the `k`
odd steps), the positions of its bad steps (at most `4 D` of them), the endpoint of each bad
step (one of the `k + 1` visited vertices, so at most `2 k`), and the labels of the visited
vertices in first-visit order. The bad steps number at most `4 D`, so their positions and
endpoints cost at most `(2 k) ^ (8 D)`. -/
theorem cellCard_le_code {n d k s t : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hst : s + t ≤ k) (hnd : 2 * k ≤ min n d) :
    cellCard n d k s t
      ≤ k.choose s * k.choose (t - 1) *
          ((4 * (k + 1 - (s + t)) + 1) * ((2 * k) ^ 8) ^ (k + 1 - (s + t))) *
          (n.descFactorial s * d.descFactorial t) := by
  sorry

/-- L2b, the labels. One unit of deficiency turns one factor `min n d` into one more row
label, at the price 2: `min n d ≤ 2 (n - k)` follows from `2 k ≤ min n d ≤ n`. -/
theorem descFac_mul_min_pow_le {n d s D k : ℕ} (hk : 2 ≤ k) (h : s + D ≤ k)
    (hnd : 2 * k ≤ min n d) :
    n.descFactorial s * (min n d) ^ D ≤ 2 ^ D * n.descFactorial (s + D) := by
  sorry

/-- L2c, the binomial step. `k.choose (t - 1) = k.choose (s + D)` by symmetry, and
`k.choose s ≤ k ^ (D - 1) * k.choose (s + D - 1)` by `D - 1` steps of
`k.choose a ≤ k * k.choose (a + 1)`. -/
theorem choose_mul_le_pow_choose {k s t D : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hD : 1 ≤ D) (hsum : s + t + D = k + 1) :
    k * (k.choose s * k.choose (t - 1))
      ≤ k ^ D * (k.choose (s + D) * k.choose (s + D - 1)) := by
  sorry

/-- L2d, the Narayana lower bound on the target cell. A cell with `s + t = k + 1` holds the
tree walks; there are `N (k, s) = C(k,s) C(k,s-1) / k` tree shapes and
`n.descFactorial s * d.descFactorial t` labelings of each. The `k ^ 3` on the right leaves a
factor `k ^ 2` of slack, so a lower bound that loses `k ^ 2` against Narayana is enough. -/
theorem tree_cell_lower {n d k s t : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hsum : s + t = k + 1) (hnd : 2 * k ≤ min n d) :
    k.choose s * k.choose (s - 1) * (n.descFactorial s * d.descFactorial t)
      ≤ k ^ 3 * cellCard n d k s t := by
  sorry

/-- L2e, the polynomial factors of the code fit in `(2 k) ^ (4 D)`: at `D = 1` and `k = 2`
the two sides are 160 and 512. -/
theorem code_pow_arith {k D : ℕ} (hk : 2 ≤ k) (hD : 1 ≤ D) :
    (4 * D + 1) * 2 ^ D * k ^ D * k ^ 3 ≤ ((2 * k) ^ 4) ^ D * k := by
  sorry

/-- X8 with the exponent 12: the count in one cell. Same hypotheses and conclusion as
`cellCard_mul_le_pos` of `Count.lean`, with `(2 k) ^ 12` in place of `(2 k) ^ 6`. -/
theorem cellCard_mul_le_pos' {n d k s t : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hst : s + t ≤ k) (hnd : 2 * k ≤ min n d) :
    cellCard n d k s t * (min n d) ^ (k + 1 - (s + t))
      ≤ ((2 * k) ^ 12) ^ (k + 1 - (s + t)) * cellCard n d k (s + (k + 1 - (s + t))) t := by
  set D := k + 1 - (s + t) with hDdef
  have hD1 : 1 ≤ D := by omega
  have hsum : (s + D) + t = k + 1 := by omega
  have hkpos : 0 < k := by omega
  have hpow : ((2 * k : ℕ) ^ 4) ^ D * ((2 * k : ℕ) ^ 8) ^ D = ((2 * k : ℕ) ^ 12) ^ D := by
    rw [← mul_pow]; congr 1; ring
  refine Nat.le_of_mul_le_mul_left ?_ hkpos
  calc k * (cellCard n d k s t * min n d ^ D)
      ≤ k * ((k.choose s * k.choose (t - 1) *
          ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          (n.descFactorial s * d.descFactorial t)) * min n d ^ D) := by
        exact Nat.mul_le_mul le_rfl
          (Nat.mul_le_mul (cellCard_le_code hk hs ht hst hnd) le_rfl)
    _ = (k * (k.choose s * k.choose (t - 1))) * ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (n.descFactorial s * min n d ^ D) := by ring
    _ ≤ (k * (k.choose s * k.choose (t - 1))) * ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (2 ^ D * n.descFactorial (s + D)) := by
        exact Nat.mul_le_mul le_rfl (descFac_mul_min_pow_le hk (by omega) hnd)
    _ ≤ (k ^ D * (k.choose (s + D) * k.choose (s + D - 1))) *
          ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (2 ^ D * n.descFactorial (s + D)) := by
        exact Nat.mul_le_mul (Nat.mul_le_mul (Nat.mul_le_mul
          (choose_mul_le_pow_choose hk hs ht hD1 (by omega)) le_rfl) le_rfl) le_rfl
    _ = ((4 * D + 1) * 2 ^ D * k ^ D) * ((2 * k) ^ 8) ^ D *
          (k.choose (s + D) * k.choose (s + D - 1) *
            (n.descFactorial (s + D) * d.descFactorial t)) := by ring
    _ ≤ ((4 * D + 1) * 2 ^ D * k ^ D) * ((2 * k) ^ 8) ^ D *
          (k ^ 3 * cellCard n d k (s + D) t) := by
        exact Nat.mul_le_mul le_rfl (tree_cell_lower hk (by omega) ht hsum hnd)
    _ = ((4 * D + 1) * 2 ^ D * k ^ D * k ^ 3) *
          (((2 * k) ^ 8) ^ D * cellCard n d k (s + D) t) := by ring
    _ ≤ (((2 * k) ^ 4) ^ D * k) * (((2 * k) ^ 8) ^ D * cellCard n d k (s + D) t) := by
        exact Nat.mul_le_mul (code_pow_arith hk hD1) le_rfl
    _ = k * (((2 * k) ^ 12) ^ D * cellCard n d k (s + D) t) := by
        rw [← hpow]; ring

end Edge
end StackedSVD
