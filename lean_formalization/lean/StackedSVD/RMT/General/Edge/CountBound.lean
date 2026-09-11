/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Code
import StackedSVD.RMT.General.Edge.Dyck

/-! # Stage 3, unit X, part 4: the Furedi-Komlos count

The assembly of the count (`notes/x8_plan.md`, section 6; Furedi and Komlos 1981, page 237).
`cellCard_mul_le_pos` puts the code bound of `Code.lean` over the tree-walk bound of
`Dyck.lean`, `cellCard_mul_le` discharges the degenerate cells, and `clsCard_mul_le` sums the
cells of one vertex count. Three arithmetic lemmas pay for the move to the target cell:
`descFac_mul_min_pow_le` (one unit of deficiency turns one factor `min n d` into one more row
label, at the price 2), `choose_mul_le_pow_choose` (the two binomial factors move from `s` to
`s + D`, at the price `k` per step) and `code_pow_arith` (the polynomial factors).

The budget of `(2 k) ^ 12` per unit of vertex deficiency `D = k + 1 - (s + t)` is `(2 k) ^ 8`
for the code (`cellCard_le_code'` of `Code.lean`) times `(2 k) ^ 4` for everything else. The
second factor is `code_pow_arith`: it absorbs the `4 D + 1` of the bad-step count, the `2 ^ D`
of the labels, the `k ^ D` of the binomial steps and the single `k ^ 3` of the Narayana ratio
(`narayana_lower` of `Dyck.lean`). -/

/-! ## Part A: the three arithmetic lemmas of the plan -/

namespace StackedSVD.Edge

/-- L2b, the labels. One unit of deficiency turns one factor `min n d` into one more row
label, at the price 2: `min n d ≤ 2 (n - k)` follows from `2 k ≤ min n d ≤ n`. -/
theorem descFac_mul_min_pow_le {n d s D k : ℕ} (h : s + D ≤ k)
    (hnd : 2 * k ≤ min n d) :
    n.descFactorial s * (min n d) ^ D ≤ 2 ^ D * n.descFactorial (s + D) := by
  have key : ∀ m : ℕ, s + m ≤ k →
      n.descFactorial s * (min n d) ^ m ≤ 2 ^ m * n.descFactorial (s + m) := by
    intro m
    induction m with
    | zero => intro _; simp
    | succ m ih =>
      intro hsm
      have hm : s + m ≤ k := by omega
      have hstep := ih hm
      have hmn : min n d ≤ n := min_le_left n d
      have hkn : 2 * k ≤ n := le_trans hnd hmn
      have hkey : min n d ≤ 2 * (n - (s + m)) := by omega
      have hprod : n.descFactorial s * (min n d) ^ m * min n d ≤
          2 ^ m * n.descFactorial (s + m) * (2 * (n - (s + m))) :=
        Nat.mul_le_mul hstep hkey
      calc n.descFactorial s * (min n d) ^ (m + 1)
          = n.descFactorial s * (min n d) ^ m * min n d := by rw [pow_succ]; ring
        _ ≤ 2 ^ m * n.descFactorial (s + m) * (2 * (n - (s + m))) := hprod
        _ = 2 ^ (m + 1) * ((n - (s + m)) * n.descFactorial (s + m)) := by
            rw [pow_succ]; ring
        _ = 2 ^ (m + 1) * n.descFactorial (s + m + 1) := by
            rw [Nat.descFactorial_succ]
  exact key D h

end StackedSVD.Edge

namespace StackedSVD.Edge

/-- L2c, the binomial step. `k.choose (t - 1) = k.choose (s + D)` by symmetry, and
`k.choose s ≤ k ^ (D - 1) * k.choose (s + D - 1)` by `D - 1` steps of
`k.choose a ≤ k * k.choose (a + 1)`. -/
theorem choose_mul_le_pow_choose {k s t D : ℕ} (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hD : 1 ≤ D) (hsum : s + t + D = k + 1) :
    k * (k.choose s * k.choose (t - 1))
      ≤ k ^ D * (k.choose (s + D) * k.choose (s + D - 1)) := by
  -- Single-step monotonicity: `k.choose a ≤ k * k.choose (a + 1)` when `a + 1 ≤ k`.
  have step : ∀ a : ℕ, a + 1 ≤ k → k.choose a ≤ k * k.choose (a + 1) := by
    intro a ha
    have heq : k.choose (a + 1) * (a + 1) = k.choose a * (k - a) :=
      Nat.choose_succ_right_eq k a
    have h1 : 1 ≤ k - a := by omega
    have e1 : k.choose a * 1 ≤ k.choose a * (k - a) :=
      mul_le_mul (le_refl _) h1 (Nat.zero_le _) (Nat.zero_le _)
    have e2 : k.choose (a + 1) * (a + 1) ≤ k.choose (a + 1) * k :=
      mul_le_mul (le_refl _) ha (Nat.zero_le _) (Nat.zero_le _)
    calc k.choose a = k.choose a * 1 := (Nat.mul_one _).symm
      _ ≤ k.choose a * (k - a) := e1
      _ = k.choose (a + 1) * (a + 1) := heq.symm
      _ ≤ k.choose (a + 1) * k := e2
      _ = k * k.choose (a + 1) := Nat.mul_comm _ _
  -- Iterate `step` `j` times: `k.choose a ≤ k ^ j * k.choose (a + j)` when `a + j ≤ k`.
  have general : ∀ j a : ℕ, a + j ≤ k → k.choose a ≤ k ^ j * k.choose (a + j) := by
    intro j
    induction j with
    | zero => intro a _; simp
    | succ j ih =>
      intro a ha
      have ha' : a + j ≤ k := by omega
      have ha'' : a + j + 1 ≤ k := by omega
      have h1 : k.choose a ≤ k ^ j * k.choose (a + j) := ih a ha'
      have h2 : k.choose (a + j) ≤ k * k.choose (a + j + 1) := step (a + j) ha''
      have h3 : k ^ j * k.choose (a + j) ≤ k ^ j * (k * k.choose (a + j + 1)) :=
        mul_le_mul (le_refl _) h2 (Nat.zero_le _) (Nat.zero_le _)
      have h4 : k ^ j * (k * k.choose (a + j + 1)) = k ^ (j + 1) * k.choose (a + j + 1) := by
        rw [pow_succ']; ring
      have h5 : a + (j + 1) = a + j + 1 := by omega
      rw [h5]
      calc k.choose a ≤ k ^ j * k.choose (a + j) := h1
        _ ≤ k ^ j * (k * k.choose (a + j + 1)) := h3
        _ = k ^ (j + 1) * k.choose (a + j + 1) := h4
  -- `t - 1 = k - (s + D)`, so `k.choose (t - 1) = k.choose (s + D)` by symmetry.
  have hSD : s + D ≤ k := by omega
  have ht1 : t - 1 = k - (s + D) := by omega
  have hsymm : k.choose (t - 1) = k.choose (s + D) := by
    rw [ht1]; exact Nat.choose_symm hSD
  -- Apply `general` with `j = D - 1`, `a = s`.
  have hgen : k.choose s ≤ k ^ (D - 1) * k.choose (s + (D - 1)) := general (D - 1) s (by omega)
  have hidx : s + (D - 1) = s + D - 1 := by omega
  rw [hidx] at hgen
  -- Combine: `k * k.choose s ≤ k ^ D * k.choose (s + D - 1)`.
  have key : k * k.choose s ≤ k ^ D * k.choose (s + D - 1) := by
    have hpow : k ^ D = k * k ^ (D - 1) := by
      conv_lhs => rw [← Nat.sub_add_cancel hD]
      rw [pow_succ']
    calc k * k.choose s
        ≤ k * (k ^ (D - 1) * k.choose (s + D - 1)) :=
          mul_le_mul (le_refl _) hgen (Nat.zero_le _) (Nat.zero_le _)
      _ = (k * k ^ (D - 1)) * k.choose (s + D - 1) := by ring
      _ = k ^ D * k.choose (s + D - 1) := by rw [hpow]
  -- Assemble the final inequality.
  rw [hsymm]
  calc k * (k.choose s * k.choose (s + D))
      = (k * k.choose s) * k.choose (s + D) := by ring
    _ ≤ (k ^ D * k.choose (s + D - 1)) * k.choose (s + D) :=
        mul_le_mul key (le_refl _) (Nat.zero_le _) (Nat.zero_le _)
    _ = k ^ D * (k.choose (s + D) * k.choose (s + D - 1)) := by ring

end StackedSVD.Edge

namespace StackedSVD.Edge

/-- L2e, the polynomial factors of the code fit in `(2 k) ^ (4 D)`: at `D = 1` and `k = 2`
the two sides are 160 and 512. -/
theorem code_pow_arith {k D : ℕ} (hk : 2 ≤ k) (hD : 1 ≤ D) :
    (4 * D + 1) * 2 ^ D * k ^ D * k ^ 3 ≤ ((2 * k) ^ 4) ^ D * k := by
  -- Step 1: `4 * D + 1 ≤ 8 ^ D` for every `D` (induction).
  have key : ∀ D : ℕ, 4 * D + 1 ≤ 8 ^ D := by
    intro D
    induction D with
    | zero => simp
    | succ n ih =>
      have h8 : (1 : ℕ) ≤ 8 ^ n := Nat.one_le_pow n 8 (by norm_num)
      have hstep : (8 : ℕ) ^ (n + 1) = 8 ^ n * 8 := pow_succ 8 n
      omega
  -- Step 2: `k ^ 3 ≤ k ^ (3 * D + 1)` since `k ≥ 1` and `3 ≤ 3 * D + 1` for `D ≥ 1`.
  have hk1 : 1 ≤ k := by omega
  have hpow3 : k ^ 3 ≤ k ^ (3 * D + 1) := by
    apply Nat.pow_le_pow_right hk1
    omega
  -- Step 3: combine, `(4 D + 1) * k ^ 3 ≤ 8 ^ D * k ^ (3 D + 1)`.
  have hcomb : (4 * D + 1) * k ^ 3 ≤ 8 ^ D * k ^ (3 * D + 1) :=
    Nat.mul_le_mul (key D) hpow3
  -- Step 4: multiply both sides by `2 ^ D * k ^ D`.
  have hfinal : (2 ^ D * k ^ D) * ((4 * D + 1) * k ^ 3)
      ≤ (2 ^ D * k ^ D) * (8 ^ D * k ^ (3 * D + 1)) :=
    Nat.mul_le_mul_left _ hcomb
  -- Step 5: the right side of `hfinal` equals the right side of the goal.
  have heq : (2 ^ D * k ^ D) * (8 ^ D * k ^ (3 * D + 1)) = ((2 * k) ^ 4) ^ D * k := by
    have e1 : (2 * k) ^ 4 = 16 * k ^ 4 := by ring
    have e3 : (16 : ℕ) = 2 * 8 := by norm_num
    rw [e1, e3, mul_pow, mul_pow, ← pow_mul,
      show (4 * D) = D + 3 * D from by ring, pow_add]
    ring
  calc (4 * D + 1) * 2 ^ D * k ^ D * k ^ 3
      = (2 ^ D * k ^ D) * ((4 * D + 1) * k ^ 3) := by ring
    _ ≤ (2 ^ D * k ^ D) * (8 ^ D * k ^ (3 * D + 1)) := hfinal
    _ = ((2 * k) ^ 4) ^ D * k := heq

end StackedSVD.Edge

/-! ## The assembly -/

namespace StackedSVD
namespace Edge

open Finset

/-- L2, the Furedi-Komlos count in one cell, the single combinatorial input of unit X. A walk
with `s` rows, `t` columns and no entry of multiplicity 1 costs at most `(2k)^12 / min(n, d)`
per unit of vertex deficiency `D = k + 1 - (s + t)` against the tree walks of the cell
`(s + D, t)`, which puts every extra vertex on the row side. The proof is the step encoding of
Furedi and Komlos 1981 page 237 (Anderson, Guionnet and Zeitouni, Lemma 2.1.23): of the `2k`
steps, `s + t - 1` are innovative and build a spanning tree, the forced steps follow the
Eulerian parity rule of the module docstring of `Count.lean`, and the bad steps are at most
`4 D` (`card_codeP_le` of `Code.lean`); a bad step is fixed by its position (at most `2k`) and
its endpoint (at most `k + 1`), so it costs `(2k)^2`, and the innovative-or-forced pattern is
part of the code.
Numerics (`scripts/stage3_count/cells.py`, 2026-09-10): 0 violations in 560 cells, `k <= 7`,
`(n, d)` from `(2k, 2k)` to `(10^12, 10^9)`; the exponent needed is at most `(2k)^1.37`,
against the `(2k)^12` available. -/
theorem cellCard_mul_le_pos {n d k s t : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
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
          (Nat.mul_le_mul (cellCard_le_code' (by omega) hs ht (by omega)) le_rfl)
    _ = (k * (k.choose s * k.choose (t - 1))) * ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (n.descFactorial s * min n d ^ D) := by ring
    _ ≤ (k * (k.choose s * k.choose (t - 1))) * ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (2 ^ D * n.descFactorial (s + D)) := by
        exact Nat.mul_le_mul le_rfl (descFac_mul_min_pow_le (by omega) hnd)
    _ ≤ (k ^ D * (k.choose (s + D) * k.choose (s + D - 1))) *
          ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          d.descFactorial t * (2 ^ D * n.descFactorial (s + D)) := by
        exact Nat.mul_le_mul (Nat.mul_le_mul (Nat.mul_le_mul
          (choose_mul_le_pow_choose hs ht hD1 (by omega)) le_rfl) le_rfl) le_rfl
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

/-- L2 with the degenerate cells discharged: a cell with no row or no column is empty. -/
private theorem cellCard_mul_le {n d k s t : ℕ} (hk : 1 ≤ k) (hst : s + t ≤ k)
    (hnd : 2 * k ≤ min n d) :
    cellCard n d k s t * (min n d) ^ (k + 1 - (s + t))
      ≤ ((2 * k) ^ 12) ^ (k + 1 - (s + t)) * cellCard n d k (s + (k + 1 - (s + t))) t := by
  rcases Nat.eq_zero_or_pos s with rfl | hs
  · rw [cellCard_eq_zero_of_zero hk (Or.inl rfl)]; simp
  · rcases Nat.eq_zero_or_pos t with rfl | ht
    · rw [cellCard_eq_zero_of_zero hk (Or.inr rfl)]; simp
    · exact cellCard_mul_le_pos (by omega) hs ht hst hnd

/-- X8, the Furedi-Komlos count in the `D`-step form, `D = k + 1 - v`: the class with `v`
vertices is at most `((2k)^12 / min(n, d))^D` times the tree class. Needs `1 ≤ k` (at `k = 0`
there is no tree walk) and `2k ≤ min(n, d)` (the assembly has `min(n, d) ≥ 2 (2k)^12` from
`fkFactor K n d k ≤ 1/2` with `1 ≤ K`). -/
theorem clsCard_mul_le {n d k v : ℕ} (hk : 1 ≤ k) (hv : v ≤ k + 1) (hnd : 2 * k ≤ min n d) :
    clsCard n d k v * (min n d) ^ (k + 1 - v)
      ≤ ((2 * k) ^ 12) ^ (k + 1 - v) * clsCard n d k (k + 1) := by
  classical
  rcases eq_or_lt_of_le hv with rfl | hlt
  · simp
  · have hvk : v ≤ k := Nat.lt_succ_iff.mp hlt
    have hstep : ∀ s ∈ range (v + 1),
        cellCard n d k s (v - s) * (min n d) ^ (k + 1 - v)
          ≤ ((2 * k) ^ 12) ^ (k + 1 - v) * cellCard n d k (s + (k + 1 - v)) (v - s) := by
      intro s hs
      have hsv : s ≤ v := Nat.lt_succ_iff.mp (Finset.mem_range.mp hs)
      have h := cellCard_mul_le (n := n) (d := d) (k := k) (s := s) (t := v - s) hk
        (by omega) hnd
      rwa [show s + (v - s) = v from by omega] at h
    have hreindex : ∀ s ∈ range (v + 1), ∀ s' ∈ range (v + 1),
        s + (k + 1 - v) = s' + (k + 1 - v) → s = s' := by
      intro s _ s' _ h
      omega
    have hsub : (range (v + 1)).image (fun s => s + (k + 1 - v)) ⊆ range (k + 2) := by
      intro z hz
      obtain ⟨s, hs, rfl⟩ := Finset.mem_image.mp hz
      rw [Finset.mem_range] at hs ⊢
      omega
    calc clsCard n d k v * (min n d) ^ (k + 1 - v)
        = ∑ s ∈ range (v + 1), cellCard n d k s (v - s) * (min n d) ^ (k + 1 - v) := by
          rw [clsCard_eq_sum_cellCard, Finset.sum_mul]
      _ ≤ ∑ s ∈ range (v + 1),
            ((2 * k) ^ 12) ^ (k + 1 - v) * cellCard n d k (s + (k + 1 - v)) (v - s) :=
          Finset.sum_le_sum hstep
      _ = ((2 * k) ^ 12) ^ (k + 1 - v) * ∑ s ∈ range (v + 1),
            cellCard n d k (s + (k + 1 - v)) (k + 1 - (s + (k + 1 - v))) := by
          rw [← Finset.mul_sum]
          refine congrArg _ (Finset.sum_congr rfl fun s hs => ?_)
          have hsv : s ≤ v := Nat.lt_succ_iff.mp (Finset.mem_range.mp hs)
          rw [show k + 1 - (s + (k + 1 - v)) = v - s from by omega]
      _ = ((2 * k) ^ 12) ^ (k + 1 - v) * ∑ z ∈ (range (v + 1)).image
            (fun s => s + (k + 1 - v)), cellCard n d k z (k + 1 - z) := by
          rw [Finset.sum_image hreindex]
      _ ≤ ((2 * k) ^ 12) ^ (k + 1 - v) * ∑ z ∈ range (k + 2), cellCard n d k z (k + 1 - z) :=
          Nat.mul_le_mul (le_refl _) (Finset.sum_le_sum_of_subset hsub)
      _ = ((2 * k) ^ 12) ^ (k + 1 - v) * clsCard n d k (k + 1) := by
          rw [clsCard_eq_sum_cellCard]

end Edge
end StackedSVD
