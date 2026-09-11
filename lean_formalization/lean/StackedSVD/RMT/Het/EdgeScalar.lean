/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.MPhet

/-!
# Campaign E, task B: the scalar bridge from the edge objective to the Silverstein API

Plan: `notes/archive/heteroedge_sharp.md`, task B. The sharp Sudakov-Fernique bound of task C
produces, for every `γ > max_i w_i²`, the finite-`N` bound
`E‖D G/√d‖² ≤ edgeObjective c_N w γ` with

  `edgeObjective c w γ = γ (1 + ∑_i c_i w_i²/(γ - w_i²))`.

This file proves the two facts the assembly needs.

* Step E6: the substitution `s = -1/γ` turns `edgeObjective` into `zfun` of
  `RMT/Het/MPhet.lean`, so the existing critical point `sStar` evaluates the objective at
  the bulk edge: `edgeObjective c w γ⋆ = bHet c w` with `γ⋆ = -1/sStar c w`. The point
  `γ⋆` is admissible: `γ⋆ > wSqMax w ≥ w_i²` for every `i`, so no denominator vanishes.
* Step E7: at the fixed `γ⋆` the objective is affine in the aspect vector `c`, so the
  finite-`N` aspect ratios `n_i/d → c_i` carry `edgeObjective c_N w γ⋆ → bHet c w`.

No new analysis enters here; every statement is algebra over `MPhet`.
-/

open Filter Topology Finset

namespace StackedSVD
namespace MPhet

open Scalars

variable {M : ℕ}

/-! ### The edge objective and its `zfun` form (step E6) -/

/-- The edge objective `γ (1 + ∑_i c_i w_i²/(γ - w_i²))` of the sharp Sudakov-Fernique
bound. It is `zfun c w (-1/γ)` for admissible `γ` (see `edgeObjective_eq_zfun`). -/
noncomputable def edgeObjective (c w : Fin M → ℝ) (γ : ℝ) : ℝ :=
  γ * (1 + ∑ i, c i * w i ^ 2 / (γ - w i ^ 2))

/-- The substitution `s = -1/γ`: the edge objective is the Silverstein map `zfun`. -/
theorem edgeObjective_eq_zfun {c w : Fin M → ℝ} {γ : ℝ} (hγ : 0 < γ)
    (hw : ∀ i, w i ^ 2 < γ) : edgeObjective c w γ = zfun c w (-1 / γ) := by
  have hγ0 : γ ≠ 0 := hγ.ne'
  unfold edgeObjective zfun
  rw [mul_add, mul_one, Finset.mul_sum]
  congr 1
  · field_simp
  · refine Finset.sum_congr rfl fun i _ => ?_
    have hd : γ - w i ^ 2 ≠ 0 := sub_ne_zero.mpr (hw i).ne'
    have h2 : 1 + w i ^ 2 * (-1 / γ) = (γ - w i ^ 2) / γ := by field_simp; ring
    rw [h2, div_div_eq_mul_div]
    field_simp

/-- The critical point is negative, so `γ⋆ = -1/s⋆` is positive. -/
theorem neg_inv_sStar_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    0 < -1 / sStar c w :=
  div_pos_of_neg_of_neg (by norm_num) (sStar_mem hc hw).2

/-- `γ⋆ = -1/s⋆ > max_i w_i²`: the objective is evaluated strictly inside its domain.
This mirrors `wSqMax_lt_gamHet` on the physical branch. -/
theorem neg_inv_sStar_gt_wSqMax {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    wSqMax w < -1 / sStar c w := by
  have hW := wSqMax_pos hw
  have hst := sStar_mem hc hw
  have h1 : -1 < sStar c w * wSqMax w := by
    have hlo := hst.1
    unfold sLo at hlo
    exact (div_lt_iff₀ hW).mp hlo
  rw [lt_div_iff_of_neg hst.2]
  linarith

/-- Every block variance is below `γ⋆`, so no denominator of the objective vanishes. -/
theorem sq_lt_neg_inv_sStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (i : Fin M) : w i ^ 2 < -1 / sStar c w :=
  lt_of_le_of_lt (le_wSqMax w i) (neg_inv_sStar_gt_wSqMax hc hw)

/-- Step E6: at the critical point the edge objective is exactly the bulk edge. -/
theorem edgeObjective_sStar_eq_bHet {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) : edgeObjective c w (-1 / sStar c w) = bHet c w := by
  have hs0 : sStar c w ≠ 0 := (sStar_mem hc hw).2.ne
  rw [edgeObjective_eq_zfun (neg_inv_sStar_pos hc hw) (sq_lt_neg_inv_sStar hc hw)]
  unfold bHet
  congr 1
  field_simp

/-! ### Aspect passage at the fixed critical point (step E7) -/

/-- Step E7: at the fixed `γ⋆ = -1/s⋆` the objective is affine in the aspect vector, so the
finite-`N` aspect ratios `n_i/d → c_i` carry the objective to the bulk edge. -/
theorem edgeObjective_tendsto {c w : Fin M → ℝ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hn : ∀ i, Tendsto (fun N => (n i N : ℝ) / d N) atTop (nhds (c i))) :
    Tendsto (fun N => edgeObjective (fun i => (n i N : ℝ) / d N) w (-1 / sStar c w))
      atTop (nhds (bHet c w)) := by
  rw [← edgeObjective_sStar_eq_bHet hc hw]
  unfold edgeObjective
  exact ((tendsto_finsetSum Finset.univ
    fun i _ => ((hn i).mul_const (w i ^ 2)).div_const
      (-1 / sStar c w - w i ^ 2)).const_add 1).const_mul _

/-! ### Sanity: the objective is positive on its domain (step E10) -/

/-- Every summand is nonnegative on the admissible range, so the objective is at least `γ`.
A cheap regression check on the sign conventions. -/
theorem lt_edgeObjective {c w : Fin M → ℝ} {γ : ℝ} (hγ : 0 < γ) (hc : ∀ i, 0 ≤ c i)
    (hw : ∀ i, w i ^ 2 < γ) : γ ≤ edgeObjective c w γ := by
  have hsum : 0 ≤ ∑ i, c i * w i ^ 2 / (γ - w i ^ 2) :=
    Finset.sum_nonneg fun i _ =>
      div_nonneg (mul_nonneg (hc i) (sq_nonneg _)) (sub_pos.mpr (hw i)).le
  unfold edgeObjective
  nlinarith

end MPhet
end StackedSVD
