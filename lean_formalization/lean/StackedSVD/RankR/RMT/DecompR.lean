/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs
import StackedSVD.Prob.GaussianMatrix

/-!
# Task U7b: the rank-`r` orthogonal decomposition

Rank-`r` twin of `StackedSVD.RMT.R6` section "Decomposition" (`RMT/R6.lean:94-175`, the block
`aOf`, `rvecOf`, `aOf_mul`, `dotProduct_rvecOf`, `eq_smul_add_rvecOf`,
`dotProduct_rvecOf_self_le`, `aOf_sq_mul_le`, `toLp_eq_smul_add`, `norm_toLp_rvecOf_le`) and of
`R6.SpikedModel.measurable_rvecOf` (`RMT/R6.lean:571`). Deliverable (part 1) of task U7 of
`notes/archive/rankr_plan_A.md`.

The rank-1 coefficient `aOf v q = (v ⬝ᵥ q) / (q ⬝ᵥ q)` becomes the least-squares coefficient
`aOfR Q v = (Qᵀ Q)⁻¹ (Qᵀ v)`, so `rvecOfR Q v = v - Q (aOfR Q v)` is the residual of the
orthogonal projection of `v` onto the column space of `Q`. Every rank-1 lemma has a direct
twin, proved from `Matrix.mulVec_mulVec`, `Matrix.mul_nonsing_inv` / `Matrix.nonsing_inv_mul`,
`Matrix.dotProduct_mulVec` and `Matrix.mulVec_transpose`, under the extra hypothesis
`IsUnit (Qᵀ * Q).det`. The rank-1 file needs no such hypothesis because its `aOf` is a
division, whose junk value at `q = 0` still satisfies the defining identity `aOf v q * (q⬝ᵥq)
= v⬝ᵥq` (`R6.aOf_mul`); `Matrix.inv` of a singular matrix is the junk value `0` instead, which
does **not** satisfy the analogous identity `(Qᵀ*Q)⁻¹ * (Qᵀ*Q) = 1`, so every orthogonality
lemma here needs `IsUnit (Qᵀ * Q).det` to rule that case out.

Imports: `StackedSVD.Defs` (`import Mathlib` in full, for the Matrix and measurability API) and
`StackedSVD.Prob.GaussianMatrix` (for `inner_euclidean_eq_dotProduct`, the bridge between the
`EuclideanSpace` norm and `⬝ᵥ`). Both are already-compiled, shared infrastructure; nothing
from `RMT/` is imported.

Namespace: `StackedSVD.DecompR`.
-/

open scoped Matrix

namespace StackedSVD
namespace DecompR

variable {d r : ℕ}

/-! ### 1. The orthogonal decomposition of `v` along the columns of `Q` -/

/-- The least-squares coefficient of `Q`'s columns in `v`: `(Qᵀ Q)⁻¹ (Qᵀ v)`. Junk value `0`
when `Qᵀ Q` is singular. -/
noncomputable def aOfR (Q : Matrix (Fin d) (Fin r) ℝ) (v : Fin d → ℝ) : Fin r → ℝ :=
  (Qᵀ * Q)⁻¹ *ᵥ (Qᵀ *ᵥ v)

/-- The residual of `v` after removing its projection onto the columns of `Q`. -/
noncomputable def rvecOfR (Q : Matrix (Fin d) (Fin r) ℝ) (v : Fin d → ℝ) : Fin d → ℝ :=
  v - Q *ᵥ aOfR Q v

theorem eq_mulVec_add_rvecOfR (Q : Matrix (Fin d) (Fin r) ℝ) (v : Fin d → ℝ) :
    v = Q *ᵥ aOfR Q v + rvecOfR Q v := by
  rw [rvecOfR]; abel

theorem transpose_mulVec_rvecOfR (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) : Qᵀ *ᵥ rvecOfR Q v = 0 := by
  have haOfR : aOfR Q v = (Qᵀ * Q)⁻¹ *ᵥ (Qᵀ *ᵥ v) := rfl
  have e1 : Qᵀ *ᵥ (Q *ᵥ aOfR Q v) = (Qᵀ * Q) *ᵥ aOfR Q v := Matrix.mulVec_mulVec _ _ _
  have e3 : (Qᵀ * Q) *ᵥ ((Qᵀ * Q)⁻¹ *ᵥ (Qᵀ *ᵥ v)) = ((Qᵀ * Q) * (Qᵀ * Q)⁻¹) *ᵥ (Qᵀ *ᵥ v) :=
    Matrix.mulVec_mulVec _ _ _
  have e4 : (Qᵀ * Q) * (Qᵀ * Q)⁻¹ = 1 := Matrix.mul_nonsing_inv _ hQ
  have e5 : Qᵀ *ᵥ (Q *ᵥ aOfR Q v) = Qᵀ *ᵥ v := by
    rw [e1, haOfR, e3, e4, Matrix.one_mulVec]
  change Qᵀ *ᵥ (v - Q *ᵥ aOfR Q v) = 0
  rw [Matrix.mulVec_sub, e5, sub_self]

theorem dotProduct_mulVec_rvecOfR (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) (a : Fin r → ℝ) : (Q *ᵥ a) ⬝ᵥ rvecOfR Q v = 0 := by
  have e1 : (Q *ᵥ a) ⬝ᵥ rvecOfR Q v = rvecOfR Q v ⬝ᵥ (Q *ᵥ a) := dotProduct_comm _ _
  have e2 : rvecOfR Q v ⬝ᵥ (Q *ᵥ a) = (rvecOfR Q v ᵥ* Q) ⬝ᵥ a := Matrix.dotProduct_mulVec _ _ _
  have e3 : rvecOfR Q v ᵥ* Q = Qᵀ *ᵥ rvecOfR Q v := (Matrix.mulVec_transpose Q (rvecOfR Q v)).symm
  have e4 : Qᵀ *ᵥ rvecOfR Q v = 0 := transpose_mulVec_rvecOfR Q hQ
  rw [e1, e2, e3, e4, zero_dotProduct]

/-- The Pythagorean identity behind both `dotProduct_rvecOfR_self_le` and `aOfR_sq_mul_le`. -/
theorem dotProduct_self_eq_add_dotProduct_self (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) :
    v ⬝ᵥ v = (Q *ᵥ aOfR Q v) ⬝ᵥ (Q *ᵥ aOfR Q v) + rvecOfR Q v ⬝ᵥ rvecOfR Q v := by
  have horth : (Q *ᵥ aOfR Q v) ⬝ᵥ rvecOfR Q v = 0 := dotProduct_mulVec_rvecOfR Q hQ (aOfR Q v)
  conv_lhs => rw [eq_mulVec_add_rvecOfR Q v]
  rw [add_dotProduct, dotProduct_add, dotProduct_add, horth,
      dotProduct_comm (rvecOfR Q v) (Q *ᵥ aOfR Q v), horth]
  ring

/-- Nonnegativity of the dot product with itself. -/
theorem dotProduct_self_nonneg' (v : Fin d → ℝ) : 0 ≤ v ⬝ᵥ v :=
  Finset.sum_nonneg fun _ _ => mul_self_nonneg _

theorem dotProduct_rvecOfR_self_le (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) : rvecOfR Q v ⬝ᵥ rvecOfR Q v ≤ v ⬝ᵥ v := by
  have h := dotProduct_self_eq_add_dotProduct_self Q hQ (v := v)
  have hnn : 0 ≤ (Q *ᵥ aOfR Q v) ⬝ᵥ (Q *ᵥ aOfR Q v) := dotProduct_self_nonneg' _
  linarith

theorem norm_toLp_rvecOfR_le (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) :
    ‖(WithLp.toLp 2 (rvecOfR Q v) : EuclideanSpace ℝ (Fin d))‖
      ≤ ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))‖ := by
  have h1 : ‖(WithLp.toLp 2 (rvecOfR Q v) : EuclideanSpace ℝ (Fin d))‖ ^ 2
      = rvecOfR Q v ⬝ᵥ rvecOfR Q v := by
    rw [← real_inner_self_eq_norm_sq]; exact inner_euclidean_eq_dotProduct _ _
  have h2 : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))‖ ^ 2 = v ⬝ᵥ v := by
    rw [← real_inner_self_eq_norm_sq]; exact inner_euclidean_eq_dotProduct _ _
  have h3 : rvecOfR Q v ⬝ᵥ rvecOfR Q v ≤ v ⬝ᵥ v := dotProduct_rvecOfR_self_le Q hQ
  rw [← Real.sqrt_sq (norm_nonneg (WithLp.toLp 2 (rvecOfR Q v) : EuclideanSpace ℝ (Fin d))),
      ← Real.sqrt_sq (norm_nonneg (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))), h1, h2]
  exact Real.sqrt_le_sqrt h3

theorem aOfR_sq_mul_le (Q : Matrix (Fin d) (Fin r) ℝ) {v : Fin d → ℝ}
    (hQ : IsUnit (Qᵀ * Q).det) {μ : ℝ}
    (hμ : ∀ y : Fin r → ℝ, μ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y)) :
    μ * (aOfR Q v ⬝ᵥ aOfR Q v) ≤ v ⬝ᵥ v := by
  have hQa : (Q *ᵥ aOfR Q v) ⬝ᵥ (Q *ᵥ aOfR Q v) = aOfR Q v ⬝ᵥ ((Qᵀ * Q) *ᵥ aOfR Q v) := by
    have e1 : (Q *ᵥ aOfR Q v) ⬝ᵥ (Q *ᵥ aOfR Q v)
        = ((Q *ᵥ aOfR Q v) ᵥ* Q) ⬝ᵥ aOfR Q v := Matrix.dotProduct_mulVec _ _ _
    have e2 : (Q *ᵥ aOfR Q v) ᵥ* Q = Qᵀ *ᵥ (Q *ᵥ aOfR Q v) :=
      (Matrix.mulVec_transpose Q (Q *ᵥ aOfR Q v)).symm
    have e3 : Qᵀ *ᵥ (Q *ᵥ aOfR Q v) = (Qᵀ * Q) *ᵥ aOfR Q v := Matrix.mulVec_mulVec _ _ _
    rw [e1, e2, e3, dotProduct_comm]
  have hbound : μ * (aOfR Q v ⬝ᵥ aOfR Q v) ≤ (Q *ᵥ aOfR Q v) ⬝ᵥ (Q *ᵥ aOfR Q v) := by
    rw [hQa]; exact hμ (aOfR Q v)
  have hpyth := dotProduct_self_eq_add_dotProduct_self Q hQ (v := v)
  have hnn : 0 ≤ rvecOfR Q v ⬝ᵥ rvecOfR Q v := dotProduct_self_nonneg' _
  linarith

theorem toLp_eq_mulVec_add (Q : Matrix (Fin d) (Fin r) ℝ) (v : Fin d → ℝ) :
    (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))
      = WithLp.toLp 2 (Q *ᵥ aOfR Q v) + WithLp.toLp 2 (rvecOfR Q v) := by
  apply WithLp.ofLp_injective
  change v = Q *ᵥ aOfR Q v + rvecOfR Q v
  exact eq_mulVec_add_rvecOfR Q v

/-! ### 2. Measurability -/

section Measurability

variable {α : Type*} [MeasurableSpace α]

/-- The entry `Q i j` of a matrix is measurable in `Q`. -/
private theorem measurable_matrix_apply {m n : ℕ} (i : Fin m) (j : Fin n) :
    Measurable fun Q : Matrix (Fin m) (Fin n) ℝ => Q i j :=
  (measurable_pi_apply j).comp (measurable_pi_apply i)

/-- Real-valued twin of `R2.measurable_det` (`RMT/R2.lean:205`, stated there for `ℂ`). -/
private theorem measurable_det_real {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℝ}
    (hM : ∀ i j, Measurable fun a => M a i j) : Measurable fun a => (M a).det := by
  simp only [Matrix.det_apply']
  refine Finset.measurable_sum _ fun σ _ => ?_
  exact measurable_const.mul (Finset.measurable_prod _ fun i _ => hM (σ i) i)

/-- Real-valued twin of `R2.measurable_adjugate` (`RMT/R2.lean:211`). -/
private theorem measurable_adjugate_real {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℝ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin n) :
    Measurable fun a => (M a).adjugate i j := by
  simp only [Matrix.adjugate_apply]
  refine measurable_det_real fun k l => ?_
  by_cases h : k = j
  · subst h
    simp only [Matrix.updateRow_self]
    exact measurable_const
  · simp only [Matrix.updateRow_ne h]
    exact hM k l

/-- Real-valued twin of `R2.measurable_inv_entry` (`RMT/R2.lean:223`): every entry of the
nonsingular inverse is a measurable (in fact polynomial) function of the matrix entries, with
no invertibility hypothesis, because `Matrix.inv` at a singular matrix is the junk value `0`. -/
private theorem measurable_inv_entry_real {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℝ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin n) :
    Measurable fun a => (M a)⁻¹ i j := by
  simp only [Matrix.inv_def, Matrix.smul_apply, smul_eq_mul, Ring.inverse_eq_inv']
  exact ((measurable_det_real hM).inv).mul (measurable_adjugate_real hM i j)

theorem measurable_transpose_mul_self_apply (k l : Fin r) :
    Measurable fun Q : Matrix (Fin d) (Fin r) ℝ => (Qᵀ * Q) k l := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply]
  exact Finset.measurable_sum _ fun i _ =>
    (measurable_matrix_apply i k).mul (measurable_matrix_apply i l)

theorem measurable_inv_transpose_mul_self_apply (k l : Fin r) :
    Measurable fun Q : Matrix (Fin d) (Fin r) ℝ => (Qᵀ * Q)⁻¹ k l :=
  measurable_inv_entry_real (fun i j => measurable_transpose_mul_self_apply i j) k l

theorem measurable_transpose_mulVec_apply (i : Fin r) :
    Measurable fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) => (q.1ᵀ *ᵥ q.2) i := by
  simp only [Matrix.mulVec_apply_eq_sum, Matrix.transpose_apply]
  exact Finset.measurable_sum _ fun k _ =>
    ((measurable_matrix_apply k i).comp measurable_fst).mul
      ((measurable_pi_apply k).comp measurable_snd)

theorem measurable_aOfR :
    Measurable (fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) => aOfR q.1 q.2) := by
  refine measurable_pi_lambda _ fun j => ?_
  simp only [aOfR, Matrix.mulVec_apply_eq_sum]
  exact Finset.measurable_sum _ fun i _ =>
    ((measurable_inv_transpose_mul_self_apply j i).comp measurable_fst).mul
      (measurable_transpose_mulVec_apply i)

theorem measurable_rvecOfR :
    Measurable (fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) => rvecOfR q.1 q.2) := by
  refine measurable_pi_lambda _ fun j => ?_
  change Measurable fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) =>
    q.2 j - (q.1 *ᵥ aOfR q.1 q.2) j
  have h1 : Measurable fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) => q.2 j :=
    (measurable_pi_apply j).comp measurable_snd
  have h2 : Measurable fun q : Matrix (Fin d) (Fin r) ℝ × (Fin d → ℝ) =>
      (q.1 *ᵥ aOfR q.1 q.2) j := by
    simp only [Matrix.mulVec_apply_eq_sum]
    exact Finset.measurable_sum _ fun k _ =>
      ((measurable_matrix_apply j k).comp measurable_fst).mul
        ((measurable_pi_apply k).comp measurable_aOfR)
  exact h1.sub h2

end Measurability

end DecompR
end StackedSVD
