/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R6het
import StackedSVD.RankR.Het.Bulk
import StackedSVD.LinAlg.Frame

/-!
# Deterministic and scalar lemmas for the one-detectable regime (F18b, units U4a, U4b, U4e.0, U4e.1)

A positive secular value caps the top eigenvalue of a rank-one update
(`lamMax_add_vecMulVec_le_of_secular_pos`); the second sorted eigenvalue of `W₀ + Q Qᵀ`
with two columns is at most the top eigenvalue of `W₀` plus one column
(`eigenvalues₀_le_of_two_cols`); a quadratic-form lower bound from entrywise closeness
(`form_le_of_close_to_form`); and `Ψ'(z) → +∞` at the edge (`PsihetDeriv_tendsto_atTop`).
Plan: `notes/archive/F18b_plan.md`, section 2 (U4) and section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace R6het

/-- U4a: a positive secular value at `z₀` caps the top eigenvalue of the rank-one update. -/
theorem lamMax_add_vecMulVec_le_of_secular_pos {p : ℕ} (hp : 0 < p)
    {W₀ : Matrix (Fin p) (Fin p) ℝ} (hW₀ : W₀.IsHermitian) (q : Fin p → ℝ)
    (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) {z₀ : ℝ}
    (h₀ : lamMax W₀ hW₀ < z₀) (hs : 0 < R4.secular W₀ q z₀) :
    lamMax (W₀ + Matrix.vecMulVec q q) hA ≤ z₀ := by
  rcases eq_or_ne q 0 with hq0 | hq
  · subst hq0
    have hWeq : W₀ + Matrix.vecMulVec (0 : Fin p → ℝ) 0 = W₀ := by
      ext i j
      simp [Matrix.add_apply]
    rw [lamMax_congr hWeq hA hW₀]
    exact h₀.le
  · by_contra hcon
    rw [not_le] at hcon
    have hlt : lamMax W₀ hW₀ < lamMax (W₀ + Matrix.vecMulVec q q) hA := lt_trans h₀ hcon
    obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hp
    have hspec : lamMax (W₀ + Matrix.vecMulVec q q) hA
        ∈ spectrum ℝ (toOp (W₀ + Matrix.vecMulVec q q)) := by
      rw [R4.spectrum_toOp]
      exact hj ▸ hA.eigenvalues_mem_spectrum_real j
    have hzero : R4.secular W₀ q (lamMax (W₀ + Matrix.vecMulVec q q) hA) = 0 :=
      (R4.secular_eq_zero_iff hW₀ hlt).mpr hspec
    have hmono := R4.secular_strictMonoOn hW₀ hq (Set.mem_Ioi.2 h₀) (Set.mem_Ioi.2 hlt) hcon
    linarith

end R6het

namespace Frame

/-- U4b.0: the Gram of a two-column matrix is the sum of the two rank-one terms. -/
theorem mul_transpose_two_cols {p : ℕ} (Q : Matrix (Fin p) (Fin 2) ℝ)
    {l₀ l₁ : Fin 2} (hne : l₀ ≠ l₁) :
    Q * Qᵀ = Matrix.vecMulVec (fun i => Q i l₀) (fun i => Q i l₀)
      + Matrix.vecMulVec (fun i => Q i l₁) (fun i => Q i l₁) := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.add_apply,
    Matrix.vecMulVec_apply, Fin.sum_univ_two]
  fin_cases l₀ <;> fin_cases l₁ <;> first
    | exact absurd rfl hne
    | (simp <;> ring)

/-- U4b: with one of the two signal columns absorbed into the bulk, every sorted eigenvalue
of index at least one sits below `τ`. -/
theorem eigenvalues₀_le_of_two_cols {p : ℕ} {W₀ S : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin 2) ℝ} (hS : S.IsHermitian)
    (hSeq : S = W₀ + Q * Qᵀ) {l₀ l₁ : Fin 2} (hne : l₀ ≠ l₁)
    (hW' : (W₀ + Matrix.vecMulVec (fun i => Q i l₁) (fun i => Q i l₁)).IsHermitian)
    {τ : ℝ} (hτ : lamMax (W₀ + Matrix.vecMulVec (fun i => Q i l₁) (fun i => Q i l₁)) hW' ≤ τ)
    (k : Fin (Fintype.card (Fin p))) (hk : 1 ≤ (k : ℕ)) : hS.eigenvalues₀ k ≤ τ := by
  set Q' : Matrix (Fin p) (Fin 1) ℝ := Matrix.of fun i (_ : Fin 1) => Q i l₀ with hQ'def
  have hQ'eq : Matrix.vecMulVec (fun i => Q i l₀) (fun i => Q i l₀) = Q' * Q'ᵀ := by
    ext i j
    simp [Matrix.mul_apply, Matrix.transpose_apply, Matrix.vecMulVec_apply, hQ'def,
      Matrix.of_apply]
  have hSeq' : S = (W₀ + Matrix.vecMulVec (fun i => Q i l₁) (fun i => Q i l₁)) + Q' * Q'ᵀ := by
    rw [hSeq, mul_transpose_two_cols Q hne, hQ'eq]
    abel
  exact Frame.eigenvalues₀_le_of_split hW' hS hSeq' hτ k hk

end Frame

/-! ### U4e helper: a form bound from closeness to a general positive matrix -/

namespace HetBulkDet

/-- U4e.0: `μ₀ - η rr` lower bound of the form from entrywise closeness to a matrix whose own
form is at least `μ₀`. The non-diagonal twin of `form_le_of_close_to_diag`. -/
theorem form_le_of_close_to_form {rr : ℕ} {G A : Matrix (Fin rr) (Fin rr) ℝ}
    {η μ₀ : ℝ} (hη : 0 ≤ η)
    (hA : ∀ y : Fin rr → ℝ, μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ (A *ᵥ y))
    (hclose : ∀ k l, |G k l - A k l| ≤ η) (y : Fin rr → ℝ) :
    (μ₀ - η * rr) * (y ⬝ᵥ y) ≤ y ⬝ᵥ (G *ᵥ y) := by
  have hE : ∀ k l, |(G - A) k l| ≤ η := fun k l => by
    rw [Matrix.sub_apply]; exact hclose k l
  have hbound : |y ⬝ᵥ ((G - A) *ᵥ y)| ≤ (rr : ℝ) * η * ∑ k, y k ^ 2 :=
    Frame.abs_dotProduct_mulVec_le hη hE y
  have hsq : y ⬝ᵥ y = ∑ k, y k ^ 2 := by
    unfold dotProduct
    exact Finset.sum_congr rfl fun k _ => by ring
  have hGA : y ⬝ᵥ ((G - A) *ᵥ y) = y ⬝ᵥ (G *ᵥ y) - y ⬝ᵥ (A *ᵥ y) := by
    rw [Matrix.sub_mulVec, dotProduct_sub]
  rw [hGA, ← hsq] at hbound
  have h1 := (abs_le.mp hbound).1
  have h2 := hA y
  have hfinal : (μ₀ - η * rr) * (y ⬝ᵥ y) = μ₀ * (y ⬝ᵥ y) - (rr : ℝ) * η * (y ⬝ᵥ y) := by ring
  rw [hfinal]
  linarith [h1, h2]

end HetBulkDet

namespace MPhet

variable {M : ℕ}

/-- U4e.1: `Ψ'(z) → +∞` at the edge, the `θ = 0` case of `FhetDeriv_tendsto_atTop`. -/
theorem PsihetDeriv_tendsto_atTop {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) :
    Tendsto (PsihetDeriv c w) (𝓝[>] bHet c w) atTop := by
  have heq : FhetDeriv (0 : Fin M → ℝ) c w = PsihetDeriv c w := by
    funext z
    have hP : PhihetDeriv (0 : Fin M → ℝ) c w z = 0 := by
      unfold PhihetDeriv
      simp
    change PhihetDeriv (0 : Fin M → ℝ) c w z + PsihetDeriv c w z = PsihetDeriv c w z
    rw [hP, zero_add]
  rw [← heq]
  exact FhetDeriv_tendsto_atTop hc hw

end MPhet

end StackedSVD
