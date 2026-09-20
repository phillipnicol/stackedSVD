/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Scalars

/-!
# The rank-one reduction of the single-weight scalar layer

Unit S3 of `notes/archive/singleweight_plan.md` section 4.1. At `r = 1`, `r_i = 1` and `R_i = 1` the
matrix secular equation of `main_paper.tex:2115` is a `1 × 1` determinant, so it is the scalar
secular equation of `lem:secular_equation`, and the limit of `main_paper.tex:2119` is the
paper's `L(w)` (`thm:stacksvd_weighted`). The two theorems here tie the new layer to the proved
rank-one tree of `StackSVDWeighted.lean`.

Numeric check: the two sides agree to 3.6e-15 over the scan of `notes/archive/singleweight_plan.md`
section 2.1.

STATUS 2026-09-05: statements only. Neither proof is written yet, both open sites have a row
in `docs/SORRIES.md`, and the statements wait for the user (`CLAUDE.md` workflow step 3).
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

/-- At `r = 1`, `r_i = 1` and `R_i = 1`, `sigMat` collapses to a `1 × 1` diagonal matrix.
Shared private helper for the two theorems below. -/
private theorem sigMat_one {M : ℕ} (θ : Fin M → ℝ) (i : Fin M) :
    sigMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) i
      = Matrix.diagonal (fun _ : Fin 1 => θ i ^ 2) := by
  unfold sigMat
  rw [Matrix.transpose_one, Matrix.one_mul, Matrix.mul_one]

/-- At `r = 1`, `r_i = 1` and `R_i = 1` the matrix secular equation of `main_paper.tex:2115`
is the scalar one of `lem:secular_equation`: `det(I₁ - secMat) = 1 - ∑_i θ_i² w_i²/(γ - w_i²)`
is `Scalars.secular θ w γ`. So a root of one is a root of the other, above `max_i w_i²`. -/
theorem isSecularRoot_one_iff {M : ℕ} (θ w : Fin M → ℝ) (γ : ℝ) :
    IsSecularRoot (r := 1) (rk := fun _ => 1) (fun i _ => θ i)
        (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ
      ↔ Scalars.IsGammaTop θ w γ := by
  have hdet0 : (1 - secMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ) 0 0
      = 1 - ∑ i : Fin M, w i ^ 2 / (γ - w i ^ 2) * θ i ^ 2 := by
    unfold secMat
    simp [Matrix.sub_apply, Matrix.one_apply_eq, Matrix.sum_apply, Matrix.smul_apply,
      sigMat_one, Matrix.diagonal_apply_eq, smul_eq_mul]
  have hterm : ∀ i : Fin M,
      w i ^ 2 / (γ - w i ^ 2) * θ i ^ 2 + θ i ^ 2 * w i ^ 2 / (w i ^ 2 - γ) = 0 := by
    intro i
    rw [show γ - w i ^ 2 = -(w i ^ 2 - γ) by ring, div_neg]
    ring
  have hsplit : (∑ i : Fin M, w i ^ 2 / (γ - w i ^ 2) * θ i ^ 2)
      + ∑ i : Fin M, θ i ^ 2 * w i ^ 2 / (w i ^ 2 - γ) = 0 := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_eq_zero (fun i _ => hterm i)
  have hdet : (1 - secMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ).det
      = Scalars.secular θ w γ := by
    rw [Matrix.det_fin_one, hdet0]
    unfold Scalars.secular
    linarith [hsplit]
  unfold IsSecularRoot Scalars.IsGammaTop
  rw [hdet]

/-- At `r = 1`, `r_i = 1` and `R_i = 1` the limit of `main_paper.tex:2119` is the paper's
`L(w)` of `thm:stacksvd_weighted` (`Scalars.Lw`, `StackSVDWeighted.lean:386`).

`h4 : Scalars.Assumption4 θ c w` is needed because `Scalars.Lw` is `0` below the
detectability threshold `eq:assumption4` while `swLimit` has no such branch. That shape
difference is modeling choice 4 of `notes/archive/singleweight_plan.md` section 4.3. -/
theorem swLimit_one {M : ℕ} (θ c w : Fin M → ℝ) {γ : ℝ}
    (z : Fin 1 → EuclideanSpace ℝ (Fin 1)) (hγ : Scalars.IsGammaTop θ w γ)
    (h4 : Scalars.Assumption4 θ c w) (hz : ‖z 0‖ = 1) :
    swLimit (r := 1) (rk := fun _ => 1) (fun i _ => θ i)
        (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w c (fun _ => γ) z
      = Scalars.Lw θ c w := by
  have hgt : Scalars.gammaTop θ w = γ := Scalars.gammaTop_eq hγ
  have hzsq : (WithLp.ofLp (z 0) 0) ^ 2 = 1 := by
    have hns : ‖z 0‖ ^ 2 = ∑ i : Fin 1, (WithLp.ofLp (z 0) i) ^ 2 :=
      EuclideanSpace.real_norm_sq_eq (z 0)
    rw [Fin.sum_univ_one, hz, one_pow] at hns
    exact hns.symm
  have hderiv0 : secDerivMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ 0 0
      = ∑ i : Fin M, θ i ^ 2 * w i ^ 2 / (γ - w i ^ 2) ^ 2 := by
    unfold secDerivMat
    simp only [Matrix.sum_apply, Matrix.smul_apply, sigMat_one, Matrix.diagonal_apply_eq,
      smul_eq_mul]
    apply Finset.sum_congr rfl
    intro i _
    ring
  have hquad : WithLp.ofLp (z 0) ⬝ᵥ
      (secDerivMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ *ᵥ
        WithLp.ofLp (z 0))
      = ∑ i : Fin M, θ i ^ 2 * w i ^ 2 / (γ - w i ^ 2) ^ 2 := by
    simp only [dotProduct, Matrix.mulVec, Fin.sum_univ_one]
    rw [show WithLp.ofLp (z 0) 0 *
        (secDerivMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ 0 0
          * WithLp.ofLp (z 0) 0)
        = secDerivMat (fun i _ => θ i) (fun _ => (1 : Matrix (Fin 1) (Fin 1) ℝ)) w γ 0 0
          * (WithLp.ofLp (z 0) 0) ^ 2 from by ring,
      hzsq, mul_one, hderiv0]
  have heta : Scalars.eta1 θ c w = 1 - ∑ i : Fin M, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by
    unfold Scalars.eta1
    rw [hgt]
  have hden : Scalars.LwDen θ w = γ * ∑ i : Fin M, θ i ^ 2 * w i ^ 2 / (γ - w i ^ 2) ^ 2 := by
    unfold Scalars.LwDen
    rw [hgt]
  unfold swLimit swTerm
  rw [Fin.sum_univ_one, hquad]
  unfold Scalars.Lw
  rw [if_pos h4, heta, hden]

end SingleWeight

end StackedSVD
