/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.EdgeDetR
import StackedSVD.RankR.RMT.DecompR
import StackedSVD.RankR.RMT.SymmetryR

/-!
# Task U7c, gap G1 glue: the deterministic layer

Deterministic half of the glue that shows the edge part of a spike overlap is small
(`notes/archive/rankr_plan_A.md` section 4, gap G1). The probabilistic half is
`RankR/RMT/EdgeGlueR.lean`, which supplies the eight good events and calls
`normSq_specProj_edge_le` below at one sample point.

Nothing here is random. Every statement is about one Hermitian split `S = W₀ + Q Qᵀ` of a
`d × d` matrix, one unit vector `v`, and the scalars of the plan.

## Content

1. Norm and projector helpers: `norm_add_sq_le_two`, `normSq_specProj_mono`,
   `inner_sq_le_normSq_specProjIdx` (Bessel at one sorted index).
2. Bilinear expansion of the resolvent forms on the columns of `Q`
   (`dotProduct_mulVec_expand`, `qform_mulVec_expand`, `qform2_mulVec_expand`).
3. Diagonal dominance (`le_sum_of_close_to_diag`): entries within `A / (2 r)` of
   `diag((λ_k + 1) A)` give the lower bound `(A / 2) ‖y‖²`.
4. Two resolvent bounds in the eigenbasis (`qform2_le_of_lamMax`, `le_qform2_of_psd`) and the
   reverse form bound `mul_neg_qform_le`.
5. `isUnit_det_of_lower_bound`: a lower bound on `yᵀ QᵀQ y` makes `QᵀQ` invertible.
6. The scaling bridge `specProjIdx_smul`: the projector at a sorted index does not change when
   the matrix is scaled by a positive number. Needed because the block of the rank-`r` split is
   `W₀ = d⁻¹ Bᵀ B` while `RankR/RMT/DelocR.lean` bounds `specProjIdx (Bᵀ B)`.
7. `normSq_specProj_edge_le`, the assembly: from the eight facts at one sample point to
   `‖P_edge v‖² ≤ 2 r / (bracket * μQ) + 2 ξ`.

The rank-one mirror of item 7 is the subcritical branch of `RMT/R6.lean`
(`align_tendstoInProb_of_subcritical`), which does the same split with `r = 1` and no
projector.
-/

open Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace EdgeGlueDetR

/-! ### 1. Norm and projector helpers -/

section Helpers

/-- `‖x + y‖² ≤ 2 ‖x‖² + 2 ‖y‖²`. -/
theorem norm_add_sq_le_two {E : Type*} [SeminormedAddCommGroup E] (x y : E) :
    ‖x + y‖ ^ 2 ≤ 2 * ‖x‖ ^ 2 + 2 * ‖y‖ ^ 2 := by
  have h := norm_add_le x y
  have h0 : (0 : ℝ) ≤ ‖x + y‖ := norm_nonneg _
  nlinarith [sq_nonneg (‖x‖ - ‖y‖), norm_nonneg x, norm_nonneg y]

variable {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}

/-- The squared norm of a spectral projection grows with the eigenvalue set. -/
theorem normSq_specProj_mono (hS : S.IsHermitian) {T T' : Set ℝ} (hTT : T ⊆ T')
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProj S T x‖ ^ 2 ≤ ‖specProj S T' x‖ ^ 2 := by
  classical
  rw [Frame.norm_sq_specProj_eq_sum hS, Frame.norm_sq_specProj_eq_sum hS]
  refine Finset.sum_le_sum fun i _ => ?_
  by_cases hi : hS.eigenvalues i ∈ T
  · rw [Set.indicator_of_mem hi, Set.indicator_of_mem (hTT hi)]
  · rw [Set.indicator_of_notMem hi]
    exact Set.indicator_apply_nonneg fun _ => sq_nonneg _

/-- **Bessel at one sorted index.** The squared coordinate on the eigenvector at the sorted
index `a` is at most the squared norm of the projection at that index. Mirror:
`overlapIdx_ge_inner_sq` (`LinAlg/SpecIdx.lean`), which is stated for a Gram matrix. -/
theorem inner_sq_le_normSq_specProjIdx (hS : S.IsHermitian)
    (a : Fin (Fintype.card (Fin p))) (x : EuclideanSpace ℝ (Fin p)) :
    ⟪hS.eigenvectorBasis (eigIdx p a), x⟫_ℝ ^ 2 ≤ ‖specProjIdx S hS (a : ℕ) x‖ ^ 2 := by
  classical
  rw [specProjIdx, Frame.norm_sq_specProj_eq_sum hS]
  have hmem : hS.eigenvalues (eigIdx p a) ∈ eigSetIdx S hS (a : ℕ) :=
    ⟨a, rfl, (eigenvalues_eigIdx hS a).symm⟩
  have hle := Finset.single_le_sum
    (f := fun i : Fin p => (eigSetIdx S hS (a : ℕ)).indicator
      (fun _ => ⟪hS.eigenvectorBasis i, x⟫_ℝ ^ 2) (hS.eigenvalues i))
    (fun i _ => Set.indicator_apply_nonneg fun _ => sq_nonneg _)
    (Finset.mem_univ (eigIdx p a))
  rwa [Set.indicator_of_mem hmem] at hle

/-- `u ⬝ᵥ (Q a) = (Qᵀ u) ⬝ᵥ a`. -/
theorem dotProduct_mulVec_transpose {D rr : ℕ} (Q : Matrix (Fin D) (Fin rr) ℝ)
    (u : Fin D → ℝ) (a : Fin rr → ℝ) : u ⬝ᵥ (Q *ᵥ a) = (Qᵀ *ᵥ u) ⬝ᵥ a := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]

/-- Cauchy-Schwarz in dot-product form. -/
theorem sq_dotProduct_le {rr : ℕ} (x y : Fin rr → ℝ) :
    (x ⬝ᵥ y) ^ 2 ≤ (x ⬝ᵥ x) * (y ⬝ᵥ y) := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin rr)) x y
  simpa [dotProduct, sq] using h

end Helpers

/-! ### 2. The bilinear expansion on the columns of `Q` -/

section Expand

variable {D rr : ℕ}

/-- A symmetric bilinear form on `Q y` expands over the columns of `Q`. -/
theorem dotProduct_mulVec_expand (M : Matrix (Fin D) (Fin D) ℝ)
    (Q : Matrix (Fin D) (Fin rr) ℝ) (y : Fin rr → ℝ) :
    (Q *ᵥ y) ⬝ᵥ (M *ᵥ (Q *ᵥ y))
      = ∑ k, ∑ l, y k * y l * ((fun i => Q i k) ⬝ᵥ (M *ᵥ fun i => Q i l)) := by
  have hQy : (Q *ᵥ y) = ∑ l, y l • (fun i => Q i l : Fin D → ℝ) := by
    funext i
    simp only [Matrix.mulVec, dotProduct, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    exact Finset.sum_congr rfl fun l _ => mul_comm _ _
  conv_lhs => rw [hQy]
  rw [Matrix.mulVec_sum, sum_dotProduct]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [dotProduct_sum]
  refine Finset.sum_congr rfl fun l _ => ?_
  rw [Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  ring

/-- `Φ_{Q y}` expands over the columns of `Q`. -/
theorem qform_mulVec_expand (W : Matrix (Fin D) (Fin D) ℝ) (z : ℝ)
    (Q : Matrix (Fin D) (Fin rr) ℝ) (y : Fin rr → ℝ) :
    R4.qform W z (Q *ᵥ y)
      = ∑ k, ∑ l, y k * y l * R4.cform W z (fun i => Q i k) (fun i => Q i l) :=
  dotProduct_mulVec_expand _ Q y

/-- `Φ²_{Q y}` expands over the columns of `Q`. -/
theorem qform2_mulVec_expand (W : Matrix (Fin D) (Fin D) ℝ) (z : ℝ)
    (Q : Matrix (Fin D) (Fin rr) ℝ) (y : Fin rr → ℝ) :
    R4.qform2 W z (Q *ᵥ y)
      = ∑ k, ∑ l, y k * y l * R4.cform2 W z (fun i => Q i k) (fun i => Q i l) :=
  dotProduct_mulVec_expand _ Q y

end Expand

/-! ### 3. Diagonal dominance -/

section Dominance

/-- **Diagonal dominance.** A matrix `F` whose entries sit within `A / (2 r)` of
`diag((λ_k + 1) A)`, with `A ≥ 0` and `λ_k ≥ 0`, gives the quadratic lower bound
`(A / 2) ‖y‖²`. Cauchy-Schwarz turns `(∑ |y_k|)²` into `r ‖y‖²`, which is where the factor
`2 r` in the accuracy is spent. -/
theorem le_sum_of_close_to_diag {rr : ℕ} (hrr : 0 < rr) {F : Fin rr → Fin rr → ℝ}
    {lam : Fin rr → ℝ} {A : ℝ} (hlam : ∀ k, 0 ≤ lam k) (hA : 0 ≤ A)
    (hclose : ∀ k l, |F k l - (if k = l then (lam k + 1) * A else 0)| ≤ A / (2 * rr))
    (y : Fin rr → ℝ) :
    (A / 2) * (y ⬝ᵥ y) ≤ ∑ k, ∑ l, y k * y l * F k l := by
  classical
  have hrrR : (0 : ℝ) < (rr : ℝ) := by exact_mod_cast hrr
  set E : Fin rr → Fin rr → ℝ := fun k l => F k l - (if k = l then (lam k + 1) * A else 0)
    with hEdef
  have hyy : y ⬝ᵥ y = ∑ k, y k ^ 2 := by simp [dotProduct, sq]
  have hsplit : ∑ k, ∑ l, y k * y l * F k l
      = (∑ k, y k ^ 2 * ((lam k + 1) * A)) + ∑ k, ∑ l, y k * y l * E k l := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    have hd : ∑ l, y k * y l * (if k = l then (lam k + 1) * A else 0)
        = y k ^ 2 * ((lam k + 1) * A) := by
      rw [Finset.sum_eq_single k]
      · rw [if_pos rfl]; ring
      · intro l _ hl
        rw [if_neg (fun hc => hl hc.symm), mul_zero]
      · intro h; exact absurd (Finset.mem_univ k) h
    have hexp : ∀ l : Fin rr, y k * y l * F k l
        = y k * y l * (if k = l then (lam k + 1) * A else 0) + y k * y l * E k l := by
      intro l
      rw [hEdef]
      ring
    rw [Finset.sum_congr rfl (fun l (_ : l ∈ Finset.univ) => hexp l),
      Finset.sum_add_distrib, hd]
  have hdiag : A * (∑ k, y k ^ 2) ≤ ∑ k, y k ^ 2 * ((lam k + 1) * A) := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun k _ => ?_
    have h1 : 0 ≤ y k ^ 2 * (lam k * A) :=
      mul_nonneg (sq_nonneg _) (mul_nonneg (hlam k) hA)
    nlinarith [h1]
  have hCS : (∑ k, |y k|) ^ 2 ≤ (rr : ℝ) * ∑ k, y k ^ 2 := by
    have h := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin rr))
      (fun _ : Fin rr => (1 : ℝ)) (fun k => |y k|)
    simpa [sq_abs] using h
  have hoff : |∑ k, ∑ l, y k * y l * E k l| ≤ (A / 2) * ∑ k, y k ^ 2 := by
    have hstep : |∑ k, ∑ l, y k * y l * E k l|
        ≤ ∑ k, ∑ l, |y k| * |y l| * (A / (2 * rr)) := by
      refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum fun k _ => ?_)
      refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum fun l _ => ?_)
      rw [abs_mul, abs_mul]
      refine mul_le_mul_of_nonneg_left (hclose k l) (by positivity)
    refine hstep.trans ?_
    have hrow : ∀ k : Fin rr, ∑ l, |y k| * |y l| * (A / (2 * (rr : ℝ)))
        = (A / (2 * (rr : ℝ))) * (|y k| * ∑ l, |y l|) := by
      intro k
      rw [Finset.mul_sum, Finset.mul_sum]
      exact Finset.sum_congr rfl fun l _ => by ring
    have hfac : ∑ k, ∑ l, |y k| * |y l| * (A / (2 * (rr : ℝ)))
        = (A / (2 * (rr : ℝ))) * (∑ k, |y k|) ^ 2 := by
      rw [Finset.sum_congr rfl (fun k (_ : k ∈ Finset.univ) => hrow k), ← Finset.mul_sum,
        ← Finset.sum_mul, sq]
    rw [hfac]
    have hA2 : (0 : ℝ) ≤ A / (2 * (rr : ℝ)) := by positivity
    have h1 : (A / (2 * (rr : ℝ))) * (∑ k, |y k|) ^ 2
        ≤ (A / (2 * (rr : ℝ))) * ((rr : ℝ) * ∑ k, y k ^ 2) := mul_le_mul_of_nonneg_left hCS hA2
    refine h1.trans (le_of_eq ?_)
    have hrne : (rr : ℝ) ≠ 0 := ne_of_gt hrrR
    field_simp
  have hlow : -((A / 2) * ∑ k, y k ^ 2) ≤ ∑ k, ∑ l, y k * y l * E k l :=
    neg_le_of_abs_le hoff
  rw [hsplit, hyy]
  linarith

end Dominance

/-! ### 4. Two resolvent bounds in the eigenbasis -/

section Resolvent

variable {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} {z : ℝ}

/-- Above the top eigenvalue, `Φ²_y z ≤ ‖y‖² / (z - λ_max)²`. -/
theorem qform2_le_of_lamMax (hW : W.IsHermitian) (hz : lamMax W hW < z) (y : Fin D → ℝ) :
    R4.qform2 W z y ≤ (y ⬝ᵥ y) / (z - lamMax W hW) ^ 2 := by
  have hgap : 0 < z - lamMax W hW := by linarith
  rw [R4.qform2_eq_sum hW hz, ← R4.dotProduct_transpose_eigU hW y]
  have hyy : ((R4.eigU hW)ᵀ *ᵥ y) ⬝ᵥ ((R4.eigU hW)ᵀ *ᵥ y)
      = ∑ a, ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 := by simp [dotProduct, sq]
  rw [hyy, Finset.sum_div]
  refine Finset.sum_le_sum fun a _ => ?_
  have hla : hW.eigenvalues a ≤ lamMax W hW := R4.eigenvalues_le_lamMax hW a
  have hpos : 0 < z - hW.eigenvalues a := by linarith
  have hkey : (hW.eigenvalues a - z)⁻¹ ^ 2 ≤ ((z - lamMax W hW) ^ 2)⁻¹ := by
    have h1 : (z - lamMax W hW) ^ 2 ≤ (hW.eigenvalues a - z) ^ 2 := by nlinarith
    have h2 : (0 : ℝ) < (z - lamMax W hW) ^ 2 := pow_pos hgap 2
    rw [inv_pow]
    exact inv_anti₀ h2 h1
  calc (hW.eigenvalues a - z)⁻¹ ^ 2 * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2
      ≤ ((z - lamMax W hW) ^ 2)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 :=
        mul_le_mul_of_nonneg_right hkey (sq_nonneg _)
    _ = ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 / (z - lamMax W hW) ^ 2 := by
        rw [div_eq_mul_inv, mul_comm]

/-- For a positive semidefinite `W` and `z` above the top eigenvalue,
`‖y‖² / z² ≤ Φ²_y z`. -/
theorem le_qform2_of_psd (hW : W.IsHermitian) (hpsd : ∀ a, 0 ≤ hW.eigenvalues a)
    (hz0 : 0 < z) (hz : lamMax W hW < z) (y : Fin D → ℝ) :
    (y ⬝ᵥ y) / z ^ 2 ≤ R4.qform2 W z y := by
  rw [R4.qform2_eq_sum hW hz, ← R4.dotProduct_transpose_eigU hW y]
  have hyy : ((R4.eigU hW)ᵀ *ᵥ y) ⬝ᵥ ((R4.eigU hW)ᵀ *ᵥ y)
      = ∑ a, ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 := by simp [dotProduct, sq]
  rw [hyy, Finset.sum_div]
  refine Finset.sum_le_sum fun a _ => ?_
  have hla : hW.eigenvalues a ≤ lamMax W hW := R4.eigenvalues_le_lamMax hW a
  have hpos : 0 < z - hW.eigenvalues a := by linarith
  have hkey : (z ^ 2)⁻¹ ≤ (hW.eigenvalues a - z)⁻¹ ^ 2 := by
    have h1 : (hW.eigenvalues a - z) ^ 2 ≤ z ^ 2 := by nlinarith [hpsd a]
    have h2 : (0 : ℝ) < (hW.eigenvalues a - z) ^ 2 := by nlinarith [mul_pos hpos hpos]
    rw [inv_pow]
    exact inv_anti₀ h2 h1
  calc ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 / z ^ 2
      = (z ^ 2)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 := by rw [div_eq_mul_inv, mul_comm]
    _ ≤ (hW.eigenvalues a - z)⁻¹ ^ 2 * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 :=
        mul_le_mul_of_nonneg_right hkey (sq_nonneg _)

/-- The reverse of `OutliersR.dot_self_le_mul_neg_qform`: the negative form is small when the
gap is wide. No positivity of `W` is needed. -/
theorem mul_neg_qform_le (hW : W.IsHermitian) (hz : lamMax W hW < z) (y : Fin D → ℝ) :
    (z - lamMax W hW) * (-R4.qform W z y) ≤ y ⬝ᵥ y := by
  have hgap : 0 < z - lamMax W hW := by linarith
  rw [R4.qform_eq_sum hW hz, ← R4.dotProduct_transpose_eigU hW y]
  have hyy : ((R4.eigU hW)ᵀ *ᵥ y) ⬝ᵥ ((R4.eigU hW)ᵀ *ᵥ y)
      = ∑ a, ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 := by simp [dotProduct, sq]
  rw [hyy, ← Finset.sum_neg_distrib, Finset.mul_sum]
  refine Finset.sum_le_sum fun a _ => ?_
  have hla : hW.eigenvalues a ≤ lamMax W hW := R4.eigenvalues_le_lamMax hW a
  have hpos : 0 < z - hW.eigenvalues a := by linarith
  have hrw : -((hW.eigenvalues a - z)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2)
      = (z - hW.eigenvalues a)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2 := by
    rw [show (hW.eigenvalues a - z)⁻¹ = -(z - hW.eigenvalues a)⁻¹ by
      rw [← neg_sub z (hW.eigenvalues a), inv_neg]]
    ring
  rw [hrw]
  have hinv : (z - hW.eigenvalues a)⁻¹ ≤ (z - lamMax W hW)⁻¹ :=
    inv_anti₀ hgap (by linarith)
  have hmul : (z - lamMax W hW) * ((z - hW.eigenvalues a)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2)
      ≤ (z - lamMax W hW) * ((z - lamMax W hW)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ y) a ^ 2) := by
    refine mul_le_mul_of_nonneg_left ?_ hgap.le
    exact mul_le_mul_of_nonneg_right hinv (sq_nonneg _)
  refine hmul.trans (le_of_eq ?_)
  rw [← mul_assoc, mul_inv_cancel₀ hgap.ne', one_mul]

end Resolvent

/-! ### 5. Invertibility of `QᵀQ` -/

section Invertible

/-- A strict lower bound on `yᵀ QᵀQ y` makes `QᵀQ` invertible. -/
theorem isUnit_det_of_lower_bound {D rr : ℕ} (Q : Matrix (Fin D) (Fin rr) ℝ) {μQ : ℝ}
    (hμ : 0 < μQ) (h : ∀ y : Fin rr → ℝ, μQ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y)) :
    IsUnit (Qᵀ * Q).det := by
  rw [isUnit_iff_ne_zero]
  intro hdet
  obtain ⟨y, hy0, hy⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  have h1 := h y
  rw [hy, dotProduct_zero] at h1
  have h2 : y ⬝ᵥ y ≤ 0 := by nlinarith [h1, hμ]
  have h3 : y = 0 := by
    have hnn : 0 ≤ y ⬝ᵥ y := DecompR.dotProduct_self_nonneg' y
    have heq : ∑ i, y i * y i = 0 := le_antisymm h2 hnn
    funext i
    have hi := (Finset.sum_eq_zero_iff_of_nonneg
      (fun i (_ : i ∈ Finset.univ) => mul_self_nonneg (y i))).mp heq i (Finset.mem_univ i)
    exact mul_self_eq_zero.mp hi
  exact hy0 h3

end Invertible

/-! ### 6. The scaling bridge -/

section Scaling

variable {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} {t : ℝ}

open Polynomial

/-- Scaling a Hermitian matrix by `t > 0` scales every sorted eigenvalue by `t`. Copy of
`Eigen.eigenvalues₀_add_one` with the shift replaced by the scaling. -/
theorem eigenvalues₀_smul (hA : A.IsHermitian) (hAt : (t • A).IsHermitian) (ht : 0 < t)
    (k : Fin (Fintype.card (Fin p))) : hAt.eigenvalues₀ k = t * hA.eigenvalues₀ k := by
  have hUU : star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 := Unitary.coe_star_mul_self _
  have hd : Matrix.diagonal (fun i => t * hA.eigenvalues i)
      = t • Matrix.diagonal hA.eigenvalues := by
    ext i j
    by_cases hij : i = j
    · subst hij; simp [Matrix.diagonal]
    · simp [Matrix.diagonal, hij]
  have hscale : t • A = (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      Matrix.diagonal (fun i => t * hA.eigenvalues i) *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) := by
    rw [hd, Matrix.mul_smul, Matrix.smul_mul, ← spectral_conj hA]
  have hchar : (t • A).charpoly = ∏ i, (X - C (t * hA.eigenvalues i)) :=
    charpoly_eq_prod_of_conj hUU hscale
  have hchar' : (t • A).charpoly
      = ∏ k : Fin (Fintype.card (Fin p)), (X - C (t * hA.eigenvalues₀ k)) := by
    rw [hchar]
    refine (Fintype.prod_equiv (eigIdx p) _ _ ?_).symm
    intro j
    rw [eigenvalues_eigIdx]
  have hanti : Antitone fun k : Fin (Fintype.card (Fin p)) => t * hA.eigenvalues₀ k :=
    fun a b hab => by
      exact mul_le_mul_of_nonneg_left (hA.eigenvalues₀_antitone hab) ht.le
  exact congrFun (eigenvalues₀_eq_of_charpoly hAt hanti hchar') k

/-- Equal subspaces give equal orthogonal projectors. The instance argument is a `Prop`, so
proof irrelevance closes the goal after the substitution. -/
theorem starProjection_congr {K K' : Submodule ℝ (EuclideanSpace ℝ (Fin p))}
    [K.HasOrthogonalProjection] [K'.HasOrthogonalProjection] (h : K = K') :
    K.starProjection = K'.starProjection := by
  subst h
  rfl

/-- The eigenspace of a scaled operator at the scaled eigenvalue. -/
theorem specSpace_smul (ht : t ≠ 0) (A : Matrix (Fin p) (Fin p) ℝ) (lam : ℝ) :
    specSpace (t • A) {t * lam} = specSpace A {lam} := by
  have hspace : ∀ (B : Matrix (Fin p) (Fin p) ℝ) (x : ℝ),
      specSpace B {x} = Module.End.eigenspace (toOp B) x := by
    intro B x
    simp only [specSpace]
    exact iSup_singleton
  rw [hspace (t • A) (t * lam), hspace A lam]
  ext x
  rw [Module.End.mem_eigenspace_iff, Module.End.mem_eigenspace_iff]
  have hop : (toOp (t • A)) x = t • (toOp A) x := by
    have hm : toOp (t • A) = t • toOp A := map_smul _ t A
    rw [hm, LinearMap.smul_apply]
  rw [hop, mul_smul]
  exact ⟨fun h => smul_right_injective (EuclideanSpace ℝ (Fin p)) ht h, fun h => by rw [h]⟩

/-- **The scaling bridge.** The projector at a sorted eigenvalue index is unchanged when the
matrix is scaled by a positive number. The rank-`r` split has `W₀ = d⁻¹ Bᵀ B` while
`lintegral_prod_normSq_specProjIdx_le` bounds the projector of `Bᵀ B`. -/
theorem specProjIdx_smul (hA : A.IsHermitian) (hAt : (t • A).IsHermitian) (ht : 0 < t)
    (k : ℕ) : specProjIdx (t • A) hAt k = specProjIdx A hA k := by
  by_cases hk : k < Fintype.card (Fin p)
  · rw [specProjIdx, specProjIdx, eigSetIdx_eq_singleton _ _ hk, eigSetIdx_eq_singleton _ _ hk,
      eigenvalues₀_smul hA hAt ht ⟨k, hk⟩]
    change (specSpace (t • A) {t * hA.eigenvalues₀ ⟨k, hk⟩}).starProjection
      = (specSpace A {hA.eigenvalues₀ ⟨k, hk⟩}).starProjection
    exact starProjection_congr (specSpace_smul ht.ne' A _)
  · rw [specProjIdx, specProjIdx, eigSetIdx_eq_empty _ _ hk, eigSetIdx_eq_empty _ _ hk]
    change (specSpace (t • A) (∅ : Set ℝ)).starProjection
      = (specSpace A (∅ : Set ℝ)).starProjection
    exact starProjection_congr (by rw [specSpace_empty, specSpace_empty])

end Scaling

/-! ### 7. The κ bridge and the assembly -/

section Assembly

variable {D rr : ℕ}

/-- The `hκ` hypothesis of `EdgeDetR.normSq_transpose_mulVec_le_of_edge` follows from a bound
on the double sum of projector norms over the sorted indices below `rr` and the columns of
`Q`. Bessel (`inner_sq_le_normSq_specProjIdx`) does the work. -/
theorem normSq_transpose_mulVec_eigvec_le {W₀ : Matrix (Fin D) (Fin D) ℝ}
    (hW₀ : W₀.IsHermitian) (Q : Matrix (Fin D) (Fin rr) ℝ) {κ : ℝ}
    (hκ : ∑ a ∈ Finset.range rr, ∑ l, ‖specProjIdx W₀ hW₀ a
        (WithLp.toLp 2 fun i => Q i l)‖ ^ 2 ≤ κ)
    (a : Fin (Fintype.card (Fin D))) (ha : (a : ℕ) < rr) :
    ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hW₀.eigenvectorBasis (eigIdx D a)))
        : EuclideanSpace ℝ (Fin rr))‖ ^ 2 ≤ κ := by
  classical
  set e : EuclideanSpace ℝ (Fin D) := hW₀.eigenvectorBasis (eigIdx D a) with hedef
  have hexp : ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp e) : EuclideanSpace ℝ (Fin rr))‖ ^ 2
      = ∑ l, ⟪e, (WithLp.toLp 2 fun i => Q i l : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2 := by
    rw [EdgeDetR.norm_toLp_sq]
    simp only [dotProduct, sq]
    refine Finset.sum_congr rfl fun l _ => ?_
    have hcoord : (Qᵀ *ᵥ WithLp.ofLp e) l
        = ⟪e, (WithLp.toLp 2 fun i => Q i l : EuclideanSpace ℝ (Fin D))⟫_ℝ := by
      rw [Frame.inner_eq_dot]
      simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
      exact Finset.sum_congr rfl fun i _ => mul_comm _ _
    rw [hcoord]
  rw [hexp]
  have hstep : ∑ l, ⟪e, (WithLp.toLp 2 fun i => Q i l : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2
      ≤ ∑ l, ‖specProjIdx W₀ hW₀ (a : ℕ) (WithLp.toLp 2 fun i => Q i l)‖ ^ 2 :=
    Finset.sum_le_sum fun l _ => inner_sq_le_normSq_specProjIdx hW₀ a _
  refine hstep.trans (le_trans ?_ hκ)
  refine Finset.single_le_sum
    (f := fun b : ℕ => ∑ l, ‖specProjIdx W₀ hW₀ b (WithLp.toLp 2 fun i => Q i l)‖ ^ 2)
    (fun b _ => Finset.sum_nonneg fun l _ => sq_nonneg _) (Finset.mem_range.mpr ha)

/-- **The deterministic assembly of gap G1.** At one sample point: the squared norm of the
projection of the unit vector `v` on the top-`rr` eigenvalues of `S` that sit at or below `τ`
is at most `2 rr / (bracket * μQ) + 2 ξ`, where `bracket = μmin - rr κ / ε₀²`.

The eight inputs are the good events of `RankR/RMT/EdgeGlueR.lean`: the edge bound `hedge`,
the two simplicity facts, the κ bound, the two form bounds `hmin` and `hQQ`, the residual
bound `hrest`, and the size condition `hrD`. -/
theorem normSq_specProj_edge_le {W₀ S : Matrix (Fin D) (Fin D) ℝ} {Q : Matrix (Fin D) (Fin rr) ℝ}
    (hW₀ : W₀.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W₀ + Q * Qᵀ)
    (hsimpleW : SimpleSpec W₀ hW₀ rr) (hsimpleS : SimpleSpec S hS rr) (hrD : rr ≤ D)
    {z₀ ε₀ τ κ μmin μQ ξ : ℝ}
    (hε₀ : 0 < ε₀) (hedge : lamMax W₀ hW₀ + ε₀ ≤ z₀) (hτ : τ ≤ z₀)
    (hκ : ∑ a ∈ Finset.range rr, ∑ l, ‖specProjIdx W₀ hW₀ a
        (WithLp.toLp 2 fun i => Q i l)‖ ^ 2 ≤ κ)
    (hmin : ∀ y : Fin rr → ℝ, μmin * (y ⬝ᵥ y) ≤ R4.qform2 W₀ z₀ (Q *ᵥ y))
    (hbr : 0 < μmin - (rr : ℝ) * κ / ε₀ ^ 2)
    (hμQ : 0 < μQ) (hQQ : ∀ y : Fin rr → ℝ, μQ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y))
    {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1)
    (hrest : ∑ k ∈ Finset.range rr,
      ‖specProjIdx S hS k (WithLp.toLp 2 (DecompR.rvecOfR Q (WithLp.ofLp v)))‖ ^ 2 ≤ ξ) :
    ‖specProj S (topEigSet S hS rr ∩ Set.Iic τ) v‖ ^ 2
      ≤ 2 * ((rr : ℝ) / ((μmin - (rr : ℝ) * κ / ε₀ ^ 2) * μQ)) + 2 * ξ := by
  classical
  set br : ℝ := μmin - (rr : ℝ) * κ / ε₀ ^ 2 with hbrdef
  set T : Set ℝ := topEigSet S hS rr ∩ Set.Iic τ with hTdef
  have hunit : IsUnit (Qᵀ * Q).det := isUnit_det_of_lower_bound Q hμQ hQQ
  have hκ0 : 0 ≤ κ :=
    le_trans (Finset.sum_nonneg fun a _ => Finset.sum_nonneg fun l _ => sq_nonneg _) hκ
  set av : Fin rr → ℝ := DecompR.aOfR Q (WithLp.ofLp v) with havdef
  set rv : Fin D → ℝ := DecompR.rvecOfR Q (WithLp.ofLp v) with hrvdef
  -- 1. the two pieces of `v`
  have hvv : WithLp.ofLp v ⬝ᵥ WithLp.ofLp v = 1 := by
    have h := EdgeDetR.norm_toLp_sq (WithLp.ofLp v)
    rw [WithLp.toLp_ofLp, hv] at h
    simpa using h.symm
  have hsplit : v = (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))
      + WithLp.toLp 2 rv := by
    have h := DecompR.toLp_eq_mulVec_add Q (WithLp.ofLp v)
    rwa [WithLp.toLp_ofLp] at h
  -- 2. the residual part
  have hres : ‖specProj S T (WithLp.toLp 2 rv)‖ ^ 2 ≤ ξ := by
    refine le_trans (normSq_specProj_mono hS Set.inter_subset_left _) ?_
    have h := normSq_specProjTop_eq_sum_of_simpleSpec hsimpleS hrD
      (WithLp.toLp 2 rv : EuclideanSpace ℝ (Fin D))
    rw [specProjTop] at h
    rw [h]
    exact hrest
  -- 3. the frame part
  have habnd : av ⬝ᵥ av ≤ 1 / μQ := by
    have h := DecompR.aOfR_sq_mul_le Q hunit hQQ (v := WithLp.ofLp v)
    rw [hvv] at h
    rw [le_div_iff₀ hμQ]
    nlinarith [h]
  set Fs : Finset (Fin D) := Finset.univ.filter fun i => hS.eigenvalues i ∈ T with hFsdef
  have hcard : Fs.card ≤ rr := by
    have hmap : ∀ i ∈ Fs, ((eigIdx D).symm i : ℕ) ∈ Finset.range rr := by
      intro i hi
      rw [hFsdef, Finset.mem_filter] at hi
      have hmem : hS.eigenvalues i ∈ topEigSet S hS rr := hi.2.1
      have hev : hS.eigenvalues₀ ((eigIdx D).symm i) = hS.eigenvalues i := by
        rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
      rw [Finset.mem_range]
      exact (mem_topEigSet_of_simpleSpec hsimpleS ((eigIdx D).symm i)).mp (hev ▸ hmem)
    have hinj : Set.InjOn (fun i => ((eigIdx D).symm i : ℕ)) Fs := by
      intro i _ i' _ h
      exact (eigIdx D).symm.injective (Fin.val_injective h)
    have := Finset.card_le_card_of_injOn _ hmap hinj
    simpa using this
  have hterm : ∀ i ∈ Fs,
      ⟪hS.eigenvectorBasis i, (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2
        ≤ 1 / (br * μQ) := by
    intro i hi
    rw [hFsdef, Finset.mem_filter] at hi
    have hmem : hS.eigenvalues i ∈ topEigSet S hS rr := hi.2.1
    have hIic : hS.eigenvalues i ≤ τ := hi.2.2
    have hev : hS.eigenvalues₀ ((eigIdx D).symm i) = hS.eigenvalues i := by
      rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
    have hidx : (((eigIdx D).symm i : Fin (Fintype.card (Fin D))) : ℕ) < rr :=
      (mem_topEigSet_of_simpleSpec hsimpleS ((eigIdx D).symm i)).mp (hev ▸ hmem)
    -- the P3 bound on `‖Qᵀ u‖²`
    have hlamMax : lamMax W₀ hW₀ < z₀ := by linarith
    have hP3 := EdgeDetR.normSq_transpose_mulVec_le_of_edge hW₀ hS hSeq hlamMax hsimpleW
      (u := hS.eigenvectorBasis i) (hS.eigenvectorBasis.orthonormal.1 i)
      (lam := hS.eigenvalues i) (Frame.toOp_eigvec hS i)
      ⟨(eigIdx D).symm i, hidx, hev⟩ (le_trans hIic hτ)
      (κ := κ) (normSq_transpose_mulVec_eigvec_le hW₀ Q hκ) hmin
    have hgz : 0 < z₀ - lamMax W₀ hW₀ := by linarith
    have hgap : ε₀ ^ 2 ≤ (z₀ - lamMax W₀ hW₀) ^ 2 := by nlinarith
    have hbr' : br ≤ μmin - (rr : ℝ) * κ / (z₀ - lamMax W₀ hW₀) ^ 2 := by
      rw [hbrdef]
      have hnum : (0 : ℝ) ≤ (rr : ℝ) * κ := by positivity
      have hden : (0 : ℝ) < ε₀ ^ 2 := pow_pos hε₀ 2
      have hinv : ((z₀ - lamMax W₀ hW₀) ^ 2)⁻¹ ≤ (ε₀ ^ 2)⁻¹ := inv_anti₀ hden hgap
      have h1 : (rr : ℝ) * κ / (z₀ - lamMax W₀ hW₀) ^ 2 ≤ (rr : ℝ) * κ / ε₀ ^ 2 := by
        rw [div_eq_mul_inv, div_eq_mul_inv]
        exact mul_le_mul_of_nonneg_left hinv hnum
      linarith
    have hnnQ : (0 : ℝ) ≤ ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
        : EuclideanSpace ℝ (Fin rr))‖ ^ 2 := sq_nonneg _
    have hQu : ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
        : EuclideanSpace ℝ (Fin rr))‖ ^ 2 ≤ 1 / br := by
      rw [le_div_iff₀ hbr]
      nlinarith [hP3, hbr', hnnQ, mul_nonneg hnnQ (sub_nonneg.mpr hbr')]
    -- Cauchy-Schwarz
    have hinner : ⟪hS.eigenvectorBasis i,
        (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ
        = (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i)) ⬝ᵥ av := by
      rw [Frame.inner_eq_dot, WithLp.ofLp_toLp, dotProduct_mulVec_transpose]
    rw [hinner]
    have hCS := sq_dotProduct_le (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i)) av
    have hnorm : (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
        ⬝ᵥ (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
        = ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
            : EuclideanSpace ℝ (Fin rr))‖ ^ 2 := (EdgeDetR.norm_toLp_sq _).symm
    rw [hnorm] at hCS
    have hnn : 0 ≤ av ⬝ᵥ av := DecompR.dotProduct_self_nonneg' av
    have hmul : ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hS.eigenvectorBasis i))
          : EuclideanSpace ℝ (Fin rr))‖ ^ 2 * (av ⬝ᵥ av) ≤ (1 / br) * (1 / μQ) :=
      mul_le_mul hQu habnd hnn (by positivity)
    have hfinal : (1 / br) * (1 / μQ) = 1 / (br * μQ) := by
      rw [div_mul_div_comm, one_mul]
    linarith [hCS, hmul, hfinal]
  have hframe : ‖specProj S T (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))‖ ^ 2
      ≤ (rr : ℝ) / (br * μQ) := by
    rw [Frame.norm_sq_specProj_eq_sum hS]
    have hg0 : ∀ i ∈ (Finset.univ : Finset (Fin D)), i ∉ Fs →
        T.indicator (fun _ =>
          ⟪hS.eigenvectorBasis i, (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2)
          (hS.eigenvalues i) = 0 := by
      intro i _ hi
      have hnot : hS.eigenvalues i ∉ T := by
        intro hc
        exact hi (by rw [hFsdef]; exact Finset.mem_filter.mpr ⟨Finset.mem_univ i, hc⟩)
      exact Set.indicator_of_notMem hnot _
    have hsum : ∑ i, T.indicator (fun _ =>
        ⟪hS.eigenvectorBasis i, (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2)
          (hS.eigenvalues i)
        = ∑ i ∈ Fs,
          ⟪hS.eigenvectorBasis i,
            (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2 := by
      rw [← Finset.sum_subset (Finset.subset_univ Fs) hg0]
      refine Finset.sum_congr rfl fun i hi => ?_
      have hmem : hS.eigenvalues i ∈ T := by
        rw [hFsdef] at hi
        exact (Finset.mem_filter.mp hi).2
      rw [Set.indicator_of_mem hmem]
    rw [hsum]
    have hbnd : ∑ i ∈ Fs,
        ⟪hS.eigenvectorBasis i, (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2
        ≤ (Fs.card : ℝ) * (1 / (br * μQ)) := by
      calc ∑ i ∈ Fs,
            ⟪hS.eigenvectorBasis i, (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))⟫_ℝ ^ 2
          ≤ ∑ _i ∈ Fs, (1 / (br * μQ)) := Finset.sum_le_sum hterm
        _ = (Fs.card : ℝ) * (1 / (br * μQ)) := by
            rw [Finset.sum_const, nsmul_eq_mul]
    refine hbnd.trans ?_
    have hcardR : (Fs.card : ℝ) ≤ (rr : ℝ) := by exact_mod_cast hcard
    have hpos : (0 : ℝ) < br * μQ := by positivity
    have heq : (rr : ℝ) / (br * μQ) = (rr : ℝ) * (1 / (br * μQ)) := by rw [mul_one_div]
    rw [heq]
    exact mul_le_mul_of_nonneg_right hcardR (by positivity)
  -- 4. assemble
  calc ‖specProj S T v‖ ^ 2
      = ‖specProj S T (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))
          + specProj S T (WithLp.toLp 2 rv)‖ ^ 2 := by
        rw [← map_add, ← hsplit]
    _ ≤ 2 * ‖specProj S T (WithLp.toLp 2 (Q *ᵥ av) : EuclideanSpace ℝ (Fin D))‖ ^ 2
          + 2 * ‖specProj S T (WithLp.toLp 2 rv)‖ ^ 2 := norm_add_sq_le_two _ _
    _ ≤ 2 * ((rr : ℝ) / (br * μQ)) + 2 * ξ := by
        have h1 := hframe
        have h2 := hres
        linarith

end Assembly

end EdgeGlueDetR

end StackedSVD
