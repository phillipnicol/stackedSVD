/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.SpecialFunctions.LogSumExp
public import Mathlib.Analysis.Calculus.Gradient.Basic
public import Mathlib.Analysis.InnerProductSpace.Calculus

set_option autoImplicit false

/-!
# Log-sum-exp smoothing of a maximum of linear forms

For a finite family `x : ι → E` of vectors of a real inner product space and `β > 0`, the
log-sum-exp function `logSumExp β x v = (1/β) log ∑ i, exp (β ⟪x i, v⟫)` is a `C¹` smoothing of
the maximum `v ↦ ⨆ i, ⟪x i, v⟫`. It is the composition of the log-sum-exp function
`Real.logSumExp β` of `COLT83.Mathlib.Analysis.SpecialFunctions.LogSumExp` with the linear map
`innerPi x : v ↦ (⟪x i, v⟫)ᵢ`, and the softmax weights `softmax β x v` are those of the vector
`innerPi x v`:
* it is sandwiched between the maximum and the maximum plus `log (card ι) / β`
  (`ciSup_inner_le_logSumExp`, `logSumExp_le_ciSup_inner_add`);
* its gradient is the convex combination `∑ i, softmax β x v i • x i` of the vectors `x i` with
  the softmax weights (`hasGradientAt_logSumExp`), hence has norm at most `sup_i ‖x i‖`
  (`norm_gradient_logSumExp_le`);
* it is smooth (`contDiff_logSumExp'`), with second derivative the softmax-weighted covariance of
  the linear forms `⟪x i, ·⟫` (`fderiv_fderiv_logSumExp_apply`), whose operator norm is at most
  `β * (sup_i ‖x i‖) ^ 2` (`norm_fderiv_fderiv_logSumExp_le`).
-/

@[expose] public section
set_option autoImplicit false

open InnerProductSpace Finset
open scoped RealInnerProductSpace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] {ι : Type*}

/-- The continuous linear map `v ↦ (⟪x i, v⟫)ᵢ` of a family of vectors `x`. -/
noncomputable def innerPi (x : ι → E) : E →L[ℝ] (ι → ℝ) :=
  ContinuousLinearMap.pi fun i ↦ innerSL ℝ (x i)

@[simp] lemma innerPi_apply (x : ι → E) (v : E) (i : ι) : innerPi x v i = ⟪x i, v⟫ := rfl

variable [Fintype ι] {β : ℝ}

/-- The softmax weights `exp (β ⟪x i, v⟫) / ∑ j, exp (β ⟪x j, v⟫)`. -/
noncomputable def softmax (β : ℝ) (x : ι → E) (v : E) (i : ι) : ℝ :=
  Real.softmax β (innerPi x v) i

/-- The log-sum-exp smoothing `(1/β) log ∑ i, exp (β ⟪x i, v⟫)` of `v ↦ max_i ⟪x i, v⟫`. -/
noncomputable def logSumExp (β : ℝ) (x : ι → E) (v : E) : ℝ :=
  Real.logSumExp β (innerPi x v)

lemma softmax_apply (β : ℝ) (x : ι → E) (v : E) (i : ι) :
    softmax β x v i = Real.exp (β * ⟪x i, v⟫) / ∑ j, Real.exp (β * ⟪x j, v⟫) := rfl

lemma logSumExp_apply (β : ℝ) (x : ι → E) (v : E) :
    logSumExp β x v = β⁻¹ * Real.log (∑ i, Real.exp (β * ⟪x i, v⟫)) := rfl

/-- `logSumExp β x` is the composition of `Real.logSumExp β` with `innerPi x`. -/
lemma logSumExp_eq_comp (β : ℝ) (x : ι → E) : logSumExp β x = Real.logSumExp β ∘ innerPi x := rfl

/-- The softmax weights are nonnegative. -/
lemma softmax_nonneg (β : ℝ) (x : ι → E) (v : E) (i : ι) : 0 ≤ softmax β x v i :=
  Real.softmax_nonneg β _ i

variable [Nonempty ι]

/-- The denominator of the softmax weights is positive. -/
lemma sum_exp_inner_pos (β : ℝ) (x : ι → E) (v : E) : 0 < ∑ i, Real.exp (β * ⟪x i, v⟫) :=
  Real.sum_exp_mul_pos β _

/-- The softmax weights sum to one. -/
lemma sum_softmax (β : ℝ) (x : ι → E) (v : E) : ∑ i, softmax β x v i = 1 :=
  Real.sum_softmax β _

/-- The log-sum-exp function is smooth. -/
lemma contDiff_logSumExp' (n : WithTop ℕ∞) (β : ℝ) (x : ι → E) :
    ContDiff ℝ n (logSumExp β x) :=
  (Real.contDiff_logSumExp n β).comp (innerPi x).contDiff

/-- The log-sum-exp function is `C¹`. -/
lemma contDiff_logSumExp (β : ℝ) (x : ι → E) : ContDiff ℝ 1 (logSumExp β x) :=
  contDiff_logSumExp' 1 β x

/-- The log-sum-exp function dominates the maximum `⨆ i, ⟪x i, v⟫`. -/
lemma ciSup_inner_le_logSumExp (hβ : 0 < β) (x : ι → E) (v : E) :
    (⨆ i, ⟪x i, v⟫) ≤ logSumExp β x v :=
  Real.ciSup_le_logSumExp hβ _

/-- The log-sum-exp function exceeds the maximum `⨆ i, ⟪x i, v⟫` by at most
`log (card ι) / β`. -/
lemma logSumExp_le_ciSup_inner_add (hβ : 0 < β) (x : ι → E) (v : E) :
    logSumExp β x v ≤ (⨆ i, ⟪x i, v⟫) + Real.log (Fintype.card ι) / β :=
  Real.logSumExp_le_ciSup_add hβ _

section fderiv

/-- The Fréchet derivative of the log-sum-exp function is the softmax-weighted average
`∑ i, softmax β x v i • ⟪x i, ·⟫` of the linear forms `⟪x i, ·⟫`. -/
lemma hasFDerivAt_logSumExp (hβ : β ≠ 0) (x : ι → E) (v : E) :
    HasFDerivAt (logSumExp β x) (∑ i, softmax β x v i • innerSL ℝ (x i)) v := by
  have h := (Real.hasFDerivAt_logSumExp hβ (innerPi x v)).comp v (innerPi x).hasFDerivAt
  refine h.congr_fderiv ?_
  ext w
  simp [softmax]

/-- The Fréchet derivative of the log-sum-exp function. -/
lemma fderiv_logSumExp (hβ : β ≠ 0) (x : ι → E) (v : E) :
    fderiv ℝ (logSumExp β x) v = ∑ i, softmax β x v i • innerSL ℝ (x i) :=
  (hasFDerivAt_logSumExp hβ x v).fderiv

/-- The Fréchet derivative of the softmax weights:
`∂ p_i = β p_i (⟪x i, ·⟫ - ∑ j, p_j ⟪x j, ·⟫)`. -/
lemma hasFDerivAt_softmax (β : ℝ) (x : ι → E) (v : E) (i : ι) :
    HasFDerivAt (fun v ↦ softmax β x v i)
      ((β * softmax β x v i) • (innerSL ℝ (x i) - ∑ j, softmax β x v j • innerSL ℝ (x j))) v := by
  have h := (Real.hasFDerivAt_softmax β (innerPi x v) i).comp v (innerPi x).hasFDerivAt
  refine h.congr_fderiv ?_
  ext w
  simp [softmax]

/-- The second derivative of the log-sum-exp function, as a bilinear form:
`D²F(v)(u, w) = β (∑ i, p_i ⟪x i, u⟫ ⟪x i, w⟫ - (∑ i, p_i ⟪x i, u⟫) (∑ i, p_i ⟪x i, w⟫))`. -/
lemma hasFDerivAt_fderiv_logSumExp (hβ : β ≠ 0) (x : ι → E) (v : E) :
    HasFDerivAt (fderiv ℝ (logSumExp β x))
      (β • (∑ i, softmax β x v i • (innerSL ℝ (x i)).smulRight (innerSL ℝ (x i)) -
        (∑ i, softmax β x v i • innerSL ℝ (x i)).smulRight
          (∑ i, softmax β x v i • innerSL ℝ (x i)))) v := by
  have h : fderiv ℝ (logSumExp β x) = fun v ↦ ∑ i, softmax β x v i • innerSL ℝ (x i) :=
    funext (fderiv_logSumExp hβ x)
  rw [h]
  refine (HasFDerivAt.fun_sum fun i _ ↦
    (hasFDerivAt_softmax β x v i).smul_const (innerSL ℝ (x i))).congr_fderiv ?_
  ext u w
  simp only [smul_apply, sub_apply, FunLike.coe_sum, Finset.sum_apply, innerSL_apply_apply,
    ContinuousLinearMap.smulRight_apply, smul_eq_mul]
  rw [mul_sub, mul_sum, mul_sum, mul_sum, ← sum_sub_distrib]
  exact sum_congr rfl fun i _ ↦ by ring

/-- The second derivative of the log-sum-exp function, applied to two vectors. -/
lemma fderiv_fderiv_logSumExp_apply (hβ : β ≠ 0) (x : ι → E) (v u w : E) :
    fderiv ℝ (fderiv ℝ (logSumExp β x)) v u w =
      β * (∑ i, softmax β x v i * (⟪x i, u⟫ * ⟪x i, w⟫) -
        (∑ i, softmax β x v i * ⟪x i, u⟫) * (∑ i, softmax β x v i * ⟪x i, w⟫)) := by
  rw [(hasFDerivAt_fderiv_logSumExp hβ x v).fderiv]
  simp only [smul_apply, sub_apply, FunLike.coe_sum, Finset.sum_apply, innerSL_apply_apply,
    ContinuousLinearMap.smulRight_apply, smul_eq_mul]

/-- If `‖x i‖ ≤ R` for all `i`, the second derivative of the log-sum-exp function has operator
norm at most `β * R ^ 2`. -/
lemma norm_fderiv_fderiv_logSumExp_le (hβ : 0 < β) {x : ι → E} {R : ℝ} (hR : ∀ i, ‖x i‖ ≤ R)
    (v : E) : ‖fderiv ℝ (fderiv ℝ (logSumExp β x)) v‖ ≤ β * R ^ 2 := by
  have hR0 : 0 ≤ R := (norm_nonneg _).trans (hR (Classical.arbitrary ι))
  refine ContinuousLinearMap.opNorm_le_bound₂ _ (by positivity) fun u w ↦ ?_
  rw [fderiv_fderiv_logSumExp_apply hβ.ne' x v u w, Real.norm_eq_abs, abs_mul, abs_of_pos hβ]
  have hbound (y : E) (i : ι) : |⟪x i, y⟫| ≤ R * ‖y‖ :=
    (abs_real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hR i) (norm_nonneg _))
  have h := abs_sum_mul_sub_mul_le (s := univ) (fun i _ ↦ softmax_nonneg β x v i)
    (sum_softmax β x v) (fun i _ ↦ hbound u i) (fun i _ ↦ hbound w i)
  calc β * |∑ i, softmax β x v i * (⟪x i, u⟫ * ⟪x i, w⟫) -
        (∑ i, softmax β x v i * ⟪x i, u⟫) * (∑ i, softmax β x v i * ⟪x i, w⟫)|
      ≤ β * (R * ‖u‖ * (R * ‖w‖)) := mul_le_mul_of_nonneg_left h hβ.le
    _ = β * R ^ 2 * ‖u‖ * ‖w‖ := by ring

end fderiv

section gradient

variable [CompleteSpace E]

/-- The gradient of the log-sum-exp function is the softmax-weighted average of the `x i`. -/
lemma hasGradientAt_logSumExp (hβ : 0 < β) (x : ι → E) (v : E) :
    HasGradientAt (logSumExp β x) (∑ i, softmax β x v i • x i) v := by
  refine hasGradientAt_iff_hasFDerivAt.2 ((hasFDerivAt_logSumExp hβ.ne' x v).congr_fderiv ?_)
  ext w
  simp

/-- If `‖x i‖ ≤ R` for all `i`, the gradient of the log-sum-exp function has norm at most `R`. -/
lemma norm_gradient_logSumExp_le (hβ : 0 < β) {x : ι → E} {R : ℝ} (hR : ∀ i, ‖x i‖ ≤ R) (v : E) :
    ‖gradient (logSumExp β x) v‖ ≤ R := by
  rw [(hasGradientAt_logSumExp hβ x v).gradient]
  calc ‖∑ i, softmax β x v i • x i‖ ≤ ∑ i, ‖softmax β x v i • x i‖ := norm_sum_le _ _
    _ = ∑ i, softmax β x v i * ‖x i‖ := by
      simp_rw [norm_smul, Real.norm_of_nonneg (softmax_nonneg _ _ _ _)]
    _ ≤ ∑ i, softmax β x v i * R :=
      sum_le_sum fun i _ ↦ mul_le_mul_of_nonneg_left (hR i) (softmax_nonneg _ _ _ _)
    _ = R := by rw [← sum_mul, sum_softmax, one_mul]

end gradient
