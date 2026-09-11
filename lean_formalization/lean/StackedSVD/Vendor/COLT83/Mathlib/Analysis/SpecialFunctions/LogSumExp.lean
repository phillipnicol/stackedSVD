/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.Calculus.ContDiff.Operations
public import Mathlib.Data.Fintype.Order
public import Mathlib.Algebra.BigOperators.Field
public import Mathlib.Analysis.Calculus.FDeriv.Prod
public import Mathlib.Analysis.SpecialFunctions.ExpDeriv
public import Mathlib.Analysis.SpecialFunctions.Log.Deriv
public import StackedSVD.Vendor.COLT83.Mathlib.Algebra.Order.BigOperators.Covariance

set_option autoImplicit false

/-!
# The log-sum-exp function and the softmax weights

For `β : ℝ` and a vector `z : ι → ℝ` with `ι` finite, the *log-sum-exp* function is
`logSumExp β z = β⁻¹ log ∑ i, exp (β z i)` and the *softmax* weights are
`softmax β z i = exp (β z i) / ∑ j, exp (β z j)`, a probability vector (`softmax_nonneg`,
`sum_softmax`). For `β > 0` the log-sum-exp function is a smoothing of the maximum:

* `ciSup_le_logSumExp`, `logSumExp_le_ciSup_add`:
  `max_i z i ≤ logSumExp β z ≤ max_i z i + log |ι| / β`;
* `logSumExp_add_const`, `softmax_add_const`: it commutes with the addition of a constant, of
  which the softmax weights are independent;
* `logSumExp_mono`, `lipschitzWith_logSumExp`: it is monotone and `1`-Lipschitz for the sup norm;
* `hasFDerivAt_logSumExp`: it is smooth (`contDiff_logSumExp`), with partial derivatives the
  softmax weights, `∂ᵢ logSumExp β z = softmax β z i`, hence `‖D(logSumExp β) z‖ ≤ 1`
  (`norm_fderiv_logSumExp_le`);
* `fderiv_fderiv_logSumExp_apply`, `fderiv_fderiv_logSumExp_single`: its Hessian is
  `β (diag p - p pᵀ)` with `p` the softmax weights, the softmax-weighted covariance of the
  coordinates; its operator norm is at most `β` (`norm_fderiv_fderiv_logSumExp_le`).

The log-sum-exp smoothing of a maximum of linear forms `v ↦ max_i ⟪x i, v⟫` on an inner product
space is the composition of `logSumExp β` with `v ↦ (⟪x i, v⟫)ᵢ`, see
`COLT83.Mathlib.Analysis.InnerProductSpace.LogSumExp`.
-/

@[expose] public section

open Finset

namespace Real

variable {ι : Type*} [Fintype ι] {β : ℝ} {z z' : ι → ℝ}

/-- The coordinate projections `ι → ℝ` to `ℝ`, as continuous linear maps. -/
local notation "proj[" i "]" => (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : ι ↦ ℝ) i)

/-- The softmax weights `exp (β z i) / ∑ j, exp (β z j)` of a vector `z`. -/
noncomputable def softmax (β : ℝ) (z : ι → ℝ) (i : ι) : ℝ := exp (β * z i) / ∑ j, exp (β * z j)

/-- The log-sum-exp function `β⁻¹ log ∑ i, exp (β z i)` of a vector `z`. -/
noncomputable def logSumExp (β : ℝ) (z : ι → ℝ) : ℝ := β⁻¹ * log (∑ i, exp (β * z i))

/-- The softmax weights are nonnegative. -/
lemma softmax_nonneg (β : ℝ) (z : ι → ℝ) (i : ι) : 0 ≤ softmax β z i :=
  div_nonneg (exp_pos _).le (sum_nonneg fun _ _ ↦ (exp_pos _).le)

/-- The softmax weights are at most `1`. -/
lemma softmax_le_one (β : ℝ) (z : ι → ℝ) (i : ι) : softmax β z i ≤ 1 :=
  div_le_one_of_le₀ (single_le_sum (fun j _ ↦ (exp_pos (β * z j)).le) (mem_univ i))
    (sum_nonneg fun _ _ ↦ (exp_pos _).le)

/-- The softmax weights do not depend on the addition of a constant to the vector. -/
lemma softmax_add_const (β : ℝ) (z : ι → ℝ) (c : ℝ) :
    softmax β (fun i ↦ z i + c) = softmax β z := by
  ext i
  simp only [softmax, mul_add, exp_add, ← sum_mul]
  rw [mul_div_mul_right _ _ (exp_pos _).ne']

section Nonempty

variable [Nonempty ι]

/-- The denominator of the softmax weights is positive. -/
lemma sum_exp_mul_pos (β : ℝ) (z : ι → ℝ) : 0 < ∑ i, exp (β * z i) :=
  sum_pos (fun _ _ ↦ exp_pos _) univ_nonempty

/-- The softmax weights are positive. -/
lemma softmax_pos (β : ℝ) (z : ι → ℝ) (i : ι) : 0 < softmax β z i :=
  div_pos (exp_pos _) (sum_exp_mul_pos β z)

/-- The softmax weights sum to one. -/
lemma sum_softmax (β : ℝ) (z : ι → ℝ) : ∑ i, softmax β z i = 1 := by
  simp only [softmax, ← Finset.sum_div]
  exact div_self (sum_exp_mul_pos β z).ne'

/-- The log-sum-exp function commutes with the addition of a constant. -/
lemma logSumExp_add_const (hβ : β ≠ 0) (z : ι → ℝ) (c : ℝ) :
    logSumExp β (fun i ↦ z i + c) = logSumExp β z + c := by
  simp only [logSumExp, mul_add, exp_add, ← sum_mul]
  rw [log_mul (sum_exp_mul_pos β z).ne' (exp_pos _).ne', log_exp, mul_add, inv_mul_cancel_left₀ hβ]

/-- The log-sum-exp function of a constant vector. -/
lemma logSumExp_const (hβ : β ≠ 0) (c : ℝ) :
    logSumExp β (fun _ : ι ↦ c) = c + log (Fintype.card ι) / β := by
  have hcard : (Fintype.card ι : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Fintype.card_ne_zero
  rw [logSumExp, sum_const, card_univ, nsmul_eq_mul, log_mul hcard (exp_pos _).ne', log_exp]
  field_simp
  ring

/-- The log-sum-exp function dominates every coordinate, for `β > 0`. -/
lemma le_logSumExp (hβ : 0 < β) (z : ι → ℝ) (i : ι) : z i ≤ logSumExp β z := by
  rw [logSumExp, ← div_eq_inv_mul, le_div_iff₀ hβ, mul_comm,
    le_log_iff_exp_le (sum_exp_mul_pos β z)]
  exact single_le_sum (fun j _ ↦ (exp_pos (β * z j)).le) (mem_univ i)

/-- The log-sum-exp function dominates the maximum, for `β > 0`. -/
lemma ciSup_le_logSumExp (hβ : 0 < β) (z : ι → ℝ) : (⨆ i, z i) ≤ logSumExp β z :=
  ciSup_le (le_logSumExp hβ z)

/-- The log-sum-exp function exceeds the maximum by at most `log |ι| / β`, for `β > 0`. -/
lemma logSumExp_le_ciSup_add (hβ : 0 < β) (z : ι → ℝ) :
    logSumExp β z ≤ (⨆ i, z i) + log (Fintype.card ι) / β := by
  have hbdd : BddAbove (Set.range z) := Finite.bddAbove_range _
  have hS : ∑ i, exp (β * z i) ≤ Fintype.card ι * exp (β * ⨆ i, z i) := by
    have h := sum_le_card_nsmul univ (fun i ↦ exp (β * z i)) (exp (β * ⨆ i, z i)) fun i _ ↦
      exp_le_exp.2 (mul_le_mul_of_nonneg_left (le_ciSup hbdd i) hβ.le)
    simpa [nsmul_eq_mul] using h
  have hcard : (0 : ℝ) < Fintype.card ι := by exact_mod_cast Fintype.card_pos
  rw [logSumExp, ← div_eq_inv_mul, div_le_iff₀ hβ]
  calc log (∑ i, exp (β * z i)) ≤ log (Fintype.card ι * exp (β * ⨆ i, z i)) :=
        log_le_log (sum_exp_mul_pos β z) hS
    _ = β * (⨆ i, z i) + log (Fintype.card ι) := by
        rw [log_mul hcard.ne' (exp_pos _).ne', log_exp, add_comm]
    _ = ((⨆ i, z i) + log (Fintype.card ι) / β) * β := by
        rw [add_mul, div_mul_cancel₀ _ hβ.ne', mul_comm]

/-- The log-sum-exp function is monotone, for `β > 0`. -/
lemma logSumExp_mono (hβ : 0 < β) (h : z ≤ z') : logSumExp β z ≤ logSumExp β z' := by
  rw [logSumExp, logSumExp]
  refine mul_le_mul_of_nonneg_left (log_le_log (sum_exp_mul_pos β z) (sum_le_sum fun i _ ↦ ?_))
    (inv_nonneg.2 hβ.le)
  exact exp_le_exp.2 (mul_le_mul_of_nonneg_left (h i) hβ.le)

/-- The log-sum-exp function is `1`-Lipschitz for the sup norm, for `β > 0`. -/
lemma logSumExp_sub_le (hβ : 0 < β) (z z' : ι → ℝ) :
    logSumExp β z - logSumExp β z' ≤ ‖z - z'‖ := by
  have h : z ≤ fun i ↦ z' i + ‖z - z'‖ := fun i ↦ by
    have := norm_le_pi_norm (z - z') i
    rw [Pi.sub_apply, Real.norm_eq_abs] at this
    linarith [(abs_le.1 this).2]
  have := logSumExp_mono hβ h
  rw [logSumExp_add_const hβ.ne'] at this
  linarith

/-- The log-sum-exp function is `1`-Lipschitz for the sup norm, for `β > 0`. -/
lemma lipschitzWith_logSumExp (hβ : 0 < β) : LipschitzWith 1 (logSumExp β : (ι → ℝ) → ℝ) := by
  refine LipschitzWith.of_dist_le_mul fun z z' ↦ ?_
  rw [NNReal.coe_one, one_mul, Real.dist_eq, dist_eq_norm]
  refine abs_sub_le_iff.2 ⟨logSumExp_sub_le hβ z z', ?_⟩
  simpa only [norm_sub_rev] using logSumExp_sub_le hβ z' z

/-! ### Derivatives -/

/-- The log-sum-exp function is smooth. -/
lemma contDiff_logSumExp (n : WithTop ℕ∞) (β : ℝ) :
    ContDiff ℝ n (logSumExp β : (ι → ℝ) → ℝ) := by
  have h : ContDiff ℝ n fun z : ι → ℝ ↦ ∑ i, exp (β * z i) :=
    ContDiff.sum fun i _ ↦ (contDiff_const.mul (contDiff_apply ℝ ℝ i)).exp
  exact contDiff_const.mul (h.log fun z ↦ (sum_exp_mul_pos β z).ne')

omit [Nonempty ι] in
/-- The Fréchet derivative of `z ↦ ∑ i, exp (β z i)`. -/
lemma hasFDerivAt_sum_exp_mul (β : ℝ) (z : ι → ℝ) :
    HasFDerivAt (fun z : ι → ℝ ↦ ∑ i, exp (β * z i))
      (∑ i, (exp (β * z i) * β) • proj[i]) z := by
  refine HasFDerivAt.fun_sum fun i _ ↦ ?_
  have h : HasFDerivAt (fun z : ι → ℝ ↦ β * z i) (β • proj[i]) z :=
    (hasFDerivAt_apply (𝕜 := ℝ) i z).const_mul β
  simpa only [smul_smul] using h.exp

/-- The Fréchet derivative of the log-sum-exp function: `∂ᵢ logSumExp β z = softmax β z i`. -/
lemma hasFDerivAt_logSumExp (hβ : β ≠ 0) (z : ι → ℝ) :
    HasFDerivAt (logSumExp β) (∑ i, softmax β z i • proj[i]) z := by
  have h := ((hasFDerivAt_sum_exp_mul β z).log (sum_exp_mul_pos β z).ne').const_mul β⁻¹
  refine h.congr_fderiv ?_
  ext w
  simp only [smul_apply, FunLike.coe_sum, Finset.sum_apply, ContinuousLinearMap.proj_apply,
    softmax, smul_eq_mul, mul_sum]
  refine sum_congr rfl fun i _ ↦ ?_
  field_simp

/-- The Fréchet derivative of the log-sum-exp function. -/
lemma fderiv_logSumExp (hβ : β ≠ 0) (z : ι → ℝ) :
    fderiv ℝ (logSumExp β) z = ∑ i, softmax β z i • proj[i] :=
  (hasFDerivAt_logSumExp hβ z).fderiv

/-- The operator norm of the derivative of the log-sum-exp function is at most `1` (for the
sup norm on `ι → ℝ`). -/
lemma norm_fderiv_logSumExp_le (hβ : β ≠ 0) (z : ι → ℝ) : ‖fderiv ℝ (logSumExp β) z‖ ≤ 1 := by
  rw [fderiv_logSumExp hβ]
  refine ContinuousLinearMap.opNorm_le_bound _ zero_le_one fun w ↦ ?_
  simp only [FunLike.coe_sum, Finset.sum_apply, smul_apply, ContinuousLinearMap.proj_apply,
    smul_eq_mul, Real.norm_eq_abs, one_mul]
  calc |∑ i, softmax β z i * w i| ≤ ∑ i, |softmax β z i * w i| := abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, softmax β z i * ‖w‖ := sum_le_sum fun i _ ↦ by
        rw [abs_mul, abs_of_nonneg (softmax_nonneg β z i)]
        exact mul_le_mul_of_nonneg_left ((Real.norm_eq_abs _).symm.trans_le (norm_le_pi_norm w i))
          (softmax_nonneg β z i)
    _ = ‖w‖ := by rw [← sum_mul, sum_softmax, one_mul]

/-- The Fréchet derivative of the softmax weights:
`D(softmax β · i) z = β pᵢ (eᵢ - ∑ j, pⱼ eⱼ)` with `p = softmax β z`. -/
lemma hasFDerivAt_softmax (β : ℝ) (z : ι → ℝ) (i : ι) :
    HasFDerivAt (fun z ↦ softmax β z i)
      ((β * softmax β z i) •
        (proj[i] - ∑ j, softmax β z j • proj[j])) z := by
  have hS := sum_exp_mul_pos β z
  have hnum : HasFDerivAt (fun z : ι → ℝ ↦ exp (β * z i))
      ((exp (β * z i) * β) • proj[i]) z := by
    have h : HasFDerivAt (fun z : ι → ℝ ↦ β * z i) (β • proj[i]) z :=
      (hasFDerivAt_apply (𝕜 := ℝ) i z).const_mul β
    simpa only [smul_smul] using h.exp
  have hden : HasFDerivAt (fun z : ι → ℝ ↦ (∑ j, exp (β * z j))⁻¹)
      ((-((∑ j, exp (β * z j)) ^ 2)⁻¹) • ∑ j, (exp (β * z j) * β) • proj[j])
      z :=
    (hasDerivAt_inv hS.ne').comp_hasFDerivAt z (hasFDerivAt_sum_exp_mul β z)
  simp_rw [softmax, div_eq_mul_inv]
  refine (hnum.fun_mul hden).congr_fderiv ?_
  ext w
  simp only [smul_apply, add_apply, sub_apply, FunLike.coe_sum, Finset.sum_apply,
    ContinuousLinearMap.proj_apply, smul_eq_mul]
  have h1 : ∑ j, exp (β * z j) * β * w j = β * ∑ j, exp (β * z j) * w j := by
    rw [mul_sum]
    exact sum_congr rfl fun j _ ↦ by ring
  have h2 : ∑ j, exp (β * z j) * (∑ j, exp (β * z j))⁻¹ * w j
      = (∑ j, exp (β * z j))⁻¹ * ∑ j, exp (β * z j) * w j := by
    rw [mul_sum]
    exact sum_congr rfl fun j _ ↦ by ring
  rw [h1, h2]
  field_simp
  ring

/-- The second derivative of the log-sum-exp function, as a bilinear form:
`D²(logSumExp β) z (u, w) = β (∑ i, pᵢ uᵢ wᵢ - (∑ i, pᵢ uᵢ) (∑ i, pᵢ wᵢ))`. -/
lemma hasFDerivAt_fderiv_logSumExp (hβ : β ≠ 0) (z : ι → ℝ) :
    HasFDerivAt (fderiv ℝ (logSumExp β))
      (β • (∑ i, softmax β z i •
          (proj[i]).smulRight (proj[i]) -
        (∑ i, softmax β z i • proj[i]).smulRight
          (∑ i, softmax β z i • proj[i]))) z := by
  have h : fderiv ℝ (logSumExp β) = fun z : ι → ℝ ↦ ∑ i, softmax β z i • proj[i] :=
    funext (fderiv_logSumExp hβ)
  rw [h]
  refine (HasFDerivAt.fun_sum fun i _ ↦
    (hasFDerivAt_softmax β z i).smul_const (proj[i])).congr_fderiv ?_
  ext u w
  simp only [smul_apply, sub_apply, FunLike.coe_sum, Finset.sum_apply,
    ContinuousLinearMap.proj_apply, ContinuousLinearMap.smulRight_apply, smul_eq_mul]
  rw [mul_sub, mul_sum, mul_sum, mul_sum, ← sum_sub_distrib]
  exact sum_congr rfl fun i _ ↦ by ring

/-- The second derivative of the log-sum-exp function, applied to two vectors. -/
lemma fderiv_fderiv_logSumExp_apply (hβ : β ≠ 0) (z u w : ι → ℝ) :
    fderiv ℝ (fderiv ℝ (logSumExp β)) z u w =
      β * (∑ i, softmax β z i * (u i * w i) -
        (∑ i, softmax β z i * u i) * (∑ i, softmax β z i * w i)) := by
  rw [(hasFDerivAt_fderiv_logSumExp hβ z).fderiv]
  simp only [smul_apply, sub_apply, FunLike.coe_sum, Finset.sum_apply,
    ContinuousLinearMap.proj_apply, ContinuousLinearMap.smulRight_apply, smul_eq_mul]

/-- The Hessian of the log-sum-exp function is `β (δᵢⱼ pᵢ - pᵢ pⱼ)` with `p = softmax β z`. -/
lemma fderiv_fderiv_logSumExp_single [DecidableEq ι] (hβ : β ≠ 0) (z : ι → ℝ) (j l : ι) :
    fderiv ℝ (fderiv ℝ (logSumExp β)) z (Pi.single j 1) (Pi.single l 1) =
      β * ((if j = l then softmax β z j else 0) - softmax β z j * softmax β z l) := by
  rw [fderiv_fderiv_logSumExp_apply hβ]
  simp only [Pi.single_apply, mul_ite, mul_one, mul_zero, sum_ite_eq', mem_univ, ite_true]
  by_cases hjl : j = l
  · subst hjl
    simp
  · simp [hjl, Ne.symm hjl]

/-- The operator norm of the Hessian of the log-sum-exp function is at most `β` (for the sup
norm on `ι → ℝ`), for `β > 0`. -/
lemma norm_fderiv_fderiv_logSumExp_le (hβ : 0 < β) (z : ι → ℝ) :
    ‖fderiv ℝ (fderiv ℝ (logSumExp β)) z‖ ≤ β := by
  refine ContinuousLinearMap.opNorm_le_bound₂ _ hβ.le fun u w ↦ ?_
  rw [fderiv_fderiv_logSumExp_apply hβ.ne' z u w, Real.norm_eq_abs, abs_mul, abs_of_pos hβ]
  have hbound (y : ι → ℝ) (i : ι) : |y i| ≤ ‖y‖ :=
    (Real.norm_eq_abs _).symm.trans_le (norm_le_pi_norm y i)
  have h := abs_sum_mul_sub_mul_le (s := univ) (fun i _ ↦ softmax_nonneg β z i)
    (sum_softmax β z) (fun i _ ↦ hbound u i) (fun i _ ↦ hbound w i)
  calc β * |∑ i, softmax β z i * (u i * w i) -
        (∑ i, softmax β z i * u i) * (∑ i, softmax β z i * w i)|
      ≤ β * (‖u‖ * ‖w‖) := mul_le_mul_of_nonneg_left h hβ.le
    _ = β * ‖u‖ * ‖w‖ := by ring

end Nonempty

end Real
