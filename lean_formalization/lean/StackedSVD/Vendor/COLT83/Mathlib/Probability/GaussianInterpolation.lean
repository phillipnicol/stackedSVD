/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinIdentity
public import Mathlib.Analysis.Calculus.ParametricIntegral
public import Mathlib.Analysis.SpecialFunctions.Sqrt

set_option autoImplicit false

/-!
# Gaussian interpolation

For linear maps `L₀ L₁ : E' → E` and `t ∈ [0, 1]`, `gaussianInterp L₀ L₁ t = √t L₁ + √(1-t) L₀`.
If `g ~ N(0, I)` on `E'`, then `gaussianInterp L₀ L₁ t g` interpolates between the Gaussian vectors
`L₀ g` (at `t = 0`) and `L₁ g` (at `t = 1`); when `L₀ g` and `L₁ g` are independent, its law is
the Gaussian measure with covariance `t Cov(L₁ g) + (1 - t) Cov(L₀ g)`.

The main result is the **Gaussian interpolation formula**
`hasDerivAt_integral_comp_gaussianInterp`: for a `C²` function `F` with bounded second derivative,
`ψ(t) := E[F(gaussianInterp L₀ L₁ t g)]` is differentiable on `(0, 1)` with
`ψ'(t) = E[∑ k, D²F(Z_t g) (Z_t (b k)) (Z_t' (b k))]` for an orthonormal basis `b` of `E'`,
where `Z_t = gaussianInterp L₀ L₁ t` and `Z_t' = gaussianInterpDeriv L₀ L₁ t` is its derivative in
`t`. It is obtained by differentiating under the integral sign and Stein's identity
(`integral_fderiv_apply_stdGaussian`). `ψ` is moreover continuous on `[0, 1]`
(`continuousOn_integral_comp_gaussianInterp`). This is the key step of the Sudakov–Fernique
inequality.
-/

@[expose] public section

open MeasureTheory Real InnerProductSpace
open scoped RealInnerProductSpace

namespace ProbabilityTheory

variable {E E' : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [NormedAddCommGroup E']
  [InnerProductSpace ℝ E']

section interp

variable (L₀ L₁ : E' →L[ℝ] E) {t : ℝ}

/-- The interpolation `√t L₁ + √(1 - t) L₀` between two linear maps. -/
noncomputable def gaussianInterp (t : ℝ) : E' →L[ℝ] E := √t • L₁ + √(1 - t) • L₀

/-- The derivative in `t` of `gaussianInterp L₀ L₁ t`, for `t ∈ (0, 1)`. -/
noncomputable def gaussianInterpDeriv (t : ℝ) : E' →L[ℝ] E :=
  (2 * √t)⁻¹ • L₁ - (2 * √(1 - t))⁻¹ • L₀

lemma gaussianInterp_apply (t : ℝ) (g : E') :
    gaussianInterp L₀ L₁ t g = √t • L₁ g + √(1 - t) • L₀ g := rfl

lemma gaussianInterpDeriv_apply (t : ℝ) (g : E') :
    gaussianInterpDeriv L₀ L₁ t g = (2 * √t)⁻¹ • L₁ g - (2 * √(1 - t))⁻¹ • L₀ g := rfl

@[simp] lemma gaussianInterp_zero : gaussianInterp L₀ L₁ 0 = L₀ := by
  simp [gaussianInterp]

@[simp] lemma gaussianInterp_one : gaussianInterp L₀ L₁ 1 = L₁ := by
  simp [gaussianInterp]

lemma continuous_gaussianInterp_apply (g : E') :
    Continuous fun t ↦ gaussianInterp L₀ L₁ t g := by
  simp_rw [gaussianInterp_apply]
  fun_prop

lemma hasDerivAt_gaussianInterp_apply (ht : t ∈ Set.Ioo 0 1) (g : E') :
    HasDerivAt (fun t ↦ gaussianInterp L₀ L₁ t g) (gaussianInterpDeriv L₀ L₁ t g) t := by
  simp_rw [gaussianInterp_apply, gaussianInterpDeriv_apply]
  have h1 : HasDerivAt (fun t : ℝ ↦ √t) (1 / (2 * √t)) t := Real.hasDerivAt_sqrt ht.1.ne'
  have h2 : HasDerivAt (fun t : ℝ ↦ √(1 - t)) (1 / (2 * √(1 - t)) * (-1)) t :=
    (Real.hasDerivAt_sqrt (sub_pos.2 ht.2).ne').comp t ((hasDerivAt_id t).const_sub 1)
  refine ((h1.smul_const (L₁ g)).add (h2.smul_const (L₀ g))).congr_deriv ?_
  simp [sub_eq_add_neg, neg_smul, one_div]

lemma norm_gaussianInterp_apply_le (ht : t ∈ Set.Icc 0 1) (g : E') :
    ‖gaussianInterp L₀ L₁ t g‖ ≤ (‖L₁‖ + ‖L₀‖) * ‖g‖ := by
  have h1 : √t ≤ 1 := Real.sqrt_le_one.2 ht.2
  have h2 : √(1 - t) ≤ 1 := Real.sqrt_le_one.2 (by linarith [ht.1])
  rw [gaussianInterp_apply]
  calc ‖√t • L₁ g + √(1 - t) • L₀ g‖ ≤ ‖√t • L₁ g‖ + ‖√(1 - t) • L₀ g‖ := norm_add_le _ _
    _ = √t * ‖L₁ g‖ + √(1 - t) * ‖L₀ g‖ := by
        rw [norm_smul, norm_smul, Real.norm_of_nonneg (Real.sqrt_nonneg _),
          Real.norm_of_nonneg (Real.sqrt_nonneg _)]
    _ ≤ 1 * (‖L₁‖ * ‖g‖) + 1 * (‖L₀‖ * ‖g‖) := by
        gcongr
        · exact L₁.le_opNorm g
        · exact L₀.le_opNorm g
    _ = (‖L₁‖ + ‖L₀‖) * ‖g‖ := by ring

lemma norm_gaussianInterpDeriv_apply_le {a b : ℝ} (ha : 0 < a) (hb : b < 1)
    (ht : t ∈ Set.Icc a b) (g : E') :
    ‖gaussianInterpDeriv L₀ L₁ t g‖ ≤ ((2 * √a)⁻¹ * ‖L₁‖ + (2 * √(1 - b))⁻¹ * ‖L₀‖) * ‖g‖ := by
  have ha' : 0 < √a := Real.sqrt_pos.2 ha
  have hb' : 0 < √(1 - b) := Real.sqrt_pos.2 (by linarith)
  have h1 : (2 * √t)⁻¹ ≤ (2 * √a)⁻¹ := by
    gcongr
    exact ht.1
  have h2 : (2 * √(1 - t))⁻¹ ≤ (2 * √(1 - b))⁻¹ := by
    gcongr
    exact ht.2
  rw [gaussianInterpDeriv_apply]
  calc ‖(2 * √t)⁻¹ • L₁ g - (2 * √(1 - t))⁻¹ • L₀ g‖
      ≤ ‖(2 * √t)⁻¹ • L₁ g‖ + ‖(2 * √(1 - t))⁻¹ • L₀ g‖ := norm_sub_le _ _
    _ = (2 * √t)⁻¹ * ‖L₁ g‖ + (2 * √(1 - t))⁻¹ * ‖L₀ g‖ := by
        rw [norm_smul, norm_smul, Real.norm_of_nonneg (by positivity),
          Real.norm_of_nonneg (by positivity)]
    _ ≤ (2 * √a)⁻¹ * (‖L₁‖ * ‖g‖) + (2 * √(1 - b))⁻¹ * (‖L₀‖ * ‖g‖) := by
        gcongr
        · exact L₁.le_opNorm g
        · exact L₀.le_opNorm g
    _ = ((2 * √a)⁻¹ * ‖L₁‖ + (2 * √(1 - b))⁻¹ * ‖L₀‖) * ‖g‖ := by ring

end interp

section integral

variable [FiniteDimensional ℝ E'] [MeasurableSpace E'] [BorelSpace E'] {L₀ L₁ : E' →L[ℝ] E}
  {F : E → ℝ} {K : ℝ} {t : ℝ}

omit [FiniteDimensional ℝ E'] [MeasurableSpace E'] [BorelSpace E'] in
/-- Quadratic domination of `F ∘ gaussianInterp L₀ L₁ t` for `t ∈ [0, 1]`. -/
lemma abs_comp_gaussianInterp_le (hF : ContDiff ℝ 2 F) (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K)
    (ht : t ∈ Set.Icc 0 1) (g : E') :
    |F (gaussianInterp L₀ L₁ t g)| ≤ |F 0| + ‖fderiv ℝ F 0‖ * (‖L₁‖ + ‖L₀‖) * ‖g‖ +
      K * (‖L₁‖ + ‖L₀‖) ^ 2 * ‖g‖ ^ 2 := by
  have hK0 : 0 ≤ K := le_trans (norm_nonneg (fderiv ℝ (fderiv ℝ F) 0)) (hK 0)
  have hz := norm_gaussianInterp_apply_le L₀ L₁ ht g
  have h := norm_le_of_norm_fderiv_fderiv_le hF hK (gaussianInterp L₀ L₁ t g)
  rw [Real.norm_eq_abs, Real.norm_eq_abs] at h
  refine h.trans ?_
  have hC : 0 ≤ (‖L₁‖ + ‖L₀‖) * ‖g‖ := by positivity
  calc |F 0| + (‖fderiv ℝ F 0‖ + K * ‖gaussianInterp L₀ L₁ t g‖) * ‖gaussianInterp L₀ L₁ t g‖
      ≤ |F 0| + (‖fderiv ℝ F 0‖ + K * ((‖L₁‖ + ‖L₀‖) * ‖g‖)) * ((‖L₁‖ + ‖L₀‖) * ‖g‖) := by
        gcongr
    _ = _ := by ring

/-- `F ∘ gaussianInterp L₀ L₁ t` is integrable under the standard Gaussian measure for
`t ∈ [0, 1]`. -/
lemma integrable_comp_gaussianInterp (hF : ContDiff ℝ 2 F)
    (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K) (ht : t ∈ Set.Icc 0 1) :
    Integrable (fun g ↦ F (gaussianInterp L₀ L₁ t g)) (stdGaussian E') :=
  IsGaussian.integrable_of_abs_le_add_mul_norm_sq
    (hF.continuous.comp (gaussianInterp L₀ L₁ t).continuous).aestronglyMeasurable
    (abs_comp_gaussianInterp_le hF hK ht)

/-- `t ↦ E[F(gaussianInterp L₀ L₁ t g)]` is continuous on `[0, 1]`. -/
lemma continuousOn_integral_comp_gaussianInterp (hF : ContDiff ℝ 2 F)
    (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K) :
    ContinuousOn (fun t ↦ ∫ g, F (gaussianInterp L₀ L₁ t g) ∂stdGaussian E') (Set.Icc 0 1) := by
  refine continuousOn_of_dominated (bound := fun g ↦ |F 0| + ‖fderiv ℝ F 0‖ * (‖L₁‖ + ‖L₀‖) * ‖g‖ +
    K * (‖L₁‖ + ‖L₀‖) ^ 2 * ‖g‖ ^ 2) (fun t ht ↦ (integrable_comp_gaussianInterp hF hK ht).1)
    (fun t ht ↦ Filter.Eventually.of_forall fun g ↦ ?_) ?_
    (Filter.Eventually.of_forall fun g ↦ ?_)
  · rw [Real.norm_eq_abs]
    exact abs_comp_gaussianInterp_le hF hK ht g
  · exact ((integrable_const _).add (IsGaussian.integrable_id.norm.const_mul _)).add
      ((IsGaussian.integrable_norm_pow 2).const_mul _)
  · exact (hF.continuous.comp (continuous_gaussianInterp_apply L₀ L₁ g)).continuousOn

variable {κ : Type*} [Fintype κ]

/-- **Gaussian interpolation formula**: for `t ∈ (0, 1)`, `t ↦ E[F(Z_t g)]` with
`Z_t = gaussianInterp L₀ L₁ t` has derivative `E[∑ k, D²F(Z_t g) (Z_t (b k)) (Z_t' (b k))]`, where
`Z_t' = gaussianInterpDeriv L₀ L₁ t` and `b` is any orthonormal basis of `E'`. -/
lemma hasDerivAt_integral_comp_gaussianInterp (hF : ContDiff ℝ 2 F)
    (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K) (b : OrthonormalBasis κ ℝ E') (ht : t ∈ Set.Ioo 0 1) :
    HasDerivAt (fun t ↦ ∫ g, F (gaussianInterp L₀ L₁ t g) ∂stdGaussian E')
      (∫ g, ∑ k, fderiv ℝ (fderiv ℝ F) (gaussianInterp L₀ L₁ t g) (gaussianInterp L₀ L₁ t (b k))
        (gaussianInterpDeriv L₀ L₁ t (b k)) ∂stdGaussian E') t := by
  have hK0 : 0 ≤ K := le_trans (norm_nonneg (fderiv ℝ (fderiv ℝ F) 0)) (hK 0)
  have hF1 : Differentiable ℝ F := hF.differentiable (by norm_num)
  set a := t / 2 with ha
  set b' := (1 + t) / 2 with hb'
  have ha0 : 0 < a := by rw [ha]; linarith [ht.1]
  have hb1 : b' < 1 := by rw [hb']; linarith [ht.2]
  have hsub : Set.Icc a b' ⊆ Set.Ioo 0 1 := fun s hs ↦ ⟨ha0.trans_le hs.1, hs.2.trans_lt hb1⟩
  set C := ‖L₁‖ + ‖L₀‖ with hC
  set C' := (2 * √a)⁻¹ * ‖L₁‖ + (2 * √(1 - b'))⁻¹ * ‖L₀‖ with hC'
  have hC0 : 0 ≤ C := by positivity
  have hC'0 : 0 ≤ C' := by positivity
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := stdGaussian E')
    (F := fun t g ↦ F (gaussianInterp L₀ L₁ t g))
    (F' := fun t g ↦ fderiv ℝ F (gaussianInterp L₀ L₁ t g) (gaussianInterpDeriv L₀ L₁ t g))
    (bound := fun g ↦ ‖fderiv ℝ F 0‖ * C' * ‖g‖ + K * C * C' * ‖g‖ ^ 2)
    (Icc_mem_nhds (by rw [ha]; linarith [ht.1]) (by rw [hb']; linarith [ht.2]))
    (Filter.Eventually.of_forall fun s ↦
      (hF.continuous.comp (gaussianInterp L₀ L₁ s).continuous).aestronglyMeasurable)
    (integrable_comp_gaussianInterp hF hK (Set.Ioo_subset_Icc_self ht))
    (((hF.continuous_fderiv (by norm_num)).comp (gaussianInterp L₀ L₁ t).continuous
      |>.clm_apply (gaussianInterpDeriv L₀ L₁ t).continuous).aestronglyMeasurable)
    (Filter.Eventually.of_forall fun g s hs ↦ ?_)
    (((IsGaussian.integrable_id.norm.const_mul _)).add
      ((IsGaussian.integrable_norm_pow 2).const_mul _))
    (Filter.Eventually.of_forall fun g s hs ↦
      (hF1 _).hasFDerivAt.comp_hasDerivAt s (hasDerivAt_gaussianInterp_apply L₀ L₁ (hsub hs) g))
  · refine key.2.congr_deriv ?_
    exact integral_fderiv_apply_stdGaussian hF hK (gaussianInterp L₀ L₁ t)
      (gaussianInterpDeriv L₀ L₁ t) b
  · have h1 := norm_fderiv_le_of_norm_fderiv_fderiv_le hF hK (gaussianInterp L₀ L₁ s g)
    have h2 := norm_gaussianInterp_apply_le L₀ L₁ (Set.Ioo_subset_Icc_self (hsub hs)) g
    have h3 := norm_gaussianInterpDeriv_apply_le L₀ L₁ ha0 hb1 hs g
    calc ‖fderiv ℝ F (gaussianInterp L₀ L₁ s g) (gaussianInterpDeriv L₀ L₁ s g)‖
        ≤ ‖fderiv ℝ F (gaussianInterp L₀ L₁ s g)‖ * ‖gaussianInterpDeriv L₀ L₁ s g‖ :=
          ContinuousLinearMap.le_opNorm _ _
      _ ≤ (‖fderiv ℝ F 0‖ + K * (C * ‖g‖)) * (C' * ‖g‖) := by
          gcongr
          exact h1.trans (by gcongr)
      _ = ‖fderiv ℝ F 0‖ * C' * ‖g‖ + K * C * C' * ‖g‖ ^ 2 := by ring

end integral

end ProbabilityTheory
