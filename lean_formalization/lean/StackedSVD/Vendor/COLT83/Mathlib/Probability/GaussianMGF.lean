/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Probability.Distributions.Gaussian.Multivariate
public import Mathlib.Probability.Distributions.Gaussian.Fernique
public import Mathlib.Probability.Moments.SubGaussian
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.Calculus.Gradient

set_option autoImplicit false

/-!
# Exponential moments of Gaussian measures

* `integral_exp_inner_stdGaussian`: the moment generating function of a linear form
  `⟪a, ·⟫` under the standard Gaussian measure is `exp (t² ‖a‖² / 2)`.
* `IsGaussian.integrable_exp_mul_norm`: `exp (c ‖x‖)` is integrable for every Gaussian measure
  on a separable Banach space (a consequence of Fernique's theorem), and therefore `exp (t F x)`
  is integrable for every function `F` of linear growth
  (`IsGaussian.integrable_exp_of_abs_le_add_mul_norm`), and functions of linear or quadratic
  growth are integrable (`IsGaussian.integrable_of_abs_le_add_mul_norm`,
  `IsGaussian.integrable_of_abs_le_add_mul_norm_sq`).
-/

@[expose] public section

open MeasureTheory ProbabilityTheory Real InnerProductSpace
open scoped RealInnerProductSpace NNReal

namespace ProbabilityTheory

section stdGaussian

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  [MeasurableSpace E] [BorelSpace E]

/-- The law of the linear form `⟪a, ·⟫` under the standard Gaussian measure is `N(0, ‖a‖²)`. -/
lemma stdGaussian_map_toDual (a : E) :
    (stdGaussian E).map (toDual ℝ E a) = gaussianReal 0 (‖a‖₊ ^ 2) := by
  rw [IsGaussian.map_eq_gaussianReal (toDual ℝ E a), integral_strongDual_stdGaussian,
    variance_dual_stdGaussian, LinearIsometryEquiv.norm_map]
  congr
  rw [← coe_nnnorm, ← NNReal.coe_pow, Real.toNNReal_coe]

lemma integrable_exp_mul_inner_stdGaussian (a : E) (t : ℝ) :
    Integrable (fun x ↦ exp (t * ⟪a, x⟫)) (stdGaussian E) := by
  have h := integrable_exp_mul_gaussianReal (μ := 0) (v := ‖a‖₊ ^ 2) t
  rw [← stdGaussian_map_toDual a, integrable_map_measure (by fun_prop) (by fun_prop)] at h
  exact h

/-- The moment generating function of `⟪a, ·⟫` under the standard Gaussian measure. -/
lemma integral_exp_mul_inner_stdGaussian (a : E) (t : ℝ) :
    ∫ x, exp (t * ⟪a, x⟫) ∂stdGaussian E = exp (t ^ 2 * ‖a‖ ^ 2 / 2) := by
  have h : ∫ x, exp (t * ⟪a, x⟫) ∂stdGaussian E =
      exp (0 * t + ((‖a‖₊ ^ 2 : ℝ≥0) : ℝ) * t ^ 2 / 2) :=
    mgf_gaussianReal (stdGaussian_map_toDual a) t
  rw [h]
  push_cast
  ring_nf

lemma integral_exp_inner_stdGaussian (a : E) :
    ∫ x, exp ⟪a, x⟫ ∂stdGaussian E = exp (‖a‖ ^ 2 / 2) := by
  simpa using integral_exp_mul_inner_stdGaussian a 1

end stdGaussian

section isGaussian

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E] {μ : Measure E} [IsGaussian μ]

/-- **Exponential moments of Gaussian norms**: `exp (c ‖x‖)` is integrable for every Gaussian
measure (Fernique's theorem). -/
lemma IsGaussian.integrable_exp_mul_norm (c : ℝ) : Integrable (fun x ↦ exp (c * ‖x‖)) μ := by
  obtain ⟨C, hC, hint⟩ := IsGaussian.exists_integrable_exp_sq μ
  refine (hint.const_mul (exp (c ^ 2 / (4 * C)))).mono'
    (continuous_const.mul continuous_norm).rexp.aestronglyMeasurable
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_of_nonneg (exp_pos _).le, ← Real.exp_add]
  refine Real.exp_le_exp.2 ?_
  rw [div_add' _ _ _ (by positivity), le_div_iff₀ (by positivity)]
  nlinarith [sq_nonneg (2 * C * ‖x‖ - c)]

/-- A measurable function of linear growth has integrable exponential moments under a Gaussian
measure. -/
lemma IsGaussian.integrable_exp_of_abs_le_add_mul_norm {F : E → ℝ} (hF : AEStronglyMeasurable F μ)
    {A B : ℝ} (hFle : ∀ x, |F x| ≤ A + B * ‖x‖) (t : ℝ) :
    Integrable (fun x ↦ exp (t * F x)) μ := by
  refine ((IsGaussian.integrable_exp_mul_norm (|t| * B)).const_mul (exp (|t| * A))).mono'
    (Real.continuous_exp.comp_aestronglyMeasurable (hF.const_mul t))
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_of_nonneg (exp_pos _).le, ← Real.exp_add]
  refine Real.exp_le_exp.2 ?_
  calc t * F x ≤ |t * F x| := le_abs_self _
    _ = |t| * |F x| := abs_mul _ _
    _ ≤ |t| * (A + B * ‖x‖) := mul_le_mul_of_nonneg_left (hFle x) (abs_nonneg t)
    _ = |t| * A + |t| * B * ‖x‖ := by ring

/-- A measurable function of linear growth is integrable under a Gaussian measure. -/
lemma IsGaussian.integrable_of_abs_le_add_mul_norm {F : E → ℝ} (hF : AEStronglyMeasurable F μ)
    {A B : ℝ} (hFle : ∀ x, |F x| ≤ A + B * ‖x‖) : Integrable F μ := by
  refine ((integrable_const A).add (IsGaussian.integrable_id.norm.const_mul B)).mono' hF
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_eq_abs]
  exact hFle x

/-- A Gaussian measure has finite moments of all orders. -/
lemma IsGaussian.integrable_norm_pow (n : ℕ) : Integrable (fun x ↦ ‖x‖ ^ n) μ :=
  (IsGaussian.memLp_id μ n (by simp)).integrable_norm_pow'

/-- A measurable function of quadratic growth is integrable under a Gaussian measure. -/
lemma IsGaussian.integrable_of_abs_le_add_mul_norm_sq {F : E → ℝ} (hF : AEStronglyMeasurable F μ)
    {A B C : ℝ} (hFle : ∀ x, |F x| ≤ A + B * ‖x‖ + C * ‖x‖ ^ 2) : Integrable F μ := by
  refine (((integrable_const A).add (IsGaussian.integrable_id.norm.const_mul B)).add
    ((IsGaussian.integrable_norm_pow 2).const_mul C)).mono' hF
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_eq_abs]
  exact hFle x

end isGaussian

section gradient

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E] {f : E → ℝ} {L : ℝ≥0}

/-- A `C¹` function with bounded gradient is integrable under a Gaussian measure, with all its
exponential moments. -/
lemma IsGaussian.integrable_exp_of_norm_gradient_le {μ : Measure E} [IsGaussian μ]
    (hf : ContDiff ℝ 1 f) (hL : ∀ x, ‖gradient f x‖ ≤ L) (t : ℝ) :
    Integrable (fun x ↦ exp (t * f x)) μ :=
  IsGaussian.integrable_exp_of_abs_le_add_mul_norm hf.continuous.aestronglyMeasurable
    (abs_le_of_norm_gradient_le (hf.differentiable one_ne_zero) hL) t

lemma IsGaussian.integrable_of_norm_gradient_le {μ : Measure E} [IsGaussian μ]
    (hf : ContDiff ℝ 1 f) (hL : ∀ x, ‖gradient f x‖ ≤ L) : Integrable f μ :=
  IsGaussian.integrable_of_abs_le_add_mul_norm hf.continuous.aestronglyMeasurable
    (abs_le_of_norm_gradient_le (hf.differentiable one_ne_zero) hL)

end gradient

section real

/-- A centered real Gaussian with variance `v` is sub-Gaussian with variance proxy `v`. -/
lemma hasSubgaussianMGF_fun_id_gaussianReal (v : ℝ≥0) :
    HasSubgaussianMGF (fun x ↦ x) v (gaussianReal 0 v) where
  integrable_exp_mul t := integrable_exp_mul_gaussianReal t
  mgf_le t := by
    rw [mgf_fun_id_gaussianReal]
    simp

end real

end ProbabilityTheory
