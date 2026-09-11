/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.Calculus.Gradient
public import Mathlib.Analysis.Convex.Integral
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv

set_option autoImplicit false

/-!
# Interpolation along a quarter circle

For a `C¹` function `f` on a real inner product space,
`f u - f v = ∫₀^{π/2} ⟪∇f (sin θ • u + cos θ • v), cos θ • u - sin θ • v⟫ dθ`
(`sub_eq_integral_quarterCircle`), together with the derivative of `θ ↦ f (sin θ • u + cos θ • v)`
and Jensen's inequality for the uniform measure on `[0, π/2]`
(`exp_mul_integral_le_integral_exp`). These are the ingredients of the Maurey–Pisier Gaussian
concentration argument.
-/

@[expose] public section

open MeasureTheory Real InnerProductSpace Set
open scoped RealInnerProductSpace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]
  {f : E → ℝ}

/-- The derivative of `θ ↦ f (sin θ • u + cos θ • v)`. -/
lemma hasDerivAt_comp_quarterCircle (hf : Differentiable ℝ f) (u v : E) (θ : ℝ) :
    HasDerivAt (fun θ ↦ f (sin θ • u + cos θ • v))
      ⟪gradient f (sin θ • u + cos θ • v), cos θ • u - sin θ • v⟫ θ := by
  have hγ : HasDerivAt (fun θ ↦ sin θ • u + cos θ • v) (cos θ • u + (-sin θ) • v) θ :=
    ((hasDerivAt_sin θ).smul_const u).add ((hasDerivAt_cos θ).smul_const v)
  have h := (hasGradientAt_iff_hasFDerivAt.1 (hf _).hasGradientAt).comp_hasDerivAt θ hγ
  exact h.congr_deriv (by rw [toDual_apply_apply, sub_eq_add_neg, neg_smul])

/-- **Interpolation along a quarter circle**:
`f u - f v = ∫₀^{π/2} ⟪∇f (sin θ • u + cos θ • v), cos θ • u - sin θ • v⟫ dθ`. -/
lemma sub_eq_integral_quarterCircle (hf : ContDiff ℝ 1 f) (u v : E) :
    f u - f v =
      ∫ θ in (0 : ℝ)..(π / 2), ⟪gradient f (sin θ • u + cos θ • v), cos θ • u - sin θ • v⟫ := by
  have hcont : Continuous fun θ : ℝ ↦
      ⟪gradient f (sin θ • u + cos θ • v), cos θ • u - sin θ • v⟫ :=
    (hf.continuous_gradient.comp (by fun_prop)).inner (by fun_prop)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt
    (fun θ _ ↦ hasDerivAt_comp_quarterCircle (hf.differentiable one_ne_zero) u v θ)
    (hcont.intervalIntegrable _ _)]
  simp

/-- **Jensen's inequality for the uniform measure on `[0, π/2]`**: for a continuous `h`,
`exp (s ∫₀^{π/2} h) ≤ (2/π) ∫₀^{π/2} exp (π s / 2 · h θ) dθ`. -/
lemma exp_mul_integral_le_integral_exp {h : ℝ → ℝ} (hc : Continuous h) (s : ℝ) :
    exp (s * ∫ θ in (0 : ℝ)..(π / 2), h θ) ≤
      2 / π * ∫ θ in (0 : ℝ)..(π / 2), exp (π * s / 2 * h θ) := by
  have hπ : 0 < π / 2 := by positivity
  set ν : Measure ℝ := volume.restrict (Ioc 0 (π / 2)) with hν
  have hν_univ : ν.real univ = π / 2 := by
    rw [hν, measureReal_restrict_apply_univ, measureReal_def, Real.volume_Ioc, sub_zero,
      ENNReal.toReal_ofReal hπ.le]
  have : NeZero ν := ⟨fun h0 ↦ by simp [h0] at hν_univ; linarith⟩
  have hint : ∀ g : ℝ → ℝ, Continuous g → Integrable g ν := fun g hg ↦
    hg.integrableOn_Icc.mono_set Ioc_subset_Icc_self
  have hJ := ConvexOn.map_average_le (μ := ν) (f := fun θ ↦ π * s / 2 * h θ) convexOn_exp
    continuous_exp.continuousOn isClosed_univ (Filter.Eventually.of_forall fun _ ↦ mem_univ _)
    (hint _ (by fun_prop)) (hint _ (by fun_prop))
  rw [average_eq, average_eq, hν_univ, smul_eq_mul, smul_eq_mul, integral_const_mul] at hJ
  rw [intervalIntegral.integral_of_le hπ.le, intervalIntegral.integral_of_le hπ.le]
  have h1 : (π / 2)⁻¹ * (π * s / 2 * ∫ θ, h θ ∂ν) = s * ∫ θ, h θ ∂ν := by
    field_simp
  rw [h1, inv_div] at hJ
  exact hJ
