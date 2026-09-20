/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Probability.Distributions.Gaussian.Real
public import Mathlib.MeasureTheory.Integral.IntegralEqImproper

set_option autoImplicit false

/-!
# Stein's identity for real Gaussian measures

For a differentiable function `f : ℝ → ℝ` with bounded derivative and `ξ ~ N(μ, v)`, Gaussian
integration by parts gives `E[(ξ - μ) f(ξ)] = v E[f'(ξ)]`.

## Main results

* `hasDerivAt_gaussianPDFReal`: the derivative of the Gaussian density `φ` is
  `φ' x = -((x - μ) / v) * φ x`.
* `integrable_sub_mul_gaussianReal`, `integrable_deriv_gaussianReal`: integrability of
  `x ↦ (x - μ) * f x` and of `deriv f` under `gaussianReal μ v` when `deriv f` is bounded.
* `integral_sub_mul_gaussianReal`: **Stein's identity** `E[(ξ - μ) f(ξ)] = v E[f'(ξ)]` for
  `ξ ~ N(μ, v)`.
* `integral_mul_gaussianReal_zero_one`: the standard case `E[ξ f(ξ)] = E[f'(ξ)]` for
  `ξ ~ N(0, 1)`, and its `HasDerivAt` form `integral_mul_gaussianReal_zero_one_of_hasDerivAt`.
-/

@[expose] public section

open MeasureTheory Real
open scoped NNReal

namespace ProbabilityTheory

/-- The derivative of the Gaussian density `φ = gaussianPDFReal μ v` is
`φ' x = -((x - μ) / v) * φ x`. Also true (as `0 = 0`) for `v = 0`. -/
lemma hasDerivAt_gaussianPDFReal (μ : ℝ) (v : ℝ≥0) (x : ℝ) :
    HasDerivAt (gaussianPDFReal μ v) (-((x - μ) / v) * gaussianPDFReal μ v x) x := by
  simp only [gaussianPDFReal]
  have h : HasDerivAt (fun y ↦ -(y - μ) ^ 2 / (2 * v)) (-(x - μ) / v) x := by
    refine ((((hasDerivAt_id' x).sub_const μ).fun_pow 2).fun_neg.div_const
      (2 * (v : ℝ))).congr_deriv ?_
    norm_num
    ring
  refine (h.exp.const_mul (√(2 * π * v))⁻¹).congr_deriv ?_
  ring

/-- For `v ≠ 0`, a function is integrable under `gaussianReal μ v` iff its product with the
Gaussian density is integrable under the Lebesgue measure. -/
lemma integrable_gaussianReal_iff {μ : ℝ} {v : ℝ≥0} (hv : v ≠ 0) {g : ℝ → ℝ} :
    Integrable g (gaussianReal μ v) ↔ Integrable (fun x ↦ gaussianPDFReal μ v x * g x) := by
  rw [gaussianReal_of_var_ne_zero _ hv, integrable_withDensity_iff_integrable_smul'
    (measurable_gaussianPDF _ _) (ae_of_all _ fun _ ↦ gaussianPDF_lt_top)]
  simp only [toReal_gaussianPDF, smul_eq_mul]

variable {f : ℝ → ℝ} {L : ℝ}

/-- A differentiable function whose derivative is bounded by `L` has linear growth:
`|f x| ≤ |f y| + L * |x - y|`. -/
lemma abs_le_add_mul_abs_sub_of_abs_deriv_le (hf : Differentiable ℝ f)
    (hL : ∀ x, |deriv f x| ≤ L) (x y : ℝ) :
    |f x| ≤ |f y| + L * |x - y| := by
  have h := convex_univ.norm_image_sub_le_of_norm_deriv_le (fun z _ ↦ hf z)
    (fun z _ ↦ (Real.norm_eq_abs _).trans_le (hL z)) (Set.mem_univ y) (Set.mem_univ x)
  simp only [Real.norm_eq_abs] at h
  calc |f x| = |f y + (f x - f y)| := by congr 1; ring
    _ ≤ |f y| + |f x - f y| := abs_add_le _ _
    _ ≤ |f y| + L * |x - y| := by gcongr

/-- A differentiable function with bounded derivative is integrable under `gaussianReal μ v`. -/
lemma integrable_gaussianReal_of_abs_deriv_le (hf : Differentiable ℝ f)
    (hL : ∀ x, |deriv f x| ≤ L) (μ : ℝ) (v : ℝ≥0) :
    Integrable f (gaussianReal μ v) := by
  have h_abs : Integrable (fun x ↦ |x - μ|) (gaussianReal μ v) :=
    (((memLp_id_gaussianReal' 1 (by simp)).integrable le_rfl).sub (integrable_const μ)).abs
  refine Integrable.mono' (g := fun x ↦ |f μ| + L * |x - μ|)
    ((integrable_const _).add (h_abs.const_mul L)) hf.continuous.aestronglyMeasurable
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_eq_abs]
  exact abs_le_add_mul_abs_sub_of_abs_deriv_le hf hL x μ

/-- `x ↦ (x - μ) * f x` is integrable under `gaussianReal μ v` when `f` has bounded derivative
(hence linear growth). -/
lemma integrable_sub_mul_gaussianReal (hf : Differentiable ℝ f) (hL : ∀ x, |deriv f x| ≤ L)
    (μ : ℝ) (v : ℝ≥0) :
    Integrable (fun x ↦ (x - μ) * f x) (gaussianReal μ v) := by
  have h_sub : MemLp (fun x ↦ x - μ) 2 (gaussianReal μ v) :=
    (memLp_id_gaussianReal' 2 (by simp)).sub (memLp_const μ)
  have h_abs : Integrable (fun x ↦ |x - μ|) (gaussianReal μ v) :=
    (h_sub.integrable (by norm_num)).abs
  have hfc := hf.continuous
  refine Integrable.mono' (g := fun x ↦ |f μ| * |x - μ| + L * (x - μ) ^ 2)
    ((h_abs.const_mul _).add (h_sub.integrable_sq.const_mul L)) (by fun_prop)
    (Filter.Eventually.of_forall fun x ↦ ?_)
  rw [Real.norm_eq_abs, abs_mul, ← sq_abs]
  calc |x - μ| * |f x| ≤ |x - μ| * (|f μ| + L * |x - μ|) := by
        gcongr
        exact abs_le_add_mul_abs_sub_of_abs_deriv_le hf hL x μ
    _ = |f μ| * |x - μ| + L * |x - μ| ^ 2 := by ring

/-- A bounded derivative is integrable under `gaussianReal μ v`. -/
lemma integrable_deriv_gaussianReal (hL : ∀ x, |deriv f x| ≤ L) (μ : ℝ) (v : ℝ≥0) :
    Integrable (deriv f) (gaussianReal μ v) :=
  Integrable.of_bound (measurable_deriv f).aestronglyMeasurable L
    (Filter.Eventually.of_forall fun x ↦ (Real.norm_eq_abs _).trans_le (hL x))

/-- **Stein's identity** for `N(μ, v)`: `E[(ξ - μ) f(ξ)] = v E[f'(ξ)]` when `f` is differentiable
with bounded derivative. Also true for `v = 0`. -/
lemma integral_sub_mul_gaussianReal (hf : Differentiable ℝ f) (hL : ∀ x, |deriv f x| ≤ L)
    (μ : ℝ) (v : ℝ≥0) :
    ∫ x, (x - μ) * f x ∂gaussianReal μ v = v * ∫ x, deriv f x ∂gaussianReal μ v := by
  by_cases hv : v = 0
  · simp [hv]
  have hv' : (v : ℝ) ≠ 0 := NNReal.coe_ne_zero.mpr hv
  -- integration by parts on `ℝ` with `u = gaussianPDFReal μ v` and `v = f`
  have h_ibp : ∫ x, gaussianPDFReal μ v x * deriv f x
      = - ∫ x, (-((x - μ) / v) * gaussianPDFReal μ v x) * f x := by
    refine integral_mul_deriv_eq_deriv_mul_of_integrable
      (fun x _ ↦ hasDerivAt_gaussianPDFReal μ v x) (fun x _ ↦ (hf x).hasDerivAt) ?_ ?_ ?_
    · exact (integrable_gaussianReal_iff hv).mp (integrable_deriv_gaussianReal hL μ v)
    · have h := ((integrable_gaussianReal_iff hv).mp
        (integrable_sub_mul_gaussianReal hf hL μ v)).const_mul (-(v : ℝ)⁻¹)
      refine h.congr (Filter.Eventually.of_forall fun x ↦ ?_)
      simp only [Pi.mul_apply]
      ring
    · exact (integrable_gaussianReal_iff hv).mp (integrable_gaussianReal_of_abs_deriv_le hf hL μ v)
  calc ∫ x, (x - μ) * f x ∂gaussianReal μ v
      = ∫ x, gaussianPDFReal μ v x * ((x - μ) * f x) := by
        simp only [integral_gaussianReal_eq_integral_smul hv, smul_eq_mul]
    _ = ∫ x, -v * ((-((x - μ) / v) * gaussianPDFReal μ v x) * f x) := by
        congr 1 with x
        field_simp
    _ = -v * ∫ x, (-((x - μ) / v) * gaussianPDFReal μ v x) * f x := integral_const_mul _ _
    _ = v * ∫ x, gaussianPDFReal μ v x * deriv f x := by rw [h_ibp]; ring
    _ = v * ∫ x, deriv f x ∂gaussianReal μ v := by
        simp only [integral_gaussianReal_eq_integral_smul hv, smul_eq_mul]

/-- **Stein's identity** for the standard Gaussian: `E[ξ f(ξ)] = E[f'(ξ)]` when `f` is
differentiable with bounded derivative. -/
lemma integral_mul_gaussianReal_zero_one (hf : Differentiable ℝ f) (hL : ∀ x, |deriv f x| ≤ L) :
    ∫ x, x * f x ∂gaussianReal 0 1 = ∫ x, deriv f x ∂gaussianReal 0 1 := by
  simpa using integral_sub_mul_gaussianReal hf hL 0 1

/-- **Stein's identity** for the standard Gaussian, `HasDerivAt` form: if `f'` is a bounded
derivative of `f`, then `E[ξ f(ξ)] = E[f'(ξ)]` for `ξ ~ N(0, 1)`. -/
lemma integral_mul_gaussianReal_zero_one_of_hasDerivAt {f' : ℝ → ℝ}
    (hf : ∀ x, HasDerivAt f (f' x) x) (hL : ∀ x, |f' x| ≤ L) :
    ∫ x, x * f x ∂gaussianReal 0 1 = ∫ x, f' x ∂gaussianReal 0 1 := by
  have hf' : deriv f = f' := funext fun x ↦ (hf x).deriv
  rw [← hf'] at hL ⊢
  exact integral_mul_gaussianReal_zero_one (fun x ↦ (hf x).differentiableAt) hL

end ProbabilityTheory
