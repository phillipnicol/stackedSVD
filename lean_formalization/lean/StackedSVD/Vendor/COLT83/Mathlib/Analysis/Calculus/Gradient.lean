/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.Calculus.Gradient.Basic
public import Mathlib.Analysis.Calculus.MeanValue
public import Mathlib.Analysis.Calculus.ContDiff.Basic
public import Mathlib.Analysis.Calculus.ContDiff.Comp

set_option autoImplicit false

/-!
# Functions with a bounded gradient

* `norm_le_of_norm_fderiv_le`, `abs_le_of_norm_fderiv_le`: a differentiable function with
  `‖Df‖ ≤ L` has linear growth `‖f x‖ ≤ ‖f 0‖ + L ‖x‖`;
* `norm_fderiv_le_of_norm_fderiv_fderiv_le`, `norm_le_of_norm_fderiv_fderiv_le`: a `C²` function
  with `‖D²f‖ ≤ K` has a derivative of linear growth and is itself of quadratic growth;
* `lipschitzWith_of_norm_gradient_le`: a differentiable function with `‖∇f‖ ≤ L` is
  `L`-Lipschitz;
* `abs_le_of_norm_gradient_le`: it has linear growth `|f x| ≤ |f 0| + L ‖x‖`;
* `abs_inner_gradient_le`: `|⟪∇f x, v⟫| ≤ L ‖v‖`;
* `ContDiff.continuous_gradient`: the gradient of a `C¹` function is continuous.

## TODO

The `fderiv` statements of the first section are proved for real normed spaces only (through
`Convex.norm_image_sub_le_of_norm_fderiv_le`); they hold over any `RCLike` field.
-/

@[expose] public section

open InnerProductSpace
open scoped RealInnerProductSpace NNReal

section fderiv

variable {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [NormedAddCommGroup F]
  [NormedSpace ℝ F] {f : E → F} {L K : ℝ}

/-- A differentiable function with `‖Df‖ ≤ L` has linear growth. -/
lemma norm_le_of_norm_fderiv_le (hf : Differentiable ℝ f) (hL : ∀ x, ‖fderiv ℝ f x‖ ≤ L)
    (x : E) : ‖f x‖ ≤ ‖f 0‖ + L * ‖x‖ := by
  have h := convex_univ.norm_image_sub_le_of_norm_fderiv_le (fun y _ ↦ hf y) (fun y _ ↦ hL y)
    (Set.mem_univ 0) (Set.mem_univ x)
  rw [sub_zero] at h
  calc ‖f x‖ = ‖f x - f 0 + f 0‖ := by rw [sub_add_cancel]
    _ ≤ ‖f x - f 0‖ + ‖f 0‖ := norm_add_le _ _
    _ ≤ ‖f 0‖ + L * ‖x‖ := by linarith

/-- A differentiable real function with `‖Df‖ ≤ L` has linear growth. -/
lemma abs_le_of_norm_fderiv_le {f : E → ℝ} (hf : Differentiable ℝ f)
    (hL : ∀ x, ‖fderiv ℝ f x‖ ≤ L) (x : E) : |f x| ≤ |f 0| + L * ‖x‖ :=
  norm_le_of_norm_fderiv_le hf hL x

/-- The derivative of a `C²` function with bounded second derivative has linear growth. -/
lemma norm_fderiv_le_of_norm_fderiv_fderiv_le (hf : ContDiff ℝ 2 f)
    (hK : ∀ x, ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ K) (x : E) :
    ‖fderiv ℝ f x‖ ≤ ‖fderiv ℝ f 0‖ + K * ‖x‖ :=
  norm_le_of_norm_fderiv_le ((hf.fderiv_right (m := 1) le_rfl).differentiable one_ne_zero) hK x

/-- A `C²` function with bounded second derivative has quadratic growth. -/
lemma norm_le_of_norm_fderiv_fderiv_le (hf : ContDiff ℝ 2 f)
    (hK : ∀ x, ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ K) (x : E) :
    ‖f x‖ ≤ ‖f 0‖ + (‖fderiv ℝ f 0‖ + K * ‖x‖) * ‖x‖ := by
  have hK0 : 0 ≤ K := le_trans (norm_nonneg (fderiv ℝ (fderiv ℝ f) 0)) (hK 0)
  have h := (convex_closedBall (0 : E) ‖x‖).norm_image_sub_le_of_norm_fderiv_le
    (fun y _ ↦ hf.differentiable (by norm_num) y)
    (fun y hy ↦ (norm_fderiv_le_of_norm_fderiv_fderiv_le hf hK y).trans
      (add_le_add_right (mul_le_mul_of_nonneg_left (mem_closedBall_zero_iff.1 hy) hK0) _))
    (Metric.mem_closedBall_self (norm_nonneg x)) (mem_closedBall_zero_iff.2 le_rfl)
  rw [sub_zero] at h
  calc ‖f x‖ = ‖f x - f 0 + f 0‖ := by rw [sub_add_cancel]
    _ ≤ ‖f x - f 0‖ + ‖f 0‖ := norm_add_le _ _
    _ ≤ ‖f 0‖ + (‖fderiv ℝ f 0‖ + K * ‖x‖) * ‖x‖ := by linarith

end fderiv

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]
  {f : E → ℝ} {L : ℝ≥0}

/-- A differentiable function with `‖∇f‖ ≤ L` is `L`-Lipschitz. -/
lemma lipschitzWith_of_norm_gradient_le (hf : Differentiable ℝ f)
    (hL : ∀ x, ‖gradient f x‖ ≤ L) : LipschitzWith L f := by
  refine lipschitzWith_of_nnnorm_fderiv_le hf fun x ↦ ?_
  rw [← toDual_gradient, LinearIsometryEquiv.nnnorm_map]
  exact_mod_cast hL x

/-- A differentiable function with `‖∇f‖ ≤ L` has linear growth. -/
lemma abs_le_of_norm_gradient_le (hf : Differentiable ℝ f) (hL : ∀ x, ‖gradient f x‖ ≤ L)
    (x : E) : |f x| ≤ |f 0| + L * ‖x‖ := by
  have h := (lipschitzWith_of_norm_gradient_le hf hL).dist_le_mul x 0
  rw [Real.dist_eq, dist_zero_right] at h
  calc |f x| = |f x - f 0 + f 0| := by ring_nf
    _ ≤ |f x - f 0| + |f 0| := abs_add_le _ _
    _ ≤ L * ‖x‖ + |f 0| := by linarith
    _ = |f 0| + L * ‖x‖ := by ring

lemma abs_inner_gradient_le (hL : ∀ x, ‖gradient f x‖ ≤ L) (x v : E) :
    |⟪gradient f x, v⟫| ≤ L * ‖v‖ :=
  (abs_real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hL x) (norm_nonneg _))

/-- The gradient of a `C¹` function is continuous. -/
lemma ContDiff.continuous_gradient (hf : ContDiff ℝ 1 f) : Continuous (gradient f) := by
  have h : gradient f = (toDual ℝ E).symm ∘ fderiv ℝ f := by
    funext x
    rw [Function.comp_apply, ← toDual_gradient, LinearIsometryEquiv.symm_apply_apply]
  rw [h]
  exact (toDual ℝ E).symm.continuous.comp (hf.continuous_fderiv one_ne_zero)
