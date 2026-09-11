/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Probability.Distributions.Gaussian.Fernique
public import Mathlib.Probability.Distributions.Gaussian.Multivariate
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.EuclideanMatrix

set_option autoImplicit false

/-!
# Transformations of multivariate Gaussian measures

We prove that the class of multivariate Gaussian measures `multivariateGaussian μ S` on
`EuclideanSpace ℝ ι` (with `S` positive semidefinite) is stable under convolution, negation,
scalar multiplication and linear maps given by matrices, and we compute the resulting means
and covariance matrices.

## Main statements

* `multivariateGaussian_conv`: the convolution of `N(μ₁, S₁)` and `N(μ₂, S₂)` is
  `N(μ₁ + μ₂, S₁ + S₂)`.
* `multivariateGaussian_map_toEuclideanCLM`: the image of `N(μ, S)` by the linear map with
  matrix `M` is `N(M μ, M S Mᵀ)`.
* `multivariateGaussian_map_smul`: the image of `N(μ, S)` by `x ↦ c • x` is `N(c • μ, c ^ 2 • S)`.
* `multivariateGaussian_map_neg`: the image of `N(μ, S)` by `x ↦ -x` is `N(-μ, S)`.

All proofs go through characteristic functions.
-/

@[expose] public section

open MeasureTheory Matrix WithLp Complex
open scoped RealInnerProductSpace

namespace ProbabilityTheory

section charFun

variable {E F : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]
  [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]
  {mE : MeasurableSpace E} [OpensMeasurableSpace E] {mF : MeasurableSpace F} [BorelSpace F]

/-- The characteristic function of the image of a measure by a continuous linear map `L` between
Hilbert spaces is the characteristic function of the measure evaluated at `L† t`. -/
lemma charFun_map_clm (μ : Measure E) (L : E →L[ℝ] F) (t : F) :
    charFun (μ.map L) t = charFun μ (ContinuousLinearMap.adjoint L t) := by
  rw [charFun_apply, charFun_apply, integral_map (by fun_prop) (by fun_prop)]
  simp_rw [ContinuousLinearMap.adjoint_inner_right]

end charFun

section multivariateGaussian

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {μ μ₁ μ₂ : EuclideanSpace ℝ ι}
  {S S₁ S₂ : Matrix ι ι ℝ}

/-- The convolution of two multivariate Gaussian measures is a multivariate Gaussian measure:
`N(μ₁, S₁) ∗ N(μ₂, S₂) = N(μ₁ + μ₂, S₁ + S₂)`. -/
lemma multivariateGaussian_conv (hS₁ : S₁.PosSemidef) (hS₂ : S₂.PosSemidef) :
    multivariateGaussian μ₁ S₁ ∗ multivariateGaussian μ₂ S₂ =
      multivariateGaussian (μ₁ + μ₂) (S₁ + S₂) := by
  apply Measure.ext_of_charFun
  ext t
  rw [charFun_conv, charFun_multivariateGaussian hS₁, charFun_multivariateGaussian hS₂,
    charFun_multivariateGaussian (hS₁.add hS₂), ← Complex.exp_add]
  congr 1
  simp only [inner_add_right, add_mulVec, dotProduct_add, ofReal_add]
  ring

/-- The image of the multivariate Gaussian measure `N(μ, S)` by the linear map with matrix `M`
is `N(M μ, M S Mᵀ)`. -/
lemma multivariateGaussian_map_toEuclideanCLM (hS : S.PosSemidef) (M : Matrix ι ι ℝ) :
    (multivariateGaussian μ S).map (toEuclideanCLM (𝕜 := ℝ) M) =
      multivariateGaussian (toEuclideanCLM (𝕜 := ℝ) M μ) (M * S * Mᵀ) := by
  apply Measure.ext_of_charFun
  ext t
  have hMSM : (M * S * Mᵀ).PosSemidef := by
    simpa [conjTranspose_eq_transpose_of_trivial] using hS.mul_mul_conjTranspose_same M
  rw [charFun_map_clm, adjoint_toEuclideanCLM, charFun_multivariateGaussian hS,
    charFun_multivariateGaussian hMSM]
  congr 3
  · rw [real_inner_comm, inner_toEuclideanCLM, inner_toEuclideanCLM,
      dotProduct_mulVec (ofLp μ) Mᵀ, vecMul_transpose, dotProduct_comm]
  · rw [ofLp_toEuclideanCLM, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec (ofLp t) M,
      ← mulVec_transpose]

/-- The image of the centered multivariate Gaussian measure `N(0, S)` by the linear map with
matrix `M` is `N(0, M S Mᵀ)`. -/
lemma multivariateGaussian_zero_map_toEuclideanCLM (hS : S.PosSemidef) (M : Matrix ι ι ℝ) :
    (multivariateGaussian 0 S).map (toEuclideanCLM (𝕜 := ℝ) M) =
      multivariateGaussian 0 (M * S * Mᵀ) := by
  simpa using multivariateGaussian_map_toEuclideanCLM (μ := 0) hS M

/-- The image of the multivariate Gaussian measure `N(μ, S)` on `EuclideanSpace ℝ ι` by the
linear map `EuclideanSpace ℝ ι → EuclideanSpace ℝ κ` with (rectangular) matrix `M` is
`N(M μ, M S Mᵀ)`. -/
lemma multivariateGaussian_map_toEuclideanLin {κ : Type*} [Fintype κ] [DecidableEq κ]
    (hS : S.PosSemidef) (M : Matrix κ ι ℝ) :
    (multivariateGaussian μ S).map (LinearMap.toContinuousLinearMap (toEuclideanLin M)) =
      multivariateGaussian (toEuclideanLin M μ) (M * S * Mᵀ) := by
  apply Measure.ext_of_charFun
  ext t
  have hMSM : (M * S * Mᵀ).PosSemidef := by
    simpa [conjTranspose_eq_transpose_of_trivial] using hS.mul_mul_conjTranspose_same M
  rw [charFun_map_clm, ← LinearMap.adjoint_toContinuousLinearMap,
    ← toEuclideanLin_conjTranspose_eq_adjoint, conjTranspose_eq_transpose_of_trivial,
    charFun_multivariateGaussian hS, charFun_multivariateGaussian hMSM]
  congr 3
  · simp only [LinearMap.coe_toContinuousLinearMap', EuclideanSpace.inner_eq_star_dotProduct,
      star_trivial, ofLp_toLpLin, toLin'_apply]
    rw [dotProduct_mulVec, vecMul_transpose]
  · simp only [LinearMap.coe_toContinuousLinearMap', ofLp_toLpLin, toLin'_apply]
    rw [← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec (ofLp t) M, ← mulVec_transpose]

/-- The image of the centered multivariate Gaussian measure `N(0, S)` on `EuclideanSpace ℝ ι` by
the linear map with (rectangular) matrix `M` is `N(0, M S Mᵀ)`. -/
lemma multivariateGaussian_zero_map_toEuclideanLin {κ : Type*} [Fintype κ] [DecidableEq κ]
    (hS : S.PosSemidef) (M : Matrix κ ι ℝ) :
    (multivariateGaussian 0 S).map (LinearMap.toContinuousLinearMap (toEuclideanLin M)) =
      multivariateGaussian 0 (M * S * Mᵀ) := by
  simpa using multivariateGaussian_map_toEuclideanLin (μ := 0) hS M

/-- The image of the multivariate Gaussian measure `N(μ, S)` by `x ↦ c • x` is
`N(c • μ, c ^ 2 • S)`. -/
lemma multivariateGaussian_map_smul (hS : S.PosSemidef) (c : ℝ) :
    (multivariateGaussian μ S).map (fun x ↦ c • x) =
      multivariateGaussian (c • μ) (c ^ 2 • S) := by
  apply Measure.ext_of_charFun
  ext t
  rw [charFun_map_smul, charFun_multivariateGaussian hS,
    charFun_multivariateGaussian (hS.smul (sq_nonneg c))]
  congr 3
  · rw [real_inner_smul_left, real_inner_smul_right]
  · congr 1
    simp only [ofLp_smul, smul_mulVec, mulVec_smul, dotProduct_smul, smul_dotProduct,
      smul_eq_mul]
    ring

/-- The image of the centered multivariate Gaussian measure `N(0, S)` by `x ↦ c • x` is
`N(0, c ^ 2 • S)`. -/
lemma multivariateGaussian_zero_map_smul (hS : S.PosSemidef) (c : ℝ) :
    (multivariateGaussian 0 S).map (fun x ↦ c • x) = multivariateGaussian 0 (c ^ 2 • S) := by
  simpa using multivariateGaussian_map_smul (μ := 0) hS c

/-- The image of the multivariate Gaussian measure `N(μ, S)` by `x ↦ -x` is `N(-μ, S)`. -/
lemma multivariateGaussian_map_neg (hS : S.PosSemidef) :
    (multivariateGaussian μ S).map (fun x ↦ -x) = multivariateGaussian (-μ) S := by
  have : (fun x : EuclideanSpace ℝ ι ↦ -x) = fun x ↦ (-1 : ℝ) • x := by simp
  rw [this, multivariateGaussian_map_smul hS]
  simp

/-- A centered multivariate Gaussian measure is symmetric: the image of `N(0, S)` by `x ↦ -x`
is `N(0, S)`. -/
lemma multivariateGaussian_zero_map_neg (hS : S.PosSemidef) :
    (multivariateGaussian 0 S).map (fun x ↦ -x) = multivariateGaussian 0 S := by
  simpa using multivariateGaussian_map_neg (μ := 0) hS

end multivariateGaussian

section translation

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Translating a centered multivariate Gaussian by `μ` gives `N(μ, S)`. -/
lemma multivariateGaussian_zero_map_const_add (μ : EuclideanSpace ℝ ι) (S : Matrix ι ι ℝ) :
    (multivariateGaussian 0 S).map (fun v ↦ μ + v) = multivariateGaussian μ S := by
  rw [multivariateGaussian, multivariateGaussian, Measure.map_map (measurable_const_add μ)
    (by fun_prop)]
  simp [Function.comp_def]

/-- `N(μ, S)` shifted by `-μ` is `N(0, S)`. -/
lemma multivariateGaussian_map_sub_const (μ : EuclideanSpace ℝ ι) (S : Matrix ι ι ℝ) :
    (multivariateGaussian μ S).map (fun v ↦ v - μ) = multivariateGaussian 0 S := by
  rw [← multivariateGaussian_zero_map_const_add μ S, Measure.map_map (by fun_prop) (by fun_prop)]
  simp [Function.comp_def]

end translation

section Moments

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {S : Matrix ι ι ℝ}

lemma memLp_two_eval_multivariateGaussian (μ : EuclideanSpace ℝ ι) (S : Matrix ι ι ℝ) (i : ι) :
    MemLp (fun x : EuclideanSpace ℝ ι ↦ x i) 2 (multivariateGaussian μ S) := by
  have h := IsGaussian.memLp_two_id
    (μ := (multivariateGaussian μ S).map (EuclideanSpace.proj (𝕜 := ℝ) i))
  rw [memLp_map_measure_iff aestronglyMeasurable_id (by fun_prop)] at h
  simpa using h

lemma integral_eval_multivariateGaussian (μ : EuclideanSpace ℝ ι) (S : Matrix ι ι ℝ) (i : ι) :
    ∫ x, x i ∂(multivariateGaussian μ S) = μ i := by
  have h := (EuclideanSpace.proj (𝕜 := ℝ) i).integral_comp_comm
    (IsGaussian.integrable_id (μ := multivariateGaussian μ S))
  simpa using h

/-- `E‖G‖² = tr S` for a centered Gaussian vector `G ~ N(0, S)`. -/
lemma integral_norm_sq_multivariateGaussian (hS : S.PosSemidef) :
    ∫ x, ‖x‖ ^ 2 ∂(multivariateGaussian 0 S) = S.trace := by
  simp_rw [EuclideanSpace.real_norm_sq_eq]
  rw [integral_finsetSum _ fun i _ ↦ (memLp_two_eval_multivariateGaussian 0 S i).integrable_sq]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Matrix.diag_apply, ← variance_eval_multivariateGaussian hS i,
    variance_of_integral_eq_zero (by fun_prop)]
  simpa using integral_eval_multivariateGaussian 0 S i

/-- `E‖G‖ ≤ √(tr S)` for a centered Gaussian vector `G ~ N(0, S)`. -/
lemma integral_norm_le_sqrt_trace_multivariateGaussian (hS : S.PosSemidef) :
    ∫ x, ‖x‖ ∂(multivariateGaussian 0 S) ≤ √S.trace := by
  rw [← integral_norm_sq_multivariateGaussian hS,
    Real.le_sqrt (integral_nonneg fun _ ↦ norm_nonneg _) (integral_nonneg fun _ ↦ sq_nonneg _)]
  have hL2 : MemLp (fun x : EuclideanSpace ℝ ι ↦ ‖x‖) 2 (multivariateGaussian 0 S) :=
    (IsGaussian.memLp_two_id (μ := multivariateGaussian 0 S)).norm
  have h := variance_eq_sub hL2
  have h0 := variance_nonneg (fun x : EuclideanSpace ℝ ι ↦ ‖x‖) (multivariateGaussian 0 S)
  simp only [Pi.pow_apply] at h
  linarith

omit [DecidableEq ι] in
/-- `E‖g‖ ≤ √d` for a standard Gaussian vector `g` in dimension `d`. -/
lemma integral_norm_stdGaussian_le :
    ∫ x, ‖x‖ ∂(stdGaussian (EuclideanSpace ℝ ι)) ≤ √(Fintype.card ι) := by
  classical
  have h := integral_norm_le_sqrt_trace_multivariateGaussian (ι := ι) Matrix.PosSemidef.one
  rwa [multivariateGaussian_zero_one, Matrix.trace_one] at h

end Moments

end ProbabilityTheory
