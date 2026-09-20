/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.InnerProductSpace.Adjoint
public import Mathlib.MeasureTheory.Group.Convolution
public import Mathlib.MeasureTheory.Group.IntegralConvolution
public import Mathlib.Probability.Distributions.Gaussian.Fernique
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.SupportFn
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.MultivariateGaussian
public import Mathlib.Probability.Distributions.Gaussian.Multivariate

set_option autoImplicit false

/-!
# Gaussian widths

For a measure `μ` on a real inner product space `E` (typically a centered Gaussian measure), the
*Gaussian width* of a set `K ⊆ E` for `μ` is `gaussianWidth K μ = ∫ ξ, supportFn K ξ ∂μ`, the
expectation of the support function `supportFn K ξ = sup_{x ∈ K} ⟪x, ξ⟫` (see
`COLT83.Mathlib.Analysis.InnerProductSpace.SupportFn`). For the standard Gaussian measure this
is the usual Gaussian width (or mean width) of `K`.

## Main results

* `integrable_supportFn`: the support function of a nonempty bounded set is integrable with
  respect to every measure with an integrable norm (Gaussian measures in particular).
* `gaussianWidth_mono`, `gaussianWidth_smul_set`, `gaussianWidth_map`, `gaussianWidth_add`:
  monotonicity, homogeneity, images by linear maps, Minkowski sums.
* `gaussianWidth_le_gaussianWidth_conv`: the Gaussian width of a compact set for `μ` is at most
  its width for `μ ∗ ν` when `ν` is centered (Jensen). For Gaussian measures this is the
  comparison by covariance domination: `N(0, S₁)` has a smaller width than `N(0, S₂)` when
  `S₁ ≤ S₂` (`gaussianWidth_multivariateGaussian_le`).
-/

@[expose] public section

open MeasureTheory Set

open scoped RealInnerProductSpace Pointwise NNReal

namespace ProbabilityTheory

section GaussianWidth

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] {mE : MeasurableSpace E}
  [OpensMeasurableSpace E] {μ ν : Measure E} {K K' : Set E} {R : ℝ}

lemma measurable_ciSup_inner {ι : Type*} [Finite ι] (y : ι → E) :
    Measurable fun v : E ↦ ⨆ i, ⟪y i, v⟫ :=
  Measurable.iSup fun _ ↦ (continuous_const.inner continuous_id).measurable

lemma measurable_ciSup_apply {κ : Type*} [Finite κ] :
    Measurable fun x : EuclideanSpace ℝ κ ↦ ⨆ k, x k :=
  Measurable.iSup fun k ↦ (PiLp.proj (𝕜 := ℝ) 2 (fun _ : κ ↦ ℝ) k).continuous.measurable

/-- Gaussian width of the set `K` with respect to the measure `μ`:
`∫ ξ, sup_{x ∈ K} ⟪x, ξ⟫ ∂μ`. -/
noncomputable def gaussianWidth (K : Set E) (μ : Measure E) : ℝ := ∫ ξ, supportFn K ξ ∂μ

/-- The support function of a nonempty bounded set is integrable for every measure with an
integrable norm. -/
lemma integrable_supportFn (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R)
    (hμ : Integrable (fun ξ ↦ ‖ξ‖) μ) : Integrable (supportFn K) μ :=
  (hμ.const_mul R).mono' (continuous_supportFn hne hK).aestronglyMeasurable
    (ae_of_all _ fun ξ ↦ by simpa [Real.norm_eq_abs] using abs_supportFn_le hne hK ξ)

lemma integrable_supportFn_of_isGaussian [BorelSpace E] [CompleteSpace E]
    [SecondCountableTopology E] [IsGaussian μ] (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) :
    Integrable (supportFn K) μ :=
  integrable_supportFn hne hK (IsGaussian.integrable_id (μ := μ)).norm

lemma abs_gaussianWidth_le (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R)
    (hμ : Integrable (fun ξ ↦ ‖ξ‖) μ) :
    |gaussianWidth K μ| ≤ R * ∫ ξ, ‖ξ‖ ∂μ := by
  rw [gaussianWidth, ← integral_const_mul]
  exact abs_integral_le_integral_abs.trans (integral_mono (integrable_supportFn hne hK hμ).abs
    (hμ.const_mul R) fun ξ ↦ abs_supportFn_le hne hK ξ)

omit [OpensMeasurableSpace E] in
lemma gaussianWidth_of_subsingleton [Subsingleton E] (K : Set E) (μ : Measure E) :
    gaussianWidth K μ = 0 := by
  simp [gaussianWidth, supportFn_of_subsingleton]

lemma gaussianWidth_mono (hne : K.Nonempty) (hK' : ∀ x ∈ K', ‖x‖ ≤ R) (hsub : K ⊆ K')
    (hμ : Integrable (fun ξ ↦ ‖ξ‖) μ) :
    gaussianWidth K μ ≤ gaussianWidth K' μ :=
  integral_mono (integrable_supportFn hne (fun x hx ↦ hK' x (hsub hx)) hμ)
    (integrable_supportFn (hne.mono hsub) hK' hμ) (supportFn_mono hne hK' hsub)

omit [OpensMeasurableSpace E] in
lemma gaussianWidth_smul_set (K : Set E) {c : ℝ} (hc : 0 ≤ c) (μ : Measure E) :
    gaussianWidth (c • K) μ = c * gaussianWidth K μ := by
  simp only [gaussianWidth, supportFn_smul_set K hc]
  exact integral_const_mul c _

lemma gaussianWidth_neg_set [BorelSpace E] (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R)
    (μ : Measure E) :
    gaussianWidth (-K) μ = gaussianWidth K (μ.map fun ξ ↦ -ξ) := by
  rw [gaussianWidth, gaussianWidth, integral_map continuous_neg.measurable.aemeasurable
    (continuous_supportFn hne hK).aestronglyMeasurable]
  simp only [supportFn_neg_set]

lemma gaussianWidth_map_smul [BorelSpace E] (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R)
    {c : ℝ} (hc : 0 ≤ c) (μ : Measure E) :
    gaussianWidth K (μ.map fun ξ ↦ c • ξ) = c * gaussianWidth K μ := by
  rw [gaussianWidth, gaussianWidth, integral_map (continuous_const_smul c).measurable.aemeasurable
    (continuous_supportFn hne hK).aestronglyMeasurable, ← integral_const_mul]
  simp only [supportFn_smul K hc]

/-- The Gaussian width of `K` for the image measure `μ.map L` is the width of the image of `K`
by the adjoint of `L` for `μ`. -/
lemma gaussianWidth_map [CompleteSpace E] {F : Type*} [NormedAddCommGroup F]
    [InnerProductSpace ℝ F] [CompleteSpace F] {mF : MeasurableSpace F} [BorelSpace F]
    {K : Set F} (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (L : E →L[ℝ] F) (μ : Measure E) :
    gaussianWidth K (μ.map L) = gaussianWidth (ContinuousLinearMap.adjoint L '' K) μ := by
  rw [gaussianWidth, gaussianWidth, integral_map L.continuous.measurable.aemeasurable
    (continuous_supportFn hne hK).aestronglyMeasurable]
  congr 1
  ext ξ
  exact (supportFn_image K (fun x ξ ↦ ContinuousLinearMap.adjoint_inner_left L ξ x) ξ).symm

lemma inner_integral_le_gaussianWidth [CompleteSpace E] (hK : ∀ x ∈ K, ‖x‖ ≤ R) {x : E} (hx : x ∈ K)
    (hμ : Integrable id μ) :
    ⟪x, ∫ ξ, ξ ∂μ⟫ ≤ gaussianWidth K μ := by
  have h := integral_inner (𝕜 := ℝ) hμ x
  simp only [id] at h
  rw [← h]
  exact integral_mono (hμ.const_inner x) (integrable_supportFn ⟨x, hx⟩ hK hμ.norm)
    fun ξ ↦ inner_le_supportFn hK hx ξ

/-- The Gaussian width of a nonempty bounded set for a centered measure is nonnegative. -/
lemma gaussianWidth_nonneg_of_integral_eq_zero [CompleteSpace E] (hne : K.Nonempty)
    (hK : ∀ x ∈ K, ‖x‖ ≤ R) (hμ : Integrable id μ) (h0 : ∫ ξ, ξ ∂μ = 0) :
    0 ≤ gaussianWidth K μ := by
  obtain ⟨x, hx⟩ := hne
  have := inner_integral_le_gaussianWidth hK hx hμ
  rwa [h0, inner_zero_right] at this

omit [OpensMeasurableSpace E] in
lemma gaussianWidth_nonneg_of_zero_mem (hK : ∀ x ∈ K, ‖x‖ ≤ R) (h0 : (0 : E) ∈ K) :
    0 ≤ gaussianWidth K μ :=
  integral_nonneg fun ξ ↦ by simpa using inner_le_supportFn hK h0 ξ

/-- The Gaussian width of a Minkowski sum is the sum of the Gaussian widths. -/
lemma gaussianWidth_add {R' : ℝ} (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (hne' : K'.Nonempty)
    (hK' : ∀ x ∈ K', ‖x‖ ≤ R') (hμ : Integrable (fun ξ ↦ ‖ξ‖) μ) :
    gaussianWidth (K + K') μ = gaussianWidth K μ + gaussianWidth K' μ := by
  simp_rw [gaussianWidth, supportFn_add hne hK hne' hK']
  exact integral_add (integrable_supportFn hne hK hμ) (integrable_supportFn hne' hK' hμ)

omit [InnerProductSpace ℝ E] in
lemma _root_.MeasureTheory.integrable_norm_conv [BorelSpace E] [SecondCountableTopology E]
    [IsFiniteMeasure μ] [IsFiniteMeasure ν] (hμ : Integrable (fun ξ ↦ ‖ξ‖) μ)
    (hν : Integrable (fun ξ ↦ ‖ξ‖) ν) :
    Integrable (fun ξ ↦ ‖ξ‖) (μ ∗ ν) := by
  rw [Measure.conv, integrable_map_measure continuous_norm.aestronglyMeasurable (by fun_prop)]
  refine ((hμ.comp_fst ν).add (hν.comp_snd μ)).mono' (by fun_prop) (ae_of_all _ fun p ↦ ?_)
  simpa using norm_add_le p.1 p.2

/-- **Comparison by convolution.** For a compact set `K` and a centered probability measure `ν`,
the Gaussian width of `K` for `μ` is at most its width for `μ ∗ ν`: adding an independent
centered noise can only increase the expected supremum (Jensen's inequality). -/
lemma gaussianWidth_le_gaussianWidth_conv [BorelSpace E] [SecondCountableTopology E]
    [CompleteSpace E] [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] (hK : IsCompact K)
    (hne : K.Nonempty) (hμ : Integrable id μ) (hν : Integrable id ν) (h0 : ∫ ξ, ξ ∂ν = 0) :
    gaussianWidth K μ ≤ gaussianWidth K (μ ∗ ν) := by
  obtain ⟨R, hR⟩ := hK.isBounded.exists_norm_le
  have hint : Integrable (supportFn K) (μ ∗ ν) :=
    integrable_supportFn hne hR (integrable_norm_conv hμ.norm hν.norm)
  rw [gaussianWidth, gaussianWidth, integral_conv hint]
  have hprod : Integrable (fun p : E × E ↦ supportFn K (p.1 + p.2)) (μ.prod ν) := by
    have h := hint
    rwa [Measure.conv, integrable_map_measure hint.1 (by fun_prop)] at h
  refine integral_mono (integrable_supportFn hne hR hμ.norm) hprod.integral_prod_left fun a ↦ ?_
  obtain ⟨x, hx, hxa⟩ := exists_supportFn_eq_inner hK hne a
  have hb : Integrable (fun b ↦ supportFn K (a + b)) ν := by
    refine (((integrable_const ‖a‖).add hν.norm).const_mul R).mono'
      ((continuous_supportFn hne hR).comp (continuous_const.add continuous_id)).aestronglyMeasurable
      (ae_of_all _ fun b ↦ ?_)
    simp only [Real.norm_eq_abs]
    exact (abs_supportFn_le hne hR _).trans
      (mul_le_mul_of_nonneg_left (norm_add_le _ _) (nonneg_of_norm_le hne hR))
  have hab : Integrable (fun b ↦ a + b) ν := (integrable_const a).add hν
  calc supportFn K a = ⟪x, a⟫ := hxa
    _ = ∫ b, ⟪x, a + b⟫ ∂ν := by
        rw [integral_inner (𝕜 := ℝ) hab, integral_add (integrable_const a) (g := fun b ↦ b) hν,
          integral_const, h0]
        simp
    _ ≤ ∫ b, supportFn K (a + b) ∂ν :=
        integral_mono (hab.const_inner x) hb fun b ↦ inner_le_supportFn hR hx _

end GaussianWidth

section Multivariate

open scoped MatrixOrder

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {K : Set (EuclideanSpace ℝ ι)}
  {S₁ S₂ : Matrix ι ι ℝ}

/-- **Comparison by covariance domination**: if `S₁ ≤ S₂` in the Loewner order, the Gaussian
width of a compact set for `N(0, S₁)` is at most its width for `N(0, S₂)`. -/
lemma gaussianWidth_multivariateGaussian_le (hK : IsCompact K) (hne : K.Nonempty)
    (hS₁ : S₁.PosSemidef) (h : S₁ ≤ S₂) :
    gaussianWidth K (multivariateGaussian 0 S₁) ≤ gaussianWidth K (multivariateGaussian 0 S₂) := by
  have h2 : multivariateGaussian 0 S₂ =
      multivariateGaussian 0 S₁ ∗ multivariateGaussian 0 (S₂ - S₁) := by
    rw [multivariateGaussian_conv hS₁ (Matrix.le_iff.1 h)]
    simp
  rw [h2]
  exact gaussianWidth_le_gaussianWidth_conv hK hne IsGaussian.integrable_id
    IsGaussian.integrable_id (by simp)

end Multivariate

end ProbabilityTheory
