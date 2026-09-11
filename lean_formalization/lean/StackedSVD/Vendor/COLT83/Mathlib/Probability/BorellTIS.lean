/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Probability.MaureyPisier
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.LogSumExp
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.SupportFnDense
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianWidth
public import StackedSVD.Vendor.COLT83.Mathlib.Matrix.Loewner
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.MultivariateGaussian

set_option autoImplicit false

/-!
# Concentration of the supremum of a linear Gaussian process (Borell–TIS)

Let `μ` be the standard Gaussian measure on a finite-dimensional inner product space `E`.

* `hasSubgaussianMGF_iSup_inner_stdGaussian`: for finitely many vectors `y i` with `‖y i‖ ≤ σ`,
  the maximum `M v = max_i ⟪y i, v⟫` satisfies `HasSubgaussianMGF (M - μ[M]) (c_G σ²) μ`
  (Maurey–Pisier applied to the log-sum-exp smoothing of `M`, then `β → ∞`).
* `hasSubgaussianMGF_supportFn_stdGaussian`: for a nonempty compact `K` with `‖x‖ ≤ R` on `K`,
  the supremum `supportFn K v = sup_{x ∈ K} ⟪x, v⟫` satisfies
  `HasSubgaussianMGF (supportFn K - μ[supportFn K]) (c_G R²) μ` (dense sequence and dominated
  convergence).
* `measureReal_supportFn_sub_integral_ge_le`, `measureReal_supportFn_sub_integral_le_le`: the
  **Borell–TIS inequality** `P(M - E M ≥ u) ≤ exp (-u² / (2 c_G R²))` and its lower-tail
  version.
-/

@[expose] public section

open MeasureTheory ProbabilityTheory Real Filter Topology
open scoped RealInnerProductSpace NNReal

namespace ProbabilityTheory

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  [MeasurableSpace E] [BorelSpace E]

section finite

variable {ι : Type*} [Finite ι] [Nonempty ι] {y : ι → E} {σ : ℝ≥0}

/-- **Concentration of a finite maximum of Gaussian linear forms**: for vectors `y i` with
`‖y i‖ ≤ σ`, the maximum `M v = max_i ⟪y i, v⟫` satisfies `HasSubgaussianMGF (M - μ[M]) (c_G σ²)`
under the standard Gaussian measure `μ`. -/
lemma hasSubgaussianMGF_iSup_inner_stdGaussian (hσ : ∀ i, ‖y i‖ ≤ σ) :
    HasSubgaussianMGF (fun v ↦ (⨆ i, ⟪y i, v⟫) - ∫ w, (⨆ i, ⟪y i, w⟫) ∂stdGaussian E)
      (gaussianConcentrationConst * σ ^ 2) (stdGaussian E) := by
  have := Fintype.ofFinite ι
  set μ := stdGaussian E with hμ
  set M : E → ℝ := fun v ↦ ⨆ i, ⟪y i, v⟫ with hM
  have hMle : ∀ v, |M v| ≤ 0 + σ * ‖v‖ := fun v ↦ by
    rw [zero_add]
    exact abs_ciSup_inner_le hσ v
  have hMm : AEStronglyMeasurable M μ := (measurable_ciSup_inner y).aestronglyMeasurable
  have hMi : Integrable M μ := IsGaussian.integrable_of_abs_le_add_mul_norm hMm hMle
  have hexpM : ∀ t, Integrable (fun v ↦ exp (t * M v)) μ :=
    IsGaussian.integrable_exp_of_abs_le_add_mul_norm hMm hMle
  refine ⟨fun t ↦ ?_, fun t ↦ ?_⟩
  · refine ((hexpM t).const_mul (exp (-(t * ∫ w, M w ∂μ)))).congr
      (Filter.Eventually.of_forall fun v ↦ ?_)
    dsimp only
    simp only [hM]
    rw [← Real.exp_add]
    ring_nf
  -- the bound for the smoothed maximum, for every `β > 0`
  set B := exp (↑(gaussianConcentrationConst * σ ^ 2) * t ^ 2 / 2) with hB
  set m : ℝ := Real.log (Fintype.card ι) with hm
  have hm0 : 0 ≤ m := Real.log_nonneg (by exact_mod_cast Fintype.card_pos)
  have key : ∀ β : ℝ, 0 < β → mgf (fun v ↦ M v - ∫ w, M w ∂μ) μ t ≤ exp (|t| * m / β) * B := by
    intro β hβ
    have hsg := hasSubgaussianMGF_sub_integral_stdGaussian (contDiff_logSumExp β y)
      (norm_gradient_logSumExp_le hβ hσ)
    have hle : ∀ v, M v ≤ logSumExp β y v := fun v ↦ ciSup_inner_le_logSumExp hβ y v
    have hge : ∀ v, logSumExp β y v ≤ M v + m / β := fun v ↦
      logSumExp_le_ciSup_inner_add hβ y v
    have hfi : Integrable (logSumExp β y) μ :=
      IsGaussian.integrable_of_norm_gradient_le (contDiff_logSumExp β y)
        (norm_gradient_logSumExp_le hβ hσ)
    have hint1 : ∫ w, M w ∂μ ≤ ∫ w, logSumExp β y w ∂μ := integral_mono hMi hfi hle
    have hint2 : ∫ w, logSumExp β y w ∂μ ≤ (∫ w, M w ∂μ) + m / β := by
      calc ∫ w, logSumExp β y w ∂μ ≤ ∫ w, (M w + m / β) ∂μ :=
            integral_mono hfi (hMi.add (integrable_const _)) hge
        _ = (∫ w, M w ∂μ) + m / β := by
            rw [integral_add hMi (integrable_const _), integral_const, probReal_univ, one_smul]
    have hpt : ∀ v, t * (M v - ∫ w, M w ∂μ) ≤
        |t| * m / β + t * (logSumExp β y v - ∫ w, logSumExp β y w ∂μ) := fun v ↦ by
      have h1 := hle v
      have h2 := hge v
      rcases le_or_gt 0 t with ht | ht
      · rw [abs_of_nonneg ht]
        have := mul_le_mul_of_nonneg_left (show M v - logSumExp β y v +
          ((∫ w, logSumExp β y w ∂μ) - ∫ w, M w ∂μ) ≤ m / β by linarith) ht
        rw [mul_div_assoc]
        linarith
      · rw [abs_of_neg ht]
        have := mul_le_mul_of_nonneg_left (show -(m / β) ≤ M v - logSumExp β y v +
          ((∫ w, logSumExp β y w ∂μ) - ∫ w, M w ∂μ) by linarith) (neg_nonneg.2 ht.le)
        have hdiv : -t * m / β = -t * (m / β) := by ring
        rw [hdiv]
        linarith
    calc mgf (fun v ↦ M v - ∫ w, M w ∂μ) μ t
        = ∫ v, exp (t * (M v - ∫ w, M w ∂μ)) ∂μ := rfl
      _ ≤ ∫ v, exp (|t| * m / β) * exp (t * (logSumExp β y v - ∫ w, logSumExp β y w ∂μ)) ∂μ := by
          refine integral_mono_of_nonneg (Filter.Eventually.of_forall fun v ↦ (exp_pos _).le)
            ((hsg.integrable_exp_mul t).const_mul _) (Filter.Eventually.of_forall fun v ↦ ?_)
          dsimp only
          rw [← Real.exp_add]
          exact exp_le_exp.2 (hpt v)
      _ = exp (|t| * m / β) * mgf (fun v ↦ logSumExp β y v - ∫ w, logSumExp β y w ∂μ) μ t := by
          rw [integral_const_mul]
          rfl
      _ ≤ exp (|t| * m / β) * B := by
          gcongr
          exact hsg.mgf_le t
  -- let `β → ∞`
  have hlim : Tendsto (fun β : ℝ ↦ exp (|t| * m / β) * B) atTop (𝓝 B) := by
    have h1 : Tendsto (fun β : ℝ ↦ |t| * m / β) atTop (𝓝 0) :=
      tendsto_const_nhds.div_atTop tendsto_id
    have h2 := (Real.continuous_exp.tendsto 0).comp h1
    rw [Function.comp_def, Real.exp_zero] at h2
    simpa using h2.mul_const B
  exact ge_of_tendsto hlim ((eventually_gt_atTop 0).mono fun β hβ ↦ key β hβ)

end finite

section compact

variable {K : Set E} {R : ℝ≥0}

/-- **Concentration of the supremum of a Gaussian linear process over a compact set**: for a
nonempty compact `K` with `‖x‖ ≤ R` on `K`, `M v = sup_{x ∈ K} ⟪x, v⟫` satisfies
`HasSubgaussianMGF (M - μ[M]) (c_G R²)` under the standard Gaussian measure `μ`. -/
theorem hasSubgaussianMGF_supportFn_stdGaussian (hK : IsCompact K) (hne : K.Nonempty)
    (hR : ∀ x ∈ K, ‖x‖ ≤ R) :
    HasSubgaussianMGF (fun v ↦ supportFn K v - ∫ w, supportFn K w ∂stdGaussian E)
      (gaussianConcentrationConst * R ^ 2) (stdGaussian E) := by
  set μ := stdGaussian E with hμ
  obtain ⟨q, hq, hsup⟩ := exists_forall_supportFn_eq_iSup hK hne
  -- the supremum and the partial maxima
  set M : E → ℝ := supportFn K with hM
  set Mn : ℕ → E → ℝ := fun n ↦ partialSupInner q n with hMn
  have hMle : ∀ v, |M v| ≤ 0 + R * ‖v‖ := fun v ↦ by
    rw [zero_add]
    exact abs_supportFn_le hne hR v
  have hMnle : ∀ n v, |Mn n v| ≤ 0 + R * ‖v‖ := fun n v ↦ by
    rw [zero_add]
    exact abs_partialSupInner_le hq hR n v
  have hMm : AEStronglyMeasurable M μ := (continuous_supportFn hne hR).aestronglyMeasurable
  have hMnm : ∀ n, AEStronglyMeasurable (Mn n) μ := fun n ↦
    (measurable_ciSup_inner (fun i : Fin (n + 1) ↦ q i)).aestronglyMeasurable
  have hMi : Integrable M μ := IsGaussian.integrable_of_abs_le_add_mul_norm hMm hMle
  have hexpM : ∀ t, Integrable (fun v ↦ exp (t * M v)) μ :=
    IsGaussian.integrable_exp_of_abs_le_add_mul_norm hMm hMle
  have hlimM : ∀ v, Tendsto (fun n ↦ Mn n v) atTop (𝓝 (M v)) := fun v ↦
    tendsto_partialSupInner hq hR hsup v
  -- the bound for each partial maximum
  have hsg : ∀ n, HasSubgaussianMGF (fun v ↦ Mn n v - ∫ w, Mn n w ∂μ)
      (gaussianConcentrationConst * R ^ 2) μ := fun n ↦
    hasSubgaussianMGF_iSup_inner_stdGaussian (y := fun i : Fin (n + 1) ↦ q i) fun i ↦ hR _ (hq i)
  refine ⟨fun t ↦ ?_, fun t ↦ ?_⟩
  · refine ((hexpM t).const_mul (exp (-(t * ∫ w, M w ∂μ)))).congr
      (Filter.Eventually.of_forall fun v ↦ ?_)
    dsimp only
    rw [← Real.exp_add]
    ring_nf
  -- dominated convergence for the means and for the exponential moments
  have hbound : Integrable (fun v : E ↦ R * ‖v‖) μ := IsGaussian.integrable_id.norm.const_mul _
  have hmean : Tendsto (fun n ↦ ∫ w, Mn n w ∂μ) atTop (𝓝 (∫ w, M w ∂μ)) := by
    refine tendsto_integral_of_dominated_convergence (fun v ↦ R * ‖v‖) hMnm hbound
      (fun n ↦ Filter.Eventually.of_forall fun v ↦ ?_) (Filter.Eventually.of_forall hlimM)
    rw [Real.norm_eq_abs]
    simpa using hMnle n v
  have hexp : Tendsto (fun n ↦ ∫ v, exp (t * Mn n v) ∂μ) atTop (𝓝 (∫ v, exp (t * M v) ∂μ)) := by
    refine tendsto_integral_of_dominated_convergence (fun v ↦ exp (|t| * (R * ‖v‖)))
      (fun n ↦ (Real.continuous_exp.comp_aestronglyMeasurable ((hMnm n).const_mul t)))
      ?_ (fun n ↦ Filter.Eventually.of_forall fun v ↦ ?_)
      (Filter.Eventually.of_forall fun v ↦ ((Real.continuous_exp.tendsto _).comp
        ((hlimM v).const_mul t)))
    · have := IsGaussian.integrable_exp_mul_norm (μ := μ) (|t| * R)
      refine this.congr (Filter.Eventually.of_forall fun v ↦ ?_)
      dsimp only
      ring_nf
    · rw [Real.norm_of_nonneg (exp_pos _).le]
      refine exp_le_exp.2 ?_
      calc t * Mn n v ≤ |t * Mn n v| := le_abs_self _
        _ = |t| * |Mn n v| := abs_mul _ _
        _ ≤ |t| * (R * ‖v‖) := by
            gcongr
            simpa using hMnle n v
  -- the moment generating functions converge
  have hmgf : Tendsto (fun n ↦ mgf (fun v ↦ Mn n v - ∫ w, Mn n w ∂μ) μ t) atTop
      (𝓝 (mgf (fun v ↦ M v - ∫ w, M w ∂μ) μ t)) := by
    have h1 : ∀ n, mgf (fun v ↦ Mn n v - ∫ w, Mn n w ∂μ) μ t =
        exp (-(t * ∫ w, Mn n w ∂μ)) * ∫ v, exp (t * Mn n v) ∂μ := fun n ↦ by
      rw [mgf, ← integral_const_mul]
      refine integral_congr_ae (Filter.Eventually.of_forall fun v ↦ ?_)
      dsimp only
      rw [← Real.exp_add]
      ring_nf
    have h2 : mgf (fun v ↦ M v - ∫ w, M w ∂μ) μ t =
        exp (-(t * ∫ w, M w ∂μ)) * ∫ v, exp (t * M v) ∂μ := by
      rw [mgf, ← integral_const_mul]
      refine integral_congr_ae (Filter.Eventually.of_forall fun v ↦ ?_)
      dsimp only
      rw [← Real.exp_add]
      ring_nf
    simp_rw [h1, h2]
    exact ((Real.continuous_exp.tendsto _).comp (hmean.const_mul t).neg).mul hexp
  exact le_of_tendsto' hmgf fun n ↦ (hsg n).mgf_le t

/-- **Borell–TIS inequality** (upper tail): for a nonempty compact `K` with `‖x‖ ≤ R` on `K`
and `M v = sup_{x ∈ K} ⟪x, v⟫`, under the standard Gaussian measure `μ`,
`μ (M - μ[M] ≥ u) ≤ exp (-u² / (2 c_G R²))` for `u ≥ 0`. -/
lemma measureReal_supportFn_sub_integral_ge_le (hK : IsCompact K) (hne : K.Nonempty)
    (hR : ∀ x ∈ K, ‖x‖ ≤ R) {u : ℝ} (hu : 0 ≤ u) :
    (stdGaussian E).real {v | u ≤ supportFn K v - ∫ w, supportFn K w ∂stdGaussian E} ≤
      exp (-u ^ 2 / (2 * (gaussianConcentrationConst * R ^ 2))) :=
  (hasSubgaussianMGF_supportFn_stdGaussian hK hne hR).measure_ge_le hu

/-- **Borell–TIS inequality** (lower tail): `μ (M - μ[M] ≤ -u) ≤ exp (-u² / (2 c_G R²))`. -/
lemma measureReal_supportFn_sub_integral_le_le (hK : IsCompact K) (hne : K.Nonempty)
    (hR : ∀ x ∈ K, ‖x‖ ≤ R) {u : ℝ} (hu : 0 ≤ u) :
    (stdGaussian E).real {v | supportFn K v - ∫ w, supportFn K w ∂stdGaussian E ≤ -u} ≤
      exp (-u ^ 2 / (2 * (gaussianConcentrationConst * R ^ 2))) := by
  have h := (hasSubgaussianMGF_supportFn_stdGaussian hK hne hR).neg.measure_ge_le hu
  refine le_of_eq_of_le ?_ h
  congr 1
  ext v
  simp only [Set.mem_ofPred_eq, Pi.neg_apply]
  exact le_neg

end compact

end ProbabilityTheory

namespace ProbabilityTheory

section multivariate

open Matrix
open scoped MatrixOrder

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {K : Set (EuclideanSpace ℝ ι)}
  {S : Matrix ι ι ℝ} {σ : ℝ≥0}

/-- **Concentration of the supremum of a Gaussian linear process**, general covariance: for a
nonempty compact `K ⊆ ℝ^d` with `xᵀ S x ≤ σ²` on `K` and `G ~ N(0, S)`,
`M = sup_{x ∈ K} ⟪x, G⟫` satisfies `HasSubgaussianMGF (M - E M) (c_G σ²)`. -/
lemma hasSubgaussianMGF_supportFn_multivariateGaussian (hK : IsCompact K) (hne : K.Nonempty)
    (hS : S.PosSemidef) (hσ : ∀ x ∈ K, WithLp.ofLp x ⬝ᵥ S *ᵥ WithLp.ofLp x ≤ σ ^ 2) :
    HasSubgaussianMGF (fun v ↦ supportFn K v - ∫ w, supportFn K w ∂multivariateGaussian 0 S)
      (gaussianConcentrationConst * σ ^ 2) (multivariateGaussian 0 S) := by
  set L := toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt S) with hL
  have hLsymm : ∀ x ξ : EuclideanSpace ℝ ι, ⟪L x, ξ⟫ = ⟪x, L ξ⟫ := fun x ξ ↦ by
    rw [← ContinuousLinearMap.adjoint_inner_left, hL, adjoint_toEuclideanCLM, transpose_sqrt]
  have hmap : multivariateGaussian 0 S = (stdGaussian (EuclideanSpace ℝ ι)).map L := by
    rw [multivariateGaussian]
    congr
    funext x
    simp only [zero_add, hL]
  have himg : ∀ v, supportFn K (L v) = supportFn (L '' K) v := fun v ↦
    (supportFn_image K hLsymm v).symm
  have hK' : IsCompact (L '' K) := hK.image L.continuous
  have hne' : (L '' K).Nonempty := hne.image L
  have hR' : ∀ y ∈ L '' K, ‖y‖ ≤ σ := by
    rintro _ ⟨x, hx, rfl⟩
    have h := norm_toEuclideanCLM_sqrt_sq hS x
    rw [← hL] at h
    exact (pow_le_pow_iff_left₀ (norm_nonneg _) σ.2 two_ne_zero).1 (h.trans_le (hσ x hx))
  obtain ⟨R, hR⟩ : ∃ R, ∀ x ∈ K, ‖x‖ ≤ R := by
    obtain ⟨r, hr⟩ := hK.isBounded.subset_closedBall (0 : EuclideanSpace ℝ ι)
    exact ⟨r, fun x hx ↦ mem_closedBall_zero_iff.1 (hr hx)⟩
  have hstd := hasSubgaussianMGF_supportFn_stdGaussian hK' hne' hR'
  have hc : ∫ w, supportFn K w ∂multivariateGaussian 0 S =
      ∫ w, supportFn (L '' K) w ∂stdGaussian (EuclideanSpace ℝ ι) := by
    rw [hmap, integral_map L.continuous.measurable.aemeasurable
      (continuous_supportFn hne hR).aestronglyMeasurable]
    simp_rw [himg]
  rw [hc, hmap]
  set c := ∫ w, supportFn (L '' K) w ∂stdGaussian (EuclideanSpace ℝ ι) with hcdef
  have hXm : Measurable fun v ↦ supportFn K v - c :=
    ((continuous_supportFn hne hR).sub continuous_const).measurable
  have hXY : HasSubgaussianMGF ((fun v ↦ supportFn K v - c) ∘ L)
      (gaussianConcentrationConst * σ ^ 2) (stdGaussian (EuclideanSpace ℝ ι)) :=
    hstd.congr (Filter.Eventually.of_forall fun v ↦ by simp [himg])
  have h1 := (HasSubgaussianMGF.id_map_iff hXY.aemeasurable).2 hXY
  rw [← Measure.map_map hXm L.continuous.measurable] at h1
  exact (HasSubgaussianMGF.id_map_iff hXm.aemeasurable).1 h1

end multivariate

end ProbabilityTheory
