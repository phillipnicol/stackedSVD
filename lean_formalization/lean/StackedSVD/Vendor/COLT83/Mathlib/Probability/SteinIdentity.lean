/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianMGF
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinReal
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.MultivariateGaussian
public import StackedSVD.Vendor.COLT83.Mathlib.Matrix.Loewner
public import Mathlib.Analysis.Calculus.FDeriv.Equiv
public import Mathlib.Analysis.Calculus.FDeriv.CompCLM
public import Mathlib.Analysis.Calculus.Deriv.Prod
public import Mathlib.Analysis.Calculus.Deriv.Comp
public import Mathlib.MeasureTheory.Integral.Prod

set_option autoImplicit false

/-!
# Stein's identity for the standard Gaussian measure

For the standard Gaussian measure on a finite-dimensional inner product space `E` and a `C¹`
function `F : E → ℝ` with bounded derivative,
`∫ ⟪a, x⟫ F x ∂(stdGaussian E) = ∫ DF(x) a ∂(stdGaussian E)`
(`integral_inner_mul_stdGaussian`, Gaussian integration by parts). It is deduced from the
one-dimensional identity (`integral_mul_gaussianReal_zero_one_of_hasDerivAt`) by Fubini on the
coordinates (`integral_apply_mul_stdGaussian`).

`integral_fderiv_apply_stdGaussian` is the form used by the Gaussian interpolation method: for
`C²` functions with bounded second derivative and linear maps `L L' : E' → E`,
`∫ DF(L g) (L' g) = ∫ ∑ k, D²F(L g) (L (b k)) (L' (b k))` for an orthonormal basis `b` of `E'`.
For a centered Gaussian vector with covariance `S`, `integral_inner_mul_multivariateGaussian`
gives `E[⟪a, X⟫ F(X)] = E[DF(X) (S a)]`.
-/

@[expose] public section

open MeasureTheory Real InnerProductSpace
open scoped RealInnerProductSpace MatrixOrder

namespace ProbabilityTheory

section integrability

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E] {μ : Measure E} [IsGaussian μ]
  {F : E → ℝ} {L : ℝ}

/-- `⟪a, x⟫ F x` is integrable under a Gaussian measure when `F` is `C¹` with bounded
derivative. -/
lemma IsGaussian.integrable_inner_mul_of_norm_fderiv_le (hF : ContDiff ℝ 1 F)
    (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (a : E) :
    Integrable (fun x ↦ ⟪a, x⟫ * F x) μ := by
  refine IsGaussian.integrable_of_abs_le_add_mul_norm_sq (A := 0) (B := ‖a‖ * |F 0|)
    (C := ‖a‖ * L) ((continuous_const.inner continuous_id).mul hF.continuous).aestronglyMeasurable
    fun x ↦ ?_
  rw [abs_mul]
  calc |⟪a, x⟫| * |F x| ≤ (‖a‖ * ‖x‖) * (|F 0| + L * ‖x‖) :=
        mul_le_mul (abs_real_inner_le_norm _ _)
          (abs_le_of_norm_fderiv_le (hF.differentiable one_ne_zero) hL x) (abs_nonneg _)
          (by positivity)
    _ = 0 + ‖a‖ * |F 0| * ‖x‖ + ‖a‖ * L * ‖x‖ ^ 2 := by ring

end integrability

section integrabilityNormed

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E] {μ : Measure E} [IsGaussian μ]
  {F : E → ℝ} {L : ℝ}

/-- `DF(x) a` is integrable under a Gaussian measure when `F` is `C¹` with bounded
derivative. -/
lemma IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le (hF : ContDiff ℝ 1 F)
    (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (a : E) :
    Integrable (fun x ↦ fderiv ℝ F x a) μ :=
  IsGaussian.integrable_of_abs_le_add_mul_norm (A := L * ‖a‖) (B := 0)
    ((hF.continuous_fderiv one_ne_zero).clm_apply continuous_const).aestronglyMeasurable fun x ↦ by
      rw [← Real.norm_eq_abs, zero_mul, add_zero]
      exact (ContinuousLinearMap.le_opNorm _ _).trans
        (mul_le_mul_of_nonneg_right (hL x) (norm_nonneg _))

end integrabilityNormed

section euclidean

variable {n : ℕ} {F : EuclideanSpace ℝ (Fin (n + 1)) → ℝ} {L : ℝ}

/-- The derivative of `t ↦ (y with `t` inserted at position `i`)`. -/
lemma hasDerivAt_toLp_insertNth (i : Fin (n + 1)) (y : Fin n → ℝ) (t : ℝ) :
    HasDerivAt (fun t : ℝ ↦ (WithLp.toLp 2 (Fin.insertNth i t y) : EuclideanSpace ℝ (Fin (n + 1))))
      (EuclideanSpace.single i 1) t := by
  have h : HasDerivAt (fun t : ℝ ↦ (Fin.insertNth i t y : Fin (n + 1) → ℝ)) (Pi.single i 1) t := by
    refine hasDerivAt_pi.2 fun j ↦ ?_
    rcases Fin.eq_self_or_eq_succAbove i j with rfl | ⟨k, rfl⟩
    · simp only [Fin.insertNth_apply_same, Pi.single_eq_same]
      exact hasDerivAt_id t
    · simp only [Fin.insertNth_apply_succAbove, Pi.single_eq_of_ne (Fin.succAbove_ne i k)]
      exact hasDerivAt_const t _
  exact (PiLp.continuousLinearEquiv 2 ℝ _).symm.hasFDerivAt.comp_hasDerivAt t h

/-- Fubini for the standard Gaussian measure on `ℝ^{n+1}`: integrate first over the `i`-th
coordinate. -/
lemma integral_stdGaussian_eq_integral_insertNth {G : EuclideanSpace ℝ (Fin (n + 1)) → ℝ}
    (hG : Integrable G (stdGaussian (EuclideanSpace ℝ (Fin (n + 1))))) (i : Fin (n + 1)) :
    ∫ x, G x ∂stdGaussian (EuclideanSpace ℝ (Fin (n + 1))) =
      ∫ y, ∫ t, G (WithLp.toLp 2 (Fin.insertNth i t y)) ∂gaussianReal 0 1
        ∂(Measure.pi fun _ : Fin n ↦ gaussianReal 0 1) := by
  have hmp := (measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) ↦ gaussianReal 0 1) i).symm
  rw [← map_pi_eq_stdGaussian] at hG ⊢
  rw [integral_map (by fun_prop) hG.1, ← hmp.integral_comp']
  have hint : Integrable (fun z ↦ G (WithLp.toLp 2 ((MeasurableEquiv.piFinSuccAbove _ i).symm z)))
      ((gaussianReal 0 1).prod (Measure.pi fun _ : Fin n ↦ gaussianReal 0 1)) := by
    rw [show (fun z ↦ G (WithLp.toLp 2 ((MeasurableEquiv.piFinSuccAbove _ i).symm z))) =
      (G ∘ WithLp.toLp 2) ∘ (MeasurableEquiv.piFinSuccAbove _ i).symm from rfl,
      hmp.integrable_comp_emb (MeasurableEquiv.measurableEmbedding _)]
    exact (integrable_map_measure hG.1 (by fun_prop)).1 hG
  rw [integral_prod_symm _ hint]
  simp only [MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv, Equiv.coe_fn_mk]

/-- **Stein's identity** for the `i`-th coordinate of a standard Gaussian vector:
`E[g_i F(g)] = E[∂_i F(g)]` for `C¹` functions with bounded derivative. -/
lemma integral_apply_mul_stdGaussian (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L)
    (i : Fin (n + 1)) :
    ∫ x, x i * F x ∂stdGaussian (EuclideanSpace ℝ (Fin (n + 1))) =
      ∫ x, fderiv ℝ F x (EuclideanSpace.single i 1)
        ∂stdGaussian (EuclideanSpace ℝ (Fin (n + 1))) := by
  have hint1 : Integrable (fun x ↦ x i * F x) (stdGaussian (EuclideanSpace ℝ (Fin (n + 1)))) := by
    have := IsGaussian.integrable_inner_mul_of_norm_fderiv_le
      (μ := stdGaussian (EuclideanSpace ℝ (Fin (n + 1)))) hF hL (EuclideanSpace.single i 1)
    simpa [EuclideanSpace.inner_single_left] using this
  have hint2 : Integrable (fun x ↦ fderiv ℝ F x (EuclideanSpace.single i 1))
      (stdGaussian (EuclideanSpace ℝ (Fin (n + 1)))) :=
    IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le hF hL _
  rw [integral_stdGaussian_eq_integral_insertNth hint1 i,
    integral_stdGaussian_eq_integral_insertNth hint2 i]
  congr 1
  ext y
  simp only [Fin.insertNth_apply_same]
  refine integral_mul_gaussianReal_zero_one_of_hasDerivAt (L := L)
    (f' := fun t ↦ fderiv ℝ F (WithLp.toLp 2 (Fin.insertNth i t y)) (EuclideanSpace.single i 1))
    (fun t ↦ ((hF.differentiable one_ne_zero) _).hasFDerivAt.comp_hasDerivAt t
      (hasDerivAt_toLp_insertNth i y t)) fun t ↦ ?_
  rw [← Real.norm_eq_abs]
  calc ‖fderiv ℝ F (WithLp.toLp 2 (Fin.insertNth i t y)) (EuclideanSpace.single i 1)‖
      ≤ ‖fderiv ℝ F (WithLp.toLp 2 (Fin.insertNth i t y))‖ * ‖(EuclideanSpace.single i 1 :
          EuclideanSpace ℝ (Fin (n + 1)))‖ := ContinuousLinearMap.le_opNorm _ _
    _ ≤ L := by rw [PiLp.norm_single, norm_one, mul_one]; exact hL _

/-- **Stein's identity** on `ℝ^{n+1}`: `E[⟪a, g⟫ F(g)] = E[DF(g) a]`. -/
lemma integral_inner_mul_stdGaussian_euclidean (hF : ContDiff ℝ 1 F)
    (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (a : EuclideanSpace ℝ (Fin (n + 1))) :
    ∫ x, ⟪a, x⟫ * F x ∂stdGaussian (EuclideanSpace ℝ (Fin (n + 1))) =
      ∫ x, fderiv ℝ F x a ∂stdGaussian (EuclideanSpace ℝ (Fin (n + 1))) := by
  have ha : a = ∑ j, a j • (EuclideanSpace.single j 1 : EuclideanSpace ℝ (Fin (n + 1))) := by
    conv_lhs => rw [← (EuclideanSpace.basisFun (Fin (n + 1)) ℝ).sum_repr a]
    simp [EuclideanSpace.basisFun_apply, EuclideanSpace.basisFun_repr]
  have h1 : ∀ x : EuclideanSpace ℝ (Fin (n + 1)), ⟪a, x⟫ * F x = ∑ j, a j * (x j * F x) := by
    intro x
    conv_lhs => rw [ha]
    simp [sum_inner, real_inner_smul_left, EuclideanSpace.inner_single_left, Finset.sum_mul,
      mul_assoc]
  have h2 : ∀ x : EuclideanSpace ℝ (Fin (n + 1)),
      fderiv ℝ F x a = ∑ j, a j * fderiv ℝ F x (EuclideanSpace.single j 1) := by
    intro x
    conv_lhs => rw [ha]
    simp [map_sum, map_smul]
  simp_rw [h1, h2]
  rw [integral_finsetSum _ fun j _ ↦ ?_, integral_finsetSum _ fun j _ ↦ ?_]
  · refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [integral_const_mul, integral_const_mul, integral_apply_mul_stdGaussian hF hL j]
  · exact (IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le hF hL _).const_mul _
  · have := IsGaussian.integrable_inner_mul_of_norm_fderiv_le
      (μ := stdGaussian (EuclideanSpace ℝ (Fin (n + 1)))) hF hL (EuclideanSpace.single j 1)
    simpa [EuclideanSpace.inner_single_left] using this.const_mul (a j)

end euclidean

section general

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  [MeasurableSpace E] [BorelSpace E] {F : E → ℝ} {L : ℝ}

/-- **Stein's identity** (Gaussian integration by parts): for the standard Gaussian measure on a
finite-dimensional inner product space `E` and a `C¹` function `F` with bounded derivative,
`E[⟪a, g⟫ F(g)] = E[DF(g) a]`. -/
lemma integral_inner_mul_stdGaussian (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L)
    (a : E) :
    ∫ x, ⟪a, x⟫ * F x ∂stdGaussian E = ∫ x, fderiv ℝ F x a ∂stdGaussian E := by
  obtain ⟨n, hn⟩ : ∃ n, Module.finrank ℝ E = n := ⟨_, rfl⟩
  rcases n with _ | n
  · have ha : a = 0 := finrank_zero_iff_forall_zero.1 hn a
    simp [ha]
  have hL0 : 0 ≤ L := le_trans (norm_nonneg _) (hL 0)
  let b : OrthonormalBasis (Fin (n + 1)) ℝ E := (stdOrthonormalBasis ℝ E).reindex (finCongr hn)
  let e : EuclideanSpace ℝ (Fin (n + 1)) ≃ₗᵢ[ℝ] E := b.repr.symm
  have h1 : ∀ y, ⟪a, e y⟫ = ⟪b.repr a, y⟫ := fun y ↦ by
    rw [← b.repr.inner_map_map, LinearIsometryEquiv.apply_symm_apply]
  have h2 : ∀ y, fderiv ℝ F (e y) a = fderiv ℝ (F ∘ e) y (b.repr a) := fun y ↦ by
    rw [show (F ∘ e) = F ∘ e.toContinuousLinearEquiv from rfl,
      e.toContinuousLinearEquiv.comp_right_fderiv]
    simp [e]
  have hFe : ContDiff ℝ 1 (F ∘ e) := hF.comp e.contDiff
  have hLe : ∀ y, ‖fderiv ℝ (F ∘ e) y‖ ≤ L := fun y ↦ by
    rw [show (F ∘ e) = F ∘ e.toContinuousLinearEquiv from rfl,
      e.toContinuousLinearEquiv.comp_right_fderiv]
    refine ContinuousLinearMap.opNorm_le_bound _ hL0 fun w ↦ ?_
    rw [ContinuousLinearMap.comp_apply]
    calc ‖fderiv ℝ F (e.toContinuousLinearEquiv y) (e.toContinuousLinearEquiv w)‖
        ≤ ‖fderiv ℝ F (e.toContinuousLinearEquiv y)‖ * ‖e.toContinuousLinearEquiv w‖ :=
          ContinuousLinearMap.le_opNorm _ _
      _ ≤ L * ‖w‖ := by
          rw [show e.toContinuousLinearEquiv w = e w from rfl, e.norm_map]
          gcongr
          exact hL _
  rw [← stdGaussian_map e, integral_map (by fun_prop) ?_, integral_map (by fun_prop) ?_]
  · simp_rw [h1, h2]
    exact integral_inner_mul_stdGaussian_euclidean hFe hLe (b.repr a)
  · exact ((hF.continuous_fderiv one_ne_zero).clm_apply continuous_const).aestronglyMeasurable
  · exact ((continuous_const.inner continuous_id).mul hF.continuous).aestronglyMeasurable

variable {E' : Type*} [NormedAddCommGroup E'] [InnerProductSpace ℝ E'] [FiniteDimensional ℝ E']
  [MeasurableSpace E'] [BorelSpace E'] {κ : Type*} [Fintype κ]

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
/-- **Stein's identity for linear images**: if `F : E → ℝ` is `C²` with bounded second
derivative and `L L' : E' → E` are linear, then for `g ~ N(0, I)` on `E'` and an orthonormal
basis `b` of `E'`, `E[DF(L g) (L' g)] = E[∑ k, D²F(L g) (L (b k)) (L' (b k))]`. -/
lemma integral_fderiv_apply_stdGaussian (hF : ContDiff ℝ 2 F) {K : ℝ}
    (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K) (L L' : E' →L[ℝ] E) (b : OrthonormalBasis κ ℝ E') :
    ∫ g, fderiv ℝ F (L g) (L' g) ∂stdGaussian E' =
      ∫ g, ∑ k, fderiv ℝ (fderiv ℝ F) (L g) (L (b k)) (L' (b k)) ∂stdGaussian E' := by
  have hK0 : 0 ≤ K := le_trans (norm_nonneg (fderiv ℝ (fderiv ℝ F) 0)) (hK 0)
  have hF1 : ContDiff ℝ 1 (fderiv ℝ F) := hF.fderiv_right (m := 1) le_rfl
  have hφd : ∀ (k : κ) (g : E'), HasFDerivAt (fun g ↦ fderiv ℝ F (L g) (L' (b k)))
      (((fderiv ℝ (fderiv ℝ F) (L g)).comp L).flip (L' (b k))) g := fun k g ↦ by
    have h1 : HasFDerivAt (fun g ↦ fderiv ℝ F (L g)) ((fderiv ℝ (fderiv ℝ F) (L g)).comp L) g :=
      ((hF1.differentiable one_ne_zero) (L g)).hasFDerivAt.comp g L.hasFDerivAt
    exact (h1.clm_apply (hasFDerivAt_const (L' (b k)) g)).congr_fderiv (by simp)
  have hφc : ∀ k : κ, ContDiff ℝ 1 (fun g ↦ fderiv ℝ F (L g) (L' (b k))) := fun k ↦
    (hF1.comp L.contDiff).clm_apply contDiff_const
  have hφL : ∀ (k : κ) (g : E'), ‖fderiv ℝ (fun g ↦ fderiv ℝ F (L g) (L' (b k))) g‖ ≤
      K * ‖L‖ * ‖L' (b k)‖ := fun k g ↦ by
    rw [(hφd k g).fderiv]
    refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun w ↦ ?_
    rw [ContinuousLinearMap.flip_apply, ContinuousLinearMap.comp_apply]
    calc ‖fderiv ℝ (fderiv ℝ F) (L g) (L w) (L' (b k))‖
        ≤ ‖fderiv ℝ (fderiv ℝ F) (L g) (L w)‖ * ‖L' (b k)‖ := ContinuousLinearMap.le_opNorm _ _
      _ ≤ ‖fderiv ℝ (fderiv ℝ F) (L g)‖ * ‖L w‖ * ‖L' (b k)‖ := by
          gcongr
          exact ContinuousLinearMap.le_opNorm _ _
      _ ≤ K * (‖L‖ * ‖w‖) * ‖L' (b k)‖ := by
          gcongr
          · exact hK _
          · exact L.le_opNorm w
      _ = K * ‖L‖ * ‖L' (b k)‖ * ‖w‖ := by ring
  have hexp : ∀ g, fderiv ℝ F (L g) (L' g) = ∑ k, ⟪b k, g⟫ * fderiv ℝ F (L g) (L' (b k)) := by
    intro g
    have hg : L' g = ∑ k, ⟪b k, g⟫ • L' (b k) := by
      conv_lhs => rw [← b.sum_repr' g]
      simp [map_sum, map_smul]
    rw [hg, map_sum]
    simp [map_smul]
  simp_rw [hexp]
  rw [integral_finsetSum _ fun k _ ↦
    IsGaussian.integrable_inner_mul_of_norm_fderiv_le (hφc k) (hφL k) (b k)]
  rw [integral_finsetSum _ fun k _ ↦ ?_]
  · refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [integral_inner_mul_stdGaussian (hφc k) (hφL k) (b k)]
    congr 1
    ext g
    rw [(hφd k g).fderiv, ContinuousLinearMap.flip_apply, ContinuousLinearMap.comp_apply]
  · have := IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le (μ := stdGaussian E') (hφc k)
      (hφL k) (b k)
    refine this.congr (Filter.Eventually.of_forall fun g ↦ ?_)
    dsimp only
    rw [(hφd k g).fderiv, ContinuousLinearMap.flip_apply, ContinuousLinearMap.comp_apply]

end general

section multivariate

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {F : EuclideanSpace ℝ ι → ℝ} {L : ℝ}

/-- **Stein's identity** for a centered Gaussian vector `X ~ N(0, S)` of `ℝ^ι`:
`E[⟪a, X⟫ F(X)] = E[DF(X) (S a)]`, i.e. `E[X_i F(X)] = ∑ j, S i j E[∂_j F(X)]`. -/
lemma integral_inner_mul_multivariateGaussian {S : Matrix ι ι ℝ} (hS : S.PosSemidef)
    (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (a : EuclideanSpace ℝ ι) :
    ∫ x, ⟪a, x⟫ * F x ∂multivariateGaussian 0 S =
      ∫ x, fderiv ℝ F x (Matrix.toEuclideanCLM (𝕜 := ℝ) S a) ∂multivariateGaussian 0 S := by
  have hL0 : 0 ≤ L := le_trans (norm_nonneg _) (hL 0)
  set M := Matrix.toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt S) with hM
  have hmap : multivariateGaussian 0 S = (stdGaussian (EuclideanSpace ℝ ι)).map M := by
    rw [← multivariateGaussian_zero_one,
      multivariateGaussian_zero_map_toEuclideanCLM Matrix.PosSemidef.one, Matrix.mul_one,
      Matrix.transpose_sqrt, hS.sqrt_mul_sqrt]
  have hsymm : ∀ x y, ⟪x, M y⟫ = ⟪M x, y⟫ := fun x y ↦ by
    rw [← ContinuousLinearMap.adjoint_inner_left, hM, Matrix.adjoint_toEuclideanCLM,
      Matrix.transpose_sqrt]
  have hMM : ∀ x, M (M x) = Matrix.toEuclideanCLM (𝕜 := ℝ) S x := fun x ↦ by
    change (M * M) x = _
    rw [hM, ← map_mul, hS.sqrt_mul_sqrt]
  have hFM : ContDiff ℝ 1 (F ∘ M) := hF.comp M.contDiff
  have hdM : ∀ g, fderiv ℝ (F ∘ M) g = (fderiv ℝ F (M g)).comp M := fun g ↦ by
    rw [fderiv_comp g (hF.differentiable one_ne_zero _) M.differentiableAt, M.fderiv]
  have hLM : ∀ g, ‖fderiv ℝ (F ∘ M) g‖ ≤ L * ‖M‖ := fun g ↦ by
    rw [hdM]
    exact (ContinuousLinearMap.opNorm_comp_le _ _).trans
      (mul_le_mul_of_nonneg_right (hL _) (norm_nonneg _))
  rw [hmap, integral_map (by fun_prop) ?_, integral_map (by fun_prop) ?_]
  · simp_rw [hsymm a]
    have h := integral_inner_mul_stdGaussian hFM hLM (M a)
    simp_rw [Function.comp_apply, hdM, ContinuousLinearMap.comp_apply, hMM] at h
    exact h
  · exact ((hF.continuous_fderiv one_ne_zero).clm_apply continuous_const).aestronglyMeasurable
  · exact ((continuous_const.inner continuous_id).mul hF.continuous).aestronglyMeasurable

end multivariate

end ProbabilityTheory
