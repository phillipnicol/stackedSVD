/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.WithDensityPi
import StackedSVD.Prob.GaussianMatrix
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.MeasureTheory.Measure.CharacteristicFunction.Basic

/-!
# The Lebesgue density of a linear image of the standard Gaussian (L5, units U1 and U3)

`gaussDensity d S` is the density of `N(0, S)` on `Fin d → ℝ`. The file proves that the
standard Gaussian product measure `Measure.pi fun _ : Fin d => gaussianReal 0 1` has the
product density (U1, from `gaussianReal_of_var_ne_zero` and `pi_withDensity`), that its
image under an invertible matrix `A` has density `gaussDensity d (A * Aᵀ)` (U3, from
`withDensity_map_mulVec`), and the characteristic function of that image (U3a, from
`charFun_pi` and `charFun_gaussianReal`). No `CFC.sqrt` and no `multivariateGaussian`
enters: the square root of the covariance is explicit in `MLEMarginal/Defs.lean`.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped ENNReal Matrix

namespace StackedSVD

/-- The density of `N(0, S)` on `Fin d → ℝ` for a positive definite `S`:
`exp(-xᵀ S⁻¹ x / 2) / √((2π)^d det S)`. -/
noncomputable def gaussDensity (d : ℕ) (S : Matrix (Fin d) (Fin d) ℝ) (x : Fin d → ℝ) : ℝ :=
  Real.exp (-(x ⬝ᵥ S⁻¹ *ᵥ x) / 2) / Real.sqrt ((2 * Real.pi) ^ d * S.det)

/-- `√(x ^ n) = (√x) ^ n` for a nonnegative `x`. Mathlib v4.33.0 has no `Real.sqrt_pow`. -/
theorem sqrt_pow_of_nonneg {x : ℝ} (hx : 0 ≤ x) (n : ℕ) :
    Real.sqrt (x ^ n) = Real.sqrt x ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [pow_succ, Real.sqrt_mul (pow_nonneg hx n), ih, pow_succ]

/-- The quadratic form of a linear image: `(M x) ⬝ᵥ (M x) = x ⬝ᵥ (Mᵀ M) x`. -/
theorem dotProduct_mulVec_self {d : ℕ} (M : Matrix (Fin d) (Fin d) ℝ) (x : Fin d → ℝ) :
    (M *ᵥ x) ⬝ᵥ (M *ᵥ x) = x ⬝ᵥ (Mᵀ * M) *ᵥ x := by
  rw [Matrix.dotProduct_mulVec, Matrix.vecMul_mulVec, Matrix.dotProduct_mulVec]

/-- The product of `d` standard Gaussian densities, in closed form. -/
theorem prod_gaussianPDFReal {d : ℕ} (w : Fin d → ℝ) :
    ∏ i, gaussianPDFReal 0 1 (w i)
      = Real.exp (-(w ⬝ᵥ w) / 2) / Real.sqrt ((2 * Real.pi) ^ d) := by
  have hval : ∀ y : ℝ, gaussianPDFReal 0 1 y
      = (Real.sqrt (2 * Real.pi))⁻¹ * Real.exp (-(y ^ 2) / 2) := by
    intro y
    simp only [gaussianPDFReal, NNReal.coe_one, mul_one, sub_zero]
  have h1 : w ⬝ᵥ w = ∑ i : Fin d, w i ^ 2 := by
    rw [dotProduct]
    exact Finset.sum_congr rfl fun i _ => (sq (w i)).symm
  have hsum : ∑ i : Fin d, (-(w i ^ 2) / 2) = -(w ⬝ᵥ w) / 2 := by
    rw [h1]
    simp only [Finset.sum_neg_distrib, ← Finset.sum_div]
  calc ∏ i, gaussianPDFReal 0 1 (w i)
      = ∏ i : Fin d, ((Real.sqrt (2 * Real.pi))⁻¹ * Real.exp (-(w i ^ 2) / 2)) :=
        Finset.prod_congr rfl fun i _ => hval (w i)
    _ = (Real.sqrt (2 * Real.pi))⁻¹ ^ d * Real.exp (∑ i : Fin d, (-(w i ^ 2) / 2)) := by
        rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.card_fin, ← Real.exp_sum]
    _ = (Real.sqrt (2 * Real.pi))⁻¹ ^ d * Real.exp (-(w ⬝ᵥ w) / 2) := by rw [hsum]
    _ = Real.exp (-(w ⬝ᵥ w) / 2) / Real.sqrt ((2 * Real.pi) ^ d) := by
        rw [sqrt_pow_of_nonneg (by positivity : (0:ℝ) ≤ 2 * Real.pi) d, inv_pow]
        ring

theorem gaussDensity_pos {d : ℕ} {S : Matrix (Fin d) (Fin d) ℝ} (hS : 0 < S.det)
    (x : Fin d → ℝ) : 0 < gaussDensity d S x := by
  have hpi : (0 : ℝ) < (2 * Real.pi) ^ d := by positivity
  exact div_pos (Real.exp_pos _) (Real.sqrt_pos.mpr (mul_pos hpi hS))

theorem measurable_gaussDensity (d : ℕ) (S : Matrix (Fin d) (Fin d) ℝ) :
    Measurable (gaussDensity d S) := by
  have hlin : Continuous fun x : Fin d → ℝ => S⁻¹ *ᵥ x :=
    (Matrix.mulVecLin S⁻¹).continuous_of_finiteDimensional
  have hq : Measurable fun x : Fin d → ℝ => x ⬝ᵥ S⁻¹ *ᵥ x := by
    have hc : Continuous fun x : Fin d → ℝ => x ⬝ᵥ S⁻¹ *ᵥ x := by
      simp only [dotProduct]
      exact continuous_finsetSum _ fun i _ =>
        (continuous_apply i).mul ((continuous_apply i).comp hlin)
    exact hc.measurable
  change Measurable fun x : Fin d → ℝ =>
    Real.exp (-(x ⬝ᵥ S⁻¹ *ᵥ x) / 2) / Real.sqrt ((2 * Real.pi) ^ d * S.det)
  exact ((hq.neg.div_const 2).exp).div_const _

/-- U1: the standard Gaussian product measure on `Fin d → ℝ` has the product density
(`gaussianReal_of_var_ne_zero`, `pi_withDensity`, `ENNReal.ofReal_prod_of_nonneg`). -/
theorem pi_gaussianReal_eq_withDensity (d : ℕ) :
    (Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (∏ i, gaussianPDFReal 0 1 (x i))) := by
  have hv : (1 : NNReal) ≠ 0 := one_ne_zero
  have hprob : ∀ _ : Fin d,
      IsProbabilityMeasure ((volume : Measure ℝ).withDensity (gaussianPDF 0 1)) := by
    intro _
    rw [← gaussianReal_of_var_ne_zero (0 : ℝ) hv]
    infer_instance
  have hfun : (fun _ : Fin d => gaussianReal 0 1)
      = fun _ : Fin d => (volume : Measure ℝ).withDensity (gaussianPDF 0 1) :=
    funext fun _ => gaussianReal_of_var_ne_zero 0 hv
  have key := pi_withDensity (fun _ : Fin d => (volume : Measure ℝ))
      (fun _ : Fin d => gaussianPDF 0 1) (fun _ => measurable_gaussianPDF 0 1)
  calc (Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = Measure.pi (fun _ : Fin d => (volume : Measure ℝ).withDensity (gaussianPDF 0 1)) := by
        rw [hfun]
    _ = (Measure.pi fun _ : Fin d => (volume : Measure ℝ)).withDensity
          (fun x => ∏ i, gaussianPDF 0 1 (x i)) := key
    _ = (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ∏ i, gaussianPDF 0 1 (x i)) := rfl
    _ = (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (∏ i, gaussianPDFReal 0 1 (x i))) := by
        congr 1
        funext x
        rw [ENNReal.ofReal_prod_of_nonneg fun i _ => gaussianPDFReal_nonneg 0 1 (x i)]
        rfl

/-- U3: the image of the standard Gaussian under an invertible matrix `A` has density
`gaussDensity d (A * Aᵀ)`. Pointwise: `|det A|⁻¹ ∏ᵢ φ((A⁻¹x)ᵢ) = gaussDensity d (A Aᵀ) x`,
since `(A⁻¹x) ⬝ᵥ (A⁻¹x) = x ⬝ᵥ (A Aᵀ)⁻¹ x` and `det (A Aᵀ) = (det A)²`. -/
theorem map_mulVec_pi_gaussianReal {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0) :
    (Measure.pi fun _ : Fin d => gaussianReal 0 1).map (fun z => A *ᵥ z)
      = (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (gaussDensity d (A * Aᵀ) x)) := by
  have hmeas : Measurable fun x : Fin d → ℝ =>
      ENNReal.ofReal (∏ i, gaussianPDFReal 0 1 (x i)) :=
    Measurable.ennreal_ofReal (Finset.measurable_prod _ fun i _ =>
      (measurable_gaussianPDFReal 0 1).comp (measurable_pi_apply i))
  have hdet : (A * Aᵀ).det = A.det ^ 2 := by
    rw [Matrix.det_mul, Matrix.det_transpose, sq]
  have hsq : Real.sqrt ((2 * Real.pi) ^ d * (A * Aᵀ).det)
      = Real.sqrt ((2 * Real.pi) ^ d) * |A.det| := by
    rw [hdet, Real.sqrt_mul (by positivity), Real.sqrt_sq_eq_abs]
  rw [pi_gaussianReal_eq_withDensity, withDensity_map_mulVec A hA hmeas]
  congr 1
  funext x
  show ENNReal.ofReal |A.det|⁻¹ * ENNReal.ofReal (∏ i, gaussianPDFReal 0 1 ((A⁻¹ *ᵥ x) i))
      = ENNReal.ofReal (gaussDensity d (A * Aᵀ) x)
  have hq : (A⁻¹ *ᵥ x) ⬝ᵥ (A⁻¹ *ᵥ x) = x ⬝ᵥ (A * Aᵀ)⁻¹ *ᵥ x := by
    rw [dotProduct_mulVec_self, Matrix.transpose_nonsing_inv, Matrix.mul_inv_rev]
  have hreal : |A.det|⁻¹ * (∏ i, gaussianPDFReal 0 1 ((A⁻¹ *ᵥ x) i))
      = gaussDensity d (A * Aᵀ) x := by
    rw [prod_gaussianPDFReal, hq]
    simp only [gaussDensity]
    rw [hsq]
    ring
  rw [← ENNReal.ofReal_mul (by positivity), hreal]

/-- U3a with the frequency written as `WithLp.toLp 2 u` for a plain vector `u`. -/
theorem charFun_map_mulVec_pi_gaussianReal' {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ)
    (u : Fin d → ℝ) :
    charFun ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
        (fun z => (WithLp.toLp 2 (A *ᵥ z) : EuclideanSpace ℝ (Fin d))))
        (WithLp.toLp 2 u : EuclideanSpace ℝ (Fin d))
      = Complex.exp (-(((u ⬝ᵥ (A * Aᵀ) *ᵥ u : ℝ) : ℂ) / 2)) := by
  have hAm : Measurable fun z : Fin d → ℝ => A *ᵥ z :=
    (Matrix.mulVecLin A).continuous_of_finiteDimensional.measurable
  have hphi : Measurable fun z : Fin d → ℝ =>
      (WithLp.toLp 2 (A *ᵥ z) : EuclideanSpace ℝ (Fin d)) :=
    (WithLp.measurable_toLp 2 (Fin d → ℝ)).comp hAm
  have hcont : ∀ s : EuclideanSpace ℝ (Fin d),
      Continuous fun y : EuclideanSpace ℝ (Fin d) =>
        Complex.exp ((inner ℝ y s : ℝ) * Complex.I) := by
    intro s
    have h1 : Continuous fun y : EuclideanSpace ℝ (Fin d) => (inner ℝ y s : ℝ) := by
      fun_prop
    exact ((Complex.continuous_ofReal.comp h1).mul continuous_const).cexp
  have hip : ∀ a b : Fin d → ℝ,
      (inner ℝ (WithLp.toLp 2 a : EuclideanSpace ℝ (Fin d))
        (WithLp.toLp 2 b : EuclideanSpace ℝ (Fin d)) : ℝ) = b ⬝ᵥ a := by
    intro a b
    simp [PiLp.inner_apply, dotProduct]
  have hwz : ∀ z : Fin d → ℝ,
      (inner ℝ (WithLp.toLp 2 (A *ᵥ z) : EuclideanSpace ℝ (Fin d))
        (WithLp.toLp 2 u : EuclideanSpace ℝ (Fin d)) : ℝ) = (Aᵀ *ᵥ u) ⬝ᵥ z := by
    intro z
    rw [hip, Matrix.dotProduct_mulVec, Matrix.mulVec_transpose]
  have hL : charFun ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
        (fun z => (WithLp.toLp 2 (A *ᵥ z) : EuclideanSpace ℝ (Fin d))))
        (WithLp.toLp 2 u : EuclideanSpace ℝ (Fin d))
      = ∫ z, Complex.exp ((((Aᵀ *ᵥ u) ⬝ᵥ z : ℝ) : ℂ) * Complex.I)
          ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1) := by
    rw [charFun_apply,
      integral_map hphi.aemeasurable (hcont (WithLp.toLp 2 u)).aestronglyMeasurable]
    exact integral_congr_ae (Filter.Eventually.of_forall fun z => by simp only [hwz])
  have hR : ∫ z, Complex.exp ((((Aᵀ *ᵥ u) ⬝ᵥ z : ℝ) : ℂ) * Complex.I)
        ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = charFun ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
          (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d)))
          (WithLp.toLp 2 (Aᵀ *ᵥ u) : EuclideanSpace ℝ (Fin d)) := by
    rw [charFun_apply, ← MeasurableEquiv.coe_toLp, integral_map_equiv]
    refine integral_congr_ae (Filter.Eventually.of_forall fun z => ?_)
    simp only [MeasurableEquiv.toLp_apply]
    rw [hip]
  have hfac : ∀ i : Fin d, charFun (gaussianReal 0 1) ((Aᵀ *ᵥ u) i)
      = Complex.exp (-(((Aᵀ *ᵥ u) i : ℂ) ^ 2 / 2)) := by
    intro i
    rw [charFun_gaussianReal]
    congr 1
    push_cast
    ring
  have hww : (Aᵀ *ᵥ u) ⬝ᵥ (Aᵀ *ᵥ u) = u ⬝ᵥ (A * Aᵀ) *ᵥ u := by
    rw [dotProduct_mulVec_self, Matrix.transpose_transpose]
  have hcast : ((u ⬝ᵥ (A * Aᵀ) *ᵥ u : ℝ) : ℂ) = ∑ i : Fin d, (((Aᵀ *ᵥ u) i : ℝ) : ℂ) ^ 2 := by
    rw [← hww, dotProduct]
    push_cast
    exact Finset.sum_congr rfl fun i _ => (sq _).symm
  rw [hL, hR, charFun_pi]
  simp only [hfac]
  rw [← Complex.exp_sum, hcast]
  congr 1
  simp only [Finset.sum_neg_distrib, ← Finset.sum_div]

/-- U3a: the characteristic function of the image of the standard Gaussian under any matrix
`A`, read on `EuclideanSpace`: `exp(-tᵀ A Aᵀ t / 2)`. -/
theorem charFun_map_mulVec_pi_gaussianReal {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ)
    (t : EuclideanSpace ℝ (Fin d)) :
    charFun ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
        (fun z => (WithLp.toLp 2 (A *ᵥ z) : EuclideanSpace ℝ (Fin d)))) t
      = Complex.exp (-(((WithLp.ofLp t ⬝ᵥ (A * Aᵀ) *ᵥ WithLp.ofLp t : ℝ) : ℂ) / 2)) :=
  charFun_map_mulVec_pi_gaussianReal' A (WithLp.ofLp t)

end StackedSVD
