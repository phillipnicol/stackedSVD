/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.GaussianMatrix
import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinIdentity
import StatsMLlib.Probability.Gaussian.LipschitzConcentration

/-!
# Gaussian adapters: Stein's identity and Lipschitz concentration on `gaussianMatrix`

Items R1, R2 and R3 need two analytic tools on the canonical Gaussian matrix law. Both exist,
on `stdGaussian` of an inner product space (COLT83) and on `stdGaussianE n` (StatsMLlib). This
file moves them to `gaussianMatrix p d` along `matrixEquivE` of
`StackedSVD/Prob/GaussianMatrix.lean`, and adds the complex-valued form of Stein that item R1
needs (the resolvent entries are complex).

Contents.

* `integral_inner_mul_stdGaussian_complex`: Stein for a `C¹` function `E → ℂ` with bounded
  derivative, by real and imaginary parts.
* `integral_apply_mul_stdGaussianE`, `integral_apply_mul_stdGaussianE_complex`: the coordinate
  form on `EuclideanSpace ℝ (Fin n)` under StatsMLlib's `stdGaussianE`.
* `integral_entry_mul_gaussianMatrix`, `integral_entry_mul_gaussianMatrix_complex`: the entry
  form on `gaussianMatrix p d`. The function is stated on `EuclideanSpace ℝ (Fin (p * d))`
  because `Matrix` carries the sup norm, not the Frobenius norm (modeling choice 2 of
  `notes/archive/rmt_R1.md`).
* `measure_ge_le_of_lipschitz`, `measure_abs_ge_le_of_lipschitz`: the one-sided and two-sided
  Gaussian concentration bounds on `gaussianMatrix p d`.

The variance corollary `Var F ≤ 4 L²` is **not** here; see the report `proof_rinfra.md`. Items
R1 and R2 use it only through Chebyshev, and the two-sided bound
`measure_abs_ge_le_of_lipschitz` gives the same conclusion directly with a better rate.
-/

open MeasureTheory ProbabilityTheory
open scoped Matrix RealInnerProductSpace ENNReal NNReal

namespace StackedSVD

/-! ### Stein's identity, complex valued -/

section ComplexStein

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  [MeasurableSpace E] [BorelSpace E] {F : E → ℂ} {L : ℝ}

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem contDiff_re_comp (hF : ContDiff ℝ 1 F) : ContDiff ℝ 1 fun x => (F x).re :=
  Complex.reCLM.contDiff.comp hF

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem contDiff_im_comp (hF : ContDiff ℝ 1 F) : ContDiff ℝ 1 fun x => (F x).im :=
  Complex.imCLM.contDiff.comp hF

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem fderiv_re_comp (hF : ContDiff ℝ 1 F) (x : E) :
    fderiv ℝ (fun y => (F y).re) x = Complex.reCLM.comp (fderiv ℝ F x) := by
  have hd : DifferentiableAt ℝ F x := (hF.differentiable one_ne_zero) x
  exact (Complex.reCLM.hasFDerivAt.comp x hd.hasFDerivAt).fderiv

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem fderiv_im_comp (hF : ContDiff ℝ 1 F) (x : E) :
    fderiv ℝ (fun y => (F y).im) x = Complex.imCLM.comp (fderiv ℝ F x) := by
  have hd : DifferentiableAt ℝ F x := (hF.differentiable one_ne_zero) x
  exact (Complex.imCLM.hasFDerivAt.comp x hd.hasFDerivAt).fderiv

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem norm_fderiv_re_comp_le (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (x : E) :
    ‖fderiv ℝ (fun y => (F y).re) x‖ ≤ L := by
  rw [fderiv_re_comp hF x]
  calc ‖Complex.reCLM.comp (fderiv ℝ F x)‖ ≤ ‖Complex.reCLM‖ * ‖fderiv ℝ F x‖ :=
        ContinuousLinearMap.opNorm_comp_le _ _
    _ ≤ L := by rw [Complex.reCLM_norm, one_mul]; exact hL x

omit [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E] in
theorem norm_fderiv_im_comp_le (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (x : E) :
    ‖fderiv ℝ (fun y => (F y).im) x‖ ≤ L := by
  rw [fderiv_im_comp hF x]
  calc ‖Complex.imCLM.comp (fderiv ℝ F x)‖ ≤ ‖Complex.imCLM‖ * ‖fderiv ℝ F x‖ :=
        ContinuousLinearMap.opNorm_comp_le _ _
    _ ≤ L := by rw [Complex.imCLM_norm, one_mul]; exact hL x

/-- **Stein's identity, complex valued.** `E[⟪a, g⟫ F(g)] = E[DF(g) a]` for a `C¹` function
`F : E → ℂ` with bounded derivative. Real and imaginary parts of the COLT83 statement
`ProbabilityTheory.integral_inner_mul_stdGaussian`. -/
theorem integral_inner_mul_stdGaussian_complex (hF : ContDiff ℝ 1 F)
    (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (a : E) :
    ∫ x, (⟪a, x⟫ : ℂ) * F x ∂(stdGaussian E) = ∫ x, fderiv ℝ F x a ∂(stdGaussian E) := by
  have hFre := contDiff_re_comp hF
  have hFim := contDiff_im_comp hF
  have hLre := norm_fderiv_re_comp_le hF hL
  have hLim := norm_fderiv_im_comp_le hF hL
  -- integrability of the four real pieces
  have hi1 : Integrable (fun x => ⟪a, x⟫ * (F x).re) (stdGaussian E) :=
    IsGaussian.integrable_inner_mul_of_norm_fderiv_le (μ := stdGaussian E) hFre hLre a
  have hi2 : Integrable (fun x => ⟪a, x⟫ * (F x).im) (stdGaussian E) :=
    IsGaussian.integrable_inner_mul_of_norm_fderiv_le (μ := stdGaussian E) hFim hLim a
  have hj1 : Integrable (fun x => fderiv ℝ (fun y => (F y).re) x a) (stdGaussian E) :=
    IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le (μ := stdGaussian E) hFre hLre a
  have hj2 : Integrable (fun x => fderiv ℝ (fun y => (F y).im) x a) (stdGaussian E) :=
    IsGaussian.integrable_fderiv_apply_of_norm_fderiv_le (μ := stdGaussian E) hFim hLim a
  -- the two real Stein identities
  have hre := ProbabilityTheory.integral_inner_mul_stdGaussian hFre hLre a
  have him := ProbabilityTheory.integral_inner_mul_stdGaussian hFim hLim a
  -- split both sides into real and imaginary parts
  have hsplitL : ∀ x : E, (⟪a, x⟫ : ℂ) * F x
      = ((⟪a, x⟫ * (F x).re : ℝ) : ℂ) + ((⟪a, x⟫ * (F x).im : ℝ) : ℂ) * Complex.I := by
    intro x
    apply Complex.ext <;> simp
  have hsplitR : ∀ x : E, fderiv ℝ F x a
      = ((fderiv ℝ (fun y => (F y).re) x a : ℝ) : ℂ)
        + ((fderiv ℝ (fun y => (F y).im) x a : ℝ) : ℂ) * Complex.I := by
    intro x
    rw [fderiv_re_comp hF x, fderiv_im_comp hF x]
    apply Complex.ext <;> simp
  have hi1C : Integrable (fun x : E => ((⟪a, x⟫ * (F x).re : ℝ) : ℂ)) (stdGaussian E) :=
    hi1.ofReal
  have hi2C : Integrable
      (fun x : E => ((⟪a, x⟫ * (F x).im : ℝ) : ℂ) * Complex.I) (stdGaussian E) :=
    (Integrable.ofReal (𝕜 := ℂ) hi2).mul_const Complex.I
  have hj1C : Integrable
      (fun x : E => ((fderiv ℝ (fun y => (F y).re) x a : ℝ) : ℂ)) (stdGaussian E) :=
    hj1.ofReal
  have hj2C : Integrable
      (fun x : E => ((fderiv ℝ (fun y => (F y).im) x a : ℝ) : ℂ) * Complex.I)
      (stdGaussian E) :=
    (Integrable.ofReal (𝕜 := ℂ) hj2).mul_const Complex.I
  have hL' : ∫ x, (⟪a, x⟫ : ℂ) * F x ∂(stdGaussian E)
      = ((∫ x, ⟪a, x⟫ * (F x).re ∂(stdGaussian E) : ℝ) : ℂ)
        + ((∫ x, ⟪a, x⟫ * (F x).im ∂(stdGaussian E) : ℝ) : ℂ) * Complex.I := by
    simp_rw [hsplitL]
    rw [integral_add hi1C hi2C, integral_mul_const, integral_complex_ofReal,
      integral_complex_ofReal]
  have hR' : ∫ x, fderiv ℝ F x a ∂(stdGaussian E)
      = ((∫ x, fderiv ℝ (fun y => (F y).re) x a ∂(stdGaussian E) : ℝ) : ℂ)
        + ((∫ x, fderiv ℝ (fun y => (F y).im) x a ∂(stdGaussian E) : ℝ) : ℂ) * Complex.I := by
    simp_rw [hsplitR]
    rw [integral_add hj1C hj2C, integral_mul_const, integral_complex_ofReal,
      integral_complex_ofReal]
  rw [hL', hR', hre, him]

end ComplexStein

/-! ### The coordinate form on `EuclideanSpace ℝ (Fin n)` -/

section Coordinates

variable {n : ℕ} {L : ℝ}

/-- **Stein's identity** in one coordinate, under StatsMLlib's `stdGaussianE`. -/
theorem integral_apply_mul_stdGaussianE {F : EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (i : Fin n) :
    ∫ x, x i * F x ∂(GaussianMeasure.stdGaussianE n)
      = ∫ x, fderiv ℝ F x (EuclideanSpace.single i 1)
          ∂(GaussianMeasure.stdGaussianE n) := by
  rw [stdGaussianE_eq]
  have h := ProbabilityTheory.integral_inner_mul_stdGaussian hF hL
    (EuclideanSpace.single i 1 : EuclideanSpace ℝ (Fin n))
  have hinner : ∀ x : EuclideanSpace ℝ (Fin n),
      ⟪(EuclideanSpace.single i 1 : EuclideanSpace ℝ (Fin n)), x⟫ = x i := by
    intro x
    rw [EuclideanSpace.inner_single_left]
    simp
  simp_rw [hinner] at h
  exact h

/-- The complex-valued form in one coordinate. -/
theorem integral_apply_mul_stdGaussianE_complex {F : EuclideanSpace ℝ (Fin n) → ℂ}
    (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (i : Fin n) :
    ∫ x, (x i : ℂ) * F x ∂(GaussianMeasure.stdGaussianE n)
      = ∫ x, fderiv ℝ F x (EuclideanSpace.single i 1)
          ∂(GaussianMeasure.stdGaussianE n) := by
  rw [stdGaussianE_eq]
  have h := integral_inner_mul_stdGaussian_complex
    (E := EuclideanSpace ℝ (Fin n)) hF hL (EuclideanSpace.single i 1)
  have hinner : ∀ x : EuclideanSpace ℝ (Fin n),
      ⟪(EuclideanSpace.single i 1 : EuclideanSpace ℝ (Fin n)), x⟫ = x i := by
    intro x
    rw [EuclideanSpace.inner_single_left]
    simp
  simp_rw [hinner] at h
  exact h

end Coordinates

/-! ### The entry form on `gaussianMatrix` -/

section Matrix

variable {p d : ℕ} {L : ℝ}

private theorem apply_matrixEquivE (Z : Matrix (Fin p) (Fin d) ℝ) (i : Fin p) (j : Fin d) :
    matrixEquivE p d Z (finProdFinEquiv (i, j)) = Z i j := by
  simp

/-- **Stein's identity for one entry of a Gaussian matrix.** `F` is stated on the flattened
space, so `F ∘ matrixEquivE p d` is the function of the matrix. -/
theorem integral_entry_mul_gaussianMatrix {F : EuclideanSpace ℝ (Fin (p * d)) → ℝ}
    (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (i : Fin p) (j : Fin d) :
    ∫ Z, Z i j * F (matrixEquivE p d Z) ∂(gaussianMatrix p d)
      = ∫ Z, fderiv ℝ F (matrixEquivE p d Z)
          (EuclideanSpace.single (finProdFinEquiv (i, j)) 1) ∂(gaussianMatrix p d) := by
  have hmp := measurePreserving_matrixEquivE p d
  have h1 : ∫ Z, Z i j * F (matrixEquivE p d Z) ∂(gaussianMatrix p d)
      = ∫ x, x (finProdFinEquiv (i, j)) * F x ∂(GaussianMeasure.stdGaussianE (p * d)) := by
    rw [← hmp.integral_comp' fun x => x (finProdFinEquiv (i, j)) * F x]
    exact integral_congr_ae (Filter.Eventually.of_forall fun Z => by
      simp only [apply_matrixEquivE])
  have h2 : ∫ Z, fderiv ℝ F (matrixEquivE p d Z)
        (EuclideanSpace.single (finProdFinEquiv (i, j)) 1) ∂(gaussianMatrix p d)
      = ∫ x, fderiv ℝ F x (EuclideanSpace.single (finProdFinEquiv (i, j)) 1)
          ∂(GaussianMeasure.stdGaussianE (p * d)) :=
    hmp.integral_comp' fun x => fderiv ℝ F x (EuclideanSpace.single (finProdFinEquiv (i, j)) 1)
  rw [h1, h2]
  exact integral_apply_mul_stdGaussianE hF hL _

/-- The complex-valued entry form, which item R1 uses on the resolvent entries. -/
theorem integral_entry_mul_gaussianMatrix_complex {F : EuclideanSpace ℝ (Fin (p * d)) → ℂ}
    (hF : ContDiff ℝ 1 F) (hL : ∀ x, ‖fderiv ℝ F x‖ ≤ L) (i : Fin p) (j : Fin d) :
    ∫ Z, (Z i j : ℂ) * F (matrixEquivE p d Z) ∂(gaussianMatrix p d)
      = ∫ Z, fderiv ℝ F (matrixEquivE p d Z)
          (EuclideanSpace.single (finProdFinEquiv (i, j)) 1) ∂(gaussianMatrix p d) := by
  have hmp := measurePreserving_matrixEquivE p d
  have h1 : ∫ Z, (Z i j : ℂ) * F (matrixEquivE p d Z) ∂(gaussianMatrix p d)
      = ∫ x, (x (finProdFinEquiv (i, j)) : ℂ) * F x
          ∂(GaussianMeasure.stdGaussianE (p * d)) := by
    rw [← hmp.integral_comp' fun x => (x (finProdFinEquiv (i, j)) : ℂ) * F x]
    exact integral_congr_ae (Filter.Eventually.of_forall fun Z => by
      simp only [apply_matrixEquivE])
  have h2 : ∫ Z, fderiv ℝ F (matrixEquivE p d Z)
        (EuclideanSpace.single (finProdFinEquiv (i, j)) 1) ∂(gaussianMatrix p d)
      = ∫ x, fderiv ℝ F x (EuclideanSpace.single (finProdFinEquiv (i, j)) 1)
          ∂(GaussianMeasure.stdGaussianE (p * d)) :=
    hmp.integral_comp' fun x => fderiv ℝ F x (EuclideanSpace.single (finProdFinEquiv (i, j)) 1)
  rw [h1, h2]
  exact integral_apply_mul_stdGaussianE_complex hF hL _

end Matrix

/-! ### Gaussian concentration on `gaussianMatrix` -/

section Concentration

variable {p d : ℕ} {F : EuclideanSpace ℝ (Fin (p * d)) → ℝ} {LL : ℝ≥0}

theorem integral_comp_matrixEquivE (F : EuclideanSpace ℝ (Fin (p * d)) → ℝ) :
    ∫ Z, F (matrixEquivE p d Z) ∂(gaussianMatrix p d)
      = ∫ x, F x ∂(GaussianMeasure.stdGaussianE (p * d)) :=
  (measurePreserving_matrixEquivE p d).integral_comp' F

/-- **One-sided Gaussian concentration on the matrix law.** For `F` `L`-Lipschitz in the
Frobenius metric (`dist_matrixEquivE_symm_sq`), the upper tail of `F` above its mean is
sub-Gaussian with variance proxy `L²`. -/
theorem measure_ge_le_of_lipschitz (hpd : 0 < p * d) (hL : 0 < LL)
    (hF : LipschitzWith LL F) {t : ℝ} (ht : 0 < t) :
    ((gaussianMatrix p d)
        {Z | (∫ Z', F (matrixEquivE p d Z') ∂(gaussianMatrix p d)) + t
              ≤ F (matrixEquivE p d Z)}).toReal
      ≤ Real.exp (-t ^ 2 / (2 * (LL : ℝ) ^ 2)) := by
  have hmean := integral_comp_matrixEquivE (p := p) (d := d) F
  set c : ℝ := ∫ x, F x ∂(GaussianMeasure.stdGaussianE (p * d)) with hc
  have hmeas : MeasurableSet {x : EuclideanSpace ℝ (Fin (p * d)) | t ≤ F x - c} := by
    have hcont : Continuous F := hF.continuous
    exact measurableSet_le measurable_const (hcont.measurable.sub measurable_const)
  have hset : {Z : Matrix (Fin p) (Fin d) ℝ |
        (∫ Z', F (matrixEquivE p d Z') ∂(gaussianMatrix p d)) + t ≤ F (matrixEquivE p d Z)}
      = (matrixEquivE p d) ⁻¹' {x | t ≤ F x - c} := by
    ext Z
    simp only [Set.mem_preimage, Set.mem_ofPred_eq, hmean]
    constructor <;> intro h <;> linarith
  rw [hset, (measurePreserving_matrixEquivE p d).measure_preimage hmeas.nullMeasurableSet]
  exact GaussianLipConcen.gaussian_lipschitz_concentration_one_sided hpd hL hF t ht

/-- **Two-sided Gaussian concentration on the matrix law.** This is what items R1 and R2 use
in place of Chebyshev. -/
theorem measure_abs_ge_le_of_lipschitz (hpd : 0 < p * d) (hL : 0 < LL)
    (hF : LipschitzWith LL F) {t : ℝ} (ht : 0 < t) :
    ((gaussianMatrix p d)
        {Z | t ≤ |F (matrixEquivE p d Z)
              - ∫ Z', F (matrixEquivE p d Z') ∂(gaussianMatrix p d)|}).toReal
      ≤ 2 * Real.exp (-t ^ 2 / (2 * (LL : ℝ) ^ 2)) := by
  have hmean := integral_comp_matrixEquivE (p := p) (d := d) F
  set c : ℝ := ∫ x, F x ∂(GaussianMeasure.stdGaussianE (p * d)) with hc
  have hmeas : MeasurableSet {x : EuclideanSpace ℝ (Fin (p * d)) | t ≤ |F x - c|} := by
    have hcont : Continuous F := hF.continuous
    exact measurableSet_le measurable_const
      ((hcont.measurable.sub measurable_const).abs)
  have hset : {Z : Matrix (Fin p) (Fin d) ℝ |
        t ≤ |F (matrixEquivE p d Z)
              - ∫ Z', F (matrixEquivE p d Z') ∂(gaussianMatrix p d)|}
      = (matrixEquivE p d) ⁻¹' {x | t ≤ |F x - c|} := by
    ext Z
    simp only [Set.mem_preimage, Set.mem_ofPred_eq, hmean]
  rw [hset, (measurePreserving_matrixEquivE p d).measure_preimage hmeas.nullMeasurableSet]
  exact GaussianLipConcen.gaussian_lipschitz_concentration hpd hL hF t ht

end Concentration

end StackedSVD
