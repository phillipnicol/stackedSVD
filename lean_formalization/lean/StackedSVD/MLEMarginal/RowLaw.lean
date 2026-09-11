/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLEMarginal.Defs

/-!
# The law of one row under the random-effects model (L5, units U3b and U4)

`reRowLaw n d θ v` is `N(0, Σ)` with `Σ = mleCov d (θ² / (n / d)) v`, the paper's
`Σ_i(v) = (1/d) I + (θ_i² / n_i) v vᵀ` (`main_paper.tex:1652`). The identification goes
through characteristic functions on `EuclideanSpace ℝ (Fin d)`:

1. `charFun_reRowLaw`: `exp(-(‖t‖²/d + θ² ⟨t, v⟩²/n) / 2)` by `integral_prod_mul`,
   `charFun_gaussianReal`, `charFun_pi`; and `tᵀ Σ t` is that quadratic form.
2. `charFun_map_mulVec_pi_gaussianReal` (`Prob/GaussianDensity.lean`) gives the same value
   for the image of the standard Gaussian under `A = sqrtMleCov` since `A Aᵀ = Σ`.
3. `Measure.ext_of_charFun`, then the measurable equivalence `toLp 2` is cancelled.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped ENNReal Matrix NNReal

namespace StackedSVD

/-! ### Two algebraic steps -/

/-- `(v vᵀ) w = (v ⬝ᵥ w) • v`, the real form of Mathlib's `Matrix.vecMulVec_mulVec`
(which carries the scalar in `MulOpposite`). -/
theorem vecMulVec_mulVec_real {d : ℕ} (v w : Fin d → ℝ) :
    Matrix.vecMulVec v v *ᵥ w = (v ⬝ᵥ w) • v := by
  funext i
  simp [Matrix.mulVec, dotProduct, Matrix.vecMulVec, Finset.mul_sum, mul_comm, mul_left_comm]

/-- The quadratic form of `mleCov`: `wᵀ Σ w = (w ⬝ᵥ w) / d + (a / d) (v ⬝ᵥ w)²`. -/
theorem dotProduct_mleCov_mulVec {d : ℕ} (a : ℝ) (v w : Fin d → ℝ) :
    w ⬝ᵥ mleCov d a v *ᵥ w = (w ⬝ᵥ w) / d + a / d * (v ⬝ᵥ w) ^ 2 := by
  rw [mleCov]
  simp only [Matrix.smul_mulVec, Matrix.add_mulVec, Matrix.one_mulVec, vecMulVec_mulVec_real,
    dotProduct_add, dotProduct_smul, smul_eq_mul]
  rw [dotProduct_comm w v]
  ring

/-- The characteristic-function integral of the standard Gaussian product measure, at a plain
vector: `∫ exp(i ⟨s, z⟩) dN(0, I) = exp(-(s ⬝ᵥ s) / 2)`. -/
theorem integral_cexp_dotProduct_pi_gaussianReal {d : ℕ} (s : Fin d → ℝ) :
    ∫ z, Complex.exp (((s ⬝ᵥ z : ℝ) : ℂ) * Complex.I)
        ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = Complex.exp (-(((s ⬝ᵥ s : ℝ) : ℂ) / 2)) := by
  have hip : ∀ a b : Fin d → ℝ,
      (inner ℝ (WithLp.toLp 2 a : EuclideanSpace ℝ (Fin d))
        (WithLp.toLp 2 b : EuclideanSpace ℝ (Fin d)) : ℝ) = b ⬝ᵥ a := by
    intro a b
    simp [PiLp.inner_apply, dotProduct]
  have hone : (fun z : Fin d → ℝ =>
        (WithLp.toLp 2 ((1 : Matrix (Fin d) (Fin d) ℝ) *ᵥ z) : EuclideanSpace ℝ (Fin d)))
      = (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d)) := by
    funext z
    rw [Matrix.one_mulVec]
  have hkey := charFun_map_mulVec_pi_gaussianReal (1 : Matrix (Fin d) (Fin d) ℝ)
      (WithLp.toLp 2 s : EuclideanSpace ℝ (Fin d))
  rw [hone] at hkey
  have hR : ∫ z, Complex.exp (((s ⬝ᵥ z : ℝ) : ℂ) * Complex.I)
        ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = charFun ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
          (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d)))
          (WithLp.toLp 2 s : EuclideanSpace ℝ (Fin d)) := by
    rw [charFun_apply, ← MeasurableEquiv.coe_toLp, integral_map_equiv]
    refine integral_congr_ae (Filter.Eventually.of_forall fun z => ?_)
    simp only [MeasurableEquiv.toLp_apply]
    rw [hip]
  rw [hR, hkey, Matrix.transpose_one, Matrix.mul_one, Matrix.one_mulVec]

/-! ### U3b: the characteristic function of one row -/

/-- U3b at a plain frequency vector `w`. -/
theorem charFun_reRowLaw' {n d : ℕ} (hn : 0 < n) (hd : 0 < d) (θ : ℝ) (v w : Fin d → ℝ) :
    charFun ((reRowLaw n d θ v).map (WithLp.toLp 2))
        (WithLp.toLp 2 w : EuclideanSpace ℝ (Fin d))
      = Complex.exp (-(((w ⬝ᵥ mleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ w : ℝ) : ℂ) / 2)) := by
  have hdR : (0 : ℝ) < d := Nat.cast_pos.mpr hd
  have hnR : (0 : ℝ) < n := Nat.cast_pos.mpr hn
  have hd0 : (d : ℝ) ≠ 0 := ne_of_gt hdR
  have hn0 : (n : ℝ) ≠ 0 := ne_of_gt hnR
  have hmap : Measurable (fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2) :=
    measurable_reRowMap d θ v
  have htoLp : Measurable (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d)) :=
    WithLp.measurable_toLp 2 (Fin d → ℝ)
  have hip : ∀ a b : Fin d → ℝ,
      (inner ℝ (WithLp.toLp 2 a : EuclideanSpace ℝ (Fin d))
        (WithLp.toLp 2 b : EuclideanSpace ℝ (Fin d)) : ℝ) = b ⬝ᵥ a := by
    intro a b
    simp [PiLp.inner_apply, dotProduct]
  have hcont : Continuous fun y : EuclideanSpace ℝ (Fin d) =>
      Complex.exp ((inner ℝ y (WithLp.toLp 2 w : EuclideanSpace ℝ (Fin d)) : ℝ)
        * Complex.I) := by
    have h1 : Continuous fun y : EuclideanSpace ℝ (Fin d) =>
        (inner ℝ y (WithLp.toLp 2 w : EuclideanSpace ℝ (Fin d)) : ℝ) := by fun_prop
    exact ((Complex.continuous_ofReal.comp h1).mul continuous_const).cexp
  have hstep : charFun ((reRowLaw n d θ v).map (WithLp.toLp 2))
        (WithLp.toLp 2 w : EuclideanSpace ℝ (Fin d))
      = ∫ p, Complex.exp (((θ * (v ⬝ᵥ w) * p.1 : ℝ) : ℂ) * Complex.I)
            * Complex.exp (((((Real.sqrt d)⁻¹ • w) ⬝ᵥ p.2 : ℝ) : ℂ) * Complex.I)
          ∂((gaussianReal 0 (n : ℝ≥0)⁻¹).prod
              (Measure.pi fun _ : Fin d => gaussianReal 0 1)) := by
    rw [reRowLaw, Measure.map_map htoLp hmap, charFun_apply,
      integral_map (htoLp.comp hmap).aemeasurable hcont.aestronglyMeasurable]
    refine integral_congr_ae (Filter.Eventually.of_forall fun p => ?_)
    have hre : w ⬝ᵥ ((θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2)
        = θ * (v ⬝ᵥ w) * p.1 + ((Real.sqrt d)⁻¹ • w) ⬝ᵥ p.2 := by
      simp only [dotProduct_add, dotProduct_smul, smul_dotProduct, smul_eq_mul]
      rw [dotProduct_comm w v]
      ring
    simp only [Function.comp_apply, hip, hre]
    push_cast
    rw [add_mul, Complex.exp_add]
  have hgauss : ∫ u, Complex.exp (((θ * (v ⬝ᵥ w) * u : ℝ) : ℂ) * Complex.I)
        ∂(gaussianReal 0 (n : ℝ≥0)⁻¹)
      = Complex.exp (-((((n : ℝ)⁻¹ * (θ * (v ⬝ᵥ w)) ^ 2 : ℝ) : ℂ) / 2)) := by
    have h1 : ∫ u, Complex.exp (((θ * (v ⬝ᵥ w) * u : ℝ) : ℂ) * Complex.I)
          ∂(gaussianReal 0 (n : ℝ≥0)⁻¹)
        = charFun (gaussianReal 0 (n : ℝ≥0)⁻¹) (θ * (v ⬝ᵥ w)) := by
      rw [charFun_apply_real]
      refine integral_congr_ae (Filter.Eventually.of_forall fun u => ?_)
      push_cast
      ring
    rw [h1, charFun_gaussianReal]
    congr 1
    push_cast
    ring
  have hdd : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = 1 / (d : ℝ) := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d), one_div]
  have hss : ((Real.sqrt d)⁻¹ • w) ⬝ᵥ ((Real.sqrt d)⁻¹ • w) = (w ⬝ᵥ w) / d := by
    rw [smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc, hdd]
    ring
  have hquad : w ⬝ᵥ mleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ w
      = (n : ℝ)⁻¹ * (θ * (v ⬝ᵥ w)) ^ 2 + ((Real.sqrt d)⁻¹ • w) ⬝ᵥ ((Real.sqrt d)⁻¹ • w) := by
    rw [dotProduct_mleCov_mulVec, hss]
    field_simp
    ring
  rw [hstep, integral_prod_mul
      (fun u : ℝ => Complex.exp (((θ * (v ⬝ᵥ w) * u : ℝ) : ℂ) * Complex.I))
      (fun z : Fin d → ℝ => Complex.exp (((((Real.sqrt d)⁻¹ • w) ⬝ᵥ z : ℝ) : ℂ) * Complex.I)),
    hgauss, integral_cexp_dotProduct_pi_gaussianReal, ← Complex.exp_add, hquad]
  congr 1
  push_cast
  ring

/-- U3b: the characteristic function of a row, read on `EuclideanSpace`, is
`exp(-tᵀ Σ t / 2)` with `Σ = mleCov d (θ² / (n / d)) v`. -/
theorem charFun_reRowLaw {n d : ℕ} (hn : 0 < n) (hd : 0 < d) (θ : ℝ) (v : Fin d → ℝ)
    (t : EuclideanSpace ℝ (Fin d)) :
    charFun ((reRowLaw n d θ v).map (WithLp.toLp 2)) t
      = Complex.exp (-(((WithLp.ofLp t ⬝ᵥ mleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ WithLp.ofLp t
          : ℝ) : ℂ) / 2)) :=
  charFun_reRowLaw' hn hd θ v (WithLp.ofLp t)

/-- U4: a row is the image of the standard Gaussian under `sqrtMleCov`. -/
theorem reRowLaw_eq_map_sqrtMleCov {n d : ℕ} (hn : 0 < n) (hd : 0 < d) (θ : ℝ)
    (v : Fin d → ℝ) :
    reRowLaw n d θ v
      = (Measure.pi fun _ : Fin d => gaussianReal 0 1).map
          (fun z => sqrtMleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ z) := by
  have ha : (0 : ℝ) ≤ θ ^ 2 / ((n : ℝ) / d) := by positivity
  have hAm : Measurable fun z : Fin d → ℝ => sqrtMleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ z :=
    (Matrix.mulVecLin (sqrtMleCov d (θ ^ 2 / ((n : ℝ) / d)) v)).continuous_of_finiteDimensional
      |>.measurable
  have htoLp : Measurable (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d)) :=
    WithLp.measurable_toLp 2 (Fin d → ℝ)
  have hp0 : IsProbabilityMeasure ((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
      (fun z => sqrtMleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ z)) :=
    Measure.isProbabilityMeasure_map hAm.aemeasurable
  refine (MeasurableEquiv.toLp 2 (Fin d → ℝ)).map_measurableEquiv_injective ?_
  rw [MeasurableEquiv.coe_toLp]
  have hp1 : IsProbabilityMeasure
      ((reRowLaw n d θ v).map (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d))) :=
    Measure.isProbabilityMeasure_map htoLp.aemeasurable
  have hp2 : IsProbabilityMeasure
      (((Measure.pi fun _ : Fin d => gaussianReal 0 1).map
          (fun z => sqrtMleCov d (θ ^ 2 / ((n : ℝ) / d)) v *ᵥ z)).map
        (WithLp.toLp 2 : (Fin d → ℝ) → EuclideanSpace ℝ (Fin d))) :=
    Measure.isProbabilityMeasure_map htoLp.aemeasurable
  refine Measure.ext_of_charFun (funext fun t => ?_)
  rw [charFun_reRowLaw hn hd θ v t, Measure.map_map htoLp hAm, Function.comp_def,
    charFun_map_mulVec_pi_gaussianReal, sqrtMleCov_mul_transpose d ha v]

/-- A row has the Gaussian density with covariance `mleCov d (θ² / (n / d)) v`. -/
theorem reRowLaw_eq_withDensity {n d : ℕ} (hn : 0 < n) (hd : 0 < d) (θ : ℝ) (v : Fin d → ℝ) :
    reRowLaw n d θ v
      = (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x)) := by
  have ha : (0 : ℝ) ≤ θ ^ 2 / ((n : ℝ) / d) := by positivity
  rw [reRowLaw_eq_map_sqrtMleCov hn hd,
    map_mulVec_pi_gaussianReal _ (det_sqrtMleCov_ne_zero hd ha v),
    sqrtMleCov_mul_transpose d ha v]

end StackedSVD
