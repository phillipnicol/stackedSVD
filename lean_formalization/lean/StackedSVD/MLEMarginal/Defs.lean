/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLE
import StackedSVD.Prob.GaussianDensity

/-!
# The random-effects model behind `app:wstacksvd_mle` (L5): definitions

The paper (`main_paper.tex:1645`) keeps `X_i = θ_i u_i vᵀ + E_i` and makes `u_i` random,
`u_i ~ N(0, I_{n_i} / n_i)`, independent of the Gaussian noise `E_i = Z_i / √d`. The review
note is `notes/archive/L5_mle_marginal.md`.

* `reRowLaw n d θ v`: the law of one row, `θ u v + z / √d` with `u ~ N(0, 1/n)`,
  `z ~ N(0, I_d)`.
* `reTableLaw n d θ v`: the law of one table, `θ u vᵀ + Z / √d` with `u ~ N(0, I_n / n)`
  and `Z ~ gaussianMatrix n d`.
* `reJointLaw n d θ v`: the `M` tables, independent.
* `lebesgueMatrix n d`: Lebesgue measure on `n × d` matrices, in the product shape of
  `gaussianMatrix`.
* `reDensity`, `reLogLik`: the joint density `∏ᵢ ∏ₖ gaussDensity d (Σ_i(v)) (X_i)_k` with
  `Σ_i(v) = mleCov d (θ_i² / (n_i / d)) v`, and its logarithm.
* `sqrtMleCov d a v = (√d)⁻¹ (I + b v vᵀ)` with `b = a / (1 + √(1 + a ‖v‖²))`, an explicit
  square root of `mleCov d a v`: `sqrtMleCov * sqrtMleCovᵀ = mleCov`.

`MultiTableModel` is untouched: its `u` stays deterministic. The random-effects law is a
separate object, and the headline theorem (`MLEMarginal/Main.lean`) evaluates its density
at the observed matrices `(m.tbl i).X N ω`.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped ENNReal Matrix NNReal

namespace StackedSVD

/-! ### An explicit square root of `mleCov` -/

/-- The scalar `b` with `(1 + b v vᵀ)² = 1 + a v vᵀ`: `b = a / (1 + √(1 + a ‖v‖²))`.
With `s = √(1 + a ‖v‖²)`: `b ‖v‖² = s - 1`, so `2b + b² ‖v‖² = b (1 + s) = a`. -/
noncomputable def sqrtCoef {d : ℕ} (a : ℝ) (v : Fin d → ℝ) : ℝ :=
  a / (1 + Real.sqrt (1 + a * (v ⬝ᵥ v)))

/-- An explicit square root of `mleCov d a v`: `(√d)⁻¹ (1 + b v vᵀ)` with `b = sqrtCoef a v`. -/
noncomputable def sqrtMleCov (d : ℕ) (a : ℝ) (v : Fin d → ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  (Real.sqrt d)⁻¹ • (1 + sqrtCoef a v • vecMulVec v v)

/-- `0 ≤ v ⬝ᵥ v` for a real vector: needed for `sqrtCoef_spec` and `det_mleCov_pos`. -/
theorem dotSq_nonneg {d : ℕ} (v : Fin d → ℝ) : 0 ≤ v ⬝ᵥ v :=
  Finset.sum_nonneg fun i _ => mul_self_nonneg (v i)

/-- `sqrtCoef` solves `2b + b² ‖v‖² = a` when `0 ≤ a`. -/
theorem sqrtCoef_spec {d : ℕ} {a : ℝ} (ha : 0 ≤ a) (v : Fin d → ℝ) :
    2 * sqrtCoef a v + sqrtCoef a v ^ 2 * (v ⬝ᵥ v) = a := by
  have hvv : 0 ≤ v ⬝ᵥ v := dotSq_nonneg v
  have harg : (0 : ℝ) ≤ 1 + a * (v ⬝ᵥ v) := by
    have h0 : (0 : ℝ) ≤ a * (v ⬝ᵥ v) := mul_nonneg ha hvv
    linarith
  rw [sqrtCoef]
  set s := Real.sqrt (1 + a * (v ⬝ᵥ v))
  have hs2 : s ^ 2 = 1 + a * (v ⬝ᵥ v) := Real.sq_sqrt harg
  have hsnn : (0 : ℝ) ≤ s := Real.sqrt_nonneg _
  have h1s : (1 : ℝ) + s ≠ 0 := by
    have hpos : (0 : ℝ) < 1 + s := by linarith
    exact hpos.ne'
  have hkey : 2 * a * (1 + s) + a ^ 2 * (v ⬝ᵥ v) = a * (1 + s) ^ 2 := by
    linear_combination (-a) * hs2
  have hexpand : 2 * (a / (1 + s)) + (a / (1 + s)) ^ 2 * (v ⬝ᵥ v)
      = (2 * a * (1 + s) + a ^ 2 * (v ⬝ᵥ v)) / (1 + s) ^ 2 := by
    field_simp
  rw [hexpand, hkey, mul_div_assoc, div_self (pow_ne_zero 2 h1s), mul_one]

theorem sqrtMleCov_mul_transpose (d : ℕ) {a : ℝ} (ha : 0 ≤ a) (v : Fin d → ℝ) :
    sqrtMleCov d a v * (sqrtMleCov d a v)ᵀ = mleCov d a v := by
  have htr : (sqrtMleCov d a v)ᵀ = sqrtMleCov d a v := by
    rw [sqrtMleCov]
    simp
  rw [htr, sqrtMleCov, mleCov, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  have hdd : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = 1 / (d : ℝ) := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d), one_div]
  rw [hdd]
  congr 1
  simp only [Matrix.add_mul, Matrix.mul_add, Matrix.one_mul, Matrix.mul_one, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.vecMulVec_mul_vecMulVec, Matrix.vecMulVec_smul, smul_smul]
  match_scalars <;> try ring1
  linear_combination sqrtCoef_spec ha v

/-- `det_mleCov` at every `v`: `det Σ = d^{-d} (1 + a ‖v‖²)`. -/
theorem det_mleCov' (d : ℕ) (a : ℝ) (v : Fin d → ℝ) :
    (mleCov d a v).det = (1 / (d : ℝ)) ^ d * (1 + a * (v ⬝ᵥ v)) := by
  rw [mleCov, Matrix.det_smul, Fintype.card_fin]
  congr 1
  rw [← Matrix.smul_vecMulVec, Matrix.vecMulVec_eq (Fin 1),
    Matrix.det_one_add_replicateCol_mul_replicateRow, dotProduct_smul]
  simp

theorem det_mleCov_pos {d : ℕ} (hd : 0 < d) {a : ℝ} (ha : 0 ≤ a) (v : Fin d → ℝ) :
    0 < (mleCov d a v).det := by
  have hvv : 0 ≤ v ⬝ᵥ v := dotSq_nonneg v
  have hdR : (0 : ℝ) < (d : ℝ) := Nat.cast_pos.mpr hd
  have h1 : (0 : ℝ) < (1 / (d : ℝ)) ^ d := pow_pos (div_pos one_pos hdR) d
  have h2 : (0 : ℝ) < 1 + a * (v ⬝ᵥ v) := by
    have h3 : (0 : ℝ) ≤ a * (v ⬝ᵥ v) := mul_nonneg ha hvv
    linarith
  rw [det_mleCov']
  exact mul_pos h1 h2

/-- `(det sqrtMleCov)² = det mleCov ≠ 0`. -/
theorem det_sqrtMleCov_ne_zero {d : ℕ} (hd : 0 < d) {a : ℝ} (ha : 0 ≤ a) (v : Fin d → ℝ) :
    (sqrtMleCov d a v).det ≠ 0 := by
  intro h
  have hdet := congrArg Matrix.det (sqrtMleCov_mul_transpose d ha v)
  rw [Matrix.det_mul, Matrix.det_transpose, h, mul_zero] at hdet
  exact (det_mleCov_pos hd ha v).ne' hdet.symm

/-! ### The random-effects laws -/

/-- One row of table `i` under the random-effects model: `θ u v + z / √d` with
`u ~ N(0, 1/n)` and `z ~ N(0, I_d)`. -/
noncomputable def reRowLaw (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) : Measure (Fin d → ℝ) :=
  ((gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
    (fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2)

/-- Table `i` under the random-effects model: `θ u vᵀ + Z / √d` with `u ~ N(0, I_n / n)` and
`Z ~ gaussianMatrix n d`. -/
noncomputable def reTableLaw (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    Measure (Matrix (Fin n) (Fin d) ℝ) :=
  ((Measure.pi fun _ : Fin n => gaussianReal 0 (n : ℝ≥0)⁻¹).prod (gaussianMatrix n d)).map
    (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
      θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2)

theorem measurable_reRowMap (d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    Measurable (fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2) := by
  refine measurable_pi_lambda _ fun j => ?_
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  have hB : Measurable (fun p : ℝ × (Fin d → ℝ) => p.2 j) :=
    (measurable_pi_apply j).comp measurable_snd
  exact ((measurable_fst.const_mul θ).mul_const (v j)).add (hB.const_mul (Real.sqrt d)⁻¹)

theorem measurable_reTableMap (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    Measurable (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
      θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.vecMulVec_apply, smul_eq_mul]
  have hA : Measurable (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ => p.1 i) :=
    (measurable_pi_apply i).comp measurable_fst
  have hB : Measurable (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ => p.2 i j) :=
    (measurable_pi_apply j).comp ((measurable_pi_apply i).comp measurable_snd)
  exact ((hA.mul_const (v j)).const_mul θ).add (hB.const_mul (Real.sqrt d)⁻¹)

instance instIsProbabilityMeasure_reRowLaw (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    IsProbabilityMeasure (reRowLaw n d θ v) :=
  Measure.isProbabilityMeasure_map (measurable_reRowMap d θ v).aemeasurable

instance instIsProbabilityMeasure_reTableLaw (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    IsProbabilityMeasure (reTableLaw n d θ v) :=
  Measure.isProbabilityMeasure_map (measurable_reTableMap n d θ v).aemeasurable

/-- The `M` tables, independent (modeling choice 3 of the note). -/
noncomputable def reJointLaw {M : ℕ} (n : Fin M → ℕ) (d : ℕ) (θ : Fin M → ℝ) (v : Fin d → ℝ) :
    Measure ((i : Fin M) → Matrix (Fin (n i)) (Fin d) ℝ) :=
  Measure.pi fun i => reTableLaw (n i) d (θ i) v

/-- Lebesgue measure on `n × d` matrices, with the product shape of `gaussianMatrix`. -/
noncomputable def lebesgueMatrix (n d : ℕ) : Measure (Matrix (Fin n) (Fin d) ℝ) :=
  Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => (volume : Measure ℝ)

/-! ### The density and the log-likelihood -/

/-- The joint Lebesgue density of the random-effects model at the data `X`:
`∏ᵢ ∏ₖ gaussDensity d (Σ_i(v)) (X_i)_k` with `Σ_i(v) = mleCov d (θ_i² / (n_i / d)) v`, the
paper's `Σ_i(v)` at `c_i = n_i / d`. -/
noncomputable def reDensity {M : ℕ} (n : Fin M → ℕ) (d : ℕ) (θ : Fin M → ℝ) (v : Fin d → ℝ)
    (X : (i : Fin M) → Matrix (Fin (n i)) (Fin d) ℝ) : ℝ :=
  ∏ i, ∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k)

/-- The marginal log-likelihood, `log reDensity`. -/
noncomputable def reLogLik {M : ℕ} (n : Fin M → ℕ) (d : ℕ) (θ : Fin M → ℝ) (v : Fin d → ℝ)
    (X : (i : Fin M) → Matrix (Fin (n i)) (Fin d) ℝ) : ℝ :=
  Real.log (reDensity n d θ v X)

end StackedSVD
