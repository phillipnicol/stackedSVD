/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R3
import StackedSVD.RMT.Het.Split
import StackedSVD.RMT.Het.MPhet

/-!
# Item H8: the Sudakov-Fernique upper edge of the heteroscedastic Wishart block

Task H8 of `notes/archive/plan_heterolaw_A.md` (section 3.6 and the H8 row of section 4). This file
proves the `edge` field of `ResolventLimitsHet` (`RMT/Het/R4het.lean`) at `b = bSF c w`:
for every `ε > 0`,

`Tendsto (fun N => μ N {ω | lamMax (W₀' N ω) _ ≤ bSF c w + ε}) atTop (𝓝 1)`,

with `W₀' = Σ^{1/2} E⊥ E⊥ᵀ Σ^{1/2}` (`Split.lean`) and
`bSF c w = (√(wSqMax w) + √(∑ c_i w_i²))²` (`MPhet.lean`).

## Mathematics

Write `P = ∑ n_i`, `p = d - 1`, `S = Σ^{1/2} = diag(w_i I_{n_i})` (signed, `Split.lean` choice
1), `L = √(wSqMax w) = max_i |w_i|`. By `exists_block_hasLaw_het`, `W₀' = d⁻¹ (S B)(S B)ᵀ`
with `B ~ gaussianMatrix P p`, so `lamMax W₀' ≤ d⁻¹ ‖S B‖²`.

1. **Sudakov-Fernique** (R3a pattern). On a finite family of unit pairs `(x_i, y_i)` compare
   `X_i = ⟪S x_i, B y_i⟫` with `Y_i = ⟪S x_i, g⟫ + ‖S x_i‖ ⟪y_i, h⟫`, `g ~ N(0, I_P)`,
   `h ~ N(0, I_p)` independent. With `a = ‖S x_i‖`, `a' = ‖S x_j‖`, `s = ⟪S x_i, S x_j⟫`,
   `t = ⟪y_i, y_j⟫`:
   `E(X_i - X_j)² = a² + a'² - 2 s t` and `E(Y_i - Y_j)² = 2a² + 2a'² - 2 s - 2 a a' t`, so
   `E(Y_i - Y_j)² - E(X_i - X_j)² = a² + a'² - 2 a a' t - 2 s (1 - t) ≥ (a - a')² ≥ 0`,
   because `s ≤ a a'` (Cauchy-Schwarz) and `1 - t ≥ 0`. Hence
   `E sup X ≤ E sup Y ≤ E ‖S g‖ + L E ‖h‖ ≤ √tr(S Sᵀ) + L √p`, by
   `⟪S x, g⟫ = ⟪x, S g⟫ ≤ ‖S g‖` (`Sᵀ = S`), `‖S x‖ ≤ L`, and
   `integral_norm_le_sqrt_trace_multivariateGaussian` for `S g ~ N(0, S Sᵀ)`.
2. **Net** (R3 `integral_opNorm_le_aux`): `E ‖S B‖ ≤ √tr(S Sᵀ) + L √p`.
3. **Concentration**: `B ↦ ‖S B‖` is `(L + 1)`-Lipschitz in the entries, so
   `μ {‖S B‖ ≥ E ‖S B‖ + κ √d} ≤ exp(-κ² d / (2 (L + 1)²))` (`measure_ge_le_of_lipschitz`).
4. **Assembly** (R3c pattern): `tr(S Sᵀ) = ∑ n_i w_i²`, so `d^{-1/2} √tr(S Sᵀ) → √(∑ c_i w_i²)`
   and `√p ≤ √d`; on the good event `lamMax W₀' ≤ (r + L + κ)² ≤ bSF + ε`.

Only `w_i²` enters (`wSqMax`, `tr(S Sᵀ)`), so the sign of `w_i` plays no role. The constant
`L + 1` in step 3 avoids a hypothesis `∃ i, w i ≠ 0`; the rate is not used.

Numeric check (a session script, `check_edge.py`, not kept; seed 20260830, `d = 400`,
`M = 3`,
`c = (0.5, 1, 2)`, `w = (1, -0.6, 0.3)`): `lamMax W₀'` mean `3.508`, max `3.613` over 20
draws, `bHet = 3.559`, `bSF = 4.080`; the increment gap of step 1 is nonnegative on 2000
random pairs (minimum `-4e-16`).

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/R3het.lean` exit 0, 0 `sorry`, no warning.
See `notes/archive/agent_reports/h8_edge.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace R3het

variable {P p : ℕ}

/-! ### Deterministic step 1: the operator norm of `S * B` -/

/-- A diagonal matrix with entries bounded by `L` has operator norm at most `L`. -/
theorem opNorm_diagonal_le (a : Fin P → ℝ) {L : ℝ} (hL : 0 ≤ L) (ha : ∀ r, |a r| ≤ L) :
    ‖Matrix.diagonal a‖ ≤ L := by
  rw [Matrix.l2_opNorm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ hL fun x => ?_
  have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap)
      (Matrix.diagonal a)) x = WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) := rfl
  rw [happ]
  have hsq : ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) :
      EuclideanSpace ℝ (Fin P))‖ ^ 2 ≤ (L * ‖x‖) ^ 2 := by
    rw [mul_pow, EuclideanSpace.real_norm_sq_eq, EuclideanSpace.real_norm_sq_eq,
      Finset.mul_sum]
    refine Finset.sum_le_sum fun r _ => ?_
    have hcoord : (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) :
        EuclideanSpace ℝ (Fin P)) r = a r * x r := by
      change (Matrix.diagonal a *ᵥ WithLp.ofLp x) r = _
      rw [Matrix.mulVec_diagonal]
    rw [hcoord, mul_pow]
    have h1 : a r ^ 2 ≤ L ^ 2 := by
      rw [← sq_abs (a r)]
      exact pow_le_pow_left₀ (abs_nonneg _) (ha r) 2
    exact mul_le_mul_of_nonneg_right h1 (sq_nonneg _)
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (mul_nonneg hL (norm_nonneg _)) two_ne_zero).1 hsq

/-- `‖S * B‖ ≤ L ‖B‖` when `‖S‖ ≤ L`. -/
theorem opNorm_mul_le_of_le (S : Matrix (Fin P) (Fin P) ℝ) {L : ℝ} (hS : ‖S‖ ≤ L)
    (B : Matrix (Fin P) (Fin p) ℝ) : ‖S * B‖ ≤ L * ‖B‖ :=
  (Matrix.l2_opNorm_mul S B).trans (mul_le_mul_of_nonneg_right hS (norm_nonneg _))

/-- `B ↦ ‖S * B‖` is `L`-Lipschitz in the entries when `‖S‖ ≤ L` (the R3 adapter
`lipschitzWith_opNorm` with the fixed left factor `S`). -/
theorem lipschitzWith_opNorm_mul (S : Matrix (Fin P) (Fin P) ℝ) {L : ℝ≥0} (hS : ‖S‖ ≤ L) :
    LipschitzWith L
      (fun w : EuclideanSpace ℝ (Fin (P * p)) => ‖S * (matrixEquivE P p).symm w‖) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  have h1 : dist ‖S * (matrixEquivE P p).symm x‖ ‖S * (matrixEquivE P p).symm y‖
      ≤ ‖S * ((matrixEquivE P p).symm x - (matrixEquivE P p).symm y)‖ := by
    rw [Real.dist_eq, Matrix.mul_sub]
    exact abs_norm_sub_norm_le _ _
  have h2 : ‖(matrixEquivE P p).symm x - (matrixEquivE P p).symm y‖ ≤ dist x y := by
    rw [← matrixEquivE_symm_sub]
    calc ‖(matrixEquivE P p).symm (x - y)‖
        ≤ Real.sqrt (∑ i, ∑ j, ((matrixEquivE P p).symm (x - y)) i j ^ 2) :=
          R3.l2_opNorm_le_frobenius _
      _ = Real.sqrt (dist x y ^ 2) := by
          have hsub : ∀ i j, (matrixEquivE P p).symm (x - y) i j
              = (matrixEquivE P p).symm x i j - (matrixEquivE P p).symm y i j := by
            intro i j; simp
          rw [dist_matrixEquivE_symm_sq]
          simp only [hsub]
      _ = dist x y := Real.sqrt_sq dist_nonneg
  calc dist ‖S * (matrixEquivE P p).symm x‖ ‖S * (matrixEquivE P p).symm y‖
      ≤ ‖S * ((matrixEquivE P p).symm x - (matrixEquivE P p).symm y)‖ := h1
    _ ≤ L * ‖(matrixEquivE P p).symm x - (matrixEquivE P p).symm y‖ :=
        opNorm_mul_le_of_le S hS _
    _ ≤ L * dist x y := mul_le_mul_of_nonneg_left h2 (NNReal.coe_nonneg L)

/-- Measurability of `B ↦ ‖S * B‖`. -/
theorem measurable_opNorm_mul (S : Matrix (Fin P) (Fin P) ℝ) :
    Measurable (fun B : Matrix (Fin P) (Fin p) ℝ => ‖S * B‖) := by
  have h1 : Measurable (fun w : EuclideanSpace ℝ (Fin (P * p)) =>
      ‖S * (matrixEquivE P p).symm w‖) :=
    (lipschitzWith_opNorm_mul S (L := ⟨‖S‖, norm_nonneg _⟩) le_rfl).continuous.measurable
  have h2 := h1.comp (matrixEquivE P p).measurable
  simpa [Function.comp_def] using h2

/-- **Rayleigh bound** for the Gram matrix on the row side, `W₀' = r • (F Fᵀ)`. -/
theorem lamMax_gramT_le_opNorm_sq {r : ℝ} (hr : 0 ≤ r) (F : Matrix (Fin P) (Fin p) ℝ)
    (hs : (r • (F * Fᵀ)).IsHermitian) : lamMax (r • (F * Fᵀ)) hs ≤ r * ‖F‖ ^ 2 := by
  have heq : r • (F * Fᵀ) = r • (Fᵀᵀ * Fᵀ) := by rw [Matrix.transpose_transpose]
  have hs' : (r • (Fᵀᵀ * Fᵀ)).IsHermitian := heq ▸ hs
  rw [lamMax_congr heq hs hs']
  have h := R3.lamMax_gram_le_opNorm_sq hr Fᵀ hs'
  have hT : ‖Fᵀ‖ = ‖F‖ := by
    rw [← Matrix.conjTranspose_eq_transpose_of_trivial]
    exact Matrix.l2_opNorm_conjTranspose F
  rwa [hT] at h

/-! ### R3a-het: Sudakov-Fernique with the comparison process `⟪S x, g⟫ + ‖S x‖ ⟪y, h⟫` -/

section Gordon

variable {ι : Type*}

/-- The matrix of the comparison process `Y_i = ⟪x_i, g⟫ + ‖x_i‖ ⟪y_i, h⟫`, for `x_i` of any
length (it is `S x̃_i` below) and `y_i` unit. -/
noncomputable def cmpMatHet (x : ι → EuclideanSpace ℝ (Fin P)) (y : ι → EuclideanSpace ℝ (Fin p)) :
    Matrix ι (Fin P ⊕ Fin p) ℝ :=
  Matrix.of fun i k =>
    Sum.elim (fun a => WithLp.ofLp (x i) a) (fun b => ‖x i‖ * WithLp.ofLp (y i) b) k

section GramSection

variable (x : ι → EuclideanSpace ℝ (Fin P)) (y : ι → EuclideanSpace ℝ (Fin p))

/-- `Cov(Y_i, Y_j) = ⟪x_i, x_j⟫ + ‖x_i‖ ‖x_j‖ ⟪y_i, y_j⟫`. -/
theorem cmpMatHet_gram (i j : ι) :
    (cmpMatHet x y * (cmpMatHet x y)ᵀ) i j
      = (∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x j) a)
        + ‖x i‖ * ‖x j‖ * ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y j) b := by
  rw [Matrix.mul_apply, Fintype.sum_sum_type, Finset.mul_sum]
  simp only [cmpMatHet, Matrix.of_apply, Matrix.transpose_apply, Sum.elim_inl, Sum.elim_inr]
  congr 1
  exact Finset.sum_congr rfl fun b _ => by ring

/-- The comparison process, split into its two blocks. -/
theorem cmpMatHet_apply (w : EuclideanSpace ℝ (Fin P ⊕ Fin p)) (i : ι) :
    WithLp.ofLp (Matrix.toEuclideanLin (cmpMatHet x y) w) i
      = (∑ a, WithLp.ofLp (x i) a * WithLp.ofLp w (Sum.inl a))
        + ‖x i‖ * ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b) := by
  change (cmpMatHet x y *ᵥ WithLp.ofLp w) i = _
  simp only [Matrix.mulVec, dotProduct, cmpMatHet, Matrix.of_apply]
  rw [Fintype.sum_sum_type, Finset.mul_sum]
  simp only [Sum.elim_inl, Sum.elim_inr]
  congr 1
  exact Finset.sum_congr rfl fun b _ => by ring

variable [Finite ι]

theorem posSemidef_cmpMatHet_gram : (cmpMatHet x y * (cmpMatHet x y)ᵀ).PosSemidef := by
  simpa using Matrix.posSemidef_self_mul_conjTranspose (cmpMatHet x y)

end GramSection

/-- `‖S v‖ ≤ ‖S‖ ‖v‖` on Euclidean vectors, in the `WithLp.toLp` form used here. -/
theorem norm_toLp_mulVec_le (S : Matrix (Fin P) (Fin P) ℝ) (v : EuclideanSpace ℝ (Fin P)) :
    ‖(WithLp.toLp 2 (S *ᵥ WithLp.ofLp v) : EuclideanSpace ℝ (Fin P))‖ ≤ ‖S‖ * ‖v‖ :=
  ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) S).le_opNorm v

/-- `E ‖S g‖ ≤ √tr(S Sᵀ)` for the first block `g` of a standard Gaussian vector of
`ℝ^{P} ⊕ ℝ^{p}`: `S g ~ N(0, S Sᵀ)` and `integral_norm_le_sqrt_trace_multivariateGaussian`. -/
theorem integral_norm_mul_projL_le (S : Matrix (Fin P) (Fin P) ℝ) :
    ∫ w, ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖
        ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) ≤ Real.sqrt (S * Sᵀ).trace := by
  have hmap : (stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (S * R3.projL P p)))
      = multivariateGaussian 0 (S * Sᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin P ⊕ Fin p),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one,
      Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (R3.projL P p),
      R3.projL_mul_transpose, Matrix.one_mul]
  have h : ∫ z, ‖z‖ ∂((stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))).map
        (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (S * R3.projL P p))))
      = ∫ w, ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖
          ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) :=
    integral_map (by fun_prop) (by fun_prop)
  rw [← h, hmap]
  have hpsd : (S * Sᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose S
  exact integral_norm_le_sqrt_trace_multivariateGaussian hpsd

/-- **R3a-het**, finite family form: Sudakov-Fernique for the process
`X_i = ⟪S x_i, B y_i⟫` on unit pairs, against `Y_i = ⟪S x_i, g⟫ + ‖S x_i‖ ⟪y_i, h⟫`. The
increment inequality is checked in the header (item 1). -/
theorem integral_iSup_bilin_le_het [Finite ι] [Nonempty ι]
    (S : Matrix (Fin P) (Fin P) ℝ) (hS : Sᵀ = S) {L : ℝ} (hSL : ‖S‖ ≤ L)
    (x : ι → EuclideanSpace ℝ (Fin P)) (y : ι → EuclideanSpace ℝ (Fin p))
    (hx : ∀ i, ‖x i‖ = 1) (hy : ∀ i, ‖y i‖ = 1) :
    ∫ B, ⨆ i, (S *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i)) ∂(gaussianMatrix P p)
      ≤ Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p := by
  classical
  have := Fintype.ofFinite ι
  set xs : ι → EuclideanSpace ℝ (Fin P) :=
    fun i => WithLp.toLp 2 (S *ᵥ WithLp.ofLp (x i)) with hxs
  have hxs_apply : ∀ i, WithLp.ofLp (xs i) = S *ᵥ WithLp.ofLp (x i) := fun i => rfl
  have hxsL : ∀ i, ‖xs i‖ ≤ L := fun i => by
    calc ‖xs i‖ ≤ ‖S‖ * ‖x i‖ := norm_toLp_mulVec_le S (x i)
      _ = ‖S‖ := by rw [hx i, mul_one]
      _ ≤ L := hSL
  have hL0 : (0 : ℝ) ≤ L := (norm_nonneg _).trans hSL
  -- the increment inequality
  have hincr : ∀ i j,
      (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) i i + (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) j j
        - 2 * (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) i j
      ≤ (cmpMatHet xs y * (cmpMatHet xs y)ᵀ) i i + (cmpMatHet xs y * (cmpMatHet xs y)ᵀ) j j
        - 2 * (cmpMatHet xs y * (cmpMatHet xs y)ᵀ) i j := by
    intro i j
    simp only [R3.bilMat_gram, cmpMatHet_gram]
    rw [R3.sum_mul_self_eq_norm_sq (xs i), R3.sum_mul_self_eq_norm_sq (xs j),
      R3.sum_mul_self_eq_norm_sq (y i), R3.sum_mul_self_eq_norm_sq (y j),
      R3.sum_mul_eq_inner (xs i) (xs j), R3.sum_mul_eq_inner (y i) (y j), hy i, hy j]
    have hs : inner ℝ (xs i) (xs j) ≤ ‖xs i‖ * ‖xs j‖ := real_inner_le_norm _ _
    have ht : |inner ℝ (y i) (y j)| ≤ 1 := by
      calc |inner ℝ (y i) (y j)| ≤ ‖y i‖ * ‖y j‖ := abs_real_inner_le_norm _ _
        _ = 1 := by rw [hy i, hy j]; norm_num
    obtain ⟨-, ht2⟩ := abs_le.mp ht
    have h1 : (0 : ℝ) ≤ 1 - inner ℝ (y i) (y j) := by linarith
    nlinarith [mul_le_mul_of_nonneg_right hs h1, sq_nonneg (‖xs i‖ - ‖xs j‖)]
  -- the left side, as an integral against `N(0, S)`
  have hSmap : (stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (R3.bilMat xs y)))
      = multivariateGaussian 0 (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin P × Fin p),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hleft : ∫ B, ⨆ i, (S *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i))
        ∂(gaussianMatrix P p)
      = ∫ z, ⨆ i, WithLp.ofLp z i
          ∂(multivariateGaussian 0 (R3.bilMat xs y * (R3.bilMat xs y)ᵀ)) := by
    have h1 : ∫ B, ⨆ i, (S *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i))
          ∂(gaussianMatrix P p)
        = ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3.bilMat xs y) w) i
            ∂(stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))) := by
      rw [← (R3.measurePreserving_matrixEquivP P p).integral_comp'
        (fun w => ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3.bilMat xs y) w) i)]
      exact integral_congr_ae (Eventually.of_forall fun B =>
        iSup_congr fun i => (R3.bilMat_apply xs y B i).symm)
    rw [h1, ← hSmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    rfl
  -- the right side
  have hTmap : (stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (cmpMatHet xs y)))
      = multivariateGaussian 0 (cmpMatHet xs y * (cmpMatHet xs y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin P ⊕ Fin p),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hright : ∫ z, ⨆ i, WithLp.ofLp z i
        ∂(multivariateGaussian 0 (cmpMatHet xs y * (cmpMatHet xs y)ᵀ))
      ≤ Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p := by
    rw [← hTmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    have hpt : ∀ w : EuclideanSpace ℝ (Fin P ⊕ Fin p),
        (⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (cmpMatHet xs y) w) i)
          ≤ ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖
            + L * ‖Matrix.toEuclideanLin (R3.projR P p) w‖ := by
      intro w
      refine ciSup_le fun i => ?_
      rw [cmpMatHet_apply]
      have h1 : ∑ a, WithLp.ofLp (xs i) a * WithLp.ofLp w (Sum.inl a)
          ≤ ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖ := by
        have heq : ∑ a, WithLp.ofLp (xs i) a * WithLp.ofLp w (Sum.inl a)
            = inner ℝ (x i) (Matrix.toEuclideanLin (S * R3.projL P p) w) := by
          rw [← R3.sum_mul_eq_inner]
          have hcoord : ∀ a, WithLp.ofLp (Matrix.toEuclideanLin (S * R3.projL P p) w) a
              = (S *ᵥ WithLp.ofLp (Matrix.toEuclideanLin (R3.projL P p) w)) a := by
            intro a
            change ((S * R3.projL P p) *ᵥ WithLp.ofLp w) a
              = (S *ᵥ (R3.projL P p *ᵥ WithLp.ofLp w)) a
            rw [Matrix.mulVec_mulVec]
          simp only [hcoord, hxs_apply]
          have hdot : ∀ u v : Fin P → ℝ, ∑ a, (S *ᵥ u) a * v a = ∑ a, u a * (S *ᵥ v) a := by
            intro u v
            have h : (S *ᵥ u) ⬝ᵥ v = u ⬝ᵥ (S *ᵥ v) := by
              rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hS]
            simpa [dotProduct] using h
          have hvec : (fun a => WithLp.ofLp w (Sum.inl a))
              = WithLp.ofLp (Matrix.toEuclideanLin (R3.projL P p) w) := by
            funext a; exact (R3.projL_apply w a).symm
          rw [hdot, hvec]
        rw [heq]
        calc inner ℝ (x i) (Matrix.toEuclideanLin (S * R3.projL P p) w)
            ≤ ‖x i‖ * ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖ := real_inner_le_norm _ _
          _ = ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖ := by rw [hx i, one_mul]
      have h2 : ‖xs i‖ * ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b)
          ≤ L * ‖Matrix.toEuclideanLin (R3.projR P p) w‖ := by
        have heq : ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b)
            = inner ℝ (y i) (Matrix.toEuclideanLin (R3.projR P p) w) := by
          rw [← R3.sum_mul_eq_inner]
          exact Finset.sum_congr rfl fun b _ => by rw [R3.projR_apply]
        have hin : inner ℝ (y i) (Matrix.toEuclideanLin (R3.projR P p) w)
            ≤ ‖Matrix.toEuclideanLin (R3.projR P p) w‖ := by
          calc inner ℝ (y i) (Matrix.toEuclideanLin (R3.projR P p) w)
              ≤ ‖y i‖ * ‖Matrix.toEuclideanLin (R3.projR P p) w‖ := real_inner_le_norm _ _
            _ = ‖Matrix.toEuclideanLin (R3.projR P p) w‖ := by rw [hy i, one_mul]
        rw [heq]
        calc ‖xs i‖ * inner ℝ (y i) (Matrix.toEuclideanLin (R3.projR P p) w)
            ≤ ‖xs i‖ * ‖Matrix.toEuclideanLin (R3.projR P p) w‖ :=
              mul_le_mul_of_nonneg_left hin (norm_nonneg _)
          _ ≤ L * ‖Matrix.toEuclideanLin (R3.projR P p) w‖ :=
              mul_le_mul_of_nonneg_right (hxsL i) (norm_nonneg _)
      linarith
    calc ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (cmpMatHet xs y) w) i
          ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p)))
        ≤ ∫ w, (‖Matrix.toEuclideanLin (S * R3.projL P p) w‖
            + L * ‖Matrix.toEuclideanLin (R3.projR P p) w‖)
            ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) := by
          refine integral_mono
            (R3.integrable_ciSup_clm
              (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (cmpMatHet xs y))) _)
            ((R3.integrable_norm_toEuclideanLin (S * R3.projL P p)).add
              ((R3.integrable_norm_toEuclideanLin (R3.projR P p)).const_mul L)) hpt
      _ = (∫ w, ‖Matrix.toEuclideanLin (S * R3.projL P p) w‖
              ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))))
            + L * ∫ w, ‖Matrix.toEuclideanLin (R3.projR P p) w‖
              ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) := by
          rw [integral_add (R3.integrable_norm_toEuclideanLin (S * R3.projL P p))
            ((R3.integrable_norm_toEuclideanLin (R3.projR P p)).const_mul L),
            integral_const_mul]
      _ ≤ Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p :=
          add_le_add (integral_norm_mul_projL_le S)
            (mul_le_mul_of_nonneg_left (R3.integral_norm_projR_le P p) hL0)
  rw [hleft]
  exact (sudakov_fernique (R3.posSemidef_bilMat_gram xs y) (posSemidef_cmpMatHet_gram xs y)
    hincr).trans hright

end Gordon

/-! ### The net argument (R3 `integral_opNorm_le_aux`, with `Z := S * B`) -/

/-- The `ε₀`-step of Gordon's bound for `S * B`. -/
theorem integral_opNorm_mul_le_aux (hP : 0 < P) (hp : 0 < p) (S : Matrix (Fin P) (Fin P) ℝ)
    (hS : Sᵀ = S) {L : ℝ} (hSL : ‖S‖ ≤ L) {e : ℝ} (he0 : 0 < e) (he : 2 * e < 1) :
    ∫ B, ‖S * B‖ ∂(gaussianMatrix P p)
      ≤ (Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p) / (1 - 2 * e) := by
  classical
  obtain ⟨Nx, hNx, hsubx⟩ := R3.exists_sphere_enet P hP he0
  obtain ⟨Ny, hNy, hsuby⟩ := R3.exists_sphere_enet p hp he0
  set Nnet : RMT.CenteredMatrixBilinearNet P p e :=
    { domainNet := Ny
      codomainNet := Nx
      domain_isNet := hNy
      codomain_isNet := hNx
      domain_subset_unitBall := hsuby.trans Metric.sphere_subset_closedBall
      codomain_subset_unitBall := hsubx.trans Metric.sphere_subset_closedBall } with hNnet
  obtain ⟨vx, hvx⟩ : Nx.Nonempty := Nnet.toMatrixBilinearNet.codomainNet_nonempty_of_pos hP
  obtain ⟨vy, hvy⟩ : Ny.Nonempty := Nnet.toMatrixBilinearNet.domainNet_nonempty_of_pos hp
  have : Nonempty (↥Nx × ↥Ny × Bool) := ⟨(⟨vx, hvx⟩, ⟨vy, hvy⟩, true)⟩
  set xf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin P) :=
    fun i => if i.2.2 then (i.1 : EuclideanSpace ℝ (Fin P))
      else -(i.1 : EuclideanSpace ℝ (Fin P)) with hxf
  set yf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin p) :=
    fun i => (i.2.1 : EuclideanSpace ℝ (Fin p)) with hyf
  have hxn : ∀ i, ‖xf i‖ = 1 := by
    intro i
    have h1 : ‖(i.1 : EuclideanSpace ℝ (Fin P))‖ = 1 :=
      RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsubx i.1.2)
    by_cases hb : i.2.2 = true <;> simp [hxf, hb, h1]
  have hyn : ∀ i, ‖yf i‖ = 1 := fun i =>
    RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsuby i.2.1.2)
  set u : Matrix (Fin P) (Fin p) ℝ → ℝ :=
    fun B => ⨆ i, (S *ᵥ WithLp.ofLp (xf i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (yf i)) with hu
  have hbdd : ∀ B, BddAbove (Set.range fun i =>
      (S *ᵥ WithLp.ofLp (xf i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (yf i))) :=
    fun B => Finite.bddAbove_range _
  -- every net value is dominated by `u B`
  have hnet : ∀ B : Matrix (Fin P) (Fin p) ℝ, ∀ a ∈ Nnet.domainNet, ∀ b ∈ Nnet.codomainNet,
      |inner ℝ (Matrix.toEuclideanLin (S * B) a) b| ≤ u B := by
    intro B a ha b hb
    have hinner : inner ℝ (Matrix.toEuclideanLin (S * B) a) b
        = (S *ᵥ WithLp.ofLp b) ⬝ᵥ (B *ᵥ WithLp.ofLp a) := by
      rw [inner_euclidean_eq_dotProduct, dotProduct_comm]
      change WithLp.ofLp b ⬝ᵥ ((S * B) *ᵥ WithLp.ofLp a) = _
      rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hS]
    have hpos := le_ciSup (hbdd B) ((⟨b, hb⟩ : ↥Nx), (⟨a, ha⟩ : ↥Ny), true)
    have hneg := le_ciSup (hbdd B) ((⟨b, hb⟩ : ↥Nx), (⟨a, ha⟩ : ↥Ny), false)
    simp only [hxf, hyf, if_true] at hpos hneg
    rw [hinner]
    refine abs_le.mpr ⟨?_, hpos⟩
    have hneg' : -((S *ᵥ WithLp.ofLp b) ⬝ᵥ (B *ᵥ WithLp.ofLp a)) ≤ u B := by
      refine le_trans (le_of_eq ?_) hneg
      simp [Matrix.mulVec_neg, neg_dotProduct]
    linarith
  have hu0 : ∀ B, 0 ≤ u B := by
    intro B
    have := hnet B vy hvy vx hvx
    exact le_trans (abs_nonneg _) this
  -- the net bound, pointwise in `B`
  have hbound : ∀ B : Matrix (Fin P) (Fin p) ℝ, ‖S * B‖ ≤ u B / (1 - 2 * e) := fun B =>
    RMT.matrixOperatorNorm_le_of_centered_bilinear_net (S * B) Nnet he0.le he (hu0 B) (hnet B)
  -- integrability of the process
  set xs : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin P) :=
    fun i => WithLp.toLp 2 (S *ᵥ WithLp.ofLp (xf i)) with hxs
  have hint : Integrable u (gaussianMatrix P p) := by
    have hg : Integrable
        (fun w : EuclideanSpace ℝ (Fin P × Fin p) =>
          ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
            (Matrix.toEuclideanLin (R3.bilMat xs yf)) w) i)
        (stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))) :=
      R3.integrable_ciSup_clm _ _
    have heq : u = (fun w : EuclideanSpace ℝ (Fin P × Fin p) =>
        ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
          (Matrix.toEuclideanLin (R3.bilMat xs yf)) w) i) ∘ (R3.matrixEquivP P p) := by
      funext B
      exact (iSup_congr fun i => R3.bilMat_apply xs yf B i).symm
    rw [heq]
    exact ((R3.measurePreserving_matrixEquivP P p).integrable_comp_emb
      (R3.matrixEquivP P p).measurableEmbedding).mpr hg
  -- integrate
  have hden : (0 : ℝ) < 1 - 2 * e := by linarith
  calc ∫ B, ‖S * B‖ ∂(gaussianMatrix P p)
      ≤ ∫ B, u B / (1 - 2 * e) ∂(gaussianMatrix P p) :=
        integral_mono_of_nonneg (Eventually.of_forall fun B => norm_nonneg _)
          (hint.div_const _) (Eventually.of_forall hbound)
    _ = (∫ B, u B ∂(gaussianMatrix P p)) / (1 - 2 * e) := integral_div _ _
    _ ≤ (Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p) / (1 - 2 * e) := by
        exact div_le_div_of_nonneg_right
          (integral_iSup_bilin_le_het S hS hSL xf yf hxn hyn) hden.le

/-- **R3a-het.** Gordon's bound in expectation for `S * B`, after `ε₀ ↓ 0`:
`E ‖S B‖ ≤ √tr(S Sᵀ) + L √p`. -/
theorem integral_opNorm_mul_le (hP : 0 < P) (hp : 0 < p) (S : Matrix (Fin P) (Fin P) ℝ)
    (hS : Sᵀ = S) {L : ℝ} (hSL : ‖S‖ ≤ L) :
    ∫ B, ‖S * B‖ ∂(gaussianMatrix P p) ≤ Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p := by
  set A : ℝ := Real.sqrt (S * Sᵀ).trace + L * Real.sqrt p with hA
  have hcont : ContinuousAt (fun e : ℝ => A / (1 - 2 * e)) 0 := by
    refine ContinuousAt.div continuousAt_const (by fun_prop) (by norm_num)
  have hf : Tendsto (fun e : ℝ => A / (1 - 2 * e)) (𝓝[>] (0 : ℝ)) (𝓝 A) := by
    have h : Tendsto (fun e : ℝ => A / (1 - 2 * e)) (𝓝[>] (0 : ℝ))
        (𝓝 (A / (1 - 2 * (0 : ℝ)))) := hcont.continuousWithinAt
    simpa using h
  refine ge_of_tendsto hf ?_
  have hsmall : ∀ᶠ e in 𝓝[>] (0 : ℝ), e ∈ Set.Iio (1 / 2 : ℝ) :=
    Filter.Eventually.filter_mono nhdsWithin_le_nhds (Iio_mem_nhds (by norm_num))
  filter_upwards [self_mem_nhdsWithin, hsmall] with e he1 he2
  have he3 : e < 1 / 2 := he2
  exact integral_opNorm_mul_le_aux hP hp S hS hSL he1 (by linarith)

/-! ### R3b-het: one-sided Gaussian concentration of `‖S * B‖` -/

/-- **R3b-het.** `measure_ge_le_of_lipschitz` at `F := fun w => ‖S * (matrixEquivE P p).symm w‖`,
which is `L`-Lipschitz when `‖S‖ ≤ L`. -/
theorem measure_opNorm_mul_ge_le (hP : 0 < P) (hp : 0 < p) (S : Matrix (Fin P) (Fin P) ℝ)
    {L : ℝ≥0} (hL : 0 < L) (hSL : ‖S‖ ≤ L) {t : ℝ} (ht : 0 < t) :
    ((gaussianMatrix P p)
        {B | (∫ B', ‖S * B'‖ ∂(gaussianMatrix P p)) + t ≤ ‖S * B‖}).toReal
      ≤ Real.exp (-t ^ 2 / (2 * (L : ℝ) ^ 2)) := by
  have h := measure_ge_le_of_lipschitz (p := P) (d := p) (LL := L)
    (Nat.mul_pos hP hp) hL (lipschitzWith_opNorm_mul S hSL) ht
  simp only [MeasurableEquiv.symm_apply_apply] at h
  simpa using h

/-! ### R3c-het: the assembly -/

/-- The scalar inequality of R3c-het: `(r + L + κ)² ≤ (L + sc)² + ε`. -/
theorem edge_arith_het {sc L κ ε r : ℝ} (hsc : 0 ≤ sc) (hL : 0 ≤ L) (hκ0 : 0 < κ)
    (hκ1 : κ ≤ 1) (hκK : κ * (4 * (L + sc) + 4) ≤ ε) (hr0 : 0 ≤ r) (hr : r ≤ sc + κ) :
    (r + L + κ) ^ 2 ≤ (L + sc) ^ 2 + ε := by
  have h1 : r + L + κ ≤ L + sc + 2 * κ := by linarith
  have h2 : (0 : ℝ) ≤ r + L + κ := by linarith
  have h3 : (r + L + κ) ^ 2 ≤ (L + sc + 2 * κ) ^ 2 := by nlinarith
  have h4 : κ ^ 2 ≤ κ := by nlinarith
  nlinarith [h3, h4, hκK, hsc, hL, hκ0.le]

/-- **R3-het (H1), block form.** The edge bound for `W₀ = d⁻¹ (S B)(S B)ᵀ` with
`B ~ gaussianMatrix P p`, `p + 1 = d`, `‖S‖ ≤ L`, `Sᵀ = S` and `tr(S Sᵀ)/d → τ`: for every
`ε > 0`, `lamMax W₀ ≤ (L + √τ)² + ε` with probability tending to `1`. The model theorem
`MultiTableModel.tendsto_measure_lamMax_W0het_le_bSF` instantiates `S = SigmaHalf`,
`L = √(wSqMax w)`, `τ = ∑ c_i w_i²`. -/
theorem tendsto_measure_lamMax_le_het
    {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] (μ : ∀ N, Measure (Ω N))
    [∀ N, IsProbabilityMeasure (μ N)] {P p d : ℕ → ℕ}
    (S : (N : ℕ) → Matrix (Fin (P N)) (Fin (P N)) ℝ) {L τ : ℝ} (hL : 0 ≤ L)
    (hS : ∀ N, (S N)ᵀ = S N) (hSL : ∀ N, ‖S N‖ ≤ L)
    (htr0 : ∀ N, 0 ≤ (S N * (S N)ᵀ).trace)
    (htr : Tendsto (fun N => (S N * (S N)ᵀ).trace / d N) atTop (𝓝 τ))
    (B : (N : ℕ) → Ω N → Matrix (Fin (P N)) (Fin (p N)) ℝ)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (P N)) (Fin (P N)) ℝ)
    (hW₀ : ∀ N ω, W₀ N ω = ((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ))
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (hlaw : ∀ N, HasLaw (B N) (gaussianMatrix (P N) (p N)) (μ N))
    (hP : ∀ N, 0 < P N) (hpd : ∀ N, p N + 1 = d N) (hdtop : Tendsto d atTop atTop)
    {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ (L + Real.sqrt τ) ^ 2 + ε})
      atTop (𝓝 1) := by
  classical
  set sc : ℝ := Real.sqrt τ with hscdef
  have hsc0 : (0 : ℝ) ≤ sc := Real.sqrt_nonneg _
  set K : ℝ := 4 * (L + sc) + 4 with hKdef
  have hK0 : (0 : ℝ) < K := by rw [hKdef]; nlinarith
  set κ : ℝ := min 1 (ε / K) with hkapdef
  have hκ0 : (0 : ℝ) < κ := lt_min one_pos (div_pos hε hK0)
  have hκ1 : κ ≤ 1 := min_le_left _ _
  have hκK : κ * K ≤ ε := by
    have hle : κ ≤ ε / K := min_le_right _ _
    have := mul_le_mul_of_nonneg_right hle hK0.le
    rwa [div_mul_cancel₀ _ hK0.ne'] at this
  have hd : ∀ N, 0 < d N := fun N => by rw [← hpd N]; exact Nat.succ_pos _
  -- the Lipschitz constant `L + 1`
  set LL : ℝ≥0 := ⟨L + 1, by positivity⟩ with hLLdef
  have hLL0 : (0 : ℝ≥0) < LL := by
    rw [← NNReal.coe_pos]; change (0 : ℝ) < L + 1; linarith
  have hLLc : (LL : ℝ) = L + 1 := rfl
  have hSLL : ∀ N, ‖S N‖ ≤ (LL : ℝ) := fun N => by rw [hLLc]; linarith [hSL N]
  -- the tail parameter and the bad event
  set t : ℕ → ℝ := fun N => κ * Real.sqrt (d N) with htdef
  have hsd : ∀ N, (0 : ℝ) < Real.sqrt (d N) := fun N =>
    Real.sqrt_pos.mpr (by exact_mod_cast hd N)
  have ht0 : ∀ N, (0 : ℝ) < t N := fun N => by
    rw [htdef]; exact mul_pos hκ0 (hsd N)
  set Bad : (N : ℕ) → Set (Matrix (Fin (P N)) (Fin (p N)) ℝ) := fun N =>
    {Z | (∫ Z', ‖S N * Z'‖ ∂(gaussianMatrix (P N) (p N))) + t N ≤ ‖S N * Z‖} with hBaddef
  have hmeasBad : ∀ N, MeasurableSet (Bad N) := fun N =>
    measurableSet_le measurable_const (measurable_opNorm_mul _)
  -- transfer of the complement to the canonical Gaussian law
  have hmeasure : ∀ N, μ N {ω | B N ω ∈ (Bad N)ᶜ}
      = (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ) := fun N =>
    (hlaw N).measure_eq (hmeasBad N).compl
  have hcompl : ∀ N, (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ)
      = 1 - (gaussianMatrix (P N) (p N)) (Bad N) := fun N =>
    prob_compl_eq_one_sub (hmeasBad N)
  -- eventually `p N > 0`
  have hd2 : ∀ᶠ N in atTop, 2 ≤ d N := hdtop.eventually (eventually_ge_atTop 2)
  have hp0 : ∀ N, 2 ≤ d N → 0 < p N := fun N h => by have := hpd N; omega
  -- the tail bound
  have hsq : ∀ N, (t N) ^ 2 = κ ^ 2 * (d N : ℝ) := fun N => by
    change (κ * Real.sqrt (d N)) ^ 2 = κ ^ 2 * (d N : ℝ)
    rw [mul_pow, Real.sq_sqrt (by positivity)]
  have hbound : ∀ᶠ N in atTop, (gaussianMatrix (P N) (p N)) (Bad N)
      ≤ ENNReal.ofReal (Real.exp (-(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2))) := by
    filter_upwards [hd2] with N hN
    have h := measure_opNorm_mul_ge_le (hP N) (hp0 N hN) (S N) hLL0 (hSLL N) (ht0 N)
    rw [hLLc, hsq N] at h
    have hfin : (gaussianMatrix (P N) (p N)) (Bad N) ≠ ⊤ := measure_ne_top _ _
    rw [← ENNReal.ofReal_toReal hfin]
    exact ENNReal.ofReal_le_ofReal h
  -- the tail bound tends to zero
  have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hdtop
  have hexp0 : Tendsto (fun N => -(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)) atTop atBot := by
    have h2 : Tendsto (fun N => (κ ^ 2 / (2 * (L + 1) ^ 2)) * (d N : ℝ)) atTop atTop :=
      Filter.Tendsto.const_mul_atTop (by positivity) hdR
    have heq : ∀ N : ℕ, -(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)
        = -((κ ^ 2 / (2 * (L + 1) ^ 2)) * (d N : ℝ)) := fun N => by ring
    simp only [heq]
    exact tendsto_neg_atTop_atBot.comp h2
  have hzero : Tendsto (fun N => (gaussianMatrix (P N) (p N)) (Bad N)) atTop (𝓝 0) := by
    have hE : Tendsto (fun N => ENNReal.ofReal
        (Real.exp (-(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)))) atTop (𝓝 0) := by
      have := ENNReal.tendsto_ofReal (Real.tendsto_exp_atBot.comp hexp0)
      simpa using this
    exact tendsto_of_tendsto_of_tendsto_of_le_of_le'
      (tendsto_const_nhds (x := (0 : ℝ≥0∞)) (f := atTop)) hE
      (Eventually.of_forall fun _ => by simp) hbound
  -- the eventual inclusion
  have hrsmall : ∀ᶠ N in atTop,
      Real.sqrt ((S N * (S N)ᵀ).trace / (d N : ℝ)) < sc + κ :=
    Filter.Tendsto.eventually_lt_const (by linarith) htr.sqrt
  have hev : ∀ᶠ N in atTop, (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N)
      ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ (L + sc) ^ 2 + ε} := by
    filter_upwards [hrsmall, hd2] with N hrN hN
    have hsub : {ω | B N ω ∈ (Bad N)ᶜ}
        ⊆ {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ (L + sc) ^ 2 + ε} := by
      intro ω hω
      simp only [hBaddef, Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω
      have hGordon := integral_opNorm_mul_le (hP N) (hp0 N hN) (S N) (hS N) (hSL N)
      -- `‖S B‖ ≤ (r + L + κ) √d` with `r = √(tr(S Sᵀ) / d)`
      set sd : ℝ := Real.sqrt (d N) with hsddef
      set r : ℝ := Real.sqrt ((S N * (S N)ᵀ).trace / (d N : ℝ)) with hrdef
      have hsd0 : (0 : ℝ) < sd := hsd N
      have hr0 : (0 : ℝ) ≤ r := Real.sqrt_nonneg _
      have hsp : r * sd = Real.sqrt (S N * (S N)ᵀ).trace := by
        rw [hrdef, hsddef, Real.sqrt_div (htr0 N), div_mul_cancel₀ _ (hsd N).ne']
      have hpd' : Real.sqrt (p N) ≤ sd := by
        have hle : p N ≤ d N := by rw [← hpd N]; exact Nat.le_succ (p N)
        rw [hsddef]
        exact Real.sqrt_le_sqrt (by exact_mod_cast hle)
      have hM0 : (0 : ℝ) ≤ ‖S N * B N ω‖ := norm_nonneg _
      have hMle : ‖S N * B N ω‖ ≤ (r + L + κ) * sd := by
        have h1 : ‖S N * B N ω‖
            < (∫ Z', ‖S N * Z'‖ ∂(gaussianMatrix (P N) (p N))) + t N := hω
        have h2 : (∫ Z', ‖S N * Z'‖ ∂(gaussianMatrix (P N) (p N)))
            ≤ Real.sqrt (S N * (S N)ᵀ).trace + L * Real.sqrt (p N) := hGordon
        have h3 : t N = κ * sd := rfl
        nlinarith [h1, h2, h3, hsp, mul_le_mul_of_nonneg_left hpd' hL]
      -- the Rayleigh bound
      have hherm : (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)).IsHermitian :=
        (hW₀ N ω) ▸ (hsymm N ω)
      have hlm : lamMax (W₀ N ω) (hsymm N ω)
          = lamMax (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)) hherm :=
        lamMax_congr (hW₀ N ω) _ _
      have hray : lamMax (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)) hherm
          ≤ ((d N : ℝ))⁻¹ * ‖S N * B N ω‖ ^ 2 :=
        lamMax_gramT_le_opNorm_sq (by positivity) _ hherm
      -- the arithmetic
      have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hsdsq : sd ^ 2 = (d N : ℝ) := Real.sq_sqrt hdpos.le
      have hrk : (0 : ℝ) ≤ r + L + κ := by linarith
      have hMsq : ‖S N * B N ω‖ ^ 2 ≤ ((r + L + κ) * sd) ^ 2 := pow_le_pow_left₀ hM0 hMle 2
      have hfinal : ((d N : ℝ))⁻¹ * ‖S N * B N ω‖ ^ 2 ≤ (r + L + κ) ^ 2 := by
        rw [inv_mul_le_iff₀ hdpos]
        calc ‖S N * B N ω‖ ^ 2 ≤ ((r + L + κ) * sd) ^ 2 := hMsq
          _ = (r + L + κ) ^ 2 * (d N : ℝ) := by rw [mul_pow, hsdsq]
          _ = (d N : ℝ) * (r + L + κ) ^ 2 := by ring
      have hedge := edge_arith_het hsc0 hL hκ0 hκ1 hκK hr0 hrN.le
      simp only [Set.mem_ofPred_eq]
      rw [hlm]
      linarith [hray, hfinal, hedge]
    calc (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N)
        = (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ) := (hcompl N).symm
      _ = μ N {ω | B N ω ∈ (Bad N)ᶜ} := (hmeasure N).symm
      _ ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ (L + sc) ^ 2 + ε} := measure_mono hsub
  -- squeeze
  have hlow : Tendsto (fun N => (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N))
      atTop (𝓝 1) := by
    have := ENNReal.Tendsto.sub (tendsto_const_nhds (x := (1 : ℝ≥0∞)) (f := atTop))
      hzero (Or.inl (by simp))
    simpa using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds hev
    (Eventually.of_forall fun N => prob_le_one)

end R3het

/-! ### The model theorem -/

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

omit [NeZero M] in
/-- `tr(Σ^{1/2} Σ^{1/2}ᵀ) = ∑ n_i w_i²`: only `w_i²` enters. -/
theorem trace_SigmaHalf_mul_transpose (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    (m.SigmaHalf w N * (m.SigmaHalf w N)ᵀ).trace = ∑ i, (n i N : ℝ) * w i ^ 2 := by
  rw [SigmaHalf, Matrix.diagonal_transpose, Matrix.diagonal_mul_diagonal, Matrix.trace_diagonal,
    ← Equiv.sum_comp finSigmaFinEquiv]
  simp only [Equiv.symm_apply_apply]
  rw [Fintype.sum_sigma]
  simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, sq]

omit [NeZero M] in
/-- `‖Σ^{1/2}‖ ≤ √(wSqMax w) = max_i |w_i|`. -/
theorem opNorm_SigmaHalf_le (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    ‖m.SigmaHalf w N‖ ≤ Real.sqrt (Scalars.wSqMax w) := by
  refine R3het.opNorm_diagonal_le _ (Real.sqrt_nonneg _) fun r => ?_
  rw [← Real.sqrt_sq_eq_abs]
  exact Real.sqrt_le_sqrt (Scalars.le_wSqMax w _)

omit [NeZero M] in
/-- `‖Σ^{1/2} B‖ ≤ √(wSqMax w) ‖B‖`, the Lipschitz constant of item 3 of the header. -/
theorem opNorm_SigmaHalf_mul_le (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    {q : ℕ} (B : Matrix (Fin (∑ i, n i N)) (Fin q) ℝ) :
    ‖m.SigmaHalf w N * B‖ ≤ Real.sqrt (Scalars.wSqMax w) * ‖B‖ :=
  R3het.opNorm_mul_le_of_le _ (m.opNorm_SigmaHalf_le w N) B

/-- **Item H8, the `edge` field of `ResolventLimitsHet` at `b = bSF c w`.** For every
`ε > 0`, `lamMax W₀' ≤ bSF c w + ε` with probability tending to `1`. Of `hreg` only
`d → ∞` and `n_i / d → c_i` are used; `hG` gives the block law through
`exists_block_hasLaw_het`. No sign or nonvanishing condition on `w` is needed. -/
theorem tendsto_measure_lamMax_W0het_le_bSF [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ MPhet.bSF c w + ε})
      atTop (𝓝 1) := by
  intro ε hε
  have hd : ∀ N, 0 < d N := (m.tbl 0).hd
  have hp : ∀ N, d N = (d N - 1) + 1 := fun N => (Nat.sub_add_cancel (hd N)).symm
  choose B hB hlaw using fun N => m.exists_block_hasLaw_het N hG (hp N)
  have hW : ∀ N ω, m.W0het w N ω = ((d N : ℝ))⁻¹ •
      ((m.SigmaHalf w N * B N ω) * (m.SigmaHalf w N * B N ω)ᵀ) :=
    fun N ω => m.W0het_eq_of_block w N ω (B N ω) (hB N ω)
  have hlawB : ∀ N, HasLaw (B N) (gaussianMatrix (∑ i, n i N) (d N - 1)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hlaw N)
  have htr : Tendsto (fun N => (m.SigmaHalf w N * (m.SigmaHalf w N)ᵀ).trace / (d N : ℝ))
      atTop (𝓝 (∑ i, c i * w i ^ 2)) := by
    simp only [trace_SigmaHalf_mul_transpose, Finset.sum_div]
    refine tendsto_finsetSum _ fun i _ => ?_
    refine ((hreg i).2.2.mul_const (w i ^ 2)).congr fun N => ?_
    ring
  have hbSF : MPhet.bSF c w
      = (Real.sqrt (Scalars.wSqMax w) + Real.sqrt (∑ i, c i * w i ^ 2)) ^ 2 := rfl
  rw [hbSF]
  exact R3het.tendsto_measure_lamMax_le_het μ (fun N => m.SigmaHalf w N) (Real.sqrt_nonneg _)
    (fun N => m.transpose_SigmaHalf w N) (fun N => m.opNorm_SigmaHalf_le w N)
    (fun N => by rw [trace_SigmaHalf_mul_transpose]; positivity) htr B (m.W0het w) hW
    (m.isHermitian_W0het w) hlawB m.stack_row_pos (fun N => Nat.sub_add_cancel (hd N))
    (hreg 0).2.1 hε

end MultiTableModel
end StackedSVD
