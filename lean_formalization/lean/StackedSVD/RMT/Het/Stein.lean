/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.SteinStep
import StackedSVD.RMT.R1

/-!
# Item H5: the Stein system for a heteroscedastic Wishart resolvent

Task H5 of `notes/archive/plan_heterolaw_A.md` (section 3.4, the H5 row of section 4, and fallback 2
of risk 2). The file is stated for an abstract row-scale profile and never mentions tables or
blocks; blocks enter in H6 through `τ_j = w_{i(j)}`.

Model. `B : Matrix (Fin p) (Fin q) ℝ` has i.i.d. `N(0,1)` entries (`gaussianMatrix p q`,
`Defs.lean`: no scale inside the law). The row scale `τ : Fin p → ℝ` is the **signed** square
root of the variance profile, `σ_j = τ_j²`, so that the identification with
`Het/Split.lean` is definitional: `W0het_eq_of_block` writes `W₀'` as
`d⁻¹ • ((SigmaHalf w N * B) * (SigmaHalf w N * B)ᵀ)` and `SigmaHalf w N = diagonal τ` for
`τ r = w (finSigmaFinEquiv.symm r).1`. With `Y = diagonal τ * B` (row `j` of `d^{-1/2} Y` has
entries of variance `τ_j² / d`; `d` is a separate scale parameter, `d = q + 1` in H0):

```
W_σ  = d⁻¹ Y Yᵀ   (p × p),   G_σ(z) = (W_σ - z)⁻¹,   g_j(z) = G_σ(z)_{jj},
W_σ' = d⁻¹ Yᵀ Y   (q × q),   s(z)   = d⁻¹ tr (W_σ' - z)⁻¹     (the `d`-side trace of 2.2).
```

Results (all for `0 < Im z`, `0 < d`).

1. `trace_identity`: `s = d⁻¹ tr G_σ - (q - p) / (d z)`, exact. Grouped by blocks this is
   the plan's `s^N = ∑ c_i^N g_i^N - (1 - ∑ c_i^N)/z` with `c_i^N = n_i / d` and the finite
   `d`-side count `q = d - 1` kept honest.
2. `stein_row`: `E[1 + z g_j] = -τ_j² z E[s g_j] + r_j` with
   `r_j = -τ_j² d⁻¹ E[g_j + z (G_σ²)_{jj}]` and `‖r_j‖ ≤ τ_j² d⁻¹ (η⁻¹ + ‖z‖ η⁻²)`, `η = Im z`.
   The residual is `O(1/(d η²))`, two powers of `η` better than the H5 row asks.
3. `stein_block`: the same identity for the average over a finite set of rows of equal
   `τ_j²`, same residual bound.
4. `variance_re_gsig_le` and friends: `Var (Re g_j) ≤ 4 lipG²` with
   `lipG = 2 √(S (η + ‖z‖) / d) / η²`, and `Var (Re s) ≤ 4 lipS²` with
   `lipS = 2 √(S p (η + ‖z‖) / d) / (d η²)`, for any `S ≥ τ_j²`; same for `Im`.

Method. The Gaussian integration by parts of `RMT/SteinStep.lean` is reused through the
linear map `lam τ : B ↦ (diagonal τ * B)ᵀ` on the flattened spaces: with `κ = d / p` and
`ζ = κ z`, `G_σ(z) = κ · Gmat ζ (lam τ B)` in the `q × p` picture of `ResolvDeriv`
(`Gsig_eq_smul_Gmat`), and every derivative is a chain rule through `lam τ`.

Numeric check: a session script (`check_stein.py`, seed 20260830; not kept).
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix NNReal

namespace StackedSVD

namespace HetStein

open ResolvDeriv SteinStep R4

variable {p q d : ℕ} {z : ℂ}

/-! ### The model -/

/-- `Y = diagonal τ * B`: row `j` of `B` scaled by `τ_j`. -/
noncomputable def Ysig (τ : Fin p → ℝ) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix (Fin p) (Fin q) ℝ :=
  Matrix.diagonal τ * B

/-- `W_σ = d⁻¹ Y Yᵀ`, the `p × p` heteroscedastic Wishart matrix (plan section 2.1, `W₀'`). -/
noncomputable def Wsig (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix (Fin p) (Fin p) ℝ :=
  (d : ℝ)⁻¹ • (Ysig τ B * (Ysig τ B)ᵀ)

/-- `W_σ' = d⁻¹ Yᵀ Y`, the `q × q` companion (the `d`-side matrix of plan section 2.2). -/
noncomputable def Wsig' (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix (Fin q) (Fin q) ℝ :=
  (d : ℝ)⁻¹ • ((Ysig τ B)ᵀ * Ysig τ B)

/-- `G_σ(z) = (W_σ - z)⁻¹`. -/
noncomputable def Gsig (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix (Fin p) (Fin p) ℂ :=
  R4C.resolvC (Wsig τ d B) z

/-- `g_j(z) = G_σ(z)_{jj}`, the row trace. -/
noncomputable def gsig (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (B : Matrix (Fin p) (Fin q) ℝ)
    (j : Fin p) : ℂ :=
  Gsig τ d z B j j

/-- `s(z) = d⁻¹ tr (W_σ' - z)⁻¹`, the `d`-side trace of plan section 2.2. -/
noncomputable def ssig (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (B : Matrix (Fin p) (Fin q) ℝ) : ℂ :=
  (d : ℂ)⁻¹ * (R4C.resolvC (Wsig' τ d B) z).trace

theorem isHermitian_Wsig (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    (Wsig τ d B).IsHermitian := by
  have h : (Ysig τ B * (Ysig τ B)ᵀ).IsHermitian := by
    simpa using Matrix.isHermitian_mul_conjTranspose_self (Ysig τ B)
  change ((d : ℝ)⁻¹ • (Ysig τ B * (Ysig τ B)ᵀ))ᴴ = (d : ℝ)⁻¹ • (Ysig τ B * (Ysig τ B)ᵀ)
  rw [Matrix.conjTranspose_smul, star_trivial, h.eq]

theorem isHermitian_Wsig' (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    (Wsig' τ d B).IsHermitian := by
  have h : ((Ysig τ B)ᵀ * Ysig τ B).IsHermitian := by
    simpa using Matrix.isHermitian_conjTranspose_mul_self (Ysig τ B)
  change ((d : ℝ)⁻¹ • ((Ysig τ B)ᵀ * Ysig τ B))ᴴ = (d : ℝ)⁻¹ • ((Ysig τ B)ᵀ * Ysig τ B)
  rw [Matrix.conjTranspose_smul, star_trivial, h.eq]

/-! ### The trace identity (plan section 2.2, finite `N`) -/

private theorem cmat_smul (c : ℝ) {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) :
    R4C.cmat (c • M) = (c : ℂ) • R4C.cmat M := by
  ext i j
  simp [R4C.cmat, Matrix.smul_apply]

private theorem cmat_Wsig (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    R4C.cmat (Wsig τ d B) = (d : ℂ)⁻¹ • (cmapR (Ysig τ B) * (cmapR (Ysig τ B))ᵀ) := by
  rw [Wsig, cmat_smul, R4C.cmat, ← cmapR_transpose, ← cmapR_mul]
  push_cast
  rfl

private theorem cmat_Wsig' (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    R4C.cmat (Wsig' τ d B) = (d : ℂ)⁻¹ • ((cmapR (Ysig τ B))ᵀ * cmapR (Ysig τ B)) := by
  rw [Wsig', cmat_smul, R4C.cmat, ← cmapR_transpose, ← cmapR_mul]
  push_cast
  rfl

/-- `Y Yᵀ G = d (1 + z G)`: the resolvent identity of `W_σ` in the form the algebra uses. -/
theorem Yc_mul_Yct_mul_Gsig (hz : z.im ≠ 0) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    cmapR (Ysig τ B) * ((cmapR (Ysig τ B))ᵀ * Gsig τ d z B)
      = (d : ℂ) • ((1 : Matrix (Fin p) (Fin p) ℂ) + z • Gsig τ d z B) := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hAG : (R4C.cmat (Wsig τ d B) - z • (1 : Matrix (Fin p) (Fin p) ℂ)) * Gsig τ d z B = 1 :=
    cmat_sub_mul_resolvC (isHermitian_Wsig τ d B) hz
  rw [cmat_Wsig] at hAG
  simp only [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul, Matrix.mul_assoc] at hAG
  have h1 : (d : ℂ)⁻¹ • (cmapR (Ysig τ B) * ((cmapR (Ysig τ B))ᵀ * Gsig τ d z B))
      = 1 + z • Gsig τ d z B := by
    rw [← sub_eq_iff_eq_add.mp hAG]
  rw [← h1, smul_smul, mul_inv_cancel₀ hdC, one_smul]

/-- The push-through identity: `(W_σ' - z)⁻¹ = -z⁻¹ (1 - d⁻¹ Yᵀ G_σ Y)`. -/
theorem resolvC_Wsig'_eq (hz : z.im ≠ 0) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    R4C.resolvC (Wsig' τ d B) z
      = -z⁻¹ • ((1 : Matrix (Fin q) (Fin q) ℂ)
          - (d : ℂ)⁻¹ • ((cmapR (Ysig τ B))ᵀ * (Gsig τ d z B * cmapR (Ysig τ B)))) := by
  have hz0 : z ≠ 0 := fun h => hz (by rw [h]; simp)
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have h2 := Yc_mul_Yct_mul_Gsig hz hd τ B
  set Yc := cmapR (Ysig τ B) with hYc
  set G := Gsig τ d z B with hG
  have h2' : Yc * (Ycᵀ * (G * Yc)) = (d : ℂ) • (Yc + z • (G * Yc)) := by
    calc Yc * (Ycᵀ * (G * Yc)) = (Yc * (Ycᵀ * G)) * Yc := by simp only [Matrix.mul_assoc]
      _ = (d : ℂ) • (Yc + z • (G * Yc)) := by
          rw [h2, Matrix.smul_mul, Matrix.add_mul, Matrix.one_mul, Matrix.smul_mul]
  have hNM : (Ycᵀ * Yc) * (Ycᵀ * (G * Yc))
      = (d : ℂ) • ((Ycᵀ * Yc) + z • (Ycᵀ * (G * Yc))) := by
    calc (Ycᵀ * Yc) * (Ycᵀ * (G * Yc)) = Ycᵀ * (Yc * (Ycᵀ * (G * Yc))) := by
          simp only [Matrix.mul_assoc]
      _ = Ycᵀ * ((d : ℂ) • (Yc + z • (G * Yc))) := by rw [h2']
      _ = (d : ℂ) • ((Ycᵀ * Yc) + z • (Ycᵀ * (G * Yc))) := by
          rw [Matrix.mul_smul, Matrix.mul_add, Matrix.mul_smul]
  refine Matrix.inv_eq_right_inv ?_
  rw [cmat_Wsig']
  set N := Ycᵀ * Yc with hN
  set M := Ycᵀ * (G * Yc) with hM
  rw [Matrix.mul_smul, Matrix.mul_sub, Matrix.sub_mul, Matrix.sub_mul, Matrix.mul_one,
    Matrix.mul_one, Matrix.smul_mul, Matrix.mul_smul, Matrix.smul_mul, Matrix.one_mul, hNM,
    smul_smul, smul_smul]
  have hdd : (d : ℂ)⁻¹ * (d : ℂ)⁻¹ * (d : ℂ) = (d : ℂ)⁻¹ := by field_simp
  rw [hdd]
  match_scalars <;> (field_simp; try ring)

/-- **(H5, `trace_identity`).** `s = d⁻¹ tr G_σ - (q - p) / (d z)`. -/
theorem trace_identity (hz : z.im ≠ 0) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    ssig τ d z B = (d : ℂ)⁻¹ * (Gsig τ d z B).trace - ((q : ℂ) - p) / ((d : ℂ) * z) := by
  have hz0 : z ≠ 0 := fun h => hz (by rw [h]; simp)
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have h2 := Yc_mul_Yct_mul_Gsig hz hd τ B
  set Yc := cmapR (Ysig τ B) with hYc
  set G := Gsig τ d z B with hG
  have htr : (Ycᵀ * (G * Yc)).trace = (d : ℂ) * ((p : ℂ) + z * G.trace) := by
    rw [← Matrix.mul_assoc, Matrix.trace_mul_comm, h2, Matrix.trace_smul, Matrix.trace_add,
      Matrix.trace_one, Matrix.trace_smul, Fintype.card_fin, smul_eq_mul, smul_eq_mul]
  rw [ssig, resolvC_Wsig'_eq hz hd, Matrix.trace_smul, Matrix.trace_sub, Matrix.trace_one,
    Matrix.trace_smul, Fintype.card_fin, htr, smul_eq_mul, smul_eq_mul]
  field_simp
  ring

/-! ### The flattening map `lam τ : B ↦ (diagonal τ * B)ᵀ`

`RMT/ResolvDeriv.lean` differentiates `Y ↦ (d⁻¹ Yᵀ Y - z)⁻¹` in the entries of `Y`. For the
`p × p` matrix `W_σ = d⁻¹ Y Yᵀ` that picture applies to `Yᵀ` (a `q × p` matrix, `gram`
normalization `p⁻¹`), and the Gaussian variable is `B`, not `Y`. The linear map below carries
the flattened `B` to the flattened `Yᵀ`; every derivative is a chain rule through it. -/

section Flatten

/-- `lam τ q` as a linear map: `(lam τ q x)_{(k, j)} = τ_j x_{(j, k)}`. -/
noncomputable def lamL (τ : Fin p → ℝ) (q : ℕ) :
    EuclideanSpace ℝ (Fin (p * q)) →ₗ[ℝ] EuclideanSpace ℝ (Fin (q * p)) where
  toFun x := WithLp.toLp 2 fun r : Fin (q * p) =>
    τ (finProdFinEquiv.symm r).2
      * x (finProdFinEquiv ((finProdFinEquiv.symm r).2, (finProdFinEquiv.symm r).1))
  map_add' x y := by
    ext r
    simp [mul_add]
  map_smul' c x := by
    ext r
    simp [mul_left_comm]

/-- `lam τ q : B ↦ (diagonal τ * B)ᵀ` on the flattened spaces (`lam_symm`). -/
noncomputable def lam (τ : Fin p → ℝ) (q : ℕ) :
    EuclideanSpace ℝ (Fin (p * q)) →L[ℝ] EuclideanSpace ℝ (Fin (q * p)) :=
  LinearMap.toContinuousLinearMap (lamL τ q)

theorem lam_apply (τ : Fin p → ℝ) (x : EuclideanSpace ℝ (Fin (p * q))) (k : Fin q)
    (j : Fin p) :
    lam τ q x (finProdFinEquiv (k, j)) = τ j * x (finProdFinEquiv (j, k)) := by
  change τ (finProdFinEquiv.symm (finProdFinEquiv (k, j))).2
      * x (finProdFinEquiv ((finProdFinEquiv.symm (finProdFinEquiv (k, j))).2,
          (finProdFinEquiv.symm (finProdFinEquiv (k, j))).1)) = _
  rw [Equiv.symm_apply_apply]

theorem lam_symm (τ : Fin p → ℝ) (x : EuclideanSpace ℝ (Fin (p * q))) :
    (matrixEquivE q p).symm (lam τ q x)
      = (Matrix.diagonal τ * (matrixEquivE p q).symm x)ᵀ := by
  ext k j
  rw [matrixEquivE_symm_apply, lam_apply, Matrix.transpose_apply, Matrix.diagonal_mul,
    matrixEquivE_symm_apply]

theorem lam_matrixEquivE (τ : Fin p → ℝ) (B : Matrix (Fin p) (Fin q) ℝ) :
    (matrixEquivE q p).symm (lam τ q (matrixEquivE p q B)) = (Ysig τ B)ᵀ := by
  rw [lam_symm, MeasurableEquiv.symm_apply_apply]
  rfl

theorem lam_single (τ : Fin p → ℝ) (j : Fin p) (k : Fin q) :
    lam τ q (EuclideanSpace.single (finProdFinEquiv (j, k)) (1 : ℝ))
      = τ j • EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ) := by
  ext r
  obtain ⟨⟨k', j'⟩, rfl⟩ := finProdFinEquiv.surjective r
  rw [lam_apply]
  simp only [PiLp.single_apply, PiLp.smul_apply, smul_eq_mul]
  by_cases h : (k', j') = (k, j)
  · obtain ⟨rfl, rfl⟩ := Prod.mk.inj h
    simp
  · have h1 : finProdFinEquiv (k', j') ≠ finProdFinEquiv (k, j) :=
      fun h' => h (finProdFinEquiv.injective h')
    have h2 : finProdFinEquiv (j', k') ≠ finProdFinEquiv (j, k) := fun h' => by
      apply h
      have := finProdFinEquiv.injective h'
      simp only [Prod.mk.injEq] at this
      exact Prod.ext this.2 this.1
    simp [h1, h2]

/-- `‖lam τ q x‖ ≤ √S ‖x‖` when `τ_j² ≤ S` for every `j`. -/
theorem norm_lam_le (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S)
    (x : EuclideanSpace ℝ (Fin (p * q))) :
    ‖lam τ q x‖ ≤ Real.sqrt S * ‖x‖ := by
  have h1 : ‖lam τ q x‖ ^ 2 ≤ S * ‖x‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq, EuclideanSpace.real_norm_sq_eq,
      ← Equiv.sum_comp finProdFinEquiv fun r => (lam τ q x r) ^ 2,
      ← Equiv.sum_comp finProdFinEquiv fun r => (x r) ^ 2,
      Fintype.sum_prod_type, Fintype.sum_prod_type, Finset.mul_sum, Finset.sum_comm]
    refine Finset.sum_le_sum fun j _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun k _ => ?_
    rw [lam_apply, mul_pow]
    exact mul_le_mul_of_nonneg_right (hS j) (sq_nonneg _)
  calc ‖lam τ q x‖ = Real.sqrt (‖lam τ q x‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt (S * ‖x‖ ^ 2) := Real.sqrt_le_sqrt h1
    _ = Real.sqrt S * ‖x‖ := by rw [Real.sqrt_mul hS0, Real.sqrt_sq (norm_nonneg _)]

theorem opNorm_lam_le (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) :
    ‖lam τ q‖ ≤ Real.sqrt S :=
  ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg S) (norm_lam_le τ hS hS0)

end Flatten

/-! ### `G_σ` in the `q × p` picture of `ResolvDeriv` -/

section Translate

/-- `κ = d / p`: `gram (Yᵀ) = p⁻¹ Y Yᵀ = κ W_σ`. -/
noncomputable def kap (p d : ℕ) : ℝ := (d : ℝ) / p

theorem kap_pos (hp : 0 < p) (hd : 0 < d) : 0 < kap p d := by
  unfold kap
  have : (0 : ℝ) < p := by exact_mod_cast hp
  have : (0 : ℝ) < d := by exact_mod_cast hd
  positivity

theorem gram_lam (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (B : Matrix (Fin p) (Fin q) ℝ) :
    gram ((matrixEquivE q p).symm (lam τ q (matrixEquivE p q B))) = kap p d • Wsig τ d B := by
  rw [lam_matrixEquivE, gram, Matrix.transpose_transpose, Wsig, smul_smul, kap]
  congr 1
  have hp' : (p : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hd' : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  field_simp

private theorem resolvC_smul (hz : z.im ≠ 0) {n : ℕ} {W : Matrix (Fin n) (Fin n) ℝ}
    (hW : W.IsHermitian) {c : ℝ} (hc : c ≠ 0) :
    R4C.resolvC (c • W) ((c : ℂ) * z) = ((c : ℂ)⁻¹) • R4C.resolvC W z := by
  have hcC : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc
  refine Matrix.inv_eq_right_inv ?_
  have h1 : R4C.cmat (c • W) - ((c : ℂ) * z) • (1 : Matrix (Fin n) (Fin n) ℂ)
      = (c : ℂ) • (R4C.cmat W - z • 1) := by
    rw [cmat_smul, smul_sub, smul_smul]
  rw [h1, Matrix.smul_mul, Matrix.mul_smul, smul_smul, mul_inv_cancel₀ hcC, one_smul,
    cmat_sub_mul_resolvC hW hz]

/-- **The translation.** `G_σ(z) = κ · Gmat (κ z) (lam τ B)`. -/
theorem Gsig_eq_smul_Gmat (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    Gsig τ d z B
      = ((kap p d : ℝ) : ℂ) • Gmat (((kap p d : ℝ) : ℂ) * z) (lam τ q (matrixEquivE p q B)) := by
  have hk : kap p d ≠ 0 := (kap_pos hp hd).ne'
  have hkC : ((kap p d : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hk
  rw [Gmat, gram_lam hp hd, resolvC_smul hz (isHermitian_Wsig τ d B) hk, smul_smul,
    mul_inv_cancel₀ hkC, one_smul]
  rfl

theorem im_kap_mul (z : ℂ) :
    (((kap p d : ℝ) : ℂ) * z).im = kap p d * z.im := by
  simp

theorem norm_kap_mul (hp : 0 < p) (hd : 0 < d) (z : ℂ) :
    ‖((kap p d : ℝ) : ℂ) * z‖ = kap p d * ‖z‖ := by
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (kap_pos hp hd)]

end Translate

/-! ### Column identities in the generic `ResolvDeriv` picture

`Y : Matrix (Fin m) (Fin n) ℝ` flattened, `G = (n⁻¹ Yᵀ Y - z)⁻¹`, `P = Y G`, `Q = Y G Yᵀ`.
`SteinStep` sums its identities over all entries; the row identity needs the sums over one
column `i`. The two helpers `transpose_Ymat_mul_Ymat'`, `cmat_gram_mul_Gmat'` are private in
`SteinStep`; the proofs are the same. -/

section Column

variable {m n : ℕ}

private theorem cmapR_smul' {a b : Type*} (c : ℝ) (M : Matrix a b ℝ) :
    cmapR (c • M) = (c : ℂ) • cmapR M := by
  ext i j
  simp only [cmapR_apply, Matrix.smul_apply, smul_eq_mul, Complex.ofReal_mul]

theorem transpose_Ymat_mul_Ymat' (hn : 0 < n) (x : EuclideanSpace ℝ (Fin (m * n))) :
    (Ymat x)ᵀ * Ymat x = (n : ℂ) • R4C.cmat (gram ((matrixEquivE m n).symm x)) := by
  have hn0 : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  set Y := (matrixEquivE m n).symm x with hY
  have h1 : Yᵀ * Y = (n : ℝ) • gram Y := by
    rw [gram, smul_smul, mul_inv_cancel₀ hn0, one_smul]
  have hcast : (((n : ℕ) : ℝ) : ℂ) = ((n : ℕ) : ℂ) := by push_cast; ring
  calc (Ymat x)ᵀ * Ymat x = cmapR Yᵀ * cmapR Y := by rw [Ymat, cmapR_transpose]
    _ = cmapR (Yᵀ * Y) := (cmapR_mul _ _).symm
    _ = cmapR ((n : ℝ) • gram Y) := by rw [h1]
    _ = (((n : ℕ) : ℝ) : ℂ) • cmapR (gram Y) := cmapR_smul' _ _
    _ = (n : ℂ) • R4C.cmat (gram Y) := by rw [hcast, cmat_eq_cmapR]

theorem cmat_gram_mul_Gmat' (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (m * n))) :
    R4C.cmat (gram ((matrixEquivE m n).symm x)) * Gmat z x = 1 + z • Gmat z x := by
  have h := cmat_sub_mul_resolvC (isHermitian_gram ((matrixEquivE m n).symm x)) (z := z) hz
  rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul] at h
  have h' : R4C.cmat (gram ((matrixEquivE m n).symm x)) * Gmat z x - z • Gmat z x = 1 := h
  rw [sub_eq_iff_eq_add] at h'
  exact h'

/-- `∑_k Y_{ki} P_{ki} = n (1 + z G_{ii})`: the left side of the row identity. -/
theorem sum_Ymat_mul_Pmat_col (hz : z.im ≠ 0) (hn : 0 < n)
    (x : EuclideanSpace ℝ (Fin (m * n))) (i : Fin n) :
    ∑ k : Fin m, Ymat x k i * Pmat z x k i = (n : ℂ) * (1 + z * Gmat z x i i) := by
  have h1 : ∑ k : Fin m, Ymat x k i * Pmat z x k i = ((Ymat x)ᵀ * Pmat z x) i i := by
    rw [Matrix.mul_apply]
    rfl
  rw [h1, Pmat, ← Matrix.mul_assoc, transpose_Ymat_mul_Ymat' hn, Matrix.smul_mul,
    cmat_gram_mul_Gmat' hz, Matrix.smul_apply, Matrix.add_apply, Matrix.one_apply_eq,
    Matrix.smul_apply, smul_eq_mul, smul_eq_mul]

/-- `∑_k P_{ki}² = n (G_{ii} + z (G²)_{ii})`. -/
theorem sum_Pmat_sq_col (hz : z.im ≠ 0) (hn : 0 < n)
    (x : EuclideanSpace ℝ (Fin (m * n))) (i : Fin n) :
    ∑ k : Fin m, Pmat z x k i * Pmat z x k i
      = (n : ℂ) * (Gmat z x i i + z * (Gmat z x * Gmat z x) i i) := by
  have h1 : ∑ k : Fin m, Pmat z x k i * Pmat z x k i = ((Pmat z x)ᵀ * Pmat z x) i i := by
    rw [Matrix.mul_apply]
    rfl
  have hPP : (Pmat z x)ᵀ * Pmat z x
      = (n : ℂ) • (Gmat z x + z • (Gmat z x * Gmat z x)) := by
    have h2 : (Pmat z x)ᵀ = Gmat z x * (Ymat x)ᵀ := by
      rw [Pmat, Matrix.transpose_mul, transpose_Gmat hz]
    rw [h2, Pmat]
    calc Gmat z x * (Ymat x)ᵀ * (Ymat x * Gmat z x)
        = Gmat z x * (((Ymat x)ᵀ * Ymat x) * Gmat z x) := by simp only [Matrix.mul_assoc]
      _ = Gmat z x * (((n : ℂ) • R4C.cmat (gram ((matrixEquivE m n).symm x))) * Gmat z x) := by
          rw [transpose_Ymat_mul_Ymat' hn]
      _ = (n : ℂ) • (Gmat z x * (R4C.cmat (gram ((matrixEquivE m n).symm x)) * Gmat z x)) := by
          rw [Matrix.smul_mul, Matrix.mul_smul]
      _ = (n : ℂ) • (Gmat z x * ((1 : Matrix (Fin n) (Fin n) ℂ) + z • Gmat z x)) := by
          rw [cmat_gram_mul_Gmat' hz]
      _ = (n : ℂ) • (Gmat z x + z • (Gmat z x * Gmat z x)) := by
          rw [Matrix.mul_add, Matrix.mul_one, Matrix.mul_smul]
  rw [h1, hPP, Matrix.smul_apply, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, smul_eq_mul]

/-- **The divergence over one column.**
`∑_k ∂F_{ki}/∂Y_{ki} = m G_{ii} - (n + z tr G) G_{ii} - (G_{ii} + z (G²)_{ii})`. -/
theorem sum_fderiv_Fentry_single_col (hz : z.im ≠ 0) (hn : 0 < n)
    (x : EuclideanSpace ℝ (Fin (m * n))) (i : Fin n) :
    ∑ k : Fin m, fderiv ℝ (Fentry (p := m) (d := n) z k i) x
        (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
      = (m : ℂ) * Gmat z x i i - ((n : ℂ) + z * (Gmat z x).trace) * Gmat z x i i
        - (Gmat z x i i + z * (Gmat z x * Gmat z x) i i) := by
  have hnC : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  have hQ : (Qmat z x).trace = ∑ k : Fin m, Qmat z x k k := by
    rw [Matrix.trace]
    exact Finset.sum_congr rfl fun k _ => rfl
  have expand : ∀ k : Fin m,
      fderiv ℝ (Fentry (p := m) (d := n) z k i) x
          (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
        = Gmat z x i i - (n : ℂ)⁻¹ * (Qmat z x k k * Gmat z x i i)
          - (n : ℂ)⁻¹ * (Pmat z x k i * Pmat z x k i) := by
    intro k
    rw [fderiv_Fentry_single hz k i x]
    ring
  rw [Finset.sum_congr rfl fun k _ => expand k, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, ← Finset.mul_sum,
    ← Finset.mul_sum, ← Finset.sum_mul, ← hQ, trace_Qmat hz hn x, sum_Pmat_sq_col hz hn x i]
  field_simp

end Column

/-! ### The Gaussian integration by parts, one row at a time -/

section SteinRow

/-- `ζ = κ z`, the spectral parameter in the `q × p` picture. -/
noncomputable def zeta (p d : ℕ) (z : ℂ) : ℂ := ((kap p d : ℝ) : ℂ) * z

theorem zeta_im_pos (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) : 0 < (zeta p d z).im := by
  rw [zeta, im_kap_mul]
  exact mul_pos (kap_pos hp hd) hz

theorem contDiff_Fsig (hζ : (zeta p d z).im ≠ 0) (τ : Fin p → ℝ) (k : Fin q) (j : Fin p) :
    ContDiff ℝ 1 (fun x : EuclideanSpace ℝ (Fin (p * q)) =>
      Fentry (zeta p d z) k j (lam τ q x)) :=
  (contDiff_Fentry hζ k j).comp (lam τ q).contDiff

theorem norm_fderiv_Fsig_le (hζ : 0 < (zeta p d z).im) (hp : 0 < p) (τ : Fin p → ℝ)
    (k : Fin q) (j : Fin p) (x : EuclideanSpace ℝ (Fin (p * q))) :
    ‖fderiv ℝ (fun x : EuclideanSpace ℝ (Fin (p * q)) => Fentry (zeta p d z) k j (lam τ q x)) x‖
      ≤ Lconst (zeta p d z) q p * ‖lam τ q‖ := by
  have hdiff := differentiableAt_Fentry hζ.ne' k j (lam τ q x)
  rw [show (fun x : EuclideanSpace ℝ (Fin (p * q)) => Fentry (zeta p d z) k j (lam τ q x))
      = Fentry (zeta p d z) k j ∘ (lam τ q) from rfl,
    fderiv_comp x hdiff (lam τ q).differentiableAt, (lam τ q).fderiv]
  exact (ContinuousLinearMap.opNorm_comp_le _ _).trans
    (mul_le_mul_of_nonneg_right (norm_fderiv_Fentry_le hζ hp k j _) (norm_nonneg _))

theorem fderiv_Fsig_single (hζ : (zeta p d z).im ≠ 0) (τ : Fin p → ℝ) (k : Fin q) (j : Fin p)
    (x : EuclideanSpace ℝ (Fin (p * q))) :
    fderiv ℝ (fun x : EuclideanSpace ℝ (Fin (p * q)) => Fentry (zeta p d z) k j (lam τ q x)) x
        (EuclideanSpace.single (finProdFinEquiv (j, k)) (1 : ℝ))
      = ((τ j : ℝ) : ℂ) * fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q x)
          (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)) := by
  have hdiff := differentiableAt_Fentry hζ k j (lam τ q x)
  rw [show (fun x : EuclideanSpace ℝ (Fin (p * q)) => Fentry (zeta p d z) k j (lam τ q x))
      = Fentry (zeta p d z) k j ∘ (lam τ q) from rfl,
    fderiv_comp x hdiff (lam τ q).differentiableAt, (lam τ q).fderiv,
    ContinuousLinearMap.comp_apply, lam_single, map_smul, Complex.real_smul]

/-- **Stein's identity for the entry `(j, k)`.** -/
theorem stein_entry (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p)
    (k : Fin q) :
    ∫ B, ((B j k : ℝ) : ℂ) * Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B))
        ∂gaussianMatrix p q
      = ∫ B, ((τ j : ℝ) : ℂ) * fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q (matrixEquivE p q B))
          (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)) ∂gaussianMatrix p q := by
  have hζ := zeta_im_pos hz hp hd
  rw [integral_entry_mul_gaussianMatrix_complex (contDiff_Fsig hζ.ne' τ k j)
    (norm_fderiv_Fsig_le hζ hp τ k j) j k]
  exact integral_congr_ae (Filter.Eventually.of_forall fun B => fderiv_Fsig_single hζ.ne' τ k j _)

private theorem integrable_coord' (r : Fin (p * q)) :
    Integrable (fun x : EuclideanSpace ℝ (Fin (p * q)) => (x r : ℝ))
      (GaussianMeasure.stdGaussianE (p * q)) := by
  rw [stdGaussianE_eq]
  have h := IsGaussian.integrable_inner_mul_of_norm_fderiv_le
    (μ := stdGaussian (EuclideanSpace ℝ (Fin (p * q))))
    (F := fun _ : EuclideanSpace ℝ (Fin (p * q)) => (1 : ℝ)) (L := 0)
    contDiff_const (fun x => by simp) (EuclideanSpace.single r (1 : ℝ))
  simp only [mul_one] at h
  have hinner : ∀ x : EuclideanSpace ℝ (Fin (p * q)),
      inner ℝ (EuclideanSpace.single r (1 : ℝ) : EuclideanSpace ℝ (Fin (p * q))) x = x r := by
    intro x
    rw [EuclideanSpace.inner_single_left]
    simp
  simpa [hinner] using h

private theorem integrable_entry' (j : Fin p) (k : Fin q) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => ((B j k : ℝ) : ℂ)) (gaussianMatrix p q) := by
  have h1 : Integrable
      (fun x : EuclideanSpace ℝ (Fin (p * q)) => ((x (finProdFinEquiv (j, k)) : ℝ) : ℂ))
      (GaussianMeasure.stdGaussianE (p * q)) := (integrable_coord' _).ofReal
  have h2 := (measurePreserving_matrixEquivE p q).integrable_comp_of_integrable h1
  simpa [Function.comp_def] using h2

theorem integrable_lhs (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p)
    (k : Fin q) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      ((B j k : ℝ) : ℂ) * Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B)))
      (gaussianMatrix p q) := by
  have hζ := zeta_im_pos hz hp hd
  have hcont : Continuous (fun x : EuclideanSpace ℝ (Fin (p * q)) =>
      Fentry (zeta p d z) k j (lam τ q x)) := (contDiff_Fsig hζ.ne' τ k j).continuous
  have hmeas : Measurable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B))) :=
    hcont.measurable.comp (matrixEquivE p q).measurable
  exact (integrable_entry' j k).mul_bdd (c := BP (zeta p d z) p) hmeas.aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => norm_Pmat_le hζ hp (lam τ q (matrixEquivE p q B)) k j)

theorem integrable_rhs (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p)
    (k : Fin q) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q (matrixEquivE p q B))
        (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ))) (gaussianMatrix p q) := by
  have hζ := zeta_im_pos hz hp hd
  have hcont : Continuous (fun y : EuclideanSpace ℝ (Fin (q * p)) =>
      fderiv ℝ (Fentry (p := q) (d := p) (zeta p d z) k j) y
        (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ))) :=
    ((contDiff_Fentry (n := 1) hζ.ne' k j).continuous_fderiv one_ne_zero).clm_apply
      continuous_const
  have hmeas := (hcont.comp (lam τ q).continuous).measurable.comp (matrixEquivE p q).measurable
  refine Integrable.mono'
    (integrable_const (Lconst (zeta p d z) q p
      * ‖EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)‖))
    (by simpa [Function.comp_def] using hmeas.aestronglyMeasurable)
    (Filter.Eventually.of_forall fun B => ?_)
  refine ((fderiv ℝ (Fentry (p := q) (d := p) (zeta p d z) k j)
    (lam τ q (matrixEquivE p q B))).le_opNorm _).trans ?_
  exact mul_le_mul_of_nonneg_right (norm_fderiv_Fentry_le hζ hp k j _) (norm_nonneg _)

/-- **The row identity in the `q × p` picture.** `E[p (1 + ζ G_{jj})]` equals `τ_j²` times
the expected divergence of the column `j`. -/
theorem stein_row_raw (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p) :
    ∫ B, (p : ℂ) * (1 + zeta p d z * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j)
        ∂gaussianMatrix p q
      = ((τ j : ℝ) : ℂ) ^ 2 * ∫ B,
          ((q : ℂ) * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
            - ((p : ℂ) + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))).trace)
                * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
            - (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
                + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))
                    * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))) j j))
          ∂gaussianMatrix p q := by
  have hζ := zeta_im_pos hz hp hd
  have hL := fun k : Fin q => integrable_lhs hz hp hd τ j k
  have hR := fun k : Fin q => integrable_rhs hz hp hd τ j k
  have e1 : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      (p : ℂ) * (1 + zeta p d z * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j)
        = ∑ k : Fin q, ((τ j : ℝ) : ℂ)
            * (((B j k : ℝ) : ℂ) * Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B))) := by
    intro B
    rw [← sum_Ymat_mul_Pmat_col hζ.ne' hp _ j]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Ymat_apply, lam_matrixEquivE, Matrix.transpose_apply, Ysig, Matrix.diagonal_mul,
      Complex.ofReal_mul, mul_assoc]
    rfl
  have e2 : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      ((q : ℂ) * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
        - ((p : ℂ) + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))).trace)
            * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
        - (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
            + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))
                * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))) j j))
      = ∑ k : Fin q, fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q (matrixEquivE p q B))
          (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)) :=
    fun B => (sum_fderiv_Fentry_single_col hζ.ne' hp _ j).symm
  calc ∫ B, (p : ℂ) * (1 + zeta p d z * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j)
        ∂gaussianMatrix p q
      = ∫ B, ∑ k : Fin q, ((τ j : ℝ) : ℂ)
            * (((B j k : ℝ) : ℂ) * Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B)))
          ∂gaussianMatrix p q := integral_congr_ae (Filter.Eventually.of_forall e1)
    _ = ∑ k : Fin q, ((τ j : ℝ) : ℂ) * ∫ B, ((B j k : ℝ) : ℂ)
            * Fentry (zeta p d z) k j (lam τ q (matrixEquivE p q B)) ∂gaussianMatrix p q := by
        rw [integral_finsetSum _ (fun k _ => (hL k).const_mul _)]
        exact Finset.sum_congr rfl fun k _ => integral_const_mul _ _
    _ = ∑ k : Fin q, ((τ j : ℝ) : ℂ) * ∫ B, ((τ j : ℝ) : ℂ)
            * fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q (matrixEquivE p q B))
              (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)) ∂gaussianMatrix p q := by
        exact Finset.sum_congr rfl fun k _ => by rw [stein_entry hz hp hd τ j k]
    _ = ((τ j : ℝ) : ℂ) ^ 2 * ∫ B, ∑ k : Fin q,
            fderiv ℝ (Fentry (zeta p d z) k j) (lam τ q (matrixEquivE p q B))
              (EuclideanSpace.single (finProdFinEquiv (k, j)) (1 : ℝ)) ∂gaussianMatrix p q := by
        rw [integral_finsetSum _ (fun k _ => hR k), Finset.mul_sum]
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [integral_const_mul]
        ring
    _ = _ := by
        rw [integral_congr_ae (Filter.Eventually.of_forall e2)]

end SteinRow

/-! ### Deterministic bounds: diagonal entries and the trace of a resolvent -/

section Bounds

variable {n : ℕ} {W : Matrix (Fin n) (Fin n) ℝ}

theorem sum_sq_eigU_row (hW : W.IsHermitian) (i : Fin n) : ∑ a, (eigU hW i a) ^ 2 = 1 := by
  have h := congrFun (congrFun (eigU_mul_transpose hW) i) i
  rw [Matrix.mul_apply, Matrix.one_apply_eq] at h
  rw [← h]
  exact Finset.sum_congr rfl fun a _ => by rw [Matrix.transpose_apply, sq]

/-- A diagonal entry of `U diag(f) Uᵀ` is bounded by `sup ‖f‖`. -/
theorem norm_conj_diag_le (hW : W.IsHermitian) (f : Fin n → ℂ) {C : ℝ}
    (hf : ∀ a, ‖f a‖ ≤ C) (i : Fin n) :
    ‖(R4C.cmat (eigU hW) * Matrix.diagonal f * (R4C.cmat (eigU hW))ᵀ) i i‖ ≤ C := by
  have hentry : (R4C.cmat (eigU hW) * Matrix.diagonal f * (R4C.cmat (eigU hW))ᵀ) i i
      = ∑ a, ((eigU hW i a : ℝ) : ℂ) * f a * ((eigU hW i a : ℝ) : ℂ) := by
    rw [Matrix.mul_apply]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Matrix.mul_diagonal, Matrix.transpose_apply]
    rfl
  rw [hentry]
  calc ‖∑ a, ((eigU hW i a : ℝ) : ℂ) * f a * ((eigU hW i a : ℝ) : ℂ)‖
      ≤ ∑ a, ‖((eigU hW i a : ℝ) : ℂ) * f a * ((eigU hW i a : ℝ) : ℂ)‖ := norm_sum_le _ _
    _ ≤ ∑ a, (eigU hW i a) ^ 2 * C := by
        refine Finset.sum_le_sum fun a _ => ?_
        rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
          show |eigU hW i a| * ‖f a‖ * |eigU hW i a| = (eigU hW i a) ^ 2 * ‖f a‖ by
            rw [← sq_abs]; ring]
        exact mul_le_mul_of_nonneg_left (hf a) (sq_nonneg _)
    _ = C := by rw [← Finset.sum_mul, sum_sq_eigU_row hW i, one_mul]

theorem norm_resolvC_diag_le (hW : W.IsHermitian) (hz : 0 < z.im) (i : Fin n) :
    ‖R4C.resolvC W z i i‖ ≤ 1 / z.im := by
  rw [R4C.resolvC_eq_conj hW hz.ne']
  refine norm_conj_diag_le hW _ (fun a => ?_) i
  rw [one_div]
  exact R4C.norm_inv_eigenvalue_sub_le hW hz a

theorem norm_resolvC_sq_diag_le (hW : W.IsHermitian) (hz : 0 < z.im) (i : Fin n) :
    ‖(R4C.resolvC W z * R4C.resolvC W z) i i‖ ≤ 1 / z.im ^ 2 := by
  rw [R4C.resolvC_mul_resolvC_eq_conj hW hz.ne']
  refine norm_conj_diag_le hW _ (fun a => ?_) i
  rw [norm_pow, one_div, ← inv_pow]
  exact pow_le_pow_left₀ (norm_nonneg _) (R4C.norm_inv_eigenvalue_sub_le hW hz a) 2

theorem norm_trace_resolvC_le (hW : W.IsHermitian) (hz : 0 < z.im) :
    ‖(R4C.resolvC W z).trace‖ ≤ (n : ℝ) / z.im := by
  rw [R4C.trace_resolvC hW hz.ne']
  calc ‖∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹‖
      ≤ ∑ _a : Fin n, (z.im)⁻¹ :=
        (norm_sum_le _ _).trans
          (Finset.sum_le_sum fun a _ => R4C.norm_inv_eigenvalue_sub_le hW hz a)
    _ = (n : ℝ) / z.im := by simp [div_eq_mul_inv]

end Bounds

/-! ### The model side: bounds, measurability, integrability -/

section ModelBounds

theorem norm_gsig_le (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ)
    (j : Fin p) : ‖gsig τ d z B j‖ ≤ 1 / z.im :=
  norm_resolvC_diag_le (isHermitian_Wsig τ d B) hz j

theorem norm_Gsig_sq_diag_le (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ)
    (B : Matrix (Fin p) (Fin q) ℝ) (j : Fin p) :
    ‖(Gsig τ d z B * Gsig τ d z B) j j‖ ≤ 1 / z.im ^ 2 :=
  norm_resolvC_sq_diag_le (isHermitian_Wsig τ d B) hz j

theorem norm_ssig_le (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    ‖ssig τ d z B‖ ≤ (q : ℝ) / ((d : ℝ) * z.im) := by
  rw [ssig, norm_mul, norm_inv, Complex.norm_natCast]
  have h := norm_trace_resolvC_le (isHermitian_Wsig' τ d B) hz
  calc ((d : ℝ))⁻¹ * ‖(R4C.resolvC (Wsig' τ d B) z).trace‖
      ≤ ((d : ℝ))⁻¹ * ((q : ℝ) / z.im) := mul_le_mul_of_nonneg_left h (by positivity)
    _ = (q : ℝ) / ((d : ℝ) * z.im) := by rw [div_eq_mul_inv, div_eq_mul_inv, mul_inv]; ring

/-- `Gmat ζ (lam τ B) = κ⁻¹ G_σ`, the translation read from the `q × p` side. -/
theorem Gmat_zeta_eq (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) = ((kap p d : ℝ) : ℂ)⁻¹ • Gsig τ d z B := by
  have hkC : ((kap p d : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (kap_pos hp hd).ne'
  rw [zeta, Gsig_eq_smul_Gmat hz hp hd, smul_smul, inv_mul_cancel₀ hkC, one_smul]

theorem measurable_Gsig_entry (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (i l : Fin p) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => Gsig τ d z B i l) := by
  have hζ := zeta_im_pos hz hp hd
  have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => Gsig τ d z B i l)
      = fun B => ((kap p d : ℝ) : ℂ)
          * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) i l := by
    funext B
    rw [Gmat_zeta_eq hz.ne' hp hd, Matrix.smul_apply, smul_eq_mul, ← mul_assoc,
      mul_inv_cancel₀ (Complex.ofReal_ne_zero.mpr (kap_pos hp hd).ne'), one_mul]
  rw [heq]
  exact measurable_const.mul
    (((contDiff_Gmat_entry (n := 0) hζ.ne' i l).continuous.comp
      (lam τ q).continuous).measurable.comp (matrixEquivE p q).measurable)

theorem measurable_gsig (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => gsig τ d z B j) :=
  measurable_Gsig_entry hz hp hd τ j j

theorem measurable_Gsig_sq_diag (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (j : Fin p) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => (Gsig τ d z B * Gsig τ d z B) j j) := by
  have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => (Gsig τ d z B * Gsig τ d z B) j j)
      = fun B => ∑ l, Gsig τ d z B j l * Gsig τ d z B l j := by
    funext B
    rw [Matrix.mul_apply]
  rw [heq]
  exact Finset.measurable_sum _ fun l _ =>
    (measurable_Gsig_entry hz hp hd τ j l).mul (measurable_Gsig_entry hz hp hd τ l j)

theorem measurable_ssig (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => ssig τ d z B) := by
  have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => ssig τ d z B)
      = fun B => (d : ℂ)⁻¹ * ∑ j, Gsig τ d z B j j - ((q : ℂ) - p) / ((d : ℂ) * z) := by
    funext B
    rw [trace_identity hz.ne' hd, Matrix.trace]
    rfl
  rw [heq]
  exact (measurable_const.mul (Finset.measurable_sum _ fun j _ =>
    measurable_Gsig_entry hz hp hd τ j j)).sub measurable_const

theorem integrable_gsig (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => gsig τ d z B j) (gaussianMatrix p q) :=
  Integrable.mono' (integrable_const (1 / z.im))
    (measurable_gsig hz hp hd τ j).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => norm_gsig_le hz τ d B j)

theorem integrable_ssig_mul_gsig (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (j : Fin p) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => ssig τ d z B * gsig τ d z B j)
      (gaussianMatrix p q) := by
  refine Integrable.mono' (integrable_const ((q : ℝ) / ((d : ℝ) * z.im) * (1 / z.im)))
    ((measurable_ssig hz hp hd τ).mul (measurable_gsig hz hp hd τ j)).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => ?_)
  rw [norm_mul]
  exact mul_le_mul (norm_ssig_le hz τ d B) (norm_gsig_le hz τ d B j) (norm_nonneg _)
    (by positivity)

theorem integrable_gsig_add (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (j : Fin p) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j) (gaussianMatrix p q) := by
  refine Integrable.mono' (integrable_const (1 / z.im + ‖z‖ * (1 / z.im ^ 2)))
    ((measurable_gsig hz hp hd τ j).add
      (measurable_const.mul (measurable_Gsig_sq_diag hz hp hd τ j))).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => ?_)
  refine (norm_add_le _ _).trans (add_le_add (norm_gsig_le hz τ d B j) ?_)
  rw [norm_mul]
  exact mul_le_mul_of_nonneg_left (norm_Gsig_sq_diag_le hz τ d B j) (norm_nonneg z)

/-- Pointwise: the left side of `stein_row_raw` in model terms. -/
theorem lhs_translate (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) (j : Fin p) :
    (p : ℂ) * (1 + zeta p d z * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j)
      = (p : ℂ) * (1 + z * gsig τ d z B j) := by
  have hkC : ((kap p d : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (kap_pos hp hd).ne'
  rw [Gmat_zeta_eq hz hp hd, Matrix.smul_apply, smul_eq_mul, zeta, gsig]
  congr 2
  field_simp

/-- Pointwise: the divergence of `stein_row_raw` in model terms, after the trace identity. -/
theorem rhs_translate (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) (j : Fin p) :
    ((q : ℂ) * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
      - ((p : ℂ) + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))).trace)
          * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
      - (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B)) j j
          + zeta p d z * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))
              * Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))) j j))
      = (p : ℂ) * (-(z * (ssig τ d z B * gsig τ d z B j))
          - (d : ℂ)⁻¹ * (gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j)) := by
  have hz0 : z ≠ 0 := fun h => hz (by rw [h]; simp)
  have hpC : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hkap : ((kap p d : ℝ) : ℂ) = (d : ℂ) / (p : ℂ) := by
    rw [kap]
    push_cast
    ring
  rw [Gmat_zeta_eq hz hp hd, trace_identity hz hd, Matrix.smul_mul, Matrix.mul_smul, smul_smul,
    Matrix.smul_apply, Matrix.smul_apply, Matrix.trace_smul, smul_eq_mul, smul_eq_mul,
    smul_eq_mul, zeta, hkap, gsig]
  field_simp
  ring

end ModelBounds

/-! ### The row identity and the block identity (plan section 3.4) -/

section RowBlock

/-- The residual `r_j = -τ_j² d⁻¹ E[g_j + z (G_σ²)_{jj}]` of `stein_row`. `q` is explicit
because the law `gaussianMatrix p q` is. -/
noncomputable def steinErr (τ : Fin p → ℝ) (q d : ℕ) (z : ℂ) (j : Fin p) : ℂ :=
  -(((τ j : ℝ) : ℂ) ^ 2 * (d : ℂ)⁻¹)
    * ∫ B, (gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j) ∂gaussianMatrix p q

/-- **(H5, `stein_row`).** `E[1 + z g_j] = -τ_j² z E[s g_j] + r_j`. -/
theorem stein_row (hz : 0 < z.im) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p) :
    ∫ B, (1 + z * gsig τ d z B j) ∂gaussianMatrix p q
      = -(((τ j : ℝ) : ℂ) ^ 2 * z) * ∫ B, ssig τ d z B * gsig τ d z B j ∂gaussianMatrix p q
        + steinErr τ q d z j := by
  have hp : 0 < p := Fin.pos j
  have hpC : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hraw := stein_row_raw (q := q) hz hp hd τ j
  rw [integral_congr_ae (Filter.Eventually.of_forall fun B => lhs_translate hz.ne' hp hd τ B j),
    integral_congr_ae (Filter.Eventually.of_forall fun B => rhs_translate hz.ne' hp hd τ B j),
    integral_const_mul, integral_const_mul,
    integral_sub (f := fun B : Matrix (Fin p) (Fin q) ℝ => -(z * (ssig τ d z B * gsig τ d z B j)))
      (g := fun B : Matrix (Fin p) (Fin q) ℝ =>
        (d : ℂ)⁻¹ * (gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j))
      (by exact ((integrable_ssig_mul_gsig hz hp hd τ j).const_mul z).neg)
      ((integrable_gsig_add hz hp hd τ j).const_mul _),
    integral_neg, integral_const_mul, integral_const_mul] at hraw
  rw [steinErr]
  apply mul_left_cancel₀ hpC
  linear_combination hraw

/-- **The residual bound.** `‖r_j‖ ≤ τ_j² d⁻¹ (η⁻¹ + ‖z‖ η⁻²)`. -/
theorem norm_steinErr_le (hz : 0 < z.im) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p) :
    ‖steinErr τ q d z j‖ ≤ (τ j) ^ 2 / d * (1 / z.im + ‖z‖ / z.im ^ 2) := by
  have hpt : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      ‖gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j‖ ≤ 1 / z.im + ‖z‖ / z.im ^ 2 := by
    intro B
    refine (norm_add_le _ _).trans (add_le_add (norm_gsig_le hz τ d B j) ?_)
    rw [norm_mul, div_eq_mul_one_div]
    exact mul_le_mul_of_nonneg_left (norm_Gsig_sq_diag_le hz τ d B j) (norm_nonneg z)
  have hb : ‖∫ B, (gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j) ∂gaussianMatrix p q‖
      ≤ 1 / z.im + ‖z‖ / z.im ^ 2 := by
    have h := norm_integral_le_of_norm_le_const (μ := gaussianMatrix p q)
      (Filter.Eventually.of_forall hpt)
    simpa using h
  rw [steinErr, norm_mul, norm_neg, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
    sq_abs, norm_inv, Complex.norm_natCast]
  calc (τ j) ^ 2 * ((d : ℝ))⁻¹
        * ‖∫ B, (gsig τ d z B j + z * (Gsig τ d z B * Gsig τ d z B) j j) ∂gaussianMatrix p q‖
      ≤ (τ j) ^ 2 * ((d : ℝ))⁻¹ * (1 / z.im + ‖z‖ / z.im ^ 2) :=
        mul_le_mul_of_nonneg_left hb (by positivity)
    _ = (τ j) ^ 2 / d * (1 / z.im + ‖z‖ / z.im ^ 2) := by ring

/-- The average of `g_j` over a set of rows `J` (a block in H6). -/
noncomputable def gsigAvg (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) : ℂ :=
  (J.card : ℂ)⁻¹ * ∑ j ∈ J, gsig τ d z B j

/-- The averaged residual. -/
noncomputable def steinErrAvg (τ : Fin p → ℝ) (q d : ℕ) (z : ℂ) (J : Finset (Fin p)) : ℂ :=
  (J.card : ℂ)⁻¹ * ∑ j ∈ J, steinErr τ q d z j

/-- **(H5, `stein_block`).** For rows of equal variance `τ_j² = σ₀`,
`E[1 + z g_J] = -σ₀ z E[s g_J] + r_J` with `g_J` the row average over `J`. -/
theorem stein_block (hz : 0 < z.im) (hd : 0 < d) (τ : Fin p → ℝ) {J : Finset (Fin p)}
    (hJ : J.Nonempty) {σ₀ : ℝ} (hσ : ∀ j ∈ J, τ j ^ 2 = σ₀) :
    ∫ B, (1 + z * gsigAvg τ d z J B) ∂gaussianMatrix p q
      = -((σ₀ : ℂ) * z) * ∫ B, ssig τ d z B * gsigAvg τ d z J B ∂gaussianMatrix p q
        + steinErrAvg τ q d z J := by
  obtain ⟨j₀, hj₀⟩ := hJ
  have hp : 0 < p := Fin.pos j₀
  have hcard : (J.card : ℂ) ≠ 0 := by
    exact_mod_cast (Finset.card_pos.mpr ⟨j₀, hj₀⟩).ne'
  have e1 : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      1 + z * gsigAvg τ d z J B = (J.card : ℂ)⁻¹ * ∑ j ∈ J, (1 + z * gsig τ d z B j) := by
    intro B
    rw [gsigAvg, Finset.sum_add_distrib, Finset.sum_const, nsmul_eq_mul, mul_one,
      ← Finset.mul_sum, mul_add, inv_mul_cancel₀ hcard]
    ring
  have e2 : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      ssig τ d z B * gsigAvg τ d z J B
        = (J.card : ℂ)⁻¹ * ∑ j ∈ J, ssig τ d z B * gsig τ d z B j := by
    intro B
    rw [gsigAvg, mul_left_comm, Finset.mul_sum]
  rw [integral_congr_ae (Filter.Eventually.of_forall e1),
    integral_congr_ae (Filter.Eventually.of_forall e2), integral_const_mul, integral_const_mul,
    integral_finsetSum (f := fun j (B : Matrix (Fin p) (Fin q) ℝ) => 1 + z * gsig τ d z B j) J
      (fun j _ => by
        exact (integrable_const (1 : ℂ)).add ((integrable_gsig hz hp hd τ j).const_mul z)),
    integral_finsetSum (f := fun j (B : Matrix (Fin p) (Fin q) ℝ) => ssig τ d z B * gsig τ d z B j)
      J (fun j _ => integrable_ssig_mul_gsig hz hp hd τ j), steinErrAvg]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j hj => ?_
  rw [stein_row hz hd τ j, ← hσ j hj]
  push_cast
  ring

/-- The averaged residual bound: `‖r_J‖ ≤ σ₀ d⁻¹ (η⁻¹ + ‖z‖ η⁻²)`. -/
theorem norm_steinErrAvg_le (hz : 0 < z.im) (hd : 0 < d) (τ : Fin p → ℝ) {J : Finset (Fin p)}
    (hJ : J.Nonempty) {σ₀ : ℝ} (hσ : ∀ j ∈ J, τ j ^ 2 = σ₀) :
    ‖steinErrAvg τ q d z J‖ ≤ σ₀ / d * (1 / z.im + ‖z‖ / z.im ^ 2) := by
  have hcardR : (J.card : ℝ) ≠ 0 := by
    exact_mod_cast (Finset.card_pos.mpr hJ).ne'
  rw [steinErrAvg, norm_mul, norm_inv, Complex.norm_natCast]
  calc ((J.card : ℝ))⁻¹ * ‖∑ j ∈ J, steinErr τ q d z j‖
      ≤ ((J.card : ℝ))⁻¹ * ∑ j ∈ J, ‖steinErr τ q d z j‖ :=
        mul_le_mul_of_nonneg_left (norm_sum_le _ _) (by positivity)
    _ ≤ ((J.card : ℝ))⁻¹ * ∑ _j ∈ J, σ₀ / d * (1 / z.im + ‖z‖ / z.im ^ 2) := by
        refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun j hj => ?_) (by positivity)
        rw [← hσ j hj]
        exact norm_steinErr_le hz hd τ j
    _ = σ₀ / d * (1 / z.im + ‖z‖ / z.im ^ 2) := by
        rw [Finset.sum_const, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hcardR, one_mul]

end RowBlock

/-! ### The derivative of a diagonal entry of `G`, generic picture

`∂G_{ii}(h) = -(2/n) ∑_{k,a} P_{ki} G_{ai} h_{ka}`, so `‖∇ G_{ii}‖ ≤ (2/n) ‖P e_i‖ ‖G e_i‖`
with `‖G e_i‖² ≤ η⁻²` and `‖P e_i‖² = n ∑_c λ_c |g_c|² U_{ic}² ≤ n (η + ‖z‖) η⁻²`. -/

section DiagDeriv

variable {m n : ℕ}

/-- `Zᵀ Z = n diag(λ)`; the `hZZ` step of `SteinStep.sum_sq_Zmat'`, exported. -/
theorem transpose_Zmat_mul_Zmat (hn : 0 < n) (x : EuclideanSpace ℝ (Fin (m * n))) :
    (Zmat x)ᵀ * Zmat x = (n : ℝ) • Matrix.diagonal (eigVal x) := by
  set Y := (matrixEquivE m n).symm x with hY
  set U := eigUx x with hU
  have hn0 : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  have hYY : Yᵀ * Y = (n : ℝ) • gram Y := by
    rw [gram, smul_smul, mul_inv_cancel₀ hn0, one_smul]
  have hUW : Uᵀ * gram Y * U = Matrix.diagonal (eigVal x) := by
    have hc := eigU_conj (isHermitian_gram Y)
    calc Uᵀ * gram Y * U
        = Uᵀ * (U * Matrix.diagonal (eigVal x) * Uᵀ) * U := by
          rw [hU, eigUx, hY]
          rw [show Matrix.diagonal (eigVal x)
              = Matrix.diagonal (isHermitian_gram ((matrixEquivE m n).symm x)).eigenvalues from
            rfl]
          rw [hc]
      _ = (Uᵀ * U) * Matrix.diagonal (eigVal x) * (Uᵀ * U) := by
          simp only [Matrix.mul_assoc]
      _ = Matrix.diagonal (eigVal x) := by
          rw [hU, eigUx, transpose_eigU_mul, Matrix.one_mul, Matrix.mul_one]
  rw [Zmat, Matrix.transpose_mul, ← hU, ← hY]
  calc Uᵀ * Yᵀ * (Y * U) = Uᵀ * (Yᵀ * Y) * U := by simp only [Matrix.mul_assoc]
    _ = Uᵀ * ((n : ℝ) • gram Y) * U := by rw [hYY]
    _ = (n : ℝ) • (Uᵀ * gram Y * U) := by rw [Matrix.mul_smul, Matrix.smul_mul]
    _ = (n : ℝ) • Matrix.diagonal (eigVal x) := by rw [hUW]

/-- `∑_k ‖∑_c w_c Z_{kc}‖² = n ∑_c λ_c ‖w_c‖²`. -/
theorem sum_normSq_Zmat_mulVec (hn : 0 < n) (x : EuclideanSpace ℝ (Fin (m * n)))
    (w : Fin n → ℂ) :
    ∑ k : Fin m, ‖∑ c : Fin n, w c * ((Zmat x k c : ℝ) : ℂ)‖ ^ 2
      = ∑ c : Fin n, (n : ℝ) * eigVal x c * ‖w c‖ ^ 2 := by
  set Z := Zmat x with hZ
  set u : Fin n → ℝ := fun c => (w c).re with hu
  set v : Fin n → ℝ := fun c => (w c).im with hv
  have hre : ∀ k, (∑ c : Fin n, w c * ((Z k c : ℝ) : ℂ)).re = (Z *ᵥ u) k := by
    intro k
    rw [Complex.re_sum, Matrix.mulVec, dotProduct]
    exact Finset.sum_congr rfl fun c _ => by simp [hu, Complex.mul_re, mul_comm]
  have him : ∀ k, (∑ c : Fin n, w c * ((Z k c : ℝ) : ℂ)).im = (Z *ᵥ v) k := by
    intro k
    rw [Complex.im_sum, Matrix.mulVec, dotProduct]
    exact Finset.sum_congr rfl fun c _ => by simp [hv, Complex.mul_im, mul_comm]
  have hnorm : ∀ k, ‖∑ c : Fin n, w c * ((Z k c : ℝ) : ℂ)‖ ^ 2
      = ((Z *ᵥ u) k) ^ 2 + ((Z *ᵥ v) k) ^ 2 := by
    intro k
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
    ring
  have hsq : ∀ t : Fin n → ℝ, ∑ k, ((Z *ᵥ t) k) ^ 2 = ∑ c, (n : ℝ) * eigVal x c * (t c) ^ 2 := by
    intro t
    have h1 : ∑ k, ((Z *ᵥ t) k) ^ 2 = (Z *ᵥ t) ⬝ᵥ (Z *ᵥ t) := by
      rw [dotProduct]
      exact Finset.sum_congr rfl fun k _ => sq _
    rw [h1, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec,
      hZ, transpose_Zmat_mul_Zmat hn, Matrix.smul_mulVec, dotProduct]
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [Pi.smul_apply, Matrix.mulVec_diagonal, smul_eq_mul]
    ring
  rw [Finset.sum_congr rfl (fun k (_ : k ∈ Finset.univ) => hnorm k), Finset.sum_add_distrib,
    hsq u, hsq v, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [Complex.sq_norm, Complex.normSq_apply, hu, hv]
  ring

/-- `∑_k |P_{ki}|² ≤ n (η + ‖z‖) / η²`. -/
theorem sum_normSq_Pmat_col_le (hz : 0 < z.im) (hn : 0 < n) (x : EuclideanSpace ℝ (Fin (m * n)))
    (i : Fin n) :
    ∑ k : Fin m, ‖Pmat z x k i‖ ^ 2 ≤ (n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) := by
  have hP : ∀ k, Pmat z x k i
      = ∑ c : Fin n, (gval z x c * ((eigUx x i c : ℝ) : ℂ)) * ((Zmat x k c : ℝ) : ℂ) := by
    intro k
    rw [Pmat_entry_eq hz.ne' x k i]
    exact Finset.sum_congr rfl fun c _ => by ring
  rw [Finset.sum_congr rfl fun k _ => by rw [hP k], sum_normSq_Zmat_mulVec hn]
  have hrow : ∑ c : Fin n, (eigUx x i c) ^ 2 = 1 := sum_sq_eigU_row (isHermitian_gram _) i
  calc ∑ c : Fin n, (n : ℝ) * eigVal x c * ‖gval z x c * ((eigUx x i c : ℝ) : ℂ)‖ ^ 2
      ≤ ∑ c : Fin n, (n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) * (eigUx x i c) ^ 2 := by
        refine Finset.sum_le_sum fun c _ => ?_
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs]
        have h := eigVal_mul_norm_gval_sq_le hz hn x c
        have hn' : (0 : ℝ) ≤ n := Nat.cast_nonneg n
        calc (n : ℝ) * eigVal x c * (‖gval z x c‖ ^ 2 * (eigUx x i c) ^ 2)
            = (n : ℝ) * (eigVal x c * ‖gval z x c‖ ^ 2) * (eigUx x i c) ^ 2 := by ring
          _ ≤ (n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) * (eigUx x i c) ^ 2 :=
              mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left h hn') (sq_nonneg _)
    _ = (n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) := by rw [← Finset.mul_sum, hrow, mul_one]

/-- `∑_a |G_{ai}|² ≤ η⁻²`. -/
theorem sum_normSq_Gmat_col_le (hz : 0 < z.im) (x : EuclideanSpace ℝ (Fin (m * n)))
    (i : Fin n) :
    ∑ a : Fin n, ‖Gmat z x a i‖ ^ 2 ≤ 1 / z.im ^ 2 := by
  have hG : ∀ a, Gmat z x a i
      = ∑ c : Fin n, (gval z x c * ((eigUx x i c : ℝ) : ℂ)) * ((eigUx x a c : ℝ) : ℂ) := by
    intro a
    rw [Gmat_entry_eq hz.ne' x a i]
    exact Finset.sum_congr rfl fun c _ => by ring
  rw [Finset.sum_congr rfl fun a _ => by rw [hG a]]
  -- orthogonal invariance, as in `ResolvDeriv.sum_normSq_mulVec` (private there)
  have hinv : ∑ a : Fin n, ‖∑ c : Fin n, (gval z x c * ((eigUx x i c : ℝ) : ℂ))
      * ((eigUx x a c : ℝ) : ℂ)‖ ^ 2 = ∑ c : Fin n, ‖gval z x c * ((eigUx x i c : ℝ) : ℂ)‖ ^ 2 := by
    set U := eigUx x with hU
    set w : Fin n → ℂ := fun c => gval z x c * ((eigUx x i c : ℝ) : ℂ) with hw
    set u : Fin n → ℝ := fun c => (w c).re with hu
    set v : Fin n → ℝ := fun c => (w c).im with hv
    have hre : ∀ a, (∑ c : Fin n, w c * ((U a c : ℝ) : ℂ)).re = (U *ᵥ u) a := by
      intro a
      rw [Complex.re_sum, Matrix.mulVec, dotProduct]
      exact Finset.sum_congr rfl fun c _ => by simp [hu, Complex.mul_re, mul_comm]
    have him : ∀ a, (∑ c : Fin n, w c * ((U a c : ℝ) : ℂ)).im = (U *ᵥ v) a := by
      intro a
      rw [Complex.im_sum, Matrix.mulVec, dotProduct]
      exact Finset.sum_congr rfl fun c _ => by simp [hv, Complex.mul_im, mul_comm]
    have hnorm : ∀ a, ‖∑ c : Fin n, w c * ((U a c : ℝ) : ℂ)‖ ^ 2
        = ((U *ᵥ u) a) ^ 2 + ((U *ᵥ v) a) ^ 2 := by
      intro a
      rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
      ring
    have hsq : ∀ t : Fin n → ℝ, ∑ a, ((U *ᵥ t) a) ^ 2 = ∑ c, (t c) ^ 2 := by
      intro t
      have h1 : ∑ a, ((U *ᵥ t) a) ^ 2 = (U *ᵥ t) ⬝ᵥ (U *ᵥ t) := by
        rw [dotProduct]
        exact Finset.sum_congr rfl fun a _ => sq _
      rw [h1, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec,
        hU, transpose_eigUx_mul, Matrix.one_mulVec, dotProduct]
      exact Finset.sum_congr rfl fun c _ => (sq _).symm
    rw [Finset.sum_congr rfl (fun a (_ : a ∈ Finset.univ) => hnorm a), Finset.sum_add_distrib,
      hsq u, hsq v, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [Complex.sq_norm, Complex.normSq_apply, hu, hv]
    ring
  rw [hinv]
  have hrow : ∑ c : Fin n, (eigUx x i c) ^ 2 = 1 := sum_sq_eigU_row (isHermitian_gram _) i
  calc ∑ c : Fin n, ‖gval z x c * ((eigUx x i c : ℝ) : ℂ)‖ ^ 2
      ≤ ∑ c : Fin n, 1 / z.im ^ 2 * (eigUx x i c) ^ 2 := by
        refine Finset.sum_le_sum fun c _ => ?_
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs]
        refine mul_le_mul_of_nonneg_right ?_ (sq_nonneg _)
        rw [one_div, ← inv_pow]
        exact pow_le_pow_left₀ (norm_nonneg _) (norm_gval_le hz x c) 2
    _ = 1 / z.im ^ 2 := by rw [← Finset.mul_sum, hrow, mul_one]

/-- **The gradient of a diagonal entry.** `∂G_{ii}(h) = -(2/n) ∑_{k,a} P_{ki} G_{ai} h_{ka}`. -/
theorem fderiv_Gmat_diag_apply (hz : z.im ≠ 0) (x h : EuclideanSpace ℝ (Fin (m * n)))
    (i : Fin n) :
    fderiv ℝ (fun y : EuclideanSpace ℝ (Fin (m * n)) => Gmat z y i i) x h
      = -(2 / (n : ℂ)) * ∑ k : Fin m, Pmat z x k i
          * ∑ a : Fin n, ((h (finProdFinEquiv (k, a)) : ℝ) : ℂ) * Gmat z x a i := by
  rw [fderiv_Gmat_entry hz i i x h]
  have hdW : dWmat m n x h
      = (n : ℂ)⁻¹ • ((Ymat x)ᵀ * Ymat h) + (n : ℂ)⁻¹ • ((Ymat h)ᵀ * Ymat x) := rfl
  set G := Gmat z x with hG
  have hGt : Gᵀ = G := transpose_Gmat hz x
  -- the second half has the same diagonal as the first (transpose)
  have hT : (G * ((Ymat h)ᵀ * Ymat x) * G) i i = (G * ((Ymat x)ᵀ * Ymat h) * G) i i := by
    have : (G * ((Ymat h)ᵀ * Ymat x) * G)ᵀ = G * ((Ymat x)ᵀ * Ymat h) * G := by
      simp only [Matrix.transpose_mul, Matrix.transpose_transpose, hGt, Matrix.mul_assoc]
    rw [← this, Matrix.transpose_apply]
  have hP : G * (Ymat x)ᵀ = (Pmat z x)ᵀ := by
    rw [Pmat, Matrix.transpose_mul, hGt]
  have hentry : (G * ((Ymat x)ᵀ * Ymat h) * G) i i
      = ∑ k : Fin m, Pmat z x k i * ∑ a : Fin n, ((h (finProdFinEquiv (k, a)) : ℝ) : ℂ)
          * Gmat z x a i := by
    rw [← Matrix.mul_assoc, hP, Matrix.mul_assoc, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Matrix.transpose_apply, Matrix.mul_apply]
    congr 1
  rw [hdW, Matrix.mul_add, Matrix.add_mul, Matrix.add_apply, Matrix.mul_smul, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.smul_mul, Matrix.smul_apply, Matrix.smul_apply, smul_eq_mul,
    smul_eq_mul, hT, hentry]
  ring

/-- **The gradient bound of a diagonal entry.** `‖∇ G_{ii}‖ ≤ √(4 (η + ‖z‖) / (n η⁴))`. -/
theorem norm_fderiv_Gmat_diag_le (hz : 0 < z.im) (hn : 0 < n)
    (x : EuclideanSpace ℝ (Fin (m * n))) (i : Fin n) :
    ‖fderiv ℝ (fun y : EuclideanSpace ℝ (Fin (m * n)) => Gmat z y i i) x‖
      ≤ Real.sqrt (4 * (z.im + ‖z‖) / ((n : ℝ) * z.im ^ 4)) := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hnn : 0 ≤ z.im + ‖z‖ := by positivity
  refine ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg _) fun h => ?_
  rw [fderiv_Gmat_diag_apply hz.ne' x h i, norm_mul, norm_neg, norm_div, Complex.norm_natCast,
    Complex.norm_ofNat]
  set A : ℝ := ∑ k : Fin m, ‖Pmat z x k i‖ ^ 2 with hA
  set Bv : ℝ := ∑ a : Fin n, ‖Gmat z x a i‖ ^ 2 with hB
  have hA' := sum_normSq_Pmat_col_le hz hn x i
  have hB' := sum_normSq_Gmat_col_le hz x i
  -- Cauchy-Schwarz over the pairs `(k, a)`
  have hcs : ‖∑ k : Fin m, Pmat z x k i
      * ∑ a : Fin n, ((h (finProdFinEquiv (k, a)) : ℝ) : ℂ) * Gmat z x a i‖
      ≤ Real.sqrt (A * Bv) * ‖h‖ := by
    have h1 : ‖∑ k : Fin m, Pmat z x k i
        * ∑ a : Fin n, ((h (finProdFinEquiv (k, a)) : ℝ) : ℂ) * Gmat z x a i‖
        ≤ ∑ k : Fin m, ∑ a : Fin n,
            (‖Pmat z x k i‖ * ‖Gmat z x a i‖) * |h (finProdFinEquiv (k, a))| := by
      refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
      rw [norm_mul]
      refine (mul_le_mul_of_nonneg_left (norm_sum_le _ _) (norm_nonneg _)).trans ?_
      rw [Finset.mul_sum]
      refine Finset.sum_le_sum fun a _ => le_of_eq ?_
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
      ring
    have h2 : ∑ k : Fin m, ∑ a : Fin n,
        (‖Pmat z x k i‖ * ‖Gmat z x a i‖) * |h (finProdFinEquiv (k, a))|
        ≤ Real.sqrt (A * Bv) * ‖h‖ := by
      rw [← Fintype.sum_prod_type']
      have hh : ‖h‖ ^ 2 = ∑ r : Fin m × Fin n, |h (finProdFinEquiv r)| ^ 2 := by
        rw [EuclideanSpace.real_norm_sq_eq,
          ← Equiv.sum_comp finProdFinEquiv fun r => (h r) ^ 2]
        exact Finset.sum_congr rfl fun r _ => (sq_abs _).symm
      have hAB : A * Bv = ∑ r : Fin m × Fin n, (‖Pmat z x r.1 i‖ * ‖Gmat z x r.2 i‖) ^ 2 := by
        rw [hA, hB, Finset.sum_mul_sum, ← Fintype.sum_prod_type']
        exact Finset.sum_congr rfl fun r _ => by ring
      have hcs' := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin m × Fin n))
        (fun r => ‖Pmat z x r.1 i‖ * ‖Gmat z x r.2 i‖) (fun r => |h (finProdFinEquiv r)|)
      have hL0 : 0 ≤ ∑ r : Fin m × Fin n,
          (‖Pmat z x r.1 i‖ * ‖Gmat z x r.2 i‖) * |h (finProdFinEquiv r)| :=
        Finset.sum_nonneg fun r _ => by positivity
      have hR0 : 0 ≤ Real.sqrt (A * Bv) * ‖h‖ := by positivity
      refine (pow_le_pow_iff_left₀ hL0 hR0 two_ne_zero).mp ?_
      rw [mul_pow, Real.sq_sqrt (by rw [hAB]; positivity), hh, hAB]
      exact hcs'
    exact h1.trans h2
  calc 2 / (n : ℝ) * ‖∑ k : Fin m, Pmat z x k i
        * ∑ a : Fin n, ((h (finProdFinEquiv (k, a)) : ℝ) : ℂ) * Gmat z x a i‖
      ≤ 2 / (n : ℝ) * (Real.sqrt (A * Bv) * ‖h‖) :=
        mul_le_mul_of_nonneg_left hcs (by positivity)
    _ ≤ Real.sqrt (4 * (z.im + ‖z‖) / ((n : ℝ) * z.im ^ 4)) * ‖h‖ := by
        rw [← mul_assoc]
        refine mul_le_mul_of_nonneg_right ?_ (norm_nonneg h)
        have hAB : A * Bv ≤ (n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) * (1 / z.im ^ 2) :=
          mul_le_mul hA' hB' (Finset.sum_nonneg fun a _ => by positivity) (by positivity)
        have hsq : 2 / (n : ℝ) * Real.sqrt (A * Bv)
            = Real.sqrt ((2 / (n : ℝ)) ^ 2 * (A * Bv)) := by
          rw [Real.sqrt_mul (show (0 : ℝ) ≤ (2 / (n : ℝ)) ^ 2 by positivity) (A * Bv),
            Real.sqrt_sq (show (0 : ℝ) ≤ 2 / (n : ℝ) by positivity)]
        rw [hsq]
        refine Real.sqrt_le_sqrt ?_
        calc (2 / (n : ℝ)) ^ 2 * (A * Bv)
            ≤ (2 / (n : ℝ)) ^ 2 * ((n : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) * (1 / z.im ^ 2)) :=
              mul_le_mul_of_nonneg_left hAB (by positivity)
          _ = 4 * (z.im + ‖z‖) / ((n : ℝ) * z.im ^ 4) := by
              field_simp
              ring

end DiagDeriv

/-! ### Lipschitz constants and variance bounds (plan section 3.4, "Lipschitz concentration")

`Var (Re g_j) ≤ 4 lipG²` and `Var (Re s) ≤ 4 lipS²` from `R1.variance_le_of_lipschitz`, with
`lipG² = 4 S (η + ‖z‖) / (d η⁴)` and `lipS² = 4 S p (η + ‖z‖) / (d³ η⁴)` for any `S ≥ max τ_j²`.
The same for `Im`. -/

section Variance

/-- `g_j` on the flattened `B`: `κ · Gmat ζ (lam τ x)_{jj}`. -/
noncomputable def gfun (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (j : Fin p)
    (x : EuclideanSpace ℝ (Fin (p * q))) : ℂ :=
  ((kap p d : ℝ) : ℂ) * Gmat (zeta p d z) (lam τ q x) j j

theorem gfun_matrixEquivE (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) (j : Fin p)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    gfun τ d z j (matrixEquivE p q B) = gsig τ d z B j := by
  have hkC : ((kap p d : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (kap_pos hp hd).ne'
  rw [gfun, Gmat_zeta_eq hz hp hd, Matrix.smul_apply, smul_eq_mul, ← mul_assoc,
    mul_inv_cancel₀ hkC, one_mul]
  rfl

theorem contDiff_gfun (hζ : (zeta p d z).im ≠ 0) (τ : Fin p → ℝ) (j : Fin p) :
    ContDiff ℝ 1 (gfun (q := q) τ d z j) :=
  contDiff_const.mul ((contDiff_Gmat_entry hζ j j).comp (lam τ q).contDiff)

/-- The Lipschitz constant of `Re g_j`: `lipG² = 4 S (η + ‖z‖) / (d η⁴)`. -/
noncomputable def lipG (z : ℂ) (d : ℕ) (S : ℝ) : ℝ :=
  Real.sqrt (4 * S * (z.im + ‖z‖) / ((d : ℝ) * z.im ^ 4))

/-- The Lipschitz constant of `Re s`: `lipS² = 4 S p (η + ‖z‖) / (d³ η⁴)`. -/
noncomputable def lipS (z : ℂ) (d p : ℕ) (S : ℝ) : ℝ :=
  Real.sqrt (4 * S * p * (z.im + ‖z‖) / ((d : ℝ) ^ 3 * z.im ^ 4))

theorem lipG_sq (hz : 0 < z.im) {S : ℝ} (hS0 : 0 ≤ S) :
    lipG z d S ^ 2 = 4 * S * (z.im + ‖z‖) / ((d : ℝ) * z.im ^ 4) :=
  Real.sq_sqrt (by have := norm_nonneg z; positivity)

theorem lipS_sq (hz : 0 < z.im) {S : ℝ} (hS0 : 0 ≤ S) :
    lipS z d p S ^ 2 = 4 * S * p * (z.im + ‖z‖) / ((d : ℝ) ^ 3 * z.im ^ 4) :=
  Real.sq_sqrt (by have := norm_nonneg z; positivity)

theorem lipG_pos (hz : 0 < z.im) (hd : 0 < d) {S : ℝ} (hS0 : 0 < S) : 0 < lipG z d S := by
  have : (0 : ℝ) < d := by exact_mod_cast hd
  have := norm_nonneg z
  exact Real.sqrt_pos.mpr (by positivity)

theorem lipS_pos (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) {S : ℝ} (hS0 : 0 < S) :
    0 < lipS z d p S := by
  have : (0 : ℝ) < d := by exact_mod_cast hd
  have : (0 : ℝ) < p := by exact_mod_cast hp
  have := norm_nonneg z
  exact Real.sqrt_pos.mpr (by positivity)

theorem norm_fderiv_gfun_le (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) {S : ℝ}
    (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) (j : Fin p) (x : EuclideanSpace ℝ (Fin (p * q))) :
    ‖fderiv ℝ (gfun τ d z j) x‖ ≤ lipG z d S := by
  have hζ := zeta_im_pos hz hp hd
  have hκ := kap_pos hp hd
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hzn := norm_nonneg z
  have hG := differentiableAt_Gmat_entry hζ.ne' j j (lam τ q x)
  have hdiff : DifferentiableAt ℝ
      (fun y : EuclideanSpace ℝ (Fin (p * q)) => Gmat (zeta p d z) (lam τ q y) j j) x :=
    hG.comp x (lam τ q).differentiableAt
  have h1 : fderiv ℝ (gfun τ d z j) x
      = ((kap p d : ℝ) : ℂ) • ((fderiv ℝ (fun y => Gmat (zeta p d z) y j j) (lam τ q x)).comp
          (lam τ q)) := by
    rw [show gfun (q := q) τ d z j
        = fun y => ((kap p d : ℝ) : ℂ) * Gmat (zeta p d z) (lam τ q y) j j from rfl,
      fderiv_const_mul hdiff,
      show (fun y : EuclideanSpace ℝ (Fin (p * q)) => Gmat (zeta p d z) (lam τ q y) j j)
        = (fun y => Gmat (zeta p d z) y j j) ∘ (lam τ q) from rfl,
      fderiv_comp x hG (lam τ q).differentiableAt, (lam τ q).fderiv]
  rw [h1, norm_smul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hκ]
  have h2 := norm_fderiv_Gmat_diag_le hζ hp (lam τ q x) j
  have h3 := opNorm_lam_le (q := q) τ hS hS0
  set X : ℝ := 4 * ((zeta p d z).im + ‖zeta p d z‖) / ((p : ℝ) * (zeta p d z).im ^ 4) with hX
  have hX0 : 0 ≤ X := by
    have := norm_nonneg (zeta p d z)
    positivity
  calc kap p d * ‖(fderiv ℝ (fun y => Gmat (zeta p d z) y j j) (lam τ q x)).comp (lam τ q)‖
      ≤ kap p d * (‖fderiv ℝ (fun y => Gmat (zeta p d z) y j j) (lam τ q x)‖ * ‖lam τ q‖) :=
        mul_le_mul_of_nonneg_left (ContinuousLinearMap.opNorm_comp_le _ _) hκ.le
    _ ≤ kap p d * (Real.sqrt X * Real.sqrt S) :=
        mul_le_mul_of_nonneg_left
          (mul_le_mul h2 h3 (norm_nonneg _) (Real.sqrt_nonneg _)) hκ.le
    _ = Real.sqrt ((kap p d) ^ 2 * (X * S)) := by
        rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hκ.le, Real.sqrt_mul hX0]
    _ = lipG z d S := by
        rw [lipG]
        congr 1
        rw [hX, zeta, im_kap_mul, norm_kap_mul hp hd, kap]
        field_simp

theorem lipschitzWith_gfun_re (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) (j : Fin p) {C : ℝ≥0}
    (hC : lipG z d S ≤ (C : ℝ)) :
    LipschitzWith C (fun x : EuclideanSpace ℝ (Fin (p * q)) => (gfun τ d z j x).re) := by
  have hζ := zeta_im_pos hz hp hd
  have hF := contDiff_gfun (q := q) hζ.ne' τ j
  refine lipschitzWith_of_nnnorm_fderiv_le
    (fun x => (contDiff_re_comp hF).differentiable_one x) fun x => ?_
  rw [← NNReal.coe_le_coe, coe_nnnorm]
  exact (norm_fderiv_re_comp_le hF (norm_fderiv_gfun_le hz hp hd τ hS hS0 j) x).trans hC

theorem lipschitzWith_gfun_im (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) (j : Fin p) {C : ℝ≥0}
    (hC : lipG z d S ≤ (C : ℝ)) :
    LipschitzWith C (fun x : EuclideanSpace ℝ (Fin (p * q)) => (gfun τ d z j x).im) := by
  have hζ := zeta_im_pos hz hp hd
  have hF := contDiff_gfun (q := q) hζ.ne' τ j
  refine lipschitzWith_of_nnnorm_fderiv_le
    (fun x => (contDiff_im_comp hF).differentiable_one x) fun x => ?_
  rw [← NNReal.coe_le_coe, coe_nnnorm]
  exact (norm_fderiv_im_comp_le hF (norm_fderiv_gfun_le hz hp hd τ hS hS0 j) x).trans hC

/-- **(H5, variance of `g_j`, real part).** `Var (Re g_j) ≤ 4 lipG²`. -/
theorem variance_re_gsig_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 < S) (j : Fin p) :
    ∫ B, ((gsig τ d z B j).re - ∫ B', (gsig τ d z B' j).re ∂gaussianMatrix p q) ^ 2
        ∂gaussianMatrix p q
      ≤ 4 * lipG z d S ^ 2 := by
  have hL := lipG_pos hz hd hS0
  have h := variance_le_of_lipschitz (Nat.mul_pos hp hq) hL
    (lipschitzWith_gfun_re (q := q) hz hp hd τ hS hS0.le j (C := ⟨lipG z d S, hL.le⟩) le_rfl)
  simpa only [gfun_matrixEquivE hz.ne' hp hd] using h

/-- **(H5, variance of `g_j`, imaginary part).** -/
theorem variance_im_gsig_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 < S) (j : Fin p) :
    ∫ B, ((gsig τ d z B j).im - ∫ B', (gsig τ d z B' j).im ∂gaussianMatrix p q) ^ 2
        ∂gaussianMatrix p q
      ≤ 4 * lipG z d S ^ 2 := by
  have hL := lipG_pos hz hd hS0
  have h := variance_le_of_lipschitz (Nat.mul_pos hp hq) hL
    (lipschitzWith_gfun_im (q := q) hz hp hd τ hS hS0.le j (C := ⟨lipG z d S, hL.le⟩) le_rfl)
  simpa only [gfun_matrixEquivE hz.ne' hp hd] using h

/-- `s = sfun ζ (lam τ B) - (q - p)/(d z)`: the `d`-side trace is the `ResolvDeriv` trace
in the `q × p` picture, up to the constant of `trace_identity`. -/
theorem ssig_eq_sfun (hz : z.im ≠ 0) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    ssig τ d z B
      = sfun (zeta p d z) (lam τ q (matrixEquivE p q B)) - ((q : ℂ) - p) / ((d : ℂ) * z) := by
  have hpC : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hs : sfun (zeta p d z) (lam τ q (matrixEquivE p q B))
      = (p : ℂ)⁻¹ * (Gmat (zeta p d z) (lam τ q (matrixEquivE p q B))).trace := by
    rw [sfun_eq]
    rfl
  rw [trace_identity hz hd, hs, Gmat_zeta_eq hz hp hd, Matrix.trace_smul, smul_eq_mul]
  congr 1
  rw [kap]
  push_cast
  field_simp

/-- `lipC ζ p · √S ≤ lipS`, the constant of `ssig` after the chain rule through `lam`. -/
theorem lipC_zeta_mul_sqrt_le (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) {S : ℝ} (hS0 : 0 ≤ S) :
    lipC (zeta p d z) p * Real.sqrt S ≤ lipS z d p S := by
  have hζ := zeta_im_pos hz hp hd
  have hκ := kap_pos hp hd
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hzn := norm_nonneg z
  have hL0 : 0 ≤ lipC (zeta p d z) p * Real.sqrt S :=
    mul_nonneg (lipC_nonneg hζ p) (Real.sqrt_nonneg S)
  refine le_of_eq ?_
  rw [lipS]
  have hsq : (lipC (zeta p d z) p * Real.sqrt S) ^ 2
      = 4 * S * p * (z.im + ‖z‖) / ((d : ℝ) ^ 3 * z.im ^ 4) := by
    have hnn : 0 ≤ (zeta p d z).im + ‖zeta p d z‖ := by
      have := norm_nonneg (zeta p d z)
      positivity
    rw [mul_pow, Real.sq_sqrt hS0, lipC, div_pow, mul_pow, Real.sq_sqrt hnn, zeta, im_kap_mul,
      norm_kap_mul hp hd, kap]
    field_simp
    ring
  rw [← hsq, Real.sqrt_sq hL0]

theorem lipschitzWith_sfun_lam_re (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) {C : ℝ≥0} (hC : lipS z d p S ≤ (C : ℝ)) :
    LipschitzWith C
      (fun x : EuclideanSpace ℝ (Fin (p * q)) => (sfun (zeta p d z) (lam τ q x)).re) := by
  have hζ := zeta_im_pos hz hp hd
  set L : ℝ≥0 := ⟨lipC (zeta p d z) p, lipC_nonneg hζ p⟩ with hLdef
  have hLc : (L : ℝ) = lipC (zeta p d z) p := rfl
  have h1 : LipschitzWith (L * ‖lam τ q‖₊)
      ((fun y : EuclideanSpace ℝ (Fin (q * p)) => (sfun (zeta p d z) y).re) ∘ (lam τ q)) :=
    (lipschitzWith_sfun_re hζ hp (C := L) le_rfl).comp (lam τ q).lipschitz
  refine h1.weaken ?_
  rw [← NNReal.coe_le_coe, NNReal.coe_mul, hLc, coe_nnnorm]
  refine le_trans ?_ hC
  refine le_trans ?_ (lipC_zeta_mul_sqrt_le hz hp hd hS0)
  exact mul_le_mul_of_nonneg_left (opNorm_lam_le τ hS hS0) (lipC_nonneg hζ p)

theorem lipschitzWith_sfun_lam_im (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) {C : ℝ≥0} (hC : lipS z d p S ≤ (C : ℝ)) :
    LipschitzWith C
      (fun x : EuclideanSpace ℝ (Fin (p * q)) => (sfun (zeta p d z) (lam τ q x)).im) := by
  have hζ := zeta_im_pos hz hp hd
  set L : ℝ≥0 := ⟨lipC (zeta p d z) p, lipC_nonneg hζ p⟩ with hLdef
  have hLc : (L : ℝ) = lipC (zeta p d z) p := rfl
  have h1 : LipschitzWith (L * ‖lam τ q‖₊)
      ((fun y : EuclideanSpace ℝ (Fin (q * p)) => (sfun (zeta p d z) y).im) ∘ (lam τ q)) :=
    (lipschitzWith_sfun_im hζ hp (C := L) le_rfl).comp (lam τ q).lipschitz
  refine h1.weaken ?_
  rw [← NNReal.coe_le_coe, NNReal.coe_mul, hLc, coe_nnnorm]
  refine le_trans ?_ hC
  refine le_trans ?_ (lipC_zeta_mul_sqrt_le hz hp hd hS0)
  exact mul_le_mul_of_nonneg_left (opNorm_lam_le τ hS hS0) (lipC_nonneg hζ p)

theorem integrable_sfun_lam_re (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      (sfun (zeta p d z) (lam τ q (matrixEquivE p q B))).re) (gaussianMatrix p q) := by
  have hζ := zeta_im_pos hz hp hd
  have hcont : Continuous (fun x : EuclideanSpace ℝ (Fin (p * q)) =>
      sfun (zeta p d z) (lam τ q x)) :=
    (contDiff_sfun (n := 0) hζ.ne').continuous.comp (lam τ q).continuous
  refine Integrable.mono' (integrable_const (1 / (zeta p d z).im))
    (Complex.measurable_re.comp
      (hcont.measurable.comp (matrixEquivE p q).measurable)).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => ?_)
  rw [Real.norm_eq_abs]
  refine (Complex.abs_re_le_norm _).trans ?_
  rw [sfun_eq]
  exact R4C.norm_stieltjesC_le (isHermitian_gram _) hζ hp

theorem integrable_sfun_lam_im (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      (sfun (zeta p d z) (lam τ q (matrixEquivE p q B))).im) (gaussianMatrix p q) := by
  have hζ := zeta_im_pos hz hp hd
  have hcont : Continuous (fun x : EuclideanSpace ℝ (Fin (p * q)) =>
      sfun (zeta p d z) (lam τ q x)) :=
    (contDiff_sfun (n := 0) hζ.ne').continuous.comp (lam τ q).continuous
  refine Integrable.mono' (integrable_const (1 / (zeta p d z).im))
    (Complex.measurable_im.comp
      (hcont.measurable.comp (matrixEquivE p q).measurable)).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => ?_)
  rw [Real.norm_eq_abs]
  refine (Complex.abs_im_le_norm _).trans ?_
  rw [sfun_eq]
  exact R4C.norm_stieltjesC_le (isHermitian_gram _) hζ hp

/-- **(H5, variance of `s`, real part).** `Var (Re s) ≤ 4 lipS²`. -/
theorem variance_re_ssig_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 < S) :
    ∫ B, ((ssig τ d z B).re - ∫ B', (ssig τ d z B').re ∂gaussianMatrix p q) ^ 2
        ∂gaussianMatrix p q
      ≤ 4 * lipS z d p S ^ 2 := by
  have hL := lipS_pos hz hp hd hS0
  set F : EuclideanSpace ℝ (Fin (p * q)) → ℝ :=
    fun x => (sfun (zeta p d z) (lam τ q x)).re with hFdef
  set c : ℝ := (((q : ℂ) - p) / ((d : ℂ) * z)).re with hcdef
  have hvar := variance_le_of_lipschitz (Nat.mul_pos hp hq) hL
    (lipschitzWith_sfun_lam_re (q := q) hz hp hd τ hS hS0.le (C := ⟨lipS z d p S, hL.le⟩) le_rfl)
  have hpt : ∀ B : Matrix (Fin p) (Fin q) ℝ, (ssig τ d z B).re = F (matrixEquivE p q B) - c := by
    intro B
    rw [ssig_eq_sfun hz.ne' hp hd, Complex.sub_re]
  have hmean : ∫ B', (ssig τ d z B').re ∂gaussianMatrix p q
      = (∫ B', F (matrixEquivE p q B') ∂gaussianMatrix p q) - c := by
    rw [integral_congr_ae (Filter.Eventually.of_forall hpt),
      integral_sub (integrable_sfun_lam_re hz hp hd τ) (integrable_const c), integral_const]
    simp [hFdef]
  refine le_trans (le_of_eq ?_) hvar
  refine integral_congr_ae (Filter.Eventually.of_forall fun B => ?_)
  dsimp only
  rw [hpt, hmean]
  ring

/-- **(H5, variance of `s`, imaginary part).** -/
theorem variance_im_ssig_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 < S) :
    ∫ B, ((ssig τ d z B).im - ∫ B', (ssig τ d z B').im ∂gaussianMatrix p q) ^ 2
        ∂gaussianMatrix p q
      ≤ 4 * lipS z d p S ^ 2 := by
  have hL := lipS_pos hz hp hd hS0
  set F : EuclideanSpace ℝ (Fin (p * q)) → ℝ :=
    fun x => (sfun (zeta p d z) (lam τ q x)).im with hFdef
  set c : ℝ := (((q : ℂ) - p) / ((d : ℂ) * z)).im with hcdef
  have hvar := variance_le_of_lipschitz (Nat.mul_pos hp hq) hL
    (lipschitzWith_sfun_lam_im (q := q) hz hp hd τ hS hS0.le (C := ⟨lipS z d p S, hL.le⟩) le_rfl)
  have hpt : ∀ B : Matrix (Fin p) (Fin q) ℝ, (ssig τ d z B).im = F (matrixEquivE p q B) - c := by
    intro B
    rw [ssig_eq_sfun hz.ne' hp hd, Complex.sub_im]
  have hmean : ∫ B', (ssig τ d z B').im ∂gaussianMatrix p q
      = (∫ B', F (matrixEquivE p q B') ∂gaussianMatrix p q) - c := by
    rw [integral_congr_ae (Filter.Eventually.of_forall hpt),
      integral_sub (integrable_sfun_lam_im hz hp hd τ) (integrable_const c), integral_const]
    simp [hFdef]
  refine le_trans (le_of_eq ?_) hvar
  refine integral_congr_ae (Filter.Eventually.of_forall fun B => ?_)
  dsimp only
  rw [hpt, hmean]
  ring

end Variance

end HetStein

end StackedSVD
