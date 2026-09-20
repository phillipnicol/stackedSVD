/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R1het
import StackedSVD.RMT.Het.R5het
import StackedSVD.RMT.R2
import StackedSVD.RMT.T

/-!
# Item H7: the isotropic resolvent forms of the heteroscedastic block

Task H7 of `notes/archive/plan_heterolaw_A.md` (sections 2.4, 3.3 and the H7 row of section 4).
Setting of `RMT/Het/Stein.lean` and `RMT/Het/R1het.lean`: `B ~ gaussianMatrix p q`, row
scale `τ_j = w_{blk j}`, `W_σ = d⁻¹ (diag τ B)(diag τ B)ᵀ`, `G_σ(z) = (W_σ - z)⁻¹`. The six
forms of `ResolventLimitsHet` (`RMT/Het/R4het.lean`) are, for a deterministic vector `y`
with block norms `∑_{j ∈ J_i} y_j² = a_i` and a standard Gaussian vector `x` independent of
`B`, with `g = d^{-1/2} diag τ x`:

* `yᵀ G_σ y → ∑ a_i g_i(z)`, `yᵀ G_σ² y → ∑ a_i g_i'(z)` (block-isotropic means and
  Lipschitz concentration in the entries of `B`, plan section 3.3);
* `gᵀ G_σ g → ∑ c_i w_i² g_i(z)`, `gᵀ G_σ² g → ∑ c_i w_i² g_i'(z)` (conditional Chebyshev
  given `B`, R2 (c));
* `yᵀ G_σ g → 0`, `yᵀ G_σ² g → 0` (a linear form in `x`, Chebyshev).

Each limit is first proved at complex `z ∈ ℂ⁺` in the norm form and then transferred to
real `x > b` by `RMT/T.lean` (`T.Gen.tendstoInProb_cform_of_complex`), with the block limits
`g_i`, `g_i'` of item H6 on the boundary (`tendsto_gC_sGlob`, `tendsto_gCDeriv_sGlob`).

**Concentration route for the fixed vector.** Plan section 3.3 asks for Lipschitz
concentration. The Lipschitz constant comes from the resolvent identity
`G - G' = G (W' - W) G'` and the two bounds `‖G y‖ ≤ ‖y‖/η` and
`‖Bᵀ diag τ G y‖² = d · (Gy)ᴴ W (Gy) ≤ d ‖y‖² (1/η + ‖z‖/η²)`, with no derivative: the map
`B ↦ yᵀ G_σ(B) y` is `2 √S ‖y‖² √(1/η + ‖z‖/η²) / (η √d)`-Lipschitz for the Frobenius norm,
where `S = max_j τ_j²`. The mean is `∑_i a_i E[g_{J_i}]` by the sign and permutation
symmetry inside each block (the R2a argument of `RMT/R2.lean`, on the rows).

**The `b` convention.** The real-axis theorems take any `b ≥ bHet c w` with the edge event
at `b` as a hypothesis and give the limit at every `x > b`. Item H11 assembles
`ResolventLimitsHet m w b Phihet Psihet PhihetDeriv PsihetDeriv` from them, at
`b = bSF w c` (item H8, `bHet_le_bSF`) and at `b = bHet c w` (`HeteroEdge`).

STATUS: see `notes/archive/agent_reports/h7_r2het.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory Finset
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace HetR2

open HetStein HetR1 MPhet ResolvDeriv R4 R2

/-! ### Squared length of a complex vector and the basic bounds -/

section Vec

variable {n : ℕ}

/-- `nsq v = ∑ ‖v_k‖²`. -/
noncomputable def nsq (v : Fin n → ℂ) : ℝ := ∑ k, ‖v k‖ ^ 2

theorem nsq_nonneg (v : Fin n → ℂ) : 0 ≤ nsq v :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

theorem nsq_cvec (y : Fin n → ℝ) : nsq (R4C.cvec y) = y ⬝ᵥ y := by
  simp only [nsq, R4C.cvec, Complex.norm_real, Real.norm_eq_abs, dotProduct, sq,
    abs_mul_abs_self]

theorem nsq_smul (a : ℂ) (v : Fin n → ℂ) : nsq (a • v) = ‖a‖ ^ 2 * nsq v := by
  simp only [nsq, Pi.smul_apply, smul_eq_mul, norm_mul, mul_pow, Finset.mul_sum]

theorem nsq_mul_le {f : Fin n → ℂ} {C : ℝ} (hf : ∀ k, ‖f k‖ ≤ C) (v : Fin n → ℂ) :
    nsq (fun k => f k * v k) ≤ C ^ 2 * nsq v := by
  simp only [nsq, Finset.mul_sum, norm_mul, mul_pow]
  refine Finset.sum_le_sum fun k _ => ?_
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (norm_nonneg _) (hf k) 2) (sq_nonneg _)

/-- `((nsq v : ℝ) : ℂ) = star v ⬝ᵥ v`. -/
theorem ofReal_nsq (v : Fin n → ℂ) : ((nsq v : ℝ) : ℂ) = star v ⬝ᵥ v := by
  simp only [nsq, dotProduct, Pi.star_apply, Complex.ofReal_sum, ← Complex.normSq_eq_norm_sq,
    Complex.normSq_eq_conj_mul_self, Complex.star_def]

/-- Cauchy-Schwarz for the bilinear (unconjugated) dot product. -/
theorem norm_dotProduct_sq_le (a b : Fin n → ℂ) : ‖a ⬝ᵥ b‖ ^ 2 ≤ nsq a * nsq b := by
  have h1 : ‖a ⬝ᵥ b‖ ≤ ∑ k, ‖a k‖ * ‖b k‖ := by
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
    rw [norm_mul]
  have h2 : (∑ k, ‖a k‖ * ‖b k‖) ^ 2 ≤ (∑ k, ‖a k‖ ^ 2) * ∑ k, ‖b k‖ ^ 2 :=
    Finset.sum_mul_sq_le_sq_mul_sq _ _ _
  exact (pow_le_pow_left₀ (norm_nonneg _) h1 2).trans h2

theorem norm_dotProduct_le (a b : Fin n → ℂ) :
    ‖a ⬝ᵥ b‖ ≤ Real.sqrt (nsq a * nsq b) :=
  Real.le_sqrt_of_sq_le (norm_dotProduct_sq_le a b)

/-- Frobenius bound: `‖A v‖² ≤ ‖A‖_F² ‖v‖²`. -/
theorem nsq_cmapR_mulVec_le {m : ℕ} (A : Matrix (Fin m) (Fin n) ℝ) (v : Fin n → ℂ) :
    nsq (cmapR A *ᵥ v) ≤ (∑ i, ∑ j, A i j ^ 2) * nsq v := by
  unfold nsq
  rw [Finset.sum_mul]
  refine Finset.sum_le_sum fun i _ => ?_
  have h1 : ‖(cmapR A *ᵥ v) i‖ ≤ ∑ j, |A i j| * ‖v j‖ := by
    simp only [Matrix.mulVec, dotProduct, cmapR_apply]
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun j _ => ?_)
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
  have h2 : (∑ j, |A i j| * ‖v j‖) ^ 2 ≤ (∑ j, |A i j| ^ 2) * ∑ j, ‖v j‖ ^ 2 :=
    Finset.sum_mul_sq_le_sq_mul_sq _ _ _
  simp only [sq_abs] at h2
  exact (pow_le_pow_left₀ (norm_nonneg _) h1 2).trans h2

/-- Diagonal scaling: `‖diag τ v‖² ≤ S ‖v‖²` when `τ_j² ≤ S`. -/
theorem nsq_diagonal_mulVec_le (τ : Fin n → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S)
    (v : Fin n → ℂ) : nsq (cmapR (Matrix.diagonal τ) *ᵥ v) ≤ S * nsq v := by
  have hk : ∀ k, (cmapR (Matrix.diagonal τ) *ᵥ v) k = ((τ k : ℝ) : ℂ) * v k := by
    intro k
    simp only [Matrix.mulVec, dotProduct, cmapR_apply, Matrix.diagonal_apply]
    rw [Finset.sum_eq_single k (fun l _ hl => by simp [Ne.symm hl]) (by simp)]
    simp
  simp only [nsq, hk, norm_mul, mul_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs,
    Finset.mul_sum]
  exact Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_right (hS k) (sq_nonneg _)

theorem cmapR_diagonal_mulVec (τ : Fin n → ℝ) (v : Fin n → ℂ) :
    cmapR (Matrix.diagonal τ) *ᵥ v = fun k => ((τ k : ℝ) : ℂ) * v k := by
  funext k
  simp only [Matrix.mulVec, dotProduct, cmapR_apply, Matrix.diagonal_apply]
  rw [Finset.sum_eq_single k (fun l _ hl => by simp [Ne.symm hl]) (by simp)]
  simp

end Vec

/-! ### Resolvent bounds in the eigenbasis, and the Frobenius-Lipschitz estimate -/

section Resolvent

variable {n p q : ℕ} {W : Matrix (Fin n) (Fin n) ℝ} {z : ℂ}

theorem conjTranspose_cmapR {m k : ℕ} (A : Matrix (Fin m) (Fin k) ℝ) :
    (cmapR A)ᴴ = cmapR Aᵀ := by
  ext i j
  simp [cmapR, Matrix.conjTranspose_apply]

theorem cmapR_sub {m k : ℕ} (A B : Matrix (Fin m) (Fin k) ℝ) :
    cmapR (A - B) = cmapR A - cmapR B := by
  ext i j; simp [cmapR]

theorem cmapR_add {m k : ℕ} (A B : Matrix (Fin m) (Fin k) ℝ) :
    cmapR (A + B) = cmapR A + cmapR B := by
  ext i j; simp [cmapR]

theorem cmapR_smul {m k : ℕ} (c : ℝ) (A : Matrix (Fin m) (Fin k) ℝ) :
    cmapR (c • A) = (c : ℂ) • cmapR A := by
  ext i j; simp [cmapR]

theorem nsq_star (u : Fin n → ℂ) : nsq (star u) = nsq u := by simp [nsq]

/-- A real orthogonal matrix preserves `nsq` of complex vectors. -/
theorem nsq_cmat_mulVec {U : Matrix (Fin n) (Fin n) ℝ} (hU : Uᵀ * U = 1) (v : Fin n → ℂ) :
    nsq (R4C.cmat U *ᵥ v) = nsq v := by
  have hH : (R4C.cmat U)ᴴ = (R4C.cmat U)ᵀ := by
    ext i j
    simp [R4C.cmat, Matrix.conjTranspose_apply]
  apply Complex.ofReal_injective
  rw [ofReal_nsq, ofReal_nsq, Matrix.star_mulVec, hH, Matrix.vecMul_transpose,
    Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec,
    R2.transpose_cmat_mul hU, Matrix.one_mulVec]

/-- `‖U diag f Uᵀ v‖² ≤ C² ‖v‖²` when `‖f_a‖ ≤ C`. -/
theorem nsq_conj_mulVec_le (hW : W.IsHermitian) (f : Fin n → ℂ) {C : ℝ}
    (hf : ∀ a, ‖f a‖ ≤ C) (v : Fin n → ℂ) :
    nsq ((R4C.cmat (eigU hW) * Matrix.diagonal f * (R4C.cmat (eigU hW))ᵀ) *ᵥ v)
      ≤ C ^ 2 * nsq v := by
  rw [R4C.conj_mulVec_gen, nsq_cmat_mulVec (transpose_eigU_mul hW)]
  have h1 : nsq ((R4C.cmat (eigU hW))ᵀ *ᵥ v) = nsq v := by
    rw [← R4C.cmat_transpose]
    exact nsq_cmat_mulVec (R2.transpose_eigU_orth hW) v
  calc nsq (fun a => f a * ((R4C.cmat (eigU hW))ᵀ *ᵥ v) a)
      ≤ C ^ 2 * nsq ((R4C.cmat (eigU hW))ᵀ *ᵥ v) := nsq_mul_le hf _
    _ = C ^ 2 * nsq v := by rw [h1]

/-- `‖G(z) v‖² ≤ ‖v‖² / (Im z)²`. -/
theorem nsq_resolvC_mulVec_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin n → ℂ) :
    nsq (R4C.resolvC W z *ᵥ v) ≤ nsq v / z.im ^ 2 := by
  rw [R4C.resolvC_eq_conj hW hz.ne', div_eq_inv_mul, ← inv_pow]
  exact nsq_conj_mulVec_le hW _ (fun a => R4C.norm_inv_eigenvalue_sub_le hW hz a) v

/-- `‖G(z)² v‖² ≤ ‖v‖² / (Im z)⁴`. -/
theorem nsq_resolvC_sq_mulVec_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin n → ℂ) :
    nsq ((R4C.resolvC W z * R4C.resolvC W z) *ᵥ v) ≤ nsq v / z.im ^ 4 := by
  rw [R4C.resolvC_mul_resolvC_eq_conj hW hz.ne', div_eq_inv_mul,
    show z.im ^ 4 = (z.im ^ 2) ^ 2 by ring, ← inv_pow]
  refine nsq_conj_mulVec_le hW _ (fun a => ?_) v
  rw [norm_pow, ← inv_pow]
  exact pow_le_pow_left₀ (norm_nonneg _) (R4C.norm_inv_eigenvalue_sub_le hW hz a) 2

/-- `diag τ B Bᵀ diag τ = d W_σ`. -/
theorem diag_mul_gram_eq (τ : Fin p → ℝ) {d : ℕ} (hd : 0 < d) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.diagonal τ * B * Bᵀ * Matrix.diagonal τ = (d : ℝ) • Wsig τ d B := by
  rw [Wsig, Ysig, smul_smul, mul_inv_cancel₀ (Nat.cast_ne_zero.mpr hd.ne'), one_smul,
    Matrix.transpose_mul, Matrix.diagonal_transpose]
  simp only [Matrix.mul_assoc]

/-- **The energy bound.** `‖Bᵀ diag τ G_σ v‖² = d · (G v)ᴴ W_σ (G v) ≤ d ‖v‖² (1/η + ‖z‖/η²)`. -/
theorem nsq_transpose_diag_Gsig_le (hz : 0 < z.im) (τ : Fin p → ℝ) {d : ℕ} (hd : 0 < d)
    (B : Matrix (Fin p) (Fin q) ℝ) (v : Fin p → ℂ) :
    nsq (cmapR Bᵀ *ᵥ (cmapR (Matrix.diagonal τ) *ᵥ (Gsig τ d z B *ᵥ v)))
      ≤ d * (nsq v * (1 / z.im + ‖z‖ / z.im ^ 2)) := by
  have hW : (Wsig τ d B).IsHermitian := isHermitian_Wsig τ d B
  set u : Fin p → ℂ := Gsig τ d z B *ᵥ v with hudef
  set y := cmapR Bᵀ *ᵥ (cmapR (Matrix.diagonal τ) *ᵥ u) with hydef
  have hsy : star y = star u ᵥ* cmapR (Matrix.diagonal τ * B) := by
    rw [hydef, Matrix.star_mulVec, Matrix.star_mulVec, conjTranspose_cmapR, conjTranspose_cmapR,
      Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.vecMul_vecMul, ← cmapR_mul]
  have hyy : star y ⬝ᵥ y = (d : ℂ) * (star u ⬝ᵥ (R4C.cmat (Wsig τ d B) *ᵥ u)) := by
    rw [hsy, ← Matrix.dotProduct_mulVec, hydef, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
      ← cmapR_mul, ← cmapR_mul, Matrix.mul_assoc (Matrix.diagonal τ) B, ← Matrix.mul_assoc,
      diag_mul_gram_eq τ hd B, cmapR_smul, Matrix.smul_mulVec, dotProduct_smul,
      smul_eq_mul, cmat_eq_cmapR, Complex.ofReal_natCast]
  have hWu : R4C.cmat (Wsig τ d B) *ᵥ u = v + z • u := by
    have h : (R4C.cmat (Wsig τ d B) - z • (1 : Matrix (Fin p) (Fin p) ℂ)) *ᵥ u = v := by
      rw [hudef, Gsig, Matrix.mulVec_mulVec, cmat_sub_mul_resolvC hW hz.ne', Matrix.one_mulVec]
    rw [Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec] at h
    rw [← h]
    abel
  have hdot : star u ⬝ᵥ (R4C.cmat (Wsig τ d B) *ᵥ u) = star u ⬝ᵥ v + z * ((nsq u : ℝ) : ℂ) := by
    rw [hWu, dotProduct_add, dotProduct_smul, smul_eq_mul, ofReal_nsq]
  have hu : nsq u ≤ nsq v / z.im ^ 2 := nsq_resolvC_mulVec_le hW hz v
  have hv0 : 0 ≤ nsq v := nsq_nonneg v
  have hηpos : 0 < z.im := hz
  have h1 : ‖star u ⬝ᵥ v‖ ≤ nsq v / z.im := by
    refine (norm_dotProduct_le _ _).trans ?_
    rw [Real.sqrt_le_left (by positivity), nsq_star]
    calc nsq u * nsq v ≤ nsq v / z.im ^ 2 * nsq v :=
          mul_le_mul_of_nonneg_right hu hv0
      _ = (nsq v / z.im) ^ 2 := by field_simp
  have h2 : ‖z * ((nsq u : ℝ) : ℂ)‖ ≤ ‖z‖ * (nsq v / z.im ^ 2) := by
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (nsq_nonneg u)]
    exact mul_le_mul_of_nonneg_left hu (norm_nonneg z)
  have hny : nsq y = ‖star y ⬝ᵥ y‖ := by
    rw [← ofReal_nsq, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (nsq_nonneg y)]
  rw [hny, hyy, norm_mul, Complex.norm_natCast, hdot]
  refine mul_le_mul_of_nonneg_left ((norm_add_le _ _).trans ?_) (Nat.cast_nonneg d)
  calc ‖star u ⬝ᵥ v‖ + ‖z * ((nsq u : ℝ) : ℂ)‖
      ≤ nsq v / z.im + ‖z‖ * (nsq v / z.im ^ 2) := add_le_add h1 h2
    _ = nsq v * (1 / z.im + ‖z‖ / z.im ^ 2) := by ring

/-- `G - G' = G (W' - W) G'`. -/
theorem resolvC_sub_resolvC {W' : Matrix (Fin n) (Fin n) ℝ} (hW : W.IsHermitian)
    (hW' : W'.IsHermitian) (hz : z.im ≠ 0) :
    R4C.resolvC W z - R4C.resolvC W' z
      = R4C.resolvC W z * (R4C.cmat W' - R4C.cmat W) * R4C.resolvC W' z := by
  have h1 : R4C.cmat W' - R4C.cmat W
      = (R4C.cmat W' - z • (1 : Matrix (Fin n) (Fin n) ℂ)) - (R4C.cmat W - z • 1) := by abel
  rw [h1, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc, cmat_sub_mul_resolvC hW' hz,
    Matrix.mul_one, resolvC_mul_cmat_sub hW hz, Matrix.one_mul]

/-- `W_σ(B') - W_σ(B) = d⁻¹ diag τ ((B' - B) B'ᵀ + B (B' - B)ᵀ) diag τ`, over `ℂ`. -/
theorem cmat_Wsig_sub (τ : Fin p → ℝ) (d : ℕ) (B B' : Matrix (Fin p) (Fin q) ℝ) :
    R4C.cmat (Wsig τ d B') - R4C.cmat (Wsig τ d B)
      = (((d : ℝ)⁻¹ : ℝ) : ℂ) • (cmapR (Matrix.diagonal τ)
          * (cmapR (B' - B) * cmapR B'ᵀ + cmapR B * cmapR (B' - B)ᵀ)
          * cmapR (Matrix.diagonal τ)) := by
  have hreal : Wsig τ d B' - Wsig τ d B
      = (d : ℝ)⁻¹ • (Matrix.diagonal τ * ((B' - B) * B'ᵀ + B * (B' - B)ᵀ)
          * Matrix.diagonal τ) := by
    simp only [Wsig, Ysig, ← smul_sub]
    congr 1
    simp only [Matrix.transpose_mul, Matrix.diagonal_transpose, Matrix.sub_mul, Matrix.mul_sub,
      Matrix.transpose_sub, Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc]
    abel
  rw [cmat_eq_cmapR, cmat_eq_cmapR, ← cmapR_sub, hreal, cmapR_smul, cmapR_mul, cmapR_mul,
    cmapR_add, cmapR_mul, cmapR_mul]

theorem sum_sq_transpose {m k : ℕ} (A : Matrix (Fin m) (Fin k) ℝ) :
    ∑ i, ∑ j, Aᵀ i j ^ 2 = ∑ i, ∑ j, A i j ^ 2 := by
  rw [Finset.sum_comm]
  rfl

theorem transpose_cmapR_diagonal (τ : Fin p → ℝ) :
    (cmapR (Matrix.diagonal τ))ᵀ = cmapR (Matrix.diagonal τ) := by
  rw [← cmapR_transpose, Matrix.diagonal_transpose]

/-- The constant of the Frobenius-Lipschitz estimate: `2 √(S na nb (1/η + ‖z‖/η²) / η² / d)`. -/
noncomputable def lipC (z : ℂ) (S na nb : ℝ) (d : ℕ) : ℝ :=
  2 * Real.sqrt (S * na * nb * (1 / z.im + ‖z‖ / z.im ^ 2) / z.im ^ 2 / d)

/-- **The Frobenius-Lipschitz estimate for the bilinear form of `G_σ`.** -/
theorem norm_dotProduct_Gsig_sub_le (hz : 0 < z.im) (τ : Fin p → ℝ) {S : ℝ}
    (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) {d : ℕ} (hd : 0 < d)
    (B B' : Matrix (Fin p) (Fin q) ℝ) {a b : Fin p → ℂ} {na nb : ℝ} (ha : nsq a ≤ na)
    (hb : nsq b ≤ nb) :
    ‖a ⬝ᵥ ((Gsig τ d z B - Gsig τ d z B') *ᵥ b)‖
      ≤ lipC z S na nb d * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2) := by
  have hW : (Wsig τ d B).IsHermitian := isHermitian_Wsig τ d B
  have hW' : (Wsig τ d B').IsHermitian := isHermitian_Wsig τ d B'
  have hna0 : 0 ≤ na := (nsq_nonneg a).trans ha
  have hnb0 : 0 ≤ nb := (nsq_nonneg b).trans hb
  have hη : 0 < z.im := hz
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  set G := Gsig τ d z B with hGdef
  set G' := Gsig τ d z B' with hG'def
  set F := ∑ i, ∑ j, (B' i j - B i j) ^ 2 with hFdef
  set K := S * na * nb * (1 / z.im + ‖z‖ / z.im ^ 2) / z.im ^ 2 with hKdef
  have hF0 : 0 ≤ F := Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => sq_nonneg _
  have hK0 : 0 ≤ K := by positivity
  set u := cmapR (Matrix.diagonal τ) *ᵥ (G *ᵥ a) with hudef
  set u' := cmapR (Matrix.diagonal τ) *ᵥ (G' *ᵥ b) with hu'def
  have hdiff : G - G' = G * (R4C.cmat (Wsig τ d B') - R4C.cmat (Wsig τ d B)) * G' :=
    resolvC_sub_resolvC hW hW' hz.ne'
  have hform : a ⬝ᵥ ((G - G') *ᵥ b)
      = (((d : ℝ)⁻¹ : ℝ) : ℂ) * (u ⬝ᵥ ((cmapR (B' - B) * cmapR B'ᵀ
          + cmapR B * cmapR (B' - B)ᵀ) *ᵥ u')) := by
    rw [hdiff, cmat_Wsig_sub τ d B B']
    simp only [← Matrix.mulVec_mulVec, Matrix.smul_mulVec, Matrix.mulVec_smul,
      dotProduct_smul, smul_eq_mul]
    congr 1
    rw [Matrix.dotProduct_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
      ← Matrix.mulVec_transpose, transpose_cmapR_diagonal, hudef, hu'def, hGdef, Gsig,
      transpose_resolvC hW hz.ne']
  have hsplit : u ⬝ᵥ ((cmapR (B' - B) * cmapR B'ᵀ + cmapR B * cmapR (B' - B)ᵀ) *ᵥ u')
      = (cmapR (B' - B)ᵀ *ᵥ u) ⬝ᵥ (cmapR B'ᵀ *ᵥ u')
        + (cmapR Bᵀ *ᵥ u) ⬝ᵥ (cmapR (B' - B)ᵀ *ᵥ u') := by
    have hbil : ∀ (X : Matrix (Fin p) (Fin q) ℝ) (Y : Matrix (Fin q) (Fin p) ℝ) (w : Fin p → ℂ),
        u ⬝ᵥ ((cmapR X * cmapR Y) *ᵥ w) = (cmapR Xᵀ *ᵥ u) ⬝ᵥ (cmapR Y *ᵥ w) := by
      intro X Y w
      rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
        cmapR_transpose]
    rw [Matrix.add_mulVec, dotProduct_add, hbil, hbil]
  have hu : nsq u ≤ S * (na / z.im ^ 2) :=
    (nsq_diagonal_mulVec_le τ hS _).trans (mul_le_mul_of_nonneg_left
      ((nsq_resolvC_mulVec_le hW hz a).trans (div_le_div_of_nonneg_right ha (by positivity)))
      hS0)
  have hu' : nsq u' ≤ S * (nb / z.im ^ 2) :=
    (nsq_diagonal_mulVec_le τ hS _).trans (mul_le_mul_of_nonneg_left
      ((nsq_resolvC_mulVec_le hW' hz b).trans (div_le_div_of_nonneg_right hb (by positivity)))
      hS0)
  have hE : nsq (cmapR Bᵀ *ᵥ u) ≤ d * (na * (1 / z.im + ‖z‖ / z.im ^ 2)) :=
    (nsq_transpose_diag_Gsig_le hz τ hd B a).trans
      (mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right ha (by positivity)) hdR.le)
  have hE' : nsq (cmapR B'ᵀ *ᵥ u') ≤ d * (nb * (1 / z.im + ‖z‖ / z.im ^ 2)) :=
    (nsq_transpose_diag_Gsig_le hz τ hd B' b).trans
      (mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hb (by positivity)) hdR.le)
  have hΔu : nsq (cmapR (B' - B)ᵀ *ᵥ u) ≤ F * (S * (na / z.im ^ 2)) := by
    refine (nsq_cmapR_mulVec_le _ _).trans ?_
    rw [sum_sq_transpose]
    exact mul_le_mul_of_nonneg_left hu hF0
  have hΔu' : nsq (cmapR (B' - B)ᵀ *ᵥ u') ≤ F * (S * (nb / z.im ^ 2)) := by
    refine (nsq_cmapR_mulVec_le _ _).trans ?_
    rw [sum_sq_transpose]
    exact mul_le_mul_of_nonneg_left hu' hF0
  have ht1 : ‖(cmapR (B' - B)ᵀ *ᵥ u) ⬝ᵥ (cmapR B'ᵀ *ᵥ u')‖ ≤ Real.sqrt (F * d * K) := by
    refine (norm_dotProduct_le _ _).trans (Real.sqrt_le_sqrt ?_)
    calc nsq (cmapR (B' - B)ᵀ *ᵥ u) * nsq (cmapR B'ᵀ *ᵥ u')
        ≤ (F * (S * (na / z.im ^ 2))) * (d * (nb * (1 / z.im + ‖z‖ / z.im ^ 2))) :=
          mul_le_mul hΔu hE' (nsq_nonneg _) (by positivity)
      _ = F * d * K := by rw [hKdef]; ring
  have ht2 : ‖(cmapR Bᵀ *ᵥ u) ⬝ᵥ (cmapR (B' - B)ᵀ *ᵥ u')‖ ≤ Real.sqrt (F * d * K) := by
    refine (norm_dotProduct_le _ _).trans (Real.sqrt_le_sqrt ?_)
    calc nsq (cmapR Bᵀ *ᵥ u) * nsq (cmapR (B' - B)ᵀ *ᵥ u')
        ≤ (d * (na * (1 / z.im + ‖z‖ / z.im ^ 2))) * (F * (S * (nb / z.im ^ 2))) :=
          mul_le_mul hE hΔu' (nsq_nonneg _) (by positivity)
      _ = F * d * K := by rw [hKdef]; ring
  have hsq : Real.sqrt (F * d * K) = d * (Real.sqrt (K / d) * Real.sqrt F) := by
    rw [show F * d * K = (d * d) * ((K / d) * F) by field_simp, Real.sqrt_mul (by positivity),
      Real.sqrt_mul_self hdR.le, Real.sqrt_mul (by positivity)]
  rw [hform, hsplit, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (inv_nonneg.mpr hdR.le), lipC, ← hKdef]
  calc ((d : ℝ))⁻¹ * ‖(cmapR (B' - B)ᵀ *ᵥ u) ⬝ᵥ (cmapR B'ᵀ *ᵥ u')
        + (cmapR Bᵀ *ᵥ u) ⬝ᵥ (cmapR (B' - B)ᵀ *ᵥ u')‖
      ≤ ((d : ℝ))⁻¹ * (Real.sqrt (F * d * K) + Real.sqrt (F * d * K)) :=
        mul_le_mul_of_nonneg_left ((norm_add_le _ _).trans (add_le_add ht1 ht2))
          (inv_nonneg.mpr hdR.le)
    _ = 2 * Real.sqrt (K / d) * Real.sqrt F := by
        rw [hsq]
        field_simp
        ring

/-- The fixed-vector form `yᵀ G_σ y` is `lipC z S ‖y‖² ‖y‖² d`-Lipschitz in the entries. -/
theorem norm_qformC_Wsig_sub_le (hz : 0 < z.im) (τ : Fin p → ℝ) {S : ℝ}
    (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) {d : ℕ} (hd : 0 < d)
    (B B' : Matrix (Fin p) (Fin q) ℝ) (y : Fin p → ℝ) {ny : ℝ} (hy : y ⬝ᵥ y ≤ ny) :
    ‖R4C.qformC (Wsig τ d B) z y - R4C.qformC (Wsig τ d B') z y‖
      ≤ lipC z S ny ny d * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2) := by
  have h : R4C.qformC (Wsig τ d B) z y - R4C.qformC (Wsig τ d B') z y
      = R4C.cvec y ⬝ᵥ ((Gsig τ d z B - Gsig τ d z B') *ᵥ R4C.cvec y) := by
    rw [Matrix.sub_mulVec, dotProduct_sub]
    rfl
  rw [h]
  exact norm_dotProduct_Gsig_sub_le hz τ hS hS0 hd B B' (by rw [nsq_cvec]; exact hy)
    (by rw [nsq_cvec]; exact hy)

/-- The second-order fixed-vector form `yᵀ G_σ² y` is `2 lipC z S (‖y‖²/η²) ‖y‖² d`-Lipschitz. -/
theorem norm_qform2C_Wsig_sub_le (hz : 0 < z.im) (τ : Fin p → ℝ) {S : ℝ}
    (hS : ∀ j, τ j ^ 2 ≤ S) (hS0 : 0 ≤ S) {d : ℕ} (hd : 0 < d)
    (B B' : Matrix (Fin p) (Fin q) ℝ) (y : Fin p → ℝ) {ny : ℝ} (hy : y ⬝ᵥ y ≤ ny) :
    ‖R4C.qform2C (Wsig τ d B) z y - R4C.qform2C (Wsig τ d B') z y‖
      ≤ 2 * lipC z S (ny / z.im ^ 2) ny d * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2) := by
  have hW : (Wsig τ d B).IsHermitian := isHermitian_Wsig τ d B
  have hW' : (Wsig τ d B').IsHermitian := isHermitian_Wsig τ d B'
  set G := Gsig τ d z B with hGdef
  set G' := Gsig τ d z B' with hG'def
  have hsplit : G * G - G' * G' = G * (G - G') + (G - G') * G' := by
    rw [Matrix.mul_sub, Matrix.sub_mul]
    abel
  have h : R4C.qform2C (Wsig τ d B) z y - R4C.qform2C (Wsig τ d B') z y
      = (G *ᵥ R4C.cvec y) ⬝ᵥ ((G - G') *ᵥ R4C.cvec y)
        + R4C.cvec y ⬝ᵥ ((G - G') *ᵥ (G' *ᵥ R4C.cvec y)) := by
    change R4C.cvec y ⬝ᵥ ((G * G) *ᵥ R4C.cvec y) - R4C.cvec y ⬝ᵥ ((G' * G') *ᵥ R4C.cvec y) = _
    rw [← dotProduct_sub, ← Matrix.sub_mulVec, hsplit, Matrix.add_mulVec, dotProduct_add,
      ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec (v := R4C.cvec y)
      (A := G), ← Matrix.mulVec_transpose, hGdef, Gsig, transpose_resolvC hW hz.ne']
  have hny : nsq (R4C.cvec y) ≤ ny := by rw [nsq_cvec]; exact hy
  have hGy : nsq (G *ᵥ R4C.cvec y) ≤ ny / z.im ^ 2 :=
    (nsq_resolvC_mulVec_le hW hz _).trans (div_le_div_of_nonneg_right hny (by positivity))
  have hG'y : nsq (G' *ᵥ R4C.cvec y) ≤ ny / z.im ^ 2 :=
    (nsq_resolvC_mulVec_le hW' hz _).trans (div_le_div_of_nonneg_right hny (by positivity))
  have h1 := norm_dotProduct_Gsig_sub_le hz τ hS hS0 hd B B' hGy hny
  have h2 := norm_dotProduct_Gsig_sub_le hz τ hS hS0 hd B B' hny hG'y
  have hsymm : lipC z S ny (ny / z.im ^ 2) d = lipC z S (ny / z.im ^ 2) ny d := by
    simp only [lipC]
    ring_nf
  rw [hsymm] at h2
  rw [h]
  calc ‖(G *ᵥ R4C.cvec y) ⬝ᵥ ((G - G') *ᵥ R4C.cvec y)
        + R4C.cvec y ⬝ᵥ ((G - G') *ᵥ (G' *ᵥ R4C.cvec y))‖
      ≤ lipC z S (ny / z.im ^ 2) ny d * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2)
        + lipC z S (ny / z.im ^ 2) ny d * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2) :=
        (norm_add_le _ _).trans (add_le_add h1 h2)
    _ = _ := by ring

theorem lipC_nonneg (z : ℂ) (S na nb : ℝ) (d : ℕ) : 0 ≤ lipC z S na nb d := by
  unfold lipC; positivity

theorem lipC_pos {z : ℂ} (hz : 0 < z.im) {S na nb : ℝ} (hS : 0 < S) (hna : 0 < na)
    (hnb : 0 < nb) {d : ℕ} (hd : 0 < d) : 0 < lipC z S na nb d := by
  unfold lipC
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  positivity

theorem tendsto_lipC {dN : ℕ → ℕ} (hd : Tendsto dN atTop atTop) (z : ℂ) (S na nb : ℝ) :
    Tendsto (fun N => lipC z S na nb (dN N)) atTop (𝓝 0) := by
  have hdR : Tendsto (fun N => ((dN N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have h1 : Tendsto (fun N => S * na * nb * (1 / z.im + ‖z‖ / z.im ^ 2) / z.im ^ 2 / (dN N : ℝ))
      atTop (𝓝 0) := tendsto_const_nhds.div_atTop hdR
  have h2 := (Real.continuous_sqrt.tendsto 0).comp h1
  simp only [Function.comp_def, Real.sqrt_zero] at h2
  simpa [lipC] using h2.const_mul 2

/-- From a Frobenius-Lipschitz bound to `LipschitzWith` through `matrixEquivE`, real and
imaginary parts. -/
theorem lipschitzWith_of_norm_sub_le {g : Matrix (Fin p) (Fin q) ℝ → ℂ} {L : ℝ} (hL : 0 ≤ L)
    (h : ∀ B B', ‖g B - g B'‖ ≤ L * Real.sqrt (∑ i, ∑ j, (B' i j - B i j) ^ 2)) :
    LipschitzWith (Real.toNNReal L)
        (fun x : EuclideanSpace ℝ (Fin (p * q)) => (g ((matrixEquivE p q).symm x)).re) ∧
      LipschitzWith (Real.toNNReal L)
        (fun x : EuclideanSpace ℝ (Fin (p * q)) => (g ((matrixEquivE p q).symm x)).im) := by
  have key : ∀ x x' : EuclideanSpace ℝ (Fin (p * q)),
      ‖g ((matrixEquivE p q).symm x) - g ((matrixEquivE p q).symm x')‖ ≤ L * dist x x' := by
    intro x x'
    refine (h _ _).trans (le_of_eq ?_)
    congr 1
    rw [← Real.sqrt_sq dist_nonneg, dist_matrixEquivE_symm_sq]
    congr 1
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [← neg_sub, neg_sq]
  constructor
  · refine LipschitzWith.of_dist_le_mul fun x x' => ?_
    rw [Real.dist_eq, Real.coe_toNNReal L hL]
    refine le_trans ?_ (key x x')
    rw [← Complex.sub_re]
    exact Complex.abs_re_le_norm _
  · refine LipschitzWith.of_dist_le_mul fun x x' => ?_
    rw [Real.dist_eq, Real.coe_toNNReal L hL]
    refine le_trans ?_ (key x x')
    rw [← Complex.sub_im]
    exact Complex.abs_im_le_norm _

end Resolvent

/-! ### The mean of the fixed-vector forms: block symmetry (the R2a argument on the rows) -/

section Mean

variable {p q : ℕ} {z : ℂ}

/-- The row action `B ↦ Pᵀ B` of a signed permutation preserves the Gaussian law. -/
theorem measurePreserving_rowSgn (σ : Equiv.Perm (Fin p)) {ε : Fin p → ℝ}
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) :
    MeasurePreserving (fun B : Matrix (Fin p) (Fin q) ℝ => (sgnPermMat σ ε)ᵀ * B)
      (gaussianMatrix p q) (gaussianMatrix p q) :=
  measurePreserving_mul_left
    (by rw [Matrix.transpose_transpose]; exact mul_eq_one_comm.mp (sgnPermMat_orth hε)) q

/-- `diag τ` commutes with a signed permutation that preserves `τ`. -/
theorem diagonal_mul_sgnPermMat_transpose (τ : Fin p → ℝ) (σ : Equiv.Perm (Fin p))
    (ε : Fin p → ℝ) (hτ : ∀ j, τ (σ j) = τ j) :
    Matrix.diagonal τ * (sgnPermMat σ ε)ᵀ = (sgnPermMat σ ε)ᵀ * Matrix.diagonal τ := by
  ext j l
  simp only [Matrix.diagonal_mul, Matrix.mul_diagonal, Matrix.transpose_apply, sgnPermMat,
    Matrix.of_apply]
  split_ifs with h
  · rw [h, hτ]; ring
  · simp

theorem Wsig_rowSgn (τ : Fin p → ℝ) (d : ℕ) (σ : Equiv.Perm (Fin p)) (ε : Fin p → ℝ)
    (hτ : ∀ j, τ (σ j) = τ j) (B : Matrix (Fin p) (Fin q) ℝ) :
    Wsig τ d ((sgnPermMat σ ε)ᵀ * B)
      = (sgnPermMat σ ε)ᵀ * Wsig τ d B * sgnPermMat σ ε := by
  simp only [Wsig, Ysig, Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  have h1 : Matrix.diagonal τ * ((sgnPermMat σ ε)ᵀ * B)
      = (sgnPermMat σ ε)ᵀ * (Matrix.diagonal τ * B) := by
    rw [← Matrix.mul_assoc, diagonal_mul_sgnPermMat_transpose τ σ ε hτ, Matrix.mul_assoc]
  rw [h1, Matrix.transpose_mul, Matrix.transpose_transpose]
  simp only [Matrix.mul_assoc]

/-- Entries of `(cmat P)ᵀ K (cmat P)` for a signed permutation `P`. -/
theorem conj_sgnPermMat_apply (K : Matrix (Fin p) (Fin p) ℂ) (σ : Equiv.Perm (Fin p))
    (ε : Fin p → ℝ) (i j : Fin p) :
    ((R4C.cmat (sgnPermMat σ ε))ᵀ * K * R4C.cmat (sgnPermMat σ ε)) i j
      = (ε i : ℂ) * (ε j : ℂ) * K (σ i) (σ j) := by
  rw [Matrix.mul_apply]
  have hinner : ∀ l : Fin p,
      ((R4C.cmat (sgnPermMat σ ε))ᵀ * K) i l * R4C.cmat (sgnPermMat σ ε) l j
        = if l = σ j then (ε i : ℂ) * K (σ i) l * (ε j : ℂ) else 0 := by
    intro l
    rw [Matrix.mul_apply]
    have houter : ∀ k : Fin p,
        (R4C.cmat (sgnPermMat σ ε))ᵀ i k * K k l
          = if k = σ i then (ε i : ℂ) * K k l else 0 := by
      intro k
      change ((if k = σ i then ε i else 0 : ℝ) : ℂ) * _ = _
      by_cases h : k = σ i <;> simp [h]
    simp only [houter]
    rw [Finset.sum_ite_eq' Finset.univ (σ i) fun k => (ε i : ℂ) * K k l]
    change (if σ i ∈ Finset.univ then (ε i : ℂ) * K (σ i) l else 0) *
      ((if l = σ j then ε j else 0 : ℝ) : ℂ) = _
    by_cases h : l = σ j <;> simp [h]
  simp only [hinner]
  rw [Finset.sum_ite_eq' Finset.univ (σ j) fun l => (ε i : ℂ) * K (σ i) l * (ε j : ℂ)]
  simp only [Finset.mem_univ, if_true]
  ring

/-- A kernel `K(B)` that is conjugated by every `τ`-preserving signed row permutation, with
measurable entries bounded by `C`. `G_σ` and `G_σ²` are the two instances. -/
structure BlockSym (τ : Fin p → ℝ) (q : ℕ)
    (K : Matrix (Fin p) (Fin q) ℝ → Matrix (Fin p) (Fin p) ℂ) (C : ℝ) : Prop where
  conj : ∀ (σ : Equiv.Perm (Fin p)) (ε : Fin p → ℝ), (∀ j, ε j = 1 ∨ ε j = -1) →
    (∀ j, τ (σ j) = τ j) → ∀ B, K ((sgnPermMat σ ε)ᵀ * B)
      = (R4C.cmat (sgnPermMat σ ε))ᵀ * K B * R4C.cmat (sgnPermMat σ ε)
  meas : ∀ j k, Measurable fun B => K B j k
  bound : ∀ B j k, ‖K B j k‖ ≤ C

/-- Entrywise mean of a kernel. -/
noncomputable def meanK (K : Matrix (Fin p) (Fin q) ℝ → Matrix (Fin p) (Fin p) ℂ) (j k : Fin p) :
    ℂ :=
  ∫ B, K B j k ∂gaussianMatrix p q

variable {τ : Fin p → ℝ} {K : Matrix (Fin p) (Fin q) ℝ → Matrix (Fin p) (Fin p) ℂ} {C : ℝ}

theorem BlockSym.integrable_entry (h : BlockSym τ q K C) (j k : Fin p) :
    Integrable (fun B => K B j k) (gaussianMatrix p q) :=
  Integrable.of_bound (h.meas j k).aestronglyMeasurable C
    (Filter.Eventually.of_forall fun B => h.bound B j k)

theorem BlockSym.meanK_sgn (h : BlockSym τ q K C) (σ : Equiv.Perm (Fin p)) (ε : Fin p → ℝ)
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) (hτ : ∀ j, τ (σ j) = τ j) (j k : Fin p) :
    meanK K j k = (ε j : ℂ) * (ε k : ℂ) * meanK K (σ j) (σ k) := by
  have hmp := measurePreserving_rowSgn (q := q) σ hε
  have hchange : ∫ B, K B j k ∂(gaussianMatrix p q)
      = ∫ B, K ((sgnPermMat σ ε)ᵀ * B) j k ∂(gaussianMatrix p q) := by
    have h' := integral_map (φ := fun B : Matrix (Fin p) (Fin q) ℝ => (sgnPermMat σ ε)ᵀ * B)
      (μ := gaussianMatrix p q) hmp.measurable.aemeasurable (f := fun B => K B j k)
      (by rw [hmp.map_eq]; exact (h.meas j k).aestronglyMeasurable)
    rw [hmp.map_eq] at h'
    exact h'
  change ∫ B, K B j k ∂(gaussianMatrix p q) = _
  rw [hchange]
  simp only [h.conj σ ε hε hτ, conj_sgnPermMat_apply]
  rw [integral_const_mul]
  rfl

theorem BlockSym.meanK_offDiag (h : BlockSym τ q K C) {j k : Fin p} (hjk : j ≠ k) :
    meanK K j k = 0 := by
  set ε : Fin p → ℝ := fun l => if l = k then -1 else 1 with hεdef
  have hε : ∀ l, ε l = 1 ∨ ε l = -1 := by
    intro l
    by_cases hl : l = k <;> simp [hεdef, hl]
  have h' := h.meanK_sgn (Equiv.refl (Fin p)) ε hε (fun _ => rfl) j k
  simp only [hεdef, Equiv.refl_apply, if_neg hjk, if_pos rfl] at h'
  push_cast at h'
  linear_combination h' / 2

theorem BlockSym.meanK_diag_eq (h : BlockSym τ q K C) {j k : Fin p} (hjk : τ j = τ k) :
    meanK K j j = meanK K k k := by
  have hτ : ∀ l, τ (Equiv.swap j k l) = τ l := by
    intro l
    rcases eq_or_ne l j with rfl | h1
    · simp [hjk]
    rcases eq_or_ne l k with rfl | h2
    · simp [hjk]
    · simp [Equiv.swap_apply_of_ne_of_ne h1 h2]
  have h' := h.meanK_sgn (Equiv.swap j k) (fun _ => 1) (fun _ => Or.inl rfl) hτ j j
  simpa [Equiv.swap_apply_left] using h'

/-- **The mean of the quadratic form of a block-symmetric kernel.** -/
theorem BlockSym.integral_qform (h : BlockSym τ q K C) (y : Fin p → ℝ) :
    ∫ B, R4C.cvec y ⬝ᵥ (K B *ᵥ R4C.cvec y) ∂gaussianMatrix p q
      = ∑ j, ((y j ^ 2 : ℝ) : ℂ) * meanK K j j := by
  have hint : ∀ (i j : Fin p), Integrable
      (fun B : Matrix (Fin p) (Fin q) ℝ => (y i : ℂ) * (K B i j * (y j : ℂ)))
      (gaussianMatrix p q) :=
    fun i j => (((h.integrable_entry i j).mul_const (y j : ℂ)).const_mul (y i : ℂ))
  have hexp : ∀ B : Matrix (Fin p) (Fin q) ℝ, R4C.cvec y ⬝ᵥ (K B *ᵥ R4C.cvec y)
      = ∑ i, ∑ j, (y i : ℂ) * (K B i j * (y j : ℂ)) := by
    intro B
    simp only [dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum]
  simp only [hexp]
  rw [integral_finsetSum _ fun i _ => integrable_finsetSum _ fun j _ => hint i j]
  have hrow : ∀ i : Fin p,
      ∫ B, ∑ j, (y i : ℂ) * (K B i j * (y j : ℂ)) ∂(gaussianMatrix p q)
        = ((y i ^ 2 : ℝ) : ℂ) * meanK K i i := by
    intro i
    rw [integral_finsetSum _ fun j _ => hint i j]
    have hterm : ∀ j : Fin p,
        ∫ B, (y i : ℂ) * (K B i j * (y j : ℂ)) ∂(gaussianMatrix p q)
          = if j = i then ((y i ^ 2 : ℝ) : ℂ) * meanK K i i else 0 := by
      intro j
      rw [integral_const_mul, MeasureTheory.integral_mul_const]
      by_cases hji : j = i
      · rw [if_pos hji, hji]
        simp only [meanK]
        push_cast
        ring
      · rw [if_neg hji]
        have hz0 : meanK K i j = 0 := h.meanK_offDiag (Ne.symm hji)
        simp only [meanK] at hz0
        rw [hz0]
        ring
    simp only [hterm]
    rw [Finset.sum_ite_eq' Finset.univ i fun _ => ((y i ^ 2 : ℝ) : ℂ) * meanK K i i]
    simp
  exact Finset.sum_congr rfl fun i _ => hrow i

/-- The block average of the diagonal means. -/
noncomputable def meanKAvg (K : Matrix (Fin p) (Fin q) ℝ → Matrix (Fin p) (Fin p) ℂ)
    (J : Finset (Fin p)) : ℂ :=
  ∫ B, (J.card : ℂ)⁻¹ * ∑ j ∈ J, K B j j ∂gaussianMatrix p q

theorem BlockSym.meanKAvg_eq (h : BlockSym τ q K C) {J : Finset (Fin p)}
    (hJ : ∀ j ∈ J, ∀ k ∈ J, τ j = τ k) {j : Fin p} (hj : j ∈ J) :
    meanKAvg K J = meanK K j j := by
  have hcard : (J.card : ℂ) ≠ 0 := by
    exact_mod_cast (Finset.card_pos.mpr ⟨j, hj⟩).ne'
  rw [meanKAvg, integral_const_mul, integral_finsetSum _ fun k _ => h.integrable_entry k k]
  have hall : ∀ k ∈ J, ∫ B, K B k k ∂(gaussianMatrix p q) = meanK K j j :=
    fun k hk => h.meanK_diag_eq (hJ k hk j hj)
  rw [Finset.sum_congr rfl hall, Finset.sum_const, nsmul_eq_mul, ← mul_assoc,
    inv_mul_cancel₀ hcard, one_mul]

/-- **The mean in block form**: `E[yᵀ K y] = ∑_i (∑_{j ∈ J_i} y_j²) E[K_{J_i}]`. -/
theorem BlockSym.integral_qform_blocks {M : ℕ} {w : Fin M → ℝ} {blk : Fin p → Fin M}
    (h : BlockSym (tauOf w blk) q K C) (y : Fin p → ℝ) :
    ∫ B, R4C.cvec y ⬝ᵥ (K B *ᵥ R4C.cvec y) ∂gaussianMatrix p q
      = ∑ i, ((∑ j ∈ blockSet blk i, y j ^ 2 : ℝ) : ℂ) * meanKAvg K (blockSet blk i) := by
  rw [h.integral_qform, ← sum_blockSet blk]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Complex.ofReal_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun j hj => ?_
  congr 1
  refine (h.meanKAvg_eq (fun a ha b hb => ?_) hj).symm
  simp only [blockSet, Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
  simp [tauOf, ha, hb]

/-- `(G²)_{jk}` as a form, and its bound. -/
theorem resolvC_sq_apply_eq_cform2C {n : ℕ} (W : Matrix (Fin n) (Fin n) ℝ) (z : ℂ)
    (i j : Fin n) :
    (R4C.resolvC W z * R4C.resolvC W z) i j
      = R4C.cform2C W z (Pi.single i 1) (Pi.single j 1) := by
  rw [cform2C_eq_bil, bil_eq_sum]
  set ei : Fin n → ℝ := Pi.single i 1 with hei
  set ej : Fin n → ℝ := Pi.single j 1 with hej
  have h1 : ∀ a : Fin n, ((ei a : ℝ) : ℂ) = if a = i then 1 else 0 := by
    intro a
    by_cases h : a = i <;> simp [hei, h]
  have h2 : ∀ b : Fin n, ((ej b : ℝ) : ℂ) = if b = j then 1 else 0 := by
    intro b
    by_cases h : b = j <;> simp [hej, h]
  simp only [h1, h2, ite_mul, one_mul, zero_mul, mul_ite, mul_one, mul_zero,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]

theorem norm_resolvC_sq_entry_le {n : ℕ} {W : Matrix (Fin n) (Fin n) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (i j : Fin n) :
    ‖(R4C.resolvC W z * R4C.resolvC W z) i j‖ ≤ 1 / z.im ^ 2 := by
  have hs : ∀ k : Fin n, (Pi.single k (1 : ℝ)) ⬝ᵥ (Pi.single k (1 : ℝ)) = 1 := by
    intro k
    simp [dotProduct, Pi.single_apply]
  have h := R4C.norm_cform2C_le hW hz (Pi.single i (1 : ℝ)) (Pi.single j (1 : ℝ))
  rw [hs, hs] at h
  rw [resolvC_sq_apply_eq_cform2C]
  simpa using h

/-- `G_σ` is block-symmetric. -/
theorem blockSym_Gsig (hz : 0 < z.im) (hp : 0 < p) {d : ℕ} (hd : 0 < d) (τ : Fin p → ℝ) :
    BlockSym τ q (Gsig τ d z) (1 / z.im) where
  conj σ ε hε hτ B := by
    rw [Gsig, Wsig_rowSgn τ d σ ε hτ B, resolvC_conj (sgnPermMat_orth hε) (isHermitian_Wsig τ d B)
      hz.ne']
    rfl
  meas j k := measurable_Gsig_entry hz hp hd τ j k
  bound B j k := norm_resolvC_entry_le (isHermitian_Wsig τ d B) hz j k

/-- `G_σ²` is block-symmetric. -/
theorem blockSym_Gsig_sq (hz : 0 < z.im) (hp : 0 < p) {d : ℕ} (hd : 0 < d) (τ : Fin p → ℝ) :
    BlockSym τ q (fun B => Gsig τ d z B * Gsig τ d z B) (1 / z.im ^ 2) where
  conj σ ε hε hτ B := by
    simp only [Gsig]
    rw [Wsig_rowSgn τ d σ ε hτ B, resolvC_sq_conj (sgnPermMat_orth hε) (isHermitian_Wsig τ d B)
      hz.ne']
  meas j k := by
    have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => (Gsig τ d z B * Gsig τ d z B) j k)
        = fun B => ∑ l, Gsig τ d z B j l * Gsig τ d z B l k := by
      funext B
      rw [Matrix.mul_apply]
    rw [heq]
    exact Finset.measurable_sum _ fun l _ =>
      (measurable_Gsig_entry hz hp hd τ j l).mul (measurable_Gsig_entry hz hp hd τ l k)
  bound B j k := norm_resolvC_sq_entry_le (isHermitian_Wsig τ d B) hz j k

end Mean

/-! ### Bounded convergence: the means of a sequence that converges in probability -/

section Bounded

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN : ℕ → ℕ}

/-- A uniformly bounded model statistic that converges in probability has convergent means. -/
theorem tendsto_integral_of_tendstoInProb (B : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (qN N)) ℝ)
    (hB : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (qN N)) (μ N))
    {g : (N : ℕ) → Matrix (Fin (pN N)) (Fin (qN N)) ℝ → ℂ} {C : ℝ} {ℓ : ℂ}
    (hmeas : ∀ᶠ N in atTop, Measurable (g N)) (hbound : ∀ N B, ‖g N B‖ ≤ C)
    (hprob : TendstoInProb μ (fun N ω => ‖g N (B N ω) - ℓ‖) 0) :
    Tendsto (fun N => ∫ B, g N B ∂gaussianMatrix (pN N) (qN N)) atTop (𝓝 ℓ) := by
  refine Metric.tendsto_nhds.mpr fun ε hε => ?_
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hbound 0 0)
  set D : ℝ := C + ‖ℓ‖ + 1 with hDdef
  have hD : 0 < D := by positivity
  have hsmall := hprob (ε / 2) (by positivity)
  have hev : ∀ᶠ N in atTop,
      μ N {ω | ε / 2 ≤ |‖g N (B N ω) - ℓ‖ - 0|} < ENNReal.ofReal (ε / (2 * D)) :=
    hsmall.eventually (gt_mem_nhds (ENNReal.ofReal_pos.mpr (by positivity)))
  filter_upwards [hev, hmeas] with N hN hmN
  set ν := gaussianMatrix (pN N) (qN N) with hνdef
  have hint : Integrable (g N) ν :=
    Integrable.of_bound hmN.aestronglyMeasurable C (Filter.Eventually.of_forall (hbound N))
  have hS : MeasurableSet {B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ | ε / 2 ≤ ‖g N B - ℓ‖} :=
    measurableSet_le measurable_const (hmN.sub_const ℓ).norm
  have h1 : ∫ B, g N B ∂ν - ℓ = ∫ B, (g N B - ℓ) ∂ν := by
    rw [integral_sub hint (integrable_const ℓ), integral_const]
    simp
  have hpt : ∀ B, ‖g N B - ℓ‖
      ≤ ε / 2 + D * ({B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ | ε / 2 ≤ ‖g N B - ℓ‖}.indicator
          (fun _ => (1 : ℝ)) B) := by
    intro B
    by_cases hb : ε / 2 ≤ ‖g N B - ℓ‖
    · rw [Set.indicator_of_mem (show B ∈ {B | ε / 2 ≤ ‖g N B - ℓ‖} from hb)]
      have := hbound N B
      calc ‖g N B - ℓ‖ ≤ ‖g N B‖ + ‖ℓ‖ := norm_sub_le _ _
        _ ≤ ε / 2 + D * 1 := by rw [hDdef]; linarith
    · rw [Set.indicator_of_notMem (show B ∉ {B | ε / 2 ≤ ‖g N B - ℓ‖} from hb)]
      rw [not_le] at hb
      linarith
  have h2 : ∫ B, ‖g N B - ℓ‖ ∂ν
      ≤ ε / 2 + D * (ν {B | ε / 2 ≤ ‖g N B - ℓ‖}).toReal := by
    calc ∫ B, ‖g N B - ℓ‖ ∂ν
        ≤ ∫ B, (ε / 2 + D * ({B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ |
            ε / 2 ≤ ‖g N B - ℓ‖}.indicator (fun _ => (1 : ℝ)) B)) ∂ν :=
          integral_mono (hint.sub (integrable_const ℓ)).norm
            ((integrable_const _).add (((integrable_const 1).indicator hS).const_mul D)) hpt
      _ = ε / 2 + D * (ν {B | ε / 2 ≤ ‖g N B - ℓ‖}).toReal := by
          rw [integral_add (integrable_const _) (((integrable_const 1).indicator hS).const_mul D),
            integral_const, integral_const_mul, integral_indicator_const _ hS]
          simp [measureReal_def]
  have h3 : (ν {B | ε / 2 ≤ ‖g N B - ℓ‖}).toReal < ε / (2 * D) := by
    have heq := (hB N).measure_eq (p := fun B => ε / 2 ≤ ‖g N B - ℓ‖) hS
    have hset : {ω | ε / 2 ≤ |‖g N (B N ω) - ℓ‖ - 0|} = {ω | ε / 2 ≤ ‖g N (B N ω) - ℓ‖} := by
      ext ω; simp
    rw [hset, heq] at hN
    exact ENNReal.toReal_lt_of_lt_ofReal hN
  rw [dist_eq_norm, h1]
  calc ‖∫ B, (g N B - ℓ) ∂ν‖ ≤ ∫ B, ‖g N B - ℓ‖ ∂ν := norm_integral_le_integral_norm _
    _ ≤ ε / 2 + D * (ν {B | ε / 2 ≤ ‖g N B - ℓ‖}).toReal := h2
    _ < ε / 2 + D * (ε / (2 * D)) := by gcongr
    _ = ε := by field_simp; ring

end Bounded

/-! ### The fixed-vector forms at complex `z` (block means and concentration) -/

section Fixed

variable {M : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

/-- `E[gsig2Avg_{J_i}]`, the second-order block mean. -/
noncomputable def meanG2 {p : ℕ} (w : Fin M → ℝ) (blk : Fin p → Fin M) (q d : ℕ) (z : ℂ)
    (i : Fin M) : ℂ :=
  ∫ B, gsig2Avg (tauOf w blk) d z (blockSet blk i) B ∂gaussianMatrix p q

theorem meanKAvg_Gsig {p q : ℕ} (w : Fin M → ℝ) (blk : Fin p → Fin M) (d : ℕ) (i : Fin M) :
    meanKAvg (q := q) (Gsig (tauOf w blk) d z) (blockSet blk i) = meanG w blk q d z i := rfl

theorem meanKAvg_Gsig_sq {p q : ℕ} (w : Fin M → ℝ) (blk : Fin p → Fin M) (d : ℕ) (i : Fin M) :
    meanKAvg (q := q) (fun B => Gsig (tauOf w blk) d z B * Gsig (tauOf w blk) d z B)
      (blockSet blk i) = meanG2 w blk q d z i := rfl

theorem measurable_gsig2Avg {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (J : Finset (Fin p)) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => gsig2Avg τ d z J B) :=
  measurable_const.mul (Finset.measurable_sum _ fun j _ => measurable_Gsig_sq_diag hz hp hd τ j)

theorem measurable_qformC_Wsig {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qformC (Wsig τ d B) z y) := by
  have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qformC (Wsig τ d B) z y)
      = fun B => ∑ i, ∑ j, (y i : ℂ) * (Gsig τ d z B i j * (y j : ℂ)) := by
    funext B
    simp only [R4C.qformC, R4C.cformC, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum]
    rfl
  rw [heq]
  exact Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ =>
    ((measurable_Gsig_entry hz hp hd τ i j).mul_const _).const_mul _

theorem measurable_qform2C_Wsig {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qform2C (Wsig τ d B) z y) := by
  have heq : (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qform2C (Wsig τ d B) z y)
      = fun B => ∑ i, ∑ j,
        (y i : ℂ) * ((∑ l, Gsig τ d z B i l * Gsig τ d z B l j) * (y j : ℂ)) := by
    funext B
    simp only [R4C.qform2C, R4C.cform2C, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum,
      Matrix.mul_apply]
    rfl
  rw [heq]
  exact Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ =>
    ((Finset.measurable_sum _ fun l _ => (measurable_Gsig_entry hz hp hd τ i l).mul
      (measurable_Gsig_entry hz hp hd τ l j)).mul_const _).const_mul _

theorem integrable_qformC_Wsig {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qformC (Wsig τ d B) z y)
      (gaussianMatrix p q) :=
  Integrable.of_bound (measurable_qformC_Wsig hz hp hd τ y).aestronglyMeasurable
    ((y ⬝ᵥ y) / z.im)
    (Filter.Eventually.of_forall fun B => R4C.norm_qformC_le (isHermitian_Wsig τ d B) hz y)

theorem integrable_qform2C_Wsig {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => R4C.qform2C (Wsig τ d B) z y)
      (gaussianMatrix p q) :=
  Integrable.of_bound (measurable_qform2C_Wsig hz hp hd τ y).aestronglyMeasurable
    ((y ⬝ᵥ y) / z.im ^ 2)
    (Filter.Eventually.of_forall fun B => R4C.norm_qform2C_le (isHermitian_Wsig τ d B) hz y)

theorem sum_blockSet_real {p : ℕ} (blk : Fin p → Fin M) (f : Fin p → ℝ) :
    ∑ i, ∑ j ∈ blockSet blk i, f j = ∑ j, f j := by
  classical
  exact Finset.sum_fiberwise Finset.univ blk f

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))
  (B : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (qN N)) ℝ)
  (hB : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (qN N)) (μ N))

include hc hw hz hd hq hn hB in
/-- The second-order block means converge: `E[gsig2Avg_i] → g_i'` (bounded convergence from
`tendstoInProb_gsig2Avg`). -/
theorem tendsto_meanG2 (i : Fin M) :
    Tendsto (fun N => meanG2 w (blk N) (qN N) (dN N) z i) atTop
      (𝓝 (gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹)) :=
  tendsto_integral_of_tendstoInProb B hB
    (g := fun N B => gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B)
    (C := 1 / z.im ^ 2)
    (by
      filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
      exact measurable_gsig2Avg hz hN.2.2.2 hN.1 _ _)
    (fun N B => norm_gsig2Avg_le hz _ _ _ B)
    (tendstoInProb_gsig2Avg hc hw hz hd hq hn B hB i)

include hc hw hz hd hq hn hB in
/-- **The fixed-vector form at complex `z`.** For deterministic `y_N` with block norms
`∑_{j ∈ J_i} y_j² = a_i`, `y_Nᵀ G_σ(z) y_N → ∑ a_i g_i(z)` in probability. -/
theorem tendstoInProb_qformC_fixed (y : ∀ N, Fin (pN N) → ℝ) (a : Fin M → ℝ)
    (hy : ∀ N i, ∑ j ∈ blockSet (blk N) i, y N j ^ 2 = a i) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (Wsig (tauOf w (blk N)) (dN N) (B N ω)) z (y N)
      - ∑ i, ((a i : ℝ) : ℂ) * gC w i z (sGlob c w z)‖) 0 := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  have hyy : ∀ N, y N ⬝ᵥ y N = ∑ i, a i := by
    intro N
    rw [← Finset.sum_congr rfl fun i _ => hy N i, sum_blockSet_real]
    simp only [dotProduct, sq]
  have ha0 : 0 ≤ ∑ i, a i := by
    rw [← hyy 0]
    exact Finset.sum_nonneg fun _ _ => mul_self_nonneg _
  refine tendstoInProb_of_lipschitz_mean B hB
    (g := fun N B => R4C.qformC (Wsig (tauOf w (blk N)) (dN N) B) z (y N))
    (L := fun N => lipC z (Scalars.wSqMax w) (∑ i, a i + 1) (∑ i, a i + 1) (dN N))
    (tendsto_lipC hd z _ _ _) ?_ ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    have hS : ∀ j, tauOf w (blk N) j ^ 2 ≤ Scalars.wSqMax w := tauOf_sq_le w (blk N)
    refine ⟨Nat.mul_pos hpN hqN, lipC_pos hz hS0 (by linarith) (by linarith) hdN,
      measurable_qformC_Wsig hz hpN hdN _ (y N), integrable_qformC_Wsig hz hpN hdN _ (y N), ?_⟩
    obtain ⟨hre, him⟩ := lipschitzWith_of_norm_sub_le (lipC_nonneg z _ _ _ (dN N))
      (fun B B' => norm_qformC_Wsig_sub_le hz (tauOf w (blk N)) hS hS0.le hdN B B' (y N)
        (ny := ∑ i, a i + 1) (by rw [hyy]; linarith))
    exact ⟨_, _, hre, him, fun B => by simp⟩
  · have hmean : ∀ᶠ N in atTop,
        ∫ B, R4C.qformC (Wsig (tauOf w (blk N)) (dN N) B) z (y N)
            ∂gaussianMatrix (pN N) (qN N)
          = ∑ i, ((a i : ℝ) : ℂ) * meanG w (blk N) (qN N) (dN N) z i := by
      filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
      obtain ⟨hdN, hqN, -, hpN⟩ := hN
      have h := (blockSym_Gsig (q := qN N) hz hpN hdN (tauOf w (blk N))).integral_qform_blocks
        (y N)
      simp only [hy, meanKAvg_Gsig] at h
      exact h
    refine (tendsto_finsetSum _ fun i _ =>
      (tendsto_meanG hc hw hz hd hq hn i).const_mul (((a i : ℝ) : ℂ))).congr' ?_
    filter_upwards [hmean] with N hN
    exact hN.symm

include hc hw hz hd hq hn hB in
/-- **The second-order fixed-vector form at complex `z`.** -/
theorem tendstoInProb_qform2C_fixed (y : ∀ N, Fin (pN N) → ℝ) (a : Fin M → ℝ)
    (hy : ∀ N i, ∑ j ∈ blockSet (blk N) i, y N j ^ 2 = a i) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) (B N ω)) z (y N)
      - ∑ i, ((a i : ℝ) : ℂ) * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹‖)
      0 := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  have hyy : ∀ N, y N ⬝ᵥ y N = ∑ i, a i := by
    intro N
    rw [← Finset.sum_congr rfl fun i _ => hy N i, sum_blockSet_real]
    simp only [dotProduct, sq]
  have ha0 : 0 ≤ ∑ i, a i := by
    rw [← hyy 0]
    exact Finset.sum_nonneg fun _ _ => mul_self_nonneg _
  refine tendstoInProb_of_lipschitz_mean B hB
    (g := fun N B => R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) B) z (y N))
    (L := fun N => 2 * lipC z (Scalars.wSqMax w) ((∑ i, a i + 1) / z.im ^ 2) (∑ i, a i + 1)
      (dN N))
    (by simpa using (tendsto_lipC hd z (Scalars.wSqMax w) ((∑ i, a i + 1) / z.im ^ 2)
      (∑ i, a i + 1)).const_mul 2) ?_ ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    have hS : ∀ j, tauOf w (blk N) j ^ 2 ≤ Scalars.wSqMax w := tauOf_sq_le w (blk N)
    have hL2 : 0 < 2 * lipC z (Scalars.wSqMax w) ((∑ i, a i + 1) / z.im ^ 2) (∑ i, a i + 1)
        (dN N) := by
      have hpos : 0 < ∑ i, a i + 1 := by linarith
      have := lipC_pos hz hS0 (div_pos hpos (pow_pos hz 2)) hpos hdN
      linarith
    refine ⟨Nat.mul_pos hpN hqN, hL2,
      measurable_qform2C_Wsig hz hpN hdN _ (y N), integrable_qform2C_Wsig hz hpN hdN _ (y N), ?_⟩
    obtain ⟨hre, him⟩ := lipschitzWith_of_norm_sub_le hL2.le
      (fun B B' => norm_qform2C_Wsig_sub_le hz (tauOf w (blk N)) hS hS0.le hdN B B' (y N)
        (ny := ∑ i, a i + 1) (by rw [hyy]; linarith))
    exact ⟨_, _, hre, him, fun B => by simp⟩
  · have hmean : ∀ᶠ N in atTop,
        ∫ B, R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) B) z (y N)
            ∂gaussianMatrix (pN N) (qN N)
          = ∑ i, ((a i : ℝ) : ℂ) * meanG2 w (blk N) (qN N) (dN N) z i := by
      filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
      obtain ⟨hdN, hqN, -, hpN⟩ := hN
      have h := (blockSym_Gsig_sq (q := qN N) hz hpN hdN
        (tauOf w (blk N))).integral_qform_blocks (y N)
      simp only [hy, meanKAvg_Gsig_sq] at h
      exact h
    refine (tendsto_finsetSum _ fun i _ =>
      (tendsto_meanG2 hc hw hz hd hq hn B hB i).const_mul (((a i : ℝ) : ℂ))).congr' ?_
    filter_upwards [hmean] with N hN
    exact hN.symm

end Fixed

/-! ### Gaussian quadratic and linear forms: second moments (R2 (c), conditional on `B`) -/

section GaussForm

variable {p : ℕ}

theorem isHermitian_of_transpose_eq {A : Matrix (Fin p) (Fin p) ℝ} (h : Aᵀ = A) :
    A.IsHermitian := by
  rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial]
  exact h

/-- `tr S = ∑ λ_a` for a real symmetric `S`. -/
theorem trace_eq_sum_eigenvalues {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian) :
    S.trace = ∑ a, hS.eigenvalues a := by
  conv_lhs => rw [← eigU_conj hS]
  rw [Matrix.trace_mul_comm, ← Matrix.mul_assoc, transpose_eigU_mul hS, Matrix.one_mul,
    Matrix.trace_diagonal]

/-- The Gaussian second moment of a real quadratic form: `E (xᵀ S x - tr S)² = varSq ∑ λ_a²`. -/
theorem integral_sq_qform_sub_trace {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian) :
    ∫ x, (x ⬝ᵥ (S *ᵥ x) - S.trace) ^ 2 ∂(piGauss p)
      = varSq * ∑ a, hS.eigenvalues a ^ 2 := by
  have hexp : ∀ x : Fin p → ℝ, x ⬝ᵥ (S *ᵥ x) - S.trace
      = ∑ a, hS.eigenvalues a * ((fun t : ℝ => t ^ 2 - 1) (((eigU hS)ᵀ *ᵥ x) a)) := by
    intro x
    have hq := R4.dotProduct_conj (eigU hS) hS.eigenvalues x
    rw [eigU_conj hS] at hq
    rw [hq, trace_eq_sum_eigenvalues hS, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun a _ => by ring
  simp only [hexp]
  rw [integral_comp_mulVec (R2.transpose_eigU_orth hS)
    (F := fun y => (∑ a, hS.eigenvalues a * ((fun t : ℝ => t ^ 2 - 1) (y a))) ^ 2)
    (centered_sq_sub_one.integrable_sq_sum _).aestronglyMeasurable]
  exact centered_sq_sub_one.integral_sq_sum _

theorem integrable_sq_qform_sub_trace {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian) :
    Integrable (fun x => (x ⬝ᵥ (S *ᵥ x) - S.trace) ^ 2) (piGauss p) := by
  have hexp : ∀ x : Fin p → ℝ, x ⬝ᵥ (S *ᵥ x) - S.trace
      = ∑ a, hS.eigenvalues a * ((fun t : ℝ => t ^ 2 - 1) (((eigU hS)ᵀ *ᵥ x) a)) := by
    intro x
    have hq := R4.dotProduct_conj (eigU hS) hS.eigenvalues x
    rw [eigU_conj hS] at hq
    rw [hq, trace_eq_sum_eigenvalues hS, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun a _ => by ring
  simp only [hexp]
  exact integrable_comp_mulVec (R2.transpose_eigU_orth hS)
    (centered_sq_sub_one.integrable_sq_sum _)

/-- Every eigenvalue is bounded by the Rayleigh bound on unit vectors. -/
theorem abs_eigenvalue_le {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian) {C : ℝ}
    (hC : ∀ u : Fin p → ℝ, u ⬝ᵥ u = 1 → |u ⬝ᵥ (S *ᵥ u)| ≤ C) (a : Fin p) :
    |hS.eigenvalues a| ≤ C := by
  have hdiag : (eigU hS)ᵀ * S * eigU hS = Matrix.diagonal hS.eigenvalues := by
    have h := congrArg (fun A => (eigU hS)ᵀ * A * eigU hS) (eigU_conj hS)
    rw [← h, show (eigU hS)ᵀ * (eigU hS * Matrix.diagonal hS.eigenvalues * (eigU hS)ᵀ) * eigU hS
        = ((eigU hS)ᵀ * eigU hS) * Matrix.diagonal hS.eigenvalues * ((eigU hS)ᵀ * eigU hS) by
        simp only [Matrix.mul_assoc], transpose_eigU_mul hS, Matrix.one_mul, Matrix.mul_one]
  set col : Fin p → ℝ := fun j => eigU hS j a with hcol
  have hunit : col ⬝ᵥ col = 1 := by
    have h := congrFun (congrFun (transpose_eigU_mul hS) a) a
    rw [Matrix.mul_apply, Matrix.one_apply_eq] at h
    rw [← h]
    exact Finset.sum_congr rfl fun j _ => by simp [hcol, Matrix.transpose_apply]
  have hentry : ((eigU hS)ᵀ * S * eigU hS) a a = col ⬝ᵥ (S *ᵥ col) := by
    simp only [Matrix.mul_apply, dotProduct, Matrix.mulVec, Matrix.transpose_apply, hcol,
      Finset.sum_mul, Finset.mul_sum]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => by ring
  have hlam : hS.eigenvalues a = col ⬝ᵥ (S *ᵥ col) := by
    rw [← hentry, hdiag, Matrix.diagonal_apply_eq]
  rw [hlam]
  exact hC col hunit

theorem norm_ofReal_add_I_mul_sq (u v : ℝ) :
    ‖(u : ℂ) + Complex.I * (v : ℂ)‖ ^ 2 = u ^ 2 + v ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im]
  ring

/-- The real and imaginary parts of a complex kernel. -/
noncomputable def reK (K : Matrix (Fin p) (Fin p) ℂ) : Matrix (Fin p) (Fin p) ℝ :=
  Matrix.of fun j k => (K j k).re

noncomputable def imK (K : Matrix (Fin p) (Fin p) ℂ) : Matrix (Fin p) (Fin p) ℝ :=
  Matrix.of fun j k => (K j k).im

theorem K_eq_re_add_im (K : Matrix (Fin p) (Fin p) ℂ) :
    K = cmapR (reK K) + Complex.I • cmapR (imK K) := by
  ext j k
  simp only [Matrix.add_apply, Matrix.smul_apply, cmapR_apply, reK, imK, Matrix.of_apply,
    smul_eq_mul]
  rw [mul_comm, Complex.re_add_im]

/-- `yᵀ K y` splits into the two real forms. -/
theorem cvec_dotProduct_K_mulVec (K : Matrix (Fin p) (Fin p) ℂ) (y : Fin p → ℝ) :
    R4C.cvec y ⬝ᵥ (K *ᵥ R4C.cvec y)
      = ((y ⬝ᵥ (reK K *ᵥ y) : ℝ) : ℂ) + Complex.I * ((y ⬝ᵥ (imK K *ᵥ y) : ℝ) : ℂ) := by
  conv_lhs => rw [K_eq_re_add_im K]
  rw [Matrix.add_mulVec, Matrix.smul_mulVec, dotProduct_add, dotProduct_smul, smul_eq_mul,
    ← cmat_eq_cmapR, ← cmat_eq_cmapR, R4C.dotProduct_cmat_mulVec, R4C.dotProduct_cmat_mulVec]

theorem transpose_reK {K : Matrix (Fin p) (Fin p) ℂ} (hK : Kᵀ = K) : (reK K)ᵀ = reK K := by
  ext j k
  simp only [Matrix.transpose_apply, reK, Matrix.of_apply]
  conv_rhs => rw [← hK]
  rfl

theorem transpose_imK {K : Matrix (Fin p) (Fin p) ℂ} (hK : Kᵀ = K) : (imK K)ᵀ = imK K := by
  ext j k
  simp only [Matrix.transpose_apply, imK, Matrix.of_apply]
  conv_rhs => rw [← hK]
  rfl

/-- The scaled vector `δ ⊙ x`. -/
def scaled (δ x : Fin p → ℝ) : Fin p → ℝ := fun j => δ j * x j

/-- `(δ ⊙ x)ᵀ A (δ ⊙ x) = xᵀ (D A D) x`. -/
theorem scaled_qform_eq (A : Matrix (Fin p) (Fin p) ℝ) (δ x : Fin p → ℝ) :
    scaled δ x ⬝ᵥ (A *ᵥ scaled δ x)
      = x ⬝ᵥ ((Matrix.diagonal δ * A * Matrix.diagonal δ) *ᵥ x) := by
  have hD : Matrix.diagonal δ *ᵥ x = scaled δ x := by
    funext j
    rw [Matrix.mulVec_diagonal]
    rfl
  have hxD : x ᵥ* Matrix.diagonal δ = scaled δ x := by
    funext j
    rw [Matrix.vecMul_diagonal, mul_comm]
    rfl
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, hD, Matrix.dotProduct_mulVec
    (v := x) (A := Matrix.diagonal δ), hxD]

theorem scaled_dotProduct_self_le (δ : Fin p → ℝ) {T : ℝ} (hT : ∀ j, δ j ^ 2 ≤ T)
    (u : Fin p → ℝ) : scaled δ u ⬝ᵥ scaled δ u ≤ T * (u ⬝ᵥ u) := by
  simp only [scaled, dotProduct, Finset.mul_sum]
  refine Finset.sum_le_sum fun j _ => ?_
  have : δ j * u j * (δ j * u j) = δ j ^ 2 * (u j * u j) := by ring
  rw [this]
  exact mul_le_mul_of_nonneg_right (hT j) (mul_self_nonneg _)

/-- **The second moment of a Gaussian quadratic form** `(δ ⊙ x)ᵀ K (δ ⊙ x)` about its mean
`∑ δ_j² K_jj`: at most `2 varSq p (T Cq)²` when `δ_j² ≤ T` and `|uᵀ K u| ≤ Cq ‖u‖²`. -/
theorem integral_normSq_qform_scaled (K : Matrix (Fin p) (Fin p) ℂ) (hKt : Kᵀ = K) {Cq : ℝ}
    (hK : ∀ u : Fin p → ℝ, ‖R4C.cvec u ⬝ᵥ (K *ᵥ R4C.cvec u)‖ ≤ Cq * (u ⬝ᵥ u))
    (δ : Fin p → ℝ) {T : ℝ} (hT : ∀ j, δ j ^ 2 ≤ T) :
    Integrable (fun x => ‖R4C.cvec (scaled δ x) ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))
        - ∑ j, ((δ j ^ 2 : ℝ) : ℂ) * K j j‖ ^ 2) (piGauss p) ∧
      ∫ x, ‖R4C.cvec (scaled δ x) ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))
        - ∑ j, ((δ j ^ 2 : ℝ) : ℂ) * K j j‖ ^ 2 ∂(piGauss p)
      ≤ 2 * varSq * (p * (T * Cq) ^ 2) := by
  set SR := Matrix.diagonal δ * reK K * Matrix.diagonal δ with hSR
  set SI := Matrix.diagonal δ * imK K * Matrix.diagonal δ with hSI
  have hSRh : SR.IsHermitian := isHermitian_of_transpose_eq (by
    rw [hSR, Matrix.transpose_mul, Matrix.transpose_mul, Matrix.diagonal_transpose,
      transpose_reK hKt, Matrix.mul_assoc])
  have hSIh : SI.IsHermitian := isHermitian_of_transpose_eq (by
    rw [hSI, Matrix.transpose_mul, Matrix.transpose_mul, Matrix.diagonal_transpose,
      transpose_imK hKt, Matrix.mul_assoc])
  -- the form and its mean, split into real and imaginary parts
  have hQ : ∀ x, R4C.cvec (scaled δ x) ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))
      = ((x ⬝ᵥ (SR *ᵥ x) : ℝ) : ℂ) + Complex.I * ((x ⬝ᵥ (SI *ᵥ x) : ℝ) : ℂ) := by
    intro x
    rw [cvec_dotProduct_K_mulVec, hSR, hSI, scaled_qform_eq, scaled_qform_eq]
  have hm : ∑ j, ((δ j ^ 2 : ℝ) : ℂ) * K j j
      = ((SR.trace : ℝ) : ℂ) + Complex.I * ((SI.trace : ℝ) : ℂ) := by
    have hdR : SR.trace = ∑ j, δ j ^ 2 * (K j j).re := by
      simp only [hSR, Matrix.trace, Matrix.diag, Matrix.mul_diagonal, Matrix.diagonal_mul, reK,
        Matrix.of_apply]
      exact Finset.sum_congr rfl fun j _ => by ring
    have hdI : SI.trace = ∑ j, δ j ^ 2 * (K j j).im := by
      simp only [hSI, Matrix.trace, Matrix.diag, Matrix.mul_diagonal, Matrix.diagonal_mul, imK,
        Matrix.of_apply]
      exact Finset.sum_congr rfl fun j _ => by ring
    rw [hdR, hdI]
    apply Complex.ext
    · simp only [Complex.re_sum, Complex.add_re, Complex.ofReal_re,
        Complex.mul_re, Complex.I_re, Complex.I_im, Complex.ofReal_im, zero_mul, one_mul,
        sub_zero, add_zero]
    · simp only [Complex.im_sum, Complex.add_im, Complex.ofReal_im,
        Complex.mul_im, Complex.I_re, Complex.I_im, Complex.ofReal_re, zero_mul, one_mul,
        zero_add, add_zero]
  have hsq : ∀ x, ‖R4C.cvec (scaled δ x) ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))
      - ∑ j, ((δ j ^ 2 : ℝ) : ℂ) * K j j‖ ^ 2
      = (x ⬝ᵥ (SR *ᵥ x) - SR.trace) ^ 2 + (x ⬝ᵥ (SI *ᵥ x) - SI.trace) ^ 2 := by
    intro x
    rw [hQ, hm, ← norm_ofReal_add_I_mul_sq]
    congr 1
    push_cast
    ring_nf
  -- the eigenvalue bounds
  have hRay : ∀ (A : Matrix (Fin p) (Fin p) ℝ) (S : Matrix (Fin p) (Fin p) ℝ),
      S = Matrix.diagonal δ * A * Matrix.diagonal δ →
      (∀ y : Fin p → ℝ, |y ⬝ᵥ (A *ᵥ y)| ≤ ‖R4C.cvec y ⬝ᵥ (K *ᵥ R4C.cvec y)‖) →
      ∀ u : Fin p → ℝ, u ⬝ᵥ u = 1 → |u ⬝ᵥ (S *ᵥ u)| ≤ T * Cq := by
    intro A S hS hA u hu
    have hCq : 0 ≤ Cq := by
      have h := hK u
      rw [hu, mul_one] at h
      exact (norm_nonneg _).trans h
    rw [hS, ← scaled_qform_eq]
    calc |scaled δ u ⬝ᵥ (A *ᵥ scaled δ u)|
        ≤ ‖R4C.cvec (scaled δ u) ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ u))‖ := hA _
      _ ≤ Cq * (scaled δ u ⬝ᵥ scaled δ u) := hK _
      _ ≤ Cq * (T * (u ⬝ᵥ u)) :=
          mul_le_mul_of_nonneg_left (scaled_dotProduct_self_le δ hT u) hCq
      _ = T * Cq := by rw [hu]; ring
  have hre : ∀ y : Fin p → ℝ, |y ⬝ᵥ (reK K *ᵥ y)| ≤ ‖R4C.cvec y ⬝ᵥ (K *ᵥ R4C.cvec y)‖ := by
    intro y
    rw [cvec_dotProduct_K_mulVec]
    have h := Complex.abs_re_le_norm (((y ⬝ᵥ (reK K *ᵥ y) : ℝ) : ℂ)
      + Complex.I * ((y ⬝ᵥ (imK K *ᵥ y) : ℝ) : ℂ))
    simpa using h
  have him : ∀ y : Fin p → ℝ, |y ⬝ᵥ (imK K *ᵥ y)| ≤ ‖R4C.cvec y ⬝ᵥ (K *ᵥ R4C.cvec y)‖ := by
    intro y
    rw [cvec_dotProduct_K_mulVec]
    have h := Complex.abs_im_le_norm (((y ⬝ᵥ (reK K *ᵥ y) : ℝ) : ℂ)
      + Complex.I * ((y ⬝ᵥ (imK K *ᵥ y) : ℝ) : ℂ))
    simpa using h
  have hsumR : ∑ a, hSRh.eigenvalues a ^ 2 ≤ p * (T * Cq) ^ 2 := by
    have : ∀ a, hSRh.eigenvalues a ^ 2 ≤ (T * Cq) ^ 2 := fun a => by
      rw [← sq_abs]
      exact pow_le_pow_left₀ (abs_nonneg _) (abs_eigenvalue_le hSRh (hRay _ _ rfl hre) a) 2
    calc ∑ a, hSRh.eigenvalues a ^ 2
        ≤ ∑ _a : Fin p, (T * Cq) ^ 2 := Finset.sum_le_sum fun a _ => this a
      _ = p * (T * Cq) ^ 2 := by simp
  have hsumI : ∑ a, hSIh.eigenvalues a ^ 2 ≤ p * (T * Cq) ^ 2 := by
    have : ∀ a, hSIh.eigenvalues a ^ 2 ≤ (T * Cq) ^ 2 := fun a => by
      rw [← sq_abs]
      exact pow_le_pow_left₀ (abs_nonneg _) (abs_eigenvalue_le hSIh (hRay _ _ rfl him) a) 2
    calc ∑ a, hSIh.eigenvalues a ^ 2
        ≤ ∑ _a : Fin p, (T * Cq) ^ 2 := Finset.sum_le_sum fun a _ => this a
      _ = p * (T * Cq) ^ 2 := by simp
  simp only [hsq]
  refine ⟨(integrable_sq_qform_sub_trace hSRh).add (integrable_sq_qform_sub_trace hSIh), ?_⟩
  rw [integral_add (integrable_sq_qform_sub_trace hSRh) (integrable_sq_qform_sub_trace hSIh),
    integral_sq_qform_sub_trace hSRh, integral_sq_qform_sub_trace hSIh]
  have hv := varSq_nonneg
  calc varSq * ∑ a, hSRh.eigenvalues a ^ 2 + varSq * ∑ a, hSIh.eigenvalues a ^ 2
      ≤ varSq * (p * (T * Cq) ^ 2) + varSq * (p * (T * Cq) ^ 2) :=
        add_le_add (mul_le_mul_of_nonneg_left hsumR hv) (mul_le_mul_of_nonneg_left hsumI hv)
    _ = 2 * varSq * (p * (T * Cq) ^ 2) := by ring

/-- **The second moment of a Gaussian linear form** `uᵀ K (δ ⊙ x)`: `≤ T ‖K u‖²`. -/
theorem integral_normSq_lin_scaled (K : Matrix (Fin p) (Fin p) ℂ) (hKt : Kᵀ = K)
    (u : Fin p → ℝ) {Cl : ℝ} (hK : nsq (K *ᵥ R4C.cvec u) ≤ Cl) (δ : Fin p → ℝ) {T : ℝ}
    (hT : ∀ j, δ j ^ 2 ≤ T) (hT0 : 0 ≤ T) :
    Integrable (fun x => ‖R4C.cvec u ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))‖ ^ 2) (piGauss p) ∧
      ∫ x, ‖R4C.cvec u ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))‖ ^ 2 ∂(piGauss p) ≤ T * Cl := by
  set cf : Fin p → ℂ := fun k => (K *ᵥ R4C.cvec u) k * (δ k : ℂ) with hcf
  have hlin : ∀ x, R4C.cvec u ⬝ᵥ (K *ᵥ R4C.cvec (scaled δ x))
      = ∑ k, cf k * (((fun t : ℝ => t) (x k) : ℝ) : ℂ) := by
    intro x
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hKt]
    simp only [dotProduct, R4C.cvec, scaled, hcf]
    push_cast
    exact Finset.sum_congr rfl fun k _ => by ring
  simp only [hlin]
  refine ⟨centered_id.integrable_normSq_sum cf, ?_⟩
  rw [centered_id.integral_normSq_sum cf]
  have h1 : ∫ t : ℝ, ((fun t : ℝ => t) t) ^ 2 ∂(gaussianReal 0 1) = 1 := integral_sq_gauss
  rw [h1, one_mul]
  calc ∑ k, ‖cf k‖ ^ 2 = ∑ k, δ k ^ 2 * ‖(K *ᵥ R4C.cvec u) k‖ ^ 2 := by
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [hcf]
        simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs]
        ring
    _ ≤ ∑ k, T * ‖(K *ᵥ R4C.cvec u) k‖ ^ 2 :=
        Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_right (hT k) (sq_nonneg _)
    _ = T * nsq (K *ᵥ R4C.cvec u) := by rw [nsq, Finset.mul_sum]
    _ ≤ T * Cl := mul_le_mul_of_nonneg_left hK hT0

end GaussForm

/-! ### The `g`-forms at complex `z`: conditional Chebyshev given `B` (R2 (c)) -/

section GaussSeq

variable {M : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

/-- The pair `(x, B)`: the Gaussian vector `√d e = Z v` and the block. -/
abbrev PairSpace (p q : ℕ) := (Fin p → ℝ) × Matrix (Fin p) (Fin q) ℝ

/-- The product law of the pair (`exists_block_hasLaw_het`). -/
noncomputable def pairLaw (p q : ℕ) : Measure (PairSpace p q) :=
  (piGauss p).prod (gaussianMatrix p q)

/-- `g = d^{-1/2} diag τ x`, in the shape of `SigmaHalf *ᵥ eHet`. -/
noncomputable def gvec {p : ℕ} (τ : Fin p → ℝ) (d : ℕ) (x : Fin p → ℝ) : Fin p → ℝ :=
  fun j => τ j * ((Real.sqrt d)⁻¹ * x j)

theorem gvec_eq_scaled {p : ℕ} (τ : Fin p → ℝ) (d : ℕ) (x : Fin p → ℝ) :
    gvec τ d x = scaled (fun j => τ j * (Real.sqrt d)⁻¹) x := by
  funext j
  simp only [gvec, scaled, mul_assoc]

theorem gvec_coef_sq_le {p : ℕ} (τ : Fin p → ℝ) {S : ℝ} (hS : ∀ j, τ j ^ 2 ≤ S) (d : ℕ)
    (j : Fin p) : (τ j * (Real.sqrt d)⁻¹) ^ 2 ≤ S / d := by
  rw [mul_pow, inv_pow, Real.sq_sqrt (Nat.cast_nonneg d), div_eq_mul_inv]
  exact mul_le_mul_of_nonneg_right (hS j) (by positivity)

/-- The centering of the `g`-form in block form:
`∑_j (τ_j²/d) K_jj = ∑_i w_i² (|J_i|/d) K_{J_i}`. -/
theorem sum_scaled_diag_eq {p : ℕ} (w : Fin M → ℝ) (blk : Fin p → Fin M) (d : ℕ)
    (K : Matrix (Fin p) (Fin p) ℂ) :
    ∑ j, (((tauOf w blk j * (Real.sqrt d)⁻¹) ^ 2 : ℝ) : ℂ) * K j j
      = ∑ i, ((w i ^ 2 * (((blockSet blk i).card : ℝ) / d) : ℝ) : ℂ)
          * (((blockSet blk i).card : ℂ)⁻¹ * ∑ j ∈ blockSet blk i, K j j) := by
  rw [← sum_blockSet blk]
  refine Finset.sum_congr rfl fun i _ => ?_
  rcases Nat.eq_zero_or_pos (blockSet blk i).card with h0 | hpos
  · rw [Finset.card_eq_zero] at h0
    simp [h0]
  · have hne : ((blockSet blk i).card : ℂ) ≠ 0 := by exact_mod_cast hpos.ne'
    have hterm : ∀ j ∈ blockSet blk i,
        (((tauOf w blk j * (Real.sqrt d)⁻¹) ^ 2 : ℝ) : ℂ) * K j j
          = ((w i ^ 2 / d : ℝ) : ℂ) * K j j := by
      intro j hj
      rw [mul_pow, inv_pow, Real.sq_sqrt (Nat.cast_nonneg d), tauOf_sq_eq w blk i hj,
        div_eq_mul_inv]
    rw [Finset.sum_congr rfl hterm, ← Finset.mul_sum, ← mul_assoc]
    congr 1
    push_cast
    field_simp

/-- **Chebyshev on the product law.** A second-moment bound in `x`, uniform in the block `B`
and tending to `0`, gives a limit in probability on the model spaces. -/
theorem tendstoInProb_of_integral_sq_le_pair
    (XB : ∀ N, Ω N → PairSpace (pN N) (qN N))
    (hXB : ∀ N, HasLaw (XB N) (pairLaw (pN N) (qN N)) (μ N))
    (F : ∀ N, (Fin (pN N) → ℝ) → Matrix (Fin (pN N)) (Fin (qN N)) ℝ → ℝ)
    (hFm : ∀ᶠ N in atTop, Measurable fun q : PairSpace (pN N) (qN N) => F N q.1 q.2)
    (hFnn : ∀ N x Y, 0 ≤ F N x Y)
    (hint : ∀ᶠ N in atTop, ∀ Y, Integrable (fun x => (F N x Y) ^ 2) (piGauss (pN N)))
    {K : ℕ → ℝ} (hK : Tendsto K atTop (𝓝 0))
    (hbnd : ∀ᶠ N in atTop, ∀ Y, ∫ x, (F N x Y) ^ 2 ∂(piGauss (pN N)) ≤ K N) :
    TendstoInProb μ (fun N ω => F N (XB N ω).1 (XB N ω).2) 0 := by
  intro ε hε
  have hK' : Tendsto (fun N => ENNReal.ofReal (K N / ε ^ 2)) atTop (𝓝 0) := by
    have := ENNReal.tendsto_ofReal (hK.div_const (ε ^ 2))
    simpa using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hK'
    (Filter.Eventually.of_forall fun _ => by simp) ?_
  filter_upwards [hbnd, hFm, hint] with N hN hFmN hintN
  have hset : {ω | ε ≤ |F N (XB N ω).1 (XB N ω).2 - 0|}
      = {ω | ε ≤ F N (XB N ω).1 (XB N ω).2} := by
    ext ω
    simp [abs_of_nonneg (hFnn N _ _)]
  have hmeasSet : MeasurableSet {q : PairSpace (pN N) (qN N) | ε ≤ F N q.1 q.2} :=
    measurableSet_le measurable_const hFmN
  rw [hset, (hXB N).measure_eq hmeasSet, pairLaw, Measure.prod_apply_symm hmeasSet]
  have hsec : ∀ Y : Matrix (Fin (pN N)) (Fin (qN N)) ℝ,
      (piGauss (pN N)) ((fun x => (x, Y)) ⁻¹' {q | ε ≤ F N q.1 q.2})
        ≤ ENNReal.ofReal (K N / ε ^ 2) := fun Y =>
    meas_ge_le_of_integral_sq _ (hFmN.comp measurable_prodMk_right) (hintN Y) hε (hN Y)
  calc ∫⁻ Y, (piGauss (pN N)) ((fun x => (x, Y)) ⁻¹' {q | ε ≤ F N q.1 q.2})
        ∂(gaussianMatrix (pN N) (qN N))
      ≤ ∫⁻ _Y : Matrix (Fin (pN N)) (Fin (qN N)) ℝ, ENNReal.ofReal (K N / ε ^ 2)
        ∂(gaussianMatrix (pN N) (qN N)) := lintegral_mono hsec
    _ = ENNReal.ofReal (K N / ε ^ 2) := by simp

/-- A deterministic sequence tending to `0` tends to `0` in probability. -/
theorem tendstoInProb_of_tendsto {r : ℕ → ℝ} (hr : Tendsto r atTop (𝓝 0)) :
    TendstoInProb μ (fun N _ => r N) 0 := by
  intro ε hε
  have hev : ∀ᶠ N in atTop, |r N - 0| < ε := by
    have := hr.eventually (Metric.ball_mem_nhds (0 : ℝ) hε)
    filter_upwards [this] with N hN
    simpa [Real.dist_eq] using hN
  refine tendsto_const_nhds.congr' ?_
  filter_upwards [hev] with N hN
  have : {ω : Ω N | ε ≤ |r N - 0|} = ∅ := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_le]
    exact hN
  rw [this, measure_empty]

/-- A finite sum of sequences tending to `0` in probability tends to `0`. -/
theorem tendstoInProb_finset_sum_zero {ι : Type*} [Fintype ι] {F : ι → ∀ N, Ω N → ℝ}
    (h : ∀ i, TendstoInProb μ (F i) 0) :
    TendstoInProb μ (fun N ω => ∑ i, F i N ω) 0 := by
  have hφ : ContinuousAt (fun u : ι → ℝ => ∑ i, u i) (fun _ => 0) :=
    (continuous_finsetSum Finset.univ fun i _ => continuous_apply i).continuousAt
  have hpi : TendstoInProbPi μ (fun N ω i => F i N ω) (fun _ => 0) := fun i => h i
  have := hpi.comp_continuous hφ
  simpa using this

theorem norm_mul_sub_mul_le (a b : ℝ) (X Y : ℂ) :
    ‖(a : ℂ) * X - (b : ℂ) * Y‖ ≤ |a - b| * ‖X‖ + |b| * ‖X - Y‖ := by
  have h : (a : ℂ) * X - (b : ℂ) * Y = ((a - b : ℝ) : ℂ) * X + (b : ℂ) * (X - Y) := by
    push_cast; ring
  rw [h]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs]

/-- The block marginal of the pair law. -/
theorem hasLaw_snd_of_pair {N : ℕ} {XB : Ω N → PairSpace (pN N) (qN N)}
    (hXB : HasLaw XB (pairLaw (pN N) (qN N)) (μ N)) :
    HasLaw (fun ω => (XB ω).2) (gaussianMatrix (pN N) (qN N)) (μ N) :=
  (MeasureTheory.measurePreserving_snd (μ := piGauss (pN N))
    (ν := gaussianMatrix (pN N) (qN N))).fun_comp_hasLaw hXB

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))
  (XB : ∀ N, Ω N → PairSpace (pN N) (qN N))
  (hXB : ∀ N, HasLaw (XB N) (pairLaw (pN N) (qN N)) (μ N))

include hd hn in
/-- `p_N (S / (d_N η))² → 0`. -/
theorem tendsto_varBound (S a : ℝ) :
    Tendsto (fun N => 2 * varSq * ((pN N : ℝ) * (S / dN N * a) ^ 2)) atTop (𝓝 0) := by
  have hpd : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 (∑ i, c i)) := by
    have := tendsto_finsetSum Finset.univ fun i _ => hn i
    refine this.congr fun N => ?_
    exact sum_cN blk dN N
  have hdR : Tendsto (fun N => ((dN N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have hinv : Tendsto (fun N => (S * a) ^ 2 / (dN N : ℝ)) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop hdR
  have hprod := (hpd.mul hinv).const_mul (2 * varSq)
  rw [mul_zero, mul_zero] at hprod
  refine hprod.congr' ?_
  filter_upwards [hd.eventually_gt_atTop 0] with N hN
  have hne : (dN N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  field_simp

theorem measurable_gvec_apply {p : ℕ} (τ : Fin p → ℝ) (d : ℕ) (j : Fin p) :
    Measurable fun x : Fin p → ℝ => gvec τ d x j :=
  ((measurable_pi_apply j).const_mul _).const_mul _

/-- Measurability of the `g`-forms on the pair space. -/
theorem measurable_qformC_gvec_pair {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) :
    Measurable fun x : PairSpace p q =>
      R4C.qformC (Wsig τ d x.2) z (gvec τ d x.1) := by
  have heq : (fun x : PairSpace p q => R4C.qformC (Wsig τ d x.2) z (gvec τ d x.1))
      = fun x => ∑ i, ∑ j, ((gvec τ d x.1 i : ℝ) : ℂ)
          * (Gsig τ d z x.2 i j * ((gvec τ d x.1 j : ℝ) : ℂ)) := by
    funext x
    simp only [R4C.qformC, R4C.cformC, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum]
    rfl
  rw [heq]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d i).comp measurable_fst)).mul
    (((measurable_Gsig_entry hz hp hd τ i j).comp measurable_snd).mul
      (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d j).comp measurable_fst)))

theorem measurable_qform2C_gvec_pair {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) :
    Measurable fun x : PairSpace p q =>
      R4C.qform2C (Wsig τ d x.2) z (gvec τ d x.1) := by
  have heq : (fun x : PairSpace p q => R4C.qform2C (Wsig τ d x.2) z (gvec τ d x.1))
      = fun x => ∑ i, ∑ j, ((gvec τ d x.1 i : ℝ) : ℂ)
          * ((∑ l, Gsig τ d z x.2 i l * Gsig τ d z x.2 l j) * ((gvec τ d x.1 j : ℝ) : ℂ)) := by
    funext x
    simp only [R4C.qform2C, R4C.cform2C, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum,
      Matrix.mul_apply]
    rfl
  rw [heq]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d i).comp measurable_fst)).mul
    ((Finset.measurable_sum _ fun l _ =>
      ((measurable_Gsig_entry hz hp hd τ i l).comp measurable_snd).mul
        ((measurable_Gsig_entry hz hp hd τ l j).comp measurable_snd)).mul
      (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d j).comp measurable_fst)))

theorem measurable_cformC_gvec_pair {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Measurable fun x : PairSpace p q =>
      R4C.cformC (Wsig τ d x.2) z y (gvec τ d x.1) := by
  have heq : (fun x : PairSpace p q => R4C.cformC (Wsig τ d x.2) z y (gvec τ d x.1))
      = fun x => ∑ i, ∑ j, (y i : ℂ) * (Gsig τ d z x.2 i j * ((gvec τ d x.1 j : ℝ) : ℂ)) := by
    funext x
    simp only [R4C.cformC, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum]
    rfl
  rw [heq]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact measurable_const.mul
    (((measurable_Gsig_entry hz hp hd τ i j).comp measurable_snd).mul
      (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d j).comp measurable_fst)))

theorem measurable_cform2C_gvec_pair {p q d : ℕ} (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (τ : Fin p → ℝ) (y : Fin p → ℝ) :
    Measurable fun x : PairSpace p q =>
      R4C.cform2C (Wsig τ d x.2) z y (gvec τ d x.1) := by
  have heq : (fun x : PairSpace p q => R4C.cform2C (Wsig τ d x.2) z y (gvec τ d x.1))
      = fun x => ∑ i, ∑ j, (y i : ℂ)
          * ((∑ l, Gsig τ d z x.2 i l * Gsig τ d z x.2 l j) * ((gvec τ d x.1 j : ℝ) : ℂ)) := by
    funext x
    simp only [R4C.cform2C, dotProduct, Matrix.mulVec, R4C.cvec, Finset.mul_sum,
      Matrix.mul_apply]
    rfl
  rw [heq]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact measurable_const.mul
    ((Finset.measurable_sum _ fun l _ =>
      ((measurable_Gsig_entry hz hp hd τ i l).comp measurable_snd).mul
        ((measurable_Gsig_entry hz hp hd τ l j).comp measurable_snd)).mul
      (Complex.measurable_ofReal.comp ((measurable_gvec_apply τ d j).comp measurable_fst)))


include hc hw hz hd hq hn hXB in
/-- **The `g`-form, centered at the random block means** (conditional Chebyshev given `B`). -/
theorem tendstoInProb_qformC_gvec_sub :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)
        - ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
            * gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (XB N ω).2‖) 0 := by
  have hcent : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
          * gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B
        = ∑ j, (((tauOf w (blk N) j * (Real.sqrt (dN N))⁻¹) ^ 2 : ℝ) : ℂ)
            * Gsig (tauOf w (blk N)) (dN N) z B j j :=
    fun N B => (sum_scaled_diag_eq w (blk N) (dN N) (Gsig (tauOf w (blk N)) (dN N) z B)).symm
  have hKq : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ) (u : Fin (pN N) → ℝ),
      ‖R4C.cvec u ⬝ᵥ (Gsig (tauOf w (blk N)) (dN N) z B *ᵥ R4C.cvec u)‖
        ≤ (1 / z.im) * (u ⬝ᵥ u) := by
    intro N B u
    have h := R4C.norm_qformC_le (isHermitian_Wsig (tauOf w (blk N)) (dN N) B) hz u
    rw [div_eq_inv_mul, ← one_div] at h
    exact h
  refine tendstoInProb_of_integral_sq_le_pair XB hXB
    (fun N x B => ‖R4C.qformC (Wsig (tauOf w (blk N)) (dN N) B) z (gvec (tauOf w (blk N)) (dN N) x)
      - ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
          * gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B‖)
    ?_ (fun _ _ _ => norm_nonneg _) ?_
    (K := fun N => 2 * varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * (1 / z.im)) ^ 2))
    (tendsto_varBound hd hn _ _) ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    exact ((measurable_qformC_gvec_pair hz hpN hdN _).sub
      (Finset.measurable_sum _ fun i _ => measurable_const.mul
        ((measurable_gsigAvg hz hpN hdN _ _).comp measurable_snd))).norm
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN B
    simp only [hcent N B, gvec_eq_scaled]
    exact (integral_normSq_qform_scaled (Gsig (tauOf w (blk N)) (dN N) z B)
      (transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne') (hKq N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N))).1
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN B
    simp only [hcent N B, gvec_eq_scaled]
    exact (integral_normSq_qform_scaled (Gsig (tauOf w (blk N)) (dN N) z B)
      (transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne') (hKq N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N))).2

include hc hw hz hd hq hn hXB in
/-- **The second-order `g`-form, centered at the random block means of `G_σ²`.** -/
theorem tendstoInProb_qform2C_gvec_sub :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)
        - ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
            * gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (XB N ω).2‖) 0 := by
  have hcent : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
          * gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B
        = ∑ j, (((tauOf w (blk N) j * (Real.sqrt (dN N))⁻¹) ^ 2 : ℝ) : ℂ)
            * (Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B) j j :=
    fun N B => (sum_scaled_diag_eq w (blk N) (dN N)
      (Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B)).symm
  have hKt : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      (Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B)ᵀ
        = Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B := by
    intro N B
    rw [Matrix.transpose_mul, Gsig, transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne']
  have hKq : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ) (u : Fin (pN N) → ℝ),
      ‖R4C.cvec u ⬝ᵥ ((Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B)
        *ᵥ R4C.cvec u)‖ ≤ (1 / z.im ^ 2) * (u ⬝ᵥ u) := by
    intro N B u
    have h := R4C.norm_qform2C_le (isHermitian_Wsig (tauOf w (blk N)) (dN N) B) hz u
    rw [div_eq_inv_mul, ← one_div] at h
    exact h
  refine tendstoInProb_of_integral_sq_le_pair XB hXB
    (fun N x B => ‖R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) B) z (gvec (tauOf w (blk N)) (dN N) x)
      - ∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ)
          * gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B‖)
    ?_ (fun _ _ _ => norm_nonneg _) ?_
    (K := fun N => 2 * varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * (1 / z.im ^ 2)) ^ 2))
    (tendsto_varBound hd hn _ _) ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    exact ((measurable_qform2C_gvec_pair hz hpN hdN _).sub
      (Finset.measurable_sum _ fun i _ => measurable_const.mul
        ((measurable_gsig2Avg hz hpN hdN _ _).comp measurable_snd))).norm
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN B
    simp only [hcent N B, gvec_eq_scaled]
    exact (integral_normSq_qform_scaled _ (hKt N B) (hKq N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N))).1
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN B
    simp only [hcent N B, gvec_eq_scaled]
    exact (integral_normSq_qform_scaled _ (hKt N B) (hKq N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N))).2

include hc hn in
omit hw hz hd hq in
/-- From the random block means to the scalar limits. -/
theorem tendstoInProb_blockMeans_sub {G : ∀ N, Ω N → Fin M → ℂ} {ℓ : Fin M → ℂ} {Cg : ℝ}
    (hG : ∀ N ω i, ‖G N ω i‖ ≤ Cg)
    (hlim : ∀ i, TendstoInProb μ (fun N ω => ‖G N ω i - ℓ i‖) 0) :
    TendstoInProb μ (fun N ω => ‖∑ i, ((w i ^ 2 * cN blk dN N i : ℝ) : ℂ) * G N ω i
      - ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ) * ℓ i‖) 0 := by
  refine TendstoInProb.of_le (g := fun N ω => ∑ i, w i ^ 2 *
      (|cN blk dN N i - c i| * Cg + c i * ‖G N ω i - ℓ i‖))
    (fun N => Filter.Eventually.of_forall fun ω => ?_) ?_
  · rw [sub_zero, abs_of_nonneg (norm_nonneg _), ← Finset.sum_sub_distrib]
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i _ => ?_)
    have h := norm_mul_sub_mul_le (w i ^ 2 * cN blk dN N i) (c i * w i ^ 2) (G N ω i) (ℓ i)
    have h1 : |w i ^ 2 * cN blk dN N i - c i * w i ^ 2| = w i ^ 2 * |cN blk dN N i - c i| := by
      rw [show w i ^ 2 * cN blk dN N i - c i * w i ^ 2 = w i ^ 2 * (cN blk dN N i - c i) by ring,
        abs_mul, abs_of_nonneg (sq_nonneg _)]
    have h2 : |c i * w i ^ 2| = c i * w i ^ 2 := abs_of_nonneg (by have := hc i; positivity)
    rw [h1, h2] at h
    refine h.trans ?_
    have hGb := hG N ω i
    have hci := (hc i).le
    calc w i ^ 2 * |cN blk dN N i - c i| * ‖G N ω i‖ + c i * w i ^ 2 * ‖G N ω i - ℓ i‖
        ≤ w i ^ 2 * |cN blk dN N i - c i| * Cg + c i * w i ^ 2 * ‖G N ω i - ℓ i‖ := by gcongr
      _ = w i ^ 2 * (|cN blk dN N i - c i| * Cg + c i * ‖G N ω i - ℓ i‖) := by ring
  · have hsum := tendstoInProb_finset_sum_zero (μ := μ)
      (F := fun i N ω => w i ^ 2 * (|cN blk dN N i - c i| * Cg + c i * ‖G N ω i - ℓ i‖))
      (fun i => ?_)
    · simpa using hsum
    · have h1 : TendstoInProb μ (fun N _ => |cN blk dN N i - c i| * Cg) 0 := by
        refine tendstoInProb_of_tendsto ?_
        have := (((hn i).sub_const (c i)).abs).mul_const Cg
        simpa using this
      have h2 := (hlim i).const_mul (c i)
      simpa using (h1.add h2).const_mul (w i ^ 2)

include hc hw hz hd hq hn hXB in
/-- **The `g`-form at complex `z`**: `gᵀ G_σ g → ∑ c_i w_i² g_i(z)`. -/
theorem tendstoInProb_qformC_gvec :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)
        - ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ) * gC w i z (sGlob c w z)‖) 0 :=
  tendstoInProb_norm_sub_trans (tendstoInProb_qformC_gvec_sub hc hw hz hd hq hn XB hXB)
    (tendstoInProb_blockMeans_sub hc hn
      (G := fun N ω i => gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (XB N ω).2)
      (fun _ _ _ => norm_gsigAvg_le hz _ _ _ _)
      (fun i => tendstoInProb_gsigAvg hc hw hz hd hq hn (fun N ω => (XB N ω).2)
        (fun N => hasLaw_snd_of_pair (hXB N)) i))

include hc hw hz hd hq hn hXB in
/-- **The second-order `g`-form at complex `z`**: `gᵀ G_σ² g → ∑ c_i w_i² g_i'(z)`. -/
theorem tendstoInProb_qform2C_gvec :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)
        - ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ)
            * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹‖) 0 :=
  tendstoInProb_norm_sub_trans (tendstoInProb_qform2C_gvec_sub hc hw hz hd hq hn XB hXB)
    (tendstoInProb_blockMeans_sub hc hn
      (G := fun N ω i => gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (XB N ω).2)
      (fun _ _ _ => norm_gsig2Avg_le hz _ _ _ _)
      (fun i => tendstoInProb_gsig2Avg hc hw hz hd hq hn (fun N ω => (XB N ω).2)
        (fun N => hasLaw_snd_of_pair (hXB N)) i))

include hc hw hz hd hq hn hXB in
/-- **The cross form at complex `z`**: `yᵀ G_σ g → 0` for bounded deterministic `y`. -/
theorem tendstoInProb_cformC_gvec (y : ∀ N, Fin (pN N) → ℝ) {Ky : ℝ}
    (hy : ∀ N, y N ⬝ᵥ y N ≤ Ky) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z (y N)
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)‖) 0 := by
  have hdR : Tendsto (fun N => ((dN N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have hS0 : ∀ N, 0 ≤ Scalars.wSqMax w / dN N :=
    fun N => div_nonneg (wSqMax_pos hw).le (Nat.cast_nonneg _)
  have hKl : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      nsq (Gsig (tauOf w (blk N)) (dN N) z B *ᵥ R4C.cvec (y N)) ≤ Ky / z.im ^ 2 :=
    fun N B => (nsq_resolvC_mulVec_le (isHermitian_Wsig _ _ B) hz _).trans
      (by rw [nsq_cvec]; exact div_le_div_of_nonneg_right (hy N) (by positivity))
  refine tendstoInProb_of_integral_sq_le_pair XB hXB
    (fun N x B => ‖R4C.cformC (Wsig (tauOf w (blk N)) (dN N) B) z (y N)
      (gvec (tauOf w (blk N)) (dN N) x)‖)
    ?_ (fun _ _ _ => norm_nonneg _) ?_
    (K := fun N => Scalars.wSqMax w / dN N * (Ky / z.im ^ 2)) ?_ ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    exact (measurable_cformC_gvec_pair hz hpN hdN _ (y N)).norm
  · filter_upwards with N B
    simp only [gvec_eq_scaled]
    exact (integral_normSq_lin_scaled _ (transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne')
      (y N) (hKl N B) _ (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N)) (hS0 N)).1
  · have := ((tendsto_const_nhds (x := Scalars.wSqMax w)).div_atTop hdR).mul_const (Ky / z.im ^ 2)
    simpa using this
  · filter_upwards with N B
    simp only [gvec_eq_scaled]
    exact (integral_normSq_lin_scaled _ (transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne')
      (y N) (hKl N B) _ (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N)) (hS0 N)).2

include hc hw hz hd hq hn hXB in
/-- **The second-order cross form at complex `z`**: `yᵀ G_σ² g → 0`. -/
theorem tendstoInProb_cform2C_gvec (y : ∀ N, Fin (pN N) → ℝ) {Ky : ℝ}
    (hy : ∀ N, y N ⬝ᵥ y N ≤ Ky) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (Wsig (tauOf w (blk N)) (dN N) (XB N ω).2) z (y N)
          (gvec (tauOf w (blk N)) (dN N) (XB N ω).1)‖) 0 := by
  have hdR : Tendsto (fun N => ((dN N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have hS0 : ∀ N, 0 ≤ Scalars.wSqMax w / dN N :=
    fun N => div_nonneg (wSqMax_pos hw).le (Nat.cast_nonneg _)
  have hKt : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      (Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B)ᵀ
        = Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B := by
    intro N B
    rw [Matrix.transpose_mul, Gsig, transpose_resolvC (isHermitian_Wsig _ _ B) hz.ne']
  have hKl : ∀ N (B : Matrix (Fin (pN N)) (Fin (qN N)) ℝ),
      nsq ((Gsig (tauOf w (blk N)) (dN N) z B * Gsig (tauOf w (blk N)) (dN N) z B)
        *ᵥ R4C.cvec (y N)) ≤ Ky / z.im ^ 4 :=
    fun N B => (nsq_resolvC_sq_mulVec_le (isHermitian_Wsig _ _ B) hz _).trans
      (by rw [nsq_cvec]; exact div_le_div_of_nonneg_right (hy N) (by positivity))
  refine tendstoInProb_of_integral_sq_le_pair XB hXB
    (fun N x B => ‖R4C.cform2C (Wsig (tauOf w (blk N)) (dN N) B) z (y N)
      (gvec (tauOf w (blk N)) (dN N) x)‖)
    ?_ (fun _ _ _ => norm_nonneg _) ?_
    (K := fun N => Scalars.wSqMax w / dN N * (Ky / z.im ^ 4)) ?_ ?_
  · filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
    obtain ⟨hdN, hqN, -, hpN⟩ := hN
    exact (measurable_cform2C_gvec_pair hz hpN hdN _ (y N)).norm
  · filter_upwards with N B
    simp only [gvec_eq_scaled]
    exact (integral_normSq_lin_scaled _ (hKt N B) (y N) (hKl N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N)) (hS0 N)).1
  · have := ((tendsto_const_nhds (x := Scalars.wSqMax w)).div_atTop hdR).mul_const (Ky / z.im ^ 4)
    simpa using this
  · filter_upwards with N B
    simp only [gvec_eq_scaled]
    exact (integral_normSq_lin_scaled _ (hKt N B) (y N) (hKl N B) _
      (gvec_coef_sq_le _ (tauOf_sq_le w (blk N)) (dN N)) (hS0 N)).2

include hw hd hn hXB in
omit hc hq in
/-- **The squared norm of `g`**: `gᵀ g → ∑ w_i² c_i` (chi-square law of large numbers). -/
theorem tendstoInProb_dotProduct_gvec :
    TendstoInProb μ (fun N ω => gvec (tauOf w (blk N)) (dN N) (XB N ω).1
      ⬝ᵥ gvec (tauOf w (blk N)) (dN N) (XB N ω).1) (∑ i, w i ^ 2 * c i) := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  have hid : ∀ N (x : Fin (pN N) → ℝ),
      gvec (tauOf w (blk N)) (dN N) x ⬝ᵥ gvec (tauOf w (blk N)) (dN N) x
        - ∑ i, w i ^ 2 * cN blk dN N i
      = ∑ j, (tauOf w (blk N) j ^ 2 / dN N) * ((fun t : ℝ => t ^ 2 - 1) (x j)) := by
    intro N x
    have hsq : ((Real.sqrt (dN N))⁻¹) ^ 2 = 1 / dN N := by
      rw [inv_pow, Real.sq_sqrt (Nat.cast_nonneg _), one_div]
    have h1 : gvec (tauOf w (blk N)) (dN N) x ⬝ᵥ gvec (tauOf w (blk N)) (dN N) x
        = ∑ j, (tauOf w (blk N) j ^ 2 / dN N) * x j ^ 2 := by
      simp only [dotProduct, gvec]
      refine Finset.sum_congr rfl fun j _ => ?_
      calc tauOf w (blk N) j * ((Real.sqrt (dN N))⁻¹ * x j)
            * (tauOf w (blk N) j * ((Real.sqrt (dN N))⁻¹ * x j))
          = tauOf w (blk N) j ^ 2 * ((Real.sqrt (dN N))⁻¹) ^ 2 * x j ^ 2 := by ring
        _ = tauOf w (blk N) j ^ 2 / dN N * x j ^ 2 := by rw [hsq]; ring
    have h2 : ∑ i, w i ^ 2 * cN blk dN N i = ∑ j, tauOf w (blk N) j ^ 2 / dN N := by
      rw [← sum_blockSet_real (blk N)]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [Finset.sum_congr rfl fun j hj => by rw [tauOf_sq_eq w (blk N) i hj], Finset.sum_const,
        nsmul_eq_mul, cN]
      ring
    rw [h1, h2, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hbase : TendstoInProb μ (fun N ω => |gvec (tauOf w (blk N)) (dN N) (XB N ω).1
      ⬝ᵥ gvec (tauOf w (blk N)) (dN N) (XB N ω).1 - ∑ i, w i ^ 2 * cN blk dN N i|) 0 := by
    refine tendstoInProb_of_integral_sq_le_pair XB hXB
      (fun N x _ => |gvec (tauOf w (blk N)) (dN N) x ⬝ᵥ gvec (tauOf w (blk N)) (dN N) x
        - ∑ i, w i ^ 2 * cN blk dN N i|)
      (Filter.Eventually.of_forall fun N => ?_) (fun _ _ _ => abs_nonneg _)
      (Filter.Eventually.of_forall fun N _ => ?_)
      (K := fun N => 2 * varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * 1) ^ 2))
      (tendsto_varBound hd hn _ _) (Filter.Eventually.of_forall fun N _ => ?_)
    · refine (Measurable.sub ?_ measurable_const).abs
      exact Finset.measurable_sum _ fun j _ =>
        ((measurable_gvec_apply _ _ j).comp measurable_fst).mul
          ((measurable_gvec_apply _ _ j).comp measurable_fst)
    · simp only [sq_abs, hid]
      exact centered_sq_sub_one.integrable_sq_sum _
    · simp only [sq_abs, hid]
      rw [centered_sq_sub_one.integral_sq_sum]
      have hv := varSq_nonneg
      have hr : ∀ j, (tauOf w (blk N) j ^ 2 / dN N) ^ 2 ≤ (Scalars.wSqMax w / dN N) ^ 2 :=
        fun j => pow_le_pow_left₀ (by positivity)
          (div_le_div_of_nonneg_right (tauOf_sq_le w (blk N) j) (Nat.cast_nonneg _)) 2
      calc (∫ t, ((fun t : ℝ => t ^ 2 - 1) t) ^ 2 ∂(gaussianReal 0 1))
            * ∑ j, (tauOf w (blk N) j ^ 2 / dN N) ^ 2
          ≤ varSq * ∑ _j : Fin (pN N), (Scalars.wSqMax w / dN N) ^ 2 :=
            mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun j _ => hr j) hv
        _ = varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * 1) ^ 2) := by simp
        _ ≤ 2 * varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * 1) ^ 2) := by
            have : 0 ≤ varSq * ((pN N : ℝ) * (Scalars.wSqMax w / dN N * 1) ^ 2) := by
              positivity
            linarith
  have hm : Tendsto (fun N => ∑ i, w i ^ 2 * cN blk dN N i) atTop (𝓝 (∑ i, w i ^ 2 * c i)) :=
    tendsto_finsetSum _ fun i _ => (hn i).const_mul _
  have hm0 : TendstoInProb μ (fun N _ => |∑ i, w i ^ 2 * cN blk dN N i - ∑ i, w i ^ 2 * c i|)
      0 := by
    refine tendstoInProb_of_tendsto ?_
    have := (hm.sub_const (∑ i, w i ^ 2 * c i)).abs
    simpa using this
  refine TendstoInProb.of_le (g := fun N ω => |gvec (tauOf w (blk N)) (dN N) (XB N ω).1
      ⬝ᵥ gvec (tauOf w (blk N)) (dN N) (XB N ω).1 - ∑ i, w i ^ 2 * cN blk dN N i|
      + |∑ i, w i ^ 2 * cN blk dN N i - ∑ i, w i ^ 2 * c i|)
    (fun N => Filter.Eventually.of_forall fun ω => ?_) (by simpa using hbase.add hm0)
  exact (abs_sub_le _ _ _)

end GaussSeq


/-! ### From complex `z` to the real axis: the T.Gen wrapper with the norm rescaling -/

section RealAxis

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {dN : ℕ → ℕ}

theorem cformC_smul_smul {n : ℕ} (W : Matrix (Fin n) (Fin n) ℝ) (z : ℂ) (s : ℝ)
    (x y : Fin n → ℝ) :
    R4C.cformC W z (s • x) (s • y) = ((s : ℂ) ^ 2) * R4C.cformC W z x y := by
  rw [R4C.cformC, R4C.cformC, R2.cvec_smul, R2.cvec_smul, Matrix.mulVec_smul, dotProduct_smul,
    smul_dotProduct, smul_eq_mul, smul_eq_mul, sq]
  ring

theorem cform2C_smul_smul {n : ℕ} (W : Matrix (Fin n) (Fin n) ℝ) (z : ℂ) (s : ℝ)
    (x y : Fin n → ℝ) :
    R4C.cform2C W z (s • x) (s • y) = ((s : ℂ) ^ 2) * R4C.cform2C W z x y := by
  rw [R4C.cform2C, R4C.cform2C, R2.cvec_smul, R2.cvec_smul, Matrix.mulVec_smul, dotProduct_smul,
    smul_dotProduct, smul_eq_mul, smul_eq_mul, sq]
  ring

theorem cform_smul_smul {n : ℕ} (W : Matrix (Fin n) (Fin n) ℝ) (x : ℝ) (s : ℝ)
    (v y : Fin n → ℝ) :
    R4.cform W x (s • v) (s • y) = s ^ 2 * R4.cform W x v y := by
  rw [R4.cform, R4.cform, Matrix.mulVec_smul, dotProduct_smul, smul_dotProduct, smul_eq_mul,
    smul_eq_mul, sq]
  ring

theorem cform2_smul_smul {n : ℕ} (W : Matrix (Fin n) (Fin n) ℝ) (x : ℝ) (s : ℝ)
    (v y : Fin n → ℝ) :
    R4.cform2 W x (s • v) (s • y) = s ^ 2 * R4.cform2 W x v y := by
  rw [R4.cform2, R4.cform2, Matrix.mulVec_smul, dotProduct_smul, smul_dotProduct, smul_eq_mul,
    smul_eq_mul, sq]
  ring

/-- The scale `s = 1/(1 + Kp + Kq)` puts both vectors inside the ball of radius `2`. -/
theorem sq_scale_mul_le {Kp Kq : ℝ} (hKp : 0 ≤ Kp) (hKq : 0 ≤ Kq) {K : ℝ} (hK : K ≤ Kp + Kq) :
    (1 / (1 + Kp + Kq)) ^ 2 * K ≤ 4 := by
  have hpos : 0 < 1 + Kp + Kq := by linarith
  rw [div_pow, one_pow, div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
  nlinarith

/-- The good event of T.Gen after rescaling by `s`, as an intersection. -/
theorem scaled_norm_event_eq {n : ℕ} (p q : Fin n → ℝ) {s : ℝ} (hs : 0 < s) :
    ((s • p) ⬝ᵥ (s • p) ≤ 4 ∧ (s • q) ⬝ᵥ (s • q) ≤ 4)
      ↔ (p ⬝ᵥ p ≤ 4 / s ^ 2 ∧ q ⬝ᵥ q ≤ 4 / s ^ 2) := by
  simp only [smul_dotProduct, dotProduct_smul, smul_eq_mul, ← mul_assoc, ← sq,
    le_div_iff₀ (pow_pos hs 2)]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨by linarith, by linarith⟩
  · rintro ⟨h1, h2⟩
    exact ⟨by linarith, by linarith⟩

/-- **T.Gen with the norm rescaling**, first order. `p` and `q` lie in the balls of radius
`√Kp`, `√Kq` with probability tending to `1`; the scale `s = 1/(1 + Kp + Kq)` puts both
inside the ball of radius `2` that `RMT/T.lean` hardwires. -/
theorem tendstoInProb_cform_of_complex_scaled [∀ N, IsProbabilityMeasure (μ N)]
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ)) {b : ℝ}
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}) atTop (𝓝 1))
    {Kp Kq : ℝ} (hKp : 0 ≤ Kp) (hKq : 0 ≤ Kq)
    (hpMeas : ∀ N r, NullMeasurableSet {ω | p N ω ⬝ᵥ p N ω ≤ r} (μ N))
    (hp : Tendsto (fun N => μ N {ω | p N ω ⬝ᵥ p N ω ≤ Kp}) atTop (𝓝 1))
    (hqMeas : ∀ N r, NullMeasurableSet {ω | q N ω ⬝ᵥ q N ω ≤ r} (μ N))
    (hq : Tendsto (fun N => μ N {ω | q N ω ⬝ᵥ q N ω ≤ Kq}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cformC (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform (W₀ N ω) x (p N ω) (q N ω)) ℓ := by
  set s : ℝ := 1 / (1 + Kp + Kq) with hsdef
  have hs : 0 < s := by rw [hsdef]; positivity
  have hs2 : s ^ 2 ≠ 0 := by positivity
  have hnormMeas : ∀ N, NullMeasurableSet
      {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4} (μ N) := by
    intro N
    have heq : {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4}
        = {ω | p N ω ⬝ᵥ p N ω ≤ 4 / s ^ 2} ∩ {ω | q N ω ⬝ᵥ q N ω ≤ 4 / s ^ 2} := by
      ext ω
      simp only [Set.mem_inter_iff, Set.mem_ofPred_eq]
      exact scaled_norm_event_eq _ _ hs
    rw [heq]
    exact (hpMeas N _).inter (hqMeas N _)
  have hnorm : Tendsto (fun N => μ N
      {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4}) atTop (𝓝 1) := by
    refine tendsto_measure_one_of_bad
      (s := fun N => {ω | p N ω ⬝ᵥ p N ω ≤ Kp}ᶜ ∪ {ω | q N ω ⬝ᵥ q N ω ≤ Kq}ᶜ)
      (fun N ω hω => ?_) (tendsto_measure_zero_union
        (tendsto_measure_compl_zero (fun N => hpMeas N Kp) hp)
        (tendsto_measure_compl_zero (fun N => hqMeas N Kq) hq))
    by_contra hcon
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_ofPred_eq, not_or, not_not] at hcon
    refine hω ?_
    change (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4
    rw [scaled_norm_event_eq _ _ hs, le_div_iff₀ (pow_pos hs 2), le_div_iff₀ (pow_pos hs 2)]
    have h1 := sq_scale_mul_le hKp hKq (K := Kp) (by linarith)
    have h2 := sq_scale_mul_le hKp hKq (K := Kq) (by linarith)
    rw [← hsdef] at h1 h2
    have hpp : 0 ≤ p N ω ⬝ᵥ p N ω := Finset.sum_nonneg fun _ _ => mul_self_nonneg _
    have hqq : 0 ≤ q N ω ⬝ᵥ q N ω := Finset.sum_nonneg fun _ _ => mul_self_nonneg _
    exact ⟨by nlinarith [sq_nonneg s, hcon.1], by nlinarith [sq_nonneg s, hcon.2]⟩
  have hcplx' : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W₀ N ω) z (s • p N ω) (s • q N ω) - ((s : ℂ) ^ 2) * L z‖) 0 := by
    intro z hz
    have h := (hcplx z hz).const_mul (s ^ 2)
    rw [mul_zero] at h
    refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) h
    rw [cformC_smul_smul, ← mul_sub, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hs, sub_zero, abs_of_nonneg (by positivity)]
  have hL' : Tendsto (fun η : ℝ => ((s : ℂ) ^ 2) * L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((s ^ 2 * ℓ : ℝ) : ℂ)) := by
    have := hL.const_mul ((s : ℂ) ^ 2)
    push_cast
    exact this
  have h1 := T.Gen.tendstoInProb_cform_of_complex W₀ hsymm (fun N ω => s • p N ω)
    (fun N ω => s • q N ω) hedgeMeas hedge hnormMeas hnorm (fun z => ((s : ℂ) ^ 2) * L z)
    (s ^ 2 * ℓ) hcplx' hx hL'
  have h2 := h1.const_mul (s ^ 2)⁻¹
  simp only [cform_smul_smul, inv_mul_cancel_left₀ hs2] at h2
  exact h2

/-- **T.Gen with the norm rescaling**, second order. -/
theorem tendstoInProb_cform2_of_complex_scaled [∀ N, IsProbabilityMeasure (μ N)]
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ)) {b : ℝ}
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}) atTop (𝓝 1))
    {Kp Kq : ℝ} (hKp : 0 ≤ Kp) (hKq : 0 ≤ Kq)
    (hpMeas : ∀ N r, NullMeasurableSet {ω | p N ω ⬝ᵥ p N ω ≤ r} (μ N))
    (hp : Tendsto (fun N => μ N {ω | p N ω ⬝ᵥ p N ω ≤ Kp}) atTop (𝓝 1))
    (hqMeas : ∀ N r, NullMeasurableSet {ω | q N ω ⬝ᵥ q N ω ≤ r} (μ N))
    (hq : Tendsto (fun N => μ N {ω | q N ω ⬝ᵥ q N ω ≤ Kq}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cform2C (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform2 (W₀ N ω) x (p N ω) (q N ω)) ℓ := by
  set s : ℝ := 1 / (1 + Kp + Kq) with hsdef
  have hs : 0 < s := by rw [hsdef]; positivity
  have hs2 : s ^ 2 ≠ 0 := by positivity
  have hnormMeas : ∀ N, NullMeasurableSet
      {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4} (μ N) := by
    intro N
    have heq : {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4}
        = {ω | p N ω ⬝ᵥ p N ω ≤ 4 / s ^ 2} ∩ {ω | q N ω ⬝ᵥ q N ω ≤ 4 / s ^ 2} := by
      ext ω
      simp only [Set.mem_inter_iff, Set.mem_ofPred_eq]
      exact scaled_norm_event_eq _ _ hs
    rw [heq]
    exact (hpMeas N _).inter (hqMeas N _)
  have hnorm : Tendsto (fun N => μ N
      {ω | (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4}) atTop (𝓝 1) := by
    refine tendsto_measure_one_of_bad
      (s := fun N => {ω | p N ω ⬝ᵥ p N ω ≤ Kp}ᶜ ∪ {ω | q N ω ⬝ᵥ q N ω ≤ Kq}ᶜ)
      (fun N ω hω => ?_) (tendsto_measure_zero_union
        (tendsto_measure_compl_zero (fun N => hpMeas N Kp) hp)
        (tendsto_measure_compl_zero (fun N => hqMeas N Kq) hq))
    by_contra hcon
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_ofPred_eq, not_or, not_not] at hcon
    refine hω ?_
    change (s • p N ω) ⬝ᵥ (s • p N ω) ≤ 4 ∧ (s • q N ω) ⬝ᵥ (s • q N ω) ≤ 4
    rw [scaled_norm_event_eq _ _ hs, le_div_iff₀ (pow_pos hs 2), le_div_iff₀ (pow_pos hs 2)]
    have h1 := sq_scale_mul_le hKp hKq (K := Kp) (by linarith)
    have h2 := sq_scale_mul_le hKp hKq (K := Kq) (by linarith)
    rw [← hsdef] at h1 h2
    have hpp : 0 ≤ p N ω ⬝ᵥ p N ω := Finset.sum_nonneg fun _ _ => mul_self_nonneg _
    have hqq : 0 ≤ q N ω ⬝ᵥ q N ω := Finset.sum_nonneg fun _ _ => mul_self_nonneg _
    exact ⟨by nlinarith [sq_nonneg s, hcon.1], by nlinarith [sq_nonneg s, hcon.2]⟩
  have hcplx' : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W₀ N ω) z (s • p N ω) (s • q N ω) - ((s : ℂ) ^ 2) * L z‖) 0 := by
    intro z hz
    have h := (hcplx z hz).const_mul (s ^ 2)
    rw [mul_zero] at h
    refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) h
    rw [cform2C_smul_smul, ← mul_sub, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hs, sub_zero, abs_of_nonneg (by positivity)]
  have hL' : Tendsto (fun η : ℝ => ((s : ℂ) ^ 2) * L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((s ^ 2 * ℓ : ℝ) : ℂ)) := by
    have := hL.const_mul ((s : ℂ) ^ 2)
    push_cast
    exact this
  have h1 := T.Gen.tendstoInProb_cform2_of_complex W₀ hsymm (fun N ω => s • p N ω)
    (fun N ω => s • q N ω) hedgeMeas hedge hnormMeas hnorm (fun z => ((s : ℂ) ^ 2) * L z)
    (s ^ 2 * ℓ) hcplx' hx hL'
  have h2 := h1.const_mul (s ^ 2)⁻¹
  simp only [cform2_smul_smul, inv_mul_cancel_left₀ hs2] at h2
  exact h2

end RealAxis

end HetR2

/-! ### The model: the six forms of the column split at complex `z` -/

namespace MultiTableModel

open HetR2 HetStein HetR1 MPhet

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

theorem d_eq_pred_succ (m : MultiTableModel μ M n d) (N : ℕ) : d N = (d N - 1) + 1 :=
  (Nat.succ_pred_eq_of_pos (m.stack.hd N)).symm

/-- The block `B` of `exists_block_hasLaw_het`, chosen once per `N`. -/
noncomputable def blockB (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise) (N : ℕ) :
    Ω N → Matrix (Fin (∑ i, n i N)) (Fin (d N - 1)) ℝ :=
  Classical.choose (m.exists_block_hasLaw_het N hG (m.d_eq_pred_succ N))

theorem blockB_gram (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise) (N : ℕ)
    (ω : Ω N) :
    m.EperpHet N ω * (m.EperpHet N ω)ᵀ
      = ((d N : ℝ))⁻¹ • (m.blockB hG N ω * (m.blockB hG N ω)ᵀ) :=
  (Classical.choose_spec (m.exists_block_hasLaw_het N hG (m.d_eq_pred_succ N))).1 ω

/-- The pair `(Z v, B)`. -/
noncomputable def pairXB (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise) (N : ℕ)
    (ω : Ω N) : PairSpace (∑ i, n i N) (d N - 1) :=
  (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N), m.blockB hG N ω)

theorem hasLaw_pairXB (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise) (N : ℕ) :
    HasLaw (m.pairXB hG N) (pairLaw (∑ i, n i N) (d N - 1)) (μ N) :=
  (Classical.choose_spec (m.exists_block_hasLaw_het N hG (m.d_eq_pred_succ N))).2

theorem hasLaw_blockB (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise) (N : ℕ) :
    HasLaw (m.blockB hG N) (gaussianMatrix (∑ i, n i N) (d N - 1)) (μ N) :=
  hasLaw_snd_of_pair (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1) (m.hasLaw_pairXB hG N)

/-- `W₀' = W_σ(B)` on the chosen block. -/
theorem W0het_eq_Wsig_blockB (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.W0het w N ω = Wsig (tauOf w (blkStack n N)) (d N) (m.blockB hG N ω) :=
  HetR1.W0het_eq_Wsig m w N ω _ (m.blockB_gram hG N ω)

/-- `Σ^{1/2} e = gvec τ d (Z v)`. -/
theorem SigmaHalf_mulVec_eHet (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.SigmaHalf w N *ᵥ m.eHet N ω
      = gvec (tauOf w (blkStack n N)) (d N) (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N)) := by
  have he : m.eHet N ω
      = (Real.sqrt (d N))⁻¹ • (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N)) := by
    change ((Real.sqrt (d N))⁻¹ • m.stack.Z N ω) *ᵥ WithLp.ofLp (m.stack.v N) = _
    rw [Matrix.smul_mulVec]
  rw [he, SigmaHalf_eq_diagonal]
  funext j
  rw [Matrix.mulVec_diagonal]
  rfl

omit [NeZero M] in
/-- Sums over a block of the stacked index, through `finSigmaFinEquiv`. -/
theorem sum_blockSet_blkStack (N : ℕ) (i : Fin M) (f : Fin (∑ i, n i N) → ℝ) :
    ∑ k ∈ blockSet (blkStack n N) i, f k = ∑ r : Fin (n i N), f (finSigmaFinEquiv ⟨i, r⟩) := by
  classical
  rw [blockSet, Finset.sum_filter]
  have h := Fintype.sum_equiv finSigmaFinEquiv
    (fun σ : Σ i', Fin (n i' N) => if σ.1 = i then f (finSigmaFinEquiv σ) else 0)
    (fun k => if blkStack n N k = i then f k else 0) (fun σ => by simp [blkStack])
  rw [← h, Fintype.sum_sigma]
  rw [Finset.sum_eq_single i (fun x _ hx => by simp [hx]) (fun h => absurd (Finset.mem_univ i) h)]
  simp

/-! ### The model: block norms of `ũ₀`, the regime parameters, and the boundary values -/

omit [NeZero M] in
theorem stackThetaSqW_nonneg (m : MultiTableModel μ M n d) (w : Fin M → ℝ) :
    0 ≤ m.stackThetaSqW w :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

/-- The block norms of `ũ₀`: `∑_{k ∈ J_i} (ũ₀)_k² = (θ_i w_i)²`. -/
theorem sum_blockSet_u0Het_sq (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (i : Fin M) :
    ∑ k ∈ blockSet (blkStack n N) i, m.u0Het w N k ^ 2 = ((m.tbl i).θ * w i) ^ 2 := by
  rw [sum_blockSet_blkStack]
  have hu : ∑ r : Fin (n i N), (WithLp.ofLp ((m.tbl i).u N) r) ^ 2 = 1 := by
    have h := R2.dotProduct_ofLp_self ((m.tbl i).hu N)
    simpa [dotProduct, sq] using h
  by_cases h0 : m.stackThetaSqW w = 0
  · have hθ : ((m.tbl i).θ * w i) ^ 2 = 0 :=
      (Finset.sum_eq_zero_iff_of_nonneg fun j _ => sq_nonneg ((m.tbl j).θ * w j)).mp h0 i
        (Finset.mem_univ i)
    have hz : m.u0Het w N = 0 := by
      change m.stackThetaW w • WithLp.ofLp (m.stackUW w N) = 0
      rw [stackThetaW, h0, Real.sqrt_zero, zero_smul]
    simp [hz, hθ]
  · have hθW : m.stackThetaW w ≠ 0 := by
      rw [stackThetaW]
      exact Real.sqrt_ne_zero'.mpr
        (lt_of_le_of_ne (m.stackThetaSqW_nonneg w) (Ne.symm h0))
    have hk : ∀ r : Fin (n i N), m.u0Het w N (finSigmaFinEquiv ⟨i, r⟩)
        = (m.tbl i).θ * w i * WithLp.ofLp ((m.tbl i).u N) r := by
      intro r
      change m.stackThetaW w * WithLp.ofLp (m.stackUW w N) (finSigmaFinEquiv ⟨i, r⟩) = _
      rw [stackUW, if_neg h0, WithLp.ofLp_toLp]
      change m.stackThetaW w * m.stackUSigmaW w N (finSigmaFinEquiv.symm (finSigmaFinEquiv ⟨i, r⟩))
        = _
      rw [Equiv.symm_apply_apply]
      simp only [stackUSigmaW]
      field_simp
    simp only [hk, mul_pow, ← Finset.mul_sum, hu, mul_one]

theorem dotProduct_u0Het_self (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    m.u0Het w N ⬝ᵥ m.u0Het w N = m.stackThetaSqW w := by
  have h : m.u0Het w N ⬝ᵥ m.u0Het w N = ∑ k, m.u0Het w N k ^ 2 := by
    simp only [dotProduct, sq]
  rw [h, ← sum_blockSet_real (blkStack n N), stackThetaSqW]
  exact Finset.sum_congr rfl fun i _ => m.sum_blockSet_u0Het_sq w N i

omit [NeZero M] in
theorem tendsto_cN_blkStack (m : MultiTableModel μ M n d) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (i : Fin M) :
    Tendsto (fun N => cN (fun N => blkStack n N) d N i) atTop (𝓝 (c i)) := by
  have h : (fun N => cN (fun N => blkStack n N) d N i) = fun N => (n i N : ℝ) / d N := by
    funext N
    simp only [cN, card_blockSet_blkStack]
  rw [h]
  exact (hreg i).2.2

theorem tendsto_pred_div {d : ℕ → ℕ} (hd : Tendsto d atTop atTop) :
    Tendsto (fun N => ((d N - 1 : ℕ) : ℝ) / d N) atTop (𝓝 1) := by
  have hdR : Tendsto (fun N => ((d N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have h1 : Tendsto (fun N => (1 : ℝ) - 1 / (d N : ℝ)) atTop (𝓝 (1 - 0)) :=
    tendsto_const_nhds.sub (tendsto_const_nhds.div_atTop hdR)
  rw [sub_zero] at h1
  refine h1.congr' ?_
  filter_upwards [hd.eventually_ge_atTop 1] with N hN
  have hne : (d N : ℝ) ≠ 0 := by exact_mod_cast (Nat.one_le_iff_ne_zero.mp hN)
  rw [Nat.cast_pred (Nat.lt_of_lt_of_le Nat.zero_lt_one hN)]
  field_simp

/-- `Φ` at complex `z`: `∑ (θ_i w_i)² g_i(z)`. -/
noncomputable def PhiC (θ c w : Fin M → ℝ) (z : ℂ) : ℂ :=
  ∑ i, (((θ i * w i) ^ 2 : ℝ) : ℂ) * gC w i z (sGlob c w z)

/-- `Ψ` at complex `z`: `∑ c_i w_i² g_i(z)`. -/
noncomputable def PsiC (c w : Fin M → ℝ) (z : ℂ) : ℂ :=
  ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ) * gC w i z (sGlob c w z)

/-- `Φ'` at complex `z`. -/
noncomputable def PhiDerivC (θ c w : Fin M → ℝ) (z : ℂ) : ℂ :=
  ∑ i, (((θ i * w i) ^ 2 : ℝ) : ℂ)
    * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹

/-- `Ψ'` at complex `z`. -/
noncomputable def PsiDerivC (c w : Fin M → ℝ) (z : ℂ) : ℂ :=
  ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ)
    * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹

omit [NeZero M] in
theorem tendsto_PhiC {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) :
    Tendsto (fun η : ℝ => PhiC θ c w ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((Phihet θ c w x : ℝ) : ℂ)) := by
  have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gC_sGlob hc hw hx i).const_mul (((θ i * w i) ^ 2 : ℝ) : ℂ)
  have hlim : ∑ i, (((θ i * w i) ^ 2 : ℝ) : ℂ) * ((ghet c w i x : ℝ) : ℂ)
      = ((Phihet θ c w x : ℝ) : ℂ) := by
    simp only [Phihet, Complex.ofReal_sum, Complex.ofReal_mul, Complex.ofReal_pow, mul_pow]
  rw [← hlim]
  exact this

omit [NeZero M] in
theorem tendsto_PsiC {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) :
    Tendsto (fun η : ℝ => PsiC c w ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((Psihet c w x : ℝ) : ℂ)) := by
  have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gC_sGlob hc hw hx i).const_mul ((c i * w i ^ 2 : ℝ) : ℂ)
  have hlim : ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ) * ((ghet c w i x : ℝ) : ℂ)
      = ((Psihet c w x : ℝ) : ℂ) := by
    simp only [Psihet, Complex.ofReal_sum, Complex.ofReal_mul, Complex.ofReal_pow]
  rw [← hlim]
  exact this

omit [NeZero M] in
theorem tendsto_PhiDerivC {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) :
    Tendsto (fun η : ℝ => PhiDerivC θ c w ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((PhihetDeriv θ c w x : ℝ) : ℂ)) := by
  have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gCDeriv_sGlob hc hw hx i).const_mul (((θ i * w i) ^ 2 : ℝ) : ℂ)
  have hlim : ∑ i, (((θ i * w i) ^ 2 : ℝ) : ℂ) * ((ghetDeriv c w i x : ℝ) : ℂ)
      = ((PhihetDeriv θ c w x : ℝ) : ℂ) := by
    simp only [PhihetDeriv, Complex.ofReal_sum, Complex.ofReal_mul, Complex.ofReal_pow, mul_pow]
  rw [← hlim]
  exact this

omit [NeZero M] in
theorem tendsto_PsiDerivC {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) :
    Tendsto (fun η : ℝ => PsiDerivC c w ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((PsihetDeriv c w x : ℝ) : ℂ)) := by
  have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gCDeriv_sGlob hc hw hx i).const_mul ((c i * w i ^ 2 : ℝ) : ℂ)
  have hlim : ∑ i, ((c i * w i ^ 2 : ℝ) : ℂ) * ((ghetDeriv c w i x : ℝ) : ℂ)
      = ((PsihetDeriv c w x : ℝ) : ℂ) := by
    simp only [PsihetDeriv, Complex.ofReal_sum, Complex.ofReal_mul, Complex.ofReal_pow]
  rw [← hlim]
  exact this

/-! ### The six forms at complex `z` in the model -/

section ComplexForms

variable (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
  (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {z : ℂ} (hz : 0 < z.im)

include hc hw hreg hG hz in
/-- `ũ₀ᵀ G₀'(z) ũ₀ → Φ(z)` at complex `z`. -/
theorem tendstoInProb_qformC_u0Het :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0het w N ω) z (m.u0Het w N)
      - PhiC (fun i => (m.tbl i).θ) c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w]
  exact tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.blockB hG) (m.hasLaw_blockB hG) (fun N => m.u0Het w N)
    (fun i => ((m.tbl i).θ * w i) ^ 2) (fun N i => m.sum_blockSet_u0Het_sq w N i)

include hc hw hreg hG hz in
/-- `ũ₀ᵀ G₀'(z)² ũ₀ → Φ'(z)` at complex `z`. -/
theorem tendstoInProb_qform2C_u0Het :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0het w N ω) z (m.u0Het w N)
      - PhiDerivC (fun i => (m.tbl i).θ) c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w]
  exact tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.blockB hG) (m.hasLaw_blockB hG) (fun N => m.u0Het w N)
    (fun i => ((m.tbl i).θ * w i) ^ 2) (fun N i => m.sum_blockSet_u0Het_sq w N i)

include hc hw hreg hG hz in
/-- `(Σ^{1/2} e)ᵀ G₀'(z) (Σ^{1/2} e) → Ψ(z)` at complex `z`. -/
theorem tendstoInProb_qformC_eHet :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω)
      - PsiC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w, m.SigmaHalf_mulVec_eHet w]
  exact tendstoInProb_qformC_gvec (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.pairXB hG) (m.hasLaw_pairXB hG)

include hc hw hreg hG hz in
/-- `(Σ^{1/2} e)ᵀ G₀'(z)² (Σ^{1/2} e) → Ψ'(z)` at complex `z`. -/
theorem tendstoInProb_qform2C_eHet :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω)
      - PsiDerivC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w, m.SigmaHalf_mulVec_eHet w]
  exact tendstoInProb_qform2C_gvec (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.pairXB hG) (m.hasLaw_pairXB hG)

include hc hw hreg hG hz in
/-- `ũ₀ᵀ G₀'(z) (Σ^{1/2} e) → 0` at complex `z`. -/
theorem tendstoInProb_cformC_u0Het_eHet :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0het w N ω) z (m.u0Het w N)
      (m.SigmaHalf w N *ᵥ m.eHet N ω)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w, m.SigmaHalf_mulVec_eHet w]
  exact tendstoInProb_cformC_gvec (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.pairXB hG) (m.hasLaw_pairXB hG) (fun N => m.u0Het w N)
    (Ky := m.stackThetaSqW w) (fun N => (m.dotProduct_u0Het_self w N).le)

include hc hw hreg hG hz in
/-- `ũ₀ᵀ G₀'(z)² (Σ^{1/2} e) → 0` at complex `z`. -/
theorem tendstoInProb_cform2C_u0Het_eHet :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0het w N ω) z (m.u0Het w N)
      (m.SigmaHalf w N *ᵥ m.eHet N ω)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0het_eq_Wsig_blockB hG w, m.SigmaHalf_mulVec_eHet w]
  exact tendstoInProb_cform2C_gvec (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_pred_div hd)
    (m.tendsto_cN_blkStack hreg) (m.pairXB hG) (m.hasLaw_pairXB hG) (fun N => m.u0Het w N)
    (Ky := m.stackThetaSqW w) (fun N => (m.dotProduct_u0Het_self w N).le)

include hw hreg hG in
omit hc in
/-- `‖Σ^{1/2} e‖² → ∑ w_i² c_i` in probability. -/
theorem tendstoInProb_SigmaHalf_eHet_norm :
    TendstoInProb μ (fun N ω => (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω))
      (∑ i, w i ^ 2 * c i) := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.SigmaHalf_mulVec_eHet w]
  exact tendstoInProb_dotProduct_gvec (pN := fun N => ∑ i, n i N) (qN := fun N => d N - 1)
    (dN := d) (blk := fun N => blkStack n N) hw hd
    (m.tendsto_cN_blkStack hreg) (m.pairXB hG) (m.hasLaw_pairXB hG)

end ComplexForms

/-! ### The six forms at real `x > b` (the fields of `ResolventLimitsHet`) -/

section RealForms

theorem measurable_eHet_apply (m : MultiTableModel μ M n d) (N : ℕ) (r : Fin (∑ i, n i N)) :
    Measurable fun ω => m.eHet N ω r := by
  have hZ : ∀ (r : Fin (∑ i, n i N)) (j : Fin (d N)), Measurable fun ω => m.stack.Z N ω r j :=
    fun r j => (measurable_pi_apply j).comp ((measurable_pi_apply r).comp (m.stack.hZ N))
  have h : (fun ω => m.eHet N ω r)
      = fun ω => ∑ j, ((Real.sqrt (d N))⁻¹ * m.stack.Z N ω r j) * WithLp.ofLp (m.stack.v N) j :=
    rfl
  rw [h]
  exact Finset.measurable_sum _ fun j _ => ((hZ r j).const_mul _).mul_const _

theorem measurable_SigmaHalf_eHet_dot (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    Measurable fun ω =>
      (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) := by
  have hcoord : ∀ r, Measurable fun ω => (m.SigmaHalf w N *ᵥ m.eHet N ω) r := by
    intro r
    have h : (fun ω => (m.SigmaHalf w N *ᵥ m.eHet N ω) r)
        = fun ω => ∑ k, m.SigmaHalf w N r k * m.eHet N ω k := rfl
    rw [h]
    exact Finset.measurable_sum _ fun k _ => (m.measurable_eHet_apply N k).const_mul _
  exact Finset.measurable_sum _ fun r _ => (hcoord r).mul (hcoord r)

/-- A set cut out by a condition that does not depend on `ω` is measurable. -/
theorem measurableSet_const_prop (N : ℕ) (P : Prop) : MeasurableSet {_ω : Ω N | P} := by
  classical
  by_cases h : P
  · simp only [h, Set.ofPred_true]
    exact MeasurableSet.univ
  · simp only [h, Set.ofPred_false]
    exact MeasurableSet.empty

variable [∀ N, IsProbabilityMeasure (μ N)] (m : MultiTableModel μ M n d) (w c : Fin M → ℝ)
  (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
  (hG : m.JointGaussianNoise) {b : ℝ} (hb : bHet c w ≤ b)
  (hedge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ b + ε}) atTop (𝓝 1))
  {x : ℝ} (hx : b < x)

include hw hreg hG in
omit hc in
/-- The norm event of `Σ^{1/2} e`: `‖Σ^{1/2} e‖² ≤ ∑ w_i² c_i + 1` with probability `→ 1`. -/
theorem tendsto_measure_SigmaHalf_eHet_norm_le :
    Tendsto (fun N => μ N {ω | (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω)
      ≤ ∑ i, w i ^ 2 * c i + 1}) atTop (𝓝 1) := by
  have h := m.tendstoInProb_SigmaHalf_eHet_norm w c hw hreg hG 1 one_pos
  refine tendsto_measure_one_of_bad (s := fun N => {ω | 1 ≤
    |(m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) - ∑ i, w i ^ 2 * c i|})
    (fun N ω hω => ?_) h
  simp only [Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω ⊢
  rw [le_abs]
  left
  linarith

include hc hw hreg hG hb hedge hx in
/-- **Field `uu` at real `x`**: `ũ₀ᵀ G₀'(x) ũ₀ → Φ(x)`. -/
theorem tendstoInProb_qform_u0Het :
    TendstoInProb μ (fun N ω => R4.qform (m.W0het w N ω) x (m.u0Het w N))
      (Phihet (fun i => (m.tbl i).θ) c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have hpMeas : ∀ N r, NullMeasurableSet {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ r} (μ N) :=
    fun N r => (measurableSet_const_prop N _).nullMeasurableSet
  have hp : Tendsto (fun N => μ N {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w})
      atTop (𝓝 1) := by
    have : ∀ N, {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w} = Set.univ :=
      fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_u0Het_self w N).le
    simp only [this, measure_univ]
    exact tendsto_const_nhds
  have h := tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N _ => m.u0Het w N) (fun N _ => m.u0Het w N)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    (m.stackThetaSqW_nonneg w) (m.stackThetaSqW_nonneg w) hpMeas hp hpMeas hp
    (PhiC (fun i => (m.tbl i).θ) c w) (Phihet (fun i => (m.tbl i).θ) c w x)
    (fun z hz => m.tendstoInProb_qformC_u0Het w c hc hw hreg hG hz) hx
    (tendsto_PhiC hc hw hbx)
  simpa only [R4.cform_self] using h

include hc hw hreg hG hb hedge hx in
/-- **Field `uu2` at real `x`**: `ũ₀ᵀ G₀'(x)² ũ₀ → Φ'(x)`. -/
theorem tendstoInProb_qform2_u0Het :
    TendstoInProb μ (fun N ω => R4.qform2 (m.W0het w N ω) x (m.u0Het w N))
      (PhihetDeriv (fun i => (m.tbl i).θ) c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have hpMeas : ∀ N r, NullMeasurableSet {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ r} (μ N) :=
    fun N r => (measurableSet_const_prop N _).nullMeasurableSet
  have hp : Tendsto (fun N => μ N {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w})
      atTop (𝓝 1) := by
    have : ∀ N, {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w} = Set.univ :=
      fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_u0Het_self w N).le
    simp only [this, measure_univ]
    exact tendsto_const_nhds
  have h := tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N _ => m.u0Het w N) (fun N _ => m.u0Het w N)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    (m.stackThetaSqW_nonneg w) (m.stackThetaSqW_nonneg w) hpMeas hp hpMeas hp
    (PhiDerivC (fun i => (m.tbl i).θ) c w) (PhihetDeriv (fun i => (m.tbl i).θ) c w x)
    (fun z hz => m.tendstoInProb_qform2C_u0Het w c hc hw hreg hG hz) hx
    (tendsto_PhiDerivC hc hw hbx)
  simpa only [R4.cform2_self] using h

include hc hw hreg hG hb hedge hx in
/-- **Field `ee` at real `x`**: `(Σ^{1/2} e)ᵀ G₀'(x) (Σ^{1/2} e) → Ψ(x)`. -/
theorem tendstoInProb_qform_eHet :
    TendstoInProb μ (fun N ω => R4.qform (m.W0het w N ω) x (m.SigmaHalf w N *ᵥ m.eHet N ω))
      (Psihet c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have hqMeas : ∀ N r, NullMeasurableSet {ω : Ω N |
      (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) ≤ r} (μ N) :=
    fun N r => (measurableSet_le (m.measurable_SigmaHalf_eHet_dot w N)
      measurable_const).nullMeasurableSet
  have hKq : 0 ≤ ∑ i, w i ^ 2 * c i + 1 := by
    have : 0 ≤ ∑ i, w i ^ 2 * c i := Finset.sum_nonneg fun i _ => by
      have := (hc i).le; positivity
    linarith
  have h := tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω) (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    hKq hKq hqMeas (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG) hqMeas
    (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG)
    (PsiC c w) (Psihet c w x)
    (fun z hz => m.tendstoInProb_qformC_eHet w c hc hw hreg hG hz) hx
    (tendsto_PsiC hc hw hbx)
  simpa only [R4.cform_self] using h

include hc hw hreg hG hb hedge hx in
/-- **Field `ee2` at real `x`**: `(Σ^{1/2} e)ᵀ G₀'(x)² (Σ^{1/2} e) → Ψ'(x)`. -/
theorem tendstoInProb_qform2_eHet :
    TendstoInProb μ (fun N ω => R4.qform2 (m.W0het w N ω) x (m.SigmaHalf w N *ᵥ m.eHet N ω))
      (PsihetDeriv c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have hqMeas : ∀ N r, NullMeasurableSet {ω : Ω N |
      (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) ≤ r} (μ N) :=
    fun N r => (measurableSet_le (m.measurable_SigmaHalf_eHet_dot w N)
      measurable_const).nullMeasurableSet
  have hKq : 0 ≤ ∑ i, w i ^ 2 * c i + 1 := by
    have : 0 ≤ ∑ i, w i ^ 2 * c i := Finset.sum_nonneg fun i _ => by
      have := (hc i).le; positivity
    linarith
  have h := tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω) (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    hKq hKq hqMeas (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG) hqMeas
    (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG)
    (PsiDerivC c w) (PsihetDeriv c w x)
    (fun z hz => m.tendstoInProb_qform2C_eHet w c hc hw hreg hG hz) hx
    (tendsto_PsiDerivC hc hw hbx)
  simpa only [R4.cform2_self] using h

include hc hw hreg hG hedge hx in
omit hb in
/-- **Field `ue` at real `x`**: `ũ₀ᵀ G₀'(x) (Σ^{1/2} e) → 0`. -/
theorem tendstoInProb_cform_u0Het_eHet :
    TendstoInProb μ (fun N ω => R4.cform (m.W0het w N ω) x (m.u0Het w N)
      (m.SigmaHalf w N *ᵥ m.eHet N ω)) 0 := by
  have hpMeas : ∀ N r, NullMeasurableSet {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ r} (μ N) :=
    fun N r => (measurableSet_const_prop N _).nullMeasurableSet
  have hp : Tendsto (fun N => μ N {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w})
      atTop (𝓝 1) := by
    have : ∀ N, {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w} = Set.univ :=
      fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_u0Het_self w N).le
    simp only [this, measure_univ]
    exact tendsto_const_nhds
  have hqMeas : ∀ N r, NullMeasurableSet {ω : Ω N |
      (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) ≤ r} (μ N) :=
    fun N r => (measurableSet_le (m.measurable_SigmaHalf_eHet_dot w N)
      measurable_const).nullMeasurableSet
  have hKq : 0 ≤ ∑ i, w i ^ 2 * c i + 1 := by
    have : 0 ≤ ∑ i, w i ^ 2 * c i := Finset.sum_nonneg fun i _ => by
      have := (hc i).le; positivity
    linarith
  have h := tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N _ => m.u0Het w N) (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    (m.stackThetaSqW_nonneg w) hKq hpMeas hp hqMeas
    (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG) (fun _ => 0) 0
    (fun z hz => by simpa using m.tendstoInProb_cformC_u0Het_eHet w c hc hw hreg hG hz) hx
    (by simp)
  exact h

include hc hw hreg hG hedge hx in
omit hb in
/-- **Field `ue2` at real `x`**: `ũ₀ᵀ G₀'(x)² (Σ^{1/2} e) → 0`. -/
theorem tendstoInProb_cform2_u0Het_eHet :
    TendstoInProb μ (fun N ω => R4.cform2 (m.W0het w N ω) x (m.u0Het w N)
      (m.SigmaHalf w N *ᵥ m.eHet N ω)) 0 := by
  have hpMeas : ∀ N r, NullMeasurableSet {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ r} (μ N) :=
    fun N r => (measurableSet_const_prop N _).nullMeasurableSet
  have hp : Tendsto (fun N => μ N {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w})
      atTop (𝓝 1) := by
    have : ∀ N, {ω : Ω N | m.u0Het w N ⬝ᵥ m.u0Het w N ≤ m.stackThetaSqW w} = Set.univ :=
      fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_u0Het_self w N).le
    simp only [this, measure_univ]
    exact tendsto_const_nhds
  have hqMeas : ∀ N r, NullMeasurableSet {ω : Ω N |
      (m.SigmaHalf w N *ᵥ m.eHet N ω) ⬝ᵥ (m.SigmaHalf w N *ᵥ m.eHet N ω) ≤ r} (μ N) :=
    fun N r => (measurableSet_le (m.measurable_SigmaHalf_eHet_dot w N)
      measurable_const).nullMeasurableSet
  have hKq : 0 ≤ ∑ i, w i ^ 2 * c i + 1 := by
    have : 0 ≤ ∑ i, w i ^ 2 * c i := Finset.sum_nonneg fun i _ => by
      have := (hc i).le; positivity
    linarith
  have h := tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0het w N ω) (fun N ω => m.isHermitian_W0het w N ω)
    (fun N _ => m.u0Het w N) (fun N ω => m.SigmaHalf w N *ᵥ m.eHet N ω)
    (fun ε _ N => (m.measurableSet_lamMax_W0het_le w N (b + ε)).nullMeasurableSet) hedge
    (m.stackThetaSqW_nonneg w) hKq hpMeas hp hqMeas
    (m.tendsto_measure_SigmaHalf_eHet_norm_le w c hw hreg hG) (fun _ => 0) 0
    (fun z hz => by simpa using m.tendstoInProb_cform2C_u0Het_eHet w c hc hw hreg hG hz) hx
    (by simp)
  exact h

end RealForms

end MultiTableModel
end StackedSVD
