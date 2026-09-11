/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Stability
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.General.Deloc
import StackedSVD.RMT.General.MPtilde
import StackedSVD.RMT.General.QuadForm

/-!
# Item G3: the Marchenko-Pastur trace law at four moments

`notes/archive/prop_single_table_general.md` section 5, unit G3. Twin of `RMT/R1.lean` for a general
noise law `ν` (mean 0, variance 1, finite fourth moment): the two limits

* `d⁻¹ tr (d⁻¹ YᵀY - z)⁻¹ → MP.mC c z` in probability,
* `d⁻¹ tr ((d⁻¹ YᵀY - z)⁻¹)² → MP.mCDeriv c z` in probability,

both in the norm form of choice 28. The Gaussian file proves the first by Gaussian
concentration and a Stein equation; neither is available here, so the route is the standard
leave-one-out argument, which needs only four moments.

This file imports none of `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` or
`Vendor/COLT83/` (choice 8 of the note).

## Route

Write `A = gram Y`, `G = resolvC A z`, `B = gramC Y`, `Gc = resolvC B z`, `g k` for the `k`-th
row of `Y` divided by `√d`, and `A_k = A - vecMulVec (g k) (g k)`, the Gram matrix of `Y` with
row `k` set to zero.

1. `qformC_sub_vecMulVec_self` and `one_add_alpha_mul`: Sherman-Morrison on the quadratic form
   of the removed row, in the product form `(1 + α k) * (1 - s k) = 1`.
2. `gram_updateRow_zero`: zeroing a row of `Y` subtracts a rank-one matrix from `gram Y`.
3. `mul_one_add_alpha_gramC_diag`: the leave-one-out formula `z (1 + α k) (Gc) k k = -1`, from
   the companion identity of unit G2 read at the `(k, k)` entry.
4. `norm_trace_resolvC_sub_le`: the deterministic trace stability `‖tr G_k - tr G‖ ≤ 1/η`.
5. `norm_quad_le`: the self-consistent equation at a finite `d`, with the residual carried by
   the row average of `‖α k - d⁻¹ tr G‖`.
6. `integral_sq_alpha_sub_le`: the four-moment step of unit G1, with the Fubini split on the
   rows of the product law.
7. `tendstoInProb_stieltjesC_general` and `tendstoInProb_stieltjes2C_general`.
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace GenRMT

variable {p d : ℕ} {z : ℂ}

/-! ### The scaled row and the leave-one-out Gram matrix -/

/-- `g k`, the `k`-th row of `Y` scaled by `d^{-1/2}`. With this scaling
`gram Y = ∑ k, vecMulVec (g k) (g k)`. -/
noncomputable def rowVec (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) : Fin d → ℝ :=
  (Real.sqrt d)⁻¹ • Y k

theorem rowVec_apply (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) (a : Fin d) :
    rowVec Y k a = (Real.sqrt d)⁻¹ * Y k a := rfl

/-- The scaled square root, twice, is `d⁻¹`, over `ℂ`. -/
private theorem sqrt_inv_sq (d : ℕ) :
    (((Real.sqrt d)⁻¹ : ℝ) : ℂ) * (((Real.sqrt d)⁻¹ : ℝ) : ℂ) = (d : ℂ)⁻¹ := by
  rw [← Complex.ofReal_mul]
  have step : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = ((d : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d)]
  rw [step, Complex.ofReal_inv]
  congr 1

/-- **Zeroing a row is a rank-one downdate of the Gram matrix.** The leave-one-out matrix
`A_k` keeps the type `Fin d`, so no reindexing is needed anywhere in this file. -/
theorem gram_updateRow_zero (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) :
    gram (Matrix.updateRow Y k 0)
      = gram Y - Matrix.vecMulVec (rowVec Y k) (rowVec Y k) := by
  have hdd : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = ((d : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d)]
  ext a b
  have hg : ∀ X : Matrix (Fin p) (Fin d) ℝ,
      gram X a b = ((d : ℝ))⁻¹ * ∑ j, X j a * X j b := by
    intro X
    simp only [gram, Matrix.smul_apply, smul_eq_mul, Matrix.mul_apply, Matrix.transpose_apply]
  have hsum : ∑ j, (Matrix.updateRow Y k 0) j a * (Matrix.updateRow Y k 0) j b
      = (∑ j, Y j a * Y j b) - Y k a * Y k b := by
    have e1 := Finset.sum_erase_add Finset.univ
      (fun j => (Matrix.updateRow Y k 0) j a * (Matrix.updateRow Y k 0) j b) (Finset.mem_univ k)
    have e2 := Finset.sum_erase_add Finset.univ
      (fun j => Y j a * Y j b) (Finset.mem_univ k)
    have e3 : ∑ j ∈ Finset.univ.erase k,
        (Matrix.updateRow Y k 0) j a * (Matrix.updateRow Y k 0) j b
        = ∑ j ∈ Finset.univ.erase k, Y j a * Y j b :=
      Finset.sum_congr rfl fun j hj => by
        rw [Matrix.updateRow_ne (Finset.ne_of_mem_erase hj)]
    have e4 : (Matrix.updateRow Y k 0) k a * (Matrix.updateRow Y k 0) k b = 0 := by
      simp
    rw [← e1, e3, e4, ← e2]
    ring
  rw [Matrix.sub_apply, hg, hg, hsum, Matrix.vecMulVec_apply, rowVec_apply, rowVec_apply]
  have hlast : ((Real.sqrt d)⁻¹ * Y k a) * ((Real.sqrt d)⁻¹ * Y k b)
      = ((d : ℝ))⁻¹ * (Y k a * Y k b) := by
    rw [← hdd]; ring
  rw [hlast]
  ring


/-! ### The two scalars of the leave-one-out step -/

/-- `s k = g_kᵀ G g_k`, the quadratic form of the full resolvent at the scaled row `k`. -/
noncomputable def sRow (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) (k : Fin p) : ℂ :=
  R4C.qformC (gram Y) z (rowVec Y k)

/-- `α k = g_kᵀ G_k g_k`, the same form at the leave-one-out resolvent. -/
noncomputable def alphaRow (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) (k : Fin p) : ℂ :=
  R4C.qformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k)

/-- `1 - s k ≠ 0`, from `R4C.qformC_ne_one` of unit G2. -/
theorem one_sub_sRow_ne_zero (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (k : Fin p) :
    1 - sRow Y z k ≠ 0 :=
  sub_ne_zero.mpr (Ne.symm (R4C.qformC_ne_one (gram_isHermitian Y) hz (rowVec Y k)))

/-- **Sherman-Morrison on the quadratic form**, in the product form that needs no division:
`(1 + α k) (1 - s k) = 1`. -/
theorem one_add_alphaRow_mul (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (k : Fin p) :
    (1 + alphaRow Y z k) * (1 - sRow Y z k) = 1 := by
  have hne := one_sub_sRow_ne_zero Y hz k
  have hsm := R4C.cformC_sub_vecMulVec (gram_isHermitian Y) hz (rowVec Y k)
    (rowVec Y k) (rowVec Y k)
  have halpha : alphaRow Y z k
      = sRow Y z k + sRow Y z k * sRow Y z k / (1 - sRow Y z k) := by
    rw [alphaRow, gram_updateRow_zero]
    exact hsm
  rw [halpha]
  field_simp
  ring


/-! ### The diagonal entry of the companion resolvent -/

/-- The `(k, k)` entry of `Y M Yᵀ` is the quadratic form of `M` at row `k` of `Y`. -/
private theorem entry_kk (Y : Matrix (Fin p) (Fin d) ℝ) (M : Matrix (Fin d) (Fin d) ℂ)
    (k : Fin p) :
    (R4C.cmat' Y * M * (R4C.cmat' Y)ᵀ) k k
      = R4C.cvec (Y k) ⬝ᵥ (M *ᵥ R4C.cvec (Y k)) := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, R4C.cmat'_apply, dotProduct,
    Matrix.mulVec, R4C.cvec, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- `s k` written on the unscaled row. -/
private theorem sRow_eq_smul (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) (k : Fin p) :
    sRow Y z k
      = (d : ℂ)⁻¹ * (R4C.cvec (Y k) ⬝ᵥ (R4C.resolvC (gram Y) z *ᵥ R4C.cvec (Y k))) := by
  have hcv : R4C.cvec (rowVec Y k) = (((Real.sqrt d)⁻¹ : ℝ) : ℂ) • R4C.cvec (Y k) := by
    funext a
    simp [R4C.cvec, rowVec_apply]
  change R4C.cvec (rowVec Y k) ⬝ᵥ (R4C.resolvC (gram Y) z *ᵥ R4C.cvec (rowVec Y k)) = _
  rw [hcv, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_smul, smul_eq_mul,
    sqrt_inv_sq]

/-- **The companion identity at one diagonal entry.** `s k = 1 + z (Gc) k k`. -/
theorem sRow_eq_one_add (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (k : Fin p) :
    sRow Y z k = 1 + z * R4C.resolvC (gramC Y) z k k := by
  have hcomp := smul_cmat'_mul_resolvC_mul_transpose Y hz hd
  have hkk : (((d : ℂ))⁻¹ • (R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ)) k k
      = ((1 : Matrix (Fin p) (Fin p) ℂ) + z • R4C.resolvC (gramC Y) z) k k := by
    rw [hcomp]
  rw [Matrix.smul_apply, smul_eq_mul, entry_kk, Matrix.add_apply, Matrix.one_apply_eq,
    Matrix.smul_apply, smul_eq_mul] at hkk
  rw [sRow_eq_smul]
  exact hkk

/-- **The leave-one-out formula**, in the product form: `z (1 + α k) (Gc) k k = -1`. -/
theorem mul_one_add_alphaRow (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (k : Fin p) :
    z * (1 + alphaRow Y z k) * R4C.resolvC (gramC Y) z k k = -1 := by
  have h1 := one_add_alphaRow_mul Y hz k
  have h2 := sRow_eq_one_add Y hz hd k
  have h3 : (1 + alphaRow Y z k) * (1 - (1 + z * R4C.resolvC (gramC Y) z k k)) = 1 := by
    rw [← h2]; exact h1
  linear_combination -h3


/-! ### Two deterministic resolvent facts -/

/-- Every diagonal entry of the resolvent has norm at most `1/η`. -/
theorem norm_resolvC_diag_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (k : Fin D) : ‖R4C.resolvC W z k k‖ ≤ 1 / z.im := by
  have hcv : R4C.cvec (Pi.single k (1 : ℝ)) = Pi.single k (1 : ℂ) := by
    funext a
    by_cases h : a = k
    · subst h; simp [R4C.cvec]
    · simp [R4C.cvec, h]
  have hq : R4C.qformC W z (Pi.single k (1 : ℝ)) = R4C.resolvC W z k k := by
    change R4C.cvec (Pi.single k (1 : ℝ)) ⬝ᵥ
      (R4C.resolvC W z *ᵥ R4C.cvec (Pi.single k (1 : ℝ))) = _
    rw [hcv, Matrix.mulVec_single, single_dotProduct]
    simp
  have hy : (Pi.single k (1 : ℝ)) ⬝ᵥ (Pi.single k (1 : ℝ)) = 1 := by
    rw [single_dotProduct]
    simp
  have h := R4C.norm_qformC_le hW hz (Pi.single k (1 : ℝ))
  rw [hq, hy] at h
  exact h

/-- The resolvent moves across a bilinear form: it is symmetric. -/
private theorem resolvC_dotProduct {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : z.im ≠ 0) (x y : Fin D → ℂ) :
    x ⬝ᵥ (R4C.resolvC W z *ᵥ y) = (R4C.resolvC W z *ᵥ x) ⬝ᵥ y := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, ResolvDeriv.transpose_resolvC hW hz]

/-- `gᵀ G² g = ‖G g‖²`, without a conjugate. -/
private theorem qform2C_eq_dotProduct {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ}
    (hW : W.IsHermitian) (hz : z.im ≠ 0) (g : Fin D → ℝ) :
    R4C.qform2C W z g
      = (R4C.resolvC W z *ᵥ R4C.cvec g) ⬝ᵥ (R4C.resolvC W z *ᵥ R4C.cvec g) := by
  change R4C.cvec g ⬝ᵥ ((R4C.resolvC W z * R4C.resolvC W z) *ᵥ R4C.cvec g) = _
  rw [← Matrix.mulVec_mulVec, resolvC_dotProduct hW hz]

/-- `‖gᵀ G² g‖` is at most the eigenvalue sum `∑ w_a / |λ_a - z|²`. -/
private theorem norm_qform2C_le_sum {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (g : Fin D → ℝ) :
    ‖R4C.qform2C W z g‖
      ≤ ∑ a, ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2 / Complex.normSq ((hW.eigenvalues a : ℂ) - z) := by
  have hz' : z.im ≠ 0 := hz.ne'
  change ‖R4C.cform2C W z g g‖ ≤ _
  rw [R4C.cform2C_eq_sum hW hz' g g]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun a _ => ?_)
  rw [norm_mul, norm_pow]
  have h1 : ‖((hW.eigenvalues a : ℂ) - z)⁻¹‖ ^ 2
      = (Complex.normSq ((hW.eigenvalues a : ℂ) - z))⁻¹ := by
    simp only [norm_inv, inv_pow, Complex.sq_norm]
  have h2 : ‖((((R4.eigU hW)ᵀ *ᵥ g) a * ((R4.eigU hW)ᵀ *ᵥ g) a : ℝ) : ℂ)‖
      = ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2 := by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (mul_self_nonneg _)]
    ring
  rw [h1, h2, div_eq_inv_mul]


/-! F37 (2026-09-09): `im_qformC_eq` used to be repeated here (`private`); `Deloc.im_qformC_eq`
(`Deloc.lean:499`, now public) is the canonical copy, used directly. -/

/-- **The denominator bound of step 5.** `η ‖gᵀ G² g‖ ≤ ‖1 - gᵀ G g‖`, for every `g` and every
Hermitian `W`. Both sides are eigenvalue sums with the same weights; the left one drops the
factor `η` that the imaginary part of the right one carries. -/
theorem im_mul_norm_qform2C_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (g : Fin D → ℝ) :
    z.im * ‖R4C.qform2C W z g‖ ≤ ‖1 - R4C.qformC W z g‖ := by
  set T : ℝ := ∑ a, ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2
      / Complex.normSq ((hW.eigenvalues a : ℂ) - z) with hT
  have hTnn : 0 ≤ T := by
    refine Finset.sum_nonneg fun a _ => ?_
    have h := Complex.normSq_nonneg ((hW.eigenvalues a : ℂ) - z)
    positivity
  have h1 : ‖R4C.qform2C W z g‖ ≤ T := norm_qform2C_le_sum hW hz g
  have h2 : (R4C.qformC W z g).im = z.im * T := Deloc.im_qformC_eq hW hz g
  have h3 : z.im * T ≤ ‖1 - R4C.qformC W z g‖ := by
    have h4 : |(1 - R4C.qformC W z g).im| ≤ ‖1 - R4C.qformC W z g‖ :=
      Complex.abs_im_le_norm _
    have h5 : (1 - R4C.qformC W z g).im = -(z.im * T) := by
      simp [h2]
    rw [h5, abs_neg, abs_of_nonneg (mul_nonneg hz.le hTnn)] at h4
    exact h4
  calc z.im * ‖R4C.qform2C W z g‖ ≤ z.im * T := mul_le_mul_of_nonneg_left h1 hz.le
    _ ≤ ‖1 - R4C.qformC W z g‖ := h3


/-! ### Step 5: the trace of the leave-one-out resolvent -/

/-- **The trace of the downdate.** `tr G_k - tr G = (1 - s k)⁻¹ g_kᵀ G² g_k`. -/
theorem trace_resolvC_sub (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (k : Fin p) :
    (R4C.resolvC (gram (Matrix.updateRow Y k 0)) z).trace
        - (R4C.resolvC (gram Y) z).trace
      = (1 - sRow Y z k)⁻¹ * R4C.qform2C (gram Y) z (rowVec Y k) := by
  have hz' : z.im ≠ 0 := hz.ne'
  have hW := gram_isHermitian Y
  have htr : (Matrix.vecMulVec (R4C.resolvC (gram Y) z *ᵥ R4C.cvec (rowVec Y k))
      (R4C.resolvC (gram Y) z *ᵥ R4C.cvec (rowVec Y k))).trace
      = R4C.qform2C (gram Y) z (rowVec Y k) := by
    rw [qform2C_eq_dotProduct hW hz' (rowVec Y k)]
    simp [Matrix.trace, Matrix.diag, Matrix.vecMulVec_apply, dotProduct]
  rw [gram_updateRow_zero, R4C.resolvC_sub_vecMulVec hW hz (rowVec Y k), Matrix.trace_add,
    Matrix.trace_smul, smul_eq_mul, htr]
  simp only [sRow]
  ring

/-- **Step 5.** `‖tr G_k - tr G‖ ≤ 1/η`, deterministically, for every row and every
realization. Step 7 uses this twice: once to replace `α k` by `d⁻¹ tr G`, once inside the
error bound. -/
theorem norm_trace_resolvC_sub_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (k : Fin p) :
    ‖(R4C.resolvC (gram (Matrix.updateRow Y k 0)) z).trace
        - (R4C.resolvC (gram Y) z).trace‖ ≤ 1 / z.im := by
  have hne : (0 : ℝ) < ‖1 - sRow Y z k‖ :=
    norm_pos_iff.mpr (one_sub_sRow_ne_zero Y hz k)
  have hbnd : z.im * ‖R4C.qform2C (gram Y) z (rowVec Y k)‖ ≤ ‖1 - sRow Y z k‖ :=
    im_mul_norm_qform2C_le (gram_isHermitian Y) hz (rowVec Y k)
  rw [trace_resolvC_sub Y hz k, norm_mul, norm_inv, inv_mul_eq_div, div_le_div_iff₀ hne hz]
  rw [mul_comm]
  linarith

/-- The same bound on the normalized traces: `‖s_k(A_k) - s_k(A)‖ ≤ 1/(d η)`. -/
theorem norm_stieltjesC_sub_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (k : Fin p) :
    ‖R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z - R4C.stieltjesC (gram Y) z‖
      ≤ 1 / ((d : ℝ) * z.im) := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hstep : R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z - R4C.stieltjesC (gram Y) z
      = (d : ℂ)⁻¹ * ((R4C.resolvC (gram (Matrix.updateRow Y k 0)) z).trace
          - (R4C.resolvC (gram Y) z).trace) := by
    simp only [R4C.stieltjesC]
    ring
  rw [hstep, norm_mul, norm_inv, Complex.norm_natCast]
  have h := norm_trace_resolvC_sub_le Y hz k
  calc ((d : ℝ))⁻¹ * ‖(R4C.resolvC (gram (Matrix.updateRow Y k 0)) z).trace
        - (R4C.resolvC (gram Y) z).trace‖
      ≤ ((d : ℝ))⁻¹ * (1 / z.im) := by
        exact mul_le_mul_of_nonneg_left h (by positivity)
    _ = 1 / ((d : ℝ) * z.im) := by
        field_simp


/-! ### Step 7: the self-consistent equation at a finite `d` -/

/-- The row average of `‖α k - d⁻¹ tr G‖`, the residual that drives the quadratic. -/
noncomputable def resid (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) : ℝ :=
  ((p : ℝ))⁻¹ * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖

theorem resid_nonneg (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) : 0 ≤ resid Y z := by
  refine mul_nonneg (by positivity) (Finset.sum_nonneg fun k _ => norm_nonneg _)

/-- **The self-consistent equation.** Averaging the leave-one-out formula over the `p` rows and
eliminating the companion transform with the finite identity of unit G5 gives the MP quadratic
at the empirical ratio `p/d`, with a residual carried by the row deviations `α k - d⁻¹ tr G`. -/
theorem quad_eq_of_gram (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hp : 0 < p)
    (hd : 0 < d) :
    MP.quad ((p : ℝ) / d) z (R4C.stieltjesC (gram Y) z)
      = -((((p : ℝ) / d : ℝ)) : ℂ) * z * (p : ℂ)⁻¹
        * ∑ k, (alphaRow Y z k - R4C.stieltjesC (gram Y) z)
            * R4C.resolvC (gramC Y) z k k := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hP : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hD : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  set s : ℂ := R4C.stieltjesC (gram Y) z with hs
  set Tc : ℂ := (R4C.resolvC (gramC Y) z).trace with hTcdef
  set S : ℂ := ∑ k, (alphaRow Y z k - s) * R4C.resolvC (gramC Y) z k k with hSdef
  have htr : ∑ k, R4C.resolvC (gramC Y) z k k = Tc := rfl
  -- the row average of the leave-one-out formula
  have hrow : z * (1 + s) * Tc + z * S = -(p : ℂ) := by
    have hexp : ∀ k : Fin p, z * (1 + s) * R4C.resolvC (gramC Y) z k k
        + z * ((alphaRow Y z k - s) * R4C.resolvC (gramC Y) z k k)
        = z * (1 + alphaRow Y z k) * R4C.resolvC (gramC Y) z k k := fun k => by ring
    have hsum : ∑ k : Fin p, (z * (1 + s) * R4C.resolvC (gramC Y) z k k
        + z * ((alphaRow Y z k - s) * R4C.resolvC (gramC Y) z k k))
        = ∑ _k : Fin p, (-1 : ℂ) := by
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [hexp k]
      exact mul_one_add_alphaRow Y hz hd k
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, htr] at hsum
    simpa [hSdef] using hsum
  -- the finite identity of unit G5
  have hfin := stieltjesC_gram_sub_gramC Y hz hd
  have hTc : Tc = (d : ℂ) * s - ((p : ℂ) - d) / z := by
    have hst : R4C.stieltjesC (gramC Y) z = (p : ℂ)⁻¹ * Tc := rfl
    have hpp : (p : ℂ) * ((p : ℂ)⁻¹ * Tc) = Tc := by field_simp
    rw [hst, ← hs, hpp] at hfin
    linear_combination -hfin
  have hzS : z * S = -(p : ℂ) - z * (1 + s) * ((d : ℂ) * s - ((p : ℂ) - d) / z) := by
    rw [← hTc]
    linear_combination hrow
  have hcast : ((((p : ℝ) / d : ℝ)) : ℂ) = (p : ℂ) / (d : ℂ) := by
    push_cast
    ring
  change z * s ^ 2 + (z + 1 - ((((p : ℝ) / d : ℝ)) : ℂ)) * s + 1 = _
  rw [hcast]
  have hgoal : -((p : ℂ) / (d : ℂ)) * z * (p : ℂ)⁻¹ * S
      = -((p : ℂ) / (d : ℂ)) * (p : ℂ)⁻¹ * (z * S) := by ring
  rw [hgoal, hzS]
  field_simp
  ring


/-- **The residual form of step 7.** The MP quadratic at the empirical ratio is at most
`(p/d) (‖z‖/η)` times the row average of `‖α k - d⁻¹ tr G‖`. -/
theorem norm_quad_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) :
    ‖MP.quad ((p : ℝ) / d) z (R4C.stieltjesC (gram Y) z)‖
      ≤ ((p : ℝ) / d) * (‖z‖ / z.im) * resid Y z := by
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hcnn : (0 : ℝ) ≤ (p : ℝ) / d := by positivity
  rw [quad_eq_of_gram Y hz hp hd]
  set S : ℂ := ∑ k, (alphaRow Y z k - R4C.stieltjesC (gram Y) z)
      * R4C.resolvC (gramC Y) z k k with hSdef
  have hS : ‖S‖ ≤ (1 / z.im) * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ := by
    rw [hSdef]
    refine (norm_sum_le _ _).trans ?_
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun k _ => ?_
    rw [norm_mul, mul_comm (1 / z.im)]
    exact mul_le_mul_of_nonneg_left
      (norm_resolvC_diag_le (gramC_isHermitian Y) hz k) (norm_nonneg _)
  have hnorm : ‖-((((p : ℝ) / d : ℝ)) : ℂ) * z * (p : ℂ)⁻¹ * S‖
      = ((p : ℝ) / d) * ‖z‖ * ((p : ℝ))⁻¹ * ‖S‖ := by
    rw [norm_mul, norm_mul, norm_mul, norm_neg, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg hcnn, norm_inv, Complex.norm_natCast]
  rw [hnorm]
  have hstep : ((p : ℝ) / d) * ‖z‖ * ((p : ℝ))⁻¹ * ‖S‖
      ≤ ((p : ℝ) / d) * ‖z‖ * ((p : ℝ))⁻¹
        * ((1 / z.im) * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖) :=
    mul_le_mul_of_nonneg_left hS (by positivity)
  refine hstep.trans (le_of_eq ?_)
  rw [resid]
  ring


/-! ### The Frobenius norm of the resolvent -/

/-- **`∑_{a,b} ‖G a b‖² ≤ D / η²`.** The sharp bound that step 4 needs; the entrywise bound
`‖G a b‖ ≤ 1/η` would lose a factor `D` and make the error term of order one. The route is the
spectral form: `G Ḡᵀ = U diag(|λ_a - z|⁻²) Uᵀ`, whose trace is the sum of the squared entries
because `U` is real orthogonal. -/
theorem sum_normSq_resolvC_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) :
    ∑ a, ∑ b, ‖R4C.resolvC W z a b‖ ^ 2 ≤ (D : ℝ) / z.im ^ 2 := by
  have hz' : z.im ≠ 0 := hz.ne'
  set V : Matrix (Fin D) (Fin D) ℂ := R4C.cmat (R4.eigU hW) with hV
  set f : Fin D → ℂ := fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹ with hf
  have hVc : V.map (starRingEnd ℂ) = V := by
    ext a b
    simp [hV, R4C.cmat]
  have hVtc : (Vᵀ).map (starRingEnd ℂ) = Vᵀ := by
    ext a b
    simp [hV, R4C.cmat]
  have hM : R4C.resolvC W z = V * Matrix.diagonal f * Vᵀ := R4C.resolvC_eq_conj hW hz'
  have hMc : ((R4C.resolvC W z).map (starRingEnd ℂ))ᵀ
      = V * Matrix.diagonal (fun a => (starRingEnd ℂ) (f a)) * Vᵀ := by
    rw [hM, Matrix.map_mul, Matrix.map_mul, Matrix.diagonal_map (by simp), hVc, hVtc,
      Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_transpose,
      Matrix.diagonal_transpose, Matrix.mul_assoc]
  have hone : Vᵀ * V = 1 := R4C.transpose_ceigU_mul hW
  have hprod : R4C.resolvC W z * ((R4C.resolvC W z).map (starRingEnd ℂ))ᵀ
      = V * Matrix.diagonal (fun a => f a * (starRingEnd ℂ) (f a)) * Vᵀ := by
    have hstep : (V * Matrix.diagonal f * Vᵀ)
        * (V * Matrix.diagonal (fun a => (starRingEnd ℂ) (f a)) * Vᵀ)
        = V * (Matrix.diagonal f * ((Vᵀ * V)
            * Matrix.diagonal (fun a => (starRingEnd ℂ) (f a)))) * Vᵀ := by
      simp only [Matrix.mul_assoc]
    rw [hMc, hM, hstep, hone, Matrix.one_mul, Matrix.diagonal_mul_diagonal]
  have hnormsq : ∀ w : ℂ, ((‖w‖ ^ 2 : ℝ) : ℂ) = w * (starRingEnd ℂ) w := by
    intro w
    rw [Complex.sq_norm, Complex.mul_conj]
  have hfrob : ((∑ a, ∑ b, ‖R4C.resolvC W z a b‖ ^ 2 : ℝ) : ℂ)
      = (R4C.resolvC W z * ((R4C.resolvC W z).map (starRingEnd ℂ))ᵀ).trace := by
    rw [Matrix.trace, Complex.ofReal_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Complex.ofReal_sum, Matrix.diag_apply, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [hnormsq, Matrix.transpose_apply, Matrix.map_apply]
  have hkey : ((∑ a, ∑ b, ‖R4C.resolvC W z a b‖ ^ 2 : ℝ) : ℂ)
      = ((∑ a, ‖f a‖ ^ 2 : ℝ) : ℂ) := by
    rw [hfrob, hprod, R4C.trace_conj _ _ hone, Matrix.trace_diagonal,
      Complex.ofReal_sum]
    exact Finset.sum_congr rfl fun a _ => (hnormsq (f a)).symm
  have hreal : ∑ a, ∑ b, ‖R4C.resolvC W z a b‖ ^ 2 = ∑ a, ‖f a‖ ^ 2 :=
    Complex.ofReal_inj.mp hkey
  rw [hreal]
  calc ∑ a, ‖f a‖ ^ 2 ≤ ∑ _a : Fin D, ((z.im)⁻¹) ^ 2 := by
        refine Finset.sum_le_sum fun a _ => ?_
        exact pow_le_pow_left₀ (norm_nonneg _) (R4C.norm_inv_eigenvalue_sub_le hW hz a) 2
    _ = (D : ℝ) / z.im ^ 2 := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, inv_pow]
        ring


/-! ### Measurability

The entries of the resolvent are rational functions of the entries of the matrix, so they are
measurable. These six lemmas copy the shape of `R2.measurable_det`, `R2.measurable_adjugate`,
`R2.measurable_inv_entry` and `R2.measurable_resolvC_entry` (`RMT/R2.lean:210` to `:247`),
which this file may not import (choice 8 of the note); they are stated for an arbitrary
measurable family of real matrices, so that both `gram Y` and `gram (updateRow Y k 0)` are
instances. -/

section Measurability

variable {α : Type*} [MeasurableSpace α]

private theorem measurable_det {D : ℕ} {M : α → Matrix (Fin D) (Fin D) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) : Measurable fun a => (M a).det := by
  simp only [Matrix.det_apply']
  refine Finset.measurable_sum _ fun σ _ => ?_
  exact measurable_const.mul (Finset.measurable_prod _ fun i _ => hM (σ i) i)

private theorem measurable_adjugate {D : ℕ} {M : α → Matrix (Fin D) (Fin D) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin D) :
    Measurable fun a => (M a).adjugate i j := by
  simp only [Matrix.adjugate_apply]
  refine measurable_det fun k l => ?_
  by_cases h : k = j
  · subst h
    simp only [Matrix.updateRow_self]
    exact measurable_const
  · simp only [Matrix.updateRow_ne h]
    exact hM k l

private theorem measurable_inv_entry {D : ℕ} {M : α → Matrix (Fin D) (Fin D) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin D) :
    Measurable fun a => (M a)⁻¹ i j := by
  simp only [Matrix.inv_def, Matrix.smul_apply, smul_eq_mul, Ring.inverse_eq_inv']
  exact ((measurable_det hM).inv).mul (measurable_adjugate hM i j)

/-- Every entry of the resolvent of a measurable family of real matrices is measurable. -/
theorem measurable_resolvC_entry {D : ℕ} {W : α → Matrix (Fin D) (Fin D) ℝ}
    (hW : ∀ i j, Measurable fun a => W a i j) (z : ℂ) (i j : Fin D) :
    Measurable fun a => R4C.resolvC (W a) z i j := by
  refine measurable_inv_entry (fun k l => ?_) i j
  simp only [Matrix.sub_apply, R4C.cmat, Matrix.map_apply, Matrix.smul_apply, smul_eq_mul]
  exact (Complex.measurable_ofReal.comp (hW k l)).sub measurable_const

/-- The normalized trace of the resolvent is measurable. -/
theorem measurable_stieltjesC {D : ℕ} {W : α → Matrix (Fin D) (Fin D) ℝ}
    (hW : ∀ i j, Measurable fun a => W a i j) (z : ℂ) :
    Measurable fun a => R4C.stieltjesC (W a) z := by
  simp only [R4C.stieltjesC, Matrix.trace, Matrix.diag]
  exact measurable_const.mul
    (Finset.measurable_sum _ fun i _ => measurable_resolvC_entry hW z i i)

/-- The normalized trace of the squared resolvent is measurable. -/
theorem measurable_stieltjes2C {D : ℕ} {W : α → Matrix (Fin D) (Fin D) ℝ}
    (hW : ∀ i j, Measurable fun a => W a i j) (z : ℂ) :
    Measurable fun a => R4C.stieltjes2C (W a) z := by
  simp only [R4C.stieltjes2C, Matrix.trace, Matrix.diag, Matrix.mul_apply]
  refine measurable_const.mul (Finset.measurable_sum _ fun i _ => ?_)
  exact Finset.measurable_sum _ fun l _ =>
    (measurable_resolvC_entry hW z i l).mul (measurable_resolvC_entry hW z l i)

/-- The quadratic form of the resolvent, at a measurable family of vectors, is measurable. -/
theorem measurable_qformC {D : ℕ} {W : α → Matrix (Fin D) (Fin D) ℝ}
    (hW : ∀ i j, Measurable fun a => W a i j) (z : ℂ) {y : α → (Fin D → ℝ)}
    (hy : ∀ i, Measurable fun a => y a i) :
    Measurable fun a => R4C.qformC (W a) z (y a) := by
  have h : ∀ a, R4C.qformC (W a) z (y a)
      = ∑ i, ∑ j, R4C.resolvC (W a) z i j * ((y a i : ℝ) : ℂ) * ((y a j : ℝ) : ℂ) := by
    intro a
    exact bil_eq_sum (R4C.resolvC (W a) z) (y a) (y a)
  simp only [h]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact ((measurable_resolvC_entry hW z i j).mul
    (Complex.measurable_ofReal.comp (hy i))).mul (Complex.measurable_ofReal.comp (hy j))

end Measurability


section MeasurableGram

variable {α : Type*} [MeasurableSpace α] {P D : ℕ}

/-- The entries of `gram` are measurable in the entries of the matrix. -/
theorem measurable_gram_entry {X : α → Matrix (Fin P) (Fin D) ℝ}
    (hX : ∀ k l, Measurable fun a => X a k l) (i j : Fin D) :
    Measurable fun a => gram (X a) i j := by
  have h : ∀ a, gram (X a) i j = ((D : ℝ))⁻¹ * ∑ l, X a l i * X a l j := by
    intro a
    simp only [gram, Matrix.smul_apply, smul_eq_mul, Matrix.mul_apply, Matrix.transpose_apply]
  simp only [h]
  exact measurable_const.mul (Finset.measurable_sum _ fun l _ => (hX l i).mul (hX l j))

theorem measurable_matrix_entry (k : Fin P) (l : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => Y k l :=
  (measurable_pi_apply l).comp (measurable_pi_apply k)

theorem measurable_updateRow_entry (k l : Fin P) (i : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => (Matrix.updateRow Y k 0) l i := by
  by_cases h : l = k
  · subst h
    have hz0 : ∀ Y : Matrix (Fin P) (Fin D) ℝ, (Matrix.updateRow Y l 0) l i = 0 := by
      intro Y; simp
    simp only [hz0]
    exact measurable_const
  · simp only [Matrix.updateRow_ne h]
    exact measurable_matrix_entry l i

theorem measurable_rowVec (k : Fin P) (i : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => rowVec Y k i := by
  simp only [rowVec_apply]
  exact measurable_const.mul (measurable_matrix_entry k i)

/-- The entries of `gram Y` and of the leave-one-out `gram (updateRow Y k 0)`. -/
theorem measurable_gram_self (i j : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => gram Y i j :=
  measurable_gram_entry (fun l i' => measurable_matrix_entry l i') i j

theorem measurable_gram_loo (k : Fin P) (i j : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => gram (Matrix.updateRow Y k 0) i j :=
  measurable_gram_entry (fun l i' => measurable_updateRow_entry k l i') i j

theorem measurable_alphaRow (z : ℂ) (k : Fin P) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => alphaRow Y z k :=
  measurable_qformC (measurable_gram_loo k) z (measurable_rowVec k)

theorem measurable_resid (z : ℂ) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => resid Y z := by
  simp only [resid]
  refine measurable_const.mul (Finset.measurable_sum _ fun k _ => ?_)
  exact ((measurable_alphaRow z k).sub (measurable_stieltjesC measurable_gram_self z)).norm

end MeasurableGram


/-! ### Step 4: the quadratic form of a row against its trace -/

/-- **The row split of the product law.** For a nonnegative measurable `F` on `(n+1) × D`
matrices, the integral against `noiseMatrix ν (n+1) D` is at most any bound `C` that holds for
the integral in row `k` alone, at every value of the other `n` rows. The statement is in
`ℝ≥0∞`, so Tonelli applies and no integrability hypothesis is needed. Public since F37
(2026-09-09), the canonical copy for `RMT/General/` (the private twins of `Iso.lean`, named
`lintegral_split_le`, and `IsoMixed.lean`, named `lintegral_split_leMx`, are dropped). -/
theorem lintegral_noiseMatrix_le {n D : ℕ} {ν : Measure ℝ} [IsProbabilityMeasure ν]
    {F : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ≥0∞} (hF : Measurable F) (k : Fin (n + 1))
    {C : ℝ≥0∞}
    (hC : ∀ r : Fin n → Fin D → ℝ,
      ∫⁻ t, F (Fin.insertNth k t r) ∂(Measure.pi fun _ : Fin D => ν) ≤ C) :
    ∫⁻ Y, F Y ∂(noiseMatrix ν (n + 1) D) ≤ C := by
  set ρ : Measure (Fin D → ℝ) := Measure.pi fun _ : Fin D => ν with hρ
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => (Fin D → ℝ)) k with he
  have hmp : MeasurePreserving e (noiseMatrix ν (n + 1) D)
      (ρ.prod (Measure.pi fun _ : Fin n => ρ)) :=
    measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) => ρ) k
  have hsym := hmp.symm e
  have h1 : ∫⁻ q, F (e.symm q) ∂(ρ.prod (Measure.pi fun _ : Fin n => ρ))
      = ∫⁻ Y, F Y ∂(noiseMatrix ν (n + 1) D) := hsym.lintegral_comp hF
  have h2 : ∫⁻ q, F (e.symm q) ∂(ρ.prod (Measure.pi fun _ : Fin n => ρ))
      = ∫⁻ r, ∫⁻ t, F (e.symm (t, r)) ∂ρ ∂(Measure.pi fun _ : Fin n => ρ) :=
    lintegral_prod_symm _ (hF.comp e.symm.measurable).aemeasurable
  rw [← h1, h2]
  calc ∫⁻ r, ∫⁻ t, F (e.symm (t, r)) ∂ρ ∂(Measure.pi fun _ : Fin n => ρ)
      ≤ ∫⁻ _r : Fin n → Fin D → ℝ, C ∂(Measure.pi fun _ : Fin n => ρ) :=
        lintegral_mono fun r => hC r
    _ = C := by simp


/-- The quadratic form at the scaled row is the `bil` of the rescaled resolvent. -/
private theorem qformC_smul_eq_bil {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ)
    (t : Fin D → ℝ) :
    R4C.qformC W z ((Real.sqrt D)⁻¹ • t) = bil ((D : ℂ)⁻¹ • R4C.resolvC W z) t t := by
  have hl : R4C.qformC W z ((Real.sqrt D)⁻¹ • t)
      = bil (R4C.resolvC W z) ((Real.sqrt D)⁻¹ • t) ((Real.sqrt D)⁻¹ • t) := rfl
  rw [hl, bil_eq_sum, bil_eq_sum]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
  have ha : (((Real.sqrt D)⁻¹ • t) a : ℝ) = (Real.sqrt D)⁻¹ * t a := rfl
  have hb : (((Real.sqrt D)⁻¹ • t) b : ℝ) = (Real.sqrt D)⁻¹ * t b := rfl
  rw [ha, hb, Matrix.smul_apply, smul_eq_mul, Complex.ofReal_mul, Complex.ofReal_mul]
  have hs := sqrt_inv_sq D
  ring_nf
  ring_nf at hs
  linear_combination (R4C.resolvC W z a b * (t a : ℂ) * (t b : ℂ)) * hs

/-- **Step 4.** The quadratic form of row `k` in the leave-one-out resolvent is within
`(ν₄ + 2)/(D η²)` of the leave-one-out trace, in mean square. Unit G1
(`GenRMT.integral_normSq_bil_sub_trace_le`) supplies the four-moment bound at a fixed matrix,
and the Fubini split of the rows supplies the independence of row `k` from that matrix. -/
theorem lintegral_normSq_alphaRow_sub_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) (k : Fin (n + 1)) :
    ∫⁻ Y, ENNReal.ofReal (‖alphaRow Y z k
        - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2)
        ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2)) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hm : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      ENNReal.ofReal (‖alphaRow Y z k
        - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2) :=
    ENNReal.measurable_ofReal.comp
      (((measurable_alphaRow z k).sub
        (measurable_stieltjesC (measurable_gram_loo k) z)).norm.pow_const 2)
  refine lintegral_noiseMatrix_le hm k fun r => ?_
  set X : Fin (n + 1) → Fin D → ℝ := Fin.insertNth k (0 : Fin D → ℝ) r with hX
  set A : Matrix (Fin D) (Fin D) ℂ := (D : ℂ)⁻¹ • R4C.resolvC (gram X) z with hA
  have hupd : ∀ t : Fin D → ℝ, Matrix.updateRow (Fin.insertNth k t r) k 0 = X := by
    intro t
    change Function.update (Fin.insertNth k t r) k (0 : Fin D → ℝ) = X
    rw [hX]
    simp
  have hrow : ∀ t : Fin D → ℝ,
      rowVec (Fin.insertNth k t r) k = (Real.sqrt D)⁻¹ • t := by
    intro t
    have h2 : (Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k t r) k = t := by
      simp
    change (Real.sqrt D)⁻¹
        • ((Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k t r) k)
      = (Real.sqrt D)⁻¹ • t
    rw [h2]
  have hpt : ∀ t : Fin D → ℝ,
      alphaRow (Fin.insertNth k t r) z k
        - R4C.stieltjesC (gram (Matrix.updateRow (Fin.insertNth k t r) k 0)) z
      = bil A t t - A.trace := by
    intro t
    have hal : alphaRow (Fin.insertNth k t r) z k
        = R4C.qformC (gram (Matrix.updateRow (Fin.insertNth k t r) k 0)) z
            (rowVec (Fin.insertNth k t r) k) := rfl
    rw [hal, hupd t, hrow t, qformC_smul_eq_bil, ← hA]
    congr 1
    rw [hA, Matrix.trace_smul, smul_eq_mul]
    rfl
  simp only [hpt]
  obtain ⟨hint, hbound⟩ := integral_normSq_bil_sub_trace_le hν A
  have hfrob : ∑ a, ∑ b, ‖A a b‖ ^ 2 ≤ 1 / ((D : ℝ) * z.im ^ 2) := by
    have hentry : ∀ a b : Fin D, ‖A a b‖ ^ 2
        = ((D : ℝ))⁻¹ ^ 2 * ‖R4C.resolvC (gram X) z a b‖ ^ 2 := by
      intro a b
      rw [hA, Matrix.smul_apply, smul_eq_mul, norm_mul, norm_inv, Complex.norm_natCast, mul_pow]
    have hsum : ∑ a, ∑ b, ‖A a b‖ ^ 2
        = ((D : ℝ))⁻¹ ^ 2 * ∑ a, ∑ b, ‖R4C.resolvC (gram X) z a b‖ ^ 2 := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun a _ => ?_
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun b _ => hentry a b
    rw [hsum]
    have h := sum_normSq_resolvC_le (gram_isHermitian X) hz
    calc ((D : ℝ))⁻¹ ^ 2 * ∑ a, ∑ b, ‖R4C.resolvC (gram X) z a b‖ ^ 2
        ≤ ((D : ℝ))⁻¹ ^ 2 * ((D : ℝ) / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h (by positivity)
      _ = 1 / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
    have := hν.one_le_integral_pow_four
    linarith
  have hfin : ∫ t, ‖bil A t t - A.trace‖ ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2) := by
    refine hbound.trans ?_
    calc ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, ‖A a b‖ ^ 2
        ≤ ((∫ x, x ^ 4 ∂ν) + 2) * (1 / ((D : ℝ) * z.im ^ 2)) :=
          mul_le_mul_of_nonneg_left hfrob hnu
      _ = ((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2) := by ring
  rw [← ofReal_integral_eq_lintegral_ofReal hint
    (Filter.Eventually.of_forall fun t => sq_nonneg _)]
  exact ENNReal.ofReal_le_ofReal hfin


/-! ### From the row bound to the residual -/

/-- Cauchy-Schwarz on the row average: the square of `resid` is at most the average of the
squares. -/
theorem sq_resid_le (Y : Matrix (Fin p) (Fin d) ℝ) (z : ℂ) (hp : 0 < p) :
    resid Y z ^ 2
      ≤ ((p : ℝ))⁻¹ * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ ^ 2 := by
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp
  have hcs : (∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖) ^ 2
      ≤ (p : ℝ) * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ ^ 2 := by
    have h := sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset (Fin p)))
      (f := fun k => ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖)
    simpa using h
  rw [resid, mul_pow]
  have hstep : ((p : ℝ))⁻¹ ^ 2
      * (∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖) ^ 2
      ≤ ((p : ℝ))⁻¹ ^ 2
        * ((p : ℝ) * ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ ^ 2) :=
    mul_le_mul_of_nonneg_left hcs (by positivity)
  refine hstep.trans (le_of_eq ?_)
  field_simp

/-- Splitting the row deviation at the leave-one-out trace: the four-moment term plus the
deterministic `1/(d η)` of step 5. -/
theorem sq_norm_row_dev_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (k : Fin p) :
    ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ ^ 2
      ≤ 2 * ‖alphaRow Y z k - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2
        + 2 * (1 / ((d : ℝ) * z.im)) ^ 2 := by
  set u : ℂ := alphaRow Y z k - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z with hu
  set v : ℂ := R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z
      - R4C.stieltjesC (gram Y) z with hv
  have hsplit : alphaRow Y z k - R4C.stieltjesC (gram Y) z = u + v := by
    rw [hu, hv]; ring
  have hvb : ‖v‖ ≤ 1 / ((d : ℝ) * z.im) := norm_stieltjesC_sub_le Y hz hd k
  have hvnn : (0 : ℝ) ≤ 1 / ((d : ℝ) * z.im) := le_trans (norm_nonneg _) hvb
  rw [hsplit]
  have h1 : ‖u + v‖ ^ 2 ≤ (‖u‖ + ‖v‖) ^ 2 :=
    pow_le_pow_left₀ (norm_nonneg _) (norm_add_le u v) 2
  have h2 : (‖u‖ + ‖v‖) ^ 2 ≤ 2 * ‖u‖ ^ 2 + 2 * ‖v‖ ^ 2 := by nlinarith [sq_nonneg (‖u‖ - ‖v‖)]
  have h3 : ‖v‖ ^ 2 ≤ (1 / ((d : ℝ) * z.im)) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hvb 2
  linarith


/-- The pointwise bound that feeds the integral: the row average of the squared deviations
against the leave-one-out traces, plus the deterministic term. -/
private theorem sq_resid_le_sum {n D : ℕ} (Y : Matrix (Fin (n + 1)) (Fin D) ℝ) (hz : 0 < z.im)
    (hD : 0 < D) :
    resid Y z ^ 2
      ≤ (∑ k, 2 * ((n : ℝ) + 1)⁻¹ * ‖alphaRow Y z k
            - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2)
        + 2 * (1 / ((D : ℝ) * z.im)) ^ 2 := by
  have hpR : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  set a : Fin (n + 1) → ℝ := fun k => ‖alphaRow Y z k
      - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2 with ha
  set e : ℝ := 2 * (1 / ((D : ℝ) * z.im)) ^ 2 with he
  have hb : ∑ k, ‖alphaRow Y z k - R4C.stieltjesC (gram Y) z‖ ^ 2
      ≤ ∑ k : Fin (n + 1), (2 * a k + e) :=
    Finset.sum_le_sum fun k _ => sq_norm_row_dev_le Y hz hD k
  have hstep := (sq_resid_le Y z (Nat.succ_pos n)).trans
    (mul_le_mul_of_nonneg_left hb (by positivity))
  refine hstep.trans (le_of_eq ?_)
  have h1 : ∑ k : Fin (n + 1), (2 * a k + e) = 2 * (∑ k, a k) + ((n : ℝ) + 1) * e := by
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul]
    push_cast
    ring
  have h2 : ∑ k : Fin (n + 1), 2 * ((n : ℝ) + 1)⁻¹ * a k
      = 2 * ((n : ℝ) + 1)⁻¹ * ∑ k, a k := (Finset.mul_sum _ _ _).symm
  rw [h1, h2]
  push_cast
  field_simp


/-- **The residual has a vanishing second moment.** Step 4 at every row, the deterministic
step 5, and the row split of the product law. The bound is `O(1/D)`, so Chebyshev sends the
residual to zero in probability. -/
theorem lintegral_sq_resid_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    (hD : 0 < D) :
    ∫⁻ Y, ENNReal.ofReal (resid Y z ^ 2) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (2 * (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2))
          + 2 * (1 / ((D : ℝ) * z.im)) ^ 2) := by
  have hprob := hν.prob
  have hpR : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
    have := integral_pow_four_nonneg ν
    linarith
  set Qr : ℝ := ((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2) with hQr
  have hQnn : 0 ≤ Qr := div_nonneg hnu (by positivity)
  set Er : ℝ := 2 * (1 / ((D : ℝ) * z.im)) ^ 2 with hEr
  have hEnn : (0 : ℝ) ≤ Er := by rw [hEr]; positivity
  set cr : ℝ := 2 * ((n : ℝ) + 1)⁻¹ with hcr
  have hcnn : (0 : ℝ) ≤ cr := by rw [hcr]; positivity
  set Hk : Matrix (Fin (n + 1)) (Fin D) ℝ → Fin (n + 1) → ℝ≥0∞ := fun Y k =>
    ENNReal.ofReal (‖alphaRow Y z k
      - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2) with hHk
  have hHkm : ∀ k, Measurable fun Y => Hk Y k := fun k =>
    ENNReal.measurable_ofReal.comp
      (((measurable_alphaRow z k).sub
        (measurable_stieltjesC (measurable_gram_loo k) z)).norm.pow_const 2)
  have hpt : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ENNReal.ofReal (resid Y z ^ 2)
      ≤ (∑ k, ENNReal.ofReal cr * Hk Y k) + ENNReal.ofReal Er := by
    intro Y
    refine (ENNReal.ofReal_le_ofReal (sq_resid_le_sum Y hz hD)).trans (le_of_eq ?_)
    rw [ENNReal.ofReal_add (Finset.sum_nonneg fun k _ => by positivity) hEnn]
    congr 1
    rw [ENNReal.ofReal_sum_of_nonneg fun k _ => by positivity]
    exact Finset.sum_congr rfl fun k _ => ENNReal.ofReal_mul hcnn
  calc ∫⁻ Y, ENNReal.ofReal (resid Y z ^ 2) ∂(noiseMatrix ν (n + 1) D)
      ≤ ∫⁻ Y, ((∑ k, ENNReal.ofReal cr * Hk Y k) + ENNReal.ofReal Er)
          ∂(noiseMatrix ν (n + 1) D) := lintegral_mono hpt
    _ = (∑ k, ENNReal.ofReal cr * ∫⁻ Y, Hk Y k ∂(noiseMatrix ν (n + 1) D))
          + ENNReal.ofReal Er := by
        rw [lintegral_add_right _ measurable_const,
          lintegral_finsetSum _ fun k _ => (hHkm k).const_mul _, lintegral_const, measure_univ,
          mul_one]
        exact congrArg (· + ENNReal.ofReal Er)
          (Finset.sum_congr rfl fun k _ => lintegral_const_mul _ (hHkm k))
    _ ≤ (∑ _k : Fin (n + 1), ENNReal.ofReal cr * ENNReal.ofReal Qr) + ENNReal.ofReal Er := by
        gcongr with k
        exact lintegral_normSq_alphaRow_sub_le hν hz hD k
    _ = ENNReal.ofReal (2 * Qr + Er) := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, ← ENNReal.ofReal_mul hcnn,
          ← ENNReal.ofReal_nsmul, ← ENNReal.ofReal_add (by positivity) hEnn]
        congr 1
        rw [nsmul_eq_mul, hcr]
        push_cast
        field_simp


/-- The Chebyshev constant of the residual, `O(1/D)`. -/
noncomputable def residConst (ν : Measure ℝ) (D : ℕ) (z : ℂ) : ℝ :=
  2 * (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2)) + 2 * (1 / ((D : ℝ) * z.im)) ^ 2

theorem residConst_nonneg (ν : Measure ℝ) (D : ℕ) {z : ℂ} (hz : 0 < z.im) :
    0 ≤ residConst ν D z := by
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
    have := integral_pow_four_nonneg ν
    linarith
  rw [residConst]
  have h1 : (0 : ℝ) ≤ ((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2) :=
    div_nonneg hnu (by positivity)
  positivity

/-- `residConst → 0` along any index with `D → ∞`. -/
theorem tendsto_residConst {dN : ℕ → ℕ} (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im)
    (hd : Tendsto dN atTop atTop) :
    Tendsto (fun N => residConst ν (dN N) z) atTop (𝓝 0) := by
  have hinv : Tendsto (fun N => ((dN N : ℝ))⁻¹) atTop (𝓝 0) :=
    (tendsto_natCast_atTop_atTop.comp hd).inv_tendsto_atTop
  have h1 : Tendsto (fun N => 2 * (((∫ x, x ^ 4 ∂ν) + 2) / ((dN N : ℝ) * z.im ^ 2)))
      atTop (𝓝 0) := by
    have := hinv.const_mul (2 * (((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2))
    rw [mul_zero] at this
    refine this.congr fun N => ?_
    rcases eq_or_ne ((dN N : ℝ)) 0 with h0 | h0
    · rw [h0]; simp
    · field_simp
  have h2 : Tendsto (fun N => 2 * (1 / ((dN N : ℝ) * z.im)) ^ 2) atTop (𝓝 0) := by
    have hsq : Tendsto (fun N => (((dN N : ℝ))⁻¹) ^ 2) atTop (𝓝 0) := by
      have := hinv.pow 2
      simpa using this
    have := hsq.const_mul (2 * (1 / z.im) ^ 2)
    rw [mul_zero] at this
    refine this.congr fun N => ?_
    rcases eq_or_ne ((dN N : ℝ)) 0 with h0 | h0
    · rw [h0]; simp
    · field_simp
  have := h1.add h2
  rw [add_zero] at this
  exact this.congr fun N => rfl

/-- **Chebyshev on the residual.** The residual exceeds `ε` on a set of measure at most
`residConst / ε²`. -/
theorem measure_resid_ge_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    (hP : 0 < P) (hD : 0 < D) {ε : ℝ} (hε : 0 < ε) :
    noiseMatrix ν P D {Y | ε ≤ resid Y z}
      ≤ ENNReal.ofReal (residConst ν D z / ε ^ 2) := by
  obtain ⟨n, rfl⟩ : ∃ n, P = n + 1 := ⟨P - 1, (Nat.succ_pred_eq_of_pos hP).symm⟩
  have hprob := hν.prob
  have hsub : {Y : Matrix (Fin (n + 1)) (Fin D) ℝ | ε ≤ resid Y z}
      ⊆ {Y | ENNReal.ofReal (ε ^ 2) ≤ ENNReal.ofReal (resid Y z ^ 2)} := fun Y hY =>
    ENNReal.ofReal_le_ofReal (pow_le_pow_left₀ hε.le hY 2)
  refine (measure_mono hsub).trans ?_
  have hme : AEMeasurable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      ENNReal.ofReal (resid Y z ^ 2)) (noiseMatrix ν (n + 1) D) :=
    (ENNReal.measurable_ofReal.comp ((measurable_resid z).pow_const 2)).aemeasurable
  have hε2 : ENNReal.ofReal (ε ^ 2) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    positivity
  refine (meas_ge_le_lintegral_div hme hε2 ENNReal.ofReal_ne_top).trans ?_
  calc (∫⁻ Y, ENNReal.ofReal (resid Y z ^ 2) ∂(noiseMatrix ν (n + 1) D))
        / ENNReal.ofReal (ε ^ 2)
      ≤ ENNReal.ofReal (residConst ν D z) / ENNReal.ofReal (ε ^ 2) := by
        gcongr
        exact lintegral_sq_resid_le hν hz hD
    _ = ENNReal.ofReal (residConst ν D z / ε ^ 2) :=
        (ENNReal.ofReal_div_of_pos (by positivity)).symm


/-! ### Step 8: from the residual to the root -/

/-- Changing `c` in the quadratic costs `|c - c'| ‖w‖`. Unit G3a does not copy
`R1.quad_sub_quad` (`RMT/R1.lean:651`); this is that one line. -/
private theorem quad_sub_quad (c c' : ℝ) (z w : ℂ) :
    MP.quad c z w - MP.quad c' z w = ((c' : ℂ) - (c : ℂ)) * w := by
  change (z * w ^ 2 + (z + 1 - (c : ℂ)) * w + 1)
      - (z * w ^ 2 + (z + 1 - (c' : ℂ)) * w + 1) = _
  ring

/-- **Steps 7 and 8, deterministic.** A small residual and a small `|p/d - c|` put the
empirical transform near the root `MP.mC c z`. The `1/2` threshold is the one
`GenRMT.norm_ge_rootLb` needs. -/
theorem norm_stieltjesC_sub_mC_le {c : ℝ} (Y : Matrix (Fin p) (Fin d) ℝ) (hc : 0 < c)
    (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d)
    (hsmall : ((p : ℝ) / d) * (‖z‖ / z.im) * resid Y z + |(p : ℝ) / d - c| * (1 / z.im)
      ≤ 1 / 2) :
    ‖R4C.stieltjesC (gram Y) z - MP.mC c z‖
      ≤ (((p : ℝ) / d) * (‖z‖ / z.im) * resid Y z + |(p : ℝ) / d - c| * (1 / z.im))
        / (‖z‖ * (z.im * rootLb c z ^ 2)) := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hδ : 0 < rootLb c z := rootLb_pos hz0
  have hDen : 0 < ‖z‖ * (z.im * rootLb c z ^ 2) := by positivity
  set s : ℂ := R4C.stieltjesC (gram Y) z with hs
  set G : ℝ := ((p : ℝ) / d) * (‖z‖ / z.im) * resid Y z + |(p : ℝ) / d - c| * (1 / z.im)
    with hG
  have hsn : ‖s‖ ≤ 1 / z.im := R4C.norm_stieltjesC_le (gram_isHermitian Y) hz hd
  have hR : ‖MP.quad c z s‖ ≤ G := by
    have hq1 : ‖MP.quad ((p : ℝ) / d) z s‖
        ≤ ((p : ℝ) / d) * (‖z‖ / z.im) * resid Y z := norm_quad_le Y hz hp hd
    have hsub : MP.quad c z s - MP.quad ((p : ℝ) / d) z s
        = ((((p : ℝ) / d : ℝ)) : ℂ) * s - ((c : ℝ) : ℂ) * s := by
      rw [quad_sub_quad]
      ring
    have hcoe : ‖((((p : ℝ) / d : ℝ)) : ℂ) - ((c : ℝ) : ℂ)‖ = |(p : ℝ) / d - c| := by
      rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
    have he : MP.quad c z s = MP.quad ((p : ℝ) / d) z s
        + (((((p : ℝ) / d : ℝ)) : ℂ) - ((c : ℝ) : ℂ)) * s := by
      have := hsub
      linear_combination this
    have h2 : ‖MP.quad c z s‖
        ≤ ‖MP.quad ((p : ℝ) / d) z s‖ + |(p : ℝ) / d - c| * ‖s‖ := by
      rw [he]
      refine (norm_add_le _ _).trans ?_
      rw [norm_mul, hcoe]
    have h3 : |(p : ℝ) / d - c| * ‖s‖ ≤ |(p : ℝ) / d - c| * (1 / z.im) :=
      mul_le_mul_of_nonneg_left hsn (abs_nonneg _)
    rw [hG]
    linarith
  set rr : ℂ := -(MP.quad c z s) with hrr
  have hq : MP.quad c z s = -rr := by rw [hrr, neg_neg]
  have hrrn : ‖rr‖ = ‖MP.quad c z s‖ := by rw [hrr, norm_neg]
  have hlb : rootLb c z ≤ ‖s‖ := norm_ge_rootLb hz0 hq (by rw [hrrn]; linarith)
  have him : z.im * ‖s‖ ^ 2 ≤ s.im := im_stieltjesC_ge (gram_isHermitian Y) hz
  have hsq : rootLb c z ^ 2 ≤ ‖s‖ ^ 2 := pow_le_pow_left₀ hδ.le hlb 2
  have h5 : z.im * rootLb c z ^ 2 ≤ z.im * ‖s‖ ^ 2 := mul_le_mul_of_nonneg_left hsq hz.le
  have h6 : 0 < z.im * rootLb c z ^ 2 := by positivity
  have himpos : 0 < s.im := by linarith
  have hbnd := norm_sub_root_le hc hz (MP.im_mC_pos hc.le hz) (MP.quad_mC hc.le hz) himpos hq
  rw [hrrn] at hbnd
  have hstep : ‖s - MP.mC c z‖ ≤ ‖MP.quad c z s‖ / (‖z‖ * (z.im * rootLb c z ^ 2)) := by
    refine hbnd.trans ?_
    refine div_le_div_of_nonneg_left (norm_nonneg _) hDen ?_
    exact mul_le_mul_of_nonneg_left (by linarith) hzn.le
  refine hstep.trans ?_
  gcongr


/-! ### R1a and R1b at four moments -/

section Limits

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {pN dN : ℕ → ℕ}

/-- **R1a at four moments.** `d⁻¹ tr (d⁻¹YᵀY - z)⁻¹ → MP.mC c z` in probability, for i.i.d.
entries of law `ν`. Twin of `R1.tendstoInProb_stieltjesC` (`RMT/R1.lean:769`), same
conclusion, `hY` at `noiseMatrix ν` in place of `gaussianMatrix`.

Proof: the leave-one-out identity turns the MP quadratic at the empirical ratio into the row
average of `‖α k - d⁻¹ tr G‖` (`norm_quad_le`); that average has second moment `O(1/d)`
(`lintegral_sq_resid_le`), so Chebyshev makes it small with probability tending to one; and on
that event the stability of the quadratic (`norm_stieltjesC_sub_mC_le`) puts the transform
within `ε` of the root. No concentration inequality and no Stein identity is used, so four
moments suffice. -/
theorem tendstoInProb_stieltjesC_general {c : ℝ} {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hc : 0 < c) (hz : 0 < z.im) (hd : Tendsto dN atTop atTop)
    (hp : ∀ N, 0 < pN N) (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c))
    (Y : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (dN N)) ℝ)
    (hY : ∀ N, HasLaw (Y N) (noiseMatrix ν (pN N) (dN N)) (μ N)) :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖) 0 := by
  intro ε hε
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hδ : 0 < rootLb c z := rootLb_pos hz0
  have hDen : 0 < ‖z‖ * (z.im * rootLb c z ^ 2) := by positivity
  have hcp1 : (0 : ℝ) < c + 1 := by linarith
  set T : ℝ := min (1 / 2) (ε * (‖z‖ * (z.im * rootLb c z ^ 2)) / 2) with hT
  have hTpos : 0 < T := lt_min (by norm_num) (by positivity)
  set t : ℝ := T * z.im / (2 * (c + 1) * ‖z‖) with ht
  have htpos : 0 < t := by rw [ht]; positivity
  have hBto : Tendsto (fun N => ENNReal.ofReal (residConst ν (dN N) z / t ^ 2)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => residConst ν (dN N) z / t ^ 2) atTop (𝓝 0) := by
      simpa using (tendsto_residConst ν hz hd).div_const (t ^ 2)
    simpa using ENNReal.tendsto_ofReal h1
  have hev1 : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  have hevc : ∀ᶠ N in atTop, |(pN N : ℝ) / dN N - c| * (1 / z.im) ≤ T / 2 := by
    have h1 : Tendsto (fun N => |(pN N : ℝ) / dN N - c|) atTop (𝓝 0) := by
      simpa using (hcN.sub_const c).abs
    have h2 : Tendsto (fun N => |(pN N : ℝ) / dN N - c| * (1 / z.im)) atTop (𝓝 0) := by
      simpa using h1.mul_const (1 / z.im)
    exact (h2.eventually (gt_mem_nhds (by positivity : (0 : ℝ) < T / 2))).mono fun N h => h.le
  have hevc2 : ∀ᶠ N in atTop, (pN N : ℝ) / dN N ≤ c + 1 :=
    (hcN.eventually (gt_mem_nhds (by linarith : c < c + 1))).mono fun N h => h.le
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hBto
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  filter_upwards [hev1, hevc, hevc2] with N hdN hcNa hcNb
  have hcNnn : (0 : ℝ) ≤ (pN N : ℝ) / dN N := by positivity
  have hmset : MeasurableSet {Yv : Matrix (Fin (pN N)) (Fin (dN N)) ℝ | t ≤ resid Yv z} :=
    measurableSet_le measurable_const (measurable_resid z)
  have hsub : {ω | ε ≤ |‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖ - 0|}
      ⊆ {ω | t ≤ resid (Y N ω) z} := by
    intro ω hω
    by_contra hcon
    have hlt : resid (Y N ω) z < t := not_le.mp hcon
    have hmem : ε ≤ ‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖ := by
      have h' : ε ≤ |‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖ - 0| := hω
      rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at h'
    have h1 : ((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
        ≤ (c + 1) * (‖z‖ / z.im) * t := by
      have hstep1 : ((pN N : ℝ) / dN N) * (‖z‖ / z.im) ≤ (c + 1) * (‖z‖ / z.im) :=
        mul_le_mul_of_nonneg_right hcNb (by positivity)
      calc ((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
          ≤ (c + 1) * (‖z‖ / z.im) * resid (Y N ω) z :=
            mul_le_mul_of_nonneg_right hstep1 (resid_nonneg (Y N ω) z)
        _ ≤ (c + 1) * (‖z‖ / z.im) * t :=
            mul_le_mul_of_nonneg_left hlt.le (by positivity)
    have h2 : (c + 1) * (‖z‖ / z.im) * t = T / 2 := by
      rw [ht]
      field_simp
    have hG : ((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
        + |(pN N : ℝ) / dN N - c| * (1 / z.im) ≤ T := by linarith
    have hsmall : ((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
        + |(pN N : ℝ) / dN N - c| * (1 / z.im) ≤ 1 / 2 :=
      hG.trans (min_le_left _ _)
    have hkey := norm_stieltjesC_sub_mC_le (Y N ω) hc hz (hp N) hdN hsmall
    have hnum : ((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
        + |(pN N : ℝ) / dN N - c| * (1 / z.im)
        ≤ ε * (‖z‖ * (z.im * rootLb c z ^ 2)) / 2 := hG.trans (min_le_right _ _)
    have hfinal : ‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖ ≤ ε / 2 := by
      refine hkey.trans ?_
      calc (((pN N : ℝ) / dN N) * (‖z‖ / z.im) * resid (Y N ω) z
            + |(pN N : ℝ) / dN N - c| * (1 / z.im)) / (‖z‖ * (z.im * rootLb c z ^ 2))
          ≤ (ε * (‖z‖ * (z.im * rootLb c z ^ 2)) / 2)
              / (‖z‖ * (z.im * rootLb c z ^ 2)) := by gcongr
        _ = ε / 2 := by field_simp
    linarith
  calc μ N {ω | ε ≤ |‖R4C.stieltjesC (gram (Y N ω)) z - MP.mC c z‖ - 0|}
      ≤ μ N {ω | t ≤ resid (Y N ω) z} := measure_mono hsub
    _ = noiseMatrix ν (pN N) (dN N) {Yv | t ≤ resid Yv z} := (hY N).measure_eq hmset
    _ ≤ ENNReal.ofReal (residConst ν (dN N) z / t ^ 2) :=
        measure_resid_ge_le hν hz (hp N) hdN htpos


/-- **R1b at four moments.** The same for `d⁻¹ tr G²` against `MP.mCDeriv c z`. Twin of
`R1.tendstoInProb_stieltjes2C` (`RMT/R1.lean:910`), whose net argument is model free and is
followed here with the `GenRMT` copies of `RMT/General/Stability.lean`.

Proof: Cauchy's estimate at radius `η/2` (`norm_stieltjes2C_sub_mCDeriv_le`) turns a uniform
bound on the circle into a bound on the derivative gap; the circle is compact, so a finite net
of radius `ρ` covers it; the gap is `32/η²`-Lipschitz in `ζ` on the disc of radius `3η/4`
(`norm_diff_sub_diff_le`); and R1a at each of the finitely many net points, with a union
bound, sends the probability to `0`. -/
theorem tendstoInProb_stieltjes2C_general {c : ℝ} {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hc : 0 < c) (hz : 0 < z.im) (hd : Tendsto dN atTop atTop)
    (hp : ∀ N, 0 < pN N) (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c))
    (Y : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (dN N)) ℝ)
    (hY : ∀ N, HasLaw (Y N) (noiseMatrix ν (pN N) (dN N)) (μ N)) :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (gram (Y N ω)) z - MP.mCDeriv c z‖) 0 := by
  intro ε hε
  set r : ℝ := z.im / 2 with hrdef
  have hrpos : (0 : ℝ) < r := by rw [hrdef]; linarith
  set Lz : ℝ := 32 / z.im ^ 2 with hLzdef
  have hLzpos : (0 : ℝ) < Lz := by rw [hLzdef]; positivity
  set Csup : ℝ := ε * r / 2 with hCdef
  have hCpos : (0 : ℝ) < Csup := by rw [hCdef]; positivity
  set ρ : ℝ := Csup / (2 * Lz) with hrhodef
  have hrhopos : (0 : ℝ) < ρ := by rw [hrhodef]; positivity
  obtain ⟨b, hbsub, hbfin, hbcov⟩ :=
    (isCompact_sphere z r).elim_finite_subcover_image
      (b := Metric.sphere z r) (c := fun ζ : ℂ => Metric.ball ζ ρ)
      (fun ζ _ => Metric.isOpen_ball)
      (fun ζ hζ => Set.mem_biUnion hζ (Metric.mem_ball_self hrhopos))
  have hyim : ∀ y ∈ b, 0 < y.im := by
    intro y hy
    have hy' : y ∈ Metric.closedBall z r := Metric.sphere_subset_closedBall (hbsub hy)
    have := im_ge_of_mem_closedBall (z := z) hy'
    rw [hrdef] at this
    linarith
  have hyball : ∀ y ∈ b, y ∈ Metric.closedBall z (3 * z.im / 4) := by
    intro y hy
    have h1 : dist y z = r := hbsub hy
    rw [Metric.mem_closedBall, h1, hrdef]
    linarith
  have hpt : ∀ y ∈ hbfin.toFinset,
      Tendsto (fun N => μ N {ω | Csup / 2 ≤ ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖})
        atTop (𝓝 0) := by
    intro y hy
    have hyb : y ∈ b := hbfin.mem_toFinset.mp hy
    have hR1 := tendstoInProb_stieltjesC_general (z := y) hν hc (hyim y hyb) hd hp hcN Y hY
      (Csup / 2) (by positivity)
    have hset : ∀ N, {ω : Ω N | Csup / 2 ≤ ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖}
        = {ω | Csup / 2 ≤ |‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖ - 0|} := by
      intro N
      ext ω
      simp [abs_of_nonneg (norm_nonneg _)]
    simpa [hset] using hR1
  have hsum : Tendsto (fun N => ∑ y ∈ hbfin.toFinset,
      μ N {ω | Csup / 2 ≤ ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖}) atTop (𝓝 0) := by
    simpa using tendsto_finsetSum hbfin.toFinset hpt
  have hev1 : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hsum
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  filter_upwards [hev1] with N hdN
  have hsubset : {ω | ε ≤ |‖R4C.stieltjes2C (gram (Y N ω)) z - MP.mCDeriv c z‖ - 0|}
      ⊆ ⋃ y ∈ hbfin.toFinset,
          {ω | Csup / 2 ≤ ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖} := by
    intro ω hω
    by_contra hcon
    simp only [Set.mem_iUnion, Set.mem_ofPred_eq, not_exists, not_le] at hcon
    have hW := gram_isHermitian (Y N ω)
    have hC : ∀ ζ ∈ Metric.sphere z r,
        ‖R4C.stieltjesC (gram (Y N ω)) ζ - MP.mC c ζ‖ ≤ Csup := by
      intro ζ hζ
      obtain ⟨y, hyb, hyball'⟩ : ∃ y ∈ b, ζ ∈ Metric.ball y ρ := by
        have := hbcov hζ
        simpa using this
      have hlt : ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖ < Csup / 2 := by
        have := hcon y (hbfin.mem_toFinset.mpr hyb)
        simpa using this
      have hζball : ζ ∈ Metric.closedBall z (3 * z.im / 4) := by
        have h1 : dist ζ z = r := hζ
        rw [Metric.mem_closedBall, h1, hrdef]
        linarith
      have hlip := norm_diff_sub_diff_le (W := gram (Y N ω)) hc hz hdN hW hζball
        (hyball y hyb)
      have hdistlt : ‖ζ - y‖ < ρ := by
        have hd' : dist ζ y < ρ := hyball'
        rwa [Complex.dist_eq] at hd'
      have htri : ‖R4C.stieltjesC (gram (Y N ω)) ζ - MP.mC c ζ‖
          ≤ ‖(R4C.stieltjesC (gram (Y N ω)) ζ - MP.mC c ζ)
              - (R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y)‖
            + ‖R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y‖ := by
        have h := norm_add_le
          ((R4C.stieltjesC (gram (Y N ω)) ζ - MP.mC c ζ)
            - (R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y))
          (R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y)
        simpa using h
      have hLzrho : Lz * ρ = Csup / 2 := by
        rw [hrhodef]
        field_simp
      have hstep : Lz * ‖ζ - y‖ ≤ Lz * ρ := mul_le_mul_of_nonneg_left hdistlt.le hLzpos.le
      have hlip' : ‖(R4C.stieltjesC (gram (Y N ω)) ζ - MP.mC c ζ)
          - (R4C.stieltjesC (gram (Y N ω)) y - MP.mC c y)‖ ≤ Lz * ‖ζ - y‖ := by
        rw [hLzdef]; exact hlip
      linarith
    have hcauchy := norm_stieltjes2C_sub_mCDeriv_le (W := gram (Y N ω)) hc hz hdN hW hC
    have hmem : ε ≤ ‖R4C.stieltjes2C (gram (Y N ω)) z - MP.mCDeriv c z‖ := by
      have hω' : ε ≤ |‖R4C.stieltjes2C (gram (Y N ω)) z - MP.mCDeriv c z‖ - 0| := hω
      rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
    have hval : Csup / (z.im / 2) = ε / 2 := by
      rw [hCdef, hrdef]
      field_simp
    rw [← hrdef, hval] at hcauchy
    linarith
  exact (measure_mono hsubset).trans (measure_biUnion_finset_le _ _)

end Limits

end GenRMT
end StackedSVD
