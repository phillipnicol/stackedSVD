/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.ResolvDeriv
import StackedSVD.RMT.MP
import StackedSVD.Prob.GaussianAdapters

/-!
# The Gaussian integration by parts of item R1 (step 2 of `notes/archive/rmt_R1.md`)

`RMT/R1.lean` states the Stein equation in the residual form `norm_quad_integral_le`. This
file supplies the algebra and the analysis behind it.

Setting (the names of `RMT/ResolvDeriv.lean`). `x` is a point of the flattened space
`EuclideanSpace ℝ (Fin (p * d))`, `Y = (matrixEquivE p d).symm x`, `W₀ = d⁻¹ Yᵀ Y`,
`G = (W₀ - z)⁻¹`. Three complex matrices carry the whole computation:

* `Ymat x = Y` over `ℂ`,
* `Pmat z x = Y G` (`p × d`), the test function of the Gaussian integration by parts,
* `Qmat z x = Y G Yᵀ` (`p × p`).

Contents.

1. Entry formulas in the eigenbasis of `W₀` and the three uniform entry bounds
   `‖G_{bi}‖ ≤ d/η`, `‖P_{ka}‖ ≤ d √(d(η+‖z‖))/η`, `‖Q_{ka}‖ ≤ d²(η+‖z‖)/η`.
2. `ContDiff ℝ n` of `F_{ki}(x) = P_{ki}`, its derivative
   `∂F_{ki}(h) = (H G)_{ki} - d⁻¹ ((Q H G)_{ki} + (P Hᵀ P)_{ki})`, and the global bound
   `‖fderiv (F_{ki}) x‖ ≤ Lconst z p d`.
3. The two trace identities `∑_{k,i} Y_{ki} P_{ki} = d (d + z tr G)` and
   `∑_{k,i} ∂F_{ki}/∂Y_{ki} = p tr G - (d + z tr G) tr G - tr G - z tr G²`.
4. Stein's identity entry by entry, summed, and the residual bound of step 2.

Everything in sections 1 to 3 is deterministic. No `axiom`, no `sorry`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix NNReal

namespace StackedSVD

namespace SteinStep

open ResolvDeriv R4

variable {p d : ℕ} {z : ℂ}

/-! ### The three matrices -/

/-- `Y` read over `ℂ`. -/
noncomputable def Ymat (x : EuclideanSpace ℝ (Fin (p * d))) : Matrix (Fin p) (Fin d) ℂ :=
  cmapR ((matrixEquivE p d).symm x)

/-- `P = Y G`, the test function of the Gaussian integration by parts. -/
noncomputable def Pmat (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin p) (Fin d) ℂ :=
  Ymat x * Gmat z x

/-- `Q = Y G Yᵀ`. -/
noncomputable def Qmat (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin p) (Fin p) ℂ :=
  Pmat z x * (Ymat x)ᵀ

theorem Ymat_apply (x : EuclideanSpace ℝ (Fin (p * d))) (k : Fin p) (j : Fin d) :
    Ymat x k j = (((matrixEquivE p d).symm x k j : ℝ) : ℂ) := rfl

/-! ### The eigenbasis: entry formulas -/

/-- `G = U diag(g) Uᵀ`, in one entry. -/
theorem Gmat_entry_eq (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) (b i : Fin d) :
    Gmat z x b i
      = ∑ c : Fin d, ((eigUx x b c : ℝ) : ℂ) * gval z x c * ((eigUx x i c : ℝ) : ℂ) := by
  have hG : Gmat z x
      = R4C.cmat (eigUx x) * Matrix.diagonal (gval z x) * (R4C.cmat (eigUx x))ᵀ :=
    R4C.resolvC_eq_conj (isHermitian_gram _) hz
  rw [hG, Matrix.mul_apply]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [Matrix.mul_diagonal, Matrix.transpose_apply]
  rfl

private theorem Pmat_eq_conj (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Pmat z x = cmapR (Zmat x) * Matrix.diagonal (gval z x) * (R4C.cmat (eigUx x))ᵀ := by
  have hG : Gmat z x
      = R4C.cmat (eigUx x) * Matrix.diagonal (gval z x) * (R4C.cmat (eigUx x))ᵀ :=
    R4C.resolvC_eq_conj (isHermitian_gram _) hz
  have hYU : Ymat (p := p) (d := d) x * R4C.cmat (eigUx x) = cmapR (Zmat x) := by
    rw [Ymat, cmat_eq_cmapR, ← cmapR_mul]
    rfl
  rw [Pmat, hG, ← Matrix.mul_assoc, ← Matrix.mul_assoc, hYU]

/-- `P = Z diag(g) Uᵀ`, in one entry, with `Z = Y U`. -/
theorem Pmat_entry_eq (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) (k : Fin p)
    (a : Fin d) :
    Pmat z x k a
      = ∑ c : Fin d, ((Zmat x k c : ℝ) : ℂ) * gval z x c * ((eigUx x a c : ℝ) : ℂ) := by
  rw [Pmat_eq_conj hz x, Matrix.mul_apply]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [Matrix.mul_diagonal, Matrix.transpose_apply]
  rfl

/-- `Q = Z diag(g) Zᵀ`, in one entry. -/
theorem Qmat_entry_eq (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) (k a : Fin p) :
    Qmat z x k a
      = ∑ c : Fin d, ((Zmat x k c : ℝ) : ℂ) * gval z x c * ((Zmat x a c : ℝ) : ℂ) := by
  have hQ : Qmat z x = cmapR (Zmat x) * Matrix.diagonal (gval z x) * (cmapR (Zmat x))ᵀ := by
    have hYU : Ymat (p := p) (d := d) x * R4C.cmat (eigUx x) = cmapR (Zmat x) := by
      rw [Ymat, cmat_eq_cmapR, ← cmapR_mul]
      rfl
    rw [Qmat, Pmat_eq_conj hz x, Matrix.mul_assoc, ← Matrix.transpose_mul, hYU]
  rw [hQ, Matrix.mul_apply]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [Matrix.mul_diagonal, Matrix.transpose_apply]
  rfl

/-! ### The eigenbasis: scalar bounds

`sum_sq_Zmat` is `private` in `RMT/ResolvDeriv.lean`; the proof below is the same. -/

theorem sum_sq_Zmat' (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) (a : Fin d) :
    ∑ k : Fin p, (Zmat x k a) ^ 2 = (d : ℝ) * eigVal x a := by
  set Y := (matrixEquivE p d).symm x with hY
  set U := eigUx x with hU
  have hd0 : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hYY : Yᵀ * Y = (d : ℝ) • gram Y := by
    rw [gram, smul_smul, mul_inv_cancel₀ hd0, one_smul]
  have hUW : Uᵀ * gram Y * U = Matrix.diagonal (eigVal x) := by
    have hc := eigU_conj (isHermitian_gram Y)
    calc Uᵀ * gram Y * U
        = Uᵀ * (U * Matrix.diagonal (eigVal x) * Uᵀ) * U := by
          rw [hU, eigUx, hY]
          rw [show Matrix.diagonal (eigVal x)
              = Matrix.diagonal (isHermitian_gram ((matrixEquivE p d).symm x)).eigenvalues from
            rfl]
          rw [hc]
      _ = (Uᵀ * U) * Matrix.diagonal (eigVal x) * (Uᵀ * U) := by
          simp only [Matrix.mul_assoc]
      _ = Matrix.diagonal (eigVal x) := by
          rw [hU, eigUx, transpose_eigU_mul, Matrix.one_mul, Matrix.mul_one]
  have hZZ : (Zmat x)ᵀ * Zmat x = (d : ℝ) • Matrix.diagonal (eigVal x) := by
    rw [Zmat, Matrix.transpose_mul, ← hU, ← hY]
    calc Uᵀ * Yᵀ * (Y * U) = Uᵀ * (Yᵀ * Y) * U := by simp only [Matrix.mul_assoc]
      _ = Uᵀ * ((d : ℝ) • gram Y) * U := by rw [hYY]
      _ = (d : ℝ) • (Uᵀ * gram Y * U) := by rw [Matrix.mul_smul, Matrix.smul_mul]
      _ = (d : ℝ) • Matrix.diagonal (eigVal x) := by rw [hUW]
  have h1 : ((Zmat x)ᵀ * Zmat x) a a = ∑ k : Fin p, (Zmat x k a) ^ 2 := by
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun k _ => by rw [Matrix.transpose_apply, sq]
  rw [← h1, hZZ, Matrix.smul_apply, Matrix.diagonal_apply_eq, smul_eq_mul]

theorem sq_Zmat_le (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) (k : Fin p) (a : Fin d) :
    (Zmat x k a) ^ 2 ≤ (d : ℝ) * eigVal x a := by
  rw [← sum_sq_Zmat' hd x a]
  exact Finset.single_le_sum (f := fun k : Fin p => (Zmat x k a) ^ 2)
    (fun j _ => sq_nonneg _) (Finset.mem_univ k)

theorem eigVal_nonneg (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) (a : Fin d) :
    0 ≤ eigVal x a := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have h0 : (0 : ℝ) ≤ ∑ k : Fin p, (Zmat x k a) ^ 2 :=
    Finset.sum_nonneg fun k _ => sq_nonneg _
  rw [sum_sq_Zmat' hd x a] at h0
  nlinarith

theorem abs_eigUx_le_one (x : EuclideanSpace ℝ (Fin (p * d))) (j c : Fin d) :
    |eigUx x j c| ≤ 1 := by
  have hrow : ∑ b : Fin d, (eigUx x j b) ^ 2 = 1 := by
    have h := eigU_mul_transpose (isHermitian_gram ((matrixEquivE p d).symm x))
    have h1 : (eigUx x * (eigUx x)ᵀ) j j = (1 : Matrix (Fin d) (Fin d) ℝ) j j := by
      rw [eigUx]; rw [h]
    rw [Matrix.mul_apply] at h1
    rw [Matrix.one_apply_eq] at h1
    rw [← h1]
    exact Finset.sum_congr rfl fun b _ => by rw [Matrix.transpose_apply, sq]
  have hle : (eigUx x j c) ^ 2 ≤ 1 := by
    rw [← hrow]
    exact Finset.single_le_sum (f := fun b : Fin d => (eigUx x j b) ^ 2)
      (fun b _ => sq_nonneg _) (Finset.mem_univ c)
  nlinarith [abs_nonneg (eigUx x j c), sq_abs (eigUx x j c)]

/-! ### The three scalar estimates on `g_c` -/

private theorem norm_sub_ge (hz : 0 < z.im) (l : ℝ) : z.im ≤ ‖(l : ℂ) - z‖ := by
  have h := Complex.abs_im_le_norm ((l : ℂ) - z)
  have him : ((l : ℂ) - z).im = -z.im := by simp
  rw [him, abs_neg, abs_of_pos hz] at h
  exact h

theorem norm_gval_le (hz : 0 < z.im) (x : EuclideanSpace ℝ (Fin (p * d))) (c : Fin d) :
    ‖gval z x c‖ ≤ (z.im)⁻¹ := by
  rw [gval, norm_inv]
  exact inv_anti₀ hz (norm_sub_ge hz _)

theorem eigVal_mul_norm_gval_le (hz : 0 < z.im) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) (c : Fin d) :
    eigVal x c * ‖gval z x c‖ ≤ (z.im + ‖z‖) / z.im := by
  set l : ℝ := eigVal x c with hl
  set D : ℝ := ‖((l : ℂ) - z)‖ with hD
  have hDge : z.im ≤ D := norm_sub_ge hz l
  have hDpos : 0 < D := lt_of_lt_of_le hz hDge
  have hgn : ‖gval z x c‖ = D⁻¹ := by rw [gval, norm_inv, ← hD]
  have hlnn : 0 ≤ l := eigVal_nonneg hd x c
  have hlle : l ≤ D + ‖z‖ := by
    have h1 : ‖(l : ℂ)‖ ≤ ‖((l : ℂ) - z)‖ + ‖z‖ := by
      have := norm_add_le ((l : ℂ) - z) z
      simpa using this
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hlnn] at h1
    exact h1
  rw [hgn]
  have hinv : D⁻¹ ≤ (z.im)⁻¹ := inv_anti₀ hz hDge
  have h2 : l * D⁻¹ ≤ (D + ‖z‖) * D⁻¹ := mul_le_mul_of_nonneg_right hlle (by positivity)
  have h3 : (D + ‖z‖) * D⁻¹ = 1 + ‖z‖ * D⁻¹ := by field_simp
  have h4 : ‖z‖ * D⁻¹ ≤ ‖z‖ * (z.im)⁻¹ := mul_le_mul_of_nonneg_left hinv (norm_nonneg z)
  have h5 : 1 + ‖z‖ * (z.im)⁻¹ = (z.im + ‖z‖) / z.im := by field_simp
  linarith [h2, h3.le, h3.ge, h4, h5.le, h5.ge]

theorem eigVal_mul_norm_gval_sq_le (hz : 0 < z.im) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) (c : Fin d) :
    eigVal x c * ‖gval z x c‖ ^ 2 ≤ (z.im + ‖z‖) / z.im ^ 2 := by
  set l : ℝ := eigVal x c with hl
  set D : ℝ := ‖((l : ℂ) - z)‖ with hD
  have hDge : z.im ≤ D := norm_sub_ge hz l
  have hDpos : 0 < D := lt_of_lt_of_le hz hDge
  have hgn : ‖gval z x c‖ = D⁻¹ := by rw [gval, norm_inv, ← hD]
  have hlnn : 0 ≤ l := eigVal_nonneg hd x c
  have hlle : l ≤ D + ‖z‖ := by
    have h1 : ‖(l : ℂ)‖ ≤ ‖((l : ℂ) - z)‖ + ‖z‖ := by
      have := norm_add_le ((l : ℂ) - z) z
      simpa using this
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hlnn] at h1
    exact h1
  rw [hgn]
  have hinv : D⁻¹ ≤ (z.im)⁻¹ := inv_anti₀ hz hDge
  have hinv2 : (D⁻¹) ^ 2 ≤ ((z.im)⁻¹) ^ 2 :=
    pow_le_pow_left₀ (by positivity) hinv 2
  have h2 : l * (D⁻¹) ^ 2 ≤ (D + ‖z‖) * (D⁻¹) ^ 2 :=
    mul_le_mul_of_nonneg_right hlle (by positivity)
  have h3 : (D + ‖z‖) * (D⁻¹) ^ 2 = D⁻¹ + ‖z‖ * (D⁻¹) ^ 2 := by field_simp
  have h4 : ‖z‖ * (D⁻¹) ^ 2 ≤ ‖z‖ * ((z.im)⁻¹) ^ 2 :=
    mul_le_mul_of_nonneg_left hinv2 (norm_nonneg z)
  have h5 : (z.im)⁻¹ + ‖z‖ * ((z.im)⁻¹) ^ 2 = (z.im + ‖z‖) / z.im ^ 2 := by field_simp
  linarith [h2, h3.le, h3.ge, h4, h5.le, h5.ge, hinv]

/-! ### The three uniform entry bounds -/

/-- `√(d(η+‖z‖))/η`, the bound on `|Z_{kc}| ‖g_c‖`. -/
noncomputable def sB (z : ℂ) (d : ℕ) : ℝ := Real.sqrt ((d : ℝ) * (z.im + ‖z‖)) / z.im

/-- The uniform bound on one entry of `G`. -/
noncomputable def BG (z : ℂ) (d : ℕ) : ℝ := (d : ℝ) / z.im

/-- The uniform bound on one entry of `P = Y G`. -/
noncomputable def BP (z : ℂ) (d : ℕ) : ℝ := (d : ℝ) * sB z d

/-- The uniform bound on one entry of `Q = Y G Yᵀ`. -/
noncomputable def BQ (z : ℂ) (d : ℕ) : ℝ := (d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im

theorem sB_nonneg (hz : 0 < z.im) (d : ℕ) : 0 ≤ sB z d := by
  have h1 : 0 ≤ (d : ℝ) * (z.im + ‖z‖) := by positivity
  unfold sB
  positivity

theorem BG_nonneg (hz : 0 < z.im) (d : ℕ) : 0 ≤ BG z d := by
  unfold BG; positivity

theorem BP_nonneg (hz : 0 < z.im) (d : ℕ) : 0 ≤ BP z d := by
  have := sB_nonneg hz d
  unfold BP; positivity

theorem BQ_nonneg (hz : 0 < z.im) (d : ℕ) : 0 ≤ BQ z d := by
  have h1 : 0 ≤ z.im + ‖z‖ := by positivity
  unfold BQ; positivity

theorem sq_sB (hz : 0 < z.im) (d : ℕ) :
    sB z d ^ 2 = (d : ℝ) * (z.im + ‖z‖) / z.im ^ 2 := by
  have h1 : (0 : ℝ) ≤ (d : ℝ) * (z.im + ‖z‖) := by positivity
  unfold sB
  rw [div_pow, Real.sq_sqrt h1]

private theorem abs_Zmat_mul_norm_gval_le (hz : 0 < z.im) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) (k : Fin p) (c : Fin d) :
    |Zmat x k c| * ‖gval z x c‖ ≤ sB z d := by
  have hdR : (0 : ℝ) ≤ (d : ℝ) := Nat.cast_nonneg d
  have hnn : 0 ≤ |Zmat x k c| * ‖gval z x c‖ := by positivity
  have hsq : (|Zmat x k c| * ‖gval z x c‖) ^ 2 ≤ sB z d ^ 2 := by
    have h1 : (|Zmat x k c| * ‖gval z x c‖) ^ 2
        = (Zmat x k c) ^ 2 * ‖gval z x c‖ ^ 2 := by
      rw [mul_pow, sq_abs]
    have h2 : (Zmat x k c) ^ 2 * ‖gval z x c‖ ^ 2
        ≤ ((d : ℝ) * eigVal x c) * ‖gval z x c‖ ^ 2 :=
      mul_le_mul_of_nonneg_right (sq_Zmat_le hd x k c) (by positivity)
    have h3 : ((d : ℝ) * eigVal x c) * ‖gval z x c‖ ^ 2
        ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) := by
      have := eigVal_mul_norm_gval_sq_le hz hd x c
      nlinarith [hdR]
    rw [h1, sq_sB hz d]
    calc (Zmat x k c) ^ 2 * ‖gval z x c‖ ^ 2
        ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im ^ 2) := h2.trans h3
      _ = (d : ℝ) * (z.im + ‖z‖) / z.im ^ 2 := by ring
  calc |Zmat x k c| * ‖gval z x c‖
      = Real.sqrt ((|Zmat x k c| * ‖gval z x c‖) ^ 2) := (Real.sqrt_sq hnn).symm
    _ ≤ Real.sqrt (sB z d ^ 2) := Real.sqrt_le_sqrt hsq
    _ = sB z d := Real.sqrt_sq (sB_nonneg hz d)

private theorem abs_Zmat_mul_abs_Zmat_le (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) (k a : Fin p) (c : Fin d) :
    |Zmat x k c| * |Zmat x a c| ≤ (d : ℝ) * eigVal x c := by
  have h1 := sq_Zmat_le hd x k c
  have h2 := sq_Zmat_le hd x a c
  nlinarith [sq_nonneg (|Zmat x k c| - |Zmat x a c|), sq_abs (Zmat x k c),
    sq_abs (Zmat x a c), abs_nonneg (Zmat x k c), abs_nonneg (Zmat x a c)]

/-- A finite sum of complex numbers with a uniform entry bound. -/
private theorem norm_sum_le_card {ι : Type*} [Fintype ι] (f : ι → ℂ) {B : ℝ}
    (hf : ∀ a, ‖f a‖ ≤ B) : ‖∑ a, f a‖ ≤ (Fintype.card ι : ℝ) * B := by
  calc ‖∑ a, f a‖ ≤ ∑ a, ‖f a‖ := norm_sum_le _ _
    _ ≤ ∑ _a : ι, B := Finset.sum_le_sum fun a _ => hf a
    _ = (Fintype.card ι : ℝ) * B := by
        rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

theorem norm_Gmat_le (hz : 0 < z.im) (x : EuclideanSpace ℝ (Fin (p * d))) (b i : Fin d) :
    ‖Gmat z x b i‖ ≤ BG z d := by
  rw [Gmat_entry_eq hz.ne' x b i]
  have hterm : ∀ c : Fin d,
      ‖((eigUx x b c : ℝ) : ℂ) * gval z x c * ((eigUx x i c : ℝ) : ℂ)‖ ≤ (z.im)⁻¹ := by
    intro c
    rw [norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
      Real.norm_eq_abs]
    have h1 := abs_eigUx_le_one x b c
    have h2 := abs_eigUx_le_one x i c
    have h3 := norm_gval_le hz x c
    have h4 : (0 : ℝ) ≤ ‖gval z x c‖ := norm_nonneg _
    have hA : |eigUx x b c| * ‖gval z x c‖ ≤ ‖gval z x c‖ := by
      nlinarith [abs_nonneg (eigUx x b c)]
    have hB : (|eigUx x b c| * ‖gval z x c‖) * |eigUx x i c| ≤ ‖gval z x c‖ * 1 :=
      mul_le_mul hA h2 (abs_nonneg _) h4
    linarith
  have h := norm_sum_le_card
    (fun c : Fin d => ((eigUx x b c : ℝ) : ℂ) * gval z x c * ((eigUx x i c : ℝ) : ℂ)) hterm
  rw [Fintype.card_fin] at h
  rw [BG]
  calc ‖∑ c : Fin d, ((eigUx x b c : ℝ) : ℂ) * gval z x c * ((eigUx x i c : ℝ) : ℂ)‖
      ≤ (d : ℝ) * (z.im)⁻¹ := h
    _ = (d : ℝ) / z.im := by field_simp

theorem norm_Pmat_le (hz : 0 < z.im) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d)))
    (k : Fin p) (a : Fin d) : ‖Pmat z x k a‖ ≤ BP z d := by
  rw [Pmat_entry_eq hz.ne' x k a]
  have hterm : ∀ c : Fin d,
      ‖((Zmat x k c : ℝ) : ℂ) * gval z x c * ((eigUx x a c : ℝ) : ℂ)‖ ≤ sB z d := by
    intro c
    rw [norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
      Real.norm_eq_abs]
    have h1 := abs_Zmat_mul_norm_gval_le hz hd x k c
    have h2 := abs_eigUx_le_one x a c
    have hB : (|Zmat x k c| * ‖gval z x c‖) * |eigUx x a c| ≤ sB z d * 1 :=
      mul_le_mul h1 h2 (abs_nonneg _) (sB_nonneg hz d)
    linarith
  have h := norm_sum_le_card
    (fun c : Fin d => ((Zmat x k c : ℝ) : ℂ) * gval z x c * ((eigUx x a c : ℝ) : ℂ)) hterm
  rw [Fintype.card_fin] at h
  exact h

theorem norm_Qmat_le (hz : 0 < z.im) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d)))
    (k a : Fin p) : ‖Qmat z x k a‖ ≤ BQ z d := by
  rw [Qmat_entry_eq hz.ne' x k a]
  have hterm : ∀ c : Fin d,
      ‖((Zmat x k c : ℝ) : ℂ) * gval z x c * ((Zmat x a c : ℝ) : ℂ)‖
        ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im) := by
    intro c
    rw [norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
      Real.norm_eq_abs]
    have h1 := abs_Zmat_mul_abs_Zmat_le hd x k a c
    have h2 := eigVal_mul_norm_gval_le hz hd x c
    have h3 : (0 : ℝ) ≤ ‖gval z x c‖ := norm_nonneg _
    have h4 : (0 : ℝ) ≤ (d : ℝ) := Nat.cast_nonneg d
    have h5 : |Zmat x k c| * ‖gval z x c‖ * |Zmat x a c|
        = (|Zmat x k c| * |Zmat x a c|) * ‖gval z x c‖ := by ring
    rw [h5]
    calc (|Zmat x k c| * |Zmat x a c|) * ‖gval z x c‖
        ≤ ((d : ℝ) * eigVal x c) * ‖gval z x c‖ :=
          mul_le_mul_of_nonneg_right h1 h3
      _ = (d : ℝ) * (eigVal x c * ‖gval z x c‖) := by ring
      _ ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im) := mul_le_mul_of_nonneg_left h2 h4
  have h := norm_sum_le_card
    (fun c : Fin d => ((Zmat x k c : ℝ) : ℂ) * gval z x c * ((Zmat x a c : ℝ) : ℂ)) hterm
  rw [Fintype.card_fin] at h
  rw [BQ]
  calc ‖∑ c : Fin d, ((Zmat x k c : ℝ) : ℂ) * gval z x c * ((Zmat x a c : ℝ) : ℂ)‖
      ≤ (d : ℝ) * ((d : ℝ) * ((z.im + ‖z‖) / z.im)) := h
    _ = (d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im := by ring

/-- One entry of `Y` is bounded by the Frobenius norm of the flattened point. -/
theorem norm_Ymat_le (h : EuclideanSpace ℝ (Fin (p * d))) (k : Fin p) (j : Fin d) :
    ‖Ymat h k j‖ ≤ ‖h‖ := by
  set H := (matrixEquivE p d).symm h with hH
  have hnorm : ‖h‖ ^ 2 = ∑ a : Fin p, ∑ b : Fin d, (H a b) ^ 2 := by
    have h1 := norm_matrixEquivE_sq (p := p) (d := d) H
    rw [hH, MeasurableEquiv.apply_symm_apply] at h1
    exact h1
  have hle : (H k j) ^ 2 ≤ ‖h‖ ^ 2 := by
    rw [hnorm]
    calc (H k j) ^ 2 ≤ ∑ b : Fin d, (H k b) ^ 2 :=
          Finset.single_le_sum (f := fun b : Fin d => (H k b) ^ 2)
            (fun b _ => sq_nonneg _) (Finset.mem_univ j)
      _ ≤ ∑ a : Fin p, ∑ b : Fin d, (H a b) ^ 2 :=
          Finset.single_le_sum (f := fun a : Fin p => ∑ b : Fin d, (H a b) ^ 2)
            (fun a _ => Finset.sum_nonneg fun b _ => sq_nonneg _) (Finset.mem_univ k)
  have hYn : ‖Ymat h k j‖ = |H k j| := by
    rw [Ymat_apply, Complex.norm_real, Real.norm_eq_abs, hH]
  rw [hYn]
  nlinarith [abs_nonneg (H k j), sq_abs (H k j), norm_nonneg h]

/-! ### The test function `F_{ki} = (Y G)_{ki}` and its derivative -/

/-- The coordinate functional of the flattened space, complexified. -/
noncomputable def coordC (p d : ℕ) (r : Fin (p * d)) :=
  Complex.ofRealCLM.comp (EuclideanSpace.proj (𝕜 := ℝ) r)

@[simp] theorem coordC_apply (r : Fin (p * d)) (x : EuclideanSpace ℝ (Fin (p * d))) :
    coordC p d r x = ((x r : ℝ) : ℂ) := rfl

theorem contDiff_coordC {n : WithTop ℕ∞} (r : Fin (p * d)) :
    ContDiff ℝ n (fun x : EuclideanSpace ℝ (Fin (p * d)) => ((x r : ℝ) : ℂ)) :=
  (coordC p d r).contDiff

theorem hasFDerivAt_coordC (r : Fin (p * d)) (x : EuclideanSpace ℝ (Fin (p * d))) :
    HasFDerivAt (fun y : EuclideanSpace ℝ (Fin (p * d)) => ((y r : ℝ) : ℂ))
      (coordC p d r) x :=
  (coordC p d r).hasFDerivAt

/-- `F_{ki}(x) = (Y G)_{ki}`, the test function of the Gaussian integration by parts. -/
noncomputable def Fentry (z : ℂ) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) : ℂ := Pmat z x k i

theorem Fentry_eq_sum (z : ℂ) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    Fentry z k i x = ∑ j : Fin d, Ymat x k j * Gmat z x j i := by
  rw [Fentry, Pmat, Matrix.mul_apply]

private theorem Fentry_eq_fun (z : ℂ) (k : Fin p) (i : Fin d) :
    Fentry (p := p) (d := d) z k i
      = fun y : EuclideanSpace ℝ (Fin (p * d)) =>
        ∑ j : Fin d, ((y (finProdFinEquiv (k, j)) : ℝ) : ℂ) * Gmat z y j i := by
  funext y
  rw [Fentry_eq_sum]
  rfl

theorem contDiff_Fentry (hz : z.im ≠ 0) {n : WithTop ℕ∞} (k : Fin p) (i : Fin d) :
    ContDiff ℝ n (Fentry (p := p) (d := d) z k i) := by
  rw [Fentry_eq_fun]
  exact ContDiff.sum fun j _ => (contDiff_coordC _).mul (contDiff_Gmat_entry hz j i)

private theorem Fentry_deriv_aux (hz : z.im ≠ 0) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (Fentry (p := p) (d := d) z k i) x ∧
      ∀ h : EuclideanSpace ℝ (Fin (p * d)),
        fderiv ℝ (Fentry (p := p) (d := d) z k i) x h
          = (Ymat h * Gmat z x) k i
            - (Ymat x * (Gmat z x * dWmat p d x h * Gmat z x)) k i := by
  have hS : HasFDerivAt (fun y : EuclideanSpace ℝ (Fin (p * d)) =>
      ∑ j : Fin d, ((y (finProdFinEquiv (k, j)) : ℝ) : ℂ) * Gmat z y j i)
      (∑ j : Fin d, (((x (finProdFinEquiv (k, j)) : ℝ) : ℂ) •
          fderiv ℝ (fun y : EuclideanSpace ℝ (Fin (p * d)) => Gmat z y j i) x
        + (Gmat z x j i) • coordC p d (finProdFinEquiv (k, j)))) x :=
    HasFDerivAt.fun_sum fun j _ =>
      (hasFDerivAt_coordC (finProdFinEquiv (k, j)) x).mul
        ((differentiableAt_Gmat_entry hz j i x).hasFDerivAt)
  rw [Fentry_eq_fun]
  refine ⟨hS.differentiableAt, fun h => ?_⟩
  rw [hS.fderiv, sum_apply]
  simp only [add_apply, smul_apply, smul_eq_mul, coordC_apply]
  rw [Matrix.mul_apply, Matrix.mul_apply, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [fderiv_Gmat_entry hz j i x h]
  simp only [Ymat_apply, matrixEquivE_symm_apply]
  ring

theorem differentiableAt_Fentry (hz : z.im ≠ 0) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (Fentry (p := p) (d := d) z k i) x :=
  (Fentry_deriv_aux hz k i x).1

/-- **The derivative of the test function.** -/
theorem fderiv_Fentry (hz : z.im ≠ 0) (k : Fin p) (i : Fin d)
    (x h : EuclideanSpace ℝ (Fin (p * d))) :
    fderiv ℝ (Fentry (p := p) (d := d) z k i) x h
      = (Ymat h * Gmat z x) k i
        - (Ymat x * (Gmat z x * dWmat p d x h * Gmat z x)) k i :=
  (Fentry_deriv_aux hz k i x).2 h

private theorem Ymat_mul_G_dW_G (x h : EuclideanSpace ℝ (Fin (p * d))) :
    Ymat x * (Gmat z x * dWmat p d x h * Gmat z x)
      = (d : ℂ)⁻¹ • (Qmat z x * Ymat h * Gmat z x)
        + (d : ℂ)⁻¹ • (Pmat z x * (Ymat h)ᵀ * Pmat z x) := by
  have hdW : dWmat p d x h
      = (d : ℂ)⁻¹ • ((Ymat x)ᵀ * Ymat h) + (d : ℂ)⁻¹ • ((Ymat h)ᵀ * Ymat x) := rfl
  rw [hdW, Qmat, Pmat]
  simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul,
    Matrix.mul_assoc]

/-- **The derivative of the test function, in `P` and `Q`.** -/
theorem fderiv_Fentry_eq (hz : z.im ≠ 0) (k : Fin p) (i : Fin d)
    (x h : EuclideanSpace ℝ (Fin (p * d))) :
    fderiv ℝ (Fentry (p := p) (d := d) z k i) x h
      = (Ymat h * Gmat z x) k i
        - (d : ℂ)⁻¹ * ((Qmat z x * Ymat h * Gmat z x) k i
            + (Pmat z x * (Ymat h)ᵀ * Pmat z x) k i) := by
  rw [fderiv_Fentry hz k i x h, Ymat_mul_G_dW_G x h]
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  ring

/-! ### The global bound on the derivative -/

/-- The global bound on `‖fderiv (F_{ki})‖`. Its value does not matter: it enters
`integral_entry_mul_gaussianMatrix_complex` only as a hypothesis. -/
noncomputable def Lconst (z : ℂ) (p d : ℕ) : ℝ :=
  (d : ℝ) * BG z d + (p : ℝ) * (d : ℝ) * (BQ z d * BG z d)
    + (p : ℝ) * (d : ℝ) * (BP z d * BP z d)

theorem Lconst_nonneg (hz : 0 < z.im) (p d : ℕ) : 0 ≤ Lconst z p d := by
  have h1 := BG_nonneg hz d
  have h2 := BP_nonneg hz d
  have h3 := BQ_nonneg hz d
  unfold Lconst
  positivity

theorem norm_fderiv_Fentry_le (hz : 0 < z.im) (hd : 0 < d) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ‖fderiv ℝ (Fentry (p := p) (d := d) z k i) x‖ ≤ Lconst z p d := by
  refine ContinuousLinearMap.opNorm_le_bound _ (Lconst_nonneg hz p d) fun h => ?_
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hbg := BG_nonneg hz (d := d)
  have hbp := BP_nonneg hz (d := d)
  have hbq := BQ_nonneg hz (d := d)
  have hnh : (0 : ℝ) ≤ ‖h‖ := norm_nonneg h
  -- the three entry bounds
  have hY : ∀ (a : Fin p) (b : Fin d), ‖Ymat h a b‖ ≤ ‖h‖ := fun a b => norm_Ymat_le h a b
  have hG : ∀ (a b : Fin d), ‖Gmat z x a b‖ ≤ BG z d := fun a b => norm_Gmat_le hz x a b
  have hP : ∀ (a : Fin p) (b : Fin d), ‖Pmat z x a b‖ ≤ BP z d :=
    fun a b => norm_Pmat_le hz hd x a b
  have hQ : ∀ a b : Fin p, ‖Qmat z x a b‖ ≤ BQ z d := fun a b => norm_Qmat_le hz hd x a b
  -- term 1
  have t1 : ‖(Ymat h * Gmat z x) k i‖ ≤ (d : ℝ) * (‖h‖ * BG z d) := by
    rw [Matrix.mul_apply]
    have hb : ∀ j : Fin d, ‖Ymat h k j * Gmat z x j i‖ ≤ ‖h‖ * BG z d := by
      intro j
      rw [norm_mul]
      exact mul_le_mul (hY k j) (hG j i) (norm_nonneg _) hnh
    have := norm_sum_le_card (fun j : Fin d => Ymat h k j * Gmat z x j i) hb
    rwa [Fintype.card_fin] at this
  -- term 2
  have t2row : ∀ a : Fin d, ‖(Qmat z x * Ymat h) k a‖ ≤ (p : ℝ) * (BQ z d * ‖h‖) := by
    intro a
    rw [Matrix.mul_apply]
    have hb : ∀ b : Fin p, ‖Qmat z x k b * Ymat h b a‖ ≤ BQ z d * ‖h‖ := by
      intro b
      rw [norm_mul]
      exact mul_le_mul (hQ k b) (hY b a) (norm_nonneg _) hbq
    have := norm_sum_le_card (fun b : Fin p => Qmat z x k b * Ymat h b a) hb
    rwa [Fintype.card_fin] at this
  have t2 : ‖(Qmat z x * Ymat h * Gmat z x) k i‖
      ≤ (d : ℝ) * ((p : ℝ) * (BQ z d * ‖h‖) * BG z d) := by
    rw [Matrix.mul_apply]
    have hb : ∀ a : Fin d, ‖(Qmat z x * Ymat h) k a * Gmat z x a i‖
        ≤ (p : ℝ) * (BQ z d * ‖h‖) * BG z d := by
      intro a
      rw [norm_mul]
      exact mul_le_mul (t2row a) (hG a i) (norm_nonneg _) (by positivity)
    have := norm_sum_le_card
      (fun a : Fin d => (Qmat z x * Ymat h) k a * Gmat z x a i) hb
    rwa [Fintype.card_fin] at this
  -- term 3
  have t3row : ∀ a : Fin p, ‖(Pmat z x * (Ymat h)ᵀ) k a‖ ≤ (d : ℝ) * (BP z d * ‖h‖) := by
    intro a
    rw [Matrix.mul_apply]
    have hb : ∀ b : Fin d, ‖Pmat z x k b * (Ymat h)ᵀ b a‖ ≤ BP z d * ‖h‖ := by
      intro b
      rw [norm_mul, Matrix.transpose_apply]
      exact mul_le_mul (hP k b) (hY a b) (norm_nonneg _) hbp
    have := norm_sum_le_card (fun b : Fin d => Pmat z x k b * (Ymat h)ᵀ b a) hb
    rwa [Fintype.card_fin] at this
  have t3 : ‖(Pmat z x * (Ymat h)ᵀ * Pmat z x) k i‖
      ≤ (p : ℝ) * ((d : ℝ) * (BP z d * ‖h‖) * BP z d) := by
    rw [Matrix.mul_apply]
    have hb : ∀ a : Fin p, ‖(Pmat z x * (Ymat h)ᵀ) k a * Pmat z x a i‖
        ≤ (d : ℝ) * (BP z d * ‖h‖) * BP z d := by
      intro a
      rw [norm_mul]
      exact mul_le_mul (t3row a) (hP a i) (norm_nonneg _) (by positivity)
    have := norm_sum_le_card
      (fun a : Fin p => (Pmat z x * (Ymat h)ᵀ) k a * Pmat z x a i) hb
    rwa [Fintype.card_fin] at this
  -- the scalar factor
  have hdinv : ‖(d : ℂ)⁻¹‖ ≤ 1 := by
    have h1 : ‖(d : ℂ)⁻¹‖ = ((d : ℝ))⁻¹ := by
      rw [norm_inv, Complex.norm_natCast]
    have h2 : (1 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
    rw [h1]
    exact inv_le_one_of_one_le₀ h2
  rw [fderiv_Fentry_eq hz.ne' k i x h]
  have hstep : ‖(Ymat h * Gmat z x) k i
      - (d : ℂ)⁻¹ * ((Qmat z x * Ymat h * Gmat z x) k i
          + (Pmat z x * (Ymat h)ᵀ * Pmat z x) k i)‖
      ≤ ‖(Ymat h * Gmat z x) k i‖
        + ‖(d : ℂ)⁻¹‖ * (‖(Qmat z x * Ymat h * Gmat z x) k i‖
            + ‖(Pmat z x * (Ymat h)ᵀ * Pmat z x) k i‖) := by
    refine (norm_sub_le _ _).trans ?_
    rw [norm_mul]
    have h2 := norm_add_le ((Qmat z x * Ymat h * Gmat z x) k i)
      ((Pmat z x * (Ymat h)ᵀ * Pmat z x) k i)
    have h4 := mul_le_mul_of_nonneg_left h2 (norm_nonneg ((d : ℂ)⁻¹))
    linarith
  have hpos : (0 : ℝ) ≤ ‖(Qmat z x * Ymat h * Gmat z x) k i‖
      + ‖(Pmat z x * (Ymat h)ᵀ * Pmat z x) k i‖ := by positivity
  have hfin : ‖(d : ℂ)⁻¹‖ * (‖(Qmat z x * Ymat h * Gmat z x) k i‖
      + ‖(Pmat z x * (Ymat h)ᵀ * Pmat z x) k i‖)
      ≤ 1 * ((d : ℝ) * ((p : ℝ) * (BQ z d * ‖h‖) * BG z d)
          + (p : ℝ) * ((d : ℝ) * (BP z d * ‖h‖) * BP z d)) := by
    refine mul_le_mul hdinv (by linarith) hpos zero_le_one
  have heq : (d : ℝ) * (‖h‖ * BG z d)
      + 1 * ((d : ℝ) * ((p : ℝ) * (BQ z d * ‖h‖) * BG z d)
          + (p : ℝ) * ((d : ℝ) * (BP z d * ‖h‖) * BP z d))
      = Lconst z p d * ‖h‖ := by
    unfold Lconst
    ring
  linarith

/-! ### The derivative in the direction of its own entry -/

private theorem Ymat_single_eq (k a : Fin p) (i b : Fin d) :
    Ymat (p := p) (d := d) (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ)) a b
      = if a = k then (if b = i then 1 else 0) else 0 := by
  rw [Ymat_apply, matrixEquivE_symm_apply, PiLp.single_apply]
  by_cases hab : a = k
  · subst hab
    by_cases hb : b = i
    · subst hb; simp
    · simp [hb]
  · simp [hab]

theorem fderiv_Fentry_single (hz : z.im ≠ 0) (k : Fin p) (i : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    fderiv ℝ (Fentry (p := p) (d := d) z k i) x
        (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
      = Gmat z x i i
        - (d : ℂ)⁻¹ * (Qmat z x k k * Gmat z x i i + Pmat z x k i * Pmat z x k i) := by
  set h : EuclideanSpace ℝ (Fin (p * d)) :=
    EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ) with hh
  have hYs : ∀ (a : Fin p) (b : Fin d),
      Ymat (p := p) (d := d) h a b = if a = k then (if b = i then 1 else 0) else 0 := by
    intro a b; rw [hh]; exact Ymat_single_eq k a i b
  have e1 : (Ymat h * Gmat z x) k i = Gmat z x i i := by
    rw [Matrix.mul_apply]
    have hterm : ∀ j : Fin d, Ymat (p := p) (d := d) h k j * Gmat z x j i
        = if j = i then Gmat z x j i else 0 := by
      intro j
      rw [hYs k j]
      by_cases hj : j = i <;> simp [hj]
    rw [Finset.sum_congr rfl fun j _ => hterm j]
    simp
  have e2row : ∀ a : Fin d, (Qmat z x * Ymat h) k a = if a = i then Qmat z x k k else 0 := by
    intro a
    rw [Matrix.mul_apply]
    have hterm : ∀ b : Fin p, Qmat z x k b * Ymat (p := p) (d := d) h b a
        = if b = k then (if a = i then Qmat z x k b else 0) else 0 := by
      intro b
      rw [hYs b a]
      by_cases hb : b = k <;> by_cases ha : a = i <;> simp [hb, ha]
    rw [Finset.sum_congr rfl fun b _ => hterm b]
    simp
  have e2 : (Qmat z x * Ymat h * Gmat z x) k i = Qmat z x k k * Gmat z x i i := by
    rw [Matrix.mul_apply]
    have hterm : ∀ a : Fin d, (Qmat z x * Ymat h) k a * Gmat z x a i
        = if a = i then Qmat z x k k * Gmat z x a i else 0 := by
      intro a
      rw [e2row a]
      by_cases ha : a = i <;> simp [ha]
    rw [Finset.sum_congr rfl fun a _ => hterm a]
    simp
  have e3row : ∀ a : Fin p,
      (Pmat z x * (Ymat h)ᵀ) k a = if a = k then Pmat z x k i else 0 := by
    intro a
    rw [Matrix.mul_apply]
    have hterm : ∀ b : Fin d, Pmat z x k b * (Ymat (p := p) (d := d) h)ᵀ b a
        = if a = k then (if b = i then Pmat z x k b else 0) else 0 := by
      intro b
      rw [Matrix.transpose_apply, hYs a b]
      by_cases ha : a = k <;> by_cases hb : b = i <;> simp [ha, hb]
    rw [Finset.sum_congr rfl fun b _ => hterm b]
    by_cases ha : a = k <;> simp [ha]
  have e3 : (Pmat z x * (Ymat h)ᵀ * Pmat z x) k i = Pmat z x k i * Pmat z x k i := by
    rw [Matrix.mul_apply]
    have hterm : ∀ a : Fin p, (Pmat z x * (Ymat h)ᵀ) k a * Pmat z x a i
        = if a = k then Pmat z x k i * Pmat z x a i else 0 := by
      intro a
      rw [e3row a]
      by_cases ha : a = k <;> simp [ha]
    rw [Finset.sum_congr rfl fun a _ => hterm a]
    simp
  rw [fderiv_Fentry_eq hz k i x h, e1, e2, e3]

/-! ### The two trace identities -/

private theorem cmapR_smul {m n : Type*} (c : ℝ) (M : Matrix m n ℝ) :
    cmapR (c • M) = (c : ℂ) • cmapR M := by
  ext i j
  simp only [cmapR_apply, Matrix.smul_apply, smul_eq_mul, Complex.ofReal_mul]

private theorem transpose_Ymat_mul_Ymat (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Ymat x)ᵀ * Ymat x = (d : ℂ) • R4C.cmat (gram ((matrixEquivE p d).symm x)) := by
  have hd0 : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  set Y := (matrixEquivE p d).symm x with hY
  have h1 : Yᵀ * Y = (d : ℝ) • gram Y := by
    rw [gram, smul_smul, mul_inv_cancel₀ hd0, one_smul]
  have hcast : (((d : ℕ) : ℝ) : ℂ) = ((d : ℕ) : ℂ) := by push_cast; ring
  calc (Ymat x)ᵀ * Ymat x = cmapR Yᵀ * cmapR Y := by rw [Ymat, cmapR_transpose]
    _ = cmapR (Yᵀ * Y) := (cmapR_mul _ _).symm
    _ = cmapR ((d : ℝ) • gram Y) := by rw [h1]
    _ = (((d : ℕ) : ℝ) : ℂ) • cmapR (gram Y) := cmapR_smul _ _
    _ = (d : ℂ) • R4C.cmat (gram Y) := by rw [hcast, cmat_eq_cmapR]

private theorem cmat_gram_mul_Gmat (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    R4C.cmat (gram ((matrixEquivE p d).symm x)) * Gmat z x = 1 + z • Gmat z x := by
  have h := cmat_sub_mul_resolvC (isHermitian_gram ((matrixEquivE p d).symm x)) (z := z) hz
  rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul] at h
  have h' : R4C.cmat (gram ((matrixEquivE p d).symm x)) * Gmat z x - z • Gmat z x = 1 := h
  rw [sub_eq_iff_eq_add] at h'
  exact h'

/-- `tr(Yᵀ Y G) = d (d + z tr G)`. -/
theorem trace_transpose_Ymat_mul_Pmat (hz : z.im ≠ 0) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ((Ymat x)ᵀ * Pmat z x).trace = (d : ℂ) * ((d : ℂ) + z * (Gmat z x).trace) := by
  calc ((Ymat x)ᵀ * Pmat z x).trace
      = (((Ymat x)ᵀ * Ymat x) * Gmat z x).trace := by rw [Pmat, Matrix.mul_assoc]
    _ = (((d : ℂ) • R4C.cmat (gram ((matrixEquivE p d).symm x))) * Gmat z x).trace := by
        rw [transpose_Ymat_mul_Ymat hd]
    _ = (d : ℂ) * (R4C.cmat (gram ((matrixEquivE p d).symm x)) * Gmat z x).trace := by
        rw [Matrix.smul_mul, Matrix.trace_smul, smul_eq_mul]
    _ = (d : ℂ) * ((1 : Matrix (Fin d) (Fin d) ℂ) + z • Gmat z x).trace := by
        rw [cmat_gram_mul_Gmat hz]
    _ = (d : ℂ) * ((d : ℂ) + z * (Gmat z x).trace) := by
        rw [Matrix.trace_add, Matrix.trace_one, Matrix.trace_smul, smul_eq_mul,
          Fintype.card_fin]

/-- `tr(Y G Yᵀ) = d (d + z tr G)`. -/
theorem trace_Qmat (hz : z.im ≠ 0) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Qmat z x).trace = (d : ℂ) * ((d : ℂ) + z * (Gmat z x).trace) := by
  rw [Qmat, Matrix.trace_mul_comm]
  exact trace_transpose_Ymat_mul_Pmat hz hd x

private theorem sum_mul_eq_trace {m n : Type*} [Fintype m] [Fintype n]
    (A B : Matrix m n ℂ) :
    ∑ k : m, ∑ i : n, A k i * B k i = (Aᵀ * B).trace := by
  rw [Matrix.trace]
  simp only [Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply]
  exact Finset.sum_comm

/-- The left side of the summed Stein identity. -/
theorem sum_Ymat_mul_Pmat (hz : z.im ≠ 0) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ∑ k : Fin p, ∑ i : Fin d, Ymat x k i * Pmat z x k i
      = (d : ℂ) * ((d : ℂ) + z * (Gmat z x).trace) := by
  rw [sum_mul_eq_trace]
  exact trace_transpose_Ymat_mul_Pmat hz hd x

/-- `∑_{k,i} (Y G)_{ki}² = d (tr G + z tr G²)`. -/
theorem sum_Pmat_sq (hz : z.im ≠ 0) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) :
    ∑ k : Fin p, ∑ i : Fin d, Pmat z x k i * Pmat z x k i
      = (d : ℂ) * ((Gmat z x).trace + z * (Gmat z x * Gmat z x).trace) := by
  rw [sum_mul_eq_trace]
  have hPP : (Pmat z x)ᵀ * Pmat z x
      = (d : ℂ) • (Gmat z x + z • (Gmat z x * Gmat z x)) := by
    have h1 : (Pmat z x)ᵀ = Gmat z x * (Ymat x)ᵀ := by
      rw [Pmat, Matrix.transpose_mul, transpose_Gmat hz]
    rw [h1, Pmat]
    calc Gmat z x * (Ymat x)ᵀ * (Ymat x * Gmat z x)
        = Gmat z x * (((Ymat x)ᵀ * Ymat x) * Gmat z x) := by simp only [Matrix.mul_assoc]
      _ = Gmat z x * (((d : ℂ) • R4C.cmat (gram ((matrixEquivE p d).symm x))) * Gmat z x) := by
          rw [transpose_Ymat_mul_Ymat hd]
      _ = (d : ℂ) • (Gmat z x * (R4C.cmat (gram ((matrixEquivE p d).symm x)) * Gmat z x)) := by
          rw [Matrix.smul_mul, Matrix.mul_smul]
      _ = (d : ℂ) • (Gmat z x * ((1 : Matrix (Fin d) (Fin d) ℂ) + z • Gmat z x)) := by
          rw [cmat_gram_mul_Gmat hz]
      _ = (d : ℂ) • (Gmat z x + z • (Gmat z x * Gmat z x)) := by
          rw [Matrix.mul_add, Matrix.mul_one, Matrix.mul_smul]
  rw [hPP, Matrix.trace_smul, Matrix.trace_add, Matrix.trace_smul, smul_eq_mul, smul_eq_mul]

/-! ### The divergence -/

private theorem trace_eq_sum_diag {n : ℕ} (A : Matrix (Fin n) (Fin n) ℂ) :
    A.trace = ∑ i : Fin n, A i i := by
  rw [Matrix.trace]
  exact Finset.sum_congr rfl fun i _ => rfl

/-- **The divergence identity.** -/
theorem sum_fderiv_Fentry_single (hz : z.im ≠ 0) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ∑ k : Fin p, ∑ i : Fin d,
        fderiv ℝ (Fentry (p := p) (d := d) z k i) x
          (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
      = (p : ℂ) * (Gmat z x).trace
        - ((d : ℂ) + z * (Gmat z x).trace) * (Gmat z x).trace
        - (Gmat z x).trace - z * (Gmat z x * Gmat z x).trace := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hG : (Gmat z x).trace = ∑ i : Fin d, Gmat z x i i := trace_eq_sum_diag _
  have hQ : (Qmat z x).trace = ∑ k : Fin p, Qmat z x k k := trace_eq_sum_diag _
  have expand : ∀ (k : Fin p) (i : Fin d),
      fderiv ℝ (Fentry (p := p) (d := d) z k i) x
          (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
        = Gmat z x i i - (d : ℂ)⁻¹ * Qmat z x k k * Gmat z x i i
          - (d : ℂ)⁻¹ * (Pmat z x k i * Pmat z x k i) := by
    intro k i
    rw [fderiv_Fentry_single hz k i x]
    ring
  rw [Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ => expand k i]
  have split : ∑ k : Fin p, ∑ i : Fin d,
      (Gmat z x i i - (d : ℂ)⁻¹ * Qmat z x k k * Gmat z x i i
        - (d : ℂ)⁻¹ * (Pmat z x k i * Pmat z x k i))
      = (∑ k : Fin p, ∑ i : Fin d, Gmat z x i i)
        - (∑ k : Fin p, ∑ i : Fin d, (d : ℂ)⁻¹ * Qmat z x k k * Gmat z x i i)
        - ∑ k : Fin p, ∑ i : Fin d, (d : ℂ)⁻¹ * (Pmat z x k i * Pmat z x k i) := by
    rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  have s1 : ∑ _k : Fin p, ∑ i : Fin d, Gmat z x i i = (p : ℂ) * (Gmat z x).trace := by
    rw [hG, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have s2 : ∑ k : Fin p, ∑ i : Fin d, (d : ℂ)⁻¹ * Qmat z x k k * Gmat z x i i
      = (d : ℂ)⁻¹ * (Qmat z x).trace * (Gmat z x).trace := by
    have inner : ∀ k : Fin p, ∑ i : Fin d, (d : ℂ)⁻¹ * Qmat z x k k * Gmat z x i i
        = ((d : ℂ)⁻¹ * Qmat z x k k) * (Gmat z x).trace := by
      intro k
      rw [hG, Finset.mul_sum]
    rw [Finset.sum_congr rfl fun k _ => inner k, ← Finset.sum_mul, ← Finset.mul_sum, ← hQ]
  have s3 : ∑ k : Fin p, ∑ i : Fin d, (d : ℂ)⁻¹ * (Pmat z x k i * Pmat z x k i)
      = (d : ℂ)⁻¹ * ∑ k : Fin p, ∑ i : Fin d, Pmat z x k i * Pmat z x k i := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.mul_sum]
  rw [split, s1, s2, s3, trace_Qmat hz hd x, sum_Pmat_sq hz hd x]
  field_simp
  ring

/-! ### The two traces as `stieltjesC` and `stieltjes2C` -/

theorem sfun_matrixEquivE (z : ℂ) (Z : Matrix (Fin p) (Fin d) ℝ) :
    sfun z (matrixEquivE p d Z) = R4C.stieltjesC (gram Z) z := by
  rw [sfun_eq, MeasurableEquiv.symm_apply_apply]

/-- `d⁻¹ tr G²`, on the flattened space. -/
noncomputable def s2fun (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) : ℂ :=
  (d : ℂ)⁻¹ * (Gmat z x * Gmat z x).trace

theorem s2fun_eq (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    s2fun z x = R4C.stieltjes2C (gram ((matrixEquivE p d).symm x)) z := rfl

theorem s2fun_matrixEquivE (z : ℂ) (Z : Matrix (Fin p) (Fin d) ℝ) :
    s2fun z (matrixEquivE p d Z) = R4C.stieltjes2C (gram Z) z := by
  rw [s2fun_eq, MeasurableEquiv.symm_apply_apply]

theorem contDiff_s2fun (hz : z.im ≠ 0) {n : WithTop ℕ∞} :
    ContDiff ℝ n (s2fun (p := p) (d := d) z) := by
  have heq : s2fun (p := p) (d := d) z
      = fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (d : ℂ)⁻¹ * ∑ i : Fin d, ∑ j : Fin d, Gmat z x i j * Gmat z x j i := by
    funext x
    rw [s2fun, trace_eq_sum_diag]
    congr 1
  rw [heq]
  exact contDiff_const.mul (ContDiff.sum fun i _ => ContDiff.sum fun j _ =>
    (contDiff_Gmat_entry hz i j).mul (contDiff_Gmat_entry hz j i))

theorem trace_Gmat_eq (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    (Gmat z (matrixEquivE p d Z)).trace = (d : ℂ) * R4C.stieltjesC (gram Z) z := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  rw [R4C.stieltjesC, Gmat, MeasurableEquiv.symm_apply_apply]
  field_simp

theorem trace_Gmat_sq_eq (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    (Gmat z (matrixEquivE p d Z) * Gmat z (matrixEquivE p d Z)).trace
      = (d : ℂ) * R4C.stieltjes2C (gram Z) z := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  rw [R4C.stieltjes2C, Gmat, MeasurableEquiv.symm_apply_apply]
  field_simp

/-! ### Measurability, bounds and integrability under `gaussianMatrix` -/

theorem measurable_sN (hz : z.im ≠ 0) :
    Measurable (fun Z : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (gram Z) z) := by
  have h : Continuous (sfun (p := p) (d := d) z) := (contDiff_sfun (n := 0) hz).continuous
  have h2 := h.measurable.comp (matrixEquivE p d).measurable
  simpa [Function.comp_def, sfun_matrixEquivE] using h2

theorem measurable_s2N (hz : z.im ≠ 0) :
    Measurable (fun Z : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjes2C (gram Z) z) := by
  have h : Continuous (s2fun (p := p) (d := d) z) := (contDiff_s2fun (n := 0) hz).continuous
  have h2 := h.measurable.comp (matrixEquivE p d).measurable
  simpa [Function.comp_def, s2fun_matrixEquivE] using h2

theorem norm_sN_le (hz : 0 < z.im) (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    ‖R4C.stieltjesC (gram Z) z‖ ≤ 1 / z.im :=
  R4C.norm_stieltjesC_le (isHermitian_gram Z) hz hd

theorem norm_s2N_le (hz : 0 < z.im) (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    ‖R4C.stieltjes2C (gram Z) z‖ ≤ 1 / z.im ^ 2 :=
  R4C.norm_stieltjes2C_le (isHermitian_gram Z) hz hd

theorem integrable_sN (hz : 0 < z.im) (hd : 0 < d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (gram Z) z)
      (gaussianMatrix p d) :=
  Integrable.mono' (integrable_const (1 / z.im)) (measurable_sN hz.ne').aestronglyMeasurable
    (Filter.Eventually.of_forall fun Z => norm_sN_le hz hd Z)

theorem integrable_sN_sq (hz : 0 < z.im) (hd : 0 < d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ => (R4C.stieltjesC (gram Z) z) ^ 2)
      (gaussianMatrix p d) := by
  refine Integrable.mono' (integrable_const ((1 / z.im) ^ 2))
    ((measurable_sN hz.ne').pow_const 2).aestronglyMeasurable
    (Filter.Eventually.of_forall fun Z => ?_)
  rw [norm_pow]
  exact pow_le_pow_left₀ (norm_nonneg _) (norm_sN_le hz hd Z) 2

theorem integrable_s2N (hz : 0 < z.im) (hd : 0 < d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjes2C (gram Z) z)
      (gaussianMatrix p d) :=
  Integrable.mono' (integrable_const (1 / z.im ^ 2)) (measurable_s2N hz.ne').aestronglyMeasurable
    (Filter.Eventually.of_forall fun Z => norm_s2N_le hz hd Z)

/-! ### The summed Stein identity -/

private theorem integrable_coord (r : Fin (p * d)) :
    Integrable (fun x : EuclideanSpace ℝ (Fin (p * d)) => (x r : ℝ))
      (GaussianMeasure.stdGaussianE (p * d)) := by
  rw [stdGaussianE_eq]
  have h := IsGaussian.integrable_inner_mul_of_norm_fderiv_le
    (μ := stdGaussian (EuclideanSpace ℝ (Fin (p * d))))
    (F := fun _ : EuclideanSpace ℝ (Fin (p * d)) => (1 : ℝ)) (L := 0)
    contDiff_const (fun x => by simp) (EuclideanSpace.single r (1 : ℝ))
  simp only [mul_one] at h
  have hinner : ∀ x : EuclideanSpace ℝ (Fin (p * d)),
      inner ℝ (EuclideanSpace.single r (1 : ℝ) : EuclideanSpace ℝ (Fin (p * d))) x = x r := by
    intro x
    rw [EuclideanSpace.inner_single_left]
    simp
  simpa [hinner] using h

private theorem integrable_entry_matrix (k : Fin p) (i : Fin d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ => ((Z k i : ℝ) : ℂ)) (gaussianMatrix p d) := by
  have h1 : Integrable
      (fun x : EuclideanSpace ℝ (Fin (p * d)) => ((x (finProdFinEquiv (k, i)) : ℝ) : ℂ))
      (GaussianMeasure.stdGaussianE (p * d)) := (integrable_coord _).ofReal
  have h2 := (measurePreserving_matrixEquivE p d).integrable_comp_of_integrable h1
  simpa [Function.comp_def] using h2

private theorem integrable_lhs_entry (hz : 0 < z.im) (hd : 0 < d) (k : Fin p) (i : Fin d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ =>
      ((Z k i : ℝ) : ℂ) * Fentry z k i (matrixEquivE p d Z)) (gaussianMatrix p d) := by
  have hcont : Continuous (Fentry (p := p) (d := d) z k i) :=
    (contDiff_Fentry (n := 0) hz.ne' k i).continuous
  have hmeas : Measurable
      (fun Z : Matrix (Fin p) (Fin d) ℝ => Fentry z k i (matrixEquivE p d Z)) :=
    hcont.measurable.comp (matrixEquivE p d).measurable
  exact (integrable_entry_matrix k i).mul_bdd (c := BP z d) hmeas.aestronglyMeasurable
    (Filter.Eventually.of_forall fun Z => norm_Pmat_le hz hd (matrixEquivE p d Z) k i)

private theorem integrable_rhs_entry (hz : 0 < z.im) (hd : 0 < d) (k : Fin p) (i : Fin d) :
    Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ =>
      fderiv ℝ (Fentry (p := p) (d := d) z k i) (matrixEquivE p d Z)
        (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))) (gaussianMatrix p d) := by
  have hcont : Continuous (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
      fderiv ℝ (Fentry (p := p) (d := d) z k i) x
        (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))) :=
    ((contDiff_Fentry (n := 1) hz.ne' k i).continuous_fderiv one_ne_zero).clm_apply
      continuous_const
  have hmeas := hcont.measurable.comp (matrixEquivE p d).measurable
  refine Integrable.mono'
    (integrable_const (Lconst z p d * ‖EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ)‖))
    (by simpa [Function.comp_def] using hmeas.aestronglyMeasurable)
    (Filter.Eventually.of_forall fun Z => ?_)
  refine ((fderiv ℝ (Fentry (p := p) (d := d) z k i)
    (matrixEquivE p d Z)).le_opNorm _).trans ?_
  exact mul_le_mul_of_nonneg_right (norm_fderiv_Fentry_le hz hd k i _) (norm_nonneg _)

theorem lhs_pointwise (hz : z.im ≠ 0) (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    ∑ k : Fin p, ∑ i : Fin d, ((Z k i : ℝ) : ℂ) * Fentry z k i (matrixEquivE p d Z)
      = (d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * R4C.stieltjesC (gram Z) z := by
  have hconv : ∀ (k : Fin p) (i : Fin d),
      ((Z k i : ℝ) : ℂ) * Fentry z k i (matrixEquivE p d Z)
        = Ymat (matrixEquivE p d Z) k i * Pmat z (matrixEquivE p d Z) k i := by
    intro k i
    rw [Ymat_apply, Fentry, MeasurableEquiv.symm_apply_apply]
  rw [Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ => hconv k i,
    sum_Ymat_mul_Pmat hz hd, trace_Gmat_eq hd]
  ring

theorem rhs_pointwise (hz : z.im ≠ 0) (hd : 0 < d) (Z : Matrix (Fin p) (Fin d) ℝ) :
    ∑ k : Fin p, ∑ i : Fin d,
        fderiv ℝ (Fentry (p := p) (d := d) z k i) (matrixEquivE p d Z)
          (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))
      = ((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * R4C.stieltjesC (gram Z) z
        - (z * (d : ℂ) ^ 2) * (R4C.stieltjesC (gram Z) z) ^ 2
        - (z * (d : ℂ)) * R4C.stieltjes2C (gram Z) z := by
  rw [sum_fderiv_Fentry_single hz hd, trace_Gmat_eq hd, trace_Gmat_sq_eq hd]
  ring

/-- **The summed Stein identity.** -/
theorem stein_summed (hz : 0 < z.im) (hd : 0 < d) :
    ∫ Z, ((d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * R4C.stieltjesC (gram Z) z) ∂gaussianMatrix p d
      = ∫ Z, (((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * R4C.stieltjesC (gram Z) z
          - (z * (d : ℂ) ^ 2) * (R4C.stieltjesC (gram Z) z) ^ 2
          - (z * (d : ℂ)) * R4C.stieltjes2C (gram Z) z) ∂gaussianMatrix p d := by
  have hL := integrable_lhs_entry hz hd (p := p) (d := d) (z := z)
  have hR := integrable_rhs_entry hz hd (p := p) (d := d) (z := z)
  have e1 : ∫ Z, ((d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * R4C.stieltjesC (gram Z) z)
        ∂gaussianMatrix p d
      = ∫ Z, (∑ k : Fin p, ∑ i : Fin d,
          ((Z k i : ℝ) : ℂ) * Fentry z k i (matrixEquivE p d Z)) ∂gaussianMatrix p d :=
    integral_congr_ae (Filter.Eventually.of_forall fun Z => (lhs_pointwise hz.ne' hd Z).symm)
  have e2 : ∫ Z, (((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * R4C.stieltjesC (gram Z) z
          - (z * (d : ℂ) ^ 2) * (R4C.stieltjesC (gram Z) z) ^ 2
          - (z * (d : ℂ)) * R4C.stieltjes2C (gram Z) z) ∂gaussianMatrix p d
      = ∫ Z, (∑ k : Fin p, ∑ i : Fin d,
          fderiv ℝ (Fentry (p := p) (d := d) z k i) (matrixEquivE p d Z)
            (EuclideanSpace.single (finProdFinEquiv (k, i)) (1 : ℝ))) ∂gaussianMatrix p d :=
    integral_congr_ae (Filter.Eventually.of_forall fun Z => (rhs_pointwise hz.ne' hd Z).symm)
  rw [e1, e2,
    integral_finsetSum _ (fun k _ => integrable_finsetSum _ fun i _ => hL k i),
    integral_finsetSum _ (fun k _ => integrable_finsetSum _ fun i _ => hR k i)]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [integral_finsetSum _ (fun i _ => hL k i), integral_finsetSum _ (fun i _ => hR k i)]
  refine Finset.sum_congr rfl fun i _ => ?_
  exact integral_entry_mul_gaussianMatrix_complex (contDiff_Fentry (n := 1) hz.ne' k i)
    (norm_fderiv_Fentry_le hz hd k i) k i

/-! ### The quadratic equation for the mean, with its residual -/

/-- **Step 2 of `notes/archive/rmt_R1.md`, as an equality.** -/
theorem quad_integral_eq (hz : 0 < z.im) (hd : 0 < d) :
    MP.quad ((p : ℝ) / d) z (∫ Z, R4C.stieltjesC (gram Z) z ∂gaussianMatrix p d)
      = -(z * (∫ Z, (R4C.stieltjesC (gram Z) z
              - ∫ Z', R4C.stieltjesC (gram Z') z ∂gaussianMatrix p d) ^ 2
            ∂gaussianMatrix p d)
          + (d : ℂ)⁻¹ * (∫ Z, R4C.stieltjesC (gram Z) z ∂gaussianMatrix p d)
          + z * (d : ℂ)⁻¹ * (∫ Z, R4C.stieltjes2C (gram Z) z ∂gaussianMatrix p d)) := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hisN := integrable_sN (p := p) (d := d) hz hd
  have hisN2 := integrable_sN_sq (p := p) (d := d) hz hd
  have hisS2 := integrable_s2N (p := p) (d := d) hz hd
  set A : ℂ := ∫ Z, R4C.stieltjesC (gram Z) z ∂gaussianMatrix p d with hAdef
  set A2 : ℂ := ∫ Z, R4C.stieltjes2C (gram Z) z ∂gaussianMatrix p d with hA2def
  set M2 : ℂ := ∫ Z, (R4C.stieltjesC (gram Z) z) ^ 2 ∂gaussianMatrix p d with hM2def
  set V : ℂ := ∫ Z, (R4C.stieltjesC (gram Z) z - A) ^ 2 ∂gaussianMatrix p d with hVdef
  have hLint : ∫ Z, ((d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * R4C.stieltjesC (gram Z) z)
      ∂gaussianMatrix p d = (d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * A := by
    rw [integral_add (integrable_const _) (hisN.const_mul _), integral_const, integral_const_mul]
    simp only [hAdef]
    simp
  have hi1 : Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ =>
      ((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * R4C.stieltjesC (gram Z) z
        - (z * (d : ℂ) ^ 2) * (R4C.stieltjesC (gram Z) z) ^ 2) (gaussianMatrix p d) :=
    (hisN.const_mul _).sub (hisN2.const_mul _)
  have hRint : ∫ Z, (((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * R4C.stieltjesC (gram Z) z
        - (z * (d : ℂ) ^ 2) * (R4C.stieltjesC (gram Z) z) ^ 2
        - (z * (d : ℂ)) * R4C.stieltjes2C (gram Z) z) ∂gaussianMatrix p d
      = ((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * A - (z * (d : ℂ) ^ 2) * M2
        - (z * (d : ℂ)) * A2 := by
    rw [integral_sub hi1 (hisS2.const_mul _),
      integral_sub (hisN.const_mul _) (hisN2.const_mul _),
      integral_const_mul, integral_const_mul, integral_const_mul]
  have heq : (d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * A
      = ((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * A - (z * (d : ℂ) ^ 2) * M2
        - (z * (d : ℂ)) * A2 := by
    rw [← hLint, ← hRint]
    exact stein_summed hz hd
  have hM : M2 = A ^ 2 + V := by
    have hexp : ∀ Z : Matrix (Fin p) (Fin d) ℝ,
        (R4C.stieltjesC (gram Z) z - A) ^ 2
          = (R4C.stieltjesC (gram Z) z) ^ 2 - (2 * A) * R4C.stieltjesC (gram Z) z + A ^ 2 :=
      fun Z => by ring
    have hi2 : Integrable (fun Z : Matrix (Fin p) (Fin d) ℝ =>
        (R4C.stieltjesC (gram Z) z) ^ 2 - (2 * A) * R4C.stieltjesC (gram Z) z)
        (gaussianMatrix p d) := hisN2.sub (hisN.const_mul _)
    have hV : V = M2 - (2 * A) * A + A ^ 2 := by
      rw [hVdef, integral_congr_ae (Filter.Eventually.of_forall hexp),
        integral_add hi2 (integrable_const _),
        integral_sub hisN2 (hisN.const_mul _), integral_const_mul, integral_const]
      simp only [hAdef, hM2def]
      simp
    rw [hV]
    ring
  have heq2 : 1 + z * A
      = ((p : ℂ) / (d : ℂ) - 1 - (d : ℂ)⁻¹) * A - z * M2 - z * (d : ℂ)⁻¹ * A2 := by
    refine mul_left_cancel₀ (pow_ne_zero 2 hdC) ?_
    calc (d : ℂ) ^ 2 * (1 + z * A) = (d : ℂ) ^ 2 + (z * (d : ℂ) ^ 2) * A := by ring
      _ = ((p : ℂ) * (d : ℂ) - (d : ℂ) ^ 2 - (d : ℂ)) * A - (z * (d : ℂ) ^ 2) * M2
            - (z * (d : ℂ)) * A2 := heq
      _ = (d : ℂ) ^ 2 * (((p : ℂ) / (d : ℂ) - 1 - (d : ℂ)⁻¹) * A - z * M2
            - z * (d : ℂ)⁻¹ * A2) := by field_simp
  rw [hM] at heq2
  have hcast : (((p : ℝ) / (d : ℝ) : ℝ) : ℂ) = (p : ℂ) / (d : ℂ) := by push_cast; ring
  change z * A ^ 2 + (z + 1 - (((p : ℝ) / (d : ℝ) : ℝ) : ℂ)) * A + 1 = _
  rw [hcast]
  linear_combination heq2

/-! ### The residual bound -/

/-- **Step 2, in the residual form `RMT/R1.lean` states.** The two variance hypotheses are
what `R1.variance_le_of_lipschitz` supplies for `Re s_N` and `Im s_N`. -/
theorem norm_quad_integral_le_aux (hz : 0 < z.im) (hd : 0 < d) {V : ℝ}
    (hre : ∫ Z, ((sfun z (matrixEquivE p d Z)).re
        - ∫ Z', (sfun z (matrixEquivE p d Z')).re ∂gaussianMatrix p d) ^ 2
        ∂gaussianMatrix p d ≤ V)
    (him : ∫ Z, ((sfun z (matrixEquivE p d Z)).im
        - ∫ Z', (sfun z (matrixEquivE p d Z')).im ∂gaussianMatrix p d) ^ 2
        ∂gaussianMatrix p d ≤ V) :
    ‖MP.quad ((p : ℝ) / d) z (∫ Z, R4C.stieltjesC (gram Z) z ∂gaussianMatrix p d)‖
      ≤ ‖z‖ * (2 * V) + 1 / ((d : ℝ) * z.im) + ‖z‖ / ((d : ℝ) * z.im ^ 2) := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hisN := integrable_sN (p := p) (d := d) hz hd
  have hisS2 := integrable_s2N (p := p) (d := d) hz hd
  have hmre : Measurable fun Z : Matrix (Fin p) (Fin d) ℝ => (R4C.stieltjesC (gram Z) z).re :=
    Complex.measurable_re.comp (measurable_sN hz.ne')
  have hmim : Measurable fun Z : Matrix (Fin p) (Fin d) ℝ => (R4C.stieltjesC (gram Z) z).im :=
    Complex.measurable_im.comp (measurable_sN hz.ne')
  set A : ℂ := ∫ Z, R4C.stieltjesC (gram Z) z ∂gaussianMatrix p d with hAdef
  set A2 : ℂ := ∫ Z, R4C.stieltjes2C (gram Z) z ∂gaussianMatrix p d with hA2def
  set V2 : ℂ := ∫ Z, (R4C.stieltjesC (gram Z) z - A) ^ 2 ∂gaussianMatrix p d with hV2def
  -- the two real integrals are the parts of `A`
  have hAre : ∫ Z, (R4C.stieltjesC (gram Z) z).re ∂gaussianMatrix p d = A.re := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.reCLM hisN
    simpa [hAdef] using h
  have hAim : ∫ Z, (R4C.stieltjesC (gram Z) z).im ∂gaussianMatrix p d = A.im := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM hisN
    simpa [hAdef] using h
  have hre' : ∫ Z, ((R4C.stieltjesC (gram Z) z).re - A.re) ^ 2 ∂gaussianMatrix p d ≤ V := by
    have h := hre
    simp only [sfun_matrixEquivE] at h
    rwa [hAre] at h
  have him' : ∫ Z, ((R4C.stieltjesC (gram Z) z).im - A.im) ^ 2 ∂gaussianMatrix p d ≤ V := by
    have h := him
    simp only [sfun_matrixEquivE] at h
    rwa [hAim] at h
  -- the complex second moment
  have hbre : ∀ Z : Matrix (Fin p) (Fin d) ℝ, |(R4C.stieltjesC (gram Z) z).re - A.re|
      ≤ 1 / z.im + ‖A‖ := by
    intro Z
    have h1 : |(R4C.stieltjesC (gram Z) z).re| ≤ 1 / z.im :=
      (Complex.abs_re_le_norm _).trans (norm_sN_le hz hd Z)
    have h2 : |A.re| ≤ ‖A‖ := Complex.abs_re_le_norm A
    calc |(R4C.stieltjesC (gram Z) z).re - A.re|
        ≤ |(R4C.stieltjesC (gram Z) z).re| + |A.re| := abs_sub _ _
      _ ≤ 1 / z.im + ‖A‖ := add_le_add h1 h2
  have hbim : ∀ Z : Matrix (Fin p) (Fin d) ℝ, |(R4C.stieltjesC (gram Z) z).im - A.im|
      ≤ 1 / z.im + ‖A‖ := by
    intro Z
    have h1 : |(R4C.stieltjesC (gram Z) z).im| ≤ 1 / z.im :=
      (Complex.abs_im_le_norm _).trans (norm_sN_le hz hd Z)
    have h2 : |A.im| ≤ ‖A‖ := Complex.abs_im_le_norm A
    calc |(R4C.stieltjesC (gram Z) z).im - A.im|
        ≤ |(R4C.stieltjesC (gram Z) z).im| + |A.im| := abs_sub _ _
      _ ≤ 1 / z.im + ‖A‖ := add_le_add h1 h2
  have hIre2 : Integrable
      (fun Z : Matrix (Fin p) (Fin d) ℝ => ((R4C.stieltjesC (gram Z) z).re - A.re) ^ 2)
      (gaussianMatrix p d) := by
    refine Integrable.mono' (integrable_const ((1 / z.im + ‖A‖) ^ 2))
      (((hmre.sub measurable_const).pow_const 2)).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Z => ?_)
    rw [Real.norm_eq_abs, abs_pow]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbre Z) 2
  have hIim2 : Integrable
      (fun Z : Matrix (Fin p) (Fin d) ℝ => ((R4C.stieltjesC (gram Z) z).im - A.im) ^ 2)
      (gaussianMatrix p d) := by
    refine Integrable.mono' (integrable_const ((1 / z.im + ‖A‖) ^ 2))
      (((hmim.sub measurable_const).pow_const 2)).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Z => ?_)
    rw [Real.norm_eq_abs, abs_pow]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbim Z) 2
  have hV2 : ‖V2‖ ≤ 2 * V := by
    have hstep : ‖V2‖ ≤ ∫ Z, ‖(R4C.stieltjesC (gram Z) z - A) ^ 2‖ ∂gaussianMatrix p d :=
      norm_integral_le_integral_norm _
    have hpt : ∀ Z : Matrix (Fin p) (Fin d) ℝ, ‖(R4C.stieltjesC (gram Z) z - A) ^ 2‖
        = ((R4C.stieltjesC (gram Z) z).re - A.re) ^ 2
          + ((R4C.stieltjesC (gram Z) z).im - A.im) ^ 2 := by
      intro Z
      rw [norm_pow, Complex.sq_norm, Complex.normSq_apply, Complex.sub_re, Complex.sub_im]
      ring
    rw [integral_congr_ae (Filter.Eventually.of_forall hpt),
      integral_add hIre2 hIim2] at hstep
    linarith
  have hAn : ‖A‖ ≤ 1 / z.im := by
    have h := norm_integral_le_of_norm_le_const (μ := gaussianMatrix p d) (C := 1 / z.im)
      (Filter.Eventually.of_forall fun Z => norm_sN_le hz hd Z)
    simpa [hAdef] using h
  have hA2n : ‖A2‖ ≤ 1 / z.im ^ 2 := by
    have h := norm_integral_le_of_norm_le_const (μ := gaussianMatrix p d) (C := 1 / z.im ^ 2)
      (Filter.Eventually.of_forall fun Z => norm_s2N_le hz hd Z)
    simpa [hA2def] using h
  have hdinv : ‖(d : ℂ)⁻¹‖ = 1 / (d : ℝ) := by
    rw [norm_inv, Complex.norm_natCast, one_div]
  rw [quad_integral_eq hz hd, norm_neg]
  have hsplit : ‖z * V2 + (d : ℂ)⁻¹ * A + z * (d : ℂ)⁻¹ * A2‖
      ≤ ‖z‖ * ‖V2‖ + ‖(d : ℂ)⁻¹‖ * ‖A‖ + ‖z‖ * ‖(d : ℂ)⁻¹‖ * ‖A2‖ := by
    refine (norm_add_le _ _).trans ?_
    have h1 := norm_add_le (z * V2) ((d : ℂ)⁻¹ * A)
    have h2 : ‖z * V2‖ = ‖z‖ * ‖V2‖ := norm_mul _ _
    have h3 : ‖(d : ℂ)⁻¹ * A‖ = ‖(d : ℂ)⁻¹‖ * ‖A‖ := norm_mul _ _
    have h4 : ‖z * (d : ℂ)⁻¹ * A2‖ = ‖z‖ * ‖(d : ℂ)⁻¹‖ * ‖A2‖ := by
      rw [norm_mul, norm_mul]
    linarith
  have hzn : (0 : ℝ) ≤ ‖z‖ := norm_nonneg z
  have hb1 : ‖z‖ * ‖V2‖ ≤ ‖z‖ * (2 * V) := mul_le_mul_of_nonneg_left hV2 hzn
  have hb2 : ‖(d : ℂ)⁻¹‖ * ‖A‖ ≤ (1 / (d : ℝ)) * (1 / z.im) := by
    rw [hdinv]
    exact mul_le_mul_of_nonneg_left hAn (by positivity)
  have hb3 : ‖z‖ * ‖(d : ℂ)⁻¹‖ * ‖A2‖ ≤ ‖z‖ * (1 / (d : ℝ)) * (1 / z.im ^ 2) := by
    rw [hdinv]
    exact mul_le_mul_of_nonneg_left hA2n (by positivity)
  have he1 : (1 / (d : ℝ)) * (1 / z.im) = 1 / ((d : ℝ) * z.im) := by
    field_simp
  have he2 : ‖z‖ * (1 / (d : ℝ)) * (1 / z.im ^ 2) = ‖z‖ / ((d : ℝ) * z.im ^ 2) := by
    field_simp
  linarith

end SteinStep

end StackedSVD
