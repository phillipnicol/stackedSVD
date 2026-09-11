/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Trace
import StackedSVD.Prob.Tensorization
import StackedSVD.Prob.Chebyshev

/-!
# Item G4: the isotropic Marchenko-Pastur law at four moments

`notes/archive/prop_single_table_general.md` section 5, unit G4. Twin of `RMT/R2.lean` for a general
noise law `ν` (mean 0, variance 1, finite fourth moment). The conclusion is a measure bound at
finite `p` and `d`, not a limit: for every pair of directions `x`, `y` inside the unit ball,

`P (ε ≤ ‖xᵀ G y - (x ⬝ y) d⁻¹ tr G‖) ≤ isoRate ν z p d / ε²`,

with `isoRate ν z p d = isoConst ν z * (1 + p/d) ^ 8 / d`. The constant `isoConst` mentions
only `∫ x ^ 4 ∂ν`, `‖z‖` and `z.im`; it is free of `x`, `y`, `p`, `d` and `ε`. That uniformity
is the content of the unit: unit G7 applies the bound to a sequence of directions that changes
with `N`, so a limit at one fixed pair would be useless.

This file imports none of `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` or
`Vendor/COLT83/` (choice 8 of the note).

## Route

Write `A = gram Y`, `G = resolvC A z`, `g k = rowVec Y k`, `A_k = gram (updateRow Y k 0)`,
`G_k = resolvC A_k z`, `s = stieltjesC A z`, `m = MP.mC (p/d) z`.

1. `z Q_u + (u ⬝ u) = ∑ k, (1 - s k) (u ⬝ g k) (g kᵀ G_k u)`, an exact identity for every
   realization, from `A G = 1 + z G` and the rank-one downdate of unit G2.
2. The row toolkit: the fourth moment of a complex linear form of one row, and the
   leave-one-out increment `Q_u - Q_u^{(k)} = -(1 + α k)⁻¹ (g kᵀ G_k u)²`.
3. The `L²` trace law `∫ ‖s - m‖² ≤ C / d`, from the residual bound of unit G3.
4. The variance of `Q_u - (u ⬝ u) s`, by Efron-Stein on the rows (`Prob/Tensorization.lean`).
5. The mean, by taking the expectation of step 1 and using the self-consistent equation. This
   step is what a bare variance bound cannot give: the law of `noiseMatrix ν p d` is not
   invariant under a sign flip of a column, so `E G` is not a multiple of the identity.
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace GenRMT

variable {p d : ℕ} {z : ℂ}

/-! ### Small helpers on the complex embedding

F37 (2026-09-09): this section used to hold `cvec_dotProduct_cvec`, `vecMulVec_self_mulVec'`
and `dotProduct_vecMulVec_mulVec'` (`private`). `Companion.lean:404`, `:73` and `:81` now
carry the public canonical copy of each (`SpikedModel.cvec_dotProduct_cvec`,
`R4C.vecMulVec_self_mulVec`, `R4C.dotProduct_vecMulVec_mulVec`); the rest of this file uses
those directly. -/


/-! ### Milestone 1: the exact identity `(*)` -/

/-- **The Gram matrix is the sum of the rank-one row terms.** -/
theorem gram_eq_sum_vecMulVec (Y : Matrix (Fin p) (Fin d) ℝ) :
    gram Y = ∑ k, Matrix.vecMulVec (rowVec Y k) (rowVec Y k) := by
  have hdd : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = ((d : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d)]
  ext a b
  rw [Matrix.sum_apply]
  simp only [gram, Matrix.smul_apply, smul_eq_mul, Matrix.mul_apply, Matrix.transpose_apply,
    Matrix.vecMulVec_apply, rowVec_apply]
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by rw [← hdd]; ring

/-- The same over `ℂ`. -/
private theorem cmat_gram_eq_sum (Y : Matrix (Fin p) (Fin d) ℝ) :
    R4C.cmat (gram Y)
      = ∑ k, Matrix.vecMulVec (R4C.cvec (rowVec Y k)) (R4C.cvec (rowVec Y k)) := by
  ext a b
  rw [Matrix.sum_apply]
  simp only [R4C.cmat, Matrix.map_apply, gram_eq_sum_vecMulVec, Matrix.sum_apply,
    Matrix.vecMulVec_apply, R4C.cvec, Complex.ofReal_sum, Complex.ofReal_mul]

/-- **The resolvent identity** `W G = 1 + z G`. -/
theorem cmat_mul_resolvC {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) :
    R4C.cmat W * R4C.resolvC W z = 1 + z • R4C.resolvC W z := by
  have hdet : IsUnit (R4C.cmat W - z • (1 : Matrix (Fin D) (Fin D) ℂ)).det :=
    ResolvDeriv.isUnit_det_cmat_sub hW hz.ne'
  have h : (R4C.cmat W - z • (1 : Matrix (Fin D) (Fin D) ℂ)) * R4C.resolvC W z = 1 :=
    Matrix.mul_nonsing_inv _ hdet
  rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul, sub_eq_iff_eq_add] at h
  exact h

/-- **The identity `(*)`, first line.** `A G = 1 + z G` read at the pair `(u, u)`, with
`A = ∑ k, g k g kᵀ` expanded on the left. -/
theorem z_mul_qformC_add_eq (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (u : Fin d → ℝ) :
    z * R4C.qformC (gram Y) z u + ((u ⬝ᵥ u : ℝ) : ℂ)
      = ∑ k, ((u ⬝ᵥ rowVec Y k : ℝ) : ℂ) * R4C.cformC (gram Y) z (rowVec Y k) u := by
  set G := R4C.resolvC (gram Y) z with hG
  have hleft : R4C.cvec u ⬝ᵥ ((R4C.cmat (gram Y) * G) *ᵥ R4C.cvec u)
      = ((u ⬝ᵥ u : ℝ) : ℂ) + z * R4C.qformC (gram Y) z u := by
    rw [hG, cmat_mul_resolvC (gram_isHermitian Y) hz, Matrix.add_mulVec, Matrix.one_mulVec,
      Matrix.smul_mulVec, dotProduct_add, dotProduct_smul, smul_eq_mul,
      SpikedModel.cvec_dotProduct_cvec]
    rfl
  have hright : R4C.cvec u ⬝ᵥ ((R4C.cmat (gram Y) * G) *ᵥ R4C.cvec u)
      = ∑ k, ((u ⬝ᵥ rowVec Y k : ℝ) : ℂ) * R4C.cformC (gram Y) z (rowVec Y k) u := by
    rw [← Matrix.mulVec_mulVec, cmat_gram_eq_sum, Matrix.sum_mulVec, dotProduct_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [R4C.dotProduct_vecMulVec_mulVec, SpikedModel.cvec_dotProduct_cvec]
    rfl
  rw [add_comm, ← hleft]
  exact hright

/-- **The downdate at the pair `(g k, u)`.** `g kᵀ G u = (1 - s k) g kᵀ G_k u`, so the sum of
`(*)` reads at the leave-one-out resolvent, which is independent of row `k`. -/
theorem cformC_rowVec_eq (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (u : Fin d → ℝ)
    (k : Fin p) :
    R4C.cformC (gram Y) z (rowVec Y k) u
      = (1 - sRow Y z k) * R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u := by
  have hne : 1 - sRow Y z k ≠ 0 := one_sub_sRow_ne_zero Y hz k
  have hkey : R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u
      = R4C.cformC (gram Y) z (rowVec Y k) u
        + sRow Y z k * R4C.cformC (gram Y) z (rowVec Y k) u / (1 - sRow Y z k) := by
    rw [gram_updateRow_zero]
    exact R4C.cformC_sub_vecMulVec (gram_isHermitian Y) hz (rowVec Y k) (rowVec Y k) u
  rw [hkey]
  field_simp
  ring

/-- **The identity `(*)`.** For every realization and every direction `u`,
`z Q_u + (u ⬝ u) = ∑ k, (1 - s k) (u ⬝ g k) (g kᵀ G_k u)`. Both the variance of milestone 4
and the mean of milestone 5 start here. -/
theorem z_mul_qformC_add_eq_loo (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (u : Fin d → ℝ) :
    z * R4C.qformC (gram Y) z u + ((u ⬝ᵥ u : ℝ) : ℂ)
      = ∑ k, (1 - sRow Y z k) * (((u ⬝ᵥ rowVec Y k : ℝ) : ℂ)
          * R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u) := by
  rw [z_mul_qformC_add_eq Y hz u]
  exact Finset.sum_congr rfl fun k _ => by rw [cformC_rowVec_eq Y hz u k]; ring

/-! ### Milestone 2: the row toolkit

The two facts the fourth-moment step needs: the sharp bound `‖G u‖² ≤ (u ⬝ u)/η²` on the
image vector, and the leave-one-out increment of the quadratic form. -/

/-- A matrix with real entries and orthogonal columns keeps the sum of the squared norms.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem sum_normSq_mulVec_real_orth {D : ℕ} (U : Matrix (Fin D) (Fin D) ℂ)
    (hU : Uᵀ * U = 1) (hUr : ∀ a b, (starRingEnd ℂ) (U a b) = U a b) (w : Fin D → ℂ) :
    ∑ a, ‖(U *ᵥ w) a‖ ^ 2 = ∑ a, ‖w a‖ ^ 2 := by
  have hc : ∀ x : ℂ, ((‖x‖ ^ 2 : ℝ) : ℂ) = x * (starRingEnd ℂ) x := fun x => by
    rw [Complex.sq_norm, Complex.mul_conj]
  refine Complex.ofReal_inj.mp ?_
  rw [Complex.ofReal_sum, Complex.ofReal_sum]
  have hL : ∀ a : Fin D, ((‖(U *ᵥ w) a‖ ^ 2 : ℝ) : ℂ)
      = ∑ b, ∑ c, (U a b * U a c) * (w b * (starRingEnd ℂ) (w c)) := by
    intro a
    have h1 : (U *ᵥ w) a = ∑ b, U a b * w b := rfl
    rw [hc, h1, map_sum]
    simp only [map_mul, hUr]
    rw [Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun c _ => by ring
  have hswap : ∑ a, ∑ b, ∑ c, (U a b * U a c) * (w b * (starRingEnd ℂ) (w c))
      = ∑ b, ∑ c, (∑ a, U a b * U a c) * (w b * (starRingEnd ℂ) (w c)) := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun c _ => (Finset.sum_mul _ _ _).symm
  have hone : ∀ b c : Fin D, (∑ a, U a b * U a c) = if b = c then (1 : ℂ) else 0 := by
    intro b c
    have h := congrArg (fun M : Matrix (Fin D) (Fin D) ℂ => M b c) hU
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.one_apply] at h
    exact h
  rw [Finset.sum_congr rfl fun a _ => hL a, hswap]
  rw [Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun c _ => by rw [hone b c]]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [Finset.sum_eq_single b (fun c _ hc' => by rw [if_neg (Ne.symm hc'), zero_mul])
    (fun hb => absurd (Finset.mem_univ b) hb), if_pos rfl, one_mul, hc]

/-- **The sharp bound on the image vector.** `∑ a, ‖(G u) a‖² ≤ (u ⬝ u) / η²`. The entrywise
bound `‖G a b‖ ≤ 1/η` would lose a factor `D` here. -/
theorem sum_normSq_resolvC_mulVec_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (u : Fin D → ℝ) :
    ∑ a, ‖(R4C.resolvC W z *ᵥ R4C.cvec u) a‖ ^ 2 ≤ (u ⬝ᵥ u) / z.im ^ 2 := by
  have hz' : z.im ≠ 0 := hz.ne'
  set V : Matrix (Fin D) (Fin D) ℂ := R4C.cmat (R4.eigU hW) with hV
  have hVr : ∀ a b, (starRingEnd ℂ) (V a b) = V a b := by
    intro a b; simp [hV, R4C.cmat]
  have hVtr : ∀ a b, (starRingEnd ℂ) (Vᵀ a b) = Vᵀ a b := by
    intro a b; simp [hV, R4C.cmat]
  set f : Fin D → ℂ := fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹ with hf
  set q : Fin D → ℂ := Vᵀ *ᵥ R4C.cvec u with hq
  have hsplit : R4C.resolvC W z *ᵥ R4C.cvec u = V *ᵥ (Matrix.diagonal f *ᵥ q) := by
    rw [R4C.resolvC_eq_conj hW hz', ← hV, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]
  have hq2 : ∑ a, ‖q a‖ ^ 2 = u ⬝ᵥ u := by
    have hVt : (Vᵀ)ᵀ * Vᵀ = 1 := by
      rw [Matrix.transpose_transpose, hV]; exact R4C.ceigU_mul_transpose hW
    rw [hq, sum_normSq_mulVec_real_orth Vᵀ hVt hVtr]
    simp only [R4C.cvec, Complex.norm_real, Real.norm_eq_abs, sq_abs, dotProduct]
    exact Finset.sum_congr rfl fun a _ => by ring
  rw [hsplit, sum_normSq_mulVec_real_orth V (R4C.transpose_ceigU_mul hW) hVr]
  have hstep : ∀ a : Fin D, ‖(Matrix.diagonal f *ᵥ q) a‖ ^ 2 ≤ (z.im)⁻¹ ^ 2 * ‖q a‖ ^ 2 := by
    intro a
    rw [Matrix.mulVec_diagonal, norm_mul, mul_pow]
    have h1 : ‖f a‖ ≤ (z.im)⁻¹ := R4C.norm_inv_eigenvalue_sub_le hW hz a
    exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (norm_nonneg _) h1 2) (by positivity)
  calc ∑ a, ‖(Matrix.diagonal f *ᵥ q) a‖ ^ 2 ≤ ∑ a, (z.im)⁻¹ ^ 2 * ‖q a‖ ^ 2 :=
        Finset.sum_le_sum fun a _ => hstep a
    _ = (z.im)⁻¹ ^ 2 * (u ⬝ᵥ u) := by rw [← Finset.mul_sum, hq2]
    _ = (u ⬝ᵥ u) / z.im ^ 2 := by rw [inv_pow]; ring

/-- `1 - s k = -z (Gc) k k`: the scalar of the downdate is a diagonal entry of the companion
resolvent, scaled by `z`. -/
theorem one_sub_sRow_eq (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d) (k : Fin p) :
    1 - sRow Y z k = -(z * R4C.resolvC (gramC Y) z k k) := by
  have h1 := one_add_alphaRow_mul Y hz k
  have h2 := mul_one_add_alphaRow Y hz hd k
  have hne : (1 + alphaRow Y z k) ≠ 0 := by
    intro h; rw [h, zero_mul] at h1; exact zero_ne_one h1
  have hz0 : (1 + alphaRow Y z k)
      * ((1 - sRow Y z k) + z * R4C.resolvC (gramC Y) z k k) = 0 := by
    linear_combination h1 + h2
  have h3 := (mul_eq_zero.mp hz0).resolve_left hne
  linear_combination h3

/-- **The uniform bound on the downdate scalar.** `‖1 - s k‖ ≤ ‖z‖/η`, for every row and every
realization. The naive bound through `g k ⬝ g k` is useless: that norm is random and
unbounded. -/
theorem norm_one_sub_sRow_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (k : Fin p) : ‖1 - sRow Y z k‖ ≤ ‖z‖ / z.im := by
  rw [one_sub_sRow_eq Y hz hd k, norm_neg, norm_mul]
  have h := norm_resolvC_diag_le (gramC_isHermitian Y) hz k
  calc ‖z‖ * ‖R4C.resolvC (gramC Y) z k k‖ ≤ ‖z‖ * (1 / z.im) :=
        mul_le_mul_of_nonneg_left h (norm_nonneg z)
    _ = ‖z‖ / z.im := by ring

/-- The downdate at the pair `(u, g k)`, the mirror of `cformC_rowVec_eq`. -/
theorem cformC_rowVec_right_eq (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (u : Fin d → ℝ)
    (k : Fin p) :
    R4C.cformC (gram Y) z u (rowVec Y k)
      = (1 - sRow Y z k) * R4C.cformC (gram (Matrix.updateRow Y k 0)) z u (rowVec Y k) := by
  have hne : 1 - sRow Y z k ≠ 0 := one_sub_sRow_ne_zero Y hz k
  have hkey : R4C.cformC (gram (Matrix.updateRow Y k 0)) z u (rowVec Y k)
      = R4C.cformC (gram Y) z u (rowVec Y k)
        + R4C.cformC (gram Y) z u (rowVec Y k) * sRow Y z k / (1 - sRow Y z k) := by
    rw [gram_updateRow_zero]
    exact R4C.cformC_sub_vecMulVec (gram_isHermitian Y) hz (rowVec Y k) u (rowVec Y k)
  rw [hkey]
  field_simp
  ring

/-- **The leave-one-out increment of the quadratic form.** `Q_u - Q_u^{(k)}` is a product of
two forms of the leave-one-out resolvent, and both are linear in row `k`. -/
theorem qformC_sub_qformC_loo (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (u : Fin d → ℝ)
    (k : Fin p) :
    R4C.qformC (gram Y) z u - R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
      = -((1 - sRow Y z k)
          * (R4C.cformC (gram (Matrix.updateRow Y k 0)) z u (rowVec Y k)
             * R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u)) := by
  have hne : 1 - sRow Y z k ≠ 0 := one_sub_sRow_ne_zero Y hz k
  have hkey : R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
      = R4C.qformC (gram Y) z u
        + R4C.cformC (gram Y) z u (rowVec Y k) * R4C.cformC (gram Y) z (rowVec Y k) u
          / (1 - sRow Y z k) := by
    rw [gram_updateRow_zero]
    exact R4C.cformC_sub_vecMulVec (gram_isHermitian Y) hz (rowVec Y k) u u
  rw [hkey, cformC_rowVec_right_eq Y hz u k, cformC_rowVec_eq Y hz u k]
  field_simp
  ring

/-! ### The fourth moment of a complex linear form of one row -/

/-- `∑ l, x l ^ 4 ≤ (∑ l, x l ^ 2) ^ 2`. -/
private theorem sum_pow_four_le {D : ℕ} (x : Fin D → ℝ) :
    ∑ l, x l ^ 4 ≤ (∑ l, x l ^ 2) ^ 2 := by
  have h := Finset.sum_sq_le_sq_sum_of_nonneg (s := (Finset.univ : Finset (Fin D)))
    (f := fun l => x l ^ 2) (fun l _ => sq_nonneg _)
  calc ∑ l, x l ^ 4 = ∑ l, (x l ^ 2) ^ 2 := Finset.sum_congr rfl fun l _ => by ring
    _ ≤ (∑ l, x l ^ 2) ^ 2 := h

/-- **The fourth moment of a complex linear form.** For `t` with i.i.d. coordinates of law `ν`
the form `∑ l, a l t l` has fourth moment at most `2 (ν₄ + 3) (∑ l, ‖a l‖²)²`. Milestone 2
uses it twice, once at the coefficient vector `G_k u` and once at `u`. -/
theorem integral_clinForm_pow_four_le {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ} (a : Fin D → ℂ) :
    Integrable (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 4)
        (Measure.pi fun _ : Fin D => ν) ∧
    ∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) * (∑ l, ‖a l‖ ^ 2) ^ 2 := by
  have hprob := hν.prob
  set ar : Fin D → ℝ := fun l => (a l).re with har
  set ai : Fin D → ℝ := fun l => (a l).im with hai
  have hRe : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).re = LinForm.linForm ar t := by
    intro t
    rw [Complex.re_sum]
    exact Finset.sum_congr rfl fun l _ => by
      simp [Complex.mul_re, har]
  have hIm : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).im = LinForm.linForm ai t := by
    intro t
    rw [Complex.im_sum]
    exact Finset.sum_congr rfl fun l _ => by
      simp [Complex.mul_im, hai]
  have hEq : ∀ t : Fin D → ℝ, ‖∑ l, a l * (t l : ℂ)‖ ^ 4
      = (LinForm.linForm ar t ^ 2 + LinForm.linForm ai t ^ 2) ^ 2 := by
    intro t
    have h2 : ‖∑ l, a l * (t l : ℂ)‖ ^ 2
        = LinForm.linForm ar t ^ 2 + LinForm.linForm ai t ^ 2 := by
      rw [Complex.sq_norm, Complex.normSq_apply, hRe t, hIm t]; ring
    calc ‖∑ l, a l * (t l : ℂ)‖ ^ 4 = (‖∑ l, a l * (t l : ℂ)‖ ^ 2) ^ 2 := by ring
      _ = _ := by rw [h2]
  have h4r := LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 ar
  have h4i := LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 ai
  have hdom : Integrable (fun t : Fin D → ℝ =>
      2 * (LinForm.linForm ar t ^ 4 + LinForm.linForm ai t ^ 4))
      (Measure.pi fun _ : Fin D => ν) := (h4r.add h4i).const_mul 2
  have hptw : ∀ t : Fin D → ℝ, ‖∑ l, a l * (t l : ℂ)‖ ^ 4
      ≤ 2 * (LinForm.linForm ar t ^ 4 + LinForm.linForm ai t ^ 4) := by
    intro t
    rw [hEq t]
    nlinarith [sq_nonneg (LinForm.linForm ar t ^ 2 - LinForm.linForm ai t ^ 2),
      sq_nonneg (LinForm.linForm ar t), sq_nonneg (LinForm.linForm ai t)]
  have hfun : (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 4)
      = fun t => (LinForm.linForm ar t ^ 2 + LinForm.linForm ai t ^ 2) ^ 2 := funext hEq
  have hmeas : AEStronglyMeasurable (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 4)
      (Measure.pi fun _ : Fin D => ν) := by
    rw [hfun]
    exact Measurable.aestronglyMeasurable
      ((((LinForm.measurable_linForm ar).pow_const 2).add
        ((LinForm.measurable_linForm ai).pow_const 2)).pow_const 2)
  have hint : Integrable (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 4)
      (Measure.pi fun _ : Fin D => ν) := by
    refine hdom.mono' hmeas (Filter.Eventually.of_forall fun t => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact hptw t
  refine ⟨hint, ?_⟩
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  have hbr := LinForm.integral_linForm_pow_four_le hν.mean hν.var hν.mom4 ar
  have hbi := LinForm.integral_linForm_pow_four_le hν.mean hν.var hν.mom4 ai
  have hsr := sum_pow_four_le ar
  have hsi := sum_pow_four_le ai
  have hnorm : ∑ l, ‖a l‖ ^ 2 = (∑ l, ar l ^ 2) + ∑ l, ai l ^ 2 := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun l _ => by
      rw [Complex.sq_norm, Complex.normSq_apply, har, hai]; ring
  have hnn1 : (0 : ℝ) ≤ ∑ l, ar l ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  have hnn2 : (0 : ℝ) ≤ ∑ l, ai l ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  have hnu4 : (0 : ℝ) ≤ ∫ x, x ^ 4 ∂ν := integral_pow_four_nonneg ν
  calc ∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ∫ t, 2 * (LinForm.linForm ar t ^ 4 + LinForm.linForm ai t ^ 4)
          ∂(Measure.pi fun _ : Fin D => ν) :=
        integral_mono hint hdom hptw
    _ = 2 * ((∫ t, LinForm.linForm ar t ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
          + ∫ t, LinForm.linForm ai t ^ 4 ∂(Measure.pi fun _ : Fin D => ν)) := by
        rw [integral_const_mul, integral_add h4r h4i]
    _ ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) * (∑ l, ‖a l‖ ^ 2) ^ 2 := by
        have e1 : (∫ x, x ^ 4 ∂ν) * ∑ l, ar l ^ 4
            ≤ (∫ x, x ^ 4 ∂ν) * (∑ l, ar l ^ 2) ^ 2 := mul_le_mul_of_nonneg_left hsr hnu4
        have e2 : (∫ x, x ^ 4 ∂ν) * ∑ l, ai l ^ 4
            ≤ (∫ x, x ^ 4 ∂ν) * (∑ l, ai l ^ 2) ^ 2 := mul_le_mul_of_nonneg_left hsi hnu4
        have e4 : (0 : ℝ)
            ≤ ((∫ x, x ^ 4 ∂ν) + 3) * ((∑ l, ar l ^ 2) * ∑ l, ai l ^ 2) :=
          mul_nonneg hnu (mul_nonneg hnn1 hnn2)
        rw [hnorm]
        nlinarith [hbr, hbi, e1, e2, e4]


/-! ### Two conversions between `lintegral` and `integral` -/

/-- A nonnegative measurable function with a finite `lintegral` is integrable.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integrable_of_lintegral_ofReal_ne_top {α : Type*} [MeasurableSpace α]
    {μ : Measure α} {f : α → ℝ} (hf : Measurable f) (hnn : ∀ a, 0 ≤ f a)
    (h : ∫⁻ a, ENNReal.ofReal (f a) ∂μ ≠ ⊤) : Integrable f μ := by
  refine ⟨hf.aestronglyMeasurable, ?_⟩
  rw [hasFiniteIntegral_iff_ofReal (Filter.Eventually.of_forall hnn)]
  exact lt_of_le_of_ne le_top h

/-- A `lintegral` bound on a nonnegative integrable function is an `integral` bound.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_le_of_lintegral_le {α : Type*} [MeasurableSpace α]
    {μ : Measure α} {f : α → ℝ} (hf : Integrable f μ) (hnn : ∀ a, 0 ≤ f a) {C : ℝ}
    (hC : 0 ≤ C) (h : ∫⁻ a, ENNReal.ofReal (f a) ∂μ ≤ ENNReal.ofReal C) : ∫ a, f a ∂μ ≤ C := by
  rw [← ofReal_integral_eq_lintegral_ofReal hf (Filter.Eventually.of_forall hnn)] at h
  exact (ENNReal.ofReal_le_ofReal_iff hC).mp h


/-! ### Milestone 3: the `L²` trace law -/

/-- The constant of the `L²` trace law: the stability constant of the MP quadratic squared,
times the residual constant of unit G3, which is `O(1/D)`. -/
noncomputable def traceSqBound (ν : Measure ℝ) (z : ℂ) (P D : ℕ) : ℝ :=
  ((P : ℝ) / D * (‖z‖ / z.im)) ^ 2
      * ((1 / (‖z‖ * (z.im * rootLb ((P : ℝ) / D) z ^ 2))) ^ 2 + 16 / z.im ^ 2)
    * residConst ν D z

/-- **The `L²` trace law.** `∫ ‖s - m‖² ≤ traceSqBound`, which is `O(1/D)` at a fixed ratio.
The deterministic stability of the quadratic (`norm_stieltjesC_sub_mC_le`) holds on the event
where the residual is small; on the complement the trivial bound `‖s - m‖ ≤ 2/η` and the
Chebyshev bound `measure_resid_ge_le` give the same order. -/
theorem lintegral_normSq_stieltjesC_sub_mC_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hP : 0 < P) (hD : 0 < D) :
    ∫⁻ Y, ENNReal.ofReal (‖R4C.stieltjesC (gram Y) z - MP.mC ((P : ℝ) / D) z‖ ^ 2)
        ∂(noiseMatrix ν P D)
      ≤ ENNReal.ofReal (traceSqBound ν z P D) := by
  obtain ⟨n, rfl⟩ : ∃ n, P = n + 1 := ⟨P - 1, (Nat.succ_pred_eq_of_pos hP).symm⟩
  have hprob := hν.prob
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hPR : (0 : ℝ) < ((n + 1 : ℕ) : ℝ) := by exact_mod_cast hP
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  set c : ℝ := ((n + 1 : ℕ) : ℝ) / D with hc
  have hcpos : 0 < c := by rw [hc]; positivity
  set kap : ℝ := c * (‖z‖ / z.im) with hkap
  have hkpos : 0 < kap := by rw [hkap]; positivity
  set th : ℝ := 1 / (‖z‖ * (z.im * rootLb c z ^ 2)) with hth
  have hδ : 0 < rootLb c z := rootLb_pos hz0
  have hthpos : 0 < th := by rw [hth]; positivity
  set e0 : ℝ := 1 / (2 * kap) with he0
  have he0pos : 0 < e0 := by rw [he0]; positivity
  set S : Set (Matrix (Fin (n + 1)) (Fin D) ℝ) := {Y | e0 ≤ resid Y z} with hS
  have hSmeas : MeasurableSet S := measurableSet_le measurable_const (measurable_resid z)
  set B : ℝ := 4 / z.im ^ 2 with hB
  have hBnn : 0 ≤ B := by rw [hB]; positivity
  have hRnn : 0 ≤ residConst ν D z := residConst_nonneg ν D hz
  have hpt : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
      ENNReal.ofReal (‖R4C.stieltjesC (gram Y) z - MP.mC c z‖ ^ 2)
        ≤ ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2)
          + S.indicator (fun _ => ENNReal.ofReal B) Y := by
    intro Y
    by_cases hY : Y ∈ S
    · refine le_trans ?_ le_add_self
      rw [Set.indicator_of_mem hY]
      refine ENNReal.ofReal_le_ofReal ?_
      have h1 : ‖R4C.stieltjesC (gram Y) z‖ ≤ 1 / z.im :=
        R4C.norm_stieltjesC_le (gram_isHermitian Y) hz hD
      have h2 : ‖MP.mC c z‖ ≤ (z.im)⁻¹ := MP.norm_mC_le hcpos.le hz
      have hhalf : (2 : ℝ) / z.im = 1 / z.im + 1 / z.im := by ring
      have h3 : ‖R4C.stieltjesC (gram Y) z - MP.mC c z‖ ≤ 2 / z.im := by
        refine (norm_sub_le _ _).trans ?_
        rw [inv_eq_one_div] at h2
        rw [hhalf]
        linarith
      have hsq : (2 / z.im) ^ 2 = 4 / z.im ^ 2 := by rw [div_pow]; norm_num
      rw [hB, ← hsq]
      exact pow_le_pow_left₀ (norm_nonneg _) h3 2
    · refine le_trans ?_ le_self_add
      refine ENNReal.ofReal_le_ofReal ?_
      have hlt : resid Y z < e0 := not_le.mp hY
      have hval : kap * e0 = 1 / 2 := by rw [he0]; field_simp
      have hklt : kap * resid Y z < 1 / 2 := by
        rw [← hval]; exact mul_lt_mul_of_pos_left hlt hkpos
      have hsmall : c * (‖z‖ / z.im) * resid Y z + |c - c| * (1 / z.im) ≤ 1 / 2 := by
        rw [sub_self, abs_zero, zero_mul, add_zero, ← hkap]
        linarith
      have hbnd := norm_stieltjesC_sub_mC_le (c := c) Y hcpos hz hP hD (by rw [← hc]; exact hsmall)
      rw [← hc, sub_self, abs_zero, zero_mul, add_zero, ← hkap] at hbnd
      have hb2 : ‖R4C.stieltjesC (gram Y) z - MP.mC c z‖ ≤ th * (kap * resid Y z) := by
        refine hbnd.trans (le_of_eq ?_)
        rw [hth]; ring
      calc ‖R4C.stieltjesC (gram Y) z - MP.mC c z‖ ^ 2
          ≤ (th * (kap * resid Y z)) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hb2 2
        _ = th ^ 2 * kap ^ 2 * resid Y z ^ 2 := by ring
  have hmr : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      ENNReal.ofReal (resid Y z ^ 2) :=
    ENNReal.measurable_ofReal.comp ((measurable_resid z).pow_const 2)
  have hmeasf : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2) :=
    ENNReal.measurable_ofReal.comp (measurable_const.mul ((measurable_resid z).pow_const 2))
  have hfirst : ∫⁻ Y, ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2)
      ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (th ^ 2 * kap ^ 2 * residConst ν D z) := by
    have hcongr : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
        ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2)
          = ENNReal.ofReal (th ^ 2 * kap ^ 2) * ENNReal.ofReal (resid Y z ^ 2) := fun Y =>
      ENNReal.ofReal_mul (by positivity)
    rw [lintegral_congr hcongr, lintegral_const_mul _ hmr,
      ENNReal.ofReal_mul (by positivity : (0 : ℝ) ≤ th ^ 2 * kap ^ 2)]
    gcongr
    exact lintegral_sq_resid_le hν hz hD
  calc ∫⁻ Y, ENNReal.ofReal (‖R4C.stieltjesC (gram Y) z - MP.mC c z‖ ^ 2)
        ∂(noiseMatrix ν (n + 1) D)
      ≤ ∫⁻ Y, (ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2)
          + S.indicator (fun _ => ENNReal.ofReal B) Y) ∂(noiseMatrix ν (n + 1) D) :=
        lintegral_mono hpt
    _ = (∫⁻ Y, ENNReal.ofReal (th ^ 2 * kap ^ 2 * resid Y z ^ 2) ∂(noiseMatrix ν (n + 1) D))
          + ENNReal.ofReal B * noiseMatrix ν (n + 1) D S := by
        rw [lintegral_add_left hmeasf, lintegral_indicator hSmeas, setLIntegral_const]
    _ ≤ ENNReal.ofReal (th ^ 2 * kap ^ 2 * residConst ν D z)
          + ENNReal.ofReal B * ENNReal.ofReal (residConst ν D z / e0 ^ 2) := by
        gcongr
        exact measure_resid_ge_le hν hz hP hD he0pos
    _ = ENNReal.ofReal (traceSqBound ν z (n + 1) D) := by
        rw [← ENNReal.ofReal_mul hBnn,
          ← ENNReal.ofReal_add (by positivity) (by positivity)]
        congr 1
        rw [traceSqBound, ← hc, ← hkap, ← hth, hB, he0]
        field_simp
        ring

/-! ### Milestone 4: the variance -/

/-- The bilinear form of the resolvent is symmetric, because the resolvent of a real symmetric
matrix is symmetric. -/
theorem cformC_comm {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian) (hz : z.im ≠ 0)
    (x y : Fin D → ℝ) : R4C.cformC W z x y = R4C.cformC W z y x := by
  change R4C.cvec x ⬝ᵥ (R4C.resolvC W z *ᵥ R4C.cvec y)
    = R4C.cvec y ⬝ᵥ (R4C.resolvC W z *ᵥ R4C.cvec x)
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, ResolvDeriv.transpose_resolvC hW hz,
    dotProduct_comm]

/-- The form at a scaled vector is a complex linear form of that vector.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem cformC_smul_left_eq {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (r : ℝ)
    (t u : Fin D → ℝ) :
    R4C.cformC W z (r • t) u
      = ∑ l, ((r : ℂ) * (R4C.resolvC W z *ᵥ R4C.cvec u) l) * (t l : ℂ) := by
  simp only [R4C.cformC, dotProduct, R4C.cvec, Pi.smul_apply, smul_eq_mul, Complex.ofReal_mul]
  exact Finset.sum_congr rfl fun l _ => by ring

/-- **The fourth moment of the row form.** For a Hermitian `W` free of the row, the form
`g kᵀ G u` has fourth moment `O(1/D²)` in the row. -/
theorem integral_normSq_cformC_row_pow_four_le {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (hz : 0 < z.im) (hD : 0 < D) {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (u : Fin D → ℝ) :
    Integrable (fun t : Fin D → ℝ => ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4)
        (Measure.pi fun _ : Fin D => ν) ∧
    ∫ t, ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) * ((u ⬝ᵥ u) / ((D : ℝ) * z.im ^ 2)) ^ 2 := by
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  set a : Fin D → ℂ := fun l => (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC W z *ᵥ R4C.cvec u) l
    with ha
  have hfun : (fun t : Fin D → ℝ => ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4)
      = fun t => ‖∑ l, a l * (t l : ℂ)‖ ^ 4 := by
    funext t
    rw [cformC_smul_left_eq]
  have hsq : ((Real.sqrt D)⁻¹ : ℝ) ^ 2 = ((D : ℝ))⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hDR.le]
  have hcoef : ∑ l, ‖a l‖ ^ 2 ≤ (u ⬝ᵥ u) / ((D : ℝ) * z.im ^ 2) := by
    have hstep : ∀ l, ‖a l‖ ^ 2
        = ((D : ℝ))⁻¹ * ‖(R4C.resolvC W z *ᵥ R4C.cvec u) l‖ ^ 2 := by
      intro l
      rw [ha]
      simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs, hsq]
    rw [Finset.sum_congr rfl fun l _ => hstep l, ← Finset.mul_sum]
    have h := sum_normSq_resolvC_mulVec_le hW hz u
    calc ((D : ℝ))⁻¹ * ∑ l, ‖(R4C.resolvC W z *ᵥ R4C.cvec u) l‖ ^ 2
        ≤ ((D : ℝ))⁻¹ * ((u ⬝ᵥ u) / z.im ^ 2) := by
          exact mul_le_mul_of_nonneg_left h (by positivity)
      _ = (u ⬝ᵥ u) / ((D : ℝ) * z.im ^ 2) := by field_simp
  obtain ⟨hint, hbnd⟩ := integral_clinForm_pow_four_le hν a
  rw [hfun]
  refine ⟨hint, hbnd.trans ?_⟩
  have hnu : (0 : ℝ) ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) := by
    have := integral_pow_four_nonneg ν; linarith
  have hnn : (0 : ℝ) ≤ ∑ l, ‖a l‖ ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hnn hcoef 2) hnu

/-- **The centered isotropic form** `Q_u - (u ⬝ u) s`, the object of statement 1. -/
noncomputable def isoF (z : ℂ) (u : Fin d → ℝ) (Y : Matrix (Fin p) (Fin d) ℝ) : ℂ :=
  R4C.qformC (gram Y) z u - ((u ⬝ᵥ u : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z

/-- `isoF` is bounded by `2 (u ⬝ u)/η` for every realization. -/
theorem norm_isoF_le (hz : 0 < z.im) (hd : 0 < d) (u : Fin d → ℝ)
    (Y : Matrix (Fin p) (Fin d) ℝ) : ‖isoF z u Y‖ ≤ 2 * (u ⬝ᵥ u) / z.im := by
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have h1 : ‖R4C.qformC (gram Y) z u‖ ≤ (u ⬝ᵥ u) / z.im :=
    R4C.norm_qformC_le (gram_isHermitian Y) hz u
  have h2 : ‖((u ⬝ᵥ u : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z‖ ≤ (u ⬝ᵥ u) / z.im := by
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg huu]
    have hst := R4C.norm_stieltjesC_le (gram_isHermitian Y) hz hd
    calc (u ⬝ᵥ u) * ‖R4C.stieltjesC (gram Y) z‖ ≤ (u ⬝ᵥ u) * (1 / z.im) :=
          mul_le_mul_of_nonneg_left hst huu
      _ = (u ⬝ᵥ u) / z.im := by ring
  refine (norm_sub_le _ _).trans ?_
  have hsplit : (2 : ℝ) * (u ⬝ᵥ u) / z.im = (u ⬝ᵥ u) / z.im + (u ⬝ᵥ u) / z.im := by ring
  rw [hsplit]
  linarith

/-- `isoF` is measurable in the matrix. -/
theorem measurable_isoF {P D : ℕ} (z : ℂ) (u : Fin D → ℝ) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => isoF z u Y := by
  simp only [isoF]
  refine Measurable.sub ?_ (measurable_const.mul (measurable_stieltjesC measurable_gram_self z))
  exact measurable_qformC (y := fun _ : Matrix (Fin P) (Fin D) ℝ => u) measurable_gram_self z
    (fun i => measurable_const)

/-- The entries of `Matrix.updateRow x k t` are measurable in `t`.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem measurable_updateRow_entry' {P D : ℕ} (x : Matrix (Fin P) (Fin D) ℝ)
    (k l : Fin P) (i : Fin D) :
    Measurable fun t : Fin D → ℝ => (Matrix.updateRow x k t) l i := by
  by_cases h : l = k
  · subst h
    have hs : ∀ t : Fin D → ℝ, (Matrix.updateRow x l t) l i = t i := fun t => by
      rw [Matrix.updateRow_self]
    simp only [hs]
    exact measurable_pi_apply i
  · have hs : ∀ t : Fin D → ℝ, (Matrix.updateRow x k t) l i = x l i := fun t => by
      rw [Matrix.updateRow_ne h]
    simp only [hs]
    exact measurable_const

/-- Updating the same row twice keeps only the last value.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem updateRow_idem {P D : ℕ} (x : Matrix (Fin P) (Fin D) ℝ) (k : Fin P)
    (t v : Fin D → ℝ) :
    Matrix.updateRow (Matrix.updateRow x k t) k v = Matrix.updateRow x k v := by
  ext l i
  by_cases h : l = k
  · subst h; rw [Matrix.updateRow_self, Matrix.updateRow_self]
  · rw [Matrix.updateRow_ne h, Matrix.updateRow_ne h, Matrix.updateRow_ne h]

/-- **The leave-one-out deviation of `isoF`**, bounded by the square of a linear form of the
row plus a deterministic `1/(D η)`. -/
private theorem norm_isoF_sub_le {P D : ℕ} (hz : 0 < z.im) (hD : 0 < D) (u : Fin D → ℝ)
    (x : Matrix (Fin P) (Fin D) ℝ) (k : Fin P) (t : Fin D → ℝ) :
    ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖
      ≤ ‖z‖ / z.im
          * ‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖ ^ 2
        + (u ⬝ᵥ u) * (1 / ((D : ℝ) * z.im)) := by
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  set Y : Matrix (Fin P) (Fin D) ℝ := Matrix.updateRow x k t with hY
  have hupd : Matrix.updateRow Y k 0 = Matrix.updateRow x k 0 := by
    rw [hY]; exact updateRow_idem x k t 0
  have hrow : rowVec Y k = (Real.sqrt D)⁻¹ • t := by
    have hYk : Y k = t := by rw [hY, Matrix.updateRow_self]
    rw [rowVec, hYk]
  have hq := qformC_sub_qformC_loo Y hz u k
  rw [hupd] at hq
  have hcomm : R4C.cformC (gram (Matrix.updateRow x k 0)) z u (rowVec Y k)
      = R4C.cformC (gram (Matrix.updateRow x k 0)) z (rowVec Y k) u :=
    cformC_comm (gram_isHermitian _) hz.ne' u (rowVec Y k)
  rw [hcomm, hrow] at hq
  have hs := norm_stieltjesC_sub_le Y hz hD k
  rw [hupd] at hs
  have hsplit : isoF z u Y - isoF z u (Matrix.updateRow x k 0)
      = (R4C.qformC (gram Y) z u - R4C.qformC (gram (Matrix.updateRow x k 0)) z u)
        - ((u ⬝ᵥ u : ℝ) : ℂ) * (R4C.stieltjesC (gram Y) z
            - R4C.stieltjesC (gram (Matrix.updateRow x k 0)) z) := by
    simp only [isoF]; ring
  rw [hsplit, hq]
  refine (norm_sub_le _ _).trans (add_le_add ?_ ?_)
  · rw [norm_neg, norm_mul, norm_mul]
    have h1 := norm_one_sub_sRow_le Y hz hD k
    calc ‖1 - sRow Y z k‖
          * (‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖
            * ‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖)
        ≤ (‖z‖ / z.im)
          * (‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖
            * ‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖) :=
          mul_le_mul_of_nonneg_right h1 (by positivity)
      _ = ‖z‖ / z.im
          * ‖R4C.cformC (gram (Matrix.updateRow x k 0)) z ((Real.sqrt D)⁻¹ • t) u‖ ^ 2 := by
          ring
  · rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg huu, norm_sub_rev]
    exact mul_le_mul_of_nonneg_left hs huu

/-- The per-row variance bound, `O(1/D²)`. -/
noncomputable def rowVarBound (ν : Measure ℝ) (z : ℂ) (D : ℕ) : ℝ :=
  4 * (‖z‖ / z.im) ^ 2 * ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) * z.im ^ 2) ^ 2
    + 2 / ((D : ℝ) * z.im) ^ 2

theorem rowVarBound_nonneg (ν : Measure ℝ) {z : ℂ} (D : ℕ) : 0 ≤ rowVarBound ν z D := by
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  rw [rowVarBound]
  have h1 : (0 : ℝ) ≤ 4 * (‖z‖ / z.im) ^ 2 * ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) * z.im ^ 2) ^ 2 :=
    div_nonneg (by positivity) (by positivity)
  have h2 : (0 : ℝ) ≤ 2 / ((D : ℝ) * z.im) ^ 2 := by positivity
  linarith

/-- **The second moment of the leave-one-out deviation in one row.** -/
private theorem integral_normSq_isoF_sub_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin P)
    (x : Matrix (Fin P) (Fin D) ℝ) :
    ∫ t, ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2
        ∂(Measure.pi fun _ : Fin D => ν)
      ≤ rowVarBound ν z D := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  set W : Matrix (Fin D) (Fin D) ℝ := gram (Matrix.updateRow x k 0) with hW
  obtain ⟨hint4, hb4⟩ := integral_normSq_cformC_row_pow_four_le hν hz hD
    (W := W) (gram_isHermitian _) u
  set A : ℝ := ‖z‖ / z.im with hA
  set Bc : ℝ := (u ⬝ᵥ u) * (1 / ((D : ℝ) * z.im)) with hBc
  have hBnn : 0 ≤ Bc := by rw [hBc]; positivity
  set g : (Fin D → ℝ) → ℝ := fun t =>
    2 * (A ^ 2 * ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4) + 2 * Bc ^ 2 with hg
  have hgint : Integrable g (Measure.pi fun _ : Fin D => ν) :=
    ((hint4.const_mul (A ^ 2)).const_mul 2).add (integrable_const _)
  have hptw : ∀ t : Fin D → ℝ,
      ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2 ≤ g t := by
    intro t
    have h := norm_isoF_sub_le hz hD u x k t
    rw [← hW, ← hA, ← hBc] at h
    have hnn : (0 : ℝ) ≤ ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ :=
      norm_nonneg _
    have hsq := pow_le_pow_left₀ hnn h 2
    refine hsq.trans (le_of_eq_of_le rfl ?_)
    rw [hg]
    nlinarith [sq_nonneg (A * ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 2 - Bc)]
  refine (integral_mono_of_nonneg (Filter.Eventually.of_forall fun t => by positivity) hgint
    (Filter.Eventually.of_forall hptw)).trans ?_
  rw [hg, integral_add ((hint4.const_mul (A ^ 2)).const_mul 2) (integrable_const _),
    integral_const_mul, integral_const_mul, integral_const]
  have hone : (Measure.pi fun _ : Fin D => ν).real Set.univ = 1 := by simp
  rw [hone, one_smul]
  have hstep : ∫ t, ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4
      ∂(Measure.pi fun _ : Fin D => ν)
      ≤ 2 * ((∫ y, y ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 := by
    refine hb4.trans ?_
    have hnu : (0 : ℝ) ≤ 2 * ((∫ y, y ^ 4 ∂ν) + 3) := by
      have := integral_pow_four_nonneg ν; linarith
    refine mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) ?_ 2) hnu
    exact div_le_div_of_nonneg_right hu (by positivity)
  have hAnn : (0 : ℝ) ≤ 2 * A ^ 2 := by positivity
  have hBc2 : Bc ^ 2 ≤ (1 / ((D : ℝ) * z.im)) ^ 2 := by
    refine pow_le_pow_left₀ hBnn ?_ 2
    rw [hBc]
    nlinarith [hu, huu, (by positivity : (0:ℝ) ≤ 1 / ((D : ℝ) * z.im))]
  rw [rowVarBound, ← hA]
  have hmain : 2 * (A ^ 2 * ∫ t, ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4
        ∂(Measure.pi fun _ : Fin D => ν))
      ≤ 4 * A ^ 2 * ((∫ y, y ^ 4 ∂ν) + 3) / ((D : ℝ) * z.im ^ 2) ^ 2 := by
    calc 2 * (A ^ 2 * ∫ t, ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4
            ∂(Measure.pi fun _ : Fin D => ν))
        = 2 * A ^ 2 * ∫ t, ‖R4C.cformC W z ((Real.sqrt D)⁻¹ • t) u‖ ^ 4
            ∂(Measure.pi fun _ : Fin D => ν) := by ring
      _ ≤ 2 * A ^ 2 * (2 * ((∫ y, y ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2) :=
          mul_le_mul_of_nonneg_left hstep hAnn
      _ = 4 * A ^ 2 * ((∫ y, y ^ 4 ∂ν) + 3) / ((D : ℝ) * z.im ^ 2) ^ 2 := by ring
  have hrest : 2 * Bc ^ 2 ≤ 2 / ((D : ℝ) * z.im) ^ 2 := by
    have : (2 : ℝ) / ((D : ℝ) * z.im) ^ 2 = 2 * (1 / ((D : ℝ) * z.im)) ^ 2 := by
      ring
    rw [this]
    linarith
  linarith

/-- The variance is at most the second moment about any fixed point.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_sq_sub_mean_le {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {g : β → ℝ} (h1 : Integrable g μ)
    (h2 : Integrable (fun t => g t ^ 2) μ) (a : ℝ) :
    ∫ t, (g t - ∫ s, g s ∂μ) ^ 2 ∂μ ≤ ∫ t, (g t - a) ^ 2 ∂μ := by
  have hone : μ.real Set.univ = 1 := by simp
  have hexp : ∀ b : ℝ, ∫ t, (g t - b) ^ 2 ∂μ
      = (∫ t, g t ^ 2 ∂μ) - 2 * b * (∫ t, g t ∂μ) + b ^ 2 := by
    intro b
    have hI2 : Integrable (fun t => -(2 * b) * g t) μ := h1.const_mul _
    have hI3 : Integrable (fun t : β => (b : ℝ) ^ 2) μ := integrable_const _
    have hI1 : Integrable (fun t => -(2 * b) * g t + b ^ 2) μ := hI2.add hI3
    have hfun : (fun t => (g t - b) ^ 2)
        = fun t => g t ^ 2 + (-(2 * b) * g t + b ^ 2) := by
      funext t; ring
    rw [hfun, integral_add h2 hI1, integral_add hI2 hI3, integral_const_mul, integral_const,
      hone, one_smul]
    ring
  rw [hexp, hexp]
  nlinarith [sq_nonneg ((∫ s, g s ∂μ) - a)]

/-- **The variance of one real component of `isoF`**, by Efron-Stein on the rows. -/
private theorem variance_component_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1)
    (φ : ℂ → ℝ) (hφm : Measurable φ) (hφs : ∀ w v : ℂ, φ w - φ v = φ (w - v))
    (hφn : ∀ w : ℂ, |φ w| ≤ ‖w‖) :
    ∫ Y, (φ (isoF z u Y) - ∫ X, φ (isoF z u X) ∂(noiseMatrix ν P D)) ^ 2
        ∂(noiseMatrix ν P D)
      ≤ (P : ℝ) * rowVarBound ν z D := by
  have hprob := hν.prob
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hf2 : ∀ {Q : ℕ} (Y : Matrix (Fin Q) (Fin D) ℝ), ‖isoF z u Y‖ ≤ 2 / z.im := by
    intro Q Y
    refine (norm_isoF_le hz hD u Y).trans ?_
    have h1 : 2 * (u ⬝ᵥ u) ≤ 2 := by linarith
    exact div_le_div_of_nonneg_right h1 hz.le
  set f : Matrix (Fin P) (Fin D) ℝ → ℝ := fun Y => φ (isoF z u Y) with hf
  have hfm : Measurable f := hφm.comp (measurable_isoF z u)
  have hbd : ∀ Y : Matrix (Fin P) (Fin D) ℝ, |f Y| ≤ 2 / z.im := fun Y =>
    (hφn _).trans (hf2 Y)
  have hfsq : Integrable (fun Y => f Y ^ 2) (noiseMatrix ν P D) := by
    refine Integrable.mono' (integrable_const ((2 / z.im) ^ 2))
      (hfm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbd Y) 2
  have hmem : MemLp f 2 (noiseMatrix ν P D) :=
    (memLp_two_iff_integrable_sq hfm.aestronglyMeasurable).mpr hfsq
  have hES : ∫ Y, (f Y - ∫ X, f X ∂(noiseMatrix ν P D)) ^ 2 ∂(noiseMatrix ν P D)
      ≤ ∑ k : Fin P, ∫ x, (∫ t, (f (Matrix.updateRow x k t)
            - ∫ s, f (Matrix.updateRow x k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
            ∂(Measure.pi fun _ : Fin D => ν)) ∂(noiseMatrix ν P D) :=
    Tensorization.variance_pi_le_sum (fun _ : Fin P => Measure.pi fun _ : Fin D => ν) f hmem
  refine hES.trans ?_
  have hinner : ∀ (k : Fin P) (x : Matrix (Fin P) (Fin D) ℝ),
      ∫ t, (f (Matrix.updateRow x k t)
          - ∫ s, f (Matrix.updateRow x k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
          ∂(Measure.pi fun _ : Fin D => ν)
        ≤ rowVarBound ν z D := by
    intro k x
    have hWm : ∀ i j, Measurable fun t : Fin D → ℝ => gram (Matrix.updateRow x k t) i j :=
      fun i j => measurable_gram_entry (fun l i' => measurable_updateRow_entry' x k l i') i j
    have hisoM : Measurable fun t : Fin D → ℝ => isoF z u (Matrix.updateRow x k t) := by
      simp only [isoF]
      refine Measurable.sub ?_ (measurable_const.mul (measurable_stieltjesC hWm z))
      exact measurable_qformC (y := fun _ : Fin D → ℝ => u) hWm z (fun i => measurable_const)
    have hgm : Measurable fun t : Fin D → ℝ => f (Matrix.updateRow x k t) := hφm.comp hisoM
    have hgbd : ∀ t : Fin D → ℝ, |f (Matrix.updateRow x k t)| ≤ 2 / z.im := fun t =>
      (hφn _).trans (hf2 _)
    have hgint : Integrable (fun t : Fin D → ℝ => f (Matrix.updateRow x k t))
        (Measure.pi fun _ : Fin D => ν) :=
      Integrable.mono' (integrable_const (2 / z.im)) hgm.aestronglyMeasurable
        (Filter.Eventually.of_forall fun t => by rw [Real.norm_eq_abs]; exact hgbd t)
    have hgsq : Integrable (fun t : Fin D → ℝ => f (Matrix.updateRow x k t) ^ 2)
        (Measure.pi fun _ : Fin D => ν) := by
      refine Integrable.mono' (integrable_const ((2 / z.im) ^ 2))
        (hgm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun t => ?_)
      rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
      exact pow_le_pow_left₀ (abs_nonneg _) (hgbd t) 2
    have hdm : Measurable fun t : Fin D → ℝ =>
        ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2 :=
      ((hisoM.sub measurable_const).norm).pow_const 2
    have hdbd : ∀ t : Fin D → ℝ,
        ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2
          ≤ (4 / z.im) ^ 2 := by
      intro t
      refine pow_le_pow_left₀ (norm_nonneg _) ?_ 2
      refine (norm_sub_le _ _).trans ?_
      have h1 := hf2 (Matrix.updateRow x k t)
      have h2 := hf2 (Matrix.updateRow x k (0 : Fin D → ℝ))
      have h3 : (4 : ℝ) / z.im = 2 / z.im + 2 / z.im := by ring
      rw [h3]; linarith
    have hdint : Integrable (fun t : Fin D → ℝ =>
        ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2)
        (Measure.pi fun _ : Fin D => ν) := by
      refine Integrable.mono' (integrable_const ((4 / z.im) ^ 2))
        hdm.aestronglyMeasurable (Filter.Eventually.of_forall fun t => ?_)
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      exact hdbd t
    refine (integral_sq_sub_mean_le hgint hgsq (f (Matrix.updateRow x k 0))).trans ?_
    refine (integral_mono_of_nonneg (Filter.Eventually.of_forall fun t => sq_nonneg _)
      hdint (Filter.Eventually.of_forall fun t => ?_)).trans
      (integral_normSq_isoF_sub_le hν hz hD hu k x)
    rw [hf, hφs]
    calc φ (isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)) ^ 2
        = |φ (isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0))| ^ 2 := by
          rw [sq_abs]
      _ ≤ ‖isoF z u (Matrix.updateRow x k t) - isoF z u (Matrix.updateRow x k 0)‖ ^ 2 :=
          pow_le_pow_left₀ (abs_nonneg _) (hφn _) 2
  have hone : (noiseMatrix ν P D).real Set.univ = 1 := by simp
  have hterm : ∀ k : Fin P,
      ∫ x, (∫ t, (f (Matrix.updateRow x k t)
          - ∫ s, f (Matrix.updateRow x k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
          ∂(Measure.pi fun _ : Fin D => ν)) ∂(noiseMatrix ν P D)
        ≤ rowVarBound ν z D := by
    intro k
    have hnn : ∀ x : Matrix (Fin P) (Fin D) ℝ, (0 : ℝ)
        ≤ ∫ t, (f (Matrix.updateRow x k t)
            - ∫ s, f (Matrix.updateRow x k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
            ∂(Measure.pi fun _ : Fin D => ν) := fun x =>
      integral_nonneg (μ := Measure.pi fun _ : Fin D => ν) fun t => sq_nonneg _
    have hmono := integral_mono_of_nonneg (Filter.Eventually.of_forall hnn)
      (integrable_const (μ := noiseMatrix ν P D) (rowVarBound ν z D))
      (Filter.Eventually.of_forall (hinner k))
    refine hmono.trans (le_of_eq ?_)
    rw [integral_const, hone, one_smul]
  calc ∑ k : Fin P, ∫ x, (∫ t, (f (Matrix.updateRow x k t)
        - ∫ s, f (Matrix.updateRow x k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
        ∂(Measure.pi fun _ : Fin D => ν)) ∂(noiseMatrix ν P D)
      ≤ ∑ _k : Fin P, rowVarBound ν z D := Finset.sum_le_sum fun k _ => hterm k
    _ = (P : ℝ) * rowVarBound ν z D := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-! ### Milestone 5: the mean

The self-consistent equation. Section 2.3 of the brief: the variance alone cannot give the
mean, because the law of `noiseMatrix ν p d` is not invariant under a sign flip of a column,
so `E G` is not a multiple of the identity. -/

/-- The scalar of the self-consistent equation: `secW c z = (1 + m)⁻¹`, written as `-z m̃` so
that the bound `‖secW‖ ≤ ‖z‖/η` is free of `c`. It is the deterministic twin of
`1 - s k = -z (Gc) k k`. -/
noncomputable def secW (c : ℝ) (z : ℂ) : ℂ := -(z * MP.mTildeC c z)

/-- `(1 + m) secW = 1`. -/
theorem one_add_mC_mul_secW {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    (1 + MP.mC c z) * secW c z = 1 := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have h1 := MP.quad_mTildeC hc hz
  have h2 := MP.c_mul_mTildeC (z := z) hc
  have h2' : (c : ℂ) * z * MP.mTildeC c z = z * MP.mC c z + (1 - (c : ℂ)) := by
    field_simp at h2
    linear_combination h2
  have key : z * MP.mTildeC c z * (1 + MP.mC c z)
      = ((c : ℂ) * z * MP.mTildeC c z ^ 2 + (z + (c : ℂ) - 1) * MP.mTildeC c z + 1) - 1 := by
    linear_combination (-(MP.mTildeC c z)) * h2'
  rw [secW]
  linear_combination -key - h1

/-- `‖secW c z‖ ≤ ‖z‖/η`, with no `c` in the bound. -/
theorem norm_secW_le {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    ‖secW c z‖ ≤ ‖z‖ / z.im := by
  rw [secW, norm_neg, norm_mul]
  have h := MP.norm_mTildeC_le hc hz
  calc ‖z‖ * ‖MP.mTildeC c z‖ ≤ ‖z‖ * (z.im)⁻¹ :=
        mul_le_mul_of_nonneg_left h (norm_nonneg z)
    _ = ‖z‖ / z.im := by rw [div_eq_mul_inv]

/-- **The algebra that makes the unit work.** `m (z - c secW) = -1`, so the coefficient of the
mean in the self-consistent equation has norm at least `η`. -/
theorem mC_mul_sub_secW {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    MP.mC c z * (z - (c : ℂ) * secW c z) = -1 := by
  have hw := one_add_mC_mul_secW hc hz
  have hne : (1 + MP.mC c z) ≠ 0 := by
    intro h; rw [h, zero_mul] at hw; exact zero_ne_one hw
  have hquad : MP.quad c z (MP.mC c z) = 0 := MP.quad_mC hc.le hz
  have hquad' : z * MP.mC c z ^ 2 + (z + 1 - (c : ℂ)) * MP.mC c z + 1 = 0 := hquad
  have hzero : (1 + MP.mC c z) * (MP.mC c z * (z - (c : ℂ) * secW c z) + 1) = 0 := by
    linear_combination hquad' + (-(c : ℂ) * MP.mC c z) * hw
  have h3 := (mul_eq_zero.mp hzero).resolve_left hne
  linear_combination h3

/-- **The second moment of a product of two linear forms**, by polarization from
`LinForm.integral_linForm_sq`. -/
private theorem integral_linForm_mul {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (a b : Fin D → ℝ) :
    Integrable (fun t => LinForm.linForm a t * LinForm.linForm b t)
        (Measure.pi fun _ : Fin D => ν) ∧
    ∫ t, LinForm.linForm a t * LinForm.linForm b t ∂(Measure.pi fun _ : Fin D => ν)
      = ∑ l, a l * b l := by
  have hprob := hν.prob
  have hsq : ∀ c : Fin D → ℝ,
      Integrable (fun t => LinForm.linForm c t ^ 2) (Measure.pi fun _ : Fin D => ν) := fun c =>
    LinForm.integrable_linForm_pow_of_four
      (LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 c) (by norm_num)
  have hval : ∀ c : Fin D → ℝ,
      ∫ t, LinForm.linForm c t ^ 2 ∂(Measure.pi fun _ : Fin D => ν) = ∑ l, c l ^ 2 := fun c =>
    LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 c
  have hadd : ∀ t : Fin D → ℝ, LinForm.linForm (a + b) t
      = LinForm.linForm a t + LinForm.linForm b t := by
    intro t
    simp only [LinForm.linForm, Pi.add_apply, add_mul]
    exact Finset.sum_add_distrib
  have hdomI : Integrable (fun t : Fin D → ℝ =>
      1 / 2 * (LinForm.linForm a t ^ 2 + LinForm.linForm b t ^ 2))
      (Measure.pi fun _ : Fin D => ν) := ((hsq a).add (hsq b)).const_mul (1 / 2)
  have hmul : Integrable (fun t => LinForm.linForm a t * LinForm.linForm b t)
      (Measure.pi fun _ : Fin D => ν) := by
    refine Integrable.mono' hdomI
      (((LinForm.measurable_linForm a).mul
        (LinForm.measurable_linForm b)).aestronglyMeasurable)
      (Filter.Eventually.of_forall fun t => ?_)
    rw [Real.norm_eq_abs, abs_mul]
    nlinarith [sq_nonneg (|LinForm.linForm a t| - |LinForm.linForm b t|),
      abs_nonneg (LinForm.linForm a t), abs_nonneg (LinForm.linForm b t),
      sq_abs (LinForm.linForm a t), sq_abs (LinForm.linForm b t)]
  refine ⟨hmul, ?_⟩
  have hfun : (fun t => LinForm.linForm (a + b) t ^ 2)
      = fun t => LinForm.linForm a t ^ 2
          + (2 * (LinForm.linForm a t * LinForm.linForm b t) + LinForm.linForm b t ^ 2) := by
    funext t; rw [hadd]; ring
  have hexp : ∫ t, LinForm.linForm (a + b) t ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      = (∫ t, LinForm.linForm a t ^ 2 ∂(Measure.pi fun _ : Fin D => ν))
        + (2 * ∫ t, LinForm.linForm a t * LinForm.linForm b t
            ∂(Measure.pi fun _ : Fin D => ν)
          + ∫ t, LinForm.linForm b t ^ 2 ∂(Measure.pi fun _ : Fin D => ν)) := by
    have hI2 : Integrable (fun t : Fin D → ℝ =>
        2 * (LinForm.linForm a t * LinForm.linForm b t))
        (Measure.pi fun _ : Fin D => ν) := hmul.const_mul 2
    have hI3 : Integrable (fun t : Fin D → ℝ =>
        2 * (LinForm.linForm a t * LinForm.linForm b t) + LinForm.linForm b t ^ 2)
        (Measure.pi fun _ : Fin D => ν) := hI2.add (hsq b)
    rw [hfun, integral_add (hsq a) hI3, integral_add hI2 (hsq b), integral_const_mul]
  rw [hval, hval, hval] at hexp
  have hsum : ∑ l, (a + b) l ^ 2
      = ∑ l, a l ^ 2 + (2 * (∑ l, a l * b l) + ∑ l, b l ^ 2) := by
    have hterm : ∀ l, (a + b) l ^ 2 = a l ^ 2 + (2 * (a l * b l) + b l ^ 2) := by
      intro l; simp only [Pi.add_apply]; ring
    rw [Finset.sum_congr rfl fun l _ => hterm l, Finset.sum_add_distrib,
      Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [hsum] at hexp
  linarith

/-- The same with a complex coefficient vector on the right. -/
private theorem integral_linForm_mul_clinForm {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (u : Fin D → ℝ) (a : Fin D → ℂ) :
    ∫ t, ((LinForm.linForm u t : ℝ) : ℂ) * (∑ l, a l * (t l : ℂ))
        ∂(Measure.pi fun _ : Fin D => ν)
      = ∑ l, (u l : ℂ) * a l := by
  have hprob := hν.prob
  set ar : Fin D → ℝ := fun l => (a l).re with har
  set ai : Fin D → ℝ := fun l => (a l).im with hai
  have hre : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).re = LinForm.linForm ar t := by
    intro t
    rw [Complex.re_sum]
    exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_re, har]
  have him : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).im = LinForm.linForm ai t := by
    intro t
    rw [Complex.im_sum]
    exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_im, hai]
  obtain ⟨hIr, hVr⟩ := integral_linForm_mul hν u ar
  obtain ⟨hIi, hVi⟩ := integral_linForm_mul hν u ai
  have hdecomp : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ))
      = ((LinForm.linForm ar t : ℝ) : ℂ)
        + Complex.I * ((LinForm.linForm ai t : ℝ) : ℂ) := by
    intro t
    apply Complex.ext <;> simp [hre t, him t]
  have hfun : (fun t : Fin D → ℝ => ((LinForm.linForm u t : ℝ) : ℂ) * (∑ l, a l * (t l : ℂ)))
      = fun t => ((LinForm.linForm u t * LinForm.linForm ar t : ℝ) : ℂ)
          + Complex.I * ((LinForm.linForm u t * LinForm.linForm ai t : ℝ) : ℂ) := by
    funext t
    rw [hdecomp t]
    push_cast
    ring
  have hIr' : Integrable (fun t : Fin D → ℝ =>
      ((LinForm.linForm u t * LinForm.linForm ar t : ℝ) : ℂ))
      (Measure.pi fun _ : Fin D => ν) := hIr.ofReal
  have hIi' : Integrable (fun t : Fin D → ℝ =>
      Complex.I * ((LinForm.linForm u t * LinForm.linForm ai t : ℝ) : ℂ))
      (Measure.pi fun _ : Fin D => ν) := hIi.ofReal.const_mul Complex.I
  rw [hfun, integral_add hIr' hIi', integral_const_mul, integral_complex_ofReal,
    integral_complex_ofReal, hVr, hVi, Complex.ofReal_sum, Complex.ofReal_sum,
    Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun l _ => ?_
  have hal : a l = ((ar l : ℝ) : ℂ) + Complex.I * ((ai l : ℝ) : ℂ) := by
    apply Complex.ext <;> simp [har, hai]
  rw [hal]
  push_cast
  ring


/-- The bilinear form of the resolvent, at measurable families, is measurable. -/
theorem measurable_cformC {α : Type*} [MeasurableSpace α] {D : ℕ}
    {W : α → Matrix (Fin D) (Fin D) ℝ} (hW : ∀ i j, Measurable fun a => W a i j) (z : ℂ)
    {x y : α → (Fin D → ℝ)} (hx : ∀ i, Measurable fun a => x a i)
    (hy : ∀ i, Measurable fun a => y a i) :
    Measurable fun a => R4C.cformC (W a) z (x a) (y a) := by
  have h : ∀ a, R4C.cformC (W a) z (x a) (y a)
      = ∑ i, ∑ j, R4C.resolvC (W a) z i j * ((x a i : ℝ) : ℂ) * ((y a j : ℝ) : ℂ) :=
    fun a => bil_eq_sum (R4C.resolvC (W a) z) (x a) (y a)
  simp only [h]
  refine Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ => ?_
  exact ((measurable_resolvC_entry hW z i j).mul
    (Complex.measurable_ofReal.comp (hx i))).mul (Complex.measurable_ofReal.comp (hy j))

/-- Cauchy-Schwarz for the real dot product. -/
theorem sq_dotProduct_le {D : ℕ} (x y : Fin D → ℝ) :
    (x ⬝ᵥ y) ^ 2 ≤ (x ⬝ᵥ x) * (y ⬝ᵥ y) := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin D)) x y
  have hx : ∑ i, x i ^ 2 = x ⬝ᵥ x := by
    rw [dotProduct]; exact Finset.sum_congr rfl fun i _ => by ring
  have hy : ∑ i, y i ^ 2 = y ⬝ᵥ y := by
    rw [dotProduct]; exact Finset.sum_congr rfl fun i _ => by ring
  rw [hx, hy] at h
  exact h

/-- **The row term of the identity `(*)`**, at the leave-one-out resolvent. -/
noncomputable def rowTerm (z : ℂ) (u : Fin d → ℝ) (Y : Matrix (Fin p) (Fin d) ℝ)
    (k : Fin p) : ℂ :=
  ((u ⬝ᵥ rowVec Y k : ℝ) : ℂ)
    * R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u

theorem measurable_rowTerm {P D : ℕ} (z : ℂ) (u : Fin D → ℝ) (k : Fin P) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => rowTerm z u Y k := by
  simp only [rowTerm]
  refine Measurable.mul ?_ ?_
  · refine Complex.measurable_ofReal.comp ?_
    simp only [dotProduct]
    exact Finset.measurable_sum _ fun i _ => measurable_const.mul (measurable_rowVec k i)
  · exact measurable_cformC (measurable_gram_loo k) z (measurable_rowVec k)
      (fun i => measurable_const)

/-- The identity `(*)` in the `rowTerm` notation. -/
theorem z_mul_qformC_add_eq_rowTerm (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im)
    (u : Fin d → ℝ) :
    z * R4C.qformC (gram Y) z u + ((u ⬝ᵥ u : ℝ) : ℂ)
      = ∑ k, (1 - sRow Y z k) * rowTerm z u Y k :=
  z_mul_qformC_add_eq_loo Y hz u

/-- The row term is bounded by the squared length of the row, which is integrable. -/
theorem norm_rowTerm_le (hz : 0 < z.im) (u : Fin d → ℝ) (hu : u ⬝ᵥ u ≤ 1)
    (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) :
    ‖rowTerm z u Y k‖ ≤ (rowVec Y k ⬝ᵥ rowVec Y k) / z.im := by
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hgg : (0 : ℝ) ≤ rowVec Y k ⬝ᵥ rowVec Y k := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hcs : (u ⬝ᵥ rowVec Y k) ^ 2 ≤ (u ⬝ᵥ u) * (rowVec Y k ⬝ᵥ rowVec Y k) :=
    sq_dotProduct_le u (rowVec Y k)
  have h3 : |u ⬝ᵥ rowVec Y k|
      ≤ Real.sqrt (u ⬝ᵥ u) * Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k) := by
    have h := Real.sqrt_le_sqrt hcs
    rwa [Real.sqrt_sq_eq_abs, Real.sqrt_mul huu] at h
  have h2 : ‖R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u‖
      ≤ Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k) * Real.sqrt (u ⬝ᵥ u) / z.im :=
    R4C.norm_cformC_le (gram_isHermitian _) hz (rowVec Y k) u
  calc ‖rowTerm z u Y k‖
      = |u ⬝ᵥ rowVec Y k| * ‖R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k) u‖ := by
        rw [rowTerm, norm_mul, Complex.norm_real, Real.norm_eq_abs]
    _ ≤ (Real.sqrt (u ⬝ᵥ u) * Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k))
          * (Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k) * Real.sqrt (u ⬝ᵥ u) / z.im) :=
        mul_le_mul h3 h2 (norm_nonneg _) (by positivity)
    _ = (Real.sqrt (u ⬝ᵥ u) * Real.sqrt (u ⬝ᵥ u))
          * (Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k) * Real.sqrt (rowVec Y k ⬝ᵥ rowVec Y k))
          / z.im := by ring
    _ = (u ⬝ᵥ u) * (rowVec Y k ⬝ᵥ rowVec Y k) / z.im := by
        rw [Real.mul_self_sqrt huu, Real.mul_self_sqrt hgg]
    _ ≤ (rowVec Y k ⬝ᵥ rowVec Y k) / z.im := by
        have hle : (u ⬝ᵥ u) * (rowVec Y k ⬝ᵥ rowVec Y k) ≤ rowVec Y k ⬝ᵥ rowVec Y k := by
          nlinarith
        exact div_le_div_of_nonneg_right hle hz.le

/-- The squared length of a row is integrable under the product law. -/
private theorem integrable_rowVec_sq {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (k : Fin P) :
    Integrable (fun Y : Matrix (Fin P) (Fin D) ℝ => rowVec Y k ⬝ᵥ rowVec Y k)
      (noiseMatrix ν P D) := by
  have hprob := hν.prob
  have h1 : Integrable (fun x : ℝ => x ^ 2) ν :=
    LinForm.integrable_pow_of_four hν.mom4 (by norm_num)
  have h2 : ∀ l : Fin D, Integrable (fun y : Fin D → ℝ => (y l) ^ 2)
      (Measure.pi fun _ : Fin D => ν) := fun l =>
    (measurePreserving_eval (fun _ : Fin D => ν) l).integrable_comp_of_integrable h1
  have h3 : ∀ l : Fin D, Integrable
      (fun Y : Matrix (Fin P) (Fin D) ℝ => (Y k l) ^ 2) (noiseMatrix ν P D) := fun l =>
    (measurePreserving_eval (fun _ : Fin P => Measure.pi fun _ : Fin D => ν) k)
      |>.integrable_comp_of_integrable (h2 l)
  have hdd : (Real.sqrt D)⁻¹ * (Real.sqrt D)⁻¹ = ((D : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg D)]
  have hfun : (fun Y : Matrix (Fin P) (Fin D) ℝ => rowVec Y k ⬝ᵥ rowVec Y k)
      = fun Y => ∑ l, ((D : ℝ))⁻¹ * (Y k l) ^ 2 := by
    funext Y
    rw [dotProduct]
    refine Finset.sum_congr rfl fun l _ => ?_
    simp only [rowVec_apply]
    rw [← hdd]
    ring
  rw [hfun]
  exact integrable_finsetSum _ fun l _ => (h3 l).const_mul _

/-- **The Fubini split on the rows**, in the Bochner form.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_noiseMatrix_split {n D : ℕ} {ν : Measure ℝ} [IsProbabilityMeasure ν]
    {F : Matrix (Fin (n + 1)) (Fin D) ℝ → ℂ} (hF : Integrable F (noiseMatrix ν (n + 1) D))
    (k : Fin (n + 1)) :
    ∫ Y, F Y ∂(noiseMatrix ν (n + 1) D)
      = ∫ r, (∫ t, F (Fin.insertNth k t r) ∂(Measure.pi fun _ : Fin D => ν))
          ∂(Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin D => ν) := by
  set ρ : Measure (Fin D → ℝ) := Measure.pi fun _ : Fin D => ν with hρ
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => (Fin D → ℝ)) k with he
  have hmp : MeasurePreserving e (noiseMatrix ν (n + 1) D)
      (ρ.prod (Measure.pi fun _ : Fin n => ρ)) :=
    measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) => ρ) k
  have hsym := hmp.symm e
  have hint : Integrable (fun q => F (e.symm q)) (ρ.prod (Measure.pi fun _ : Fin n => ρ)) :=
    hsym.integrable_comp_of_integrable hF
  have h1 : ∫ q, F (e.symm q) ∂(ρ.prod (Measure.pi fun _ : Fin n => ρ))
      = ∫ Y, F Y ∂(noiseMatrix ν (n + 1) D) := hsym.integral_comp' F
  rw [← h1, integral_prod_symm _ hint]
  rfl

/-- Zeroing the inserted row gives back the matrix with a zero row.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem updateRow_insertNth {n D : ℕ} (k : Fin (n + 1)) (r : Fin n → Fin D → ℝ)
    (t : Fin D → ℝ) :
    Matrix.updateRow (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k 0
      = (Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r
          : Matrix (Fin (n + 1)) (Fin D) ℝ) := by
  change Function.update (Fin.insertNth k t r) k (0 : Fin D → ℝ)
    = Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r
  simp

/-- The scaled inserted row. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem rowVec_insertNth {n D : ℕ} (k : Fin (n + 1)) (r : Fin n → Fin D → ℝ)
    (t : Fin D → ℝ) :
    rowVec (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k = (Real.sqrt D)⁻¹ • t := by
  have h2 : (Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k t r) k = t := by simp
  change (Real.sqrt D)⁻¹ • ((Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k t r) k)
    = (Real.sqrt D)⁻¹ • t
  rw [h2]

/-- `rowTerm` is integrable under the product law. -/
private theorem integrable_rowTerm {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin P) :
    Integrable (fun Y : Matrix (Fin P) (Fin D) ℝ => rowTerm z u Y k) (noiseMatrix ν P D) := by
  have hprob := hν.prob
  refine Integrable.mono' ((integrable_rowVec_sq hν k).div_const z.im)
    (measurable_rowTerm z u k).aestronglyMeasurable
    (Filter.Eventually.of_forall fun Y => ?_)
  exact norm_rowTerm_le hz u hu Y k

/-- The leave-one-out quadratic form is bounded, hence integrable. -/
private theorem integrable_qformC_loo {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin P) :
    Integrable (fun Y : Matrix (Fin P) (Fin D) ℝ =>
      R4C.qformC (gram (Matrix.updateRow Y k 0)) z u) (noiseMatrix ν P D) := by
  have hprob := hν.prob
  have hm : Measurable fun Y : Matrix (Fin P) (Fin D) ℝ =>
      R4C.qformC (gram (Matrix.updateRow Y k 0)) z u :=
    measurable_qformC (y := fun _ : Matrix (Fin P) (Fin D) ℝ => u)
      (measurable_gram_loo k) z (fun i => measurable_const)
  refine Integrable.mono' (integrable_const (μ := noiseMatrix ν P D) (1 / z.im))
    hm.aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
  refine (R4C.norm_qformC_le (gram_isHermitian _) hz u).trans ?_
  exact div_le_div_of_nonneg_right hu hz.le

/-- **Step 1 of section 2.3: the conditional mean is exact.** The row term integrates to
`D⁻¹` times the leave-one-out quadratic form, with no error at all. -/
private theorem integral_rowTerm_eq {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin (n + 1)) :
    ∫ Y, rowTerm z u Y k ∂(noiseMatrix ν (n + 1) D)
      = ((D : ℂ))⁻¹ * ∫ Y, R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
          ∂(noiseMatrix ν (n + 1) D) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hone : (Measure.pi fun _ : Fin D => ν).real Set.univ = 1 := by simp
  have hsq : ((Real.sqrt D)⁻¹ : ℝ) * ((Real.sqrt D)⁻¹ : ℝ) = ((D : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt hDR.le]
  rw [integral_noiseMatrix_split (integrable_rowTerm hν hz hu k) k,
    integral_noiseMatrix_split (integrable_qformC_loo hν hz hu k) k, ← integral_const_mul]
  refine integral_congr_ae (Filter.Eventually.of_forall fun r => ?_)
  set X : Matrix (Fin (n + 1)) (Fin D) ℝ :=
    Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r with hX
  set a : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l with ha
  have hterm : ∀ t : Fin D → ℝ,
      rowTerm z u (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k
        = (((Real.sqrt D)⁻¹ : ℝ) : ℂ)
          * (((LinForm.linForm u t : ℝ) : ℂ) * (∑ l, a l * (t l : ℂ))) := by
    intro t
    simp only [rowTerm]
    rw [updateRow_insertNth k r t, ← hX, rowVec_insertNth k r t, cformC_smul_left_eq]
    simp only [ha]
    have hdot : (u ⬝ᵥ ((Real.sqrt D)⁻¹ • t) : ℝ)
        = ((Real.sqrt D)⁻¹ : ℝ) * LinForm.linForm u t := by
      rw [dotProduct, LinForm.linForm, Finset.mul_sum]
      exact Finset.sum_congr rfl fun l _ => by
        simp only [Pi.smul_apply, smul_eq_mul]; ring
    rw [hdot, Complex.ofReal_mul]
    ring
  have hq : ∀ t : Fin D → ℝ,
      R4C.qformC (gram (Matrix.updateRow
          (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k 0)) z u
        = R4C.qformC (gram X) z u := by
    intro t
    rw [updateRow_insertNth k r t, ← hX]
  simp only [hterm, hq]
  rw [integral_const_mul, integral_linForm_mul_clinForm hν u a, integral_const, hone, one_smul]
  have hqsum : R4C.qformC (gram X) z u
      = ∑ l, (u l : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l := rfl
  rw [hqsum]
  simp only [ha]
  rw [Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun l _ => ?_
  have hcast : (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (((Real.sqrt D)⁻¹ : ℝ) : ℂ) = ((D : ℂ))⁻¹ := by
    rw [← Complex.ofReal_mul, hsq, Complex.ofReal_inv]
    norm_num
  calc (((Real.sqrt D)⁻¹ : ℝ) : ℂ)
        * ((u l : ℂ) * ((((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l))
      = ((((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (((Real.sqrt D)⁻¹ : ℝ) : ℂ))
          * ((u l : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l) := by ring
    _ = ((D : ℂ))⁻¹ * ((u l : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l) := by rw [hcast]


/-- **The second moment of a complex linear form** is the squared length of its coefficient
vector. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_normSq_clinForm {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (a : Fin D → ℂ) :
    ∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 2 ∂(Measure.pi fun _ : Fin D => ν) = ∑ l, ‖a l‖ ^ 2 := by
  have hprob := hν.prob
  set ar : Fin D → ℝ := fun l => (a l).re with har
  set ai : Fin D → ℝ := fun l => (a l).im with hai
  have hre : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).re = LinForm.linForm ar t := by
    intro t
    rw [Complex.re_sum]
    exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_re, har]
  have him : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ)).im = LinForm.linForm ai t := by
    intro t
    rw [Complex.im_sum]
    exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_im, hai]
  have hsq : ∀ c : Fin D → ℝ,
      Integrable (fun t => LinForm.linForm c t ^ 2) (Measure.pi fun _ : Fin D => ν) := fun c =>
    LinForm.integrable_linForm_pow_of_four
      (LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 c) (by norm_num)
  have hfun : (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 2)
      = fun t => LinForm.linForm ar t ^ 2 + LinForm.linForm ai t ^ 2 := by
    funext t
    rw [Complex.sq_norm, Complex.normSq_apply, hre t, him t]; ring
  rw [hfun, integral_add (hsq ar) (hsq ai),
    LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 ar,
    LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 ai, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun l _ => by
    rw [Complex.sq_norm, Complex.normSq_apply, har, hai]; ring

/-- The squared norm of a complex linear form is integrable.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integrable_normSq_clinForm {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (a : Fin D → ℂ) :
    Integrable (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 2)
      (Measure.pi fun _ : Fin D => ν) := by
  have hprob := hν.prob
  have hsq : ∀ c : Fin D → ℝ,
      Integrable (fun t => LinForm.linForm c t ^ 2) (Measure.pi fun _ : Fin D => ν) := fun c =>
    LinForm.integrable_linForm_pow_of_four
      (LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 c) (by norm_num)
  have hfun : (fun t : Fin D → ℝ => ‖∑ l, a l * (t l : ℂ)‖ ^ 2)
      = fun t => LinForm.linForm (fun l => (a l).re) t ^ 2
          + LinForm.linForm (fun l => (a l).im) t ^ 2 := by
    funext t
    have hre : (∑ l, a l * (t l : ℂ)).re = LinForm.linForm (fun l => (a l).re) t := by
      rw [Complex.re_sum]
      exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_re]
    have him : (∑ l, a l * (t l : ℂ)).im = LinForm.linForm (fun l => (a l).im) t := by
      rw [Complex.im_sum]
      exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_im]
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]; ring
  rw [hfun]
  exact (hsq _).add (hsq _)

/-- The `lintegral` form of the second moment of a complex linear form.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem lintegral_normSq_clinForm {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (a : Fin D → ℂ) :
    ∫⁻ t, ENNReal.ofReal (‖∑ l, a l * (t l : ℂ)‖ ^ 2) ∂(Measure.pi fun _ : Fin D => ν)
      = ENNReal.ofReal (∑ l, ‖a l‖ ^ 2) := by
  rw [← ofReal_integral_eq_lintegral_ofReal (integrable_normSq_clinForm hν a)
    (Filter.Eventually.of_forall fun t => sq_nonneg _), integral_normSq_clinForm hν a]

/-- **Step 3 of section 2.3.** The leave-one-out change of the quadratic form is `O(1/D)` in
`L¹`. -/
private theorem lintegral_norm_qformC_sub_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin (n + 1)) :
    ∫⁻ Y, ENNReal.ofReal (‖R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
        - R4C.qformC (gram Y) z u‖) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (‖z‖ / z.im * (1 / ((D : ℝ) * z.im ^ 2))) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hsq : ((Real.sqrt D)⁻¹ : ℝ) ^ 2 = ((D : ℝ))⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hDR.le]
  have hmF : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      ENNReal.ofReal (‖R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
        - R4C.qformC (gram Y) z u‖) :=
    ENNReal.measurable_ofReal.comp
      (((measurable_qformC (y := fun _ : Matrix (Fin (n + 1)) (Fin D) ℝ => u)
          (measurable_gram_loo k) z (fun i => measurable_const)).sub
        (measurable_qformC (y := fun _ : Matrix (Fin (n + 1)) (Fin D) ℝ => u)
          measurable_gram_self z (fun i => measurable_const))).norm)
  refine lintegral_noiseMatrix_le hmF k fun r => ?_
  set X : Matrix (Fin (n + 1)) (Fin D) ℝ :=
    Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r with hX
  set a : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l with ha
  have hcoef : ∑ l, ‖a l‖ ^ 2 ≤ 1 / ((D : ℝ) * z.im ^ 2) := by
    have hstep : ∀ l, ‖a l‖ ^ 2
        = ((D : ℝ))⁻¹ * ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l‖ ^ 2 := by
      intro l
      rw [ha]
      simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs, hsq]
    rw [Finset.sum_congr rfl fun l _ => hstep l, ← Finset.mul_sum]
    have h := sum_normSq_resolvC_mulVec_le (gram_isHermitian X) hz u
    calc ((D : ℝ))⁻¹ * ∑ l, ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l‖ ^ 2
        ≤ ((D : ℝ))⁻¹ * ((u ⬝ᵥ u) / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h (by positivity)
      _ ≤ ((D : ℝ))⁻¹ * (1 / z.im ^ 2) := by
          have : (u ⬝ᵥ u) / z.im ^ 2 ≤ 1 / z.im ^ 2 :=
            div_le_div_of_nonneg_right hu (by positivity)
          exact mul_le_mul_of_nonneg_left this (by positivity)
      _ = 1 / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hpt : ∀ t : Fin D → ℝ,
      ENNReal.ofReal (‖R4C.qformC (gram (Matrix.updateRow
          (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k 0)) z u
          - R4C.qformC (gram (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ)) z u‖)
        ≤ ENNReal.ofReal (‖z‖ / z.im) * ENNReal.ofReal (‖∑ l, a l * (t l : ℂ)‖ ^ 2) := by
    intro t
    set Y : Matrix (Fin (n + 1)) (Fin D) ℝ := Fin.insertNth k t r with hY
    have hupd : Matrix.updateRow Y k 0 = X := by rw [hY, updateRow_insertNth k r t, ← hX]
    have hrow : rowVec Y k = (Real.sqrt D)⁻¹ • t := by rw [hY]; exact rowVec_insertNth k r t
    have hq := qformC_sub_qformC_loo Y hz u k
    rw [hupd] at hq
    have hcomm : R4C.cformC (gram X) z u (rowVec Y k) = R4C.cformC (gram X) z (rowVec Y k) u :=
      cformC_comm (gram_isHermitian X) hz.ne' u (rowVec Y k)
    rw [hcomm, hrow, cformC_smul_left_eq] at hq
    rw [hupd]
    simp only [ha]
    set L : ℂ := ∑ l, (((Real.sqrt D)⁻¹ : ℝ) : ℂ)
      * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l * (t l : ℂ) with hL
    have hval : ‖R4C.qformC (gram X) z u - R4C.qformC (gram Y) z u‖
        = ‖1 - sRow Y z k‖ * ‖L‖ ^ 2 := by
      have hneg : R4C.qformC (gram X) z u - R4C.qformC (gram Y) z u
          = (1 - sRow Y z k) * (L * L) := by linear_combination -hq
      rw [hneg, norm_mul, norm_mul]
      ring
    rw [hval, ← ENNReal.ofReal_mul (by positivity)]
    refine ENNReal.ofReal_le_ofReal ?_
    exact mul_le_mul_of_nonneg_right (norm_one_sub_sRow_le Y hz hD k) (by positivity)
  calc ∫⁻ t, ENNReal.ofReal (‖R4C.qformC (gram (Matrix.updateRow
          (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k 0)) z u
          - R4C.qformC (gram (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ)) z u‖)
        ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ∫⁻ t, ENNReal.ofReal (‖z‖ / z.im) * ENNReal.ofReal (‖∑ l, a l * (t l : ℂ)‖ ^ 2)
        ∂(Measure.pi fun _ : Fin D => ν) := lintegral_mono hpt
    _ = ENNReal.ofReal (‖z‖ / z.im) * ENNReal.ofReal (∑ l, ‖a l‖ ^ 2) := by
        have hmeasL : Measurable fun t : Fin D → ℝ =>
            ENNReal.ofReal (‖∑ l, a l * (t l : ℂ)‖ ^ 2) := by
          refine ENNReal.measurable_ofReal.comp (Measurable.pow_const (Measurable.norm ?_) 2)
          exact Finset.measurable_sum _ fun l _ =>
            measurable_const.mul (Complex.measurable_ofReal.comp (measurable_pi_apply l))
        rw [lintegral_const_mul _ hmeasL, lintegral_normSq_clinForm hν a]
    _ ≤ ENNReal.ofReal (‖z‖ / z.im * (1 / ((D : ℝ) * z.im ^ 2))) := by
        rw [← ENNReal.ofReal_mul (by positivity)]
        exact ENNReal.ofReal_le_ofReal
          (mul_le_mul_of_nonneg_left hcoef (by positivity))

/-- Additivity of `∫⁻ ∘ ofReal` on nonnegative measurable functions.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem lintegral_ofReal_add {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f g : α → ℝ} (hf : Measurable f) (hfn : ∀ a, 0 ≤ f a) (hgn : ∀ a, 0 ≤ g a) :
    ∫⁻ a, ENNReal.ofReal (f a + g a) ∂μ
      = (∫⁻ a, ENNReal.ofReal (f a) ∂μ) + ∫⁻ a, ENNReal.ofReal (g a) ∂μ := by
  have hcongr : ∀ a, ENNReal.ofReal (f a + g a)
      = ENNReal.ofReal (f a) + ENNReal.ofReal (g a) := fun a => ENNReal.ofReal_add (hfn a) (hgn a)
  rw [lintegral_congr hcongr]
  exact lintegral_add_left (ENNReal.measurable_ofReal.comp hf) _

/-- A constant comes out of `∫⁻ ∘ ofReal`. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem lintegral_ofReal_const_mul {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f : α → ℝ} (hf : Measurable f) {c : ℝ} (hc : 0 ≤ c) :
    ∫⁻ a, ENNReal.ofReal (c * f a) ∂μ
      = ENNReal.ofReal c * ∫⁻ a, ENNReal.ofReal (f a) ∂μ := by
  have hcongr : ∀ a, ENNReal.ofReal (c * f a)
      = ENNReal.ofReal c * ENNReal.ofReal (f a) := fun a => ENNReal.ofReal_mul hc
  rw [lintegral_congr hcongr]
  exact lintegral_const_mul _ (ENNReal.measurable_ofReal.comp hf)

/-- The per-row second moment of the row term, `O(1/D²)`. -/
noncomputable def rowTermSqBound (ν : Measure ℝ) (z : ℂ) (D : ℕ) : ℝ :=
  ((∫ x, x ^ 4 ∂ν) + 3) / (2 * (D : ℝ) ^ 2)
    + ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) ^ 2 * z.im ^ 4)

theorem rowTermSqBound_nonneg (ν : Measure ℝ) (z : ℂ) (D : ℕ) : 0 ≤ rowTermSqBound ν z D := by
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  rw [rowTermSqBound]
  have h1 : (0 : ℝ) ≤ ((∫ x, x ^ 4 ∂ν) + 3) / (2 * (D : ℝ) ^ 2) :=
    div_nonneg hnu (by positivity)
  have h2 : (0 : ℝ) ≤ ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) ^ 2 * z.im ^ 4) :=
    div_nonneg hnu (by positivity)
  linarith

/-- **The second moment of the row term.** It is `O(1/D²)`, which is what the weighted
Cauchy step of section 2.3 needs. -/
private theorem lintegral_normSq_rowTerm_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin (n + 1)) :
    ∫⁻ Y, ENNReal.ofReal (‖rowTerm z u Y k‖ ^ 2) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (rowTermSqBound ν z D) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have husum : ∑ l, u l ^ 2 = u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_congr rfl fun l _ => by ring
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  have hnu4 : (0 : ℝ) ≤ ∫ x, x ^ 4 ∂ν := integral_pow_four_nonneg ν
  have hsqi : ((Real.sqrt D)⁻¹ : ℝ) ^ 2 = ((D : ℝ))⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hDR.le]
  refine lintegral_noiseMatrix_le
    (ENNReal.measurable_ofReal.comp ((measurable_rowTerm z u k).norm.pow_const 2)) k fun r => ?_
  set X : Matrix (Fin (n + 1)) (Fin D) ℝ :=
    Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r with hX
  set a : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l with ha
  have hcoef : ∑ l, ‖a l‖ ^ 2 ≤ 1 / ((D : ℝ) * z.im ^ 2) := by
    have hstep : ∀ l, ‖a l‖ ^ 2
        = ((D : ℝ))⁻¹ * ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l‖ ^ 2 := by
      intro l
      rw [ha]
      simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs, hsqi]
    rw [Finset.sum_congr rfl fun l _ => hstep l, ← Finset.mul_sum]
    have h := sum_normSq_resolvC_mulVec_le (gram_isHermitian X) hz u
    calc ((D : ℝ))⁻¹ * ∑ l, ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l‖ ^ 2
        ≤ ((D : ℝ))⁻¹ * ((u ⬝ᵥ u) / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h (by positivity)
      _ ≤ ((D : ℝ))⁻¹ * (1 / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left (div_le_div_of_nonneg_right hu (by positivity))
            (by positivity)
      _ = 1 / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hcnn : (0 : ℝ) ≤ ∑ l, ‖a l‖ ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  have hterm : ∀ t : Fin D → ℝ,
      ‖rowTerm z u (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k‖ ^ 2
        = ((D : ℝ))⁻¹
          * (LinForm.linForm u t ^ 2 * ‖∑ l, a l * (t l : ℂ)‖ ^ 2) := by
    intro t
    have hasum : (∑ l, ((((Real.sqrt D)⁻¹ : ℝ) : ℂ)
        * (R4C.resolvC (gram X) z *ᵥ R4C.cvec u) l) * (t l : ℂ))
        = ∑ l, a l * (t l : ℂ) := by simp only [ha]
    simp only [rowTerm]
    rw [updateRow_insertNth k r t, ← hX, rowVec_insertNth k r t, cformC_smul_left_eq, hasum]
    have hdot : (u ⬝ᵥ ((Real.sqrt D)⁻¹ • t) : ℝ)
        = ((Real.sqrt D)⁻¹ : ℝ) * LinForm.linForm u t := by
      rw [dotProduct, LinForm.linForm, Finset.mul_sum]
      exact Finset.sum_congr rfl fun l _ => by
        simp only [Pi.smul_apply, smul_eq_mul]; ring
    rw [hdot, Complex.ofReal_mul, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
      Complex.norm_real, Real.norm_eq_abs]
    rw [mul_pow, mul_pow, sq_abs, sq_abs, hsqi]
    ring
  set g : (Fin D → ℝ) → ℝ := fun t =>
    LinForm.linForm u t ^ 4 / (2 * (D : ℝ) ^ 2) + ‖∑ l, a l * (t l : ℂ)‖ ^ 4 / 2 with hg
  have hptw : ∀ t : Fin D → ℝ,
      ‖rowTerm z u (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k‖ ^ 2 ≤ g t := by
    intro t
    rw [hterm t]
    simp only [hg]
    set A : ℝ := LinForm.linForm u t
    set B : ℝ := ‖∑ l, a l * (t l : ℂ)‖
    have hkey : A ^ 4 / (2 * (D : ℝ) ^ 2) + B ^ 4 / 2 - ((D : ℝ))⁻¹ * (A ^ 2 * B ^ 2)
        = (A ^ 2 - (D : ℝ) * B ^ 2) ^ 2 / (2 * (D : ℝ) ^ 2) := by
      field_simp
      ring
    have hnn : (0 : ℝ) ≤ (A ^ 2 - (D : ℝ) * B ^ 2) ^ 2 / (2 * (D : ℝ) ^ 2) := by positivity
    linarith
  obtain ⟨h4int, h4bd⟩ := integral_clinForm_pow_four_le hν a
  have hu4int : Integrable (fun t : Fin D → ℝ => LinForm.linForm u t ^ 4)
      (Measure.pi fun _ : Fin D => ν) :=
    LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 u
  have hgint : Integrable g (Measure.pi fun _ : Fin D => ν) := by
    rw [hg]
    exact (hu4int.div_const _).add (h4int.div_const 2)
  have hgnn : ∀ t, 0 ≤ g t := by
    intro t; rw [hg]; positivity
  have hu4bd : ∫ t, LinForm.linForm u t ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ (∫ x, x ^ 4 ∂ν) + 3 := by
    refine (LinForm.integral_linForm_pow_four_le hν.mean hν.var hν.mom4 u).trans ?_
    have h1 : ∑ l, u l ^ 4 ≤ 1 := by
      refine (sum_pow_four_le u).trans ?_
      rw [husum]
      nlinarith
    have h2 : (∑ l, u l ^ 2) ^ 2 ≤ 1 := by rw [husum]; nlinarith
    nlinarith [h1, h2, hnu4]
  have hgbd : ∫ t, g t ∂(Measure.pi fun _ : Fin D => ν) ≤ rowTermSqBound ν z D := by
    rw [hg, integral_add (hu4int.div_const _) (h4int.div_const 2), integral_div, integral_div,
      rowTermSqBound]
    have hA : (∫ t, LinForm.linForm u t ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
        / (2 * (D : ℝ) ^ 2) ≤ ((∫ x, x ^ 4 ∂ν) + 3) / (2 * (D : ℝ) ^ 2) :=
      div_le_div_of_nonneg_right hu4bd (by positivity)
    have hB : (∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)) / 2
        ≤ ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) ^ 2 * z.im ^ 4) := by
      have hstep : ∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
          ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 := by
        refine h4bd.trans ?_
        exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hcnn hcoef 2) (by linarith)
      have heq : 2 * ((∫ x, x ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 / 2
          = ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) ^ 2 * z.im ^ 4) := by
        field_simp
      calc (∫ t, ‖∑ l, a l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)) / 2
          ≤ 2 * ((∫ x, x ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 / 2 :=
            div_le_div_of_nonneg_right hstep (by norm_num)
        _ = ((∫ x, x ^ 4 ∂ν) + 3) / ((D : ℝ) ^ 2 * z.im ^ 4) := heq
    linarith
  calc ∫⁻ t, ENNReal.ofReal (‖rowTerm z u
        (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k‖ ^ 2)
        ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ∫⁻ t, ENNReal.ofReal (g t) ∂(Measure.pi fun _ : Fin D => ν) :=
        lintegral_mono fun t => ENNReal.ofReal_le_ofReal (hptw t)
    _ = ENNReal.ofReal (∫ t, g t ∂(Measure.pi fun _ : Fin D => ν)) :=
        (ofReal_integral_eq_lintegral_ofReal hgint
          (Filter.Eventually.of_forall hgnn)).symm
    _ ≤ ENNReal.ofReal (rowTermSqBound ν z D) := ENNReal.ofReal_le_ofReal hgbd

/-- The second moment of the row scalar against the deterministic root, `O(1/D)`. -/
noncomputable def alphaSqBound (ν : Measure ℝ) (z : ℂ) (P D : ℕ) : ℝ :=
  3 * (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2))
    + (3 * (1 / ((D : ℝ) * z.im)) ^ 2 + 3 * traceSqBound ν z P D)

theorem alphaSqBound_nonneg (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im) (P D : ℕ) :
    0 ≤ alphaSqBound ν z P D := by
  have hnu : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
    have := integral_pow_four_nonneg ν; linarith
  have hR : 0 ≤ traceSqBound ν z P D := by
    rw [traceSqBound]
    exact mul_nonneg (by positivity) (residConst_nonneg ν D hz)
  rw [alphaSqBound]
  have h1 : (0 : ℝ) ≤ 3 * (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2)) := by
    have : (0 : ℝ) ≤ ((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2) :=
      div_nonneg hnu (by positivity)
    linarith
  have h2 : (0 : ℝ) ≤ 3 * (1 / ((D : ℝ) * z.im)) ^ 2 := by positivity
  linarith

/-- **The `L²` bound on `α k - m`.** The row scalar is within `O(D^{-1/2})` of the MP root, in
mean square: the four-moment step of unit G3, the deterministic trace stability, and the `L²`
trace law of milestone 3. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem lintegral_normSq_alphaRow_sub_mC_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) (k : Fin (n + 1)) :
    ∫⁻ Y, ENNReal.ofReal (‖alphaRow Y z k - MP.mC (((n + 1 : ℕ) : ℝ) / D) z‖ ^ 2)
        ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (alphaSqBound ν z (n + 1) D) := by
  have hprob := hν.prob
  have hP : 0 < n + 1 := Nat.succ_pos n
  set m : ℂ := MP.mC (((n + 1 : ℕ) : ℝ) / D) z with hm
  set f1 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y =>
    ‖alphaRow Y z k - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ ^ 2 with hf1
  set f3 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y =>
    ‖R4C.stieltjesC (gram Y) z - m‖ ^ 2 with hf3
  set B2 : ℝ := 3 * (1 / ((D : ℝ) * z.im)) ^ 2 with hB2
  have hf1m : Measurable f1 :=
    (((measurable_alphaRow z k).sub
      (measurable_stieltjesC (measurable_gram_loo k) z)).norm.pow_const 2)
  have hf3m : Measurable f3 :=
    (((measurable_stieltjesC measurable_gram_self z).sub measurable_const).norm.pow_const 2)
  have hf1n : ∀ Y, 0 ≤ f1 Y := fun Y => by rw [hf1]; positivity
  have hf3n : ∀ Y, 0 ≤ f3 Y := fun Y => by rw [hf3]; positivity
  have hpt : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
      ‖alphaRow Y z k - m‖ ^ 2 ≤ 3 * f1 Y + (B2 + 3 * f3 Y) := by
    intro Y
    set A : ℝ := ‖alphaRow Y z k - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z‖ with hA
    set Bm : ℝ := ‖R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z
      - R4C.stieltjesC (gram Y) z‖ with hBm
    set C : ℝ := ‖R4C.stieltjesC (gram Y) z - m‖ with hC
    have htri : ‖alphaRow Y z k - m‖ ≤ A + Bm + C := by
      have h1 : alphaRow Y z k - m
          = (alphaRow Y z k - R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z)
            + ((R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z
                - R4C.stieltjesC (gram Y) z) + (R4C.stieltjesC (gram Y) z - m)) := by ring
      rw [h1]
      calc ‖_ + _‖ ≤ A + ‖(R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z
              - R4C.stieltjesC (gram Y) z) + (R4C.stieltjesC (gram Y) z - m)‖ :=
            norm_add_le _ _
        _ ≤ A + (Bm + C) := by
            have := norm_add_le (R4C.stieltjesC (gram (Matrix.updateRow Y k 0)) z
              - R4C.stieltjesC (gram Y) z) (R4C.stieltjesC (gram Y) z - m)
            linarith
        _ = A + Bm + C := by ring
    have hB : Bm ≤ 1 / ((D : ℝ) * z.im) := norm_stieltjesC_sub_le Y hz hD k
    have hBn : 0 ≤ Bm := norm_nonneg _
    have hAn : 0 ≤ A := norm_nonneg _
    have hCn : 0 ≤ C := norm_nonneg _
    have hnn : 0 ≤ ‖alphaRow Y z k - m‖ := norm_nonneg _
    have hsq : ‖alphaRow Y z k - m‖ ^ 2 ≤ (A + Bm + C) ^ 2 := pow_le_pow_left₀ hnn htri 2
    have hBsq : Bm ^ 2 ≤ (1 / ((D : ℝ) * z.im)) ^ 2 := pow_le_pow_left₀ hBn hB 2
    have hf1v : f1 Y = A ^ 2 := rfl
    have hf3v : f3 Y = C ^ 2 := rfl
    rw [hf1v, hf3v, hB2]
    nlinarith [hsq, hBsq, sq_nonneg (A - Bm), sq_nonneg (Bm - C), sq_nonneg (A - C)]
  have hstep1 := lintegral_normSq_alphaRow_sub_le (ν := ν) hν hz hD k
  have hstep3 := lintegral_normSq_stieltjesC_sub_mC_le (P := n + 1) (D := D) hν hz hP hD
  have hnu2 : (0 : ℝ) ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
    have := integral_pow_four_nonneg ν; linarith
  have hRnn : (0 : ℝ) ≤ traceSqBound ν z (n + 1) D := by
    rw [traceSqBound]
    exact mul_nonneg (by positivity) (residConst_nonneg ν D hz)
  calc ∫⁻ Y, ENNReal.ofReal (‖alphaRow Y z k - m‖ ^ 2) ∂(noiseMatrix ν (n + 1) D)
      ≤ ∫⁻ Y, ENNReal.ofReal (3 * f1 Y + (B2 + 3 * f3 Y)) ∂(noiseMatrix ν (n + 1) D) :=
        lintegral_mono fun Y => ENNReal.ofReal_le_ofReal (hpt Y)
    _ = (∫⁻ Y, ENNReal.ofReal (3 * f1 Y) ∂(noiseMatrix ν (n + 1) D))
          + ∫⁻ Y, ENNReal.ofReal (B2 + 3 * f3 Y) ∂(noiseMatrix ν (n + 1) D) :=
        lintegral_ofReal_add (hf1m.const_mul 3) (fun Y => by
          have := hf1n Y; linarith)
          (fun Y => by have := hf3n Y; rw [hB2]; positivity)
    _ = (ENNReal.ofReal 3 * ∫⁻ Y, ENNReal.ofReal (f1 Y) ∂(noiseMatrix ν (n + 1) D))
          + ((∫⁻ _Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ENNReal.ofReal B2
              ∂(noiseMatrix ν (n + 1) D))
            + ENNReal.ofReal 3 * ∫⁻ Y, ENNReal.ofReal (f3 Y) ∂(noiseMatrix ν (n + 1) D)) := by
        rw [lintegral_ofReal_const_mul hf1m (by norm_num),
          lintegral_ofReal_add measurable_const (fun _ => by rw [hB2]; positivity)
            (fun Y => by have := hf3n Y; linarith),
          lintegral_ofReal_const_mul hf3m (by norm_num)]
    _ ≤ (ENNReal.ofReal 3 * ENNReal.ofReal (((∫ x, x ^ 4 ∂ν) + 2) / ((D : ℝ) * z.im ^ 2)))
          + (ENNReal.ofReal B2
            + ENNReal.ofReal 3 * ENNReal.ofReal (traceSqBound ν z (n + 1) D)) := by
        have hconst : (∫⁻ _Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ENNReal.ofReal B2
            ∂(noiseMatrix ν (n + 1) D)) = ENNReal.ofReal B2 := by
          rw [lintegral_const, measure_univ, mul_one]
        rw [hconst]
        gcongr
    _ = ENNReal.ofReal (alphaSqBound ν z (n + 1) D) := by
        rw [← ENNReal.ofReal_mul (by norm_num), ← ENNReal.ofReal_mul (by norm_num),
          ← ENNReal.ofReal_add (by rw [hB2]; positivity) (by positivity),
          ← ENNReal.ofReal_add (by positivity) (by rw [hB2]; positivity), alphaSqBound, hB2]

/-- The downdate scalar against its deterministic limit.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem sub_secW_eq {c : ℝ} (hc : 0 < c) (hz : 0 < z.im)
    (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) :
    (1 - sRow Y z k) - secW c z
      = (1 - sRow Y z k) * (MP.mC c z - alphaRow Y z k) * secW c z := by
  have h1 := one_add_alphaRow_mul Y hz k
  have h2 := one_add_mC_mul_secW hc hz
  linear_combination (secW c z) * h1 - (1 - sRow Y z k) * h2

/-- The `L∞`-to-`L²` bound on that difference. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem norm_sub_secW_le {c : ℝ} (hc : 0 < c) (hz : 0 < z.im) (hd : 0 < d)
    (Y : Matrix (Fin p) (Fin d) ℝ) (k : Fin p) :
    ‖(1 - sRow Y z k) - secW c z‖
      ≤ (‖z‖ / z.im) ^ 2 * ‖alphaRow Y z k - MP.mC c z‖ := by
  rw [sub_secW_eq hc hz Y k, norm_mul, norm_mul, norm_sub_rev (MP.mC c z)]
  have h1 := norm_one_sub_sRow_le Y hz hd k
  have h2 := norm_secW_le hc hz
  have h3 : (0 : ℝ) ≤ ‖alphaRow Y z k - MP.mC c z‖ := norm_nonneg _
  have h4 : (0 : ℝ) ≤ ‖z‖ / z.im := by positivity
  calc ‖1 - sRow Y z k‖ * ‖alphaRow Y z k - MP.mC c z‖ * ‖secW c z‖
      = (‖1 - sRow Y z k‖ * ‖secW c z‖) * ‖alphaRow Y z k - MP.mC c z‖ := by ring
    _ ≤ ((‖z‖ / z.im) * (‖z‖ / z.im)) * ‖alphaRow Y z k - MP.mC c z‖ :=
        mul_le_mul_of_nonneg_right (mul_le_mul h1 h2 (norm_nonneg _) h4) h3
    _ = (‖z‖ / z.im) ^ 2 * ‖alphaRow Y z k - MP.mC c z‖ := by ring

/-- **Step 2 of section 2.3.** Replacing the random downdate scalar by the deterministic
`secW` costs `O(D^{-3/2})` per row: the weighted Cauchy step at the weight `√D`. -/
private theorem norm_integral_dev_rowTerm_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) (k : Fin (n + 1)) :
    ‖∫ Y, ((1 - sRow Y z k) - secW (((n + 1 : ℕ) : ℝ) / D) z) * rowTerm z u Y k
        ∂(noiseMatrix ν (n + 1) D)‖
      ≤ 1 / (2 * Real.sqrt D) * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D)
        + Real.sqrt D / 2 * rowTermSqBound ν z D := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hlam : (0 : ℝ) < Real.sqrt D := Real.sqrt_pos.mpr hDR
  set c : ℝ := ((n + 1 : ℕ) : ℝ) / D with hc
  have hcpos : 0 < c := by rw [hc]; positivity
  set dev : Matrix (Fin (n + 1)) (Fin D) ℝ → ℂ := fun Y =>
    (1 - sRow Y z k) - secW c z with hdev
  have hdevm : Measurable dev :=
    (measurable_const.sub
      (measurable_qformC measurable_gram_self z (measurable_rowVec k))).sub measurable_const
  set F1 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => ‖dev Y‖ ^ 2 with hF1
  set F2 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => ‖rowTerm z u Y k‖ ^ 2 with hF2
  have hF1m : Measurable F1 := hdevm.norm.pow_const 2
  have hF2m : Measurable F2 := (measurable_rowTerm z u k).norm.pow_const 2
  have hF1nn : ∀ Y, 0 ≤ F1 Y := fun Y => by rw [hF1]; positivity
  have hF2nn : ∀ Y, 0 ≤ F2 Y := fun Y => by rw [hF2]; positivity
  have hKnn : (0 : ℝ) ≤ (‖z‖ / z.im) ^ 4 := by positivity
  have hAnn : (0 : ℝ) ≤ alphaSqBound ν z (n + 1) D := alphaSqBound_nonneg ν hz (n + 1) D
  have hRnn : (0 : ℝ) ≤ rowTermSqBound ν z D := rowTermSqBound_nonneg ν z D
  -- the two second moments
  have hL1 : ∫⁻ Y, ENNReal.ofReal (F1 Y) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D) := by
    have hpt : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
        ENNReal.ofReal (F1 Y)
          ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * ‖alphaRow Y z k - MP.mC c z‖ ^ 2) := by
      intro Y
      refine ENNReal.ofReal_le_ofReal ?_
      have h := norm_sub_secW_le hcpos hz hD Y k
      have h2 := pow_le_pow_left₀ (norm_nonneg (dev Y)) h 2
      calc F1 Y ≤ ((‖z‖ / z.im) ^ 2 * ‖alphaRow Y z k - MP.mC c z‖) ^ 2 := h2
        _ = (‖z‖ / z.im) ^ 4 * ‖alphaRow Y z k - MP.mC c z‖ ^ 2 := by ring
    calc ∫⁻ Y, ENNReal.ofReal (F1 Y) ∂(noiseMatrix ν (n + 1) D)
        ≤ ∫⁻ Y, ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * ‖alphaRow Y z k - MP.mC c z‖ ^ 2)
          ∂(noiseMatrix ν (n + 1) D) := lintegral_mono hpt
      _ = ENNReal.ofReal ((‖z‖ / z.im) ^ 4)
            * ∫⁻ Y, ENNReal.ofReal (‖alphaRow Y z k - MP.mC c z‖ ^ 2)
              ∂(noiseMatrix ν (n + 1) D) :=
          lintegral_ofReal_const_mul
            (((measurable_alphaRow z k).sub measurable_const).norm.pow_const 2) hKnn
      _ ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4) * ENNReal.ofReal (alphaSqBound ν z (n + 1) D) := by
          gcongr
          exact lintegral_normSq_alphaRow_sub_mC_le hν hz hD k
      _ = ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D) :=
          (ENNReal.ofReal_mul hKnn).symm
  have hL2 := lintegral_normSq_rowTerm_le hν hz hD hu k
  have hF1int : Integrable F1 (noiseMatrix ν (n + 1) D) :=
    integrable_of_lintegral_ofReal_ne_top hF1m hF1nn
      (ne_top_of_le_ne_top ENNReal.ofReal_ne_top hL1)
  have hF2int : Integrable F2 (noiseMatrix ν (n + 1) D) :=
    integrable_of_lintegral_ofReal_ne_top hF2m hF2nn
      (ne_top_of_le_ne_top ENNReal.ofReal_ne_top hL2)
  have hF1bd : ∫ Y, F1 Y ∂(noiseMatrix ν (n + 1) D)
      ≤ (‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D :=
    integral_le_of_lintegral_le hF1int hF1nn (mul_nonneg hKnn hAnn) hL1
  have hF2bd : ∫ Y, F2 Y ∂(noiseMatrix ν (n + 1) D) ≤ rowTermSqBound ν z D :=
    integral_le_of_lintegral_le hF2int hF2nn hRnn hL2
  -- the weighted Cauchy step
  set G : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y =>
    1 / (2 * Real.sqrt D) * F1 Y + Real.sqrt D / 2 * F2 Y with hG
  have hGint : Integrable G (noiseMatrix ν (n + 1) D) :=
    (hF1int.const_mul _).add (hF2int.const_mul _)
  have hprod : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ‖dev Y * rowTerm z u Y k‖ ≤ G Y := by
    intro Y
    rw [norm_mul]
    simp only [hG, hF1, hF2]
    set x : ℝ := ‖dev Y‖
    set y : ℝ := ‖rowTerm z u Y k‖
    have hkey : 1 / (2 * Real.sqrt D) * x ^ 2 + Real.sqrt D / 2 * y ^ 2 - x * y
        = (x - Real.sqrt D * y) ^ 2 / (2 * Real.sqrt D) := by
      field_simp
      ring
    have hnn : (0 : ℝ) ≤ (x - Real.sqrt D * y) ^ 2 / (2 * Real.sqrt D) := by positivity
    linarith
  have hpint : Integrable (fun Y => dev Y * rowTerm z u Y k) (noiseMatrix ν (n + 1) D) :=
    Integrable.mono' hGint (hdevm.mul (measurable_rowTerm z u k)).aestronglyMeasurable
      (Filter.Eventually.of_forall hprod)
  calc ‖∫ Y, dev Y * rowTerm z u Y k ∂(noiseMatrix ν (n + 1) D)‖
      ≤ ∫ Y, ‖dev Y * rowTerm z u Y k‖ ∂(noiseMatrix ν (n + 1) D) :=
        norm_integral_le_integral_norm _
    _ ≤ ∫ Y, G Y ∂(noiseMatrix ν (n + 1) D) :=
        integral_mono hpint.norm hGint hprod
    _ = 1 / (2 * Real.sqrt D) * (∫ Y, F1 Y ∂(noiseMatrix ν (n + 1) D))
          + Real.sqrt D / 2 * ∫ Y, F2 Y ∂(noiseMatrix ν (n + 1) D) := by
        rw [hG, integral_add (hF1int.const_mul _) (hF2int.const_mul _), integral_const_mul,
          integral_const_mul]
    _ ≤ 1 / (2 * Real.sqrt D) * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D)
          + Real.sqrt D / 2 * rowTermSqBound ν z D := by
        gcongr

/-- The total error of the self-consistent equation, `O(P D^{-3/2})`. -/
noncomputable def meanErr (ν : Measure ℝ) (z : ℂ) (P D : ℕ) : ℝ :=
  (P : ℝ) * ((1 / (2 * Real.sqrt D) * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z P D)
      + Real.sqrt D / 2 * rowTermSqBound ν z D)
    + (‖z‖ / z.im) * (1 / (D : ℝ)) * ((‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2))))

/-- **The mean of section 2.3.** The expected quadratic form is within `O(P D^{-3/2}/η)` of
`(u ⬝ u) m`. This is the step a variance bound cannot give. -/
private theorem norm_integral_qformC_sub_mC_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) :
    ‖(∫ Y, R4C.qformC (gram Y) z u ∂(noiseMatrix ν (n + 1) D))
        - ((u ⬝ᵥ u : ℝ) : ℂ) * MP.mC (((n + 1 : ℕ) : ℝ) / D) z‖
      ≤ 1 / z.im * meanErr ν z (n + 1) D := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hlam : (0 : ℝ) < Real.sqrt D := Real.sqrt_pos.mpr hDR
  set μ : Measure (Matrix (Fin (n + 1)) (Fin D) ℝ) := noiseMatrix ν (n + 1) D with hμ
  set c : ℝ := ((n + 1 : ℕ) : ℝ) / D with hc
  have hcpos : 0 < c := by rw [hc]; positivity
  set w : ℂ := secW c z with hw
  set m : ℂ := MP.mC c z with hm
  set Λ : ℂ := ∫ Y, R4C.qformC (gram Y) z u ∂μ with hΛ
  set dev : Matrix (Fin (n + 1)) (Fin D) ℝ → Fin (n + 1) → ℂ := fun Y k =>
    (1 - sRow Y z k) - w with hdev
  -- integrability
  have hQint : Integrable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      R4C.qformC (gram Y) z u) μ := by
    have hm2 : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
        R4C.qformC (gram Y) z u :=
      measurable_qformC (y := fun _ : Matrix (Fin (n + 1)) (Fin D) ℝ => u)
        measurable_gram_self z (fun i => measurable_const)
    refine Integrable.mono' (integrable_const (μ := μ) (1 / z.im))
      hm2.aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    exact (R4C.norm_qformC_le (gram_isHermitian Y) hz u).trans
      (div_le_div_of_nonneg_right hu hz.le)
  have hrint : ∀ k, Integrable (fun Y => rowTerm z u Y k) μ := fun k =>
    integrable_rowTerm hν hz hu k
  have hdevbd : ∀ (Y : Matrix (Fin (n + 1)) (Fin D) ℝ) (k), ‖dev Y k‖ ≤ 2 * (‖z‖ / z.im) := by
    intro Y k
    refine (norm_sub_le _ _).trans ?_
    have h1 := norm_one_sub_sRow_le Y hz hD k
    have h2 := norm_secW_le hcpos hz
    rw [hw]
    linarith
  have hdevm : ∀ k, Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => dev Y k := fun k =>
    (measurable_const.sub
      (measurable_qformC measurable_gram_self z (measurable_rowVec k))).sub measurable_const
  have hdevint : ∀ k, Integrable (fun Y => dev Y k * rowTerm z u Y k) μ := by
    intro k
    refine Integrable.mono' (((hrint k).norm.const_mul (2 * (‖z‖ / z.im))))
      ((hdevm k).mul (measurable_rowTerm z u k)).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => ?_)
    rw [norm_mul]
    exact mul_le_mul_of_nonneg_right (hdevbd Y k) (norm_nonneg _)
  have hfint : ∀ k, Integrable (fun Y => (1 - sRow Y z k) * rowTerm z u Y k) μ := by
    intro k
    have hsplit : (fun Y => (1 - sRow Y z k) * rowTerm z u Y k)
        = fun Y => dev Y k * rowTerm z u Y k + w * rowTerm z u Y k := by
      funext Y; rw [hdev]; ring
    rw [hsplit]
    exact (hdevint k).add ((hrint k).const_mul w)
  -- the integrated identity
  have hLHS : ∫ Y, (z * R4C.qformC (gram Y) z u + ((u ⬝ᵥ u : ℝ) : ℂ)) ∂μ
      = z * Λ + ((u ⬝ᵥ u : ℝ) : ℂ) := by
    rw [integral_add (hQint.const_mul z) (integrable_const _), integral_const_mul, integral_const]
    have hone : μ.real Set.univ = 1 := by rw [hμ]; simp
    rw [hone, one_smul, hΛ]
  have hRHS : ∫ Y, (∑ k, (1 - sRow Y z k) * rowTerm z u Y k) ∂μ
      = ∑ k, ∫ Y, (1 - sRow Y z k) * rowTerm z u Y k ∂μ :=
    integral_finsetSum _ fun k _ => hfint k
  have hEq : z * Λ + ((u ⬝ᵥ u : ℝ) : ℂ)
      = ∑ k, ∫ Y, (1 - sRow Y z k) * rowTerm z u Y k ∂μ := by
    rw [← hLHS, ← hRHS]
    exact integral_congr_ae (Filter.Eventually.of_forall fun Y =>
      z_mul_qformC_add_eq_loo Y hz u)
  -- each row integral
  have hrow : ∀ k, ∫ Y, (1 - sRow Y z k) * rowTerm z u Y k ∂μ
      = (∫ Y, dev Y k * rowTerm z u Y k ∂μ)
        + w * (((D : ℂ))⁻¹ * (Λ + ∫ Y, (R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
            - R4C.qformC (gram Y) z u) ∂μ)) := by
    intro k
    have hsplit : (fun Y => (1 - sRow Y z k) * rowTerm z u Y k)
        = fun Y => dev Y k * rowTerm z u Y k + w * rowTerm z u Y k := by
      funext Y; rw [hdev]; ring
    rw [hsplit, integral_add (hdevint k) ((hrint k).const_mul w), integral_const_mul,
      integral_rowTerm_eq hν hz hD hu k]
    congr 2
    rw [integral_sub (integrable_qformC_loo hν hz hu k) hQint, hΛ]
    ring
  rw [Finset.sum_congr rfl fun k _ => hrow k] at hEq
  -- solve for Λ
  set E : Fin (n + 1) → ℂ := fun k => ∫ Y, dev Y k * rowTerm z u Y k ∂μ with hE
  set Δ : Fin (n + 1) → ℂ := fun k => ∫ Y, (R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
      - R4C.qformC (gram Y) z u) ∂μ with hΔ
  set R : ℂ := (∑ k, E k) + w * ((D : ℂ))⁻¹ * ∑ k, Δ k with hR
  have hcC : ((c : ℝ) : ℂ) = ((n + 1 : ℕ) : ℂ) * ((D : ℂ))⁻¹ := by
    rw [hc]
    push_cast
    ring
  have hsum : ∑ k, (E k + w * (((D : ℂ))⁻¹ * (Λ + Δ k)))
      = R + ((c : ℝ) : ℂ) * w * Λ := by
    have hexp : ∀ k, E k + w * (((D : ℂ))⁻¹ * (Λ + Δ k))
        = (E k + w * ((D : ℂ))⁻¹ * Δ k) + w * ((D : ℂ))⁻¹ * Λ := fun k => by ring
    have hA : ∑ k, (E k + w * ((D : ℂ))⁻¹ * Δ k)
        = (∑ k, E k) + w * ((D : ℂ))⁻¹ * ∑ k, Δ k := by
      rw [Finset.sum_add_distrib, Finset.mul_sum]
    have hB : (∑ _k : Fin (n + 1), w * ((D : ℂ))⁻¹ * Λ)
        = ((n + 1 : ℕ) : ℂ) * ((D : ℂ))⁻¹ * w * Λ := by
      rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
      push_cast
      ring
    rw [Finset.sum_congr rfl fun k _ => hexp k, Finset.sum_add_distrib, hA, hB, hR, hcC]
  rw [hsum] at hEq
  have heq2 : (z - ((c : ℝ) : ℂ) * w) * Λ + ((u ⬝ᵥ u : ℝ) : ℂ) = R := by
    linear_combination hEq
  have hkey : m * (z - ((c : ℝ) : ℂ) * w) = -1 := mC_mul_sub_secW hcpos hz
  have hfinal : Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m = -(m * R) := by
    linear_combination Λ * hkey - m * heq2
  -- the size of R
  have hEbd : ∀ k, ‖E k‖ ≤ 1 / (2 * Real.sqrt D) * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D)
      + Real.sqrt D / 2 * rowTermSqBound ν z D := fun k =>
    norm_integral_dev_rowTerm_le hν hz hD hu k
  have hΔbd : ∀ k, ‖Δ k‖ ≤ (‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2)) := by
    intro k
    have hint : Integrable (fun Y => R4C.qformC (gram (Matrix.updateRow Y k 0)) z u
        - R4C.qformC (gram Y) z u) μ := (integrable_qformC_loo hν hz hu k).sub hQint
    refine (norm_integral_le_integral_norm _).trans ?_
    refine integral_le_of_lintegral_le hint.norm (fun Y => norm_nonneg _) (by positivity) ?_
    exact lintegral_norm_qformC_sub_le hν hz hD hu k
  have hRbd : ‖R‖ ≤ meanErr ν z (n + 1) D := by
    have hwn : ‖w‖ ≤ ‖z‖ / z.im := norm_secW_le hcpos hz
    have hDn : ‖((D : ℂ))⁻¹‖ = 1 / (D : ℝ) := by
      rw [norm_inv, Complex.norm_natCast, one_div]
    have h1 : ‖∑ k, E k‖ ≤ (((n + 1 : ℕ) : ℝ))
        * (1 / (2 * Real.sqrt D) * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D)
          + Real.sqrt D / 2 * rowTermSqBound ν z D) := by
      refine (norm_sum_le _ _).trans ?_
      calc ∑ k, ‖E k‖ ≤ ∑ _k : Fin (n + 1), (1 / (2 * Real.sqrt D)
              * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D)
            + Real.sqrt D / 2 * rowTermSqBound ν z D) := Finset.sum_le_sum fun k _ => hEbd k
        _ = _ := by
            rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    have h2 : ‖∑ k, Δ k‖ ≤ (((n + 1 : ℕ) : ℝ)) * ((‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2))) := by
      refine (norm_sum_le _ _).trans ?_
      calc ∑ k, ‖Δ k‖ ≤ ∑ _k : Fin (n + 1), ((‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2))) :=
            Finset.sum_le_sum fun k _ => hΔbd k
        _ = _ := by rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    have hnn1 : (0 : ℝ) ≤ ‖z‖ / z.im := by positivity
    have hnn2 : (0 : ℝ) ≤ (‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2)) := by positivity
    rw [hR, meanErr]
    refine (norm_add_le _ _).trans ?_
    have h3 : ‖w * ((D : ℂ))⁻¹ * ∑ k, Δ k‖
        ≤ (‖z‖ / z.im) * (1 / (D : ℝ))
          * ((((n + 1 : ℕ) : ℝ)) * ((‖z‖ / z.im) * (1 / ((D : ℝ) * z.im ^ 2)))) := by
      rw [norm_mul, norm_mul, hDn]
      refine mul_le_mul (mul_le_mul hwn le_rfl (by positivity) hnn1) h2 (norm_nonneg _)
        (by positivity)
    have hPnn : (0 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) := by positivity
    nlinarith [h1, h3, hPnn]
  rw [hfinal, norm_neg, norm_mul]
  have hmn : ‖m‖ ≤ 1 / z.im := by
    have h := MP.norm_mC_le hcpos.le hz
    rw [← one_div] at h
    exact h
  have hRnn : (0 : ℝ) ≤ ‖R‖ := norm_nonneg _
  have hMnn : (0 : ℝ) ≤ meanErr ν z (n + 1) D := le_trans hRnn hRbd
  calc ‖m‖ * ‖R‖ ≤ (1 / z.im) * ‖R‖ := mul_le_mul_of_nonneg_right hmn hRnn
    _ ≤ 1 / z.im * meanErr ν z (n + 1) D :=
        mul_le_mul_of_nonneg_left hRbd (by positivity)

/-! ### The second moment of `isoF`, and statements 1 and 2 -/

/-- The second moment about a fixed point. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_sq_sub_const {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {g : β → ℝ} (h1 : Integrable g μ)
    (h2 : Integrable (fun t => g t ^ 2) μ) (b : ℝ) :
    ∫ t, (g t - b) ^ 2 ∂μ = (∫ t, g t ^ 2 ∂μ) - 2 * b * (∫ t, g t ∂μ) + b ^ 2 := by
  have hone : μ.real Set.univ = 1 := by simp
  have hI2 : Integrable (fun t => -(2 * b) * g t) μ := h1.const_mul _
  have hI3 : Integrable (fun _t : β => (b : ℝ) ^ 2) μ := integrable_const _
  have hI1 : Integrable (fun t => -(2 * b) * g t + b ^ 2) μ := hI2.add hI3
  have hfun : (fun t => (g t - b) ^ 2)
      = fun t => g t ^ 2 + (-(2 * b) * g t + b ^ 2) := by funext t; ring
  rw [hfun, integral_add h2 hI1, integral_add hI2 hI3, integral_const_mul, integral_const,
    hone, one_smul]
  ring

/-- Jensen on a probability measure. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem sq_integral_le {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {g : β → ℝ} (h1 : Integrable g μ)
    (h2 : Integrable (fun t => g t ^ 2) μ) :
    (∫ t, g t ∂μ) ^ 2 ≤ ∫ t, g t ^ 2 ∂μ := by
  have h := integral_sq_sub_const h1 h2 (∫ t, g t ∂μ)
  have hnn : (0 : ℝ) ≤ ∫ t, (g t - ∫ s, g s ∂μ) ^ 2 ∂μ :=
    integral_nonneg (μ := μ) fun t => sq_nonneg _
  linarith

/-- The second moment splits into the variance and the squared mean.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem integral_sq_eq_var_add {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {g : β → ℝ} (h1 : Integrable g μ)
    (h2 : Integrable (fun t => g t ^ 2) μ) :
    ∫ t, g t ^ 2 ∂μ
      = (∫ t, (g t - ∫ s, g s ∂μ) ^ 2 ∂μ) + (∫ t, g t ∂μ) ^ 2 := by
  rw [integral_sq_sub_const h1 h2 (∫ t, g t ∂μ)]
  ring

/-- The second-moment bound on `isoF`. It mentions `ν`, `z`, `p` and `d` and nothing else:
no direction and no `ε`. Every component is `O(1/d)` at a fixed ratio `p/d`. -/
noncomputable def isoBase (ν : Measure ℝ) (z : ℂ) (p d : ℕ) : ℝ :=
  2 * (p : ℝ) * rowVarBound ν z d
    + 2 * (1 / z.im * meanErr ν z p d) ^ 2 + 2 * traceSqBound ν z p d

/-- **The rate of the isotropic bound.** Eight times `isoBase`, so that the polarization of
statement 2 carries the same rate as statement 1. It is free of the two directions and of
`ε`. -/
noncomputable def isoRate (ν : Measure ℝ) (z : ℂ) (p d : ℕ) : ℝ := 8 * isoBase ν z p d

/-- **The second moment of the isotropic form.** The variance of milestone 4, the mean of
milestone 5 and the `L²` trace law of milestone 3, put together. -/
theorem integral_normSq_isoF_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    (hP : 0 < P) (hD : 0 < D) {u : Fin D → ℝ} (hu : u ⬝ᵥ u ≤ 1) :
    ∫ Y, ‖isoF z u Y‖ ^ 2 ∂(noiseMatrix ν P D) ≤ isoBase ν z P D := by
  obtain ⟨n, rfl⟩ : ∃ n, P = n + 1 := ⟨P - 1, (Nat.succ_pred_eq_of_pos hP).symm⟩
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  set μ : Measure (Matrix (Fin (n + 1)) (Fin D) ℝ) := noiseMatrix ν (n + 1) D with hμ
  set c : ℝ := ((n + 1 : ℕ) : ℝ) / D with hc
  have hcpos : 0 < c := by rw [hc]; positivity
  set m : ℂ := MP.mC c z with hm
  have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  -- integrability of `isoF` and of its two components
  have hFm : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => isoF z u Y :=
    measurable_isoF z u
  have hFbd : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ‖isoF z u Y‖ ≤ 2 / z.im := by
    intro Y
    refine (norm_isoF_le hz hD u Y).trans ?_
    have h1 : 2 * (u ⬝ᵥ u) ≤ 2 := by linarith
    exact div_le_div_of_nonneg_right h1 hz.le
  have hFint : Integrable (fun Y => isoF z u Y) μ :=
    Integrable.mono' (integrable_const (μ := μ) (2 / z.im)) hFm.aestronglyMeasurable
      (Filter.Eventually.of_forall hFbd)
  have hcomp : ∀ (φ : ℂ → ℝ), Measurable φ → (∀ w : ℂ, |φ w| ≤ ‖w‖) →
      Integrable (fun Y => φ (isoF z u Y)) μ
        ∧ Integrable (fun Y => φ (isoF z u Y) ^ 2) μ := by
    intro φ hφm hφn
    have hm2 : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => φ (isoF z u Y) :=
      hφm.comp hFm
    have hbd : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, |φ (isoF z u Y)| ≤ 2 / z.im := fun Y =>
      (hφn _).trans (hFbd Y)
    refine ⟨Integrable.mono' (integrable_const (μ := μ) (2 / z.im)) hm2.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by rw [Real.norm_eq_abs]; exact hbd Y), ?_⟩
    refine Integrable.mono' (integrable_const (μ := μ) ((2 / z.im) ^ 2))
      (hm2.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbd Y) 2
  obtain ⟨hReI, hReS⟩ := hcomp Complex.re Complex.measurable_re Complex.abs_re_le_norm
  obtain ⟨hImI, hImS⟩ := hcomp Complex.im Complex.measurable_im Complex.abs_im_le_norm
  -- split the squared norm
  have hsplit : ∫ Y, ‖isoF z u Y‖ ^ 2 ∂μ
      = (∫ Y, (isoF z u Y).re ^ 2 ∂μ) + ∫ Y, (isoF z u Y).im ^ 2 ∂μ := by
    have hfun : (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => ‖isoF z u Y‖ ^ 2)
        = fun Y => (isoF z u Y).re ^ 2 + (isoF z u Y).im ^ 2 := by
      funext Y
      rw [Complex.sq_norm, Complex.normSq_apply]; ring
    rw [hfun, integral_add hReS hImS]
  -- the two variances
  have hvRe := variance_component_le (P := n + 1) hν hz hD hu Complex.re Complex.measurable_re
    (fun w v => (Complex.sub_re w v).symm) Complex.abs_re_le_norm
  have hvIm := variance_component_le (P := n + 1) hν hz hD hu Complex.im Complex.measurable_im
    (fun w v => (Complex.sub_im w v).symm) Complex.abs_im_le_norm
  -- the mean
  have hmeanEq : (∫ Y, (isoF z u Y).re ∂μ) ^ 2 + (∫ Y, (isoF z u Y).im ∂μ) ^ 2
      = ‖∫ Y, isoF z u Y ∂μ‖ ^ 2 := by
    have hre : (∫ Y, isoF z u Y ∂μ).re = ∫ Y, (isoF z u Y).re ∂μ := (integral_re hFint).symm
    have him : (∫ Y, isoF z u Y ∂μ).im = ∫ Y, (isoF z u Y).im ∂μ := (integral_im hFint).symm
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
    ring
  rw [← hμ] at hvRe hvIm
  -- the mean of `isoF`
  have hQint : Integrable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      R4C.qformC (gram Y) z u) μ := by
    have hm2 : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
        R4C.qformC (gram Y) z u :=
      measurable_qformC (y := fun _ : Matrix (Fin (n + 1)) (Fin D) ℝ => u)
        measurable_gram_self z (fun i => measurable_const)
    refine Integrable.mono' (integrable_const (μ := μ) (1 / z.im))
      hm2.aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    exact (R4C.norm_qformC_le (gram_isHermitian Y) hz u).trans
      (div_le_div_of_nonneg_right hu hz.le)
  have hSint : Integrable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
      R4C.stieltjesC (gram Y) z) μ := by
    refine Integrable.mono' (integrable_const (μ := μ) (1 / z.im))
      (measurable_stieltjesC measurable_gram_self z).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => ?_)
    exact R4C.norm_stieltjesC_le (gram_isHermitian Y) hz hD
  set S : ℂ := ∫ Y, R4C.stieltjesC (gram Y) z ∂μ with hS
  set Λ : ℂ := ∫ Y, R4C.qformC (gram Y) z u ∂μ with hΛ
  have hone : μ.real Set.univ = 1 := by rw [hμ]; simp
  have hisoI : ∫ Y, isoF z u Y ∂μ = Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * S := by
    simp only [isoF]
    rw [integral_sub hQint (hSint.const_mul _), integral_const_mul, hΛ, hS]
  set g : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y =>
    ‖R4C.stieltjesC (gram Y) z - m‖ with hgdef
  have hgm : Measurable g := ((measurable_stieltjesC measurable_gram_self z).sub
    measurable_const).norm
  have hgbd : ∀ Y, g Y ≤ 2 / z.im := by
    intro Y
    rw [hgdef]
    refine (norm_sub_le _ _).trans ?_
    have h1 := R4C.norm_stieltjesC_le (gram_isHermitian Y) hz hD
    have h2 : ‖m‖ ≤ 1 / z.im := by
      have h := MP.norm_mC_le hcpos.le hz
      rw [← one_div] at h
      exact h
    have h3 : (2 : ℝ) / z.im = 1 / z.im + 1 / z.im := by ring
    rw [h3]
    linarith
  have hgint : Integrable g μ :=
    Integrable.mono' (integrable_const (μ := μ) (2 / z.im)) hgm.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by
        rw [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]; exact hgbd Y)
  have hg2int : Integrable (fun Y => g Y ^ 2) μ := by
    refine Integrable.mono' (integrable_const (μ := μ) ((2 / z.im) ^ 2))
      (hgm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact pow_le_pow_left₀ (norm_nonneg _) (hgbd Y) 2
  have hg2bd : ∫ Y, g Y ^ 2 ∂μ ≤ traceSqBound ν z (n + 1) D := by
    refine integral_le_of_lintegral_le hg2int (fun Y => sq_nonneg _) ?_ ?_
    · rw [traceSqBound]
      exact mul_nonneg (by positivity) (residConst_nonneg ν D hz)
    · rw [hμ]
      exact lintegral_normSq_stieltjesC_sub_mC_le hν hz (Nat.succ_pos n) hD
  have hSm : (∫ Y, (R4C.stieltjesC (gram Y) z - m) ∂μ) = S - m := by
    rw [integral_sub hSint (integrable_const _), integral_const, hone, one_smul, hS]
  have hSmbd : ‖S - m‖ ^ 2 ≤ traceSqBound ν z (n + 1) D := by
    have h1 : ‖S - m‖ ≤ ∫ Y, g Y ∂μ := by
      rw [← hSm]
      exact norm_integral_le_integral_norm _
    have h2 : (∫ Y, g Y ∂μ) ^ 2 ≤ ∫ Y, g Y ^ 2 ∂μ := sq_integral_le hgint hg2int
    have h3 : (0 : ℝ) ≤ ‖S - m‖ := norm_nonneg _
    nlinarith
  have hmeanbd : ‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ ≤ 1 / z.im * meanErr ν z (n + 1) D := by
    rw [hΛ, hm, hμ]
    exact norm_integral_qformC_sub_mC_le hν hz hD hu
  have hMnn : (0 : ℝ) ≤ 1 / z.im * meanErr ν z (n + 1) D :=
    le_trans (norm_nonneg _) hmeanbd
  have hfinal : ‖∫ Y, isoF z u Y ∂μ‖ ^ 2
      ≤ 2 * (1 / z.im * meanErr ν z (n + 1) D) ^ 2
        + 2 * traceSqBound ν z (n + 1) D := by
    have hdecomp : ∫ Y, isoF z u Y ∂μ
        = (Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m) - ((u ⬝ᵥ u : ℝ) : ℂ) * (S - m) := by
      rw [hisoI]; ring
    have hb2 : ‖((u ⬝ᵥ u : ℝ) : ℂ) * (S - m)‖ ≤ ‖S - m‖ := by
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg huu]
      nlinarith [norm_nonneg (S - m)]
    have htri : ‖∫ Y, isoF z u Y ∂μ‖
        ≤ ‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ + ‖S - m‖ := by
      rw [hdecomp]
      exact (norm_sub_le _ _).trans (by linarith)
    have hnn1 : (0 : ℝ) ≤ ‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ := norm_nonneg _
    have hnn2 : (0 : ℝ) ≤ ‖S - m‖ := norm_nonneg _
    have hsq : ‖∫ Y, isoF z u Y ∂μ‖ ^ 2
        ≤ (‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ + ‖S - m‖) ^ 2 :=
      pow_le_pow_left₀ (norm_nonneg _) htri 2
    have hsq1 : ‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ ^ 2
        ≤ (1 / z.im * meanErr ν z (n + 1) D) ^ 2 :=
      pow_le_pow_left₀ hnn1 hmeanbd 2
    nlinarith [hsq, hsq1, hSmbd, sq_nonneg (‖Λ - ((u ⬝ᵥ u : ℝ) : ℂ) * m‖ - ‖S - m‖)]
  rw [hsplit, integral_sq_eq_var_add hReI hReS, integral_sq_eq_var_add hImI hImS, isoBase]
  linarith [hvRe, hvIm, hfinal, hmeanEq]

theorem isoBase_nonneg (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im) (p d : ℕ) :
    0 ≤ isoBase ν z p d := by
  have h1 : (0 : ℝ) ≤ 2 * (p : ℝ) * rowVarBound ν z d :=
    mul_nonneg (by positivity) (rowVarBound_nonneg ν d)
  have h2 : (0 : ℝ) ≤ 2 * (1 / z.im * meanErr ν z p d) ^ 2 := by positivity
  have h3 : (0 : ℝ) ≤ 2 * traceSqBound ν z p d := by
    have h4 : (0 : ℝ) ≤ traceSqBound ν z p d := by
      rw [traceSqBound]
      exact mul_nonneg (by positivity) (residConst_nonneg ν d hz)
    linarith
  rw [isoBase]
  linarith

/-- **Chebyshev on `isoF`.** -/
theorem measure_isoF_ge_le {ν : Measure ℝ} (hν : NoiseLaw ν) {z : ℂ} (hz : 0 < z.im)
    {p d : ℕ} (hp : 0 < p) (hd : 0 < d) (u : Fin d → ℝ) (hu : u ⬝ᵥ u ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    noiseMatrix ν p d {Y | ε ≤ ‖isoF z u Y‖} ≤ ENNReal.ofReal (isoBase ν z p d / ε ^ 2) := by
  have hprob := hν.prob
  have hFm : Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => ‖isoF z u Y‖ :=
    (measurable_isoF z u).norm
  have hFbd : ∀ Y : Matrix (Fin p) (Fin d) ℝ, ‖isoF z u Y‖ ≤ 2 / z.im := by
    intro Y
    refine (norm_isoF_le hz hd u Y).trans ?_
    have huu : (0 : ℝ) ≤ u ⬝ᵥ u := by
      rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
    have h1 : 2 * (u ⬝ᵥ u) ≤ 2 := by linarith
    exact div_le_div_of_nonneg_right h1 hz.le
  have hint : Integrable (fun Y : Matrix (Fin p) (Fin d) ℝ => ‖isoF z u Y‖ ^ 2)
      (noiseMatrix ν p d) := by
    refine Integrable.mono' (integrable_const (μ := noiseMatrix ν p d) ((2 / z.im) ^ 2))
      (hFm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact pow_le_pow_left₀ (norm_nonneg _) (hFbd Y) 2
  exact Cheb.meas_ge_le_of_integral_sq _ hFm hint hε
    (integral_normSq_isoF_le hν hz hp hd hu)

/-- **Statement 1: the quadratic form.** For every direction `u` in the unit ball, the
quadratic form of the resolvent is within `ε` of `(u ⬝ u) d⁻¹ tr G` outside a set of measure at
most `isoRate ν z p d / ε²`. The rate is free of `u` and of `ε`. -/
theorem measure_qformC_sub_ge_le {ν : Measure ℝ} (hν : NoiseLaw ν) {z : ℂ} (hz : 0 < z.im)
    {p d : ℕ} (hp : 0 < p) (hd : 0 < d) (u : Fin d → ℝ) (hu : u ⬝ᵥ u ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    noiseMatrix ν p d {Y | ε ≤ ‖R4C.qformC (gram Y) z u
        - ((u ⬝ᵥ u : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z‖}
      ≤ ENNReal.ofReal (isoRate ν z p d / ε ^ 2) := by
  have hset : {Y : Matrix (Fin p) (Fin d) ℝ | ε ≤ ‖R4C.qformC (gram Y) z u
      - ((u ⬝ᵥ u : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z‖}
      = {Y | ε ≤ ‖isoF z u Y‖} := rfl
  rw [hset]
  refine (measure_isoF_ge_le hν hz hp hd u hu hε).trans (ENNReal.ofReal_le_ofReal ?_)
  have h0 : 0 ≤ isoBase ν z p d := isoBase_nonneg ν hz p d
  rw [isoRate]
  refine div_le_div_of_nonneg_right ?_ (by positivity)
  linarith

/-! ### Statement 2: the bilinear form, by polarization -/

/-- The form is additive in its left argument. Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem cformC_add_left {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ)
    (a b c : Fin D → ℝ) :
    R4C.cformC W z (a + b) c = R4C.cformC W z a c + R4C.cformC W z b c := by
  simp only [R4C.cformC, R4C.cvec, dotProduct, Pi.add_apply, Complex.ofReal_add]
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

private theorem cformC_sub_left {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ)
    (a b c : Fin D → ℝ) :
    R4C.cformC W z (a - b) c = R4C.cformC W z a c - R4C.cformC W z b c := by
  simp only [R4C.cformC, R4C.cvec, dotProduct, Pi.sub_apply, Complex.ofReal_sub]
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- The form is homogeneous in its left argument.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem cformC_smulR_left {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (r : ℝ)
    (a b : Fin D → ℝ) :
    R4C.cformC W z (r • a) b = (r : ℂ) * R4C.cformC W z a b := by
  simp only [R4C.cformC, R4C.cvec, dotProduct, Pi.smul_apply, smul_eq_mul, Complex.ofReal_mul]
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **Polarization.** The bilinear form is a difference of two quadratic forms at directions of
squared length at most one. -/
theorem isoF_polarization {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im)
    (x y : Fin d → ℝ) :
    isoF z ((2⁻¹ : ℝ) • (x + y)) Y - isoF z ((2⁻¹ : ℝ) • (x - y)) Y
      = R4C.cformC (gram Y) z x y - ((x ⬝ᵥ y : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z := by
  have hcomm : ∀ a b : Fin d → ℝ, R4C.cformC (gram Y) z a b = R4C.cformC (gram Y) z b a :=
    fun a b => cformC_comm (gram_isHermitian Y) hz.ne' a b
  have hAR : ∀ a b c : Fin d → ℝ, R4C.cformC (gram Y) z a (b + c)
      = R4C.cformC (gram Y) z a b + R4C.cformC (gram Y) z a c := by
    intro a b c
    rw [hcomm a (b + c), cformC_add_left, hcomm b a, hcomm c a]
  have hSR : ∀ a b c : Fin d → ℝ, R4C.cformC (gram Y) z a (b - c)
      = R4C.cformC (gram Y) z a b - R4C.cformC (gram Y) z a c := by
    intro a b c
    rw [hcomm a (b - c), cformC_sub_left, hcomm b a, hcomm c a]
  have hMR : ∀ (r : ℝ) (a b : Fin d → ℝ), R4C.cformC (gram Y) z a (r • b)
      = (r : ℂ) * R4C.cformC (gram Y) z a b := by
    intro r a b
    rw [hcomm a (r • b), cformC_smulR_left, hcomm b a]
  have hq : ∀ (r : ℝ) (w : Fin d → ℝ), R4C.qformC (gram Y) z (r • w)
      = (r : ℂ) * (r : ℂ) * R4C.cformC (gram Y) z w w := by
    intro r w
    change R4C.cformC (gram Y) z (r • w) (r • w) = _
    rw [cformC_smulR_left, hMR]
    ring
  have hqp : R4C.qformC (gram Y) z ((2⁻¹ : ℝ) • (x + y))
      = ((2⁻¹ : ℝ) : ℂ) * ((2⁻¹ : ℝ) : ℂ)
        * (R4C.cformC (gram Y) z x x + R4C.cformC (gram Y) z x y
          + (R4C.cformC (gram Y) z y x + R4C.cformC (gram Y) z y y)) := by
    rw [hq, cformC_add_left, hAR, hAR]
  have hqm : R4C.qformC (gram Y) z ((2⁻¹ : ℝ) • (x - y))
      = ((2⁻¹ : ℝ) : ℂ) * ((2⁻¹ : ℝ) : ℂ)
        * (R4C.cformC (gram Y) z x x - R4C.cformC (gram Y) z x y
          - (R4C.cformC (gram Y) z y x - R4C.cformC (gram Y) z y y)) := by
    rw [hq, cformC_sub_left, hSR, hSR]
  have hdotp : ((2⁻¹ : ℝ) • (x + y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x + y))
      - ((2⁻¹ : ℝ) • (x - y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x - y)) = x ⬝ᵥ y := by
    simp only [dotProduct, Pi.smul_apply, Pi.add_apply, Pi.sub_apply, smul_eq_mul]
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  simp only [isoF]
  rw [hqp, hqm, hcomm y x]
  have hcast : ((2⁻¹ : ℝ) : ℂ) = 2⁻¹ := by norm_num
  have hdotC : (((2⁻¹ : ℝ) • (x + y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x + y)) : ℝ)
      - (((2⁻¹ : ℝ) • (x - y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x - y)) : ℝ) = x ⬝ᵥ y := hdotp
  have hdotCC : ((((2⁻¹ : ℝ) • (x + y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x + y)) : ℝ) : ℂ)
      - ((((2⁻¹ : ℝ) • (x - y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x - y)) : ℝ) : ℂ) = ((x ⬝ᵥ y : ℝ) : ℂ) := by
    rw [← Complex.ofReal_sub, hdotC]
  rw [hcast]
  linear_combination (-(R4C.stieltjesC (gram Y) z)) * hdotCC

private theorem dotProduct_add_self {D : ℕ} (a b : Fin D → ℝ) :
    (a + b) ⬝ᵥ (a + b) = a ⬝ᵥ a + 2 * (a ⬝ᵥ b) + b ⬝ᵥ b := by
  simp only [dotProduct, Pi.add_apply]
  rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

private theorem dotProduct_sub_self {D : ℕ} (a b : Fin D → ℝ) :
    (a - b) ⬝ᵥ (a - b) = a ⬝ᵥ a - 2 * (a ⬝ᵥ b) + b ⬝ᵥ b := by
  simp only [dotProduct, Pi.sub_apply]
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

private theorem dotProduct_smul_self {D : ℕ} (r : ℝ) (a : Fin D → ℝ) :
    (r • a) ⬝ᵥ (r • a) = r * r * (a ⬝ᵥ a) := by
  simp only [dotProduct, Pi.smul_apply, smul_eq_mul]
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **Statement 2: the bilinear form.** The same bound at two directions in the unit ball, by
polarization; unit G6 and unit G7 consume this form. -/
theorem measure_cformC_sub_ge_le {ν : Measure ℝ} (hν : NoiseLaw ν) {z : ℂ} (hz : 0 < z.im)
    {p d : ℕ} (hp : 0 < p) (hd : 0 < d) (x y : Fin d → ℝ) (hx : x ⬝ᵥ x ≤ 1) (hy : y ⬝ᵥ y ≤ 1)
    {ε : ℝ} (hε : 0 < ε) :
    noiseMatrix ν p d {Y | ε ≤ ‖R4C.cformC (gram Y) z x y
        - ((x ⬝ᵥ y : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z‖}
      ≤ ENNReal.ofReal (isoRate ν z p d / ε ^ 2) := by
  have hnnx : (0 : ℝ) ≤ x ⬝ᵥ x := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hnny : (0 : ℝ) ≤ y ⬝ᵥ y := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hnnp : (0 : ℝ) ≤ (x + y) ⬝ᵥ (x + y) := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hnnm : (0 : ℝ) ≤ (x - y) ⬝ᵥ (x - y) := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rw [dotProduct_add_self] at hnnp
  rw [dotProduct_sub_self] at hnnm
  have hvp : ((2⁻¹ : ℝ) • (x + y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x + y)) ≤ 1 := by
    rw [dotProduct_smul_self, dotProduct_add_self]
    nlinarith
  have hvm : ((2⁻¹ : ℝ) • (x - y)) ⬝ᵥ ((2⁻¹ : ℝ) • (x - y)) ≤ 1 := by
    rw [dotProduct_smul_self, dotProduct_sub_self]
    nlinarith
  have hhalf : (0 : ℝ) < ε / 2 := by positivity
  have hsub : {Y : Matrix (Fin p) (Fin d) ℝ | ε ≤ ‖R4C.cformC (gram Y) z x y
        - ((x ⬝ᵥ y : ℝ) : ℂ) * R4C.stieltjesC (gram Y) z‖}
      ⊆ {Y | ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x + y)) Y‖}
        ∪ {Y | ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x - y)) Y‖} := by
    intro Y hY
    by_cases h1 : ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x + y)) Y‖
    · exact Or.inl h1
    by_cases h2 : ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x - y)) Y‖
    · exact Or.inr h2
    exfalso
    rw [not_le] at h1 h2
    have hpol := isoF_polarization Y hz x y
    have hlt : ‖isoF z ((2⁻¹ : ℝ) • (x + y)) Y - isoF z ((2⁻¹ : ℝ) • (x - y)) Y‖ < ε := by
      refine lt_of_le_of_lt (norm_sub_le _ _) ?_
      linarith
    rw [hpol] at hlt
    exact absurd hY (not_le.mpr hlt)
  refine (measure_mono hsub).trans ?_
  refine (measure_union_le _ _).trans ?_
  have hb1 := measure_isoF_ge_le hν hz hp hd _ hvp hhalf
  have hb2 := measure_isoF_ge_le hν hz hp hd _ hvm hhalf
  have h0 : (0 : ℝ) ≤ isoBase ν z p d := isoBase_nonneg ν hz p d
  calc noiseMatrix ν p d {Y | ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x + y)) Y‖}
        + noiseMatrix ν p d {Y | ε / 2 ≤ ‖isoF z ((2⁻¹ : ℝ) • (x - y)) Y‖}
      ≤ ENNReal.ofReal (isoBase ν z p d / (ε / 2) ^ 2)
          + ENNReal.ofReal (isoBase ν z p d / (ε / 2) ^ 2) := add_le_add hb1 hb2
    _ = ENNReal.ofReal (isoRate ν z p d / ε ^ 2) := by
        rw [← ENNReal.ofReal_add (by positivity) (by positivity), isoRate]
        congr 1
        field_simp
        ring

/-! ### The rate vanishes

Unit G6 and unit G7 need one statement that the rate goes to zero, so that they do not unfold
the five component definitions. -/

/-- The factor of `traceSqBound` that depends on the ratio only.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
noncomputable def traceCoef (z : ℂ) (r : ℝ) : ℝ :=
  (r * (‖z‖ / z.im)) ^ 2 * ((1 / (‖z‖ * (z.im * rootLb r z ^ 2))) ^ 2 + 16 / z.im ^ 2)

/-- `traceSqBound` unfolds to `traceCoef` times `residConst`.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem traceSqBound_eq (ν : Measure ℝ) (z : ℂ) (P D : ℕ) :
    traceSqBound ν z P D = traceCoef z ((P : ℝ) / D) * residConst ν D z := rfl

/-- `traceCoef` is continuous in the ratio: `rootLb` is.
Public since F36 (2026-09-09) for `IsoMixed.lean`. -/
theorem tendsto_traceCoef {rN : ℕ → ℝ} {cc : ℝ} {z : ℂ} (hz : 0 < z.im)
    (hr : Tendsto rN atTop (𝓝 cc)) :
    Tendsto (fun N => traceCoef z (rN N)) atTop (𝓝 (traceCoef z cc)) := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hcoe : Tendsto (fun N => ((rN N : ℝ) : ℂ)) atTop (𝓝 ((cc : ℝ) : ℂ)) :=
    (Complex.continuous_ofReal.tendsto cc).comp hr
  have hnorm : Tendsto (fun N => ‖z + 1 - ((rN N : ℝ) : ℂ)‖) atTop
      (𝓝 ‖z + 1 - ((cc : ℝ) : ℂ)‖) := (tendsto_const_nhds.sub hcoe).norm
  have hden : (0 : ℝ) < 2 * (‖z‖ + ‖z + 1 - ((cc : ℝ) : ℂ)‖) := by
    have h := norm_nonneg (z + 1 - ((cc : ℝ) : ℂ))
    linarith
  have hinner : Tendsto (fun N => 1 / (2 * (‖z‖ + ‖z + 1 - ((rN N : ℝ) : ℂ)‖))) atTop
      (𝓝 (1 / (2 * (‖z‖ + ‖z + 1 - ((cc : ℝ) : ℂ)‖)))) :=
    tendsto_const_nhds.div (tendsto_const_nhds.mul (tendsto_const_nhds.add hnorm)) hden.ne'
  have hroot : Tendsto (fun N => rootLb (rN N) z) atTop (𝓝 (rootLb cc z)) := by
    simp only [rootLb]
    exact Tendsto.min tendsto_const_nhds hinner
  have hrpos : 0 < rootLb cc z := rootLb_pos hz0
  have hden2 : ‖z‖ * (z.im * rootLb cc z ^ 2) ≠ 0 :=
    (mul_pos hzn (mul_pos hz (pow_pos hrpos 2))).ne'
  simp only [traceCoef]
  refine Tendsto.mul ((hr.mul tendsto_const_nhds).pow 2) (Tendsto.add ?_ tendsto_const_nhds)
  exact Tendsto.pow (tendsto_const_nhds.div
    (tendsto_const_nhds.mul (tendsto_const_nhds.mul (hroot.pow 2))) hden2) 2

/-- **The rate vanishes.** At a ratio `p/d` with a limit, `isoRate` is `O(1/d)`, so it goes to
zero. The shape of the proof copies `tendsto_residConst` (`RMT/General/Trace.lean:905`). -/
theorem tendsto_isoRate {pN dN : ℕ → ℕ} (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im) {c : ℝ}
    (hd : Tendsto dN atTop atTop)
    (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c)) :
    Tendsto (fun N => isoRate ν z (pN N) (dN N)) atTop (𝓝 0) := by
  have hzi : z.im ≠ 0 := hz.ne'
  have hev : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  have ht : Tendsto (fun N => ((dN N : ℝ))⁻¹) atTop (𝓝 0) :=
    (tendsto_natCast_atTop_atTop.comp hd).inv_tendsto_atTop
  have hs : Tendsto (fun N => (Real.sqrt (dN N))⁻¹) atTop (𝓝 0) :=
    (Real.tendsto_sqrt_atTop.comp (tendsto_natCast_atTop_atTop.comp hd)).inv_tendsto_atTop
  have hsmall : ∀ (a : ℝ), Tendsto (fun N => a / ((dN N : ℝ) * z.im ^ 2)) atTop (𝓝 0) := by
    intro a
    have h := ht.const_mul (a / z.im ^ 2)
    rw [mul_zero] at h
    refine h.congr fun N => ?_
    rcases eq_or_ne ((dN N : ℝ)) 0 with h0 | h0
    · rw [h0]; simp
    · field_simp
  -- `D * residConst` converges
  have hresD : Tendsto (fun N => ((dN N : ℝ)) * residConst ν (dN N) z) atTop
      (𝓝 (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2)) := by
    have hnice : Tendsto (fun N => 2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
        + 2 / ((dN N : ℝ) * z.im ^ 2)) atTop
        (𝓝 (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2)) := by
      have h : Tendsto (fun N : ℕ => 2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
          + 2 / ((dN N : ℝ) * z.im ^ 2)) atTop
          (𝓝 (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2 + 0)) :=
        tendsto_const_nhds.add (hsmall 2)
      rwa [add_zero] at h
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [residConst]
    field_simp
  -- the traces
  have hcoefT : Tendsto (fun N => traceCoef z ((pN N : ℝ) / dN N)) atTop (𝓝 (traceCoef z c)) :=
    tendsto_traceCoef hz hcN
  have htrace : Tendsto (fun N => traceSqBound ν z (pN N) (dN N)) atTop (𝓝 0) := by
    have h := hcoefT.mul (tendsto_residConst ν hz hd)
    rw [mul_zero] at h
    exact h.congr fun N => (traceSqBound_eq ν z (pN N) (dN N)).symm
  have htraceD : Tendsto (fun N => ((dN N : ℝ)) * traceSqBound ν z (pN N) (dN N)) atTop
      (𝓝 (traceCoef z c * (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2))) := by
    refine (hcoefT.mul hresD).congr fun N => ?_
    rw [traceSqBound_eq]
    ring
  -- `D * alphaSqBound` converges
  have halphaD : Tendsto (fun N => ((dN N : ℝ)) * alphaSqBound ν z (pN N) (dN N)) atTop
      (𝓝 (3 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
        + 3 * (traceCoef z c * (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2)))) := by
    have hnice : Tendsto (fun N => (3 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
          + 3 / ((dN N : ℝ) * z.im ^ 2))
        + 3 * (((dN N : ℝ)) * traceSqBound ν z (pN N) (dN N))) atTop
        (𝓝 (3 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
          + 3 * (traceCoef z c * (2 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2)))) := by
      have h1 : Tendsto (fun N : ℕ => 3 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2
          + 3 / ((dN N : ℝ) * z.im ^ 2)) atTop
          (𝓝 (3 * ((∫ x, x ^ 4 ∂ν) + 2) / z.im ^ 2 + 0)) :=
        tendsto_const_nhds.add (hsmall 3)
      rw [add_zero] at h1
      exact h1.add (htraceD.const_mul 3)
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    set T : ℝ := traceSqBound ν z (pN N) (dN N) with hT
    rw [alphaSqBound, ← hT]
    field_simp
    ring
  -- `D² rowTermSqBound` and `D² rowVarBound` are constants
  have hrtD : Tendsto (fun N => ((dN N : ℝ)) ^ 2 * rowTermSqBound ν z (dN N)) atTop
      (𝓝 (((∫ x, x ^ 4 ∂ν) + 3) / 2 + ((∫ x, x ^ 4 ∂ν) + 3) / z.im ^ 4)) := by
    refine Tendsto.congr' ?_ tendsto_const_nhds
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [rowTermSqBound]
    field_simp
  have hrvD : Tendsto (fun N => ((dN N : ℝ)) ^ 2 * rowVarBound ν z (dN N)) atTop
      (𝓝 (4 * (‖z‖ / z.im) ^ 2 * ((∫ x, x ^ 4 ∂ν) + 3) / z.im ^ 4 + 2 / z.im ^ 2)) := by
    refine Tendsto.congr' ?_ tendsto_const_nhds
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [rowVarBound]
    field_simp
  -- the variance term
  have hPrv : Tendsto (fun N => (pN N : ℝ) * rowVarBound ν z (dN N)) atTop (𝓝 0) := by
    have h := (hcN.mul ht).mul hrvD
    rw [mul_zero, zero_mul] at h
    refine Tendsto.congr' ?_ h
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    field_simp
  -- the mean term
  have hmean : Tendsto (fun N => meanErr ν z (pN N) (dN N)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => (‖z‖ / z.im) ^ 4 / 2
        * ((Real.sqrt (dN N))⁻¹ * (((dN N : ℝ)) * alphaSqBound ν z (pN N) (dN N))))
        atTop (𝓝 0) := by
      have h := (hs.mul halphaD).const_mul ((‖z‖ / z.im) ^ 4 / 2)
      rw [zero_mul, mul_zero] at h
      exact h
    have h2 : Tendsto (fun N => 1 / 2
        * ((Real.sqrt (dN N))⁻¹ * (((dN N : ℝ)) ^ 2 * rowTermSqBound ν z (dN N))))
        atTop (𝓝 0) := by
      have h := (hs.mul hrtD).const_mul (1 / 2 : ℝ)
      rw [zero_mul, mul_zero] at h
      exact h
    have h3 : Tendsto (fun N => (‖z‖ / z.im) ^ 2 / z.im ^ 2 * ((dN N : ℝ))⁻¹)
        atTop (𝓝 0) := by
      have h := ht.const_mul ((‖z‖ / z.im) ^ 2 / z.im ^ 2)
      rw [mul_zero] at h
      exact h
    have hsum := hcN.mul ((h1.add h2).add h3)
    rw [add_zero, add_zero, mul_zero] at hsum
    refine Tendsto.congr' ?_ hsum
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [meanErr]
    set S : ℝ := Real.sqrt ((dN N : ℝ)) with hSdef
    have hSpos : 0 < S := by rw [hSdef]; exact Real.sqrt_pos.mpr hDR
    have hS2 : S * S = ((dN N : ℝ)) := by rw [hSdef]; exact Real.mul_self_sqrt hDR.le
    set A : ℝ := alphaSqBound ν z (pN N) (dN N) with hA
    set Rt : ℝ := rowTermSqBound ν z (dN N) with hRt
    rw [← hS2]
    field_simp
  -- assemble
  have hbase : Tendsto (fun N => isoBase ν z (pN N) (dN N)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => 2 * (pN N : ℝ) * rowVarBound ν z (dN N)) atTop (𝓝 0) := by
      have h := hPrv.const_mul (2 : ℝ)
      rw [mul_zero] at h
      exact h.congr fun N => by ring
    have h2 : Tendsto (fun N => 2 * (1 / z.im * meanErr ν z (pN N) (dN N)) ^ 2)
        atTop (𝓝 0) := by
      have h := ((hmean.const_mul (1 / z.im)).pow 2).const_mul (2 : ℝ)
      rw [mul_zero, zero_pow (two_ne_zero), mul_zero] at h
      exact h
    have h3 : Tendsto (fun N => 2 * traceSqBound ν z (pN N) (dN N)) atTop (𝓝 0) := by
      have h := htrace.const_mul (2 : ℝ)
      rw [mul_zero] at h
      exact h
    have h := (h1.add h2).add h3
    rw [add_zero, add_zero] at h
    exact h.congr fun N => by rw [isoBase]
  have h := hbase.const_mul (8 : ℝ)
  rw [mul_zero] at h
  exact h.congr fun N => by rw [isoRate]


end GenRMT

end StackedSVD
