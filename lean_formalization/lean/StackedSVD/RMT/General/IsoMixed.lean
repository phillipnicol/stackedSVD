/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Iso

/-!
# Item G4, extension: the mixed isotropic form

`notes/archive/prop_single_table_general.md` section 8 item 1, and section 2 of the brief
`brief_G4_mixed.md`. Unit G6 needs the field `vgC` of `ResolventFormsC`, which is
`vᵀ G₀ g` with `g = m.gvec`; Sherman-Morrison turns it into the **mixed** form

`F = xᵀ G Eᵀ y`,   `E = (√d)⁻¹ Y`,   `x` in the right space, `y` in the left space,

over a bounded denominator. None of the statements of `RMT/General/Iso.lean` gives this form:
`E` and `G` share the same matrix, so no direction is independent of the resolvent.

The conclusion is again a measure bound at finite `p` and `d`, with a rate free of `x`, `y`
and `ε`, because unit G7 applies it to directions that change with `N`.

## Route

Write `g k = rowVec Y k`, `A_k = gram (updateRow Y k 0)`, `G_k = resolvC A_k z`,
`s k = sRow Y z k`, `ξ k = alphaRow Y z k`.

1. `F = ∑ k, y k * (xᵀ G g k)` and `xᵀ G g k = (1 - s k) * (xᵀ G_k g k)`, so every row term is
   a bounded scalar times a linear form of that row.
2. The variance, by Efron-Stein on the rows, with the leave-one-out center
   `F_k = xᵀ G_k E_kᵀ y`. The new input is the operator bound
   `‖G_k E_kᵀ y‖² ≤ ‖z‖ (y ⬝ y) / η²`, which comes from the companion identity
   `E G Eᵀ = 1 + z Gc` and the resolvent identity `‖G w‖² = Im (wᵀ G w) / η`.
3. The mean: conditionally on the other rows, `xᵀ G_k g k` is a linear form with mean zero, so
   the deterministic `secW` may be subtracted from `1 - s k` for free, and what is left is
   `O(1/d)` per row by the `L²` bound on `ξ k - m`.
4. Chebyshev, as statements 1 and 2 do.

F36/F37 (2026-09-09): the private helpers of `RMT/General/Iso.lean` used to be invisible here,
so the ones this file needed were repeated with the suffix `Mx`, each a verbatim copy. Those
copies are now dropped: `Iso.lean` made its own twins public, and `Companion.lean` carries the
public canonical copy of the few that were shared with `FormsBridge.lean` too. This file uses
the public copies directly.
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace GenRMT

variable {p d : ℕ} {z : ℂ}

/-! ### F36/F37 (2026-09-09)

This section used to hold 23 private helpers, each named with the suffix `Mx` (for example
`integrable_of_lintegral_ofReal_ne_topMx`), verbatim copies of the private helpers of
`RMT/General/Iso.lean`. `Iso.lean` made each twin public (dropping `private`); this file uses
those public copies directly (unqualified, since both files sit in `namespace GenRMT`), and
the 23 local copies are dropped.
-/
/-! ### Milestone 1: the identity and the operator bound -/

/-- The mixed direction is the `y`-combination of the scaled rows. -/
private theorem cvec_mix (Y : Matrix (Fin p) (Fin d) ℝ) (y : Fin p → ℝ) :
    R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y)) = ∑ k, (y k : ℂ) • R4C.cvec (rowVec Y k) := by
  funext a
  simp only [R4C.cvec, Pi.smul_apply, smul_eq_mul, Matrix.mulVec, Matrix.transpose_apply,
    dotProduct, Finset.sum_apply, rowVec_apply, Complex.ofReal_sum, Complex.ofReal_mul,
    Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by push_cast; ring

/-- **The mixed form is the `y`-combination of the row forms.** -/
theorem cformC_mixed_eq_sum (Y : Matrix (Fin p) (Fin d) ℝ) (x : Fin d → ℝ) (y : Fin p → ℝ) :
    R4C.cformC (gram Y) z x ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))
      = ∑ k, (y k : ℂ) * R4C.cformC (gram Y) z x (rowVec Y k) := by
  change R4C.cvec x ⬝ᵥ (R4C.resolvC (gram Y) z *ᵥ R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))) = _
  rw [cvec_mix, Matrix.mulVec_sum, dotProduct_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul]
  rfl

/-- Cauchy-Schwarz for a complex dot product, in the squared form. -/
private theorem normSq_dotProduct_le {D : ℕ} (a b : Fin D → ℂ) :
    ‖a ⬝ᵥ b‖ ^ 2 ≤ (∑ i, ‖a i‖ ^ 2) * ∑ i, ‖b i‖ ^ 2 := by
  have h1 : ‖a ⬝ᵥ b‖ ≤ ∑ i, ‖a i‖ * ‖b i‖ := by
    refine (norm_sum_le _ _).trans (le_of_eq ?_)
    exact Finset.sum_congr rfl fun i _ => norm_mul _ _
  have h2 : (∑ i, ‖a i‖ * ‖b i‖) ^ 2 ≤ (∑ i, ‖a i‖ ^ 2) * ∑ i, ‖b i‖ ^ 2 :=
    Finset.sum_mul_sq_le_sq_mul_sq _ _ _
  have h3 : (0 : ℝ) ≤ ∑ i, ‖a i‖ * ‖b i‖ :=
    Finset.sum_nonneg fun i _ => mul_nonneg (norm_nonneg _) (norm_nonneg _)
  exact (pow_le_pow_left₀ (norm_nonneg _) h1 2).trans h2

/-- **The resolvent identity in the squared form.** `∑ a, ‖(G u) a‖² = Im (uᵀ G u) / η`. -/
theorem sum_normSq_resolvC_mulVec_eq {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (u : Fin D → ℝ) :
    ∑ a, ‖(R4C.resolvC W z *ᵥ R4C.cvec u) a‖ ^ 2 = (R4C.qformC W z u).im / z.im := by
  have hz' : z.im ≠ 0 := hz.ne'
  set V : Matrix (Fin D) (Fin D) ℂ := R4C.cmat (R4.eigU hW) with hV
  have hVr : ∀ a b, (starRingEnd ℂ) (V a b) = V a b := by
    intro a b; simp [hV, R4C.cmat]
  have hVtr : ∀ a b, (starRingEnd ℂ) (Vᵀ a b) = Vᵀ a b := by
    intro a b; simp [hV, R4C.cmat]
  set f : Fin D → ℂ := fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹ with hf
  set q : Fin D → ℂ := Vᵀ *ᵥ R4C.cvec u with hq
  have hqval : ∀ a, q a = ((((R4.eigU hW)ᵀ *ᵥ u) a : ℝ) : ℂ) := by
    intro a
    simp only [hq, hV, Matrix.mulVec, dotProduct, Matrix.transpose_apply, R4C.cmat,
      Matrix.map_apply, R4C.cvec, Complex.ofReal_sum, Complex.ofReal_mul]
  have hsplit : R4C.resolvC W z *ᵥ R4C.cvec u = V *ᵥ (Matrix.diagonal f *ᵥ q) := by
    rw [R4C.resolvC_eq_conj hW hz', ← hV, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]
  have hne : ∀ a, ((hW.eigenvalues a : ℂ) - z) ≠ 0 := fun a =>
    R4C.eigenvalue_sub_ne_zero hW hz' a
  rw [hsplit, sum_normSq_mulVec_real_orth V (R4C.transpose_ceigU_mul hW) hVr]
  have hterm : ∀ a, ‖(Matrix.diagonal f *ᵥ q) a‖ ^ 2
      = ((R4.eigU hW)ᵀ *ᵥ u) a ^ 2 / Complex.normSq ((hW.eigenvalues a : ℂ) - z) := by
    intro a
    rw [Matrix.mulVec_diagonal, norm_mul, mul_pow, hqval a, Complex.norm_real,
      Real.norm_eq_abs, sq_abs, hf]
    rw [norm_inv, ← Complex.sq_norm, inv_pow]
    ring
  rw [Finset.sum_congr rfl fun a _ => hterm a]
  have him : (R4C.qformC W z u).im
      = ∑ a, z.im * (((R4.eigU hW)ᵀ *ᵥ u) a ^ 2
          / Complex.normSq ((hW.eigenvalues a : ℂ) - z)) := by
    change (R4C.cformC W z u u).im = _
    rw [R4C.cformC_eq_sum hW hz' u u, Complex.im_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, mul_zero, zero_add,
      Complex.inv_im]
    have hIm : ((hW.eigenvalues a : ℂ) - z).im = -z.im := by simp
    rw [hIm]
    field_simp
  rw [him, ← Finset.mul_sum]
  field_simp

-- F37 (2026-09-09): `cvec_smulMx`, `cmat'_transpose_mulVec_cvecMx`, `cvec_dotProduct_cvecMx`
-- and `dotProduct_mulVec_eq_rectMx` used to be repeated here (`private`). `Companion.lean:404`
-- to `:425` now carries the public canonical copy of each, in `namespace SpikedModel`; this
-- file uses those directly (qualified, as `SpikedModel.cvec_smul` and so on).


/-- **The companion identity at the mixed direction.** `(Eᵀ y)ᵀ G (Eᵀ y) = (y ⬝ y) + z yᵀ Gc y`,
the twin of `SpikedModel.qformC_gvec_eq` (`RMT/General/Companion.lean:429`) at a general `y`. -/
theorem qformC_gram_mix_eq (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (hd : 0 < d)
    (y : Fin p → ℝ) :
    R4C.qformC (gram Y) z ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))
      = ((y ⬝ᵥ y : ℝ) : ℂ) + z * R4C.qformC (gramC Y) z y := by
  have hscal : (((Real.sqrt d)⁻¹ : ℝ) : ℂ) * (((Real.sqrt d)⁻¹ : ℝ) : ℂ) = ((d : ℕ) : ℂ)⁻¹ := by
    rw [← Complex.ofReal_mul]
    have step : (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = ((d : ℝ))⁻¹ := by
      rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d)]
    rw [step, Complex.ofReal_inv]
    congr 1
  have hcvw : R4C.cvec (Yᵀ *ᵥ y) = (R4C.cmat' Y)ᵀ *ᵥ R4C.cvec y :=
    (SpikedModel.cmat'_transpose_mulVec_cvec Y y).symm
  have hbridge : R4C.cvec y
      ⬝ᵥ ((R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ) *ᵥ R4C.cvec y)
      = ((R4C.cmat' Y)ᵀ *ᵥ R4C.cvec y)
        ⬝ᵥ (R4C.resolvC (gram Y) z *ᵥ ((R4C.cmat' Y)ᵀ *ᵥ R4C.cvec y)) := by
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
      SpikedModel.dotProduct_mulVec_eq_rect (R4C.cmat' Y) (R4C.cvec y)
        (R4C.resolvC (gram Y) z *ᵥ ((R4C.cmat' Y)ᵀ *ᵥ R4C.cvec y))]
  have hq1 : R4C.qformC (gram Y) z ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))
      = ((d : ℕ) : ℂ)⁻¹ * (R4C.cvec y
          ⬝ᵥ ((R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ) *ᵥ R4C.cvec y)) := by
    change R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))
      ⬝ᵥ (R4C.resolvC (gram Y) z *ᵥ R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))) = _
    rw [SpikedModel.cvec_smul, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_smul,
      smul_eq_mul, hscal, hcvw, hbridge]
  rw [hq1]
  have h := smul_cmat'_mul_resolvC_mul_transpose Y hz hd
  have hswap : ((d : ℕ) : ℂ)⁻¹ * (R4C.cvec y
        ⬝ᵥ ((R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ) *ᵥ R4C.cvec y))
      = R4C.cvec y ⬝ᵥ ((1 + z • R4C.resolvC (gramC Y) z) *ᵥ R4C.cvec y) := by
    have h2 : R4C.cvec y
        ⬝ᵥ ((((d : ℕ) : ℂ)⁻¹
          • (R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ)) *ᵥ R4C.cvec y)
        = R4C.cvec y ⬝ᵥ ((1 + z • R4C.resolvC (gramC Y) z) *ᵥ R4C.cvec y) := by rw [h]
    rwa [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h2
  rw [hswap, Matrix.add_mulVec, dotProduct_add, Matrix.one_mulVec, Matrix.smul_mulVec,
    dotProduct_smul, smul_eq_mul, SpikedModel.cvec_dotProduct_cvec]
  rfl

/-- **The operator bound.** `‖G Eᵀ y‖² ≤ ‖z‖ (y ⬝ y)/η²`, uniform in `Y`: the naive bound
`‖Eᵀ y‖²/η²` is useless, because `‖Eᵀ y‖` is random and unbounded. -/
theorem sum_normSq_resolvC_mulVec_mix_le (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im)
    (hd : 0 < d) (y : Fin p → ℝ) :
    ∑ a, ‖(R4C.resolvC (gram Y) z *ᵥ R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))) a‖ ^ 2
      ≤ ‖z‖ * (y ⬝ᵥ y) / z.im ^ 2 := by
  have hyy : (0 : ℝ) ≤ y ⬝ᵥ y := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rw [sum_normSq_resolvC_mulVec_eq (gram_isHermitian Y) hz, qformC_gram_mix_eq Y hz hd y]
  have h1 : (((y ⬝ᵥ y : ℝ) : ℂ) + z * R4C.qformC (gramC Y) z y).im
      = (z * R4C.qformC (gramC Y) z y).im := by simp
  rw [h1]
  have h2 : (z * R4C.qformC (gramC Y) z y).im ≤ ‖z * R4C.qformC (gramC Y) z y‖ :=
    le_trans (le_abs_self _) (Complex.abs_im_le_norm _)
  have h3 : ‖z * R4C.qformC (gramC Y) z y‖ ≤ ‖z‖ * ((y ⬝ᵥ y) / z.im) := by
    rw [norm_mul]
    exact mul_le_mul_of_nonneg_left (R4C.norm_qformC_le (gramC_isHermitian Y) hz y)
      (norm_nonneg z)
  have h4 : (z * R4C.qformC (gramC Y) z y).im ≤ ‖z‖ * ((y ⬝ᵥ y) / z.im) := le_trans h2 h3
  rw [div_le_div_iff₀ hz (by positivity : (0:ℝ) < z.im ^ 2)]
  have h5 : ‖z‖ * ((y ⬝ᵥ y) / z.im) * z.im ^ 2 = ‖z‖ * (y ⬝ᵥ y) * z.im := by
    field_simp
  nlinarith [h4, hz, sq_nonneg z.im]

/-- **The mixed form** `xᵀ G Eᵀ y`, with `Eᵀ y = (√d)⁻¹ • (Yᵀ *ᵥ y)`, the shape of `gvec`. -/
noncomputable def mixF (z : ℂ) (x : Fin d → ℝ) (y : Fin p → ℝ)
    (Y : Matrix (Fin p) (Fin d) ℝ) : ℂ :=
  R4C.cformC (gram Y) z x ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))

/-- The mixed form is bounded, uniformly in the realization. -/
theorem normSq_mixF_le (hz : 0 < z.im) (hd : 0 < d) {x : Fin d → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin p → ℝ} (hy : y ⬝ᵥ y ≤ 1) (Y : Matrix (Fin p) (Fin d) ℝ) :
    ‖mixF z x y Y‖ ^ 2 ≤ ‖z‖ / z.im ^ 2 := by
  have hyy : (0 : ℝ) ≤ y ⬝ᵥ y := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hxx : (0 : ℝ) ≤ x ⬝ᵥ x := by
    rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  have hcs := normSq_dotProduct_le (R4C.cvec x)
    (R4C.resolvC (gram Y) z *ᵥ R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y)))
  have hxsum : ∑ a, ‖R4C.cvec x a‖ ^ 2 = x ⬝ᵥ x := by
    simp only [R4C.cvec, Complex.norm_real, Real.norm_eq_abs, sq_abs, dotProduct]
    exact Finset.sum_congr rfl fun a _ => by ring
  rw [hxsum] at hcs
  have hbnd := sum_normSq_resolvC_mulVec_mix_le Y hz hd y
  have hsnn : (0 : ℝ) ≤ ∑ a, ‖(R4C.resolvC (gram Y) z
      *ᵥ R4C.cvec ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))) a‖ ^ 2 :=
    Finset.sum_nonneg fun a _ => sq_nonneg _
  have hzz : (0 : ℝ) < z.im ^ 2 := by positivity
  have hznn : (0 : ℝ) ≤ ‖z‖ := norm_nonneg z
  have hstep : ‖mixF z x y Y‖ ^ 2 ≤ (x ⬝ᵥ x) * (‖z‖ * (y ⬝ᵥ y) / z.im ^ 2) := by
    refine le_trans hcs ?_
    exact mul_le_mul_of_nonneg_left hbnd hxx
  refine hstep.trans ?_
  have hfrac : ‖z‖ * (y ⬝ᵥ y) / z.im ^ 2 ≤ ‖z‖ / z.im ^ 2 := by
    apply div_le_div_of_nonneg_right _ (le_of_lt hzz)
    nlinarith
  nlinarith [hfrac, div_nonneg (mul_nonneg hznn hyy) hzz.le]


/-! ### Milestone 2: the leave-one-out decomposition -/

/-- The mixed direction splits off the row `k` part. -/
theorem mix_split (Y : Matrix (Fin p) (Fin d) ℝ) (y : Fin p → ℝ) (k : Fin p) :
    (Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y)
      = (Real.sqrt d)⁻¹ • ((Matrix.updateRow Y k 0)ᵀ *ᵥ y) + (y k) • rowVec Y k := by
  funext a
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Matrix.mulVec, Matrix.transpose_apply,
    dotProduct, rowVec_apply]
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_erase_add _ _ (Finset.mem_univ k)]
  have hz0 : ∀ j ∈ Finset.univ.erase k,
      (Real.sqrt d)⁻¹ * ((Matrix.updateRow Y k 0) j a * y j)
        = (Real.sqrt d)⁻¹ * (Y j a * y j) := by
    intro j hj
    rw [Matrix.updateRow_ne (Finset.ne_of_mem_erase hj)]
  rw [← Finset.sum_erase_add _ (fun j => (Real.sqrt d)⁻¹ * ((Matrix.updateRow Y k 0) j a * y j))
    (Finset.mem_univ k), Finset.sum_congr rfl hz0, Matrix.updateRow_self]
  simp
  ring

/-- **The leave-one-out decomposition of the mixed form.** `F - F_k` is a bounded scalar times a
product of two linear forms of row `k`. -/
theorem mixF_sub_loo (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (x : Fin d → ℝ)
    (y : Fin p → ℝ) (k : Fin p) :
    mixF z x y Y - mixF z x y (Matrix.updateRow Y k 0)
      = (1 - sRow Y z k) * (R4C.cformC (gram (Matrix.updateRow Y k 0)) z x (rowVec Y k)
          * (((y k : ℝ) : ℂ)
            - R4C.cformC (gram (Matrix.updateRow Y k 0)) z (rowVec Y k)
                ((Real.sqrt d)⁻¹ • ((Matrix.updateRow Y k 0)ᵀ *ᵥ y)))) := by
  have hne : 1 - sRow Y z k ≠ 0 := one_sub_sRow_ne_zero Y hz k
  have hcomm : ∀ a b : Fin d → ℝ, R4C.cformC (gram Y) z a b = R4C.cformC (gram Y) z b a :=
    fun a b => cformC_comm (gram_isHermitian Y) hz.ne' a b
  set Xk : Matrix (Fin p) (Fin d) ℝ := Matrix.updateRow Y k 0 with hXk
  set wk : Fin d → ℝ := (Real.sqrt d)⁻¹ • (Xkᵀ *ᵥ y) with hwk
  set gk : Fin d → ℝ := rowVec Y k with hgk
  -- `mixF Y` splits
  have hsplitF : mixF z x y Y
      = R4C.cformC (gram Y) z x wk + ((y k : ℝ) : ℂ) * R4C.cformC (gram Y) z x gk := by
    rw [mixF, mix_split Y y k, ← hXk, ← hwk, ← hgk]
    rw [hcomm x (wk + (y k) • gk), cformC_add_left,
      cformC_smulR_left (gram Y) z (y k) gk x, hcomm wk x, hcomm gk x]
  -- the downdate at `(x, wk)`
  have hdown : R4C.cformC (gram Xk) z x wk
      = R4C.cformC (gram Y) z x wk
        + R4C.cformC (gram Y) z x gk * R4C.cformC (gram Y) z gk wk / (1 - sRow Y z k) := by
    rw [hXk, gram_updateRow_zero]
    exact R4C.cformC_sub_vecMulVec (gram_isHermitian Y) hz (rowVec Y k) x wk
  have hxg : R4C.cformC (gram Y) z x gk
      = (1 - sRow Y z k) * R4C.cformC (gram Xk) z x gk := cformC_rowVec_right_eq Y hz x k
  have hgw : R4C.cformC (gram Y) z gk wk
      = (1 - sRow Y z k) * R4C.cformC (gram Xk) z gk wk := cformC_rowVec_eq Y hz wk k
  have hmixk : mixF z x y Xk = R4C.cformC (gram Xk) z x wk := rfl
  rw [hsplitF, hmixk]
  rw [hxg] at hdown ⊢
  rw [hgw] at hdown
  field_simp at hdown ⊢
  linear_combination -hdown


/-- The per-row second moment of the leave-one-out deviation. -/
noncomputable def mixRowBound (ν : Measure ℝ) (z : ℂ) (D : ℕ) (t : ℝ) : ℝ :=
  2 * (‖z‖ / z.im) ^ 2 * (t ^ 2 / ((D : ℝ) * z.im ^ 2))
    + 2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2)
      / ((D : ℝ) ^ 2 * z.im ^ 4)

/-- **The second moment of the leave-one-out deviation in one row.** -/
private theorem integral_normSq_mixF_sub_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin P → ℝ} (hy : y ⬝ᵥ y ≤ 1) (k : Fin P) (x0 : Matrix (Fin P) (Fin D) ℝ) :
    ∫ t, ‖mixF z x y (Matrix.updateRow x0 k t) - mixF z x y (Matrix.updateRow x0 k 0)‖ ^ 2
        ∂(Measure.pi fun _ : Fin D => ν)
      ≤ mixRowBound ν z D (y k) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hnu : (0 : ℝ) ≤ (∫ w, w ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  have hsqi : ((Real.sqrt D)⁻¹ : ℝ) ^ 2 = ((D : ℝ))⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hDR.le]
  set Xk : Matrix (Fin P) (Fin D) ℝ := Matrix.updateRow x0 k 0 with hXk
  set wk : Fin D → ℝ := (Real.sqrt D)⁻¹ • (Xkᵀ *ᵥ y) with hwk
  set aa : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram Xk) z *ᵥ R4C.cvec x) l with haa
  set bb : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram Xk) z *ᵥ R4C.cvec wk) l with hbb
  have hcoef : ∀ (v : Fin D → ℝ) (S : ℝ),
      (∑ l, ‖(R4C.resolvC (gram Xk) z *ᵥ R4C.cvec v) l‖ ^ 2) ≤ S →
      ∑ l, ‖(((Real.sqrt D)⁻¹ : ℝ) : ℂ)
          * (R4C.resolvC (gram Xk) z *ᵥ R4C.cvec v) l‖ ^ 2 ≤ ((D : ℝ))⁻¹ * S := by
    intro v S hS
    have hstep : ∀ l, ‖(((Real.sqrt D)⁻¹ : ℝ) : ℂ)
        * (R4C.resolvC (gram Xk) z *ᵥ R4C.cvec v) l‖ ^ 2
        = ((D : ℝ))⁻¹ * ‖(R4C.resolvC (gram Xk) z *ᵥ R4C.cvec v) l‖ ^ 2 := by
      intro l
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs, hsqi]
    rw [Finset.sum_congr rfl fun l _ => hstep l, ← Finset.mul_sum]
    exact mul_le_mul_of_nonneg_left hS (by positivity)
  have hA : ∑ l, ‖aa l‖ ^ 2 ≤ 1 / ((D : ℝ) * z.im ^ 2) := by
    have h : ∑ l, ‖aa l‖ ^ 2 ≤ ((D : ℝ))⁻¹ * ((x ⬝ᵥ x) / z.im ^ 2) :=
      hcoef x ((x ⬝ᵥ x) / z.im ^ 2)
        (sum_normSq_resolvC_mulVec_le (gram_isHermitian Xk) hz x)
    refine h.trans ?_
    have h2 : (x ⬝ᵥ x) / z.im ^ 2 ≤ 1 / z.im ^ 2 :=
      div_le_div_of_nonneg_right hx (by positivity)
    calc ((D : ℝ))⁻¹ * ((x ⬝ᵥ x) / z.im ^ 2) ≤ ((D : ℝ))⁻¹ * (1 / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h2 (by positivity)
      _ = 1 / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hB : ∑ l, ‖bb l‖ ^ 2 ≤ ‖z‖ / ((D : ℝ) * z.im ^ 2) := by
    have h : ∑ l, ‖bb l‖ ^ 2 ≤ ((D : ℝ))⁻¹ * (‖z‖ * (y ⬝ᵥ y) / z.im ^ 2) :=
      hcoef wk (‖z‖ * (y ⬝ᵥ y) / z.im ^ 2)
        (sum_normSq_resolvC_mulVec_mix_le Xk hz hD y)
    refine h.trans ?_
    have hyy : (0 : ℝ) ≤ y ⬝ᵥ y := by
      rw [dotProduct]; exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
    have h2 : ‖z‖ * (y ⬝ᵥ y) / z.im ^ 2 ≤ ‖z‖ / z.im ^ 2 := by
      refine div_le_div_of_nonneg_right ?_ (by positivity)
      nlinarith [norm_nonneg z]
    calc ((D : ℝ))⁻¹ * (‖z‖ * (y ⬝ᵥ y) / z.im ^ 2) ≤ ((D : ℝ))⁻¹ * (‖z‖ / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h2 (by positivity)
      _ = ‖z‖ / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hAnn : (0 : ℝ) ≤ ∑ l, ‖aa l‖ ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  have hBnn : (0 : ℝ) ≤ ∑ l, ‖bb l‖ ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  set K : ℝ := ‖z‖ / z.im with hK
  set g : (Fin D → ℝ) → ℝ := fun t =>
    2 * K ^ 2 * ((y k) ^ 2 * ‖∑ l, aa l * (t l : ℂ)‖ ^ 2)
      + K ^ 2 * (‖∑ l, aa l * (t l : ℂ)‖ ^ 4 + ‖∑ l, bb l * (t l : ℂ)‖ ^ 4) with hg
  have hpt : ∀ t : Fin D → ℝ,
      ‖mixF z x y (Matrix.updateRow x0 k t) - mixF z x y (Matrix.updateRow x0 k 0)‖ ^ 2
        ≤ g t := by
    intro t
    have hupd : Matrix.updateRow (Matrix.updateRow x0 k t) k 0 = Xk := by
      rw [hXk]; exact updateRow_idem x0 k t 0
    have hrow : rowVec (Matrix.updateRow x0 k t) k = (Real.sqrt D)⁻¹ • t := by
      have hYk : (Matrix.updateRow x0 k t) k = t := by rw [Matrix.updateRow_self]
      rw [rowVec, hYk]
    have hdec := mixF_sub_loo (Matrix.updateRow x0 k t) hz x y k
    rw [hupd, hrow] at hdec
    have hPeq : R4C.cformC (gram Xk) z x ((Real.sqrt D)⁻¹ • t) = ∑ l, aa l * (t l : ℂ) := by
      rw [cformC_comm (gram_isHermitian Xk) hz.ne', cformC_smul_left_eq]
    have hQeq : R4C.cformC (gram Xk) z ((Real.sqrt D)⁻¹ • t) wk = ∑ l, bb l * (t l : ℂ) := by
      rw [cformC_smul_left_eq]
    rw [hdec, hPeq, ← hwk, hQeq]
    set Pv : ℂ := ∑ l, aa l * (t l : ℂ)
    set Qv : ℂ := ∑ l, bb l * (t l : ℂ)
    have hs : ‖1 - sRow (Matrix.updateRow x0 k t) z k‖ ≤ K :=
      norm_one_sub_sRow_le _ hz hD k
    have hknn : (0 : ℝ) ≤ K := by rw [hK]; positivity
    have hb1 : ‖(1 - sRow (Matrix.updateRow x0 k t) z k) * (Pv * (((y k : ℝ) : ℂ) - Qv))‖
        ≤ K * (‖Pv‖ * (|y k| + ‖Qv‖)) := by
      rw [norm_mul, norm_mul]
      refine mul_le_mul hs ?_ (by positivity) hknn
      refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
      refine (norm_sub_le _ _).trans (le_of_eq ?_)
      rw [Complex.norm_real, Real.norm_eq_abs]
    have hsq := pow_le_pow_left₀ (norm_nonneg _) hb1 2
    refine hsq.trans ?_
    rw [hg]
    have hk2 : (0 : ℝ) ≤ K ^ 2 := sq_nonneg K
    have e1 : (|y k| + ‖Qv‖) ^ 2 ≤ 2 * ((y k) ^ 2 + ‖Qv‖ ^ 2) := by
      nlinarith [sq_nonneg (|y k| - ‖Qv‖), sq_abs (y k)]
    have e2 : 2 * (‖Pv‖ ^ 2 * ‖Qv‖ ^ 2) ≤ ‖Pv‖ ^ 4 + ‖Qv‖ ^ 4 := by
      nlinarith [sq_nonneg (‖Pv‖ ^ 2 - ‖Qv‖ ^ 2)]
    calc (K * (‖Pv‖ * (|y k| + ‖Qv‖))) ^ 2
        = K ^ 2 * ‖Pv‖ ^ 2 * (|y k| + ‖Qv‖) ^ 2 := by ring
      _ ≤ K ^ 2 * ‖Pv‖ ^ 2 * (2 * ((y k) ^ 2 + ‖Qv‖ ^ 2)) :=
          mul_le_mul_of_nonneg_left e1 (by positivity)
      _ = 2 * K ^ 2 * ((y k) ^ 2 * ‖Pv‖ ^ 2) + K ^ 2 * (2 * (‖Pv‖ ^ 2 * ‖Qv‖ ^ 2)) := by ring
      _ ≤ 2 * K ^ 2 * ((y k) ^ 2 * ‖Pv‖ ^ 2) + K ^ 2 * (‖Pv‖ ^ 4 + ‖Qv‖ ^ 4) := by
          have h := mul_le_mul_of_nonneg_left e2 hk2
          linarith
  obtain ⟨h4a, hb4a⟩ := integral_clinForm_pow_four_le hν aa
  obtain ⟨h4b, hb4b⟩ := integral_clinForm_pow_four_le hν bb
  have h2a := integrable_normSq_clinForm hν aa
  have hI1 : Integrable (fun t : Fin D → ℝ =>
      2 * K ^ 2 * ((y k) ^ 2 * ‖∑ l, aa l * (t l : ℂ)‖ ^ 2))
      (Measure.pi fun _ : Fin D => ν) := (h2a.const_mul ((y k) ^ 2)).const_mul (2 * K ^ 2)
  have hI2 : Integrable (fun t : Fin D → ℝ =>
      K ^ 2 * (‖∑ l, aa l * (t l : ℂ)‖ ^ 4 + ‖∑ l, bb l * (t l : ℂ)‖ ^ 4))
      (Measure.pi fun _ : Fin D => ν) := (h4a.add h4b).const_mul (K ^ 2)
  have hgint : Integrable g (Measure.pi fun _ : Fin D => ν) := by
    simp only [hg]
    exact hI1.add hI2
  refine (integral_mono_of_nonneg (Filter.Eventually.of_forall fun t => sq_nonneg _) hgint
    (Filter.Eventually.of_forall hpt)).trans ?_
  simp only [hg]
  rw [integral_add hI1 hI2, integral_const_mul, integral_const_mul, integral_const_mul,
    integral_add h4a h4b, integral_normSq_clinForm hν aa]
  have hka : (0 : ℝ) ≤ 2 * K ^ 2 * (y k) ^ 2 := by positivity
  have hstepA : 2 * K ^ 2 * ((y k) ^ 2 * ∑ l, ‖aa l‖ ^ 2)
      ≤ 2 * (‖z‖ / z.im) ^ 2 * ((y k) ^ 2 / ((D : ℝ) * z.im ^ 2)) := by
    calc 2 * K ^ 2 * ((y k) ^ 2 * ∑ l, ‖aa l‖ ^ 2)
        = 2 * K ^ 2 * (y k) ^ 2 * ∑ l, ‖aa l‖ ^ 2 := by ring
      _ ≤ 2 * K ^ 2 * (y k) ^ 2 * (1 / ((D : ℝ) * z.im ^ 2)) :=
          mul_le_mul_of_nonneg_left hA hka
      _ = 2 * (‖z‖ / z.im) ^ 2 * ((y k) ^ 2 / ((D : ℝ) * z.im ^ 2)) := by rw [hK]; ring
  have hstepB : K ^ 2 * ((∫ t, ‖∑ l, aa l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
        + ∫ t, ‖∑ l, bb l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
      ≤ 2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2)
        / ((D : ℝ) ^ 2 * z.im ^ 4) := by
    have hA2 : (∑ l, ‖aa l‖ ^ 2) ^ 2 ≤ (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 :=
      pow_le_pow_left₀ hAnn hA 2
    have hB2 : (∑ l, ‖bb l‖ ^ 2) ^ 2 ≤ (‖z‖ / ((D : ℝ) * z.im ^ 2)) ^ 2 :=
      pow_le_pow_left₀ hBnn hB 2
    have hfa : ∫ t, ‖∑ l, aa l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
        ≤ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2 :=
      hb4a.trans (mul_le_mul_of_nonneg_left hA2 (by linarith))
    have hfb : ∫ t, ‖∑ l, bb l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
        ≤ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (‖z‖ / ((D : ℝ) * z.im ^ 2)) ^ 2 :=
      hb4b.trans (mul_le_mul_of_nonneg_left hB2 (by linarith))
    have hknn2 : (0 : ℝ) ≤ K ^ 2 := sq_nonneg K
    have hcomb : K ^ 2 * ((∫ t, ‖∑ l, aa l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
          + ∫ t, ‖∑ l, bb l * (t l : ℂ)‖ ^ 4 ∂(Measure.pi fun _ : Fin D => ν))
        ≤ K ^ 2 * (2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 / ((D : ℝ) * z.im ^ 2)) ^ 2
          + 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (‖z‖ / ((D : ℝ) * z.im ^ 2)) ^ 2) :=
      mul_le_mul_of_nonneg_left (by linarith) hknn2
    refine hcomb.trans (le_of_eq ?_)
    rw [hK]
    field_simp
  rw [mixRowBound]
  linarith


/-! ### Milestone 3: the variance -/

/-- One coordinate of the mixed direction `(√d)⁻¹ (Yᵀ y)` is measurable in `Y`. Confirmed
canonical for `RMT/General/` by F37 (2026-09-09): the private twins `measurable_mixDir`
(`FormsGeneral.lean`) and `measurable_mixDir'` (`DelocUniform.lean`) are dropped in its
favor. -/
theorem measurable_mix_dir {P D : ℕ} (y : Fin P → ℝ) (i : Fin D) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ =>
      ((Real.sqrt D)⁻¹ • (Yᵀ *ᵥ y)) i := by
  have h : ∀ Y : Matrix (Fin P) (Fin D) ℝ, ((Real.sqrt D)⁻¹ • (Yᵀ *ᵥ y)) i
      = (Real.sqrt D)⁻¹ * ∑ k, Y k i * y k := by
    intro Y
    simp only [Pi.smul_apply, smul_eq_mul, Matrix.mulVec, Matrix.transpose_apply, dotProduct]
  simp only [h]
  exact measurable_const.mul
    (Finset.measurable_sum _ fun k _ => (measurable_matrix_entry k i).mul measurable_const)

theorem measurable_mixF {P D : ℕ} (z : ℂ) (x : Fin D → ℝ) (y : Fin P → ℝ) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => mixF z x y Y := by
  simp only [mixF]
  exact measurable_cformC measurable_gram_self z (fun i => measurable_const)
    (fun i => measurable_mix_dir y i)

private theorem measurable_mixF_row {P D : ℕ} (z : ℂ) (x : Fin D → ℝ) (y : Fin P → ℝ)
    (k : Fin P) (x0 : Matrix (Fin P) (Fin D) ℝ) :
    Measurable fun t : Fin D → ℝ => mixF z x y (Matrix.updateRow x0 k t) := by
  have hWm : ∀ i j, Measurable fun t : Fin D → ℝ => gram (Matrix.updateRow x0 k t) i j :=
    fun i j => measurable_gram_entry (fun l i' => measurable_updateRow_entry' x0 k l i') i j
  have hvm : ∀ i, Measurable fun t : Fin D → ℝ =>
      ((Real.sqrt D)⁻¹ • ((Matrix.updateRow x0 k t)ᵀ *ᵥ y)) i := by
    intro i
    have h : ∀ t : Fin D → ℝ, ((Real.sqrt D)⁻¹ • ((Matrix.updateRow x0 k t)ᵀ *ᵥ y)) i
        = (Real.sqrt D)⁻¹ * ∑ j, (Matrix.updateRow x0 k t) j i * y j := by
      intro t
      simp only [Pi.smul_apply, smul_eq_mul, Matrix.mulVec, Matrix.transpose_apply, dotProduct]
    simp only [h]
    exact measurable_const.mul (Finset.measurable_sum _ fun j _ =>
      (measurable_updateRow_entry' x0 k j i).mul measurable_const)
  simp only [mixF]
  exact measurable_cformC hWm z (fun i => measurable_const) hvm

/-- The mixed form is bounded by `√(‖z‖)/η`. -/
theorem norm_mixF_le (hz : 0 < z.im) (hd : 0 < d) {x : Fin d → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin p → ℝ} (hy : y ⬝ᵥ y ≤ 1) (Y : Matrix (Fin p) (Fin d) ℝ) :
    ‖mixF z x y Y‖ ≤ Real.sqrt (‖z‖ / z.im ^ 2) := by
  have h := normSq_mixF_le hz hd hx hy Y
  have hnn : (0 : ℝ) ≤ ‖mixF z x y Y‖ := norm_nonneg _
  rw [show ‖mixF z x y Y‖ = Real.sqrt (‖mixF z x y Y‖ ^ 2) from (Real.sqrt_sq hnn).symm]
  exact Real.sqrt_le_sqrt h

/-- The row sum of the per-row bound. -/
noncomputable def mixVarBound (ν : Measure ℝ) (z : ℂ) (P D : ℕ) : ℝ :=
  2 * (‖z‖ / z.im) ^ 2 * (1 / ((D : ℝ) * z.im ^ 2))
    + (P : ℝ) * (2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2)
      / ((D : ℝ) ^ 2 * z.im ^ 4))

private theorem sum_mixRowBound_le {P D : ℕ} (ν : Measure ℝ) (z : ℂ) (hz : 0 < z.im)
    {y : Fin P → ℝ} (hy : y ⬝ᵥ y ≤ 1) :
    ∑ k, mixRowBound ν z D (y k) ≤ mixVarBound ν z P D := by
  have hnu : (0 : ℝ) ≤ (∫ w, w ^ 4 ∂ν) + 3 := by
    have := integral_pow_four_nonneg ν; linarith
  have hysum : ∑ k, (y k) ^ 2 = y ⬝ᵥ y := by
    rw [dotProduct]; exact Finset.sum_congr rfl fun k _ => by ring
  have hsplit : ∑ k, mixRowBound ν z D (y k)
      = 2 * (‖z‖ / z.im) ^ 2 * ((∑ k, (y k) ^ 2) / ((D : ℝ) * z.im ^ 2))
        + (P : ℝ) * (2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2)
          / ((D : ℝ) ^ 2 * z.im ^ 4)) := by
    simp only [mixRowBound]
    rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul]
    congr 1
    rw [Finset.sum_div, Finset.mul_sum]
  rw [hsplit, hysum, mixVarBound]
  have hcoef : (0 : ℝ) ≤ 2 * (‖z‖ / z.im) ^ 2 := by positivity
  have hstep : 2 * (‖z‖ / z.im) ^ 2 * ((y ⬝ᵥ y) / ((D : ℝ) * z.im ^ 2))
      ≤ 2 * (‖z‖ / z.im) ^ 2 * (1 / ((D : ℝ) * z.im ^ 2)) := by
    refine mul_le_mul_of_nonneg_left ?_ hcoef
    exact div_le_div_of_nonneg_right hy (by positivity)
  linarith

/-- **The variance of one real component of the mixed form**, by Efron-Stein on the rows. -/
private theorem variance_component_mix_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin P → ℝ} (hy : y ⬝ᵥ y ≤ 1)
    (φ : ℂ → ℝ) (hφm : Measurable φ) (hφs : ∀ w v : ℂ, φ w - φ v = φ (w - v))
    (hφn : ∀ w : ℂ, |φ w| ≤ ‖w‖) :
    ∫ Y, (φ (mixF z x y Y) - ∫ X, φ (mixF z x y X) ∂(noiseMatrix ν P D)) ^ 2
        ∂(noiseMatrix ν P D)
      ≤ mixVarBound ν z P D := by
  have hprob := hν.prob
  set Bd : ℝ := Real.sqrt (‖z‖ / z.im ^ 2) with hBd
  have hBdnn : (0 : ℝ) ≤ Bd := Real.sqrt_nonneg _
  have hf2 : ∀ Y : Matrix (Fin P) (Fin D) ℝ, ‖mixF z x y Y‖ ≤ Bd := fun Y =>
    norm_mixF_le hz hD hx hy Y
  set f : Matrix (Fin P) (Fin D) ℝ → ℝ := fun Y => φ (mixF z x y Y) with hf
  have hfm : Measurable f := hφm.comp (measurable_mixF z x y)
  have hbd : ∀ Y : Matrix (Fin P) (Fin D) ℝ, |f Y| ≤ Bd := fun Y => (hφn _).trans (hf2 Y)
  have hfsq : Integrable (fun Y => f Y ^ 2) (noiseMatrix ν P D) := by
    refine Integrable.mono' (integrable_const (μ := noiseMatrix ν P D) (Bd ^ 2))
      (hfm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbd Y) 2
  have hmem : MemLp f 2 (noiseMatrix ν P D) :=
    (memLp_two_iff_integrable_sq hfm.aestronglyMeasurable).mpr hfsq
  have hES : ∫ Y, (f Y - ∫ X, f X ∂(noiseMatrix ν P D)) ^ 2 ∂(noiseMatrix ν P D)
      ≤ ∑ k : Fin P, ∫ x0, (∫ t, (f (Matrix.updateRow x0 k t)
            - ∫ s, f (Matrix.updateRow x0 k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
            ∂(Measure.pi fun _ : Fin D => ν)) ∂(noiseMatrix ν P D) :=
    Tensorization.variance_pi_le_sum
      (fun _ : Fin P => Measure.pi fun _ : Fin D => ν) f hmem
  refine hES.trans ?_
  have hone : (noiseMatrix ν P D).real Set.univ = 1 := by simp
  have hinner : ∀ (k : Fin P) (x0 : Matrix (Fin P) (Fin D) ℝ),
      ∫ t, (f (Matrix.updateRow x0 k t)
          - ∫ s, f (Matrix.updateRow x0 k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
          ∂(Measure.pi fun _ : Fin D => ν)
        ≤ mixRowBound ν z D (y k) := by
    intro k x0
    have hmixm : Measurable fun t : Fin D → ℝ => mixF z x y (Matrix.updateRow x0 k t) :=
      measurable_mixF_row z x y k x0
    have hgm : Measurable fun t : Fin D → ℝ => f (Matrix.updateRow x0 k t) := hφm.comp hmixm
    have hgbd : ∀ t : Fin D → ℝ, |f (Matrix.updateRow x0 k t)| ≤ Bd := fun t =>
      (hφn _).trans (hf2 _)
    have hgint : Integrable (fun t : Fin D → ℝ => f (Matrix.updateRow x0 k t))
        (Measure.pi fun _ : Fin D => ν) :=
      Integrable.mono' (integrable_const Bd) hgm.aestronglyMeasurable
        (Filter.Eventually.of_forall fun t => by rw [Real.norm_eq_abs]; exact hgbd t)
    have hgsq : Integrable (fun t : Fin D → ℝ => f (Matrix.updateRow x0 k t) ^ 2)
        (Measure.pi fun _ : Fin D => ν) := by
      refine Integrable.mono' (integrable_const (Bd ^ 2))
        (hgm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun t => ?_)
      rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
      exact pow_le_pow_left₀ (abs_nonneg _) (hgbd t) 2
    have hdm : Measurable fun t : Fin D → ℝ =>
        ‖mixF z x y (Matrix.updateRow x0 k t) - mixF z x y (Matrix.updateRow x0 k 0)‖ ^ 2 :=
      ((hmixm.sub measurable_const).norm).pow_const 2
    have hdint : Integrable (fun t : Fin D → ℝ =>
        ‖mixF z x y (Matrix.updateRow x0 k t) - mixF z x y (Matrix.updateRow x0 k 0)‖ ^ 2)
        (Measure.pi fun _ : Fin D => ν) := by
      refine Integrable.mono' (integrable_const ((2 * Bd) ^ 2))
        hdm.aestronglyMeasurable (Filter.Eventually.of_forall fun t => ?_)
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      refine pow_le_pow_left₀ (norm_nonneg _) ?_ 2
      refine (norm_sub_le _ _).trans ?_
      have h1 := hf2 (Matrix.updateRow x0 k t)
      have h2 := hf2 (Matrix.updateRow x0 k (0 : Fin D → ℝ))
      linarith
    refine (integral_sq_sub_mean_le hgint hgsq (f (Matrix.updateRow x0 k 0))).trans ?_
    refine (integral_mono_of_nonneg (Filter.Eventually.of_forall fun t => sq_nonneg _)
      hdint (Filter.Eventually.of_forall fun t => ?_)).trans
      (integral_normSq_mixF_sub_le hν hz hD hx hy k x0)
    rw [hf, hφs]
    calc φ (mixF z x y (Matrix.updateRow x0 k t) - mixF z x y (Matrix.updateRow x0 k 0)) ^ 2
        = |φ (mixF z x y (Matrix.updateRow x0 k t)
            - mixF z x y (Matrix.updateRow x0 k 0))| ^ 2 := by rw [sq_abs]
      _ ≤ ‖mixF z x y (Matrix.updateRow x0 k t)
            - mixF z x y (Matrix.updateRow x0 k 0)‖ ^ 2 :=
          pow_le_pow_left₀ (abs_nonneg _) (hφn _) 2
  have hterm : ∀ k : Fin P,
      ∫ x0, (∫ t, (f (Matrix.updateRow x0 k t)
          - ∫ s, f (Matrix.updateRow x0 k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
          ∂(Measure.pi fun _ : Fin D => ν)) ∂(noiseMatrix ν P D)
        ≤ mixRowBound ν z D (y k) := by
    intro k
    have hnn : ∀ x0 : Matrix (Fin P) (Fin D) ℝ, (0 : ℝ)
        ≤ ∫ t, (f (Matrix.updateRow x0 k t)
            - ∫ s, f (Matrix.updateRow x0 k s) ∂(Measure.pi fun _ : Fin D => ν)) ^ 2
            ∂(Measure.pi fun _ : Fin D => ν) := fun x0 =>
      integral_nonneg (μ := Measure.pi fun _ : Fin D => ν) fun t => sq_nonneg _
    have hmono := integral_mono_of_nonneg (Filter.Eventually.of_forall hnn)
      (integrable_const (μ := noiseMatrix ν P D) (mixRowBound ν z D (y k)))
      (Filter.Eventually.of_forall (hinner k))
    refine hmono.trans (le_of_eq ?_)
    rw [integral_const, hone, one_smul]
  refine le_trans (Finset.sum_le_sum fun k _ => hterm k) ?_
  exact sum_mixRowBound_le ν z hz hy


/-! ### Milestone 4: the mean -/

/-- A complex linear form of one row has mean zero. -/
private theorem integral_clinForm_eq_zero {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (a : Fin D → ℂ) :
    ∫ t, (∑ l, a l * (t l : ℂ)) ∂(Measure.pi fun _ : Fin D => ν) = 0 := by
  have hprob := hν.prob
  set ar : Fin D → ℝ := fun l => (a l).re with har
  set ai : Fin D → ℝ := fun l => (a l).im with hai
  have hdecomp : ∀ t : Fin D → ℝ, (∑ l, a l * (t l : ℂ))
      = ((LinForm.linForm ar t : ℝ) : ℂ)
        + Complex.I * ((LinForm.linForm ai t : ℝ) : ℂ) := by
    intro t
    have hre : (∑ l, a l * (t l : ℂ)).re = LinForm.linForm ar t := by
      rw [Complex.re_sum]
      exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_re, har]
    have him : (∑ l, a l * (t l : ℂ)).im = LinForm.linForm ai t := by
      rw [Complex.im_sum]
      exact Finset.sum_congr rfl fun l _ => by simp [Complex.mul_im, hai]
    apply Complex.ext <;> simp [hre, him]
  have hIr : Integrable (fun t : Fin D → ℝ => LinForm.linForm ar t)
      (Measure.pi fun _ : Fin D => ν) := by
    have h := LinForm.integrable_linForm_pow_of_four (a := ar)
      (LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 ar) (k := 1) (by norm_num)
    simpa using h
  have hIi : Integrable (fun t : Fin D → ℝ => LinForm.linForm ai t)
      (Measure.pi fun _ : Fin D => ν) := by
    have h := LinForm.integrable_linForm_pow_of_four (a := ai)
      (LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 ai) (k := 1) (by norm_num)
    simpa using h
  have hIr' : Integrable (fun t : Fin D → ℝ => ((LinForm.linForm ar t : ℝ) : ℂ))
      (Measure.pi fun _ : Fin D => ν) := hIr.ofReal
  have hIi' : Integrable (fun t : Fin D → ℝ =>
      Complex.I * ((LinForm.linForm ai t : ℝ) : ℂ))
      (Measure.pi fun _ : Fin D => ν) := hIi.ofReal.const_mul Complex.I
  rw [integral_congr_ae (Filter.Eventually.of_forall hdecomp),
    integral_add hIr' hIi', integral_const_mul, integral_complex_ofReal,
    integral_complex_ofReal, LinForm.integral_linForm hν.mean hν.var hν.mom4 ar,
    LinForm.integral_linForm hν.mean hν.var hν.mom4 ai]
  simp

/-- The row form of the mixed statement: `xᵀ G_k g k`. -/
noncomputable def rowForm (z : ℂ) (x : Fin d → ℝ) (Y : Matrix (Fin p) (Fin d) ℝ)
    (k : Fin p) : ℂ :=
  R4C.cformC (gram (Matrix.updateRow Y k 0)) z x (rowVec Y k)

theorem measurable_rowForm {P D : ℕ} (z : ℂ) (x : Fin D → ℝ) (k : Fin P) :
    Measurable fun Y : Matrix (Fin P) (Fin D) ℝ => rowForm z x Y k := by
  simp only [rowForm]
  exact measurable_cformC (measurable_gram_loo k) z (fun i => measurable_const)
    (measurable_rowVec k)

/-- The mixed form is the `y`-combination of `(1 - s k) * rowForm`. -/
theorem mixF_eq_sum_rowForm (Y : Matrix (Fin p) (Fin d) ℝ) (hz : 0 < z.im) (x : Fin d → ℝ)
    (y : Fin p → ℝ) :
    mixF z x y Y = ∑ k, (y k : ℂ) * ((1 - sRow Y z k) * rowForm z x Y k) := by
  rw [mixF, cformC_mixed_eq_sum]
  exact Finset.sum_congr rfl fun k _ => by rw [cformC_rowVec_right_eq Y hz x k]; rfl

/-- The row form has an `L²` bound `1/(D η²)` and mean zero. -/
private theorem lintegral_normSq_rowForm_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1) (k : Fin (n + 1)) :
    ∫⁻ Y, ENNReal.ofReal (‖rowForm z x Y k‖ ^ 2) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal (1 / ((D : ℝ) * z.im ^ 2)) := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hsqi : ((Real.sqrt D)⁻¹ : ℝ) ^ 2 = ((D : ℝ))⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hDR.le]
  refine lintegral_noiseMatrix_le
    (ENNReal.measurable_ofReal.comp ((measurable_rowForm z x k).norm.pow_const 2)) k fun r => ?_
  set X : Matrix (Fin (n + 1)) (Fin D) ℝ :=
    Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r with hX
  set aa : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec x) l with haa
  have hcoef : ∑ l, ‖aa l‖ ^ 2 ≤ 1 / ((D : ℝ) * z.im ^ 2) := by
    have hstep : ∀ l, ‖aa l‖ ^ 2
        = ((D : ℝ))⁻¹ * ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec x) l‖ ^ 2 := by
      intro l
      rw [haa]
      simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs, hsqi]
    rw [Finset.sum_congr rfl fun l _ => hstep l, ← Finset.mul_sum]
    have h := sum_normSq_resolvC_mulVec_le (gram_isHermitian X) hz x
    calc ((D : ℝ))⁻¹ * ∑ l, ‖(R4C.resolvC (gram X) z *ᵥ R4C.cvec x) l‖ ^ 2
        ≤ ((D : ℝ))⁻¹ * ((x ⬝ᵥ x) / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left h (by positivity)
      _ ≤ ((D : ℝ))⁻¹ * (1 / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left (div_le_div_of_nonneg_right hx (by positivity))
            (by positivity)
      _ = 1 / ((D : ℝ) * z.im ^ 2) := by field_simp
  have hpt : ∀ t : Fin D → ℝ,
      rowForm z x (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k
        = ∑ l, aa l * (t l : ℂ) := by
    intro t
    simp only [rowForm]
    rw [updateRow_insertNth k r t, ← hX, rowVec_insertNth k r t,
      cformC_comm (gram_isHermitian X) hz.ne', cformC_smul_left_eq]
  calc ∫⁻ t, ENNReal.ofReal (‖rowForm z x
        (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k‖ ^ 2)
        ∂(Measure.pi fun _ : Fin D => ν)
      = ∫⁻ t, ENNReal.ofReal (‖∑ l, aa l * (t l : ℂ)‖ ^ 2)
        ∂(Measure.pi fun _ : Fin D => ν) := lintegral_congr fun t => by rw [hpt t]
    _ = ENNReal.ofReal (∑ l, ‖aa l‖ ^ 2) := lintegral_normSq_clinForm hν aa
    _ ≤ ENNReal.ofReal (1 / ((D : ℝ) * z.im ^ 2)) := ENNReal.ofReal_le_ofReal hcoef

private theorem integral_rowForm_eq_zero {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (_hx : x ⬝ᵥ x ≤ 1) (k : Fin (n + 1))
    (hint : Integrable (fun Y => rowForm z x Y k) (noiseMatrix ν (n + 1) D)) :
    ∫ Y, rowForm z x Y k ∂(noiseMatrix ν (n + 1) D) = 0 := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  rw [integral_noiseMatrix_split hint k]
  refine integral_eq_zero_of_ae (Filter.Eventually.of_forall fun r => ?_)
  set X : Matrix (Fin (n + 1)) (Fin D) ℝ :=
    Fin.insertNth (α := fun _ : Fin (n + 1) => (Fin D → ℝ)) k (0 : Fin D → ℝ) r with hX
  set aa : Fin D → ℂ := fun l =>
    (((Real.sqrt D)⁻¹ : ℝ) : ℂ) * (R4C.resolvC (gram X) z *ᵥ R4C.cvec x) l with haa
  have hpt : ∀ t : Fin D → ℝ,
      rowForm z x (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k
        = ∑ l, aa l * (t l : ℂ) := by
    intro t
    simp only [rowForm]
    rw [updateRow_insertNth k r t, ← hX, rowVec_insertNth k r t,
      cformC_comm (gram_isHermitian X) hz.ne', cformC_smul_left_eq]
  calc ∫ t, rowForm z x (Fin.insertNth k t r : Matrix (Fin (n + 1)) (Fin D) ℝ) k
        ∂(Measure.pi fun _ : Fin D => ν)
      = ∫ t, (∑ l, aa l * (t l : ℂ)) ∂(Measure.pi fun _ : Fin D => ν) :=
        integral_congr_ae (Filter.Eventually.of_forall hpt)
    _ = 0 := integral_clinForm_eq_zero hν aa


/-- The mean error of one row, `O(1/D)`. -/
noncomputable def mixMeanBound (ν : Measure ℝ) (z : ℂ) (P D : ℕ) : ℝ :=
  1 / 2 * ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z P D) + 1 / 2 * (1 / ((D : ℝ) * z.im ^ 2))

theorem mixMeanBound_nonneg (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im) (P D : ℕ) :
    0 ≤ mixMeanBound ν z P D := by
  have h1 : (0 : ℝ) ≤ (‖z‖ / z.im) ^ 4 * alphaSqBound ν z P D :=
    mul_nonneg (by positivity) (alphaSqBound_nonneg ν hz P D)
  have h2 : (0 : ℝ) ≤ 1 / ((D : ℝ) * z.im ^ 2) := by positivity
  rw [mixMeanBound]
  linarith

/-- **The mean of one row term.** The linear form has conditional mean zero, so the
deterministic `secW` may be subtracted for free, and what is left is `O(1/D)`. -/
private theorem norm_integral_mixRow_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1) (k : Fin (n + 1)) :
    ‖∫ Y, (1 - sRow Y z k) * rowForm z x Y k ∂(noiseMatrix ν (n + 1) D)‖
      ≤ mixMeanBound ν z (n + 1) D := by
  have hprob := hν.prob
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  set c : ℝ := ((n + 1 : ℕ) : ℝ) / D with hc
  have hcpos : 0 < c := by rw [hc]; positivity
  set w : ℂ := secW c z with hw
  set F2 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => ‖rowForm z x Y k‖ ^ 2 with hF2
  have hF2m : Measurable F2 := (measurable_rowForm z x k).norm.pow_const 2
  have hF2nn : ∀ Y, 0 ≤ F2 Y := fun Y => by rw [hF2]; positivity
  have hL2 := lintegral_normSq_rowForm_le hν hz hD hx k
  have hF2int : Integrable F2 (noiseMatrix ν (n + 1) D) :=
    integrable_of_lintegral_ofReal_ne_top hF2m hF2nn
      (ne_top_of_le_ne_top ENNReal.ofReal_ne_top hL2)
  have hF2bd : ∫ Y, F2 Y ∂(noiseMatrix ν (n + 1) D) ≤ 1 / ((D : ℝ) * z.im ^ 2) :=
    integral_le_of_lintegral_le hF2int hF2nn (by positivity) hL2
  have hdom : Integrable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => 1 / 2 * (F2 Y + 1))
      (noiseMatrix ν (n + 1) D) := (hF2int.add (integrable_const 1)).const_mul (1 / 2)
  have hrnorm : Integrable (fun Y => ‖rowForm z x Y k‖) (noiseMatrix ν (n + 1) D) := by
    refine Integrable.mono' hdom (measurable_rowForm z x k).norm.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
    have hF2v : F2 Y = ‖rowForm z x Y k‖ ^ 2 := rfl
    rw [hF2v]
    nlinarith [sq_nonneg (‖rowForm z x Y k‖ - 1), norm_nonneg (rowForm z x Y k)]
  have hrint : Integrable (fun Y => rowForm z x Y k) (noiseMatrix ν (n + 1) D) :=
    Integrable.mono' hrnorm (measurable_rowForm z x k).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => le_rfl)
  set F1 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => ‖(1 - sRow Y z k) - w‖ ^ 2 with hF1
  have hdevm : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => (1 - sRow Y z k) - w :=
    (measurable_const.sub
      (measurable_qformC measurable_gram_self z (measurable_rowVec k))).sub measurable_const
  have hF1m : Measurable F1 := hdevm.norm.pow_const 2
  have hF1nn : ∀ Y, 0 ≤ F1 Y := fun Y => by rw [hF1]; positivity
  have hKnn : (0 : ℝ) ≤ (‖z‖ / z.im) ^ 4 := by positivity
  have hAnn : (0 : ℝ) ≤ alphaSqBound ν z (n + 1) D := alphaSqBound_nonneg ν hz (n + 1) D
  have hL1 : ∫⁻ Y, ENNReal.ofReal (F1 Y) ∂(noiseMatrix ν (n + 1) D)
      ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D) := by
    have hpt : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
        ENNReal.ofReal (F1 Y)
          ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * ‖alphaRow Y z k - MP.mC c z‖ ^ 2) := by
      intro Y
      refine ENNReal.ofReal_le_ofReal ?_
      have h := norm_sub_secW_le hcpos hz hD Y k
      have h2 := pow_le_pow_left₀ (norm_nonneg ((1 - sRow Y z k) - w)) h 2
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
      _ ≤ ENNReal.ofReal ((‖z‖ / z.im) ^ 4)
            * ENNReal.ofReal (alphaSqBound ν z (n + 1) D) := by
          gcongr
          exact lintegral_normSq_alphaRow_sub_mC_le hν hz hD k
      _ = ENNReal.ofReal ((‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D) :=
          (ENNReal.ofReal_mul hKnn).symm
  have hF1int : Integrable F1 (noiseMatrix ν (n + 1) D) :=
    integrable_of_lintegral_ofReal_ne_top hF1m hF1nn
      (ne_top_of_le_ne_top ENNReal.ofReal_ne_top hL1)
  have hF1bd : ∫ Y, F1 Y ∂(noiseMatrix ν (n + 1) D)
      ≤ (‖z‖ / z.im) ^ 4 * alphaSqBound ν z (n + 1) D :=
    integral_le_of_lintegral_le hF1int hF1nn (mul_nonneg hKnn hAnn) hL1
  set G : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => 1 / 2 * F1 Y + 1 / 2 * F2 Y with hG
  have hGint : Integrable G (noiseMatrix ν (n + 1) D) :=
    (hF1int.const_mul _).add (hF2int.const_mul _)
  have hpbd : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ,
      ‖((1 - sRow Y z k) - w) * rowForm z x Y k‖ ≤ G Y := by
    intro Y
    rw [norm_mul, hG, hF1, hF2]
    nlinarith [sq_nonneg (‖(1 - sRow Y z k) - w‖ - ‖rowForm z x Y k‖),
      norm_nonneg ((1 - sRow Y z k) - w), norm_nonneg (rowForm z x Y k)]
  have hpint : Integrable (fun Y => ((1 - sRow Y z k) - w) * rowForm z x Y k)
      (noiseMatrix ν (n + 1) D) :=
    Integrable.mono' hGint (hdevm.mul (measurable_rowForm z x k)).aestronglyMeasurable
      (Filter.Eventually.of_forall hpbd)
  have hsplitint : ∫ Y, (1 - sRow Y z k) * rowForm z x Y k ∂(noiseMatrix ν (n + 1) D)
      = ∫ Y, ((1 - sRow Y z k) - w) * rowForm z x Y k ∂(noiseMatrix ν (n + 1) D) := by
    have hfun : (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ =>
        (1 - sRow Y z k) * rowForm z x Y k)
        = fun Y => ((1 - sRow Y z k) - w) * rowForm z x Y k + w * rowForm z x Y k := by
      funext Y; ring
    rw [hfun, integral_add hpint (hrint.const_mul w), integral_const_mul,
      integral_rowForm_eq_zero hν hz hD hx k hrint, mul_zero, add_zero]
  rw [hsplitint]
  calc ‖∫ Y, ((1 - sRow Y z k) - w) * rowForm z x Y k ∂(noiseMatrix ν (n + 1) D)‖
      ≤ ∫ Y, ‖((1 - sRow Y z k) - w) * rowForm z x Y k‖ ∂(noiseMatrix ν (n + 1) D) :=
        norm_integral_le_integral_norm _
    _ ≤ ∫ Y, G Y ∂(noiseMatrix ν (n + 1) D) := integral_mono hpint.norm hGint hpbd
    _ = 1 / 2 * (∫ Y, F1 Y ∂(noiseMatrix ν (n + 1) D))
          + 1 / 2 * ∫ Y, F2 Y ∂(noiseMatrix ν (n + 1) D) := by
        rw [hG, integral_add (hF1int.const_mul _) (hF2int.const_mul _), integral_const_mul,
          integral_const_mul]
    _ ≤ mixMeanBound ν z (n + 1) D := by
        rw [mixMeanBound]
        linarith


private theorem integrable_norm_rowForm {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1) (k : Fin (n + 1)) :
    Integrable (fun Y => ‖rowForm z x Y k‖) (noiseMatrix ν (n + 1) D) := by
  have hprob := hν.prob
  set F2 : Matrix (Fin (n + 1)) (Fin D) ℝ → ℝ := fun Y => ‖rowForm z x Y k‖ ^ 2 with hF2
  have hF2m : Measurable F2 := (measurable_rowForm z x k).norm.pow_const 2
  have hF2nn : ∀ Y, 0 ≤ F2 Y := fun Y => by rw [hF2]; positivity
  have hF2int : Integrable F2 (noiseMatrix ν (n + 1) D) :=
    integrable_of_lintegral_ofReal_ne_top hF2m hF2nn
      (ne_top_of_le_ne_top ENNReal.ofReal_ne_top (lintegral_normSq_rowForm_le hν hz hD hx k))
  have hdom : Integrable (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => 1 / 2 * (F2 Y + 1))
      (noiseMatrix ν (n + 1) D) := (hF2int.add (integrable_const 1)).const_mul (1 / 2)
  refine Integrable.mono' hdom (measurable_rowForm z x k).norm.aestronglyMeasurable
    (Filter.Eventually.of_forall fun Y => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
  have hF2v : F2 Y = ‖rowForm z x Y k‖ ^ 2 := rfl
  rw [hF2v]
  nlinarith [sq_nonneg (‖rowForm z x Y k‖ - 1), norm_nonneg (rowForm z x Y k)]

private theorem integrable_mixRow {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1) (k : Fin (n + 1)) :
    Integrable (fun Y => (1 - sRow Y z k) * rowForm z x Y k) (noiseMatrix ν (n + 1) D) := by
  have hprob := hν.prob
  have hdevm : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => 1 - sRow Y z k :=
    measurable_const.sub (measurable_qformC measurable_gram_self z (measurable_rowVec k))
  refine Integrable.mono' ((integrable_norm_rowForm hν hz hD hx k).const_mul (‖z‖ / z.im))
    (hdevm.mul (measurable_rowForm z x k)).aestronglyMeasurable
    (Filter.Eventually.of_forall fun Y => ?_)
  rw [norm_mul]
  exact mul_le_mul_of_nonneg_right (norm_one_sub_sRow_le Y hz hD k) (norm_nonneg _)

/-- **The rate of the mixed bound.** It mentions `ν`, `z`, `p` and `d` and nothing else: no
direction and no `ε`. -/
noncomputable def isoRateMixed (ν : Measure ℝ) (z : ℂ) (p d : ℕ) : ℝ :=
  2 * mixVarBound ν z p d + (p : ℝ) * (mixMeanBound ν z p d) ^ 2

/-- **The mean of the mixed form.** -/
private theorem normSq_integral_mixF_le {n D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hz : 0 < z.im) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin (n + 1) → ℝ} (hy : y ⬝ᵥ y ≤ 1) :
    ‖∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)‖ ^ 2
      ≤ ((n + 1 : ℕ) : ℝ) * (mixMeanBound ν z (n + 1) D) ^ 2 := by
  have hprob := hν.prob
  have hMnn : (0 : ℝ) ≤ mixMeanBound ν z (n + 1) D := mixMeanBound_nonneg ν hz (n + 1) D
  have hint : ∀ k, Integrable (fun Y => ((y k : ℝ) : ℂ) * ((1 - sRow Y z k) * rowForm z x Y k))
      (noiseMatrix ν (n + 1) D) := fun k => (integrable_mixRow hν hz hD hx k).const_mul _
  have hEq : ∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)
      = ∑ k, ((y k : ℝ) : ℂ)
        * ∫ Y, (1 - sRow Y z k) * rowForm z x Y k ∂(noiseMatrix ν (n + 1) D) := by
    rw [integral_congr_ae (Filter.Eventually.of_forall fun Y => mixF_eq_sum_rowForm Y hz x y),
      integral_finsetSum _ fun k _ => hint k]
    exact Finset.sum_congr rfl fun k _ => integral_const_mul _ _
  have hbd : ‖∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)‖
      ≤ (∑ k, |y k|) * mixMeanBound ν z (n + 1) D := by
    rw [hEq]
    refine (norm_sum_le _ _).trans ?_
    rw [Finset.sum_mul]
    refine Finset.sum_le_sum fun k _ => ?_
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
    exact mul_le_mul_of_nonneg_left (norm_integral_mixRow_le hν hz hD hx k) (abs_nonneg _)
  have hcard : (∑ k, |y k|) ^ 2 ≤ ((n + 1 : ℕ) : ℝ) * ∑ k, |y k| ^ 2 := by
    have h := sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset (Fin (n + 1))))
      (f := fun k => |y k|)
    rwa [Finset.card_univ, Fintype.card_fin] at h
  have hysum : ∑ k, |y k| ^ 2 = y ⬝ᵥ y := by
    rw [dotProduct]
    exact Finset.sum_congr rfl fun k _ => by rw [sq_abs]; ring
  rw [hysum] at hcard
  have hsn : (0 : ℝ) ≤ ∑ k, |y k| := Finset.sum_nonneg fun k _ => abs_nonneg _
  have hsq := pow_le_pow_left₀ (norm_nonneg _) hbd 2
  refine hsq.trans ?_
  rw [mul_pow]
  have h2 : (∑ k, |y k|) ^ 2 ≤ ((n + 1 : ℕ) : ℝ) := by
    refine hcard.trans ?_
    nlinarith [Nat.cast_nonneg (α := ℝ) (n + 1)]
  exact mul_le_mul_of_nonneg_right h2 (sq_nonneg _)

/-- **The second moment of the mixed form.** -/
theorem integral_normSq_mixF_le {P D : ℕ} {ν : Measure ℝ} (hν : NoiseLaw ν) (hz : 0 < z.im)
    (hP : 0 < P) (hD : 0 < D) {x : Fin D → ℝ} (hx : x ⬝ᵥ x ≤ 1)
    {y : Fin P → ℝ} (hy : y ⬝ᵥ y ≤ 1) :
    ∫ Y, ‖mixF z x y Y‖ ^ 2 ∂(noiseMatrix ν P D) ≤ isoRateMixed ν z P D := by
  obtain ⟨n, rfl⟩ : ∃ n, P = n + 1 := ⟨P - 1, (Nat.succ_pred_eq_of_pos hP).symm⟩
  have hprob := hν.prob
  set Bd : ℝ := Real.sqrt (‖z‖ / z.im ^ 2) with hBd
  have hFm : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => mixF z x y Y :=
    measurable_mixF z x y
  have hFbd : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, ‖mixF z x y Y‖ ≤ Bd := fun Y =>
    norm_mixF_le hz hD hx hy Y
  have hFint : Integrable (fun Y => mixF z x y Y) (noiseMatrix ν (n + 1) D) :=
    Integrable.mono' (integrable_const (μ := noiseMatrix ν (n + 1) D) Bd)
      hFm.aestronglyMeasurable (Filter.Eventually.of_forall hFbd)
  have hcomp : ∀ (φ : ℂ → ℝ), Measurable φ → (∀ w : ℂ, |φ w| ≤ ‖w‖) →
      Integrable (fun Y => φ (mixF z x y Y)) (noiseMatrix ν (n + 1) D)
        ∧ Integrable (fun Y => φ (mixF z x y Y) ^ 2) (noiseMatrix ν (n + 1) D) := by
    intro φ hφm hφn
    have hm2 : Measurable fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => φ (mixF z x y Y) :=
      hφm.comp hFm
    have hbd : ∀ Y : Matrix (Fin (n + 1)) (Fin D) ℝ, |φ (mixF z x y Y)| ≤ Bd := fun Y =>
      (hφn _).trans (hFbd Y)
    refine ⟨Integrable.mono' (integrable_const (μ := noiseMatrix ν (n + 1) D) Bd)
      hm2.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by rw [Real.norm_eq_abs]; exact hbd Y), ?_⟩
    refine Integrable.mono' (integrable_const (μ := noiseMatrix ν (n + 1) D) (Bd ^ 2))
      (hm2.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), ← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) (hbd Y) 2
  obtain ⟨hReI, hReS⟩ := hcomp Complex.re Complex.measurable_re Complex.abs_re_le_norm
  obtain ⟨hImI, hImS⟩ := hcomp Complex.im Complex.measurable_im Complex.abs_im_le_norm
  have hsplit : ∫ Y, ‖mixF z x y Y‖ ^ 2 ∂(noiseMatrix ν (n + 1) D)
      = (∫ Y, (mixF z x y Y).re ^ 2 ∂(noiseMatrix ν (n + 1) D))
        + ∫ Y, (mixF z x y Y).im ^ 2 ∂(noiseMatrix ν (n + 1) D) := by
    have hfun : (fun Y : Matrix (Fin (n + 1)) (Fin D) ℝ => ‖mixF z x y Y‖ ^ 2)
        = fun Y => (mixF z x y Y).re ^ 2 + (mixF z x y Y).im ^ 2 := by
      funext Y
      rw [Complex.sq_norm, Complex.normSq_apply]; ring
    rw [hfun, integral_add hReS hImS]
  have hvRe := variance_component_mix_le hν hz hD hx hy Complex.re Complex.measurable_re
    (fun w v => (Complex.sub_re w v).symm) Complex.abs_re_le_norm
  have hvIm := variance_component_mix_le hν hz hD hx hy Complex.im Complex.measurable_im
    (fun w v => (Complex.sub_im w v).symm) Complex.abs_im_le_norm
  have hmeanEq : (∫ Y, (mixF z x y Y).re ∂(noiseMatrix ν (n + 1) D)) ^ 2
      + (∫ Y, (mixF z x y Y).im ∂(noiseMatrix ν (n + 1) D)) ^ 2
      = ‖∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)‖ ^ 2 := by
    have hre : (∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)).re
        = ∫ Y, (mixF z x y Y).re ∂(noiseMatrix ν (n + 1) D) := (integral_re hFint).symm
    have him : (∫ Y, mixF z x y Y ∂(noiseMatrix ν (n + 1) D)).im
        = ∫ Y, (mixF z x y Y).im ∂(noiseMatrix ν (n + 1) D) := (integral_im hFint).symm
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
    ring
  have hmean := normSq_integral_mixF_le hν hz hD hx hy
  rw [hsplit, integral_sq_eq_var_add hReI hReS, integral_sq_eq_var_add hImI hImS,
    isoRateMixed]
  linarith [hvRe, hvIm, hmean, hmeanEq]

/-- **Statement 5: the mixed form.** For `x` in the unit ball of the right space and `y` in the
unit ball of the left space, `xᵀ G Eᵀ y` is within `ε` of `0` outside a set of measure at most
`isoRateMixed ν z p d / ε²`, and the rate is free of `x`, `y` and `ε`. -/
theorem measure_cformC_mixed_ge_le {ν : Measure ℝ} (hν : NoiseLaw ν) {z : ℂ} (hz : 0 < z.im)
    {p d : ℕ} (hp : 0 < p) (hd : 0 < d) (x : Fin d → ℝ) (hx : x ⬝ᵥ x ≤ 1)
    (y : Fin p → ℝ) (hy : y ⬝ᵥ y ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    noiseMatrix ν p d {Y | ε ≤ ‖R4C.cformC (gram Y) z x ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))‖}
      ≤ ENNReal.ofReal (isoRateMixed ν z p d / ε ^ 2) := by
  have hprob := hν.prob
  have hFm : Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => ‖mixF z x y Y‖ :=
    (measurable_mixF z x y).norm
  have hFbd : ∀ Y : Matrix (Fin p) (Fin d) ℝ,
      ‖mixF z x y Y‖ ≤ Real.sqrt (‖z‖ / z.im ^ 2) := fun Y => norm_mixF_le hz hd hx hy Y
  have hint : Integrable (fun Y : Matrix (Fin p) (Fin d) ℝ => ‖mixF z x y Y‖ ^ 2)
      (noiseMatrix ν p d) := by
    refine Integrable.mono' (integrable_const (μ := noiseMatrix ν p d)
      (Real.sqrt (‖z‖ / z.im ^ 2) ^ 2)) (hFm.pow_const 2).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact pow_le_pow_left₀ (norm_nonneg _) (hFbd Y) 2
  have hset : {Y : Matrix (Fin p) (Fin d) ℝ
      | ε ≤ ‖R4C.cformC (gram Y) z x ((Real.sqrt d)⁻¹ • (Yᵀ *ᵥ y))‖}
      = {Y | ε ≤ ‖mixF z x y Y‖} := rfl
  rw [hset]
  exact Cheb.meas_ge_le_of_integral_sq _ hFm hint hε
    (integral_normSq_mixF_le hν hz hp hd hx hy)

/-! ### The rate vanishes -/

/-- **The mixed rate vanishes.** Same shape as `tendsto_isoRate` (`RMT/General/Iso.lean`). -/
theorem tendsto_isoRateMixed {pN dN : ℕ → ℕ} (ν : Measure ℝ) {z : ℂ} (hz : 0 < z.im) {c : ℝ}
    (hd : Tendsto dN atTop atTop)
    (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c)) :
    Tendsto (fun N => isoRateMixed ν z (pN N) (dN N)) atTop (𝓝 0) := by
  have hzi : z.im ≠ 0 := hz.ne'
  have hev : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  have ht : Tendsto (fun N => ((dN N : ℝ))⁻¹) atTop (𝓝 0) :=
    (tendsto_natCast_atTop_atTop.comp hd).inv_tendsto_atTop
  have hsmall : ∀ a : ℝ, Tendsto (fun N => a / ((dN N : ℝ) * z.im ^ 2)) atTop (𝓝 0) := by
    intro a
    have h := ht.const_mul (a / z.im ^ 2)
    rw [mul_zero] at h
    refine h.congr fun N => ?_
    rcases eq_or_ne ((dN N : ℝ)) 0 with h0 | h0
    · rw [h0]; simp
    · field_simp
  -- `D * residConst` converges
  have hresD : Tendsto (fun N => ((dN N : ℝ)) * residConst ν (dN N) z) atTop
      (𝓝 (2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2)) := by
    have hnice : Tendsto (fun N : ℕ => 2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2
        + 2 / ((dN N : ℝ) * z.im ^ 2)) atTop
        (𝓝 (2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2 + 0)) :=
      tendsto_const_nhds.add (hsmall 2)
    rw [add_zero] at hnice
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [residConst]
    field_simp
  have hcoefT : Tendsto (fun N => traceCoef z ((pN N : ℝ) / dN N)) atTop
      (𝓝 (traceCoef z c)) := tendsto_traceCoef hz hcN
  have htraceD : Tendsto (fun N => ((dN N : ℝ)) * traceSqBound ν z (pN N) (dN N)) atTop
      (𝓝 (traceCoef z c * (2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2))) := by
    refine (hcoefT.mul hresD).congr fun N => ?_
    rw [traceSqBound_eq]
    ring
  have halphaD : Tendsto (fun N => ((dN N : ℝ)) * alphaSqBound ν z (pN N) (dN N)) atTop
      (𝓝 (3 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2
        + 3 * (traceCoef z c * (2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2)))) := by
    have h1 : Tendsto (fun N : ℕ => 3 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2
        + 3 / ((dN N : ℝ) * z.im ^ 2)) atTop
        (𝓝 (3 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2 + 0)) :=
      tendsto_const_nhds.add (hsmall 3)
    rw [add_zero] at h1
    have hnice := h1.add (htraceD.const_mul 3)
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    set T : ℝ := traceSqBound ν z (pN N) (dN N) with hT
    rw [alphaSqBound, ← hT]
    field_simp
    ring
  -- `D * mixMeanBound` converges
  have hmeanD : Tendsto (fun N => ((dN N : ℝ)) * mixMeanBound ν z (pN N) (dN N)) atTop
      (𝓝 (1 / 2 * ((‖z‖ / z.im) ^ 4 * (3 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2
          + 3 * (traceCoef z c * (2 * ((∫ w, w ^ 4 ∂ν) + 2) / z.im ^ 2))))
        + 1 / 2 * (1 / z.im ^ 2))) := by
    have hnice := ((halphaD.const_mul ((‖z‖ / z.im) ^ 4)).const_mul (1 / 2 : ℝ)).add
      (tendsto_const_nhds (x := (1 / 2 : ℝ) * (1 / z.im ^ 2)) (f := atTop (α := ℕ)))
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    set A : ℝ := alphaSqBound ν z (pN N) (dN N) with hA
    rw [mixMeanBound, ← hA]
    field_simp
  -- the two halves
  have hvar : Tendsto (fun N => mixVarBound ν z (pN N) (dN N)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => 2 * (‖z‖ / z.im) ^ 2 / z.im ^ 2 * ((dN N : ℝ))⁻¹)
        atTop (𝓝 0) := by
      have h := ht.const_mul (2 * (‖z‖ / z.im) ^ 2 / z.im ^ 2)
      rw [mul_zero] at h
      exact h
    have h2 : Tendsto (fun N => ((pN N : ℝ) / dN N) * (((dN N : ℝ))⁻¹
        * (2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2) / z.im ^ 4)))
        atTop (𝓝 0) := by
      have h := hcN.mul (ht.mul (tendsto_const_nhds
        (x := 2 * (‖z‖ / z.im) ^ 2 * ((∫ w, w ^ 4 ∂ν) + 3) * (1 + ‖z‖ ^ 2) / z.im ^ 4)))
      rw [zero_mul, mul_zero] at h
      exact h
    have hnice := h1.add h2
    rw [add_zero] at hnice
    refine Tendsto.congr' ?_ hnice
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    rw [mixVarBound]
    field_simp
  have hmean : Tendsto (fun N => ((pN N : ℝ)) * (mixMeanBound ν z (pN N) (dN N)) ^ 2)
      atTop (𝓝 0) := by
    have h := (hcN.mul ht).mul (hmeanD.mul hmeanD)
    rw [mul_zero, zero_mul] at h
    refine Tendsto.congr' ?_ h
    filter_upwards [hev] with N hN
    have hDR : (0 : ℝ) < (dN N : ℝ) := by exact_mod_cast hN
    set M : ℝ := mixMeanBound ν z (pN N) (dN N) with hM
    field_simp
  have hall := (hvar.const_mul (2 : ℝ)).add hmean
  rw [mul_zero, add_zero] at hall
  exact hall.congr fun N => by rw [isoRateMixed]

end GenRMT

end StackedSVD
