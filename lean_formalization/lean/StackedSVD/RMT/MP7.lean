/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.MP

/-!
# Item MP7: the holomorphic branch of the MP transform on the upper half plane

Specification: `notes/archive/rmt_MP.md` (MP6), `notes/archive/rmt_R1.md` step 6 and modeling choice
5, `notes/archive/rmt_T.md` (T3). Numeric confirmation: two session scripts, not kept
(`check_mp7.py`, seed 20260829, and `check_mp7b.py`, seed 20260831).

For a fixed `c > 0` this file builds `mC c : ℂ → ℂ`, the branch of the MP Stieltjes
transform with `Im > 0` on the upper half plane, and proves the four facts that items R1
and T need.

1. `quad_mC`, `im_mC_pos`: for `0 < z.im` the value `mC c z` is a root of `MP.quad c z`
   with positive imaginary part, hence *the* root of `MP.existsUnique_root_im_pos`
   (`mC_eq_root`).
2. `hasDerivAt_mC`: `HasDerivAt (mC c) (mCDeriv c z) z`, with the closed form
   `mCDeriv c z = -(mC c z * (mC c z + 1)) / sqrtDisc c z`, the same shape as the real
   `MP.mDeriv`. `mCDeriv_eq_secular` gives the form `-(m² + m)/(2 z m + z + 1 - c)` that
   `notes/archive/rmt_R1.md` writes.
3. `norm_mC_le : ‖mC c z‖ ≤ (z.im)⁻¹` and `norm_mCDeriv_le : ‖mCDeriv c z‖ ≤ (z.im)⁻¹ ^ 2`,
   the bounds R1b assumes. Both constants are `1`, and both are sharp (the grid maximum is
   0.999999999 at `z = 10000 i`).
4. `tendsto_mC`, `tendsto_mCDeriv`: continuity at the real axis from above, above the edge.

The square root is `csqrt w = exp (log w / 2)` applied to **each factor** of the
discriminant, `sqrtDisc c z = csqrt (z - b) * csqrt (z - b')`. The product form
`exp (log ((z-b)(z-b')) / 2)` of `notes/archive/rmt_T.md` is *not* the branch with `Im > 0`: at
`c = 1`, `z = i` it returns `-1.3002 - 0.6248 i` while the root with `Im > 0` is
`0.3002 + 0.6248 i` (7776 of 20000 grid points fail). The factored form is correct because
`Im z > 0` puts each factor in the upper half plane, so each square root lies in the open
first quadrant and the product has `Im > 0`.
-/

open Filter Topology

namespace StackedSVD
namespace MP

/-! ### A holomorphic square root on the slit plane -/

/-- Principal square root, `exp (log w / 2)`. Junk at `w = 0`, where it returns `1`. -/
noncomputable def csqrt (w : ℂ) : ℂ := Complex.exp (Complex.log w / 2)

theorem csqrt_ne_zero (w : ℂ) : csqrt w ≠ 0 := Complex.exp_ne_zero _

theorem csqrt_mul_self {w : ℂ} (hw : w ≠ 0) : csqrt w * csqrt w = w := by
  unfold csqrt
  rw [← Complex.exp_add, show Complex.log w / 2 + Complex.log w / 2 = Complex.log w by ring]
  exact Complex.exp_log hw

theorem csqrt_sq {w : ℂ} (hw : w ≠ 0) : csqrt w ^ 2 = w := by
  rw [sq]; exact csqrt_mul_self hw

theorem im_log_div_two (w : ℂ) : (Complex.log w / 2).im = Complex.arg w / 2 := by
  rw [Complex.div_im, Complex.log_im, Complex.log_re]
  norm_num [Complex.normSq_apply]
  ring

/-- On the upper half plane the argument is in `(0, π)`. -/
theorem arg_pos_of_im_pos {w : ℂ} (hw : 0 < w.im) : 0 < Complex.arg w := by
  rcases lt_or_eq_of_le (Complex.arg_nonneg_iff.mpr hw.le) with h | h
  · exact h
  · exact absurd (Complex.arg_eq_zero_iff.mp h.symm).2 hw.ne'

theorem arg_lt_pi_of_im_pos {w : ℂ} (hw : 0 < w.im) : Complex.arg w < Real.pi :=
  Complex.arg_lt_pi_iff.mpr (Or.inr hw.ne')

/-- `Im w > 0` puts `csqrt w` in the open first quadrant. -/
theorem re_csqrt_pos {w : ℂ} (hw : 0 < w.im) : 0 < (csqrt w).re := by
  have h1 := arg_pos_of_im_pos hw
  have h2 := arg_lt_pi_of_im_pos hw
  unfold csqrt
  rw [Complex.exp_re, im_log_div_two]
  refine mul_pos (Real.exp_pos _) (Real.cos_pos_of_mem_Ioo ⟨?_, ?_⟩)
  · have := Real.pi_pos; linarith
  · linarith

theorem im_csqrt_pos {w : ℂ} (hw : 0 < w.im) : 0 < (csqrt w).im := by
  have h1 := arg_pos_of_im_pos hw
  have h2 := arg_lt_pi_of_im_pos hw
  unfold csqrt
  rw [Complex.exp_im, im_log_div_two]
  refine mul_pos (Real.exp_pos _) (Real.sin_pos_of_pos_of_lt_pi (by linarith) ?_)
  have := Real.pi_pos; linarith

theorem csqrt_ofReal {r : ℝ} (hr : 0 < r) : csqrt (r : ℂ) = ((Real.sqrt r : ℝ) : ℂ) := by
  unfold csqrt
  rw [← Complex.ofReal_log hr.le,
    show ((Real.log r : ℝ) : ℂ) / 2 = (((Real.log r / 2 : ℝ)) : ℂ) by push_cast; ring,
    ← Complex.ofReal_exp]
  congr 1
  rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hr]
  congr 1
  ring

theorem hasDerivAt_csqrt {w : ℂ} (hw : w ∈ Complex.slitPlane) :
    HasDerivAt csqrt (csqrt w / (2 * w)) w := by
  have hne : w ≠ 0 := by
    intro h
    rw [h] at hw
    exact Complex.zero_notMem_slitPlane hw
  have h2 : HasDerivAt (fun u : ℂ => Complex.log u / 2) (w⁻¹ / 2) w :=
    (Complex.hasDerivAt_log hw).div_const 2
  have h3 : HasDerivAt (fun u : ℂ => Complex.exp (Complex.log u / 2))
      (Complex.exp (Complex.log w / 2) * (w⁻¹ / 2)) w := h2.cexp
  have hval : csqrt w / (2 * w) = Complex.exp (Complex.log w / 2) * (w⁻¹ / 2) := by
    unfold csqrt
    field_simp
  rw [hval]
  exact h3

/-- Chain rule for `csqrt`, in the form that avoids `HasDerivAt.comp` (whose `ℂ`-over-`ℂ`
module instance does not match the one `HasDerivAt` picks here). -/
theorem hasDerivAt_csqrt_comp {f : ℂ → ℂ} {f' z : ℂ} (hf : HasDerivAt f f' z)
    (h : f z ∈ Complex.slitPlane) :
    HasDerivAt (fun w : ℂ => csqrt (f w)) (csqrt (f z) * (f' / f z / 2)) z :=
  ((hf.clog h).div_const 2).cexp

theorem hasDerivAt_csqrt_sub {z r : ℂ} (h : z - r ∈ Complex.slitPlane) :
    HasDerivAt (fun w : ℂ => csqrt (w - r)) (csqrt (z - r) / (2 * (z - r))) z := by
  have hne : z - r ≠ 0 := by
    intro hzero
    rw [hzero] at h
    exact Complex.zero_notMem_slitPlane h
  have hf := hasDerivAt_csqrt_comp ((hasDerivAt_id' (x := z)).sub_const r) h
  have heq : csqrt (z - r) / (2 * (z - r)) = csqrt (z - r) * (1 / (z - r) / 2) := by
    field_simp
  rw [heq]
  exact hf

/-! ### The discriminant and its square root -/

/-- `√((z - b) (z - b'))`, as the product of the two principal square roots. -/
noncomputable def sqrtDisc (c : ℝ) (z : ℂ) : ℂ :=
  csqrt (z - (bulkEdge c : ℂ)) * csqrt (z - (bulkEdgeLo c : ℂ))

section Disc

variable {c : ℝ} {z : ℂ}

theorem bulkEdge_add_bulkEdgeLo (hc : 0 ≤ c) : bulkEdge c + bulkEdgeLo c = 2 + 2 * c := by
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 ≤ s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_nonneg c, (Real.sq_sqrt hc).symm⟩
  unfold bulkEdge bulkEdgeLo
  rw [Real.sqrt_sq hs0]
  ring

/-- MP2 over `ℂ`. -/
theorem discr_eqC (hc : 0 ≤ c) (z : ℂ) :
    (z + 1 - (c : ℂ)) ^ 2 - 4 * z = (z - (bulkEdge c : ℂ)) * (z - (bulkEdgeLo c : ℂ)) := by
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 ≤ s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_nonneg c, (Real.sq_sqrt hc).symm⟩
  unfold bulkEdge bulkEdgeLo
  rw [Real.sqrt_sq hs0]
  push_cast
  ring

theorem sub_bulkEdge_ne_zero (hz : 0 < z.im) : z - (bulkEdge c : ℂ) ≠ 0 := by
  intro h
  have h2 : (z - (bulkEdge c : ℂ)).im = 0 := by rw [h]; simp
  simp only [Complex.sub_im, Complex.ofReal_im, sub_zero] at h2
  linarith

theorem sub_bulkEdgeLo_ne_zero (hz : 0 < z.im) : z - (bulkEdgeLo c : ℂ) ≠ 0 := by
  intro h
  have h2 : (z - (bulkEdgeLo c : ℂ)).im = 0 := by rw [h]; simp
  simp only [Complex.sub_im, Complex.ofReal_im, sub_zero] at h2
  linarith

theorem sqrtDisc_ne_zero : sqrtDisc c z ≠ 0 :=
  mul_ne_zero (csqrt_ne_zero _) (csqrt_ne_zero _)

theorem sqrtDisc_sq (hc : 0 ≤ c) (hz : 0 < z.im) :
    sqrtDisc c z ^ 2 = (z + 1 - (c : ℂ)) ^ 2 - 4 * z := by
  rw [sqrtDisc, mul_pow, csqrt_sq (sub_bulkEdge_ne_zero hz),
    csqrt_sq (sub_bulkEdgeLo_ne_zero hz), discr_eqC hc]

/-- The branch fact: `Im √((z-b)(z-b')) ≥ Im z` on the upper half plane. Both square roots
lie in the open first quadrant and `2 (Re s) (Im s) = Im z` for each of them, so the claim
is AM-GM. -/
theorem im_le_im_sqrtDisc (hz : 0 < z.im) : z.im ≤ (sqrtDisc c z).im := by
  have hb : 0 < (z - (bulkEdge c : ℂ)).im := by
    simp only [Complex.sub_im, Complex.ofReal_im, sub_zero]; exact hz
  have hb' : 0 < (z - (bulkEdgeLo c : ℂ)).im := by
    simp only [Complex.sub_im, Complex.ofReal_im, sub_zero]; exact hz
  have hne1 := sub_bulkEdge_ne_zero (c := c) hz
  have hne2 := sub_bulkEdgeLo_ne_zero (c := c) hz
  have e1 : 2 * ((csqrt (z - (bulkEdge c : ℂ))).re * (csqrt (z - (bulkEdge c : ℂ))).im)
      = z.im := by
    have h := congrArg Complex.im (csqrt_mul_self hne1)
    simp only [Complex.mul_im, Complex.sub_im, Complex.ofReal_im, sub_zero] at h
    linear_combination h
  have e2 : 2 * ((csqrt (z - (bulkEdgeLo c : ℂ))).re * (csqrt (z - (bulkEdgeLo c : ℂ))).im)
      = z.im := by
    have h := congrArg Complex.im (csqrt_mul_self hne2)
    simp only [Complex.mul_im, Complex.sub_im, Complex.ofReal_im, sub_zero] at h
    linear_combination h
  have h1r : 0 < (csqrt (z - (bulkEdge c : ℂ))).re := re_csqrt_pos hb
  have h1i : 0 < (csqrt (z - (bulkEdge c : ℂ))).im := im_csqrt_pos hb
  have h2r : 0 < (csqrt (z - (bulkEdgeLo c : ℂ))).re := re_csqrt_pos hb'
  have h2i : 0 < (csqrt (z - (bulkEdgeLo c : ℂ))).im := im_csqrt_pos hb'
  rw [sqrtDisc, Complex.mul_im]
  set a := (csqrt (z - (bulkEdge c : ℂ))).re
  set b := (csqrt (z - (bulkEdge c : ℂ))).im
  set p := (csqrt (z - (bulkEdgeLo c : ℂ))).re
  set q := (csqrt (z - (bulkEdgeLo c : ℂ))).im
  have hprod : (2 * (a * b)) * (2 * (p * q)) = z.im * z.im := by rw [e1, e2]
  have hsum : 0 < a * q + b * p := add_pos (mul_pos h1r h2i) (mul_pos h1i h2r)
  have key : z.im * z.im ≤ (a * q + b * p) ^ 2 := by
    nlinarith [sq_nonneg (a * q - b * p), hprod]
  nlinarith [key, hsum, hz]

theorem im_sqrtDisc_pos (hz : 0 < z.im) : 0 < (sqrtDisc c z).im :=
  lt_of_lt_of_le hz (im_le_im_sqrtDisc hz)

/-- The algebra behind `hasDerivAt_sqrtDisc`: `s₁' s₂ + s₁ s₂' = (z-1-c)/(s₁ s₂)`. -/
private theorem sqrtDisc_deriv_alg {A1 A2 z b1 b2 cc : ℂ} (h1 : A1 * A1 = z - b1)
    (h2 : A2 * A2 = z - b2) (hA1 : A1 ≠ 0) (hA2 : A2 ≠ 0) (hb : b1 + b2 = 2 + 2 * cc) :
    (z - 1 - cc) / (A1 * A2)
      = A1 / (2 * (z - b1)) * A2 + A1 * (A2 / (2 * (z - b2))) := by
  have key : A1 / (2 * (A1 * A1)) * A2 + A1 * (A2 / (2 * (A2 * A2)))
      = (A1 * A1 + A2 * A2) / (2 * (A1 * A2)) := by
    field_simp
    ring
  rw [← h1, ← h2, key, h1, h2, div_eq_div_iff (mul_ne_zero hA1 hA2)
    (mul_ne_zero (two_ne_zero) (mul_ne_zero hA1 hA2))]
  linear_combination (A1 * A2) * hb

theorem hasDerivAt_sqrtDisc (hc : 0 ≤ c) (hz : 0 < z.im) :
    HasDerivAt (sqrtDisc c) ((z - 1 - (c : ℂ)) / sqrtDisc c z) z := by
  have hb : 0 < (z - (bulkEdge c : ℂ)).im := by
    simp only [Complex.sub_im, Complex.ofReal_im, sub_zero]; exact hz
  have hb' : 0 < (z - (bulkEdgeLo c : ℂ)).im := by
    simp only [Complex.sub_im, Complex.ofReal_im, sub_zero]; exact hz
  have hne1 := sub_bulkEdge_ne_zero (c := c) hz
  have hne2 := sub_bulkEdgeLo_ne_zero (c := c) hz
  have hs1 : (z - (bulkEdge c : ℂ)) ∈ Complex.slitPlane :=
    Complex.mem_slitPlane_iff.mpr (Or.inr hb.ne')
  have hs2 : (z - (bulkEdgeLo c : ℂ)) ∈ Complex.slitPlane :=
    Complex.mem_slitPlane_iff.mpr (Or.inr hb'.ne')
  have d1 := hasDerivAt_csqrt_sub (z := z) (r := (bulkEdge c : ℂ)) hs1
  have d2 := hasDerivAt_csqrt_sub (z := z) (r := (bulkEdgeLo c : ℂ)) hs2
  have hsumC : (bulkEdge c : ℂ) + (bulkEdgeLo c : ℂ) = 2 + 2 * (c : ℂ) := by
    have h := bulkEdge_add_bulkEdgeLo hc
    have h2 : ((bulkEdge c + bulkEdgeLo c : ℝ) : ℂ) = ((2 + 2 * c : ℝ) : ℂ) := by rw [h]
    push_cast at h2
    exact h2
  have hval : (z - 1 - (c : ℂ)) / sqrtDisc c z
      = csqrt (z - (bulkEdge c : ℂ)) / (2 * (z - (bulkEdge c : ℂ)))
          * csqrt (z - (bulkEdgeLo c : ℂ))
        + csqrt (z - (bulkEdge c : ℂ))
          * (csqrt (z - (bulkEdgeLo c : ℂ)) / (2 * (z - (bulkEdgeLo c : ℂ)))) := by
    rw [sqrtDisc]
    exact sqrtDisc_deriv_alg (csqrt_mul_self hne1) (csqrt_mul_self hne2)
      (csqrt_ne_zero _) (csqrt_ne_zero _) hsumC
  rw [hval]
  exact d1.mul d2

end Disc

/-! ### The branch `mC` and its derivative -/

/-- The MP transform on the upper half plane, the root of `MP.quad c z` with `Im > 0`.
Junk outside the upper half plane: the same formula, with `csqrt` on the principal branch,
and the value `0` at `z = 0`. -/
noncomputable def mC (c : ℝ) (z : ℂ) : ℂ := (-(z + 1 - (c : ℂ)) + sqrtDisc c z) / (2 * z)

/-- The derivative of `mC`, in closed form; the complex twin of `MP.mDeriv`. -/
noncomputable def mCDeriv (c : ℝ) (z : ℂ) : ℂ := -(mC c z * (mC c z + 1)) / sqrtDisc c z

section Branch

variable {c : ℝ} {z : ℂ}

theorem ne_zero_of_im_pos (hz : 0 < z.im) : z ≠ 0 := by
  intro h; rw [h] at hz; simp at hz

/-- `mC c z` solves the MP quadratic. -/
theorem quad_mC (hc : 0 ≤ c) (hz : 0 < z.im) : quad c z (mC c z) = 0 := by
  change quad c z ((-(z + 1 - (c : ℂ)) + sqrtDisc c z) / (2 * z)) = 0
  exact quad_root_of_sq (ne_zero_of_im_pos hz) (sqrtDisc_sq hc hz)

/-- `2 z m + z + 1 - c = √((z-b)(z-b'))`, the complex twin of `MP.two_mul_mul_m_add`. -/
theorem two_mul_mul_mC_add (hz : z ≠ 0) :
    2 * z * mC c z + z + 1 - (c : ℂ) = sqrtDisc c z := by
  unfold mC
  field_simp
  ring

theorem im_add_pos (hz : 0 < z.im) : 0 < (sqrtDisc c z + (z + 1 - (c : ℂ))).im := by
  have h := im_sqrtDisc_pos (c := c) hz
  simp only [Complex.add_im, Complex.sub_im, Complex.one_im, Complex.ofReal_im, add_zero,
    sub_zero]
  linarith

theorem add_ne_zero_of_im_pos (hz : 0 < z.im) : sqrtDisc c z + (z + 1 - (c : ℂ)) ≠ 0 := by
  intro h
  have h2 := im_add_pos (c := c) hz
  rw [h] at h2
  simp at h2

/-- Rationalized form, the complex twin of `MP.m_eq_neg_two_div`. -/
theorem mC_eq (hc : 0 ≤ c) (hz : 0 < z.im) :
    mC c z = -2 / (sqrtDisc c z + (z + 1 - (c : ℂ))) := by
  have h2z : (2 : ℂ) * z ≠ 0 := mul_ne_zero two_ne_zero (ne_zero_of_im_pos hz)
  unfold mC
  rw [div_eq_div_iff h2z (add_ne_zero_of_im_pos (c := c) hz)]
  linear_combination sqrtDisc_sq hc hz

/-- The branch has positive imaginary part. -/
theorem im_mC_pos (hc : 0 ≤ c) (hz : 0 < z.im) : 0 < (mC c z).im := by
  have hWim := im_add_pos (c := c) hz
  have hWne := add_ne_zero_of_im_pos (c := c) hz
  have hns : 0 < Complex.normSq (sqrtDisc c z + (z + 1 - (c : ℂ))) :=
    Complex.normSq_pos.mpr hWne
  have hinv : ((sqrtDisc c z + (z + 1 - (c : ℂ)))⁻¹).im < 0 := by
    rw [Complex.inv_im]
    exact div_neg_of_neg_of_pos (by linarith) hns
  rw [mC_eq hc hz, div_eq_mul_inv, Complex.mul_im]
  simp only [Complex.neg_re, Complex.neg_im, Complex.re_ofNat, Complex.im_ofNat, neg_zero,
    zero_mul, add_zero]
  linarith

/-- `mC c z` is the root of `MP.existsUnique_root_im_pos`. -/
theorem mC_eq_root (hc : 0 < c) (hz : 0 < z.im) {u : ℂ} (hu : 0 < u.im)
    (hq : quad c z u = 0) : u = mC c z :=
  (existsUnique_root_im_pos hc hz).unique ⟨hu, hq⟩ ⟨im_mC_pos hc.le hz, quad_mC hc.le hz⟩

/-- The form of the derivative that `notes/archive/rmt_R1.md` writes. -/
theorem mCDeriv_eq_secular (hz : 0 < z.im) :
    mCDeriv c z = -(mC c z * (mC c z + 1)) / (2 * z * mC c z + z + 1 - (c : ℂ)) := by
  rw [mCDeriv, two_mul_mul_mC_add (ne_zero_of_im_pos hz)]

theorem mCDeriv_eq_div (hc : 0 ≤ c) (hz : 0 < z.im) :
    mCDeriv c z = 2 * (sqrtDisc c z + z - 1 - (c : ℂ))
      / (sqrtDisc c z * ((sqrtDisc c z + z - 1 - (c : ℂ)) + 2) ^ 2) := by
  have hSne : sqrtDisc c z ≠ 0 := sqrtDisc_ne_zero
  have hWne := add_ne_zero_of_im_pos (c := c) hz
  have hWT : sqrtDisc c z + (z + 1 - (c : ℂ)) = (sqrtDisc c z + z - 1 - (c : ℂ)) + 2 := by ring
  have hmC : mC c z = -2 / ((sqrtDisc c z + z - 1 - (c : ℂ)) + 2) := by
    rw [mC_eq hc hz, hWT]
  rw [hWT] at hWne
  unfold mCDeriv
  rw [hmC]
  field_simp
  ring

theorem hasDerivAt_mC (hc : 0 ≤ c) (hz : 0 < z.im) :
    HasDerivAt (mC c) (mCDeriv c z) z := by
  have hSne : sqrtDisc c z ≠ 0 := sqrtDisc_ne_zero
  have hWne := add_ne_zero_of_im_pos (c := c) hz
  have hd1 : HasDerivAt (fun w : ℂ => sqrtDisc c w + (w + 1 - (c : ℂ)))
      ((z - 1 - (c : ℂ)) / sqrtDisc c z + 1) z :=
    (hasDerivAt_sqrtDisc hc hz).add
      (((hasDerivAt_id' (x := z)).add_const (1 : ℂ)).sub_const ((c : ℂ)))
  have hd2 : HasDerivAt (fun w : ℂ => (-2 : ℂ) / (sqrtDisc c w + (w + 1 - (c : ℂ))))
      ((0 * (sqrtDisc c z + (z + 1 - (c : ℂ)))
        - (-2 : ℂ) * ((z - 1 - (c : ℂ)) / sqrtDisc c z + 1))
        / (sqrtDisc c z + (z + 1 - (c : ℂ))) ^ 2) z :=
    (hasDerivAt_const z (-2 : ℂ)).div hd1 hWne
  have heq : mC c =ᶠ[𝓝 z] fun w : ℂ => (-2 : ℂ) / (sqrtDisc c w + (w + 1 - (c : ℂ))) := by
    have hopen : IsOpen {w : ℂ | 0 < w.im} := isOpen_lt continuous_const Complex.continuous_im
    filter_upwards [hopen.mem_nhds hz] with w hw
    exact mC_eq hc hw
  have hval : mCDeriv c z
      = (0 * (sqrtDisc c z + (z + 1 - (c : ℂ)))
        - (-2 : ℂ) * ((z - 1 - (c : ℂ)) / sqrtDisc c z + 1))
        / (sqrtDisc c z + (z + 1 - (c : ℂ))) ^ 2 := by
    unfold mCDeriv
    rw [mC_eq hc hz]
    field_simp
    ring
  rw [hval]
  exact hd2.congr_of_eventuallyEq heq

theorem deriv_mC (hc : 0 ≤ c) (hz : 0 < z.im) : deriv (mC c) z = mCDeriv c z :=
  (hasDerivAt_mC hc hz).deriv

end Branch

/-! ### The two norm bounds -/

/-- Scalar core of `norm_mCDeriv_le`. With `nt = ‖T‖`, `ns = ‖S‖`, `nu = ‖T+2‖`,
`ti = Im T` and `η = Im z`. -/
private theorem norm_deriv_scalar {η nt ns nu ti cc : ℝ} (hη : 0 < η) (hnt : 0 < nt)
    (hns : 0 ≤ ns) (h1 : 2 * η * nt ^ 2 = ti * (nt ^ 2 - 4 * cc))
    (h2 : nt ^ 2 ≤ 2 * ns * nt + 4 * cc) (h3 : ti ≤ nu) (h4 : 2 * η ≤ ti) :
    2 * nt * η ^ 2 ≤ ns * nu ^ 2 := by
  have hti : 0 < ti := by linarith
  have hkey : nt ^ 2 - 4 * cc ≤ 2 * ns * nt := by linarith
  have hstep : η * nt ≤ ti * ns := by
    nlinarith [mul_le_mul_of_nonneg_left hkey hti.le, hnt]
  have hc2 : 2 * η * (η * nt) ≤ 2 * η * (ti * ns) :=
    mul_le_mul_of_nonneg_left hstep (by linarith)
  have hc3 : 2 * η * (ti * ns) ≤ ti * (ti * ns) := by
    have h0 : 0 ≤ ti * ns := mul_nonneg hti.le hns
    nlinarith [h4, h0]
  have hc4 : ti * ti * ns ≤ nu * nu * ns :=
    mul_le_mul_of_nonneg_right (mul_self_le_mul_self hti.le h3) hns
  nlinarith [hc2, hc3, hc4]

section Bounds

variable {c : ℝ} {z : ℂ}

/-- `‖m(z)‖ ≤ 1 / Im z`, the Stieltjes bound. Sharp. -/
theorem norm_mC_le (hc : 0 ≤ c) (hz : 0 < z.im) : ‖mC c z‖ ≤ (z.im)⁻¹ := by
  have hSim := im_le_im_sqrtDisc (c := c) hz
  have hWim : 2 * z.im ≤ (sqrtDisc c z + (z + 1 - (c : ℂ))).im := by
    simp only [Complex.add_im, Complex.sub_im, Complex.one_im, Complex.ofReal_im, add_zero,
      sub_zero]
    linarith
  have hWnorm : 2 * z.im ≤ ‖sqrtDisc c z + (z + 1 - (c : ℂ))‖ :=
    le_trans hWim (Complex.im_le_norm _)
  have hWpos : 0 < ‖sqrtDisc c z + (z + 1 - (c : ℂ))‖ := by linarith
  have hnum : ‖mC c z‖ = 2 / ‖sqrtDisc c z + (z + 1 - (c : ℂ))‖ := by
    rw [mC_eq hc hz, norm_div, norm_neg, Complex.norm_two]
  rw [hnum, div_le_iff₀ hWpos, inv_mul_eq_div, le_div_iff₀ hz]
  exact hWnorm

/-- `‖m'(z)‖ ≤ (Im z)⁻²`, the bound R1b assumes. Sharp. -/
theorem norm_mCDeriv_le (hc : 0 < c) (hz : 0 < z.im) : ‖mCDeriv c z‖ ≤ (z.im)⁻¹ ^ 2 := by
  have hSne : sqrtDisc c z ≠ 0 := sqrtDisc_ne_zero
  have hSim := im_le_im_sqrtDisc (c := c) hz
  have hS2 := sqrtDisc_sq hc.le hz
  have hDeq := mCDeriv_eq_div hc.le hz
  -- the two algebraic identities in `T = S + z - 1 - c`
  have hidT : (sqrtDisc c z + z - 1 - (c : ℂ)) * (sqrtDisc c z + z - 1 - (c : ℂ))
      = sqrtDisc c z * (sqrtDisc c z + z - 1 - (c : ℂ))
        + sqrtDisc c z * (sqrtDisc c z + z - 1 - (c : ℂ)) + ((4 * c : ℝ) : ℂ) := by
    push_cast
    linear_combination -hS2
  have hidz : 2 * z * (sqrtDisc c z + z - 1 - (c : ℂ))
      = (sqrtDisc c z + z - 1 - (c : ℂ)) * (sqrtDisc c z + z - 1 - (c : ℂ))
        + 2 * (1 + (c : ℂ)) * (sqrtDisc c z + z - 1 - (c : ℂ)) + 4 * (c : ℂ) := by
    linear_combination -hS2
  have hTim : (sqrtDisc c z + z - 1 - (c : ℂ)).im = (sqrtDisc c z).im + z.im := by
    simp only [Complex.sub_im, Complex.add_im, Complex.one_im, Complex.ofReal_im, sub_zero]
  have hTim2 : 2 * z.im ≤ (sqrtDisc c z + z - 1 - (c : ℂ)).im := by rw [hTim]; linarith
  set T := sqrtDisc c z + z - 1 - (c : ℂ) with hTdef
  have hTpos : 0 < T.im := by linarith
  have hTne : T ≠ 0 := by intro h; rw [h] at hTpos; simp at hTpos
  have hntpos : 0 < ‖T‖ := norm_pos_iff.mpr hTne
  have hT2im : (T + 2).im = T.im := by simp
  have hT2ne : T + 2 ≠ 0 := by
    intro h
    rw [h] at hT2im
    simp only [Complex.zero_im] at hT2im
    linarith
  have hnupos : 0 < ‖T + 2‖ := norm_pos_iff.mpr hT2ne
  have hnspos : 0 < ‖sqrtDisc c z‖ := norm_pos_iff.mpr hSne
  -- h1: the real identity linking `Im z`, `Im T` and `‖T‖`
  have h1 : 2 * z.im * ‖T‖ ^ 2 = T.im * (‖T‖ ^ 2 - 4 * c) := by
    have hre := congrArg Complex.re hidz
    have him := congrArg Complex.im hidz
    simp only [Complex.mul_re, Complex.mul_im, Complex.add_re, Complex.add_im,
      Complex.ofReal_re, Complex.ofReal_im, Complex.one_re, Complex.one_im,
      Complex.re_ofNat, Complex.im_ofNat] at hre him
    rw [Complex.sq_norm, Complex.normSq_apply]
    linear_combination T.re * him - T.im * hre
  -- h2: the triangle inequality on `T * T = 2 S T + 4c`
  have h2 : ‖T‖ ^ 2 ≤ 2 * ‖sqrtDisc c z‖ * ‖T‖ + 4 * c := by
    have hsq : ‖T‖ ^ 2 = ‖T * T‖ := by rw [norm_mul]; ring
    rw [hsq, hidT]
    calc ‖sqrtDisc c z * T + sqrtDisc c z * T + ((4 * c : ℝ) : ℂ)‖
        ≤ ‖sqrtDisc c z * T + sqrtDisc c z * T‖ + ‖((4 * c : ℝ) : ℂ)‖ := norm_add_le _ _
      _ ≤ ‖sqrtDisc c z * T‖ + ‖sqrtDisc c z * T‖ + ‖((4 * c : ℝ) : ℂ)‖ := by
          gcongr
          exact norm_add_le _ _
      _ = 2 * ‖sqrtDisc c z‖ * ‖T‖ + 4 * c := by
          rw [norm_mul, Complex.norm_real, Real.norm_eq_abs,
            abs_of_pos (by linarith : (0 : ℝ) < 4 * c)]
          ring
  have h3 : T.im ≤ ‖T + 2‖ := by
    rw [← hT2im]; exact Complex.im_le_norm _
  -- assemble
  have hden : 0 < ‖sqrtDisc c z‖ * ‖T + 2‖ ^ 2 := by positivity
  have hz2 : 0 < z.im ^ 2 := by positivity
  have hnorm : ‖mCDeriv c z‖ = 2 * ‖T‖ / (‖sqrtDisc c z‖ * ‖T + 2‖ ^ 2) := by
    rw [hDeq, norm_div, norm_mul, norm_mul, norm_pow, Complex.norm_two]
  rw [hnorm, inv_pow, ← one_div, div_le_div_iff₀ hden hz2, one_mul]
  exact norm_deriv_scalar hz hntpos hnspos.le h1 h2 h3 hTim2

/-- The `1 / (Im z)^2` shape that `notes/archive/rmt_R1.md` R1b states. -/
theorem norm_mCDeriv_le' (hc : 0 < c) (hz : 0 < z.im) : ‖mCDeriv c z‖ ≤ 1 / z.im ^ 2 := by
  rw [one_div, ← inv_pow]
  exact norm_mCDeriv_le hc hz

end Bounds

/-! ### T3: the real axis from above -/

section RealAxis

variable {c x : ℝ}

theorem sqrtDisc_ofReal (hc : 0 < c) (hx : bulkEdge c < x) :
    sqrtDisc c (x : ℂ)
      = ((Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)) : ℝ) : ℂ) := by
  have h1 : 0 < x - bulkEdge c := by linarith
  have h2 : 0 < x - bulkEdgeLo c := by
    have := bulkEdgeLo_lt_bulkEdge hc; linarith
  have e1 : (x : ℂ) - (bulkEdge c : ℂ) = ((x - bulkEdge c : ℝ) : ℂ) := by push_cast; ring
  have e2 : (x : ℂ) - (bulkEdgeLo c : ℂ) = ((x - bulkEdgeLo c : ℝ) : ℂ) := by push_cast; ring
  rw [sqrtDisc, e1, e2, csqrt_ofReal h1, csqrt_ofReal h2, ← Complex.ofReal_mul,
    ← Real.sqrt_mul h1.le]

theorem mC_ofReal (hc : 0 < c) (hx : bulkEdge c < x) : mC c (x : ℂ) = ((m c x : ℝ) : ℂ) := by
  rw [mC, sqrtDisc_ofReal hc hx, m]
  push_cast
  ring

theorem mCDeriv_ofReal (hc : 0 < c) (hx : bulkEdge c < x) :
    mCDeriv c (x : ℂ) = ((mDeriv c x : ℝ) : ℂ) := by
  rw [mCDeriv, mC_ofReal hc hx, sqrtDisc_ofReal hc hx, mDeriv]
  push_cast
  ring

theorem continuousAt_sqrtDisc_ofReal (hc : 0 < c) (hx : bulkEdge c < x) :
    ContinuousAt (sqrtDisc c) (x : ℂ) := by
  have h1 : 0 < x - bulkEdge c := by linarith
  have h2 : 0 < x - bulkEdgeLo c := by
    have := bulkEdgeLo_lt_bulkEdge hc; linarith
  have hs1 : ((x : ℂ) - (bulkEdge c : ℂ)) ∈ Complex.slitPlane := by
    refine Complex.mem_slitPlane_iff.mpr (Or.inl ?_)
    simp only [Complex.sub_re, Complex.ofReal_re]
    linarith
  have hs2 : ((x : ℂ) - (bulkEdgeLo c : ℂ)) ∈ Complex.slitPlane := by
    refine Complex.mem_slitPlane_iff.mpr (Or.inl ?_)
    simp only [Complex.sub_re, Complex.ofReal_re]
    linarith
  have d1 := (hasDerivAt_csqrt_sub (z := (x : ℂ)) (r := (bulkEdge c : ℂ)) hs1).continuousAt
  have d2 := (hasDerivAt_csqrt_sub (z := (x : ℂ)) (r := (bulkEdgeLo c : ℂ)) hs2).continuousAt
  exact d1.mul d2

theorem continuousAt_mC_ofReal (hc : 0 < c) (hx : bulkEdge c < x) :
    ContinuousAt (mC c) (x : ℂ) := by
  have hx0 : (0 : ℝ) < x := lt_trans (bulkEdge_pos hc.le) hx
  have hden : (2 : ℂ) * (x : ℂ) ≠ 0 :=
    mul_ne_zero two_ne_zero (Complex.ofReal_ne_zero.mpr hx0.ne')
  have hnum : ContinuousAt (fun w : ℂ => -(w + 1 - (c : ℂ)) + sqrtDisc c w) (x : ℂ) :=
    (by fun_prop : ContinuousAt (fun w : ℂ => -(w + 1 - (c : ℂ))) (x : ℂ)).add
      (continuousAt_sqrtDisc_ofReal hc hx)
  exact hnum.div (by fun_prop) hden

theorem continuousAt_mCDeriv_ofReal (hc : 0 < c) (hx : bulkEdge c < x) :
    ContinuousAt (mCDeriv c) (x : ℂ) := by
  have hm := continuousAt_mC_ofReal hc hx
  exact ((hm.mul (hm.add continuousAt_const)).neg).div
    (continuousAt_sqrtDisc_ofReal hc hx) sqrtDisc_ne_zero

theorem tendsto_ofReal_add_mul_I (x : ℝ) :
    Tendsto (fun η : ℝ => (x : ℂ) + (η : ℂ) * Complex.I) (𝓝[>] (0 : ℝ)) (𝓝 ((x : ℂ))) := by
  have hcont : Continuous (fun η : ℝ => (x : ℂ) + (η : ℂ) * Complex.I) := by fun_prop
  have h0 : Tendsto (fun η : ℝ => (x : ℂ) + (η : ℂ) * Complex.I) (𝓝 (0 : ℝ))
      (𝓝 ((x : ℂ) + (((0 : ℝ) : ℂ)) * Complex.I)) := hcont.tendsto 0
  have h1 := h0.mono_left (nhdsWithin_le_nhds (s := Set.Ioi (0 : ℝ)))
  simpa using h1

/-- **T3**, scalar part: the branch converges to the real transform above the edge. -/
theorem tendsto_mC (hc : 0 < c) (hx : bulkEdge c < x) :
    Tendsto (fun η : ℝ => mC c ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 ((m c x : ℝ) : ℂ)) := by
  have h := (continuousAt_mC_ofReal hc hx).tendsto.comp (tendsto_ofReal_add_mul_I x)
  rw [mC_ofReal hc hx] at h
  exact h

/-- **T3**, derivative part. -/
theorem tendsto_mCDeriv (hc : 0 < c) (hx : bulkEdge c < x) :
    Tendsto (fun η : ℝ => mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 ((mDeriv c x : ℝ) : ℂ)) := by
  have h := (continuousAt_mCDeriv_ofReal hc hx).tendsto.comp (tendsto_ofReal_add_mul_I x)
  rw [mCDeriv_ofReal hc hx] at h
  exact h

/-- **T3**, derivative part in the `deriv` form that `notes/archive/L2_STATEMENTS.md` states. -/
theorem tendsto_deriv_mC (hc : 0 < c) (hx : bulkEdge c < x) :
    Tendsto (fun η : ℝ => deriv (mC c) ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 ((mDeriv c x : ℝ) : ℂ)) := by
  refine Filter.Tendsto.congr' ?_ (tendsto_mCDeriv hc hx)
  filter_upwards [self_mem_nhdsWithin] with η hη
  have hη' : 0 < ((x : ℂ) + (η : ℂ) * Complex.I).im := by
    simpa using Set.mem_Ioi.mp hη
  exact (deriv_mC hc.le hη').symm

end RealAxis

end MP
end StackedSVD
