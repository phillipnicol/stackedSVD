/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.ResolvDeriv
import StackedSVD.RMT.MP7

/-!
# Item G5: the companion Stieltjes transform

`notes/archive/prop_single_table_general.md` section 5, unit G5. This file gives the companion
Stieltjes transform `MP.mTildeC`, the transform of the `p x p` block `d⁻¹ Y Yᵀ` normalized
by `p`, next to `MP.mC`, the transform of the `d x d` block `d⁻¹ Yᵀ Y` normalized by `d`. The
two matrices share their nonzero eigenvalues, so the two traces differ only by the zero
eigenvalues; the difference is exact at every finite `p` and `d`
(`GenRMT.trace_resolvC_gram_sub_gramC`). Unit G3 runs its self-consistent equation on the
companion side and needs both this finite identity and the limit object; unit G4 states the
companion limit with the name `MP.mTildeC`.

Does not import `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` or any `Vendor/COLT83/`
file (choice 8 of the note).

## Content

1. `MP.mTildeC`, its closed form `mTildeC_eq`, the defining relation `c_mul_mTildeC`, and the
   companion quadratic `quad_mTildeC`.
2. The duality `c_mul_mTildeC_eq_mC_inv`: the companion transform at ratio `c` is the transform
   at the reciprocal ratio and the rescaled argument. `im_mTildeC_pos` and `norm_mTildeC_le`
   follow from it, through the private helper `sqrtDisc_dual`.
3. The finite identity `GenRMT.trace_resolvC_gram_sub_gramC` and its normalized form
   `GenRMT.stieltjesC_gram_sub_gramC`, from a private intertwining lemma `Y * gram Y =
   gramC Y * Y`. Unit G2 proves the public twin `GenRMT.cmat'_mul_resolvC_gram` in
   `RMT/General/Companion.lean` (wave 2, in flight, not read here); follow-up item F34 covers
   the dedup.
-/

open scoped Matrix

namespace StackedSVD

namespace MP

/-- **The companion Stieltjes transform.** `c m̃(z) = m(z) - (c - 1)/z`, the transform of the
`p × p` block `d⁻¹ Y Yᵀ` normalized by `p`, against that of the `d × d` block `d⁻¹ YᵀY`
normalized by `d`, with `c = lim p/d`. -/
noncomputable def mTildeC (c : ℝ) (z : ℂ) : ℂ := (MP.mC c z + (1 - (c : ℂ)) / z) / c

/-- The closed form, the twin of `MP.mC` (`RMT/MP7.lean:272`) with `c` and `1` exchanged in
the linear term and `2 c z` in the denominator. -/
theorem mTildeC_eq {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    mTildeC c z = (-(z + (c : ℂ) - 1) + sqrtDisc c z) / (2 * (c : ℂ) * z) := by
  have hzne : z ≠ 0 := ne_zero_of_im_pos hz
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  unfold mTildeC mC
  field_simp
  ring

/-- The defining relation, in the form the finite-`d` identity of section 2 takes. -/
theorem c_mul_mTildeC {c : ℝ} {z : ℂ} (hc : 0 < c) :
    (c : ℂ) * mTildeC c z = mC c z + (1 - (c : ℂ)) / z := by
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  unfold mTildeC
  field_simp

/-- **The companion quadratic.** `c z m̃² + (z + c - 1) m̃ + 1 = 0`, the twin of
`MP.quad_mC` (`RMT/MP7.lean:285`). -/
theorem quad_mTildeC {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    (c : ℂ) * z * mTildeC c z ^ 2 + (z + (c : ℂ) - 1) * mTildeC c z + 1 = 0 := by
  have hzne : z ≠ 0 := ne_zero_of_im_pos hz
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  have hmC0 : mC c z = (c : ℂ) * mTildeC c z - (1 - (c : ℂ)) / z := by
    have h := c_mul_mTildeC (c := c) (z := z) hc
    rw [h]; ring
  have hkey : (c : ℂ) * ((c : ℂ) * z * mTildeC c z ^ 2 + (z + (c : ℂ) - 1) * mTildeC c z + 1)
      = z * mC c z ^ 2 + (z + 1 - (c : ℂ)) * mC c z + 1 := by
    rw [hmC0]
    field_simp
    ring
  have hquad0 : z * mC c z ^ 2 + (z + 1 - (c : ℂ)) * mC c z + 1 = 0 := quad_mC hc.le hz
  rw [hquad0] at hkey
  exact (mul_eq_zero.mp hkey).resolve_left hcne

/-- **The duality of the discriminant.** `c * sqrtDisc c⁻¹ (z/c) = sqrtDisc c z`. Private
helper: the one fact `c_mul_mTildeC_eq_mC_inv` needs. Both sides square to the same real
quadratic in `z` and `c` (`sqrtDisc_sq` at `(c, z)` and at `(c⁻¹, z/c)`), and both sides have
positive imaginary part, so the difference of squares forces them equal. -/
private theorem sqrtDisc_dual {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    (c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ)) = sqrtDisc c z := by
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  have hzcim : 0 < (z / (c : ℂ)).im := by
    rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have hsq1 : sqrtDisc c z ^ 2 = (z + 1 - (c : ℂ)) ^ 2 - 4 * z := sqrtDisc_sq hc.le hz
  have hsq2 : sqrtDisc c⁻¹ (z / (c : ℂ)) ^ 2
      = (z / (c : ℂ) + 1 - ((c⁻¹ : ℝ) : ℂ)) ^ 2 - 4 * (z / (c : ℂ)) :=
    sqrtDisc_sq (inv_nonneg.mpr hc.le) hzcim
  have hxsq : ((c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ))) ^ 2 = sqrtDisc c z ^ 2 := by
    rw [mul_pow, hsq2, hsq1, Complex.ofReal_inv]
    field_simp
    ring
  have hyim : 0 < (sqrtDisc c z).im := im_sqrtDisc_pos hz
  have hxim0 : 0 < (sqrtDisc c⁻¹ (z / (c : ℂ))).im := im_sqrtDisc_pos hzcim
  have hximtot : 0 < ((c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ))).im := by
    rw [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero]
    exact mul_pos hc hxim0
  have hsumim : 0 < ((c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ)) + sqrtDisc c z).im := by
    rw [Complex.add_im]; linarith
  have hsumne : (c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ)) + sqrtDisc c z ≠ 0 := by
    intro h
    rw [h] at hsumim
    simp at hsumim
  have hfactor : ((c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ)) - sqrtDisc c z)
      * ((c : ℂ) * sqrtDisc c⁻¹ (z / (c : ℂ)) + sqrtDisc c z) = 0 := by
    linear_combination hxsq
  exact sub_eq_zero.mp ((mul_eq_zero.mp hfactor).resolve_right hsumne)

/-- **The duality.** `c m̃_c(z) = m_{1/c}(z/c)`: the companion transform at ratio `c` is the
transform at the reciprocal ratio and the rescaled argument. Everything the two bounds below
need comes from this identity. -/
theorem c_mul_mTildeC_eq_mC_inv {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    (c : ℂ) * mTildeC c z = mC c⁻¹ (z / (c : ℂ)) := by
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  have hzne : z ≠ 0 := ne_zero_of_im_pos hz
  have hzcim : 0 < (z / (c : ℂ)).im := by
    rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have hzcne : z / (c : ℂ) ≠ 0 := ne_zero_of_im_pos hzcim
  have hdual := sqrtDisc_dual hc hz
  simp only [mTildeC, mC, Complex.ofReal_inv]
  rw [← hdual]
  field_simp
  ring

theorem im_mTildeC_pos {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    0 < (mTildeC c z).im := by
  have hzcim : 0 < (z / (c : ℂ)).im := by
    rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have h := im_mC_pos (inv_nonneg.mpr hc.le) hzcim
  rw [← c_mul_mTildeC_eq_mC_inv hc hz, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
    zero_mul, add_zero] at h
  exact (mul_pos_iff_of_pos_left hc).mp h

theorem norm_mTildeC_le {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    ‖mTildeC c z‖ ≤ (z.im)⁻¹ := by
  have hzcim : 0 < (z / (c : ℂ)).im := by
    rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have h := norm_mC_le (inv_nonneg.mpr hc.le) hzcim
  rw [← c_mul_mTildeC_eq_mC_inv hc hz, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos hc, Complex.div_ofReal_im, inv_div] at h
  have h2 : c * ‖mTildeC c z‖ ≤ c * (z.im)⁻¹ := by
    rw [← div_eq_mul_inv]; exact h
  exact le_of_mul_le_mul_left h2 hc

end MP

namespace GenRMT

/-! ### The finite identity: private helpers

Two small facts about `R4C.cmat'` that `RMT/General/Defs.lean` does not carry (it only has
`cmat'_apply`, `cmat'_transpose` and `cmat'_eq_cmat`): `cmat'` is multiplicative on rectangular
products and commutes with a real scalar, in the shape of `R4C.cmat_mul` (`RMT/R4C.lean:171`).
F37 (2026-09-09): the multiplicative fact is now `R4C.cmat'_mul`, the public canonical copy of
`Companion.lean:49` (this file's own private twin is dropped); only the scalar fact stays
private here, since it has no twin. -/

private theorem cmat'_smul {p q : ℕ} (r : ℝ) (A : Matrix (Fin p) (Fin q) ℝ) :
    R4C.cmat' (r • A) = (r : ℂ) • R4C.cmat' A := by
  ext i j
  simp only [R4C.cmat'_apply, Matrix.smul_apply, smul_eq_mul]
  push_cast
  ring

/-! F37 (2026-09-09): the intertwining `cmat' Y * G = Gc * cmat' Y` used to be repeated here
as `cmat'_mul_resolvC_gram_priv` (`private`). `GenRMT.cmat'_mul_resolvC_gram`
(`Companion.lean:299`, already public) is the canonical copy; this file uses it directly.
Follow-up item F34 tracked this dedup and is closed by it. -/

/-- `r * (r⁻¹ * tr M) = tr M`, for any `r : ℕ` (including `r = 0`, where the index type is
empty and `tr M = 0` too). Private helper for `stieltjesC_gram_sub_gramC`, so that the
statement carries no positivity hypothesis on `p`. -/
private theorem trace_cancel {r : ℕ} (M : Matrix (Fin r) (Fin r) ℂ) :
    (r : ℂ) * ((r : ℂ)⁻¹ * M.trace) = M.trace := by
  rcases Nat.eq_zero_or_pos r with hr | hr
  · subst hr
    simp
  · rw [← mul_assoc, mul_inv_cancel₀ (Nat.cast_ne_zero.mpr hr.ne'), one_mul]

/-- **The finite identity.** The two Gram matrices share their nonzero spectrum, so the two
traces differ by the zero eigenvalues only. Exact at every `p` and `d`. -/
theorem trace_resolvC_gram_sub_gramC {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) {z : ℂ}
    (hz : 0 < z.im) (hd : 0 < d) :
    (R4C.resolvC (gram Y) z).trace - (R4C.resolvC (gramC Y) z).trace = ((p : ℂ) - d) / z := by
  have _hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hzne : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hint : R4C.cmat' Y * R4C.resolvC (gram Y) z
      = R4C.resolvC (gramC Y) z * R4C.cmat' Y := cmat'_mul_resolvC_gram Y hz
  have hcgram : R4C.cmat (gram Y) = ((d : ℂ))⁻¹ • ((R4C.cmat' Y)ᵀ * R4C.cmat' Y) := by
    ext i j
    simp only [R4C.cmat, Matrix.map_apply, gram, Matrix.smul_apply, Matrix.mul_apply,
      Matrix.transpose_apply, R4C.cmat'_apply, smul_eq_mul]
    push_cast
    ring
  have hcgramC : R4C.cmat (gramC Y) = ((d : ℂ))⁻¹ • (R4C.cmat' Y * (R4C.cmat' Y)ᵀ) := by
    ext i j
    simp only [R4C.cmat, Matrix.map_apply, gramC, Matrix.smul_apply, Matrix.mul_apply,
      Matrix.transpose_apply, R4C.cmat'_apply, smul_eq_mul]
    push_cast
    ring
  have hAeq : R4C.cmat (gram Y) * R4C.resolvC (gram Y) z
      = 1 + z • R4C.resolvC (gram Y) z := by
    have hA := ResolvDeriv.cmat_sub_mul_resolvC (gram_isHermitian Y) hz.ne'
    rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul] at hA
    exact sub_eq_iff_eq_add.mp hA
  have hBeq : R4C.cmat (gramC Y) * R4C.resolvC (gramC Y) z
      = 1 + z • R4C.resolvC (gramC Y) z := by
    have hB := ResolvDeriv.cmat_sub_mul_resolvC (gramC_isHermitian Y) hz.ne'
    rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul] at hB
    exact sub_eq_iff_eq_add.mp hB
  have htrG : (R4C.cmat (gram Y) * R4C.resolvC (gram Y) z).trace
      = (d : ℂ) + z * (R4C.resolvC (gram Y) z).trace := by
    rw [hAeq, Matrix.trace_add, Matrix.trace_smul, smul_eq_mul, Matrix.trace_one,
      Fintype.card_fin]
  have htrGc : (R4C.cmat (gramC Y) * R4C.resolvC (gramC Y) z).trace
      = (p : ℂ) + z * (R4C.resolvC (gramC Y) z).trace := by
    rw [hBeq, Matrix.trace_add, Matrix.trace_smul, smul_eq_mul, Matrix.trace_one,
      Fintype.card_fin]
  have hsame : (R4C.cmat (gram Y) * R4C.resolvC (gram Y) z).trace
      = (R4C.cmat (gramC Y) * R4C.resolvC (gramC Y) z).trace := by
    rw [hcgram, hcgramC, Matrix.smul_mul, Matrix.smul_mul, Matrix.trace_smul,
      Matrix.trace_smul, smul_eq_mul, smul_eq_mul]
    congr 1
    calc ((R4C.cmat' Y)ᵀ * R4C.cmat' Y * R4C.resolvC (gram Y) z).trace
        = ((R4C.cmat' Y)ᵀ * (R4C.cmat' Y * R4C.resolvC (gram Y) z)).trace := by
          rw [Matrix.mul_assoc]
      _ = ((R4C.cmat' Y)ᵀ * (R4C.resolvC (gramC Y) z * R4C.cmat' Y)).trace := by rw [hint]
      _ = ((R4C.cmat' Y)ᵀ * R4C.resolvC (gramC Y) z * R4C.cmat' Y).trace := by
          rw [Matrix.mul_assoc]
      _ = (R4C.cmat' Y * (R4C.cmat' Y)ᵀ * R4C.resolvC (gramC Y) z).trace :=
          Matrix.trace_mul_cycle _ _ _
  have hcore : (d : ℂ) + z * (R4C.resolvC (gram Y) z).trace
      = (p : ℂ) + z * (R4C.resolvC (gramC Y) z).trace := by
    rw [← htrG, ← htrGc, hsame]
  have hzG : z * ((R4C.resolvC (gram Y) z).trace - (R4C.resolvC (gramC Y) z).trace)
      = (p : ℂ) - d := by linear_combination hcore
  field_simp
  linear_combination hzG

/-- The same identity on the normalized traces, the form unit G3 divides by `d`. `R4C.stieltjesC`
normalizes by the size of its own argument, so `R4C.stieltjesC (gram Y) z = (d:ℂ)⁻¹ * tr G` and
`R4C.stieltjesC (gramC Y) z = (p:ℂ)⁻¹ * tr Gc`; this is the finite identity multiplied out.
Careful at `p = 0`: `trace_cancel` handles it without a `0 < p` hypothesis. -/
theorem stieltjesC_gram_sub_gramC {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) {z : ℂ}
    (hz : 0 < z.im) (hd : 0 < d) :
    (d : ℂ) * R4C.stieltjesC (gram Y) z - (p : ℂ) * R4C.stieltjesC (gramC Y) z
      = ((p : ℂ) - d) / z := by
  have hfin := trace_resolvC_gram_sub_gramC Y hz hd
  simp only [R4C.stieltjesC]
  rw [trace_cancel, trace_cancel]
  exact hfin

end GenRMT
end StackedSVD
