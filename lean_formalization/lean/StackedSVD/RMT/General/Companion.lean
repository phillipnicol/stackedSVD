/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.ResolvDeriv

/-!
# Companion identities and the rank-one downdate, general noise law (Stage 1, unit G2)

`notes/archive/prop_single_table_general.md` section 5, unit G2. Two deterministic identities that
replace the Gaussian block-rotation route of item R0 (choice 5 of the note): the resolvent
form of `E G Eᵀ = 1 + z Gc` (the companion identity between the Wishart resolvent and the
companion resolvent) and Sherman-Morrison for the rank-one downdate `W₀ = W - g gᵀ`. Every
later item that needs a form of the noise vector `g = Eᵀ u` rewrites it through these two
identities into forms of the companion resolvent at `u` and of the full-noise resolvent at
`W = Eᵀ E`.

This file does not import `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` or any
`Vendor/COLT83/` file (choice 8 of the note). Section 3.3 of the unit brief supplies a
`private` copy of the Sherman-Morrison formula
(`Vendor/COLT83/Mathlib/Matrix/Loewner.lean:57`, Rémy Degenne, Apache 2.0) for that reason;
follow-up item F34 tracks its dedup against the original.

## Content

1. `R4C.qformC_ne_one`: the quadratic form of the resolvent never reaches `1` above the real
   axis, because its imaginary part is strictly positive off the zero vector.
2. `R4C.resolvC_sub_vecMulVec`, `R4C.cformC_sub_vecMulVec`, `R4C.cform2C_sub_vecMulVec`:
   Sherman-Morrison for the downdate `W - g gᵀ`, on the resolvent and on its two bilinear
   forms.
3. `GenRMT.cmat'_mul_resolvC_gram`, `GenRMT.smul_cmat'_mul_resolvC_mul_transpose`,
   `GenRMT.smul_cmat'_mul_resolvC_sq_mul_transpose`: the intertwining relation
   `Y G = Gc Y` and the companion identity `d⁻¹ Y G Yᵀ = 1 + z Gc`, at `G` and at `G²`.
4. `SpikedModel.W0_eq_gram_sub`, `SpikedModel.qformC_gvec_eq`: the two model corollaries,
   `W₀ = W - g gᵀ` at `W = gram Z` and the companion identity read at `g = Eᵀ u`.
-/

open MeasureTheory Filter Topology
open scoped Matrix

namespace StackedSVD

namespace R4C

/-- Rectangular twin of `cmat_mul` (`RMT/R4C.lean:171`): the cast to `ℂ` is multiplicative
across a rectangular-times-rectangular product. Public since F37 (2026-09-09), the canonical
copy for `RMT/General/` (the private twin of `MPtilde.lean` is dropped). -/
theorem cmat'_mul {p q r : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (B : Matrix (Fin q) (Fin r) ℝ) : cmat' (A * B) = cmat' A * cmat' B := by
  ext i j
  simp only [cmat'_apply, Matrix.mul_apply]
  push_cast
  rfl

/-- A real scalar passes through the square cast `cmat`. Public since F37 (2026-09-09), the
canonical copy for `RMT/General/` (the private twin of `FormsBridge.lean` is dropped). -/
theorem cmat_smul {d : ℕ} (r : ℝ) (A : Matrix (Fin d) (Fin d) ℝ) :
    cmat (r • A) = (r : ℂ) • cmat A := by
  ext i j
  simp only [cmat, Matrix.map_apply, Matrix.smul_apply, smul_eq_mul]
  push_cast
  ring

/-- The cast of the rank-one downdate, in the shape of `cmat_sub_smul` (`RMT/R4C.lean:186`). -/
private theorem cmat_sub_vecMulVec {d : ℕ} (W : Matrix (Fin d) (Fin d) ℝ) (g : Fin d → ℝ) :
    cmat (W - Matrix.vecMulVec g g) = cmat W - Matrix.vecMulVec (cvec g) (cvec g) := by
  ext i j
  simp only [cmat, cvec, Matrix.map_apply, Matrix.sub_apply, Matrix.vecMulVec_apply]
  push_cast
  ring

/-- `vecMulVec w w *ᵥ y = (w ⬝ᵥ y) • w`, the rank-one matrix applied to a vector. Public since
F37 (2026-09-09), the canonical copy for `RMT/General/` (the private twins of `Deloc.lean` and
`Iso.lean` are dropped). -/
theorem vecMulVec_self_mulVec {d : ℕ} (w y : Fin d → ℂ) :
    Matrix.vecMulVec w w *ᵥ y = (w ⬝ᵥ y) • w := by
  funext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul]
  rw [Finset.sum_mul]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- The bilinear form of a rank-one matrix factors. Public since F37 (2026-09-09), the
canonical copy for `RMT/General/` (the private twin of `Iso.lean` is dropped). -/
theorem dotProduct_vecMulVec_mulVec {d : ℕ} (x w y : Fin d → ℂ) :
    x ⬝ᵥ (Matrix.vecMulVec w w *ᵥ y) = (x ⬝ᵥ w) * (w ⬝ᵥ y) := by
  rw [vecMulVec_self_mulVec, dotProduct_smul, smul_eq_mul]
  ring

/-- Moving a matrix across a bilinear form onto the transpose, in the shape of
`R4.dotProduct_mulVec_comm` (`RMT/R4C.lean:64`) but for a possibly non-symmetric complex
matrix. -/
private theorem dotProduct_mulVec_eq {d : ℕ} (A : Matrix (Fin d) (Fin d) ℂ) (x z : Fin d → ℂ) :
    x ⬝ᵥ (A *ᵥ z) = (Aᵀ *ᵥ x) ⬝ᵥ z := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]

/-- The resolvent is symmetric, so it moves freely across a bilinear form. -/
private theorem resolvC_dotProduct_mulVec_eq {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ}
    (hW : W.IsHermitian) {z : ℂ} (hz : z.im ≠ 0) (x y : Fin d → ℂ) :
    x ⬝ᵥ (resolvC W z *ᵥ y) = (resolvC W z *ᵥ x) ⬝ᵥ y := by
  rw [dotProduct_mulVec_eq, ResolvDeriv.transpose_resolvC hW hz]

/-- **Sherman-Morrison formula**, copied from
`Vendor/COLT83/Mathlib/Matrix/Loewner.lean:57` (Rémy Degenne, Apache 2.0 license) because
`RMT/General/` may not import `Vendor/COLT83/` (choice 8 of
`notes/archive/prop_single_table_general.md`). Follow-up item F34 tracks the dedup of this copy
against the original. -/
private theorem inv_add_vecMulVec {d : ℕ} {A : Matrix (Fin d) (Fin d) ℂ} (hA : IsUnit A.det)
    (u v : Fin d → ℂ) (h : 1 + v ⬝ᵥ A⁻¹ *ᵥ u ≠ 0) :
    (A + Matrix.vecMulVec u v)⁻¹ =
      A⁻¹ - (1 + v ⬝ᵥ A⁻¹ *ᵥ u)⁻¹ • Matrix.vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹) := by
  set s := v ⬝ᵥ A⁻¹ *ᵥ u with hs
  refine Matrix.inv_eq_right_inv ?_
  have h1 : A * Matrix.vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹) = Matrix.vecMulVec u (v ᵥ* A⁻¹) := by
    rw [Matrix.mul_vecMulVec, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv A hA,
      Matrix.one_mulVec]
  have h2 : Matrix.vecMulVec u v * Matrix.vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹)
      = s • Matrix.vecMulVec u (v ᵥ* A⁻¹) := by
    rw [Matrix.vecMulVec_mul_vecMulVec, Matrix.vecMulVec_smul]
  have h3 : Matrix.vecMulVec u v * A⁻¹ = Matrix.vecMulVec u (v ᵥ* A⁻¹) :=
    Matrix.vecMulVec_mul _ _ _
  rw [Matrix.add_mul, Matrix.mul_sub, Matrix.mul_sub, Matrix.mul_nonsing_inv A hA,
    Matrix.mul_smul, Matrix.mul_smul, h1, h2, h3, smul_smul]
  have h4 : (1 + s)⁻¹ * s = 1 - (1 + s)⁻¹ := by field_simp; ring
  rw [h4, sub_smul, one_smul]
  abel

/-- The quadratic form of the resolvent never reaches 1 above the real axis: it is 0 at
`g = 0`, and its imaginary part is positive otherwise. This discharges the denominator of
the Sherman-Morrison downdate once and for all. -/
theorem qformC_ne_one {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian) {z : ℂ}
    (hz : 0 < z.im) (g : Fin d → ℝ) :
    qformC W z g ≠ 1 := by
  change cformC W z g g ≠ 1
  rcases eq_or_ne g 0 with hg | hg
  · subst hg
    have hcvec : cvec (0 : Fin d → ℝ) = 0 := by funext a; simp [cvec]
    simp [cformC, hcvec]
  · have hz' : z.im ≠ 0 := hz.ne'
    have hsum := cformC_eq_sum hW hz' g g
    have him : (cformC W z g g).im
        = ∑ a, (z.im / Complex.normSq ((hW.eigenvalues a : ℂ) - z))
            * (((R4.eigU hW)ᵀ *ᵥ g) a) ^ 2 := by
      rw [hsum, Complex.im_sum]
      refine Finset.sum_congr rfl fun a _ => ?_
      rw [Complex.mul_im]
      have hre0 : ((((R4.eigU hW)ᵀ *ᵥ g) a * ((R4.eigU hW)ᵀ *ᵥ g) a : ℝ) : ℂ).im = 0 := by
        simp
      have hre1 : ((((R4.eigU hW)ᵀ *ᵥ g) a * ((R4.eigU hW)ᵀ *ᵥ g) a : ℝ) : ℂ).re
          = ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2 := by rw [Complex.ofReal_re]; ring
      rw [hre0, hre1, mul_zero, zero_add, Complex.inv_im]
      have hIm : ((hW.eigenvalues a : ℂ) - z).im = -z.im := by simp
      rw [hIm]
      ring
    have hUg : (R4.eigU hW)ᵀ *ᵥ g ≠ 0 := fun h => hg (R4.eq_zero_of_transpose_eigU_mulVec hW h)
    obtain ⟨a0, ha0⟩ := Function.ne_iff.mp hUg
    have hnesq0 : ((hW.eigenvalues a0 : ℂ) - z) ≠ 0 := eigenvalue_sub_ne_zero hW hz' a0
    have hnormsq0 : 0 < Complex.normSq ((hW.eigenvalues a0 : ℂ) - z) :=
      Complex.normSq_pos.mpr hnesq0
    have hpos0 : 0 < (z.im / Complex.normSq ((hW.eigenvalues a0 : ℂ) - z))
        * (((R4.eigU hW)ᵀ *ᵥ g) a0) ^ 2 := by
      have h1 : 0 < z.im / Complex.normSq ((hW.eigenvalues a0 : ℂ) - z) := div_pos hz hnormsq0
      have h2 : 0 < (((R4.eigU hW)ᵀ *ᵥ g) a0) ^ 2 := sq_pos_of_ne_zero ha0
      exact mul_pos h1 h2
    have hnonneg : ∀ a, 0 ≤ (z.im / Complex.normSq ((hW.eigenvalues a : ℂ) - z))
        * (((R4.eigU hW)ᵀ *ᵥ g) a) ^ 2 := by
      intro a
      have hnesq : ((hW.eigenvalues a : ℂ) - z) ≠ 0 := eigenvalue_sub_ne_zero hW hz' a
      have hnormsq : 0 < Complex.normSq ((hW.eigenvalues a : ℂ) - z) :=
        Complex.normSq_pos.mpr hnesq
      positivity
    have hImPos : 0 < (cformC W z g g).im := by
      rw [him]
      exact Finset.sum_pos' (fun a _ => hnonneg a) ⟨a0, Finset.mem_univ a0, hpos0⟩
    intro hcon
    rw [hcon] at hImPos
    simp at hImPos

/-- **Sherman-Morrison for the rank-one downdate** `W - g gᵀ`. -/
theorem resolvC_sub_vecMulVec {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z : ℂ} (hz : 0 < z.im) (g : Fin d → ℝ) :
    resolvC (W - Matrix.vecMulVec g g) z
      = resolvC W z
        + (1 - qformC W z g)⁻¹ •
            Matrix.vecMulVec (resolvC W z *ᵥ cvec g)
              (resolvC W z *ᵥ cvec g) := by
  have hz' : z.im ≠ 0 := hz.ne'
  set A : Matrix (Fin d) (Fin d) ℂ := cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ) with hA
  have hAdet : IsUnit A.det := ResolvDeriv.isUnit_det_cmat_sub hW hz'
  have hAinv : A⁻¹ = resolvC W z := rfl
  have hvA : cvec g ᵥ* A⁻¹ = resolvC W z *ᵥ cvec g := by
    rw [hAinv, ← Matrix.mulVec_transpose, ResolvDeriv.transpose_resolvC hW hz']
  have hden0 : cvec g ⬝ᵥ A⁻¹ *ᵥ (-cvec g) = -qformC W z g := by
    rw [Matrix.mulVec_neg, dotProduct_neg, hAinv]
    rfl
  have hden : (1 : ℂ) + cvec g ⬝ᵥ A⁻¹ *ᵥ (-cvec g) ≠ 0 := by
    rw [hden0, ← sub_eq_add_neg]
    exact sub_ne_zero.mpr (Ne.symm (qformC_ne_one hW hz g))
  have hlhs : A + Matrix.vecMulVec (-cvec g) (cvec g)
      = cmat (W - Matrix.vecMulVec g g) - z • (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [hA, Matrix.neg_vecMulVec, cmat_sub_vecMulVec]
    abel
  have hkey := inv_add_vecMulVec hAdet (-cvec g) (cvec g) hden
  rw [hlhs] at hkey
  rw [hden0, ← sub_eq_add_neg, hvA, hAinv, Matrix.mulVec_neg, Matrix.neg_vecMulVec, smul_neg,
    sub_neg_eq_add] at hkey
  exact hkey

/-- The downdate on the bilinear form. -/
theorem cformC_sub_vecMulVec {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z : ℂ} (hz : 0 < z.im) (g x y : Fin d → ℝ) :
    cformC (W - Matrix.vecMulVec g g) z x y
      = cformC W z x y
        + cformC W z x g * cformC W z g y / (1 - qformC W z g) := by
  have hz' : z.im ≠ 0 := hz.ne'
  change cvec x ⬝ᵥ (resolvC (W - Matrix.vecMulVec g g) z *ᵥ cvec y) = _
  rw [resolvC_sub_vecMulVec hW hz g, Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec,
    dotProduct_smul, smul_eq_mul, dotProduct_vecMulVec_mulVec]
  have hxy : cvec x ⬝ᵥ (resolvC W z *ᵥ cvec y) = cformC W z x y := rfl
  have hxg : cvec x ⬝ᵥ (resolvC W z *ᵥ cvec g) = cformC W z x g := rfl
  have hgy : (resolvC W z *ᵥ cvec g) ⬝ᵥ cvec y = cformC W z g y :=
    (resolvC_dotProduct_mulVec_eq hW hz' (cvec g) (cvec y)).symm
  rw [hxy, hxg, hgy]
  ring

/-- The downdate on the squared bilinear form. -/
theorem cform2C_sub_vecMulVec {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z : ℂ} (hz : 0 < z.im) (g x y : Fin d → ℝ) :
    cform2C (W - Matrix.vecMulVec g g) z x y
      = cform2C W z x y
        + (cform2C W z x g * cformC W z g y
            + cformC W z x g * cform2C W z g y) / (1 - qformC W z g)
        + cformC W z x g * qform2C W z g * cformC W z g y
            / (1 - qformC W z g) ^ 2 := by
  have hz' : z.im ≠ 0 := hz.ne'
  change cvec x ⬝ᵥ ((resolvC (W - Matrix.vecMulVec g g) z
      * resolvC (W - Matrix.vecMulVec g g) z) *ᵥ cvec y) = _
  rw [resolvC_sub_vecMulVec hW hz g]
  set w : Fin d → ℂ := resolvC W z *ᵥ cvec g with hw
  set α : ℂ := (1 - qformC W z g)⁻¹ with hα
  have e1 : cvec x ⬝ᵥ (resolvC W z *ᵥ w) = cform2C W z x g := by
    rw [hw, Matrix.mulVec_mulVec]; rfl
  have e2 : w ⬝ᵥ cvec y = cformC W z g y := by
    rw [hw]
    exact (resolvC_dotProduct_mulVec_eq hW hz' (cvec g) (cvec y)).symm
  have e3 : cvec x ⬝ᵥ w = cformC W z x g := by
    rw [hw]
    rfl
  have e4 : w ⬝ᵥ (resolvC W z *ᵥ cvec y) = cform2C W z g y := by
    rw [hw, (resolvC_dotProduct_mulVec_eq hW hz' (cvec g) (resolvC W z *ᵥ cvec y)).symm,
      Matrix.mulVec_mulVec]; rfl
  have e5 : w ⬝ᵥ w = qform2C W z g := by
    rw [hw, (resolvC_dotProduct_mulVec_eq hW hz' (cvec g) (resolvC W z *ᵥ cvec g)).symm,
      Matrix.mulVec_mulVec]; rfl
  have hexpand : (resolvC W z + α • Matrix.vecMulVec w w)
      * (resolvC W z + α • Matrix.vecMulVec w w)
      = (resolvC W z * resolvC W z) + α • (resolvC W z * Matrix.vecMulVec w w)
        + α • (Matrix.vecMulVec w w * resolvC W z)
        + (α * α * (w ⬝ᵥ w)) • Matrix.vecMulVec w w := by
    rw [Matrix.add_mul, Matrix.mul_add, Matrix.mul_add, Matrix.mul_smul, Matrix.smul_mul,
      Matrix.smul_mul, Matrix.mul_smul, Matrix.vecMulVec_mul_vecMulVec, Matrix.vecMulVec_smul,
      smul_smul, smul_smul]
    abel
  rw [hexpand]
  simp only [Matrix.add_mulVec, dotProduct_add]
  have t0 : cvec x ⬝ᵥ ((resolvC W z * resolvC W z) *ᵥ cvec y) = cform2C W z x y := rfl
  have t1 : cvec x ⬝ᵥ ((α • (resolvC W z * Matrix.vecMulVec w w)) *ᵥ cvec y)
      = cform2C W z x g * cformC W z g y * α := by
    rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, ← Matrix.mulVec_mulVec,
      vecMulVec_self_mulVec, Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul, e2, e1]
    ring
  have t2 : cvec x ⬝ᵥ ((α • (Matrix.vecMulVec w w * resolvC W z)) *ᵥ cvec y)
      = cformC W z x g * cform2C W z g y * α := by
    rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, ← Matrix.mulVec_mulVec,
      vecMulVec_self_mulVec, dotProduct_smul, smul_eq_mul, e4, e3]
    ring
  have t3 : cvec x ⬝ᵥ (((α * α * (w ⬝ᵥ w)) • Matrix.vecMulVec w w) *ᵥ cvec y)
      = cformC W z x g * qform2C W z g * cformC W z g y * (α * α) := by
    rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, dotProduct_vecMulVec_mulVec, e3, e2, e5]
    ring
  rw [t0, t1, t2, t3, hα]
  field_simp
  ring

end R4C

namespace GenRMT


/-- A generic "intertwining implies conjugation of inverses" fact: if `X a = b X` with `a`,
`b` right/left invertible by `G`, `Gc` respectively, then `X G = Gc X`. Isolates the matrix
associativity bookkeeping of the intertwining relation from the `gram`/`gramC` instance. -/
private theorem intertwine_of_conj {n p : ℕ} {X : Matrix (Fin p) (Fin n) ℂ}
    {a G : Matrix (Fin n) (Fin n) ℂ} {b Gc : Matrix (Fin p) (Fin p) ℂ}
    (haG : a * G = 1) (hGcb : Gc * b = 1) (hshift : X * a = b * X) :
    X * G = Gc * X := by
  have e1 : Gc * X = Gc * (X * a * G) := by
    rw [Matrix.mul_assoc X a G, haG, Matrix.mul_one]
  rw [e1, hshift, ← Matrix.mul_assoc Gc (b * X) G, ← Matrix.mul_assoc Gc b X, hGcb,
    Matrix.one_mul]

/-- **The intertwining relation** `Y G = Gc Y`, from `Y (d⁻¹ YᵀY) = (d⁻¹ Y Yᵀ) Y`. -/
theorem cmat'_mul_resolvC_gram {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) {z : ℂ}
    (hz : 0 < z.im) :
    R4C.cmat' Y * R4C.resolvC (gram Y) z = R4C.resolvC (gramC Y) z * R4C.cmat' Y := by
  have hz' : z.im ≠ 0 := hz.ne'
  have hAherm := gram_isHermitian Y
  have hCherm := gramC_isHermitian Y
  have hreal : Y * gram Y = gramC Y * Y := by
    change Y * (((d : ℝ))⁻¹ • (Yᵀ * Y)) = (((d : ℝ))⁻¹ • (Y * Yᵀ)) * Y
    rw [Matrix.mul_smul, Matrix.smul_mul]
    congr 1
    rw [Matrix.mul_assoc]
  have h1 : R4C.cmat' (Y * gram Y) = R4C.cmat' Y * R4C.cmat (gram Y) := by
    rw [R4C.cmat'_mul, R4C.cmat'_eq_cmat]
  have h2 : R4C.cmat' (gramC Y * Y) = R4C.cmat (gramC Y) * R4C.cmat' Y := by
    rw [R4C.cmat'_mul, R4C.cmat'_eq_cmat]
  have hshift : R4C.cmat' Y * (R4C.cmat (gram Y) - z • (1 : Matrix (Fin d) (Fin d) ℂ))
      = (R4C.cmat (gramC Y) - z • (1 : Matrix (Fin p) (Fin p) ℂ)) * R4C.cmat' Y := by
    rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one,
      Matrix.one_mul, ← h1, ← h2, hreal]
  exact intertwine_of_conj (ResolvDeriv.cmat_sub_mul_resolvC hAherm hz')
    (ResolvDeriv.resolvC_mul_cmat_sub hCherm hz') hshift

/-- **The companion identity** `d⁻¹ Y G Yᵀ = 1 + z Gc`. -/
theorem smul_cmat'_mul_resolvC_mul_transpose {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ)
    {z : ℂ} (hz : 0 < z.im) (hd : 0 < d) :
    ((d : ℂ))⁻¹ • (R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ)
      = 1 + z • R4C.resolvC (gramC Y) z := by
  have hz' : z.im ≠ 0 := hz.ne'
  have hCherm := gramC_isHermitian Y
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hdR : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hAAt : R4C.cmat' Y * (R4C.cmat' Y)ᵀ = (d : ℂ) • R4C.cmat (gramC Y) := by
    have hYYt : Y * Yᵀ = (d : ℝ) • gramC Y := by
      change Y * Yᵀ = (d : ℝ) • (((d : ℝ))⁻¹ • (Y * Yᵀ))
      rw [smul_smul, mul_inv_cancel₀ hdR, one_smul]
    have h1 : R4C.cmat' (Y * Yᵀ) = R4C.cmat' Y * (R4C.cmat' Y)ᵀ := by
      rw [R4C.cmat'_mul, R4C.cmat'_transpose]
    have h2 : R4C.cmat' ((d : ℝ) • gramC Y) = (d : ℂ) • R4C.cmat (gramC Y) := by
      rw [R4C.cmat'_eq_cmat, R4C.cmat_smul]
      congr 1
    rw [← h1, hYYt, h2]
  have hGcmat : R4C.resolvC (gramC Y) z * R4C.cmat (gramC Y)
      = 1 + z • R4C.resolvC (gramC Y) z := by
    have h := ResolvDeriv.resolvC_mul_cmat_sub hCherm hz'
    rwa [Matrix.mul_sub, Matrix.mul_smul, Matrix.mul_one, sub_eq_iff_eq_add] at h
  calc ((d : ℂ))⁻¹ • (R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ)
      = ((d : ℂ))⁻¹ • (R4C.resolvC (gramC Y) z * R4C.cmat' Y * (R4C.cmat' Y)ᵀ) := by
        rw [cmat'_mul_resolvC_gram Y hz]
    _ = ((d : ℂ))⁻¹ • (R4C.resolvC (gramC Y) z * (R4C.cmat' Y * (R4C.cmat' Y)ᵀ)) := by
        rw [Matrix.mul_assoc]
    _ = ((d : ℂ))⁻¹ • (R4C.resolvC (gramC Y) z * ((d : ℂ) • R4C.cmat (gramC Y))) := by
        rw [hAAt]
    _ = ((d : ℂ))⁻¹ • ((d : ℂ) • (R4C.resolvC (gramC Y) z * R4C.cmat (gramC Y))) := by
        rw [Matrix.mul_smul]
    _ = R4C.resolvC (gramC Y) z * R4C.cmat (gramC Y) := by
        rw [smul_smul, inv_mul_cancel₀ hdC, one_smul]
    _ = 1 + z • R4C.resolvC (gramC Y) z := hGcmat

/-- **The squared companion identity** `d⁻¹ Y G² Yᵀ = Gc + z Gc²`. -/
theorem smul_cmat'_mul_resolvC_sq_mul_transpose {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ)
    {z : ℂ} (hz : 0 < z.im) (hd : 0 < d) :
    ((d : ℂ))⁻¹ • (R4C.cmat' Y * (R4C.resolvC (gram Y) z * R4C.resolvC (gram Y) z)
        * (R4C.cmat' Y)ᵀ)
      = R4C.resolvC (gramC Y) z + z • (R4C.resolvC (gramC Y) z * R4C.resolvC (gramC Y) z) := by
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hAGAt : R4C.cmat' Y * R4C.resolvC (gram Y) z * (R4C.cmat' Y)ᵀ
      = (d : ℂ) • (1 + z • R4C.resolvC (gramC Y) z) := by
    have h := smul_cmat'_mul_resolvC_mul_transpose Y hz hd
    rw [← h, smul_smul, mul_inv_cancel₀ hdC, one_smul]
  have e1 : R4C.cmat' Y * (R4C.resolvC (gram Y) z * R4C.resolvC (gram Y) z)
      = R4C.resolvC (gramC Y) z * (R4C.cmat' Y * R4C.resolvC (gram Y) z) := by
    have step := cmat'_mul_resolvC_gram Y hz
    calc R4C.cmat' Y * (R4C.resolvC (gram Y) z * R4C.resolvC (gram Y) z)
        = R4C.cmat' Y * R4C.resolvC (gram Y) z * R4C.resolvC (gram Y) z := by
          rw [← Matrix.mul_assoc]
      _ = R4C.resolvC (gramC Y) z * R4C.cmat' Y * R4C.resolvC (gram Y) z := by rw [step]
      _ = R4C.resolvC (gramC Y) z * (R4C.cmat' Y * R4C.resolvC (gram Y) z) := by
          rw [Matrix.mul_assoc]
  have hstep : R4C.cmat' Y * (R4C.resolvC (gram Y) z * R4C.resolvC (gram Y) z) * (R4C.cmat' Y)ᵀ
      = (d : ℂ) • (R4C.resolvC (gramC Y) z
          + z • (R4C.resolvC (gramC Y) z * R4C.resolvC (gramC Y) z)) := by
    rw [e1, Matrix.mul_assoc, hAGAt, Matrix.mul_smul, Matrix.mul_add, Matrix.mul_one,
      Matrix.mul_smul]
  rw [hstep, smul_smul, inv_mul_cancel₀ hdC, one_smul]

end GenRMT

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **`W₀ = W - g gᵀ`**, the downdate the stage runs on. -/
theorem W0_eq_gram_sub (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    m.W0 N ω = GenRMT.gram (m.Z N ω) - Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω) := by
  have hdR : (d N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (m.hd N).ne'
  have hw : (m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N) = Real.sqrt (d N) • m.gvec N ω :=
    (m.sqrt_smul_gvec N ω).symm
  have hvv : Matrix.vecMulVec ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))
      ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))
      = (d N : ℝ) • Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω) := by
    rw [hw, Matrix.smul_vecMulVec, Matrix.vecMulVec_smul, smul_smul,
      Real.mul_self_sqrt (Nat.cast_nonneg (d N))]
  rw [m.W0_eq_smul N ω, hvv, smul_sub, smul_smul, inv_mul_cancel₀ hdR, one_smul]
  rfl

/-- `cvec` commutes with a real scalar. Public since F37 (2026-09-09), the canonical copy for
`RMT/General/` (the private twins of `DelocLimits.lean`, `FormsBridge.lean` and
`IsoMixed.lean` are dropped). -/
theorem cvec_smul {k : ℕ} (c : ℝ) (v : Fin k → ℝ) :
    R4C.cvec (c • v) = (c : ℂ) • R4C.cvec v := by
  funext a
  simp [R4C.cvec]

/-- Moves the transposed rectangular cast across `cvec`. Public since F37 (2026-09-09), the
canonical copy for `RMT/General/` (the private twins of `FormsBridge.lean` and
`IsoMixed.lean` are dropped). -/
theorem cmat'_transpose_mulVec_cvec {p k : ℕ} (Y : Matrix (Fin p) (Fin k) ℝ)
    (x : Fin p → ℝ) : (R4C.cmat' Y)ᵀ *ᵥ R4C.cvec x = R4C.cvec (Yᵀ *ᵥ x) := by
  funext a
  simp only [R4C.cvec, Matrix.mulVec, dotProduct, Matrix.transpose_apply, R4C.cmat'_apply]
  push_cast
  rfl

/-- `cvec` preserves the dot product. Public since F37 (2026-09-09), the canonical copy for
`RMT/General/` (the private twins of `FormsBridge.lean`, `Iso.lean` and `IsoMixed.lean` are
dropped). -/
theorem cvec_dotProduct_cvec {k : ℕ} (x y : Fin k → ℝ) :
    R4C.cvec x ⬝ᵥ R4C.cvec y = ((x ⬝ᵥ y : ℝ) : ℂ) := by
  have h := R4C.dotProduct_cmat_mulVec (1 : Matrix (Fin k) (Fin k) ℝ) x y
  rwa [R4C.cmat_one, Matrix.one_mulVec, Matrix.one_mulVec] at h

/-- Rectangular twin of `R4C.dotProduct_mulVec_eq`: moves a (possibly non-square) matrix
across a bilinear form onto its transpose. Public since F37 (2026-09-09), the canonical copy
for `RMT/General/` (the private twins of `FormsBridge.lean` and `IsoMixed.lean` are
dropped). -/
theorem dotProduct_mulVec_eq_rect {p k : ℕ} (M : Matrix (Fin p) (Fin k) ℂ)
    (x : Fin p → ℂ) (z : Fin k → ℂ) : x ⬝ᵥ (M *ᵥ z) = (Mᵀ *ᵥ x) ⬝ᵥ z := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]

/-- `gᵀ G g = 1 + z uᵀ Gc u`, the companion identity at `u`. Choice 5 of the plan note:
the denominator `1 - gᵀ G g` of the downdate is `-z uᵀ Gc u`. -/
theorem qformC_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im) :
    R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = 1 + z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) := by
  set Z := m.Z N ω with hZ
  set u := WithLp.ofLp (m.u N) with hu
  set G := R4C.resolvC (GenRMT.gram Z) z with hG
  set Gc := R4C.resolvC (GenRMT.gramC Z) z with hGc
  set A := R4C.cmat' Z with hA
  have hgw : R4C.cvec (m.gvec N ω) = (((Real.sqrt (d N))⁻¹ : ℝ) : ℂ) • R4C.cvec (Zᵀ *ᵥ u) := by
    rw [m.gvec_eq_smul N ω, cvec_smul]
  have hcvw : R4C.cvec (Zᵀ *ᵥ u) = Aᵀ *ᵥ R4C.cvec u := (cmat'_transpose_mulVec_cvec Z u).symm
  have hscal : (((Real.sqrt (d N))⁻¹ : ℝ) : ℂ) * (((Real.sqrt (d N))⁻¹ : ℝ) : ℂ)
      = (d N : ℂ)⁻¹ := by
    rw [← Complex.ofReal_mul]
    have step : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = (d N : ℝ)⁻¹ := by
      rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg (d N))]
    rw [step, Complex.ofReal_inv]
    congr 1
  have hbridge : R4C.cvec u ⬝ᵥ ((A * G * Aᵀ) *ᵥ R4C.cvec u)
      = (Aᵀ *ᵥ R4C.cvec u) ⬝ᵥ (G *ᵥ (Aᵀ *ᵥ R4C.cvec u)) := by
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
      dotProduct_mulVec_eq_rect A (R4C.cvec u) (G *ᵥ (Aᵀ *ᵥ R4C.cvec u))]
  have hq1 : R4C.qformC (GenRMT.gram Z) z (m.gvec N ω)
      = (d N : ℂ)⁻¹ * (R4C.cvec u ⬝ᵥ ((A * G * Aᵀ) *ᵥ R4C.cvec u)) := by
    change R4C.cvec (m.gvec N ω) ⬝ᵥ (G *ᵥ R4C.cvec (m.gvec N ω)) = _
    rw [hgw, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_smul, smul_eq_mul, hscal,
      hcvw, hbridge]
  rw [hq1]
  have hswap : (d N : ℂ)⁻¹ * (R4C.cvec u ⬝ᵥ ((A * G * Aᵀ) *ᵥ R4C.cvec u))
      = R4C.cvec u ⬝ᵥ ((1 + z • Gc) *ᵥ R4C.cvec u) := by
    have h := GenRMT.smul_cmat'_mul_resolvC_mul_transpose Z hz (m.hd N)
    have h2 : R4C.cvec u ⬝ᵥ (((d N : ℂ)⁻¹ • (A * G * Aᵀ)) *ᵥ R4C.cvec u)
        = R4C.cvec u ⬝ᵥ ((1 + z • Gc) *ᵥ R4C.cvec u) := by rw [h]
    rwa [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h2
  rw [hswap, Matrix.add_mulVec, dotProduct_add, Matrix.one_mulVec, Matrix.smul_mulVec,
    dotProduct_smul, smul_eq_mul]
  have hreal : u ⬝ᵥ u = 1 := m.dotProduct_u_self N
  have huu : R4C.cvec u ⬝ᵥ R4C.cvec u = 1 := by
    rw [cvec_dotProduct_cvec, hreal]
    norm_num
  rw [huu]
  rfl

end SpikedModel
end StackedSVD
