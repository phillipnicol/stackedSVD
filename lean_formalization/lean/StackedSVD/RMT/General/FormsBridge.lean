/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.General.MPtilde

/-!
# The deterministic bridge from `gram Z` to `W₀`, the scalar identities of `mTildeC`, and the
companion rescaling (Stage 1, unit G6a)

`notes/archive/prop_single_table_general.md` section 5, unit G6 (the deterministic half, unit G6a,
2026-09-09). Every statement here is a pointwise algebraic identity or inequality: it holds
for every `ω` and every `z` with `0 < z.im`, with no probability and no `TendstoInProb`. Unit
G6 (a later agent) rewrites each of the six fields of `ResolventFormsC`
(`RMT/General/Defs.lean`) through these identities and then passes to the limit.

## Content

1. Section A (namespace `StackedSVD.SpikedModel`): the downdate `W₀ = W - g gᵀ` read through
   the companion identity `a = 1 + z b` (`Companion.lean:429`, `a := qformC W z g`,
   `b := qformC Wc z u`), so the denominator `1 - a` of the Sherman-Morrison downdate becomes
   `-(z b)`.
2. Section B (namespace `StackedSVD.MP`): the derivative of `mTildeC`, its relation to
   `mCDeriv`, and `mC = -1 - (z m̃)⁻¹`, the scalar twin of item A2.
3. Section C (namespace `StackedSVD.GenRMT`): the companion resolvent as a rescaled `gram`
   resolvent at the transposed matrix, and the Lipschitz bound of `cformC`/`cform2C` in `z`.

Imports only `RMT/General/Companion` and `RMT/General/MPtilde` (choice 8 of the note): no
`RMT/General/Iso`, `RMT/General/Trace`, `RMT/General/R3minus`.
-/

open MeasureTheory Filter Topology
open scoped Matrix

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **A0.** `z b ≠ 0`, `b` the companion quadratic form at `u`. From `qformC_ne_one` on `W` at
`g` (`R4C.qformC_ne_one`) and the companion identity `a = 1 + z b`
(`SpikedModel.qformC_gvec_eq`, `Companion.lean:429`): if `z b = 0` then `a = 1`, contradicting
`qformC_ne_one`. -/
theorem z_mul_qformC_gramC_ne_zero (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ}
    (hz : 0 < z.im) :
    z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) ≠ 0 := by
  intro hcon
  have ha := m.qformC_gvec_eq N ω hz
  rw [hcon, add_zero] at ha
  exact R4C.qformC_ne_one (GenRMT.gram_isHermitian (m.Z N ω)) hz (m.gvec N ω) ha

/-- Algebra helper for the A-series rewrites: `1 - a = -(z b)` turns a downdate denominator
`p / (1 - a)` into `-(p / (z b))`, absorbing the sign flip once and for all. -/
private theorem div_one_sub_eq_neg_div (a z b p : ℂ) (ha : a = 1 + z * b) :
    p / (1 - a) = -(p / (z * b)) := by
  have hsub : (1 : ℂ) - a = -(z * b) := by rw [ha]; ring
  rw [hsub, div_neg]

/-- Squared twin of `div_one_sub_eq_neg_div`, for the `A5` to `A7` rewrites: the sign cancels
under the square, `(1 - a) ^ 2 = (z b) ^ 2`. -/
private theorem div_one_sub_sq_eq_div (a z b p : ℂ) (ha : a = 1 + z * b) :
    p / (1 - a) ^ 2 = p / (z * b) ^ 2 := by
  have hsub : (1 : ℂ) - a = -(z * b) := by rw [ha]; ring
  rw [hsub, neg_sq]

/-- **A1.** The downdate on `qformC` at a general `x`, with the denominator `1 - a` already
rewritten as `-(z b)` through the companion identity. From `SpikedModel.W0_eq_gram_sub`
(`Companion.lean:391`), `R4C.cformC_sub_vecMulVec` (`Companion.lean:206`) and
`div_one_sub_eq_neg_div`. -/
theorem qformC_W0_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im)
    (x : Fin (d N) → ℝ) :
    R4C.qformC (m.W0 N ω) z x
      = R4C.qformC (GenRMT.gram (m.Z N ω)) z x
        - R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
            * R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) x
            / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) := by
  have hW := GenRMT.gram_isHermitian (m.Z N ω)
  have ha := m.qformC_gvec_eq N ω hz
  have hq : R4C.qformC (m.W0 N ω) z x = R4C.cformC (m.W0 N ω) z x x := rfl
  have hWq : R4C.qformC (GenRMT.gram (m.Z N ω)) z x
      = R4C.cformC (GenRMT.gram (m.Z N ω)) z x x := rfl
  rw [hq, m.W0_eq_gram_sub N ω, R4C.cformC_sub_vecMulVec hW hz (m.gvec N ω) x x, hWq,
    div_one_sub_eq_neg_div _ _ _ _ ha]
  ring

/-- **A2.** The downdate on `qformC` read at `g` itself: `a + a² / (1 - a) = -1 - (z b)⁻¹`. -/
theorem qformC_W0_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im) :
    R4C.qformC (m.W0 N ω) z (m.gvec N ω)
      = -1 - (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))⁻¹ := by
  have hW := GenRMT.gram_isHermitian (m.Z N ω)
  have ha := m.qformC_gvec_eq N ω hz
  have hzb := z_mul_qformC_gramC_ne_zero m N ω hz
  have hzne : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hbne : R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) ≠ 0 :=
    right_ne_zero_of_mul hzb
  have hsub : (1 : ℂ) - R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = -(z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) := by
    rw [ha]; ring
  have hq : R4C.qformC (m.W0 N ω) z (m.gvec N ω)
      = R4C.cformC (m.W0 N ω) z (m.gvec N ω) (m.gvec N ω) := rfl
  have haa : R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  rw [hq, m.W0_eq_gram_sub N ω,
    R4C.cformC_sub_vecMulVec hW hz (m.gvec N ω) (m.gvec N ω) (m.gvec N ω), haa, hsub, ha]
  field_simp
  ring

/-- **A3.** The downdate on `cformC x g`. From `R4C.cformC_sub_vecMulVec` at `y = g`. -/
theorem cformC_W0_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im)
    (x : Fin (d N) → ℝ) :
    R4C.cformC (m.W0 N ω) z x (m.gvec N ω)
      = -(R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
          / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))) := by
  have hW := GenRMT.gram_isHermitian (m.Z N ω)
  have ha := m.qformC_gvec_eq N ω hz
  have hzb := z_mul_qformC_gramC_ne_zero m N ω hz
  have hzne : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hbne : R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) ≠ 0 :=
    right_ne_zero_of_mul hzb
  have hsub : (1 : ℂ) - R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = -(z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) := by
    rw [ha]; ring
  have hgg : R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  rw [m.W0_eq_gram_sub N ω,
    R4C.cformC_sub_vecMulVec hW hz (m.gvec N ω) x (m.gvec N ω), hgg, hsub, ha]
  field_simp
  ring

end SpikedModel

namespace MP

/-- **B1.** The derivative of `mTildeC c` in `z`. -/
noncomputable def mTildeCDeriv (c : ℝ) (z : ℂ) : ℂ := (mCDeriv c z - (1 - (c : ℂ)) / z ^ 2) / c

/-- **B2.** `mTildeC c` is differentiable at every `z` with `0 < z.im`, with derivative
`mTildeCDeriv c z`. `mTildeC c = fun w => (mC c w + (1 - c) / w) / c` (`MPtilde.lean:48`) is a
sum and a constant quotient of `mC c` (`MP7.lean:354`) and `fun w => (1 - c) / w`, so the
derivative is built directly, with no need for `HasDerivAt.congr_of_eventuallyEq`. -/
theorem hasDerivAt_mTildeC {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    HasDerivAt (mTildeC c) (mTildeCDeriv c z) z := by
  have hzne : z ≠ 0 := ne_zero_of_im_pos hz
  have h1 : HasDerivAt (mC c) (mCDeriv c z) z := hasDerivAt_mC hc.le hz
  have h2 : HasDerivAt (fun w : ℂ => (1 - (c : ℂ)) / w)
      ((0 * z - (1 - (c : ℂ)) * 1) / z ^ 2) z :=
    (hasDerivAt_const z (1 - (c : ℂ))).div (hasDerivAt_id' (x := z)) hzne
  have h3 : HasDerivAt (fun w : ℂ => mC c w + (1 - (c : ℂ)) / w)
      (mCDeriv c z + (0 * z - (1 - (c : ℂ)) * 1) / z ^ 2) z := h1.add h2
  have h4 : HasDerivAt (fun w : ℂ => (mC c w + (1 - (c : ℂ)) / w) / (c : ℂ))
      ((mCDeriv c z + (0 * z - (1 - (c : ℂ)) * 1) / z ^ 2) / (c : ℂ)) z := h3.div_const (c : ℂ)
  have hfun : (fun w : ℂ => (mC c w + (1 - (c : ℂ)) / w) / (c : ℂ)) = mTildeC c := rfl
  rw [hfun] at h4
  convert h4 using 1
  unfold mTildeCDeriv
  ring

/-- **B3.** `z m̃ ≠ 0`, both factors having positive imaginary part
(`MPtilde.lean:134`, `im_mTildeC_pos`). -/
theorem z_mul_mTildeC_ne_zero {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    z * mTildeC c z ≠ 0 :=
  mul_ne_zero (ne_zero_of_im_pos hz) (ne_zero_of_im_pos (im_mTildeC_pos hc hz))

/-- **B4.** `mC = -1 - (z m̃)⁻¹`, the scalar twin of item A2. From `quad_mTildeC`
(`MPtilde.lean:69`) and `c_mul_mTildeC` (`MPtilde.lean:61`). -/
theorem mC_eq_neg_one_sub_inv {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    mC c z = -1 - (z * mTildeC c z)⁻¹ := by
  have hzne : z ≠ 0 := ne_zero_of_im_pos hz
  have hmne : mTildeC c z ≠ 0 := ne_zero_of_im_pos (im_mTildeC_pos hc hz)
  have hmc' : (c : ℂ) * mTildeC c z - (1 - (c : ℂ)) / z = mC c z := by
    linear_combination c_mul_mTildeC (c := c) (z := z) hc
  rw [← hmc']
  have hQ := quad_mTildeC hc hz
  field_simp
  linear_combination hQ

/-- **B5.** The derivative of `mC` read through `mTildeC`, by uniqueness of the derivative:
`fun w => -1 - (w m̃(w))⁻¹` has derivative `(m̃ + z m̃') / (z m̃)²` at `z` (built from B2, B3
and the standard `HasDerivAt` combinators) and agrees with `mC c` on the open upper half
plane (B4), so its derivative at `z` is also `mCDeriv c z` (`hasDerivAt_mC`,
`MP7.lean:354`). -/
theorem mCDeriv_eq_of_mTildeC {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    mCDeriv c z = (mTildeC c z + z * mTildeCDeriv c z) / (z * mTildeC c z) ^ 2 := by
  have hzm := z_mul_mTildeC_ne_zero hc hz
  have hid : HasDerivAt (fun w : ℂ => w) 1 z := hasDerivAt_id' (x := z)
  have hmt : HasDerivAt (mTildeC c) (mTildeCDeriv c z) z := hasDerivAt_mTildeC hc hz
  have hprod : HasDerivAt (fun w : ℂ => w * mTildeC c w)
      (1 * mTildeC c z + z * mTildeCDeriv c z) z := hid.mul hmt
  have hinv : HasDerivAt (fun w : ℂ => (w * mTildeC c w)⁻¹)
      (-(1 * mTildeC c z + z * mTildeCDeriv c z) / (z * mTildeC c z) ^ 2) z :=
    hprod.inv hzm
  have hfull : HasDerivAt (fun w : ℂ => -1 - (w * mTildeC c w)⁻¹)
      (-(-(1 * mTildeC c z + z * mTildeCDeriv c z) / (z * mTildeC c z) ^ 2)) z :=
    hinv.const_sub (-1 : ℂ)
  have heq : mC c =ᶠ[𝓝 z] fun w : ℂ => -1 - (w * mTildeC c w)⁻¹ := by
    have hopen : IsOpen {w : ℂ | 0 < w.im} := isOpen_lt continuous_const Complex.continuous_im
    filter_upwards [hopen.mem_nhds hz] with w hw
    exact mC_eq_neg_one_sub_inv hc hw
  have hresult := hfull.congr_of_eventuallyEq heq
  have hu := hresult.unique (hasDerivAt_mC hc.le hz)
  rw [← hu]
  ring

/-- **B6.** `m̃'` read through `mCDeriv` at the dual ratio `c⁻¹` and the rescaled point `z/c`,
by the same uniqueness method as B5 applied to `c_mul_mTildeC_eq_mC_inv`
(`MPtilde.lean:121`). -/
theorem mTildeCDeriv_eq_mCDeriv_inv {c : ℝ} {z : ℂ} (hc : 0 < c) (hz : 0 < z.im) :
    mTildeCDeriv c z = ((c : ℂ)⁻¹) ^ 2 * mCDeriv c⁻¹ (z / c) := by
  have hcne : (c : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hc.ne'
  have hzcim : 0 < (z / (c : ℂ)).im := by
    rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have hmt : HasDerivAt (mTildeC c) (mTildeCDeriv c z) z := hasDerivAt_mTildeC hc hz
  have hL : HasDerivAt (fun w : ℂ => (c : ℂ) * mTildeC c w) ((c : ℂ) * mTildeCDeriv c z) z :=
    hmt.const_mul (c : ℂ)
  have hdivc : HasDerivAt (fun w : ℂ => w / (c : ℂ)) (1 / (c : ℂ)) z :=
    (hasDerivAt_id' (x := z)).div_const (c : ℂ)
  have hmcinv : HasDerivAt (mC c⁻¹) (mCDeriv c⁻¹ (z / (c : ℂ))) (z / (c : ℂ)) :=
    hasDerivAt_mC (inv_nonneg.mpr hc.le) hzcim
  have hR : HasDerivAt (mC c⁻¹ ∘ fun w : ℂ => w / (c : ℂ))
      (mCDeriv c⁻¹ (z / (c : ℂ)) * (1 / (c : ℂ))) z :=
    HasDerivAt.comp_of_eq z hmcinv hdivc rfl
  have heq : (fun w : ℂ => (c : ℂ) * mTildeC c w) =ᶠ[𝓝 z] mC c⁻¹ ∘ fun w : ℂ => w / (c : ℂ) := by
    have hopen : IsOpen {w : ℂ | 0 < w.im} := isOpen_lt continuous_const Complex.continuous_im
    filter_upwards [hopen.mem_nhds hz] with w hw
    exact c_mul_mTildeC_eq_mC_inv hc hw
  have hresult := hL.congr_of_eventuallyEq heq.symm
  have hu := hresult.unique hR
  -- Abstract `mCDeriv c⁻¹ (z / c)` to a local constant first, so `field_simp` never rewrites
  -- its argument `c⁻¹` to `1 / c` (it would otherwise, splitting one atom into two).
  set K := mCDeriv c⁻¹ (z / (c : ℂ)) with hK
  field_simp at hu ⊢
  linear_combination hu

end MP

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! F37 (2026-09-09): `cvec_smul`, `cmat'_transpose_mulVec_cvec`, `cvec_dotProduct_cvec` and
`dotProduct_mulVec_eq_rect` used to be repeated here (`private`, since `Companion.lean`'s
copies were `private` too). `Companion.lean:404` to `:425` now carries the public canonical
copy of each; this file uses those directly. -/

/-- **A4.** The derivative twin of the companion identity at `u`: `Φ²_g(z) = b + z b₂`, with
`b = Φ_u(z)` and `b₂ = Φ²_u(z)` on the companion side. Template:
`SpikedModel.qformC_gvec_eq` (`Companion.lean:429` to `:472`), with
`GenRMT.smul_cmat'_mul_resolvC_sq_mul_transpose` (`Companion.lean:358`) in place of
`smul_cmat'_mul_resolvC_mul_transpose`. -/
theorem qform2C_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im) :
    R4C.qform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
        + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) := by
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
  have hbridge : R4C.cvec u ⬝ᵥ ((A * (G * G) * Aᵀ) *ᵥ R4C.cvec u)
      = (Aᵀ *ᵥ R4C.cvec u) ⬝ᵥ ((G * G) *ᵥ (Aᵀ *ᵥ R4C.cvec u)) := by
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
      dotProduct_mulVec_eq_rect A (R4C.cvec u) ((G * G) *ᵥ (Aᵀ *ᵥ R4C.cvec u))]
  have hq1 : R4C.qform2C (GenRMT.gram Z) z (m.gvec N ω)
      = (d N : ℂ)⁻¹ * (R4C.cvec u ⬝ᵥ ((A * (G * G) * Aᵀ) *ᵥ R4C.cvec u)) := by
    change R4C.cvec (m.gvec N ω) ⬝ᵥ ((G * G) *ᵥ R4C.cvec (m.gvec N ω)) = _
    rw [hgw, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_smul, smul_eq_mul, hscal,
      hcvw, hbridge]
  rw [hq1]
  have hswap : (d N : ℂ)⁻¹ * (R4C.cvec u ⬝ᵥ ((A * (G * G) * Aᵀ) *ᵥ R4C.cvec u))
      = R4C.cvec u ⬝ᵥ ((Gc + z • (Gc * Gc)) *ᵥ R4C.cvec u) := by
    have h := GenRMT.smul_cmat'_mul_resolvC_sq_mul_transpose Z hz (m.hd N)
    have h2 : R4C.cvec u ⬝ᵥ (((d N : ℂ)⁻¹ • (A * (G * G) * Aᵀ)) *ᵥ R4C.cvec u)
        = R4C.cvec u ⬝ᵥ ((Gc + z • (Gc * Gc)) *ᵥ R4C.cvec u) := by rw [h]
    rwa [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h2
  rw [hswap, Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]
  rfl

/-- **A5.** The downdate on `cform2C` at general `x`, `y`, with `qform2C W z g` already
rewritten to `b + z b₂` by A4 and the denominators rewritten to powers of `z b` by the
companion identity. From `R4C.cform2C_sub_vecMulVec` (`Companion.lean:223`). -/
theorem cform2C_W0_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im)
    (x y : Fin (d N) → ℝ) :
    R4C.cform2C (m.W0 N ω) z x y
      = R4C.cform2C (GenRMT.gram (m.Z N ω)) z x y
        - (R4C.cform2C (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
              * R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) y
            + R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
              * R4C.cform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) y)
            / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
        + R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
            * (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
                + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
            * R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) y
            / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) ^ 2 := by
  have hW := GenRMT.gram_isHermitian (m.Z N ω)
  have ha := m.qformC_gvec_eq N ω hz
  have hA4 := m.qform2C_gvec_eq N ω hz
  have hsub : (1 : ℂ) - R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = -(z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) := by
    rw [ha]; ring
  rw [m.W0_eq_gram_sub N ω, R4C.cform2C_sub_vecMulVec hW hz (m.gvec N ω) x y, hA4, hsub,
    div_neg, neg_sq]
  ring

/-- **A6.** The downdate on `qform2C` read at `g` itself: the three terms of A5 collapse
because `(1 - a)² + 2 a (1 - a) + a² = ((1 - a) + a)² = 1`. -/
theorem qform2C_W0_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im) :
    R4C.qform2C (m.W0 N ω) z (m.gvec N ω)
      = (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
          + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
          / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) ^ 2 := by
  have hW := GenRMT.gram_isHermitian (m.Z N ω)
  have ha := m.qformC_gvec_eq N ω hz
  have hA4 := m.qform2C_gvec_eq N ω hz
  have hzb := z_mul_qformC_gramC_ne_zero m N ω hz
  have hzne : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hbne : R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) ≠ 0 :=
    right_ne_zero_of_mul hzb
  have hsub : (1 : ℂ) - R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω)
      = -(z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) := by
    rw [ha]; ring
  have hq2 : R4C.qform2C (m.W0 N ω) z (m.gvec N ω)
      = R4C.cform2C (m.W0 N ω) z (m.gvec N ω) (m.gvec N ω) := rfl
  have hcgg : R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  have hc2gg : R4C.cform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  rw [hq2, m.W0_eq_gram_sub N ω,
    R4C.cform2C_sub_vecMulVec hW hz (m.gvec N ω) (m.gvec N ω) (m.gvec N ω), hcgg, hc2gg, hsub,
    ha, hA4]
  field_simp
  ring

/-- **A7.** The downdate on `cform2C x g`. A5 at `y = g`, with `cformC W z g g = a = 1 + z b`
and `cform2C W z g g = qform2C W z g = b + z b₂` (A4). -/
theorem cform2C_W0_gvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z : ℂ} (hz : 0 < z.im)
    (x : Fin (d N) → ℝ) :
    R4C.cform2C (m.W0 N ω) z x (m.gvec N ω)
      = R4C.cform2C (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
        - (R4C.cform2C (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
              * (1 + z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
            + R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
              * (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
                  + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))))
            / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
        + R4C.cformC (GenRMT.gram (m.Z N ω)) z x (m.gvec N ω)
            * (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
                + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
            * (1 + z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
            / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) ^ 2 := by
  have ha := m.qformC_gvec_eq N ω hz
  have hA4 := m.qform2C_gvec_eq N ω hz
  have hA5 := cform2C_W0_eq m N ω hz x (m.gvec N ω)
  have hcgg : R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  have hc2gg : R4C.cform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (m.gvec N ω)
      = R4C.qform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) := rfl
  rw [hA5, hcgg, hc2gg, ha, hA4]

end SpikedModel

namespace GenRMT

/-! F37 (2026-09-09): `cmat_smul` used to be repeated here (`private`); `Companion.lean:57`
now carries the public canonical copy, used directly as `R4C.cmat_smul`. -/

/-- The rescaled point `z d / p` has positive imaginary part when `z` does and `p`, `d` are
positive naturals: `z d / p = z (d / p : ℝ)`, a real positive multiple of `z`. Reused by `C1`
and `C2`. Public since F37 (2026-09-09), the canonical copy for `RMT/General/` (the private
twin of `FormsLimits.lean` is dropped). -/
theorem im_mul_natCast_div_natCast_pos {z : ℂ} (hz : 0 < z.im) {p d : ℕ}
    (hp : 0 < p) (hd : 0 < d) : 0 < (z * (d : ℂ) / (p : ℂ)).im := by
  have hr : z * (d : ℂ) / (p : ℂ) = z * (((d : ℝ) / (p : ℝ) : ℝ) : ℂ) := by
    push_cast
    ring
  rw [hr, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, mul_zero, zero_add]
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hdR : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd
  exact mul_pos hz (div_pos hdR hpR)

/-- **C1.** The companion resolvent of `Y` is a rescaled `gram` resolvent of `Yᵀ` at the
rescaled point: `Gc(z) = (d / p) G_{Yᵀ}(z d / p)`, from `gramC_eq_smul_gram_transpose`
(`Defs.lean:94`) and the two-sided inverse characterization of the resolvent
(`ResolvDeriv.cmat_sub_mul_resolvC`, `ResolvDeriv.lean:196`). -/
theorem resolvC_gramC_eq {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (hp : 0 < p) (hd : 0 < d)
    {z : ℂ} (hz : 0 < z.im) :
    R4C.resolvC (gramC Y) z = ((d : ℂ) / p) • R4C.resolvC (gram Yᵀ) (z * d / p) := by
  have hpC : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hzdp : 0 < (z * (d : ℂ) / (p : ℂ)).im := im_mul_natCast_div_natCast_pos hz hp hd
  have hAB : R4C.cmat (gramC Y) = ((p : ℂ) / d) • R4C.cmat (gram Yᵀ) := by
    rw [gramC_eq_smul_gram_transpose Y hp, R4C.cmat_smul]
    push_cast
    rfl
  have hzw : z = ((p : ℂ) / d) * (z * d / p) := by
    field_simp
  have hfact : R4C.cmat (gramC Y) - z • (1 : Matrix (Fin p) (Fin p) ℂ)
      = ((p : ℂ) / d) • (R4C.cmat (gram Yᵀ) - (z * d / p) • (1 : Matrix (Fin p) (Fin p) ℂ)) := by
    rw [hAB, smul_sub, smul_smul, ← hzw]
  have hone : (R4C.cmat (gramC Y) - z • (1 : Matrix (Fin p) (Fin p) ℂ)) *
      (((d : ℂ) / p) • R4C.resolvC (gram Yᵀ) (z * d / p)) = 1 := by
    rw [hfact, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    have hcancel : (p : ℂ) / d * ((d : ℂ) / p) = 1 := by field_simp
    rw [hcancel, one_smul]
    exact ResolvDeriv.cmat_sub_mul_resolvC (gram_isHermitian Yᵀ) hzdp.ne'
  exact Matrix.inv_eq_right_inv hone

/-- **C2, first half.** The companion `cformC` as a rescaled `gram` form of `Yᵀ`, from C1. -/
theorem cformC_gramC_eq {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (hp : 0 < p) (hd : 0 < d)
    {z : ℂ} (hz : 0 < z.im) (x y : Fin p → ℝ) :
    R4C.cformC (gramC Y) z x y = ((d : ℂ) / p) * R4C.cformC (gram Yᵀ) (z * d / p) x y := by
  have hq : R4C.cformC (gramC Y) z x y
      = R4C.cvec x ⬝ᵥ (R4C.resolvC (gramC Y) z *ᵥ R4C.cvec y) := rfl
  have hq2 : R4C.cvec x ⬝ᵥ (R4C.resolvC (gram Yᵀ) (z * d / p) *ᵥ R4C.cvec y)
      = R4C.cformC (gram Yᵀ) (z * d / p) x y := rfl
  rw [hq, resolvC_gramC_eq Y hp hd hz, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, hq2]

/-- **C2, second half.** The companion `cform2C` as a rescaled `gram` form of `Yᵀ`, from C1. -/
theorem cform2C_gramC_eq {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (hp : 0 < p) (hd : 0 < d)
    {z : ℂ} (hz : 0 < z.im) (x y : Fin p → ℝ) :
    R4C.cform2C (gramC Y) z x y
      = ((d : ℂ) / p) ^ 2 * R4C.cform2C (gram Yᵀ) (z * d / p) x y := by
  have hq : R4C.cform2C (gramC Y) z x y
      = R4C.cvec x ⬝ᵥ ((R4C.resolvC (gramC Y) z * R4C.resolvC (gramC Y) z) *ᵥ R4C.cvec y) := rfl
  have hq2 : R4C.cvec x ⬝ᵥ
      ((R4C.resolvC (gram Yᵀ) (z * d / p) * R4C.resolvC (gram Yᵀ) (z * d / p)) *ᵥ R4C.cvec y)
      = R4C.cform2C (gram Yᵀ) (z * d / p) x y := rfl
  rw [hq, resolvC_gramC_eq Y hp hd hz, Matrix.smul_mul, Matrix.mul_smul, smul_smul,
    Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, hq2]
  ring

/-- **C3.** The Lipschitz bound of `cformC` in `z`, uniform in the Hermitian matrix. Template
`R4C.norm_cformC_le` (`RMT/R4C.lean:358` to `:376`): two `cformC_eq_sum`, the termwise identity
`(λ - z₁)⁻¹ - (λ - z₂)⁻¹ = (z₁ - z₂) ((λ - z₁)⁻¹ (λ - z₂)⁻¹)`, `R4C.norm_inv_eigenvalue_sub_le`
on each factor, then `R4C.sum_abs_coords_le`. -/
theorem norm_cformC_sub_cformC_le {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z₁ z₂ : ℂ} (hz₁ : 0 < z₁.im) (hz₂ : 0 < z₂.im) (x y : Fin d → ℝ) :
    ‖R4C.cformC W z₁ x y - R4C.cformC W z₂ x y‖
      ≤ Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) * ‖z₁ - z₂‖ / (z₁.im * z₂.im) := by
  have hz₁' : z₁.im ≠ 0 := ne_of_gt hz₁
  have hz₂' : z₂.im ≠ 0 := ne_of_gt hz₂
  rw [R4C.cformC_eq_sum hW hz₁', R4C.cformC_eq_sum hW hz₂', ← Finset.sum_sub_distrib]
  have hterm : ∀ a : Fin d,
      ((hW.eigenvalues a : ℂ) - z₁)⁻¹
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)
        - ((hW.eigenvalues a : ℂ) - z₂)⁻¹
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)
      = (z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ) := by
    intro a
    have hne1 : ((hW.eigenvalues a : ℂ) - z₁) ≠ 0 := R4C.eigenvalue_sub_ne_zero hW hz₁' a
    have hne2 : ((hW.eigenvalues a : ℂ) - z₂) ≠ 0 := R4C.eigenvalue_sub_ne_zero hW hz₂' a
    field_simp
    ring
  rw [Finset.sum_congr rfl fun a _ => hterm a]
  calc ‖∑ a, (z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖
      ≤ ∑ a, ‖(z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖ := norm_sum_le _ _
    _ ≤ ∑ a, ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹)
          * (|((R4.eigU hW)ᵀ *ᵥ x) a| * |((R4.eigU hW)ᵀ *ᵥ y) a|) := by
        refine Finset.sum_le_sum fun a _ => ?_
        rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul]
        have hb1 : ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹‖ ≤ (z₁.im)⁻¹ :=
          R4C.norm_inv_eigenvalue_sub_le hW hz₁ a
        have hb2 : ‖((hW.eigenvalues a : ℂ) - z₂)⁻¹‖ ≤ (z₂.im)⁻¹ :=
          R4C.norm_inv_eigenvalue_sub_le hW hz₂ a
        have hprodbound : ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹‖
            ≤ (z₁.im)⁻¹ * (z₂.im)⁻¹ := by
          rw [norm_mul]
          exact mul_le_mul hb1 hb2 (norm_nonneg _) (by positivity)
        exact mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_left hprodbound (norm_nonneg _)) (by positivity)
    _ = ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹)
          * ∑ a, |((R4.eigU hW)ᵀ *ᵥ x) a| * |((R4.eigU hW)ᵀ *ᵥ y) a| := by
        rw [Finset.mul_sum]
    _ ≤ ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹) * (Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y)) := by
        exact mul_le_mul_of_nonneg_left (R4C.sum_abs_coords_le hW x y) (by positivity)
    _ = Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) * ‖z₁ - z₂‖ / (z₁.im * z₂.im) := by
        field_simp

/-- **C4.** The Lipschitz bound of `cform2C` in `z`. From `cform2C_eq_sum`, the termwise
identity `r₁² - r₂² = (r₁ - r₂)(r₁ + r₂)` with C3's factor `r₁ - r₂ = (z₁ - z₂) r₁ r₂`, and
`norm_add_le` on `r₁ + r₂`. -/
theorem norm_cform2C_sub_cform2C_le {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z₁ z₂ : ℂ} (hz₁ : 0 < z₁.im) (hz₂ : 0 < z₂.im) (x y : Fin d → ℝ) :
    ‖R4C.cform2C W z₁ x y - R4C.cform2C W z₂ x y‖
      ≤ Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) * ‖z₁ - z₂‖ * ((z₁.im)⁻¹ + (z₂.im)⁻¹)
          / (z₁.im * z₂.im) := by
  have hz₁' : z₁.im ≠ 0 := ne_of_gt hz₁
  have hz₂' : z₂.im ≠ 0 := ne_of_gt hz₂
  rw [R4C.cform2C_eq_sum hW hz₁', R4C.cform2C_eq_sum hW hz₂', ← Finset.sum_sub_distrib]
  have hterm : ∀ a : Fin d,
      (((hW.eigenvalues a : ℂ) - z₁)⁻¹) ^ 2
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)
        - (((hW.eigenvalues a : ℂ) - z₂)⁻¹) ^ 2
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)
      = (z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ + ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ) := by
    intro a
    have hne1 : ((hW.eigenvalues a : ℂ) - z₁) ≠ 0 := R4C.eigenvalue_sub_ne_zero hW hz₁' a
    have hne2 : ((hW.eigenvalues a : ℂ) - z₂) ≠ 0 := R4C.eigenvalue_sub_ne_zero hW hz₂' a
    field_simp
    ring
  rw [Finset.sum_congr rfl fun a _ => hterm a]
  calc ‖∑ a, (z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ + ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖
      ≤ ∑ a, ‖(z₁ - z₂) * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((hW.eigenvalues a : ℂ) - z₁)⁻¹ + ((hW.eigenvalues a : ℂ) - z₂)⁻¹)
          * (((((R4.eigU hW)ᵀ *ᵥ x) a * ((R4.eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖ := norm_sum_le _ _
    _ ≤ ∑ a, ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹) * ((z₁.im)⁻¹ + (z₂.im)⁻¹)
          * (|((R4.eigU hW)ᵀ *ᵥ x) a| * |((R4.eigU hW)ᵀ *ᵥ y) a|) := by
        refine Finset.sum_le_sum fun a _ => ?_
        have hb1 : ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹‖ ≤ (z₁.im)⁻¹ :=
          R4C.norm_inv_eigenvalue_sub_le hW hz₁ a
        have hb2 : ‖((hW.eigenvalues a : ℂ) - z₂)⁻¹‖ ≤ (z₂.im)⁻¹ :=
          R4C.norm_inv_eigenvalue_sub_le hW hz₂ a
        have hsumbound : ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹ + ((hW.eigenvalues a : ℂ) - z₂)⁻¹‖
            ≤ (z₁.im)⁻¹ + (z₂.im)⁻¹ := le_trans (norm_add_le _ _) (add_le_add hb1 hb2)
        have hprodbound : ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹‖
            ≤ (z₁.im)⁻¹ * (z₂.im)⁻¹ := by
          rw [norm_mul]; exact mul_le_mul hb1 hb2 (norm_nonneg _) (by positivity)
        have hABbound : ‖z₁ - z₂‖
              * ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹ * ((hW.eigenvalues a : ℂ) - z₂)⁻¹‖
              * ‖((hW.eigenvalues a : ℂ) - z₁)⁻¹ + ((hW.eigenvalues a : ℂ) - z₂)⁻¹‖
            ≤ ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹) * ((z₁.im)⁻¹ + (z₂.im)⁻¹) :=
          mul_le_mul (mul_le_mul_of_nonneg_left hprodbound (norm_nonneg _)) hsumbound
            (norm_nonneg _) (by positivity)
        rw [norm_mul, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul]
        exact mul_le_mul hABbound (le_refl _) (by positivity) (by positivity)
    _ = ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹) * ((z₁.im)⁻¹ + (z₂.im)⁻¹)
          * ∑ a, |((R4.eigU hW)ᵀ *ᵥ x) a| * |((R4.eigU hW)ᵀ *ᵥ y) a| := by
        rw [Finset.mul_sum]
    _ ≤ ‖z₁ - z₂‖ * ((z₁.im)⁻¹ * (z₂.im)⁻¹) * ((z₁.im)⁻¹ + (z₂.im)⁻¹)
          * (Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y)) := by
        exact mul_le_mul_of_nonneg_left (R4C.sum_abs_coords_le hW x y) (by positivity)
    _ = Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) * ‖z₁ - z₂‖ * ((z₁.im)⁻¹ + (z₂.im)⁻¹)
          / (z₁.im * z₂.im) := by
        field_simp

end GenRMT

end StackedSVD
