/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Spectral
import StackedSVD.RMT.R4C
import StackedSVD.RMT.ResolvDeriv
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.Companion

/-!
# Item Deloc at a general noise law: the deterministic layer

This file carries the deterministic (model-free, measure-free) input of item R6' at a
general noise law, the subcritical `align` field
(`SpikedModel.align_tendstoInProb_of_subcritical_general`), the frozen statement of the plan
(`notes/archive/prop_single_table_general.md`, section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35).

The Gaussian proof of item R6' (`RMT/R6.lean`) splits `v = a q + r` with `q` the spike plus
noise vector of the R0 split, bounds the `q` part by the secular equation, and bounds the
`r` part by rotation invariance. A general law has no rotation invariance, so the `r` part
is bounded instead by a **spectral window**: the top eigenvalue of the Gram matrix sits in a
window `[x - η, x + η]` around the bulk edge, and the squared norm of the spectral projector
on that window is at most `2 η Im (w ᵀ G(x + i η) w)`, a resolvent quantity that the general
law controls.

## Content

1. The model-free sections of `RMT/R6.lean`, copied here because `RMT/R6.lean` is Gaussian
   and `RMT/General/` may not import it: `overlap` is homogeneous and subadditive
   (`overlap_smul`, `overlap_smul_add_le`), the decomposition `v = a q + r`
   (`aOf`, `rvecOf`, `overlap_le_decomp`) and the secular bound on the `q` part
   (`overlap_q_le`).
2. `normSq_specProj_mono`: the squared norm of a spectral projector grows with the set.
3. `overlap_le_specProj_Icc`: the overlap is at most the window projector at any window that
   contains the top eigenvalue.
4. `normSq_specProj_eq_sum_eigU`: the window projector in the coordinates of `R4.eigU`.
5. `normSq_specProj_Icc_le`: the window bound
   `‖P_[x-η,x+η] w‖² ≤ 2 η Im (w ᵀ G(x + i η) w)`.
6. `cformC_add_vecMulVec_eq`: the rank-one **update** twin of `cformC_sub_vecMulVec`
   (`RMT/General/Companion.lean:206`), in the form that leaves the updated resolvent on the
   right side.

Everything here is deterministic. Nothing is asymptotic, and no measure appears.
-/

open Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD
namespace GenRMT
namespace Deloc

/-! ### 1. Elementary facts about `overlap`

Copied verbatim from `RMT/R6.lean:50-82` (section `Elementary`), which this file may not
import. -/

section Elementary

variable {n d : ℕ}

/-- `overlap` is homogeneous of degree two in the test direction.
Copy of `R6.overlap_smul` (`RMT/R6.lean:53`). -/
theorem overlap_smul (X : Matrix (Fin n) (Fin d) ℝ) (a : ℝ) (w : EuclideanSpace ℝ (Fin d)) :
    overlap X (a • w) = a ^ 2 * overlap X w := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) (a • w)‖ ^ 2
      = a ^ 2 * ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2
  rw [map_smul, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]

/-- Copy of `R6.overlap_zero` (`RMT/R6.lean:60`). -/
theorem overlap_zero (X : Matrix (Fin n) (Fin d) ℝ) :
    overlap X (0 : EuclideanSpace ℝ (Fin d)) = 0 := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) 0‖ ^ 2 = 0
  rw [map_zero, norm_zero]
  norm_num

/-- The triangle inequality for the top projector, in squared form.
Copy of `R6.overlap_smul_add_le` (`RMT/R6.lean:68`). -/
theorem overlap_smul_add_le (X : Matrix (Fin n) (Fin d) ℝ) (a : ℝ)
    (y r : EuclideanSpace ℝ (Fin d)) :
    overlap X (a • y + r) ≤ 2 * (a ^ 2 * overlap X y) + 2 * overlap X r := by
  set P := topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) with hP
  have h1 : ‖P (a • y + r)‖ ≤ |a| * ‖P y‖ + ‖P r‖ := by
    rw [map_add, map_smul]
    refine (norm_add_le _ _).trans ?_
    rw [norm_smul, Real.norm_eq_abs]
  change ‖P (a • y + r)‖ ^ 2 ≤ 2 * (a ^ 2 * ‖P y‖ ^ 2) + 2 * ‖P r‖ ^ 2
  calc ‖P (a • y + r)‖ ^ 2 ≤ (|a| * ‖P y‖ + ‖P r‖) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ 2 * (a ^ 2 * ‖P y‖ ^ 2) + 2 * ‖P r‖ ^ 2 := by
        nlinarith [sq_nonneg (|a| * ‖P y‖ - ‖P r‖), sq_abs a]

end Elementary

/-! ### 2. The orthogonal decomposition of `v` along `q`

Copied verbatim from `RMT/R6.lean:86-167` (section `Decomposition`). -/

section Decomposition

variable {d : ℕ}

/-- The coefficient of `q` in `v`. Junk value `0` at `q = 0`.
Copy of `R6.aOf` (`RMT/R6.lean:94`). -/
noncomputable def aOf (v q : Fin d → ℝ) : ℝ := (v ⬝ᵥ q) / (q ⬝ᵥ q)

/-- The part of `v` orthogonal to `q`. Copy of `R6.rvecOf` (`RMT/R6.lean:97`). -/
noncomputable def rvecOf (v q : Fin d → ℝ) : Fin d → ℝ := v - aOf v q • q

/-- Nonnegativity of the dot product with itself.
Copy of `R6.dotProduct_self_nonneg'` (`RMT/R6.lean:101`). -/
theorem dotProduct_self_nonneg' (q : Fin d → ℝ) : 0 ≤ q ⬝ᵥ q :=
  Finset.sum_nonneg fun _ _ => mul_self_nonneg _

/-- Copy of `R6.aOf_mul` (`RMT/R6.lean:104`). -/
theorem aOf_mul (v q : Fin d → ℝ) : aOf v q * (q ⬝ᵥ q) = v ⬝ᵥ q := by
  rcases eq_or_ne (q ⬝ᵥ q) 0 with h0 | h0
  · have hq : q = 0 := dotProduct_self_eq_zero.mp h0
    subst hq
    simp [aOf]
  · rw [aOf, div_mul_cancel₀ _ h0]

/-- Copy of `R6.dotProduct_rvecOf` (`RMT/R6.lean:112`). -/
theorem dotProduct_rvecOf (v q : Fin d → ℝ) : rvecOf v q ⬝ᵥ q = 0 := by
  rw [rvecOf, sub_dotProduct, smul_dotProduct, smul_eq_mul, aOf_mul, sub_self]

/-- Copy of `R6.eq_smul_add_rvecOf` (`RMT/R6.lean:115`). -/
theorem eq_smul_add_rvecOf (v q : Fin d → ℝ) : v = aOf v q • q + rvecOf v q := by
  rw [rvecOf]
  abel

/-- Copy of `R6.dotProduct_rvecOf_self_le` (`RMT/R6.lean:119`). -/
theorem dotProduct_rvecOf_self_le (v q : Fin d → ℝ) :
    rvecOf v q ⬝ᵥ rvecOf v q ≤ v ⬝ᵥ v := by
  have hexp : rvecOf v q ⬝ᵥ rvecOf v q
      = v ⬝ᵥ v - 2 * (aOf v q * (v ⬝ᵥ q)) + aOf v q ^ 2 * (q ⬝ᵥ q) := by
    rw [rvecOf]
    simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul]
    rw [dotProduct_comm q v]
    ring
  have hac : aOf v q * (v ⬝ᵥ q) = aOf v q ^ 2 * (q ⬝ᵥ q) := by
    rw [← aOf_mul v q]; ring
  rw [hexp, hac]
  nlinarith [mul_nonneg (sq_nonneg (aOf v q)) (dotProduct_self_nonneg' q)]

/-- Copy of `R6.aOf_sq_mul_le` (`RMT/R6.lean:133`). -/
theorem aOf_sq_mul_le (v q : Fin d → ℝ) : aOf v q ^ 2 * (q ⬝ᵥ q) ≤ v ⬝ᵥ v := by
  have h := dotProduct_rvecOf_self_le v q
  have hnn : 0 ≤ rvecOf v q ⬝ᵥ rvecOf v q := dotProduct_self_nonneg' _
  have hexp : rvecOf v q ⬝ᵥ rvecOf v q
      = v ⬝ᵥ v - aOf v q ^ 2 * (q ⬝ᵥ q) := by
    have hexp' : rvecOf v q ⬝ᵥ rvecOf v q
        = v ⬝ᵥ v - 2 * (aOf v q * (v ⬝ᵥ q)) + aOf v q ^ 2 * (q ⬝ᵥ q) := by
      rw [rvecOf]
      simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul]
      rw [dotProduct_comm q v]
      ring
    have hac : aOf v q * (v ⬝ᵥ q) = aOf v q ^ 2 * (q ⬝ᵥ q) := by
      rw [← aOf_mul v q]; ring
    rw [hexp', hac]; ring
  linarith [hexp ▸ hnn]

/-- Copy of `R6.toLp_eq_smul_add` (`RMT/R6.lean:151`). -/
theorem toLp_eq_smul_add (v q : Fin d → ℝ) :
    (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))
      = aOf v q • WithLp.toLp 2 q + WithLp.toLp 2 (rvecOf v q) := by
  apply WithLp.ofLp_injective
  change v = aOf v q • q + rvecOf v q
  exact eq_smul_add_rvecOf v q

/-- Copy of `R6.norm_toLp_rvecOf_le` (`RMT/R6.lean:158`). -/
theorem norm_toLp_rvecOf_le (v q : Fin d → ℝ) (hv : v ⬝ᵥ v ≤ 1) :
    ‖(WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d))‖ ≤ 1 := by
  have h1 : ‖(WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d))‖ ^ 2
      = rvecOf v q ⬝ᵥ rvecOf v q := by
    rw [← real_inner_self_eq_norm_sq]
    exact inner_euclidean_eq_dotProduct _ _
  nlinarith [norm_nonneg (WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d)),
    dotProduct_rvecOf_self_le v q, h1]

/-- **Step 1.** The decomposition bound. Copy of `R6.overlap_le_decomp`
(`RMT/R6.lean:167`). -/
theorem overlap_le_decomp {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) (v q : Fin d → ℝ) :
    overlap X (WithLp.toLp 2 v)
      ≤ 2 * (aOf v q ^ 2 * overlap X (WithLp.toLp 2 q))
        + 2 * overlap X (WithLp.toLp 2 (rvecOf v q)) := by
  rw [toLp_eq_smul_add v q]
  exact overlap_smul_add_le X _ _ _

end Decomposition

/-! ### 3. The deterministic bound on `overlap X q`

Copied from `RMT/R6.lean:171-274` (section `QBound`). Two names that the Gaussian file takes
from banned modules are reproved here privately: `R2.dotProduct_ofLp_self` and
`Symmetry.lamMax_mem_eigSet`. -/

section QBound

variable {n d : ℕ}

/-- `‖x‖ = 1` reads as `x ⬝ᵥ x = 1`. Private copy of `R2.dotProduct_ofLp_self`
(`RMT/R2.lean`), which `RMT/General/` may not import. -/
private theorem dotProduct_ofLp_self {x : EuclideanSpace ℝ (Fin d)} (hx : ‖x‖ = 1) :
    WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = 1 := by
  have h : ⟪x, x⟫_ℝ = WithLp.ofLp x ⬝ᵥ WithLp.ofLp x := inner_euclidean_eq_dotProduct x x
  rw [real_inner_self_eq_norm_sq, hx] at h
  simpa using h.symm

/-- `lamMax` is attained by an eigenvector. Private copy of `Symmetry.lamMax_mem_eigSet`
(`RMT/Symmetry.lean:166`), stated without the `eigSet` wrapper. -/
private theorem exists_eigvec_lamMax {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hd : 0 < d) :
    ∃ x : EuclideanSpace ℝ (Fin d), x ≠ 0 ∧ toOp A x = lamMax A hA • x := by
  have hcard : (0 : ℕ) < Fintype.card (Fin d) := by simpa using hd
  set i : Fin (Fintype.card (Fin d)) := ⟨0, hcard⟩ with hi
  refine ⟨(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, ?_, ?_⟩
  · have h1 : ‖(symmOp hA).eigenvectorBasis (finrank_euclideanSpace (ι := Fin d)) i‖ = 1 :=
      ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).orthonormal.1 i
    intro h
    rw [h, norm_zero] at h1
    exact zero_ne_one h1
  · rw [apply_eigvec hA i, lamMax, dif_pos hd]

/-- **Step 2, the secular case.** At a secular root above `lamMax W₀` the overlap with `q`
is the inverse of the squared resolvent form. Copy of `R6.overlap_q_eq_inv`
(`RMT/R6.lean:177`). -/
theorem overlap_q_eq_inv {W₀ : Matrix (Fin d) (Fin d) ℝ} (hW₀ : W₀.IsHermitian)
    {q : Fin d → ℝ} (hq : q ≠ 0) {lam : ℝ} (hlam : lamMax W₀ hW₀ < lam)
    (e : R4.secular W₀ q lam = 0) (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    ‖topProj (W₀ + Matrix.vecMulVec q q) hA (WithLp.toLp 2 q)‖ ^ 2
      = 1 / R4.qform2 W₀ lam q := by
  rw [R4.topProj_norm_sq hW₀ hq hlam e hA (WithLp.toLp 2 q)]
  have h2 : q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q) = -1 := by
    have he : 1 + R4.qform W₀ lam q = 0 := e
    have : R4.qform W₀ lam q = q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q) := rfl
    linarith [this ▸ he]
  change (q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q)) ^ 2 / _ = _
  rw [h2]
  norm_num
  rfl

/-- **Step 2.** On the simple event the overlap with `q` is at most the inverse of the
squared resolvent form at any `z₀` above the top eigenvalue of `Xᵀ X`. Copy of
`R6.overlap_q_le` (`RMT/R6.lean:200`). -/
theorem overlap_q_le (X : Matrix (Fin n) (Fin d) ℝ) {W₀ : Matrix (Fin d) (Fin d) ℝ}
    (hW₀ : W₀.IsHermitian) {q : Fin d → ℝ}
    (hgram : Xᵀ * X = W₀ + Matrix.vecMulVec q q) (hd : 0 < d)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {z₀ : ℝ} (hz₀ : gramLamMax X ≤ z₀) (hpos : 0 < R4.qform2 W₀ z₀ q) :
    overlap X (WithLp.toLp 2 q) ≤ 1 / R4.qform2 W₀ z₀ q := by
  have hA : (W₀ + Matrix.vecMulVec q q).IsHermitian := by
    rw [← hgram]; exact isHermitian_transpose_mul_self X
  have hgl : gramLamMax X = lamMax (W₀ + Matrix.vecMulVec q q) hA :=
    lamMax_congr hgram (isHermitian_transpose_mul_self X) hA
  have hle : lamMax W₀ hW₀ ≤ gramLamMax X := by
    rw [hgl]; exact R4.lamMax_le_lamMax_vecMulVec hW₀ hA
  rcases lt_or_eq_of_le hle with hlt | heq
  · -- the secular case
    have hq : q ≠ 0 := by
      intro h0
      have hz : Matrix.vecMulVec q q = 0 := by
        ext i j; rw [h0]; simp
      have : gramLamMax X = lamMax W₀ hW₀ := by
        rw [hgl]
        exact (lamMax_congr (by rw [hz, add_zero]) hA hW₀)
      linarith
    have hspec : gramLamMax X ∈ spectrum ℝ (toOp (W₀ + Matrix.vecMulVec q q)) := by
      obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hd
      rw [R4.spectrum_toOp, hgl]
      exact hj ▸ hA.eigenvalues_mem_spectrum_real j
    have hsec : R4.secular W₀ q (gramLamMax X) = 0 :=
      (R4.secular_eq_zero_iff hW₀ hlt).mpr hspec
    have hov : overlap X (WithLp.toLp 2 q) = 1 / R4.qform2 W₀ (gramLamMax X) q := by
      change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) (WithLp.toLp 2 q)‖ ^ 2 = _
      rw [topProj_congr hgram (isHermitian_transpose_mul_self X) hA (WithLp.toLp 2 q)]
      exact overlap_q_eq_inv hW₀ hq hlt hsec hA
    rw [hov]
    have hanti := R4.qform2_antitoneOn hW₀ q (Set.mem_Ioi.2 hlt)
      (Set.mem_Ioi.2 (lt_of_lt_of_le hlt hz₀)) hz₀
    exact one_div_le_one_div_of_le hpos hanti
  · -- the degenerate case: the top eigenvector of `W₀` is orthogonal to `q`
    obtain ⟨x, hx0, hxe⟩ := exists_eigvec_lamMax hW₀ hd
    have hxnorm : ‖x‖ ≠ 0 := norm_ne_zero_iff.mpr hx0
    set φ : EuclideanSpace ℝ (Fin d) := ‖x‖⁻¹ • x with hφdef
    have hφn : ‖φ‖ = 1 := by
      rw [hφdef, norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hxnorm]
    have hxm : W₀ *ᵥ WithLp.ofLp x = lamMax W₀ hW₀ • WithLp.ofLp x :=
      (R4.mem_eigenspace_iff' W₀ (lamMax W₀ hW₀) x).mp
        (Module.End.mem_eigenspace_iff.mpr hxe)
    have hφe : W₀ *ᵥ WithLp.ofLp φ = lamMax W₀ hW₀ • WithLp.ofLp φ := by
      change W₀ *ᵥ (‖x‖⁻¹ • WithLp.ofLp x) = lamMax W₀ hW₀ • (‖x‖⁻¹ • WithLp.ofLp x)
      rw [Matrix.mulVec_smul, hxm, smul_comm]
    have hφφ : WithLp.ofLp φ ⬝ᵥ WithLp.ofLp φ = 1 := dotProduct_ofLp_self hφn
    have hquad : WithLp.ofLp φ ⬝ᵥ ((W₀ + Matrix.vecMulVec q q) *ᵥ WithLp.ofLp φ)
        = lamMax W₀ hW₀ + (q ⬝ᵥ WithLp.ofLp φ) ^ 2 := by
      rw [Matrix.add_mulVec, dotProduct_add, hφe, R4.vecMulVec_mulVec]
      simp only [dotProduct_smul, smul_eq_mul]
      rw [hφφ, dotProduct_comm (WithLp.ofLp φ) q]
      ring
    have hle2 : WithLp.ofLp φ ⬝ᵥ ((W₀ + Matrix.vecMulVec q q) *ᵥ WithLp.ofLp φ)
        ≤ lamMax (W₀ + Matrix.vecMulVec q q) hA * (WithLp.ofLp φ ⬝ᵥ WithLp.ofLp φ) :=
      R4.dotProduct_mulVec_le_lamMax hA _
    have hqφ : q ⬝ᵥ WithLp.ofLp φ = 0 := by
      rw [hquad, hφφ, mul_one, ← hgl, ← heq] at hle2
      nlinarith [sq_nonneg (q ⬝ᵥ WithLp.ofLp φ)]
    have hAφ : (Xᵀ * X) *ᵥ WithLp.ofLp φ = gramLamMax X • WithLp.ofLp φ := by
      rw [hgram, Matrix.add_mulVec, hφe, R4.vecMulVec_mulVec, hqφ, zero_smul, add_zero, heq]
    have hmem : φ ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
      change φ ∈ specSpace (Xᵀ * X) {lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)}
      rw [specSpace]
      simp only [Set.mem_singleton_iff, iSup_iSup_eq_left]
      exact (R4.mem_eigenspace_iff' _ _ _).mpr hAφ
    have hzero : overlap X (WithLp.toLp 2 q) = 0 := by
      rw [overlap_eq_inner_sq X (WithLp.toLp 2 q) hsimple hmem hφn,
        inner_euclidean_eq_dotProduct]
      change (WithLp.ofLp φ ⬝ᵥ q) ^ 2 = 0
      rw [dotProduct_comm, hqφ]
      norm_num
    rw [hzero]
    exact div_nonneg zero_le_one hpos.le

end QBound

/-! ### 4. The rank-one update of the complex resolvent

`RMT/General/Companion.lean:178-215` does this algebra for the **downdate** `W - g gᵀ`
through Sherman-Morrison. The **update** `W + q qᵀ` needs no Sherman-Morrison and no
denominator: the two resolvent identities `G_A (A - z) = 1` and `(W - z) G = 1` give
`G_A = G - G_A q qᵀ G` directly, at the price of leaving `G_A` on the right side. -/

section Update

variable {d : ℕ}

/-- The cast of the rank-one update, the `+` twin of `Companion.cmat_sub_vecMulVec`
(`RMT/General/Companion.lean:65`). -/
private theorem cmat_add_vecMulVec (W : Matrix (Fin d) (Fin d) ℝ) (g : Fin d → ℝ) :
    R4C.cmat (W + Matrix.vecMulVec g g)
      = R4C.cmat W + Matrix.vecMulVec (R4C.cvec g) (R4C.cvec g) := by
  ext i j
  simp only [R4C.cmat, R4C.cvec, Matrix.map_apply, Matrix.add_apply, Matrix.vecMulVec_apply]
  push_cast
  ring

/-! F37 (2026-09-09): `vecMulVec_self_mulVec` used to be repeated here (`private`);
`R4C.vecMulVec_self_mulVec` (`Companion.lean:73`, now public) is the canonical copy, used
directly. -/

/-- A rank-one update keeps a real symmetric matrix symmetric. -/
theorem isHermitian_add_vecMulVec {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    (q : Fin d → ℝ) : (W + Matrix.vecMulVec q q).IsHermitian := by
  refine hW.add ?_
  ext i j
  simp [Matrix.vecMulVec_apply, mul_comm]

/-- **The rank-one update identity for the bilinear form of the resolvent.** With
`A = W + q qᵀ`, `G_A = G - G_A q qᵀ G`, so
`x ᵀ G_A y = x ᵀ G y - (x ᵀ G_A q) (q ᵀ G y)`.
This is the deterministic input that transfers a limit for the resolvent of the Wishart
block `W₀` to the resolvent of the full Gram matrix `Xᵀ X = W₀ + q qᵀ` of the R0 split. -/
theorem cformC_add_vecMulVec_eq {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    {z : ℂ} (hz : 0 < z.im) (q x y : Fin d → ℝ) :
    R4C.cformC (W + Matrix.vecMulVec q q) z x y
      = R4C.cformC W z x y
        - R4C.cformC (W + Matrix.vecMulVec q q) z x q * R4C.cformC W z q y := by
  have hz' : z.im ≠ 0 := hz.ne'
  have hA : (W + Matrix.vecMulVec q q).IsHermitian := isHermitian_add_vecMulVec hW q
  have hkey : R4C.resolvC (W + Matrix.vecMulVec q q) z
      = R4C.resolvC W z
        - R4C.resolvC (W + Matrix.vecMulVec q q) z
            * Matrix.vecMulVec (R4C.cvec q) (R4C.cvec q) * R4C.resolvC W z := by
    have h1 : R4C.resolvC (W + Matrix.vecMulVec q q) z
        * (R4C.cmat (W + Matrix.vecMulVec q q) - z • (1 : Matrix (Fin d) (Fin d) ℂ)) = 1 :=
      ResolvDeriv.resolvC_mul_cmat_sub hA hz'
    have h2 : (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) * R4C.resolvC W z = 1 :=
      ResolvDeriv.cmat_sub_mul_resolvC hW hz'
    have h3 : R4C.cmat (W + Matrix.vecMulVec q q) - z • (1 : Matrix (Fin d) (Fin d) ℂ)
        = (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ))
          + Matrix.vecMulVec (R4C.cvec q) (R4C.cvec q) := by
      rw [cmat_add_vecMulVec]
      abel
    rw [h3, Matrix.mul_add] at h1
    have h4 : R4C.resolvC (W + Matrix.vecMulVec q q) z
        * (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ))
          = 1 - R4C.resolvC (W + Matrix.vecMulVec q q) z
              * Matrix.vecMulVec (R4C.cvec q) (R4C.cvec q) := by
      rw [eq_sub_iff_add_eq]
      exact h1
    calc R4C.resolvC (W + Matrix.vecMulVec q q) z
        = R4C.resolvC (W + Matrix.vecMulVec q q) z
            * ((R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) * R4C.resolvC W z) := by
          rw [h2, Matrix.mul_one]
      _ = R4C.resolvC (W + Matrix.vecMulVec q q) z
            * (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) * R4C.resolvC W z := by
          rw [Matrix.mul_assoc]
      _ = R4C.resolvC W z
            - R4C.resolvC (W + Matrix.vecMulVec q q) z
                * Matrix.vecMulVec (R4C.cvec q) (R4C.cvec q) * R4C.resolvC W z := by
          rw [h4, Matrix.sub_mul, Matrix.one_mul]
  have hexp : R4C.cformC (W + Matrix.vecMulVec q q) z x y
      = R4C.cvec x ⬝ᵥ (R4C.resolvC W z *ᵥ R4C.cvec y)
        - R4C.cvec x ⬝ᵥ ((R4C.resolvC (W + Matrix.vecMulVec q q) z
            * Matrix.vecMulVec (R4C.cvec q) (R4C.cvec q) * R4C.resolvC W z) *ᵥ R4C.cvec y) := by
    change R4C.cvec x ⬝ᵥ (R4C.resolvC (W + Matrix.vecMulVec q q) z *ᵥ R4C.cvec y) = _
    conv_lhs => rw [hkey]
    rw [Matrix.sub_mulVec, dotProduct_sub]
  rw [hexp]
  congr 1
  have hqy : R4C.cvec q ⬝ᵥ (R4C.resolvC W z *ᵥ R4C.cvec y) = R4C.cformC W z q y := rfl
  have hxq : R4C.cvec x ⬝ᵥ (R4C.resolvC (W + Matrix.vecMulVec q q) z *ᵥ R4C.cvec q)
      = R4C.cformC (W + Matrix.vecMulVec q q) z x q := rfl
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, R4C.vecMulVec_self_mulVec,
    Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul, hqy, hxq]
  ring

end Update

/-! ### 5. The spectral window

The general-law replacement for the rotation-invariance step of `RMT/R6.lean`: the overlap
with the top eigenvalue is at most the projector on any window that contains it, and the
window projector is bounded by the imaginary part of the resolvent at the window's center,
lifted by `η`. -/

section Window

variable {d : ℕ}

open scoped Classical in
/-- The squared norm of a spectral projector grows with the set: the eigenbasis sum of
`Spectral.normSq_specProj` is termwise monotone in `S`. -/
theorem normSq_specProj_mono (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    {S T : Set ℝ} (hST : S ⊆ T) (w : EuclideanSpace ℝ (Fin d)) :
    ‖specProj A S w‖ ^ 2 ≤ ‖specProj A T w‖ ^ 2 := by
  rw [normSq_specProj A hA S w, normSq_specProj A hA T w]
  refine Finset.sum_le_sum fun i _ => ?_
  by_cases h : hA.eigenvalues₀ i ∈ S
  · rw [if_pos h, if_pos (hST h)]
  · rw [if_neg h, show (0 : ℝ) ^ 2 = 0 from by norm_num]
    exact sq_nonneg _

/-- **The window step.** If the top eigenvalue of `Xᵀ X` is within `η` of `x`, then the
overlap is at most the squared norm of the projector on the window `[x - η, x + η]`. -/
theorem overlap_le_specProj_Icc {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) {x η : ℝ}
    (h : |gramLamMax X - x| ≤ η) (w : EuclideanSpace ℝ (Fin d)) :
    overlap X w ≤ ‖specProj (Xᵀ * X) (Set.Icc (x - η) (x + η)) w‖ ^ 2 := by
  have hsub : ({gramLamMax X} : Set ℝ) ⊆ Set.Icc (x - η) (x + η) := by
    intro t ht
    rw [Set.mem_singleton_iff] at ht
    subst ht
    rw [abs_le] at h
    exact ⟨by linarith [h.1], by linarith [h.2]⟩
  exact normSq_specProj_mono (Xᵀ * X) (isHermitian_transpose_mul_self X) hsub w
/-- The index equiv of `Matrix.IsHermitian.eigenvalues`: `Fin (Fintype.card (Fin d)) ≃ Fin d`.
-/
private noncomputable def eqvIdx (d : ℕ) : Fin (Fintype.card (Fin d)) ≃ Fin d :=
  Fintype.equivOfCardEq (Fintype.card_fin (Fintype.card (Fin d)))

/-- `eigenvalues₀` reindexed. Private copy of the private `Spectral.eigenvalues₀_eq_eigenvalues`
(`Spectral.lean:289`). -/
private theorem eigenvalues₀_eq_eigenvalues {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin d))) : hA.eigenvalues₀ k = hA.eigenvalues (eqvIdx d k) := by
  simp [Matrix.IsHermitian.eigenvalues, eqvIdx]

/-- The eigenvector basis of `Spectral.normSq_specProj`, reindexed to the basis of
`Matrix.IsHermitian.eigenvectorBasis`. -/
private theorem eigenvectorBasis_eqvIdx {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin d))) :
    (symmOp hA).eigenvectorBasis finrank_euclideanSpace k
      = hA.eigenvectorBasis (eqvIdx d k) := by
  simp [Matrix.IsHermitian.eigenvectorBasis, OrthonormalBasis.reindex_apply, eqvIdx]

/-- The `a`-th coordinate in the eigenbasis, as an inner product. -/
private theorem transpose_eigU_mulVec_apply {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (w : Fin d → ℝ) (a : Fin d) :
    ((R4.eigU hA)ᵀ *ᵥ w) a
      = ⟪hA.eigenvectorBasis a, (WithLp.toLp 2 w : EuclideanSpace ℝ (Fin d))⟫_ℝ := by
  rw [inner_euclidean_eq_dotProduct, WithLp.ofLp_toLp]
  simp [R4.eigU, Matrix.mulVec, dotProduct, Matrix.transpose_apply]

open scoped Classical in
/-- **The basis bridge.** The eigenbasis sum of `Spectral.normSq_specProj`, written in the
coordinates `(eigU hA)ᵀ w` and the eigenvalues `hA.eigenvalues` that the resolvent formulas
of `RMT/R4C.lean` use. -/
theorem normSq_specProj_eq_sum_eigU (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (S : Set ℝ) (w : Fin d → ℝ) :
    ‖specProj A S (WithLp.toLp 2 w)‖ ^ 2
      = ∑ a, (if hA.eigenvalues a ∈ S then ((R4.eigU hA)ᵀ *ᵥ w) a else 0) ^ 2 := by
  rw [normSq_specProj A hA S (WithLp.toLp 2 w)]
  refine Fintype.sum_equiv (eqvIdx d) _ _ fun k => ?_
  rw [eigenvalues₀_eq_eigenvalues hA k, eigenvectorBasis_eqvIdx hA k,
    transpose_eigU_mulVec_apply hA w (eqvIdx d k)]
/-- `Im (g ᵀ G(z) g) = Im z ∑ w_a / |λ_a - z|²`. Public since F37 (2026-09-09), the canonical
copy for `RMT/General/` (the private twin of `Trace.lean` is dropped; `Trace.lean` imports
this file to reach it, since this file cannot import `Trace.lean`). -/
theorem im_qformC_eq {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian) {z : ℂ}
    (hz : 0 < z.im) (g : Fin d → ℝ) :
    (R4C.qformC W z g).im
      = z.im * ∑ a, ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2
          / Complex.normSq ((hW.eigenvalues a : ℂ) - z) := by
  have hz' : z.im ≠ 0 := hz.ne'
  have hstep : (R4C.cformC W z g g).im
      = ∑ a, (z.im / Complex.normSq ((hW.eigenvalues a : ℂ) - z))
          * (((R4.eigU hW)ᵀ *ᵥ g) a) ^ 2 := by
    rw [R4C.cformC_eq_sum hW hz' g g, Complex.im_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Complex.mul_im]
    have hre0 : ((((R4.eigU hW)ᵀ *ᵥ g) a * ((R4.eigU hW)ᵀ *ᵥ g) a : ℝ) : ℂ).im = 0 := by simp
    have hre1 : ((((R4.eigU hW)ᵀ *ᵥ g) a * ((R4.eigU hW)ᵀ *ᵥ g) a : ℝ) : ℂ).re
        = ((R4.eigU hW)ᵀ *ᵥ g) a ^ 2 := by rw [Complex.ofReal_re]; ring
    rw [hre0, hre1, mul_zero, zero_add, Complex.inv_im]
    have hIm : ((hW.eigenvalues a : ℂ) - z).im = -z.im := by simp
    rw [hIm]
    ring
  change (R4C.cformC W z g g).im = _
  rw [hstep, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => by ring

open scoped Classical in
/-- The window bound at a window centered on `Re z` of half-width `Im z`. -/
private theorem normSq_specProj_Icc_le_aux (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    {z : ℂ} (hz : 0 < z.im) (w : Fin d → ℝ) :
    ‖specProj A (Set.Icc (z.re - z.im) (z.re + z.im)) (WithLp.toLp 2 w)‖ ^ 2
      ≤ 2 * z.im * (R4C.qformC A z w).im := by
  have hD : ∀ lam : ℝ, Complex.normSq ((lam : ℂ) - z) = (lam - z.re) ^ 2 + z.im ^ 2 := by
    intro lam
    have hre' : ((lam : ℂ) - z).re = lam - z.re := by simp
    have him' : ((lam : ℂ) - z).im = -z.im := by simp
    rw [Complex.normSq_apply, hre', him']
    ring
  rw [normSq_specProj_eq_sum_eigU A hA _ w, im_qformC_eq hA hz w]
  calc _ ≤ ∑ a, 2 * z.im * (z.im * (((R4.eigU hA)ᵀ *ᵥ w) a ^ 2
          / Complex.normSq ((hA.eigenvalues a : ℂ) - z))) := by
        refine Finset.sum_le_sum fun a _ => ?_
        rw [hD (hA.eigenvalues a)]
        have hDpos : 0 < (hA.eigenvalues a - z.re) ^ 2 + z.im ^ 2 := by positivity
        have hcnn : 0 ≤ ((R4.eigU hA)ᵀ *ᵥ w) a ^ 2 / ((hA.eigenvalues a - z.re) ^ 2 + z.im ^ 2) :=
          div_nonneg (sq_nonneg _) hDpos.le
        by_cases hmem : hA.eigenvalues a ∈ Set.Icc (z.re - z.im) (z.re + z.im)
        · rw [if_pos hmem]
          obtain ⟨h1, h2⟩ := hmem
          have hwin : (hA.eigenvalues a - z.re) ^ 2 + z.im ^ 2 ≤ 2 * z.im ^ 2 := by
            nlinarith [h1, h2]
          have hexp : 2 * z.im * (z.im * (((R4.eigU hA)ᵀ *ᵥ w) a ^ 2
                / ((hA.eigenvalues a - z.re) ^ 2 + z.im ^ 2)))
              = 2 * z.im ^ 2 * ((R4.eigU hA)ᵀ *ᵥ w) a ^ 2
                / ((hA.eigenvalues a - z.re) ^ 2 + z.im ^ 2) := by
            field_simp
          rw [hexp, le_div_iff₀ hDpos]
          nlinarith [mul_le_mul_of_nonneg_left hwin (sq_nonneg (((R4.eigU hA)ᵀ *ᵥ w) a))]
        · rw [if_neg hmem, show (0 : ℝ) ^ 2 = 0 from by norm_num]
          exact mul_nonneg (by linarith : (0 : ℝ) ≤ 2 * z.im) (mul_nonneg hz.le hcnn)
    _ = 2 * z.im * (z.im * ∑ a, ((R4.eigU hA)ᵀ *ᵥ w) a ^ 2
          / Complex.normSq ((hA.eigenvalues a : ℂ) - z)) := by
        rw [Finset.mul_sum, Finset.mul_sum]

open scoped Classical in
/-- **The window bound.** The squared norm of the projector on `[x - η, x + η]` is at most
`2 η Im (w ᵀ G(x + i η) w)`. This is the general-law replacement for the rotation-invariance
step of `RMT/R6.lean`: it caps the whole spectral window at once, so it needs neither a
density nor a simple top eigenvalue. -/
theorem normSq_specProj_Icc_le (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    {x η : ℝ} (hη : 0 < η) (w : Fin d → ℝ) :
    ‖specProj A (Set.Icc (x - η) (x + η)) (WithLp.toLp 2 w)‖ ^ 2
      ≤ 2 * η * (R4C.qformC A (x + η * Complex.I) w).im := by
  have hre : ((x : ℂ) + (η : ℂ) * Complex.I).re = x := by simp
  have him : ((x : ℂ) + (η : ℂ) * Complex.I).im = η := by simp
  have h := normSq_specProj_Icc_le_aux A hA
    (z := (x : ℂ) + (η : ℂ) * Complex.I) (by rw [him]; exact hη) w
  rw [hre, him] at h
  exact h

end Window

end Deloc
end GenRMT
end StackedSVD
