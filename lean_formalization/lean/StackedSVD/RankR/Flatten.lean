/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.General
import StackedSVD.RankR.Weighted

/-!
# The flatten bridge: general `r_i` block objects as `r_i = 1` objects on `Fin (rtot rk)`

STATUS 2026-09-01: proved, 0 `sorry` (task B0 of `notes/archive/rankr_plan_B.md`). The plan itself
is the review note; the audit is `notes/archive/audit_rankr_plan_B_2026-09-01.md`.

`RankR/Defs.lean` and `RankR/Weighted.lean` state `BR`, `Dmat`, `AbetaR`, `limitR`, `AbetaRW`,
`BRW`, `limitRW`, `optWR`, `abetaRSqrt`, `DcongR`, `limitROpt` for an arbitrary table count
`Fin M'`. Read at `M' := rtot rk`, `β := betaFlat β`, `R := Rcol R`, each of these objects is
literally the general-`r_i` block object of `RankR/General.lean`, so every bridge below closes
by `rfl`. `Rcol R p` is the flattened form of the family
`R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ`: row `k` of column `p` reads entry `k` of
`R (blk p).1` in column `(blk p).2`.

The two lemmas that are not `rfl`, `norm_Rcol` and `inner_Rcol`, turn `‖R i‖ = 1` and the
paper's inner product form into the shape the block bookkeeping of `RankR/General.lean`
(`ABlock_intra` and its neighbors) consumes.
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

variable {M r : ℕ} {rk : Fin M → ℕ}

/-- Column `p` of the flattened alignment family: row `k` is `R (blk p).1 k (blk p).2`. The
flattened form of `R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ`, read by the `r_i = 1`
objects of `RankR/Defs.lean` and `RankR/Weighted.lean` at table count `rtot rk`. -/
noncomputable def Rcol (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (p : Fin (rtot rk)) :
    EuclideanSpace ℝ (Fin r) :=
  WithLp.toLp 2 fun k => R (blk p).1 k (blk p).2

/-- `B_R` at general `r_i` is `BR` of `RankR/Defs.lean`, read at table count `rtot rk`. -/
theorem BBlock_eq_BR (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    BBlock β R = BR (betaFlat β) (Rcol R) := rfl

/-- `D` at general `r_i` is `Dmat` of `RankR/Defs.lean`, read at table count `rtot rk`. -/
theorem DBlock_eq_Dmat (β : (i : Fin M) → Fin (rk i) → ℝ) :
    DBlock (rk := rk) β = Dmat (betaFlat β) := rfl

/-- `A_{β,R}` at general `r_i` is `AbetaR` of `RankR/Defs.lean`. -/
theorem ABlock_eq_AbetaR (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    ABlock β R = AbetaR (betaFlat β) (Rcol R) := rfl

/-- `W A_{β,R} Wᵀ` at general `r_i` is `AbetaRW` of `RankR/Weighted.lean`. -/
theorem ABlockW_eq_AbetaRW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    ABlockW W β R = AbetaRW W (betaFlat β) (Rcol R) := rfl

/-- `W B_R` at general `r_i` is `BRW` of `RankR/Weighted.lean`. -/
theorem BBlockW_eq_BRW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    BBlockW W β R = BRW W (betaFlat β) (Rcol R) := rfl

/-- The unweighted limit at general `r_i` is `limitR` of `RankR/Defs.lean`. -/
theorem limitRG_eq_limitR (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    limitRG β R = limitR (betaFlat β) (Rcol R) := rfl

/-- The weighted limit at general `r_i` is `limitRW` of `RankR/Weighted.lean`. -/
theorem limitRGW_eq_limitRW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    limitRGW W β R = limitRW W (betaFlat β) (Rcol R) := rfl

/-- The optimal weight `D^{-1/2}` at general `r_i` is `optWR` of `RankR/Weighted.lean`. -/
theorem optWG_eq_optWR (β : (i : Fin M) → Fin (rk i) → ℝ) :
    optWG (rk := rk) β = optWR (betaFlat β) := rfl

/-- `A_{β,R}^{1/2}` at general `r_i` is `abetaRSqrt` of `RankR/Weighted.lean`. -/
theorem ABlockSqrt_eq_abetaRSqrt (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    ABlockSqrt β R = abetaRSqrt (betaFlat β) (Rcol R) := rfl

/-- `A_{β,R}^{-1/2} D A_{β,R}^{-1/2}` at general `r_i` is `DcongR` of `RankR/Weighted.lean`. -/
theorem DcongG_eq_DcongR (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    DcongG β R = DcongR (betaFlat β) (Rcol R) := rfl

/-- `L⋆` at general `r_i` is `limitROpt` of `RankR/Weighted.lean`. -/
theorem limitOptG_eq_limitROpt (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ)
    (hr : r ≤ Fintype.card (Fin (rtot rk))) :
    limitOptG β R hr = limitROpt (betaFlat β) (Rcol R) hr := rfl

/-- The inner product of two flattened columns is the paper's `(R_i)ᵀ R_{i'}` read at the two
block indices, the form `abetaR_apply` and `ABlock_apply` both use. -/
theorem inner_Rcol (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (p q : Fin (rtot rk)) :
    ⟪Rcol R p, Rcol R q⟫_ℝ = ((R (blk p).1)ᵀ * R (blk q).1) (blk p).2 (blk q).2 := by
  rw [real_inner_eq_dotProduct, Matrix.mul_apply]
  simp [Rcol, dotProduct, Matrix.transpose_apply]

/-- A column of `Rcol R` is a unit vector once its table's alignment matrix has orthonormal
columns, the hypothesis `UnalignedModelR.hR` supplies at every table. -/
theorem norm_Rcol (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (p : Fin (rtot rk))
    (hR : (R (blk p).1)ᵀ * R (blk p).1 = 1) : ‖Rcol R p‖ = 1 := by
  have h : ‖Rcol R p‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, inner_Rcol, hR, Matrix.one_apply_eq]
  rw [← Real.sqrt_sq (norm_nonneg (Rcol R p)), h, Real.sqrt_one]

end StackedSVD
