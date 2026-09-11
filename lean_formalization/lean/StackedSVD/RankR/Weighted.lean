/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Defs
import StackedSVD.LinAlg.Eigen
import StackedSVD.LinAlg.KyFan

/-!
# `thm:gen_rank_weight_svdstak`: weighted svdstack in the unaligned setting, at `r_i = 1`

STATUS 2026-08-30: proved, 0 `sorry` (task T5b of `notes/RANK_R_PLAN.md`, proofs under
decisions D15 and D16 of `notes/FLAGGED.md`; the Layer 2 corollaries of section 7 under task
`rank_r_weighted_gaussian`). The review note is `notes/archive/thm_gen_rank_weight_svdstak.md`; the
reports are `notes/archive/agent_reports/rank_r_t5b_proofs.md` and
`notes/archive/agent_reports/rank_r_weighted_gaussian.md`.

`main_paper.tex:889` defines `V̂_svdstack(W)` as the top `r` right singular vectors of
`W Ṽ`, with a **matrix** weight `W ∈ ℝ^{r̃ × r̃}`, and `thm:gen_rank_weight_svdstak`
(`main_paper.tex:893`) says that the optimum over all such `W` is attained by the diagonal
`W⋆ = diag(1/√(1 - β_ij²)) = D^{-1/2}`, with value
`L⋆ = r - ∑_{ℓ=1}^{r} λ_{r̃+1-ℓ}(A_{β,R}^{-1/2} D A_{β,R}^{-1/2})`.

This file takes the same `r_i = 1` slice as `RankR/Defs.lean`, so `r̃ = ∑_i r_i = M`, each
`R_i` is a unit vector of `ℝ^r`, and the per-table input is the proved rank-one
`SingleTableLaw`.

## Content

1. `AbetaRW W β R = W A_{β,R} Wᵀ`, `BRW W β R = W B_R`, and the weighted limit
   `limitRW W β R = tr((W B_R)ᵀ (specInvTop (W A_{β,R} Wᵀ) r) (W B_R))`. At `W = 1` these are
   the unweighted objects of `RankR/Defs.lean` (`abetaRW_one`, `limitRW_one`).
2. `BR β R = diag(β) R_stack`, with `Rstack` and the rank bridge
   `rank(∑_i R_i R_iᵀ) = rank(R_stack)` in `LinAlg/Eigen.lean`. The paper's rank condition
   `Rank(∑_i R_i R_iᵀ) = r` is `(Rstack R).rank = r`.
3. `optWR β = D^{-1/2} = diagonal (optW β)`, the paper's `W⋆`. Its defining property is
   `optWR_mul_Dmat_mul_optWR : W⋆ D W⋆ᵀ = 1`, and its entries are the rank-one optimal
   weights `optW` of `SVDStack/Defs.lean` (`main_paper.tex:521`).
4. `abetaRSqrt β R = A_{β,R}^{1/2}` from `CFC.sqrt`, `DcongR β R = A^{-1/2} D A^{-1/2}`, and
   the optimal value `limitROpt β R hr = r - ∑_{k < r} λ_{rev k}(A^{-1/2} D A^{-1/2})`.
   `Fin.rev` is the index form of `KyFan.lean`: `k = ℓ - 1` gives the paper's `λ_{r̃+1-ℓ}`.
5. `UnalignedModel.VtW`, `gramW`, `VtVW` and the weighted performance `perfRW`, the trace
   form of `‖V̂_svdstack(W)ᵀ V‖_F²` (the paper's own first display at `main_paper.tex:2016`,
   with `X^{(d)} = Wᵀ Q_r Λ_r^{-1/2}` and `X^{(d)} (X^{(d)})ᵀ = Wᵀ specInvTop(W Ṽ Ṽᵀ Wᵀ) W`).
6. The results: `limitRW_le_opt` (step 1 of the paper's proof, Ky Fan), `limitRW_optWR`
   (step 2, attainment, with no eigengap hypothesis), `thm_gen_rank_weight_svdstak_general`
   (the limit of `perfRW W` for one admissible `W`), `thm_gen_rank_weight_svdstak` (the
   theorem at `W⋆`), `thm_gen_rank_weight_svdstak_opt` (no admissible `W` beats `L⋆`) and
   `thm_gen_rank_weight_svdstak_max` (both halves in one declaration, with the convergence
   at every admissible `W`, which is what the paper claims). `topGap_optWR_of_rankBR` carries
   the eigengap check for `W⋆`, and
   `topGap_optWR_of_rank` states it in the paper's own hypotheses.
   `thm_gen_rank_weight_svdstak_gaussian` and `thm_gen_rank_weight_svdstak_max_gaussian`
   (section 7) are the Layer 2 forms of the last two: the `SingleTableLaw` hypothesis is
   replaced by the proportional regime `Regime (c i)` of each table.
7. The general spectral, rank and frame facts that Mathlib lacks on this pin moved to
   `LinAlg/Eigen.lean` in cleanup wave 2: `eigenvalues₀_inv`, `eigenvalues₀_add_one`,
   `eigenvalues₀_eq_of_charpoly_eq`, `eigenvalues₀_pos_of_lt_rank`, `exists_topFrame`,
   `trace_specInvTop_add_one`, `Rstack` and `rank_sum_vecMulVec`.

## Route of the proofs (T5b)

Step 1, `limitRW_le_opt`. Under `hgapW` and `hposW` the paper's frame
`X = Wᵀ Q_r(A_W) Λ_r^{-1/2}(A_W)` satisfies `Xᵀ A_{β,R} X = 1` and
`limitRW W β R = tr(Xᵀ B_R B_Rᵀ X) = r - tr(Xᵀ D X)`, by `A_{β,R} = B_R B_Rᵀ + D`
(the definition of `AbetaR`). Then `kyFan_min_congr` of `LinAlg/KyFan.lean`, with
`S := abetaRSqrt β R`, `hS := transpose_abetaRSqrt`, `hSA := sq_abetaRSqrt` and
`hSu := isUnit_det_abetaRSqrt`, gives
`∑_{k < r} λ_{rev k}(S⁻¹ D S⁻¹) ≤ tr(Xᵀ D X)`, which is the claim.

Step 2, `limitRW_optWR`. No eigengap hypothesis is needed. `optWR_mul_Dmat_mul_optWR` gives
`W⋆ A_{β,R} W⋆ᵀ = W⋆ B_R B_Rᵀ W⋆ᵀ + 1`, so the top `r` eigenvectors of `A_{W⋆}` are those of
`W⋆ B_R B_Rᵀ W⋆ᵀ` and `tr(Λ^{-1/2} Qᵀ W⋆ B_R B_Rᵀ W⋆ᵀ Q Λ^{-1/2}) = tr(1_r - Λ^{-1})`. The
eigenvalues of `A_{W⋆}` and of `A^{-1/2} D A^{-1/2}` are reciprocal in reverse order, because
`D^{1/2} A^{-1} D^{1/2} = (D^{1/2}A^{-1/2})(D^{1/2}A^{-1/2})ᵀ` has the spectrum of
`(D^{1/2}A^{-1/2})ᵀ(D^{1/2}A^{-1/2}) = A^{-1/2} D A^{-1/2}` (`main_paper.tex:2120`).

Steps 3 and 4, the two probabilistic statements. `gramW_eq` and `VtVW_eq` write the two
random matrices as `W (Ṽ Ṽᵀ) Wᵀ` and `W (Ṽ V)`, whose entries are polynomials in the entries
of `Ṽ Ṽᵀ` and `Ṽ V`. `gramR` and `VtV_tendsto` of `RankR/Defs.lean` give the entrywise
limits, and `continuousAt_trace_specInvTop` of `LinAlg/SpecProjPerturb.lean` composed with
that polynomial map is the continuous function that `TendstoInProbPi.comp_continuous`
consumes. `thm_gen_rank_weight_svdstak` is the general form at `W⋆`, with
`topGap_optWR_of_rank` for the gap, `limitRW_optWR` for the value, and
`abetaRW_optWR_posDef` for `hposW`.

`thm_gen_rank_weight_svdstak_general` through `_max_gaussian` moved to
`RankR/WeightedMain.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

/-! ### 1. The weighted matrices and the weighted limit

This file is the slice `r_i = 1`, where `r̃ = ∑_i r_i = M` (`main_paper.tex:766`), so the
paper's eigenvalue index `r̃ + 1 - ℓ` is `Fin.rev` at `ℓ = k + 1` on
`Fin (Fintype.card (Fin M))`, which is the index form of `LinAlg/KyFan.lean`. -/

variable {M r : ℕ}

/-- `A_{β,R,W} = W A_{β,R} Wᵀ`, the limit of the weighted Gram matrix `(W Ṽ)(W Ṽ)ᵀ`. The
paper writes `W A_{β,R} Wᵀ` throughout the proof (`main_paper.tex:2011`). -/
noncomputable def AbetaRW (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : Matrix (Fin M) (Fin M) ℝ :=
  W * AbetaR β R * Wᵀ

theorem isHermitian_AbetaRW (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : (AbetaRW W β R).IsHermitian := by
  have h := Matrix.isHermitian_mul_mul_conjTranspose (A := AbetaR β R) W (isHermitian_AbetaR β R)
  rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h

theorem abetaRW_one (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    AbetaRW 1 β R = AbetaR β R := by
  rw [AbetaRW, Matrix.one_mul, Matrix.transpose_one, Matrix.mul_one]

/-- `W B_R`, the limit of `(W Ṽ) V`. -/
noncomputable def BRW (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : Matrix (Fin M) (Fin r) ℝ :=
  W * BR β R

theorem brW_one (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    BRW 1 β R = BR β R := Matrix.one_mul _

/-- The limit of the weighted performance, in the trace form of `limitR`:
`tr((W B_R)ᵀ (specInvTop (W A_{β,R} Wᵀ) r) (W B_R))`. It is the paper's `‖Xᵀ B_R‖_F²` at
`X = Wᵀ Q_r(A_W) Λ_r^{-1/2}(A_W)` (`main_paper.tex:2029`). -/
noncomputable def limitRW (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : ℝ :=
  Matrix.trace ((BRW W β R)ᵀ * specInvTop (AbetaRW W β R) (isHermitian_AbetaRW W β R) r *
    BRW W β R)

theorem limitRW_one (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    limitRW 1 β R = limitR β R := by
  rw [limitRW, limitR, brW_one,
    specInvTop_congr_mat (isHermitian_AbetaRW 1 β R) (isHermitian_AbetaR β R)
      (abetaRW_one β R) r]

/-! ### 2. `B_R` and the paper's rank condition

`Rstack`, `Rstack_transpose_mul` and `rank_sum_vecMulVec` moved to `LinAlg/Eigen.lean` in
cleanup wave 2; only the bridge to `BR` stays here. -/

/-- `B_R = diag(β) R_stack` (`main_paper.tex:1951`). -/
theorem BR_eq_diagonal_mul_Rstack (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    BR β R = Matrix.diagonal β * Rstack R := by
  ext i k
  rw [BR, Matrix.of_apply, Matrix.diagonal_mul, Rstack, Matrix.of_apply]

/-! ### 3. The paper's weights `W⋆ = D^{-1/2}` -/

/-- `W⋆ = D^{-1/2} = diag(1/√(1 - β_i²))` (`main_paper.tex:895`). The diagonal entries are the
rank-one optimal weights `optW` of `SVDStack/Defs.lean` (`main_paper.tex:521`), which is the
paper's statement that the optimal weighting does not depend on `R_i`. -/
noncomputable def optWR (β : Fin M → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.diagonal (optW β)

theorem transpose_optWR (β : Fin M → ℝ) : (optWR β)ᵀ = optWR β :=
  Matrix.diagonal_transpose _

theorem optW_mul_one_sub_mul {β : Fin M → ℝ} {i : Fin M} (h : β i ^ 2 < 1) :
    optW β i * (1 - β i ^ 2) * optW β i = 1 := by
  have hpos : (0 : ℝ) < 1 - β i ^ 2 := by linarith
  have hsq : Real.sqrt (1 - β i ^ 2) * Real.sqrt (1 - β i ^ 2) = 1 - β i ^ 2 :=
    Real.mul_self_sqrt hpos.le
  have key : optW β i * (1 - β i ^ 2) * optW β i
      = (1 - β i ^ 2) / (Real.sqrt (1 - β i ^ 2) * Real.sqrt (1 - β i ^ 2)) := by
    rw [optW]
    ring
  rw [key, hsq, div_self hpos.ne']

/-- The defining property of `W⋆`: `D^{-1/2} D D^{-1/2} = I`. Step 2 of the paper's proof uses
it as `W⋆ A_{β,R} W⋆ᵀ = W⋆ B_R B_Rᵀ W⋆ᵀ + I` (`main_paper.tex:2103`). -/
theorem optWR_mul_Dmat_mul_optWR {β : Fin M → ℝ} (hβ : ∀ i, β i ^ 2 < 1) :
    optWR β * Dmat β * (optWR β)ᵀ = 1 := by
  rw [transpose_optWR, optWR, Dmat, Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal,
    ← Matrix.diagonal_one]
  congr 1
  funext i
  simpa using optW_mul_one_sub_mul (hβ i)

/-! ### 4. `A_{β,R}^{1/2}`, the congruence `A^{-1/2} D A^{-1/2}` and the optimal value

`Matrix.PosDef.sqrt` and `Matrix.PosSemidef.sqrt` are gone on this Mathlib pin. The square
root is `CFC.sqrt`, which needs the scoped order instances of `MatrixOrder`. `KyFan.lean`
takes the root as an argument `S` with `Sᵀ = S`, `S * S = A` and `IsUnit S.det`; the three
lemmas below supply exactly those three facts for `S = abetaRSqrt β R`. -/

/-- `A_{β,R}^{1/2}`, the positive semidefinite square root from the continuous functional
calculus. -/
noncomputable def abetaRSqrt (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Matrix (Fin M) (Fin M) ℝ :=
  CFC.sqrt (AbetaR β R)

theorem posSemidef_abetaRSqrt (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (abetaRSqrt β R).PosSemidef :=
  Matrix.nonneg_iff_posSemidef.mp (CFC.sqrt_nonneg (AbetaR β R))

/-- `S := A_{β,R}^{1/2}` is symmetric: the hypothesis `hS` of `kyFan_min_congr`. -/
theorem transpose_abetaRSqrt (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (abetaRSqrt β R)ᵀ = abetaRSqrt β R := by
  have h := (posSemidef_abetaRSqrt β R).isHermitian
  rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- `S * S = A_{β,R}`: the hypothesis `hSA` of `kyFan_min_congr`. It needs `A_{β,R} ⪰ 0`,
which `abetaR_posDef` gives from `0 ≤ β_i < 1`. -/
theorem sq_abetaRSqrt {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) :
    abetaRSqrt β R * abetaRSqrt β R = AbetaR β R := by
  have hA : (0 : Matrix (Fin M) (Fin M) ℝ) ≤ AbetaR β R :=
    (abetaR_posDef R h0 h1).posSemidef.nonneg
  have h := CFC.sq_sqrt (AbetaR β R) hA
  rwa [pow_two] at h

/-- `IsUnit S.det`: the hypothesis `hSu` of `kyFan_min_congr`. `det S * det S = det A > 0`. -/
theorem isUnit_det_abetaRSqrt {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) : IsUnit (abetaRSqrt β R).det := by
  have hdet : 0 < (AbetaR β R).det := (abetaR_posDef R h0 h1).det_pos
  have hmul : (abetaRSqrt β R).det * (abetaRSqrt β R).det = (AbetaR β R).det := by
    rw [← Matrix.det_mul, sq_abetaRSqrt R h0 h1]
  refine isUnit_iff_ne_zero.mpr fun h => ?_
  rw [h, zero_mul] at hmul
  exact hdet.ne' hmul.symm

/-- `A_{β,R}^{-1/2} D A_{β,R}^{-1/2}`, the matrix whose bottom `r` eigenvalues give the paper's
`L⋆` (`main_paper.tex:897`). -/
noncomputable def DcongR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Matrix (Fin M) (Fin M) ℝ :=
  (abetaRSqrt β R)⁻¹ * Dmat β * (abetaRSqrt β R)⁻¹

theorem isHermitian_Dmat (β : Fin M → ℝ) : (Dmat β).IsHermitian :=
  Matrix.isHermitian_diagonal _

theorem isHermitian_DcongR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (DcongR β R).IsHermitian :=
  isHermitian_inv_mul_mul_inv (isHermitian_Dmat β) (transpose_abetaRSqrt β R)

/-- `L⋆ = r - ∑_{ℓ=1}^{r} λ_{r̃+1-ℓ}(A_{β,R}^{-1/2} D A_{β,R}^{-1/2})`
(`main_paper.tex:897`). The paper's index `r̃ + 1 - ℓ` at `ℓ = k + 1` is `Fin.rev` on
`Fin (Fintype.card (Fin M))`, the index form of `kyFan_min` in `LinAlg/KyFan.lean`. The
argument `hr` is a `Prop`, so the value does not depend on which proof is supplied. -/
noncomputable def limitROpt (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (hr : r ≤ Fintype.card (Fin M)) : ℝ :=
  (r : ℝ) - ∑ k : Fin r, (isHermitian_DcongR β R).eigenvalues₀ (Fin.castLE hr k).rev

/-! ### 5. The weighted estimator on an `UnalignedModel` -/

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- `Ṽ_W = W Ṽ`, the matrix whose top `r` right singular vectors are `V̂_svdstack(W)`
(`main_paper.tex:889`). -/
noncomputable def VtW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin M) (Fin (d N)) ℝ :=
  W * m.Vt N ω

/-- `Ṽ_W Ṽ_Wᵀ = W (Ṽ Ṽᵀ) Wᵀ`. -/
noncomputable def gramW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin M) (Fin M) ℝ :=
  m.VtW W N ω * (m.VtW W N ω)ᵀ

theorem isHermitian_gramW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (N : ℕ) (ω : Ω N) : (m.gramW W N ω).IsHermitian :=
  isHermitian_mul_transpose_self (m.VtW W N ω)

/-- `Ṽ_W V = W (Ṽ V)`. -/
noncomputable def VtVW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin M) (Fin r) ℝ :=
  m.VtW W N ω * m.V N

/-- Performance of weighted svdstack: `tr((Ṽ_W V)ᵀ (specInvTop (Ṽ_W Ṽ_Wᵀ) r) (Ṽ_W V))`, which
equals `‖V̂_svdstack(W)ᵀ V‖_F²` whenever `λ_r(Ṽ_W Ṽ_Wᵀ) > 0`. The paper's own first display
(`main_paper.tex:2016`) writes `V̂_svdstack(W) = Ṽᵀ X^{(d)}` with
`X^{(d)} = Wᵀ Q_r(W Ṽ Ṽᵀ Wᵀ) Λ_r^{-1/2}(W Ṽ Ṽᵀ Wᵀ)`, and
`X^{(d)} (X^{(d)})ᵀ = Wᵀ (specInvTop (W Ṽ Ṽᵀ Wᵀ) r) W`, so this is the same number. At `W = 1`
it is `perfR` (`perfRW_one`). -/
noncomputable def perfRW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVW W N ω)ᵀ * specInvTop (m.gramW W N ω) (m.isHermitian_gramW W N ω) r *
    m.VtVW W N ω)

theorem gramW_eq (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : m.gramW W N ω = W * m.gram N ω * Wᵀ := by
  rw [UnalignedModel.gramW, UnalignedModel.VtW, UnalignedModel.gram, Matrix.transpose_mul,
    Matrix.mul_assoc, ← Matrix.mul_assoc (m.Vt N ω), ← Matrix.mul_assoc]

theorem VtVW_eq (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : m.VtVW W N ω = W * m.VtV N ω := by
  rw [UnalignedModel.VtVW, UnalignedModel.VtW, UnalignedModel.VtV, Matrix.mul_assoc]

theorem gramW_one (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.gramW 1 N ω = m.gram N ω := by
  rw [m.gramW_eq, Matrix.one_mul, Matrix.transpose_one, Matrix.mul_one]

theorem VtVW_one (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.VtVW 1 N ω = m.VtV N ω := by
  rw [m.VtVW_eq, Matrix.one_mul]

theorem perfRW_one (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.perfRW 1 N ω = m.perfR N ω := by
  rw [UnalignedModel.perfRW, UnalignedModel.perfR,
    specInvTop_congr_mat (m.isHermitian_gramW 1 N ω) (m.isHermitian_gram N ω)
      (m.gramW_one N ω) r, m.VtVW_one N ω]

end UnalignedModel

/-! ### 6. The statements

`hgapW` is the paper's admissibility condition `λ_r(W A_{β,R} Wᵀ) > λ_{r+1}(W A_{β,R} Wᵀ)`
(`eq:weighted_eigengap`, `main_paper.tex:2010`, and the footnote at `main_paper.tex:890`
which excludes every other `W`). `hposW` says that the top `r` eigenvalues of `W A_{β,R} Wᵀ`
are positive, which the paper needs for `Λ_r^{-1/2}(W A_{β,R} Wᵀ)` to exist.

Cleanup wave 2 kept `hposW` on `limitRW_le_opt` and on
`thm_gen_rank_weight_svdstak_general`, against the audit's request to drop it. It excludes
only `rank W < r`, where the bound still holds numerically
(`notes/archive/audit_rank_r_weighted_post_2026-08-30.md`, finding 4.2, worst excess -8.0e-4 over
4000 draws), but both proofs read it, not only the statement: `limitRW_le_opt` builds the
frame `X = Wᵀ Q_r Λ_r^{-1/2}` through `exists_topFrame`, whose `Λ_r^{-1/2}` needs
`0 < λ_{r-1}`, and `thm_gen_rank_weight_svdstak_general` needs the same bound inside
`continuousAt_trace_specInvTop`. A route without `Λ_r^{-1/2}` is a separate task. -/

/-- `W⋆ A_{β,R} W⋆ᵀ = (W⋆ B_R)(W⋆ B_R)ᵀ + I` (`main_paper.tex:2103`), the identity that step 2
of the paper's proof starts from. -/
theorem abetaRW_optWR_eq {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h1 : ∀ i, β i ^ 2 < 1) :
    AbetaRW (optWR β) β R = BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1 := by
  have hD := optWR_mul_Dmat_mul_optWR (β := β) h1
  rw [Matrix.mul_assoc] at hD
  simp only [AbetaRW, AbetaR, BRW, Matrix.mul_add, Matrix.add_mul, Matrix.transpose_mul,
    Matrix.mul_assoc, hD]

/-- `W⋆ A_{β,R} W⋆ᵀ ⪰ I`, which gives the hypothesis `hposW` at `W⋆` for free. -/
theorem abetaRW_optWR_posDef {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) : (AbetaRW (optWR β) β R).PosDef := by
  rw [abetaRW_optWR_eq R (fun i => by nlinarith [h0 i, h1 i])]
  exact posSemidef_add_posDef
    (by simpa using Matrix.posSemidef_self_mul_conjTranspose (BRW (optWR β) β R))
    Matrix.PosDef.one

/-- **Step 1 of `thm:gen_rank_weight_svdstak`** (`main_paper.tex:2062`): no admissible weight
matrix beats `L⋆`. The paper's supremum over `{X : Xᵀ A_{β,R} X = I_r}` of `‖Xᵀ B_R‖_F²` is
`r - inf tr(Xᵀ D X)`, and the infimum is Ky Fan's minimum principle in the congruence form
`kyFan_min_congr`. -/
theorem limitRW_le_opt (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1)
    (hr : 0 < r) (hrM : r ≤ Fintype.card (Fin M))
    (hgapW : TopGap (AbetaRW W β R) (isHermitian_AbetaRW W β R) r)
    (hposW : 0 < (isHermitian_AbetaRW W β R).eigenvalues₀ ⟨r - 1, by omega⟩) :
    limitRW W β R ≤ limitROpt β R hrM := by
  classical
  obtain ⟨Q, L, hLpos, hQQ, hQAQ, hspec⟩ :=
    exists_topFrame (isHermitian_AbetaRW W β R) hr hrM hgapW hposW
  set X : Matrix (Fin M) (Fin r) ℝ :=
    Wᵀ * Q * Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) with hXdef
  have hXt : Xᵀ = Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) * Qᵀ * W := by
    rw [hXdef, Matrix.transpose_mul, Matrix.transpose_mul, Matrix.diagonal_transpose,
      Matrix.transpose_transpose]
    simp only [Matrix.mul_assoc]
  have hsqrtne : ∀ k, Real.sqrt (L k) ≠ 0 := fun k => (Real.sqrt_pos.mpr (hLpos k)).ne'
  have hsq : ∀ k, (Real.sqrt (L k))⁻¹ * L k * (Real.sqrt (L k))⁻¹ = 1 := by
    intro k
    have hs : Real.sqrt (L k) * Real.sqrt (L k) = L k := Real.mul_self_sqrt (hLpos k).le
    calc (Real.sqrt (L k))⁻¹ * L k * (Real.sqrt (L k))⁻¹
        = L k * ((Real.sqrt (L k))⁻¹ * (Real.sqrt (L k))⁻¹) := by ring
      _ = L k * (Real.sqrt (L k) * Real.sqrt (L k))⁻¹ := by rw [mul_inv]
      _ = L k * (L k)⁻¹ := by rw [hs]
      _ = 1 := mul_inv_cancel₀ (hLpos k).ne'
  have hss : Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) *
      Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹)
      = Matrix.diagonal (fun k => (L k)⁻¹) := by
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext k
    have hs : Real.sqrt (L k) * Real.sqrt (L k) = L k := Real.mul_self_sqrt (hLpos k).le
    rw [← mul_inv, hs]
  have hXA : Xᵀ * AbetaR β R * X = 1 := by
    have hkey : Xᵀ * AbetaR β R * X
        = Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) *
          (Qᵀ * AbetaRW W β R * Q) * Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) := by
      rw [hXt, hXdef, AbetaRW]
      simp only [Matrix.mul_assoc]
    rw [hkey, hQAQ, Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal,
      ← Matrix.diagonal_one]
    congr 1
    funext k
    simpa using hsq k
  have hXXt : X * Xᵀ
      = Wᵀ * specInvTop (AbetaRW W β R) (isHermitian_AbetaRW W β R) r * W := by
    calc X * Xᵀ
        = Wᵀ * Q * (Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹) *
            Matrix.diagonal (fun k => (Real.sqrt (L k))⁻¹)) * Qᵀ * W := by
          rw [hXdef, hXt]; simp only [Matrix.mul_assoc]
      _ = Wᵀ * (Q * Matrix.diagonal (fun k => (L k)⁻¹) * Qᵀ) * W := by
          rw [hss]; simp only [Matrix.mul_assoc]
      _ = _ := by rw [hspec]
  have hlim : limitRW W β R = Matrix.trace (Xᵀ * (BR β R * (BR β R)ᵀ) * X) := by
    have hA1 : (BRW W β R)ᵀ * specInvTop (AbetaRW W β R) (isHermitian_AbetaRW W β R) r *
          BRW W β R
        = (BR β R)ᵀ * (X * Xᵀ) * BR β R := by
      rw [hXXt, BRW, Matrix.transpose_mul]
      simp only [Matrix.mul_assoc]
    have hA2 : (BR β R)ᵀ * (X * Xᵀ) * BR β R = ((BR β R)ᵀ * X) * (Xᵀ * BR β R) := by
      simp only [Matrix.mul_assoc]
    have hA3 : Xᵀ * (BR β R * (BR β R)ᵀ) * X = (Xᵀ * BR β R) * ((BR β R)ᵀ * X) := by
      simp only [Matrix.mul_assoc]
    rw [limitRW, hA1, hA2, hA3]
    exact Matrix.trace_mul_comm _ _
  have hsplit : Matrix.trace (Xᵀ * (BR β R * (BR β R)ᵀ) * X)
      = (r : ℝ) - Matrix.trace (Xᵀ * Dmat β * X) := by
    have hBB : BR β R * (BR β R)ᵀ = AbetaR β R - Dmat β := by
      rw [AbetaR, add_sub_cancel_right]
    rw [hBB, Matrix.mul_sub, Matrix.sub_mul, Matrix.trace_sub, hXA, Matrix.trace_one]
    simp
  have hkf : ∑ k : Fin r, (isHermitian_DcongR β R).eigenvalues₀ (Fin.castLE hrM k).rev
      ≤ Matrix.trace (Xᵀ * Dmat β * X) :=
    (kyFan_min_congr (transpose_abetaRSqrt β R) (sq_abetaRSqrt R h0 h1)
      (isUnit_det_abetaRSqrt R h0 h1) (isHermitian_DcongR β R) hrM).2 ⟨X, hXA, rfl⟩
  rw [limitROpt, hlim, hsplit]
  linarith

/-- **Step 2 of `thm:gen_rank_weight_svdstak`** (`main_paper.tex:2101`): `W⋆ = D^{-1/2}`
attains `L⋆`. The eigengap of `W⋆ A_{β,R} W⋆ᵀ` is **not** needed: `W⋆ A_{β,R} W⋆ᵀ` has the
shape `C Cᵀ + I` with `C` of `r` columns, so a boundary tie sits at the eigenvalue `1` and
contributes nothing (`trace_specInvTop_add_one`). This drops the hypothesis `hgap` of the
first statement of this file (`notes/archive/audit_rank_r_weighted_2026-08-30.md`, finding 5.2;
numeric row H5 of `notes/archive/agent_reports/rank_r_t5b_proofs.md`). -/
theorem limitRW_optWR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) (hrM : r ≤ Fintype.card (Fin M)) :
    limitRW (optWR β) β R = limitROpt β R hrM := by
  classical
  have hsq : ∀ i, β i ^ 2 < 1 := fun i => by nlinarith [h0 i, h1 i]
  have heq : AbetaRW (optWR β) β R
      = BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1 := abetaRW_optWR_eq R hsq
  have hA' : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1).IsHermitian :=
    (isHermitian_mul_transpose_self _).add Matrix.isHermitian_one
  have hval : limitRW (optWR β) β R
      = (r : ℝ) - ∑ k : Fin r, (hA'.eigenvalues₀ (Fin.castLE hrM k))⁻¹ := by
    rw [limitRW, specInvTop_congr_mat (isHermitian_AbetaRW (optWR β) β R) hA' heq r]
    exact trace_specInvTop_add_one _ hrM hA'
  have hHW : (AbetaRW (optWR β) β R).IsHermitian := isHermitian_AbetaRW (optWR β) β R
  have hAWpd : (AbetaRW (optWR β) β R).PosDef := abetaRW_optWR_posDef R h0 h1
  have hAWpos : ∀ k, 0 < hHW.eigenvalues₀ k := by
    intro k
    rw [← eigenvalues_eigIdx hHW]
    exact hAWpd.eigenvalues_pos _
  have hWt : (AbetaRW (optWR β) β R)ᵀ = AbetaRW (optWR β) β R := by
    have h := hHW
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have hinvH : ((AbetaRW (optWR β) β R)⁻¹).IsHermitian := by
    rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial,
      Matrix.transpose_nonsing_inv, hWt]
  have hAdet : IsUnit (AbetaR β R).det :=
    isUnit_iff_ne_zero.mpr (abetaR_posDef R h0 h1).det_pos.ne'
  have hDsq : Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) *
      Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) = Dmat β := by
    rw [Matrix.diagonal_mul_diagonal, Dmat]
    congr 1
    funext i
    exact Real.mul_self_sqrt (by linarith [hsq i])
  have hWoD : optWR β * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) = 1 := by
    rw [optWR, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    funext i
    have hp : (0 : ℝ) < 1 - β i ^ 2 := by linarith [hsq i]
    have hs : Real.sqrt (1 - β i ^ 2) ≠ 0 := (Real.sqrt_pos.mpr hp).ne'
    rw [optW]
    field_simp
  have hSS : (abetaRSqrt β R)⁻¹ * (abetaRSqrt β R)⁻¹ = (AbetaR β R)⁻¹ := by
    rw [← sq_abetaRSqrt R h0 h1, Matrix.mul_inv_rev]
  have hAWinv : (AbetaRW (optWR β) β R)⁻¹
      = Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (AbetaR β R)⁻¹ *
        Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) := by
    refine Matrix.inv_eq_right_inv ?_
    rw [AbetaRW, transpose_optWR]
    calc optWR β * AbetaR β R * optWR β *
          (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (AbetaR β R)⁻¹ *
            Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)))
        = optWR β * AbetaR β R *
            (optWR β * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2))) *
            ((AbetaR β R)⁻¹ * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2))) := by
          simp only [Matrix.mul_assoc]
      _ = 1 := by
          rw [hWoD, Matrix.mul_one, Matrix.mul_assoc, ← Matrix.mul_assoc (AbetaR β R),
            Matrix.mul_nonsing_inv _ hAdet, Matrix.one_mul, hWoD]
  have hGt : (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹)ᵀ
      = (abetaRSqrt β R)⁻¹ * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) := by
    rw [Matrix.transpose_mul, Matrix.diagonal_transpose, Matrix.transpose_nonsing_inv,
      transpose_abetaRSqrt]
  have hGtG : (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹)ᵀ *
      (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹)
      = DcongR β R := by
    rw [hGt, DcongR]
    calc (abetaRSqrt β R)⁻¹ * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) *
          (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹)
        = (abetaRSqrt β R)⁻¹ * (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) *
            Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2))) * (abetaRSqrt β R)⁻¹ := by
          simp only [Matrix.mul_assoc]
      _ = (abetaRSqrt β R)⁻¹ * Dmat β * (abetaRSqrt β R)⁻¹ := by rw [hDsq]
  have hGGt : (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹) *
      (Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹)ᵀ
      = (AbetaRW (optWR β) β R)⁻¹ := by
    rw [hGt, hAWinv]
    calc Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) * (abetaRSqrt β R)⁻¹ *
          ((abetaRSqrt β R)⁻¹ * Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)))
        = Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) *
            ((abetaRSqrt β R)⁻¹ * (abetaRSqrt β R)⁻¹) *
            Matrix.diagonal (fun i => Real.sqrt (1 - β i ^ 2)) := by
          simp only [Matrix.mul_assoc]
      _ = _ := by rw [hSS]
  have hchar : (DcongR β R).charpoly = ((AbetaRW (optWR β) β R)⁻¹).charpoly := by
    rw [← hGtG, ← hGGt]
    exact (Matrix.charpoly_mul_comm _ _).symm
  have hDeq := eigenvalues₀_eq_of_charpoly_eq (isHermitian_DcongR β R) hinvH hchar
  have hkey : ∀ j, (isHermitian_DcongR β R).eigenvalues₀ j = (hA'.eigenvalues₀ j.rev)⁻¹ := by
    intro j
    rw [congrFun hDeq j, eigenvalues₀_inv hHW hAWpos hinvH j,
      congrFun (eigenvalues₀_congr_mat hHW hA' heq) j.rev]
  rw [hval, limitROpt]
  congr 1
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [hkey, Fin.rev_rev]

/-! ### 6b. The eigengap at `W⋆`

The paper's check (`main_paper.tex:2103`): `W⋆ A_{β,R} W⋆ᵀ = (W⋆ B_R)(W⋆ B_R)ᵀ + I`, and the
positive semidefinite `(W⋆ B_R)(W⋆ B_R)ᵀ` has exactly `r` nonzero eigenvalues when `B_R` has
rank `r`. `rank B_R = r` is what the proof needs; the paper's `β_ij > 0` plus
`Rank(∑ R_i R_iᵀ) = r` is one sufficient condition for it (`rank_BR_of_ne_zero`), and it is
strictly stronger: at `M = 3`, `r = 2`, `R = (e₁, e₂, e₁)` and `β = (0.5, 0, 0.5)` the paper's
rank condition holds, `β ≠ 0` fails for one index, and the gap is exactly `0` (numeric row H3,
`notes/archive/agent_reports/rank_r_t5b_proofs.md`). -/

/-- `rank B_R = r` from `β_i ≠ 0` and the paper's rank condition on `R_stack`. -/
theorem rank_BR_of_ne_zero {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (hβ : ∀ i, β i ≠ 0) (hrank : (Rstack R).rank = r) : (BR β R).rank = r := by
  have hdet : IsUnit (Matrix.diagonal β).det := by
    rw [Matrix.det_diagonal]
    exact isUnit_iff_ne_zero.mpr (Finset.prod_ne_zero_iff.mpr fun i _ => hβ i)
  rw [BR_eq_diagonal_mul_Rstack, Matrix.rank_mul_eq_right_of_isUnit_det _ _ hdet]
  exact hrank

/-- `W⋆` is invertible, so it does not change the rank of `B_R`. -/
theorem rank_BRW_optWR {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h1 : ∀ i, β i ^ 2 < 1) : (BRW (optWR β) β R).rank = (BR β R).rank := by
  have hne : ∀ i, optW β i ≠ 0 := by
    intro i
    have hp : (0 : ℝ) < 1 - β i ^ 2 := by linarith [h1 i]
    have : 0 < optW β i := by
      rw [optW]
      positivity
    exact this.ne'
  have hdet : IsUnit (optWR β).det := by
    rw [optWR, Matrix.det_diagonal]
    exact isUnit_iff_ne_zero.mpr (Finset.prod_ne_zero_iff.mpr fun i _ => hne i)
  rw [BRW, Matrix.rank_mul_eq_right_of_isUnit_det _ _ hdet]

/-- **The eigengap of step 2** in the form the proof needs: `rank B_R = r`. -/
theorem topGap_optWR_of_rankBR {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) (hrankB : (BR β R).rank = r) :
    TopGap (AbetaRW (optWR β) β R) (isHermitian_AbetaRW (optWR β) β R) r := by
  have hsq : ∀ i, β i ^ 2 < 1 := fun i => by nlinarith [h0 i, h1 i]
  have hPrank : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).rank = r := by
    rw [Matrix.rank_self_mul_transpose, rank_BRW_optWR R hsq, hrankB]
  have hP : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).IsHermitian :=
    isHermitian_mul_transpose_self _
  have hPSD : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose (BRW (optWR β) β R)
  have hnn : ∀ k, 0 ≤ hP.eigenvalues₀ k := by
    intro k
    rw [← eigenvalues_eigIdx hP]
    exact hPSD.eigenvalues_nonneg _
  have hP1 : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1).IsHermitian :=
    hP.add (Matrix.isHermitian_one)
  have heq : AbetaRW (optWR β) β R
      = BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1 := abetaRW_optWR_eq R hsq
  have hev : ∀ j, (isHermitian_AbetaRW (optWR β) β R).eigenvalues₀ j = hP.eigenvalues₀ j + 1 := by
    intro j
    rw [eigenvalues₀_congr_mat (isHermitian_AbetaRW (optWR β) β R) hP1 heq]
    exact eigenvalues₀_add_one hP hP1 j
  intro k l hk hl
  rw [hev k, hev l]
  have hzero : hP.eigenvalues₀ l = 0 :=
    eigenvalues₀_eq_zero_of_rank_le hP hnn (le_of_eq hPrank) hl
  have hposk : 0 < hP.eigenvalues₀ k :=
    eigenvalues₀_pos_of_lt_rank hP hnn (le_of_eq hPrank.symm) hk
  linarith

/-- The eigengap check of step 2 (`main_paper.tex:2103`) in the paper's own hypotheses:
`diag(β)` is invertible (`hβpos`, the paper's `β_ij > 0`) and `R_stack` has rank `r` (`hrank`,
the rank condition of `assum:unaligned`). -/
theorem topGap_optWR_of_rank (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (hβpos : ∀ i, 0 < β i) (h1 : ∀ i, β i < 1) (hrank : (Rstack R).rank = r) :
    TopGap (AbetaRW (optWR β) β R) (isHermitian_AbetaRW (optWR β) β R) r :=
  topGap_optWR_of_rankBR R (fun i => (hβpos i).le) h1
    (rank_BR_of_ne_zero R (fun i => (hβpos i).ne') hrank)

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- Two entrywise helpers: `W G Wᵀ` and `W Y` are fixed real linear combinations of the
entries of `G` and of `Y`, which is what makes the weighted limits a continuous image of the
unweighted ones. -/
theorem mul_mul_transpose_apply (W G : Matrix (Fin M) (Fin M) ℝ) (i j : Fin M) :
    (W * G * Wᵀ) i j = ∑ b : Fin M, ∑ a : Fin M, W i a * G a b * W j b := by
  rw [Matrix.mul_apply]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [Matrix.mul_apply, Matrix.transpose_apply, Finset.sum_mul]

end UnalignedModel

end StackedSVD
