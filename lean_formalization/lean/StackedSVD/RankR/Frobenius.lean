/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SubspaceMain
import StackedSVD.RankR.WeightedMain

/-!
# The paper's Frobenius performance `‖V̂ᵀ V‖_F²` at rank `r`

STATUS 2026-08-31: see `notes/archive/agent_reports/polish_frobenius_secular.md` and
`notes/archive/agent_reports/fix_global_audit.md` (finding 1 of the global audit: the frame
hypothesis of the two corollaries below is now quantified with probability tending to one, not over
every `ω`).

Section 7 of the paper measures both methods by `‖V̂ᵀ V‖_F²` with `V̂` a matrix of `r`
orthonormal columns. `RankR/Defs.lean` and `RankR/Subspace.lean` state the two propositions
with `perfR` (a trace against `specInvTop`) and `perfStackR` (a sum of squared projector
norms), and the identification with the paper's display lives only in a doc-comment. This
file removes that gap: it defines the estimator objects and proves the two identities, so
that each proposition gets a corollary stated exactly on `‖V̂ᵀ V‖_F²`.

## Content

1. `frobSq X = ∑_{i,j} X_{ij}²`, with `frobSq_eq_trace`. Mathlib's Frobenius norm is a scoped
   instance on `Matrix` that clashes with `Matrix.Norms.L2Operator`, which `LinAlg/` opens, so
   the square of the norm is written out.
2. `IsTopFrame A hA Y`: the columns of `Y` are orthonormal and they span `specTop A hA r`.
   This is "any orthonormal basis of the top-`r` eigenspace". No choice of eigenvector
   appears, and `frobenius_eq_perf` shows the performance does not see the choice. A frame
   exists exactly when `dim (specTop A hA r) = r`, that is when no eigenvalue of `A` ties
   across the index `r` boundary.
3. `frobenius_eq_perf`: `‖Yᵀ V‖_F² = ∑_k ‖P_top(A, r) (V e_k)‖²` for every top frame `Y`.
   With `A` the stacksvd Gram this is exactly `perfStackR`, so
   `prop_stacksvd_subspace_frobenius` is a corollary of `prop_stacksvd_subspace`.
4. `IsTopEigFrame A hA Y lam`: a top frame together with its eigenvalues, `A Y = Y diag(λ)`.
   This is the paper's pair `(Q_r, Λ_r)`. `toOp_specInvTop_apply` gives
   `specInvTop A hA r x = ∑_j λ_j⁻¹ ⟪y_j, x⟫ y_j`, and `frobSq_eq_trace_of_eigFrame` turns the
   trace form into a Frobenius norm.
5. `topEigMat` and `topEigVal`, the canonical choice of `(Q_r, Λ_r)`: the first `r` vectors of
   the sorted eigenbasis, in the index convention of `eigenvalues₀` and of `specProjTop`.
   `isTopEigFrame_topEigMat` proves it is a top-`r` eigenframe under `TopGap A hA r`, which is
   the existence statement that the corollaries need.
6. `UnalignedModel.vhatSvdstack` and `vhatSvdstackW`: the paper's `V̂_svdstack = Ṽᵀ Q Λ^{-1/2}`
   (`main_paper.tex:1982` and `2016`). `frobSq_vhatSvdstack` and `frobSq_vhatSvdstackW` prove
   `‖V̂ᵀ V‖_F² = perfR` and `= perfRW W`, so the two svdstack propositions get corollaries on
   the paper's own quantity.
7. `topGap_gram_whp`: the top-`r` gap of `Ṽ Ṽᵀ` with probability tending to one, from `gramR`
   and `specTop_simple_whp_of_tendsto`. With item 5 it discharges the frame hypothesis of the
   unweighted svdstack corollary (`..._frobenius_eig`).

## What is proved for svdstack, exactly

`perfR` is `tr((Ṽ V)ᵀ (specInvTop (Ṽ Ṽᵀ) r) (Ṽ V))`, a trace on the `M × M` side. The paper's
`V̂_svdstack` is a `d × r` matrix, and the paper's own first display writes it as
`V̂_svdstack = Ṽᵀ Q_r Λ_r^{-1/2}` with `(Q_r, Λ_r)` the top-`r` eigenpairs of `Ṽ Ṽᵀ`. This file
follows that display: `vhatSvdstack` is that matrix, and `frobSq_vhatSvdstack` proves
`‖V̂_svdstackᵀ V‖_F² = perfR` for **any** top-`r` eigenframe `(Q, Λ)` of `Ṽ Ṽᵀ` with
`λ_j ≥ 0`. Two further lemmas justify the name: `vhatSvdstack_transpose_mul_self` shows the
columns are orthonormal (when `λ_j > 0`), and `gram_mulVec_vhatSvdstack` shows each column is
an eigenvector of `Ṽᵀ Ṽ` with eigenvalue `λ_j`. So `V̂_svdstack` is a matrix of `r`
orthonormal right singular vectors of `Ṽ` for the eigenvalues `λ_1, …, λ_r`.

The one statement this file does **not** prove is that these `r` right singular vectors are a
frame of `specTop (Ṽᵀ Ṽ) r`, that is that the top-`r` index sets of `Ṽ Ṽᵀ` and of `Ṽᵀ Ṽ`
select the same eigenvalues. That needs the sorted eigenvalue lists of `B Bᵀ` and `Bᵀ B` to
agree on the nonzero part, which `topEigSet` states through `eigenvalues₀` indices and which
no lemma in the project supplies. It is not needed: the paper's display defines `V̂_svdstack`
from `(Q_r, Λ_r)`, and that is what is used here.

## Where the frame hypothesis holds

A top-`r` frame of `A` exists exactly when `dim (specTop A hA r) = r`, that is when no
eigenvalue of `A` ties across the index `r` boundary. A Gaussian model has such an `ω` at
every `N` when `r < M` (svdstack: all `v̂_i` orthogonal gives `Ṽ Ṽᵀ = I_M`) or when
`r < d N` (stacksvd: `E = -θ u vᵀ` gives `X_stack = 0`). So a frame hypothesis quantified over
every `ω` is false there, and the corollary that carries it is vacuous (global audit
2026-08-31, finding 1). The two corollaries of section 4c ask the frame only on an event whose
probability tends to one and transport the conclusion off the bad event
(`TendstoInProb.of_tendsto_measure_ne_of_tendsto`), as `SVDStack/Main.lean:301` does for the
rank-one tie. Their `_eig` companions fix the canonical frame of item 5 and carry no frame
hypothesis: for unweighted svdstack the hypotheses are exactly those of
`prop_general_rank_unweighted_svdstack` plus `0 < r`; for stacksvd the top-`r` gap of
`X_stackᵀ X_stack` with probability tending to one stays an explicit hypothesis, because
`SubspaceLaw` has no gap field (decision D18). `thm_gen_rank_weight_svdstak_frobenius` keeps
the `∀ N ω` frame; it is used at `r = M`, where `TopGap` holds for every matrix
(`topGap_of_card_le`) and item 5 gives the frame at every `ω`.

## Measurability

`TendstoInProb μ f a` is `∀ ε > 0, Tendsto (fun N => μ N {ω | ε ≤ |f N ω - a|}) atTop (𝓝 0)`,
and `μ N` is an outer measure, defined on every subset. So the selections `Q N ω`, `lam N ω`
and `Y N ω` below need no measurability: the corollaries transport a limit in probability
along an event, or along a pointwise equality of functions (`TendstoInProb.congr`), and
neither reads a measurable set. The same remark applies to `perfR` itself, which is why the
propositions of `RankR/Defs.lean` carry no measurability hypothesis either.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. The squared Frobenius norm -/

section Frob

variable {p q r : ℕ}

/-- Column `j` of a matrix, as a vector of `EuclideanSpace`. -/
noncomputable def frameCol (Y : Matrix (Fin p) (Fin r) ℝ) (j : Fin r) :
    EuclideanSpace ℝ (Fin p) := WithLp.toLp 2 fun i => Y i j

@[simp]
theorem frameCol_apply (Y : Matrix (Fin p) (Fin r) ℝ) (j : Fin r) (i : Fin p) :
    frameCol Y j i = Y i j := rfl

/-- The squared Frobenius norm `∑_{i,j} X_{ij}²`. -/
noncomputable def frobSq (X : Matrix (Fin p) (Fin q) ℝ) : ℝ := ∑ i, ∑ j, X i j ^ 2

theorem frobSq_eq_trace (X : Matrix (Fin p) (Fin q) ℝ) :
    frobSq X = Matrix.trace (Xᵀ * X) := by
  rw [Matrix.trace, frobSq, Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Matrix.diag_apply, Matrix.mul_apply]
  exact Finset.sum_congr rfl fun i _ => by rw [Matrix.transpose_apply]; exact sq _

/-- `⟪y_j, x_k⟫` is the `(j, k)` entry of `Yᵀ X`. -/
theorem inner_frameCol (Y : Matrix (Fin p) (Fin r) ℝ) (X : Matrix (Fin p) (Fin q) ℝ)
    (j : Fin r) (k : Fin q) : ⟪frameCol Y j, frameCol X k⟫_ℝ = (Yᵀ * X) j k := by
  rw [real_inner_eq_dotProduct, Matrix.mul_apply, dotProduct]
  exact Finset.sum_congr rfl fun i _ => rfl

/-- Orthonormal columns, read on the inner products. -/
theorem inner_frameCol_of_ortho {Y : Matrix (Fin p) (Fin r) ℝ} (h : Yᵀ * Y = 1) (i j : Fin r) :
    ⟪frameCol Y i, frameCol Y j⟫_ℝ = if i = j then 1 else 0 := by
  rw [inner_frameCol, h, Matrix.one_apply]

end Frob

/-! ### 2. Top-`r` frames -/

section Frame

variable {p q r : ℕ}

/-- `Y` is a top-`r` frame of `A`: its `r` columns are orthonormal and they span the top-`r`
eigenspace. This is the paper's "matrix of the top `r` eigenvectors", written without a choice
of eigenvector: two frames differ by an `r × r` orthogonal matrix, and `frobenius_eq_perf`
shows the performance does not see the difference. -/
structure IsTopFrame (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (Y : Matrix (Fin p) (Fin r) ℝ) : Prop where
  /-- the columns are orthonormal -/
  ortho : Yᵀ * Y = 1
  /-- the columns span the top-`r` eigenspace -/
  span_eq : Submodule.span ℝ (Set.range fun j => frameCol Y j) = specTop A hA r

variable {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {Y : Matrix (Fin p) (Fin r) ℝ}

/-- A frame column lies in the top-`r` eigenspace. -/
theorem IsTopFrame.mem (hY : IsTopFrame A hA Y) (j : Fin r) :
    frameCol Y j ∈ specTop A hA r := by
  rw [← hY.span_eq]
  exact Submodule.subset_span ⟨j, rfl⟩

/-- The projector onto the top-`r` eigenspace, expanded in a frame. -/
theorem specProjTop_eq_frame_sum (hY : IsTopFrame A hA Y) (x : EuclideanSpace ℝ (Fin p)) :
    specProjTop A hA r x = ∑ j, ⟪frameCol Y j, x⟫_ℝ • frameCol Y j := by
  have hproj : specProjTop A hA r x = (specTop A hA r).starProjection x := rfl
  rw [hproj]
  refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero
    (Submodule.sum_mem _ fun j _ => Submodule.smul_mem _ _ (hY.mem j)) fun w hw => ?_
  rw [← hY.span_eq] at hw
  refine Submodule.span_induction ?_ ?_ ?_ ?_ hw
  · rintro _ ⟨k, rfl⟩
    rw [inner_sub_left, sum_inner]
    have hterm : ∀ j : Fin r,
        ⟪⟪frameCol Y j, x⟫_ℝ • frameCol Y j, frameCol Y k⟫_ℝ
          = if j = k then ⟪frameCol Y k, x⟫_ℝ else 0 := by
      intro j
      rw [real_inner_smul_left, inner_frameCol_of_ortho hY.ortho]
      by_cases h : j = k
      · subst h; simp
      · simp [h]
    rw [Finset.sum_congr rfl fun j (_ : j ∈ Finset.univ) => hterm j,
      Finset.sum_ite_eq' Finset.univ k fun _ => ⟪frameCol Y k, x⟫_ℝ]
    simp [real_inner_comm x (frameCol Y k)]
  · simp
  · intro u v _ _ hu hv
    rw [inner_add_right, hu, hv, add_zero]
  · intro a u _ hu
    rw [real_inner_smul_right, hu, mul_zero]

/-- The squared norm of the projection, in a frame. -/
theorem norm_sq_specProjTop_frame (hY : IsTopFrame A hA Y) (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProjTop A hA r x‖ ^ 2 = ∑ j, ⟪frameCol Y j, x⟫_ℝ ^ 2 := by
  rw [← real_inner_self_eq_norm_sq, specProjTop_eq_frame_sum hY, sum_inner]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [real_inner_smul_left, inner_sum]
  have hterm : ∀ i : Fin r, ⟪frameCol Y j, ⟪frameCol Y i, x⟫_ℝ • frameCol Y i⟫_ℝ
      = if i = j then ⟪frameCol Y j, x⟫_ℝ else 0 := by
    intro i
    rw [real_inner_smul_right, inner_frameCol_of_ortho hY.ortho]
    by_cases h : i = j
    · subst h; simp
    · simp [h, Ne.symm h]
  rw [Finset.sum_congr rfl fun i (_ : i ∈ Finset.univ) => hterm i,
    Finset.sum_ite_eq' Finset.univ j fun _ => ⟪frameCol Y j, x⟫_ℝ]
  simp [sq]

/-- **The Frobenius identity.** For every top-`r` frame `Y` of `A`,
`‖Yᵀ V‖_F² = ∑_k ‖P_top(A, r) (V e_k)‖²`. The right hand side reads no frame, so the left
hand side does not depend on the choice of frame. -/
theorem frobenius_eq_perf (hY : IsTopFrame A hA Y) (V : Matrix (Fin p) (Fin q) ℝ) :
    frobSq (Yᵀ * V) = ∑ k, ‖specProjTop A hA r (frameCol V k)‖ ^ 2 := by
  rw [frobSq, Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [norm_sq_specProjTop_frame hY]
  exact Finset.sum_congr rfl fun j _ => by rw [inner_frameCol]

end Frame

/-! ### 3. Top-`r` eigenframes and `specInvTop` -/

section EigFrame

variable {p q r : ℕ}

/-- The matrix eigenbasis diagonalizes the operator. `LinAlg/SpecProjPerturb.lean` has the
same lemma, private. -/
private theorem apply_eigvecF {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (i : Fin p) :
    toOp A (hA.eigenvectorBasis i) = hA.eigenvalues i • hA.eigenvectorBasis i := by
  apply WithLp.ofLp_injective
  simpa using hA.mulVec_eigenvectorBasis i

/-- `toEuclideanCLM` is `toOp`, as a function. -/
private theorem toCLM_eq_toOp (B : Matrix (Fin p) (Fin p) ℝ) (x : EuclideanSpace ℝ (Fin p)) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) B x = toOp B x := rfl

/-- `specInvTop` is symmetric: it is a real linear combination of the rank-one matrices
`b_i b_iᵀ`. -/
theorem isHermitian_specInvTop {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (r : ℕ) :
    (specInvTop A hA r).IsHermitian := by
  have hsym : ∀ i j, specInvTop A hA r i j = specInvTop A hA r j i := by
    intro i j
    simp only [specInvTop, Matrix.sum_apply, Matrix.smul_apply, Matrix.vecMulVec_apply,
      smul_eq_mul]
    exact Finset.sum_congr rfl fun k _ => by ring
  ext i j
  simpa [Matrix.conjTranspose_apply] using (hsym j i)

/-- A top-`r` frame together with the eigenvalues of its columns, `A Y = Y diag(λ)`. This is
the paper's pair `(Q_r, Λ_r)` (`main_paper.tex:1982`). -/
structure IsTopEigFrame (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (Y : Matrix (Fin p) (Fin r) ℝ) (lam : Fin r → ℝ) : Prop where
  /-- the columns are an orthonormal basis of the top-`r` eigenspace -/
  frame : IsTopFrame A hA Y
  /-- column `j` is an eigenvector with eigenvalue `lam j` -/
  eig : A * Y = Y * Matrix.diagonal lam

variable {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {Y : Matrix (Fin p) (Fin r) ℝ}
  {lam : Fin r → ℝ}

set_option linter.deprecated false in
/-- The eigenvector equation, column by column. -/
theorem toOp_frameCol (h : A * Y = Y * Matrix.diagonal lam) (j : Fin r) :
    toOp A (frameCol Y j) = lam j • frameCol Y j := by
  apply WithLp.ofLp_injective
  have hmv : A *ᵥ WithLp.ofLp (frameCol Y j) = lam j • WithLp.ofLp (frameCol Y j) := by
    funext i
    have hij := congrFun (congrFun h i) j
    have hL : (A * Y) i j = (A *ᵥ WithLp.ofLp (frameCol Y j)) i := by
      rw [Matrix.mul_apply, Matrix.mulVec, dotProduct]
      exact Finset.sum_congr rfl fun k _ => rfl
    have hR : (Y * Matrix.diagonal lam) i j = lam j * Y i j := by
      rw [Matrix.mul_apply]
      simp [Matrix.diagonal_apply, mul_comm]
    rw [hL, hR] at hij
    simpa using hij
  simpa [Matrix.toEuclideanLin_apply] using hmv

/-- Every frame eigenvalue is a top-`r` eigenvalue of `A`. Otherwise its eigenvector would be
orthogonal to the whole top-`r` eigenspace, which contains it, so it would be zero. -/
theorem IsTopEigFrame.lam_mem (hY : IsTopEigFrame A hA Y lam) (j : Fin r) :
    lam j ∈ topEigSet A hA r := by
  by_contra hmem
  have hle : specSpace A (topEigSet A hA r) ≤ (ℝ ∙ frameCol Y j)ᗮ := by
    refine iSup_le fun t => iSup_le fun ht z hz => ?_
    rw [Module.End.mem_eigenspace_iff] at hz
    refine Submodule.mem_orthogonal_singleton_iff_inner_right.mpr ?_
    exact inner_eq_zero_of_ne hA (fun h => hmem (by rw [h]; exact ht))
      (toOp_frameCol hY.eig j) hz
  have hzero : ⟪frameCol Y j, frameCol Y j⟫_ℝ = 0 :=
    Submodule.mem_orthogonal_singleton_iff_inner_right.mp (hle (hY.frame.mem j))
  rw [inner_frameCol_of_ortho hY.frame.ortho, if_pos rfl] at hzero
  exact one_ne_zero hzero

/-- `specInvTop` inverts `A` on a frame column. -/
theorem toOp_specInvTop_frameCol (hY : IsTopEigFrame A hA Y lam) (j : Fin r) :
    toOp (specInvTop A hA r) (frameCol Y j) = (lam j)⁻¹ • frameCol Y j := by
  have hmem := hY.lam_mem j
  have hexp := (hA.eigenvectorBasis).sum_repr' (frameCol Y j)
  calc toOp (specInvTop A hA r) (frameCol Y j)
      = toOp (specInvTop A hA r)
          (∑ i, ⟪hA.eigenvectorBasis i, frameCol Y j⟫_ℝ • hA.eigenvectorBasis i) := by
        rw [hexp]
    _ = ∑ i, ⟪hA.eigenvectorBasis i, frameCol Y j⟫_ℝ •
          (invCoef A hA r i • hA.eigenvectorBasis i) := by
        rw [map_sum]
        exact Finset.sum_congr rfl fun i _ => by
          rw [map_smul, toOp_specInvTop_eigvec]
    _ = ∑ i, ((lam j)⁻¹ * ⟪hA.eigenvectorBasis i, frameCol Y j⟫_ℝ) • hA.eigenvectorBasis i := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [smul_smul]
        by_cases hev : hA.eigenvalues i = lam j
        · rw [invCoef_of_mem (by rw [hev]; exact hmem), hev, mul_comm]
        · have h0 : ⟪hA.eigenvectorBasis i, frameCol Y j⟫_ℝ = 0 :=
            inner_eq_zero_of_ne hA hev (apply_eigvecF hA i) (toOp_frameCol hY.eig j)
          simp [h0]
    _ = (lam j)⁻¹ • ∑ i, ⟪hA.eigenvectorBasis i, frameCol Y j⟫_ℝ • hA.eigenvectorBasis i := by
        rw [Finset.smul_sum]
        exact Finset.sum_congr rfl fun i _ => mul_smul _ _ _
    _ = (lam j)⁻¹ • frameCol Y j := by rw [hexp]

/-- **`specInvTop` in a top-`r` eigenframe**: `specInvTop A r x = ∑_j λ_j⁻¹ ⟪y_j, x⟫ y_j`,
that is `specInvTop A r = Y diag(λ⁻¹) Yᵀ`. -/
theorem toOp_specInvTop_apply (hY : IsTopEigFrame A hA Y lam) (x : EuclideanSpace ℝ (Fin p)) :
    toOp (specInvTop A hA r) x
      = ∑ j, ((lam j)⁻¹ * ⟪frameCol Y j, x⟫_ℝ) • frameCol Y j := by
  have hproj : specProjTop A hA r (toOp (specInvTop A hA r) x) = toOp (specInvTop A hA r) x := by
    simpa only [toCLM_eq_toOp] using proj_comp_specInvTop hA (r := r) x
  rw [← hproj, specProjTop_eq_frame_sum hY.frame]
  refine Finset.sum_congr rfl fun j _ => ?_
  congr 1
  have hsymm := symmOp (isHermitian_specInvTop hA r) (frameCol Y j) x
  rw [toOp_specInvTop_frameCol hY, real_inner_smul_left] at hsymm
  exact hsymm.symm

set_option linter.deprecated false in
/-- The trace form of `specInvTop`, read in an eigenframe. -/
theorem trace_specInvTop_eigFrame (hY : IsTopEigFrame A hA Y lam)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Bᵀ * specInvTop A hA r * B)
      = ∑ j, ∑ k, (lam j)⁻¹ * ((Yᵀ * B) j k) ^ 2 := by
  have hstep : ∀ k : Fin q, (Bᵀ * specInvTop A hA r * B) k k
      = ∑ j, (lam j)⁻¹ * ((Yᵀ * B) j k) ^ 2 := by
    intro k
    have hop : ⟪frameCol B k, toOp (specInvTop A hA r) (frameCol B k)⟫_ℝ
        = (Bᵀ * specInvTop A hA r * B) k k := by
      rw [real_inner_eq_dotProduct]
      have hlp : WithLp.ofLp (toOp (specInvTop A hA r) (frameCol B k))
          = specInvTop A hA r *ᵥ WithLp.ofLp (frameCol B k) := by
        simp [Matrix.toEuclideanLin_apply]
      rw [hlp, dotProduct, Matrix.mul_apply]
      have hL : ∀ b : Fin p, (Bᵀ * specInvTop A hA r) k b * B b k
          = ∑ a, B a k * specInvTop A hA r a b * B b k := by
        intro b
        rw [Matrix.mul_apply, Finset.sum_mul]
        exact Finset.sum_congr rfl fun a _ => rfl
      have hR : ∀ a : Fin p, WithLp.ofLp (frameCol B k) a *
            (specInvTop A hA r *ᵥ WithLp.ofLp (frameCol B k)) a
          = ∑ b, B a k * specInvTop A hA r a b * B b k := by
        intro a
        rw [Matrix.mulVec, dotProduct, Finset.mul_sum]
        exact Finset.sum_congr rfl fun b _ => by rw [← mul_assoc]; rfl
      rw [Finset.sum_congr rfl fun a (_ : a ∈ Finset.univ) => hR a,
        Finset.sum_congr rfl fun b (_ : b ∈ Finset.univ) => hL b, Finset.sum_comm]
    rw [← hop, toOp_specInvTop_apply hY, inner_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [real_inner_smul_right, inner_frameCol, real_inner_comm, inner_frameCol]
    ring
  have htr : Matrix.trace (Bᵀ * specInvTop A hA r * B)
      = ∑ k, ∑ j, (lam j)⁻¹ * ((Yᵀ * B) j k) ^ 2 := by
    rw [Matrix.trace]
    exact Finset.sum_congr rfl fun k _ => by rw [Matrix.diag_apply]; exact hstep k
  rw [htr, Finset.sum_comm]

/-- **The svdstack Frobenius identity.** With `(Y, λ)` a top-`r` eigenframe of `A` and
`λ_j ≥ 0`, the trace form of the performance is a squared Frobenius norm:
`tr(Bᵀ (specInvTop A r) B) = ‖Λ^{-1/2} Yᵀ B‖_F²`. -/
theorem frobSq_eq_trace_of_eigFrame (hY : IsTopEigFrame A hA Y lam) (hlam : ∀ j, 0 ≤ lam j)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    frobSq (Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) * (Yᵀ * B))
      = Matrix.trace (Bᵀ * specInvTop A hA r * B) := by
  rw [trace_specInvTop_eigFrame hY, frobSq]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  rw [Matrix.diagonal_mul, mul_pow, ← Real.sqrt_inv,
    Real.sq_sqrt (inv_nonneg.mpr (hlam j))]

/-- A top-`r` eigenframe transports along an equality of matrices. -/
theorem IsTopEigFrame.congr_mat {B : Matrix (Fin p) (Fin p) ℝ} (hB : B.IsHermitian)
    (h : A = B) (hY : IsTopEigFrame B hB Y lam) : IsTopEigFrame A hA Y lam := by
  subst h; exact hY

end EigFrame

/-! ### 3b. The canonical top-`r` eigenframe

`IsTopFrame A hA Y` asks `r` orthonormal columns to span `specTop A hA r`. That subspace has
dimension above `r` at a matrix whose eigenvalues tie across the index `r`, so no frame exists
there, and a hypothesis "`Y N ω` is a frame at every `ω`" is false on a Gaussian model (global
audit 2026-08-31, finding 1). Under `TopGap A hA r` no tie crosses the boundary and the first
`r` sorted eigenvectors are a frame. `topEigMat` is that canonical choice, with the index
convention of `eigenvalues₀`, of `topEigSet` and of `specProjTop`.
`specTop_simple_whp_of_tendsto` (`LinAlg/SpecProjPerturb.lean`) makes the gap hold with
probability tending to one, which is the form section 4c uses. -/

section Canonical

variable {p r : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}

/-- The sorted eigenvector of index `k`. -/
noncomputable def sortedEigvec (hA : A.IsHermitian) (k : Fin (Fintype.card (Fin p))) :
    EuclideanSpace ℝ (Fin p) :=
  (Matrix.isSymmetric_toEuclideanLin_iff.mpr hA).eigenvectorBasis finrank_euclideanSpace k

/-- The sorted eigenbasis diagonalizes the operator, with the eigenvalue read as
`eigenvalues₀`. -/
theorem toOp_sortedEigvec (hA : A.IsHermitian) (k : Fin (Fintype.card (Fin p))) :
    toOp A (sortedEigvec hA k) = hA.eigenvalues₀ k • sortedEigvec hA k :=
  (Matrix.isSymmetric_toEuclideanLin_iff.mpr hA).apply_eigenvectorBasis finrank_euclideanSpace k

/-- The sorted eigenbasis is orthonormal. -/
theorem inner_sortedEigvec (hA : A.IsHermitian) (k l : Fin (Fintype.card (Fin p))) :
    ⟪sortedEigvec hA k, sortedEigvec hA l⟫_ℝ = if k = l then 1 else 0 :=
  orthonormal_iff_ite.mp
    ((Matrix.isSymmetric_toEuclideanLin_iff.mpr hA).eigenvectorBasis
      finrank_euclideanSpace).orthonormal k l

/-- The paper's `Q_r`: the `r` sorted eigenvectors of index below `r`, as columns. -/
noncomputable def topEigMat (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) :
    Matrix (Fin p) (Fin r) ℝ :=
  Matrix.of fun i j => sortedEigvec hA (Fin.castLE hrp j) i

/-- The paper's `Λ_r`: the `r` largest eigenvalues. -/
noncomputable def topEigVal (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) :
    Fin r → ℝ := fun j => hA.eigenvalues₀ (Fin.castLE hrp j)

@[simp]
theorem frameCol_topEigMat (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) (j : Fin r) :
    frameCol (topEigMat hA hrp) j = sortedEigvec hA (Fin.castLE hrp j) := by
  ext i
  rfl

/-- `eigenvalues₀` and `eigenvalues` are the same list under the index equivalence. -/
private theorem eigenvalues₀_eq_eigF (hA : A.IsHermitian) (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k = hA.eigenvalues (Fintype.equivOfCardEq (Fintype.card_fin _) k) := by
  rw [Matrix.IsHermitian.eigenvalues, Equiv.symm_apply_apply]

/-- Under the top-`r` gap the set `topEigSet` selects exactly the indices below `r`.
`LinAlg/SpecProjPerturb.lean` has the same lemma, private. -/
private theorem mem_topEigSetF (hA : A.IsHermitian) (hgap : TopGap A hA r)
    (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k ∈ topEigSet A hA r ↔ (k : ℕ) < r := by
  constructor
  · rintro ⟨j, hj, hjk⟩
    by_contra hk
    have hlt := hgap j k hj (by omega)
    rw [hjk] at hlt
    exact lt_irrefl _ hlt
  · intro hk
    exact ⟨k, hk, rfl⟩

/-- The columns of `topEigMat` are orthonormal. -/
theorem topEigMat_ortho (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) :
    (topEigMat hA hrp)ᵀ * topEigMat hA hrp = 1 := by
  ext j k
  rw [← inner_frameCol, frameCol_topEigMat, frameCol_topEigMat, inner_sortedEigvec,
    Matrix.one_apply]
  simp [Fin.castLE_inj]

set_option linter.deprecated false in
/-- Column `j` of `topEigMat` is an eigenvector with eigenvalue `topEigVal j`. -/
theorem topEigMat_eig (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) :
    A * topEigMat hA hrp = topEigMat hA hrp * Matrix.diagonal (topEigVal hA hrp) := by
  ext i j
  have hmv : A *ᵥ WithLp.ofLp (sortedEigvec hA (Fin.castLE hrp j))
      = hA.eigenvalues₀ (Fin.castLE hrp j) •
        WithLp.ofLp (sortedEigvec hA (Fin.castLE hrp j)) := by
    simpa [Matrix.toEuclideanLin_apply] using
      congrArg WithLp.ofLp (toOp_sortedEigvec hA (Fin.castLE hrp j))
  have hL : (A * topEigMat hA hrp) i j
      = (A *ᵥ WithLp.ofLp (sortedEigvec hA (Fin.castLE hrp j))) i := by
    rw [Matrix.mul_apply, Matrix.mulVec, dotProduct]
    rfl
  have hR : (topEigMat hA hrp * Matrix.diagonal (topEigVal hA hrp)) i j
      = topEigVal hA hrp j * topEigMat hA hrp i j := by
    rw [Matrix.mul_apply]
    simp [Matrix.diagonal_apply, mul_comm]
  rw [hL, hmv, hR]
  rfl

/-- The columns of `topEigMat` span at most the top-`r` eigenspace. No gap is needed. -/
theorem span_topEigMat_le (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) :
    Submodule.span ℝ (Set.range fun j => frameCol (topEigMat hA hrp) j) ≤ specTop A hA r := by
  refine Submodule.span_le.mpr ?_
  rintro _ ⟨j, rfl⟩
  simp only [frameCol_topEigMat]
  refine Submodule.mem_iSup_of_mem (hA.eigenvalues₀ (Fin.castLE hrp j)) ?_
  refine Submodule.mem_iSup_of_mem
    (⟨Fin.castLE hrp j, by simp, rfl⟩ :
      hA.eigenvalues₀ (Fin.castLE hrp j) ∈ topEigSet A hA r) ?_
  exact Module.End.mem_eigenspace_iff.mpr (toOp_sortedEigvec hA _)

/-- Under the top-`r` gap the columns of `topEigMat` span the whole top-`r` eigenspace. -/
theorem specTop_le_span_topEigMat (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) :
    specTop A hA r ≤ Submodule.span ℝ (Set.range fun j => frameCol (topEigMat hA hrp) j) := by
  refine iSup_le fun t => iSup_le fun ht x hx => ?_
  rw [Module.End.mem_eigenspace_iff] at hx
  have hrepr : ∑ k, ⟪sortedEigvec hA k, x⟫_ℝ • sortedEigvec hA k = x :=
    ((Matrix.isSymmetric_toEuclideanLin_iff.mpr hA).eigenvectorBasis
      finrank_euclideanSpace).sum_repr' x
  rw [← hrepr]
  refine Submodule.sum_mem _ fun k _ => ?_
  by_cases hk : (k : ℕ) < r
  · refine Submodule.smul_mem _ _ (Submodule.subset_span ⟨⟨(k : ℕ), hk⟩, ?_⟩)
    simp only [frameCol_topEigMat]
    congr 1
  · have hne : hA.eigenvalues₀ k ≠ t := by
      intro h
      exact hk ((mem_topEigSetF hA hgap k).mp (by rw [h]; exact ht))
    have hz : ⟪sortedEigvec hA k, x⟫_ℝ = 0 :=
      inner_eq_zero_of_ne hA hne (toOp_sortedEigvec hA k) hx
    rw [hz, zero_smul]
    exact Submodule.zero_mem _

/-- **The canonical top-`r` frame exists as soon as `A` has a top-`r` gap.** -/
theorem isTopFrame_topEigMat (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) : IsTopFrame A hA (topEigMat hA hrp) :=
  ⟨topEigMat_ortho hA hrp,
    le_antisymm (span_topEigMat_le hA hrp) (specTop_le_span_topEigMat hA hrp hgap)⟩

/-- **The canonical top-`r` eigenframe**, the pair `(Q_r, Λ_r)` of the paper. -/
theorem isTopEigFrame_topEigMat (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) :
    IsTopEigFrame A hA (topEigMat hA hrp) (topEigVal hA hrp) :=
  ⟨isTopFrame_topEigMat hA hrp hgap, topEigMat_eig hA hrp⟩

/-- The eigenvalues of the canonical frame are nonnegative at a positive semidefinite `A`. -/
theorem topEigVal_nonneg {hA : A.IsHermitian} (hrp : r ≤ Fintype.card (Fin p))
    (hpsd : A.PosSemidef) (j : Fin r) : 0 ≤ topEigVal hA hrp j := by
  rw [topEigVal, eigenvalues₀_eq_eigF hA]
  exact hpsd.eigenvalues_nonneg _

end Canonical

/-! ### 4. The paper's estimators and the Frobenius corollaries -/

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! #### 4a. `V̂_svdstack`, unweighted -/

/-- `V̂_svdstack` at rank `r`, written as the paper writes it (`main_paper.tex:1982`):
`V̂_svdstack = Ṽᵀ Q_r Λ_r^{-1/2}` with `(Q, λ)` a top-`r` eigenframe of `Ṽ Ṽᵀ`. It is a
`d × r` matrix, and `svdGram_mul_vhatSvdstack` shows its columns are eigenvectors of `Ṽᵀ Ṽ`,
that is right singular vectors of `Ṽ`. -/
noncomputable def vhatSvdstack (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin M) (Fin r) ℝ) (lam : Fin r → ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  (m.Vt N ω)ᵀ * Q * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹

/-- The paper's display: `V̂_svdstackᵀ V = Λ^{-1/2} Q_rᵀ Ṽ V`. -/
theorem vhatSvdstack_transpose_mul_V (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin M) (Fin r) ℝ) (lam : Fin r → ℝ) :
    (m.vhatSvdstack N ω Q lam)ᵀ * m.V N
      = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) * (Qᵀ * m.VtV N ω) := by
  simp only [vhatSvdstack, UnalignedModel.VtV, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.mul_assoc]

/-- **`perfR` is the paper's `‖V̂_svdstackᵀ V‖_F²`.** This is the identification that
`RankR/Defs.lean` states in a doc-comment. -/
theorem frobSq_vhatSvdstack (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {Q : Matrix (Fin M) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gram N ω) (m.isHermitian_gram N ω) Q lam) (hlam : ∀ j, 0 ≤ lam j) :
    frobSq ((m.vhatSvdstack N ω Q lam)ᵀ * m.V N) = m.perfR N ω := by
  rw [m.vhatSvdstack_transpose_mul_V N ω Q lam, UnalignedModel.perfR]
  exact frobSq_eq_trace_of_eigFrame hQ hlam (m.VtV N ω)

/-- The `r` columns of `V̂_svdstack` are orthonormal. -/
theorem vhatSvdstack_transpose_mul_self (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {Q : Matrix (Fin M) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gram N ω) (m.isHermitian_gram N ω) Q lam) (hlam : ∀ j, 0 < lam j) :
    (m.vhatSvdstack N ω Q lam)ᵀ * m.vhatSvdstack N ω Q lam = 1 := by
  have hone : (fun j => (Real.sqrt (lam j))⁻¹ * (lam j * (Real.sqrt (lam j))⁻¹))
      = fun _ : Fin r => (1 : ℝ) := by
    funext j
    have hs : 0 < Real.sqrt (lam j) := Real.sqrt_pos.mpr (hlam j)
    have hss : Real.sqrt (lam j) * Real.sqrt (lam j) = lam j :=
      Real.mul_self_sqrt (hlam j).le
    calc (Real.sqrt (lam j))⁻¹ * (lam j * (Real.sqrt (lam j))⁻¹)
        = (Real.sqrt (lam j))⁻¹ *
            (Real.sqrt (lam j) * Real.sqrt (lam j) * (Real.sqrt (lam j))⁻¹) := by rw [hss]
      _ = ((Real.sqrt (lam j))⁻¹ * Real.sqrt (lam j)) *
            (Real.sqrt (lam j) * (Real.sqrt (lam j))⁻¹) := by ring
      _ = 1 := by
          rw [inv_mul_cancel₀ (ne_of_gt hs), mul_inv_cancel₀ (ne_of_gt hs), one_mul]
  calc (m.vhatSvdstack N ω Q lam)ᵀ * m.vhatSvdstack N ω Q lam
      = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) *
          (Qᵀ * (m.gram N ω * (Q * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹))) := by
        simp only [vhatSvdstack, UnalignedModel.gram, Matrix.transpose_mul,
          Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.mul_assoc]
    _ = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) *
          (Matrix.diagonal lam * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹) := by
        rw [← Matrix.mul_assoc (m.gram N ω), hQ.eig, Matrix.mul_assoc,
          ← Matrix.mul_assoc Qᵀ, hQ.frame.ortho, Matrix.one_mul]
    _ = 1 := by
        rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal, hone]
        simp

/-- Each column of `V̂_svdstack` is an eigenvector of `Ṽᵀ Ṽ` with eigenvalue `λ_j`: the columns
are right singular vectors of `Ṽ` for the singular values `√λ_j`. -/
theorem svdGram_mul_vhatSvdstack (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {Q : Matrix (Fin M) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gram N ω) (m.isHermitian_gram N ω) Q lam) :
    (m.Vt N ω)ᵀ * m.Vt N ω * m.vhatSvdstack N ω Q lam
      = m.vhatSvdstack N ω Q lam * Matrix.diagonal lam := by
  have hcomm : Matrix.diagonal lam * (Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹)
      = (Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹) * Matrix.diagonal lam := by
    rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
    exact congrArg Matrix.diagonal (funext fun j => mul_comm _ _)
  calc (m.Vt N ω)ᵀ * m.Vt N ω * m.vhatSvdstack N ω Q lam
      = (m.Vt N ω)ᵀ * (m.gram N ω * Q) *
          Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) := by
        simp only [vhatSvdstack, UnalignedModel.gram, Matrix.mul_assoc]
    _ = (m.Vt N ω)ᵀ * Q *
          (Matrix.diagonal lam * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹) := by
        rw [hQ.eig]
        simp only [Matrix.mul_assoc]
    _ = m.vhatSvdstack N ω Q lam * Matrix.diagonal lam := by
        rw [hcomm]
        simp only [vhatSvdstack, Matrix.mul_assoc]

/-! #### 4b. `V̂_svdstack(W)`, weighted -/

/-- `V̂_svdstack(W)` at rank `r`, as the paper writes it (`main_paper.tex:2016`):
`V̂ = Ṽ_Wᵀ Q_r Λ_r^{-1/2}` with `(Q, λ)` a top-`r` eigenframe of `Ṽ_W Ṽ_Wᵀ`. -/
noncomputable def vhatSvdstackW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (N : ℕ) (ω : Ω N) (Q : Matrix (Fin M) (Fin r) ℝ) (lam : Fin r → ℝ) :
    Matrix (Fin (d N)) (Fin r) ℝ :=
  (m.VtW W N ω)ᵀ * Q * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹

theorem vhatSvdstackW_transpose_mul_V (m : UnalignedModel μ M n d r)
    (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ) (ω : Ω N) (Q : Matrix (Fin M) (Fin r) ℝ)
    (lam : Fin r → ℝ) :
    (m.vhatSvdstackW W N ω Q lam)ᵀ * m.V N
      = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) * (Qᵀ * m.VtVW W N ω) := by
  simp only [vhatSvdstackW, UnalignedModel.VtVW, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.mul_assoc]

/-- **`perfRW W` is the paper's `‖V̂_svdstack(W)ᵀ V‖_F²`.** -/
theorem frobSq_vhatSvdstackW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (N : ℕ) (ω : Ω N) {Q : Matrix (Fin M) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gramW W N ω) (m.isHermitian_gramW W N ω) Q lam)
    (hlam : ∀ j, 0 ≤ lam j) :
    frobSq ((m.vhatSvdstackW W N ω Q lam)ᵀ * m.V N) = m.perfRW W N ω := by
  rw [m.vhatSvdstackW_transpose_mul_V W N ω Q lam, UnalignedModel.perfRW]
  exact frobSq_eq_trace_of_eigFrame hQ hlam (m.VtVW W N ω)

/-! #### 4c. The Frobenius corollaries

Each one is the corresponding proposition composed with the identity of section 4a or 4b, so
the only added hypothesis is the frame itself. The selections `Q N ω`, `lam N ω` and `Y N ω`
need no measurability: `TendstoInProb` is a statement about `μ N` of an arbitrary set (see the
header).

**Where the frame hypothesis holds.** A top-`r` frame of `A` exists exactly when
`dim (specTop A hA r) = r`, that is when no eigenvalue of `A` ties across the index `r`
boundary. On a Gaussian model with `r < M` (svdstack) or `r < d N` (stacksvd) there are `ω`
with such a tie at every `N`, so a hypothesis quantified over every `ω` is false and the
corollary says nothing (global audit 2026-08-31, finding 1). The two corollaries below
therefore ask the frame only on an event whose probability tends to one, and transport the
conclusion off the bad event with `TendstoInProb.of_tendsto_measure_ne_of_tendsto`, as
`SVDStack/Main.lean:301` does for the rank-one tie. The `_eig` companions take the canonical
frame `topEigMat` of section 3b and carry no frame hypothesis at all: for svdstack the gap of
`Ṽ Ṽᵀ` follows from `hgap` and `gramR` through `specTop_simple_whp_of_tendsto`, and for
stacksvd the gap of `X_stackᵀ X_stack` stays an explicit hypothesis, because `SubspaceLaw` has
no gap field (decision D18).

Cleanup wave 3 gives `thm_gen_rank_weight_svdstak_frobenius` the same treatment, since the
same defect appears there at `r < M`: the frame of `Ṽ_W Ṽ_Wᵀ` is asked with probability
tending to one, and the `_eig` companion reads the canonical frame, whose gap comes from
`hgapW` through `gramRW`. -/

/-- The Gram matrix `Ṽ Ṽᵀ` is positive semidefinite, so its eigenvalues are nonnegative. -/
theorem posSemidef_gram (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    (m.gram N ω).PosSemidef := by
  simpa [UnalignedModel.gram] using Matrix.posSemidef_self_mul_conjTranspose (m.Vt N ω)

/-- The top-`r` gap of `Ṽ Ṽᵀ` holds with probability tending to one. The entries converge to
`A_{β,R}` (`gramR`), which has a gap by `hgap`; `specTop_simple_whp_of_tendsto` does the rest.
This is what makes the canonical frame of section 3b available with probability tending to
one. At `r = 0` the set is empty: `TopGap A hA 0` asks nothing, because no index is below
`0`. The root lemma `specTop_simple_whp_of_tendsto` still takes `0 < r`. -/
theorem topGap_gram_whp (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ Fintype.card (Fin M))
    (hgap : TopGap (AbetaR β m.R) (isHermitian_AbetaR β m.R) r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.gram N ω) (m.isHermitian_gram N ω) r}) atTop (𝓝 0) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · have hempty : ∀ N, {ω : Ω N | ¬ TopGap (m.gram N ω) (m.isHermitian_gram N ω) 0} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact fun k _ hk _ => absurd hk (Nat.not_lt_zero _)
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds
  · exact specTop_simple_whp_of_tendsto (isHermitian_AbetaR β m.R)
      (fun N ω => m.isHermitian_gram N ω)
      (fun i j => m.gramR c β hβdef law hG i j) hr hrM hgap

/-- `prop:general_rank_unweighted_svdstack` (`main_paper.tex:799`) at `r_i = 1`, on the
paper's own quantity `‖V̂_svdstackᵀ V‖_F²`. `(Q N ω, λ N ω)` is any selection that is a
top-`r` eigenframe of `Ṽ Ṽᵀ` with nonnegative eigenvalues on an event whose probability tends
to one. Quantifying the frame over every `ω` would be a false hypothesis (see the note above
this declaration). -/
theorem prop_general_rank_unweighted_svdstack_frobenius (m : UnalignedModel μ M n d r)
    (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M) (hgap : TopGap (AbetaR β m.R) (isHermitian_AbetaR β m.R) r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise)
    (Q : ∀ N, Ω N → Matrix (Fin M) (Fin r) ℝ) (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gram N ω) (m.isHermitian_gram N ω)
      (Q N ω) (lam N ω) ∧ ∀ j, 0 ≤ lam N ω j)}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstack N ω (Q N ω) (lam N ω))ᵀ * m.V N))
      (limitR β m.R) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.prop_general_rank_unweighted_svdstack c β hc hβdef hrM hgap law hG)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hQ (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgood
  exact hω (m.frobSq_vhatSvdstack N ω hgood.1 hgood.2)

/-- The same at the canonical frame `(Q_r, Λ_r) = (topEigMat, topEigVal)` of `Ṽ Ṽᵀ`, with **no
frame hypothesis**: the hypotheses are exactly those of
`prop_general_rank_unweighted_svdstack`. So the paper's `‖V̂_svdstackᵀ V‖_F²` converges for a
definite estimator, and the corollary above is not vacuous. Second cleanup pass, 2026-09-02:
`0 < r` is gone; `topGap_gram_whp` now covers `r = 0`. -/
theorem prop_general_rank_unweighted_svdstack_frobenius_eig (m : UnalignedModel μ M n d r)
    (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ Fintype.card (Fin M))
    (hgap : TopGap (AbetaR β m.R) (isHermitian_AbetaR β m.R) r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstack N ω
        (topEigMat (m.isHermitian_gram N ω) hrM)
        (topEigVal (m.isHermitian_gram N ω) hrM))ᵀ * m.V N))
      (limitR β m.R) := by
  refine m.prop_general_rank_unweighted_svdstack_frobenius c β hc hβdef (by simpa using hrM)
    hgap law hG _ _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topGap_gram_whp c β hβdef hrM hgap law hG) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω ⟨isTopEigFrame_topEigMat (m.isHermitian_gram N ω) hrM hgapN,
    fun j => topEigVal_nonneg hrM (m.posSemidef_gram N ω) j⟩

/-- The weighted Gram matrix `Ṽ_W Ṽ_Wᵀ` is positive semidefinite. -/
theorem posSemidef_gramW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (N : ℕ) (ω : Ω N) : (m.gramW W N ω).PosSemidef := by
  simpa [UnalignedModel.gramW] using Matrix.posSemidef_self_mul_conjTranspose (m.VtW W N ω)

/-- `Ṽ_W Ṽ_Wᵀ → W A_{β,R} Wᵀ` entrywise in probability. Each entry is a fixed real linear
combination of the entries of `Ṽ Ṽᵀ` (`mul_mul_transpose_apply`), so `gramR` and the
continuous mapping theorem give it. Same step as inside
`thm_gen_rank_weight_svdstak_general`, named here for `topGap_gramW_whp`. -/
theorem gramRW (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (W : Matrix (Fin M) (Fin M) ℝ) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) (i j : Fin M) :
    TendstoInProb μ (fun N ω => m.gramW W N ω i j) (AbetaRW W β m.R i j) := by
  have hconv0 : TendstoInProbPi μ
      (fun N ω => fun t : Fin M × Fin M => m.gram N ω t.1 t.2)
      (fun t : Fin M × Fin M => AbetaR β m.R t.1 t.2) :=
    fun t => m.gramR c β hβdef law hG t.1 t.2
  have hcont : Continuous fun z : Fin M × Fin M → ℝ =>
      ∑ b : Fin M, ∑ a : Fin M, W i a * z (a, b) * W j b :=
    continuous_finsetSum _ fun b _ => continuous_finsetSum _ fun a _ =>
      (continuous_const.mul (continuous_apply _)).mul continuous_const
  have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
  have hlim : AbetaRW W β m.R i j
      = ∑ b : Fin M, ∑ a : Fin M, W i a * AbetaR β m.R a b * W j b :=
    mul_mul_transpose_apply W (AbetaR β m.R) i j
  rw [hlim]
  refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
  change ∑ b : Fin M, ∑ a : Fin M, W i a * m.gram N ω a b * W j b = m.gramW W N ω i j
  rw [m.gramW_eq]
  exact (mul_mul_transpose_apply W (m.gram N ω) i j).symm

/-- The top-`r` gap of `Ṽ_W Ṽ_Wᵀ` holds with probability tending to one, from the gap of
`W A_{β,R} Wᵀ`. Weighted twin of `topGap_gram_whp`, `r = 0` included. -/
theorem topGap_gramW_whp (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (W : Matrix (Fin M) (Fin M) ℝ) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ Fintype.card (Fin M))
    (hgapW : TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.gramW W N ω) (m.isHermitian_gramW W N ω) r}) atTop
      (𝓝 0) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · have hempty : ∀ N,
        {ω : Ω N | ¬ TopGap (m.gramW W N ω) (m.isHermitian_gramW W N ω) 0} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact fun k _ hk _ => absurd hk (Nat.not_lt_zero _)
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds
  · exact specTop_simple_whp_of_tendsto (isHermitian_AbetaRW W β m.R)
      (fun N ω => m.isHermitian_gramW W N ω) (fun i j => m.gramRW c β W hβdef law hG i j)
      hr hrM hgapW

/-- `thm:gen_rank_weight_svdstak` for one admissible weight matrix `W`
(`main_paper.tex:893`) at `r_i = 1`, on the paper's own quantity. `(Q N ω, λ N ω)` is any
selection that is a top-`r` eigenframe of `Ṽ_W Ṽ_Wᵀ` with nonnegative eigenvalues on an event
whose probability tends to one. Quantifying the frame over every `ω` would be a false
hypothesis at `r < M` (see the note above this declaration; cleanup wave 3). -/
theorem thm_gen_rank_weight_svdstak_frobenius (m : UnalignedModel μ M n d r)
    (c β : Fin M → ℝ) (W : Matrix (Fin M) (Fin M) ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hr : 0 < r) (hrM : r ≤ M)
    (hgapW : TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r)
    (hposW : 0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise)
    (Q : ∀ N, Ω N → Matrix (Fin M) (Fin r) ℝ) (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gramW W N ω)
      (m.isHermitian_gramW W N ω) (Q N ω) (lam N ω) ∧ ∀ j, 0 ≤ lam N ω j)}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackW W N ω (Q N ω) (lam N ω))ᵀ * m.V N))
      (limitRW W β m.R) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.thm_gen_rank_weight_svdstak_general c β W hβdef hr hrM hgapW hposW law hG)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hQ (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgood
  exact hω (m.frobSq_vhatSvdstackW W N ω hgood.1 hgood.2)

/-- The same at the canonical frame `(Q_r, Λ_r) = (topEigMat, topEigVal)` of `Ṽ_W Ṽ_Wᵀ`, with
**no frame hypothesis**: the hypotheses are those of `thm_gen_rank_weight_svdstak_general`.
The gap of `Ṽ_W Ṽ_Wᵀ` comes from `hgapW` through `topGap_gramW_whp`, and the eigenvalues are
nonnegative by `posSemidef_gramW`. -/
theorem thm_gen_rank_weight_svdstak_frobenius_eig (m : UnalignedModel μ M n d r)
    (c β : Fin M → ℝ) (W : Matrix (Fin M) (Fin M) ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hr : 0 < r)
    (hrM : r ≤ Fintype.card (Fin M))
    (hgapW : TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r)
    (hposW : 0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀ ⟨r - 1, by omega⟩)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackW W N ω
        (topEigMat (m.isHermitian_gramW W N ω) hrM)
        (topEigVal (m.isHermitian_gramW W N ω) hrM))ᵀ * m.V N))
      (limitRW W β m.R) := by
  refine m.thm_gen_rank_weight_svdstak_frobenius c β W hβdef hr (by simpa using hrM) hgapW
    hposW law hG _ _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topGap_gramW_whp c β W hβdef hrM hgapW law hG) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω ⟨isTopEigFrame_topEigMat (m.isHermitian_gramW W N ω) hrM hgapN,
    fun j => topEigVal_nonneg hrM (m.posSemidef_gramW W N ω) j⟩

/-- `prop:stacksvd_subspace` (`main_paper.tex:834`) at `r_i = 1`, on the paper's own quantity:
`Y N ω` is any selection that is an orthonormal basis of the top-`r` eigenspace of
`X_stackᵀ X_stack` on an event whose probability tends to one, and then
`‖V̂_stacksvdᵀ V‖_F² → ∑_j β²(√λ_j(C), ‖c‖₁)`. Here `V̂_stacksvd = Y` itself: the stacksvd
estimator is a matrix of top eigenvectors of the `d × d` Gram matrix, so no `Λ^{-1/2}`
appears. -/
theorem prop_stacksvd_subspace_frobenius (m : UnalignedModel μ M n d r) (c : Fin M → ℝ)
    (law : m.SubspaceLaw (∑ i, c i))
    (Y : ∀ N, Ω N → Matrix (Fin (d N)) (Fin r) ℝ)
    (hY : Tendsto (fun N => μ N {ω | ¬ IsTopFrame (m.stackGram N ω)
      (m.isHermitian_stackGram N ω) (Y N ω)}) atTop (𝓝 0)) :
    TendstoInProb μ (fun N ω => frobSq ((Y N ω)ᵀ * m.V N))
      (limitStackR (fun i => (m.tbl i).θ) m.R c) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ (m.prop_stacksvd_subspace c law)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hY (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hframe
  refine hω ?_
  rw [frobenius_eq_perf hframe (m.V N)]
  rfl

/-- The same at the canonical frame `topEigMat` of `X_stackᵀ X_stack`. The frame hypothesis
becomes the top-`r` gap of the stacksvd Gram matrix with probability tending to one, which is
the field that `SubspaceLaw` does not carry (decision D18). -/
theorem prop_stacksvd_subspace_frobenius_eig (m : UnalignedModel μ M n d r) (c : Fin M → ℝ)
    (law : m.SubspaceLaw (∑ i, c i)) (hrd : ∀ N, r ≤ Fintype.card (Fin (d N)))
    (hgap : Tendsto (fun N => μ N {ω | ¬ TopGap (m.stackGram N ω)
      (m.isHermitian_stackGram N ω) r}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((topEigMat (m.isHermitian_stackGram N ω) (hrd N))ᵀ * m.V N))
      (limitStackR (fun i => (m.tbl i).θ) m.R c) := by
  refine m.prop_stacksvd_subspace_frobenius c law _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hgap (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω (isTopFrame_topEigMat (m.isHermitian_stackGram N ω) (hrd N) hgapN)

end UnalignedModel

end StackedSVD
