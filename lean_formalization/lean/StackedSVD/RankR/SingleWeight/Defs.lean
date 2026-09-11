/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.StackGamma

/-!
# Single-weight stacksvd on the model: the estimator and its performance

Unit M1 of `notes/archive/singleweight_plan.md` section 4.2. Namespace `StackedSVD.UnalignedModelR`.

The paper (`main_paper.tex:2089`) defines `V̂_stacksvd(w)` as the top `r` eigenvectors of
`∑_i w_i² X_iᵀ X_i`, and measures it by `‖V̂_stacksvd(w)ᵀ V‖_F²`. The model side needs no new
object: `m.stackXW w` (`RankR/StackGamma.lean:258`) already scales block row `i` by the single
scalar `w i`, and `m.stackGramW w` is its Gram matrix `∑_i w_i² X_iᵀ X_i`
(`stackGramW_eq_sum`).

## The two definitions

* `perfSW` is `‖V̂_stacksvd(w)ᵀ V‖_F²` in **index projector** form,
  `∑_l ∑_k overlapIdx (X_stack(w)) l v_k`. A top-`r` frame of the eigenspace does not exist at
  a tie, so the projector sum is the total object; it equals the paper's quantity as soon as
  each of the `r` eigenvalues is simple at its own index (`perfSW_eq_inner_sq`). Same choice
  as `frobSqStackR` (`RankR/StackGamma.lean:422`).
* `vhatSW w l` is the `l`-th right singular vector of the weighted stack, with the paper's
  sign convention `⟪v̂_l, v_l⟫ ≥ 0`. Mirror: `vhatStackR` (`RankR/StackGamma.lean:364`).

STATUS 2026-09-05: the two definitions are real definitions; the two theorems are statements
only, with a row each in `docs/SORRIES.md`.
-/

open MeasureTheory Filter Topology

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The performance of single-weight stacksvd (`main_paper.tex:2089`): `‖V̂_stacksvd(w)ᵀ V‖_F²`
in index projector form, `∑_l ∑_k overlapIdx (X_stack(w)) l v_k`. At a simple `l`-th
eigenvalue each summand is `(v̂_lᵀ v_k)²` (`perfSW_eq_inner_sq`). -/
noncomputable def perfSW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  ∑ l : Fin r, ∑ k : Fin r, overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k)

/-- `v̂_l`, the `l`-th right singular vector of the weighted stack `X_stack(w)`, with the sign
fixed by the paper's convention `⟪v̂_l, v_l⟫ ≥ 0` exactly as `vhatStackR` fixes it
(`RankR/StackGamma.lean:364`). -/
noncomputable def vhatSW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (l : Fin r)
    (N : ℕ) (ω : Ω N) : EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ),
      m.colVecG N l⟫_ℝ then
    vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)
  else
    -vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)

/-- `perfSW` is a sum of squared norms, so it is nonnegative. -/
theorem perfSW_nonneg (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    0 ≤ m.perfSW w N ω :=
  Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => overlapIdx_nonneg _ _ _

/-- The projector form and the paper's own form agree once each of the `r` eigenvalues of the
weighted stack Gram matrix is simple at its own index. Mirror: `stackOverlapJ_eq_inner_sq`
(`RankR/StackGamma.lean:412`), which reads one index; here the hypothesis is `SimpleIdx` at
every `l < r`, because `perfSW` reads all of them. -/
theorem perfSW_eq_inner_sq (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N)
    (hs : ∀ l : Fin r,
      SimpleIdx (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)) :
    m.perfSW w N ω
      = ∑ l : Fin r, ∑ k : Fin r, ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2 := by
  rw [perfSW]
  refine Finset.sum_congr rfl fun l _ => Finset.sum_congr rfl fun k _ => ?_
  rw [overlapIdx, vhatSW]
  split_ifs with h
  · exact normSq_specProjIdx_eq_inner_sq (hs l) _
  · rw [inner_neg_left, neg_sq]
    exact normSq_specProjIdx_eq_inner_sq (hs l) _

end UnalignedModelR

end StackedSVD
