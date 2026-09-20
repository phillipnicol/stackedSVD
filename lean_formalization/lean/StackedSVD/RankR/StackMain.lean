/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.StackGamma

/-!
# `thm:rank_r_stacksvd`: rank-`r` weighted stacksvd, Layer 1

Track D item D2 of `notes/RANK_R_PLAN.md`. The review note is
`notes/archive/rankr_D1_statement.md`; `RankR/StackGamma.lean` (item D1) carries the scalars
`γ_j`, `θ̃_jk²`, `ℓ_j` and the model objects `X_stack^(j)`, `v̂_j`, `(v_jᵀ v̂_j)²`,
`‖Vᵀ V̂‖_F²`.

## The paper

`main_paper.tex:2337`, Corollary `thm:rank_r_stacksvd`. On the exactly aligned rank-`r`
model, rank-`r` weighted stacksvd satisfies

```
(v_jᵀ v̂_{j,stacksvd})² →p γ_j     for j = 1, …, r,
‖Vᵀ V̂_stacksvd‖_F²     →p ∑_j γ_j.
```

## Layer 1

The limit laws are hypotheses, as everywhere in this development. `HeteroLawR` is the
rank-`r` twin of `MultiTableModel.HeteroLaw` (`StackSVDWeighted.lean:1267`). It has one
field per thing the corollary reads:

* `align`: the per-component projector overlap of stack `j` tends to `γ_j`. This is the
  first display, in the form that carries no simplicity side condition.
* `crossProj`: every off-diagonal projector overlap `overlapIdx X_stack^(j) ℓ_j v_k`,
  `k ≠ j`, tends to `0`. The paper states this as "the columns of `V̂` will be
  asymptotically orthogonal" (`main_paper.tex:2318`). It does **not** follow from `align`,
  and the second display needs it.
* `simpleIdxJ`: the `ℓ_j`-th eigenvalue of the Gram matrix of `X_stack^(j)` is almost surely
  simple **at that one index** (`SimpleIdx`), which turns the projector overlap into the
  paper's `⟪v̂_j, v_j⟫²`.

`HeteroLawR` has no discharge in the tree today. The rank-1 `HeteroLaw` is discharged from
`SingleTableLaw` only at `w = 1` (`heteroLaw_one_of_singleTableLaw`); at general `w` it is
black box 2 of `notes/README.md`. `HeteroLawR` is strictly harder, because
`X_stack^(j)` carries `r` spikes and `M` block-row variances `w_ij²`, so a Gaussian facade
needs a multi-spike heteroscedastic BBP law (`liu2023asymptotic` Theorem 2 with `r` spikes,
plus `hong2023optimally`).

## Hypotheses that are absent on purpose

The draft of `notes/archive/rankr_D1_statement.md` carried `hc : ∀ i, 0 < c i`, `hR : ∀ i, m.R i =
1` and the paper's separation `θ̃_jj ≠ θ̃_jk`. None of the three is read by the implication below,
so rule 5 of `CLAUDE.md` (hypothesis necessity) removes them (coordinator decision,
`notes/FLAGGED.md` item 16). They are what makes `HeteroLawR` satisfiable, not what the implication
needs, and they belong on the future Gaussian facade that discharges `HeteroLawR`. The rank-1
precedent is `hθ` on `thm_stacksvd_weighted`, which is kept only as documentation.

## Contents

1. `TendstoInProb.finsum` (`Prob/TendstoInProb.lean`, public since the 2026-09-02 dedupe)
   and `tendstoInProb_finsetSum`: a finite sum of limits in probability.
   `tendstoInProb_finsetSum` stays local: it sums over an arbitrary `Finset`, not over
   `univ`.
2. `UnalignedModelR.HeteroLawR`: the hypothesis structure.
3. `UnalignedModelR.HeteroLawR.cross`: the inner-product form of `crossProj`, by Bessel.
4. `thm_rank_r_stacksvd_proj`: the first display in projector form (`law.align`).
5. `thm_rank_r_stacksvd_inner`: the first display in the paper's inner-product form.
6. `thm_rank_r_stacksvd_frobenius`: the second display.
7. `thm_rank_r_stacksvd`: the paper-facing conjunction of 5 and 6.

No `sorry`, no new axiom.

`structure HeteroLawR` moved to `RankR/StackGamma.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. A finite sum of limits in probability -/

section FinSum

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}

/-- The same over an arbitrary `Finset`, which is what the off-diagonal sum
`∑ k ∈ Finset.univ.erase j` of `frobSqStackR_diag` needs. Terms outside `s` are replaced by
the constant `0`, and `Finset.sum_ite_mem` folds the indicator back. -/
private theorem tendstoInProb_finsetSum {ι : Type*} [Finite ι]
    (s : Finset ι) {F : ι → ∀ N, Ω N → ℝ} {v : ι → ℝ}
    (h : ∀ i ∈ s, TendstoInProb μ (F i) (v i)) :
    TendstoInProb μ (fun N ω => ∑ i ∈ s, F i N ω) (∑ i ∈ s, v i) := by
  classical
  have : Fintype ι := Fintype.ofFinite ι
  have h' : ∀ i : ι, TendstoInProb μ (fun N ω => if i ∈ s then F i N ω else 0)
      (if i ∈ s then v i else 0) := by
    intro i
    by_cases hi : i ∈ s
    · simpa only [hi, if_true] using h i hi
    · simpa only [hi, if_false] using TendstoInProb.const μ (0 : ℝ)
  have hkey := TendstoInProb.finsum h'
  have e2 : (∑ i, if i ∈ s then v i else 0) = ∑ i ∈ s, v i := by
    rw [Finset.sum_ite_mem, Finset.univ_inter]
  have e1 : ∀ (N : ℕ) (ω : Ω N),
      (∑ i, if i ∈ s then F i N ω else 0) = ∑ i ∈ s, F i N ω := by
    intro N ω
    rw [Finset.sum_ite_mem, Finset.univ_inter]
  rw [e2] at hkey
  simpa only [e1] using hkey

end FinSum

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The inner-product form of `crossProj`: every off-diagonal overlap `⟪v̂_j, v_k⟫²`,
`k ≠ j`, tends to `0`. Bessel (`overlapIdx_ge_inner_sq`, `LinAlg/SpecIdx.lean:270`) bounds
the inner form by the projector form with no simplicity hypothesis, and
`TendstoInProb.of_le` (`Prob/TendstoInProb.lean:247`) transfers the limit. This was a field
of `HeteroLawR` before the Track E audit. -/
theorem HeteroLawR.cross {m : UnalignedModelR μ M n d r (alignedRk M r)} {c : Fin M → ℝ}
    (law : m.HeteroLawR c) (j k : Fin r) (hk : k ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2) 0 := by
  refine TendstoInProb.of_le (g := fun N ω => overlapIdx (m.stackXJ c j N ω)
    (Scalars.ellR m.thetaAligned c j) (m.colVecG N k)) (fun N => ?_) (law.crossProj j k hk)
  filter_upwards with ω
  rw [m.inner_vhatStackR_sq c j k N ω, sub_zero, abs_of_nonneg (sq_nonneg _)]
  exact overlapIdx_ge_inner_sq _ _ _

/-! ### 3. The first display -/

/-- The first display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) in projector form.
This one is `law.align` and needs no simplicity. Rank-1 mirror:
`thm_stacksvd_weighted_general` (`StackSVDWeighted.lean:1277`). -/
theorem thm_rank_r_stacksvd_proj (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (law : m.HeteroLawR c) (j : Fin r) :
    TendstoInProb μ (fun N ω => m.stackOverlapJ c j N ω)
      (Scalars.gammaR m.thetaAligned c j) :=
  law.align j

/-- The first display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) in the paper's
inner-product form: `(v_jᵀ v̂_{j,stacksvd})² →p γ_j`. Route:
`stackOverlapJ_eq_inner_sq` (`RankR/StackGamma.lean:412`) on the almost sure event
`law.simpleIdxJ j N`, transferred by `TendstoInProb.congr`. Rank-1 mirror:
`thm_stacksvd_weighted_inner` (`StackSVD/Weighted.lean:64`). -/
theorem thm_rank_r_stacksvd_inner (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (law : m.HeteroLawR c) (j : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
      (Scalars.gammaR m.thetaAligned c j) := by
  refine (law.align j).congr fun N => ?_
  filter_upwards [law.simpleIdxJ j N] with ω hω
  exact m.stackOverlapJ_eq_inner_sq c j N ω hω

/-! ### 4. The second display -/

/-- The second display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`):
`‖Vᵀ V̂_stacksvd‖_F² →p ∑_j γ_j`. Route: `frobSqStackR_diag` (`RankR/StackGamma.lean:433`)
on the almost sure event where every `ℓ_j`-th eigenvalue is simple at its index (`ae_all_iff`
over the countable `Fin r`), then the `r` diagonal limits (`law.align`) plus the `r(r-1)`
off-diagonal limits (`HeteroLawR.cross`, limit `0`). -/
theorem thm_rank_r_stacksvd_frobenius (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (law : m.HeteroLawR c) :
    TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
      (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) := by
  classical
  have hdiag : TendstoInProb μ (fun N ω => ∑ j : Fin r, m.stackOverlapJ c j N ω)
      (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) :=
    TendstoInProb.finsum fun j => law.align j
  have hoffj : ∀ j : Fin r, TendstoInProb μ
      (fun N ω => ∑ k ∈ Finset.univ.erase j,
        ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2) 0 := by
    intro j
    have h0 : (0 : ℝ) = ∑ _k ∈ Finset.univ.erase j, (0 : ℝ) := by simp
    rw [h0]
    exact tendstoInProb_finsetSum _ fun k hk => law.cross j k (Finset.ne_of_mem_erase hk)
  have hoff : TendstoInProb μ
      (fun N ω => ∑ j : Fin r, ∑ k ∈ Finset.univ.erase j,
        ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2) 0 := by
    have h0 : (0 : ℝ) = ∑ _j : Fin r, (0 : ℝ) := by simp
    rw [h0]
    exact TendstoInProb.finsum fun j => hoffj j
  have hsum := hdiag.add hoff
  rw [add_zero] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards [ae_all_iff.mpr fun j : Fin r => law.simpleIdxJ j N] with ω hω
  exact (m.frobSqStackR_diag c N ω hω).symm

/-! ### 5. The corollary -/

/-- **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`), Layer 1 form: under `HeteroLawR` the
rank-`r` weighted stacksvd estimator satisfies both displays of the corollary,
`(v_jᵀ v̂_j)² →p γ_j` for every `j` and `‖Vᵀ V̂‖_F² →p ∑_j γ_j`.

`hc : ∀ i, 0 < c i`, `hR : ∀ i, m.R i = 1` and the separation `θ̃_jj ≠ θ̃_jk` are not
hypotheses here: the implication does not read them (`notes/FLAGGED.md` item 16). They are
conditions for `HeteroLawR` to hold and belong on a Gaussian facade.

Scope. The spectral index of the statement is `Scalars.ellR m.thetaAligned c j`, the paper's
`ℓ_j`, everywhere: in `stackOverlapJ` and `vhatStackR` (`RankR/StackGamma.lean:346, :362`)
and in the `crossProj` and `simpleIdxJ` fields of `HeteroLawR`. It equals the outlier index
`j` by `ellR_thetaAligned` (`RankR/StackGamma.lean:323`), which needs `hc : ∀ i, 0 < c i` and
`hex : ∃ i, 0 < θ_ij`. This theorem assumes neither, so read the index of the statement as
the paper's `ℓ_j`. `SpikedModelR.hθnn` and `hθanti` order the spikes strictly in every table
and keep them nonnegative, so on that family the paper's index is `ℓ_j = j`
(`main_paper.tex:2369`) at every component that some table carries (F8). The paper
(`main_paper.tex:2306`, with the condition and the conclusion at `:2368` and `:2369`) allows
an unordered `θ`, where `ℓ_j` can differ from the outlier index; that case needs a
permutation `R_i`, or a model class without `hθanti` (`notes/FLAGGED.md` item 17 (7)). More
than generality is at stake there: with the paper's own `ℓ_j` the first display is false on
an explicit `M = 2`, `r = 2` instance, where index 0 carries `0.571 ± 0.008` and the paper's
index 1 carries `0.011 ± 0.003` at `d = 1200` (`paper_edits.md` finding E7). -/
theorem thm_rank_r_stacksvd (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (law : m.HeteroLawR c) :
    (∀ j : Fin r, TendstoInProb μ
        (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR m.thetaAligned c j)) ∧
      TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
        (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) :=
  ⟨fun j => m.thm_rank_r_stacksvd_inner c law j, m.thm_rank_r_stacksvd_frobenius c law⟩

end UnalignedModelR

end StackedSVD
