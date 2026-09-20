/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Defs
import StackedSVD.RankR.SingleWeight.Scalars

/-!
# `prop:gen_rank_stacksvd_singleweight`, Layer 1

Unit M2 of `notes/archive/singleweight_plan.md` section 4.2. `RankR/SingleWeight/Defs.lean` carries
the estimator `vhatSW` and the performance `perfSW`; `RankR/SingleWeight/Scalars.lean` carries
the scalar layer `secMat`, `swTerm`, `swLimit` and `EigSep`.

## The paper

`main_paper.tex:2112`, `prop:gen_rank_stacksvd_singleweight`. Under
`assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2104`),

```
‖V̂_stacksvd(w)ᵀ V‖_F² →p ∑_ℓ (1 - ∑_i c_i w_i⁴/(γ_ℓ - w_i²)²)
                          / (γ_ℓ z_ℓᵀ [∑_i w_i²/(w_i² - γ_ℓ)² R_i Θ_i² R_iᵀ] z_ℓ).
```

## Layer 1

The limit law is a hypothesis, as everywhere in this development (`CLAUDE.md` hard rule 2).
`SingleWeightLaw` is the single-weight twin of `HeteroLawR` (`RankR/StackGamma.lean:473`), with
one field per thing the proposition reads. `EigSep` is not a field of it: `EigSep` is
deterministic, and it is the hypothesis under which the law holds.

STATUS 2026-09-05: both theorems are proved (Layer 1, user OK on the statements the same day;
`notes/archive/singleweight_plan.md` section 6 records the split). The Gaussian discharge of
`SingleWeightLaw` is Track G (`notes/archive/trackG_plan.md`).
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The RMT input of `prop:gen_rank_stacksvd_singleweight` (`main_paper.tex:2112`), in the
shape of `HeteroLawR` (`RankR/StackGamma.lean:473`). One field per thing the proposition reads.

* `align` is the per-pair limit of the display at `main_paper.tex:2153`:
  `(v̂_lᵀ v_k)² → η_l (z_l)_k² / (γ_l z_lᵀ secDerivMat z_l)`, that is
  `swTerm … * (z_l)_k²`. It is in **projector** form, so it needs no simplicity, and it is
  stronger than the inner form by Bessel (`overlapIdx_ge_inner_sq`,
  `LinAlg/SpecIdx.lean:270`).
* `simpleIdx` is what turns the projector overlap into the paper's `(v̂_lᵀ v_k)²`. It asks
  simplicity at the one index `l` (`SimpleIdx`, `LinAlg/SpecIdxPerturb.lean:484`), which is
  the form the Gaussian discharge supplies.

`EigSep` (`RankR/SingleWeight/Scalars.lean`) is not a field: it is deterministic, and it is
the condition under which these limits hold. -/
structure SingleWeightLaw (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : Prop where
  /-- per pair `(l, k)`: the squared overlap of `v_k` with the `l`-th right singular subspace
  of `X_stack(w)` tends to `swTerm … (γ l) (z l) * (z l)_k²` -/
  align : ∀ l k : Fin r, TendstoInProb μ
    (fun N ω => overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k))
    (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
      (WithLp.ofLp (z l) k) ^ 2)
  /-- the `l`-th eigenvalue of the weighted stack Gram matrix is almost surely simple at that
  index -/
  simpleIdx : ∀ (l : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N),
    SimpleIdx (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)

/-- **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`), second half, Layer 1
form: under `EigSep` and the limit law the single-weight stacksvd performance converges in
probability to the paper's sum at `main_paper.tex:2119`.

The first half of the proposition, the existence of the unit eigenvector `z_ℓ`, is the
deterministic theorem `SingleWeight.exists_unit_eigvec_secMat`
(`RankR/SingleWeight/Scalars.lean`); it carries no probability and is stated on its own.

Proof plan (3 lines). Sum `law.align` over the `r²` pairs with
`StackedSVD.TendstoInProb.finsum` (`Prob/TendstoInProb.lean:396`); the inner sum over `k` of
`(z l)_k²` is `‖z l‖² = 1` by `hsep.eigvec`; so the total is `∑_l swTerm … = swLimit`. Same
shape as `thm_rank_r_stacksvd_frobenius` (`RankR/StackMain.lean:165`). -/
theorem prop_gen_rank_stacksvd_singleweight
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z)
    (law : m.SingleWeightLaw w c γ z) :
    TendstoInProb μ (fun N ω => m.perfSW w N ω)
      (SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z) := by
  have hz : ∀ l : Fin r, ∑ k : Fin r, (WithLp.ofLp (z l) k) ^ 2 = 1 := fun l => by
    have h := EuclideanSpace.real_norm_sq_eq (z l)
    rw [(hsep.eigvec l).1, one_pow] at h
    exact h.symm
  have hsum : TendstoInProb μ (fun N ω => ∑ l : Fin r, ∑ k : Fin r,
      overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k))
      (∑ l : Fin r, ∑ k : Fin r,
        SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
          (WithLp.ofLp (z l) k) ^ 2) :=
    TendstoInProb.finsum fun l => TendstoInProb.finsum fun k => law.align l k
  have hval : (∑ l : Fin r, ∑ k : Fin r,
        SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
          (WithLp.ofLp (z l) k) ^ 2)
      = SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z := by
    rw [SingleWeight.swLimit]
    exact Finset.sum_congr rfl fun l _ => by rw [← Finset.mul_sum, hz l, mul_one]
  rw [hval] at hsum
  exact hsum

/-- **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`) in the paper's own
inner-product form, one pair at a time: `(v̂_lᵀ v_k)² →p swTerm … (γ_l) (z_l) (z_l)_k²`.
Mirror: `thm_rank_r_stacksvd_inner` (`RankR/StackMain.lean:150`). Route: `law.align` on the
almost sure event `law.simpleIdx l N`, through `normSq_specProjIdx_eq_inner_sq`. No `EigSep`
hypothesis: the pairwise form needs only the law (hypothesis scan of 2026-09-05, after the
user's OK; the mirror carries none either). -/
theorem prop_gen_rank_stacksvd_singleweight_inner
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (law : m.SingleWeightLaw w c γ z) (l k : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2)
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
        (WithLp.ofLp (z l) k) ^ 2) := by
  refine (law.align l k).congr fun N => ?_
  filter_upwards [law.simpleIdx l N] with ω hω
  rw [overlapIdx, vhatSW]
  split_ifs with h
  · exact normSq_specProjIdx_eq_inner_sq hω _
  · rw [inner_neg_left, neg_sq]
    exact normSq_specProjIdx_eq_inner_sq hω _

end UnalignedModelR

end StackedSVD
