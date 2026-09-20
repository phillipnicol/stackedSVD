/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.RankR.Weighted
import StackedSVD.RankR.Unweighted

/-!
# `thm:gen_rank_weight_svdstak`: the weighted rank-`r` svdstack limit

Moved out of `RankR/Weighted.lean` (F28, 2026-09-08):
`thm_gen_rank_weight_svdstak_general` through `_max_gaussian`, the Layer 1 and
Gaussian forms of the weighted general-rank svdstack limit. No proof changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

variable {M r : ℕ}

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- `thm:gen_rank_weight_svdstak` for one admissible weight matrix `W`, at `r_i = 1`. Same
hypotheses as `prop_general_rank_unweighted_svdstack`, with the gap and the positivity now
read on `W A_{β,R} Wᵀ`. At `W = 1` it is `prop_general_rank_unweighted_svdstack`
(`perfRW_one`, `limitRW_one`, `abetaRW_one`). The hypothesis `hc` of the first statement of
this file is not needed: `hgapW` and `hposW` carry everything the proof reads off `β`. -/
theorem thm_gen_rank_weight_svdstak_general (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (W : Matrix (Fin M) (Fin M) ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hr : 0 < r) (hrM : r ≤ M)
    (hgapW : TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r)
    (hposW : 0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW W N ω) (limitRW W β m.R) := by
  have hA : (AbetaRW W β m.R).IsHermitian := isHermitian_AbetaRW W β m.R
  have hrp : r ≤ Fintype.card (Fin M) := by simpa using hrM
  have hconv0 : TendstoInProbPi μ
      (fun N ω => Sum.elim (fun t : Fin M × Fin M => m.gram N ω t.1 t.2)
        (fun t : Fin M × Fin r => m.VtV N ω t.1 t.2))
      (Sum.elim (fun t : Fin M × Fin M => AbetaR β m.R t.1 t.2)
        (fun t : Fin M × Fin r => BR β m.R t.1 t.2)) := by
    rintro (t | t)
    · exact m.gramR c β hβdef law hG t.1 t.2
    · exact m.VtV_tendsto c β hβdef law t.1 t.2
  have hconv : TendstoInProbPi μ
      (fun N ω => Sum.elim (fun t : Fin M × Fin M => m.gramW W N ω t.1 t.2)
        (fun t : Fin M × Fin r => m.VtVW W N ω t.1 t.2))
      (Sum.elim (fun t : Fin M × Fin M => AbetaRW W β m.R t.1 t.2)
        (fun t : Fin M × Fin r => BRW W β m.R t.1 t.2)) := by
    rintro (t | t)
    · obtain ⟨i, j⟩ := t
      have hcont : Continuous fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
          ∑ b : Fin M, ∑ a : Fin M, W i a * z (Sum.inl (a, b)) * W j b :=
        continuous_finsetSum _ fun b _ => continuous_finsetSum _ fun a _ =>
          (continuous_const.mul (continuous_apply _)).mul continuous_const
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inl] at h
      have hlim : AbetaRW W β m.R i j
          = ∑ b : Fin M, ∑ a : Fin M, W i a * AbetaR β m.R a b * W j b :=
        mul_mul_transpose_apply W (AbetaR β m.R) i j
      rw [Sum.elim_inl, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inl]
      rw [m.gramW_eq, mul_mul_transpose_apply]
    · obtain ⟨i, k⟩ := t
      have hcont : Continuous fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
          ∑ a : Fin M, W i a * z (Sum.inr (a, k)) :=
        continuous_finsetSum _ fun a _ => continuous_const.mul (continuous_apply _)
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inr] at h
      have hlim : BRW W β m.R i k = ∑ a : Fin M, W i a * BR β m.R a k := Matrix.mul_apply
      rw [Sum.elim_inr, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inr]
      rw [m.VtVW_eq]
      exact (Matrix.mul_apply).symm
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hr hrp hgapW hposW (BRW W β m.R)) hconv
  rw [traceFun_eq hA (BRW W β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gramW W N ω) (m.VtVW W N ω)

/-- **`thm:gen_rank_weight_svdstak`** (`main_paper.tex:893`) at `r_i = 1`: with the weights
`W⋆ = D^{-1/2}` the performance of weighted svdstack tends to
`L⋆ = r - ∑_{ℓ=1}^{r} λ_{r̃+1-ℓ}(A_{β,R}^{-1/2} D A_{β,R}^{-1/2})`. The optimality half of the
paper's claim is `thm_gen_rank_weight_svdstak_opt`, and the two are put together in
`thm_gen_rank_weight_svdstak_max`.

The hypothesis is `rank B_R = r`, which is what the eigengap at `W⋆` needs. The paper's own
`β_ij > 0` plus `Rank(∑ R_i R_iᵀ) = r` implies it (`thm_gen_rank_weight_svdstak_paper`) and is
strictly stronger; at `r = 1` the hypothesis reads `∃ k, β_k ≠ 0`, which is the hypothesis of
the rank-one `thm_svdstack_weighted`. There is no `0 < r` hypothesis: at `r = 0` both sides
are `0` (cleanup wave 2, audit T7 item 7). -/
theorem thm_gen_rank_weight_svdstak (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M)
    (hrankB : (BR β m.R).rank = r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
      (limitROpt β m.R (by simpa using hrM)) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · -- `r = 0`: `perfRW` is the trace of a `0 × 0` matrix and `L⋆ = 0 - 0`
    have hzero : ∀ (N : ℕ) (ω : Ω N), m.perfRW (optWR β) N ω = 0 := fun N ω => by
      simp [UnalignedModel.perfRW, Matrix.trace]
    have hlim : limitROpt β m.R (Nat.zero_le _) = 0 := by simp [limitROpt]
    rw [hlim]
    exact (TendstoInProb.const μ 0).congr fun N =>
      Filter.Eventually.of_forall fun ω => (hzero N ω).symm
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have h0 : ∀ i, 0 ≤ β i := fun i => (hβ01 i).1
  have h1 : ∀ i, β i < 1 := fun i => (hβ01 i).2
  have hrp : r ≤ Fintype.card (Fin M) := by simpa using hrM
  have hgap := topGap_optWR_of_rankBR m.R h0 h1 hrankB
  have hpd := abetaRW_optWR_posDef m.R h0 h1
  have hpos : 0 < (isHermitian_AbetaRW (optWR β) β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ := by
    rw [← eigenvalues_eigIdx]
    exact hpd.eigenvalues_pos _
  have h := m.thm_gen_rank_weight_svdstak_general c β (optWR β) hβdef hr hrM hgap hpos law hG
  rwa [limitRW_optWR β m.R h0 h1 hrp] at h

/-- `thm:gen_rank_weight_svdstak` with the paper's own hypotheses: `β_ij > 0` for every `i, j`
and the rank condition of `assum:unaligned`. -/
theorem thm_gen_rank_weight_svdstak_paper (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hβpos : ∀ i, 0 < β i) (hrM : r ≤ M)
    (hrank : (Rstack m.R).rank = r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
      (limitROpt β m.R (by simpa using hrM)) :=
  m.thm_gen_rank_weight_svdstak c β hc hβdef hrM
    (rank_BR_of_ne_zero m.R (fun i => (hβpos i).ne') hrank) law hG

/-- **The optimality half of `thm:gen_rank_weight_svdstak`**: no admissible weight matrix beats
`L⋆`. This is `limitRW_le_opt` with the quantifier over `W` in front, so that the file states
the paper's "the optimal weights are `W⋆`" and not only the limit at `W⋆`. -/
theorem thm_gen_rank_weight_svdstak_opt (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hr : 0 < r) (hrM : r ≤ M) :
    ∀ W : Matrix (Fin M) (Fin M) ℝ,
      TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r →
      0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
        ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
      limitRW W β m.R ≤ limitROpt β m.R (by simpa using hrM) := by
  intro W hgapW hposW
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  exact limitRW_le_opt W β m.R (fun i => (hβ01 i).1) (fun i => (hβ01 i).2) hr
    (by simpa using hrM) hgapW hposW

/-- **`thm:gen_rank_weight_svdstak`, both halves, in one declaration.** The performance of
weighted svdstack at `W⋆ = D^{-1/2}` tends to `L⋆`; and for every admissible weight matrix `W`
the performance at `W` tends to `limitRW W`, which is at most `L⋆`. So the maximum over `W` of
the limit of `‖V̂_svdstack(W)ᵀ V‖_F²` is attained at `W⋆`, which is what the paper claims
(`main_paper.tex:893` with the definition of `W\opt` at `main_paper.tex:889`).

Cleanup wave 2 added the convergence half of the second conjunct. Before it the statement
bounded `limitRW W`, a defined quantity, and never said that `limitRW W` is the limit of the
performance at `W`; that link sat in `thm_gen_rank_weight_svdstak_general` and needed `law`
and `hG`, which this theorem already carries
(`notes/archive/audit_rank_r_weighted_post_2026-08-30.md`, finding 4.1). -/
theorem thm_gen_rank_weight_svdstak_max (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hr : 0 < r) (hrM : r ≤ M)
    (hrankB : (BR β m.R).rank = r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
        (limitROpt β m.R (by simpa using hrM)) ∧
      ∀ W : Matrix (Fin M) (Fin M) ℝ,
        TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r →
        0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRW W N ω) (limitRW W β m.R) ∧
          limitRW W β m.R ≤ limitROpt β m.R (by simpa using hrM) :=
  ⟨m.thm_gen_rank_weight_svdstak c β hc hβdef hrM hrankB law hG,
    fun W hgapW hposW =>
      ⟨m.thm_gen_rank_weight_svdstak_general c β W hβdef hr hrM hgapW hposW law hG,
        m.thm_gen_rank_weight_svdstak_opt c β hc hβdef hr hrM W hgapW hposW⟩⟩

/-! ### 7. Layer 2: the Gaussian corollaries

`thm_gen_rank_weight_svdstak` and `thm_gen_rank_weight_svdstak_max` take
`SingleTableLaw (c i)` of every table as a hypothesis. Under the proportional regime
`Regime (c i)` and the joint Gaussian law that hypothesis is a theorem
(`SpikedModel.singleTableLaw_of_gaussian` of `RMT/Full.lean`, whose Gaussian marginal is
`gaussianNoise_of_joint`). The two corollaries below carry the same conclusions with the
`SingleTableLaw` hypothesis discharged, exactly as
`prop_general_rank_unweighted_svdstack_gaussian` does for the unweighted proposition
(`RankR/Defs.lean`). The weight matrix `W⋆` needs no extra condition: its admissibility
(`TopGap` and `0 < λ_{r-1}` at `W⋆ A_{β,R} W⋆ᵀ`) comes from `hrankB` inside
`thm_gen_rank_weight_svdstak`. The `∀ W` half keeps the two admissibility hypotheses,
because an arbitrary weight matrix can be singular. -/

/-- **`thm:gen_rank_weight_svdstak`** (`main_paper.tex:893`) at `r_i = 1`, Layer 2 form: the
same conclusion as `thm_gen_rank_weight_svdstak` from the proportional regime of each table
and the joint Gaussian law, with no `SingleTableLaw` hypothesis. -/
theorem thm_gen_rank_weight_svdstak_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M)
    (hrankB : (BR β m.R).rank = r)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
      (limitROpt β m.R (by simpa using hrM)) :=
  m.thm_gen_rank_weight_svdstak c β hc hβdef hrM hrankB
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG

/-- **`thm:gen_rank_weight_svdstak`, both halves, Layer 2 form.** The performance of weighted
svdstack at `W⋆ = D^{-1/2}` tends to `L⋆`; the performance at any admissible `W` tends to
`limitRW W ≤ L⋆`; and the only hypotheses are the proportional regime and the joint Gaussian
law. Same statement as `thm_gen_rank_weight_svdstak_max` with `SingleTableLaw` discharged. -/
theorem thm_gen_rank_weight_svdstak_max_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hr : 0 < r) (hrM : r ≤ M)
    (hrankB : (BR β m.R).rank = r)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
        (limitROpt β m.R (by simpa using hrM)) ∧
      ∀ W : Matrix (Fin M) (Fin M) ℝ,
        TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r →
        0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRW W N ω) (limitRW W β m.R) ∧
          limitRW W β m.R ≤ limitROpt β m.R (by simpa using hrM) :=
  ⟨m.thm_gen_rank_weight_svdstak_gaussian c β hc hβdef hrM hrankB hreg hG,
    fun W hgapW hposW =>
      ⟨m.thm_gen_rank_weight_svdstak_general c β W hβdef hr hrM hgapW hposW
        (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
          (m.gaussianNoise_of_joint hG i)) hG,
        m.thm_gen_rank_weight_svdstak_opt c β hc hβdef hr hrM W hgapW hposW⟩⟩


end UnalignedModel

end StackedSVD
