/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.RankR.Defs
import StackedSVD.SVDStack.Main

/-!
# The unweighted rank-`r` svdstack limit

Moved out of `RankR/Defs.lean` (F28, 2026-09-08): `align_inner` through
`prop_general_rank_unweighted_svdstack_gaussian`, the Layer 1 and Gaussian forms of
the general-rank unweighted svdstack limit. No proof changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! ### 3d. `align` against a deterministic direction -/

/-- `⟪v̂_i, v_i⟫ → β_i`, the signed form of `SingleTableLaw.align`. The sign convention of
`vhat` turns the square root of the overlap into the inner product itself. -/
theorem align_inner (m : UnalignedModel μ M n d r) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).SingleTableLaw ci) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, (m.tbl i).v N⟫_ℝ) (beta (m.tbl i).θ ci) := by
  have h1 : TendstoInProb μ
      (fun N ω => Real.sqrt (overlap ((m.tbl i).X N ω) ((m.tbl i).v N)))
      (Real.sqrt (betaSq (m.tbl i).θ ci)) :=
    law.align.comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
  refine h1.congr fun N => ?_
  filter_upwards [law.topSimple N] with ω hω
  rw [overlap_eq_inner_sq _ _ hω (m.mem_topSpace_vhat i N ω) (m.norm_vhat i N ω),
    Real.sqrt_sq_eq_abs, abs_of_nonneg (m.inner_vhat_nonneg i N ω)]

/-- `⟪v̂_i, q_N⟫ → 0` for a deterministic family `q` orthogonal to `v_i` with `‖q_N‖ ≤ 1`.
The normalized `q_N` is a unit vector orthogonal to `v_i`, so `SingleTableLaw.delocUniform`
bounds the probability; the case `q_N = 0` gives an empty event. -/
theorem inner_vhat_perp_tendsto (m : UnalignedModel μ M n d r) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).SingleTableLaw ci) (q : (N : ℕ) → EuclideanSpace ℝ (Fin (d N)))
    (hq : ∀ N, ⟪q N, (m.tbl i).v N⟫_ℝ = 0) (hn : ∀ N, ‖q N‖ ≤ 1) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, q N⟫_ℝ) 0 := by
  intro ε hε
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (law.delocUniform (ε ^ 2) (by positivity)) (fun _ => zero_le) (fun N => ?_)
  by_cases h0 : q N = 0
  · have hempty : {ω : Ω N | ε ≤ |⟪m.vhat i N ω, q N⟫_ℝ - 0|} = ∅ := by
      ext ω
      simp only [h0, inner_zero_right, sub_zero, abs_zero, Set.mem_ofPred_eq,
        Set.mem_empty_iff_false, iff_false, not_le]
      exact hε
    rw [hempty]
    simp
  · have hqpos : 0 < ‖q N‖ := norm_pos_iff.mpr h0
    have hmem : (‖q N‖⁻¹ • q N) ∈ (m.tbl i).orthUnit N := by
      refine ⟨?_, ?_⟩
      · rw [norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
        exact inv_mul_cancel₀ hqpos.ne'
      · rw [real_inner_smul_left, hq N, mul_zero]
    refine le_trans (measure_mono ?_) (le_iSup₂ (f := fun w (_ : w ∈ (m.tbl i).orthUnit N) =>
      μ N {ω | ε ^ 2 ≤ overlap ((m.tbl i).X N ω) w}) _ hmem)
    intro ω hω
    have h1 : ε ≤ |⟪m.vhat i N ω, q N⟫_ℝ| := by
      have h := hω
      rw [Set.mem_ofPred_eq, sub_zero] at h
      exact h
    have hqe : q N = ‖q N‖ • (‖q N‖⁻¹ • q N) := by
      rw [smul_smul, mul_inv_cancel₀ hqpos.ne', one_smul]
    have hval : ⟪m.vhat i N ω, q N⟫_ℝ
        = ‖q N‖ * ⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ := by
      conv_lhs => rw [hqe]
      rw [real_inner_smul_right]
    rw [hval, abs_mul, abs_of_nonneg (norm_nonneg _)] at h1
    have h2 : ε ≤ |⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ| := by
      nlinarith [abs_nonneg (⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ), hn N, hqpos]
    have h3 : ⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ ^ 2
        ≤ overlap ((m.tbl i).X N ω) (‖q N‖⁻¹ • q N) :=
      overlap_ge_inner_sq _ _ (m.mem_topSpace_vhat i N ω) (m.norm_vhat i N ω)
    have h4 : ε ^ 2 ≤ overlap ((m.tbl i).X N ω) (‖q N‖⁻¹ • q N) := by
      nlinarith [sq_abs (⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ), abs_nonneg
        (⟪m.vhat i N ω, ‖q N‖⁻¹ • q N⟫_ℝ)]
    exact h4

/-- `⟪v̂_i, y_N⟫ → β_i a` for a deterministic family `y` of norm at most one with
`⟪v_i, y_N⟫ = a` for every `N`. The split is the paper's
`y = (v_i v_iᵀ) y + (I - v_i v_iᵀ) y` (`main_paper.tex:1961`): the first half is
`align_inner`, the second is `inner_vhat_perp_tendsto`. -/
theorem align_inner_det (m : UnalignedModel μ M n d r) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).SingleTableLaw ci) (y : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))) (a : ℝ)
    (ha : ∀ N, ⟪(m.tbl i).v N, y N⟫_ℝ = a) (hy : ∀ N, ‖y N‖ ≤ 1) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, y N⟫_ℝ) (beta (m.tbl i).θ ci * a) := by
  have hA := (m.align_inner law).mul_const a
  have hB : TendstoInProb μ
      (fun N ω => ⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (y N)⟫_ℝ) 0 :=
    m.inner_vhat_perp_tendsto law (fun N => perpOf ((m.tbl i).v N) (y N))
      (fun N => by
        rw [real_inner_comm]
        exact inner_v_perpOf ((m.tbl i).hv N) (y N))
      (fun N => le_trans (norm_perpOf_le ((m.tbl i).hv N) (y N)) (hy N))
  have hsum := hA.add hB
  rw [add_zero] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  simp only [perpOf, inner_sub_right, real_inner_smul_right, ha N]
  ring

/-! ### 4. The statements -/

/-- `lem:general_rank_delocalization` (`main_paper.tex:1948`), off-diagonal half, signed: for
`i ≠ j`, `⟪v̂_i, v̂_j⟫ → β_i β_j ⟪R_i, R_j⟫`. The rank-`r` twin of
`MultiTableModel.lem_delocalization`. The paper splits
`v̂_iᵀ v̂_j = v̂_iᵀ v_i v_iᵀ v̂_j + v̂_iᵀ (I - v_i v_iᵀ) v̂_j`; the first term needs
`⟪v̂_j, v_i⟫ → β_j ⟪v_j, v_i⟫ = β_j ⟪R_j, R_i⟫`, which is `VtV_tendsto` in disguise, and the
second needs `delocUniform` of table `i` at the random direction of table `j`, hence Fubini
over the product law and the hypothesis `hG`. -/
theorem lem_general_rank_delocalization (m : UnalignedModel μ M n d r) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j) * ⟪m.R i, m.R j⟫_ℝ) := by
  have hAj : TendstoInProb μ (fun N ω => ⟪m.vhat j N ω, (m.tbl i).v N⟫_ℝ)
      (beta (m.tbl j).θ (c j) * ⟪m.R j, m.R i⟫_ℝ) :=
    m.align_inner_det (law j) (fun N => (m.tbl i).v N) _ (fun N => m.inner_v_v j i N)
      (fun N => le_of_eq ((m.tbl i).hv N))
  have hA := (m.align_inner (law i)).mul hAj
  have hB : TendstoInProb μ
      (fun N ω => ⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ) 0 := by
    refine TendstoInProb.of_le
      (g := fun N ω => Real.sqrt (overlap ((m.tbl i).X N ω)
        (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)))) (fun N => ?_) ?_
    · filter_upwards [(law j).topSimple N] with ω hω
      rw [sub_zero]
      have h1 : |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|
          ≤ |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ| :=
        abs_inner_perpOf_le ((m.tbl i).hv N) _ hω (m.mem_topSpace_vhat j N ω)
          (m.norm_vhat j N ω) _
      have h2 : |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ|
          ≤ Real.sqrt (overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))) := by
        rw [← Real.sqrt_sq_eq_abs]
        exact Real.sqrt_le_sqrt (overlap_ge_inner_sq _ _ (m.mem_topSpace_vhat i N ω)
          (m.norm_vhat i N ω))
      linarith
    · intro ε hε
      have hK := (law i).delocUniform (ε ^ 2) (by positivity)
      refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hK
        (fun _ => zero_le) (fun N => ?_)
      have hset : {ω : Ω N | ε ≤ |Real.sqrt (overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))) - 0|}
          = {ω : Ω N | ε ^ 2 ≤ overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))} := by
        ext ω
        simp only [sub_zero, Set.mem_ofPred_eq, abs_of_nonneg (Real.sqrt_nonneg _)]
        exact Real.le_sqrt hε.le (overlap_nonneg _ _)
      rw [hset]
      exact m.measure_deloc_le hG hij N (by positivity)
  have hsum := hA.add hB
  rw [add_zero] at hsum
  have hlim : beta (m.tbl i).θ (c i) * (beta (m.tbl j).θ (c j) * ⟪m.R j, m.R i⟫_ℝ)
      = beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j) * ⟪m.R i, m.R j⟫_ℝ := by
    rw [real_inner_comm]
    ring
  rw [← hlim]
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  simp only [perpOf, inner_sub_right, real_inner_smul_right]
  rw [real_inner_comm (m.vhat j N ω) ((m.tbl i).v N)]
  ring

/-- `lem:general_rank_delocalization`, second display: `Ṽ Ṽᵀ → A_{β,R}` entrywise and signed.
The diagonal is exact, not a limit: `v̂_i` is a unit vector, so `(Ṽ Ṽᵀ)_{ii} = 1`, and
`abetaR_diag` gives `(A_{β,R})_{ii} = 1`. The off-diagonal is
`lem_general_rank_delocalization`. -/
theorem gramR (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) (i j : Fin M) :
    TendstoInProb μ (fun N ω => m.gram N ω i j) (AbetaR β m.R i j) := by
  rcases eq_or_ne i j with rfl | hij
  · rw [abetaR_diag β m.R (m.hR i)]
    refine (TendstoInProb.const μ 1).congr fun N => ?_
    filter_upwards with ω
    rw [m.gram_eq_inner N ω i i, real_inner_self_eq_norm_sq, m.norm_vhat i N ω, one_pow]
  · have hd : AbetaR β m.R i j = β i * β j * ⟪m.R i, m.R j⟫_ℝ := by
      rw [abetaR_apply, if_neg hij, add_zero]
    rw [hd, hβdef i, hβdef j]
    refine (m.lem_general_rank_delocalization c law hG hij).congr fun N => ?_
    filter_upwards with ω
    exact (m.gram_eq_inner N ω i j).symm

/-- `lem:general_rank_delocalization`, first display: `Ṽ V → B_R` entrywise and signed,
`⟪v̂_i, V e_k⟫ → β_i (R_i)_k`. Route: `SingleTableLaw.align` of table `i` gives the component
along `v_i = V R_i`, with the sign fixed by the convention of `vhat`, and
`SingleTableLaw.deloc_seq` at the deterministic unit vector
`(V e_k - (R_i)_k v_i) / ‖V e_k - (R_i)_k v_i‖` kills the rest. Split off the degenerate case
`V e_k = ± v_i`, where that vector is `0`. Then `(V R_i)ᵀ V = R_iᵀ` by `hV`. This needs no
independence across tables, so `hG` does not appear. -/
theorem VtV_tendsto (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (i : Fin M) (k : Fin r) :
    TendstoInProb μ (fun N ω => m.VtV N ω i k) (BR β m.R i k) := by
  have hBR : BR β m.R i k = beta (m.tbl i).θ (c i) * m.R i k := by
    rw [BR, Matrix.of_apply, hβdef i]
  rw [hBR]
  have h := m.align_inner_det (law i) (fun N => m.colVec N k) (m.R i k)
    (fun N => m.inner_v_colVec i k N) (fun N => le_of_eq (m.norm_colVec N k))
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact (m.VtV_eq_inner N ω i k).symm

/-- `prop:general_rank_unweighted_svdstack` (`main_paper.tex:799`) at `r_i = 1`, Layer 1 form.
The paper's hypothesis that each `Θ_i` has `r_i` distinct entries is vacuous here. `hgap` is
the paper's `λ_r(A_{β,R}) - λ_{r+1}(A_{β,R}) > 0` with its convention `λ_{r̃+1} := -∞`, which
`TopGap` encodes by quantifying over the indices. `hc` gives `0 ≤ β_i < 1`, so
`A_{β,R} ⪰ D ≻ 0` and the inverse in `specInvTop` never meets a zero eigenvalue. There is no
`0 < r` hypothesis: at `r = 0` both sides are the trace of a `0 × 0` matrix (cleanup wave 2,
audit T7 item 7).

Route: `gramR` and `VtV_tendsto` give the two entrywise limits, jointly on the index type
`(Fin M × Fin M) ⊕ (Fin M × Fin r)`; `perfR` is `traceFun r` of that family at every `ω`; and
`continuousAt_trace_specInvTop` with `TendstoInProbPi.comp_continuous` transports the limit.
No good event is needed: the gap and `0 < λ_{r-1}` are hypotheses on the deterministic limit
`A_{β,R}`, not on `Ṽ Ṽᵀ`. -/
theorem prop_general_rank_unweighted_svdstack (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M)
    (hgap : TopGap (AbetaR β m.R) (isHermitian_AbetaR β m.R) r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω) (limitR β m.R) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · -- `r = 0`: `perfR` and `limitR` are traces of `0 × 0` matrices, so both are `0`
    have hzero : ∀ (N : ℕ) (ω : Ω N), m.perfR N ω = 0 := fun N ω => by
      simp [UnalignedModel.perfR, Matrix.trace]
    have hlim : limitR β m.R = 0 := by simp [limitR, Matrix.trace]
    rw [hlim]
    exact (TendstoInProb.const μ 0).congr fun N =>
      Filter.Eventually.of_forall fun ω => (hzero N ω).symm
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have hA : (AbetaR β m.R).IsHermitian := isHermitian_AbetaR β m.R
  have hpd : (AbetaR β m.R).PosDef :=
    abetaR_posDef m.R (fun k => (hβ01 k).1) (fun k => (hβ01 k).2)
  have hrp : r ≤ Fintype.card (Fin M) := by simpa using hrM
  have hrm : r - 1 < Fintype.card (Fin M) := by
    simp only [Fintype.card_fin]
    omega
  have hpos : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ := by
    set j : Fin M := Fintype.equivOfCardEq (Fintype.card_fin _) ⟨r - 1, hrm⟩ with hj
    have hev : hA.eigenvalues j = hA.eigenvalues₀ ⟨r - 1, hrm⟩ := by
      rw [Matrix.IsHermitian.eigenvalues, hj, Equiv.symm_apply_apply]
    rw [← hev]
    exact hpd.eigenvalues_pos j
  have hconv : TendstoInProbPi μ
      (fun N ω => Sum.elim (fun t : Fin M × Fin M => m.gram N ω t.1 t.2)
        (fun t : Fin M × Fin r => m.VtV N ω t.1 t.2))
      (Sum.elim (fun t : Fin M × Fin M => AbetaR β m.R t.1 t.2)
        (fun t : Fin M × Fin r => BR β m.R t.1 t.2)) := by
    rintro (t | t)
    · exact m.gramR c β hβdef law hG t.1 t.2
    · exact m.VtV_tendsto c β hβdef law t.1 t.2
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hr hrp hgap hpos (BR β m.R)) hconv
  rw [traceFun_eq hA (BR β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gram N ω) (m.VtV N ω)

/-- `prop:general_rank_unweighted_svdstack` at `r_i = 1`, Layer 2 form: the same conclusion
from the proportional regime of each table and the joint Gaussian law, with no
`SingleTableLaw` hypothesis. The proof is `prop_general_rank_unweighted_svdstack` applied to
`SpikedModel.singleTableLaw_of_gaussian` of each table, whose Gaussian marginal comes from
`gaussianNoise_of_joint`. -/
theorem prop_general_rank_unweighted_svdstack_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M)
    (hgap : TopGap (AbetaR β m.R) (isHermitian_AbetaR β m.R) r)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω) (limitR β m.R) :=
  m.prop_general_rank_unweighted_svdstack c β hc hβdef hrM hgap
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG


end UnalignedModel

end StackedSVD
