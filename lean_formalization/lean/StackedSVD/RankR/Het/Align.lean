/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.Outliers
import StackedSVD.RankR.Het.Duality
import StackedSVD.RankR.Het.Scalars
import StackedSVD.RankR.Het.Edge
import StackedSVD.LinAlg.SpecWindow

/-!
# Track E, task E8a: the supercritical `align` and `cross` limits of the weighted stack

Paper: `thm:rank_r_stacksvd` (the weighted rank-`r` stacked SVD, `main_paper.tex:2337`), read
through `thm:stacksvd_weighted` and `eq:assumption4`. This file lands the overlap limit of the
sorted right singular direction `ellSup k` of the weighted stack `X_W` with the spike
directions `v_l` at a supercritical component `k`: the limit is `Scalars.Lw θ_k c w` at
`l = k` (the `align` field) and `0` at `l ≠ k` (the `cross` field).

The Track C mirror is `SpikedModelR.align_cross_of_gaussian_supercritical_aux`
(`RankR/RMT/TableAlignSup.lean`); the rank-one mirror is
`MultiTableModel.align_tendstoInProb_het` (`RMT/Het/R5het.lean`). Two things are new here.

1. The overlap `overlapIdx (stackXW) ℓ (colVecG N l)` reads the `d`-side matrix `X_Wᵀ X_W`,
   while the column split and the limits of task E7 live on the `n` side `X_W X_Wᵀ`. The
   duality of task E2 (`Het.overlapIdx_eq_normSq_specProjIdx_div`, `RankR/Het/Duality.lean`)
   crosses the gap at the price of a division by the sorted eigenvalue `λ_ℓ(X_W X_Wᵀ)`, so
   this file also proves that this eigenvalue tends to `rhoHet θ_k c w`
   (`UnalignedModelR.tendstoInProb_eigVal_het`), and the final limit is a quotient.
2. The count filter is the `Assumption4 θ_l c w ∧ τ < rhoHet θ_l c w` shape of task E7, and
   the sorted index is `Scalars.ellSup m.thetaAligned c w k` (task E0), the number of
   supercritical outliers strictly above `rhoHet θ_k c w`.

## The route

The sandwich of `notes/archive/rankr_TrackE_plan.md` step 7. A margin `δ > 0` separates
`ρ_k = rhoHet θ_k c w` from the bulk edge `bHet c w` and from every other supercritical
outlier (`exists_margin_het`). At the thresholds `τ∓ = ρ_k ∓ δ'` the count of supercritical
outliers above `τ₋` is `ellSup k + 1` and above `τ₊` it is `ellSup k`
(`card_filter_sub_of_gap`, `card_filter_add_of_gap`), so the count theorem of task E7 traps the
sorted eigenvalue `λ_{ellSup k}` in `(τ₋, τ₊]` with probability tending to 1. On that event
the window lemmas `specProjIdx_eq_specProj_Ioc` and `normSq_specProj_Ioc`
(`LinAlg/SpecWindow.lean`) write the index projector norm as the difference of two half-line
norms, and the half-line limits of task E7 give `1 / nuHet θ_k c w` at `l = k` and `0`
otherwise. The quotient by `λ_{ellSup k} → ρ_k` and the identity `1 / (ρ ν) = L(w)`
(`MPhet.overlap_identity_nuHet`) finish.

A tie between two supercritical outliers puts two eigenvalues in the window and breaks both
counts. The general lemma `exists_margin_het` therefore takes a separation hypothesis
`hsep : ∀ l, l ≠ k → Assumption4 θ_l c w → rhoHet θ_l c w ≠ rhoHet θ_k c w`. The two model
theorems do not: `SpikedModelR.hθnn` and `hθanti` order the spikes strictly in every table,
so `Scalars.rhoHet_lt_of_lt` derives `hsep` from `hc` and `hw`. The junk value `rhoHet = 0` of
a subcritical component never enters: every filter and every hypothesis reads `rhoHet θ_l`
under `Assumption4 θ_l`.

Numeric check: `check_sandwich.py` (agent scratch, reproduced in the report), seed 20260902; see
`notes/archive/agent_reports/d32_E8a_hetalign.md`.

Every declaration is proved in full.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. The margin around a supercritical outlier -/

/-- **The margin.** At a supercritical component `k` with no tie among the supercritical
outliers, a `δ > 0` separates `ρ_k = rhoHet θ_k c w` from the bulk edge `bHet c w` and from
every other supercritical outlier by `2 δ`. The two thresholds of the sandwich are `ρ_k ∓ δ'`
for any `0 < δ' ≤ δ`. The Track C mirror is `SpikedModelR.exists_margin`
(`RankR/RMT/TableAlignSup.lean`); there the outliers are ordered by the index, here the order
is free and only the separation `hsep` enters. A subcritical `l` costs nothing: its filter
entry is `false` and its `rhoHet` is never read. -/
theorem exists_margin_het (m : UnalignedModelR μ M n d r (alignedRk M r)) {c w : Fin M → ℝ}
    (hc : ∀ i, 0 < c i) {k : Fin r}
    (hk : Scalars.Assumption4 (fun i => m.thetaAligned i k) c w)
    (hsep : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      MPhet.rhoHet (fun i => m.thetaAligned i l) c w
        ≠ MPhet.rhoHet (fun i => m.thetaAligned i k) c w) :
    ∃ δ : ℝ, 0 < δ ∧
      MPhet.bHet c w + 2 * δ < MPhet.rhoHet (fun i => m.thetaAligned i k) c w ∧
      ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
        (MPhet.rhoHet (fun i => m.thetaAligned i l) c w + 2 * δ
            < MPhet.rhoHet (fun i => m.thetaAligned i k) c w
          ∨ MPhet.rhoHet (fun i => m.thetaAligned i k) c w + 2 * δ
            < MPhet.rhoHet (fun i => m.thetaAligned i l) c w) := by
  classical
  have hbe : MPhet.bHet c w < MPhet.rhoHet (fun i => m.thetaAligned i k) c w :=
    MPhet.bHet_lt_rhoHet hc hk
  have hg : ∀ l : Fin r, 0 <
      (if l = k ∨ ¬ Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
        then MPhet.rhoHet (fun i => m.thetaAligned i k) c w - MPhet.bHet c w
        else |MPhet.rhoHet (fun i => m.thetaAligned i l) c w
          - MPhet.rhoHet (fun i => m.thetaAligned i k) c w|) := by
    intro l
    by_cases h : l = k ∨ ¬ Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
    · rw [if_pos h]; linarith
    · rw [if_neg h, abs_pos, sub_ne_zero]
      push Not at h
      exact hsep l h.1 h.2
  obtain ⟨L, hL, hLle⟩ := ScalarsC.exists_pos_lower_bound _ hg
  refine ⟨L / 3, by linarith, ?_, fun l hlk hl => ?_⟩
  · have h1 := hLle k
    rw [if_pos (Or.inl rfl)] at h1
    linarith
  · have h1 := hLle l
    rw [if_neg (not_or.mpr ⟨hlk, not_not.mpr hl⟩)] at h1
    rcases lt_or_gt_of_ne (hsep l hlk hl) with h | h
    · left
      rw [abs_of_neg (by linarith)] at h1
      linarith
    · right
      rw [abs_of_pos (by linarith)] at h1
      linarith

/-! ### 2. The separation hypothesis of task E7 at the two thresholds -/

/-- The separation hypothesis of task E7 at the lower threshold `τ = ρ_k - δ'`, with margin
`δ'`: the outlier `ρ_k` sits at `τ + δ'` (right branch), and the gap of `exists_margin_het`
gives the branch of every other supercritical outlier. -/
theorem sep_of_gap_sub (m : UnalignedModelR μ M n d r (alignedRk M r)) {c w : Fin M → ℝ}
    {k : Fin r} {δ δ' τ : ℝ} (hδ' : 0 < δ') (hδδ : δ' ≤ δ)
    (hgap : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i k) c w
        ∨ MPhet.rhoHet (fun i => m.thetaAligned i k) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i l) c w))
    (hτ : τ = MPhet.rhoHet (fun i => m.thetaAligned i k) c w - δ') :
    ∀ l : Fin r, Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w < τ
        ∨ τ + δ' ≤ MPhet.rhoHet (fun i => m.thetaAligned i l) c w) := by
  intro l hl
  rcases eq_or_ne l k with rfl | hlk
  · right; rw [hτ]; linarith
  · rcases hgap l hlk hl with h | h
    · left; rw [hτ]; linarith
    · right; rw [hτ]; linarith

/-- The separation hypothesis of task E7 at the upper threshold `τ = ρ_k + δ'`, with margin
`δ'`: the outlier `ρ_k` sits at `τ - δ'` (left branch). -/
theorem sep_of_gap_add (m : UnalignedModelR μ M n d r (alignedRk M r)) {c w : Fin M → ℝ}
    {k : Fin r} {δ δ' τ : ℝ} (hδ' : 0 < δ') (hδδ : δ' ≤ δ)
    (hgap : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i k) c w
        ∨ MPhet.rhoHet (fun i => m.thetaAligned i k) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i l) c w))
    (hτ : τ = MPhet.rhoHet (fun i => m.thetaAligned i k) c w + δ') :
    ∀ l : Fin r, Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w < τ
        ∨ τ + δ' ≤ MPhet.rhoHet (fun i => m.thetaAligned i l) c w) := by
  intro l hl
  rcases eq_or_ne l k with rfl | hlk
  · left; rw [hτ]; linarith
  · rcases hgap l hlk hl with h | h
    · left; rw [hτ]; linarith
    · right; rw [hτ]; linarith

/-! ### 3. The two counts -/

open Classical in
/-- **The count above the lower threshold is `ellSup k + 1`.** The filter of task E7 at
`τ = ρ_k - δ'` holds the supercritical outliers strictly above `ρ_k` (the filter of
`Scalars.ellSup`) and the component `k` itself. Both filters use the `Classical` instance of
`Scalars.ellSup`. The Track C mirror is the count `hcm` inside
`SpikedModelR.align_cross_of_gaussian_supercritical_aux`. -/
theorem card_filter_sub_of_gap (m : UnalignedModelR μ M n d r (alignedRk M r))
    {c w : Fin M → ℝ} {k : Fin r} (hk : Scalars.Assumption4 (fun i => m.thetaAligned i k) c w)
    {δ δ' τ : ℝ} (hδ' : 0 < δ') (hδδ : δ' ≤ δ)
    (hgap : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i k) c w
        ∨ MPhet.rhoHet (fun i => m.thetaAligned i k) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i l) c w))
    (hτ : τ = MPhet.rhoHet (fun i => m.thetaAligned i k) c w - δ') :
    (Finset.univ.filter fun l : Fin r =>
        Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
          ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w).card
      = Scalars.ellSup m.thetaAligned c w k + 1 := by
  have hfilt : (Finset.univ.filter fun l : Fin r =>
        Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
          ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w)
      = insert k (Finset.univ.filter fun l : Fin r =>
          Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
            ∧ MPhet.rhoHet (fun i => m.thetaAligned i k) c w
              < MPhet.rhoHet (fun i => m.thetaAligned i l) c w) := by
    ext l
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert]
    constructor
    · rintro ⟨hl, hlt⟩
      rcases eq_or_ne l k with rfl | hlk
      · exact Or.inl rfl
      · right
        refine ⟨hl, ?_⟩
        rcases hgap l hlk hl with h | h
        · rw [hτ] at hlt; linarith
        · linarith
    · rintro (rfl | ⟨hl, hlt⟩)
      · exact ⟨hk, by rw [hτ]; linarith⟩
      · exact ⟨hl, by rw [hτ]; linarith⟩
  have hnot : k ∉ (Finset.univ.filter fun l : Fin r =>
      Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
        ∧ MPhet.rhoHet (fun i => m.thetaAligned i k) c w
          < MPhet.rhoHet (fun i => m.thetaAligned i l) c w) := by
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, not_and, not_lt]
    intro _
    exact le_rfl
  rw [hfilt, Finset.card_insert_of_notMem hnot]
  rfl

open Classical in
/-- **The count above the upper threshold is `ellSup k`.** The filter of task E7 at
`τ = ρ_k + δ'` is the filter of `Scalars.ellSup`: a supercritical outlier strictly above
`ρ_k` is above `ρ_k + 2 δ` by the gap. The Track C mirror is the count `hcp` inside
`SpikedModelR.align_cross_of_gaussian_supercritical_aux`. -/
theorem card_filter_add_of_gap (m : UnalignedModelR μ M n d r (alignedRk M r))
    {c w : Fin M → ℝ} {k : Fin r} {δ δ' τ : ℝ} (hδ' : 0 < δ') (hδδ : δ' ≤ δ)
    (hgap : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i l) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i k) c w
        ∨ MPhet.rhoHet (fun i => m.thetaAligned i k) c w + 2 * δ
          < MPhet.rhoHet (fun i => m.thetaAligned i l) c w))
    (hτ : τ = MPhet.rhoHet (fun i => m.thetaAligned i k) c w + δ') :
    (Finset.univ.filter fun l : Fin r =>
        Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
          ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w).card
      = Scalars.ellSup m.thetaAligned c w k := by
  rw [Scalars.ellSup]
  congr 1
  refine Finset.filter_congr fun l _ => ?_
  constructor
  · rintro ⟨hl, hlt⟩
    exact ⟨hl, by rw [hτ] at hlt; linarith⟩
  · rintro ⟨hl, hlt⟩
    refine ⟨hl, ?_⟩
    have hlk : l ≠ k := by
      rintro rfl
      exact lt_irrefl _ hlt
    rcases hgap l hlk hl with h | h
    · linarith
    · rw [hτ]; linarith

/-! ### 4. The index is in range -/

/-- `ellSup k < ∑ i, n i N`: the sorted index of the outlier sits inside the `n`-side
spectrum at every `N`, because `r ≤ n_0 N` (`SpikedModelR.rk_le_n`) and `ellSup k < r`
(`Scalars.ellSup_lt`). So `eigVal` at this index is the sorted eigenvalue, never the junk `0`.
The pattern is `UnalignedModelR.simpleIdxJ_of_gaussian` (`RankR/Het/Simplicity.lean`). -/
theorem ellSup_lt_sum_n [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c w : Fin M → ℝ) (k : Fin r) (N : ℕ) :
    Scalars.ellSup m.thetaAligned c w k < ∑ i, n i N := by
  have hrn : r ≤ n ⟨0, NeZero.pos M⟩ N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_n N
  have hsum : r ≤ ∑ i, n i N :=
    hrn.trans (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
      (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  exact lt_of_lt_of_le (Scalars.ellSup_lt m.thetaAligned c w k) hsum

/-- `ellSup k < min (∑ i, n i N) (d N)`, the range condition of the duality of task E2
(`Het.overlapIdx_eq_normSq_specProjIdx_div`). -/
theorem ellSup_lt_min [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c w : Fin M → ℝ) (k : Fin r) (N : ℕ) :
    Scalars.ellSup m.thetaAligned c w k < min (∑ i, n i N) (d N) := by
  have hrd : r ≤ d N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_d N
  exact lt_min (m.ellSup_lt_sum_n c w k N)
    (lt_of_lt_of_le (Scalars.ellSup_lt m.thetaAligned c w k) hrd)

/-! ### 5. The count event is measurable -/

/-- The count event of task E7 is a finite intersection of level sets of the sorted
eigenvalues of `X_W X_Wᵀ`, so it is measurable. The Track C mirror is
`RankRStack.measurableSet_count_Ioi` (`RankR/RMT/TableAlignSup.lean`); the level sets are
measurable by `UnalignedModelR.measurable_gramHet_eigenvalues₀`. -/
theorem measurableSet_count_Ioi_het (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (τ : ℝ) (u : ℕ) :
    MeasurableSet {ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q ↔ (q : ℕ) < u)} := by
  have hset : {ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q ↔ (q : ℕ) < u)}
      = ⋂ q : Fin (Fintype.card (Fin (∑ i, n i N))),
          {ω : Ω N | τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
            ↔ (q : ℕ) < u} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iInter]
  rw [hset]
  refine MeasurableSet.iInter fun q => ?_
  by_cases hq : (q : ℕ) < u
  · have he : {ω : Ω N | τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < u}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q)
            ⁻¹' Set.Ioi τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Ioi, hq, iff_true]
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N q measurableSet_Ioi
  · have he : {ω : Ω N | τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < u}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q)
            ⁻¹' Set.Iic τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic, hq, iff_false, not_lt]
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N q measurableSet_Iic

/-! ### 6. The sorted eigenvalue at the outlier index converges to `ρ_k` -/

/-- **The eigenvalue limit.** At a supercritical component `k` with no tie among the
supercritical outliers, the sorted eigenvalue of `X_W X_Wᵀ` at the index `ellSup k` tends in
`hsep` is not a hypothesis: inside the model class a tie is impossible. `SpikedModelR.hθnn`
and `hθanti` order the spikes strictly in every table, so `Scalars.rhoHet_lt_of_lt` gives
`rhoHet θ_l ≠ rhoHet θ_k` for every supercritical `l ≠ k` when every weight is nonzero
(`hw`).

Route: at `δ' = min δ (ε / 2)` the two count events of
task E7 (`tendsto_measure_count_Ioi_tau_het`) at `ρ_k ∓ δ'` trap the eigenvalue in
`(ρ_k - δ', ρ_k + δ']`, so `|λ - ρ_k| ≤ δ' < ε` outside the union of their complements. The
index is in range by `ellSup_lt_sum_n`, so `eigVal` is the sorted eigenvalue. The rank-one
mirror is the eigenvalue limit inside `MultiTableModel.align_tendstoInProb_het`
(`RMT/Het/R5het.lean`); Track C has no separate statement, because its duality does not
divide by the eigenvalue. -/
theorem tendstoInProb_eigVal_het [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {k : Fin r} (hk : Scalars.Assumption4 (fun i => m.thetaAligned i k) c w) :
    TendstoInProb μ
      (fun N ω => eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (m.stackXW w N ω))
        (Scalars.ellSup m.thetaAligned c w k))
      (MPhet.rhoHet (fun i => m.thetaAligned i k) c w) := by
  have hsep : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      MPhet.rhoHet (fun i => m.thetaAligned i l) c w
        ≠ MPhet.rhoHet (fun i => m.thetaAligned i k) c w := by
    intro l hlk h4l
    have hθnn : ∀ i j, 0 ≤ m.thetaAligned i j := fun i j => (m.tbl i).hθnn j
    have hθanti : ∀ i, StrictAnti (m.thetaAligned i) := fun i => (m.tbl i).hθanti
    rcases lt_or_gt_of_ne hlk with hlt | hgt
    · exact ne_of_gt
        (Scalars.rhoHet_lt_of_lt hc hθnn hθanti (NeZero.pos M) hw hlt hk)
    · exact ne_of_lt
        (Scalars.rhoHet_lt_of_lt hc hθnn hθanti (NeZero.pos M) hw hgt h4l)
  obtain ⟨δ, hδ, hgapb, hgap⟩ := m.exists_margin_het hc hk hsep
  have hw' : ∃ i, w i ≠ 0 := ⟨⟨0, NeZero.pos M⟩, hw _⟩
  intro ε hε
  obtain ⟨δ', hδ'def⟩ : ∃ t : ℝ, t = min δ (ε / 2) := ⟨_, rfl⟩
  have hδ'pos : 0 < δ' := by rw [hδ'def]; exact lt_min hδ (by linarith)
  have hδ'le : δ' ≤ δ := by rw [hδ'def]; exact min_le_left _ _
  have hδ'ε : δ' ≤ ε / 2 := by rw [hδ'def]; exact min_le_right _ _
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ,
    t = MPhet.rhoHet (fun i => m.thetaAligned i k) c w - δ' := ⟨_, rfl⟩
  obtain ⟨τp, hτpdef⟩ : ∃ t : ℝ,
    t = MPhet.rhoHet (fun i => m.thetaAligned i k) c w + δ' := ⟨_, rfl⟩
  have hτmb : MPhet.bHet c w < τm := by rw [hτmdef]; linarith
  have hτpb : MPhet.bHet c w < τp := by rw [hτpdef]; linarith
  have hsepm := m.sep_of_gap_sub hδ'pos hδ'le hgap hτmdef
  have hsepp := m.sep_of_gap_add hδ'pos hδ'le hgap hτpdef
  have hcm := m.card_filter_sub_of_gap hk hδ'pos hδ'le hgap hτmdef
  have hcp := m.card_filter_add_of_gap hδ'pos hδ'le hgap hτpdef
  -- the two count events of task E7, and their complements
  have hAm := m.tendsto_measure_count_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge hδ'pos hτmb
    hsepm
  have hAp := m.tendsto_measure_count_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge hδ'pos hτpb
    hsepp
  rw [hcm] at hAm
  rw [hcp] at hAp
  have hBm : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_count_Ioi_het w N τm
      (Scalars.ellSup m.thetaAligned c w k + 1)).nullMeasurableSet) hAm
  have hBp : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_count_Ioi_het w N τp
      (Scalars.ellSup m.thetaAligned c w k)).nullMeasurableSet) hAp
  -- the union bound
  refine tendsto_measure_zero_of_subset
    (t := fun N => ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)})ᶜ
      ∪ ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)})ᶜ)
    ?_ (tendsto_measure_zero_union hBm hBp)
  intro N ω hω
  by_cases hlow : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)
  · by_cases hhigh : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)
    · exfalso
      have hℓ : Scalars.ellSup m.thetaAligned c w k < Fintype.card (Fin (∑ i, n i N)) := by
        rw [Fintype.card_fin]; exact m.ellSup_lt_sum_n c w k N
      have hω' : ε ≤ |eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω))
          (Scalars.ellSup m.thetaAligned c w k)
          - MPhet.rhoHet (fun i => m.thetaAligned i k) c w| := hω
      rw [eigVal_eq _ _ hℓ] at hω'
      have h1 : τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨Scalars.ellSup m.thetaAligned c w k, hℓ⟩ :=
        (hlow ⟨Scalars.ellSup m.thetaAligned c w k, hℓ⟩).mpr (Nat.lt_succ_self _)
      have h2 : ¬ τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨Scalars.ellSup m.thetaAligned c w k, hℓ⟩ :=
        fun h => lt_irrefl _ ((hhigh ⟨Scalars.ellSup m.thetaAligned c w k, hℓ⟩).mp h)
      rw [not_lt] at h2
      rw [hτmdef] at h1
      rw [hτpdef] at h2
      have habs : |(isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨Scalars.ellSup m.thetaAligned c w k, hℓ⟩
          - MPhet.rhoHet (fun i => m.thetaAligned i k) c w| ≤ δ' :=
        abs_le.mpr ⟨by linarith, by linarith⟩
      linarith
    · exact Set.mem_union_right _ hhigh
  · exact Set.mem_union_left _ hlow

/-! ### 7. The target: `align` and `cross` at the supercritical index -/

/-- **Task E8a, the target** (the outlier branch of `thm:rank_r_stacksvd`,
`main_paper.tex:2337`; the rank-one mirror is `thm:stacksvd_weighted`, `:463`). At a
supercritical component `k` (`eq:assumption4` at the
weights `w`) with no tie among the supercritical outliers, the overlap of the sorted right
singular direction `ellSup k` of the weighted stack `X_W` with the spike direction `v_l`
tends in probability to `Scalars.Lw θ_k c w` at `l = k` and to `0` at `l ≠ k`. At the stack
weights `w = wStackR θ c k` the value `Lw` is `gammaR θ c k` (`Scalars.Lw_wStackR_eq_gammaR`),
which task E9 reads; this statement stops at `Lw`.

Route (the sandwich, `notes/archive/rankr_TrackE_plan.md` step 7): the two thresholds `ρ_k ∓ δ` of
`exists_margin_het`, the two count events of task E7 at both, the eigenvalue limit
`tendstoInProb_eigVal_het`, the two half-line limits of task E7
(`tendstoInProb_normSq_specProj_Ioi_tau_het`) and their difference, the quotient, and on the
intersection of the count events the duality of task E2
(`Het.overlapIdx_eq_normSq_specProjIdx_div`) with the window lemmas
`specProjIdx_eq_specProj_Ioc` and `normSq_specProj_Ioc`. The Track C mirror is
`SpikedModelR.align_cross_of_gaussian_supercritical_aux` (`RankR/RMT/TableAlignSup.lean`); the
rank-one mirror is `MultiTableModel.align_tendstoInProb_het` (`RMT/Het/R5het.lean`).

`hsep` is not a hypothesis: inside the model class a tie is impossible. `SpikedModelR.hθnn`
and `hθanti` order the spikes strictly in every table, so `Scalars.rhoHet_lt_of_lt` gives
`rhoHet θ_l ≠ rhoHet θ_k` for every supercritical `l ≠ k` when every weight is nonzero
(`hw`). -/
theorem align_cross_het_of_gaussian_sup [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {k : Fin r} (hk : Scalars.Assumption4 (fun i => m.thetaAligned i k) c w)
    (l : Fin r) :
    TendstoInProb μ
      (fun N ω => overlapIdx (m.stackXW w N ω) (Scalars.ellSup m.thetaAligned c w k)
        (m.colVecG N l))
      (if l = k then Scalars.Lw (fun i => m.thetaAligned i k) c w else 0) := by
  have hsep : ∀ l : Fin r, l ≠ k → Scalars.Assumption4 (fun i => m.thetaAligned i l) c w →
      MPhet.rhoHet (fun i => m.thetaAligned i l) c w
        ≠ MPhet.rhoHet (fun i => m.thetaAligned i k) c w := by
    intro l hlk h4l
    have hθnn : ∀ i j, 0 ≤ m.thetaAligned i j := fun i j => (m.tbl i).hθnn j
    have hθanti : ∀ i, StrictAnti (m.thetaAligned i) := fun i => (m.tbl i).hθanti
    rcases lt_or_gt_of_ne hlk with hlt | hgt
    · exact ne_of_gt
        (Scalars.rhoHet_lt_of_lt hc hθnn hθanti (NeZero.pos M) hw hlt hk)
    · exact ne_of_lt
        (Scalars.rhoHet_lt_of_lt hc hθnn hθanti (NeZero.pos M) hw hgt h4l)
  obtain ⟨δ, hδ, hgapb, hgap⟩ := m.exists_margin_het hc hk hsep
  have hw' : ∃ i, w i ≠ 0 := ⟨⟨0, NeZero.pos M⟩, hw _⟩
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ,
    t = MPhet.rhoHet (fun i => m.thetaAligned i k) c w - δ := ⟨_, rfl⟩
  obtain ⟨τp, hτpdef⟩ : ∃ t : ℝ,
    t = MPhet.rhoHet (fun i => m.thetaAligned i k) c w + δ := ⟨_, rfl⟩
  have hτmb : MPhet.bHet c w < τm := by rw [hτmdef]; linarith
  have hτpb : MPhet.bHet c w < τp := by rw [hτpdef]; linarith
  have hab : τm ≤ τp := by rw [hτmdef, hτpdef]; linarith
  have hsepm := m.sep_of_gap_sub hδ le_rfl hgap hτmdef
  have hsepp := m.sep_of_gap_add hδ le_rfl hgap hτpdef
  have hcm := m.card_filter_sub_of_gap hk hδ le_rfl hgap hτmdef
  have hcp := m.card_filter_add_of_gap hδ le_rfl hgap hτpdef
  -- step 1: the two count events of task E7, and their complements
  have hAm := m.tendsto_measure_count_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge hδ hτmb hsepm
  have hAp := m.tendsto_measure_count_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge hδ hτpb hsepp
  rw [hcm] at hAm
  rw [hcp] at hAp
  have hBm : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_count_Ioi_het w N τm
      (Scalars.ellSup m.thetaAligned c w k + 1)).nullMeasurableSet) hAm
  have hBp : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_count_Ioi_het w N τp
      (Scalars.ellSup m.thetaAligned c w k)).nullMeasurableSet) hAp
  -- step 2: the two half-line limits of task E7 and the eigenvalue limit
  have hg1 := m.tendstoInProb_normSq_specProj_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge l hδ
    hτmb hsepm
  have hg2 := m.tendstoInProb_normSq_specProj_Ioi_tau_het w c hc hw' hR hreg hG hpd hedge l hδ
    hτpb hsepp
  have hlam := m.tendstoInProb_eigVal_het w c hc hw hR hreg hG hpd hedge hk
  have hρpos : 0 < MPhet.rhoHet (fun i => m.thetaAligned i k) c w :=
    lt_trans (MPhet.bHet_pos hc hw') (MPhet.bHet_lt_rhoHet hc hk)
  -- step 3: the quotient limit, by cases on `l = k`
  have hq : TendstoInProb μ
      (fun N ω => (‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τm)
          (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2
        - ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τp)
          (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2)
        / eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω))
          (Scalars.ellSup m.thetaAligned c w k))
      (if l = k then Scalars.Lw (fun i => m.thetaAligned i k) c w else 0) := by
    rcases eq_or_ne l k with rfl | hlk
    · rw [if_pos rfl]
      rw [if_pos ⟨hk, by rw [hτmdef]; linarith⟩] at hg1
      rw [if_neg (fun h => by rw [hτpdef] at h; linarith [h.2])] at hg2
      have hLw : 1 / MPhet.nuHet (fun i => m.thetaAligned i l) c w
          / MPhet.rhoHet (fun i => m.thetaAligned i l) c w
          = Scalars.Lw (fun i => m.thetaAligned i l) c w := by
        rw [div_div, mul_comm (MPhet.nuHet _ c w) (MPhet.rhoHet _ c w),
          MPhet.overlap_identity_nuHet hc hk]
      have := (hg1.sub hg2).div hlam hρpos.ne'
      rwa [sub_zero, hLw] at this
    · rw [if_neg hlk]
      by_cases hl : Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
      · rcases hgap l hlk hl with h | h
        · rw [if_neg (fun h' => by rw [hτmdef] at h'; linarith [h'.2])] at hg1
          rw [if_neg (fun h' => by rw [hτpdef] at h'; linarith [h'.2])] at hg2
          have := (hg1.sub hg2).div hlam hρpos.ne'
          rwa [sub_zero, zero_div] at this
        · rw [if_pos ⟨hl, by rw [hτmdef]; linarith⟩] at hg1
          rw [if_pos ⟨hl, by rw [hτpdef]; linarith⟩] at hg2
          have := (hg1.sub hg2).div hlam hρpos.ne'
          rwa [sub_self, zero_div] at this
      · rw [if_neg (fun h' => hl h'.1)] at hg1
        rw [if_neg (fun h' => hl h'.1)] at hg2
        have := (hg1.sub hg2).div hlam hρpos.ne'
        rwa [sub_zero, zero_div] at this
  -- step 4: the union bound, with the duality of task E2 on the intersection
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hq
  refine tendsto_measure_zero_of_subset
    (t := fun N => ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)})ᶜ
      ∪ ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)})ᶜ)
    ?_ (tendsto_measure_zero_union hBm hBp)
  intro N ω hω
  by_cases hlow : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k + 1)
  · by_cases hhigh : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < Scalars.ellSup m.thetaAligned c w k)
    · refine absurd ?_ hω
      have hℓ : Scalars.ellSup m.thetaAligned c w k < min (∑ i, n i N) (d N) :=
        m.ellSup_lt_min c w k N
      have hℓ' : Scalars.ellSup m.thetaAligned c w k < Fintype.card (Fin (∑ i, n i N)) := by
        rw [Fintype.card_fin]; exact m.ellSup_lt_sum_n c w k N
      have hpos' : 0 < eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω))
          (Scalars.ellSup m.thetaAligned c w k) := by
        rw [eigVal_eq _ _ hℓ']
        exact lt_trans (lt_trans (MPhet.bHet_pos hc hw') hτmb)
          ((hlow ⟨Scalars.ellSup m.thetaAligned c w k, hℓ'⟩).mpr (Nat.lt_succ_self _))
      have hpos : 0 < eigVal ((m.stackXW w N ω)ᵀ * m.stackXW w N ω)
          (isHermitian_transpose_mul_self (m.stackXW w N ω))
          (Scalars.ellSup m.thetaAligned c w k) := by
        rw [Het.eigVal_gram_comm (m.stackXW w N ω) hℓ]
        exact hpos'
      change overlapIdx (m.stackXW w N ω) (Scalars.ellSup m.thetaAligned c w k) (m.colVecG N l)
        = (‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τm)
            (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2
          - ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τp)
            (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2)
          / eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
            (isHermitian_mul_transpose_self (m.stackXW w N ω))
            (Scalars.ellSup m.thetaAligned c w k)
      rw [Het.overlapIdx_eq_normSq_specProjIdx_div (m.stackXW w N ω) hpos (m.colVecG N l),
        m.stackXW_mulVec_colVecG w N ω l,
        specProjIdx_eq_specProj_Ioc (isHermitian_mul_transpose_self (m.stackXW w N ω)) hlow
          hhigh,
        normSq_specProj_Ioc (isHermitian_mul_transpose_self (m.stackXW w N ω)) hab]
    · exact Set.mem_union_right _ hhigh
  · exact Set.mem_union_left _ hlow

end UnalignedModelR

end StackedSVD
