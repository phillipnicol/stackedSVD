/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.AlignTauR
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.LinAlg.SpecWindow
import StackedSVD.RankR.RMT.ShiftR
import StackedSVD.RankR.RMT.EdgeTauR

/-!
# Task C5: `align` and `cross` at a supercritical index of one rank-`rk` table

Task C5 of `notes/archive/rankr_plan_C.md` section 2. The rank-1 mirror is
`SpikedModel.align_of_gaussian_supercritical` (`RMT/Full.lean`); the subcritical twin of this
file is `SpikedModelR.align_cross_of_gaussian_subcritical` (`RankR/RMT/R6R.lean`, task C9).

One theorem covers both fields of `SpikedModelR.TableLawR` at a supercritical index `k`: the
limit is `betaSq (θ k) c` at `l = k` (the `align` field) and `0` at `l ≠ k` (the `cross`
field).

The route is the sandwich of the plan. Two thresholds `τ∓ = ρ_k ∓ δ` trap the sorted
eigenvalue `λ_k`, because the count of outliers above `τ₋` is `k + 1` and the count above
`τ₊` is `k`. On the event where both counts hold at every sorted index (task C3.2, an event
of probability tending to 1) the window lemma
`overlapIdx_eq_normSq_specProj_Ioi_sub` (`LinAlg/SpecWindow.lean`, task C4) writes the
overlap as a difference of two half-line overlaps, and task C3.1 gives the limit of each.

Contents:

1. `SpikedModelR.rhoSq_lt_rhoSq_of_lt`, `rhoSq_lt_rhoSq_of_gt`, `exists_margin`: the outliers
   `ρ_j = rhoSq (θ j) c` are strictly ordered around a supercritical index `k`, with a
   positive margin `δ` that also clears the bulk edge. A subcritical index costs nothing:
   its outlier sits at `bulkEdge c`, below `ρ_k`.
2. `SpikedModelR.card_filter_rhoSq_coreEig`: the count of outliers above a threshold, read on
   the stack indices, equals the count read on the model indices. The transport is the
   permutation `σ` of task C1.3.
3. `RankRStack.measurableSet_count_Ioi`: the count event of task C3.2 is measurable, so
   `tendsto_measure_compl_zero` turns its probability limit `1` into a bad event of
   probability limit `0`.
4. `SpikedModelR.align_cross_of_gaussian_supercritical_aux` and
   `SpikedModelR.align_cross_of_gaussian_supercritical`: the target, with and without the
   block side condition `∀ N, n N = rk + pp N`. The tail shift of task C1 removes it.

The supercriticality hypothesis is **per index**: `hsup : c < m.θ k ^ 4` for the index `k` of
the statement, not for every index. The mixed case costs nothing, because a subcritical spike
has `rhoSq = bulkEdge c` and stays below `τ₋` on its own. Task C10 is then a plain case split
against `SpikedModelR.align_cross_of_gaussian_subcritical` at one index.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 1. The count event is measurable -/

/-- The count event of task C3.2 is a finite intersection of level sets of the sorted
eigenvalues, so it is measurable. `RankRStack.measurableSet_eigenvalues₀_le`
(`RankR/RMT/EdgeR.lean`) is the same statement for the edge event; both stand on
`RankRStack.measurable_gram_eigenvalues₀`. Task C5 needs it because C3.2 gives a probability
limit `1` and the assembly needs a bad event of probability limit `0`. -/
theorem measurableSet_count_Ioi (s : RankRStack μ ns d r) (N : ℕ) (τ : ℝ) (u : ℕ) :
    MeasurableSet {ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
      (τ < (s.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < u)} := by
  have hset : {ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
        (τ < (s.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < u)}
      = ⋂ q : Fin (Fintype.card (Fin (d N))),
          {ω : Ω N | τ < (s.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < u} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iInter]
  rw [hset]
  refine MeasurableSet.iInter fun q => ?_
  by_cases hq : (q : ℕ) < u
  · have he : {ω : Ω N | τ < (s.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < u}
        = (fun ω => (s.isHermitian_gram N ω).eigenvalues₀ q) ⁻¹' Set.Ioi τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Ioi, hq, iff_true]
    rw [he]
    exact s.measurable_gram_eigenvalues₀ N q measurableSet_Ioi
  · have he : {ω : Ω N | τ < (s.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < u}
        = (fun ω => (s.isHermitian_gram N ω).eigenvalues₀ q) ⁻¹' Set.Iic τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic, hq, iff_false, not_lt]
    rw [he]
    exact s.measurable_gram_eigenvalues₀ N q measurableSet_Iic

end RankRStack

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-! ### 2. The outliers around a supercritical index -/

/-- Above a supercritical index the outlier is larger: `j < k` gives `θ k < θ j`, hence
`c < θ j ^ 4`, hence `ρ_k < ρ_j` by `ScalarsC.rhoSq_lt_rhoSq`. -/
theorem rhoSq_lt_rhoSq_of_lt (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) {k : Fin rk}
    (hsup : c < m.θ k ^ 4) {j : Fin rk} (hjk : j < k) :
    rhoSq (m.θ k) c < rhoSq (m.θ j) c :=
  ScalarsC.rhoSq_lt_rhoSq hc (m.θ_pos_of_sup hc hsup) hsup (m.hθanti hjk)

/-- Below a supercritical index the outlier is smaller. A subcritical spike sits at the bulk
edge (`rhoSq = bulkEdge c`), which is below `ρ_k` by `MP.bulkEdge_lt_rhoSq`; a supercritical
one is below `ρ_k` by `ScalarsC.rhoSq_lt_rhoSq`. This is where the per-index form of the
supercriticality hypothesis costs nothing. -/
theorem rhoSq_lt_rhoSq_of_gt (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) {k : Fin rk}
    (hsup : c < m.θ k ^ 4) {j : Fin rk} (hkj : k < j) :
    rhoSq (m.θ j) c < rhoSq (m.θ k) c := by
  have hθlt : m.θ j < m.θ k := m.hθanti hkj
  rcases lt_or_ge c (m.θ j ^ 4) with hsj | hsj
  · exact ScalarsC.rhoSq_lt_rhoSq hc (m.θ_pos_of_sup hc hsj) hsj hθlt
  · have h0 : rhoSq (m.θ j) c = bulkEdge c := by
      rw [rhoSq, if_neg (not_lt.mpr hsj)]
    rw [h0]
    exact MP.bulkEdge_lt_rhoSq hc (m.θ_pos_of_sup hc hsup) hsup

/-- **Step 3 of the C5 route.** A margin `δ > 0` that separates the outlier `ρ_k` of a
supercritical index from every other outlier and from the bulk edge, by `2 δ` on each side.
The two thresholds of the sandwich are `ρ_k ∓ δ`. -/
theorem exists_margin (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) {k : Fin rk}
    (hsup : c < m.θ k ^ 4) :
    ∃ δ : ℝ, 0 < δ ∧
      (∀ j : Fin rk, j < k → rhoSq (m.θ k) c + 2 * δ < rhoSq (m.θ j) c) ∧
      (∀ j : Fin rk, k < j → rhoSq (m.θ j) c + 2 * δ < rhoSq (m.θ k) c) ∧
      bulkEdge c + 2 * δ < rhoSq (m.θ k) c := by
  have hbe : bulkEdge c < rhoSq (m.θ k) c := MP.bulkEdge_lt_rhoSq hc (m.θ_pos_of_sup hc hsup) hsup
  have hg : ∀ j : Fin rk, 0 < (if j = k then rhoSq (m.θ k) c - bulkEdge c
      else |rhoSq (m.θ j) c - rhoSq (m.θ k) c|) := by
    intro j
    by_cases hjk : j = k
    · rw [if_pos hjk]; linarith
    · rw [if_neg hjk, abs_pos, sub_ne_zero]
      rcases lt_or_gt_of_ne hjk with h | h
      · exact (m.rhoSq_lt_rhoSq_of_lt hc hsup h).ne'
      · exact (m.rhoSq_lt_rhoSq_of_gt hc hsup h).ne
  obtain ⟨L, hL, hLle⟩ := ScalarsC.exists_pos_lower_bound _ hg
  refine ⟨L / 3, by linarith, fun j hjk => ?_, fun j hkj => ?_, ?_⟩
  · have hlt := m.rhoSq_lt_rhoSq_of_lt hc hsup hjk
    have h1 := hLle j
    rw [if_neg (ne_of_lt hjk), abs_of_pos (by linarith)] at h1
    linarith
  · have hlt := m.rhoSq_lt_rhoSq_of_gt hc hsup hkj
    have h1 := hLle j
    rw [if_neg (ne_of_gt hkj), abs_of_neg (by linarith)] at h1
    linarith
  · have h1 := hLle k
    rw [if_pos rfl] at h1
    linarith

/-! ### 3. The count transported over the permutation of task C1.3 -/

/-- The number of outliers above a threshold does not depend on the index type: the stack
index `i` and the model index `j` are matched by the permutation `σ` of
`SpikedModelR.exists_perm_coreEig`, and `√(coreEig (σ j)) = θ j` because `θ j` is nonnegative. -/
theorem card_filter_rhoSq_coreEig (m : SpikedModelR μ n d rk) {σ : Equiv.Perm (Fin rk)}
    (hσ : ∀ j, m.toStack.coreEig (σ j) = m.θ j ^ 2) (c t : ℝ) :
    (Finset.univ.filter fun i : Fin rk =>
        t < rhoSq (Real.sqrt (m.toStack.coreEig i)) c).card
      = (Finset.univ.filter fun j : Fin rk => t < rhoSq (m.θ j) c).card := by
  refine (Finset.card_equiv σ fun j => ?_).symm
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  rw [hσ j, Real.sqrt_sq (m.hθnn j)]

/-- The sorted spectrum of the stack Gram matrix is the sorted spectrum of `Xᵀ X`: the two
Hermitian witnesses are the same proof of the same statement. This is the bridge between the
count event of task C3.2 and the window lemma of task C4. -/
theorem eigenvalues₀_toStack_gram (m : SpikedModelR μ n d rk) (N : ℕ) (ω : Ω N) :
    (m.toStack.isHermitian_gram N ω).eigenvalues₀
      = (isHermitian_transpose_mul_self (m.X N ω)).eigenvalues₀ := rfl

/-! ### 4. The limit at a supercritical index -/

/-- **Task C5.1.** At a supercritical index `k` (`c < θ_k^4`) of a Gaussian rank-`rk` table the
overlap of the `k`-th singular subspace with the spike direction `l` converges in probability
to `betaSq (θ k) c` when `l = k` and to `0` when `l ≠ k`. This is the `align` field and the
`cross` field of `SpikedModelR.TableLawR` at one index.

The side condition `∀ N, n N = rk + pp N` with `0 < pp N` is the block-nonempty hypothesis of
the Track A chain; `align_cross_of_gaussian_supercritical` removes it by the tail shift.

Route: the two thresholds `ρ_k ∓ δ` of `exists_margin`, the count event of task C3.2 at both,
the window lemma `overlapIdx_eq_normSq_specProj_Ioi_sub` of task C4 on that event, and the
two half-line limits of task C3.1. -/
theorem align_cross_of_gaussian_supercritical_aux [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop) (hns : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, n N = rk + pp N) (hp : ∀ N, 0 < pp N)
    {k : Fin rk} (hsup : c < m.θ k ^ 4) (l : Fin rk) :
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l))
      (if k = l then betaSq (m.θ k) c else 0) := by
  obtain ⟨σ, hσ, hσv⟩ := m.exists_perm_coreEig
  obtain ⟨δ, hδ, hgap1, hgap2, hgapb⟩ := m.exists_margin hc hsup
  have hρσ : ∀ j : Fin rk,
      rhoSq (Real.sqrt (m.toStack.coreEig (σ j))) c = rhoSq (m.θ j) c := by
    intro j
    rw [hσ j, Real.sqrt_sq (m.hθnn j)]
  have hβσ : ∀ j : Fin rk,
      betaSq (Real.sqrt (m.toStack.coreEig (σ j))) c = betaSq (m.θ j) c := by
    intro j
    rw [hσ j, Real.sqrt_sq (m.hθnn j)]
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ, t = rhoSq (m.θ k) c - δ := ⟨_, rfl⟩
  obtain ⟨τp, hτpdef⟩ : ∃ t : ℝ, t = rhoSq (m.θ k) c + δ := ⟨_, rfl⟩
  have hτmb : bulkEdge c < τm := by rw [hτmdef]; linarith
  have hτpb : bulkEdge c < τp := by rw [hτpdef]; linarith
  have hab : τm ≤ τp := by rw [hτmdef, hτpdef]; linarith
  -- step 3: the separation hypothesis of task C3 at both thresholds
  have hsepm : ∀ i : Fin rk, rhoSq (Real.sqrt (m.toStack.coreEig i)) c < τm
      ∨ τm + δ ≤ rhoSq (Real.sqrt (m.toStack.coreEig i)) c := by
    intro i
    obtain ⟨j, rfl⟩ : ∃ j, σ j = i := ⟨σ.symm i, σ.apply_symm_apply i⟩
    rw [hρσ j, hτmdef]
    rcases lt_trichotomy j k with h | h | h
    · exact Or.inr (by have := hgap1 j h; linarith)
    · subst h; exact Or.inr (by linarith)
    · exact Or.inl (by have := hgap2 j h; linarith)
  have hsepp : ∀ i : Fin rk, rhoSq (Real.sqrt (m.toStack.coreEig i)) c < τp
      ∨ τp + δ ≤ rhoSq (Real.sqrt (m.toStack.coreEig i)) c := by
    intro i
    obtain ⟨j, rfl⟩ : ∃ j, σ j = i := ⟨σ.symm i, σ.apply_symm_apply i⟩
    rw [hρσ j, hτpdef]
    rcases lt_trichotomy j k with h | h | h
    · exact Or.inr (by have := hgap1 j h; linarith)
    · subst h; exact Or.inl (by linarith)
    · exact Or.inl (by have := hgap2 j h; linarith)
  -- step 4: the two counts
  have hcm : (Finset.univ.filter fun i : Fin rk =>
      τm < rhoSq (Real.sqrt (m.toStack.coreEig i)) c).card = (k : ℕ) + 1 := by
    rw [m.card_filter_rhoSq_coreEig hσ c τm]
    have hfilt : (Finset.univ.filter fun j : Fin rk => τm < rhoSq (m.θ j) c)
        = Finset.univ.filter fun j : Fin rk => (j : ℕ) < (k : ℕ) + 1 := by
      refine Finset.filter_congr fun j _ => ?_
      rw [hτmdef]
      constructor
      · intro h
        by_contra hcon
        have hkj : k < j := by rw [Fin.lt_def]; omega
        have := hgap2 j hkj
        linarith
      · intro h
        rcases lt_or_eq_of_le (Nat.lt_succ_iff.mp h) with h' | h'
        · have hjk : j < k := by rw [Fin.lt_def]; exact h'
          have := hgap1 j hjk
          linarith
        · have hjk : j = k := Fin.ext h'
          subst hjk
          linarith
    rw [hfilt, card_filter_lt k.isLt]
  have hcp : (Finset.univ.filter fun i : Fin rk =>
      τp < rhoSq (Real.sqrt (m.toStack.coreEig i)) c).card = (k : ℕ) := by
    rw [m.card_filter_rhoSq_coreEig hσ c τp]
    have hfilt : (Finset.univ.filter fun j : Fin rk => τp < rhoSq (m.θ j) c)
        = Finset.univ.filter fun j : Fin rk => (j : ℕ) < (k : ℕ) := by
      refine Finset.filter_congr fun j _ => ?_
      rw [hτpdef]
      constructor
      · intro h
        by_contra hcon
        rcases eq_or_lt_of_le (not_lt.mp hcon) with h' | h'
        · have hjk : j = k := Fin.ext h'.symm
          subst hjk
          linarith
        · have hkj : k < j := by rw [Fin.lt_def]; exact h'
          have := hgap2 j hkj
          linarith
      · intro h
        have hjk : j < k := by rw [Fin.lt_def]; exact h
        have := hgap1 j hjk
        linarith
    rw [hfilt, card_filter_lt (le_of_lt k.isLt)]
  -- step 5: the two count events, and their complements
  have hAm := m.toStack.tendsto_measure_count_Ioi_tau (m.gaussianNoise_toStack hG) hc hdtop
    hns hn hp hδ hτmb hsepm
  have hAp := m.toStack.tendsto_measure_count_Ioi_tau (m.gaussianNoise_toStack hG) hc hdtop
    hns hn hp hδ hτpb hsepp
  simp only [hcm] at hAm
  simp only [hcp] at hAp
  have hBm : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
      (τm < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ) + 1)})ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.toStack.measurableSet_count_Ioi N τm ((k : ℕ) + 1)).nullMeasurableSet) hAm
  have hBp : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
      (τp < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ))})ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.toStack.measurableSet_count_Ioi N τp (k : ℕ)).nullMeasurableSet) hAp
  -- step 7: the two half-line limits, and their difference
  have hg1 := m.toStack.tendstoInProb_normSq_specProj_Ioi_tau (m.gaussianNoise_toStack hG) hc
    hdtop hns hn hp (σ l) hδ hτmb hsepm
  have hg2 := m.toStack.tendstoInProb_normSq_specProj_Ioi_tau (m.gaussianNoise_toStack hG) hc
    hdtop hns hn hp (σ l) hδ hτpb hsepp
  simp only [hρσ, hβσ] at hg1 hg2
  have hval : (if τm < rhoSq (m.θ l) c then betaSq (m.θ l) c else 0)
      - (if τp < rhoSq (m.θ l) c then betaSq (m.θ l) c else 0)
      = if k = l then betaSq (m.θ k) c else 0 := by
    by_cases hkl : k = l
    · subst hkl
      rw [if_pos (show τm < rhoSq (m.θ k) c by rw [hτmdef]; linarith),
        if_neg (show ¬ τp < rhoSq (m.θ k) c by rw [hτpdef]; linarith), if_pos rfl, sub_zero]
    · rw [if_neg hkl]
      rcases lt_or_gt_of_ne (fun h : l = k => hkl h.symm) with h | h
      · have hga := hgap1 l h
        rw [if_pos (show τm < rhoSq (m.θ l) c by rw [hτmdef]; linarith),
          if_pos (show τp < rhoSq (m.θ l) c by rw [hτpdef]; linarith), sub_self]
      · have hga := hgap2 l h
        rw [if_neg (show ¬ τm < rhoSq (m.θ l) c by rw [hτmdef]; linarith),
          if_neg (show ¬ τp < rhoSq (m.θ l) c by rw [hτpdef]; linarith), sub_zero]
  have hg := hg1.sub hg2
  rw [hval] at hg
  -- step 8: the union bound
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hg
  refine tendsto_measure_zero_of_subset
    (t := fun N => ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
        (τm < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ) + 1)})ᶜ
      ∪ ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (d N))),
        (τp < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ))})ᶜ)
    ?_ (tendsto_measure_zero_union hBm hBp)
  intro N ω hω
  by_cases hlow : ∀ q : Fin (Fintype.card (Fin (d N))),
      (τm < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ) + 1)
  · by_cases hhigh : ∀ q : Fin (Fintype.card (Fin (d N))),
        (τp < (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ↔ (q : ℕ) < (k : ℕ))
    · refine absurd ?_ hω
      show overlapIdx (m.X N ω) (k : ℕ) (m.col N l)
        = ‖specProj (m.toStack.gram N ω) (Set.Ioi τm) (m.toStack.spikeVec (σ l) N)‖ ^ 2
          - ‖specProj (m.toStack.gram N ω) (Set.Ioi τp) (m.toStack.spikeVec (σ l) N)‖ ^ 2
      rw [m.normSq_specProj_spikeVec_eq hσv l N ω (Set.Ioi τm),
        m.normSq_specProj_spikeVec_eq hσv l N ω (Set.Ioi τp)]
      exact overlapIdx_eq_normSq_specProj_Ioi_sub (m.X N ω) hab hlow hhigh (m.col N l)
    · exact Set.mem_union_right _ hhigh
  · exact Set.mem_union_left _ hlow

/-- **Task C5.2, the target.** The same limit with no side condition on `n N`: the regime
gives a shift after which `rk < n (N + j)`, and `SpikedModelR.tendstoInProb_overlapIdx_of_shift`
carries the limit back. The shifted table has the same `θ`, so `hsup` transfers as it stands.
This is the supercritical half of task C10; the other half is
`SpikedModelR.align_cross_of_gaussian_subcritical` (`RankR/RMT/R6R.lean`, task C9). -/
theorem align_cross_of_gaussian_supercritical [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hreg : m.Regime c) {k : Fin rk} (hsup : c < m.θ k ^ 4) (l : Fin rk) :
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l))
      (if k = l then betaSq (m.θ k) c else 0) := by
  obtain ⟨j, pp, hn, hp⟩ := RankRStack.exists_shift_lt rk hreg.1
  have : ∀ N, IsProbabilityMeasure (μ (N + j)) := SpikedModel.isProbabilityMeasure_shift j
  exact m.tendstoInProb_overlapIdx_of_shift j (k : ℕ) l
    ((m.shift j).align_cross_of_gaussian_supercritical_aux (m.shift_GaussianNoise j hG) hc
      (RankRStack.shift_hdtop j hreg.2.1) (RankRStack.shift_hns j hreg.2.2) hn hp hsup l)

end SpikedModelR

end StackedSVD
