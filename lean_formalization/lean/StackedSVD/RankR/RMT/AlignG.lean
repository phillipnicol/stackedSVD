/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.AlignOutG
import StackedSVD.RankR.RMT.EdgeGlueR

/-!
# Task A2b: the rank-`r` spiked law of a Gaussian stack, every spike, mixed regime

The rank-1 mirror is `RankRStack.align_of_gaussian_supercritical`
(`RankR/RMT/Outliers.lean`), which needs every spike supercritical. Here the spikes may be
subcritical as well, and the limit `betaSq (√λ_j) c` is `0` at a subcritical spike.

Two earlier tasks supply the two halves of the projector. Task A2a
(`RankRStack.tendstoInProb_normSq_specProj_Ioi`, `RankR/RMT/AlignOutG.lean`) gives the part
above the threshold `τ = bulkEdge c + ε₁`; task U7c
(`RankRStack.tendsto_measure_normSq_specProj_edge_gt`, `RankR/RMT/EdgeGlueR.lean`) gives the
edge part, the top-`r` eigenvalues at most `τ`. Task U7a
(`RankRStack.tendsto_measure_eigenvalues₀_gt`, through `RankR/RMT/EdgeGlueR.lean`) says that
every sorted eigenvalue of index at least `r` is at most `τ`, which is the hypothesis of the
split `Frame.norm_sq_specProjTop_split`.

Contents:

1. `RankRStack.align_of_gaussian_aux`: the shared-interface statement.
2. `UnalignedModel.subspaceLaw_of_gaussian_aux`: the `r_i = 1` facade,
   `prop:single_table` in rank `r`.

The `_aux` suffix marks the side condition `hn`/`hp` (that is, `r < ns N` at every `N`).
Task A3 removes it by the tail shift and owns the unsuffixed names.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- **Task A2b.** For a `RankRStack` with Gaussian noise and aspect ratio `c`, the top-`r`
eigenspace of the Gram matrix overlaps each spike direction `V q_j` by `betaSq (√λ_j) c` in
probability. No spike has to be supercritical: this is
`RankRStack.align_of_gaussian_supercritical` (`RankR/RMT/Outliers.lean`) without `hsup`.

Ties among the `λ_j(C)` are allowed. The side condition `hn`/`hp` is task A3. -/
theorem align_of_gaussian_aux [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) :
    ∀ j : Fin r, TendstoInProb μ
      (fun N ω => ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
        (s.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (s.coreEig j)) c) := by
  classical
  obtain ⟨U, hU, hsig, hgram, h⟩ := s.resolventLimitsR_of_gaussian hG hc hdtop hns hn hp
  intro j ε hε
  have hδ : (0 : ℝ) < ε / 2 := by linarith
  obtain ⟨ε₀, hε₀, hedge⟩ :=
    s.tendsto_measure_normSq_specProj_edge_gt hG hc hdtop hns hn hp j (ε / 2) hδ
  obtain ⟨ε₀', hε₀', hmar⟩ := s.exists_margin_le_rhoSq hc
  set ε₁ : ℝ := min ε₀ ε₀' with hε₁def
  have hε₁ : 0 < ε₁ := lt_min hε₀ hε₀'
  have hε₁a : ε₁ ≤ ε₀ := min_le_left _ _
  have hε₁b : ε₁ ≤ ε₀' := min_le_right _ _
  -- the three bad families
  have hBcnt : Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (d N))),
      r ≤ (k : ℕ) → (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε₁}ᶜ)
      atTop (𝓝 0) :=
    s.tendsto_measure_eigenvalues₀_gt hc h hgram (Equiv.emptySum (Fin 0) (Fin r))
      (fun a => a.elim0) ε₁ hε₁
  have hBedge : Tendsto (fun N => μ N {ω | ε / 2 < ‖specProj (s.gram N ω)
      (topEigSet (s.gram N ω) (s.isHermitian_gram N ω) r ∩ Set.Iic (bulkEdge c + ε₁))
      (s.spikeVec j N)‖ ^ 2}) atTop (𝓝 0) := hedge ε₁ hε₁ hε₁a
  have hBout : Tendsto (fun N => μ N {ω | ε / 2 ≤ |‖specProj (s.gram N ω)
      (Set.Ioi (bulkEdge c + ε₁)) (s.spikeVec j N)‖ ^ 2
      - betaSq (Real.sqrt (s.coreEig j)) c|}) atTop (𝓝 0) :=
    s.tendstoInProb_normSq_specProj_Ioi hG hc hdtop hns hn hp j hε₁
      (fun k hk => hmar ε₁ hε₁ hε₁b k hk) (ε / 2) hδ
  -- the bad set sits inside the union of the three
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | ∀ k : Fin (Fintype.card (Fin (d N))),
          r ≤ (k : ℕ) → (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε₁}ᶜ
        ∪ ({ω | ε / 2 < ‖specProj (s.gram N ω)
              (topEigSet (s.gram N ω) (s.isHermitian_gram N ω) r
                ∩ Set.Iic (bulkEdge c + ε₁)) (s.spikeVec j N)‖ ^ 2}
          ∪ {ω | ε / 2 ≤ |‖specProj (s.gram N ω) (Set.Ioi (bulkEdge c + ε₁))
              (s.spikeVec j N)‖ ^ 2 - betaSq (Real.sqrt (s.coreEig j)) c|})) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_ofPred_eq,
      not_or, not_not, not_lt, not_le] at hbad
    obtain ⟨hcnt, hbe, hbo⟩ := hbad
    have hω' : ε ≤ |‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
        (s.spikeVec j N)‖ ^ 2 - betaSq (Real.sqrt (s.coreEig j)) c| := hω
    rw [Frame.norm_sq_specProjTop_split (s.isHermitian_gram N ω) hcnt (s.spikeVec j N)] at hω'
    set a : ℝ := ‖specProj (s.gram N ω) (Set.Ioi (bulkEdge c + ε₁)) (s.spikeVec j N)‖ ^ 2
      with hadef
    set b : ℝ := ‖specProj (s.gram N ω)
        (topEigSet (s.gram N ω) (s.isHermitian_gram N ω) r ∩ Set.Iic (bulkEdge c + ε₁))
        (s.spikeVec j N)‖ ^ 2 with hbdef
    have hb0 : (0 : ℝ) ≤ b := by rw [hbdef]; positivity
    have habs : |a + b - betaSq (Real.sqrt (s.coreEig j)) c|
        ≤ |a - betaSq (Real.sqrt (s.coreEig j)) c| + b := by
      have heq : a + b - betaSq (Real.sqrt (s.coreEig j)) c
          = (a - betaSq (Real.sqrt (s.coreEig j)) c) + b := by ring
      rw [heq]
      calc |(a - betaSq (Real.sqrt (s.coreEig j)) c) + b|
          ≤ |a - betaSq (Real.sqrt (s.coreEig j)) c| + |b| := abs_add_le _ _
        _ = |a - betaSq (Real.sqrt (s.coreEig j)) c| + b := by rw [abs_of_nonneg hb0]
    linarith
  · exact tendsto_measure_zero_union hBcnt (tendsto_measure_zero_union hBedge hBout)

end RankRStack

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **`prop:single_table` in rank `r`** for an `UnalignedModel` with joint Gaussian noise, in
the mixed regime. This is `UnalignedModel.subspaceLaw_of_gaussian_supercritical`
(`RankR/RMT/Outliers.lean`) without `hsup`, so it is the exact hypothesis of
`prop_stacksvd_subspace` and the two compose with no glue.

Ties among the `λ_j(C)` are allowed. The side condition `hn`/`hp` is task A3. -/
theorem subspaceLaw_of_gaussian_aux [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i)
    {pp : ℕ → ℕ} (hn : ∀ N, ∑ i, n i N = r + pp N) (hp : ∀ N, 0 < pp N) :
    m.SubspaceLaw (∑ i, cc i) :=
  ⟨m.toStack.align_of_gaussian_aux (m.gaussianNoise_toStack hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTable.stack_regime cc hreg).2.2 hn hp⟩

end UnalignedModel

end StackedSVD
