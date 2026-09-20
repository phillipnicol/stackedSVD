/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.RankR.RMT.EdgeGlueR
import StackedSVD.RankR.RMT.EdgeR

/-!
# Task C9: the subcritical index of one rank-`rk` table

Task C9 of `notes/archive/rankr_plan_C.md`, in the form the audit of 2026-09-02 gives it
(`notes/archive/audit_rankr_plan_C_2026-09-02.md`, attack 1). At a **subcritical** spike index
`k` (`θ_k^4 ≤ c`) the `k`-th singular subspace of `X` has no overlap with any spike
direction of the table: the limit is `0` for every `l`, so this one statement covers the
`align` field of `SpikedModelR.TableLawR` at a subcritical index (`betaSq θ c = 0` there)
and the `cross` field at such an index. The paper statements are `prop:single_table` (rank
1) and `lem:general_rank_delocalization` (`main_paper.tex:1948`, general rank); the
rank-1 mirror is `SpikedModel.align_tendstoInProb_of_subcritical` (**R6'**, `RMT/R6.lean`).

## The route

The plan's original route (a decomposition along the column span of `Q`, a port of
`RMT/R6.lean`) is not used, and its step 2 is unsound: at a subcritical index the
eigenvalue sits at or below the threshold, so `specProj (Set.Ioi τ)` does not dominate
`specProjIdx k`. The proved rank-`r` template is task U7c-edge
(`RankRStack.tendsto_measure_normSq_specProj_edge_gt`, `RankR/RMT/EdgeGlueR.lean`), which
bounds the projection of a spike direction on `topEigSet r ∩ Set.Iic (bulkEdge c + ε₁)`
with **no** supercriticality. This file is the assembly of `RankR/RMT/AlignG.lean` with two
bad events instead of three and the target `0`:

1. `normSq_specProjIdx_le_edge`: the deterministic domination. When the sorted eigenvalue at
   the index `k` is at most `τ` and `k < r`, the index projector is dominated by the edge
   projector, through `EdgeGlueDetR.normSq_specProj_mono`.
2. `SpikedModelR.exists_split_subcritical`: `hθanti` makes every model index at or above a
   subcritical `k` subcritical too, so the spikes split as
   `Fin (rk - k) ⊕ Fin k ≃ Fin rk` with the subcritical block on the left. The count event
   `RankRStack.tendsto_measure_eigenvalues₀_gt` (task U7a, `RankR/RMT/EdgeR.lean`) then puts
   every sorted eigenvalue of index `k` or above at or below `bulkEdge c + ε₁`.
3. `SpikedModelR.align_cross_of_gaussian_subcritical_aux`: the union of the two bad events,
   with `C1.4` (`normSq_specProj_spikeVec_eq`) to move `col N l` to `± spikeVec (σ l)`.
4. `SpikedModelR.align_cross_of_gaussian_subcritical`: the tail shift of
   `RankR/RMT/ShiftR.lean` and `RankR/RMT/TableStack.lean` removes the side condition
   `∀ N, n N = rk + pp N`.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 1. The deterministic domination -/

/-- **The domination step of task C9** (audit of `notes/archive/rankr_plan_C.md`, attack 1). The
eigenvalue set at the sorted index `k` is a subset of `topEigSet r ∩ Set.Iic τ` as soon as
`k < r` and the sorted eigenvalue at that index is at most `τ`, so the index projector is
dominated by the edge projector U7c-edge bounds. -/
theorem normSq_specProjIdx_le_edge {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}
    (hS : S.IsHermitian) {r k : ℕ} (hk : k < r) {τ : ℝ}
    (hle : ∀ q : Fin (Fintype.card (Fin p)), (q : ℕ) = k → hS.eigenvalues₀ q ≤ τ)
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx S hS k x‖ ^ 2
      ≤ ‖specProj S (topEigSet S hS r ∩ Set.Iic τ) x‖ ^ 2 := by
  rw [specProjIdx]
  refine EdgeGlueDetR.normSq_specProj_mono hS ?_ x
  rintro t ⟨q, hq, rfl⟩
  exact ⟨⟨q, hq ▸ hk, rfl⟩, hle q hq⟩

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-! ### 2. The split at a subcritical model index -/

/-- **The split.** `hθanti` is strict, so every model index at or above a subcritical index
`k` is subcritical. The spikes therefore split as `Fin (rk - k) ⊕ Fin k ≃ Fin rk` with the
subcritical block on the left, which is the shape
`RankRStack.tendsto_measure_eigenvalues₀_le` and its complement form take. The right block
has exactly `k` members, so the count event reaches the index `k` itself. -/
theorem exists_split_subcritical (m : SpikedModelR μ n d rk) {c : ℝ}
    {σ : Equiv.Perm (Fin rk)} (hσ : ∀ j, m.toStack.coreEig (σ j) = m.θ j ^ 2)
    {k : Fin rk} (hsub : m.θ k ^ 4 ≤ c) :
    ∃ e : Fin (rk - (k : ℕ)) ⊕ Fin (k : ℕ) ≃ Fin rk,
      ∀ a : Fin (rk - (k : ℕ)), m.toStack.coreEig (e (Sum.inl a)) ^ 2 ≤ c := by
  have hkk : (k : ℕ) + (rk - (k : ℕ)) = rk := by
    have := k.isLt
    omega
  refine ⟨((Equiv.sumComm (Fin (rk - (k : ℕ))) (Fin (k : ℕ))).trans finSumFinEquiv).trans
    ((finCongr hkk).trans σ), fun a => ?_⟩
  set i : Fin rk := finCongr hkk (Fin.natAdd (k : ℕ) a) with hidef
  have hidx : (((Equiv.sumComm (Fin (rk - (k : ℕ))) (Fin (k : ℕ))).trans finSumFinEquiv).trans
      ((finCongr hkk).trans σ)) (Sum.inl a) = σ i := rfl
  have hival : (i : ℕ) = (k : ℕ) + (a : ℕ) := rfl
  have hki : k ≤ i := Fin.le_def.mpr (by rw [hival]; omega)
  have hθle : m.θ i ≤ m.θ k := m.hθanti.antitone hki
  rw [hidx, hσ i]
  have h4 : m.θ i ^ 4 ≤ m.θ k ^ 4 := pow_le_pow_left₀ (m.hθnn i) hθle 4
  calc (m.θ i ^ 2) ^ 2 = m.θ i ^ 4 := by ring
    _ ≤ m.θ k ^ 4 := h4
    _ ≤ c := hsub

/-! ### 3. The limit, with the side condition of the Gaussian chain -/

/-- **Task C9 with the side condition.** At a subcritical index `k` the overlap of the
`k`-th singular subspace of `X` with any spike direction `col N l` tends to `0` in
probability. The `_aux` suffix marks the side condition `hn`/`hp` (that is, `rk < n N` at
every `N`), which `align_cross_of_gaussian_subcritical` removes by the tail shift.

The two bad events are the count event of task U7a
(`RankRStack.tendsto_measure_eigenvalues₀_gt`) and the edge event of task U7c-edge
(`RankRStack.tendsto_measure_normSq_specProj_edge_gt`). No spike has to be supercritical
and no margin between the spikes is assumed. -/
theorem align_cross_of_gaussian_subcritical_aux [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, n N = rk + pp N) (hp : ∀ N, 0 < pp N)
    {k : Fin rk} (hsub : m.θ k ^ 4 ≤ c) (l : Fin rk) :
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l)) 0 := by
  classical
  obtain ⟨σ, hσeig, hσvec⟩ := m.exists_perm_coreEig
  obtain ⟨U, hU, hsig, hgram, h⟩ :=
    m.toStack.resolventLimitsR_of_gaussian (m.gaussianNoise_toStack hG) hc hdtop hns hn hp
  obtain ⟨e, hesub⟩ := m.exists_split_subcritical hσeig hsub
  intro ε hε
  obtain ⟨ε₀, hε₀, hedge⟩ := m.toStack.tendsto_measure_normSq_specProj_edge_gt
    (m.gaussianNoise_toStack hG) hc hdtop hns hn hp (σ l) (ε / 2) (by linarith)
  have hBcnt := m.toStack.tendsto_measure_eigenvalues₀_gt hc h hgram e hesub ε₀ hε₀
  have hBedge := hedge ε₀ hε₀ le_rfl
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | ∀ q : Fin (Fintype.card (Fin (d N))), (k : ℕ) ≤ (q : ℕ) →
          (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ≤ bulkEdge c + ε₀}ᶜ
        ∪ {ω | ε / 2 < ‖specProj (m.toStack.gram N ω)
            (topEigSet (m.toStack.gram N ω) (m.toStack.isHermitian_gram N ω) rk
              ∩ Set.Iic (bulkEdge c + ε₀)) (m.toStack.spikeVec (σ l) N)‖ ^ 2})
    ?_ (tendsto_measure_zero_union hBcnt hBedge)
  intro N ω hω
  by_cases hcnt : ∀ q : Fin (Fintype.card (Fin (d N))), (k : ℕ) ≤ (q : ℕ) →
      (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ≤ bulkEdge c + ε₀
  · -- on the count event the index projector is dominated by the edge projector
    refine Set.mem_union_right _ ?_
    have hle : ∀ q : Fin (Fintype.card (Fin (d N))), (q : ℕ) = (k : ℕ) →
        (m.toStack.isHermitian_gram N ω).eigenvalues₀ q ≤ bulkEdge c + ε₀ :=
      fun q hq => hcnt q (le_of_eq hq.symm)
    have hdom : overlapIdx (m.X N ω) (k : ℕ) (m.col N l)
        ≤ ‖specProj (m.toStack.gram N ω)
            (topEigSet (m.toStack.gram N ω) (m.toStack.isHermitian_gram N ω) rk
              ∩ Set.Iic (bulkEdge c + ε₀)) (m.toStack.spikeVec (σ l) N)‖ ^ 2 := by
      rw [m.normSq_specProj_spikeVec_eq hσvec l N ω _,
        m.overlapIdx_toStack N ω (k : ℕ) (m.col N l)]
      exact normSq_specProjIdx_le_edge (m.toStack.isHermitian_gram N ω) k.isLt hle _
    have hnn : 0 ≤ overlapIdx (m.X N ω) (k : ℕ) (m.col N l) := by
      rw [overlapIdx]
      positivity
    have hωv : ε ≤ |overlapIdx (m.X N ω) (k : ℕ) (m.col N l) - 0| := hω
    rw [sub_zero, abs_of_nonneg hnn] at hωv
    exact lt_of_lt_of_le (by linarith : ε / 2 < overlapIdx (m.X N ω) (k : ℕ) (m.col N l)) hdom
  · exact Set.mem_union_left _ hcnt

/-! ### 4. The limit -/

/-- **Task C9.** At a subcritical index `k` (`θ_k^4 ≤ c`) the `k`-th singular subspace of a
Gaussian rank-`rk` table has vanishing overlap with every spike direction of the table. This
is the `align` field of `SpikedModelR.TableLawR` at a subcritical index, because
`betaSq (m.θ k) c = 0` there, and the `cross` field at such an index. The rank-1 mirror is
`SpikedModel.align_tendstoInProb_of_subcritical` (**R6'**, `RMT/R6.lean`).

The only hypotheses are the Gaussian noise law, the regime and `θ_k^4 ≤ c`, the last with
equality allowed. Task C10 pairs this with the supercritical facade of task C8. -/
theorem align_cross_of_gaussian_subcritical [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hreg : m.Regime c) {k : Fin rk} (hsub : m.θ k ^ 4 ≤ c) (l : Fin rk) :
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l)) 0 := by
  obtain ⟨j, pp, hn, hp⟩ := RankRStack.exists_shift_lt rk hreg.1
  have : ∀ N, IsProbabilityMeasure (μ (N + j)) := SpikedModel.isProbabilityMeasure_shift j
  exact m.tendstoInProb_overlapIdx_of_shift j (k : ℕ) l
    ((m.shift j).align_cross_of_gaussian_subcritical_aux (m.shift_GaussianNoise j hG) hc
      (RankRStack.shift_hdtop j hreg.2.1) (RankRStack.shift_hns j hreg.2.2) hn hp hsub l)

end SpikedModelR

end StackedSVD
