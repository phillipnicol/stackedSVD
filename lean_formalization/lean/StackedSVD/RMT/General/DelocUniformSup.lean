/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.DelocAlign
import StackedSVD.RMT.General.DelocUniform
import StackedSVD.RMT.General.R3minus
import StackedSVD.RMT.R5
import StackedSVD.RMT.MP7

/-!
# Item Sym (unit G7) at a general law: the `ResolventLimits` wrapper

`delocUniform_of_general` (the plan note, item Sym) is the Stage 1 target: it takes the edge
hypothesis `hedge` directly. This file builds `delocUniform_of_general'`, the same
conclusion from the seven-field interface `m.ResolventLimits c` instead, in the shape that
`RMT/General/Sup.lean`'s `singleTableLaw_of_general` already consumes for its `align` and
`lamMax` fields.

The route is one `by_cases` on `c < m.θ ^ 4`, exactly as `RMT/Full.lean` and
`RMT/General/Sup.lean`:

* Supercritical: `lamMax_tendstoInProb` (`RMT/R5.lean:711`) gives the Gram top eigenvalue's
  limit `rhoSq m.θ c`, strictly above `bulkEdge c` (`MP.bulkEdge_lt_rhoSq`,
  `RMT/MP.lean:399`); `GenRMT.mC_above_small` supplies the shrinking-window bound there from
  `MP.tendsto_mC` (`RMT/MP7.lean:567`).
* Subcritical: `lamMax_tendstoInProb_of_subcritical_general` (`RMT/General/R3minus.lean:511`)
  gives the same limit at `rhoSq m.θ c = bulkEdge c`; `GenRMT.mC_edge_small` supplies the
  window bound there directly, from the crude uniform bound `GenRMT.norm_mC_edge_le`
  (`RMT/General/DelocAlign.lean:49`).

Both branches close through `delocUniform_of_lamMax` (`RMT/General/DelocUniform.lean`, the
core of unit G7), which takes the `lamMax` limit and the window bound as hypotheses and
produces `delocUniform`. `RMT/General/Sup.lean` consumes `delocUniform_of_general'` in both
regimes; the target statement `delocUniform_of_general` of the plan note follows from it
with `H` and `RL` built from `hedge` there.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

namespace GenRMT

/-- **The window bound at the bulk edge.** `GenRMT.norm_mC_edge_le` caps `‖mC c‖` on the
segment `bulkEdge c + i η`, `0 < η ≤ 1`, by the constant `mCbound c`; taking
`η = min 1 (δ / (2 (mCbound c + 1)))` forces `η ‖mC c (bulkEdge c + i η)‖` below any `δ`. -/
theorem mC_edge_small {c : ℝ} (hc : 0 < c) : ∀ δ > 0, ∃ η : ℝ, 0 < η ∧ η ≤ 1 ∧
    η * ‖MP.mC c ((bulkEdge c : ℂ) + (η : ℂ) * Complex.I)‖ < δ := by
  intro δ hδ
  have hcb : 0 ≤ mCbound c := mCbound_nonneg hc
  have hden : (0 : ℝ) < 2 * (mCbound c + 1) := by linarith
  have hfrac : (0 : ℝ) < δ / (2 * (mCbound c + 1)) := div_pos hδ hden
  set η : ℝ := min 1 (δ / (2 * (mCbound c + 1))) with hηdef
  have hη1 : η ≤ 1 := min_le_left _ _
  have hη0 : 0 < η := lt_min one_pos hfrac
  have hηb : η ≤ δ / (2 * (mCbound c + 1)) := min_le_right _ _
  have hbound : ‖MP.mC c ((bulkEdge c : ℂ) + (η : ℂ) * Complex.I)‖ ≤ mCbound c :=
    norm_mC_edge_le hc hη0 hη1
  have hlt1 : mCbound c / (2 * (mCbound c + 1)) < 1 := by
    rw [div_lt_one hden]; linarith
  refine ⟨η, hη0, hη1, ?_⟩
  calc η * ‖MP.mC c ((bulkEdge c : ℂ) + (η : ℂ) * Complex.I)‖
      ≤ η * mCbound c := mul_le_mul_of_nonneg_left hbound hη0.le
    _ ≤ (δ / (2 * (mCbound c + 1))) * mCbound c := mul_le_mul_of_nonneg_right hηb hcb
    _ = δ * (mCbound c / (2 * (mCbound c + 1))) := by ring
    _ < δ * 1 := mul_lt_mul_of_pos_left hlt1 hδ
    _ = δ := by ring

/-- **The window bound strictly above the edge.** `η ‖mC c (x + i η)‖ → 0 ‖m c x‖ = 0` along
`η ↓ 0` (`MP.tendsto_mC`), so for every `δ` some `η ∈ (0, 1]` makes it small. -/
theorem mC_above_small {c : ℝ} (hc : 0 < c) {x : ℝ} (hx : bulkEdge c < x) : ∀ δ > 0, ∃ η :
    ℝ, 0 < η ∧ η ≤ 1 ∧ η * ‖MP.mC c ((x : ℂ) + (η : ℂ) * Complex.I)‖ < δ := by
  intro δ hδ
  have htend : Tendsto (fun η : ℝ => η * ‖MP.mC c ((x : ℂ) + (η : ℂ) * Complex.I)‖)
      (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) := by
    have h1 : Tendsto (fun η : ℝ => MP.mC c ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
        (𝓝 ((MP.m c x : ℝ) : ℂ)) := MP.tendsto_mC hc hx
    have h3 : Tendsto (fun η : ℝ => η) (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
      tendsto_nhdsWithin_of_tendsto_nhds tendsto_id
    simpa using h3.mul h1.norm
  rw [Metric.tendsto_nhdsWithin_nhds] at htend
  obtain ⟨δ', hδ', hwit⟩ := htend δ hδ
  set η : ℝ := min (δ' / 2) 1 with hηdef
  have hη0 : 0 < η := lt_min (by linarith) one_pos
  have hη1 : η ≤ 1 := min_le_right _ _
  have hηδ' : η < δ' := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  have hmem : η ∈ Set.Ioi (0 : ℝ) := hη0
  have hdist : dist η (0 : ℝ) < δ' := by
    rw [Real.dist_eq, sub_zero, abs_of_pos hη0]; exact hηδ'
  have hfin := hwit hmem hdist
  rw [Real.dist_eq, sub_zero] at hfin
  exact ⟨η, hη0, hη1, (abs_lt.mp hfin).2⟩

end GenRMT

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **Item Sym at a general law, `ResolventLimits` form.** Same conclusion as
`delocUniform_of_general`, from the seven-field interface instead of the bare edge. Both
regimes, from `delocUniform_of_lamMax` through the window bounds `GenRMT.mC_above_small`
(supercritical) and `GenRMT.mC_edge_small` (subcritical, after `rhoSq m.θ c = bulkEdge c`). -/
theorem delocUniform_of_general' [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (H : m.ResolventFormsC c) (RL : m.ResolventLimits c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) :
    ∀ ε > 0, Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w})
      atTop (𝓝 0) := by
  by_cases hθ : c < m.θ ^ 4
  · have hlam := lamMax_tendstoInProb RL hc hθ
    have hθ0 : 0 < m.θ := by
      rcases m.hθ.lt_or_eq with h | h
      · exact h
      · exfalso
        have h4 : m.θ ^ 4 = 0 := by rw [← h]; norm_num
        linarith
    have hbρ : bulkEdge c < rhoSq m.θ c := MP.bulkEdge_lt_rhoSq hc hθ0 hθ
    exact delocUniform_of_lamMax hc m H hreg hν hG hlam (GenRMT.mC_above_small hc hbρ)
  · have hθ' : m.θ ^ 4 ≤ c := not_lt.mp hθ
    have hlam := lamMax_tendstoInProb_of_subcritical_general RL hc hreg hν hG hθ'
    have hrho : rhoSq m.θ c = bulkEdge c := by rw [rhoSq, if_neg (not_lt.mpr hθ')]
    rw [hrho] at hlam
    exact delocUniform_of_lamMax hc m H hreg hν hG hlam (GenRMT.mC_edge_small hc)

end SpikedModel

end StackedSVD
