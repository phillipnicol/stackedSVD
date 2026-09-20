/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.TendstoInProb

/-!
# Chebyshev, moved out of the Gaussian tree

Two lemmas that `RMT/R2.lean` and `ThetaEst.lean` each prove once, for their own file only.
Stage 1 (the general-noise single-table law,
`notes/archive/prop_single_table_general.md`) needs both again for the `RMT/General/` files.
Those files must not import `RMT/R2.lean` or
`ThetaEst.lean` (choice 8 of the note), so this file gives one shared copy. `RMT/R2.lean` and
`ThetaEst.lean` keep their own copies unchanged, so that the Gaussian tree does not rebuild.
The dedup of the three copies into one is a later item.

`meas_ge_le_of_integral_sq` is a word-for-word copy of `R2.meas_ge_le_of_integral_sq`
(`RMT/R2.lean:795`). `tendstoInProb_of_meas_le` is a word-for-word copy of
`ThetaEst.tendstoInProb_of_meas_le` (`ThetaEst.lean:83`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped ENNReal

namespace StackedSVD
namespace Cheb

/-- Chebyshev for a nonnegative random variable with a bounded second moment. Copy of
`R2.meas_ge_le_of_integral_sq` (`RMT/R2.lean:795`), placed here so that the general-law
files (`RMT/General/`) do not import `RMT/R2.lean`. -/
theorem meas_ge_le_of_integral_sq {α : Type*} [MeasurableSpace α] (ν : Measure α)
    {F : α → ℝ} (hFm : Measurable F)
    (hint : Integrable (fun a => F a ^ 2) ν) {ε C : ℝ} (hε : 0 < ε)
    (hC : ∫ a, F a ^ 2 ∂ν ≤ C) :
    ν {a | ε ≤ F a} ≤ ENNReal.ofReal (C / ε ^ 2) := by
  have hsub : {a | ε ≤ F a} ⊆ {a | ENNReal.ofReal (ε ^ 2) ≤ ENNReal.ofReal (F a ^ 2)} := by
    intro a ha
    exact ENNReal.ofReal_le_ofReal (pow_le_pow_left₀ hε.le ha 2)
  refine (measure_mono hsub).trans ?_
  have hme : AEMeasurable (fun a => ENNReal.ofReal (F a ^ 2)) ν :=
    (ENNReal.measurable_ofReal.comp (hFm.pow_const 2)).aemeasurable
  have hε2 : ENNReal.ofReal (ε ^ 2) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    positivity
  refine (meas_ge_le_lintegral_div hme hε2 ENNReal.ofReal_ne_top).trans ?_
  have hli : ∫⁻ a, ENNReal.ofReal (F a ^ 2) ∂ν = ENNReal.ofReal (∫ a, F a ^ 2 ∂ν) :=
    (ofReal_integral_eq_lintegral_ofReal hint
      (Filter.Eventually.of_forall fun a => sq_nonneg _)).symm
  rw [hli, ← ENNReal.ofReal_div_of_pos (by positivity)]
  exact ENNReal.ofReal_le_ofReal (by gcongr)

/-- Chebyshev bounds with a vanishing constant give a limit in probability. Copy of
`ThetaEst.tendstoInProb_of_meas_le` (`ThetaEst.lean:83`). -/
theorem tendstoInProb_of_meas_le {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {f : ∀ N, Ω N → ℝ} {ctr : ℕ → ℝ} {K : ℕ → ℝ}
    (hK : Tendsto K atTop (𝓝 0))
    (hb : ∀ ε : ℝ, 0 < ε → ∀ N, μ N {ω | ε ≤ |f N ω - ctr N|}
      ≤ ENNReal.ofReal (K N / ε ^ 2)) :
    TendstoInProb μ (fun N ω => f N ω - ctr N) 0 := by
  intro ε hε
  have h1 : Tendsto (fun N => ENNReal.ofReal (K N / ε ^ 2)) atTop (𝓝 0) := by
    have h2 : Tendsto (fun N => K N / ε ^ 2) atTop (𝓝 0) := by
      simpa using hK.div_const (ε ^ 2)
    simpa using ENNReal.tendsto_ofReal h2
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds h1
    (fun _ => zero_le) fun N => ?_
  have hset : {ω : Ω N | ε ≤ |f N ω - ctr N - 0|} = {ω | ε ≤ |f N ω - ctr N|} := by
    simp
  rw [hset]
  exact hb ε hε N

end Cheb
end StackedSVD
