/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# RMT limit laws as hypothesis structures (v2)

STATUS 2026-08-29: compiled on Mathlib v4.33.0 (lake env lean, server).

`SingleTableLaw` is `prop:single_table` (both parts) in projector form, plus the two fields
the rest of the paper consumes (top eigenvalue, a.s. simplicity). Layer 1 theorems take it
as a hypothesis. Layer 2 proves `singleTableLaw_of_gaussian` (`notes/archive/rmt_roadmap.md`,
`notes/archive/SINGLE_TABLE_PLAN.md`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- Unit vectors orthogonal to `v_N` (the test directions of `prop:single_table`, part 2). -/
def orthUnit (m : SpikedModel μ n d) (N : ℕ) : Set (EuclideanSpace ℝ (Fin (d N))) :=
  {w | ‖w‖ = 1 ∧ ⟪w, m.v N⟫_ℝ = 0}

/-- `prop:single_table` and its companions.
* `align`: `overlap(X_N, v_N) → β²`.
* `delocUniform`: `sup_{w ⊥ v_N, ‖w‖ = 1} μ_N{overlap(X_N, w) ≥ ε} → 0`. This implies the
  paper's sequential form (`deloc_seq`, to prove) and is what the cross-table lemma
  `lem:delocalization` uses.
* `lamMax`: `λ_max(X_Nᵀ X_N) → ρ²(θ, c)` (used by `thm:theta_est`).
* `topSimple`: the top eigenvalue is simple a.s. (so `overlap = ⟨v̂, ·⟩²` a.s.). -/
structure SingleTableLaw (m : SpikedModel μ n d) (c : ℝ) : Prop where
  align : TendstoInProb μ (fun N ω => overlap (m.X N ω) (m.v N)) (betaSq m.θ c)
  delocUniform : ∀ ε > 0,
    Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w}) atTop (𝓝 0)
  lamMax : TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) (rhoSq m.θ c)
  topSimple : ∀ N, ∀ᵐ ω ∂(μ N),
    TopSimple ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω))

/-- The paper's sequential form of part 2, derived from `delocUniform`. -/
theorem SingleTableLaw.deloc_seq {m : SpikedModel μ n d} {c : ℝ} (law : m.SingleTableLaw c)
    (w : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))) (hw : ∀ N, w N ∈ m.orthUnit N) :
    TendstoInProb μ (fun N ω => overlap (m.X N ω) (w N)) 0 := by
  intro ε hε
  have hnn : ∀ (N : ℕ) (ω : Ω N), 0 ≤ overlap (m.X N ω) (w N) := by
    intro N ω
    unfold overlap
    positivity
  have hset : ∀ N : ℕ, {ω : Ω N | ε ≤ |overlap (m.X N ω) (w N) - 0|}
      = {ω : Ω N | ε ≤ overlap (m.X N ω) (w N)} := by
    intro N
    ext ω
    simp [abs_of_nonneg (hnn N ω)]
  simp only [hset]
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds (law.delocUniform ε hε)
    (fun N => by simp) (fun N => ?_)
  exact le_iSup₂ (f := fun v (_ : v ∈ m.orthUnit N) => μ N {ω | ε ≤ overlap (m.X N ω) v})
    (w N) (hw N)

-- The Layer 2 target `singleTableLaw_of_gaussian` now lives in `RMT/Full.lean`, which
-- imports `RMT/TailShift.lean` and `RMT/R6.lean`; this file stays upstream of both.

end SpikedModel

end StackedSVD

