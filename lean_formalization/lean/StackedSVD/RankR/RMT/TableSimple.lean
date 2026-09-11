/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.RankR.RMT.SimplicityAffineR
import StackedSVD.RankR.RMT.EdgeGlueR

/-!
# Task C7: `simple` for the Gaussian table at rank `rk`

Task C7 of `notes/archive/rankr_plan_C.md`, section "C7 and C8". This file proves the `simple` field
of `SpikedModelR.TableLawR` (`RankR/General.lean:163-168`): for Gaussian noise, almost every draw of
the table gives a Gram matrix with a simple top-`rk` spectrum. Rank-`rk` port of
`singleTableLaw_topSimple_of_gaussian` (`RMT/Simplicity.lean:374-397`).

## Route

`m.X N ω = A + t • m.Z N ω` with `A := m.U N * diagonal m.θ * (m.V N)ᵀ` and
`t := (√(d N))⁻¹` (`SpikedModelR.X`, `SpikedModelR.E`, both `def`, unfolded by defeq).
`simpleSpec_ae_affine (m.hn N) (m.hd N) A ht rk (le_min (m.rk_le_n N) (m.rk_le_d N))`
(`RankR/RMT/SimplicityAffineR.lean:76`) gives the almost-sure statement for `A + t • Z` under
`gaussianMatrix`, transported along the law of `m.Z N` by `((hG N).ae_iff hpm).mpr`
(`HasLaw.ae_iff`, `Mathlib.Probability.HasLaw`).

## Deviation from the plan

The task brief expects `hpm` to be rebuilt from `EdgeGlueR.measurableSet_simpleSpec_gram` and a
local restatement of the `private` `measurable_affine` (`RMT/Simplicity.lean`), the way the
rank-1 file does it. `RankR/RMT/EdgeGlueR.lean` already carries that exact composite at rank
`rk`, public, under `RankRStack.measurable_simpleSpec_affine`, built from the same
`measurableSet_simpleSpec_gram` (also `EdgeGlueR.lean`) and the public `measurable_affine_matrix`
(`LinAlg/SpecIdxMeas.lean`; cleanup pass 2 moved it there from `EdgeGlueR.lean`). Since this
file already imports `EdgeGlueR` and the composite lemma has the right shape with no side
condition, this file calls it directly instead of re-deriving it, which drops the local `private`
restatement the brief carried. No `sorry`, no `axiom`, no edits to any existing file.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- **Item C7.** The `simple` field of `TableLawR`, for Gaussian noise. Rank-`rk` port of
`singleTableLaw_topSimple_of_gaussian` (`RMT/Simplicity.lean:374-397`). -/
theorem simple_of_gaussian (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) (N : ℕ) :
    ∀ᵐ ω ∂(μ N),
      SimpleSpec ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω)) rk := by
  set A : Matrix (Fin (n N)) (Fin (d N)) ℝ := m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ with hA
  set t : ℝ := (Real.sqrt (d N))⁻¹ with hts
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.hd N
  have ht : t ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  have hae : ∀ᵐ Z ∂(gaussianMatrix (n N) (d N)),
      SimpleSpec ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) rk :=
    simpleSpec_ae_affine (m.hn N) (m.hd N) A ht rk (le_min (m.rk_le_n N) (m.rk_le_d N))
  have htrans : ∀ᵐ ω ∂(μ N), SimpleSpec ((A + t • m.Z N ω)ᵀ * (A + t • m.Z N ω))
      (isHermitian_transpose_mul_self (A + t • m.Z N ω)) rk :=
    ((hG N).ae_iff (RankRStack.measurable_simpleSpec_affine A t rk)).mpr hae
  filter_upwards [htrans] with ω hω
  exact hω

end SpikedModelR

end StackedSVD
