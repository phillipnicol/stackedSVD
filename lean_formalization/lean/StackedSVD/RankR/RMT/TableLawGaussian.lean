/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.TableAlignSup
import StackedSVD.RankR.RMT.R6R
import StackedSVD.RankR.RMT.DelocAffineR
import StackedSVD.RankR.RMT.TableSimple
import StackedSVD.RankR.GramR

/-!
# Task C10: the Gaussian `TableLawR` at general rank, proved from the noise law

Task C10 of `notes/archive/rankr_plan_C.md` section 2. This file closes Track C: under Gaussian
noise and the proportional regime, a rank-`rk` table satisfies `SpikedModelR.TableLawR c`. The black
box of Section 7 is therefore discharged, not assumed, at every `r_i` and at every spike strength.
The rank-1 mirror is `SpikedModelR.tableLawR_of_gaussian` (`RankR/GeneralFrob.lean:123`), which
stays in place as the `rk = 1` facade; the name here carries the suffix `_rk` so that a file can
import both.

The four fields come from the four Track C results:

1. `align` and `cross` split per index on `c < θ_k ^ 4`. The supercritical branch is
   `SpikedModelR.align_cross_of_gaussian_supercritical` (task C5,
   `RankR/RMT/TableAlignSup.lean`), the subcritical branch is
   `SpikedModelR.align_cross_of_gaussian_subcritical` (task C9, `RankR/RMT/R6R.lean`). At a
   subcritical index `betaSq (θ k) c` is `0` by the `else` branch of `betaSq`, so the two
   branches state the same limit.
2. `delocUniform` is `SpikedModelR.delocUniform_of_gaussian` (task C6,
   `RankR/RMT/DelocAffineR.lean`), which needs only `d → ∞`.
3. `simple` is `SpikedModelR.simple_of_gaussian` (task C7, `RankR/RMT/TableSimple.lean`).

`UnalignedModelR.tableLawR_of_gaussian_rk` lifts the table statement to every table of a
multi-table model through `gaussianNoise_of_joint_piR` (`RankR/GramR.lean`), the marginal of
the product law. It is the twin of `UnalignedModelR.tableLawR_of_gaussian`
(`RankR/GeneralFrob.lean:190`) with the rank restriction `rk = 1` removed.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- **`lem:general_rank_delocalization` under Gaussian noise** (`main_paper.tex:1948`, the
general-rank form of `prop:single_table`, `:278`): the black box `TableLawR` is discharged,
not assumed. The hypotheses are the proportional regime
`hreg` with limit `c > 0` and the Gaussian noise law `hG`; no assumption on the spikes beyond
the model fields `hθnn` and `hθanti`, and in particular no supercriticality.

Rank-1 mirror: `SpikedModelR.tableLawR_of_gaussian` (`RankR/GeneralFrob.lean:123`), which
routes through `SpikedModel.singleTableLaw_of_gaussian` instead.

Route: `align` and `cross` case split per index on `c < θ_k ^ 4` (task C5 above the
threshold, task C9 at or below it); `delocUniform` is task C6, `simple` is task C7. -/
theorem tableLawR_of_gaussian_rk [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) (hreg : m.Regime c)
    (hG : m.GaussianNoise) : m.TableLawR c := by
  refine { align := ?_, cross := ?_, delocUniform := ?_, simple := ?_ }
  · -- `align`: the limit `betaSq (θ k) c` on the supercritical branch, `0 = betaSq (θ k) c`
    -- on the other one.
    intro k
    by_cases h : c < m.θ k ^ 4
    · have hk := m.align_cross_of_gaussian_supercritical hG hc hreg h k
      rwa [if_pos rfl] at hk
    · -- `betaSq` (`Defs.lean:41`) branches on `θ ^ 4 > c`, which is `c < θ ^ 4`, so the
      -- limit is `0` here and the subcritical theorem gives it.
      rw [betaSq, if_neg h]
      exact m.align_cross_of_gaussian_subcritical hG hc hreg (not_lt.mp h) k
  · -- `cross`: the same split, with `if_neg` on the supercritical branch.
    intro k l hkl
    by_cases h : c < m.θ k ^ 4
    · have hk := m.align_cross_of_gaussian_supercritical hG hc hreg h l
      rwa [if_neg hkl] at hk
    · exact m.align_cross_of_gaussian_subcritical hG hc hreg (not_lt.mp h) l
  · -- `delocUniform`: task C6, which reads only `d → ∞` from the regime.
    exact fun k => m.delocUniform_of_gaussian hG hreg.2.1 k
  · -- `simple`: task C7.
    exact fun N => m.simple_of_gaussian hG N

end SpikedModelR

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- **Every table of a multi-table model satisfies `TableLawR` under Gaussian noise, at general
`r_i`.** The twin of `UnalignedModelR.tableLawR_of_gaussian` (`RankR/GeneralFrob.lean:190`)
without the restriction `rk = fun _ => 1`. `gaussianNoise_of_joint_piR` (`RankR/GramR.lean`)
takes the marginal of the product law, so joint Gaussian noise gives Gaussian noise on each
table. This is the input every Layer 1 theorem of Section 7 asks for. -/
theorem tableLawR_of_gaussian_rk [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∀ i, (m.tbl i).TableLawR (c i) :=
  fun i => (m.tbl i).tableLawR_of_gaussian_rk (hc i) (hreg i)
    (gaussianNoise_of_joint_piR m.tbl hG i)

end UnalignedModelR

end StackedSVD
