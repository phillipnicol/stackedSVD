/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SubspaceG
import StackedSVD.RankR.RMT.AlignG

/-!
# Task A2b, part 2: the general-`r_i` stack as `RankRStack`, and `prop:single_table` on it

`RankR/SubspaceG.lean` (task PA) built every object a `RankRStack`
(`RankR/RMT/Stack.lean`) asks of the unit-weight stack of an `UnalignedModelR`, except two:
the law of the stacked noise and its measurability. This file adds those two, assembles
`UnalignedModelR.toStackG`, and reads task A2b (`RankRStack.align_of_gaussian_aux`) through
it. The rank-1 mirror is `UnalignedModel.toStack` (`RankR/RMT/Split.lean`) with
`UnalignedModel.subspaceLaw_of_gaussian_aux` (`RankR/RMT/AlignG.lean`).

The law of `stackZG` comes from `MultiTableModel.stack_law`, which is stated for a
`MultiTableModel` of rank-one tables. `tblShell` is the rank-one shell of table `i`: the same
noise, a zero spike and coordinate singular vectors. Only its noise field is read, so the
shell carries the law with no extra work (decision recorded in `notes/FLAGGED.md` item 8).

Contents:

1. `UnalignedModelR.tblShell`, `toMultiTableShell`, `stackZG_eq_shell`: the shell.
2. `UnalignedModelR.hasLaw_stackZG`, `measurable_stackZG`: the two missing data.
3. `UnalignedModelR.toStackG` and its `rfl` bridges.
4. `UnalignedModelR.subspaceLawG_of_gaussian_aux`: `prop:single_table` in rank `r` at general
   `r_i`, in the mixed regime.

The `_aux` suffix marks the side condition `hn`/`hp` (that is, `r < ∑ i, n i N` at every `N`).
Task A3 removes it.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. The rank-one shell -/

/-- A rank-one shell of table `i`: the same noise, a zero spike, unit coordinate vectors as
the singular vectors. Only its noise is read; it exists so that `MultiTableModel.stack_law`
gives the law of `stackZG`. -/
noncomputable def tblShell (m : UnalignedModelR μ M n d r rk) (i : Fin M) :
    SpikedModel μ (n i) d where
  θ := 0
  u := fun N => EuclideanSpace.single ⟨0, (m.tbl i).hn N⟩ 1
  v := fun N => EuclideanSpace.single ⟨0, (m.tbl i).hd N⟩ 1
  Z := (m.tbl i).Z
  hθ := le_rfl
  hn := (m.tbl i).hn
  hd := (m.tbl i).hd
  hu := fun N => by simp
  hv := fun N => by simp
  hZ := (m.tbl i).hZ

/-- The shell as a `MultiTableModel`; the shared right singular vector is the same coordinate
vector in every table. -/
noncomputable def toMultiTableShell (m : UnalignedModelR μ M n d r rk) :
    MultiTableModel μ M n d where
  tbl := m.tblShell
  hv := fun _ _ _ => rfl

/-- The stacked unscaled noise of the model is the stacked noise of the shell. -/
theorem stackZG_eq_shell (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.stackZG N ω = m.toMultiTableShell.stackZ N ω := rfl

/-! ### 2. The law and the measurability of the stacked noise -/

/-- **The law of the stacked noise.** `MultiTableModel.stack_law` on the shell. `NeZero M` is
the side condition of `MultiTableModel.stack`, which reads table `0`. -/
theorem hasLaw_stackZG [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) (N : ℕ) :
    HasLaw (m.stackZG N) (gaussianMatrix (∑ i, n i N) (d N)) (μ N) :=
  m.toMultiTableShell.stack_law hG N

/-- The stacked unscaled noise is measurable, entry by entry from `SpikedModelR.hZ`. The
rank-one mirror is `UnalignedModel.measurable_stackZu` (`RankR/RMT/Split.lean`). -/
theorem measurable_stackZG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Measurable (m.stackZG N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => m.stackZG N ω q l)
      = fun ω => (m.tbl (finSigmaFinEquiv.symm q).1).Z N ω
          (finSigmaFinEquiv.symm q).2 l := rfl
  rw [h]
  exact (measurable_pi_apply l).comp
    ((measurable_pi_apply _).comp ((m.tbl (finSigmaFinEquiv.symm q).1).hZ N))

/-! ### 3. The stack as `RankRStack` -/

/-- **The general-`r_i` stack as `RankRStack` data.** Every field is a field or a proved
identity of `RankR/SubspaceG.lean`, so each derived object of `RankRStack` is definitionally
the matching object of `UnalignedModelR` and the bridge lemmas below are `rfl`. `NeZero M`
supplies `0 < d N`, which the model carries only inside a table. -/
noncomputable def toStackG [NeZero M] (m : UnalignedModelR μ M n d r rk) :
    RankRStack μ (fun N => ∑ i, n i N) d r where
  V := m.V
  A := m.signalFactorG
  core := m.coreG
  Zu := m.stackZG
  E := m.stackEG
  X := m.stackXG
  hV := m.hV
  hcore := m.signalFactorG_transpose_mul_self
  hE := m.stackE_eqG
  hX := m.stackX_eqG
  hd := fun N => (m.tbl ⟨0, NeZero.pos M⟩).hd N
  hZmeas := m.measurable_stackZG

/-- Joint Gaussian noise of the tables gives Gaussian noise of the stack data. -/
theorem gaussianNoise_toStackG [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) : m.toStackG.GaussianNoise :=
  fun N => m.hasLaw_stackZG hG N

theorem toStackG_gram [NeZero M] (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.toStackG.gram N ω = m.stackGramG N ω := rfl

theorem toStackG_core [NeZero M] (m : UnalignedModelR μ M n d r rk) :
    m.toStackG.core = m.coreG := rfl

theorem toStackG_coreEig [NeZero M] (m : UnalignedModelR μ M n d r rk) (j : Fin r) :
    m.toStackG.coreEig j = m.coreEigG j := rfl

theorem toStackG_spikeVec [NeZero M] (m : UnalignedModelR μ M n d r rk) (j : Fin r)
    (N : ℕ) : m.toStackG.spikeVec j N = m.spikeVecG j N := rfl

/-! ### 4. `prop:single_table` in rank `r` at general `r_i` -/

/-- **`prop:single_table` in rank `r`** for the general-`r_i` model with joint Gaussian noise,
in the mixed regime. It is the exact hypothesis of
`UnalignedModelR.prop_stacksvd_subspace_general` (`RankR/SubspaceG.lean`), so the two compose
with no glue. The rank-1 mirror is `UnalignedModel.subspaceLaw_of_gaussian_aux`
(`RankR/RMT/AlignG.lean`).

Ties among the `λ_j(C)` are allowed, and no spike has to be supercritical. The side condition
`hn`/`hp` is task A3. -/
theorem subspaceLawG_of_gaussian_aux [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i)
    {pp : ℕ → ℕ} (hn : ∀ N, ∑ i, n i N = r + pp N) (hp : ∀ N, 0 < pp N) :
    m.SubspaceLawG (∑ i, cc i) :=
  ⟨m.toStackG.align_of_gaussian_aux (m.gaussianNoise_toStackG hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTableShell.stack_regime cc hreg).2.2 hn hp⟩

end UnalignedModelR

end StackedSVD
