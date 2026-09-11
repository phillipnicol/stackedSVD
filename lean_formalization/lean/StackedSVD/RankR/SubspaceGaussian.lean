/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SubspaceGStack
import StackedSVD.RankR.RMT.ShiftR
import StackedSVD.RankR.Frobenius
import StackedSVD.RankR.RMT.SimplicityAffineR
import StackedSVD.RankR.SubspaceMain

/-!
# Task A3b: `prop:stacksvd_subspace` for Gaussian noise, with no side condition

Task A2b proves the rank-`r` alignment of a Gaussian stack under the extra hypothesis
`ns N = r + pp N`, `0 < pp N` (the `_aux` names of `RankR/RMT/AlignG.lean` and
`RankR/SubspaceGStack.lean`). Task A3a builds the tail shift `RankRStack.shift`
(`RankR/RMT/ShiftR.lean`), which produces that side condition from the regime limit alone.
This file joins the two, so every statement here takes the regime and Gaussian noise and
nothing else.

## Content

1. `RankRStack.align_of_gaussian`: the `_aux` statement without `hn` and `hp`.
2. `UnalignedModel.subspaceLaw_of_gaussian` and `UnalignedModelR.subspaceLawG_of_gaussian`:
   `prop:single_table` in rank `r`, at `r_i = 1` and at general `r_i`.
3. `UnalignedModel.prop_stacksvd_subspace_gaussian` and
   `UnalignedModelR.prop_stacksvd_subspace_general_gaussian`: `prop:stacksvd_subspace`
   (`main_paper.tex:834`) with the law discharged.
4. The Frobenius form of the same proposition at the canonical eigenframe:
   `UnalignedModel.r_le_d`, `UnalignedModel.topGap_stackGram_whp_of_gaussian` and
   `UnalignedModel.prop_stacksvd_subspace_frobenius_eig_gaussian`. The top-`r` gap of the
   stack Gram matrix holds almost surely at every large `N`, by `simpleSpec_ae_affine`
   (`RankR/RMT/SimplicityAffineR.lean`) transported through the law of the stacked noise.
4b. The general-`r_i` twin of section 4, all six declarations in this file (unlike section 4,
   whose two Layer 1 corollaries live in `RankR/Frobenius.lean`): `UnalignedModelR.r_le_d`,
   `UnalignedModelR.prop_stacksvd_subspace_general_frobenius`,
   `UnalignedModelR.prop_stacksvd_subspace_general_frobenius_eig`,
   `UnalignedModelR.ae_topGap_stackGramG`, `UnalignedModelR.topGap_stackGramG_whp_of_gaussian`
   and `UnalignedModelR.prop_stacksvd_subspace_general_frobenius_eig_gaussian`, each the same
   proof as its section 4 counterpart with the `G`-suffixed stack objects.

5. `RankR.Example.example_perfStackR_tendsto_gaussian`: the stacksvd half of the worked
   example of Section 7 (`eq:perf_rankr_ex_stacksvd`), with its law hypothesis discharged.

The rank-1 mirror of section 3 is `MultiTableModel.prop_stacksvd_general_gaussian`
(`StackSVD/Main.lean:91`).

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 0. A simple top-`rk` spectrum gives the top-`rk` gap -/

/-- `SimpleSpec A hA rk` gives `TopGap A hA rk`. This is the last four lines of
`topGap_ae_affine` (`RankR/RMT/SimplicityAffineR.lean:98`) read as a deterministic step: the
sorted eigenvalues are antitone, so a strict inequality follows from `≤` and `≠`. -/
theorem topGap_of_simpleSpec {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {rk : ℕ} (h : SimpleSpec A hA rk) : TopGap A hA rk := by
  intro k l hk hl
  have hval : (k : ℕ) < (l : ℕ) := lt_of_lt_of_le hk hl
  have hne : k ≠ l := by
    intro hkl
    rw [hkl] at hval
    exact lt_irrefl _ hval
  exact lt_of_le_of_ne (hA.eigenvalues₀_antitone (Fin.le_def.mpr hval.le)) (h k l hk hne)

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 1. The alignment of every spike, with no side condition -/

/-- **Task A3b, the unconditional alignment.** For a `RankRStack` with Gaussian noise and
aspect ratio `c`, the top-`r` eigenspace of the Gram matrix overlaps the spike direction
`V q_j` by `betaSq (√λ_j) c` in probability. No spike has to be supercritical, ties among the
`λ_j(C)` are allowed, and no lower bound on `ns N` is assumed.

This is `RankRStack.align_of_gaussian_aux` (`RankR/RMT/AlignG.lean:54`) with `hn` and `hp`
discharged by the tail shift of `RankR/RMT/ShiftR.lean`. The rank-1 mirror is
`SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean:49`). -/
theorem align_of_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)) (j : Fin r) :
    TendstoInProb μ (fun N ω => ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
      (s.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (s.coreEig j)) c) := by
  obtain ⟨k, pp, hn, hp⟩ := exists_shift_lt r (tendsto_ns_atTop hc hdtop hns)
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := SpikedModel.isProbabilityMeasure_shift k
  exact s.tendstoInProb_normSq_specProjTop_of_shift k j
    ((s.shift k).align_of_gaussian_aux (s.shift_GaussianNoise k hG) hc
      (shift_hdtop k hdtop) (shift_hns k hns) hn hp j)

end RankRStack

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! ### 2. `prop:single_table` in rank `r`, at `r_i = 1` -/

/-- **`prop:single_table` in rank `r`** for an `UnalignedModel` with joint Gaussian noise, in
the mixed regime, with no side condition on `∑_i n_i`. This is
`UnalignedModel.subspaceLaw_of_gaussian_aux` (`RankR/RMT/AlignG.lean:134`) with `hn` and `hp`
removed by the tail shift, so it is the exact hypothesis of `prop_stacksvd_subspace`
(`RankR/Subspace.lean:498`) and the two compose with no glue. -/
theorem subspaceLaw_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    m.SubspaceLaw (∑ i, cc i) :=
  ⟨fun j => m.toStack.align_of_gaussian (m.gaussianNoise_toStack hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTable.stack_regime cc hreg).2.2 j⟩

/-! ### 3. `prop:stacksvd_subspace` at `r_i = 1`, unconditional -/

/-- **`prop:stacksvd_subspace`** (`main_paper.tex:834`) at `r_i = 1` for joint Gaussian noise:
the stacksvd subspace performance tends to one rank-one performance per spike of the core
matrix. The Layer 1 form `prop_stacksvd_subspace` (`RankR/Subspace.lean:498`) takes the
rank-`r` spiked law as a hypothesis; here that law is proved, so the only hypotheses left are
the model, the Gaussian noise and the regime. -/
theorem prop_stacksvd_subspace_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω)
      (limitStackR (fun i => (m.tbl i).θ) m.R cc) :=
  m.prop_stacksvd_subspace cc (m.subspaceLaw_of_gaussian hG hreg hc)

/-! ### 4. The Frobenius form at the canonical eigenframe -/

/-- The model field `hV` forces `r ≤ d N`: the `r` columns of `V N` are independent in
`ℝ^{d N}`. Stated with `Fintype.card (Fin (d N))` in place of `d N`, because that is the form
the frame constructor `topEigMat` and `prop_stacksvd_subspace_frobenius_eig` take. Mirror:
`SpikedModelR.rk_le_d` (`RankR/GeneralMain.lean:114`). -/
theorem r_le_d (m : UnalignedModel μ M n d r) (N : ℕ) : r ≤ Fintype.card (Fin (d N)) := by
  have h1 : ((m.V N)ᵀ * m.V N).rank = r := by
    rw [m.hV N]
    simp
  have h2 : ((m.V N)ᵀ * m.V N).rank ≤ (m.V N).rank := Matrix.rank_mul_le_right _ _
  have h3 : (m.V N).rank ≤ Fintype.card (Fin (d N)) := Matrix.rank_le_card_height (m.V N)
  omega

/-- **The top-`r` gap of the stack Gram matrix holds almost surely**, at every `N` with
`r ≤ min (∑_i n_i N) (d N)`. The stack is the affine family `X_stack = A Vᵀ + d^{-1/2} Z_stack`
(`stackX_eq`, `stackE_eq_smul`) and `Z_stack` has law `gaussianMatrix (∑_i n_i N) (d N)`
(`hasLaw_stackZu`), so `simpleSpec_ae_affine` (`RankR/RMT/SimplicityAffineR.lean:76`) applies
and `topGap_of_simpleSpec` reads the gap off the simple spectrum. Mirror:
`SpikedModel.singleTableLaw_topSimple_of_gaussian` (`RMT/Simplicity.lean:374`). -/
theorem ae_topGap_stackGram [NeZero M] (m : UnalignedModel μ M n d r)
    (hG : m.JointGaussianNoise) (N : ℕ) (hr : r ≤ min (∑ i, n i N) (d N)) :
    ∀ᵐ ω ∂(μ N), TopGap (m.stackGram N ω) (m.isHermitian_stackGram N ω) r := by
  classical
  have hdN : 0 < d N := (m.tbl ⟨0, NeZero.pos M⟩).hd N
  have hnN : 0 < ∑ i, n i N :=
    lt_of_lt_of_le ((m.tbl ⟨0, NeZero.pos M⟩).hn N)
      (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
        (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hdN
  have ht : (Real.sqrt (d N))⁻¹ ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  have hae : ∀ᵐ ω ∂(μ N), SimpleSpec
      ((m.signalPart N + (Real.sqrt (d N))⁻¹ • m.stackZu N ω)ᵀ *
        (m.signalPart N + (Real.sqrt (d N))⁻¹ • m.stackZu N ω))
      (isHermitian_transpose_mul_self _) r :=
    ((m.hasLaw_stackZu hG N).ae_iff
        (RankRStack.measurable_simpleSpec_affine (m.signalPart N)
          ((Real.sqrt (d N))⁻¹) r)).mpr
      (simpleSpec_ae_affine hnN hdN (m.signalPart N) ht r hr)
  filter_upwards [hae] with ω hω
  have hX : m.stackX N ω = m.signalPart N + (Real.sqrt (d N))⁻¹ • m.stackZu N ω := by
    rw [m.stackX_eq N ω, m.stackE_eq_smul N ω]
  have hgram : (m.signalPart N + (Real.sqrt (d N))⁻¹ • m.stackZu N ω)ᵀ *
      (m.signalPart N + (Real.sqrt (d N))⁻¹ • m.stackZu N ω) = m.stackGram N ω := by
    rw [← hX]
    rfl
  exact topGap_of_simpleSpec (m.isHermitian_stackGram N ω)
    (EdgeGlueR.simpleSpec_congr hgram _ (m.isHermitian_stackGram N ω) hω)

/-- **The gap hypothesis of `prop_stacksvd_subspace_frobenius_eig`, discharged.** Under the
regime the two size conditions `r ≤ ∑_i n_i N` and `r ≤ d N` hold at every large `N`, and
there the bad set is null by `ae_topGap_stackGram`, so its measure is eventually `0`. Only the
regime of the tables is used, not `0 < ∑_i c_i`. -/
theorem topGap_stackGram_whp_of_gaussian [NeZero M] (m : UnalignedModel μ M n d r)
    (hG : m.JointGaussianNoise) {cc : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.stackGram N ω) (m.isHermitian_stackGram N ω) r})
      atTop (𝓝 0) := by
  refine tendsto_const_nhds.congr' ?_
  have hstack : Tendsto (fun N => ∑ i, n i N) atTop atTop :=
    (m.toMultiTable.stack_regime cc hreg).1
  filter_upwards [hstack.eventually_ge_atTop r] with N hN
  have hrd : r ≤ d N := by simpa using m.r_le_d N
  exact (ae_iff.mp (m.ae_topGap_stackGram hG N (le_min hN hrd))).symm

/-- **`prop:stacksvd_subspace`** (`main_paper.tex:834`) at `r_i = 1` for joint Gaussian noise,
on the paper's own quantity `‖V̂_stacksvdᵀ V‖_F²`, with `V̂_stacksvd` the canonical top-`r`
eigenframe of `X_stackᵀ X_stack`. Both side conditions of
`prop_stacksvd_subspace_frobenius_eig` (`RankR/Frobenius.lean:886`) are proved above, so the
hypotheses are again the model, the Gaussian noise and the regime. -/
theorem prop_stacksvd_subspace_frobenius_eig_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModel μ M n d r)
    (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ
      (fun N ω => frobSq ((topEigMat (m.isHermitian_stackGram N ω) (m.r_le_d N))ᵀ * m.V N))
      (limitStackR (fun i => (m.tbl i).θ) m.R cc) :=
  m.prop_stacksvd_subspace_frobenius_eig cc (m.subspaceLaw_of_gaussian hG hreg hc)
    (fun N => m.r_le_d N) (m.topGap_stackGram_whp_of_gaussian hG hreg)

end UnalignedModel

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- **`prop:single_table` in rank `r`** for the general-`r_i` model with joint Gaussian noise,
in the mixed regime, with no side condition on `∑_i n_i`. This is
`UnalignedModelR.subspaceLawG_of_gaussian_aux` (`RankR/SubspaceGStack.lean:147`) with `hn` and
`hp` removed by the tail shift. -/
theorem subspaceLawG_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    m.SubspaceLawG (∑ i, cc i) :=
  ⟨fun j => m.toStackG.align_of_gaussian (m.gaussianNoise_toStackG hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTableShell.stack_regime cc hreg).2.2 j⟩

/-- **`prop:stacksvd_subspace`** (`main_paper.tex:834`) at general `r_i` for joint Gaussian
noise. The general-`r_i` twin of `UnalignedModel.prop_stacksvd_subspace_gaussian`; the Layer 1
form is `prop_stacksvd_subspace_general` (`RankR/SubspaceMain.lean:193`).

Scope (Track A audit F1, 2026-09-02): each table is a `SpikedModelR`, whose fields `hθnn` and
`hθanti : StrictAnti θ` require distinct nonnegative `θ_ij` inside a table. The paper makes no
ordering assumption (`main_paper.tex:767`) and assumes distinct entries for svdstack only
(`:768`), so this statement is narrower than the paper's at a table with a repeated `θ_ij`.
The stacksvd proofs read neither field; the order is a model field, not a hypothesis. -/
theorem prop_stacksvd_subspace_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackRG N ω)
      (limitStackRG (fun i => (m.tbl i).θ) m.R cc) :=
  m.prop_stacksvd_subspace_general cc (m.subspaceLawG_of_gaussian hG hreg hc)

/-! ### 4b. The Frobenius form at general `r_i` -/

/-- The general-`r_i` twin of `UnalignedModel.r_le_d` (`RankR/SubspaceGaussian.lean:134`): the
model field `hV` forces `r ≤ d N`, read through `Fintype.card (Fin (d N))`. Same proof. -/
theorem r_le_d (m : UnalignedModelR μ M n d r rk) (N : ℕ) : r ≤ Fintype.card (Fin (d N)) := by
  have h1 : ((m.V N)ᵀ * m.V N).rank = r := by
    rw [m.hV N]
    simp
  have h2 : ((m.V N)ᵀ * m.V N).rank ≤ (m.V N).rank := Matrix.rank_mul_le_right _ _
  have h3 : (m.V N).rank ≤ Fintype.card (Fin (d N)) := Matrix.rank_le_card_height (m.V N)
  omega

/-- The general-`r_i` twin of `UnalignedModel.prop_stacksvd_subspace_frobenius`
(`RankR/Frobenius.lean:871`): `prop:stacksvd_subspace` (`main_paper.tex:834`) on the paper's
own quantity `‖V̂_stacksvdᵀ V‖_F²`, for any selection `Y` that is a top-`r` frame of the stack
Gram matrix with probability tending to one. Same proof, through
`prop_stacksvd_subspace_general` (`RankR/SubspaceMain.lean:193`) in place of
`prop_stacksvd_subspace`. -/
theorem prop_stacksvd_subspace_general_frobenius (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (law : m.SubspaceLawG (∑ i, c i))
    (Y : ∀ N, Ω N → Matrix (Fin (d N)) (Fin r) ℝ)
    (hY : Tendsto (fun N => μ N {ω | ¬ IsTopFrame (m.stackGramG N ω)
      (m.isHermitian_stackGramG N ω) (Y N ω)}) atTop (𝓝 0)) :
    TendstoInProb μ (fun N ω => frobSq ((Y N ω)ᵀ * m.V N))
      (limitStackRG (fun i => (m.tbl i).θ) m.R c) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.prop_stacksvd_subspace_general c law)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hY (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hframe
  refine hω ?_
  rw [frobenius_eq_perf hframe (m.V N)]
  rfl

/-- The general-`r_i` twin of `UnalignedModel.prop_stacksvd_subspace_frobenius_eig`
(`RankR/Frobenius.lean:889`): the same statement at the canonical top-`r` eigenframe
`topEigMat` of the stack Gram matrix, with the top-`r` gap as an explicit hypothesis. Same
proof. -/
theorem prop_stacksvd_subspace_general_frobenius_eig (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (law : m.SubspaceLawG (∑ i, c i)) (hrd : ∀ N, r ≤ Fintype.card (Fin (d N)))
    (hgap : Tendsto (fun N => μ N {ω | ¬ TopGap (m.stackGramG N ω)
      (m.isHermitian_stackGramG N ω) r}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((topEigMat (m.isHermitian_stackGramG N ω) (hrd N))ᵀ * m.V N))
      (limitStackRG (fun i => (m.tbl i).θ) m.R c) := by
  refine m.prop_stacksvd_subspace_general_frobenius c law _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hgap (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω (isTopFrame_topEigMat (m.isHermitian_stackGramG N ω) (hrd N) hgapN)

/-- The general-`r_i` twin of `UnalignedModel.ae_topGap_stackGram`
(`RankR/SubspaceGaussian.lean:148`): the top-`r` gap of the stack Gram matrix holds almost
surely, at every `N` with `r ≤ min (∑_i n_i N) (d N)`. Same proof, through the `G`-suffixed
stack objects (`stackX_eqG`, `stackE_eqG`, `signalPartG`, `stackZG`) and `hasLaw_stackZG`
(`RankR/SubspaceGStack.lean:82`) in place of `hasLaw_stackZu`. -/
theorem ae_topGap_stackGramG [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) (N : ℕ) (hr : r ≤ min (∑ i, n i N) (d N)) :
    ∀ᵐ ω ∂(μ N), TopGap (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r := by
  classical
  have hdN : 0 < d N := (m.tbl ⟨0, NeZero.pos M⟩).hd N
  have hnN : 0 < ∑ i, n i N :=
    lt_of_lt_of_le ((m.tbl ⟨0, NeZero.pos M⟩).hn N)
      (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
        (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hdN
  have ht : (Real.sqrt (d N))⁻¹ ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  have hae : ∀ᵐ ω ∂(μ N), SimpleSpec
      ((m.signalPartG N + (Real.sqrt (d N))⁻¹ • m.stackZG N ω)ᵀ *
        (m.signalPartG N + (Real.sqrt (d N))⁻¹ • m.stackZG N ω))
      (isHermitian_transpose_mul_self _) r :=
    ((m.hasLaw_stackZG hG N).ae_iff
        (RankRStack.measurable_simpleSpec_affine (m.signalPartG N)
          ((Real.sqrt (d N))⁻¹) r)).mpr
      (simpleSpec_ae_affine hnN hdN (m.signalPartG N) ht r hr)
  filter_upwards [hae] with ω hω
  have hX : m.stackXG N ω = m.signalPartG N + (Real.sqrt (d N))⁻¹ • m.stackZG N ω := by
    rw [m.stackX_eqG N ω, m.stackE_eqG N ω]
  have hgram : (m.signalPartG N + (Real.sqrt (d N))⁻¹ • m.stackZG N ω)ᵀ *
      (m.signalPartG N + (Real.sqrt (d N))⁻¹ • m.stackZG N ω) = m.stackGramG N ω := by
    rw [← hX]
    rfl
  exact topGap_of_simpleSpec (m.isHermitian_stackGramG N ω)
    (EdgeGlueR.simpleSpec_congr hgram _ (m.isHermitian_stackGramG N ω) hω)

/-- The general-`r_i` twin of `UnalignedModel.topGap_stackGram_whp_of_gaussian`
(`RankR/SubspaceGaussian.lean:181`): the gap hypothesis of
`prop_stacksvd_subspace_general_frobenius_eig`, discharged under the regime. Same proof,
through `toMultiTableShell.stack_regime` (`RankR/SubspaceGStack.lean:69`) in place of
`toMultiTable.stack_regime`. -/
theorem topGap_stackGramG_whp_of_gaussian [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r})
      atTop (𝓝 0) := by
  refine tendsto_const_nhds.congr' ?_
  have hstack : Tendsto (fun N => ∑ i, n i N) atTop atTop :=
    (m.toMultiTableShell.stack_regime cc hreg).1
  filter_upwards [hstack.eventually_ge_atTop r] with N hN
  have hrd : r ≤ d N := by simpa using m.r_le_d N
  exact (ae_iff.mp (m.ae_topGap_stackGramG hG N (le_min hN hrd))).symm

/-- The general-`r_i` twin of `UnalignedModel.prop_stacksvd_subspace_frobenius_eig_gaussian`
(`RankR/SubspaceGaussian.lean:198`): `prop:stacksvd_subspace` (`main_paper.tex:834`) on the
paper's own quantity, with the law and the gap both discharged, so the hypotheses are the
model, the Gaussian noise and the regime. Same proof, through `subspaceLawG_of_gaussian` in
place of `subspaceLaw_of_gaussian`. -/
theorem prop_stacksvd_subspace_general_frobenius_eig_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ
      (fun N ω => frobSq ((topEigMat (m.isHermitian_stackGramG N ω) (m.r_le_d N))ᵀ * m.V N))
      (limitStackRG (fun i => (m.tbl i).θ) m.R cc) :=
  m.prop_stacksvd_subspace_general_frobenius_eig cc (m.subspaceLawG_of_gaussian hG hreg hc)
    (fun N => m.r_le_d N) (m.topGap_stackGramG_whp_of_gaussian hG hreg)

end UnalignedModelR

end StackedSVD

namespace StackedSVD.RankR.Example

open MeasureTheory Filter Topology

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ}

/-! ### 5. The worked example of Section 7, stacksvd half, unconditional -/

/-- **`eq:perf_rankr_ex_stacksvd`** (`main_paper.tex:848`) for joint Gaussian noise. This is
`example_perfStackR_tendsto` (`RankR/SubspaceMain.lean:156`) with its hypothesis `law` discharged
by `UnalignedModel.subspaceLaw_of_gaussian`, so the only hypotheses left are the two tables of
`eq:psi_equation`, the common regime `c` and the Gaussian noise. The svdstack twin is
`example_perfR_tendsto_gaussian` (`RankR/Example.lean:976`) and the two together are the
comparison of `main_paper.tex:864`. -/
theorem example_perfStackR_tendsto_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω) (stacksvdEx θ c (Real.sin ψ)) := by
  have hc2 : (∑ i : Fin 2, (fun _ : Fin 2 => c) i) = 2 * c := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    simp [two_mul]
  have hcsum : 0 < ∑ i : Fin 2, (fun _ : Fin 2 => c) i := by
    rw [hc2]
    linarith
  have law : m.SubspaceLaw (2 * c) := by
    have h := m.subspaceLaw_of_gaussian hG (cc := fun _ : Fin 2 => c) hreg hcsum
    rwa [hc2] at h
  exact example_perfStackR_tendsto m θ c ψ hθ hR law

end StackedSVD.RankR.Example
