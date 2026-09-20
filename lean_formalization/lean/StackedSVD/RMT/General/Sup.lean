/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.R3minus
import StackedSVD.RMT.General.Simplicity
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.FormsGeneral
import StackedSVD.RMT.General.DelocAlign
import StackedSVD.RMT.General.DelocUniformSup

/-! # Item Sup at a general law: the Stage 1 endpoint

`SpikedModel.resolventLimits_of_general` builds the seven-field `ResolventLimits` interface
(item R5, `RMT/R5.lean`) for a spiked model with i.i.d. noise of a fixed law `ν` (`NoiseLaw ν`)
in the proportional regime, from the edge hypothesis `hedge` and the six complex forms of
`resolventFormsC_of_general'` (`RMT/General/FormsGeneral.lean`, unit F3). No
`hn2 : ∀ N, 2 ≤ n N` is needed: item R0's block split does not appear at a general law.

`SpikedModel.singleTableLaw_of_general` is the Stage 1 endpoint, `prop:single_table` at
`assum:general_noise` (`docs/THEOREMS.md` 3.1), both regimes. Route: one `by_cases` on
`c < θ⁴`, as `RMT/Full.lean`. The supercritical branch reads `align` and `lamMax` off
`ResolventLimits` through the model-free `RMT/R5.lean`; the subcritical branch reads them off
unit G9 (`RMT/General/R3minus.lean`). `delocUniform` comes from unit G7
(`delocUniform_of_general`) in both regimes and `topSimple` from unit G8
(`singleTableLaw_topSimple_of_general`, `RMT/General/Simplicity.lean`), every `N`, no
threshold.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! ### (b) the seven-field interface, with the edge as a hypothesis -/

/-- **(b).** `ResolventLimits` from the six complex forms and the edge, taken as the
hypothesis `hedge` in the shape of the `edge` field of the structure (the sharp constant
`bulkEdge c`). Item T does the transfer from complex `z` to real `z > bulkEdge c`. -/
theorem resolventLimits_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)) :
    m.ResolventLimits c := by
  classical
  have hd : Tendsto d atTop atTop := hreg.2.1
  have hvdot : ∀ N, WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 := m.dotProduct_v_self
  have hedge0 : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0) :=
    fun ε hε => tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet) (hedge ε hε)
  have H : m.ResolventFormsC c := m.resolventFormsC_of_general' hc hreg hν hG
  have hgdot : TendstoInProb μ (fun N ω => m.gvec N ω ⬝ᵥ m.gvec N ω) 1 :=
    m.tendstoInProb_dotProduct_gvec_general hν hG hd
  have hgbad :
      Tendsto (fun N => μ N {ω : Ω N | m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω : Ω N | (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|})
      (fun N ω hω => ?_) (hgdot 3 (by norm_num))
    have hnot : ¬ (m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    have h4 : (4 : ℝ) < m.gvec N ω ⬝ᵥ m.gvec N ω := not_le.mp hnot
    change (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|
    rw [abs_of_nonneg (by linarith)]
    linarith
  have hnormvv : Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
        WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4}ᶜ) atTop (𝓝 0) := by
    have hset : ∀ N : ℕ, {ω : Ω N |
        WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
          WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4}ᶜ = (∅ : Set (Ω N)) := by
      intro N
      ext ω
      simp [hvdot N]
    simp only [hset, measure_empty]
    exact tendsto_const_nhds
  have hnormgg : Tendsto (fun N => μ N {ω : Ω N |
      m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4 ∧ m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) hgbad
    have hnot : ¬ (m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4 ∧ m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    exact fun hg => hnot ⟨hg, hg⟩
  have hnormvg : Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
        m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) hgbad
    have hnot : ¬ (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
      m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    refine fun hg => hnot ⟨?_, hg⟩
    rw [hvdot N]
    norm_num
  have hzeroL : ∀ x : ℝ, Tendsto
      (fun η : ℝ => (fun _ : ℂ => (0 : ℂ)) ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 (((0 : ℝ) : ℂ))) := by
    intro x
    simp only [Complex.ofReal_zero]
    exact tendsto_const_nhds
  refine ⟨hedge, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) (fun N _ => WithLp.ofLp (m.v N)) hedge0 hnormvv
      (MP.mC c) (MP.m c x) (fun z hz => H.vvC z hz) hx (MP.tendsto_mC hc hx)
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0 m.gvec m.gvec
      hedge0 hnormgg (MP.mC c) (MP.m c x) (fun z hz => H.ggC z hz) hx
      (MP.tendsto_mC hc hx)
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) m.gvec hedge0 hnormvg (fun _ => (0 : ℂ)) 0
      (fun z hz => by simpa only [sub_zero] using H.vgC z hz) hx (hzeroL x)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) (fun N _ => WithLp.ofLp (m.v N)) hedge0 hnormvv
      (MP.mCDeriv c) (MP.mDeriv c x) (fun z hz => H.vv2C z hz) hx
      (MP.tendsto_mCDeriv hc hx)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0 m.gvec m.gvec
      hedge0 hnormgg (MP.mCDeriv c) (MP.mDeriv c x) (fun z hz => H.gg2C z hz) hx
      (MP.tendsto_mCDeriv hc hx)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) m.gvec hedge0 hnormvg (fun _ => (0 : ℂ)) 0
      (fun z hz => by simpa only [sub_zero] using H.vg2C z hz) hx (hzeroL x)

/-! ### (c) the Stage 1 target -/

/-- **Item Sym at a general law, the target statement of the plan note (unit G7).** The
`delocUniform` field, both regimes, uniform over the unit sphere of `v ⊥`, from the edge
hypothesis `hedge`: the six complex forms (`resolventFormsC_of_general'`, unit F3) and the
seven-field interface (`resolventLimits_of_general`, (b) above) feed
`delocUniform_of_general'` (`RMT/General/DelocUniformSup.lean`), whose core is the window
bound `delocUniform_of_lamMax` (`RMT/General/DelocUniform.lean`). No density hypothesis:
the window bound caps the whole top eigenspace at once, so it does not need a simple top
eigenvalue. -/
theorem delocUniform_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)) :
    ∀ ε > 0, Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w})
      atTop (𝓝 0) :=
  m.delocUniform_of_general' hc (m.resolventFormsC_of_general' hc hreg hν hG)
    (m.resolventLimits_of_general hc hreg hν hG hedge) hreg hν hG

/-- **(c). `prop:single_table` at `assum:general_noise`, given the edge.** Both regimes, no
`hn2 : ∀ N, 2 ≤ n N`. Against `singleTableLaw_of_gaussian` (`RMT/Full.lean`) this replaces
`hG : m.GaussianNoise` by the three hypotheses `hν`, `hac`, `hG` and adds `hedge`. -/
theorem singleTableLaw_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)) :
    m.SingleTableLaw c := by
  by_cases hθ : c < m.θ ^ 4
  · have RL := m.resolventLimits_of_general hc hreg hν hG hedge
    have := hν.prob
    exact ⟨align_tendstoInProb RL hc hθ, m.delocUniform_of_general hc hreg hν hG hedge,
      lamMax_tendstoInProb RL hc hθ, singleTableLaw_topSimple_of_general m hac hG⟩
  · have hθ' : m.θ ^ 4 ≤ c := not_lt.mp hθ
    have H := m.resolventFormsC_of_general' hc hreg hν hG
    have RL := m.resolventLimits_of_general hc hreg hν hG hedge
    have := hν.prob
    exact ⟨m.align_tendstoInProb_of_subcritical_general' hc H RL hreg hν hac hG hθ',
      m.delocUniform_of_general hc hreg hν hG hedge,
      m.lamMax_tendstoInProb_of_subcritical_general RL hc hreg hν hG hθ',
      singleTableLaw_topSimple_of_general m hac hG⟩

end SpikedModel

end StackedSVD
