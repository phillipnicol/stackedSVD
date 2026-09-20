/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Sup
import StackedSVD.RMT.General.Edge.Sup
import StackedSVD.SVDStack.Gram
import StackedSVD.SVDStack.Main
import StackedSVD.StackSVD
import StackedSVD.StackSVD.Main
import StackedSVD.ThetaEst

/-!
# Stage 1 of the non-Gaussian extension: the Layer 1 statements

The Layer 1 half of the 12 target statements of `notes/archive/prop_single_table_general.md` (status
`user OK`, 2026-09-09): the two bridges and the four representative corollaries, proved by
unit G11, 2026-09-09. The single-table half was `RMT/General/Statements.lean`, retired by
follow-up item F35 (2026-09-09); its structure `ResolventFormsC` now lives in
`RMT/General/Defs.lean`.

This file sits outside `RMT/` on purpose. Gate 8 (`scripts/check_core_imports.py`) puts every
module under `StackedSVD.RMT.` into the Gaussian RMT core, and a core module may not import
`ThetaEst`, `StackSVD` or `SVDStack`, which these six statements need.

Since Stage 3 (`notes/stage3_edge.md`, the sharp upper edge at a general law,
`RMT/General/Edge/Sup.lean`) the file also holds the four corollaries with the edge
hypothesis discharged, named `*_of_moments`: the paper's noise class alone
(`assum:general_noise`: i.i.d. entries of one law with mean `0`, variance `1` and a finite
fourth moment), plus the Lebesgue density `hac` that the simplicity of the top eigenvalue
uses.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### the two bridges the Layer 1 corollaries need -/

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- Independent tables with i.i.d. entries of law `ν` are independent tables
(`IndepNoise`, `SVDStack/Gram.lean`). Two lines, the twin of
`JointGaussianNoise.indepNoise`. -/
theorem JointGeneralNoise.indepNoise {m : MultiTableModel μ M n d} {ν : Measure ℝ}
    [IsProbabilityMeasure ν] (hG : m.JointGeneralNoise ν) : m.IndepNoise :=
  ⟨fun N i => noiseMatrix ν (n i N) (d N), fun _ _ => inferInstance, hG⟩

set_option linter.unusedSectionVars false in
/-- The stacking map on plain product types, at a general law `ν`: the twin of
`measurePreserving_stackPi` (`StackSVD.lean:240`), a `Measure.pi_eq` computation that never
uses the Gaussian law. -/
private theorem measurePreserving_stackPi_general (ν : Measure ℝ) [SigmaFinite ν] (N : ℕ) :
    MeasurePreserving
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
        Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2)
      (Measure.pi fun i : Fin M =>
        Measure.pi fun _ : Fin (n i N) => Measure.pi fun _ : Fin (d N) => ν)
      (Measure.pi fun _ : Fin (∑ i, n i N) =>
        Measure.pi fun _ : Fin (d N) => ν) := by
  have hmeas : Measurable
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
        Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2) :=
    measurable_pi_lambda _ fun r =>
      (measurable_pi_apply (finSigmaFinEquiv.symm r).2).comp
        (measurable_pi_apply (finSigmaFinEquiv.symm r).1)
  refine ⟨hmeas, ?_⟩
  refine (Measure.pi_eq fun s hs => ?_).symm
  rw [Measure.map_apply hmeas (MeasurableSet.univ_pi hs)]
  have hpre :
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
          Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2) ⁻¹' Set.univ.pi s
        = Set.univ.pi fun i : Fin M =>
            Set.univ.pi fun j : Fin (n i N) => s (finSigmaFinEquiv ⟨i, j⟩) := by
    ext Zs
    simp only [Set.mem_preimage, Set.mem_univ_pi]
    constructor
    · intro hZ i j
      have h1 := hZ (finSigmaFinEquiv ⟨i, j⟩)
      rw [Equiv.symm_apply_apply] at h1
      exact h1
    · intro hZ r
      obtain ⟨p, rfl⟩ : ∃ p, finSigmaFinEquiv p = r :=
        ⟨finSigmaFinEquiv.symm r, Equiv.apply_symm_apply _ _⟩
      obtain ⟨i, j⟩ := p
      rw [Equiv.symm_apply_apply]
      exact hZ i j
  have hRHS : ∏ r : Fin (∑ i, n i N),
        (Measure.pi fun _ : Fin (d N) => ν) (s r)
      = ∏ i : Fin M, ∏ j : Fin (n i N),
          (Measure.pi fun _ : Fin (d N) => ν) (s (finSigmaFinEquiv ⟨i, j⟩)) := by
    rw [← Equiv.prod_comp finSigmaFinEquiv
      (fun r => (Measure.pi fun _ : Fin (d N) => ν) (s r))]
    exact Fintype.prod_sigma _
  rw [hpre, hRHS]
  simp only [Measure.pi_pi]

set_option linter.unusedSectionVars false in
/-- The same statement in the `Matrix` types of `Defs.lean`, at a general law `ν`: the twin of
`measurePreserving_stackFun` (`StackSVD.lean:287`). -/
private theorem measurePreserving_stackFun_general (ν : Measure ℝ) [SigmaFinite ν] (N : ℕ) :
    MeasurePreserving
      (fun (Zs : (i : Fin M) → Matrix (Fin (n i N)) (Fin (d N)) ℝ) =>
        Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
          (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => Zs p.1 p.2 k))
      (Measure.pi fun i : Fin M => noiseMatrix ν (n i N) (d N))
      (noiseMatrix ν (∑ i, n i N) (d N)) :=
  measurePreserving_stackPi_general ν N

/-- The stacked noise has i.i.d. entries of law `ν`: the twin of `stack_law`
(`StackSVD.lean:298`), same proof with `ν` in place of `gaussianReal 0 1`. -/
theorem stack_law_general [NeZero M] (m : MultiTableModel μ M n d) {ν : Measure ℝ}
    [SigmaFinite ν] (h : m.JointGeneralNoise ν) : m.stack.GeneralNoise ν := by
  intro N
  exact (measurePreserving_stackFun_general ν N).fun_comp_hasLaw (h N)

/-! ### four representative Layer 1 corollaries -/

/-- `prop:stacksvd_general` (`Main.prop_stacksvd_general`) at four moments: one edge
hypothesis, at the stack. -/
theorem prop_stacksvd_general_of_general [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.stack.W0 N ω) (m.stack.isHermitian_W0 N ω)
        ≤ bulkEdge (∑ i, c i) + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  have hne : Nonempty (Fin M) := ⟨0⟩
  have hcpos : 0 < ∑ i, c i := Finset.sum_pos (fun i _ => hc i) Finset.univ_nonempty
  have := hν.prob
  exact m.prop_stacksvd_general c
    (SpikedModel.singleTableLaw_of_general hcpos m.stack (m.stack_regime c hreg) hν hac
      (m.stack_law_general hG) hedge)

/-- `lem:delocalization` (`Main.lem_delocalization`) at four moments: one edge hypothesis per
table. -/
theorem lem_delocalization_of_general [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ (i : Fin M), ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge (c i) + ε}) atTop (𝓝 1))
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j)) := by
  have := hν.prob
  exact m.lem_delocalization c
    (fun k => SpikedModel.singleTableLaw_of_general (hc k) (m.tbl k) (hreg k) hν hac
      (m.generalNoise_of_joint hG k) (hedge k))
    hG.indepNoise hij

/-- `thm:svd_stack_general` (`Main.thm_svd_stack_general`) at four moments. -/
theorem thm_svd_stack_general_of_general [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ (i : Fin M), ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge (c i) + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β) := by
  have := hν.prob
  exact m.thm_svd_stack_general c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_general (hc i) (m.tbl i) (hreg i) hν hac
      (m.generalNoise_of_joint hG i) (hedge i))
    hG.indepNoise

/-- `thm:theta_est` (`Main.thm_theta_est`) at four moments: `ThetaEstLaw` is Stage 0, the
single-table law of table `i` is Stage 1, so only the edge of table `i` stays open. -/
theorem thm_theta_est_of_general [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge ci + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ := by
  have := hν.prob
  exact m.thm_theta_est hci.le hthr
    (SpikedModel.singleTableLaw_of_general hci (m.tbl i) hregi hν hac
      (m.generalNoise_of_joint hG i) hedge)
    (m.thetaEstLaw_of_general hν hij hG hregj)

/-! ### the same four corollaries with the edge discharged (Stage 3) -/

/-- `prop:stacksvd_general` at `assum:general_noise`: `prop_stacksvd_general_of_general` with
its edge hypothesis supplied by `SpikedModel.lamMax_W0_edge_of_general` at the stack. -/
theorem prop_stacksvd_general_of_moments [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  have hne : Nonempty (Fin M) := ⟨0⟩
  have hcpos : 0 < ∑ i, c i := Finset.sum_pos (fun i _ => hc i) Finset.univ_nonempty
  have := hν.prob
  exact m.prop_stacksvd_general_of_general c hc hreg hν hac hG
    (m.stack.lamMax_W0_edge_of_general hcpos (m.stack_regime c hreg) hν
      (m.stack_law_general hG))

/-- `lem:delocalization` at `assum:general_noise`: `lem_delocalization_of_general` with the
edge of every table supplied by `SpikedModel.lamMax_W0_edge_of_general`. -/
theorem lem_delocalization_of_moments [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j)) := by
  have := hν.prob
  exact m.lem_delocalization_of_general c hc hreg hν hac hG
    (fun k => (m.tbl k).lamMax_W0_edge_of_general (hc k) (hreg k) hν
      (m.generalNoise_of_joint hG k)) hij

/-- `thm:svd_stack_general` at `assum:general_noise`: `thm_svd_stack_general_of_general`
with the edge of every table supplied by `SpikedModel.lamMax_W0_edge_of_general`. -/
theorem thm_svd_stack_general_of_moments [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β) := by
  have := hν.prob
  exact m.thm_svd_stack_general_of_general c β hc hβdef hthr hreg hν hac hG
    (fun i => (m.tbl i).lamMax_W0_edge_of_general (hc i) (hreg i) hν
      (m.generalNoise_of_joint hG i))

/-- `thm:theta_est` at `assum:general_noise`: `thm_theta_est_of_general` with the edge of
table `i` supplied by `SpikedModel.lamMax_W0_edge_of_general`. -/
theorem thm_theta_est_of_moments [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ := by
  have := hν.prob
  exact m.thm_theta_est_of_general hij hci hthr hν hac hG hregi hregj
    ((m.tbl i).lamMax_W0_edge_of_general hci hregi hν (m.generalNoise_of_joint hG i))

end MultiTableModel

end StackedSVD
