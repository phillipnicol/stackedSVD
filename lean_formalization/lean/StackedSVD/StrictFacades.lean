/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Sup
import StackedSVD.StackSVD.Weighted

/-!
# Strict facades: `thm:stacksvd_binary_optimal_svd_stack` and `prop:dominance`

`notes/SERVER_TODO.md`, "After the third external review", item P2. The paper states two
comparisons with strict inequalities under extra hypotheses:

* `thm:stacksvd_binary_optimal_svd_stack` (`main_paper.tex:626`): binary-weighted stacksvd on
  the detectable set beats optimally weighted svdstack strictly once `β₂ > 0`, that is, once
  at least two tables clear their own detection threshold.
* `prop:dominance` (`main_paper.tex:634`): optimally weighted stacksvd beats optimally
  weighted svdstack strictly once at least two tables carry signal, and beats unweighted
  stacksvd strictly once `θ_i²/c_i` is not constant across tables.

The scalar strict inequalities (`Scalars.svdstackOpt_lt_binary`,
`Scalars.svdstackOpt_lt_stackSVDLimitW`, `Scalars.stackSVDLimit_lt_stackSVDLimitW`) already
exist in `Scalars.lean`. This file only assembles them with the model-side convergence
already proved for the weak (`≤`) facades
(`MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian`,
`StackSVD/Weighted.lean:214`; `MultiTableModel.prop_dominance_gaussian`,
`RMT/Het/Sup.lean:410`): the random matrix theory content does not change, only the scalar
comparison. No new mathematics is proved here.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### `thm:stacksvd_binary_optimal_svd_stack`, strict half -/

/-- **`thm:stacksvd_binary_optimal_svd_stack`** (`main_paper.tex:626`), strict half, Gaussian
noise: with every kept table at `c_i ≤ 1` and at least two tables above their own detection
threshold `c_i < θ_i⁴` (`hcard`, the paper's hypothesis `β₂ > 0`), the binary weighting on the
detectable set `S` converges to `Scalars.binaryStackSVDLimit S θ c`, and that limit is
*strictly* above optimally weighted svdstack. `hcard` replaces the weak facade's `hdet`
(existence of one detectable table): `2 ≤ S.card` already forces `S` nonempty. Convergence
comes from `stackPerfW_binary_tendsto_gaussian`; the strict inequality from
`Scalars.svdstackOpt_lt_binary`. -/
theorem thm_stacksvd_binary_optimal_svd_stack_gaussian_strict
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1)
    (hcard : 2 ≤ (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4).card)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW
        (fun i => if i ∈ Finset.univ.filter (fun i => c i < (m.tbl i).θ ^ 4) then 1 else 0)
        N ω)
      (Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
        (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        < Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
            (fun i => (m.tbl i).θ) c := by
  classical
  set S : Finset (Fin M) := Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4 with hSdef
  have hSne : S.Nonempty := Finset.card_pos.mp (by omega)
  have : NeZero S.card := ⟨(Finset.card_pos.mpr hSne).ne'⟩
  exact ⟨m.stackPerfW_binary_tendsto_gaussian S c hc hreg hG,
    Scalars.svdstackOpt_lt_binary hc hc1 hcard⟩

/-! ### `prop:dominance`, strict halves -/

/-- **`prop:dominance`** (`main_paper.tex:634`), both strict halves, Gaussian noise: above the
recovery threshold `hthr` (`∑ θ_i⁴/c_i > 1`, the condition of `thm:stacksvd_weighted` at
`main_paper.tex:471`; `eq:assumption4` at line 1368 is the general-`w` form, equivalent to
it only at the optimal weights, line 1602), with at least two tables carrying signal
(`htwo`, the paper's hypothesis for the svdstack comparison) and `θ_i²/c_i` not constant
across tables (`hnc`, the paper's hypothesis for the unweighted stacksvd
comparison), optimally weighted stacksvd converges to `Scalars.stackSVDLimitW` and strictly
dominates both optimally weighted svdstack and unweighted stacksvd. The weak facade's `hθ`
(some `θ_i ≠ 0`) is derived from `hthr` via `Scalars.exists_ne_zero_of_thr`, since
`∑ θ_i⁴/c_i > 0` already forces some `θ_i ≠ 0`. -/
theorem prop_dominance_gaussian_strict
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hthr : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (htwo : ∃ i j, i ≠ j ∧ (m.tbl i).θ ≠ 0 ∧ (m.tbl j).θ ≠ 0)
    (hnc : ∃ i j, (m.tbl i).θ ^ 2 * c j ≠ (m.tbl j).θ ^ 2 * c i) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c := by
  have hθ : ∃ i, (m.tbl i).θ ≠ 0 := Scalars.exists_ne_zero_of_thr hthr
  refine ⟨(m.prop_dominance_gaussian c hc hθ hreg hG).1,
    Scalars.svdstackOpt_lt_stackSVDLimitW hc hthr htwo,
    Scalars.stackSVDLimit_lt_stackSVDLimitW hc hthr hnc⟩

/-- **`prop:dominance`**, svdstack half only, Gaussian noise: same as
`prop_dominance_gaussian_strict` but keeps the unweighted-stacksvd comparison weak, so the
strict svdstack claim carries only its own hypothesis `htwo` and not `hnc`. -/
theorem prop_dominance_gaussian_strict_svdstack
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hthr : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (htwo : ∃ i j, i ≠ j ∧ (m.tbl i).θ ≠ 0 ∧ (m.tbl j).θ ≠ 0) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c := by
  have hθ : ∃ i, (m.tbl i).θ ≠ 0 := Scalars.exists_ne_zero_of_thr hthr
  have hbase := m.prop_dominance_gaussian c hc hθ hreg hG
  exact ⟨hbase.1, Scalars.svdstackOpt_lt_stackSVDLimitW hc hthr htwo, hbase.2.2⟩

/-- **`prop:dominance`**, unweighted-stacksvd half only, Gaussian noise: same as
`prop_dominance_gaussian_strict` but keeps the svdstack comparison weak, so the strict
unweighted-stacksvd claim carries only its own hypothesis `hnc` and not `htwo`. -/
theorem prop_dominance_gaussian_strict_unweighted
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hthr : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (hnc : ∃ i j, (m.tbl i).θ ^ 2 * c j ≠ (m.tbl j).θ ^ 2 * c i) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c := by
  have hθ : ∃ i, (m.tbl i).θ ≠ 0 := Scalars.exists_ne_zero_of_thr hthr
  have hbase := m.prop_dominance_gaussian c hc hθ hreg hG
  exact ⟨hbase.1, hbase.2.1, Scalars.stackSVDLimit_lt_stackSVDLimitW hc hthr hnc⟩

end MultiTableModel

end StackedSVD
