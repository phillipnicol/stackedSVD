/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.StackSVD
import StackedSVD.Scalars

/-!
# `prop:stacksvd_general`: the unweighted stacksvd limit

Moved out of `StackSVD.lean` (F28, 2026-09-08): the Layer 1 and Gaussian forms of
`prop:stacksvd_general` and `thm:simple_thm1`, stacksvd half. F32 (2026-09-08) drops the local
closed form `stackSVDLimit_const_eq`; the theorem now cites the public
`Scalars.simple_thm1_stacksvd` directly, so the file imports `StackedSVD.Scalars`. No proof
changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section Stack

-- One table at least. `Fin M` must be nonempty: the stack needs a row and a shared `v`.
variable [NeZero M]

/-! ### `prop:stacksvd_general` -/

omit [NeZero M] in
/-- `prop:stacksvd_general`, Layer 1 form: the single-table law for the stack at
`(‖θ‖₂, ‖c‖₁)` gives the stacksvd limit. The proof is `law.align` rewritten by
`betaSq_sqrt`.

**On the measure.** This statement takes no `IsProbabilityMeasure` instance, so it holds for
every measure family `μ`, the zero measure included; there `TendstoInProb` is trivial and the
statement carries no content (mechanical audit 2026-08-31, findings 1 to 3; probe P1). That
is a feature of a Layer 1 implication, not a gap: the hypothesis `law` is what a caller must
produce, and the Gaussian corollary `prop_stacksvd_general_gaussian` pins the measure through
`GaussianNoise`, which forces `IsProbabilityMeasure (μ N)`. The same reading applies to
`thm_stacksvd_weighted` in `StackSVDWeighted.lean`. -/
theorem prop_stacksvd_general [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : m.stack.SingleTableLaw (∑ i, c i)) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  have hT : (0 : ℝ) ≤ ∑ i, (m.tbl i).θ ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  have halign := law.align
  rw [stack_theta, betaSq_sqrt hT] at halign
  exact halign

omit [NeZero M] in
/-- `prop:stacksvd_general`, Layer 2 form: the same conclusion from the proportional regime of
each table and the joint Gaussian law, with no `SingleTableLaw` hypothesis. The proof is
`prop_stacksvd_general` applied to `SpikedModel.singleTableLaw_of_gaussian` of the stack,
whose regime is `stack_regime` and whose noise law is `stack_law`; the aspect ratio
`∑ i, c i` is positive because `hc` is positive and `Fin M` is nonempty. This is the missing
link of `notes/archive/audit_scope_2026-08-29.md` section 3.2 item 1. -/
theorem prop_stacksvd_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  have hne : Nonempty (Fin M) := ⟨0⟩
  have hcpos : 0 < ∑ i, c i := Finset.sum_pos (fun i _ => hc i) Finset.univ_nonempty
  exact m.prop_stacksvd_general c
    (SpikedModel.singleTableLaw_of_gaussian hcpos m.stack (m.stack_regime c hreg)
      (m.stack_law hG))

set_option linter.unusedVariables false in
/-- `prop:stacksvd_general` in the paper's own form `|⟨v̂_stacksvd, v⟩|²`, for any selection
`vhat` of a unit top eigenvector of the stack Gram matrix `X_stackᵀ X_stack`. The projector
form `overlap` equals that square on the event where the top eigenvalue is simple
(`overlap_eq_inner_sq`), which `SingleTableLaw.topSimple` of the stack gives almost everywhere
at each `N`. No measurability of `vhat` is needed: `TendstoInProb` bounds `μ N` of an arbitrary
set. This mirrors `thm_svd_stack_general_inner` (`SVDStack/Main.lean`) and
`thm_stacksvd_weighted_inner` (`StackSVDWeighted.lean`). -/
theorem prop_stacksvd_general_inner (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : m.stack.SingleTableLaw (∑ i, c i))
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace ((m.stack.X N ω)ᵀ * m.stack.X N ω)
      (isHermitian_transpose_mul_self (m.stack.X N ω)))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, m.stack.v N⟫_ℝ ^ 2)
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  intro ε hε
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.prop_stacksvd_general c law ε hε) (fun _ => zero_le) (fun N => ?_)
  refine measure_mono_ae ?_
  filter_upwards [law.topSimple N] with ω hω hmem'
  have heq : overlap (m.stack.X N ω) (m.stack.v N) = ⟪vhat N ω, m.stack.v N⟫_ℝ ^ 2 :=
    overlap_eq_inner_sq _ _ hω (hmem N ω) (hnorm N ω)
  have hfin : ε ≤ |overlap (m.stack.X N ω) (m.stack.v N)
      - stackSVDLimit (fun i => (m.tbl i).θ) c| := by
    rw [heq]
    exact hmem'
  exact hfin

/-- The Layer 2 form of `prop_stacksvd_general_inner`: the paper's inner product, from the
proportional regime of each table and the joint Gaussian law alone. -/
theorem prop_stacksvd_general_inner_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace ((m.stack.X N ω)ᵀ * m.stack.X N ω)
      (isHermitian_transpose_mul_self (m.stack.X N ω)))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, m.stack.v N⟫_ℝ ^ 2)
      (stackSVDLimit (fun i => (m.tbl i).θ) c) := by
  have hne : Nonempty (Fin M) := ⟨0⟩
  have hcpos : 0 < ∑ i, c i := Finset.sum_pos (fun i _ => hc i) Finset.univ_nonempty
  exact m.prop_stacksvd_general_inner c
    (SpikedModel.singleTableLaw_of_gaussian hcpos m.stack (m.stack_regime c hreg)
      (m.stack_law hG)) vhat hmem hnorm

omit [NeZero M] in
/-- `thm:simple_thm1`, stacksvd half, with every random matrix theory hypothesis discharged.
The tables share signal `θ₀` and regime `c₀`; below the recovery threshold `c₀ < M θ₀⁴` the
limit is the paper's closed form `1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²)`, and `0` otherwise. The proof is
`prop_stacksvd_general_gaussian` at constant `θ` and `c`, rewritten with the closed form
`Scalars.simple_thm1_stacksvd` (F32, 2026-09-08 drops the local `stackSVDLimit_const_eq`). -/
theorem thm_simple_thm1_stacksvd_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (if c₀ < (M : ℝ) * θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) else 0) := by
  have hM : 0 < M := NeZero.pos M
  have h := m.prop_stacksvd_general_gaussian (fun _ => c₀) (fun _ => hc) hreg hG
  rw [show (fun i => (m.tbl i).θ) = (fun _ : Fin M => θ₀) from funext hθ,
    Scalars.simple_thm1_stacksvd hM hc] at h
  exact h


end Stack

end MultiTableModel

end StackedSVD
