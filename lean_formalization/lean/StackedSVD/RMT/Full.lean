/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.TailShift
import StackedSVD.RMT.R6

/-!
# Item Full: the Layer 2 target, both regimes

`SpikedModel.singleTableLaw_of_gaussian` (`prop:single_table` for a Gaussian spiked model in
the proportional regime) used to sit in `RMT.lean` with a `sorry`. It moved here because its
proof needs `RMT/TailShift.lean` and `RMT/R6.lean`, and both of those import `RMT.lean`.
`RMT.lean` keeps the statement of `SingleTableLaw` and the two-line pointer to this file.

The theorem is one `by_cases` on `c < m.θ ^ 4`.

* Supercritical (`c < θ⁴`): `singleTableLaw_of_gaussian_supercritical'` (`RMT/TailShift.lean`),
  which is milestone L2a of `RMT/Sup.lean` after the index shift of decision D11.
* Subcritical (`θ⁴ ≤ c`): the four fields one by one. `align` is item R6'
  (`SpikedModel.align_tendstoInProb_of_subcritical`), `lamMax` is item R3⁻ through
  `lamMax_tendstoInProb_of_subcritical'`, `delocUniform` is item Sym
  (`delocUniform_of_gaussian`, both regimes, `d N → ∞` only), and `topSimple` is item S
  (`singleTableLaw_topSimple_of_gaussian`, every `N`, no threshold). The two `ResolventLimits`
  consumers take the same interface, built once by `resolventLimits_of_gaussian'`.

Every hypothesis of the target is used: `hc` and `hreg` reach both branches, `hG` reaches all
four fields, and the case split supplies `hθ` or `hθ'`. No `hn2 : ∀ N, 2 ≤ n N` appears; the
primed theorems of `RMT/TailShift.lean` remove it by a tail shift.

The doc-comment of the theorem is the one that stood in `RMT.lean`, kept word for word so
that the statement diff is empty. Its sentence "Tracked in `SORRIES.md`" no longer applies:
the `sorry` is gone and the `docs/SORRIES.md` row went with it.

STATUS 2026-08-31: see `notes/archive/agent_reports/proof_full.md`. This file has no `sorry` and no
`axiom`. It inherited a `sorryAx` from `RMT/R1.lean` until 2026-08-30 13:00;
`RMT/SteinStep.lean` closed that, so the theorem now depends on the three standard axioms
only.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- Layer 2 target: the single-table law holds under Gaussian noise in the proportional
regime. Proved 2026-08-30 from items R0, MP, MP7, R4, R4C, S, Sym, R3, R3⁻, T, R1 (with
`ResolvDeriv`), R2, R5, R6' and the tail shift of `notes/archive/rmt_roadmap.md`; the case split is
on `c < θ ^ 4`. -/
theorem singleTableLaw_of_gaussian [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) :
    m.SingleTableLaw c := by
  by_cases hθ : c < m.θ ^ 4
  · exact singleTableLaw_of_gaussian_supercritical' hc m hreg hG hθ
  · have hθ' : m.θ ^ 4 ≤ c := not_lt.mp hθ
    have H : m.ResolventLimits c := resolventLimits_of_gaussian' hc m hreg hG
    exact ⟨m.align_tendstoInProb_of_subcritical H hc hreg hG hθ',
      delocUniform_of_gaussian m hG hreg.2.1,
      m.lamMax_tendstoInProb_of_subcritical' H hc hreg hG hθ',
      singleTableLaw_topSimple_of_gaussian m hG⟩

end SpikedModel

end StackedSVD
