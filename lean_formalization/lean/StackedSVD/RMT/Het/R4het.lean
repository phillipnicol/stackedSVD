/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Duality

/-!
# Item H3: the heteroscedastic resolvent interface, and the H4, H11, H12 targets

Task H3 of `notes/archive/plan_heterolaw_A.md` (section 3.1 and the H3, H4, H11, H12 rows of
section 4). `ResolventLimitsHet` is the heteroscedastic twin of `SpikedModel.ResolventLimits`
(`RMT/R5.lean:311`): the same seven fields, on the `n`-side objects `W₀'` and `q` of
`RMT/Het/Split.lean`, at every real `z > b`. The deterministic layer it feeds is `RMT/R4.lean`
unchanged, because the column split is a rank-one positive update (plan section 3.2).

**Scalars are parameters here.** Item H1 (`RMT/Het/MPhet.lean`, written by another agent)
owns `zfun`, `sStar`, `bHet`, `sPhys`, `rhoHet` and `bSF`. Every scalar of the limit
theory enters as an explicit parameter; the substitution (made in `R5het.lean`) is:

| Parameter | TODO(H1): substitute |
|---|---|
| `b` | `Scalars.bSF w c` (proved edge) or `Scalars.bHet c w` (hypothesis `HeteroEdge`) |
| `rho` | `Scalars.rhoHet (fun i => (m.tbl i).θ) c w` |
| `Phi z` | `∑ i, θ i ^ 2 * w i ^ 2 * gHet c w i z` |
| `Psi z` | `∑ i, c i * w i ^ 2 * gHet c w i z` |
| `Phi2 z`, `Psi2 z` | the `z`-derivatives of `Phi` and `Psi` |

`HetScalarFacts` bundles the seven scalar statements that the analytic core needs about
them. Item H1 landed (`RMT/Het/MPhet.lean`); `MPhet.hetScalarFacts_of_assumption4` in
`RMT/Het/R5het.lean` proves every field for the closed forms `Phihet`, `Psihet`,
`PhihetDeriv`, `PsihetDeriv` at `rho = rhoHet` and any `b` with `bHet ≤ b < rhoHet`. The
parameters stay so that the analytic core is model-free, as `R5.lean` is.

Paper: `main_paper.tex` lines 1409 to 1440 (`thm:stacksvd_weighted`, `eq:assumption4`).

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/R4het.lean` exit 0; the interface
only, no `sorry`. The H4 targets are in `R5het.lean`; the H11 targets
`heteroLaw_of_gaussian_margin` and `heteroLaw_of_gaussian` are in `RMT/Het/Sup.lean`; the
H12 target `align_tendstoInProb_het_subcritical` is proved in `RMT/Het/R6het.lean`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD
namespace MultiTableModel

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### The interface (plan section 3.1) -/

/-- **(H1) and (H2) for the column split**, in the seven-field shape of
`SpikedModel.ResolventLimits`. `W₀'` is `m.W0het`, the deterministic vector `ũ₀ = m.u0Het`
plays `v`, and the Gaussian vector `Σ^{1/2} e = m.SigmaHalf *ᵥ m.eHet` plays `g`.

The limits are `Phi`, `Psi`, `0`, `Phi2`, `Psi2`, `0`. They are
`Φ(z) = ∑ θ_i² w_i² g_i(z)`, `Ψ(z) = ∑ c_i w_i² g_i(z)` and their derivatives
(`MPhet.Phihet`, `MPhet.Psihet` and the two derivatives of `RMT/Het/MPhet.lean`). The
aspect ratios `c` are not a parameter of this structure: they enter only through `Phi`,
`Psi`, `Phi2`, `Psi2` and through `b` (cleanup wave 4, audit finding 5). -/
structure ResolventLimitsHet (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (b : ℝ)
    (Phi Psi Phi2 Psi2 : ℝ → ℝ) : Prop where
  /-- (H1), the edge bound: item H8 at `b = bSF`, hypothesis `HeteroEdge` at `b = bHet`. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ b + ε}) atTop (𝓝 1)
  /-- (H2), `ũ₀ᵀ G₀' ũ₀ → Φ`. -/
  uu : ∀ z > b, TendstoInProb μ
    (fun N ω => qform (m.W0het w N ω) z (m.u0Het w N)) (Phi z)
  /-- (H2), `eᵀ Σ^{1/2} G₀' Σ^{1/2} e → Ψ`. -/
  ee : ∀ z > b, TendstoInProb μ
    (fun N ω => qform (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω)) (Psi z)
  /-- (H2), the cross form `→ 0`. -/
  ue : ∀ z > b, TendstoInProb μ
    (fun N ω => cform (m.W0het w N ω) z (m.u0Het w N) (m.SigmaHalf w N *ᵥ m.eHet N ω)) 0
  /-- (H2), `ũ₀ᵀ G₀'² ũ₀ → Φ'`. -/
  uu2 : ∀ z > b, TendstoInProb μ
    (fun N ω => qform2 (m.W0het w N ω) z (m.u0Het w N)) (Phi2 z)
  /-- (H2), `eᵀ Σ^{1/2} G₀'² Σ^{1/2} e → Ψ'`. -/
  ee2 : ∀ z > b, TendstoInProb μ
    (fun N ω => qform2 (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω)) (Psi2 z)
  /-- (H2), the squared cross form `→ 0`. -/
  ue2 : ∀ z > b, TendstoInProb μ
    (fun N ω => cform2 (m.W0het w N ω) z (m.u0Het w N) (m.SigmaHalf w N *ᵥ m.eHet N ω)) 0

/-- **The heteroscedastic edge**, as a hypothesis structure for a general noise law. The
field says `λ_max(W₀') ≤ bHet + ε` with probability tending to `1`, for every `ε > 0`.
`bHet` is a parameter until H1 lands.

Plan section 3.6 named this the single residual black box. Campaign E proves it for
Gaussian noise (`MultiTableModel.heteroEdge_of_gaussian`, `RMT/Het/EdgeSharp.lean`), so
`heteroLaw_of_gaussian` (`RMT/Het/Sup.lean`) needs no hypothesis of this kind. The structure
stays as the general interface for any other noise law: Sudakov-Fernique proves the weaker
bound at `b = bSF w c` (item H8), which is why the margin theorem needs no hypothesis of
this kind either. -/
structure HeteroEdge (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (bHet : ℝ) : Prop where
  /-- `λ_max(W₀') ≤ bHet + ε` with probability tending to `1`, for every `ε > 0`. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ bHet + ε}) atTop (𝓝 1)

/-- The seven scalar statements of plan sections 2.3, 2.4 and 3.5 that item H1 proves about
the heteroscedastic branch; `MPhet.hetScalarFacts_of_assumption4` (`RMT/Het/R5het.lean`)
is the concrete instance. `neg` and `pos` say that the limiting secular function
`1 + Φ + Ψ` changes sign at `rho`; `overlap` is the identity `1/(ρ F'(ρ)) = L(w)` checked
numerically in plan section 2.4 (`0.639273` both ways, difference `1.6e-11`). -/
structure HetScalarFacts (θ c w : Fin M → ℝ) (b rho : ℝ) (Phi Psi Phi2 Psi2 : ℝ → ℝ) :
    Prop where
  /-- The outlier sits above the edge. -/
  lt : b < rho
  /-- `1 + F < 0` below the root. -/
  neg : ∀ z, b < z → z < rho → 1 + Phi z + Psi z < 0
  /-- `1 + F > 0` above the root. -/
  pos : ∀ z, rho < z → 0 < 1 + Phi z + Psi z
  /-- Continuity of `F` on `(b, ∞)`, for the bracket of R5 step 8. -/
  cont : ContinuousOn (fun z => Phi z + Psi z) (Set.Ioi b)
  /-- Continuity of `F'` on `(b, ∞)`. -/
  cont2 : ContinuousOn (fun z => Phi2 z + Psi2 z) (Set.Ioi b)
  /-- The denominator of the overlap is positive at the root. -/
  den_pos : 0 < Phi2 rho + Psi2 rho
  /-- `1/(ρ F'(ρ)) = L(w)`. -/
  overlap : 1 / (rho * (Phi2 rho + Psi2 rho)) = Scalars.Lw θ c w

/-! ### Targets of item H4

The three H4 conclusions `lamMax_tendstoInProb_het`, `tendsto_measure_topSimple_het` and
`align_tendstoInProb_het` are proved in `RMT/Het/R5het.lean` (2026-08-31), over the
interface above, with `[∀ N, IsProbabilityMeasure (μ N)]` (audit item F1). -/

/-! ### Item H12

Item H11 (`RMT/Het/Sup.lean`, 2026-08-30) proved the two theorems that used to stand here,
`heteroLaw_of_gaussian_margin` and `heteroLaw_of_gaussian`. Item H12
(`RMT/Het/R6het.lean`, 2026-08-30) proved the subcritical half of `HeteroLaw.align`,
`MultiTableModel.align_tendstoInProb_het_subcritical`, over the interface above at
`b = MPhet.bHet c w`; with it `heteroLaw_of_gaussian` covers both regimes under `HeteroEdge`
alone, and no `sorry` is left in this layer. -/

end MultiTableModel
end StackedSVD
