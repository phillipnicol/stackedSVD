/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R2het
import StackedSVD.RMT.Het.R3het
import StackedSVD.RMT.Het.Simplicity
import StackedSVD.RMT.Het.R6het
import StackedSVD.RMT.Het.EdgeSharp
import StackedSVD.StackSVD.Weighted

/-!
# Item H11: the Gaussian assembly of `HeteroLaw`

Task H11 of `notes/archive/plan_heterolaw_A.md` (sections 1, 3.6, 3.7 and the H11 row of section 4).
Nothing new is proved about random matrices here. The file plugs the six isotropic limits of
`RMT/Het/R2het.lean` (item H7), the Sudakov-Fernique edge of `RMT/Het/R3het.lean` (item H8)
and the almost sure simplicity of `RMT/Het/Simplicity.lean` (item H10) into the analytic core
of `RMT/Het/R5het.lean` (item H4), and so discharges the black box
`MultiTableModel.HeteroLaw` for Gaussian noise.

## Content

1. `resolventLimitsHet_of_edge_bound`: the seven-field interface `ResolventLimitsHet` at any
   edge bound `b ≥ bHet` whose edge event is given. Two instances follow:
   `resolventLimitsHet_of_gaussian` at `b = bSF c w` (proved, item H8) and
   `resolventLimitsHet_of_edge` at `b = bHet c w` (hypothesis `HeteroEdge`).
2. `heteroLaw_of_gaussian_margin`: `HeteroLaw` with no black box, above `eq:assumption4` and
   above the margin `bSF c w < rhoHet θ c w` of plan section 3.6. Campaign E later proves the
   exact edge (item 3), so this theorem is a weaker corollary, kept as a record of the
   interim result.
3. `heteroLaw_of_gaussian`: `HeteroLaw` in both regimes at the exact edge `bHet`
   (supercritical half from item H4, subcritical half from item H12, `RMT/Het/R6het.lean`),
   with the edge itself proved in `RMT/Het/EdgeSharp.lean` (campaign E).
4. The paper-level corollaries at the exact edge, both regimes, model hypotheses only:
   `thm_stacksvd_weighted_gaussian`, `_inner`, `_smul` (at `t • optWstack` for every
   `t ≠ 0`, the paper's `w ∝ θ_i/√(θ_i² + c_i)`), `_opt` (the optimal weights beat every
   other nonzero `w`) and `prop_dominance_gaussian` (`prop:dominance`,
   `main_paper.tex:634`: optimal stacksvd beats optimal svdstack and unweighted stacksvd).
   The margin corollaries `thm_stacksvd_weighted_gaussian_margin`, `_inner`, `_smul` are
   weaker versions of the same statements, above `bSF c w < rhoHet θ c w` instead of the
   exact edge; they are kept as a record of the interim result, not because any later
   theorem needs the weaker form.
   None of the theorems above takes `[∀ N, IsProbabilityMeasure (μ N)]`: `hG` gives it
   through `(hG N).isProbabilityMeasure` (`notes/INTERFACES.md`).
5. `stackPerfW_binary_inad_gaussian` and `stackPerfW_opt_inad_gaussian`:
   `prop:binarystacksvd_inadmissable` (`main_paper.tex:660`) for Gaussian noise, on the
   instance with every `θ_i = 1`. `heteroLaw_of_gaussian` (item 3) supplies the law they
   need, which is why these two theorems live here and not next to the rest of
   `prop:binarystacksvd_inadmissable` in `StackSVDWeighted.lean`. Unlike the theorems of
   item 4, both take `[∀ N, IsProbabilityMeasure (μ N)]` explicitly, matching the
   convention of the base theorems they wrap (`StackSVDWeighted.lean:1443,1665`).
6. A satisfiability check at `M = 2`, `θ_i = √5`, `c_i = 1`: both optimal weights are
   `√5/√6 ≠ 0`, `1 < ∑ θ_i⁴/c_i = 50`, and `bSF < 9 < 11 = rhoHet`.

## The shape of `HeteroLaw.align`, and where `h4` stays a hypothesis

`HeteroLaw.align` (`StackSVDWeighted.lean:1267`) is one statement for both regimes: the
overlap tends to `Scalars.Lw θ c w`, which is the closed form under `Scalars.Assumption4`
and `0` when `Assumption4` fails. Item H4 proves the supercritical half
(`align_tendstoInProb_het_of_assumption4` takes `h4`) at any edge bound `b` with
`bHet ≤ b < rhoHet`; item H12 (`align_tendstoInProb_het_subcritical`, `RMT/Het/R6het.lean`,
2026-08-30) proves the subcritical half at the exact edge `b = bHet` only, because its
argument needs `F'(z) → ∞` as `z ↓ b`. So the margin theorems, which know the edge only up
to `bSF ≥ bHet`, keep `h4` (at `w = optWstack θ c` it is the paper's detectability threshold
`1 < ∑_i θ_i⁴/c_i`, `Scalars.assumption4_optW_iff`), while `heteroLaw_of_gaussian` and
`thm_stacksvd_weighted_gaussian` run at the exact edge `bHet` and cover both regimes.
No field of `HeteroLaw` is weakened in either case.

Paper: `main_paper.tex` line 462 (`thm:stacksvd_weighted`), lines 1409 to 1440
(`eq:assumption4`, the proof), line 1602 (the threshold `∑_j θ_j⁴/c_j > 1`).

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/Sup.lean` exit 0; no `sorry`,
no `axiom`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace Scalars

/-- Audit item F5 (`notes/archive/audit_het_skeleton_2026-08-31.md`): the optimal weights are not
all zero as soon as one table has signal, which is what the simplicity input `hw` of
`RMT/Het/Simplicity.lean` needs. -/
theorem optWstack_ne_zero {M : ℕ} {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, θ i ≠ 0) : ∃ i, optWstack θ c i ≠ 0 := by
  obtain ⟨i, hi⟩ := hθ
  refine ⟨i, ?_⟩
  have hpos : 0 < θ i ^ 2 + c i := by
    have := hc i; positivity
  have hs : Real.sqrt (θ i ^ 2 + c i) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hpos)
  simpa [optWstack] using div_ne_zero hi hs

end Scalars

namespace MultiTableModel

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### The interface, from the six limits of item H7 -/

/-- **The assembly of `ResolventLimitsHet`** (plan section 3.1, the "For H11" paragraph of
`notes/archive/agent_reports/h7_r2het.md`). The six theorems of `RMT/Het/R2het.lean` hold at every
real `x > b` for any `b` above `bHet`, once the edge event at `b` is known.

`[∀ N, IsProbabilityMeasure (μ N)]` is not an argument of any theorem of this file: `hG`
gives it, by `(hG N).isProbabilityMeasure` (`notes/INTERFACES.md`, cleanup wave 4). -/
theorem resolventLimitsHet_of_edge_bound
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {b : ℝ} (hb : MPhet.bHet c w ≤ b)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ b + ε}) atTop (𝓝 1)) :
    m.ResolventLimitsHet w b (MPhet.Phihet (fun i => (m.tbl i).θ) c w) (MPhet.Psihet c w)
      (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w) (MPhet.PsihetDeriv c w) := by
  have : ∀ N, IsProbabilityMeasure (μ N) := fun N => (hG N).isProbabilityMeasure
  exact
  ⟨hedge, fun _z hz => m.tendstoInProb_qform_u0Het w c hc hw hreg hG hb hedge hz,
    fun _z hz => m.tendstoInProb_qform_eHet w c hc hw hreg hG hb hedge hz,
    fun _z hz => m.tendstoInProb_cform_u0Het_eHet w c hc hw hreg hG hedge hz,
    fun _z hz => m.tendstoInProb_qform2_u0Het w c hc hw hreg hG hb hedge hz,
    fun _z hz => m.tendstoInProb_qform2_eHet w c hc hw hreg hG hb hedge hz,
    fun _z hz => m.tendstoInProb_cform2_u0Het_eHet w c hc hw hreg hG hedge hz⟩

/-- **The interface with no black box**, at the Sudakov-Fernique edge `bSF c w` of plan
section 3.6. The edge event is `tendsto_measure_lamMax_W0het_le_bSF` (item H8) and
`bHet ≤ bSF` is `MPhet.bHet_le_bSF`. -/
theorem resolventLimitsHet_of_gaussian
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.ResolventLimitsHet w (MPhet.bSF c w)
      (MPhet.Phihet (fun i => (m.tbl i).θ) c w) (MPhet.Psihet c w)
      (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w) (MPhet.PsihetDeriv c w) := by
  have : ∀ N, IsProbabilityMeasure (μ N) := fun N => (hG N).isProbabilityMeasure
  exact m.resolventLimitsHet_of_edge_bound w c hc hw hreg hG (MPhet.bHet_le_bSF hc hw)
    (m.tendsto_measure_lamMax_W0het_le_bSF w c hreg hG)

/-- **The interface at the exact edge**, from the single black box `HeteroEdge` of plan
section 3.6 (`RMT/Het/R4het.lean`), which is the `edge` field at `b = MPhet.bHet c w`. -/
theorem resolventLimitsHet_of_edge
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hedge : m.HeteroEdge w c (MPhet.bHet c w)) :
    m.ResolventLimitsHet w (MPhet.bHet c w)
      (MPhet.Phihet (fun i => (m.tbl i).θ) c w) (MPhet.Psihet c w)
      (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w) (MPhet.PsihetDeriv c w) :=
  m.resolventLimitsHet_of_edge_bound w c hc hw hreg hG le_rfl hedge.edge

/-! ### `HeteroLaw` for Gaussian noise -/

/-- **Item H11, the unconditional half.** Above `eq:assumption4` and above the
Sudakov-Fernique margin `bSF c w < rhoHet θ c w`, the Gaussian model satisfies the black box
`HeteroLaw`, with no random matrix theory input left over.

`align` is `align_tendstoInProb_het_of_assumption4` at `b = bSF c w` (item H4 over items H5
to H8); `topSimple` is `heteroLaw_topSimple_of_gaussian` (item H10), which needs `w ≠ 0`.
That is not a hypothesis: `h4.1` is the existence of the outlier root, and
`MPhet.exists_w_ne_zero_of_root` reads `∃ i, w i ≠ 0` off it (cleanup wave 4, audit
finding 2).

The margin is not vacuous and not the paper's condition: it is implied by no hypothesis of
the paper. The coverage figure of plan section 3.6 (82 % of supercritical draws,
`notes/archive/audit_het_skeleton_2026-08-31.md` section 3, 36286 draws, seed 20260834) is measured
at `w = Scalars.optWstack θ c` only; this theorem is stated for a general `w`, where the
margin can fail more often. Campaign E later proves `HeteroEdge` for Gaussian noise
(`heteroEdge_of_gaussian`, `RMT/Het/EdgeSharp.lean`), so `heteroLaw_of_gaussian` below covers
the band `bHet < bSF` and the subcritical regime with no margin hypothesis at all. This
theorem is kept as a weaker corollary, a record of the interim result proved before item
H12 closed the subcritical case.

Paper: `main_paper.tex` line 462 (`thm:stacksvd_weighted`), line 1409 (`eq:assumption4`). -/
theorem heteroLaw_of_gaussian_margin
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (h4 : Scalars.Assumption4 (fun i => (m.tbl i).θ) c w)
    (hmargin : MPhet.bSF c w < MPhet.rhoHet (fun i => (m.tbl i).θ) c w) :
    m.HeteroLaw w c := by
  have : ∀ N, IsProbabilityMeasure (μ N) := fun N => (hG N).isProbabilityMeasure
  have hw : ∃ i, w i ≠ 0 := MPhet.exists_w_ne_zero_of_root h4.1
  exact
  { align := align_tendstoInProb_het_of_assumption4 hc h4 (MPhet.bHet_le_bSF hc hw) hmargin
      (m.resolventLimitsHet_of_gaussian w c hc hw hreg hG)
    topSimple := m.heteroLaw_topSimple_of_gaussian w hG hw }

omit [NeZero M] in
/-- **Item H11 under the exact edge, both regimes.** `HeteroLaw` with `bSF` replaced by the
true edge `bHet`, with no black box and no condition on the parameters: above
`eq:assumption4` the limit is the closed form (item H4, `MPhet.bHet_lt_rhoHet` gives the
margin for free), below it the limit is `0` (item H12,
`align_tendstoInProb_het_subcritical` in `RMT/Het/R6het.lean`).

The edge `HeteroEdge` is no longer a hypothesis. Campaign E proves it for Gaussian noise by
the sharp Sudakov-Fernique bound (`MultiTableModel.heteroEdge_of_gaussian`,
`RMT/Het/EdgeSharp.lean`).

`hw : ∃ i, w i ≠ 0` stays a hypothesis here, unlike in `heteroLaw_of_gaussian_margin`: this
theorem has no `h4` to read it off, and at `w = 0` the statement is false, because the
weighted stack is the zero matrix and `topSimple` fails as soon as `2 ≤ d N`. -/
theorem heteroLaw_of_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroLaw w c := by
  have : ∀ N, IsProbabilityMeasure (μ N) := fun N => (hG N).isProbabilityMeasure
  have hedge : m.HeteroEdge w c (MPhet.bHet c w) :=
    m.heteroEdge_of_gaussian w c hc hw hreg hG
  refine ⟨?_, m.heteroLaw_topSimple_of_gaussian w hG hw⟩
  by_cases h4 : Scalars.Assumption4 (fun i => (m.tbl i).θ) c w
  · exact align_tendstoInProb_het_of_assumption4 hc h4 le_rfl (MPhet.bHet_lt_rhoHet hc h4)
      (m.resolventLimitsHet_of_edge w c hc hw hreg hG hedge)
  · exact m.align_tendstoInProb_het_subcritical w c hc hw hG
      (m.resolventLimitsHet_of_edge w c hc hw hreg hG hedge) h4

/-! ### `thm:stacksvd_weighted` for Gaussian noise (`main_paper.tex:462`) -/

/-- **`thm:stacksvd_weighted` with no black box** (`main_paper.tex:462`): for Gaussian noise,
above the detectability threshold `1 < ∑_i θ_i⁴/c_i` (`main_paper.tex:1602`) and above the
Sudakov-Fernique margin, the performance of stacksvd at the paper's optimal weights
`w_i⋆ ∝ θ_i/√(θ_i² + c_i)` tends in probability to `γ_opt`, the unique root in `(0,1)` of
`∑_i θ_i⁴(1-x)/(c_i + xθ_i²) = 1`.

The threshold enters as `Scalars.assumption4_optW_iff`, so the hypothesis is the paper's own
condition at the optimal weights. Every hypothesis is a hypothesis on the model.

`∃ i, θ_i ≠ 0` is not a hypothesis: `hdet` gives it, because every `θ_i = 0` makes the sum
`0` (cleanup wave 4, audit finding 1). -/
theorem thm_stacksvd_weighted_gaussian_margin
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) (hdet : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (hmargin : MPhet.bSF c (Scalars.optWstack (fun i => (m.tbl i).θ) c) <
      MPhet.rhoHet (fun i => (m.tbl i).θ) c
        (Scalars.optWstack (fun i => (m.tbl i).θ) c)) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  have hθ : ∃ i, (m.tbl i).θ ≠ 0 := by
    by_contra hcon
    have hz : ∑ i, (m.tbl i).θ ^ 4 / c i = 0 :=
      Finset.sum_eq_zero fun i _ => by
        rw [not_not.mp (fun h : (m.tbl i).θ ≠ 0 => hcon ⟨i, h⟩)]; norm_num
    rw [hz] at hdet
    exact absurd hdet (by norm_num)
  exact m.thm_stacksvd_weighted c hc hθ
    (m.heteroLaw_of_gaussian_margin _ c hc hreg hG
      ((Scalars.assumption4_optW_iff hc).mpr hdet) hmargin)

/-- The paper's inner-product form of the corollary above: for any selected unit top
eigenvector `v̂` of the weighted stack Gram matrix, `(v̂ᵀ v)² → γ_opt`
(`main_paper.tex:462`). The simplicity field of `HeteroLaw` is what makes the projector form
and the inner-product form agree. `∃ i, θ_i ≠ 0` follows from `hdet`, as above. -/
theorem thm_stacksvd_weighted_gaussian_margin_inner
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) (hdet : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (hmargin : MPhet.bSF c (Scalars.optWstack (fun i => (m.tbl i).θ) c) <
      MPhet.rhoHet (fun i => (m.tbl i).θ) c
        (Scalars.optWstack (fun i => (m.tbl i).θ) c))
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈
      topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
        (m.isHermitian_stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  have hθ : ∃ i, (m.tbl i).θ ≠ 0 := by
    by_contra hcon
    have hz : ∑ i, (m.tbl i).θ ^ 4 / c i = 0 :=
      Finset.sum_eq_zero fun i _ => by
        rw [not_not.mp (fun h : (m.tbl i).θ ≠ 0 => hcon ⟨i, h⟩)]; norm_num
    rw [hz] at hdet
    exact absurd hdet (by norm_num)
  exact m.thm_stacksvd_weighted_inner c hc hθ
    (m.heteroLaw_of_gaussian_margin _ c hc hreg hG
      ((Scalars.assumption4_optW_iff hc).mpr hdet) hmargin) vhat hmem hnorm

omit [NeZero M] in
/-- **`thm:stacksvd_weighted` under the exact edge, both regimes** (`main_paper.tex:462`):
with no black box, the performance at the optimal weights tends to `γ_opt` above the
threshold `1 < ∑_i θ_i⁴/c_i` and to `0` below it (`Scalars.stackSVDLimitW` carries both
branches), with no condition on the parameters. The exact edge comes from
`MultiTableModel.heteroEdge_of_gaussian` (campaign E, `RMT/Het/EdgeSharp.lean`).

`hθ : ∃ i, θ_i ≠ 0` stays a hypothesis here: this theorem has no detectability threshold to
read it off, and below the threshold the conclusion at `θ = 0` is false (audit finding 1
covers the two margin corollaries only). -/
theorem thm_stacksvd_weighted_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) :=
  m.thm_stacksvd_weighted c hc hθ
    (m.heteroLaw_of_gaussian _ c hc (Scalars.optWstack_ne_zero hc hθ) hreg hG)

/-- The inner-product form of `thm_stacksvd_weighted_gaussian`, under the exact edge and
both regimes (`main_paper.tex:462`). Mirrors `thm_stacksvd_weighted_gaussian_margin_inner`,
with `heteroLaw_of_gaussian` in place of the margin route: every random matrix theory
hypothesis is discharged for Gaussian noise. -/
theorem thm_stacksvd_weighted_gaussian_inner
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈
      topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
        (m.isHermitian_stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) :=
  m.thm_stacksvd_weighted_inner c hc hθ
    (m.heteroLaw_of_gaussian _ c hc (Scalars.optWstack_ne_zero hc hθ) hreg hG) vhat hmem hnorm

/-! ### The optimal weights are fixed up to a scale (`w ∝ θ_i/√(θ_i² + c_i)`)

`main_paper.tex:462` writes the optimal weights as a proportionality. The performance of
stacksvd sees only the direction of `w`, because `t • w` multiplies the whole weighted stack
by `t` and `Spectral.overlap_smul` is scale invariant. The two lemmas below sit here, and
not in `StackSVDWeighted.lean`, so that no olean outside `RMT/Het/` is invalidated; move
them up when that file is next edited.
-/

/-- The weighted stack at `t • w` is `t` times the weighted stack at `w`. The constant-weight
case is `stackW_const_X` (`StackSVDWeighted.lean`). -/
theorem stackW_smul_X (m : MultiTableModel μ M n d) (t : ℝ) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.stackW (t • w)).X N ω = t • (m.stackW w).X N ω := by
  ext r k
  obtain ⟨p, rfl⟩ : ∃ p, finSigmaFinEquiv p = r :=
    ⟨finSigmaFinEquiv.symm r, Equiv.apply_symm_apply _ _⟩
  obtain ⟨i, j⟩ := p
  rw [m.stackW_X (t • w) N ω i j k, Matrix.smul_apply, m.stackW_X w N ω i j k]
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

/-- **`stackPerfW` is scale invariant**: `(v̂ᵀ v)²` at the weighting `t • w` equals `(v̂ᵀ v)²`
at `w`, for every `t ≠ 0`. -/
theorem stackPerfW_smul (m : MultiTableModel μ M n d) {t : ℝ} (ht : t ≠ 0) (w : Fin M → ℝ)
    (N : ℕ) (ω : Ω N) : m.stackPerfW (t • w) N ω = m.stackPerfW w N ω := by
  change overlap ((m.stackW (t • w)).X N ω) ((m.tbl 0).v N)
      = overlap ((m.stackW w).X N ω) ((m.tbl 0).v N)
  rw [m.stackW_smul_X t w N ω, overlap_smul ht]

/-- **`thm:stacksvd_weighted` at every weighting proportional to `w⋆`**
(`main_paper.tex:462`, `w_i⋆ ∝ θ_i/√(θ_i² + c_i)`): the conclusion of
`thm_stacksvd_weighted_gaussian_margin` holds at `t • w⋆` for every `t ≠ 0`. The scalar
limit does not move, because `Scalars.stackSVDLimitW` does not depend on `w`; the random
side does not move either, by `stackPerfW_smul`. -/
theorem thm_stacksvd_weighted_gaussian_margin_smul
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) (hdet : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (hmargin : MPhet.bSF c (Scalars.optWstack (fun i => (m.tbl i).θ) c) <
      MPhet.rhoHet (fun i => (m.tbl i).θ) c
        (Scalars.optWstack (fun i => (m.tbl i).θ) c))
    {t : ℝ} (ht : t ≠ 0) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (t • Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  have h := m.thm_stacksvd_weighted_gaussian_margin c hc hreg hG hdet hmargin
  simpa only [m.stackPerfW_smul ht] using h

/-- **`thm_stacksvd_weighted_gaussian` at every weighting proportional to `w⋆`**
(`main_paper.tex:462`, `w_i⋆ ∝ θ_i/√(θ_i² + c_i)`), under the exact edge and both regimes.
Mirrors `thm_stacksvd_weighted_gaussian_margin_smul`: every random matrix theory hypothesis
is discharged for Gaussian noise. -/
theorem thm_stacksvd_weighted_gaussian_smul
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {t : ℝ} (ht : t ≠ 0) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (t • Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  have h := m.thm_stacksvd_weighted_gaussian c hc hθ hreg hG
  simpa only [m.stackPerfW_smul ht] using h

/-! ### The optimal weights against every other weighting, and against svdstack and
unweighted stacksvd, under the exact edge -/

omit [NeZero M] in
/-- **`thm:stacksvd_weighted`, bundled with `L_le_opt`** (`main_paper.tex:462`): the
optimal weights `optWstack` attain `Scalars.stackSVDLimitW`, and every other nonzero
weighting `w` converges to a value `Scalars.Lw θ c w` that is at most `stackSVDLimitW`. The
first half is `thm_stacksvd_weighted_gaussian`; the second is `heteroLaw_of_gaussian.align`
for the convergence and `Scalars.L_le_opt` for the simplex optimization, so every random
matrix theory hypothesis is discharged for Gaussian noise. -/
theorem thm_stacksvd_weighted_gaussian_opt [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ ∀ w : Fin M → ℝ, (∃ i, w i ≠ 0) →
        TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
          (Scalars.Lw (fun i => (m.tbl i).θ) c w)
        ∧ Scalars.Lw (fun i => (m.tbl i).θ) c w
            ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c :=
  ⟨m.thm_stacksvd_weighted_gaussian c hc hθ hreg hG,
    fun w hw => ⟨(m.heteroLaw_of_gaussian w c hc hw hreg hG).align, Scalars.L_le_opt hc w⟩⟩

omit [NeZero M] in
/-- **`prop:dominance`** (`main_paper.tex:634`), with the Gaussian estimator convergence at
the optimal weights attached: optimally weighted stacksvd tends to `Scalars.stackSVDLimitW`
in probability, and that limit dominates both optimally weighted svdstack
(`Scalars.svdstackOpt_le_stackSVDLimitW`) and unweighted stacksvd
(`Scalars.stackSVDLimit_le_stackSVDLimitW`). Every random matrix theory hypothesis is
discharged for Gaussian noise. -/
theorem prop_dominance_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c :=
  ⟨m.thm_stacksvd_weighted_gaussian c hc hθ hreg hG,
    Scalars.svdstackOpt_le_stackSVDLimitW hc, Scalars.stackSVDLimit_le_stackSVDLimitW hc⟩

/-! ### `prop:binarystacksvd_inadmissable`, Gaussian noise (`main_paper.tex:660`)

The scalar content and the model half over a generic law are in `StackSVDWeighted.lean`.
These two Gaussian corollaries need `heteroLaw_of_gaussian`, so they land here rather than
in `StackSVDWeighted.lean`, which `RMT/Het/Sup.lean` imports (`grep '^import'` on both
files; see `notes/archive/agent_reports/facade_H.md`). -/

omit [NeZero M] in
/-- On the instance of `prop:binarystacksvd_inadmissable`, every binary weighting of
Gaussian stacksvd tends to `0` in probability. `S` is any nonempty subset, so this covers
the optimal binary weighting. Route: `stackPerfW_binary_tendsto_gaussian` at
`c = Scalars.inadC M`, then `Scalars.inad_binaryStackSVDLimit_eq_zero`. -/
theorem stackPerfW_binary_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1) (S : Finset (Fin M))
    [NeZero S.card] (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0 := by
  have hfun : (fun i => (m.tbl i).θ) = Scalars.inadTheta M := funext fun i => hθ i
  have h := m.stackPerfW_binary_tendsto_gaussian S (Scalars.inadC M) (Scalars.inadC_pos M)
    hreg hG
  rw [hfun, Scalars.inad_binaryStackSVDLimit_eq_zero S] at h
  exact h

omit [NeZero M] in
/-- On the same instance, optimally weighted Gaussian stacksvd tends to a value above
`1 - ε`, as soon as `M ≥ e^{-γ} exp(2/ε)`. The Gaussian law of `HeteroLaw` at `optWstack` is
`heteroLaw_of_gaussian`, discharged by `w ≠ 0` (`Scalars.optWstack_ne_zero`, at table `0`
where `θ = 1 ≠ 0`), `hreg` and `hG`. -/
theorem stackPerfW_opt_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
        (fun N ω => m.stackPerfW
          (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
        (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M) := by
  have hw : ∃ i, Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M) i ≠ 0 :=
    Scalars.optWstack_ne_zero (Scalars.inadC_pos M) ⟨0, by rw [hθ 0]; exact one_ne_zero⟩
  have law : m.HeteroLaw (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M))
      (Scalars.inadC M) :=
    m.heteroLaw_of_gaussian _ _ (Scalars.inadC_pos M) hw hreg hG
  exact m.stackPerfW_opt_inad hθ hε hε1 hM law

end MultiTableModel

/-! ### Satisfiability of the margin (deliverable 5 of item H11)

`M = 2`, `θ_1 = θ_2 = √5`, `c_1 = c_2 = 1`. Then `w⋆_i = √5/√6` for both tables,
`∑ θ_i⁴/c_i = 50 > 1`, `γ₁ = 55/6`, `rhoHet = 11` and `bSF = (√(5/6) + √(5/3))² < 9`, so
every hypothesis of `thm_stacksvd_weighted_gaussian_margin` holds at once. Numerically
(a session script, not kept; no random input for the scalars): `bSF = bHet = 4.857023`, so
the weights are equal across the tables and Sudakov-Fernique is sharp here;
`L(w⋆) = 49/55 = 0.890909`. Simulation `d = 400`, `n_i = 400`, seed 20260830, 6 draws:
`λ_max(W₀') = 4.803 ± 0.050`, below the edge. -/

namespace HetSupExample

open Scalars

/-- The two signal strengths of the example. -/
noncomputable def thetaEx : Fin 2 → ℝ := fun _ => Real.sqrt 5

/-- The two aspect ratios of the example. -/
def cEx : Fin 2 → ℝ := fun _ => 1

theorem cEx_pos : ∀ i, 0 < cEx i := fun _ => one_pos

theorem thetaEx_ne_zero : ∀ i, thetaEx i ≠ 0 := fun _ =>
  ne_of_gt (Real.sqrt_pos.mpr (by norm_num))

theorem thetaEx_sq (i : Fin 2) : thetaEx i ^ 2 = 5 :=
  Real.sq_sqrt (by norm_num)

/-- The optimal weights are the constant `√5/√6`, and both are nonzero. -/
theorem optW_eq : optWstack thetaEx cEx = fun _ => Real.sqrt 5 / Real.sqrt 6 := by
  funext i
  have h : thetaEx i ^ 2 + cEx i = 6 := by rw [thetaEx_sq]; norm_num [cEx]
  change thetaEx i / Real.sqrt (thetaEx i ^ 2 + cEx i) = Real.sqrt 5 / Real.sqrt 6
  rw [h]
  rfl

theorem optW_sq (i : Fin 2) : optWstack thetaEx cEx i ^ 2 = 5 / 6 := by
  rw [optW_eq]
  rw [div_pow, Real.sq_sqrt (by norm_num : (5:ℝ) ≥ 0).le,
    Real.sq_sqrt (by norm_num : (6:ℝ) ≥ 0).le]

/-- The detectability threshold of `main_paper.tex:1602`: `∑ θ_i⁴/c_i = 50 > 1`. -/
theorem det_ex : 1 < ∑ i, thetaEx i ^ 4 / cEx i := by
  have h : ∀ i : Fin 2, thetaEx i ^ 4 / cEx i = 25 := by
    intro i
    have : thetaEx i ^ 4 = (thetaEx i ^ 2) ^ 2 := by ring
    rw [this, thetaEx_sq]
    norm_num [cEx]
  rw [Finset.sum_congr rfl fun i _ => h i]
  norm_num

/-- `γ₁ = t²(1 + ∑ θ_i²) = (5/6)·11`. -/
theorem gammaTop_ex : gammaTop thetaEx (optWstack thetaEx cEx) = 55 / 6 := by
  have hne : Real.sqrt 5 / Real.sqrt 6 ≠ 0 :=
    div_ne_zero (ne_of_gt (Real.sqrt_pos.mpr (by norm_num)))
      (ne_of_gt (Real.sqrt_pos.mpr (by norm_num)))
  have hw : optWstack thetaEx cEx
      = fun i => (Real.sqrt 5 / Real.sqrt 6) * (1 : Fin 2 → ℝ) i := by
    rw [optW_eq]; funext i; simp
  have hsum : ∑ i, thetaEx i ^ 2 = 10 := by
    rw [Fin.sum_univ_two, thetaEx_sq, thetaEx_sq]; norm_num
  have ht : (Real.sqrt 5 / Real.sqrt 6) ^ 2 = 5 / 6 := by
    rw [div_pow, Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 5),
      Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 6)]
  rw [hw, gammaTop_smul, gammaTop_one ⟨0, thetaEx_ne_zero 0⟩, hsum, ht]
  norm_num

/-- `rhoHet = γ₁(1 + ∑ c_i w_i²/(γ₁ - w_i²)) = 11`. -/
theorem rhoHet_ex : MPhet.rhoHet thetaEx cEx (optWstack thetaEx cEx) = 11 := by
  rw [MPhet.rhoHet, gammaTop_ex]
  rw [Finset.sum_congr rfl fun i (_ : i ∈ Finset.univ) => by rw [optW_sq i]]
  norm_num [cEx]

/-- `bSF = (√(5/6) + √(5/3))² < 9`, by the crude bounds `√(5/6) < 1` and `√(5/3) < 2`. -/
theorem bSF_ex_lt : MPhet.bSF cEx (optWstack thetaEx cEx) < 9 := by
  have hmax : wSqMax (optWstack thetaEx cEx) = 5 / 6 := by
    have : (fun i => optWstack thetaEx cEx i ^ 2) = fun _ : Fin 2 => (5 : ℝ) / 6 :=
      funext fun i => optW_sq i
    rw [wSqMax, this, ciSup_const]
  have hsum : ∑ i, cEx i * optWstack thetaEx cEx i ^ 2 = 5 / 3 := by
    rw [Fin.sum_univ_two, optW_sq, optW_sq]
    norm_num [cEx]
  have h1 : Real.sqrt (5 / 6) < 1 := by
    rw [show (1 : ℝ) = Real.sqrt 1 by simp]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have h2 : Real.sqrt (5 / 3) < 2 := by
    rw [show (2 : ℝ) = Real.sqrt 4 by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have hn1 : 0 ≤ Real.sqrt (5 / 6) := Real.sqrt_nonneg _
  have hn2 : 0 ≤ Real.sqrt (5 / 3) := Real.sqrt_nonneg _
  rw [MPhet.bSF, hmax, hsum]
  nlinarith

/-- **The margin holds** at the example: `bSF < 9 < 11 = rhoHet`. -/
theorem margin_ex : MPhet.bSF cEx (optWstack thetaEx cEx) <
    MPhet.rhoHet thetaEx cEx (optWstack thetaEx cEx) := by
  rw [rhoHet_ex]
  exact lt_trans bSF_ex_lt (by norm_num)

open MultiTableModel in
/-- Every hypothesis of `thm_stacksvd_weighted_gaussian_margin` is satisfiable at once: on
any Gaussian two-table model in the regime `c = (1, 1)` with `θ = (√5, √5)`, the theorem
applies and gives the conclusion of `thm:stacksvd_weighted` (`main_paper.tex:462`). -/
example {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
    {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ}
    (m : MultiTableModel μ 2 n d) (hθv : ∀ i, (m.tbl i).θ = Real.sqrt 5)
    (hreg : ∀ i, (m.tbl i).Regime (cEx i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (optWstack (fun i => (m.tbl i).θ) cEx) N ω)
      (stackSVDLimitW (fun i => (m.tbl i).θ) cEx) := by
  have hfun : (fun i => (m.tbl i).θ) = thetaEx := funext hθv
  refine m.thm_stacksvd_weighted_gaussian_margin cEx cEx_pos hreg hG ?_ ?_
  · have h := det_ex
    simpa [thetaEx, hθv] using h
  · rw [hfun]; exact margin_ex

end HetSupExample

end StackedSVD
