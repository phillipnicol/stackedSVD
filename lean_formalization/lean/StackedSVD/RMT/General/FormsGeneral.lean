/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.IsoMixed
import StackedSVD.RMT.General.FormsBridge
import StackedSVD.RMT.General.Forms
import StackedSVD.RMT.General.FormsSup
import StackedSVD.RMT.General.FormsLimits

/-!
# Item (a) of Stage 1: the six complex resolvent forms at a general law

The endpoint of item (a) of `notes/archive/prop_single_table_general.md`, the frozen
statement `SpikedModel.resolventFormsC_of_general` of the plan (section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35): at i.i.d. noise
of a fixed law `ν` with mean zero, variance one and a finite fourth moment (`NoiseLaw ν`), the
six
isotropic forms of `SpikedModel.ResolventFormsC` converge in probability at every `z` in the
open upper half plane.

Two things happen here.

1. **L3**, `SpikedModel.cformC_gram_v_gvec_tendsto`: the mixed form `vᵀ G g` goes to zero in
   probability. It is the one scalar limit of the six that no unit upstream of this file
   produces, because `g = Eᵀ u` is a function of the same noise matrix as `G`. The measure
   bound at finite `n` and `d` is statement 5 of `RMT/General/IsoMixed.lean`
   (`GenRMT.measure_cformC_mixed_ge_le`); this file transports it through the law
   `hG N : HasLaw (m.Z N) (noiseMatrix ν (n N) (d N)) (μ N)` and lets the rate vanish.
2. **The chaining**, `SpikedModel.resolventFormsC_of_general'`: the six scalar limits L1, L2,
   L3, L1', L2', L3' feed `SpikedModel.resolventFormsC_of_limits`.

## Where the six limits come from

| limit | statement | file |
|---|---|---|
| L1  | `qformC_gram_v_tendsto`        | `RMT/General/FormsLimits.lean` (unit F2) |
| L2  | `qformC_gramC_u_tendsto`       | `RMT/General/FormsLimits.lean` (unit F2) |
| L3  | `cformC_gram_v_gvec_tendsto`   | this file |
| L1' | `qform2C_gram_v_tendsto`       | `RMT/General/Forms.lean` (unit F1a) |
| L2' | `qform2C_gramC_u_tendsto`      | `RMT/General/Forms.lean` (unit F1a) |
| L3' | `cform2C_gram_v_gvec_tendsto`  | `RMT/General/Forms.lean` (unit F1a) |

and the assembly `resolventFormsC_of_limits` is `RMT/General/FormsSup.lean` (unit F1b).
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal

namespace StackedSVD

/-! ### L3: the mixed form `vᵀ G g` -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! F37 (2026-09-09): `measurable_mixDir` used to be repeated here (`private`, to avoid a name
clash with the public `GenRMT.measurable_mix_dir` of `RMT/General/IsoMixed.lean:858`, which
this file already imports). That clash does not arise (the two names differ, and Lean's
`GenRMT` and `SpikedModel` namespaces do not collide), so this file now uses
`GenRMT.measurable_mix_dir` directly. -/

/-- **L3 of item (a).** The mixed form `vᵀ G g` goes to zero in probability, at i.i.d. noise of
a general law `ν` with four moments. `G` is the resolvent of `gram Z` at a complex `z` with
`Im z > 0` and `g = m.gvec = Eᵀ u = (√d)⁻¹ (Zᵀ u)`, so the direction on the right is a function
of the same noise matrix as the resolvent: this is the one limit of the six that the isotropic
law for a fixed direction (`GenRMT.measure_qformC_sub_ge_le`) cannot give. The finite `(n, d)`
bound is statement 5 of unit G4 (`GenRMT.measure_cformC_mixed_ge_le`,
`RMT/General/IsoMixed.lean`); this proof transports it from `noiseMatrix ν (n N) (d N)` to
`μ N` through the law of `m.Z N` and lets the rate vanish. -/
theorem cformC_gram_v_gvec_tendsto [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ}
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0 := by
  have hK : Tendsto (fun N => GenRMT.isoRateMixed ν z (n N) (d N)) atTop (𝓝 0) :=
    GenRMT.tendsto_isoRateMixed ν hz hreg.2.1 hreg.2.2
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N))
          (m.gvec N ω)‖ - 0|}
        ≤ ENNReal.ofReal (GenRMT.isoRateMixed ν z (n N) (d N) / ε ^ 2) := by
    intro ε hε N
    have hxx : WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 1 :=
      le_of_eq (m.dotProduct_v_self N)
    have hyy : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) ≤ 1 :=
      le_of_eq (m.dotProduct_u_self N)
    have hmc : Measurable fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        R4C.cformC (GenRMT.gram Y) z (WithLp.ofLp (m.v N))
          ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N))) :=
      GenRMT.measurable_cformC (W := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => GenRMT.gram Y)
        GenRMT.measurable_gram_self z
        (x := fun _ : Matrix (Fin (n N)) (Fin (d N)) ℝ => WithLp.ofLp (m.v N))
        (y := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
          (Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))
        (fun _ => measurable_const) (fun i => GenRMT.measurable_mix_dir (WithLp.ofLp (m.u N)) i)
    have hmf : Measurable fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        ‖R4C.cformC (GenRMT.gram Y) z (WithLp.ofLp (m.v N))
          ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖ := hmc.norm
    have hmset : MeasurableSet {Y : Matrix (Fin (n N)) (Fin (d N)) ℝ |
        ε ≤ ‖R4C.cformC (GenRMT.gram Y) z (WithLp.ofLp (m.v N))
          ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖} :=
      measurableSet_le measurable_const hmf
    have hlaw := (hG N).measure_eq (p := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ε ≤ ‖R4C.cformC (GenRMT.gram Y) z (WithLp.ofLp (m.v N))
        ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖) hmset
    have hpt : ∀ ω : Ω N,
        |‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖ - 0|
          = ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N))
              ((Real.sqrt (d N))⁻¹ • ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)))‖ := by
      intro ω
      rw [sub_zero, abs_of_nonneg (norm_nonneg _), m.gvec_eq_smul N ω]
    have hseteq : {ω : Ω N |
        ε ≤ |‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖ - 0|}
        = {ω | ε ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N))
            ((Real.sqrt (d N))⁻¹ • ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)))‖} := by
      ext ω
      simp only [Set.mem_ofPred_eq, hpt ω]
    rw [hseteq, hlaw]
    exact GenRMT.measure_cformC_mixed_ge_le hν hz (m.hn N) (m.hd N) _ hxx _ hyy hε
  have h0 := Cheb.tendstoInProb_of_meas_le (ctr := fun _ => (0 : ℝ)) hK hbnd
  simpa using h0

/-! ### Item (a): the six forms -/

/-- **Item (a) of `notes/archive/prop_single_table_general.md`**, the frozen statement
`SpikedModel.resolventFormsC_of_general` of the plan (section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35). At i.i.d. noise
of a fixed law `ν` with mean zero, variance one and a finite fourth moment, the six isotropic
resolvent forms of `SpikedModel.ResolventFormsC` converge in probability at every `z` in the
open upper half plane. There is no edge hypothesis and no `hn2 : ∀ N, 2 ≤ n N`.

The six scalar limits come from `SpikedModel.qformC_gram_v_tendsto` and
`qformC_gramC_u_tendsto` (unit F2, `RMT/General/FormsLimits.lean`),
`cformC_gram_v_gvec_tendsto` above, and the three derivative transfers
`qform2C_gram_v_tendsto`, `qform2C_gramC_u_tendsto`, `cform2C_gram_v_gvec_tendsto` (unit F1a,
`RMT/General/Forms.lean`); `resolventFormsC_of_limits` (unit F1b,
`RMT/General/FormsSup.lean`) turns them into the six fields by algebra. -/
theorem resolventFormsC_of_general' [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) :
    m.ResolventFormsC c := by
  have hL1 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qformC
      (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0 :=
    fun z hz => qformC_gram_v_tendsto hc m hreg hν hG z hz
  have hL2 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qformC
      (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) - MP.mTildeC c z‖) 0 :=
    fun z hz => qformC_gramC_u_tendsto hc m hreg hν hG z hz
  have hL3 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.cformC
      (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0 :=
    fun z hz => cformC_gram_v_gvec_tendsto m hreg hν hG z hz
  exact resolventFormsC_of_limits hc m hL1 hL2 hL3
    (fun z hz => qform2C_gram_v_tendsto hc m hreg hν hG hL1 z hz)
    (fun z hz => qform2C_gramC_u_tendsto hc m hreg hν hG hL2 z hz)
    (fun z hz => cform2C_gram_v_gvec_tendsto hc m hreg hν hG hL3 z hz)

end SpikedModel

end StackedSVD
