/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Sup
import StackedSVD.RMT.R3minus

/-!
# Item TailShift: an index shift removes `hn2 : ∀ N, 2 ≤ n N`

Decision D11 of `notes/FLAGGED.md`: item R0 splits one row off `n N`, so R1, R2, R3 and R3⁻
need `2 ≤ n N` for every `N`, while the target `singleTableLaw_of_gaussian` (`RMT.lean`) only
assumes `m.Regime c`, which gives `2 ≤ n N` for large `N` alone. This file closes the gap by
a tail shift. No statement of another file changes.

## Content

1. `SpikedModel.shift m k` reindexes a model by `N ↦ N + k`, together with `shift_Regime`,
   `shift_GaussianNoise`, the projection lemmas (`shift_X` and friends, all `rfl`) and the
   transport of `IsProbabilityMeasure`.
2. Tail invariance. `tendsto_of_shift` and `tendstoInProb_of_shift` turn a limit of the
   shifted model into the limit of the model, and `tendsto_shift`, `tendstoInProb_shift` go
   the other way; `Filter.tendsto_add_atTop_iff_nat` is the one fact behind all four.
   `SingleTableLaw.of_shift` and `ResolventLimits.of_shift` (with `ResolventLimits.shift`)
   apply them field by field. The `topSimple` field is not a limit, so `of_shift` takes it as
   an argument; `singleTableLaw_topSimple_of_gaussian` (item S) proves it for every `N`
   without `hn2`.
3. `exists_shift_two_le`: `Tendsto n atTop atTop` gives a `k` with `2 ≤ n (N + k)` for all `N`.
4. The primed theorems. Every model-level theorem of `RMT/Sup.lean` and `RMT/R3minus.lean`
   that carries `hn2` gets a copy without it, same name plus `'`, same hypotheses otherwise.
   Items R1, R2, R3 state their own results on an abstract block with `hp : ∀ N, 0 < p N`,
   never on a model, so they need no primed copy.

`singleTableLaw_of_gaussian_supercritical'` differs from the target
`singleTableLaw_of_gaussian` by the single hypothesis `hθ : c < m.θ ^ 4`.

STATUS: see `notes/archive/agent_reports/proof_tailshift.md`. The `sorryAx` this file inherited,
until 2026-08-30 13:00, then closed by `RMT/SteinStep.lean`: now standard axioms only; formerly
inherits comes only from the three tracked `sorry` of `RMT/R1.lean`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! ### 1. The shifted model -/

/-- The model reindexed by `N ↦ N + k`: same `θ`, every sequence evaluated at `N + k`. The
probability spaces move with the index, so `Ω` becomes `fun N => Ω (N + k)` and `μ` becomes
`fun N => μ (N + k)`. -/
def shift (m : SpikedModel μ n d) (k : ℕ) :
    SpikedModel (fun N => μ (N + k)) (fun N => n (N + k)) (fun N => d (N + k)) where
  θ := m.θ
  u := fun N => m.u (N + k)
  v := fun N => m.v (N + k)
  Z := fun N => m.Z (N + k)
  hθ := m.hθ
  hn := fun N => m.hn (N + k)
  hd := fun N => m.hd (N + k)
  hu := fun N => m.hu (N + k)
  hv := fun N => m.hv (N + k)
  hZ := fun N => m.hZ (N + k)

section Projections

variable (m : SpikedModel μ n d) (k N : ℕ)

@[simp] theorem shift_θ : (m.shift k).θ = m.θ := rfl

@[simp] theorem shift_u : (m.shift k).u N = m.u (N + k) := rfl

@[simp] theorem shift_v : (m.shift k).v N = m.v (N + k) := rfl

@[simp] theorem shift_Z : (m.shift k).Z N = m.Z (N + k) := rfl

@[simp] theorem shift_E (ω : Ω (N + k)) : (m.shift k).E N ω = m.E (N + k) ω := rfl

@[simp] theorem shift_X (ω : Ω (N + k)) : (m.shift k).X N ω = m.X (N + k) ω := rfl

@[simp] theorem shift_orthUnit : (m.shift k).orthUnit N = m.orthUnit (N + k) := rfl

@[simp] theorem shift_gvec (ω : Ω (N + k)) : (m.shift k).gvec N ω = m.gvec (N + k) ω := rfl

@[simp] theorem shift_W0 (ω : Ω (N + k)) : (m.shift k).W0 N ω = m.W0 (N + k) ω := rfl

end Projections

/-- The regime is a tail property. -/
theorem shift_Regime {c : ℝ} (m : SpikedModel μ n d) (k : ℕ) (hreg : m.Regime c) :
    (m.shift k).Regime c :=
  ⟨(tendsto_add_atTop_iff_nat (f := n) k).mpr hreg.1,
    (tendsto_add_atTop_iff_nat (f := d) k).mpr hreg.2.1,
    (tendsto_add_atTop_iff_nat (f := fun N => (n N : ℝ) / (d N : ℝ)) k).mpr hreg.2.2⟩

/-- The Gaussian noise law is stated index by index, so it shifts by evaluation. -/
theorem shift_GaussianNoise (m : SpikedModel μ n d) (k : ℕ) (hG : m.GaussianNoise) :
    (m.shift k).GaussianNoise := fun N => hG (N + k)

/-- Transport of the `IsProbabilityMeasure` instance to the shifted spaces. -/
theorem isProbabilityMeasure_shift [∀ N, IsProbabilityMeasure (μ N)] (k : ℕ) :
    ∀ N, IsProbabilityMeasure (μ (N + k)) := fun _ => inferInstance

/-! ### 2. Tail invariance of the limits -/

/-- A limit along `atTop` is a tail property (the direction that removes the shift). -/
theorem tendsto_of_shift {α : Type*} {l : Filter α} {F : ℕ → α} (k : ℕ)
    (h : Tendsto (fun N => F (N + k)) atTop l) : Tendsto F atTop l :=
  (tendsto_add_atTop_iff_nat k).mp h

/-- A limit along `atTop` is a tail property (the direction that adds the shift). -/
theorem tendsto_shift {α : Type*} {l : Filter α} {F : ℕ → α} (k : ℕ)
    (h : Tendsto F atTop l) : Tendsto (fun N => F (N + k)) atTop l :=
  (tendsto_add_atTop_iff_nat k).mpr h

/-- Convergence in probability is a tail property (shift removed). -/
theorem tendstoInProb_of_shift {f : ∀ N, Ω N → ℝ} {a : ℝ} (k : ℕ)
    (h : TendstoInProb (fun N => μ (N + k)) (fun N => f (N + k)) a) :
    TendstoInProb μ f a := fun ε hε =>
  (tendsto_add_atTop_iff_nat (f := fun N => μ N {ω | ε ≤ |f N ω - a|}) k).mp (h ε hε)

/-- Convergence in probability is a tail property (shift added). -/
theorem tendstoInProb_shift {f : ∀ N, Ω N → ℝ} {a : ℝ} (k : ℕ) (h : TendstoInProb μ f a) :
    TendstoInProb (fun N => μ (N + k)) (fun N => f (N + k)) a := fun ε hε =>
  (tendsto_add_atTop_iff_nat (f := fun N => μ N {ω | ε ≤ |f N ω - a|}) k).mpr (h ε hε)

/-- **The single-table law is a tail property, except for `topSimple`.** The three limit
fields come from the shifted model; `topSimple` is a statement about one `N` at a time, so it
enters as the second argument (item S proves it for every `N` with no `hn2`). -/
theorem SingleTableLaw.of_shift {m : SpikedModel μ n d} {c : ℝ} (k : ℕ)
    (law : (m.shift k).SingleTableLaw c)
    (hsimple : ∀ N, ∀ᵐ ω ∂(μ N), TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω))) :
    m.SingleTableLaw c := by
  refine ⟨tendstoInProb_of_shift k law.align, fun ε hε => ?_,
    tendstoInProb_of_shift k law.lamMax, hsimple⟩
  exact tendsto_of_shift k (law.delocUniform ε hε)

/-- The seven fields of `ResolventLimits` are limits, so they shift back. -/
theorem ResolventLimits.of_shift {m : SpikedModel μ n d} {c : ℝ} (k : ℕ)
    (H : (m.shift k).ResolventLimits c) : m.ResolventLimits c := by
  refine ⟨fun ε hε => tendsto_of_shift k (H.edge ε hε), fun z hz => ?_, fun z hz => ?_,
    fun z hz => ?_, fun z hz => ?_, fun z hz => ?_, fun z hz => ?_⟩
  · exact tendstoInProb_of_shift k (H.vv z hz)
  · exact tendstoInProb_of_shift k (H.gg z hz)
  · exact tendstoInProb_of_shift k (H.vg z hz)
  · exact tendstoInProb_of_shift k (H.vv2 z hz)
  · exact tendstoInProb_of_shift k (H.gg2 z hz)
  · exact tendstoInProb_of_shift k (H.vg2 z hz)

/-- The same seven fields, in the direction that adds the shift. Every `f` is named: in this
direction the goal is the shifted statement, and `fun N => ?f (N + k)` against it is not a
pattern for the unifier. -/
theorem ResolventLimits.shift {m : SpikedModel μ n d} {c : ℝ} (H : m.ResolventLimits c)
    (k : ℕ) : (m.shift k).ResolventLimits c := by
  refine ⟨fun ε hε => tendsto_shift k (H.edge ε hε), fun z hz => ?_, fun z hz => ?_,
    fun z hz => ?_, fun z hz => ?_, fun z hz => ?_, fun z hz => ?_⟩
  · exact tendstoInProb_shift
      (f := fun N ω => R4.qform (m.W0 N ω) z (WithLp.ofLp (m.v N))) k (H.vv z hz)
  · exact tendstoInProb_shift
      (f := fun N ω => R4.qform (m.W0 N ω) z (m.gvec N ω)) k (H.gg z hz)
  · exact tendstoInProb_shift
      (f := fun N ω => R4.cform (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)) k
      (H.vg z hz)
  · exact tendstoInProb_shift
      (f := fun N ω => R4.qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N))) k (H.vv2 z hz)
  · exact tendstoInProb_shift
      (f := fun N ω => R4.qform2 (m.W0 N ω) z (m.gvec N ω)) k (H.gg2 z hz)
  · exact tendstoInProb_shift
      (f := fun N ω => R4.cform2 (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)) k
      (H.vg2 z hz)

/-! ### 3. The shift that makes the block nonempty -/

/-- `n N → ∞` gives one `k` with `2 ≤ n (N + k)` for every `N`. This is the whole content of
the tail shift: the hypothesis `hn2` of items R1, R2, R3 and R3⁻ holds for the shifted model
under `m.Regime c` alone. -/
theorem exists_shift_two_le {c : ℝ} (m : SpikedModel μ n d) (hreg : m.Regime c) :
    ∃ k, ∀ N, 2 ≤ n (N + k) := by
  obtain ⟨k, hk⟩ := eventually_atTop.mp (hreg.1.eventually_ge_atTop 2)
  exact ⟨k, fun N => hk (N + k) (Nat.le_add_left k N)⟩

/-! ### 4. The primed theorems: the same statements without `hn2` -/

section Primed

variable [∀ N, IsProbabilityMeasure (μ N)]

/-- `Sup.resolventLimits_of_gaussian` without `hn2`. The binder order copies the unprimed
theorem: only `hn2` is gone. -/
theorem resolventLimits_of_gaussian' {c : ℝ} (hc : 0 < c) (m : SpikedModel μ n d)
    (hreg : m.Regime c) (hG : m.GaussianNoise) : m.ResolventLimits c := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact ResolventLimits.of_shift k
    (resolventLimits_of_gaussian hc (m.shift k) (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk)

/-- **`Sup.singleTableLaw_of_gaussian_supercritical` without `hn2`.** Against the target
`singleTableLaw_of_gaussian` (`RMT.lean`) this keeps one extra hypothesis, `hθ : c < m.θ ^ 4`.
-/
theorem singleTableLaw_of_gaussian_supercritical' {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) (hθ : c < m.θ ^ 4) :
    m.SingleTableLaw c := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact SingleTableLaw.of_shift k
    (singleTableLaw_of_gaussian_supercritical hc (m.shift k) (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk hθ)
    (singleTableLaw_topSimple_of_gaussian m hG)

section Model

variable {m : SpikedModel μ n d} {c : ℝ}

omit [∀ N, IsProbabilityMeasure (μ N)] in
/-- `R3minus.tendstoInProb_stieltjes2C_W0` without `hn2`. -/
theorem tendstoInProb_stieltjes2C_W0' (hc : 0 < c) (hreg : m.Regime c) (hG : m.GaussianNoise)
    {z : ℂ} (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (m.W0 N ω) z - MP.mCDeriv c z‖) 0 := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  exact tendstoInProb_of_shift k
    (tendstoInProb_stieltjes2C_W0 (m := m.shift k) hc (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk hz)

/-- `R3minus.tendstoInProb_trace2_W0` without `hn2`. -/
theorem tendstoInProb_trace2_W0' (H : m.ResolventLimits c) (hc : 0 < c) (hreg : m.Regime c)
    (hG : m.GaussianNoise) {x : ℝ} (hx : bulkEdge c < x) :
    TendstoInProb μ (fun N ω =>
      (d N : ℝ)⁻¹ * (R4.resolv (m.W0 N ω) x * R4.resolv (m.W0 N ω) x).trace)
      (MP.mDeriv c x) := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact tendstoInProb_of_shift k
    (tendstoInProb_trace2_W0 (H.shift k) hc (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk hx)

/-- `R3minus.tendsto_measure_lamMax_W0_le_sub` without `hn2`. -/
theorem tendsto_measure_lamMax_W0_le_sub' (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c - δ})
      atTop (𝓝 0) := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact tendsto_of_shift k
    (tendsto_measure_lamMax_W0_le_sub (H.shift k) hc (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk hδ)

/-- `R3minus.tendsto_measure_gramLamMax_le_sub` without `hn2`. -/
theorem tendsto_measure_gramLamMax_le_sub' (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | gramLamMax (m.X N ω) ≤ bulkEdge c - δ}) atTop (𝓝 0) := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact tendsto_of_shift k
    (tendsto_measure_gramLamMax_le_sub (H.shift k) hc (m.shift_Regime k hreg)
      (m.shift_GaussianNoise k hG) hk hδ)

/-- **`R3minus.lamMax_tendstoInProb_of_subcritical` without `hn2`.** This is the subcritical
`lamMax` field of `SingleTableLaw` under `m.Regime c` alone. -/
theorem lamMax_tendstoInProb_of_subcritical' (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) (hθ : m.θ ^ 4 ≤ c) :
    TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) (rhoSq m.θ c) := by
  obtain ⟨k, hk⟩ := m.exists_shift_two_le hreg
  have : ∀ N : ℕ, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  have h := lamMax_tendstoInProb_of_subcritical (H.shift k) hc (m.shift_Regime k hreg)
    (m.shift_GaussianNoise k hG) hk hθ
  exact tendstoInProb_of_shift (f := fun N ω => gramLamMax (m.X N ω)) k h

end Model

end Primed

end SpikedModel

end StackedSVD
