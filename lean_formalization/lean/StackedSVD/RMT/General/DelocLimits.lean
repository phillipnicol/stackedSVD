/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.General.Iso
import StackedSVD.RMT.General.ProbC
import StackedSVD.RMT.R5
import StackedSVD.RMT.R4C

/-!
# Probabilistic limits for item R6' at a general noise law (Stage 1, unit D1b)

The probabilistic inputs of `align_tendstoInProb_of_subcritical_general`, the subcritical
`align` field of `prop:single_table` at a general noise law. This was the frozen statement of
the plan (`notes/archive/prop_single_table_general.md`, section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35). The assembly
that consumes them is `RMT/General/DelocAlign.lean` (unit D1c).

## Content

1. `tendstoInProb_qform2_qvec`: the squared real resolvent form at `q = θ v + g` converges to
   `(θ² + 1) m'(z)` above the bulk edge. A model-free copy of
   `SpikedModel.tendstoInProb_qform2_qvec` (`RMT/R6.lean:856`), which lives in the Gaussian
   file `RMT/R6.lean` only because that file is where item R6 was proved.
2. `aOf_tendsto`: the scalar `a = (v ⬝ q) / (q ⬝ q)` converges to `θ / (θ² + 1)`, from the two
   scalar limits of `RMT/General/Defs.lean`.
3. `qformC_W0_rdir_tendsto` and `cformC_W0_qvec_rdir_tendsto`: the complex resolvent forms in
   the direction `r = v - a q`, from the six forms of `SpikedModel.ResolventFormsC`.

Every theorem takes the model `m` explicitly. The complex limits are stated in the norm form
`TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0` of `ResolventFormsC`. The closure
lemmas for that form come from `RMT/General/ProbC.lean` (unit F1a). The bilinearity of
`R4C.cformC` stays a private copy here, with the suffix `_dl`: its right-hand twins live in
the private block around `RMT/General/Iso.lean:2258` and have no public home to import.
F37 (2026-09-09): `cformC_comm_dl` is dropped, since `GenRMT.cformC_comm`
(`RMT/General/Iso.lean:513`) is already public; this file now imports `Iso`, which the
earlier note here wrongly said it could not (checked: no import cycle results, `Iso.lean`
does not transitively depend on this file).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace GenRMT

namespace Deloc

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}
  [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ}

set_option linter.unusedSectionVars false

/-! ### Part 1: the squared real form at `q` -/

/-- **(H2) at `q = θ v + g`.** The squared resolvent form at the rank-one direction converges
to `(θ² + 1) m'(z)` at every real `z` above the bulk edge. Model-free copy of
`SpikedModel.tendstoInProb_qform2_qvec` (`RMT/R6.lean:856-877`) with the model explicit: the
proof uses only the `vv2`, `vg2` and `gg2` fields of `ResolventLimits` and the bilinearity of
`R4.qform2`. It serves item R6' of `prop:single_table`. -/
theorem tendstoInProb_qform2_qvec (m : SpikedModel μ n d) (RL : m.ResolventLimits c) {z : ℝ}
    (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.qform2 (m.W0 N ω) z (m.qvec N ω))
      ((m.θ ^ 2 + 1) * MP.mDeriv c z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), R4.qform2 (m.W0 N ω) z (m.qvec N ω)
      = m.θ ^ 2 * R4.qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N))
        + 2 * m.θ * R4.cform2 (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
        + R4.qform2 (m.W0 N ω) z (m.gvec N ω) := by
    intro N ω
    change R4.qform2 (m.W0 N ω) z (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
    rw [R4.qform2_smul_add (m.isHermitian_W0 N ω)]
  have h := (((RL.vv2 z hz).const_mul (m.θ ^ 2)).add
    ((RL.vg2 z hz).const_mul (2 * m.θ))).add (RL.gg2 z hz)
  have hlim : m.θ ^ 2 * MP.mDeriv c z + 2 * m.θ * 0 + MP.mDeriv c z
      = (m.θ ^ 2 + 1) * MP.mDeriv c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-! ### Part 2: the scalar `a = (v ⬝ q) / (q ⬝ q)` -/

/-- **The projection coefficient.** `a = (v ⬝ q) / (q ⬝ q) → θ / (θ² + 1)` for i.i.d. noise of
a fixed law `ν` with mean 0, variance 1 and a finite fourth moment. The numerator is
`θ (v ⬝ v) + v ⬝ g = θ + v ⬝ g → θ` and the denominator is `q ⬝ q → θ² + 1 > 0`.

The quotient is `R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω)` of the Gaussian file `RMT/R6.lean`
and `Deloc.aOf` of unit D1a (`RMT/General/Deloc.lean`). This file may not import either, so
the statement writes the quotient out. -/
theorem aOf_tendsto (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) :
    TendstoInProb μ
      (fun N ω => (WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω))
      (m.θ / (m.θ ^ 2 + 1)) := by
  have hnum : TendstoInProb μ (fun N ω => WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) m.θ := by
    have hfun : ∀ (N : ℕ) (ω : Ω N), WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω
        = m.θ + WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω := by
      intro N ω
      change WithLp.ofLp (m.v N) ⬝ᵥ (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
      rw [dotProduct_add, dotProduct_smul, smul_eq_mul, m.dotProduct_v_self N, mul_one]
    have h := (TendstoInProb.const μ m.θ).add
      (m.tendstoInProb_dotProduct_v_gvec_general hν hG hreg.2.1)
    rw [add_zero] at h
    simpa only [hfun] using h
  exact hnum.div (m.tendstoInProb_dotProduct_qvec_general hν hG hreg.2.1) (by positivity)

/-! ### The bilinearity of the complex resolvent forms

Private copies of the right-hand twins of the private block at `RMT/General/Iso.lean:2258`,
which has no public home to import. F37 (2026-09-09): `cvec_smul_dl` and `cformC_comm_dl` are
dropped; `GenRMT.cvec_smul` (`Companion.lean:404`, public) and `GenRMT.cformC_comm`
(`Iso.lean:513`, public) are used directly instead.
-/

private theorem cvec_add_dl {D : ℕ} (x y : Fin D → ℝ) :
    R4C.cvec (x + y) = R4C.cvec x + R4C.cvec y := by
  funext i; simp [R4C.cvec]

/-- The form is additive in its right argument. -/
private theorem cformC_add_right_dl {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ)
    (w x y : Fin D → ℝ) :
    R4C.cformC W z w (x + y) = R4C.cformC W z w x + R4C.cformC W z w y := by
  simp only [R4C.cformC, cvec_add_dl, Matrix.mulVec_add, dotProduct_add]

/-- The form is homogeneous in its right argument. -/
private theorem cformC_smul_right_dl {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (r : ℝ)
    (w x : Fin D → ℝ) :
    R4C.cformC W z w (r • x) = (r : ℂ) * R4C.cformC W z w x := by
  simp only [R4C.cformC, SpikedModel.cvec_smul, Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul]

/-- The form on a two-term combination in its right argument. -/
private theorem cformC_lin2_right_dl {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (s t : ℝ)
    (w x y : Fin D → ℝ) :
    R4C.cformC W z w (s • x + t • y)
      = (s : ℂ) * R4C.cformC W z w x + (t : ℂ) * R4C.cformC W z w y := by
  rw [cformC_add_right_dl, cformC_smul_right_dl, cformC_smul_right_dl]

/-- The quadratic form on a two-term combination. -/
private theorem qformC_lin2_dl {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} (hW : W.IsHermitian)
    {z : ℂ} (hz : z.im ≠ 0) (s t : ℝ) (x y : Fin D → ℝ) :
    R4C.qformC W z (s • x + t • y)
      = (s : ℂ) ^ 2 * R4C.qformC W z x + 2 * ((s : ℂ) * (t : ℂ)) * R4C.cformC W z x y
        + (t : ℂ) ^ 2 * R4C.qformC W z y := by
  simp only [R4C.qformC]
  rw [cformC_lin2_right_dl, cformC_comm hW hz (s • x + t • y) x,
    cformC_comm hW hz (s • x + t • y) y, cformC_lin2_right_dl, cformC_lin2_right_dl,
    cformC_comm hW hz y x]
  ring

/-! ### Part 3: the complex forms in the direction `r = v - a q` -/

/-- `v ᵀ G₀(z) q → θ m(z)`, from the `vvC` and `vgC` fields and `q = θ v + g`. -/
private theorem cformC_v_qvec_tendsto (m : SpikedModel μ n d) (H : m.ResolventFormsC c)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.qvec N ω)
      - (m.θ : ℂ) * MP.mC c z‖) 0 := by
  have hfun : ∀ (N : ℕ) (ω : Ω N),
      R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.qvec N ω)
        = (m.θ : ℂ) * R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N))
          + R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω) := by
    intro N ω
    have hq : m.qvec N ω = m.θ • WithLp.ofLp (m.v N) + (1 : ℝ) • m.gvec N ω := by
      change m.θ • WithLp.ofLp (m.v N) + m.gvec N ω = _
      rw [one_smul]
    rw [hq, cformC_lin2_right_dl]
    simp only [R4C.qformC]
    push_cast
    ring
  have hvg : TendstoInProb μ
      (fun N ω => ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω) - 0‖) 0 := by
    simpa using H.vgC z hz
  have h := tendstoInProbC_add (tendstoInProbC_const_mul (m.θ : ℂ) (H.vvC z hz)) hvg
  have hlim : (m.θ : ℂ) * MP.mC c z + 0 = (m.θ : ℂ) * MP.mC c z := add_zero _
  rw [hlim] at h
  simpa only [hfun] using h

/-- `q ᵀ G₀(z) q → (θ² + 1) m(z)`, from the `vvC`, `ggC` and `vgC` fields. -/
private theorem qformC_qvec_tendsto (m : SpikedModel μ n d) (H : m.ResolventFormsC c)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0 N ω) z (m.qvec N ω)
      - ((m.θ : ℂ) ^ 2 + 1) * MP.mC c z‖) 0 := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), R4C.qformC (m.W0 N ω) z (m.qvec N ω)
      = (m.θ : ℂ) ^ 2 * R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N))
        + 2 * (m.θ : ℂ) * R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
        + R4C.qformC (m.W0 N ω) z (m.gvec N ω) := by
    intro N ω
    have hq : m.qvec N ω = m.θ • WithLp.ofLp (m.v N) + (1 : ℝ) • m.gvec N ω := by
      change m.θ • WithLp.ofLp (m.v N) + m.gvec N ω = _
      rw [one_smul]
    rw [hq, qformC_lin2_dl (m.isHermitian_W0 N ω) hz.ne']
    push_cast
    ring
  have hvg : TendstoInProb μ
      (fun N ω => ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω) - 0‖) 0 := by
    simpa using H.vgC z hz
  have h := tendstoInProbC_add
    (tendstoInProbC_add (tendstoInProbC_const_mul ((m.θ : ℂ) ^ 2) (H.vvC z hz))
      (tendstoInProbC_const_mul (2 * (m.θ : ℂ)) hvg)) (H.ggC z hz)
  have hlim : (m.θ : ℂ) ^ 2 * MP.mC c z + 2 * (m.θ : ℂ) * 0 + MP.mC c z
      = ((m.θ : ℂ) ^ 2 + 1) * MP.mC c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-- `(θ² + 1) ≠ 0` over `ℂ`. -/
private theorem thetaSq_add_one_ne_zero_dl (m : SpikedModel μ n d) :
    ((m.θ : ℂ) ^ 2 + 1) ≠ 0 := by
  have h1 : (0 : ℝ) < m.θ ^ 2 + 1 := by positivity
  have h2 : ((m.θ ^ 2 + 1 : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr h1.ne'
  push_cast at h2
  exact h2

/-- **The quadratic form in the direction `r = v - a q`.** With `a = (v ⬝ q) / (q ⬝ q)` the
complex resolvent form at `r` converges to `m(z) / (θ² + 1)` at every `z` with `Im z > 0`.

Write `r = (1 - a θ) v - a g`. The expansion in `v` and `q`,
`Φ_r = Φ_v - 2 a Ψ_{v,q} + a² Φ_q`, together with `a → θ / (θ² + 1)`, `Ψ_{v,q} → θ m` and
`Φ_q → (θ² + 1) m`, gives `m (1 - θ² / (θ² + 1)) = m / (θ² + 1)`. It serves item R6' of
`prop:single_table` at a general noise law. -/
theorem qformC_W0_rdir_tendsto (m : SpikedModel μ n d) (H : m.ResolventFormsC c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N)
        - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)) • m.qvec N ω)
      - MP.mC c z / (m.θ ^ 2 + 1)‖) 0 := by
  have hA := tendstoInProbC_ofReal (aOf_tendsto m hreg hν hG)
  have hfun : ∀ (N : ℕ) (ω : Ω N), R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N)
        - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)) • m.qvec N ω)
      = R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N))
        + (-2 : ℂ) * (((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω)
              / (m.qvec N ω ⬝ᵥ m.qvec N ω) : ℝ) : ℂ)
            * R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.qvec N ω)
        + ((((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω) : ℝ) : ℂ)
            * ((((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω)
                / (m.qvec N ω ⬝ᵥ m.qvec N ω) : ℝ) : ℂ)))
            * R4C.qformC (m.W0 N ω) z (m.qvec N ω) := by
    intro N ω
    have hr : WithLp.ofLp (m.v N)
          - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)) • m.qvec N ω
        = (1 : ℝ) • WithLp.ofLp (m.v N)
          + (-((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)))
              • m.qvec N ω := by
      rw [one_smul, neg_smul, ← sub_eq_add_neg]
    rw [hr, qformC_lin2_dl (m.isHermitian_W0 N ω) hz.ne']
    push_cast
    ring
  have h := tendstoInProbC_add
    (tendstoInProbC_add (H.vvC z hz)
      (tendstoInProbC_mul (tendstoInProbC_const_mul (-2 : ℂ) hA)
        (cformC_v_qvec_tendsto m H z hz)))
    (tendstoInProbC_mul (tendstoInProbC_mul hA hA) (qformC_qvec_tendsto m H z hz))
  have hne := thetaSq_add_one_ne_zero_dl m
  have hcast : ((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ) = (m.θ : ℂ) / ((m.θ : ℂ) ^ 2 + 1) := by
    push_cast
    ring
  have hlim : MP.mC c z
        + (-2 : ℂ) * ((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ) * ((m.θ : ℂ) * MP.mC c z)
        + ((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ) * ((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ)
            * (((m.θ : ℂ) ^ 2 + 1) * MP.mC c z)
      = MP.mC c z / ((m.θ : ℂ) ^ 2 + 1) := by
    rw [hcast]
    field_simp
    ring
  rw [hlim] at h
  simpa only [hfun] using h

/-- **The cross form of `q` with `r = v - a q`.** It vanishes: `r` is the component of `v`
orthogonal to `q` in the metric of the resolvent, because
`Ψ_{q,r} = Ψ_{v,q} - a Φ_q → θ m - θ m = 0`. It serves item R6' of `prop:single_table` at a
general noise law. -/
theorem cformC_W0_qvec_rdir_tendsto (m : SpikedModel μ n d) (H : m.ResolventFormsC c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0 N ω) z (m.qvec N ω) (WithLp.ofLp (m.v N)
      - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω))
        • m.qvec N ω)‖) 0 := by
  have hA := tendstoInProbC_ofReal (aOf_tendsto m hreg hν hG)
  have hfun : ∀ (N : ℕ) (ω : Ω N), R4C.cformC (m.W0 N ω) z (m.qvec N ω) (WithLp.ofLp (m.v N)
        - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)) • m.qvec N ω)
      = R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.qvec N ω)
        + -(((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω) : ℝ) : ℂ)
            * R4C.qformC (m.W0 N ω) z (m.qvec N ω) := by
    intro N ω
    have hr : WithLp.ofLp (m.v N)
          - ((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)) • m.qvec N ω
        = (1 : ℝ) • WithLp.ofLp (m.v N)
          + (-((WithLp.ofLp (m.v N) ⬝ᵥ m.qvec N ω) / (m.qvec N ω ⬝ᵥ m.qvec N ω)))
              • m.qvec N ω := by
      rw [one_smul, neg_smul, ← sub_eq_add_neg]
    rw [hr, cformC_lin2_right_dl,
      cformC_comm (m.isHermitian_W0 N ω) hz.ne' (m.qvec N ω) (WithLp.ofLp (m.v N))]
    simp only [R4C.qformC]
    push_cast
    ring
  have h := tendstoInProbC_add (cformC_v_qvec_tendsto m H z hz)
    (tendstoInProbC_mul (tendstoInProbC_neg hA) (qformC_qvec_tendsto m H z hz))
  have hne := thetaSq_add_one_ne_zero_dl m
  have hcast : ((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ) = (m.θ : ℂ) / ((m.θ : ℂ) ^ 2 + 1) := by
    push_cast
    ring
  have hlim : (m.θ : ℂ) * MP.mC c z
      + -((m.θ / (m.θ ^ 2 + 1) : ℝ) : ℂ) * (((m.θ : ℂ) ^ 2 + 1) * MP.mC c z) = 0 := by
    rw [hcast]
    field_simp
    ring
  rw [hlim] at h
  simpa only [hfun, sub_zero] using h

end Deloc

end GenRMT

end StackedSVD
