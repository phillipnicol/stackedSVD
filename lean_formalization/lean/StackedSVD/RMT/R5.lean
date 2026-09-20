/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RMT.R4C
import StackedSVD.RMT.MP
import StackedSVD.RMT.T
import StackedSVD.Spectral

/-!
# Item R5: the analytic core

See `notes/archive/rmt_R5.md` (restructured 2026-08-30, decision D6). Nothing here is Gaussian and
nothing here uses independence. R5 turns two probabilistic inputs, an edge bound (H1) and six
resolvent limits at every fixed real `z > bulkEdge c` (H2), into three conclusions:

* `lamMax_tendstoInProb`  : `λ_max(Xᵀ X) → ρ² = rhoSq θ c` in probability (`c < θ⁴`);
* `align_tendstoInProb`   : `overlap(X, v) → β² = betaSq θ c` in probability (`c < θ⁴`);

Cleanup wave 3 (2026-08-31) deleted two consumer-less lemmas of this file,
`tendsto_measure_localized` (step 4) and `tendsto_measure_topSimple` (the w.h.p. simplicity;
the almost sure form that the project uses is item S, `RMT/Simplicity.lean`). Both follow in
five lines from `exists_localization_bad` and `localize` if a later task needs them.
* `tendsto_measure_lamMax_le_of_subcritical` : the upper edge bound for `θ⁴ ≤ c`.

The two inputs are bundled in the structure `SpikedModel.ResolventLimits` (7 fields, D6).
The route, with `q = θ v + g`, `W₀ = E⊥ᵀ E⊥` and `Xᵀ X = W₀ + q qᵀ` (item R0):

1. `R4.qform_smul_add` expands the secular function of `q` in the six (H2) forms, so
   `secular W₀ q z → 1 + (θ² + 1) m(z)` at every fixed `z > b` (`tendstoInProb_secular`).
2. Above threshold `1 + (θ²+1) m` is negative at `zlo = ρ - δ` and positive at `zhi = ρ + δ`,
   so with probability tending to `1` the intermediate value theorem puts the unique secular
   root, which is `λ_max(Xᵀ X)` by `R4.lamMax_eq`, inside `(zlo, zhi)` (`localize`).
3. `R4.qform_le_qform` (increasing) and `R4.qform2_antitoneOn` (decreasing) bracket every form
   at the random root between its values at `zlo` and `zhi`; the cross forms follow by
   polarization. This gives `overlap_in_bracket`, with the four random endpoints
   `numLoR`, `numHiR`, `denLoR`, `denHiR`.
4. Those endpoints converge to the deterministic `R5.numLoD`, ..., `R5.denHiD`, whose brackets
   shrink to `betaSq θ c` as `δ ↓ 0` (`R5.exists_delta`, from `MP.continuousOn_m`,
   `MP.continuousOn_mDeriv`, `MP.overlap_identity`).

STATUS: see `notes/archive/agent_reports/proof_r5.md`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### Transport along an equality of matrices

`lamMax_congr`, `topProj_congr` and `topSimple_congr` (`Spectral.lean`) move `lamMax`,
`topProj` and `TopSimple` across the identity `Xᵀ X = W₀ + q qᵀ` of item R0. This file held
private copies of them until cleanup wave 2. -/

namespace R4

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} {z θ : ℝ}

/-! ### Bilinearity of the four deterministic forms

`R4C.lean` proves the two polarization identities at every `z`. These are the remaining
scaling and additivity facts that R5 needs to expand `q = θ v + g`. -/

theorem qform_add (hW : W.IsHermitian) (z : ℝ) (x y : Fin d → ℝ) :
    qform W z (x + y) = qform W z x + 2 * cform W z x y + qform W z y := by
  rw [cform_eq_polarization hW z x y]; ring

theorem qform2_add (hW : W.IsHermitian) (z : ℝ) (x y : Fin d → ℝ) :
    qform2 W z (x + y) = qform2 W z x + 2 * cform2 W z x y + qform2 W z y := by
  rw [cform2_eq_polarization hW z x y]; ring

theorem qform_smul (z θ : ℝ) (x : Fin d → ℝ) : qform W z (θ • x) = θ ^ 2 * qform W z x := by
  simp only [qform, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul]
  ring

theorem qform2_smul (z θ : ℝ) (x : Fin d → ℝ) :
    qform2 W z (θ • x) = θ ^ 2 * qform2 W z x := by
  simp only [qform2, Matrix.mulVec_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul]
  ring

theorem cform_smul_left (z θ : ℝ) (x y : Fin d → ℝ) :
    cform W z (θ • x) y = θ * cform W z x y := by
  simp only [cform, smul_dotProduct, smul_eq_mul]

theorem cform2_smul_left (z θ : ℝ) (x y : Fin d → ℝ) :
    cform2 W z (θ • x) y = θ * cform2 W z x y := by
  simp only [cform2, smul_dotProduct, smul_eq_mul]

theorem cform_add_right (z : ℝ) (x y y' : Fin d → ℝ) :
    cform W z x (y + y') = cform W z x y + cform W z x y' := by
  simp only [cform, Matrix.mulVec_add, dotProduct_add]

theorem cform_smul_right (z θ : ℝ) (x y : Fin d → ℝ) :
    cform W z x (θ • y) = θ * cform W z x y := by
  simp only [cform, Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul]

/-- The numerator of the overlap formula, expanded in `v` and `g`. -/
theorem cform_smul_add (z θ : ℝ) (x y : Fin d → ℝ) :
    cform W z x (θ • x + y) = θ * qform W z x + cform W z x y := by
  rw [cform_add_right, cform_smul_right]
  rfl

/-- The secular function of `q = θ v + g`, expanded in `v` and `g`. -/
theorem qform_smul_add (hW : W.IsHermitian) (z θ : ℝ) (x y : Fin d → ℝ) :
    qform W z (θ • x + y)
      = θ ^ 2 * qform W z x + 2 * θ * cform W z x y + qform W z y := by
  rw [qform_add hW, qform_smul, cform_smul_left]; ring

/-- The denominator of the overlap formula, expanded in `v` and `g`. -/
theorem qform2_smul_add (hW : W.IsHermitian) (z θ : ℝ) (x y : Fin d → ℝ) :
    qform2 W z (θ • x + y)
      = θ ^ 2 * qform2 W z x + 2 * θ * cform2 W z x y + qform2 W z y := by
  rw [qform2_add hW, qform2_smul, cform2_smul_left]; ring

/-! ### The deterministic sandwich (step 6) -/

/-- Monotonicity of `Φ` brackets the numerator at the random root between its two endpoint
values. -/
theorem num_bracket (hW : W.IsHermitian) {zlo zhi lam : ℝ} (hθ : 0 ≤ θ)
    (h₀ : lamMax W hW < zlo) (h₁ : zlo ≤ lam) (h₂ : lam ≤ zhi) (V G : Fin d → ℝ) :
    θ * qform W zlo V + (qform W zlo (V + G) - qform W zhi V - qform W zhi G) / 2
        ≤ cform W lam V (θ • V + G) ∧
      cform W lam V (θ • V + G)
        ≤ θ * qform W zhi V + (qform W zhi (V + G) - qform W zlo V - qform W zlo G) / 2 := by
  have hlam : lamMax W hW < lam := lt_of_lt_of_le h₀ h₁
  have e1 : cform W lam V (θ • V + G)
      = θ * qform W lam V + (qform W lam (V + G) - qform W lam V - qform W lam G) / 2 := by
    rw [cform_smul_add, cform_eq_polarization hW]
  have hV₁ : qform W zlo V ≤ qform W lam V := qform_le_qform hW h₀ h₁
  have hV₂ : qform W lam V ≤ qform W zhi V := qform_le_qform hW hlam h₂
  have hG₁ : qform W zlo G ≤ qform W lam G := qform_le_qform hW h₀ h₁
  have hG₂ : qform W lam G ≤ qform W zhi G := qform_le_qform hW hlam h₂
  have hS₁ : qform W zlo (V + G) ≤ qform W lam (V + G) := qform_le_qform hW h₀ h₁
  have hS₂ : qform W lam (V + G) ≤ qform W zhi (V + G) := qform_le_qform hW hlam h₂
  have hm₁ : θ * qform W zlo V ≤ θ * qform W lam V := mul_le_mul_of_nonneg_left hV₁ hθ
  have hm₂ : θ * qform W lam V ≤ θ * qform W zhi V := mul_le_mul_of_nonneg_left hV₂ hθ
  rw [e1]
  constructor <;> linarith

/-- Antitonicity of `Φ²` brackets the denominator at the random root. -/
theorem den_bracket (hW : W.IsHermitian) {zlo zhi lam : ℝ} (hθ : 0 ≤ θ)
    (h₀ : lamMax W hW < zlo) (h₁ : zlo ≤ lam) (h₂ : lam ≤ zhi) (V G : Fin d → ℝ) :
    θ ^ 2 * qform2 W zhi V + θ * (qform2 W zhi (V + G) - qform2 W zlo V - qform2 W zlo G)
          + qform2 W zhi G ≤ qform2 W lam (θ • V + G) ∧
      qform2 W lam (θ • V + G) ≤ θ ^ 2 * qform2 W zlo V
        + θ * (qform2 W zlo (V + G) - qform2 W zhi V - qform2 W zhi G) + qform2 W zlo G := by
  have hlam : lamMax W hW < lam := lt_of_lt_of_le h₀ h₁
  have hmlo : zlo ∈ Set.Ioi (lamMax W hW) := h₀
  have hmL : lam ∈ Set.Ioi (lamMax W hW) := hlam
  have hmhi : zhi ∈ Set.Ioi (lamMax W hW) := lt_of_lt_of_le hlam h₂
  have hV₁ : qform2 W lam V ≤ qform2 W zlo V := qform2_antitoneOn hW V hmlo hmL h₁
  have hV₂ : qform2 W zhi V ≤ qform2 W lam V := qform2_antitoneOn hW V hmL hmhi h₂
  have hG₁ : qform2 W lam G ≤ qform2 W zlo G := qform2_antitoneOn hW G hmlo hmL h₁
  have hG₂ : qform2 W zhi G ≤ qform2 W lam G := qform2_antitoneOn hW G hmL hmhi h₂
  have hS₁ : qform2 W lam (V + G) ≤ qform2 W zlo (V + G) :=
    qform2_antitoneOn hW (V + G) hmlo hmL h₁
  have hS₂ : qform2 W zhi (V + G) ≤ qform2 W lam (V + G) :=
    qform2_antitoneOn hW (V + G) hmL hmhi h₂
  have e1 : qform2 W lam (θ • V + G) = θ ^ 2 * qform2 W lam V
      + θ * (qform2 W lam (V + G) - qform2 W lam V - qform2 W lam G) + qform2 W lam G := by
    rw [qform2_smul_add hW, cform2_eq_polarization hW]; ring
  have hsq : (0 : ℝ) ≤ θ ^ 2 := sq_nonneg θ
  have hm₁ : θ ^ 2 * qform2 W zhi V ≤ θ ^ 2 * qform2 W lam V :=
    mul_le_mul_of_nonneg_left hV₂ hsq
  have hm₂ : θ ^ 2 * qform2 W lam V ≤ θ ^ 2 * qform2 W zlo V :=
    mul_le_mul_of_nonneg_left hV₁ hsq
  have hm₃ : θ * (qform2 W zhi (V + G) - qform2 W zlo V - qform2 W zlo G)
      ≤ θ * (qform2 W lam (V + G) - qform2 W lam V - qform2 W lam G) :=
    mul_le_mul_of_nonneg_left (by linarith) hθ
  have hm₄ : θ * (qform2 W lam (V + G) - qform2 W lam V - qform2 W lam G)
      ≤ θ * (qform2 W zlo (V + G) - qform2 W zhi V - qform2 W zhi G) :=
    mul_le_mul_of_nonneg_left (by linarith) hθ
  rw [e1]
  constructor <;> linarith

end R4

namespace R5

/-! ### The scalar bracket -/

/-- The ratio `x²/y` of a negative bracketed numerator and a positive bracketed denominator
is bracketed. Pure real arithmetic. -/
theorem div_bracket {nl nh x dl dh y : ℝ} (h1 : nl ≤ x) (h2 : x ≤ nh) (hnh : nh < 0)
    (h3 : dl ≤ y) (h4 : y ≤ dh) (hdl : 0 < dl) :
    nh * nh / dh ≤ x * x / y ∧ x * x / y ≤ nl * nl / dl := by
  have hy : 0 < y := lt_of_lt_of_le hdl h3
  have hdh : 0 < dh := lt_of_lt_of_le hy h4
  have hx : x < 0 := lt_of_le_of_lt h2 hnh
  have hsq1 : nh * nh ≤ x * x := by nlinarith
  have hsq2 : x * x ≤ nl * nl := by nlinarith
  have hxx : (0 : ℝ) ≤ x * x := mul_self_nonneg x
  have hnn : (0 : ℝ) ≤ nh * nh := mul_self_nonneg nh
  constructor
  · rw [div_le_div_iff₀ hdh hy]
    have k1 : nh * nh * y ≤ x * x * y := mul_le_mul_of_nonneg_right hsq1 hy.le
    have k2 : x * x * y ≤ x * x * dh := mul_le_mul_of_nonneg_left h4 hxx
    linarith
  · rw [div_le_div_iff₀ hy hdl]
    have k1 : x * x * dl ≤ nl * nl * dl := mul_le_mul_of_nonneg_right hsq2 hdl.le
    have k2 : nl * nl * dl ≤ nl * nl * y := mul_le_mul_of_nonneg_left h3 (mul_self_nonneg nl)
    linarith

/-! ### The deterministic bracket endpoints and the choice of `δ` -/

/-- Lower endpoint of the numerator bracket, in the limit. -/
noncomputable def numLoD (θ c ρ δ : ℝ) : ℝ :=
  θ * MP.m c (ρ - δ) + (MP.m c (ρ - δ) - MP.m c (ρ + δ))

/-- Upper endpoint of the numerator bracket, in the limit. -/
noncomputable def numHiD (θ c ρ δ : ℝ) : ℝ :=
  θ * MP.m c (ρ + δ) + (MP.m c (ρ + δ) - MP.m c (ρ - δ))

/-- Lower endpoint of the denominator bracket, in the limit. -/
noncomputable def denLoD (θ c ρ δ : ℝ) : ℝ :=
  (θ ^ 2 + 1) * MP.mDeriv c (ρ + δ) + 2 * θ * (MP.mDeriv c (ρ + δ) - MP.mDeriv c (ρ - δ))

/-- Upper endpoint of the denominator bracket, in the limit. -/
noncomputable def denHiD (θ c ρ δ : ℝ) : ℝ :=
  (θ ^ 2 + 1) * MP.mDeriv c (ρ - δ) + 2 * θ * (MP.mDeriv c (ρ - δ) - MP.mDeriv c (ρ + δ))

/-- **How small `δ` must be** (`notes/archive/rmt_R5.md`, step 7 and the numeric section). The
denominator bracket is vacuous for large `δ`, so the proof fixes `δ` before it divides. For
every target accuracy `ε` there is a `δ` with a positive denominator bracket, a negative upper
numerator endpoint, and both ends of the overlap bracket within `ε` of the limit `β`. -/
theorem exists_delta {θ c ρ β : ℝ} (hc : 0 < c) (hθ : 0 < θ) (hbρ : bulkEdge c < ρ)
    (hmneg : MP.m c ρ < 0)
    (hβ : θ * MP.m c ρ * (θ * MP.m c ρ) / ((θ ^ 2 + 1) * MP.mDeriv c ρ) = β)
    {ε : ℝ} (hε : 0 < ε) :
    ∃ δ : ℝ, 0 < δ ∧ δ < ρ - bulkEdge c ∧ numHiD θ c ρ δ < 0 ∧
      0 < denLoD θ c ρ δ ∧ 0 < denHiD θ c ρ δ ∧
      |numHiD θ c ρ δ * numHiD θ c ρ δ / denHiD θ c ρ δ - β| < ε ∧
      |numLoD θ c ρ δ * numLoD θ c ρ δ / denLoD θ c ρ δ - β| < ε := by
  have hdρ : 0 < MP.mDeriv c ρ := MP.mDeriv_pos hc hbρ
  have hden0 : 0 < (θ ^ 2 + 1) * MP.mDeriv c ρ := mul_pos (by positivity) hdρ
  have hwin : Set.Ioo (0 : ℝ) (ρ - bulkEdge c) ∈ 𝓝[>] (0 : ℝ) :=
    Ioo_mem_nhdsGT (by linarith)
  have hev : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < δ ∧ δ < ρ - bulkEdge c := hwin
  have hsub : Tendsto (fun δ : ℝ => ρ - δ) (𝓝[>] (0 : ℝ)) (𝓝 ρ) := by
    have h : Tendsto (fun δ : ℝ => ρ - δ) (𝓝 (0 : ℝ)) (𝓝 (ρ - 0)) :=
      tendsto_const_nhds.sub tendsto_id
    simpa using h.mono_left nhdsWithin_le_nhds
  have hadd : Tendsto (fun δ : ℝ => ρ + δ) (𝓝[>] (0 : ℝ)) (𝓝 ρ) := by
    have h : Tendsto (fun δ : ℝ => ρ + δ) (𝓝 (0 : ℝ)) (𝓝 (ρ + 0)) :=
      tendsto_const_nhds.add tendsto_id
    simpa using h.mono_left nhdsWithin_le_nhds
  have tmlo : Tendsto (fun δ => MP.m c (ρ - δ)) (𝓝[>] (0 : ℝ)) (𝓝 (MP.m c ρ)) := by
    refine Tendsto.comp ((MP.continuousOn_m hc) ρ (Set.mem_Ici.2 hbρ.le)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hsub, hev.mono fun δ h => Set.mem_Ici.2 (by linarith [h.2])⟩
  have tmhi : Tendsto (fun δ => MP.m c (ρ + δ)) (𝓝[>] (0 : ℝ)) (𝓝 (MP.m c ρ)) := by
    refine Tendsto.comp ((MP.continuousOn_m hc) ρ (Set.mem_Ici.2 hbρ.le)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hadd, hev.mono fun δ h => Set.mem_Ici.2 (by linarith [h.1])⟩
  have tdlo : Tendsto (fun δ => MP.mDeriv c (ρ - δ)) (𝓝[>] (0 : ℝ)) (𝓝 (MP.mDeriv c ρ)) := by
    refine Tendsto.comp ((MP.continuousOn_mDeriv hc) ρ (Set.mem_Ioi.2 hbρ)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hsub, hev.mono fun δ h => Set.mem_Ioi.2 (by linarith [h.2])⟩
  have tdhi : Tendsto (fun δ => MP.mDeriv c (ρ + δ)) (𝓝[>] (0 : ℝ)) (𝓝 (MP.mDeriv c ρ)) := by
    refine Tendsto.comp ((MP.continuousOn_mDeriv hc) ρ (Set.mem_Ioi.2 hbρ)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hadd, hev.mono fun δ h => Set.mem_Ioi.2 (by linarith [h.1])⟩
  have tnl : Tendsto (fun δ => numLoD θ c ρ δ) (𝓝[>] (0 : ℝ)) (𝓝 (θ * MP.m c ρ)) := by
    have h := (tmlo.const_mul θ).add (tmlo.sub tmhi)
    simp only [sub_self, add_zero] at h
    simpa only [numLoD] using h
  have tnh : Tendsto (fun δ => numHiD θ c ρ δ) (𝓝[>] (0 : ℝ)) (𝓝 (θ * MP.m c ρ)) := by
    have h := (tmhi.const_mul θ).add (tmhi.sub tmlo)
    simp only [sub_self, add_zero] at h
    simpa only [numHiD] using h
  have tdl : Tendsto (fun δ => denLoD θ c ρ δ) (𝓝[>] (0 : ℝ))
      (𝓝 ((θ ^ 2 + 1) * MP.mDeriv c ρ)) := by
    have h := (tdhi.const_mul (θ ^ 2 + 1)).add ((tdhi.sub tdlo).const_mul (2 * θ))
    simp only [sub_self, mul_zero, add_zero] at h
    simpa only [denLoD] using h
  have tdh : Tendsto (fun δ => denHiD θ c ρ δ) (𝓝[>] (0 : ℝ))
      (𝓝 ((θ ^ 2 + 1) * MP.mDeriv c ρ)) := by
    have h := (tdlo.const_mul (θ ^ 2 + 1)).add ((tdlo.sub tdhi).const_mul (2 * θ))
    simp only [sub_self, mul_zero, add_zero] at h
    simpa only [denHiD] using h
  have tlo : Tendsto (fun δ => numHiD θ c ρ δ * numHiD θ c ρ δ / denHiD θ c ρ δ)
      (𝓝[>] (0 : ℝ)) (𝓝 β) := by
    have h := (tnh.mul tnh).div tdh hden0.ne'
    rwa [hβ] at h
  have thi : Tendsto (fun δ => numLoD θ c ρ δ * numLoD θ c ρ δ / denLoD θ c ρ δ)
      (𝓝[>] (0 : ℝ)) (𝓝 β) := by
    have h := (tnl.mul tnl).div tdl hden0.ne'
    rwa [hβ] at h
  have hnegLim : θ * MP.m c ρ < 0 := mul_neg_of_pos_of_neg hθ hmneg
  have e1 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), numHiD θ c ρ δ < 0 := tnh.eventually_lt_const hnegLim
  have e2 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < denLoD θ c ρ δ := tdl.eventually_const_lt hden0
  have e3 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < denHiD θ c ρ δ := tdh.eventually_const_lt hden0
  have e4 : ∀ᶠ δ in 𝓝[>] (0 : ℝ),
      numHiD θ c ρ δ * numHiD θ c ρ δ / denHiD θ c ρ δ ∈ Metric.ball β ε :=
    tlo.eventually (Metric.ball_mem_nhds β hε)
  have e5 : ∀ᶠ δ in 𝓝[>] (0 : ℝ),
      numLoD θ c ρ δ * numLoD θ c ρ δ / denLoD θ c ρ δ ∈ Metric.ball β ε :=
    thi.eventually (Metric.ball_mem_nhds β hε)
  obtain ⟨δ, ⟨hδ0, hδb⟩, h1, h2, h3, h4, h5⟩ :=
    (hev.and (e1.and (e2.and (e3.and (e4.and e5))))).exists
  refine ⟨δ, hδ0, hδb, h1, h2, h3, ?_, ?_⟩
  · rw [Metric.mem_ball, Real.dist_eq] at h4; exact h4
  · rw [Metric.mem_ball, Real.dist_eq] at h5; exact h5

end R5

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

open R4

/-- **(H1) and (H2) in seven fields** (decision D6, 2026-08-30): the `edge` bound of item R3
and the six resolvent limits of items R1, R2 and T, at every real `z > bulkEdge c`. The six
`w` fields and `TendstoUnifOrth` of the 2026-08-29 draft are deleted; item Sym proves
`delocUniform` instead. Every field is about `m.W0`, `m.v` and `m.gvec` of item R0. -/
structure ResolventLimits (m : SpikedModel μ n d) (c : ℝ) : Prop where
  /-- (H1), from item R3 (`tendsto_measure_lamMax_le`). -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)
  /-- (H2), `vᵀ G₀ v → m`. -/
  vv : ∀ z > bulkEdge c,
    TendstoInProb μ (fun N ω => qform (m.W0 N ω) z (WithLp.ofLp (m.v N))) (MP.m c z)
  /-- (H2), `gᵀ G₀ g → m`. -/
  gg : ∀ z > bulkEdge c,
    TendstoInProb μ (fun N ω => qform (m.W0 N ω) z (m.gvec N ω)) (MP.m c z)
  /-- (H2), `vᵀ G₀ g → 0`. -/
  vg : ∀ z > bulkEdge c,
    TendstoInProb μ
      (fun N ω => cform (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)) 0
  /-- (H2), `vᵀ G₀² v → m'`. -/
  vv2 : ∀ z > bulkEdge c,
    TendstoInProb μ (fun N ω => qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N))) (MP.mDeriv c z)
  /-- (H2), `gᵀ G₀² g → m'`. -/
  gg2 : ∀ z > bulkEdge c,
    TendstoInProb μ (fun N ω => qform2 (m.W0 N ω) z (m.gvec N ω)) (MP.mDeriv c z)
  /-- (H2), `vᵀ G₀² g → 0`. -/
  vg2 : ∀ z > bulkEdge c,
    TendstoInProb μ
      (fun N ω => cform2 (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)) 0

/-! ### The four random bracket endpoints (step 6) -/

/-- Lower endpoint of the numerator bracket at the random root. -/
noncomputable def numLoR (m : SpikedModel μ n d) (zlo zhi : ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  m.θ * qform (m.W0 N ω) zlo (WithLp.ofLp (m.v N))
    + (qform (m.W0 N ω) zlo (WithLp.ofLp (m.v N) + m.gvec N ω)
        - qform (m.W0 N ω) zhi (WithLp.ofLp (m.v N))
        - qform (m.W0 N ω) zhi (m.gvec N ω)) / 2

/-- Upper endpoint of the numerator bracket at the random root. -/
noncomputable def numHiR (m : SpikedModel μ n d) (zlo zhi : ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  m.θ * qform (m.W0 N ω) zhi (WithLp.ofLp (m.v N))
    + (qform (m.W0 N ω) zhi (WithLp.ofLp (m.v N) + m.gvec N ω)
        - qform (m.W0 N ω) zlo (WithLp.ofLp (m.v N))
        - qform (m.W0 N ω) zlo (m.gvec N ω)) / 2

/-- Lower endpoint of the denominator bracket at the random root. -/
noncomputable def denLoR (m : SpikedModel μ n d) (zlo zhi : ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  m.θ ^ 2 * qform2 (m.W0 N ω) zhi (WithLp.ofLp (m.v N))
    + m.θ * (qform2 (m.W0 N ω) zhi (WithLp.ofLp (m.v N) + m.gvec N ω)
        - qform2 (m.W0 N ω) zlo (WithLp.ofLp (m.v N))
        - qform2 (m.W0 N ω) zlo (m.gvec N ω))
    + qform2 (m.W0 N ω) zhi (m.gvec N ω)

/-- Upper endpoint of the denominator bracket at the random root. -/
noncomputable def denHiR (m : SpikedModel μ n d) (zlo zhi : ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  m.θ ^ 2 * qform2 (m.W0 N ω) zlo (WithLp.ofLp (m.v N))
    + m.θ * (qform2 (m.W0 N ω) zlo (WithLp.ofLp (m.v N) + m.gvec N ω)
        - qform2 (m.W0 N ω) zhi (WithLp.ofLp (m.v N))
        - qform2 (m.W0 N ω) zhi (m.gvec N ω))
    + qform2 (m.W0 N ω) zlo (m.gvec N ω)

/-! ### Measurability of the edge event

The `edge` field is stated with `→ 𝓝 1`, the shape item R3 proves. R5 needs the complement,
so it needs the event to be measurable. `W₀ = E⊥ᵀ E⊥` makes `λ_max(W₀)` the `gramLamMax` of
`E⊥`, which is a measurable function of the noise (`m.hZ`). -/

theorem measurable_Eperp (m : SpikedModel μ n d) (N : ℕ) : Measurable (m.Eperp N) := by
  have hZ : ∀ (r : Fin (n N)) (j : Fin (d N)), Measurable fun ω => m.Z N ω r j := by
    intro r j
    have h1 : Measurable fun ω => m.Z N ω r := (measurable_pi_apply r).comp (m.hZ N)
    exact (measurable_pi_apply j).comp h1
  have hE : ∀ (r : Fin (n N)) (j : Fin (d N)), Measurable fun ω => m.E N ω r j := by
    intro r j
    have h : (fun ω => m.E N ω r j) = fun ω => (Real.sqrt (d N))⁻¹ * m.Z N ω r j := rfl
    rw [h]
    exact (hZ r j).const_mul _
  have hg : ∀ j : Fin (d N), Measurable fun ω => m.gvec N ω j := by
    intro j
    have h : (fun ω => m.gvec N ω j)
        = fun ω => ∑ k : Fin (n N), m.E N ω k j * WithLp.ofLp (m.u N) k := rfl
    rw [h]
    exact Finset.measurable_sum _ fun k _ => (hE k j).mul_const _
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.Eperp N ω r j)
      = fun ω => m.E N ω r j - WithLp.ofLp (m.u N) r * m.gvec N ω j := rfl
  rw [h]
  exact (hE r j).sub ((hg j).const_mul _)

theorem measurableSet_lamMax_W0_le (m : SpikedModel μ n d) (N : ℕ) (r : ℝ) :
    MeasurableSet {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ r} := by
  have h : {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ r}
      = (fun ω => gramLamMax (m.Eperp N ω)) ⁻¹' Set.Iic r := rfl
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_Eperp N)) measurableSet_Iic

/-! ### Localization at the secular root (step 4) -/

/-- **Step 4, one `ω`.** If `λ_max(W₀) < zlo`, the secular function is negative at `zlo` and
positive at `zhi`, then `q ≠ 0`, the top eigenvalue of `Xᵀ X` is the secular root inside
`(zlo, zhi)`, and it is simple. -/
theorem localize (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {zlo zhi : ℝ}
    (h₀ : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) < zlo) (hz : zlo < zhi)
    (hslo : secular (m.W0 N ω) (m.qvec N ω) zlo < 0)
    (hshi : 0 < secular (m.W0 N ω) (m.qvec N ω) zhi) :
    m.qvec N ω ≠ 0 ∧ zlo < gramLamMax (m.X N ω) ∧ gramLamMax (m.X N ω) < zhi ∧
      secular (m.W0 N ω) (m.qvec N ω) (gramLamMax (m.X N ω)) = 0 ∧
      TopSimple ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω)) := by
  have hW := m.isHermitian_W0 N ω
  have hq : m.qvec N ω ≠ 0 := by
    intro h0
    have h1 : secular (m.W0 N ω) (m.qvec N ω) zlo = 1 := by
      change 1 + qform (m.W0 N ω) zlo (m.qvec N ω) = 1
      rw [h0]
      simp [qform]
    linarith
  have hcont : ContinuousOn (fun t => secular (m.W0 N ω) (m.qvec N ω) t) (Set.Icc zlo zhi) := by
    intro s hs
    have hs' : lamMax (m.W0 N ω) hW < s := lt_of_lt_of_le h₀ hs.1
    exact (continuousAt_const.add
      (hasDerivAt_qform (y := m.qvec N ω) hW hs').continuousAt).continuousWithinAt
  obtain ⟨lam, hlmem, hlz⟩ :=
    intermediate_value_Ioo hz.le hcont (Set.mem_Ioo.2 ⟨hslo, hshi⟩)
  have hlamlt : lamMax (m.W0 N ω) hW < lam := lt_trans h₀ hlmem.1
  have hA' : (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)).IsHermitian := by
    rw [← m.gram_eq N ω]; exact isHermitian_transpose_mul_self _
  have hgl : gramLamMax (m.X N ω)
      = lamMax (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hA' :=
    lamMax_congr (m.gram_eq N ω) _ hA'
  have hlamMax : gramLamMax (m.X N ω) = lam := by
    rw [hgl, lamMax_eq hW hq hlamlt hlz hA']
  refine ⟨hq, ?_, ?_, ?_, ?_⟩
  · rw [hlamMax]; exact hlmem.1
  · rw [hlamMax]; exact hlmem.2
  · rw [hlamMax]; exact hlz
  · exact topSimple_congr (m.gram_eq N ω) _ hA' (topSimple hW hq hlamlt hlz hA')

/-- **Step 7, one `ω`.** On the localization event the overlap sits inside the random
bracket built from the forms at `zlo` and `zhi`. -/
theorem overlap_in_bracket (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {zlo zhi : ℝ}
    (h₀ : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) < zlo) (hz : zlo < zhi)
    (hslo : secular (m.W0 N ω) (m.qvec N ω) zlo < 0)
    (hshi : 0 < secular (m.W0 N ω) (m.qvec N ω) zhi)
    (hnh : m.numHiR zlo zhi N ω < 0) (hdl : 0 < m.denLoR zlo zhi N ω) :
    m.numHiR zlo zhi N ω * m.numHiR zlo zhi N ω / m.denHiR zlo zhi N ω
        ≤ overlap (m.X N ω) (m.v N) ∧
      overlap (m.X N ω) (m.v N)
        ≤ m.numLoR zlo zhi N ω * m.numLoR zlo zhi N ω / m.denLoR zlo zhi N ω := by
  obtain ⟨hq, hlo, hhi, hsec, -⟩ := m.localize N ω h₀ hz hslo hshi
  have hW := m.isHermitian_W0 N ω
  have hA' : (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)).IsHermitian := by
    rw [← m.gram_eq N ω]; exact isHermitian_transpose_mul_self _
  have hlamlt : lamMax (m.W0 N ω) hW < gramLamMax (m.X N ω) := lt_trans h₀ hlo
  have hqdef : m.qvec N ω = m.θ • WithLp.ofLp (m.v N) + m.gvec N ω := rfl
  -- the overlap formula of R4d.7
  have hov : overlap (m.X N ω) (m.v N)
      = cform (m.W0 N ω) (gramLamMax (m.X N ω)) (WithLp.ofLp (m.v N)) (m.qvec N ω)
          * cform (m.W0 N ω) (gramLamMax (m.X N ω)) (WithLp.ofLp (m.v N)) (m.qvec N ω)
        / qform2 (m.W0 N ω) (gramLamMax (m.X N ω)) (m.qvec N ω) := by
    have h1 : overlap (m.X N ω) (m.v N)
        = ‖topProj (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hA'
            (m.v N)‖ ^ 2 := by
      change ‖topProj ((m.X N ω)ᵀ * m.X N ω)
        (isHermitian_transpose_mul_self (m.X N ω)) (m.v N)‖ ^ 2 = _
      rw [topProj_congr (m.gram_eq N ω) (isHermitian_transpose_mul_self (m.X N ω)) hA'
        (m.v N)]
    rw [h1, topProj_norm_sq hW hq hlamlt hsec hA' (m.v N)]
    simp only [cform, qform2]
    ring
  rw [hov, hqdef]
  refine R5.div_bracket ?_ ?_ ?_ ?_ ?_ hdl
  · exact (num_bracket hW m.hθ h₀ hlo.le hhi.le (WithLp.ofLp (m.v N)) (m.gvec N ω)).1
  · exact (num_bracket hW m.hθ h₀ hlo.le hhi.le (WithLp.ofLp (m.v N)) (m.gvec N ω)).2
  · exact hnh
  · exact (den_bracket hW m.hθ h₀ hlo.le hhi.le (WithLp.ofLp (m.v N)) (m.gvec N ω)).1
  · exact (den_bracket hW m.hθ h₀ hlo.le hhi.le (WithLp.ofLp (m.v N)) (m.gvec N ω)).2

/-- **Step 9, one `ω`.** If the secular function is positive at `z₀ > λ_max(W₀)` then the top
eigenvalue of `Xᵀ X` is at most `z₀`. -/
theorem gramLamMax_le_of_secular_pos (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {z₀ : ℝ}
    (h₀ : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) < z₀)
    (hs : 0 < secular (m.W0 N ω) (m.qvec N ω) z₀) : gramLamMax (m.X N ω) ≤ z₀ := by
  by_contra hcon
  rw [not_le] at hcon
  have hW := m.isHermitian_W0 N ω
  have hA' : (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)).IsHermitian := by
    rw [← m.gram_eq N ω]; exact isHermitian_transpose_mul_self _
  have hgl : gramLamMax (m.X N ω)
      = lamMax (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hA' :=
    lamMax_congr (m.gram_eq N ω) _ hA'
  have hlt : lamMax (m.W0 N ω) hW < gramLamMax (m.X N ω) := lt_trans h₀ hcon
  obtain ⟨j, hj⟩ := exists_eigenvalues_eq_lamMax hA' (m.hd N)
  have hspec : lamMax (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hA'
      ∈ spectrum ℝ (toOp (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω))) := by
    rw [spectrum_toOp]
    exact hj ▸ hA'.eigenvalues_mem_spectrum_real j
  have hzero : secular (m.W0 N ω) (m.qvec N ω) (gramLamMax (m.X N ω)) = 0 := by
    rw [hgl]
    exact (secular_eq_zero_iff hW (hgl ▸ hlt)).mpr hspec
  rcases eq_or_ne (m.qvec N ω) 0 with hq0 | hq
  · have h1 : secular (m.W0 N ω) (m.qvec N ω) (gramLamMax (m.X N ω)) = 1 := by
      change 1 + qform (m.W0 N ω) (gramLamMax (m.X N ω)) (m.qvec N ω) = 1
      rw [hq0]
      simp [qform]
    rw [h1] at hzero
    norm_num at hzero
  · have hmono := secular_strictMonoOn hW hq (Set.mem_Ioi.2 h₀) (Set.mem_Ioi.2 hlt) hcon
    linarith

/-! ### Measure helpers

The four helpers this file used to hold (`tendsto_measure_zero_of_subset`,
`tendsto_measure_zero_union`, `tendsto_measure_compl_zero`,
`tendsto_measure_one_of_bad`) are public in `Prob/TendstoInProb.lean` since cleanup
wave 1 (2026-08-30). `tendsto_measure_compl_zero` there takes `NullMeasurableSet`,
so a `MeasurableSet` proof reaches it through `MeasurableSet.nullMeasurableSet`. -/

variable [∀ N, IsProbabilityMeasure (μ N)] {m : SpikedModel μ n d} {c : ℝ}

set_option linter.unusedSectionVars false

/-! ### Step 2 and step 3: limits of the random forms -/

theorem tendstoInProb_qform_add (H : m.ResolventLimits c) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ
      (fun N ω => qform (m.W0 N ω) z (WithLp.ofLp (m.v N) + m.gvec N ω)) (2 * MP.m c z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N),
      qform (m.W0 N ω) z (WithLp.ofLp (m.v N) + m.gvec N ω)
        = qform (m.W0 N ω) z (WithLp.ofLp (m.v N))
            + 2 * cform (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
          + qform (m.W0 N ω) z (m.gvec N ω) :=
    fun N ω => qform_add (m.isHermitian_W0 N ω) z _ _
  have h := ((H.vv z hz).add ((H.vg z hz).const_mul 2)).add (H.gg z hz)
  have hlim : MP.m c z + 2 * (0 : ℝ) + MP.m c z = 2 * MP.m c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

theorem tendstoInProb_qform2_add (H : m.ResolventLimits c) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ
      (fun N ω => qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N) + m.gvec N ω))
      (2 * MP.mDeriv c z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N),
      qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N) + m.gvec N ω)
        = qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N))
            + 2 * cform2 (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
          + qform2 (m.W0 N ω) z (m.gvec N ω) :=
    fun N ω => qform2_add (m.isHermitian_W0 N ω) z _ _
  have h := ((H.vv2 z hz).add ((H.vg2 z hz).const_mul 2)).add (H.gg2 z hz)
  have hlim : MP.mDeriv c z + 2 * (0 : ℝ) + MP.mDeriv c z = 2 * MP.mDeriv c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-- **Step 3.** The random secular function converges at every fixed `z > b`. -/
theorem tendstoInProb_secular (H : m.ResolventLimits c) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => secular (m.W0 N ω) (m.qvec N ω) z)
      (1 + (m.θ ^ 2 + 1) * MP.m c z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N),
      secular (m.W0 N ω) (m.qvec N ω) z
        = 1 + m.θ ^ 2 * qform (m.W0 N ω) z (WithLp.ofLp (m.v N))
            + 2 * m.θ * cform (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
          + qform (m.W0 N ω) z (m.gvec N ω) := by
    intro N ω
    change 1 + qform (m.W0 N ω) z (m.qvec N ω) = _
    rw [show m.qvec N ω = m.θ • WithLp.ofLp (m.v N) + m.gvec N ω from rfl,
      qform_smul_add (m.isHermitian_W0 N ω)]
    ring
  have h := (((TendstoInProb.const μ (1 : ℝ)).add ((H.vv z hz).const_mul (m.θ ^ 2))).add
    ((H.vg z hz).const_mul (2 * m.θ))).add (H.gg z hz)
  have hlim : (1 : ℝ) + m.θ ^ 2 * MP.m c z + 2 * m.θ * 0 + MP.m c z
      = 1 + (m.θ ^ 2 + 1) * MP.m c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-! ### Limits of the four random bracket endpoints -/

theorem tendstoInProb_numLoR (H : m.ResolventLimits c) {ρ δ : ℝ}
    (hlo : bulkEdge c < ρ - δ) (hhi : bulkEdge c < ρ + δ) :
    TendstoInProb μ (m.numLoR (ρ - δ) (ρ + δ)) (R5.numLoD m.θ c ρ δ) := by
  have h := ((H.vv (ρ - δ) hlo).const_mul m.θ).add
    ((((tendstoInProb_qform_add H hlo).sub (H.vv (ρ + δ) hhi)).sub (H.gg (ρ + δ) hhi)).div
      (TendstoInProb.const μ (2 : ℝ)) two_ne_zero)
  have hlim : m.θ * MP.m c (ρ - δ)
      + (2 * MP.m c (ρ - δ) - MP.m c (ρ + δ) - MP.m c (ρ + δ)) / 2
      = R5.numLoD m.θ c ρ δ := by unfold R5.numLoD; ring
  rw [hlim] at h
  exact h

theorem tendstoInProb_numHiR (H : m.ResolventLimits c) {ρ δ : ℝ}
    (hlo : bulkEdge c < ρ - δ) (hhi : bulkEdge c < ρ + δ) :
    TendstoInProb μ (m.numHiR (ρ - δ) (ρ + δ)) (R5.numHiD m.θ c ρ δ) := by
  have h := ((H.vv (ρ + δ) hhi).const_mul m.θ).add
    ((((tendstoInProb_qform_add H hhi).sub (H.vv (ρ - δ) hlo)).sub (H.gg (ρ - δ) hlo)).div
      (TendstoInProb.const μ (2 : ℝ)) two_ne_zero)
  have hlim : m.θ * MP.m c (ρ + δ)
      + (2 * MP.m c (ρ + δ) - MP.m c (ρ - δ) - MP.m c (ρ - δ)) / 2
      = R5.numHiD m.θ c ρ δ := by unfold R5.numHiD; ring
  rw [hlim] at h
  exact h

theorem tendstoInProb_denLoR (H : m.ResolventLimits c) {ρ δ : ℝ}
    (hlo : bulkEdge c < ρ - δ) (hhi : bulkEdge c < ρ + δ) :
    TendstoInProb μ (m.denLoR (ρ - δ) (ρ + δ)) (R5.denLoD m.θ c ρ δ) := by
  have h := (((H.vv2 (ρ + δ) hhi).const_mul (m.θ ^ 2)).add
    ((((tendstoInProb_qform2_add H hhi).sub (H.vv2 (ρ - δ) hlo)).sub
      (H.gg2 (ρ - δ) hlo)).const_mul m.θ)).add (H.gg2 (ρ + δ) hhi)
  have hlim : m.θ ^ 2 * MP.mDeriv c (ρ + δ)
      + m.θ * (2 * MP.mDeriv c (ρ + δ) - MP.mDeriv c (ρ - δ) - MP.mDeriv c (ρ - δ))
      + MP.mDeriv c (ρ + δ) = R5.denLoD m.θ c ρ δ := by unfold R5.denLoD; ring
  rw [hlim] at h
  exact h

theorem tendstoInProb_denHiR (H : m.ResolventLimits c) {ρ δ : ℝ}
    (hlo : bulkEdge c < ρ - δ) (hhi : bulkEdge c < ρ + δ) :
    TendstoInProb μ (m.denHiR (ρ - δ) (ρ + δ)) (R5.denHiD m.θ c ρ δ) := by
  have h := (((H.vv2 (ρ - δ) hlo).const_mul (m.θ ^ 2)).add
    ((((tendstoInProb_qform2_add H hlo).sub (H.vv2 (ρ + δ) hhi)).sub
      (H.gg2 (ρ + δ) hhi)).const_mul m.θ)).add (H.gg2 (ρ - δ) hlo)
  have hlim : m.θ ^ 2 * MP.mDeriv c (ρ - δ)
      + m.θ * (2 * MP.mDeriv c (ρ - δ) - MP.mDeriv c (ρ + δ) - MP.mDeriv c (ρ + δ))
      + MP.mDeriv c (ρ - δ) = R5.denHiD m.θ c ρ δ := by unfold R5.denHiD; ring
  rw [hlim] at h
  exact h

/-! ### The supercritical conclusions -/

private theorem theta_pos (hc : 0 < c) (hθ : c < m.θ ^ 4) : 0 < m.θ := by
  rcases m.hθ.lt_or_eq with h | h
  · exact h
  · exfalso
    have h4 : m.θ ^ 4 = 0 := by rw [← h]; norm_num
    linarith

/-- **The localization bad set.** Its measure tends to `0`, and off it the three hypotheses of
`localize` hold at `zlo = ρ - δ` and `zhi = ρ + δ`. Every supercritical conclusion uses this
one lemma; (H1) enters here and nowhere else. -/
theorem exists_localization_bad (H : m.ResolventLimits c) (hc : 0 < c) (hθ : c < m.θ ^ 4)
    {δ : ℝ} (hδ : 0 < δ) (hδ' : δ < rhoSq m.θ c - bulkEdge c) :
    ∃ B : ∀ N, Set (Ω N), Tendsto (fun N => μ N (B N)) atTop (𝓝 0) ∧
      ∀ (N : ℕ) (ω : Ω N), ω ∉ B N →
        lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) < rhoSq m.θ c - δ ∧
        secular (m.W0 N ω) (m.qvec N ω) (rhoSq m.θ c - δ) < 0 ∧
        0 < secular (m.W0 N ω) (m.qvec N ω) (rhoSq m.θ c + δ) := by
  have hθ0 : 0 < m.θ := theta_pos hc hθ
  have hbρ : bulkEdge c < rhoSq m.θ c := MP.bulkEdge_lt_rhoSq hc hθ0 hθ
  have hlo : bulkEdge c < rhoSq m.θ c - δ := by linarith
  have hhi : bulkEdge c < rhoSq m.θ c + δ := by linarith
  have hpos : (0 : ℝ) < m.θ ^ 2 + 1 := by positivity
  have hzero : 1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c) = 0 :=
    (MP.secular_eq_zero_iff hc hθ0 hθ hbρ).mpr rfl
  have hFlo : 1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c - δ) < 0 := by
    have hmono := MP.m_strictMonoOn hc (Set.mem_Ici.2 hlo.le) (Set.mem_Ici.2 hbρ.le)
      (by linarith)
    have := mul_lt_mul_of_pos_left hmono hpos
    linarith
  have hFhi : 0 < 1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c + δ) := by
    have hmono := MP.m_strictMonoOn hc (Set.mem_Ici.2 hbρ.le) (Set.mem_Ici.2 hhi.le)
      (by linarith)
    have := mul_lt_mul_of_pos_left hmono hpos
    linarith
  have hε₀ : (0 : ℝ) < (rhoSq m.θ c - δ - bulkEdge c) / 2 := by linarith
  have hB1 : Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω)
        ≤ bulkEdge c + (rhoSq m.θ c - δ - bulkEdge c) / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet)
      (H.edge _ hε₀)
  have hB2 := tendstoInProb_secular H hlo
    (-(1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c - δ)) / 2) (by linarith)
  have hB3 := tendstoInProb_secular H hhi
    ((1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c + δ)) / 2) (by linarith)
  refine ⟨_, tendsto_measure_zero_union hB1 (tendsto_measure_zero_union hB2 hB3), ?_⟩
  intro N ω hω
  simp only [Set.mem_union, not_or] at hω
  obtain ⟨g1, g2, g3⟩ := hω
  have hlam : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω)
      ≤ bulkEdge c + (rhoSq m.θ c - δ - bulkEdge c) / 2 := by
    by_contra hx; exact g1 (Set.mem_compl hx)
  refine ⟨by linarith, ?_, ?_⟩
  · by_contra hx
    rw [not_lt] at hx
    refine g2 ?_
    change -(1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c - δ)) / 2
      ≤ |secular (m.W0 N ω) (m.qvec N ω) (rhoSq m.θ c - δ)
        - (1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c - δ))|
    rw [abs_of_nonneg (by linarith)]
    linarith
  · by_contra hx
    rw [not_lt] at hx
    refine g3 ?_
    change (1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c + δ)) / 2
      ≤ |secular (m.W0 N ω) (m.qvec N ω) (rhoSq m.θ c + δ)
        - (1 + (m.θ ^ 2 + 1) * MP.m c (rhoSq m.θ c + δ))|
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- **Conclusion 1.** `SingleTableLaw.lamMax`, supercritical branch. -/
theorem lamMax_tendstoInProb (H : m.ResolventLimits c) (hc : 0 < c) (hθ : c < m.θ ^ 4) :
    TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) (rhoSq m.θ c) := by
  have hθ0 : 0 < m.θ := theta_pos hc hθ
  have hbρ : bulkEdge c < rhoSq m.θ c := MP.bulkEdge_lt_rhoSq hc hθ0 hθ
  intro ε hε
  obtain ⟨δ, hδ, hδ', hδε⟩ :
      ∃ δ : ℝ, 0 < δ ∧ δ < rhoSq m.θ c - bulkEdge c ∧ δ ≤ ε :=
    ⟨min ε ((rhoSq m.θ c - bulkEdge c) / 2), lt_min hε (by linarith),
      lt_of_le_of_lt (min_le_right _ _) (by linarith), min_le_left _ _⟩
  obtain ⟨B, hB, hBgood⟩ := exists_localization_bad H hc hθ hδ hδ'
  refine tendsto_measure_zero_of_subset ?_ hB
  intro N ω hω
  by_contra hnb
  obtain ⟨hlam, hslo, hshi⟩ := hBgood N ω hnb
  obtain ⟨-, hlo, hhi, -, -⟩ := m.localize N ω hlam (by linarith) hslo hshi
  have h2 : ε ≤ |gramLamMax (m.X N ω) - rhoSq m.θ c| := hω
  have h3 : |gramLamMax (m.X N ω) - rhoSq m.θ c| < δ := by
    rw [abs_lt]; constructor <;> linarith
  linarith

/-- **Conclusion 3.** `SingleTableLaw.align`. The `ε`-`δ` order matters: `δ` is fixed first,
small enough that the deterministic bracket is inside `ε/2` of `betaSq θ c` and that its
denominator end is positive (`R5.exists_delta`); only then is the limit in `N` taken, at the
two points `ρ ∓ δ`. -/
theorem align_tendstoInProb (H : m.ResolventLimits c) (hc : 0 < c) (hθ : c < m.θ ^ 4) :
    TendstoInProb μ (fun N ω => overlap (m.X N ω) (m.v N)) (betaSq m.θ c) := by
  have hθ0 : 0 < m.θ := theta_pos hc hθ
  have hbρ : bulkEdge c < rhoSq m.θ c := MP.bulkEdge_lt_rhoSq hc hθ0 hθ
  have hmneg : MP.m c (rhoSq m.θ c) < 0 := MP.m_neg hc hbρ
  have hβ : m.θ * MP.m c (rhoSq m.θ c) * (m.θ * MP.m c (rhoSq m.θ c))
      / ((m.θ ^ 2 + 1) * MP.mDeriv c (rhoSq m.θ c)) = betaSq m.θ c := by
    rw [← MP.overlap_identity hc hθ0 hθ]; ring
  intro ε hε
  obtain ⟨δ, hδ, hδ', hnhD, hdlD, hdhD, hLoD, hHiD⟩ :=
    R5.exists_delta hc hθ0 hbρ hmneg hβ (half_pos hε)
  have hlo : bulkEdge c < rhoSq m.θ c - δ := by linarith
  have hhi : bulkEdge c < rhoSq m.θ c + δ := by linarith
  obtain ⟨B, hB, hBgood⟩ := exists_localization_bad H hc hθ hδ hδ'
  have hB4 := tendstoInProb_numHiR H hlo hhi
    (-R5.numHiD m.θ c (rhoSq m.θ c) δ / 2) (by linarith)
  have hB5 := tendstoInProb_denLoR H hlo hhi
    (R5.denLoD m.θ c (rhoSq m.θ c) δ / 2) (by linarith)
  have hB6 := ((tendstoInProb_numHiR H hlo hhi).mul (tendstoInProb_numHiR H hlo hhi)).div
    (tendstoInProb_denHiR H hlo hhi) hdhD.ne' (ε / 2) (half_pos hε)
  have hB7 := ((tendstoInProb_numLoR H hlo hhi).mul (tendstoInProb_numLoR H hlo hhi)).div
    (tendstoInProb_denLoR H hlo hhi) hdlD.ne' (ε / 2) (half_pos hε)
  refine tendsto_measure_zero_of_subset ?_
    (tendsto_measure_zero_union hB (tendsto_measure_zero_union hB4
      (tendsto_measure_zero_union hB5 (tendsto_measure_zero_union hB6 hB7))))
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g0, g4, g5, g6, g7⟩ := hcon
  obtain ⟨hlam, hslo, hshi⟩ := hBgood N ω g0
  have g4' : |m.numHiR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
      - R5.numHiD m.θ c (rhoSq m.θ c) δ| < -R5.numHiD m.θ c (rhoSq m.θ c) δ / 2 :=
    not_le.mp g4
  have g5' : |m.denLoR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
      - R5.denLoD m.θ c (rhoSq m.θ c) δ| < R5.denLoD m.θ c (rhoSq m.θ c) δ / 2 :=
    not_le.mp g5
  have g6' : |m.numHiR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
        * m.numHiR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
        / m.denHiR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
      - R5.numHiD m.θ c (rhoSq m.θ c) δ * R5.numHiD m.θ c (rhoSq m.θ c) δ
        / R5.denHiD m.θ c (rhoSq m.θ c) δ| < ε / 2 := not_le.mp g6
  have g7' : |m.numLoR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
        * m.numLoR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
        / m.denLoR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω
      - R5.numLoD m.θ c (rhoSq m.θ c) δ * R5.numLoD m.θ c (rhoSq m.θ c) δ
        / R5.denLoD m.θ c (rhoSq m.θ c) δ| < ε / 2 := not_le.mp g7
  -- the two sign conditions of the sandwich
  have hnh : m.numHiR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω < 0 := by
    have := (abs_lt.mp g4').2
    linarith
  have hdl : 0 < m.denLoR (rhoSq m.θ c - δ) (rhoSq m.θ c + δ) N ω := by
    have := (abs_lt.mp g5').1
    linarith
  obtain ⟨hlow, hhigh⟩ :=
    m.overlap_in_bracket N ω hlam (by linarith) hslo hshi hnh hdl
  have h6 := abs_lt.mp g6'
  have h7 := abs_lt.mp g7'
  have h8 := abs_lt.mp hLoD
  have h9 := abs_lt.mp hHiD
  have hfin : |overlap (m.X N ω) (m.v N) - betaSq m.θ c| < ε := by
    rw [abs_lt]
    constructor <;> linarith
  exact absurd (show ε ≤ |overlap (m.X N ω) (m.v N) - betaSq m.θ c| from hω)
    (not_le.mpr hfin)

/-- **Subcritical.** The upper edge bound; the lower bound is R3⁻ and the overlap is R6'. -/
theorem tendsto_measure_lamMax_le_of_subcritical (H : m.ResolventLimits c) (hc : 0 < c)
    (hθ : m.θ ^ 4 ≤ c) : ∀ ε > 0,
    Tendsto (fun N => μ N {ω | gramLamMax (m.X N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  intro ε hε
  have hz₀ : bulkEdge c < bulkEdge c + ε := by linarith
  have hF : 0 < 1 + (m.θ ^ 2 + 1) * MP.m c (bulkEdge c + ε) :=
    MP.secular_pos_of_subcritical hc m.hθ hθ hz₀
  have hε₀ : (0 : ℝ) < ε / 2 := by linarith
  have hB1 : Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet)
      (H.edge (ε / 2) hε₀)
  have hB2 := tendstoInProb_secular H hz₀
    ((1 + (m.θ ^ 2 + 1) * MP.m c (bulkEdge c + ε)) / 2) (by linarith)
  refine tendsto_measure_one_of_bad ?_ (tendsto_measure_zero_union hB1 hB2)
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g1, g2⟩ := hcon
  have hlam : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε / 2 := by
    by_contra hx; exact g1 (Set.mem_compl hx)
  have hlam' : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) < bulkEdge c + ε := by linarith
  have hs : 0 < secular (m.W0 N ω) (m.qvec N ω) (bulkEdge c + ε) := by
    by_contra hx
    rw [not_lt] at hx
    refine g2 ?_
    change (1 + (m.θ ^ 2 + 1) * MP.m c (bulkEdge c + ε)) / 2
      ≤ |secular (m.W0 N ω) (m.qvec N ω) (bulkEdge c + ε)
        - (1 + (m.θ ^ 2 + 1) * MP.m c (bulkEdge c + ε))|
    rw [abs_of_nonpos (by linarith)]
    linarith
  exact hω (m.gramLamMax_le_of_secular_pos N ω hlam' hs)

end SpikedModel

end StackedSVD
