/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# Tensorization of the variance

This file proves the Efron-Stein tensorization inequality for the variance of a square
integrable function on a finite product measure: the variance of `f` is at most the sum,
over the coordinates, of the mean of the variance in that coordinate alone. Mathlib has no
Efron-Stein inequality and no martingale-difference `L²` identity, so this file gives the
substitute needed for the leave-one-out step of the trace law (Bai and Silverstein, section
3.3), by direct induction instead of a martingale argument.

The proof is one induction on the number of coordinates, in the style of
`Prob/LinFormMoments.lean`. The step splits off one coordinate with the measure preserving
equivalence `MeasureTheory.measurePreserving_piFinSuccAbove`, applies the one-coordinate law
of total variance in Fubini form, and recurses on the mean over the split coordinate. The two
ingredients of the recursion step are proved once, for a general two-factor product measure:
`total_var_eq` (the law of total variance, an equality) and `var_mean_le_mean_var` (the
variance of a conditional mean is at most the mean of the conditional variance, a form of
Jensen's inequality applied twice). The induction itself is carried out for `ι = Fin D` and
then transported to a general `Fintype` by the reindexing equivalence
`MeasureTheory.measurePreserving_piCongrLeft`, in the style of
`MeasureTheory.Integrable.fintype_prod_dep`.

The file ends with the bounded-difference corollary: if replacing one coordinate moves `f` by
at most a fixed amount, the variance is at most the sum of the squares of those amounts. The
constant proved is `1`, not the sharper `1/2` of the classical bounded-differences inequality;
the sharper constant needs recentering at an independent copy rather than at a fixed point,
which is not used here.
-/

open MeasureTheory ProbabilityTheory Filter

namespace StackedSVD
namespace Tensorization

/-! ### Variance algebra on a single probability space -/

/-- `Var(f) = E[f²] - (E f)²`, from bare integrability of `f` and `f²` (not `MemLp`). -/
private theorem integral_sq_sub_sq {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {f : β → ℝ} (hf1 : Integrable f μ)
    (hf2 : Integrable (fun x => f x ^ 2) μ) :
    ∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ = (∫ x, f x ^ 2 ∂μ) - (∫ x, f x ∂μ) ^ 2 := by
  set m := ∫ y, f y ∂μ with hm_def
  have hc : Integrable (fun x : β => 2 * m * f x) μ := hf1.const_mul _
  have hstep : ∫ x, (f x - m) ^ 2 ∂μ = ∫ x, (f x ^ 2 - 2 * m * f x + m ^ 2) ∂μ :=
    integral_congr_ae (Filter.Eventually.of_forall fun x => by ring)
  have h1 : ∫ x, (f x ^ 2 - 2 * m * f x + m ^ 2) ∂μ
      = (∫ x, (f x ^ 2 - 2 * m * f x) ∂μ) + ∫ x, m ^ 2 ∂μ := by
    simpa using integral_add (hf2.sub hc) (integrable_const (m ^ 2))
  have h2 : ∫ x, (f x ^ 2 - 2 * m * f x) ∂μ = (∫ x, f x ^ 2 ∂μ) - ∫ x, 2 * m * f x ∂μ := by
    simpa using integral_sub hf2 hc
  have h3 : ∫ x, 2 * m * f x ∂μ = 2 * m * ∫ x, f x ∂μ := integral_const_mul (2 * m) f
  have h4 : ∫ x, m ^ 2 ∂μ = m ^ 2 := by simp
  rw [hstep, h1, h2, h3, h4]
  ring

/-- The same, from `MemLp`. -/
private theorem memLp_sq_sub_sq {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {f : β → ℝ} (hf : MemLp f 2 μ) :
    ∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ = (∫ x, f x ^ 2 ∂μ) - (∫ x, f x ∂μ) ^ 2 :=
  integral_sq_sub_sq (hf.integrable (by norm_num)) hf.integrable_sq

/-- Jensen: `(E f)² ≤ E[f²]`. -/
private theorem sq_integral_le {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {f : β → ℝ} (hf1 : Integrable f μ)
    (hf2 : Integrable (fun x => f x ^ 2) μ) :
    (∫ x, f x ∂μ) ^ 2 ≤ ∫ x, f x ^ 2 ∂μ := by
  have h0 : 0 ≤ ∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ := integral_nonneg (fun x => sq_nonneg _)
  rw [integral_sq_sub_sq hf1 hf2] at h0
  linarith

/-- The mean minimizes the mean square deviation: `E[(f - a)²] = Var(f) + (E f - a)²` for any
fixed point `a`. -/
private theorem integral_sq_sub_const_eq {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {f : β → ℝ} (hf1 : Integrable f μ)
    (hf2 : Integrable (fun x => f x ^ 2) μ) (a : ℝ) :
    ∫ x, (f x - a) ^ 2 ∂μ = (∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ) + (∫ y, f y ∂μ - a) ^ 2 := by
  have hca : Integrable (fun x : β => 2 * a * f x) μ := hf1.const_mul _
  have hstep : ∫ x, (f x - a) ^ 2 ∂μ = ∫ x, (f x ^ 2 - 2 * a * f x + a ^ 2) ∂μ :=
    integral_congr_ae (Filter.Eventually.of_forall fun x => by ring)
  have hexp1 : ∫ x, (f x ^ 2 - 2 * a * f x + a ^ 2) ∂μ
      = (∫ x, (f x ^ 2 - 2 * a * f x) ∂μ) + ∫ x, a ^ 2 ∂μ := by
    simpa using integral_add (hf2.sub hca) (integrable_const (a ^ 2))
  have hexp2 : ∫ x, (f x ^ 2 - 2 * a * f x) ∂μ = (∫ x, f x ^ 2 ∂μ) - ∫ x, 2 * a * f x ∂μ := by
    simpa using integral_sub hf2 hca
  have hexp3 : ∫ x, 2 * a * f x ∂μ = 2 * a * ∫ x, f x ∂μ := integral_const_mul (2 * a) f
  have hexp4 : ∫ x, a ^ 2 ∂μ = a ^ 2 := by simp
  rw [hstep, hexp1, hexp2, hexp3, hexp4, integral_sq_sub_sq hf1 hf2]
  ring

/-- The variance is at most the mean square deviation from any fixed point `a`, which is in
turn bounded by `c²` when `f` stays within `c` of `a`. Only `AEStronglyMeasurable` is needed;
the integrability of `f` and `f²` is derived from the bound itself. -/
private theorem integral_sq_sub_le_of_bounded {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {f : β → ℝ} (hfm : AEStronglyMeasurable f μ) (a c : ℝ)
    (hc : ∀ x, |f x - a| ≤ c) :
    ∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ ≤ c ^ 2 := by
  have hpt : ∀ x, (f x - a) ^ 2 ≤ c ^ 2 := fun x => by
    have h := hc x; rw [abs_le] at h; nlinarith [h.1, h.2]
  have hint : Integrable (fun x => (f x - a) ^ 2) μ := by
    apply Integrable.mono' (integrable_const (c ^ 2))
      ((continuous_pow 2).comp_aestronglyMeasurable (hfm.sub aestronglyMeasurable_const))
    filter_upwards with x
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact hpt x
  have hf1 : Integrable f μ := by
    apply Integrable.mono' (integrable_const (c + |a|)) hfm
    filter_upwards with x
    have h := hc x
    rw [abs_le] at h
    rw [Real.norm_eq_abs, abs_le]
    constructor <;> nlinarith [h.1, h.2, le_abs_self a, neg_abs_le a]
  have hf2 : Integrable (fun x => f x ^ 2) μ := by
    have heq : (fun x => f x ^ 2)
        = fun x => (f x - a) ^ 2 + (2 * a * f x - a ^ 2) := by funext x; ring
    rw [heq]
    exact hint.add ((hf1.const_mul (2 * a)).sub (integrable_const (a ^ 2)))
  have hstep4 : ∫ x, (f x - a) ^ 2 ∂μ ≤ c ^ 2 := by
    calc ∫ x, (f x - a) ^ 2 ∂μ ≤ ∫ x, c ^ 2 ∂μ := integral_mono hint (integrable_const _) hpt
      _ = c ^ 2 := by simp [integral_const]
  have hexpand1 := integral_sq_sub_const_eq hf1 hf2 a
  nlinarith [hstep4, hexpand1, sq_nonneg (∫ x, f x ∂μ - a)]

/-- The independent-copy identity: the variance is half the mean squared difference between
two independent draws of `f`. -/
private theorem var_eq_half_double {β : Type*} [MeasurableSpace β] {μ : Measure β}
    [IsProbabilityMeasure μ] {W : β → ℝ} (hW1 : Integrable W μ)
    (hW2 : Integrable (fun t => W t ^ 2) μ) :
    ∫ t, (W t - ∫ s, W s ∂μ) ^ 2 ∂μ = (1 / 2) * ∫ t, ∫ s, (W t - W s) ^ 2 ∂μ ∂μ := by
  have hinner : ∀ t, ∫ s, (W t - W s) ^ 2 ∂μ
      = (∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) + (∫ u, W u ∂μ - W t) ^ 2 := by
    intro t
    have hcongr : ∫ s, (W t - W s) ^ 2 ∂μ = ∫ s, (W s - W t) ^ 2 ∂μ :=
      integral_congr_ae (Filter.Eventually.of_forall fun s => by ring)
    rw [hcongr]
    exact integral_sq_sub_const_eq hW1 hW2 (W t)
  have houter : ∫ t, (∫ s, (W t - W s) ^ 2 ∂μ) ∂μ
      = ∫ t, ((∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) + (∫ u, W u ∂μ - W t) ^ 2) ∂μ :=
    integral_congr_ae (Filter.Eventually.of_forall hinner)
  have hVarW := integral_sq_sub_sq hW1 hW2
  have hconst : Integrable (fun _ : β => ∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) μ := integrable_const _
  have hsq0 : Integrable (fun t => W t ^ 2 - 2 * (∫ u, W u ∂μ) * W t + (∫ u, W u ∂μ) ^ 2) μ :=
    (hW2.sub (hW1.const_mul (2 * ∫ u, W u ∂μ))).add (integrable_const _)
  have hsq : Integrable (fun t => (∫ u, W u ∂μ - W t) ^ 2) μ :=
    hsq0.congr (Filter.Eventually.of_forall fun t => by ring)
  have hsplit : ∫ t, ((∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) + (∫ u, W u ∂μ - W t) ^ 2) ∂μ
      = (∫ t, (∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) ∂μ) + ∫ t, (∫ u, W u ∂μ - W t) ^ 2 ∂μ := by
    simpa using integral_add hconst hsq
  have hc1 : ∫ t, (∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ) ∂μ = ∫ s, (W s - ∫ u, W u ∂μ) ^ 2 ∂μ := by
    simp
  have hc2 : ∫ t, (∫ u, W u ∂μ - W t) ^ 2 ∂μ = ∫ t, (W t - ∫ u, W u ∂μ) ^ 2 ∂μ :=
    integral_congr_ae (Filter.Eventually.of_forall fun t => by ring)
  rw [houter, hsplit, hc1, hc2]
  linarith [hVarW]

/-! ### The law of total variance and a conditional Jensen bound, for a two-factor product -/

/-- **Law of total variance**, in Fubini form: the variance of `G` on `μA.prod μB` equals the
mean, over `μB`, of the variance in the `A`-coordinate, plus the variance over `μB` of the
`A`-coordinate mean. -/
private theorem total_var_eq {A B : Type*} [MeasurableSpace A] [MeasurableSpace B]
    {μA : Measure A} {μB : Measure B} [IsProbabilityMeasure μA] [IsProbabilityMeasure μB]
    {G : A × B → ℝ} (hG1 : Integrable G (μA.prod μB))
    (hG2 : Integrable (fun p => G p ^ 2) (μA.prod μB)) :
    ∫ p, (G p - ∫ q, G q ∂(μA.prod μB)) ^ 2 ∂(μA.prod μB)
      = (∫ b, ∫ a, (G (a, b) - ∫ a', G (a', b) ∂μA) ^ 2 ∂μA ∂μB)
        + ∫ b, ((∫ a, G (a, b) ∂μA) - ∫ q, G q ∂(μA.prod μB)) ^ 2 ∂μB := by
  set g : B → ℝ := fun b => ∫ a, G (a, b) ∂μA with hg_def
  have hg1 : Integrable g μB := hG1.integral_prod_right
  have hg2 : Integrable (fun b => g b ^ 2) μB := by
    apply Integrable.mono' hG2.integral_prod_right
      ((continuous_pow 2).comp_aestronglyMeasurable hg1.aestronglyMeasurable)
    filter_upwards [hG1.prod_left_ae, hG2.prod_left_ae] with b h1 h2
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact sq_integral_le h1 h2
  have e1 : ∫ p, G p ∂(μA.prod μB) = ∫ b, g b ∂μB := integral_prod_symm G hG1
  have e2 : ∫ p, G p ^ 2 ∂(μA.prod μB) = ∫ b, ∫ a, G (a, b) ^ 2 ∂μA ∂μB :=
    integral_prod_symm (fun p => G p ^ 2) hG2
  have hterm1 : (∫ b, ∫ a, (G (a, b) - g b) ^ 2 ∂μA ∂μB)
      = (∫ b, ∫ a, G (a, b) ^ 2 ∂μA ∂μB) - ∫ b, g b ^ 2 ∂μB := by
    have hcongr : ∫ b, (∫ a, (G (a, b) - g b) ^ 2 ∂μA) ∂μB
        = ∫ b, ((∫ a, G (a, b) ^ 2 ∂μA) - g b ^ 2) ∂μB := by
      apply integral_congr_ae
      filter_upwards [hG1.prod_left_ae, hG2.prod_left_ae] with b h1 h2
      exact integral_sq_sub_sq h1 h2
    rw [hcongr]
    simpa using integral_sub hG2.integral_prod_right hg2
  have hterm2 : ∫ b, (g b - ∫ b', g b' ∂μB) ^ 2 ∂μB = (∫ b, g b ^ 2 ∂μB) - (∫ b, g b ∂μB) ^ 2 :=
    integral_sq_sub_sq hg1 hg2
  rw [integral_sq_sub_sq hG1 hG2, e1, e2, hterm1, hterm2]
  ring

/-- **Conditional Jensen.** The variance, over `μB`, of the `A`-coordinate mean of `G` is at
most the mean, over `μA`, of the variance of `G` in the `B`-coordinate. Proved through the
independent-copy identity `var_eq_half_double`, applied once to `g` and once (a.e. in `a`) to
`G (a, ·)`, then compared pointwise by Jensen on the difference of two independent copies. -/
private theorem var_mean_le_mean_var {A B : Type*} [MeasurableSpace A] [MeasurableSpace B]
    {μA : Measure A} {μB : Measure B} [IsProbabilityMeasure μA] [IsProbabilityMeasure μB]
    {G : A × B → ℝ} (hG1 : Integrable G (μA.prod μB))
    (hG2 : Integrable (fun p => G p ^ 2) (μA.prod μB)) :
    ∫ b, ((∫ a, G (a, b) ∂μA) - ∫ b', ∫ a, G (a, b') ∂μA ∂μB) ^ 2 ∂μB
      ≤ ∫ a, ∫ b, (G (a, b) - ∫ b', G (a, b') ∂μB) ^ 2 ∂μB ∂μA := by
  set g : B → ℝ := fun b => ∫ a, G (a, b) ∂μA with hg_def
  have hg1 : Integrable g μB := hG1.integral_prod_right
  have hg2 : Integrable (fun b => g b ^ 2) μB := by
    apply Integrable.mono' hG2.integral_prod_right
      ((continuous_pow 2).comp_aestronglyMeasurable hg1.aestronglyMeasurable)
    filter_upwards [hG1.prod_left_ae, hG2.prod_left_ae] with b h1 h2
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact sq_integral_le h1 h2
  have hgmem : MemLp g 2 μB := (memLp_two_iff_integrable_sq hg1.aestronglyMeasurable).mpr hg2
  have hGmem : MemLp G 2 (μA.prod μB) :=
    (memLp_two_iff_integrable_sq hG1.aestronglyMeasurable).mpr hG2
  -- transport g to B × B via the two projections
  have hgmem1 : MemLp (g ∘ Prod.fst) 2 (μB.prod μB) :=
    hgmem.comp_measurePreserving measurePreserving_fst
  have hgmem2 : MemLp (g ∘ Prod.snd) 2 (μB.prod μB) :=
    hgmem.comp_measurePreserving measurePreserving_snd
  -- transport G to A × (B × B) via id on A and the two projections on B × B
  have hid : MeasurePreserving (id : A → A) μA μA := ⟨measurable_id, Measure.map_id⟩
  have hπ1 : MeasurePreserving (Prod.map id (Prod.fst : B × B → B)) (μA.prod (μB.prod μB))
      (μA.prod μB) := hid.prod measurePreserving_fst
  have hπ2 : MeasurePreserving (Prod.map id (Prod.snd : B × B → B)) (μA.prod (μB.prod μB))
      (μA.prod μB) := hid.prod measurePreserving_snd
  have hGmem1 : MemLp (G ∘ Prod.map id (Prod.fst : B × B → B)) 2 (μA.prod (μB.prod μB)) :=
    hGmem.comp_measurePreserving hπ1
  have hGmem2 : MemLp (G ∘ Prod.map id (Prod.snd : B × B → B)) 2 (μA.prod (μB.prod μB)) :=
    hGmem.comp_measurePreserving hπ2
  set Φ : A × (B × B) → ℝ := fun p => (G (p.1, p.2.1) - G (p.1, p.2.2)) ^ 2 with hΦ_def
  have hΦeq : Φ
      = fun p => ((G ∘ Prod.map id (Prod.fst : B × B → B)) p
          - (G ∘ Prod.map id (Prod.snd : B × B → B)) p) ^ 2 := by
    funext p
    obtain ⟨a, b, b'⟩ := p
    simp [hΦ_def]
  have hΦint : Integrable Φ (μA.prod (μB.prod μB)) := by
    rw [hΦeq]; exact (hGmem1.sub hGmem2).integrable_sq
  set Φg : B × B → ℝ := fun q => (g q.1 - g q.2) ^ 2 with hΦg_def
  have hΦgeq : Φg = fun q => ((g ∘ Prod.fst) q - (g ∘ Prod.snd) q) ^ 2 := rfl
  have hΦgint : Integrable Φg (μB.prod μB) := by
    rw [hΦgeq]; exact (hgmem1.sub hgmem2).integrable_sq
  have hRHSint : Integrable (fun q : B × B => ∫ a, Φ (a, q) ∂μA) (μB.prod μB) :=
    hΦint.integral_prod_right
  -- the pointwise Jensen bound, a.e. on B × B
  have hJoint : ∀ᵐ q : B × B ∂(μB.prod μB),
      Integrable (fun a => G (a, q.1)) μA ∧ Integrable (fun a => G (a, q.1) ^ 2) μA ∧
      Integrable (fun a => G (a, q.2)) μA ∧ Integrable (fun a => G (a, q.2) ^ 2) μA := by
    set mpfst := measurePreserving_fst (μ := μB) (ν := μB) with hmpfst_def
    set mpsnd := measurePreserving_snd (μ := μB) (ν := μB) with hmpsnd_def
    have hae1 := mpfst.quasiMeasurePreserving.tendsto_ae.eventually
      (p := fun b => Integrable (fun a => G (a, b)) μA) hG1.prod_left_ae
    have hae2 := mpfst.quasiMeasurePreserving.tendsto_ae.eventually
      (p := fun b => Integrable (fun a => G (a, b) ^ 2) μA) hG2.prod_left_ae
    have hae3 := mpsnd.quasiMeasurePreserving.tendsto_ae.eventually
      (p := fun b => Integrable (fun a => G (a, b)) μA) hG1.prod_left_ae
    have hae4 := mpsnd.quasiMeasurePreserving.tendsto_ae.eventually
      (p := fun b => Integrable (fun a => G (a, b) ^ 2) μA) hG2.prod_left_ae
    filter_upwards [hae1, hae2, hae3, hae4] with q h1 h2 h3 h4 using ⟨h1, h2, h3, h4⟩
  have hptJensen : ∀ᵐ q : B × B ∂(μB.prod μB),
      (g q.1 - g q.2) ^ 2 ≤ ∫ a, (G (a, q.1) - G (a, q.2)) ^ 2 ∂μA := by
    filter_upwards [hJoint] with q hjq
    obtain ⟨hb1, hb2, hb1', hb2'⟩ := hjq
    have hcross : Integrable (fun a => G (a, q.1) * G (a, q.2)) μA :=
      (((memLp_two_iff_integrable_sq hb1'.aestronglyMeasurable).mpr hb2').mul
        ((memLp_two_iff_integrable_sq hb1.aestronglyMeasurable).mpr hb2)).integrable le_rfl
    have hdiff1 : Integrable (fun a => G (a, q.1) - G (a, q.2)) μA := hb1.sub hb1'
    have hdiff2 : Integrable (fun a => (G (a, q.1) - G (a, q.2)) ^ 2) μA := by
      have heq : (fun a => (G (a, q.1) - G (a, q.2)) ^ 2)
          = fun a => G (a, q.1) ^ 2 - 2 * (G (a, q.1) * G (a, q.2)) + G (a, q.2) ^ 2 := by
        funext a; ring
      rw [heq]; exact (hb2.sub (hcross.const_mul 2)).add hb2'
    have hgeq : g q.1 - g q.2 = ∫ a, (G (a, q.1) - G (a, q.2)) ∂μA := (integral_sub hb1 hb1').symm
    rw [hgeq]
    exact sq_integral_le hdiff1 hdiff2
  have hcompare : ∫ q : B × B, Φg q ∂(μB.prod μB)
      ≤ ∫ q : B × B, (∫ a, Φ (a, q) ∂μA) ∂(μB.prod μB) := by
    apply integral_mono_ae hΦgint hRHSint
    filter_upwards [hptJensen] with q hq
    exact hq
  -- left side: Φg integrates to 2 * Var_B(g)
  have hleft : ∫ q : B × B, Φg q ∂(μB.prod μB) = 2 * ∫ b, (g b - ∫ b', g b' ∂μB) ^ 2 ∂μB := by
    have e1 : ∫ q : B × B, Φg q ∂(μB.prod μB) = ∫ b, ∫ b', Φg (b, b') ∂μB ∂μB :=
      integral_prod Φg hΦgint
    rw [e1, var_eq_half_double hg1 hg2]
    ring
  -- right side: relate to the target's RHS via a Fubini swap of Φ and var_eq_half_double at a
  have hright : ∫ q : B × B, (∫ a, Φ (a, q) ∂μA) ∂(μB.prod μB)
      = 2 * ∫ a, ∫ b, (G (a, b) - ∫ b', G (a, b') ∂μB) ^ 2 ∂μB ∂μA := by
    have e2 : ∫ q : B × B, (∫ a, Φ (a, q) ∂μA) ∂(μB.prod μB)
        = ∫ p : A × (B × B), Φ p ∂(μA.prod (μB.prod μB)) := (integral_prod_symm Φ hΦint).symm
    have e3 : ∫ p : A × (B × B), Φ p ∂(μA.prod (μB.prod μB))
        = ∫ a, ∫ q, Φ (a, q) ∂(μB.prod μB) ∂μA := integral_prod Φ hΦint
    have hae : ∀ᵐ a ∂μA,
        Integrable (fun b => G (a, b)) μB ∧ Integrable (fun b => G (a, b) ^ 2) μB := by
      filter_upwards [hG1.prod_right_ae, hG2.prod_right_ae] with a h1 h2 using ⟨h1, h2⟩
    have e4 : ∀ᵐ a ∂μA, ∫ q, Φ (a, q) ∂(μB.prod μB)
        = 2 * ∫ b, (G (a, b) - ∫ b', G (a, b') ∂μB) ^ 2 ∂μB := by
      filter_upwards [hae] with a ha
      obtain ⟨hWa1, hWa2⟩ := ha
      have hΦa : Integrable (fun q : B × B => Φ (a, q)) (μB.prod μB) := by
        have heq : (fun q : B × B => Φ (a, q)) = fun q => (G (a, q.1) - G (a, q.2)) ^ 2 := rfl
        rw [heq]
        have hcross : Integrable (fun q : B × B => G (a, q.1) * G (a, q.2)) (μB.prod μB) := by
          have h1' : MemLp (fun q : B × B => G (a, q.1)) 2 (μB.prod μB) :=
            ((memLp_two_iff_integrable_sq hWa1.aestronglyMeasurable).mpr
              hWa2).comp_measurePreserving measurePreserving_fst
          have h2' : MemLp (fun q : B × B => G (a, q.2)) 2 (μB.prod μB) :=
            ((memLp_two_iff_integrable_sq hWa1.aestronglyMeasurable).mpr
              hWa2).comp_measurePreserving measurePreserving_snd
          exact (h2'.mul h1').integrable le_rfl
        have hd1sq : Integrable (fun q : B × B => G (a, q.1) ^ 2) (μB.prod μB) :=
          measurePreserving_fst.integrable_comp_of_integrable hWa2
        have hd2sq : Integrable (fun q : B × B => G (a, q.2) ^ 2) (μB.prod μB) :=
          measurePreserving_snd.integrable_comp_of_integrable hWa2
        have hexp : (fun q : B × B => (G (a, q.1) - G (a, q.2)) ^ 2)
            = fun q => G (a, q.1) ^ 2 - 2 * (G (a, q.1) * G (a, q.2)) + G (a, q.2) ^ 2 := by
          funext q; ring
        rw [hexp]
        exact (hd1sq.sub (hcross.const_mul 2)).add hd2sq
      have e5 : ∫ q : B × B, Φ (a, q) ∂(μB.prod μB) = ∫ b, ∫ b', Φ (a, (b, b')) ∂μB ∂μB :=
        integral_prod (fun q => Φ (a, q)) hΦa
      rw [e5, var_eq_half_double hWa1 hWa2]
      ring
    rw [e2, e3, integral_congr_ae e4, integral_const_mul]
  rw [hleft, hright] at hcompare
  linarith [hcompare]

/-- Change of variables along a (not necessarily injective) measure preserving map, without
the `MeasurableEmbedding` side condition of `MeasurePreserving.integral_comp`. -/
private theorem integral_comp_of_measurePreserving {A B : Type*} [MeasurableSpace A]
    [MeasurableSpace B] {μ : Measure A} {ν : Measure B} {T : A → B} (hT : MeasurePreserving T μ ν)
    (g : B → ℝ) (hg : AEStronglyMeasurable g ν) :
    ∫ x, g (T x) ∂μ = ∫ y, g y ∂ν := by
  rw [← hT.map_eq] at hg ⊢
  exact (integral_map hT.measurable.aemeasurable hg).symm

/-- For any coordinate `i` of a `Fin (n + 1)`-indexed product, the mean of the variance in that
coordinate is an integrable function of the point, and its integral equals the two-factor form
obtained by splitting off coordinate `i` with `measurePreserving_piFinSuccAbove`. Reused for the
split-off coordinate `0` and, in the recursive step, for each coordinate `j.succ`. -/
private theorem integral_var_term_eq {n : ℕ} {α : Fin (n + 1) → Type*}
    [∀ k, MeasurableSpace (α k)] (ν : ∀ k, Measure (α k)) [∀ k, IsProbabilityMeasure (ν k)]
    {f : (∀ k, α k) → ℝ} (hf1 : Integrable f (Measure.pi ν))
    (hf2 : Integrable (fun x => f x ^ 2) (Measure.pi ν)) (i : Fin (n + 1)) :
    Integrable (fun x => ∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) (Measure.pi ν) ∧
    ∫ x, (∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν)
      = ∫ r, ∫ y, (f (Fin.insertNth i y r) - ∫ y', f (Fin.insertNth i y' r) ∂(ν i)) ^ 2 ∂(ν i)
          ∂(Measure.pi fun k : Fin n => ν (i.succAbove k)) := by
  set e : (∀ k, α k) ≃ᵐ α i × ∀ k : Fin n, α (i.succAbove k) := MeasurableEquiv.piFinSuccAbove α i
    with he_def
  have hmp : MeasurePreserving e (Measure.pi ν)
      ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) := measurePreserving_piFinSuccAbove ν i
  set F : α i × (∀ k : Fin n, α (i.succAbove k)) → ℝ := f ∘ e.symm with hF_def
  have hFint1 : Integrable F ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) :=
    hmp.symm.integrable_comp_of_integrable hf1
  have hFint2 : Integrable (fun p => F p ^ 2)
      ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) :=
    hmp.symm.integrable_comp_of_integrable hf2
  have hg1 : Integrable (fun r => ∫ y, F (y, r) ∂(ν i))
      (Measure.pi fun k => ν (i.succAbove k)) := hFint1.integral_prod_right
  have hg2 : Integrable (fun r => (∫ y, F (y, r) ∂(ν i)) ^ 2)
      (Measure.pi fun k => ν (i.succAbove k)) := by
    apply Integrable.mono' hFint2.integral_prod_right
      ((continuous_pow 2).comp_aestronglyMeasurable hg1.aestronglyMeasurable)
    filter_upwards [hFint1.prod_left_ae, hFint2.prod_left_ae] with r h1 h2
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact sq_integral_le h1 h2
  set H : (∀ k : Fin n, α (i.succAbove k)) → ℝ :=
    fun r => ∫ y, (F (y, r) - ∫ y', F (y', r) ∂(ν i)) ^ 2 ∂(ν i) with hH_def
  have hHeq : ∀ᵐ r ∂(Measure.pi fun k => ν (i.succAbove k)),
      H r = (∫ y, F (y, r) ^ 2 ∂(ν i)) - (∫ y, F (y, r) ∂(ν i)) ^ 2 := by
    filter_upwards [hFint1.prod_left_ae, hFint2.prod_left_ae] with r h1 h2
    exact integral_sq_sub_sq h1 h2
  have hHint : Integrable H (Measure.pi fun k => ν (i.succAbove k)) :=
    (hFint2.integral_prod_right.sub hg2).congr (hHeq.mono fun r hr => hr.symm)
  have hK : Integrable ((H ∘ Prod.snd) ∘ e) (Measure.pi ν) :=
    hmp.integrable_comp_of_integrable
      ((measurePreserving_snd (μ := ν i)).integrable_comp_of_integrable hHint)
  have hupdate : ∀ (x : ∀ k, α k) (t : α i), Function.update x i t = e.symm (t, (e x).2) := by
    intro x t
    change Function.update x i t = Fin.insertNth i t (Fin.removeNth i x)
    have hx : x = Fin.insertNth i (x i) (Fin.removeNth i x) :=
      (Fin.insertNth_self_removeNth i x).symm
    conv_lhs => rw [hx]
    exact Fin.update_insertNth i (x i) t (Fin.removeNth i x)
  have hFeq : ∀ (x : ∀ k, α k) (t : α i), F (t, (e x).2) = f (Function.update x i t) := by
    intro x t; rw [hupdate x t]; rfl
  have hpt : ∀ x : ∀ k, α k, ((H ∘ Prod.snd) ∘ e) x = ∫ t, (f (Function.update x i t)
      - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i) := by
    intro x
    change (∫ y, (F (y, (e x).2) - ∫ y', F (y', (e x).2) ∂(ν i)) ^ 2 ∂(ν i))
        = ∫ t, (f (Function.update x i t) - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)
    simp_rw [hFeq x]
  refine ⟨hK.congr (Filter.Eventually.of_forall hpt), ?_⟩
  have hstep1 : ∫ x, (∫ t, (f (Function.update x i t)
      - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν)
      = ∫ x, H ((e x).2) ∂(Measure.pi ν) := by
    apply integral_congr_ae
    filter_upwards with x
    exact (hpt x).symm
  have hsnde : MeasurePreserving (Prod.snd ∘ e) (Measure.pi ν)
      (Measure.pi fun k => ν (i.succAbove k)) := (measurePreserving_snd (μ := ν i)).comp hmp
  have hstep2 : ∫ x, H ((e x).2) ∂(Measure.pi ν)
      = ∫ r, H r ∂(Measure.pi fun k => ν (i.succAbove k)) := by
    simpa using integral_comp_of_measurePreserving hsnde H hHint.aestronglyMeasurable
  have hstep3 : ∫ r, H r ∂(Measure.pi fun k => ν (i.succAbove k))
      = ∫ r, ∫ y, (f (Fin.insertNth i y r) - ∫ y', f (Fin.insertNth i y' r) ∂(ν i)) ^ 2 ∂(ν i)
          ∂(Measure.pi fun k => ν (i.succAbove k)) := by
    apply integral_congr_ae
    filter_upwards with r
    rfl
  rw [hstep1, hstep2, hstep3]

/-! ### The induction on `Fin D` -/

/-- The tensorization inequality for a `Fin D`-indexed product, by induction on `D`. The step
splits off coordinate `0` with `measurePreserving_piFinSuccAbove`, applies `total_var_eq` to
the resulting two-factor product, and bounds the recursion term on the remaining `D`
coordinates with the induction hypothesis and `var_mean_le_mean_var`. -/
private theorem variance_pi_le_sum_fin :
    ∀ (D : ℕ) {α : Fin D → Type*} [∀ i, MeasurableSpace (α i)] (ν : ∀ i, Measure (α i))
      [∀ i, IsProbabilityMeasure (ν i)] (f : (∀ i, α i) → ℝ), MemLp f 2 (Measure.pi ν) →
      ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
        ≤ ∑ i, ∫ x, (∫ t, (f (Function.update x i t)
            - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν) := by
  intro D
  induction D with
  | zero =>
    intro α _ ν _ f _
    have hmean : ∀ x, ∫ y, f y ∂(Measure.pi ν) = f x := fun x => by
      have heq : (fun y => f y) = fun _ => f x := funext fun y => by
        congr 1; exact Subsingleton.elim y x
      rw [heq]; simp
    have hz : ∀ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 = 0 := fun x => by rw [hmean x]; ring
    simp [hz]
  | succ D ih =>
    intro α _ ν _ f hf
    classical
    set μ0 : Measure (α 0) := ν 0 with hμ0_def
    set μR : Measure (∀ j : Fin D, α j.succ) := Measure.pi (fun j => ν j.succ) with hμR_def
    set e : (∀ i, α i) ≃ᵐ α 0 × ∀ j : Fin D, α j.succ := MeasurableEquiv.piFinSuccAbove α 0
      with he_def
    have hmp : MeasurePreserving e (Measure.pi ν) (μ0.prod μR) :=
      measurePreserving_piFinSuccAbove ν 0
    set F : α 0 × (∀ j : Fin D, α j.succ) → ℝ := f ∘ e.symm with hF_def
    have hFmem : MemLp F 2 (μ0.prod μR) := hf.comp_measurePreserving hmp.symm
    have hF1 : Integrable F (μ0.prod μR) := hFmem.integrable (by norm_num)
    have hF2 : Integrable (fun p => F p ^ 2) (μ0.prod μR) := hFmem.integrable_sq
    have hFe : ∀ x, F (e x) = f x := fun x => by simp [hF_def]
    have hcenter : ∫ y, f y ∂(Measure.pi ν) = ∫ q, F q ∂(μ0.prod μR) := by
      rw [← hmp.integral_comp' F]
      exact integral_congr_ae (Filter.Eventually.of_forall fun x => (hFe x).symm)
    have hLHSeq : ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
        = ∫ p, (F p - ∫ q, F q ∂(μ0.prod μR)) ^ 2 ∂(μ0.prod μR) := by
      rw [← hmp.integral_comp' (fun p => (F p - ∫ q, F q ∂(μ0.prod μR)) ^ 2)]
      apply integral_congr_ae
      filter_upwards with x
      rw [hFe x, hcenter]
    have htv := total_var_eq (μA := μ0) (μB := μR) (G := F) hF1 hF2
    have hf1 : Integrable f (Measure.pi ν) := hf.integrable (by norm_num)
    have hf2 : Integrable (fun x => f x ^ 2) (Measure.pi ν) := hf.integrable_sq
    have hi0 := integral_var_term_eq ν hf1 hf2 (0 : Fin (D + 1))
    have hvmm := var_mean_le_mean_var (μA := μ0) (μB := μR) (G := F) hF1 hF2
    have haeF : ∀ᵐ y ∂μ0,
        Integrable (fun r => F (y, r)) μR ∧ Integrable (fun r => F (y, r) ^ 2) μR := by
      filter_upwards [hF1.prod_right_ae, hF2.prod_right_ae] with y h1 h2 using ⟨h1, h2⟩
    have hstep2 : ∀ᵐ y ∂μ0, ∫ r, (F (y, r) - ∫ r', F (y, r') ∂μR) ^ 2 ∂μR
        ≤ ∑ j, ∫ r, (∫ t, (F (y, Function.update r j t)
            - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR := by
      filter_upwards [haeF] with y hy
      exact ih (fun j => ν j.succ) (fun r => F (y, r))
        ((memLp_two_iff_integrable_sq hy.1.aestronglyMeasurable).mpr hy.2)
    have hkey : ∀ (y : α 0) (r : ∀ k : Fin D, α k.succ) (j : Fin D) (t : α j.succ),
        Function.update (e.symm (y, r)) j.succ t = e.symm (y, Function.update r j t) :=
      fun y r j t => (Fin.insertNth_update 0 y j t r).symm
    have hXieq : ∀ (j : Fin D) (y : α 0) (r : ∀ k : Fin D, α k.succ),
        (∫ t, (F (y, Function.update r j t)
            - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
          = ∫ t, (f (Function.update (e.symm (y, r)) j.succ t)
              - ∫ s, f (Function.update (e.symm (y, r)) j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ) := by
      intro j y r
      have hFeq : ∀ t, F (y, Function.update r j t)
          = f (Function.update (e.symm (y, r)) j.succ t) := fun t => by rw [hkey y r j t]; rfl
      simp_rw [hFeq]
    have hLj : ∀ j : Fin D, Integrable (fun x => ∫ t, (f (Function.update x j.succ t)
        - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) (Measure.pi ν) :=
      fun j => (integral_var_term_eq ν hf1 hf2 j.succ).1
    have hpull : ∀ j : Fin D, Integrable (fun p : α 0 × (∀ k : Fin D, α k.succ) =>
        ∫ t, (F (p.1, Function.update p.2 j t)
            - ∫ s, F (p.1, Function.update p.2 j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) (μ0.prod μR) := by
      intro j
      have h0 := hmp.symm.integrable_comp_of_integrable (hLj j)
      refine h0.congr (Filter.Eventually.of_forall fun p => ?_)
      exact (hXieq j p.1 p.2).symm
    have hpj : ∀ j : Fin D, Integrable (fun y => ∫ r, (∫ t, (F (y, Function.update r j t)
        - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR) μ0 :=
      fun j => (hpull j).integral_prod_left
    have hm1 : Integrable (fun y => ∫ r, F (y, r) ∂μR) μ0 := hF1.integral_prod_left
    have hm2 : Integrable (fun y => (∫ r, F (y, r) ∂μR) ^ 2) μ0 := by
      apply Integrable.mono' hF2.integral_prod_left
        ((continuous_pow 2).comp_aestronglyMeasurable hm1.aestronglyMeasurable)
      filter_upwards [hF1.prod_right_ae, hF2.prod_right_ae] with y h1 h2
      rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
      exact sq_integral_le h1 h2
    have hΨeq : ∀ᵐ y ∂μ0, (∫ r, (F (y, r) - ∫ r', F (y, r') ∂μR) ^ 2 ∂μR)
        = (∫ r, F (y, r) ^ 2 ∂μR) - (∫ r, F (y, r) ∂μR) ^ 2 := by
      filter_upwards [hF1.prod_right_ae, hF2.prod_right_ae] with y h1 h2
      exact integral_sq_sub_sq h1 h2
    have hΨint : Integrable (fun y => ∫ r, (F (y, r) - ∫ r', F (y, r') ∂μR) ^ 2 ∂μR) μ0 :=
      (hF2.integral_prod_left.sub hm2).congr (hΨeq.mono fun y hy => hy.symm)
    have hstep3 : ∫ y, ∫ r, (F (y, r) - ∫ r', F (y, r') ∂μR) ^ 2 ∂μR ∂μ0
        ≤ ∫ y, ∑ j, ∫ r, (∫ t, (F (y, Function.update r j t)
            - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0 := by
      apply integral_mono_ae hΨint ?_ hstep2
      exact integrable_finsetSum Finset.univ fun j _ => hpj j
    have hswap : ∫ y, ∑ j, ∫ r, (∫ t, (F (y, Function.update r j t)
        - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0
        = ∑ j, ∫ y, ∫ r, (∫ t, (F (y, Function.update r j t)
            - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0 :=
      integral_finsetSum Finset.univ fun j _ => hpj j
    have hstep4 : ∀ j : Fin D, ∫ y, ∫ r, (∫ t, (F (y, Function.update r j t)
        - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0
        = ∫ x, (∫ t, (f (Function.update x j.succ t)
            - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
            ∂(Measure.pi ν) := by
      intro j
      have e1 : ∫ y, ∫ r, (∫ t, (F (y, Function.update r j t)
          - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0
          = ∫ p : α 0 × (∀ k : Fin D, α k.succ), (∫ t, (F (p.1, Function.update p.2 j t)
              - ∫ s, F (p.1, Function.update p.2 j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂(μ0.prod μR) :=
        (integral_prod _ (hpull j)).symm
      have e2 : ∫ x, (∫ t, (f (Function.update x j.succ t)
          - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂(Measure.pi ν)
          = ∫ p : α 0 × (∀ k : Fin D, α k.succ), (∫ t, (f (Function.update (e.symm p) j.succ t)
              - ∫ s, f (Function.update (e.symm p) j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
              ∂(μ0.prod μR) := by
        rw [← hmp.integral_comp' (fun p => ∫ t, (f (Function.update (e.symm p) j.succ t)
            - ∫ s, f (Function.update (e.symm p) j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))]
        apply integral_congr_ae
        filter_upwards with x
        rw [MeasurableEquiv.symm_apply_apply]
      rw [e1, e2]
      apply integral_congr_ae
      filter_upwards with p
      exact hXieq j p.1 p.2
    have hgcenter : ∫ q, F q ∂(μ0.prod μR) = ∫ r', ∫ y, F (y, r') ∂μ0 ∂μR :=
      integral_prod_symm F hF1
    have hi0term : ∫ x, (∫ t, (f (Function.update x 0 t)
        - ∫ s, f (Function.update x 0 s) ∂(ν 0)) ^ 2 ∂(ν 0)) ∂(Measure.pi ν)
        = ∫ r, ∫ y, (F (y, r) - ∫ y', F (y', r) ∂μ0) ^ 2 ∂μ0 ∂μR := hi0.2
    have hbound : ∫ r, ((∫ y, F (y, r) ∂μ0) - ∫ q, F q ∂(μ0.prod μR)) ^ 2 ∂μR
        ≤ ∑ j : Fin D, ∫ x, (∫ t, (f (Function.update x j.succ t)
            - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
            ∂(Measure.pi ν) := by
      rw [hgcenter]
      calc ∫ r, ((∫ y, F (y, r) ∂μ0) - ∫ r', ∫ y, F (y, r') ∂μ0 ∂μR) ^ 2 ∂μR
          ≤ ∫ y, ∫ r, (F (y, r) - ∫ r', F (y, r') ∂μR) ^ 2 ∂μR ∂μ0 := hvmm
        _ ≤ ∫ y, ∑ j, ∫ r, (∫ t, (F (y, Function.update r j t)
              - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0 :=
          hstep3
        _ = ∑ j, ∫ y, ∫ r, (∫ t, (F (y, Function.update r j t)
              - ∫ s, F (y, Function.update r j s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂μR ∂μ0 :=
          hswap
        _ = ∑ j : Fin D, ∫ x, (∫ t, (f (Function.update x j.succ t)
              - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
              ∂(Measure.pi ν) := Finset.sum_congr rfl fun j _ => hstep4 j
    have hsplit : ∑ i : Fin (D + 1), ∫ x, (∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν)
        = (∫ x, (∫ t, (f (Function.update x 0 t)
            - ∫ s, f (Function.update x 0 s) ∂(ν 0)) ^ 2 ∂(ν 0)) ∂(Measure.pi ν))
          + ∑ j : Fin D, ∫ x, (∫ t, (f (Function.update x j.succ t)
              - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ)) ∂(Measure.pi ν) :=
      Fin.sum_univ_succ (fun i => ∫ x, (∫ t, (f (Function.update x i t)
          - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν))
    rw [hsplit]
    calc ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
        = (∫ r, ∫ y, (F (y, r) - ∫ y', F (y', r) ∂μ0) ^ 2 ∂μ0 ∂μR)
          + ∫ r, ((∫ y, F (y, r) ∂μ0) - ∫ q, F q ∂(μ0.prod μR)) ^ 2 ∂μR := by rw [hLHSeq, htv]
      _ = (∫ x, (∫ t, (f (Function.update x 0 t)
              - ∫ s, f (Function.update x 0 s) ∂(ν 0)) ^ 2 ∂(ν 0)) ∂(Measure.pi ν))
          + ∫ r, ((∫ y, F (y, r) ∂μ0) - ∫ q, F q ∂(μ0.prod μR)) ^ 2 ∂μR := by rw [hi0term]
      _ ≤ (∫ x, (∫ t, (f (Function.update x 0 t)
              - ∫ s, f (Function.update x 0 s) ∂(ν 0)) ^ 2 ∂(ν 0)) ∂(Measure.pi ν))
          + ∑ j : Fin D, ∫ x, (∫ t, (f (Function.update x j.succ t)
              - ∫ s, f (Function.update x j.succ s) ∂(ν j.succ)) ^ 2 ∂(ν j.succ))
              ∂(Measure.pi ν) := by linarith [hbound]

/-! ### Transport to a general `Fintype` -/

/-- Reindexing along an equivalence commutes with `Function.update`, in the sense that updating
the transported point at a transported index matches transporting the update. -/
private theorem piCongrLeft_update {ι β : Type*} [DecidableEq ι] [DecidableEq β] {α : ι → Type*}
    (e : β ≃ ι) (z : ∀ b, α (e b)) (b₀ : β) (t : α (e b₀)) :
    Equiv.piCongrLeft α e (Function.update z b₀ t)
      = Function.update (Equiv.piCongrLeft α e z) (e b₀) t := by
  funext i
  obtain ⟨b, rfl⟩ := e.surjective i
  rw [Equiv.piCongrLeft_apply_apply]
  rcases eq_or_ne b b₀ with rfl | hb
  · simp
  · rw [Function.update_of_ne hb, Function.update_of_ne (e.injective.ne hb),
      Equiv.piCongrLeft_apply_apply]

/-- **Tensorization of the variance (Efron-Stein) on a finite product measure.** For a
square integrable `f` on `Measure.pi ν`, the variance is at most the sum over the
coordinates of the mean of the variance in that coordinate. Proved for `Fin D` in
`variance_pi_le_sum_fin`, then transported by the reindexing equivalence
`e : Fin (Fintype.card ι) ≃ ι` and `MeasureTheory.measurePreserving_piCongrLeft`, in the style
of `MeasureTheory.Integrable.fintype_prod_dep`. -/
theorem variance_pi_le_sum {ι : Type*} [Fintype ι] [DecidableEq ι] {α : ι → Type*}
    [∀ i, MeasurableSpace (α i)] (ν : ∀ i, Measure (α i)) [∀ i, IsProbabilityMeasure (ν i)]
    (f : (∀ i, α i) → ℝ) (hf : MemLp f 2 (Measure.pi ν)) :
    ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
      ≤ ∑ i, ∫ x, (∫ t, (f (Function.update x i t)
          - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν) := by
  set e : Fin (Fintype.card ι) ≃ ι := (Fintype.equivFin ι).symm with he_def
  set E : (∀ b : Fin (Fintype.card ι), α (e b)) ≃ᵐ ∀ i, α i := MeasurableEquiv.piCongrLeft α e
    with hE_def
  set μ' : ∀ b : Fin (Fintype.card ι), Measure (α (e b)) := fun b => ν (e b) with hμ'_def
  have hmp : MeasurePreserving E (Measure.pi μ') (Measure.pi ν) := measurePreserving_piCongrLeft ν e
  set F : (∀ b : Fin (Fintype.card ι), α (e b)) → ℝ := f ∘ E with hF_def
  have hFmem : MemLp F 2 (Measure.pi μ') := hf.comp_measurePreserving hmp
  have main := variance_pi_le_sum_fin (Fintype.card ι) μ' F hFmem
  have hcenter : ∫ y, f y ∂(Measure.pi ν) = ∫ z, F z ∂(Measure.pi μ') := (hmp.integral_comp' f).symm
  have hLHSeq : ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
      = ∫ z, (F z - ∫ z', F z' ∂(Measure.pi μ')) ^ 2 ∂(Measure.pi μ') := by
    rw [← hmp.integral_comp' (fun x => (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2)]
    apply integral_congr_ae
    filter_upwards with z
    rw [hcenter]; rfl
  have hupdate : ∀ (z : ∀ b, α (e b)) (b₀ : Fin (Fintype.card ι)) (t : α (e b₀)),
      E (Function.update z b₀ t) = Function.update (E z) (e b₀) t :=
    fun z b₀ t => piCongrLeft_update e z b₀ t
  have hterm : ∀ b : Fin (Fintype.card ι), ∫ z, (∫ t, (F (Function.update z b t)
      - ∫ s, F (Function.update z b s) ∂(μ' b)) ^ 2 ∂(μ' b)) ∂(Measure.pi μ')
      = ∫ x, (∫ t, (f (Function.update x (e b) t)
          - ∫ s, f (Function.update x (e b) s) ∂(ν (e b))) ^ 2 ∂(ν (e b)))
          ∂(Measure.pi ν) := by
    intro b
    rw [← hmp.integral_comp' (fun x => ∫ t, (f (Function.update x (e b) t)
        - ∫ s, f (Function.update x (e b) s) ∂(ν (e b))) ^ 2 ∂(ν (e b)))]
    apply integral_congr_ae
    filter_upwards with z
    have hFeq : ∀ t' : α (e b), F (Function.update z b t')
        = f (Function.update (E z) (e b) t') := fun t' => by
      change f (E (Function.update z b t')) = f (Function.update (E z) (e b) t')
      rw [hupdate z b t']
    simp_rw [hFeq]
    rfl
  have hsum : ∑ i, ∫ x, (∫ t, (f (Function.update x i t)
      - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν)
      = ∑ b : Fin (Fintype.card ι), ∫ z, (∫ t, (F (Function.update z b t)
          - ∫ s, F (Function.update z b s) ∂(μ' b)) ^ 2 ∂(μ' b)) ∂(Measure.pi μ') := by
    rw [← Equiv.sum_comp e (fun i => ∫ x, (∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν))]
    exact Finset.sum_congr rfl fun b _ => (hterm b).symm
  rw [hsum, hLHSeq]
  exact main

/-! ### Bounded differences -/

/-- The `Fin (n + 1)` form of one summand of the bounded-difference bound: if replacing
coordinate `i` moves `f` by at most `c`, that coordinate's variance term is at most `c ^ 2`.
Reuses the same coordinate split as `integral_var_term_eq`. -/
private theorem integral_var_term_le {n : ℕ} {α : Fin (n + 1) → Type*}
    [∀ k, MeasurableSpace (α k)] (ν : ∀ k, Measure (α k)) [∀ k, IsProbabilityMeasure (ν k)]
    {f : (∀ k, α k) → ℝ} (hf1 : Integrable f (Measure.pi ν))
    (hf2 : Integrable (fun x => f x ^ 2) (Measure.pi ν)) (i : Fin (n + 1)) (c : ℝ)
    (hc : ∀ x t, |f (Function.update x i t) - f x| ≤ c) :
    ∫ x, (∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν) ≤ c ^ 2 := by
  set e : (∀ k, α k) ≃ᵐ α i × ∀ k : Fin n, α (i.succAbove k) := MeasurableEquiv.piFinSuccAbove α i
    with he_def
  have hmp : MeasurePreserving e (Measure.pi ν)
      ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) := measurePreserving_piFinSuccAbove ν i
  set F : α i × (∀ k : Fin n, α (i.succAbove k)) → ℝ := f ∘ e.symm with hF_def
  have hFint1 : Integrable F ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) :=
    hmp.symm.integrable_comp_of_integrable hf1
  have hFint2 : Integrable (fun p => F p ^ 2)
      ((ν i).prod (Measure.pi fun k => ν (i.succAbove k))) :=
    hmp.symm.integrable_comp_of_integrable hf2
  obtain ⟨y₀⟩ := nonempty_of_isProbabilityMeasure (ν i)
  have hupdate : ∀ (x : ∀ k, α k) (t : α i), Function.update x i t = e.symm (t, (e x).2) := by
    intro x t
    change Function.update x i t = Fin.insertNth i t (Fin.removeNth i x)
    have hx : x = Fin.insertNth i (x i) (Fin.removeNth i x) :=
      (Fin.insertNth_self_removeNth i x).symm
    conv_lhs => rw [hx]
    exact Fin.update_insertNth i (x i) t (Fin.removeNth i x)
  have hcF : ∀ (r : ∀ k : Fin n, α (i.succAbove k)) (t : α i), |F (t, r) - F (y₀, r)| ≤ c := by
    intro r t
    have he : e (e.symm (y₀, r)) = (y₀, r) := MeasurableEquiv.apply_symm_apply e (y₀, r)
    have h1 : F (t, r) = f (Function.update (e.symm (y₀, r)) i t) := by
      rw [hupdate (e.symm (y₀, r)) t, he]; rfl
    have h2 : F (y₀, r) = f (e.symm (y₀, r)) := rfl
    rw [h1, h2]
    exact hc (e.symm (y₀, r)) t
  have hbound : ∀ᵐ r ∂(Measure.pi fun k => ν (i.succAbove k)),
      ∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i) ≤ c ^ 2 := by
    filter_upwards [hFint1.prod_left_ae] with r hr
    exact integral_sq_sub_le_of_bounded hr.aestronglyMeasurable (F (y₀, r)) c (hcF r)
  have hHmeas : AEStronglyMeasurable (fun r => ∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i))
      (Measure.pi fun k => ν (i.succAbove k)) := by
    have h1 : AEStronglyMeasurable (fun r => ∫ t, F (t, r) ^ 2 ∂(ν i))
        (Measure.pi fun k => ν (i.succAbove k)) := hFint2.integral_prod_right.aestronglyMeasurable
    have h2 : AEStronglyMeasurable (fun r => ∫ t, F (t, r) ∂(ν i))
        (Measure.pi fun k => ν (i.succAbove k)) := hFint1.integral_prod_right.aestronglyMeasurable
    have hHeq : ∀ᵐ r ∂(Measure.pi fun k => ν (i.succAbove k)),
        (∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i))
          = (∫ t, F (t, r) ^ 2 ∂(ν i)) - (∫ t, F (t, r) ∂(ν i)) ^ 2 := by
      filter_upwards [hFint1.prod_left_ae, hFint2.prod_left_ae] with r hr1 hr2
      exact integral_sq_sub_sq hr1 hr2
    exact (h1.sub ((continuous_pow 2).comp_aestronglyMeasurable h2)).congr
      (hHeq.mono fun r hr => hr.symm)
  have hHint : Integrable (fun r => ∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i))
      (Measure.pi fun k => ν (i.succAbove k)) := by
    apply Integrable.mono' (integrable_const (c ^ 2)) hHmeas
    filter_upwards [hbound] with r hr
    rwa [Real.norm_eq_abs, abs_of_nonneg (integral_nonneg (fun _ => sq_nonneg _))]
  calc ∫ x, (∫ t, (f (Function.update x i t)
        - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν)
      = ∫ x, (∫ t, (F (t, (e x).2) - ∫ s, F (s, (e x).2) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν) := by
        apply integral_congr_ae
        filter_upwards with x
        simp_rw [hupdate x]
        rfl
    _ = ∫ r, ∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i)
          ∂(Measure.pi fun k => ν (i.succAbove k)) := by
        simpa using integral_comp_of_measurePreserving
          ((measurePreserving_snd (μ := ν i)).comp hmp)
          (fun r => ∫ t, (F (t, r) - ∫ s, F (s, r) ∂(ν i)) ^ 2 ∂(ν i)) hHmeas
    _ ≤ ∫ _r, c ^ 2 ∂(Measure.pi fun k => ν (i.succAbove k)) :=
        integral_mono_ae hHint (integrable_const _) hbound
    _ = c ^ 2 := by simp

/-- **Bounded differences, `Fin D` form.** If replacing coordinate `i` moves `f` by at most
`c i`, the variance is at most `∑ i, (c i) ^ 2`. This is the constant-`1` form of the classical
bounded-difference inequality; the sharper constant `1 / 2` recenters each coordinate's
variance at an independent copy rather than at `f x`, which is not used here. Stated for
`Fin D`: the general `Fintype` case needs the same reindexing transport as
`variance_pi_le_sum`, composed with the coordinate split inside `integral_var_term_le`, which
needs `Fintype.card ι` in the literal form `n + 1`; that composition is not carried out here. -/
theorem variance_pi_le_of_bounded_diff_fin {D : ℕ} {α : Fin D → Type*}
    [∀ i, MeasurableSpace (α i)] (ν : ∀ i, Measure (α i)) [∀ i, IsProbabilityMeasure (ν i)]
    (f : (∀ i, α i) → ℝ) (hf : MemLp f 2 (Measure.pi ν)) (c : Fin D → ℝ)
    (hc : ∀ i x t, |f (Function.update x i t) - f x| ≤ c i) :
    ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν) ≤ ∑ i, c i ^ 2 := by
  have hf1 : Integrable f (Measure.pi ν) := hf.integrable (by norm_num)
  have hf2 : Integrable (fun x => f x ^ 2) (Measure.pi ν) := hf.integrable_sq
  cases D with
  | zero =>
    have hmean : ∀ x, ∫ y, f y ∂(Measure.pi ν) = f x := fun x => by
      have heq : (fun y => f y) = fun _ => f x := funext fun y => by
        congr 1; exact Subsingleton.elim y x
      rw [heq]; simp
    have hz : ∀ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 = 0 := fun x => by rw [hmean x]; ring
    simp [hz]
  | succ n =>
    calc ∫ x, (f x - ∫ y, f y ∂(Measure.pi ν)) ^ 2 ∂(Measure.pi ν)
        ≤ ∑ i, ∫ x, (∫ t, (f (Function.update x i t)
            - ∫ s, f (Function.update x i s) ∂(ν i)) ^ 2 ∂(ν i)) ∂(Measure.pi ν) :=
      variance_pi_le_sum_fin (n + 1) ν f hf
      _ ≤ ∑ i, c i ^ 2 :=
        Finset.sum_le_sum fun i _ => integral_var_term_le ν hf1 hf2 i (c i) (hc i)

end Tensorization
end StackedSVD
