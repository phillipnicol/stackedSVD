/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# Convergence in probability to a constant, one space per `N`

STATUS 2026-08-29: rewritten for the per-`N` definition of `Defs.lean` v2; compiled on
Mathlib v4.33.0 (lake env lean, server). Every proof is complete. The file is clean under
the Mathlib standard linter set (`weak.linter.mathlibStandardSet=true`).

This file proves the closure lemmas of `StackedSVD.TendstoInProb` of `Defs.lean`. That
definition gives one probability space `(Ω N, μ N)` for each index `N`. No proof can
therefore compare two values of `N` on one space, and no proof can use an almost sure limit.

## Content

1. Measure helpers (public, section `MeasureHelpers`): a monotone squeeze, an almost
   everywhere squeeze, a union bound for two sets and for a finite family, the complement of
   a full family, and the good-set form. One copy for the whole project.
2. A private predicate `TendstoInProbDist` for a limit in a pseudometric space. It carries
   the continuous mapping theorem for `ℝ`, for `ℝ × ℝ`, and for `ι → ℝ`.
3. The public lemmas: `const`, `congr`, `of_le`, `add`, `neg`, `sub`, `mul`, `const_mul`,
   `mul_const`, `inv`, `div`, `comp_continuous`, `comp_continuous₂`,
   `of_tendsto_measure_ne_of_tendsto`, and `sup_le`.

## Measurability

Almost every lemma here needs no measurable set and no measurable function. The proofs use
monotonicity, subadditivity, or almost everywhere monotonicity of a measure. These hold for
arbitrary sets.

`tendsto_measure_compl_zero` is the exception. It needs the set to be null measurable.
Reason: `μ N s` is an outer measure on an arbitrary set. A saturated nonmeasurable set `s`
and its complement can both have measure 1 in a probability space. Then `μ N s = 1` gives no
bound on the measure of the complement. Every route from a full measure event therefore goes
through `of_tendsto_measure_ne`, which states the fact on the bad event `{ω | f N ω ≠ a}` and
needs no hypothesis. Cleanup wave 3 deleted the consumer-less `of_tendsto_prob_one`, which
was the only lemma here that started from `𝓝 1` (mechanical audit 2026-08-31, finding 7).
-/

open MeasureTheory Filter Topology Set
open scoped ENNReal

namespace StackedSVD

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
variable {f f' g : ∀ N, Ω N → ℝ} {a b : ℝ}

/-! ### Measure helpers

These eight lemmas are the measure facts that the whole project uses about a family of events
that vanishes, or that fills the space, as `N → ∞`. They are public and they live here so
that one copy serves every consumer: before cleanup wave 1 (2026-08-30) the same statements
sat private in this file, public in `namespace SpikedModel` of `RMT/R5.lean`, and private in
`RMT/T.lean`.

None of them needs a measurable set, except `tendsto_measure_compl_zero`. That one needs
null measurability because `μ N` is an outer measure on an arbitrary set: a saturated
nonmeasurable set and its complement can both have measure `1`. -/

section MeasureHelpers

variable {s t : ∀ N, Set (Ω N)}

/-- Squeeze: a family of sets inside a null tending family is null tending. -/
theorem tendsto_measure_zero_of_subset (hst : ∀ N, s N ⊆ t N)
    (ht : Tendsto (fun N => μ N (t N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (s N)) atTop (𝓝 0) :=
  tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds ht (fun _ => zero_le)
    fun N => measure_mono (hst N)

/-- Same as `tendsto_measure_zero_of_subset`, with an almost everywhere inclusion. -/
theorem tendsto_measure_zero_of_subset_ae
    (hst : ∀ N, s N ≤ᵐ[μ N] t N) (ht : Tendsto (fun N => μ N (t N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (s N)) atTop (𝓝 0) :=
  tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds ht (fun _ => zero_le)
    fun N => measure_mono_ae (hst N)

/-- A union bound for two families of sets. -/
theorem tendsto_measure_zero_union
    (hs : Tendsto (fun N => μ N (s N)) atTop (𝓝 0))
    (ht : Tendsto (fun N => μ N (t N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (s N ∪ t N)) atTop (𝓝 0) := by
  have hsum : Tendsto (fun N => μ N (s N) + μ N (t N)) atTop (𝓝 0) := by
    simpa using hs.add ht
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum (fun _ => zero_le)
    fun N => measure_union_le _ _

/-- A union bound over a finite index type. -/
theorem tendsto_measure_zero_iUnion {ι : Type*} [Finite ι] {u : ι → ∀ N, Set (Ω N)}
    (hs : ∀ i, Tendsto (fun N => μ N (u i N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (⋃ i, u i N)) atTop (𝓝 0) := by
  have : Fintype ι := Fintype.ofFinite ι
  have hsum : Tendsto (fun N => ∑ i, μ N (u i N)) atTop (𝓝 0) := by
    have := tendsto_finsetSum (Finset.univ : Finset ι) fun i (_ : i ∈ Finset.univ) => hs i
    simpa using this
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum (fun _ => zero_le)
    fun N => measure_iUnion_fintype_le _ _

/-- The complement of a null measurable family whose measure tends to `1` is null tending.
A `MeasurableSet` hypothesis reaches this form through `MeasurableSet.nullMeasurableSet`. -/
theorem tendsto_measure_compl_zero [∀ N, IsProbabilityMeasure (μ N)]
    (hm : ∀ N, NullMeasurableSet (s N) (μ N))
    (h : Tendsto (fun N => μ N (s N)) atTop (𝓝 1)) :
    Tendsto (fun N => μ N (s N)ᶜ) atTop (𝓝 0) := by
  have key : ∀ N, μ N (s N)ᶜ = 1 - μ N (s N) := fun N => prob_compl_eq_one_sub₀ (hm N)
  simp only [key]
  have h1 : Tendsto (fun N => (1 : ℝ≥0∞) - μ N (s N)) atTop (𝓝 (1 - 1)) :=
    ENNReal.Tendsto.sub tendsto_const_nhds h (Or.inl ENNReal.one_ne_top)
  simpa using h1

/-- From a bad family of vanishing measure to a good family of measure tending to `1`. No
measurability of the good family is needed. -/
theorem tendsto_measure_one_of_bad [∀ N, IsProbabilityMeasure (μ N)]
    (hst : ∀ N, (t N)ᶜ ⊆ s N) (hs : Tendsto (fun N => μ N (s N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (t N)) atTop (𝓝 1) := by
  have hone : Tendsto (fun _ : ℕ => (1 : ℝ≥0∞)) atTop (𝓝 1) := tendsto_const_nhds
  have hg : Tendsto (fun N => 1 - μ N (s N)) atTop (𝓝 1) := by
    have h1 : Tendsto (fun N => (1 : ℝ≥0∞) - μ N (s N)) atTop (𝓝 (1 - 0)) :=
      ENNReal.Tendsto.sub tendsto_const_nhds hs (Or.inl ENNReal.one_ne_top)
    simpa using h1
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le hg hone (fun N => ?_)
    (fun N => prob_le_one)
  · refine tsub_le_iff_right.mpr ?_
    calc (1 : ℝ≥0∞) = μ N Set.univ := measure_univ.symm
      _ = μ N (t N ∪ (t N)ᶜ) := by rw [Set.union_compl_self]
      _ ≤ μ N (t N) + μ N ((t N)ᶜ) := measure_union_le _ _
      _ ≤ μ N (t N) + μ N (s N) := add_le_add le_rfl (measure_mono (hst N))

/-! ### Skeletons for `TendstoInProb`

Two `δ`-by-`δ` skeletons: the target set of a limit in probability sits inside a union of two
or of three families of vanishing measure. Only monotonicity and subadditivity enter, so no
measurability is needed. Cleanup wave 2 moved them here from `RMT/T.lean`
(`tendstoInProb_of_subset_union₃`, private) and `RMT/R3minus.lean`
(`tendstoInProb_of_subset_union₂`, private). -/

/-- A target set inside a union of two null tending families. -/
theorem tendstoInProb_of_subset_union₂ {f : ∀ N, Ω N → ℝ} {ℓ : ℝ}
    (h : ∀ δ > 0, ∃ s t : ∀ N, Set (Ω N),
      (∀ N, {ω | δ ≤ |f N ω - ℓ|} ⊆ s N ∪ t N) ∧
        Tendsto (fun N => μ N (s N)) atTop (𝓝 0) ∧
        Tendsto (fun N => μ N (t N)) atTop (𝓝 0)) :
    TendstoInProb μ f ℓ := by
  intro δ hδ
  obtain ⟨s, t, hsub, hs, ht⟩ := h δ hδ
  have hsum : Tendsto (fun N => μ N (s N) + μ N (t N)) atTop (𝓝 0) := by
    simpa using hs.add ht
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum
    (fun _ => zero_le) fun N => ?_
  calc μ N {ω | δ ≤ |f N ω - ℓ|} ≤ μ N (s N ∪ t N) := measure_mono (hsub N)
    _ ≤ μ N (s N) + μ N (t N) := measure_union_le _ _

/-- A target set inside a union of three null tending families. -/
theorem tendstoInProb_of_subset_union₃ {f : ∀ N, Ω N → ℝ} {ℓ : ℝ}
    (h : ∀ δ > 0, ∃ s t u : ∀ N, Set (Ω N),
      (∀ N, {ω | δ ≤ |f N ω - ℓ|} ⊆ s N ∪ t N ∪ u N) ∧
      Tendsto (fun N => μ N (s N)) atTop (𝓝 0) ∧
      Tendsto (fun N => μ N (t N)) atTop (𝓝 0) ∧
      Tendsto (fun N => μ N (u N)) atTop (𝓝 0)) :
    TendstoInProb μ f ℓ := by
  intro δ hδ
  obtain ⟨s, t, u, hsub, hs, ht, hu⟩ := h δ hδ
  have hsum : Tendsto (fun N => μ N (s N) + μ N (t N) + μ N (u N)) atTop (𝓝 0) := by
    simpa using (hs.add ht).add hu
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum
    (fun _ => zero_le) fun N => ?_
  calc μ N {ω | δ ≤ |f N ω - ℓ|} ≤ μ N (s N ∪ t N ∪ u N) := measure_mono (hsub N)
    _ ≤ μ N (s N ∪ t N) + μ N (u N) := measure_union_le _ _
    _ ≤ μ N (s N) + μ N (t N) + μ N (u N) := add_le_add (measure_union_le _ _) le_rfl

end MeasureHelpers


/-! ### An auxiliary limit in a pseudometric space -/

/-- Convergence in probability to a constant `c` of a sequence with values in a pseudometric
space. Auxiliary: it carries the continuous mapping theorem for `ℝ × ℝ` and for `ι → ℝ`. -/
private def TendstoInProbDist {E : Type*} [PseudoMetricSpace E] (μ : ∀ N, Measure (Ω N))
    (F : ∀ N, Ω N → E) (c : E) : Prop :=
  ∀ ε > 0, Tendsto (fun N => μ N {ω | ε ≤ dist (F N ω) c}) atTop (𝓝 0)

/-- On `ℝ` the auxiliary predicate is the definition. -/
private theorem tendstoInProbDist_real_iff : TendstoInProbDist μ f a ↔ TendstoInProb μ f a := by
  simp only [TendstoInProbDist, TendstoInProb, Real.dist_eq]

/-- Continuous mapping theorem, generic form. -/
private theorem tendstoInProb_comp {E : Type*} [PseudoMetricSpace E] {F : ∀ N, Ω N → E} {c : E}
    {φ : E → ℝ} (hφ : ContinuousAt φ c) (hF : TendstoInProbDist μ F c) :
    TendstoInProb μ (fun N ω => φ (F N ω)) (φ c) := by
  intro ε hε
  obtain ⟨δ, hδ, hδφ⟩ := Metric.continuousAt_iff.mp hφ ε hε
  refine tendsto_measure_zero_of_subset (t := fun N => {ω | δ ≤ dist (F N ω) c}) ?_ (hF δ hδ)
  intro N ω hω
  have hω' : ε ≤ |φ (F N ω) - φ c| := hω
  change δ ≤ dist (F N ω) c
  by_contra hlt
  rw [not_le] at hlt
  have hd : dist (φ (F N ω)) (φ c) < ε := hδφ hlt
  rw [Real.dist_eq] at hd
  exact absurd hω' (not_le.mpr hd)

/-- A pair of convergent sequences converges in the product metric. -/
private theorem tendstoInProbDist_prod (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) :
    TendstoInProbDist μ (fun N ω => (f N ω, g N ω)) (a, b) := by
  intro ε hε
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ε ≤ |f N ω - a|} ∪ {ω | ε ≤ |g N ω - b|}) ?_
    (tendsto_measure_zero_union (hf ε hε) (hg ε hε))
  intro N ω hω
  have hω' : ε ≤ dist (f N ω, g N ω) (a, b) := hω
  rw [Prod.dist_eq, le_max_iff, Real.dist_eq, Real.dist_eq] at hω'
  exact hω'

/-! ### Basic lemmas -/

/-- A constant sequence converges in probability to that constant. -/
theorem TendstoInProb.const (μ : ∀ N, Measure (Ω N)) (a : ℝ) :
    TendstoInProb μ (fun _ _ => a) a := by
  intro ε hε
  have hempty : ∀ N, {ω : Ω N | ε ≤ |a - a|} = (∅ : Set (Ω N)) := by
    intro N
    ext ω
    simp [not_le.mpr hε]
  simp only [hempty, measure_empty]
  exact tendsto_const_nhds

/-- Transfer along an almost everywhere equality for each `N`. -/
theorem TendstoInProb.congr (hf : TendstoInProb μ f a) (h : ∀ N, f N =ᵐ[μ N] f' N) :
    TendstoInProb μ f' a := by
  intro ε hε
  refine tendsto_measure_zero_of_subset_ae (t := fun N => {ω | ε ≤ |f N ω - a|}) ?_ (hf ε hε)
  intro N
  filter_upwards [h N] with ω hω
  intro hmem
  have hmem' : ε ≤ |f' N ω - a| := hmem
  change ε ≤ |f N ω - a|
  rw [hω]
  exact hmem'

/-- Domination: if `|f N ω - a| ≤ g N ω` almost everywhere and `g` tends to `0` in
probability, then `f` tends to `a` in probability. -/
theorem TendstoInProb.of_le (hfg : ∀ N, ∀ᵐ ω ∂(μ N), |f N ω - a| ≤ g N ω)
    (hg : TendstoInProb μ g 0) : TendstoInProb μ f a := by
  intro ε hε
  refine tendsto_measure_zero_of_subset_ae (t := fun N => {ω | ε ≤ |g N ω - 0|}) ?_ (hg ε hε)
  intro N
  filter_upwards [hfg N] with ω hω
  intro hmem
  have hmem' : ε ≤ |f N ω - a| := hmem
  change ε ≤ |g N ω - 0|
  rw [sub_zero]
  exact hmem'.trans (hω.trans (le_abs_self _))

/-- If the probability that `f N` differs from `a` tends to `0`, then `f` tends to `a` in
probability. No measurability and no finiteness of `μ N` is needed. -/
theorem TendstoInProb.of_tendsto_measure_ne
    (h : Tendsto (fun N => μ N {ω | f N ω ≠ a}) atTop (𝓝 0)) : TendstoInProb μ f a := by
  intro ε hε
  refine tendsto_measure_zero_of_subset (t := fun N => {ω | f N ω ≠ a}) ?_ h
  intro N ω hω
  have hω' : ε ≤ |f N ω - a| := hω
  change f N ω ≠ a
  intro heq
  rw [heq, sub_self, abs_zero] at hω'
  exact absurd hω' (not_le.mpr hε)

/-- Transfer along an event of vanishing probability. If `f N` and `g N` agree outside a set
whose measure tends to `0`, and `g` tends to `a` in probability, then so does `f`. This is the
form that a "with probability tending to one" identity gives; `congr` needs an almost
everywhere identity for each `N`, which such an event cannot supply. No measurability is
needed: `{ε ≤ |f N - a|} ⊆ {f N ≠ g N} ∪ {ε ≤ |g N - a|}` and the union bound. -/
theorem TendstoInProb.of_tendsto_measure_ne_of_tendsto
    (h : Tendsto (fun N => μ N {ω | f N ω ≠ g N ω}) atTop (𝓝 0))
    (hg : TendstoInProb μ g a) : TendstoInProb μ f a := by
  intro ε hε
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | f N ω ≠ g N ω} ∪ {ω | ε ≤ |g N ω - a|}) ?_
    (tendsto_measure_zero_union h (hg ε hε))
  intro N ω hω
  have hω' : ε ≤ |f N ω - a| := hω
  by_cases hfg : f N ω = g N ω
  · right
    change ε ≤ |g N ω - a|
    rw [← hfg]
    exact hω'
  · left
    exact hfg

/-! ### Composition with a continuous function -/

/-- Continuous mapping theorem for a constant limit, one variable. -/
theorem TendstoInProb.comp_continuous {φ : ℝ → ℝ} (hφ : ContinuousAt φ a)
    (hf : TendstoInProb μ f a) : TendstoInProb μ (fun N ω => φ (f N ω)) (φ a) :=
  tendstoInProb_comp hφ (tendstoInProbDist_real_iff.mpr hf)

/-- Continuous mapping theorem for a constant limit, two variables. -/
theorem TendstoInProb.comp_continuous₂ {φ : ℝ × ℝ → ℝ} (hφ : ContinuousAt φ (a, b))
    (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) :
    TendstoInProb μ (fun N ω => φ (f N ω, g N ω)) (φ (a, b)) :=
  tendstoInProb_comp hφ (tendstoInProbDist_prod hf hg)

/-! ### Algebra -/

/-- Sums. -/
theorem TendstoInProb.add (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) :
    TendstoInProb μ (fun N ω => f N ω + g N ω) (a + b) :=
  TendstoInProb.comp_continuous₂ (φ := fun p => p.1 + p.2)
    (continuous_fst.add continuous_snd).continuousAt hf hg

/-- Negation. -/
theorem TendstoInProb.neg (hf : TendstoInProb μ f a) :
    TendstoInProb μ (fun N ω => -f N ω) (-a) :=
  TendstoInProb.comp_continuous (φ := fun x => -x) continuous_neg.continuousAt hf

/-- Differences. -/
theorem TendstoInProb.sub (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) :
    TendstoInProb μ (fun N ω => f N ω - g N ω) (a - b) :=
  TendstoInProb.comp_continuous₂ (φ := fun p => p.1 - p.2)
    (continuous_fst.sub continuous_snd).continuousAt hf hg

/-- Products. -/
theorem TendstoInProb.mul (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) :
    TendstoInProb μ (fun N ω => f N ω * g N ω) (a * b) :=
  TendstoInProb.comp_continuous₂ (φ := fun p => p.1 * p.2)
    (continuous_fst.mul continuous_snd).continuousAt hf hg

/-- A constant factor on the left. -/
theorem TendstoInProb.const_mul (c : ℝ) (hf : TendstoInProb μ f a) :
    TendstoInProb μ (fun N ω => c * f N ω) (c * a) :=
  TendstoInProb.comp_continuous (φ := fun x => c * x)
    (continuous_const.mul continuous_id).continuousAt hf

/-- A constant factor on the right. -/
theorem TendstoInProb.mul_const (c : ℝ) (hf : TendstoInProb μ f a) :
    TendstoInProb μ (fun N ω => f N ω * c) (a * c) :=
  TendstoInProb.comp_continuous (φ := fun x => x * c)
    (continuous_id.mul continuous_const).continuousAt hf

/-- Inverses, for a nonzero limit. -/
theorem TendstoInProb.inv (ha : a ≠ 0) (hf : TendstoInProb μ f a) :
    TendstoInProb μ (fun N ω => (f N ω)⁻¹) a⁻¹ :=
  TendstoInProb.comp_continuous (φ := fun x => x⁻¹) (continuousAt_inv₀ ha) hf

/-- Quotients, for a nonzero limit in the denominator. -/
theorem TendstoInProb.div (hf : TendstoInProb μ f a) (hg : TendstoInProb μ g b) (hb : b ≠ 0) :
    TendstoInProb μ (fun N ω => f N ω / g N ω) (a / b) :=
  TendstoInProb.comp_continuous₂ (φ := fun p => p.1 / p.2)
    (continuousAt_fst.div continuousAt_snd (by simpa using hb)) hf hg

/-! ### Finitely many coordinates -/

/-- Coordinatewise convergence in probability of a sequence of vectors indexed by `ι`. -/
def TendstoInProbPi {ι : Type*} (μ : ∀ N, Measure (Ω N)) (F : ∀ N, Ω N → ι → ℝ) (v : ι → ℝ) :
    Prop :=
  ∀ i, TendstoInProb μ (fun N ω => F N ω i) (v i)

section Pi

variable {ι : Type*} {F : ∀ N, Ω N → ι → ℝ} {v : ι → ℝ}

/-- Coordinatewise convergence gives convergence for the sup metric on `ι → ℝ`. -/
private theorem tendstoInProbDist_pi [Fintype ι] (hF : TendstoInProbPi μ F v) :
    TendstoInProbDist μ F v := by
  intro ε hε
  refine tendsto_measure_zero_of_subset (t := fun N => ⋃ i, {ω | ε ≤ |F N ω i - v i|}) ?_
    (tendsto_measure_zero_iUnion fun i => hF i ε hε)
  intro N ω hω
  have hω' : ε ≤ dist (F N ω) v := hω
  by_contra hnot
  simp only [Set.mem_iUnion, Set.mem_ofPred_eq, not_exists, not_le] at hnot
  have hlt : dist (F N ω) v < ε := by
    refine (dist_pi_lt_iff hε).mpr fun i => ?_
    rw [Real.dist_eq]
    exact hnot i
  exact absurd hω' (not_le.mpr hlt)

/-- Continuous mapping theorem for a constant limit, finitely many variables. -/
theorem TendstoInProbPi.comp_continuous [Finite ι] {φ : (ι → ℝ) → ℝ} (hφ : ContinuousAt φ v)
    (hF : TendstoInProbPi μ F v) : TendstoInProb μ (fun N ω => φ (F N ω)) (φ v) := by
  have : Fintype ι := Fintype.ofFinite ι
  exact tendstoInProb_comp hφ (tendstoInProbDist_pi hF)

end Pi

/-! ### A finite sum of sequences -/

/-- **A finite sum of limits in probability is the limit of the sum.** Induction over the
index `Finset` with `TendstoInProb.add` at each step. Four files kept a private copy of this
statement until the 2026-09-02 dedupe (`RankR/Subspace.lean`, `RankR/SubspaceG.lean`,
`RankR/GeneralMain.lean`, `RankR/StackMain.lean`). -/
theorem TendstoInProb.finsum {ι : Type*} [Fintype ι] {F : ι → ∀ N, Ω N → ℝ} {v : ι → ℝ}
    (h : ∀ i, TendstoInProb μ (F i) (v i)) :
    TendstoInProb μ (fun N ω => ∑ i, F i N ω) (∑ i, v i) := by
  classical
  have key : ∀ s : Finset ι,
      TendstoInProb μ (fun N ω => ∑ i ∈ s, F i N ω) (∑ i ∈ s, v i) := by
    intro s
    refine Finset.induction_on s ?_ ?_
    · simpa using TendstoInProb.const μ 0
    · intro a t ha ih
      simp only [Finset.sum_insert ha]
      exact (h a).add ih
  exact key Finset.univ

/-! ### The maximum of finitely many sequences -/

/-- A uniform bound on the coordinates bounds the distance of the two suprema. -/
private theorem abs_ciSup_sub_ciSup_le {ι : Type*} [Finite ι] [Nonempty ι] {u v : ι → ℝ}
    {ε : ℝ} (h : ∀ i, |u i - v i| ≤ ε) : |(⨆ i, u i) - ⨆ i, v i| ≤ ε := by
  have key : ∀ x y : ι → ℝ, BddAbove (Set.range y) → (∀ i, x i - y i ≤ ε) →
      (⨆ i, x i) - ⨆ i, y i ≤ ε := by
    intro x y hy hxy
    have hle : (⨆ i, x i) ≤ (⨆ i, y i) + ε := by
      refine ciSup_le fun i => ?_
      have h1 : x i ≤ y i + ε := by linarith [hxy i]
      have h2 : y i ≤ ⨆ j, y j := le_ciSup hy i
      linarith
    linarith
  have hu : BddAbove (Set.range u) := Finite.bddAbove_range u
  have hv : BddAbove (Set.range v) := Finite.bddAbove_range v
  rw [abs_le]
  constructor
  · have hle := key v u hu fun i => by
      have hi := h i
      rw [abs_le] at hi
      linarith [hi.1]
    linarith
  · exact key u v hv fun i => by
      have hi := h i
      rw [abs_le] at hi
      linarith [hi.2]

/-- The supremum over a finite nonempty index type is continuous. -/
private theorem continuousAt_ciSup {ι : Type*} [Finite ι] [Nonempty ι] (v : ι → ℝ) :
    ContinuousAt (fun w : ι → ℝ => ⨆ i, w i) v := by
  have : Fintype ι := Fintype.ofFinite ι
  refine Metric.continuousAt_iff.mpr fun ε hε => ⟨ε / 2, by positivity, fun {w} hw => ?_⟩
  rw [Real.dist_eq]
  refine lt_of_le_of_lt (abs_ciSup_sub_ciSup_le (ε := ε / 2) fun i => ?_) (by linarith)
  have hi := (dist_pi_lt_iff (show (0 : ℝ) < ε / 2 by positivity)).mp hw i
  rw [Real.dist_eq] at hi
  exact hi.le

/-- The maximum of finitely many sequences converges to the maximum of the limits. -/
theorem TendstoInProb.sup_le {ι : Type*} [Finite ι] {F : ι → ∀ N, Ω N → ℝ} {v : ι → ℝ}
    (hF : ∀ i, TendstoInProb μ (F i) (v i)) :
    TendstoInProb μ (fun N ω => ⨆ i, F i N ω) (⨆ i, v i) := by
  rcases isEmpty_or_nonempty ι with hι | hι
  · have := hι
    simp only [Real.iSup_of_isEmpty]
    exact TendstoInProb.const μ 0
  · exact TendstoInProbPi.comp_continuous (φ := fun w : ι → ℝ => ⨆ i, w i)
      (continuousAt_ciSup v) fun i => hF i

end StackedSVD
