/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Defs
import StackedSVD.RMT.R3
import StackedSVD.Prob.TendstoInProb

/-! # Stage 3, unit L: the discarded part

With probability tending to 1 the discarded part `R = discardMat T Z` has operator norm at
most `ε √d`: no entry above `ε √d` (the fourth moment, by dominated convergence), and no row
or column holds two discarded entries (a union bound at `T = d^(1/2 - a)`, `a < 1/8`); a
matrix with at most one nonzero entry per row and per column has operator norm at most its
largest absolute entry (`notes/stage3_edge.md`, route C; `notes/STAGE3_CAMPAIGN.md`). -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace Edge

/-! ### L1: the deterministic bound at one nonzero entry per line -/

/-- A finite family with at most one nonzero member: the square of the sum is the sum of the
squares. -/
private theorem sq_sum_of_atMostOne {ι : Type*} [Fintype ι] (f : ι → ℝ)
    (h : ∀ i i', f i ≠ 0 → f i' ≠ 0 → i = i') :
    (∑ i, f i) ^ 2 = ∑ i, f i ^ 2 := by
  by_cases hall : ∀ i, f i = 0
  · simp [hall]
  · obtain ⟨i₀, hi₀⟩ := not_forall.mp hall
    have hz : ∀ i, i ≠ i₀ → f i = 0 := by
      intro i hi
      by_contra hfi
      exact hi (h i i₀ hfi hi₀)
    have h1 : ∑ i, f i = f i₀ :=
      Finset.sum_eq_single i₀ (fun i _ hi => hz i hi) fun hc => absurd (Finset.mem_univ _) hc
    have h2 : ∑ i, f i ^ 2 = f i₀ ^ 2 :=
      Finset.sum_eq_single i₀ (fun i _ hi => by rw [hz i hi]; ring)
        fun hc => absurd (Finset.mem_univ _) hc
    rw [h1, h2]

/-- A finite family with at most one nonzero member, every member at most `B`, sums to at
most `B`. -/
private theorem sum_le_of_atMostOne {ι : Type*} [Fintype ι] {f : ι → ℝ} {B : ℝ} (hB : 0 ≤ B)
    (h : ∀ i i', f i ≠ 0 → f i' ≠ 0 → i = i') (hb : ∀ i, f i ≤ B) :
    ∑ i, f i ≤ B := by
  by_cases hall : ∀ i, f i = 0
  · simpa [hall] using hB
  · obtain ⟨i₀, hi₀⟩ := not_forall.mp hall
    have hz : ∀ i, i ≠ i₀ → f i = 0 := by
      intro i hi
      by_contra hfi
      exact hi (h i i₀ hfi hi₀)
    have h1 : ∑ i, f i = f i₀ :=
      Finset.sum_eq_single i₀ (fun i _ hi => hz i hi) fun hc => absurd (Finset.mem_univ _) hc
    rw [h1]
    exact hb i₀

/-- L1, deterministic. A matrix with at most one nonzero entry per row and per column has
operator norm at most the largest absolute entry. Mathlib v4.33.0 has no Schur test for the
`l2` operator norm, so this is written directly. -/
theorem opNorm_le_of_sparse {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) {b : ℝ} (hb : 0 ≤ b)
    (hrow : ∀ i j j', A i j ≠ 0 → A i j' ≠ 0 → j = j')
    (hcol : ∀ i i' j, A i j ≠ 0 → A i' j ≠ 0 → i = i')
    (hentry : ∀ i j, |A i j| ≤ b) :
    ‖A‖ ≤ b := by
  rw [Matrix.l2_opNorm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ hb fun x => ?_
  have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) A) x
      = WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) := rfl
  rw [happ]
  have hsq : ‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖ ^ 2
      ≤ b ^ 2 * ‖x‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq, EuclideanSpace.real_norm_sq_eq, Finset.mul_sum]
    have hrowsq : ∀ i : Fin n,
        (WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n)) i ^ 2
          = ∑ j, (A i j * x j) ^ 2 := by
      intro i
      have hco : (WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n)) i
          = ∑ j, A i j * x j := rfl
      rw [hco]
      refine sq_sum_of_atMostOne _ fun j j' hj hj' => ?_
      exact hrow i j j' (left_ne_zero_of_mul hj) (left_ne_zero_of_mul hj')
    calc ∑ i, (WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n)) i ^ 2
        = ∑ i, ∑ j, (A i j * x j) ^ 2 := Finset.sum_congr rfl fun i _ => hrowsq i
      _ = ∑ j, ∑ i, (A i j * x j) ^ 2 := Finset.sum_comm
      _ ≤ ∑ j, b ^ 2 * x j ^ 2 := by
          refine Finset.sum_le_sum fun j _ => ?_
          refine sum_le_of_atMostOne (by positivity) (fun i i' hi hi' => ?_) fun i => ?_
          · refine hcol i i' j ?_ ?_
            · exact left_ne_zero_of_mul (pow_ne_zero_iff (n := 2) (by norm_num) |>.mp hi)
            · exact left_ne_zero_of_mul (pow_ne_zero_iff (n := 2) (by norm_num) |>.mp hi')
          · have hA : A i j ^ 2 ≤ b ^ 2 := by
              have := hentry i j
              nlinarith [abs_nonneg (A i j), sq_abs (A i j)]
            have hexp : (A i j * x j) ^ 2 = A i j ^ 2 * x j ^ 2 := by ring
            rw [hexp]
            exact mul_le_mul_of_nonneg_right hA (by positivity)
  calc ‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖
      = Real.sqrt (‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖ ^ 2) :=
        (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt (b ^ 2 * ‖x‖ ^ 2) := Real.sqrt_le_sqrt hsq
    _ = b * ‖x‖ := by
        rw [Real.sqrt_mul (by positivity), Real.sqrt_sq hb, Real.sqrt_sq (norm_nonneg _)]

/-- L1b. The discarded part of a matrix is measurable. -/
theorem measurable_discardMat (T : ℝ) (n d : ℕ) :
    Measurable (discardMat T : Matrix (Fin n) (Fin d) ℝ → Matrix (Fin n) (Fin d) ℝ) := by
  have hmap : Measurable (discardMap T) :=
    Measurable.ite (measurableSet_le measurable_id.abs measurable_const)
      measurable_const measurable_id
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  exact hmap.comp ((measurable_pi_apply j).comp (measurable_pi_apply i))

/-! ### L2a: the fourth moment at the entry level -/

/-- The tail set `{x | s < |x|}` is measurable. -/
private theorem measurableSet_tail (s : ℝ) : MeasurableSet {x : ℝ | s < |x|} :=
  measurableSet_lt measurable_const measurable_id.abs

/-- Markov at order 4, with the truncated integral: `s⁴ ν{|x| > s} ≤ ∫ x⁴ 1{|x| > s} ∂ν`. -/
private theorem measure_tail_mul_le {ν : Measure ℝ} [IsFiniteMeasure ν]
    (hint : Integrable (fun x : ℝ => x ^ 4) ν) {s : ℝ} (hs : 0 ≤ s) :
    (ν {x : ℝ | s < |x|}).toReal * s ^ 4
      ≤ ∫ x, Set.indicator {x : ℝ | s < |x|} (fun x => x ^ 4) x ∂ν := by
  have hS := measurableSet_tail s
  have h1 : ∫ x, Set.indicator {x : ℝ | s < |x|} (fun _ => s ^ 4) x ∂ν
      = (ν {x : ℝ | s < |x|}).toReal * s ^ 4 := by
    rw [integral_indicator_const _ hS, smul_eq_mul, measureReal_def]
  rw [← h1]
  refine integral_mono ((integrable_const (s ^ 4)).indicator hS) (hint.indicator hS) fun x => ?_
  by_cases hx : x ∈ {x : ℝ | s < |x|}
  · simp only [Set.indicator_of_mem hx]
    have hxx : s < |x| := hx
    calc s ^ 4 ≤ |x| ^ 4 := by gcongr
      _ = x ^ 4 := Even.pow_abs (by norm_num) x
  · simp only [Set.indicator_of_notMem hx]
    positivity

/-- Markov at order 4: `ν {|x| > s} ≤ ν₄ / s⁴` for `s > 0`. -/
private theorem measure_tail_le {ν : Measure ℝ} (hν : NoiseLaw ν) {s : ℝ} (hs : 0 < s) :
    (ν {x : ℝ | s < |x|}).toReal ≤ (∫ x, x ^ 4 ∂ν) / s ^ 4 := by
  have := hν.prob
  have hS := measurableSet_tail s
  have h1 := measure_tail_mul_le hν.mom4 hs.le
  have h2 : ∫ x, Set.indicator {x : ℝ | s < |x|} (fun x => x ^ 4) x ∂ν ≤ ∫ x, x ^ 4 ∂ν := by
    refine integral_mono (hν.mom4.indicator hS) hν.mom4 fun x => ?_
    by_cases hx : x ∈ {x : ℝ | s < |x|}
    · simp only [Set.indicator_of_mem hx, le_refl]
    · simp only [Set.indicator_of_notMem hx]
      positivity
  rw [le_div_iff₀ (by positivity)]
  linarith

/-- L2a. The fourth moment is spent here: `n d P(|Z| > ε √d) → 0` by dominated convergence.
Plain Markov gives only `O(1)`. -/
theorem tendsto_entry_tail {ν : Measure ℝ} (hν : NoiseLaw ν) {n d : ℕ → ℕ} {c : ℝ}
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => ((n N : ℝ) * d N) * (ν {x | ε * Real.sqrt (d N) < |x|}).toReal)
      atTop (𝓝 0) := by
  have := hν.prob
  set G : ℕ → ℝ → ℝ := fun N =>
    Set.indicator {x : ℝ | ε * Real.sqrt (d N) < |x|} fun x => x ^ 4 with hGdef
  -- the level `ε √d_N` tends to infinity
  have hlev : Tendsto (fun N => ε * Real.sqrt (d N)) atTop atTop :=
    Filter.Tendsto.const_mul_atTop hε
      (Real.tendsto_sqrt_atTop.comp (tendsto_natCast_atTop_atTop.comp hdtop))
  -- the truncated fourth moment tends to 0 by dominated convergence
  have hGtend : Tendsto (fun N => ∫ x, G N x ∂ν) atTop (𝓝 0) := by
    have hlim : ∀ᵐ x ∂ν, Tendsto (fun N => G N x) atTop (𝓝 ((fun _ : ℝ => (0 : ℝ)) x)) := by
      filter_upwards with x
      refine Tendsto.congr' ?_ tendsto_const_nhds
      filter_upwards [hlev.eventually_gt_atTop |x|] with N hN
      have hx : x ∉ {y : ℝ | ε * Real.sqrt (d N) < |y|} := fun hc =>
        absurd (hc : ε * Real.sqrt (d N) < |x|) (not_lt.mpr hN.le)
      rw [hGdef]
      exact (Set.indicator_of_notMem hx _).symm
    have hdom := tendsto_integral_of_dominated_convergence (F := G) (f := fun _ : ℝ => (0 : ℝ))
      (bound := fun x : ℝ => x ^ 4) (μ := ν)
      (fun N => (hν.mom4.indicator (measurableSet_tail _)).aestronglyMeasurable)
      hν.mom4
      (fun N => by
        filter_upwards with x
        rw [hGdef, Real.norm_eq_abs]
        by_cases hx : x ∈ {y : ℝ | ε * Real.sqrt (d N) < |y|}
        · simp only [Set.indicator_of_mem hx]
          rw [abs_of_nonneg (by positivity)]
        · simp only [Set.indicator_of_notMem hx, abs_zero]
          positivity)
      hlim
    simpa using hdom
  -- the squeeze
  have hbound : Tendsto (fun N => (n N : ℝ) / d N * (∫ x, G N x ∂ν) / ε ^ 4) atTop (𝓝 0) := by
    have h := (hcN.mul hGtend).div_const (ε ^ 4)
    simpa using h
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hbound
    (Eventually.of_forall fun N => by positivity) ?_
  filter_upwards [hdtop.eventually_ge_atTop 1] with N hN
  have hD : (1 : ℝ) ≤ (d N : ℝ) := by exact_mod_cast hN
  have hDpos : (0 : ℝ) < (d N : ℝ) := by linarith
  have hs : 0 < ε * Real.sqrt (d N) := by positivity
  have hm := measure_tail_mul_le (ν := ν) hν.mom4 hs.le
  have hsq : Real.sqrt (d N) ^ 4 = (d N : ℝ) ^ 2 := by
    have h2 : Real.sqrt (d N) ^ 2 = (d N : ℝ) := Real.sq_sqrt hDpos.le
    calc Real.sqrt (d N) ^ 4 = (Real.sqrt (d N) ^ 2) ^ 2 := by ring
      _ = (d N : ℝ) ^ 2 := by rw [h2]
  have hsq4 : (ε * Real.sqrt (d N)) ^ 4 = ε ^ 4 * (d N : ℝ) ^ 2 := by rw [mul_pow, hsq]
  rw [hsq4] at hm
  have hq : (0 : ℝ) ≤ (ν {x : ℝ | ε * Real.sqrt (d N) < |x|}).toReal := ENNReal.toReal_nonneg
  have hnn : (0 : ℝ) ≤ (n N : ℝ) := Nat.cast_nonneg _
  have hprod := mul_le_mul_of_nonneg_left hm hnn
  rw [div_mul_eq_mul_div, div_div, le_div_iff₀ (by positivity)]
  nlinarith [hprod, hq, hnn, hDpos]

/-! ### L2: the three bad events and the assembly -/

/-- The law of one entry under `noiseMatrix` is `ν`. -/
private theorem measure_entry_eq {ν : Measure ℝ} [IsProbabilityMeasure ν] {n d : ℕ}
    (i : Fin n) (j : Fin d) {S : Set ℝ} (hS : MeasurableSet S) :
    noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ | Z i j ∈ S} = ν S := by
  have h1 : MeasurePreserving (fun Z : Matrix (Fin n) (Fin d) ℝ => Z i)
      (noiseMatrix ν n d) (Measure.pi fun _ : Fin d => ν) :=
    measurePreserving_eval (fun _ : Fin n => Measure.pi fun _ : Fin d => ν) i
  have h2 : MeasurePreserving (fun y : Fin d → ℝ => y j)
      (Measure.pi fun _ : Fin d => ν) ν := measurePreserving_eval _ j
  exact (h2.comp h1).measure_preimage hS.nullMeasurableSet

/-- Two entries of one row are independent: the two-entry event has measure `(ν S)²`. -/
private theorem measure_row_pair {ν : Measure ℝ} [IsProbabilityMeasure ν] {n d : ℕ}
    (i : Fin n) {j j' : Fin d} (hjj : j ≠ j') {S : Set ℝ} (hS : MeasurableSet S) :
    noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ | Z i j ∈ S ∧ Z i j' ∈ S} = ν S * ν S := by
  have hrow : MeasurePreserving (fun Z : Matrix (Fin n) (Fin d) ℝ => Z i)
      (noiseMatrix ν n d) (Measure.pi fun _ : Fin d => ν) :=
    measurePreserving_eval (fun _ : Fin n => Measure.pi fun _ : Fin d => ν) i
  have hmeasF : MeasurableSet {y : Fin d → ℝ | y j ∈ S ∧ y j' ∈ S} :=
    ((measurable_pi_apply j) hS).inter ((measurable_pi_apply j') hS)
  have hinner : (Measure.pi fun _ : Fin d => ν) {y : Fin d → ℝ | y j ∈ S ∧ y j' ∈ S}
      = ν S * ν S := by
    have h := (PiLaw.indepFun_coord ν hjj).measure_inter_preimage_eq_mul S S hS hS
    rw [(PiLaw.measurePreserving_coord ν j).measure_preimage hS.nullMeasurableSet,
      (PiLaw.measurePreserving_coord ν j').measure_preimage hS.nullMeasurableSet] at h
    exact h
  exact (hrow.measure_preimage hmeasF.nullMeasurableSet).trans hinner

/-- Two entries of one column are independent: the two-entry event has measure `(ν S)²`. -/
private theorem measure_col_pair {ν : Measure ℝ} [IsProbabilityMeasure ν] {n d : ℕ}
    {i i' : Fin n} (hii : i ≠ i') (j : Fin d) {S : Set ℝ} (hS : MeasurableSet S) :
    noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ | Z i j ∈ S ∧ Z i' j ∈ S} = ν S * ν S := by
  have hS₀ : MeasurableSet ((fun y : Fin d → ℝ => y j) ⁻¹' S) := (measurable_pi_apply j) hS
  have h := IndepFun.measure_inter_preimage_eq_mul
    (PiLaw.indepFun_coord (Measure.pi fun _ : Fin d => ν) hii) _ _ hS₀ hS₀
  have e : ∀ k : Fin n, (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => ν)
      ((fun Z : Fin n → Fin d → ℝ => Z k) ⁻¹' ((fun y : Fin d → ℝ => y j) ⁻¹' S)) = ν S := by
    intro k
    have h1 : MeasurePreserving (fun Z : Fin n → Fin d → ℝ => Z k)
        (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => ν)
        (Measure.pi fun _ : Fin d => ν) :=
      measurePreserving_eval (fun _ : Fin n => Measure.pi fun _ : Fin d => ν) k
    rw [h1.measure_preimage hS₀.nullMeasurableSet]
    exact (PiLaw.measurePreserving_coord ν j).measure_preimage hS.nullMeasurableSet
  rw [e i, e i'] at h
  exact h

/-- A union over a finite index type, every event of measure at most `m`. -/
private theorem measure_iUnion_le_card {α : Type*} [MeasurableSpace α] (μ : Measure α)
    {ι : Type*} [Fintype ι] (u : ι → Set α) {m : ℝ≥0∞} (hu : ∀ i, μ (u i) ≤ m) :
    μ (⋃ i, u i) ≤ (Fintype.card ι : ℝ≥0∞) * m := by
  refine le_trans (measure_iUnion_fintype_le _ _) ?_
  calc ∑ i, μ (u i) ≤ ∑ _i : ι, m := Finset.sum_le_sum fun i _ => hu i
    _ = (Fintype.card ι : ℝ≥0∞) * m := by
        rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

/-- A cast times a finite measure value, against a real bound. -/
private theorem natCast_mul_le_ofReal {k : ℕ} {m : ℝ≥0∞} (hm : m ≠ ∞) {r : ℝ}
    (h : (k : ℝ) * m.toReal ≤ r) : (k : ℝ≥0∞) * m ≤ ENNReal.ofReal r := by
  rw [← ENNReal.ofReal_natCast k, ← ENNReal.ofReal_toReal hm,
    ← ENNReal.ofReal_mul (Nat.cast_nonneg k)]
  exact ENNReal.ofReal_le_ofReal h

/-- A family of events of measure at most `ENNReal.ofReal (g N)`, with `g` tending to `0`,
is null tending. -/
private theorem tendsto_measure_zero_of_le_ofReal {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {s : ∀ N, Set (Ω N)} {g : ℕ → ℝ}
    (hle : ∀ N, μ N (s N) ≤ ENNReal.ofReal (g N)) (hg : Tendsto g atTop (𝓝 0)) :
    Tendsto (fun N => μ N (s N)) atTop (𝓝 0) := by
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds ?_ (fun _ => zero_le) hle
  simpa [Function.comp_def] using (ENNReal.continuous_ofReal.tendsto 0).comp hg

/-- Bad event 1: some entry is above `b` in absolute value. -/
private def bigEntrySet (b : ℝ) (n d : ℕ) : Set (Matrix (Fin n) (Fin d) ℝ) :=
  ⋃ p : Fin n × Fin d, {Z | b < |Z p.1 p.2|}

/-- Bad event 2: one row holds two entries above `T` in absolute value. -/
private def rowPairSet (T : ℝ) (n d : ℕ) : Set (Matrix (Fin n) (Fin d) ℝ) :=
  ⋃ p : Fin n × Fin d × Fin d, {Z | T < |Z p.1 p.2.1| ∧ T < |Z p.1 p.2.2| ∧ p.2.1 ≠ p.2.2}

/-- Bad event 3: one column holds two entries above `T` in absolute value. -/
private def colPairSet (T : ℝ) (n d : ℕ) : Set (Matrix (Fin n) (Fin d) ℝ) :=
  ⋃ p : Fin n × Fin n × Fin d, {Z | T < |Z p.1 p.2.2| ∧ T < |Z p.2.1 p.2.2| ∧ p.1 ≠ p.2.1}

/-- The union bound on bad event 1: `n d` entries, each above `b` with probability
`ν {|x| > b}`. -/
private theorem measure_bigEntrySet_le {ν : Measure ℝ} [IsProbabilityMeasure ν] (n d : ℕ)
    (b : ℝ) :
    noiseMatrix ν n d (bigEntrySet b n d)
      ≤ ENNReal.ofReal ((n : ℝ) * d * (ν {x : ℝ | b < |x|}).toReal) := by
  have hb : ∀ p : Fin n × Fin d,
      noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ | b < |Z p.1 p.2|}
        ≤ ν {x : ℝ | b < |x|} :=
    fun p => le_of_eq (measure_entry_eq p.1 p.2 (measurableSet_tail b))
  simp only [bigEntrySet]
  refine le_trans (measure_iUnion_le_card _ _ hb) (natCast_mul_le_ofReal (measure_ne_top _ _) ?_)
  simp only [Fintype.card_prod, Fintype.card_fin]
  push_cast
  ring_nf
  exact le_rfl

/-- The union bound on bad event 2: `n d²` ordered row pairs, each of probability `p²` by the
independence of two entries of one row. -/
private theorem measure_rowPairSet_le {ν : Measure ℝ} [IsProbabilityMeasure ν] (n d : ℕ)
    (T : ℝ) :
    noiseMatrix ν n d (rowPairSet T n d)
      ≤ ENNReal.ofReal ((n : ℝ) * d * d * (ν {x : ℝ | T < |x|}).toReal ^ 2) := by
  have hS := measurableSet_tail T
  have hb : ∀ p : Fin n × Fin d × Fin d,
      noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ |
          T < |Z p.1 p.2.1| ∧ T < |Z p.1 p.2.2| ∧ p.2.1 ≠ p.2.2}
        ≤ ν {x : ℝ | T < |x|} * ν {x : ℝ | T < |x|} := by
    intro p
    by_cases hjj : p.2.1 = p.2.2
    · have he : {Z : Matrix (Fin n) (Fin d) ℝ |
          T < |Z p.1 p.2.1| ∧ T < |Z p.1 p.2.2| ∧ p.2.1 ≠ p.2.2} = ∅ := by
        ext Z
        simp [hjj]
      rw [he]
      simp
    · have he : {Z : Matrix (Fin n) (Fin d) ℝ |
            T < |Z p.1 p.2.1| ∧ T < |Z p.1 p.2.2| ∧ p.2.1 ≠ p.2.2}
          = {Z : Matrix (Fin n) (Fin d) ℝ |
            Z p.1 p.2.1 ∈ {x : ℝ | T < |x|} ∧ Z p.1 p.2.2 ∈ {x : ℝ | T < |x|}} := by
        ext Z
        simp [hjj]
      rw [he, measure_row_pair p.1 hjj hS]
  simp only [rowPairSet]
  refine le_trans (measure_iUnion_le_card _ _ hb)
    (natCast_mul_le_ofReal (ENNReal.mul_ne_top (measure_ne_top _ _) (measure_ne_top _ _))
      (le_of_eq ?_))
  rw [ENNReal.toReal_mul]
  simp only [Fintype.card_prod, Fintype.card_fin]
  push_cast
  ring

/-- The union bound on bad event 3: `n² d` ordered column pairs, each of probability `p²` by
the independence of two entries of one column. -/
private theorem measure_colPairSet_le {ν : Measure ℝ} [IsProbabilityMeasure ν] (n d : ℕ)
    (T : ℝ) :
    noiseMatrix ν n d (colPairSet T n d)
      ≤ ENNReal.ofReal ((n : ℝ) * n * d * (ν {x : ℝ | T < |x|}).toReal ^ 2) := by
  have hS := measurableSet_tail T
  have hb : ∀ p : Fin n × Fin n × Fin d,
      noiseMatrix ν n d {Z : Matrix (Fin n) (Fin d) ℝ |
          T < |Z p.1 p.2.2| ∧ T < |Z p.2.1 p.2.2| ∧ p.1 ≠ p.2.1}
        ≤ ν {x : ℝ | T < |x|} * ν {x : ℝ | T < |x|} := by
    intro p
    by_cases hii : p.1 = p.2.1
    · have he : {Z : Matrix (Fin n) (Fin d) ℝ |
          T < |Z p.1 p.2.2| ∧ T < |Z p.2.1 p.2.2| ∧ p.1 ≠ p.2.1} = ∅ := by
        ext Z
        simp [hii]
      rw [he]
      simp
    · have he : {Z : Matrix (Fin n) (Fin d) ℝ |
            T < |Z p.1 p.2.2| ∧ T < |Z p.2.1 p.2.2| ∧ p.1 ≠ p.2.1}
          = {Z : Matrix (Fin n) (Fin d) ℝ |
            Z p.1 p.2.2 ∈ {x : ℝ | T < |x|} ∧ Z p.2.1 p.2.2 ∈ {x : ℝ | T < |x|}} := by
        ext Z
        simp [hii]
      rw [he, measure_col_pair hii p.2.2 hS]
  simp only [colPairSet]
  refine le_trans (measure_iUnion_le_card _ _ hb)
    (natCast_mul_le_ofReal (ENNReal.mul_ne_top (measure_ne_top _ _) (measure_ne_top _ _))
      (le_of_eq ?_))
  rw [ENNReal.toReal_mul]
  simp only [Fintype.card_prod, Fintype.card_fin]
  push_cast
  ring

/-- The Markov bound at the truncation level, in the form the two pair bounds consume:
`D³ q² ≤ ν₄² D^(8a-1)` at `q = ν{|x| > T_D}` and `T_D = D^(1/2-a)`, since
`D³ = (T_D)⁸ D^(8a-1)`. -/
private theorem cube_mul_tail_sq_le {ν : Measure ℝ} (hν : NoiseLaw ν) {a : ℝ} {D : ℕ}
    (hD : 0 < D) :
    (D : ℝ) ^ 3 * (ν {x : ℝ | truncLevel a D < |x|}).toReal ^ 2
      ≤ (∫ x, x ^ 4 ∂ν) ^ 2 * (D : ℝ) ^ (8 * a - 1) := by
  have hDpos : (0 : ℝ) < (D : ℝ) := by exact_mod_cast hD
  have hT : 0 < truncLevel a D := truncLevel_pos hD
  have hq : (0 : ℝ) ≤ (ν {x : ℝ | truncLevel a D < |x|}).toReal := ENNReal.toReal_nonneg
  have hq4 : (ν {x : ℝ | truncLevel a D < |x|}).toReal * truncLevel a D ^ 4 ≤ ∫ x, x ^ 4 ∂ν :=
    (le_div_iff₀ (by positivity)).mp (measure_tail_le hν hT)
  have hkey : (D : ℝ) ^ 3 = truncLevel a D ^ 8 * (D : ℝ) ^ (8 * a - 1) := by
    simp only [truncLevel]
    rw [← Real.rpow_natCast ((D : ℝ) ^ ((1 : ℝ) / 2 - a)) 8, ← Real.rpow_mul hDpos.le,
      ← Real.rpow_add hDpos, ← Real.rpow_natCast (D : ℝ) 3]
    congr 1
    push_cast
    ring
  have hE : (0 : ℝ) ≤ (D : ℝ) ^ (8 * a - 1) := Real.rpow_nonneg hDpos.le _
  have hsq : ((ν {x : ℝ | truncLevel a D < |x|}).toReal * truncLevel a D ^ 4) ^ 2
      ≤ (∫ x, x ^ 4 ∂ν) ^ 2 := by
    have h0 : 0 ≤ (ν {x : ℝ | truncLevel a D < |x|}).toReal * truncLevel a D ^ 4 := by
      positivity
    nlinarith [hq4, h0]
  calc (D : ℝ) ^ 3 * (ν {x : ℝ | truncLevel a D < |x|}).toReal ^ 2
      = (truncLevel a D ^ 8 * (D : ℝ) ^ (8 * a - 1))
          * (ν {x : ℝ | truncLevel a D < |x|}).toReal ^ 2 := by rw [hkey]
    _ = (D : ℝ) ^ (8 * a - 1)
          * ((ν {x : ℝ | truncLevel a D < |x|}).toReal * truncLevel a D ^ 4) ^ 2 := by ring
    _ ≤ (D : ℝ) ^ (8 * a - 1) * (∫ x, x ^ 4 ∂ν) ^ 2 := mul_le_mul_of_nonneg_left hsq hE
    _ = (∫ x, x ^ 4 ∂ν) ^ 2 * (D : ℝ) ^ (8 * a - 1) := by ring

set_option linter.unusedVariables false in
/-- L2. With probability tending to 1 the discarded part has operator norm at most `ε √d`:
no entry above `ε √d` (L2a) and no row or column with two nonzeros (a union bound at
`a < 1/8`), then L1.

The proof reads neither `hn` nor `ha`: the sparsity comes from the union bound at the level
`T`, and the entry bound from L2a at the level `ε √d`, so `a > 0` and `0 < n N` are slack
here. Both stay in the statement, which the review gate froze, so the unused-variable linter
is off for this declaration. -/
theorem tendsto_measure_discard_opNorm {ν : Measure ℝ} (hν : NoiseLaw ν) {n d : ℕ → ℕ} {c : ℝ}
    (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N) (hdtop : Tendsto d atTop atTop)
    (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha8 : a < 1 / 8) {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => (noiseMatrix ν (n N) (d N))
      {Z | ‖discardMat (truncLevel a (d N)) Z‖ ≤ ε * Real.sqrt (d N)}) atTop (𝓝 1) := by
  have hprob := hν.prob
  -- `d^(8a-1) → 0`, the place where `a < 1/8` enters
  have hEt : Tendsto (fun N => (d N : ℝ) ^ (8 * a - 1)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hdtop
    have h2 := (tendsto_rpow_neg_atTop (y := 1 - 8 * a) (by linarith)).comp h1
    simpa [Function.comp_def] using h2
  have hnu2 : Tendsto (fun N => (∫ x, x ^ 4 ∂ν) ^ 2 * (d N : ℝ) ^ (8 * a - 1)) atTop (𝓝 0) := by
    simpa using tendsto_const_nhds.mul hEt
  refine tendsto_measure_one_of_bad
    (μ := fun N => noiseMatrix ν (n N) (d N))
    (t := fun N => {Z : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      ‖discardMat (truncLevel a (d N)) Z‖ ≤ ε * Real.sqrt (d N)})
    (s := fun N => bigEntrySet (ε * Real.sqrt (d N)) (n N) (d N) ∪
      (rowPairSet (truncLevel a (d N)) (n N) (d N) ∪
        colPairSet (truncLevel a (d N)) (n N) (d N))) ?_ ?_
  · -- the good event contains the complement of the three bad events, by L1
    intro N Z hZ
    by_contra hns
    refine hZ ?_
    have hA : ∀ i j, |Z i j| ≤ ε * Real.sqrt (d N) := by
      intro i j
      by_contra hc
      exact hns (Or.inl (Set.mem_iUnion.2 ⟨(i, j), not_le.1 hc⟩))
    have hB : ∀ i j j', truncLevel a (d N) < |Z i j| → truncLevel a (d N) < |Z i j'| →
        j = j' := by
      intro i j j' h1 h2
      by_contra hjj
      exact hns (Or.inr (Or.inl (Set.mem_iUnion.2 ⟨(i, j, j'), h1, h2, hjj⟩)))
    have hC : ∀ i i' j, truncLevel a (d N) < |Z i j| → truncLevel a (d N) < |Z i' j| →
        i = i' := by
      intro i i' j h1 h2
      by_contra hii
      exact hns (Or.inr (Or.inr (Set.mem_iUnion.2 ⟨(i, i', j), h1, h2, hii⟩)))
    have hne : ∀ x : ℝ, discardMap (truncLevel a (d N)) x ≠ 0 →
        truncLevel a (d N) < |x| := by
      intro x hx
      by_contra hcc
      apply hx
      unfold discardMap
      rw [if_pos (not_lt.1 hcc)]
    refine opNorm_le_of_sparse _ (mul_nonneg hε.le (Real.sqrt_nonneg _))
      (fun i j j' h1 h2 => hB i j j' (hne _ h1) (hne _ h2))
      (fun i i' j h1 h2 => hC i i' j (hne _ h1) (hne _ h2)) fun i j => ?_
    change |discardMap (truncLevel a (d N)) (Z i j)| ≤ ε * Real.sqrt (d N)
    unfold discardMap
    by_cases hcc : |Z i j| ≤ truncLevel a (d N)
    · rw [if_pos hcc, abs_zero]
      exact mul_nonneg hε.le (Real.sqrt_nonneg _)
    · rw [if_neg hcc]
      exact hA i j
  · -- the three bad events are null tending
    refine tendsto_measure_zero_union ?_ (tendsto_measure_zero_union ?_ ?_)
    · exact tendsto_measure_zero_of_le_ofReal
        (fun N => measure_bigEntrySet_le (n N) (d N) (ε * Real.sqrt (d N)))
        (tendsto_entry_tail hν hdtop hcN hε)
    · refine tendsto_measure_zero_of_le_ofReal
        (g := fun N => (n N : ℝ) / d N * ((∫ x, x ^ 4 ∂ν) ^ 2 * (d N : ℝ) ^ (8 * a - 1)))
        (fun N => ?_) (by simpa using hcN.mul hnu2)
      refine le_trans (measure_rowPairSet_le (n N) (d N) _) (ENNReal.ofReal_le_ofReal ?_)
      have hDpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hM : (0 : ℝ) ≤ (n N : ℝ) := Nat.cast_nonneg _
      have hcube := cube_mul_tail_sq_le (a := a) hν (hd N)
      rw [div_mul_eq_mul_div, le_div_iff₀ hDpos]
      calc (n N : ℝ) * d N * d N * (ν {x : ℝ | truncLevel a (d N) < |x|}).toReal ^ 2 * d N
          = (n N : ℝ) * ((d N : ℝ) ^ 3
              * (ν {x : ℝ | truncLevel a (d N) < |x|}).toReal ^ 2) := by ring
        _ ≤ (n N : ℝ) * ((∫ x, x ^ 4 ∂ν) ^ 2 * (d N : ℝ) ^ (8 * a - 1)) :=
            mul_le_mul_of_nonneg_left hcube hM
    · refine tendsto_measure_zero_of_le_ofReal
        (g := fun N => ((n N : ℝ) / d N) ^ 2 * ((∫ x, x ^ 4 ∂ν) ^ 2 * (d N : ℝ) ^ (8 * a - 1)))
        (fun N => ?_) (by simpa using (hcN.pow 2).mul hnu2)
      refine le_trans (measure_colPairSet_le (n N) (d N) _) (ENNReal.ofReal_le_ofReal ?_)
      have hDpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hM : (0 : ℝ) ≤ (n N : ℝ) ^ 2 := by positivity
      have hcube := cube_mul_tail_sq_le (a := a) hν (hd N)
      rw [div_pow, div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
      calc (n N : ℝ) * n N * d N * (ν {x : ℝ | truncLevel a (d N) < |x|}).toReal ^ 2
              * (d N : ℝ) ^ 2
          = (n N : ℝ) ^ 2 * ((d N : ℝ) ^ 3
              * (ν {x : ℝ | truncLevel a (d N) < |x|}).toReal ^ 2) := by ring
        _ ≤ (n N : ℝ) ^ 2 * ((∫ x, x ^ 4 ∂ν) ^ 2 * (d N : ℝ) ^ (8 * a - 1)) :=
            mul_le_mul_of_nonneg_left hcube hM

end Edge
end StackedSVD
