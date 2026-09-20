/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Defs
import StackedSVD.RMT.R3

/-! # Stage 3, unit T: the truncation

The one-entry facts about the truncated and recentered law `truncLaw ν T` and the split
`Z = Ẑ + R + m_T J` (`notes/stage3_edge.md`, route C; the unit table is in
`notes/STAGE3_CAMPAIGN.md`). -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace Edge

/-! ### The kept part of one entry

The helper lemmas below are about `x ↦ x 1{|x| ≤ T}`, the truncation before the recentering.
`truncMean ν T` is its integral and `truncMap ν T` is it minus that integral. -/

/-- The kept part of one entry, `x 1{|x| ≤ T}`, is measurable. -/
private theorem measurable_keep (T : ℝ) :
    Measurable (fun x : ℝ => if |x| ≤ T then x else 0) :=
  Measurable.ite (measurableSet_le continuous_abs.measurable measurable_const) measurable_id
    measurable_const

/-- The kept part of one entry is at most `T` in absolute value. -/
private theorem abs_keep_le {T : ℝ} (hT : 0 ≤ T) (x : ℝ) :
    |if |x| ≤ T then x else 0| ≤ T := by
  split_ifs with h
  · exact h
  · simpa using hT

/-- The kept part of one entry is integrable against a probability measure. -/
private theorem integrable_keep (ν : Measure ℝ) [IsProbabilityMeasure ν] {T : ℝ} (hT : 0 ≤ T) :
    Integrable (fun x : ℝ => if |x| ≤ T then x else 0) ν :=
  Integrable.mono' (integrable_const T) (measurable_keep T).aestronglyMeasurable
    (Eventually.of_forall fun x => by
      rw [Real.norm_eq_abs]; exact abs_keep_le hT x)

/-- The square of the kept part of one entry is integrable against a probability measure. -/
private theorem integrable_keep_sq (ν : Measure ℝ) [IsProbabilityMeasure ν] {T : ℝ}
    (hT : 0 ≤ T) : Integrable (fun x : ℝ => (if |x| ≤ T then x else 0) ^ 2) ν :=
  Integrable.mono' (integrable_const (T ^ 2))
    ((measurable_keep T).pow_const 2).aestronglyMeasurable
    (Eventually.of_forall fun x => by
      rw [Real.norm_eq_abs, abs_pow]
      exact pow_le_pow_left₀ (abs_nonneg _) (abs_keep_le hT x) 2)

/-- The mean the truncation removes is at most the level `T` in absolute value. -/
private theorem abs_truncMean_le_level (ν : Measure ℝ) [IsProbabilityMeasure ν] {T : ℝ}
    (hT : 0 ≤ T) : |truncMean ν T| ≤ T := by
  have h1 : |truncMean ν T| ≤ ∫ x, |if |x| ≤ T then x else 0| ∂ν :=
    abs_integral_le_integral_abs
  refine h1.trans ?_
  have h2 : ∫ x, |if |x| ≤ T then x else 0| ∂ν ≤ ∫ _x : ℝ, T ∂ν :=
    integral_mono_ae (integrable_keep ν hT).abs (integrable_const T)
      (Eventually.of_forall fun x => abs_keep_le hT x)
  simpa using h2

/-! ### T0 to T4: the truncated law -/

/-- T0. The truncated and recentered entry map is measurable. -/
theorem measurable_truncMap (ν : Measure ℝ) (T : ℝ) : Measurable (truncMap ν T) :=
  (measurable_keep T).sub_const _

/-- T1. `truncLaw` is a probability measure. -/
instance isProbabilityMeasure_truncLaw (ν : Measure ℝ) [IsProbabilityMeasure ν] (T : ℝ) :
    IsProbabilityMeasure (truncLaw ν T) :=
  Measure.isProbabilityMeasure_map (measurable_truncMap ν T).aemeasurable

/-- T2 to T4, bundled. The truncated law is centered, has variance at most 1, and is
supported in `[-2T, 2T]`. -/
theorem truncNoiseLaw_truncLaw {ν : Measure ℝ} (hν : NoiseLaw ν) {T : ℝ} (hT : 0 < T) :
    TruncNoiseLaw (truncLaw ν T) (2 * T) := by
  have hprob := hν.prob
  have hmap : AEMeasurable (truncMap ν T) ν := (measurable_truncMap ν T).aemeasurable
  have hkeep : Integrable (fun x : ℝ => if |x| ≤ T then x else 0) ν := integrable_keep ν hT.le
  have hkeep2 : Integrable (fun x : ℝ => (if |x| ≤ T then x else 0) ^ 2) ν :=
    integrable_keep_sq ν hT.le
  have hkm : ∫ x, (if |x| ≤ T then x else 0) ∂ν = truncMean ν T := rfl
  -- T2: the truncated law is centered.
  have hmean : ∫ x, x ∂(truncLaw ν T) = 0 := by
    have h1 : ∫ x, x ∂(truncLaw ν T) = ∫ x, truncMap ν T x ∂ν :=
      integral_map hmap (measurable_id.aestronglyMeasurable)
    rw [h1]
    simp only [truncMap]
    rw [integral_sub hkeep (integrable_const _), integral_const, hkm]
    simp
  -- T3: the truncated law has variance at most one.
  have hvar : ∫ x, x ^ 2 ∂(truncLaw ν T) ≤ 1 := by
    have h1 : ∫ x, x ^ 2 ∂(truncLaw ν T) = ∫ x, (truncMap ν T x) ^ 2 ∂ν :=
      integral_map hmap ((measurable_id.pow_const 2).aestronglyMeasurable)
    have hexp : ∀ x : ℝ, (truncMap ν T x) ^ 2
        = ((if |x| ≤ T then x else 0) ^ 2 - 2 * truncMean ν T * (if |x| ≤ T then x else 0))
          + truncMean ν T ^ 2 := by
      intro x; simp only [truncMap]; ring
    have hmul : Integrable
        (fun x : ℝ => 2 * truncMean ν T * (if |x| ≤ T then x else 0)) ν :=
      hkeep.const_mul (2 * truncMean ν T)
    have hsub : Integrable (fun x : ℝ => (if |x| ≤ T then x else 0) ^ 2
        - 2 * truncMean ν T * (if |x| ≤ T then x else 0)) ν := hkeep2.sub hmul
    rw [h1]
    simp only [hexp]
    rw [integral_add hsub (integrable_const _), integral_sub hkeep2 hmul,
      integral_const_mul, integral_const, hkm]
    have hle : ∫ x, (if |x| ≤ T then x else 0) ^ 2 ∂ν ≤ 1 := by
      have h2 : Integrable (fun x : ℝ => x ^ 2) ν := hν.integrable_pow (by norm_num)
      have h3 := integral_mono_ae hkeep2 h2 (Eventually.of_forall fun x => by
        change (if |x| ≤ T then x else 0) ^ 2 ≤ x ^ 2
        split_ifs with h
        · exact le_rfl
        · simpa using sq_nonneg x)
      rwa [hν.var] at h3
    simp only [probReal_univ, one_smul]
    nlinarith [sq_nonneg (truncMean ν T)]
  -- T4: the truncated law is supported in `[-2T, 2T]`.
  have hbdd : ∀ᵐ x ∂(truncLaw ν T), |x| ≤ 2 * T := by
    have hset : MeasurableSet {y : ℝ | |y| ≤ 2 * T} :=
      measurableSet_le continuous_abs.measurable measurable_const
    have hb : ∀ x : ℝ, |truncMap ν T x| ≤ 2 * T := by
      intro x
      have hsplit : |truncMap ν T x|
          ≤ |if |x| ≤ T then x else 0| + |truncMean ν T| := by
        simp only [truncMap, sub_eq_add_neg]
        simpa [abs_neg] using abs_add_le (if |x| ≤ T then x else 0) (-truncMean ν T)
      have := add_le_add (abs_keep_le hT.le x) (abs_truncMean_le_level ν hT.le)
      linarith
    exact (ae_map_iff hmap hset).mpr (Eventually.of_forall hb)
  exact { prob := inferInstance, mean := hmean, var_le := hvar, bdd := hbdd }

/-- T5. The absolute moments of a bounded, centered law of variance at most 1:
`∫ |x|^p ∂ρ ≤ K^(p-2) ∫ x² ∂ρ ≤ K^(p-2)` for `p ≥ 2`. Unit X consumes this. -/
theorem integral_abs_pow_le {ρ : Measure ℝ} {K : ℝ} (hK : 0 ≤ K) (hρ : TruncNoiseLaw ρ K)
    {p : ℕ} (hp : 2 ≤ p) :
    ∫ x, |x| ^ p ∂ρ ≤ K ^ (p - 2) * ∫ x, x ^ 2 ∂ρ := by
  have hprob := hρ.prob
  have hKp : (0 : ℝ) ≤ K ^ (p - 2) := pow_nonneg hK _
  have hsq : Integrable (fun x : ℝ => x ^ 2) ρ := by
    refine Integrable.mono' (integrable_const (K ^ 2))
      ((measurable_id.pow_const 2).aestronglyMeasurable) ?_
    filter_upwards [hρ.bdd] with x hx
    rw [Real.norm_eq_abs, abs_pow]
    exact pow_le_pow_left₀ (abs_nonneg x) hx 2
  have hkey : ∀ᵐ x ∂ρ, |x| ^ p ≤ K ^ (p - 2) * x ^ 2 := by
    filter_upwards [hρ.bdd] with x hx
    have hsplit : |x| ^ p = |x| ^ (p - 2) * |x| ^ 2 := by
      rw [← pow_add]
      congr 1
      omega
    rw [hsplit, sq_abs]
    exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (abs_nonneg x) hx _) (sq_nonneg x)
  have habs : Integrable (fun x : ℝ => |x| ^ p) ρ := by
    refine Integrable.mono' (hsq.const_mul (K ^ (p - 2)))
      (((continuous_abs.measurable.comp measurable_id).pow_const p).aestronglyMeasurable) ?_
    filter_upwards [hkey] with x hx
    rw [Real.norm_eq_abs, abs_pow, abs_abs]
    exact hx
  calc ∫ x, |x| ^ p ∂ρ ≤ ∫ x, K ^ (p - 2) * x ^ 2 ∂ρ :=
        integral_mono_ae habs (hsq.const_mul _) hkey
    _ = K ^ (p - 2) * ∫ x, x ^ 2 ∂ρ := integral_const_mul _ _

/-- T6. The mean the truncation removes is at most `ν₄ T^(-3)`. -/
theorem abs_truncMean_le {ν : Measure ℝ} (hν : NoiseLaw ν) {T : ℝ} (hT : 0 < T) :
    |truncMean ν T| ≤ (∫ x, x ^ 4 ∂ν) / T ^ 3 := by
  have hprob := hν.prob
  have hid : Integrable (fun x : ℝ => x) ν := by
    simpa using hν.integrable_pow (k := 1) (by norm_num)
  have hkeep : Integrable (fun x : ℝ => if |x| ≤ T then x else 0) ν := integrable_keep ν hT.le
  have hdis : Integrable (fun x : ℝ => x - (if |x| ≤ T then x else 0)) ν := hid.sub hkeep
  have hval : ∫ x, (x - (if |x| ≤ T then x else 0)) ∂ν = -truncMean ν T := by
    rw [integral_sub hid hkeep, hν.mean]
    simp [truncMean]
  have h1 : |truncMean ν T| = |∫ x, (x - (if |x| ≤ T then x else 0)) ∂ν| := by
    rw [hval, abs_neg]
  have hbnd : ∀ᵐ x ∂ν, |x - (if |x| ≤ T then x else 0)| ≤ x ^ 4 / T ^ 3 := by
    refine Eventually.of_forall fun x => ?_
    have hT3 : (0 : ℝ) < T ^ 3 := by positivity
    split_ifs with h
    · simp only [sub_self, abs_zero]
      positivity
    · simp only [sub_zero]
      rw [le_div_iff₀ hT3]
      have h4 : x ^ 4 = |x| ^ 4 := (Even.pow_abs (by norm_num) x).symm
      have h3 : T ^ 3 ≤ |x| ^ 3 := pow_le_pow_left₀ hT.le (not_le.mp h).le 3
      have hpow : |x| ^ 4 = |x| * |x| ^ 3 := by ring
      rw [h4, hpow]
      exact mul_le_mul_of_nonneg_left h3 (abs_nonneg x)
  have hfin : ∫ x, x ^ 4 / T ^ 3 ∂ν = (∫ x, x ^ 4 ∂ν) / T ^ 3 := integral_div _ _
  rw [h1, ← hfin]
  exact abs_integral_le_integral_abs.trans
    (integral_mono_ae hdis.abs (hν.mom4.div_const _) hbnd)

/-! ### T7 to T8: the matrix split and the law of the truncated matrix -/

/-- T7. The deterministic split `Z = Ẑ + R + m_T J`. -/
theorem mat_split (ν : Measure ℝ) (T : ℝ) {n d : ℕ} (Z : Matrix (Fin n) (Fin d) ℝ) :
    Z = truncMat ν T Z + discardMat T Z + meanMat ν T n d := by
  ext i j
  simp only [Matrix.add_apply, truncMat, discardMat, meanMat, Matrix.of_apply, truncMap,
    discardMap]
  split_ifs with h <;> ring

/-- T7b. The entrywise truncation of a matrix is measurable. -/
theorem measurable_truncMat (ν : Measure ℝ) (T : ℝ) (n d : ℕ) :
    Measurable (truncMat ν T : Matrix (Fin n) (Fin d) ℝ → Matrix (Fin n) (Fin d) ℝ) := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  exact (measurable_truncMap ν T).comp ((measurable_pi_apply j).comp (measurable_pi_apply i))

/-- T8. The law of the truncated matrix: the entrywise map transports `noiseMatrix`. -/
theorem map_truncMat (ν : Measure ℝ) [IsProbabilityMeasure ν] (T : ℝ) (n d : ℕ) :
    Measure.map (truncMat ν T) (noiseMatrix ν n d) = noiseMatrix (truncLaw ν T) n d := by
  have hlaw : IsProbabilityMeasure (Measure.map (truncMap ν T) ν) :=
    Measure.isProbabilityMeasure_map (measurable_truncMap ν T).aemeasurable
  have hrowmeas : Measurable (fun (r : Fin d → ℝ) (j : Fin d) => truncMap ν T (r j)) :=
    measurable_pi_lambda _ fun j => (measurable_truncMap ν T).comp (measurable_pi_apply j)
  have hrow : Measure.map (fun (r : Fin d → ℝ) (j : Fin d) => truncMap ν T (r j))
        (Measure.pi fun _ : Fin d => ν)
      = Measure.pi fun _ : Fin d => truncLaw ν T :=
    Measure.pi_map_pi (fun _ => (measurable_truncMap ν T).aemeasurable)
  have hrowprob : IsProbabilityMeasure
      (Measure.map (fun (r : Fin d → ℝ) (j : Fin d) => truncMap ν T (r j))
        (Measure.pi fun _ : Fin d => ν)) :=
    Measure.isProbabilityMeasure_map hrowmeas.aemeasurable
  have houter := Measure.pi_map_pi (ι := Fin n) (X := fun _ : Fin n => Fin d → ℝ)
    (Y := fun _ : Fin n => Fin d → ℝ)
    (μ := fun _ : Fin n => Measure.pi fun _ : Fin d => ν)
    (f := fun _ : Fin n => fun (r : Fin d → ℝ) (j : Fin d) => truncMap ν T (r j))
    (fun _ => hrowmeas.aemeasurable)
  rw [hrow] at houter
  exact houter

/-- T8b. The measure transport K8 uses: an event about `Ẑ` under `noiseMatrix ν` is the same
event about `Y` under `noiseMatrix (truncLaw ν T)`. Proved from T8. -/
theorem measure_truncMat_preimage (ν : Measure ℝ) [IsProbabilityMeasure ν] (T : ℝ) (n d : ℕ)
    {S : Set (Matrix (Fin n) (Fin d) ℝ)} (hS : MeasurableSet S) :
    (noiseMatrix ν n d) (truncMat ν T ⁻¹' S) = (noiseMatrix (truncLaw ν T) n d) S := by
  rw [← map_truncMat ν T n d,
    Measure.map_apply (measurable_truncMat ν T n d) hS]

/-! ### T9: the mean shift is negligible -/

/-- T9a. The operator norm of the constant matrix: `‖m_T J‖ = |m_T| √(n d)`, bounded here. -/
theorem opNorm_meanMat_le (ν : Measure ℝ) (T : ℝ) (n d : ℕ) :
    ‖meanMat ν T n d‖ ≤ |truncMean ν T| * Real.sqrt ((n : ℝ) * d) := by
  have hsum : ∑ i : Fin n, ∑ j : Fin d, (meanMat ν T n d) i j ^ 2
      = ((n : ℝ) * d) * truncMean ν T ^ 2 := by
    simp only [meanMat, Matrix.of_apply, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul]
    ring
  calc ‖meanMat ν T n d‖ ≤ Real.sqrt (∑ i, ∑ j, (meanMat ν T n d) i j ^ 2) :=
        R3.l2_opNorm_le_frobenius _
    _ = |truncMean ν T| * Real.sqrt ((n : ℝ) * d) := by
        rw [hsum, Real.sqrt_mul (by positivity), Real.sqrt_sq_eq_abs, mul_comm]

/-- One `N` of T9b: at the level `T = D^(1/2-a)` and `n ≤ B D` the mean shift is at most
`ν₄ √B D^(3a-1) √D`. -/
private theorem meanMat_opNorm_le_pow {ν : Measure ℝ} (hν : NoiseLaw ν) (a : ℝ)
    {n D : ℕ} (hD : 0 < D) {B : ℝ} (hB : 0 < B) (hnB : (n : ℝ) ≤ B * D) :
    ‖meanMat ν (truncLevel a D) n D‖
      ≤ (∫ x, x ^ 4 ∂ν) * Real.sqrt B * ((D : ℝ) ^ (3 * a - 1)) * Real.sqrt D := by
  have hDR : (0 : ℝ) < (D : ℝ) := by exact_mod_cast hD
  have hT : (0 : ℝ) < truncLevel a D := truncLevel_pos hD
  have hstep1 := opNorm_meanMat_le ν (truncLevel a D) n D
  have hstep2 := abs_truncMean_le hν hT
  have hsq : Real.sqrt ((n : ℝ) * D) ≤ Real.sqrt B * (D : ℝ) := by
    have hnd : (n : ℝ) * (D : ℝ) ≤ B * (D : ℝ) ^ 2 := by nlinarith [hDR.le]
    calc Real.sqrt ((n : ℝ) * D) ≤ Real.sqrt (B * (D : ℝ) ^ 2) := Real.sqrt_le_sqrt hnd
      _ = Real.sqrt B * (D : ℝ) := by rw [Real.sqrt_mul hB.le, Real.sqrt_sq hDR.le]
  have hcomb : ‖meanMat ν (truncLevel a D) n D‖
      ≤ ((∫ x, x ^ 4 ∂ν) / (truncLevel a D) ^ 3) * (Real.sqrt B * (D : ℝ)) :=
    hstep1.trans (mul_le_mul hstep2 hsq (Real.sqrt_nonneg _) (by positivity))
  refine hcomb.trans (le_of_eq ?_)
  have hTL : truncLevel a D = (D : ℝ) ^ ((1 : ℝ) / 2 - a) := rfl
  have hT3 : (truncLevel a D) ^ 3 = (D : ℝ) ^ ((3 : ℝ) / 2 - 3 * a) := by
    have he : ((1 : ℝ) / 2 - a) * ((3 : ℕ) : ℝ) = (3 : ℝ) / 2 - 3 * a := by push_cast; ring
    rw [hTL, ← Real.rpow_natCast ((D : ℝ) ^ ((1 : ℝ) / 2 - a)) 3, ← Real.rpow_mul hDR.le, he]
  have hP : (0 : ℝ) < (D : ℝ) ^ ((3 : ℝ) / 2 - 3 * a) := Real.rpow_pos_of_pos hDR _
  have hDeq : (D : ℝ)
      = (D : ℝ) ^ (3 * a - 1) * Real.sqrt D * ((D : ℝ) ^ ((3 : ℝ) / 2 - 3 * a)) := by
    rw [Real.sqrt_eq_rpow, ← Real.rpow_add hDR, ← Real.rpow_add hDR,
      show 3 * a - 1 + (1 : ℝ) / 2 + ((3 : ℝ) / 2 - 3 * a) = 1 by ring, Real.rpow_one]
  have halg : ∀ P Q R Db nb sb : ℝ, P ≠ 0 → Db = Q * R * P →
      nb / P * (sb * Db) = nb * sb * Q * R := by
    intro P Q R Db nb sb hP0 hDb
    subst hDb
    field_simp
  rw [hT3]
  exact halg _ _ _ _ _ _ (ne_of_gt hP) hDeq

/-- T9b. The mean shift is negligible against `√d` at `T = d^(1/2-a)` with `a < 1/6`. -/
theorem meanMat_opNorm_small {ν : Measure ℝ} (hν : NoiseLaw ν) {n d : ℕ → ℕ} {c : ℝ}
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha6 : a < 1 / 6) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ N in atTop,
      ‖meanMat ν (truncLevel a (d N)) (n N) (d N)‖ ≤ ε * Real.sqrt (d N) := by
  have hprob := hν.prob
  have hBpos : (0 : ℝ) < |c| + 1 := by positivity
  have hd1 : ∀ᶠ N in atTop, 1 ≤ d N := hdtop.eventually_ge_atTop 1
  have hnB : ∀ᶠ N in atTop, (n N : ℝ) ≤ (|c| + 1) * d N := by
    have h1 : ∀ᶠ N in atTop, (n N : ℝ) / d N < c + 1 :=
      hcN.eventually_lt_const (by linarith)
    filter_upwards [h1, hd1] with N hN hd
    have hdpos : (0 : ℝ) < (d N : ℝ) := by
      have hd0 : 0 < d N := hd
      exact_mod_cast hd0
    rw [div_lt_iff₀ hdpos] at hN
    nlinarith [le_abs_self c]
  have hkey : ∀ᶠ N in atTop,
      (∫ x, x ^ 4 ∂ν) * Real.sqrt (|c| + 1) * ((d N : ℝ) ^ (3 * a - 1)) ≤ ε := by
    have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hdtop
    have h1 : Tendsto (fun N => ((d N : ℝ)) ^ (3 * a - 1)) atTop (𝓝 0) := by
      have h2 : (0 : ℝ) < 1 - 3 * a := by linarith
      have h3 := (tendsto_rpow_neg_atTop h2).comp hdR
      simp only [Function.comp_def] at h3
      rwa [show -(1 - 3 * a) = 3 * a - 1 by ring] at h3
    have h0 : Tendsto
        (fun N => (∫ x, x ^ 4 ∂ν) * Real.sqrt (|c| + 1) * ((d N : ℝ) ^ (3 * a - 1)))
        atTop (𝓝 0) := by
      simpa using h1.const_mul ((∫ x, x ^ 4 ∂ν) * Real.sqrt (|c| + 1))
    exact (h0.eventually_lt_const hε).mono fun N hN => hN.le
  filter_upwards [hd1, hnB, hkey] with N hd hnb hk
  have hdpos : 0 < d N := hd
  exact (meanMat_opNorm_le_pow hν a hdpos hBpos hnb).trans
    (mul_le_mul_of_nonneg_right hk (Real.sqrt_nonneg _))

end Edge
end StackedSVD
