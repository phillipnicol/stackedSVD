/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Trace

/-! # Stage 3, unit C: the comparison with the Gaussian

The expected trace moment at a `TruncNoiseLaw` is at most the Gaussian one plus the weighted
sum over the excess walks (`notes/stage3_edge.md`, route C; `notes/STAGE3_CAMPAIGN.md`). The
walks split into three classes: `A` (every multiplicity `0` or `2`; weight `(∫x²∂ρ)^k ≤ 1`
at `ρ`, weight `1` at the Gaussian), `O` (some multiplicity `1`; weight `0` at every centered
law) and `B` (`IsExcess`; the Gaussian weight is nonnegative and drops). -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace Edge

/-- The zeroth moment of a probability law is `1`. -/
private theorem integral_pow_zero {μ : Measure ℝ} [IsProbabilityMeasure μ] :
    ∫ x : ℝ, x ^ 0 ∂μ = 1 := by simp

/-- A law with bounded support has every moment: `|x| ^ p ≤ K ^ p` almost everywhere. -/
private theorem TruncNoiseLaw.integrable_pow {ρ : Measure ℝ} {K : ℝ}
    (hρ : TruncNoiseLaw ρ K) (p : ℕ) : Integrable (fun x : ℝ => x ^ p) ρ := by
  have := hρ.prob
  refine Integrable.mono' (integrable_const (K ^ p))
    (measurable_id.pow_const p).aestronglyMeasurable ?_
  filter_upwards [hρ.bdd] with x hx
  rw [Real.norm_eq_abs, abs_pow]
  exact pow_le_pow_left₀ (abs_nonneg x) hx p

/-- Every moment of the standard Gaussian exists, so unit E2 applies to it. -/
private theorem integrable_pow_gaussian (p : ℕ) :
    Integrable (fun x : ℝ => x ^ p) (gaussianReal 0 1) := by
  have hm : MemLp (id : ℝ → ℝ) ((p : ℕ) : ENNReal) (gaussianReal 0 1) :=
    memLp_id_gaussianReal' _ (by simp)
  have h := hm.integrable_norm_pow'
  rw [← integrable_norm_iff (by fun_prop)]
  simpa [norm_pow] using h

/-- Every moment of the standard Gaussian is nonnegative: an even power is nonnegative
pointwise, an odd moment is `0` because the law is symmetric (`gaussianReal_map_neg`). -/
private theorem integral_pow_gaussian_nonneg (m : ℕ) :
    0 ≤ ∫ x, x ^ m ∂(gaussianReal 0 1) := by
  rcases Nat.even_or_odd m with hm | hm
  · exact integral_nonneg fun x => hm.pow_nonneg x
  · have hsymm : Measure.map (fun x : ℝ => -x) (gaussianReal 0 1) = gaussianReal 0 1 := by
      rw [gaussianReal_map_neg, neg_zero]
    have h : ∫ x, x ^ m ∂(gaussianReal 0 1) = ∫ x, (-x) ^ m ∂(gaussianReal 0 1) := by
      conv_lhs => rw [← hsymm]
      exact integral_map measurable_neg.aemeasurable (by fun_prop)
    simp only [hm.neg_pow, integral_neg] at h
    linarith

/-- The paired class and the excess class are disjoint: a paired walk has no multiplicity
above `2`. -/
private theorem not_isExcess_of_isPaired {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (hP : IsPaired i j) : ¬ IsExcess i j := by
  intro hE
  obtain ⟨e, he⟩ := hE.2
  rcases hP e with h | h <;> omega

/-- A signed moment is at most the matching absolute moment in absolute value. -/
private theorem abs_integral_pow_le {ρ : Measure ℝ} (m : ℕ) :
    |∫ x, x ^ m ∂ρ| ≤ ∫ x, |x| ^ m ∂ρ := by
  have h := norm_integral_le_integral_norm (μ := ρ) (fun x : ℝ => x ^ m)
  simpa [Real.norm_eq_abs, abs_pow] using h

set_option linter.unusedVariables false in
/-- C. Split the walks into `A` (every multiplicity `0` or `2`), `O` (some multiplicity `1`,
weight `0` at any centered law) and `B` (`IsExcess`). On `A` the weight is `(∫x²∂ρ)^k ≤ 1` and
the Gaussian weight is `1`; on `B` the Gaussian weight is nonnegative.

`hK : 0 ≤ K` is not needed: `hρ.bdd` gives `|x| ≤ K` almost everywhere, and `0 ≤ |x|` then
supplies the sign. The hypothesis stays because the plan note and every call site carry it;
the linter warning about it is switched off for this declaration only. -/
theorem trace_le_gaussian_add_excess {ρ : Measure ℝ} {K : ℝ} (hK : 0 ≤ K)
    (hρ : TruncNoiseLaw ρ K) {n d k : ℕ} (hk : 1 ≤ k) :
    ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(noiseMatrix ρ n d)
      ≤ (∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix n d)) + excessSum ρ n d k := by
  have := hρ.prob
  rw [integral_trace_gram_pow (fun p _ => hρ.integrable_pow p) hk,
    gaussianMatrix_eq_noiseMatrix,
    integral_trace_gram_pow (ρ := gaussianReal 0 1) (fun p _ => integrable_pow_gaussian p) hk]
  simp only [excessSum]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_le_sum fun i _ => ?_
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_le_sum fun j _ => ?_
  by_cases h1 : ∃ e : Fin n × Fin d, walkMult i j e = 1
  · -- Class `O`: one factor is the mean, which is `0` at both laws.
    obtain ⟨e₀, he₀⟩ := h1
    have hLz : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ρ) = 0 :=
      Finset.prod_eq_zero (Finset.mem_univ e₀) (by rw [he₀]; simpa using hρ.mean)
    have hGz : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂(gaussianReal 0 1)) = 0 :=
      Finset.prod_eq_zero (Finset.mem_univ e₀) (by rw [he₀]; simp)
    have hEx : ¬ IsExcess i j := fun h => h.1 e₀ he₀
    simp only [hLz, hGz, if_neg hEx, add_zero]
    exact le_rfl
  · have h1 : ∀ e : Fin n × Fin d, walkMult i j e ≠ 1 := fun e he => h1 ⟨e, he⟩
    by_cases hP : IsPaired i j
    · -- Class `A`: the Gaussian weight is `1` and the `ρ` weight is at most `1`.
      have hG1 : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂(gaussianReal 0 1)) = 1 := by
        refine Finset.prod_eq_one fun e _ => ?_
        rcases hP e with h | h
        · rw [h]; exact integral_pow_zero
        · rw [h]; exact noiseLaw_gaussian.var
      have hL1 : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ρ) ≤ 1 := by
        refine Finset.prod_le_one (fun e _ => ?_) fun e _ => ?_
        · rcases hP e with h | h
          · rw [h, integral_pow_zero]; exact zero_le_one
          · rw [h]; exact integral_nonneg fun x => by positivity
        · rcases hP e with h | h
          · rw [h, integral_pow_zero]
          · rw [h]; exact hρ.var_le
      simp only [hG1, if_neg (not_isExcess_of_isPaired hP), add_zero]
      exact hL1
    · -- Class `B`: the Gaussian weight is nonnegative and the `ρ` weight is at most the
      -- absolute weight, which is the `excessSum` summand.
      have hEx : IsExcess i j := by
        refine ⟨h1, ?_⟩
        simp only [IsPaired, not_forall, not_or] at hP
        obtain ⟨e, he0, he2⟩ := hP
        exact ⟨e, by have := h1 e; omega⟩
      simp only [if_pos hEx]
      have hGnn : (0 : ℝ)
          ≤ ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂(gaussianReal 0 1) :=
        Finset.prod_nonneg fun e _ => integral_pow_gaussian_nonneg _
      have hAbs : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ρ)
          ≤ ∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ := by
        refine le_trans (le_abs_self _) ?_
        rw [Finset.abs_prod]
        exact Finset.prod_le_prod (fun e _ => abs_nonneg _) fun e _ => abs_integral_pow_le _
      linarith

end Edge
end StackedSVD
