/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Defs

/-! # Stage 3, unit E: the trace expansion

`trace ((YᵀY)^k)` as a sum over the closed walks of length `2k` on the complete bipartite
graph `Fin n × Fin d`, and its expectation under `noiseMatrix ρ n d` as a sum of products of
entry moments (`notes/stage3_edge.md`, route C; `notes/STAGE3_CAMPAIGN.md`). The closed form
at `k = 2` is the check of the expansion. -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace Edge

/-- The cyclic successor of a non-final index is the ordinary successor. -/
private theorem cycSucc_castSucc {p : ℕ} (t : Fin p) : cycSucc t.castSucc = t.succ := by
  have ht : t.1 + 1 < p + 1 := by omega
  simp [cycSucc, Fin.ext_iff, Nat.mod_eq_of_lt ht]

/-- The cyclic successor of the final index is `0`. -/
private theorem cycSucc_last (p : ℕ) : cycSucc (Fin.last p) = 0 := by
  simp [cycSucc]

/-- A sum over `(p+1)`-tuples splits into the first entry and the remaining `p` entries. -/
private theorem sum_cons_eq {d p : ℕ} (F : (Fin (p + 1) → Fin d) → ℝ) :
    ∑ v : Fin (p + 1) → Fin d, F v = ∑ c : Fin d, ∑ w : Fin p → Fin d, F (Fin.cons c w) := by
  have h := (Fintype.sum_equiv (Fin.consEquiv fun _ : Fin (p + 1) => Fin d)
    (fun x : Fin d × (Fin p → Fin d) => F (Fin.cons x.1 x.2)) F fun _ => rfl).symm
  rw [h, Fintype.sum_prod_type]

/-- The entry of a matrix power as a sum over the open walks of length `p+1`. -/
private theorem pow_apply_cons {d : ℕ} (M : Matrix (Fin d) (Fin d) ℝ) (p : ℕ) (a b : Fin d) :
    (M ^ (p + 1)) a b
      = ∑ w : Fin p → Fin d,
          (∏ t : Fin p, M ((Fin.cons a w : Fin (p + 1) → Fin d) t.castSucc)
              ((Fin.cons a w : Fin (p + 1) → Fin d) t.succ))
            * M ((Fin.cons a w : Fin (p + 1) → Fin d) (Fin.last p)) b := by
  induction p generalizing a b with
  | zero => simp
  | succ q ih =>
      have hstep : (M ^ (q + 1 + 1)) a b = ∑ c : Fin d, M a c * (M ^ (q + 1)) c b := by
        rw [pow_succ', Matrix.mul_apply]
      rw [hstep, sum_cons_eq]
      refine Finset.sum_congr rfl fun c _ => ?_
      rw [ih c b, Finset.mul_sum]
      refine Finset.sum_congr rfl fun w _ => ?_
      rw [Fin.prod_univ_succ]
      simp only [Fin.castSucc_zero, Fin.cons_zero, Fin.cons_succ, ← Fin.succ_castSucc,
        ← Fin.succ_last]
      ring

/-- The trace of a power of a square matrix as a sum over its closed walks. -/
private theorem trace_pow_eq_sum_cyc {d : ℕ} (M : Matrix (Fin d) (Fin d) ℝ) (p : ℕ) :
    Matrix.trace (M ^ (p + 1))
      = ∑ j : Fin (p + 1) → Fin d, ∏ t : Fin (p + 1), M (j t) (j (cycSucc t)) := by
  rw [sum_cons_eq, Matrix.trace]
  simp only [Matrix.diag_apply]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [pow_apply_cons M p a a]
  refine Finset.sum_congr rfl fun w _ => ?_
  rw [Fin.prod_univ_castSucc]
  congr 1
  · exact Finset.prod_congr rfl fun t _ => by rw [cycSucc_castSucc]
  · rw [cycSucc_last, Fin.cons_zero]

/-- E1, deterministic. The trace of a power of the Gram matrix as a sum over closed walks
on the complete bipartite graph. -/
theorem trace_gram_pow_eq {n d k : ℕ} (hk : 1 ≤ k) (Y : Matrix (Fin n) (Fin d) ℝ) :
    Matrix.trace ((Yᵀ * Y) ^ k)
      = ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t)) := by
  obtain ⟨p, rfl⟩ : ∃ p, k = p + 1 := ⟨k - 1, by omega⟩
  rw [trace_pow_eq_sum_cyc, Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  have h1 : ∀ t : Fin (p + 1), (Yᵀ * Y) (j t) (j (cycSucc t))
      = ∑ a : Fin n, Y a (j t) * Y a (j (cycSucc t)) := by
    intro t; rw [Matrix.mul_apply]; rfl
  simp only [h1]
  rw [Fintype.prod_sum]

/-- Every walk multiplicity is at most `2k`. -/
private theorem walkMult_le {n d k : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (e : Fin n × Fin d) : walkMult i j e ≤ 2 * k := by
  have h1 : ({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k)).card ≤ k := by
    simpa using Finset.card_filter_le (Finset.univ : Finset (Fin k)) _
  have h2 : ({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k)).card ≤ k := by
    simpa using Finset.card_filter_le (Finset.univ : Finset (Fin k)) _
  have : walkMult i j e
      = ({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k)).card
        + ({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k)).card := rfl
  omega

/-- One walk product, regrouped by the entry it uses. -/
private theorem prod_walk_eq {n d k : ℕ} (Y : Matrix (Fin n) (Fin d) ℝ)
    (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t)))
      = ∏ e : Fin n × Fin d, Y e.1 e.2 ^ walkMult i j e := by
  have key : ∀ g : Fin k → Fin d, (∏ t : Fin k, Y (i t) (g t))
      = ∏ e : Fin n × Fin d,
          Y e.1 e.2 ^ ({t : Fin k | i t = e.1 ∧ g t = e.2} : Finset (Fin k)).card := by
    intro g
    rw [← Finset.prod_fiberwise' (Finset.univ : Finset (Fin k)) (fun t => (i t, g t))
      (fun e : Fin n × Fin d => Y e.1 e.2)]
    refine Finset.prod_congr rfl fun e _ => ?_
    rw [Finset.prod_const]
    congr 1
    congr 1
    ext t
    simp [Prod.ext_iff]
  rw [Finset.prod_mul_distrib, key j, key (fun t => j (cycSucc t)), ← Finset.prod_mul_distrib]
  exact Finset.prod_congr rfl fun e _ => (pow_add _ _ _).symm


/-- The integral of a product of entry powers factors into the entry moments. -/
private theorem integral_prod_entry_pow {n d : ℕ} {ρ : Measure ℝ} [IsProbabilityMeasure ρ]
    (m : Fin n × Fin d → ℕ) :
    ∫ Y : Matrix (Fin n) (Fin d) ℝ, (∏ e : Fin n × Fin d, Y e.1 e.2 ^ m e) ∂(noiseMatrix ρ n d)
      = ∏ e : Fin n × Fin d, ∫ x, x ^ m e ∂ρ := by
  have hrow : ∀ Y : Matrix (Fin n) (Fin d) ℝ,
      (∏ e : Fin n × Fin d, Y e.1 e.2 ^ m e)
        = ∏ r : Fin n, ∏ c : Fin d, Y r c ^ m (r, c) := fun Y => Fintype.prod_prod_type _
  simp only [hrow]
  have h1 : ∫ Y : Matrix (Fin n) (Fin d) ℝ, (∏ r : Fin n, ∏ c : Fin d, Y r c ^ m (r, c))
        ∂(noiseMatrix ρ n d)
      = ∏ r : Fin n, ∫ row : Fin d → ℝ, (∏ c : Fin d, row c ^ m (r, c))
          ∂(Measure.pi fun _ : Fin d => ρ) :=
    MeasureTheory.integral_fintype_prod_eq_prod
      (fun (r : Fin n) (row : Fin d → ℝ) => ∏ c : Fin d, row c ^ m (r, c))
  have h2 : ∀ r : Fin n, ∫ row : Fin d → ℝ, (∏ c : Fin d, row c ^ m (r, c))
        ∂(Measure.pi fun _ : Fin d => ρ)
      = ∏ c : Fin d, ∫ x : ℝ, x ^ m (r, c) ∂ρ :=
    fun r => MeasureTheory.integral_fintype_prod_eq_prod (fun (c : Fin d) (x : ℝ) => x ^ m (r, c))
  rw [h1]
  simp only [h2]
  exact (Fintype.prod_prod_type (fun e : Fin n × Fin d => ∫ x : ℝ, x ^ m e ∂ρ)).symm

/-- A product of entry powers is integrable when every power up to `2k` is. -/
private theorem integrable_prod_entry_pow {n d k : ℕ} {ρ : Measure ℝ} [IsProbabilityMeasure ρ]
    (hmom : ∀ p : ℕ, p ≤ 2 * k → Integrable (fun x => x ^ p) ρ)
    (m : Fin n × Fin d → ℕ) (hm : ∀ e, m e ≤ 2 * k) :
    Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ => ∏ e : Fin n × Fin d, Y e.1 e.2 ^ m e)
      (noiseMatrix ρ n d) := by
  have hrow : (fun Y : Matrix (Fin n) (Fin d) ℝ => ∏ e : Fin n × Fin d, Y e.1 e.2 ^ m e)
      = fun Y : Matrix (Fin n) (Fin d) ℝ => ∏ r : Fin n, (∏ c : Fin d, Y r c ^ m (r, c)) := by
    funext Y; exact Fintype.prod_prod_type _
  rw [hrow]
  refine MeasureTheory.Integrable.fintype_prod (f := fun (r : Fin n) (row : Fin d → ℝ) =>
    ∏ c : Fin d, row c ^ m (r, c)) fun r => ?_
  exact MeasureTheory.Integrable.fintype_prod (f := fun (c : Fin d) (x : ℝ) => x ^ m (r, c))
    fun c => hmom _ (hm _)

/-- E2. The expectation of one walk factors into the entry moments of `ρ`.
`hmom` is the weakest hypothesis the proof needs: every power up to `2k` is integrable. -/
theorem integral_trace_gram_pow {n d k : ℕ} {ρ : Measure ℝ} [IsProbabilityMeasure ρ]
    (hmom : ∀ p : ℕ, p ≤ 2 * k → Integrable (fun x => x ^ p) ρ) (hk : 1 ≤ k) :
    ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(noiseMatrix ρ n d)
      = ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
          ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ρ := by
  have hfun : ∀ (i : Fin k → Fin n) (j : Fin k → Fin d),
      (fun Y : Matrix (Fin n) (Fin d) ℝ =>
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t)))
        = fun Y : Matrix (Fin n) (Fin d) ℝ => ∏ e : Fin n × Fin d, Y e.1 e.2 ^ walkMult i j e :=
    fun i j => funext fun Y => prod_walk_eq Y i j
  have hint : ∀ (i : Fin k → Fin n) (j : Fin k → Fin d),
      Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ =>
        ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))) (noiseMatrix ρ n d) := by
    intro i j
    rw [hfun i j]
    exact integrable_prod_entry_pow hmom _ (walkMult_le i j)
  calc ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(noiseMatrix ρ n d)
      = ∫ Y : Matrix (Fin n) (Fin d) ℝ, (∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))) ∂(noiseMatrix ρ n d) := by
        simp only [trace_gram_pow_eq hk]
    _ = ∑ i : Fin k → Fin n, ∫ Y : Matrix (Fin n) (Fin d) ℝ, (∑ j : Fin k → Fin d,
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))) ∂(noiseMatrix ρ n d) :=
        integral_finsetSum _ fun i _ => integrable_finsetSum _ fun j _ => hint i j
    _ = ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d, ∫ Y : Matrix (Fin n) (Fin d) ℝ,
          (∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))) ∂(noiseMatrix ρ n d) :=
        Finset.sum_congr rfl fun i _ => integral_finsetSum _ fun j _ => hint i j
    _ = ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
          ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ρ := by
        refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
        rw [hfun i j]
        exact integral_prod_entry_pow _

/-- The cyclic successor of `0` in `Fin 2`. -/
private theorem cycSucc_zero_two : cycSucc (0 : Fin 2) = 1 := rfl

/-- The cyclic successor of `1` in `Fin 2`. -/
private theorem cycSucc_one_two : cycSucc (1 : Fin 2) = 0 := rfl

/-- The walk multiplicity at `k = 2` as a sum of four indicators. -/
private theorem walkMult_two {n d : ℕ} (i : Fin 2 → Fin n) (j : Fin 2 → Fin d)
    (e : Fin n × Fin d) :
    walkMult i j e
      = ((if i 0 = e.1 ∧ j 0 = e.2 then 1 else 0)
          + (if i 1 = e.1 ∧ j 1 = e.2 then 1 else 0))
        + ((if i 0 = e.1 ∧ j 1 = e.2 then 1 else 0)
          + (if i 1 = e.1 ∧ j 0 = e.2 then 1 else 0)) := by
  simp only [walkMult, Finset.card_filter, Fin.sum_univ_two, cycSucc_zero_two, cycSucc_one_two]
  rfl

/-- A sum over pairs of an equality test. -/
private theorem sum_fin_two_ite {m : ℕ} (A B : ℝ) :
    (∑ v : Fin 2 → Fin m, if v 0 = v 1 then A else B)
      = (m : ℝ) * A + ((m : ℝ) * m - m) * B := by
  have he : (∑ v : Fin 2 → Fin m, if v 0 = v 1 then A else B)
      = ∑ p : Fin m × Fin m, if p.1 = p.2 then A else B :=
    Fintype.sum_equiv (finTwoArrowEquiv (Fin m)) _ _ fun v => by simp
  rw [he, Fintype.sum_prod_type]
  have hin : ∀ a : Fin m, (∑ b : Fin m, if a = b then A else B) = (A - B) + (m : ℝ) * B := by
    intro a
    have hb : ∀ b : Fin m, (if a = b then A else B) = (if a = b then A - B else 0) + B := by
      intro b; by_cases h : a = b <;> simp [h]
    rw [Finset.sum_congr rfl fun b _ => hb b, Finset.sum_add_distrib, Finset.sum_ite_eq,
      Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    simp
  rw [Finset.sum_congr rfl fun a _ => hin a, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul]
  ring

/-- E3. The closed form at `k = 2`, the check of E2: `∫ trace ((ZᵀZ)²) = n d (n + d + ν₄ - 2)`. -/
theorem integral_trace_gram_sq {ν : Measure ℝ} (hν : NoiseLaw ν) (n d : ℕ) :
    ∫ Y, Matrix.trace ((Yᵀ * Y) ^ 2) ∂(noiseMatrix ν n d)
      = (n : ℝ) * d * ((n : ℝ) + d + (∫ x, x ^ 4 ∂ν) - 2) := by
  have hprob := hν.prob
  rw [integral_trace_gram_pow (n := n) (d := d) (k := 2) (ρ := ν)
    (fun p hp => hν.integrable_pow (by omega)) (by norm_num)]
  have hzero : (∫ x, x ^ (0 : ℕ) ∂ν) = 1 := by simp
  have hone : (∫ x, x ^ (1 : ℕ) ∂ν) = 0 := by simpa using hν.mean
  have htwo : (∫ x, x ^ (2 : ℕ) ∂ν) = 1 := hν.var
  have hval : ∀ (i : Fin 2 → Fin n) (j : Fin 2 → Fin d),
      (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ν)
        = if i 0 = i 1 then (if j 0 = j 1 then (∫ x, x ^ 4 ∂ν) else 1)
            else (if j 0 = j 1 then 1 else 0) := by
    intro i j
    by_cases hi : i 0 = i 1 <;> by_cases hj : j 0 = j 1
    · rw [if_pos hi, if_pos hj]
      have hm : ∀ e : Fin n × Fin d, walkMult i j e = if e = (i 0, j 0) then 4 else 0 := by
        intro e
        rw [walkMult_two, ← hi, ← hj]
        by_cases h : e = (i 0, j 0)
        · subst h; simp
        · have h' : ¬(i 0 = e.1 ∧ j 0 = e.2) := by
            rintro ⟨h1, h2⟩
            exact h (Prod.ext h1.symm h2.symm)
          simp [h', h]
      have hsingle : (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂ν)
          = ∫ x, x ^ (walkMult i j (i 0, j 0)) ∂ν :=
        Finset.prod_eq_single_of_mem _ (Finset.mem_univ _) fun b _ hb => by
          rw [hm b, if_neg hb]; simp
      rw [hsingle, hm (i 0, j 0), if_pos rfl]
    · rw [if_pos hi, if_neg hj]
      refine Finset.prod_eq_one fun e _ => ?_
      have hm : walkMult i j e = 0 ∨ walkMult i j e = 2 := by
        rw [walkMult_two, ← hi]
        by_cases hie : i 0 = e.1
        · by_cases h0 : j 0 = e.2
          · have h1 : ¬j 1 = e.2 := fun h => hj (h0.trans h.symm)
            right; simp [hie, h0, h1]
          · by_cases h1 : j 1 = e.2
            · right; simp [hie, h0, h1]
            · left; simp [hie, h0, h1]
        · left; simp [hie]
      rcases hm with h | h <;> rw [h]
      · simp
      · exact htwo
    · rw [if_neg hi, if_pos hj]
      refine Finset.prod_eq_one fun e _ => ?_
      have hm : walkMult i j e = 0 ∨ walkMult i j e = 2 := by
        rw [walkMult_two, ← hj]
        by_cases hje : j 0 = e.2
        · by_cases h0 : i 0 = e.1
          · have h1 : ¬i 1 = e.1 := fun h => hi (h0.trans h.symm)
            right; simp [hje, h0, h1]
          · by_cases h1 : i 1 = e.1
            · right; simp [hje, h0, h1]
            · left; simp [hje, h0, h1]
        · left; simp [hje]
      rcases hm with h | h <;> rw [h]
      · simp
      · exact htwo
    · rw [if_neg hi, if_neg hj]
      refine Finset.prod_eq_zero (Finset.mem_univ (i 0, j 0)) ?_
      have hm : walkMult i j (i 0, j 0) = 1 := by
        have c1 : i 0 = ((i 0, j 0) : Fin n × Fin d).1 ∧ j 0 = ((i 0, j 0) : Fin n × Fin d).2 :=
          ⟨rfl, rfl⟩
        have c2 : ¬(i 1 = ((i 0, j 0) : Fin n × Fin d).1
            ∧ j 1 = ((i 0, j 0) : Fin n × Fin d).2) := fun h => hi h.1.symm
        have c3 : ¬(i 0 = ((i 0, j 0) : Fin n × Fin d).1
            ∧ j 1 = ((i 0, j 0) : Fin n × Fin d).2) := fun h => hj h.2.symm
        have c4 : ¬(i 1 = ((i 0, j 0) : Fin n × Fin d).1
            ∧ j 0 = ((i 0, j 0) : Fin n × Fin d).2) := fun h => hi h.1.symm
        rw [walkMult_two, if_pos c1, if_neg c2, if_neg c3, if_neg c4]
      rw [hm]
      simpa using hone
  simp only [hval]
  have hinner : ∀ i : Fin 2 → Fin n,
      (∑ j : Fin 2 → Fin d, if i 0 = i 1 then (if j 0 = j 1 then (∫ x, x ^ 4 ∂ν) else 1)
          else (if j 0 = j 1 then 1 else 0))
        = if i 0 = i 1 then ((d : ℝ) * (∫ x, x ^ 4 ∂ν) + ((d : ℝ) * d - d) * 1)
            else ((d : ℝ) * 1 + ((d : ℝ) * d - d) * 0) := by
    intro i
    by_cases h : i 0 = i 1
    · simp only [h, if_true]
      exact sum_fin_two_ite _ _
    · simp only [h, if_false]
      exact sum_fin_two_ite _ _
  rw [Finset.sum_congr rfl fun i _ => hinner i, sum_fin_two_ite]
  ring

/-- E4. The trace moment is nonnegative at every matrix: `(YᵀY)^k` is positive semidefinite. -/
theorem trace_gram_pow_nonneg {n d k : ℕ} (Y : Matrix (Fin n) (Fin d) ℝ) :
    0 ≤ Matrix.trace ((Yᵀ * Y) ^ k) := by
  have h : (Yᵀ * Y).PosSemidef := by
    have h0 := Matrix.posSemidef_conjTranspose_mul_self Y
    rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h0
  exact (h.pow k).trace_nonneg

/-- E6. The trace moment is a measurable function of the matrix. Unit K8 needs it to name
the bad event. -/
theorem measurable_trace_gram_pow (n d k : ℕ) :
    Measurable (fun Y : Matrix (Fin n) (Fin d) ℝ => Matrix.trace ((Yᵀ * Y) ^ k)) := by
  have hentry : ∀ a b : Fin d,
      Measurable (fun Y : Matrix (Fin n) (Fin d) ℝ => ((Yᵀ * Y) ^ k) a b) := by
    induction k with
    | zero => intro a b; simp only [pow_zero]; fun_prop
    | succ m ih =>
        intro a b
        have hs : ∀ Y : Matrix (Fin n) (Fin d) ℝ,
            ((Yᵀ * Y) ^ (m + 1)) a b = ∑ c : Fin d, ((Yᵀ * Y) ^ m) a c * (Yᵀ * Y) c b := by
          intro Y; rw [pow_succ, Matrix.mul_apply]
        simp only [hs]
        refine Finset.measurable_sum _ fun c _ => (ih a c).mul ?_
        have he : ∀ Y : Matrix (Fin n) (Fin d) ℝ, (Yᵀ * Y) c b = ∑ r : Fin n, Y r c * Y r b := by
          intro Y; rw [Matrix.mul_apply]; rfl
        simp only [he]
        exact Finset.measurable_sum _ fun r _ => by fun_prop
  simp only [Matrix.trace, Matrix.diag_apply]
  exact Finset.measurable_sum _ fun a _ => hentry a a

/-- Under a bounded law every entry of the noise matrix is bounded almost everywhere. -/
private theorem ae_entry_le {n d : ℕ} {ρ : Measure ℝ} {K : ℝ} (hρ : TruncNoiseLaw ρ K) :
    ∀ᵐ Y ∂(noiseMatrix ρ n d), ∀ (i : Fin n) (j : Fin d), |Y i j| ≤ K := by
  have := hρ.prob
  rw [MeasureTheory.ae_all_iff]
  intro i
  rw [MeasureTheory.ae_all_iff]
  intro j
  have h1 : MeasurePreserving (fun Y : Matrix (Fin n) (Fin d) ℝ => Y i)
      (noiseMatrix ρ n d) (Measure.pi fun _ : Fin d => ρ) :=
    measurePreserving_eval (fun _ : Fin n => Measure.pi fun _ : Fin d => ρ) i
  have h2 : MeasurePreserving (fun row : Fin d → ℝ => row j)
      (Measure.pi fun _ : Fin d => ρ) ρ :=
    measurePreserving_eval (fun _ : Fin d => ρ) j
  exact (h2.comp h1).quasiMeasurePreserving.ae hρ.bdd

/-- E5. The trace moment is integrable under a bounded law. Unit K8 needs it for Markov. -/
theorem integrable_trace_gram_pow {n d k : ℕ} {ρ : Measure ℝ} {K : ℝ} (hK : 0 ≤ K)
    (hρ : TruncNoiseLaw ρ K) :
    Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ => Matrix.trace ((Yᵀ * Y) ^ k))
      (noiseMatrix ρ n d) := by
  have := hρ.prob
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp only [pow_zero, Matrix.trace_one]
    exact integrable_const _
  refine Integrable.mono' (g := fun _ : Matrix (Fin n) (Fin d) ℝ =>
      (n : ℝ) ^ k * (d : ℝ) ^ k * (K ^ 2) ^ k) (integrable_const _)
    (measurable_trace_gram_pow n d k).aestronglyMeasurable ?_
  filter_upwards [ae_entry_le (n := n) (d := d) hρ] with Y hY
  rw [Real.norm_eq_abs, trace_gram_pow_eq hk Y]
  have hterm : ∀ (i : Fin k → Fin n) (j : Fin k → Fin d),
      |∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))| ≤ (K ^ 2) ^ k := by
    intro i j
    rw [Finset.abs_prod]
    calc ∏ t : Fin k, |Y (i t) (j t) * Y (i t) (j (cycSucc t))|
        ≤ ∏ _t : Fin k, K ^ 2 := by
          refine Finset.prod_le_prod (fun t _ => abs_nonneg _) fun t _ => ?_
          rw [abs_mul, sq]
          exact mul_le_mul (hY _ _) (hY _ _) (abs_nonneg _) hK
      _ = (K ^ 2) ^ k := by simp
  calc |∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))|
      ≤ ∑ i : Fin k → Fin n, |∑ j : Fin k → Fin d,
          ∏ t : Fin k, Y (i t) (j t) * Y (i t) (j (cycSucc t))| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i : Fin k → Fin n, ∑ _j : Fin k → Fin d, (K ^ 2) ^ k := by
        refine Finset.sum_le_sum fun i _ => ?_
        exact (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun j _ => hterm i j)
    _ = (n : ℝ) ^ k * (d : ℝ) ^ k * (K ^ 2) ^ k := by
        simp [Finset.sum_const, mul_assoc]

end Edge
end StackedSVD
