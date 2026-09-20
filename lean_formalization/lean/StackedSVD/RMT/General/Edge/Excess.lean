/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Compare
import StackedSVD.RMT.General.Edge.Trunc
import StackedSVD.RMT.General.Edge.CountBound

/-! # Stage 3, unit X: the excess count

The Furedi-Komlos count: the weighted sum over the excess walks (some entry of multiplicity
at least 3, none of multiplicity 1) is at most `2 fkFactor K n d k` times the Gaussian trace
moment, once `fkFactor K n d k ≤ 1/2` (`notes/stage3_edge.md`, route C;
`notes/STAGE3_CAMPAIGN.md`). The bound is relative to the tree-walk class, whose size the
Gaussian moment dominates; an absolute count loses `4^k` against `(1 + √c)^(2k)` and does
not close the assembly. -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace Edge

/-! ## The private helpers

Every helper of this unit lives in `ExcessAux`. `Count.lean` proves several of the same walk
facts under the same short names in `StackedSVD.Edge`, and a private declaration whose short
name matches an imported public one in the same namespace is a hard error, so the helpers are
one namespace deeper. `notes/INTERFACES.md` records which of them belong in `Count.lean`. -/

namespace ExcessAux

/-! ### The walk arithmetic -/

/-- X2. The `2k` factors of a closed walk distribute over the entries: the multiplicities sum
to `2k`. Each half is `Finset.card_eq_sum_card_fiberwise` for the map `t ↦ (i t, j t)`,
respectively `t ↦ (i t, j (cycSucc t))`. -/
private theorem sum_walkMult {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    ∑ e : Fin n × Fin d, walkMult i j e = 2 * k := by
  have key : ∀ g : Fin k → Fin d, ∑ e : Fin n × Fin d,
      (Finset.univ.filter (fun t : Fin k => i t = e.1 ∧ g t = e.2)).card = k := by
    intro g
    have h := Finset.card_eq_sum_card_fiberwise
      (f := fun t : Fin k => (i t, g t)) (s := (Finset.univ : Finset (Fin k)))
      (t := (Finset.univ : Finset (Fin n × Fin d))) (by intro x _; simp)
    simp only [Finset.card_univ, Fintype.card_fin] at h
    refine Eq.trans ?_ h.symm
    refine Finset.sum_congr rfl fun e _ => ?_
    congr 1
    ext t
    simp [Prod.ext_iff]
  have hsplit : ∀ e : Fin n × Fin d, walkMult i j e
      = (Finset.univ.filter (fun t : Fin k => i t = e.1 ∧ j t = e.2)).card
        + (Finset.univ.filter (fun t : Fin k => i t = e.1 ∧ j (cycSucc t) = e.2)).card :=
    fun _ => rfl
  simp only [hsplit]
  rw [Finset.sum_add_distrib, key j, key (fun t => j (cycSucc t))]
  omega

/-- Membership in the edge set is a nonzero multiplicity. -/
private theorem mem_walkEdges {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    {e : Fin n × Fin d} : e ∈ walkEdges i j ↔ walkMult i j e ≠ 0 := by
  simp [walkEdges]

/-- The multiplicities sum to `2k` over the edge set alone. -/
private theorem sum_walkMult_edges {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    ∑ e ∈ walkEdges i j, walkMult i j e = 2 * k := by
  rw [← sum_walkMult i j]
  refine Finset.sum_subset (Finset.subset_univ _) fun e _ he => ?_
  by_contra h
  exact he (mem_walkEdges.mpr h)

/-- With no entry of multiplicity 1 every edge carries multiplicity at least 2. -/
private theorem two_le_walkMult {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) {e : Fin n × Fin d} (he : e ∈ walkEdges i j) : 2 ≤ walkMult i j e := by
  have h0 := mem_walkEdges.mp he
  have h1 := h e
  omega

/-- Twice the edge count is at most the multiplicity sum, when no multiplicity is 1. -/
private theorem two_mul_edges_card_le {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) : 2 * (walkEdges i j).card ≤ 2 * k := by
  rw [← sum_walkMult_edges i j]
  calc 2 * (walkEdges i j).card = ∑ _e ∈ walkEdges i j, 2 := by
        rw [Finset.sum_const, smul_eq_mul, mul_comm]
    _ ≤ ∑ e ∈ walkEdges i j, walkMult i j e :=
        Finset.sum_le_sum fun e he => two_le_walkMult h he

/-- X4. A walk with no entry of multiplicity 1 has at most `k` edges. -/
private theorem edges_card_le {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) : (walkEdges i j).card ≤ k :=
  Nat.le_of_mul_le_mul_left (two_mul_edges_card_le h) (by norm_num)

/-- One edge of multiplicity at least 3 makes the edge count drop below `k`. -/
private theorem edges_card_lt {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) {e₀ : Fin n × Fin d} (he₀ : 3 ≤ walkMult i j e₀) :
    (walkEdges i j).card < k := by
  have hmem : e₀ ∈ walkEdges i j := mem_walkEdges.mpr (by omega)
  have hlt : 2 * (walkEdges i j).card < ∑ e ∈ walkEdges i j, walkMult i j e := by
    calc 2 * (walkEdges i j).card = ∑ _e ∈ walkEdges i j, 2 := by
          rw [Finset.sum_const, smul_eq_mul, mul_comm]
      _ < ∑ e ∈ walkEdges i j, walkMult i j e :=
          Finset.sum_lt_sum (fun e he => two_le_walkMult h he) ⟨e₀, hmem, by omega⟩
  rw [sum_walkMult_edges i j] at hlt
  omega

/-! ### The geometric tail -/

/-- X9. The tail of a geometric series with ratio at most `1/2` is at most twice the ratio. -/
private theorem geom_tail {f : ℝ} (hf0 : 0 ≤ f) (hf : f ≤ 1 / 2) (m : ℕ) :
    ∑ D ∈ Finset.Icc 1 m, f ^ D ≤ 2 * f := by
  have key : ∀ p : ℕ, ∑ D ∈ Finset.Icc 1 p, f ^ D ≤ 2 * f - 2 * f ^ (p + 1) := by
    intro p
    induction p with
    | zero => simp
    | succ q ih =>
      rw [Finset.sum_Icc_succ_top (by omega)]
      have hpow : (0 : ℝ) ≤ f ^ (q + 1) := pow_nonneg hf0 _
      have hstep : 2 * f ^ (q + 1 + 1) ≤ f ^ (q + 1) := by
        have hrw : f ^ (q + 1 + 1) = f ^ (q + 1) * f := by ring
        rw [hrw]
        nlinarith
      linarith
  have h1 := key m
  have h2 : (0 : ℝ) ≤ 2 * f ^ (m + 1) := by positivity
  linarith

/-! ### The first-visit injection -/

/-- The least index of a nonempty index set, `0` on the empty set. -/
private def firstIdx {k : ℕ} [NeZero k] (s : Finset (Fin k)) : Fin k :=
  if h : s.Nonempty then s.min' h else 0

/-- The least index belongs to a nonempty index set. -/
private theorem firstIdx_mem {k : ℕ} [NeZero k] {s : Finset (Fin k)} (h : s.Nonempty) :
    firstIdx s ∈ s := by
  rw [firstIdx, dif_pos h]
  exact s.min'_mem h

/-- The least index is at most every index of the set. -/
private theorem firstIdx_le {k : ℕ} [NeZero k] {s : Finset (Fin k)} {t : Fin k} (ht : t ∈ s) :
    firstIdx s ≤ t := by
  rw [firstIdx, dif_pos ⟨t, ht⟩]
  exact s.min'_le t ht

/-- The first time the walk visits the row `r`. -/
private def rowTime {k n : ℕ} [NeZero k] (i : Fin k → Fin n) (r : Fin n) : Fin k :=
  firstIdx (Finset.univ.filter (fun t : Fin k => i t = r))

/-- The first time the walk visits the column `c`. -/
private def colTime {k d : ℕ} [NeZero k] (j : Fin k → Fin d) (c : Fin d) : Fin k :=
  firstIdx (Finset.univ.filter (fun t : Fin k => j t = c))

/-- The walk is at the row `r` at the time `rowTime i r`, once `r` is visited at all. -/
private theorem rowTime_spec {k n : ℕ} [NeZero k] {i : Fin k → Fin n} {r : Fin n} {t : Fin k}
    (ht : i t = r) : i (rowTime i r) = r := by
  have h : (Finset.univ.filter (fun t : Fin k => i t = r)).Nonempty := ⟨t, by simp [ht]⟩
  have hm := firstIdx_mem h
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hm
  exact hm

/-- No earlier time visits the row `r`. -/
private theorem rowTime_le {k n : ℕ} [NeZero k] {i : Fin k → Fin n} {r : Fin n} {t : Fin k}
    (ht : i t = r) : rowTime i r ≤ t :=
  firstIdx_le (by simp [ht])

/-- The walk is at the column `c` at the time `colTime j c`, once `c` is visited at all. -/
private theorem colTime_spec {k d : ℕ} [NeZero k] {j : Fin k → Fin d} {c : Fin d} {t : Fin k}
    (ht : j t = c) : j (colTime j c) = c := by
  have h : (Finset.univ.filter (fun t : Fin k => j t = c)).Nonempty := ⟨t, by simp [ht]⟩
  have hm := firstIdx_mem h
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hm
  exact hm

/-- No earlier time visits the column `c`. -/
private theorem colTime_le {k d : ℕ} [NeZero k] {j : Fin k → Fin d} {c : Fin d} {t : Fin k}
    (ht : j t = c) : colTime j c ≤ t :=
  firstIdx_le (by simp [ht])

/-- The predecessor of an index, `0` at `0`. -/
private def prevIdx {k : ℕ} (t : Fin k) : Fin k :=
  ⟨t.1 - 1, Nat.lt_of_le_of_lt (Nat.sub_le _ _) t.2⟩

/-- One cyclic step after the predecessor of a nonzero index is the index itself. -/
private theorem cycSucc_prevIdx {k : ℕ} {t : Fin k} (h : 1 ≤ t.1) : cycSucc (prevIdx t) = t := by
  have hlt := t.isLt
  have h1 : t.1 - 1 + 1 = t.1 := by omega
  ext
  simp only [cycSucc, prevIdx, h1]
  exact Nat.mod_eq_of_lt hlt

/-- The edge of the first visit to a row is an edge of the walk. -/
private theorem rowEdge_mem {k n d : ℕ} [NeZero k] {i : Fin k → Fin n} {j : Fin k → Fin d}
    {r : Fin n} {t : Fin k} (ht : i t = r) :
    (r, j (rowTime i r)) ∈ walkEdges i j := by
  rw [mem_walkEdges]
  have hmem : rowTime i r ∈
      Finset.univ.filter (fun s : Fin k => i s = r ∧ j s = j (rowTime i r)) := by
    simp [rowTime_spec ht]
  have hpos : 0 < (Finset.univ.filter
      (fun s : Fin k => i s = r ∧ j s = j (rowTime i r))).card :=
    Finset.card_pos.mpr ⟨_, hmem⟩
  have hw : walkMult i j (r, j (rowTime i r))
      = (Finset.univ.filter (fun s : Fin k => i s = r ∧ j s = j (rowTime i r))).card
        + (Finset.univ.filter (fun s : Fin k =>
            i s = r ∧ j (cycSucc s) = j (rowTime i r))).card := rfl
  omega

/-- The first visit to a column other than the start happens after time `0`. -/
private theorem one_le_colTime {k d : ℕ} [NeZero k] {j : Fin k → Fin d} {c : Fin d} {t : Fin k}
    (ht : j t = c) (hne : c ≠ j 0) : 1 ≤ (colTime j c).1 := by
  by_contra hcon
  have hz : colTime j c = 0 := by
    have : (colTime j c).1 = (0 : Fin k).1 := by
      simp only [Fin.val_zero]
      omega
    exact Fin.val_injective this
  exact hne (by rw [← colTime_spec ht, hz])

/-- The edge just before the first visit to a column is an edge of the walk. -/
private theorem colEdge_mem {k n d : ℕ} [NeZero k] {i : Fin k → Fin n} {j : Fin k → Fin d}
    {c : Fin d} {t : Fin k} (ht : j t = c) (hne : c ≠ j 0) :
    (i (prevIdx (colTime j c)), c) ∈ walkEdges i j := by
  rw [mem_walkEdges]
  have hstep : j (cycSucc (prevIdx (colTime j c))) = c := by
    rw [cycSucc_prevIdx (one_le_colTime ht hne)]
    exact colTime_spec ht
  have hmem : prevIdx (colTime j c) ∈ Finset.univ.filter
      (fun s : Fin k => i s = i (prevIdx (colTime j c)) ∧ j (cycSucc s) = c) := by
    simp [hstep]
  have hpos : 0 < (Finset.univ.filter
      (fun s : Fin k => i s = i (prevIdx (colTime j c)) ∧ j (cycSucc s) = c)).card :=
    Finset.card_pos.mpr ⟨_, hmem⟩
  have hw : walkMult i j (i (prevIdx (colTime j c)), c)
      = (Finset.univ.filter (fun s : Fin k =>
            i s = i (prevIdx (colTime j c)) ∧ j s = c)).card
        + (Finset.univ.filter (fun s : Fin k =>
            i s = i (prevIdx (colTime j c)) ∧ j (cycSucc s) = c)).card := rfl
  omega

/-- X3. A walk has at most one more vertex than it has edges. Every row goes to the edge of
its first visit, and every column except the start `j 0` to the edge just before its first
visit; the two families are disjoint and injective, so they hold `v - 1` distinct edges. -/
private theorem walkVerts_le_edges_succ {k n d : ℕ} [NeZero k] (i : Fin k → Fin n)
    (j : Fin k → Fin d) : walkVerts i j ≤ (walkEdges i j).card + 1 := by
  set R : Finset (Fin n) := Finset.univ.image i with hRdef
  set C : Finset (Fin d) := Finset.univ.image j with hCdef
  set fr : Fin n → Fin n × Fin d := fun r => (r, j (rowTime i r)) with hfr
  set fc : Fin d → Fin n × Fin d := fun c => (i (prevIdx (colTime j c)), c) with hfc
  set E1 : Finset (Fin n × Fin d) := R.image fr with hE1
  set E2 : Finset (Fin n × Fin d) := (C.erase (j 0)).image fc with hE2
  have hj0 : j 0 ∈ C := Finset.mem_image_of_mem j (Finset.mem_univ _)
  have hCpos : 1 ≤ C.card := Finset.card_pos.mpr ⟨j 0, hj0⟩
  have hmemR : ∀ r ∈ R, ∃ t : Fin k, i t = r := by
    intro r hr
    rw [hRdef, Finset.mem_image] at hr
    obtain ⟨t, _, ht⟩ := hr
    exact ⟨t, ht⟩
  have hmemC : ∀ c ∈ C, ∃ t : Fin k, j t = c := by
    intro c hc
    rw [hCdef, Finset.mem_image] at hc
    obtain ⟨t, _, ht⟩ := hc
    exact ⟨t, ht⟩
  have hinjr : Function.Injective fr := by
    intro a b hab
    have := congrArg Prod.fst hab
    simpa [hfr] using this
  have hinjc : Function.Injective fc := by
    intro a b hab
    have := congrArg Prod.snd hab
    simpa [hfc] using this
  have hE1card : E1.card = R.card := by
    rw [hE1, Finset.card_image_of_injective _ hinjr]
  have hE2card : E2.card = C.card - 1 := by
    rw [hE2, Finset.card_image_of_injective _ hinjc, Finset.card_erase_of_mem hj0]
  have hsub : E1 ∪ E2 ⊆ walkEdges i j := by
    intro e he
    rcases Finset.mem_union.mp he with h | h
    · rw [hE1, Finset.mem_image] at h
      obtain ⟨r, hr, hre⟩ := h
      obtain ⟨t, ht⟩ := hmemR r hr
      rw [← hre, hfr]
      exact rowEdge_mem ht
    · rw [hE2, Finset.mem_image] at h
      obtain ⟨c, hc, hce⟩ := h
      have hne : c ≠ j 0 := (Finset.mem_erase.mp hc).1
      obtain ⟨t, ht⟩ := hmemC c (Finset.mem_erase.mp hc).2
      rw [← hce, hfc]
      exact colEdge_mem ht hne
  have hdisj : Disjoint E1 E2 := by
    refine Finset.disjoint_left.mpr fun e he1 he2 => ?_
    rw [hE1, Finset.mem_image] at he1
    rw [hE2, Finset.mem_image] at he2
    obtain ⟨r, hr, hre⟩ := he1
    obtain ⟨c, hc, hce⟩ := he2
    have hne : c ≠ j 0 := (Finset.mem_erase.mp hc).1
    obtain ⟨tr0, htr0⟩ := hmemR r hr
    obtain ⟨tc0, htc0⟩ := hmemC c (Finset.mem_erase.mp hc).2
    have heq : fr r = fc c := by rw [hre, hce]
    have h1 : r = i (prevIdx (colTime j c)) := by
      have := congrArg Prod.fst heq
      simpa [hfr, hfc] using this
    have h2 : j (rowTime i r) = c := by
      have := congrArg Prod.snd heq
      simpa [hfr, hfc] using this
    have hle1 : rowTime i r ≤ prevIdx (colTime j c) := rowTime_le h1.symm
    have hle2 : colTime j c ≤ rowTime i r := colTime_le h2
    have hpos := one_le_colTime htc0 hne
    have hv1 := (Fin.le_def).mp hle1
    have hv2 := (Fin.le_def).mp hle2
    have hprev : (prevIdx (colTime j c)).1 = (colTime j c).1 - 1 := rfl
    omega
  have hcard : R.card + (C.card - 1) ≤ (walkEdges i j).card := by
    calc R.card + (C.card - 1) = E1.card + E2.card := by rw [hE1card, hE2card]
      _ = (E1 ∪ E2).card := (Finset.card_union_of_disjoint hdisj).symm
      _ ≤ (walkEdges i j).card := Finset.card_le_card hsub
  have hv : walkVerts i j = R.card + C.card := rfl
  omega

/-- X5. An excess walk has at most `k` vertices: one entry of multiplicity at least 3 forces
the edge count below `k`, and the vertex count is at most the edge count plus one. -/
private theorem walkVerts_le_of_excess {k n d : ℕ} [NeZero k] {i : Fin k → Fin n}
    {j : Fin k → Fin d} (h : IsExcess i j) : walkVerts i j ≤ k := by
  obtain ⟨h1, e₀, he₀⟩ := h
  have h2 := walkVerts_le_edges_succ i j
  have h3 := edges_card_lt h1 he₀
  omega

/-- A walk with no entry of multiplicity 1 and `k + 1` vertices is a tree walk: `k + 1 ≤ |E| + 1`
and `|E| ≤ k` pin the edge count to `k`, and then no entry can carry multiplicity 3. -/
private theorem isPaired_of_verts {k n d : ℕ} [NeZero k] {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) (hv : walkVerts i j = k + 1) : IsPaired i j := by
  intro e
  by_cases he : walkMult i j e = 0
  · exact Or.inl he
  · refine Or.inr ?_
    by_contra hne
    have h3 : 3 ≤ walkMult i j e := by have := h e; omega
    have h4 := edges_card_lt h h3
    have h5 := walkVerts_le_edges_succ i j
    omega

/-! ### The weight of one class-`B` walk -/

/-- The weight of a walk with no entry of multiplicity 1 is at most `K^(2(k - |E|))`. Off the
edge set the factor is `∫ |x|^0 ∂ρ = 1`; on it `∫ |x|^m ∂ρ ≤ K^(m-2) ∫ x^2 ∂ρ ≤ K^(m-2)`, and
the exponents sum to `2k - 2|E|`. -/
private theorem weight_le {ρ : Measure ℝ} {K : ℝ} (hK : 1 ≤ K) (hρ : TruncNoiseLaw ρ K)
    {n d k : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d} (h : NoSingle i j) :
    (∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ)
      ≤ K ^ (2 * (k - (walkEdges i j).card)) := by
  have : IsProbabilityMeasure ρ := hρ.prob
  have hK0 : (0 : ℝ) ≤ K := by linarith
  have hstep : (∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ)
      ≤ ∏ e : Fin n × Fin d,
          (if walkMult i j e ≠ 0 then K ^ (walkMult i j e - 2) else 1) := by
    refine Finset.prod_le_prod (fun e _ => integral_nonneg fun x => by positivity) ?_
    intro e _
    by_cases he : walkMult i j e = 0
    · simp [he]
    · rw [if_pos he]
      have h2 : 2 ≤ walkMult i j e := two_le_walkMult h (mem_walkEdges.mpr he)
      calc ∫ x, |x| ^ (walkMult i j e) ∂ρ ≤ K ^ (walkMult i j e - 2) * ∫ x, x ^ 2 ∂ρ :=
            integral_abs_pow_le hK0 hρ h2
        _ ≤ K ^ (walkMult i j e - 2) * 1 :=
            mul_le_mul_of_nonneg_left hρ.var_le (pow_nonneg hK0 _)
        _ = K ^ (walkMult i j e - 2) := mul_one _
  have hfilt : (∏ e : Fin n × Fin d,
        (if walkMult i j e ≠ 0 then K ^ (walkMult i j e - 2) else 1))
      = ∏ e ∈ walkEdges i j, K ^ (walkMult i j e - 2) := (Finset.prod_filter _ _).symm
  have hexp : ∑ e ∈ walkEdges i j, (walkMult i j e - 2) = 2 * (k - (walkEdges i j).card) := by
    have hid : ∑ e ∈ walkEdges i j, walkMult i j e
        = ∑ e ∈ walkEdges i j, ((walkMult i j e - 2) + 2) :=
      Finset.sum_congr rfl fun e he => by have := two_le_walkMult h he; omega
    rw [Finset.sum_add_distrib, Finset.sum_const, smul_eq_mul, sum_walkMult_edges i j] at hid
    have hEk : (walkEdges i j).card ≤ k := edges_card_le h
    omega
  calc (∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ)
      ≤ ∏ e : Fin n × Fin d, (if walkMult i j e ≠ 0 then K ^ (walkMult i j e - 2) else 1) :=
        hstep
    _ = ∏ e ∈ walkEdges i j, K ^ (walkMult i j e - 2) := hfilt
    _ = K ^ (∑ e ∈ walkEdges i j, (walkMult i j e - 2)) := Finset.prod_pow_eq_pow_sum _ _ _
    _ = K ^ (2 * (k - (walkEdges i j).card)) := by rw [hexp]

/-- The vertex deficiency dominates the multiplicity excess, so the weight of a class-`B` walk
is at most `K^(2(k + 1 - v))`. This is the one step that needs `1 ≤ K`. -/
private theorem weight_le_deficiency {ρ : Measure ℝ} {K : ℝ} (hK : 1 ≤ K)
    (hρ : TruncNoiseLaw ρ K) {n d k : ℕ} [NeZero k] {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) :
    (∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ)
      ≤ K ^ (2 * (k + 1 - walkVerts i j)) := by
  refine (weight_le hK hρ h).trans (pow_le_pow_right₀ hK ?_)
  have h1 := walkVerts_le_edges_succ i j
  omega

open scoped Classical in
/-- X6. Grade the class-`B` walks by the vertex deficiency `D = k + 1 - v`, which lies in
`[1, k + 1]` by X5. Each such walk has weight at most `K^(2D)`, and the walks of a given `D`
sit inside the class `clsCard n d k (k + 1 - D)`. -/
private theorem excessSum_le_graded {ρ : Measure ℝ} {K : ℝ} (hK : 1 ≤ K)
    (hρ : TruncNoiseLaw ρ K) {n d k : ℕ} [NeZero k] :
    excessSum ρ n d k
      ≤ ∑ D ∈ Finset.Icc 1 (k + 1), K ^ (2 * D) * (clsCard n d k (k + 1 - D) : ℝ) := by
  have hK0 : (0 : ℝ) ≤ K := by linarith
  have hmaps : ∀ p ∈ Finset.univ.filter
      (fun p : (Fin k → Fin n) × (Fin k → Fin d) => IsExcess p.1 p.2),
      (k + 1 - walkVerts p.1 p.2) ∈ Finset.Icc 1 (k + 1) := by
    intro p hp
    have hEx : IsExcess p.1 p.2 := (Finset.mem_filter.mp hp).2
    have hv := walkVerts_le_of_excess hEx
    simp only [Finset.mem_Icc]
    omega
  have hprod : excessSum ρ n d k
      = ∑ p : (Fin k → Fin n) × (Fin k → Fin d),
          (if IsExcess p.1 p.2 then
            ∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult p.1 p.2 e) ∂ρ else 0) := by
    rw [excessSum, Fintype.sum_prod_type]
  calc excessSum ρ n d k
      = ∑ p : (Fin k → Fin n) × (Fin k → Fin d),
          (if IsExcess p.1 p.2 then
            ∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult p.1 p.2 e) ∂ρ else 0) := hprod
    _ ≤ ∑ p : (Fin k → Fin n) × (Fin k → Fin d),
          (if IsExcess p.1 p.2 then K ^ (2 * (k + 1 - walkVerts p.1 p.2)) else 0) := by
        refine Finset.sum_le_sum fun p _ => ?_
        by_cases hp : IsExcess p.1 p.2
        · rw [if_pos hp, if_pos hp]
          exact weight_le_deficiency hK hρ hp.1
        · rw [if_neg hp, if_neg hp]
    _ = ∑ p ∈ Finset.univ.filter
          (fun p : (Fin k → Fin n) × (Fin k → Fin d) => IsExcess p.1 p.2),
          K ^ (2 * (k + 1 - walkVerts p.1 p.2)) := (Finset.sum_filter _ _).symm
    _ = ∑ D ∈ Finset.Icc 1 (k + 1), ∑ p ∈ (Finset.univ.filter
          (fun p : (Fin k → Fin n) × (Fin k → Fin d) => IsExcess p.1 p.2)).filter
          (fun p => k + 1 - walkVerts p.1 p.2 = D),
          K ^ (2 * (k + 1 - walkVerts p.1 p.2)) :=
        (Finset.sum_fiberwise_of_maps_to hmaps _).symm
    _ ≤ ∑ D ∈ Finset.Icc 1 (k + 1), K ^ (2 * D) * (clsCard n d k (k + 1 - D) : ℝ) := by
        refine Finset.sum_le_sum fun D _ => ?_
        have hconst : ∀ p ∈ (Finset.univ.filter
            (fun p : (Fin k → Fin n) × (Fin k → Fin d) => IsExcess p.1 p.2)).filter
            (fun p => k + 1 - walkVerts p.1 p.2 = D),
            K ^ (2 * (k + 1 - walkVerts p.1 p.2)) = K ^ (2 * D) := by
          intro p hp
          rw [(Finset.mem_filter.mp hp).2]
        have hcls : clsCard n d k (k + 1 - D)
            = (Finset.univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
                NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = k + 1 - D)).card := rfl
        have hcard : ((Finset.univ.filter
              (fun p : (Fin k → Fin n) × (Fin k → Fin d) => IsExcess p.1 p.2)).filter
              (fun p => k + 1 - walkVerts p.1 p.2 = D)).card
            ≤ clsCard n d k (k + 1 - D) := by
          rw [hcls]
          refine Finset.card_le_card fun p hp => ?_
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hp ⊢
          obtain ⟨hEx, hDp⟩ := hp
          have hv := walkVerts_le_of_excess hEx
          exact ⟨hEx.1, by omega⟩
        rw [Finset.sum_congr rfl hconst, Finset.sum_const, nsmul_eq_mul, mul_comm]
        refine mul_le_mul_of_nonneg_left ?_ (pow_nonneg hK0 _)
        exact_mod_cast Nat.cast_le.mpr hcard

/-! ### The tree class against the Gaussian trace -/

/-- Every moment of the standard Gaussian exists, so unit E2 applies to it. A copy of the
private lemma of `Compare.lean`, which this file cannot reach. -/
private theorem integrable_pow_gauss (p : ℕ) :
    Integrable (fun x : ℝ => x ^ p) (gaussianReal 0 1) := by
  have hm : MemLp (id : ℝ → ℝ) ((p : ℕ) : ENNReal) (gaussianReal 0 1) :=
    memLp_id_gaussianReal' _ (by simp)
  have h := hm.integrable_norm_pow'
  rw [← integrable_norm_iff (by fun_prop)]
  simpa [norm_pow] using h

/-- Every moment of the standard Gaussian is nonnegative: an even power is nonnegative
pointwise, an odd moment is `0` by symmetry. A copy of the private lemma of `Compare.lean`. -/
private theorem integral_pow_gauss_nonneg (m : ℕ) :
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

open scoped Classical in
/-- X7. The tree class is at most the Gaussian trace moment. Unit E2 at the Gaussian law writes
the trace moment as a sum of products of Gaussian moments; every summand is nonnegative, and a
walk with no entry of multiplicity 1 and `k + 1` vertices contributes exactly `1`. -/
private theorem clsCard_le_gaussian {n d k : ℕ} [NeZero k] (hk : 1 ≤ k) :
    (clsCard n d k (k + 1) : ℝ)
      ≤ ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix n d) := by
  rw [gaussianMatrix_eq_noiseMatrix,
    integral_trace_gram_pow (ρ := gaussianReal 0 1) (fun p _ => integrable_pow_gauss p) hk]
  have hsum : ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
        (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult i j e) ∂(gaussianReal 0 1))
      = ∑ p : (Fin k → Fin n) × (Fin k → Fin d),
        ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult p.1 p.2 e) ∂(gaussianReal 0 1) := by
    rw [Fintype.sum_prod_type]
  rw [hsum]
  have hcls : clsCard n d k (k + 1)
      = (Finset.univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
          NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = k + 1)).card := rfl
  have hone : ∀ p ∈ Finset.univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
        NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = k + 1),
      (∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult p.1 p.2 e) ∂(gaussianReal 0 1)) = 1 := by
    intro p hp
    have hp' := (Finset.mem_filter.mp hp).2
    have hpair : IsPaired p.1 p.2 := isPaired_of_verts hp'.1 hp'.2
    refine Finset.prod_eq_one fun e _ => ?_
    rcases hpair e with hm | hm
    · rw [hm]; simp
    · rw [hm]; exact noiseLaw_gaussian.var
  calc (clsCard n d k (k + 1) : ℝ)
      = ∑ _p ∈ Finset.univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
          NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = k + 1), (1 : ℝ) := by
        rw [hcls]; simp
    _ = ∑ p ∈ Finset.univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
          NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = k + 1),
          ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult p.1 p.2 e) ∂(gaussianReal 0 1) :=
        Finset.sum_congr rfl fun p hp => (hone p hp).symm
    _ ≤ ∑ p : (Fin k → Fin n) × (Fin k → Fin d),
          ∏ e : Fin n × Fin d, ∫ x, x ^ (walkMult p.1 p.2 e) ∂(gaussianReal 0 1) :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) fun p _ _ =>
          Finset.prod_nonneg fun e _ => integral_pow_gauss_nonneg _

end ExcessAux

/-- X. The weighted sum over `B` is at most `2 f` times the Gaussian trace, with
`f = fkFactor K n d k`. The bound is RELATIVE to the tree-walk class: the class with `k+1`
distinct vertices sits inside `A`, every Gaussian weight there is 1, so its size is at most
`∫ trace ((Y_GᵀY_G)^k)`, which unit G bounds sharply. An absolute count of the shape
`2^k k^(6j) n^s d^t` loses `4^k` against `(1+√c)^(2k)` and does not close the assembly. -/
theorem excessSum_le {ρ : Measure ℝ} {K : ℝ} (hK : 1 ≤ K) (hρ : TruncNoiseLaw ρ K)
    {n d k : ℕ} (hk : 1 ≤ k) (hn : 0 < n) (hd : 0 < d)
    (hf : fkFactor K n d k ≤ 1 / 2) :
    excessSum ρ n d k
      ≤ 2 * fkFactor K n d k * ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix n d) := by
  have : NeZero k := ⟨by omega⟩
  have hK0 : (0 : ℝ) < K := by linarith
  have hMpos : (0 : ℝ) < min (n : ℝ) (d : ℝ) :=
    lt_min (by exact_mod_cast hn) (by exact_mod_cast hd)
  have hfnn : 0 ≤ fkFactor K n d k := by
    rw [fkFactor]
    exact div_nonneg (by positivity) hMpos.le
  -- `hf` with `1 ≤ K` and `1 ≤ k` forces `min (n, d) ≥ 2 (2k)^12 ≥ 2k`.
  have hk1 : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
  have hx6 : (2 * (k : ℝ)) ^ 1 ≤ (2 * (k : ℝ)) ^ 12 :=
    pow_le_pow_right₀ (by linarith) (by norm_num)
  rw [pow_one] at hx6
  have hK2 : (1 : ℝ) ≤ K ^ 2 := one_le_pow₀ hK
  have hnd : 2 * k ≤ min n d := by
    have hf' : K ^ 2 * (2 * (k : ℝ)) ^ 12 / min (n : ℝ) (d : ℝ) ≤ 1 / 2 := hf
    rw [div_le_iff₀ hMpos] at hf'
    have hpow6 : (0 : ℝ) < (2 * (k : ℝ)) ^ 12 := by positivity
    have h1 : (2 * (k : ℝ)) ^ 12 ≤ K ^ 2 * (2 * (k : ℝ)) ^ 12 := by nlinarith
    have h3 : ((2 * k : ℕ) : ℝ) ≤ ((min n d : ℕ) : ℝ) := by push_cast; linarith
    exact_mod_cast h3
  have hCnn : (0 : ℝ) ≤ (clsCard n d k (k + 1) : ℝ) := Nat.cast_nonneg _
  have hterm : ∀ D ∈ Finset.Icc 1 (k + 1),
      K ^ (2 * D) * (clsCard n d k (k + 1 - D) : ℝ)
        ≤ fkFactor K n d k ^ D * (clsCard n d k (k + 1) : ℝ) := by
    intro D hD
    rw [Finset.mem_Icc] at hD
    have hcount := clsCard_mul_le (n := n) (d := d) (k := k) (v := k + 1 - D) hk
      (by omega) hnd
    have hexp : k + 1 - (k + 1 - D) = D := by omega
    rw [hexp] at hcount
    have hcountR : (clsCard n d k (k + 1 - D) : ℝ) * min (n : ℝ) (d : ℝ) ^ D
        ≤ ((2 * (k : ℝ)) ^ 12) ^ D * (clsCard n d k (k + 1) : ℝ) := by
      have hc := (Nat.cast_le (α := ℝ)).mpr hcount
      push_cast at hc
      exact hc
    have hMD : (0 : ℝ) < min (n : ℝ) (d : ℝ) ^ D := pow_pos hMpos D
    have hfD : fkFactor K n d k ^ D
        = K ^ (2 * D) * ((2 * (k : ℝ)) ^ 12) ^ D / min (n : ℝ) (d : ℝ) ^ D := by
      rw [fkFactor, div_pow, mul_pow, pow_mul]
    rw [hfD, div_mul_eq_mul_div, le_div_iff₀ hMD]
    calc K ^ (2 * D) * (clsCard n d k (k + 1 - D) : ℝ) * min (n : ℝ) (d : ℝ) ^ D
        = K ^ (2 * D) * ((clsCard n d k (k + 1 - D) : ℝ) * min (n : ℝ) (d : ℝ) ^ D) := by ring
      _ ≤ K ^ (2 * D) * (((2 * (k : ℝ)) ^ 12) ^ D * (clsCard n d k (k + 1) : ℝ)) :=
          mul_le_mul_of_nonneg_left hcountR (by positivity)
      _ = K ^ (2 * D) * ((2 * (k : ℝ)) ^ 12) ^ D * (clsCard n d k (k + 1) : ℝ) := by ring
  calc excessSum ρ n d k
      ≤ ∑ D ∈ Finset.Icc 1 (k + 1), K ^ (2 * D) * (clsCard n d k (k + 1 - D) : ℝ) :=
        ExcessAux.excessSum_le_graded hK hρ
    _ ≤ ∑ D ∈ Finset.Icc 1 (k + 1), fkFactor K n d k ^ D * (clsCard n d k (k + 1) : ℝ) :=
        Finset.sum_le_sum hterm
    _ = (∑ D ∈ Finset.Icc 1 (k + 1), fkFactor K n d k ^ D) * (clsCard n d k (k + 1) : ℝ) := by
        rw [Finset.sum_mul]
    _ ≤ 2 * fkFactor K n d k * (clsCard n d k (k + 1) : ℝ) :=
        mul_le_mul_of_nonneg_right (ExcessAux.geom_tail hfnn hf _) hCnn
    _ ≤ 2 * fkFactor K n d k * ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix n d) :=
        mul_le_mul_of_nonneg_left (ExcessAux.clsCard_le_gaussian hk) (by linarith)

end Edge
end StackedSVD
