/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Defs

/-! # Stage 3, unit X, part 1: the walk classes and the Furedi-Komlos count

The combinatorial input of the excess bound (`notes/stage3_edge.md`, route C;
`notes/STAGE3_CAMPAIGN.md`; the plan of unit X in the campaign note). A closed walk
`j 0 → i 0 → j 1 → i 1 → ... → j (k-1) → i (k-1) → j 0` on the complete bipartite graph
`Fin n × Fin d` is graded by its vertex count `walkVerts i j` (distinct rows plus distinct
columns). `clsCard n d k v` counts the walks with no entry of multiplicity 1 and exactly `v`
vertices. The walks with `k + 1` vertices are the tree walks (every entry has multiplicity
exactly 2, so they lie in class `A` of `Defs.lean`).

`clsCard_mul_le` (`CountBound.lean`) is the count: one unit of vertex deficiency
`D = k + 1 - v` costs at most `(2k)^12 / min(n, d)` relative to the tree class, in the
`D`-step form the assembly reads (Furedi and Komlos 1981, page 237; Anderson, Guionnet and
Zeitouni, Lemma 2.1.23). This file holds the definitions and the elementary walk facts the
count reads, and three modules on top of it prove the count. `Code.lean` encodes a walk by its
steps: each of the `2k` steps is innovative (a first visit; there are `v - 1` of these and
they build a spanning tree), forced (an old step along the one entry at the current vertex
with an odd number of traversals so far, the Eulerian parity rule) or bad, and the code
determines the walk. `Dyck.lean` counts the tree walks of the target cell from below, through
the contour walk of a Dyck word and the Narayana lower bound. `CountBound.lean` assembles the
two halves and sums the cells of one vertex count. A bad step is determined by its position
(at most `2k` values) and its endpoint (at most `k + 1` values), and the number of bad steps
is at most `4 D` (`card_codeP_le` of `Code.lean`). The code must also carry the
innovative-or-forced pattern: the tree and the bad steps alone do not determine the walk
(`x c1 x c2 x c1 w c1 x` and `x c1 w c1 x c1 x c2 x` share both, 2026-09-10). Caution: the
bracket rule (a forced step retraces the innovative step the bracket order matches) is wrong;
the `2L`-cycle traversed twice has `D = 1` and `k` bad steps under it. The row and column
counts of the tree skeleton equal those of the walk, so the count is cell by cell in `(s, t)`
and no `4^k` is lost against `(1 + √c)^(2k)`.

Numerical evidence (the plan of unit X, `verify2.py`, 2026-09-10): the statement has 0
violations for `k ≤ 7` at `(n, d)` from `(k + 1, k + 1)` to `(10^12, 10^9)`, and it still holds
with `(2k)^2` in place of `(2k)^12`. The hypotheses are needed: at `k = 0` the empty walk has
`v = 0 < k + 1` and no tree walk exists; at `n = d = 3`, `k = 6`, `clsCard 3 3 6 7 = 0` while
`clsCard 3 3 6 6 > 0`, so `2k ≤ min(n, d)` cannot be dropped. -/

open Finset

namespace StackedSVD
namespace Edge

/-- The edge set of a walk: the entries of multiplicity at least 1. -/
def walkEdges {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin n × Fin d) :=
  univ.filter (fun e => walkMult i j e ≠ 0)

/-- The vertex count of a walk: the distinct rows plus the distinct columns. -/
def walkVerts {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : ℕ :=
  (univ.image i).card + (univ.image j).card

/-- No entry of multiplicity 1 (the union of the classes `A` and `B` of `Defs.lean`). -/
def NoSingle {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Prop :=
  ∀ e : Fin n × Fin d, walkMult i j e ≠ 1

open scoped Classical in
/-- The number of walks of length `2k` on `Fin n × Fin d` with no entry of multiplicity 1 and
exactly `v` vertices. -/
noncomputable def clsCard (n d k v : ℕ) : ℕ :=
  (univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
    NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = v)).card

/-! ## The elementary walk facts -/

/-- The entry the walk traverses at the step `j t → i t` is an edge of the walk. -/
theorem mem_walkEdges_fst {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d} (t : Fin k) :
    (i t, j t) ∈ walkEdges i j := by
  classical
  rw [walkEdges, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, ?_⟩
  rw [walkMult]
  have h : (({s : Fin k | i s = (i t, j t).1 ∧ j s = (i t, j t).2} : Finset (Fin k))).Nonempty :=
    ⟨t, by simp⟩
  have := Finset.card_pos.mpr h
  omega

/-- The entry the walk traverses at the step `i t → j (cycSucc t)` is an edge of the walk. -/
theorem mem_walkEdges_snd {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d} (t : Fin k) :
    (i t, j (cycSucc t)) ∈ walkEdges i j := by
  classical
  rw [walkEdges, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, ?_⟩
  rw [walkMult]
  have h : (({s : Fin k | i s = (i t, j (cycSucc t)).1 ∧
      j (cycSucc s) = (i t, j (cycSucc t)).2} : Finset (Fin k))).Nonempty := ⟨t, by simp⟩
  have := Finset.card_pos.mpr h
  omega

/-- The value of the cyclic successor. -/
theorem cycSucc_val {k : ℕ} (x : Fin k) : (cycSucc x).1 = (x.1 + 1) % k := rfl

/-- A walk has no entry of multiplicity 1 as soon as every entry its first step family
traverses is also traversed by its second family, and conversely. -/
theorem noSingle_of_families {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h1 : ∀ x : Fin k, ∃ y : Fin k, i y = i x ∧ j (cycSucc y) = j x)
    (h2 : ∀ x : Fin k, ∃ y : Fin k, i y = i x ∧ j y = j (cycSucc x)) : NoSingle i j := by
  classical
  intro e he
  rw [walkMult] at he
  by_cases hAne : (({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k))).Nonempty
  · obtain ⟨x, hxA⟩ := hAne
    have hx := hxA
    rw [Finset.mem_filter] at hx
    obtain ⟨y, hy1, hy2⟩ := h1 x
    have hyB : y ∈ ({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k)) := by
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ _, by rw [hy1]; exact hx.2.1, by rw [hy2]; exact hx.2.2⟩
    have hcA := Finset.card_pos.mpr ⟨x, hxA⟩
    have hcB := Finset.card_pos.mpr ⟨y, hyB⟩
    omega
  · rw [Finset.not_nonempty_iff_eq_empty] at hAne
    rw [hAne, Finset.card_empty] at he
    have hBne : (({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k))).Nonempty :=
      Finset.card_pos.mp (by omega)
    obtain ⟨x, hxB⟩ := hBne
    have hx := hxB
    rw [Finset.mem_filter] at hx
    obtain ⟨y, hy1, hy2⟩ := h2 x
    have hyA : y ∈ ({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k)) := by
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ _, by rw [hy1]; exact hx.2.1, by rw [hy2]; exact hx.2.2⟩
    rw [hAne] at hyA
    exact absurd hyA (Finset.notMem_empty y)
/-- The fibers of a map out of `Fin k` have `k` elements in total. -/
private theorem sum_card_fiber {k n d : ℕ} (f : Fin k → Fin n × Fin d) :
    ∑ e : Fin n × Fin d, (univ.filter fun t => f t = e).card = k := by
  classical
  exact (Finset.card_eq_sum_card_fiberwise (fun x _ => mem_univ (f x))).symm.trans
    (Finset.card_fin k)

/-- A fiber written with a conjunction is a fiber of the pair map. -/
private theorem filter_pair_eq {k n d : ℕ} (a : Fin k → Fin n) (b : Fin k → Fin d)
    (e : Fin n × Fin d) :
    ({t : Fin k | a t = e.1 ∧ b t = e.2} : Finset (Fin k))
      = univ.filter fun t => (a t, b t) = e := by
  apply Finset.filter_congr
  intro t _
  simp [Prod.ext_iff]

/-- X2: the total multiplicity of a walk is its number of steps, `2k`. -/
theorem sum_walkMult {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    ∑ e : Fin n × Fin d, walkMult i j e = 2 * k := by
  classical
  simp only [walkMult, filter_pair_eq, Finset.sum_add_distrib, sum_card_fiber]
  ring

/-- The multiplicity of an edge of the walk is at least 2 when no entry has multiplicity 1. -/
theorem two_le_walkMult {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) {e : Fin n × Fin d} (he : e ∈ walkEdges i j) : 2 ≤ walkMult i j e := by
  classical
  rw [walkEdges, Finset.mem_filter] at he
  have := h e
  omega

/-- X4: a walk with no entry of multiplicity 1 has at most `k` edges. -/
theorem card_walkEdges_le {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) : (walkEdges i j).card ≤ k := by
  classical
  have hlow : (walkEdges i j).card * 2 ≤ ∑ e ∈ walkEdges i j, walkMult i j e := by
    calc (walkEdges i j).card * 2 = ∑ _e ∈ walkEdges i j, 2 := by
          rw [Finset.sum_const, smul_eq_mul]
      _ ≤ _ := Finset.sum_le_sum fun e he => two_le_walkMult h he
  have hsub : ∑ e ∈ walkEdges i j, walkMult i j e ≤ 2 * k := by
    rw [← sum_walkMult i j]
    exact Finset.sum_le_sum_of_subset (Finset.subset_univ _)
  omega

/-! ## The first-visit injection -/

/-- The least element of a `Finset (Fin k)`, with a default value for the empty set. -/
private noncomputable def firstVisit {k : ℕ} (s : Finset (Fin k)) (t0 : Fin k) : Fin k :=
  if h : s.Nonempty then s.min' h else t0

/-- The first visit lies in the set when the set is nonempty. -/
private theorem firstVisit_mem {k : ℕ} {s : Finset (Fin k)} {t0 : Fin k} (h : s.Nonempty) :
    firstVisit s t0 ∈ s := by
  rw [firstVisit, dif_pos h]; exact Finset.min'_mem _ _

/-- The first visit is the least element of the set. -/
private theorem firstVisit_le {k : ℕ} {s : Finset (Fin k)} {t0 t : Fin k} (ht : t ∈ s) :
    firstVisit s t0 ≤ t := by
  rw [firstVisit, dif_pos ⟨t, ht⟩]; exact Finset.min'_le _ _ ht

/-- The predecessor on `Fin k`, with `0` fixed. -/
private def predFin {k : ℕ} (t : Fin k) : Fin k :=
  ⟨t.1 - 1, Nat.lt_of_le_of_lt (Nat.sub_le _ _) t.2⟩

/-- X3: a walk has at most one more vertex than it has edges. Every vertex except the column
`j 0` where the walk starts is charged to the edge of its first visit, and the charge is
injective: a row and a column can only collide against the minimality of both first visits. -/
theorem walkVerts_le_edges_succ {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    walkVerts i j ≤ (walkEdges i j).card + 1 := by
  classical
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp [walkVerts]
  · set t0 : Fin k := ⟨0, hk⟩ with ht0
    obtain ⟨tR, htRmem, htRle⟩ : ∃ f : Fin n → Fin k,
        (∀ r ∈ univ.image i, i (f r) = r) ∧ (∀ r t, i t = r → f r ≤ t) := by
      refine ⟨fun r => firstVisit (univ.filter fun t => i t = r) t0, ?_, ?_⟩
      · intro r hr
        obtain ⟨t, -, ht⟩ := Finset.mem_image.mp hr
        have := firstVisit_mem (t0 := t0) (s := univ.filter fun t => i t = r) ⟨t, by simp [ht]⟩
        simpa using this
      · intro r t ht
        exact firstVisit_le (by simp [ht])
    obtain ⟨tC, htCmem, htCle⟩ : ∃ f : Fin d → Fin k,
        (∀ c ∈ univ.image j, j (f c) = c) ∧ (∀ c t, j t = c → f c ≤ t) := by
      refine ⟨fun c => firstVisit (univ.filter fun t => j t = c) t0, ?_, ?_⟩
      · intro c hc
        obtain ⟨t, -, ht⟩ := Finset.mem_image.mp hc
        have := firstVisit_mem (t0 := t0) (s := univ.filter fun t => j t = c) ⟨t, by simp [ht]⟩
        simpa using this
      · intro c t ht
        exact firstVisit_le (by simp [ht])
    set Φ : Fin n ⊕ Fin d → Fin n × Fin d :=
      Sum.elim (fun r => (r, j (tR r))) (fun c => (i (predFin (tC c)), c)) with hΦ
    set S : Finset (Fin n ⊕ Fin d) :=
      (univ.image i).disjSum ((univ.image j).erase (j t0)) with hS
    have hposC : ∀ c ∈ (univ.image j).erase (j t0), 0 < (tC c).1 := by
      intro c hc
      have hcne : c ≠ j t0 := (Finset.mem_erase.mp hc).1
      have hjc : j (tC c) = c := htCmem c (Finset.mem_erase.mp hc).2
      rcases Nat.eq_zero_or_pos (tC c).1 with h0 | h0
      · exact absurd (by rw [← hjc]; congr 1; exact Fin.ext (by simp [ht0, h0])) hcne
      · exact h0
    have hstep : ∀ c ∈ (univ.image j).erase (j t0), cycSucc (predFin (tC c)) = tC c := by
      intro c hc
      have hpos := hposC c hc
      apply Fin.ext
      change ((predFin (tC c)).1 + 1) % k = (tC c).1
      have hval : (predFin (tC c)).1 + 1 = (tC c).1 := by simp only [predFin]; omega
      rw [hval]
      exact Nat.mod_eq_of_lt (tC c).2
    have hmaps : ∀ x ∈ S, Φ x ∈ walkEdges i j := by
      rintro (r | c) hx
      · have hr : r ∈ univ.image i := by rw [hS, Finset.inl_mem_disjSum] at hx; exact hx
        have hir : i (tR r) = r := htRmem r hr
        simpa [hΦ, hir] using mem_walkEdges_fst (i := i) (j := j) (tR r)
      · have hc : c ∈ (univ.image j).erase (j t0) := by
          rw [hS, Finset.inr_mem_disjSum] at hx; exact hx
        have hjc : j (tC c) = c := htCmem c (Finset.mem_erase.mp hc).2
        have hmem := mem_walkEdges_snd (i := i) (j := j) (predFin (tC c))
        rwa [hstep c hc, hjc] at hmem
    have hinj : Set.InjOn Φ ↑S := by
      rintro (r | c) hx (r' | c') hy hEq
      · simp only [hΦ, Sum.elim_inl, Prod.mk.injEq] at hEq
        rw [hEq.1]
      · exfalso
        simp only [hΦ, Sum.elim_inl, Sum.elim_inr, Prod.mk.injEq] at hEq
        have hc' : c' ∈ (univ.image j).erase (j t0) := by
          rw [Finset.mem_coe, hS, Finset.inr_mem_disjSum] at hy; exact hy
        have h1 : tC c' ≤ tR r := htCle c' (tR r) hEq.2
        have h2 : tR r ≤ predFin (tC c') := htRle r (predFin (tC c')) hEq.1.symm
        have hpos := hposC c' hc'
        rw [Fin.le_def] at h1 h2
        simp only [predFin] at h2
        omega
      · exfalso
        simp only [hΦ, Sum.elim_inl, Sum.elim_inr, Prod.mk.injEq] at hEq
        have hc : c ∈ (univ.image j).erase (j t0) := by
          rw [Finset.mem_coe, hS, Finset.inr_mem_disjSum] at hx; exact hx
        have h1 : tC c ≤ tR r' := htCle c (tR r') hEq.2.symm
        have h2 : tR r' ≤ predFin (tC c) := htRle r' (predFin (tC c)) hEq.1
        have hpos := hposC c hc
        rw [Fin.le_def] at h1 h2
        simp only [predFin] at h2
        omega
      · simp only [hΦ, Sum.elim_inr, Prod.mk.injEq] at hEq
        rw [hEq.2]
    have hcard := Finset.card_le_card_of_injOn Φ hmaps hinj
    have hj0 : j t0 ∈ univ.image j := Finset.mem_image.mpr ⟨t0, Finset.mem_univ _, rfl⟩
    rw [hS, Finset.card_disjSum, Finset.card_erase_of_mem hj0] at hcard
    have hj0card : 1 ≤ (univ.image j).card := Finset.card_pos.mpr ⟨j t0, hj0⟩
    simp only [walkVerts]
    omega

/-- A walk of positive length visits at least one row and at least one column. -/
theorem two_le_walkVerts {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    2 ≤ walkVerts i j := by
  have h1 : 1 ≤ (univ.image i).card :=
    Finset.card_pos.mpr ⟨i ⟨0, hk⟩, Finset.mem_image.mpr ⟨⟨0, hk⟩, Finset.mem_univ _, rfl⟩⟩
  have h2 : 1 ≤ (univ.image j).card :=
    Finset.card_pos.mpr ⟨j ⟨0, hk⟩, Finset.mem_image.mpr ⟨⟨0, hk⟩, Finset.mem_univ _, rfl⟩⟩
  simp only [walkVerts]
  omega

/-- A walk with no entry of multiplicity 1 has at most `k + 1` vertices. -/
theorem walkVerts_le_of_noSingle {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) : walkVerts i j ≤ k + 1 := by
  have h1 := walkVerts_le_edges_succ i j
  have h2 := card_walkEdges_le h
  omega

/-- X5: a walk of class `B` has at most `k` vertices. The entry of multiplicity at least 3
costs one edge, and the vertex count is at most the edge count plus one. -/
theorem walkVerts_le_of_excess {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : IsExcess i j) : walkVerts i j ≤ k := by
  classical
  obtain ⟨hns, e0, he0⟩ := h
  have he0mem : e0 ∈ walkEdges i j := by
    rw [walkEdges, Finset.mem_filter]
    exact ⟨Finset.mem_univ _, by omega⟩
  have hsplit : walkMult i j e0 + ∑ e ∈ (walkEdges i j).erase e0, walkMult i j e
      = ∑ e ∈ walkEdges i j, walkMult i j e := Finset.add_sum_erase _ _ he0mem
  have hlow : ((walkEdges i j).erase e0).card * 2
      ≤ ∑ e ∈ (walkEdges i j).erase e0, walkMult i j e := by
    calc ((walkEdges i j).erase e0).card * 2 = ∑ _e ∈ (walkEdges i j).erase e0, 2 := by
          rw [Finset.sum_const, smul_eq_mul]
      _ ≤ _ := Finset.sum_le_sum fun e he =>
          two_le_walkMult hns (Finset.mem_of_mem_erase he)
  have hsub : ∑ e ∈ walkEdges i j, walkMult i j e ≤ 2 * k := by
    rw [← sum_walkMult i j]
    exact Finset.sum_le_sum_of_subset (Finset.subset_univ _)
  have hcard : ((walkEdges i j).erase e0).card = (walkEdges i j).card - 1 :=
    Finset.card_erase_of_mem he0mem
  have hpos : 1 ≤ (walkEdges i j).card := Finset.card_pos.mpr ⟨e0, he0mem⟩
  have hv := walkVerts_le_edges_succ i j
  omega

/-- A walk with no entry of multiplicity 1 and `k + 1` vertices is a tree walk: every entry
multiplicity is `0` or `2`, so the walk lies in class `A`. -/
theorem isPaired_of_walkVerts_eq {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    (h : NoSingle i j) (hv : walkVerts i j = k + 1) : IsPaired i j := by
  intro e
  rcases Nat.lt_or_ge (walkMult i j e) 3 with hlt | hge
  · have := h e; omega
  · exact absurd (walkVerts_le_of_excess ⟨h, e, hge⟩) (by omega)

/-! ## The Furedi-Komlos count -/

open scoped Classical in
/-- The walks of length `2k` on `Fin n × Fin d` with no entry of multiplicity 1 and exactly `v`
vertices. `clsCard n d k v` is its cardinality. -/
noncomputable def walkSet (n d k v : ℕ) : Finset ((Fin k → Fin n) × (Fin k → Fin d)) :=
  univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
    NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = v)

/-- `clsCard` counts `walkSet`. -/
theorem clsCard_eq_card_walkSet (n d k v : ℕ) : clsCard n d k v = (walkSet n d k v).card := rfl

/-- Membership in `walkSet`. -/
theorem mem_walkSet {n d k v : ℕ} {p : (Fin k → Fin n) × (Fin k → Fin d)} :
    p ∈ walkSet n d k v ↔ NoSingle p.1 p.2 ∧ walkVerts p.1 p.2 = v := by
  classical
  rw [walkSet, Finset.mem_filter]
  exact ⟨fun h => h.2, fun h => ⟨Finset.mem_univ _, h⟩⟩

open scoped Classical in
/-- The `(rows, columns)` cell: the walks with no entry of multiplicity 1, exactly `s` distinct
rows and exactly `t` distinct columns. The Furedi-Komlos count is cell by cell, since a bound
that forgets the cell loses a factor `4 ^ k` against `(1 + sqrt c) ^ (2 k)`. -/
noncomputable def cellSet (n d k s t : ℕ) : Finset ((Fin k → Fin n) × (Fin k → Fin d)) :=
  univ.filter (fun p : (Fin k → Fin n) × (Fin k → Fin d) =>
    NoSingle p.1 p.2 ∧ (univ.image p.1).card = s ∧ (univ.image p.2).card = t)

/-- The number of walks in one `(rows, columns)` cell. -/
noncomputable def cellCard (n d k s t : ℕ) : ℕ := (cellSet n d k s t).card

/-- Membership in `cellSet`. -/
theorem mem_cellSet {n d k s t : ℕ} {p : (Fin k → Fin n) × (Fin k → Fin d)} :
    p ∈ cellSet n d k s t ↔
      NoSingle p.1 p.2 ∧ (univ.image p.1).card = s ∧ (univ.image p.2).card = t := by
  classical
  rw [cellSet, Finset.mem_filter]
  exact ⟨fun h => h.2, fun h => ⟨Finset.mem_univ _, h⟩⟩

/-- A cell with no row, or with no column, is empty at a walk of positive length. -/
theorem cellCard_eq_zero_of_zero {n d k s t : ℕ} (hk : 1 ≤ k) (h : s = 0 ∨ t = 0) :
    cellCard n d k s t = 0 := by
  classical
  rw [cellCard, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
  rintro p hp
  rw [mem_cellSet] at hp
  have h1 : 1 ≤ (univ.image p.1).card :=
    Finset.card_pos.mpr ⟨p.1 ⟨0, hk⟩, Finset.mem_image.mpr ⟨⟨0, hk⟩, Finset.mem_univ _, rfl⟩⟩
  have h2 : 1 ≤ (univ.image p.2).card :=
    Finset.card_pos.mpr ⟨p.2 ⟨0, hk⟩, Finset.mem_image.mpr ⟨⟨0, hk⟩, Finset.mem_univ _, rfl⟩⟩
  rcases h with h | h <;> omega

/-- The image of a family that takes exactly the first `s` values has `s` elements. -/
private theorem card_image_val {k n s : ℕ} (f : Fin k → Fin n) (hlt : ∀ x, (f x).1 < s)
    (hsurj : ∀ u, u < s → ∃ x : Fin k, (f x).1 = u) : (univ.image f).card = s := by
  classical
  have h1 : ((univ.image f).image Fin.val).card = (univ.image f).card :=
    Finset.card_image_of_injective _ Fin.val_injective
  have h2 : (univ.image f).image Fin.val = Finset.range s := by
    ext u
    simp only [Finset.mem_image, Finset.mem_range, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨e, ⟨x, rfl⟩, rfl⟩
      exact hlt x
    · intro hu
      obtain ⟨x, hx⟩ := hsurj u hu
      exact ⟨f x, ⟨x, rfl⟩, hx⟩
  rw [← h1, h2, Finset.card_range]

/-- The row index of the caterpillar walk: `x + 1` while `x + 2 ≤ s`, then `0`. -/
private def catRow (s x : ℕ) : ℕ := if x + 2 ≤ s then x + 1 else 0

/-- The column index of the caterpillar walk: `0` while `x + 1 ≤ s`, then `x + 1 - s`. -/
private def catCol (s x : ℕ) : ℕ := if x + 1 ≤ s then 0 else x + 1 - s

/-- Every cell of the tree class is nonempty. The caterpillar walk hangs the rows
`r 1, ..., r (s-1)` on the column `c 0` and the columns `c 1, ..., c (t-1)` on the row `r 0`,
and walks it once: `c 0 → r 1 → c 0 → ... → r (s-1) → c 0 → r 0 → c 1 → r 0 → ... →
c (t-1) → r 0 → c 0`. It has `s` rows, `t` columns and no entry of multiplicity 1. -/
theorem one_le_cellCard {n d k s t : ℕ} (hs : 1 ≤ s) (ht : 1 ≤ t) (hst : s + t = k + 1)
    (hn : s ≤ n) (hd : t ≤ d) : 1 ≤ cellCard n d k s t := by
  classical
  have hk : 1 ≤ k := by omega
  have hrow : ∀ x : Fin k, catRow s x.1 < n := by
    intro x; simp only [catRow]; split_ifs <;> omega
  have hcol : ∀ x : Fin k, catCol s x.1 < d := by
    intro x; have := x.2; simp only [catCol]; split_ifs <;> omega
  set i : Fin k → Fin n := fun x => ⟨catRow s x.1, hrow x⟩ with hi
  set j : Fin k → Fin d := fun x => ⟨catCol s x.1, hcol x⟩ with hj
  have hival : ∀ x : Fin k, (i x).1 = catRow s x.1 := fun x => rfl
  have hjval : ∀ x : Fin k, (j x).1 = catCol s x.1 := fun x => rfl
  have h1 : ∀ x : Fin k, ∃ y : Fin k, i y = i x ∧ j (cycSucc y) = j x := by
    intro x
    have hxk := x.2
    by_cases hA : x.1 + 2 ≤ s
    · refine ⟨x, rfl, Fin.ext ?_⟩
      rw [hjval, hjval, cycSucc_val, Nat.mod_eq_of_lt (by omega)]
      simp only [catCol]
      split_ifs <;> omega
    · by_cases hB : x.1 + 1 ≤ s
      · refine ⟨⟨k - 1, by omega⟩, Fin.ext ?_, Fin.ext ?_⟩
        · rw [hival, hival]
          simp only [catRow]
          split_ifs <;> omega
        · rw [hjval, hjval, cycSucc_val, Fin.val_mk,
            show (k - 1 + 1) % k = 0 from by rw [show k - 1 + 1 = k from by omega, Nat.mod_self]]
          simp only [catCol]
          split_ifs <;> omega
      · refine ⟨⟨x.1 - 1, by omega⟩, Fin.ext ?_, Fin.ext ?_⟩
        · rw [hival, hival]
          simp only [catRow]
          split_ifs <;> omega
        · rw [hjval, hjval, cycSucc_val, Fin.val_mk,
            show (x.1 - 1 + 1) % k = x.1 % k from by rw [show x.1 - 1 + 1 = x.1 from by omega],
            Nat.mod_eq_of_lt hxk]
  have h2 : ∀ x : Fin k, ∃ y : Fin k, i y = i x ∧ j y = j (cycSucc x) := by
    intro x
    have hxk := x.2
    by_cases hA : x.1 + 2 ≤ s
    · refine ⟨x, rfl, Fin.ext ?_⟩
      rw [hjval, hjval, cycSucc_val, Nat.mod_eq_of_lt (by omega)]
      simp only [catCol]
      split_ifs <;> omega
    · by_cases hL : x.1 + 1 = k
      · refine ⟨⟨s - 1, by omega⟩, Fin.ext ?_, Fin.ext ?_⟩
        · rw [hival, hival]
          simp only [catRow]
          split_ifs <;> omega
        · rw [hjval, hjval, cycSucc_val, Fin.val_mk,
            show (x.1 + 1) % k = 0 from by rw [hL, Nat.mod_self]]
          simp only [catCol]
          split_ifs <;> omega
      · refine ⟨⟨x.1 + 1, by omega⟩, Fin.ext ?_, Fin.ext ?_⟩
        · rw [hival, hival]
          simp only [catRow]
          split_ifs <;> omega
        · rw [hjval, hjval, cycSucc_val, Fin.val_mk, Nat.mod_eq_of_lt (by omega)]
  have hrows : (univ.image i).card = s := by
    refine card_image_val i (fun x => ?_) (fun u hu => ?_)
    · rw [hival]; simp only [catRow]; split_ifs <;> omega
    · rcases Nat.eq_zero_or_pos u with rfl | hu0
      · exact ⟨⟨s - 1, by omega⟩, by
          rw [hival]; simp only [catRow]; split_ifs <;> omega⟩
      · exact ⟨⟨u - 1, by omega⟩, by
          rw [hival]; simp only [catRow]; split_ifs <;> omega⟩
  have hcols : (univ.image j).card = t := by
    refine card_image_val j (fun x => ?_) (fun u hu => ?_)
    · rw [hjval]; have := x.2; simp only [catCol]; split_ifs <;> omega
    · rcases Nat.eq_zero_or_pos u with rfl | hu0
      · exact ⟨⟨0, by omega⟩, by
          rw [hjval]; simp only [catCol]; split_ifs <;> omega⟩
      · exact ⟨⟨s + u - 1, by omega⟩, by
          rw [hjval]; simp only [catCol]; split_ifs <;> omega⟩
  rw [cellCard]
  exact Finset.card_pos.mpr ⟨(i, j), mem_cellSet.mpr ⟨noSingle_of_families h1 h2, hrows, hcols⟩⟩

/-- L1: the class with `v` vertices is the disjoint union of its `(rows, columns)` cells. -/
theorem clsCard_eq_sum_cellCard (n d k v : ℕ) :
    clsCard n d k v = ∑ s ∈ range (v + 1), cellCard n d k s (v - s) := by
  classical
  have hfib : ∀ p ∈ walkSet n d k v, (univ.image p.1).card ∈ range (v + 1) := by
    intro p hp
    rw [mem_walkSet] at hp
    have hv := hp.2
    simp only [walkVerts] at hv
    rw [Finset.mem_range]
    omega
  rw [clsCard_eq_card_walkSet, Finset.card_eq_sum_card_fiberwise hfib]
  refine Finset.sum_congr rfl fun s hs => ?_
  have hsv : s ≤ v := Nat.lt_succ_iff.mp (Finset.mem_range.mp hs)
  rw [cellCard]
  congr 1
  ext p
  rw [Finset.mem_filter, mem_walkSet, mem_cellSet]
  simp only [walkVerts]
  constructor
  · rintro ⟨⟨hns, hv⟩, hA⟩
    exact ⟨hns, hA, by omega⟩
  · rintro ⟨hns, hA, hB⟩
    exact ⟨⟨hns, by omega⟩, hA⟩

/-- A class above the tree class is empty. -/
theorem clsCard_eq_zero_of_gt {n d k v : ℕ} (h : k + 1 < v) : clsCard n d k v = 0 := by
  classical
  rw [clsCard_eq_card_walkSet, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
  intro p hp
  rw [mem_walkSet] at hp
  have := walkVerts_le_of_noSingle hp.1
  omega

/-- The tree class is nonempty: the star walk that hangs `k` rows on one column lies in it. -/
theorem one_le_clsCard_tree {n d k : ℕ} (hk : 1 ≤ k) (hnd : 2 * k ≤ min n d) :
    1 ≤ clsCard n d k (k + 1) := by
  have hn : k ≤ n := by have := min_le_left n d; omega
  have hd : 1 ≤ d := by have := min_le_right n d; omega
  have hcell : 1 ≤ cellCard n d k k 1 := one_le_cellCard hk le_rfl (by omega) hn hd
  have hsum := Finset.single_le_sum (f := fun s => cellCard n d k s (k + 1 - s))
    (fun s _ => Nat.zero_le _) (Finset.mem_range.mpr (show k < k + 1 + 1 by omega))
  simp only [show k + 1 - k = 1 from by omega] at hsum
  rw [clsCard_eq_sum_cellCard]
  omega

end Edge
end StackedSVD
