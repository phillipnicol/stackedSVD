/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Count

/-! # Stage 3, unit X, part 2: the step code of a walk

The upper half of the Furedi-Komlos count (`notes/x8_plan.md`, sections 1 to 5; Furedi and
Komlos 1981, page 237). A closed walk of length `2 k` on `Fin n × Fin d` is read one step at
a time. A step is innovative when it reaches a vertex for the first time, forced when it
leaves the current vertex along the one entry of odd prefix count (the Eulerian parity rule),
and bad otherwise.

`code` bundles six parts: the innovative even steps, the innovative odd steps, the bad steps,
the endpoint index of every bad step in first-visit order, and the rows and the columns in
first-visit order. `code_injective` proves that the code determines the walk.

`card_codeP_le` bounds the bad steps of the cell `(s, t)` by `4 D`, with the vertex
deficiency `D = k + 1 - (s + t)`. Every step opens its entry or closes it, so a bad step is a
bad open or a bad close. A bad open gives one credit to the vertex it reaches; a bad close
takes one debt from the vertex it leaves. The invariant `credit_invariant` holds the debt of
a vertex at or below its credit, so the bad closes are at most the bad opens.

`cellCard_le_code'` is the corollary the assembly reads: the cell `(s, t)` holds at most
`C(k, s) * C(k, t - 1) * ((4 D + 1) * ((2 k) ^ 8) ^ D)` times the labels
`n.descFactorial s * d.descFactorial t`. -/

/-! ## Part U: the upper chain (U1, U3, U4, U5) -/

/-! ## X8 unit U1: the step types of a walk, and the innovative-step counts

Sections 1 and 2 of `notes/x8_plan.md`. A closed walk of length `2 k` is the pair
`(i : Fin k → Fin n, j : Fin k → Fin d)`; its vertex sequence is `u_(2t) = column (j t)`,
`u_(2t+1) = row (i t)`, cyclic. Step `s : Fin (2 k)` runs from `u_s` to `u_(s+1)` and
traverses `(i t, j t)` at `s = 2 t` and `(i t, j (cycSucc t))` at `s = 2 t + 1`, so the `2 k`
steps are exactly the `2 k` factors of `walkMult` (`walkMult_eq_card_stepEntry`).

Each step is innovative (`isInnov`: the endpoint is new), forced (`isForced`: the endpoint is
old and `oddSet` has exactly one entry at the current vertex, the one the step traverses) or
bad (`isBad`). `card_innov_even` counts one innovative even step per distinct row and
`card_innov_odd` one innovative odd step per distinct column other than the root `j 0`. -/

open Finset

namespace StackedSVD
namespace Edge

/-! ## 1. The vertex sequence and the entry of a step -/

/-- The walk index `⌊s / 2⌋` of the step `s`: the steps `2 t` and `2 t + 1` both read the
row `i t`. -/
def walkIdx {k : ℕ} (s : Fin (2 * k)) : Fin k := ⟨s.val / 2, by have := s.isLt; omega⟩

/-- The even step `2 t`, from the column `j t` to the row `i t`. -/
def evenStep {k : ℕ} (t : Fin k) : Fin (2 * k) := ⟨2 * t.val, by have := t.isLt; omega⟩

/-- The odd step `2 t + 1`, from the row `i t` to the column `j (cycSucc t)`. -/
def oddStep {k : ℕ} (t : Fin k) : Fin (2 * k) := ⟨2 * t.val + 1, by have := t.isLt; omega⟩

@[simp] theorem walkIdx_val {k : ℕ} (s : Fin (2 * k)) : (walkIdx s).val = s.val / 2 := rfl
@[simp] theorem evenStep_val {k : ℕ} (t : Fin k) : (evenStep t).val = 2 * t.val := rfl
@[simp] theorem oddStep_val {k : ℕ} (t : Fin k) : (oddStep t).val = 2 * t.val + 1 := rfl

@[simp] theorem walkIdx_evenStep {k : ℕ} (t : Fin k) : walkIdx (evenStep t) = t :=
  Fin.ext (by simp)

@[simp] theorem walkIdx_oddStep {k : ℕ} (t : Fin k) : walkIdx (oddStep t) = t :=
  Fin.ext (by simp; omega)

/-- An even step is `2 t` for its own walk index. -/
theorem evenStep_walkIdx {k : ℕ} {s : Fin (2 * k)} (h : s.val % 2 = 0) :
    evenStep (walkIdx s) = s := Fin.ext (by simp; omega)

/-- An odd step is `2 t + 1` for its own walk index. -/
theorem oddStep_walkIdx {k : ℕ} {s : Fin (2 * k)} (h : s.val % 2 = 1) :
    oddStep (walkIdx s) = s := Fin.ext (by simp; omega)

/-- The `s`-th vertex of the walk: `u_(2t)` is the column `j t`, `u_(2t+1)` is the row `i t`. -/
def vert {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) :
    Fin n ⊕ Fin d :=
  if s.val % 2 = 0 then Sum.inr (j (walkIdx s)) else Sum.inl (i (walkIdx s))

/-- The vertex at an even step is a column. -/
theorem vert_of_even {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) {s : Fin (2 * k)}
    (h : s.val % 2 = 0) : vert i j s = Sum.inr (j (walkIdx s)) := if_pos h

/-- The vertex at an odd step is a row. -/
theorem vert_of_odd {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) {s : Fin (2 * k)}
    (h : s.val % 2 = 1) : vert i j s = Sum.inl (i (walkIdx s)) := if_neg (by omega)

@[simp] theorem vert_evenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (t : Fin k) :
    vert i j (evenStep t) = Sum.inr (j t) := by
  rw [vert_of_even i j (by simp), walkIdx_evenStep]

@[simp] theorem vert_oddStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (t : Fin k) :
    vert i j (oddStep t) = Sum.inl (i t) := by
  rw [vert_of_odd i j (by simp), walkIdx_oddStep]

/-- The cyclic successor of an even step is the odd step of the same walk index. -/
@[simp] theorem cycSucc_evenStep {k : ℕ} (t : Fin k) : cycSucc (evenStep t) = oddStep t := by
  have ht := t.isLt
  refine Fin.ext ?_
  change (2 * t.val + 1) % (2 * k) = 2 * t.val + 1
  exact Nat.mod_eq_of_lt (by omega)

/-- The cyclic successor of an odd step is the even step of the next walk index. -/
@[simp] theorem cycSucc_oddStep {k : ℕ} (t : Fin k) :
    cycSucc (oddStep t) = evenStep (cycSucc t) := by
  refine Fin.ext ?_
  change (2 * t.val + 1 + 1) % (2 * k) = 2 * ((t.val + 1) % k)
  have h : 2 * t.val + 1 + 1 = 2 * (t.val + 1) := by ring
  rw [h, Nat.mul_mod_mul_left]

/-- The entry of `Fin n × Fin d` that the step `s` traverses: `(i t, j t)` at `s = 2 t` and
`(i t, j (cycSucc t))` at `s = 2 t + 1`. -/
def stepEntry {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) :
    Fin n × Fin d :=
  if s.val % 2 = 0 then (i (walkIdx s), j (walkIdx s))
  else (i (walkIdx s), j (cycSucc (walkIdx s)))

@[simp] theorem stepEntry_evenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (t : Fin k) : stepEntry i j (evenStep t) = (i t, j t) := by
  rw [stepEntry, if_pos (by simp), walkIdx_evenStep]

@[simp] theorem stepEntry_oddStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (t : Fin k) : stepEntry i j (oddStep t) = (i t, j (cycSucc t)) := by
  rw [stepEntry, if_neg (by simp), walkIdx_oddStep]

/-! ## 2. The bridge to `walkMult` -/

/-- The even steps that traverse `e` are the walk indices `t` with `(i t, j t) = e`. -/
private theorem card_even_stepEntry {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (e : Fin n × Fin d) :
    (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 0)).card
      = ({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k)).card := by
  refine Finset.card_bij' (fun s _ => walkIdx s) (fun t _ => evenStep t) ?_ ?_ ?_ ?_
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs ⊢
    rw [← evenStep_walkIdx hs.2, stepEntry_evenStep] at hs
    exact ⟨congrArg Prod.fst hs.1, congrArg Prod.snd hs.1⟩
  · intro t ht
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ht ⊢
    exact ⟨by rw [stepEntry_evenStep]; exact Prod.ext ht.1 ht.2, by simp⟩
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    exact evenStep_walkIdx hs.2
  · intro t _; exact walkIdx_evenStep t

/-- The odd steps that traverse `e` are the walk indices `t` with `(i t, j (cycSucc t)) = e`. -/
private theorem card_odd_stepEntry {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (e : Fin n × Fin d) :
    (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 1)).card
      = ({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k)).card := by
  refine Finset.card_bij' (fun s _ => walkIdx s) (fun t _ => oddStep t) ?_ ?_ ?_ ?_
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs ⊢
    rw [← oddStep_walkIdx hs.2, stepEntry_oddStep] at hs
    exact ⟨congrArg Prod.fst hs.1, congrArg Prod.snd hs.1⟩
  · intro t ht
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ht ⊢
    refine ⟨by rw [stepEntry_oddStep]; exact Prod.ext ht.1 ht.2, by simp⟩
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    exact oddStep_walkIdx hs.2
  · intro t _; exact walkIdx_oddStep t

/-- The bridge to `Count.lean`: the multiplicity of an entry is the number of steps that
traverse it. The `2 k` steps are exactly the `2 k` factors of `walkMult`. -/
theorem walkMult_eq_card_stepEntry {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (e : Fin n × Fin d) :
    walkMult i j e = (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e)).card := by
  have hunion : (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e))
      = (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 0))
        ∪ (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 1)) := by
    ext x
    simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · intro h
      have h2 : x.val % 2 = 0 ∨ x.val % 2 = 1 := by omega
      rcases h2 with h2 | h2
      · exact Or.inl ⟨h, h2⟩
      · exact Or.inr ⟨h, h2⟩
    · rintro (⟨h, -⟩ | ⟨h, -⟩) <;> exact h
  have hdisj : Disjoint
      (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 0))
      (univ.filter (fun s : Fin (2 * k) => stepEntry i j s = e ∧ s.val % 2 = 1)) := by
    rw [Finset.disjoint_left]
    intro x hx hx'
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx hx'
    have h0 := hx.2
    have h1 := hx'.2
    omega
  rw [hunion, Finset.card_union_of_disjoint hdisj, card_even_stepEntry, card_odd_stepEntry]
  rfl

/-! ## 3. The odd set, and the three step types -/

/-- The number of steps among `0, ..., s - 1` that traverse the entry `e`. -/
def prefixMult {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ)
    (e : Fin n × Fin d) : ℕ :=
  (univ.filter (fun r : Fin (2 * k) => r.val < s ∧ stepEntry i j r = e)).card

/-- `O_s`, the entries traversed an odd number of times by the steps `0, ..., s - 1`. The
Eulerian parity rule reads this set at the current vertex. -/
def oddSet {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ) :
    Finset (Fin n × Fin d) :=
  univ.filter (fun e => Odd (prefixMult i j s e))

/-- Membership in the odd set. -/
theorem mem_oddSet {k n d : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d} {s : ℕ}
    {e : Fin n × Fin d} : e ∈ oddSet i j s ↔ Odd (prefixMult i j s e) := by
  simp [oddSet]

/-- The vertex `v` is an endpoint of the entry `e`. -/
def Incident {n d : ℕ} (v : Fin n ⊕ Fin d) (e : Fin n × Fin d) : Prop :=
  v = Sum.inl e.1 ∨ v = Sum.inr e.2

instance decidableIncident {n d : ℕ} (v : Fin n ⊕ Fin d) : DecidablePred (Incident v) :=
  fun e => inferInstanceAs (Decidable (v = Sum.inl e.1 ∨ v = Sum.inr e.2))

/-- Step `s` is **innovative**: its endpoint `u_(s+1)` is not among `u_0, ..., u_s`. -/
def isInnov {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) : Prop :=
  ∀ r : Fin (2 * k), r.val ≤ s.val → vert i j (cycSucc s) ≠ vert i j r

instance decidableIsInnov {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    DecidablePred (isInnov i j) := fun s =>
  inferInstanceAs (Decidable (∀ r : Fin (2 * k), r.val ≤ s.val →
    vert i j (cycSucc s) ≠ vert i j r))

/-- Step `s` is **forced**: its endpoint is old, and the odd set `O_s` has exactly one entry
at the current vertex `u_s`, namely the entry that step `s` traverses. -/
def isForced {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) : Prop :=
  ¬ isInnov i j s ∧
    (oddSet i j s.val).filter (Incident (vert i j s)) = {stepEntry i j s}

instance decidableIsForced {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    DecidablePred (isForced i j) := fun s =>
  inferInstanceAs (Decidable (¬ isInnov i j s ∧
    (oddSet i j s.val).filter (Incident (vert i j s)) = {stepEntry i j s}))

/-- Step `s` is **bad**: neither innovative nor forced. -/
def isBad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) : Prop :=
  ¬ isInnov i j s ∧ ¬ isForced i j s

instance decidableIsBad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    DecidablePred (isBad i j) := fun s =>
  inferInstanceAs (Decidable (¬ isInnov i j s ∧ ¬ isForced i j s))

/-- Every step is innovative, forced or bad. -/
theorem isInnov_or_isForced_or_isBad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) : isInnov i j s ∨ isForced i j s ∨ isBad i j s := by
  by_cases h : isInnov i j s
  · exact Or.inl h
  · exact Or.inr (if hf : isForced i j s then Or.inl hf else Or.inr ⟨h, hf⟩)

/-! ## 4. The first-visit order -/

/-- The visited vertices in first-visit order, without repetition. `List.dedup` keeps the
last occurrence, so the list is reversed on both sides to keep the first one. -/
def firstVisitList {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    List (Fin n ⊕ Fin d) :=
  ((List.ofFn (vert i j)).reverse.dedup).reverse

/-- The index of a vertex in the first-visit order. The decoder of section 2 of the plan
reads the endpoint of a bad step through this index. This is `List.idxOf` of the
`DecidableEq`-derived `BEq`, written with `List.findIdx` so that no `LawfulBEq` instance on
`Fin n ⊕ Fin d` is needed. -/
def firstVisitIdx {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (v : Fin n ⊕ Fin d) : ℕ := (firstVisitList i j).findIdx (fun x => decide (x = v))

/-- The first-visit list has no repetition. -/
theorem firstVisitList_nodup {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (firstVisitList i j).Nodup :=
  List.nodup_reverse.mpr (List.nodup_dedup _)

/-- The first-visit list holds exactly the visited vertices. -/
theorem mem_firstVisitList {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (v : Fin n ⊕ Fin d) :
    v ∈ firstVisitList i j ↔ ∃ s : Fin (2 * k), vert i j s = v := by
  rw [firstVisitList, List.mem_reverse, List.mem_dedup, List.mem_reverse, List.mem_ofFn]

/-- A visited vertex has an index inside the first-visit list. -/
theorem firstVisitIdx_lt {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) :
    firstVisitIdx i j (vert i j s) < (firstVisitList i j).length :=
  List.findIdx_lt_length_of_exists
    ⟨vert i j s, (mem_firstVisitList i j _).mpr ⟨s, rfl⟩, by simp⟩

/-! ## 5. The innovative-step counts -/

/-- The indices at which `f` takes a value for the first time. -/
def firstOcc {k : ℕ} {α : Type*} [DecidableEq α] (f : Fin k → α) : Finset (Fin k) :=
  univ.filter (fun t => ∀ r : Fin k, r.val < t.val → f r ≠ f t)

/-- Every value of `f` has a first occurrence. -/
private theorem exists_first_occ {k : ℕ} {α : Type*} [DecidableEq α] (f : Fin k → α)
    (a : Fin k) : ∃ m : Fin k, f m = f a ∧ m ∈ firstOcc f := by
  classical
  have hne : (univ.filter (fun r : Fin k => f r = f a)).Nonempty := ⟨a, by simp⟩
  have hmem := Finset.min'_mem _ hne
  have hfa : f ((univ.filter (fun r : Fin k => f r = f a)).min' hne) = f a := by simpa using hmem
  refine ⟨_, hfa, ?_⟩
  simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and]
  intro r hr hfr
  have hrmem : r ∈ univ.filter (fun r : Fin k => f r = f a) := by
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    rw [hfr, hfa]
  have hle : ((univ.filter (fun r : Fin k => f r = f a)).min' hne).val ≤ r.val :=
    Finset.min'_le _ r hrmem
  omega

/-- One first occurrence per distinct value. -/
theorem card_firstOcc {k : ℕ} {α : Type*} [DecidableEq α] (f : Fin k → α) :
    (firstOcc f).card = (univ.image f).card := by
  refine Finset.card_bij (fun t _ => f t) ?_ ?_ ?_
  · intro a _; exact Finset.mem_image_of_mem f (Finset.mem_univ a)
  · intro a1 h1 a2 h2 heq
    simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and] at h1 h2
    by_contra hne
    rcases Nat.lt_trichotomy a1.val a2.val with h | h | h
    · exact h2 a1 h heq
    · exact hne (Fin.ext h)
    · exact h1 a2 h heq.symm
  · intro b hb
    obtain ⟨a, -, rfl⟩ := Finset.mem_image.mp hb
    obtain ⟨m, hm, hmem⟩ := exists_first_occ f a
    exact ⟨m, hmem, hm⟩

/-- An even step is innovative exactly when its walk index is the first occurrence of its
row label. -/
theorem isInnov_evenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (t : Fin k) :
    isInnov i j (evenStep t) ↔ ∀ r : Fin k, r.val < t.val → i r ≠ i t := by
  rw [isInnov, cycSucc_evenStep, vert_oddStep]
  constructor
  · intro h q hq hiq
    have ht := t.isLt
    exact h (oddStep q) (by simp; omega) (by rw [vert_oddStep]; exact congrArg Sum.inl hiq.symm)
  · intro h r hr
    have hrle : r.val ≤ 2 * t.val := hr
    by_cases hpar : r.val % 2 = 0
    · rw [vert_of_even i j hpar]; simp
    · rw [vert_of_odd i j (by omega)]
      have hlt : (walkIdx r).val < t.val := by simp only [walkIdx_val]; omega
      exact fun hc => h (walkIdx r) hlt (Sum.inl.inj hc).symm

/-- An odd step is innovative exactly when no column of index at most its walk index carries
the label of its endpoint. -/
theorem isInnov_oddStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (t : Fin k) :
    isInnov i j (oddStep t) ↔ ∀ r : Fin k, r.val ≤ t.val → j r ≠ j (cycSucc t) := by
  rw [isInnov, cycSucc_oddStep, vert_evenStep]
  constructor
  · intro h q hq hjq
    have ht := t.isLt
    exact h (evenStep q) (by simp; omega)
      (by rw [vert_evenStep]; exact congrArg Sum.inr hjq.symm)
  · intro h r hr
    have hrle : r.val ≤ 2 * t.val + 1 := hr
    by_cases hpar : r.val % 2 = 0
    · rw [vert_of_even i j hpar]
      have hle : (walkIdx r).val ≤ t.val := by simp only [walkIdx_val]; omega
      exact fun hc => h (walkIdx r) hle (Sum.inr.inj hc).symm
    · rw [vert_of_odd i j (by omega)]; simp

/-- Each distinct row is created by exactly one even step, its first visit. -/
theorem card_innov_even {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (univ.filter (fun s : Fin (2 * k) => Even s.val ∧ isInnov i j s)).card
      = (univ.image i).card := by
  rw [← card_firstOcc i]
  refine Finset.card_bij' (fun s _ => walkIdx s) (fun t _ => evenStep t) ?_ ?_ ?_ ?_
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    obtain ⟨hev, hin⟩ := hs
    have hpar : s.val % 2 = 0 := Nat.even_iff.mp hev
    rw [← evenStep_walkIdx hpar, isInnov_evenStep] at hin
    simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and]
    exact hin
  · intro t ht
    simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and] at ht
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨by simp, (isInnov_evenStep i j t).mpr ht⟩
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    exact evenStep_walkIdx (Nat.even_iff.mp hs.1)
  · intro t _; exact walkIdx_evenStep t

/-- A nonzero index of `Fin k` has a positive value. -/
private theorem one_le_val_of_ne_zero {k : ℕ} {hk : 1 ≤ k} {t : Fin k}
    (hne : t ≠ (⟨0, hk⟩ : Fin k)) : 1 ≤ t.val := by
  rcases Nat.eq_zero_or_pos t.val with h | h
  · exact absurd (Fin.ext h) hne
  · exact h

/-- The cyclic successor of `t - 1` is `t` at a positive `t`. -/
private theorem cycSucc_pred {k : ℕ} {t : Fin k} (htv : 1 ≤ t.val) :
    cycSucc (⟨t.val - 1, by have := t.isLt; omega⟩ : Fin k) = t := by
  have ht := t.isLt
  refine Fin.ext ?_
  change (t.val - 1 + 1) % k = t.val
  rw [show t.val - 1 + 1 = t.val from by omega]
  exact Nat.mod_eq_of_lt ht

/-- An innovative odd step does not close the cycle, so the index of its endpoint is the
successor of its walk index. -/
private theorem innov_oddStep_succ {k n d : ℕ} (hk : 1 ≤ k) {i : Fin k → Fin n}
    {j : Fin k → Fin d} {t : Fin k} (h : isInnov i j (oddStep t)) :
    (cycSucc t).val = t.val + 1 := by
  rw [isInnov_oddStep] at h
  have ht := t.isLt
  have htk : t.val + 1 < k := by
    rcases Nat.lt_or_ge (t.val + 1) k with hlt | hge
    · exact hlt
    · have hc : cycSucc t = (⟨0, hk⟩ : Fin k) := Fin.ext (by
        change (t.val + 1) % k = 0
        rw [show t.val + 1 = k from by omega, Nat.mod_self])
      have hbad := h ⟨0, hk⟩ (Nat.zero_le _)
      rw [hc] at hbad
      exact absurd rfl hbad
  change (t.val + 1) % k = t.val + 1
  exact Nat.mod_eq_of_lt htk

/-- Each distinct column other than the root `j 0` is created by exactly one odd step. -/
theorem card_innov_odd {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (univ.filter (fun s : Fin (2 * k) => Odd s.val ∧ isInnov i j s)).card + 1
      = (univ.image j).card := by
  have hzero : (⟨0, hk⟩ : Fin k) ∈ firstOcc j := by
    simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and]
    intro r hr; omega
  rw [← card_firstOcc j, ← Finset.card_erase_add_one hzero]
  congr 1
  refine Finset.card_bij' (fun s _ => cycSucc (walkIdx s))
    (fun t _ => oddStep ⟨t.val - 1, by have := t.isLt; omega⟩) ?_ ?_ ?_ ?_
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    have hin := hs.2
    rw [← oddStep_walkIdx (Nat.odd_iff.mp hs.1)] at hin
    have hsucc := innov_oddStep_succ hk hin
    rw [isInnov_oddStep] at hin
    refine Finset.mem_erase.mpr ⟨?_, ?_⟩
    · intro hc
      rw [hc] at hsucc
      have h0 : (0 : ℕ) = (walkIdx s).val + 1 := hsucc
      omega
    · simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and]
      exact fun r hr => hin r (by omega)
  · intro t ht
    obtain ⟨hne, hmem⟩ := Finset.mem_erase.mp ht
    simp only [firstOcc, Finset.mem_filter, Finset.mem_univ, true_and] at hmem
    have htv : 1 ≤ t.val := one_le_val_of_ne_zero hne
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    refine ⟨by simp, ?_⟩
    rw [isInnov_oddStep, cycSucc_pred htv]
    intro r hr
    have hr' : r.val ≤ t.val - 1 := hr
    exact hmem r (by omega)
  · intro s hs
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
    have hin := hs.2
    rw [← oddStep_walkIdx (Nat.odd_iff.mp hs.1)] at hin
    have hsucc := innov_oddStep_succ hk hin
    have hpar : s.val % 2 = 1 := Nat.odd_iff.mp hs.1
    have hw : (walkIdx s).val = s.val / 2 := rfl
    refine Fin.ext ?_
    change 2 * ((cycSucc (walkIdx s)).val - 1) + 1 = s.val
    omega
  · intro t ht
    have htv : 1 ≤ t.val := one_le_val_of_ne_zero (Finset.mem_erase.mp ht).1
    rw [walkIdx_oddStep]
    exact cycSucc_pred htv


/-! ## 6. The code of a walk (U4, deliverable 1)

The six components of the code of section 1 of `x8_plan.md`. `codeA` and `codeB` are the
innovative even and odd steps, `codeP` the bad steps, `codeE` the endpoint index of each bad
step in first-visit order, and `rowList`, `colList` the rows and the columns in first-visit
order. -/

/-- `A`: the innovative even steps, one per distinct row. -/
def codeA {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => Even s.val ∧ isInnov i j s)

/-- `B`: the innovative odd steps, one per distinct column other than the root `j 0`. -/
def codeB {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => Odd s.val ∧ isInnov i j s)

/-- `P`: the bad steps. -/
def codeP {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => isBad i j s)

/-- `e`: the index of the endpoint of a bad step in first-visit order, `0` at every other
step. -/
def codeE {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) : ℕ :=
  if isBad i j s then firstVisitIdx i j (vert i j (cycSucc s)) else 0

/-- The rows in first-visit order. -/
def rowList {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : List (Fin n) :=
  (firstVisitList i j).filterMap Sum.getLeft?

/-- The columns in first-visit order. -/
def colList {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : List (Fin d) :=
  (firstVisitList i j).filterMap Sum.getRight?

/-- The code of a walk: the innovative even steps `A`, the innovative odd steps `B`, the bad
steps `P`, the endpoint index of each bad step in first-visit order (`0` elsewhere), and the
rows and the columns in first-visit order. -/
def code {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    Finset (Fin (2 * k)) × Finset (Fin (2 * k)) × Finset (Fin (2 * k)) ×
      (Fin (2 * k) → ℕ) × List (Fin n) × List (Fin d) :=
  (codeA i j, codeB i j, codeP i j, codeE i j, rowList i j, colList i j)

section CodeParts

variable {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d}

theorem codeA_eq (h : code i j = code i' j') : codeA i j = codeA i' j' :=
  congrArg (fun c => c.1) h

theorem codeB_eq (h : code i j = code i' j') : codeB i j = codeB i' j' :=
  congrArg (fun c => c.2.1) h

theorem codeP_eq (h : code i j = code i' j') : codeP i j = codeP i' j' :=
  congrArg (fun c => c.2.2.1) h

theorem codeE_eq (h : code i j = code i' j') : codeE i j = codeE i' j' :=
  congrArg (fun c => c.2.2.2.1) h

theorem rowList_eq (h : code i j = code i' j') : rowList i j = rowList i' j' :=
  congrArg (fun c => c.2.2.2.2.1) h

theorem colList_eq (h : code i j = code i' j') : colList i j = colList i' j' :=
  congrArg (fun c => c.2.2.2.2.2) h

end CodeParts

/-! ## 7. The prefix-dependence lemmas

`stepEntry` at the step `s` reads `vert` at `s` and at `cycSucc s` only, so `oddSet i j m`
and the first-visit list of the prefix `u_0, ..., u_(m-1)` are functions of `vert` below `m`.
These are the workhorses of the induction of section 9. -/

/-- The endpoint of an even step is its row. -/
theorem vert_cycSucc_of_even {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : s.val % 2 = 0) :
    vert i j (cycSucc s) = Sum.inl (i (walkIdx s)) := by
  obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = evenStep t := ⟨walkIdx s, (evenStep_walkIdx h).symm⟩
  simp

/-- The endpoint of an odd step is the next column. -/
theorem vert_cycSucc_of_odd {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : s.val % 2 = 1) :
    vert i j (cycSucc s) = Sum.inr (j (cycSucc (walkIdx s))) := by
  obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = oddStep t := ⟨walkIdx s, (oddStep_walkIdx h).symm⟩
  simp

/-- The entry of a step is read from its two ends. -/
theorem stepEntry_eq_of_vert_eq {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d}
    {s : Fin (2 * k)} (h1 : vert i j s = vert i' j' s)
    (h2 : vert i j (cycSucc s) = vert i' j' (cycSucc s)) :
    stepEntry i j s = stepEntry i' j' s := by
  by_cases hpar : s.val % 2 = 0
  · obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = evenStep t := ⟨walkIdx s, (evenStep_walkIdx hpar).symm⟩
    rw [vert_evenStep, vert_evenStep] at h1
    rw [cycSucc_evenStep, vert_oddStep, vert_oddStep] at h2
    rw [stepEntry_evenStep, stepEntry_evenStep]
    exact Prod.ext (Sum.inl.inj h2) (Sum.inr.inj h1)
  · obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = oddStep t :=
      ⟨walkIdx s, (oddStep_walkIdx (by omega)).symm⟩
    rw [vert_oddStep, vert_oddStep] at h1
    rw [cycSucc_oddStep, vert_evenStep, vert_evenStep] at h2
    rw [stepEntry_oddStep, stepEntry_oddStep]
    exact Prod.ext (Sum.inl.inj h1) (Sum.inr.inj h2)

/-- The far end of a step is read from its entry. -/
theorem vert_cycSucc_eq_of_stepEntry_eq {k n d : ℕ} {i i' : Fin k → Fin n}
    {j j' : Fin k → Fin d} {s : Fin (2 * k)}
    (h : stepEntry i j s = stepEntry i' j' s) :
    vert i j (cycSucc s) = vert i' j' (cycSucc s) := by
  by_cases hpar : s.val % 2 = 0
  · obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = evenStep t := ⟨walkIdx s, (evenStep_walkIdx hpar).symm⟩
    rw [stepEntry_evenStep, stepEntry_evenStep] at h
    rw [cycSucc_evenStep, vert_oddStep, vert_oddStep]
    exact congrArg Sum.inl (congrArg Prod.fst h)
  · obtain ⟨t, rfl⟩ : ∃ t : Fin k, s = oddStep t :=
      ⟨walkIdx s, (oddStep_walkIdx (by omega)).symm⟩
    rw [stepEntry_oddStep, stepEntry_oddStep] at h
    rw [cycSucc_oddStep, vert_evenStep, vert_evenStep]
    exact congrArg Sum.inr (congrArg Prod.snd h)

/-- `oddSet` at `m` reads only the entries of the steps below `m`. -/
theorem oddSet_congr {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d} {m : ℕ}
    (hstep : ∀ p : Fin (2 * k), p.val < m → stepEntry i j p = stepEntry i' j' p) :
    oddSet i j m = oddSet i' j' m := by
  have hpm : ∀ e, prefixMult i j m e = prefixMult i' j' m e := by
    intro e
    unfold prefixMult
    congr 1
    ext p
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨h1, by rw [← hstep p h1]; exact h2⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1, by rw [hstep p h1]; exact h2⟩
  unfold oddSet
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, hpm e]

/-- The visited vertices of the prefix `u_0, ..., u_(m-1)`, in first-visit order. The decoder
of section 2 of the plan reads the endpoint of a bad step out of this list. -/
def visitListUpTo {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (m : ℕ) :
    List (Fin n ⊕ Fin d) :=
  ((((List.ofFn (vert i j)).take m).reverse).dedup).reverse

/-- The whole walk gives the whole first-visit list. -/
theorem visitListUpTo_self {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    visitListUpTo i j (2 * k) = firstVisitList i j := by
  unfold visitListUpTo firstVisitList
  rw [List.take_of_length_le (by simp)]

/-- The first-visit list of a prefix reads only `vert` below `m`. -/
theorem visitListUpTo_congr {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d} {m : ℕ}
    (hv : ∀ p : Fin (2 * k), p.val < m → vert i j p = vert i' j' p) :
    visitListUpTo i j m = visitListUpTo i' j' m := by
  have hteq : ((List.ofFn (vert i j)).take m) = ((List.ofFn (vert i' j')).take m) := by
    refine List.ext_getElem? ?_
    intro p
    rcases Nat.lt_or_ge p m with hlt | hge
    · rw [List.getElem?_take_of_lt hlt, List.getElem?_take_of_lt hlt]
      rcases Nat.lt_or_ge p (2 * k) with h2 | h2
      · rw [List.getElem?_eq_getElem (by simpa using h2),
            List.getElem?_eq_getElem (by simpa using h2)]
        simp only [List.getElem_ofFn]
        exact congrArg some (hv ⟨p, h2⟩ hlt)
      · rw [List.getElem?_eq_none (by simpa using h2), List.getElem?_eq_none (by simpa using h2)]
    · rw [List.getElem?_eq_none (by simp; omega), List.getElem?_eq_none (by simp; omega)]
  unfold visitListUpTo
  rw [hteq]

/-- The first-visit list of a prefix is a prefix of the first-visit list. -/
theorem visitListUpTo_prefix {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (m : ℕ) :
    visitListUpTo i j m <+: firstVisitList i j := by
  unfold visitListUpTo firstVisitList
  rw [List.reverse_prefix]
  have hsplit : (List.ofFn (vert i j)).reverse
      = ((List.ofFn (vert i j)).drop m).reverse ++ ((List.ofFn (vert i j)).take m).reverse := by
    conv_lhs => rw [← List.take_append_drop m (List.ofFn (vert i j))]
    rw [List.reverse_append]
  rw [hsplit, List.dedup_append]
  exact List.suffix_union_right _ _

/-- Every vertex of the prefix is in the first-visit list of the prefix. -/
theorem mem_visitListUpTo_of_lt {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) {m : ℕ}
    {p : Fin (2 * k)} (hp : p.val < m) : vert i j p ∈ visitListUpTo i j m := by
  unfold visitListUpTo
  rw [List.mem_reverse, List.mem_dedup, List.mem_reverse]
  refine List.mem_of_getElem? (i := p.val) ?_
  rw [List.getElem?_take_of_lt hp, List.getElem?_eq_getElem (by simp [p.isLt])]
  simp

/-- The first-visit list of a prefix holds exactly the vertices of the prefix. -/
theorem mem_visitListUpTo {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) {m : ℕ}
    {v : Fin n ⊕ Fin d} :
    v ∈ visitListUpTo i j m ↔ ∃ p : Fin (2 * k), p.val < m ∧ vert i j p = v := by
  constructor
  · intro h
    unfold visitListUpTo at h
    rw [List.mem_reverse, List.mem_dedup, List.mem_reverse, List.mem_iff_getElem?] at h
    obtain ⟨q, hq⟩ := h
    have hqm : q < m := by
      by_contra hc
      rw [List.getElem?_eq_none (by simp only [List.length_take, List.length_ofFn]; omega)] at hq
      exact absurd hq (by simp)
    rw [List.getElem?_take_of_lt hqm] at hq
    have hq2 : q < 2 * k := by
      by_contra hc
      rw [List.getElem?_eq_none (by simp only [List.length_ofFn]; omega)] at hq
      exact absurd hq (by simp)
    refine ⟨⟨q, hq2⟩, hqm, ?_⟩
    rw [List.getElem?_eq_getElem (by simp only [List.length_ofFn]; exact hq2)] at hq
    simpa using hq
  · rintro ⟨p, hp, rfl⟩
    exact mem_visitListUpTo_of_lt i j hp

/-- A new vertex is appended to the first-visit list of the prefix. -/
theorem visitListUpTo_succ_of_notMem {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {m : ℕ} (hm : m < 2 * k) (h : vert i j ⟨m, hm⟩ ∉ visitListUpTo i j m) :
    visitListUpTo i j (m + 1) = visitListUpTo i j m ++ [vert i j ⟨m, hm⟩] := by
  have hL : (List.ofFn (vert i j)).take (m + 1)
      = (List.ofFn (vert i j)).take m ++ [vert i j ⟨m, hm⟩] := by
    rw [List.take_add_one,
      List.getElem?_eq_getElem (by simp only [List.length_ofFn]; exact hm)]
    simp
  have hnot : vert i j ⟨m, hm⟩ ∉ ((List.ofFn (vert i j)).take m).reverse := by
    intro hc
    refine h ?_
    unfold visitListUpTo
    rw [List.mem_reverse, List.mem_dedup]
    exact hc
  unfold visitListUpTo
  rw [hL, List.reverse_append]
  simp only [List.reverse_cons, List.reverse_nil, List.nil_append, List.cons_append]
  rw [List.dedup_cons_of_notMem hnot, List.reverse_cons]

/-- An old vertex does not change the first-visit list of the prefix. -/
theorem visitListUpTo_succ_of_mem {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {m : ℕ} (hm : m < 2 * k) (h : vert i j ⟨m, hm⟩ ∈ visitListUpTo i j m) :
    visitListUpTo i j (m + 1) = visitListUpTo i j m := by
  have hL : (List.ofFn (vert i j)).take (m + 1)
      = (List.ofFn (vert i j)).take m ++ [vert i j ⟨m, hm⟩] := by
    rw [List.take_add_one,
      List.getElem?_eq_getElem (by simp only [List.length_ofFn]; exact hm)]
    simp
  have hmem : vert i j ⟨m, hm⟩ ∈ ((List.ofFn (vert i j)).take m).reverse := by
    rw [List.mem_reverse]
    have := h
    unfold visitListUpTo at this
    rw [List.mem_reverse, List.mem_dedup, List.mem_reverse] at this
    exact this
  unfold visitListUpTo
  rw [hL, List.reverse_append]
  simp only [List.reverse_cons, List.reverse_nil, List.nil_append, List.cons_append]
  rw [List.dedup_cons_of_mem hmem]

/-- An innovative step appends its endpoint to the first-visit list of the prefix, so the
whole first-visit list splits at that point. -/
theorem firstVisitList_split_innov {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hin : isInnov i j s) :
    ∃ R : List (Fin n ⊕ Fin d), firstVisitList i j
      = visitListUpTo i j (s.val + 1) ++ vert i j (cycSucc s) :: R := by
  have hs2 : s.val + 1 < 2 * k := by
    rcases Nat.lt_or_ge (s.val + 1) (2 * k) with hlt | hge
    · exact hlt
    · have hs := s.isLt
      have h0 : (0 : ℕ) < 2 * k := by omega
      have hc : cycSucc s = (⟨0, h0⟩ : Fin (2 * k)) := by
        refine Fin.ext ?_
        change (s.val + 1) % (2 * k) = 0
        rw [show s.val + 1 = 2 * k from by omega, Nat.mod_self]
      exact absurd (by rw [hc]) (hin ⟨0, h0⟩ (Nat.zero_le _))
  have hcyc : cycSucc s = (⟨s.val + 1, hs2⟩ : Fin (2 * k)) := by
    refine Fin.ext ?_
    change (s.val + 1) % (2 * k) = s.val + 1
    exact Nat.mod_eq_of_lt hs2
  have hnot : vert i j (⟨s.val + 1, hs2⟩ : Fin (2 * k)) ∉ visitListUpTo i j (s.val + 1) := by
    intro hc
    obtain ⟨p, hp, hpv⟩ := (mem_visitListUpTo i j).mp hc
    exact hin p (by omega) (by rw [hcyc]; exact hpv.symm)
  have hgrow := visitListUpTo_succ_of_notMem i j hs2 hnot
  obtain ⟨R, hR⟩ := visitListUpTo_prefix i j (s.val + 1 + 1)
  refine ⟨R, ?_⟩
  rw [← hR, hgrow, hcyc, List.append_assoc]
  rfl

/-- Reading one element out of a `filterMap` at a known split. -/
private theorem getElem?_filterMap_of_split {α β : Type*} {Q Pre R : List α} {v : α} {x : β}
    {f : α → Option β} (hQ : Q = Pre ++ v :: R) (hv : f v = some x) :
    (Q.filterMap f)[(Pre.filterMap f).length]? = some x := by
  rw [hQ, List.filterMap_append, List.filterMap_cons_some hv,
    List.getElem?_append_right (Nat.le_refl _), Nat.sub_self]
  simp

/-! ## 8. The three step types agree when the codes agree -/

section StepTypes

variable {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d}

/-- `A` and `B` name the innovative steps. -/
theorem isInnov_congr (hA : codeA i j = codeA i' j') (hB : codeB i j = codeB i' j')
    (s : Fin (2 * k)) : isInnov i j s ↔ isInnov i' j' s := by
  by_cases hpar : s.val % 2 = 0
  · have h1 : s ∈ codeA i j ↔ s ∈ codeA i' j' := by rw [hA]
    simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and] at h1
    have hev : Even s.val := Nat.even_iff.mpr hpar
    exact ⟨fun hi => (h1.mp ⟨hev, hi⟩).2, fun hi => (h1.mpr ⟨hev, hi⟩).2⟩
  · have h1 : s ∈ codeB i j ↔ s ∈ codeB i' j' := by rw [hB]
    simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and] at h1
    have hod : Odd s.val := Nat.odd_iff.mpr (by omega)
    exact ⟨fun hi => (h1.mp ⟨hod, hi⟩).2, fun hi => (h1.mpr ⟨hod, hi⟩).2⟩

/-- `P` names the bad steps. -/
theorem isBad_congr (hP : codeP i j = codeP i' j') (s : Fin (2 * k)) :
    isBad i j s ↔ isBad i' j' s := by
  have h1 : s ∈ codeP i j ↔ s ∈ codeP i' j' := by rw [hP]
  simpa only [codeP, Finset.mem_filter, Finset.mem_univ, true_and] using h1

end StepTypes

/-! ## 9. The three sub-lemmas of the induction

Each one says where the decoder of section 2 of the plan reads the endpoint of one step.
They are about one walk, not about a pair of walks; the induction of section 10 applies each
to both walks of the pair and compares. -/

/-- `List.dedup` keeps the last occurrence of each value, so it keeps the last element. -/
private theorem getLast?_dedup {α : Type*} [DecidableEq α] (l : List α) :
    l.dedup.getLast? = l.getLast? := by
  induction l with
  | nil => simp
  | cons a t ih =>
    rcases t with _ | ⟨b, u⟩
    · simp
    · rw [List.getLast?_cons_cons]
      by_cases hm : a ∈ (b :: u)
      · rw [List.dedup_cons_of_mem hm]; exact ih
      · rw [List.dedup_cons_of_notMem hm]
        have hne : (b :: u).dedup ≠ [] := by
          intro hc
          have hnil := (List.dedup_eq_nil (b :: u)).mp hc
          simp at hnil
        obtain ⟨c, v, hcv⟩ := List.exists_cons_of_ne_nil hne
        rw [hcv, List.getLast?_cons_cons, ← hcv]
        exact ih

/-- The first-visit list starts at the first vertex of the walk. -/
theorem head_firstVisitList {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) :
    (firstVisitList i j).head? = some (vert i j ⟨0, hk⟩) := by
  unfold firstVisitList
  rw [List.head?_reverse, getLast?_dedup, List.getLast?_reverse, List.head?_eq_getElem?,
    List.getElem?_eq_getElem (by simpa using hk)]
  simp

/-- Base: the first vertex of the walk is the first column of `colList`. -/
theorem head_colList {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (colList i j).head? = Sum.getRight? (vert i j ⟨0, hk⟩) := by
  have hv0 : vert i j (⟨0, hk⟩ : Fin (2 * k)) = Sum.inr (j (walkIdx ⟨0, hk⟩)) :=
    vert_of_even i j rfl
  have hhead := head_firstVisitList hk i j
  obtain ⟨L, hL⟩ : ∃ L, firstVisitList i j = vert i j ⟨0, hk⟩ :: L := by
    rcases hFV : firstVisitList i j with _ | ⟨a, L⟩
    · rw [hFV] at hhead; simp at hhead
    · refine ⟨L, ?_⟩
      rw [hFV] at hhead
      simp only [List.head?_cons, Option.some.injEq] at hhead
      rw [hhead]
  unfold colList
  rw [hL, hv0]
  simp

/-- COUNTING, still open: one innovative even step per distinct row of the prefix
`u_0, ..., u_s`. The induction is on `s`: `visitListUpTo` grows by one vertex exactly at an
innovative step (`visitListUpTo_succ_of_notMem`), and the new vertex is a row exactly when
the step is even. -/
private theorem card_filter_lt_succ_of_mem {k : ℕ} {S : Finset (Fin (2 * k))} {m : ℕ}
    (hm : m < 2 * k) (h : (⟨m, hm⟩ : Fin (2 * k)) ∈ S) :
    (S.filter (fun q => q.val < m + 1)).card = (S.filter (fun q => q.val < m)).card + 1 := by
  have hset : S.filter (fun q => q.val < m + 1)
      = insert (⟨m, hm⟩ : Fin (2 * k)) (S.filter (fun q => q.val < m)) := by
    ext q
    simp only [Finset.mem_filter, Finset.mem_insert]
    constructor
    · rintro ⟨hq, hlt⟩
      rcases Nat.lt_or_ge q.val m with h1 | h1
      · exact Or.inr ⟨hq, h1⟩
      · exact Or.inl (Fin.ext (show q.val = m by omega))
    · rintro (rfl | ⟨hq, h1⟩)
      · exact ⟨h, Nat.lt_succ_self m⟩
      · exact ⟨hq, by omega⟩
  rw [hset, Finset.card_insert_of_notMem (by simp)]

private theorem card_filter_lt_succ_of_notMem {k : ℕ} {S : Finset (Fin (2 * k))} {m : ℕ}
    (hm : m < 2 * k) (h : (⟨m, hm⟩ : Fin (2 * k)) ∉ S) :
    (S.filter (fun q => q.val < m + 1)).card = (S.filter (fun q => q.val < m)).card := by
  congr 1
  ext q
  simp only [Finset.mem_filter]
  constructor
  · rintro ⟨hq, hlt⟩
    refine ⟨hq, ?_⟩
    rcases Nat.lt_or_ge q.val m with h1 | h1
    · exact h1
    · exact absurd (show q = (⟨m, hm⟩ : Fin (2 * k)) from Fin.ext (show q.val = m by omega))
        fun hc => h (hc ▸ hq)
  · rintro ⟨hq, h1⟩
    exact ⟨hq, by omega⟩

/-- The counting step: the first-visit list of the prefix `u_0, ..., u_m` holds one row per
innovative even step below `m` and one column per innovative odd step below `m`, plus the
root column `u_0`. -/
theorem rowCol_count_upTo {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    ∀ m : ℕ, m < 2 * k →
      ((visitListUpTo i j (m + 1)).filterMap Sum.getLeft?).length
          = ((codeA i j).filter (fun q => q.val < m)).card ∧
        ((visitListUpTo i j (m + 1)).filterMap Sum.getRight?).length
          = ((codeB i j).filter (fun q => q.val < m)).card + 1 := by
  intro m
  induction m with
  | zero =>
    intro hm
    have h0 : visitListUpTo i j 0 = ([] : List (Fin n ⊕ Fin d)) := by simp [visitListUpTo]
    have h1 : visitListUpTo i j 1 = [vert i j ⟨0, hm⟩] := by
      have hstep := visitListUpTo_succ_of_notMem i j hm (by rw [h0]; simp)
      rw [hstep, h0, List.nil_append]
    have hv : vert i j (⟨0, hm⟩ : Fin (2 * k)) = Sum.inr (j (walkIdx ⟨0, hm⟩)) :=
      vert_of_even i j rfl
    rw [h1, hv]
    exact ⟨by simp, by simp⟩
  | succ m ihm =>
    intro hm1
    have hm : m < 2 * k := by omega
    obtain ⟨ihA, ihB⟩ := ihm hm
    have hcyc : cycSucc (⟨m, hm⟩ : Fin (2 * k)) = (⟨m + 1, hm1⟩ : Fin (2 * k)) := by
      refine Fin.ext ?_
      change (m + 1) % (2 * k) = m + 1
      exact Nat.mod_eq_of_lt hm1
    have hinn : isInnov i j (⟨m, hm⟩ : Fin (2 * k))
        ↔ vert i j (⟨m + 1, hm1⟩ : Fin (2 * k)) ∉ visitListUpTo i j (m + 1) := by
      rw [mem_visitListUpTo]
      constructor
      · rintro h ⟨p, hp, hpv⟩
        exact h p (Nat.lt_succ_iff.mp hp) (by rw [hcyc]; exact hpv.symm)
      · intro h p hp hpv
        rw [hcyc] at hpv
        exact h ⟨p, Nat.lt_succ_of_le hp, hpv.symm⟩
    by_cases hnew : vert i j (⟨m + 1, hm1⟩ : Fin (2 * k)) ∈ visitListUpTo i j (m + 1)
    · have hnot : ¬ isInnov i j (⟨m, hm⟩ : Fin (2 * k)) := fun h => (hinn.mp h) hnew
      have hgrow : visitListUpTo i j (m + 1 + 1) = visitListUpTo i j (m + 1) :=
        visitListUpTo_succ_of_mem i j hm1 hnew
      have hA : (⟨m, hm⟩ : Fin (2 * k)) ∉ codeA i j := by
        simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and]
        exact fun h => hnot h.2
      have hB : (⟨m, hm⟩ : Fin (2 * k)) ∉ codeB i j := by
        simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and]
        exact fun h => hnot h.2
      rw [hgrow, card_filter_lt_succ_of_notMem hm hA, card_filter_lt_succ_of_notMem hm hB]
      exact ⟨ihA, ihB⟩
    · have hin : isInnov i j (⟨m, hm⟩ : Fin (2 * k)) := hinn.mpr hnew
      have hgrow : visitListUpTo i j (m + 1 + 1)
          = visitListUpTo i j (m + 1) ++ [vert i j ⟨m + 1, hm1⟩] :=
        visitListUpTo_succ_of_notMem i j hm1 hnew
      by_cases hpar : (m + 1) % 2 = 0
      · have hv : vert i j (⟨m + 1, hm1⟩ : Fin (2 * k))
            = Sum.inr (j (walkIdx ⟨m + 1, hm1⟩)) := vert_of_even i j hpar
        have hAnot : (⟨m, hm⟩ : Fin (2 * k)) ∉ codeA i j := by
          simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and]
          rintro ⟨hev, -⟩
          rw [Nat.even_iff] at hev
          omega
        have hBmem : (⟨m, hm⟩ : Fin (2 * k)) ∈ codeB i j := by
          simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and]
          exact ⟨Nat.odd_iff.mpr (by omega), hin⟩
        rw [hgrow, hv, card_filter_lt_succ_of_notMem hm hAnot,
          card_filter_lt_succ_of_mem hm hBmem]
        exact ⟨by simp [ihA], by simp [ihB]⟩
      · have hv : vert i j (⟨m + 1, hm1⟩ : Fin (2 * k))
            = Sum.inl (i (walkIdx ⟨m + 1, hm1⟩)) :=
          vert_of_odd i j (by change (m + 1) % 2 = 1; omega)
        have hAmem : (⟨m, hm⟩ : Fin (2 * k)) ∈ codeA i j := by
          simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and]
          exact ⟨Nat.even_iff.mpr (by omega), hin⟩
        have hBnot : (⟨m, hm⟩ : Fin (2 * k)) ∉ codeB i j := by
          simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and]
          rintro ⟨hod, -⟩
          rw [Nat.odd_iff] at hod
          omega
        rw [hgrow, hv, card_filter_lt_succ_of_mem hm hAmem,
          card_filter_lt_succ_of_notMem hm hBnot]
        exact ⟨by simp [ihA], by simp [ihB]⟩

theorem card_codeA_lt_eq_rowCount {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) :
    ((codeA i j).filter (fun q => q.val < s.val)).card
      = ((visitListUpTo i j (s.val + 1)).filterMap Sum.getLeft?).length :=
  (rowCol_count_upTo i j s.val s.isLt).1.symm

/-- COUNTING, still open: one innovative odd step per distinct column of the prefix
`u_0, ..., u_s` other than the root `j 0`, which is `u_0`. -/
theorem card_codeB_lt_eq_colCount {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) :
    ((codeB i j).filter (fun q => q.val < s.val)).card + 1
      = ((visitListUpTo i j (s.val + 1)).filterMap Sum.getRight?).length :=
  (rowCol_count_upTo i j s.val s.isLt).2.symm

/-- Innovative even step: its endpoint is the row of rank `|A ∩ [0, s)|` in first-visit
order. -/
theorem rowList_getElem_innov_even {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hpar : s.val % 2 = 0) (hin : isInnov i j s) :
    (rowList i j)[((codeA i j).filter (fun q => q.val < s.val)).card]?
      = some (i (walkIdx s)) := by
  obtain ⟨R, hR⟩ := firstVisitList_split_innov i j hin
  rw [vert_cycSucc_of_even i j hpar] at hR
  rw [card_codeA_lt_eq_rowCount i j s]
  unfold rowList
  exact getElem?_filterMap_of_split hR rfl

/-- Innovative odd step: its endpoint is the column of rank `|B ∩ [0, s)| + 1` in first-visit
order, the `+ 1` for the root column `j 0`. -/
theorem colList_getElem_innov_odd {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hpar : s.val % 2 = 1) (hin : isInnov i j s) :
    (colList i j)[((codeB i j).filter (fun q => q.val < s.val)).card + 1]?
      = some (j (cycSucc (walkIdx s))) := by
  obtain ⟨R, hR⟩ := firstVisitList_split_innov i j hin
  rw [vert_cycSucc_of_odd i j hpar] at hR
  rw [card_codeB_lt_eq_colCount i j s]
  unfold colList
  exact getElem?_filterMap_of_split hR rfl

/-- Bad step: its endpoint is old, so the index `codeE i j s` of the endpoint in the whole
first-visit list already points inside the first-visit list of the prefix `u_0, ..., u_s`,
and it points at the endpoint. -/
theorem visitListUpTo_getElem_bad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hbad : isBad i j s) :
    (visitListUpTo i j (s.val + 1))[firstVisitIdx i j (vert i j (cycSucc s))]?
      = some (vert i j (cycSucc s)) := by
  set v := vert i j (cycSucc s) with hvdef
  set P := visitListUpTo i j (s.val + 1) with hPdef
  set pr : (Fin n ⊕ Fin d) → Bool := fun x => decide (x = v) with hprdef
  -- the endpoint of a bad step is old, so it is a vertex of the prefix
  have hold : ∃ r : Fin (2 * k), r.val ≤ s.val ∧ v = vert i j r := by
    by_contra hc
    exact hbad.1 (fun r hr hveq => hc ⟨r, hr, hveq⟩)
  obtain ⟨r, hr, hrv⟩ := hold
  have hmem : v ∈ P := by
    rw [hrv]
    exact mem_visitListUpTo_of_lt i j (by omega)
  have hlt : P.findIdx pr < P.length :=
    List.findIdx_lt_length_of_exists ⟨v, hmem, by simp [hprdef]⟩
  -- the index in the whole list is the index in the prefix
  obtain ⟨R, hR⟩ := visitListUpTo_prefix i j (s.val + 1)
  have hfind : (firstVisitList i j).findIdx pr = P.findIdx pr := by
    rw [← hR, List.findIdx_append, if_pos hlt]
  have hval : P[P.findIdx pr]'hlt = v := by
    have h1 : pr (P[P.findIdx pr]'hlt) = true := List.findIdx_getElem (w := hlt)
    simpa [hprdef] using h1
  change P[firstVisitIdx i j v]? = some v
  unfold firstVisitIdx
  rw [hfind, List.getElem?_eq_getElem hlt, hval]

/-! ## 10. Injectivity (U4, deliverable 2) -/

section Injective

variable {k n d : ℕ} {i i' : Fin k → Fin n} {j j' : Fin k → Fin d}

/-- One step of the decoder: the far end of the step `s` is a function of the code and of the
prefix `u_0, ..., u_s`. -/
theorem vert_cycSucc_eq_of_code_eq (h : code i j = code i' j') (s : Fin (2 * k))
    (hpre : ∀ p : Fin (2 * k), p.val ≤ s.val → vert i j p = vert i' j' p) :
    vert i j (cycSucc s) = vert i' j' (cycSucc s) := by
  have hA := codeA_eq h
  have hB := codeB_eq h
  have hP := codeP_eq h
  by_cases hin : isInnov i j s
  · -- innovative: the endpoint is the next unused label of its side
    have hin' : isInnov i' j' s := (isInnov_congr hA hB s).mp hin
    by_cases hpar : s.val % 2 = 0
    · have e1 := rowList_getElem_innov_even i j hpar hin
      have e2 := rowList_getElem_innov_even i' j' hpar hin'
      rw [rowList_eq h, hA] at e1
      have hrow : i (walkIdx s) = i' (walkIdx s) := Option.some.inj (e1.symm.trans e2)
      rw [vert_cycSucc_of_even i j hpar, vert_cycSucc_of_even i' j' hpar, hrow]
    · have hpar1 : s.val % 2 = 1 := by omega
      have e1 := colList_getElem_innov_odd i j hpar1 hin
      have e2 := colList_getElem_innov_odd i' j' hpar1 hin'
      rw [colList_eq h, hB] at e1
      have hcol : j (cycSucc (walkIdx s)) = j' (cycSucc (walkIdx s)) :=
        Option.some.inj (e1.symm.trans e2)
      rw [vert_cycSucc_of_odd i j hpar1, vert_cycSucc_of_odd i' j' hpar1, hcol]
  · by_cases hbad : isBad i j s
    · -- bad: the endpoint is the `e s`-th vertex of the first-visit list of the prefix
      have hbad' : isBad i' j' s := (isBad_congr hP s).mp hbad
      have hlist : visitListUpTo i j (s.val + 1) = visitListUpTo i' j' (s.val + 1) :=
        visitListUpTo_congr (fun p hp => hpre p (by omega))
      have hidx : firstVisitIdx i j (vert i j (cycSucc s))
          = firstVisitIdx i' j' (vert i' j' (cycSucc s)) := by
        have := congrFun (codeE_eq h) s
        simpa only [codeE, if_pos hbad, if_pos hbad'] using this
      have e1 := visitListUpTo_getElem_bad i j hbad
      have e2 := visitListUpTo_getElem_bad i' j' hbad'
      rw [hlist, hidx] at e1
      exact Option.some.inj (e1.symm.trans e2)
    · -- forced: the endpoint is the far end of the unique open entry at `u_s`
      have hfor : isForced i j s := by
        rcases isInnov_or_isForced_or_isBad i j s with h1 | h1 | h1
        · exact absurd h1 hin
        · exact h1
        · exact absurd h1 hbad
      have hfor' : isForced i' j' s := by
        rcases isInnov_or_isForced_or_isBad i' j' s with h1 | h1 | h1
        · exact absurd h1 ((isInnov_congr hA hB s).not.mp hin)
        · exact h1
        · exact absurd h1 ((isBad_congr hP s).not.mp hbad)
      have hstep : ∀ p : Fin (2 * k), p.val < s.val →
          stepEntry i j p = stepEntry i' j' p := by
        intro p hp
        refine stepEntry_eq_of_vert_eq (hpre p (by omega)) ?_
        have hcyc : (cycSucc p).val = p.val + 1 := by
          have hp2 := p.isLt
          change (p.val + 1) % (2 * k) = p.val + 1
          exact Nat.mod_eq_of_lt (by omega)
        exact hpre (cycSucc p) (by omega)
      have hos : oddSet i j s.val = oddSet i' j' s.val := oddSet_congr hstep
      have hv := hpre s le_rfl
      have hsing : ({stepEntry i j s} : Finset (Fin n × Fin d)) = {stepEntry i' j' s} := by
        rw [← hfor.2, ← hfor'.2, hos, hv]
      exact vert_cycSucc_eq_of_stepEntry_eq (Finset.singleton_injective hsing)

/-- The vertex sequence is a function of the code. -/
theorem vert_eq_of_code_eq (h : code i j = code i' j') (s : Fin (2 * k)) :
    vert i j s = vert i' j' s := by
  have key : ∀ m : ℕ, ∀ p : Fin (2 * k), p.val ≤ m → vert i j p = vert i' j' p := by
    intro m
    induction m with
    | zero =>
      intro p hp
      have hp0 : p.val = 0 := Nat.le_zero.mp hp
      have hk2 : 0 < 2 * k := by have := p.isLt; omega
      have hpe : p = (⟨0, hk2⟩ : Fin (2 * k)) := Fin.ext hp0
      subst hpe
      have hz : ((⟨0, hk2⟩ : Fin (2 * k))).val % 2 = 0 := rfl
      have h1 := head_colList hk2 i j
      have h2 := head_colList hk2 i' j'
      rw [colList_eq h] at h1
      have h3 : Sum.getRight? (vert i j (⟨0, hk2⟩ : Fin (2 * k)))
          = Sum.getRight? (vert i' j' (⟨0, hk2⟩ : Fin (2 * k))) := h1.symm.trans h2
      rw [vert_of_even i j hz, vert_of_even i' j' hz] at h3 ⊢
      exact congrArg Sum.inr (Option.some.inj (by simpa using h3))
    | succ q ihq =>
      intro p hp
      rcases Nat.lt_or_ge p.val (q + 1) with hlt | hge
      · exact ihq p (by omega)
      · have hpq : p.val = q + 1 := by omega
        have hq : q < 2 * k := by have := p.isLt; omega
        have hcyc : cycSucc (⟨q, hq⟩ : Fin (2 * k)) = p := by
          refine Fin.ext ?_
          change (q + 1) % (2 * k) = p.val
          rw [hpq]
          exact Nat.mod_eq_of_lt (by omega)
        rw [← hcyc]
        exact vert_cycSucc_eq_of_code_eq h ⟨q, hq⟩ (fun r hr => ihq r hr)
  exact key s.val s le_rfl

/-- **The code is injective.** Two walks with the same code are the same walk. -/
theorem code_injective (h : code i j = code i' j') : i = i' ∧ j = j' := by
  have hv := vert_eq_of_code_eq h
  constructor
  · funext t
    have ht := hv (oddStep t)
    rw [vert_oddStep, vert_oddStep] at ht
    exact Sum.inl.inj ht
  · funext t
    have ht := hv (evenStep t)
    rw [vert_evenStep, vert_evenStep] at ht
    exact Sum.inr.inj ht

end Injective


/-! ## 11. U5: the count of the codes

Section 6 of `notes/x8_plan.md`, row U5. The code map of section 6 is injective on a cell
(`code_injective`), so the cell is at most as large as the set of codes. The code lands in a
five-fold product: the innovative even steps (a subset of the `k` even steps of size `s`),
the innovative odd steps (a subset of the `k` odd steps of size `t - 1`), the pair (bad
steps, endpoint index), the rows in first-visit order and the columns in first-visit order.

The bad-step bound `#codeP ≤ 4 D` is unit U3 and is not in Lean yet, so it enters as the
hypothesis `hbad`. -/

section U5

/-- The `k` even steps. -/
def evenSteps (k : ℕ) : Finset (Fin (2 * k)) := univ.filter (fun x => Even x.val)

/-- The `k` odd steps. -/
def oddSteps (k : ℕ) : Finset (Fin (2 * k)) := univ.filter (fun x => Odd x.val)

theorem card_evenSteps (k : ℕ) : (evenSteps k).card = k := by
  have h : (evenSteps k).card = (univ : Finset (Fin k)).card := by
    refine Finset.card_bij' (fun x _ => walkIdx x) (fun t _ => evenStep t) ?_ ?_ ?_ ?_
    · intro x _; exact Finset.mem_univ _
    · intro y _
      simp only [evenSteps, Finset.mem_filter, Finset.mem_univ, true_and, evenStep_val,
        Nat.even_iff]
      omega
    · intro x hx
      simp only [evenSteps, Finset.mem_filter, Finset.mem_univ, true_and, Nat.even_iff] at hx
      exact evenStep_walkIdx hx
    · intro y _; exact walkIdx_evenStep y
  simpa using h

theorem card_oddSteps (k : ℕ) : (oddSteps k).card = k := by
  have h : (oddSteps k).card = (univ : Finset (Fin k)).card := by
    refine Finset.card_bij' (fun x _ => walkIdx x) (fun t _ => oddStep t) ?_ ?_ ?_ ?_
    · intro x _; exact Finset.mem_univ _
    · intro y _
      simp only [oddSteps, Finset.mem_filter, Finset.mem_univ, true_and, oddStep_val,
        Nat.odd_iff]
      omega
    · intro x hx
      simp only [oddSteps, Finset.mem_filter, Finset.mem_univ, true_and, Nat.odd_iff] at hx
      exact oddStep_walkIdx hx
    · intro y _; exact walkIdx_oddStep y
  simpa using h

/-- The pairs `(P, e)` a bad-step code can take: `P` has at most `4 D` steps, `e` is `0` off
`P` and at most `k` on `P`. -/
def codePE (k D : ℕ) : Finset (Finset (Fin (2 * k)) × (Fin (2 * k) → ℕ)) :=
  ((univ : Finset (Finset (Fin (2 * k)))).filter (fun P => P.card ≤ 4 * D)).biUnion
    (fun P => ({P} : Finset (Finset (Fin (2 * k)))) ×ˢ
      Fintype.piFinset (fun x => if x ∈ P then Finset.range (k + 1) else ({0} : Finset ℕ)))

/-- The number of bad-step codes. `Σ_{b ≤ 4 D} C(2k, b) (k+1)^b ≤ (4D+1) (2k)^(4D) (k+1)^(4D)`
and `k + 1 ≤ 2 k` for `k ≥ 1`. Open here (U5 step 3 of the plan). -/
theorem card_codePE (k D : ℕ) (hk : 1 ≤ k) :
    (codePE k D).card ≤ (4 * D + 1) * ((2 * k) ^ 8) ^ D := by
  set S : Finset (Finset (Fin (2 * k))) :=
    (univ : Finset (Finset (Fin (2 * k)))).filter (fun P => P.card ≤ 4 * D) with hS
  have hfib : ∀ P : Finset (Fin (2 * k)),
      (({P} : Finset (Finset (Fin (2 * k)))) ×ˢ
          Fintype.piFinset
            (fun x => if x ∈ P then Finset.range (k + 1) else ({0} : Finset ℕ))).card
        = (k + 1) ^ P.card := by
    intro P
    rw [Finset.card_product, Finset.card_singleton, one_mul, Fintype.card_piFinset]
    have hc : ∀ x : Fin (2 * k),
        (if x ∈ P then Finset.range (k + 1) else ({0} : Finset ℕ)).card
          = if x ∈ P then (k + 1) else 1 := by
      intro x; by_cases hx : x ∈ P <;> simp [hx]
    rw [Finset.prod_congr rfl (fun x _ => hc x), Finset.prod_ite_mem, Finset.univ_inter,
      Finset.prod_const]
  have h1 : (codePE k D).card ≤ ∑ P ∈ S, (k + 1) ^ P.card := by
    have hcp : codePE k D = S.biUnion (fun P => ({P} : Finset (Finset (Fin (2 * k)))) ×ˢ
        Fintype.piFinset
          (fun x => if x ∈ P then Finset.range (k + 1) else ({0} : Finset ℕ))) := by
      rw [hS]; rfl
    rw [hcp]
    exact le_trans Finset.card_biUnion_le
      (le_of_eq (Finset.sum_congr rfl (fun P _ => hfib P)))
  have h2 : ∑ P ∈ S, (k + 1) ^ P.card ≤ S.card * (k + 1) ^ (4 * D) := by
    calc ∑ P ∈ S, (k + 1) ^ P.card ≤ ∑ _P ∈ S, (k + 1) ^ (4 * D) := by
          refine Finset.sum_le_sum ?_
          intro P hP
          rw [hS, Finset.mem_filter] at hP
          exact Nat.pow_le_pow_right (by omega) hP.2
      _ = S.card * (k + 1) ^ (4 * D) := by rw [Finset.sum_const, smul_eq_mul]
  have h3 : S.card ≤ (4 * D + 1) * (2 * k) ^ (4 * D) := by
    have hsub : S ⊆ (Finset.range (4 * D + 1)).biUnion
        (fun b => Finset.powersetCard b (univ : Finset (Fin (2 * k)))) := by
      intro P hP
      rw [hS, Finset.mem_filter] at hP
      rw [Finset.mem_biUnion]
      exact ⟨P.card, Finset.mem_range.mpr (by omega),
        Finset.mem_powersetCard.mpr ⟨Finset.subset_univ _, rfl⟩⟩
    refine le_trans (Finset.card_le_card hsub) (le_trans Finset.card_biUnion_le ?_)
    have hb : ∀ b ∈ Finset.range (4 * D + 1),
        (Finset.powersetCard b (univ : Finset (Fin (2 * k)))).card ≤ (2 * k) ^ (4 * D) := by
      intro b hb
      have hble : b ≤ 4 * D := by have := Finset.mem_range.mp hb; omega
      rw [Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]
      exact le_trans (Nat.choose_le_pow _ _) (Nat.pow_le_pow_right (by omega) hble)
    refine le_trans (Finset.sum_le_sum hb) ?_
    rw [Finset.sum_const, Finset.card_range, smul_eq_mul]
  have h4 : (k + 1) ^ (4 * D) ≤ (2 * k) ^ (4 * D) := Nat.pow_le_pow_left (by omega) _
  have h5 : (2 * k) ^ (4 * D) * (2 * k) ^ (4 * D) = ((2 * k) ^ 8) ^ D := by
    rw [← pow_add, ← pow_mul]
    congr 1
    omega
  calc (codePE k D).card ≤ ∑ P ∈ S, (k + 1) ^ P.card := h1
    _ ≤ S.card * (k + 1) ^ (4 * D) := h2
    _ ≤ ((4 * D + 1) * (2 * k) ^ (4 * D)) * (2 * k) ^ (4 * D) := Nat.mul_le_mul h3 h4
    _ = (4 * D + 1) * ((2 * k) ^ (4 * D) * (2 * k) ^ (4 * D)) := by ring
    _ = (4 * D + 1) * ((2 * k) ^ 8) ^ D := by rw [h5]

/-- The number of injective maps `Fin a → Fin b`. -/
theorem card_injFun (a b : ℕ) :
    ((univ : Finset (Fin a → Fin b)).filter (fun g => Function.Injective g)).card
      = b.descFactorial a := by
  classical
  rw [← Fintype.card_subtype]
  rw [Fintype.card_congr (Equiv.subtypeInjectiveEquivEmbedding (Fin a) (Fin b))]
  simp [Fintype.card_embedding_eq]

/-- A `Nodup` list of length `m` gives an injective map out of `Fin m`. -/
private theorem list_getD_inj {α : Type*} {l : List α} (hl : l.Nodup) {m : ℕ}
    (hlen : l.length = m) (a : α) :
    Function.Injective (fun x : Fin m => l.getD x.val a) := by
  intro x y hxy
  have hx : x.val < l.length := by rw [hlen]; exact x.isLt
  have hy : y.val < l.length := by rw [hlen]; exact y.isLt
  simp only [List.getD_eq_getElem l a hx, List.getD_eq_getElem l a hy] at hxy
  exact Fin.ext (List.Nodup.getElem_inj_iff hl |>.mp hxy)

/-- Two lists of the same length with the same `getD` map are equal. -/
private theorem list_eq_of_getD_eq {α : Type*} {l l' : List α} {m : ℕ}
    (h : l.length = m) (h' : l'.length = m) (a : α)
    (heq : (fun x : Fin m => l.getD x.val a) = fun x : Fin m => l'.getD x.val a) :
    l = l' := by
  refine List.ext_getElem (by rw [h, h']) ?_
  intro q hq hq'
  have hlt : q < m := by rw [← h]; exact hq
  have hc := congrFun heq ⟨q, hlt⟩
  simpa only [List.getD_eq_getElem l a hq, List.getD_eq_getElem l' a hq'] using hc

/-- The row list has no repetition. -/
theorem rowList_nodup {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (rowList i j).Nodup := by
  refine List.Nodup.filterMap ?_ (firstVisitList_nodup i j)
  intro x y b hx hy
  cases x with
  | inl x => cases y with
    | inl y =>
      simp only [Sum.getLeft?, Option.mem_def, Option.some.injEq] at hx hy
      rw [hx, hy]
    | inr y => simp [Sum.getLeft?] at hy
  | inr x => simp [Sum.getLeft?] at hx

/-- The column list has no repetition. -/
theorem colList_nodup {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (colList i j).Nodup := by
  refine List.Nodup.filterMap ?_ (firstVisitList_nodup i j)
  intro x y b hx hy
  cases x with
  | inr x => cases y with
    | inr y =>
      simp only [Sum.getRight?, Option.mem_def, Option.some.injEq] at hx hy
      rw [hx, hy]
    | inl y => simp [Sum.getRight?] at hy
  | inl x => simp [Sum.getRight?] at hx

/-- The last step of a closed walk is never innovative: it returns to the root column. -/
theorem not_isInnov_last {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    ¬ isInnov i j ⟨2 * k - 1, by omega⟩ := by
  intro hin
  have hz : (2 * k - 1 : ℕ) < 2 * k := by omega
  have hcyc : cycSucc (⟨2 * k - 1, hz⟩ : Fin (2 * k)) = (⟨0, by omega⟩ : Fin (2 * k)) := by
    refine Fin.ext ?_
    rw [cycSucc_val]
    have : 2 * k - 1 + 1 = 2 * k := by omega
    rw [this, Nat.mod_self]
  exact hin ⟨0, by omega⟩ (by simp) (by rw [hcyc])

/-- The row list has one entry per distinct row. -/
theorem rowList_length {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (rowList i j).length = (univ.image i).card := by
  have hlt : 2 * k - 1 < 2 * k := by omega
  have h := (rowCol_count_upTo i j (2 * k - 1) hlt).1
  have hsucc : 2 * k - 1 + 1 = 2 * k := by omega
  rw [hsucc, visitListUpTo_self] at h
  rw [rowList, h, ← card_innov_even i j]
  congr 1
  refine Finset.filter_true_of_mem ?_
  intro x hx
  simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and, Nat.even_iff] at hx
  have := x.isLt
  omega

/-- The column list has one entry per distinct column. -/
theorem colList_length {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (colList i j).length = (univ.image j).card := by
  have hlt : 2 * k - 1 < 2 * k := by omega
  have h := (rowCol_count_upTo i j (2 * k - 1) hlt).2
  have hsucc : 2 * k - 1 + 1 = 2 * k := by omega
  rw [hsucc, visitListUpTo_self] at h
  rw [colList, h, ← card_innov_odd hk i j]
  congr 2
  refine Finset.filter_true_of_mem ?_
  intro x hx
  simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and] at hx
  rcases Nat.lt_or_ge x.val (2 * k - 1) with hx1 | hx1
  · exact hx1
  · exfalso
    have hxv : x.val = 2 * k - 1 := by have := x.isLt; omega
    refine not_isInnov_last hk i j ?_
    have : (⟨2 * k - 1, hlt⟩ : Fin (2 * k)) = x := Fin.ext hxv.symm
    rw [this]
    exact hx.2

/-- Every element of a list of a sum type is a row or a column. -/
private theorem length_sum_split {α β : Type*} (l : List (α ⊕ β)) :
    (l.filterMap Sum.getLeft?).length + (l.filterMap Sum.getRight?).length = l.length := by
  induction l with
  | nil => simp
  | cons a l ih =>
    cases a with
    | inl x =>
      have h1 : Sum.getLeft? (Sum.inl x : α ⊕ β) = some x := rfl
      have h2 : Sum.getRight? (Sum.inl x : α ⊕ β) = none := rfl
      simp only [List.filterMap_cons, h1, h2, List.length_cons]
      omega
    | inr y =>
      have h1 : Sum.getLeft? (Sum.inr y : α ⊕ β) = none := rfl
      have h2 : Sum.getRight? (Sum.inr y : α ⊕ β) = some y := rfl
      simp only [List.filterMap_cons, h1, h2, List.length_cons]
      omega

/-- The first-visit list holds one entry per distinct row and one per distinct column. -/
theorem firstVisitList_length {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) :
    (firstVisitList i j).length = (univ.image i).card + (univ.image j).card := by
  rw [← rowList_length hk i j, ← colList_length hk i j]
  exact (length_sum_split (firstVisitList i j)).symm

/-- The endpoint index of a bad step points inside the first-visit list. -/
theorem codeE_lt_length {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (hpos : 0 < (firstVisitList i j).length) (x : Fin (2 * k)) :
    codeE i j x < (firstVisitList i j).length := by
  unfold codeE
  split_ifs with h
  · exact firstVisitIdx_lt i j (cycSucc x)
  · exact hpos

/-- The endpoint index is `0` off the bad steps. -/
theorem codeE_eq_zero_of_notMem {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {x : Fin (2 * k)} (hx : x ∉ codeP i j) : codeE i j x = 0 := by
  have hb : ¬ isBad i j x := by
    simp only [codeP, Finset.mem_filter, Finset.mem_univ, true_and] at hx
    exact hx
  unfold codeE
  rw [if_neg hb]

/-- **U5**: the cell is at most as large as the set of codes. The bad-step bound of unit U3
(`#codeP ≤ 4 D`) is the hypothesis `hbad`; with it in place, this is
`cellCard_le_code` of `scripts/stage3_count/X8_skeleton.lean`. -/
theorem cellCard_le_code_of_bad {n d k s t : ℕ} (hk : 1 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hst : s + t ≤ k + 1)
    (hbad : ∀ p ∈ cellSet n d k s t, (codeP p.1 p.2).card ≤ 4 * (k + 1 - (s + t))) :
    cellCard n d k s t
      ≤ k.choose s * k.choose (t - 1) *
          ((4 * (k + 1 - (s + t)) + 1) * ((2 * k) ^ 8) ^ (k + 1 - (s + t))) *
          (n.descFactorial s * d.descFactorial t) := by
  classical
  rcases Nat.eq_zero_or_pos n with hn0 | hn
  · subst hn0
    have hempty : cellSet 0 d k s t = ∅ := by
      rw [Finset.eq_empty_iff_forall_notMem]
      intro p _
      exact (p.1 ⟨0, hk⟩).elim0
    simp [cellCard, hempty]
  rcases Nat.eq_zero_or_pos d with hd0 | hd
  · subst hd0
    have hempty : cellSet n 0 k s t = ∅ := by
      rw [Finset.eq_empty_iff_forall_notMem]
      intro p _
      exact (p.2 ⟨0, hk⟩).elim0
    simp [cellCard, hempty]
  set D := k + 1 - (s + t) with hDdef
  set r0 : Fin n := ⟨0, hn⟩ with hr0
  set c0 : Fin d := ⟨0, hd⟩ with hc0
  -- the two list lengths on the cell
  have rlen : ∀ q ∈ cellSet n d k s t, (rowList q.1 q.2).length = s := by
    intro q hq
    rw [rowList_length hk q.1 q.2, (mem_cellSet.mp hq).2.1]
  have clen : ∀ q ∈ cellSet n d k s t, (colList q.1 q.2).length = t := by
    intro q hq
    rw [colList_length hk q.1 q.2, (mem_cellSet.mp hq).2.2]
  set TGT := ((evenSteps k).powersetCard s) ×ˢ ((oddSteps k).powersetCard (t - 1)) ×ˢ
      (codePE k D) ×ˢ
      ((univ : Finset (Fin s → Fin n)).filter (fun g => Function.Injective g)) ×ˢ
      ((univ : Finset (Fin t → Fin d)).filter (fun g => Function.Injective g)) with hTGT
  have key : cellCard n d k s t ≤ TGT.card := by
    rw [cellCard]
    refine Finset.card_le_card_of_injOn
      (fun p => (codeA p.1 p.2, codeB p.1 p.2, (codeP p.1 p.2, codeE p.1 p.2),
        (fun x : Fin s => (rowList p.1 p.2).getD x.val r0),
        (fun x : Fin t => (colList p.1 p.2).getD x.val c0))) ?_ ?_
    · intro p hp
      rw [Finset.mem_coe] at hp
      obtain ⟨-, hrows, hcols⟩ := mem_cellSet.mp hp
      rw [Finset.mem_coe, hTGT]
      simp only [Finset.mem_product]
      refine ⟨?_, ?_, ?_, ?_, ?_⟩
      · rw [Finset.mem_powersetCard]
        refine ⟨?_, ?_⟩
        · intro x hx
          simp only [codeA, Finset.mem_filter, Finset.mem_univ, true_and] at hx
          simp only [evenSteps, Finset.mem_filter, Finset.mem_univ, true_and]
          exact hx.1
        · simp only [codeA]
          rw [card_innov_even p.1 p.2, hrows]
      · rw [Finset.mem_powersetCard]
        refine ⟨?_, ?_⟩
        · intro x hx
          simp only [codeB, Finset.mem_filter, Finset.mem_univ, true_and] at hx
          simp only [oddSteps, Finset.mem_filter, Finset.mem_univ, true_and]
          exact hx.1
        · have hb := card_innov_odd hk p.1 p.2
          rw [hcols] at hb
          simp only [codeB]
          omega
      · rw [hTGT] at *
        rw [codePE, Finset.mem_biUnion]
        refine ⟨codeP p.1 p.2, ?_, ?_⟩
        · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
          exact hbad p hp
        · rw [Finset.mem_product]
          refine ⟨Finset.mem_singleton_self _, ?_⟩
          rw [Fintype.mem_piFinset]
          intro x
          have hlen : (firstVisitList p.1 p.2).length = s + t := by
            rw [firstVisitList_length hk p.1 p.2, hrows, hcols]
          by_cases hxP : x ∈ codeP p.1 p.2
          · rw [if_pos hxP, Finset.mem_range]
            have hlt2 := codeE_lt_length p.1 p.2 (by omega) x
            change codeE p.1 p.2 x < k + 1
            omega
          · rw [if_neg hxP, Finset.mem_singleton]
            exact codeE_eq_zero_of_notMem p.1 p.2 hxP
      · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        exact list_getD_inj (rowList_nodup p.1 p.2) (rlen p hp) r0
      · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        exact list_getD_inj (colList_nodup p.1 p.2) (clen p hp) c0
    · intro p hp p' hp' hEq
      rw [Finset.mem_coe] at hp hp'
      have h1 := congrArg (fun c => c.1) hEq
      have h2 := congrArg (fun c => c.2.1) hEq
      have h3 := congrArg (fun c => c.2.2.1.1) hEq
      have h4 := congrArg (fun c => c.2.2.1.2) hEq
      have h5 := congrArg (fun c => c.2.2.2.1) hEq
      have h6 := congrArg (fun c => c.2.2.2.2) hEq
      simp only at h1 h2 h3 h4 h5 h6
      have hrl : rowList p.1 p.2 = rowList p'.1 p'.2 :=
        list_eq_of_getD_eq (rlen p hp) (rlen p' hp') r0 h5
      have hcl : colList p.1 p.2 = colList p'.1 p'.2 :=
        list_eq_of_getD_eq (clen p hp) (clen p' hp') c0 h6
      have hcode : code p.1 p.2 = code p'.1 p'.2 := by
        simp only [code, h1, h2, h3, h4, hrl, hcl]
      obtain ⟨e1, e2⟩ := code_injective hcode
      exact Prod.ext e1 e2
  have hcard : TGT.card = k.choose s * (k.choose (t - 1) *
      ((codePE k D).card * (n.descFactorial s * d.descFactorial t))) := by
    rw [hTGT, Finset.card_product, Finset.card_product, Finset.card_product,
      Finset.card_product, Finset.card_powersetCard, Finset.card_powersetCard,
      card_evenSteps, card_oddSteps, card_injFun, card_injFun]
  calc cellCard n d k s t ≤ TGT.card := key
    _ = k.choose s * (k.choose (t - 1) *
          ((codePE k D).card * (n.descFactorial s * d.descFactorial t))) := hcard
    _ ≤ k.choose s * (k.choose (t - 1) *
          (((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
            (n.descFactorial s * d.descFactorial t))) := by
        gcongr
        exact card_codePE k D hk
    _ = k.choose s * k.choose (t - 1) * ((4 * D + 1) * ((2 * k) ^ 8) ^ D) *
          (n.descFactorial s * d.descFactorial t) := by ring

end U5

/-! ## U3: the bad-step bound `#bad ≤ 4 D`

Section 6 of `scripts/stage3_count/x8_scratch/U3/u3_proof.md`. A step **opens** when its entry has
an even prefix count, and **closes** otherwise. Every innovative step opens and every forced step
closes, so the bad steps split into bad opens and bad closes. The bad opens are counted exactly
(`card_badOpen_eq`, Proposition 7) and the bad closes are bounded by them
(`card_badClose_le_card_badOpen`, Corollary 15) through the credit-debt invariant
(`credit_invariant`, Theorem 14). -/

section U3

/-- Step `s` **opens**: its entry is outside the odd set `O_s`, so the step puts the entry in.
A step that is not open **closes**: it takes its entry out of `O_s`. -/
def isOpenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : Fin (2 * k)) : Prop :=
  stepEntry i j s ∉ oddSet i j s.val

instance decidableIsOpenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    DecidablePred (isOpenStep i j) := fun s =>
  inferInstanceAs (Decidable (stepEntry i j s ∉ oddSet i j s.val))

/-- The **level** `l_z(s)` of the vertex `z` at time `s`: the number of entries of `O_s` that
have `z` as an endpoint. -/
def levelAt {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ)
    (z : Fin n ⊕ Fin d) : ℕ :=
  ((oddSet i j s).filter (Incident z)).card

/-- The vertex `u_s` of the closed walk at any natural time, `u_(2k) = u_0`. -/
def vertC {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ) :
    Fin n ⊕ Fin d :=
  vert i j ⟨s % (2 * k), Nat.mod_lt _ hk⟩

/-- The **credit** of `z` at time `s`: the bad opening steps before `s` with endpoint `z`. -/
def credAt {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ)
    (z : Fin n ⊕ Fin d) : ℕ :=
  (univ.filter (fun r : Fin (2 * k) =>
    r.val < s ∧ isBad i j r ∧ isOpenStep i j r ∧ vert i j (cycSucc r) = z)).card

/-- The **debt** of `z` at time `s`: the bad closing steps before `s` that leave `z`. -/
def debtAt {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (s : ℕ)
    (z : Fin n ⊕ Fin d) : ℕ :=
  (univ.filter (fun r : Fin (2 * k) =>
    r.val < s ∧ isBad i j r ∧ ¬ isOpenStep i j r ∧ vert i j r = z)).card

/-- The opening steps. -/
def openSteps {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => isOpenStep i j s)

/-- The closing steps. -/
def closeSteps {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => ¬ isOpenStep i j s)

/-- The bad opening steps. -/
def badOpen {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => isBad i j s ∧ isOpenStep i j s)

/-- The bad closing steps. -/
def badClose {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => isBad i j s ∧ ¬ isOpenStep i j s)

/-- The innovative steps. -/
def innovSteps {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin (2 * k)) :=
  univ.filter (fun s => isInnov i j s)

/-- `Oend`, the entries of odd multiplicity. -/
def oddEnd {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Finset (Fin n × Fin d) :=
  oddSet i j (2 * k)

/-! ### 1. One step of the odd set -/

/-- Counting the steps below `s + 1` that satisfy `P`: one more than below `s` exactly when
`s` itself satisfies `P`. -/
theorem card_lt_succ_split {k : ℕ} (P : Fin (2 * k) → Prop) [DecidablePred P]
    (s : Fin (2 * k)) :
    (univ.filter (fun r : Fin (2 * k) => r.val < s.val + 1 ∧ P r)).card
      = (univ.filter (fun r : Fin (2 * k) => r.val < s.val ∧ P r)).card
        + (if P s then 1 else 0) := by
  classical
  have hsplit :
      (univ.filter (fun r : Fin (2 * k) => r.val < s.val + 1 ∧ P r))
        = (univ.filter (fun r : Fin (2 * k) => r.val < s.val ∧ P r))
          ∪ (univ.filter (fun r : Fin (2 * k) => r = s ∧ P r)) := by
    ext r
    simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨h1, h2⟩
      rcases Nat.lt_or_ge r.val s.val with h | h
      · exact Or.inl ⟨h, h2⟩
      · exact Or.inr ⟨Fin.ext (by omega), h2⟩
    · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
      · exact ⟨by omega, h2⟩
      · subst h1; exact ⟨by omega, h2⟩
  have hdisj :
      Disjoint (univ.filter (fun r : Fin (2 * k) => r.val < s.val ∧ P r))
        (univ.filter (fun r : Fin (2 * k) => r = s ∧ P r)) := by
    rw [Finset.disjoint_left]
    intro a ha ha'
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha ha'
    have : a = s := ha'.1
    subst this
    omega
  have hlast :
      (univ.filter (fun r : Fin (2 * k) => r = s ∧ P r)).card = (if P s then 1 else 0) := by
    by_cases h : P s
    · rw [if_pos h]
      have hone : (univ.filter (fun r : Fin (2 * k) => r = s ∧ P r)) = {s} := by
        ext r
        simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
        exact ⟨fun h1 => h1.1, fun h1 => by subst h1; exact ⟨rfl, h⟩⟩
      rw [hone, Finset.card_singleton]
    · rw [if_neg h]
      have hnil : (univ.filter (fun r : Fin (2 * k) => r = s ∧ P r)) = ∅ := by
        rw [Finset.eq_empty_iff_forall_notMem]
        intro r hr
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
        obtain ⟨rfl, h2⟩ := hr
        exact h h2
      rw [hnil, Finset.card_empty]
  rw [hsplit, Finset.card_union_of_disjoint hdisj, hlast]

/-- The prefix count grows by one at a step that traverses the entry, and not at all
otherwise. -/
theorem prefixMult_succ {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) (e : Fin n × Fin d) :
    prefixMult i j (s.val + 1) e =
      prefixMult i j s.val e + (if stepEntry i j s = e then 1 else 0) := by
  classical
  have hsplit :
      (univ.filter (fun r : Fin (2 * k) => r.val < s.val + 1 ∧ stepEntry i j r = e))
        = (univ.filter (fun r : Fin (2 * k) => r.val < s.val ∧ stepEntry i j r = e))
          ∪ (univ.filter (fun r : Fin (2 * k) => r = s ∧ stepEntry i j r = e)) := by
    ext r
    simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨h1, h2⟩
      rcases Nat.lt_or_ge r.val s.val with h | h
      · exact Or.inl ⟨h, h2⟩
      · exact Or.inr ⟨Fin.ext (by omega), h2⟩
    · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
      · exact ⟨by omega, h2⟩
      · subst h1; exact ⟨by omega, h2⟩
  have hdisj :
      Disjoint (univ.filter (fun r : Fin (2 * k) => r.val < s.val ∧ stepEntry i j r = e))
        (univ.filter (fun r : Fin (2 * k) => r = s ∧ stepEntry i j r = e)) := by
    rw [Finset.disjoint_left]
    intro a ha ha'
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha ha'
    have : a = s := ha'.1
    subst this
    omega
  have hlast :
      (univ.filter (fun r : Fin (2 * k) => r = s ∧ stepEntry i j r = e)).card
        = (if stepEntry i j s = e then 1 else 0) := by
    by_cases h : stepEntry i j s = e
    · rw [if_pos h]
      have : (univ.filter (fun r : Fin (2 * k) => r = s ∧ stepEntry i j r = e)) = {s} := by
        ext r
        simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
        constructor
        · rintro ⟨h1, _⟩; exact h1
        · rintro rfl; exact ⟨rfl, h⟩
      rw [this, Finset.card_singleton]
    · rw [if_neg h]
      have : (univ.filter (fun r : Fin (2 * k) => r = s ∧ stepEntry i j r = e)) = ∅ := by
        rw [Finset.eq_empty_iff_forall_notMem]
        intro r hr
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
        obtain ⟨rfl, h2⟩ := hr
        exact h h2
      rw [this, Finset.card_empty]
  simp only [prefixMult]
  rw [hsplit, Finset.card_union_of_disjoint hdisj, hlast]

/-- An opening step inserts its entry into the odd set. -/
theorem oddSet_succ_of_open {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : isOpenStep i j s) :
    oddSet i j (s.val + 1) = insert (stepEntry i j s) (oddSet i j s.val) := by
  classical
  ext e
  simp only [mem_oddSet, Finset.mem_insert, prefixMult_succ]
  constructor
  · intro he
    by_cases hes : stepEntry i j s = e
    · exact Or.inl hes.symm
    · rw [if_neg hes, Nat.add_zero] at he
      exact Or.inr he
  · rintro (rfl | he)
    · rw [if_pos rfl]
      have := h
      rw [isOpenStep, mem_oddSet] at this
      rcases Nat.even_or_odd (prefixMult i j s.val (stepEntry i j s)) with he' | he'
      · rw [Nat.odd_iff]
        rw [Nat.even_iff] at he'
        omega
      · exact absurd he' this
    · by_cases hes : stepEntry i j s = e
      · exfalso
        apply h
        rw [mem_oddSet, hes]
        exact he
      · rw [if_neg hes, Nat.add_zero]
        exact he

/-- A closing step removes its entry from the odd set. -/
theorem oddSet_succ_of_close {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : ¬ isOpenStep i j s) :
    oddSet i j (s.val + 1) = (oddSet i j s.val).erase (stepEntry i j s) := by
  classical
  rw [isOpenStep, not_not, mem_oddSet, Nat.odd_iff] at h
  ext e
  simp only [mem_oddSet, Finset.mem_erase, prefixMult_succ]
  constructor
  · intro he
    by_cases hes : stepEntry i j s = e
    · exfalso
      subst hes
      rw [if_pos rfl, Nat.odd_iff] at he
      omega
    · rw [if_neg hes, Nat.add_zero] at he
      exact ⟨fun hc => hes hc.symm, he⟩
  · rintro ⟨hne, he⟩
    rw [if_neg (fun hc => hne hc.symm), Nat.add_zero]
    exact he

/-! ### 2. The two ends of a step -/

/-- The two ends of a step are distinct: one is a row, the other is a column. -/
theorem vert_ne_vert_cycSucc {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) : vert i j s ≠ vert i j (cycSucc s) := by
  rcases Nat.even_or_odd s.val with h | h
  · rw [Nat.even_iff] at h
    obtain ⟨t, rfl⟩ : ∃ t, s = evenStep t := ⟨walkIdx s, (evenStep_walkIdx h).symm⟩
    rw [vert_evenStep, cycSucc_evenStep, vert_oddStep]
    simp
  · rw [Nat.odd_iff] at h
    obtain ⟨t, rfl⟩ : ∃ t, s = oddStep t := ⟨walkIdx s, (oddStep_walkIdx h).symm⟩
    rw [vert_oddStep, cycSucc_oddStep, vert_evenStep]
    simp

/-- A vertex is an endpoint of the entry of a step exactly when it is one of the two ends of
the step. -/
theorem incident_stepEntry_iff {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) (z : Fin n ⊕ Fin d) :
    Incident z (stepEntry i j s) ↔ z = vert i j s ∨ z = vert i j (cycSucc s) := by
  rcases Nat.even_or_odd s.val with h | h
  · rw [Nat.even_iff] at h
    obtain ⟨t, rfl⟩ : ∃ t, s = evenStep t := ⟨walkIdx s, (evenStep_walkIdx h).symm⟩
    rw [stepEntry_evenStep, vert_evenStep, cycSucc_evenStep, vert_oddStep]
    exact Or.comm
  · rw [Nat.odd_iff] at h
    obtain ⟨t, rfl⟩ : ∃ t, s = oddStep t := ⟨walkIdx s, (oddStep_walkIdx h).symm⟩
    rw [stepEntry_oddStep, vert_oddStep, cycSucc_oddStep, vert_evenStep]
    exact Iff.rfl

/-- The current vertex is an endpoint of the entry of the step. -/
theorem incident_stepEntry_self {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) : Incident (vert i j s) (stepEntry i j s) :=
  (incident_stepEntry_iff i j s _).mpr (Or.inl rfl)

/-- The endpoint of the step is an endpoint of the entry of the step. -/
theorem incident_stepEntry_cycSucc {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) : Incident (vert i j (cycSucc s)) (stepEntry i j s) :=
  (incident_stepEntry_iff i j s _).mpr (Or.inr rfl)

/-- The walk vertex at a step index. -/
theorem vertC_val {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) : vertC hk i j s.val = vert i j s := by
  have h : (⟨s.val % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = s :=
    Fin.ext (Nat.mod_eq_of_lt s.isLt)
  rw [vertC, h]

/-- The root of the walk. -/
theorem vertC_zero {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    vertC hk i j 0 = vert i j ⟨0, hk⟩ := by
  have h : (⟨0 % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = ⟨0, hk⟩ :=
    Fin.ext (Nat.zero_mod _)
  rw [vertC, h]

/-- An entry with a positive prefix count is traversed by an earlier step. -/
theorem exists_stepEntry_of_prefixMult_pos {k n d : ℕ} (i : Fin k → Fin n)
    (j : Fin k → Fin d) {m : ℕ} {e : Fin n × Fin d} (h : 0 < prefixMult i j m e) :
    ∃ t : Fin (2 * k), t.val < m ∧ stepEntry i j t = e := by
  rw [prefixMult, Finset.card_pos] at h
  obtain ⟨t, ht⟩ := h
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ht
  exact ⟨t, ht.1, ht.2⟩

/-- A vertex that is an endpoint of an entry of `O_s` is already visited at time `s`. -/
theorem exists_vert_of_incident_oddSet {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} {z : Fin n ⊕ Fin d} {e : Fin n × Fin d}
    (he : e ∈ oddSet i j s.val) (hz : Incident z e) :
    ∃ r : Fin (2 * k), r.val ≤ s.val ∧ vert i j r = z := by
  rw [mem_oddSet, Nat.odd_iff] at he
  have hpos : 0 < prefixMult i j s.val e := by omega
  obtain ⟨t, htlt, hte⟩ := exists_stepEntry_of_prefixMult_pos i j hpos
  subst hte
  rw [incident_stepEntry_iff] at hz
  have hs2 : s.val < 2 * k := s.isLt
  rcases hz with rfl | rfl
  · exact ⟨t, by omega, rfl⟩
  · refine ⟨cycSucc t, ?_, rfl⟩
    have hc : (cycSucc t).val = (t.val + 1) % (2 * k) := rfl
    rw [hc, Nat.mod_eq_of_lt (by omega)]
    omega

/-! ### 3. One step of the level -/

/-- An opening step raises the level of each of its two ends by one. -/
theorem levelAt_succ_open_incident {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : isOpenStep i j s) {z : Fin n ⊕ Fin d}
    (hz : Incident z (stepEntry i j s)) :
    levelAt i j (s.val + 1) z = levelAt i j s.val z + 1 := by
  classical
  rw [levelAt, levelAt, oddSet_succ_of_open i j h, Finset.filter_insert, if_pos hz]
  rw [Finset.card_insert_of_notMem]
  simp only [Finset.mem_filter]
  intro hc
  exact h hc.1

/-- A closing step lowers the level of each of its two ends by one. -/
theorem levelAt_succ_close_incident {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : ¬ isOpenStep i j s) {z : Fin n ⊕ Fin d}
    (hz : Incident z (stepEntry i j s)) :
    levelAt i j (s.val + 1) z + 1 = levelAt i j s.val z := by
  classical
  rw [isOpenStep, not_not] at h
  rw [levelAt, levelAt, oddSet_succ_of_close i j (by rw [isOpenStep, not_not]; exact h),
    Finset.filter_erase]
  rw [Finset.card_erase_of_mem (Finset.mem_filter.mpr ⟨h, hz⟩)]
  have : 1 ≤ ((oddSet i j s.val).filter (Incident z)).card :=
    Finset.card_pos.mpr ⟨_, Finset.mem_filter.mpr ⟨h, hz⟩⟩
  omega

/-- A step that misses a vertex does not change its level. -/
theorem levelAt_succ_of_not_incident {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : Fin (2 * k)) {z : Fin n ⊕ Fin d} (hz : ¬ Incident z (stepEntry i j s)) :
    levelAt i j (s.val + 1) z = levelAt i j s.val z := by
  classical
  by_cases h : isOpenStep i j s
  · rw [levelAt, levelAt, oddSet_succ_of_open i j h, Finset.filter_insert, if_neg hz]
  · rw [levelAt, levelAt, oddSet_succ_of_close i j h, Finset.filter_erase,
      Finset.erase_eq_of_notMem]
    simp only [Finset.mem_filter]
    intro hc
    exact hz hc.2

/-! ### 4. The parity of the level (Lemma 1) -/

/-- **Lemma 1, the parity of the level.** `l_z(s) = [z = u_0] + [z = u_s] (mod 2)`.

Proof: Induction on `s` from `levelAt_succ_open_incident`, `levelAt_succ_close_incident` and
`levelAt_succ_of_not_incident`, with the case split `z = u_s`, `z = u_(s+1)`, neither. The
step needs `vert_ne_vert_cycSucc` (the two ends differ, so no vertex is both), and
`vertC hk i j (s + 1) = vert i j (cycSucc ⟨s, _⟩)` for `s < 2 k`. -/
theorem level_parity {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : ℕ) (hs : s ≤ 2 * k) (z : Fin n ⊕ Fin d) :
    levelAt i j s z % 2 =
      ((if z = vertC hk i j 0 then 1 else 0) + (if z = vertC hk i j s then 1 else 0)) % 2 := by
  classical
  revert hs
  induction s with
  | zero =>
      intro _
      have h0 : levelAt i j 0 z = 0 := by
        rw [levelAt, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
        intro e he
        rw [Finset.mem_filter, mem_oddSet] at he
        have hz : prefixMult i j 0 e = 0 := by
          rw [prefixMult, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
          intro r hr
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
          omega
        rw [hz, Nat.odd_iff] at he
        omega
      rw [h0]
      by_cases hz0 : z = vertC hk i j 0 <;> simp [hz0]
  | succ m ih =>
      intro hs
      have hmlt : m < 2 * k := by omega
      have ihm := ih (by omega)
      set σ : Fin (2 * k) := ⟨m, hmlt⟩ with hσ
      have hσv : σ.val = m := rfl
      have hvm : vertC hk i j m = vert i j σ := by
        have hfin : (⟨m % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = σ :=
          Fin.ext (Nat.mod_eq_of_lt hmlt)
        rw [vertC, hfin]
      have hvm1 : vertC hk i j (m + 1) = vert i j (cycSucc σ) := by
        have hfin : (⟨(m + 1) % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = cycSucc σ := rfl
        rw [vertC, hfin]
      have hne := vert_ne_vert_cycSucc i j σ
      set a : ℕ := (if z = vertC hk i j 0 then 1 else 0) with ha
      by_cases hinc : Incident z (stepEntry i j σ)
      · have hflip : levelAt i j (m + 1) z % 2 = (levelAt i j m z + 1) % 2 := by
          by_cases hop : isOpenStep i j σ
          · have hh := levelAt_succ_open_incident i j hop hinc
            rw [hσv] at hh
            rw [hh]
          · have hh := levelAt_succ_close_incident i j hop hinc
            rw [hσv] at hh
            omega
        rcases (incident_stepEntry_iff i j σ z).mp hinc with hA | hB
        · have hz2 : z ≠ vert i j (cycSucc σ) := by rw [hA]; exact hne
          rw [hvm, if_pos hA] at ihm
          rw [hvm1, if_neg hz2]
          omega
        · have hz1 : z ≠ vert i j σ := by rw [hB]; exact fun hc => hne hc.symm
          rw [hvm, if_neg hz1] at ihm
          rw [hvm1, if_pos hB]
          omega
      · have hsame : levelAt i j (m + 1) z = levelAt i j m z := by
          have hh := levelAt_succ_of_not_incident i j σ hinc
          rw [hσv] at hh
          exact hh
        have hz1 : z ≠ vert i j σ := fun hc =>
          hinc ((incident_stepEntry_iff i j σ z).mpr (Or.inl hc))
        have hz2 : z ≠ vert i j (cycSucc σ) := fun hc =>
          hinc ((incident_stepEntry_iff i j σ z).mpr (Or.inr hc))
        rw [hvm, if_neg hz1] at ihm
        rw [hvm1, if_neg hz2]
        omega

/-- Corollary 1(a): at a nonroot current vertex the level is odd, so at least one. -/
theorem one_le_level_self {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) (s : ℕ) (hs : s ≤ 2 * k) (hz : vertC hk i j s ≠ vertC hk i j 0) :
    1 ≤ levelAt i j s (vertC hk i j s) := by
  have h := level_parity hk i j s hs (vertC hk i j s)
  rw [if_neg hz, if_pos rfl] at h
  omega

/-! ### 5. The step types against opening and closing (Lemmas 2, 3 and 4) -/

/-- **Lemma 2(a).** Every innovative step opens.

Proof: If the entry were in `O_s` its prefix count is at least one, so an earlier step `t < s`
traverses it; the endpoint `u_(s+1)` is then one of `u_t`, `u_(t+1)`, both at index at most
`s`, against `isInnov`. -/
theorem innov_isOpenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : isInnov i j s) : isOpenStep i j s := by
  intro hmem
  obtain ⟨r, hrle, hrv⟩ :=
    exists_vert_of_incident_oddSet i j hmem (incident_stepEntry_cycSucc i j s)
  exact h r hrle hrv.symm

/-- **Lemma 2(b).** Every forced step closes. -/
theorem forced_not_isOpenStep {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (h : isForced i j s) : ¬ isOpenStep i j s := by
  intro hopen
  apply hopen
  have := h.2
  have hmem : stepEntry i j s ∈ (oddSet i j s.val).filter (Incident (vert i j s)) := by
    rw [this]
    exact Finset.mem_singleton_self _
  exact (Finset.mem_filter.mp hmem).1

/-- **Lemma 3, one half.** A closing step whose current vertex has level at least two is
bad.

Proof: The step is not innovative by `innov_isOpenStep`. It is not forced: the filter of `O_s`
at `u_s` has at least two elements, so it is not the singleton `{e_s}`
(`Finset.card_le_one` or `Finset.card_eq_one`). -/
theorem close_level_ge_two_isBad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hclose : ¬ isOpenStep i j s)
    (hlv : 2 ≤ levelAt i j s.val (vert i j s)) : isBad i j s := by
  refine ⟨fun hinnov => hclose (innov_isOpenStep i j hinnov), fun hforced => ?_⟩
  rw [levelAt, hforced.2, Finset.card_singleton] at hlv
  omega

/-- **Lemma 3, the other half.** A closing step whose current vertex has level one is
forced, so it is not bad.

Proof: `Finset.card_eq_one` turns the level into a singleton `{e}`; the entry `e_s` lies in
that filter because the step closes and `u_s` is one of its ends, so `e = e_s`. -/
theorem close_level_one_not_isBad {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    {s : Fin (2 * k)} (hclose : ¬ isOpenStep i j s)
    (hlv : levelAt i j s.val (vert i j s) = 1) : ¬ isBad i j s := by
  intro hbad
  apply hbad.2
  refine ⟨hbad.1, ?_⟩
  rw [levelAt, Finset.card_eq_one] at hlv
  obtain ⟨e, he⟩ := hlv
  have hmem : stepEntry i j s ∈ (oddSet i j s.val).filter (Incident (vert i j s)) :=
    Finset.mem_filter.mpr ⟨not_not.mp hclose, incident_stepEntry_self i j s⟩
  rw [he, Finset.mem_singleton] at hmem
  rw [he, hmem]

/-- **Lemma 4.** An opening step whose endpoint is the root or already has a positive level
is bad.

Proof: It is enough to show that the step is not innovative. If `u_(s+1) = u_0` that is
immediate. If the level of `u_(s+1)` is positive, an entry of `O_s` at `u_(s+1)` has a
positive prefix count, so an earlier step traverses it and `u_(s+1)` is already visited. -/
theorem open_arrival_isBad {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) {s : Fin (2 * k)} (hopen : isOpenStep i j s)
    (harr : vert i j (cycSucc s) = vertC hk i j 0 ∨ 1 ≤ levelAt i j s.val (vert i j (cycSucc s))) :
    isBad i j s := by
  refine ⟨?_, fun hf => forced_not_isOpenStep i j hf hopen⟩
  intro hinnov
  rcases harr with hroot | hlv
  · exact hinnov ⟨0, hk⟩ (Nat.zero_le _) (by rw [hroot, vertC_zero])
  · have hlv' : 0 < ((oddSet i j s.val).filter (Incident (vert i j (cycSucc s)))).card := by
      rw [levelAt] at hlv; omega
    rw [Finset.card_pos] at hlv'
    obtain ⟨e, hemem⟩ := hlv'
    rw [Finset.mem_filter] at hemem
    obtain ⟨r, hrle, hrv⟩ := exists_vert_of_incident_oddSet i j hemem.1 hemem.2
    exact hinnov r hrle hrv.symm

/-! ### 6. The counts (Lemmas 5 and 6, Propositions 7 and 8) -/

/-- Every step opens or closes. -/
theorem card_openSteps_add_card_closeSteps {k n d : ℕ} (i : Fin k → Fin n)
    (j : Fin k → Fin d) : (openSteps i j).card + (closeSteps i j).card = 2 * k := by
  classical
  rw [openSteps, closeSteps, Finset.card_filter_add_card_filter_not,
    Finset.card_univ, Fintype.card_fin]

/-- **Lemma 5.** `#open - #close = |Oend|`, in the form without subtraction.

Proof: The telescoping sum of `|O_s|` over `Finset.range (2 k)`: `|O_(s+1)| = |O_s| + 1` on an
opening step and `|O_(s+1)| + 1 = |O_s|` on a closing step, from `oddSet_succ_of_open` and
`oddSet_succ_of_close`. Induction on the number of steps with the statement
`|O_m| + #{closing steps below m} = #{opening steps below m}`. -/
theorem card_closeSteps_add_card_oddEnd {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (closeSteps i j).card + (oddEnd i j).card = (openSteps i j).card := by
  classical
  have key : ∀ m : ℕ, m ≤ 2 * k →
      (univ.filter (fun r : Fin (2 * k) => r.val < m ∧ ¬ isOpenStep i j r)).card
        + (oddSet i j m).card
        = (univ.filter (fun r : Fin (2 * k) => r.val < m ∧ isOpenStep i j r)).card := by
    intro m
    induction m with
    | zero =>
        intro _
        have h0 : oddSet i j 0 = (∅ : Finset (Fin n × Fin d)) := by
          rw [Finset.eq_empty_iff_forall_notMem]
          intro e he
          rw [mem_oddSet] at he
          have hz : prefixMult i j 0 e = 0 := by
            rw [prefixMult, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
            intro r hr
            simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
            omega
          rw [hz, Nat.odd_iff] at he
          omega
        have he1 : (univ.filter (fun r : Fin (2 * k) => r.val < 0 ∧ ¬ isOpenStep i j r))
            = ∅ := by
          rw [Finset.eq_empty_iff_forall_notMem]
          intro r hr
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
          omega
        have he2 : (univ.filter (fun r : Fin (2 * k) => r.val < 0 ∧ isOpenStep i j r))
            = ∅ := by
          rw [Finset.eq_empty_iff_forall_notMem]
          intro r hr
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
          omega
        rw [h0, he1, he2]
        simp
    | succ m ih =>
        intro hm
        have hmlt : m < 2 * k := by omega
        have ihm := ih (by omega)
        set s : Fin (2 * k) := ⟨m, hmlt⟩ with hsdef
        have hsv : s.val = m := rfl
        have h1 := card_lt_succ_split (fun r : Fin (2 * k) => ¬ isOpenStep i j r) s
        have h2 := card_lt_succ_split (fun r : Fin (2 * k) => isOpenStep i j r) s
        rw [hsv] at h1 h2
        by_cases hop : isOpenStep i j s
        · have hnm : stepEntry i j s ∉ oddSet i j m := hop
          have ho : oddSet i j (m + 1) = insert (stepEntry i j s) (oddSet i j m) := by
            have hh := oddSet_succ_of_open i j hop
            rwa [hsv] at hh
          rw [ho, Finset.card_insert_of_notMem hnm, h1, h2, if_neg (not_not_intro hop),
            if_pos hop]
          omega
        · have hmem : stepEntry i j s ∈ oddSet i j m := not_not.mp hop
          have ho : oddSet i j (m + 1) = (oddSet i j m).erase (stepEntry i j s) := by
            have hh := oddSet_succ_of_close i j hop
            rwa [hsv] at hh
          have hpos : 1 ≤ (oddSet i j m).card := Finset.card_pos.mpr ⟨_, hmem⟩
          rw [ho, Finset.card_erase_of_mem hmem, h1, h2, if_pos hop, if_neg hop]
          omega
  have h := key (2 * k) le_rfl
  have e1 : (univ.filter (fun r : Fin (2 * k) => r.val < 2 * k ∧ ¬ isOpenStep i j r))
      = closeSteps i j := by
    ext r
    simp only [closeSteps, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨fun hh => hh.2, fun hh => ⟨r.isLt, hh⟩⟩
  have e2 : (univ.filter (fun r : Fin (2 * k) => r.val < 2 * k ∧ isOpenStep i j r))
      = openSteps i j := by
    ext r
    simp only [openSteps, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨fun hh => hh.2, fun hh => ⟨r.isLt, hh⟩⟩
  rw [e1, e2] at h
  exact h

/-- **Lemma 6.** `#innov = V - 1`, without subtraction.

Proof: The innovative steps split by parity into `codeA` and `codeB`; `card_innov_even` counts
the first as the distinct rows and `card_innov_odd` the second as the distinct columns less
one. `Finset.filter_card_add_filter_neg_card_eq_card` on the parity of `s.val` splits the
filter. -/
theorem card_innovSteps {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (innovSteps i j).card + 1 = walkVerts i j := by
  classical
  have hA := card_innov_even i j
  have hB := card_innov_odd hk i j
  have e1 : (innovSteps i j).filter (fun s : Fin (2 * k) => Even s.val)
      = univ.filter (fun s : Fin (2 * k) => Even s.val ∧ isInnov i j s) := by
    ext x
    simp only [innovSteps, Finset.mem_filter, Finset.mem_univ, true_and]
    tauto
  have e2 : (innovSteps i j).filter (fun s : Fin (2 * k) => ¬ Even s.val)
      = univ.filter (fun s : Fin (2 * k) => Odd s.val ∧ isInnov i j s) := by
    ext x
    simp only [innovSteps, Finset.mem_filter, Finset.mem_univ, true_and,
      Nat.not_even_iff_odd]
    tauto
  have hsum := Finset.card_filter_add_card_filter_not
    (s := innovSteps i j) (p := fun s : Fin (2 * k) => Even s.val)
  rw [e1, e2] at hsum
  rw [walkVerts]
  omega

/-- **Proposition 7.** The opening steps are the innovative steps and the bad opening steps.

Proof: `innov_isOpenStep` gives the inclusion; `isInnov_or_isForced_or_isBad` with
`forced_not_isOpenStep` gives that an opening step is innovative or bad. Then
`Finset.filter_card_add_filter_neg_card_eq_card` inside `openSteps`. -/
theorem card_innovSteps_add_card_badOpen {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (innovSteps i j).card + (badOpen i j).card = (openSteps i j).card := by
  classical
  have e1 : (openSteps i j).filter (fun s : Fin (2 * k) => isInnov i j s) = innovSteps i j := by
    ext x
    simp only [openSteps, innovSteps, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨fun hh => hh.2, fun hh => ⟨innov_isOpenStep i j hh, hh⟩⟩
  have e2 : (openSteps i j).filter (fun s : Fin (2 * k) => ¬ isInnov i j s)
      = badOpen i j := by
    ext x
    simp only [openSteps, badOpen, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨hop, hni⟩
      have hb : isBad i j x := ⟨hni, fun hf => forced_not_isOpenStep i j hf hop⟩
      exact ⟨hb, hop⟩
    · rintro ⟨hb, hop⟩
      exact ⟨hop, hb.1⟩
  have hsum := Finset.card_filter_add_card_filter_not
    (s := openSteps i j) (p := fun s : Fin (2 * k) => isInnov i j s)
  rw [e1, e2] at hsum
  exact hsum

/-- The bad steps are the bad opening steps and the bad closing steps. -/
theorem card_codeP_eq {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) :
    (codeP i j).card = (badOpen i j).card + (badClose i j).card := by
  classical
  rw [badOpen, badClose, codeP]
  rw [← Finset.card_filter_add_card_filter_not (p := fun s => isOpenStep i j s)
    (s := univ.filter (fun s : Fin (2 * k) => isBad i j s))]
  congr 1
  · rw [Finset.filter_filter]
  · rw [Finset.filter_filter]

/-- **Proposition 8, first half.** `|Oend| ≤ 2 X`, in the form `|Oend| + 2 E ≤ 2 k`.

Proof: Every entry of `Oend` has odd multiplicity and `NoSingle` gives multiplicity at least
two, so at least three; then
`|Oend| + 2 |walkEdges| ≤ Σ_(e ∈ walkEdges) walkMult e = 2 k` by `sum_walkMult`. -/
theorem card_oddEnd_add_card_walkEdges_le {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d)
    (hns : NoSingle i j) :
    (oddEnd i j).card + 2 * (walkEdges i j).card ≤ 2 * k := by
  classical
  have hpm : ∀ e : Fin n × Fin d, prefixMult i j (2 * k) e = walkMult i j e := by
    intro e
    rw [walkMult_eq_card_stepEntry, prefixMult]
    congr 1
    ext r
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨fun h => h.2, fun h => ⟨r.isLt, h⟩⟩
  have hodd : ∀ e ∈ oddEnd i j, 3 ≤ walkMult i j e := by
    intro e he
    rw [oddEnd, mem_oddSet, hpm, Nat.odd_iff] at he
    have hne := hns e
    omega
  have hsub : oddEnd i j ⊆ walkEdges i j := by
    intro e he
    have h3 := hodd e he
    rw [walkEdges, Finset.mem_filter]
    exact ⟨Finset.mem_univ _, by omega⟩
  have hstep : ∀ e ∈ walkEdges i j,
      2 + (if e ∈ oddEnd i j then 1 else 0) ≤ walkMult i j e := by
    intro e he
    by_cases h : e ∈ oddEnd i j
    · rw [if_pos h]; exact hodd e h
    · rw [if_neg h]
      have := two_le_walkMult hns he
      omega
  have hsum := Finset.sum_le_sum hstep
  have hL : ∑ e ∈ walkEdges i j, (2 + (if e ∈ oddEnd i j then 1 else 0))
      = 2 * (walkEdges i j).card + (oddEnd i j).card := by
    rw [Finset.sum_add_distrib, Finset.sum_const, smul_eq_mul, Finset.sum_ite_mem,
      Finset.inter_eq_right.mpr hsub, Finset.sum_const, smul_eq_mul]
    ring
  have hR : ∑ e ∈ walkEdges i j, walkMult i j e ≤ 2 * k := by
    rw [← sum_walkMult i j]
    exact Finset.sum_le_sum_of_subset (Finset.subset_univ _)
  omega

/-- **Proposition 8, second half.** `V ≤ E + 1`. `Count.lean` already proves it. -/
theorem walkVerts_le_card_walkEdges_succ {k n d : ℕ} (_hk : 1 ≤ k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) : walkVerts i j ≤ (walkEdges i j).card + 1 :=
  walkVerts_le_edges_succ i j

/-! ### 7. The credit-debt invariant (Theorem 14) -/

/-- **Theorem 14, the invariant.** At every time `s` and every vertex `z`,
`debt_z(s) ≤ cred_z(s)` and `⌊l_z(s) / 2⌋ ≤ cred_z(s) - debt_z(s) + [z ∉ {u_0, u_s}]`.

Proof: Induction on `s`, with the case split of section 6 of
`scripts/stage3_count/x8_scratch/U3/u3_proof.md`: a vertex outside `{u_s, u_(s+1)}` (nothing
changes), the current vertex `u_s` (only the debt can grow, by `close_level_ge_two_isBad` and
`close_level_one_not_isBad`), the endpoint `u_(s+1)` (only the credit can grow, by
`open_arrival_isBad`). Every case closes with `omega` once the parity of the level (`level_parity`)
and the type of the step are in context. -/
theorem credit_invariant {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n) (j : Fin k → Fin d)
    (s : ℕ) (hs : s ≤ 2 * k) (z : Fin n ⊕ Fin d) :
    debtAt i j s z ≤ credAt i j s z ∧
      levelAt i j s z / 2 + debtAt i j s z
        ≤ credAt i j s z +
            (if z ≠ vertC hk i j 0 ∧ z ≠ vertC hk i j s then 1 else 0) := by
  classical
  revert hs
  induction s with
  | zero =>
      intro _
      have hc0 : credAt i j 0 z = 0 := by
        rw [credAt, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
        intro r hr
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
        exact absurd hr.1 (Nat.not_lt_zero _)
      have hd0 : debtAt i j 0 z = 0 := by
        rw [debtAt, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
        intro r hr
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
        exact absurd hr.1 (Nat.not_lt_zero _)
      have hl0 : levelAt i j 0 z = 0 := by
        rw [levelAt, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
        intro e he
        rw [Finset.mem_filter, mem_oddSet] at he
        have hz : prefixMult i j 0 e = 0 := by
          rw [prefixMult, Finset.card_eq_zero, Finset.eq_empty_iff_forall_notMem]
          intro r hr
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hr
          exact absurd hr.1 (Nat.not_lt_zero _)
        have he1 := he.1
        rw [hz, Nat.odd_iff] at he1
        omega
      rw [hc0, hd0, hl0]
      omega
  | succ m ih =>
      intro hs
      have hmlt : m < 2 * k := by omega
      obtain ⟨ihA, ihB⟩ := ih (by omega)
      set σ : Fin (2 * k) := ⟨m, hmlt⟩ with hσ
      have hσv : σ.val = m := rfl
      have hvm : vertC hk i j m = vert i j σ := by
        have hfin : (⟨m % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = σ :=
          Fin.ext (Nat.mod_eq_of_lt hmlt)
        rw [vertC, hfin]
      have hvm1 : vertC hk i j (m + 1) = vert i j (cycSucc σ) := by
        have hfin : (⟨(m + 1) % (2 * k), Nat.mod_lt _ hk⟩ : Fin (2 * k)) = cycSucc σ := rfl
        rw [vertC, hfin]
      have hne := vert_ne_vert_cycSucc i j σ
      have hcred : credAt i j (m + 1) z = credAt i j m z
          + (if isBad i j σ ∧ isOpenStep i j σ ∧ vert i j (cycSucc σ) = z then 1 else 0) := by
        have hh := card_lt_succ_split
          (fun r : Fin (2 * k) => isBad i j r ∧ isOpenStep i j r ∧ vert i j (cycSucc r) = z) σ
        rw [hσv] at hh
        exact hh
      have hdebt : debtAt i j (m + 1) z = debtAt i j m z
          + (if isBad i j σ ∧ ¬ isOpenStep i j σ ∧ vert i j σ = z then 1 else 0) := by
        have hh := card_lt_succ_split
          (fun r : Fin (2 * k) => isBad i j r ∧ ¬ isOpenStep i j r ∧ vert i j r = z) σ
        rw [hσv] at hh
        exact hh
      have hpar := level_parity hk i j m (by omega) z
      rw [hvm] at hpar ihB
      rw [hvm1]
      by_cases hx : z = vert i j σ
      · -- the current vertex: only the debt can grow
        have hy : z ≠ vert i j (cycSucc σ) := by rw [hx]; exact hne
        have hinc : Incident z (stepEntry i j σ) :=
          (incident_stepEntry_iff i j σ z).mpr (Or.inl hx)
        rw [if_neg (fun hcon => hcon.2 hx)] at ihB
        rw [if_pos hx] at hpar
        rw [if_neg (fun hcon => hy hcon.2.2.symm)] at hcred
        by_cases hop : isOpenStep i j σ
        · rw [if_neg (fun hcon => hcon.2.1 hop)] at hdebt
          have hlev : levelAt i j (m + 1) z = levelAt i j m z + 1 := by
            have hh := levelAt_succ_open_incident i j hop hinc
            rw [hσv] at hh
            exact hh
          by_cases hz0 : z = vertC hk i j 0
          · rw [if_pos hz0] at hpar
            rw [if_neg (fun hcon => hcon.1 hz0)]
            omega
          · rw [if_neg hz0] at hpar
            rw [if_pos ⟨hz0, hy⟩]
            omega
        · have hmem : stepEntry i j σ ∈ oddSet i j m := not_not.mp hop
          have hLpos : 1 ≤ levelAt i j m z := by
            rw [levelAt]
            exact Finset.card_pos.mpr ⟨stepEntry i j σ, Finset.mem_filter.mpr ⟨hmem, hinc⟩⟩
          have hlev : levelAt i j (m + 1) z + 1 = levelAt i j m z := by
            have hh := levelAt_succ_close_incident i j hop hinc
            rw [hσv] at hh
            exact hh
          by_cases hbad : isBad i j σ
          · rw [if_pos ⟨hbad, hop, hx.symm⟩] at hdebt
            have hL2 : 2 ≤ levelAt i j m z := by
              by_contra hcon
              have h1 : levelAt i j σ.val (vert i j σ) = 1 := by
                rw [hσv, ← hx]; omega
              exact close_level_one_not_isBad i j hop h1 hbad
            by_cases hz0 : z = vertC hk i j 0
            · rw [if_pos hz0] at hpar
              rw [if_neg (fun hcon => hcon.1 hz0)]
              omega
            · rw [if_neg hz0] at hpar
              rw [if_pos ⟨hz0, hy⟩]
              omega
          · rw [if_neg (fun hcon => hbad hcon.1)] at hdebt
            by_cases hz0 : z = vertC hk i j 0
            · rw [if_pos hz0] at hpar
              rw [if_neg (fun hcon => hcon.1 hz0)]
              omega
            · rw [if_neg hz0] at hpar
              rw [if_pos ⟨hz0, hy⟩]
              omega
      · rw [if_neg hx] at hpar
        rw [if_neg (fun hcon => hx hcon.2.2.symm)] at hdebt
        by_cases hy : z = vert i j (cycSucc σ)
        · -- the endpoint: only the credit can grow
          have hinc : Incident z (stepEntry i j σ) :=
            (incident_stepEntry_iff i j σ z).mpr (Or.inr hy)
          rw [if_neg (fun hcon => hcon.2 hy)]
          by_cases hop : isOpenStep i j σ
          · have hlev : levelAt i j (m + 1) z = levelAt i j m z + 1 := by
              have hh := levelAt_succ_open_incident i j hop hinc
              rw [hσv] at hh
              exact hh
            by_cases hz0 : z = vertC hk i j 0
            · have hbad : isBad i j σ :=
                open_arrival_isBad hk i j hop (Or.inl (hy.symm.trans hz0))
              rw [if_pos ⟨hbad, hop, hy.symm⟩] at hcred
              rw [if_pos hz0] at hpar
              rw [if_neg (fun hcon => hcon.1 hz0)] at ihB
              omega
            · rw [if_neg hz0] at hpar
              rw [if_pos ⟨hz0, hx⟩] at ihB
              by_cases hLpos : 1 ≤ levelAt i j m z
              · have hbad : isBad i j σ :=
                  open_arrival_isBad hk i j hop (Or.inr (by rw [hσv, ← hy]; exact hLpos))
                rw [if_pos ⟨hbad, hop, hy.symm⟩] at hcred
                omega
              · omega
          · have hmem : stepEntry i j σ ∈ oddSet i j m := not_not.mp hop
            have hLpos : 1 ≤ levelAt i j m z := by
              rw [levelAt]
              exact Finset.card_pos.mpr ⟨stepEntry i j σ, Finset.mem_filter.mpr ⟨hmem, hinc⟩⟩
            have hlev : levelAt i j (m + 1) z + 1 = levelAt i j m z := by
              have hh := levelAt_succ_close_incident i j hop hinc
              rw [hσv] at hh
              exact hh
            rw [if_neg (fun hcon => hop hcon.2.1)] at hcred
            by_cases hz0 : z = vertC hk i j 0
            · rw [if_pos hz0] at hpar
              rw [if_neg (fun hcon => hcon.1 hz0)] at ihB
              omega
            · rw [if_neg hz0] at hpar
              rw [if_pos ⟨hz0, hx⟩] at ihB
              omega
        · -- a vertex the step misses
          have hinc : ¬ Incident z (stepEntry i j σ) := by
            intro hc
            rcases (incident_stepEntry_iff i j σ z).mp hc with h | h
            · exact hx h
            · exact hy h
          have hlev : levelAt i j (m + 1) z = levelAt i j m z := by
            have hh := levelAt_succ_of_not_incident i j σ hinc
            rw [hσv] at hh
            exact hh
          rw [if_neg (fun hcon => hy hcon.2.2.symm)] at hcred
          by_cases hz0 : z = vertC hk i j 0
          · rw [if_neg (fun hcon => hcon.1 hz0)]
            rw [if_neg (fun hcon => hcon.1 hz0)] at ihB
            omega
          · rw [if_pos ⟨hz0, hy⟩]
            rw [if_pos ⟨hz0, hx⟩] at ihB
            omega

/-- **Corollary 15.** `#badclose ≤ #badopen`.

Proof: Sum part (i) of `credit_invariant` at `s = 2 k` over the vertex type `Fin n ⊕ Fin d`.
Both sides are a sum over `z` of a filtered card: `badClose` is partitioned by `vert i j r`
and `badOpen` by `vert i j (cycSucc r)`, so `Finset.card_eq_sum_card_fiberwise` on each side
turns the inequality into `Finset.sum_le_sum`. -/
theorem card_badClose_le_card_badOpen {k n d : ℕ} (hk : 0 < 2 * k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) : (badClose i j).card ≤ (badOpen i j).card := by
  classical
  have hd : ∀ z : Fin n ⊕ Fin d,
      ((badClose i j).filter (fun r => vert i j r = z)).card = debtAt i j (2 * k) z := by
    intro z
    rw [debtAt, badClose, Finset.filter_filter]
    congr 1
    ext r
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨⟨hb, hcl⟩, hv⟩
      exact ⟨r.isLt, hb, hcl, hv⟩
    · rintro ⟨-, hb, hcl, hv⟩
      exact ⟨⟨hb, hcl⟩, hv⟩
  have hc : ∀ z : Fin n ⊕ Fin d,
      ((badOpen i j).filter (fun r => vert i j (cycSucc r) = z)).card
        = credAt i j (2 * k) z := by
    intro z
    rw [credAt, badOpen, Finset.filter_filter]
    congr 1
    ext r
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨⟨hb, hop⟩, hv⟩
      exact ⟨r.isLt, hb, hop, hv⟩
    · rintro ⟨-, hb, hop, hv⟩
      exact ⟨⟨hb, hop⟩, hv⟩
  have h1 : (badClose i j).card = ∑ z : Fin n ⊕ Fin d, debtAt i j (2 * k) z := by
    rw [Finset.card_eq_sum_card_fiberwise
      (f := fun r : Fin (2 * k) => vert i j r) (t := (univ : Finset (Fin n ⊕ Fin d)))
      (fun x _ => Finset.mem_univ _)]
    exact Finset.sum_congr rfl (fun z _ => hd z)
  have h2 : (badOpen i j).card = ∑ z : Fin n ⊕ Fin d, credAt i j (2 * k) z := by
    rw [Finset.card_eq_sum_card_fiberwise
      (f := fun r : Fin (2 * k) => vert i j (cycSucc r)) (t := (univ : Finset (Fin n ⊕ Fin d)))
      (fun x _ => Finset.mem_univ _)]
    exact Finset.sum_congr rfl (fun z _ => hc z)
  rw [h1, h2]
  exact Finset.sum_le_sum (fun z _ => (credit_invariant hk i j (2 * k) le_rfl z).1)

/-! ### 8. The theorem (Theorem 16) -/

/-- **Theorem 16.** A walk with no entry of multiplicity one has at most `4 D` bad steps,
where `D = k + 1 - V` and `V = walkVerts i j`. -/
theorem card_codeP_le_of_noSingle {k n d : ℕ} (hk : 1 ≤ k) (i : Fin k → Fin n)
    (j : Fin k → Fin d) (hns : NoSingle i j) :
    (codeP i j).card ≤ 4 * (k + 1 - walkVerts i j) := by
  have hk2 : 0 < 2 * k := by omega
  have h1 := card_codeP_eq i j
  have h2 := card_badClose_le_card_badOpen hk2 i j
  have h3 := card_innovSteps_add_card_badOpen i j
  have h4 := card_closeSteps_add_card_oddEnd i j
  have h5 := card_openSteps_add_card_closeSteps i j
  have h6 := card_innovSteps hk i j
  have h7 := card_oddEnd_add_card_walkEdges_le i j hns
  have h8 := walkVerts_le_card_walkEdges_succ hk i j
  have h9 := walkVerts_le_of_noSingle hns
  omega

/-- **The target.** On the cell `(s, t)` the number of bad steps is at most
`4 (k + 1 - (s + t))`. -/
theorem card_codeP_le {n d k s t : ℕ} (hk : 1 ≤ k)
    {p : (Fin k → Fin n) × (Fin k → Fin d)} (hp : p ∈ cellSet n d k s t) :
    (codeP p.1 p.2).card ≤ 4 * (k + 1 - (s + t)) := by
  rw [mem_cellSet] at hp
  obtain ⟨hns, hrows, hcols⟩ := hp
  have h := card_codeP_le_of_noSingle hk p.1 p.2 hns
  rwa [walkVerts, hrows, hcols] at h

/-- The corollary the assembly reads: the cell count with the bad-step bound discharged. -/
theorem cellCard_le_code' {n d k s t : ℕ} (hk : 1 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hst : s + t ≤ k + 1) :
    cellCard n d k s t
      ≤ k.choose s * k.choose (t - 1) *
          ((4 * (k + 1 - (s + t)) + 1) * ((2 * k) ^ 8) ^ (k + 1 - (s + t))) *
          (n.descFactorial s * d.descFactorial t) :=
  cellCard_le_code_of_bad hk hs ht hst (fun _ hp => card_codeP_le hk hp)

end U3

end Edge
end StackedSVD
