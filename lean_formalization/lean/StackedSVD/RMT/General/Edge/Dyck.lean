/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Count
import Mathlib.Combinatorics.Enumerative.DyckWord

/-! # Stage 3, unit X, part 3: the tree walks of a cell and the Narayana lower bound

The lower half of the Furedi-Komlos count (`notes/x8_plan.md`, section 5b and items L1 to L4
of section 6). The target cell of the count has `s + t = k + 1`, so its walks are tree walks,
and the assembly needs a lower bound on how many there are.

`dyckWalk` builds the contour walk of a plane tree, presented as a Dyck word of semilength
`k` (`Mathlib.Combinatorics.Enumerative.DyckWord`). An up-letter descends to a new vertex and
a down-letter returns to the parent, so the walk is closed and every entry has multiplicity
2. `dyckSet k s` collects the words whose walk has `s` rows, hence `t = k + 1 - s` columns.

`dyck_label_card_le` labels the vertices by an embedding of the rows in `Fin n` and of the
columns in `Fin d`, and proves that the labeled walks inject into the cell, so the cell holds
at least `#(dyckSet k s) * n.descFactorial s * d.descFactorial t` walks.

`narayana_lower` gives `C(k, s) * C(k, s - 1) ≤ k ^ 3 * #(dyckSet k s)`. The proof is the
double rotation: a pair of `s`-subsets of `range k` names the even and the odd up-letters,
some rotation of the pair is a Dyck word, and a rotation orbit has at most `k ^ 2` members.
`tree_cell_lower` combines the two facts for the assembly. -/

/-! ## Part L: the lower chain (L1, L2, L3, L4) -/

/-! ## X8, item L1: the contour walk of a Dyck word

The plan is `notes/x8_plan.md`, section 5b and row L1 of section 6. The target of the item is
`tree_cell_lower` of `scripts/stage3_count/X8_skeleton.lean`: a cell with `s + t = k + 1`
holds at least `C(k,s) C(k,s-1) / k` labeled tree walks. This file builds the walk of one
plane tree, presented as a Dyck word `w` of semilength `k`.

## The construction

The contour of the plane tree visits, after `s` letters of `w`, the vertex at depth `h s`
(the height of the word after `s` letters). Depth parity equals letter parity, so the even
positions carry the columns and the odd positions the rows; the root is the column `j 0`.

A vertex is named by the position of the up-letter that creates it, so `dstack w.toList s` is
the list of the open positions of the current vertex and of its ancestors (an up-letter pushes
the position, a down-letter pops it) and `vid w.toList s` is `p + 1` at the vertex opened at
`p`, and `0` at the root. Then `vid` is even at a column and odd at a row, `vid L s ≤ s`, and
the label of the vertex is `vid L s / 2`, which is below `k`. So `dyckWalk` sends `w` to the
walk `(i, j)` with `i x` the row at position `2 x + 1` and `j x` the column at position `2 x`.

## The pairing

`NoSingle` holds because every tree edge is traversed exactly twice. The down-letter at `b`
and the up-letter at `a = (dstack L b).head` traverse the same edge (`matched_edge`), and the
two positions have opposite parity, so one of the two traversals is an even step and the other
an odd step. That is exactly the pair `noSingle_of_families` of `Count.lean` asks for.
`pop_partner` supplies the up-letter of a down-letter and `push_partner` the down-letter of an
up-letter; the second one uses that the stack is empty again after `2 k` letters. -/

open Finset DyckStep

namespace StackedSVD
namespace Edge

/-! ## The stack of open positions -/

/-- The letter of `L` at the position `s`, read as a down-letter past the end of `L`. -/
def letterAt (L : List DyckStep) (s : ℕ) : DyckStep := (L[s]?).getD D

/-- The list of open positions after `s` letters: the position that created the current
vertex, then the ones that created its ancestors. An up-letter pushes, a down-letter pops. -/
def dstack (L : List DyckStep) : ℕ → List ℕ
  | 0 => []
  | s + 1 => if letterAt L s = U then s :: dstack L s else (dstack L s).tail

/-- The name of a vertex: `p + 1` at the vertex opened by the up-letter at `p`, `0` at the
root. -/
def headId : List ℕ → ℕ
  | [] => 0
  | p :: _ => p + 1

/-- The name of the vertex the contour visits after `s` letters. -/
def vid (L : List DyckStep) (s : ℕ) : ℕ := headId (dstack L s)

@[simp] theorem dstack_zero (L : List DyckStep) : dstack L 0 = [] := rfl

/-- An up-letter pushes its own position. -/
theorem dstack_succ_U {L : List DyckStep} {s : ℕ} (h : letterAt L s = U) :
    dstack L (s + 1) = s :: dstack L s := by
  simp [dstack, h]

/-- A down-letter pops the top of the stack. -/
theorem dstack_succ_D {L : List DyckStep} {s : ℕ} (h : letterAt L s = D) :
    dstack L (s + 1) = (dstack L s).tail := by
  simp [dstack, h]

@[simp] theorem vid_zero (L : List DyckStep) : vid L 0 = 0 := rfl

/-- The name of the vertex of an empty stack is `0` (the root). -/
theorem vid_of_nil {L : List DyckStep} {s : ℕ} (h : dstack L s = []) : vid L s = 0 := by
  simp [vid, h, headId]

/-- The name of the vertex whose stack has the top `p` is `p + 1`. -/
theorem vid_of_cons {L : List DyckStep} {s p : ℕ} {rest : List ℕ}
    (h : dstack L s = p :: rest) : vid L s = p + 1 := by
  simp [vid, h, headId]

/-! ## The elementary stack invariants -/

/-- Every open position is below the current position. -/
theorem dstack_mem_lt (L : List DyckStep) : ∀ (s : ℕ), ∀ p ∈ dstack L s, p < s := by
  intro s
  induction s with
  | zero => simp
  | succ s ih =>
    intro p hp
    by_cases hU : letterAt L s = U
    · rw [dstack_succ_U hU] at hp
      rcases List.mem_cons.1 hp with rfl | hp
      · omega
      · exact Nat.lt_succ_of_lt (ih p hp)
    · rw [dstack_succ_D (by cases h : letterAt L s with
        | U => exact absurd h hU
        | D => rfl)] at hp
      exact Nat.lt_succ_of_lt (ih p (List.mem_of_mem_tail hp))

/-- Every open position carries an up-letter. -/
theorem dstack_mem_U (L : List DyckStep) : ∀ (s : ℕ), ∀ p ∈ dstack L s, letterAt L p = U := by
  intro s
  induction s with
  | zero => simp
  | succ s ih =>
    intro p hp
    by_cases hU : letterAt L s = U
    · rw [dstack_succ_U hU] at hp
      rcases List.mem_cons.1 hp with rfl | hp
      · exact hU
      · exact ih p hp
    · rw [dstack_succ_D (by cases h : letterAt L s with
        | U => exact absurd h hU
        | D => rfl)] at hp
      exact ih p (List.mem_of_mem_tail hp)

/-- The name is at most the position: `vid L s ≤ s`. This is what puts every label below `k`. -/
theorem vid_le (L : List DyckStep) (s : ℕ) : vid L s ≤ s := by
  cases h : dstack L s with
  | nil => simp [vid_of_nil h]
  | cons p rest =>
    rw [vid_of_cons h]
    have := dstack_mem_lt L s p (by rw [h]; exact List.mem_cons_self ..)
    omega

/-- The stack below the top is the stack at the position of the top. This is the fact that
makes the two traversals of a tree edge land on the same pair of vertices. -/
theorem dstack_tail_eq (L : List DyckStep) :
    ∀ (s : ℕ) (p : ℕ) (rest : List ℕ), dstack L s = p :: rest → rest = dstack L p := by
  intro s
  induction s using Nat.strong_induction_on with
  | _ s ih =>
    match s with
    | 0 => intro p rest h; simp at h
    | s + 1 =>
      intro p rest h
      by_cases hU : letterAt L s = U
      · rw [dstack_succ_U hU] at h
        obtain ⟨rfl, rfl⟩ := List.cons.inj h
        rfl
      · rw [dstack_succ_D (by cases hc : letterAt L s with
          | U => exact absurd hc hU
          | D => rfl)] at h
        cases hs : dstack L s with
        | nil => rw [hs] at h; simp at h
        | cons q qs =>
          rw [hs, List.tail_cons] at h
          have hq : qs = dstack L q := ih s (by omega) q qs hs
          have hqlt : q < s := dstack_mem_lt L s q (by rw [hs]; exact List.mem_cons_self ..)
          exact ih q (by omega) p rest (by rw [← hq]; exact h)

/-! ## The height of the stack -/

/-- Every letter is an up-letter or a down-letter. -/
private theorem count_U_add_count_D (l : List DyckStep) :
    l.count U + l.count D = l.length := by
  induction l with
  | nil => simp
  | cons a t ih => cases a <;> simp <;> omega

/-- The letter read from the list agrees with `letterAt` inside the list. -/
private theorem getElem?_eq_letterAt {L : List DyckStep} {b : ℕ} (hb : b < L.length) :
    L[b]? = some (letterAt L b) := by
  rw [letterAt, List.getElem?_eq_getElem hb]; rfl

/-- One more letter appends that letter to the prefix. -/
private theorem take_succ_eq (L : List DyckStep) {s : ℕ} {x : DyckStep} (h : L[s]? = some x) :
    L.take (s + 1) = L.take s ++ [x] := by
  rw [List.take_add_one, h]; rfl

private theorem count_U_singleton_U : ([U] : List DyckStep).count U = 1 := by decide
private theorem count_D_singleton_U : ([U] : List DyckStep).count D = 0 := by decide
private theorem count_U_singleton_D : ([D] : List DyckStep).count U = 0 := by decide
private theorem count_D_singleton_D : ([D] : List DyckStep).count D = 1 := by decide

/-- The stack height is the number of up-letters minus the number of down-letters read so
far. The balance condition of a Dyck word makes every pop legal. -/
theorem dstack_length (w : DyckWord) (s : ℕ) :
    (dstack w.toList s).length + (w.toList.take s).count D = (w.toList.take s).count U := by
  induction s with
  | zero => simp
  | succ s ih =>
    rcases hL : w.toList[s]? with _ | x
    · have hlen : w.toList.length ≤ s := by
        by_contra hc
        exact absurd hL (by simp [List.getElem?_eq_getElem (Nat.lt_of_not_le hc)])
      have hall : w.toList.take s = w.toList := List.take_of_length_le hlen
      have htake : w.toList.take (s + 1) = w.toList.take s := by
        rw [hall, List.take_of_length_le (by omega)]
      have hz : (dstack w.toList s).length = 0 := by
        rw [hall] at ih
        have := w.count_U_eq_count_D
        omega
      have hnil : dstack w.toList s = [] := List.length_eq_zero_iff.1 hz
      have hD : letterAt w.toList s = D := by simp [letterAt, hL]
      rw [dstack_succ_D hD, hnil, htake]
      simpa [hz] using ih
    · have hx : letterAt w.toList s = x := by simp [letterAt, hL]
      have htake : w.toList.take (s + 1) = w.toList.take s ++ [x] := take_succ_eq _ hL
      cases x with
      | U =>
        rw [dstack_succ_U hx, htake, List.count_append, List.count_append,
          count_U_singleton_U, count_D_singleton_U, List.length_cons]
        omega
      | D =>
        have hbal := w.count_D_le_count_U (s + 1)
        rw [htake, List.count_append, List.count_append, count_U_singleton_D,
          count_D_singleton_D] at hbal
        have hpos : 1 ≤ (dstack w.toList s).length := by omega
        rw [dstack_succ_D hx, htake, List.count_append, List.count_append,
          count_U_singleton_D, count_D_singleton_D, List.length_tail]
        omega

/-- The stack is empty again after the whole word. -/
theorem dstack_end (w : DyckWord) : dstack w.toList (2 * w.semilength) = [] := by
  have h1 := dstack_length w (2 * w.semilength)
  have hall : w.toList.take (2 * w.semilength) = w.toList :=
    List.take_of_length_le (by rw [w.two_mul_semilength_eq_length])
  rw [hall, w.count_U_eq_count_D] at h1
  exact List.length_eq_zero_iff.1 (by omega)

/-- The root is the vertex at the position `2 k`, as it is at the position `0`. -/
theorem vid_end (w : DyckWord) : vid w.toList (2 * w.semilength) = 0 :=
  vid_of_nil (dstack_end w)

/-- The depth parity equals the position parity: the even positions carry the columns and
the odd positions the rows. -/
theorem dstack_length_parity (w : DyckWord) {s : ℕ} (hs : s ≤ 2 * w.semilength) :
    (dstack w.toList s).length % 2 = s % 2 := by
  have h1 := dstack_length w s
  have h2 := count_U_add_count_D (w.toList.take s)
  have h3 : (w.toList.take s).length = s := by
    rw [List.length_take]
    have := w.two_mul_semilength_eq_length
    omega
  omega

/-- The stack is not empty at a down-letter: the pop is legal. -/
theorem dstack_ne_nil_of_D (w : DyckWord) {b : ℕ} (hb : b < 2 * w.semilength)
    (hD : letterAt w.toList b = D) : dstack w.toList b ≠ [] := by
  have hblen : b < w.toList.length := by
    have := w.two_mul_semilength_eq_length
    omega
  have hL : w.toList[b]? = some D := by rw [getElem?_eq_letterAt hblen, hD]
  have htake : w.toList.take (b + 1) = w.toList.take b ++ [D] := take_succ_eq _ hL
  have hbal := w.count_D_le_count_U (b + 1)
  rw [htake, List.count_append, List.count_append, count_U_singleton_D,
    count_D_singleton_D] at hbal
  have h1 := dstack_length w b
  intro hnil
  rw [hnil] at h1
  simp at h1
  omega

/-! ## The two traversals of one tree edge -/

/-- The down-letter at `b` and the up-letter at `a`, the top of the stack at `b`, traverse the
same tree edge in opposite directions: `vid L b = vid L (a + 1)` is the child and
`vid L (b + 1) = vid L a` is the parent. The stack is one deeper at `b` than at `a`. -/
theorem matched_edge {L : List DyckStep} {b a : ℕ} {rest : List ℕ}
    (hD : letterAt L b = D) (hb : dstack L b = a :: rest) :
    a < b ∧ letterAt L a = U ∧ vid L b = vid L (a + 1) ∧ vid L (b + 1) = vid L a ∧
      (dstack L b).length = (dstack L a).length + 1 := by
  have hmem : a ∈ dstack L b := by rw [hb]; exact List.mem_cons_self ..
  have halt : a < b := dstack_mem_lt L b a hmem
  have haU : letterAt L a = U := dstack_mem_U L b a hmem
  have hrest : rest = dstack L a := dstack_tail_eq L b a rest hb
  have hba : dstack L b = a :: dstack L a := by rw [hb, hrest]
  refine ⟨halt, haU, ?_, ?_, ?_⟩
  · rw [vid_of_cons hba, vid_of_cons (dstack_succ_U haU)]
  · change headId (dstack L (b + 1)) = headId (dstack L a)
    rw [dstack_succ_D hD, hba, List.tail_cons]
  · rw [hba, List.length_cons]

/-- Every up-letter inside the word is popped later by a down-letter, because the stack is
empty again after the whole word. -/
theorem push_partner (w : DyckWord) {a : ℕ} (ha : a < 2 * w.semilength)
    (haU : letterAt w.toList a = U) :
    ∃ b, a < b ∧ b < 2 * w.semilength ∧ letterAt w.toList b = D ∧
      dstack w.toList b = a :: dstack w.toList a := by
  set L := w.toList with hL
  have hex : ∃ m, a ∉ dstack L (a + 1 + m) := by
    refine ⟨2 * w.semilength - (a + 1), ?_⟩
    have : a + 1 + (2 * w.semilength - (a + 1)) = 2 * w.semilength := by omega
    rw [this, dstack_end]
    simp
  classical
  set m := Nat.find hex with hm
  have hmspec : a ∉ dstack L (a + 1 + m) := Nat.find_spec hex
  have hmpos : 0 < m := by
    rcases Nat.eq_zero_or_pos m with h0 | h; swap
    · exact h
    · exfalso
      rw [h0] at hmspec
      exact hmspec (by rw [dstack_succ_U haU]; exact List.mem_cons_self ..)
  have hmle : m ≤ 2 * w.semilength - (a + 1) := Nat.find_min' hex (by
    have : a + 1 + (2 * w.semilength - (a + 1)) = 2 * w.semilength := by omega
    rw [this, dstack_end]; simp)
  have hprev : a ∈ dstack L (a + 1 + (m - 1)) := by
    by_contra hc
    have := Nat.find_min' hex hc
    omega
  refine ⟨a + m, by omega, by omega, ?_, ?_⟩
  · by_contra hcU
    have hU : letterAt L (a + m) = U := by
      cases hc : letterAt L (a + m) with
      | U => rfl
      | D => exact absurd hc hcU
    have hstep : a + 1 + m = (a + m) + 1 := by omega
    rw [hstep, dstack_succ_U hU] at hmspec
    exact hmspec (List.mem_cons_of_mem _ (by
      have : a + 1 + (m - 1) = a + m := by omega
      rwa [this] at hprev))
  · have hD : letterAt L (a + m) = D := by
      by_contra hcU
      have hU : letterAt L (a + m) = U := by
        cases hc : letterAt L (a + m) with
        | U => rfl
        | D => exact absurd hc hcU
      have hstep : a + 1 + m = (a + m) + 1 := by omega
      rw [hstep, dstack_succ_U hU] at hmspec
      exact hmspec (List.mem_cons_of_mem _ (by
        have : a + 1 + (m - 1) = a + m := by omega
        rwa [this] at hprev))
    have hin : a ∈ dstack L (a + m) := by
      have : a + 1 + (m - 1) = a + m := by omega
      rwa [this] at hprev
    have hout : a ∉ (dstack L (a + m)).tail := by
      have hstep : a + 1 + m = (a + m) + 1 := by omega
      rw [hstep, dstack_succ_D hD] at hmspec
      exact hmspec
    cases hs : dstack L (a + m) with
    | nil => rw [hs] at hin; simp at hin
    | cons q rest =>
      have hrest : rest = dstack L q := dstack_tail_eq L (a + m) q rest hs
      rw [hs] at hin hout
      rw [List.tail_cons] at hout
      have hqa : q = a := by
        rcases List.mem_cons.1 hin with h | h
        · exact h.symm
        · exact absurd h hout
      rw [hrest, hqa]

/-- Every step of the contour has a partner step of the opposite parity that traverses the
same tree edge in the opposite direction. This is the whole content of `NoSingle`. -/
theorem exists_partner (w : DyckWord) {p : ℕ} (hp : p < 2 * w.semilength) :
    ∃ q, q < 2 * w.semilength ∧ q % 2 ≠ p % 2 ∧
      vid w.toList p = vid w.toList (q + 1) ∧ vid w.toList (p + 1) = vid w.toList q := by
  cases hc : letterAt w.toList p with
  | D =>
    obtain ⟨a, rest, hb⟩ : ∃ a rest, dstack w.toList p = a :: rest := by
      cases hs : dstack w.toList p with
      | nil => exact absurd hs (dstack_ne_nil_of_D w hp hc)
      | cons a rest => exact ⟨a, rest, rfl⟩
    obtain ⟨halt, haU, h1, h2, h3⟩ := matched_edge hc hb
    have hpar1 := dstack_length_parity w (s := p) (by omega)
    have hpar2 := dstack_length_parity w (s := a) (by omega)
    exact ⟨a, by omega, by omega, h1, h2⟩
  | U =>
    obtain ⟨b, hab, hb2k, hbD, hbs⟩ := push_partner w hp hc
    obtain ⟨halt, haU, h1, h2, h3⟩ := matched_edge hbD hbs
    have hpar1 := dstack_length_parity w (s := p) (by omega)
    have hpar2 := dstack_length_parity w (s := b) (by omega)
    exact ⟨b, hb2k, by omega, h2.symm, h1.symm⟩

/-! ## The walk -/

set_option linter.unusedVariables false in
/-- The contour walk of a Dyck word `w` of semilength `k`, as a closed walk of length `2 k`
on `Fin n × Fin d`. The row `i x` is the vertex at the position `2 x + 1` and the column
`j x` the vertex at the position `2 x`; both labels are `vid / 2`, which is below `k`. -/
def dyckWalk {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n) (hd : k ≤ d) :
    (Fin k → Fin n) × (Fin k → Fin d) :=
  (fun x => ⟨vid w.toList (2 * x.1 + 1) / 2, by
      have h := vid_le w.toList (2 * x.1 + 1); have := x.2; omega⟩,
   fun x => ⟨vid w.toList (2 * x.1) / 2, by
      have h := vid_le w.toList (2 * x.1); have := x.2; omega⟩)

/-- The row label of the contour walk. -/
@[simp] theorem dyckWalk_fst_val {k n d : ℕ} {w : DyckWord} {hw : w.semilength = k}
    {hn : k ≤ n} {hd : k ≤ d} (x : Fin k) :
    (((dyckWalk w hw hn hd).1 x : Fin n) : ℕ) = vid w.toList (2 * x.1 + 1) / 2 := rfl

/-- The column label of the contour walk. -/
@[simp] theorem dyckWalk_snd_val {k n d : ℕ} {w : DyckWord} {hw : w.semilength = k}
    {hn : k ≤ n} {hd : k ≤ d} (x : Fin k) :
    (((dyckWalk w hw hn hd).2 x : Fin d) : ℕ) = vid w.toList (2 * x.1) / 2 := rfl

/-- The cyclic wrap of the column index: at `y + 1 = k` the walk returns to the root. -/
theorem vid_cyc {k : ℕ} {w : DyckWord} (hw : w.semilength = k) {y : ℕ} (hy : y < k) :
    vid w.toList (2 * ((y + 1) % k)) = vid w.toList (2 * y + 2) := by
  rcases Nat.lt_or_ge (y + 1) k with h | h
  · rw [Nat.mod_eq_of_lt h]; ring_nf
  · have hyk : y + 1 = k := by omega
    rw [hyk, Nat.mod_self]
    have : 2 * y + 2 = 2 * w.semilength := by omega
    rw [this, vid_end]
    simp

/-- L1, the main deliverable: the contour walk of a Dyck word has no entry of multiplicity 1.
Every tree edge is traversed exactly twice, once by an even step and once by an odd step. -/
theorem noSingle_dyckWalk {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) : NoSingle (dyckWalk w hw hn hd).1 (dyckWalk w hw hn hd).2 := by
  subst hw
  refine noSingle_of_families ?_ ?_
  · intro x
    obtain ⟨q, hq2k, hqpar, h1, h2⟩ := exists_partner w (p := 2 * x.1) (by omega)
    have hqodd : q % 2 = 1 := by omega
    have hylt : q / 2 < w.semilength := by omega
    refine ⟨⟨q / 2, hylt⟩, ?_, ?_⟩
    · apply Fin.ext
      rw [dyckWalk_fst_val, dyckWalk_fst_val]
      have : 2 * (q / 2) + 1 = q := by omega
      rw [this, ← h2]
    · apply Fin.ext
      rw [dyckWalk_snd_val, dyckWalk_snd_val]
      have hc : ((cycSucc (⟨q / 2, hylt⟩ : Fin w.semilength)) : Fin w.semilength).1
          = (q / 2 + 1) % w.semilength := rfl
      rw [hc, vid_cyc rfl hylt]
      have : 2 * (q / 2) + 2 = q + 1 := by omega
      rw [this, ← h1]
  · intro x
    obtain ⟨q, hq2k, hqpar, h1, h2⟩ := exists_partner w (p := 2 * x.1 + 1) (by omega)
    have hqeven : q % 2 = 0 := by omega
    have hylt : q / 2 < w.semilength := by omega
    refine ⟨⟨q / 2, hylt⟩, ?_, ?_⟩
    · apply Fin.ext
      rw [dyckWalk_fst_val, dyckWalk_fst_val]
      have : 2 * (q / 2) + 1 = q + 1 := by omega
      rw [this, ← h1]
    · apply Fin.ext
      rw [dyckWalk_snd_val, dyckWalk_snd_val]
      have hc : ((cycSucc x) : Fin w.semilength).1 = (x.1 + 1) % w.semilength := rfl
      rw [hc, vid_cyc rfl x.2]
      have h2' : 2 * (q / 2) = q := by omega
      rw [h2']
      have : 2 * x.1 + 2 = 2 * x.1 + 1 + 1 := by omega
      rw [this, h2]

/-! ## The row and column counts -/

/-- The number of rows of the contour walk: an up-letter at an even position opens a vertex
of odd depth, that is a row. -/
def dyckRows (w : DyckWord) : ℕ :=
  ((range (2 * w.semilength)).filter (fun s => s % 2 = 0 ∧ letterAt w.toList s = U)).card

/-- The number of columns of the contour walk: the root, plus one column for every up-letter
at an odd position. -/
def dyckCols (w : DyckWord) : ℕ :=
  ((range (2 * w.semilength)).filter (fun s => s % 2 = 1 ∧ letterAt w.toList s = U)).card + 1

/-- The row count is the number of up-letters at even positions, by definition. -/
theorem dyckRows_eq_card_even_ups (w : DyckWord) :
    dyckRows w
      = ((range (2 * w.semilength)).filter (fun s => s % 2 = 0 ∧ letterAt w.toList s = U)).card :=
  rfl

/-- The positions of a letter, counted over the index range. -/
private theorem card_filter_letterAt (L : List DyckStep) (x : DyckStep) :
    ((range L.length).filter (fun s => letterAt L s = x)).card = L.count x := by
  classical
  induction L with
  | nil => simp
  | cons a t ih =>
    have hshift : ∀ s, letterAt (a :: t) (s + 1) = letterAt t s := fun s => by simp [letterAt]
    have h0 : letterAt (a :: t) 0 = a := by simp [letterAt]
    rw [Finset.card_filter] at ih ⊢
    rw [List.length_cons, Finset.sum_range_succ']
    simp only [hshift, h0, ih, List.count_cons]
    by_cases h : a = x <;> simp [h]

/-- The word has `k` up-letters inside its `2 k` positions. -/
private theorem card_filter_U (w : DyckWord) :
    ((range (2 * w.semilength)).filter (fun s => letterAt w.toList s = U)).card
      = w.semilength := by
  have hlen : 2 * w.semilength = w.toList.length := w.two_mul_semilength_eq_length
  rw [hlen, card_filter_letterAt w.toList U]
  rfl

/-- The tree of a Dyck word of semilength `k` has `k + 1` vertices: `dyckRows w` rows, at the
even up-letters, and `dyckCols w` columns, at the root and the odd up-letters. -/
theorem dyckRows_add_dyckCols (w : DyckWord) {k : ℕ} (hw : w.semilength = k) :
    dyckRows w + dyckCols w = k + 1 := by
  classical
  subst hw
  have hU : ∑ s ∈ range (2 * w.semilength), (if letterAt w.toList s = U then 1 else 0)
      = w.semilength := by
    rw [← Finset.card_filter]; exact card_filter_U w
  have hpoint : ∀ s ∈ range (2 * w.semilength),
      ((if (s % 2 = 0 ∧ letterAt w.toList s = U) then 1 else 0)
        + (if (s % 2 = 1 ∧ letterAt w.toList s = U) then 1 else 0))
      = (if letterAt w.toList s = U then 1 else 0) := by
    intro s _
    have hp : s % 2 = 0 ∨ s % 2 = 1 := by omega
    by_cases h : letterAt w.toList s = U
    · rcases hp with h2 | h2 <;> simp [h, h2]
    · simp [h]
  have hsum := Finset.sum_congr rfl hpoint
  rw [Finset.sum_add_distrib, hU] at hsum
  rw [dyckRows, dyckCols, Finset.card_filter, Finset.card_filter]
  omega

/-- The top of a nonempty stack has the parity opposite to the position, so the vertex name
`vid` has the parity of the position: even at a column, odd at a row. -/
theorem top_parity (w : DyckWord) {s p : ℕ} {rest : List ℕ} (hs : s ≤ 2 * w.semilength)
    (h : dstack w.toList s = p :: rest) : (p + 1) % 2 = s % 2 := by
  have hrest : rest = dstack w.toList p := dstack_tail_eq _ s p rest h
  have hplt : p < s := dstack_mem_lt _ s p (by rw [h]; exact List.mem_cons_self ..)
  have h1 := dstack_length_parity w (s := s) hs
  have h2 := dstack_length_parity w (s := p) (by omega)
  rw [h, hrest, List.length_cons] at h1
  omega

/-- At an odd position the stack is not empty: the walk is at a row, never at the root. -/
theorem dstack_ne_nil_of_odd (w : DyckWord) {s : ℕ} (hs : s ≤ 2 * w.semilength)
    (hodd : s % 2 = 1) : dstack w.toList s ≠ [] := by
  intro h
  have h1 := dstack_length_parity w hs
  rw [h] at h1
  simp only [List.length_nil] at h1
  omega

/-- The row labels of the walk, read as natural numbers. -/
theorem card_image_fst_val {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) :
    (univ.image (dyckWalk w hw hn hd).1).card
      = (univ.image (fun x : Fin k => vid w.toList (2 * x.1 + 1) / 2)).card := by
  classical
  rw [← Finset.card_image_of_injective (univ.image (dyckWalk w hw hn hd).1) Fin.val_injective,
    Finset.image_image]
  rfl

/-- The column labels of the walk, read as natural numbers. -/
theorem card_image_snd_val {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) :
    (univ.image (dyckWalk w hw hn hd).2).card
      = (univ.image (fun x : Fin k => vid w.toList (2 * x.1) / 2)).card := by
  classical
  rw [← Finset.card_image_of_injective (univ.image (dyckWalk w hw hn hd).2) Fin.val_injective,
    Finset.image_image]
  rfl

/-- The rows of the contour walk are exactly the vertices opened by an up-letter at an even
position, so the walk has `dyckRows w` distinct rows. -/
theorem card_image_dyckWalk_fst {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k)
    (hn : k ≤ n) (hd : k ≤ d) :
    (univ.image (dyckWalk w hw hn hd).1).card = dyckRows w := by
  classical
  rw [card_image_fst_val w hw hn hd]
  subst hw
  have himg : (univ.image (fun x : Fin w.semilength => vid w.toList (2 * x.1 + 1) / 2))
      = ((range (2 * w.semilength)).filter
          (fun p => p % 2 = 0 ∧ letterAt w.toList p = U)).image (fun p => p / 2) := by
    ext u
    simp only [Finset.mem_image, Finset.mem_univ, true_and, Finset.mem_filter, Finset.mem_range]
    constructor
    · rintro ⟨x, rfl⟩
      have hx := x.2
      obtain ⟨p, rest, hp⟩ : ∃ p rest, dstack w.toList (2 * x.1 + 1) = p :: rest := by
        cases hc : dstack w.toList (2 * x.1 + 1) with
        | nil => exact absurd hc (dstack_ne_nil_of_odd w (by omega) (by omega))
        | cons p rest => exact ⟨p, rest, rfl⟩
      have hmem : p ∈ dstack w.toList (2 * x.1 + 1) := by rw [hp]; exact List.mem_cons_self ..
      have hpar := top_parity w (s := 2 * x.1 + 1) (by omega) hp
      have hplt : p < 2 * x.1 + 1 := dstack_mem_lt _ _ p hmem
      have hpU : letterAt w.toList p = U := dstack_mem_U _ _ p hmem
      exact ⟨p, ⟨by omega, by omega, hpU⟩, by rw [vid_of_cons hp]; omega⟩
    · rintro ⟨p, ⟨hplt, hpev, hpU⟩, rfl⟩
      refine ⟨⟨p / 2, by omega⟩, ?_⟩
      have h1 : 2 * (p / 2) + 1 = p + 1 := by omega
      simp only []
      rw [h1, vid_of_cons (dstack_succ_U hpU)]
      omega
  have hinj : Set.InjOn (fun p : ℕ => p / 2)
      (((range (2 * w.semilength)).filter
        (fun p => p % 2 = 0 ∧ letterAt w.toList p = U) : Finset ℕ) : Set ℕ) := by
    intro a ha b hb hab
    simp only [Finset.coe_filter, Set.mem_ofPred_eq] at ha hb
    simp only at hab
    omega
  rw [himg, Finset.card_image_of_injOn hinj, dyckRows]

/-- The columns of the contour walk are the root and the vertices opened by an up-letter at an
odd position, so the walk has `dyckCols w` distinct columns. The hypothesis `1 ≤ k` is needed:
at `k = 0` the walk is empty and has no column, while `dyckCols w = 1`. -/
theorem card_image_dyckWalk_snd {k n d : ℕ} (w : DyckWord) (hk : 1 ≤ k) (hw : w.semilength = k)
    (hn : k ≤ n) (hd : k ≤ d) :
    (univ.image (dyckWalk w hw hn hd).2).card = dyckCols w := by
  classical
  rw [card_image_snd_val w hw hn hd]
  subst hw
  have himg : (univ.image (fun x : Fin w.semilength => vid w.toList (2 * x.1) / 2))
      = insert 0 (((range (2 * w.semilength)).filter
          (fun p => p % 2 = 1 ∧ letterAt w.toList p = U)).image (fun p => (p + 1) / 2)) := by
    ext u
    simp only [Finset.mem_image, Finset.mem_univ, true_and, Finset.mem_filter, Finset.mem_range,
      Finset.mem_insert]
    constructor
    · rintro ⟨x, rfl⟩
      have hx := x.2
      cases hp : dstack w.toList (2 * x.1) with
      | nil => exact Or.inl (by rw [vid_of_nil hp])
      | cons p rest =>
        have hmem : p ∈ dstack w.toList (2 * x.1) := by rw [hp]; exact List.mem_cons_self ..
        have hpar := top_parity w (s := 2 * x.1) (by omega) hp
        have hplt : p < 2 * x.1 := dstack_mem_lt _ _ p hmem
        have hpU : letterAt w.toList p = U := dstack_mem_U _ _ p hmem
        exact Or.inr ⟨p, ⟨by omega, by omega, hpU⟩, by rw [vid_of_cons hp]⟩
    · rintro (rfl | ⟨p, ⟨hplt, hpodd, hpU⟩, rfl⟩)
      · exact ⟨⟨0, hk⟩, by simp⟩
      · obtain ⟨b, hab, hb2k, -, -⟩ := push_partner w hplt hpU
        refine ⟨⟨(p + 1) / 2, by omega⟩, ?_⟩
        have h1 : 2 * ((p + 1) / 2) = p + 1 := by omega
        simp only []
        rw [h1, vid_of_cons (dstack_succ_U hpU)]
  have hzero : (0 : ℕ) ∉ (((range (2 * w.semilength)).filter
      (fun p => p % 2 = 1 ∧ letterAt w.toList p = U)).image (fun p => (p + 1) / 2)) := by
    simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_range, not_exists]
    rintro p ⟨⟨-, hpodd, -⟩, hp0⟩
    omega
  have hinj : Set.InjOn (fun p : ℕ => (p + 1) / 2)
      (((range (2 * w.semilength)).filter
        (fun p => p % 2 = 1 ∧ letterAt w.toList p = U) : Finset ℕ) : Set ℕ) := by
    intro a ha b hb hab
    simp only [Finset.coe_filter, Set.mem_ofPred_eq] at ha hb
    simp only at hab
    omega
  rw [himg, Finset.card_insert_of_notMem hzero, Finset.card_image_of_injOn hinj, dyckCols]

/-- L1: the contour walk of a Dyck word of semilength `k ≥ 1` visits `k + 1` vertices, so it
is a tree walk, and it lands in the cell `(dyckRows w, dyckCols w)`. -/
theorem walkVerts_dyckWalk' {k n d : ℕ} (w : DyckWord) (hk : 1 ≤ k) (hw : w.semilength = k)
    (hn : k ≤ n) (hd : k ≤ d) :
    walkVerts (dyckWalk w hw hn hd).1 (dyckWalk w hw hn hd).2 = k + 1 := by
  rw [walkVerts, card_image_dyckWalk_fst w hw hn hd, card_image_dyckWalk_snd w hk hw hn hd,
    dyckRows_add_dyckCols w hw]

/-! ## What L2 reads -/

/-- The first letter of a nonempty Dyck word is an up-letter, so the walk opens a row at the
position `0`. -/
theorem letterAt_zero (w : DyckWord) (hk : 1 ≤ w.semilength) : letterAt w.toList 0 = U := by
  have hbal := w.count_D_le_count_U 1
  cases hc : letterAt w.toList 0 with
  | U => rfl
  | D =>
    exfalso
    have hlen : 0 < w.toList.length := by
      have := w.two_mul_semilength_eq_length; omega
    have hL : w.toList[0]? = some D := by rw [getElem?_eq_letterAt hlen, hc]
    have htake : w.toList.take 1 = w.toList.take 0 ++ [D] := take_succ_eq _ hL
    rw [htake, List.count_append, List.count_append, count_U_singleton_D,
      count_D_singleton_D] at hbal
    simp at hbal

/-- A nonempty word has at least one row. -/
theorem one_le_dyckRows (w : DyckWord) (hk : 1 ≤ w.semilength) : 1 ≤ dyckRows w := by
  rw [dyckRows]
  exact Finset.card_pos.mpr ⟨0, Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (by omega),
    ⟨rfl, letterAt_zero w hk⟩⟩⟩

/-- Every word has at least one column, the root. -/
theorem one_le_dyckCols (w : DyckWord) : 1 ≤ dyckCols w := by
  rw [dyckCols]; omega

/-- The row count of a word of semilength `k` is between `1` and `k`. -/
theorem dyckRows_le (w : DyckWord) {k : ℕ} (hw : w.semilength = k) : dyckRows w ≤ k := by
  have h := dyckRows_add_dyckCols w hw
  have := one_le_dyckCols w
  omega

/-- The column count of a word of semilength `k` is between `1` and `k`. -/
theorem dyckCols_le (w : DyckWord) {k : ℕ} (hk : 1 ≤ k) (hw : w.semilength = k) :
    dyckCols w ≤ k := by
  have h := dyckRows_add_dyckCols w hw
  have := one_le_dyckRows w (by omega)
  omega

/-- L1, assembled for L2: the contour walk of a Dyck word of semilength `k ≥ 1` is a walk of
the cell `(dyckRows w, dyckCols w)`, a cell of the tree class since the two counts add up to
`k + 1`. -/
theorem dyckWalk_mem_cellSet {k n d : ℕ} (w : DyckWord) (hk : 1 ≤ k) (hw : w.semilength = k)
    (hn : k ≤ n) (hd : k ≤ d) :
    dyckWalk w hw hn hd ∈ cellSet n d k (dyckRows w) (dyckCols w) :=
  mem_cellSet.mpr ⟨noSingle_dyckWalk w hw hn hd, card_image_dyckWalk_fst w hw hn hd,
    card_image_dyckWalk_snd w hk hw hn hd⟩

/-! ## What L2 reads: the word is recoverable from the walk -/

/-- The name of a vertex has the parity of its position, so `vid` is recoverable from the
label `vid / 2` and the parity of the position. -/
theorem vid_parity (w : DyckWord) {s : ℕ} (hs : s ≤ 2 * w.semilength) :
    vid w.toList s % 2 = s % 2 := by
  cases h : dstack w.toList s with
  | nil =>
    have h1 := dstack_length_parity w hs
    rw [h] at h1
    simp only [List.length_nil] at h1
    rw [vid_of_nil h]
    omega
  | cons p rest =>
    rw [vid_of_cons h]
    exact top_parity w hs h

/-- The vertex name from its label: `vid L s = 2 * (label) + (s % 2)`. -/
theorem vid_eq_two_mul_label (w : DyckWord) {s : ℕ} (hs : s ≤ 2 * w.semilength) :
    vid w.toList s = 2 * (vid w.toList s / 2) + s % 2 := by
  have := vid_parity w hs
  omega

/-- A position carries an up-letter exactly when it opens a new vertex. With
`vid_eq_two_mul_label` this recovers the word from the labeled walk, which is what the
injection of L2 needs. -/
theorem letterAt_eq_U_iff (L : List DyckStep) (s : ℕ) :
    letterAt L s = U ↔ vid L (s + 1) = s + 1 := by
  constructor
  · intro h
    rw [vid_of_cons (dstack_succ_U h)]
  · intro h
    by_contra hc
    have hD : letterAt L s = D := by
      cases hx : letterAt L s with
      | U => exact absurd hx hc
      | D => rfl
    have hlt : vid L (s + 1) ≤ s := by
      cases hx : dstack L (s + 1) with
      | nil => rw [vid_of_nil hx]; omega
      | cons p rest =>
        rw [vid_of_cons hx]
        have hmem : p ∈ (dstack L s).tail := by
          rw [← dstack_succ_D hD, hx]; exact List.mem_cons_self ..
        have := dstack_mem_lt L s p (List.mem_of_mem_tail hmem)
        omega
    omega


/-! ## X8, item L2: the labeled tree walks inject into the cell

The plan is `notes/x8_plan.md`, section 5b and row L2 of section 6. L1 sends a Dyck word `w`
of semilength `k` to a walk `dyckWalk w` of the cell `(dyckRows w, dyckCols w)`. L2 labels
that walk: an injection of the rows the walk uses into `Fin n` and an injection of the
columns into `Fin d` give a new walk of the same cell, and the map
`(w, row labeling, column labeling) -> walk` is injective.

## The two halves

1. **The image is in the cell.** A relabeling by an injection changes no repetition pattern,
   so `noSingle_of_pattern` carries `NoSingle` across, and the row count and the column count
   are unchanged.
2. **The map is injective.** The word is recoverable from the pattern of the walk. A position
   `p` of `w` carries an up-letter exactly when the walk is at a first visit just after `p`
   (`letterAt_even_iff`, `letterAt_odd_iff`), and "first visit" is a statement about which
   earlier steps repeat the current vertex (`row_first_visit`, `col_first_visit`), so it only
   reads the pattern. The last position is a down-letter (`letterAt_last_D`). Once the word is
   recovered, the two labelings are recovered by evaluation.

The labelings of one word `w` are the injections out of the row set and the column set of
`dyckWalk w`, two types that depend on `w`, so the source of the injection is a
`Finset.sigma` and not a product. That is what makes the count exact: the fiber over `w` has
`n.descFactorial s * d.descFactorial t` elements by `Fintype.card_embedding_eq`. -/

/-! ## The elementary bounds on the labels -/

/-- The row label at the step `x` is at most `x`: the vertex was opened at or before `2 x`. -/
theorem dyckWalk_fst_le {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (x : Fin k) : ((dyckWalk w hw hn hd).1 x).1 ≤ x.1 := by
  have h := vid_le w.toList (2 * x.1 + 1)
  simp only [dyckWalk_fst_val]
  omega

/-- The column label at the step `x` is at most `x`. -/
theorem dyckWalk_snd_le {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (x : Fin k) : ((dyckWalk w hw hn hd).2 x).1 ≤ x.1 := by
  have h := vid_le w.toList (2 * x.1)
  simp only [dyckWalk_snd_val]
  omega

/-- The vertex name at the odd position `2 x + 1`, read from the row label. -/
theorem vid_odd_eq {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n) (hd : k ≤ d)
    (x : Fin k) : vid w.toList (2 * x.1 + 1) = 2 * ((dyckWalk w hw hn hd).1 x).1 + 1 := by
  have hx := x.2
  have hs : 2 * x.1 + 1 ≤ 2 * w.semilength := by omega
  have h := vid_eq_two_mul_label w hs
  simp only [dyckWalk_fst_val]
  omega

/-- The vertex name at the even position `2 x`, read from the column label. -/
theorem vid_even_eq {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n) (hd : k ≤ d)
    (x : Fin k) : vid w.toList (2 * x.1) = 2 * ((dyckWalk w hw hn hd).2 x).1 := by
  have hx := x.2
  have hs : 2 * x.1 ≤ 2 * w.semilength := by omega
  have h := vid_eq_two_mul_label w hs
  simp only [dyckWalk_snd_val]
  omega

/-! ## First visits -/

/-- The row of the step `x` is new exactly when its label is `x`. The label of a row is the
index of the step that opens it, so a repeated row carries a strictly smaller label. -/
theorem row_first_visit {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (x : Fin k) :
    ((dyckWalk w hw hn hd).1 x).1 = x.1 ↔
      ∀ y : Fin k, y.1 < x.1 → (dyckWalk w hw hn hd).1 y ≠ (dyckWalk w hw hn hd).1 x := by
  constructor
  · intro h y hy hcon
    have hle := dyckWalk_fst_le w hw hn hd y
    rw [hcon, h] at hle
    omega
  · intro h
    by_contra hne
    have hle := dyckWalk_fst_le w hw hn hd x
    set m := ((dyckWalk w hw hn hd).1 x).1 with hm
    have hmx : m < x.1 := by omega
    have hmk : m < k := by have := x.2; omega
    have hvid := vid_odd_eq w hw hn hd x
    rcases hc : dstack w.toList (2 * x.1 + 1) with _ | ⟨p, rest⟩
    · rw [vid_of_nil hc] at hvid; omega
    · have hp1 : p + 1 = 2 * m + 1 := by rw [← vid_of_cons hc]; exact hvid
      have hpU : letterAt w.toList p = U :=
        dstack_mem_U w.toList (2 * x.1 + 1) p (by rw [hc]; exact List.mem_cons_self ..)
      have hvm : vid w.toList (p + 1) = p + 1 := (letterAt_eq_U_iff _ _).mp hpU
      have hpm : p = 2 * m := by omega
      rw [hpm] at hvm
      have h2 := vid_odd_eq w hw hn hd (⟨m, hmk⟩ : Fin k)
      simp only [] at h2
      refine h ⟨m, hmk⟩ hmx (Fin.ext ?_)
      omega

/-- The column of the step `x` is new exactly when its label is `x`. -/
theorem col_first_visit {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (x : Fin k) :
    ((dyckWalk w hw hn hd).2 x).1 = x.1 ↔
      ∀ y : Fin k, y.1 < x.1 → (dyckWalk w hw hn hd).2 y ≠ (dyckWalk w hw hn hd).2 x := by
  constructor
  · intro h y hy hcon
    have hle := dyckWalk_snd_le w hw hn hd y
    rw [hcon, h] at hle
    omega
  · intro h
    by_contra hne
    have hle := dyckWalk_snd_le w hw hn hd x
    set m := ((dyckWalk w hw hn hd).2 x).1 with hm
    have hmx : m < x.1 := by omega
    have hmk : m < k := by have := x.2; omega
    have hvid := vid_even_eq w hw hn hd x
    rcases Nat.eq_zero_or_pos m with hm0 | hm1
    · have hk0 : 0 < k := by omega
      have h0 := vid_even_eq w hw hn hd (⟨0, hk0⟩ : Fin k)
      have hz : (⟨0, hk0⟩ : Fin k).1 = 0 := rfl
      rw [hz, Nat.mul_zero, vid_zero] at h0
      refine h ⟨0, hk0⟩ (by rw [hz]; omega) (Fin.ext ?_)
      omega
    · rcases hc : dstack w.toList (2 * x.1) with _ | ⟨p, rest⟩
      · rw [vid_of_nil hc] at hvid; omega
      · have hp1 : p + 1 = 2 * m := by rw [← vid_of_cons hc]; exact hvid
        have hpU : letterAt w.toList p = U :=
          dstack_mem_U w.toList (2 * x.1) p (by rw [hc]; exact List.mem_cons_self ..)
        have hvm : vid w.toList (p + 1) = p + 1 := (letterAt_eq_U_iff _ _).mp hpU
        have hpm : p + 1 = 2 * m := hp1
        rw [hpm] at hvm
        have h2 := vid_even_eq w hw hn hd (⟨m, hmk⟩ : Fin k)
        simp only [] at h2
        refine h ⟨m, hmk⟩ hmx (Fin.ext ?_)
        omega

/-! ## The word, read back from the walk -/

/-- An even position carries an up-letter exactly when the row of that step is new. -/
theorem letterAt_even_iff {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (x : Fin k) :
    letterAt w.toList (2 * x.1) = U ↔ ((dyckWalk w hw hn hd).1 x).1 = x.1 := by
  rw [letterAt_eq_U_iff]
  have h := vid_odd_eq w hw hn hd x
  omega

/-- An odd position `2 y - 1` carries an up-letter exactly when the column of the step `y` is
new. -/
theorem letterAt_odd_iff {k n d : ℕ} (w : DyckWord) (hw : w.semilength = k) (hn : k ≤ n)
    (hd : k ≤ d) (y : Fin k) (hy : 1 ≤ y.1) :
    letterAt w.toList (2 * y.1 - 1) = U ↔ ((dyckWalk w hw hn hd).2 y).1 = y.1 := by
  rw [letterAt_eq_U_iff]
  have he : 2 * y.1 - 1 + 1 = 2 * y.1 := by omega
  rw [he]
  have h := vid_even_eq w hw hn hd y
  omega

/-- The last position of a word of semilength `k ≥ 1` carries a down-letter. -/
theorem letterAt_last_D {k : ℕ} (w : DyckWord) (hk : 1 ≤ k) (hw : w.semilength = k) :
    letterAt w.toList (2 * k - 1) = D := by
  cases hc : letterAt w.toList (2 * k - 1) with
  | D => rfl
  | U =>
    exfalso
    have h := (letterAt_eq_U_iff w.toList (2 * k - 1)).mp hc
    have he := vid_end w
    rw [hw] at he
    have h2 : 2 * k - 1 + 1 = 2 * k := by omega
    rw [h2, he] at h
    omega

/-- The word is recoverable from the repetition pattern of its contour walk. Two words of
semilength `k` whose walks repeat the same steps are equal. -/
theorem dyckWord_eq_of_pattern {k n d : ℕ} {w₁ w₂ : DyckWord}
    (hw₁ : w₁.semilength = k) (hw₂ : w₂.semilength = k) (hn : k ≤ n) (hd : k ≤ d)
    (hi : ∀ x y : Fin k, (dyckWalk w₁ hw₁ hn hd).1 x = (dyckWalk w₁ hw₁ hn hd).1 y ↔
        (dyckWalk w₂ hw₂ hn hd).1 x = (dyckWalk w₂ hw₂ hn hd).1 y)
    (hj : ∀ x y : Fin k, (dyckWalk w₁ hw₁ hn hd).2 x = (dyckWalk w₁ hw₁ hn hd).2 y ↔
        (dyckWalk w₂ hw₂ hn hd).2 x = (dyckWalk w₂ hw₂ hn hd).2 y) : w₁ = w₂ := by
  have hlen₁ : w₁.toList.length = 2 * k := by rw [← w₁.two_mul_semilength_eq_length, hw₁]
  have hlen₂ : w₂.toList.length = 2 * k := by rw [← w₂.two_mul_semilength_eq_length, hw₂]
  have hU : ∀ p, p < 2 * k → (letterAt w₁.toList p = U ↔ letterAt w₂.toList p = U) := by
    intro p hp
    have hpar : p % 2 = 0 ∨ p % 2 = 1 := by omega
    rcases hpar with hpar | hpar
    · have hxk : p / 2 < k := by omega
      have hpx : 2 * (⟨p / 2, hxk⟩ : Fin k).1 = p := by simp only []; omega
      rw [← hpx, letterAt_even_iff w₁ hw₁ hn hd, letterAt_even_iff w₂ hw₂ hn hd,
        row_first_visit w₁ hw₁ hn hd, row_first_visit w₂ hw₂ hn hd]
      exact ⟨fun hA y hy hcon => hA y hy ((hi y _).mpr hcon),
        fun hA y hy hcon => hA y hy ((hi y _).mp hcon)⟩
    · rcases Nat.lt_or_ge ((p + 1) / 2) k with hyk | hyk
      · have hpy : 2 * (⟨(p + 1) / 2, hyk⟩ : Fin k).1 - 1 = p := by simp only []; omega
        have hy1 : 1 ≤ (⟨(p + 1) / 2, hyk⟩ : Fin k).1 := by simp only []; omega
        rw [← hpy, letterAt_odd_iff w₁ hw₁ hn hd _ hy1, letterAt_odd_iff w₂ hw₂ hn hd _ hy1,
          col_first_visit w₁ hw₁ hn hd, col_first_visit w₂ hw₂ hn hd]
        exact ⟨fun hA y hy hcon => hA y hy ((hj y _).mpr hcon),
          fun hA y hy hcon => hA y hy ((hj y _).mp hcon)⟩
      · have hk1 : 1 ≤ k := by omega
        have hplast : p = 2 * k - 1 := by omega
        rw [hplast, letterAt_last_D w₁ hk1 hw₁, letterAt_last_D w₂ hk1 hw₂]
  rw [DyckWord.ext_iff]
  refine List.ext_getElem (by omega) ?_
  intro p h₁ h₂
  have e₁ : w₁.toList[p] = letterAt w₁.toList p := by
    simp [letterAt, List.getElem?_eq_getElem h₁]
  have e₂ : w₂.toList[p] = letterAt w₂.toList p := by
    simp [letterAt, List.getElem?_eq_getElem h₂]
  have hiff := hU p (by omega)
  rw [e₁, e₂]
  cases h1 : letterAt w₁.toList p with
  | U => rw [h1] at hiff; exact (hiff.mp rfl).symm
  | D =>
    cases h2 : letterAt w₂.toList p with
    | D => rfl
    | U => rw [h1, h2] at hiff; exact absurd (hiff.mpr rfl) (by simp)

/-! ## A relabeling keeps the walk in its cell -/

/-- Two walks that repeat the same steps have the same multiplicity profile, so a walk with no
entry of multiplicity 1 stays one under any relabeling that is injective on the labels the
walk uses. -/
private theorem pattern_filter_eq {k n d n' d' : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    {i' : Fin k → Fin n'} {j' : Fin k → Fin d'}
    (hi : ∀ x y : Fin k, i' x = i' y ↔ i x = i y)
    (hj : ∀ x y : Fin k, j' x = j' y ↔ j x = j y)
    {e : Fin n' × Fin d'} {x u : Fin k} (hx : i' x = e.1) (hu : j' u = e.2) :
    ({t : Fin k | i' t = e.1 ∧ j' t = e.2} : Finset (Fin k))
        = ({t : Fin k | i t = i x ∧ j t = j u} : Finset (Fin k))
      ∧ ({t : Fin k | i' t = e.1 ∧ j' (cycSucc t) = e.2} : Finset (Fin k))
        = ({t : Fin k | i t = i x ∧ j (cycSucc t) = j u} : Finset (Fin k)) := by
  classical
  constructor
  · ext t
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, ← hx, ← hu]
    rw [hi t x, hj t u]
  · ext t
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, ← hx, ← hu]
    rw [hi t x, hj (cycSucc t) u]

/-- `NoSingle` only reads the repetition pattern of a walk. -/
theorem noSingle_of_pattern {k n d n' d' : ℕ} {i : Fin k → Fin n} {j : Fin k → Fin d}
    {i' : Fin k → Fin n'} {j' : Fin k → Fin d'} (h : NoSingle i j)
    (hi : ∀ x y : Fin k, i' x = i' y ↔ i x = i y)
    (hj : ∀ x y : Fin k, j' x = j' y ↔ j x = j y) : NoSingle i' j' := by
  classical
  intro e he
  rw [walkMult] at he
  by_cases hA : (({t : Fin k | i' t = e.1 ∧ j' t = e.2} : Finset (Fin k))).Nonempty
  · obtain ⟨x, hx⟩ := hA
    rw [Finset.mem_filter] at hx
    obtain ⟨-, hx1, hx2⟩ := hx
    obtain ⟨E1, E2⟩ := pattern_filter_eq hi hj hx1 hx2
    refine h (i x, j x) ?_
    simp only [walkMult, ← E1, ← E2]
    exact he
  · rw [Finset.not_nonempty_iff_eq_empty] at hA
    rw [hA, Finset.card_empty] at he
    have hB : (({t : Fin k | i' t = e.1 ∧ j' (cycSucc t) = e.2} : Finset (Fin k))).Nonempty :=
      Finset.card_pos.mp (by omega)
    obtain ⟨x, hx⟩ := hB
    rw [Finset.mem_filter] at hx
    obtain ⟨-, hx1, hx2⟩ := hx
    obtain ⟨E1, E2⟩ := pattern_filter_eq hi hj hx1 hx2
    refine h (i x, j (cycSucc x)) ?_
    simp only [walkMult, ← E1, ← E2]
    rw [hA, Finset.card_empty]
    exact he

/-! ## The labeled walk -/

/-- The rows the contour walk of `w` uses. -/
abbrev dyckRowType {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d)
    (w : {w : DyckWord // w.semilength = k}) :=
  {a : Fin n // a ∈ univ.image (dyckWalk w.1 w.2 hn hd).1}

/-- The columns the contour walk of `w` uses. -/
abbrev dyckColType {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d)
    (w : {w : DyckWord // w.semilength = k}) :=
  {b : Fin d // b ∈ univ.image (dyckWalk w.1 w.2 hn hd).2}

/-- A word of semilength `k`, an injective labeling of the rows of its contour walk by
`Fin n`, and an injective labeling of the columns by `Fin d`. -/
abbrev DyckLabel {k : ℕ} (n d : ℕ) (hn : k ≤ n) (hd : k ≤ d) :=
  Σ w : {w : DyckWord // w.semilength = k},
    (dyckRowType hn hd w ↪ Fin n) × (dyckColType hn hd w ↪ Fin d)

/-- The relabeled contour walk. -/
def labelWalk {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d) (p : DyckLabel n d hn hd) :
    (Fin k → Fin n) × (Fin k → Fin d) :=
  (fun x => p.2.1 ⟨(dyckWalk p.1.1 p.1.2 hn hd).1 x, mem_image_of_mem _ (mem_univ x)⟩,
   fun x => p.2.2 ⟨(dyckWalk p.1.1 p.1.2 hn hd).2 x, mem_image_of_mem _ (mem_univ x)⟩)

/-- The relabeling repeats exactly the rows the walk repeats. -/
theorem labelWalk_fst_iff {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d) (p : DyckLabel n d hn hd)
    (x y : Fin k) : (labelWalk hn hd p).1 x = (labelWalk hn hd p).1 y ↔
      (dyckWalk p.1.1 p.1.2 hn hd).1 x = (dyckWalk p.1.1 p.1.2 hn hd).1 y := by
  constructor
  · intro h
    exact congrArg Subtype.val (p.2.1.injective h)
  · intro h
    exact congrArg (fun a : dyckRowType hn hd p.1 => p.2.1 a) (Subtype.ext h)

/-- The relabeling repeats exactly the columns the walk repeats. -/
theorem labelWalk_snd_iff {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d) (p : DyckLabel n d hn hd)
    (x y : Fin k) : (labelWalk hn hd p).2 x = (labelWalk hn hd p).2 y ↔
      (dyckWalk p.1.1 p.1.2 hn hd).2 x = (dyckWalk p.1.1 p.1.2 hn hd).2 y := by
  constructor
  · intro h
    exact congrArg Subtype.val (p.2.2.injective h)
  · intro h
    exact congrArg (fun b : dyckColType hn hd p.1 => p.2.2 b) (Subtype.ext h)

/-- The image of a relabeling has as many rows as the walk. -/
theorem card_image_labelWalk_fst {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d)
    (p : DyckLabel n d hn hd) :
    (univ.image (labelWalk hn hd p).1).card
      = (univ.image (dyckWalk p.1.1 p.1.2 hn hd).1).card := by
  classical
  have hg : (univ : Finset (Fin k)).image
      (fun x => (⟨(dyckWalk p.1.1 p.1.2 hn hd).1 x, mem_image_of_mem _ (mem_univ x)⟩ :
        dyckRowType hn hd p.1)) = univ := by
    apply Finset.eq_univ_iff_forall.mpr
    rintro ⟨a, ha⟩
    obtain ⟨x, -, hx⟩ := Finset.mem_image.mp ha
    exact Finset.mem_image.mpr ⟨x, Finset.mem_univ x, Subtype.ext hx⟩
  have : (univ.image (labelWalk hn hd p).1)
      = ((univ : Finset (Fin k)).image
          (fun x => (⟨(dyckWalk p.1.1 p.1.2 hn hd).1 x, mem_image_of_mem _ (mem_univ x)⟩ :
            dyckRowType hn hd p.1))).image p.2.1 := by
    rw [Finset.image_image]; rfl
  rw [this, hg, Finset.card_image_of_injective _ p.2.1.injective, Finset.card_univ,
    Fintype.card_coe]

/-- The image of a relabeling has as many columns as the walk. -/
theorem card_image_labelWalk_snd {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d)
    (p : DyckLabel n d hn hd) :
    (univ.image (labelWalk hn hd p).2).card
      = (univ.image (dyckWalk p.1.1 p.1.2 hn hd).2).card := by
  classical
  have hg : (univ : Finset (Fin k)).image
      (fun x => (⟨(dyckWalk p.1.1 p.1.2 hn hd).2 x, mem_image_of_mem _ (mem_univ x)⟩ :
        dyckColType hn hd p.1)) = univ := by
    apply Finset.eq_univ_iff_forall.mpr
    rintro ⟨b, hb⟩
    obtain ⟨x, -, hx⟩ := Finset.mem_image.mp hb
    exact Finset.mem_image.mpr ⟨x, Finset.mem_univ x, Subtype.ext hx⟩
  have : (univ.image (labelWalk hn hd p).2)
      = ((univ : Finset (Fin k)).image
          (fun x => (⟨(dyckWalk p.1.1 p.1.2 hn hd).2 x, mem_image_of_mem _ (mem_univ x)⟩ :
            dyckColType hn hd p.1))).image p.2.2 := by
    rw [Finset.image_image]; rfl
  rw [this, hg, Finset.card_image_of_injective _ p.2.2.injective, Finset.card_univ,
    Fintype.card_coe]

/-- The relabeled contour walk is a walk of the cell `(dyckRows w, dyckCols w)`. -/
theorem labelWalk_mem_cellSet {k n d s t : ℕ} (hk : 1 ≤ k) (hn : k ≤ n) (hd : k ≤ d)
    (p : DyckLabel n d hn hd) (hs : dyckRows p.1.1 = s) (ht : dyckCols p.1.1 = t) :
    labelWalk hn hd p ∈ cellSet n d k s t := by
  refine mem_cellSet.mpr ⟨?_, ?_, ?_⟩
  · exact noSingle_of_pattern (noSingle_dyckWalk p.1.1 p.1.2 hn hd)
      (labelWalk_fst_iff hn hd p) (labelWalk_snd_iff hn hd p)
  · rw [card_image_labelWalk_fst, card_image_dyckWalk_fst, hs]
  · rw [card_image_labelWalk_snd, card_image_dyckWalk_snd p.1.1 hk p.1.2 hn hd, ht]

/-- The labeled contour walk determines the word and the two labelings. -/
theorem labelWalk_injective {k n d : ℕ} (hn : k ≤ n) (hd : k ≤ d) :
    Function.Injective (labelWalk (k := k) hn hd) := by
  rintro ⟨w₁, σ₁, τ₁⟩ ⟨w₂, σ₂, τ₂⟩ heq
  have h1 : (labelWalk hn hd ⟨w₁, σ₁, τ₁⟩).1 = (labelWalk hn hd ⟨w₂, σ₂, τ₂⟩).1 :=
    congrArg Prod.fst heq
  have h2 : (labelWalk hn hd ⟨w₁, σ₁, τ₁⟩).2 = (labelWalk hn hd ⟨w₂, σ₂, τ₂⟩).2 :=
    congrArg Prod.snd heq
  have hw : w₁.1 = w₂.1 := by
    refine dyckWord_eq_of_pattern w₁.2 w₂.2 hn hd ?_ ?_
    · intro x y
      rw [← labelWalk_fst_iff hn hd ⟨w₁, σ₁, τ₁⟩, ← labelWalk_fst_iff hn hd ⟨w₂, σ₂, τ₂⟩,
        h1]
    · intro x y
      rw [← labelWalk_snd_iff hn hd ⟨w₁, σ₁, τ₁⟩, ← labelWalk_snd_iff hn hd ⟨w₂, σ₂, τ₂⟩,
        h2]
  obtain rfl : w₁ = w₂ := Subtype.ext hw
  have hσ : σ₁ = σ₂ := by
    apply Function.Embedding.ext
    rintro ⟨a, ha⟩
    obtain ⟨x, -, hx⟩ := Finset.mem_image.mp ha
    have := congrFun h1 x
    simp only [labelWalk] at this
    have hax : (⟨a, ha⟩ : dyckRowType hn hd w₁)
        = ⟨(dyckWalk w₁.1 w₁.2 hn hd).1 x, mem_image_of_mem _ (mem_univ x)⟩ :=
      Subtype.ext hx.symm
    rw [hax]
    exact this
  have hτ : τ₁ = τ₂ := by
    apply Function.Embedding.ext
    rintro ⟨b, hb⟩
    obtain ⟨x, -, hx⟩ := Finset.mem_image.mp hb
    have := congrFun h2 x
    simp only [labelWalk] at this
    have hbx : (⟨b, hb⟩ : dyckColType hn hd w₁)
        = ⟨(dyckWalk w₁.1 w₁.2 hn hd).2 x, mem_image_of_mem _ (mem_univ x)⟩ :=
      Subtype.ext hx.symm
    rw [hbx]
    exact this
  rw [hσ, hτ]

/-! ## The count -/

/-- The Dyck words of semilength `k` whose contour walk has `s` rows, inside the `Fintype` of
the words of semilength `k`. -/
def dyckFinset (k s : ℕ) : Finset {w : DyckWord // w.semilength = k} :=
  univ.filter (fun w => dyckRows w.1 = s)

/-- The Dyck words of semilength `k` with `s` rows. -/
def dyckSet (k s : ℕ) : Finset DyckWord := (dyckFinset k s).image Subtype.val

/-- The two presentations of the word set have the same size. -/
theorem card_dyckSet (k s : ℕ) : (dyckSet k s).card = (dyckFinset k s).card :=
  Finset.card_image_of_injective _ Subtype.val_injective

/-- The source of the injection: every word with `s` rows, with every pair of labelings. -/
noncomputable def dyckLabelSet {k : ℕ} (n d : ℕ) (hn : k ≤ n) (hd : k ≤ d) (s : ℕ) :
    Finset (DyckLabel n d hn hd) :=
  (dyckFinset k s).sigma (fun _ => univ)

/-- The source has `#(Dyck words with s rows) * n.descFactorial s * d.descFactorial t`
elements: the fiber over one word is the pair of embeddings, counted by
`Fintype.card_embedding_eq`. -/
theorem card_dyckLabelSet {k n d s t : ℕ} (hk : 1 ≤ k) (hn : k ≤ n) (hd : k ≤ d)
    (hsum : s + t = k + 1) :
    (dyckLabelSet n d hn hd s).card
      = (dyckSet k s).card * (n.descFactorial s * d.descFactorial t) := by
  classical
  rw [dyckLabelSet, Finset.card_sigma, card_dyckSet]
  have hfib : ∀ w ∈ dyckFinset k s,
      (univ : Finset ((dyckRowType hn hd w ↪ Fin n) × (dyckColType hn hd w ↪ Fin d))).card
        = n.descFactorial s * d.descFactorial t := by
    intro w hw
    rw [dyckFinset, Finset.mem_filter] at hw
    have hrows : dyckRows w.1 = s := hw.2
    have hcols : dyckCols w.1 = t := by
      have := dyckRows_add_dyckCols w.1 w.2
      omega
    rw [Finset.card_univ, Fintype.card_prod, Fintype.card_embedding_eq,
      Fintype.card_embedding_eq, Fintype.card_fin, Fintype.card_fin, Fintype.card_coe,
      Fintype.card_coe, card_image_dyckWalk_fst, card_image_dyckWalk_snd w.1 hk w.2 hn hd,
      hrows, hcols]
  rw [Finset.sum_congr rfl hfib, Finset.sum_const, smul_eq_mul]

/-- L2: the labeled contour walks of the Dyck words with `s` rows inject into the cell
`(s, t)`, so the cell holds at least `#(dyckSet k s) * n.descFactorial s * d.descFactorial t`
walks. -/
theorem dyck_label_card_le {n d k s t : ℕ} (hk : 1 ≤ k) (hsum : s + t = k + 1)
    (hn : k ≤ n) (hd : k ≤ d) :
    (dyckSet k s).card * (n.descFactorial s * d.descFactorial t) ≤ cellCard n d k s t := by
  classical
  rw [← card_dyckLabelSet hk hn hd hsum, cellCard]
  refine Finset.card_le_card_of_injOn (labelWalk hn hd) ?_ (labelWalk_injective hn hd).injOn
  intro p hp
  rw [Finset.mem_coe, dyckLabelSet, Finset.mem_sigma, dyckFinset, Finset.mem_filter] at hp
  have hrows : dyckRows p.1.1 = s := hp.1.2
  have hcols : dyckCols p.1.1 = t := by
    have := dyckRows_add_dyckCols p.1.1 p.1.2
    omega
  exact Finset.mem_coe.mpr (labelWalk_mem_cellSet hk hn hd p hrows hcols)

end Edge
end StackedSVD

/-! ## X8, item L3: the Narayana lower bound on `dyckSet k s`

The proof is `scripts/stage3_count/x8_scratch/L3/l3_proof.md`. Route A, the double rotation:

* `parityPairs k s` is the set of pairs `(x, y)`, where `x ⊆ range k` names the even positions
  that carry an up-letter and `y ⊆ range k` names the odd ones. Its size is `C(k,s)^2`.
* `pairWord k x y a b` interleaves the two cyclic words, `x` rotated by `a` at the even
  positions and `y` rotated by `b` at the odd positions. It always has `s` up-letters at even
  positions, so it lies in `dyckSet k s` as soon as it is a Dyck word.
* `exists_rotation_isDyck`: for `1 ≤ s ≤ k` some rotation `(a, b)` gives a Dyck word. The proof
  is the parity argument of section 4 of `l3_proof.md`: with
  `Psi d n = PX (n+1) + PY (n+d)`, the word of `(a,b)` is Dyck iff `x a = U` and `Psi (b-a) a`
  is the minimum of `Psi (b-a)`; if no minimum ever sat on an up-letter then
  `mu (d+1) = mu d + 1` for every `d`, while `mu (d+k) = mu d + k - 2 s`, forcing `s = 0`.
* A rotation orbit has at most `k^2` members, so `Finset.card_le_mul_card_image` gives
  `C(k,s)^2 ≤ k^2 * (dyckSet k s).card`.
* `C(k,s-1) ≤ k * C(k,s)` then gives the target with the `k^3` the plan allows.
-/

namespace StackedSVD
namespace Edge

open Finset DyckStep

/-- The letter of the cyclic word named by `x ⊆ range k` at the index `i`, read modulo `k`. -/
def rotIn (k : ℕ) (x : Finset ℕ) (i : ℕ) : DyckStep := if i % k ∈ x then U else D

/-- The word of the pair `(x, y)` rotated by `(a, b)`: the even positions carry `x` from `a`,
the odd positions carry `y` from `b`. -/
def pairWord (k : ℕ) (x y : Finset ℕ) (a b : ℕ) : List DyckStep :=
  (List.range k).flatMap (fun m => [rotIn k x (a + m), rotIn k y (b + m)])

/-- The pairs of up-position sets: `s` of the `k` even positions, `k - s` of the odd ones. -/
def parityPairs (k s : ℕ) : Finset (Finset ℕ × Finset ℕ) :=
  (Finset.powersetCard s (range k)) ×ˢ (Finset.powersetCard (k - s) (range k))

/-- The source of the map has `C(k,s)^2` elements. -/
theorem card_parityPairs {k s : ℕ} (hsk : s ≤ k) :
    (parityPairs k s).card = k.choose s ^ 2 := by
  rw [parityPairs, Finset.card_product, Finset.card_powersetCard, Finset.card_powersetCard,
    Finset.card_range, Nat.choose_symm hsk, sq]

/-- The word of a pair has length `2 * k`. -/
theorem pairWord_length (k : ℕ) (x y : Finset ℕ) (a b : ℕ) :
    (pairWord k x y a b).length = 2 * k := by
  simp [pairWord, Nat.mul_comm]

/-! ### The letters of an interleaved word -/

/-- The length of a word interleaved from two letter functions. -/
private theorem length_interleave (f g : ℕ → DyckStep) (K : ℕ) :
    ((List.range K).flatMap (fun t => [f t, g t])).length = 2 * K := by
  simp [Nat.mul_comm]

/-- The even positions of an interleaved word carry the first letter function. -/
private theorem letterAt_interleave_even (f g : ℕ → DyckStep) :
    ∀ K m : ℕ, m < K →
      letterAt ((List.range K).flatMap (fun t => [f t, g t])) (2 * m) = f m := by
  intro K
  induction K with
  | zero => intro m hm; exact absurd hm (Nat.not_lt_zero m)
  | succ K ih =>
    intro m hm
    have hlen : ((List.range K).flatMap (fun t => [f t, g t])).length = 2 * K :=
      length_interleave f g K
    rw [List.range_succ, List.flatMap_append]
    rcases Nat.lt_or_ge m K with h | h
    · have h2 : 2 * m < ((List.range K).flatMap (fun t => [f t, g t])).length := by omega
      simp only [letterAt, List.getElem?_append_left h2]
      exact ih m h
    · obtain rfl : m = K := by omega
      have h2 : ((List.range m).flatMap (fun t => [f t, g t])).length ≤ 2 * m := by omega
      simp only [letterAt, List.getElem?_append_right h2, hlen]
      simp

/-- The odd positions of an interleaved word carry the second letter function. -/
private theorem letterAt_interleave_odd (f g : ℕ → DyckStep) :
    ∀ K m : ℕ, m < K →
      letterAt ((List.range K).flatMap (fun t => [f t, g t])) (2 * m + 1) = g m := by
  intro K
  induction K with
  | zero => intro m hm; exact absurd hm (Nat.not_lt_zero m)
  | succ K ih =>
    intro m hm
    have hlen : ((List.range K).flatMap (fun t => [f t, g t])).length = 2 * K :=
      length_interleave f g K
    rw [List.range_succ, List.flatMap_append]
    rcases Nat.lt_or_ge m K with h | h
    · have h2 : 2 * m + 1 < ((List.range K).flatMap (fun t => [f t, g t])).length := by omega
      simp only [letterAt, List.getElem?_append_left h2]
      exact ih m h
    · obtain rfl : m = K := by omega
      have h2 : ((List.range m).flatMap (fun t => [f t, g t])).length ≤ 2 * m + 1 := by omega
      simp only [letterAt, List.getElem?_append_right h2, hlen]
      simp

/-- The even position `2 m` of the word of a pair carries the rotated `x`. -/
theorem letterAt_pairWord_even (k : ℕ) (x y : Finset ℕ) (a b : ℕ) {m : ℕ} (hm : m < k) :
    letterAt (pairWord k x y a b) (2 * m) = rotIn k x (a + m) :=
  letterAt_interleave_even _ _ k m hm

/-- The odd position `2 m + 1` of the word of a pair carries the rotated `y`. -/
theorem letterAt_pairWord_odd (k : ℕ) (x y : Finset ℕ) (a b : ℕ) {m : ℕ} (hm : m < k) :
    letterAt (pairWord k x y a b) (2 * m + 1) = rotIn k y (b + m) :=
  letterAt_interleave_odd _ _ k m hm

/-- A rotation permutes `range k`, so the rotated membership test counts `x` itself. -/
theorem card_filter_rotIn {k : ℕ} (x : Finset ℕ) (hx : x ⊆ range k) (a : ℕ) :
    ((range k).filter (fun m => rotIn k x (a + m) = U)).card = x.card := by
  classical
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk
    rw [Finset.range_zero, Finset.subset_empty] at hx
    subst hx
    simp
  · have hU : ∀ m : ℕ, (rotIn k x (a + m) = U) ↔ ((a + m) % k ∈ x) := by
      intro m
      by_cases h : (a + m) % k ∈ x <;> simp [rotIn, h]
    have hinj : ∀ p q : ℕ, p < k → q < k → (a + p) % k = (a + q) % k → p = q := by
      intro p q hp hq hpq
      have h1 : p ≡ q [MOD k] := Nat.ModEq.add_left_cancel' a hpq
      have h2 : p % k = q % k := h1
      rwa [Nat.mod_eq_of_lt hp, Nat.mod_eq_of_lt hq] at h2
    have himg : (range k).image (fun m => (a + m) % k) = range k := by
      apply Finset.eq_of_subset_of_card_le
      · intro j hj
        rw [Finset.mem_image] at hj
        obtain ⟨m, -, rfl⟩ := hj
        exact Finset.mem_range.mpr (Nat.mod_lt _ hk)
      · refine le_of_eq (Finset.card_image_of_injOn ?_).symm
        intro p hp q hq hpq
        exact hinj p q (by simpa using hp) (by simpa using hq) hpq
    refine Finset.card_bij (fun m _ => (a + m) % k) ?_ ?_ ?_
    · intro m hm
      rw [Finset.mem_filter] at hm
      exact (hU m).mp hm.2
    · intro p hp q hq hpq
      rw [Finset.mem_filter, Finset.mem_range] at hp hq
      exact hinj p q hp.1 hq.1 hpq
    · intro j hj
      have hjk : j ∈ (range k).image (fun m => (a + m) % k) := by rw [himg]; exact hx hj
      rw [Finset.mem_image] at hjk
      obtain ⟨m, hm, hmj⟩ := hjk
      refine ⟨m, ?_, hmj⟩
      rw [Finset.mem_filter]
      exact ⟨hm, (hU m).mpr (by rw [hmj]; exact hj)⟩

/-- The word of a pair has `s` up-letters at even positions, whatever the rotation: the even
positions of `pairWord k x y a b` carry one full period of the cyclic word `x`. -/
theorem pairWord_rows {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (a b : ℕ) :
    ((range (2 * k)).filter
        (fun i => i % 2 = 0 ∧ letterAt (pairWord k x y a b) i = U)).card = s := by
  classical
  have step1 : ((range (2 * k)).filter
      (fun i => i % 2 = 0 ∧ letterAt (pairWord k x y a b) i = U)).card
      = ((range k).filter (fun m => rotIn k x (a + m) = U)).card := by
    refine Finset.card_bij (fun i _ => i / 2) ?_ ?_ ?_
    · intro i hi
      rw [Finset.mem_filter, Finset.mem_range] at hi
      obtain ⟨hi1, hi2, hi3⟩ := hi
      rw [Finset.mem_filter, Finset.mem_range]
      refine ⟨by omega, ?_⟩
      have he : 2 * (i / 2) = i := by omega
      rw [← letterAt_pairWord_even k x y a b (m := i / 2) (by omega), he]
      exact hi3
    · intro p hp q hq hpq
      rw [Finset.mem_filter, Finset.mem_range] at hp hq
      obtain ⟨-, hp2, -⟩ := hp
      obtain ⟨-, hq2, -⟩ := hq
      omega
    · intro m hm
      rw [Finset.mem_filter, Finset.mem_range] at hm
      obtain ⟨hm1, hm2⟩ := hm
      refine ⟨2 * m, ?_, by omega⟩
      rw [Finset.mem_filter, Finset.mem_range]
      refine ⟨by omega, by omega, ?_⟩
      rw [letterAt_pairWord_even k x y a b hm1]
      exact hm2
  rw [step1, ← hxc]
  exact card_filter_rotIn x hx a

/-! ### The height of a prefix, and the Dyck test -/

/-- The height of `L` after `i` letters: the up-letters minus the down-letters. -/
def hgt (L : List DyckStep) (i : ℕ) : ℤ :=
  ((L.take i).count U : ℤ) - ((L.take i).count D : ℤ)

/-- The height starts at `0`. -/
theorem hgt_zero (L : List DyckStep) : hgt L 0 = 0 := by simp [hgt]

/-- An up-letter raises the height by one. -/
theorem hgt_succ_U (L : List DyckStep) {i : ℕ} (hi : i < L.length) (h : letterAt L i = U) :
    hgt L (i + 1) = hgt L i + 1 := by
  have htake : L.take (i + 1) = L.take i ++ [U] := by
    rw [take_succ_eq _ (getElem?_eq_letterAt hi), h]
  rw [hgt, hgt, htake, List.count_append, List.count_append, count_U_singleton_U,
    count_D_singleton_U]
  push_cast
  ring

/-- A down-letter lowers the height by one. -/
theorem hgt_succ_D (L : List DyckStep) {i : ℕ} (hi : i < L.length) (h : letterAt L i = D) :
    hgt L (i + 1) = hgt L i - 1 := by
  have htake : L.take (i + 1) = L.take i ++ [D] := by
    rw [take_succ_eq _ (getElem?_eq_letterAt hi), h]
  rw [hgt, hgt, htake, List.count_append, List.count_append, count_U_singleton_D,
    count_D_singleton_D]
  push_cast
  ring

/-- Past the end of the list the height does not move. -/
theorem hgt_of_length_le (L : List DyckStep) {i : ℕ} (hi : L.length ≤ i) :
    hgt L i = hgt L L.length := by
  rw [hgt, hgt, List.take_of_length_le hi, List.take_of_length_le (le_refl L.length)]

/-- A list with nonnegative prefix heights and total height `0` is a Dyck word. -/
theorem exists_dyckWord_of_hgt {L : List DyckStep} (hend : hgt L L.length = 0)
    (hnn : ∀ i, 0 ≤ hgt L i) : ∃ w : DyckWord, w.toList = L := by
  refine ⟨⟨L, ?_, ?_⟩, rfl⟩
  · rw [hgt, List.take_of_length_le (le_refl L.length)] at hend
    omega
  · intro i
    have h := hnn i
    rw [hgt] at h
    omega

/-- The word of a pair has `k - s` up-letters at odd positions, whatever the rotation. -/
theorem pairWord_cols {k s : ℕ} {y : Finset ℕ} (hy : y ⊆ range k) (hyc : y.card = k - s)
    (x : Finset ℕ) (a b : ℕ) :
    ((range (2 * k)).filter
        (fun i => i % 2 = 1 ∧ letterAt (pairWord k x y a b) i = U)).card = k - s := by
  classical
  have step1 : ((range (2 * k)).filter
      (fun i => i % 2 = 1 ∧ letterAt (pairWord k x y a b) i = U)).card
      = ((range k).filter (fun m => rotIn k y (b + m) = U)).card := by
    refine Finset.card_bij (fun i _ => i / 2) ?_ ?_ ?_
    · intro i hi
      rw [Finset.mem_filter, Finset.mem_range] at hi
      obtain ⟨hi1, hi2, hi3⟩ := hi
      rw [Finset.mem_filter, Finset.mem_range]
      refine ⟨by omega, ?_⟩
      have he : 2 * (i / 2) + 1 = i := by omega
      rw [← letterAt_pairWord_odd k x y a b (m := i / 2) (by omega), he]
      exact hi3
    · intro p hp q hq hpq
      rw [Finset.mem_filter, Finset.mem_range] at hp hq
      obtain ⟨-, hp2, -⟩ := hp
      obtain ⟨-, hq2, -⟩ := hq
      omega
    · intro m hm
      rw [Finset.mem_filter, Finset.mem_range] at hm
      obtain ⟨hm1, hm2⟩ := hm
      refine ⟨2 * m + 1, ?_, by omega⟩
      rw [Finset.mem_filter, Finset.mem_range]
      refine ⟨by omega, by omega, ?_⟩
      rw [letterAt_pairWord_odd k x y a b hm1]
      exact hm2
  rw [step1, ← hyc]
  exact card_filter_rotIn y hy b

/-- The word of a pair has `k` up-letters in all. -/
theorem pairWord_count_U {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) (hsk : s ≤ k) (a b : ℕ) :
    (pairWord k x y a b).count U = k := by
  classical
  have hlen : (pairWord k x y a b).length = 2 * k := pairWord_length k x y a b
  have hcount := card_filter_letterAt (pairWord k x y a b) U
  rw [hlen] at hcount
  have hpoint : ∀ i ∈ range (2 * k),
      ((if (i % 2 = 0 ∧ letterAt (pairWord k x y a b) i = U) then 1 else 0)
        + (if (i % 2 = 1 ∧ letterAt (pairWord k x y a b) i = U) then 1 else 0))
      = (if letterAt (pairWord k x y a b) i = U then 1 else 0) := by
    intro i _
    have hp : i % 2 = 0 ∨ i % 2 = 1 := by omega
    by_cases h : letterAt (pairWord k x y a b) i = U
    · rcases hp with h2 | h2 <;> simp [h, h2]
    · simp [h]
  have hsum := Finset.sum_congr rfl hpoint
  rw [Finset.sum_add_distrib] at hsum
  have hr := pairWord_rows (y := y) hx hxc a b
  have hc := pairWord_cols hy hyc x a b
  simp only [Finset.card_filter] at hr hc hcount
  omega

/-- The word of a pair ends at height `0`. -/
theorem pairWord_hgt_end {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) (hsk : s ≤ k) (a b : ℕ) :
    hgt (pairWord k x y a b) (2 * k) = 0 := by
  have hlen : (pairWord k x y a b).length = 2 * k := pairWord_length k x y a b
  have hU : (pairWord k x y a b).count U = k := pairWord_count_U hx hxc hy hyc hsk a b
  have hUD := count_U_add_count_D (pairWord k x y a b)
  rw [hgt, List.take_of_length_le (by omega)]
  omega

/-! ### The partial sums of the two periodic step words -/

/-- The step of the cyclic word named by `x`: `+1` at an up-letter, `-1` at a down-letter. -/
def rotStep (k : ℕ) (x : Finset ℕ) (i : ℕ) : ℤ := if i % k ∈ x then 1 else -1

/-- `PX` of `l3_proof.md`: the partial sum of the steps of the cyclic word named by `x`. -/
def rotSum (k : ℕ) (x : Finset ℕ) (n : ℕ) : ℤ := ∑ i ∈ range n, rotStep k x i

theorem rotSum_zero (k : ℕ) (x : Finset ℕ) : rotSum k x 0 = 0 := by simp [rotSum]

theorem rotSum_succ (k : ℕ) (x : Finset ℕ) (n : ℕ) :
    rotSum k x (n + 1) = rotSum k x n + rotStep k x n := Finset.sum_range_succ _ _

/-- An up-letter has the step `+1`. -/
theorem rotStep_of_U {k : ℕ} {x : Finset ℕ} {i : ℕ} (h : rotIn k x i = U) :
    rotStep k x i = 1 := by
  by_cases hm : i % k ∈ x
  · simp [rotStep, hm]
  · rw [rotIn, if_neg hm] at h; exact absurd h (by decide)

/-- A down-letter has the step `-1`. -/
theorem rotStep_of_D {k : ℕ} {x : Finset ℕ} {i : ℕ} (h : rotIn k x i = D) :
    rotStep k x i = -1 := by
  by_cases hm : i % k ∈ x
  · rw [rotIn, if_pos hm] at h; exact absurd h (by decide)
  · simp [rotStep, hm]

/-- The step word is `k`-periodic. -/
theorem rotStep_add_k (k : ℕ) (x : Finset ℕ) (i : ℕ) :
    rotStep k x (i + k) = rotStep k x i := by
  simp [rotStep, Nat.add_mod_right]

/-- One period of the step word adds `2 |x| - k`. -/
theorem rotSum_period {k : ℕ} {x : Finset ℕ} (hx : x ⊆ range k) :
    rotSum k x k = 2 * (x.card : ℤ) - k := by
  classical
  have h1 : rotSum k x k = ∑ i ∈ range k, (2 * (if i ∈ x then (1 : ℤ) else 0) - 1) := by
    rw [rotSum]
    refine Finset.sum_congr rfl ?_
    intro i hi
    rw [Finset.mem_range] at hi
    rw [rotStep, Nat.mod_eq_of_lt hi]
    by_cases h : i ∈ x <;> simp [h]
  have h2 : ∑ i ∈ range k, (if i ∈ x then (1 : ℤ) else 0) = (x.card : ℤ) := by
    rw [← Finset.sum_filter, Finset.filter_mem_eq_inter, Finset.inter_eq_right.mpr hx]
    simp
  rw [h1, Finset.sum_sub_distrib, ← Finset.mul_sum, h2]
  simp

/-- The partial sum grows by `2 |x| - k` over one period. -/
theorem rotSum_add_k {k : ℕ} {x : Finset ℕ} (hx : x ⊆ range k) (n : ℕ) :
    rotSum k x (n + k) = rotSum k x n + (2 * (x.card : ℤ) - k) := by
  induction n with
  | zero => simpa [rotSum_zero] using rotSum_period hx
  | succ n ih =>
    have h1 : n + 1 + k = (n + k) + 1 := by omega
    rw [h1, rotSum_succ, ih, rotSum_succ, rotStep_add_k]
    ring

/-! ### The height of the word of a pair, read from the partial sums -/

/-- The height after an even prefix. -/
theorem hgt_pairWord_even {k : ℕ} (x y : Finset ℕ) (a b : ℕ) :
    ∀ m : ℕ, m ≤ k →
      hgt (pairWord k x y a b) (2 * m)
        = (rotSum k x (a + m) - rotSum k x a) + (rotSum k y (b + m) - rotSum k y b) := by
  intro m
  induction m with
  | zero => intro _; simp [hgt]
  | succ m ih =>
    intro hm
    have hmk : m < k := by omega
    have hlen : (pairWord k x y a b).length = 2 * k := pairWord_length k x y a b
    have h1 : hgt (pairWord k x y a b) (2 * m + 1)
        = hgt (pairWord k x y a b) (2 * m) + rotStep k x (a + m) := by
      have hlt : 2 * m < (pairWord k x y a b).length := by omega
      have he := letterAt_pairWord_even k x y a b hmk
      rcases (rotIn k x (a + m)).dichotomy with h | h
      · rw [hgt_succ_U _ hlt (he.trans h), rotStep_of_U h]
      · rw [hgt_succ_D _ hlt (he.trans h), rotStep_of_D h]; ring
    have h2 : hgt (pairWord k x y a b) (2 * m + 1 + 1)
        = hgt (pairWord k x y a b) (2 * m + 1) + rotStep k y (b + m) := by
      have hlt : 2 * m + 1 < (pairWord k x y a b).length := by omega
      have ho := letterAt_pairWord_odd k x y a b hmk
      rcases (rotIn k y (b + m)).dichotomy with h | h
      · rw [hgt_succ_U _ hlt (ho.trans h), rotStep_of_U h]
      · rw [hgt_succ_D _ hlt (ho.trans h), rotStep_of_D h]; ring
    have h3 : 2 * (m + 1) = 2 * m + 1 + 1 := by ring
    rw [h3, h2, h1, ih (by omega)]
    have e1 : a + (m + 1) = (a + m) + 1 := by omega
    have e2 : b + (m + 1) = (b + m) + 1 := by omega
    rw [e1, e2, rotSum_succ k x (a + m), rotSum_succ k y (b + m)]
    ring

/-- The height after an odd prefix. -/
theorem hgt_pairWord_odd {k : ℕ} (x y : Finset ℕ) (a b : ℕ) {m : ℕ} (hm : m < k) :
    hgt (pairWord k x y a b) (2 * m + 1)
      = (rotSum k x (a + m + 1) - rotSum k x a) + (rotSum k y (b + m) - rotSum k y b) := by
  have hlen : (pairWord k x y a b).length = 2 * k := pairWord_length k x y a b
  have hlt : 2 * m < (pairWord k x y a b).length := by omega
  have he := letterAt_pairWord_even k x y a b hm
  have hev := hgt_pairWord_even (k := k) x y a b m (by omega)
  have h1 : hgt (pairWord k x y a b) (2 * m + 1)
      = hgt (pairWord k x y a b) (2 * m) + rotStep k x (a + m) := by
    rcases (rotIn k x (a + m)).dichotomy with h | h
    · rw [hgt_succ_U _ hlt (he.trans h), rotStep_of_U h]
    · rw [hgt_succ_D _ hlt (he.trans h), rotStep_of_D h]; ring
  rw [h1, hev, rotSum_succ k x (a + m)]
  ring

/-- The Dyck test only reads the odd positions: a word of length `2 k` that ends at height `0`
and has height at least `1` after every odd prefix has nonnegative heights everywhere. -/
theorem hgt_nonneg_of_odd {L : List DyckStep} {k : ℕ} (hlen : L.length = 2 * k)
    (hend : hgt L (2 * k) = 0) (hodd : ∀ m, m < k → 1 ≤ hgt L (2 * m + 1)) :
    ∀ i, 0 ≤ hgt L i := by
  intro i
  rcases Nat.lt_or_ge i (2 * k) with hi | hi
  · obtain ⟨m, hm⟩ : ∃ m, i = 2 * m ∨ i = 2 * m + 1 := ⟨i / 2, by omega⟩
    rcases hm with rfl | rfl
    · have hmk : m < k := by omega
      have h1 := hodd m hmk
      have hlt : 2 * m < L.length := by omega
      rcases (letterAt L (2 * m)).dichotomy with h | h
      · have h2 := hgt_succ_U L hlt h
        omega
      · have h2 := hgt_succ_D L hlt h
        omega
    · have hmk : m < k := by omega
      have h1 := hodd m hmk
      omega
  · have h1 : hgt L i = hgt L L.length := hgt_of_length_le L (by omega)
    rw [hlen] at h1
    omega

/-! ### `Psi` and its minimum, sections 2 and 3 of `l3_proof.md`

The remaining combinatorial content of L3 is sections 2 to 4 of `l3_proof.md`: for
`1 ≤ s ≤ k` some double rotation raises every odd prefix height to at least `1`. With `PX`,
`PY` the partial sums of the two periodic step words and `Psi d n = PX (n+1) + PY (n+d)`, the
rotation `(a, a + d)` works iff `x a = U` and `Psi d a` is the minimum of `Psi d`. -/

/-- `Psi d n` of `l3_proof.md`: the value the Dyck test of the rotation `(a, a + d)` compares. -/
def psi (k : ℕ) (x y : Finset ℕ) (d n : ℕ) : ℤ := rotSum k x (n + 1) + rotSum k y (n + d)

/-- The minimum of `Psi d` over one period. The window carries `0` as well, so that the
`Finset` is nonempty even at `k = 0`; at `k ≥ 1` the window is exactly `range k`. -/
noncomputable def psiMin (k : ℕ) (x y : Finset ℕ) (d : ℕ) : ℤ :=
  ((insert 0 (range k)).image (fun n => psi k x y d n)).min' (by simp)

theorem psiMin_def (k : ℕ) (x y : Finset ℕ) (d : ℕ) :
    psiMin k x y d = ((insert 0 (range k)).image (fun n => psi k x y d n)).min' (by simp) := rfl

/-- Every step of `Psi d` is a sum of two steps, so it is even. -/
theorem rotStep_eq (k : ℕ) (x : Finset ℕ) (i : ℕ) :
    rotStep k x i = 1 ∨ rotStep k x i = -1 := by
  by_cases h : i % k ∈ x <;> simp [rotStep, h]

/-- Lemma 1(1) of `l3_proof.md`: `Psi d` is `k`-periodic. -/
theorem psi_add_k {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) (hsk : s ≤ k) (d n : ℕ) :
    psi k x y d (n + k) = psi k x y d n := by
  have h1 : rotSum k x (n + 1 + k) = rotSum k x (n + 1) + (2 * (x.card : ℤ) - k) :=
    rotSum_add_k hx _
  have h2 : rotSum k y (n + d + k) = rotSum k y (n + d) + (2 * (y.card : ℤ) - k) :=
    rotSum_add_k hy _
  have hyz : (y.card : ℤ) = (k : ℤ) - s := by rw [hyc, Nat.cast_sub hsk]
  have e1 : n + k + 1 = n + 1 + k := by omega
  have e2 : n + k + d = n + d + k := by omega
  rw [psi, psi, e1, e2, h1, h2, hxc, hyz]
  ring

/-- The period, iterated. -/
theorem psi_add_mul_k {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) (hsk : s ≤ k) (d r q : ℕ) :
    psi k x y d (r + k * q) = psi k x y d r := by
  induction q with
  | zero => simp
  | succ q ih =>
    have e : r + k * (q + 1) = (r + k * q) + k := by ring
    rw [e, psi_add_k hx hxc hy hyc hsk, ih]

/-- The minimum over one period is a minimum over all of `ℕ`. -/
theorem psiMin_le {k s : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) (hsk : s ≤ k) (hk : 0 < k) (d n : ℕ) :
    psiMin k x y d ≤ psi k x y d n := by
  have h1 : psi k x y d (n % k + k * (n / k)) = psi k x y d (n % k) :=
    psi_add_mul_k hx hxc hy hyc hsk d _ _
  rw [Nat.mod_add_div] at h1
  rw [h1, psiMin_def]
  refine Finset.min'_le _ _ ?_
  rw [Finset.mem_image]
  exact ⟨n % k, Finset.mem_insert_of_mem (Finset.mem_range.mpr (Nat.mod_lt _ hk)), rfl⟩

/-- The minimum is attained inside the window. -/
theorem exists_psi_eq_psiMin (k : ℕ) (x y : Finset ℕ) (hk : 0 < k) (d : ℕ) :
    ∃ n, n < k ∧ psi k x y d n = psiMin k x y d := by
  have hmem := Finset.min'_mem ((insert 0 (range k)).image (fun n => psi k x y d n)) (by simp)
  rw [Finset.mem_image] at hmem
  obtain ⟨n, hn, hval⟩ := hmem
  refine ⟨n, ?_, ?_⟩
  · rw [Finset.mem_insert, Finset.mem_range] at hn
    rcases hn with rfl | hn
    · exact hk
    · exact hn
  · rw [psiMin_def]
    exact hval

/-- Lemma 2(1) of `l3_proof.md`: how `Psi` moves when `d` moves. -/
theorem psi_succ_d (k : ℕ) (x y : Finset ℕ) (d n : ℕ) :
    psi k x y (d + 1) n = psi k x y d (n + 1) - rotStep k x (n + 1) := by
  have e : n + (d + 1) = n + 1 + d := by omega
  rw [psi, psi, e, rotSum_succ k x (n + 1)]
  ring

/-- Lemma 2(2) of `l3_proof.md`: one period of `d` lowers `Psi` by `2 s - k`. -/
theorem psi_add_k_d {k s : ℕ} {x y : Finset ℕ} (hy : y ⊆ range k) (hyc : y.card = k - s)
    (hsk : s ≤ k) (d n : ℕ) :
    psi k x y (d + k) n = psi k x y d n - (2 * (s : ℤ) - k) := by
  have h2 : rotSum k y (n + d + k) = rotSum k y (n + d) + (2 * (y.card : ℤ) - k) :=
    rotSum_add_k hy _
  have hyz : (y.card : ℤ) = (k : ℤ) - s := by rw [hyc, Nat.cast_sub hsk]
  have e : n + (d + k) = n + d + k := by omega
  rw [psi, psi, e, h2, hyz]
  ring

/-- Lemma 1(2) of `l3_proof.md`: every value of `Psi d` has the same parity. -/
theorem psi_parity (k : ℕ) (x y : Finset ℕ) (d n : ℕ) :
    psi k x y d n % 2 = psi k x y d 0 % 2 := by
  induction n with
  | zero => rfl
  | succ n ih =>
    have e : n + 1 + d = n + d + 1 := by omega
    have h1 : psi k x y d (n + 1)
        = psi k x y d n + rotStep k x (n + 1) + rotStep k y (n + d) := by
      rw [psi, psi, e, rotSum_succ k x (n + 1), rotSum_succ k y (n + d)]
      ring
    rcases rotStep_eq k x (n + 1) with h2 | h2 <;> rcases rotStep_eq k y (n + d) with h3 | h3 <;>
      omega

/-- Theorem 3 of `l3_proof.md`, the one open step: for `1 ≤ s` some shift `d` has a minimizer
of `Psi d` that sits on an up-letter of `x`. The proof is the parity argument of section 4:
if every minimizer of every `Psi d` sat on a down-letter then `psiMin (d+1) = psiMin d + 1`
for every `d` (Lemma 2(1), `psi_succ_d`, plus `psi_parity`), so `psiMin (0 + k) = psiMin 0 + k`,
while `psi_add_k_d` gives `psiMin (0 + k) = psiMin 0 - (2 s - k)`; hence `s = 0`. -/
theorem exists_psi_min_at_U {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k)
    {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) :
    ∃ d a : ℕ, rotStep k x a = 1 ∧ psi k x y d a = psiMin k x y d := by
  have hk : 0 < k := lt_of_lt_of_le hs hsk
  by_contra hcon
  push Not at hcon
  have hH : ∀ d n : ℕ, psi k x y d n = psiMin k x y d → rotStep k x n = -1 := by
    intro d n hmin
    rcases rotStep_eq k x n with h | h
    · exact absurd hmin (hcon d n h)
    · exact h
  have hpar : ∀ d n : ℕ, psi k x y d n % 2 = psiMin k x y d % 2 := by
    intro d n
    obtain ⟨n1, -, h1⟩ := exists_psi_eq_psiMin k x y hk d
    rw [← h1, psi_parity k x y d n, psi_parity k x y d n1]
  have hgap : ∀ d n : ℕ, psi k x y d n ≠ psiMin k x y d →
      psiMin k x y d + 2 ≤ psi k x y d n := by
    intro d n hne
    have h1 := psiMin_le hx hxc hy hyc hsk hk d n
    have h2 := hpar d n
    omega
  have hstep : ∀ d : ℕ, psiMin k x y (d + 1) = psiMin k x y d + 1 := by
    intro d
    have hge : ∀ n : ℕ, psiMin k x y d + 1 ≤ psi k x y (d + 1) n := by
      intro n
      have he := psi_succ_d k x y d n
      by_cases hcase : psi k x y d (n + 1) = psiMin k x y d
      · have hU := hH d (n + 1) hcase
        omega
      · have h1 := hgap d (n + 1) hcase
        rcases rotStep_eq k x (n + 1) with h2 | h2 <;> omega
    have hle : psiMin k x y (d + 1) ≤ psiMin k x y d + 1 := by
      obtain ⟨n1, hn1, h1⟩ := exists_psi_eq_psiMin k x y hk d
      obtain ⟨n, hn⟩ : ∃ n, n + 1 = n1 + k := ⟨n1 + k - 1, by omega⟩
      have hper : psi k x y d (n1 + k) = psi k x y d n1 := psi_add_k hx hxc hy hyc hsk d n1
      have hmin' : psi k x y d (n + 1) = psiMin k x y d := by rw [hn, hper, h1]
      have hU := hH d (n + 1) hmin'
      have he := psi_succ_d k x y d n
      have hml := psiMin_le hx hxc hy hyc hsk hk (d + 1) n
      omega
    obtain ⟨n2, -, h2⟩ := exists_psi_eq_psiMin k x y hk (d + 1)
    have h3 := hge n2
    omega
  have hiter : ∀ j : ℕ, psiMin k x y j = psiMin k x y 0 + j := by
    intro j
    induction j with
    | zero => simp
    | succ j ih => rw [hstep j, ih]; push_cast; ring
  have hshift : psiMin k x y k = psiMin k x y 0 - (2 * (s : ℤ) - k) := by
    have hA : ∀ n : ℕ, psi k x y k n = psi k x y 0 n - (2 * (s : ℤ) - k) := by
      intro n
      have h := psi_add_k_d (x := x) hy hyc hsk 0 n
      rwa [Nat.zero_add] at h
    obtain ⟨n1, -, h1⟩ := exists_psi_eq_psiMin k x y hk 0
    obtain ⟨n2, -, h2⟩ := exists_psi_eq_psiMin k x y hk k
    have hml1 := psiMin_le hx hxc hy hyc hsk hk k n1
    have hml2 := psiMin_le hx hxc hy hyc hsk hk 0 n2
    have e1 := hA n1
    have e2 := hA n2
    omega
  have hfin := hiter k
  omega

/-- The rotation `(a, a + d)` of Theorem 3 satisfies the Dyck test. -/
theorem exists_min_rotation {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k)
    {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) :
    ∃ a b : ℕ, ∀ m, m < k →
      rotSum k x a + rotSum k y b + 1 ≤ rotSum k x (a + m + 1) + rotSum k y (b + m) := by
  have hk : 0 < k := lt_of_lt_of_le hs hsk
  obtain ⟨d, a, hU, hmin⟩ := exists_psi_min_at_U hs hsk hx hxc hy hyc
  refine ⟨a, a + d, fun m hm => ?_⟩
  have h1 : psi k x y d a = rotSum k x (a + 1) + rotSum k y (a + d) := rfl
  have h2 : psi k x y d (a + m) = rotSum k x (a + m + 1) + rotSum k y (a + m + d) := rfl
  have h3 : psiMin k x y d ≤ psi k x y d (a + m) := psiMin_le hx hxc hy hyc hsk hk d _
  have h4 : rotSum k x (a + 1) = rotSum k x a + rotStep k x a := rotSum_succ k x a
  have e : a + d + m = a + m + d := by omega
  rw [e]
  omega

/-- Sections 2 to 4 of `l3_proof.md`, read back as a statement about the word. -/
theorem exists_rotation_hgt_odd {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k)
    {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) :
    ∃ a b : ℕ, ∀ m, m < k → 1 ≤ hgt (pairWord k x y a b) (2 * m + 1) := by
  obtain ⟨a, b, hab⟩ := exists_min_rotation hs hsk hx hxc hy hyc
  refine ⟨a, b, fun m hm => ?_⟩
  rw [hgt_pairWord_odd x y a b hm]
  have h := hab m hm
  omega

/-- Sections 2 to 4 of `l3_proof.md`, in the form the Dyck test wants. -/
theorem exists_rotation_hgt_nonneg {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k)
    {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) :
    ∃ a b : ℕ, ∀ i : ℕ, 0 ≤ hgt (pairWord k x y a b) i := by
  obtain ⟨a, b, hab⟩ := exists_rotation_hgt_odd hs hsk hx hxc hy hyc
  exact ⟨a, b, hgt_nonneg_of_odd (pairWord_length k x y a b)
    (pairWord_hgt_end hx hxc hy hyc hsk a b) hab⟩

/-- Theorem 3 of `l3_proof.md`, the existence step: for `1 ≤ s ≤ k` some double rotation of the
pair is a Dyck word. This is the only combinatorial content of L3. -/
theorem exists_rotation_isDyck {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k)
    {x y : Finset ℕ} (hx : x ⊆ range k) (hxc : x.card = s)
    (hy : y ⊆ range k) (hyc : y.card = k - s) :
    ∃ a b : ℕ, ∃ w : DyckWord, w.toList = pairWord k x y a b := by
  obtain ⟨a, b, hab⟩ := exists_rotation_hgt_nonneg hs hsk hx hxc hy hyc
  have hend : hgt (pairWord k x y a b) (pairWord k x y a b).length = 0 := by
    rw [pairWord_length k x y a b]
    exact pairWord_hgt_end hx hxc hy hyc hsk a b
  obtain ⟨w, hw⟩ := exists_dyckWord_of_hgt hend hab
  exact ⟨a, b, w, hw⟩

/-! ### The orbit bound -/

/-- A rotation permutes `range k`. -/
private theorem image_rot_range (k a : ℕ) :
    (range k).image (fun m => (a + m) % k) = range k := by
  classical
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk; simp
  · have hinj : ∀ p q : ℕ, p < k → q < k → (a + p) % k = (a + q) % k → p = q := by
      intro p q hp hq hpq
      have h1 : p ≡ q [MOD k] := Nat.ModEq.add_left_cancel' a hpq
      have h2 : p % k = q % k := h1
      rwa [Nat.mod_eq_of_lt hp, Nat.mod_eq_of_lt hq] at h2
    apply Finset.eq_of_subset_of_card_le
    · intro j hj
      rw [Finset.mem_image] at hj
      obtain ⟨m, -, rfl⟩ := hj
      exact Finset.mem_range.mpr (Nat.mod_lt _ hk)
    · refine le_of_eq (Finset.card_image_of_injOn ?_).symm
      intro p hp q hq hpq
      exact hinj p q (by simpa using hp) (by simpa using hq) hpq

/-- The up-positions of the rotated word, rotated back, give `x` again. -/
theorem image_filter_rotIn {k : ℕ} {x : Finset ℕ} (hx : x ⊆ range k) (a : ℕ) :
    ((range k).filter (fun m => rotIn k x (a + m) = U)).image (fun m => (a + m) % k) = x := by
  classical
  have hU : ∀ m : ℕ, (rotIn k x (a + m) = U) ↔ ((a + m) % k ∈ x) := by
    intro m
    by_cases h : (a + m) % k ∈ x <;> simp [rotIn, h]
  ext j
  simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_range]
  constructor
  · rintro ⟨m, ⟨-, hm2⟩, rfl⟩
    exact (hU m).mp hm2
  · intro hj
    have hjr : j ∈ (range k).image (fun m => (a + m) % k) := by
      rw [image_rot_range]; exact hx hj
    rw [Finset.mem_image] at hjr
    obtain ⟨m, hm, hmj⟩ := hjr
    exact ⟨m, ⟨Finset.mem_range.mp hm, (hU m).mpr (by rw [hmj]; exact hj)⟩, hmj⟩

/-- The rotation map only reads `a` modulo `k`. -/
private theorem rot_fun_mod (k a a' : ℕ) (h : a % k = a' % k) :
    (fun m => (a + m) % k) = (fun m => (a' + m) % k) := by
  funext m
  rw [Nat.add_mod, Nat.add_mod a' m, h]

/-- The pair `(x, y)` is recovered from the word and the rotation. -/
theorem recover_pairWord {k : ℕ} {x y : Finset ℕ} (hx : x ⊆ range k) (hy : y ⊆ range k)
    (a b : ℕ) {L : List DyckStep} (hL : L = pairWord k x y a b) :
    ((range k).filter (fun m => letterAt L (2 * m) = U)).image (fun m => (a + m) % k) = x
      ∧ ((range k).filter (fun m => letterAt L (2 * m + 1) = U)).image
          (fun m => (b + m) % k) = y := by
  classical
  subst hL
  constructor
  · have hfx : (range k).filter (fun m => letterAt (pairWord k x y a b) (2 * m) = U)
        = (range k).filter (fun m => rotIn k x (a + m) = U) := by
      apply Finset.filter_congr
      intro m hm
      rw [letterAt_pairWord_even k x y a b (Finset.mem_range.mp hm)]
    rw [hfx]
    exact image_filter_rotIn hx a
  · have hfy : (range k).filter (fun m => letterAt (pairWord k x y a b) (2 * m + 1) = U)
        = (range k).filter (fun m => rotIn k y (b + m) = U) := by
      apply Finset.filter_congr
      intro m hm
      rw [letterAt_pairWord_odd k x y a b (Finset.mem_range.mp hm)]
    rw [hfy]
    exact image_filter_rotIn hy b

/-- A word of semilength `k` with `s` rows lies in `dyckSet k s`. -/
theorem mem_dyckSet_of {k s : ℕ} {w : DyckWord} (hw : w.semilength = k) (hr : dyckRows w = s) :
    w ∈ dyckSet k s := by
  classical
  rw [dyckSet, Finset.mem_image]
  exact ⟨⟨w, hw⟩, by simp [dyckFinset, hr], rfl⟩

/-- The membership test of `parityPairs`, unpacked. -/
theorem mem_parityPairs {k s : ℕ} {p : Finset ℕ × Finset ℕ} (hp : p ∈ parityPairs k s) :
    p.1 ⊆ range k ∧ p.1.card = s ∧ p.2 ⊆ range k ∧ p.2.card = k - s := by
  rw [parityPairs, Finset.mem_product, Finset.mem_powersetCard, Finset.mem_powersetCard] at hp
  exact ⟨hp.1.1, hp.1.2, hp.2.1, hp.2.2⟩

/-- Theorem 5 of `l3_proof.md`: the fiber of the map `pair ↦ (its first Dyck rotation)` sits
inside one rotation orbit, which has at most `k^2` members. -/
theorem choose_sq_le_card_dyckSet {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k) :
    k.choose s ^ 2 ≤ k ^ 2 * (dyckSet k s).card := by
  classical
  have hk : 0 < k := lt_of_lt_of_le hs hsk
  have key : ∀ p : Finset ℕ × Finset ℕ, ∃ q : (ℕ × ℕ) × DyckWord,
      p ∈ parityPairs k s → q.2.toList = pairWord k p.1 p.2 q.1.1 q.1.2 := by
    intro p
    by_cases hp : p ∈ parityPairs k s
    · obtain ⟨hxp, hxcp, hyp, hycp⟩ := mem_parityPairs hp
      obtain ⟨a, b, w, hw⟩ := exists_rotation_isDyck hs hsk hxp hxcp hyp hycp
      exact ⟨((a, b), w), fun _ => hw⟩
    · exact ⟨((0, 0), 0), fun h => absurd h hp⟩
  choose F hF using key
  have hmaps : Set.MapsTo (fun p => ((F p).2, ((F p).1.1 % k, (F p).1.2 % k)))
      (parityPairs k s : Set (Finset ℕ × Finset ℕ))
      (((dyckSet k s) ×ˢ ((range k) ×ˢ (range k)) : Finset (DyckWord × ℕ × ℕ)) :
        Set (DyckWord × ℕ × ℕ)) := by
    intro p hp
    simp only [Finset.mem_coe] at hp ⊢
    obtain ⟨hxp, hxcp, hyp, hycp⟩ := mem_parityPairs hp
    have hwp := hF p hp
    have hlen : (F p).2.toList.length = 2 * k := by
      rw [hwp]; exact pairWord_length k p.1 p.2 _ _
    have hsemi : (F p).2.semilength = k := by
      have h2 := (F p).2.two_mul_semilength_eq_length
      omega
    have hrows : dyckRows (F p).2 = s := by
      rw [dyckRows_eq_card_even_ups, hsemi, hwp]
      exact pairWord_rows hxp hxcp _ _
    rw [Finset.mem_product]
    refine ⟨mem_dyckSet_of hsemi hrows, ?_⟩
    rw [Finset.mem_product]
    exact ⟨Finset.mem_range.mpr (Nat.mod_lt _ hk), Finset.mem_range.mpr (Nat.mod_lt _ hk)⟩
  have hinj : Set.InjOn (fun p => ((F p).2, ((F p).1.1 % k, (F p).1.2 % k)))
      (parityPairs k s : Set (Finset ℕ × Finset ℕ)) := by
    intro p hp q hq hpq
    simp only [Finset.mem_coe] at hp hq
    obtain ⟨hxp, hxcp, hyp, hycp⟩ := mem_parityPairs hp
    obtain ⟨hxq, hxcq, hyq, hycq⟩ := mem_parityPairs hq
    have hwp := hF p hp
    have hwq := hF q hq
    simp only [Prod.mk.injEq] at hpq
    obtain ⟨hw, ha, hb⟩ := hpq
    obtain ⟨rp1, rp2⟩ := recover_pairWord hxp hyp _ _ hwp
    obtain ⟨rq1, rq2⟩ := recover_pairWord hxq hyq _ _ hwq
    have h1 : p.1 = q.1 := by
      rw [← rp1, ← rq1, hw, rot_fun_mod k _ _ ha]
    have h2 : p.2 = q.2 := by
      rw [← rp2, ← rq2, hw, rot_fun_mod k _ _ hb]
    exact Prod.ext h1 h2
  have hcard := Finset.card_le_card_of_injOn _ hmaps hinj
  rw [card_parityPairs hsk, Finset.card_product, Finset.card_product, Finset.card_range] at hcard
  calc k.choose s ^ 2 ≤ (dyckSet k s).card * (k * k) := hcard
    _ = k ^ 2 * (dyckSet k s).card := by ring

/-- Lemma 6 of `l3_proof.md`: `C(k,s) * s = C(k,s-1) * (k-s+1)` and `1 ≤ k - s + 1`. -/
theorem choose_pred_le_mul {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k) :
    k.choose (s - 1) ≤ k * k.choose s := by
  obtain ⟨r, rfl⟩ : ∃ r, s = r + 1 := ⟨s - 1, by omega⟩
  have h := Nat.choose_succ_right_eq k r
  have h1 : 1 ≤ k - r := by omega
  have hr : r + 1 - 1 = r := by omega
  rw [hr]
  calc k.choose r = k.choose r * 1 := (mul_one _).symm
    _ ≤ k.choose r * (k - r) := Nat.mul_le_mul_left _ h1
    _ = k.choose (r + 1) * (r + 1) := h.symm
    _ ≤ k.choose (r + 1) * k := Nat.mul_le_mul_left _ (by omega)
    _ = k * k.choose (r + 1) := Nat.mul_comm _ _

/-- L3, Theorem 7 of `l3_proof.md`: the Narayana lower bound with the `k^2` slack the plan
allows. The exact count is `C(k,s) C(k,s-1) = k * (dyckSet k s).card`. -/
theorem narayana_lower {k s : ℕ} (hs : 1 ≤ s) (hsk : s ≤ k) :
    k.choose s * k.choose (s - 1) ≤ k ^ 3 * (dyckSet k s).card := by
  calc k.choose s * k.choose (s - 1) ≤ k.choose s * (k * k.choose s) :=
        Nat.mul_le_mul_left _ (choose_pred_le_mul hs hsk)
    _ = k * k.choose s ^ 2 := by ring
    _ ≤ k * (k ^ 2 * (dyckSet k s).card) :=
        Nat.mul_le_mul_left _ (choose_sq_le_card_dyckSet hs hsk)
    _ = k ^ 3 * (dyckSet k s).card := by ring

/-- L4, `tree_cell_lower` of `X8_skeleton.lean`: L3 and `dyck_label_card_le` of L2. -/
theorem tree_cell_lower {n d k s t : ℕ} (hk : 2 ≤ k) (hs : 1 ≤ s) (ht : 1 ≤ t)
    (hsum : s + t = k + 1) (hnd : 2 * k ≤ min n d) :
    k.choose s * k.choose (s - 1) * (n.descFactorial s * d.descFactorial t)
      ≤ k ^ 3 * cellCard n d k s t := by
  have hmn : min n d ≤ n := min_le_left n d
  have hmd : min n d ≤ d := min_le_right n d
  have hn : k ≤ n := by omega
  have hd : k ≤ d := by omega
  have hsk : s ≤ k := by omega
  calc k.choose s * k.choose (s - 1) * (n.descFactorial s * d.descFactorial t)
      ≤ k ^ 3 * (dyckSet k s).card * (n.descFactorial s * d.descFactorial t) :=
        Nat.mul_le_mul_right _ (narayana_lower hs hsk)
    _ = k ^ 3 * ((dyckSet k s).card * (n.descFactorial s * d.descFactorial t)) := by ring
    _ ≤ k ^ 3 * cellCard n d k s t :=
        Nat.mul_le_mul_left _ (dyck_label_card_le (by omega) hsum hn hd)

end Edge
end StackedSVD
