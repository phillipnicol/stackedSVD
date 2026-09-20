/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Frame
import StackedSVD.LinAlg.SpecIdx

/-!
# The window projector at one eigenvalue index

Task C4 of `notes/archive/rankr_plan_C.md` section 2, plus the down-set lemma that the audit
`notes/archive/audit_rankr_plan_C_2026-09-02.md` (attack 2) asks for.

Section 7 at general `r_i` reads the overlap at **one** sorted eigenvalue index `j`
(`overlapIdx`, `LinAlg/SpecIdx.lean`). Every limit law of Track A and Track C is stated for a
projector on a half line `Set.Ioi τ`. This file joins the two: two thresholds `a < b` that
trap the sorted eigenvalue `λ_j` turn the index projector into the window projector on
`Set.Ioc a b`, and Pythagoras turns the window norm into a difference of two half-line norms.
The rank-1 mirror is `Frame.specProjTop_eq_specProj_Ioi` (`LinAlg/Frame.lean:611`) at `r = 1`.

The window is half open. With `j` eigenvalues above `b` and `j + 1` above `a`, the sorted
eigenvalue `λ_j` sits in `(a, b]`, not in `(a, b)`.

## Content

1. The down-set lemma, `Frame.lt_eigenvalues₀_iff_of_card`. A count of the sorted indices
   above `τ` gives the `↔` shape `τ < λ_k ↔ k < s`, because `eigenvalues₀` is antitone, so
   `{k | τ < λ_k}` is a down-set of `Fin (Fintype.card (Fin p))`, and a down-set of
   cardinality `s` is `{k | k < s}`. `Frame.card_filter_of_count` (`Frame.lean:598`) is the
   converse. `Frame.lt_eigenvalues₀_iff_of_card_eigenvalues` takes the same count on the
   unsorted matrix index, through `Frame.card_filter_eigenvalues` (`Frame.lean:542`).
2. C4.1, `normSq_specProj_Ioc`: Pythagoras on `Set.Ioi a = Set.Ioc a b ∪ Set.Ioi b`.
3. C4.2, `specProjIdx_eq_specProj_Ioc`: the two counts make the eigenvalue sets agree, so
   `Frame.specProj_congr_of_iff` identifies the two projectors. No range condition on `j` is
   needed: out of range both sets hold no eigenvalue.
4. C4.3, `overlapIdx_eq_normSq_specProj_Ioc`, and the sandwich
   `overlapIdx_eq_normSq_specProj_Ioi_sub` that task C5 consumes (steps 5 and 6 of its
   route).

## Route

`Frame.norm_sq_specProj_union` (`Frame.lean:369`) is Pythagoras for two disjoint eigenvalue
sets. `Frame.specProj_congr_of_iff` (`Frame.lean:348`) needs only that the two sets hold the
same eigenvalues, so a repeated eigenvalue costs nothing: membership depends on the value,
not on the index.
-/

open Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. The down-set lemma -/

namespace Frame

/-- A count of the sorted eigenvalue indices above `τ` gives the `↔` shape. The set
`{k | τ < λ_k}` is a down-set because `Matrix.IsHermitian.eigenvalues₀` is antitone, and a
down-set of `Fin (Fintype.card (Fin p))` with cardinality `s` is `{k | (k : ℕ) < s}`.
Converse of `Frame.card_filter_of_count`. -/
theorem lt_eigenvalues₀_iff_of_card {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}
    (hS : S.IsHermitian) {τ : ℝ} {s : ℕ}
    (hcard : (Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) =>
      τ < hS.eigenvalues₀ k).card = s) (k : Fin (Fintype.card (Fin p))) :
    τ < hS.eigenvalues₀ k ↔ (k : ℕ) < s := by
  have hdown : ∀ l m : Fin (Fintype.card (Fin p)), l ≤ m → τ < hS.eigenvalues₀ m →
      τ < hS.eigenvalues₀ l := fun l m hlm h =>
    lt_of_lt_of_le h (hS.eigenvalues₀_antitone hlm)
  constructor
  · intro h
    have hsub : (Finset.univ.filter fun l : Fin (Fintype.card (Fin p)) =>
        (l : ℕ) < (k : ℕ) + 1)
        ⊆ Finset.univ.filter fun l => τ < hS.eigenvalues₀ l := by
      intro l hl
      rw [Finset.mem_filter] at hl ⊢
      exact ⟨hl.1, hdown l k (Fin.le_def.mpr (by omega)) h⟩
    have hle := Finset.card_le_card hsub
    rw [card_filter_lt k.isLt, hcard] at hle
    omega
  · intro h
    by_contra hnk
    rw [not_lt] at hnk
    have hsub : (Finset.univ.filter fun l : Fin (Fintype.card (Fin p)) =>
        τ < hS.eigenvalues₀ l)
        ⊆ Finset.univ.filter fun l => (l : ℕ) < (k : ℕ) := by
      intro l hl
      rw [Finset.mem_filter] at hl ⊢
      refine ⟨hl.1, ?_⟩
      by_contra hge
      rw [not_lt] at hge
      exact absurd (hdown k l (Fin.le_def.mpr hge) hl.2) (not_lt.mpr hnk)
    have hle := Finset.card_le_card hsub
    rw [hcard, card_filter_lt (le_of_lt k.isLt)] at hle
    omega

/-- The same statement from the count on the unsorted matrix index, which is the shape that
`OutliersR.count_eq_of_forms_of_edge` delivers. -/
theorem lt_eigenvalues₀_iff_of_card_eigenvalues {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}
    (hS : S.IsHermitian) {τ : ℝ} {s : ℕ}
    (hcard : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    (k : Fin (Fintype.card (Fin p))) :
    τ < hS.eigenvalues₀ k ↔ (k : ℕ) < s := by
  refine lt_eigenvalues₀_iff_of_card hS ?_ k
  rw [← card_filter_eigenvalues hS (fun t => τ < t)]
  exact hcard

end Frame

/-! ### 2. C4.1: the window norm is a difference of two half-line norms -/

/-- Pythagoras on `Set.Ioi a = Set.Ioc a b ∪ Set.Ioi b`. Every limit law of Track A and
Track C is stated on a half line, so this is how a window enters. -/
theorem normSq_specProj_Ioc {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    {a b : ℝ} (hab : a ≤ b) (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProj S (Set.Ioc a b) x‖ ^ 2
      = ‖specProj S (Set.Ioi a) x‖ ^ 2 - ‖specProj S (Set.Ioi b) x‖ ^ 2 := by
  have hdisj : Disjoint (Set.Ioc a b) (Set.Ioi b) := by
    rw [Set.disjoint_left]
    rintro t ⟨-, ht2⟩ ht3
    exact absurd ht3 (not_lt.mpr ht2)
  have hunion : Set.Ioc a b ∪ Set.Ioi b = Set.Ioi a := by
    ext t
    simp only [Set.mem_union, Set.mem_Ioc, Set.mem_Ioi]
    constructor
    · rintro (⟨h1, -⟩ | h1)
      · exact h1
      · exact lt_of_le_of_lt hab h1
    · intro h1
      rcases le_or_gt t b with h2 | h2
      · exact Or.inl ⟨h1, h2⟩
      · exact Or.inr h2
  have h := Frame.norm_sq_specProj_union hS hdisj x
  rw [hunion] at h
  linarith

/-! ### 3. C4.2: two counts identify the index projector with the window projector -/

/-- With `j + 1` sorted eigenvalues above `a` and `j` above `b`, the projector at the sorted
index `j` is the projector on the half-open window `(a, b]`. The window is half open because
the count at `b` is a strict inequality, so `λ_j ≤ b`. -/
theorem specProjIdx_eq_specProj_Ioc {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}
    (hS : S.IsHermitian) {a b : ℝ} {j : ℕ}
    (hlow : ∀ k : Fin (Fintype.card (Fin p)), a < hS.eigenvalues₀ k ↔ (k : ℕ) < j + 1)
    (hhigh : ∀ k : Fin (Fintype.card (Fin p)), b < hS.eigenvalues₀ k ↔ (k : ℕ) < j) :
    specProjIdx S hS j = specProj S (Set.Ioc a b) := by
  rw [specProjIdx]
  refine Frame.specProj_congr_of_iff hS fun i => ?_
  have hi : hS.eigenvalues i = hS.eigenvalues₀ ((eigIdx p).symm i) := by
    rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
  rw [hi]
  constructor
  · rintro ⟨q, hq, hqk⟩
    have h1 : a < hS.eigenvalues₀ q := (hlow q).mpr (by omega)
    have h2 : hS.eigenvalues₀ q ≤ b := by
      by_contra hb
      exact absurd ((hhigh q).mp (not_le.mp hb)) (by omega)
    rw [← hqk]
    exact ⟨h1, h2⟩
  · rintro ⟨h1, h2⟩
    have hk1 : ((eigIdx p).symm i : ℕ) < j + 1 := (hlow _).mp h1
    have hk2 : ¬ (((eigIdx p).symm i : ℕ) < j) := fun h =>
      absurd ((hhigh _).mpr h) (not_lt.mpr h2)
    exact ⟨(eigIdx p).symm i, by omega, rfl⟩

/-! ### 4. C4.3: the same at the Gram matrix, and the sandwich -/

/-- The overlap at the sorted index `j` is the squared norm of the window projector. This is
`specProjIdx_eq_specProj_Ioc` at `S = Xᵀ * X`, which is what `overlapIdx` reads. -/
theorem overlapIdx_eq_normSq_specProj_Ioc {q p : ℕ} (X : Matrix (Fin q) (Fin p) ℝ)
    {a b : ℝ} {j : ℕ}
    (hlow : ∀ k : Fin (Fintype.card (Fin p)),
      a < (isHermitian_transpose_mul_self X).eigenvalues₀ k ↔ (k : ℕ) < j + 1)
    (hhigh : ∀ k : Fin (Fintype.card (Fin p)),
      b < (isHermitian_transpose_mul_self X).eigenvalues₀ k ↔ (k : ℕ) < j)
    (w : EuclideanSpace ℝ (Fin p)) :
    overlapIdx X j w = ‖specProj (Xᵀ * X) (Set.Ioc a b) w‖ ^ 2 := by
  rw [overlapIdx, specProjIdx_eq_specProj_Ioc (isHermitian_transpose_mul_self X) hlow hhigh]

/-- The sandwich task C5 consumes: the overlap at the sorted index `j` is the difference of
the two half-line overlaps, at the thresholds that trap `λ_j`. Steps 5 and 6 of the C5
route in one statement. -/
theorem overlapIdx_eq_normSq_specProj_Ioi_sub {q p : ℕ} (X : Matrix (Fin q) (Fin p) ℝ)
    {a b : ℝ} (hab : a ≤ b) {j : ℕ}
    (hlow : ∀ k : Fin (Fintype.card (Fin p)),
      a < (isHermitian_transpose_mul_self X).eigenvalues₀ k ↔ (k : ℕ) < j + 1)
    (hhigh : ∀ k : Fin (Fintype.card (Fin p)),
      b < (isHermitian_transpose_mul_self X).eigenvalues₀ k ↔ (k : ℕ) < j)
    (w : EuclideanSpace ℝ (Fin p)) :
    overlapIdx X j w
      = ‖specProj (Xᵀ * X) (Set.Ioi a) w‖ ^ 2 - ‖specProj (Xᵀ * X) (Set.Ioi b) w‖ ^ 2 := by
  rw [overlapIdx_eq_normSq_specProj_Ioc X hlow hhigh w,
    normSq_specProj_Ioc (isHermitian_transpose_mul_self X) hab w]

end StackedSVD
