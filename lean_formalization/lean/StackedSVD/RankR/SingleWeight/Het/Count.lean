/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Frame

/-!
# Unit G2c: the eigenvalue count at an intermediate gap point

Track G, unit G2c of `notes/archive/trackG_plan.md` (2026-09-05), risk 3. Namespace
`StackedSVD.Frame`, next to `Frame.count_of_split_of_frame` (`LinAlg/Frame.lean:977`), the
file this one mirrors. Nothing here is random: every declaration is a statement about one
real symmetric matrix, one threshold, and one finite family of approximate eigenvectors.

## Why a new count

`Frame.count_of_split_of_frame` counts the eigenvalues above a threshold `τ` that sits below
**all** `r` outliers of `S = W + Q Qᵀ`. `OutliersR.count_eq_of_frame_of_edge`
(`RankR/RMT/OutliersG.lean:507`) allows a partial frame of `s < r` vectors, but it takes the
edge bound `∀ k, s ≤ k → eigenvalues₀ k ≤ τ` ready made, and Track E supplies that from a
sub-block edge over a column subset of `Q`. At general `r_i` the outliers are not attached to
columns of `Q`, so there is no column subset to take and that route is closed.

`count_eq_of_frame_of_split` replaces it. It counts at an intermediate `τ`: `s` of the `r`
approximate eigenvalues `ρ k` sit above `τ + mg`, the other `r - s` sit at `ρ k + 2 δ ≤ τ`,
and every `ρ k` sits above `lamMax W + 2 δ`. The lower bound `s ≤ #{i | τ < λ_i}` is
`Frame.le_card_filter_of_frame` on the first `s` vectors. The upper bound is the new part:
the `r - s` low vectors localize `r - s` **distinct** eigenvalues inside `(lamMax W, τ]`
(`exists_eigenvalue_near_two` for the localization, the separation `hsep` for the
distinctness), while `Frame.card_filter_le_of_split` caps the whole window `(lamMax W, ∞)` at
`r` eigenvalues. So at most `s` eigenvalues are left above `τ`.

The near radius is `2 δ`, not `δ`, because the frame is not exactly unit norm. The consumer
(unit G3) builds its frame from `OutliersR.yhatv`, which is only Gram close to orthonormal,
so `‖y k‖ = 1` is not available. Under `δ' ≤ 1 / 2` the diagonal of `hgram` gives
`1 / 2 ≤ ‖y k‖ ^ 2`, and `exists_eigenvalue_near_two` then localizes within `2 δ`.

## Content

0. `exists_eigenvalue_near_two`: `Frame.exists_eigenvalue_near` (`LinAlg/Frame.lean:401`)
   without exact unit norm. A vector with `1 / 2 ≤ ‖y‖ ^ 2` and residual at most `δ` at `ρ`
   puts an eigenvalue of `S` within `2 δ` of `ρ`.
1. `count_eq_of_frame_of_split`: the count in the matrix index,
   `#{i | τ < hS.eigenvalues i} = s`.
2. `count_of_split_of_frame_gap`: the same count in the sorted index,
   `∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < s`. This is the form the consumers state
   (compare `Frame.count_of_split_of_frame`). The sorted set is a down set because
   `eigenvalues₀` is antitone, so the cardinality fixes it.

## Hypotheses

`hsr : s ≤ r` is not in the plan and is not implied by the rest. Without it the statement is
false: take `r = 0`, `s = 1`, `p = 1`, `W = 0`, `τ = 0`. Every hypothesis quantified over
`Fin r` is then vacuous, `hsmall` holds at `δ = δ' = 0` and `mg = 1`, and the claim reads
`#{i | 0 < λ_i(W)} = 1` against `lamMax W ≤ 0`.

`hδ'h : δ' ≤ 1 / 2` replaces the earlier `hyn : ∀ k, ‖y k‖ = 1`. It is the only thing the
localization needs about the norms, and every consumer that controls the Gram already has it.
No rescaling of the frame is required.

`hedge : ∀ k, lamMax W hW + 2 δ < ρ k` is not in the plan. It is what puts every localized
eigenvalue inside the window `(lamMax W, ∞)` that `card_filter_le_of_split` counts. It is
necessary, not only convenient: with `p = 40`, `r = 2`, `s = 1`, two strong columns in `Q`
and `y 1` a bulk eigenvector of `S` at `ρ 1 = 1.80 < lamMax W = 2.386`, every other
hypothesis holds at `δ = 1e-10`, `τ = lamMax W + 1e-6`, and the count is `2`, not `s = 1`
(seed 20260905).

`hsep` is weaker than the plan's `∀ k l, k ≠ l → 2 δ < |ρ k - ρ l|` in the index range and
stronger in the radius: the proof separates only the **low** points, so ties among the `s`
detected outliers are allowed, as in `Frame.count_of_split_of_frame`, but the radius is
`4 δ` because each low point is localized within `2 δ`. Whether the low separation itself
can be dropped is open; a two sided window frame argument would replace it, and that is a
new spectral lemma.

`hτ : lamMax W hW ≤ τ` is not implied by `hedge` and `hbelow` when `s = r`, so it stays.

Numeric satisfiability check (seed 20260905, `p = 40`, `n = 120`, `r = 4`): four instances
(`s = 2` with four detectable outliers, `s = 2` with a tied high pair, `s = r`, `s = 1` with
three low points) satisfy every hypothesis in the `2 δ` and `4 δ` form at `δ ≈ 6e-2`,
`δ' ≈ 0.2` and give the claimed count. The frames there are deliberately not unit norm:
`‖y k‖ ^ 2` runs from 0.83 to 1.20.
-/

open Matrix Finset
open scoped InnerProductSpace Matrix

namespace StackedSVD
namespace Frame

/-- `Frame.exists_eigenvalue_near` (`LinAlg/Frame.lean:401`) without exact unit norm: a vector
with `1 / 2 ≤ ‖y‖ ^ 2` and residual at most `δ` puts an eigenvalue of `S` within `2 * δ` of
`ρ`. The mass lost by the norm defect costs a factor `2` in the radius. -/
theorem exists_eigenvalue_near_two {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    {ρ δ : ℝ} {y : EuclideanSpace ℝ (Fin p)} (hy : 1 / 2 ≤ ‖y‖ ^ 2)
    (hres : ‖toOp S y - ρ • y‖ ≤ δ) : ∃ i, |hS.eigenvalues i - ρ| ≤ 2 * δ := by
  by_contra hcon
  push Not at hcon
  have hδ : 0 ≤ δ := le_trans (norm_nonneg _) hres
  have hsq : ‖toOp S y - ρ • y‖ ^ 2 ≤ δ ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hres 2
  rw [norm_sq_residual hS] at hsq
  have hone : ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 = ‖y‖ ^ 2 := (norm_sq_eq_sum_coord hS y).symm
  have hhalf : 1 / 2 ≤ ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 := by
    rw [hone]
    exact hy
  obtain ⟨i₀, -, hi₀⟩ := Finset.exists_ne_zero_of_sum_ne_zero
    (show ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 ≠ 0 by
      intro hz
      rw [hz] at hhalf
      linarith)
  have hlt : 4 * δ ^ 2 * ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2
      < ∑ i, (hS.eigenvalues i - ρ) ^ 2 * ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 := by
    rw [Finset.mul_sum]
    refine Finset.sum_lt_sum (fun i _ => ?_) ⟨i₀, Finset.mem_univ _, ?_⟩
    · have h2 : (2 * δ) ^ 2 ≤ |hS.eigenvalues i - ρ| ^ 2 :=
        pow_le_pow_left₀ (by linarith) (hcon i).le 2
      rw [sq_abs] at h2
      have h3 : 4 * δ ^ 2 ≤ (hS.eigenvalues i - ρ) ^ 2 := by nlinarith [h2]
      exact mul_le_mul_of_nonneg_right h3 (sq_nonneg _)
    · have h2 : (2 * δ) ^ 2 < |hS.eigenvalues i₀ - ρ| ^ 2 :=
        pow_lt_pow_left₀ (hcon i₀) (by linarith) (by norm_num)
      rw [sq_abs] at h2
      have h3 : 4 * δ ^ 2 < (hS.eigenvalues i₀ - ρ) ^ 2 := by nlinarith [h2]
      exact mul_lt_mul_of_pos_right h3 ((sq_nonneg _).lt_of_ne hi₀.symm)
  nlinarith [hsq, hlt, hhalf, sq_nonneg δ,
    mul_nonneg (mul_nonneg (by norm_num : (0 : ℝ) ≤ 4) (sq_nonneg δ))
      (by linarith : (0 : ℝ) ≤ (∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2) - 1 / 2)]

/-- The count at an intermediate gap point. `S = W + Q Qᵀ` with `Q` of `r` columns, `τ` above
`lamMax W`, and `r` approximate eigenvectors `y k` with residual `δ` at the points `ρ k`. The
first `s` points sit at least `mg` above `τ`, the last `r - s` sit at least `2 δ` below `τ`,
every point sits more than `2 δ` above `lamMax W`, and distinct low points are more than
`4 δ` apart. The frame is not assumed unit norm; `δ' ≤ 1 / 2` replaces that. Then exactly `s`
eigenvalues of `S` exceed `τ`. -/
theorem count_eq_of_frame_of_split {p r s : ℕ} {S W : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) (hsr : s ≤ r) {τ mg δ δ' : ℝ}
    (hτ : lamMax W hW ≤ τ) (hmg : 0 < mg) (hδ' : 0 ≤ δ') (hδ'h : δ' ≤ 1 / 2)
    {y : Fin r → EuclideanSpace ℝ (Fin p)} {ρ : Fin r → ℝ}
    (hedge : ∀ k, lamMax W hW + 2 * δ < ρ k)
    (hbelow : ∀ k : Fin r, s ≤ (k : ℕ) → ρ k + 2 * δ ≤ τ)
    (habove : ∀ k : Fin r, (k : ℕ) < s → τ + mg ≤ ρ k)
    (hsep : ∀ k l : Fin r, s ≤ (k : ℕ) → s ≤ (l : ℕ) → k ≠ l → 4 * δ < |ρ k - ρ l|)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : (s : ℝ) * (δ' + (δ / mg) ^ 2) < 1) :
    (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s := by
  classical
  -- the lower bound: the `s` high vectors are a frame above `τ + mg`
  have hlow : s ≤ (Finset.univ.filter fun i => τ < hS.eigenvalues i).card := by
    have hres' : ∀ k : Fin s,
        ∑ i, (hS.eigenvalues i - ρ (Fin.castLE hsr k)) ^ 2
          * (fun k : Fin s => fun i => ⟪hS.eigenvectorBasis i, y (Fin.castLE hsr k)⟫_ℝ)
            k i ^ 2 ≤ δ ^ 2 := by
      intro k
      have h := pow_le_pow_left₀ (norm_nonneg _) (hres (Fin.castLE hsr k)) 2
      rwa [norm_sq_residual hS] at h
    have hgram' : ∀ k l : Fin s,
        |(fun k : Fin s => fun i => ⟪hS.eigenvectorBasis i, y (Fin.castLE hsr k)⟫_ℝ) k
            ⬝ᵥ (fun k : Fin s => fun i => ⟪hS.eigenvectorBasis i, y (Fin.castLE hsr k)⟫_ℝ) l
          - if k = l then 1 else 0| ≤ δ' := by
      intro k l
      have h : (fun k : Fin s => fun i => ⟪hS.eigenvectorBasis i, y (Fin.castLE hsr k)⟫_ℝ) k
          ⬝ᵥ (fun k : Fin s => fun i => ⟪hS.eigenvectorBasis i, y (Fin.castLE hsr k)⟫_ℝ) l
          = ⟪y (Fin.castLE hsr k), y (Fin.castLE hsr l)⟫_ℝ := by
        rw [inner_eq_sum_coord hS]
        rfl
      rw [h]
      have h2 := hgram (Fin.castLE hsr k) (Fin.castLE hsr l)
      rcases eq_or_ne k l with rfl | hne
      · simpa using h2
      · rw [if_neg hne]
        rw [if_neg fun hc => hne (Fin.castLE_injective hsr hc)] at h2
        exact h2
    have hρ' : ∀ k : Fin s, τ + mg ≤ ρ (Fin.castLE hsr k) := by
      intro k
      refine habove (Fin.castLE hsr k) ?_
      rw [Fin.val_castLE]
      exact k.isLt
    exact le_card_filter_of_frame hmg hδ' hρ' hres' hgram' hsmall
  -- the upper bound: the `r - s` low vectors use up `r - s` of the `r` slots above `lamMax W`
  have hup : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card ≤ s := by
    have hAB : (Finset.univ.filter fun i => τ < hS.eigenvalues i)
        ⊆ Finset.univ.filter fun i => lamMax W hW < hS.eigenvalues i := by
      intro i hi
      rw [Finset.mem_filter] at hi ⊢
      exact ⟨hi.1, lt_of_le_of_lt hτ hi.2⟩
    have hBr : (Finset.univ.filter fun i => lamMax W hW < hS.eigenvalues i).card ≤ r :=
      card_filter_le_of_split hW hS hSeq le_rfl
    have hyhalf : ∀ k : Fin r, 1 / 2 ≤ ‖y k‖ ^ 2 := by
      intro k
      have hk := hgram k k
      rw [if_pos rfl, real_inner_self_eq_norm_sq] at hk
      have hk2 := (abs_le.mp hk).1
      linarith
    have hnear : ∀ k : Fin r, ∃ i, |hS.eigenvalues i - ρ k| ≤ 2 * δ := fun k =>
      exists_eigenvalue_near_two hS (hyhalf k) (hres k)
    choose ι hι using hnear
    have hmaps : ∀ k ∈ Finset.univ.filter fun k : Fin r => ¬ ((k : ℕ) < s),
        ι k ∈ (Finset.univ.filter fun i => lamMax W hW < hS.eigenvalues i)
          \ (Finset.univ.filter fun i => τ < hS.eigenvalues i) := by
      intro k hk
      rw [Finset.mem_filter] at hk
      have hks : s ≤ (k : ℕ) := not_lt.mp hk.2
      obtain ⟨h1, h2⟩ := abs_le.mp (hι k)
      have he := hedge k
      have hb := hbelow k hks
      rw [Finset.mem_sdiff, Finset.mem_filter, Finset.mem_filter]
      refine ⟨⟨Finset.mem_univ _, by linarith⟩, ?_⟩
      rintro ⟨-, hlt⟩
      linarith
    have hinj : Set.InjOn ι
        ((Finset.univ.filter fun k : Fin r => ¬ ((k : ℕ) < s) : Finset (Fin r)) : Set (Fin r)) := by
      intro k hkC l hlC hkl
      rw [Finset.mem_coe] at hkC hlC
      by_contra hne
      obtain ⟨h1, h2⟩ := abs_le.mp (hι k)
      obtain ⟨h3, h4⟩ := abs_le.mp (hι l)
      rw [hkl] at h1 h2
      have habs : |ρ k - ρ l| ≤ 4 * δ := by
        rw [abs_le]
        constructor <;> linarith
      exact absurd (hsep k l (not_lt.mp (Finset.mem_filter.mp hkC).2)
        (not_lt.mp (Finset.mem_filter.mp hlC).2) hne) (not_lt.mpr habs)
    have hCcard := Finset.card_le_card_of_injOn ι hmaps hinj
    have hsplit : (Finset.univ.filter fun k : Fin r => (k : ℕ) < s).card
        + (Finset.univ.filter fun k : Fin r => ¬ ((k : ℕ) < s)).card
        = (Finset.univ : Finset (Fin r)).card :=
      Finset.card_filter_add_card_filter_not _
    rw [card_filter_lt hsr, Finset.card_univ, Fintype.card_fin] at hsplit
    have hBA : ((Finset.univ.filter fun i => lamMax W hW < hS.eigenvalues i)
        \ (Finset.univ.filter fun i => τ < hS.eigenvalues i)).card
        + (Finset.univ.filter fun i => τ < hS.eigenvalues i).card
        = (Finset.univ.filter fun i => lamMax W hW < hS.eigenvalues i).card :=
      Finset.card_sdiff_add_card_eq_card hAB
    omega
  exact le_antisymm hup hlow

/-- The sorted form of `count_eq_of_frame_of_split`, the shape the consumer takes (compare
`Frame.count_of_split_of_frame` and `OutliersR.count_eq_of_frame_of_edge`): the sorted
eigenvalue indices above `τ` are exactly the top `s`. -/
theorem count_of_split_of_frame_gap {p r s : ℕ} {S W : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) (hsr : s ≤ r) {τ mg δ δ' : ℝ}
    (hτ : lamMax W hW ≤ τ) (hmg : 0 < mg) (hδ' : 0 ≤ δ') (hδ'h : δ' ≤ 1 / 2)
    {y : Fin r → EuclideanSpace ℝ (Fin p)} {ρ : Fin r → ℝ}
    (hedge : ∀ k, lamMax W hW + 2 * δ < ρ k)
    (hbelow : ∀ k : Fin r, s ≤ (k : ℕ) → ρ k + 2 * δ ≤ τ)
    (habove : ∀ k : Fin r, (k : ℕ) < s → τ + mg ≤ ρ k)
    (hsep : ∀ k l : Fin r, s ≤ (k : ℕ) → s ≤ (l : ℕ) → k ≠ l → 4 * δ < |ρ k - ρ l|)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : (s : ℝ) * (δ' + (δ / mg) ^ 2) < 1) :
    ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < s := by
  classical
  have hcard := count_eq_of_frame_of_split hW hS hSeq hsr hτ hmg hδ' hδ'h hedge hbelow habove
    hsep hres hgram hsmall
  have hsp : s ≤ p := by
    rw [← hcard]
    calc (Finset.univ.filter fun i => τ < hS.eigenvalues i).card
        ≤ (Finset.univ : Finset (Fin p)).card := Finset.card_filter_le _ _
      _ = p := by simp
  have hcard0 : (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k).card = s := by
    rw [← card_filter_eigenvalues hS (fun t => τ < t)]
    exact hcard
  -- the sorted set is a down set, so the cardinality fixes it
  have hsub : (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k)
      ⊆ Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < s := by
    intro k hk
    rw [Finset.mem_filter] at hk ⊢
    refine ⟨hk.1, ?_⟩
    have hdown : (Finset.univ.filter
        fun l : Fin (Fintype.card (Fin p)) => (l : ℕ) < (k : ℕ) + 1)
        ⊆ Finset.univ.filter fun l => τ < hS.eigenvalues₀ l := by
      intro l hl
      rw [Finset.mem_filter] at hl ⊢
      refine ⟨hl.1, lt_of_lt_of_le hk.2 (hS.eigenvalues₀_antitone ?_)⟩
      rw [Fin.le_def]
      omega
    have h1 := Finset.card_le_card hdown
    rw [card_filter_lt (by omega : (k : ℕ) + 1 ≤ Fintype.card (Fin p)), hcard0] at h1
    omega
  have hcards : (Finset.univ.filter
      fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < s).card = s :=
    card_filter_lt (by simpa using hsp)
  have heq := Finset.eq_of_subset_of_card_le hsub (le_of_eq (by rw [hcards, hcard0]))
  intro k
  constructor
  · intro h
    have hmem : k ∈ Finset.univ.filter fun k => τ < hS.eigenvalues₀ k :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, h⟩
    rw [heq] at hmem
    exact (Finset.mem_filter.mp hmem).2
  · intro h
    have hmem : k ∈ Finset.univ.filter
        fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < s :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, h⟩
    rw [← heq] at hmem
    exact (Finset.mem_filter.mp hmem).2

end Frame
end StackedSVD
