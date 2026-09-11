/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.SpecProjPerturb

/-!
# General eigenvalue, rank and frame lemmas

STATUS 2026-08-30: proved, 0 `sorry`. Cleanup wave 2 moved every declaration here from
sections 0, 2 and 5b of `StackedSVD/RankR/Weighted.lean`, where the header already asked for
the move. No statement changed; only the home file did. `RankR/Weighted.lean` imports this
file, and so does anything else that needs the sorted spectrum of an inverse, of `A + I`, or
the top-`r` eigenframe.

Mathlib v4.33.0 has none of these on this pin. The three groups are:

1. **Sorted eigenvalues.** `eigenvalues₀_eq_of_charpoly` reads the sorted list off the
   characteristic polynomial and carries every other fact of the group:
   `eigenvalues₀_eq_of_charpoly_eq` (equal charpoly gives equal spectrum, which with
   `Matrix.charpoly_mul_comm` is the `B Bᵀ` versus `Bᵀ B` rule), `eigenvalues₀_inv` (the
   spectrum of `A⁻¹` is the reciprocal set in reverse order, `main_paper.tex:2120`) and
   `eigenvalues₀_add_one` (`A + I` shifts every sorted eigenvalue by one).
2. **Rank and the spectrum.** `card_nz_eigenvalues₀`, `eigenvalues₀_eq_zero_of_rank_le` and
   `eigenvalues₀_pos_of_lt_rank`: a positive semidefinite matrix of rank `ρ` has `λ_k > 0`
   exactly for `k < ρ`. `rank_sum_vecMulVec` is the companion rank bridge
   `Rank(∑_i R_i R_iᵀ) = Rank(R_stack)`, with `Rstack` and `Rstack_transpose_mul`.
3. **The top-`r` eigenframe and its trace forms.** `exists_topFrame` builds the paper's
   `Q_r`, `Λ_r` with `Qᵀ Q = I`, `Qᵀ A Q = Λ_r` and `specInvTop A r = Q Λ_r⁻¹ Qᵀ`;
   `trace_specInvTop_conj` reads `tr(Bᵀ (specInvTop A r) B)` in the eigenbasis; and
   `trace_specInvTop_add_one` evaluates it at a matrix of the shape `C Cᵀ + I`, where every
   sorted eigenvalue at or above index `r` is exactly `1`, so **no eigengap is needed**.
4. **The spectrum of a Gram matrix with orthogonal columns** (2026-09-02, item 2.3 of
   `notes/archive/rankr_D4_plan.md`). `eigenvalues₀_mul_transpose_of_transpose_mul_diagonal`: when
   `Cᵀ C = diag S` with `S` antitone and nonnegative and `q ≤ p`, the sorted spectrum of
   `C Cᵀ` is `S` padded by zeros. `antitone_padSpec` and `prod_X_sub_C_padSpec` support it.

`eigIdx`, `eigenvalues_eigIdx`, `spectral_conj`, `charpoly_eq_prod_of_conj`,
`eigenvalues₀_congr_mat`, `inv_conj`, the three `Fin` index helpers and
`trace_conj_vecMulVec` and `sum_sq_transpose_mulVec` support them. The local copy
`mem_topEigSet_iff'` is gone (cleanup wave 3): `LinAlg/SpecProjPerturb.mem_topEigSet_iff` is
now public and this file uses it.
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. Sorted eigenvalues from the characteristic polynomial -/

section SpectrumHelpers

open Polynomial Unitary

variable {p : ℕ}

/-- The canonical index equiv between the sorted index type and the matrix index type. -/
noncomputable def eigIdx (p : ℕ) : Fin (Fintype.card (Fin p)) ≃ Fin p :=
  Fintype.equivOfCardEq (Fintype.card_fin _)

theorem eigenvalues_eigIdx {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin p))) : hA.eigenvalues (eigIdx p k) = hA.eigenvalues₀ k := by
  rw [Matrix.IsHermitian.eigenvalues, eigIdx, Equiv.symm_apply_apply]

/-- The spectral theorem in the plain matrix form `A = U diag(λ) U⋆`. -/
theorem spectral_conj {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) :
    A = (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) * Matrix.diagonal hA.eigenvalues *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) := by
  conv_lhs => rw [hA.spectral_theorem]
  simp

/-- `U diag(f) U⋆` has the characteristic polynomial of `diag f`. -/
theorem charpoly_eq_prod_of_conj {U B : Matrix (Fin p) (Fin p) ℝ} {f : Fin p → ℝ}
    (hU : star U * U = 1) (hB : B = U * Matrix.diagonal f * star U) :
    B.charpoly = ∏ i, (X - C (f i)) := by
  rw [hB, Matrix.charpoly_mul_comm, ← Matrix.mul_assoc, hU, Matrix.one_mul,
    Matrix.charpoly_diagonal]

/-- A Hermitian matrix whose characteristic polynomial is `∏ (X - g k)` with `g` antitone has
`g` for its sorted eigenvalue list. -/
theorem eigenvalues₀_eq_of_charpoly {B : Matrix (Fin p) (Fin p) ℝ} (hB : B.IsHermitian)
    {g : Fin (Fintype.card (Fin p)) → ℝ} (hg : Antitone g)
    (hchar : B.charpoly = ∏ k, (X - C (g k))) : hB.eigenvalues₀ = g := by
  have hprod : (∏ k, (X - C (g k))) ≠ 0 := by
    simp [Finset.prod_ne_zero_iff, Polynomial.X_sub_C_ne_zero]
  have hroots : B.charpoly.roots = Multiset.map g Finset.univ.val := by
    rw [hchar, Polynomial.roots_prod _ _ hprod]
    simp
  have hsorted : (B.charpoly.roots.map RCLike.re).sort (· ≥ ·) = List.ofFn g := by
    rw [hroots]
    simp_rw [Fin.univ_val_map, Multiset.map_coe, List.map_ofFn, Function.comp_def,
      Multiset.coe_sort]
    apply List.mergeSort_of_pairwise
    simp_rw [decide_eq_true_eq, ← List.sortedGE_iff_pairwise]
    exact hg.sortedGE_ofFn
  have h := hB.sort_roots_charpoly_eq_eigenvalues₀
  rw [hsorted] at h
  exact List.ofFn_inj.mp h.symm

/-- `eigenvalues₀` transports along a matrix equality. -/
theorem eigenvalues₀_congr_mat {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) : hA.eigenvalues₀ = hB.eigenvalues₀ := by
  subst h
  rfl

/-- Two Hermitian matrices with the same characteristic polynomial have the same sorted
eigenvalues. With `Matrix.charpoly_mul_comm` this is the `B Bᵀ` versus `Bᵀ B` rule that step 2
of the paper's proof needs. -/
theorem eigenvalues₀_eq_of_charpoly_eq {B C : Matrix (Fin p) (Fin p) ℝ} (hB : B.IsHermitian)
    (hC : C.IsHermitian) (h : B.charpoly = C.charpoly) : hB.eigenvalues₀ = hC.eigenvalues₀ := by
  rw [← List.ofFn_inj, ← hB.sort_roots_charpoly_eq_eigenvalues₀,
    ← hC.sort_roots_charpoly_eq_eigenvalues₀, h]

/-- `(U D U⋆)⁻¹ = U E U⋆` when `D E = 1` and `U` is unitary. -/
theorem inv_conj {U D E : Matrix (Fin p) (Fin p) ℝ} (hUU : star U * U = 1)
    (hUU' : U * star U = 1) (hDE : D * E = 1) : (U * D * star U)⁻¹ = U * E * star U := by
  refine Matrix.inv_eq_right_inv ?_
  calc U * D * star U * (U * E * star U)
      = U * (D * (star U * U) * E) * star U := by simp only [Matrix.mul_assoc]
    _ = 1 := by rw [hUU, Matrix.mul_one, hDE, Matrix.mul_one, hUU']

/-- The sorted eigenvalues of the inverse of a positive definite matrix are the reciprocals in
reverse order (`main_paper.tex:2120`). Mathlib has no such rule on this pin. -/
theorem eigenvalues₀_inv {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hApos : ∀ k, 0 < hA.eigenvalues₀ k) (hAi : (A⁻¹).IsHermitian) (k) :
    hAi.eigenvalues₀ k = (hA.eigenvalues₀ k.rev)⁻¹ := by
  have hUU : star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 := Unitary.coe_star_mul_self _
  have hUU' : (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 := Unitary.coe_mul_star_self _
  have hne : ∀ i, hA.eigenvalues i ≠ 0 := by
    intro i
    have h := hApos ((eigIdx p).symm i)
    rw [Matrix.IsHermitian.eigenvalues]
    exact ne_of_gt h
  have hD : Matrix.diagonal hA.eigenvalues *
      Matrix.diagonal (fun i => (hA.eigenvalues i)⁻¹) = 1 := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    funext i
    exact mul_inv_cancel₀ (hne i)
  have hinv : A⁻¹ = (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      Matrix.diagonal (fun i => (hA.eigenvalues i)⁻¹) *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) := by
    conv_lhs => rw [spectral_conj hA]
    exact inv_conj hUU hUU' hD
  have hchar : (A⁻¹).charpoly = ∏ i, (X - C (hA.eigenvalues i)⁻¹) :=
    charpoly_eq_prod_of_conj hUU hinv
  have hchar' : (A⁻¹).charpoly = ∏ k : Fin (Fintype.card (Fin p)),
      (X - C (hA.eigenvalues₀ k.rev)⁻¹) := by
    rw [hchar]
    refine (Fintype.prod_equiv ((Fin.revPerm).trans (eigIdx p)) _ _ ?_).symm
    intro j
    simp only [Equiv.trans_apply, Fin.revPerm_apply]
    rw [eigenvalues_eigIdx]
  have hanti : Antitone fun k : Fin (Fintype.card (Fin p)) => (hA.eigenvalues₀ k.rev)⁻¹ := by
    intro a b hab
    have hrev : b.rev ≤ a.rev := Fin.rev_le_rev.mpr hab
    have h1 : hA.eigenvalues₀ a.rev ≤ hA.eigenvalues₀ b.rev := hA.eigenvalues₀_antitone hrev
    exact inv_anti₀ (hApos a.rev) h1
  exact congrFun (eigenvalues₀_eq_of_charpoly hAi hanti hchar') k

/-- Adding the identity shifts every sorted eigenvalue by one. -/
theorem eigenvalues₀_add_one {P : Matrix (Fin p) (Fin p) ℝ} (hP : P.IsHermitian)
    (hP1 : (P + 1).IsHermitian) (k) : hP1.eigenvalues₀ k = hP.eigenvalues₀ k + 1 := by
  have hUU : star (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 := Unitary.coe_star_mul_self _
  have hUU' : (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      star (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 := Unitary.coe_mul_star_self _
  have hshift : P + 1 = (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      Matrix.diagonal (fun i => hP.eigenvalues i + 1) *
      star (hP.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) := by
    have hd : Matrix.diagonal (fun i => hP.eigenvalues i + 1)
        = Matrix.diagonal hP.eigenvalues + 1 := by
      rw [← Matrix.diagonal_one, ← Matrix.diagonal_add]
    rw [hd, Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hUU', ← spectral_conj hP]
  have hchar : (P + 1).charpoly = ∏ i, (X - C (hP.eigenvalues i + 1)) :=
    charpoly_eq_prod_of_conj hUU hshift
  have hchar' : (P + 1).charpoly = ∏ k : Fin (Fintype.card (Fin p)),
      (X - C (hP.eigenvalues₀ k + 1)) := by
    rw [hchar]
    refine (Fintype.prod_equiv (eigIdx p) _ _ ?_).symm
    intro j
    rw [eigenvalues_eigIdx]
  have hanti : Antitone fun k : Fin (Fintype.card (Fin p)) => hP.eigenvalues₀ k + 1 :=
    fun a b hab => by simpa using hP.eigenvalues₀_antitone hab
  exact congrFun (eigenvalues₀_eq_of_charpoly hP1 hanti hchar') k

end SpectrumHelpers

section FinIndex

/-- Membership in the image of `Fin.castLEEmb`. -/
theorem mem_map_castLEEmb {n s : ℕ} (h : s ≤ n) (k : Fin n) :
    k ∈ Finset.univ.map (Fin.castLEEmb h) ↔ (k : ℕ) < s := by
  constructor
  · intro hk
    obtain ⟨i, -, rfl⟩ := Finset.mem_map.mp hk
    simp
  · intro hk
    exact Finset.mem_map.mpr ⟨⟨(k : ℕ), hk⟩, Finset.mem_univ _, by ext; rfl⟩

theorem card_filter_lt {n s : ℕ} (h : s ≤ n) :
    (Finset.univ.filter (fun j : Fin n => (j : ℕ) < s)).card = s := by
  classical
  have hset : (Finset.univ.filter (fun j : Fin n => (j : ℕ) < s))
      = Finset.univ.map (Fin.castLEEmb h) := by
    ext j
    rw [Finset.mem_filter, mem_map_castLEEmb h]
    simp
  rw [hset, Finset.card_map, Finset.card_univ, Fintype.card_fin]

/-- A sum whose terms vanish above `s` is a sum over `Fin s`. -/
theorem sum_eq_sum_castLE {α : Type*} [AddCommMonoid α] {n s : ℕ} (h : s ≤ n) (f : Fin n → α)
    (hzero : ∀ k : Fin n, s ≤ (k : ℕ) → f k = 0) :
    ∑ k : Fin n, f k = ∑ k : Fin s, f (Fin.castLE h k) := by
  classical
  have h1 : ∑ k ∈ Finset.univ.map (Fin.castLEEmb h), f k = ∑ k : Fin n, f k := by
    refine Finset.sum_subset (Finset.subset_univ _) ?_
    intro x _ hx
    exact hzero x (by simpa [mem_map_castLEEmb h] using hx)
  rw [← h1, Finset.sum_map]
  rfl

end FinIndex

section RankEigen

variable {p : ℕ}

/-- The number of nonzero sorted eigenvalues is the rank. -/
theorem card_nz_eigenvalues₀ {P : Matrix (Fin p) (Fin p) ℝ} (hP : P.IsHermitian) :
    (Finset.univ.filter (fun k => hP.eigenvalues₀ k ≠ 0)).card = P.rank := by
  classical
  have hcard : Fintype.card {i : Fin p // hP.eigenvalues i ≠ 0}
      = Fintype.card {k : Fin (Fintype.card (Fin p)) // hP.eigenvalues₀ k ≠ 0} :=
    Fintype.card_congr (Equiv.subtypeEquiv (eigIdx p)
      (fun k => by rw [eigenvalues_eigIdx])).symm
  rw [hP.rank_eq_card_non_zero_eigs, hcard, Fintype.card_subtype]

/-- A positive semidefinite matrix has a zero eigenvalue at every sorted index at or above its
rank. -/
theorem eigenvalues₀_eq_zero_of_rank_le {P : Matrix (Fin p) (Fin p) ℝ} (hP : P.IsHermitian)
    (hnn : ∀ k, 0 ≤ hP.eigenvalues₀ k) {s : ℕ} (hrank : P.rank ≤ s)
    {l : Fin (Fintype.card (Fin p))} (hl : s ≤ (l : ℕ)) : hP.eigenvalues₀ l = 0 := by
  classical
  by_contra hne
  have hpos : 0 < hP.eigenvalues₀ l := lt_of_le_of_ne (hnn l) (Ne.symm hne)
  have hsub : Finset.univ.filter (fun j : Fin (Fintype.card (Fin p)) => (j : ℕ) < (l : ℕ) + 1)
      ⊆ Finset.univ.filter (fun k => hP.eigenvalues₀ k ≠ 0) := by
    intro j hj
    rw [Finset.mem_filter] at hj ⊢
    refine ⟨Finset.mem_univ _, ?_⟩
    have hjl : j ≤ l := Fin.le_def.mpr (by omega)
    have h2 := hP.eigenvalues₀_antitone hjl
    intro hz
    rw [hz] at h2
    linarith
  have hc1 := Finset.card_le_card hsub
  rw [card_filter_lt (by omega), card_nz_eigenvalues₀] at hc1
  omega

/-- A positive semidefinite matrix has a positive eigenvalue at every sorted index below its
rank. -/
theorem eigenvalues₀_pos_of_lt_rank {P : Matrix (Fin p) (Fin p) ℝ} (hP : P.IsHermitian)
    (hnn : ∀ k, 0 ≤ hP.eigenvalues₀ k) {s : ℕ} (hrank : s ≤ P.rank)
    {k : Fin (Fintype.card (Fin p))} (hk : (k : ℕ) < s) : 0 < hP.eigenvalues₀ k := by
  classical
  rcases lt_or_eq_of_le (hnn k) with h | h
  · exact h
  exfalso
  have hsub : Finset.univ.filter (fun j => hP.eigenvalues₀ j ≠ 0)
      ⊆ Finset.univ.filter (fun j : Fin (Fintype.card (Fin p)) => (j : ℕ) < (k : ℕ)) := by
    intro j hj
    rw [Finset.mem_filter] at hj ⊢
    refine ⟨Finset.mem_univ _, ?_⟩
    by_contra hjk
    have hjk' : (k : ℕ) ≤ (j : ℕ) := not_lt.mp hjk
    have h2 := hP.eigenvalues₀_antitone (Fin.le_def.mpr hjk')
    rw [← h] at h2
    exact hj.2 (le_antisymm h2 (hnn j))
  have hc1 := Finset.card_le_card hsub
  rw [card_nz_eigenvalues₀, card_filter_lt (le_of_lt k.isLt)] at hc1
  omega

end RankEigen

section TraceForm

variable {p q : ℕ}

/-- `tr(Bᵀ (u uᵀ) B) = ‖Bᵀ u‖²`. -/
theorem trace_conj_vecMulVec (u : Fin p → ℝ) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Bᵀ * Matrix.vecMulVec u u * B) = ∑ j, (Bᵀ *ᵥ u) j ^ 2 := by
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.vecMulVec_apply,
    Matrix.transpose_apply, Matrix.mulVec, dotProduct, sq]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_mul_sum, Finset.sum_comm]
  simp only [Finset.sum_mul]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- `‖Bᵀ u‖² = uᵀ (B Bᵀ) u`. -/
theorem sum_sq_transpose_mulVec (B : Matrix (Fin p) (Fin q) ℝ) (u : Fin p → ℝ) :
    ∑ j, (Bᵀ *ᵥ u) j ^ 2 = ((B * Bᵀ) *ᵥ u) ⬝ᵥ u := by
  simp only [Matrix.mulVec, dotProduct, Matrix.mul_apply, Matrix.transpose_apply, sq,
    Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun j _ => by ring

/-- `tr(Bᵀ (specInvTop A r) B)` read in the eigenbasis of `A`. -/
theorem trace_specInvTop_conj (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Bᵀ * specInvTop A hA r * B)
      = ∑ i : Fin p, Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues i) *
          ∑ j : Fin q, (Bᵀ *ᵥ (WithLp.ofLp (hA.eigenvectorBasis i))) j ^ 2 := by
  rw [specInvTop]
  simp only [Matrix.mul_sum, Matrix.sum_mul, Matrix.trace_sum, Matrix.mul_smul,
    Matrix.smul_mul, Matrix.trace_smul, smul_eq_mul]
  exact Finset.sum_congr rfl fun i _ => by rw [trace_conj_vecMulVec]

end TraceForm

/-! ### 2. `R_stack` and the rank bridge

`Rstack` is the vertical concatenation of the `R_iᵀ` (`main_paper.tex:797`). It carries no
model and no probability, so it lives here with the rank lemma that reads it. -/

section RankBridge

variable {M r : ℕ}

/-- `R_stack ∈ ℝ^{r̃ × r}`, the vertical concatenation of the `R_iᵀ` (`main_paper.tex:797`).
At `r_i = 1` row `i` is the unit vector `R_i`. -/
noncomputable def Rstack (R : Fin M → EuclideanSpace ℝ (Fin r)) : Matrix (Fin M) (Fin r) ℝ :=
  Matrix.of fun i k => R i k

/-- `R_stackᵀ R_stack = ∑_i R_i R_iᵀ`, the matrix of the paper's rank condition
`Rank(∑_i R_i R_iᵀ) = r` (`assum:unaligned`, `main_paper.tex:757`). Over `ℝ` that rank equals
`(Rstack R).rank`, which is the form the statements below use. -/
theorem Rstack_transpose_mul (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (Rstack R)ᵀ * Rstack R
      = ∑ i, Matrix.vecMulVec (WithLp.ofLp (R i)) (WithLp.ofLp (R i)) := by
  ext k l
  simp [Rstack, Matrix.mul_apply, Matrix.sum_apply, Matrix.vecMulVec_apply]

/-- The rank bridge: the paper's `Rank(∑_i R_i R_iᵀ) = r` (`assum:unaligned`,
`main_paper.tex:757`) and the Lean form `(Rstack R).rank = r` are the same condition, by
`Matrix.rank_transpose_mul_self`. -/
theorem rank_sum_vecMulVec (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (∑ i, Matrix.vecMulVec (WithLp.ofLp (R i)) (WithLp.ofLp (R i))).rank = (Rstack R).rank := by
  rw [← Rstack_transpose_mul, Matrix.rank_transpose_mul_self]

end RankBridge


/-! ### 3. The top-`r` eigenframe

Step 1 of the proof (`main_paper.tex:2029`) needs, for a weight matrix with a top-`r` gap and a
positive `λ_{r-1}`, an `r`-column frame `Q` of eigenvectors with `Qᵀ Q = I`,
`Qᵀ A_W Q = Λ_r` and `specInvTop A_W r = Q Λ_r^{-1} Qᵀ`. `exists_topFrame` builds it once.
-/

section Frame

variable {p : ℕ}

/-- The top-`r` eigenframe of a Hermitian matrix with a top-`r` gap and a positive
`λ_{r-1}`: `Q Λ_r Qᵀ` is the paper's `Q_r Λ_r Q_rᵀ` and `specInvTop A r = Q Λ_r^{-1} Qᵀ`. -/
theorem exists_topFrame {r : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p)) (hgap : TopGap A hA r)
    (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩) :
    ∃ (Q : Matrix (Fin p) (Fin r) ℝ) (L : Fin r → ℝ), (∀ k, 0 < L k) ∧ Qᵀ * Q = 1 ∧
      Qᵀ * A * Q = Matrix.diagonal L ∧
      specInvTop A hA r = Q * Matrix.diagonal (fun k => (L k)⁻¹) * Qᵀ := by
  classical
  set u : Fin r → Fin p → ℝ :=
    fun k => WithLp.ofLp (hA.eigenvectorBasis (eigIdx p (Fin.castLE hrp k))) with hu
  set L : Fin r → ℝ := fun k => hA.eigenvalues₀ (Fin.castLE hrp k) with hL
  set Q : Matrix (Fin p) (Fin r) ℝ := Matrix.of fun i k => u k i with hQ
  have hLpos : ∀ k, 0 < L k := by
    intro k
    have hk := k.isLt
    have hle : Fin.castLE hrp k ≤ (⟨r - 1, by omega⟩ : Fin (Fintype.card (Fin p))) := by
      rw [Fin.le_def]
      simp only [Fin.val_castLE]
      omega
    exact lt_of_lt_of_le hpos (hA.eigenvalues₀_antitone hle)
  have hQQ : Qᵀ * Q = 1 := by
    ext k l
    have hinner : ∑ i, Qᵀ k i * Q i l
        = ⟪hA.eigenvectorBasis (eigIdx p (Fin.castLE hrp k)),
           hA.eigenvectorBasis (eigIdx p (Fin.castLE hrp l))⟫_ℝ := by
      rw [real_inner_eq_dotProduct, dotProduct]
      exact Finset.sum_congr rfl fun i _ => rfl
    rw [Matrix.mul_apply, Matrix.one_apply, hinner]
    by_cases h : k = l
    · subst h
      rw [if_pos rfl, real_inner_self_eq_norm_sq, hA.eigenvectorBasis.orthonormal.1, one_pow]
    · rw [if_neg h]
      refine hA.eigenvectorBasis.orthonormal.2 fun hc => h ?_
      exact Fin.castLE_inj.mp ((eigIdx p).injective hc)
  have hAQ : A * Q = Q * Matrix.diagonal L := by
    ext i k
    have hev := congrFun (hA.mulVec_eigenvectorBasis (eigIdx p (Fin.castLE hrp k))) i
    simp only [Matrix.mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at hev
    rw [Matrix.mul_diagonal, Matrix.mul_apply, hL]
    simp only [hQ, hu, Matrix.of_apply]
    rw [← eigenvalues_eigIdx hA, mul_comm]
    exact hev
  have hQAQ : Qᵀ * A * Q = Matrix.diagonal L := by
    rw [Matrix.mul_assoc, hAQ, ← Matrix.mul_assoc, hQQ, Matrix.one_mul]
  have hspec : specInvTop A hA r = Q * Matrix.diagonal (fun k => (L k)⁻¹) * Qᵀ := by
    ext i j
    have hrhs : (Q * Matrix.diagonal (fun k => (L k)⁻¹) * Qᵀ) i j
        = ∑ k : Fin r, (L k)⁻¹ * (u k i * u k j) := by
      rw [Matrix.mul_assoc, Matrix.mul_apply]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [Matrix.diagonal_mul]
      simp only [hQ, hu, Matrix.of_apply, Matrix.transpose_apply]
      ring
    have hre : ∑ a : Fin p,
          (Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues a) *
            (WithLp.ofLp (hA.eigenvectorBasis a) i * WithLp.ofLp (hA.eigenvectorBasis a) j))
        = ∑ k : Fin (Fintype.card (Fin p)),
          (Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues (eigIdx p k)) *
            (WithLp.ofLp (hA.eigenvectorBasis (eigIdx p k)) i *
              WithLp.ofLp (hA.eigenvectorBasis (eigIdx p k)) j)) :=
      (Fintype.sum_equiv (eigIdx p) _ _ fun k => rfl).symm
    have hzero : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) →
        (Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues (eigIdx p k)) *
          (WithLp.ofLp (hA.eigenvectorBasis (eigIdx p k)) i *
            WithLp.ofLp (hA.eigenvectorBasis (eigIdx p k)) j)) = 0 := by
      intro k hk
      have hnm : hA.eigenvalues₀ k ∉ topEigSet A hA r := by
        rw [mem_topEigSet_iff hA hgap]
        omega
      rw [eigenvalues_eigIdx hA, Set.indicator_of_notMem hnm, zero_mul]
    rw [hrhs, specInvTop]
    simp only [Matrix.sum_apply, Matrix.smul_apply, Matrix.vecMulVec_apply, smul_eq_mul]
    rw [hre, sum_eq_sum_castLE hrp _ hzero]
    refine Finset.sum_congr rfl fun k _ => ?_
    have hk := k.isLt
    have hmem : hA.eigenvalues₀ (Fin.castLE hrp k) ∈ topEigSet A hA r := by
      rw [mem_topEigSet_iff hA hgap]
      simp
    rw [eigenvalues_eigIdx hA, Set.indicator_of_mem hmem]
  exact ⟨Q, L, hLpos, hQQ, hQAQ, hspec⟩


/-- The trace form at a matrix of the shape `C Cᵀ + I` with `C` of `r` columns. Every sorted
eigenvalue at or above index `r` is exactly `1`, so it contributes nothing and **no eigengap is
needed**: this is what lets `limitRW_optWR` drop the hypothesis `hgap`
(`notes/archive/audit_rank_r_weighted_2026-08-30.md`, finding 5.2). -/
theorem trace_specInvTop_add_one {r : ℕ} (C : Matrix (Fin p) (Fin r) ℝ)
    (hrp : r ≤ Fintype.card (Fin p)) (hA : (C * Cᵀ + 1).IsHermitian) :
    Matrix.trace (Cᵀ * specInvTop (C * Cᵀ + 1) hA r * C)
      = (r : ℝ) - ∑ k : Fin r, (hA.eigenvalues₀ (Fin.castLE hrp k))⁻¹ := by
  classical
  have hP : (C * Cᵀ).IsHermitian := by simpa using Matrix.isHermitian_mul_conjTranspose_self C
  have hPSD : (C * Cᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose C
  have hnn : ∀ k, 0 ≤ hP.eigenvalues₀ k := by
    intro k
    rw [← eigenvalues_eigIdx hP]
    exact hPSD.eigenvalues_nonneg _
  have hrank : (C * Cᵀ).rank ≤ r := by
    rw [Matrix.rank_self_mul_transpose]
    simpa using Matrix.rank_le_card_width C
  have hev : ∀ j, hA.eigenvalues₀ j = hP.eigenvalues₀ j + 1 := eigenvalues₀_add_one hP hA
  have hone : ∀ j : Fin (Fintype.card (Fin p)), r ≤ (j : ℕ) → hA.eigenvalues₀ j = 1 := by
    intro j hj
    rw [hev j, eigenvalues₀_eq_zero_of_rank_le hP hnn hrank hj, zero_add]
  have hge : ∀ j, 1 ≤ hA.eigenvalues₀ j := by
    intro j
    rw [hev j]
    linarith [hnn j]
  have hnorm : ∀ i : Fin p,
      ∑ j, (Cᵀ *ᵥ (WithLp.ofLp (hA.eigenvectorBasis i))) j ^ 2 = hA.eigenvalues i - 1 := by
    intro i
    rw [sum_sq_transpose_mulVec]
    have huu : (WithLp.ofLp (hA.eigenvectorBasis i)) ⬝ᵥ
        (WithLp.ofLp (hA.eigenvectorBasis i)) = 1 := by
      have h := real_inner_self_eq_norm_sq (hA.eigenvectorBasis i)
      rw [real_inner_eq_dotProduct] at h
      rw [h, hA.eigenvectorBasis.orthonormal.1, one_pow]
    have hmv := hA.mulVec_eigenvectorBasis i
    rw [Matrix.add_mulVec, Matrix.one_mulVec] at hmv
    rw [eq_sub_of_add_eq hmv, sub_dotProduct, smul_dotProduct, huu]
    simp
  rw [trace_specInvTop_conj]
  simp only [hnorm]
  have hre : ∑ i : Fin p,
        (Set.indicator (topEigSet (C * Cᵀ + 1) hA r) (fun t => t⁻¹) (hA.eigenvalues i) *
          (hA.eigenvalues i - 1))
      = ∑ k : Fin (Fintype.card (Fin p)),
        (Set.indicator (topEigSet (C * Cᵀ + 1) hA r) (fun t => t⁻¹) (hA.eigenvalues₀ k) *
          (hA.eigenvalues₀ k - 1)) := by
    refine (Fintype.sum_equiv (eigIdx p) _ _ fun k => ?_).symm
    rw [eigenvalues_eigIdx hA]
  have hzero : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) →
      (Set.indicator (topEigSet (C * Cᵀ + 1) hA r) (fun t => t⁻¹) (hA.eigenvalues₀ k) *
        (hA.eigenvalues₀ k - 1)) = 0 := by
    intro k hk
    rw [hone k hk, sub_self, mul_zero]
  rw [hre, sum_eq_sum_castLE hrp _ hzero]
  have hterm : ∀ k : Fin r,
      (Set.indicator (topEigSet (C * Cᵀ + 1) hA r) (fun t => t⁻¹)
          (hA.eigenvalues₀ (Fin.castLE hrp k)) * (hA.eigenvalues₀ (Fin.castLE hrp k) - 1))
        = 1 - (hA.eigenvalues₀ (Fin.castLE hrp k))⁻¹ := by
    intro k
    have hk := k.isLt
    have hmem : hA.eigenvalues₀ (Fin.castLE hrp k) ∈ topEigSet (C * Cᵀ + 1) hA r :=
      ⟨Fin.castLE hrp k, by simp, rfl⟩
    have hpos : (0 : ℝ) < hA.eigenvalues₀ (Fin.castLE hrp k) := by linarith [hge (Fin.castLE hrp k)]
    rw [Set.indicator_of_mem hmem]
    field_simp
  simp only [hterm]
  rw [Finset.sum_sub_distrib]
  simp

end Frame

/-! ### 4. The sorted spectrum of `C Cᵀ` from `Cᵀ C = diag S`

Item 2.3 of `notes/archive/rankr_D4_plan.md`. When the columns of `C ∈ ℝ^{p × q}` are orthogonal
with `Cᵀ C = diag(S)`, `S` antitone and nonnegative, and `q ≤ p`, the sorted spectrum of `C Cᵀ` is
`S` padded by zeros. Route: `Matrix.charpoly_mul_comm_of_le` gives `(C Cᵀ).charpoly = X^{p-q} (Cᵀ
C).charpoly`, `Matrix.charpoly_diagonal` evaluates the second factor, and
`eigenvalues₀_eq_of_charpoly` reads the sorted list off the product. The aligned model of
`thm:rank_r_svdstack` consumes it at `C = W⋆ B_R` (`RankR/AlignedComponent.lean`). -/

section GramSpectrum

variable {p q : ℕ}

/-- The padded list `S_0, …, S_{q-1}, 0, …, 0` is antitone when `S` is antitone and
nonnegative. -/
theorem antitone_padSpec {S : Fin q → ℝ} (hanti : Antitone S) (hnn : ∀ k, 0 ≤ S k) {n : ℕ} :
    Antitone fun k : Fin n => if h : (k : ℕ) < q then S ⟨k, h⟩ else 0 := by
  intro a b hab
  have hab' : (a : ℕ) ≤ (b : ℕ) := hab
  change (if h : (b : ℕ) < q then S ⟨(b : ℕ), h⟩ else 0)
      ≤ (if h : (a : ℕ) < q then S ⟨(a : ℕ), h⟩ else 0)
  by_cases hb : (b : ℕ) < q
  · have ha : (a : ℕ) < q := lt_of_le_of_lt hab' hb
    rw [dif_pos ha, dif_pos hb]
    exact hanti (show (⟨(a : ℕ), ha⟩ : Fin q) ≤ ⟨(b : ℕ), hb⟩ from hab')
  · rw [dif_neg hb]
    by_cases ha : (a : ℕ) < q
    · rw [dif_pos ha]
      exact hnn _
    · rw [dif_neg ha]

/-- The characteristic polynomial of the padded list splits as `X^{n-q}` times the product
over `Fin q`. -/
theorem prod_X_sub_C_padSpec {S : Fin q → ℝ} {n : ℕ} (hqn : q ≤ n) :
    (∏ k : Fin n, (Polynomial.X - Polynomial.C
        (if h : (k : ℕ) < q then S ⟨k, h⟩ else 0)) : Polynomial ℝ)
      = Polynomial.X ^ (n - q) * ∏ j : Fin q, (Polynomial.X - Polynomial.C (S j)) := by
  classical
  have hsplit : q + (n - q) = n := Nat.add_sub_cancel' hqn
  set F : Fin n → Polynomial ℝ :=
    fun k => Polynomial.X - Polynomial.C (if h : (k : ℕ) < q then S ⟨k, h⟩ else 0) with hF
  have hcongr : (∏ i : Fin (q + (n - q)), F (Fin.cast hsplit i)) = ∏ k : Fin n, F k :=
    Fin.prod_congr' F hsplit
  rw [← hcongr, Fin.prod_univ_add]
  have hleft : ∀ i : Fin q, F (Fin.cast hsplit (Fin.castAdd (n - q) i))
      = Polynomial.X - Polynomial.C (S i) := by
    intro i
    have hlt : ((Fin.cast hsplit (Fin.castAdd (n - q) i) : Fin n) : ℕ) < q := i.isLt
    rw [hF]
    simp only [dif_pos hlt]
    rfl
  have hright : ∀ i : Fin (n - q), F (Fin.cast hsplit (Fin.natAdd q i)) = Polynomial.X := by
    intro i
    have hnlt : ¬ ((Fin.cast hsplit (Fin.natAdd q i) : Fin n) : ℕ) < q := by
      have hval : ((Fin.cast hsplit (Fin.natAdd q i) : Fin n) : ℕ) = q + (i : ℕ) := rfl
      omega
    rw [hF]
    simp only [dif_neg hnlt, map_zero, sub_zero]
  rw [Finset.prod_congr rfl fun i _ => hleft i, Finset.prod_congr rfl fun i _ => hright i,
    Finset.prod_const, Finset.card_univ, Fintype.card_fin, mul_comm]

/-- **The sorted spectrum of `C Cᵀ` when `Cᵀ C` is diagonal** (item 2.3 of
`notes/archive/rankr_D4_plan.md`). With `q ≤ p`, `Cᵀ C = diag S`, `S` antitone and nonnegative, the
sorted eigenvalues of `C Cᵀ` are `S` on the first `q` indices and `0` beyond. -/
theorem eigenvalues₀_mul_transpose_of_transpose_mul_diagonal (C : Matrix (Fin p) (Fin q) ℝ)
    {S : Fin q → ℝ} (hC : Cᵀ * C = Matrix.diagonal S) (hanti : Antitone S)
    (hnn : ∀ k, 0 ≤ S k) (hqp : q ≤ p) (hH : (C * Cᵀ).IsHermitian)
    (k : Fin (Fintype.card (Fin p))) :
    hH.eigenvalues₀ k = if h : (k : ℕ) < q then S ⟨k, h⟩ else 0 := by
  classical
  have hqn : q ≤ Fintype.card (Fin p) := by simpa using hqp
  have hcard : Fintype.card (Fin p) - Fintype.card (Fin q) = Fintype.card (Fin p) - q := by
    simp
  have hchar : (C * Cᵀ).charpoly
      = ∏ l : Fin (Fintype.card (Fin p)), (Polynomial.X - Polynomial.C
          (if h : (l : ℕ) < q then S ⟨l, h⟩ else 0)) := by
    rw [prod_X_sub_C_padSpec hqn,
      Matrix.charpoly_mul_comm_of_le C Cᵀ (by simpa using hqp), hcard, hC,
      Matrix.charpoly_diagonal]
  exact congrFun (eigenvalues₀_eq_of_charpoly hH (antitone_padSpec hanti hnn) hchar) k

end GramSpectrum

end StackedSVD
