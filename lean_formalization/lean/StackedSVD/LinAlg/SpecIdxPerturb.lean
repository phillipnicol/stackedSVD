/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.SpecIdx
import StackedSVD.LinAlg.Frame

/-!
# Continuity of the spectral projector at one eigenvalue index

Items 3.1 to 3.4 and 4.3 of `notes/archive/rankr_D4_plan.md`. Everything here is deterministic
linear algebra, except the last section, which carries the one probabilistic statement that item 4.3
consumes.

`LinAlg/TopProjPerturb.lean` states Weyl and Davis-Kahan at the single index `0`;
`LinAlg/SpecProjPerturb.lean` states them for the projector on the top `r` eigenvalues. A
component-labelled result needs the projector at **one** free index `j`, which is
`specProjIdx` of `LinAlg/SpecIdx.lean`. This file connects the two: at a matrix with a top-`j`
gap and a top-`(j+1)` gap,

    ‖P_j y‖² = ‖P_{top (j+1)} y‖² - ‖P_{top j} y‖²,

and each term on the right is continuous in the entries of the matrix and in `y`. No window
projector and no threshold appears.

## Content

1. `eigVal`, the sorted eigenvalue at a free index, total (junk value `0` out of range).
   Matrix twin of `gramEig` (`LinAlg/SpecIdxMeas.lean`), which reads `Xᵀ X`. With
   `eigVal_eq`, `eigVal_of_le`, `eigVal_congr_mat`, the Weyl bound `abs_eigVal_sub_le` and the
   entrywise continuity `continuous_eigVal_symMat` (item 3.1).
2. The splitting of the index projector (item 3.2). `topEigSet_succ` writes
   `topEigSet A hA (j+1)` as the disjoint union of `topEigSet A hA j` and `eigSetIdx A hA j`,
   and `normSq_specProjIdx_eq_sub` is Pythagoras on that union
   (`Frame.norm_sq_specProj_union`). Only the gap at `j` is needed, not the gap at `j+1`.
3. The continuity of the top-`k` projector in the entries and in the vector (item 3.3).
   `normSq_specProjTop_continuousAt_fixed` is the fixed-vector form, a copy of
   `norm_topProj_sq_continuousAt` (`LinAlg/TopProjPerturb.lean`) with `specProjTop_perturb` in
   place of `topProj_perturb`; `normSq_specProjTop_continuousAt` lets the vector move.
   `eventually_topGap` says that a top-`k` gap survives a small entrywise perturbation.
4. The quotient functional of item 3.4. `compFun j jc` reads
   `‖P_j (column jc)‖² / λ_j` off a joint entry family on
   `(Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ`, the shape `TendstoInProbPi.comp_continuous`
   consumes; `continuousAt_compFun` is its continuity at a matrix with the two gaps and a
   nonzero `λ_j`. Mirror: `traceFun` and `continuousAt_traceFun` (`RankR/Defs.lean`).
5. Simplicity at one index (item 4.3). `SimpleIdx A hA j` says that no other sorted index
   carries `λ_j`; the two gaps give it (`simpleIdx_of_topGap`), it collapses the projector to
   a rank-one form (`normSq_specProjIdx_eq_inner_sq`), and it holds with probability tending
   to one under an entrywise limit in probability (`specIdxSimple_whp_of_tendsto`). Mirror:
   `specTop_simple_whp_of_tendsto` (`LinAlg/SpecProjPerturb.lean`).

## Two restatements

`topGap_of_gap` and the singleton rank-one collapse are private in `LinAlg/SpecProjPerturb.lean`
and `LinAlg/SpecIdx.lean`. This file does not edit either: `topGap_of_gap` is not needed (the
`TopGap` of the perturbed matrix comes out of `specProjTop_perturb` itself), and the rank-one
collapse is re-derived from the public `Frame.norm_sq_specProj_eq_sum`, which gives the squared
norm directly and is the only form the consumers need.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix Matrix.Norms.L2Operator

namespace StackedSVD

variable {p q : ℕ}

/-! ### 1. The sorted eigenvalue at a free index -/

/-- `λ_j(A)` at a free sorted index, total: junk value `0` out of range. Matrix twin of
`gramEig` (`LinAlg/SpecIdxMeas.lean`), which reads `Xᵀ X`, and of `vEig`
(`LinAlg/SpecIdx.lean`), which reads the eigenvector. -/
noncomputable def eigVal (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (j : ℕ) : ℝ :=
  if h : j < Fintype.card (Fin p) then hA.eigenvalues₀ ⟨j, h⟩ else 0

/-- In range, `eigVal` is the sorted eigenvalue. -/
theorem eigVal_eq (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {j : ℕ}
    (h : j < Fintype.card (Fin p)) : eigVal A hA j = hA.eigenvalues₀ ⟨j, h⟩ := dif_pos h

/-- Out of range, `eigVal` is `0`. -/
theorem eigVal_of_le (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {j : ℕ}
    (h : Fintype.card (Fin p) ≤ j) : eigVal A hA j = 0 := dif_neg (by omega)

/-- `eigVal` transports along a matrix equality. -/
theorem eigVal_congr_mat {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (j : ℕ) : eigVal A hA j = eigVal B hB j := by
  subst h; rfl

/-- Weyl at a free index: `eigVal` is `1`-Lipschitz in the `L2` operator norm. -/
theorem abs_eigVal_sub_le (A B : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (j : ℕ) : |eigVal A hA j - eigVal B hB j| ≤ ‖A - B‖ := by
  by_cases h : j < Fintype.card (Fin p)
  · rw [eigVal_eq A hA h, eigVal_eq B hB h]
    exact abs_eigenvalues₀_sub_le A B hA hB ⟨j, h⟩
  · rw [eigVal_of_le A hA (by omega), eigVal_of_le B hB (by omega), sub_zero, abs_zero]
    exact norm_nonneg _

/-- Value of the symmetrized `eigVal` at a symmetric matrix. -/
theorem eigVal_symMat {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (j : ℕ) :
    eigVal (symMat fun t => A t.1 t.2) (isHermitian_symMat _) j = eigVal A hA j :=
  eigVal_congr_mat _ hA (symMat_entries hA) j

/-- **Item 3.1.** `w ↦ λ_j(symMat w)` is continuous on `Fin p × Fin p → ℝ` (Weyl at a free
index, entrywise form). Mirror: `continuous_lamMax_symMat` (`LinAlg/TopProjPerturb.lean`),
which is this statement at `j = 0`. -/
theorem continuous_eigVal_symMat (j : ℕ) :
    Continuous (fun w : Fin p × Fin p → ℝ => eigVal (symMat w) (isHermitian_symMat w) j) := by
  refine Metric.continuous_iff.mpr fun w' ε hε =>
    ⟨ε / ((p : ℝ) + 1), by positivity, fun w hw => ?_⟩
  rw [Real.dist_eq]
  exact lt_of_le_of_lt (abs_eigVal_sub_le _ _ _ _ j) (norm_symMat_sub_lt hε hw)

/-- **Item 3.1**, in the `ContinuousAt` form `TendstoInProbPi.comp_continuous` consumes. -/
theorem continuousAt_eigVal_symMat (j : ℕ) (w : Fin p × Fin p → ℝ) :
    ContinuousAt (fun w : Fin p × Fin p → ℝ => eigVal (symMat w) (isHermitian_symMat w) j) w :=
  (continuous_eigVal_symMat j).continuousAt

/-! ### 2. The index projector as a difference of two top projectors -/

/-- `TopGap` at index `0` holds for every matrix: no pair of indices exists. -/
theorem topGap_zero (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) : TopGap A hA 0 :=
  fun _ _ hk _ => absurd hk (Nat.not_lt_zero _)

/-- The top-`0` eigenvalue set is empty. -/
theorem topEigSet_zero (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) :
    topEigSet A hA 0 = ∅ := by
  ext t
  simp only [topEigSet, Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
  rintro ⟨k, hk, -⟩
  exact absurd hk (Nat.not_lt_zero _)

/-- The top-`0` projector sends every vector to `0`. -/
theorem specProjTop_apply_zero (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (y : EuclideanSpace ℝ (Fin p)) : specProjTop A hA 0 y = 0 := by
  have h : specSpace A (topEigSet A hA 0) = ⊥ := by
    rw [topEigSet_zero A hA]
    exact specSpace_empty A
  change (specSpace A (topEigSet A hA 0)).starProjection y = 0
  rw [h, Submodule.starProjection_bot]
  simp

/-- The top `j + 1` eigenvalues are the top `j` together with the one at index `j`. No gap is
needed: out of range both sides collapse to the same set. -/
theorem topEigSet_succ (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (j : ℕ) :
    topEigSet A hA (j + 1) = topEigSet A hA j ∪ eigSetIdx A hA j := by
  ext t
  simp only [Set.mem_union]
  constructor
  · rintro ⟨k, hk, rfl⟩
    rcases Nat.lt_succ_iff_lt_or_eq.mp hk with h | h
    · exact Or.inl ⟨k, h, rfl⟩
    · exact Or.inr ⟨k, h, rfl⟩
  · rintro (⟨k, hk, rfl⟩ | ⟨k, hk, rfl⟩)
    · exact ⟨k, by omega, rfl⟩
    · exact ⟨k, by omega, rfl⟩

/-- Under the top-`j` gap the eigenvalue at index `j` is not one of the top `j`. -/
theorem disjoint_topEigSet_eigSetIdx {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {j : ℕ} (hgap : TopGap A hA j) : Disjoint (topEigSet A hA j) (eigSetIdx A hA j) := by
  rw [Set.disjoint_right]
  rintro t ⟨k, hk, rfl⟩ ⟨l, hl, hlt⟩
  have h := hgap l k hl (by omega)
  rw [hlt] at h
  exact lt_irrefl _ h

/-- Pythagoras at one index: the top-`(j+1)` projector splits into the top-`j` projector and
the projector at index `j`. -/
theorem normSq_specProjTop_succ {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {j : ℕ}
    (hgap : TopGap A hA j) (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProjTop A hA (j + 1) y‖ ^ 2
      = ‖specProjTop A hA j y‖ ^ 2 + ‖specProjIdx A hA j y‖ ^ 2 := by
  simp only [specProjTop, specProjIdx]
  rw [topEigSet_succ A hA j]
  exact Frame.norm_sq_specProj_union hA (disjoint_topEigSet_eigSetIdx hA hgap) y

/-- **Item 3.2.** The squared norm of the index-`j` projector is a difference of two top
projector norms. Only the gap at `j` is needed. -/
theorem normSq_specProjIdx_eq_sub {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {j : ℕ}
    (hgap : TopGap A hA j) (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx A hA j y‖ ^ 2
      = ‖specProjTop A hA (j + 1) y‖ ^ 2 - ‖specProjTop A hA j y‖ ^ 2 := by
  rw [normSq_specProjTop_succ hA hgap y]
  ring

/-! ### 3. Continuity of the top-`k` projector, entrywise and in the vector -/

/-- `specProjTop` transports along a matrix equality. -/
theorem specProjTop_congr_mat {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (k : ℕ) : specProjTop A hA k = specProjTop B hB k := by
  subst h; rfl

/-- `specProjIdx` transports along a matrix equality. -/
theorem specProjIdx_congr_mat {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (k : ℕ) : specProjIdx A hA k = specProjIdx B hB k := by
  subst h; rfl

/-- An orthogonal projector is a contraction. -/
theorem norm_specProj_le (A : Matrix (Fin p) (Fin p) ℝ) (S : Set ℝ)
    (y : EuclideanSpace ℝ (Fin p)) : ‖specProj A S y‖ ≤ ‖y‖ :=
  Submodule.norm_starProjection_apply_le _ y

/-- The difference of two squared norms, by the reverse triangle inequality. -/
theorem abs_normSq_sub_normSq_le {E : Type*} [NormedAddCommGroup E] (a b : E) :
    |‖a‖ ^ 2 - ‖b‖ ^ 2| ≤ ‖a - b‖ * (‖a‖ + ‖b‖) := by
  have h1 : |‖a‖ - ‖b‖| ≤ ‖a - b‖ := abs_norm_sub_norm_le a b
  have h2 : ‖a‖ ^ 2 - ‖b‖ ^ 2 = (‖a‖ - ‖b‖) * (‖a‖ + ‖b‖) := by ring
  rw [h2, abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ ‖a‖ + ‖b‖)]
  exact mul_le_mul_of_nonneg_right h1 (by positivity)

/-- A top-`k` gap survives a small entrywise perturbation. This replaces the private
`topGap_of_gap` of `LinAlg/SpecProjPerturb.lean`: the gap of the perturbed matrix comes out of
`specProjTop_perturb` itself. -/
theorem eventually_topGap {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {k : ℕ}
    (hgap : TopGap A hA k) :
    ∀ᶠ w in 𝓝 (fun t : Fin p × Fin p => A t.1 t.2),
      TopGap (symMat w) (isHermitian_symMat w) k := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · exact Eventually.of_forall fun w => topGap_zero _ _
  by_cases hkp : k < Fintype.card (Fin p)
  · have hrm : k - 1 < Fintype.card (Fin p) := by omega
    have hγ : 0 < hA.eigenvalues₀ ⟨k - 1, hrm⟩ - hA.eigenvalues₀ ⟨k, hkp⟩ := by
      have h := hgap ⟨k - 1, hrm⟩ ⟨k, hkp⟩ (show k - 1 < k by omega) le_rfl
      linarith
    set γ := hA.eigenvalues₀ ⟨k - 1, hrm⟩ - hA.eigenvalues₀ ⟨k, hkp⟩ with hγdef
    have hopen : IsOpen {w : Fin p × Fin p → ℝ | ‖symMat w - A‖ < γ / 3} :=
      isOpen_lt (continuous_norm_symMat_sub A) continuous_const
    have hmem : (fun t : Fin p × Fin p => A t.1 t.2) ∈
        {w : Fin p × Fin p → ℝ | ‖symMat w - A‖ < γ / 3} := by
      change ‖symMat (fun t : Fin p × Fin p => A t.1 t.2) - A‖ < γ / 3
      rw [symMat_entries hA, sub_self, norm_zero]
      linarith
    filter_upwards [hopen.mem_nhds hmem] with w hw
    have hw' : ‖symMat w - A‖ < γ / 3 := hw
    obtain ⟨-, h, -⟩ := specProjTop_perturb hA (isHermitian_symMat w) hk (le_of_lt hkp) hγ
      (fun _ => le_rfl) (le_refl ‖symMat w - A‖) (by linarith)
    exact h
  · exact Eventually.of_forall fun w => topGap_of_card_le _ (by omega)

/-- A projector moves the squared norm by at most `‖u - v‖ (‖u - v‖ + 2 ‖v‖)`. -/
theorem abs_normSq_specProj_sub_le (A : Matrix (Fin p) (Fin p) ℝ) (S : Set ℝ)
    (u v : EuclideanSpace ℝ (Fin p)) :
    |‖specProj A S u‖ ^ 2 - ‖specProj A S v‖ ^ 2| ≤ ‖u - v‖ * (‖u - v‖ + 2 * ‖v‖) := by
  have h0 : specProj A S u - specProj A S v = specProj A S (u - v) := (map_sub _ _ _).symm
  have h1 : ‖specProj A S u - specProj A S v‖ ≤ ‖u - v‖ := by
    rw [h0]
    exact norm_specProj_le _ _ _
  have h2 : ‖specProj A S u‖ ≤ ‖u‖ := norm_specProj_le _ _ _
  have h3 : ‖specProj A S v‖ ≤ ‖v‖ := norm_specProj_le _ _ _
  have h4 : ‖u‖ ≤ ‖u - v‖ + ‖v‖ := by simpa using norm_add_le (u - v) v
  calc |‖specProj A S u‖ ^ 2 - ‖specProj A S v‖ ^ 2|
      ≤ ‖specProj A S u - specProj A S v‖ * (‖specProj A S u‖ + ‖specProj A S v‖) :=
        abs_normSq_sub_normSq_le _ _
    _ ≤ ‖u - v‖ * (‖u - v‖ + 2 * ‖v‖) :=
        mul_le_mul h1 (by linarith) (by positivity) (norm_nonneg _)

/-- **Item 3.3**, fixed vector. `w ↦ ‖P_{top k}(symMat w) x‖²` is continuous at the entries of
a matrix with a top-`k` gap. Copy of `norm_topProj_sq_continuousAt`
(`LinAlg/TopProjPerturb.lean`) with `specProjTop_perturb` in place of `topProj_perturb`. -/
theorem normSq_specProjTop_continuousAt_fixed {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {k : ℕ} (hgap : TopGap A hA k) (x : EuclideanSpace ℝ (Fin p)) :
    ContinuousAt (fun w : Fin p × Fin p → ℝ =>
      ‖specProjTop (symMat w) (isHermitian_symMat w) k x‖ ^ 2)
      (fun t => A t.1 t.2) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · have hconst : (fun w : Fin p × Fin p → ℝ =>
        ‖specProjTop (symMat w) (isHermitian_symMat w) 0 x‖ ^ 2) = fun _ => (0 : ℝ) := by
      funext w
      rw [specProjTop_apply_zero, norm_zero]
      norm_num
    rw [hconst]
    exact continuousAt_const
  by_cases hkp : k < Fintype.card (Fin p)
  · have hrm : k - 1 < Fintype.card (Fin p) := by omega
    have hγ : 0 < hA.eigenvalues₀ ⟨k - 1, hrm⟩ - hA.eigenvalues₀ ⟨k, hkp⟩ := by
      have h := hgap ⟨k - 1, hrm⟩ ⟨k, hkp⟩ (show k - 1 < k by omega) le_rfl
      linarith
    set γ := hA.eigenvalues₀ ⟨k - 1, hrm⟩ - hA.eigenvalues₀ ⟨k, hkp⟩ with hγdef
    refine Metric.continuousAt_iff.mpr fun ε hε => ?_
    set δ := min (γ / 3) (ε * γ / (8 * (‖x‖ ^ 2 + 1))) with hδdef
    have hδ0 : 0 < δ := lt_min (by linarith) (by positivity)
    refine ⟨δ / ((p : ℝ) + 1), by positivity, fun {w} hw => ?_⟩
    have hop : ‖symMat w - symMat fun t : Fin p × Fin p => A t.1 t.2‖ < δ :=
      norm_symMat_sub_lt hδ0 hw
    rw [symMat_entries hA] at hop
    have hd0 : (0 : ℝ) ≤ ‖symMat w - A‖ := norm_nonneg _
    have hδ1 : ‖symMat w - A‖ < γ / 3 := lt_of_lt_of_le hop (min_le_left _ _)
    have hδ2 : ‖symMat w - A‖ < ε * γ / (8 * (‖x‖ ^ 2 + 1)) :=
      lt_of_lt_of_le hop (min_le_right _ _)
    obtain ⟨-, -, hproj⟩ := specProjTop_perturb hA (isHermitian_symMat w) hk (le_of_lt hkp) hγ
      (fun _ => le_rfl) (le_refl ‖symMat w - A‖) (by linarith)
    rw [Real.dist_eq,
      show ‖specProjTop (symMat fun t : Fin p × Fin p => A t.1 t.2)
          (isHermitian_symMat _) k x‖ ^ 2 = ‖specProjTop A hA k x‖ ^ 2 by
        rw [specProjTop_congr_mat (isHermitian_symMat _) hA (symMat_entries hA) k]]
    set P := specProjTop (symMat w) (isHermitian_symMat w) k with hP
    set Q := specProjTop A hA k with hQ
    have hbound : |‖P x‖ ^ 2 - ‖Q x‖ ^ 2| ≤ ‖P x - Q x‖ * (‖P x‖ + ‖Q x‖) :=
      abs_normSq_sub_normSq_le _ _
    have h1 : ‖P x - Q x‖ ≤ ‖P - Q‖ * ‖x‖ := by
      have hsub : (P - Q) x = P x - Q x := rfl
      rw [← hsub]
      exact ContinuousLinearMap.le_opNorm _ _
    have h2 : ‖P x‖ ≤ ‖x‖ := by rw [hP]; exact norm_specProj_le _ _ _
    have h3 : ‖Q x‖ ≤ ‖x‖ := by rw [hQ]; exact norm_specProj_le _ _ _
    have hA1 : ‖P x - Q x‖ ≤ 4 * ‖symMat w - A‖ / γ * ‖x‖ :=
      le_trans h1 (mul_le_mul_of_nonneg_right hproj (norm_nonneg x))
    have hprod : ‖P x - Q x‖ * (‖P x‖ + ‖Q x‖)
        ≤ 4 * ‖symMat w - A‖ / γ * ‖x‖ * (2 * ‖x‖) :=
      mul_le_mul hA1 (by linarith) (by positivity) (by positivity)
    have heq : 4 * ‖symMat w - A‖ / γ * ‖x‖ * (2 * ‖x‖)
        = 8 * ‖symMat w - A‖ * ‖x‖ ^ 2 / γ := by
      field_simp
      ring
    rw [heq] at hprod
    have hfrac : 8 * ‖symMat w - A‖ * (‖x‖ ^ 2 + 1) < ε * γ := by
      have hpos : (0 : ℝ) < 8 * (‖x‖ ^ 2 + 1) := by positivity
      calc 8 * ‖symMat w - A‖ * (‖x‖ ^ 2 + 1)
          = ‖symMat w - A‖ * (8 * (‖x‖ ^ 2 + 1)) := by ring
        _ < ε * γ / (8 * (‖x‖ ^ 2 + 1)) * (8 * (‖x‖ ^ 2 + 1)) :=
            mul_lt_mul_of_pos_right hδ2 hpos
        _ = ε * γ := by field_simp
    have hlt : 8 * ‖symMat w - A‖ * ‖x‖ ^ 2 / γ < ε := by
      rw [div_lt_iff₀ hγ]
      nlinarith [sq_nonneg ‖x‖, hd0]
    linarith
  · have hle : Fintype.card (Fin p) ≤ k := by omega
    have hconst : (fun w : Fin p × Fin p → ℝ =>
        ‖specProjTop (symMat w) (isHermitian_symMat w) k x‖ ^ 2) = fun _ => ‖x‖ ^ 2 := by
      funext w
      rw [specProjTop_eq_id (isHermitian_symMat w) hle]
      simp
    rw [hconst]
    exact continuousAt_const

/-- **Item 3.3.** `(w, y) ↦ ‖P_{top k}(symMat w) y‖²` is continuous at the entries of a matrix
with a top-`k` gap and at any vector. -/
theorem normSq_specProjTop_continuousAt {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {k : ℕ} (hgap : TopGap A hA k) (y₀ : Fin p → ℝ) :
    ContinuousAt (fun z : (Fin p × Fin p → ℝ) × (Fin p → ℝ) =>
        ‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k (WithLp.toLp 2 z.2)‖ ^ 2)
      ((fun t : Fin p × Fin p => A t.1 t.2), y₀) := by
  set x₀ : EuclideanSpace ℝ (Fin p) := WithLp.toLp 2 y₀ with hx₀
  have hx0nn : (0 : ℝ) ≤ ‖x₀‖ := norm_nonneg _
  have hfix := normSq_specProjTop_continuousAt_fixed hA hgap x₀
  have hcontLp : Continuous (fun y : Fin p → ℝ =>
      (WithLp.toLp 2 y : EuclideanSpace ℝ (Fin p))) := PiLp.continuous_toLp _ _
  refine Metric.continuousAt_iff.mpr fun ε hε => ?_
  obtain ⟨δ₁, hδ₁, h1⟩ := Metric.continuousAt_iff.mp hfix (ε / 2) (by linarith)
  set σ := min 1 (ε / (2 * (2 * ‖x₀‖ + 1))) with hσdef
  have hσ0 : 0 < σ := lt_min one_pos (by positivity)
  have hσ1 : σ ≤ 1 := min_le_left _ _
  have hσ2 : σ ≤ ε / (2 * (2 * ‖x₀‖ + 1)) := min_le_right _ _
  have hσε : σ * (2 * ‖x₀‖ + 1) ≤ ε / 2 := by
    calc σ * (2 * ‖x₀‖ + 1) ≤ ε / (2 * (2 * ‖x₀‖ + 1)) * (2 * ‖x₀‖ + 1) :=
          mul_le_mul_of_nonneg_right hσ2 (by positivity)
      _ = ε / 2 := by field_simp
  obtain ⟨δ₂, hδ₂, h2⟩ := Metric.continuous_iff.mp hcontLp y₀ σ hσ0
  refine ⟨min δ₁ δ₂, lt_min hδ₁ hδ₂, fun {z} hz => ?_⟩
  rw [Prod.dist_eq, max_lt_iff] at hz
  have hz1 : dist z.1 (fun t : Fin p × Fin p => A t.1 t.2) < δ₁ :=
    lt_of_lt_of_le hz.1 (min_le_left _ _)
  have hz2 : dist z.2 y₀ < δ₂ := lt_of_lt_of_le hz.2 (min_le_right _ _)
  have hstep1 := h1 hz1
  rw [Real.dist_eq] at hstep1 ⊢
  have he : ‖(WithLp.toLp 2 z.2 : EuclideanSpace ℝ (Fin p)) - x₀‖ < σ := by
    rw [← dist_eq_norm]
    exact h2 _ hz2
  have he0 : (0 : ℝ) ≤ ‖(WithLp.toLp 2 z.2 : EuclideanSpace ℝ (Fin p)) - x₀‖ := norm_nonneg _
  have hfirst : |‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k (WithLp.toLp 2 z.2)‖ ^ 2
      - ‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k x₀‖ ^ 2| < ε / 2 := by
    have hb := abs_normSq_specProj_sub_le (symMat z.1)
      (topEigSet (symMat z.1) (isHermitian_symMat z.1) k) (WithLp.toLp 2 z.2) x₀
    have hb' : |‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k (WithLp.toLp 2 z.2)‖ ^ 2
        - ‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k x₀‖ ^ 2|
        ≤ ‖(WithLp.toLp 2 z.2 : EuclideanSpace ℝ (Fin p)) - x₀‖ *
          (‖(WithLp.toLp 2 z.2 : EuclideanSpace ℝ (Fin p)) - x₀‖ + 2 * ‖x₀‖) := hb
    nlinarith [hb', he, he0, hσ1, hσε, hx0nn, hσ0]
  have htri := abs_sub_le
    (‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k (WithLp.toLp 2 z.2)‖ ^ 2)
    (‖specProjTop (symMat z.1) (isHermitian_symMat z.1) k x₀‖ ^ 2)
    (‖specProjTop (symMat fun t : Fin p × Fin p => A t.1 t.2) (isHermitian_symMat _) k x₀‖ ^ 2)
  linarith

/-! ### 4. The quotient functional on a joint entry family -/

/-- `z ↦ ‖P_{top k}(symMat z_inl) (column jc of z_inr)‖²`, read off the joint entry family. -/
noncomputable def projTopFun (k : ℕ) (jc : Fin q)
    (z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ) : ℝ :=
  ‖specProjTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) k
      (WithLp.toLp 2 fun i => z (Sum.inr (i, jc)))‖ ^ 2

/-- `z ↦ ‖P_j(symMat z_inl) (column jc of z_inr)‖²`, read off the joint entry family. -/
noncomputable def projIdxFun (j : ℕ) (jc : Fin q)
    (z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ) : ℝ :=
  ‖specProjIdx (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) j
      (WithLp.toLp 2 fun i => z (Sum.inr (i, jc)))‖ ^ 2

/-- `z ↦ λ_j(symMat z_inl)`, read off the joint entry family. -/
noncomputable def eigValFun (j : ℕ) (z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ) : ℝ :=
  eigVal (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) j

/-- **The functional of item 3.4:** `‖P_j (column jc)‖² / λ_j`, read off the joint entry family
`z` on `(Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ`. Same domain shape as `traceFun`
(`RankR/Defs.lean`), which is what `TendstoInProbPi.comp_continuous` consumes. -/
noncomputable def compFun (j : ℕ) (jc : Fin q)
    (z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ) : ℝ :=
  projIdxFun j jc z / eigValFun j z

/-- `projTopFun` at the entry family of a symmetric `C` and an arbitrary `Y`. -/
theorem projTopFun_eq {C : Matrix (Fin p) (Fin p) ℝ} (hC : C.IsHermitian)
    (Y : Matrix (Fin p) (Fin q) ℝ) (k : ℕ) (jc : Fin q) :
    projTopFun (p := p) k jc (Sum.elim (fun t : Fin p × Fin p => C t.1 t.2)
        (fun t : Fin p × Fin q => Y t.1 t.2))
      = ‖specProjTop C hC k (WithLp.toLp 2 fun i => Y i jc)‖ ^ 2 := by
  simp only [projTopFun, Sum.elim_inl, Sum.elim_inr]
  rw [specProjTop_congr_mat (isHermitian_symMat _) hC (symMat_entries hC) k]

/-- `compFun` at the entry family of a symmetric `C` and an arbitrary `Y`. -/
theorem compFun_eq {C : Matrix (Fin p) (Fin p) ℝ} (hC : C.IsHermitian)
    (Y : Matrix (Fin p) (Fin q) ℝ) (j : ℕ) (jc : Fin q) :
    compFun (p := p) j jc (Sum.elim (fun t : Fin p × Fin p => C t.1 t.2)
        (fun t : Fin p × Fin q => Y t.1 t.2))
      = ‖specProjIdx C hC j (WithLp.toLp 2 fun i => Y i jc)‖ ^ 2 / eigVal C hC j := by
  simp only [compFun, projIdxFun, eigValFun, Sum.elim_inl, Sum.elim_inr]
  rw [specProjIdx_congr_mat (isHermitian_symMat _) hC (symMat_entries hC) j,
    eigVal_congr_mat (isHermitian_symMat _) hC (symMat_entries hC) j]

/-- `projTopFun k jc` is continuous at the joint entry family of a matrix with a top-`k` gap. -/
theorem continuousAt_projTopFun {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {k : ℕ}
    (hgap : TopGap A hA k) (B : Matrix (Fin p) (Fin q) ℝ) (jc : Fin q) :
    ContinuousAt (projTopFun (p := p) k jc)
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2)) := by
  have hsplit : Continuous (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
      ((fun t : Fin p × Fin p => z (Sum.inl t)), (fun i : Fin p => z (Sum.inr (i, jc))))) :=
    (continuous_pi fun t => continuous_apply (Sum.inl t)).prodMk
      (continuous_pi fun i => continuous_apply (Sum.inr (i, jc)))
  exact ContinuousAt.comp (x := Sum.elim (fun t : Fin p × Fin p => A t.1 t.2)
      (fun t : Fin p × Fin q => B t.1 t.2))
    (normSq_specProjTop_continuousAt hA hgap fun i => B i jc) hsplit.continuousAt

/-- `projIdxFun j jc` is continuous at the joint entry family of a matrix with gaps at `j` and
at `j + 1`. -/
theorem continuousAt_projIdxFun {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {j : ℕ}
    (hgapj : TopGap A hA j) (hgapj1 : TopGap A hA (j + 1)) (B : Matrix (Fin p) (Fin q) ℝ)
    (jc : Fin q) :
    ContinuousAt (projIdxFun (p := p) j jc)
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2)) := by
  have hinl : Continuous (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
      (fun t : Fin p × Fin p => z (Sum.inl t))) := continuous_pi fun t => continuous_apply _
  have hev := (hinl.continuousAt (x := Sum.elim (fun t : Fin p × Fin p => A t.1 t.2)
      (fun t : Fin p × Fin q => B t.1 t.2))).eventually (eventually_topGap hA hgapj)
  refine ((continuousAt_projTopFun hA hgapj1 B jc).sub
    (continuousAt_projTopFun hA hgapj B jc)).congr ?_
  filter_upwards [hev] with z hz
  simp only [projIdxFun]
  exact (normSq_specProjIdx_eq_sub (isHermitian_symMat _) hz _).symm

/-- **Item 3.4.** The quotient functional is continuous at the joint entry family of a matrix
with gaps at `j` and at `j + 1` and a nonzero eigenvalue at index `j`. -/
theorem continuousAt_compFun {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {j : ℕ}
    (hgapj : TopGap A hA j) (hgapj1 : TopGap A hA (j + 1)) (hpos : 0 < eigVal A hA j)
    (B : Matrix (Fin p) (Fin q) ℝ) (jc : Fin q) :
    ContinuousAt (compFun (p := p) j jc)
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2)) := by
  have hinl : Continuous (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
      (fun t : Fin p × Fin p => z (Sum.inl t))) := continuous_pi fun t => continuous_apply _
  have hden : ContinuousAt (eigValFun (p := p) (q := q) j)
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2)) :=
    (continuous_eigVal_symMat j).continuousAt.comp hinl.continuousAt
  have hne : eigValFun (p := p) (q := q) j
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2))
      ≠ 0 := by
    change eigVal (symMat fun t : Fin p × Fin p => A t.1 t.2) (isHermitian_symMat _) j ≠ 0
    rw [eigVal_symMat hA]
    exact ne_of_gt hpos
  exact (continuousAt_projIdxFun hA hgapj hgapj1 B jc).div hden hne

/-! ### 5. Simplicity at one index -/

/-- The eigenvalue at the sorted index `j` is simple: no other sorted index carries it. The
one-index twin of `SimpleSpec` (`LinAlg/SpecIdx.lean`), which asks the same of every index
below `rk`. -/
def SimpleIdx (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (j : ℕ) : Prop :=
  ∀ k l : Fin (Fintype.card (Fin p)), (k : ℕ) = j → k ≠ l →
    hA.eigenvalues₀ l ≠ hA.eigenvalues₀ k

/-- Gaps at `j` and at `j + 1` make the eigenvalue at index `j` simple. -/
theorem simpleIdx_of_topGap {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {j : ℕ}
    (hgapj : TopGap A hA j) (hgapj1 : TopGap A hA (j + 1)) : SimpleIdx A hA j := by
  intro k l hk hkl hcon
  rcases lt_or_ge (l : ℕ) j with hl | hl
  · have h := hgapj l k hl (by omega)
    rw [hcon] at h
    exact lt_irrefl _ h
  · have hlj : j + 1 ≤ (l : ℕ) := by
      rcases eq_or_lt_of_le hl with h | h
      · exact absurd (Fin.val_injective (show (k : ℕ) = (l : ℕ) by omega)) hkl
      · omega
    have h := hgapj1 k l (by omega) hlj
    rw [hcon] at h
    exact lt_irrefl _ h

/-- Under simplicity at index `j` the projector collapses to the rank-one form: the paper's
`(v_jᵀ y)²`. The sign of `vEig` is absorbed by the square. Mirror: `overlapIdx_eq_inner_sq`
(`LinAlg/SpecIdx.lean`), which asks for `SimpleSpec`. -/
theorem normSq_specProjIdx_eq_inner_sq {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {j : ℕ} (hsimple : SimpleIdx A hA j) (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx A hA j y‖ ^ 2 = ⟪vEig A hA j, y⟫_ℝ ^ 2 := by
  classical
  by_cases hj : j < Fintype.card (Fin p)
  · have hveq : vEig A hA j = hA.eigenvectorBasis (eigIdx p ⟨j, hj⟩) := by
      rw [vEig, dif_pos hj]
    have hval : ∀ i : Fin p, hA.eigenvalues i = hA.eigenvalues₀ ((eigIdx p).symm i) := by
      intro i
      conv_lhs => rw [← Equiv.apply_symm_apply (eigIdx p) i]
      exact eigenvalues_eigIdx hA _
    have hmem : ∀ i : Fin p,
        hA.eigenvalues i = hA.eigenvalues₀ ⟨j, hj⟩ ↔ i = eigIdx p ⟨j, hj⟩ := by
      intro i
      constructor
      · intro hi
        rw [hval i] at hi
        by_contra hne
        exact hsimple ⟨j, hj⟩ ((eigIdx p).symm i) rfl
          (fun hc => hne (by rw [hc, Equiv.apply_symm_apply])) hi
      · rintro rfl
        exact eigenvalues_eigIdx hA _
    simp only [specProjIdx]
    rw [eigSetIdx_eq_singleton A hA hj, Frame.norm_sq_specProj_eq_sum hA,
      Finset.sum_eq_single_of_mem (eigIdx p ⟨j, hj⟩) (Finset.mem_univ _)]
    · rw [Set.indicator_of_mem (Set.mem_singleton_iff.mpr ((hmem _).mpr rfl)), hveq]
    · intro i _ hi
      exact Set.indicator_of_notMem
        (fun hc => hi ((hmem i).mp (Set.mem_singleton_iff.mp hc))) _
  · rw [specProjIdx_apply_of_not_lt A hA hj, vEig_of_not_lt A hA hj, inner_zero_left,
      norm_zero]

section EntrywiseLimit

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
variable {A : Matrix (Fin p) (Fin p) ℝ} {G : (N : ℕ) → Ω N → Matrix (Fin p) (Fin p) ℝ}

/-- `specTop_simple_whp_of_tendsto` (`LinAlg/SpecProjPerturb.lean`) with the hypothesis
`0 < k` removed: at `k = 0` the gap holds for every matrix. -/
theorem topGap_whp_of_tendsto (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    {k : ℕ} (hkp : k ≤ Fintype.card (Fin p)) (hgap : TopGap A hA k) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (G N ω) (hsymm N ω) k}) atTop (𝓝 0) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · have hempty : ∀ N, {ω : Ω N | ¬ TopGap (G N ω) (hsymm N ω) 0} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact topGap_zero _ _
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds
  · exact specTop_simple_whp_of_tendsto hA hsymm hconv hk hkp hgap

/-- **Item 4.3.** The free-index mirror of `specTop_simple_whp_of_tendsto`
(`LinAlg/SpecProjPerturb.lean`): if the entries of a random symmetric matrix converge in
probability to those of `A`, and `A` has gaps at `j` and at `j + 1`, then the eigenvalue at
sorted index `j` of the random matrix is simple with probability tending to one. -/
theorem specIdxSimple_whp_of_tendsto (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    {j : ℕ} (hjp : j < Fintype.card (Fin p)) (hgapj : TopGap A hA j)
    (hgapj1 : TopGap A hA (j + 1)) :
    Tendsto (fun N => μ N {ω | ¬ SimpleIdx (G N ω) (hsymm N ω) j}) atTop (𝓝 0) := by
  refine tendsto_measure_zero_of_subset (t := fun N =>
    {ω | ¬ TopGap (G N ω) (hsymm N ω) j} ∪ {ω | ¬ TopGap (G N ω) (hsymm N ω) (j + 1)})
    (fun N ω hω => ?_)
    (tendsto_measure_zero_union (topGap_whp_of_tendsto hA hsymm hconv (by omega) hgapj)
      (topGap_whp_of_tendsto hA hsymm hconv (by omega) hgapj1))
  by_contra hnot
  simp only [Set.mem_union, Set.mem_ofPred_eq, not_or, not_not] at hnot
  exact hω (simpleIdx_of_topGap (hsymm N ω) hnot.1 hnot.2)

end EntrywiseLimit

end StackedSVD
