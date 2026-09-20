/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Eigen
import StackedSVD.SVDStack.Defs

/-!
# The spectral projector at one eigenvalue index

`Defs.lean` gives `topProj A hA = specProj A {lamMax A hA}` and `overlap X w = ‖topProj w‖²`;
`LinAlg/SpecProjPerturb.lean` gives the top-`r` subspace projector `specProjTop A hA r`. A
rank-`r_i` table needs the projector at **one** sorted eigenvalue index `k`, because the rows
of `Ṽ` are the individual singular vectors. This file states that projector and proves the
interface that Section 7 at general `r_i` consumes.

Task B1 of `notes/archive/rankr_plan_B.md` section 2. The six definitions (`eigSetIdx`,
`specProjIdx`, `overlapIdx`, `vEig`, `SimpleSpec`, `vEig_zero`) moved here verbatim from
`RankR/General.lean` on 2026-09-01; `General.lean` now imports this file.

## Content

1. The two regimes of the index. `eigSetIdx A hA k` is the singleton `{λ_k(A)}` for
   `k < card`, and `∅` otherwise, so `specProjIdx A hA k` is `0` out of range and no side
   condition is carried (`eigSetIdx_eq_singleton`, `eigSetIdx_eq_empty`,
   `specProjIdx_apply_of_not_lt`).
2. The `vEig` bridge. `vEig` reads `Matrix.IsHermitian.eigenvectorBasis`, while the lemmas of
   `Spectral.lean` read the eigenbasis of the symmetric operator `symmOp hA`. Mathlib defines
   the first as the second reindexed by `eigIdx`, so the two agree at matching indices with no
   sign correction (`eigenvectorBasis_eigIdx`, `vEig_eq`).
3. The unit eigenvector at an index: `norm_vEig`, `mem_specSpace_vEig`. Mirrors: `norm_vMax`
   and `mem_topSpace_vMax` (`SVDStack/Defs.lean`).
4. The rank-one collapse under `SimpleSpec`: `specProjIdx_eq_rankOne` and
   `overlapIdx_eq_inner_sq`. Mirrors: `topProj_eq_rankOne` (`SVDStack/Defs.lean`) and
   `overlap_eq_inner_sq` (`Spectral.lean`). Without simplicity the one-sided bound
   `overlapIdx_ge_inner_sq` still holds. Mirror: `overlap_ge_inner_sq`.
5. The two reductions to index `0`: `vEig_zero` and `overlapIdx_zero`.
6. `simpleSpec_one_of_topSimple`: `TopSimple A hA → SimpleSpec A hA 1`. Mirror: the step
   `huniq` inside `specInvTop_one_mulVec` (`RankR/Defs.lean`). The converse is not used.

## Route

`normSq_specProj` of `Spectral.lean` reads `‖specProj A S w‖²` in the sorted eigenbasis for an
arbitrary set `S`. The rank-one collapse then expands `w` in that basis
(`OrthonormalBasis.sum_repr'`) and drops every term except the one index that carries the
eigenvalue, so no `finrank` argument appears.
-/

open Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

section EigIdx

variable {p : ℕ}

/-! ### 1. Definitions (moved verbatim from `RankR/General.lean`) -/

/-- The eigenvalue at the sorted index `k`, as a set: `{λ_k(A)}` when `k < p`, else `∅`.
`eigenvalues₀` is antitone, so `k = 0` gives `{lamMax A hA}`. -/
def eigSetIdx (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (k : ℕ) : Set ℝ :=
  {t | ∃ q : Fin (Fintype.card (Fin p)), (q : ℕ) = k ∧ hA.eigenvalues₀ q = t}

/-- Orthogonal projector onto the eigenspace of the `k`-th sorted eigenvalue. At `k = 0` it is
`topProj`. -/
noncomputable def specProjIdx (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (k : ℕ) :
    EuclideanSpace ℝ (Fin p) →L[ℝ] EuclideanSpace ℝ (Fin p) :=
  specProj A (eigSetIdx A hA k)

/-- Squared overlap of the `k`-th right singular subspace of `X` with `w`. At `k = 0` it is
`overlap X w` of `Defs.lean`. -/
noncomputable def overlapIdx {q : ℕ} (X : Matrix (Fin q) (Fin p) ℝ) (k : ℕ)
    (w : EuclideanSpace ℝ (Fin p)) : ℝ :=
  ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) k w‖ ^ 2

/-- `overlapIdx` is a squared norm, so it is nonnegative. The tree used to carry three copies
of this line: `EdgeGlueR.overlapIdx_nonneg`, and a private `overlapIdx_nonneg'` in each of
`RankR/RMT/DelocAffineR.lean` and `RankR/GeneralMain.lean`. Second cleanup pass, 2026-09-02. -/
theorem overlapIdx_nonneg {q : ℕ} (X : Matrix (Fin q) (Fin p) ℝ) (k : ℕ)
    (w : EuclideanSpace ℝ (Fin p)) : 0 ≤ overlapIdx X k w := sq_nonneg _

/-- A unit eigenvector at the sorted index `k`, from `Matrix.IsHermitian.eigenvectorBasis`
through the index equivalence `eigIdx` of `LinAlg/Eigen.lean`. Junk value `0` when `p ≤ k`.
The sign is arbitrary; `UnalignedModelR.vhatG` fixes it by the paper's convention. -/
noncomputable def vEig (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (k : ℕ) :
    EuclideanSpace ℝ (Fin p) :=
  if h : k < Fintype.card (Fin p) then hA.eigenvectorBasis (eigIdx p ⟨k, h⟩) else 0

/-- Each of the top `rk` eigenvalues of `A` is simple. The rank-`rk` twin of `TopSimple`: it is
what makes `specProjIdx A hA k` a rank-one projector for `k < rk`, hence
`overlapIdx X k w = ⟪v̂_k, w⟫²`. -/
def SimpleSpec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (rk : ℕ) : Prop :=
  ∀ k l : Fin (Fintype.card (Fin p)), (k : ℕ) < rk → k ≠ l →
    hA.eigenvalues₀ l ≠ hA.eigenvalues₀ k

/-- `vEig A hA 0` is the top eigenvector `vMax A hA` of `SVDStack/Defs.lean`. -/
theorem vEig_zero (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) :
    vEig A hA 0 = vMax A hA := by
  rw [vEig, vMax, eigIdx]
  by_cases h : 0 < p
  · rw [dif_pos (by simpa using h), dif_pos h]
  · rw [dif_neg (by simpa using h), dif_neg h]

/-! ### 2. The two regimes of the index -/

/-- In range, the eigenvalue set at an index is a singleton. Mirror: the step `hsetiff` inside
`specInvTop_one_mulVec` (`RankR/Defs.lean`). -/
theorem eigSetIdx_eq_singleton (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : k < Fintype.card (Fin p)) : eigSetIdx A hA k = {hA.eigenvalues₀ ⟨k, hk⟩} := by
  ext t
  constructor
  · rintro ⟨q, hq, rfl⟩
    exact congrArg hA.eigenvalues₀ (Fin.val_injective hq)
  · rintro rfl
    exact ⟨⟨k, hk⟩, rfl, rfl⟩

/-- Out of range, the eigenvalue set at an index is empty. -/
theorem eigSetIdx_eq_empty (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : ¬ k < Fintype.card (Fin p)) : eigSetIdx A hA k = ∅ := by
  ext t
  simp only [eigSetIdx, Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
  rintro ⟨q, hq, -⟩
  exact hk (by rw [← hq]; exact q.isLt)

/-- The spectral subspace at the empty set of eigenvalues is trivial. -/
theorem specSpace_empty (A : Matrix (Fin p) (Fin p) ℝ) : specSpace A (∅ : Set ℝ) = ⊥ := by
  simp [specSpace]

/-- Out of range, the projector at an index sends every vector to `0`. -/
theorem specProjIdx_apply_of_not_lt (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    {k : ℕ} (hk : ¬ k < Fintype.card (Fin p)) (w : EuclideanSpace ℝ (Fin p)) :
    specProjIdx A hA k w = 0 := by
  have h : specSpace A (eigSetIdx A hA k) = ⊥ := by
    rw [eigSetIdx_eq_empty A hA hk]
    exact specSpace_empty A
  change (specSpace A (eigSetIdx A hA k)).starProjection w = 0
  rw [h, Submodule.starProjection_bot]
  simp

/-- Out of range, `vEig` takes its junk value `0`. -/
theorem vEig_of_not_lt (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : ¬ k < Fintype.card (Fin p)) : vEig A hA k = 0 := by
  rw [vEig, dif_neg hk]

/-! ### 3. `vEig` against the eigenbasis of the symmetric operator -/

/-- Mathlib defines `Matrix.IsHermitian.eigenvectorBasis` as the eigenbasis of the symmetric
operator, reindexed by `eigIdx`. The two therefore agree at matching indices, with no sign
correction. This is the bridge that every lemma of `Spectral.lean` needs, because those lemmas
read the operator eigenbasis and `vEig` reads the matrix one. -/
theorem eigenvectorBasis_eigIdx {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (i : Fin (Fintype.card (Fin p))) :
    hA.eigenvectorBasis (eigIdx p i)
      = (symmOp hA).eigenvectorBasis finrank_euclideanSpace i := by
  change ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).reindex (eigIdx p) (eigIdx p i)
      = _
  rw [OrthonormalBasis.reindex_apply, Equiv.symm_apply_apply]

/-- `vEig A hA k` in the eigenbasis of the symmetric operator. -/
theorem vEig_eq (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : k < Fintype.card (Fin p)) :
    vEig A hA k = (symmOp hA).eigenvectorBasis finrank_euclideanSpace ⟨k, hk⟩ := by
  rw [vEig, dif_pos hk, eigenvectorBasis_eigIdx]

/-- `vEig` is a unit vector in range. Mirror: `norm_vMax` (`SVDStack/Defs.lean`). -/
theorem norm_vEig (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : k < Fintype.card (Fin p)) : ‖vEig A hA k‖ = 1 := by
  rw [vEig, dif_pos hk]
  exact hA.eigenvectorBasis.norm_eq_one _

/-- `vEig` lies in the spectral subspace at its own index. Mirror: `mem_topSpace_vMax`
(`SVDStack/Defs.lean`). -/
theorem mem_specSpace_vEig (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) {k : ℕ}
    (hk : k < Fintype.card (Fin p)) : vEig A hA k ∈ specSpace A (eigSetIdx A hA k) := by
  rw [eigSetIdx_eq_singleton A hA hk, specSpace_singleton, vEig_eq A hA hk]
  exact Module.End.mem_eigenspace_iff.mpr (apply_eigvec hA ⟨k, hk⟩)

/-! ### 4. Index `0` -/

/-- `overlapIdx X 0 w` is `overlap X w` of `Defs.lean`: at a positive dimension the set
`eigSetIdx A hA 0` is the singleton `{lamMax A hA}`. -/
theorem overlapIdx_zero {q : ℕ} (hp : 0 < p) (X : Matrix (Fin q) (Fin p) ℝ)
    (w : EuclideanSpace ℝ (Fin p)) : overlapIdx X 0 w = overlap X w := by
  have hcard : 0 < Fintype.card (Fin p) := by simpa using hp
  have hset : eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) 0
      = {lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)} := by
    rw [eigSetIdx_eq_singleton _ _ hcard, lamMax, dif_pos hp]
  rw [overlapIdx, overlap, specProjIdx, hset]
  rfl

/-! ### 5. The rank-one collapse under `SimpleSpec` -/

/-- On a singleton spectral set whose eigenvalue is carried by exactly one sorted index, the
projector is the rank-one projector on that eigenvector. Proved by expanding `w` in the sorted
eigenbasis, so no `finrank` argument is needed. -/
private theorem specProj_singleton_eq_rankOne {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) (i₀ : Fin (Fintype.card (Fin p)))
    (huniq : ∀ i, hA.eigenvalues₀ i = hA.eigenvalues₀ i₀ → i = i₀)
    (w : EuclideanSpace ℝ (Fin p)) :
    specProj A {hA.eigenvalues₀ i₀} w
      = ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i₀, w⟫_ℝ •
        (symmOp hA).eigenvectorBasis finrank_euclideanSpace i₀ := by
  classical
  set b := (symmOp hA).eigenvectorBasis finrank_euclideanSpace with hb
  have hterm : ∀ i, specProj A {hA.eigenvalues₀ i₀} (b i) = if i = i₀ then b i else 0 := by
    intro i
    rw [hb, specProj_eigvec hA _ i, ← hb]
    by_cases hi : i = i₀
    · subst hi
      rw [if_pos (Set.mem_singleton_iff.mpr rfl), if_pos rfl]
    · rw [if_neg hi, if_neg (fun hc => hi (huniq i (Set.mem_singleton_iff.mp hc)))]
  calc specProj A {hA.eigenvalues₀ i₀} w
      = ∑ i, ⟪b i, w⟫_ℝ • specProj A {hA.eigenvalues₀ i₀} (b i) := by
        conv_lhs => rw [← b.sum_repr' w]
        rw [map_sum]
        exact Finset.sum_congr rfl fun i _ => map_smul _ _ _
    _ = ⟪b i₀, w⟫_ℝ • b i₀ := by
        rw [Finset.sum_eq_single i₀]
        · rw [hterm i₀, if_pos rfl]
        · intro i _ hi
          rw [hterm i, if_neg hi, smul_zero]
        · intro h
          exact absurd (Finset.mem_univ i₀) h

/-- Under `SimpleSpec A hA rk`, the projector at an index `k < rk` is the rank-one projector on
`vEig A hA k`. Mirror: `topProj_eq_rankOne` (`SVDStack/Defs.lean`). Out of range both sides are
`0`, so the statement carries no bound of `rk` by the dimension. -/
theorem specProjIdx_eq_rankOne {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {rk k : ℕ}
    (hsimple : SimpleSpec A hA rk) (hk : k < rk) (w : EuclideanSpace ℝ (Fin p)) :
    specProjIdx A hA k w = ⟪vEig A hA k, w⟫_ℝ • vEig A hA k := by
  by_cases hkc : k < Fintype.card (Fin p)
  · have huniq : ∀ i, hA.eigenvalues₀ i = hA.eigenvalues₀ ⟨k, hkc⟩ → i = ⟨k, hkc⟩ := by
      intro i hi
      by_contra hne
      exact hsimple ⟨k, hkc⟩ i hk (fun hc => hne hc.symm) hi
    rw [specProjIdx, eigSetIdx_eq_singleton A hA hkc,
      specProj_singleton_eq_rankOne hA ⟨k, hkc⟩ huniq w, ← vEig_eq A hA hkc]
  · rw [specProjIdx_apply_of_not_lt A hA hkc, vEig_of_not_lt A hA hkc]
    simp

/-- Under `SimpleSpec`, the projector form and the paper's inner-product form agree at every
index of the top `rk`. The sign ambiguity of `vEig` is absorbed by the square. Mirror:
`overlap_eq_inner_sq` (`Spectral.lean`). -/
theorem overlapIdx_eq_inner_sq {n : ℕ} (X : Matrix (Fin n) (Fin p) ℝ) {rk k : ℕ}
    (hsimple : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) rk) (hk : k < rk)
    (w : EuclideanSpace ℝ (Fin p)) :
    overlapIdx X k w
      = ⟪vEig (Xᵀ * X) (isHermitian_transpose_mul_self X) k, w⟫_ℝ ^ 2 := by
  rw [overlapIdx, specProjIdx_eq_rankOne hsimple hk, norm_smul, mul_pow, Real.norm_eq_abs,
    sq_abs]
  by_cases hkc : k < Fintype.card (Fin p)
  · rw [norm_vEig _ _ hkc, one_pow, mul_one]
  · rw [vEig_of_not_lt _ _ hkc]
    simp

/-- Without simplicity the projector form still dominates the inner-product form. This
transfers every `overlapIdx → 0` statement to the selected eigenvector. Mirror:
`overlap_ge_inner_sq` (`Spectral.lean`). -/
private theorem inner_sq_le_normSq_specProj {A : Matrix (Fin p) (Fin p) ℝ} (S : Set ℝ)
    {e : EuclideanSpace ℝ (Fin p)} (he : e ∈ specSpace A S) (hne : ‖e‖ = 1)
    (w : EuclideanSpace ℝ (Fin p)) : ⟪e, w⟫_ℝ ^ 2 ≤ ‖specProj A S w‖ ^ 2 := by
  set K := specSpace A S with hK
  have h1 : ⟪e, w⟫_ℝ = ⟪e, K.starProjection w⟫_ℝ := by
    rw [← Submodule.inner_starProjection_left_eq_right,
      Submodule.starProjection_eq_self_iff.mpr he]
  have h2 : |⟪e, w⟫_ℝ| ≤ ‖K.starProjection w‖ := by
    rw [h1]
    calc |⟪e, K.starProjection w⟫_ℝ| ≤ ‖e‖ * ‖K.starProjection w‖ := abs_real_inner_le_norm _ _
      _ = ‖K.starProjection w‖ := by rw [hne, one_mul]
  change _ ≤ ‖K.starProjection w‖ ^ 2
  nlinarith [abs_nonneg (⟪e, w⟫_ℝ), sq_abs (⟪e, w⟫_ℝ)]

/-- The one-sided bound, with no simplicity hypothesis. -/
theorem overlapIdx_ge_inner_sq {n : ℕ} (X : Matrix (Fin n) (Fin p) ℝ) (k : ℕ)
    (w : EuclideanSpace ℝ (Fin p)) :
    ⟪vEig (Xᵀ * X) (isHermitian_transpose_mul_self X) k, w⟫_ℝ ^ 2 ≤ overlapIdx X k w := by
  by_cases hkc : k < Fintype.card (Fin p)
  · exact inner_sq_le_normSq_specProj _ (mem_specSpace_vEig _ _ hkc) (norm_vEig _ _ hkc) w
  · rw [vEig_of_not_lt _ _ hkc, inner_zero_left, overlapIdx,
      zero_pow (by norm_num : (2 : ℕ) ≠ 0)]
    positivity

/-! ### 6. `TopSimple` gives `SimpleSpec` at `rk = 1` -/

/-- A simple top eigenvalue is `SimpleSpec` at `rk = 1`. Mirror: the step `huniq` inside
`specInvTop_one_mulVec` (`RankR/Defs.lean`). The converse is not used. -/
theorem simpleSpec_one_of_topSimple (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (h : TopSimple A hA) : SimpleSpec A hA 1 := by
  intro k l hk hkl hcon
  have hcard : 0 < Fintype.card (Fin p) := lt_of_le_of_lt (Nat.zero_le _) k.isLt
  have hp : 0 < p := by simpa using hcard
  have hk0 : (k : ℕ) = 0 := Nat.lt_one_iff.mp hk
  have hlam : hA.eigenvalues₀ k = lamMax A hA := by
    rw [lamMax, dif_pos hp]
    exact congrArg hA.eigenvalues₀ (Fin.val_injective hk0)
  have hmem : ∀ i : Fin (Fintype.card (Fin p)), hA.eigenvalues₀ i = lamMax A hA →
      (symmOp hA).eigenvectorBasis finrank_euclideanSpace i ∈ topSpace A hA := by
    intro i hi
    rw [topSpace_eq_eigenspace]
    refine Module.End.mem_eigenspace_iff.mpr ?_
    rw [← hi]
    exact apply_eigvec hA i
  have hnorm : ∀ i : Fin (Fintype.card (Fin p)),
      ‖(symmOp hA).eigenvectorBasis finrank_euclideanSpace i‖ = 1 :=
    fun i => ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).orthonormal.1 i
  have hproj := topProj_eq_rankOne h (hmem k hlam) (hnorm k)
  have h1 : topProj A hA ((symmOp hA).eigenvectorBasis finrank_euclideanSpace l)
      = (symmOp hA).eigenvectorBasis finrank_euclideanSpace l :=
    Submodule.starProjection_eq_self_iff.mpr (hmem l (hcon.trans hlam))
  have h2 := hproj ((symmOp hA).eigenvectorBasis finrank_euclideanSpace l)
  rw [((symmOp hA).eigenvectorBasis finrank_euclideanSpace).orthonormal.2 hkl, zero_smul] at h2
  have h3 := hnorm l
  rw [← h1, h2, norm_zero] at h3
  exact absurd h3 (by norm_num)

end EigIdx

end StackedSVD
