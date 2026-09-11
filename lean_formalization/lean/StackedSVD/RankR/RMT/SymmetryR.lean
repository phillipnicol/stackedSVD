/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.DelocR

/-!
# Task U6: the top-`r` overlap of a direction orthogonal to the signal rows

Task U6 of `notes/archive/rankr_plan_A.md` section 3, the rank-`r` twin of `RMT/R6.lean` sections 2
and 4. The consumer is task U7 part 2, which applies this bound at the component of `V q_j`
orthogonal to the column space of `Q`.

The headline, for the `(r + p) × d` matrix with `r` fixed rows `Y` on top and a scaled
Gaussian block `t B` below:

```
∫⁻ B, ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B) _ r x‖ ^ 2) ∂(gaussianMatrix p d)
  ≤ ofReal (r / (d - r))
```

for every `x` with `‖x‖ ≤ 1` that is orthogonal to every row of `Y`. At `r = 1` the constant is
`1 / (d - 1)`, the rank-one constant of `R6.lintegral_overlap_Ymat_le`.

## The argument

1. **The model matrix** (section 1). `YmatR Y t B` is `Matrix.fromRows Y (t • B)`, so the row
   index is `Fin r ⊕ Fin p`. The Gram matrix is `Yᵀ Y + t² Bᵀ B` (`gram_YmatR`), right
   multiplication distributes over the two blocks (`YmatR_mul`), and the map `B ↦ YmatR Y t B`
   is measurable (`measurable_YmatR`, which U7 needs). A row reindexing does not change the
   Gram matrix, so the sum index type costs the consumer nothing.
2. **Conjugation** (section 2). `topEigSet` is defined from `eigenvalues₀`, which
   `DelocR.eigenvalues₀_conj` shows is invariant under an orthogonal conjugation, so the top-`r`
   projector conjugates (`specProjTop_conj`) and a right rotation moves the test direction
   (`normSq_specProjTop_mul_right`). Twin of `DelocR.normSq_specProjIdx_mul_right`.
3. **Bessel at rank `r`** (section 3). Under `SimpleSpec A hA r` the top-`r` projector splits
   as the sum of the `r` projectors at the sorted indices `0` to `r - 1`
   (`normSq_specProjTop_eq_sum_of_simpleSpec`), so
   `DelocR.sum_normSq_specProjIdx_le_one` applies `r` times and the sum over an orthonormal
   family is at most `r` (`sum_normSq_specProjTop_le_of_simpleSpec`). This is the one new step
   against `R6.lean`, where `sum_overlap_le_one` gives the bound `1`.
4. **Exchangeability** (section 4). The Householder reflection with axis `x - x'` fixes every
   row of `Y`, so it fixes `Y * O` and moves the Gaussian block only. Its law is invariant
   (`Symmetry.measurePreserving_mul_right`), and the change of variables runs through the
   measurable embedding of `DelocR` section `MeasEmb`, which asks nothing of the integrand.
5. **The count** (section 5). The family is an orthonormal basis of the orthogonal complement
   of the span of the rows of `Y`, which has at least `d - r` elements, so the common value of
   the integral is at most `r / (d - r)`.

`p = 0` needs no exclusion, and `r = 0` is covered (both sides are then `0`).

## Deviation from the task brief

`Defs.isHermitian_transpose_mul_self` asks for a row index of the form `Fin n`, and the row
index here is `Fin r ⊕ Fin p`. This file therefore states the general-index twin
`isHermitian_transpose_mul_self'` and uses it in every signature. No existing file changes.

No `sorry`, no `axiom`.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-- The Gram matrix `Xᵀ X` is Hermitian, for an arbitrary finite row index. The twin of
`Defs.isHermitian_transpose_mul_self`, which fixes the row index to `Fin n`. -/
theorem isHermitian_transpose_mul_self' {m : Type*} [Fintype m] {d : ℕ}
    (X : Matrix m (Fin d) ℝ) : (Xᵀ * X).IsHermitian := by
  simpa using Matrix.isHermitian_conjTranspose_mul_self X

/-! ### 1. The conditional model matrix at rank `r` -/

section YmatR

variable {p d r : ℕ}

/-- The `(r + p) × d` matrix with fixed rows `Y` on top and `t • B` below. Its Gram matrix is
`Yᵀ Y + t² Bᵀ B`. Rank-`r` twin of `R6.Ymat`; the row index is `Fin r ⊕ Fin p`, which makes the
Gram identity one rewrite. -/
noncomputable def YmatR (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) :
    Matrix (Fin r ⊕ Fin p) (Fin d) ℝ :=
  Matrix.fromRows Y (t • B)

@[simp]
theorem YmatR_apply_inl (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ)
    (k : Fin r) : YmatR Y t B (Sum.inl k) = Y k := rfl

@[simp]
theorem YmatR_apply_inr (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ)
    (i : Fin p) : YmatR Y t B (Sum.inr i) = t • B i := rfl

/-- **The Gram matrix of the model.** Rank-`r` twin of `R6.gram_Ymat`, with `Yᵀ Y` in place of
the rank one term `vecMulVec y y`. -/
theorem gram_YmatR (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) :
    (YmatR Y t B)ᵀ * YmatR Y t B = Yᵀ * Y + t ^ 2 • (Bᵀ * B) := by
  rw [YmatR, Matrix.transpose_fromRows, Matrix.fromCols_mul_fromRows, Matrix.transpose_smul,
    Matrix.smul_mul, Matrix.mul_smul, smul_smul, sq]

/-- **Right multiplication distributes over the two blocks.** Rank-`r` twin of `R6.Ymat_mul`,
stated with no hypothesis on `O`; the fixed-rows case is `Y * O = Y`. -/
theorem YmatR_mul (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ)
    (O : Matrix (Fin d) (Fin d) ℝ) : YmatR Y t B * O = YmatR (Y * O) t (B * O) := by
  rw [YmatR, YmatR, Matrix.fromRows_mul, Matrix.smul_mul]

/-- **The model matrix is a measurable function of the Gaussian block.** Rank-`r` twin of
`R6.measurable_Ymat`; task U7 needs it. -/
theorem measurable_YmatR (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ) :
    Measurable (fun B : Matrix (Fin p) (Fin d) ℝ => YmatR Y t B) := by
  refine measurable_pi_lambda _ fun i => ?_
  cases i with
  | inl k =>
      change Measurable fun _ : Matrix (Fin p) (Fin d) ℝ => Y k
      exact measurable_const
  | inr i' =>
      refine measurable_pi_lambda _ fun j => ?_
      change Measurable fun B : Matrix (Fin p) (Fin d) ℝ => t * B i' j
      have h1 : Measurable fun B : Matrix (Fin p) (Fin d) ℝ => B i' := measurable_pi_apply i'
      exact ((measurable_pi_apply j).comp h1).const_mul t

end YmatR

/-! ### 2. The top-`r` projector under an orthogonal conjugation -/

section Conj

variable {d : ℕ} {O : Matrix (Fin d) (Fin d) ℝ}

/-- The top-`r` eigenvalue set is invariant under conjugation by an orthogonal matrix, because
it is defined from `eigenvalues₀`. Twin of `DelocR.eigSetIdx_conj`. -/
theorem topEigSet_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) (r : ℕ) :
    topEigSet (Oᵀ * A * O) hB r = topEigSet A hA r := by
  unfold topEigSet
  rw [eigenvalues₀_conj hO hA hB]

/-- Equal matrices have the same top-`r` projector; the two `IsHermitian` proofs are
irrelevant. Twin of `DelocR.specProjIdx_congr`. -/
theorem specProjTop_congr {A B : Matrix (Fin d) (Fin d) ℝ} (h : A = B) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (r : ℕ) : specProjTop A hA r = specProjTop B hB r := by
  subst h
  rfl

/-- **The top-`r` projector conjugates.** Twin of `DelocR.specProjIdx_conj`. -/
theorem specProjTop_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) (r : ℕ) (x : EuclideanSpace ℝ (Fin d)) :
    specProjTop (Oᵀ * A * O) hB r x =
      rotIso Oᵀ (transpose_orth hO) (specProjTop A hA r (rotIso O hO x)) := by
  rw [specProjTop, specProjTop, topEigSet_conj hO hA hB r, specProj_conj hO]

end Conj

/-- **Right multiplication by an orthogonal matrix moves the test direction.** Twin of
`DelocR.normSq_specProjIdx_mul_right`, with an arbitrary finite row index so that it applies to
`YmatR`. -/
theorem normSq_specProjTop_mul_right {m : Type*} [Fintype m] {d : ℕ} (Z : Matrix m (Fin d) ℝ)
    {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1) (r : ℕ) (w : EuclideanSpace ℝ (Fin d)) :
    ‖specProjTop ((Z * O)ᵀ * (Z * O)) (isHermitian_transpose_mul_self' (Z * O)) r w‖
      = ‖specProjTop (Zᵀ * Z) (isHermitian_transpose_mul_self' Z) r (rotIso O hO w)‖ := by
  have hgram : (Z * O)ᵀ * (Z * O) = Oᵀ * (Zᵀ * Z) * O := by
    simp [Matrix.transpose_mul, Matrix.mul_assoc]
  have hA : (Zᵀ * Z).IsHermitian := isHermitian_transpose_mul_self' Z
  have hB : (Oᵀ * (Zᵀ * Z) * O).IsHermitian := by
    rw [← hgram]
    exact isHermitian_transpose_mul_self' _
  rw [specProjTop_congr hgram (isHermitian_transpose_mul_self' (Z * O)) hB r,
    specProjTop_conj hO hA hB r w, LinearIsometryEquiv.norm_map]

/-! ### 3. Bessel at rank `r`: the top-`r` projector splits into `r` index projectors -/

section Bessel

variable {d : ℕ} {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian} {r : ℕ}

/-- Under `SimpleSpec A hA r` the top-`r` eigenvalue set selects exactly the sorted indices
below `r`. Unlike `mem_topEigSet_iff`, which asks for `TopGap`, simplicity of the top `r`
eigenvalues is enough. -/
theorem mem_topEigSet_of_simpleSpec (hsimple : SimpleSpec A hA r)
    (i : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ i ∈ topEigSet A hA r ↔ (i : ℕ) < r := by
  constructor
  · rintro ⟨k, hk, hEq⟩
    by_cases hik : i = k
    · rw [hik]; exact hk
    · exact absurd hEq.symm (hsimple k i hk (Ne.symm hik))
  · intro hi
    exact ⟨i, hi, rfl⟩

/-- Under `SimpleSpec A hA r` the eigenvalue set at a sorted index `k < r` selects exactly the
index `k`. -/
theorem mem_eigSetIdx_of_simpleSpec (hsimple : SimpleSpec A hA r) {k : ℕ} (hk : k < r)
    (i : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ i ∈ eigSetIdx A hA k ↔ (i : ℕ) = k := by
  constructor
  · rintro ⟨q, hq, hEq⟩
    by_cases hiq : i = q
    · rw [hiq]; exact hq
    · exact absurd hEq.symm (hsimple q i (by rw [hq]; exact hk) (Ne.symm hiq))
  · intro hi
    exact ⟨i, hi, rfl⟩

/-- **The top-`r` projector splits.** Under `SimpleSpec A hA r` the squared norm of the top-`r`
projection is the sum of the squared norms of the `r` projections at the sorted indices `0` to
`r - 1`. Both sides expand in the sorted eigenbasis through `normSq_specProj`. -/
theorem normSq_specProjTop_eq_sum_of_simpleSpec (hsimple : SimpleSpec A hA r) (hr : r ≤ d)
    (x : EuclideanSpace ℝ (Fin d)) :
    ‖specProjTop A hA r x‖ ^ 2 = ∑ k ∈ Finset.range r, ‖specProjIdx A hA k x‖ ^ 2 := by
  classical
  set c : ℕ → ℝ := fun k =>
    if h : k < Fintype.card (Fin d) then
      ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace ⟨k, h⟩, x⟫_ℝ else 0 with hc
  have hcval : ∀ (k : ℕ) (h : k < Fintype.card (Fin d)),
      c k = ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace ⟨k, h⟩, x⟫_ℝ := by
    intro k h
    rw [hc]
    simp only [dif_pos h]
  have hcvalFin : ∀ i : Fin (Fintype.card (Fin d)),
      c (i : ℕ) = ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ := fun i =>
    hcval (i : ℕ) i.isLt
  -- the left side, in the sorted eigenbasis
  have hLHS : ‖specProjTop A hA r x‖ ^ 2
      = ∑ i : Fin (Fintype.card (Fin d)), (if (i : ℕ) < r then c (i : ℕ) else 0) ^ 2 := by
    rw [specProjTop, normSq_specProj A hA]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [hcvalFin i]
    by_cases h : (i : ℕ) < r
    · rw [if_pos ((mem_topEigSet_of_simpleSpec hsimple i).mpr h), if_pos h]
    · rw [if_neg (fun hm => h ((mem_topEigSet_of_simpleSpec hsimple i).mp hm)), if_neg h]
  -- each index term, in the same basis
  have hRHS : ∀ k, k < r → ‖specProjIdx A hA k x‖ ^ 2 = c k ^ 2 := by
    intro k hk
    have hkc : k < Fintype.card (Fin d) := by
      rw [Fintype.card_fin]; omega
    rw [specProjIdx, normSq_specProj A hA,
      Finset.sum_eq_single (⟨k, hkc⟩ : Fin (Fintype.card (Fin d)))]
    · rw [if_pos ((mem_eigSetIdx_of_simpleSpec hsimple hk ⟨k, hkc⟩).mpr rfl), hcval k hkc]
    · intro i _ hik
      rw [if_neg (fun hm => hik (Fin.ext ((mem_eigSetIdx_of_simpleSpec hsimple hk i).mp hm)))]
      norm_num
    · intro h
      exact absurd (Finset.mem_univ _) h
  rw [hLHS, Finset.sum_congr rfl (fun k hk => hRHS k (Finset.mem_range.mp hk))]
  have hsq : ∀ i : Fin (Fintype.card (Fin d)),
      (if (i : ℕ) < r then c (i : ℕ) else 0) ^ 2
        = if (i : ℕ) < r then c (i : ℕ) ^ 2 else 0 := by
    intro i
    by_cases h : (i : ℕ) < r <;> simp [h]
  rw [Finset.sum_congr rfl (fun i _ => hsq i),
    Fin.sum_univ_eq_sum_range (fun k => if k < r then c k ^ 2 else 0) (Fintype.card (Fin d)),
    ← Finset.sum_filter]
  congr 1
  ext k
  simp only [Finset.mem_filter, Finset.mem_range, Fintype.card_fin]
  omega

/-- **The Bessel bound at rank `r`.** Under `SimpleSpec A hA r` the squared norms of the top-`r`
projection of an orthonormal family sum to at most `r`. Rank-`r` twin of
`Symmetry.sum_overlap_le_one`, whose bound is `1`. -/
theorem sum_normSq_specProjTop_le_of_simpleSpec {ι : Type*} [Fintype ι]
    (hsimple : SimpleSpec A hA r) (hr : r ≤ d)
    {e : ι → EuclideanSpace ℝ (Fin d)} (he : Orthonormal ℝ e) :
    ∑ j, ‖specProjTop A hA r (e j)‖ ^ 2 ≤ (r : ℝ) := by
  calc ∑ j, ‖specProjTop A hA r (e j)‖ ^ 2
      = ∑ j, ∑ k ∈ Finset.range r, ‖specProjIdx A hA k (e j)‖ ^ 2 :=
        Finset.sum_congr rfl fun j _ =>
          normSq_specProjTop_eq_sum_of_simpleSpec hsimple hr (e j)
    _ = ∑ k ∈ Finset.range r, ∑ j, ‖specProjIdx A hA k (e j)‖ ^ 2 := Finset.sum_comm
    _ ≤ ∑ _k ∈ Finset.range r, (1 : ℝ) := by
        refine Finset.sum_le_sum fun k hk => ?_
        have hkr : k < r := Finset.mem_range.mp hk
        exact sum_normSq_specProjIdx_le_one (simpleSpec_mono (by omega) hsimple) (by omega) he
    _ = (r : ℝ) := by simp

end Bessel

/-! ### 4. Exchangeability of the test direction -/

/-- **Exchangeability.** Conditionally on `Y` the mean squared top-`r` overlap is the same for
every unit direction orthogonal to every row of `Y`: the Householder reflection with axis
`x - x'` fixes `Y` and `B ↦ B O` preserves the Gaussian law. Rank-`r` twin of
`R6.lintegral_overlap_Ymat_eq`; the change of variables runs through the measurable embedding of
`DelocR`, so the integrand needs no measurability. -/
theorem lintegral_normSq_specProjTop_YmatR_eq {p d r : ℕ} (Y : Matrix (Fin r) (Fin d) ℝ)
    (t : ℝ) {x x' : EuclideanSpace ℝ (Fin d)} (hx : ‖x‖ = 1) (hx' : ‖x'‖ = 1)
    (hxY : ∀ k, Y k ⬝ᵥ WithLp.ofLp x = 0) (hx'Y : ∀ k, Y k ⬝ᵥ WithLp.ofLp x' = 0) :
    ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d)
      = ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r x'‖ ^ 2) ∂(gaussianMatrix p d) := by
  by_cases hxx : x = x'
  · rw [hxx]
  have hdot : ∀ u v : EuclideanSpace ℝ (Fin d),
      WithLp.ofLp u ⬝ᵥ WithLp.ofLp v = ⟪u, v⟫_ℝ := fun u v =>
    (inner_euclidean_eq_dotProduct u v).symm
  have hx1 : WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hx]
    norm_num
  have hx'1 : WithLp.ofLp x' ⬝ᵥ WithLp.ofLp x' = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hx']
    norm_num
  have hne : WithLp.ofLp x ≠ WithLp.ofLp x' := fun h => hxx (by
    have := congrArg (WithLp.toLp 2) h
    simpa using this)
  set a : Fin d → ℝ := WithLp.ofLp x - WithLp.ofLp x' with ha
  set O : Matrix (Fin d) (Fin d) ℝ := householder a with hOdef
  have ha0 : a ⬝ᵥ a ≠ 0 := fun h => (sub_ne_zero.mpr hne) (dotProduct_self_eq_zero.mp h)
  have hO : Oᵀ * O = 1 := householder_orth ha0
  have haY : ∀ k, a ⬝ᵥ Y k = 0 := by
    intro k
    have h1 : WithLp.ofLp x ⬝ᵥ Y k = 0 := by rw [dotProduct_comm]; exact hxY k
    have h2 : WithLp.ofLp x' ⬝ᵥ Y k = 0 := by rw [dotProduct_comm]; exact hx'Y k
    rw [ha, sub_dotProduct, h1, h2, sub_zero]
  have hrow : ∀ k, Oᵀ *ᵥ Y k = Y k := by
    intro k
    rw [hOdef, householder_transpose]
    exact householder_apply_of_orth (haY k)
  have hYO : Y * O = Y := by
    ext k j
    have h := congrFun (hrow k) j
    change ∑ l, Y k l * O l j = Y k j
    rw [← h]
    simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
    exact Finset.sum_congr rfl fun l _ => mul_comm _ _
  have hOx : rotIso O hO x = x' := by
    rw [rotIso_apply]
    exact congrArg (WithLp.toLp 2) (householder_sub_apply hx1 hx'1 hne)
  have key : ∀ B : Matrix (Fin p) (Fin d) ℝ,
      ‖specProjTop ((YmatR Y t (B * O))ᵀ * YmatR Y t (B * O))
          (isHermitian_transpose_mul_self' (YmatR Y t (B * O))) r x‖
        = ‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
          (isHermitian_transpose_mul_self' (YmatR Y t B)) r x'‖ := by
    intro B
    have hmul : YmatR Y t B * O = YmatR Y t (B * O) := by
      rw [YmatR_mul, hYO]
    rw [← hmul, normSq_specProjTop_mul_right (YmatR Y t B) hO r x, hOx]
  calc ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
          (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d)
      = ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t (B * O))ᵀ * YmatR Y t (B * O))
          (isHermitian_transpose_mul_self' (YmatR Y t (B * O))) r x‖ ^ 2)
          ∂(gaussianMatrix p d) :=
        ((measurePreserving_mul_right hO p).lintegral_comp_emb
          (measurableEmbedding_mul_right hO p) _).symm
    _ = _ := by simp only [key]

/-! ### 5. The uniform bound `r / (d - r)` -/

/-- **Task U6, the unit case.** Bessel at rank `r` over an orthonormal basis of the orthogonal
complement of the row space of `Y`, which has at least `d - r` elements, together with
exchangeability. Rank-`r` twin of `R6.lintegral_overlap_Ymat_le`, whose constant is
`1 / (d - 1)`. -/
theorem lintegral_normSq_specProjTop_YmatR_le_of_unit {p d r : ℕ} (hd : r < d)
    (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ)
    (hsimple : ∀ᵐ B ∂(gaussianMatrix p d),
      SimpleSpec ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r)
    {x : EuclideanSpace ℝ (Fin d)} (hx : ‖x‖ = 1) (hxY : ∀ k, Y k ⬝ᵥ WithLp.ofLp x = 0) :
    ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal ((r : ℝ) / ((d : ℝ) - r)) := by
  classical
  set K₀ : Submodule ℝ (EuclideanSpace ℝ (Fin d)) :=
    Submodule.span ℝ (Set.range fun k : Fin r =>
      (WithLp.toLp 2 (Y k) : EuclideanSpace ℝ (Fin d))) with hK₀
  set K : Submodule ℝ (EuclideanSpace ℝ (Fin d)) := K₀ᗮ with hK
  have hrk : d - r ≤ Module.finrank ℝ K := by
    have h2 : Module.finrank ℝ K₀ + Module.finrank ℝ K
        = Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) := by
      rw [hK]
      exact Submodule.finrank_add_finrank_orthogonal (𝕜 := ℝ) (K := K₀)
    have h1 : Module.finrank ℝ K₀ ≤ r := by
      have h3 := finrank_range_le_card (R := ℝ)
        (fun k : Fin r => (WithLp.toLp 2 (Y k) : EuclideanSpace ℝ (Fin d)))
      simpa [Set.finrank, hK₀] using h3
    rw [finrank_euclideanSpace_fin] at h2
    omega
  set b := stdOrthonormalBasis ℝ K with hb
  set e : Fin (Module.finrank ℝ K) → EuclideanSpace ℝ (Fin d) :=
    fun j => (b j : EuclideanSpace ℝ (Fin d)) with he
  have hon : Orthonormal ℝ e :=
    (K.subtypeₗᵢ.orthonormal_comp_iff (v := fun j => b j)).mpr b.orthonormal
  have henorm : ∀ j, ‖e j‖ = 1 := fun j => hon.1 j
  have heY : ∀ j, ∀ k, Y k ⬝ᵥ WithLp.ofLp (e j) = 0 := by
    intro j k
    have hmem : (b j : EuclideanSpace ℝ (Fin d)) ∈ K₀ᗮ := (b j).2
    have hsp : (WithLp.toLp 2 (Y k) : EuclideanSpace ℝ (Fin d)) ∈ K₀ :=
      Submodule.subset_span (Set.mem_range_self k)
    have hzero := (Submodule.mem_orthogonal K₀ _).mp hmem _ hsp
    rw [inner_euclidean_eq_dotProduct] at hzero
    exact hzero
  set I : ℝ≥0∞ := ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
      (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d) with hI
  have hsame : ∀ j, ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
      (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
      ∂(gaussianMatrix p d) = I := fun j =>
    (lintegral_normSq_specProjTop_YmatR_eq Y t hx (henorm j) hxY (heY j)).symm
  have hae : ∀ᵐ B ∂(gaussianMatrix p d),
      ∑ j, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
        ≤ ENNReal.ofReal (r : ℝ) := by
    filter_upwards [hsimple] with B hB
    have hbes := sum_normSq_specProjTop_le_of_simpleSpec hB hd.le hon
    calc ∑ j, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
        = ENNReal.ofReal (∑ j, ‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2) :=
          (ENNReal.ofReal_sum_of_nonneg fun j _ => by positivity).symm
      _ ≤ ENNReal.ofReal (r : ℝ) := ENNReal.ofReal_le_ofReal hbes
  have hsum : ∑ j, ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
      (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
      ∂(gaussianMatrix p d) ≤ ENNReal.ofReal (r : ℝ) := by
    calc ∑ j, ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
            ∂(gaussianMatrix p d)
        ≤ ∫⁻ B, ∑ j, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r (e j)‖ ^ 2)
            ∂(gaussianMatrix p d) := sum_lintegral_le _ _
      _ ≤ ∫⁻ _, ENNReal.ofReal (r : ℝ) ∂(gaussianMatrix p d) := lintegral_mono_ae hae
      _ = ENNReal.ofReal (r : ℝ) := by simp
  simp only [hsame, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
  have hpos : (0 : ℝ) < (d : ℝ) - r := by
    have : (r : ℝ) < (d : ℝ) := by exact_mod_cast hd
    linarith
  have hcard : ENNReal.ofReal ((d : ℝ) - r) ≤ (Module.finrank ℝ K : ℝ≥0∞) := by
    have h1 : (d : ℝ) - r ≤ (Module.finrank ℝ K : ℝ) := by
      have hle : ((d - r : ℕ) : ℝ) ≤ (Module.finrank ℝ K : ℝ) := by exact_mod_cast hrk
      have hcast : ((d - r : ℕ) : ℝ) = (d : ℝ) - r := by
        have hdr : r ≤ d := hd.le
        push_cast [Nat.cast_sub hdr]
        ring
      linarith [hcast ▸ hle]
    calc ENNReal.ofReal ((d : ℝ) - r)
        ≤ ENNReal.ofReal ((Module.finrank ℝ K : ℝ)) := ENNReal.ofReal_le_ofReal h1
      _ = _ := by rw [ENNReal.ofReal_natCast]
  have hstep : ENNReal.ofReal ((d : ℝ) - r) * I ≤ ENNReal.ofReal (r : ℝ) :=
    le_trans (by gcongr) hsum
  have hne0 : ENNReal.ofReal ((d : ℝ) - r) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    exact hpos
  have hdiv : I ≤ ENNReal.ofReal (r : ℝ) / ENNReal.ofReal ((d : ℝ) - r) := by
    rw [ENNReal.le_div_iff_mul_le (Or.inl hne0) (Or.inl ENNReal.ofReal_ne_top), mul_comm]
    exact hstep
  exact hdiv.trans (le_of_eq (ENNReal.ofReal_div_of_pos hpos).symm)

/-- **Task U6.** The same bound for a direction of norm at most one, by homogeneity of the
projector. Rank-`r` twin of `R6.lintegral_overlap_Ymat_le'`. -/
theorem lintegral_normSq_specProjTop_YmatR_le {p d r : ℕ} (hd : r < d)
    (Y : Matrix (Fin r) (Fin d) ℝ) (t : ℝ)
    (hsimple : ∀ᵐ B ∂(gaussianMatrix p d),
      SimpleSpec ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r)
    {x : EuclideanSpace ℝ (Fin d)} (hx : ‖x‖ ≤ 1) (hxY : ∀ k, Y k ⬝ᵥ WithLp.ofLp x = 0) :
    ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal ((r : ℝ) / ((d : ℝ) - r)) := by
  rcases eq_or_ne x 0 with h0 | h0
  · subst h0
    simp only [map_zero, norm_zero]
    norm_num
  · have hn0 : ‖x‖ ≠ 0 := norm_ne_zero_iff.mpr h0
    set u : EuclideanSpace ℝ (Fin d) := ‖x‖⁻¹ • x with hu
    have hun : ‖u‖ = 1 := by
      rw [hu, norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hn0]
    have huY : ∀ k, Y k ⬝ᵥ WithLp.ofLp u = 0 := by
      intro k
      change Y k ⬝ᵥ (‖x‖⁻¹ • WithLp.ofLp x) = 0
      rw [dotProduct_smul, hxY k, smul_eq_mul, mul_zero]
    have hxu : x = ‖x‖ • u := by
      rw [hu, smul_smul, mul_inv_cancel₀ hn0, one_smul]
    have hmono : ∀ B : Matrix (Fin p) (Fin d) ℝ,
        ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2)
          ≤ ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r u‖ ^ 2) := by
      intro B
      refine ENNReal.ofReal_le_ofReal ?_
      conv_lhs => rw [hxu]
      rw [map_smul, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]
      have h1 : ‖x‖ ^ 2 ≤ 1 := by nlinarith [norm_nonneg x]
      nlinarith [sq_nonneg ‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
        (isHermitian_transpose_mul_self' (YmatR Y t B)) r u‖]
    calc ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r x‖ ^ 2) ∂(gaussianMatrix p d)
        ≤ ∫⁻ B, ENNReal.ofReal (‖specProjTop ((YmatR Y t B)ᵀ * YmatR Y t B)
            (isHermitian_transpose_mul_self' (YmatR Y t B)) r u‖ ^ 2) ∂(gaussianMatrix p d) :=
          lintegral_mono hmono
      _ ≤ ENNReal.ofReal ((r : ℝ) / ((d : ℝ) - r)) :=
          lintegral_normSq_specProjTop_YmatR_le_of_unit hd Y t hsimple hun huY

end StackedSVD
