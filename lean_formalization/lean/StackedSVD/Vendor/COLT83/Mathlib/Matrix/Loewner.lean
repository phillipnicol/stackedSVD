/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Matrix.LoewnerInv
public import Mathlib.Analysis.CStarAlgebra.Matrix
public import Mathlib.LinearAlgebra.Matrix.SchurComplement

set_option autoImplicit false

/-!
# The Loewner order, quadratic forms and square roots

Elementary facts about the Loewner order `A ≤ B` on matrices (`Matrix.le_iff : A ≤ B ↔
(B - A).PosSemidef`, scoped in `MatrixOrder`) that are stated in terms of quadratic forms:

* `A ≤ B` implies `x* A x ≤ x* B x` for every `x`, and conversely for Hermitian matrices
  (`dotProduct_mulVec_le_of_le`, `le_iff_dotProduct_mulVec_le`);
* the order is preserved by congruence `A ↦ M A M*` (`mul_mul_conjTranspose_le_of_le`);
* a matrix above a positive definite matrix is positive definite (`PosDef.of_le`);
* the scalar form of the antitonicity of inversion: `c • A ≤ B` with `A ≻ 0`, `c > 0` gives
  `B⁻¹ ≤ c⁻¹ • A⁻¹`, hence `x* B⁻¹ x ≤ c⁻¹ x* A⁻¹ x` (`PosDef.inv_le_inv_smul`,
  `PosDef.dotProduct_inv_mulVec_le_of_smul_le`).

We also record the matrix determinant lemma in `vecMulVec` form (`det_add_vecMulVec`) and facts
about the positive square root `CFC.sqrt A` of a positive definite matrix: `CFC.sqrt A⁻¹ * A *
CFC.sqrt A⁻¹ = 1` (`PosDef.sqrt_inv_mul_mul_sqrt_inv`), and for the associated linear maps on
`EuclideanSpace ℝ n`, `‖√A x‖ ^ 2 = x ⬝ᵥ A *ᵥ x` (`norm_toEuclideanCLM_sqrt_sq`) and
`⟪√A u, √A⁻¹ v⟫ = ⟪u, v⟫` (`inner_toEuclideanCLM_sqrt_sqrt_inv`).
-/

@[expose] public section

open scoped MatrixOrder ComplexOrder RealInnerProductSpace

namespace Matrix

section det

variable {R n : Type*} [CommRing R] [Fintype n] [DecidableEq n]

/-- The **matrix determinant lemma** for a rank-one update `A + u vᵀ`. -/
lemma det_add_vecMulVec {A : Matrix n n R} (hA : IsUnit A.det) (u v : n → R) :
    (A + vecMulVec u v).det = A.det * (1 + v ⬝ᵥ A⁻¹ *ᵥ u) := by
  have h : A + vecMulVec u v = A * (1 + vecMulVec (A⁻¹ *ᵥ u) v) := by
    rw [Matrix.mul_add, Matrix.mul_one, mul_vecMulVec, mulVec_mulVec, mul_nonsing_inv A hA,
      one_mulVec]
  rw [h, det_mul, vecMulVec_eq Unit, det_one_add_replicateCol_mul_replicateRow]

/-- The **Sherman–Morrison formula**: the inverse of a rank-one update `A + u vᵀ` of an
invertible matrix `A` with `1 + vᵀ A⁻¹ u ≠ 0`. -/
lemma inv_add_vecMulVec {K : Type*} [Field K] {A : Matrix n n K} (hA : IsUnit A.det)
    (u v : n → K) (h : 1 + v ⬝ᵥ A⁻¹ *ᵥ u ≠ 0) :
    (A + vecMulVec u v)⁻¹ =
      A⁻¹ - (1 + v ⬝ᵥ A⁻¹ *ᵥ u)⁻¹ • vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹) := by
  set s := v ⬝ᵥ A⁻¹ *ᵥ u with hs
  refine inv_eq_right_inv ?_
  have h1 : A * vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹) = vecMulVec u (v ᵥ* A⁻¹) := by
    rw [mul_vecMulVec, mulVec_mulVec, mul_nonsing_inv A hA, one_mulVec]
  have h2 : vecMulVec u v * vecMulVec (A⁻¹ *ᵥ u) (v ᵥ* A⁻¹) = s • vecMulVec u (v ᵥ* A⁻¹) := by
    rw [vecMulVec_mul_vecMulVec, vecMulVec_smul]
  have h3 : vecMulVec u v * A⁻¹ = vecMulVec u (v ᵥ* A⁻¹) := vecMulVec_mul _ _ _
  rw [Matrix.add_mul, Matrix.mul_sub, Matrix.mul_sub, mul_nonsing_inv A hA, Matrix.mul_smul,
    Matrix.mul_smul, h1, h2, h3, smul_smul]
  have h4 : (1 + s)⁻¹ * s = 1 - (1 + s)⁻¹ := by field_simp; ring
  rw [h4, sub_smul, one_smul]
  abel

/-- The trace form of the Sherman–Morrison formula for a symmetric update `A + u uᵀ` of an
invertible symmetric matrix `A`. -/
lemma trace_inv_add_vecMulVec_self {K : Type*} [Field K] {A : Matrix n n K} (hA : IsUnit A.det)
    (hAt : Aᵀ = A) (u : n → K) (h : 1 + u ⬝ᵥ A⁻¹ *ᵥ u ≠ 0) :
    (A + vecMulVec u u)⁻¹.trace = A⁻¹.trace - (u ⬝ᵥ (A⁻¹ * A⁻¹) *ᵥ u) / (1 + u ⬝ᵥ A⁻¹ *ᵥ u) := by
  have key : (A⁻¹ *ᵥ u) ⬝ᵥ (u ᵥ* A⁻¹) = u ⬝ᵥ (A⁻¹ * A⁻¹) *ᵥ u := by
    rw [← mulVec_transpose, transpose_nonsing_inv, hAt, dotProduct_mulVec, ← mulVec_transpose,
      transpose_nonsing_inv, hAt, mulVec_mulVec, dotProduct_comm]
  rw [inv_add_vecMulVec hA u u h, trace_sub, trace_smul, trace_vecMulVec, smul_eq_mul,
    div_eq_inv_mul, key]

end det

section order

variable {𝕜 m n : Type*} [RCLike 𝕜] [Fintype n] {A B : Matrix n n 𝕜}

/-- The Loewner order implies the order of the quadratic forms. -/
lemma dotProduct_mulVec_le_of_le (h : A ≤ B) (x : n → 𝕜) :
    star x ⬝ᵥ A *ᵥ x ≤ star x ⬝ᵥ B *ᵥ x := by
  have := (le_iff.mp h).dotProduct_mulVec_nonneg x
  rwa [sub_mulVec, dotProduct_sub, sub_nonneg] at this

/-- For Hermitian matrices, the Loewner order is the order of the quadratic forms. -/
lemma le_iff_dotProduct_mulVec_le (hA : A.IsHermitian) (hB : B.IsHermitian) :
    A ≤ B ↔ ∀ x, star x ⬝ᵥ A *ᵥ x ≤ star x ⬝ᵥ B *ᵥ x := by
  refine ⟨dotProduct_mulVec_le_of_le, fun h ↦ le_iff.mpr ?_⟩
  refine PosSemidef.of_dotProduct_mulVec_nonneg (hB.sub hA) fun x ↦ ?_
  rw [sub_mulVec, dotProduct_sub, sub_nonneg]
  exact h x

/-- Congruence `A ↦ M A M*` preserves the Loewner order. -/
lemma mul_mul_conjTranspose_le_of_le [Finite m] (h : A ≤ B) (M : Matrix m n 𝕜) :
    M * A * Mᴴ ≤ M * B * Mᴴ := by
  have := Fintype.ofFinite m
  rw [le_iff] at h ⊢
  simpa [Matrix.mul_sub, Matrix.sub_mul] using h.mul_mul_conjTranspose_same M

/-- Congruence `A ↦ M* A M` preserves the Loewner order. -/
lemma conjTranspose_mul_mul_le_of_le [Finite m] (h : A ≤ B) (M : Matrix n m 𝕜) :
    Mᴴ * A * M ≤ Mᴴ * B * M := by
  have := Fintype.ofFinite m
  rw [le_iff] at h ⊢
  simpa [Matrix.mul_sub, Matrix.sub_mul] using h.conjTranspose_mul_mul_same M

omit [Fintype n] in
/-- A matrix above a positive definite matrix in the Loewner order is positive definite. -/
lemma PosDef.of_le (hA : A.PosDef) (hAB : A ≤ B) : B.PosDef := by
  simpa using hA.add_posSemidef (le_iff.mp hAB)

variable [DecidableEq n]

/-- If `c • A ≤ B` with `A` positive definite and `c > 0`, then `B⁻¹ ≤ c⁻¹ • A⁻¹`. -/
lemma PosDef.inv_le_inv_smul (hA : A.PosDef) {c : 𝕜} (hc : 0 < c) (h : c • A ≤ B) :
    B⁻¹ ≤ c⁻¹ • A⁻¹ := by
  have hA' : IsUnit A.det := (isUnit_iff_isUnit_det A).mp hA.isUnit
  have h_inv : (c • A)⁻¹ = c⁻¹ • A⁻¹ := by
    refine inv_eq_left_inv ?_
    rw [smul_mul_smul_comm, nonsing_inv_mul A hA', inv_mul_cancel₀ hc.ne', one_smul]
  rw [← h_inv]
  exact (hA.smul hc).inv_le_inv h

/-- If `c • A ≤ B` with `A` positive definite and `c > 0`, then
`x* B⁻¹ x ≤ c⁻¹ * x* A⁻¹ x` for every `x`. -/
lemma PosDef.dotProduct_inv_mulVec_le_of_smul_le (hA : A.PosDef) {c : 𝕜} (hc : 0 < c)
    (h : c • A ≤ B) (x : n → 𝕜) :
    star x ⬝ᵥ B⁻¹ *ᵥ x ≤ c⁻¹ * (star x ⬝ᵥ A⁻¹ *ᵥ x) := by
  have := dotProduct_mulVec_le_of_le (hA.inv_le_inv_smul hc h) x
  rwa [smul_mulVec, dotProduct_smul, smul_eq_mul] at this

end order

section sqrt

variable {𝕜 n : Type*} [RCLike 𝕜] [Fintype n] [DecidableEq n] {A : Matrix n n 𝕜}

-- The `DecidableEq n` instance enters the type through the continuous functional calculus
-- instance on matrices; the linter does not see it.
set_option linter.unusedDecidableInType false in
/-- The square root of a positive definite matrix is positive definite. -/
lemma PosDef.sqrt (hA : A.PosDef) : (CFC.sqrt A).PosDef := by
  have hnn : (CFC.sqrt A).PosSemidef := nonneg_iff_posSemidef.mp (CFC.sqrt_nonneg A)
  refine hnn.posDef_iff_isUnit.mpr ?_
  rw [CFC.isUnit_sqrt_iff A hA.posSemidef.nonneg]
  exact hA.isUnit

set_option linter.unusedDecidableInType false in
/-- `√A * √A = A` for a positive semidefinite matrix. -/
lemma PosSemidef.sqrt_mul_sqrt (hA : A.PosSemidef) : CFC.sqrt A * CFC.sqrt A = A :=
  CFC.sqrt_mul_sqrt_self A hA.nonneg

/-- `√A * √(A⁻¹) = 1` for a positive definite matrix. -/
lemma PosDef.sqrt_mul_sqrt_inv (hA : A.PosDef) : CFC.sqrt A * CFC.sqrt A⁻¹ = 1 := by
  rw [← hA.posSemidef.inv_sqrt, mul_nonsing_inv _ ((isUnit_iff_isUnit_det _).mp hA.sqrt.isUnit)]

/-- `√(A⁻¹) * √A = 1` for a positive definite matrix. -/
lemma PosDef.sqrt_inv_mul_sqrt (hA : A.PosDef) : CFC.sqrt A⁻¹ * CFC.sqrt A = 1 := by
  rw [← hA.posSemidef.inv_sqrt, nonsing_inv_mul _ ((isUnit_iff_isUnit_det _).mp hA.sqrt.isUnit)]

/-- `√(A⁻¹) * A * √(A⁻¹) = 1` for a positive definite matrix: `√(A⁻¹)` whitens `A`. -/
lemma PosDef.sqrt_inv_mul_mul_sqrt_inv (hA : A.PosDef) : CFC.sqrt A⁻¹ * A * CFC.sqrt A⁻¹ = 1 := by
  calc CFC.sqrt A⁻¹ * A * CFC.sqrt A⁻¹
      = CFC.sqrt A⁻¹ * (CFC.sqrt A * CFC.sqrt A) * CFC.sqrt A⁻¹ := by
        rw [hA.posSemidef.sqrt_mul_sqrt]
    _ = 1 := by
        rw [← mul_assoc, mul_assoc _ (CFC.sqrt A), hA.sqrt_mul_sqrt_inv, mul_one,
          hA.sqrt_inv_mul_sqrt]

/-- `√A * A⁻¹ * √A = 1` for a positive definite matrix. -/
lemma PosDef.sqrt_mul_inv_mul_sqrt (hA : A.PosDef) : CFC.sqrt A * A⁻¹ * CFC.sqrt A = 1 := by
  calc CFC.sqrt A * A⁻¹ * CFC.sqrt A
      = CFC.sqrt A * (CFC.sqrt A⁻¹ * CFC.sqrt A⁻¹) * CFC.sqrt A := by
        rw [hA.inv.posSemidef.sqrt_mul_sqrt]
    _ = 1 := by
        rw [← mul_assoc, mul_assoc _ (CFC.sqrt A⁻¹), hA.sqrt_inv_mul_sqrt, mul_one,
          hA.sqrt_mul_sqrt_inv]

set_option linter.unusedDecidableInType false in
/-- The square root of a real matrix is symmetric (`CFC.sqrt A` is positive semidefinite for
every `A`, being `0` when `A` is not positive semidefinite). -/
lemma transpose_sqrt (A : Matrix n n ℝ) : (CFC.sqrt A)ᵀ = CFC.sqrt A := by
  have h : (CFC.sqrt A).PosSemidef := nonneg_iff_posSemidef.1 (CFC.sqrt_nonneg _)
  simpa using h.isHermitian.eq

end sqrt

section euclidean

variable {n : Type*} [Fintype n] [DecidableEq n] {A : Matrix n n ℝ}

/-- `‖√A x‖ ^ 2 = xᵀ A x` for a positive semidefinite matrix `A`. -/
lemma norm_toEuclideanCLM_sqrt_sq (hA : A.PosSemidef) (x : EuclideanSpace ℝ n) :
    ‖toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt A) x‖ ^ 2 = WithLp.ofLp x ⬝ᵥ A *ᵥ WithLp.ofLp x := by
  rw [← real_inner_self_eq_norm_sq, EuclideanSpace.inner_eq_star_dotProduct, star_trivial,
    ofLp_toEuclideanCLM, dotProduct_mulVec, ← mulVec_transpose, transpose_sqrt, mulVec_mulVec,
    hA.sqrt_mul_sqrt, dotProduct_comm]

/-- `‖√(A⁻¹) x‖ ^ 2 = xᵀ A⁻¹ x` for a positive definite matrix `A`. -/
lemma norm_toEuclideanCLM_sqrt_inv_sq (hA : A.PosDef) (x : EuclideanSpace ℝ n) :
    ‖toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt A⁻¹) x‖ ^ 2 = WithLp.ofLp x ⬝ᵥ A⁻¹ *ᵥ WithLp.ofLp x :=
  norm_toEuclideanCLM_sqrt_sq hA.inv.posSemidef x

/-- `⟪√A u, √(A⁻¹) v⟫ = ⟪u, v⟫` for a positive definite matrix `A`. -/
lemma inner_toEuclideanCLM_sqrt_sqrt_inv (hA : A.PosDef) (u v : EuclideanSpace ℝ n) :
    ⟪toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt A) u, toEuclideanCLM (𝕜 := ℝ) (CFC.sqrt A⁻¹) v⟫ =
      ⟪u, v⟫ := by
  rw [EuclideanSpace.inner_eq_star_dotProduct, EuclideanSpace.inner_eq_star_dotProduct,
    star_trivial, star_trivial, ofLp_toEuclideanCLM, ofLp_toEuclideanCLM, dotProduct_mulVec,
    ← mulVec_transpose, transpose_sqrt, mulVec_mulVec, hA.sqrt_mul_sqrt_inv, one_mulVec,
    dotProduct_comm]

end euclidean

section convex

variable {𝕜 n : Type*} [RCLike 𝕜]

/-- The cone of positive semidefinite matrices is convex. -/
lemma convex_posSemidef : Convex ℝ {A : Matrix n n 𝕜 | A.PosSemidef} := by
  intro A hA B hB a b ha hb _
  exact (hA.smul ha).add (hB.smul hb)

end convex

end Matrix
