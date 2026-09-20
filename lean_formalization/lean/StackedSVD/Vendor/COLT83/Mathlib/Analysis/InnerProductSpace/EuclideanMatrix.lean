/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.CStarAlgebra.Matrix
public import Mathlib.Analysis.InnerProductSpace.Adjoint
public import Mathlib.LinearAlgebra.Matrix.PosDef
public import Mathlib.Topology.Instances.Matrix

set_option autoImplicit false

/-!
# Matrices as linear maps on Euclidean spaces

Small lemmas about `Matrix.toEuclideanLin` and `Matrix.toEuclideanCLM` for real matrices:
coordinates of `toEuclideanLin M v`, products, adjoints, inverses, spans of images, and the
continuity of the quadratic form `x ↦ xᵀ M x`.
-/

@[expose] public section

open scoped RealInnerProductSpace

namespace Matrix

lemma ofLp_toEuclideanLin {ι κ : Type*} [Fintype κ] [DecidableEq κ] (M : Matrix ι κ ℝ)
    (v : EuclideanSpace ℝ κ) :
    WithLp.ofLp (Matrix.toEuclideanLin M v) = M *ᵥ WithLp.ofLp v := rfl

lemma toEuclideanLin_mul_apply {ι : Type*} [Fintype ι] [DecidableEq ι] (A B : Matrix ι ι ℝ)
    (v : EuclideanSpace ℝ ι) :
    Matrix.toEuclideanLin (A * B) v = Matrix.toEuclideanLin A (Matrix.toEuclideanLin B v) := by
  apply WithLp.ofLp_injective
  change (A * B) *ᵥ WithLp.ofLp v = A *ᵥ (B *ᵥ WithLp.ofLp v)
  rw [Matrix.mulVec_mulVec]

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {𝒳 : Set (EuclideanSpace ℝ ι)}

/-- The adjoint of the continuous linear map associated to a real matrix `M` is the map
associated to `Mᵀ`. -/
lemma adjoint_toEuclideanCLM (M : Matrix ι ι ℝ) :
    ContinuousLinearMap.adjoint (toEuclideanCLM (𝕜 := ℝ) M) = toEuclideanCLM (𝕜 := ℝ) Mᵀ := by
  rw [← ContinuousLinearMap.star_eq_adjoint, ← map_star, star_eq_conjTranspose,
    conjTranspose_eq_transpose_of_trivial]

lemma toEuclideanCLM_inv_apply_toEuclideanCLM {S : Matrix ι ι ℝ} (hS : S.PosDef)
    (θ : EuclideanSpace ℝ ι) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) S⁻¹ (Matrix.toEuclideanCLM (𝕜 := ℝ) S θ) = θ := by
  rw [← mul_apply_eq_comp, ← map_mul,
    Matrix.nonsing_inv_mul _ ((Matrix.isUnit_iff_isUnit_det _).1 hS.isUnit), map_one,
    one_apply_eq_self]

/-- The image of a spanning set by an invertible linear map is spanning. -/
lemma span_image_toEuclideanCLM_eq_top {M : Matrix ι ι ℝ} (hM : IsUnit M)
    (hspan : Submodule.span ℝ 𝒳 = ⊤) :
    Submodule.span ℝ (Matrix.toEuclideanCLM (𝕜 := ℝ) M '' 𝒳) = ⊤ := by
  have hdet : IsUnit M.det := (Matrix.isUnit_iff_isUnit_det M).1 hM
  have h1 : Matrix.toEuclideanCLM (𝕜 := ℝ) M * Matrix.toEuclideanCLM (𝕜 := ℝ) M⁻¹ = 1 := by
    rw [← map_mul, Matrix.mul_nonsing_inv _ hdet, map_one]
  have hsurj : Function.Surjective (Matrix.toEuclideanCLM (𝕜 := ℝ) M) := fun y ↦
    ⟨Matrix.toEuclideanCLM (𝕜 := ℝ) M⁻¹ y, by
      simpa using congr_fun (congrArg DFunLike.coe h1) y⟩
  rw [← ContinuousLinearMap.coe_coe, Submodule.span_image, hspan, Submodule.map_top,
    LinearMap.range_eq_top]
  exact hsurj

omit [DecidableEq ι] in
lemma continuous_dotProduct_mulVec (M : Matrix ι ι ℝ) :
    Continuous fun x : EuclideanSpace ℝ ι ↦ WithLp.ofLp x ⬝ᵥ M *ᵥ WithLp.ofLp x :=
  (PiLp.continuous_ofLp _ _).dotProduct
    (continuous_const.matrix_mulVec (PiLp.continuous_ofLp _ _))

omit [DecidableEq ι] in
lemma bddAbove_range_dotProduct_mulVec (h𝒳 : IsCompact 𝒳) (M : Matrix ι ι ℝ) :
    BddAbove (Set.range fun x : 𝒳 ↦
      WithLp.ofLp (x : EuclideanSpace ℝ ι) ⬝ᵥ M *ᵥ WithLp.ofLp (x : EuclideanSpace ℝ ι)) := by
  refine (h𝒳.image (continuous_dotProduct_mulVec M)).bddAbove.mono ?_
  rintro _ ⟨x, rfl⟩
  exact ⟨x, x.2, rfl⟩

end Matrix
