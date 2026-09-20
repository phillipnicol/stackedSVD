/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.Matrix.Order

set_option autoImplicit false

/-!
# Inversion is antitone for the Loewner order

For positive definite matrices `A ≤ B` (Loewner order), we have `B⁻¹ ≤ A⁻¹`.

Mathlib proves this for units of a unital C⋆-algebra (`CStarAlgebra.inv_le_inv`), which applies
to `Matrix n n ℂ` but not to `Matrix n n ℝ` (there is no `CStarAlgebra` instance on real
matrices). Here we give a direct proof over any `RCLike 𝕜`, using only the real continuous
functional calculus on matrices: writing `S := √(A⁻¹)` (so that `S * A * S = 1`), the inequality
`A ≤ B` becomes `1 ≤ S * B * S`; a positive definite matrix `C` with `1 ≤ C` has `C⁻¹ ≤ 1`
(a spectral statement), and conjugating back by `S` gives `B⁻¹ ≤ A⁻¹`.
-/

@[expose] public section

open scoped MatrixOrder ComplexOrder

namespace Matrix

variable {𝕜 n : Type*} [RCLike 𝕜] [Fintype n] [DecidableEq n]

/-- A positive definite matrix `C` with `1 ≤ C` satisfies `C⁻¹ ≤ 1`. -/
lemma PosDef.inv_le_one_of_one_le {C : Matrix n n 𝕜} (hC : C.PosDef) (h : 1 ≤ C) : C⁻¹ ≤ 1 := by
  have hC' : IsSelfAdjoint C := hC.isHermitian.isSelfAdjoint
  have h_spec : ∀ x ∈ spectrum ℝ C, 1 ≤ x := (CFC.one_le_iff (R := ℝ) C hC').mp h
  have h_cont : ContinuousOn (fun x : ℝ ↦ x⁻¹) (spectrum ℝ C) :=
    continuousOn_id.inv₀ fun x hx ↦ (zero_lt_one.trans_le (h_spec x hx)).ne'
  rw [nonsing_inv_eq_ringInverse, ← cfc_ringInverse_id (R := ℝ) C hC.isUnit hC',
    cfc_le_one_iff _ C h_cont hC']
  exact fun x hx ↦ inv_le_one_of_one_le₀ (h_spec x hx)

/-- Matrix inversion is antitone on positive definite matrices for the Loewner order. -/
lemma PosDef.inv_le_inv {A B : Matrix n n 𝕜} (hA : A.PosDef) (hAB : A ≤ B) : B⁻¹ ≤ A⁻¹ := by
  -- `S` is the positive square root of `A⁻¹`
  set S := CFC.sqrt A⁻¹ with hS_def
  have hS_nonneg : 0 ≤ S := CFC.sqrt_nonneg _
  have hSS : S * S = A⁻¹ := CFC.sqrt_mul_sqrt_self A⁻¹ hA.inv.posSemidef.nonneg
  have hA_det : IsUnit A.det := (isUnit_iff_isUnit_det A).mp hA.isUnit
  have hS_det : IsUnit S.det := by
    have : IsUnit (S * S).det := by rw [hSS]; exact (isUnit_iff_isUnit_det _).mp hA.inv.isUnit
    rw [det_mul] at this
    exact (IsUnit.mul_iff.mp this).1
  have hS_unit : IsUnit S := (isUnit_iff_isUnit_det S).mpr hS_det
  -- `S * A * S = 1`
  have hSAS : S * A * S = 1 := by
    have hA' : A = S⁻¹ * S⁻¹ := by rw [← mul_inv_rev, hSS, nonsing_inv_nonsing_inv A hA_det]
    calc S * A * S = (S * S⁻¹) * (S⁻¹ * S) := by rw [hA']; simp only [mul_assoc]
      _ = 1 := by rw [mul_nonsing_inv S hS_det, nonsing_inv_mul S hS_det, one_mul]
  -- `C := S * B * S` satisfies `1 ≤ C` and is positive definite
  have hC : 1 ≤ S * B * S := hSAS ▸ conjugate_le_conjugate_of_nonneg hAB hS_nonneg
  have hC_pd : (S * B * S).PosDef := by
    have hB : B.PosDef := by
      simpa using hA.add_posSemidef (le_iff.mp hAB)
    exact (nonneg_iff_posSemidef.mp (zero_le_one.trans hC)).posDef_iff_isUnit.mpr
      ((hS_unit.mul hB.isUnit).mul hS_unit)
  -- conjugate `C⁻¹ ≤ 1` back by `S`
  have hSCS : S * (S * B * S)⁻¹ * S = B⁻¹ := by
    rw [mul_inv_rev, mul_inv_rev]
    calc S * (S⁻¹ * (B⁻¹ * S⁻¹)) * S = (S * S⁻¹) * B⁻¹ * (S⁻¹ * S) := by simp only [mul_assoc]
      _ = B⁻¹ := by rw [mul_nonsing_inv S hS_det, nonsing_inv_mul S hS_det, one_mul, mul_one]
  calc B⁻¹ = S * (S * B * S)⁻¹ * S := hSCS.symm
    _ ≤ S * 1 * S := conjugate_le_conjugate_of_nonneg (hC_pd.inv_le_one_of_one_le hC) hS_nonneg
    _ = A⁻¹ := by rw [mul_one, hSS]

end Matrix
