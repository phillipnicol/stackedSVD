/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Algebra.Order.BigOperators.Ring.Finset
public import Mathlib.Algebra.Order.Ring.Abs

set_option autoImplicit false

/-!
# A bound on weighted covariances

For a probability vector `p` on a finset `s` (`0 ≤ p i` and `∑ i ∈ s, p i = 1`) and two families
`a`, `b` bounded by `A` and `B` in absolute value on `s`, the weighted covariance
`∑ i ∈ s, p i * (a i * b i) - (∑ i ∈ s, p i * a i) * ∑ i ∈ s, p i * b i` is at most `A * B` in
absolute value (`Finset.abs_sum_mul_sub_mul_le`): it is the sum of the centered products
(`Finset.sum_mul_sub_mul_eq`), the weighted centered second moments are at most `A ^ 2` and
`B ^ 2` (`Finset.sum_mul_sq_sub_le`), and the Cauchy–Schwarz inequality concludes. Everything
holds in a linearly ordered commutative ring.
-/

@[expose] public section

namespace Finset

variable {ι R : Type*} [CommRing R] [LinearOrder R] [IsStrictOrderedRing R] {s : Finset ι}
  {p a b : ι → R} {A B : R}

/-- The weighted centered second moment of a family bounded by `C` is at most `C ^ 2`. -/
lemma sum_mul_sq_sub_le (hp : ∀ i ∈ s, 0 ≤ p i) (hp1 : ∑ i ∈ s, p i = 1) {c : ι → R} {C : R}
    (hc : ∀ i ∈ s, |c i| ≤ C) :
    ∑ i ∈ s, p i * (c i - ∑ j ∈ s, p j * c j) ^ 2 ≤ C ^ 2 := by
  have h1 : ∑ i ∈ s, p i * (c i - ∑ j ∈ s, p j * c j) ^ 2
      = ∑ i ∈ s, p i * c i ^ 2 - (∑ j ∈ s, p j * c j) ^ 2 := by
    have : ∀ i, p i * (c i - ∑ j ∈ s, p j * c j) ^ 2
        = p i * c i ^ 2 - (2 * ∑ j ∈ s, p j * c j) * (p i * c i) + (∑ j ∈ s, p j * c j) ^ 2 * p i :=
      fun i ↦ by ring
    simp_rw [this]
    rw [sum_add_distrib, sum_sub_distrib, ← mul_sum, ← mul_sum, hp1]
    ring
  calc ∑ i ∈ s, p i * (c i - ∑ j ∈ s, p j * c j) ^ 2
      ≤ ∑ i ∈ s, p i * c i ^ 2 := by rw [h1]; exact sub_le_self _ (sq_nonneg _)
    _ ≤ ∑ i ∈ s, p i * C ^ 2 := sum_le_sum fun i hi ↦ mul_le_mul_of_nonneg_left
        (sq_le_sq' (abs_le.1 (hc i hi)).1 (abs_le.1 (hc i hi)).2) (hp i hi)
    _ = C ^ 2 := by rw [← sum_mul, hp1, one_mul]

omit [LinearOrder R] [IsStrictOrderedRing R] in
/-- The weighted covariance is the weighted sum of the centered products. -/
lemma sum_mul_sub_mul_eq (hp1 : ∑ i ∈ s, p i = 1) :
    ∑ i ∈ s, p i * (a i * b i) - (∑ i ∈ s, p i * a i) * (∑ i ∈ s, p i * b i)
      = ∑ i ∈ s, p i * ((a i - ∑ j ∈ s, p j * a j) * (b i - ∑ j ∈ s, p j * b j)) := by
  have : ∀ i, p i * ((a i - ∑ j ∈ s, p j * a j) * (b i - ∑ j ∈ s, p j * b j))
      = p i * (a i * b i) - (∑ j ∈ s, p j * b j) * (p i * a i) - (∑ j ∈ s, p j * a j) * (p i * b i)
        + (∑ j ∈ s, p j * a j) * (∑ j ∈ s, p j * b j) * p i :=
    fun i ↦ by ring
  simp_rw [this]
  rw [sum_add_distrib, sum_sub_distrib, sum_sub_distrib, ← mul_sum, ← mul_sum, ← mul_sum, hp1]
  ring

/-- **Bound on a weighted covariance**: if `p` is a probability vector on `s` and `|a i| ≤ A`,
`|b i| ≤ B` on `s`, then
`|∑ i ∈ s, p i * (a i * b i) - (∑ i ∈ s, p i * a i) * ∑ i ∈ s, p i * b i| ≤ A * B`. -/
lemma abs_sum_mul_sub_mul_le (hp : ∀ i ∈ s, 0 ≤ p i) (hp1 : ∑ i ∈ s, p i = 1)
    (ha : ∀ i ∈ s, |a i| ≤ A) (hb : ∀ i ∈ s, |b i| ≤ B) :
    |∑ i ∈ s, p i * (a i * b i) - (∑ i ∈ s, p i * a i) * (∑ i ∈ s, p i * b i)| ≤ A * B := by
  obtain ⟨i₀, hi₀⟩ := nonempty_of_sum_ne_zero (hp1.trans_ne one_ne_zero)
  have hA : 0 ≤ A := (abs_nonneg _).trans (ha i₀ hi₀)
  have hB : 0 ≤ B := (abs_nonneg _).trans (hb i₀ hi₀)
  rw [sum_mul_sub_mul_eq hp1]
  refine abs_le_of_sq_le_sq ?_ (mul_nonneg hA hB)
  rw [mul_pow]
  refine (sum_sq_le_sum_mul_sum_of_sq_le_mul s
    (f := fun i ↦ p i * (a i - ∑ j ∈ s, p j * a j) ^ 2)
    (g := fun i ↦ p i * (b i - ∑ j ∈ s, p j * b j) ^ 2)
    (fun i hi ↦ mul_nonneg (hp i hi) (sq_nonneg _))
    (fun i hi ↦ mul_nonneg (hp i hi) (sq_nonneg _)) fun i _ ↦ le_of_eq (by ring)).trans ?_
  exact mul_le_mul (sum_mul_sq_sub_le hp hp1 ha) (sum_mul_sq_sub_le hp hp1 hb)
    (sum_nonneg fun i hi ↦ mul_nonneg (hp i hi) (sq_nonneg _)) (sq_nonneg _)

end Finset
