/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import Mathlib.Analysis.InnerProductSpace.Basic
public import Mathlib.Analysis.InnerProductSpace.Continuous
public import Mathlib.Data.Fintype.Order
public import Mathlib.Analysis.Convex.Function
public import Mathlib.Analysis.Normed.Group.Pointwise
public import Mathlib.Topology.Order.Compact

set_option autoImplicit false

/-!
# The support function of a set in an inner product space

For a set `K` in a real inner product space `E`, the *support function* of `K` is
`supportFn K ξ = sup_{x ∈ K} ⟪x, ξ⟫`. For a nonempty set bounded by `R` it is `R`-Lipschitz,
convex and positively homogeneous, and it is attained on a compact set.

## Main results

* `inner_le_supportFn`, `supportFn_le`, `abs_supportFn_le`: the defining inequalities;
* `lipschitzWith_supportFn`, `convexOn_supportFn`, `supportFn_smul`, `supportFn_image`,
  `supportFn_add`, `supportFn_mono`: the basic properties of the support function;
* `exists_supportFn_eq_inner`: on a compact set, the support function is attained.

## TODO

The support function makes sense for a set of a normed space, as a function on the dual, and
for a set of the dual as a function on the space (`supportFn K L = sup_{x ∈ K} L x`); the
inner-product formulation is the one used by the Gaussian-width developments.
-/

@[expose] public section

open Set

open scoped RealInnerProductSpace Pointwise NNReal

section SupportFn

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The support function of a set `K`: `supportFn K ξ = sup_{x ∈ K} ⟪x, ξ⟫`. -/
noncomputable def supportFn (K : Set E) (ξ : E) : ℝ := ⨆ x : K, ⟪(x : E), ξ⟫

variable {K K' : Set E} {R : ℝ} {ξ : E}

lemma supportFn_eq_sSup (K : Set E) (ξ : E) : supportFn K ξ = sSup ((fun x ↦ ⟪x, ξ⟫) '' K) := by
  rw [supportFn, iSup, image_eq_range]

lemma supportFn_of_subsingleton [Subsingleton E] (K : Set E) (ξ : E) : supportFn K ξ = 0 := by
  simp [supportFn, Subsingleton.elim ξ 0, Real.iSup_const_zero]

lemma bddAbove_range_inner (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ : E) :
    BddAbove (range fun x : K ↦ ⟪(x : E), ξ⟫) := by
  refine ⟨R * ‖ξ‖, ?_⟩
  rintro _ ⟨x, rfl⟩
  exact (real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hK x x.2) (norm_nonneg _))

lemma inner_le_supportFn (hK : ∀ x ∈ K, ‖x‖ ≤ R) {x : E} (hx : x ∈ K) (ξ : E) :
    ⟪x, ξ⟫ ≤ supportFn K ξ :=
  le_ciSup (bddAbove_range_inner hK ξ) ⟨x, hx⟩

lemma supportFn_le (hne : K.Nonempty) {c : ℝ} (h : ∀ x ∈ K, ⟪x, ξ⟫ ≤ c) : supportFn K ξ ≤ c := by
  have : Nonempty K := hne.to_subtype
  exact ciSup_le fun x ↦ h x x.2

lemma supportFn_le_mul_norm (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ : E) :
    supportFn K ξ ≤ R * ‖ξ‖ :=
  supportFn_le hne fun x hx ↦
    (real_inner_le_norm x ξ).trans (mul_le_mul_of_nonneg_right (hK x hx) (norm_nonneg _))

lemma neg_mul_norm_le_supportFn (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ : E) :
    -(R * ‖ξ‖) ≤ supportFn K ξ := by
  obtain ⟨x, hx⟩ := hne
  refine le_trans ?_ (inner_le_supportFn hK hx ξ)
  have h1 := abs_real_inner_le_norm x ξ
  have h2 : ‖x‖ * ‖ξ‖ ≤ R * ‖ξ‖ := mul_le_mul_of_nonneg_right (hK x hx) (norm_nonneg _)
  linarith [(abs_le.1 h1).1]

lemma abs_supportFn_le (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ : E) :
    |supportFn K ξ| ≤ R * ‖ξ‖ :=
  abs_le.2 ⟨neg_mul_norm_le_supportFn hne hK ξ, supportFn_le_mul_norm hne hK ξ⟩

lemma supportFn_sub_le (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ ξ' : E) :
    supportFn K ξ - supportFn K ξ' ≤ R * ‖ξ - ξ'‖ := by
  rw [sub_le_iff_le_add']
  refine supportFn_le hne fun x hx ↦ ?_
  calc ⟪x, ξ⟫ = ⟪x, ξ'⟫ + ⟪x, ξ - ξ'⟫ := by rw [← inner_add_right, add_sub_cancel]
    _ ≤ supportFn K ξ' + R * ‖ξ - ξ'‖ :=
      add_le_add (inner_le_supportFn hK hx ξ') ((real_inner_le_norm _ _).trans
        (mul_le_mul_of_nonneg_right (hK x hx) (norm_nonneg _)))

omit [InnerProductSpace ℝ E] in
lemma nonneg_of_norm_le (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) : 0 ≤ R :=
  hne.elim fun x hx ↦ (norm_nonneg x).trans (hK x hx)

lemma lipschitzWith_supportFn (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) :
    LipschitzWith R.toNNReal (supportFn K) := by
  refine LipschitzWith.of_dist_le_mul fun ξ ξ' ↦ ?_
  rw [Real.dist_eq, Real.coe_toNNReal R (nonneg_of_norm_le hne hK), dist_eq_norm]
  refine abs_sub_le_iff.2 ⟨supportFn_sub_le hne hK ξ ξ', ?_⟩
  simpa [norm_sub_rev] using supportFn_sub_le hne hK ξ' ξ

lemma continuous_supportFn (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) :
    Continuous (supportFn K) :=
  (lipschitzWith_supportFn hne hK).continuous

lemma supportFn_add_le (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (ξ ξ' : E) :
    supportFn K (ξ + ξ') ≤ supportFn K ξ + supportFn K ξ' :=
  supportFn_le hne fun x hx ↦ by
    rw [inner_add_right]
    exact add_le_add (inner_le_supportFn hK hx ξ) (inner_le_supportFn hK hx ξ')

/-- The support function is positively homogeneous (this holds for every set `K`, by the
convention `sSup ∅ = 0` of `ℝ`). -/
lemma supportFn_smul (K : Set E) {c : ℝ} (hc : 0 ≤ c) (ξ : E) :
    supportFn K (c • ξ) = c * supportFn K ξ := by
  simp only [supportFn, inner_smul_right]
  rw [Real.mul_iSup_of_nonneg hc]

lemma convexOn_supportFn (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) :
    ConvexOn ℝ univ (supportFn K) := by
  refine ⟨convex_univ, fun ξ _ ξ' _ a b ha hb _ ↦ ?_⟩
  calc supportFn K (a • ξ + b • ξ') ≤ supportFn K (a • ξ) + supportFn K (b • ξ') :=
        supportFn_add_le hne hK _ _
    _ = a • supportFn K ξ + b • supportFn K ξ' := by
        rw [supportFn_smul K ha, supportFn_smul K hb, smul_eq_mul, smul_eq_mul]

/-- The support function of an image `L '' K` is the support function of `K` composed with the
adjoint `L'` of `L`. -/
lemma supportFn_image {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F]
    (K : Set E) {L : E → F} {L' : F → E} (h : ∀ x ξ, ⟪L x, ξ⟫ = ⟪x, L' ξ⟫) (ξ : F) :
    supportFn (L '' K) ξ = supportFn K (L' ξ) := by
  rw [supportFn_eq_sSup, supportFn_eq_sSup, image_image]
  simp_rw [h]

lemma supportFn_neg_set (K : Set E) (ξ : E) : supportFn (-K) ξ = supportFn K (-ξ) := by
  rw [← image_neg_eq_neg]
  exact supportFn_image K (fun _ _ ↦ by rw [inner_neg_left, inner_neg_right]) ξ

lemma supportFn_smul_set (K : Set E) {c : ℝ} (hc : 0 ≤ c) (ξ : E) :
    supportFn (c • K) ξ = c * supportFn K ξ := by
  rw [← image_smul, supportFn_image K (L' := fun ξ ↦ c • ξ) (fun _ _ ↦ by
    simp [inner_smul_left, inner_smul_right]), supportFn_smul K hc]

lemma supportFn_mono (hne : K.Nonempty) (hK' : ∀ x ∈ K', ‖x‖ ≤ R) (hsub : K ⊆ K') (ξ : E) :
    supportFn K ξ ≤ supportFn K' ξ :=
  supportFn_le hne fun _ hx ↦ inner_le_supportFn hK' (hsub hx) ξ

/-- The support function of a Minkowski sum is the sum of the support functions. -/
lemma supportFn_add {R' : ℝ} (hne : K.Nonempty) (hK : ∀ x ∈ K, ‖x‖ ≤ R) (hne' : K'.Nonempty)
    (hK' : ∀ x ∈ K', ‖x‖ ≤ R') (ξ : E) :
    supportFn (K + K') ξ = supportFn K ξ + supportFn K' ξ := by
  have hKK' : ∀ z ∈ K + K', ‖z‖ ≤ R + R' := by
    rintro _ ⟨x, hx, y, hy, rfl⟩
    exact (norm_add_le _ _).trans (add_le_add (hK x hx) (hK' y hy))
  refine le_antisymm (supportFn_le (hne.add hne') ?_) ?_
  · rintro _ ⟨x, hx, y, hy, rfl⟩
    rw [inner_add_left]
    exact add_le_add (inner_le_supportFn hK hx ξ) (inner_le_supportFn hK' hy ξ)
  · have h : supportFn K' ξ ≤ supportFn (K + K') ξ - supportFn K ξ := by
      refine supportFn_le hne' fun y hy ↦ ?_
      rw [le_sub_comm]
      refine supportFn_le hne fun x hx ↦ ?_
      rw [le_sub_iff_add_le, ← inner_add_left]
      exact inner_le_supportFn hKK' (add_mem_add hx hy) ξ
    linarith

/-- On a compact set, the support function is attained. -/
lemma exists_supportFn_eq_inner (hK : IsCompact K) (hne : K.Nonempty) (ξ : E) :
    ∃ x ∈ K, supportFn K ξ = ⟪x, ξ⟫ := by
  have hc : Continuous fun x : E ↦ ⟪x, ξ⟫ := continuous_id.inner continuous_const
  obtain ⟨x, hx, hmax⟩ := hK.exists_isMaxOn hne hc.continuousOn
  obtain ⟨R, hR⟩ := hK.isBounded.exists_norm_le
  exact ⟨x, hx, le_antisymm (supportFn_le hne fun y hy ↦ hmax hy) (inner_le_supportFn hR hx ξ)⟩

section finite

variable {ι : Type*} [Finite ι] [Nonempty ι] {y : ι → E} {σ : ℝ≥0}

lemma abs_ciSup_inner_le (hσ : ∀ i, ‖y i‖ ≤ σ) (v : E) :
    |⨆ i, ⟪y i, v⟫| ≤ σ * ‖v‖ := by
  have hb : ∀ i, |⟪y i, v⟫| ≤ σ * ‖v‖ := fun i ↦
    (abs_real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hσ i) (norm_nonneg _))
  refine abs_le.2 ⟨?_, ciSup_le fun i ↦ (le_abs_self _).trans (hb i)⟩
  obtain ⟨i⟩ := ‹Nonempty ι›
  exact (neg_le.2 ((neg_le_abs _).trans (hb i))).trans
    (le_ciSup (Finite.bddAbove_range fun i ↦ ⟪y i, v⟫) i)

end finite

end SupportFn
