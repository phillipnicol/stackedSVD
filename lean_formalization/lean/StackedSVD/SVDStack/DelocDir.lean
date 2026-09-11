/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.SVDStack.Defs
import StackedSVD.RMT
import StackedSVD.Prob.GaussianMatrix

/-!
# The `delocDir` machinery of `lem:delocalization`

Moved out of `SVDStack/Gram.lean` (F28, 2026-09-08): the sections `DelocDir`,
`PiMarginal` and `TableFamily`, consumer-free (`perpOf` through
`measure_deloc_le_of_pi`). No proof changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

section DelocDir

variable {D : ℕ}

/-- Component of `x` orthogonal to the unit vector `v`. -/
noncomputable def perpOf (v x : EuclideanSpace ℝ (Fin D)) : EuclideanSpace ℝ (Fin D) :=
  x - ⟪v, x⟫_ℝ • v

/-- `perpOf v` is homogeneous of degree one. -/
theorem perpOf_smul (v : EuclideanSpace ℝ (Fin D)) (a : ℝ)
    (x : EuclideanSpace ℝ (Fin D)) : perpOf v (a • x) = a • perpOf v x := by
  simp only [perpOf, real_inner_smul_right, smul_sub, smul_smul]

/-- `perpOf v` is self-adjoint: `⟪perpOf v x, y⟫ = ⟪x, perpOf v y⟫`. -/
theorem inner_perpOf_comm (v x y : EuclideanSpace ℝ (Fin D)) :
    ⟪perpOf v x, y⟫_ℝ = ⟪x, perpOf v y⟫_ℝ := by
  simp only [perpOf, inner_sub_left, inner_sub_right, real_inner_smul_left,
    real_inner_smul_right]
  rw [real_inner_comm x v]
  ring

/-- `perpOf v x` is orthogonal to the unit vector `v`. -/
theorem inner_v_perpOf {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1)
    (x : EuclideanSpace ℝ (Fin D)) : ⟪v, perpOf v x⟫_ℝ = 0 := by
  simp only [perpOf, inner_sub_right, real_inner_smul_right, real_inner_self_eq_norm_sq, hv]
  ring

/-- `perpOf v` does not increase the norm. -/
theorem norm_perpOf_le {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1)
    (x : EuclideanSpace ℝ (Fin D)) : ‖perpOf v x‖ ≤ ‖x‖ := by
  have hx : x = perpOf v x + ⟪v, x⟫_ℝ • v := by
    simp only [perpOf]; abel
  have horth : ⟪perpOf v x, ⟪v, x⟫_ℝ • v⟫_ℝ = 0 := by
    rw [real_inner_smul_right, real_inner_comm v (perpOf v x), inner_v_perpOf hv, mul_zero]
  have hpy : ‖x‖ ^ 2 = ‖perpOf v x‖ ^ 2 + ‖⟪v, x⟫_ℝ • v‖ ^ 2 := by
    conv_lhs => rw [hx]
    rw [norm_add_sq_real, horth]
    ring
  nlinarith [norm_nonneg (perpOf v x), norm_nonneg x, sq_nonneg (‖⟪v, x⟫_ℝ • v‖)]

/-- `perpOf v` is continuous. -/
theorem continuous_perpOf (v : EuclideanSpace ℝ (Fin D)) : Continuous (perpOf v) := by
  unfold perpOf
  fun_prop

/-- Polarization: the coordinates of the top projector are differences of overlaps. -/
theorem topProj_apply_eq {r : ℕ} (X : Matrix (Fin r) (Fin D) ℝ)
    (z : EuclideanSpace ℝ (Fin D)) (k : Fin D) :
    (topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) z) k
      = (overlap X (z + EuclideanSpace.single k 1)
          - overlap X (z - EuclideanSpace.single k 1)) / 4 := by
  set K := topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) with hKdef
  set u : EuclideanSpace ℝ (Fin D) := EuclideanSpace.single k 1 with hu
  have hP : ∀ w : EuclideanSpace ℝ (Fin D),
      topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w = K.starProjection w := fun _ => rfl
  have hov : ∀ w : EuclideanSpace ℝ (Fin D), overlap X w = ‖K.starProjection w‖ ^ 2 :=
    fun _ => rfl
  have hinner : ⟪K.starProjection z, u⟫_ℝ = (K.starProjection z) k := by
    rw [hu, EuclideanSpace.inner_single_right]
    simp
  have hidem : K.starProjection (K.starProjection u) = K.starProjection u :=
    Submodule.starProjection_eq_self_iff.mpr (K.starProjection_apply_mem u)
  have h3 : ⟪K.starProjection z, K.starProjection u⟫_ℝ = ⟪K.starProjection z, u⟫_ℝ := by
    rw [Submodule.inner_starProjection_left_eq_right, hidem,
      ← Submodule.inner_starProjection_left_eq_right]
  have h1 : overlap X (z + u)
      = ‖K.starProjection z‖ ^ 2 + 2 * ⟪K.starProjection z, K.starProjection u⟫_ℝ
        + ‖K.starProjection u‖ ^ 2 := by
    rw [hov, map_add, norm_add_sq_real]
  have h2 : overlap X (z - u)
      = ‖K.starProjection z‖ ^ 2 - 2 * ⟪K.starProjection z, K.starProjection u⟫_ℝ
        + ‖K.starProjection u‖ ^ 2 := by
    rw [hov, map_sub, norm_sub_sq_real]
  rw [hP, h1, h2, h3, hinner]
  ring

/-- `X ↦ P_top(Xᵀ X) z` is measurable, by the polarization identity `topProj_apply_eq`. -/
theorem measurable_topProj_apply {r : ℕ} (z : EuclideanSpace ℝ (Fin D)) :
    Measurable (fun X : Matrix (Fin r) (Fin D) ℝ =>
      topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) z) := by
  have hrw : (fun X : Matrix (Fin r) (Fin D) ℝ =>
      topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) z)
      = fun X : Matrix (Fin r) (Fin D) ℝ => (WithLp.toLp 2
        (fun k => (topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) z) k) :
          EuclideanSpace ℝ (Fin D)) := rfl
  rw [hrw]
  refine (WithLp.measurable_toLp 2 _).comp (measurable_pi_lambda _ fun k => ?_)
  simp only [topProj_apply_eq]
  exact ((measurable_overlap _).sub (measurable_overlap _)).div_const 4

/-- `Pperp (P_top(X) (Pperp e_k))`: parallel to `Pperp v̂` on the simple event. -/
noncomputable def perpTopVec {r : ℕ} (v : EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin r) (Fin D) ℝ) (k : Fin D) : EuclideanSpace ℝ (Fin D) :=
  perpOf v (topProj (Xᵀ * X) (isHermitian_transpose_mul_self X)
    (perpOf v (EuclideanSpace.single k 1)))

/-- `X ↦ perpTopVec v X k` is measurable. -/
theorem measurable_perpTopVec {r : ℕ} (v : EuclideanSpace ℝ (Fin D)) (k : Fin D) :
    Measurable (fun X : Matrix (Fin r) (Fin D) ℝ => perpTopVec v X k) :=
  (continuous_perpOf v).measurable.comp (measurable_topProj_apply _)

open scoped Classical in
/-- A measurable unit vector orthogonal to `v`, parallel to `Pperp v̂(X)` on the simple event. -/
noncomputable def delocDir {r : ℕ} (v : EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin r) (Fin D) ℝ) : EuclideanSpace ℝ (Fin D) :=
  ∑ k : Fin D, if perpTopVec v X k ≠ 0 ∧ ∀ l, l < k → perpTopVec v X l = 0
    then ‖perpTopVec v X k‖⁻¹ • perpTopVec v X k else 0

/-- `delocDir v X` is either `0`, and then every `perpTopVec v X k` is `0`, or the normalized
`perpTopVec v X k` at the first index `k` where that vector is not `0`. -/
theorem delocDir_spec {r : ℕ} (v : EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin r) (Fin D) ℝ) :
    (delocDir v X = 0 ∧ ∀ k, perpTopVec v X k = 0) ∨
      ∃ k, perpTopVec v X k ≠ 0 ∧
        delocDir v X = ‖perpTopVec v X k‖⁻¹ • perpTopVec v X k := by
  classical
  by_cases h : ∃ k, perpTopVec v X k ≠ 0
  · right
    obtain ⟨k, hk⟩ := h
    obtain ⟨k₀, hk0mem, hlt⟩ : ∃ k₀ : Fin D, perpTopVec v X k₀ ≠ 0 ∧
        ∀ l, l < k₀ → perpTopVec v X l = 0 := by
      have hSne : (Finset.univ.filter (fun k => perpTopVec v X k ≠ 0)).Nonempty :=
        ⟨k, Finset.mem_filter.mpr ⟨Finset.mem_univ k, hk⟩⟩
      refine ⟨(Finset.univ.filter (fun k => perpTopVec v X k ≠ 0)).min' hSne, ?_, ?_⟩
      · exact (Finset.mem_filter.mp (Finset.min'_mem _ hSne)).2
      · intro l hl
        by_contra hne
        exact absurd (Finset.min'_le _ l (Finset.mem_filter.mpr ⟨Finset.mem_univ l, hne⟩))
          (not_le.mpr hl)
    refine ⟨k₀, hk0mem, ?_⟩
    have key : delocDir v X
        = if perpTopVec v X k₀ ≠ 0 ∧ ∀ l, l < k₀ → perpTopVec v X l = 0
          then ‖perpTopVec v X k₀‖⁻¹ • perpTopVec v X k₀ else 0 := by
      simp only [delocDir]
      refine Finset.sum_eq_single k₀ (fun b _ hb => ?_)
        (fun hb => absurd (Finset.mem_univ k₀) hb)
      rcases lt_or_gt_of_ne hb with hblt | hbgt
      · rw [if_neg]
        rintro ⟨h1, -⟩
        exact h1 (hlt b hblt)
      · rw [if_neg]
        rintro ⟨-, h2⟩
        exact hk0mem (h2 k₀ hbgt)
    rw [key, if_pos ⟨hk0mem, hlt⟩]
  · left
    simp only [not_exists, not_not] at h
    refine ⟨?_, h⟩
    simp only [delocDir]
    refine Finset.sum_eq_zero fun k _ => ?_
    rw [if_neg]
    rintro ⟨h1, -⟩
    exact h1 (h k)

/-- The zero set of `X ↦ perpTopVec v X l` is measurable. -/
theorem measurableSet_perpTopVec_eq_zero {r : ℕ} (v : EuclideanSpace ℝ (Fin D))
    (l : Fin D) :
    MeasurableSet {X : Matrix (Fin r) (Fin D) ℝ | perpTopVec v X l = 0} :=
  (measurable_perpTopVec v l) (measurableSet_singleton 0)

/-- `X ↦ delocDir v X` is measurable. This is what makes the random test direction of
`lem:delocalization` a function of one table alone. -/
theorem measurable_delocDir {r : ℕ} (v : EuclideanSpace ℝ (Fin D)) :
    Measurable (fun X : Matrix (Fin r) (Fin D) ℝ => delocDir v X) := by
  classical
  have hcond : ∀ k : Fin D, MeasurableSet {X : Matrix (Fin r) (Fin D) ℝ |
      perpTopVec v X k ≠ 0 ∧ ∀ l, l < k → perpTopVec v X l = 0} := by
    intro k
    have h1 : MeasurableSet {X : Matrix (Fin r) (Fin D) ℝ | perpTopVec v X k ≠ 0} :=
      (measurableSet_perpTopVec_eq_zero v k).compl
    have h2 : MeasurableSet {X : Matrix (Fin r) (Fin D) ℝ |
        ∀ l, l < k → perpTopVec v X l = 0} := by
      have hset : {X : Matrix (Fin r) (Fin D) ℝ | ∀ l, l < k → perpTopVec v X l = 0}
          = ⋂ l : {l : Fin D // l < k}, {X : Matrix (Fin r) (Fin D) ℝ |
            perpTopVec v X l.1 = 0} := by
        ext X
        simp
      rw [hset]
      exact MeasurableSet.iInter fun l => measurableSet_perpTopVec_eq_zero v l.1
    exact h1.inter h2
  simp only [delocDir]
  refine Finset.measurable_sum _ fun k _ => ?_
  refine Measurable.ite (hcond k) ?_ measurable_const
  exact ((measurable_perpTopVec v k).norm.inv).smul (measurable_perpTopVec v k)

/-- On the event that the top eigenvalue of `Xᵀ X` is simple with unit top eigenvector `e`,
`perpTopVec v X k` is a multiple of `perpOf v e`. -/
theorem perpTopVec_eq {r : ℕ} {v : EuclideanSpace ℝ (Fin D)}
    (X : Matrix (Fin r) (Fin D) ℝ)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {e : EuclideanSpace ℝ (Fin D)}
    (he : e ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) (hen : ‖e‖ = 1)
    (k : Fin D) : perpTopVec v X k = (perpOf v e k) • perpOf v e := by
  rw [perpTopVec, topProj_eq_rankOne hsimple he hen, perpOf_smul]
  congr 1
  rw [← inner_perpOf_comm, EuclideanSpace.inner_single_right]
  simp

/-- On the simple event `delocDir v X` is a unit vector parallel to `perpOf v e`, so it
carries the whole perpendicular part of any inner product with `e`. -/
theorem abs_inner_perpOf_le {r : ℕ} {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1)
    (X : Matrix (Fin r) (Fin D) ℝ)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {e : EuclideanSpace ℝ (Fin D)}
    (he : e ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) (hen : ‖e‖ = 1)
    (y : EuclideanSpace ℝ (Fin D)) :
    |⟪y, perpOf v e⟫_ℝ| ≤ |⟪y, delocDir v X⟫_ℝ| := by
  by_cases hq0 : perpOf v e = 0
  · rw [hq0]
    simp
  · have hqk : ∃ k, perpOf v e k ≠ 0 := by
      by_contra hc
      simp only [not_exists, not_not] at hc
      exact hq0 (by apply PiLp.ext; intro k; simpa using hc k)
    obtain ⟨k1, hk1⟩ := hqk
    rcases delocDir_spec v X with ⟨-, hall⟩ | ⟨k, hkne, hdir⟩
    · exact absurd (hall k1)
        (by rw [perpTopVec_eq X hsimple he hen k1]; exact smul_ne_zero hk1 hq0)
    · rw [perpTopVec_eq X hsimple he hen k] at hkne hdir
      have hqk0 : perpOf v e k ≠ 0 := fun h => hkne (by rw [h, zero_smul])
      have hnq : 0 < ‖perpOf v e‖ := norm_pos_iff.mpr hq0
      have hle1 : ‖perpOf v e‖ ≤ 1 := hen ▸ norm_perpOf_le hv e
      rw [hdir, real_inner_smul_right, real_inner_smul_right, norm_smul, Real.norm_eq_abs]
      rw [abs_mul, abs_inv, abs_mul, abs_mul, abs_abs, abs_of_pos hnq]
      rw [inv_mul_eq_div, le_div_iff₀ (mul_pos (abs_pos.mpr hqk0) hnq)]
      nlinarith [mul_nonneg (mul_nonneg (abs_nonneg (⟪y, perpOf v e⟫_ℝ))
        (abs_nonneg ((perpOf v e) k))) (sub_nonneg.mpr hle1)]

/-- When `delocDir v X` is not `0` it is a unit vector orthogonal to `v`. -/
theorem norm_delocDir {r : ℕ} {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1)
    (X : Matrix (Fin r) (Fin D) ℝ) (h : delocDir v X ≠ 0) :
    ‖delocDir v X‖ = 1 ∧ ⟪delocDir v X, v⟫_ℝ = 0 := by
  rcases delocDir_spec v X with ⟨h0, -⟩ | ⟨k, hkne, hdir⟩
  · exact absurd h0 h
  · constructor
    · rw [hdir, norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
      exact inv_mul_cancel₀ (norm_ne_zero_iff.mpr hkne)
    · rw [hdir, real_inner_smul_left, perpTopVec, real_inner_comm v, inner_v_perpOf hv,
        mul_zero]

/-- Fubini over a product measure: a set that reads only two independent coordinates. -/
theorem pi_pair_le {ι : Type*} [Fintype ι] {α : ι → Type*}
    [∀ k, MeasurableSpace (α k)] (P : ∀ k, Measure (α k)) (hP : ∀ k, IsProbabilityMeasure (P k))
    {i j : ι} (hij : i ≠ j) {S : Set (α i × α j)} (hS : MeasurableSet S) (K : ℝ≥0∞)
    (hK : ∀ zj : α j, P i {zi | (zi, zj) ∈ S} ≤ K) :
    Measure.pi P {Zs | (Zs i, Zs j) ∈ S} ≤ K := by
  classical
  have := hP
  have hmeas : Measurable fun Zs : (∀ k, α k) => (Zs i, Zs j) :=
    (measurable_pi_apply i).prodMk (measurable_pi_apply j)
  have hpush : (Measure.pi P).map (fun Zs => (Zs i, Zs j)) = (P i).prod (P j) := by
    refine (Measure.prod_eq fun s t hs ht => ?_).symm
    rw [Measure.map_apply hmeas (hs.prod ht)]
    have hpre : (fun Zs : (∀ k, α k) => (Zs i, Zs j)) ⁻¹' (s ×ˢ t)
        = Set.univ.pi (Function.update (Function.update
            (fun k => (Set.univ : Set (α k))) i s) j t) := by
      ext Zs
      simp only [Set.mem_preimage, Set.mem_prod, Set.mem_univ_pi]
      constructor
      · rintro ⟨h1, h2⟩ k
        rcases eq_or_ne k j with rfl | hkj
        · rw [Function.update_self]; exact h2
        · rw [Function.update_of_ne hkj]
          rcases eq_or_ne k i with rfl | hki
          · rw [Function.update_self]; exact h1
          · rw [Function.update_of_ne hki]; trivial
      · intro h
        refine ⟨?_, ?_⟩
        · have hi := h i
          rwa [Function.update_of_ne hij, Function.update_self] at hi
        · have hj := h j
          rwa [Function.update_self] at hj
    rw [hpre, Measure.pi_pi]
    rw [Finset.prod_eq_mul_of_mem i j (Finset.mem_univ i) (Finset.mem_univ j) hij ?_]
    · rw [Function.update_of_ne hij, Function.update_self, Function.update_self]
    · rintro k - ⟨hki, hkj⟩
      rw [Function.update_of_ne hkj, Function.update_of_ne hki]
      simp
  calc Measure.pi P {Zs | (Zs i, Zs j) ∈ S}
      = ((Measure.pi P).map (fun Zs => (Zs i, Zs j))) S := (Measure.map_apply hmeas hS).symm
    _ = ((P i).prod (P j)) S := by rw [hpush]
    _ = ∫⁻ zj, P i ((fun zi => (zi, zj)) ⁻¹' S) ∂(P j) := Measure.prod_apply_symm hS
    _ ≤ ∫⁻ _, K ∂(P j) := lintegral_mono hK
    _ = K := by simp

/-- `overlap X 0 = 0`. -/
theorem overlap_zero {r : ℕ} (X : Matrix (Fin r) (Fin D) ℝ) :
    overlap X (0 : EuclideanSpace ℝ (Fin D)) = 0 := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) 0‖ ^ 2 = 0
  rw [map_zero, norm_zero]
  ring

end DelocDir

section PiMarginal

/-! ### The marginal of a product law

`JointGaussianNoise` of `MultiTableModel` and of `UnalignedModel` both say that the family of
noise matrices has the product law `Measure.pi P`. The law of one coordinate is then `P i`.
Stated once here, for `Measure.pi`, so that neither model needs its own copy. -/

/-- One coordinate of a random point with a product law has the corresponding factor as its
law. `MeasureTheory.measurePreserving_eval` gives the push-forward of `Measure.pi P` under the
evaluation at `i`, and `HasLaw` composes with a measure-preserving map. -/
theorem hasLaw_eval_of_hasLaw_pi {Ω : Type*} [MeasurableSpace Ω] {ν : Measure Ω} {ι : Type*}
    [Fintype ι] {α : ι → Type*} [∀ k, MeasurableSpace (α k)] (P : ∀ k, Measure (α k))
    [∀ k, IsProbabilityMeasure (P k)] {f : Ω → ∀ k, α k}
    (hf : HasLaw f (Measure.pi P) ν) (i : ι) :
    HasLaw (fun ω => f ω i) (P i) ν :=
  (measurePreserving_eval P i).fun_comp_hasLaw hf

end PiMarginal

/-! ### 3b. The noise-to-data map and the Fubini step, model-free

`SpikedModel.dataOf` writes one table as a function of its noise matrix. The Fubini bound
`measure_deloc_le_of_pi` and the marginal `gaussianNoise_of_joint_pi` read only the family of
tables and the product law, never the rest of a model. Both `MultiTableModel` (below) and
`UnalignedModel` (`RankR/Defs.lean`) are one line each on top of them; before cleanup wave 2
`RankR/Defs.lean` repeated the whole argument, 108 lines. -/

section TableFamily

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

namespace SpikedModel

variable {n₁ : ℕ → ℕ}

/-- One table as a function of its noise matrix: `X = θ u vᵀ + d^{-1/2} Z`. -/
noncomputable def dataOf (t : SpikedModel μ n₁ d) (N : ℕ)
    (Z : Matrix (Fin (n₁ N)) (Fin (d N)) ℝ) : Matrix (Fin (n₁ N)) (Fin (d N)) ℝ :=
  t.θ • Matrix.vecMulVec (WithLp.ofLp (t.u N)) (WithLp.ofLp (t.v N))
    + (Real.sqrt (d N))⁻¹ • Z

/-- `dataOf` at the model noise is the data matrix of that table. -/
theorem dataOf_eq (t : SpikedModel μ n₁ d) (N : ℕ) (ω : Ω N) :
    t.dataOf N (t.Z N ω) = t.X N ω := rfl

/-- The noise-to-data map of one table is measurable. -/
theorem measurable_dataOf (t : SpikedModel μ n₁ d) (N : ℕ) : Measurable (t.dataOf N) := by
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun l => ?_
  simp only [SpikedModel.dataOf, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  have hz : Measurable fun Z : Matrix (Fin (n₁ N)) (Fin (d N)) ℝ => Z r l :=
    (measurable_pi_apply l).comp (measurable_pi_apply r)
  exact measurable_const.add (hz.const_mul _)

end SpikedModel

/-- Marginal of the product law: every table gets its own `SpikedModel.GaussianNoise`. This is
what lets a Layer 2 corollary discharge `SingleTableLaw` table by table with
`SpikedModel.singleTableLaw_of_gaussian`. -/
theorem gaussianNoise_of_joint_pi (tbl : (i : Fin M) → SpikedModel μ (n i) d)
    (hG : ∀ N, HasLaw (fun ω i => (tbl i).Z N ω)
      (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N)) (i : Fin M) :
    (tbl i).GaussianNoise := by
  intro N
  have : ∀ k : Fin M, IsProbabilityMeasure (gaussianMatrix (n k N) (d N)) := fun _ =>
    inferInstance
  exact hasLaw_eval_of_hasLaw_pi (fun k : Fin M => gaussianMatrix (n k N) (d N)) (hG N) i

/-- The Fubini step of `lem:delocalization`, model-free and for an arbitrary product law.
The probability that table `i` has a large overlap with the random direction
`delocDir v_i X_j` of table `j` is at most the supremum of the overlap over deterministic unit
directions orthogonal to `v_i`, which `SingleTableLaw.delocUniform` sends to `0`. Independence
across tables enters through `hind`; the laws `ν N i` of the single tables are arbitrary
probability measures, so Gaussian noise is one case of it (L4, 2026-09-02; mirror of
`measure_deloc_le_of_piR`, `RankR/GramR.lean`). -/
theorem measure_deloc_le_of_pi_indep (tbl : (i : Fin M) → SpikedModel μ (n i) d)
    {ν : (N : ℕ) → (i : Fin M) → Measure (Matrix (Fin (n i N)) (Fin (d N)) ℝ)}
    [hν : ∀ N i, IsProbabilityMeasure (ν N i)]
    (hind : ∀ N, HasLaw (fun ω i => (tbl i).Z N ω) (Measure.pi (ν N)) (μ N))
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) {η : ℝ} (hη : 0 < η) :
    μ N {ω | η ≤ overlap ((tbl i).X N ω)
        (delocDir ((tbl i).v N) ((tbl j).X N ω))}
      ≤ ⨆ w ∈ (tbl i).orthUnit N, μ N {ω | η ≤ overlap ((tbl i).X N ω) w} := by
  classical
  have hprob : ∀ k : Fin M, IsProbabilityMeasure (ν N k) := fun k => hν N k
  have hinst := hprob
  -- the law of one table's noise
  have hmarg : ∀ (k : Fin M) (p : Matrix (Fin (n k N)) (Fin (d N)) ℝ → Prop),
      MeasurableSet {z | p z} →
      μ N {ω | p ((tbl k).Z N ω)} = ν N k {z | p z} := by
    intro k p hp
    have h1 : μ N {ω | p ((tbl k).Z N ω)}
        = (Measure.pi fun l : Fin M => ν N l) {Zs | p (Zs k)} :=
      (hind N).measure_eq (p := fun Zs => p (Zs k)) ((measurable_pi_apply k) hp)
    rw [h1]
    exact (measurePreserving_eval (fun l : Fin M => ν N l)
      k).measure_preimage hp.nullMeasurableSet
  -- the bad set, as a set of pairs of noise matrices
  have hSmeas : MeasurableSet {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
      Matrix (Fin (n j N)) (Fin (d N)) ℝ |
      η ≤ overlap ((tbl i).dataOf N q.1) (delocDir ((tbl i).v N) ((tbl j).dataOf N q.2))} := by
    have h1 : Measurable fun q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ =>
        ((tbl i).dataOf N q.1, delocDir ((tbl i).v N) ((tbl j).dataOf N q.2)) :=
      (((tbl i).measurable_dataOf N).comp measurable_fst).prodMk
        ((measurable_delocDir ((tbl i).v N)).comp (((tbl j).measurable_dataOf N).comp
          measurable_snd))
    have hf := measurable_overlap₂.comp h1
    exact hf measurableSet_Ici
  have hpairmeas : MeasurableSet {Zs : (k : Fin M) → Matrix (Fin (n k N)) (Fin (d N)) ℝ |
      (Zs i, Zs j) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ |
        η ≤ overlap ((tbl i).dataOf N q.1) (delocDir ((tbl i).v N) ((tbl j).dataOf N q.2))}} :=
    ((measurable_pi_apply i).prodMk (measurable_pi_apply j)) hSmeas
  have key : μ N {ω | η ≤ overlap ((tbl i).X N ω)
        (delocDir ((tbl i).v N) ((tbl j).X N ω))}
      = (Measure.pi fun k : Fin M => ν N k)
        {Zs | (Zs i, Zs j) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n j N)) (Fin (d N)) ℝ |
          η ≤ overlap ((tbl i).dataOf N q.1) (delocDir ((tbl i).v N) ((tbl j).dataOf N q.2))}} :=
    (hind N).measure_eq hpairmeas
  rw [key]
  refine pi_pair_le (fun k : Fin M => ν N k) hprob hij hSmeas _ ?_
  intro zj
  by_cases hw : delocDir ((tbl i).v N) ((tbl j).dataOf N zj) = 0
  · have hempty : {zi : Matrix (Fin (n i N)) (Fin (d N)) ℝ |
        (zi, zj) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n j N)) (Fin (d N)) ℝ |
          η ≤ overlap ((tbl i).dataOf N q.1) (delocDir ((tbl i).v N) ((tbl j).dataOf N q.2))}}
        = (∅ : Set (Matrix (Fin (n i N)) (Fin (d N)) ℝ)) := by
      ext zi
      simp only [Set.mem_ofPred_eq, hw, overlap_zero, Set.mem_empty_iff_false, iff_false, not_le]
      exact hη
    rw [hempty]
    simp
  · have hmem : delocDir ((tbl i).v N) ((tbl j).dataOf N zj) ∈ (tbl i).orthUnit N :=
      norm_delocDir ((tbl i).hv N) _ hw
    have hsl : (ν N i) {zi : Matrix (Fin (n i N)) (Fin (d N)) ℝ |
        (zi, zj) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n j N)) (Fin (d N)) ℝ |
          η ≤ overlap ((tbl i).dataOf N q.1) (delocDir ((tbl i).v N) ((tbl j).dataOf N q.2))}}
        = μ N {ω | η ≤ overlap ((tbl i).X N ω)
          (delocDir ((tbl i).v N) ((tbl j).dataOf N zj))} := by
      refine (hmarg i (fun z => η ≤ overlap ((tbl i).dataOf N z)
        (delocDir ((tbl i).v N) ((tbl j).dataOf N zj))) ?_).symm
      have hf := (measurable_overlap
        (delocDir ((tbl i).v N) ((tbl j).dataOf N zj))).comp ((tbl i).measurable_dataOf N)
      exact hf measurableSet_Ici
    rw [hsl]
    exact le_iSup₂ (f := fun w (_ : w ∈ (tbl i).orthUnit N) =>
      μ N {ω | η ≤ overlap ((tbl i).X N ω) w}) _ hmem

/-- The Fubini step of `lem:delocalization` for Gaussian tables: `measure_deloc_le_of_pi_indep`
at `ν N k = gaussianMatrix (n k N) (d N)`. Kept in this form because `RankR/Defs.lean`
(`UnalignedModel.measure_deloc_le`) calls it. -/
theorem measure_deloc_le_of_pi (tbl : (i : Fin M) → SpikedModel μ (n i) d)
    (hG : ∀ N, HasLaw (fun ω i => (tbl i).Z N ω)
      (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N))
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) {η : ℝ} (hη : 0 < η) :
    μ N {ω | η ≤ overlap ((tbl i).X N ω)
        (delocDir ((tbl i).v N) ((tbl j).X N ω))}
      ≤ ⨆ w ∈ (tbl i).orthUnit N, μ N {ω | η ≤ overlap ((tbl i).X N ω) w} :=
  measure_deloc_le_of_pi_indep tbl (ν := fun N k => gaussianMatrix (n k N) (d N)) hG hij N hη

end TableFamily

end StackedSVD
