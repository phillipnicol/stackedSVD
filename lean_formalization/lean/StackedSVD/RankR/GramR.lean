/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.General
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# The per-table delocalization machinery at general `r_i`

Task B2b of `notes/archive/rankr_plan_B.md` section 2. Every result here is the rank-`r_i` twin of a
result of `StackedSVD/SVDStack/Gram.lean`, with three replacements.

1. The single spike direction `v` becomes the whole signal frame `col 1, …, col rk` of the
   table, so `perpOf v` becomes `perpFrame v`, the component orthogonal to the span of the
   frame. The frame is orthonormal by `SpikedModelR.hV`.
2. `topProj` at the top eigenvalue becomes `specProjIdx` at the eigenvalue index `j`, and
   `overlap` becomes `overlapIdx`; the measurability chain runs on
   `measurable_overlapIdx₂ j` and `measurable_overlapIdx j w` (`LinAlg/SpecIdxMeas.lean`).
3. `TopSimple` becomes `SimpleSpec … rk'`, through `specProjIdx_eq_rankOne`
   (`LinAlg/SpecIdx.lean`).

## Content

1. `perpFrame` and its five algebra lemmas (mirror `Gram.lean:47-86`).
2. `specProjIdx_apply_eq`, `measurable_specProjIdx_apply` (mirror `:88-131`).
3. `perpTopVecIdx`, `measurable_perpTopVecIdx` (mirror `:133-143`).
4. `delocDirIdx`, `delocDirIdx_spec`, `measurableSet_perpTopVecIdx_eq_zero`,
   `measurable_delocDirIdx` (mirror `:145-226`).
5. `specProjIdx_eq_rankOne_of_mem`, `perpTopVecIdx_eq`, `abs_inner_perpFrame_le`,
   `norm_delocDirIdx` (mirror `:228-281`).
6. `SpikedModelR.perpSpan`, `inner_col`, `dataOf`, `dataOf_eq`, `measurable_dataOf`
   (mirror `:372-392`).
7. `gaussianNoise_of_joint_piR` and `measure_deloc_le_of_piR`, the Fubini step
   (mirror `:394-490`).
8. `UnalignedModelR.IndepNoise`, the paper's `assum:general_noise`, and the fact that
   `JointGaussianNoise` implies it. The Fubini step needs only independence across tables,
   so every Layer 1 theorem of Track B takes `IndepNoise`.

`pi_pair_le` (`Gram.lean:282`) and `hasLaw_eval_of_hasLaw_pi` (`:345`) are model-free and
public, so they are reused unchanged; no primed copy is needed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

section DelocDirIdx

variable {D rk : ℕ}

/-! ### 1. The component orthogonal to a frame -/

/-- Component of `x` orthogonal to the span of the frame `v`. At `rk = 1` it is `perpOf (v 0)`
of `SVDStack/Gram.lean`. -/
noncomputable def perpFrame (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (x : EuclideanSpace ℝ (Fin D)) : EuclideanSpace ℝ (Fin D) :=
  x - ∑ k, ⟪v k, x⟫_ℝ • v k

/-- `perpFrame v` is homogeneous of degree one. Mirror: `perpOf_smul`. -/
theorem perpFrame_smul (v : Fin rk → EuclideanSpace ℝ (Fin D)) (a : ℝ)
    (x : EuclideanSpace ℝ (Fin D)) : perpFrame v (a • x) = a • perpFrame v x := by
  simp only [perpFrame, real_inner_smul_right, smul_sub, Finset.smul_sum, smul_smul]

/-- `perpFrame v` is self-adjoint. Mirror: `inner_perpOf_comm`. No orthonormality is needed. -/
theorem inner_perpFrame_comm (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (x y : EuclideanSpace ℝ (Fin D)) :
    ⟪perpFrame v x, y⟫_ℝ = ⟪x, perpFrame v y⟫_ℝ := by
  simp only [perpFrame, inner_sub_left, inner_sub_right, sum_inner, inner_sum,
    real_inner_smul_left, real_inner_smul_right]
  congr 1
  exact Finset.sum_congr rfl fun k _ => by rw [real_inner_comm x (v k)]; ring

/-- `perpFrame v x` is orthogonal to every member of an orthonormal frame `v`.
Mirror: `inner_v_perpOf`. -/
theorem inner_v_perpFrame {v : Fin rk → EuclideanSpace ℝ (Fin D)}
    (hv : ∀ k l, ⟪v k, v l⟫_ℝ = if k = l then 1 else 0) (l : Fin rk)
    (x : EuclideanSpace ℝ (Fin D)) : ⟪v l, perpFrame v x⟫_ℝ = 0 := by
  simp only [perpFrame, inner_sub_right, inner_sum, real_inner_smul_right, hv]
  simp

/-- `perpFrame v` does not increase the norm. Mirror: `norm_perpOf_le`. -/
theorem norm_perpFrame_le {v : Fin rk → EuclideanSpace ℝ (Fin D)}
    (hv : ∀ k l, ⟪v k, v l⟫_ℝ = if k = l then 1 else 0)
    (x : EuclideanSpace ℝ (Fin D)) : ‖perpFrame v x‖ ≤ ‖x‖ := by
  have hx : x = perpFrame v x + ∑ k, ⟪v k, x⟫_ℝ • v k := by
    simp only [perpFrame]; abel
  have horth : ⟪perpFrame v x, ∑ k, ⟪v k, x⟫_ℝ • v k⟫_ℝ = 0 := by
    rw [inner_sum]
    refine Finset.sum_eq_zero fun k _ => ?_
    rw [real_inner_smul_right, real_inner_comm (v k) (perpFrame v x),
      inner_v_perpFrame hv k x, mul_zero]
  have hpy : ‖x‖ ^ 2 = ‖perpFrame v x‖ ^ 2 + ‖∑ k, ⟪v k, x⟫_ℝ • v k‖ ^ 2 := by
    conv_lhs => rw [hx]
    rw [norm_add_sq_real, horth]
    ring
  nlinarith [norm_nonneg (perpFrame v x), norm_nonneg x,
    sq_nonneg ‖∑ k, ⟪v k, x⟫_ℝ • v k‖]

/-- `perpFrame v` is continuous. Mirror: `continuous_perpOf`. -/
theorem continuous_perpFrame (v : Fin rk → EuclideanSpace ℝ (Fin D)) :
    Continuous (perpFrame v) := by
  unfold perpFrame
  fun_prop

/-! ### 2. The projector at an eigenvalue index, coordinate by coordinate -/

/-- Polarization: the coordinates of `specProjIdx` are differences of `overlapIdx` values.
Mirror: `topProj_apply_eq`. -/
theorem specProjIdx_apply_eq {q : ℕ} (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ)
    (z : EuclideanSpace ℝ (Fin D)) (k : Fin D) :
    (specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j z) k
      = (overlapIdx X j (z + EuclideanSpace.single k 1)
          - overlapIdx X j (z - EuclideanSpace.single k 1)) / 4 := by
  set K := specSpace (Xᵀ * X)
    (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j) with hKdef
  set u : EuclideanSpace ℝ (Fin D) := EuclideanSpace.single k 1 with hu
  have hP : ∀ w : EuclideanSpace ℝ (Fin D),
      specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j w = K.starProjection w :=
    fun _ => rfl
  have hov : ∀ w : EuclideanSpace ℝ (Fin D),
      overlapIdx X j w = ‖K.starProjection w‖ ^ 2 := fun _ => rfl
  have hinner : ⟪K.starProjection z, u⟫_ℝ = (K.starProjection z) k := by
    rw [hu, EuclideanSpace.inner_single_right]
    simp
  have hidem : K.starProjection (K.starProjection u) = K.starProjection u :=
    Submodule.starProjection_eq_self_iff.mpr (K.starProjection_apply_mem u)
  have h3 : ⟪K.starProjection z, K.starProjection u⟫_ℝ = ⟪K.starProjection z, u⟫_ℝ := by
    rw [Submodule.inner_starProjection_left_eq_right, hidem,
      ← Submodule.inner_starProjection_left_eq_right]
  have h1 : overlapIdx X j (z + u)
      = ‖K.starProjection z‖ ^ 2 + 2 * ⟪K.starProjection z, K.starProjection u⟫_ℝ
        + ‖K.starProjection u‖ ^ 2 := by
    rw [hov, map_add, norm_add_sq_real]
  have h2 : overlapIdx X j (z - u)
      = ‖K.starProjection z‖ ^ 2 - 2 * ⟪K.starProjection z, K.starProjection u⟫_ℝ
        + ‖K.starProjection u‖ ^ 2 := by
    rw [hov, map_sub, norm_sub_sq_real]
  rw [hP, h1, h2, h3, hinner]
  ring

/-- `X ↦ specProjIdx (Xᵀ X) j z` is measurable, by `specProjIdx_apply_eq` and
`measurable_overlapIdx`. Mirror: `measurable_topProj_apply`. -/
theorem measurable_specProjIdx_apply {q : ℕ} (j : ℕ) (z : EuclideanSpace ℝ (Fin D)) :
    Measurable (fun X : Matrix (Fin q) (Fin D) ℝ =>
      specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j z) := by
  have hrw : (fun X : Matrix (Fin q) (Fin D) ℝ =>
      specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j z)
      = fun X : Matrix (Fin q) (Fin D) ℝ => (WithLp.toLp 2
        (fun k => (specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j z) k) :
          EuclideanSpace ℝ (Fin D)) := rfl
  rw [hrw]
  refine (WithLp.measurable_toLp 2 _).comp (measurable_pi_lambda _ fun k => ?_)
  simp only [specProjIdx_apply_eq]
  exact ((measurable_overlapIdx j _).sub (measurable_overlapIdx j _)).div_const 4

/-! ### 3. The perpendicular part of the index-`j` projector -/

/-- `Pperp (specProjIdx (XᵀX) j (Pperp e_k))`: parallel to `Pperp v̂_j` on the simple event.
Mirror: `perpTopVec`. -/
noncomputable def perpTopVecIdx {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ) (k : Fin D) : EuclideanSpace ℝ (Fin D) :=
  perpFrame v (specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j
    (perpFrame v (EuclideanSpace.single k 1)))

/-- `X ↦ perpTopVecIdx v X j k` is measurable. Mirror: `measurable_perpTopVec`. -/
theorem measurable_perpTopVecIdx {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D)) (j : ℕ)
    (k : Fin D) :
    Measurable (fun X : Matrix (Fin q) (Fin D) ℝ => perpTopVecIdx v X j k) :=
  (continuous_perpFrame v).measurable.comp (measurable_specProjIdx_apply j _)

/-! ### 4. The measurable random test direction -/

open scoped Classical in
/-- A measurable unit vector orthogonal to the frame `v`, parallel to `Pperp v̂_j(X)` on the
simple event. Mirror: `delocDir`. -/
noncomputable def delocDirIdx {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ) : EuclideanSpace ℝ (Fin D) :=
  ∑ k : Fin D, if perpTopVecIdx v X j k ≠ 0 ∧ ∀ l, l < k → perpTopVecIdx v X j l = 0
    then ‖perpTopVecIdx v X j k‖⁻¹ • perpTopVecIdx v X j k else 0

/-- `delocDirIdx v X j` is either `0`, and then every `perpTopVecIdx v X j k` is `0`, or the
normalized `perpTopVecIdx v X j k` at the first index `k` where that vector is not `0`.
Mirror: `delocDir_spec`. -/
theorem delocDirIdx_spec {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ) :
    (delocDirIdx v X j = 0 ∧ ∀ k, perpTopVecIdx v X j k = 0) ∨
      ∃ k, perpTopVecIdx v X j k ≠ 0 ∧
        delocDirIdx v X j = ‖perpTopVecIdx v X j k‖⁻¹ • perpTopVecIdx v X j k := by
  classical
  by_cases h : ∃ k, perpTopVecIdx v X j k ≠ 0
  · right
    obtain ⟨k, hk⟩ := h
    obtain ⟨k₀, hk0mem, hlt⟩ : ∃ k₀ : Fin D, perpTopVecIdx v X j k₀ ≠ 0 ∧
        ∀ l, l < k₀ → perpTopVecIdx v X j l = 0 := by
      have hSne : (Finset.univ.filter (fun k => perpTopVecIdx v X j k ≠ 0)).Nonempty :=
        ⟨k, Finset.mem_filter.mpr ⟨Finset.mem_univ k, hk⟩⟩
      refine ⟨(Finset.univ.filter (fun k => perpTopVecIdx v X j k ≠ 0)).min' hSne, ?_, ?_⟩
      · exact (Finset.mem_filter.mp (Finset.min'_mem _ hSne)).2
      · intro l hl
        by_contra hne
        exact absurd (Finset.min'_le _ l (Finset.mem_filter.mpr ⟨Finset.mem_univ l, hne⟩))
          (not_le.mpr hl)
    refine ⟨k₀, hk0mem, ?_⟩
    have key : delocDirIdx v X j
        = if perpTopVecIdx v X j k₀ ≠ 0 ∧ ∀ l, l < k₀ → perpTopVecIdx v X j l = 0
          then ‖perpTopVecIdx v X j k₀‖⁻¹ • perpTopVecIdx v X j k₀ else 0 := by
      simp only [delocDirIdx]
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
    simp only [delocDirIdx]
    refine Finset.sum_eq_zero fun k _ => ?_
    rw [if_neg]
    rintro ⟨h1, -⟩
    exact h1 (h k)

/-- The zero set of `X ↦ perpTopVecIdx v X j l` is measurable.
Mirror: `measurableSet_perpTopVec_eq_zero`. -/
theorem measurableSet_perpTopVecIdx_eq_zero {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D))
    (j : ℕ) (l : Fin D) :
    MeasurableSet {X : Matrix (Fin q) (Fin D) ℝ | perpTopVecIdx v X j l = 0} :=
  (measurable_perpTopVecIdx v j l) (measurableSet_singleton 0)

/-- `X ↦ delocDirIdx v X j` is measurable. This is what makes the random test direction of
`lem:general_rank_delocalization` a function of one table alone.
Mirror: `measurable_delocDir`. -/
theorem measurable_delocDirIdx {q : ℕ} (v : Fin rk → EuclideanSpace ℝ (Fin D)) (j : ℕ) :
    Measurable (fun X : Matrix (Fin q) (Fin D) ℝ => delocDirIdx v X j) := by
  classical
  have hcond : ∀ k : Fin D, MeasurableSet {X : Matrix (Fin q) (Fin D) ℝ |
      perpTopVecIdx v X j k ≠ 0 ∧ ∀ l, l < k → perpTopVecIdx v X j l = 0} := by
    intro k
    have h1 : MeasurableSet {X : Matrix (Fin q) (Fin D) ℝ | perpTopVecIdx v X j k ≠ 0} :=
      (measurableSet_perpTopVecIdx_eq_zero v j k).compl
    have h2 : MeasurableSet {X : Matrix (Fin q) (Fin D) ℝ |
        ∀ l, l < k → perpTopVecIdx v X j l = 0} := by
      have hset : {X : Matrix (Fin q) (Fin D) ℝ | ∀ l, l < k → perpTopVecIdx v X j l = 0}
          = ⋂ l : {l : Fin D // l < k}, {X : Matrix (Fin q) (Fin D) ℝ |
            perpTopVecIdx v X j l.1 = 0} := by
        ext X
        simp
      rw [hset]
      exact MeasurableSet.iInter fun l => measurableSet_perpTopVecIdx_eq_zero v j l.1
    exact h1.inter h2
  simp only [delocDirIdx]
  refine Finset.measurable_sum _ fun k _ => ?_
  refine Measurable.ite (hcond k) ?_ measurable_const
  exact ((measurable_perpTopVecIdx v j k).norm.inv).smul (measurable_perpTopVecIdx v j k)

/-! ### 5. The rank-one collapse on the simple event -/

/-- Under `SimpleSpec`, the projector at an index `j < rk'` is the rank-one projector on **any**
unit vector of the index-`j` spectral subspace, not only on `vEig`. This is what lets a
consumer feed the sign-corrected `v̂_ij` in place of `vEig`.
Mirror: `topProj_eq_rankOne` (`SVDStack/Defs.lean:152`). -/
theorem specProjIdx_eq_rankOne_of_mem {A : Matrix (Fin D) (Fin D) ℝ} {hA : A.IsHermitian}
    {rk' j : ℕ} (hsimple : SimpleSpec A hA rk') (hj : j < rk')
    {e : EuclideanSpace ℝ (Fin D)} (he : e ∈ specSpace A (eigSetIdx A hA j)) (hen : ‖e‖ = 1)
    (w : EuclideanSpace ℝ (Fin D)) : specProjIdx A hA j w = ⟪e, w⟫_ℝ • e := by
  have h1 : specProjIdx A hA j e = ⟪vEig A hA j, e⟫_ℝ • vEig A hA j :=
    specProjIdx_eq_rankOne hsimple hj e
  have h2 : specProjIdx A hA j e = e := by
    change (specSpace A (eigSetIdx A hA j)).starProjection e = e
    exact Submodule.starProjection_eq_self_iff.mpr he
  have he' : e = ⟪vEig A hA j, e⟫_ℝ • vEig A hA j := by rw [← h1, h2]
  have hqne : vEig A hA j ≠ 0 := by
    intro h0
    rw [h0, smul_zero] at he'
    rw [he', norm_zero] at hen
    norm_num at hen
  have hjc : j < Fintype.card (Fin D) := by
    by_contra hc
    exact hqne (vEig_of_not_lt A hA hc)
  have hqn : ‖vEig A hA j‖ = 1 := norm_vEig A hA hjc
  obtain ⟨a, hea⟩ : ∃ a : ℝ, e = a • vEig A hA j := ⟨_, he'⟩
  have hasq : a ^ 2 = 1 := by
    have hnorm : ‖e‖ = |a| * ‖vEig A hA j‖ := by
      rw [hea, norm_smul, Real.norm_eq_abs]
    rw [hen, hqn, mul_one] at hnorm
    nlinarith [sq_abs a]
  rw [specProjIdx_eq_rankOne hsimple hj w, hea, real_inner_smul_left, smul_smul]
  congr 1
  linear_combination (-⟪vEig A hA j, w⟫_ℝ) * hasq

/-- On the event that the index-`j` eigenvalue of `Xᵀ X` is simple with unit eigenvector `e`,
`perpTopVecIdx v X j k` is a multiple of `perpFrame v e`. Mirror: `perpTopVec_eq`. -/
theorem perpTopVecIdx_eq {q : ℕ} {v : Fin rk → EuclideanSpace ℝ (Fin D)}
    (X : Matrix (Fin q) (Fin D) ℝ) {rk' j : ℕ}
    (hsimple : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) rk') (hj : j < rk')
    {e : EuclideanSpace ℝ (Fin D)}
    (he : e ∈ specSpace (Xᵀ * X)
      (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j)) (hen : ‖e‖ = 1) (k : Fin D) :
    perpTopVecIdx v X j k = (perpFrame v e k) • perpFrame v e := by
  rw [perpTopVecIdx, specProjIdx_eq_rankOne_of_mem hsimple hj he hen, perpFrame_smul]
  congr 1
  rw [← inner_perpFrame_comm, EuclideanSpace.inner_single_right]
  simp

/-- On the simple event `delocDirIdx v X j` is a unit vector parallel to `perpFrame v e`, so it
carries the whole perpendicular part of any inner product with `e`.
Mirror: `abs_inner_perpOf_le`. -/
theorem abs_inner_perpFrame_le {q : ℕ} {v : Fin rk → EuclideanSpace ℝ (Fin D)}
    (hv : ∀ k l, ⟪v k, v l⟫_ℝ = if k = l then 1 else 0)
    (X : Matrix (Fin q) (Fin D) ℝ) {rk' j : ℕ}
    (hsimple : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) rk') (hj : j < rk')
    {e : EuclideanSpace ℝ (Fin D)}
    (he : e ∈ specSpace (Xᵀ * X)
      (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j)) (hen : ‖e‖ = 1)
    (y : EuclideanSpace ℝ (Fin D)) :
    |⟪y, perpFrame v e⟫_ℝ| ≤ |⟪y, delocDirIdx v X j⟫_ℝ| := by
  by_cases hq0 : perpFrame v e = 0
  · rw [hq0]
    simp
  · have hqk : ∃ k, perpFrame v e k ≠ 0 := by
      by_contra hc
      simp only [not_exists, not_not] at hc
      exact hq0 (by apply PiLp.ext; intro k; simpa using hc k)
    obtain ⟨k1, hk1⟩ := hqk
    rcases delocDirIdx_spec v X j with ⟨-, hall⟩ | ⟨k, hkne, hdir⟩
    · exact absurd (hall k1)
        (by rw [perpTopVecIdx_eq X hsimple hj he hen k1]; exact smul_ne_zero hk1 hq0)
    · rw [perpTopVecIdx_eq X hsimple hj he hen k] at hkne hdir
      have hqk0 : perpFrame v e k ≠ 0 := fun h => hkne (by rw [h, zero_smul])
      have hnq : 0 < ‖perpFrame v e‖ := norm_pos_iff.mpr hq0
      have hle1 : ‖perpFrame v e‖ ≤ 1 := hen ▸ norm_perpFrame_le hv e
      rw [hdir, real_inner_smul_right, real_inner_smul_right, norm_smul, Real.norm_eq_abs]
      rw [abs_mul, abs_inv, abs_mul, abs_mul, abs_abs, abs_of_pos hnq]
      rw [inv_mul_eq_div, le_div_iff₀ (mul_pos (abs_pos.mpr hqk0) hnq)]
      nlinarith [mul_nonneg (mul_nonneg (abs_nonneg (⟪y, perpFrame v e⟫_ℝ))
        (abs_nonneg ((perpFrame v e) k))) (sub_nonneg.mpr hle1)]

/-- When `delocDirIdx v X j` is not `0` it is a unit vector orthogonal to the whole frame `v`,
that is, a member of `orthUnitR`. Mirror: `norm_delocDir`. -/
theorem norm_delocDirIdx {q : ℕ} {v : Fin rk → EuclideanSpace ℝ (Fin D)}
    (hv : ∀ k l, ⟪v k, v l⟫_ℝ = if k = l then 1 else 0)
    (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ) (h : delocDirIdx v X j ≠ 0) :
    ‖delocDirIdx v X j‖ = 1 ∧ ∀ l, ⟪delocDirIdx v X j, v l⟫_ℝ = 0 := by
  rcases delocDirIdx_spec v X j with ⟨h0, -⟩ | ⟨k, hkne, hdir⟩
  · exact absurd h0 h
  · refine ⟨?_, fun l => ?_⟩
    · rw [hdir, norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
      exact inv_mul_cancel₀ (norm_ne_zero_iff.mpr hkne)
    · rw [hdir, real_inner_smul_left, perpTopVecIdx, real_inner_comm,
        inner_v_perpFrame hv l, mul_zero]

/-- `overlapIdx X j 0 = 0`. Mirror: `overlap_zero` (`Gram.lean:325`). -/
theorem overlapIdx_apply_zero {q : ℕ} (X : Matrix (Fin q) (Fin D) ℝ) (j : ℕ) :
    overlapIdx X j (0 : EuclideanSpace ℝ (Fin D)) = 0 := by
  change ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) j 0‖ ^ 2 = 0
  rw [map_zero, norm_zero]
  ring

end DelocDirIdx

/-! ### 6. The table: the signal frame, the perpendicular part, and the noise-to-data map -/

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- The columns of `V N` form an orthonormal frame; this is the model field `hV` in the shape
that every `perpFrame` lemma of section 1 takes. -/
theorem inner_col (t : SpikedModelR μ n d rk) (N : ℕ) (k l : Fin rk) :
    ⟪t.col N k, t.col N l⟫_ℝ = if k = l then 1 else 0 := by
  have hVV := t.hV N
  have h : ∑ m, (t.V N) m k * (t.V N) m l = ((t.V N)ᵀ * t.V N) k l := by
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun m _ => by rw [Matrix.transpose_apply]
  rw [real_inner_eq_dotProduct]
  change ∑ m, (t.V N) m k * (t.V N) m l = _
  rw [h, hVV, Matrix.one_apply]

/-- The component of `x` orthogonal to the whole signal span `V R_i` of the table.
Mirror: `perpOf (m.v N)`. -/
noncomputable def perpSpan (t : SpikedModelR μ n d rk) (N : ℕ)
    (x : EuclideanSpace ℝ (Fin (d N))) : EuclideanSpace ℝ (Fin (d N)) :=
  perpFrame (t.col N) x

theorem perpSpan_eq (t : SpikedModelR μ n d rk) (N : ℕ) (x : EuclideanSpace ℝ (Fin (d N))) :
    t.perpSpan N x = x - ∑ k, ⟪t.col N k, x⟫_ℝ • t.col N k := rfl

/-- One table as a function of its noise matrix: `X = U Θ Vᵀ + d^{-1/2} Z`.
Mirror: `SpikedModel.dataOf`. -/
noncomputable def dataOf (t : SpikedModelR μ n d rk) (N : ℕ)
    (Z : Matrix (Fin (n N)) (Fin (d N)) ℝ) : Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  t.U N * Matrix.diagonal t.θ * (t.V N)ᵀ + (Real.sqrt (d N))⁻¹ • Z

/-- `dataOf` at the model noise is the data matrix of that table. Mirror: `dataOf_eq`. -/
theorem dataOf_eq (t : SpikedModelR μ n d rk) (N : ℕ) (ω : Ω N) :
    t.dataOf N (t.Z N ω) = t.X N ω := rfl

/-- The noise-to-data map of one table is measurable. Mirror: `measurable_dataOf`. -/
theorem measurable_dataOf (t : SpikedModelR μ n d rk) (N : ℕ) : Measurable (t.dataOf N) := by
  refine measurable_pi_lambda _ fun a => measurable_pi_lambda _ fun l => ?_
  simp only [SpikedModelR.dataOf, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  have hz : Measurable fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ => Z a l :=
    (measurable_pi_apply l).comp (measurable_pi_apply a)
  exact measurable_const.add (hz.const_mul _)

end SpikedModelR

/-! ### 7. The marginal of the product law and the Fubini step -/

section TableFamilyR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {rk : Fin M → ℕ}

/-- Marginal of the product law: every table gets its own `SpikedModelR.GaussianNoise`.
Mirror: `gaussianNoise_of_joint_pi` (`Gram.lean:394`), with `hasLaw_eval_of_hasLaw_pi` reused
unchanged. -/
theorem gaussianNoise_of_joint_piR (tbl : (i : Fin M) → SpikedModelR μ (n i) d (rk i))
    (hG : ∀ N, HasLaw (fun ω i => (tbl i).Z N ω)
      (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N)) (i : Fin M) :
    (tbl i).GaussianNoise := by
  intro N
  have : ∀ k : Fin M, IsProbabilityMeasure (gaussianMatrix (n k N) (d N)) := fun _ =>
    inferInstance
  exact hasLaw_eval_of_hasLaw_pi (fun k : Fin M => gaussianMatrix (n k N) (d N)) (hG N) i

/-- The Fubini step of `lem:general_rank_delocalization`, model-free. The probability that the
index-`j` singular subspace of table `i` has a large overlap with the random direction
`delocDirIdx (col_i) (X_{i'}) j'` of table `i'` is at most the supremum of that overlap over
deterministic unit directions orthogonal to the whole signal span of table `i`, which
`TableLawR.delocUniform` sends to `0`. Independence across tables enters through `hind`; the
laws `ν N i` of the single tables are arbitrary probability measures, so Gaussian noise is one
case of it. Mirror: `measure_deloc_le_of_pi` (`Gram.lean:407`). -/
theorem measure_deloc_le_of_piR (tbl : (i : Fin M) → SpikedModelR μ (n i) d (rk i))
    {ν : (N : ℕ) → (i : Fin M) → Measure (Matrix (Fin (n i N)) (Fin (d N)) ℝ)}
    [hν : ∀ N i, IsProbabilityMeasure (ν N i)]
    (hind : ∀ N, HasLaw (fun ω i => (tbl i).Z N ω) (Measure.pi (ν N)) (μ N))
    {i i' : Fin M} (hii : i ≠ i') (j j' : ℕ) (N : ℕ) {η : ℝ} (hη : 0 < η) :
    μ N {ω | η ≤ overlapIdx ((tbl i).X N ω) j
        (delocDirIdx ((tbl i).col N) ((tbl i').X N ω) j')}
      ≤ ⨆ w ∈ (tbl i).orthUnitR N, μ N {ω | η ≤ overlapIdx ((tbl i).X N ω) j w} := by
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
      Matrix (Fin (n i' N)) (Fin (d N)) ℝ |
      η ≤ overlapIdx ((tbl i).dataOf N q.1) j
        (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j')} := by
    have h1 : Measurable fun q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i' N)) (Fin (d N)) ℝ =>
        ((tbl i).dataOf N q.1,
          delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j') :=
      (((tbl i).measurable_dataOf N).comp measurable_fst).prodMk
        ((measurable_delocDirIdx ((tbl i).col N) j').comp
          (((tbl i').measurable_dataOf N).comp measurable_snd))
    have hf := (measurable_overlapIdx₂ j).comp h1
    exact hf measurableSet_Ici
  have hpairmeas : MeasurableSet {Zs : (k : Fin M) → Matrix (Fin (n k N)) (Fin (d N)) ℝ |
      (Zs i, Zs i') ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i' N)) (Fin (d N)) ℝ |
        η ≤ overlapIdx ((tbl i).dataOf N q.1) j
          (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j')}} :=
    ((measurable_pi_apply i).prodMk (measurable_pi_apply i')) hSmeas
  have key : μ N {ω | η ≤ overlapIdx ((tbl i).X N ω) j
        (delocDirIdx ((tbl i).col N) ((tbl i').X N ω) j')}
      = (Measure.pi fun k : Fin M => ν N k)
        {Zs | (Zs i, Zs i') ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n i' N)) (Fin (d N)) ℝ |
          η ≤ overlapIdx ((tbl i).dataOf N q.1) j
            (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j')}} :=
    (hind N).measure_eq hpairmeas
  rw [key]
  refine pi_pair_le (fun k : Fin M => ν N k) hprob hii hSmeas _ ?_
  intro zj
  by_cases hw : delocDirIdx ((tbl i).col N) ((tbl i').dataOf N zj) j' = 0
  · have hempty : {zi : Matrix (Fin (n i N)) (Fin (d N)) ℝ |
        (zi, zj) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n i' N)) (Fin (d N)) ℝ |
          η ≤ overlapIdx ((tbl i).dataOf N q.1) j
            (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j')}}
        = (∅ : Set (Matrix (Fin (n i N)) (Fin (d N)) ℝ)) := by
      ext zi
      simp only [Set.mem_ofPred_eq, hw, overlapIdx_apply_zero, Set.mem_empty_iff_false,
        iff_false, not_le]
      exact hη
    rw [hempty]
    simp
  · have hmem : delocDirIdx ((tbl i).col N) ((tbl i').dataOf N zj) j'
        ∈ (tbl i).orthUnitR N :=
      norm_delocDirIdx ((tbl i).inner_col N) _ j' hw
    have hsl : (ν N i) {zi : Matrix (Fin (n i N)) (Fin (d N)) ℝ |
        (zi, zj) ∈ {q : Matrix (Fin (n i N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n i' N)) (Fin (d N)) ℝ |
          η ≤ overlapIdx ((tbl i).dataOf N q.1) j
            (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N q.2) j')}}
        = μ N {ω | η ≤ overlapIdx ((tbl i).X N ω) j
          (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N zj) j')} := by
      refine (hmarg i (fun z => η ≤ overlapIdx ((tbl i).dataOf N z) j
        (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N zj) j')) ?_).symm
      have hf := (measurable_overlapIdx j
        (delocDirIdx ((tbl i).col N) ((tbl i').dataOf N zj) j')).comp
          ((tbl i).measurable_dataOf N)
      exact hf measurableSet_Ici
    rw [hsl]
    exact le_iSup₂ (f := fun w (_ : w ∈ (tbl i).orthUnitR N) =>
      μ N {ω | η ≤ overlapIdx ((tbl i).X N ω) j w}) _ hmem

end TableFamilyR

/-! ### 8. Independent noise across tables: the paper's `assum:general_noise` -/

section IndepNoiseR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- `assum:general_noise` as the paper states it: the tables are independent, each with some
probability law on its noise matrix. `JointGaussianNoise` is the Gaussian case. Every Layer 1
theorem of Track B reads the noise only through this predicate, because the single place that
needs it is the Fubini step `measure_deloc_le_of_piR`. -/
def UnalignedModelR.IndepNoise (m : UnalignedModelR μ M n d r rk) : Prop :=
  ∃ ν : (N : ℕ) → (i : Fin M) → Measure (Matrix (Fin (n i N)) (Fin (d N)) ℝ),
    (∀ N i, IsProbabilityMeasure (ν N i)) ∧
      ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω) (Measure.pi (ν N)) (μ N)

/-- Independent Gaussian tables are independent tables. -/
theorem UnalignedModelR.JointGaussianNoise.indepNoise {m : UnalignedModelR μ M n d r rk}
    (hG : m.JointGaussianNoise) : m.IndepNoise :=
  ⟨fun N i => gaussianMatrix (n i N) (d N), fun _ _ => inferInstance, hG⟩

end IndepNoiseR

end StackedSVD
