/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Spectral
import StackedSVD.Prob.TendstoInProb
import StatsMLlib.LinearAlgebra.Matrix.Perturbation

/-!
# Weyl and Davis-Kahan for `topProj`: the deterministic core of step (iv)

This file has the deterministic perturbation facts that `lem_entrywise_conv_eigenvec` and
`thm_svd_stack_general` of `SVDStack.lean` need. Nothing here is probabilistic.

The matrix norm is the `L2` operator norm of `Mathlib.Analysis.CStarAlgebra.Matrix`
(`open scoped Matrix.Norms.L2Operator`), which is `‖toEuclideanCLM A‖` by definition.

## Content

1. `abs_eigenvalues₀_sub_le`, `lamMax_lipschitz`, `abs_eigenvalues₀_one_sub_le`: Weyl.
   Every sorted eigenvalue is `1`-Lipschitz in the operator norm.
2. `topSimple_of_gap`: a positive gap `λ₁ - λ₂ > 0` gives `TopSimple`.
3. `topProj_apply_of_topSimple`, `norm_topProj_sq_eq_inner_sq`: the plain-matrix form of
   `overlap_eq_inner_sq` of `Spectral.lean` (that one is stated for a Gram matrix).
4. `topProj_perturb`: Davis-Kahan in projector form. A gap `γ` for `A` and `‖B - A‖ ≤ δ`
   with `2 δ < γ` give `TopSimple B` and `‖topProj B - topProj A‖ ≤ 4 δ / γ`.
5. `l2_opNorm_le_of_entries`: `‖C‖ ≤ d * max |C i j|`, the bridge from entrywise
   convergence to operator-norm convergence.
6. `norm_topProj_sq_continuous_at` (operator norm) and `norm_topProj_sq_continuousAt`
   (entrywise, a `ContinuousAt` on `Fin d × Fin d → ℝ`): `‖topProj B x‖²` is continuous in
   `B` at a matrix with a simple top eigenvalue. The second form is the one
   `TendstoInProbPi.comp_continuous` takes.
7. `opNorm_sub_tendstoInProb`: entrywise convergence in probability of a family of symmetric
   random matrices gives `‖G_N - A‖ → 0` in probability. It came from
   `SVDStack/Deterministic.lean` on 2026-08-30 so that `LinAlg/SpecProjPerturb.lean`, which
   does not import `SVDStack`, can use it. This is the one probabilistic statement in the
   file, and the reason for the `StackedSVD.Prob.TendstoInProb` import.

## The constant of `topProj_perturb`

StatsMLlib delivers `davisKahan_eigenvector_projection_of_gap_le`, that is
`‖P_{(span v)ᗮ} u‖ ≤ 2 ‖A - B‖ / γ` with `u`, `v` the unit top eigenvectors and `γ` the gap
of `A`. Write `z = u - ⟪v, u⟫ v`, so `‖z‖ = ‖P_{(span v)ᗮ} u‖`. Two facts finish it:
`‖v - ⟪v, u⟫ u‖ = ‖z‖` (both squares equal `1 - ⟪v, u⟫²`), and the identity
`P_B x - P_A x = ⟪v, x⟫ (v - ⟪v, u⟫ u) - ⟪z, x⟫ u`. The triangle inequality then gives
`‖P_B - P_A‖ ≤ 2 ‖z‖ ≤ 4 δ / γ`. So the constant is `4`, not the `2` of the sine-theta bound.
-/

open Filter Topology
open scoped InnerProductSpace Matrix Matrix.Norms.L2Operator

namespace StackedSVD

variable {d : ℕ}

/-! ### Bridges to StatsMLlib and to the eigenbasis -/

section Bridge

/-- `Matrix.IsHermitian.eigenvalues₀` is the operator eigenvalue list of StatsMLlib. -/
private theorem eigenvalues₀_eq_symm {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (i : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ i = (symmOp hA).eigenvalues finrank_euclideanSpace i := rfl

/-- The `L2` operator norm of a matrix difference is the operator norm of the difference of
the two operators. -/
private theorem opNorm_toCLM_sub (A B : Matrix (Fin d) (Fin d) ℝ) :
    ‖(toOp A - toOp B).toContinuousLinearMap‖ = ‖A - B‖ := by
  rw [show toOp A - toOp B = toOp (A - B) from (map_sub Matrix.toEuclideanLin A B).symm]
  rfl

private theorem pos_card (h1 : 1 < Fintype.card (Fin d)) : 0 < Fintype.card (Fin d) :=
  lt_trans Nat.zero_lt_one h1

private theorem pos_dim (hcard : 0 < Fintype.card (Fin d)) : 0 < d := by simpa using hcard

/-- `lamMax` is the sorted eigenvalue of index `0`. -/
private theorem lamMax_eq_eigenvalues₀ {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hcard : 0 < Fintype.card (Fin d)) : lamMax A hA = hA.eigenvalues₀ ⟨0, hcard⟩ := by
  rw [lamMax, dif_pos (pos_dim hcard)]

/-- The top eigenvector, in StatsMLlib's index convention. -/
private noncomputable def topEigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hcard : 0 < Fintype.card (Fin d)) : EuclideanSpace ℝ (Fin d) :=
  (symmOp hA).eigenvectorBasis finrank_euclideanSpace ⟨0, hcard⟩

private theorem norm_topEigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hcard : 0 < Fintype.card (Fin d)) : ‖topEigvec hA hcard‖ = 1 :=
  ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).norm_eq_one _

/-- The top eigenspace is the eigenspace at `lamMax`. Public since cleanup wave 2
(2026-08-30): `RankR/Defs.lean` re-derived it as `topSpace_eq_eigenspace'`. -/
theorem topSpace_eq_eigenspace (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    topSpace A hA = Module.End.eigenspace (toOp A) (lamMax A hA) := by
  simp [topSpace, specSpace]

private theorem apply_topEigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hcard : 0 < Fintype.card (Fin d)) :
    toOp A (topEigvec hA hcard) = lamMax A hA • topEigvec hA hcard := by
  rw [lamMax_eq_eigenvalues₀ hA hcard]
  exact apply_eigvec hA ⟨0, hcard⟩

private theorem mem_topSpace_topEigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hcard : 0 < Fintype.card (Fin d)) : topEigvec hA hcard ∈ topSpace A hA := by
  rw [topSpace_eq_eigenspace]
  exact Module.End.mem_eigenspace_iff.mpr (apply_topEigvec hA hcard)

/-- Below the top index every eigenvalue is at most the second one. -/
private theorem eigenvalues₀_le_second {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d)) {i : Fin (Fintype.card (Fin d))}
    (hi : i ≠ ⟨0, pos_card h1⟩) : hA.eigenvalues₀ i ≤ hA.eigenvalues₀ ⟨1, h1⟩ := by
  have hival : 1 ≤ i.val := by
    have hne : i.val ≠ 0 := fun h => hi (Fin.ext h)
    omega
  exact hA.eigenvalues₀_antitone (by rw [Fin.le_def]; exact hival)

end Bridge

/-! ### 1. Weyl: every sorted eigenvalue is 1-Lipschitz -/

/-- Weyl's inequality, in the index convention of `Matrix.IsHermitian.eigenvalues₀`. -/
theorem abs_eigenvalues₀_sub_le (A B : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (i : Fin (Fintype.card (Fin d))) :
    |hA.eigenvalues₀ i - hB.eigenvalues₀ i| ≤ ‖A - B‖ := by
  have key := (symmOp hA).abs_eigenvalues_sub_le_opNorm (symmOp hB) finrank_euclideanSpace i
  rwa [opNorm_toCLM_sub] at key

/-- `lamMax` is 1-Lipschitz for the `L2` operator norm. -/
theorem lamMax_lipschitz (A B : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) : |lamMax A hA - lamMax B hB| ≤ ‖A - B‖ := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    rw [lamMax, lamMax, dif_neg (by omega), dif_neg (by omega)]
    simp
  · have hcard : 0 < Fintype.card (Fin d) := by simpa using hd
    rw [lamMax_eq_eigenvalues₀ hA hcard, lamMax_eq_eigenvalues₀ hB hcard]
    exact abs_eigenvalues₀_sub_le A B hA hB ⟨0, hcard⟩

/-- The second eigenvalue is 1-Lipschitz too. -/
theorem abs_eigenvalues₀_one_sub_le (A B : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h1 : 1 < Fintype.card (Fin d)) :
    |hA.eigenvalues₀ ⟨1, h1⟩ - hB.eigenvalues₀ ⟨1, h1⟩| ≤ ‖A - B‖ :=
  abs_eigenvalues₀_sub_le A B hA hB ⟨1, h1⟩

/-! ### 2. A positive gap makes the top eigenvalue simple -/

/-- With a positive gap the top eigenspace is the line of the top eigenvector. -/
private theorem topSpace_eq_span_topEigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d))
    (hgap : 0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) :
    topSpace A hA = ℝ ∙ topEigvec hA (pos_card h1) := by
  have hc : 0 < Fintype.card (Fin d) := pos_card h1
  refine le_antisymm ?_ ?_
  · intro x hx
    have hxe : toOp A x = lamMax A hA • x :=
      Module.End.mem_eigenspace_iff.mp (by rwa [topSpace_eq_eigenspace] at hx)
    have hzero : ∀ i : Fin (Fintype.card (Fin d)), i ≠ ⟨0, hc⟩ →
        ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr x i = 0 := by
      intro i hi
      have hle : hA.eigenvalues₀ i ≤ hA.eigenvalues₀ ⟨1, h1⟩ :=
        eigenvalues₀_le_second hA h1 hi
      have hne : hA.eigenvalues₀ i ≠ lamMax A hA := by
        intro h
        rw [h] at hle
        linarith
      rw [OrthonormalBasis.repr_apply_apply]
      exact inner_eq_zero_of_ne hA hne (apply_eigvec hA i) hxe
    have hexp : ∑ i, ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr x i •
          ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) i
        = ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr x ⟨0, hc⟩ •
          ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) ⟨0, hc⟩ :=
      Finset.sum_eq_single _ (fun i _ hi => by rw [hzero i hi, zero_smul])
        (fun h => absurd (Finset.mem_univ _) h)
    have hx2 : x = ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr x ⟨0, hc⟩ •
        topEigvec hA hc := by
      calc x = ∑ i, ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr x i •
            ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) i :=
          (((symmOp hA).eigenvectorBasis finrank_euclideanSpace).sum_repr x).symm
        _ = _ := hexp
    rw [hx2]
    exact Submodule.smul_mem _ _ (Submodule.mem_span_singleton_self _)
  · rw [Submodule.span_singleton_le_iff_mem]
    exact mem_topSpace_topEigvec hA hc

/-- A positive gap `λ₁ - λ₂ > 0` gives a simple top eigenvalue. -/
theorem topSimple_of_gap {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d))
    (hgap : 0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) : TopSimple A hA := by
  have hne : topEigvec hA (pos_card h1) ≠ 0 := by
    intro h
    have hn := norm_topEigvec hA (pos_card h1)
    rw [h, norm_zero] at hn
    exact zero_ne_one hn
  rw [TopSimple, topSpace_eq_span_topEigvec hA h1 hgap, finrank_span_singleton hne]

/-! ### 3. The rank-one form of `topProj` -/

/-- Under simplicity the top projector is the rank-one projector on a unit top eigenvector.
Same statement as the private `topProj_eq_rankOne` of `SVDStack.lean`, under a different name
so that the two do not clash when `SVDStack.lean` imports this file. -/
theorem topProj_apply_of_topSimple {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian}
    (hsimple : TopSimple A hA) {e : EuclideanSpace ℝ (Fin d)} (he : e ∈ topSpace A hA)
    (hne : ‖e‖ = 1) (w : EuclideanSpace ℝ (Fin d)) : topProj A hA w = ⟪e, w⟫_ℝ • e := by
  have he0 : e ≠ 0 := by
    intro h
    rw [h, norm_zero] at hne
    exact zero_ne_one hne
  have hle : (ℝ ∙ e) ≤ topSpace A hA := (Submodule.span_singleton_le_iff_mem _ _).2 he
  have hrank : Module.finrank ℝ (ℝ ∙ e) = Module.finrank ℝ (topSpace A hA) := by
    rw [finrank_span_singleton he0, hsimple]
  have hspan : (ℝ ∙ e) = topSpace A hA := Submodule.eq_of_le_of_finrank_eq hle hrank
  change (topSpace A hA).starProjection w = _
  rw [← hspan, Submodule.starProjection_unit_singleton ℝ hne]

/-- The plain-matrix form of `overlap_eq_inner_sq` (which is stated for a Gram matrix):
`‖topProj A x‖² = ⟪e, x⟫²` for a unit top eigenvector `e`. -/
theorem norm_topProj_sq_eq_inner_sq {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian}
    (hsimple : TopSimple A hA) {e : EuclideanSpace ℝ (Fin d)} (he : e ∈ topSpace A hA)
    (hne : ‖e‖ = 1) (w : EuclideanSpace ℝ (Fin d)) :
    ‖topProj A hA w‖ ^ 2 = ⟪e, w⟫_ℝ ^ 2 := by
  rw [topProj_apply_of_topSimple hsimple he hne, norm_smul, hne, mul_one, Real.norm_eq_abs, sq_abs]

/-- `⟪P w, w⟫ = ‖P w‖²` for the top projector. -/
theorem inner_topProj_self_eq_norm_sq (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (w : EuclideanSpace ℝ (Fin d)) : ⟪topProj A hA w, w⟫_ℝ = ‖topProj A hA w‖ ^ 2 := by
  have h0 : ⟪(topSpace A hA).starProjection w, w - (topSpace A hA).starProjection w⟫_ℝ = 0 :=
    (topSpace A hA).sub_starProjection_mem_orthogonal w _
      ((topSpace A hA).starProjection_apply_mem w)
  have h1 : ⟪(topSpace A hA).starProjection w, w⟫_ℝ -
      ⟪(topSpace A hA).starProjection w, (topSpace A hA).starProjection w⟫_ℝ = 0 := by
    rw [← inner_sub_right]
    exact h0
  change ⟪(topSpace A hA).starProjection w, w⟫_ℝ = ‖(topSpace A hA).starProjection w‖ ^ 2
  rw [← real_inner_self_eq_norm_sq]
  linarith

/-! ### 4. Davis-Kahan in projector form -/

/-- Davis-Kahan for the top spectral projector. If `A` has gap `γ > 0` between its first and
second eigenvalue, and `‖B - A‖ ≤ δ` with `2 δ < γ`, then both top eigenvalues are simple and
the two projectors are `4 δ / γ` apart in operator norm. -/
theorem topProj_perturb {A B : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h1 : 1 < Fintype.card (Fin d)) {γ δ : ℝ} (hγ : 0 < γ)
    (hgapA : γ ≤ lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) (hδ : ‖B - A‖ ≤ δ) (hlt : 2 * δ < γ) :
    TopSimple A hA ∧ TopSimple B hB ∧ ‖topProj B hB - topProj A hA‖ ≤ 4 * δ / γ := by
  have hc : 0 < Fintype.card (Fin d) := pos_card h1
  have hnorm : ‖A - B‖ ≤ δ := by rwa [norm_sub_rev]
  have hδ0 : 0 ≤ δ := le_trans (norm_nonneg _) hδ
  have hsA : TopSimple A hA := topSimple_of_gap hA h1 (by linarith)
  -- the gap of `B`, by Weyl at the indices `0` and `1`
  have hw0 := abs_le.mp (lamMax_lipschitz A B hA hB)
  have hw1 := abs_le.mp (abs_eigenvalues₀_one_sub_le A B hA hB h1)
  have hgapB : 0 < lamMax B hB - hB.eigenvalues₀ ⟨1, h1⟩ := by
    have h01 := hw0.1
    have h02 := hw0.2
    have h11 := hw1.1
    have h12 := hw1.2
    linarith
  have hsB : TopSimple B hB := topSimple_of_gap hB h1 hgapB
  refine ⟨hsA, hsB, ?_⟩
  -- the gap hypothesis of Davis-Kahan
  have hgapDK : ∀ i : Fin (Fintype.card (Fin d)), i ≠ ⟨0, hc⟩ →
      γ ≤ |hA.eigenvalues₀ ⟨0, hc⟩ - hA.eigenvalues₀ i| := by
    intro i hi
    have hle : hA.eigenvalues₀ i ≤ hA.eigenvalues₀ ⟨1, h1⟩ := eigenvalues₀_le_second hA h1 hi
    have hlam : lamMax A hA = hA.eigenvalues₀ ⟨0, hc⟩ := lamMax_eq_eigenvalues₀ hA hc
    rw [abs_of_nonneg (by linarith)]
    linarith
  -- the two unit top eigenvectors, with the Davis-Kahan bound
  obtain ⟨u, v, hun, hvn, humem, hvmem, hDK⟩ :
      ∃ u v : EuclideanSpace ℝ (Fin d), ‖u‖ = 1 ∧ ‖v‖ = 1 ∧ u ∈ topSpace A hA ∧
        v ∈ topSpace B hB ∧ ‖((ℝ ∙ v)ᗮ).starProjection u‖ ≤ 2 * ‖A - B‖ / γ :=
    ⟨topEigvec hA hc, topEigvec hB hc, norm_topEigvec hA hc, norm_topEigvec hB hc,
      mem_topSpace_topEigvec hA hc, mem_topSpace_topEigvec hB hc, by
        have h := (symmOp hA).davisKahan_eigenvector_projection_of_gap_le (symmOp hB)
          finrank_euclideanSpace ⟨0, hc⟩ hγ hgapDK
        rwa [opNorm_toCLM_sub] at h⟩
  set a := ⟪v, u⟫_ℝ with ha
  set z := u - a • v with hz
  have hzproj : ((ℝ ∙ v)ᗮ).starProjection u = z := by
    rw [Submodule.starProjection_orthogonal_val, Submodule.starProjection_unit_singleton ℝ hvn,
      ← ha, hz]
  have hzbound : ‖z‖ ≤ 2 * δ / γ := by
    rw [hzproj] at hDK
    refine le_trans hDK ?_
    rw [div_eq_mul_inv, div_eq_mul_inv]
    have hinv : (0 : ℝ) ≤ γ⁻¹ := le_of_lt (inv_pos.mpr hγ)
    nlinarith [hnorm]
  -- the mirror vector has the same norm
  have huv : ⟪u, v⟫_ℝ = a := by rw [ha, real_inner_comm]
  have hsq1 : ‖z‖ ^ 2 = 1 - a ^ 2 := by
    rw [hz, norm_sub_sq_real, real_inner_smul_right, norm_smul, hun, huv, Real.norm_eq_abs, hvn]
    nlinarith [sq_abs a]
  have hsq2 : ‖v - a • u‖ ^ 2 = 1 - a ^ 2 := by
    rw [norm_sub_sq_real, real_inner_smul_right, norm_smul, hun, hvn, Real.norm_eq_abs, ← ha]
    nlinarith [sq_abs a]
  have hmirror : ‖v - a • u‖ = ‖z‖ := by
    rw [← Real.sqrt_sq (norm_nonneg (v - a • u)), ← Real.sqrt_sq (norm_nonneg z), hsq1, hsq2]
  -- the operator bound
  have hop : ‖topProj B hB - topProj A hA‖ ≤ 2 * ‖z‖ := by
    refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun x => ?_
    have hzx : ⟪z, x⟫_ℝ = ⟪u, x⟫_ℝ - a * ⟪v, x⟫_ℝ := by
      rw [hz, inner_sub_left, real_inner_smul_left]
    have hsub : (topProj B hB - topProj A hA) x = topProj B hB x - topProj A hA x := rfl
    have hkey : (topProj B hB - topProj A hA) x
        = ⟪v, x⟫_ℝ • (v - a • u) - ⟪z, x⟫_ℝ • u := by
      rw [hsub, topProj_apply_of_topSimple hsB hvmem hvn,
        topProj_apply_of_topSimple hsA humem hun, hzx]
      module
    rw [hkey]
    have hb1 : |⟪v, x⟫_ℝ| ≤ ‖x‖ := by
      have h := abs_real_inner_le_norm v x
      rwa [hvn, one_mul] at h
    have hb2 : |⟪z, x⟫_ℝ| ≤ ‖z‖ * ‖x‖ := abs_real_inner_le_norm z x
    calc ‖⟪v, x⟫_ℝ • (v - a • u) - ⟪z, x⟫_ℝ • u‖
        ≤ ‖⟪v, x⟫_ℝ • (v - a • u)‖ + ‖⟪z, x⟫_ℝ • u‖ := norm_sub_le _ _
      _ = |⟪v, x⟫_ℝ| * ‖z‖ + |⟪z, x⟫_ℝ| * 1 := by
          rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs, hmirror, hun]
      _ ≤ 2 * ‖z‖ * ‖x‖ := by
          nlinarith [mul_le_mul_of_nonneg_right hb1 (norm_nonneg z), hb2]
  have hfin : (4 : ℝ) * δ / γ = 2 * (2 * δ / γ) := by ring
  rw [hfin]
  linarith

/-! ### 5. From entrywise size to operator norm -/

private theorem sum_sq_mulVec_le (C : Matrix (Fin d) (Fin d) ℝ) {K : ℝ}
    (hK : ∀ i j, |C i j| ≤ K) (v : Fin d → ℝ) :
    ∑ i, ((C *ᵥ v) i) ^ 2 ≤ ((d : ℝ) * K) ^ 2 * ∑ j, (v j) ^ 2 := by
  have hcs : (∑ j, |v j|) ^ 2 ≤ (d : ℝ) * ∑ j, (v j) ^ 2 := by
    have h := sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset (Fin d)))
      (f := fun j => |v j|)
    have hcard : ((Finset.univ : Finset (Fin d)).card : ℝ) = (d : ℝ) := by simp
    rw [hcard] at h
    refine h.trans (le_of_eq ?_)
    congr 1
    exact Finset.sum_congr rfl fun j _ => sq_abs (v j)
  have hnn : 0 ≤ ∑ j, |v j| := Finset.sum_nonneg fun j _ => abs_nonneg _
  have hrow : ∀ i, ((C *ᵥ v) i) ^ 2 ≤ K ^ 2 * ((d : ℝ) * ∑ j, (v j) ^ 2) := by
    intro i
    have hK0 : 0 ≤ K := le_trans (abs_nonneg _) (hK i i)
    have hmv : (C *ᵥ v) i = ∑ j, C i j * v j := by simp [Matrix.mulVec, dotProduct]
    have habs : |(C *ᵥ v) i| ≤ K * ∑ j, |v j| := by
      rw [hmv]
      calc |∑ j, C i j * v j| ≤ ∑ j, |C i j * v j| := Finset.abs_sum_le_sum_abs _ _
        _ ≤ ∑ j, K * |v j| := by
            refine Finset.sum_le_sum fun j _ => ?_
            rw [abs_mul]
            exact mul_le_mul_of_nonneg_right (hK i j) (abs_nonneg _)
        _ = K * ∑ j, |v j| := by rw [Finset.mul_sum]
    have h2 : ((C *ᵥ v) i) ^ 2 ≤ (K * ∑ j, |v j|) ^ 2 := by
      nlinarith [abs_nonneg ((C *ᵥ v) i), sq_abs ((C *ᵥ v) i), habs, mul_nonneg hK0 hnn]
    calc ((C *ᵥ v) i) ^ 2 ≤ (K * ∑ j, |v j|) ^ 2 := h2
      _ = K ^ 2 * (∑ j, |v j|) ^ 2 := by ring
      _ ≤ K ^ 2 * ((d : ℝ) * ∑ j, (v j) ^ 2) := mul_le_mul_of_nonneg_left hcs (by positivity)
  calc ∑ i, ((C *ᵥ v) i) ^ 2
      ≤ ∑ _i : Fin d, K ^ 2 * ((d : ℝ) * ∑ j, (v j) ^ 2) := Finset.sum_le_sum fun i _ => hrow i
    _ = (d : ℝ) * (K ^ 2 * ((d : ℝ) * ∑ j, (v j) ^ 2)) := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    _ = ((d : ℝ) * K) ^ 2 * ∑ j, (v j) ^ 2 := by ring

/-- The `L2` operator norm is at most `d` times the largest entry in absolute value. The
constant is not sharp; it only has to be finite. -/
theorem l2_opNorm_le_of_entries (C : Matrix (Fin d) (Fin d) ℝ) {K : ℝ}
    (hK : ∀ i j, |C i j| ≤ K) : ‖C‖ ≤ d * K := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    have hC : C = 0 := by
      ext i j
      exact i.elim0
    simp [hC]
  have hK0 : 0 ≤ K := le_trans (abs_nonneg _) (hK ⟨0, hd⟩ ⟨0, hd⟩)
  rw [Matrix.cstar_norm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun x => ?_
  have hxsq : ‖x‖ ^ 2 = ∑ j, (WithLp.ofLp x j) ^ 2 := EuclideanSpace.real_norm_sq_eq x
  have hysq : ‖(Matrix.toEuclideanCLM (𝕜 := ℝ) C) x‖ ^ 2
      = ∑ i, ((C *ᵥ WithLp.ofLp x) i) ^ 2 := EuclideanSpace.real_norm_sq_eq _
  have hmain := sum_sq_mulVec_le C hK (WithLp.ofLp x)
  rw [← hxsq, ← hysq] at hmain
  nlinarith [norm_nonneg ((Matrix.toEuclideanCLM (𝕜 := ℝ) C) x), norm_nonneg x, hK0,
    mul_nonneg (mul_nonneg (Nat.cast_nonneg (α := ℝ) d) hK0) (norm_nonneg x)]

/-! ### 6. Continuity of `‖topProj B x‖²` at a matrix with a simple top eigenvalue -/

/-- Operator-norm form: the `ε`-`δ` statement that the continuous mapping theorem for
convergence in probability needs. -/
theorem norm_topProj_sq_continuous_at {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d)) (hgap : 0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩)
    (x : EuclideanSpace ℝ (Fin d)) {ε : ℝ} (hε : 0 < ε) :
    ∃ δ > 0, ∀ (B : Matrix (Fin d) (Fin d) ℝ) (hB : B.IsHermitian), ‖B - A‖ < δ →
      |‖topProj B hB x‖ ^ 2 - ‖topProj A hA x‖ ^ 2| < ε := by
  set γ := lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩ with hγdef
  refine ⟨min (γ / 3) (ε * γ / (4 * (‖x‖ ^ 2 + 1))), by positivity, ?_⟩
  intro B hB hBA
  have hlt1 : ‖B - A‖ < γ / 3 := lt_of_lt_of_le hBA (min_le_left _ _)
  have hlt2 : ‖B - A‖ < ε * γ / (4 * (‖x‖ ^ 2 + 1)) := lt_of_lt_of_le hBA (min_le_right _ _)
  have hd0 : 0 ≤ ‖B - A‖ := norm_nonneg _
  have hx2 : (0 : ℝ) ≤ ‖x‖ ^ 2 := by positivity
  obtain ⟨-, -, hproj⟩ :=
    topProj_perturb hA hB h1 hgap le_rfl (le_refl ‖B - A‖) (by linarith)
  -- the projector bound gives a bound on the two quadratic forms
  have hquad : |‖topProj B hB x‖ ^ 2 - ‖topProj A hA x‖ ^ 2|
      ≤ ‖topProj B hB - topProj A hA‖ * ‖x‖ ^ 2 := by
    have hsub : (topProj B hB - topProj A hA) x = topProj B hB x - topProj A hA x := rfl
    have hdiff : ‖topProj B hB x‖ ^ 2 - ‖topProj A hA x‖ ^ 2
        = ⟪(topProj B hB - topProj A hA) x, x⟫_ℝ := by
      rw [hsub, inner_sub_left, inner_topProj_self_eq_norm_sq, inner_topProj_self_eq_norm_sq]
    rw [hdiff]
    have hb1 : |⟪(topProj B hB - topProj A hA) x, x⟫_ℝ|
        ≤ ‖(topProj B hB - topProj A hA) x‖ * ‖x‖ := abs_real_inner_le_norm _ _
    have hb2 : ‖(topProj B hB - topProj A hA) x‖ ≤ ‖topProj B hB - topProj A hA‖ * ‖x‖ :=
      ContinuousLinearMap.le_opNorm _ _
    nlinarith [mul_le_mul_of_nonneg_right hb2 (norm_nonneg x), hb1]
  -- and that bound is below `ε`
  have hkey : 4 * ‖B - A‖ * (‖x‖ ^ 2 + 1) < ε * γ := by
    have hpos : (0 : ℝ) < 4 * (‖x‖ ^ 2 + 1) := by positivity
    calc 4 * ‖B - A‖ * (‖x‖ ^ 2 + 1) = ‖B - A‖ * (4 * (‖x‖ ^ 2 + 1)) := by ring
      _ < ε * γ / (4 * (‖x‖ ^ 2 + 1)) * (4 * (‖x‖ ^ 2 + 1)) := mul_lt_mul_of_pos_right hlt2 hpos
      _ = ε * γ := by field_simp
  have hmul : 4 * ‖B - A‖ / γ * ‖x‖ ^ 2 * γ < ε * γ := by
    have hrw : 4 * ‖B - A‖ / γ * ‖x‖ ^ 2 * γ = 4 * ‖B - A‖ * ‖x‖ ^ 2 := by
      field_simp
    rw [hrw]
    nlinarith [hkey, hd0]
  have hfrac : 4 * ‖B - A‖ / γ * ‖x‖ ^ 2 < ε :=
    lt_of_mul_lt_mul_right hmul (le_of_lt hgap)
  have hstep : ‖topProj B hB - topProj A hA‖ * ‖x‖ ^ 2 ≤ 4 * ‖B - A‖ / γ * ‖x‖ ^ 2 :=
    mul_le_mul_of_nonneg_right hproj hx2
  linarith

/-! ### 6b. The same statement entrywise, as a `ContinuousAt` -/

/-- Symmetrization of a function of index pairs, read as a matrix. It is symmetric for every
input, so `topProj` applies to it with no side condition. -/
noncomputable def symMat (w : Fin d × Fin d → ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  Matrix.of fun i j => (w (i, j) + w (j, i)) / 2

theorem isHermitian_symMat (w : Fin d × Fin d → ℝ) : (symMat w).IsHermitian := by
  change (symMat w)ᴴ = symMat w
  ext i j
  simp only [Matrix.conjTranspose_apply, symMat, Matrix.of_apply, star_trivial]
  ring

theorem symMat_entries {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) :
    symMat (fun p => A p.1 p.2) = A := by
  ext i j
  have h : A j i = A i j := by simpa using hA.apply i j
  simp only [symMat, Matrix.of_apply, h]
  ring

/-- Value of the symmetrized map at a symmetric matrix. -/
theorem norm_topProj_sq_symMat {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (x : EuclideanSpace ℝ (Fin d)) :
    ‖topProj (symMat fun p => A p.1 p.2) (isHermitian_symMat _) x‖ ^ 2
      = ‖topProj A hA x‖ ^ 2 := by
  rw [topProj_congr (symMat_entries hA) _ hA x]

/-- Entrywise continuity of `‖topProj B x‖²` at a matrix with a simple top eigenvalue, in the
form `TendstoInProbPi.comp_continuous` takes: a `ContinuousAt` for a function of the
`Fin d × Fin d` family of entries. -/
theorem norm_topProj_sq_continuousAt {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d)) (hgap : 0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩)
    (x : EuclideanSpace ℝ (Fin d)) :
    ContinuousAt (fun w : Fin d × Fin d → ℝ =>
      ‖topProj (symMat w) (isHermitian_symMat w) x‖ ^ 2) (fun p => A p.1 p.2) := by
  refine Metric.continuousAt_iff.mpr fun ε hε => ?_
  obtain ⟨δ, hδ0, hδ⟩ := norm_topProj_sq_continuous_at hA h1 hgap x hε
  have hdpos : (0 : ℝ) < (d : ℝ) + 1 := by positivity
  refine ⟨δ / ((d : ℝ) + 1), by positivity, fun {w} hw => ?_⟩
  have hAji : ∀ i j, A j i = A i j := fun i j => by simpa using hA.apply i j
  have hentry : ∀ i j, |(symMat w - A) i j| ≤ δ / ((d : ℝ) + 1) := by
    intro i j
    have hij := (dist_pi_lt_iff (show (0 : ℝ) < δ / ((d : ℝ) + 1) by positivity)).mp hw (i, j)
    have hji := (dist_pi_lt_iff (show (0 : ℝ) < δ / ((d : ℝ) + 1) by positivity)).mp hw (j, i)
    simp only [Real.dist_eq] at hij hji
    have hval : (symMat w - A) i j
        = ((w (i, j) - A i j) + (w (j, i) - A i j)) / 2 := by
      simp only [Matrix.sub_apply, symMat, Matrix.of_apply]
      ring
    rw [hval, abs_div, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 2)]
    have habs : |w (i, j) - A i j + (w (j, i) - A i j)|
        ≤ |w (i, j) - A i j| + |w (j, i) - A i j| := abs_add_le _ _
    rw [hAji i j] at hji
    linarith
  have hop : ‖symMat w - A‖ < δ := by
    have h := l2_opNorm_le_of_entries (symMat w - A) hentry
    have hlt : (d : ℝ) * (δ / ((d : ℝ) + 1)) < δ := by
      have heq : (d : ℝ) * (δ / ((d : ℝ) + 1)) = δ - δ / ((d : ℝ) + 1) := by
        field_simp
        ring
      rw [heq]
      linarith [div_pos hδ0 hdpos]
    exact lt_of_le_of_lt h hlt
  have hres := hδ (symMat w) (isHermitian_symMat w) hop
  rw [Real.dist_eq, norm_topProj_sq_symMat hA x]
  exact hres

/-! ### 7. Additions for the assembly of `thm:svd_stack_general` (2026-08-29)

Appended after the first version of the file. Content:

* `finrank_eigenspace_transpose_mul_eq`, `topSimple_transpose_mul_iff`: for `t ≠ 0` the
  eigenspaces of `Aᵀ A` and `A Aᵀ` at `t` have the same dimension (the maps `x ↦ A x` and
  `y ↦ Aᵀ y` are injective between them), so simplicity of the top eigenvalue transfers
  between the two Gram matrices when that eigenvalue is positive and shared.
* `gap_pos_of_topSimple`: the converse of `topSimple_of_gap` for `1 < d`.
* `topProj_apply_of_card_le_one`: for `d ≤ 1` the top projector is the identity.
* `norm_topProj_sq_continuousAt_of_topSimple`: `norm_topProj_sq_continuousAt` with the gap
  hypothesis replaced by `TopSimple A`, which is the hypothesis of
  `lem_entrywise_conv_eigenvec`.
* `norm_symMat_sub_lt`, `lamMax_symMat`, `continuous_lamMax_symMat`,
  `continuous_norm_symMat_sub`: the entrywise form of the Lipschitz bound `lamMax_lipschitz`
  and of the operator norm, for `TendstoInProbPi.comp_continuous`.
-/

section Assembly

/-- The operator of a matrix product is the composition of the operators. -/
private theorem toEuclideanLin_mul_apply {p q r : ℕ} (B : Matrix (Fin p) (Fin q) ℝ)
    (A : Matrix (Fin q) (Fin r) ℝ) (x : EuclideanSpace ℝ (Fin r)) :
    Matrix.toEuclideanLin (B * A) x = Matrix.toEuclideanLin B (Matrix.toEuclideanLin A x) := by
  apply WithLp.ofLp_injective
  change (B * A) *ᵥ WithLp.ofLp x = B *ᵥ (A *ᵥ WithLp.ofLp x)
  exact (Matrix.mulVec_mulVec _ _ _).symm

/-- `x ↦ B x` is injective from the `t`-eigenspace of `C B` to that of `B C` when `t ≠ 0`. -/
private theorem finrank_eigenspace_le_of_mul {p q : ℕ} (B : Matrix (Fin p) (Fin q) ℝ)
    (C : Matrix (Fin q) (Fin p) ℝ) {t : ℝ} (ht : t ≠ 0) :
    Module.finrank ℝ (Module.End.eigenspace (Matrix.toEuclideanLin (C * B)) t)
      ≤ Module.finrank ℝ (Module.End.eigenspace (Matrix.toEuclideanLin (B * C)) t) := by
  have hmaps : ∀ x ∈ Module.End.eigenspace (Matrix.toEuclideanLin (C * B)) t,
      Matrix.toEuclideanLin B x ∈ Module.End.eigenspace (Matrix.toEuclideanLin (B * C)) t := by
    intro x hx
    rw [Module.End.mem_eigenspace_iff] at hx ⊢
    rw [toEuclideanLin_mul_apply, ← toEuclideanLin_mul_apply C B x, hx, map_smul]
  have hinj : Function.Injective ((Matrix.toEuclideanLin B).restrict hmaps) := by
    intro x y hxy
    have h : Matrix.toEuclideanLin B (x : EuclideanSpace ℝ (Fin q))
        = Matrix.toEuclideanLin B (y : EuclideanSpace ℝ (Fin q)) := by
      have := congrArg Subtype.val hxy
      simpa [LinearMap.restrict_apply] using this
    have hx := Module.End.mem_eigenspace_iff.mp x.2
    have hy := Module.End.mem_eigenspace_iff.mp y.2
    have hst : t • (x : EuclideanSpace ℝ (Fin q)) = t • (y : EuclideanSpace ℝ (Fin q)) := by
      rw [← hx, ← hy, toEuclideanLin_mul_apply, toEuclideanLin_mul_apply, h]
    exact Subtype.ext (smul_right_injective _ ht hst)
  exact LinearMap.finrank_le_finrank_of_injective hinj

/-- For `t ≠ 0` the eigenspaces of `Aᵀ A` and `A Aᵀ` at `t` have the same dimension. -/
theorem finrank_eigenspace_transpose_mul_eq {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ) {t : ℝ}
    (ht : t ≠ 0) :
    Module.finrank ℝ (Module.End.eigenspace (toOp (Aᵀ * A)) t)
      = Module.finrank ℝ (Module.End.eigenspace (toOp (A * Aᵀ)) t) :=
  le_antisymm (finrank_eigenspace_le_of_mul A Aᵀ ht) (finrank_eigenspace_le_of_mul Aᵀ A ht)

/-- Simplicity of the top eigenvalue transfers between `Aᵀ A` and `A Aᵀ` when the two top
eigenvalues agree and are positive. -/
theorem topSimple_transpose_mul_iff {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (hA : (Aᵀ * A).IsHermitian) (hB : (A * Aᵀ).IsHermitian)
    (hlam : lamMax (Aᵀ * A) hA = lamMax (A * Aᵀ) hB) (hpos : 0 < lamMax (A * Aᵀ) hB) :
    TopSimple (Aᵀ * A) hA ↔ TopSimple (A * Aᵀ) hB := by
  unfold TopSimple
  rw [topSpace_eq_eigenspace, topSpace_eq_eigenspace, hlam,
    finrank_eigenspace_transpose_mul_eq A hpos.ne']

/-- Converse of `topSimple_of_gap`: a simple top eigenvalue has a positive gap below it. -/
theorem gap_pos_of_topSimple {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (h1 : 1 < Fintype.card (Fin d)) (hsimple : TopSimple A hA) :
    0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩ := by
  have hc : 0 < Fintype.card (Fin d) := pos_card h1
  by_contra hle
  rw [not_lt] at hle
  have hle' : hA.eigenvalues₀ ⟨1, h1⟩ ≤ lamMax A hA := by
    rw [lamMax_eq_eigenvalues₀ hA hc]
    exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])
  have heq : hA.eigenvalues₀ ⟨1, h1⟩ = lamMax A hA := le_antisymm hle' (by linarith)
  set b := (symmOp hA).eigenvectorBasis finrank_euclideanSpace with hb
  set f : Fin 2 → Fin (Fintype.card (Fin d)) := fun i => ⟨i.val, by omega⟩ with hf
  have hfinj : Function.Injective f := fun i j hij => Fin.ext (by simpa [hf] using hij)
  have hmem : ∀ i : Fin 2, b (f i) ∈ topSpace A hA := by
    intro i
    rw [topSpace_eq_eigenspace]
    refine Module.End.mem_eigenspace_iff.mpr ?_
    rw [apply_eigvec hA (f i)]
    congr 1
    fin_cases i
    · exact (lamMax_eq_eigenvalues₀ hA hc).symm
    · exact heq
  have hli : LinearIndependent ℝ (fun i : Fin 2 => (⟨b (f i), hmem i⟩ : topSpace A hA)) := by
    refine LinearIndependent.of_comp (topSpace A hA).subtype ?_
    exact (b.orthonormal.comp f hfinj).linearIndependent
  have hcard := hli.fintype_card_le_finrank
  rw [hsimple, Fintype.card_fin] at hcard
  omega

/-- For `d ≤ 1` the top projector is the identity: the top eigenspace is the whole space. -/
theorem topProj_apply_of_card_le_one (hd : Fintype.card (Fin d) ≤ 1)
    {B : Matrix (Fin d) (Fin d) ℝ} (hB : B.IsHermitian) (x : EuclideanSpace ℝ (Fin d)) :
    topProj B hB x = x := by
  rcases Nat.eq_zero_or_pos (Fintype.card (Fin d)) with hc | hc
  · have hd0 : d = 0 := by simpa using hc
    subst hd0
    apply PiLp.ext
    intro i
    exact i.elim0
  · have hE : Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) = 1 := by
      rw [finrank_euclideanSpace]
      omega
    have hne : topEigvec hB hc ≠ 0 := by
      intro h
      have hn := norm_topEigvec hB hc
      rw [h, norm_zero] at hn
      exact zero_ne_one hn
    have hspan : (ℝ ∙ topEigvec hB hc) ≤ topSpace B hB :=
      (Submodule.span_singleton_le_iff_mem _ _).2 (mem_topSpace_topEigvec hB hc)
    have h1 : 1 ≤ Module.finrank ℝ (topSpace B hB) := by
      rw [← finrank_span_singleton (K := ℝ) hne]
      exact Submodule.finrank_mono hspan
    have h2 : Module.finrank ℝ (topSpace B hB) ≤ 1 := hE ▸ Submodule.finrank_le _
    have htop : topSpace B hB = ⊤ := Submodule.eq_top_of_finrank_eq (by omega)
    change (topSpace B hB).starProjection x = x
    rw [Submodule.starProjection_eq_self_iff, htop]
    exact Submodule.mem_top

/-- `norm_topProj_sq_continuousAt` under the hypothesis `TopSimple A` instead of the gap. For
`1 < d` the gap follows from simplicity (`gap_pos_of_topSimple`); for `d ≤ 1` the function is
the constant `‖x‖²`. -/
theorem norm_topProj_sq_continuousAt_of_topSimple {A : Matrix (Fin d) (Fin d) ℝ}
    (hA : A.IsHermitian) (hsimple : TopSimple A hA) (x : EuclideanSpace ℝ (Fin d)) :
    ContinuousAt (fun w : Fin d × Fin d → ℝ =>
      ‖topProj (symMat w) (isHermitian_symMat w) x‖ ^ 2) (fun p => A p.1 p.2) := by
  by_cases h1 : 1 < Fintype.card (Fin d)
  · exact norm_topProj_sq_continuousAt hA h1 (gap_pos_of_topSimple hA h1 hsimple) x
  · rw [not_lt] at h1
    have hconst : (fun w : Fin d × Fin d → ℝ =>
        ‖topProj (symMat w) (isHermitian_symMat w) x‖ ^ 2) = fun _ => ‖x‖ ^ 2 := by
      funext w
      rw [topProj_apply_of_card_le_one h1]
    rw [hconst]
    exact continuousAt_const

/-- Entrywise closeness in the sup metric gives operator-norm closeness of the symmetrized
matrices, with the constant `d + 1`. Same computation as inside
`norm_topProj_sq_continuousAt`. -/
theorem norm_symMat_sub_lt {w w' : Fin d × Fin d → ℝ} {δ : ℝ} (hδ : 0 < δ)
    (hw : dist w w' < δ / ((d : ℝ) + 1)) : ‖symMat w - symMat w'‖ < δ := by
  have hdpos : (0 : ℝ) < (d : ℝ) + 1 := by positivity
  have hentry : ∀ i j, |(symMat w - symMat w') i j| ≤ δ / ((d : ℝ) + 1) := by
    intro i j
    have hij := (dist_pi_lt_iff (show (0 : ℝ) < δ / ((d : ℝ) + 1) by positivity)).mp hw (i, j)
    have hji := (dist_pi_lt_iff (show (0 : ℝ) < δ / ((d : ℝ) + 1) by positivity)).mp hw (j, i)
    simp only [Real.dist_eq] at hij hji
    have hval : (symMat w - symMat w') i j
        = ((w (i, j) - w' (i, j)) + (w (j, i) - w' (j, i))) / 2 := by
      simp only [Matrix.sub_apply, symMat, Matrix.of_apply]
      ring
    rw [hval, abs_div, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 2)]
    have habs : |w (i, j) - w' (i, j) + (w (j, i) - w' (j, i))|
        ≤ |w (i, j) - w' (i, j)| + |w (j, i) - w' (j, i)| := abs_add_le _ _
    linarith
  have h := l2_opNorm_le_of_entries (symMat w - symMat w') hentry
  have hlt : (d : ℝ) * (δ / ((d : ℝ) + 1)) < δ := by
    have heq : (d : ℝ) * (δ / ((d : ℝ) + 1)) = δ - δ / ((d : ℝ) + 1) := by
      field_simp
      ring
    rw [heq]
    linarith [div_pos hδ hdpos]
  exact lt_of_le_of_lt h hlt

/-- Value of the symmetrized `lamMax` at a symmetric matrix. -/
theorem lamMax_symMat {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) :
    lamMax (symMat fun p => A p.1 p.2) (isHermitian_symMat _) = lamMax A hA :=
  lamMax_congr (symMat_entries hA) _ hA

/-- The symmetrized `lamMax` is continuous on `Fin d × Fin d → ℝ` (Weyl, entrywise form). -/
theorem continuous_lamMax_symMat :
    Continuous (fun w : Fin d × Fin d → ℝ => lamMax (symMat w) (isHermitian_symMat w)) := by
  refine Metric.continuous_iff.mpr fun w' ε hε =>
    ⟨ε / ((d : ℝ) + 1), by positivity, fun w hw => ?_⟩
  rw [Real.dist_eq]
  exact lt_of_le_of_lt (lamMax_lipschitz _ _ _ _) (norm_symMat_sub_lt hε hw)

/-- `w ↦ ‖symMat w - A‖` is continuous on `Fin d × Fin d → ℝ`. -/
theorem continuous_norm_symMat_sub (A : Matrix (Fin d) (Fin d) ℝ) :
    Continuous (fun w : Fin d × Fin d → ℝ => ‖symMat w - A‖) := by
  refine Metric.continuous_iff.mpr fun w' ε hε =>
    ⟨ε / ((d : ℝ) + 1), by positivity, fun w hw => ?_⟩
  rw [Real.dist_eq]
  refine lt_of_le_of_lt ((abs_norm_sub_norm_le _ _).trans (le_of_eq ?_))
    (norm_symMat_sub_lt hε hw)
  rw [sub_sub_sub_cancel_right]

end Assembly

section EntrywiseOpNorm

/-! ### Entrywise convergence in probability gives operator-norm convergence

`G N ω` is a symmetric random `M × M` matrix whose entries converge in probability to those of
a symmetric `A`. This section states the one consequence that the perturbation lemmas of this
file and of `LinAlg/SpecProjPerturb.lean` consume. -/

open MeasureTheory

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {M : ℕ}
variable {A : Matrix (Fin M) (Fin M) ℝ} {G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ}

/-- Entrywise convergence gives convergence of the `L2` operator norm of `G_N - A`. -/
theorem opNorm_sub_tendstoInProb (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j)) :
    TendstoInProb μ (fun N ω => ‖G N ω - A‖) 0 := by
  have hentries : TendstoInProbPi μ (fun N ω => fun p : Fin M × Fin M => G N ω p.1 p.2)
      (fun p => A p.1 p.2) := fun p => hconv p.1 p.2
  have h := TendstoInProbPi.comp_continuous (continuous_norm_symMat_sub A).continuousAt hentries
  rw [symMat_entries hA, sub_self, norm_zero] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  rw [symMat_entries (hsymm N ω)]

end EntrywiseOpNorm

end StackedSVD
