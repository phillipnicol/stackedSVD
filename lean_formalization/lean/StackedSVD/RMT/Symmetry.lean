/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT
import StackedSVD.RMT.Simplicity
import StackedSVD.Prob.GaussianMatrix

/-!
# `delocUniform` from right-rotation invariance (decision D6, symmetry route)

This file proves the `delocUniform` field of `SpikedModel.SingleTableLaw` for Gaussian noise,
in both regimes and with no hypothesis on `θ` and no random matrix input. It is the route of
`notes/archive/L2_CHOICES_ANSWERS.md`, top decision 1 (`notes/FLAGGED.md`, D6), and it is the
paper's own remark (`main_paper.tex:350`, Loffler et al. Lemma 4.4).

## The argument

Write `X = θ u vᵀ + t Z` with `Z` a canonical Gaussian matrix and `t = d^{-1/2}`.

1. **Right rotation is a covariance of the top projector** (deterministic). For `Oᵀ O = 1`,
   `(X O)ᵀ (X O) = Oᵀ (Xᵀ X) O`, the two matrices have the same eigenvalues, and the top
   eigenspace of the conjugate is the image of the top eigenspace under `Oᵀ`. Hence
   `overlap (X O) w = overlap X (O w)` (`overlap_mul_right`).
2. **Right rotation preserves the law.** `Z ↦ Z O` preserves `gaussianMatrix`
   (`measurePreserving_mul_right`, two transposes and `gaussianMatrix_map_mul`), and
   `(θ u vᵀ) O = θ u (Oᵀ v)ᵀ = θ u vᵀ` when `Oᵀ v = v`.
3. **Exchangeability.** For two unit vectors `w, w'` orthogonal to `v` the Householder
   reflection `H = 1 - 2 a aᵀ / (a ⬝ a)` with `a = w - w'` fixes `v` and sends `w` to `w'`.
   Steps 1 and 2 then give `E[overlap X w] = E[overlap X w']` (`lintegral_overlap_eq_of_orth`).
4. **Bessel.** On the almost sure event that the top eigenvalue is simple
   (`topSimple_ae_affine`, item S), `overlap X e = ⟪v̂, e⟫²` for a unit top eigenvector `v̂`, so
   the sum over any orthonormal family is at most `‖v̂‖² = 1` (`sum_overlap_le_one`).
5. Taking the family to be an orthonormal basis of `v^⊥`, which has `d - 1` elements, steps 3
   and 4 give `E[overlap X w] ≤ 1/(d-1)` (`lintegral_overlap_le`).
6. **Markov** turns this into `μ {overlap X w ≥ ε} ≤ 1/(ε (d-1))`, uniformly in
   `w ∈ orthUnit N`, and `d N → ∞` finishes (`delocUniform_of_gaussian`).

Nothing here needs `Regime` beyond `d N → ∞`, and nothing needs a threshold on `θ`.

## Statement choices

See `notes/archive/agent_reports/proof_symmetry.md`. The two that matter: the integral bound is
stated as a `lintegral` of `ENNReal.ofReal (overlap ...)`, so that no integrability side condition
appears; and the signal invariance enters the canonical lemmas as the hypothesis `∀ O, Oᵀ *ᵥ v = v →
A * O = A`, which `vecMulVec_mul_of_transpose_mulVec_eq` discharges for the rank one signal.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### 1. Orthogonal matrices and the rotation isometry -/

section Orth

variable {d : ℕ} {O : Matrix (Fin d) (Fin d) ℝ}

/-- A one sided inverse of a square matrix is a two sided inverse. -/
theorem mul_transpose_of_orth (hO : Oᵀ * O = 1) : O * Oᵀ = 1 := mul_eq_one_comm.mp hO

/-- The transpose of an orthogonal matrix is orthogonal. -/
theorem transpose_orth (hO : Oᵀ * O = 1) : (Oᵀ)ᵀ * Oᵀ = 1 := by
  rw [Matrix.transpose_transpose]
  exact mul_transpose_of_orth hO

/-- `Oᵀ` undoes `O` on `EuclideanSpace`. -/
theorem rot_left_inv (hO : Oᵀ * O = 1) (x : EuclideanSpace ℝ (Fin d)) :
    rotIso Oᵀ (transpose_orth hO) (rotIso O hO x) = x := by
  simp only [rotIso_apply, Matrix.mulVec_mulVec, hO, Matrix.one_mulVec, WithLp.toLp_ofLp]

/-- `O` undoes `Oᵀ` on `EuclideanSpace`. -/
theorem rot_right_inv (hO : Oᵀ * O = 1) (x : EuclideanSpace ℝ (Fin d)) :
    rotIso O hO (rotIso Oᵀ (transpose_orth hO) x) = x := by
  simp only [rotIso_apply, Matrix.mulVec_mulVec, mul_transpose_of_orth hO, Matrix.one_mulVec,
    WithLp.toLp_ofLp]

/-- The inverse of the isometry of `Oᵀ` is the isometry of `O`. -/
theorem rotIso_transpose_symm (hO : Oᵀ * O = 1) (x : EuclideanSpace ℝ (Fin d)) :
    (rotIso Oᵀ (transpose_orth hO)).symm x = rotIso O hO x := by
  have h := congrArg (rotIso Oᵀ (transpose_orth hO)).symm (rot_left_inv hO x)
  rw [LinearIsometryEquiv.symm_apply_apply] at h
  exact h.symm

end Orth

/-! ### 2. Conjugation of the spectral data -/

section Conj

variable {d : ℕ} {O : Matrix (Fin d) (Fin d) ℝ}

/-- The operator of `Oᵀ A O` is the operator of `A`, conjugated by the two isometries. -/
theorem toOp_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ)
    (x : EuclideanSpace ℝ (Fin d)) :
    toOp (Oᵀ * A * O) x = rotIso Oᵀ (transpose_orth hO) (toOp A (rotIso O hO x)) := by
  simp only [rotIso_apply]
  change WithLp.toLp 2 ((Oᵀ * A * O) *ᵥ WithLp.ofLp x)
    = WithLp.toLp 2 (Oᵀ *ᵥ (A *ᵥ (O *ᵥ WithLp.ofLp x)))
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]

/-- Eigenspaces of a conjugate are the images of the eigenspaces. -/
theorem eigenspace_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ) (t : ℝ) :
    Module.End.eigenspace (toOp (Oᵀ * A * O)) t =
      (Module.End.eigenspace (toOp A) t).map
        ((rotIso Oᵀ (transpose_orth hO)).toLinearEquiv :
          EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d)) := by
  ext x
  rw [Module.End.mem_eigenspace_iff, Submodule.mem_map]
  constructor
  · intro h
    refine ⟨rotIso O hO x, ?_, rot_left_inv hO x⟩
    rw [Module.End.mem_eigenspace_iff]
    have h2 := congrArg (rotIso O hO) h
    rw [toOp_conj hO A x, rot_right_inv hO, map_smul] at h2
    exact h2
  · rintro ⟨y, hy, rfl⟩
    rw [Module.End.mem_eigenspace_iff] at hy
    change toOp (Oᵀ * A * O) (rotIso Oᵀ (transpose_orth hO) y)
      = t • rotIso Oᵀ (transpose_orth hO) y
    rw [toOp_conj hO A _, rot_right_inv hO, hy, map_smul]

/-- Spectral subspaces of a conjugate are the images of the spectral subspaces. -/
theorem specSpace_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ) (S : Set ℝ) :
    specSpace (Oᵀ * A * O) S =
      (specSpace A S).map ((rotIso Oᵀ (transpose_orth hO)).toLinearEquiv :
        EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d)) := by
  unfold specSpace
  simp only [Submodule.map_iSup, eigenspace_conj hO]

/-! #### `lamMax` is invariant

`lamMax` is the greatest point of the set of eigenvalues, and conjugation by an invertible
matrix is a bijection on eigenvectors, so the two sets are equal. -/

/-- The set of eigenvalues of a real matrix, as a set of reals. -/
def eigSet (A : Matrix (Fin d) (Fin d) ℝ) : Set ℝ :=
  {t | ∃ x : EuclideanSpace ℝ (Fin d), x ≠ 0 ∧ toOp A x = t • x}

/-- Every eigenvalue is at most `lamMax`. -/
theorem le_lamMax_of_mem_eigSet {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d)
    {t : ℝ} (ht : t ∈ eigSet A) : t ≤ lamMax A hA := by
  obtain ⟨x, hx0, hxe⟩ := ht
  set b := (symmOp hA).eigenvectorBasis (finrank_euclideanSpace (ι := Fin d)) with hb
  -- some coordinate of `x` in the eigenbasis is nonzero
  have hex : ∃ i, ⟪b i, x⟫_ℝ ≠ 0 := by
    by_contra hcon
    simp only [not_exists, ne_eq, not_not] at hcon
    have hzero : ‖x‖ ^ 2 = 0 := by
      rw [← b.sum_sq_norm_inner_right x]
      exact Finset.sum_eq_zero fun i _ => by simp [hcon i]
    exact hx0 (by simpa using pow_eq_zero_iff (n := 2) (by norm_num) |>.mp hzero)
  obtain ⟨i, hi⟩ := hex
  have h1 : ⟪b i, toOp A x⟫_ℝ = hA.eigenvalues₀ i * ⟪b i, x⟫_ℝ := by
    rw [← (symmOp hA) (b i) x, apply_eigvec hA i, real_inner_smul_left]
  have h2 : ⟪b i, toOp A x⟫_ℝ = t * ⟪b i, x⟫_ℝ := by
    rw [hxe, real_inner_smul_right]
  have h3 : hA.eigenvalues₀ i = t := mul_right_cancel₀ hi (h1.symm.trans h2)
  rw [lamMax, dif_pos hd, ← h3]
  exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

/-- `lamMax` is an eigenvalue. -/
theorem lamMax_mem_eigSet {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d) :
    lamMax A hA ∈ eigSet A := by
  have hcard : (0 : ℕ) < Fintype.card (Fin d) := by simpa using hd
  set i : Fin (Fintype.card (Fin d)) := ⟨0, hcard⟩ with hi
  refine ⟨(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, ?_, ?_⟩
  · have : ‖(symmOp hA).eigenvectorBasis (finrank_euclideanSpace (ι := Fin d)) i‖ = 1 :=
      ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).orthonormal.1 i
    intro h
    rw [h, norm_zero] at this
    exact zero_ne_one this
  · rw [apply_eigvec hA i, lamMax, dif_pos hd]

/-- `lamMax` is the greatest eigenvalue. -/
theorem lamMax_isGreatest_eigSet {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d) :
    IsGreatest (eigSet A) (lamMax A hA) :=
  ⟨lamMax_mem_eigSet hA hd, fun _ ht => le_lamMax_of_mem_eigSet hA hd ht⟩

/-- Conjugation by an orthogonal matrix does not change the set of eigenvalues. -/
theorem eigSet_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ) :
    eigSet (Oᵀ * A * O) = eigSet A := by
  ext t
  constructor
  · rintro ⟨x, hx0, hxe⟩
    refine ⟨rotIso O hO x, ?_, ?_⟩
    · intro h
      apply hx0
      have := congrArg (rotIso Oᵀ (transpose_orth hO)) h
      rwa [rot_left_inv hO, map_zero] at this
    · have h2 := congrArg (rotIso O hO) hxe
      rw [toOp_conj hO A x, rot_right_inv hO, map_smul] at h2
      exact h2
  · rintro ⟨y, hy0, hye⟩
    refine ⟨rotIso Oᵀ (transpose_orth hO) y, ?_, ?_⟩
    · intro h
      apply hy0
      have := congrArg (rotIso O hO) h
      rwa [rot_right_inv hO, map_zero] at this
    · rw [toOp_conj hO A _, rot_right_inv hO, hye, map_smul]

/-- **`lamMax` is invariant under conjugation by an orthogonal matrix.** -/
theorem lamMax_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) : lamMax (Oᵀ * A * O) hB = lamMax A hA := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    rw [lamMax, lamMax, dif_neg (by omega), dif_neg (by omega)]
  · exact IsGreatest.unique ((eigSet_conj hO A) ▸ lamMax_isGreatest_eigSet hB hd)
      (lamMax_isGreatest_eigSet hA hd)

/-- The top eigenspace of a conjugate is the image of the top eigenspace. -/
theorem topSpace_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) :
    topSpace (Oᵀ * A * O) hB =
      (topSpace A hA).map ((rotIso Oᵀ (transpose_orth hO)).toLinearEquiv :
        EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d)) := by
  unfold topSpace
  rw [lamMax_conj hO hA hB]
  exact specSpace_conj hO A _

/-- **The top projector of a conjugate is the conjugated top projector.** -/
theorem topProj_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) (x : EuclideanSpace ℝ (Fin d)) :
    topProj (Oᵀ * A * O) hB x =
      rotIso Oᵀ (transpose_orth hO) (topProj A hA (rotIso O hO x)) := by
  change (topSpace (Oᵀ * A * O) hB).starProjection x = _
  rw [topSpace_conj hO hA hB, Submodule.starProjection_map_apply, rotIso_transpose_symm hO]
  rfl

end Conj

/-! ### 3. `overlap` under a right rotation -/

/-- Right multiplication by an orthogonal matrix moves the test vector by that matrix. -/
theorem overlap_mul_right {n d : ℕ} (X : Matrix (Fin n) (Fin d) ℝ)
    {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1) (w : EuclideanSpace ℝ (Fin d)) :
    overlap (X * O) w = overlap X (WithLp.toLp 2 (O *ᵥ WithLp.ofLp w)) := by
  have hgram : (X * O)ᵀ * (X * O) = Oᵀ * (Xᵀ * X) * O := by
    simp [Matrix.transpose_mul, Matrix.mul_assoc]
  have hA : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  have hB : (Oᵀ * (Xᵀ * X) * O).IsHermitian := by
    rw [← hgram]
    exact isHermitian_transpose_mul_self _
  change ‖topProj ((X * O)ᵀ * (X * O)) (isHermitian_transpose_mul_self (X * O)) w‖ ^ 2 = _
  rw [topProj_congr hgram _ hB, topProj_conj hO hA hB w, LinearIsometryEquiv.norm_map]
  rfl

/-! ### 4. Householder reflections -/

section Householder

variable {d : ℕ}

/-- `vecMulVec a b` applied to `x` is `(b ⬝ x) • a`. -/
theorem vecMulVec_mulVec' (a b x : Fin d → ℝ) :
    Matrix.vecMulVec a b *ᵥ x = (b ⬝ᵥ x) • a := by
  ext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- `vecMulVec a b * vecMulVec c e = (b ⬝ c) • vecMulVec a e`. -/
theorem vecMulVec_mul_vecMulVec (a b c e : Fin d → ℝ) :
    Matrix.vecMulVec a b * Matrix.vecMulVec c e = (b ⬝ᵥ c) • Matrix.vecMulVec a e := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.vecMulVec_apply, Matrix.smul_apply, smul_eq_mul, dotProduct]
  rw [Finset.sum_mul]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The Householder reflection in the hyperplane orthogonal to `a`. -/
noncomputable def householder (a : Fin d → ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  1 - (2 / (a ⬝ᵥ a)) • Matrix.vecMulVec a a

theorem householder_transpose (a : Fin d → ℝ) : (householder a)ᵀ = householder a := by
  simp [householder, Matrix.transpose_sub, Matrix.transpose_smul, Matrix.transpose_vecMulVec]

theorem householder_mulVec (a x : Fin d → ℝ) :
    householder a *ᵥ x = x - ((2 / (a ⬝ᵥ a)) * (a ⬝ᵥ x)) • a := by
  rw [householder, Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.smul_mulVec,
    vecMulVec_mulVec', smul_smul]

/-- A Householder reflection fixes every vector orthogonal to its axis. -/
theorem householder_apply_of_orth {a x : Fin d → ℝ} (h : a ⬝ᵥ x = 0) :
    householder a *ᵥ x = x := by
  rw [householder_mulVec, h, mul_zero, zero_smul, sub_zero]

/-- A Householder reflection is an involution. -/
theorem householder_mul_self {a : Fin d → ℝ} (ha : a ⬝ᵥ a ≠ 0) :
    householder a * householder a = 1 := by
  have hS : Matrix.vecMulVec a a * Matrix.vecMulVec a a
      = (a ⬝ᵥ a) • Matrix.vecMulVec a a := vecMulVec_mul_vecMulVec a a a a
  have hcc : 2 / (a ⬝ᵥ a) * (2 / (a ⬝ᵥ a)) * (a ⬝ᵥ a) = 2 * (2 / (a ⬝ᵥ a)) := by
    field_simp
  rw [householder]
  simp only [Matrix.sub_mul, Matrix.mul_sub, Matrix.one_mul, Matrix.mul_one, Matrix.smul_mul,
    Matrix.mul_smul, smul_smul, hS]
  match_scalars <;> (field_simp; try ring)

theorem householder_orth {a : Fin d → ℝ} (ha : a ⬝ᵥ a ≠ 0) :
    (householder a)ᵀ * householder a = 1 := by
  rw [householder_transpose]
  exact householder_mul_self ha

/-- The Householder reflection with axis `w - w'` sends `w` to `w'`, for two distinct unit
vectors. -/
theorem householder_sub_apply {w w' : Fin d → ℝ} (hw : w ⬝ᵥ w = 1) (hw' : w' ⬝ᵥ w' = 1)
    (hne : w ≠ w') : householder (w - w') *ᵥ w = w' := by
  have ha0 : w - w' ≠ 0 := sub_ne_zero.mpr hne
  have haa : (w - w') ⬝ᵥ (w - w') ≠ 0 := fun h => ha0 (dotProduct_self_eq_zero.mp h)
  have hexp : (w - w') ⬝ᵥ (w - w') = 2 - 2 * (w ⬝ᵥ w') := by
    rw [sub_dotProduct, dotProduct_sub, dotProduct_sub, hw, hw', dotProduct_comm w' w]
    ring
  have hexpw : (w - w') ⬝ᵥ w = 1 - (w ⬝ᵥ w') := by
    rw [sub_dotProduct, hw, dotProduct_comm w' w]
  have h1 : (1 : ℝ) - w ⬝ᵥ w' ≠ 0 := by
    intro h
    apply haa
    rw [hexp]
    linarith
  have hscale : (2 / ((w - w') ⬝ᵥ (w - w'))) * ((w - w') ⬝ᵥ w) = 1 := by
    rw [hexp, hexpw, show (2 : ℝ) - 2 * (w ⬝ᵥ w') = 2 * (1 - w ⬝ᵥ w') by ring]
    field_simp
  rw [householder_mulVec, hscale, one_smul]
  abel

end Householder

/-! ### 5. Right rotation preserves the Gaussian matrix law -/

/-- **Right multiplication by an orthogonal matrix preserves `gaussianMatrix`.** -/
theorem measurePreserving_mul_right {d : ℕ} {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1)
    (n : ℕ) :
    MeasurePreserving (fun Z : Matrix (Fin n) (Fin d) ℝ => Z * O)
      (gaussianMatrix n d) (gaussianMatrix n d) := by
  have hcomp := (measurePreserving_transpose d n).comp
    ((measurePreserving_mul_left (U := Oᵀ) (transpose_orth hO) n).comp
      (measurePreserving_transpose n d))
  have hfun : (fun Z : Matrix (Fin n) (Fin d) ℝ => Z * O) =
      (Matrix.transpose ∘ (fun Y : Matrix (Fin d) (Fin n) ℝ => Oᵀ * Y) ∘
        (Matrix.transpose : Matrix (Fin n) (Fin d) ℝ → Matrix (Fin d) (Fin n) ℝ)) := by
    funext Z
    change Z * O = (Oᵀ * Zᵀ)ᵀ
    rw [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.transpose_transpose]
  rw [hfun]
  exact hcomp

/-- A rank one signal is fixed by a right rotation that fixes its right factor. -/
theorem vecMulVec_mul_of_transpose_mulVec_eq {n d : ℕ} (u : Fin n → ℝ) (v : Fin d → ℝ)
    {O : Matrix (Fin d) (Fin d) ℝ} (hv : Oᵀ *ᵥ v = v) :
    Matrix.vecMulVec u v * O = Matrix.vecMulVec u v := by
  have hgen : ∀ P : Matrix (Fin d) (Fin d) ℝ,
      Matrix.vecMulVec u v * P = Matrix.vecMulVec u (Pᵀ *ᵥ v) := by
    intro P
    ext i k
    simp only [Matrix.mul_apply, Matrix.vecMulVec_apply, Matrix.mulVec, dotProduct,
      Matrix.transpose_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hgen O, hv]

/-! ### 6. Measurability of the affine family -/

private theorem measurable_affine' {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A + t • Z := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  change Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A i j + t * Z i j
  have hij : Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => Z i j :=
    (measurable_pi_apply j).comp (measurable_pi_apply i)
  exact (hij.const_mul t).const_add (A i j)

private theorem measurable_overlap_affine {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ)
    (w : EuclideanSpace ℝ (Fin d)) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => ENNReal.ofReal (overlap (A + t • Z) w) :=
  ENNReal.measurable_ofReal.comp ((measurable_overlap w).comp (measurable_affine' A t))

/-! ### 7. Exchangeability of the test direction -/

/-- **The law of `overlap X w` is the same for every unit `w` orthogonal to `v`.**
Stated as equality of the two `lintegral`s, which is all that the bound needs. -/
theorem lintegral_overlap_eq_of_orth {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ)
    (v : EuclideanSpace ℝ (Fin d))
    (hA : ∀ O : Matrix (Fin d) (Fin d) ℝ,
      Oᵀ *ᵥ WithLp.ofLp v = WithLp.ofLp v → A * O = A)
    {w w' : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) (hw' : ‖w'‖ = 1)
    (hwv : ⟪w, v⟫_ℝ = 0) (hw'v : ⟪w', v⟫_ℝ = 0) :
    ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w) ∂(gaussianMatrix n d)
      = ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w') ∂(gaussianMatrix n d) := by
  by_cases hww : w = w'
  · rw [hww]
  -- dot products of the underlying plain vectors
  have hdot : ∀ x y : EuclideanSpace ℝ (Fin d),
      WithLp.ofLp x ⬝ᵥ WithLp.ofLp y = ⟪x, y⟫_ℝ := fun x y =>
    (inner_euclidean_eq_dotProduct x y).symm
  have hww1 : WithLp.ofLp w ⬝ᵥ WithLp.ofLp w = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw]
    norm_num
  have hww'1 : WithLp.ofLp w' ⬝ᵥ WithLp.ofLp w' = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw']
    norm_num
  have hne : WithLp.ofLp w ≠ WithLp.ofLp w' := fun h => hww (by
    have := congrArg (WithLp.toLp 2) h
    simpa using this)
  set a : Fin d → ℝ := WithLp.ofLp w - WithLp.ofLp w' with ha
  set O : Matrix (Fin d) (Fin d) ℝ := householder a with hOdef
  have ha0 : a ⬝ᵥ a ≠ 0 := fun h => (sub_ne_zero.mpr hne) (dotProduct_self_eq_zero.mp h)
  have hO : Oᵀ * O = 1 := householder_orth ha0
  have hav : a ⬝ᵥ WithLp.ofLp v = 0 := by
    rw [ha, sub_dotProduct, hdot, hdot, hwv, hw'v, sub_zero]
  have hOv : Oᵀ *ᵥ WithLp.ofLp v = WithLp.ofLp v := by
    rw [hOdef, householder_transpose]
    exact householder_apply_of_orth hav
  have hOw : O *ᵥ WithLp.ofLp w = WithLp.ofLp w' :=
    householder_sub_apply hww1 hww'1 hne
  -- the change of variables
  have key : ∀ Z : Matrix (Fin n) (Fin d) ℝ,
      overlap (A + t • (Z * O)) w = overlap (A + t • Z) w' := by
    intro Z
    have hmul : (A + t • Z) * O = A + t • (Z * O) := by
      rw [Matrix.add_mul, hA O hOv, Matrix.smul_mul]
    rw [← hmul, overlap_mul_right _ hO, hOw]
  have hmeas : Measurable
      (fun Z : Matrix (Fin n) (Fin d) ℝ => ENNReal.ofReal (overlap (A + t • Z) w)) :=
    measurable_overlap_affine A t w
  calc ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w) ∂(gaussianMatrix n d)
      = ∫⁻ Z, ENNReal.ofReal (overlap (A + t • (Z * O)) w) ∂(gaussianMatrix n d) :=
        ((measurePreserving_mul_right hO n).lintegral_comp hmeas).symm
    _ = ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w') ∂(gaussianMatrix n d) := by
        simp only [key]

/-! ### 8. Bessel: the sum over an orthonormal family is at most one -/

/-- On the event that the top eigenvalue is simple, `overlap X ·` is `⟪v̂, ·⟫²`, so Bessel's
inequality bounds the sum over any orthonormal family by `1`. -/
theorem sum_overlap_le_one {n d : ℕ} {ι : Type*} [Fintype ι] (X : Matrix (Fin n) (Fin d) ℝ)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {e : ι → EuclideanSpace ℝ (Fin d)} (he : Orthonormal ℝ e) :
    ∑ k, overlap X (e k) ≤ 1 := by
  set K := topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) with hK
  have hpos : 0 < Module.finrank ℝ K := by rw [hsimple]; norm_num
  have hnt : Nontrivial K := Module.finrank_pos_iff.mp hpos
  obtain ⟨y, hy⟩ := exists_ne (0 : K)
  have hy0 : (y : EuclideanSpace ℝ (Fin d)) ≠ 0 := fun h => hy (Subtype.ext h)
  set f : EuclideanSpace ℝ (Fin d) := ‖(y : EuclideanSpace ℝ (Fin d))‖⁻¹ • (y : _) with hf
  have hfK : f ∈ K := K.smul_mem _ y.2
  have hfn : ‖f‖ = 1 := by
    rw [hf, norm_smul, norm_inv, norm_norm]
    field_simp
  have hval : ∀ k, overlap X (e k) = ⟪f, e k⟫_ℝ ^ 2 := fun k =>
    overlap_eq_inner_sq X (e k) hsimple hfK hfn
  have hcomm : ∀ k, overlap X (e k) = ‖⟪e k, f⟫_ℝ‖ ^ 2 := by
    intro k
    rw [hval k, real_inner_comm, Real.norm_eq_abs, sq_abs]
  calc ∑ k, overlap X (e k) = ∑ k, ‖⟪e k, f⟫_ℝ‖ ^ 2 :=
        Finset.sum_congr rfl fun k _ => hcomm k
    _ ≤ ‖f‖ ^ 2 := Orthonormal.sum_inner_products_le f he
    _ = 1 := by rw [hfn]; norm_num

/-! ### 9. The uniform bound `E[overlap X w] ≤ 1/(d-1)` -/

/-- **Step 3 of the task.** For a signal `A` fixed by every right rotation fixing `v`, and for
every unit `w` orthogonal to `v`, the mean overlap is at most `1/(d-1)`. -/
theorem lintegral_overlap_le {n d : ℕ} (hn : 0 < n) (hd : 2 ≤ d)
    (A : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0)
    {v : EuclideanSpace ℝ (Fin d)} (hv : ‖v‖ = 1)
    (hA : ∀ O : Matrix (Fin d) (Fin d) ℝ,
      Oᵀ *ᵥ WithLp.ofLp v = WithLp.ofLp v → A * O = A)
    {w : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) (hwv : ⟪w, v⟫_ℝ = 0) :
    ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w) ∂(gaussianMatrix n d)
      ≤ ENNReal.ofReal (1 / ((d : ℝ) - 1)) := by
  have hd0 : 0 < d := by omega
  -- an orthonormal basis of `v^⊥`
  have hv0 : v ≠ 0 := fun h => by rw [h, norm_zero] at hv; exact zero_ne_one hv
  have hrk : Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) = d - 1 := by
    have h1 : Module.finrank ℝ (ℝ ∙ v) = 1 := finrank_span_singleton hv0
    have h2 := Submodule.finrank_add_finrank_orthogonal (𝕜 := ℝ) (K := (ℝ ∙ v))
    rw [h1, finrank_euclideanSpace_fin] at h2
    omega
  set b := stdOrthonormalBasis ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) with hb
  set e : Fin (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d)))) →
      EuclideanSpace ℝ (Fin d) := fun k => (b k : EuclideanSpace ℝ (Fin d)) with he
  have hon : Orthonormal ℝ e :=
    (((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))).subtypeₗᵢ.orthonormal_comp_iff
      (v := fun k => b k)).mpr b.orthonormal
  have henorm : ∀ k, ‖e k‖ = 1 := fun k => hon.1 k
  have hev : ∀ k, ⟪e k, v⟫_ℝ = 0 := by
    intro k
    have : (b k : EuclideanSpace ℝ (Fin d)) ∈ (ℝ ∙ v)ᗮ := (b k).2
    rw [Submodule.mem_orthogonal_singleton_iff_inner_right] at this
    rw [real_inner_comm]
    exact this
  -- every direction has the same mean overlap
  set I : ℝ≥0∞ := ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w) ∂(gaussianMatrix n d) with hI
  have hsame : ∀ k, ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) (e k)) ∂(gaussianMatrix n d) = I :=
    fun k => (lintegral_overlap_eq_of_orth A t v hA hw (henorm k) hwv (hev k)).symm
  -- the sum over the basis is at most one, almost surely
  have hae : ∀ᵐ Z ∂(gaussianMatrix n d),
      ENNReal.ofReal (∑ k, overlap (A + t • Z) (e k)) ≤ 1 := by
    filter_upwards [topSimple_ae_affine hn hd0 A ht] with Z hZ
    exact ENNReal.ofReal_le_one.mpr (sum_overlap_le_one _ hZ hon)
  have hsum : ∑ k, ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) (e k)) ∂(gaussianMatrix n d)
      ≤ 1 := by
    rw [← lintegral_finsetSum _ fun k _ => measurable_overlap_affine A t (e k)]
    calc ∫⁻ Z, ∑ k, ENNReal.ofReal (overlap (A + t • Z) (e k)) ∂(gaussianMatrix n d)
        = ∫⁻ Z, ENNReal.ofReal (∑ k, overlap (A + t • Z) (e k)) ∂(gaussianMatrix n d) := by
          refine lintegral_congr fun Z => ?_
          rw [ENNReal.ofReal_sum_of_nonneg fun k _ => overlap_nonneg _ _]
      _ ≤ ∫⁻ _, (1 : ℝ≥0∞) ∂(gaussianMatrix n d) := lintegral_mono_ae hae
      _ = 1 := by simp
  -- conclude
  simp only [hsame, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
  have hcard : (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) : ℝ≥0∞)
      = ENNReal.ofReal ((d : ℝ) - 1) := by
    rw [hrk]
    rw [show ((d - 1 : ℕ) : ℝ≥0∞) = ENNReal.ofReal ((d - 1 : ℕ) : ℝ) by
      rw [ENNReal.ofReal_natCast]]
    congr 1
    have : (1 : ℕ) ≤ d := by omega
    push_cast [Nat.cast_sub this]
    ring
  rw [hcard] at hsum
  have hpos : (0 : ℝ) < (d : ℝ) - 1 := by
    have : (2 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
    linarith
  have hne0 : ENNReal.ofReal ((d : ℝ) - 1) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    exact hpos
  have hnetop : ENNReal.ofReal ((d : ℝ) - 1) ≠ ⊤ := ENNReal.ofReal_ne_top
  have hdiv : I ≤ 1 / ENNReal.ofReal ((d : ℝ) - 1) := by
    rw [ENNReal.le_div_iff_mul_le (Or.inl hne0) (Or.inl hnetop), mul_comm]
    exact hsum
  refine hdiv.trans (le_of_eq ?_)
  rw [ENNReal.ofReal_div_of_pos hpos, ENNReal.ofReal_one]

/-! ### 10. Transfer to the spiked model, Markov, and `delocUniform` -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- The mean overlap bound for the model, for one `N`. -/
theorem lintegral_overlap_model_le (m : SpikedModel μ n d) (hG : m.GaussianNoise) (N : ℕ)
    (hd : 2 ≤ d N) {w : EuclideanSpace ℝ (Fin (d N))} (hw : w ∈ m.orthUnit N) :
    ∫⁻ ω, ENNReal.ofReal (overlap (m.X N ω) w) ∂(μ N)
      ≤ ENNReal.ofReal (1 / ((d N : ℝ) - 1)) := by
  set A : Matrix (Fin (n N)) (Fin (d N)) ℝ :=
    m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) with hAdef
  set t : ℝ := (Real.sqrt (d N))⁻¹ with htdef
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.hd N
  have ht : t ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  have hAinv : ∀ O : Matrix (Fin (d N)) (Fin (d N)) ℝ,
      Oᵀ *ᵥ WithLp.ofLp (m.v N) = WithLp.ofLp (m.v N) → A * O = A := by
    intro O hO
    rw [hAdef, Matrix.smul_mul, vecMulVec_mul_of_transpose_mulVec_eq _ _ hO]
  have hX : ∀ ω, m.X N ω = A + t • m.Z N ω := fun ω => rfl
  have hfun : (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ENNReal.ofReal (overlap (A + t • Z) w)) ∘ (m.Z N) =
      fun ω => ENNReal.ofReal (overlap (m.X N ω) w) := by
    funext ω
    rw [Function.comp_apply, hX ω]
  calc ∫⁻ ω, ENNReal.ofReal (overlap (m.X N ω) w) ∂(μ N)
      = ∫⁻ Z, ENNReal.ofReal (overlap (A + t • Z) w) ∂(gaussianMatrix (n N) (d N)) := by
        rw [← hfun]
        exact (hG N).lintegral_comp (measurable_overlap_affine A t w).aemeasurable
    _ ≤ ENNReal.ofReal (1 / ((d N : ℝ) - 1)) :=
        lintegral_overlap_le (m.hn N) hd A ht (m.hv N) hAinv hw.1 hw.2

/-- Markov's inequality turns the mean bound into a tail bound. -/
theorem measure_overlap_ge_le (m : SpikedModel μ n d) (hG : m.GaussianNoise) (N : ℕ)
    (hd : 2 ≤ d N) {ε : ℝ} (hε : 0 < ε) {w : EuclideanSpace ℝ (Fin (d N))}
    (hw : w ∈ m.orthUnit N) :
    μ N {ω | ε ≤ overlap (m.X N ω) w} ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1))) := by
  have hmeas : AEMeasurable
      (fun ω => ENNReal.ofReal (overlap (m.X N ω) w)) (μ N) := by
    have h1 : Measurable fun ω => ENNReal.ofReal (overlap (m.X N ω) w) := by
      have : (fun ω => ENNReal.ofReal (overlap (m.X N ω) w)) =
          (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
            ENNReal.ofReal (overlap ((m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N))
              (WithLp.ofLp (m.v N))) + (Real.sqrt (d N))⁻¹ • Z) w)) ∘ (m.Z N) := rfl
      rw [this]
      exact (measurable_overlap_affine _ _ w).comp (m.hZ N)
    exact h1.aemeasurable
  have hset : {ω | ε ≤ overlap (m.X N ω) w}
      = {ω | ENNReal.ofReal ε ≤ ENNReal.ofReal (overlap (m.X N ω) w)} := by
    ext ω
    simp [ENNReal.ofReal_le_ofReal_iff (overlap_nonneg _ _)]
  have hεne : ENNReal.ofReal ε ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hε]
  have hmk := meas_ge_le_lintegral_div hmeas hεne ENNReal.ofReal_ne_top
  rw [hset]
  refine hmk.trans ?_
  have hbound := lintegral_overlap_model_le m hG N hd hw
  refine (ENNReal.div_le_div_right hbound _).trans (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hε]
  congr 1
  field_simp

/-- **The `delocUniform` field of `SingleTableLaw`, for Gaussian noise.** Both regimes, no
hypothesis on `θ`, and uniform over `w ∈ orthUnit N` with no extra work. -/
theorem delocUniform_of_gaussian (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) :
    ∀ ε > 0, Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w})
      atTop (𝓝 0) := by
  intro ε hε
  have hle : ∀ᶠ N in atTop,
      (⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w})
        ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1))) := by
    filter_upwards [hd (Filter.eventually_ge_atTop 2)] with N hN
    exact iSup₂_le fun w hw => measure_overlap_ge_le m hG N hN hε hw
  have hreal : Tendsto (fun N => 1 / (ε * ((d N : ℝ) - 1))) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => ε * ((d N : ℝ) - 1)) atTop atTop := by
      refine Filter.Tendsto.const_mul_atTop hε ?_
      exact (tendsto_natCast_atTop_atTop.comp hd).atTop_add tendsto_const_nhds
    exact Filter.Tendsto.congr (fun N => (one_div _).symm) h1.inv_tendsto_atTop
  have htend : Tendsto (fun N => ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1)))) atTop (𝓝 0) := by
    have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
    rw [ENNReal.ofReal_zero] at h3
    exact h3
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend ?_ hle
  exact Filter.Eventually.of_forall fun _ => by simp


end SpikedModel

end StackedSVD
