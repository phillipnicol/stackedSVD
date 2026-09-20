/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R2het
import StackedSVD.RankR.RMT.DelocR
import StackedSVD.RankR.RMT.EdgeGlueR
import StackedSVD.RankR.Het.Split
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.Het.Duality
import StackedSVD.RankR.GramR

/-!
# Heteroscedastic delocalization of the index projector (Track E, item E6)

Model-free layer (`StackedSVD.HetDeloc`). For `Wsig τ d B = d⁻¹ Ysig Ysigᵀ` with
`Ysig τ B = diagonal τ * B` and `B` a `p × q` Gaussian matrix, the spectral projector at
index `a` has a matrix `projMat` whose entries are bounded by 1, whose law is symmetric under
signed permutations inside each weight block (`blockSym_specProjIdx`), and whose diagonal
trace is at most 1 at a simple index. The block-symmetry lemma
`BlockSym.integral_qform_blocks` of `R2het` turns this into the mean bound

  `E ‖P_a y‖² ≤ ∑ i, ‖y_i‖² / |block i|`

(`lintegral_normSq_specProjIdx_Wsig_le`), the weaker corollary with `‖y‖² ∑ i |block i|⁻¹`,
and the product form with a first component of arbitrary law.

Model layer (`StackedSVD.UnalignedModelR`). `W0hetR` is `Wsig` of the Gaussian block from
`exists_block_hasLaw_hetR`, the columns of `QmatHetR` are deterministic functions of the
other Gaussian block, and Markov gives `measure_kappa_ge_le_het`, the hetero twin of
`EdgeGlueR.measure_kappa_ge_le`.

Simplicity route. The a.s. simplicity of the `p × p` Gram `Ysig Ysigᵀ` at indices below
`min p q` uses the polynomial witnesses of `RMT/Het/Simplicity.lean` (`gramRowPolyLeft`
when `p ≤ q`, and a new `gramRowPolyRightDet` when `q ≤ p`, transferred through
`simpleSpec_p_of_gram_left`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace HetDeloc

open HetStein HetR1 HetR2 R2

variable {p q d : ℕ}

/-! ## The matrix of the index projector -/

/-- The matrix of the spectral projector at index `a` of the symmetric matrix `A`. -/
noncomputable def projMat (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ) :
    Matrix (Fin p) (Fin p) ℝ :=
  Matrix.toEuclideanLin.symm
    (specProjIdx A hA a : EuclideanSpace ℝ (Fin p) →ₗ[ℝ] EuclideanSpace ℝ (Fin p))

theorem toEuclideanLin_projMat (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ) :
    Matrix.toEuclideanLin (projMat A hA a) =
      (specProjIdx A hA a : EuclideanSpace ℝ (Fin p) →ₗ[ℝ] EuclideanSpace ℝ (Fin p)) :=
  LinearEquiv.apply_symm_apply _ _

theorem projMat_mulVec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ)
    (y : EuclideanSpace ℝ (Fin p)) :
    projMat A hA a *ᵥ WithLp.ofLp y = WithLp.ofLp (specProjIdx A hA a y) := by
  have h : Matrix.toEuclideanLin (projMat A hA a) y = specProjIdx A hA a y := by
    rw [toEuclideanLin_projMat]; rfl
  rw [← h]; rfl

theorem projMat_apply (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ)
    (j k : Fin p) :
    projMat A hA a j k = (specProjIdx A hA a (EuclideanSpace.single k 1)) j := by
  have h := projMat_mulVec A hA a (EuclideanSpace.single k 1)
  rw [PiLp.ofLp_single, Matrix.mulVec_single_one] at h
  have h' := congrFun h j
  rw [Matrix.col_apply] at h'
  exact h'

theorem projMat_congr {A B : Matrix (Fin p) (Fin p) ℝ} (h : A = B) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (a : ℕ) : projMat A hA a = projMat B hB a := by
  subst h; rfl

theorem normSq_specProj_eq_inner (A : Matrix (Fin p) (Fin p) ℝ) (S : Set ℝ)
    (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProj A S y‖ ^ 2 = ⟪y, specProj A S y⟫_ℝ := by
  rw [← real_inner_self_eq_norm_sq, specProj, Submodule.inner_starProjection_left_eq_right]
  congr 1
  exact Submodule.starProjection_eq_self_iff.mpr (Submodule.starProjection_apply_mem _ _)

theorem normSq_specProjIdx_eq_dot (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ)
    (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx A hA a y‖ ^ 2 = WithLp.ofLp y ⬝ᵥ (projMat A hA a *ᵥ WithLp.ofLp y) := by
  rw [projMat_mulVec, ← inner_euclidean_eq_dotProduct]
  exact normSq_specProj_eq_inner A _ y

theorem projMat_diag_eq (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ)
    (j : Fin p) :
    projMat A hA a j j = ‖specProjIdx A hA a (EuclideanSpace.single j 1)‖ ^ 2 := by
  rw [normSq_specProjIdx_eq_dot, PiLp.ofLp_single, Matrix.mulVec_single_one,
    single_dotProduct, one_mul, Matrix.col_apply]

theorem abs_projMat_le_one (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (a : ℕ)
    (j k : Fin p) : |projMat A hA a j k| ≤ 1 := by
  rw [projMat_apply, ← Real.norm_eq_abs]
  refine (PiLp.norm_apply_le _ _).trans ?_
  refine (norm_specProj_le A _ _).trans ?_
  simp [PiLp.norm_single]

theorem projMat_conj {O : Matrix (Fin p) (Fin p) ℝ} (hO : Oᵀ * O = 1)
    {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (hB : (Oᵀ * A * O).IsHermitian)
    (a : ℕ) : projMat (Oᵀ * A * O) hB a = Oᵀ * projMat A hA a * O := by
  apply Matrix.toEuclideanLin.injective
  rw [toEuclideanLin_projMat]
  refine LinearMap.ext fun y => ?_
  change specProjIdx (Oᵀ * A * O) hB a y =
    WithLp.toLp 2 ((Oᵀ * projMat A hA a * O) *ᵥ WithLp.ofLp y)
  rw [specProjIdx_conj hO hA hB a y, rotIso_apply, rotIso_apply, ← Matrix.mulVec_mulVec,
    ← Matrix.mulVec_mulVec]
  congr 1
  congr 1
  exact (projMat_mulVec _ _ _ _).symm

theorem sum_projMat_diag_le_one {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {a : ℕ}
    (hsimple : SimpleSpec A hA (a + 1)) (ha : a < p) :
    ∑ j, projMat A hA a j j ≤ 1 := by
  have h := sum_normSq_specProjIdx_le_one hsimple ha
    (EuclideanSpace.basisFun (Fin p) ℝ).orthonormal
  simp only [EuclideanSpace.basisFun_apply] at h
  simpa only [projMat_diag_eq] using h

/-! ## `Wsig` as a Gram matrix -/

theorem specProjIdx_Wsig_eq (τ : Fin p → ℝ) (hd : 0 < d) (B : Matrix (Fin p) (Fin q) ℝ)
    (a : ℕ) :
    specProjIdx (Wsig τ d B) (isHermitian_Wsig τ d B) a =
      specProjIdx ((Ysig τ B)ᵀᵀ * (Ysig τ B)ᵀ) (isHermitian_transpose_mul_self (Ysig τ B)ᵀ)
        a := by
  have hpos : (0 : ℝ) < ((d : ℝ))⁻¹ := inv_pos.mpr (by exact_mod_cast hd)
  calc specProjIdx (Wsig τ d B) (isHermitian_Wsig τ d B) a
      = specProjIdx (Ysig τ B * (Ysig τ B)ᵀ) (isHermitian_mul_transpose_self (Ysig τ B)) a :=
        EdgeGlueDetR.specProjIdx_smul (isHermitian_mul_transpose_self (Ysig τ B))
          (isHermitian_Wsig τ d B) hpos a
    _ = _ := specProjIdx_congr (by rw [Matrix.transpose_transpose]) _ _ a

theorem measurable_Ysig_transpose (τ : Fin p → ℝ) :
    Measurable fun B : Matrix (Fin p) (Fin q) ℝ => (Ysig τ B)ᵀ := by
  refine measurable_pi_lambda _ fun k => measurable_pi_lambda _ fun j => ?_
  change Measurable fun B : Matrix (Fin p) (Fin q) ℝ => (Matrix.diagonal τ * B) j k
  simp only [Matrix.diagonal_mul]
  have h1 : Measurable fun B : Matrix (Fin p) (Fin q) ℝ => B j k :=
    (measurable_pi_apply k).comp (measurable_pi_apply j)
  exact h1.const_mul (τ j)

theorem measurable_specProjIdx_Wsig (τ : Fin p → ℝ) (hd : 0 < d) (a : ℕ)
    (y : EuclideanSpace ℝ (Fin p)) :
    Measurable fun B : Matrix (Fin p) (Fin q) ℝ =>
      specProjIdx (Wsig τ d B) (isHermitian_Wsig τ d B) a y := by
  have h : (fun B : Matrix (Fin p) (Fin q) ℝ =>
      specProjIdx (Wsig τ d B) (isHermitian_Wsig τ d B) a y) =
      (fun X : Matrix (Fin q) (Fin p) ℝ =>
        specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) a y) ∘
        (fun B => (Ysig τ B)ᵀ) := by
    funext B
    simp only [Function.comp]
    rw [specProjIdx_Wsig_eq τ hd]
  rw [h]
  exact (measurable_specProjIdx_apply a y).comp (measurable_Ysig_transpose τ)

theorem measurable_projMat_Wsig (τ : Fin p → ℝ) (hd : 0 < d) (a : ℕ) (j k : Fin p) :
    Measurable fun B : Matrix (Fin p) (Fin q) ℝ =>
      projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j k := by
  simp_rw [projMat_apply]
  exact (measurable_pi_apply j).comp
    ((WithLp.measurable_ofLp 2 _).comp (measurable_specProjIdx_Wsig τ hd a _))

theorem normSq_specProjIdx_Wsig_eq_overlapIdx (τ : Fin p → ℝ) (hd : 0 < d)
    (B : Matrix (Fin p) (Fin q) ℝ) (a : ℕ) (y : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx (Wsig τ d B) (isHermitian_Wsig τ d B) a y‖ ^ 2 =
      overlapIdx (Ysig τ B)ᵀ a y := by
  rw [overlapIdx, specProjIdx_Wsig_eq τ hd]

/-! ## Block symmetry -/

theorem blockSym_specProjIdx (τ : Fin p → ℝ) (hd : 0 < d) (a : ℕ) :
    BlockSym τ q (fun B : Matrix (Fin p) (Fin q) ℝ =>
      R4C.cmat (projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a)) 1 where
  conj σ ε hε hτ B := by
    have hW := Wsig_rowSgn τ d σ ε hτ B
    have hH : ((sgnPermMat σ ε)ᵀ * Wsig τ d B * sgnPermMat σ ε).IsHermitian := by
      rw [← hW]; exact isHermitian_Wsig τ d _
    rw [projMat_congr hW _ hH a, projMat_conj (sgnPermMat_orth hε) (isHermitian_Wsig τ d B) hH a,
      R4C.cmat_mul, R4C.cmat_mul, R4C.cmat_transpose]
  meas j k := Complex.measurable_ofReal.comp (measurable_projMat_Wsig τ hd a j k)
  bound B j k := by
    simp only [R4C.cmat, Matrix.map_apply, Complex.norm_real, Real.norm_eq_abs]
    exact abs_projMat_le_one _ _ _ _ _

/-! ## Almost sure simplicity of the `p × p` Gram of `Ysig` -/

theorem simpleSpec_smul {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {t : ℝ}
    (hAt : (t • A).IsHermitian) (ht : 0 < t) {rk : ℕ} (hs : SimpleSpec A hA rk) :
    SimpleSpec (t • A) hAt rk := by
  intro k l hk hkl heq
  rw [EdgeGlueDetR.eigenvalues₀_smul hA hAt ht, EdgeGlueDetR.eigenvalues₀_smul hA hAt ht] at heq
  exact hs k l hk hkl (mul_left_cancel₀ ht.ne' heq)

/-- The separability witness times the determinant of the `q × q` Gram of the affine row
scaling `A + rowScale s Y`. It vanishes at `Y` exactly when the Gram has a repeated
eigenvalue or is singular. -/
noncomputable def gramRowPolyRightDet (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ) :
    MvPolynomial (Fin p × Fin q) ℝ :=
  gramRowPolyRight A s * ((affMatRow A s)ᵀ * affMatRow A s).det

theorem eval_gramRowPolyRightDet (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    MvPolynomial.eval (fun rk => Y rk.1 rk.2) (gramRowPolyRightDet A s) =
      charRes ((A + rowScale s Y)ᵀ * (A + rowScale s Y)) *
        ((A + rowScale s Y)ᵀ * (A + rowScale s Y)).det := by
  simp only [gramRowPolyRightDet, gramRowPolyRight, map_mul, charRes_map, RingHom.map_det,
    RingHom.mapMatrix_apply, Matrix.map_mul, Matrix.transpose_map, affMatRow_map]

theorem gramRowPolyRightDet_ne_zero (hq : 0 < q) (hqp : q ≤ p) (A : Matrix (Fin p) (Fin q) ℝ)
    {s : Fin p → ℝ} (hs : ∀ r, s r ≠ 0) : gramRowPolyRightDet A s ≠ 0 := by
  intro h0
  have hev := eval_gramRowPolyRightDet A s
    (Matrix.of fun r k => ((witMat p q) r k - A r k) / s r)
  rw [h0, rowScale_hits hs A (witMat p q), witMat_gram hqp] at hev
  simp only [map_zero] at hev
  exact mul_ne_zero (charRes_diagonal_ne_zero hq wit_inj) wit_det_ne_zero hev.symm

theorem injective_eigenvalues₀_congr {A A' : Matrix (Fin p) (Fin p) ℝ} (h : A = A')
    (hA : A.IsHermitian) (hA' : A'.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues₀) : Function.Injective hA'.eigenvalues₀ := by
  subst h; exact hinj

theorem Ysig_eq_rowScale (τ : Fin p → ℝ) (B : Matrix (Fin p) (Fin q) ℝ) :
    (0 : Matrix (Fin p) (Fin q) ℝ) + rowScale τ B = Ysig τ B := by
  ext j k
  simp [Ysig, Matrix.diagonal_mul]

theorem simpleSpec_ae_Ysig_gram (hp : 0 < p) (hq : 0 < q) {τ : Fin p → ℝ}
    (hτ : ∀ j, τ j ≠ 0) (rk : ℕ) (hrk : rk ≤ min p q) :
    ∀ᵐ B ∂gaussianMatrix p q,
      SimpleSpec (Ysig τ B * (Ysig τ B)ᵀ) (isHermitian_mul_transpose_self (Ysig τ B)) rk := by
  rcases le_or_gt p q with hpq | hpq
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix _ (gramRowPolyLeft_ne_zero hp hpq 0 hτ)]
      with B hB
    rw [eval_gramRowPolyLeft, Ysig_eq_rowScale τ B] at hB
    exact simpleSpec_of_injective _ (injective_eigenvalues₀_of_separable _
      ((charRes_ne_zero_iff_separable hp _).mp (left_ne_zero_of_mul hB))) rk
  · have hqp : q ≤ p := hpq.le
    filter_upwards [ae_eval_ne_zero_gaussianMatrix _ (gramRowPolyRightDet_ne_zero hq hqp 0 hτ)]
      with B hB
    rw [eval_gramRowPolyRightDet, Ysig_eq_rowScale τ B] at hB
    have hinjB : Function.Injective
        (isHermitian_mul_transpose_self (Ysig τ B)ᵀ).eigenvalues₀ :=
      injective_eigenvalues₀_congr (by rw [Matrix.transpose_transpose]) _ _
        (injective_eigenvalues₀_of_separable (isHermitian_transpose_mul_self (Ysig τ B))
          ((charRes_ne_zero_iff_separable hq _).mp (left_ne_zero_of_mul hB)))
    have hdet : ((Ysig τ B)ᵀ * (Ysig τ B)ᵀᵀ).det ≠ 0 := by
      rw [Matrix.transpose_transpose]; exact right_ne_zero_of_mul hB
    have hS := simpleSpec_p_of_gram_left hqp (Ysig τ B)ᵀ hinjB hdet
    exact simpleSpec_mono (hrk.trans (min_le_right _ _))
      (EdgeGlueR.simpleSpec_congr (by rw [Matrix.transpose_transpose]) _ _ hS)

theorem simpleSpec_ae_Wsig (hp : 0 < p) (hq : 0 < q) (hd : 0 < d) {τ : Fin p → ℝ}
    (hτ : ∀ j, τ j ≠ 0) (rk : ℕ) (hrk : rk ≤ min p q) :
    ∀ᵐ B ∂gaussianMatrix p q, SimpleSpec (Wsig τ d B) (isHermitian_Wsig τ d B) rk := by
  filter_upwards [simpleSpec_ae_Ysig_gram hp hq hτ rk hrk] with B hB
  exact simpleSpec_smul (isHermitian_mul_transpose_self (Ysig τ B)) (isHermitian_Wsig τ d B)
    (inv_pos.mpr (by exact_mod_cast hd)) hB

/-! ## The mean bound at a simple index -/

theorem integrable_projMat_Wsig (τ : Fin p → ℝ) (hd : 0 < d) (a : ℕ) (j k : Fin p) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j k) (gaussianMatrix p q) :=
  Integrable.of_bound (measurable_projMat_Wsig τ hd a j k).aestronglyMeasurable 1
    (Eventually.of_forall fun B => by
      rw [Real.norm_eq_abs]; exact abs_projMat_le_one _ _ _ _ _)

theorem integral_projMat_diag_nonneg (τ : Fin p → ℝ) (a : ℕ) (j : Fin p) :
    0 ≤ ∫ B, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j ∂gaussianMatrix p q :=
  integral_nonneg fun B => by rw [projMat_diag_eq]; exact sq_nonneg _

/-- At a simple index the mean diagonal of the projector matrix sums to at most 1. -/
theorem sum_integral_projMat_diag_le_one (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    {τ : Fin p → ℝ} (hτ : ∀ j, τ j ≠ 0) {a : ℕ} (ha : a < min p q) :
    ∑ j, ∫ B, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j ∂gaussianMatrix p q ≤ 1 := by
  have hae : ∀ᵐ B ∂gaussianMatrix p q,
      ∑ j, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j ≤ 1 := by
    filter_upwards [simpleSpec_ae_Wsig hp hq hd hτ (a + 1) (Nat.succ_le_of_lt ha)] with B hB
    exact sum_projMat_diag_le_one hB (lt_of_lt_of_le ha (min_le_left _ _))
  rw [← integral_finsetSum _ (fun j _ => integrable_projMat_Wsig τ hd a j j)]
  calc ∫ B, ∑ j, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j ∂gaussianMatrix p q
      ≤ ∫ _, (1 : ℝ) ∂gaussianMatrix p q :=
        integral_mono_ae (integrable_finsetSum _ (fun j _ => integrable_projMat_Wsig τ hd a j j))
          (integrable_const 1) hae
    _ = 1 := by simp

theorem meanKAvg_projMat (τ : Fin p → ℝ) (hd : 0 < d) (a : ℕ) (J : Finset (Fin p)) :
    meanKAvg (fun B : Matrix (Fin p) (Fin q) ℝ =>
        R4C.cmat (projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a)) J =
      ((((J.card : ℝ))⁻¹ * ∑ j ∈ J, ∫ B, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j
        ∂gaussianMatrix p q : ℝ) : ℂ) := by
  have hpt : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      (J.card : ℂ)⁻¹ * ∑ j ∈ J, R4C.cmat (projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a) j j =
      ((((J.card : ℝ))⁻¹ * ∑ j ∈ J, projMat (Wsig τ d B) (isHermitian_Wsig τ d B) a j j : ℝ) :
        ℂ) := by
    intro B
    simp only [R4C.cmat, Matrix.map_apply]
    push_cast
    rfl
  rw [meanKAvg, integral_congr_ae (Eventually.of_forall hpt), integral_complex_ofReal,
    integral_const_mul, integral_finsetSum J (fun j _ => integrable_projMat_Wsig τ hd a j j)]

/-- The heteroscedastic delocalization bound. At an index `a < min p q` the mean squared
projection of a fixed vector `y` on the index eigenvector of `Wsig (tauOf w blk) d B` is at
most `∑ i, ‖y_i‖² / |block i|`, where `y_i` is the restriction of `y` to block `i`. Blocks of
size 0 contribute 0 on both sides. The weights may have any sign; only `w i ≠ 0` is used. -/
theorem lintegral_normSq_specProjIdx_Wsig_le {M : ℕ} (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0)
    (blk : Fin p → Fin M) (hd : 0 < d) {a : ℕ} (ha : a < min p q)
    (y : EuclideanSpace ℝ (Fin p)) :
    ∫⁻ B, ENNReal.ofReal
        (‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2)
        ∂gaussianMatrix p q ≤
      ENNReal.ofReal (∑ i, (∑ j ∈ blockSet blk i, y j ^ 2) / (blockSet blk i).card) := by
  have hp : 0 < p := lt_of_le_of_lt (Nat.zero_le a) (lt_of_lt_of_le ha (min_le_left _ _))
  have hq : 0 < q := lt_of_le_of_lt (Nat.zero_le a) (lt_of_lt_of_le ha (min_le_right _ _))
  have hτ : ∀ j, tauOf w blk j ≠ 0 := fun j => hw (blk j)
  have hK : BlockSym (tauOf w blk) q (fun B : Matrix (Fin p) (Fin q) ℝ =>
      R4C.cmat (projMat (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a)) 1 :=
    blockSym_specProjIdx (tauOf w blk) hd a
  have hfmeas : Measurable fun B : Matrix (Fin p) (Fin q) ℝ =>
      ‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2 :=
    ((measurable_specProjIdx_Wsig (tauOf w blk) hd a y).norm).pow_const 2
  have hfbound : ∀ B : Matrix (Fin p) (Fin q) ℝ,
      ‖‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2‖
        ≤ ‖y‖ ^ 2 := fun B => by
    rw [Real.norm_of_nonneg (sq_nonneg _)]
    exact pow_le_pow_left₀ (norm_nonneg _) (norm_specProj_le _ _ _) 2
  have hfint : Integrable (fun B : Matrix (Fin p) (Fin q) ℝ =>
      ‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2)
      (gaussianMatrix p q) :=
    Integrable.of_bound hfmeas.aestronglyMeasurable (‖y‖ ^ 2) (Eventually.of_forall hfbound)
  rw [← ofReal_integral_eq_lintegral_ofReal hfint (Eventually.of_forall fun B => sq_nonneg _)]
  refine ENNReal.ofReal_le_ofReal ?_
  have hcplx := hK.integral_qform_blocks (WithLp.ofLp y)
  have hL : ∫ B, R4C.cvec (WithLp.ofLp y) ⬝ᵥ
      (R4C.cmat (projMat (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a) *ᵥ
        R4C.cvec (WithLp.ofLp y)) ∂gaussianMatrix p q =
      ((∫ B, ‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2
        ∂gaussianMatrix p q : ℝ) : ℂ) := by
    rw [← integral_complex_ofReal]
    refine integral_congr_ae (Eventually.of_forall fun B => ?_)
    beta_reduce
    rw [R4C.dotProduct_cmat_mulVec, normSq_specProjIdx_eq_dot]
  have hreal : ∫ B, ‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B)
      a y‖ ^ 2 ∂gaussianMatrix p q =
      ∑ i, (∑ j ∈ blockSet blk i, y j ^ 2) *
        ((((blockSet blk i).card : ℝ))⁻¹ * ∑ j ∈ blockSet blk i,
          ∫ B, projMat (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a j j
            ∂gaussianMatrix p q) := by
    apply Complex.ofReal_injective
    rw [← hL, hcplx]
    simp only [meanKAvg_projMat (tauOf w blk) hd a]
    push_cast
    rfl
  have hS : ∀ i, ∑ j ∈ blockSet blk i,
      ∫ B, projMat (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a j j
        ∂gaussianMatrix p q ≤ 1 := fun i =>
    (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      (fun j _ _ => integral_projMat_diag_nonneg (tauOf w blk) a j)).trans
      (sum_integral_projMat_diag_le_one hp hq hd hτ ha)
  rw [hreal]
  refine Finset.sum_le_sum fun i _ => ?_
  rw [div_eq_mul_inv]
  refine mul_le_mul_of_nonneg_left ?_ (Finset.sum_nonneg fun j _ => sq_nonneg _)
  exact mul_le_of_le_one_right (inv_nonneg.mpr (Nat.cast_nonneg _)) (hS i)

/-- The weaker form with `‖y‖²` in front. -/
theorem lintegral_normSq_specProjIdx_Wsig_le' {M : ℕ} (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0)
    (blk : Fin p → Fin M) (hd : 0 < d) {a : ℕ} (ha : a < min p q)
    (y : EuclideanSpace ℝ (Fin p)) :
    ∫⁻ B, ENNReal.ofReal
        (‖specProjIdx (Wsig (tauOf w blk) d B) (isHermitian_Wsig (tauOf w blk) d B) a y‖ ^ 2)
        ∂gaussianMatrix p q ≤
      ENNReal.ofReal (‖y‖ ^ 2 * ∑ i, (((blockSet blk i).card : ℝ))⁻¹) := by
  refine (lintegral_normSq_specProjIdx_Wsig_le w hw blk hd ha y).trans
    (ENNReal.ofReal_le_ofReal ?_)
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ => ?_
  rw [div_eq_mul_inv]
  refine mul_le_mul_of_nonneg_right ?_ (inv_nonneg.mpr (Nat.cast_nonneg _))
  rw [EuclideanSpace.real_norm_sq_eq]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) fun j _ _ => sq_nonneg _

/-- The product form. The first component has an arbitrary law `ν`, the Gaussian block is the
second component. `MeasureTheory.lintegral_prod_le` needs no measurability and no finiteness
of `ν`. -/
theorem lintegral_prod_normSq_specProjIdx_Wsig_le {M : ℕ} (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0)
    (blk : Fin p → Fin M) (hd : 0 < d) {a : ℕ} (ha : a < min p q)
    {Ω : Type*} [MeasurableSpace Ω] (ν : Measure Ω) (y : Ω → EuclideanSpace ℝ (Fin p)) :
    ∫⁻ ξ : Ω × Matrix (Fin p) (Fin q) ℝ, ENNReal.ofReal
        (‖specProjIdx (Wsig (tauOf w blk) d ξ.2) (isHermitian_Wsig (tauOf w blk) d ξ.2) a
          (y ξ.1)‖ ^ 2) ∂(ν.prod (gaussianMatrix p q)) ≤
      (∫⁻ ω, ENNReal.ofReal (‖y ω‖ ^ 2) ∂ν) *
        ENNReal.ofReal (∑ i, (((blockSet blk i).card : ℝ))⁻¹) := by
  refine (lintegral_prod_le _).trans ?_
  rw [← lintegral_mul_const' _ _ ENNReal.ofReal_ne_top]
  refine lintegral_mono fun ω => ?_
  refine (lintegral_normSq_specProjIdx_Wsig_le' w hw blk hd ha (y ω)).trans (le_of_eq ?_)
  rw [ENNReal.ofReal_mul (sq_nonneg _)]

end HetDeloc

/-! ## Model layer: `W0hetR` and the truncated columns of `QmatHetR` -/

namespace UnalignedModelR

open HetDeloc HetStein HetR1

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- `W₀'` is `Wsig` of the Gaussian block `B` with the stacked weights, once the Gram
identity `E⊥ E⊥ᵀ = d⁻¹ B Bᵀ` holds. Rank-one mirror: `HetR1.W0het_eq_Wsig`. -/
theorem W0hetR_eq_Wsig (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {p : ℕ} (B : Matrix (Fin (∑ i, n i N)) (Fin p) ℝ)
    (hB : m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ = ((d N : ℝ))⁻¹ • (B * Bᵀ)) :
    m.W0hetR w N ω = Wsig (tauOf w (blkStack n N)) (d N) B := by
  rw [m.W0hetR_eq_of_block w N ω B hB]; rfl

/-- Column `l` of `Q = Ũ + Σ^{1/2} E V` as a vector of `EuclideanSpace`. -/
noncomputable def qColHet (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) : EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  WithLp.toLp 2 fun q => m.QmatHetR w N ω q l

/-- The column, truncated to `0` when its squared norm exceeds `C`. -/
noncomputable def qColTruncHet (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) (C : ℝ) : EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  if ‖m.qColHet w N ω l‖ ^ 2 ≤ C then m.qColHet w N ω l else 0

theorem norm_qColTruncHet_sq_le (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) {C : ℝ} (hC : 0 ≤ C) : ‖m.qColTruncHet w N ω l C‖ ^ 2 ≤ C := by
  rw [qColTruncHet]
  split_ifs with h
  · exact h
  · simpa using hC

theorem qColTruncHet_eq (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) {C : ℝ} (h : ‖m.qColHet w N ω l‖ ^ 2 ≤ C) :
    m.qColTruncHet w N ω l C = m.qColHet w N ω l :=
  if_pos h

theorem measurable_qColHet (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (l : Fin r) : Measurable fun ω => m.qColHet w N ω l := by
  have hentry : ∀ q : Fin (∑ i, n i N),
      Measurable fun X : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ => X q l := fun q =>
    (measurable_pi_apply l).comp (measurable_pi_apply q)
  exact (WithLp.measurable_toLp 2 _).comp
    (measurable_pi_lambda _ fun q => (hentry q).comp (m.measurable_QmatHetR w N))

/-- The deterministic column map: column `l` of `Ũ + Σ^{1/2} (√d)⁻¹ ZV` as a function of
the Gaussian block `ZV = Z_stack V`. -/
noncomputable def qColOfZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (l : Fin r) :
    EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  WithLp.toLp 2 fun q => (m.UtildeR w N + m.SigmaHalfR w N * ((Real.sqrt (d N))⁻¹ • ZV)) q l

/-- The deterministic truncated column map. -/
noncomputable def qColTruncOfZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (l : Fin r) (C : ℝ) :
    EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  if ‖m.qColOfZV w N ZV l‖ ^ 2 ≤ C then m.qColOfZV w N ZV l else 0

theorem qColHet_eq_ofZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) :
    m.qColHet w N ω l = m.qColOfZV w N (m.stackZG N ω * m.V N) l := by
  simp only [qColHet, qColOfZV, QmatHetR, m.stackE_eqG N ω, Matrix.smul_mul]

theorem qColTruncHet_eq_ofZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (l : Fin r) (C : ℝ) :
    m.qColTruncHet w N ω l C = m.qColTruncOfZV w N (m.stackZG N ω * m.V N) l C := by
  simp only [qColTruncHet, qColTruncOfZV, qColHet_eq_ofZV]

theorem norm_qColTruncOfZV_sq_le (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (l : Fin r) {C : ℝ} (hC : 0 ≤ C) :
    ‖m.qColTruncOfZV w N ZV l C‖ ^ 2 ≤ C := by
  rw [qColTruncOfZV]
  split_ifs with h
  · exact h
  · simpa using hC

theorem measurable_qColOfZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (l : Fin r) : Measurable fun ZV => m.qColOfZV w N ZV l := by
  refine (WithLp.measurable_toLp 2 _).comp (measurable_pi_lambda _ fun q => ?_)
  have h : (fun ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ =>
      (m.UtildeR w N + m.SigmaHalfR w N * ((Real.sqrt (d N))⁻¹ • ZV)) q l) =
      fun ZV => m.UtildeR w N q l +
        ∑ k, m.SigmaHalfR w N q k * ((Real.sqrt (d N))⁻¹ * ZV k l) := by
    funext ZV
    simp only [Matrix.add_apply, Matrix.mul_apply, Matrix.smul_apply, smul_eq_mul]
  rw [h]
  refine measurable_const.add (Finset.measurable_sum _ fun k _ => ?_)
  have hk : Measurable fun ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ => ZV k l :=
    (measurable_pi_apply l).comp (measurable_pi_apply k)
  exact (hk.const_mul _).const_mul _

theorem measurable_qColTruncOfZV (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (l : Fin r) (C : ℝ) : Measurable fun ZV => m.qColTruncOfZV w N ZV l C := by
  have hset : MeasurableSet {ZV : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ |
      ‖m.qColOfZV w N ZV l‖ ^ 2 ≤ C} :=
    ((measurable_norm.comp (m.measurable_qColOfZV w N l)).pow_const 2) measurableSet_Iic
  exact Measurable.ite hset (m.measurable_qColOfZV w N l) measurable_const

/-- **Heteroscedastic (G4), Markov form.** The double sum of the squared projections of the
truncated columns of `Q` on the top-`r` eigenvectors of `W₀'` exceeds `κ` with probability
at most `r (∑ C_l) (∑ i 1/n_i) / κ`. The weights may have any sign. -/
theorem measure_kappa_ge_le_het [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0)
    (hG : m.JointGaussianNoise) (N : ℕ) {p : ℕ} (hpd : d N = p + r)
    (hrp : r ≤ min (∑ i, n i N) p) {CC : Fin r → ℝ} (hCC : ∀ l, 0 ≤ CC l)
    {κ : ℝ} (hκ : 0 < κ) :
    μ N {ω | κ ≤ ∑ a ∈ Finset.range r, ∑ l,
        ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
          (m.qColTruncHet w N ω l (CC l))‖ ^ 2}
      ≤ ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹) / κ) := by
  classical
  obtain ⟨B, hB, hlaw⟩ := m.exists_block_hasLaw_hetR hG N hpd
  have hdN : 0 < d N := (m.tbl ⟨0, NeZero.pos M⟩).hd N
  have hSblk : ∑ i, (((blockSet (blkStack n N) i).card : ℝ))⁻¹ = ∑ i, ((n i N : ℝ))⁻¹ := by
    simp only [card_blockSet_blkStack]
  have hSnn : 0 ≤ ∑ i, ((n i N : ℝ))⁻¹ :=
    Finset.sum_nonneg fun i _ => inv_nonneg.mpr (Nat.cast_nonneg _)
  set F : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ × Matrix (Fin (∑ i, n i N)) (Fin p) ℝ → ℝ≥0∞ :=
    fun ξ => ∑ a ∈ Finset.range r, ∑ l, ENNReal.ofReal
      (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
        (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
        (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2) with hFdef
  -- measurability of `F`
  have hterm : ∀ (a : ℕ) (l : Fin r), Measurable fun ξ : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ ×
      Matrix (Fin (∑ i, n i N)) (Fin p) ℝ => ENNReal.ofReal
        (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
          (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
          (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2) := by
    intro a l
    simp_rw [normSq_specProjIdx_Wsig_eq_overlapIdx (tauOf w (blkStack n N)) hdN]
    have hg : Measurable fun ξ : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ ×
        Matrix (Fin (∑ i, n i N)) (Fin p) ℝ =>
        ((Ysig (tauOf w (blkStack n N)) ξ.2)ᵀ, m.qColTruncOfZV w N ξ.1 l (CC l)) :=
      ((measurable_Ysig_transpose _).comp measurable_snd).prodMk
        ((m.measurable_qColTruncOfZV w N l (CC l)).comp measurable_fst)
    have h2 := (measurable_overlapIdx₂ (n := p) (p := ∑ i, n i N) a).comp hg
    exact ENNReal.measurable_ofReal.comp h2
  have hFm : Measurable F := by
    rw [hFdef]
    exact Finset.measurable_sum _ fun a _ => Finset.measurable_sum _ fun l _ => hterm a l
  -- the pointwise identification
  have hpt : ∀ ω, ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
        (m.qColTruncHet w N ω l (CC l))‖ ^ 2)
      = F (m.stackZG N ω * m.V N, B ω) := by
    intro ω
    rw [hFdef]
    rw [ENNReal.ofReal_sum_of_nonneg
      (fun a _ => Finset.sum_nonneg fun l _ => sq_nonneg _)]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [ENNReal.ofReal_sum_of_nonneg (fun l _ => sq_nonneg _)]
    refine Finset.sum_congr rfl fun l _ => ?_
    congr 1
    rw [m.qColTruncHet_eq_ofZV w N ω l (CC l),
      specProjIdx_congr (m.W0hetR_eq_Wsig w N ω (B ω) (hB ω)) (m.isHermitian_W0hetR w N ω)
        (isHermitian_Wsig _ _ _) a]
  -- Markov
  have hmeas : AEMeasurable (fun ω => ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
        (m.qColTruncHet w N ω l (CC l))‖ ^ 2)) (μ N) := by
    refine AEMeasurable.congr ?_ (Filter.Eventually.of_forall fun ω => (hpt ω).symm)
    exact hFm.aemeasurable.comp_aemeasurable hlaw.aemeasurable
  have hset : {ω | κ ≤ ∑ a ∈ Finset.range r, ∑ l,
        ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
          (m.qColTruncHet w N ω l (CC l))‖ ^ 2}
      = {ω | ENNReal.ofReal κ ≤ ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
          ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
            (m.qColTruncHet w N ω l (CC l))‖ ^ 2)} := by
    ext ω
    simp only [Set.mem_ofPred_eq]
    rw [ENNReal.ofReal_le_ofReal_iff
      (Finset.sum_nonneg fun a _ => Finset.sum_nonneg fun l _ => sq_nonneg _)]
  have hκne : ENNReal.ofReal κ ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hκ]
  rw [hset]
  refine (meas_ge_le_lintegral_div hmeas hκne ENNReal.ofReal_ne_top).trans ?_
  -- the integral over the product measure
  have hint : ∫⁻ ω, ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
        (m.qColTruncHet w N ω l (CC l))‖ ^ 2) ∂(μ N)
      = ∫⁻ ξ, F ξ ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p)) := by
    rw [lintegral_congr hpt]
    exact hlaw.lintegral_comp hFm.aemeasurable
  have hbound : ∫⁻ ξ, F ξ ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p))
      ≤ ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹)) := by
    rw [hFdef]
    rw [lintegral_finsetSum _ (fun a _ => Finset.measurable_sum _ fun l _ => hterm a l)]
    have hin : ∀ a ∈ Finset.range r,
        ∫⁻ ξ, ∑ l, ENNReal.ofReal
            (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
              (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
              (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2)
            ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p))
          ≤ ENNReal.ofReal ((∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹)) := by
      intro a ha
      rw [lintegral_finsetSum _ (fun l _ => hterm a l)]
      have hak : a < min (∑ i, n i N) p := lt_of_lt_of_le (Finset.mem_range.mp ha) hrp
      have hone : ∀ l : Fin r,
          ∫⁻ ξ, ENNReal.ofReal
              (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
                (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
                (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2)
              ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p))
            ≤ ENNReal.ofReal (CC l * ∑ i, ((n i N : ℝ))⁻¹) := by
        intro l
        refine (lintegral_prod_normSq_specProjIdx_Wsig_le w hw (blkStack n N) hdN hak
          (gaussianMatrix (∑ i, n i N) r)
          (fun α => m.qColTruncOfZV w N α l (CC l))).trans ?_
        rw [hSblk]
        have hqb : ∫⁻ α, ENNReal.ofReal (‖m.qColTruncOfZV w N α l (CC l)‖ ^ 2)
            ∂(gaussianMatrix (∑ i, n i N) r) ≤ ENNReal.ofReal (CC l) := by
          calc ∫⁻ α, ENNReal.ofReal (‖m.qColTruncOfZV w N α l (CC l)‖ ^ 2)
                ∂(gaussianMatrix (∑ i, n i N) r)
              ≤ ∫⁻ _, ENNReal.ofReal (CC l) ∂(gaussianMatrix (∑ i, n i N) r) :=
                lintegral_mono fun α => ENNReal.ofReal_le_ofReal
                  (m.norm_qColTruncOfZV_sq_le w N α l (hCC l))
            _ = ENNReal.ofReal (CC l) := by simp
        calc (∫⁻ α, ENNReal.ofReal (‖m.qColTruncOfZV w N α l (CC l)‖ ^ 2)
                ∂(gaussianMatrix (∑ i, n i N) r)) * ENNReal.ofReal (∑ i, ((n i N : ℝ))⁻¹)
            ≤ ENNReal.ofReal (CC l) * ENNReal.ofReal (∑ i, ((n i N : ℝ))⁻¹) :=
              mul_le_mul' hqb le_rfl
          _ = ENNReal.ofReal (CC l * ∑ i, ((n i N : ℝ))⁻¹) := (ENNReal.ofReal_mul (hCC l)).symm
      calc ∑ l, ∫⁻ ξ, ENNReal.ofReal
              (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
                (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
                (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2)
              ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p))
          ≤ ∑ l, ENNReal.ofReal (CC l * ∑ i, ((n i N : ℝ))⁻¹) :=
            Finset.sum_le_sum fun l _ => hone l
        _ = ENNReal.ofReal ((∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹)) := by
            rw [Finset.sum_mul, ENNReal.ofReal_sum_of_nonneg
              (fun l _ => mul_nonneg (hCC l) hSnn)]
    calc ∑ a ∈ Finset.range r, ∫⁻ ξ, ∑ l, ENNReal.ofReal
            (‖specProjIdx (Wsig (tauOf w (blkStack n N)) (d N) ξ.2)
              (isHermitian_Wsig (tauOf w (blkStack n N)) (d N) ξ.2) a
              (m.qColTruncOfZV w N ξ.1 l (CC l))‖ ^ 2)
            ∂((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p))
        ≤ ∑ _a ∈ Finset.range r, ENNReal.ofReal ((∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹)) :=
          Finset.sum_le_sum hin
      _ = ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) * (∑ i, ((n i N : ℝ))⁻¹)) := by
          rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul,
            ← ENNReal.ofReal_natCast r, ← ENNReal.ofReal_mul (Nat.cast_nonneg r), mul_assoc]
  rw [hint]
  refine (ENNReal.div_le_div_right hbound _).trans (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hκ]

/-- In the regime the size condition `r ≤ min n_tot (d - r)` holds eventually. -/
theorem eventually_r_le_min [NeZero M] (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    ∀ᶠ N in atTop, r ≤ min (∑ i, n i N) (d N - r) := by
  filter_upwards [(hreg ⟨0, NeZero.pos M⟩).1.eventually_ge_atTop r,
    (hreg ⟨0, NeZero.pos M⟩).2.1.eventually_ge_atTop (2 * r)] with N hN1 hN2
  exact le_min (hN1.trans (Finset.single_le_sum (f := fun i => n i N)
    (fun i _ => Nat.zero_le _) (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))) (by omega)

end UnalignedModelR

end StackedSVD
