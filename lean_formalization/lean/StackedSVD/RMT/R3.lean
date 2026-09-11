/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RMT.R4
import StackedSVD.Prob.GaussianAdapters
import StackedSVD.Vendor.COLT83.Mathlib.Probability.SudakovFernique
import StatsMLlib.Probability.RandomMatrix.Basic

/-!
# Item R3: the upper edge of the Wishart block

Review note: `notes/archive/rmt_R3.md` (choices 9 to 16 of `notes/archive/L2_STATEMENTS.md`, all
`accept` in `notes/archive/L2_CHOICES_ANSWERS.md`). This file proves (H1), the `edge` field of
`ResolventLimits`: for every `ε > 0`,

`Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) _ ≤ bulkEdge c + ε}) atTop (𝓝 1)`.

Contents, in the order of the note.

* Deterministic: `l2_opNorm_le_frobenius`, `lipschitzWith_opNorm`, `lamMax_le_of_dotProduct_le`,
  `lamMax_gram_le_opNorm_sq`.
* R3a: `integral_iSup_bilin_le` (Sudakov-Fernique on a finite family of unit pairs) and
  `integral_opNorm_le` (Gordon, `∫ ‖Z‖ ≤ √p + √d`), through an `ε₀`-net of the two unit
  spheres and `ε₀ ↓ 0`.
* R3b: `measure_opNorm_ge_le`, one-sided Gaussian concentration at `L = 1`.
* R3c: `tendsto_measure_lamMax_le`, the assembly with `t = κ √d`.

Item R3⁻ (the lower edge) is **not** in this file: it needs R1 and item T.

The `Matrix` type carries the sup norm of the pi instance, so the operator norm enters through
the scoped instance `Matrix.Norms.L2Operator`; `‖A‖` below is always that norm.

STATUS 2026-08-30: `lean -j 3 -R . StackedSVD/RMT/R3.lean` exit 0, 0 `sorry`, no warning.
See `notes/archive/agent_reports/proof_r3.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace R3

variable {p d : ℕ}

/-! ### Deterministic step 1: the operator norm against the Frobenius norm -/

/-- **No `l2_opNorm_le_frobenius` in Mathlib v4.33.0.** Row-wise Cauchy-Schwarz gives the
bound directly. The two norms are separate scoped instances on `Matrix`, so the Frobenius
side is written as the explicit sum. -/
theorem l2_opNorm_le_frobenius (A : Matrix (Fin p) (Fin d) ℝ) :
    ‖A‖ ≤ Real.sqrt (∑ i, ∑ j, A i j ^ 2) := by
  rw [Matrix.l2_opNorm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg _) fun x => ?_
  have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) A) x
      = WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) := rfl
  rw [happ]
  have hsq : ‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin p))‖ ^ 2
      ≤ (∑ i, ∑ j, A i j ^ 2) * ‖x‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq, Finset.sum_mul]
    refine Finset.sum_le_sum fun i _ => ?_
    have hcoord : (WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin p)) i
        = ∑ j, A i j * x j := rfl
    rw [hcoord, EuclideanSpace.real_norm_sq_eq]
    exact Finset.sum_mul_sq_le_sq_mul_sq _ _ _
  have hnn : (0 : ℝ) ≤ ∑ i, ∑ j, A i j ^ 2 := by positivity
  calc ‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin p))‖
      = Real.sqrt (‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin p))‖ ^ 2) :=
        (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt ((∑ i, ∑ j, A i j ^ 2) * ‖x‖ ^ 2) := Real.sqrt_le_sqrt hsq
    _ = Real.sqrt (∑ i, ∑ j, A i j ^ 2) * ‖x‖ := by
        rw [Real.sqrt_mul hnn, Real.sqrt_sq (norm_nonneg _)]

/-- **Adapter 2.** `Z ↦ ‖Z‖_op` is `1`-Lipschitz in the entries: `|‖A‖ - ‖B‖| ≤ ‖A - B‖_op`
and `‖M‖_op ≤ ‖M‖_F`, the last being `l2_opNorm_le_frobenius` and `dist_matrixEquivE_symm_sq`. -/
theorem lipschitzWith_opNorm :
    LipschitzWith 1
      (fun w : EuclideanSpace ℝ (Fin (p * d)) => ‖(matrixEquivE p d).symm w‖) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  have h1 : dist ‖(matrixEquivE p d).symm x‖ ‖(matrixEquivE p d).symm y‖
      ≤ ‖(matrixEquivE p d).symm x - (matrixEquivE p d).symm y‖ := by
    rw [Real.dist_eq]
    exact abs_norm_sub_norm_le _ _
  have h2 : ‖(matrixEquivE p d).symm x - (matrixEquivE p d).symm y‖ ≤ dist x y := by
    rw [← matrixEquivE_symm_sub]
    calc ‖(matrixEquivE p d).symm (x - y)‖
        ≤ Real.sqrt (∑ i, ∑ j, ((matrixEquivE p d).symm (x - y)) i j ^ 2) :=
          l2_opNorm_le_frobenius _
      _ = Real.sqrt (dist x y ^ 2) := by
          have hsub : ∀ i j, (matrixEquivE p d).symm (x - y) i j
              = (matrixEquivE p d).symm x i j - (matrixEquivE p d).symm y i j := by
            intro i j; simp
          rw [dist_matrixEquivE_symm_sq]
          simp only [hsub]
      _ = dist x y := Real.sqrt_sq dist_nonneg
  simpa using h1.trans h2

/-- Measurability of the operator norm on the matrix space (the pi sigma-algebra). -/
theorem measurable_opNorm (p d : ℕ) :
    Measurable (fun Z : Matrix (Fin p) (Fin d) ℝ => ‖Z‖) := by
  have h1 : Measurable (fun w : EuclideanSpace ℝ (Fin (p * d)) =>
      ‖(matrixEquivE p d).symm w‖) := (lipschitzWith_opNorm).continuous.measurable
  have h2 := h1.comp (matrixEquivE p d).measurable
  simpa [Function.comp_def] using h2

/-! ### Deterministic step 1, continued: the Rayleigh bound -/

/-- `lamMax M ≤ b` from the Rayleigh quotient bound. `lamMax` is attained at a unit
eigenvector (`R4.exists_eigenvalues_eq_lamMax`). -/
theorem lamMax_le_of_dotProduct_le {D : ℕ} {M : Matrix (Fin D) (Fin D) ℝ} (hM : M.IsHermitian)
    {b : ℝ} (hb : 0 ≤ b) (h : ∀ x : Fin D → ℝ, x ⬝ᵥ (M *ᵥ x) ≤ b * (x ⬝ᵥ x)) :
    lamMax M hM ≤ b := by
  rcases Nat.eq_zero_or_pos D with hD | hD
  · subst hD; simpa [lamMax] using hb
  obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hM hD
  set y : Fin D → ℝ := R4.eigU hM *ᵥ Pi.single j 1 with hy
  have hcoord : (R4.eigU hM)ᵀ *ᵥ y = Pi.single j 1 := by
    rw [hy, Matrix.mulVec_mulVec, R4.transpose_eigU_mul, Matrix.one_mulVec]
  have hyy : y ⬝ᵥ y = 1 := by
    rw [← R4.dotProduct_transpose_eigU hM y, hcoord]
    simp [dotProduct, Pi.single_apply]
  have hquad : y ⬝ᵥ (M *ᵥ y) = lamMax M hM := by
    conv_lhs => rw [← R4.eigU_conj hM]
    rw [R4.dotProduct_conj, hcoord]
    simp [Pi.single_apply, hj]
  have hy2 := h y
  rw [hquad, hyy, mul_one] at hy2
  exact hy2

/-- **Rayleigh bound**, in the scaled form `W₀ = r • (Fᵀ F)` with `r = 1/d`. -/
theorem lamMax_gram_le_opNorm_sq {r : ℝ} (hr : 0 ≤ r) (F : Matrix (Fin p) (Fin d) ℝ)
    (hs : (r • (Fᵀ * F)).IsHermitian) : lamMax (r • (Fᵀ * F)) hs ≤ r * ‖F‖ ^ 2 := by
  refine lamMax_le_of_dotProduct_le hs (by positivity) fun x => ?_
  have hq : x ⬝ᵥ ((r • (Fᵀ * F)) *ᵥ x) = r * ((F *ᵥ x) ⬝ᵥ (F *ᵥ x)) := by
    rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, ← Matrix.mulVec_mulVec,
      Matrix.dotProduct_mulVec, Matrix.vecMul_transpose]
  have hop : ‖(WithLp.toLp 2 (F *ᵥ x) : EuclideanSpace ℝ (Fin p))‖
      ≤ ‖F‖ * ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin d))‖ :=
    ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) F).le_opNorm
      (WithLp.toLp 2 x)
  have hL : (F *ᵥ x) ⬝ᵥ (F *ᵥ x)
      = ‖(WithLp.toLp 2 (F *ᵥ x) : EuclideanSpace ℝ (Fin p))‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
    simp [dotProduct, sq]
  have hR : x ⬝ᵥ x = ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin d))‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
    simp [dotProduct, sq]
  rw [hq, hL, hR]
  have hnn : (0 : ℝ) ≤ ‖(WithLp.toLp 2 (F *ᵥ x) : EuclideanSpace ℝ (Fin p))‖ := norm_nonneg _
  have hnn2 : (0 : ℝ) ≤ ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin d))‖ := norm_nonneg _
  have hsq : ‖(WithLp.toLp 2 (F *ᵥ x) : EuclideanSpace ℝ (Fin p))‖ ^ 2
      ≤ (‖F‖ * ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin d))‖) ^ 2 := by
    exact pow_le_pow_left₀ hnn hop 2
  nlinarith [hsq, mul_pow ‖F‖ ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin d))‖ 2]

/-! ### R3a: Gordon's bound in expectation

The route of `notes/archive/rmt_R3.md`: Sudakov-Fernique against the comparison process
`Y_{x,y} = ⟪g, x⟫ + ⟪h, y⟫` on a finite family of unit pairs, then an `ε₀`-net of the two unit
spheres and `ε₀ ↓ 0` (choices 11, 12). -/

section Gordon

/-! #### Coordinate suprema of Gaussian images -/

/-- A coordinate of a Euclidean vector is at most its norm in absolute value. -/
theorem abs_apply_le_norm {n : Type*} [Fintype n] (z : EuclideanSpace ℝ n) (i : n) :
    |WithLp.ofLp z i| ≤ ‖z‖ := by
  have h1 : (WithLp.ofLp z i) ^ 2 ≤ ‖z‖ ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
    exact Finset.single_le_sum (fun j _ => sq_nonneg _) (Finset.mem_univ i)
  rw [← Real.sqrt_sq_eq_abs]
  calc Real.sqrt ((WithLp.ofLp z i) ^ 2) ≤ Real.sqrt (‖z‖ ^ 2) := Real.sqrt_le_sqrt h1
    _ = ‖z‖ := Real.sqrt_sq (norm_nonneg _)

/-- The coordinate supremum of a Euclidean vector is at most its norm in absolute value. -/
theorem abs_ciSup_apply_le {n : Type*} [Fintype n] [Nonempty n] (z : EuclideanSpace ℝ n) :
    |⨆ i, WithLp.ofLp z i| ≤ ‖z‖ := by
  obtain ⟨i₀⟩ := ‹Nonempty n›
  refine abs_le.mpr ⟨?_, ciSup_le fun i => (le_abs_self _).trans (abs_apply_le_norm z i)⟩
  have h0 : -‖z‖ ≤ WithLp.ofLp z i₀ := (abs_le.mp (abs_apply_le_norm z i₀)).1
  exact h0.trans (le_ciSup (Finite.bddAbove_range fun i => WithLp.ofLp z i) i₀)

/-- The coordinate supremum of the image of a Gaussian measure by a linear map is integrable:
it is at most `‖L‖ ‖w‖`, and a Gaussian measure integrates linear growth. -/
theorem integrable_ciSup_clm {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    [FiniteDimensional ℝ E] [MeasurableSpace E] [BorelSpace E]
    {n : Type*} [Finite n] [Nonempty n] (L : E →L[ℝ] EuclideanSpace ℝ n)
    (ν : Measure E) [IsGaussian ν] :
    Integrable (fun w => ⨆ i, WithLp.ofLp (L w) i) ν := by
  have := Fintype.ofFinite n
  have hm : Measurable (fun w => ⨆ i, WithLp.ofLp (L w) i) :=
    measurable_ciSup_apply.comp L.continuous.measurable
  refine IsGaussian.integrable_of_abs_le_add_mul_norm hm.aestronglyMeasurable
    (A := 0) (B := ‖L‖) fun w => ?_
  calc |⨆ i, WithLp.ofLp (L w) i| ≤ ‖L w‖ := abs_ciSup_apply_le _
    _ ≤ ‖L‖ * ‖w‖ := L.le_opNorm w
    _ = 0 + ‖L‖ * ‖w‖ := by ring

/-! #### The flattening with index type `Fin p × Fin d` -/

/-- The matrix space flattened to `EuclideanSpace ℝ (Fin p × Fin d)`. R3 needs this index type,
not the `Fin (p * d)` of `matrixEquivE`, because the covariance of the bilinear process factors
over the two blocks. -/
def matrixEquivP (p d : ℕ) : Matrix (Fin p) (Fin d) ℝ ≃ᵐ EuclideanSpace ℝ (Fin p × Fin d) :=
  (matrixUncurry p d).trans (MeasurableEquiv.toLp 2 (Fin p × Fin d → ℝ))

@[simp]
theorem matrixEquivP_apply (Z : Matrix (Fin p) (Fin d) ℝ) (q : Fin p × Fin d) :
    WithLp.ofLp (matrixEquivP p d Z) q = Z q.1 q.2 := rfl

theorem measurePreserving_matrixEquivP (p d : ℕ) :
    MeasurePreserving (matrixEquivP p d) (gaussianMatrix p d)
      (stdGaussian (EuclideanSpace ℝ (Fin p × Fin d))) := by
  have h1 : MeasurePreserving
      (WithLp.toLp 2 : (Fin p × Fin d → ℝ) → EuclideanSpace ℝ (Fin p × Fin d))
      (Measure.pi fun _ : Fin p × Fin d => gaussianReal 0 1)
      (stdGaussian (EuclideanSpace ℝ (Fin p × Fin d))) :=
    ⟨WithLp.measurable_toLp _ _, map_pi_eq_stdGaussian⟩
  exact h1.comp (measurePreserving_matrixUncurry p d)

/-! #### The two covariance matrices -/

variable {ι : Type*}

/-- The matrix of the bilinear process `X_i = x_iᵀ Z y_i` as a linear form in the entries. -/
def bilMat (x : ι → EuclideanSpace ℝ (Fin p)) (y : ι → EuclideanSpace ℝ (Fin d)) :
    Matrix ι (Fin p × Fin d) ℝ :=
  Matrix.of fun i q => WithLp.ofLp (x i) q.1 * WithLp.ofLp (y i) q.2

/-- The matrix of the comparison process `Y_i = ⟪g, x_i⟫ + ⟪h, y_i⟫`. -/
def cmpMat (x : ι → EuclideanSpace ℝ (Fin p)) (y : ι → EuclideanSpace ℝ (Fin d)) :
    Matrix ι (Fin p ⊕ Fin d) ℝ :=
  Matrix.of fun i k => Sum.elim (fun a => WithLp.ofLp (x i) a) (fun b => WithLp.ofLp (y i) b) k

section GramSection

variable (x : ι → EuclideanSpace ℝ (Fin p)) (y : ι → EuclideanSpace ℝ (Fin d))

/-- `∑ a, u a v a` is the real inner product. -/
theorem sum_mul_eq_inner {n : ℕ} (u v : EuclideanSpace ℝ (Fin n)) :
    ∑ a, WithLp.ofLp u a * WithLp.ofLp v a = inner ℝ u v := by
  rw [inner_euclidean_eq_dotProduct]
  rfl

theorem sum_mul_self_eq_norm_sq {n : ℕ} (u : EuclideanSpace ℝ (Fin n)) :
    ∑ a, WithLp.ofLp u a * WithLp.ofLp u a = ‖u‖ ^ 2 := by
  rw [EuclideanSpace.real_norm_sq_eq]
  exact Finset.sum_congr rfl fun a _ => by ring

/-- `Cov(X_i, X_j) = ⟪x_i, x_j⟫ ⟪y_i, y_j⟫`. -/
theorem bilMat_gram (i j : ι) :
    (bilMat x y * (bilMat x y)ᵀ) i j
      = (∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x j) a)
        * (∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y j) b) := by
  rw [Matrix.mul_apply, Fintype.sum_mul_sum, Fintype.sum_prod_type]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by
    simp only [bilMat, Matrix.of_apply, Matrix.transpose_apply]; ring

/-- `Cov(Y_i, Y_j) = ⟪x_i, x_j⟫ + ⟪y_i, y_j⟫`. -/
theorem cmpMat_gram (i j : ι) :
    (cmpMat x y * (cmpMat x y)ᵀ) i j
      = (∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x j) a)
        + (∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y j) b) := by
  rw [Matrix.mul_apply, Fintype.sum_sum_type]
  simp [cmpMat]

/-- The bilinear process, realized on the flattened matrix space. -/
theorem bilMat_apply (Z : Matrix (Fin p) (Fin d) ℝ) (i : ι) :
    WithLp.ofLp (Matrix.toEuclideanLin (bilMat x y) (matrixEquivP p d Z)) i
      = WithLp.ofLp (x i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (y i)) := by
  change (bilMat x y *ᵥ WithLp.ofLp (matrixEquivP p d Z)) i
      = WithLp.ofLp (x i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (y i))
  simp only [Matrix.mulVec, dotProduct, bilMat, Matrix.of_apply, matrixEquivP_apply]
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun b _ => by ring

/-- The comparison process, split into its two blocks. -/
theorem cmpMat_apply (w : EuclideanSpace ℝ (Fin p ⊕ Fin d)) (i : ι) :
    WithLp.ofLp (Matrix.toEuclideanLin (cmpMat x y) w) i
      = (∑ a, WithLp.ofLp (x i) a * WithLp.ofLp w (Sum.inl a))
        + ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b) := by
  change (cmpMat x y *ᵥ WithLp.ofLp w) i = _
  simp only [Matrix.mulVec, dotProduct, cmpMat, Matrix.of_apply]
  rw [Fintype.sum_sum_type]
  simp

variable [Finite ι]

theorem posSemidef_bilMat_gram : (bilMat x y * (bilMat x y)ᵀ).PosSemidef := by
  simpa using Matrix.posSemidef_self_mul_conjTranspose (bilMat x y)

theorem posSemidef_cmpMat_gram : (cmpMat x y * (cmpMat x y)ᵀ).PosSemidef := by
  simpa using Matrix.posSemidef_self_mul_conjTranspose (cmpMat x y)

end GramSection

/-! #### The two coordinate blocks of `ℝ^{p} ⊕ ℝ^{d}` -/

/-- The first block projection, as a rectangular matrix. -/
def projL (p d : ℕ) : Matrix (Fin p) (Fin p ⊕ Fin d) ℝ := Matrix.fromCols 1 0

/-- The second block projection, as a rectangular matrix. -/
def projR (p d : ℕ) : Matrix (Fin d) (Fin p ⊕ Fin d) ℝ := Matrix.fromCols 0 1

theorem projL_mul_transpose (p d : ℕ) : projL p d * (projL p d)ᵀ = 1 := by
  rw [projL, Matrix.transpose_fromCols, Matrix.fromCols_mul_fromRows]
  simp

theorem projR_mul_transpose (p d : ℕ) : projR p d * (projR p d)ᵀ = 1 := by
  rw [projR, Matrix.transpose_fromCols, Matrix.fromCols_mul_fromRows]
  simp

theorem projL_apply (w : EuclideanSpace ℝ (Fin p ⊕ Fin d)) (a : Fin p) :
    WithLp.ofLp (Matrix.toEuclideanLin (projL p d) w) a = WithLp.ofLp w (Sum.inl a) := by
  change (projL p d *ᵥ WithLp.ofLp w) a = _
  simp only [Matrix.mulVec, dotProduct, projL]
  rw [Fintype.sum_sum_type]
  simp [Matrix.one_apply]

theorem projR_apply (w : EuclideanSpace ℝ (Fin p ⊕ Fin d)) (b : Fin d) :
    WithLp.ofLp (Matrix.toEuclideanLin (projR p d) w) b = WithLp.ofLp w (Sum.inr b) := by
  change (projR p d *ᵥ WithLp.ofLp w) b = _
  simp only [Matrix.mulVec, dotProduct, projR]
  rw [Fintype.sum_sum_type]
  simp [Matrix.one_apply]

theorem integrable_norm_toEuclideanLin {n : ℕ} (M : Matrix (Fin n) (Fin p ⊕ Fin d) ℝ) :
    Integrable (fun w : EuclideanSpace ℝ (Fin p ⊕ Fin d) => ‖Matrix.toEuclideanLin M w‖)
      (stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) := by
  set L := LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin M) with hL
  refine IsGaussian.integrable_of_abs_le_add_mul_norm
    (L.continuous.norm.measurable.aestronglyMeasurable) (A := 0) (B := ‖L‖) fun w => ?_
  rw [abs_of_nonneg (norm_nonneg _), zero_add]
  exact L.le_opNorm w

/-- `E ‖g‖ ≤ √p` for the first block of a standard Gaussian vector of `ℝ^{p} ⊕ ℝ^{d}`. -/
theorem integral_norm_projL_le (p d : ℕ) :
    ∫ w, ‖Matrix.toEuclideanLin (projL p d) w‖
        ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) ≤ Real.sqrt p := by
  have hmap : (stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (projL p d)))
      = stdGaussian (EuclideanSpace ℝ (Fin p)) := by
    rw [← multivariateGaussian_zero_one (ι := Fin p ⊕ Fin d),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one,
      projL_mul_transpose, multivariateGaussian_zero_one]
  have h : ∫ z, ‖z‖ ∂((stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))).map
        (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (projL p d))))
      = ∫ w, ‖Matrix.toEuclideanLin (projL p d) w‖
          ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) :=
    integral_map (by fun_prop) (by fun_prop)
  rw [← h, hmap]
  simpa using integral_norm_stdGaussian_le (ι := Fin p)

/-- `E ‖h‖ ≤ √d` for the second block. -/
theorem integral_norm_projR_le (p d : ℕ) :
    ∫ w, ‖Matrix.toEuclideanLin (projR p d) w‖
        ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) ≤ Real.sqrt d := by
  have hmap : (stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (projR p d)))
      = stdGaussian (EuclideanSpace ℝ (Fin d)) := by
    rw [← multivariateGaussian_zero_one (ι := Fin p ⊕ Fin d),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one,
      projR_mul_transpose, multivariateGaussian_zero_one]
  have h : ∫ z, ‖z‖ ∂((stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))).map
        (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (projR p d))))
      = ∫ w, ‖Matrix.toEuclideanLin (projR p d) w‖
          ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) :=
    integral_map (by fun_prop) (by fun_prop)
  rw [← h, hmap]
  simpa using integral_norm_stdGaussian_le (ι := Fin d)

/-! #### R3a proper -/

/-- **R3a**, finite family form: Sudakov-Fernique on unit pairs. -/
theorem integral_iSup_bilin_le {ι : Type*} [Finite ι] [Nonempty ι]
    (x : ι → EuclideanSpace ℝ (Fin p)) (y : ι → EuclideanSpace ℝ (Fin d))
    (hx : ∀ i, ‖x i‖ = 1) (hy : ∀ i, ‖y i‖ = 1) :
    ∫ Z, ⨆ i, WithLp.ofLp (x i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (y i)) ∂(gaussianMatrix p d)
      ≤ Real.sqrt p + Real.sqrt d := by
  classical
  have := Fintype.ofFinite ι
  -- the increment inequality `2 - 2 a b ≤ 4 - 2 a - 2 b`
  have hincr : ∀ i j, (bilMat x y * (bilMat x y)ᵀ) i i + (bilMat x y * (bilMat x y)ᵀ) j j
      - 2 * (bilMat x y * (bilMat x y)ᵀ) i j
      ≤ (cmpMat x y * (cmpMat x y)ᵀ) i i + (cmpMat x y * (cmpMat x y)ᵀ) j j
        - 2 * (cmpMat x y * (cmpMat x y)ᵀ) i j := by
    intro i j
    simp only [bilMat_gram, cmpMat_gram]
    have hxi : ∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x i) a = 1 := by
      rw [sum_mul_self_eq_norm_sq, hx i]; norm_num
    have hxj : ∑ a, WithLp.ofLp (x j) a * WithLp.ofLp (x j) a = 1 := by
      rw [sum_mul_self_eq_norm_sq, hx j]; norm_num
    have hyi : ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y i) b = 1 := by
      rw [sum_mul_self_eq_norm_sq, hy i]; norm_num
    have hyj : ∑ b, WithLp.ofLp (y j) b * WithLp.ofLp (y j) b = 1 := by
      rw [sum_mul_self_eq_norm_sq, hy j]; norm_num
    have hxij : |∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x j) a| ≤ 1 := by
      rw [sum_mul_eq_inner]
      calc |inner ℝ (x i) (x j)| ≤ ‖x i‖ * ‖x j‖ := abs_real_inner_le_norm _ _
        _ = 1 := by rw [hx i, hx j]; norm_num
    have hyij : |∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y j) b| ≤ 1 := by
      rw [sum_mul_eq_inner]
      calc |inner ℝ (y i) (y j)| ≤ ‖y i‖ * ‖y j‖ := abs_real_inner_le_norm _ _
        _ = 1 := by rw [hy i, hy j]; norm_num
    rw [hxi, hxj, hyi, hyj]
    obtain ⟨h1, h2⟩ := abs_le.mp hxij
    obtain ⟨h3, h4⟩ := abs_le.mp hyij
    nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 - ∑ a, WithLp.ofLp (x i) a * WithLp.ofLp (x j) a)
      (by linarith : (0:ℝ) ≤ 1 - ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp (y j) b)]
  -- the left side, as an integral against `N(0, S)`
  have hSmap : (stdGaussian (EuclideanSpace ℝ (Fin p × Fin d))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (bilMat x y)))
      = multivariateGaussian 0 (bilMat x y * (bilMat x y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin p × Fin d),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hleft : ∫ Z, ⨆ i, WithLp.ofLp (x i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (y i)) ∂(gaussianMatrix p d)
      = ∫ z, ⨆ i, WithLp.ofLp z i ∂(multivariateGaussian 0 (bilMat x y * (bilMat x y)ᵀ)) := by
    have h1 : ∫ Z, ⨆ i, WithLp.ofLp (x i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (y i)) ∂(gaussianMatrix p d)
        = ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (bilMat x y) w) i
            ∂(stdGaussian (EuclideanSpace ℝ (Fin p × Fin d))) := by
      rw [← (measurePreserving_matrixEquivP p d).integral_comp'
        (fun w => ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (bilMat x y) w) i)]
      exact integral_congr_ae (Eventually.of_forall fun Z =>
        iSup_congr fun i => (bilMat_apply x y Z i).symm)
    rw [h1, ← hSmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    rfl
  -- the right side
  have hTmap : (stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (cmpMat x y)))
      = multivariateGaussian 0 (cmpMat x y * (cmpMat x y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin p ⊕ Fin d),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hright : ∫ z, ⨆ i, WithLp.ofLp z i
        ∂(multivariateGaussian 0 (cmpMat x y * (cmpMat x y)ᵀ)) ≤ Real.sqrt p + Real.sqrt d := by
    rw [← hTmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    have hpt : ∀ w : EuclideanSpace ℝ (Fin p ⊕ Fin d),
        (⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (cmpMat x y) w) i)
          ≤ ‖Matrix.toEuclideanLin (projL p d) w‖
            + ‖Matrix.toEuclideanLin (projR p d) w‖ := by
      intro w
      refine ciSup_le fun i => ?_
      rw [cmpMat_apply]
      have h1 : ∑ a, WithLp.ofLp (x i) a * WithLp.ofLp w (Sum.inl a)
          ≤ ‖Matrix.toEuclideanLin (projL p d) w‖ := by
        have heq : ∑ a, WithLp.ofLp (x i) a * WithLp.ofLp w (Sum.inl a)
            = inner ℝ (x i) (Matrix.toEuclideanLin (projL p d) w) := by
          rw [← sum_mul_eq_inner]
          exact Finset.sum_congr rfl fun a _ => by rw [projL_apply]
        rw [heq]
        calc inner ℝ (x i) (Matrix.toEuclideanLin (projL p d) w)
            ≤ ‖x i‖ * ‖Matrix.toEuclideanLin (projL p d) w‖ := real_inner_le_norm _ _
          _ = ‖Matrix.toEuclideanLin (projL p d) w‖ := by rw [hx i, one_mul]
      have h2 : ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b)
          ≤ ‖Matrix.toEuclideanLin (projR p d) w‖ := by
        have heq : ∑ b, WithLp.ofLp (y i) b * WithLp.ofLp w (Sum.inr b)
            = inner ℝ (y i) (Matrix.toEuclideanLin (projR p d) w) := by
          rw [← sum_mul_eq_inner]
          exact Finset.sum_congr rfl fun b _ => by rw [projR_apply]
        rw [heq]
        calc inner ℝ (y i) (Matrix.toEuclideanLin (projR p d) w)
            ≤ ‖y i‖ * ‖Matrix.toEuclideanLin (projR p d) w‖ := real_inner_le_norm _ _
          _ = ‖Matrix.toEuclideanLin (projR p d) w‖ := by rw [hy i, one_mul]
      linarith
    calc ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (cmpMat x y) w) i
          ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d)))
        ≤ ∫ w, (‖Matrix.toEuclideanLin (projL p d) w‖
            + ‖Matrix.toEuclideanLin (projR p d) w‖)
            ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) := by
          refine integral_mono
            (integrable_ciSup_clm
              (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (cmpMat x y))) _)
            (((integrable_norm_toEuclideanLin (projL p d)).add
              (integrable_norm_toEuclideanLin (projR p d)))) hpt
      _ = (∫ w, ‖Matrix.toEuclideanLin (projL p d) w‖
              ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))))
            + ∫ w, ‖Matrix.toEuclideanLin (projR p d) w‖
              ∂(stdGaussian (EuclideanSpace ℝ (Fin p ⊕ Fin d))) :=
          integral_add (integrable_norm_toEuclideanLin (projL p d))
            (integrable_norm_toEuclideanLin (projR p d))
      _ ≤ Real.sqrt p + Real.sqrt d :=
          add_le_add (integral_norm_projL_le p d) (integral_norm_projR_le p d)
  rw [hleft]
  exact (sudakov_fernique (posSemidef_bilMat_gram x y) (posSemidef_cmpMat_gram x y) hincr).trans
    hright

/-! #### The net argument -/

/-- A finite `e`-net of the unit sphere, whose points lie **on** the sphere (choice 11). -/
theorem exists_sphere_enet (n : ℕ) (hn : 0 < n) {e : ℝ} (he : 0 < e) :
    ∃ t : Finset (EuclideanSpace ℝ (Fin n)),
      IsENet t e (RMT.euclideanUnitSphere n) ∧
        (↑t : Set (EuclideanSpace ℝ (Fin n))) ⊆ RMT.euclideanUnitSphere n := by
  obtain ⟨t, ht, hsub, -⟩ := exists_enet_subset_from_half he
    ((isCompact_sphere (0 : EuclideanSpace ℝ (Fin n)) 1).totallyBounded)
    (RMT.euclideanUnitSphere_nonempty_of_pos hn)
  exact ⟨t, ht, hsub⟩

/-- The `ε₀`-step of Gordon's bound. -/
theorem integral_opNorm_le_aux (hp : 0 < p) (hd : 0 < d) {e : ℝ} (he0 : 0 < e)
    (he : 2 * e < 1) :
    ∫ Z, ‖Z‖ ∂(gaussianMatrix p d) ≤ (Real.sqrt p + Real.sqrt d) / (1 - 2 * e) := by
  classical
  obtain ⟨Nx, hNx, hsubx⟩ := exists_sphere_enet p hp he0
  obtain ⟨Ny, hNy, hsuby⟩ := exists_sphere_enet d hd he0
  set Nnet : RMT.CenteredMatrixBilinearNet p d e :=
    { domainNet := Ny
      codomainNet := Nx
      domain_isNet := hNy
      codomain_isNet := hNx
      domain_subset_unitBall := hsuby.trans Metric.sphere_subset_closedBall
      codomain_subset_unitBall := hsubx.trans Metric.sphere_subset_closedBall } with hNnet
  obtain ⟨vx, hvx⟩ : Nx.Nonempty := Nnet.toMatrixBilinearNet.codomainNet_nonempty_of_pos hp
  obtain ⟨vy, hvy⟩ : Ny.Nonempty := Nnet.toMatrixBilinearNet.domainNet_nonempty_of_pos hd
  have : Nonempty (↥Nx × ↥Ny × Bool) := ⟨(⟨vx, hvx⟩, ⟨vy, hvy⟩, true)⟩
  set xf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin p) :=
    fun i => if i.2.2 then (i.1 : EuclideanSpace ℝ (Fin p))
      else -(i.1 : EuclideanSpace ℝ (Fin p)) with hxf
  set yf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin d) :=
    fun i => (i.2.1 : EuclideanSpace ℝ (Fin d)) with hyf
  have hxn : ∀ i, ‖xf i‖ = 1 := by
    intro i
    have h1 : ‖(i.1 : EuclideanSpace ℝ (Fin p))‖ = 1 :=
      RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsubx i.1.2)
    by_cases hb : i.2.2 = true <;> simp [hxf, hb, h1]
  have hyn : ∀ i, ‖yf i‖ = 1 := fun i =>
    RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsuby i.2.1.2)
  set u : Matrix (Fin p) (Fin d) ℝ → ℝ :=
    fun Z => ⨆ i, WithLp.ofLp (xf i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (yf i)) with hu
  have hbdd : ∀ Z, BddAbove (Set.range fun i => WithLp.ofLp (xf i) ⬝ᵥ (Z *ᵥ WithLp.ofLp (yf i))) :=
    fun Z => Finite.bddAbove_range _
  -- every net value is dominated by `u Z`
  have hnet : ∀ Z : Matrix (Fin p) (Fin d) ℝ, ∀ a ∈ Nnet.domainNet, ∀ b ∈ Nnet.codomainNet,
      |inner ℝ (Matrix.toEuclideanLin Z a) b| ≤ u Z := by
    intro Z a ha b hb
    have hinner : inner ℝ (Matrix.toEuclideanLin Z a) b
        = WithLp.ofLp b ⬝ᵥ (Z *ᵥ WithLp.ofLp a) := by
      rw [inner_euclidean_eq_dotProduct, dotProduct_comm]
      rfl
    have hpos := le_ciSup (hbdd Z) ((⟨b, hb⟩ : ↥Nx), (⟨a, ha⟩ : ↥Ny), true)
    have hneg := le_ciSup (hbdd Z) ((⟨b, hb⟩ : ↥Nx), (⟨a, ha⟩ : ↥Ny), false)
    simp only [hxf, hyf, if_true] at hpos hneg
    rw [hinner]
    refine abs_le.mpr ⟨?_, hpos⟩
    have hneg' : -(WithLp.ofLp b ⬝ᵥ (Z *ᵥ WithLp.ofLp a)) ≤ u Z := by
      refine le_trans (le_of_eq ?_) hneg
      simp [neg_dotProduct]
    linarith
  have hu0 : ∀ Z, 0 ≤ u Z := by
    intro Z
    have := hnet Z vy hvy vx hvx
    exact le_trans (abs_nonneg _) this
  -- the net bound, pointwise in `Z`
  have hbound : ∀ Z : Matrix (Fin p) (Fin d) ℝ, ‖Z‖ ≤ u Z / (1 - 2 * e) := fun Z =>
    RMT.matrixOperatorNorm_le_of_centered_bilinear_net Z Nnet he0.le he (hu0 Z) (hnet Z)
  -- integrability of the process
  have hint : Integrable u (gaussianMatrix p d) := by
    have hg : Integrable
        (fun w : EuclideanSpace ℝ (Fin p × Fin d) =>
          ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
            (Matrix.toEuclideanLin (bilMat xf yf)) w) i)
        (stdGaussian (EuclideanSpace ℝ (Fin p × Fin d))) :=
      integrable_ciSup_clm _ _
    have heq : u = (fun w : EuclideanSpace ℝ (Fin p × Fin d) =>
        ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
          (Matrix.toEuclideanLin (bilMat xf yf)) w) i) ∘ (matrixEquivP p d) := by
      funext Z
      exact (iSup_congr fun i => bilMat_apply xf yf Z i).symm
    rw [heq]
    exact ((measurePreserving_matrixEquivP p d).integrable_comp_emb
      (matrixEquivP p d).measurableEmbedding).mpr hg
  -- integrate
  have hden : (0 : ℝ) < 1 - 2 * e := by linarith
  calc ∫ Z, ‖Z‖ ∂(gaussianMatrix p d)
      ≤ ∫ Z, u Z / (1 - 2 * e) ∂(gaussianMatrix p d) :=
        integral_mono_of_nonneg (Eventually.of_forall fun Z => norm_nonneg _)
          (hint.div_const _) (Eventually.of_forall hbound)
    _ = (∫ Z, u Z ∂(gaussianMatrix p d)) / (1 - 2 * e) := integral_div _ _
    _ ≤ (Real.sqrt p + Real.sqrt d) / (1 - 2 * e) := by
        exact div_le_div_of_nonneg_right (integral_iSup_bilin_le xf yf hxn hyn) hden.le

/-- **R3a.** Gordon's bound in expectation, after `ε₀ ↓ 0` (choice 12). -/
theorem integral_opNorm_le (hp : 0 < p) (hd : 0 < d) :
    ∫ Z, ‖Z‖ ∂(gaussianMatrix p d) ≤ Real.sqrt p + Real.sqrt d := by
  have hcont : ContinuousAt
      (fun e : ℝ => (Real.sqrt p + Real.sqrt d) / (1 - 2 * e)) 0 := by
    refine ContinuousAt.div continuousAt_const (by fun_prop) (by norm_num)
  have hf : Tendsto (fun e : ℝ => (Real.sqrt p + Real.sqrt d) / (1 - 2 * e))
      (𝓝[>] (0 : ℝ)) (𝓝 (Real.sqrt p + Real.sqrt d)) := by
    have h : Tendsto (fun e : ℝ => (Real.sqrt p + Real.sqrt d) / (1 - 2 * e))
        (𝓝[>] (0 : ℝ)) (𝓝 ((Real.sqrt p + Real.sqrt d) / (1 - 2 * (0 : ℝ)))) :=
      hcont.continuousWithinAt
    simpa using h
  refine ge_of_tendsto hf ?_
  have hsmall : ∀ᶠ e in 𝓝[>] (0 : ℝ), e ∈ Set.Iio (1 / 2 : ℝ) :=
    Filter.Eventually.filter_mono nhdsWithin_le_nhds (Iio_mem_nhds (by norm_num))
  filter_upwards [self_mem_nhdsWithin, hsmall] with e he1 he2
  have he3 : e < 1 / 2 := he2
  exact integral_opNorm_le_aux hp hd he1 (by linarith)

end Gordon

/-! ### R3b: one-sided Gaussian concentration -/

/-- **R3b.** One-sided Gaussian concentration of the operator norm. It is
`measure_ge_le_of_lipschitz` at `LL = 1` and `F := fun w => ‖(matrixEquivE p d).symm w‖`,
for which `F (matrixEquivE p d Z) = ‖Z‖`. -/
theorem measure_opNorm_ge_le (hp : 0 < p) (hd : 0 < d) {t : ℝ} (ht : 0 < t) :
    ((gaussianMatrix p d)
        {Z | (∫ Z', ‖Z'‖ ∂(gaussianMatrix p d)) + t ≤ ‖Z‖}).toReal
      ≤ Real.exp (-t ^ 2 / 2) := by
  have h := measure_ge_le_of_lipschitz (p := p) (d := d) (LL := 1)
    (Nat.mul_pos hp hd) (by norm_num) (lipschitzWith_opNorm (p := p) (d := d)) ht
  simp only [MeasurableEquiv.symm_apply_apply] at h
  simpa using h

/-! ### R3c: the assembly -/

/-- The scalar inequality of R3c: `(r + 1 + κ)² ≤ (1 + √c)² + ε`. -/
theorem edge_arith {sc κ ε r : ℝ} (hsc : 0 ≤ sc) (hκ0 : 0 < κ) (hκ1 : κ ≤ 1)
    (hκK : κ * (4 * (1 + sc) + 4) ≤ ε) (hr0 : 0 ≤ r) (hr : r ≤ sc + κ) :
    (r + 1 + κ) ^ 2 ≤ (1 + sc) ^ 2 + ε := by
  have h1 : r + 1 + κ ≤ sc + 1 + 2 * κ := by linarith
  have h2 : (0 : ℝ) ≤ r + 1 + κ := by linarith
  have h3 : (r + 1 + κ) ^ 2 ≤ (sc + 1 + 2 * κ) ^ 2 := by nlinarith
  have h4 : κ ^ 2 ≤ κ := by nlinarith
  nlinarith [h3, h4, hκK, hsc, hκ0.le]

/-- **R3 (H1).** The form `notes/archive/rmt_R5.md` consumes, with the block `B` of
`notes/archive/rmt_R0.md` (`exists_block_hasLaw_snd`: `W₀ = d⁻¹ Bᵀ B`, `B ~ gaussianMatrix p d`). -/
theorem tendsto_measure_lamMax_le
    {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] (μ : ∀ N, Measure (Ω N))
    [∀ N, IsProbabilityMeasure (μ N)] {p d : ℕ → ℕ} {c : ℝ} (hc : 0 < c)
    (B : (N : ℕ) → Ω N → Matrix (Fin (p N)) (Fin (d N)) ℝ)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (d N)) (Fin (d N)) ℝ)
    (hW₀ : ∀ N ω, W₀ N ω = ((d N : ℝ))⁻¹ • ((B N ω)ᵀ * B N ω))
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (hlaw : ∀ N, HasLaw (B N) (gaussianMatrix (p N) (d N)) (μ N))
    (hp : ∀ N, 0 < p N) (hd : ∀ N, 0 < d N) (hdtop : Tendsto d atTop atTop)
    (hcN : Tendsto (fun N => (p N : ℝ) / d N) atTop (𝓝 c)) {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  classical
  set sc : ℝ := Real.sqrt c with hscdef
  have hsc0 : (0 : ℝ) ≤ sc := Real.sqrt_nonneg _
  set K : ℝ := 4 * (1 + sc) + 4 with hKdef
  have hK0 : (0 : ℝ) < K := by rw [hKdef]; nlinarith
  set κ : ℝ := min 1 (ε / K) with hkapdef
  have hκ0 : (0 : ℝ) < κ := lt_min one_pos (div_pos hε hK0)
  have hκ1 : κ ≤ 1 := min_le_left _ _
  have hκK : κ * K ≤ ε := by
    have hle : κ ≤ ε / K := min_le_right _ _
    have := mul_le_mul_of_nonneg_right hle hK0.le
    rwa [div_mul_cancel₀ _ hK0.ne'] at this
  -- the tail parameter and the bad event
  set t : ℕ → ℝ := fun N => κ * Real.sqrt (d N) with htdef
  have hsd : ∀ N, (0 : ℝ) < Real.sqrt (d N) := fun N =>
    Real.sqrt_pos.mpr (by exact_mod_cast hd N)
  have ht0 : ∀ N, (0 : ℝ) < t N := fun N => by
    rw [htdef]; exact mul_pos hκ0 (hsd N)
  set Bad : (N : ℕ) → Set (Matrix (Fin (p N)) (Fin (d N)) ℝ) := fun N =>
    {Z | (∫ Z', ‖Z'‖ ∂(gaussianMatrix (p N) (d N))) + t N ≤ ‖Z‖} with hBaddef
  have hmeasBad : ∀ N, MeasurableSet (Bad N) := fun N =>
    measurableSet_le measurable_const (measurable_opNorm _ _)
  -- transfer of the complement to the canonical Gaussian law
  have hmeasure : ∀ N, μ N {ω | B N ω ∈ (Bad N)ᶜ}
      = (gaussianMatrix (p N) (d N)) ((Bad N)ᶜ) := fun N =>
    (hlaw N).measure_eq (hmeasBad N).compl
  have hcompl : ∀ N, (gaussianMatrix (p N) (d N)) ((Bad N)ᶜ)
      = 1 - (gaussianMatrix (p N) (d N)) (Bad N) := fun N =>
    prob_compl_eq_one_sub (hmeasBad N)
  -- the tail bound
  have hbound : ∀ N, (gaussianMatrix (p N) (d N)) (Bad N)
      ≤ ENNReal.ofReal (Real.exp (-(κ ^ 2 * (d N : ℝ)) / 2)) := by
    intro N
    have h := measure_opNorm_ge_le (p := p N) (d := d N) (hp N) (hd N) (ht0 N)
    have hsq : (t N) ^ 2 = κ ^ 2 * (d N : ℝ) := by
      change (κ * Real.sqrt (d N)) ^ 2 = κ ^ 2 * (d N : ℝ)
      rw [mul_pow, Real.sq_sqrt (by positivity)]
    have hfin : (gaussianMatrix (p N) (d N)) (Bad N) ≠ ⊤ := measure_ne_top _ _
    rw [← ENNReal.ofReal_toReal hfin]
    refine ENNReal.ofReal_le_ofReal ?_
    rw [← hsq]
    exact h
  -- the tail bound tends to zero
  have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hdtop
  have hexp0 : Tendsto (fun N => -(κ ^ 2 * (d N : ℝ)) / 2) atTop atBot := by
    have h2 : Tendsto (fun N => (κ ^ 2 / 2) * (d N : ℝ)) atTop atTop :=
      Filter.Tendsto.const_mul_atTop (by positivity) hdR
    have heq : ∀ N : ℕ, -(κ ^ 2 * (d N : ℝ)) / 2 = -((κ ^ 2 / 2) * (d N : ℝ)) := fun N => by
      ring
    simp only [heq]
    exact tendsto_neg_atTop_atBot.comp h2
  have hzero : Tendsto (fun N => (gaussianMatrix (p N) (d N)) (Bad N)) atTop (𝓝 0) := by
    have hE : Tendsto (fun N => ENNReal.ofReal (Real.exp (-(κ ^ 2 * (d N : ℝ)) / 2)))
        atTop (𝓝 0) := by
      have := ENNReal.tendsto_ofReal (Real.tendsto_exp_atBot.comp hexp0)
      simpa using this
    exact tendsto_of_tendsto_of_tendsto_of_le_of_le'
      (tendsto_const_nhds (x := (0 : ℝ≥0∞)) (f := atTop)) hE
      (Eventually.of_forall fun _ => by simp) (Eventually.of_forall hbound)
  -- the eventual inclusion
  have hrsmall : ∀ᶠ N in atTop, Real.sqrt ((p N : ℝ) / (d N : ℝ)) < sc + κ :=
    Filter.Tendsto.eventually_lt_const (by linarith) hcN.sqrt
  have hev : ∀ᶠ N in atTop, (1 : ℝ≥0∞) - (gaussianMatrix (p N) (d N)) (Bad N)
      ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε} := by
    filter_upwards [hrsmall] with N hrN
    have hsub : {ω | B N ω ∈ (Bad N)ᶜ}
        ⊆ {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε} := by
      intro ω hω
      simp only [hBaddef, Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω
      have hGordon := integral_opNorm_le (p := p N) (d := d N) (hp N) (hd N)
      -- `‖B N ω‖ ≤ (r + 1 + κ) √(d N)` with `r = √(p N / d N)`
      set sd : ℝ := Real.sqrt (d N) with hsddef
      set r : ℝ := Real.sqrt ((p N : ℝ) / (d N : ℝ)) with hrdef
      have hsd0 : (0 : ℝ) < sd := hsd N
      have hr0 : (0 : ℝ) ≤ r := Real.sqrt_nonneg _
      have hsp : r * sd = Real.sqrt (p N) := by
        rw [hrdef, hsddef, Real.sqrt_div (by positivity : (0 : ℝ) ≤ (p N : ℝ)),
          div_mul_cancel₀ _ (hsd N).ne']
      have hM0 : (0 : ℝ) ≤ ‖B N ω‖ := norm_nonneg _
      have hMle : ‖B N ω‖ ≤ (r + 1 + κ) * sd := by
        have h1 : ‖B N ω‖ < (∫ Z', ‖Z'‖ ∂(gaussianMatrix (p N) (d N))) + t N := hω
        have h2 : (∫ Z', ‖Z'‖ ∂(gaussianMatrix (p N) (d N))) ≤ Real.sqrt (p N) + sd :=
          hGordon
        have h3 : t N = κ * sd := rfl
        nlinarith [h1, h2, h3, hsp]
      -- the Rayleigh bound
      have hherm : (((d N : ℝ))⁻¹ • ((B N ω)ᵀ * B N ω)).IsHermitian := (hW₀ N ω) ▸ (hsymm N ω)
      have hlm : lamMax (W₀ N ω) (hsymm N ω)
          = lamMax (((d N : ℝ))⁻¹ • ((B N ω)ᵀ * B N ω)) hherm :=
        lamMax_congr (hW₀ N ω) _ _
      have hray : lamMax (((d N : ℝ))⁻¹ • ((B N ω)ᵀ * B N ω)) hherm
          ≤ ((d N : ℝ))⁻¹ * ‖B N ω‖ ^ 2 :=
        lamMax_gram_le_opNorm_sq (by positivity) _ hherm
      -- the arithmetic
      have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hsdsq : sd ^ 2 = (d N : ℝ) := Real.sq_sqrt hdpos.le
      have hrk : (0 : ℝ) ≤ r + 1 + κ := by linarith
      have hMsq : ‖B N ω‖ ^ 2 ≤ ((r + 1 + κ) * sd) ^ 2 := pow_le_pow_left₀ hM0 hMle 2
      have hfinal : ((d N : ℝ))⁻¹ * ‖B N ω‖ ^ 2 ≤ (r + 1 + κ) ^ 2 := by
        rw [inv_mul_le_iff₀ hdpos]
        calc ‖B N ω‖ ^ 2 ≤ ((r + 1 + κ) * sd) ^ 2 := hMsq
          _ = (r + 1 + κ) ^ 2 * (d N : ℝ) := by rw [mul_pow, hsdsq]
          _ = (d N : ℝ) * (r + 1 + κ) ^ 2 := by ring
      have hedge := edge_arith hsc0 hκ0 hκ1 hκK hr0 hrN.le
      simp only [Set.mem_ofPred_eq, bulkEdge, ← hscdef]
      rw [hlm]
      linarith [hray, hfinal, hedge]
    calc (1 : ℝ≥0∞) - (gaussianMatrix (p N) (d N)) (Bad N)
        = (gaussianMatrix (p N) (d N)) ((Bad N)ᶜ) := (hcompl N).symm
      _ = μ N {ω | B N ω ∈ (Bad N)ᶜ} := (hmeasure N).symm
      _ ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε} := measure_mono hsub
  -- squeeze
  have hlow : Tendsto (fun N => (1 : ℝ≥0∞) - (gaussianMatrix (p N) (d N)) (Bad N))
      atTop (𝓝 1) := by
    have := ENNReal.Tendsto.sub (tendsto_const_nhds (x := (1 : ℝ≥0∞)) (f := atTop))
      hzero (Or.inl (by simp))
    simpa using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds hev
    (Eventually.of_forall fun N => prob_le_one)

end R3
end StackedSVD
