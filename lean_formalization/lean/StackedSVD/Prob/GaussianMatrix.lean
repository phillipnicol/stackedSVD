/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.PolynomialNull
import StatsMLlib.Probability.Gaussian.Basic

/-!
# The canonical Gaussian matrix law: flattening, transpose, rotation invariance

Shared infrastructure for items R0, R1, R2 and R3 of `notes/archive/rmt_roadmap.md`.
`gaussianMatrix p d` of `Defs.lean` is a `Measure.pi` of a `Measure.pi`. Every analytic tool
that Layer 2 uses (Stein's identity, Gaussian concentration) lives on an inner product space,
so the first job of this file is one measurable equivalence

  `matrixEquivE p d : Matrix (Fin p) (Fin d) ℝ ≃ᵐ EuclideanSpace ℝ (Fin (p * d))`

that carries `gaussianMatrix p d` to `stdGaussianE (p * d)` and that is an `l2` isometry for
the Frobenius norm (`norm_matrixEquivE_sq`, `dist_matrixEquivE_symm_sq`). The `Matrix` type
carries the sup norm of the pi instance, not the Frobenius norm, so every Lipschitz statement
of items R1, R2 and R3 is stated on `EuclideanSpace ℝ (Fin (p * d))` and moved here.

Contents.

* `matrixUncurry`, `piCongrEquiv`, `matrixEquivFun`, `matrixEquivE`: the flattening chain
  `Matrix (Fin p) (Fin d) ℝ ≃ᵐ (Fin p × Fin d → ℝ) ≃ᵐ (Fin (p * d) → ℝ) ≃ᵐ EuclideanSpace ...`.
* `stdGaussianE_eq`: StatsMLlib's `stdGaussianE n` is Mathlib's `stdGaussian` of
  `EuclideanSpace ℝ (Fin n)`.
* `measurePreserving_matrixEquivE`, `map_matrixEquivE_symm`: the two transport statements.
* `gaussianMatrix_map_transpose`, `measurePreserving_transpose`.
* `gaussianMatrix_map_mul`: invariance under left multiplication by an orthogonal matrix.
* `instIsProbabilityMeasureGaussianMatrix`.

Nothing in this file is asymptotic.
-/

open MeasureTheory ProbabilityTheory
open scoped Matrix RealInnerProductSpace

namespace StackedSVD

variable {p d : ℕ}

/-! ### `gaussianMatrix` is a probability measure -/

instance instIsProbabilityMeasureGaussianMatrix (p d : ℕ) :
    IsProbabilityMeasure (gaussianMatrix p d) := by
  have h : IsProbabilityMeasure
      (Measure.pi fun _ : Fin p => Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
    inferInstance
  exact h

/-! ### The flattening chain -/

/-- Reading a matrix as a function on index pairs, as a measurable equivalence. -/
def matrixUncurry (p d : ℕ) : Matrix (Fin p) (Fin d) ℝ ≃ᵐ (Fin p × Fin d → ℝ) where
  toEquiv := (Equiv.curry (Fin p) (Fin d) ℝ).symm
  measurable_toFun := measurable_pi_lambda _ fun q =>
    (measurable_pi_apply q.2).comp (measurable_pi_apply q.1)
  measurable_invFun :=
    measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => measurable_pi_apply (i, j)

@[simp]
theorem matrixUncurry_apply (Z : Matrix (Fin p) (Fin d) ℝ) (q : Fin p × Fin d) :
    matrixUncurry p d Z q = Z q.1 q.2 := rfl

@[simp]
theorem matrixUncurry_symm_apply (x : Fin p × Fin d → ℝ) (i : Fin p) (j : Fin d) :
    (matrixUncurry p d).symm x i j = x (i, j) := rfl

theorem measurePreserving_matrixUncurry (p d : ℕ) :
    MeasurePreserving (matrixUncurry p d) (gaussianMatrix p d)
      (Measure.pi fun _ : Fin p × Fin d => gaussianReal 0 1) :=
  measurePreserving_uncurry (gaussianReal 0 1)

/-- Reindexing the coordinates of a function space, as a measurable equivalence. -/
def piCongrEquiv {ι κ : Type*} (e : κ ≃ ι) : (ι → ℝ) ≃ᵐ (κ → ℝ) where
  toEquiv :=
    { toFun := fun x => x ∘ e
      invFun := fun y => y ∘ e.symm
      left_inv := fun x => by funext i; simp
      right_inv := fun y => by funext k; simp }
  measurable_toFun := measurable_pi_lambda _ fun k => measurable_pi_apply (e k)
  measurable_invFun := measurable_pi_lambda _ fun i => measurable_pi_apply (e.symm i)

@[simp]
theorem piCongrEquiv_apply {ι κ : Type*} (e : κ ≃ ι) (x : ι → ℝ) (k : κ) :
    piCongrEquiv e x k = x (e k) := rfl

@[simp]
theorem piCongrEquiv_symm_apply {ι κ : Type*} (e : κ ≃ ι) (y : κ → ℝ) (i : ι) :
    (piCongrEquiv e).symm y i = y (e.symm i) := rfl

/-- Reindexing a product of copies of one measure preserves it. -/
theorem measurePreserving_piCongrEquiv {ι κ : Type*} [Fintype ι] [Fintype κ] (e : κ ≃ ι) :
    MeasurePreserving (piCongrEquiv e) (Measure.pi fun _ : ι => gaussianReal 0 1)
      (Measure.pi fun _ : κ => gaussianReal 0 1) := by
  refine ⟨(piCongrEquiv e).measurable, ?_⟩
  refine (Measure.pi_eq fun s hs => ?_).symm
  rw [Measure.map_apply (piCongrEquiv e).measurable (MeasurableSet.univ_pi hs)]
  have hpre : (piCongrEquiv e) ⁻¹' Set.univ.pi s = Set.univ.pi fun i => s (e.symm i) := by
    ext x
    simp only [Set.mem_preimage, Set.mem_univ_pi, piCongrEquiv_apply]
    constructor
    · intro h i
      have hi := h (e.symm i)
      rwa [Equiv.apply_symm_apply] at hi
    · intro h k
      have hk := h (e k)
      rwa [Equiv.symm_apply_apply] at hk
  rw [hpre, Measure.pi_pi]
  exact Equiv.prod_comp e.symm fun k => gaussianReal 0 1 (s k)

/-- The matrix space, flattened to `Fin (p * d) → ℝ`. -/
def matrixEquivFun (p d : ℕ) : Matrix (Fin p) (Fin d) ℝ ≃ᵐ (Fin (p * d) → ℝ) :=
  (matrixUncurry p d).trans (piCongrEquiv finProdFinEquiv.symm)

/-- **The flattening equivalence.** `Matrix (Fin p) (Fin d) ℝ ≃ᵐ EuclideanSpace ℝ (Fin (p * d))`,
an `l2` isometry for the Frobenius norm (`norm_matrixEquivE_sq`). -/
def matrixEquivE (p d : ℕ) : Matrix (Fin p) (Fin d) ℝ ≃ᵐ EuclideanSpace ℝ (Fin (p * d)) :=
  (matrixEquivFun p d).trans (MeasurableEquiv.toLp 2 (Fin (p * d) → ℝ))

@[simp]
theorem matrixEquivE_apply (Z : Matrix (Fin p) (Fin d) ℝ) (r : Fin (p * d)) :
    matrixEquivE p d Z r =
      Z (finProdFinEquiv.symm r).1 (finProdFinEquiv.symm r).2 := rfl

@[simp]
theorem matrixEquivE_symm_apply (x : EuclideanSpace ℝ (Fin (p * d))) (i : Fin p) (j : Fin d) :
    (matrixEquivE p d).symm x i j = x (finProdFinEquiv (i, j)) := rfl

theorem matrixEquivE_symm_sub (x y : EuclideanSpace ℝ (Fin (p * d))) :
    (matrixEquivE p d).symm (x - y) =
      (matrixEquivE p d).symm x - (matrixEquivE p d).symm y := by
  funext i j
  simp

/-! ### `stdGaussianE` is Mathlib's `stdGaussian` -/

/-- StatsMLlib's `stdGaussianE n` is Mathlib's `stdGaussian` on `EuclideanSpace ℝ (Fin n)`. -/
theorem stdGaussianE_eq (n : ℕ) :
    GaussianMeasure.stdGaussianE n = stdGaussian (EuclideanSpace ℝ (Fin n)) := by
  rw [GaussianMeasure.stdGaussianE, GaussianMeasure.stdGaussianPi, ← map_pi_eq_stdGaussian]
  rfl

/-- The standard Gaussian on `EuclideanSpace` is the image of the product measure. -/
theorem measurePreserving_toLp_stdGaussianE (n : ℕ) :
    MeasurePreserving (WithLp.toLp 2 : (Fin n → ℝ) → EuclideanSpace ℝ (Fin n))
      (Measure.pi fun _ : Fin n => gaussianReal 0 1) (GaussianMeasure.stdGaussianE n) :=
  ⟨WithLp.measurable_toLp _ _, by rw [stdGaussianE_eq]; exact map_pi_eq_stdGaussian⟩

/-! ### The two transport statements -/

/-- **`matrixEquivE` carries the Gaussian matrix law to the standard Gaussian.** -/
theorem measurePreserving_matrixEquivE (p d : ℕ) :
    MeasurePreserving (matrixEquivE p d) (gaussianMatrix p d)
      (GaussianMeasure.stdGaussianE (p * d)) :=
  (measurePreserving_toLp_stdGaussianE (p * d)).comp
    ((measurePreserving_piCongrEquiv finProdFinEquiv.symm).comp
      (measurePreserving_matrixUncurry p d))

/-- The inverse direction, in the shape item R3 asks for. -/
theorem measurePreserving_matrixEquivE_symm (p d : ℕ) :
    MeasurePreserving (matrixEquivE p d).symm (GaussianMeasure.stdGaussianE (p * d))
      (gaussianMatrix p d) :=
  (measurePreserving_matrixEquivE p d).symm (matrixEquivE p d)

theorem map_matrixEquivE_symm (p d : ℕ) :
    (GaussianMeasure.stdGaussianE (p * d)).map (matrixEquivE p d).symm = gaussianMatrix p d :=
  (measurePreserving_matrixEquivE_symm p d).map_eq

/-! ### `matrixEquivE` is a Frobenius isometry -/

/-- The `l2` norm of the flattened matrix is the Frobenius norm. -/
theorem norm_matrixEquivE_sq (Z : Matrix (Fin p) (Fin d) ℝ) :
    ‖matrixEquivE p d Z‖ ^ 2 = ∑ i, ∑ j, Z i j ^ 2 := by
  rw [EuclideanSpace.real_norm_sq_eq,
    ← Equiv.sum_comp finProdFinEquiv fun r => matrixEquivE p d Z r ^ 2, Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  simp

theorem norm_matrixEquivE (Z : Matrix (Fin p) (Fin d) ℝ) :
    ‖matrixEquivE p d Z‖ = Real.sqrt (∑ i, ∑ j, Z i j ^ 2) := by
  rw [← norm_matrixEquivE_sq, Real.sqrt_sq (norm_nonneg _)]

/-- The distance of two points of `EuclideanSpace ℝ (Fin (p * d))` is the Frobenius distance of
the two matrices they name. Items R1, R2 and R3 use this to read a bound on the entries as a
`LipschitzWith` statement. -/
theorem dist_matrixEquivE_symm_sq (x y : EuclideanSpace ℝ (Fin (p * d))) :
    dist x y ^ 2 = ∑ i, ∑ j,
      ((matrixEquivE p d).symm x i j - (matrixEquivE p d).symm y i j) ^ 2 := by
  have h : x - y = matrixEquivE p d ((matrixEquivE p d).symm x - (matrixEquivE p d).symm y) := by
    rw [← matrixEquivE_symm_sub]
    exact ((matrixEquivE p d).apply_symm_apply (x - y)).symm
  rw [dist_eq_norm, h, norm_matrixEquivE_sq]
  rfl

/-! ### Transpose -/

/-- **The transpose of a canonical Gaussian matrix is a canonical Gaussian matrix.** Not in
Mathlib v4.33.0 (`Equiv.piComm` carries no measure theory). -/
theorem measurePreserving_transpose (p d : ℕ) :
    MeasurePreserving (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin d) (Fin p) ℝ)
      (gaussianMatrix p d) (gaussianMatrix d p) := by
  have h1 := measurePreserving_matrixUncurry p d
  have h2 := measurePreserving_piCongrEquiv (Equiv.prodComm (Fin d) (Fin p))
  have h3 := (measurePreserving_matrixUncurry d p).symm (matrixUncurry d p)
  have hcomp := h3.comp (h2.comp h1)
  have hfun : ((matrixUncurry d p).symm ∘ (piCongrEquiv (Equiv.prodComm (Fin d) (Fin p))) ∘
      (matrixUncurry p d)) = (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → _) := rfl
  rwa [hfun] at hcomp

theorem gaussianMatrix_map_transpose (p d : ℕ) :
    (gaussianMatrix p d).map
        (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin d) (Fin p) ℝ)
      = gaussianMatrix d p :=
  (measurePreserving_transpose p d).map_eq

/-! ### Rotation invariance -/

section Rotation

variable {U : Matrix (Fin p) (Fin p) ℝ}

/-- The real inner product of `EuclideanSpace` is the dot product. -/
theorem inner_euclidean_eq_dotProduct (a b : EuclideanSpace ℝ (Fin p)) :
    ⟪a, b⟫ = WithLp.ofLp a ⬝ᵥ WithLp.ofLp b := by
  simp [PiLp.inner_apply, dotProduct, mul_comm]

/-- An orthogonal matrix preserves the dot product. -/
theorem dotProduct_mulVec_mulVec (hU : Uᵀ * U = 1) (x y : Fin p → ℝ) :
    (U *ᵥ x) ⬝ᵥ (U *ᵥ y) = x ⬝ᵥ y := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec, hU,
    Matrix.one_mulVec]

theorem inner_toEuclideanLin (hU : Uᵀ * U = 1) (x y : EuclideanSpace ℝ (Fin p)) :
    ⟪Matrix.toEuclideanLin U x, Matrix.toEuclideanLin U y⟫ = ⟪x, y⟫ := by
  rw [inner_euclidean_eq_dotProduct, inner_euclidean_eq_dotProduct]
  change (U *ᵥ WithLp.ofLp x) ⬝ᵥ (U *ᵥ WithLp.ofLp y) = _
  rw [dotProduct_mulVec_mulVec hU]

/-- A real matrix with `Uᵀ U = 1` acts on `EuclideanSpace ℝ (Fin p)` as a linear isometry
equivalence. -/
noncomputable def rotIso (U : Matrix (Fin p) (Fin p) ℝ) (hU : Uᵀ * U = 1) :
    EuclideanSpace ℝ (Fin p) ≃ₗᵢ[ℝ] EuclideanSpace ℝ (Fin p) :=
  LinearIsometry.toLinearIsometryEquiv
    ((Matrix.toEuclideanLin U).isometryOfInner (inner_toEuclideanLin hU)) rfl

@[simp]
theorem rotIso_apply (hU : Uᵀ * U = 1) (x : EuclideanSpace ℝ (Fin p)) :
    rotIso U hU x = WithLp.toLp 2 (U *ᵥ WithLp.ofLp x) := rfl

/-- **Rotation invariance for one Gaussian vector.** -/
theorem measurePreserving_mulVec (hU : Uᵀ * U = 1) :
    MeasurePreserving (fun y : Fin p → ℝ => U *ᵥ y)
      (Measure.pi fun _ : Fin p => gaussianReal 0 1)
      (Measure.pi fun _ : Fin p => gaussianReal 0 1) := by
  set g : Measure (Fin p → ℝ) := Measure.pi fun _ : Fin p => gaussianReal 0 1 with hg
  have hmeas : Measurable (fun y : Fin p → ℝ => U *ᵥ y) :=
    (Matrix.mulVecLin U).continuous_of_finiteDimensional.measurable
  refine ⟨hmeas, ?_⟩
  set e : (Fin p → ℝ) ≃ᵐ EuclideanSpace ℝ (Fin p) := MeasurableEquiv.toLp 2 (Fin p → ℝ) with he
  have h1 : g.map e = stdGaussian (EuclideanSpace ℝ (Fin p)) := map_pi_eq_stdGaussian
  have hcomp : (⇑e ∘ fun y : Fin p → ℝ => U *ᵥ y) = (rotIso U hU) ∘ ⇑e := rfl
  have h3 : (g.map (fun y : Fin p → ℝ => U *ᵥ y)).map e = g.map e := by
    rw [Measure.map_map e.measurable hmeas, hcomp,
      ← Measure.map_map (rotIso U hU).continuous.measurable e.measurable, h1,
      stdGaussian_map (rotIso U hU)]
  calc g.map (fun y : Fin p → ℝ => U *ᵥ y)
      = ((g.map (fun y : Fin p → ℝ => U *ᵥ y)).map e).map e.symm :=
        (MeasurableEquiv.map_symm_map e).symm
    _ = (g.map e).map e.symm := by rw [h3]
    _ = g := MeasurableEquiv.map_symm_map e

/-- Row by row action on a Gaussian matrix. -/
theorem measurePreserving_rowMulVec (hU : Uᵀ * U = 1) (d : ℕ) :
    MeasurePreserving (fun (Y : Matrix (Fin d) (Fin p) ℝ) (i : Fin d) => U *ᵥ Y i)
      (gaussianMatrix d p) (gaussianMatrix d p) :=
  measurePreserving_pi _ _ fun _ => measurePreserving_mulVec hU

/-- **Rotation invariance.** Left multiplication by an orthogonal matrix preserves the
canonical Gaussian matrix law. -/
theorem measurePreserving_mul_left (hU : Uᵀ * U = 1) (d : ℕ) :
    MeasurePreserving (fun Z : Matrix (Fin p) (Fin d) ℝ => U * Z)
      (gaussianMatrix p d) (gaussianMatrix p d) := by
  have hcomp := (measurePreserving_transpose d p).comp
    ((measurePreserving_rowMulVec hU d).comp (measurePreserving_transpose p d))
  have hfun : (fun Z : Matrix (Fin p) (Fin d) ℝ => U * Z) =
      (Matrix.transpose ∘ (fun (Y : Matrix (Fin d) (Fin p) ℝ) (i : Fin d) => U *ᵥ Y i) ∘
        (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin d) (Fin p) ℝ)) := by
    funext Z
    ext k i
    change (U * Z) k i = (U *ᵥ Zᵀ i) k
    simp [Matrix.mul_apply, Matrix.mulVec, dotProduct]
  rw [hfun]
  exact hcomp

theorem gaussianMatrix_map_mul (hU : Uᵀ * U = 1) (d : ℕ) :
    (gaussianMatrix p d).map (fun Z : Matrix (Fin p) (Fin d) ℝ => U * Z) = gaussianMatrix p d :=
  (measurePreserving_mul_left hU d).map_eq

end Rotation

end StackedSVD
