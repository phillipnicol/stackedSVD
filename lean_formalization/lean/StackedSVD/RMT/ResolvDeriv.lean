/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R4C
import StackedSVD.Prob.GaussianMatrix

/-!
# The derivative of the complex resolvent in the entries of the Gaussian block

Item R1 (`notes/archive/rmt_R1.md`, steps 2 and 3) needs one analytic fact that `RMT/R4C.lean` does
not carry: the map `Y ↦ (cmat (d⁻¹ Yᵀ Y) - z)⁻¹` is `C^∞` in the entries of `Y`, with the
derivative `-G (dW) G`, and the trace `s(Y) = d⁻¹ tr G` has the global gradient bound
`‖∇ Re s‖_F ≤ 2 √(Im z + ‖z‖) / (d (Im z)²)`.

`Matrix` carries **no** norm instance in Mathlib v4.33.0 (every candidate in
`Mathlib/Analysis/Matrix/Normed.lean` is a scoped non-instance), and the scoped ones clash
with `instTopologicalSpaceMatrix`, so `E →L[ℝ] Matrix ...` is not a normed space. The
calculus therefore happens in the operator algebra `Op d = (Fin d → ℂ) →L[ℝ] (Fin d → ℂ)`,
which is a complete normed `ℝ`-algebra, and `toOpₗ` carries a matrix there. Matrices come
back through the `ℝ`-linear functionals `entryOp` and `traceOp`.

Everything here is deterministic. Nothing is asymptotic.
-/

open Filter Topology
open scoped Matrix NNReal RealInnerProductSpace

namespace StackedSVD

namespace ResolvDeriv

open R4

variable {p d : ℕ} {z : ℂ} {W : Matrix (Fin d) (Fin d) ℝ}

/-! ### The Wishart block and the complexification -/

/-- `W₀(Y) = d⁻¹ Yᵀ Y`. The same term as `R1.W0`. -/
noncomputable def gram (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  (d : ℝ)⁻¹ • (Yᵀ * Y)

theorem isHermitian_gram (Y : Matrix (Fin p) (Fin d) ℝ) : (gram Y).IsHermitian := by
  have h := isHermitian_transpose_mul_self Y
  change ((d : ℝ)⁻¹ • (Yᵀ * Y))ᴴ = (d : ℝ)⁻¹ • (Yᵀ * Y)
  rw [Matrix.conjTranspose_smul, star_trivial, h.eq]

/-- A real matrix read over `ℂ`; the rectangular form of `R4C.cmat`. -/
noncomputable def cmapR {m n : Type*} (M : Matrix m n ℝ) : Matrix m n ℂ :=
  M.map (fun a => (a : ℂ))

@[simp] theorem cmapR_apply {m n : Type*} (M : Matrix m n ℝ) (i : m) (j : n) :
    cmapR M i j = (M i j : ℂ) := rfl

theorem cmat_eq_cmapR (W : Matrix (Fin d) (Fin d) ℝ) : R4C.cmat W = cmapR W := rfl

theorem cmapR_mul {l m n : Type*} [Fintype m] (A : Matrix l m ℝ) (B : Matrix m n ℝ) :
    cmapR (A * B) = cmapR A * cmapR B := by
  ext i j
  simp only [cmapR_apply, Matrix.mul_apply]
  push_cast
  rfl

theorem cmapR_transpose {m n : Type*} (A : Matrix m n ℝ) : cmapR Aᵀ = (cmapR A)ᵀ := rfl

/-! ### Matrices as operators on `Fin d → ℂ` -/

/-- A complex matrix as a real-linear continuous operator on `Fin d → ℂ`. -/
noncomputable def toOpₗ (d : ℕ) :
    Matrix (Fin d) (Fin d) ℂ →ₗ[ℝ] ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) where
  toFun M := LinearMap.toContinuousLinearMap ((M.mulVecLin).restrictScalars ℝ)
  map_add' := by
    intro M N
    ext v i
    simp
  map_smul' := by
    intro c M
    ext v i
    simp

@[simp] theorem toOpₗ_apply (M : Matrix (Fin d) (Fin d) ℂ) (v : Fin d → ℂ) :
    toOpₗ d M v = M *ᵥ v := rfl

theorem toOpₗ_mul (M N : Matrix (Fin d) (Fin d) ℂ) :
    toOpₗ d (M * N) = toOpₗ d M * toOpₗ d N := by
  ext v i
  change ((M * N) *ᵥ v) i = (M *ᵥ (N *ᵥ v)) i
  rw [Matrix.mulVec_mulVec]

theorem toOpₗ_one : toOpₗ d (1 : Matrix (Fin d) (Fin d) ℂ) = 1 := by
  ext v i
  change ((1 : Matrix (Fin d) (Fin d) ℂ) *ᵥ v) i = v i
  rw [Matrix.one_mulVec]

/-- One entry of an operator, as a real-linear functional. -/
noncomputable def entryOpₗ (d : ℕ) (i j : Fin d) :
    ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) →ₗ[ℝ] ℂ where
  toFun T := T (Pi.single j 1) i
  map_add' := by intro T S; simp
  map_smul' := by intro c T; simp

theorem entryOpₗ_toOpₗ (M : Matrix (Fin d) (Fin d) ℂ) (i j : Fin d) :
    entryOpₗ d i j (toOpₗ d M) = M i j := by
  change (M *ᵥ Pi.single j (1 : ℂ)) i = M i j
  simp

/-- The trace of an operator, as a real-linear functional. -/
noncomputable def traceOpₗ (d : ℕ) : ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) →ₗ[ℝ] ℂ :=
  ∑ i, entryOpₗ d i i

theorem traceOpₗ_toOpₗ (M : Matrix (Fin d) (Fin d) ℂ) :
    traceOpₗ d (toOpₗ d M) = M.trace := by
  rw [traceOpₗ, LinearMap.sum_apply, Matrix.trace]
  exact Finset.sum_congr rfl fun i _ => entryOpₗ_toOpₗ M i i

/-! ### Bilinear maps in finite dimension -/

/-- A bilinear map between finite-dimensional real normed spaces, as a continuous bilinear
map. -/
noncomputable def bilCLM {E F G : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] [NormedAddCommGroup F] [NormedSpace ℝ F] [FiniteDimensional ℝ F]
    [NormedAddCommGroup G] [NormedSpace ℝ G] (B : E →ₗ[ℝ] F →ₗ[ℝ] G) : E →L[ℝ] F →L[ℝ] G :=
  LinearMap.toContinuousLinearMap
    ((LinearMap.toContinuousLinearMap :
        (F →ₗ[ℝ] G) ≃ₗ[ℝ] (F →L[ℝ] G)).toLinearMap.comp B)

@[simp] theorem bilCLM_apply {E F G : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] [NormedAddCommGroup F] [NormedSpace ℝ F] [FiniteDimensional ℝ F]
    [NormedAddCommGroup G] [NormedSpace ℝ G] (B : E →ₗ[ℝ] F →ₗ[ℝ] G) (u : E) (v : F) :
    bilCLM B u v = B u v := rfl

/-! ### The Gram bilinear map on the flattened space -/

theorem cmapR_symm_add (u v : EuclideanSpace ℝ (Fin (p * d))) :
    cmapR ((matrixEquivE p d).symm (u + v))
      = cmapR ((matrixEquivE p d).symm u) + cmapR ((matrixEquivE p d).symm v) := by
  ext i j
  simp [cmapR]

theorem cmapR_symm_smul (c : ℝ) (u : EuclideanSpace ℝ (Fin (p * d))) :
    cmapR ((matrixEquivE p d).symm (c • u)) = c • cmapR ((matrixEquivE p d).symm u) := by
  ext i j
  simp [cmapR, Complex.real_smul]

/-- `(u, v) ↦ d⁻¹ Uᵀ V` over `ℂ`, on the flattened space, as a bilinear map. -/
noncomputable def gramBilₗ (p d : ℕ) :
    EuclideanSpace ℝ (Fin (p * d)) →ₗ[ℝ] EuclideanSpace ℝ (Fin (p * d)) →ₗ[ℝ]
      Matrix (Fin d) (Fin d) ℂ :=
  LinearMap.mk₂ ℝ
    (fun u v => (d : ℂ)⁻¹ •
      ((cmapR ((matrixEquivE p d).symm u))ᵀ * cmapR ((matrixEquivE p d).symm v)))
    (by
      intro u₁ u₂ v
      rw [cmapR_symm_add, Matrix.transpose_add, Matrix.add_mul, smul_add])
    (by
      intro c u v
      rw [cmapR_symm_smul, Matrix.transpose_smul, Matrix.smul_mul, smul_comm])
    (by
      intro u v₁ v₂
      rw [cmapR_symm_add, Matrix.mul_add, smul_add])
    (by
      intro c u v
      rw [cmapR_symm_smul, Matrix.mul_smul, smul_comm])

theorem gramBilₗ_apply (u v : EuclideanSpace ℝ (Fin (p * d))) :
    gramBilₗ p d u v = (d : ℂ)⁻¹ •
      ((cmapR ((matrixEquivE p d).symm u))ᵀ * cmapR ((matrixEquivE p d).symm v)) := rfl

theorem gramBilₗ_self (x : EuclideanSpace ℝ (Fin (p * d))) :
    gramBilₗ p d x x = R4C.cmat (gram ((matrixEquivE p d).symm x)) := by
  ext a b
  rw [gramBilₗ_apply]
  simp only [Matrix.smul_apply, Matrix.mul_apply, Matrix.transpose_apply, cmapR_apply,
    smul_eq_mul, cmat_eq_cmapR, cmapR_apply, gram, Matrix.smul_apply, smul_eq_mul]
  push_cast
  rfl

/-- The Gram map into the operator algebra, as a continuous bilinear map. -/
noncomputable def gramOpₗ (p d : ℕ) :
    EuclideanSpace ℝ (Fin (p * d)) →ₗ[ℝ] EuclideanSpace ℝ (Fin (p * d)) →ₗ[ℝ]
      ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) :=
  LinearMap.mk₂ ℝ (fun u v => toOpₗ d (gramBilₗ p d u v))
    (by intro u₁ u₂ v; rw [map_add, LinearMap.add_apply, map_add])
    (by intro c u v; rw [map_smul, LinearMap.smul_apply, map_smul])
    (by intro u v₁ v₂; rw [map_add, map_add])
    (by intro c u v; rw [map_smul, map_smul])

noncomputable def gramOp (p d : ℕ) := bilCLM (gramOpₗ p d)

theorem gramOp_apply (u v : EuclideanSpace ℝ (Fin (p * d))) :
    gramOp p d u v = toOpₗ d (gramBilₗ p d u v) := rfl

/-! ### The resolvent is a two-sided inverse -/

theorem cmat_sub_mul_resolvC (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) * R4C.resolvC W z = 1 := by
  have hne := R4C.eigenvalue_sub_ne_zero hW hz
  have hDD : Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
      Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [Matrix.diagonal_mul_diagonal,
      show (fun a => ((hW.eigenvalues a : ℂ) - z) * ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = fun _ => (1 : ℂ) from funext fun a => mul_inv_cancel₀ (hne a)]
    exact Matrix.diagonal_one
  rw [R4C.resolvC_eq_conj hW hz, R4C.cmat_sub_smul_one_eq_conj hW]
  have hassoc : R4C.cmat (eigU hW) * Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
      (R4C.cmat (eigU hW))ᵀ *
      (R4C.cmat (eigU hW) * Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
        (R4C.cmat (eigU hW))ᵀ)
      = R4C.cmat (eigU hW) * (Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
          ((R4C.cmat (eigU hW))ᵀ * R4C.cmat (eigU hW)) *
          Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)) *
        (R4C.cmat (eigU hW))ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, R4C.transpose_ceigU_mul, Matrix.mul_one, hDD, Matrix.mul_one,
    R4C.ceigU_mul_transpose]

theorem isUnit_det_cmat_sub (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    IsUnit (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)).det :=
  Matrix.isUnit_det_of_right_inverse (cmat_sub_mul_resolvC hW hz)

theorem resolvC_mul_cmat_sub (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    R4C.resolvC W z * (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) = 1 :=
  Matrix.nonsing_inv_mul _ (isUnit_det_cmat_sub hW hz)

/-- The complex resolvent of a real symmetric matrix is symmetric. -/
theorem transpose_resolvC (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    (R4C.resolvC W z)ᵀ = R4C.resolvC W z := by
  rw [R4C.resolvC_eq_conj hW hz]
  simp only [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.diagonal_transpose,
    Matrix.mul_assoc]

/-! ### `A(x)` and `G(x)` in the operator algebra -/

/-- `A(x) = cmat (W₀ Y) - z`. -/
noncomputable def Amat (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin d) (Fin d) ℂ :=
  R4C.cmat (gram ((matrixEquivE p d).symm x)) - z • 1

/-- `G(x) = (cmat (W₀ Y) - z)⁻¹`. -/
noncomputable def Gmat (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin d) (Fin d) ℂ :=
  R4C.resolvC (gram ((matrixEquivE p d).symm x)) z

theorem Amat_mul_Gmat (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Amat z x * Gmat z x = 1 :=
  cmat_sub_mul_resolvC (isHermitian_gram _) hz

theorem Gmat_mul_Amat (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Gmat z x * Amat z x = 1 :=
  resolvC_mul_cmat_sub (isHermitian_gram _) hz

noncomputable def Aop (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Fin d → ℂ) →L[ℝ] (Fin d → ℂ) :=
  gramOp p d x x - toOpₗ d (z • 1)

noncomputable def Gop (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Fin d → ℂ) →L[ℝ] (Fin d → ℂ) :=
  toOpₗ d (Gmat z x)

theorem Aop_eq (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Aop z x = toOpₗ d (Amat z x) := by
  rw [Aop, Amat, gramOp_apply, gramBilₗ_self, map_sub]

/-- `A(x)` as a unit of the operator algebra. -/
noncomputable def AopUnit (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ))ˣ where
  val := Aop z x
  inv := Gop z x
  val_inv := by rw [Aop_eq, Gop, ← toOpₗ_mul, Amat_mul_Gmat hz, toOpₗ_one]
  inv_val := by rw [Aop_eq, Gop, ← toOpₗ_mul, Gmat_mul_Amat hz, toOpₗ_one]

private theorem coe_AopUnit_inv (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (((AopUnit hz x)⁻¹ : ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ))ˣ) :
        (Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) = Gop z x := rfl

theorem ringInverse_Aop (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Ring.inverse (Aop z x) = Gop z x :=
  Ring.inverse_unit (AopUnit hz x)

theorem Gop_eq_comp (hz : z.im ≠ 0) :
    (fun y => Ring.inverse (Aop (p := p) (d := d) z y)) = Gop z :=
  funext fun y => ringInverse_Aop hz y

/-! ### The derivative -/

/-- The derivative of `x ↦ d⁻¹ Yᵀ Y` at `x`, as a matrix: `H ↦ d⁻¹ (Yᵀ H + Hᵀ Y)`. -/
noncomputable def dWmat (p d : ℕ) (x h : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin d) (Fin d) ℂ :=
  gramBilₗ p d x h + gramBilₗ p d h x

noncomputable def dAop (p d : ℕ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    EuclideanSpace ℝ (Fin (p * d)) →L[ℝ] ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) :=
  ContinuousLinearMap.comp (gramOp p d x)
      (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin (p * d))))
    + (gramOp p d).flip x

theorem dAop_apply (x h : EuclideanSpace ℝ (Fin (p * d))) :
    dAop p d x h = toOpₗ d (dWmat p d x h) := by
  change gramOp p d x h + gramOp p d h x = _
  rw [gramOp_apply, gramOp_apply, dWmat, map_add]

theorem hasFDerivAt_Aop (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    HasFDerivAt (Aop (p := p) (d := d) z) (dAop p d x) x := by
  have h0 : HasFDerivAt (fun y : EuclideanSpace ℝ (Fin (p * d)) => gramOp p d y y)
      (dAop p d x) x := ((gramOp p d).hasFDerivAt (x := x)).clm_apply (hasFDerivAt_id x)
  exact h0.sub_const _

theorem contDiff_Aop (z : ℂ) {n : WithTop ℕ∞} :
    ContDiff ℝ n (Aop (p := p) (d := d) z) := by
  have hc : ContDiff ℝ n (⇑(gramOp p d)) := (gramOp p d).contDiff
  have h1 : ContDiff ℝ n (fun y : EuclideanSpace ℝ (Fin (p * d)) => gramOp p d y y) :=
    hc.clm_apply contDiff_id
  exact h1.sub contDiff_const

theorem contDiff_Gop (hz : z.im ≠ 0) {n : WithTop ℕ∞} :
    ContDiff ℝ n (Gop (p := p) (d := d) z) := by
  rw [contDiff_iff_contDiffAt]
  intro x
  have h1 : ContDiffAt ℝ n Ring.inverse (Aop z x) := contDiffAt_ringInverse ℝ (AopUnit hz x)
  have h2 := h1.comp x (contDiff_Aop (p := p) (d := d) (n := n) z).contDiffAt
  rw [Function.comp_def, Gop_eq_comp hz] at h2
  exact h2

/-! ### The trace and its derivative -/

/-- The trace of an operator, as a continuous real-linear functional. -/
noncomputable def traceOp (d : ℕ) : ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) →L[ℝ] ℂ :=
  LinearMap.toContinuousLinearMap (traceOpₗ d)

@[simp] theorem traceOp_apply (T : (Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) :
    traceOp d T = traceOpₗ d T := rfl

/-- `s(x) = d⁻¹ tr G(x)`, on the flattened space. -/
noncomputable def sfun (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) : ℂ :=
  (d : ℂ)⁻¹ * traceOp d (Gop z x)

theorem sfun_eq (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    sfun z x = R4C.stieltjesC (gram ((matrixEquivE p d).symm x)) z := by
  rw [sfun, traceOp_apply, Gop, traceOpₗ_toOpₗ]
  rfl

/-! ### The derivative of the trace

The derivative of `Ring.inverse` carries `ContinuousLinearMap.mulLeftRight ℝ R`, whose
instance arguments do not match the ones Lean finds for a hand-written
`(Fin d → ℂ) →L[ℝ] (Fin d → ℂ)`. Every statement below therefore names no continuous linear
map of its own: Mathlib produces the derivative and only its values are stated. -/

private theorem sfun_deriv_aux (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (sfun (p := p) (d := d) z) x ∧
      ∀ h, fderiv ℝ (sfun (p := p) (d := d) z) x h
        = -((d : ℂ)⁻¹ * (Gmat z x * dWmat p d x h * Gmat z x).trace) := by
  have hinv := hasFDerivAt_ringInverse (𝕜 := ℝ) (AopUnit hz x)
  have hG0 := hinv.comp x (hasFDerivAt_Aop (p := p) (d := d) z x)
  rw [Function.comp_def, Gop_eq_comp hz] at hG0
  have hT := ((traceOp d).hasFDerivAt).comp x hG0
  rw [Function.comp_def] at hT
  have hS : HasFDerivAt (sfun (p := p) (d := d) z) _ x := hT.const_mul ((d : ℂ)⁻¹)
  refine ⟨hS.differentiableAt, fun h => ?_⟩
  rw [hS.fderiv]
  simp only [smul_apply, ContinuousLinearMap.coe_comp, Function.comp_apply, neg_apply,
    ContinuousLinearMap.mulLeftRight_apply, smul_eq_mul, traceOp_apply, dAop_apply]
  rw [coe_AopUnit_inv, Gop, ← toOpₗ_mul, ← toOpₗ_mul, map_neg, traceOpₗ_toOpₗ]
  ring

theorem differentiableAt_sfun (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (sfun (p := p) (d := d) z) x := (sfun_deriv_aux hz x).1

theorem differentiable_sfun (hz : z.im ≠ 0) :
    Differentiable ℝ (sfun (p := p) (d := d) z) := fun x => differentiableAt_sfun hz x

/-- The gradient matrix `Y G²` of step 3 of `notes/archive/rmt_R1.md`. -/
noncomputable def gradMat (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    Matrix (Fin p) (Fin d) ℂ :=
  cmapR ((matrixEquivE p d).symm x) * Gmat z x * Gmat z x

theorem transpose_Gmat (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Gmat z x)ᵀ = Gmat z x := transpose_resolvC (isHermitian_gram _) hz

theorem transpose_Gmat_sq (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    (Gmat z x * Gmat z x)ᵀ = Gmat z x * Gmat z x := by
  rw [Matrix.transpose_mul, transpose_Gmat hz]

/-- Both halves of `∂W` contribute the same trace. -/
private theorem trace_swap (hz : z.im ≠ 0) (x h : EuclideanSpace ℝ (Fin (p * d))) :
    (Gmat z x * Gmat z x * (cmapR ((matrixEquivE p d).symm x))ᵀ *
        cmapR ((matrixEquivE p d).symm h)).trace
      = (Gmat z x * Gmat z x * (cmapR ((matrixEquivE p d).symm h))ᵀ *
        cmapR ((matrixEquivE p d).symm x)).trace := by
  set S := Gmat z x * Gmat z x with hS
  set Yc := cmapR ((matrixEquivE p d).symm x) with hYc
  set Hc := cmapR ((matrixEquivE p d).symm h) with hHc
  have hSt : Sᵀ = S := transpose_Gmat_sq hz x
  calc (S * Ycᵀ * Hc).trace = ((S * Ycᵀ * Hc)ᵀ).trace := (Matrix.trace_transpose _).symm
    _ = (Hcᵀ * Yc * S).trace := by
        rw [Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_transpose, hSt,
          ← Matrix.mul_assoc]
    _ = (S * Hcᵀ * Yc).trace := Matrix.trace_mul_cycle _ _ _

private theorem trace_eq_double_sum (_hz : z.im ≠ 0) (x h : EuclideanSpace ℝ (Fin (p * d))) :
    (Gmat z x * Gmat z x * (cmapR ((matrixEquivE p d).symm h))ᵀ *
        cmapR ((matrixEquivE p d).symm x)).trace
      = ∑ k : Fin p, ∑ j : Fin d, gradMat z x k j * (h (finProdFinEquiv (k, j)) : ℂ) := by
  set S := Gmat z x * Gmat z x with hS
  set Yc := cmapR ((matrixEquivE p d).symm x) with hYc
  set Hc := cmapR ((matrixEquivE p d).symm h) with hHc
  have h1 : (S * Hcᵀ * Yc).trace = (Yc * S * Hcᵀ).trace := Matrix.trace_mul_cycle _ _ _
  rw [h1, Matrix.trace]
  refine Finset.sum_congr rfl fun k _ => ?_
  simp only [Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hg : gradMat z x k j = (Yc * S) k j := by
    rw [gradMat, hYc, hS, Matrix.mul_assoc]
  rw [hg, Matrix.mul_apply, hHc]
  rfl

theorem fderiv_sfun_apply (hz : z.im ≠ 0) (x h : EuclideanSpace ℝ (Fin (p * d))) :
    fderiv ℝ (sfun (p := p) (d := d) z) x h
      = -(2 / (d : ℂ) ^ 2) *
        ∑ k : Fin p, ∑ j : Fin d, gradMat z x k j * (h (finProdFinEquiv (k, j)) : ℂ) := by
  rw [(sfun_deriv_aux hz x).2 h]
  have hcyc : (Gmat z x * dWmat p d x h * Gmat z x).trace
      = (Gmat z x * Gmat z x * dWmat p d x h).trace := Matrix.trace_mul_cycle _ _ _
  have hdW : dWmat p d x h
      = (d : ℂ)⁻¹ • ((cmapR ((matrixEquivE p d).symm x))ᵀ * cmapR ((matrixEquivE p d).symm h))
        + (d : ℂ)⁻¹ • ((cmapR ((matrixEquivE p d).symm h))ᵀ *
            cmapR ((matrixEquivE p d).symm x)) := rfl
  rw [hcyc, hdW, Matrix.mul_add, Matrix.trace_add, Matrix.mul_smul, Matrix.mul_smul,
    Matrix.trace_smul, Matrix.trace_smul, ← Matrix.mul_assoc, ← Matrix.mul_assoc,
    trace_swap hz x h, trace_eq_double_sum hz x h]
  simp only [smul_eq_mul]
  ring

/-! ### The Frobenius norm of `Y G²`

`‖Y G²‖_F² = d ∑_a λ_a |g_a|⁴ ≤ d² (η + ‖z‖) / η⁴`, in the eigenbasis of `W₀(Y)`. -/

/-- The eigenbasis of `W₀(Y)`. -/
noncomputable def eigUx (x : EuclideanSpace ℝ (Fin (p * d))) : Matrix (Fin d) (Fin d) ℝ :=
  eigU (isHermitian_gram ((matrixEquivE p d).symm x))

/-- The eigenvalues of `W₀(Y)`. -/
noncomputable def eigVal (x : EuclideanSpace ℝ (Fin (p * d))) (a : Fin d) : ℝ :=
  (isHermitian_gram ((matrixEquivE p d).symm x)).eigenvalues a

/-- `Z = Y U`, the block in the eigenbasis. -/
noncomputable def Zmat (x : EuclideanSpace ℝ (Fin (p * d))) : Matrix (Fin p) (Fin d) ℝ :=
  (matrixEquivE p d).symm x * eigUx x

/-- `g_a = (λ_a - z)⁻¹`. -/
noncomputable def gval (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) (a : Fin d) : ℂ :=
  ((eigVal x a : ℂ) - z)⁻¹

theorem transpose_eigUx_mul (x : EuclideanSpace ℝ (Fin (p * d))) :
    (eigUx x)ᵀ * eigUx x = 1 := transpose_eigU_mul _

/-- Orthogonal invariance of the `l²` norm of a complex vector. -/
private theorem sum_normSq_mulVec {V : Matrix (Fin d) (Fin d) ℝ} (hV : Vᵀ * V = 1)
    (w : Fin d → ℂ) :
    ∑ j : Fin d, ‖∑ a : Fin d, w a * ((V j a : ℝ) : ℂ)‖ ^ 2 = ∑ a : Fin d, ‖w a‖ ^ 2 := by
  set u : Fin d → ℝ := fun a => (w a).re with hu
  set v : Fin d → ℝ := fun a => (w a).im with hv
  have hre : ∀ j, (∑ a : Fin d, w a * ((V j a : ℝ) : ℂ)).re = (V *ᵥ u) j := by
    intro j
    rw [Complex.re_sum, Matrix.mulVec, dotProduct]
    exact Finset.sum_congr rfl fun a _ => by simp [hu, Complex.mul_re, mul_comm]
  have him : ∀ j, (∑ a : Fin d, w a * ((V j a : ℝ) : ℂ)).im = (V *ᵥ v) j := by
    intro j
    rw [Complex.im_sum, Matrix.mulVec, dotProduct]
    exact Finset.sum_congr rfl fun a _ => by simp [hv, Complex.mul_im, mul_comm]
  have hnorm : ∀ j, ‖∑ a : Fin d, w a * ((V j a : ℝ) : ℂ)‖ ^ 2
      = ((V *ᵥ u) j) ^ 2 + ((V *ᵥ v) j) ^ 2 := by
    intro j
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
    ring
  have hsq : ∀ (t : Fin d → ℝ), ∑ j, ((V *ᵥ t) j) ^ 2 = ∑ a, (t a) ^ 2 := by
    intro t
    have hd := dotProduct_mulVec_mulVec hV t t
    simp only [dotProduct] at hd
    calc ∑ j, ((V *ᵥ t) j) ^ 2 = ∑ j, (V *ᵥ t) j * (V *ᵥ t) j :=
          Finset.sum_congr rfl fun j _ => sq _
      _ = ∑ a, t a * t a := hd
      _ = ∑ a, (t a) ^ 2 := Finset.sum_congr rfl fun a _ => (sq _).symm
  rw [Finset.sum_congr rfl (fun j (_ : j ∈ Finset.univ) => hnorm j), Finset.sum_add_distrib,
    hsq u, hsq v, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun a _ => by
    rw [Complex.sq_norm, Complex.normSq_apply, hu, hv]; ring

/-- `Y G²` in the eigenbasis. -/
private theorem gradMat_eq_sum (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d)))
    (k : Fin p) (j : Fin d) :
    gradMat z x k j
      = ∑ a : Fin d, ((Zmat x k a : ℂ) * (gval z x a) ^ 2) * ((eigUx x j a : ℝ) : ℂ) := by
  have hconj := R4C.resolvC_mul_resolvC_eq_conj (isHermitian_gram ((matrixEquivE p d).symm x)) hz
  have hG : Gmat z x * Gmat z x
      = R4C.cmat (eigUx x) *
        Matrix.diagonal (fun a => (gval z x a) ^ 2) * (R4C.cmat (eigUx x))ᵀ := hconj
  rw [gradMat, Matrix.mul_assoc, hG, ← Matrix.mul_assoc, ← Matrix.mul_assoc]
  rw [Matrix.mul_apply]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Matrix.mul_diagonal, Matrix.transpose_apply]
  have h1 : (cmapR ((matrixEquivE p d).symm x) * R4C.cmat (eigUx x)) k a
      = ((Zmat x k a : ℝ) : ℂ) := by
    rw [cmat_eq_cmapR, ← cmapR_mul]
    rfl
  rw [h1]
  rfl

/-- The column norms of `Z = Y U`. -/
private theorem sum_sq_Zmat (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) (a : Fin d) :
    ∑ k : Fin p, (Zmat x k a) ^ 2 = (d : ℝ) * eigVal x a := by
  set Y := (matrixEquivE p d).symm x with hY
  set U := eigUx x with hU
  have hd0 : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hYY : Yᵀ * Y = (d : ℝ) • gram Y := by
    rw [gram, smul_smul, mul_inv_cancel₀ hd0, one_smul]
  have hUW : Uᵀ * gram Y * U = Matrix.diagonal (eigVal x) := by
    have hc := eigU_conj (isHermitian_gram Y)
    calc Uᵀ * gram Y * U
        = Uᵀ * (U * Matrix.diagonal (eigVal x) * Uᵀ) * U := by
          rw [hU, eigUx, hY]
          rw [show Matrix.diagonal (eigVal x)
              = Matrix.diagonal (isHermitian_gram ((matrixEquivE p d).symm x)).eigenvalues from rfl]
          rw [hc]
      _ = (Uᵀ * U) * Matrix.diagonal (eigVal x) * (Uᵀ * U) := by
          simp only [Matrix.mul_assoc]
      _ = Matrix.diagonal (eigVal x) := by
          rw [hU, eigUx, transpose_eigU_mul, Matrix.one_mul, Matrix.mul_one]
  have hZZ : (Zmat x)ᵀ * Zmat x = (d : ℝ) • Matrix.diagonal (eigVal x) := by
    rw [Zmat, Matrix.transpose_mul, ← hU, ← hY]
    calc Uᵀ * Yᵀ * (Y * U) = Uᵀ * (Yᵀ * Y) * U := by simp only [Matrix.mul_assoc]
      _ = Uᵀ * ((d : ℝ) • gram Y) * U := by rw [hYY]
      _ = (d : ℝ) • (Uᵀ * gram Y * U) := by
          rw [Matrix.mul_smul, Matrix.smul_mul]
      _ = (d : ℝ) • Matrix.diagonal (eigVal x) := by rw [hUW]
  have h1 : ((Zmat x)ᵀ * Zmat x) a a = ∑ k : Fin p, (Zmat x k a) ^ 2 := by
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun k _ => by rw [Matrix.transpose_apply, sq]
  rw [← h1, hZZ, Matrix.smul_apply, Matrix.diagonal_apply_eq, smul_eq_mul]

/-- **The Frobenius identity of step 3.** -/
theorem sum_sq_norm_gradMat (hz : z.im ≠ 0) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ∑ k : Fin p, ∑ j : Fin d, ‖gradMat z x k j‖ ^ 2
      = ∑ a : Fin d, ((d : ℝ) * eigVal x a) * ‖gval z x a‖ ^ 4 := by
  have hrow : ∀ k : Fin p, ∑ j : Fin d, ‖gradMat z x k j‖ ^ 2
      = ∑ a : Fin d, (Zmat x k a) ^ 2 * ‖gval z x a‖ ^ 4 := by
    intro k
    rw [Finset.sum_congr rfl (fun j (_ : j ∈ Finset.univ) => by
      rw [gradMat_eq_sum hz x k j] :
      ∀ j ∈ Finset.univ, ‖gradMat z x k j‖ ^ 2
        = ‖∑ a : Fin d, ((Zmat x k a : ℂ) * (gval z x a) ^ 2) * ((eigUx x j a : ℝ) : ℂ)‖ ^ 2)]
    rw [sum_normSq_mulVec (V := eigUx x) (transpose_eigUx_mul x)
      (fun a => (Zmat x k a : ℂ) * (gval z x a) ^ 2)]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, norm_pow, mul_pow, sq_abs]
    ring
  rw [Finset.sum_congr rfl (fun k (_ : k ∈ Finset.univ) => hrow k), Finset.sum_comm]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [← Finset.sum_mul, sum_sq_Zmat hd x a]

/-- **The Frobenius bound of step 3.** -/
theorem sum_sq_norm_gradMat_le (hz : 0 < z.im) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    ∑ k : Fin p, ∑ j : Fin d, ‖gradMat z x k j‖ ^ 2
      ≤ (d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im ^ 4 := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  rw [sum_sq_norm_gradMat hz.ne' hd x]
  have hterm : ∀ a : Fin d, ((d : ℝ) * eigVal x a) * ‖gval z x a‖ ^ 4
      ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im ^ 4) := by
    intro a
    set l : ℝ := eigVal x a with hl
    set D : ℝ := ‖((l : ℂ) - z)‖ with hD
    have hDge : z.im ≤ D := by
      have h := Complex.abs_im_le_norm ((l : ℂ) - z)
      have him : ((l : ℂ) - z).im = -z.im := by simp
      rw [him, abs_neg, abs_of_pos hz] at h
      exact h
    have hDpos : 0 < D := lt_of_lt_of_le hz hDge
    have hgn : ‖gval z x a‖ = D⁻¹ := by rw [gval, norm_inv, ← hD]
    have hlnn : 0 ≤ l := by
      have hsum : ∑ k : Fin p, (Zmat x k a) ^ 2 = (d : ℝ) * l := sum_sq_Zmat hd x a
      have h0 : (0 : ℝ) ≤ ∑ k : Fin p, (Zmat x k a) ^ 2 :=
        Finset.sum_nonneg fun k _ => sq_nonneg _
      rw [hsum] at h0
      nlinarith
    have hlle : l ≤ D + ‖z‖ := by
      have h1 : ‖(l : ℂ)‖ ≤ ‖((l : ℂ) - z)‖ + ‖z‖ := by
        have := norm_add_le ((l : ℂ) - z) z
        simpa using this
      rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hlnn] at h1
      exact h1
    have hDinv4 : (D⁻¹) ^ 4 = (D ^ 4)⁻¹ := by rw [inv_pow]
    have hkey : l * (D⁻¹) ^ 4 ≤ (z.im + ‖z‖) / z.im ^ 4 := by
      have h2 : l * (D⁻¹) ^ 4 ≤ (D + ‖z‖) * (D⁻¹) ^ 4 :=
        mul_le_mul_of_nonneg_right hlle (by positivity)
      have h3 : (D + ‖z‖) * (D⁻¹) ^ 4 = D⁻¹ ^ 3 + ‖z‖ * D⁻¹ ^ 4 := by
        field_simp
      have h4 : D⁻¹ ^ 3 ≤ (z.im)⁻¹ ^ 3 :=
        pow_le_pow_left₀ (by positivity) (by exact inv_anti₀ hz hDge) 3
      have h5 : ‖z‖ * D⁻¹ ^ 4 ≤ ‖z‖ * (z.im)⁻¹ ^ 4 :=
        mul_le_mul_of_nonneg_left
          (pow_le_pow_left₀ (by positivity) (by exact inv_anti₀ hz hDge) 4) (norm_nonneg z)
      have h6 : (z.im)⁻¹ ^ 3 + ‖z‖ * (z.im)⁻¹ ^ 4 = (z.im + ‖z‖) / z.im ^ 4 := by
        field_simp
      linarith [h2, h3.le, h3.ge, h4, h5, h6.le, h6.ge]
    rw [hgn, hDinv4, ← hDinv4]
    calc (d : ℝ) * l * (D⁻¹) ^ 4 = (d : ℝ) * (l * (D⁻¹) ^ 4) := by ring
      _ ≤ (d : ℝ) * ((z.im + ‖z‖) / z.im ^ 4) :=
        mul_le_mul_of_nonneg_left hkey hdR.le
  calc ∑ a : Fin d, ((d : ℝ) * eigVal x a) * ‖gval z x a‖ ^ 4
      ≤ ∑ _a : Fin d, (d : ℝ) * ((z.im + ‖z‖) / z.im ^ 4) :=
        Finset.sum_le_sum fun a _ => hterm a
    _ = (d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im ^ 4 := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
        ring

/-! ### The gradient and the Lipschitz bound -/

/-- The Lipschitz constant of step 3 of `notes/archive/rmt_R1.md`; the same term as `R1.lipConst`.
-/
noncomputable def lipC (z : ℂ) (d : ℕ) : ℝ := 2 * Real.sqrt (z.im + ‖z‖) / (d * z.im ^ 2)

/-- `∇_Y Re s = -(2/d²) Re (Y G²)`, on the flattened space. -/
noncomputable def gradReVec (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    EuclideanSpace ℝ (Fin (p * d)) :=
  matrixEquivE p d (Matrix.of fun k j => -(2 / (d : ℝ) ^ 2) * (gradMat z x k j).re)

/-- `∇_Y Im s = -(2/d²) Im (Y G²)`. -/
noncomputable def gradImVec (z : ℂ) (x : EuclideanSpace ℝ (Fin (p * d))) :
    EuclideanSpace ℝ (Fin (p * d)) :=
  matrixEquivE p d (Matrix.of fun k j => -(2 / (d : ℝ) ^ 2) * (gradMat z x k j).im)

private theorem inner_matrixEquivE (M : Matrix (Fin p) (Fin d) ℝ)
    (h : EuclideanSpace ℝ (Fin (p * d))) :
    innerSL ℝ (matrixEquivE p d M) h
      = ∑ k : Fin p, ∑ j : Fin d, M k j * h (finProdFinEquiv (k, j)) := by
  rw [innerSL_apply_apply, inner_euclidean_eq_dotProduct, dotProduct,
    ← Equiv.sum_comp finProdFinEquiv
      (fun r => WithLp.ofLp (matrixEquivE p d M) r * WithLp.ofLp h r),
    Fintype.sum_prod_type]
  exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun j _ => by simp

theorem hasFDerivAt_sfun_re (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    HasFDerivAt (fun y => (sfun (p := p) (d := d) z y).re) (innerSL ℝ (gradReVec z x)) x := by
  have hS := (differentiableAt_sfun hz x).hasFDerivAt
  have h2 := Complex.reCLM.hasFDerivAt.comp x hS
  rw [Function.comp_def] at h2
  simp only [Complex.reCLM_apply] at h2
  have heq : Complex.reCLM.comp (fderiv ℝ (sfun (p := p) (d := d) z) x)
      = innerSL ℝ (gradReVec z x) := by
    ext h
    rw [ContinuousLinearMap.coe_comp, Function.comp_apply, Complex.reCLM_apply,
      fderiv_sfun_apply hz x h, gradReVec, inner_matrixEquivE,
      show (-(2 / (d : ℂ) ^ 2)) = (((-(2 / (d : ℝ) ^ 2)) : ℝ) : ℂ) by push_cast; ring]
    simp only [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, zero_mul, sub_zero,
      Complex.re_sum, Matrix.of_apply]
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [← heq]
  exact h2

theorem hasFDerivAt_sfun_im (hz : z.im ≠ 0) (x : EuclideanSpace ℝ (Fin (p * d))) :
    HasFDerivAt (fun y => (sfun (p := p) (d := d) z y).im) (innerSL ℝ (gradImVec z x)) x := by
  have hS := (differentiableAt_sfun hz x).hasFDerivAt
  have h2 := Complex.imCLM.hasFDerivAt.comp x hS
  rw [Function.comp_def] at h2
  simp only [Complex.imCLM_apply] at h2
  have heq : Complex.imCLM.comp (fderiv ℝ (sfun (p := p) (d := d) z) x)
      = innerSL ℝ (gradImVec z x) := by
    ext h
    rw [ContinuousLinearMap.coe_comp, Function.comp_apply, Complex.imCLM_apply,
      fderiv_sfun_apply hz x h, gradImVec, inner_matrixEquivE,
      show (-(2 / (d : ℂ) ^ 2)) = (((-(2 / (d : ℝ) ^ 2)) : ℝ) : ℂ) by push_cast; ring]
    simp only [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero,
      Complex.im_sum, Matrix.of_apply]
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [← heq]
  exact h2

theorem lipC_nonneg (hz : 0 < z.im) (d : ℕ) : 0 ≤ lipC z d := by
  have h1 : 0 ≤ z.im + ‖z‖ := by positivity
  unfold lipC
  positivity

private theorem norm_grad_le_aux (hz : 0 < z.im) (hd : 0 < d)
    (x : EuclideanSpace ℝ (Fin (p * d))) {f : ℂ → ℝ} (hf : ∀ c : ℂ, |f c| ≤ ‖c‖) :
    ‖matrixEquivE p d
        (Matrix.of fun k j => -(2 / (d : ℝ) ^ 2) * f (gradMat z x k j))‖ ≤ lipC z d := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hnn : 0 ≤ z.im + ‖z‖ := by positivity
  set v := matrixEquivE p d
    (Matrix.of fun k j => -(2 / (d : ℝ) ^ 2) * f (gradMat z x k j)) with hv
  have h1 : ‖v‖ ^ 2
      = ∑ k : Fin p, ∑ j : Fin d, (-(2 / (d : ℝ) ^ 2) * f (gradMat z x k j)) ^ 2 :=
    norm_matrixEquivE_sq _
  have h2 : ∀ (k : Fin p) (j : Fin d),
      (-(2 / (d : ℝ) ^ 2) * f (gradMat z x k j)) ^ 2
        ≤ (2 / (d : ℝ) ^ 2) ^ 2 * ‖gradMat z x k j‖ ^ 2 := by
    intro k j
    have hb : |f (gradMat z x k j)| ≤ ‖gradMat z x k j‖ := hf _
    have hsq : (f (gradMat z x k j)) ^ 2 ≤ ‖gradMat z x k j‖ ^ 2 := by
      have := sq_abs (f (gradMat z x k j))
      nlinarith [abs_nonneg (f (gradMat z x k j)), norm_nonneg (gradMat z x k j)]
    calc (-(2 / (d : ℝ) ^ 2) * f (gradMat z x k j)) ^ 2
        = (2 / (d : ℝ) ^ 2) ^ 2 * (f (gradMat z x k j)) ^ 2 := by ring
      _ ≤ (2 / (d : ℝ) ^ 2) ^ 2 * ‖gradMat z x k j‖ ^ 2 :=
        mul_le_mul_of_nonneg_left hsq (by positivity)
  have h3 : ‖v‖ ^ 2 ≤ (2 / (d : ℝ) ^ 2) ^ 2 *
      ∑ k : Fin p, ∑ j : Fin d, ‖gradMat z x k j‖ ^ 2 := by
    rw [h1, Finset.mul_sum]
    refine Finset.sum_le_sum fun k _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_le_sum fun j _ => h2 k j
  have h4 : ‖v‖ ^ 2 ≤ (2 / (d : ℝ) ^ 2) ^ 2 * ((d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im ^ 4) :=
    h3.trans (mul_le_mul_of_nonneg_left (sum_sq_norm_gradMat_le hz hd x) (by positivity))
  have hd0 : (d : ℝ) ≠ 0 := ne_of_gt hdR
  have hη0 : z.im ≠ 0 := ne_of_gt hz
  have h5 : (2 / (d : ℝ) ^ 2) ^ 2 * ((d : ℝ) ^ 2 * (z.im + ‖z‖) / z.im ^ 4)
      = lipC z d ^ 2 := by
    have hs : Real.sqrt (z.im + ‖z‖) ^ 2 = z.im + ‖z‖ := Real.sq_sqrt hnn
    have hlip : lipC z d ^ 2 = 4 * (z.im + ‖z‖) / ((d : ℝ) ^ 2 * z.im ^ 4) := by
      unfold lipC
      rw [div_pow, mul_pow, mul_pow, hs]
      ring
    rw [hlip]
    field_simp
    ring
  rw [h5] at h4
  calc ‖v‖ = Real.sqrt (‖v‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt (lipC z d ^ 2) := Real.sqrt_le_sqrt h4
    _ = lipC z d := Real.sqrt_sq (lipC_nonneg hz d)

theorem norm_gradReVec_le (hz : 0 < z.im) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) :
    ‖gradReVec z x‖ ≤ lipC z d :=
  norm_grad_le_aux hz hd x (fun c => Complex.abs_re_le_norm c)

theorem norm_gradImVec_le (hz : 0 < z.im) (hd : 0 < d) (x : EuclideanSpace ℝ (Fin (p * d))) :
    ‖gradImVec z x‖ ≤ lipC z d :=
  norm_grad_le_aux hz hd x (fun c => Complex.abs_im_le_norm c)

/-- **Step 3 of `notes/archive/rmt_R1.md`, real part.** -/
theorem lipschitzWith_sfun_re (hz : 0 < z.im) (hd : 0 < d) {C : ℝ≥0}
    (hC : lipC z d ≤ (C : ℝ)) :
    LipschitzWith C (fun x : EuclideanSpace ℝ (Fin (p * d)) => (sfun z x).re) := by
  refine lipschitzWith_of_nnnorm_fderiv_le
    (fun x => (hasFDerivAt_sfun_re hz.ne' x).differentiableAt) ?_
  intro x
  rw [(hasFDerivAt_sfun_re hz.ne' x).fderiv, ← NNReal.coe_le_coe, coe_nnnorm,
    innerSL_apply_norm]
  exact (norm_gradReVec_le hz hd x).trans hC

/-- **Step 3 of `notes/archive/rmt_R1.md`, imaginary part.** -/
theorem lipschitzWith_sfun_im (hz : 0 < z.im) (hd : 0 < d) {C : ℝ≥0}
    (hC : lipC z d ≤ (C : ℝ)) :
    LipschitzWith C (fun x : EuclideanSpace ℝ (Fin (p * d)) => (sfun z x).im) := by
  refine lipschitzWith_of_nnnorm_fderiv_le
    (fun x => (hasFDerivAt_sfun_im hz.ne' x).differentiableAt) ?_
  intro x
  rw [(hasFDerivAt_sfun_im hz.ne' x).fderiv, ← NNReal.coe_le_coe, coe_nnnorm,
    innerSL_apply_norm]
  exact (norm_gradImVec_le hz hd x).trans hC

/-! ### The entries of `G`, for the Stein step

`Prob/GaussianAdapters.lean`'s `integral_entry_mul_gaussianMatrix_complex` takes
`ContDiff ℝ 1 F` and `∀ x, ‖fderiv ℝ F x‖ ≤ L` for `F : EuclideanSpace ℝ (Fin (p*d)) → ℂ`.
The first hypothesis is `contDiff_Gmat_entry` below; the entry derivative is
`fderiv_Gmat_entry`. -/

/-- One entry of an operator, as a continuous real-linear functional. -/
noncomputable def entryOp (d : ℕ) (i j : Fin d) :
    ((Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) →L[ℝ] ℂ :=
  LinearMap.toContinuousLinearMap (entryOpₗ d i j)

@[simp] theorem entryOp_apply (i j : Fin d) (T : (Fin d → ℂ) →L[ℝ] (Fin d → ℂ)) :
    entryOp d i j T = entryOpₗ d i j T := rfl

theorem contDiff_Gmat_entry (hz : z.im ≠ 0) {n : WithTop ℕ∞} (i j : Fin d) :
    ContDiff ℝ n (fun x : EuclideanSpace ℝ (Fin (p * d)) => Gmat z x i j) := by
  have h1 : ContDiff ℝ n (Gop (p := p) (d := d) z) := contDiff_Gop hz
  have h2 : ContDiff ℝ n (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
      (Gop (p := p) (d := d) z x) (Pi.single j (1 : ℂ))) := h1.clm_apply contDiff_const
  have h3 := (contDiff_pi.mp h2) i
  have heq : (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
      (Gop (p := p) (d := d) z x) (Pi.single j (1 : ℂ)) i) = fun x => Gmat z x i j := by
    funext x
    exact entryOpₗ_toOpₗ (Gmat z x) i j
  rwa [heq] at h3

private theorem Gmat_entry_deriv_aux (hz : z.im ≠ 0) (i j : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (fun y : EuclideanSpace ℝ (Fin (p * d)) => Gmat z y i j) x ∧
      ∀ h, fderiv ℝ (fun y : EuclideanSpace ℝ (Fin (p * d)) => Gmat z y i j) x h
        = -((Gmat z x * dWmat p d x h * Gmat z x) i j) := by
  have hinv := hasFDerivAt_ringInverse (𝕜 := ℝ) (AopUnit hz x)
  have hG0 := hinv.comp x (hasFDerivAt_Aop (p := p) (d := d) z x)
  rw [Function.comp_def, Gop_eq_comp hz] at hG0
  have hT := ((entryOp d i j).hasFDerivAt).comp x hG0
  rw [Function.comp_def] at hT
  have hE : (fun y : EuclideanSpace ℝ (Fin (p * d)) => entryOp d i j (Gop z y))
      = fun y => Gmat z y i j := by
    funext y
    exact entryOpₗ_toOpₗ (Gmat z y) i j
  rw [hE] at hT
  refine ⟨hT.differentiableAt, fun h => ?_⟩
  rw [hT.fderiv]
  simp only [ContinuousLinearMap.coe_comp, Function.comp_apply, neg_apply,
    ContinuousLinearMap.mulLeftRight_apply, entryOp_apply, dAop_apply]
  rw [coe_AopUnit_inv, Gop, ← toOpₗ_mul, ← toOpₗ_mul, map_neg, entryOpₗ_toOpₗ]

theorem differentiableAt_Gmat_entry (hz : z.im ≠ 0) (i j : Fin d)
    (x : EuclideanSpace ℝ (Fin (p * d))) :
    DifferentiableAt ℝ (fun y : EuclideanSpace ℝ (Fin (p * d)) => Gmat z y i j) x :=
  (Gmat_entry_deriv_aux hz i j x).1

/-- **The entry derivative** `∂G = -G (∂W) G`. -/
theorem fderiv_Gmat_entry (hz : z.im ≠ 0) (i j : Fin d)
    (x h : EuclideanSpace ℝ (Fin (p * d))) :
    fderiv ℝ (fun y : EuclideanSpace ℝ (Fin (p * d)) => Gmat z y i j) x h
      = -((Gmat z x * dWmat p d x h * Gmat z x) i j) :=
  (Gmat_entry_deriv_aux hz i j x).2 h

theorem contDiff_sfun (hz : z.im ≠ 0) {n : WithTop ℕ∞} :
    ContDiff ℝ n (sfun (p := p) (d := d) z) := by
  have heq : sfun (p := p) (d := d) z
      = fun x : EuclideanSpace ℝ (Fin (p * d)) => (d : ℂ)⁻¹ * ∑ i : Fin d, Gmat z x i i := by
    funext x
    rw [sfun, traceOp_apply, Gop, traceOpₗ_toOpₗ]
    rfl
  rw [heq]
  exact contDiff_const.mul (ContDiff.sum fun i _ => contDiff_Gmat_entry hz i i)

end ResolvDeriv








end StackedSVD

