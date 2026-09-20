/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RMT.R4C
import StackedSVD.RMT.MP7
import StackedSVD.RMT.Symmetry
import StackedSVD.Prob.GaussianAdapters
import StackedSVD.Prob.TendstoInProb

/-!
# Item R2: isotropic resolvent forms at complex `z`

Review note: `notes/archive/rmt_R2.md` (restructured 2026-08-30, decisions D6, D7). This file proves
the six `TendstoInProb` limits of R2d, which item T consumes.

## Route

`notes/archive/rmt_R2.md` proves the two `v`-forms with a gradient bound on `Y ↦ vᵀ G(d⁻¹ YᵀY) v`
and Gaussian Lipschitz concentration. This file replaces that step by a **symmetry** step,
which needs no matrix calculus (see `Statement changes` in
`notes/archive/agent_reports/proof_r2.md`):

1. `resolvC_conj`: for `Oᵀ O = 1`, `G(Oᵀ W O) = Oᵀ G(W) O`, so every form conjugates and the
   trace is invariant. With `W0 (Y O) = Oᵀ (W0 Y) O` and `measurePreserving_mul_right`
   (item Sym) the law of `‖vᵀ G v - c‖` is the same for every unit `v`.
2. The Householder reflection of item Sym maps `v` to `u := x/‖x‖` for every `x ≠ 0`, so the
   law of `‖vᵀ G v - c‖` is the law of `‖uᵀ G u - c‖` with `u` the direction of an
   independent Gaussian vector, that is of `‖(gᵀ G g)/(g ⬝ᵥ g) - c‖` (`meas_qform_eq`).
3. The `g`-forms are handled by Chebyshev in the eigenbasis of `W0`: with `y = Uᵀ x` again
   standard Gaussian (`measurePreserving_mulVec`),
   `gᵀ G g - s_N = d⁻¹ ∑ (λ_a - z)⁻¹ (y_a² - 1)` and `vᵀ G g = d^{-1/2} ∑ (λ_a - z)⁻¹ k_a y_a`,
   whose second moments are `varSq/d²·∑|w_a|²` and `d⁻¹ ∑ |w_a k_a|²`. Only
   `∫ y y' = δ` and `∫ (y²-1)(y'²-1) = varSq δ` enter, where
   `varSq = ∫ (t²-1)² ∂N(0,1)` stays an abstract finite constant (its value `2` is never
   needed).
4. `g ⬝ᵥ g → 1` in probability comes from the same second moment at `w = 1`, and closes
   step 2.

R2a (`E[G] = E[s] I` by signed permutations) is proved as the note states it, from the same
conjugation lemma; the assembly does not use it, since step 3 centers at the random `s_N`
exactly.

STATUS: see `notes/archive/agent_reports/proof_r2.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace R2

variable {p d : ℕ} {z : ℂ} {W : Matrix (Fin d) (Fin d) ℝ}

/-! ### Definitions -/

/-- The Wishart block `W₀ = d⁻¹ Yᵀ Y`. Same definition as `R1.W0`; R2 does not import R1
(the two files are written in parallel), see `Statement changes`. -/
noncomputable def W0 (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  ((d : ℝ))⁻¹ • (Yᵀ * Y)

/-- `g = d^{-1/2} x`, the scaled Gaussian vector. -/
noncomputable def gOf (x : Fin d → ℝ) : Fin d → ℝ := fun j => (Real.sqrt d)⁻¹ * x j

/-- The noise pair `(x, Y)` of `SpikedModel.exists_block_hasLaw`. -/
abbrev NoiseSpace (p d : ℕ) := (Fin d → ℝ) × Matrix (Fin p) (Fin d) ℝ

/-- The product law of the pair. -/
noncomputable def noiseLaw (p d : ℕ) : Measure (NoiseSpace p d) :=
  (Measure.pi fun _ : Fin d => gaussianReal 0 1).prod (gaussianMatrix p d)

/-- Bilinear form of a general complex matrix in real vectors. -/
noncomputable def bil (B : Matrix (Fin d) (Fin d) ℂ) (x y : Fin d → ℝ) : ℂ :=
  R4C.cvec x ⬝ᵥ (B *ᵥ R4C.cvec y)

theorem bil_eq_sum (B : Matrix (Fin d) (Fin d) ℂ) (x y : Fin d → ℝ) :
    bil B x y = ∑ i, (x i : ℂ) * ∑ j, B i j * (y j : ℂ) := rfl

theorem cformC_eq_bil (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x y : Fin d → ℝ) :
    R4C.cformC W z x y = bil (R4C.resolvC W z) x y := rfl

theorem cform2C_eq_bil (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x y : Fin d → ℝ) :
    R4C.cform2C W z x y = bil (R4C.resolvC W z * R4C.resolvC W z) x y := rfl

theorem isHermitian_W0 (Y : Matrix (Fin p) (Fin d) ℝ) : (W0 Y).IsHermitian :=
  (isHermitian_transpose_mul_self Y).smul (star_trivial _)

theorem W0_apply (Y : Matrix (Fin p) (Fin d) ℝ) (i j : Fin d) :
    W0 Y i j = ((d : ℝ))⁻¹ * ∑ k, Y k i * Y k j := rfl

/-- Right rotation conjugates the block. -/
theorem W0_mul (Y : Matrix (Fin p) (Fin d) ℝ) (O : Matrix (Fin d) (Fin d) ℝ) :
    W0 (Y * O) = Oᵀ * W0 Y * O := by
  rw [W0, W0, Matrix.transpose_mul, Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  simp only [Matrix.mul_assoc]

/-! ### The complex resolvent under conjugation -/

theorem cmat_mulVec_cvec (O : Matrix (Fin d) (Fin d) ℝ) (x : Fin d → ℝ) :
    R4C.cmat O *ᵥ R4C.cvec x = R4C.cvec (O *ᵥ x) := by
  funext a
  simp only [R4C.cvec, R4C.cmat, Matrix.mulVec, dotProduct, Matrix.map_apply]
  push_cast
  rfl

/-- `G(z)` is a right inverse of `W - z`. -/
theorem mul_resolvC (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    (R4C.cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)) * R4C.resolvC W z = 1 := by
  have hne := R4C.eigenvalue_sub_ne_zero hW hz
  have hDD : Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
      Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [Matrix.diagonal_mul_diagonal,
      show (fun a => ((hW.eigenvalues a : ℂ) - z) * ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = fun _ => (1 : ℂ) from funext fun a => mul_inv_cancel₀ (hne a)]
    exact Matrix.diagonal_one
  rw [R4C.cmat_sub_smul_one_eq_conj hW, R4C.resolvC_eq_conj hW hz]
  have hassoc : R4C.cmat (R4.eigU hW) *
        Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) * (R4C.cmat (R4.eigU hW))ᵀ *
      (R4C.cmat (R4.eigU hW) *
        Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
          (R4C.cmat (R4.eigU hW))ᵀ)
      = R4C.cmat (R4.eigU hW) * (Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
          ((R4C.cmat (R4.eigU hW))ᵀ * R4C.cmat (R4.eigU hW)) *
          Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)) *
        (R4C.cmat (R4.eigU hW))ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, R4C.transpose_ceigU_mul hW, Matrix.mul_one, hDD, Matrix.mul_one,
    R4C.ceigU_mul_transpose hW]

section Conj

variable {O : Matrix (Fin d) (Fin d) ℝ}

theorem transpose_cmat_mul (hO : Oᵀ * O = 1) :
    (R4C.cmat O)ᵀ * R4C.cmat O = 1 := by
  rw [← R4C.cmat_transpose, ← R4C.cmat_mul, hO, R4C.cmat_one]

theorem cmat_mul_transpose (hO : Oᵀ * O = 1) :
    R4C.cmat O * (R4C.cmat O)ᵀ = 1 := by
  rw [← R4C.cmat_transpose, ← R4C.cmat_mul, mul_transpose_of_orth hO, R4C.cmat_one]

/-- **The conjugation lemma.** `G(Oᵀ W O) = Oᵀ G(W) O` for `Oᵀ O = 1`. -/
theorem resolvC_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    R4C.resolvC (Oᵀ * W * O) z = (R4C.cmat O)ᵀ * R4C.resolvC W z * R4C.cmat O := by
  have hcO := transpose_cmat_mul hO
  have hcOO := cmat_mul_transpose hO
  refine Matrix.inv_eq_right_inv ?_
  have hsub : R4C.cmat (Oᵀ * W * O) - z • (1 : Matrix (Fin d) (Fin d) ℂ)
      = (R4C.cmat O)ᵀ * (R4C.cmat W - z • 1) * R4C.cmat O := by
    rw [R4C.cmat_mul, R4C.cmat_mul, R4C.cmat_transpose, Matrix.mul_sub, Matrix.sub_mul,
      Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul, hcO]
  rw [hsub]
  have key : (R4C.cmat O)ᵀ * (R4C.cmat W - z • 1) * R4C.cmat O *
        ((R4C.cmat O)ᵀ * R4C.resolvC W z * R4C.cmat O)
      = (R4C.cmat O)ᵀ * ((R4C.cmat W - z • 1) *
          (R4C.cmat O * (R4C.cmat O)ᵀ) * R4C.resolvC W z) * R4C.cmat O := by
    simp only [Matrix.mul_assoc]
  rw [key, hcOO, Matrix.mul_one, mul_resolvC hW hz, Matrix.mul_one, hcO]

theorem resolvC_sq_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    R4C.resolvC (Oᵀ * W * O) z * R4C.resolvC (Oᵀ * W * O) z
      = (R4C.cmat O)ᵀ * (R4C.resolvC W z * R4C.resolvC W z) * R4C.cmat O := by
  rw [resolvC_conj hO hW hz]
  have key : (R4C.cmat O)ᵀ * R4C.resolvC W z * R4C.cmat O *
        ((R4C.cmat O)ᵀ * R4C.resolvC W z * R4C.cmat O)
      = (R4C.cmat O)ᵀ * (R4C.resolvC W z *
          (R4C.cmat O * (R4C.cmat O)ᵀ) * R4C.resolvC W z) * R4C.cmat O := by
    simp only [Matrix.mul_assoc]
  rw [key, cmat_mul_transpose hO, Matrix.mul_one]

/-- A conjugated kernel gives the form at the rotated vectors. -/
theorem bil_conj (O : Matrix (Fin d) (Fin d) ℝ) (B : Matrix (Fin d) (Fin d) ℂ)
    (x y : Fin d → ℝ) :
    bil ((R4C.cmat O)ᵀ * B * R4C.cmat O) x y = bil B (O *ᵥ x) (O *ᵥ y) := by
  have hmv : ((R4C.cmat O)ᵀ * B * R4C.cmat O) *ᵥ R4C.cvec y
      = (R4C.cmat O)ᵀ *ᵥ (B *ᵥ R4C.cvec (O *ᵥ y)) := by
    rw [← cmat_mulVec_cvec, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]
  rw [bil, bil, hmv, Matrix.dotProduct_mulVec, Matrix.vecMul_transpose, cmat_mulVec_cvec]

theorem cformC_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin d → ℝ) :
    R4C.cformC (Oᵀ * W * O) z x y = R4C.cformC W z (O *ᵥ x) (O *ᵥ y) := by
  rw [cformC_eq_bil, cformC_eq_bil, resolvC_conj hO hW hz, bil_conj]

theorem cform2C_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin d → ℝ) :
    R4C.cform2C (Oᵀ * W * O) z x y = R4C.cform2C W z (O *ᵥ x) (O *ᵥ y) := by
  rw [cform2C_eq_bil, cform2C_eq_bil, resolvC_sq_conj hO hW hz, bil_conj]

theorem stieltjesC_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    R4C.stieltjesC (Oᵀ * W * O) z = R4C.stieltjesC W z := by
  rw [R4C.stieltjesC, R4C.stieltjesC, resolvC_conj hO hW hz, Matrix.trace_mul_comm,
    ← Matrix.mul_assoc, cmat_mul_transpose hO, Matrix.one_mul]

theorem stieltjes2C_conj (hO : Oᵀ * O = 1) (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    R4C.stieltjes2C (Oᵀ * W * O) z = R4C.stieltjes2C W z := by
  rw [R4C.stieltjes2C, R4C.stieltjes2C, resolvC_sq_conj hO hW hz, Matrix.trace_mul_comm,
    ← Matrix.mul_assoc, cmat_mul_transpose hO, Matrix.one_mul]

end Conj

/-! ### Measurability of the resolvent entries in the block -/

section Measurability

variable {α : Type*} [MeasurableSpace α]

theorem measurable_det {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) : Measurable fun a => (M a).det := by
  simp only [Matrix.det_apply']
  refine Finset.measurable_sum _ fun σ _ => ?_
  exact measurable_const.mul (Finset.measurable_prod _ fun i _ => hM (σ i) i)

theorem measurable_adjugate {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin n) :
    Measurable fun a => (M a).adjugate i j := by
  simp only [Matrix.adjugate_apply]
  refine measurable_det fun k l => ?_
  by_cases h : k = j
  · subst h
    simp only [Matrix.updateRow_self]
    exact measurable_const
  · simp only [Matrix.updateRow_ne h]
    exact hM k l

theorem measurable_inv_entry {n : ℕ} {M : α → Matrix (Fin n) (Fin n) ℂ}
    (hM : ∀ i j, Measurable fun a => M a i j) (i j : Fin n) :
    Measurable fun a => (M a)⁻¹ i j := by
  simp only [Matrix.inv_def, Matrix.smul_apply, smul_eq_mul, Ring.inverse_eq_inv']
  exact ((measurable_det hM).inv).mul (measurable_adjugate hM i j)

theorem measurable_matrix_entry (k : Fin p) (l : Fin d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => Y k l :=
  (measurable_pi_apply l).comp (measurable_pi_apply k)

theorem measurable_W0_entry (i j : Fin d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => W0 Y i j := by
  simp only [W0_apply]
  refine measurable_const.mul (Finset.measurable_sum _ fun k _ => ?_)
  exact (measurable_matrix_entry k i).mul (measurable_matrix_entry k j)

theorem measurable_resolvC_entry (z : ℂ) (i j : Fin d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.resolvC (W0 Y) z i j := by
  refine measurable_inv_entry (fun k l => ?_) i j
  simp only [Matrix.sub_apply, R4C.cmat, Matrix.map_apply, Matrix.smul_apply, smul_eq_mul]
  exact (Complex.measurable_ofReal.comp (measurable_W0_entry k l)).sub measurable_const

theorem measurable_resolvC_sq_entry (z : ℂ) (i j : Fin d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ =>
      (R4C.resolvC (W0 Y) z * R4C.resolvC (W0 Y) z) i j := by
  simp only [Matrix.mul_apply]
  exact Finset.measurable_sum _ fun k _ =>
    (measurable_resolvC_entry z i k).mul (measurable_resolvC_entry z k j)

/-- Every form of a kernel with measurable entries, at measurable vector families, is
measurable. -/
theorem measurable_bil {K : α → Matrix (Fin d) (Fin d) ℂ}
    (hK : ∀ i j, Measurable fun a => K a i j)
    {x y : α → (Fin d → ℝ)} (hx : ∀ i, Measurable fun a => x a i)
    (hy : ∀ j, Measurable fun a => y a j) :
    Measurable fun a => bil (K a) (x a) (y a) := by
  simp only [bil_eq_sum]
  refine Finset.measurable_sum _ fun i _ => ?_
  refine (Complex.measurable_ofReal.comp (hx i)).mul (Finset.measurable_sum _ fun j _ => ?_)
  exact (hK i j).mul (Complex.measurable_ofReal.comp (hy j))

theorem measurable_stieltjesC (z : ℂ) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (W0 Y) z := by
  simp only [R4C.stieltjesC, Matrix.trace, Matrix.diag]
  exact measurable_const.mul
    (Finset.measurable_sum _ fun i _ => measurable_resolvC_entry z i i)

theorem measurable_stieltjes2C (z : ℂ) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjes2C (W0 Y) z := by
  simp only [R4C.stieltjes2C, Matrix.trace, Matrix.diag]
  exact measurable_const.mul
    (Finset.measurable_sum _ fun i _ => measurable_resolvC_sq_entry z i i)

end Measurability

/-! ### Moments of the coordinate product measure -/

section Moments

/-- The law of the standard Gaussian vector `x` of `notes/archive/rmt_R0.md`. -/
noncomputable abbrev piGauss (d : ℕ) : Measure (Fin d → ℝ) :=
  Measure.pi fun _ : Fin d => gaussianReal 0 1

theorem measurePreserving_coord (j : Fin d) :
    MeasurePreserving (fun x : Fin d → ℝ => x j) (piGauss d) (gaussianReal 0 1) :=
  measurePreserving_eval _ j

theorem integral_comp_coord {f : ℝ → ℝ} (hf : AEStronglyMeasurable f (gaussianReal 0 1))
    (j : Fin d) : ∫ x, f (x j) ∂(piGauss d) = ∫ t, f t ∂(gaussianReal 0 1) := by
  have hmp := measurePreserving_coord (d := d) j
  have h := integral_map (φ := fun x : Fin d → ℝ => x j) (μ := piGauss d)
    (measurable_pi_apply j).aemeasurable (f := f) (by rw [hmp.map_eq]; exact hf)
  rw [hmp.map_eq] at h
  exact h.symm

theorem integrable_comp_coord {f : ℝ → ℝ} (hf : Integrable f (gaussianReal 0 1)) (j : Fin d) :
    Integrable (fun x : Fin d → ℝ => f (x j)) (piGauss d) :=
  memLp_one_iff_integrable.mp
    ((memLp_one_iff_integrable.mpr hf).comp_measurePreserving (measurePreserving_coord j))

theorem indepFun_coord {a b : Fin d} (hab : a ≠ b) :
    IndepFun (fun x : Fin d → ℝ => x a) (fun x : Fin d → ℝ => x b) (piGauss d) :=
  (iIndepFun_pi (X := fun _ : Fin d => (id : ℝ → ℝ)) fun _ => aemeasurable_id).indepFun hab

/-- Every power of the coordinate is integrable under the one dimensional Gaussian. -/
theorem integrable_pow_gauss (n : ℕ) : Integrable (fun t : ℝ => t ^ n) (gaussianReal 0 1) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  have hm : MemLp (id : ℝ → ℝ) (n : ℝ≥0∞) (gaussianReal 0 1) :=
    memLp_id_gaussianReal' _ (by simp)
  have h := hm.integrable_norm_pow hn.ne'
  rw [← integrable_norm_iff (by fun_prop)]
  simpa [norm_pow] using h

theorem integrable_sq_sub_one :
    Integrable (fun t : ℝ => (t ^ 2 - 1) ^ 2) (gaussianReal 0 1) := by
  have h4 := integrable_pow_gauss 4
  have h2 := integrable_pow_gauss 2
  have heq : (fun t : ℝ => (t ^ 2 - 1) ^ 2) = fun t : ℝ => t ^ 4 - 2 * t ^ 2 + 1 := by
    funext t; ring
  rw [heq]
  exact (h4.sub (h2.const_mul 2)).add (integrable_const 1)

/-- `∫ (t² - 1)² ∂N(0,1)`. Its value is `2`; only that it is a finite constant is used. -/
noncomputable def varSq : ℝ := ∫ t, (t ^ 2 - 1) ^ 2 ∂(gaussianReal 0 1)

theorem varSq_nonneg : 0 ≤ varSq :=
  integral_nonneg fun _ => sq_nonneg _

theorem integral_id_gauss : ∫ t : ℝ, t ∂(gaussianReal 0 1) = 0 := by
  simp

theorem integral_sq_gauss : ∫ t : ℝ, t ^ 2 ∂(gaussianReal 0 1) = 1 := by
  have hm : MemLp (id : ℝ → ℝ) 2 (gaussianReal 0 1) := memLp_id_gaussianReal' 2 (by simp)
  have h := variance_eq_sub (μ := gaussianReal 0 1) hm
  rw [variance_id_gaussianReal] at h
  simp only [Pi.pow_apply, id_eq] at h
  rw [integral_id_gaussianReal] at h
  simpa using h.symm

theorem integral_sub_one_gauss : ∫ t : ℝ, (t ^ 2 - 1) ∂(gaussianReal 0 1) = 0 := by
  rw [integral_sub (integrable_pow_gauss 2) (integrable_const 1), integral_sq_gauss]
  simp

/-! #### One centered function of one coordinate, summed against weights -/

/-- The hypotheses shared by the two kernels `h = id` and `h = t² - 1`. -/
structure Centered (h : ℝ → ℝ) : Prop where
  meas : Measurable h
  int : Integrable h (gaussianReal 0 1)
  sqInt : Integrable (fun t => h t ^ 2) (gaussianReal 0 1)
  zero : ∫ t, h t ∂(gaussianReal 0 1) = 0

theorem centered_id : Centered (fun t : ℝ => t) where
  meas := measurable_id
  int := by simpa using integrable_pow_gauss 1
  sqInt := by simpa using integrable_pow_gauss 2
  zero := integral_id_gauss

theorem centered_sq_sub_one : Centered (fun t : ℝ => t ^ 2 - 1) where
  meas := (measurable_id.pow_const 2).sub measurable_const
  int := (integrable_pow_gauss 2).sub (integrable_const 1)
  sqInt := integrable_sq_sub_one
  zero := integral_sub_one_gauss

variable {h : ℝ → ℝ}

theorem Centered.integrable_pair (hh : Centered h) (a b : Fin d) :
    Integrable (fun x : Fin d → ℝ => h (x a) * h (x b)) (piGauss d) := by
  rcases eq_or_ne a b with rfl | hab
  · have hI := integrable_comp_coord hh.sqInt a
    simpa [pow_two] using hI
  · have hind := (indepFun_coord hab).comp hh.meas hh.meas
    simp only [Function.comp_def] at hind
    exact hind.integrable_mul (integrable_comp_coord hh.int a)
      (integrable_comp_coord hh.int b)

theorem Centered.integral_pair (hh : Centered h) (a b : Fin d) :
    ∫ x : Fin d → ℝ, h (x a) * h (x b) ∂(piGauss d)
      = if a = b then ∫ t, h t ^ 2 ∂(gaussianReal 0 1) else 0 := by
  rcases eq_or_ne a b with rfl | hab
  · have hI := integral_comp_coord hh.sqInt.aestronglyMeasurable a
    simpa [pow_two] using hI
  · rw [if_neg hab]
    have hind := (indepFun_coord hab).comp hh.meas hh.meas
    simp only [Function.comp_def] at hind
    have hmul := hind.integral_fun_mul_eq_mul_integral
      (integrable_comp_coord hh.int a).aestronglyMeasurable
      (integrable_comp_coord hh.int b).aestronglyMeasurable
    rw [hmul, integral_comp_coord hh.int.aestronglyMeasurable a, hh.zero, zero_mul]

private theorem sq_sum_expand (r u : Fin d → ℝ) :
    (∑ a, r a * u a) ^ 2 = ∑ a, ∑ b, r a * r b * (u a * u b) := by
  rw [sq, Finset.sum_mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

theorem Centered.integrable_sq_sum (hh : Centered h) (r : Fin d → ℝ) :
    Integrable (fun x : Fin d → ℝ => (∑ a, r a * h (x a)) ^ 2) (piGauss d) := by
  simp only [sq_sum_expand]
  refine integrable_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ => ?_
  exact ((hh.integrable_pair a b).const_mul (r a * r b))

theorem Centered.integral_sq_sum (hh : Centered h) (r : Fin d → ℝ) :
    ∫ x : Fin d → ℝ, (∑ a, r a * h (x a)) ^ 2 ∂(piGauss d)
      = (∫ t, h t ^ 2 ∂(gaussianReal 0 1)) * ∑ a, r a ^ 2 := by
  simp only [sq_sum_expand]
  rw [integral_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ =>
    ((hh.integrable_pair a b).const_mul (r a * r b))]
  have hb : ∀ a : Fin d, ∫ x : Fin d → ℝ, ∑ b, r a * r b * (h (x a) * h (x b)) ∂(piGauss d)
      = (∫ t, h t ^ 2 ∂(gaussianReal 0 1)) * r a ^ 2 := by
    intro a
    rw [integral_finsetSum _ fun b _ => ((hh.integrable_pair a b).const_mul (r a * r b))]
    have : ∀ b : Fin d, ∫ x : Fin d → ℝ, r a * r b * (h (x a) * h (x b)) ∂(piGauss d)
        = r a * r b * (if a = b then ∫ t, h t ^ 2 ∂(gaussianReal 0 1) else 0) := by
      intro b
      rw [integral_const_mul, hh.integral_pair a b]
    simp only [this, mul_ite, mul_zero]
    rw [Finset.sum_ite_eq Finset.univ a fun b => r a * r b * ∫ t, h t ^ 2 ∂(gaussianReal 0 1)]
    simp only [Finset.mem_univ, if_true]
    ring
  rw [Finset.sum_congr rfl fun a _ => hb a, ← Finset.mul_sum]

/-! #### The complex version -/

theorem normSq_sum_eq (c : Fin d → ℂ) (u : Fin d → ℝ) :
    ‖∑ a, c a * (u a : ℂ)‖ ^ 2
      = (∑ a, (c a).re * u a) ^ 2 + (∑ a, (c a).im * u a) ^ 2 := by
  have hre : (∑ a, c a * (u a : ℂ)).re = ∑ a, (c a).re * u a := by
    rw [Complex.re_sum]
    exact Finset.sum_congr rfl fun a _ => by simp [Complex.mul_re]
  have him : (∑ a, c a * (u a : ℂ)).im = ∑ a, (c a).im * u a := by
    rw [Complex.im_sum]
    exact Finset.sum_congr rfl fun a _ => by simp [Complex.mul_im]
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, hre, him]
  ring

theorem Centered.integrable_normSq_sum (hh : Centered h) (c : Fin d → ℂ) :
    Integrable (fun x : Fin d → ℝ => ‖∑ a, c a * ((h (x a) : ℝ) : ℂ)‖ ^ 2) (piGauss d) := by
  simp only [normSq_sum_eq]
  exact (hh.integrable_sq_sum fun a => (c a).re).add (hh.integrable_sq_sum fun a => (c a).im)

theorem Centered.integral_normSq_sum (hh : Centered h) (c : Fin d → ℂ) :
    ∫ x : Fin d → ℝ, ‖∑ a, c a * ((h (x a) : ℝ) : ℂ)‖ ^ 2 ∂(piGauss d)
      = (∫ t, h t ^ 2 ∂(gaussianReal 0 1)) * ∑ a, ‖c a‖ ^ 2 := by
  simp only [normSq_sum_eq]
  rw [integral_add (hh.integrable_sq_sum fun a => (c a).re)
    (hh.integrable_sq_sum fun a => (c a).im), hh.integral_sq_sum, hh.integral_sq_sum,
    ← mul_add, ← Finset.sum_add_distrib]
  congr 1
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]
  ring

end Moments

/-! ### The forms in `g`, in the eigenbasis of `W₀` -/

section Forms

variable {W : Matrix (Fin d) (Fin d) ℝ}

theorem integral_comp_mulVec {V : Matrix (Fin d) (Fin d) ℝ} (hV : Vᵀ * V = 1)
    {F : (Fin d → ℝ) → ℝ} (hF : AEStronglyMeasurable F (piGauss d)) :
    ∫ x, F (V *ᵥ x) ∂(piGauss d) = ∫ y, F y ∂(piGauss d) := by
  have hmp := measurePreserving_mulVec (U := V) hV
  have h := integral_map (φ := fun x : Fin d → ℝ => V *ᵥ x) (μ := piGauss d)
    hmp.measurable.aemeasurable (f := F) (by rw [hmp.map_eq]; exact hF)
  rw [hmp.map_eq] at h
  exact h.symm

theorem integrable_comp_mulVec {V : Matrix (Fin d) (Fin d) ℝ} (hV : Vᵀ * V = 1)
    {F : (Fin d → ℝ) → ℝ} (hF : Integrable F (piGauss d)) :
    Integrable (fun x => F (V *ᵥ x)) (piGauss d) :=
  memLp_one_iff_integrable.mp
    ((memLp_one_iff_integrable.mpr hF).comp_measurePreserving (measurePreserving_mulVec hV))

theorem integral_normSq_comp {h : ℝ → ℝ} (hh : Centered h) {V : Matrix (Fin d) (Fin d) ℝ}
    (hV : Vᵀ * V = 1) (c : Fin d → ℂ) :
    ∫ x, ‖∑ a, c a * ((h ((V *ᵥ x) a) : ℝ) : ℂ)‖ ^ 2 ∂(piGauss d)
      = (∫ t, h t ^ 2 ∂(gaussianReal 0 1)) * ∑ a, ‖c a‖ ^ 2 := by
  rw [integral_comp_mulVec hV (F := fun y => ‖∑ a, c a * ((h (y a) : ℝ) : ℂ)‖ ^ 2)
    (hh.integrable_normSq_sum c).aestronglyMeasurable]
  exact hh.integral_normSq_sum c

theorem integrable_normSq_comp {h : ℝ → ℝ} (hh : Centered h) {V : Matrix (Fin d) (Fin d) ℝ}
    (hV : Vᵀ * V = 1) (c : Fin d → ℂ) :
    Integrable (fun x => ‖∑ a, c a * ((h ((V *ᵥ x) a) : ℝ) : ℂ)‖ ^ 2) (piGauss d) :=
  integrable_comp_mulVec hV (hh.integrable_normSq_sum c)

theorem sum_normSq_le {c : Fin d → ℂ} {M : ℝ} {r : Fin d → ℝ}
    (h : ∀ a, ‖c a‖ ≤ M * |r a|) : ∑ a, ‖c a‖ ^ 2 ≤ M ^ 2 * ∑ a, r a ^ 2 := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun a _ => ?_
  have h1 : ‖c a‖ ^ 2 ≤ (M * |r a|) ^ 2 :=
    pow_le_pow_left₀ (norm_nonneg _) (h a) 2
  calc ‖c a‖ ^ 2 ≤ (M * |r a|) ^ 2 := h1
    _ = M ^ 2 * r a ^ 2 := by rw [mul_pow, sq_abs]

theorem sqrt_inv_mul_sqrt_inv (d : ℕ) :
    (Real.sqrt d)⁻¹ * (Real.sqrt d)⁻¹ = ((d : ℝ))⁻¹ := by
  rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg d)]

theorem gOf_eq_smul (x : Fin d → ℝ) : gOf x = (Real.sqrt d)⁻¹ • x := rfl

theorem mulVec_gOf (V : Matrix (Fin d) (Fin d) ℝ) (x : Fin d → ℝ) :
    V *ᵥ gOf x = gOf (V *ᵥ x) := by
  rw [gOf_eq_smul, gOf_eq_smul, Matrix.mulVec_smul]

theorem sum_sq_transpose_eigU (hW : W.IsHermitian) (v : Fin d → ℝ) :
    ∑ a, (((R4.eigU hW)ᵀ *ᵥ v) a) ^ 2 = v ⬝ᵥ v := by
  rw [← R4.dotProduct_transpose_eigU hW v, dotProduct]
  exact Finset.sum_congr rfl fun a _ => sq _

theorem transpose_eigU_orth (hW : W.IsHermitian) :
    ((R4.eigU hW)ᵀ)ᵀ * (R4.eigU hW)ᵀ = 1 := by
  rw [Matrix.transpose_transpose]
  exact R4.eigU_mul_transpose hW

private theorem cast_mul_inv_sq (d : ℕ) (η : ℝ) :
    (d : ℝ) * (((d : ℝ))⁻¹ * η⁻¹) ^ 2 = (η ^ 2)⁻¹ * ((d : ℝ))⁻¹ := by
  rcases Nat.eq_zero_or_pos d with rfl | hd
  · simp
  · have hne : ((d : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
    field_simp

/-! #### The four expansions -/

theorem qformC_gOf_sub (hW : W.IsHermitian) (hz : z.im ≠ 0) (x : Fin d → ℝ) :
    R4C.qformC W z (gOf x) - R4C.stieltjesC W z
      = ∑ a, ((d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹) *
          (((fun t : ℝ => t ^ 2 - 1) (((R4.eigU hW)ᵀ *ᵥ x) a) : ℝ) : ℂ) := by
  have hs : R4C.stieltjesC W z = ∑ a, (d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹ := by
    rw [R4C.stieltjesC, R4C.trace_resolvC hW hz, Finset.mul_sum]
  have hq : R4C.qformC W z (gOf x)
      = ∑ a, ((d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹) *
          ((((R4.eigU hW)ᵀ *ᵥ x) a ^ 2 : ℝ) : ℂ) := by
    rw [R4C.qformC, R4C.cformC_eq_sum hW hz]
    refine Finset.sum_congr rfl fun a _ => ?_
    have hy : ((R4.eigU hW)ᵀ *ᵥ gOf x) a = (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) := by
      rw [mulVec_gOf]; rfl
    have hsq : (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) *
        ((Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a))
        = ((d : ℝ))⁻¹ * ((((R4.eigU hW)ᵀ *ᵥ x) a) ^ 2) := by
      rw [← sqrt_inv_mul_sqrt_inv d]; ring
    rw [hy, hsq]
    push_cast
    ring
  rw [hq, hs, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  push_cast
  ring

theorem qform2C_gOf_sub (hW : W.IsHermitian) (hz : z.im ≠ 0) (x : Fin d → ℝ) :
    R4C.qform2C W z (gOf x) - R4C.stieltjes2C W z
      = ∑ a, ((d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2) *
          (((fun t : ℝ => t ^ 2 - 1) (((R4.eigU hW)ᵀ *ᵥ x) a) : ℝ) : ℂ) := by
  have hs : R4C.stieltjes2C W z
      = ∑ a, (d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 := by
    rw [R4C.stieltjes2C, R4C.trace_resolvC_sq hW hz, Finset.mul_sum]
  have hq : R4C.qform2C W z (gOf x)
      = ∑ a, ((d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2) *
          ((((R4.eigU hW)ᵀ *ᵥ x) a ^ 2 : ℝ) : ℂ) := by
    rw [R4C.qform2C, R4C.cform2C_eq_sum hW hz]
    refine Finset.sum_congr rfl fun a _ => ?_
    have hy : ((R4.eigU hW)ᵀ *ᵥ gOf x) a = (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) := by
      rw [mulVec_gOf]; rfl
    have hsq : (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) *
        ((Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a))
        = ((d : ℝ))⁻¹ * ((((R4.eigU hW)ᵀ *ᵥ x) a) ^ 2) := by
      rw [← sqrt_inv_mul_sqrt_inv d]; ring
    rw [hy, hsq]
    push_cast
    ring
  rw [hq, hs, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  push_cast
  ring

theorem cformC_v_gOf (hW : W.IsHermitian) (hz : z.im ≠ 0) (v x : Fin d → ℝ) :
    R4C.cformC W z v (gOf x)
      = ∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹ *
            (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ)) *
          (((fun t : ℝ => t) (((R4.eigU hW)ᵀ *ᵥ x) a) : ℝ) : ℂ) := by
  rw [R4C.cformC_eq_sum hW hz]
  refine Finset.sum_congr rfl fun a _ => ?_
  have hy : ((R4.eigU hW)ᵀ *ᵥ gOf x) a = (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) := by
    rw [mulVec_gOf]; rfl
  rw [hy]
  push_cast
  ring

theorem cform2C_v_gOf (hW : W.IsHermitian) (hz : z.im ≠ 0) (v x : Fin d → ℝ) :
    R4C.cform2C W z v (gOf x)
      = ∑ a, ((((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
            (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ)) *
          (((fun t : ℝ => t) (((R4.eigU hW)ᵀ *ᵥ x) a) : ℝ) : ℂ) := by
  rw [R4C.cform2C_eq_sum hW hz]
  refine Finset.sum_congr rfl fun a _ => ?_
  have hy : ((R4.eigU hW)ᵀ *ᵥ gOf x) a = (Real.sqrt d)⁻¹ * (((R4.eigU hW)ᵀ *ᵥ x) a) := by
    rw [mulVec_gOf]; rfl
  rw [hy]
  push_cast
  ring

/-! #### The four second moment bounds -/

theorem norm_coef_quad_le (hW : W.IsHermitian) (hz : 0 < z.im) (a : Fin d) :
    ‖(d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹‖ ≤ (((d : ℝ))⁻¹ * (z.im)⁻¹) * |(1 : ℝ)| := by
  rw [norm_mul, abs_one, mul_one]
  have h1 : ‖(d : ℂ)⁻¹‖ = ((d : ℝ))⁻¹ := by simp
  rw [h1]
  exact mul_le_mul_of_nonneg_left (R4C.norm_inv_eigenvalue_sub_le hW hz a) (by positivity)

theorem sum_one_sq (d : ℕ) : ∑ _a : Fin d, (1 : ℝ) ^ 2 = (d : ℝ) := by simp

theorem integral_normSq_qformC_le (hW : W.IsHermitian) (hz : 0 < z.im) :
    ∫ x, ‖R4C.qformC W z (gOf x) - R4C.stieltjesC W z‖ ^ 2 ∂(piGauss d)
      ≤ (varSq / z.im ^ 2) / d := by
  simp only [fun x => qformC_gOf_sub hW hz.ne' x]
  rw [integral_normSq_comp centered_sq_sub_one (transpose_eigU_orth hW)]
  have hsum : ∑ a, ‖(d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹‖ ^ 2
      ≤ (((d : ℝ))⁻¹ * (z.im)⁻¹) ^ 2 * (d : ℝ) := by
    have := sum_normSq_le (c := fun a => (d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹)
      (M := ((d : ℝ))⁻¹ * (z.im)⁻¹) (r := fun _ => (1 : ℝ)) (norm_coef_quad_le hW hz)
    rwa [sum_one_sq] at this
  calc (∫ t, ((fun t : ℝ => t ^ 2 - 1) t) ^ 2 ∂(gaussianReal 0 1)) *
        ∑ a, ‖(d : ℂ)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹‖ ^ 2
      ≤ varSq * ((((d : ℝ))⁻¹ * (z.im)⁻¹) ^ 2 * (d : ℝ)) :=
        mul_le_mul_of_nonneg_left hsum varSq_nonneg
    _ = (varSq / z.im ^ 2) / d := by
        rw [mul_comm ((((d : ℝ))⁻¹ * (z.im)⁻¹) ^ 2) (d : ℝ), cast_mul_inv_sq]
        field_simp

theorem integrable_normSq_qformC (hW : W.IsHermitian) (hz : 0 < z.im) :
    Integrable (fun x => ‖R4C.qformC W z (gOf x) - R4C.stieltjesC W z‖ ^ 2) (piGauss d) := by
  simp only [fun x => qformC_gOf_sub hW hz.ne' x]
  exact integrable_normSq_comp centered_sq_sub_one (transpose_eigU_orth hW) _

theorem norm_coef_quad2_le (hW : W.IsHermitian) (hz : 0 < z.im) (a : Fin d) :
    ‖(d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖
      ≤ (((d : ℝ))⁻¹ * ((z.im) ^ 2)⁻¹) * |(1 : ℝ)| := by
  rw [norm_mul, abs_one, mul_one]
  have h1 : ‖(d : ℂ)⁻¹‖ = ((d : ℝ))⁻¹ := by simp
  rw [h1, norm_pow]
  refine mul_le_mul_of_nonneg_left ?_ (by positivity)
  rw [← inv_pow]
  exact pow_le_pow_left₀ (norm_nonneg _) (R4C.norm_inv_eigenvalue_sub_le hW hz a) 2

theorem integral_normSq_qform2C_le (hW : W.IsHermitian) (hz : 0 < z.im) :
    ∫ x, ‖R4C.qform2C W z (gOf x) - R4C.stieltjes2C W z‖ ^ 2 ∂(piGauss d)
      ≤ (varSq / z.im ^ 4) / d := by
  simp only [fun x => qform2C_gOf_sub hW hz.ne' x]
  rw [integral_normSq_comp centered_sq_sub_one (transpose_eigU_orth hW)]
  have hsum : ∑ a, ‖(d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖ ^ 2
      ≤ (((d : ℝ))⁻¹ * ((z.im) ^ 2)⁻¹) ^ 2 * (d : ℝ) := by
    have := sum_normSq_le (c := fun a => (d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2)
      (M := ((d : ℝ))⁻¹ * ((z.im) ^ 2)⁻¹) (r := fun _ => (1 : ℝ)) (norm_coef_quad2_le hW hz)
    rwa [sum_one_sq] at this
  have hzz : (0 : ℝ) < z.im ^ 2 := by positivity
  calc (∫ t, ((fun t : ℝ => t ^ 2 - 1) t) ^ 2 ∂(gaussianReal 0 1)) *
        ∑ a, ‖(d : ℂ)⁻¹ * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖ ^ 2
      ≤ varSq * ((((d : ℝ))⁻¹ * ((z.im) ^ 2)⁻¹) ^ 2 * (d : ℝ)) :=
        mul_le_mul_of_nonneg_left hsum varSq_nonneg
    _ = (varSq / z.im ^ 4) / d := by
        rw [mul_comm ((((d : ℝ))⁻¹ * ((z.im) ^ 2)⁻¹) ^ 2) (d : ℝ), cast_mul_inv_sq]
        have : ((z.im ^ 2) ^ 2)⁻¹ = (z.im ^ 4)⁻¹ := by
          rw [← pow_mul]
        rw [this]
        field_simp

theorem integrable_normSq_qform2C (hW : W.IsHermitian) (hz : 0 < z.im) :
    Integrable (fun x => ‖R4C.qform2C W z (gOf x) - R4C.stieltjes2C W z‖ ^ 2) (piGauss d) := by
  simp only [fun x => qform2C_gOf_sub hW hz.ne' x]
  exact integrable_normSq_comp centered_sq_sub_one (transpose_eigU_orth hW) _

theorem norm_coef_lin_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) (a : Fin d) :
    ‖((hW.eigenvalues a : ℂ) - z)⁻¹ *
        (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ)‖
      ≤ ((z.im)⁻¹ * (Real.sqrt d)⁻¹) * |((R4.eigU hW)ᵀ *ᵥ v) a| := by
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul,
    abs_of_nonneg (inv_nonneg.2 (Real.sqrt_nonneg _))]
  rw [← mul_assoc]
  exact mul_le_mul_of_nonneg_right
    (mul_le_mul_of_nonneg_right (R4C.norm_inv_eigenvalue_sub_le hW hz a)
      (inv_nonneg.2 (Real.sqrt_nonneg _))) (abs_nonneg _)

theorem integral_normSq_cformC_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) :
    ∫ x, ‖R4C.cformC W z v (gOf x)‖ ^ 2 ∂(piGauss d) ≤ ((v ⬝ᵥ v) / z.im ^ 2) / d := by
  simp only [fun x => cformC_v_gOf hW hz.ne' v x]
  rw [integral_normSq_comp centered_id (transpose_eigU_orth hW)]
  have hsum := sum_normSq_le
    (c := fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹ *
      (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ))
    (M := (z.im)⁻¹ * (Real.sqrt d)⁻¹) (r := fun a => ((R4.eigU hW)ᵀ *ᵥ v) a)
    (norm_coef_lin_le hW hz v)
  rw [sum_sq_transpose_eigU hW v] at hsum
  have hid : (∫ t, ((fun t : ℝ => t) t) ^ 2 ∂(gaussianReal 0 1)) = 1 := integral_sq_gauss
  rw [hid, one_mul]
  refine hsum.trans (le_of_eq ?_)
  rw [mul_pow, inv_pow, inv_pow, Real.sq_sqrt (Nat.cast_nonneg d)]
  field_simp

theorem integrable_normSq_cformC (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) :
    Integrable (fun x => ‖R4C.cformC W z v (gOf x)‖ ^ 2) (piGauss d) := by
  simp only [fun x => cformC_v_gOf hW hz.ne' v x]
  exact integrable_normSq_comp centered_id (transpose_eigU_orth hW) _

theorem norm_coef_lin2_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) (a : Fin d) :
    ‖(((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
        (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ)‖
      ≤ (((z.im) ^ 2)⁻¹ * (Real.sqrt d)⁻¹) * |((R4.eigU hW)ᵀ *ᵥ v) a| := by
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul,
    abs_of_nonneg (inv_nonneg.2 (Real.sqrt_nonneg _)), ← mul_assoc]
  refine mul_le_mul_of_nonneg_right ?_ (abs_nonneg _)
  refine mul_le_mul_of_nonneg_right ?_ (inv_nonneg.2 (Real.sqrt_nonneg _))
  rw [norm_pow, ← inv_pow]
  exact pow_le_pow_left₀ (norm_nonneg _) (R4C.norm_inv_eigenvalue_sub_le hW hz a) 2

theorem integral_normSq_cform2C_le (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) :
    ∫ x, ‖R4C.cform2C W z v (gOf x)‖ ^ 2 ∂(piGauss d) ≤ ((v ⬝ᵥ v) / z.im ^ 4) / d := by
  simp only [fun x => cform2C_v_gOf hW hz.ne' v x]
  rw [integral_normSq_comp centered_id (transpose_eigU_orth hW)]
  have hsum := sum_normSq_le
    (c := fun a => (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
      (((Real.sqrt d)⁻¹ * ((R4.eigU hW)ᵀ *ᵥ v) a : ℝ) : ℂ))
    (M := ((z.im) ^ 2)⁻¹ * (Real.sqrt d)⁻¹) (r := fun a => ((R4.eigU hW)ᵀ *ᵥ v) a)
    (norm_coef_lin2_le hW hz v)
  rw [sum_sq_transpose_eigU hW v] at hsum
  have hid : (∫ t, ((fun t : ℝ => t) t) ^ 2 ∂(gaussianReal 0 1)) = 1 := integral_sq_gauss
  rw [hid, one_mul]
  refine hsum.trans (le_of_eq ?_)
  rw [mul_pow, inv_pow, inv_pow, Real.sq_sqrt (Nat.cast_nonneg d), ← pow_mul]
  field_simp
  ring

theorem integrable_normSq_cform2C (hW : W.IsHermitian) (hz : 0 < z.im) (v : Fin d → ℝ) :
    Integrable (fun x => ‖R4C.cform2C W z v (gOf x)‖ ^ 2) (piGauss d) := by
  simp only [fun x => cform2C_v_gOf hW hz.ne' v x]
  exact integrable_normSq_comp centered_id (transpose_eigU_orth hW) _

/-! #### The norm of `g` -/

theorem dotProduct_gOf_sub (hd : 0 < d) (x : Fin d → ℝ) :
    gOf x ⬝ᵥ gOf x - 1 = ∑ a, ((d : ℝ))⁻¹ * ((fun t : ℝ => t ^ 2 - 1) (x a)) := by
  have hne : ((d : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have h1 : gOf x ⬝ᵥ gOf x = ∑ a, ((d : ℝ))⁻¹ * (x a) ^ 2 := by
    rw [dotProduct]
    refine Finset.sum_congr rfl fun a _ => ?_
    change ((Real.sqrt d)⁻¹ * x a) * ((Real.sqrt d)⁻¹ * x a) = _
    rw [← sqrt_inv_mul_sqrt_inv d]
    ring
  rw [h1]
  simp only [mul_sub, mul_one]
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  rw [mul_inv_cancel₀ hne]

theorem integral_sq_dotProduct_gOf_le (hd : 0 < d) :
    ∫ x, (gOf x ⬝ᵥ gOf x - 1) ^ 2 ∂(piGauss d) ≤ varSq / d := by
  simp only [fun x => dotProduct_gOf_sub hd x]
  rw [Centered.integral_sq_sum centered_sq_sub_one]
  have hne : ((d : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hsum : ∑ _a : Fin d, (((d : ℝ))⁻¹) ^ 2 = ((d : ℝ))⁻¹ := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    field_simp
  rw [hsum]
  simp only [varSq]
  rw [div_eq_mul_inv]

theorem integrable_sq_dotProduct_gOf :
    Integrable (fun x => (gOf x ⬝ᵥ gOf x - 1) ^ 2) (piGauss d) := by
  rcases Nat.eq_zero_or_pos d with rfl | hd
  · have hconst : (fun x : Fin 0 → ℝ => (gOf x ⬝ᵥ gOf x - 1) ^ 2) = fun _ => 1 := by
      funext x
      simp [dotProduct]
    rw [hconst]
    exact integrable_const 1
  · simp only [fun x => dotProduct_gOf_sub hd x]
    exact Centered.integrable_sq_sum centered_sq_sub_one _

end Forms

/-! ### Chebyshev, and the transfer to the model spaces -/

section Transfer

/-- Chebyshev for a nonnegative random variable with a bounded second moment. -/
theorem meas_ge_le_of_integral_sq {α : Type*} [MeasurableSpace α] (ν : Measure α)
    {F : α → ℝ} (hFm : Measurable F)
    (hint : Integrable (fun a => F a ^ 2) ν) {ε C : ℝ} (hε : 0 < ε)
    (hC : ∫ a, F a ^ 2 ∂ν ≤ C) :
    ν {a | ε ≤ F a} ≤ ENNReal.ofReal (C / ε ^ 2) := by
  have hsub : {a | ε ≤ F a} ⊆ {a | ENNReal.ofReal (ε ^ 2) ≤ ENNReal.ofReal (F a ^ 2)} := by
    intro a ha
    exact ENNReal.ofReal_le_ofReal (pow_le_pow_left₀ hε.le ha 2)
  refine (measure_mono hsub).trans ?_
  have hme : AEMeasurable (fun a => ENNReal.ofReal (F a ^ 2)) ν :=
    (ENNReal.measurable_ofReal.comp (hFm.pow_const 2)).aemeasurable
  have hε2 : ENNReal.ofReal (ε ^ 2) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    positivity
  refine (meas_ge_le_lintegral_div hme hε2 ENNReal.ofReal_ne_top).trans ?_
  have hli : ∫⁻ a, ENNReal.ofReal (F a ^ 2) ∂ν = ENNReal.ofReal (∫ a, F a ^ 2 ∂ν) :=
    (ofReal_integral_eq_lintegral_ofReal hint
      (Filter.Eventually.of_forall fun a => sq_nonneg _)).symm
  rw [hli, ← ENNReal.ofReal_div_of_pos (by positivity)]
  exact ENNReal.ofReal_le_ofReal (by gcongr)

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN dN : ℕ → ℕ}

theorem tendsto_ofReal_div_atTop {K : ℝ} (hd : Tendsto dN atTop atTop) (ε : ℝ) :
    Tendsto (fun N => ENNReal.ofReal (K / dN N / ε ^ 2)) atTop (𝓝 0) := by
  have hdR : Tendsto (fun N => ((dN N : ℝ))) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have h1 : Tendsto (fun N => K / (dN N : ℝ)) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop hdR
  have h2 : Tendsto (fun N => K / (dN N : ℝ) / ε ^ 2) atTop (𝓝 0) := by
    simpa using h1.div_const (ε ^ 2)
  simpa using ENNReal.tendsto_ofReal h2

/-- **The transfer step.** A second moment bound in `x`, uniform in the block `Y`, gives a
limit in probability on the model spaces. `Measure.prod_apply_symm` integrates the sections
against the block law (modeling choice 4). -/
theorem tendstoInProb_of_integral_sq_le
    (ZZ : ∀ N, Ω N → NoiseSpace (pN N) (dN N))
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (F : ∀ N, (Fin (dN N) → ℝ) → Matrix (Fin (pN N)) (Fin (dN N)) ℝ → ℝ)
    (hFm : ∀ N, Measurable fun q : NoiseSpace (pN N) (dN N) => F N q.1 q.2)
    (hFnn : ∀ N x Y, 0 ≤ F N x Y)
    (hint : ∀ N Y, Integrable (fun x => (F N x Y) ^ 2) (piGauss (dN N)))
    {K : ℝ}
    (hbnd : ∀ᶠ N in atTop, ∀ Y, ∫ x, (F N x Y) ^ 2 ∂(piGauss (dN N)) ≤ K / dN N)
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω => F N (ZZ N ω).1 (ZZ N ω).2) 0 := by
  intro ε hε
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds
    (tendsto_ofReal_div_atTop (K := K) hd ε)
    (Filter.Eventually.of_forall fun _ => by simp) ?_
  filter_upwards [hbnd] with N hN
  have hset : {ω | ε ≤ |F N (ZZ N ω).1 (ZZ N ω).2 - 0|}
      = {ω | ε ≤ F N (ZZ N ω).1 (ZZ N ω).2} := by
    ext ω
    simp [abs_of_nonneg (hFnn N _ _)]
  have hmeasSet : MeasurableSet {q : NoiseSpace (pN N) (dN N) | ε ≤ F N q.1 q.2} :=
    measurableSet_le measurable_const (hFm N)
  rw [hset, (hZZ N).measure_eq hmeasSet, noiseLaw, Measure.prod_apply_symm hmeasSet]
  have hsec : ∀ Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ,
      (piGauss (dN N)) ((fun x => (x, Y)) ⁻¹' {q | ε ≤ F N q.1 q.2})
        ≤ ENNReal.ofReal (K / dN N / ε ^ 2) := by
    intro Y
    exact meas_ge_le_of_integral_sq _ ((hFm N).comp measurable_prodMk_right)
      (hint N Y) hε (hN Y)
  calc ∫⁻ Y, (piGauss (dN N)) ((fun x => (x, Y)) ⁻¹' {q | ε ≤ F N q.1 q.2})
        ∂(gaussianMatrix (pN N) (dN N))
      ≤ ∫⁻ _Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ, ENNReal.ofReal (K / dN N / ε ^ 2)
        ∂(gaussianMatrix (pN N) (dN N)) := lintegral_mono hsec
    _ = ENNReal.ofReal (K / dN N / ε ^ 2) := by simp

end Transfer

/-! ### The `v`-forms by symmetry -/

section Symmetry

/-- The direction of the Gaussian vector. -/
noncomputable def uOf (x : Fin d → ℝ) : Fin d → ℝ :=
  (Real.sqrt (gOf x ⬝ᵥ gOf x))⁻¹ • gOf x

theorem dotProduct_gOf (x : Fin d → ℝ) : gOf x ⬝ᵥ gOf x = ((d : ℝ))⁻¹ * (x ⬝ᵥ x) := by
  rw [dotProduct, dotProduct, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  change ((Real.sqrt d)⁻¹ * x a) * ((Real.sqrt d)⁻¹ * x a) = _
  rw [← sqrt_inv_mul_sqrt_inv d]
  ring

theorem dotProduct_uOf {x : Fin d → ℝ} (hs : gOf x ⬝ᵥ gOf x ≠ 0) : uOf x ⬝ᵥ uOf x = 1 := by
  have hnn : 0 ≤ gOf x ⬝ᵥ gOf x := by
    rw [dotProduct]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rw [uOf, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc,
    ← Real.sqrt_inv, ← Real.sqrt_mul (le_of_lt (inv_pos.2 (lt_of_le_of_ne hnn (Ne.symm hs)))),
    Real.sqrt_mul_self (le_of_lt (inv_pos.2 (lt_of_le_of_ne hnn (Ne.symm hs))))]
  exact inv_mul_cancel₀ hs

theorem cvec_smul (r : ℝ) (x : Fin d → ℝ) : R4C.cvec (r • x) = (r : ℂ) • R4C.cvec x := by
  funext a
  simp [R4C.cvec]

theorem bil_smul (B : Matrix (Fin d) (Fin d) ℂ) (r : ℝ) (x y : Fin d → ℝ) :
    bil B (r • x) (r • y) = ((r ^ 2 : ℝ) : ℂ) * bil B x y := by
  rw [bil, bil, cvec_smul, cvec_smul, Matrix.mulVec_smul, dotProduct_smul, smul_dotProduct,
    smul_eq_mul, smul_eq_mul]
  push_cast
  ring

theorem qformC_smul (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (r : ℝ) (y : Fin d → ℝ) :
    R4C.qformC W z (r • y) = ((r ^ 2 : ℝ) : ℂ) * R4C.qformC W z y :=
  bil_smul _ r y y

theorem qform2C_smul (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (r : ℝ) (y : Fin d → ℝ) :
    R4C.qform2C W z (r • y) = ((r ^ 2 : ℝ) : ℂ) * R4C.qform2C W z y :=
  bil_smul _ r y y

theorem qformC_uOf (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x : Fin d → ℝ) :
    R4C.qformC W z (uOf x) = (((gOf x ⬝ᵥ gOf x)⁻¹ : ℝ) : ℂ) * R4C.qformC W z (gOf x) := by
  have hnn : 0 ≤ gOf x ⬝ᵥ gOf x := by
    rw [dotProduct]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rw [uOf, qformC_smul]
  congr 2
  rw [← Real.sqrt_inv, Real.sq_sqrt (inv_nonneg.2 hnn)]

theorem qform2C_uOf (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x : Fin d → ℝ) :
    R4C.qform2C W z (uOf x) = (((gOf x ⬝ᵥ gOf x)⁻¹ : ℝ) : ℂ) * R4C.qform2C W z (gOf x) := by
  have hnn : 0 ≤ gOf x ⬝ᵥ gOf x := by
    rw [dotProduct]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rw [uOf, qform2C_smul]
  congr 2
  rw [← Real.sqrt_inv, Real.sq_sqrt (inv_nonneg.2 hnn)]

/-- **Exchangeability of the unit vector.** -/
theorem meas_form_eq {v u : Fin d → ℝ} (hv : v ⬝ᵥ v = 1) (hu : u ⬝ᵥ u = 1)
    (Ψ : Matrix (Fin d) (Fin d) ℝ → (Fin d → ℝ) → ℝ)
    (hΨ : ∀ O : Matrix (Fin d) (Fin d) ℝ, Oᵀ * O = 1 → ∀ (Y : Matrix (Fin p) (Fin d) ℝ) y,
      Ψ (W0 (Y * O)) y = Ψ (W0 Y) (O *ᵥ y))
    (hmeas : ∀ y, Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => Ψ (W0 Y) y) (ε : ℝ) :
    (gaussianMatrix p d) {Y | ε ≤ Ψ (W0 Y) u} = (gaussianMatrix p d) {Y | ε ≤ Ψ (W0 Y) v} := by
  rcases eq_or_ne u v with rfl | hne
  · rfl
  have hne0 : (u - v) ⬝ᵥ (u - v) ≠ 0 := fun h =>
    hne (sub_eq_zero.mp (dotProduct_self_eq_zero.mp h))
  have hO : (householder (u - v))ᵀ * householder (u - v) = 1 := householder_orth hne0
  have hOu : householder (u - v) *ᵥ u = v := householder_sub_apply hu hv hne
  have hS : MeasurableSet {Y : Matrix (Fin p) (Fin d) ℝ | ε ≤ Ψ (W0 Y) u} :=
    measurableSet_le measurable_const (hmeas u)
  have hpre := (measurePreserving_mul_right hO p).measure_preimage hS.nullMeasurableSet
  rw [← hpre]
  congr 1
  ext Y
  simp only [Set.mem_preimage, Set.mem_ofPred_eq, hΨ _ hO Y u, hOu]

/-- The set where the Gaussian vector vanishes is null. -/
theorem piGauss_dotProduct_eq_zero (hd : 0 < d) :
    (piGauss d) {x : Fin d → ℝ | gOf x ⬝ᵥ gOf x = 0} = 0 := by
  set j : Fin d := ⟨0, hd⟩ with hjdef
  have hsub : {x : Fin d → ℝ | gOf x ⬝ᵥ gOf x = 0} ⊆ (fun x : Fin d → ℝ => x j) ⁻¹' {0} := by
    intro x hx
    have hx0 : x ⬝ᵥ x = 0 := by
      rcases eq_or_ne ((d : ℝ))⁻¹ 0 with h | h
      · exact absurd (Nat.cast_ne_zero.mpr hd.ne' : ((d : ℝ)) ≠ 0) (by simpa using h)
      · have hx' : ((d : ℝ))⁻¹ * (x ⬝ᵥ x) = 0 := by
          rw [← dotProduct_gOf]
          exact hx
        exact (mul_eq_zero.mp hx').resolve_left h
    have : x = 0 := dotProduct_self_eq_zero.mp hx0
    simp [this]
  refine measure_mono_null hsub ?_
  have hmp := (measurePreserving_coord (d := d) j).measure_preimage
    (measurableSet_singleton (0 : ℝ)).nullMeasurableSet
  rw [hmp]
  exact gaussianReal_absolutelyContinuous 0 (by simp) (by simp)

end Symmetry

/-! ### Assembly: the six limits of (H2) -/

section Assembly

/-! #### Measurability of the six forms on the product space -/

theorem measurable_gOf_coord (i : Fin d) :
    Measurable fun q : NoiseSpace p d => gOf q.1 i :=
  measurable_const.mul ((measurable_pi_apply i).comp measurable_fst)

theorem measurable_dotProduct_gOf :
    Measurable fun q : NoiseSpace p d => gOf q.1 ⬝ᵥ gOf q.1 := by
  simp only [dotProduct]
  exact Finset.measurable_sum _ fun a _ =>
    (measurable_gOf_coord a).mul (measurable_gOf_coord a)

theorem measurable_uOf_coord (i : Fin d) :
    Measurable fun q : NoiseSpace p d => uOf q.1 i :=
  ((Real.continuous_sqrt.measurable.comp measurable_dotProduct_gOf).inv).mul
    (measurable_gOf_coord i)

theorem measurable_resolvC_snd (z : ℂ) (i j : Fin d) :
    Measurable fun q : NoiseSpace p d => R4C.resolvC (W0 q.2) z i j :=
  (measurable_resolvC_entry z i j).comp measurable_snd

theorem measurable_resolvC_sq_snd (z : ℂ) (i j : Fin d) :
    Measurable fun q : NoiseSpace p d =>
      (R4C.resolvC (W0 q.2) z * R4C.resolvC (W0 q.2) z) i j :=
  (measurable_resolvC_sq_entry z i j).comp measurable_snd

theorem measurable_qformC_gOf (z : ℂ) :
    Measurable fun q : NoiseSpace p d => R4C.qformC (W0 q.2) z (gOf q.1) :=
  measurable_bil (measurable_resolvC_snd z) measurable_gOf_coord measurable_gOf_coord

theorem measurable_qform2C_gOf (z : ℂ) :
    Measurable fun q : NoiseSpace p d => R4C.qform2C (W0 q.2) z (gOf q.1) :=
  measurable_bil (measurable_resolvC_sq_snd z) measurable_gOf_coord measurable_gOf_coord

theorem measurable_qformC_uOf (z : ℂ) :
    Measurable fun q : NoiseSpace p d => R4C.qformC (W0 q.2) z (uOf q.1) :=
  measurable_bil (measurable_resolvC_snd z) measurable_uOf_coord measurable_uOf_coord

theorem measurable_qform2C_uOf (z : ℂ) :
    Measurable fun q : NoiseSpace p d => R4C.qform2C (W0 q.2) z (uOf q.1) :=
  measurable_bil (measurable_resolvC_sq_snd z) measurable_uOf_coord measurable_uOf_coord

theorem measurable_cformC_gOf (z : ℂ) (v : Fin d → ℝ) :
    Measurable fun q : NoiseSpace p d => R4C.cformC (W0 q.2) z v (gOf q.1) :=
  measurable_bil (measurable_resolvC_snd z) (fun _ => measurable_const) measurable_gOf_coord

theorem measurable_cform2C_gOf (z : ℂ) (v : Fin d → ℝ) :
    Measurable fun q : NoiseSpace p d => R4C.cform2C (W0 q.2) z v (gOf q.1) :=
  measurable_bil (measurable_resolvC_sq_snd z) (fun _ => measurable_const) measurable_gOf_coord

theorem measurable_qformC_fixed (z : ℂ) (v : Fin d → ℝ) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.qformC (W0 Y) z v :=
  measurable_bil (fun i j => measurable_resolvC_entry z i j) (fun _ => measurable_const)
    fun _ => measurable_const

theorem measurable_qform2C_fixed (z : ℂ) (v : Fin d → ℝ) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.qform2C (W0 Y) z v :=
  measurable_bil (fun i j => measurable_resolvC_sq_entry z i j) (fun _ => measurable_const)
    fun _ => measurable_const

/-! #### Convergence in probability, general steps -/

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN dN : ℕ → ℕ} {c : ℝ}

/-- Two sequences with the same law at every level have the same limits in probability. -/
theorem tendstoInProb_of_meas_eq {f g : ∀ N, Ω N → ℝ} {a : ℝ}
    (h : ∀ N, ∀ ε > 0, μ N {ω | ε ≤ |f N ω - a|} = μ N {ω | ε ≤ |g N ω - a|})
    (hg : TendstoInProb μ g a) : TendstoInProb μ f a := by
  intro ε hε
  simp only [h _ ε hε]
  exact hg ε hε

/-- Centering: the gap to a random center plus the gap of the center. -/
theorem tendstoInProb_norm_sub_trans {A B : ∀ N, Ω N → ℂ} {m : ℂ}
    (h1 : TendstoInProb μ (fun N ω => ‖A N ω - B N ω‖) 0)
    (h2 : TendstoInProb μ (fun N ω => ‖B N ω - m‖) 0) :
    TendstoInProb μ (fun N ω => ‖A N ω - m‖) 0 := by
  refine TendstoInProb.of_le (g := fun N ω => ‖A N ω - B N ω‖ + ‖B N ω - m‖) ?_ ?_
  · intro N
    filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    calc ‖A N ω - m‖ = ‖(A N ω - B N ω) + (B N ω - m)‖ := by ring_nf
      _ ≤ ‖A N ω - B N ω‖ + ‖B N ω - m‖ := norm_add_le _ _
  · simpa using h1.add h2

theorem dotProduct_ofLp_self {D : ℕ} {v : EuclideanSpace ℝ (Fin D)} (hv : ‖v‖ = 1) :
    WithLp.ofLp v ⬝ᵥ WithLp.ofLp v = 1 := by
  have h := inner_euclidean_eq_dotProduct v v
  rw [real_inner_self_eq_norm_sq, hv] at h
  simpa using h.symm

/-! #### The four forms in `g` -/

variable (ZZ : ∀ N, Ω N → NoiseSpace (pN N) (dN N))

theorem tendstoInProb_qformC_gOf_sub (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1)
      - R4C.stieltjesC (W0 (ZZ N ω).2) z‖) 0 :=
  tendstoInProb_of_integral_sq_le ZZ hZZ
    (fun _ x Y => ‖R4C.qformC (W0 Y) z (gOf x) - R4C.stieltjesC (W0 Y) z‖)
    (fun _ => ((measurable_qformC_gOf z).sub
      ((measurable_stieltjesC z).comp measurable_snd)).norm)
    (fun _ _ _ => norm_nonneg _)
    (fun _ Y => integrable_normSq_qformC (isHermitian_W0 Y) hz)
    (K := varSq / z.im ^ 2)
    (Filter.Eventually.of_forall fun _ Y => integral_normSq_qformC_le (isHermitian_W0 Y) hz) hd

theorem tendstoInProb_qform2C_gOf_sub (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1)
      - R4C.stieltjes2C (W0 (ZZ N ω).2) z‖) 0 :=
  tendstoInProb_of_integral_sq_le ZZ hZZ
    (fun _ x Y => ‖R4C.qform2C (W0 Y) z (gOf x) - R4C.stieltjes2C (W0 Y) z‖)
    (fun _ => ((measurable_qform2C_gOf z).sub
      ((measurable_stieltjes2C z).comp measurable_snd)).norm)
    (fun _ _ _ => norm_nonneg _)
    (fun _ Y => integrable_normSq_qform2C (isHermitian_W0 Y) hz)
    (K := varSq / z.im ^ 4)
    (Filter.Eventually.of_forall fun _ Y => integral_normSq_qform2C_le (isHermitian_W0 Y) hz) hd

theorem tendstoInProb_cformC_gOf (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (v : (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))) (hv : ∀ N, ‖v N‖ = 1)
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) (gOf (ZZ N ω).1)‖) 0 := by
  refine tendstoInProb_of_integral_sq_le ZZ hZZ
    (fun N x Y => ‖R4C.cformC (W0 Y) z (WithLp.ofLp (v N)) (gOf x)‖)
    (fun N => (measurable_cformC_gOf z (WithLp.ofLp (v N))).norm)
    (fun _ _ _ => norm_nonneg _)
    (fun N Y => integrable_normSq_cformC (isHermitian_W0 Y) hz _)
    (K := 1 / z.im ^ 2) (Filter.Eventually.of_forall fun N Y => ?_) hd
  have h := integral_normSq_cformC_le (isHermitian_W0 Y) hz (WithLp.ofLp (v N))
  rwa [dotProduct_ofLp_self (hv N)] at h

theorem tendstoInProb_cform2C_gOf (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (v : (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))) (hv : ∀ N, ‖v N‖ = 1)
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) (gOf (ZZ N ω).1)‖) 0 := by
  refine tendstoInProb_of_integral_sq_le ZZ hZZ
    (fun N x Y => ‖R4C.cform2C (W0 Y) z (WithLp.ofLp (v N)) (gOf x)‖)
    (fun N => (measurable_cform2C_gOf z (WithLp.ofLp (v N))).norm)
    (fun _ _ _ => norm_nonneg _)
    (fun N Y => integrable_normSq_cform2C (isHermitian_W0 Y) hz _)
    (K := 1 / z.im ^ 4) (Filter.Eventually.of_forall fun N Y => ?_) hd
  have h := integral_normSq_cform2C_le (isHermitian_W0 Y) hz (WithLp.ofLp (v N))
  rwa [dotProduct_ofLp_self (hv N)] at h

theorem tendstoInProb_dotProduct_gOf
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hd : Tendsto dN atTop atTop) :
    TendstoInProb μ (fun N ω => gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1) 1 := by
  have hbase : TendstoInProb μ
      (fun N ω => |gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1 - 1|) 0 := by
    refine tendstoInProb_of_integral_sq_le ZZ hZZ
      (fun _ x _ => |gOf x ⬝ᵥ gOf x - 1|)
      (fun _ => (measurable_dotProduct_gOf.sub measurable_const).abs)
      (fun _ _ _ => abs_nonneg _) (fun N Y => ?_) (K := varSq) ?_ hd
    · simpa [sq_abs] using integrable_sq_dotProduct_gOf (d := dN N)
    · filter_upwards [hd.eventually_gt_atTop 0] with N hN Y
      simpa [sq_abs] using integral_sq_dotProduct_gOf_le (d := dN N) hN
  refine TendstoInProb.of_le
    (g := fun N ω => |gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1 - 1|) (fun N => ?_) hbase
  filter_upwards with ω
  exact le_rfl

/-! #### The two forms at `v`, by symmetry -/

theorem qformC_W0_mul (hz : z.im ≠ 0) {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1)
    (Y : Matrix (Fin p) (Fin d) ℝ) (y : Fin d → ℝ) :
    R4C.qformC (W0 (Y * O)) z y = R4C.qformC (W0 Y) z (O *ᵥ y) := by
  rw [W0_mul, R4C.qformC, R4C.qformC, cformC_conj hO (isHermitian_W0 Y) hz]

theorem qform2C_W0_mul (hz : z.im ≠ 0) {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1)
    (Y : Matrix (Fin p) (Fin d) ℝ) (y : Fin d → ℝ) :
    R4C.qform2C (W0 (Y * O)) z y = R4C.qform2C (W0 Y) z (O *ᵥ y) := by
  rw [W0_mul, R4C.qform2C, R4C.qform2C, cform2C_conj hO (isHermitian_W0 Y) hz]

theorem pos_of_dotProduct_one {D : ℕ} {v : Fin D → ℝ} (hv : v ⬝ᵥ v = 1) : 0 < D := by
  rcases Nat.eq_zero_or_pos D with rfl | h
  · simp [dotProduct] at hv
  · exact h

/-- **The measure identity of the symmetry route.** -/
theorem noiseLaw_meas_eq (hd : 0 < d) {v : Fin d → ℝ} (hv : v ⬝ᵥ v = 1)
    (Ψ : Matrix (Fin d) (Fin d) ℝ → (Fin d → ℝ) → ℝ)
    (hΨ : ∀ O : Matrix (Fin d) (Fin d) ℝ, Oᵀ * O = 1 → ∀ (Y : Matrix (Fin p) (Fin d) ℝ) y,
      Ψ (W0 (Y * O)) y = Ψ (W0 Y) (O *ᵥ y))
    (hmeasv : ∀ y, Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => Ψ (W0 Y) y)
    (hmeasq : Measurable fun q : NoiseSpace p d => Ψ (W0 q.2) (uOf q.1)) (ε : ℝ) :
    noiseLaw p d {q | ε ≤ Ψ (W0 q.2) v} = noiseLaw p d {q | ε ≤ Ψ (W0 q.2) (uOf q.1)} := by
  have hset1 : MeasurableSet {q : NoiseSpace p d | ε ≤ Ψ (W0 q.2) v} :=
    measurableSet_le measurable_const ((hmeasv v).comp measurable_snd)
  have hset2 : MeasurableSet {q : NoiseSpace p d | ε ≤ Ψ (W0 q.2) (uOf q.1)} :=
    measurableSet_le measurable_const hmeasq
  rw [noiseLaw, Measure.prod_apply hset1, Measure.prod_apply hset2]
  refine lintegral_congr_ae ?_
  have hae : ∀ᵐ x ∂(piGauss d), gOf x ⬝ᵥ gOf x ≠ 0 := by
    rw [ae_iff]
    simpa using piGauss_dotProduct_eq_zero (d := d) hd
  filter_upwards [hae] with x hx
  have hu : uOf x ⬝ᵥ uOf x = 1 := dotProduct_uOf hx
  exact (meas_form_eq hv hu Ψ hΨ hmeasv ε).symm

theorem norm_ofReal_mul_sub_le (s : ℝ) (Q m : ℂ) :
    ‖((s : ℝ) : ℂ) * Q - m‖ ≤ |s| * ‖Q - m‖ + |s - 1| * ‖m‖ := by
  have h : ((s : ℝ) : ℂ) * Q - m = ((s : ℝ) : ℂ) * (Q - m) + (((s - 1 : ℝ)) : ℂ) * m := by
    push_cast; ring
  rw [h]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs]

theorem tendstoInProb_qformC_uOf (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hd : Tendsto dN atTop atTop)
    (hR1a : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (W0 (ZZ N ω).2) z - MP.mC c z‖) 0) :
    TendstoInProb μ
      (fun N ω => ‖R4C.qformC (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mC c z‖) 0 := by
  have hgg : TendstoInProb μ
      (fun N ω => ‖R4C.qformC (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mC c z‖) 0 :=
    tendstoInProb_norm_sub_trans (tendstoInProb_qformC_gOf_sub ZZ hz hZZ hd) hR1a
  have hs := tendstoInProb_dotProduct_gOf ZZ hZZ hd
  have hinv : TendstoInProb μ (fun N ω => (gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹) 1 := by
    simpa using hs.inv one_ne_zero
  have habs : TendstoInProb μ (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹|) 1 := by
    simpa using hinv.comp_continuous (φ := fun t : ℝ => |t|) continuous_abs.continuousAt
  have h1 : TendstoInProb μ
      (fun N ω => (gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1) 0 := by
    simpa using hinv.sub (TendstoInProb.const μ 1)
  have h2 : TendstoInProb μ
      (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1| * ‖MP.mC c z‖) 0 := by
    have h3 : TendstoInProb μ
        (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1|) 0 := by
      simpa using h1.comp_continuous (φ := fun t : ℝ => |t|) continuous_abs.continuousAt
    simpa using h3.mul_const ‖MP.mC c z‖
  refine TendstoInProb.of_le (g := fun N ω =>
    |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹| *
        ‖R4C.qformC (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mC c z‖
      + |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1| * ‖MP.mC c z‖) (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _), qformC_uOf]
    exact norm_ofReal_mul_sub_le _ _ _
  · simpa using (habs.mul hgg).add h2

theorem tendstoInProb_qform2C_uOf (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hd : Tendsto dN atTop atTop)
    (hR1b : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0) :
    TendstoInProb μ
      (fun N ω => ‖R4C.qform2C (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mCDeriv c z‖) 0 := by
  have hgg : TendstoInProb μ
      (fun N ω => ‖R4C.qform2C (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mCDeriv c z‖) 0 :=
    tendstoInProb_norm_sub_trans (tendstoInProb_qform2C_gOf_sub ZZ hz hZZ hd) hR1b
  have hs := tendstoInProb_dotProduct_gOf ZZ hZZ hd
  have hinv : TendstoInProb μ (fun N ω => (gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹) 1 := by
    simpa using hs.inv one_ne_zero
  have habs : TendstoInProb μ (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹|) 1 := by
    simpa using hinv.comp_continuous (φ := fun t : ℝ => |t|) continuous_abs.continuousAt
  have h1 : TendstoInProb μ
      (fun N ω => (gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1) 0 := by
    simpa using hinv.sub (TendstoInProb.const μ 1)
  have h2 : TendstoInProb μ
      (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1| * ‖MP.mCDeriv c z‖) 0 := by
    have h3 : TendstoInProb μ
        (fun N ω => |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1|) 0 := by
      simpa using h1.comp_continuous (φ := fun t : ℝ => |t|) continuous_abs.continuousAt
    simpa using h3.mul_const ‖MP.mCDeriv c z‖
  refine TendstoInProb.of_le (g := fun N ω =>
    |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹| *
        ‖R4C.qform2C (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mCDeriv c z‖
      + |(gOf (ZZ N ω).1 ⬝ᵥ gOf (ZZ N ω).1)⁻¹ - 1| * ‖MP.mCDeriv c z‖) (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _), qform2C_uOf]
    exact norm_ofReal_mul_sub_le _ _ _
  · simpa using (habs.mul hgg).add h2

theorem tendstoInProb_qformC_v (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (v : (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))) (hv : ∀ N, ‖v N‖ = 1)
    (hd : Tendsto dN atTop atTop)
    (hR1a : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (W0 (ZZ N ω).2) z - MP.mC c z‖) 0) :
    TendstoInProb μ
      (fun N ω => ‖R4C.qformC (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mC c z‖) 0 := by
  refine tendstoInProb_of_meas_eq (g := fun N ω =>
    ‖R4C.qformC (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mC c z‖) (fun N ε hε => ?_)
    (tendstoInProb_qformC_uOf ZZ hz hZZ hd hR1a)
  set Ψ : Matrix (Fin (dN N)) (Fin (dN N)) ℝ → (Fin (dN N) → ℝ) → ℝ :=
    fun W y => ‖R4C.qformC W z y - MP.mC c z‖ with hΨdef
  have hVunit : WithLp.ofLp (v N) ⬝ᵥ WithLp.ofLp (v N) = 1 := dotProduct_ofLp_self (hv N)
  have hdpos : 0 < dN N := pos_of_dotProduct_one hVunit
  have hconj : ∀ O : Matrix (Fin (dN N)) (Fin (dN N)) ℝ, Oᵀ * O = 1 →
      ∀ (Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ) y,
      Ψ (W0 (Y * O)) y = Ψ (W0 Y) (O *ᵥ y) := by
    intro O hO Y y
    simp only [hΨdef, qformC_W0_mul hz.ne' hO]
  have hmv : ∀ y, Measurable fun Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ => Ψ (W0 Y) y :=
    fun y => ((measurable_qformC_fixed z y).sub measurable_const).norm
  have hmq : Measurable fun q : NoiseSpace (pN N) (dN N) => Ψ (W0 q.2) (uOf q.1) :=
    ((measurable_qformC_uOf z).sub measurable_const).norm
  have hset1 : {ω | ε ≤ |‖R4C.qformC (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mC c z‖ - 0|}
      = {ω | ε ≤ Ψ (W0 (ZZ N ω).2) (WithLp.ofLp (v N))} := by
    ext ω
    simp [hΨdef, abs_of_nonneg (norm_nonneg _)]
  have hset2 : {ω | ε ≤ |‖R4C.qformC (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mC c z‖ - 0|}
      = {ω | ε ≤ Ψ (W0 (ZZ N ω).2) (uOf (ZZ N ω).1)} := by
    ext ω
    simp [hΨdef, abs_of_nonneg (norm_nonneg _)]
  have e1 := (hZZ N).measure_eq
    (p := fun q : NoiseSpace (pN N) (dN N) => ε ≤ Ψ (W0 q.2) (WithLp.ofLp (v N)))
    (measurableSet_le measurable_const ((hmv _).comp measurable_snd))
  have e2 := (hZZ N).measure_eq
    (p := fun q : NoiseSpace (pN N) (dN N) => ε ≤ Ψ (W0 q.2) (uOf q.1))
    (measurableSet_le measurable_const hmq)
  rw [hset1, hset2, e1, e2]
  exact noiseLaw_meas_eq hdpos hVunit Ψ hconj hmv hmq ε

theorem tendstoInProb_qform2C_v (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (v : (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))) (hv : ∀ N, ‖v N‖ = 1)
    (hd : Tendsto dN atTop atTop)
    (hR1b : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0) :
    TendstoInProb μ
      (fun N ω => ‖R4C.qform2C (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mCDeriv c z‖) 0 := by
  refine tendstoInProb_of_meas_eq (g := fun N ω =>
    ‖R4C.qform2C (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mCDeriv c z‖) (fun N ε hε => ?_)
    (tendstoInProb_qform2C_uOf ZZ hz hZZ hd hR1b)
  set Ψ : Matrix (Fin (dN N)) (Fin (dN N)) ℝ → (Fin (dN N) → ℝ) → ℝ :=
    fun W y => ‖R4C.qform2C W z y - MP.mCDeriv c z‖ with hΨdef
  have hVunit : WithLp.ofLp (v N) ⬝ᵥ WithLp.ofLp (v N) = 1 := dotProduct_ofLp_self (hv N)
  have hdpos : 0 < dN N := pos_of_dotProduct_one hVunit
  have hconj : ∀ O : Matrix (Fin (dN N)) (Fin (dN N)) ℝ, Oᵀ * O = 1 →
      ∀ (Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ) y,
      Ψ (W0 (Y * O)) y = Ψ (W0 Y) (O *ᵥ y) := by
    intro O hO Y y
    simp only [hΨdef, qform2C_W0_mul hz.ne' hO]
  have hmv : ∀ y, Measurable fun Y : Matrix (Fin (pN N)) (Fin (dN N)) ℝ => Ψ (W0 Y) y :=
    fun y => ((measurable_qform2C_fixed z y).sub measurable_const).norm
  have hmq : Measurable fun q : NoiseSpace (pN N) (dN N) => Ψ (W0 q.2) (uOf q.1) :=
    ((measurable_qform2C_uOf z).sub measurable_const).norm
  have hset1 :
      {ω | ε ≤ |‖R4C.qform2C (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mCDeriv c z‖ - 0|}
      = {ω | ε ≤ Ψ (W0 (ZZ N ω).2) (WithLp.ofLp (v N))} := by
    ext ω
    simp [hΨdef, abs_of_nonneg (norm_nonneg _)]
  have hset2 :
      {ω | ε ≤ |‖R4C.qform2C (W0 (ZZ N ω).2) z (uOf (ZZ N ω).1) - MP.mCDeriv c z‖ - 0|}
      = {ω | ε ≤ Ψ (W0 (ZZ N ω).2) (uOf (ZZ N ω).1)} := by
    ext ω
    simp [hΨdef, abs_of_nonneg (norm_nonneg _)]
  have e1 := (hZZ N).measure_eq
    (p := fun q : NoiseSpace (pN N) (dN N) => ε ≤ Ψ (W0 q.2) (WithLp.ofLp (v N)))
    (measurableSet_le measurable_const ((hmv _).comp measurable_snd))
  have e2 := (hZZ N).measure_eq
    (p := fun q : NoiseSpace (pN N) (dN N) => ε ≤ Ψ (W0 q.2) (uOf q.1))
    (measurableSet_le measurable_const hmq)
  rw [hset1, hset2, e1, e2]
  exact noiseLaw_meas_eq hdpos hVunit Ψ hconj hmv hmq ε

end Assembly

/-! ### R2d: the six limits, as item T consumes them -/

section Main

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN dN : ℕ → ℕ} {c : ℝ}

set_option linter.unusedVariables false in
/-- **R2d.** The six forms item T consumes, each a plain `TendstoInProb` in the norm form of
choice 28. `hR1a` and `hR1b` are R1a and R1b (mismatch M4): R2 proves the gap to the random
`s_N(z)` and R1 the gap from `s_N(z)` to the scalar limit. No supremum and no `orth N`
appears: `delocUniform` is item Sym (decision D6). The hypotheses `hc`, `hp`, `hcN` and
`IsProbabilityMeasure` are not used by this proof; they stay for the interface of
`notes/archive/rmt_R2.md`. -/
theorem tendsto_forms [∀ N, IsProbabilityMeasure (μ N)]
    (hc : 0 < c) (hz : 0 < z.im)
    (hd : Tendsto dN atTop atTop) (hp : ∀ N, 0 < pN N)
    (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c))
    (v : (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))) (hv : ∀ N, ‖v N‖ = 1)
    (ZZ : ∀ N, Ω N → NoiseSpace (pN N) (dN N))
    (hZZ : ∀ N, HasLaw (ZZ N) (noiseLaw (pN N) (dN N)) (μ N))
    (hR1a : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (W0 (ZZ N ω).2) z - MP.mC c z‖) 0)
    (hR1b : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0) :
    -- 1. `vᵀ G v → mᶜ`
    TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mC c z‖) 0 ∧
    -- 2. `gᵀ G g → mᶜ`
    TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mC c z‖) 0 ∧
    -- 3. `vᵀ G g → 0`
    TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) (gOf (ZZ N ω).1)‖) 0 ∧
    -- 4. `vᵀ G² v → mᶜ'`
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) - MP.mCDeriv c z‖) 0 ∧
    -- 5. `gᵀ G² g → mᶜ'`
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (W0 (ZZ N ω).2) z (gOf (ZZ N ω).1) - MP.mCDeriv c z‖) 0 ∧
    -- 6. `vᵀ G² g → 0`
    TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W0 (ZZ N ω).2) z (WithLp.ofLp (v N)) (gOf (ZZ N ω).1)‖) 0 :=
  ⟨tendstoInProb_qformC_v ZZ hz hZZ v hv hd hR1a,
    tendstoInProb_norm_sub_trans (tendstoInProb_qformC_gOf_sub ZZ hz hZZ hd) hR1a,
    tendstoInProb_cformC_gOf ZZ hz hZZ v hv hd,
    tendstoInProb_qform2C_v ZZ hz hZZ v hv hd hR1b,
    tendstoInProb_norm_sub_trans (tendstoInProb_qform2C_gOf_sub ZZ hz hZZ hd) hR1b,
    tendstoInProb_cform2C_gOf ZZ hz hZZ v hv hd⟩

end Main

/-! ### R2a: the mean of the resolvent is a multiple of the identity

The assembly above does not use this section: step 3 centers the `g`-forms at the random
`s_N(z)` exactly, so no `E[s_N]` appears. R2a is the statement of `notes/archive/rmt_R2.md`, proved
from the same conjugation lemma with signed permutations. -/

section MeanResolvent

variable {p d : ℕ} {z : ℂ}

/-- The signed permutation matrix `D`. -/
noncomputable def sgnPermMat (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ) :
    Matrix (Fin d) (Fin d) ℝ :=
  Matrix.of fun l j => if l = σ j then ε j else 0

/-- `Y ↦ Y D` in coordinates. -/
noncomputable def sgnPerm (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ)
    (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin p) (Fin d) ℝ :=
  Matrix.of fun k j => ε j * Y k (σ j)

theorem sgnPerm_eq_mul (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ)
    (Y : Matrix (Fin p) (Fin d) ℝ) : sgnPerm σ ε Y = Y * sgnPermMat σ ε := by
  ext k j
  change ε j * Y k (σ j) = ∑ l, Y k l * (if l = σ j then ε j else 0)
  rw [Finset.sum_eq_single (σ j) (fun l _ hl => by simp [hl]) (by simp)]
  simp [mul_comm]

theorem sgnPermMat_orth {σ : Equiv.Perm (Fin d)} {ε : Fin d → ℝ}
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) : (sgnPermMat σ ε)ᵀ * sgnPermMat σ ε = 1 := by
  ext i j
  rw [Matrix.mul_apply]
  rcases eq_or_ne i j with rfl | hij
  · have hone : ∀ l : Fin d,
        (sgnPermMat σ ε)ᵀ i l * sgnPermMat σ ε l i
          = if l = σ i then ε i * ε i else 0 := by
      intro l
      change (if l = σ i then ε i else 0) * (if l = σ i then ε i else 0) = _
      by_cases h : l = σ i <;> simp [h]
    simp only [hone]
    rw [Finset.sum_ite_eq' Finset.univ (σ i) fun _ => ε i * ε i]
    have : ε i * ε i = 1 := by rcases hε i with h | h <;> simp [h]
    simp [this]
  · have hzero : ∀ l : Fin d, (sgnPermMat σ ε)ᵀ i l * sgnPermMat σ ε l j = 0 := by
      intro l
      change (if l = σ i then ε i else 0) * (if l = σ j then ε j else 0) = 0
      by_cases h1 : l = σ i
      · have h2 : l ≠ σ j := by
          rw [h1]
          exact fun h => hij (σ.injective h)
        simp [h2]
      · simp [h1]
    rw [Matrix.one_apply_ne hij]
    exact Finset.sum_eq_zero fun l _ => hzero l

theorem measurePreserving_sgnPerm (σ : Equiv.Perm (Fin d)) {ε : Fin d → ℝ}
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) :
    MeasurePreserving (sgnPerm σ ε : Matrix (Fin p) (Fin d) ℝ → _)
      (gaussianMatrix p d) (gaussianMatrix p d) := by
  have h := measurePreserving_mul_right (sgnPermMat_orth (σ := σ) hε) p
  simpa only [← sgnPerm_eq_mul] using h

theorem gaussianMatrix_map_sgnPerm (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ)
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) :
    (gaussianMatrix p d).map (sgnPerm σ ε) = gaussianMatrix p d :=
  (measurePreserving_sgnPerm σ hε).map_eq

theorem resolvC_W0_sgnPerm (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ)
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) (hz : 0 < z.im) (Y : Matrix (Fin p) (Fin d) ℝ)
    (i j : Fin d) :
    R4C.resolvC (W0 (sgnPerm σ ε Y)) z i j
      = (ε i : ℂ) * (ε j : ℂ) * R4C.resolvC (W0 Y) z (σ i) (σ j) := by
  have hD := sgnPermMat_orth (σ := σ) hε
  rw [sgnPerm_eq_mul, W0_mul, resolvC_conj hD (isHermitian_W0 Y) hz.ne']
  rw [Matrix.mul_apply]
  have hinner : ∀ l : Fin d,
      ((R4C.cmat (sgnPermMat σ ε))ᵀ * R4C.resolvC (W0 Y) z) i l *
          R4C.cmat (sgnPermMat σ ε) l j
        = if l = σ j then
            (ε i : ℂ) * R4C.resolvC (W0 Y) z (σ i) l * (ε j : ℂ) else 0 := by
    intro l
    rw [Matrix.mul_apply]
    have houter : ∀ k : Fin d,
        (R4C.cmat (sgnPermMat σ ε))ᵀ i k * R4C.resolvC (W0 Y) z k l
          = if k = σ i then (ε i : ℂ) * R4C.resolvC (W0 Y) z k l else 0 := by
      intro k
      change ((if k = σ i then ε i else 0 : ℝ) : ℂ) * _ = _
      by_cases h : k = σ i <;> simp [h]
    simp only [houter]
    rw [Finset.sum_ite_eq' Finset.univ (σ i)
      fun k => (ε i : ℂ) * R4C.resolvC (W0 Y) z k l]
    change (if σ i ∈ Finset.univ then (ε i : ℂ) * R4C.resolvC (W0 Y) z (σ i) l else 0) *
      ((if l = σ j then ε j else 0 : ℝ) : ℂ) = _
    by_cases h : l = σ j <;> simp [h]
  simp only [hinner]
  rw [Finset.sum_ite_eq' Finset.univ (σ j)
    fun l => (ε i : ℂ) * R4C.resolvC (W0 Y) z (σ i) l * (ε j : ℂ)]
  simp only [Finset.mem_univ, if_true]
  ring

/-- Entrywise mean of the resolvent (modeling choice 2: no Bochner integral of matrices). -/
noncomputable def meanResolv (p d : ℕ) (z : ℂ) : Matrix (Fin d) (Fin d) ℂ :=
  Matrix.of fun i j => ∫ Y, R4C.resolvC (W0 Y) z i j ∂gaussianMatrix p d

theorem resolvC_apply_eq_cformC (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (i j : Fin d) :
    R4C.resolvC W z i j = R4C.cformC W z (Pi.single i 1) (Pi.single j 1) := by
  rw [cformC_eq_bil, bil_eq_sum]
  set ei : Fin d → ℝ := Pi.single i 1 with hei
  set ej : Fin d → ℝ := Pi.single j 1 with hej
  have h1 : ∀ a : Fin d, ((ei a : ℝ) : ℂ) = if a = i then 1 else 0 := by
    intro a
    by_cases h : a = i <;> simp [hei, h]
  have h2 : ∀ b : Fin d, ((ej b : ℝ) : ℂ) = if b = j then 1 else 0 := by
    intro b
    by_cases h : b = j <;> simp [hej, h]
  simp only [h1, h2, ite_mul, one_mul, zero_mul, mul_ite, mul_one, mul_zero,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]

theorem norm_resolvC_entry_le {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (i j : Fin d) : ‖R4C.resolvC W z i j‖ ≤ 1 / z.im := by
  have hs : ∀ k : Fin d, (Pi.single k (1 : ℝ)) ⬝ᵥ (Pi.single k (1 : ℝ)) = 1 := by
    intro k
    simp [dotProduct, Pi.single_apply]
  have h := R4C.norm_cformC_le hW hz (Pi.single i (1 : ℝ)) (Pi.single j (1 : ℝ))
  rw [hs, hs] at h
  rw [resolvC_apply_eq_cformC]
  simpa using h

theorem integrable_resolvC_entry (hz : 0 < z.im) (i j : Fin d) :
    Integrable (fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.resolvC (W0 Y) z i j)
      (gaussianMatrix p d) :=
  Integrable.of_bound (measurable_resolvC_entry z i j).aestronglyMeasurable (1 / z.im)
    (Filter.Eventually.of_forall fun Y => norm_resolvC_entry_le (isHermitian_W0 Y) hz i j)

theorem meanResolv_sgnPerm (hz : 0 < z.im) (σ : Equiv.Perm (Fin d)) (ε : Fin d → ℝ)
    (hε : ∀ j, ε j = 1 ∨ ε j = -1) (i j : Fin d) :
    meanResolv p d z i j = (ε i : ℂ) * (ε j : ℂ) * meanResolv p d z (σ i) (σ j) := by
  have hmp := measurePreserving_sgnPerm (p := p) σ hε
  have hchange : ∫ Y, R4C.resolvC (W0 Y) z i j ∂(gaussianMatrix p d)
      = ∫ Y, R4C.resolvC (W0 (sgnPerm σ ε Y)) z i j ∂(gaussianMatrix p d) := by
    have h := integral_map (φ := (sgnPerm σ ε : Matrix (Fin p) (Fin d) ℝ → _))
      (μ := gaussianMatrix p d) hmp.measurable.aemeasurable
      (f := fun Y => R4C.resolvC (W0 Y) z i j)
      (by rw [hmp.map_eq]; exact (measurable_resolvC_entry z i j).aestronglyMeasurable)
    rw [hmp.map_eq] at h
    exact h
  change ∫ Y, R4C.resolvC (W0 Y) z i j ∂(gaussianMatrix p d) = _
  rw [hchange]
  simp only [resolvC_W0_sgnPerm σ ε hε hz]
  rw [integral_const_mul]
  rfl

theorem meanResolv_offDiag (hz : 0 < z.im) {i j : Fin d} (hij : i ≠ j) :
    meanResolv p d z i j = 0 := by
  set ε : Fin d → ℝ := fun l => if l = j then -1 else 1 with hεdef
  have hε : ∀ l, ε l = 1 ∨ ε l = -1 := by
    intro l
    by_cases h : l = j <;> simp [hεdef, h]
  have h := meanResolv_sgnPerm (p := p) hz (Equiv.refl (Fin d)) ε hε i j
  simp only [hεdef, Equiv.refl_apply, if_neg hij, if_pos rfl] at h
  push_cast at h
  linear_combination h / 2

theorem meanResolv_diag_eq (hz : 0 < z.im) (i i' : Fin d) :
    meanResolv p d z i i = meanResolv p d z i' i' := by
  have h := meanResolv_sgnPerm (p := p) hz (Equiv.swap i i') (fun _ => 1) (fun _ => Or.inl rfl)
    i i
  simpa [Equiv.swap_apply_left] using h

theorem meanResolv_diag (hz : 0 < z.im) (hd : 0 < d) (i : Fin d) :
    meanResolv p d z i i = ∫ Y, R4C.stieltjesC (W0 Y) z ∂gaussianMatrix p d := by
  have hne : ((d : ℂ)) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have hsum : ∫ Y, R4C.stieltjesC (W0 Y) z ∂(gaussianMatrix p d)
      = (d : ℂ)⁻¹ * ∑ k, meanResolv p d z k k := by
    have hb : ∀ Y : Matrix (Fin p) (Fin d) ℝ,
        R4C.stieltjesC (W0 Y) z = (d : ℂ)⁻¹ * ∑ k, R4C.resolvC (W0 Y) z k k := by
      intro Y
      rw [R4C.stieltjesC, Matrix.trace]
      rfl
    simp only [hb]
    rw [integral_const_mul, integral_finsetSum _ fun k _ => integrable_resolvC_entry hz k k]
    rfl
  rw [hsum]
  have hall : ∀ k : Fin d, meanResolv p d z k k = meanResolv p d z i i :=
    fun k => meanResolv_diag_eq hz k i
  rw [Finset.sum_congr rfl fun k _ => hall k, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul]
  field_simp

/-- **R2a.** Exact centering, for every pair of vectors. -/
theorem integral_cformC_eq (hz : 0 < z.im) (hd : 0 < d) (v w : Fin d → ℝ) :
    ∫ Y, R4C.cformC (W0 Y) z v w ∂gaussianMatrix p d
      = ((v ⬝ᵥ w : ℝ) : ℂ) * ∫ Y, R4C.stieltjesC (W0 Y) z ∂gaussianMatrix p d := by
  have hint : ∀ (i j : Fin d), Integrable
      (fun Y : Matrix (Fin p) (Fin d) ℝ =>
        (v i : ℂ) * (R4C.resolvC (W0 Y) z i j * (w j : ℂ))) (gaussianMatrix p d) :=
    fun i j => (((integrable_resolvC_entry hz i j).mul_const (w j : ℂ)).const_mul (v i : ℂ))
  have hexp : ∀ Y : Matrix (Fin p) (Fin d) ℝ, R4C.cformC (W0 Y) z v w
      = ∑ i, ∑ j, (v i : ℂ) * (R4C.resolvC (W0 Y) z i j * (w j : ℂ)) := by
    intro Y
    rw [cformC_eq_bil, bil_eq_sum]
    exact Finset.sum_congr rfl fun i _ => by rw [Finset.mul_sum]
  simp only [hexp]
  rw [integral_finsetSum _ fun i _ => integrable_finsetSum _ fun j _ => hint i j]
  have hrow : ∀ i : Fin d,
      ∫ Y, ∑ j, (v i : ℂ) * (R4C.resolvC (W0 Y) z i j * (w j : ℂ)) ∂(gaussianMatrix p d)
        = (v i : ℂ) * (w i : ℂ) * meanResolv p d z i i := by
    intro i
    rw [integral_finsetSum _ fun j _ => hint i j]
    have hterm : ∀ j : Fin d,
        ∫ Y, (v i : ℂ) * (R4C.resolvC (W0 Y) z i j * (w j : ℂ)) ∂(gaussianMatrix p d)
          = if j = i then (v i : ℂ) * (w i : ℂ) * meanResolv p d z i i else 0 := by
      intro j
      rw [integral_const_mul, MeasureTheory.integral_mul_const]
      by_cases hji : j = i
      · rw [if_pos hji, hji]
        simp only [meanResolv, Matrix.of_apply]
        ring
      · rw [if_neg hji]
        have hz0 : (meanResolv p d z i j) = 0 := meanResolv_offDiag hz (Ne.symm hji)
        simp only [meanResolv, Matrix.of_apply] at hz0
        rw [hz0]
        ring
    simp only [hterm]
    rw [Finset.sum_ite_eq' Finset.univ i fun _ => (v i : ℂ) * (w i : ℂ) * meanResolv p d z i i]
    simp
  simp only [hrow]
  simp only [meanResolv_diag hz hd]
  rw [← Finset.sum_mul]
  congr 1
  rw [dotProduct]
  push_cast
  rfl

end MeanResolvent

end R2
end StackedSVD
