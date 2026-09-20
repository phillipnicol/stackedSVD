/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R4
import StackedSVD.RMT.MP

/-!
# R4C: the complex resolvent, and the missing deterministic facts of items R5 and T

This file owns two things that several Layer 2 notes duplicate.

1. **Additions to the real layer `R4`** that item R5 needs (`notes/archive/rmt_R5.md`, modeling
   choice 1): the squared quadratic form `qform2`, the cross forms `cform` and `cform2`, the
   two polarization identities, and `qform2_antitoneOn`, the fact that `z ↦ y ⬝ᵥ G₀(z)² y`
   falls above the top eigenvalue. `RMT/R4.lean` proves that `qform` rises, never that
   `qform2` falls.
2. **The complex resolvent** `resolvC W z = (W - z I)⁻¹` over `ℂ`, its spectral form (the
   complex twin of `R4.resolv_eq_conj`), the deterministic bounds `‖·‖ ≤ 1/Im z` for the
   bilinear forms and the trace, the identification with `R4.resolv` at real `z` above the
   spectrum, and the two scalar bounds of item T. `notes/archive/rmt_R1.md` and
   `notes/archive/rmt_T.md` each define `resolvC` and the complex `resolv_eq_conj`; both should
   import them from here.

No probability and no asymptotics. Every bound holds for every realization.

Naming for the coordinator: `R1.resolvC`, `R1.stieltjes`, `R1.stieltjes₂`, `T.resolvC`,
`T.cmat` and `T.bilC` become `R4C.resolvC`, `R4C.stieltjesC`, `R4C.stieltjes2C`, `R4C.cmat`
and `R4C.cformC` / `R4C.cform2C`. `R4.qform2`, `R4.cform`, `R4.cform2`,
`R4.cform_eq_polarization`, `R4.cform2_eq_polarization`, `R4.qform2_antitoneOn` keep the names
of `notes/archive/rmt_R5.md`.
-/

open Filter Topology
open scoped Matrix

namespace StackedSVD

namespace R4

variable {d : ℕ} {W₀ M : Matrix (Fin d) (Fin d) ℝ} {z z₁ z₂ : ℝ} {x y : Fin d → ℝ}

/-! ### The three extra deterministic forms of item R5 -/

/-- `Φ²_y z = y ⬝ᵥ (G₀ z)² y`. -/
noncomputable def qform2 (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (y : Fin d → ℝ) : ℝ :=
  y ⬝ᵥ ((resolv W z * resolv W z) *ᵥ y)

/-- `Ψ_{x,y} z = x ⬝ᵥ G₀ z y`. -/
noncomputable def cform (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (x y : Fin d → ℝ) : ℝ :=
  x ⬝ᵥ (resolv W z *ᵥ y)

/-- `Ψ²_{x,y} z = x ⬝ᵥ (G₀ z)² y`. -/
noncomputable def cform2 (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (x y : Fin d → ℝ) : ℝ :=
  x ⬝ᵥ ((resolv W z * resolv W z) *ᵥ y)

theorem cform_self (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (y : Fin d → ℝ) :
    cform W z y y = qform W z y := rfl

theorem cform2_self (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (y : Fin d → ℝ) :
    cform2 W z y y = qform2 W z y := rfl

/-- A symmetric kernel gives a symmetric bilinear form. -/
theorem dotProduct_mulVec_comm (hM : Mᵀ = M) (x y : Fin d → ℝ) :
    x ⬝ᵥ (M *ᵥ y) = y ⬝ᵥ (M *ᵥ x) := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hM, dotProduct_comm]

/-- **Polarization.** A symmetric bilinear form is a combination of three quadratic forms. -/
theorem polarization (hM : Mᵀ = M) (x y : Fin d → ℝ) :
    x ⬝ᵥ (M *ᵥ y) =
      ((x + y) ⬝ᵥ (M *ᵥ (x + y)) - x ⬝ᵥ (M *ᵥ x) - y ⬝ᵥ (M *ᵥ y)) / 2 := by
  have h := dotProduct_mulVec_comm hM x y
  rw [Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add]
  linarith

/-- Polarization for `Ψ_{x,y}`, the identity `notes/archive/rmt_R5.md` names. -/
theorem cform_eq_polarization (hW₀ : W₀.IsHermitian) (z : ℝ) (x y : Fin d → ℝ) :
    cform W₀ z x y = (qform W₀ z (x + y) - qform W₀ z x - qform W₀ z y) / 2 :=
  polarization (transpose_resolv hW₀) x y

/-- Polarization for `Ψ²_{x,y}`. -/
theorem cform2_eq_polarization (hW₀ : W₀.IsHermitian) (z : ℝ) (x y : Fin d → ℝ) :
    cform2 W₀ z x y = (qform2 W₀ z (x + y) - qform2 W₀ z x - qform2 W₀ z y) / 2 := by
  refine polarization ?_ x y
  rw [Matrix.transpose_mul, transpose_resolv hW₀]

/-- `Φ²_y` in the eigenbasis. -/
theorem qform2_eq_sum (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (y : Fin d → ℝ) :
    qform2 W₀ z y = ∑ a, ((hW₀.eigenvalues a - z)⁻¹) ^ 2 * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 := by
  rw [qform2, resolv_mul_resolv_eq_conj hW₀ hz, dotProduct_conj]

theorem qform2_nonneg (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (y : Fin d → ℝ) :
    0 ≤ qform2 W₀ z y := by
  rw [qform2_eq_sum hW₀ hz]
  exact Finset.sum_nonneg fun a _ => by positivity

/-- **The new deterministic fact of item R5.** `Φ²_y` falls above the top eigenvalue. -/
theorem qform2_antitoneOn (hW₀ : W₀.IsHermitian) (y : Fin d → ℝ) :
    AntitoneOn (fun t => qform2 W₀ t y) (Set.Ioi (lamMax W₀ hW₀)) := by
  intro a ha b hb hab
  simp only [Set.mem_Ioi] at ha hb
  simp only
  rw [qform2_eq_sum hW₀ hb, qform2_eq_sum hW₀ ha]
  refine Finset.sum_le_sum fun i _ => ?_
  have hle : hW₀.eigenvalues i ≤ lamMax W₀ hW₀ := eigenvalues_le_lamMax hW₀ i
  have h1 : 0 < a - hW₀.eigenvalues i := by linarith
  have h2 : 0 < b - hW₀.eigenvalues i := by linarith
  have e1 : ((hW₀.eigenvalues i - b)⁻¹) ^ 2 = ((b - hW₀.eigenvalues i) ^ 2)⁻¹ := by
    rw [inv_pow]
    congr 1
    ring
  have e2 : ((hW₀.eigenvalues i - a)⁻¹) ^ 2 = ((a - hW₀.eigenvalues i) ^ 2)⁻¹ := by
    rw [inv_pow]
    congr 1
    ring
  have hmono : ((b - hW₀.eigenvalues i) ^ 2)⁻¹ ≤ ((a - hW₀.eigenvalues i) ^ 2)⁻¹ := by
    have hsq : (a - hW₀.eigenvalues i) ^ 2 ≤ (b - hW₀.eigenvalues i) ^ 2 := by nlinarith
    have hpos : (0 : ℝ) < (a - hW₀.eigenvalues i) ^ 2 := by positivity
    exact inv_anti₀ hpos hsq
  rw [e1, e2]
  exact mul_le_mul_of_nonneg_right hmono (sq_nonneg _)

end R4

/-! ### The complex resolvent -/

namespace R4C

open R4

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} {z : ℂ} {t : ℝ}

/-- A real matrix read over `ℂ`. -/
noncomputable def cmat (W : Matrix (Fin d) (Fin d) ℝ) : Matrix (Fin d) (Fin d) ℂ :=
  W.map (fun a => (a : ℂ))

/-- A real vector read over `ℂ`. -/
def cvec (v : Fin d → ℝ) : Fin d → ℂ := fun a => (v a : ℂ)

/-- **The complex resolvent** `G(z) = (W - z I)⁻¹`. Junk value `0` on the spectrum, which is
never reached when `Im z ≠ 0`. -/
noncomputable def resolvC (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) : Matrix (Fin d) (Fin d) ℂ :=
  (cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ))⁻¹

/-- `Ψ_{x,y}(z) = x ⬝ᵥ G(z) y` for real `x`, `y`. -/
noncomputable def cformC (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x y : Fin d → ℝ) : ℂ :=
  cvec x ⬝ᵥ (resolvC W z *ᵥ cvec y)

/-- `Ψ²_{x,y}(z) = x ⬝ᵥ G(z)² y` for real `x`, `y`. -/
noncomputable def cform2C (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (x y : Fin d → ℝ) : ℂ :=
  cvec x ⬝ᵥ ((resolvC W z * resolvC W z) *ᵥ cvec y)

/-- `Φ_y(z) = y ⬝ᵥ G(z) y`. -/
noncomputable def qformC (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (y : Fin d → ℝ) : ℂ :=
  cformC W z y y

/-- `Φ²_y(z) = y ⬝ᵥ G(z)² y`. -/
noncomputable def qform2C (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) (y : Fin d → ℝ) : ℂ :=
  cform2C W z y y

/-- `s(z) = d⁻¹ tr G(z)`, the Stieltjes transform of the empirical spectral law. -/
noncomputable def stieltjesC (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) : ℂ :=
  (d : ℂ)⁻¹ * (resolvC W z).trace

/-- `d⁻¹ tr G(z)²`. -/
noncomputable def stieltjes2C (W : Matrix (Fin d) (Fin d) ℝ) (z : ℂ) : ℂ :=
  (d : ℂ)⁻¹ * (resolvC W z * resolvC W z).trace

/-! #### `cmat` is a ring map -/

theorem cmat_mul (A B : Matrix (Fin d) (Fin d) ℝ) : cmat (A * B) = cmat A * cmat B := by
  ext i j
  simp only [cmat, Matrix.map_apply, Matrix.mul_apply]
  push_cast
  rfl

theorem cmat_one : cmat (1 : Matrix (Fin d) (Fin d) ℝ) = 1 :=
  Matrix.map_one _ (by simp) (by simp)

theorem cmat_transpose (A : Matrix (Fin d) (Fin d) ℝ) : cmat Aᵀ = (cmat A)ᵀ := rfl

theorem cmat_diagonal (f : Fin d → ℝ) :
    cmat (Matrix.diagonal f) = Matrix.diagonal fun a => (f a : ℂ) :=
  Matrix.diagonal_map (by simp)

theorem cmat_sub_smul (W : Matrix (Fin d) (Fin d) ℝ) (t : ℝ) :
    cmat (W - t • (1 : Matrix (Fin d) (Fin d) ℝ)) = cmat W - (t : ℂ) • 1 := by
  ext i j
  by_cases h : i = j
  · subst h; simp [cmat, Matrix.one_apply_eq]
  · simp [cmat, Matrix.one_apply_ne h]

theorem dotProduct_cmat_mulVec (A : Matrix (Fin d) (Fin d) ℝ) (x y : Fin d → ℝ) :
    cvec x ⬝ᵥ (cmat A *ᵥ cvec y) = ((x ⬝ᵥ (A *ᵥ y) : ℝ) : ℂ) := by
  simp only [cvec, cmat, dotProduct, Matrix.mulVec, Matrix.map_apply]
  push_cast
  rfl

theorem transpose_cmat_mulVec_cvec (V : Matrix (Fin d) (Fin d) ℝ) (x : Fin d → ℝ) :
    (cmat V)ᵀ *ᵥ cvec x = cvec (Vᵀ *ᵥ x) := by
  funext a
  simp only [cvec, cmat, Matrix.mulVec, dotProduct, Matrix.transpose_apply, Matrix.map_apply]
  push_cast
  rfl

/-! #### The spectral form -/

theorem transpose_ceigU_mul (hW : W.IsHermitian) :
    (cmat (eigU hW))ᵀ * cmat (eigU hW) = 1 := by
  rw [← cmat_transpose, ← cmat_mul, transpose_eigU_mul hW, cmat_one]

theorem ceigU_mul_transpose (hW : W.IsHermitian) :
    cmat (eigU hW) * (cmat (eigU hW))ᵀ = 1 := by
  rw [← cmat_transpose, ← cmat_mul, eigU_mul_transpose hW, cmat_one]

theorem cmat_sub_smul_one_eq_conj (hW : W.IsHermitian) (z : ℂ) :
    cmat W - z • (1 : Matrix (Fin d) (Fin d) ℂ)
      = cmat (eigU hW) * Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z)
        * (cmat (eigU hW))ᵀ := by
  have hconj : cmat W
      = cmat (eigU hW) * Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ))
        * (cmat (eigU hW))ᵀ := by
    conv_lhs => rw [← eigU_conj hW]
    rw [cmat_mul, cmat_mul, cmat_diagonal, cmat_transpose]
  have hdiag : Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z)
      = Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ)) - z • 1 := by
    ext i j
    rcases eq_or_ne i j with h | h
    · subst h; simp
    · simp [Matrix.diagonal_apply_ne _ h, Matrix.one_apply_ne h]
  rw [hdiag, Matrix.mul_sub, Matrix.sub_mul, ← hconj]
  congr 1
  rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one, ceigU_mul_transpose]

theorem eigenvalue_sub_ne_zero (hW : W.IsHermitian) (hz : z.im ≠ 0) (a : Fin d) :
    (hW.eigenvalues a : ℂ) - z ≠ 0 := by
  intro hcon
  apply hz
  have h : ((hW.eigenvalues a : ℂ) - z).im = 0 := by rw [hcon]; simp
  simpa using h

/-- **Master lemma, complex form.** `G(z) = U diag((λ_a - z)⁻¹) Uᵀ` off the real spectrum. -/
theorem resolvC_eq_conj (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    resolvC W z
      = cmat (eigU hW) * Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)
        * (cmat (eigU hW))ᵀ := by
  have hne := eigenvalue_sub_ne_zero hW hz
  have hDD : Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
      Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [Matrix.diagonal_mul_diagonal,
      show (fun a => ((hW.eigenvalues a : ℂ) - z) * ((hW.eigenvalues a : ℂ) - z)⁻¹)
        = fun _ => (1 : ℂ) from funext fun a => mul_inv_cancel₀ (hne a)]
    exact Matrix.diagonal_one
  refine Matrix.inv_eq_right_inv ?_
  rw [cmat_sub_smul_one_eq_conj hW]
  have hassoc : cmat (eigU hW) * Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
      (cmat (eigU hW))ᵀ *
      (cmat (eigU hW) * Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
        (cmat (eigU hW))ᵀ)
      = cmat (eigU hW) * (Matrix.diagonal (fun a => (hW.eigenvalues a : ℂ) - z) *
          ((cmat (eigU hW))ᵀ * cmat (eigU hW)) *
          Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)) * (cmat (eigU hW))ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, transpose_ceigU_mul, Matrix.mul_one, hDD, Matrix.mul_one, ceigU_mul_transpose]

theorem resolvC_mul_resolvC_eq_conj (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    resolvC W z * resolvC W z
      = cmat (eigU hW) * Matrix.diagonal (fun a => (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2)
        * (cmat (eigU hW))ᵀ := by
  rw [resolvC_eq_conj hW hz]
  have hassoc : cmat (eigU hW) * Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
      (cmat (eigU hW))ᵀ *
      (cmat (eigU hW) * Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
        (cmat (eigU hW))ᵀ)
      = cmat (eigU hW) * (Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹) *
          ((cmat (eigU hW))ᵀ * cmat (eigU hW)) *
          Matrix.diagonal (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹)) * (cmat (eigU hW))ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, transpose_ceigU_mul, Matrix.mul_one, Matrix.diagonal_mul_diagonal,
    show (fun a => ((hW.eigenvalues a : ℂ) - z)⁻¹ * ((hW.eigenvalues a : ℂ) - z)⁻¹)
      = fun a => (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 from funext fun a => (sq _).symm]

/-! #### The bilinear forms in the eigenbasis -/

theorem conj_mulVec_gen {R : Type*} [CommRing R] (V : Matrix (Fin d) (Fin d) R)
    (f y : Fin d → R) :
    (V * Matrix.diagonal f * Vᵀ) *ᵥ y = V *ᵥ fun a => f a * (Vᵀ *ᵥ y) a := by
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
  congr 1
  funext a
  simp [Matrix.mulVec_diagonal]

theorem dotProduct_conj_gen {R : Type*} [CommRing R] (V : Matrix (Fin d) (Fin d) R)
    (f x y : Fin d → R) :
    x ⬝ᵥ ((V * Matrix.diagonal f * Vᵀ) *ᵥ y) = ∑ a, f a * ((Vᵀ *ᵥ x) a * (Vᵀ *ᵥ y) a) := by
  rw [conj_mulVec_gen, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]
  simp only [dotProduct]
  exact Finset.sum_congr rfl fun a _ => by ring

theorem cformC_eq_sum (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin d → ℝ) :
    cformC W z x y = ∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹ *
      (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ) := by
  rw [cformC, resolvC_eq_conj hW hz, dotProduct_conj_gen, transpose_cmat_mulVec_cvec,
    transpose_cmat_mulVec_cvec]
  refine Finset.sum_congr rfl fun a _ => ?_
  simp only [cvec]
  push_cast
  ring

theorem cform2C_eq_sum (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin d → ℝ) :
    cform2C W z x y = ∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
      (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ) := by
  rw [cform2C, resolvC_mul_resolvC_eq_conj hW hz, dotProduct_conj_gen,
    transpose_cmat_mulVec_cvec, transpose_cmat_mulVec_cvec]
  refine Finset.sum_congr rfl fun a _ => ?_
  simp only [cvec]
  push_cast
  ring

/-! #### The deterministic bounds -/

theorem norm_inv_eigenvalue_sub_le (hW : W.IsHermitian) (hz : 0 < z.im) (a : Fin d) :
    ‖((hW.eigenvalues a : ℂ) - z)⁻¹‖ ≤ (z.im)⁻¹ := by
  have him : z.im ≤ ‖(hW.eigenvalues a : ℂ) - z‖ := by
    have h := Complex.abs_im_le_norm ((hW.eigenvalues a : ℂ) - z)
    have him' : ((hW.eigenvalues a : ℂ) - z).im = -z.im := by simp
    rw [him', abs_neg, abs_of_pos hz] at h
    exact h
  rw [norm_inv]
  exact inv_anti₀ hz him

/-- Cauchy-Schwarz in the eigenbasis. -/
theorem sum_abs_coords_le (hW : W.IsHermitian) (x y : Fin d → ℝ) :
    ∑ a, |((eigU hW)ᵀ *ᵥ x) a| * |((eigU hW)ᵀ *ᵥ y) a|
      ≤ Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) := by
  set c : Fin d → ℝ := (eigU hW)ᵀ *ᵥ x with hc
  set l : Fin d → ℝ := (eigU hW)ᵀ *ᵥ y with hl
  have hcx : ∑ a, c a ^ 2 = x ⬝ᵥ x := by
    have h := dotProduct_transpose_eigU hW x
    rw [← h, dotProduct]
    exact Finset.sum_congr rfl fun a _ => (sq (c a))
  have hly : ∑ a, l a ^ 2 = y ⬝ᵥ y := by
    have h := dotProduct_transpose_eigU hW y
    rw [← h, dotProduct]
    exact Finset.sum_congr rfl fun a _ => (sq (l a))
  have hcs : (∑ a, |c a| * |l a|) ^ 2 ≤ (∑ a, c a ^ 2) * ∑ a, l a ^ 2 := by
    have h := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun a => |c a|) fun a => |l a|
    simpa [sq_abs] using h
  have hnn : 0 ≤ ∑ a, |c a| * |l a| := Finset.sum_nonneg fun a _ => by positivity
  have hx0 : 0 ≤ ∑ a, c a ^ 2 := Finset.sum_nonneg fun a _ => sq_nonneg _
  calc ∑ a, |c a| * |l a| = Real.sqrt ((∑ a, |c a| * |l a|) ^ 2) := (Real.sqrt_sq hnn).symm
    _ ≤ Real.sqrt ((∑ a, c a ^ 2) * ∑ a, l a ^ 2) := Real.sqrt_le_sqrt hcs
    _ = Real.sqrt (∑ a, c a ^ 2) * Real.sqrt (∑ a, l a ^ 2) := Real.sqrt_mul hx0 _
    _ = Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) := by rw [hcx, hly]

/-- **The first order bound.** `|x ⬝ᵥ G(z) y| ≤ ‖x‖ ‖y‖ / Im z`. -/
theorem norm_cformC_le (hW : W.IsHermitian) (hz : 0 < z.im) (x y : Fin d → ℝ) :
    ‖cformC W z x y‖ ≤ Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) / z.im := by
  rw [cformC_eq_sum hW (ne_of_gt hz)]
  calc ‖∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹ *
          (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖
      ≤ ∑ a, ‖((hW.eigenvalues a : ℂ) - z)⁻¹ *
          (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖ := norm_sum_le _ _
    _ ≤ ∑ a, (z.im)⁻¹ *
          (|((eigU hW)ᵀ *ᵥ x) a| * |((eigU hW)ᵀ *ᵥ y) a|) := by
        refine Finset.sum_le_sum fun a _ => ?_
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul]
        exact mul_le_mul_of_nonneg_right (norm_inv_eigenvalue_sub_le hW hz a) (by positivity)
    _ = (z.im)⁻¹ * ∑ a, |((eigU hW)ᵀ *ᵥ x) a| * |((eigU hW)ᵀ *ᵥ y) a| := by
        rw [Finset.mul_sum]
    _ ≤ (z.im)⁻¹ * (Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y)) := by
        exact mul_le_mul_of_nonneg_left (sum_abs_coords_le hW x y) (by positivity)
    _ = Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) / z.im := by
        field_simp

/-- **The second order bound.** `|x ⬝ᵥ G(z)² y| ≤ ‖x‖ ‖y‖ / (Im z)²`. -/
theorem norm_cform2C_le (hW : W.IsHermitian) (hz : 0 < z.im) (x y : Fin d → ℝ) :
    ‖cform2C W z x y‖ ≤ Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) / z.im ^ 2 := by
  rw [cform2C_eq_sum hW (ne_of_gt hz)]
  calc ‖∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
          (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖
      ≤ ∑ a, ‖(((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 *
          (((((eigU hW)ᵀ *ᵥ x) a * ((eigU hW)ᵀ *ᵥ y) a : ℝ)) : ℂ)‖ := norm_sum_le _ _
    _ ≤ ∑ a, ((z.im)⁻¹) ^ 2 *
          (|((eigU hW)ᵀ *ᵥ x) a| * |((eigU hW)ᵀ *ᵥ y) a|) := by
        refine Finset.sum_le_sum fun a _ => ?_
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul, norm_pow]
        refine mul_le_mul_of_nonneg_right ?_ (by positivity)
        exact pow_le_pow_left₀ (norm_nonneg _) (norm_inv_eigenvalue_sub_le hW hz a) 2
    _ = ((z.im)⁻¹) ^ 2 * ∑ a, |((eigU hW)ᵀ *ᵥ x) a| * |((eigU hW)ᵀ *ᵥ y) a| := by
        rw [Finset.mul_sum]
    _ ≤ ((z.im)⁻¹) ^ 2 * (Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y)) := by
        exact mul_le_mul_of_nonneg_left (sum_abs_coords_le hW x y) (by positivity)
    _ = Real.sqrt (x ⬝ᵥ x) * Real.sqrt (y ⬝ᵥ y) / z.im ^ 2 := by
        rw [inv_pow]
        ring

theorem norm_qformC_le (hW : W.IsHermitian) (hz : 0 < z.im) (y : Fin d → ℝ) :
    ‖qformC W z y‖ ≤ (y ⬝ᵥ y) / z.im := by
  have h := norm_cformC_le hW hz y y
  have hy : 0 ≤ y ⬝ᵥ y := by
    rw [dotProduct]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rwa [Real.mul_self_sqrt hy] at h

theorem norm_qform2C_le (hW : W.IsHermitian) (hz : 0 < z.im) (y : Fin d → ℝ) :
    ‖qform2C W z y‖ ≤ (y ⬝ᵥ y) / z.im ^ 2 := by
  have h := norm_cform2C_le hW hz y y
  have hy : 0 ≤ y ⬝ᵥ y := by
    rw [dotProduct]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  rwa [Real.mul_self_sqrt hy] at h

/-! #### The trace bounds -/

theorem trace_conj (V D : Matrix (Fin d) (Fin d) ℂ) (hV : Vᵀ * V = 1) :
    (V * D * Vᵀ).trace = D.trace := by
  rw [Matrix.trace_mul_comm, ← Matrix.mul_assoc, hV, Matrix.one_mul]

theorem trace_resolvC (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    (resolvC W z).trace = ∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹ := by
  rw [resolvC_eq_conj hW hz, trace_conj _ _ (transpose_ceigU_mul hW), Matrix.trace_diagonal]

theorem trace_resolvC_sq (hW : W.IsHermitian) (hz : z.im ≠ 0) :
    (resolvC W z * resolvC W z).trace = ∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 := by
  rw [resolvC_mul_resolvC_eq_conj hW hz, trace_conj _ _ (transpose_ceigU_mul hW),
    Matrix.trace_diagonal]

theorem norm_stieltjesC_le (hW : W.IsHermitian) (hz : 0 < z.im) (hd : 0 < d) :
    ‖stieltjesC W z‖ ≤ 1 / z.im := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  rw [stieltjesC, trace_resolvC hW (ne_of_gt hz), norm_mul, norm_inv]
  have h1 : ‖∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹‖ ≤ (d : ℝ) * (z.im)⁻¹ := by
    calc ‖∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹‖
        ≤ ∑ _a : Fin d, (z.im)⁻¹ :=
          (norm_sum_le _ _).trans
            (Finset.sum_le_sum fun a _ => norm_inv_eigenvalue_sub_le hW hz a)
      _ = (d : ℝ) * (z.im)⁻¹ := by simp
  have hnd : ‖(d : ℂ)‖ = (d : ℝ) := by simp
  rw [hnd]
  calc ((d : ℝ))⁻¹ * ‖∑ a, ((hW.eigenvalues a : ℂ) - z)⁻¹‖
      ≤ ((d : ℝ))⁻¹ * ((d : ℝ) * (z.im)⁻¹) :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = 1 / z.im := by field_simp

theorem norm_stieltjes2C_le (hW : W.IsHermitian) (hz : 0 < z.im) (hd : 0 < d) :
    ‖stieltjes2C W z‖ ≤ 1 / z.im ^ 2 := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  rw [stieltjes2C, trace_resolvC_sq hW (ne_of_gt hz), norm_mul, norm_inv]
  have h1 : ‖∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖ ≤ (d : ℝ) * ((z.im)⁻¹) ^ 2 := by
    calc ‖∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖
        ≤ ∑ _a : Fin d, ((z.im)⁻¹) ^ 2 :=
          (norm_sum_le _ _).trans (Finset.sum_le_sum fun a _ => by
            rw [norm_pow]
            exact pow_le_pow_left₀ (norm_nonneg _) (norm_inv_eigenvalue_sub_le hW hz a) 2)
      _ = (d : ℝ) * ((z.im)⁻¹) ^ 2 := by simp
  have hnd : ‖(d : ℂ)‖ = (d : ℝ) := by simp
  rw [hnd]
  calc ((d : ℝ))⁻¹ * ‖∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2‖
      ≤ ((d : ℝ))⁻¹ * ((d : ℝ) * ((z.im)⁻¹) ^ 2) :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = 1 / z.im ^ 2 := by
        rw [inv_pow]
        field_simp

/-! #### At a real point above the spectrum, `resolvC` is `R4.resolv` -/

theorem resolvC_ofReal (hW : W.IsHermitian) (ht : lamMax W hW < t) :
    resolvC W (t : ℂ) = cmat (resolv W t) := by
  refine Matrix.inv_eq_right_inv ?_
  rw [← cmat_sub_smul, ← cmat_mul, mul_resolv hW ht, cmat_one]

theorem cformC_ofReal (hW : W.IsHermitian) (ht : lamMax W hW < t) (x y : Fin d → ℝ) :
    cformC W (t : ℂ) x y = ((cform W t x y : ℝ) : ℂ) := by
  rw [cformC, resolvC_ofReal hW ht, dotProduct_cmat_mulVec]
  rfl

theorem cform2C_ofReal (hW : W.IsHermitian) (ht : lamMax W hW < t) (x y : Fin d → ℝ) :
    cform2C W (t : ℂ) x y = ((cform2 W t x y : ℝ) : ℂ) := by
  rw [cform2C, resolvC_ofReal hW ht, ← cmat_mul, dotProduct_cmat_mulVec]
  rfl

/-! ### Item T: the two termwise scalar bounds

`A = λ - x` is the real gap, `B = λ - x - iη` its complex shift, and `g` a lower bound on
`x - λ`. The constants are the ones of `notes/archive/rmt_T.md`; the note multiplies them by `4` and
`16` because its hypotheses give the gap `ε/2`. -/

private theorem norm_shift_sq {a η : ℝ} :
    ‖((-a : ℝ) : ℂ) - η * Complex.I‖ ^ 2 = a ^ 2 + η ^ 2 := by
  have hre : (((-a : ℝ) : ℂ) - η * Complex.I).re = -a := by simp
  have him : (((-a : ℝ) : ℂ) - η * Complex.I).im = -η := by simp
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, hre, him]
  ring

private theorem norm_sum_shift_sq {a η : ℝ} :
    ‖((-a : ℝ) : ℂ) + (((-a : ℝ) : ℂ) - η * Complex.I)‖ ^ 2 = 4 * a ^ 2 + η ^ 2 := by
  have hre : (((-a : ℝ) : ℂ) + (((-a : ℝ) : ℂ) - η * Complex.I)).re = -(2 * a) := by
    simp; ring
  have him : (((-a : ℝ) : ℂ) + (((-a : ℝ) : ℂ) - η * Complex.I)).im = -η := by simp
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, hre, him]
  ring

/-- **Item T, first order.** `|(λ - x)⁻¹ - (λ - x - iη)⁻¹| ≤ η / g²` when `x - λ ≥ g > 0`. -/
theorem norm_inv_sub_inv_le {a η g : ℝ} (hη : 0 < η) (hg : 0 < g) (hga : g ≤ a) :
    ‖((-a : ℝ) : ℂ)⁻¹ - (((-a : ℝ) : ℂ) - η * Complex.I)⁻¹‖ ≤ η / g ^ 2 := by
  have ha : 0 < a := lt_of_lt_of_le hg hga
  set A : ℂ := ((-a : ℝ) : ℂ) with hA
  set B : ℂ := A - η * Complex.I with hB
  have hnA : ‖A‖ = a := by
    rw [hA, Complex.norm_real, Real.norm_eq_abs, abs_neg, abs_of_pos ha]
  have hnB2 : ‖B‖ ^ 2 = a ^ 2 + η ^ 2 := norm_shift_sq
  have hnBnn : (0 : ℝ) ≤ ‖B‖ := norm_nonneg _
  have hnBge : a ≤ ‖B‖ := by nlinarith [sq_nonneg η]
  have hposB : (0 : ℝ) < ‖B‖ := lt_of_lt_of_le ha hnBge
  have hA0 : A ≠ 0 := by
    rw [hA, ne_eq, Complex.ofReal_eq_zero, neg_eq_zero]
    exact ha.ne'
  have hB0 : B ≠ 0 := norm_pos_iff.mp hposB
  have hBA : B - A = -(η * Complex.I) := by rw [hB]; ring
  have hnBA : ‖B - A‖ = η := by
    rw [hBA, norm_neg, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos hη]
  have hden : g ^ 2 ≤ a * ‖B‖ := by nlinarith
  rw [inv_sub_inv hA0 hB0, norm_div, norm_mul, hnBA, hnA]
  exact div_le_div_of_nonneg_left hη.le (by positivity) hden

/-- **Item T, second order.** `|(λ-x)⁻² - (λ-x-iη)⁻²| ≤ 2η / g³` when `x - λ ≥ g > 0`. -/
theorem norm_inv_sq_sub_inv_sq_le {a η g : ℝ} (hη : 0 < η) (hg : 0 < g) (hga : g ≤ a) :
    ‖(((-a : ℝ) : ℂ)⁻¹) ^ 2 - ((((-a : ℝ) : ℂ) - η * Complex.I)⁻¹) ^ 2‖ ≤ 2 * η / g ^ 3 := by
  have ha : 0 < a := lt_of_lt_of_le hg hga
  set A : ℂ := ((-a : ℝ) : ℂ) with hA
  set B : ℂ := A - η * Complex.I with hB
  have hnA : ‖A‖ = a := by
    rw [hA, Complex.norm_real, Real.norm_eq_abs, abs_neg, abs_of_pos ha]
  have hnB2 : ‖B‖ ^ 2 = a ^ 2 + η ^ 2 := norm_shift_sq
  have hnBnn : (0 : ℝ) ≤ ‖B‖ := norm_nonneg _
  have hnBge : a ≤ ‖B‖ := by nlinarith [sq_nonneg η]
  have hposB : (0 : ℝ) < ‖B‖ := lt_of_lt_of_le ha hnBge
  have hA0 : A ≠ 0 := by
    rw [hA, ne_eq, Complex.ofReal_eq_zero, neg_eq_zero]
    exact ha.ne'
  have hB0 : B ≠ 0 := norm_pos_iff.mp hposB
  have hBA : B - A = -(η * Complex.I) := by rw [hB]; ring
  have hnBA : ‖B - A‖ = η := by
    rw [hBA, norm_neg, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos hη]
  have hnS2 : ‖A + B‖ ^ 2 = 4 * a ^ 2 + η ^ 2 := norm_sum_shift_sq
  have hnSnn : (0 : ℝ) ≤ ‖A + B‖ := norm_nonneg _
  have hsqle : (a * ‖A + B‖) ^ 2 ≤ (2 * (a ^ 2 + η ^ 2)) ^ 2 := by
    have hexp : (a * ‖A + B‖) ^ 2 = a ^ 2 * ‖A + B‖ ^ 2 := by ring
    rw [hexp, hnS2]
    nlinarith [sq_nonneg η, sq_nonneg a, sq_nonneg (a * η)]
  have hXnn : (0 : ℝ) ≤ a * ‖A + B‖ := by positivity
  have hYnn : (0 : ℝ) ≤ 2 * (a ^ 2 + η ^ 2) := by positivity
  have key : a * ‖A + B‖ ≤ 2 * (a ^ 2 + η ^ 2) := by nlinarith [hsqle, hXnn, hYnn]
  have hfac : (A⁻¹) ^ 2 - (B⁻¹) ^ 2 = (A⁻¹ - B⁻¹) * (A⁻¹ + B⁻¹) := by ring
  rw [hfac, norm_mul, inv_sub_inv hA0 hB0, inv_add_inv hA0 hB0, norm_div, norm_div,
    norm_mul, hnBA, hnA, div_mul_div_comm,
    div_le_div_iff₀ (by positivity) (by positivity)]
  have hBB : a * ‖B‖ * (a * ‖B‖) = a ^ 2 * (a ^ 2 + η ^ 2) := by
    have hsq : ‖B‖ * ‖B‖ = a ^ 2 + η ^ 2 := by rw [← sq]; exact hnB2
    nlinarith [hsq]
  rw [hBB]
  have hg3 : g ^ 3 ≤ a ^ 3 := pow_le_pow_left₀ hg.le hga 3
  calc η * ‖A + B‖ * g ^ 3 ≤ η * ‖A + B‖ * a ^ 3 :=
        mul_le_mul_of_nonneg_left hg3 (by positivity)
    _ = η * a ^ 2 * (a * ‖A + B‖) := by ring
    _ ≤ η * a ^ 2 * (2 * (a ^ 2 + η ^ 2)) := mul_le_mul_of_nonneg_left key (by positivity)
    _ = 2 * η * (a ^ 2 * (a ^ 2 + η ^ 2)) := by ring

end R4C

end StackedSVD
