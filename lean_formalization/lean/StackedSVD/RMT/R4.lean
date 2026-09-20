/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Spectral

/-!
# R4: the deterministic finite-rank layer

Review note: `notes/archive/rmt_R4.md` (status `user OK` for proofs, 2026-08-29).
Nothing in this file is random and nothing is asymptotic. Item R5 turns it into a limit law.

Contents, in the order of the note:

* `resolv`, `qform`, `secular`: the resolvent `G₀ z = (W - z I)⁻¹`, the quadratic form
  `y ⬝ᵥ G₀ z y`, and the rank-one secular function `1 + q ⬝ᵥ G₀ z q`.
* R4a `det_sub_smul_one`: the finite-rank determinant identity.
* R4b `resolv_sub_resolv`, `hasDerivAt_resolv`: the resolvent identity and `G₀' = G₀²`.
* R4c `qform_nonpos`, `qform_strictMonoOn`, `hasDerivAt_qform`, `qform_tendsto_zero`,
  `posSemidef_sub_qform`.
* R4d `secular_strictMonoOn`, `secular_tendsto_one`, `secular_eq_zero_iff`,
  `eigenvalue_above_unique`, `topSpace_eq_span`, `lamMax_eq`, `topSimple`, `topProj_norm_sq`.
* R4e `lamMax_le_lamMax`, `opNorm_spiked_le`.

Everything rests on one master lemma, `resolv_eq_conj`:
`G₀ z = U diag((λ_a - z)⁻¹) Uᵀ`, with `U` the eigenvector matrix of `W`. The scalar facts then
come from the explicit sum `qform W z y = ∑ a, (λ_a - z)⁻¹ (Uᵀ y)_a ²`.
-/

open Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD
namespace R4

variable {d k : ℕ}

/-- `G₀ z = (W - z I)⁻¹`. Junk value `0` on the spectrum, because `Matrix.nonsing_inv` of a
singular matrix is `0`. -/
noncomputable def resolv (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) :
    Matrix (Fin d) (Fin d) ℝ :=
  (W - z • (1 : Matrix (Fin d) (Fin d) ℝ))⁻¹

/-- `Φ_y z = y ⬝ᵥ G₀ z y`. -/
noncomputable def qform (W : Matrix (Fin d) (Fin d) ℝ) (z : ℝ) (y : Fin d → ℝ) : ℝ :=
  y ⬝ᵥ (resolv W z *ᵥ y)

/-- The rank-one secular function `f z = 1 + q ⬝ᵥ G₀ z q`. -/
noncomputable def secular (W : Matrix (Fin d) (Fin d) ℝ) (q : Fin d → ℝ) (z : ℝ) : ℝ :=
  1 + qform W z q

variable {W₀ M : Matrix (Fin d) (Fin d) ℝ} {Q : Matrix (Fin d) (Fin k) ℝ}
  {q y : Fin d → ℝ} {z z₁ z₂ lam : ℝ}

/-! ### Spectral toolkit

`lamMax` of `Defs.lean` is `eigenvalues₀ 0`, the largest eigenvalue. These lemmas convert
between `eigenvalues₀`, `eigenvalues` and `lamMax`, and set up the orthogonal conjugation
`M = U diag(λ) Uᵀ` that every later proof uses.
-/

/-- Every eigenvalue is at most `lamMax`. -/
theorem eigenvalues_le_lamMax (hM : M.IsHermitian) (i : Fin d) :
    hM.eigenvalues i ≤ lamMax M hM := by
  have hd : 0 < d := i.pos
  rw [lamMax, dif_pos hd]
  exact hM.eigenvalues₀_antitone (Fin.le_def.mpr (Nat.zero_le _))

/-- `lamMax` is attained, provided the dimension is positive. -/
theorem exists_eigenvalues_eq_lamMax (hM : M.IsHermitian) (hd : 0 < d) :
    ∃ i, hM.eigenvalues i = lamMax M hM := by
  rw [lamMax, dif_pos hd]
  refine ⟨(Fintype.equivOfCardEq (Fintype.card_fin _)) ⟨0, by simpa using hd⟩, ?_⟩
  simp only [Matrix.IsHermitian.eigenvalues, Equiv.symm_apply_apply]

/-- The orthogonal matrix of eigenvectors of a real symmetric matrix. -/
noncomputable def eigU (hM : M.IsHermitian) : Matrix (Fin d) (Fin d) ℝ :=
  (hM.eigenvectorUnitary : Matrix (Fin d) (Fin d) ℝ)

theorem transpose_eigU_mul (hM : M.IsHermitian) : (eigU hM)ᵀ * eigU hM = 1 := by
  have h := hM.eigenvectorUnitary.2
  rw [Matrix.mem_unitaryGroup_iff'] at h
  rw [eigU, ← Matrix.conjTranspose_eq_transpose_of_trivial (α := ℝ)]
  exact h

theorem eigU_mul_transpose (hM : M.IsHermitian) : eigU hM * (eigU hM)ᵀ = 1 := by
  have h := hM.eigenvectorUnitary.2
  rw [Matrix.mem_unitaryGroup_iff] at h
  rw [eigU, ← Matrix.conjTranspose_eq_transpose_of_trivial (α := ℝ)]
  exact h

/-- Spectral theorem in the form this file uses. -/
theorem eigU_conj (hM : M.IsHermitian) :
    eigU hM * Matrix.diagonal hM.eigenvalues * (eigU hM)ᵀ = M := by
  conv_rhs => rw [hM.spectral_theorem]
  simp [eigU, Unitary.conjStarAlgAut_apply, Matrix.star_eq_conjTranspose,
    Matrix.conjTranspose_eq_transpose_of_trivial]

/-- Entry of a conjugated diagonal matrix. -/
theorem conj_apply (V : Matrix (Fin d) (Fin d) ℝ) (f : Fin d → ℝ) (i j : Fin d) :
    (V * Matrix.diagonal f * Vᵀ) i j = ∑ a, f a * (V i a * V j a) := by
  rw [Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, Matrix.transpose_apply]
  exact Finset.sum_congr rfl fun a _ => by ring

/-- Action of a conjugated diagonal matrix on a vector. -/
theorem conj_mulVec (V : Matrix (Fin d) (Fin d) ℝ) (f y : Fin d → ℝ) :
    (V * Matrix.diagonal f * Vᵀ) *ᵥ y = V *ᵥ (fun a => f a * (Vᵀ *ᵥ y) a) := by
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
  congr 1
  funext a
  simp [Matrix.mulVec_diagonal]

/-- Quadratic form of a conjugated diagonal matrix. -/
theorem dotProduct_conj (V : Matrix (Fin d) (Fin d) ℝ) (f y : Fin d → ℝ) :
    y ⬝ᵥ ((V * Matrix.diagonal f * Vᵀ) *ᵥ y) = ∑ a, f a * (Vᵀ *ᵥ y) a ^ 2 := by
  rw [conj_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]
  simp only [dotProduct]
  exact Finset.sum_congr rfl fun a _ => by ring

/-- An orthogonal change of coordinates preserves the squared length. -/
theorem dotProduct_transpose_eigU (hM : M.IsHermitian) (y : Fin d → ℝ) :
    ((eigU hM)ᵀ *ᵥ y) ⬝ᵥ ((eigU hM)ᵀ *ᵥ y) = y ⬝ᵥ y := by
  rw [Matrix.dotProduct_mulVec, Matrix.vecMul_transpose, Matrix.mulVec_mulVec,
    eigU_mul_transpose, Matrix.one_mulVec]

/-- The coordinates in the eigenbasis vanish only for the zero vector. -/
theorem eq_zero_of_transpose_eigU_mulVec (hM : M.IsHermitian) (h : (eigU hM)ᵀ *ᵥ y = 0) :
    y = 0 := by
  have hy : eigU hM *ᵥ ((eigU hM)ᵀ *ᵥ y) = y := by
    rw [Matrix.mulVec_mulVec, eigU_mul_transpose, Matrix.one_mulVec]
  rw [← hy, h, Matrix.mulVec_zero]

/-- Rayleigh upper bound: `x ⬝ᵥ M x ≤ lamMax M ‖x‖²`. -/
theorem dotProduct_mulVec_le_lamMax (hM : M.IsHermitian) (x : Fin d → ℝ) :
    x ⬝ᵥ (M *ᵥ x) ≤ lamMax M hM * (x ⬝ᵥ x) := by
  conv_lhs => rw [← eigU_conj hM]
  rw [dotProduct_conj, ← dotProduct_transpose_eigU hM x]
  simp only [dotProduct, Finset.mul_sum]
  refine Finset.sum_le_sum fun a _ => ?_
  have h1 : hM.eigenvalues a ≤ lamMax M hM := eigenvalues_le_lamMax hM a
  nlinarith [sq_nonneg (((eigU hM)ᵀ *ᵥ x) a)]

/-! ### The resolvent -/

theorem sub_smul_one_eq_conj (hM : M.IsHermitian) (z : ℝ) :
    M - z • (1 : Matrix (Fin d) (Fin d) ℝ)
      = eigU hM * Matrix.diagonal (fun a => hM.eigenvalues a - z) * (eigU hM)ᵀ := by
  have hdiag : Matrix.diagonal (fun a => hM.eigenvalues a - z)
      = Matrix.diagonal hM.eigenvalues - z • (1 : Matrix (Fin d) (Fin d) ℝ) := by
    ext i j
    rcases eq_or_ne i j with h | h
    · subst h; simp
    · simp [Matrix.diagonal_apply_ne _ h, Matrix.one_apply_ne h]
  rw [hdiag, Matrix.mul_sub, Matrix.sub_mul, eigU_conj]
  congr 1
  rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one, eigU_mul_transpose]

theorem det_sub_smul_one_eq_prod (hM : M.IsHermitian) (z : ℝ) :
    (M - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det = ∏ a, (hM.eigenvalues a - z) := by
  rw [sub_smul_one_eq_conj hM, Matrix.det_mul, Matrix.det_mul, Matrix.det_diagonal]
  have h : (eigU hM).det * ((eigU hM)ᵀ).det = 1 := by
    rw [← Matrix.det_mul, eigU_mul_transpose, Matrix.det_one]
  calc (eigU hM).det * (∏ a, (hM.eigenvalues a - z)) * ((eigU hM)ᵀ).det
      = ((eigU hM).det * ((eigU hM)ᵀ).det) * ∏ a, (hM.eigenvalues a - z) := by ring
    _ = ∏ a, (hM.eigenvalues a - z) := by rw [h, one_mul]

/-- Above the top eigenvalue every eigenvalue shift is strictly negative. -/
theorem eigenvalues_sub_neg (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (a : Fin d) :
    hW₀.eigenvalues a - z < 0 := by
  have := eigenvalues_le_lamMax hW₀ a
  linarith

/-- Above the top eigenvalue the shifted matrix is invertible. -/
theorem isUnit_det_sub (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    IsUnit (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det := by
  rw [isUnit_iff_ne_zero, det_sub_smul_one_eq_prod hW₀]
  exact Finset.prod_ne_zero_iff.mpr fun a _ => ne_of_lt (eigenvalues_sub_neg hW₀ hz a)

/-- **Master lemma.** `G₀ z = U diag((λ_a - z)⁻¹) Uᵀ` above the top eigenvalue. -/
theorem resolv_eq_conj (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    resolv W₀ z
      = eigU hW₀ * Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹) * (eigU hW₀)ᵀ := by
  have hne : ∀ a, hW₀.eigenvalues a - z ≠ 0 := fun a =>
    ne_of_lt (eigenvalues_sub_neg hW₀ hz a)
  have hDD : Matrix.diagonal (fun a => hW₀.eigenvalues a - z) *
      Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹)
        = (1 : Matrix (Fin d) (Fin d) ℝ) := by
    rw [Matrix.diagonal_mul_diagonal,
      show (fun a => (hW₀.eigenvalues a - z) * (hW₀.eigenvalues a - z)⁻¹) = fun _ => (1 : ℝ) from
        funext fun a => mul_inv_cancel₀ (hne a)]
    exact Matrix.diagonal_one
  refine Matrix.inv_eq_right_inv ?_
  rw [sub_smul_one_eq_conj hW₀]
  have hassoc : eigU hW₀ * Matrix.diagonal (fun a => hW₀.eigenvalues a - z) * (eigU hW₀)ᵀ *
      (eigU hW₀ * Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹) * (eigU hW₀)ᵀ)
      = eigU hW₀ * (Matrix.diagonal (fun a => hW₀.eigenvalues a - z) *
          ((eigU hW₀)ᵀ * eigU hW₀) *
          Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹)) * (eigU hW₀)ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, transpose_eigU_mul, Matrix.mul_one, hDD, Matrix.mul_one, eigU_mul_transpose]

theorem mul_resolv (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)) * resolv W₀ z = 1 :=
  Matrix.mul_nonsing_inv _ (isUnit_det_sub hW₀ hz)

theorem resolv_mul (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    resolv W₀ z * (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)) = 1 :=
  Matrix.nonsing_inv_mul _ (isUnit_det_sub hW₀ hz)

/-- The resolvent is symmetric, whatever `z` is (junk value `0` included). -/
theorem resolv_isHermitian (hW₀ : W₀.IsHermitian) : (resolv W₀ z).IsHermitian := by
  have h : (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ))ᴴ = W₀ - z • 1 := by
    rw [Matrix.conjTranspose_sub, Matrix.conjTranspose_smul, Matrix.conjTranspose_one, hW₀.eq]
    simp
  change (resolv W₀ z)ᴴ = resolv W₀ z
  rw [resolv, Matrix.conjTranspose_nonsing_inv, h]

theorem transpose_resolv (hW₀ : W₀.IsHermitian) : (resolv W₀ z)ᵀ = resolv W₀ z := by
  have h := (resolv_isHermitian (W₀ := W₀) (z := z) hW₀)
  rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- The square of the resolvent, in conjugated form. -/
theorem resolv_mul_resolv_eq_conj (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    resolv W₀ z * resolv W₀ z
      = eigU hW₀ * Matrix.diagonal (fun a => ((hW₀.eigenvalues a - z)⁻¹) ^ 2) * (eigU hW₀)ᵀ := by
  rw [resolv_eq_conj hW₀ hz]
  have hassoc : eigU hW₀ * Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹) * (eigU hW₀)ᵀ *
      (eigU hW₀ * Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹) * (eigU hW₀)ᵀ)
      = eigU hW₀ * (Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹) *
          ((eigU hW₀)ᵀ * eigU hW₀) *
          Matrix.diagonal (fun a => (hW₀.eigenvalues a - z)⁻¹)) * (eigU hW₀)ᵀ := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, transpose_eigU_mul, Matrix.mul_one, Matrix.diagonal_mul_diagonal,
    show (fun a => (hW₀.eigenvalues a - z)⁻¹ * (hW₀.eigenvalues a - z)⁻¹)
      = fun a => ((hW₀.eigenvalues a - z)⁻¹) ^ 2 from funext fun a => (sq _).symm]

/-! ### R4a: the determinant identity -/

/-- **R4a.** The finite-rank determinant identity. -/
theorem det_sub_smul_one (hz : IsUnit (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det) :
    (W₀ + Q * Qᵀ - z • 1).det
      = (W₀ - z • 1).det * (1 + Qᵀ * resolv W₀ z * Q).det := by
  have h : W₀ + Q * Qᵀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)
      = (W₀ - z • 1) + Q * Qᵀ := by abel
  rw [h, Matrix.det_add_mul _ _ hz]
  rfl

/-! ### R4b: the resolvent identity and `G₀' = G₀²` -/

/-- **R4b.** `G₀ z₂ - G₀ z₁ = (z₂ - z₁) G₀ z₁ G₀ z₂`. -/
theorem resolv_sub_resolv
    (h₁ : IsUnit (W₀ - z₁ • (1 : Matrix (Fin d) (Fin d) ℝ)).det)
    (h₂ : IsUnit (W₀ - z₂ • (1 : Matrix (Fin d) (Fin d) ℝ)).det) :
    resolv W₀ z₂ - resolv W₀ z₁ = (z₂ - z₁) • (resolv W₀ z₁ * resolv W₀ z₂) := by
  have e1 : resolv W₀ z₁ * (W₀ - z₁ • (1 : Matrix (Fin d) (Fin d) ℝ)) = 1 :=
    Matrix.nonsing_inv_mul _ h₁
  have e2 : (W₀ - z₂ • (1 : Matrix (Fin d) (Fin d) ℝ)) * resolv W₀ z₂ = 1 :=
    Matrix.mul_nonsing_inv _ h₂
  have key : resolv W₀ z₁ * ((W₀ - z₁ • (1 : Matrix (Fin d) (Fin d) ℝ))
      - (W₀ - z₂ • (1 : Matrix (Fin d) (Fin d) ℝ))) * resolv W₀ z₂
      = resolv W₀ z₂ - resolv W₀ z₁ := by
    rw [Matrix.mul_sub, Matrix.sub_mul, e1, Matrix.one_mul, Matrix.mul_assoc, e2, Matrix.mul_one]
  rw [← key]
  have hsub : (W₀ - z₁ • (1 : Matrix (Fin d) (Fin d) ℝ))
      - (W₀ - z₂ • (1 : Matrix (Fin d) (Fin d) ℝ))
      = (z₂ - z₁) • (1 : Matrix (Fin d) (Fin d) ℝ) := by
    rw [sub_smul]
    abel
  rw [hsub, Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul]

/-- Entry of the resolvent, as an explicit real function of `z`. -/
theorem resolv_apply_eq_sum (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (i j : Fin d) :
    resolv W₀ z i j
      = ∑ a, (hW₀.eigenvalues a - z)⁻¹ * (eigU hW₀ i a * eigU hW₀ j a) := by
  rw [resolv_eq_conj hW₀ hz, conj_apply]

/-- `d/dz (c - z)⁻¹ = ((c - z)⁻¹)²`. -/
theorem hasDerivAt_inv_sub {c : ℝ} (hc : c - z ≠ 0) :
    HasDerivAt (fun t : ℝ => (c - t)⁻¹) (((c - z)⁻¹) ^ 2) z := by
  have h1 : HasDerivAt (fun t : ℝ => c - t) (-1) z := by
    simpa using (hasDerivAt_id z).const_sub c
  have h2 := h1.inv hc
  have h3 : -(-1 : ℝ) / (c - z) ^ 2 = ((c - z)⁻¹) ^ 2 := by
    rw [neg_neg, one_div, inv_pow]
  rwa [h3] at h2

/-- **R4b, derivative form.** `G₀' = G₀²`. -/
theorem hasDerivAt_resolv (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    HasDerivAt (resolv W₀) (resolv W₀ z * resolv W₀ z) z := by
  have hne : ∀ a, hW₀.eigenvalues a - z ≠ 0 := fun a =>
    ne_of_lt (eigenvalues_sub_neg hW₀ hz a)
  refine hasDerivAt_pi.mpr fun i => hasDerivAt_pi.mpr fun j => ?_
  have hev : (fun t => resolv W₀ t i j)
      =ᶠ[𝓝 z] fun t => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * (eigU hW₀ i a * eigU hW₀ j a) := by
    filter_upwards [Ioi_mem_nhds hz] with t ht
    exact resolv_apply_eq_sum hW₀ ht i j
  have hderiv : HasDerivAt
      (fun t => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * (eigU hW₀ i a * eigU hW₀ j a))
      (∑ a, ((hW₀.eigenvalues a - z)⁻¹) ^ 2 * (eigU hW₀ i a * eigU hW₀ j a)) z :=
    HasDerivAt.fun_sum fun a _ =>
      (hasDerivAt_inv_sub (c := hW₀.eigenvalues a) (z := z) (hne a)).mul_const _
  have hval : (resolv W₀ z * resolv W₀ z) i j
      = ∑ a, ((hW₀.eigenvalues a - z)⁻¹) ^ 2 * (eigU hW₀ i a * eigU hW₀ j a) := by
    rw [resolv_mul_resolv_eq_conj hW₀ hz, conj_apply]
  rw [hval]
  exact hderiv.congr_of_eventuallyEq hev

/-! ### R4c: the monotone quadratic form -/

/-- `qform` in the eigenbasis: `Φ_y z = ∑ a (λ_a - z)⁻¹ c_a²`, `c = Uᵀ y`. -/
theorem qform_eq_sum (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (y : Fin d → ℝ) :
    qform W₀ z y = ∑ a, (hW₀.eigenvalues a - z)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 := by
  rw [qform, resolv_eq_conj hW₀ hz, dotProduct_conj]

private theorem inv_sub_lt_inv_sub {u a b : ℝ} (hua : u < a) (hab : a < b) :
    (u - a)⁻¹ < (u - b)⁻¹ := by
  have h1 : u - a < 0 := by linarith
  have h2 : u - b < 0 := by linarith
  have h1' : u - a ≠ 0 := ne_of_lt h1
  have h2' : u - b ≠ 0 := ne_of_lt h2
  have hprod : 0 < (u - a) * (u - b) := mul_pos_of_neg_of_neg h1 h2
  rw [← sub_pos]
  have hrw : (u - b)⁻¹ - (u - a)⁻¹ = (b - a) / ((u - a) * (u - b)) := by
    field_simp
    ring
  rw [hrw]
  exact div_pos (by linarith) hprod

private theorem inv_sub_le_inv_sub {u a b : ℝ} (hua : u < a) (hab : a ≤ b) :
    (u - a)⁻¹ ≤ (u - b)⁻¹ := by
  rcases eq_or_lt_of_le hab with h | h
  · rw [h]
  · exact (inv_sub_lt_inv_sub hua h).le

/-- **R4c.1.** `Φ_y ≤ 0` above the top eigenvalue. -/
theorem qform_nonpos (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    qform W₀ z y ≤ 0 := by
  rw [qform_eq_sum hW₀ hz]
  refine Finset.sum_nonpos fun a _ => ?_
  have h1 : hW₀.eigenvalues a - z < 0 := eigenvalues_sub_neg hW₀ hz a
  have h2 : (hW₀.eigenvalues a - z)⁻¹ < 0 := inv_neg''.mpr h1
  nlinarith [sq_nonneg (((eigU hW₀)ᵀ *ᵥ y) a)]

/-- **R4c.1, strict form.** `Φ_y < 0` above the top eigenvalue when `y ≠ 0`. -/
theorem qform_neg (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (hy : y ≠ 0) :
    qform W₀ z y < 0 := by
  obtain ⟨j, hj⟩ : ∃ j, ((eigU hW₀)ᵀ *ᵥ y) j ≠ 0 := by
    by_contra hcon
    push Not at hcon
    exact hy (eq_zero_of_transpose_eigU_mulVec hW₀ (funext hcon))
  have hstep : ∀ a : Fin d, (hW₀.eigenvalues a - z)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 ≤ 0 := by
    intro a
    have h2 : (hW₀.eigenvalues a - z)⁻¹ < 0 := inv_neg''.mpr (eigenvalues_sub_neg hW₀ hz a)
    nlinarith [sq_nonneg (((eigU hW₀)ᵀ *ᵥ y) a)]
  have hlt : ∑ a, (hW₀.eigenvalues a - z)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2
      < ∑ _a : Fin d, (0 : ℝ) := by
    refine Finset.sum_lt_sum (fun i _ => hstep i) ⟨j, Finset.mem_univ j, ?_⟩
    have h2 : (hW₀.eigenvalues j - z)⁻¹ < 0 := inv_neg''.mpr (eigenvalues_sub_neg hW₀ hz j)
    have hpos : 0 < ((eigU hW₀)ᵀ *ᵥ y) j ^ 2 := by positivity
    nlinarith
  rw [qform_eq_sum hW₀ hz]
  simpa using hlt

/-- **R4c.** `-G₀ z` is positive definite above the top eigenvalue. -/
theorem posDef_neg_resolv (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    (-resolv W₀ z).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos (resolv_isHermitian hW₀).neg fun x hx => ?_
  have hstar : star x = x := rfl
  rw [hstar, Matrix.neg_mulVec, dotProduct_neg, neg_pos]
  exact qform_neg hW₀ hz hx

/-- Monotonicity of `Φ_y`, non-strict version. -/
theorem qform_le_qform (hW₀ : W₀.IsHermitian) (h₁ : lamMax W₀ hW₀ < z₁) (h₂ : z₁ ≤ z₂) :
    qform W₀ z₁ y ≤ qform W₀ z₂ y := by
  have h₂' : lamMax W₀ hW₀ < z₂ := lt_of_lt_of_le h₁ h₂
  rw [qform_eq_sum hW₀ h₁, qform_eq_sum hW₀ h₂']
  refine Finset.sum_le_sum fun a _ => ?_
  have hle : hW₀.eigenvalues a ≤ lamMax W₀ hW₀ := eigenvalues_le_lamMax hW₀ a
  have hmono := inv_sub_le_inv_sub (u := hW₀.eigenvalues a) (a := z₁) (b := z₂)
    (lt_of_le_of_lt hle h₁) h₂
  nlinarith [sq_nonneg (((eigU hW₀)ᵀ *ᵥ y) a)]

/-- **R4c.2.** `Φ_y` is strictly increasing above the top eigenvalue when `y ≠ 0`. -/
theorem qform_strictMonoOn (hW₀ : W₀.IsHermitian) (hy : y ≠ 0) :
    StrictMonoOn (fun t => qform W₀ t y) (Set.Ioi (lamMax W₀ hW₀)) := by
  intro a ha b hb hab
  simp only [Set.mem_Ioi] at ha hb
  obtain ⟨j, hj⟩ : ∃ j, ((eigU hW₀)ᵀ *ᵥ y) j ≠ 0 := by
    by_contra hcon
    push Not at hcon
    exact hy (eq_zero_of_transpose_eigU_mulVec hW₀ (funext hcon))
  simp only
  rw [qform_eq_sum hW₀ ha, qform_eq_sum hW₀ hb]
  refine Finset.sum_lt_sum (fun i _ => ?_) ⟨j, Finset.mem_univ j, ?_⟩
  · have hle : hW₀.eigenvalues i ≤ lamMax W₀ hW₀ := eigenvalues_le_lamMax hW₀ i
    have hmono := inv_sub_le_inv_sub (u := hW₀.eigenvalues i) (a := a) (b := b)
      (lt_of_le_of_lt hle ha) hab.le
    nlinarith [sq_nonneg (((eigU hW₀)ᵀ *ᵥ y) i)]
  · have hle : hW₀.eigenvalues j ≤ lamMax W₀ hW₀ := eigenvalues_le_lamMax hW₀ j
    have hlt := inv_sub_lt_inv_sub (u := hW₀.eigenvalues j) (a := a) (b := b)
      (lt_of_le_of_lt hle ha) hab
    have hpos : 0 < ((eigU hW₀)ᵀ *ᵥ y) j ^ 2 := by positivity
    exact mul_lt_mul_of_pos_right hlt hpos

/-- **R4c.3.** `d/dz Φ_y z = y ⬝ᵥ G₀ z² y`. -/
theorem hasDerivAt_qform (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    HasDerivAt (fun t => qform W₀ t y) (y ⬝ᵥ ((resolv W₀ z * resolv W₀ z) *ᵥ y)) z := by
  have hne : ∀ a, hW₀.eigenvalues a - z ≠ 0 := fun a =>
    ne_of_lt (eigenvalues_sub_neg hW₀ hz a)
  have hev : (fun t => qform W₀ t y)
      =ᶠ[𝓝 z] fun t => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 := by
    filter_upwards [Ioi_mem_nhds hz] with t ht
    exact qform_eq_sum hW₀ ht y
  have hderiv : HasDerivAt
      (fun t => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2)
      (∑ a, ((hW₀.eigenvalues a - z)⁻¹) ^ 2 * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2) z :=
    HasDerivAt.fun_sum fun a _ =>
      (hasDerivAt_inv_sub (c := hW₀.eigenvalues a) (z := z) (hne a)).mul_const _
  have hval : y ⬝ᵥ ((resolv W₀ z * resolv W₀ z) *ᵥ y)
      = ∑ a, ((hW₀.eigenvalues a - z)⁻¹) ^ 2 * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 := by
    rw [resolv_mul_resolv_eq_conj hW₀ hz, dotProduct_conj]
  rw [hval]
  exact hderiv.congr_of_eventuallyEq hev

/-- The derivative of `Φ_y` is the squared length of `G₀ z y`; in particular it is `≥ 0`. -/
theorem dotProduct_resolv_sq_eq_norm_sq (hW₀ : W₀.IsHermitian) :
    y ⬝ᵥ ((resolv W₀ z * resolv W₀ z) *ᵥ y)
      = ‖(WithLp.toLp 2 (resolv W₀ z *ᵥ y) : EuclideanSpace ℝ (Fin d))‖ ^ 2 := by
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    transpose_resolv hW₀, EuclideanSpace.norm_eq, Real.sq_sqrt (by positivity)]
  simp [dotProduct, sq]

/-- **R4c.4.** `Φ_y z → 0` as `z → ∞`. -/
theorem qform_tendsto_zero (hW₀ : W₀.IsHermitian) :
    Tendsto (fun t => qform W₀ t y) atTop (𝓝 0) := by
  have hev : (fun t => qform W₀ t y)
      =ᶠ[atTop] fun t => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2 := by
    filter_upwards [eventually_gt_atTop (lamMax W₀ hW₀)] with t ht
    exact qform_eq_sum hW₀ ht y
  refine Tendsto.congr' hev.symm ?_
  have hsum : Tendsto (fun t : ℝ => ∑ a, (hW₀.eigenvalues a - t)⁻¹ * ((eigU hW₀)ᵀ *ᵥ y) a ^ 2)
      atTop (𝓝 (∑ _a : Fin d, (0 : ℝ))) := by
    refine tendsto_finsetSum _ fun a _ => ?_
    have h1 : Tendsto (fun t : ℝ => hW₀.eigenvalues a - t) atTop atBot :=
      tendsto_atBot_add_const_left _ _ tendsto_neg_atTop_atBot
    have h2 : Tendsto (fun t : ℝ => (hW₀.eigenvalues a - t)⁻¹) atTop (𝓝 0) :=
      h1.inv_tendsto_atBot
    simpa using h2.mul_const (((eigU hW₀)ᵀ *ᵥ y) a ^ 2)
  simpa using hsum

/-- **R4c, Loewner form (order-free).** `z ↦ Qᵀ G₀ z Q` increases. -/
theorem posSemidef_sub_qform (hW₀ : W₀.IsHermitian)
    (h₁ : lamMax W₀ hW₀ < z₁) (h₂ : z₁ < z₂) :
    (Qᵀ * resolv W₀ z₂ * Q - Qᵀ * resolv W₀ z₁ * Q).PosSemidef := by
  have hHerm : ∀ w : ℝ, (Qᵀ * resolv W₀ w * Q).IsHermitian := by
    intro w
    have h := Matrix.isHermitian_conjTranspose_mul_mul (A := resolv W₀ w) Q
      (resolv_isHermitian hW₀)
    rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have hquad : ∀ (w : ℝ) (x : Fin k → ℝ),
      x ⬝ᵥ ((Qᵀ * resolv W₀ w * Q) *ᵥ x) = qform W₀ w (Q *ᵥ x) := by
    intro w x
    rw [qform, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
      Matrix.vecMul_transpose]
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ((hHerm z₂).sub (hHerm z₁)) fun x => ?_
  have hstar : star x = x := rfl
  rw [hstar, Matrix.sub_mulVec, dotProduct_sub, hquad, hquad, sub_nonneg]
  exact qform_le_qform hW₀ h₁ h₂.le

/-! ### R4d: the rank-one secular function -/

/-- **R4d.1.** -/
theorem secular_strictMonoOn (hW₀ : W₀.IsHermitian) (hq : q ≠ 0) :
    StrictMonoOn (secular W₀ q) (Set.Ioi (lamMax W₀ hW₀)) := by
  intro a ha b hb hab
  have h := qform_strictMonoOn hW₀ hq ha hb hab
  change 1 + qform W₀ a q < 1 + qform W₀ b q
  simp only at h
  linarith

/-- **R4d.2.** -/
theorem secular_tendsto_one (hW₀ : W₀.IsHermitian) :
    Tendsto (secular W₀ q) atTop (𝓝 1) := by
  have h := (qform_tendsto_zero (W₀ := W₀) (y := q) hW₀).const_add 1
  change Tendsto (fun t => 1 + qform W₀ t q) atTop (𝓝 1)
  simpa using h

/-- `vecMulVec a b` acts as `x ↦ (b ⬝ᵥ x) • a`. -/
theorem vecMulVec_mulVec {m : ℕ} (a : Fin m → ℝ) (b x : Fin d → ℝ) :
    Matrix.vecMulVec a b *ᵥ x = (b ⬝ᵥ x) • a := by
  funext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- `vecMulVec q q = Q Qᵀ` for the `d × 1` matrix `Q` with column `q`. -/
theorem vecMulVec_eq_replicateCol (q : Fin d → ℝ) :
    Matrix.vecMulVec q q
      = Matrix.replicateCol (Fin 1) q * (Matrix.replicateCol (Fin 1) q)ᵀ := by
  ext i j
  simp [Matrix.vecMulVec_apply, Matrix.mul_apply]

/-- The `1 × 1` determinant of R4a is the secular function. -/
theorem det_one_add_col (G : Matrix (Fin d) (Fin d) ℝ) (q : Fin d → ℝ) :
    ((1 : Matrix (Fin 1) (Fin 1) ℝ) + (Matrix.replicateCol (Fin 1) q)ᵀ * G *
        Matrix.replicateCol (Fin 1) q).det
      = 1 + q ⬝ᵥ (G *ᵥ q) := by
  rw [Matrix.det_unique]
  simp only [Matrix.add_apply, Matrix.one_apply_eq, Matrix.mul_apply, Matrix.transpose_apply,
    Matrix.replicateCol_apply, dotProduct, Matrix.mulVec, Finset.sum_mul, Finset.mul_sum]
  refine congrArg _ ?_
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

set_option linter.deprecated false in
/-- The spectrum of `toOp N` is the spectrum of `N`. -/
theorem spectrum_toOp (N : Matrix (Fin d) (Fin d) ℝ) :
    spectrum ℝ (toOp N) = spectrum ℝ N :=
  Matrix.spectrum_toEuclideanLin

/-- **R4d.3.** `z` above `lamMax W₀` is an eigenvalue of `W₀ + q qᵀ` iff `f z = 0`. -/
theorem secular_eq_zero_iff (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    secular W₀ q z = 0 ↔ z ∈ spectrum ℝ (toOp (W₀ + Matrix.vecMulVec q q)) := by
  have hdet : (W₀ + Matrix.vecMulVec q q - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det
      = (W₀ - z • 1).det * secular W₀ q z := by
    rw [vecMulVec_eq_replicateCol,
      det_sub_smul_one (Q := Matrix.replicateCol (Fin 1) q) (isUnit_det_sub hW₀ hz),
      det_one_add_col]
    rfl
  have hne : (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det ≠ 0 :=
    isUnit_iff_ne_zero.mp (isUnit_det_sub hW₀ hz)
  rw [spectrum_toOp, Matrix.mem_spectrum_iff_not_isUnit_eval_charpoly, Matrix.eval_charpoly,
    isUnit_iff_ne_zero, not_ne_iff]
  have hscalar : (Matrix.scalar (Fin d)) z = z • (1 : Matrix (Fin d) (Fin d) ℝ) := by
    rw [Matrix.smul_one_eq_diagonal]
    rfl
  have hneg : (z • (1 : Matrix (Fin d) (Fin d) ℝ) - (W₀ + Matrix.vecMulVec q q))
      = -(W₀ + Matrix.vecMulVec q q - z • 1) := by abel
  rw [hscalar, hneg, Matrix.det_neg, hdet]
  constructor
  · intro h
    simp [h]
  · intro h
    have h' : (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)).det * secular W₀ q z = 0 := by
      rcases mul_eq_zero.mp h with h'' | h''
      · exact absurd h'' (pow_ne_zero _ (by norm_num))
      · exact h''
    rcases mul_eq_zero.mp h' with h'' | h''
    · exact absurd h'' hne
    · exact h''

/-- **R4d.4.** At most one eigenvalue above `lamMax W₀`. -/
theorem eigenvalue_above_unique (hW₀ : W₀.IsHermitian) (hq : q ≠ 0)
    (h₁ : lamMax W₀ hW₀ < z₁) (h₂ : lamMax W₀ hW₀ < z₂)
    (e₁ : secular W₀ q z₁ = 0) (e₂ : secular W₀ q z₂ = 0) : z₁ = z₂ :=
  (secular_strictMonoOn hW₀ hq).injOn h₁ h₂ (e₁.trans e₂.symm)

theorem resolv_mulVec_ne_zero (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) (hq : q ≠ 0) :
    resolv W₀ z *ᵥ q ≠ 0 := by
  intro h
  apply hq
  have hbase : (W₀ - z • (1 : Matrix (Fin d) (Fin d) ℝ)) *ᵥ (resolv W₀ z *ᵥ q) = q := by
    rw [Matrix.mulVec_mulVec, mul_resolv hW₀ hz, Matrix.one_mulVec]
  rw [h, Matrix.mulVec_zero] at hbase
  exact hbase.symm

set_option linter.deprecated false in
/-- Membership in an eigenspace of `toOp`, in matrix language. -/
theorem mem_eigenspace_iff' (N : Matrix (Fin d) (Fin d) ℝ) (t : ℝ)
    (x : EuclideanSpace ℝ (Fin d)) :
    x ∈ Module.End.eigenspace (toOp N) t ↔ N *ᵥ (WithLp.ofLp x) = t • (WithLp.ofLp x) := by
  rw [Module.End.mem_eigenspace_iff]
  constructor
  · intro h
    have h' := congrArg WithLp.ofLp h
    simpa [Matrix.toEuclideanLin_apply] using h'
  · intro h
    apply WithLp.ofLp_injective
    simpa [Matrix.toEuclideanLin_apply] using h

set_option linter.unusedVariables false in
/-- **R4d.5 and R4d.6.** The eigenspace at a secular root is the line through `G₀(λ) q`.
`hq` is kept for interface stability with `notes/archive/rmt_R4.md`; the proof does not use it,
because `secular W₀ 0 lam = 1 ≠ 0` already makes `e` unsatisfiable when `q = 0`. -/
theorem eigenspace_eq_span (hW₀ : W₀.IsHermitian) (hq : q ≠ 0)
    (hlam : lamMax W₀ hW₀ < lam) (e : secular W₀ q lam = 0) :
    Module.End.eigenspace (toOp (W₀ + Matrix.vecMulVec q q)) lam
      = Submodule.span ℝ
          {(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d))} := by
  have hqx : q ⬝ᵥ (resolv W₀ lam *ᵥ q) = -1 := by
    have hs : secular W₀ q lam = 1 + q ⬝ᵥ (resolv W₀ lam *ᵥ q) := rfl
    rw [hs] at e
    linarith
  have hinv : ∀ w : Fin d → ℝ,
      resolv W₀ lam *ᵥ ((W₀ - lam • (1 : Matrix (Fin d) (Fin d) ℝ)) *ᵥ w) = w := by
    intro w
    rw [Matrix.mulVec_mulVec, resolv_mul hW₀ hlam, Matrix.one_mulVec]
  have hbase : (W₀ - lam • (1 : Matrix (Fin d) (Fin d) ℝ)) *ᵥ (resolv W₀ lam *ᵥ q) = q := by
    rw [Matrix.mulVec_mulVec, mul_resolv hW₀ hlam, Matrix.one_mulVec]
  have hfwd : ∀ w : Fin d → ℝ,
      (W₀ + Matrix.vecMulVec q q) *ᵥ w = lam • w
        ↔ (W₀ - lam • (1 : Matrix (Fin d) (Fin d) ℝ)) *ᵥ w = -((q ⬝ᵥ w) • q) := by
    intro w
    rw [Matrix.add_mulVec, vecMulVec_mulVec, Matrix.sub_mulVec, Matrix.smul_mulVec,
      Matrix.one_mulVec]
    constructor
    · intro h
      rw [← h]
      abel
    · intro h
      rw [sub_eq_iff_eq_add] at h
      rw [h]
      abel
  ext x
  rw [mem_eigenspace_iff', Submodule.mem_span_singleton]
  constructor
  · intro hx
    have h1 := (hfwd (WithLp.ofLp x)).mp hx
    have h2 := hinv (WithLp.ofLp x)
    rw [h1, Matrix.mulVec_neg, Matrix.mulVec_smul] at h2
    refine ⟨-(q ⬝ᵥ WithLp.ofLp x), ?_⟩
    apply WithLp.ofLp_injective
    simpa using h2
  · rintro ⟨c, rfl⟩
    have hcx : WithLp.ofLp
        (c • (WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)))
        = c • (resolv W₀ lam *ᵥ q) := by
      simp
    rw [hcx, hfwd, Matrix.mulVec_smul, hbase, dotProduct_smul, hqx]
    simp

/-- **R4d.4 to R4d.6.** With a secular root above `lamMax W₀`, that root is `lamMax` of `A`. -/
theorem lamMax_eq (hW₀ : W₀.IsHermitian) (hq : q ≠ 0) (hlam : lamMax W₀ hW₀ < lam)
    (e : secular W₀ q lam = 0) (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    lamMax (W₀ + Matrix.vecMulVec q q) hA = lam := by
  have hd : 0 < d := by
    rcases Nat.eq_zero_or_pos d with h | h
    · exact absurd (funext fun i => absurd (h ▸ i.pos) (by omega)) hq
    · exact h
  have hmem : lam ∈ spectrum ℝ (W₀ + Matrix.vecMulVec q q) := by
    have h := (secular_eq_zero_iff (q := q) hW₀ hlam).mp e
    rwa [spectrum_toOp] at h
  obtain ⟨i, hi⟩ : ∃ i, hA.eigenvalues i = lam := by
    rw [hA.spectrum_real_eq_range_eigenvalues] at hmem
    exact hmem
  have hle : lam ≤ lamMax (W₀ + Matrix.vecMulVec q q) hA := hi ▸ eigenvalues_le_lamMax hA i
  obtain ⟨j, hj⟩ := exists_eigenvalues_eq_lamMax hA hd
  have hjspec : lamMax (W₀ + Matrix.vecMulVec q q) hA ∈ spectrum ℝ (W₀ + Matrix.vecMulVec q q) :=
    hj ▸ hA.eigenvalues_mem_spectrum_real j
  have hjop : lamMax (W₀ + Matrix.vecMulVec q q) hA
      ∈ spectrum ℝ (toOp (W₀ + Matrix.vecMulVec q q)) := by
    rwa [spectrum_toOp]
  have hjlt : lamMax W₀ hW₀ < lamMax (W₀ + Matrix.vecMulVec q q) hA := lt_of_lt_of_le hlam hle
  have hjsec : secular W₀ q (lamMax (W₀ + Matrix.vecMulVec q q) hA) = 0 :=
    (secular_eq_zero_iff hW₀ hjlt).mpr hjop
  exact eigenvalue_above_unique hW₀ hq hjlt hlam hjsec e

/-- **R4d.5 and R4d.6.** -/
theorem topSpace_eq_span (hW₀ : W₀.IsHermitian) (hq : q ≠ 0)
    (hlam : lamMax W₀ hW₀ < lam) (e : secular W₀ q lam = 0)
    (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    topSpace (W₀ + Matrix.vecMulVec q q) hA
      = Submodule.span ℝ
          {(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d))} := by
  rw [topSpace, specSpace, lamMax_eq hW₀ hq hlam e hA]
  simp only [Set.mem_singleton_iff, iSup_iSup_eq_left]
  exact eigenspace_eq_span hW₀ hq hlam e

theorem resolv_toLp_ne_zero (hW₀ : W₀.IsHermitian) (hlam : lamMax W₀ hW₀ < lam) (hq : q ≠ 0) :
    (WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)) ≠ 0 := by
  intro h
  exact resolv_mulVec_ne_zero hW₀ hlam hq (by simpa using congrArg WithLp.ofLp h)

/-- **R4d.6.** The top eigenvalue is simple. -/
theorem topSimple (hW₀ : W₀.IsHermitian) (hq : q ≠ 0) (hlam : lamMax W₀ hW₀ < lam)
    (e : secular W₀ q lam = 0) (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    TopSimple (W₀ + Matrix.vecMulVec q q) hA := by
  rw [TopSimple, topSpace_eq_span hW₀ hq hlam e hA]
  exact finrank_span_singleton (resolv_toLp_ne_zero hW₀ hlam hq)

/-- **R4d.7.** The overlap formula. `overlap X w` of `Defs.lean` is `‖topProj (Xᵀ X) w‖²`, so
this is stated on `topProj` of the abstract `A`. -/
theorem topProj_norm_sq (hW₀ : W₀.IsHermitian) (hq : q ≠ 0) (hlam : lamMax W₀ hW₀ < lam)
    (e : secular W₀ q lam = 0) (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian)
    (w : EuclideanSpace ℝ (Fin d)) :
    ‖topProj (W₀ + Matrix.vecMulVec q q) hA w‖ ^ 2
      = (WithLp.ofLp w ⬝ᵥ (resolv W₀ lam *ᵥ q)) ^ 2
        / (q ⬝ᵥ ((resolv W₀ lam * resolv W₀ lam) *ᵥ q)) := by
  have hx : (WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)) ≠ 0 :=
    resolv_toLp_ne_zero hW₀ hlam hq
  have hnx : ‖(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d))‖ ≠ 0 :=
    norm_ne_zero_iff.mpr hx
  have hproj0 : topProj (W₀ + Matrix.vecMulVec q q) hA
      = (topSpace (W₀ + Matrix.vecMulVec q q) hA).starProjection := rfl
  have hproj : topProj (W₀ + Matrix.vecMulVec q q) hA w
      = (⟪(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)), w⟫_ℝ /
          ‖(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d))‖ ^ 2) •
        (WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)) := by
    rw [hproj0, topSpace_eq_span hW₀ hq hlam e hA]
    exact Submodule.starProjection_singleton ℝ w
  have hden : ‖(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d))‖ ^ 2
      = q ⬝ᵥ ((resolv W₀ lam * resolv W₀ lam) *ᵥ q) :=
    (dotProduct_resolv_sq_eq_norm_sq hW₀).symm
  have hnum : ⟪(WithLp.toLp 2 (resolv W₀ lam *ᵥ q) : EuclideanSpace ℝ (Fin d)), w⟫_ℝ
      = WithLp.ofLp w ⬝ᵥ (resolv W₀ lam *ᵥ q) := by
    rw [EuclideanSpace.inner_eq_star_dotProduct]
    simp
  rw [hproj, norm_smul, mul_pow, ← hden, hnum, Real.norm_eq_abs, sq_abs]
  field_simp

/-! ### R4e: crude bounds -/

/-- **R4e.1, general form.** A perturbation with a nonnegative quadratic form does not lower
the top eigenvalue. -/
theorem lamMax_le_lamMax_of_nonneg (hW₀ : W₀.IsHermitian) (P : Matrix (Fin d) (Fin d) ℝ)
    (hP : ∀ x : Fin d → ℝ, 0 ≤ x ⬝ᵥ (P *ᵥ x)) (hA : (W₀ + P).IsHermitian) :
    lamMax W₀ hW₀ ≤ lamMax (W₀ + P) hA := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    simp [lamMax]
  obtain ⟨j, hj⟩ := exists_eigenvalues_eq_lamMax hW₀ hd
  set y : Fin d → ℝ := eigU hW₀ *ᵥ Pi.single j 1 with hy
  have hcoord : (eigU hW₀)ᵀ *ᵥ y = Pi.single j 1 := by
    rw [hy, Matrix.mulVec_mulVec, transpose_eigU_mul, Matrix.one_mulVec]
  have hyy : y ⬝ᵥ y = 1 := by
    rw [← dotProduct_transpose_eigU hW₀ y, hcoord]
    simp [dotProduct, Pi.single_apply]
  have hquad : y ⬝ᵥ (W₀ *ᵥ y) = lamMax W₀ hW₀ := by
    conv_lhs => rw [← eigU_conj hW₀]
    rw [dotProduct_conj, hcoord]
    simp [Pi.single_apply, hj]
  have hupper := dotProduct_mulVec_le_lamMax hA y
  rw [Matrix.add_mulVec, dotProduct_add, hquad, hyy, mul_one] at hupper
  linarith [hP y]

/-- **R4e.1.** `lamMax W₀ ≤ lamMax (W₀ + Q Qᵀ)`. -/
theorem lamMax_le_lamMax (hW₀ : W₀.IsHermitian) (hA : (W₀ + Q * Qᵀ).IsHermitian) :
    lamMax W₀ hW₀ ≤ lamMax (W₀ + Q * Qᵀ) hA := by
  refine lamMax_le_lamMax_of_nonneg hW₀ _ (fun x => ?_) hA
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]
  simp only [dotProduct]
  exact Finset.sum_nonneg fun a _ => mul_self_nonneg _

/-- **R4e.1, rank one.** `lamMax W₀ ≤ lamMax (W₀ + q qᵀ)`. -/
theorem lamMax_le_lamMax_vecMulVec (hW₀ : W₀.IsHermitian)
    (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    lamMax W₀ hW₀ ≤ lamMax (W₀ + Matrix.vecMulVec q q) hA := by
  refine lamMax_le_lamMax_of_nonneg hW₀ _ (fun x => ?_) hA
  rw [vecMulVec_mulVec, dotProduct_smul, smul_eq_mul, dotProduct_comm]
  exact mul_self_nonneg _

section OpNorm
open scoped Matrix.Norms.L2Operator

/-- **R4e.2.** `‖θ u vᵀ + E‖ ≤ ‖E‖ + θ` in the `l2` operator norm. -/
theorem opNorm_spiked_le {n : ℕ} (θ : ℝ) (hθ : 0 ≤ θ)
    (u : EuclideanSpace ℝ (Fin n)) (v : EuclideanSpace ℝ (Fin d))
    (hu : ‖u‖ = 1) (hv : ‖v‖ = 1) (E : Matrix (Fin n) (Fin d) ℝ) :
    ‖θ • Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v) + E‖ ≤ ‖E‖ + θ := by
  have hrank : ‖Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v)‖ ≤ 1 := by
    rw [Matrix.l2_opNorm_def]
    refine ContinuousLinearMap.opNorm_le_bound _ zero_le_one fun x => ?_
    have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap)
        (Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v))) x
        = WithLp.toLp 2 (Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v) *ᵥ WithLp.ofLp x) :=
      rfl
    rw [happ, vecMulVec_mulVec]
    have hsm : (WithLp.toLp 2 ((WithLp.ofLp v ⬝ᵥ WithLp.ofLp x) • WithLp.ofLp u)
        : EuclideanSpace ℝ (Fin n)) = (WithLp.ofLp v ⬝ᵥ WithLp.ofLp x) • u := by
      simp
    rw [hsm, norm_smul, hu, mul_one, Real.norm_eq_abs, one_mul]
    have hinner : ⟪v, x⟫_ℝ = WithLp.ofLp v ⬝ᵥ WithLp.ofLp x := by
      rw [EuclideanSpace.inner_eq_star_dotProduct]
      simp [dotProduct_comm]
    rw [← hinner]
    calc |⟪v, x⟫_ℝ| ≤ ‖v‖ * ‖x‖ := abs_real_inner_le_norm v x
      _ = ‖x‖ := by rw [hv, one_mul]
  calc ‖θ • Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v) + E‖
      ≤ ‖θ • Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v)‖ + ‖E‖ := norm_add_le _ _
    _ = θ * ‖Matrix.vecMulVec (WithLp.ofLp u) (WithLp.ofLp v)‖ + ‖E‖ := by
        rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg hθ]
    _ ≤ θ * 1 + ‖E‖ := by nlinarith
    _ = ‖E‖ + θ := by ring

end OpNorm

end R4
end StackedSVD

