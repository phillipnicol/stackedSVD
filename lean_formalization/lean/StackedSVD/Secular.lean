/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R4
import StackedSVD.StackSVDWeighted

/-!
# `lem:secular_equation`: the spectrum of `R = ũ₀ ũ₀ᵀ + Σ`

STATUS 2026-08-30: see `notes/archive/agent_reports/polish_frobenius_secular.md`.

`main_paper.tex:1378` states, for `Σ = diag(w_1², …, w_1², …, w_M², …, w_M²)` with block `i`
repeated `n_i` times and `ũ₀` the block vector `(θ_i w_i u_i)_i` with `‖u_i‖ = 1`:

1. at least `n_i - 1` eigenvalues of `R` equal `w_i²`;
2. the remaining eigenvalues are the roots of `f(λ) = 1 + ∑_j θ_j² w_j²/(w_j² - λ)`;
3. under `eq:assumption4` the largest eigenvalue `γ₁` is unique, with eigenvector
   `ξ₁ ∝ (Σ - γ₁ I)⁻¹ ũ₀`.

`StackSVDWeighted.lean` proved the scalar half of this: `Scalars.secular` is `f`,
`Scalars.IsGammaTop` says "root above `max_i w_i²`", `existsUnique_gammaTop` gives the root and
`sum_secular_eq_neg_one` is the form the limit `L(w)` consumes. Nothing there mentioned the
matrix `R`. This file adds the matrix half, so that the name `γ₁` is tied to the eigenvalue it
denotes in the paper.

## Content

1. Section 1, a general diagonal plus a rank one: `Rmat σ q = diag(σ) + q qᵀ` and
   `secularDiag σ q λ = 1 + ∑_i q_i²/(σ_i - λ)`. `det_Rmat_sub` is the paper's own display
   `det(ũ₀ũ₀ᵀ + Σ - λI) = (1 + ũ₀ᵀ(Σ - λI)⁻¹ũ₀) det(Σ - λI)`, from `R4.det_sub_smul_one`;
   `mem_spectrum_iff_secularDiag` is claim 2, "for `λ ∉ {σ_i}`, `λ` is an eigenvalue of `R` if
   and only if `f(λ) = 0`"; `mulVec_eq_of_const_on` is claim 1 in the geometric form the paper
   proves (`span(e_{block i}) ∩ ũ₀^⊥` sits inside the eigenspace at `w_i²`), and
   `card_sub_one_le_finrank_eigenspace` is the dimension count the paper reads off it, "at
   least `n_i - 1` eigenvalues equal `w_i²`".
2. Section 2, above `max_i σ_i`: `lamMax_diagonal`, `root_above_unique`, `lamMax_Rmat_eq`,
   `topSimple_Rmat` and `topSpace_Rmat_eq_span`, which is claim 3, `ξ₁ ∝ (Σ - γ₁ I)⁻¹ ũ₀`.
   These are `R4.eigenvalue_above_unique`, `R4.lamMax_eq`, `R4.topSimple` and
   `R4.topSpace_eq_span` at `W₀ = diag(σ)`.
3. Section 3, the paper's block matrix: `sigmaStack`, `u0Stack`, `secularDiag_stack` (the block
   sum collapses to `Scalars.secular` because `‖u_i‖ = 1`), `iSup_sigmaStack` (the largest
   diagonal entry is `Scalars.wSqMax w`) and `isGammaTop_iff_stack`. The conclusion is
   `lamMax_Rstack_eq_gammaTop`: under `eq:assumption4` the largest eigenvalue of `R` is
   `Scalars.gammaTop θ w` and it is simple.

## Hypotheses

`det_Rmat_sub` and `mem_spectrum_iff_secularDiag` need only `λ ∉ {σ_i}`; positivity of the
`σ_i` is not used, and `q = 0` is allowed. The results of section 2 need `q ≠ 0`, which is the
paper's implicit "some table carries signal": `Scalars.exists_signal_of_root` derives it from
the existence of the root, and `u0Stack_ne_zero` turns it into `ũ₀ ≠ 0`. Section 3 needs every
block nonempty (`0 < n_i`), since an empty block contributes no diagonal entry `w_i²` and
`Scalars.wSqMax` would then be too large.

## What this file does not prove

The tie to a model is not formalized. No declaration here says that
`Rmat (sigmaStack ν w) (u0Stack θ w u)` is `E[X_stack X_stackᵀ]` of a `MultiTableModel`, and
nothing connects the paper's `η₁ (ũ₀ᵀ ξ₁)²/γ₁` (`main_paper.tex:1414`) to `Scalars.Lw`. The
file is a self-contained spectral statement about a diagonal matrix plus a rank one, with
`σ` and `q` supplied by the caller; `sigmaStack` and `u0Stack` only put them in the paper's
block shape. Nothing outside this file reads it (audit of the unaudited pieces, 2026-08-31,
finding 2). Closing the gap needs the expectation of the weighted stack Gram matrix, which
`StackSVDWeighted.lean` does not define today.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace Secular

/-! ### 1. A diagonal matrix plus a rank one -/

section Diag

variable {nn : ℕ}

/-- `R = ũ₀ ũ₀ᵀ + Σ` of `lem:secular_equation`, with `Σ = diag(σ)` and `ũ₀ = q`. The summands
are written in the order `R4` uses; `Rmat_eq` gives the paper's order. -/
noncomputable def Rmat (σ q : Fin nn → ℝ) : Matrix (Fin nn) (Fin nn) ℝ :=
  Matrix.diagonal σ + Matrix.vecMulVec q q

/-- The paper's order, `ũ₀ ũ₀ᵀ + Σ`. -/
theorem Rmat_eq (σ q : Fin nn → ℝ) :
    Rmat σ q = Matrix.vecMulVec q q + Matrix.diagonal σ := add_comm _ _

theorem isHermitian_Rmat (σ q : Fin nn → ℝ) : (Rmat σ q).IsHermitian := by
  refine Matrix.IsHermitian.add (Matrix.isHermitian_diagonal _) ?_
  ext i j
  simp [Matrix.vecMulVec_apply, Matrix.conjTranspose_apply, mul_comm]

/-- The secular function of `lem:secular_equation`, written on a general diagonal:
`f(λ) = 1 + ∑_i q_i²/(σ_i - λ)`. -/
noncomputable def secularDiag (σ q : Fin nn → ℝ) (lam : ℝ) : ℝ :=
  1 + ∑ i, q i ^ 2 / (σ i - lam)

/-- `Σ - λ I = diag(σ - λ)`. -/
theorem diagonal_sub_smul_one (σ : Fin nn → ℝ) (z : ℝ) :
    Matrix.diagonal σ - z • (1 : Matrix (Fin nn) (Fin nn) ℝ)
      = Matrix.diagonal fun i => σ i - z := by
  rw [Matrix.smul_one_eq_diagonal, ← Matrix.diagonal_sub]

/-- `det(Σ - λ I) = ∏_i (σ_i - λ)`. -/
theorem det_diagonal_sub (σ : Fin nn → ℝ) (z : ℝ) :
    (Matrix.diagonal σ - z • (1 : Matrix (Fin nn) (Fin nn) ℝ)).det = ∏ i, (σ i - z) := by
  rw [diagonal_sub_smul_one, Matrix.det_diagonal]

/-- `(Σ - λ I)⁻¹ = diag((σ - λ)⁻¹)` away from the diagonal entries. -/
theorem resolv_diagonal {σ : Fin nn → ℝ} {z : ℝ} (hz : ∀ i, σ i ≠ z) :
    R4.resolv (Matrix.diagonal σ) z = Matrix.diagonal fun i => (σ i - z)⁻¹ := by
  rw [R4.resolv, diagonal_sub_smul_one]
  refine Matrix.inv_eq_right_inv ?_
  rw [Matrix.diagonal_mul_diagonal]
  have hfun : (fun i => (σ i - z) * (σ i - z)⁻¹) = fun _ : Fin nn => (1 : ℝ) := by
    funext i
    exact mul_inv_cancel₀ (sub_ne_zero.mpr (hz i))
  rw [hfun]
  simp

/-- `R4.secular` at a diagonal matrix is `secularDiag`. -/
theorem R4secular_diagonal {σ q : Fin nn → ℝ} {z : ℝ} (hz : ∀ i, σ i ≠ z) :
    R4.secular (Matrix.diagonal σ) q z = secularDiag σ q z := by
  rw [R4.secular, R4.qform, resolv_diagonal hz, secularDiag]
  congr 1
  rw [dotProduct]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.mulVec_diagonal]
  rw [div_eq_inv_mul, sq]
  ring

/-- **`lem:secular_equation`, the determinant identity** (`main_paper.tex:1392`):
`det(ũ₀ũ₀ᵀ + Σ - λI) = det(Σ - λI) (1 + ũ₀ᵀ(Σ - λI)⁻¹ũ₀)`. -/
theorem det_Rmat_sub (σ q : Fin nn → ℝ) {lam : ℝ} (hlam : ∀ i, σ i ≠ lam) :
    (Rmat σ q - lam • (1 : Matrix (Fin nn) (Fin nn) ℝ)).det
      = (∏ i, (σ i - lam)) * secularDiag σ q lam := by
  have hunit : IsUnit (Matrix.diagonal σ - lam • (1 : Matrix (Fin nn) (Fin nn) ℝ)).det := by
    rw [det_diagonal_sub]
    exact isUnit_iff_ne_zero.mpr
      (Finset.prod_ne_zero_iff.mpr fun i _ => sub_ne_zero.mpr (hlam i))
  rw [Rmat, R4.vecMulVec_eq_replicateCol,
    R4.det_sub_smul_one (Q := Matrix.replicateCol (Fin 1) q) hunit, R4.det_one_add_col,
    det_diagonal_sub]
  congr 1
  rw [← R4secular_diagonal (q := q) hlam]
  rfl

/-- **`lem:secular_equation`, claim 2.** A number that is not a diagonal entry of `Σ` is an
eigenvalue of `R` if and only if it is a root of the secular equation. -/
theorem mem_spectrum_iff_secularDiag (σ q : Fin nn → ℝ) {lam : ℝ} (hlam : ∀ i, σ i ≠ lam) :
    lam ∈ spectrum ℝ (Rmat σ q) ↔ secularDiag σ q lam = 0 := by
  have hprod : (∏ i, (σ i - lam)) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr fun i _ => sub_ne_zero.mpr (hlam i)
  have hscalar : (Matrix.scalar (Fin nn)) lam = lam • (1 : Matrix (Fin nn) (Fin nn) ℝ) := by
    rw [Matrix.smul_one_eq_diagonal]
    rfl
  have hneg : (lam • (1 : Matrix (Fin nn) (Fin nn) ℝ) - Rmat σ q)
      = -(Rmat σ q - lam • 1) := by abel
  rw [Matrix.mem_spectrum_iff_not_isUnit_eval_charpoly, Matrix.eval_charpoly,
    isUnit_iff_ne_zero, not_ne_iff, hscalar, hneg, Matrix.det_neg, det_Rmat_sub σ q hlam]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h' | h'
    · exact absurd h' (pow_ne_zero _ (by norm_num))
    · rcases mul_eq_zero.mp h' with h'' | h''
      · exact absurd h'' hprod
      · exact h''
  · intro h
    rw [h, mul_zero, mul_zero]

/-- **`lem:secular_equation`, claim 1**, in the geometric form the paper proves: a vector
supported where `Σ` is constant and orthogonal to `ũ₀` is an eigenvector of `R` at that
constant. The paper reads "at least `n_i - 1` eigenvalues equal `w_i²`" off this, because
`span(e_{block i}) ∩ ũ₀^⊥` has dimension at least `n_i - 1`. -/
theorem mulVec_eq_of_const_on (σ q : Fin nn → ℝ) {t : ℝ} {S : Finset (Fin nn)}
    (hS : ∀ i ∈ S, σ i = t) {x : Fin nn → ℝ} (hsupp : ∀ i, i ∉ S → x i = 0)
    (hperp : q ⬝ᵥ x = 0) : Rmat σ q *ᵥ x = t • x := by
  have hzero : Matrix.vecMulVec q q *ᵥ x = 0 := by
    rw [R4.vecMulVec_mulVec, hperp, zero_smul]
  funext i
  rw [Rmat, Matrix.add_mulVec, hzero]
  simp only [Pi.add_apply, Pi.zero_apply, add_zero, Matrix.mulVec_diagonal, Pi.smul_apply,
    smul_eq_mul]
  by_cases hi : i ∈ S
  · rw [hS i hi]
  · rw [hsupp i hi, mul_zero, mul_zero]

/-- **`lem:secular_equation`, claim 1, dimension form.** If `Σ` is constant equal to `t` on a
set `S` of indices, the eigenspace of `R` at `t` has dimension at least `|S| - 1`. At the
paper's block `Σ`, with `S` the `i`-th block, this is "at least `n_i - 1` eigenvalues of `R`
equal `w_i²`". The proof is the paper's: `span(e_S)` has dimension `|S|`, one linear condition
(`⟪ũ₀, x⟫ = 0`) cuts it by at most one, and `mulVec_eq_of_const_on` puts what is left inside
the eigenspace. -/
theorem card_sub_one_le_finrank_eigenspace (σ q : Fin nn → ℝ) {t : ℝ} {S : Finset (Fin nn)}
    (hS : ∀ i ∈ S, σ i = t) :
    S.card - 1 ≤ Module.finrank ℝ (Module.End.eigenspace (toOp (Rmat σ q)) t) := by
  classical
  have hb : ∀ i : Fin nn,
      (EuclideanSpace.basisFun (Fin nn) ℝ).toBasis i = EuclideanSpace.single i (1 : ℝ) := by
    intro i
    rw [OrthonormalBasis.coe_toBasis]
    simp
  have hli : LinearIndependent ℝ
      (fun i : S => (EuclideanSpace.single (i : Fin nn) (1 : ℝ) : EuclideanSpace ℝ (Fin nn))) := by
    have h0 := (EuclideanSpace.basisFun (Fin nn) ℝ).toBasis.linearIndependent
    have hinj : Function.Injective (fun i : S => (i : Fin nn)) := fun a b h => Subtype.ext h
    have heq : ((EuclideanSpace.basisFun (Fin nn) ℝ).toBasis ∘ fun i : S => (i : Fin nn))
        = fun i : S => (EuclideanSpace.single (i : Fin nn) (1 : ℝ) : EuclideanSpace ℝ (Fin nn)) :=
      funext fun i => hb _
    rw [← heq]
    exact h0.comp _ hinj
  set W : Submodule ℝ (EuclideanSpace ℝ (Fin nn)) :=
    Submodule.span ℝ (Set.range fun i : S =>
      (EuclideanSpace.single (i : Fin nn) (1 : ℝ) : EuclideanSpace ℝ (Fin nn))) with hWdef
  set L : EuclideanSpace ℝ (Fin nn) →ₗ[ℝ] ℝ :=
    (innerSL ℝ (WithLp.toLp 2 q : EuclideanSpace ℝ (Fin nn))).toLinearMap with hLdef
  have hWcard : Module.finrank ℝ W = S.card := by
    rw [hWdef, finrank_span_eq_card hli, Fintype.card_coe]
  have hsupp : ∀ x ∈ W, ∀ i, i ∉ S → x i = 0 := by
    intro x hx i hi
    rw [hWdef] at hx
    refine Submodule.span_induction ?_ ?_ ?_ ?_ hx
    · rintro _ ⟨k, rfl⟩
      have hne : i ≠ (k : Fin nn) := fun h => hi (h ▸ k.2)
      simp [hne]
    · simp
    · intro u v _ _ hu hv
      simp [hu, hv]
    · intro a u _ hu
      simp [hu]
  have hE : Module.finrank ℝ (EuclideanSpace ℝ (Fin nn)) = nn := finrank_euclideanSpace_fin
  have hker : nn ≤ Module.finrank ℝ (LinearMap.ker L) + 1 := by
    have h := LinearMap.finrank_range_add_finrank_ker L
    rw [hE] at h
    have hr : Module.finrank ℝ (LinearMap.range L) ≤ 1 := by
      have h1 := Submodule.finrank_le (LinearMap.range L)
      simpa using h1
    omega
  have hinf := Submodule.finrank_sup_add_finrank_inf_eq W (LinearMap.ker L)
  have hsuple : Module.finrank ℝ ↥(W ⊔ LinearMap.ker L) ≤ nn := by
    have h := Submodule.finrank_le (W ⊔ LinearMap.ker L)
    rwa [hE] at h
  have hfinal : S.card - 1 ≤ Module.finrank ℝ ↥(W ⊓ LinearMap.ker L) := by
    rw [hWcard] at hinf
    omega
  have hle : W ⊓ LinearMap.ker L ≤ Module.End.eigenspace (toOp (Rmat σ q)) t := by
    rintro x ⟨hxW, hxK⟩
    rw [R4.mem_eigenspace_iff']
    refine mulVec_eq_of_const_on σ q hS (fun i hi => hsupp x hxW i hi) ?_
    have hz : ⟪(WithLp.toLp 2 q : EuclideanSpace ℝ (Fin nn)), x⟫_ℝ = 0 :=
      LinearMap.mem_ker.mp hxK
    rwa [real_inner_eq_dotProduct] at hz
  exact le_trans hfinal (Submodule.finrank_mono hle)

end Diag

/-! ### 2. The largest eigenvalue -/

section Top

variable {nn : ℕ}

/-- The largest eigenvalue of a diagonal matrix is the largest diagonal entry. -/
theorem lamMax_diagonal (σ : Fin nn → ℝ) (hnn : 0 < nn) :
    lamMax (Matrix.diagonal σ) (Matrix.isHermitian_diagonal σ) = ⨆ i, σ i := by
  have : Nonempty (Fin nn) := Fin.pos_iff_nonempty.mp hnn
  set hD := Matrix.isHermitian_diagonal σ with hDdef
  have hspec : spectrum ℝ (Matrix.diagonal σ) = Set.range σ := _root_.spectrum_diagonal σ
  have hrange : Set.range hD.eigenvalues = Set.range σ := by
    rw [← hD.spectrum_real_eq_range_eigenvalues, hspec]
  have hle : ∀ i, σ i ≤ lamMax (Matrix.diagonal σ) hD := by
    intro i
    have hmem : σ i ∈ Set.range hD.eigenvalues := by rw [hrange]; exact ⟨i, rfl⟩
    obtain ⟨k, hk⟩ := hmem
    rw [← hk]
    exact R4.eigenvalues_le_lamMax hD k
  have hge : lamMax (Matrix.diagonal σ) hD ≤ ⨆ i, σ i := by
    obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hD hnn
    have hmem : hD.eigenvalues j ∈ Set.range σ := by rw [← hrange]; exact ⟨j, rfl⟩
    obtain ⟨i, hi⟩ := hmem
    rw [← hj, ← hi]
    exact le_ciSup (Set.Finite.bddAbove (Set.finite_range _)) i
  exact le_antisymm hge (ciSup_le hle)

variable {σ q : Fin nn → ℝ}

/-- Above the largest diagonal entry, the diagonal entries are not met. -/
private theorem ne_of_lt_iSup {z : ℝ} (hz : (⨆ i, σ i) < z) (hnn : 0 < nn) : ∀ i, σ i ≠ z := by
  have : Nonempty (Fin nn) := Fin.pos_iff_nonempty.mp hnn
  intro i hi
  exact absurd (hi ▸ le_ciSup (Set.Finite.bddAbove (Set.finite_range _)) i) (not_le.mpr hz)

/-- The hypothesis `R4` consumes, in terms of the diagonal entries. -/
private theorem lamMax_lt {z : ℝ} (hz : (⨆ i, σ i) < z) (hnn : 0 < nn) :
    lamMax (Matrix.diagonal σ) (Matrix.isHermitian_diagonal σ) < z := by
  rwa [lamMax_diagonal σ hnn]

/-- A root of the secular equation forces `q ≠ 0`, since `f` is the constant `1` at `q = 0`. -/
theorem ne_zero_of_secularDiag_eq_zero {lam : ℝ} (e : secularDiag σ q lam = 0) : q ≠ 0 := by
  rintro rfl
  rw [secularDiag] at e
  simp at e

/-- A root of the secular equation forces `0 < nn`, since `q ≠ 0` needs a nonempty index. -/
theorem pos_of_secularDiag_eq_zero {lam : ℝ} (e : secularDiag σ q lam = 0) : 0 < nn := by
  rcases Nat.eq_zero_or_pos nn with h | h
  · subst h
    exact absurd (by ext i; exact i.elim0) (ne_zero_of_secularDiag_eq_zero e)
  · exact h

/-- **`lem:secular_equation`, claim 3, uniqueness.** At most one root of the secular equation
lies above `max_i σ_i`. This is the Bunch-Nielsen-Sorensen interlacing statement the paper
cites; here it is the strict monotonicity of `f` on `(max_i σ_i, ∞)`. -/
theorem root_above_unique {z₁ z₂ : ℝ}
    (h₁ : (⨆ i, σ i) < z₁) (h₂ : (⨆ i, σ i) < z₂)
    (e₁ : secularDiag σ q z₁ = 0) (e₂ : secularDiag σ q z₂ = 0) : z₁ = z₂ := by
  have hq : q ≠ 0 := ne_zero_of_secularDiag_eq_zero e₁
  have hnn : 0 < nn := pos_of_secularDiag_eq_zero e₁
  refine R4.eigenvalue_above_unique (Matrix.isHermitian_diagonal σ) hq
    (lamMax_lt h₁ hnn) (lamMax_lt h₂ hnn) ?_ ?_
  · rw [R4secular_diagonal (ne_of_lt_iSup h₁ hnn)]; exact e₁
  · rw [R4secular_diagonal (ne_of_lt_iSup h₂ hnn)]; exact e₂

/-- **`lem:secular_equation`, claim 3.** A root of the secular equation above `max_i σ_i` is
the largest eigenvalue `γ₁` of `R`. -/
theorem lamMax_Rmat_eq {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) :
    lamMax (Rmat σ q) (isHermitian_Rmat σ q) = lam := by
  have hq : q ≠ 0 := ne_zero_of_secularDiag_eq_zero e
  have hnn : 0 < nn := pos_of_secularDiag_eq_zero e
  refine R4.lamMax_eq (Matrix.isHermitian_diagonal σ) hq (lamMax_lt hlam hnn) ?_
    (isHermitian_Rmat σ q)
  rw [R4secular_diagonal (ne_of_lt_iSup hlam hnn)]
  exact e

/-- **`lem:secular_equation`, claim 3, simplicity.** `γ₁` is a simple eigenvalue of `R`. -/
theorem topSimple_Rmat {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) : TopSimple (Rmat σ q) (isHermitian_Rmat σ q) := by
  have hq : q ≠ 0 := ne_zero_of_secularDiag_eq_zero e
  have hnn : 0 < nn := pos_of_secularDiag_eq_zero e
  refine R4.topSimple (Matrix.isHermitian_diagonal σ) hq (lamMax_lt hlam hnn) ?_
    (isHermitian_Rmat σ q)
  rw [R4secular_diagonal (ne_of_lt_iSup hlam hnn)]
  exact e

/-- **`lem:secular_equation`, claim 3, the eigenvector** `ξ₁ ∝ (Σ - γ₁ I)⁻¹ ũ₀`: the top
eigenspace of `R` is the line through the vector with entries `q_i/(σ_i - γ₁)`. -/
theorem topSpace_Rmat_eq_span {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) :
    topSpace (Rmat σ q) (isHermitian_Rmat σ q)
      = Submodule.span ℝ
          {(WithLp.toLp 2 fun i => q i / (σ i - lam) : EuclideanSpace ℝ (Fin nn))} := by
  have hq : q ≠ 0 := ne_zero_of_secularDiag_eq_zero e
  have hnn : 0 < nn := pos_of_secularDiag_eq_zero e
  have hne := ne_of_lt_iSup hlam hnn
  have hres : R4.resolv (Matrix.diagonal σ) lam *ᵥ q = fun i => q i / (σ i - lam) := by
    rw [resolv_diagonal hne]
    funext i
    rw [Matrix.mulVec_diagonal, div_eq_inv_mul]
  rw [← hres]
  refine R4.topSpace_eq_span (Matrix.isHermitian_diagonal σ) hq (lamMax_lt hlam hnn) ?_
    (isHermitian_Rmat σ q)
  rw [R4secular_diagonal hne]
  exact e

end Top

/-! ### 3. The paper's block matrix and `Scalars.gammaTop` -/

section Stack

variable {M : ℕ} {ν : Fin M → ℕ}

/-- A sum over the stacked index is the double sum over the block index. `RankR/Subspace.lean`
has the same lemma, private. -/
private theorem sum_stack_indexS {α : Type*} [AddCommMonoid α]
    (F : (i : Fin M) → Fin (ν i) → α) :
    ∑ p : Fin (∑ i, ν i), F (finSigmaFinEquiv.symm p).1 (finSigmaFinEquiv.symm p).2
      = ∑ i, ∑ j, F i j := by
  rw [Equiv.sum_comp finSigmaFinEquiv.symm (fun p : (i : Fin M) × Fin (ν i) => F p.1 p.2)]
  exact Fintype.sum_sigma' F

/-- `Σ` of `lem:secular_equation`: the diagonal with `w_i²` repeated `n_i` times. -/
noncomputable def sigmaStack (ν : Fin M → ℕ) (w : Fin M → ℝ) : Fin (∑ i, ν i) → ℝ :=
  fun p => w (finSigmaFinEquiv.symm p).1 ^ 2

/-- `ũ₀` of `lem:secular_equation`: block `i` is `θ_i w_i u_i`. -/
noncomputable def u0Stack (θ w : Fin M → ℝ) (u : (i : Fin M) → Fin (ν i) → ℝ) :
    Fin (∑ i, ν i) → ℝ :=
  fun p => θ (finSigmaFinEquiv.symm p).1 * w (finSigmaFinEquiv.symm p).1 *
    u (finSigmaFinEquiv.symm p).1 (finSigmaFinEquiv.symm p).2

/-- The block sum collapses: `∑_{k} (θ_i w_i u_{ik})²/(w_i² - λ) = θ_i² w_i²/(w_i² - λ)`,
because `‖u_i‖ = 1`. So the secular function of the stacked `(Σ, ũ₀)` is the scalar
`Scalars.secular θ w` of `StackSVDWeighted.lean`. -/
theorem secularDiag_stack (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (lam : ℝ) :
    secularDiag (sigmaStack ν w) (u0Stack θ w u) lam = Scalars.secular θ w lam := by
  rw [secularDiag, Scalars.secular]
  congr 1
  have hkey := sum_stack_indexS (ν := ν)
    (fun (i : Fin M) (k : Fin (ν i)) => (θ i * w i * u i k) ^ 2 / (w i ^ 2 - lam))
  calc ∑ p : Fin (∑ i, ν i), u0Stack θ w u p ^ 2 / (sigmaStack ν w p - lam)
      = ∑ i, ∑ k, (θ i * w i * u i k) ^ 2 / (w i ^ 2 - lam) := hkey
    _ = ∑ j, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - lam) := by
        refine Finset.sum_congr rfl fun i _ => ?_
        have hnum : ∑ k, (θ i * w i * u i k) ^ 2 = θ i ^ 2 * w i ^ 2 := by
          have hfac : (∑ k, (θ i * w i * u i k) ^ 2)
              = (θ i ^ 2 * w i ^ 2) * ∑ k, u i k ^ 2 := by
            rw [Finset.mul_sum]
            exact Finset.sum_congr rfl fun k _ => by ring
          rw [hfac, hu i, mul_one]
        rw [← Finset.sum_div, hnum]

/-- The largest diagonal entry of the stacked `Σ` is `max_i w_i²`, provided every block is
nonempty. -/
theorem iSup_sigmaStack (w : Fin M → ℝ) (hν : ∀ i, 0 < ν i) (hM : 0 < M) :
    (⨆ p, sigmaStack ν w p) = Scalars.wSqMax w := by
  have : Nonempty (Fin M) := Fin.pos_iff_nonempty.mp hM
  have : Nonempty (Fin (∑ i, ν i)) := by
    refine Fin.pos_iff_nonempty.mp ?_
    exact Finset.sum_pos (fun i _ => hν i) Finset.univ_nonempty
  refine le_antisymm (ciSup_le fun p => Scalars.le_wSqMax w _) (ciSup_le fun i => ?_)
  have hk : Nonempty (Fin (ν i)) := Fin.pos_iff_nonempty.mp (hν i)
  obtain ⟨k⟩ := hk
  have hval : sigmaStack ν w (finSigmaFinEquiv ⟨i, k⟩) = w i ^ 2 := by
    simp [sigmaStack]
  rw [← hval]
  exact le_ciSup (Set.Finite.bddAbove (Set.finite_range _)) _

/-- `Scalars.IsGammaTop` is exactly "root of the stacked secular equation above the largest
diagonal entry of `Σ`". -/
theorem isGammaTop_iff_stack (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (hν : ∀ i, 0 < ν i) (hM : 0 < M) (g : ℝ) :
    Scalars.IsGammaTop θ w g ↔
      ((⨆ p, sigmaStack ν w p) < g ∧ secularDiag (sigmaStack ν w) (u0Stack θ w u) g = 0) := by
  rw [Scalars.IsGammaTop, iSup_sigmaStack w hν hM, secularDiag_stack θ w hu]

/-- A table with nonzero signal and nonzero weight makes `ũ₀ ≠ 0`. -/
theorem u0Stack_ne_zero (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) {j : Fin M} (hj : θ j * w j ≠ 0) :
    u0Stack θ w u ≠ (0 : Fin (∑ i, ν i) → ℝ) := by
  obtain ⟨k, hk⟩ : ∃ k, u j k ^ 2 ≠ 0 := by
    by_contra hcon
    have hz : ∀ k, u j k ^ 2 = 0 := fun k => not_not.mp (not_exists.mp hcon k)
    have h1 := hu j
    rw [Finset.sum_congr rfl fun k (_ : k ∈ Finset.univ) => hz k] at h1
    simp at h1
  have hP : finSigmaFinEquiv.symm (finSigmaFinEquiv (⟨j, k⟩ : (i : Fin M) × Fin (ν i)))
      = ⟨j, k⟩ := Equiv.symm_apply_apply _ _
  have hval : u0Stack θ w u (finSigmaFinEquiv ⟨j, k⟩) = θ j * w j * u j k :=
    congrArg (fun X : (i : Fin M) × Fin (ν i) => θ X.1 * w X.1 * u X.1 X.2) hP
  intro h0
  have h1 := congrFun h0 (finSigmaFinEquiv ⟨j, k⟩)
  rw [hval] at h1
  simp only [Pi.zero_apply] at h1
  rcases mul_eq_zero.mp h1 with h | h
  · exact hj h
  · exact hk (by rw [h]; ring)

/-- **`lem:secular_equation`, the conclusion.** Under `eq:assumption4` (of which the part used
here is the existence of the outlier root) the largest eigenvalue of `R = ũ₀ũ₀ᵀ + Σ` is
`Scalars.gammaTop θ w`, and it is simple. This is what ties the scalar `γ₁` of
`StackSVDWeighted.lean` to the eigenvalue the paper denotes by `γ₁`. -/
theorem lamMax_Rstack_eq_gammaTop (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (hν : ∀ i, 0 < ν i) (hM : 0 < M)
    (hg : ∃ g, Scalars.IsGammaTop θ w g) :
    lamMax (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u)) = Scalars.gammaTop θ w
      ∧ TopSimple (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u)) := by
  obtain ⟨hlt, he⟩ :=
    (isGammaTop_iff_stack θ w hu hν hM _).mp (Scalars.gammaTop_of_exists hg)
  exact ⟨lamMax_Rmat_eq hlt he, topSimple_Rmat hlt he⟩

/-- The eigenvector of `γ₁`, in block form: `ξ₁ ∝ (Σ - γ₁ I)⁻¹ ũ₀`, whose block `i` is
`θ_i w_i u_i/(w_i² - γ₁)`. This is the paper's second display after `lem:secular_equation`. -/
theorem topSpace_Rstack_eq_span (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (hν : ∀ i, 0 < ν i) (hM : 0 < M)
    (hg : ∃ g, Scalars.IsGammaTop θ w g) :
    topSpace (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u))
      = Submodule.span ℝ
          {(WithLp.toLp 2 fun p => u0Stack θ w u p /
              (sigmaStack ν w p - Scalars.gammaTop θ w) :
            EuclideanSpace ℝ (Fin (∑ i, ν i)))} := by
  obtain ⟨hlt, he⟩ :=
    (isGammaTop_iff_stack θ w hu hν hM _).mp (Scalars.gammaTop_of_exists hg)
  exact topSpace_Rmat_eq_span hlt he

end Stack

end Secular

end StackedSVD
