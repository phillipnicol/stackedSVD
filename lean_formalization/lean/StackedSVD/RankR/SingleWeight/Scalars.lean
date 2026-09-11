/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVDWeighted

/-!
# `prop:gen_rank_stacksvd_singleweight`: the deterministic scalar layer

Units S1 and S2 of `notes/archive/singleweight_plan.md` section 4.1. Namespace
`StackedSVD.SingleWeight`. Nothing here is random: every declaration is an identity or an
existence claim about finitely many real matrices.

## The paper

`main_paper.tex:2112`, `prop:gen_rank_stacksvd_singleweight`. Stacksvd with one weight `w_i`
per table has performance controlled by the spectrum of `G = A Aᵀ + Σ`, whose outliers
`γ_1 > … > γ_r` solve the matrix secular equation

```
det( I_r - ∑_i w_i²/(γ - w_i²) R_i Θ_i² R_iᵀ ) = 0,      γ > sup_i w_i²    (main_paper.tex:2115)
```

and whose eigenvector data enters the limit

```
‖V̂_stacksvd(w)ᵀ V‖_F² →p ∑_ℓ  (1 - ∑_i c_i w_i⁴/(γ_ℓ - w_i²)²)
                             / (γ_ℓ z_ℓᵀ [∑_i w_i²/(w_i² - γ_ℓ)² R_i Θ_i² R_iᵀ] z_ℓ)
                                                                            (main_paper.tex:2119)
```

## Content

1. `sigMat`, `secMat`, `secDerivMat`: the three matrices of the two displays.
2. `IsSecularRoot`: `γ` is a root above `sup_i w_i²`.
3. `exists_unit_eigvec_secMat`: the first half of the proposition, the unit eigenvector `z_ℓ`.
4. `swTerm`, `swLimit`: one summand and the whole right side of `main_paper.tex:2119`.
5. `EigSep`: `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2104`) on the data `(γ, z)`.
6. `isHermitian_secMat`, `posSemidef_secMat`, `posSemidef_secDerivMat`: the three structural
   facts that the downstream proofs read.

STATUS 2026-09-05: proved (Layer 1, user OK on the statements the same day). The Gaussian
discharge of `SingleWeightLaw` is Track G (`notes/archive/trackG_plan.md`).
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

/-! ### 1. The three matrices of the proposition -/

/-- `S_i = R_i Θ_i² R_iᵀ`, the signal covariance of table `i` in the shared basis. It is the
matrix that carries `Θ_i` into the displays of `main_paper.tex:2115` and `:2119`. -/
noncomputable def sigMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) :
    Matrix (Fin r) (Fin r) ℝ :=
  R i * Matrix.diagonal (fun j => θ i j ^ 2) * (R i)ᵀ

/-- `∑_i w_i²/(γ - w_i²) R_i Θ_i² R_iᵀ`, the matrix of `main_paper.tex:2115`. -/
noncomputable def secMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  ∑ i, (w i ^ 2 / (γ - w i ^ 2)) • sigMat θ R i

/-- `∑_i w_i²/(γ - w_i²)² R_i Θ_i² R_iᵀ`, the matrix inside the denominator of
`main_paper.tex:2119`. It is `- d/dγ (secMat θ R w γ)`. The paper writes the denominator
`(w_i² - γ)²`, which is the same number. -/
noncomputable def secDerivMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  ∑ i, (w i ^ 2 / (γ - w i ^ 2) ^ 2) • sigMat θ R i

/-- `γ` is a root of the rank-`r` secular equation above `max_i w_i²`
(`main_paper.tex:2115`). The rank-one twin is `Scalars.IsGammaTop`
(`StackSVDWeighted.lean:107`). -/
def IsSecularRoot {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) : Prop :=
  Scalars.wSqMax w < γ ∧ (1 - secMat θ R w γ).det = 0

/-! ### 2. The first half of the proposition -/

/-- The first half of **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`): a
root of the secular equation gives the matrix of `main_paper.tex:2115` a unit eigenvector
`z_ℓ` with eigenvalue `1`. Route: `1 - secMat` is singular, so its kernel is nonzero
(`Matrix.exists_mulVec_eq_zero_iff`); normalize any kernel vector. -/
theorem exists_unit_eigvec_secMat {M r : ℕ} {rk : Fin M → ℕ}
    {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w : Fin M → ℝ} {γ : ℝ}
    (h : IsSecularRoot θ R w γ) :
    ∃ z : EuclideanSpace ℝ (Fin r), ‖z‖ = 1 ∧
      secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z := by
  obtain ⟨v, hv0, hveq⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr h.2
  rw [Matrix.sub_mulVec, Matrix.one_mulVec, sub_eq_zero] at hveq
  have hkey : secMat θ R w γ *ᵥ v = v := hveq.symm
  have he0 : (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin r)) ≠ 0 := by simpa using hv0
  have hnorm0 : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin r))‖ ≠ 0 := norm_ne_zero_iff.mpr he0
  refine ⟨‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin r))‖⁻¹ •
      (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin r)), ?_, ?_⟩
  · rw [norm_smul, norm_inv, Real.norm_eq_abs,
      abs_of_nonneg (norm_nonneg (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin r))),
      inv_mul_cancel₀ hnorm0]
  · rw [WithLp.ofLp_smul, WithLp.ofLp_toLp, Matrix.mulVec_smul, hkey]

/-! ### 3. The limit of the proposition -/

/-- One summand of the limit of `main_paper.tex:2119`. -/
noncomputable def swTerm {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ) (γ : ℝ)
    (z : EuclideanSpace ℝ (Fin r)) : ℝ :=
  (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) /
    (γ * (WithLp.ofLp z ⬝ᵥ (secDerivMat θ R w γ *ᵥ WithLp.ofLp z)))

/-- The right side of `main_paper.tex:2119`, the limit of
`prop:gen_rank_stacksvd_singleweight`. -/
noncomputable def swLimit {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : ℝ :=
  ∑ l, swTerm θ R w c (γ l) (z l)

/-! ### 4. `assum:gen_rank_stacksvd_eig_sep` -/

/-- `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2104`) as a structure on the data
`(γ, z)`. `sorted` replaces the paper's "the top `r` eigenvalues of `G`, which are distinct":
`γ` lists the roots in strictly decreasing order, so `γ l` is the `l`-th one and the sorted
index of the matching outlier of the stack Gram matrix is `l`.

Scope. A tie is excluded by `sorted`, and dropping it is not free: the limit then depends on
which basis of the tied eigenspace the `z_ℓ` form. The tie tolerant replacement is the
`secDerivMat`-orthogonality of the `z_ℓ`, which is vacuous under `sorted`
(`notes/archive/singleweight_plan.md` section 4.1). -/
structure EigSep {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : Prop where
  /-- each `γ_ℓ` is a root of the secular equation above `sup_i w_i²` -/
  root : ∀ l, IsSecularRoot θ R w (γ l)
  /-- the roots are listed in strictly decreasing order, so no two of them are tied -/
  sorted : StrictAnti γ
  /-- `z_ℓ` is a unit eigenvector of the matrix of `main_paper.tex:2115` at eigenvalue `1` -/
  eigvec : ∀ l, ‖z l‖ = 1 ∧ secMat θ R w (γ l) *ᵥ WithLp.ofLp (z l) = WithLp.ofLp (z l)
  /-- the detectability threshold of `main_paper.tex:2106` -/
  thresh : ∀ l, ∑ i, c i * w i ^ 4 / (γ l - w i ^ 2) ^ 2 < 1

/-! ### 5. The three structural facts -/

/-- `sigMat` is positive semidefinite: it is `R_i` times a diagonal matrix of squares times
`R_iᵀ`. Private helper for the three theorems below; its Hermitian half (`.1`) also gives
`isHermitian_secMat`. -/
private theorem posSemidef_sigMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) :
    (sigMat θ R i).PosSemidef := by
  unfold sigMat
  rw [← Matrix.conjTranspose_eq_transpose_of_trivial (R i)]
  exact (Matrix.posSemidef_diagonal_iff.2 (fun j => sq_nonneg (θ i j))).mul_mul_conjTranspose_same
    (R i)

/-- `secMat` is symmetric: every summand is a scalar multiple of `R_i Θ_i² R_iᵀ`. -/
theorem isHermitian_secMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    (secMat θ R w γ).IsHermitian := by
  unfold secMat
  refine Finset.sum_induction _ Matrix.IsHermitian (fun a b ha hb => ha.add hb)
    Matrix.isHermitian_zero (fun i _ => ?_)
  exact (posSemidef_sigMat θ R i).1.smul (IsSelfAdjoint.all (w i ^ 2 / (γ - w i ^ 2)))

/-- Above `max_i w_i²` every coefficient `w_i²/(γ - w_i²)` is nonnegative, so `secMat` is
positive semidefinite. -/
theorem posSemidef_secMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) {γ : ℝ}
    (hγ : Scalars.wSqMax w < γ) :
    (secMat θ R w γ).PosSemidef := by
  unfold secMat
  refine Matrix.posSemidef_sum _ (fun i _ => ?_)
  have hpos : 0 < γ - w i ^ 2 := by linarith [Scalars.le_wSqMax w i]
  exact (posSemidef_sigMat θ R i).smul (div_nonneg (sq_nonneg (w i)) hpos.le)

/-- `secDerivMat` is positive semidefinite for every `γ`: each coefficient
`w_i²/(γ - w_i²)²` is a ratio of two squares. This is what makes the denominator of
`swTerm` nonnegative. -/
theorem posSemidef_secDerivMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    (secDerivMat θ R w γ).PosSemidef := by
  unfold secDerivMat
  refine Matrix.posSemidef_sum _ (fun i _ => ?_)
  exact (posSemidef_sigMat θ R i).smul (by positivity)

end SingleWeight

end StackedSVD
