/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Eigen

/-!
# `specInvTop` under a simultaneous relabeling of the rows and the columns

STATUS 2026-09-01: proved, 0 `sorry` (task B4 of `notes/archive/rankr_plan_B.md`).

`specInvTop A hA r` is written with a **choice** of eigenbasis
(`Matrix.IsHermitian.eigenvectorBasis`), so it is not obvious that a relabeling of the index
type commutes with it: the relabeled matrix carries its own eigenbasis, which is not the
relabeled one. The value does not depend on the choice, because the coefficient of the sum
reads the eigenvalue and not the index. This file proves that, and it then moves `specInvTop`
and the trace functional `tr(Bᵀ (specInvTop A r) B)` along an arbitrary `Equiv` of the index
type.

## Content

1. `specFun A hA f = ∑_i f(λ_i) q_i q_iᵀ`, the spectral function of `f` at `A`. `specInvTop`
   is the case `f = 1{· ∈ topEigSet A hA r} (·)⁻¹` (`specInvTop_eq_specFun`, `rfl`).
2. `specFun_eq_of_eigenfamily`: every family `u` of eigenvectors of `A` that resolves the
   identity gives the same matrix, `∑_i f(lam_i) u_i u_iᵀ = specFun A hA f`. This is the basis
   independence. The proof reads both matrices on the canonical eigenbasis: an eigenvector
   `u_i` that is not orthogonal to `q_j` carries the eigenvalue `λ_j`, so its coefficient is
   `f(λ_j)`, and the family then resolves `q_j` itself.
3. `specFun_submatrix`, `topEigSet_submatrix`, `specInvTop_submatrix`: the relabeled
   eigenbasis resolves the identity, which gives the first; `Matrix.charpoly_reindex` with
   `eigenvalues₀_eq_of_charpoly_eq` gives the second; the third is their composition.
4. `trace_specInvTop_submatrix` and the congruence form `trace_specInvTop_congr_submatrix`,
   which `RankR/General.lean` consumes for `limitRG_one_eq_limitR`.

No statement here reads an eigengap or a positivity hypothesis. Every one holds for an
arbitrary real symmetric matrix, an arbitrary `r` and an arbitrary `f`.
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. The spectral function of a real symmetric matrix -/

section SpecFun

variable {p : ℕ}

/-- `f(A) = ∑_i f(λ_i) q_i q_iᵀ`, the spectral function of `f` at `A`, written in the
eigenbasis of `A`. -/
noncomputable def specFun (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (f : ℝ → ℝ) :
    Matrix (Fin p) (Fin p) ℝ :=
  ∑ i : Fin p, f (hA.eigenvalues i) •
    Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis i))
      (WithLp.ofLp (hA.eigenvectorBasis i))

/-- `specInvTop` is the spectral function of `1{· ∈ topEigSet A hA r} (·)⁻¹`. -/
theorem specInvTop_eq_specFun (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) :
    specInvTop A hA r = specFun A hA (Set.indicator (topEigSet A hA r) fun t => t⁻¹) := rfl

/-- Two matrices with the same action on every vector are equal. -/
private theorem eq_of_mulVec_eq {S T : Matrix (Fin p) (Fin p) ℝ}
    (h : ∀ y : Fin p → ℝ, S *ᵥ y = T *ᵥ y) : S = T := by
  ext a b
  have h1 := congrFun (h (Pi.single b 1)) a
  simpa [Matrix.mulVec, dotProduct, Pi.single_apply] using h1

/-- The action of `∑_i c_i u_i u_iᵀ` on a vector. -/
private theorem sum_vecMulVec_mulVec {ι : Type*} [Fintype ι] (c : ι → ℝ)
    (u : ι → EuclideanSpace ℝ (Fin p)) (y : Fin p → ℝ) :
    (∑ i, c i • Matrix.vecMulVec (WithLp.ofLp (u i)) (WithLp.ofLp (u i))) *ᵥ y
      = ∑ i, (c i * (WithLp.ofLp (u i) ⬝ᵥ y)) • WithLp.ofLp (u i) := by
  rw [Matrix.sum_mulVec]
  refine Finset.sum_congr rfl fun i _ => ?_
  funext j
  simp only [Matrix.mulVec, dotProduct, Matrix.smul_apply, Matrix.vecMulVec_apply,
    Pi.smul_apply, smul_eq_mul, Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- A real symmetric matrix moves across a dot product. -/
private theorem dotProduct_mulVec_comm {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (x y : Fin p → ℝ) : (A *ᵥ x) ⬝ᵥ y = x ⬝ᵥ (A *ᵥ y) := by
  have hAT : Aᵀ = A := by
    have h := hA
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
  rw [Matrix.dotProduct_mulVec]
  conv_rhs => rw [← hAT]
  rw [Matrix.vecMul_transpose]

/-- The canonical eigenbasis resolves the identity: `∑_i q_i (q_iᵀ y) = y`. -/
private theorem eigenvectorBasis_resolves {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (y : Fin p → ℝ) :
    ∑ i, (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ y) • WithLp.ofLp (hA.eigenvectorBasis i)
      = y := by
  have h := hA.eigenvectorBasis.sum_repr' (WithLp.toLp 2 y)
  have h2 : ∑ i, WithLp.ofLp
      (⟪hA.eigenvectorBasis i, WithLp.toLp 2 y⟫_ℝ • hA.eigenvectorBasis i) = y := by
    rw [← WithLp.ofLp_sum, h]
  simpa only [WithLp.ofLp_smul, real_inner_eq_dotProduct, WithLp.ofLp_toLp] using h2

/-- `∑_i f(lam_i) u_i u_iᵀ` sends the `j`-th canonical eigenvector to `f(λ_j)` times itself,
for every family `u` of eigenvectors of `A` that resolves the identity. -/
private theorem sum_vecMulVec_apply_eigvec {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (f : ℝ → ℝ) {ι : Type*} [Fintype ι] (u : ι → EuclideanSpace ℝ (Fin p)) (lam : ι → ℝ)
    (hres : ∀ y : Fin p → ℝ, ∑ i, (WithLp.ofLp (u i) ⬝ᵥ y) • WithLp.ofLp (u i) = y)
    (heig : ∀ i, A *ᵥ WithLp.ofLp (u i) = lam i • WithLp.ofLp (u i)) (j : Fin p) :
    (∑ i, f (lam i) • Matrix.vecMulVec (WithLp.ofLp (u i)) (WithLp.ofLp (u i))) *ᵥ
        WithLp.ofLp (hA.eigenvectorBasis j)
      = f (hA.eigenvalues j) • WithLp.ofLp (hA.eigenvectorBasis j) := by
  have hw : A *ᵥ WithLp.ofLp (hA.eigenvectorBasis j)
      = hA.eigenvalues j • WithLp.ofLp (hA.eigenvectorBasis j) :=
    hA.mulVec_eigenvectorBasis j
  have hcoef : ∀ i : ι,
      (f (lam i) * (WithLp.ofLp (u i) ⬝ᵥ WithLp.ofLp (hA.eigenvectorBasis j))) •
          WithLp.ofLp (u i)
        = f (hA.eigenvalues j) •
            ((WithLp.ofLp (u i) ⬝ᵥ WithLp.ofLp (hA.eigenvectorBasis j)) •
              WithLp.ofLp (u i)) := by
    intro i
    rw [smul_smul]
    rcases eq_or_ne (WithLp.ofLp (u i) ⬝ᵥ WithLp.ofLp (hA.eigenvectorBasis j)) 0 with h | h
    · rw [h, mul_zero, mul_zero]
    · have h1 := dotProduct_mulVec_comm hA (WithLp.ofLp (u i))
        (WithLp.ofLp (hA.eigenvectorBasis j))
      rw [heig i, hw, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul] at h1
      rw [mul_right_cancel₀ h h1]
  rw [sum_vecMulVec_mulVec, Finset.sum_congr rfl fun i (_ : i ∈ Finset.univ) => hcoef i,
    ← Finset.smul_sum, hres]

/-- **Basis independence of a spectral function.** Every family `u` of eigenvectors of `A`
that resolves the identity produces `specFun A hA f`, whatever the eigenvalue list `lam` and
the index type are. -/
theorem specFun_eq_of_eigenfamily {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (f : ℝ → ℝ) {ι : Type*} [Fintype ι] (u : ι → EuclideanSpace ℝ (Fin p)) (lam : ι → ℝ)
    (hres : ∀ y : Fin p → ℝ, ∑ i, (WithLp.ofLp (u i) ⬝ᵥ y) • WithLp.ofLp (u i) = y)
    (heig : ∀ i, A *ᵥ WithLp.ofLp (u i) = lam i • WithLp.ofLp (u i)) :
    ∑ i, f (lam i) • Matrix.vecMulVec (WithLp.ofLp (u i)) (WithLp.ofLp (u i))
      = specFun A hA f := by
  simp only [specFun]
  refine eq_of_mulVec_eq fun y => ?_
  have hy := eigenvectorBasis_resolves hA y
  have hL := sum_vecMulVec_apply_eigvec hA f u lam hres heig
  have hR := sum_vecMulVec_apply_eigvec hA f (fun i => hA.eigenvectorBasis i) hA.eigenvalues
    (eigenvectorBasis_resolves hA) hA.mulVec_eigenvectorBasis
  rw [← hy, Matrix.mulVec_sum, Matrix.mulVec_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Matrix.mulVec_smul, Matrix.mulVec_smul, hL j, hR j]

end SpecFun

/-! ### 2. A simultaneous relabeling of the rows and the columns -/

section Reindex

variable {p q : ℕ}

/-- The vector `x` relabeled along `e`: entry `a` reads entry `e a` of `x`. -/
private def tvec (e : Fin q ≃ Fin p) (x : EuclideanSpace ℝ (Fin p)) :
    EuclideanSpace ℝ (Fin q) :=
  WithLp.toLp 2 fun a => WithLp.ofLp x (e a)

private theorem ofLp_tvec (e : Fin q ≃ Fin p) (x : EuclideanSpace ℝ (Fin p)) (a : Fin q) :
    WithLp.ofLp (tvec e x) a = WithLp.ofLp x (e a) := rfl

/-- The relabeled eigenbasis still resolves the identity. -/
private theorem tvec_resolves {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (e : Fin q ≃ Fin p) (y : Fin q → ℝ) :
    ∑ i : Fin p, (WithLp.ofLp (tvec e (hA.eigenvectorBasis i)) ⬝ᵥ y) •
      WithLp.ofLp (tvec e (hA.eigenvectorBasis i)) = y := by
  funext a
  have hz := congrFun (eigenvectorBasis_resolves hA fun c => y (e.symm c)) (e a)
  have hdot : ∀ i : Fin p,
      WithLp.ofLp (tvec e (hA.eigenvectorBasis i)) ⬝ᵥ y
        = WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ fun c => y (e.symm c) := by
    intro i
    change ∑ b : Fin q, WithLp.ofLp (tvec e (hA.eigenvectorBasis i)) b * y b
      = ∑ c : Fin p, WithLp.ofLp (hA.eigenvectorBasis i) c * y (e.symm c)
    rw [← Equiv.sum_comp e fun c => WithLp.ofLp (hA.eigenvectorBasis i) c * y (e.symm c)]
    exact Finset.sum_congr rfl fun b _ => by rw [ofLp_tvec, Equiv.symm_apply_apply]
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, ofLp_tvec] at hz ⊢
  simp only [hdot]
  rw [hz, Equiv.symm_apply_apply]

/-- The relabeled eigenvectors are eigenvectors of the relabeled matrix, with the same
eigenvalues. -/
private theorem tvec_eig {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (e : Fin q ≃ Fin p) (i : Fin p) :
    A.submatrix e e *ᵥ WithLp.ofLp (tvec e (hA.eigenvectorBasis i))
      = hA.eigenvalues i • WithLp.ofLp (tvec e (hA.eigenvectorBasis i)) := by
  funext a
  have h := congrFun (hA.mulVec_eigenvectorBasis i) (e a)
  have hsum : ∑ b : Fin q, A (e a) (e b) * WithLp.ofLp (hA.eigenvectorBasis i) (e b)
      = ∑ c : Fin p, A (e a) c * WithLp.ofLp (hA.eigenvectorBasis i) c :=
    Equiv.sum_comp e fun c => A (e a) c * WithLp.ofLp (hA.eigenvectorBasis i) c
  simp only [Matrix.mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at h ⊢
  simp only [Matrix.submatrix_apply, ofLp_tvec]
  rw [hsum, h]

/-- The spectral function commutes with a simultaneous relabeling of the rows and the
columns. -/
theorem specFun_submatrix (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (e : Fin q ≃ Fin p) (hA' : (A.submatrix e e).IsHermitian) (f : ℝ → ℝ) :
    specFun (A.submatrix e e) hA' f = (specFun A hA f).submatrix e e := by
  rw [← specFun_eq_of_eigenfamily hA' f (fun i => tvec e (hA.eigenvectorBasis i))
    hA.eigenvalues (tvec_resolves hA e) (tvec_eig hA e)]
  ext a b
  simp only [specFun, Matrix.sum_apply, Matrix.smul_apply, Matrix.vecMulVec_apply,
    Matrix.submatrix_apply, smul_eq_mul, ofLp_tvec]

/-- The top-`r` eigenvalue set does not see a relabeling: the characteristic polynomial does
not, so the sorted eigenvalue list does not. -/
theorem topEigSet_submatrix (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (e : Fin q ≃ Fin p) (hA' : (A.submatrix e e).IsHermitian) (r : ℕ) :
    topEigSet (A.submatrix e e) hA' r = topEigSet A hA r := by
  have hqp : q = p := by simpa using Fintype.card_congr e
  subst hqp
  have hchar : (A.submatrix e e).charpoly = A.charpoly := by
    have h := Matrix.charpoly_reindex e.symm A
    rwa [Matrix.reindex_apply, Equiv.symm_symm] at h
  have heq : hA'.eigenvalues₀ = hA.eigenvalues₀ :=
    eigenvalues₀_eq_of_charpoly_eq hA' hA hchar
  simp only [topEigSet, heq]

/-- **`specInvTop` under a simultaneous relabeling.** -/
theorem specInvTop_submatrix (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (e : Fin q ≃ Fin p) (hA' : (A.submatrix e e).IsHermitian) (r : ℕ) :
    specInvTop (A.submatrix e e) hA' r = (specInvTop A hA r).submatrix e e := by
  rw [specInvTop_eq_specFun (A.submatrix e e) hA' r, specInvTop_eq_specFun A hA r,
    topEigSet_submatrix A hA e hA' r, specFun_submatrix A hA e hA']

/-- The trace functional `tr(Bᵀ (specInvTop A r) B)` does not see a simultaneous relabeling of
the rows and the columns of `A` together with the rows of `B`. The matrix arguments come with
their equations, so that the caller does not have to rewrite under the `IsHermitian` proof of
`specInvTop`. -/
theorem trace_specInvTop_congr_submatrix {s : ℕ} (A : Matrix (Fin p) (Fin p) ℝ)
    (hA : A.IsHermitian) (B : Matrix (Fin p) (Fin s) ℝ) (e : Fin q ≃ Fin p)
    (A' : Matrix (Fin q) (Fin q) ℝ) (hA' : A'.IsHermitian) (B' : Matrix (Fin q) (Fin s) ℝ)
    (hAeq : A' = A.submatrix e e) (hBeq : B' = B.submatrix e id) (r : ℕ) :
    Matrix.trace (B'ᵀ * specInvTop A' hA' r * B')
      = Matrix.trace (Bᵀ * specInvTop A hA r * B) := by
  subst hAeq
  subst hBeq
  rw [specInvTop_submatrix A hA e hA' r]
  simp only [Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv, Matrix.submatrix_id_id]

/-- `trace_specInvTop_congr_submatrix` with the two equations discharged by `rfl`. -/
theorem trace_specInvTop_submatrix {s : ℕ} (A : Matrix (Fin p) (Fin p) ℝ)
    (hA : A.IsHermitian) (B : Matrix (Fin p) (Fin s) ℝ) (e : Fin q ≃ Fin p)
    (hA' : (A.submatrix e e).IsHermitian) (r : ℕ) :
    Matrix.trace ((B.submatrix e id)ᵀ * specInvTop (A.submatrix e e) hA' r *
        B.submatrix e id)
      = Matrix.trace (Bᵀ * specInvTop A hA r * B) :=
  trace_specInvTop_congr_submatrix A hA B e _ hA' _ rfl rfl r

end Reindex

end StackedSVD
