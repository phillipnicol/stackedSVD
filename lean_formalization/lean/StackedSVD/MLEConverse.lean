/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLE

/-!
# `app:wstacksvd_mle`: the converse, every maximizer sits in the top eigenspace

`StackedSVD/MLE.lean` proves one direction of the paper's `v̂_MLE = v_max(...)`
(`main_paper.tex:1700`): a unit vector of the top eigenspace of the weighted stack Gram matrix
maximizes the marginal log-likelihood `ℓ` (`MultiTableModel.mleLogLik_le_of_mem_topSpace`).
Modeling choice 5 of `notes/archive/mle_identity.md` left the converse open. This file adds it.

## Content

1. `inner_toOp_self_eq_lamMax_iff`: the equality case of the Rayleigh bound
   `⟪A x, x⟫ ≤ λ_max ‖x‖²`. For a real symmetric matrix `A` the equality holds exactly when
   `x` lies in `topSpace A`. This lemma takes no model and belongs in `Spectral.lean`; it
   lives here because the audit of that file runs in parallel.
2. `MultiTableModel.mem_topSpace_of_mleLogLik_max`: every unit maximizer of `ℓ` lies in the
   top eigenspace of `stackGramW (optWstack θ c)`.
3. `MultiTableModel.mleLogLik_max_iff_mem_topSpace`: the two-sided form. A unit vector
   maximizes `ℓ` over unit vectors if and only if it lies in the top eigenspace. Together with
   the projector convention of `notes/INTERFACES.md` this is the paper's `argmax = v_max`.

The route of item 1 reads `x` in the sorted eigenbasis of `A`. Equality forces
`∑_i (λ_max − λ_i) a_i² = 0` with every term nonnegative, so `a_i = 0` at each `λ_i < λ_max`.
Item 2 turns the maximizer hypothesis into the Rayleigh equality with `vMax A` as the test
vector, then applies item 1. No probability and no `N → ∞` enters.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

open Matrix

/-! ### The equality case of the Rayleigh bound -/

section RayleighEquality

variable {d : ℕ}

/-- Every eigenvalue is at most `lamMax`. Copy of the private `eigenvalues₀_le_lamMax` of
`Spectral.lean`; the index type is nonempty, so `0 < d`. -/
private theorem eigenvalues₀_le_lamMax' (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (i : Fin (Fintype.card (Fin d))) : hA.eigenvalues₀ i ≤ lamMax A hA := by
  have hd : 0 < d := by
    have h1 : (0 : ℕ) < Fintype.card (Fin d) := Nat.lt_of_le_of_lt (Nat.zero_le i.val) i.isLt
    simpa using h1
  rw [lamMax, dif_pos hd]
  exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

/-- The quadratic form in the sorted eigenbasis: `⟪A x, x⟫ = ∑_i λ_i ⟪e_i, x⟫²`. This is the
public twin of the private `quadForm_pow` of `Spectral.lean` at the power `1`. -/
private theorem quadForm_eigen (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (x : EuclideanSpace ℝ (Fin d)) :
    ⟪toOp A x, x⟫_ℝ
      = ∑ i, hA.eigenvalues₀ i
          * ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ ^ 2 := by
  rw [real_inner_comm,
    inner_eq_sum_inner ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) x (toOp A x)]
  refine Finset.sum_congr rfl fun i _ => ?_
  have h1 : ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, toOp A x⟫_ℝ
      = ⟪toOp A ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i), x⟫_ℝ :=
    (symmOp hA _ _).symm
  rw [h1, apply_eigvec hA i, real_inner_smul_left]
  ring

/-- **Equality case of the Rayleigh bound.** For a real symmetric matrix `A`, the bound
`⟪A x, x⟫ ≤ λ_max ‖x‖²` of `inner_toOp_self_le` is an equality exactly when `x` lies in the top
eigenspace. Route: expand `x` in the sorted eigenbasis; equality gives
`∑_i (λ_max − λ_i) ⟪e_i, x⟫² = 0` with every term nonnegative, so the coordinate vanishes at
every eigenvalue below `λ_max`. This lemma takes no model and belongs in `Spectral.lean`. -/
theorem inner_toOp_self_eq_lamMax_iff (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (x : EuclideanSpace ℝ (Fin d)) :
    ⟪toOp A x, x⟫_ℝ = lamMax A hA * ‖x‖ ^ 2 ↔ x ∈ topSpace A hA := by
  constructor
  · intro heq
    have hsum : ∑ i, (lamMax A hA - hA.eigenvalues₀ i)
        * ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ ^ 2 = 0 := by
      have hexp : ∑ i, (lamMax A hA - hA.eigenvalues₀ i)
          * ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ ^ 2
          = lamMax A hA * ‖x‖ ^ 2 - ⟪toOp A x, x⟫_ℝ := by
        rw [norm_sq_eq_sum_inner ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) x,
          quadForm_eigen A hA x, Finset.mul_sum, ← Finset.sum_sub_distrib]
        exact Finset.sum_congr rfl fun i _ => by ring
      rw [hexp, heq]
      ring
    have hnn : ∀ i ∈ Finset.univ, 0 ≤ (lamMax A hA - hA.eigenvalues₀ i)
        * ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ ^ 2 := by
      intro i _
      have h1 := eigenvalues₀_le_lamMax' A hA i
      have h2 : 0 ≤ lamMax A hA - hA.eigenvalues₀ i := by linarith
      positivity
    have hterm := (Finset.sum_eq_zero_iff_of_nonneg hnn).mp hsum
    rw [← OrthonormalBasis.sum_repr' ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) x]
    refine Submodule.sum_mem _ fun i _ => ?_
    by_cases hli : hA.eigenvalues₀ i = lamMax A hA
    · refine Submodule.smul_mem _ _ ?_
      rw [topSpace_eq_eigenspace, Module.End.mem_eigenspace_iff, apply_eigvec hA i, hli]
    · have hzero := hterm i (Finset.mem_univ i)
      have hcoef : ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, x⟫_ℝ = 0 := by
        rcases mul_eq_zero.mp hzero with h | h
        · exact absurd (sub_eq_zero.mp h).symm hli
        · exact pow_eq_zero_iff (two_ne_zero) |>.mp h
      rw [hcoef, zero_smul]
      exact Submodule.zero_mem _
  · intro hx
    rw [toOp_of_mem_topSpace hx, real_inner_smul_left, real_inner_self_eq_norm_sq]

end RayleighEquality

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-- `app:wstacksvd_mle`, the converse of `mleLogLik_le_of_mem_topSpace`: every unit maximizer
of the marginal log-likelihood lies in the top eigenspace of the weighted stack Gram matrix at
the paper's optimal weights (`main_paper.tex:1700`). Route: `mleLogLik_le_iff` turns the
maximizer hypothesis into the Rayleigh equality (the top eigenvector `vMax` is the test vector
for the lower bound, `inner_toOp_self_le` for the upper one), and
`inner_toOp_self_eq_lamMax_iff` reads the equality case. -/
theorem mem_topSpace_of_mleLogLik_max (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (N : ℕ) (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1)
    (hmax : ∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
      m.mleLogLik c N ω v ≤ m.mleLogLik c N ω (WithLp.ofLp x)) :
    x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (m.isHermitian_stackGramW _ N ω) := by
  set A := m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω with hA
  have hAh : A.IsHermitian := m.isHermitian_stackGramW _ N ω
  have hinner : ∀ y : EuclideanSpace ℝ (Fin (d N)),
      ⟪toOp A y, y⟫_ℝ = WithLp.ofLp y ⬝ᵥ A *ᵥ WithLp.ofLp y := fun y => by
    rw [real_inner_eq_dotProduct, dotProduct_comm]
    rfl
  have hxdot : WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = 1 := by
    have h := real_inner_eq_dotProduct x x
    rw [real_inner_self_eq_norm_sq, hx] at h
    simpa using h.symm
  have hd : 0 < d N := pos_of_dotProduct_self_eq_one hxdot
  have hudot : WithLp.ofLp (vMax A hAh) ⬝ᵥ WithLp.ofLp (vMax A hAh) = 1 := by
    have h := real_inner_eq_dotProduct (vMax A hAh) (vMax A hAh)
    rw [real_inner_self_eq_norm_sq, norm_vMax hd A hAh] at h
    simpa using h.symm
  have huray : WithLp.ofLp (vMax A hAh) ⬝ᵥ A *ᵥ WithLp.ofLp (vMax A hAh) = lamMax A hAh := by
    rw [← hinner (vMax A hAh), toOp_of_mem_topSpace (mem_topSpace_vMax hd A hAh),
      real_inner_smul_left, real_inner_self_eq_norm_sq, norm_vMax hd A hAh]
    ring
  have hge : lamMax A hAh ≤ WithLp.ofLp x ⬝ᵥ A *ᵥ WithLp.ofLp x := by
    have h := hmax (WithLp.ofLp (vMax A hAh)) hudot
    rw [m.mleLogLik_le_iff c hc N ω (WithLp.ofLp (vMax A hAh)) (WithLp.ofLp x) hudot hxdot,
      ← hA, huray] at h
    exact h
  have hle : WithLp.ofLp x ⬝ᵥ A *ᵥ WithLp.ofLp x ≤ lamMax A hAh := by
    have h := inner_toOp_self_le A hAh x
    rw [hinner x, hx] at h
    simpa using h
  refine (inner_toOp_self_eq_lamMax_iff A hAh x).mp ?_
  rw [hinner x, hx]
  simp only [one_pow, mul_one]
  linarith

omit [NeZero M] in
/-- `app:wstacksvd_mle`, the two-sided statement: a unit vector maximizes the marginal
log-likelihood over unit vectors if and only if it lies in the top eigenspace of
`stackGramW (optWstack θ c)`. The forward direction is `mem_topSpace_of_mleLogLik_max`, the
backward one is `mleLogLik_le_of_mem_topSpace` of `MLE.lean`. Read with the projector
convention of `notes/INTERFACES.md`, this is the paper's
`v̂_MLE = v_max(∑_i θ_i²/(c_i+θ_i²) X_iᵀ X_i)`. -/
theorem mleLogLik_max_iff_mem_topSpace [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (N : ℕ) (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
        m.mleLogLik c N ω v ≤ m.mleLogLik c N ω (WithLp.ofLp x))
      ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
          (m.isHermitian_stackGramW _ N ω) :=
  ⟨fun hmax => m.mem_topSpace_of_mleLogLik_max c hc N ω x hx hmax,
    fun hxtop v hv => m.mleLogLik_le_of_mem_topSpace c hc N ω v hv x hx hxtop⟩

end MultiTableModel

end StackedSVD
