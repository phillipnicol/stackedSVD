/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVDWeighted

/-!
# `app:wstacksvd_mle`: the marginal log-likelihood is the weighted stackSVD objective

The review note is `notes/archive/mle_identity.md` (task D5). The paper's appendix
`app:wstacksvd_mle` (`main_paper.tex:1630`) marginalizes over `u_i`, gets rows of `X_i` that
are independent Gaussian vectors with covariance `Σ_i(v) = (1/d)(I + (θ_i²/c_i) v vᵀ)`, and
then reduces the marginal log-likelihood to the weighted stackSVD objective.

## Content

The Gaussian marginalization itself is not formalized (modeling choice 1 of the note): `ℓ`
enters as the definition `MultiTableModel.mleLogLik`, the paper's displayed formula. What the
file proves is the paper's chain of equalities that follows it.

1. `mleCov d a v = (1/d)(I + a v vᵀ)`, the paper's `Σ_i(v)` at `a = θ_i²/c_i`.
2. `det_mleCov`: `det Σ = d^{-d}(1 + a)` at `‖v‖ = 1`, so the determinant does not see the
   direction of `v`.
3. `inv_mleCov`: Sherman-Morrison, `Σ⁻¹ = d(I - (a/(1+a)) v vᵀ)`.
4. `trace_inv_mleCov_mul`: `tr(Σ⁻¹ G) = d tr G - d (a/(1+a)) vᵀ G v` for every square `G`.
5. `MultiTableModel.mleLogLik_eq`: on unit vectors `ℓ(v)` is `mleConst` plus `d/2` times the
   Rayleigh form of `stackGramW (optWstack θ c)`. The weight identity behind it is
   `optWstack θ c i ^ 2 = θ_i²/(θ_i² + c_i) = a_i/(1 + a_i)`.
6. `MultiTableModel.mleLogLik_le_iff` and `mleLogLik_le_of_mem_topSpace`: over unit vectors
   `ℓ` and that Rayleigh form have the same order, so every unit vector of the top eigenspace
   of the weighted stack Gram matrix maximizes `ℓ`. This is the paper's
   `v̂_MLE = v_max(∑_i θ_i²/(c_i+θ_i²) X_iᵀ X_i)`, the weighted stackSVD estimator.

No probability and no `N → ∞` enters: every statement is at a fixed `N` and a fixed `ω`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

open Matrix

/-! ### The one-table marginal covariance -/

/-- `Σ(a, v) = (1/d)(I + a v vᵀ)`, the paper's `Σ_i(v)` with `a = θ_i²/c_i`
(`main_paper.tex:1654`). -/
noncomputable def mleCov (d : ℕ) (a : ℝ) (v : Fin d → ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  (1 / (d : ℝ)) • (1 + a • vecMulVec v v)

/-- A unit vector needs a nonempty index type: `v ⬝ᵥ v = 1` forces `0 < d`. -/
theorem pos_of_dotProduct_self_eq_one {d : ℕ} {v : Fin d → ℝ} (hv : v ⬝ᵥ v = 1) : 0 < d := by
  rcases Nat.eq_zero_or_pos d with rfl | h
  · simp [dotProduct] at hv
  · exact h

/-- The determinant of `Σ(a, v)` at a unit `v` is `d^{-d}(1 + a)`, which does not see the
direction of `v` (`main_paper.tex:1668`). -/
theorem det_mleCov (d : ℕ) (a : ℝ) (v : Fin d → ℝ) (hv : v ⬝ᵥ v = 1) :
    (mleCov d a v).det = (1 / (d : ℝ)) ^ d * (1 + a) := by
  rw [mleCov, Matrix.det_smul, Fintype.card_fin]
  congr 1
  rw [← Matrix.smul_vecMulVec, Matrix.vecMulVec_eq (Fin 1),
    Matrix.det_one_add_replicateCol_mul_replicateRow]
  rw [dotProduct_smul, hv]
  simp

/-- Sherman-Morrison for `Σ(a, v)` at a unit `v`: `Σ⁻¹ = d(I - (a/(1+a)) v vᵀ)`
(`main_paper.tex:1674`). The hypothesis `0 ≤ a` gives `1 + a ≠ 0`. -/
theorem inv_mleCov (d : ℕ) (a : ℝ) (v : Fin d → ℝ) (hv : v ⬝ᵥ v = 1) (ha : 0 ≤ a) :
    (mleCov d a v)⁻¹ = (d : ℝ) • (1 - (a / (1 + a)) • vecMulVec v v) := by
  have h1a : (1 : ℝ) + a ≠ 0 := by positivity
  have hd : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (pos_of_dotProduct_self_eq_one hv).ne'
  refine Matrix.inv_eq_right_inv ?_
  rw [mleCov, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [show (1 / (d : ℝ)) * (d : ℝ) = 1 by field_simp, one_smul]
  rw [Matrix.add_mul, Matrix.one_mul, Matrix.mul_sub, Matrix.mul_one, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.vecMulVec_mul_vecMulVec, hv, one_smul, smul_smul]
  match_scalars <;> (field_simp; try ring)

/-- The trace that the log-likelihood needs: `tr(Σ(a, v)⁻¹ G) = d tr G - d (a/(1+a)) vᵀ G v`
for every square `G` (`main_paper.tex:1680`). -/
theorem trace_inv_mleCov_mul (d : ℕ) (a : ℝ) (v : Fin d → ℝ) (hv : v ⬝ᵥ v = 1) (ha : 0 ≤ a)
    (G : Matrix (Fin d) (Fin d) ℝ) :
    ((mleCov d a v)⁻¹ * G).trace
      = (d : ℝ) * G.trace - (d : ℝ) * (a / (1 + a)) * (v ⬝ᵥ G *ᵥ v) := by
  have hPG : (vecMulVec v v * G).trace = v ⬝ᵥ G *ᵥ v := by
    rw [Matrix.trace_mul_comm, Matrix.mul_vecMulVec, Matrix.trace_vecMulVec, dotProduct_comm]
  rw [inv_mleCov d a v hv ha, Matrix.smul_mul, Matrix.trace_smul, Matrix.sub_mul,
    Matrix.one_mul, Matrix.trace_sub, Matrix.smul_mul, Matrix.trace_smul, hPG]
  simp only [smul_eq_mul]
  ring

namespace Scalars

/-- The paper's weight identity: `w_i⋆² = θ_i²/(θ_i² + c_i) = a_i/(1 + a_i)` at
`a_i = θ_i²/c_i`. It turns the coefficient of `trace_inv_mleCov_mul` into the square of the
optimal weight of `thm:stacksvd_weighted`. -/
theorem optWstack_sq {M : ℕ} (θ c : Fin M → ℝ) (hc : ∀ i, 0 < c i) (i : Fin M) :
    optWstack θ c i ^ 2 = (θ i ^ 2 / c i) / (1 + θ i ^ 2 / c i) := by
  have hci := hc i
  have hs : (0 : ℝ) < θ i ^ 2 + c i := add_pos_of_nonneg_of_pos (sq_nonneg _) hci
  have hcne : c i ≠ 0 := hci.ne'
  have hsne : θ i ^ 2 + c i ≠ 0 := hs.ne'
  rw [optWstack, div_pow, Real.sq_sqrt hs.le]
  field_simp
  ring

end Scalars

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### The weighted stack Gram matrix as a weighted sum of table Gram matrices -/

/-- The Gram matrix of the weighted stack is `∑_i w_i² X_iᵀ X_i`, the matrix of
`eq:weightedStackSVD_form`. This is `stackGramW_apply` read as a matrix identity. -/
theorem stackGramW_eq_sum (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.stackGramW w N ω = ∑ i, w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) := by
  ext k l
  change (((m.stackW w).X N ω)ᵀ * (m.stackW w).X N ω) k l = _
  rw [m.stackGramW_apply w N ω k l]
  simp [Matrix.sum_apply, Matrix.mul_apply]

/-- Rayleigh form of the weighted stack Gram matrix: `vᵀ (∑_i w_i² X_iᵀ X_i) v`. -/
theorem dotProduct_stackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (v : Fin (d N) → ℝ) :
    v ⬝ᵥ m.stackGramW w N ω *ᵥ v
      = ∑ i, w i ^ 2 * (v ⬝ᵥ (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) *ᵥ v) := by
  rw [m.stackGramW_eq_sum w N ω, Matrix.sum_mulVec, dotProduct_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]

/-! ### The marginal log-likelihood -/

/-- The paper's marginal log-likelihood `ℓ(v)` (`main_paper.tex:1662`), taken as a definition:
`-(1/2) ∑_i (n_i log det Σ_i(v) + tr(Σ_i(v)⁻¹ X_iᵀX_i))` with `Σ_i(v) = mleCov d (θ_i²/c_i) v`.
The Gaussian marginalization that produces this formula is `thm_wstacksvd_mle_marginal`
(`MLEMarginal/Main.lean`, 2026-09-03; modeling choice 1 of `notes/archive/mle_identity.md`
recorded the gap before that). -/
noncomputable def mleLogLik (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (v : Fin (d N) → ℝ) : ℝ :=
  -(1 / 2) * ∑ i,
    ((n i N : ℝ) * Real.log (mleCov (d N) ((m.tbl i).θ ^ 2 / c i) v).det
      + ((mleCov (d N) ((m.tbl i).θ ^ 2 / c i) v)⁻¹
          * (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)).trace)

/-- The `v`-free part of `ℓ`. -/
noncomputable def mleConst (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    ℝ :=
  -(1 / 2) * ∑ i,
    ((n i N : ℝ) * Real.log ((1 / (d N : ℝ)) ^ (d N) * (1 + (m.tbl i).θ ^ 2 / c i))
      + (d N : ℝ) * (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω).trace)

omit [NeZero M] in
/-- `app:wstacksvd_mle`, the identity: on unit vectors `ℓ(v)` is a constant plus `d/2` times
the Rayleigh form of the weighted stack Gram matrix at the paper's optimal weights
(`main_paper.tex:1690`). -/
theorem mleLogLik_eq [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (N : ℕ) (ω : Ω N) (v : Fin (d N) → ℝ) (hv : v ⬝ᵥ v = 1) :
    m.mleLogLik c N ω v
      = m.mleConst c N ω
        + (d N : ℝ) / 2
          * (v ⬝ᵥ m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω *ᵥ v) := by
  set w := Scalars.optWstack (fun i => (m.tbl i).θ) c with hw
  have hwsq : ∀ i, w i ^ 2 = ((m.tbl i).θ ^ 2 / c i) / (1 + (m.tbl i).θ ^ 2 / c i) := fun i =>
    Scalars.optWstack_sq (fun i => (m.tbl i).θ) c hc i
  have hsum : ∑ i,
      ((n i N : ℝ) * Real.log (mleCov (d N) ((m.tbl i).θ ^ 2 / c i) v).det
        + ((mleCov (d N) ((m.tbl i).θ ^ 2 / c i) v)⁻¹
            * (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)).trace)
      = (∑ i,
          ((n i N : ℝ) * Real.log ((1 / (d N : ℝ)) ^ (d N) * (1 + (m.tbl i).θ ^ 2 / c i))
            + (d N : ℝ) * (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω).trace))
        - (d N : ℝ)
            * ∑ i, w i ^ 2 * (v ⬝ᵥ (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) *ᵥ v) := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have ha : (0 : ℝ) ≤ (m.tbl i).θ ^ 2 / c i := div_nonneg (sq_nonneg _) (hc i).le
    rw [det_mleCov (d N) _ v hv, trace_inv_mleCov_mul (d N) _ v hv ha, hwsq i]
    ring
  rw [mleLogLik, mleConst, dotProduct_stackGramW, hsum]
  ring

/-- `app:wstacksvd_mle`, the order statement: on unit vectors `ℓ` and the Rayleigh form of
`stackGramW (optWstack θ c)` have the same order. -/
theorem mleLogLik_le_iff (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (N : ℕ) (ω : Ω N) (v v' : Fin (d N) → ℝ) (hv : v ⬝ᵥ v = 1) (hv' : v' ⬝ᵥ v' = 1) :
    m.mleLogLik c N ω v ≤ m.mleLogLik c N ω v'
      ↔ v ⬝ᵥ m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω *ᵥ v
        ≤ v' ⬝ᵥ m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω *ᵥ v' := by
  have hd : (0 : ℝ) < (d N : ℝ) / 2 := by
    have := pos_of_dotProduct_self_eq_one hv
    positivity
  rw [m.mleLogLik_eq c hc N ω v hv, m.mleLogLik_eq c hc N ω v' hv', add_le_add_iff_left]
  exact ⟨fun h => le_of_mul_le_mul_left h hd, fun h => mul_le_mul_of_nonneg_left h hd.le⟩

/-- `app:wstacksvd_mle`, the conclusion: every unit vector of the top eigenspace of the
weighted stack Gram matrix maximizes `ℓ` over unit vectors. So the weighted stackSVD
estimator at the optimal weights is a marginal MLE (`main_paper.tex:1700`). -/
theorem mleLogLik_le_of_mem_topSpace (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (N : ℕ) (ω : Ω N) (v : Fin (d N) → ℝ) (hv : v ⬝ᵥ v = 1)
    (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1)
    (hxtop : x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (m.isHermitian_stackGramW _ N ω)) :
    m.mleLogLik c N ω v ≤ m.mleLogLik c N ω (WithLp.ofLp x) := by
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
  have hxray : WithLp.ofLp x ⬝ᵥ A *ᵥ WithLp.ofLp x = lamMax A hAh := by
    rw [← hinner x, toOp_of_mem_topSpace hxtop, real_inner_smul_left,
      real_inner_self_eq_norm_sq, hx]
    ring
  have hvray : v ⬝ᵥ A *ᵥ v ≤ lamMax A hAh := by
    have hy := inner_toOp_self_le A hAh (WithLp.toLp 2 v)
    rw [hinner (WithLp.toLp 2 v)] at hy
    have hnorm : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin (d N)))‖ ^ 2 = 1 := by
      rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct]
      exact hv
    rw [hnorm, mul_one] at hy
    exact hy
  rw [m.mleLogLik_le_iff c hc N ω v (WithLp.ofLp x) hv hxdot, ← hA, hxray]
  exact hvray

end MultiTableModel

end StackedSVD
