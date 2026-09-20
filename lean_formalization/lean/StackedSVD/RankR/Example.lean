/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Defs
import StackedSVD.RankR.Unweighted

/-!
# The worked example of Section 7: `eq:psi_equation` and its two closed forms

Task T6 of `notes/RANK_R_PLAN.md`. Class (A): every result here is an identity or an
inequality between finitely many real numbers, or an identity between explicit `2 × 2`
matrices. No probability appears and no random matrix theory input is used.

The paper's example (`main_paper.tex:848`, `eq:psi_equation`) is

```
R_1 = (1, 0)ᵀ,  R_2 = (sin ψ, cos ψ)ᵀ,  θ_1 = θ_2 = θ ≥ 1,  c_1 = c_2 = 1,  ψ ∈ [0, π/2),
```

with `M = 2`, `r = 2` and `r_i = 1`, so `r̃ = M = r = 2`. Write `s := sin ψ`. The paper then
states the two limits

```
‖V̂_svdstackᵀ V‖_F² → β²(1+s)/(1+sβ²) + β²(1-s)/(1-sβ²)                (eq:perf_rankr_ex_svdstack)
‖V̂_stacksvdᵀ V‖_F² → h(θ²(1+s)) + h(θ²(1-s)),  h(λ) = (λ²-2c)/(λ(λ+1)) 1{λ > √(2c)}
                                                                     (eq:perf_rankr_ex_stacksvd)
```

and analyzes both in `app:unaligned_example_proofs` (`main_paper.tex:2198`).

## Content

1. `specInvTop_eq_inv`: at `r ≥ p` and `det A ≠ 0` the top-`r` partial inverse of
   `LinAlg/SpecProjPerturb.lean` is the ordinary inverse; `specInvTop_two_eq_inv` is the
   `2 × 2` case. `limitR_eq_trace_inv` is the consequence for the limit of
   `prop:general_rank_unweighted_svdstack`: at `M ≤ r` it is `tr(B_Rᵀ A_{β,R}⁻¹ B_R)`. That
   is the `r = M` case that task T4 also needs.
2. `Rex`, `abetaR_example`: the explicit `2 × 2` matrix
   `A_{β,R} = [[1, β₁β₂ sin ψ], [β₁β₂ sin ψ, 1]]`, its characteristic polynomial
   (`det_sub_smul_abetaR_example`, roots `1 ± β₁β₂ sin ψ`) and its eigenvectors `(1, ±1)`
   (`abetaR_example_mulVec_add`, `abetaR_example_mulVec_sub`).
3. `limitR_example`: `eq:perf_rankr_ex_svdstack`, `limitR` in closed form. `svdstackEx` is
   the right hand side; `hasDerivAt_svdstackEx` and `svdstackEx_strictAntiOn` are the
   appendix's derivative and monotonicity claims (`main_paper.tex:2205`), and
   `svdstackEx_zero`, `svdstackEx_one` are its two degenerate cases.
4. `Cex`: the core matrix `C = ∑_i R_i Θ_i² R_iᵀ` of `prop:stacksvd_subspace`, with the
   eigenpairs `θ²(1 ± s)` at `R_1 ± R_2` (`main_paper.tex:2216`).
5. `hPerf`, `stacksvdEx`: `eq:perf_rankr_ex_stacksvd` (`stacksvdEx_eq_paper`), the link
   `h(λ) = betaSq (√λ) (2c)` to `prop:single_table` (`hPerf_eq_betaSq`), the kinks `sStar`
   and `sDag` with their locations, the four piecewise closed forms, the derivatives `h'`
   and `h''`, the two monotonicity claims, the continuity `continuousOn_stacksvdEx`, and the
   slope jump `2θ²/(1 + √(2c))` (`slope_jump_at_kink`).
6. The endpoint `s = 0`: `svdstackEx_zero_s` gives `2β²`, `stacksvdEx_zero_s` gives
   `2 betaSq θ (2c)`, and `stacksvdEx_lt_svdstackEx_zero` is the paper's strict comparison
   (`main_paper.tex:866`), which needs `c < θ⁴` (below that threshold both sides are `0`).
7. `example_perfR_tendsto` and `example_perfR_tendsto_gaussian`: `eq:perf_rankr_ex_svdstack`
   as a limit theorem, `perfR →p svdstackEx β s`. This composes
   `prop:general_rank_unweighted_svdstack` (`RankR/Defs.lean`) with item 3. At `r = M = 2` the
   top-`r` eigengap is vacuous (`topGap_of_card_le`), so the two theorems carry no condition
   on `sin ψ` and none on `β`.

## What is not claimed

The stacksvd limit `eq:perf_rankr_ex_stacksvd` is **not** proved here. `stacksvdEx` is a
scalar definition and `stacksvdEx_eq_paper` an identity between two closed forms. The
convergence is `example_perfStackR_tendsto` of `RankR/Subspace.lean`, which holds modulo the
hypothesis structure `SubspaceLaw` (`notes/archive/prop_stacksvd_subspace.md`).

The svdstack limit is proved (item 7) with `perfR` of `RankR/Defs.lean` as the performance
metric. `perfR` is the right side of the paper's own first display (`main_paper.tex:1982`);
`RankR/Frobenius.lean` ties it to the left side `‖V̂_svdstackᵀ V‖_F²` through the definition
`vhatSvdstack` and the identity `frobSq_vhatSvdstack`.

Every closed form was checked with sympy (34 identities, all exact) and against a direct
numpy evaluation of `limitR` at 20 random `(θ, ψ)` points (worst error 4.4e-16, seed
20260830); see `notes/archive/agent_reports/rank_r_t6.md`.
-/

open Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD.RankR.Example

/-! ### 1. `specInvTop` at `r ≥ p` is the ordinary inverse

`prop:general_rank_unweighted_svdstack` carries the top-`r` eigengap of `A_{β,R}`, which is
vacuous when `r = r̃`. The paper's own example is that case (`M = r = r̃ = 2`), and there the
partial inverse `specInvTop` is the ordinary matrix inverse. -/

variable {p : ℕ}

/-- Completeness of the eigenbasis: `∑_i v_i v_iᵀ = I`. The columns of
`Matrix.IsHermitian.eigenvectorUnitary` are the eigenvectors, and the unitary group gives
`U Uᵀ = I`. -/
private theorem sum_vecMulVec_eigenvectorBasis {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) :
    ∑ i : Fin p, Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis i))
      (WithLp.ofLp (hA.eigenvectorBasis i)) = (1 : Matrix (Fin p) (Fin p) ℝ) := by
  have hU : (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 :=
    Unitary.coe_mul_star_self _
  ext j k
  have h := congrFun (congrFun hU j) k
  simp only [Matrix.mul_apply, Matrix.star_apply, Matrix.IsHermitian.eigenvectorUnitary_apply,
    star_trivial] at h
  simpa [Matrix.sum_apply, Matrix.vecMulVec_apply] using h

/-- When `r` is at least the dimension and `A` is invertible, `specInvTop A hA r = A⁻¹`.
At `r ≥ p` the set `topEigSet A hA r` holds every eigenvalue, so every coefficient of
`specInvTop` is `λ_i⁻¹` and the sum is the spectral form of `A⁻¹`. -/
theorem specInvTop_eq_inv {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (hpr : p ≤ r) (hdet : A.det ≠ 0) : specInvTop A hA r = A⁻¹ := by
  have heig : ∀ i : Fin p, hA.eigenvalues i ≠ 0 := by
    intro i h
    exact hdet (by
      rw [hA.det_eq_prod_eigenvalues]
      exact Finset.prod_eq_zero (Finset.mem_univ i) (by simpa using h))
  have hmem : ∀ i : Fin p, hA.eigenvalues i ∈ topEigSet A hA r := by
    intro i
    obtain ⟨k, hk⟩ : ∃ k : Fin (Fintype.card (Fin p)), hA.eigenvalues i = hA.eigenvalues₀ k :=
      ⟨_, rfl⟩
    exact ⟨k, lt_of_lt_of_le (by simpa using k.isLt) hpr, hk.symm⟩
  refine (Matrix.inv_eq_right_inv ?_).symm
  have key : A * specInvTop A hA r = ∑ i : Fin p,
      Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis i))
        (WithLp.ofLp (hA.eigenvectorBasis i)) := by
    rw [specInvTop_eq_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Matrix.mul_smul, Matrix.mul_vecMulVec, hA.mulVec_eigenvectorBasis,
      Matrix.smul_vecMulVec, smul_smul, invCoef_of_mem (hmem i),
      inv_mul_cancel₀ (heig i), one_smul]
  rw [key, sum_vecMulVec_eigenvectorBasis]

/-- The `2 × 2` case of `specInvTop_eq_inv`, which is the paper's worked example
(`M = r = r̃ = 2`). A symmetric `2 × 2` matrix with two nonzero eigenvalues has
`det ≠ 0`, so the hypothesis `hdet` is the paper's "two positive distinct eigenvalues"
weakened to what the proof uses. -/
theorem specInvTop_two_eq_inv {A : Matrix (Fin 2) (Fin 2) ℝ} (hA : A.IsHermitian)
    (hdet : A.det ≠ 0) : specInvTop A hA 2 = A⁻¹ :=
  specInvTop_eq_inv hA le_rfl hdet

/-- The limit of `prop:general_rank_unweighted_svdstack` at `r = r̃`, that is `M ≤ r`:
`limitR β R = tr(B_Rᵀ A_{β,R}⁻¹ B_R)`. -/
theorem limitR_eq_trace_inv {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (hMr : M ≤ r) (hdet : (AbetaR β R).det ≠ 0) :
    limitR β R = Matrix.trace ((BR β R)ᵀ * (AbetaR β R)⁻¹ * BR β R) := by
  rw [limitR, specInvTop_eq_inv (isHermitian_AbetaR β R) hMr hdet]

/-! ### 2. The example's alignment vectors and the matrix `A_{β,R}` -/

/-- `eq:psi_equation` (`main_paper.tex:848`): `R_1 = (1, 0)ᵀ` and `R_2 = (sin ψ, cos ψ)ᵀ`.
Both are unit vectors of `ℝ²` for every `ψ`, and `⟪R_1, R_2⟫ = sin ψ =: s`. -/
noncomputable def Rex (ψ : ℝ) : Fin 2 → EuclideanSpace ℝ (Fin 2) :=
  ![WithLp.toLp 2 ![1, 0], WithLp.toLp 2 ![Real.sin ψ, Real.cos ψ]]

theorem inner_Rex (ψ : ℝ) (i j : Fin 2) :
    ⟪Rex ψ i, Rex ψ j⟫_ℝ = if i = j then 1 else Real.sin ψ := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  rw [real_inner_eq_dotProduct]
  fin_cases i <;> fin_cases j <;>
    simp only [Rex, dotProduct, Fin.sum_univ_two] <;> norm_num
  nlinarith [hpy]

/-- Every `R_i` of the example is a unit vector, which is the field `hR` of
`UnalignedModel`. -/
theorem norm_Rex (ψ : ℝ) (i : Fin 2) : ‖Rex ψ i‖ = 1 := by
  have h := real_inner_self_eq_norm_sq (Rex ψ i)
  rw [inner_Rex ψ i i, if_pos rfl] at h
  nlinarith [norm_nonneg (Rex ψ i)]

/-- The explicit `2 × 2` matrix `A_{β,R}` of the example: diagonal `1`, off diagonal
`β₁ β₂ sin ψ`. -/
theorem abetaR_example (β : Fin 2 → ℝ) (ψ : ℝ) :
    AbetaR β (Rex ψ) =
      !![1, β 0 * β 1 * Real.sin ψ; β 0 * β 1 * Real.sin ψ, 1] := by
  ext i j
  rw [abetaR_apply, inner_Rex]
  fin_cases i <;> fin_cases j <;> norm_num [mul_comm] <;> ring

/-- The explicit `2 × 2` matrix `B_R = diag(β) R_stack` of the example. -/
theorem BR_example (β : Fin 2 → ℝ) (ψ : ℝ) :
    BR β (Rex ψ) = !![β 0, 0; β 1 * Real.sin ψ, β 1 * Real.cos ψ] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [BR, Rex]

/-- `B_R B_Rᵀ` of the example. This is the one place where `sin² + cos² = 1` enters. -/
theorem BBt_example (β : Fin 2 → ℝ) (ψ : ℝ) :
    BR β (Rex ψ) * (BR β (Rex ψ))ᵀ =
      !![β 0 ^ 2, β 0 * β 1 * Real.sin ψ; β 0 * β 1 * Real.sin ψ, β 1 ^ 2] := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  rw [BR_example]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two] <;> nlinarith [hpy]

theorem det_abetaR_example (β : Fin 2 → ℝ) (ψ : ℝ) :
    (AbetaR β (Rex ψ)).det = 1 - (β 0 * β 1 * Real.sin ψ) ^ 2 := by
  rw [abetaR_example, Matrix.det_fin_two_of]; ring

/-- The eigenvalues of `A_{β,R}` in the example are `1 ± β₁ β₂ sin ψ`, stated as the
factorization of the characteristic polynomial `det(t I - A_{β,R})`. -/
theorem det_sub_smul_abetaR_example (β : Fin 2 → ℝ) (ψ t : ℝ) :
    (t • (1 : Matrix (Fin 2) (Fin 2) ℝ) - AbetaR β (Rex ψ)).det =
      (t - (1 + β 0 * β 1 * Real.sin ψ)) * (t - (1 - β 0 * β 1 * Real.sin ψ)) := by
  rw [abetaR_example, Matrix.det_fin_two]
  simp
  ring

/-- The eigenvector of `A_{β,R}` at the eigenvalue `1 + β₁ β₂ sin ψ` is `(1, 1)`. -/
theorem abetaR_example_mulVec_add (β : Fin 2 → ℝ) (ψ : ℝ) :
    AbetaR β (Rex ψ) *ᵥ ![1, 1] = (1 + β 0 * β 1 * Real.sin ψ) • ![1, 1] := by
  rw [abetaR_example]
  funext i
  fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  ring

/-- The eigenvector of `A_{β,R}` at the eigenvalue `1 - β₁ β₂ sin ψ` is `(1, -1)`. -/
theorem abetaR_example_mulVec_sub (β : Fin 2 → ℝ) (ψ : ℝ) :
    AbetaR β (Rex ψ) *ᵥ ![1, -1] = (1 - β 0 * β 1 * Real.sin ψ) • ![1, -1] := by
  rw [abetaR_example]
  funext i
  fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;> ring

private theorem inv_two_sym (q : ℝ) (hq : (1 : ℝ) - q ^ 2 ≠ 0) :
    (!![1, q; q, 1] : Matrix (Fin 2) (Fin 2) ℝ)⁻¹ = (1 - q ^ 2)⁻¹ • !![1, -q; -q, 1] := by
  refine Matrix.inv_eq_right_inv ?_
  rw [Matrix.mul_smul]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp only [Matrix.cons_mul, Nat.succ_eq_add_one, Nat.reduceAdd, Matrix.vecMul_cons,
      Matrix.head_cons, one_smul, Matrix.tail_cons, Matrix.smul_cons, smul_eq_mul, mul_neg,
      mul_one, Matrix.smul_empty, Matrix.empty_vecMul, add_zero, Matrix.add_cons, neg_add_cancel,
      Matrix.empty_add_empty, add_neg_cancel, Matrix.empty_mul, Equiv.symm_apply_apply,
      Fin.zero_eta, Fin.isValue, Matrix.smul_apply, Matrix.of_apply, Matrix.cons_val',
      Matrix.cons_val_zero, Matrix.cons_val_fin_one, Matrix.one_apply_eq, Fin.mk_one,
      Matrix.cons_val_one, mul_zero, ne_eq, zero_ne_one, not_false_eq_true, Matrix.one_apply_ne,
      one_ne_zero] <;>
    (field_simp; ring)

private theorem trace_two (a d q k : ℝ) :
    Matrix.trace ((!![a, q; q, d] : Matrix (Fin 2) (Fin 2) ℝ) * (k • !![1, -q; -q, 1]))
      = k * (a + d - 2 * q ^ 2) := by
  rw [Matrix.trace_fin_two]
  simp [Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-! ### 3. `eq:perf_rankr_ex_svdstack` -/

/-- `limitR` of the example with two different `β_i`, in closed form. The paper only states
the common-`β` case; this is the same computation with `β_1 ≠ β_2` allowed. -/
theorem limitR_example_gen (β : Fin 2 → ℝ) (ψ : ℝ)
    (hq : (β 0 * β 1 * Real.sin ψ) ^ 2 < 1) :
    limitR β (Rex ψ) =
      (β 0 ^ 2 + β 1 ^ 2 - 2 * (β 0 * β 1 * Real.sin ψ) ^ 2) /
        (1 - (β 0 * β 1 * Real.sin ψ) ^ 2) := by
  have hne : (1 : ℝ) - (β 0 * β 1 * Real.sin ψ) ^ 2 ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; exact absurd h.symm (ne_of_lt hq)
  have hdet : (AbetaR β (Rex ψ)).det ≠ 0 := by rw [det_abetaR_example]; exact hne
  rw [limitR_eq_trace_inv β (Rex ψ) le_rfl hdet, abetaR_example, inv_two_sym _ hne,
    Matrix.trace_mul_comm, ← Matrix.mul_assoc, BBt_example, trace_two]
  field_simp

/-- The right hand side of `eq:perf_rankr_ex_svdstack` (`main_paper.tex:859`), as a function
of `β` and `s = sin ψ`. -/
noncomputable def svdstackEx (b s : ℝ) : ℝ :=
  b ^ 2 * (1 + s) / (1 + s * b ^ 2) + b ^ 2 * (1 - s) / (1 - s * b ^ 2)

/-- `eq:perf_rankr_ex_svdstack` (`main_paper.tex:859`): with `β_1 = β_2 = β`, the limit
`limitR` of `prop:general_rank_unweighted_svdstack` is
`β²(1+s)/(1+sβ²) + β²(1-s)/(1-sβ²)`, `s = sin ψ`. -/
theorem limitR_example (b ψ : ℝ) (hb : b ^ 2 < 1) :
    limitR (fun _ => b) (Rex ψ) = svdstackEx b (Real.sin ψ) := by
  have hs1 : Real.sin ψ ≤ 1 := Real.sin_le_one ψ
  have hs2 : -1 ≤ Real.sin ψ := Real.neg_one_le_sin ψ
  have hb0 : (0 : ℝ) ≤ b ^ 2 := sq_nonneg b
  have h1 : (0 : ℝ) < 1 + Real.sin ψ * b ^ 2 := by nlinarith
  have h2 : (0 : ℝ) < 1 - Real.sin ψ * b ^ 2 := by nlinarith
  have hq : (b * b * Real.sin ψ) ^ 2 < 1 := by nlinarith
  rw [limitR_example_gen (fun _ => b) ψ hq, svdstackEx]
  have h1' : (1 : ℝ) + b ^ 2 * Real.sin ψ ≠ 0 := by nlinarith
  have h2' : (1 : ℝ) - b ^ 2 * Real.sin ψ ≠ 0 := by nlinarith
  have hne : (1 : ℝ) - b ^ 4 * Real.sin ψ ^ 2 ≠ 0 := by nlinarith
  field_simp
  ring

/-- The single-fraction form of `eq:perf_rankr_ex_svdstack`:
`2β²(1 - s²β²)/(1 - s²β⁴)`. -/
theorem svdstackEx_eq_ratio (b s : ℝ) (h1 : 1 + s * b ^ 2 ≠ 0) (h2 : 1 - s * b ^ 2 ≠ 0) :
    svdstackEx b s = 2 * b ^ 2 * (1 - s ^ 2 * b ^ 2) / (1 - s ^ 2 * b ^ 4) := by
  have h3 : (1 : ℝ) - s ^ 2 * b ^ 4 ≠ 0 := by
    have h : (1 : ℝ) - s ^ 2 * b ^ 4 = (1 + s * b ^ 2) * (1 - s * b ^ 2) := by ring
    rw [h]; exact mul_ne_zero h1 h2
  have h1' : (1 : ℝ) + b ^ 2 * s ≠ 0 := by rw [mul_comm]; exact h1
  have h2' : (1 : ℝ) - b ^ 2 * s ≠ 0 := by rw [mul_comm]; exact h2
  have h3' : (1 : ℝ) - b ^ 4 * s ^ 2 ≠ 0 := by
    intro h; exact h3 (by nlinarith [h])
  rw [svdstackEx]
  field_simp
  ring

/-- Degenerate case `β = 0` of `main_paper.tex:2209`: the svdstack limit is identically `0`. -/
theorem svdstackEx_zero (s : ℝ) : svdstackEx 0 s = 0 := by simp [svdstackEx]

/-- Degenerate case `β = 1` of `main_paper.tex:2209`: the svdstack limit is identically `2`
on `s ∈ (-1, 1)`. -/
theorem svdstackEx_one {s : ℝ} (h1 : s ≠ -1) (h2 : s ≠ 1) : svdstackEx 1 s = 2 := by
  have ha : (1 : ℝ) + s ≠ 0 := fun h => h1 (by linarith)
  have hb : (1 : ℝ) - s ≠ 0 := fun h => h2 (by linarith)
  rw [svdstackEx]
  field_simp
  ring

/-- `main_paper.tex:2205`: the derivative of `eq:perf_rankr_ex_svdstack` in `s` is
`β²(1-β²)[(1+sβ²)^{-2} - (1-sβ²)^{-2}]`. -/
theorem hasDerivAt_svdstackEx (b : ℝ) {s : ℝ} (h1 : 1 + s * b ^ 2 ≠ 0)
    (h2 : 1 - s * b ^ 2 ≠ 0) :
    HasDerivAt (svdstackEx b)
      (b ^ 2 * (1 - b ^ 2) * (1 / (1 + s * b ^ 2) ^ 2 - 1 / (1 - s * b ^ 2) ^ 2)) s := by
  have hid : HasDerivAt (fun t : ℝ => t) 1 s := hasDerivAt_id s
  have hn1 : HasDerivAt (fun t : ℝ => b ^ 2 * (1 + t)) (b ^ 2) s := by
    simpa using (hid.const_add (1 : ℝ)).const_mul (b ^ 2)
  have hd1 : HasDerivAt (fun t : ℝ => 1 + t * b ^ 2) (b ^ 2) s := by
    simpa using (hid.mul_const (b ^ 2)).const_add (1 : ℝ)
  have hn2 : HasDerivAt (fun t : ℝ => b ^ 2 * (1 - t)) (-b ^ 2) s := by
    simpa using (hid.const_sub (1 : ℝ)).const_mul (b ^ 2)
  have hd2 : HasDerivAt (fun t : ℝ => 1 - t * b ^ 2) (-b ^ 2) s := by
    simpa using (hid.mul_const (b ^ 2)).const_sub (1 : ℝ)
  have hfun : svdstackEx b = fun t : ℝ =>
      b ^ 2 * (1 + t) / (1 + t * b ^ 2) + b ^ 2 * (1 - t) / (1 - t * b ^ 2) := rfl
  have h1' : ((1 : ℝ) + s * b ^ 2) ^ 2 ≠ 0 := pow_ne_zero 2 h1
  have h2' : ((1 : ℝ) - s * b ^ 2) ^ 2 ≠ 0 := pow_ne_zero 2 h2
  rw [hfun]
  refine ((hn1.div hd1 h1).add (hn2.div hd2 h2)).congr_deriv ?_
  field_simp
  ring

/-- `main_paper.tex:2205`: for `β ∈ (0, 1)` the svdstack limit of the example is strictly
decreasing in `s` on `[0, 1)`. -/
theorem svdstackEx_strictAntiOn {b : ℝ} (hb0 : 0 < b) (hb1 : b < 1) :
    StrictAntiOn (svdstackEx b) (Set.Ico (0 : ℝ) 1) := by
  have hbsq : (0 : ℝ) < b ^ 2 := by positivity
  have hbsq1 : b ^ 2 < 1 := by nlinarith
  have hden : ∀ t : ℝ, 0 ≤ t → t < 1 → (1 + t * b ^ 2 ≠ 0) ∧ (1 - t * b ^ 2 ≠ 0) := by
    intro t ht0 ht1
    constructor <;> [skip; skip] <;> intro h <;> nlinarith
  refine strictAntiOn_of_deriv_neg (convex_Ico 0 1) ?_ ?_
  · intro t ht
    obtain ⟨ht0, ht1⟩ := ht
    obtain ⟨ha, hb'⟩ := hden t ht0 ht1
    exact ((hasDerivAt_svdstackEx b ha hb').continuousAt).continuousWithinAt
  · intro t ht
    rw [interior_Ico] at ht
    obtain ⟨ht0, ht1⟩ := ht
    obtain ⟨ha, hb'⟩ := hden t (le_of_lt ht0) ht1
    rw [(hasDerivAt_svdstackEx b ha hb').deriv]
    have hpos : (0 : ℝ) < 1 - t * b ^ 2 := by nlinarith
    have hlt : (1 - t * b ^ 2) ^ 2 < (1 + t * b ^ 2) ^ 2 := by nlinarith
    have hd : 1 / (1 + t * b ^ 2) ^ 2 - 1 / (1 - t * b ^ 2) ^ 2 < 0 := by
      have := one_div_lt_one_div_of_lt (by positivity : (0:ℝ) < (1 - t * b ^ 2) ^ 2) hlt
      linarith
    exact mul_neg_of_pos_of_neg (by nlinarith) hd

/-! ### 4. The core matrix `C` of `prop:stacksvd_subspace`

`C = ∑_i R_i Θ_i² R_iᵀ` (`main_paper.tex:826`). In the example both `Θ_i` are the scalar `θ`,
so `C = θ² (R_1 R_1ᵀ + R_2 R_2ᵀ)`. The paper computes its spectrum at
`main_paper.tex:2216`. -/

/-- The core matrix of the example. -/
noncomputable def Cex (θ ψ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  θ ^ 2 • ∑ i : Fin 2,
    Matrix.vecMulVec (WithLp.ofLp (Rex ψ i)) (WithLp.ofLp (Rex ψ i))

/-- `main_paper.tex:2216`: `C = θ² [[1 + s², s cos ψ], [s cos ψ, cos²ψ]]`, with
`cos ψ = √(1 - s²)` on `[0, π/2)`. -/
theorem cex_eq (θ ψ : ℝ) :
    Cex θ ψ = θ ^ 2 • !![1 + Real.sin ψ ^ 2, Real.sin ψ * Real.cos ψ;
      Real.sin ψ * Real.cos ψ, Real.cos ψ ^ 2] := by
  rw [Cex]
  congr 1
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Rex, Matrix.sum_apply, Fin.sum_univ_two, Matrix.vecMulVec_apply] <;> ring

/-- The eigenvector of `C` at the eigenvalue `θ²(1 + s)` is `R_1 + R_2 = (1 + s, cos ψ)`. -/
theorem cex_mulVec_add (θ ψ : ℝ) :
    Cex θ ψ *ᵥ ![1 + Real.sin ψ, Real.cos ψ]
      = (θ ^ 2 * (1 + Real.sin ψ)) • ![1 + Real.sin ψ, Real.cos ψ] := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  rw [cex_eq]
  funext i
  fin_cases i
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    linear_combination (θ ^ 2 * Real.sin ψ) * hpy
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    linear_combination (θ ^ 2 * Real.cos ψ) * hpy

/-- The eigenvector of `C` at the eigenvalue `θ²(1 - s)` is `R_1 - R_2 = (1 - s, -cos ψ)`. -/
theorem cex_mulVec_sub (θ ψ : ℝ) :
    Cex θ ψ *ᵥ ![1 - Real.sin ψ, -Real.cos ψ]
      = (θ ^ 2 * (1 - Real.sin ψ)) • ![1 - Real.sin ψ, -Real.cos ψ] := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  rw [cex_eq]
  funext i
  fin_cases i
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    linear_combination (-(θ ^ 2 * Real.sin ψ)) * hpy
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    linear_combination (-(θ ^ 2 * Real.cos ψ)) * hpy

/-- The eigenvalues of the core matrix are `θ²(1 ± s)`, as the factorization of
`det(t I - C)`. -/
theorem det_sub_smul_cex (θ ψ t : ℝ) :
    (t • (1 : Matrix (Fin 2) (Fin 2) ℝ) - Cex θ ψ).det =
      (t - θ ^ 2 * (1 + Real.sin ψ)) * (t - θ ^ 2 * (1 - Real.sin ψ)) := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  rw [cex_eq, Matrix.det_fin_two]
  simp
  linear_combination (θ ^ 4 - t * θ ^ 2) * hpy

/-! ### 5. `eq:perf_rankr_ex_stacksvd`: the function `h`, the limit `F`, and the kinks

`main_paper.tex:2219` defines `h(λ) = (λ² - 2c)/(λ(λ+1)) 1{λ > √(2c)}` and writes the
stacksvd limit of the example as `F(s) = h(λ_+(s)) + h(λ_-(s))` with `λ_±(s) = θ²(1 ± s)`,
the two eigenvalues of the core matrix. -/

/-- The smooth branch of `h`: `(λ² - 2c)/(λ(λ+1))`. -/
noncomputable def hSmooth (c lam : ℝ) : ℝ := (lam ^ 2 - 2 * c) / (lam * (lam + 1))

/-- `h(λ) = (λ² - 2c)/(λ(λ+1)) 1{λ > √(2c)}` (`main_paper.tex:2219`). -/
noncomputable def hPerf (c lam : ℝ) : ℝ :=
  if Real.sqrt (2 * c) < lam then hSmooth c lam else 0

private theorem two_mul_le_sq_of_sqrt_le {c lam : ℝ} (h : Real.sqrt (2 * c) ≤ lam) :
    2 * c ≤ lam ^ 2 := by
  have h0 : 0 ≤ Real.sqrt (2 * c) := Real.sqrt_nonneg _
  have hs : Real.sqrt (2 * c) ^ 2 = max (2 * c) 0 := Real.sq_sqrt'
  nlinarith [le_max_left (2 * c) 0]

private theorem hSmooth_den_pos {lam : ℝ} (hlam : 0 < lam) : 0 < lam * (lam + 1) :=
  mul_pos hlam (by linarith)

private theorem hSmooth_nonneg {c lam : ℝ} (hlam : 0 < lam) (h : 2 * c ≤ lam ^ 2) :
    0 ≤ hSmooth c lam := div_nonneg (by linarith) (hSmooth_den_pos hlam).le

private theorem hSmooth_nonpos {c lam : ℝ} (hlam : 0 < lam) (h : lam ^ 2 ≤ 2 * c) :
    hSmooth c lam ≤ 0 := by
  rw [hSmooth, div_le_iff₀ (hSmooth_den_pos hlam)]
  simpa using by linarith

/-- `h` written with the indicator in the squared form of `eq:perf_rankr_ex_stacksvd`. -/
theorem hPerf_apply {c lam : ℝ} (hlam : 0 < lam) :
    hPerf c lam = if 2 * c < lam ^ 2 then (lam ^ 2 - 2 * c) / (lam ^ 2 + lam) else 0 := by
  have hiff : Real.sqrt (2 * c) < lam ↔ 2 * c < lam ^ 2 := Real.sqrt_lt' hlam
  rw [hPerf]
  by_cases h : 2 * c < lam ^ 2
  · rw [if_pos (hiff.mpr h), if_pos h, hSmooth]
    congr 1
    ring
  · rw [if_neg (fun hh => h (hiff.mp hh)), if_neg h]

/-- `h(λ) = betaSq (√λ) (2c)`: each term of `eq:perf_rankr_ex_stacksvd` is the rank-one
performance of `prop:single_table` at spike `√λ` and aspect ratio `‖c‖₁ = 2c`
(`main_paper.tex:838`). -/
theorem hPerf_eq_betaSq {c lam : ℝ} (hlam : 0 < lam) :
    hPerf c lam = betaSq (Real.sqrt lam) (2 * c) := by
  have h2 : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hlam.le
  have h4 : Real.sqrt lam ^ 4 = lam ^ 2 := by
    rw [show (4 : ℕ) = 2 * 2 from rfl, pow_mul, h2]
  simp only [hPerf_apply hlam, betaSq, h4, h2, gt_iff_lt]

/-- `h(√(2c)) = 0` on the smooth branch as well: the performance vanishes continuously at the
detection threshold (`main_paper.tex:2268`). -/
theorem hSmooth_sqrt {c : ℝ} (hc : 0 < c) : hSmooth c (Real.sqrt (2 * c)) = 0 := by
  rw [hSmooth, Real.sq_sqrt (by linarith : (0 : ℝ) ≤ 2 * c), sub_self, zero_div]

/-- The indicator in `h` is the positive part: for `λ > 0`, `h(λ) = max(hSmooth c λ, 0)`.
This is why `F` is continuous even where an indicator switches. -/
theorem hPerf_eq_max {c lam : ℝ} (hlam : 0 < lam) :
    hPerf c lam = max (hSmooth c lam) 0 := by
  have hiff : Real.sqrt (2 * c) < lam ↔ 2 * c < lam ^ 2 := Real.sqrt_lt' hlam
  rw [hPerf]
  by_cases h : Real.sqrt (2 * c) < lam
  · rw [if_pos h, max_eq_left (hSmooth_nonneg hlam (hiff.mp h).le)]
  · rw [if_neg h,
      max_eq_right (hSmooth_nonpos hlam (not_lt.mp fun hh => h (hiff.mpr hh)))]

/-- `F(s) = h(θ²(1 + s)) + h(θ²(1 - s))`, the right hand side of
`eq:perf_rankr_ex_stacksvd` (`main_paper.tex:860`). -/
noncomputable def stacksvdEx (θ c s : ℝ) : ℝ :=
  hPerf c (θ ^ 2 * (1 + s)) + hPerf c (θ ^ 2 * (1 - s))

/-- `eq:perf_rankr_ex_stacksvd` verbatim: `F` is the paper's two-term expression with the two
indicators `θ⁴(1 ± s)² > 2c`. -/
theorem stacksvdEx_eq_paper {θ c s : ℝ} (hp : 0 < θ ^ 2 * (1 + s))
    (hm : 0 < θ ^ 2 * (1 - s)) :
    stacksvdEx θ c s =
      (if 2 * c < θ ^ 4 * (1 + s) ^ 2 then
        (θ ^ 4 * (1 + s) ^ 2 - 2 * c) / (θ ^ 4 * (1 + s) ^ 2 + θ ^ 2 * (1 + s)) else 0) +
      (if 2 * c < θ ^ 4 * (1 - s) ^ 2 then
        (θ ^ 4 * (1 - s) ^ 2 - 2 * c) / (θ ^ 4 * (1 - s) ^ 2 + θ ^ 2 * (1 - s)) else 0) := by
  have e1 : (θ ^ 2 * (1 + s)) ^ 2 = θ ^ 4 * (1 + s) ^ 2 := by ring
  have e2 : (θ ^ 2 * (1 - s)) ^ 2 = θ ^ 4 * (1 - s) ^ 2 := by ring
  simp only [stacksvdEx, hPerf_apply hp, hPerf_apply hm, e1, e2]

/-! #### The kinks -/

/-- `s⋆ = 1 - √(2c)/θ²`, the kink of the strong signal regime (`eq:kink_location`,
`main_paper.tex:2242`). -/
noncomputable def sStar (θ c : ℝ) : ℝ := 1 - Real.sqrt (2 * c) / θ ^ 2

/-- `s_† = √(2c)/θ² - 1`, the kink of the weak signal regime (`main_paper.tex:2281`). -/
noncomputable def sDag (θ c : ℝ) : ℝ := Real.sqrt (2 * c) / θ ^ 2 - 1

/-- The kink `s⋆` is where component `2` leaves detection: `λ_-(s⋆) = √(2c)`. -/
theorem lam_neg_sStar {θ : ℝ} (hθ : θ ≠ 0) (c : ℝ) :
    θ ^ 2 * (1 - sStar θ c) = Real.sqrt (2 * c) := by
  have h : (θ : ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 hθ
  have e : (1 : ℝ) - sStar θ c = Real.sqrt (2 * c) / θ ^ 2 := by rw [sStar]; ring
  rw [e]
  field_simp

/-- The kink `s_†` is where component `1` enters detection: `λ_+(s_†) = √(2c)`. -/
theorem lam_pos_sDag {θ : ℝ} (hθ : θ ≠ 0) (c : ℝ) :
    θ ^ 2 * (1 + sDag θ c) = Real.sqrt (2 * c) := by
  have h : (θ : ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 hθ
  have e : (1 : ℝ) + sDag θ c = Real.sqrt (2 * c) / θ ^ 2 := by rw [sDag]; ring
  rw [e]
  field_simp

private theorem sqrt_two_mul_lt_sq {θ c : ℝ} (hθ : 0 < θ) (h : 2 * c < θ ^ 4) :
    Real.sqrt (2 * c) < θ ^ 2 :=
  (Real.sqrt_lt' (pow_pos hθ 2)).mpr (by nlinarith)

private theorem sq_lt_sqrt_two_mul {θ c : ℝ} (h : θ ^ 4 < 2 * c) :
    θ ^ 2 < Real.sqrt (2 * c) :=
  (Real.lt_sqrt (sq_nonneg θ)).mpr (by nlinarith)

/-- In the strong signal regime `θ⁴ > 2c` the kink `s⋆` is positive
(`main_paper.tex:2240`). -/
theorem sStar_pos {θ c : ℝ} (hθ : 0 < θ) (hstrong : 2 * c < θ ^ 4) : 0 < sStar θ c := by
  rw [sStar, sub_pos, div_lt_one (pow_pos hθ 2)]
  exact sqrt_two_mul_lt_sq hθ hstrong

/-- The kink `s⋆` is below `1` whenever `c > 0`. -/
theorem sStar_lt_one {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c) : sStar θ c < 1 := by
  have h1 : 0 < Real.sqrt (2 * c) := Real.sqrt_pos.mpr (by linarith)
  rw [sStar, sub_lt_self_iff]
  exact div_pos h1 (pow_pos hθ 2)

/-- In the weak signal regime `θ⁴ < 2c` the kink `s_†` is positive
(`main_paper.tex:2281`). -/
theorem sDag_pos {θ c : ℝ} (hθ : 0 < θ) (hweak : θ ^ 4 < 2 * c) : 0 < sDag θ c := by
  rw [sDag, sub_pos, lt_div_iff₀ (pow_pos hθ 2), one_mul]
  exact sq_lt_sqrt_two_mul hweak

/-- `s_† < 1` when `θ⁴ > c/2` (`main_paper.tex:2284`). -/
theorem sDag_lt_one {θ c : ℝ} (hθ : 0 < θ) (h : c / 2 < θ ^ 4) : sDag θ c < 1 := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  rw [sDag, sub_lt_iff_lt_add, div_lt_iff₀ hθ2]
  exact (Real.sqrt_lt' (by linarith)).mpr (by nlinarith)

/-- `θ⁴ ≤ c/2` puts the kink at or above `1`, and then the stacksvd limit is `0` on all of
`[0, 1)` (`main_paper.tex:2284`). -/
theorem one_le_sDag {θ c : ℝ} (hθ : 0 < θ) (h : θ ^ 4 ≤ c / 2) : 1 ≤ sDag θ c := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  rw [sDag, le_sub_iff_add_le, le_div_iff₀ hθ2]
  exact (Real.le_sqrt' (by linarith)).mpr (by nlinarith)

/-! #### The piecewise closed form -/

/-- Right of the kink in the strong signal regime `θ⁴ > 2c`: component `2` is at or below its
detection threshold, so only `λ_+` contributes (`main_paper.tex:2251`). -/
theorem stacksvdEx_of_sStar_le {θ c s : ℝ} (hθ : 0 < θ) (hstrong : 2 * c < θ ^ 4)
    (hs : sStar θ c ≤ s) : stacksvdEx θ c s = hSmooth c (θ ^ 2 * (1 + s)) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hstar : 0 < sStar θ c := sStar_pos hθ hstrong
  have hs0 : 0 < s := lt_of_lt_of_le hstar hs
  have hlt : Real.sqrt (2 * c) < θ ^ 2 := sqrt_two_mul_lt_sq hθ hstrong
  have hplus : Real.sqrt (2 * c) < θ ^ 2 * (1 + s) := by nlinarith
  have hkink : θ ^ 2 * (1 - sStar θ c) = Real.sqrt (2 * c) := lam_neg_sStar (ne_of_gt hθ) c
  have hminus : ¬ Real.sqrt (2 * c) < θ ^ 2 * (1 - s) := by
    refine not_lt.mpr ?_
    rw [← hkink]
    nlinarith [mul_nonneg hθ2.le (sub_nonneg.mpr hs)]
  simp only [stacksvdEx, hPerf, if_pos hplus, if_neg hminus, add_zero]

/-- Left of the kink in the strong signal regime: both components are detected. At `s = s⋆`
the second term is `h(√(2c)) = 0 = hSmooth c (√(2c))`, so the two branches agree there and
the curve is continuous (`main_paper.tex:2268`). -/
theorem stacksvdEx_of_le_sStar {θ c s : ℝ} (hθ : 0 < θ) (hc : 0 < c)
    (hstrong : 2 * c < θ ^ 4) (hs0 : 0 ≤ s) (hs : s ≤ sStar θ c) :
    stacksvdEx θ c s = hSmooth c (θ ^ 2 * (1 + s)) + hSmooth c (θ ^ 2 * (1 - s)) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hs1 : s < 1 := lt_of_le_of_lt hs (sStar_lt_one hθ hc)
  have hp : 0 < θ ^ 2 * (1 + s) := by nlinarith
  have hm : 0 < θ ^ 2 * (1 - s) := by nlinarith
  have hkink : θ ^ 2 * (1 - sStar θ c) = Real.sqrt (2 * c) := lam_neg_sStar (ne_of_gt hθ) c
  have hgem : Real.sqrt (2 * c) ≤ θ ^ 2 * (1 - s) := by
    rw [← hkink]
    nlinarith [mul_nonneg hθ2.le (sub_nonneg.mpr hs)]
  have hlt : Real.sqrt (2 * c) < θ ^ 2 := sqrt_two_mul_lt_sq hθ hstrong
  have hgep : Real.sqrt (2 * c) ≤ θ ^ 2 * (1 + s) := by nlinarith
  rw [stacksvdEx, hPerf_eq_max hp, hPerf_eq_max hm,
    max_eq_left (hSmooth_nonneg hp (two_mul_le_sq_of_sqrt_le hgep)),
    max_eq_left (hSmooth_nonneg hm (two_mul_le_sq_of_sqrt_le hgem))]

/-- Below the weak-regime kink `s_†` nothing is detected and `F ≡ 0`
(`main_paper.tex:2281`). -/
theorem stacksvdEx_of_le_sDag {θ c s : ℝ} (hθ : 0 < θ) (hs0 : 0 ≤ s) (hs : s ≤ sDag θ c) :
    stacksvdEx θ c s = 0 := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hkink : θ ^ 2 * (1 + sDag θ c) = Real.sqrt (2 * c) := lam_pos_sDag (ne_of_gt hθ) c
  have hple : θ ^ 2 * (1 + s) ≤ Real.sqrt (2 * c) := by
    rw [← hkink]
    nlinarith [mul_nonneg hθ2.le (sub_nonneg.mpr hs)]
  have hp : ¬ Real.sqrt (2 * c) < θ ^ 2 * (1 + s) := not_lt.mpr hple
  have hm : ¬ Real.sqrt (2 * c) < θ ^ 2 * (1 - s) := by
    refine not_lt.mpr (le_trans ?_ hple)
    nlinarith
  simp only [stacksvdEx, hPerf, if_neg hp, if_neg hm, add_zero]

/-- Above the weak-regime kink only component `1` is detected (`main_paper.tex:2282`). -/
theorem stacksvdEx_of_sDag_le {θ c s : ℝ} (hθ : 0 < θ) (hweak : θ ^ 4 < 2 * c)
    (hs0 : 0 ≤ s) (hs : sDag θ c ≤ s) (hs1 : s < 1) :
    stacksvdEx θ c s = hSmooth c (θ ^ 2 * (1 + s)) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hlt : θ ^ 2 < Real.sqrt (2 * c) := sq_lt_sqrt_two_mul hweak
  have hp0 : 0 < θ ^ 2 * (1 + s) := by nlinarith
  have hsub : (0 : ℝ) < 1 - s := sub_pos.mpr hs1
  have hmneg : ¬ Real.sqrt (2 * c) < θ ^ 2 * (1 - s) := not_lt.mpr (by nlinarith [hsub])
  have hkink : θ ^ 2 * (1 + sDag θ c) = Real.sqrt (2 * c) := lam_pos_sDag (ne_of_gt hθ) c
  have hgep : Real.sqrt (2 * c) ≤ θ ^ 2 * (1 + s) := by
    rw [← hkink]
    nlinarith [mul_nonneg hθ2.le (sub_nonneg.mpr hs)]
  rw [stacksvdEx, hPerf_eq_max hp0,
    max_eq_left (hSmooth_nonneg hp0 (two_mul_le_sq_of_sqrt_le hgep)), hPerf, if_neg hmneg,
    add_zero]

/-- In the lowest signal regime `θ⁴ ≤ c/2` no method has nonzero performance
(`main_paper.tex:2285`). -/
theorem stacksvdEx_eq_zero_of_very_weak {θ c s : ℝ} (hθ : 0 < θ) (h : θ ^ 4 ≤ c / 2)
    (hs0 : 0 ≤ s) (hs1 : s < 1) : stacksvdEx θ c s = 0 :=
  stacksvdEx_of_le_sDag hθ hs0 (le_trans hs1.le (one_le_sDag hθ h))

/-! #### The calculus of `h` -/

/-- `h'(λ) = (λ² + 4cλ + 2c)/(λ²(λ+1)²)` (`main_paper.tex:2246`). -/
noncomputable def hSmoothD (c lam : ℝ) : ℝ :=
  (lam ^ 2 + 4 * c * lam + 2 * c) / (lam ^ 2 * (lam + 1) ^ 2)

theorem hasDerivAt_hSmooth (c : ℝ) {lam : ℝ} (hlam : 0 < lam) :
    HasDerivAt (hSmooth c) (hSmoothD c lam) lam := by
  have hfun : hSmooth c = fun x : ℝ => (x ^ 2 - 2 * c) / (x * (x + 1)) := rfl
  have hd : lam * (lam + 1) ≠ 0 := ne_of_gt (hSmooth_den_pos hlam)
  have hnum : HasDerivAt (fun x : ℝ => x ^ 2 - 2 * c) (2 * lam) lam := by
    simpa using ((hasDerivAt_id lam).pow 2).sub_const (2 * c)
  have hden : HasDerivAt (fun x : ℝ => x * (x + 1)) (2 * lam + 1) lam :=
    ((hasDerivAt_id lam).mul ((hasDerivAt_id lam).add_const (1 : ℝ))).congr_deriv (by simp; ring)
  rw [hfun]
  refine (hnum.div hden hd).congr_deriv ?_
  rw [hSmoothD]
  field_simp
  ring

/-- `h''(λ) = -2(λ³ + 6cλ² + 6cλ + 2c)/(λ³(λ+1)³)` (`main_paper.tex:2247`). -/
theorem hasDerivAt_hSmoothD (c : ℝ) {lam : ℝ} (hlam : 0 < lam) :
    HasDerivAt (hSmoothD c)
      (-2 * (lam ^ 3 + 6 * c * lam ^ 2 + 6 * c * lam + 2 * c) / (lam ^ 3 * (lam + 1) ^ 3))
      lam := by
  have hfun : hSmoothD c =
      fun x : ℝ => (x ^ 2 + 4 * c * x + 2 * c) / (x ^ 2 * (x + 1) ^ 2) := rfl
  have hd : lam ^ 2 * (lam + 1) ^ 2 ≠ 0 :=
    ne_of_gt (mul_pos (pow_pos hlam 2) (pow_pos (by linarith) 2))
  have hnum : HasDerivAt (fun x : ℝ => x ^ 2 + 4 * c * x + 2 * c) (2 * lam + 4 * c) lam := by
    simpa using
      (((hasDerivAt_id lam).pow 2).add ((hasDerivAt_id lam).const_mul (4 * c))).add_const (2 * c)
  have hden : HasDerivAt (fun x : ℝ => x ^ 2 * (x + 1) ^ 2)
      (2 * lam * (lam + 1) ^ 2 + lam ^ 2 * (2 * (lam + 1))) lam :=
    (((hasDerivAt_id lam).pow 2).mul
      (((hasDerivAt_id lam).add_const (1 : ℝ)).pow 2)).congr_deriv
      (by simp only [id_eq, Pi.pow_apply]; push_cast; ring)
  rw [hfun]
  refine (hnum.div hden hd).congr_deriv ?_
  field_simp
  ring

/-- `h' > 0`: `h` is strictly increasing (`main_paper.tex:2246`). -/
theorem hSmoothD_pos {c lam : ℝ} (hc : 0 < c) (hlam : 0 < lam) : 0 < hSmoothD c lam :=
  div_pos (by nlinarith) (mul_pos (pow_pos hlam 2) (pow_pos (by linarith) 2))

/-- `h'' < 0`: `h` is strictly concave (`main_paper.tex:2247`). -/
theorem hSmoothDD_neg {c lam : ℝ} (hc : 0 < c) (hlam : 0 < lam) :
    -2 * (lam ^ 3 + 6 * c * lam ^ 2 + 6 * c * lam + 2 * c) / (lam ^ 3 * (lam + 1) ^ 3) < 0 :=
  div_neg_of_neg_of_pos (by nlinarith) (mul_pos (pow_pos hlam 3) (pow_pos (by linarith) 3))

/-- `h` is strictly increasing on `(0, ∞)`. -/
theorem hSmooth_strictMonoOn {c : ℝ} (hc : 0 < c) : StrictMonoOn (hSmooth c) (Set.Ioi 0) := by
  refine strictMonoOn_of_deriv_pos (convex_Ioi 0) ?_ ?_
  · exact fun x hx => ((hasDerivAt_hSmooth c hx).continuousAt).continuousWithinAt
  · intro x hx
    rw [interior_Ioi] at hx
    rw [(hasDerivAt_hSmooth c hx).deriv]
    exact hSmoothD_pos hc hx

/-- `h'` is strictly decreasing on `(0, ∞)`, the concavity that drives the left branch
(`main_paper.tex:2249`). -/
theorem hSmoothD_strictAntiOn {c : ℝ} (hc : 0 < c) :
    StrictAntiOn (hSmoothD c) (Set.Ioi 0) := by
  refine strictAntiOn_of_deriv_neg (convex_Ioi 0) ?_ ?_
  · exact fun x hx => ((hasDerivAt_hSmoothD c hx).continuousAt).continuousWithinAt
  · intro x hx
    rw [interior_Ioi] at hx
    rw [(hasDerivAt_hSmoothD c hx).deriv]
    exact hSmoothDD_neg hc hx

/-! #### The two branches of `F`, its monotonicity, its continuity and the slope jump -/

private theorem hasDerivAt_lamPos (θ s : ℝ) :
    HasDerivAt (fun t : ℝ => θ ^ 2 * (1 + t)) (θ ^ 2) s := by
  simpa using ((hasDerivAt_id s).const_add (1 : ℝ)).const_mul (θ ^ 2)

private theorem hasDerivAt_lamNeg (θ s : ℝ) :
    HasDerivAt (fun t : ℝ => θ ^ 2 * (1 - t)) (-θ ^ 2) s := by
  simpa using ((hasDerivAt_id s).const_sub (1 : ℝ)).const_mul (θ ^ 2)

/-- Slope of the right branch: `d/ds h(λ_+(s)) = θ² h'(λ_+(s))` (`main_paper.tex:2255`). -/
private theorem hasDerivAt_Fright (θ c : ℝ) {s : ℝ} (hp : 0 < θ ^ 2 * (1 + s)) :
    HasDerivAt (fun t : ℝ => hSmooth c (θ ^ 2 * (1 + t)))
      (θ ^ 2 * hSmoothD c (θ ^ 2 * (1 + s))) s := by
  have h : HasDerivAt (hSmooth c ∘ fun t : ℝ => θ ^ 2 * (1 + t))
      (hSmoothD c (θ ^ 2 * (1 + s)) * θ ^ 2) s :=
    (hasDerivAt_hSmooth c hp).comp s (hasDerivAt_lamPos θ s)
  rw [Function.comp_def] at h
  exact h.congr_deriv (by ring)

/-- Slope of the left branch: `F'(s) = θ²[h'(λ_+(s)) - h'(λ_-(s))]` (`main_paper.tex:2262`). -/
private theorem hasDerivAt_Fleft (θ c : ℝ) {s : ℝ} (hp : 0 < θ ^ 2 * (1 + s))
    (hm : 0 < θ ^ 2 * (1 - s)) :
    HasDerivAt (fun t : ℝ => hSmooth c (θ ^ 2 * (1 + t)) + hSmooth c (θ ^ 2 * (1 - t)))
      (θ ^ 2 * (hSmoothD c (θ ^ 2 * (1 + s)) - hSmoothD c (θ ^ 2 * (1 - s)))) s := by
  have h1 : HasDerivAt (hSmooth c ∘ fun t : ℝ => θ ^ 2 * (1 + t))
      (hSmoothD c (θ ^ 2 * (1 + s)) * θ ^ 2) s :=
    (hasDerivAt_hSmooth c hp).comp s (hasDerivAt_lamPos θ s)
  have h2 : HasDerivAt (hSmooth c ∘ fun t : ℝ => θ ^ 2 * (1 - t))
      (hSmoothD c (θ ^ 2 * (1 - s)) * -θ ^ 2) s :=
    (hasDerivAt_hSmooth c hm).comp s (hasDerivAt_lamNeg θ s)
  rw [Function.comp_def] at h1 h2
  exact (h1.add h2).congr_deriv (by ring)

/-- Right of the kink `F` has the slope `θ² h'(λ_+(s))` (`main_paper.tex:2255`). -/
theorem hasDerivWithinAt_stacksvdEx_right {θ c s : ℝ} (hθ : 0 < θ) (hstrong : 2 * c < θ ^ 4)
    (hs : sStar θ c ≤ s) :
    HasDerivWithinAt (stacksvdEx θ c) (θ ^ 2 * hSmoothD c (θ ^ 2 * (1 + s)))
      (Set.Ici (sStar θ c)) s := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hstar : 0 < sStar θ c := sStar_pos hθ hstrong
  have hp : ∀ t : ℝ, sStar θ c ≤ t → 0 < θ ^ 2 * (1 + t) := fun t ht => by nlinarith
  exact ((hasDerivAt_Fright θ c (hp s hs)).hasDerivWithinAt).congr
    (fun t ht => stacksvdEx_of_sStar_le hθ hstrong ht)
    (stacksvdEx_of_sStar_le hθ hstrong hs)

/-- Left of the kink `F` has the slope `θ²[h'(λ_+) - h'(λ_-)]` (`main_paper.tex:2262`). -/
theorem hasDerivWithinAt_stacksvdEx_left {θ c s : ℝ} (hθ : 0 < θ) (hc : 0 < c)
    (hstrong : 2 * c < θ ^ 4) (hs0 : 0 ≤ s) (hs : s ≤ sStar θ c) :
    HasDerivWithinAt (stacksvdEx θ c)
      (θ ^ 2 * (hSmoothD c (θ ^ 2 * (1 + s)) - hSmoothD c (θ ^ 2 * (1 - s))))
      (Set.Icc 0 (sStar θ c)) s := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hstar1 : sStar θ c < 1 := sStar_lt_one hθ hc
  have hpos : ∀ t : ℝ, 0 ≤ t → t ≤ sStar θ c → 0 < θ ^ 2 * (1 + t) ∧ 0 < θ ^ 2 * (1 - t) :=
    fun t h0 h1 => ⟨by nlinarith, by nlinarith⟩
  exact ((hasDerivAt_Fleft θ c (hpos s hs0 hs).1 (hpos s hs0 hs).2).hasDerivWithinAt).congr
    (fun t ht => stacksvdEx_of_le_sStar hθ hc hstrong ht.1 ht.2)
    (stacksvdEx_of_le_sStar hθ hc hstrong hs0 hs)

/-- At `s = 0` the two core eigenvalues coincide and the slope vanishes
(`main_paper.tex:2276`). -/
theorem hasDerivWithinAt_stacksvdEx_zero {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c)
    (hstrong : 2 * c < θ ^ 4) :
    HasDerivWithinAt (stacksvdEx θ c) 0 (Set.Icc 0 (sStar θ c)) 0 := by
  simpa using hasDerivWithinAt_stacksvdEx_left hθ hc hstrong le_rfl (sStar_pos hθ hstrong).le

/-- Right of the kink `F` is strictly increasing (`main_paper.tex:2257`). -/
theorem stacksvdEx_strictMonoOn {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c)
    (hstrong : 2 * c < θ ^ 4) :
    StrictMonoOn (stacksvdEx θ c) (Set.Ici (sStar θ c)) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hstar : 0 < sStar θ c := sStar_pos hθ hstrong
  intro x hx y hy hxy
  have hx' : sStar θ c ≤ x := hx
  have hy' : sStar θ c ≤ y := hy
  rw [stacksvdEx_of_sStar_le hθ hstrong hx', stacksvdEx_of_sStar_le hθ hstrong hy']
  exact hSmooth_strictMonoOn hc (Set.mem_Ioi.mpr (by nlinarith))
    (Set.mem_Ioi.mpr (by nlinarith)) (by nlinarith)

/-- Left of the kink `F` is strictly decreasing: `λ_+ > λ_-` and `h'` is strictly decreasing
by concavity (`main_paper.tex:2262`). -/
theorem stacksvdEx_strictAntiOn {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c)
    (hstrong : 2 * c < θ ^ 4) :
    StrictAntiOn (stacksvdEx θ c) (Set.Icc 0 (sStar θ c)) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hstar1 : sStar θ c < 1 := sStar_lt_one hθ hc
  have hpos : ∀ t : ℝ, 0 ≤ t → t ≤ sStar θ c → 0 < θ ^ 2 * (1 + t) ∧ 0 < θ ^ 2 * (1 - t) :=
    fun t h0 h1 => ⟨by nlinarith, by nlinarith⟩
  have hAnti : StrictAntiOn
      (fun t : ℝ => hSmooth c (θ ^ 2 * (1 + t)) + hSmooth c (θ ^ 2 * (1 - t)))
      (Set.Icc 0 (sStar θ c)) := by
    refine strictAntiOn_of_deriv_neg (convex_Icc 0 (sStar θ c)) ?_ ?_
    · intro t ht
      exact ((hasDerivAt_Fleft θ c (hpos t ht.1 ht.2).1
        (hpos t ht.1 ht.2).2).continuousAt).continuousWithinAt
    · intro t ht
      rw [interior_Icc] at ht
      have h0 : (0 : ℝ) ≤ t := le_of_lt ht.1
      have h1 : t ≤ sStar θ c := le_of_lt ht.2
      rw [(hasDerivAt_Fleft θ c (hpos t h0 h1).1 (hpos t h0 h1).2).deriv]
      have hlt : θ ^ 2 * (1 - t) < θ ^ 2 * (1 + t) := by nlinarith [ht.1]
      have hmono := hSmoothD_strictAntiOn hc (Set.mem_Ioi.mpr (hpos t h0 h1).2)
        (Set.mem_Ioi.mpr (hpos t h0 h1).1) hlt
      exact mul_neg_of_pos_of_neg hθ2 (by linarith)
  intro x hx y hy hxy
  rw [stacksvdEx_of_le_sStar hθ hc hstrong hx.1 hx.2,
    stacksvdEx_of_le_sStar hθ hc hstrong hy.1 hy.2]
  exact hAnti hx hy hxy

/-- `F` is continuous on `(-1, 1)`, kinks included: an indicator switches exactly where the
term it multiplies vanishes (`main_paper.tex:2268`). -/
theorem continuousOn_stacksvdEx {θ c : ℝ} (hθ : 0 < θ) :
    ContinuousOn (stacksvdEx θ c) (Set.Ioo (-1 : ℝ) 1) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have hp : ∀ t ∈ Set.Ioo (-1 : ℝ) 1, 0 < θ ^ 2 * (1 + t) := fun t ht => by nlinarith [ht.1]
  have hm : ∀ t ∈ Set.Ioo (-1 : ℝ) 1, 0 < θ ^ 2 * (1 - t) := fun t ht => by nlinarith [ht.2]
  have c1 : ContinuousOn (fun t : ℝ => hSmooth c (θ ^ 2 * (1 + t))) (Set.Ioo (-1 : ℝ) 1) :=
    fun t ht => ((hasDerivAt_Fright θ c (hp t ht)).continuousAt).continuousWithinAt
  have c2 : ContinuousOn (fun t : ℝ => hSmooth c (θ ^ 2 * (1 - t))) (Set.Ioo (-1 : ℝ) 1) := by
    intro t ht
    have h : HasDerivAt (hSmooth c ∘ fun u : ℝ => θ ^ 2 * (1 - u))
        (hSmoothD c (θ ^ 2 * (1 - t)) * -θ ^ 2) t :=
      (hasDerivAt_hSmooth c (hm t ht)).comp t (hasDerivAt_lamNeg θ t)
    rw [Function.comp_def] at h
    exact h.continuousAt.continuousWithinAt
  have hmax : ContinuousOn
      (fun t : ℝ => max (hSmooth c (θ ^ 2 * (1 + t))) 0 + max (hSmooth c (θ ^ 2 * (1 - t))) 0)
      (Set.Ioo (-1 : ℝ) 1) :=
    (c1.sup continuousOn_const).add (c2.sup continuousOn_const)
  refine hmax.congr fun t ht => ?_
  rw [stacksvdEx, hPerf_eq_max (hp t ht), hPerf_eq_max (hm t ht)]

/-- `h'(√(2c)) = 2/(1 + √(2c))`. -/
theorem hSmoothD_sqrt {c : ℝ} (hc : 0 < c) :
    hSmoothD c (Real.sqrt (2 * c)) = 2 / (1 + Real.sqrt (2 * c)) := by
  have h0 : 0 < Real.sqrt (2 * c) := Real.sqrt_pos.mpr (by linarith)
  have hsq : Real.sqrt (2 * c) ^ 2 = 2 * c := Real.sq_sqrt (by linarith)
  have h1 : Real.sqrt (2 * c) + 1 ≠ 0 := by linarith
  have h2 : (1 : ℝ) + Real.sqrt (2 * c) ≠ 0 := by linarith
  have h3 : (2 : ℝ) * c ≠ 0 := by linarith
  rw [hSmoothD, hsq]
  field_simp
  ring

/-- The slope jump at the kink (`main_paper.tex:2270`):
`F'(s⋆+) - F'(s⋆-) = θ² h'(√(2c)) = 2θ²/(1 + √(2c))`. -/
theorem slope_jump_at_kink {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c) :
    θ ^ 2 * hSmoothD c (θ ^ 2 * (1 + sStar θ c))
        - θ ^ 2 * (hSmoothD c (θ ^ 2 * (1 + sStar θ c))
          - hSmoothD c (θ ^ 2 * (1 - sStar θ c)))
      = 2 * θ ^ 2 / (1 + Real.sqrt (2 * c)) := by
  have hq : θ ^ 2 * (1 - sStar θ c) = Real.sqrt (2 * c) := lam_neg_sStar (ne_of_gt hθ) c
  have hcollapse : θ ^ 2 * hSmoothD c (θ ^ 2 * (1 + sStar θ c))
      - θ ^ 2 * (hSmoothD c (θ ^ 2 * (1 + sStar θ c))
        - hSmoothD c (θ ^ 2 * (1 - sStar θ c)))
      = θ ^ 2 * hSmoothD c (θ ^ 2 * (1 - sStar θ c)) := by ring
  rw [hcollapse, hq, hSmoothD_sqrt hc]
  ring

/-- The jump is strictly positive, so the slope really does jump upward at `s⋆`. -/
theorem slope_jump_pos {θ c : ℝ} (hθ : 0 < θ) :
    0 < 2 * θ ^ 2 / (1 + Real.sqrt (2 * c)) :=
  div_pos (by nlinarith [pow_pos hθ 2]) (by positivity)

/-! #### The endpoint `s = 0` -/

theorem svdstackEx_zero_s (b : ℝ) : svdstackEx b 0 = 2 * b ^ 2 := by
  rw [svdstackEx]
  norm_num
  ring

theorem stacksvdEx_zero_s {θ c : ℝ} (hθ : 0 < θ) :
    stacksvdEx θ c 0 = 2 * betaSq θ (2 * c) := by
  have hθ2 : (0 : ℝ) < θ ^ 2 := pow_pos hθ 2
  have h1 : θ ^ 2 * (1 + 0) = θ ^ 2 := by ring
  have h2 : θ ^ 2 * (1 - 0) = θ ^ 2 := by ring
  rw [stacksvdEx, h1, h2, hPerf_eq_betaSq hθ2, Real.sqrt_sq hθ.le]
  ring

private theorem betaSq_nonneg' (θ c : ℝ) : 0 ≤ betaSq θ c := by
  rw [betaSq]
  split_ifs with h
  · exact div_nonneg (by linarith) (by positivity)
  · exact le_rfl

private theorem beta_sq' (θ c : ℝ) : beta θ c ^ 2 = betaSq θ c :=
  Real.sq_sqrt (betaSq_nonneg' θ c)

/-- Doubling the aspect ratio strictly lowers the rank-one performance above the threshold.
This is the paper's `2β²(θ, 2c) < 2β²(θ, c)` at `s = 0` (`main_paper.tex:866`). -/
theorem betaSq_two_mul_lt {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c) (h : c < θ ^ 4) :
    betaSq θ (2 * c) < betaSq θ c := by
  have hden : (0 : ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith [pow_pos hθ 4, pow_pos hθ 2]
  simp only [betaSq, gt_iff_lt]
  rw [if_pos h]
  split_ifs with h2
  · have key : (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) - (θ ^ 4 - 2 * c) / (θ ^ 4 + θ ^ 2)
        = c / (θ ^ 4 + θ ^ 2) := by field_simp; ring
    linarith [div_pos hc hden, key]
  · exact div_pos (by linarith) hden

/-- `main_paper.tex:866`: at `s = 0` svdstack strictly dominates stacksvd. svdstack analyzes
each table on its own and attains `2β²(θ, c)`; stacksvd dilutes the signal with the other
table's rows and attains only `2β²(θ, 2c)`. -/
theorem stacksvdEx_lt_svdstackEx_zero {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c) (h : c < θ ^ 4) :
    stacksvdEx θ c 0 < svdstackEx (beta θ c) 0 := by
  rw [stacksvdEx_zero_s hθ, svdstackEx_zero_s, beta_sq']
  linarith [betaSq_two_mul_lt hθ hc h]

/-- The paper's own parametrization of `eq:perf_rankr_ex_svdstack`: `β = β(θ, c)` of
`prop:single_table`, which lies in `[0, 1)` for every `c > 0`. -/
theorem limitR_example_beta {θ c ψ : ℝ} (hc : 0 < c) :
    limitR (fun _ => beta θ c) (Rex ψ) = svdstackEx (beta θ c) (Real.sin ψ) := by
  obtain ⟨h0, h1⟩ := beta_mem_Ico (θ := θ) hc
  exact limitR_example _ _ (by nlinarith)

/-! ### 7. The example as a limit theorem

Sections 1 to 6 are identities between real numbers. This section joins the svdstack half to
`prop:general_rank_unweighted_svdstack` and gives `eq:perf_rankr_ex_svdstack`
(`main_paper.tex:859`) as a statement about the data: `perfR →p svdstackEx β s`.

At `M = 2` and `r = 2` the model has `r̃ = M = r`, so the top-`r` eigengap of `A_{β,R}` is
vacuous. It is the paper's own convention `λ_{r̃+1} := -∞`, which
`LinAlg/SpecProjPerturb.lean` proves as `topGap_of_card_le`. The two theorems below therefore
need no condition on `sin ψ` and no condition `β > 0`; the eigenvalues `1 ± β₁ β₂ sin ψ` of
`abetaR_example` are not used. A gap between those two eigenvalues is what a **rank-one**
reading of the same matrix would need, not this one.

The stacksvd half is `example_perfStackR_tendsto` of `RankR/Subspace.lean`, proved modulo
the hypothesis structure `SubspaceLaw`.
-/

section Convergence

open MeasureTheory

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ}

/-- `eq:perf_rankr_ex_svdstack` (`main_paper.tex:859`) as a limit theorem, Layer 1 form.
Take the paper's example `eq:psi_equation` (`main_paper.tex:848`): `M = 2`, `r = r_i = 2`,
`R = Rex ψ`, `θ_1 = θ_2 = θ` and `c_1 = c_2 = c`. Then the svdstack performance `perfR`
converges in probability to `β²(1+s)/(1+sβ²) + β²(1-s)/(1-sβ²)` with `β = β(θ, c)` and
`s = sin ψ`.

This is `prop_general_rank_unweighted_svdstack` at `r = M = 2`, where the eigengap is free,
followed by the closed form `limitR_example_beta`. The paper writes `c = 1`; the statement
holds for every `c > 0`, and it does not need `ψ ∈ [0, π/2)` or `θ ≥ 1`. -/
theorem example_perfR_tendsto (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ)
    (law : ∀ i, (m.tbl i).SingleTableLaw c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω) (svdstackEx (beta θ c) (Real.sin ψ)) := by
  have hβ : ∀ i : Fin 2, (fun _ : Fin 2 => beta θ c) i = beta (m.tbl i).θ c := by
    intro i
    rw [hθ i]
  have h := m.prop_general_rank_unweighted_svdstack (fun _ => c) (fun _ => beta θ c)
    (fun _ => hc) hβ le_rfl
    (topGap_of_card_le _ (by simp)) law hG
  rwa [hR, limitR_example_beta hc] at h

/-- `eq:perf_rankr_ex_svdstack` as a limit theorem, Layer 2 form: the same conclusion from the
proportional regime of each table and the joint Gaussian noise law, with no `SingleTableLaw`
hypothesis. The single-table black box is discharged inside
`prop_general_rank_unweighted_svdstack_gaussian` through `singleTableLaw_of_gaussian`. -/
theorem example_perfR_tendsto_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω) (svdstackEx (beta θ c) (Real.sin ψ)) := by
  have hβ : ∀ i : Fin 2, (fun _ : Fin 2 => beta θ c) i = beta (m.tbl i).θ c := by
    intro i
    rw [hθ i]
  have h := m.prop_general_rank_unweighted_svdstack_gaussian (fun _ => c) (fun _ => beta θ c)
    (fun _ => hc) hβ le_rfl
    (topGap_of_card_le _ (by simp)) hreg hG
  rwa [hR, limitR_example_beta hc] at h

end Convergence

end StackedSVD.RankR.Example
