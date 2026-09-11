/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.WeightedMain
import StackedSVD.RankR.Frobenius

/-!
# E2: the weight-free upper bound of `thm:gen_rank_weight_svdstak`, at `r_i = 1`

`RankR/Weighted.lean` proves `thm:gen_rank_weight_svdstak` (`main_paper.tex:893`) under
`hrankB : rank B_R = r`, and bounds `limitRW W` only for an admissible `W` (a top-`r` eigengap
of `W A_{β,R} Wᵀ` and `0 < λ_{r-1}`). This file closes both gaps with one deterministic
inequality, which is item E2 of `paper_edits.md` and the rank-`r` twin of `SVDStack/Rayleigh.lean`.

The mathematics. Write `Ṽ` for the `M × d` matrix of the per-table estimates, `g = Ṽ V` and
`G = Ṽ Ṽᵀ`. For every weight matrix `W`, at every `N` and every `ω` with `det G ≠ 0`,

    tr((W g)ᵀ (specInvTop (W G Wᵀ) r) (W g)) ≤ tr(gᵀ G⁻¹ g).

The right side does not read `W`. It converges in probability to
`tr(B_Rᵀ A_{β,R}⁻¹ B_R)`, which is the paper's `L⋆`. The bound gives the optimality half of
the paper's claim over **all** weight matrices, data dependent ones included, and the
truncated performance at `W⋆` gives the attainment half with no hypothesis on `rank B_R`.

## Content

1. `specInvTop_mul_mul_self` (`P A P = P`) and `trace_specInvTop_conj_le_inv`, the
   deterministic inequality; `UnalignedModel.rowBoundR` and `perfRW_le_rowBoundR`.
2. `rbPhiR`, `rowBoundR_tendsto`, `tendsto_measure_det_gram_eq_zero`: the bound converges.
3. `limitRTrace` and `limitRTrace_eq_limitROpt`: the trace form of `L⋆`.
4. `perfRWk`, `limitRWk`, `perfRWk_tendsto`, `perfRWk_le_perfRW` and
   `thm_gen_rank_weight_svdstak_norank`: attainment at `W⋆` without `rank B_R = r`.
5. `perfRW_uniform_bound`, `perfRW_le_opt_whp`, `thm_gen_rank_weight_svdstak_full` and its
   Gaussian facade.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

/-! ### 0. Copies of two lemmas that live downstream of `RankR/Weighted.lean`

`isHermitian_specInvTop` is in `RankR/Frobenius.lean` and `specInvTop_eq_inv` is in
`RankR/Example.lean`; both files import `RankR/Weighted.lean`, so this file cannot reach them.
The two proofs are copied under primed names, exactly as `SVDStack/Rayleigh.lean` copies the
private helpers of `SVDStack/Defs.lean`. -/

section Copies

variable {p : ℕ}

/-- Copy of `isHermitian_specInvTop` (`RankR/Frobenius.lean`): `specInvTop` is symmetric,
because it is a real linear combination of the rank-one matrices `b_i b_iᵀ`. -/
theorem isHermitian_specInvTop' {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (r : ℕ) :
    (specInvTop A hA r).IsHermitian := by
  have hsym : ∀ i j, specInvTop A hA r i j = specInvTop A hA r j i := by
    intro i j
    simp only [specInvTop, Matrix.sum_apply, Matrix.smul_apply, Matrix.vecMulVec_apply,
      smul_eq_mul]
    exact Finset.sum_congr rfl fun k _ => by ring
  ext i j
  simpa [Matrix.conjTranspose_apply] using (hsym j i)

/-- Copy of the private `sum_vecMulVec_eigenvectorBasis` of `RankR/Example.lean`: `U Uᵀ = I`
read on the eigenvector basis. -/
private theorem sum_vecMulVec_eigenvectorBasis' {A : Matrix (Fin p) (Fin p) ℝ}
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

/-- Copy of `specInvTop_eq_inv` (`RankR/Example.lean`): at `p ≤ r` and `det A ≠ 0` the
truncated inverse is the full inverse. -/
theorem specInvTop_eq_inv' {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
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
  rw [key, sum_vecMulVec_eigenvectorBasis' hA]

end Copies

/-! ### 1. The deterministic inequality

`specInvTop A hA k` is the pseudo inverse of `A` on its top-`k` eigenspace, so `P A P = P`.
The trace bound then follows from a positive semidefinite square: with `T = Wᵀ P W` and
`K = G⁻¹ - T` one has `T G T = T`, `K G K = K` and `gᵀ K g = (K g)ᵀ G (K g) ⪰ 0`. -/

section Deterministic

variable {p : ℕ}

/-- The scalar identity behind `P A P = P`: `c λ c = c` for `c = λ⁻¹` or `c = 0`. It holds at
`λ = 0` too, because `0⁻¹ = 0` in Lean. -/
theorem invCoef_mul_eigenvalue_mul_invCoef (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (k : ℕ) (i : Fin p) :
    invCoef A hA k i * hA.eigenvalues i * invCoef A hA k i = invCoef A hA k i := by
  by_cases hmem : hA.eigenvalues i ∈ topEigSet A hA k
  · rw [invCoef_of_mem hmem]
    rcases eq_or_ne (hA.eigenvalues i) 0 with h | h
    · rw [h]
      simp
    · rw [inv_mul_cancel₀ h, one_mul]
  · rw [invCoef_of_notMem hmem]
    ring

/-- The eigenvectors are orthonormal, read as a dot product. -/
private theorem dotProduct_eigvec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (a b : Fin p) :
    WithLp.ofLp (hA.eigenvectorBasis a) ⬝ᵥ WithLp.ofLp (hA.eigenvectorBasis b)
      = if a = b then 1 else 0 := by
  rw [← real_inner_eq_dotProduct]
  by_cases hab : a = b
  · subst hab
    have h := real_inner_self_eq_norm_sq (hA.eigenvectorBasis a)
    rw [if_pos rfl, h, hA.eigenvectorBasis.orthonormal.1, one_pow]
  · rw [if_neg hab, hA.eigenvectorBasis.orthonormal.2 hab]

/-- `u_i ⬝ᵥ (P y) = c_i (u_i ⬝ᵥ y)`, the diagonal action of `specInvTop` read against an
eigenvector. -/
private theorem dotProduct_eigvec_specInvTop_mulVec (A : Matrix (Fin p) (Fin p) ℝ)
    (hA : A.IsHermitian) (k : ℕ) (i : Fin p) (y : Fin p → ℝ) :
    WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ (specInvTop A hA k *ᵥ y)
      = invCoef A hA k i * (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ y) := by
  classical
  rw [specInvTop_mulVec, dotProduct_sum, Finset.sum_eq_single i]
  · rw [dotProduct_smul, smul_eq_mul, dotProduct_eigvec A hA i i, if_pos rfl, mul_one]
  · intro j _ hj
    rw [dotProduct_smul, smul_eq_mul, dotProduct_eigvec A hA i j, if_neg (Ne.symm hj),
      mul_zero]
  · intro h
    exact absurd (Finset.mem_univ i) h

/-- `u_i ⬝ᵥ (A z) = λ_i (u_i ⬝ᵥ z)` for a symmetric `A`. -/
private theorem dotProduct_eigvec_mulVec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (i : Fin p) (z : Fin p → ℝ) :
    WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ (A *ᵥ z)
      = hA.eigenvalues i * (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ z) := by
  have hAt : Aᵀ = A := by
    have h := hA
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have hvm : WithLp.ofLp (hA.eigenvectorBasis i) ᵥ* A
      = hA.eigenvalues i • WithLp.ofLp (hA.eigenvectorBasis i) := by
    rw [← Matrix.mulVec_transpose, hAt]
    exact hA.mulVec_eigenvectorBasis i
  rw [Matrix.dotProduct_mulVec, hvm, smul_dotProduct, smul_eq_mul]

/-- `P A P = P` for `P = specInvTop A hA k`: on the eigenvector `u_i` the coefficient chain is
`c_i λ_i c_i = c_i`. -/
theorem specInvTop_mul_mul_self (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (k : ℕ) :
    specInvTop A hA k * A * specInvTop A hA k = specInvTop A hA k := by
  classical
  refine Matrix.ext_iff_mulVec.mpr fun y => ?_
  have hkey : ∀ i : Fin p, WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ
      (A *ᵥ (specInvTop A hA k *ᵥ y))
      = hA.eigenvalues i * (invCoef A hA k i *
        (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ y)) := by
    intro i
    rw [dotProduct_eigvec_mulVec A hA i, dotProduct_eigvec_specInvTop_mulVec A hA k i y]
  have hassoc : (specInvTop A hA k * A * specInvTop A hA k) *ᵥ y
      = specInvTop A hA k *ᵥ (A *ᵥ (specInvTop A hA k *ᵥ y)) := by
    rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]
  rw [hassoc, specInvTop_mulVec A hA k (A *ᵥ (specInvTop A hA k *ᵥ y))]
  simp only [hkey]
  rw [specInvTop_mulVec A hA k y]
  refine Finset.sum_congr rfl fun i _ => ?_
  congr 1
  have hc := invCoef_mul_eigenvalue_mul_invCoef A hA k i
  linear_combination (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ y) * hc

end Deterministic

section TraceBound

variable {M r : ℕ}

/-- **The weight-free trace bound** (`paper_edits.md`, item E2). For every `W` and every
positive semidefinite invertible `G`,
`tr((W g)ᵀ (specInvTop (W G Wᵀ) k) (W g)) ≤ tr(gᵀ G⁻¹ g)`. The right side does not read `W`.
Route: `T = Wᵀ P W` satisfies `T G T = T`, so `K = G⁻¹ - T` satisfies `K G K = K` and
`gᵀ K g = (K g)ᵀ G (K g)` is positive semidefinite. -/
theorem trace_specInvTop_conj_le_inv {G : Matrix (Fin M) (Fin M) ℝ}
    (hG : G.PosSemidef) (hdet : IsUnit G.det) (W : Matrix (Fin M) (Fin M) ℝ)
    (hS : (W * G * Wᵀ).IsHermitian) (g : Matrix (Fin M) (Fin r) ℝ) (k : ℕ) :
    Matrix.trace ((W * g)ᵀ * specInvTop (W * G * Wᵀ) hS k * (W * g))
      ≤ Matrix.trace (gᵀ * G⁻¹ * g) := by
  classical
  set P : Matrix (Fin M) (Fin M) ℝ := specInvTop (W * G * Wᵀ) hS k with hPdef
  set T : Matrix (Fin M) (Fin M) ℝ := Wᵀ * P * W with hTdef
  have hGt : Gᵀ = G := by
    have h := hG.isHermitian
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have hPt : Pᵀ = P := by
    have h := isHermitian_specInvTop' hS k
    rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
    rw [hPdef]
    exact h
  have hTt : Tᵀ = T := by
    simp only [hTdef, Matrix.transpose_mul, Matrix.transpose_transpose, hPt, Matrix.mul_assoc]
  have hTGT : T * G * T = T := by
    have hsq := specInvTop_mul_mul_self (W * G * Wᵀ) hS k
    rw [← hPdef] at hsq
    calc T * G * T = Wᵀ * (P * (W * G * Wᵀ) * P) * W := by
          simp only [hTdef, Matrix.mul_assoc]
      _ = Wᵀ * P * W := by rw [hsq]
      _ = T := hTdef.symm
  set K : Matrix (Fin M) (Fin M) ℝ := G⁻¹ - T with hKdef
  have hKt : Kᵀ = K := by
    rw [hKdef, Matrix.transpose_sub, hTt, Matrix.transpose_nonsing_inv, hGt]
  have hKGK : K * G * K = K := by
    have h1 : K * G = 1 - T * G := by
      rw [hKdef, Matrix.sub_mul, Matrix.nonsing_inv_mul _ hdet]
    have h2 : T * G * G⁻¹ = T := by
      rw [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hdet, Matrix.mul_one]
    rw [h1]
    nth_rewrite 1 [hKdef]
    rw [Matrix.sub_mul, Matrix.one_mul, Matrix.mul_sub, h2, hTGT, hKdef]
    abel
  have hpsd : (gᵀ * K * g).PosSemidef := by
    have h := hG.conjTranspose_mul_mul_same (K * g)
    rw [Matrix.conjTranspose_eq_transpose_of_trivial] at h
    have heq : (K * g)ᵀ * G * (K * g) = gᵀ * K * g := by
      rw [Matrix.transpose_mul, hKt]
      calc gᵀ * K * G * (K * g) = gᵀ * (K * G * K) * g := by simp only [Matrix.mul_assoc]
        _ = gᵀ * K * g := by rw [hKGK]
    rwa [heq] at h
  have hsplit : Matrix.trace (gᵀ * K * g)
      = Matrix.trace (gᵀ * G⁻¹ * g) - Matrix.trace (gᵀ * T * g) := by
    rw [hKdef, Matrix.mul_sub, Matrix.sub_mul, Matrix.trace_sub]
  have hLHS : (W * g)ᵀ * P * (W * g) = gᵀ * T * g := by
    rw [Matrix.transpose_mul, hTdef]
    simp only [Matrix.mul_assoc]
  rw [hLHS]
  linarith [hpsd.trace_nonneg, hsplit]

end TraceBound

/-! ### 1b. The bound on an `UnalignedModel`

`rowBoundR` is the right side of the trace bound at `G = Ṽ Ṽᵀ` and `g = Ṽ V`. It reads no
weight matrix, so `perfRW_le_rowBoundR` covers data dependent weights. -/

namespace UnalignedModel

variable {M r : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- The weight-free upper bound `tr((Ṽ V)ᵀ (Ṽ Ṽᵀ)⁻¹ (Ṽ V))` (`paper_edits.md`, item E2). The
rank-one mirror is `MultiTableModel.rowBoundClosed` of `SVDStack/Rayleigh.lean`. -/
noncomputable def rowBoundR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtV N ω)ᵀ * (m.gram N ω)⁻¹ * m.VtV N ω)

/-- **The uniform bound, deterministic form.** At every `N` and every `ω` where `Ṽ Ṽᵀ` is
invertible, no weight matrix beats `rowBoundR` (`paper_edits.md`, item E2). -/
theorem perfRW_le_rowBoundR (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (N : ℕ) (ω : Ω N) (hdet : IsUnit (m.gram N ω).det) :
    m.perfRW W N ω ≤ m.rowBoundR N ω := by
  have hS : (W * m.gram N ω * Wᵀ).IsHermitian := by
    have h := Matrix.isHermitian_mul_mul_conjTranspose (A := m.gram N ω) W
      (m.isHermitian_gram N ω)
    rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have h := trace_specInvTop_conj_le_inv (m.posSemidef_gram N ω) hdet W hS (m.VtV N ω) r
  rw [UnalignedModel.perfRW, UnalignedModel.rowBoundR,
    specInvTop_congr_mat (m.isHermitian_gramW W N ω) hS (m.gramW_eq W N ω) r,
    m.VtVW_eq W N ω]
  exact h

end UnalignedModel

/-! ### 2. `L⋆` in trace form and the limit of the bound

`limitRTrace β R = tr(B_Rᵀ A_{β,R}⁻¹ B_R)` is the value the upper bound converges to. Section 3
shows that it is the paper's `L⋆` (`limitROpt`). The continuous functional `rbPhiR` writes
`tr(Yᵀ G⁻¹ Y)` as `(det G)⁻¹` times a polynomial in the entries, which is the rank-`r` twin of
`rbPhi` of `SVDStack/Rayleigh.lean`. -/

section LimitTrace

variable {M r : ℕ}

/-- `L⋆` in the trace form the proof produces, `tr(B_Rᵀ A_{β,R}⁻¹ B_R)` (`paper_edits.md`,
item E2, "Identity"). -/
noncomputable def limitRTrace (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) : ℝ :=
  Matrix.trace ((BR β R)ᵀ * (AbetaR β R)⁻¹ * BR β R)

/-- `tr(Yᵀ C Y)` as a triple sum over the entries. -/
theorem trace_conj_eq_sum (C : Matrix (Fin M) (Fin M) ℝ) (Y : Matrix (Fin M) (Fin r) ℝ) :
    Matrix.trace (Yᵀ * C * Y)
      = ∑ k : Fin r, ∑ i : Fin M, ∑ j : Fin M, Y i k * C i j * Y j k := by
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [Finset.sum_comm]

/-- The `M × M` block of a point of `(Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ`. -/
def rbMatR (z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.of fun i j => z (Sum.inl (i, j))

/-- `tr(Yᵀ G⁻¹ Y)` as an explicit polynomial divided by `det G`. Rank-`r` twin of `rbPhi`
(`SVDStack/Rayleigh.lean`). -/
noncomputable def rbPhiR (z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ) : ℝ :=
  ((rbMatR z).det)⁻¹ * ∑ k : Fin r, ∑ i : Fin M, ∑ j : Fin M,
    z (Sum.inr (i, k)) * (rbMatR z).adjugate i j * z (Sum.inr (j, k))

/-- `rbPhiR` reads the trace form at any point that carries `G` and `Y`. No invertibility is
needed: `Matrix.inv_def` is unconditional. -/
theorem rbPhiR_eq (G : Matrix (Fin M) (Fin M) ℝ) (Y : Matrix (Fin M) (Fin r) ℝ)
    (z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ)
    (hG : ∀ i j, z (Sum.inl (i, j)) = G i j) (hY : ∀ i k, z (Sum.inr (i, k)) = Y i k) :
    rbPhiR z = Matrix.trace (Yᵀ * G⁻¹ * Y) := by
  have hmat : rbMatR z = G := by
    ext i j
    exact hG i j
  have hinv : ∀ i j, G⁻¹ i j = (G.det)⁻¹ * G.adjugate i j := by
    intro i j
    rw [Matrix.inv_def, Ring.inverse_eq_inv, Matrix.smul_apply, smul_eq_mul]
  rw [rbPhiR, hmat, trace_conj_eq_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [hinv i j, hY i k, hY j k]
  ring

/-- The `M × M` block is continuous in the point. -/
theorem continuous_rbMatR :
    Continuous (fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ => rbMatR z) :=
  continuous_matrix fun _ _ => continuous_apply _

/-- `rbPhiR` is continuous at every point whose matrix block is invertible. -/
theorem continuousAt_rbPhiR (u : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ)
    (hdet : (rbMatR u).det ≠ 0) : ContinuousAt (rbPhiR (M := M) (r := r)) u := by
  have hadj : Continuous (fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
      (rbMatR z).adjugate) := continuous_rbMatR.matrix_adjugate
  have hnum : Continuous (fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
      ∑ k : Fin r, ∑ i : Fin M, ∑ j : Fin M,
        z (Sum.inr (i, k)) * (rbMatR z).adjugate i j * z (Sum.inr (j, k))) :=
    continuous_finsetSum _ fun k _ => continuous_finsetSum _ fun i _ =>
      continuous_finsetSum _ fun j _ =>
        ((continuous_apply _).mul (hadj.matrix_elem i j)).mul (continuous_apply _)
  have hdetc : Continuous (fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ => (rbMatR z).det) :=
    continuous_rbMatR.matrix_det
  change ContinuousAt (fun z => ((rbMatR z).det)⁻¹ * ∑ k : Fin r, ∑ i : Fin M, ∑ j : Fin M,
    z (Sum.inr (i, k)) * (rbMatR z).adjugate i j * z (Sum.inr (j, k))) u
  exact (hdetc.continuousAt.inv₀ hdet).mul hnum.continuousAt

/-- The limit point: `A_{β,R}` on the left block, `B_R` on the right block. -/
noncomputable def rbLimR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ :=
  Sum.elim (fun t : Fin M × Fin M => AbetaR β R t.1 t.2) (fun t : Fin M × Fin r => BR β R t.1 t.2)

theorem rbMatR_rbLimR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    rbMatR (rbLimR β R) = AbetaR β R := rfl

theorem rbPhiR_rbLimR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    rbPhiR (rbLimR β R) = limitRTrace β R :=
  rbPhiR_eq (AbetaR β R) (BR β R) (rbLimR β R) (fun _ _ => rfl) (fun _ _ => rfl)

end LimitTrace

namespace UnalignedModel

variable {M r : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- The random point: `Ṽ Ṽᵀ` on the left block, `Ṽ V` on the right block. -/
noncomputable def rbFR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ :=
  Sum.elim (fun t : Fin M × Fin M => m.gram N ω t.1 t.2)
    (fun t : Fin M × Fin r => m.VtV N ω t.1 t.2)

theorem rbMatR_rbFR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    rbMatR (m.rbFR N ω) = m.gram N ω := rfl

theorem rbPhiR_rbFR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    rbPhiR (m.rbFR N ω) = m.rowBoundR N ω :=
  rbPhiR_eq (m.gram N ω) (m.VtV N ω) (m.rbFR N ω) (fun _ _ => rfl) (fun _ _ => rfl)

/-- `gramR` and `VtV_tendsto` in one family (the six lines inside
`thm_gen_rank_weight_svdstak_general`). -/
theorem tendstoInProbPi_rbFR (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProbPi μ (fun N ω => m.rbFR N ω) (rbLimR β m.R) := by
  rintro (t | t)
  · exact m.gramR c β hβdef law hG t.1 t.2
  · exact m.VtV_tendsto c β hβdef law t.1 t.2

/-- The determinant of `Ṽ Ṽᵀ` converges in probability to `det A_{β,R} > 0`. -/
theorem det_gram_tendsto (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => (m.gram N ω).det) ((AbetaR β m.R).det) := by
  have hcont : ContinuousAt
      (fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ => (rbMatR z).det) (rbLimR β m.R) :=
    continuous_rbMatR.matrix_det.continuousAt
  have h := TendstoInProbPi.comp_continuous hcont (m.tendstoInProbPi_rbFR c β hβdef law hG)
  rw [rbMatR_rbLimR] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  rw [m.rbMatR_rbFR N ω]

/-- The null determinant event vanishes: `det (Ṽ Ṽᵀ) →p det A_{β,R} > 0`. -/
theorem tendsto_measure_det_gram_eq_zero (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    Tendsto (fun N => μ N {ω | ¬ IsUnit (m.gram N ω).det}) atTop (𝓝 0) := by
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have hpos : 0 < (AbetaR β m.R).det :=
    (abetaR_posDef m.R (fun i => (hβ01 i).1) (fun i => (hβ01 i).2)).det_pos
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | (AbetaR β m.R).det ≤ |(m.gram N ω).det - (AbetaR β m.R).det|}) ?_
    ((m.det_gram_tendsto c β hβdef law hG) _ hpos)
  intro N ω hω
  have hω' : (m.gram N ω).det = 0 := by
    by_contra hne
    exact hω (isUnit_iff_ne_zero.mpr hne)
  change (AbetaR β m.R).det ≤ |(m.gram N ω).det - (AbetaR β m.R).det|
  rw [hω', zero_sub, abs_neg, abs_of_pos hpos]

/-- **The bound converges** (`paper_edits.md`, item E2). Rank-`r` twin of `rowBound_tendsto`
(`SVDStack/Rayleigh.lean`). No eigengap and no rank hypothesis: `rbPhiR` is continuous at the
limit point because `A_{β,R}` is positive definite. -/
theorem rowBoundR_tendsto (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.rowBoundR N ω) (limitRTrace β m.R) := by
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have hpos : 0 < (AbetaR β m.R).det :=
    (abetaR_posDef m.R (fun i => (hβ01 i).1) (fun i => (hβ01 i).2)).det_pos
  have hphi : ContinuousAt (rbPhiR (M := M) (r := r)) (rbLimR β m.R) := by
    refine continuousAt_rbPhiR (rbLimR β m.R) ?_
    rw [rbMatR_rbLimR]
    exact hpos.ne'
  have h := TendstoInProbPi.comp_continuous hphi (m.tendstoInProbPi_rbFR c β hβdef law hG)
  rw [rbPhiR_rbLimR] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact m.rbPhiR_rbFR N ω

end UnalignedModel

/-! ### 3. `limitRTrace = L⋆`

The truncation index of `specInvTop` is free in this section. At a matrix of the shape
`C Cᵀ + I` every sorted eigenvalue at index `rank (C Cᵀ)` or above equals `1`, so the trace
form does not see the truncation index once the index is at least that rank. Taking the index
`M` gives the full inverse (`specInvTop_eq_inv'`), and taking it `r` gives `limitRW`. -/

section TruncationFree

variable {p s : ℕ}

/-- **The truncation index does not matter above the rank.** At `S = C Cᵀ + I`, the trace
`tr(Cᵀ (specInvTop S k) C)` is the same for every `k` at or above `rank (C Cᵀ)`: an eigenvalue
outside the top `k` equals `1` and its term carries the factor `λ - 1 = 0`. Generalization of
`trace_specInvTop_add_one` (`LinAlg/Eigen.lean`), where the index is the column count. -/
theorem trace_specInvTop_add_one_congr (C : Matrix (Fin p) (Fin s) ℝ)
    (hA : (C * Cᵀ + 1).IsHermitian) {k l : ℕ}
    (hk : (C * Cᵀ).rank ≤ k) (hl : (C * Cᵀ).rank ≤ l) :
    Matrix.trace (Cᵀ * specInvTop (C * Cᵀ + 1) hA k * C)
      = Matrix.trace (Cᵀ * specInvTop (C * Cᵀ + 1) hA l * C) := by
  classical
  have hP : (C * Cᵀ).IsHermitian := by simpa using Matrix.isHermitian_mul_conjTranspose_self C
  have hPSD : (C * Cᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose C
  have hnn : ∀ j, 0 ≤ hP.eigenvalues₀ j := by
    intro j
    rw [← eigenvalues_eigIdx hP]
    exact hPSD.eigenvalues_nonneg _
  have hev : ∀ j, hA.eigenvalues₀ j = hP.eigenvalues₀ j + 1 := eigenvalues₀_add_one hP hA
  have hnorm : ∀ i : Fin p,
      ∑ j, (Cᵀ *ᵥ (WithLp.ofLp (hA.eigenvectorBasis i))) j ^ 2 = hA.eigenvalues i - 1 := by
    intro i
    rw [sum_sq_transpose_mulVec]
    have huu : (WithLp.ofLp (hA.eigenvectorBasis i)) ⬝ᵥ
        (WithLp.ofLp (hA.eigenvectorBasis i)) = 1 := by
      have h := real_inner_self_eq_norm_sq (hA.eigenvectorBasis i)
      rw [real_inner_eq_dotProduct] at h
      rw [h, hA.eigenvectorBasis.orthonormal.1, one_pow]
    have hmv := hA.mulVec_eigenvectorBasis i
    rw [Matrix.add_mulVec, Matrix.one_mulVec] at hmv
    rw [eq_sub_of_add_eq hmv, sub_dotProduct, smul_dotProduct, huu]
    simp
  rw [trace_specInvTop_conj, trace_specInvTop_conj]
  simp only [hnorm]
  refine Finset.sum_congr rfl fun i _ => ?_
  obtain ⟨j, hj⟩ : ∃ j : Fin (Fintype.card (Fin p)), hA.eigenvalues i = hA.eigenvalues₀ j :=
    ⟨_, rfl⟩
  by_cases hjr : (j : ℕ) < (C * Cᵀ).rank
  · have hmemk : hA.eigenvalues i ∈ topEigSet (C * Cᵀ + 1) hA k := ⟨j, by omega, hj.symm⟩
    have hmeml : hA.eigenvalues i ∈ topEigSet (C * Cᵀ + 1) hA l := ⟨j, by omega, hj.symm⟩
    rw [Set.indicator_of_mem hmemk, Set.indicator_of_mem hmeml]
  · have hz : hP.eigenvalues₀ j = 0 :=
      eigenvalues₀_eq_zero_of_rank_le hP hnn le_rfl (by omega)
    have hone : hA.eigenvalues i = 1 := by rw [hj, hev j, hz, zero_add]
    rw [hone]
    simp

end TruncationFree

/-! ### 3b. The truncated weighted limit and the identity `L⋆ = tr(B_Rᵀ A⁻¹ B_R)` -/

section OptIdentity

variable {M r : ℕ}

/-- The weighted limit with the truncation index `k` free. At `k = r` it is `limitRW`. -/
noncomputable def limitRWk (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) (k : ℕ) : ℝ :=
  Matrix.trace ((BRW W β R)ᵀ * specInvTop (AbetaRW W β R) (isHermitian_AbetaRW W β R) k *
    BRW W β R)

theorem limitRWk_r (W : Matrix (Fin M) (Fin M) ℝ) (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : limitRWk W β R r = limitRW W β R := rfl

/-- `rank B_R ≤ r`, because `B_R` has `r` columns. -/
theorem rank_BR_le (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (BR β R).rank ≤ r := by
  simpa using Matrix.rank_le_card_width (BR β R)

/-- `W⋆` is invertible (the determinant computation inside `rank_BRW_optWR`). -/
theorem isUnit_det_optWR {β : Fin M → ℝ} (h1 : ∀ i, β i ^ 2 < 1) : IsUnit (optWR β).det := by
  have hne : ∀ i, optW β i ≠ 0 := by
    intro i
    have hp : (0 : ℝ) < 1 - β i ^ 2 := by linarith [h1 i]
    have hpos : 0 < optW β i := by
      rw [optW]
      positivity
    exact hpos.ne'
  rw [optWR, Matrix.det_diagonal]
  exact isUnit_iff_ne_zero.mpr (Finset.prod_ne_zero_iff.mpr fun i _ => hne i)

/-- `A_{β,R}⁻¹ = W⋆ (W⋆ A_{β,R} W⋆ᵀ)⁻¹ W⋆`, which is what turns `tr(B_Rᵀ A⁻¹ B_R)` into a
statement about `W⋆`. -/
theorem optWR_mul_inv_abetaRW_mul_optWR {β : Fin M → ℝ}
    (R : Fin M → EuclideanSpace ℝ (Fin r)) (h1 : ∀ i, β i ^ 2 < 1) :
    optWR β * (AbetaRW (optWR β) β R)⁻¹ * optWR β = (AbetaR β R)⁻¹ := by
  have hWu : IsUnit (optWR β).det := isUnit_det_optWR h1
  have hW1 : optWR β * (optWR β)⁻¹ = 1 := Matrix.mul_nonsing_inv _ hWu
  have hW2 : (optWR β)⁻¹ * optWR β = 1 := Matrix.nonsing_inv_mul _ hWu
  rw [AbetaRW, transpose_optWR, Matrix.mul_inv_rev, Matrix.mul_inv_rev]
  calc optWR β * ((optWR β)⁻¹ * ((AbetaR β R)⁻¹ * (optWR β)⁻¹)) * optWR β
      = optWR β * (optWR β)⁻¹ * (AbetaR β R)⁻¹ * ((optWR β)⁻¹ * optWR β) := by
        simp only [Matrix.mul_assoc]
    _ = (AbetaR β R)⁻¹ := by rw [hW1, hW2, Matrix.one_mul, Matrix.mul_one]

/-- **The identity of item E2**: `tr(B_Rᵀ A_{β,R}⁻¹ B_R)` is the truncated weighted limit at
`W⋆` for **every** truncation index at or above `rank B_R`. No rank hypothesis. -/
theorem limitRTrace_eq_limitRWk_optWR (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1)
    {k : ℕ} (hk : (BR β R).rank ≤ k) :
    limitRTrace β R = limitRWk (optWR β) β R k := by
  have hsq : ∀ i, β i ^ 2 < 1 := fun i => by nlinarith [h0 i, h1 i]
  have hH : (AbetaRW (optWR β) β R).IsHermitian := isHermitian_AbetaRW (optWR β) β R
  have heq : AbetaRW (optWR β) β R
      = BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1 := abetaRW_optWR_eq R hsq
  have hA' : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1).IsHermitian :=
    (isHermitian_mul_transpose_self _).add Matrix.isHermitian_one
  have hrankZ : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).rank = (BR β R).rank := by
    rw [Matrix.rank_self_mul_transpose, rank_BRW_optWR R hsq]
  have hkZ : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).rank ≤ k := by rw [hrankZ]; exact hk
  have hMZ : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).rank ≤ M := by
    simpa using Matrix.rank_le_card_width
      (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ)
  have hAWpd : (AbetaRW (optWR β) β R).PosDef := abetaRW_optWR_posDef R h0 h1
  have hspec : specInvTop (AbetaRW (optWR β) β R) hH M = (AbetaRW (optWR β) β R)⁻¹ :=
    specInvTop_eq_inv' hH le_rfl hAWpd.det_pos.ne'
  have hLT : limitRTrace β R
      = Matrix.trace ((BRW (optWR β) β R)ᵀ * (AbetaRW (optWR β) β R)⁻¹ *
        BRW (optWR β) β R) := by
    rw [limitRTrace, ← optWR_mul_inv_abetaRW_mul_optWR R hsq, BRW, Matrix.transpose_mul,
      transpose_optWR]
    simp only [Matrix.mul_assoc]
  rw [hLT, ← hspec, limitRWk,
    specInvTop_congr_mat hH hA' heq M, specInvTop_congr_mat hH hA' heq k]
  exact trace_specInvTop_add_one_congr _ hA' hMZ hkZ

/-- **`L⋆` in trace form** (`paper_edits.md`, item E2, "Identity"):
`tr(B_Rᵀ A_{β,R}⁻¹ B_R) = r - ∑_{ℓ=1}^{r} λ_{r̃+1-ℓ}(A^{-1/2} D A^{-1/2})`, for every rank of
`B_R`. -/
theorem limitRTrace_eq_limitROpt (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) (hrM : r ≤ Fintype.card (Fin M)) :
    limitRTrace β R = limitROpt β R hrM := by
  rw [limitRTrace_eq_limitRWk_optWR β R h0 h1 (rank_BR_le β R), limitRWk_r]
  exact limitRW_optWR β R h0 h1 hrM

end OptIdentity

/-! ### 4. Attainment at `W⋆` without `rank B_R = r`

The truncated performance `perfRWk W k` uses `k` eigenvectors instead of `r`. It is monotone
in `k` and it converges whenever `W A_{β,R} Wᵀ` has a gap at the index `k`. At `W⋆` the gap at
`k = rank B_R` holds with no hypothesis, because `W⋆ A_{β,R} W⋆ᵀ = Z Zᵀ + I` and `Z Zᵀ` has
exactly `rank B_R` nonzero eigenvalues. -/

section TruncationMonotone

variable {p q : ℕ}

/-- `specInvTop A hA 0 = 0`: the set `topEigSet A hA 0` is empty. -/
theorem specInvTop_zero (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) :
    specInvTop A hA 0 = 0 := by
  rw [specInvTop_eq_sum]
  refine Finset.sum_eq_zero fun i _ => ?_
  rw [invCoef_of_notMem, zero_smul]
  rintro ⟨j, hj, -⟩
  omega

/-- The trace form is monotone in the truncation index at a positive semidefinite matrix:
`topEigSet A hA k ⊆ topEigSet A hA l` for `k ≤ l`, and the extra coefficients `λ⁻¹` are
nonnegative. -/
theorem trace_specInvTop_conj_mono {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    (hPSD : S.PosSemidef) {k l : ℕ} (hkl : k ≤ l) (Y : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Yᵀ * specInvTop S hS k * Y)
      ≤ Matrix.trace (Yᵀ * specInvTop S hS l * Y) := by
  rw [trace_specInvTop_conj, trace_specInvTop_conj]
  refine Finset.sum_le_sum fun i _ => ?_
  have hs : (0 : ℝ) ≤ ∑ j, (Yᵀ *ᵥ (WithLp.ofLp (hS.eigenvectorBasis i))) j ^ 2 :=
    Finset.sum_nonneg fun j _ => sq_nonneg _
  have hnn : 0 ≤ hS.eigenvalues i := hPSD.eigenvalues_nonneg i
  by_cases hmem : hS.eigenvalues i ∈ topEigSet S hS k
  · obtain ⟨j, hj, hje⟩ := hmem
    have hmemk : hS.eigenvalues i ∈ topEigSet S hS k := ⟨j, hj, hje⟩
    have hmeml : hS.eigenvalues i ∈ topEigSet S hS l := ⟨j, lt_of_lt_of_le hj hkl, hje⟩
    rw [Set.indicator_of_mem hmemk, Set.indicator_of_mem hmeml]
  · rw [Set.indicator_of_notMem hmem, zero_mul]
    by_cases hmem2 : hS.eigenvalues i ∈ topEigSet S hS l
    · rw [Set.indicator_of_mem hmem2]
      exact mul_nonneg (inv_nonneg.mpr hnn) hs
    · rw [Set.indicator_of_notMem hmem2, zero_mul]

end TruncationMonotone

section GapAtOpt

variable {M r : ℕ}

/-- The top-`k` eigengap of `W⋆ A_{β,R} W⋆ᵀ` at `k = rank B_R`, with **no** hypothesis on the
rank. Generalization of `topGap_optWR_of_rankBR` (`RankR/Weighted.lean`), whose index is `r`
and which needs `rank B_R = r`. -/
theorem topGap_optWR_rank {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) :
    TopGap (AbetaRW (optWR β) β R) (isHermitian_AbetaRW (optWR β) β R) (BR β R).rank := by
  have hsq : ∀ i, β i ^ 2 < 1 := fun i => by nlinarith [h0 i, h1 i]
  have hPrank : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).rank = (BR β R).rank := by
    rw [Matrix.rank_self_mul_transpose, rank_BRW_optWR R hsq]
  have hP : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).IsHermitian :=
    isHermitian_mul_transpose_self _
  have hPSD : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose (BRW (optWR β) β R)
  have hnn : ∀ k, 0 ≤ hP.eigenvalues₀ k := by
    intro k
    rw [← eigenvalues_eigIdx hP]
    exact hPSD.eigenvalues_nonneg _
  have hP1 : (BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1).IsHermitian :=
    hP.add Matrix.isHermitian_one
  have heq : AbetaRW (optWR β) β R
      = BRW (optWR β) β R * (BRW (optWR β) β R)ᵀ + 1 := abetaRW_optWR_eq R hsq
  have hev : ∀ j, (isHermitian_AbetaRW (optWR β) β R).eigenvalues₀ j
      = hP.eigenvalues₀ j + 1 := by
    intro j
    rw [eigenvalues₀_congr_mat (isHermitian_AbetaRW (optWR β) β R) hP1 heq]
    exact eigenvalues₀_add_one hP hP1 j
  intro k l hk hl
  rw [hev k, hev l]
  have hzero : hP.eigenvalues₀ l = 0 :=
    eigenvalues₀_eq_zero_of_rank_le hP hnn (le_of_eq hPrank) hl
  have hposk : 0 < hP.eigenvalues₀ k :=
    eigenvalues₀_pos_of_lt_rank hP hnn (le_of_eq hPrank.symm) hk
  linarith

/-- Every eigenvalue of `W⋆ A_{β,R} W⋆ᵀ` is positive: the matrix is positive definite
(`abetaRW_optWR_posDef`), so the hypothesis `hposW` holds at every index. -/
theorem pos_eigenvalues₀_optWR {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) (j : Fin (Fintype.card (Fin M))) :
    0 < (isHermitian_AbetaRW (optWR β) β R).eigenvalues₀ j := by
  rw [← eigenvalues_eigIdx]
  exact (abetaRW_optWR_posDef R h0 h1).eigenvalues_pos _

end GapAtOpt

namespace UnalignedModel

variable {M r : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- The performance of weighted svdstack truncated at `k` eigenvectors instead of `r`. At
`k = r` it is `perfRW`. -/
noncomputable def perfRWk (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    (k : ℕ) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVW W N ω)ᵀ * specInvTop (m.gramW W N ω) (m.isHermitian_gramW W N ω) k *
    m.VtVW W N ω)

theorem perfRWk_r (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ) (N : ℕ)
    (ω : Ω N) : m.perfRWk W r N ω = m.perfRW W N ω := rfl

/-- The truncated performance is monotone in the truncation index. -/
theorem perfRWk_le_perfRW (m : UnalignedModel μ M n d r) (W : Matrix (Fin M) (Fin M) ℝ)
    {k : ℕ} (hk : k ≤ r) (N : ℕ) (ω : Ω N) : m.perfRWk W k N ω ≤ m.perfRW W N ω :=
  trace_specInvTop_conj_mono (m.isHermitian_gramW W N ω) (m.posSemidef_gramW W N ω) hk
    (m.VtVW W N ω)

/-- `thm_gen_rank_weight_svdstak_general` with the truncation index `k` free. Same proof; only
the index inside `continuousAt_traceFun` changes. -/
theorem perfRWk_tendsto (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (W : Matrix (Fin M) (Fin M) ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    {k : ℕ} (hk : 0 < k) (hkM : k ≤ M)
    (hgapW : TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) k)
    (hposW : 0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
      ⟨k - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRWk W k N ω) (limitRWk W β m.R k) := by
  have hA : (AbetaRW W β m.R).IsHermitian := isHermitian_AbetaRW W β m.R
  have hkp : k ≤ Fintype.card (Fin M) := by simpa using hkM
  have hconv0 : TendstoInProbPi μ
      (fun N ω => Sum.elim (fun t : Fin M × Fin M => m.gram N ω t.1 t.2)
        (fun t : Fin M × Fin r => m.VtV N ω t.1 t.2))
      (Sum.elim (fun t : Fin M × Fin M => AbetaR β m.R t.1 t.2)
        (fun t : Fin M × Fin r => BR β m.R t.1 t.2)) := by
    rintro (t | t)
    · exact m.gramR c β hβdef law hG t.1 t.2
    · exact m.VtV_tendsto c β hβdef law t.1 t.2
  have hconv : TendstoInProbPi μ
      (fun N ω => Sum.elim (fun t : Fin M × Fin M => m.gramW W N ω t.1 t.2)
        (fun t : Fin M × Fin r => m.VtVW W N ω t.1 t.2))
      (Sum.elim (fun t : Fin M × Fin M => AbetaRW W β m.R t.1 t.2)
        (fun t : Fin M × Fin r => BRW W β m.R t.1 t.2)) := by
    rintro (t | t)
    · obtain ⟨i, j⟩ := t
      have hcont : Continuous fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
          ∑ b : Fin M, ∑ a : Fin M, W i a * z (Sum.inl (a, b)) * W j b :=
        continuous_finsetSum _ fun b _ => continuous_finsetSum _ fun a _ =>
          (continuous_const.mul (continuous_apply _)).mul continuous_const
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inl] at h
      have hlim : AbetaRW W β m.R i j
          = ∑ b : Fin M, ∑ a : Fin M, W i a * AbetaR β m.R a b * W j b :=
        mul_mul_transpose_apply W (AbetaR β m.R) i j
      rw [Sum.elim_inl, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inl]
      rw [m.gramW_eq, mul_mul_transpose_apply]
    · obtain ⟨i, kk⟩ := t
      have hcont : Continuous fun z : (Fin M × Fin M) ⊕ (Fin M × Fin r) → ℝ =>
          ∑ a : Fin M, W i a * z (Sum.inr (a, kk)) :=
        continuous_finsetSum _ fun a _ => continuous_const.mul (continuous_apply _)
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inr] at h
      have hlim : BRW W β m.R i kk = ∑ a : Fin M, W i a * BR β m.R a kk := Matrix.mul_apply
      rw [Sum.elim_inr, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inr]
      rw [m.VtVW_eq]
      exact (Matrix.mul_apply).symm
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hk hkp hgapW hposW (BRW W β m.R)) hconv
  rw [traceFun_eq hA (BRW W β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gramW W N ω) (m.VtVW W N ω)

/-- **`thm:gen_rank_weight_svdstak`, attainment half, with no rank hypothesis**
(`paper_edits.md`, item E2, "Attainment"). The performance at `W⋆` is squeezed between the
truncated performance at `k = rank B_R`, which converges because the gap at that index is
automatic, and `rowBoundR`, which converges with no gap at all. -/
theorem thm_gen_rank_weight_svdstak_norank (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hrM : r ≤ M)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
      (limitROpt β m.R (by simpa using hrM)) := by
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i, 0 ≤ β i := fun i => (hβ01 i).1
  have h1 : ∀ i, β i < 1 := fun i => (hβ01 i).2
  have hrp : r ≤ Fintype.card (Fin M) := by simpa using hrM
  have hLeq : limitROpt β m.R hrp = limitRTrace β m.R :=
    (limitRTrace_eq_limitROpt β m.R h0 h1 hrp).symm
  rw [hLeq]
  have hlow : TendstoInProb μ (fun N ω => m.perfRWk (optWR β) (BR β m.R).rank N ω)
      (limitRTrace β m.R) := by
    rcases Nat.eq_zero_or_pos (BR β m.R).rank with hk0 | hkpos
    · have hzero : ∀ (N : ℕ) (ω : Ω N),
          m.perfRWk (optWR β) (BR β m.R).rank N ω = 0 := by
        intro N ω
        rw [UnalignedModel.perfRWk, hk0, specInvTop_zero]
        simp
      have hlim : limitRTrace β m.R = 0 := by
        rw [limitRTrace_eq_limitRWk_optWR β m.R h0 h1 (k := 0) (le_of_eq hk0), limitRWk,
          specInvTop_zero]
        simp
      rw [hlim]
      exact (TendstoInProb.const μ 0).congr fun N =>
        Filter.Eventually.of_forall fun ω => (hzero N ω).symm
    · have hkM : (BR β m.R).rank ≤ M := Matrix.rank_le_height (BR β m.R)
      have h := m.perfRWk_tendsto c β (optWR β) hβdef hkpos hkM
        (topGap_optWR_rank m.R h0 h1) (pos_eigenvalues₀_optWR m.R h0 h1 _) law hG
      rwa [← limitRTrace_eq_limitRWk_optWR β m.R h0 h1 le_rfl] at h
  have hup : TendstoInProb μ (fun N ω => m.rowBoundR N ω) (limitRTrace β m.R) :=
    m.rowBoundR_tendsto c β hc hβdef law hG
  have hnull := m.tendsto_measure_det_gram_eq_zero c β hc hβdef law hG
  refine tendstoInProb_of_subset_union₃ ?_
  intro δ hδ
  refine ⟨fun N => {ω | δ ≤ |m.perfRWk (optWR β) (BR β m.R).rank N ω - limitRTrace β m.R|},
    fun N => {ω | δ ≤ |m.rowBoundR N ω - limitRTrace β m.R|},
    fun N => {ω | ¬ IsUnit (m.gram N ω).det}, ?_, hlow δ hδ, hup δ hδ, hnull⟩
  intro N ω hω
  have hωd : δ ≤ |m.perfRW (optWR β) N ω - limitRTrace β m.R| := hω
  simp only [Set.mem_union]
  by_cases hdet : IsUnit (m.gram N ω).det
  · by_cases hA1 : δ ≤ |m.perfRWk (optWR β) (BR β m.R).rank N ω - limitRTrace β m.R|
    · exact Or.inl (Or.inl hA1)
    · by_cases hA2 : δ ≤ |m.rowBoundR N ω - limitRTrace β m.R|
      · exact Or.inl (Or.inr hA2)
      · exfalso
        have hA1' := not_le.mp hA1
        have hA2' := not_le.mp hA2
        have hle1 : m.perfRWk (optWR β) (BR β m.R).rank N ω ≤ m.perfRW (optWR β) N ω :=
          m.perfRWk_le_perfRW (optWR β) (rank_BR_le β m.R) N ω
        have hle2 : m.perfRW (optWR β) N ω ≤ m.rowBoundR N ω :=
          m.perfRW_le_rowBoundR (optWR β) N ω hdet
        rw [abs_lt] at hA1' hA2'
        have hlt : |m.perfRW (optWR β) N ω - limitRTrace β m.R| < δ := by
          rw [abs_lt]
          exact ⟨by linarith [hA1'.1], by linarith [hA2'.2]⟩
        linarith
  · exact Or.inr hdet

end UnalignedModel

/-! ### 5. The statements of item E2

`perfRW_uniform_bound` is the sentence proposed for the paper: "the bound
`‖V̂(W)ᵀ V‖_F² ≤ tr((ṼV)ᵀ(ṼṼᵀ)⁻¹ ṼV)` holds for every `W` at every `N`, so `W⋆` is optimal
among all weightings, including those for which `W A_{β,R} Wᵀ` has no eigengap and
`‖V̂(W)ᵀ V‖_F²` does not converge to a constant." The set of the bad `ω` is not assumed
measurable; only monotonicity and subadditivity of `μ N` enter. -/

namespace UnalignedModel

variable {M r : ℕ} {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- **E2, the uniform bound.** With probability tending to one no weight matrix, data
dependent or not, beats `L⋆ + ε`. Rank-`r` twin of `svdstackPerfW_uniform_bound`
(`SVDStack/Rayleigh.lean`). -/
theorem perfRW_uniform_bound (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hrM : r ≤ M)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin M) (Fin M) ℝ,
      limitROpt β m.R (by simpa using hrM) + ε ≤ m.perfRW W N ω}) atTop (𝓝 0) := by
  intro ε hε
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have hrp : r ≤ Fintype.card (Fin M) := by simpa using hrM
  have hLeq : limitROpt β m.R hrp = limitRTrace β m.R :=
    (limitRTrace_eq_limitROpt β m.R (fun i => (hβ01 i).1) (fun i => (hβ01 i).2) hrp).symm
  have hrb := m.rowBoundR_tendsto c β hc hβdef law hG
  have hnull := m.tendsto_measure_det_gram_eq_zero c β hc hβdef law hG
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ¬ IsUnit (m.gram N ω).det} ∪
      {ω | ε ≤ |m.rowBoundR N ω - limitRTrace β m.R|}) ?_
    (tendsto_measure_zero_union hnull (hrb ε hε))
  intro N ω hω
  obtain ⟨W, hW⟩ := hω
  simp only [Set.mem_union]
  by_cases hdet : IsUnit (m.gram N ω).det
  · right
    have hb := m.perfRW_le_rowBoundR W N ω hdet
    have hW' : limitRTrace β m.R + ε ≤ m.perfRW W N ω := by
      rw [← hLeq]
      exact hW
    change ε ≤ |m.rowBoundR N ω - limitRTrace β m.R|
    have hle : ε ≤ m.rowBoundR N ω - limitRTrace β m.R := by linarith
    exact le_trans hle (le_abs_self _)
  · left
    exact hdet

/-- The per-`W` form of the uniform bound, for one fixed weight matrix. -/
theorem perfRW_le_opt_whp (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hrM : r ≤ M)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise)
    (W : Matrix (Fin M) (Fin M) ℝ) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | limitROpt β m.R (by simpa using hrM) + ε ≤ m.perfRW W N ω}) atTop (𝓝 0) := by
  intro ε hε
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ∃ W' : Matrix (Fin M) (Fin M) ℝ,
      limitROpt β m.R (by simpa using hrM) + ε ≤ m.perfRW W' N ω}) ?_
    (m.perfRW_uniform_bound c β hc hβdef hrM law hG ε hε)
  intro N ω hω
  exact ⟨W, hω⟩

/-- **`thm:gen_rank_weight_svdstak` in the form of item E2.** Three conclusions in one
declaration and **no hypothesis on `rank B_R`**:

1. the performance at `W⋆ = D^{-1/2}` tends to `L⋆`;
2. for every admissible `W` the performance tends to `limitRW W ≤ L⋆`;
3. with probability tending to one no weight matrix at all beats `L⋆ + ε`.

Conjunct 3 replaces the paper's footnote (`main_paper.tex:888`), which excludes every `W`
without an eigengap. Conjunct 1 drops the hypothesis `β_ij > 0` of `main_paper.tex:893` and
the removal remark after it. -/
theorem thm_gen_rank_weight_svdstak_full (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hr : 0 < r) (hrM : r ≤ M)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
        (limitROpt β m.R (by simpa using hrM)) ∧
      (∀ W : Matrix (Fin M) (Fin M) ℝ,
        TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r →
        0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRW W N ω) (limitRW W β m.R) ∧
          limitRW W β m.R ≤ limitROpt β m.R (by simpa using hrM)) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin M) (Fin M) ℝ,
        limitROpt β m.R (by simpa using hrM) + ε ≤ m.perfRW W N ω}) atTop (𝓝 0) :=
  ⟨m.thm_gen_rank_weight_svdstak_norank c β hc hβdef hrM law hG,
    fun W hgapW hposW =>
      ⟨m.thm_gen_rank_weight_svdstak_general c β W hβdef hr hrM hgapW hposW law hG,
        m.thm_gen_rank_weight_svdstak_opt c β hc hβdef hr hrM W hgapW hposW⟩,
    m.perfRW_uniform_bound c β hc hβdef hrM law hG⟩

/-- **Item E2, Layer 2 form.** The same three conclusions from the proportional regime of
each table and the joint Gaussian law, with no `SingleTableLaw` hypothesis. Same discharge as
`thm_gen_rank_weight_svdstak_gaussian` (`RankR/Weighted.lean`). -/
theorem thm_gen_rank_weight_svdstak_full_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hr : 0 < r) (hrM : r ≤ M)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
        (limitROpt β m.R (by simpa using hrM)) ∧
      (∀ W : Matrix (Fin M) (Fin M) ℝ,
        TopGap (AbetaRW W β m.R) (isHermitian_AbetaRW W β m.R) r →
        0 < (isHermitian_AbetaRW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRW W N ω) (limitRW W β m.R) ∧
          limitRW W β m.R ≤ limitROpt β m.R (by simpa using hrM)) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin M) (Fin M) ℝ,
        limitROpt β m.R (by simpa using hrM) + ε ≤ m.perfRW W N ω}) atTop (𝓝 0) :=
  m.thm_gen_rank_weight_svdstak_full c β hc hβdef hr hrM
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG

end UnalignedModel

end StackedSVD
