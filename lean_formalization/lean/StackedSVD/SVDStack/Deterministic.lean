/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Defs

/-!
# `thm:svd_stack_general`: the deterministic and model-free steps

The second module of `StackedSVD.SVDStack`. Nothing here reads the model
`MultiTableModel`: the results are facts about `A_β`, or limit lemmas for any sequence of
random symmetric `M × M` matrices with an entrywise limit.

## Content

1. `beta_mem_Ico` (`0 ≤ β < 1`) and `one_le_lamMax_Abeta` (`λ_max(A_β) ≥ 1`).
2. `one_lt_card_fin_of_ne` and the gap bound `abeta_gap` (`λ₁ - λ₂ ≥ β_i β_j`,
   `main_paper.tex:1195`), step (iv) of the paper's proof.
3. Section `EntrywiseEigenvec`: `lamMax_tendstoInProb`, `topProj_overlap_tendsto`,
   `topSimple_whp_of_tendsto`, `lamMax_gt_half_whp_of_tendsto`, and the special case
   `lem_entrywise_conv_eigenvec`. These take any symmetric limit matrix, so the weighted
   svdstack theorem (limit `W A_β W`) can reuse them. `opNorm_sub_tendstoInProb`, which they
   all rest on, moved to `LinAlg/TopProjPerturb.lean` on 2026-08-30 so that
   `LinAlg/SpecProjPerturb.lean` can use it; the name and the statement did not change.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

section SpectralHelpers

variable {d : ℕ}

private theorem inner_toOp_left (A : Matrix (Fin d) (Fin d) ℝ)
    (x y : EuclideanSpace ℝ (Fin d)) :
    ⟪toOp A x, y⟫_ℝ = A.mulVec (WithLp.ofLp x) ⬝ᵥ WithLp.ofLp y := by
  rw [real_inner_eq_dotProduct]
  rfl

/-- Rayleigh form of a rank-one plus diagonal matrix:
`wᵀ (a aᵀ + diag D) w = (aᵀ w)² + ∑_k D_k w_k²`. Both `A_β` and `A_{β,w}` have this shape. -/
theorem inner_vecMulVec_add_diagonal {p : ℕ} (a D : Fin p → ℝ)
    (w : EuclideanSpace ℝ (Fin p)) :
    ⟪toOp (Matrix.vecMulVec a a + Matrix.diagonal D) w, w⟫_ℝ =
      (a ⬝ᵥ WithLp.ofLp w) ^ 2 + ∑ k, D k * WithLp.ofLp w k ^ 2 := by
  rw [inner_toOp_left, Matrix.add_mulVec, add_dotProduct]
  congr 1
  · simp only [dotProduct, Matrix.mulVec, Matrix.vecMulVec_apply, pow_two, Finset.sum_mul,
      Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun l _ => by ring
  · simp only [dotProduct, Matrix.mulVec_diagonal, pow_two]
    exact Finset.sum_congr rfl fun k _ => by ring

/-- Rayleigh form of `A_β`: `wᵀ A_β w = (βᵀ w)² + ∑_k (1 - β_k²) w_k²`. -/
theorem inner_Abeta {p : ℕ} (β : Fin p → ℝ) (w : EuclideanSpace ℝ (Fin p)) :
    ⟪toOp (Abeta β) w, w⟫_ℝ =
      (β ⬝ᵥ WithLp.ofLp w) ^ 2 + ∑ k, (1 - β k ^ 2) * WithLp.ofLp w k ^ 2 := by
  unfold Abeta
  exact inner_vecMulVec_add_diagonal β _ w

/-- `0 ≤ β < 1` for `β = beta θ c` with `c > 0`: `betaSq < 1` because `θ⁴ - c < θ⁴ + θ²`. -/
theorem beta_mem_Ico {θ c : ℝ} (hc : 0 < c) : 0 ≤ beta θ c ∧ beta θ c < 1 := by
  refine ⟨Real.sqrt_nonneg _, ?_⟩
  have h : betaSq θ c < 1 := by
    unfold betaSq
    split_ifs with hθ
    · have hpos : 0 < θ ^ 4 + θ ^ 2 :=
        add_pos_of_pos_of_nonneg (lt_trans hc hθ) (by positivity)
      rw [div_lt_one hpos]
      linarith [sq_nonneg θ]
    · exact zero_lt_one
  rw [beta, Real.sqrt_lt' zero_lt_one, one_pow]
  exact h

/-- `λ_max(A_β) ≥ 1`: the Rayleigh quotient at `e_i` is `(A_β)_{ii} = 1`. -/
theorem one_le_lamMax_Abeta {p : ℕ} (β : Fin p → ℝ) (i : Fin p) :
    1 ≤ lamMax (Abeta β) (isHermitian_Abeta β) := by
  set u : Fin p → ℝ := Pi.single i 1 with hudef
  set w : EuclideanSpace ℝ (Fin p) := WithLp.toLp 2 u with hwdef
  have hofLp : WithLp.ofLp w = u := rfl
  have hsum : ∀ f : Fin p → ℝ, ∑ k, f k * u k ^ 2 = f i := by
    intro f
    rw [Finset.sum_eq_single i]
    · simp [hudef]
    · intro k _ hk
      simp [hudef, Pi.single_eq_of_ne hk]
    · intro h
      exact absurd (Finset.mem_univ i) h
  have hnormw : ‖w‖ ^ 2 = 1 := by
    have h := hsum (fun _ => 1)
    simp only [one_mul] at h
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, hofLp]
    simp only [dotProduct, ← pow_two, h]
  have hRay : ⟪toOp (Abeta β) w, w⟫_ℝ = 1 := by
    rw [inner_Abeta, hofLp]
    have h1 : β ⬝ᵥ u = β i := by simp [hudef, dotProduct_single]
    rw [h1, hsum (fun k => 1 - β k ^ 2)]
    ring
  have h := inner_toOp_self_le (Abeta β) (isHermitian_Abeta β) w
  rw [hRay, hnormw] at h
  linarith

/-! ### The weighted limit matrix `A_{β,w} = W A_β W`

`AbetaW` and `isHermitian_AbetaW` are in `SVDStack/Defs.lean`. The three facts below are the
weighted twins of `inner_Abeta` and `one_le_lamMax_Abeta`; `SVDStack/Weighted.lean` uses them
for `abetaW_gap` and for the positivity of `λ_max(A_{β,w})`. -/

/-- `A_{β,w} = (Wβ)(Wβ)ᵀ + diag(w_l² (1 - β_l²))`, the rank-one plus diagonal form. -/
theorem abetaW_eq {p : ℕ} (w β : Fin p → ℝ) :
    AbetaW w β = Matrix.vecMulVec (w * β) (w * β)
      + Matrix.diagonal fun l => w l ^ 2 * (1 - β l ^ 2) := by
  ext i j
  have hAW : AbetaW w β i j = w i * Abeta β i j * w j := by
    simp [AbetaW, Matrix.mul_apply, Matrix.diagonal, Matrix.of_apply, Finset.sum_ite_eq,
      Finset.sum_ite_eq', mul_comm, mul_left_comm]
  rw [hAW, Abeta]
  by_cases h : i = j
  · subst h
    simp only [Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply_eq, Pi.mul_apply]
    ring
  · simp only [Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply_ne _ h,
      add_zero, Pi.mul_apply]
    ring

/-- Rayleigh form of `A_{β,w}`: `yᵀ A_{β,w} y = ((Wβ)ᵀ y)² + ∑_l w_l²(1 - β_l²) y_l²`. -/
theorem inner_AbetaW {p : ℕ} (w β : Fin p → ℝ) (y : EuclideanSpace ℝ (Fin p)) :
    ⟪toOp (AbetaW w β) y, y⟫_ℝ = ((w * β) ⬝ᵥ WithLp.ofLp y) ^ 2 +
      ∑ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp y l ^ 2 := by
  rw [abetaW_eq]
  exact inner_vecMulVec_add_diagonal _ _ y

/-- `λ_max(A_{β,w}) ≥ w_k²`: the Rayleigh quotient at `e_k` is `(A_{β,w})_{kk} = w_k²`. The
weighted twin of `one_le_lamMax_Abeta`, which is the case `w = 1`. -/
theorem sq_le_lamMax_AbetaW {p : ℕ} (w β : Fin p → ℝ) (k : Fin p) :
    w k ^ 2 ≤ lamMax (AbetaW w β) (isHermitian_AbetaW w β) := by
  set u : Fin p → ℝ := Pi.single k 1 with hudef
  set y : EuclideanSpace ℝ (Fin p) := WithLp.toLp 2 u with hydef
  have hofLp : WithLp.ofLp y = u := rfl
  have hsum : ∀ f : Fin p → ℝ, ∑ l, f l * u l ^ 2 = f k := by
    intro f
    rw [Finset.sum_eq_single k]
    · simp [hudef]
    · intro l _ hl
      simp [hudef, Pi.single_eq_of_ne hl]
    · intro h
      exact absurd (Finset.mem_univ k) h
  have hnormy : ‖y‖ ^ 2 = 1 := by
    have h := hsum (fun _ => 1)
    simp only [one_mul] at h
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, hofLp]
    simp only [dotProduct, ← pow_two, h]
  have hRay : ⟪toOp (AbetaW w β) y, y⟫_ℝ = w k ^ 2 := by
    rw [inner_AbetaW, hofLp]
    have h1 : (w * β) ⬝ᵥ u = w k * β k := by simp [hudef, dotProduct_single]
    rw [h1, hsum fun l => w l ^ 2 * (1 - β l ^ 2)]
    ring
  have h := inner_toOp_self_le (AbetaW w β) (isHermitian_AbetaW w β) y
  rw [hRay, hnormy] at h
  linarith

end SpectralHelpers

/-! ### Spectral facts about `A_β` -/

/-- Two distinct indices of `Fin M` force `2 ≤ M`, in the form `eigenvalues₀` needs. -/
theorem one_lt_card_fin_of_ne {M : ℕ} {i j : Fin M} (hij : i ≠ j) :
    1 < Fintype.card (Fin M) :=
  Fintype.one_lt_card_iff_nontrivial.2 ⟨⟨i, j, hij⟩⟩

-- `hβ01` is not used by the proof below: the bound `λ₁ - λ₂ ≥ β_i β_j` needs only
-- `(A_β)_{kk} = 1` and `(A_β)_{kl} = β_k β_l`. The hypothesis stays in the statement because
-- the note asks for it (step-8 choice 4).
set_option linter.unusedVariables false in
/-- The spectral gap of `A_β` (`main_paper.tex:1195`, step (iv) of the proof of
`thm:svd_stack_general`): `λ₁(A_β) - λ₂(A_β) ≥ β_i β_j` for any two distinct indices with
`β_i > 0` and `β_j > 0`. The paper states it for the two largest entries of a sorted `β`; the
unsorted form here follows because `β₁ β₂` is the largest pairwise product of a nonnegative
sorted `β`. The second eigenvalue is `eigenvalues₀ ⟨1, _⟩`, the second entry of Mathlib's
antitone eigenvalue list, and `i ≠ j` supplies the bound `1 < M` that the index needs. -/
theorem abeta_gap {M : ℕ} (β : Fin M → ℝ) (hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1)
    {i j : Fin M} (hij : i ≠ j) (hi : 0 < β i) (hj : 0 < β j) :
    β i * β j ≤ lamMax (Abeta β) (isHermitian_Abeta β) -
      (isHermitian_Abeta β).eigenvalues₀ ⟨1, one_lt_card_fin_of_ne hij⟩ := by
  have hcard : 1 < Fintype.card (Fin M) := one_lt_card_fin_of_ne hij
  have hcardM : Fintype.card (Fin M) = M := Fintype.card_fin M
  -- step 1: `λ₁ ≥ 1 + β_i β_j`, by the Rayleigh quotient at `e_i + e_j`
  have hlam1 : 1 + β i * β j ≤ lamMax (Abeta β) (isHermitian_Abeta β) := by
    set u : Fin M → ℝ := Pi.single i 1 + Pi.single j 1 with hudef
    set w : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 u with hwdef
    have hofLp : WithLp.ofLp w = u := rfl
    have hsum : ∀ f : Fin M → ℝ, ∑ k, f k * u k ^ 2 = f i + f j := by
      intro f
      have hterm : ∀ k, f k * u k ^ 2
          = (if k = i then f i else 0) + (if k = j then f j else 0) := by
        intro k
        by_cases h1 : k = i
        · subst h1; simp [hudef, hij]
        · by_cases h2 : k = j
          · subst h2; simp [hudef, h1]
          · simp [hudef, h1, h2]
      simp [hterm, Finset.sum_add_distrib]
    have hnormw : ‖w‖ ^ 2 = 2 := by
      have h := hsum (fun _ => 1)
      simp only [one_mul] at h
      rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, hofLp]
      simp only [dotProduct, ← pow_two, h]
      norm_num
    have hRay : ⟪toOp (Abeta β) w, w⟫_ℝ = 2 + 2 * (β i * β j) := by
      rw [inner_Abeta, hofLp]
      have h1 : β ⬝ᵥ u = β i + β j := by simp [hudef, dotProduct_add, dotProduct_single]
      rw [h1, hsum (fun k => 1 - β k ^ 2)]
      ring
    have h := inner_toOp_self_le (Abeta β) (isHermitian_Abeta β) w
    rw [hRay, hnormw] at h
    linarith
  -- step 2: `λ₂ ≤ 1`, by a Rayleigh bound on a nonzero vector of the top 2-dimensional
  -- eigenspace that is orthogonal to `β`
  have hlam2 : (isHermitian_Abeta β).eigenvalues₀ ⟨1, one_lt_card_fin_of_ne hij⟩ ≤ 1 := by
    have hT : (toOp (Abeta β)).IsSymmetric :=
      Matrix.isSymmetric_toEuclideanLin_iff.mpr (isHermitian_Abeta β)
    have hn : Module.finrank ℝ (EuclideanSpace ℝ (Fin M)) = Fintype.card (Fin M) :=
      finrank_euclideanSpace
    set i1 : Fin (Fintype.card (Fin M)) := ⟨1, one_lt_card_fin_of_ne hij⟩ with hi1
    set b : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 β with hbdef
    have hb0 : b ≠ 0 := by
      intro h
      have hbi : β i = 0 := congrArg (fun y : EuclideanSpace ℝ (Fin M) => WithLp.ofLp y i) h
      linarith
    have horth : Module.finrank ℝ (ℝ ∙ b) + Module.finrank ℝ ((ℝ ∙ b)ᗮ)
        = Module.finrank ℝ (EuclideanSpace ℝ (Fin M)) :=
      Submodule.finrank_add_finrank_orthogonal _
    rw [finrank_span_singleton hb0, hn, hcardM] at horth
    have hU : Module.finrank ℝ (hT.leadingEigenSubspace hn (Nat.succ_le_of_lt i1.2))
        = (i1 : ℕ) + 1 := hT.finrank_leadingEigenSubspace hn _
    have hi1val : (i1 : ℕ) = 1 := rfl
    obtain ⟨z, hzU, hzW, hz0⟩ := Submodule.exists_ne_zero_mem_inf_of_finrank_lt_add_finrank
      (hT.leadingEigenSubspace hn (Nat.succ_le_of_lt i1.2)) ((ℝ ∙ b)ᗮ) (by
        rw [hU, hn]
        omega)
    have h1 : (isHermitian_Abeta β).eigenvalues₀ i1 ≤ (toOp (Abeta β)).rayleighQuotient z :=
      hT.eigenvalues_le_rayleighQuotient_of_mem_leadingEigenSubspace hn i1 hzU hz0
    have hbz : β ⬝ᵥ WithLp.ofLp z = 0 := by
      have hz := (Submodule.mem_orthogonal _ _).mp hzW b (Submodule.mem_span_singleton_self b)
      rw [real_inner_eq_dotProduct] at hz
      simpa [hbdef] using hz
    have h3 : ∑ k, WithLp.ofLp z k ^ 2 = ‖z‖ ^ 2 := by
      rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, dotProduct]
      exact Finset.sum_congr rfl fun k _ => by ring
    have hRz : ⟪toOp (Abeta β) z, z⟫_ℝ ≤ ‖z‖ ^ 2 := by
      rw [inner_Abeta, hbz, ← h3]
      have h2 : ∑ k, (1 - β k ^ 2) * WithLp.ofLp z k ^ 2 ≤ ∑ k, WithLp.ofLp z k ^ 2 :=
        Finset.sum_le_sum fun k _ => by nlinarith [sq_nonneg (β k), sq_nonneg (WithLp.ofLp z k)]
      linarith
    have hzn : (0 : ℝ) < ‖z‖ ^ 2 := pow_pos (norm_pos_iff.mpr hz0) 2
    have h2 : (toOp (Abeta β)).rayleighQuotient z ≤ 1 := by
      rw [LinearMap.rayleighQuotient, div_le_one hzn]
      simpa using hRz
    rw [hi1] at h1
    linarith
  linarith

section EntrywiseEigenvec

/-! ### Random symmetric matrices with an entrywise limit

General forms, for a symmetric limit `A` that need not be `A_β` (the weighted svdstack
theorem has the limit `W A_β W`). `G N ω` is a symmetric random `M × M` matrix whose entries
converge in probability to those of `A`. -/

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {M : ℕ}
variable {A : Matrix (Fin M) (Fin M) ℝ} {G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ}

/-- Entrywise convergence, as `TendstoInProbPi` on the family of entries. -/
private theorem tendstoInProbPi_entries
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j)) :
    TendstoInProbPi μ (fun N ω => fun p : Fin M × Fin M => G N ω p.1 p.2)
      (fun p => A p.1 p.2) :=
  fun p => hconv p.1 p.2

/-- Step (iii) of the paper's proof in general form: `λ_max(G_N) → λ_max(A)` (Weyl). -/
theorem lamMax_tendstoInProb (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j)) :
    TendstoInProb μ (fun N ω => lamMax (G N ω) (hsymm N ω)) (lamMax A hA) := by
  have h := TendstoInProbPi.comp_continuous
    (φ := fun w : Fin M × Fin M → ℝ => lamMax (symMat w) (isHermitian_symMat w))
    continuous_lamMax_symMat.continuousAt (tendstoInProbPi_entries hconv)
  rw [lamMax_symMat hA] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact lamMax_symMat (hsymm N ω)

/-- Step (iv) of the paper's proof in general form. If the top eigenvalue of the symmetric
limit `A` is simple, then for each fixed `x` the overlap `‖P_top(G_N) x‖²` tends to
`⟪v_max(A), x⟫²`. Route: `norm_topProj_sq_continuousAt_of_topSimple` (Davis-Kahan, entrywise
form) and the continuous mapping theorem `TendstoInProbPi.comp_continuous`. -/
theorem topProj_overlap_tendsto (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    (hsimple : TopSimple A hA) (x : EuclideanSpace ℝ (Fin M)) :
    TendstoInProb μ (fun N ω => ‖topProj (G N ω) (hsymm N ω) x‖ ^ 2)
      (⟪vMax A hA, x⟫_ℝ ^ 2) := by
  have hM : 0 < M := by
    rcases Nat.eq_zero_or_pos M with h0 | h0
    · exfalso
      subst h0
      have hle : Module.finrank ℝ (topSpace A hA)
          ≤ Module.finrank ℝ (EuclideanSpace ℝ (Fin 0)) := Submodule.finrank_le _
      rw [finrank_euclideanSpace, Fintype.card_fin] at hle
      rw [TopSimple] at hsimple
      omega
    · exact h0
  have h := TendstoInProbPi.comp_continuous
    (norm_topProj_sq_continuousAt_of_topSimple hA hsimple x) (tendstoInProbPi_entries hconv)
  rw [norm_topProj_sq_symMat hA x,
    norm_topProj_sq_eq_inner_sq hsimple (mem_topSpace_vMax hM _ _) (norm_vMax hM _ _)] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact norm_topProj_sq_symMat (hsymm N ω) x

/-- With probability tending to one the top eigenvalue of `G_N` is simple, when `A` has a
positive gap `λ₁(A) - λ₂(A)`. Route: `‖G_N - A‖ < γ / 3` w.h.p. (`opNorm_sub_tendstoInProb`)
and Davis-Kahan (`topProj_perturb`). -/
theorem topSimple_whp_of_tendsto (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    (h1 : 1 < Fintype.card (Fin M)) (hgap : 0 < lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) :
    Tendsto (fun N => μ N {ω | ¬ TopSimple (G N ω) (hsymm N ω)}) atTop (𝓝 0) := by
  have hop := opNorm_sub_tendstoInProb hA hsymm hconv
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (hop ((lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) / 3) (by positivity)) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : ¬ TopSimple (G N ω) (hsymm N ω) := hω
  change (lamMax A hA - hA.eigenvalues₀ ⟨1, h1⟩) / 3 ≤ |‖G N ω - A‖ - 0|
  rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
  by_contra hlt
  rw [not_le] at hlt
  obtain ⟨-, hs, -⟩ := topProj_perturb hA (hsymm N ω) h1 hgap le_rfl hlt.le (by linarith)
  exact hω' hs

/-- With probability tending to one `λ_max(G_N) > λ_max(A) / 2`, when `λ_max(A) > 0`. This is
the positivity of the top eigenvalue that the transfer between `Ṽᵀ Ṽ` and `Ṽ Ṽᵀ` needs. -/
theorem lamMax_gt_half_whp_of_tendsto (hA : A.IsHermitian)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    (hpos : 0 < lamMax A hA) :
    Tendsto (fun N => μ N {ω | lamMax (G N ω) (hsymm N ω) ≤ lamMax A hA / 2}) atTop (𝓝 0) := by
  have hl := lamMax_tendstoInProb hA hsymm hconv
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (hl (lamMax A hA / 2) (by positivity)) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : lamMax (G N ω) (hsymm N ω) ≤ lamMax A hA / 2 := hω
  change lamMax A hA / 2 ≤ |lamMax (G N ω) (hsymm N ω) - lamMax A hA|
  rw [abs_sub_comm, abs_of_nonneg (by linarith)]
  linarith

/-- `lem:entrywise_conv_eigenvec`, step (iv) of the paper's proof, in the form the proof of
`thm:svd_stack_general` uses. If a sequence of symmetric `M × M` random matrices converges
entrywise in probability to `A_β`, and the top eigenvalue of `A_β` is simple, then the top
eigenvector converges up to sign. The sign is not stated: the conclusion is the overlap
`‖P_top(G_N) x‖² → ⟪v_max(A_β), x⟫²` for each fixed `x`, the same projector form as `overlap`
in `Defs.lean`, which squares the sign away. The route is Davis-Kahan, whose gap hypothesis is
`abeta_gap`. Special case `A = A_β` of `topProj_overlap_tendsto`. -/
theorem lem_entrywise_conv_eigenvec (β : Fin M → ℝ)
    (G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (Abeta β i j))
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β))
    (x : EuclideanSpace ℝ (Fin M)) :
    TendstoInProb μ (fun N ω => ‖topProj (G N ω) (hsymm N ω) x‖ ^ 2)
      (⟪vMax (Abeta β) (isHermitian_Abeta β), x⟫_ℝ ^ 2) :=
  topProj_overlap_tendsto (isHermitian_Abeta β) hsymm hconv hsimple x

end EntrywiseEigenvec

end StackedSVD
