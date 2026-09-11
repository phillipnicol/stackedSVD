/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Main

/-!
# `thm:svdstack_weighted`: weighted svdstack

Fifth module of `StackedSVD.SVDStack`. It proves the weighted svdstack results of the paper
(`main_paper.tex:503`, proof at `main_paper.tex:1279`). See `notes/archive/thm_svdstack_weighted.md`
for the review note and `notes/archive/audit_weighted_2026-08-30.md` for the audit that fixed the
statements.

Weighted svdstack stacks the scaled rows `w_i v̂_iᵀ` and returns the top right singular vector
of `W Ṽ` with `W = diag(w)` (`main_paper.tex:309`, `svdstack.def`). Its Gram matrix is
`(W Ṽ)ᵀ (W Ṽ) = ∑_i w_i² P_i`, so the whole unweighted machinery applies with the limit
matrix `A_{β,w} = W A_β W` in place of `A_β`.

## Content

1. `abetaW_apply`, `abetaW_one`: entries of `A_{β,w} = W A_β W` (the matrix itself and
   `isHermitian_AbetaW` are in `SVDStack/Defs.lean`, and the three deterministic facts
   `abetaW_eq`, `inner_AbetaW`, `sq_le_lamMax_AbetaW` are in `SVDStack/Deterministic.lean`).
2. `svdstackLimitW`, the limit value, and `paperW`, the paper's displayed weights
   `eq:svdstack.weight`. `optW`, `Sval` and `svdstackLimitOpt` are in `SVDStack/Defs.lean`.
3. `abetaW_gap`, `abetaW_optW_eq`, `svdstackLimitW_eq_opt`, `svdstackLimitW_optW`,
   `svdstackLimitW_paperW`, `svdstackLimitW_le_opt`: the deterministic `M × M` facts.
4. `MultiTableModel.VtW`, `svdstackGramW`, `svdstackPerfW`, `svdstackGramW_eq`,
   `gramEntriesW`, `goodEventW`, `tendsto_measure_not_goodEventW`, `topSimple_svdstackGramW`.
5. `thm_svdstack_weighted_general` (any weights, any simple `W A_β W`),
   `thm_svdstack_weighted` and `thm_svdstack_weighted_inner` (the paper's theorem at `w⋆`),
   `thm_svdstack_weighted_paper` (at the literal weights `eq:svdstack.weight`) and
   `thm_svdstack_weighted_zero` (every `β_i = 0`, `sec:svdstack_threshold` case 2).
6. `thm_svdstack_weighted_gaussian`, the Layer 2 corollary of `thm_svdstack_weighted`: the
   same conclusion from the proportional regime and the joint Gaussian law alone. It reads
   `SpikedModel.singleTableLaw_of_gaussian` through `StackedSVD.RMT.Full`, which
   `SVDStack/Main.lean` imports, and `MultiTableModel.gaussianNoise_of_joint`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

/-! ### The weighted limit matrix `A_{β,w} = W A_β W` -/

section Weighted

variable {M : ℕ}

/-- `∑_l y_l² = ‖y‖²` for a Euclidean vector. -/
private theorem sum_sq_ofLp {p : ℕ} (y : EuclideanSpace ℝ (Fin p)) :
    ∑ l, WithLp.ofLp y l ^ 2 = ‖y‖ ^ 2 := by
  rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, dotProduct]
  exact Finset.sum_congr rfl fun l _ => by ring

/-- Entries of `A_{β,w}`: `(W A_β W)_{ij} = w_i (A_β)_{ij} w_j`. -/
theorem abetaW_apply (w β : Fin M → ℝ) (i j : Fin M) :
    AbetaW w β i j = w i * Abeta β i j * w j := by
  simp [AbetaW, Matrix.mul_apply, Matrix.diagonal, Matrix.of_apply, Finset.sum_ite_eq,
    Finset.sum_ite_eq', mul_comm, mul_left_comm]

/-- At `w = 1` the weighted objects are the unweighted ones. -/
theorem abetaW_one (β : Fin M → ℝ) : AbetaW 1 β = Abeta β := by
  simp [AbetaW]

/-- Limit value of `thm:svdstack_weighted` for a weight vector `w`
(`main_paper.tex:1286`): `(v_max(W A_β W)ᵀ W β)² / λ_max(W A_β W)`. The vector `W β` is the
pointwise product `w * β`. -/
noncomputable def svdstackLimitW (w β : Fin M → ℝ) : ℝ :=
  ((w * β) ⬝ᵥ WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β))) ^ 2
    / lamMax (AbetaW w β) (isHermitian_AbetaW w β)

/-- The paper's displayed weights `eq:svdstack.weight` (`main_paper.tex:505`):
`w_i = θ_i √((θ_i²+1)/(θ_i²+c_i)) 1{θ_i⁴ > c_i}`. Below the threshold this is `0`, while
`optW β i = 1`; both give the same limit (`main_paper.tex:1330`, `svdstackLimitW_eq_opt`). -/
noncomputable def paperW (θ c : ℝ) : ℝ :=
  if c < θ ^ 4 then θ * Real.sqrt ((θ ^ 2 + 1) / (θ ^ 2 + c)) else 0

theorem paperW_nonneg {θ c : ℝ} (hθ : 0 ≤ θ) : 0 ≤ paperW θ c := by
  rw [paperW]
  split_ifs with h
  · positivity
  · exact le_rfl

/-- `optW` is the paper's weight `eq:svdstack.weight` (`main_paper.tex:505`) on a detectable
table: `1/√(1 - β_i²) = θ_i √((θ_i² + 1)/(θ_i² + c_i))` when `θ_i⁴ > c_i`
(`main_paper.tex:521`). Below the threshold `β_i = 0`, so `optW β i = 1` while the paper's
weight is `0`; `svdstackLimitW_eq_opt` shows both give the same limit. -/
theorem optW_eq_paper_weight {θ c : ℝ} (hθ : 0 ≤ θ) (hc : 0 < c) (hthr : c < θ ^ 4) :
    1 / Real.sqrt (1 - beta θ c ^ 2) = θ * Real.sqrt ((θ ^ 2 + 1) / (θ ^ 2 + c)) := by
  have hθ2 : 0 < θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
  have hden : (0 : ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith
  have hnum : (0 : ℝ) < θ ^ 2 + c := by linarith
  have hnn : (0 : ℝ) ≤ betaSq θ c := betaSq_nonneg θ c
  have hb : beta θ c ^ 2 = (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) := by
    rw [beta, Real.sq_sqrt hnn]
    unfold betaSq
    rw [if_pos hthr]
  have h1 : 1 - beta θ c ^ 2 = (θ ^ 2 + c) / (θ ^ 4 + θ ^ 2) := by
    rw [hb, eq_div_iff hden.ne', sub_mul, div_mul_cancel₀ _ hden.ne']
    ring
  have hθe : θ * Real.sqrt ((θ ^ 2 + 1) / (θ ^ 2 + c))
      = Real.sqrt (θ ^ 2 * ((θ ^ 2 + 1) / (θ ^ 2 + c))) := by
    rw [Real.sqrt_mul (by positivity), Real.sqrt_sq hθ]
  rw [h1, one_div, ← Real.sqrt_inv, inv_div, hθe]
  congr 1
  field_simp

/-- Spectral gap of `A_{β,w}`, the weighted form of `abeta_gap`. `A_{β,w} = z zᵀ + D'` with
`z = W β` and `D' = diag(w_l² (1 - β_l²))`, so the Rayleigh quotient at `e_k` gives
`λ₁ ≥ w_k²` (`sq_le_lamMax_AbetaW`) and the orthogonal-complement argument of `abeta_gap`
gives `λ₂ ≤ max_l D'_ll`.

`hmax` asks that the index `k` also attains the maximum of `w_l² (1 - β_l²)`; it is a real
restriction, not a generalization of `abeta_gap`. At the optimal weights `w⋆` every `D'_ll`
equals `1`, so `hmax` holds for every `k`, and one detectable table is enough (`β_1 > 0` in
the paper). With `w = 1` and `β_l = 0` for some `l`, `hmax` fails and the correct route is
`abeta_gap`. The audit (`notes/archive/audit_weighted_2026-08-30.md`, section 5.3) removed the three
hypotheses `0 < w i`, `0 ≤ β i < 1` and `0 < β k` that the proof never uses. -/
theorem abetaW_gap (w β : Fin M → ℝ) (h1 : 1 < Fintype.card (Fin M)) {k : Fin M}
    (hmax : ∀ l, w l ^ 2 * (1 - β l ^ 2) ≤ w k ^ 2 * (1 - β k ^ 2)) :
    w k ^ 2 * β k ^ 2 ≤ lamMax (AbetaW w β) (isHermitian_AbetaW w β) -
      (isHermitian_AbetaW w β).eigenvalues₀ ⟨1, h1⟩ := by
  have hA : (AbetaW w β).IsHermitian := isHermitian_AbetaW w β
  have hcardM : Fintype.card (Fin M) = M := Fintype.card_fin M
  -- step 1: `λ₁ ≥ w_k²`, the Rayleigh quotient at `e_k`
  have hlam1 : w k ^ 2 ≤ lamMax (AbetaW w β) hA := sq_le_lamMax_AbetaW w β k
  -- step 2: `λ₂ ≤ max_l w_l²(1 - β_l²) = w_k²(1 - β_k²)`
  have hlam2 : hA.eigenvalues₀ ⟨1, h1⟩ ≤ w k ^ 2 * (1 - β k ^ 2) := by
    have hT : (toOp (AbetaW w β)).IsSymmetric :=
      Matrix.isSymmetric_toEuclideanLin_iff.mpr hA
    have hn : Module.finrank ℝ (EuclideanSpace ℝ (Fin M)) = Fintype.card (Fin M) :=
      finrank_euclideanSpace
    set i1 : Fin (Fintype.card (Fin M)) := ⟨1, h1⟩ with hi1
    set a : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 (w * β) with hadef
    have hi1val : (i1 : ℕ) = 1 := rfl
    have hU : Module.finrank ℝ (hT.leadingEigenSubspace hn (Nat.succ_le_of_lt i1.2))
        = (i1 : ℕ) + 1 := hT.finrank_leadingEigenSubspace hn _
    have hperp : Fintype.card (Fin M) ≤ 1 + Module.finrank ℝ ((ℝ ∙ a)ᗮ) := by
      by_cases ha0 : a = 0
      · rw [ha0, Submodule.span_zero_singleton, Submodule.bot_orthogonal_eq_top,
          finrank_top, hn]
        omega
      · have h := Submodule.finrank_add_finrank_orthogonal (K := (ℝ ∙ a))
        rw [finrank_span_singleton ha0, hn] at h
        omega
    obtain ⟨z, hzU, hzW, hz0⟩ := Submodule.exists_ne_zero_mem_inf_of_finrank_lt_add_finrank
      (hT.leadingEigenSubspace hn (Nat.succ_le_of_lt i1.2)) ((ℝ ∙ a)ᗮ) (by
        rw [hU, hn]
        omega)
    have hray : hA.eigenvalues₀ i1 ≤ (toOp (AbetaW w β)).rayleighQuotient z :=
      hT.eigenvalues_le_rayleighQuotient_of_mem_leadingEigenSubspace hn i1 hzU hz0
    have haz : (w * β) ⬝ᵥ WithLp.ofLp z = 0 := by
      have hz := (Submodule.mem_orthogonal _ _).mp hzW a (Submodule.mem_span_singleton_self a)
      rw [real_inner_eq_dotProduct] at hz
      simpa [hadef] using hz
    have hzn : (0 : ℝ) < ‖z‖ ^ 2 := pow_pos (norm_pos_iff.mpr hz0) 2
    have hRz : ⟪toOp (AbetaW w β) z, z⟫_ℝ ≤ w k ^ 2 * (1 - β k ^ 2) * ‖z‖ ^ 2 := by
      rw [inner_AbetaW, haz]
      have h2 : ∑ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp z l ^ 2
          ≤ ∑ l, w k ^ 2 * (1 - β k ^ 2) * WithLp.ofLp z l ^ 2 :=
        Finset.sum_le_sum fun l _ => mul_le_mul_of_nonneg_right (hmax l) (sq_nonneg _)
      rw [← Finset.mul_sum, sum_sq_ofLp z] at h2
      simpa using h2
    have h2 : (toOp (AbetaW w β)).rayleighQuotient z ≤ w k ^ 2 * (1 - β k ^ 2) := by
      rw [LinearMap.rayleighQuotient, div_le_iff₀ hzn]
      simpa using hRz
    rw [hi1] at hray
    linarith
  nlinarith [hlam1, hlam2]

/-- The Sherman-Morrison step of `main_paper.tex:1311`: at the optimal weights the diagonal
part of `A_{β,w}` is the identity, `A_{β,w⋆} = z zᵀ + I` with `z = w⋆ * β = β/√(1-β²)`. -/
theorem abetaW_optW_eq (β : Fin M → ℝ) (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) :
    AbetaW (optW β) β = Matrix.vecMulVec (optW β * β) (optW β * β) + 1 := by
  have hpos : ∀ k : Fin M, (0 : ℝ) < 1 - β k ^ 2 := by
    intro k
    obtain ⟨h0, h1⟩ := hβ k
    nlinarith
  have hsq : ∀ k : Fin M, optW β k ^ 2 * (1 - β k ^ 2) = 1 := by
    intro k
    rw [optW, div_pow, one_pow, Real.sq_sqrt (hpos k).le, one_div,
      inv_mul_cancel₀ (hpos k).ne']
  ext i j
  rw [abetaW_apply, Abeta]
  by_cases h : i = j
  · subst h
    simp only [Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply_eq,
      Matrix.one_apply_eq, Pi.mul_apply]
    ring_nf
    nlinarith [hsq i]
  · simp only [Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply_ne _ h,
      Matrix.one_apply_ne h, add_zero, Pi.mul_apply]
    ring

/-- The key deterministic computation of `main_paper.tex:1311` to `:1318`, in the form the
paper's closing remark (`main_paper.tex:1330`) needs. Write `z = W β` and
`D_l = w_l²(1 - β_l²)`, so that `A_{β,w} = z zᵀ + diag D`. If every `D_l ≤ 1`, with equality
on the detectable tables (`β_l ≠ 0`), then `λ_max = ‖z‖² + 1 = S + 1`, the top eigenvector is
`z/‖z‖`, and the weighted limit is `S/(S+1)`.

Both the optimal weights `w⋆` (every `D_l = 1`) and the paper's displayed weights
`eq:svdstack.weight` (`D_l = 1` above the threshold, `0` below it) satisfy the two
hypotheses, which is why they give the same limit. -/
theorem svdstackLimitW_eq_opt (w β : Fin M → ℝ)
    (hle : ∀ l, w l ^ 2 * (1 - β l ^ 2) ≤ 1)
    (heq : ∀ l, β l ≠ 0 → w l ^ 2 * (1 - β l ^ 2) = 1) :
    svdstackLimitW w β = svdstackLimitOpt β := by
  have hA : (AbetaW w β).IsHermitian := isHermitian_AbetaW w β
  set zE : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 (w * β) with hzEdef
  have hzdot : WithLp.ofLp zE = w * β := rfl
  have hzl : ∀ l, WithLp.ofLp zE l = w l * β l := fun _ => rfl
  -- `‖z‖² = S`
  have hnorm2 : ‖zE‖ ^ 2 = Sval β := by
    rw [← sum_sq_ofLp zE, Sval]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [hzl l]
    by_cases hb : β l = 0
    · rw [hb]
      simp
    · have h1 := heq l hb
      have hne : 1 - β l ^ 2 ≠ 0 := by
        intro h0
        rw [h0, mul_zero] at h1
        exact zero_ne_one h1
      rw [eq_div_iff hne]
      nlinarith [h1]
  -- degenerate case `z = 0`
  by_cases hz0 : zE = 0
  · have hS : Sval β = 0 := by rw [← hnorm2, hz0, norm_zero]; ring
    have hwb : (w * β : Fin M → ℝ) = 0 := by rw [← hzdot, hz0]; rfl
    rw [svdstackLimitW, svdstackLimitOpt, hwb, hS, zero_dotProduct]
    norm_num
  have hM : 0 < M := by
    rcases Nat.eq_zero_or_pos M with h | h
    · subst h
      exact absurd (by ext i; exact i.elim0) hz0
    · exact h
  have hzn2 : (0 : ℝ) < ‖zE‖ ^ 2 := pow_pos (norm_pos_iff.mpr hz0) 2
  set u : EuclideanSpace ℝ (Fin M) := vMax (AbetaW w β) hA with hudef
  have hu1 : ‖u‖ = 1 := norm_vMax hM _ hA
  have hAu : toOp (AbetaW w β) u = lamMax (AbetaW w β) hA • u :=
    toOp_of_mem_topSpace (mem_topSpace_vMax hM _ hA)
  -- `λ_max = ‖z‖² + 1`
  have hDz : ∀ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp zE l ^ 2 = WithLp.ofLp zE l ^ 2 := by
    intro l
    by_cases hb : β l = 0
    · rw [hzl l, hb, mul_zero]
      ring
    · rw [heq l hb, one_mul]
  have hlamge : ‖zE‖ ^ 2 + 1 ≤ lamMax (AbetaW w β) hA := by
    have hray : ⟪toOp (AbetaW w β) zE, zE⟫_ℝ = ‖zE‖ ^ 2 * ‖zE‖ ^ 2 + ‖zE‖ ^ 2 := by
      rw [inner_AbetaW]
      have h1 : (w * β) ⬝ᵥ WithLp.ofLp zE = ‖zE‖ ^ 2 := by
        rw [← hzdot, ← real_inner_eq_dotProduct, real_inner_self_eq_norm_sq]
      rw [h1, Finset.sum_congr rfl fun l _ => hDz l, sum_sq_ofLp zE]
      ring
    have h := inner_toOp_self_le (AbetaW w β) hA zE
    rw [hray] at h
    by_contra hcon
    rw [not_le] at hcon
    nlinarith [mul_lt_mul_of_pos_right hcon hzn2]
  have hlamle : lamMax (AbetaW w β) hA ≤ ‖zE‖ ^ 2 + 1 := by
    have hquad : ⟪toOp (AbetaW w β) u, u⟫_ℝ = lamMax (AbetaW w β) hA := by
      rw [hAu, real_inner_smul_left, real_inner_self_eq_norm_sq, hu1]
      ring
    rw [← hquad, inner_AbetaW]
    have h1 : (w * β) ⬝ᵥ WithLp.ofLp u = ⟪zE, u⟫_ℝ := by
      rw [real_inner_eq_dotProduct, hzdot]
    rw [h1]
    have hcs : ⟪zE, u⟫_ℝ ^ 2 ≤ ‖zE‖ ^ 2 := by
      have h := abs_real_inner_le_norm zE u
      rw [hu1, mul_one] at h
      nlinarith [abs_nonneg (⟪zE, u⟫_ℝ), sq_abs (⟪zE, u⟫_ℝ), norm_nonneg zE]
    have hD : ∑ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l ^ 2 ≤ 1 := by
      calc ∑ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l ^ 2
          ≤ ∑ l, 1 * WithLp.ofLp u l ^ 2 :=
            Finset.sum_le_sum fun l _ => mul_le_mul_of_nonneg_right (hle l) (sq_nonneg _)
        _ = ‖u‖ ^ 2 := by simpa using sum_sq_ofLp u
        _ = 1 := by rw [hu1]; norm_num
    linarith
  have hlam : lamMax (AbetaW w β) hA = ‖zE‖ ^ 2 + 1 := le_antisymm hlamle hlamge
  -- the top eigenvector: `⟪z, u⟫ z = ‖z‖² u` componentwise
  have hopc : ∀ l : Fin M, ⟪zE, u⟫_ℝ * (w l * β l)
      + w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l
      = lamMax (AbetaW w β) hA * WithLp.ofLp u l := by
    intro l
    have hstep : WithLp.ofLp (toOp (AbetaW w β) u) l
        = ⟪zE, u⟫_ℝ * (w l * β l) + w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l := by
      have hmv : WithLp.ofLp (toOp (AbetaW w β) u) l
          = ((AbetaW w β).mulVec (WithLp.ofLp u)) l := rfl
      rw [hmv, abetaW_eq, Matrix.add_mulVec, Pi.add_apply, Matrix.mulVec_diagonal,
        real_inner_eq_dotProduct, hzdot]
      congr 1
      simp only [Matrix.mulVec, dotProduct, Matrix.vecMulVec_apply, Pi.mul_apply,
        Finset.sum_mul]
      exact Finset.sum_congr rfl fun j _ => by ring
    rw [← hstep, hAu]
    rfl
  have hcomp : ∀ l : Fin M, ⟪zE, u⟫_ℝ * (w l * β l) = ‖zE‖ ^ 2 * WithLp.ofLp u l := by
    intro l
    have hr := hopc l
    rw [hlam] at hr
    by_cases hb : β l = 0
    · have hz : w l * β l = 0 := by rw [hb, mul_zero]
      rw [hz, mul_zero, zero_add] at hr
      have hlt : w l ^ 2 * (1 - β l ^ 2) - (‖zE‖ ^ 2 + 1) < 0 := by
        have := hle l
        linarith
      have h0 : (w l ^ 2 * (1 - β l ^ 2) - (‖zE‖ ^ 2 + 1)) * WithLp.ofLp u l = 0 := by
        rw [sub_mul]
        linarith
      rcases mul_eq_zero.mp h0 with h | h
      · exact absurd h (ne_of_lt hlt)
      · rw [hz, h, mul_zero, mul_zero]
    · rw [heq l hb, one_mul] at hr
      linarith
  -- `⟪z, u⟫² = ‖z‖²`
  have hinner : ⟪zE, u⟫_ℝ = (w * β) ⬝ᵥ WithLp.ofLp u := by
    rw [real_inner_eq_dotProduct, hzdot]
  have hkey : ⟪zE, u⟫_ℝ ^ 2 = ‖zE‖ ^ 2 := by
    have e1 : ∑ l, ⟪zE, u⟫_ℝ * (w l * β l) * WithLp.ofLp u l
        = ⟪zE, u⟫_ℝ * ⟪zE, u⟫_ℝ := by
      have hfac : ∑ l, ⟪zE, u⟫_ℝ * (w l * β l) * WithLp.ofLp u l
          = ⟪zE, u⟫_ℝ * ∑ l, w l * β l * WithLp.ofLp u l := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun l _ => by ring
      rw [hfac]
      congr 1
      rw [hinner, dotProduct]
      exact Finset.sum_congr rfl fun l _ => rfl
    have e2 : ∑ l, ‖zE‖ ^ 2 * WithLp.ofLp u l * WithLp.ofLp u l = ‖zE‖ ^ 2 := by
      have hfac : ∑ l, ‖zE‖ ^ 2 * WithLp.ofLp u l * WithLp.ofLp u l
          = ‖zE‖ ^ 2 * ∑ l, WithLp.ofLp u l ^ 2 := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun l _ => by ring
      rw [hfac, sum_sq_ofLp u, hu1]
      norm_num
    have e3 : ∑ l, ⟪zE, u⟫_ℝ * (w l * β l) * WithLp.ofLp u l
        = ∑ l, ‖zE‖ ^ 2 * WithLp.ofLp u l * WithLp.ofLp u l :=
      Finset.sum_congr rfl fun l _ => by rw [hcomp l]
    rw [e1] at e3
    rw [e2] at e3
    nlinarith [e3]
  -- assemble
  rw [svdstackLimitW, svdstackLimitOpt, ← hudef, ← hnorm2, hlam, ← hinner, hkey]

/-- The value of the weighted limit at the optimal weights: `S / (S + 1)`
(`main_paper.tex:1328`). Deterministic. No detectable table is needed: when `β = 0` both
sides are `0`. -/
theorem svdstackLimitW_optW (β : Fin M → ℝ) (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) :
    svdstackLimitW (optW β) β = svdstackLimitOpt β := by
  have hpos : ∀ k : Fin M, (0 : ℝ) < 1 - β k ^ 2 := by
    intro k
    obtain ⟨h0, h1⟩ := hβ k
    nlinarith
  have hsq : ∀ k : Fin M, optW β k ^ 2 * (1 - β k ^ 2) = 1 := by
    intro k
    rw [optW, div_pow, one_pow, Real.sq_sqrt (hpos k).le, one_div,
      inv_mul_cancel₀ (hpos k).ne']
  exact svdstackLimitW_eq_opt (optW β) β (fun l => (hsq l).le) (fun l _ => hsq l)

/-! ### The paper's displayed weights `eq:svdstack.weight` -/

/-- Below the detectability threshold `β = 0`. -/
theorem beta_eq_zero_of_not_thr {θ c : ℝ} (h : ¬ c < θ ^ 4) : beta θ c = 0 := by
  rw [beta, betaSq, if_neg (by simpa using h), Real.sqrt_zero]

/-- Above the detectability threshold `β > 0`. -/
theorem beta_pos_of_thr {θ c : ℝ} (hc : 0 < c) (hthr : c < θ ^ 4) : 0 < beta θ c := by
  have hθ2 : 0 < θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
  have hden : (0 : ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith
  rw [beta, Real.sqrt_pos, betaSq, if_pos hthr]
  exact div_pos (by linarith) hden

/-- The diagonal part of `A_{β,w}` at the paper's weights is `1` above the threshold. -/
theorem paperW_sq_mul_one_sub {θ c : ℝ} (hθ : 0 ≤ θ) (hc : 0 < c) (hthr : c < θ ^ 4) :
    paperW θ c ^ 2 * (1 - beta θ c ^ 2) = 1 := by
  obtain ⟨hb0, hb1⟩ := beta_mem_Ico (θ := θ) hc
  have hpos : (0 : ℝ) < 1 - beta θ c ^ 2 := by nlinarith
  rw [paperW, if_pos hthr, ← optW_eq_paper_weight hθ hc hthr, div_pow, one_pow,
    Real.sq_sqrt hpos.le, one_div, inv_mul_cancel₀ hpos.ne']

/-- The paper's displayed weights `eq:svdstack.weight` give the limit `S/(S+1)`, as the
paper's closing remark says (`main_paper.tex:1330`). Below the threshold the weight is `0`
and `optW β i = 1`, so the two weight vectors differ, but both satisfy the hypotheses of
`svdstackLimitW_eq_opt`. -/
theorem svdstackLimitW_paperW (θ c β : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (θ i) (c i)) :
    svdstackLimitW (fun i => paperW (θ i) (c i)) β = svdstackLimitOpt β := by
  have hcase : ∀ l, (c l < θ l ^ 4 ∧ paperW (θ l) (c l) ^ 2 * (1 - β l ^ 2) = 1)
      ∨ (β l = 0 ∧ paperW (θ l) (c l) = 0) := by
    intro l
    by_cases h : c l < θ l ^ 4
    · exact Or.inl ⟨h, by rw [hβdef l]; exact paperW_sq_mul_one_sub (hθ l) (hc l) h⟩
    · exact Or.inr ⟨by rw [hβdef l]; exact beta_eq_zero_of_not_thr h, by rw [paperW, if_neg h]⟩
  refine svdstackLimitW_eq_opt _ β (fun l => ?_) (fun l hb => ?_)
  · rcases hcase l with ⟨-, h⟩ | ⟨-, h⟩
    · exact h.le
    · rw [h]
      norm_num
  · rcases hcase l with ⟨-, h⟩ | ⟨h, -⟩
    · exact h
    · exact absurd h hb

/-! ### Optimality of the weights `w⋆` -/

-- `hw` is not used by the proof: every quantity sees `w` through `w_i²`. The hypothesis
-- stays because the paper takes `w ∈ R^M_{≥0}` (`main_paper.tex:300`).
set_option linter.unusedVariables false in
/-- Optimality of `w⋆` (`main_paper.tex:1296`): every weight vector gives a limit at most
`S / (S + 1)`. The weighted limit is `(xᵀ β)²` for `x = W v_max(A_{β,w}) / √λ_max(A_{β,w})`,
which satisfies `xᵀ A_β x = 1`; Cauchy-Schwarz in the `A_β` inner product then bounds
`(xᵀ β)²` by `βᵀ A_β⁻¹ β = S/(S+1)`. Together with `svdstackLimitW_optW` this is the paper's
claim that `w⋆` is an optimal weighting. The audit dropped `hwne : ∃ k, w k ≠ 0`: at `w = 0`
the left side is Lean's junk value `0/0 = 0`, which is still at most `S/(S+1)`. -/
theorem svdstackLimitW_le_opt (w β : Fin M → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) :
    svdstackLimitW w β ≤ svdstackLimitOpt β := by
  have hpos : ∀ l : Fin M, (0 : ℝ) < 1 - β l ^ 2 := by
    intro l
    obtain ⟨h0, h1⟩ := hβ l
    nlinarith
  have hS0 : 0 ≤ Sval β :=
    Finset.sum_nonneg fun l _ => div_nonneg (sq_nonneg _) (hpos l).le
  have hRHS : 0 ≤ svdstackLimitOpt β := div_nonneg hS0 (by linarith)
  have hA : (AbetaW w β).IsHermitian := isHermitian_AbetaW w β
  rcases le_or_gt (lamMax (AbetaW w β) hA) 0 with hlam | hlam
  · refine le_trans ?_ hRHS
    rw [svdstackLimitW]
    exact div_nonpos_of_nonneg_of_nonpos (sq_nonneg _) hlam
  have hM : 0 < M := by
    rcases Nat.eq_zero_or_pos M with h | h
    · exfalso
      subst h
      rw [lamMax, dif_neg (by omega)] at hlam
      exact lt_irrefl 0 hlam
    · exact h
  set lam := lamMax (AbetaW w β) hA with hlamdef
  set u : EuclideanSpace ℝ (Fin M) := vMax (AbetaW w β) hA with hudef
  have hu1 : ‖u‖ = 1 := norm_vMax hM _ hA
  have hAu : toOp (AbetaW w β) u = lam • u := toOp_of_mem_topSpace (mem_topSpace_vMax hM _ hA)
  have hquad : ⟪toOp (AbetaW w β) u, u⟫_ℝ = lam := by
    rw [hAu, real_inner_smul_left, real_inner_self_eq_norm_sq, hu1]
    ring
  have hsq : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hlam.le
  have hsqpos : 0 < Real.sqrt lam := Real.sqrt_pos.mpr hlam
  set x : Fin M → ℝ := fun l => w l * WithLp.ofLp u l / Real.sqrt lam with hxdef
  have hdot : β ⬝ᵥ x = ((w * β) ⬝ᵥ WithLp.ofLp u) / Real.sqrt lam := by
    rw [dotProduct, dotProduct, Finset.sum_div]
    refine Finset.sum_congr rfl fun l _ => ?_
    simp only [hxdef, Pi.mul_apply]
    ring
  have hxsq : ∀ l, (1 - β l ^ 2) * x l ^ 2
      = w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l ^ 2 / lam := by
    intro l
    simp only [hxdef]
    rw [div_pow, hsq]
    ring
  have hray : ((w * β) ⬝ᵥ WithLp.ofLp u) ^ 2
      + ∑ l, w l ^ 2 * (1 - β l ^ 2) * WithLp.ofLp u l ^ 2 = lam := by
    rw [← hquad, inner_AbetaW]
  have hxA : (β ⬝ᵥ x) ^ 2 + ∑ l, (1 - β l ^ 2) * x l ^ 2 = 1 := by
    rw [hdot, div_pow, hsq, Finset.sum_congr rfl fun l _ => hxsq l, ← Finset.sum_div,
      ← add_div, hray, div_self hlam.ne']
  have hCS : (β ⬝ᵥ x) ^ 2 ≤ (∑ l, (1 - β l ^ 2) * x l ^ 2) * Sval β := by
    have h := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ
      (fun l => Real.sqrt (1 - β l ^ 2) * x l) (fun l => β l / Real.sqrt (1 - β l ^ 2))
    have e1 : ∑ l, Real.sqrt (1 - β l ^ 2) * x l * (β l / Real.sqrt (1 - β l ^ 2))
        = β ⬝ᵥ x := by
      rw [dotProduct]
      refine Finset.sum_congr rfl fun l _ => ?_
      have hne : Real.sqrt (1 - β l ^ 2) ≠ 0 := (Real.sqrt_pos.mpr (hpos l)).ne'
      field_simp
    have e2 : ∑ l, (Real.sqrt (1 - β l ^ 2) * x l) ^ 2 = ∑ l, (1 - β l ^ 2) * x l ^ 2 :=
      Finset.sum_congr rfl fun l _ => by
        rw [mul_pow, Real.sq_sqrt (hpos l).le]
    have e3 : ∑ l, (β l / Real.sqrt (1 - β l ^ 2)) ^ 2 = Sval β := by
      rw [Sval]
      exact Finset.sum_congr rfl fun l _ => by
        rw [div_pow, Real.sq_sqrt (hpos l).le]
    rw [← e1, ← e2, ← e3]
    exact h
  have hLHS : svdstackLimitW w β = (β ⬝ᵥ x) ^ 2 := by
    have hnum2 : (w * β) ⬝ᵥ WithLp.ofLp u = (β ⬝ᵥ x) * Real.sqrt lam := by
      rw [hdot, div_mul_cancel₀ _ hsqpos.ne']
    rw [svdstackLimitW, ← hudef, ← hlamdef, hnum2, mul_pow, hsq,
      mul_div_assoc, div_self hlam.ne', mul_one]
  rw [hLHS, svdstackLimitOpt, le_div_iff₀ (by linarith : (0 : ℝ) < Sval β + 1)]
  nlinarith [hCS, hxA]

end Weighted

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section SVDStackWeighted

/-! ### The weighted estimator -/

/-- `Ṽ_w = W Ṽ`, the `M × d` matrix whose row `i` is `w_i v̂_iᵀ` (`main_paper.tex:515`). -/
noncomputable def VtW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin M) (Fin (d N)) ℝ :=
  Matrix.diagonal w * m.Vt N ω

/-- Entries of `Ṽ_w`: row `i` is `w_i` times row `i` of `Ṽ`. -/
theorem vtW_apply (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (i : Fin M) (p : Fin (d N)) : m.VtW w N ω i p = w i * m.vhat i N ω p := by
  rw [MultiTableModel.VtW, Matrix.diagonal_mul]
  rfl

/-- `∑_i w_i² P_i`, the Gram matrix of the weighted estimator. It equals `Ṽ_wᵀ Ṽ_w` on the
almost sure event that every table has a simple top eigenvalue (`svdstackGramW_eq`), and its
top eigenvector is `v̂_svdstack(w)`. The projector form makes the object independent of the
sign convention of `vhat`, exactly as `svdstackGram` is. -/
noncomputable def svdstackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  ∑ i, (w i ^ 2) • m.P i N ω

theorem isHermitian_svdstackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.svdstackGramW w N ω).IsHermitian := by
  change (∑ i, (w i ^ 2) • m.P i N ω)ᴴ = ∑ i, (w i ^ 2) • m.P i N ω
  rw [Matrix.conjTranspose_sum]
  exact Finset.sum_congr rfl fun i _ =>
    (m.isHermitian_P i N ω).smul (IsSelfAdjoint.all (w i ^ 2))

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-- Performance of weighted svdstack in projector form, the weighted `svdstackPerf`. -/
noncomputable def svdstackPerfW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  ‖topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω) ((m.tbl 0).v N)‖ ^ 2

end SVDStackWeighted

section SVDStackWeightedLemmas

/-! ### Consequences of the single-table laws -/

/-- `∑_i w_i² P_i = Ṽ_wᵀ Ṽ_w` on the almost sure event that every table has a simple top
eigenvalue. Same content as `svdstackGram_eq`: there `P_i = v̂_i v̂_iᵀ` (`P_eq_vecMulVec`), so
the weighted sum is `Ṽᵀ W² Ṽ = (W Ṽ)ᵀ (W Ṽ)`. -/
theorem svdstackGramW_eq (m : MultiTableModel μ M n d) (c w : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω := by
  have hall : ∀ᵐ ω ∂(μ N), ∀ i : Fin M,
      TopSimple (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) :=
    ae_all_iff.mpr fun i => (law i).topSimple N
  filter_upwards [hall] with ω hω
  ext k l
  rw [MultiTableModel.svdstackGramW, Matrix.sum_apply, Matrix.mul_apply]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.smul_apply, m.P_eq_vecMulVec i N ω (hω i), Matrix.vecMulVec_apply,
    Matrix.transpose_apply, m.vtW_apply w N ω i k, m.vtW_apply w N ω i l]
  change w i ^ 2 * (m.vhat i N ω k * m.vhat i N ω l)
      = w i * m.vhat i N ω k * (w i * m.vhat i N ω l)
  ring

/-- The weighted `M × M` Gram matrix converges entrywise to `A_{β,w}`:
`(Ṽ_w Ṽ_wᵀ)_{ij} = w_i w_j (Ṽ Ṽᵀ)_{ij} → w_i w_j (A_β)_{ij}`. Consequence of `gramEntries`
and `TendstoInProb.const_mul`. -/
theorem gramEntriesW (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) (i j : Fin M) :
    TendstoInProb μ (fun N ω => (m.VtW w N ω * (m.VtW w N ω)ᵀ) i j) (AbetaW w β i j) := by
  have h := (m.gramEntries c β hβdef law hI i j).const_mul (w i * w j)
  have hlim : w i * w j * Abeta β i j = AbetaW w β i j := by
    rw [abetaW_apply]
    ring
  rw [hlim] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  show w i * w j * ((m.Vt N ω * (m.Vt N ω)ᵀ) i j) = (m.VtW w N ω * (m.VtW w N ω)ᵀ) i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, m.vtW_apply, Finset.mul_sum]
  exact Finset.sum_congr rfl fun p _ => by
    change w i * w j * (m.Vt N ω i p * m.Vt N ω j p)
        = w i * m.vhat i N ω p * (w j * m.vhat j N ω p)
    change w i * w j * (m.vhat i N ω p * m.vhat j N ω p)
        = w i * m.vhat i N ω p * (w j * m.vhat j N ω p)
    ring

/-! ### The good event of the weighted assembly -/

/-- The weighted good event: `∑_i w_i² P_i = Ṽ_wᵀ Ṽ_w`, and `Ṽ_w Ṽ_wᵀ` has a simple positive
top eigenvalue. Weighted twin of `goodEvent`. -/
def goodEventW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) : Prop :=
  m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω ∧
    TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω)) ∧
    0 < lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))

/-- The complement of the weighted good event has vanishing probability. Route: `gramEntriesW`
gives the entrywise limit `Ṽ_w Ṽ_wᵀ → A_{β,w}`; `topSimple_whp_of_tendsto` with the gap of
`hsimple` (`gap_pos_of_topSimple`) and `lamMax_gt_half_whp_of_tendsto` with
`λ_max(A_{β,w}) ≥ w_k² > 0` (`sq_le_lamMax_AbetaW`, `hwne`) make the top eigenvalue of
`Ṽ_w Ṽ_wᵀ` simple and positive with probability tending to one; `svdstackGramW_eq` is a.s.
for each `N`. Cleanup wave 3 dropped the unused `hw : 0 ≤ w i`: every quantity here sees `w`
through `w_i²` (mechanical audit 2026-08-31, finding 8). -/
theorem tendsto_measure_not_goodEventW (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hwne : ∃ k, w k ≠ 0) (h1 : 1 < Fintype.card (Fin M))
    (hsimple : TopSimple (AbetaW w β) (isHermitian_AbetaW w β))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ m.goodEventW w N ω}) atTop (𝓝 0) := by
  obtain ⟨k, hk⟩ := hwne
  have hA : (AbetaW w β).IsHermitian := isHermitian_AbetaW w β
  have hwk : (0 : ℝ) < w k ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hApos : 0 < lamMax (AbetaW w β) hA :=
    lt_of_lt_of_le hwk (sq_le_lamMax_AbetaW w β k)
  have hγ : 0 < lamMax (AbetaW w β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    gap_pos_of_topSimple hA h1 hsimple
  have hsymm : ∀ (N : ℕ) (ω : Ω N), (m.VtW w N ω * (m.VtW w N ω)ᵀ).IsHermitian :=
    fun N ω => isHermitian_mul_transpose_self (m.VtW w N ω)
  have hS := topSimple_whp_of_tendsto hA hsymm (m.gramEntriesW c β w hβdef law hI) h1 hγ
  have hL := lamMax_gt_half_whp_of_tendsto hA hsymm (m.gramEntriesW c β w hβdef law hI) hApos
  have hzero : ∀ N, μ N {ω | ¬ m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω} = 0 :=
    fun N => ae_iff.mp (m.svdstackGramW_eq c w law N)
  have hsum : Tendsto (fun N => μ N {ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)}
      + μ N {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω) ≤ lamMax (AbetaW w β) hA / 2})
      atTop (𝓝 0) := by
    simpa using hS.add hL
  have hsub : ∀ N, {ω | ¬ m.goodEventW w N ω} ⊆
      {ω | ¬ m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω} ∪
        ({ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
            ≤ lamMax (AbetaW w β) hA / 2}) := by
    intro N ω hω
    have hω' : ¬ m.goodEventW w N ω := hω
    by_cases ha : m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω
    · right
      by_cases hb : TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
      · right
        change lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω) ≤ lamMax (AbetaW w β) hA / 2
        by_contra hc'
        rw [not_le] at hc'
        exact hω' ⟨ha, hb, by linarith⟩
      · left
        exact hb
    · left
      exact ha
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum
    (fun _ => zero_le) (fun N => ?_)
  calc μ N {ω | ¬ m.goodEventW w N ω}
      ≤ μ N ({ω | ¬ m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω} ∪
        ({ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
            ≤ lamMax (AbetaW w β) hA / 2})) := measure_mono (hsub N)
    _ ≤ μ N {ω | ¬ m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω} +
        μ N ({ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
            ≤ lamMax (AbetaW w β) hA / 2}) := measure_union_le _ _
    _ ≤ μ N {ω | ¬ m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω} +
        (μ N {ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)} +
          μ N {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
            ≤ lamMax (AbetaW w β) hA / 2}) := add_le_add le_rfl (measure_union_le _ _)
    _ = μ N {ω | ¬ TopSimple (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)} +
          μ N {ω | lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (hsymm N ω)
            ≤ lamMax (AbetaW w β) hA / 2} := by rw [hzero N, zero_add]

/-- The top eigenvalue of `∑_i w_i² P_i` is simple with probability tending to one. Weighted
form of `topSimple_svdstackGram`; it turns the projector form `svdstackPerfW` into the
paper's `|⟨v̂_svdstack(w), v⟩|²`. The audits dropped the unused `hc : 0 < c i` and the unused
`hw : 0 ≤ w i`. -/
theorem topSimple_svdstackGramW (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hwne : ∃ k, w k ≠ 0) (h1 : 1 < Fintype.card (Fin M))
    (hsimple : TopSimple (AbetaW w β) (isHermitian_AbetaW w β))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopSimple (m.svdstackGramW w N ω)
      (m.isHermitian_svdstackGramW w N ω)}) atTop (𝓝 0) := by
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.tendsto_measure_not_goodEventW c β w hβdef hwne h1 hsimple law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : ¬ TopSimple (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω) := hω
  change ¬ m.goodEventW w N ω
  rintro ⟨heq, hsimpleN, hpos⟩
  refine hω' ?_
  rw [topSimple_congr_iff heq _ (isHermitian_transpose_mul_self (m.VtW w N ω))]
  exact (topSimple_transpose_mul_iff (m.VtW w N ω) _ _
    (lamMax_mul_transpose_self_eq (m.VtW w N ω) hpos) hpos).mpr hsimpleN

/-! ### The theorems -/

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-- Step (v) of the paper's proof on the weighted good event. Weighted twin of the private
`svdstackPerf_eq_of_goodEvent` of `Main.lean`. -/
private theorem svdstackPerfW_eq_of_goodEventW (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (N : ℕ) (ω : Ω N) (h : m.goodEventW w N ω) :
    m.svdstackPerfW w N ω =
      ‖topProj (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))
        (WithLp.toLp 2 ((m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2 /
      lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω)) := by
  obtain ⟨heq, hsimple, hpos⟩ := h
  set Vt := m.VtW w N ω with hVt
  set v := (m.tbl 0).v N with hv
  have hA := isHermitian_transpose_mul_self Vt
  have hB := isHermitian_mul_transpose_self Vt
  have hlam : lamMax (Vtᵀ * Vt) hA = lamMax (Vt * Vtᵀ) hB :=
    lamMax_mul_transpose_self_eq Vt hpos
  have hposA : 0 < lamMax (Vtᵀ * Vt) hA := by
    rw [hlam]
    exact hpos
  have hsimpleA : TopSimple (Vtᵀ * Vt) hA :=
    (topSimple_transpose_mul_iff Vt hA hB hlam hpos).mpr hsimple
  have hperf : m.svdstackPerfW w N ω = ‖topProj (Vtᵀ * Vt) hA v‖ ^ 2 := by
    change ‖topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω) v‖ ^ 2 = _
    rw [topProj_congr heq (m.isHermitian_svdstackGramW w N ω) hA]
  rw [hperf, svdstackPerf_eq_closed Vt v hsimpleA hposA, svdstackPerfClosed]
  have hM : 0 < M := NeZero.pos M
  set x := vMax (Vt * Vtᵀ) hB with hxdef
  have hxn : ‖x‖ = 1 := norm_vMax hM _ hB
  have hxmem : x ∈ topSpace (Vt * Vtᵀ) hB := mem_topSpace_vMax hM _ hB
  have hBx : toOp (Vt * Vtᵀ) x = lamMax (Vt * Vtᵀ) hB • x := toOp_of_mem_topSpace hxmem
  have hden : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x) = lamMax (Vt * Vtᵀ) hB := by
    have hd : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x)
        = ⟪x, toOp (Vt * Vtᵀ) x⟫_ℝ := by
      rw [real_inner_eq_dotProduct]
      rfl
    rw [hd, hBx, real_inner_smul_right, real_inner_self_eq_norm_sq, hxn]
    ring
  have hnum : WithLp.ofLp x ⬝ᵥ Vt.mulVec (WithLp.ofLp v)
      = ⟪x, (WithLp.toLp 2 (Vt.mulVec (WithLp.ofLp v)) : EuclideanSpace ℝ (Fin M))⟫_ℝ := by
    rw [real_inner_eq_dotProduct]
  rw [hden, hnum, norm_topProj_sq_eq_inner_sq hsimple hxmem hxn]

omit [NeZero M] in
/-- The matrix `P_i` is the orthogonal projector onto the top eigenspace of `X_iᵀ X_i`, so its
own top eigenspace is that same subspace. The eigenvalues of an orthogonal projector are `0`
and `1`; `1` is attained because the top eigenspace of `X_iᵀ X_i` is nonzero. Used by the
`M = 1` branch of `thm_svdstack_weighted_paper` and by the same branch of
`thm_simple_thm1_svdstack_gaussian_full` (`SVDStack/Simple.lean`). -/
theorem topSpace_P_eq (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    topSpace (m.P i N ω) (m.isHermitian_P i N ω)
      = topSpace (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) := by
  have hd : 0 < d N := (m.tbl i).hd N
  have hP := m.isHermitian_P i N ω
  have hGm := m.isHermitian_tableGram i N ω
  have htoOp : ∀ x, toOp (m.P i N ω) x
      = (topSpace (m.tableGram i N ω) hGm).starProjection x := by
    intro x
    have h : Matrix.toEuclideanLin (m.P i N ω)
        = (topProj (m.tableGram i N ω) hGm).toLinearMap := by
      rw [MultiTableModel.P, LinearEquiv.apply_symm_apply]
    rw [show toOp (m.P i N ω) = Matrix.toEuclideanLin (m.P i N ω) from rfl, h]
    rfl
  have hu : vMax (m.tableGram i N ω) hGm ∈ topSpace (m.tableGram i N ω) hGm :=
    mem_topSpace_vMax hd _ _
  have hnu : ‖vMax (m.tableGram i N ω) hGm‖ = 1 := norm_vMax hd _ _
  have hge : (1 : ℝ) ≤ lamMax (m.P i N ω) hP := by
    have h := inner_toOp_self_le (m.P i N ω) hP (vMax (m.tableGram i N ω) hGm)
    rw [htoOp, Submodule.starProjection_eq_self_iff.mpr hu, real_inner_self_eq_norm_sq, hnu,
      one_pow] at h
    linarith
  have hle : lamMax (m.P i N ω) hP ≤ 1 := by
    have hmem := mem_topSpace_vMax hd (m.P i N ω) hP
    have hnorm := norm_vMax hd (m.P i N ω) hP
    have hray : toOp (m.P i N ω) (vMax (m.P i N ω) hP)
        = lamMax (m.P i N ω) hP • vMax (m.P i N ω) hP := toOp_of_mem_topSpace hmem
    have hcs : ⟪toOp (m.P i N ω) (vMax (m.P i N ω) hP), vMax (m.P i N ω) hP⟫_ℝ ≤ 1 := by
      rw [htoOp]
      calc ⟪(topSpace (m.tableGram i N ω) hGm).starProjection (vMax (m.P i N ω) hP),
            vMax (m.P i N ω) hP⟫_ℝ
          ≤ ‖(topSpace (m.tableGram i N ω) hGm).starProjection (vMax (m.P i N ω) hP)‖ *
              ‖vMax (m.P i N ω) hP‖ := real_inner_le_norm _ _
        _ ≤ ‖vMax (m.P i N ω) hP‖ * ‖vMax (m.P i N ω) hP‖ := by
            gcongr
            exact Submodule.norm_starProjection_apply_le _ _
        _ = 1 := by rw [hnorm]; ring
    rw [hray, real_inner_smul_left, real_inner_self_eq_norm_sq, hnorm, one_pow, mul_one] at hcs
    exact hcs
  have hlam : lamMax (m.P i N ω) hP = 1 := le_antisymm hle hge
  rw [topSpace_eq_eigenspace, hlam]
  ext x
  rw [Module.End.mem_eigenspace_iff, htoOp, one_smul]
  exact Submodule.starProjection_eq_self_iff


/-! #### The case `M = 1`

At `M = 1` the general route is unavailable: `abetaW_gap`, `topSimple_of_gap` and
`topSimple_whp_of_tendsto` all read the eigenvalue of index `1`. The paper has no `2 ≤ M`
condition and neither do the theorems below (global audit 2026-08-31, finding 9; cleanup wave
3 extends the `M = 1` split from `thm_svdstack_weighted_paper` to the other four forms). At
`M = 1` the statement is `prop:single_table`: the weighted Gram matrix is `w_0² P_0`, its top
eigenspace is that of `X_0ᵀ X_0`, and the limit is `β_0²` for every nonzero weight. -/

/-- At `M = 1` the weighted svdstack Gram matrix `∑_i w_i² P_i` is `w_0² P_0`, so its top
eigenspace is the top eigenspace of `X_0ᵀ X_0` (`topSpace_P_eq` and `topSpace_smul`). -/
private theorem topSpace_svdstackGramW_one (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (h1 : ¬ 1 < Fintype.card (Fin M)) (hwne : ∃ k, w k ≠ 0) (N : ℕ) (ω : Ω N) :
    topSpace (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω)
      = topSpace (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω) := by
  obtain ⟨k, hk⟩ := hwne
  have hM : M = 1 := by
    have h0 := NeZero.pos M
    simp only [Fintype.card_fin, not_lt] at h1
    omega
  subst hM
  have hk0 : k = 0 := Subsingleton.elim k 0
  subst hk0
  have hw2 : (0 : ℝ) < w 0 ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hEq : m.svdstackGramW w N ω = w 0 ^ 2 • m.P 0 N ω := by
    rw [MultiTableModel.svdstackGramW, Fin.sum_univ_one]
  have hsm : (w 0 ^ 2 • m.P 0 N ω).IsHermitian := isHermitian_smul (m.isHermitian_P 0 N ω) _
  rw [topSpace_congr hEq (m.isHermitian_svdstackGramW w N ω) hsm,
    topSpace_smul hw2 (m.P 0 N ω) (m.isHermitian_P 0 N ω) hsm, m.topSpace_P_eq 0 N ω]

/-- At `M = 1` the weighted svdstack performance is the single-table overlap. -/
private theorem svdstackPerfW_eq_overlap_one (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (h1 : ¬ 1 < Fintype.card (Fin M)) (hwne : ∃ k, w k ≠ 0) (N : ℕ) (ω : Ω N) :
    m.svdstackPerfW w N ω = overlap ((m.tbl 0).X N ω) ((m.tbl 0).v N) := by
  have hsp := m.topSpace_svdstackGramW_one w h1 hwne N ω
  have hproj : topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω)
        ((m.tbl 0).v N)
      = topProj (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω) ((m.tbl 0).v N) := by
    change (topSpace (m.svdstackGramW w N ω)
        (m.isHermitian_svdstackGramW w N ω)).starProjection ((m.tbl 0).v N)
      = (topSpace (m.tableGram 0 N ω)
        (m.isHermitian_tableGram 0 N ω)).starProjection ((m.tbl 0).v N)
    rw [hsp]
  change ‖topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω)
    ((m.tbl 0).v N)‖ ^ 2 = _
  rw [hproj]
  rfl

/-- At `M = 1` the weighted limit is `β_0²` for every nonzero weight vector: `A_{β,w}` is the
`1 × 1` matrix `w_0²`, its unit top eigenvector has `u_0² = 1`, and the numerator is
`(w_0 β_0 u_0)²`. The weight cancels. -/
private theorem svdstackLimitW_one (w β : Fin M → ℝ)
    (h1 : ¬ 1 < Fintype.card (Fin M)) (hwne : ∃ k, w k ≠ 0) :
    svdstackLimitW w β = β 0 ^ 2 := by
  obtain ⟨k, hk⟩ := hwne
  have hM : M = 1 := by
    have h0 := NeZero.pos M
    simp only [Fintype.card_fin, not_lt] at h1
    omega
  subst hM
  have hk0 : k = 0 := Subsingleton.elim k 0
  subst hk0
  have hw2 : (0 : ℝ) < w 0 ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hA00 : AbetaW w β 0 0 = w 0 ^ 2 := by
    rw [abetaW_apply, Abeta, Matrix.add_apply, Matrix.vecMulVec_apply,
      Matrix.diagonal_apply_eq]
    ring
  have hu1 : ‖vMax (AbetaW w β) (isHermitian_AbetaW w β)‖ = 1 :=
    norm_vMax one_pos _ (isHermitian_AbetaW w β)
  have hAu : toOp (AbetaW w β) (vMax (AbetaW w β) (isHermitian_AbetaW w β))
      = lamMax (AbetaW w β) (isHermitian_AbetaW w β) •
        vMax (AbetaW w β) (isHermitian_AbetaW w β) :=
    toOp_of_mem_topSpace (mem_topSpace_vMax one_pos _ (isHermitian_AbetaW w β))
  have hu0 : WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 *
      WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 = 1 := by
    have h := real_inner_self_eq_norm_sq (vMax (AbetaW w β) (isHermitian_AbetaW w β))
    rw [hu1, one_pow, real_inner_eq_dotProduct, dotProduct, Fin.sum_univ_one] at h
    exact h
  have hlam : lamMax (AbetaW w β) (isHermitian_AbetaW w β) = w 0 ^ 2 := by
    have h0 := congrFun (congrArg WithLp.ofLp hAu) 0
    have hl : toOp (AbetaW w β) (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0
        = AbetaW w β 0 0 * WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 := by
      change ∑ j, AbetaW w β 0 j *
        WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) j = _
      rw [Fin.sum_univ_one]
    rw [hl, hA00] at h0
    have hune : WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 ≠ 0 := by
      intro h
      rw [h, mul_zero] at hu0
      exact zero_ne_one hu0
    exact (mul_right_cancel₀ hune h0).symm
  have hdot : (w * β) ⬝ᵥ WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β))
      = w 0 * β 0 * WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 := by
    rw [dotProduct, Fin.sum_univ_one]
    rfl
  have hu02 : WithLp.ofLp (vMax (AbetaW w β) (isHermitian_AbetaW w β)) 0 ^ 2 = 1 := by
    rw [sq]; exact hu0
  rw [svdstackLimitW, hdot, hlam, div_eq_iff hw2.ne', mul_pow, mul_pow, hu02]
  ring

/-- `thm:svdstack_weighted` at `M = 1`, which is `prop:single_table`: the weighted svdstack
estimator is single-table SVD and its performance tends to `β_0² = betaSq θ_0 c_0`, the value
`svdstackLimitW w β` takes for every nonzero weight vector. -/
private theorem thm_svdstack_weighted_of_card_le_one (m : MultiTableModel μ M n d)
    (c β w : Fin M → ℝ) (h1 : ¬ 1 < Fintype.card (Fin M))
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hwne : ∃ k, w k ≠ 0)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β) := by
  rw [svdstackLimitW_one w β h1 hwne, hβdef 0, beta, Real.sq_sqrt (betaSq_nonneg _ _)]
  exact ((law 0).align).congr fun N =>
    Filter.Eventually.of_forall fun ω => (m.svdstackPerfW_eq_overlap_one w h1 hwne N ω).symm

/-- `thm:svdstack_weighted` for a general weight vector, the weighted `thm_svd_stack_general`.
`hsimple` is the paper's proviso "provided that `W A_β W` has a unique largest eigenvalue"
(`main_paper.tex:1284`); `abetaW_gap` and `topSimple_of_gap` discharge it at the optimal
weights. `hwne` gives `λ_max(A_{β,w}) ≥ w_k² > 0`, which the transfer between `Ṽ_wᵀ Ṽ_w` and
`Ṽ_w Ṽ_wᵀ` needs. The audit dropped the unused `hc : 0 < c i`.

Cleanup wave 3 dropped two more hypotheses. `hw : 0 ≤ w i` was unused: every quantity sees
`w` through `w_i²`, so the theorem holds for any real weight vector, while the paper takes
`w ∈ R^M_{≥0}` (`main_paper.tex:300`; mechanical audit 2026-08-31, finding 8). `h1 : 2 ≤ M`
is gone by the `M = 1` split above. -/
theorem thm_svdstack_weighted_general (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hwne : ∃ k, w k ≠ 0)
    (hsimple : TopSimple (AbetaW w β) (isHermitian_AbetaW w β))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β) := by
  by_cases h1 : 1 < Fintype.card (Fin M)
  case neg => exact m.thm_svdstack_weighted_of_card_le_one c β w h1 hβdef hwne law
  obtain ⟨k, hk⟩ := id hwne
  have hA : (AbetaW w β).IsHermitian := isHermitian_AbetaW w β
  have hwk : (0 : ℝ) < w k ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hApos : 0 < lamMax (AbetaW w β) hA := lt_of_lt_of_le hwk (sq_le_lamMax_AbetaW w β k)
  set bvec : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 (w * β) with hbvec
  have hsymm : ∀ (N : ℕ) (ω : Ω N), (m.VtW w N ω * (m.VtW w N ω)ᵀ).IsHermitian :=
    fun N ω => isHermitian_mul_transpose_self (m.VtW w N ω)
  -- (a) `λ_max(G_N) → λ_max(A_{β,w})`, by Weyl
  have hlam : TendstoInProb μ (fun N ω =>
      lamMax (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω)))
      (lamMax (AbetaW w β) hA) :=
    lamMax_tendstoInProb hA hsymm (m.gramEntriesW c β w hβdef law hI)
  -- (b) `‖P_top(G_N) (W β)‖² → ⟪v_max(A_{β,w}), W β⟫²`
  have hnumβ : TendstoInProb μ (fun N ω =>
      ‖topProj (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))
        bvec‖ ^ 2) (⟪vMax (AbetaW w β) hA, bvec⟫_ℝ ^ 2) :=
    topProj_overlap_tendsto hA hsymm (m.gramEntriesW c β w hβdef law hI) hsimple bvec
  -- (c) `Ṽ_w v → W β`
  have hy : TendstoInProbPi μ
      (fun N ω => (m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) (w * β) := by
    intro l
    have h := (m.align c law l).const_mul (w l)
    rw [← hβdef l] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    show w l * ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N)) l)
        = (m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N)) l
    rw [MultiTableModel.VtW, ← Matrix.mulVec_mulVec, Matrix.mulVec_diagonal]
  have hyβ : TendstoInProb μ (fun N ω =>
      ‖(WithLp.toLp 2 ((m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) :
          EuclideanSpace ℝ (Fin M)) - bvec‖ *
        (‖(WithLp.toLp 2 ((m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) :
          EuclideanSpace ℝ (Fin M))‖ + ‖bvec‖)) 0 := by
    have hφ : Continuous (fun y : Fin M → ℝ =>
        ‖(WithLp.toLp 2 y : EuclideanSpace ℝ (Fin M)) - bvec‖ *
          (‖(WithLp.toLp 2 y : EuclideanSpace ℝ (Fin M))‖ + ‖bvec‖)) := by
      fun_prop
    have h := TendstoInProbPi.comp_continuous hφ.continuousAt hy
    rw [hbvec, sub_self, norm_zero, zero_mul] at h
    exact h
  have hdiff : TendstoInProb μ (fun N ω =>
      ‖topProj (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))
        (WithLp.toLp 2 ((m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2 -
      ‖topProj (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))
        bvec‖ ^ 2) 0 := by
    refine TendstoInProb.of_le (fun N => ?_) hyβ
    filter_upwards with ω
    rw [sub_zero]
    exact abs_norm_sq_topProj_sub_le _ _ _ _
  have hnum : TendstoInProb μ (fun N ω =>
      ‖topProj (m.VtW w N ω * (m.VtW w N ω)ᵀ) (isHermitian_mul_transpose_self (m.VtW w N ω))
        (WithLp.toLp 2 ((m.VtW w N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2)
      (⟪vMax (AbetaW w β) hA, bvec⟫_ℝ ^ 2) := by
    have h := hdiff.add hnumβ
    rw [zero_add] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    ring
  -- (d) the quotient, and the identification of the limit
  have hg := hnum.div hlam hApos.ne'
  have hlimit : svdstackLimitW w β = ⟪vMax (AbetaW w β) hA, bvec⟫_ℝ ^ 2
      / lamMax (AbetaW w β) hA := by
    rw [svdstackLimitW, real_inner_eq_dotProduct, dotProduct_comm]
  rw [hlimit]
  -- (e) transfer from the good event
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hg
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.tendsto_measure_not_goodEventW c β w hβdef hwne h1 hsimple law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  change ¬ m.goodEventW w N ω
  intro hgood
  exact hω (m.svdstackPerfW_eq_of_goodEventW w N ω hgood)

omit [NeZero M] in
/-- `thm:svdstack_weighted` (`main_paper.tex:503`), Layer 1 form (`law` and `hI` in the
signature: the single-table law structure for each table, and the independence of the
tables; L4 of 2026-09-02). At the optimal weights `w_i⋆ = 1/√(1 - β_i²)` the performance of svdstack
tends to `S/(S+1)`. One detectable table is enough (`hthr`), unlike the unweighted
`thm_svd_stack_general`, which needs two. There is no `2 ≤ M` hypothesis: at `M = 1` the
statement is `prop:single_table` (cleanup wave 3). -/
theorem thm_svdstack_weighted [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β) := by
  obtain ⟨k, hk⟩ := hthr
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have hpos : ∀ l : Fin M, (0 : ℝ) < 1 - β l ^ 2 := by
    intro l
    obtain ⟨h0, h2⟩ := hβ01 l
    nlinarith
  have hsq : ∀ l : Fin M, optW β l ^ 2 * (1 - β l ^ 2) = 1 := by
    intro l
    rw [optW, div_pow, one_pow, Real.sq_sqrt (hpos l).le, one_div,
      inv_mul_cancel₀ (hpos l).ne']
  have hoptpos : ∀ l : Fin M, 0 < optW β l := fun l => by
    rw [optW]
    exact div_pos one_pos (Real.sqrt_pos.mpr (hpos l))
  by_cases h1 : 1 < Fintype.card (Fin M)
  case neg =>
    have h := m.thm_svdstack_weighted_of_card_le_one c β (optW β) h1 hβdef
      ⟨k, (hoptpos k).ne'⟩ law
    rwa [svdstackLimitW_optW β hβ01] at h
  have hA : (AbetaW (optW β) β).IsHermitian := isHermitian_AbetaW (optW β) β
  have hgap := abetaW_gap (optW β) β h1 (k := k) (fun l => by rw [hsq l, hsq k])
  have hγ : 0 < lamMax (AbetaW (optW β) β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    lt_of_lt_of_le (mul_pos (pow_pos (hoptpos k) 2) (pow_pos hk 2)) hgap
  have hsimple : TopSimple (AbetaW (optW β) β) hA := topSimple_of_gap hA h1 hγ
  have h := m.thm_svdstack_weighted_general c β (optW β) hβdef
    ⟨k, (hoptpos k).ne'⟩ hsimple law hI
  rwa [svdstackLimitW_optW β hβ01] at h

omit [NeZero M] in
/-- `thm:svdstack_weighted` at the paper's displayed weights `eq:svdstack.weight`
(`main_paper.tex:505`), the literal statement of the theorem. The weight of an undetectable
table is `0`, not `optW β i = 1`, and the limit is the same
(`main_paper.tex:1330`, `svdstackLimitW_paperW`).

The paper has no `2 ≤ M` condition and neither does this statement (global audit 2026-08-31,
finding 9). At `M = 1` the general route is unavailable, because `abetaW_gap` and
`topSimple_whp_of_tendsto` read the eigenvalue of index `1`. The `M = 1` case is
`prop:single_table` instead: the weighted svdstack Gram matrix is `w_0² P_0`, whose top
projector is the top projector of `X_0ᵀ X_0` (`topSpace_P_eq` and `topProj_smul`), and
`S/(S+1) = β_0²` there. -/
theorem thm_svdstack_weighted_paper [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (fun i => paperW (m.tbl i).θ (c i)) N ω)
      (svdstackLimitOpt β) := by
  by_cases h1 : 1 < Fintype.card (Fin M)
  case neg =>
    -- `M = 1`: the estimator is single-table SVD and the limit is `β_0²`
    obtain ⟨k, hk⟩ := hthr
    have hM : M = 1 := by
      have hk2 := k.isLt
      simp only [Fintype.card_fin, not_lt] at h1
      omega
    subst hM
    have hk0 : k = 0 := Subsingleton.elim k 0
    subst hk0
    set pw : Fin 1 → ℝ := fun i => paperW (m.tbl i).θ (c i) with hpwdef
    have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
      rw [hβdef i]
      exact beta_mem_Ico (hc i)
    have hposβ : (0 : ℝ) < 1 - β 0 ^ 2 := by
      obtain ⟨h0, h2⟩ := hβ01 0
      nlinarith
    have hkthr : c 0 < (m.tbl 0).θ ^ 4 := by
      by_contra h
      rw [hβdef 0, beta_eq_zero_of_not_thr h] at hk
      exact lt_irrefl 0 hk
    have hone : pw 0 ^ 2 * (1 - β 0 ^ 2) = 1 := by
      rw [hpwdef, hβdef 0]
      exact paperW_sq_mul_one_sub (m.tbl 0).hθ (hc 0) hkthr
    have hpw2 : (0 : ℝ) < pw 0 ^ 2 := by nlinarith
    -- the limit
    have hB : β 0 ^ 2 = betaSq (m.tbl 0).θ (c 0) := by
      rw [hβdef 0, beta, Real.sq_sqrt (betaSq_nonneg _ _)]
    have hden : β 0 ^ 2 / (1 - β 0 ^ 2) + 1 = 1 / (1 - β 0 ^ 2) := by
      field_simp
      ring
    have hlim : svdstackLimitOpt β = betaSq (m.tbl 0).θ (c 0) := by
      have hS : Sval β = β 0 ^ 2 / (1 - β 0 ^ 2) := by rw [Sval, Fin.sum_univ_one]
      rw [svdstackLimitOpt, hS, hden, ← hB, div_eq_iff (one_div_ne_zero hposβ.ne'),
        mul_one_div]
    -- the estimator
    have hpwne0 : pw 0 ≠ 0 := by
      intro h0
      rw [h0] at hpw2
      norm_num at hpw2
    rw [hlim]
    exact (law 0).align.congr fun N => Filter.Eventually.of_forall fun ω =>
      (m.svdstackPerfW_eq_overlap_one pw h1 ⟨0, hpwne0⟩ N ω).symm
  obtain ⟨k, hk⟩ := hthr
  set pw : Fin M → ℝ := fun i => paperW (m.tbl i).θ (c i) with hpwdef
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have hle : ∀ l, pw l ^ 2 * (1 - β l ^ 2) ≤ 1 := by
    intro l
    by_cases h : c l < (m.tbl l).θ ^ 4
    · rw [hpwdef, hβdef l]
      exact (paperW_sq_mul_one_sub (m.tbl l).hθ (hc l) h).le
    · rw [hpwdef]
      simp only [paperW, if_neg h]
      norm_num
  have hkthr : c k < (m.tbl k).θ ^ 4 := by
    by_contra h
    rw [hβdef k, beta_eq_zero_of_not_thr h] at hk
    exact lt_irrefl 0 hk
  have hkone : pw k ^ 2 * (1 - β k ^ 2) = 1 := by
    rw [hpwdef, hβdef k]
    exact paperW_sq_mul_one_sub (m.tbl k).hθ (hc k) hkthr
  have hposk : (0 : ℝ) < 1 - β k ^ 2 := by
    obtain ⟨h0, h2⟩ := hβ01 k
    nlinarith
  have hpwk : (0 : ℝ) < pw k ^ 2 := by nlinarith [hkone, hposk]
  have hA : (AbetaW pw β).IsHermitian := isHermitian_AbetaW pw β
  have hgap := abetaW_gap pw β h1 (k := k) (fun l => by rw [hkone]; exact hle l)
  have hγ : 0 < lamMax (AbetaW pw β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    lt_of_lt_of_le (mul_pos hpwk (pow_pos hk 2)) hgap
  have hsimple : TopSimple (AbetaW pw β) hA := topSimple_of_gap hA h1 hγ
  have hpwne : pw k ≠ 0 := by
    intro h0
    rw [h0] at hpwk
    norm_num at hpwk
  have h := m.thm_svdstack_weighted_general c β pw hβdef ⟨k, hpwne⟩ hsimple law hI
  rwa [svdstackLimitW_paperW (fun i => (m.tbl i).θ) c β (fun i => (m.tbl i).hθ) hc hβdef] at h

/-- `thm:svdstack_weighted` in the paper's form `|⟨v̂_svdstack(w⋆), v⟩|² → S/(S+1)`, for any
selection `vhat` of a unit top eigenvector of `∑_i (w_i⋆)² P_i`. Weighted form of
`thm_svd_stack_general_inner`; the transfer uses `topSimple_svdstackGramW`, and at `M = 1`
the almost sure `SingleTableLaw.topSimple` of table `0` (cleanup wave 3, no `2 ≤ M`). -/
theorem thm_svdstack_weighted_inner (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace (m.svdstackGramW (optW β) N ω)
      (m.isHermitian_svdstackGramW (optW β) N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2) (svdstackLimitOpt β) := by
  obtain ⟨k, hk⟩ := id hthr
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have hpos : ∀ l : Fin M, (0 : ℝ) < 1 - β l ^ 2 := by
    intro l
    obtain ⟨h0, h2⟩ := hβ01 l
    nlinarith
  have hsq : ∀ l : Fin M, optW β l ^ 2 * (1 - β l ^ 2) = 1 := by
    intro l
    rw [optW, div_pow, one_pow, Real.sq_sqrt (hpos l).le, one_div,
      inv_mul_cancel₀ (hpos l).ne']
  have hoptpos : ∀ l : Fin M, 0 < optW β l := fun l => by
    rw [optW]
    exact div_pos one_pos (Real.sqrt_pos.mpr (hpos l))
  by_cases h1 : 1 < Fintype.card (Fin M)
  case neg =>
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
      (m.thm_svdstack_weighted c β hc hβdef hthr law hI)
    have hzero : ∀ N, μ N {ω | ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2
        ≠ m.svdstackPerfW (optW β) N ω} = 0 := by
      intro N
      refine measure_mono_null (fun ω hω => ?_) (ae_iff.mp ((law 0).topSimple N))
      intro hsimpleN
      have hsp := m.topSpace_svdstackGramW_one (optW β) h1 ⟨k, (hoptpos k).ne'⟩ N ω
      have hs : TopSimple (m.svdstackGramW (optW β) N ω)
          (m.isHermitian_svdstackGramW (optW β) N ω) := by
        change Module.finrank ℝ (topSpace (m.svdstackGramW (optW β) N ω)
          (m.isHermitian_svdstackGramW (optW β) N ω)) = 1
        rw [hsp]
        exact hsimpleN
      exact hω (norm_topProj_sq_eq_inner_sq hs (hmem N ω) (hnorm N ω) _).symm
    simp only [hzero]
    exact tendsto_const_nhds
  have hA : (AbetaW (optW β) β).IsHermitian := isHermitian_AbetaW (optW β) β
  have hgap := abetaW_gap (optW β) β h1 (k := k) (fun l => by rw [hsq l, hsq k])
  have hγ : 0 < lamMax (AbetaW (optW β) β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    lt_of_lt_of_le (mul_pos (pow_pos (hoptpos k) 2) (pow_pos hk 2)) hgap
  have hsimple : TopSimple (AbetaW (optW β) β) hA := topSimple_of_gap hA h1 hγ
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.thm_svdstack_weighted c β hc hβdef hthr law hI)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topSimple_svdstackGramW c β (optW β) hβdef
      ⟨k, (hoptpos k).ne'⟩ h1 hsimple law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2 ≠ m.svdstackPerfW (optW β) N ω := hω
  change ¬ TopSimple (m.svdstackGramW (optW β) N ω) (m.isHermitian_svdstackGramW (optW β) N ω)
  intro hsimpleN
  exact hω' (norm_topProj_sq_eq_inner_sq hsimpleN (hmem N ω) (hnorm N ω) _).symm

/-! ### The case `β = 0` (`sec:svdstack_threshold`, case 2) -/

/-- Deterministic bound behind the weighted zero case. The top eigenvalue of `∑_i w_i² P_i` is
at least `w_k²` for every `k`, because `∑_i w_i² P_i ⪰ w_k² P_k` and `P_k` is a nonzero
projector. Writing `u = P_top(∑_i w_i² P_i) v`, the eigenvalue relation gives
`λ ‖u‖² = ⟪u, ∑_i w_i² P_i v⟫ ≤ ‖u‖ ∑_i w_i² ‖P_i v‖`, so `‖u‖ ≤ (∑_i w_i² ‖P_i v‖)/w_k²`.
The unweighted `svdstackPerf_le_sum` uses `λ ≥ 1`, which fails for small weights. -/
private theorem svdstackPerfW_le_sum (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {k : Fin M} (hk : w k ≠ 0) :
    m.svdstackPerfW w N ω ≤
      ((∑ i, w i ^ 2 * ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
        ((m.tbl 0).v N)‖) / w k ^ 2) ^ 2 := by
  set v := (m.tbl 0).v N with hvdef
  set G := m.svdstackGramW w N ω with hGdef
  set hGh := m.isHermitian_svdstackGramW w N ω with hGhdef
  set Q : Fin M → EuclideanSpace ℝ (Fin (d N)) →L[ℝ] EuclideanSpace ℝ (Fin (d N)) :=
    fun i => topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) with hQdef
  set S : ℝ := ∑ i, w i ^ 2 * ‖Q i v‖ with hSdef
  have hwk : (0 : ℝ) < w k ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hSnn : 0 ≤ S := Finset.sum_nonneg fun i _ => by positivity
  have hop : ∀ x, toOp G x = ∑ i, w i ^ 2 • Q i x := by
    have hmap : Matrix.toEuclideanLin G = ∑ i, w i ^ 2 • ((Q i : _ →L[ℝ] _) : _ →ₗ[ℝ] _) := by
      rw [hGdef]
      change Matrix.toEuclideanLin (∑ i, (w i ^ 2) • m.P i N ω) = _
      rw [map_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [map_smul]
      simp only [MultiTableModel.P, hQdef, LinearEquiv.apply_symm_apply]
    intro x
    rw [show toOp G x = Matrix.toEuclideanLin G x from rfl, hmap]
    simp
  set u := topProj G hGh v with hudef
  have humem : u ∈ topSpace G hGh := Submodule.starProjection_apply_mem _ v
  have hGu : toOp G u = lamMax G hGh • u := toOp_of_mem_topSpace humem
  have hsym : (toOp G).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hGh
  have hd0 : 0 < d N := (m.tbl 0).hd N
  -- `λ ≥ w_k²`
  have hone : w k ^ 2 ≤ lamMax G hGh := by
    set y := vMax (m.tableGram k N ω) (m.isHermitian_tableGram k N ω) with hydef
    have hyn : ‖y‖ = 1 := norm_vMax hd0 _ _
    have hQk : Q k y = y :=
      Submodule.starProjection_eq_self_iff.mpr (mem_topSpace_vMax hd0 _ _)
    have hsum : ⟪toOp G y, y⟫_ℝ = ∑ i, w i ^ 2 * ‖Q i y‖ ^ 2 := by
      rw [hop y, sum_inner]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [real_inner_smul_left, inner_topProj_self_eq_norm_sq]
    have h1 : w k ^ 2 ≤ ∑ i, w i ^ 2 * ‖Q i y‖ ^ 2 := by
      have h := Finset.single_le_sum (f := fun i => w i ^ 2 * ‖Q i y‖ ^ 2)
        (fun i _ => by positivity) (Finset.mem_univ k)
      rwa [hQk, hyn, one_pow, mul_one] at h
    have h2 := inner_toOp_self_le G hGh y
    rw [hsum, hyn] at h2
    nlinarith
  -- `λ ‖u‖² ≤ ‖u‖ S`
  have key : lamMax G hGh * ‖u‖ ^ 2 ≤ ‖u‖ * S := by
    have h1 : ⟪toOp G u, v⟫_ℝ = lamMax G hGh * ‖u‖ ^ 2 := by
      rw [hGu, real_inner_smul_left]
      congr 1
      exact inner_topProj_self_eq_norm_sq G hGh v
    have h2 : ⟪toOp G u, v⟫_ℝ = ∑ i, w i ^ 2 * ⟪u, Q i v⟫_ℝ := by
      rw [hsym u v, hop v, inner_sum]
      exact Finset.sum_congr rfl fun i _ => real_inner_smul_right _ _ _
    have h3 : ∑ i, w i ^ 2 * ⟪u, Q i v⟫_ℝ ≤ ∑ i, w i ^ 2 * (‖u‖ * ‖Q i v‖) :=
      Finset.sum_le_sum fun i _ =>
        mul_le_mul_of_nonneg_left (real_inner_le_norm _ _) (sq_nonneg _)
    rw [← h1, h2]
    calc ∑ i, w i ^ 2 * ⟪u, Q i v⟫_ℝ ≤ ∑ i, w i ^ 2 * (‖u‖ * ‖Q i v‖) := h3
      _ = ‖u‖ * S := by
          rw [hSdef, Finset.mul_sum]
          exact Finset.sum_congr rfl fun i _ => by ring
  change ‖u‖ ^ 2 ≤ (S / w k ^ 2) ^ 2
  rcases eq_or_lt_of_le (norm_nonneg u) with hu0 | hu0
  · rw [← hu0]
    simpa using sq_nonneg (S / w k ^ 2)
  · have hstep : w k ^ 2 * ‖u‖ ^ 2 ≤ lamMax G hGh * ‖u‖ ^ 2 :=
      mul_le_mul_of_nonneg_right hone (sq_nonneg _)
    have h5 : w k ^ 2 * ‖u‖ ≤ S := by
      have h : (w k ^ 2 * ‖u‖) * ‖u‖ ≤ S * ‖u‖ := by nlinarith [hstep, key]
      exact le_of_mul_le_mul_right h hu0
    have h6 : ‖u‖ ≤ S / w k ^ 2 := by
      rw [le_div_iff₀ hwk]
      linarith
    nlinarith [h6, norm_nonneg u]

-- `hI` is not used by the proof: the bound `svdstackPerfW ≤ (∑_i w_i² ‖P_i v‖ / w_k²)²` is
-- deterministic, and every term of the sum tends to `0` by `SingleTableLaw.align` alone. The
-- hypothesis stays in the statement to match `thm_svdstack_weighted_general`.
set_option linter.unusedVariables false in
/-- `sec:svdstack_threshold` case 2 for weighted svdstack: when every `β_i = 0` the weighted
performance tends to `0`, for every weight vector that is not identically zero. Weighted twin
of `thm_svd_stack_general_zero`; the audit found this paper claim had no Lean statement
(`notes/archive/audit_weighted_2026-08-30.md`, summary item 6). -/
theorem thm_svdstack_weighted_zero (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0) (hwne : ∃ k, w k ≠ 0)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) 0 := by
  obtain ⟨k, hk⟩ := hwne
  have hb : ∀ i, betaSq (m.tbl i).θ (c i) = 0 := by
    intro i
    have h : (0 : ℝ) = Real.sqrt (betaSq (m.tbl i).θ (c i)) := by
      rw [← hβ0 i, hβdef i, beta]
    have hle : betaSq (m.tbl i).θ (c i) ≤ 0 := Real.sqrt_eq_zero'.mp h.symm
    linarith [betaSq_nonneg (m.tbl i).θ (c i)]
  set F : (N : ℕ) → Ω N → Fin M → ℝ := fun N ω i =>
    w i ^ 2 * Real.sqrt (overlap ((m.tbl i).X N ω) ((m.tbl 0).v N)) with hFdef
  have hFconv : TendstoInProbPi μ F fun _ => (0 : ℝ) := by
    intro i
    have h2 : TendstoInProb μ (fun N ω => overlap ((m.tbl i).X N ω) ((m.tbl 0).v N))
        (betaSq (m.tbl i).θ (c i)) := by
      have h := (law i).align
      simpa [m.hv i 0] using h
    rw [hb i] at h2
    have h3 := h2.comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
    have h4 := h3.const_mul (w i ^ 2)
    simpa [hFdef] using h4
  have hlim : TendstoInProb μ (fun N ω => ((∑ i, F N ω i) / w k ^ 2) ^ 2) 0 := by
    have h := TendstoInProbPi.comp_continuous
      (φ := fun y : Fin M → ℝ => ((∑ i, y i) / w k ^ 2) ^ 2)
      (Continuous.continuousAt (by fun_prop)) hFconv
    simpa using h
  refine TendstoInProb.of_le (g := fun N ω => ((∑ i, F N ω i) / w k ^ 2) ^ 2)
    (fun N => ?_) hlim
  filter_upwards with ω
  have hle := m.svdstackPerfW_le_sum w N ω hk
  have hEq : ∀ i : Fin M, F N ω i = w i ^ 2 *
      ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) ((m.tbl 0).v N)‖ := by
    intro i
    simp only [hFdef]
    rw [show overlap ((m.tbl i).X N ω) ((m.tbl 0).v N)
        = ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) ((m.tbl 0).v N)‖ ^ 2
      from rfl, Real.sqrt_sq (norm_nonneg _)]
  have hnn : (0 : ℝ) ≤ m.svdstackPerfW w N ω := by
    change (0 : ℝ) ≤
      ‖topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω)
        ((m.tbl 0).v N)‖ ^ 2
    positivity
  rw [sub_zero, abs_of_nonneg hnn]
  calc m.svdstackPerfW w N ω
      ≤ ((∑ i, w i ^ 2 * ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
          ((m.tbl 0).v N)‖) / w k ^ 2) ^ 2 := hle
    _ = ((∑ i, F N ω i) / w k ^ 2) ^ 2 := by
        congr 2
        exact Finset.sum_congr rfl fun i _ => (hEq i).symm


omit [NeZero M] in
/-- `thm:svdstack_weighted` (`main_paper.tex:503`) with every random matrix theory hypothesis
discharged. The tables are Gaussian and independent (`hG`) and each is in the proportional
regime (`hreg`); no `SingleTableLaw` is assumed. Same proof shape as
`thm_svd_stack_general_gaussian`: apply `SpikedModel.singleTableLaw_of_gaussian` to each
table, whose Gaussian marginal comes from `gaussianNoise_of_joint`. -/
theorem thm_svdstack_weighted_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β) :=
  m.thm_svdstack_weighted c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

/-- `thm:svdstack_weighted` for a general weight vector, with every random matrix theory
hypothesis discharged: the tables are Gaussian and independent (`hG`) and each is in the
proportional regime (`hreg`); no `SingleTableLaw` is assumed. `hsimple` stays as a hypothesis:
it is a condition on the deterministic limit matrix `A_{β,w}`, not a law structure, so
`singleTableLaw_of_gaussian` cannot discharge it. The wrapper adds `hc : ∀ i, 0 < c i`, which
the base theorem does not need but which `singleTableLaw_of_gaussian` does. -/
theorem thm_svdstack_weighted_general_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hwne : ∃ k, w k ≠ 0)
    (hsimple : TopSimple (AbetaW w β) (isHermitian_AbetaW w β))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β) :=
  m.thm_svdstack_weighted_general c β w hβdef hwne hsimple
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

/-- `thm:svdstack_weighted` at the paper's displayed weights `eq:svdstack.weight`
(`main_paper.tex:505`), with every random matrix theory hypothesis discharged: the tables are
Gaussian and independent (`hG`) and each is in the proportional regime (`hreg`); no
`SingleTableLaw` is assumed. -/
theorem thm_svdstack_weighted_paper_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (fun i => paperW (m.tbl i).θ (c i)) N ω)
      (svdstackLimitOpt β) :=
  m.thm_svdstack_weighted_paper c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

/-- `thm:svdstack_weighted` in the paper's form `|⟨v̂_svdstack(w⋆), v⟩|² → S/(S+1)`, with every
random matrix theory hypothesis discharged: the tables are Gaussian and independent (`hG`) and
each is in the proportional regime (`hreg`); no `SingleTableLaw` is assumed. -/
theorem thm_svdstack_weighted_inner_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace (m.svdstackGramW (optW β) N ω)
      (m.isHermitian_svdstackGramW (optW β) N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2) (svdstackLimitOpt β) :=
  m.thm_svdstack_weighted_inner c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise vhat hmem hnorm

omit [NeZero M] in
/-- `sec:svdstack_threshold` case 2 for weighted svdstack, with every random matrix theory
hypothesis discharged: the tables are Gaussian and independent (`hG`) and each is in the
proportional regime (`hreg`); no `SingleTableLaw` is assumed. The wrapper adds
`hc : ∀ i, 0 < c i`, which the base theorem does not need but which
`singleTableLaw_of_gaussian` does. -/
theorem thm_svdstack_weighted_zero_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0) (hwne : ∃ k, w k ≠ 0)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) 0 :=
  m.thm_svdstack_weighted_zero c β w hβdef hβ0 hwne
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

omit [NeZero M] in
/-- Bundle of `thm:svdstack_weighted` and its optimality (`main_paper.tex:1296`), with every
random matrix theory hypothesis discharged. At the optimal weights `optW β` the performance of
svdstack tends to `svdstackLimitOpt β`. For every other weight vector `w` with a nonzero entry,
nonnegative entries, and a simple top eigenvalue of `A_{β,w}`, the performance tends to
`svdstackLimitW w β`, and that limit is at most `svdstackLimitOpt β`. -/
theorem thm_svdstack_weighted_gaussian_opt [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)
    ∧ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → (∀ i, 0 ≤ w i) →
        TopSimple (AbetaW w β) (isHermitian_AbetaW w β) →
        TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β)
        ∧ svdstackLimitW w β ≤ svdstackLimitOpt β := by
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  refine ⟨m.thm_svdstack_weighted_gaussian c β hc hβdef hthr hreg hG, ?_⟩
  intro w hwne hw hsimple
  exact ⟨m.thm_svdstack_weighted_general_gaussian c β w hc hβdef hwne hsimple hreg hG,
    svdstackLimitW_le_opt w β hw hβ01⟩

end SVDStackWeightedLemmas

end MultiTableModel

end StackedSVD
