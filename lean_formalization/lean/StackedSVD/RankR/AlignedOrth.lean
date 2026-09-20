/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.GeneralGaussian

/-!
# Item F2: the svdstack columns are asymptotically orthogonal

In the exactly aligned rank-`r` model at the optimal weights `W⋆ = optWG β`, the matrix
`V̂_svdstack(W⋆)ᵀ V` converges to a diagonal matrix: every off-diagonal entry squared tends to
`0` in probability. The paper states asymptotic orthogonality for the stacksvd columns only
(`main_paper.tex:926`) and says nothing for svdstack. This file closes that asymmetry.

## Route

Write `q_a` for the sorted eigenvector of `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ` at index `a`, `λ_a` for the sorted
eigenvalue there, and `y_j` for column `j` of `Ṽ_{W⋆} V`. At the canonical frame,

`((V̂_svdstack(W⋆)ᵀ V) a b)² = ⟪q_a, y_b⟫² / λ_a`

(`vhatSvdstackGW_entry_sq_gen`, the off-diagonal twin of `vhatSvdstackGW_entry_sq`). Fix
`k ≠ j`. Bessel at the two distinct sorted indices gives

`⟪q_k, y_j⟫² ≤ ‖y_j‖² - ⟪q_j, y_j⟫² ≤ ‖y_j‖² - (diagonal entry j)² λ_j`,

where the second step is `(a / λ) λ ≤ a`, which needs no positivity of `λ`. Three limits close
it: `‖y_j‖² → S_j` (item 2.2 of `notes/archive/rankr_D4_plan.md`, on the columns of `C = W⋆ B_R`),
the component clause `(diagonal entry j)² → S_j / (S_j + 1)`, and `λ_a → 1 + S_a` (item 2.4). The
numerator tends to `S_j - S_j = 0` and the denominator to `1 + S_k ≥ 1`.

So the argument is a squeeze against the component clause alone. It needs neither the
aggregate clause `thm_rank_r_svdstack_aggregate` nor the Frobenius identity
`frobSq_vhatSvdstackGW`, and hence no top-`r` eigengap of `A⋆`: that gap holds only when every
`S_j` is positive, and the component clause carries no such hypothesis (modeling choice 3 of
`notes/archive/F_batch_2026-09-05.md`).

## Statements

* `UnalignedModelR.thm_rank_r_svdstack_offdiag_of_sep`, one off-diagonal entry, Layer 1.
* `UnalignedModelR.thm_rank_r_svdstack_offdiag_of_model_gaussian`, the Gaussian facade with
  the hypotheses of `thm_rank_r_svdstack_component_of_model_gaussian`.
* `UnalignedModelR.thm_rank_r_svdstack_offdiag_frob_of_sep` and its Gaussian facade: the whole
  off-diagonal Frobenius mass tends to `0`.

No `sorry`, no `axiom`, no edit to any existing file.
-/

open MeasureTheory Filter Topology

open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. Two deterministic inputs -/

section Bessel

variable {p : ℕ}

/-- **Bessel at two distinct sorted eigenvector indices.** The sorted eigenvectors `vEig A hA a`
and `vEig A hA b` are two distinct members of one orthonormal eigenbasis, so the squared
inner products of any vector with them sum to at most the squared norm. -/
theorem inner_vEig_sq_add_le {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {a b : ℕ}
    (ha : a < Fintype.card (Fin p)) (hb : b < Fintype.card (Fin p)) (hab : a ≠ b)
    (y : EuclideanSpace ℝ (Fin p)) :
    ⟪vEig A hA a, y⟫_ℝ ^ 2 + ⟪vEig A hA b, y⟫_ℝ ^ 2 ≤ ‖y‖ ^ 2 := by
  classical
  have hne : (⟨a, ha⟩ : Fin (Fintype.card (Fin p))) ≠ ⟨b, hb⟩ := fun h =>
    hab (congrArg Fin.val h)
  have hbes := Orthonormal.sum_inner_products_le (𝕜 := ℝ)
    (s := ({⟨a, ha⟩, ⟨b, hb⟩} : Finset (Fin (Fintype.card (Fin p))))) y
    ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).orthonormal
  rw [Finset.sum_pair hne] at hbes
  rw [vEig_eq A hA ha, vEig_eq A hA hb]
  simpa [Real.norm_eq_abs, sq_abs] using hbes

end Bessel

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- **The frame quantity at a general entry.** The `(a, b)` entry of `V̂_svdstack(W)ᵀ V` is
`(√λ_a)^{-1} ⟪q_a, y_b⟫`, so its square is `⟪q_a, y_b⟫² / λ_a`. This is
`vhatSvdstackGW_entry_sq` (`RankR/AlignedMain.lean`) with the row index and the column index
kept apart; the proof is the same. The identity also holds at `λ_a = 0`, where both sides are
`0`. -/
theorem vhatSvdstackGW_entry_sq_gen (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) (a b : Fin r)
    (hlam : 0 ≤ lam a) :
    (((m.vhatSvdstackGW W N ω Q lam)ᵀ * m.V N) a b) ^ 2
      = ⟪frameCol Q a, m.VtVWGcol W N ω b⟫_ℝ ^ 2 / lam a := by
  have hinv : (Real.sqrt (lam a))⁻¹ ^ 2 = (lam a)⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hlam]
  rw [m.vhatSvdstackGW_transpose_mul_V W N ω Q lam, Matrix.diagonal_mul,
    ← m.frameCol_VtVWG W N ω b, inner_frameCol, mul_pow, hinv, inv_mul_eq_div]

end UnalignedModelR

/-! ### 2. The off-diagonal clause -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **Item F2: the svdstack columns are asymptotically orthogonal.** In the exactly aligned
rank-`r` model at `W⋆`, the off-diagonal entry `(k, j)` of `V̂_svdstack(W⋆)ᵀ V` squared tends to
`0` in probability, for every `k ≠ j`. The hypotheses are those of the component clause
`thm_rank_r_svdstack_component_eig_of_sep`; in particular there is no hypothesis `0 < S_j`. -/
theorem thm_rank_r_svdstack_offdiag_of_sep (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β) (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    {k j : Fin r} (hkj : k ≠ j) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) k j) ^ 2) 0 := by
  have hβ01 : ∀ i a, 0 ≤ β i a ∧ β i a < 1 := fun i a => by
    rw [hβdef i a]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i a, 0 ≤ β i a := fun i a => (hβ01 i a).1
  have h1 : ∀ i a, β i a < 1 := fun i a => (hβ01 i a).2
  have hsq : ∀ i a, β i a ^ 2 < 1 := fun i a => by nlinarith [h0 i a, h1 i a]
  have hrp : r ≤ Fintype.card (Fin (rtot (alignedRk M r))) := by
    simpa using le_rtot_alignedRk (M := M) (r := r) hM
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  have hkc : (k : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM k.isLt
  have hkjval : (k : ℕ) ≠ (j : ℕ) := fun h => hkj (Fin.val_injective h)
  have hgram := m.gramWG_tendsto_aligned c β hβdef hR law hG
  have hVtV := m.VtVWG_tendsto_aligned c β hβdef hR law
  -- the eigenvalue at a component index tends to `1 + S_a`
  have hlamOf : ∀ a : Fin r, TendstoInProb μ
      (fun N ω => eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (a : ℕ))
      (1 + Sagg β a) := by
    intro a
    have hpi : TendstoInProbPi μ
        (fun N ω => fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
          m.gramWG (optWG β) N ω t.1 t.2)
        (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
          ABlockW (optWG β) β (alignedR M r) t.1 t.2) := fun t => hgram t.1 t.2
    have h := TendstoInProbPi.comp_continuous
      (continuous_eigVal_symMat (p := rtot (alignedRk M r)) (a : ℕ)).continuousAt hpi
    rw [eigVal_symMat (isHermitian_ABlockW (optWG β) β (alignedR M r)),
      eigVal_ABlockW_aligned β h0 h1 hsep.1 hM a] at h
    refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
    exact eigVal_symMat (m.isHermitian_gramWG (optWG β) N ω) (a : ℕ)
  -- the squared norm of column `j` of `Ṽ_{W⋆} V` tends to `S_j`
  have hy : TendstoInProb μ (fun N ω => ‖m.VtVWGcol (optWG β) N ω j‖ ^ 2) (Sagg β j) := by
    have hpi : TendstoInProbPi μ
        (fun N ω => fun q : Fin (rtot (alignedRk M r)) => m.VtVWG (optWG β) N ω q j)
        (fun q : Fin (rtot (alignedRk M r)) => BBlockW (optWG β) β (alignedR M r) q j) :=
      fun q => hVtV q j
    have hcont : Continuous fun z : Fin (rtot (alignedRk M r)) → ℝ => ∑ q, z q ^ 2 :=
      continuous_finsetSum _ fun q _ => (continuous_apply q).pow 2
    have h := TendstoInProbPi.comp_continuous hcont.continuousAt hpi
    have hlim : ∑ q : Fin (rtot (alignedRk M r)),
        BBlockW (optWG β) β (alignedR M r) q j ^ 2 = Sagg β j := by
      have hc0 : ‖BBlockWcol (optWG β) β (alignedR M r) j‖ ^ 2 = Sagg β j :=
        norm_sq_BBlockWcol_aligned β hsq j
      rw [← hc0]
      exact (normSq_toLp fun q => BBlockW (optWG β) β (alignedR M r) q j).symm
    rw [hlim] at h
    refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
    exact (normSq_toLp fun q => m.VtVWG (optWG β) N ω q j).symm
  -- the entry identity at the canonical frame
  have hentry : ∀ (a b : Fin r) (N : ℕ) (ω : Ω N),
      (((m.vhatSvdstackGW (optWG β) N ω (topEigMat (m.isHermitian_gramWG (optWG β) N ω) hrp)
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω) hrp))ᵀ * m.V N) a b) ^ 2
        = ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (a : ℕ),
            m.VtVWGcol (optWG β) N ω b⟫_ℝ ^ 2
          / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (a : ℕ) := by
    intro a b N ω
    rw [m.vhatSvdstackGW_entry_sq_gen (optWG β) N ω _ _ a b
        (topEigVal_nonneg hrp (m.posSemidef_gramWG (optWG β) N ω) a),
      frameCol_topEigMat_eq_vEig, topEigVal_eq_eigVal]
  -- the eigenvalues of the random Gram matrix are nonnegative
  have hLnn : ∀ (a : Fin r) (N : ℕ) (ω : Ω N),
      0 ≤ eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (a : ℕ) := by
    intro a N ω
    rw [← topEigVal_eq_eigVal (m.isHermitian_gramWG (optWG β) N ω) hrp a]
    exact topEigVal_nonneg hrp (m.posSemidef_gramWG (optWG β) N ω) a
  -- the component clause, read on the quotient
  have hcomp : TendstoInProb μ
      (fun N ω => ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ),
          m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
        / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ))
      (Sagg β j / (Sagg β j + 1)) :=
    (m.thm_rank_r_svdstack_component_eig_of_sep c β hc hβdef hR hM hsep law hG j).congr
      fun N => Filter.Eventually.of_forall fun ω => hentry j j N ω
  -- the dominating sequence tends to `0`
  have hSj : (0 : ℝ) ≤ Sagg β j := Sagg_nonneg hsq j
  have hSk : (0 : ℝ) ≤ Sagg β k := Sagg_nonneg hsq k
  have hSjne : Sagg β j + 1 ≠ 0 := ne_of_gt (by linarith)
  have hSkne : 1 + Sagg β k ≠ 0 := ne_of_gt (by linarith)
  have hnum := hy.sub (hcomp.mul (hlamOf j))
  have hzero : Sagg β j - Sagg β j / (Sagg β j + 1) * (1 + Sagg β j) = 0 := by
    field_simp
    ring
  rw [hzero] at hnum
  have hquot := hnum.div (hlamOf k) hSkne
  rw [zero_div] at hquot
  -- the squeeze
  suffices h : TendstoInProb μ
      (fun N ω => ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ),
          m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
        / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ)) 0 from
    h.congr fun N => Filter.Eventually.of_forall fun ω => (hentry k j N ω).symm
  refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hquot
  have hLk := hLnn k N ω
  have hLj := hLnn j N ω
  have hbes := inner_vEig_sq_add_le (m.isHermitian_gramWG (optWG β) N ω) hkc hjc hkjval
    (m.VtVWGcol (optWG β) N ω j)
  have hdm : ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ),
        m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
      / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
      * eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
      ≤ ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ),
        m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2 := by
    rcases eq_or_lt_of_le hLj with hz | hz
    · rw [← hz, mul_zero]
      positivity
    · rw [div_mul_eq_mul_div, mul_div_assoc, div_self (ne_of_gt hz), mul_one]
  rw [sub_zero, abs_of_nonneg (div_nonneg (sq_nonneg _) hLk)]
  exact div_le_div_of_nonneg_right (by linarith) hLk

/-- **Item F2, the whole off-diagonal Frobenius mass.** The sum of the squares of every
off-diagonal entry of `V̂_svdstack(W⋆)ᵀ V` tends to `0` in probability, so the matrix converges
to a diagonal one. A finite sum of the entry clause. -/
theorem thm_rank_r_svdstack_offdiag_frob_of_sep (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β) (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => ∑ t : Fin r × Fin r, if t.1 = t.2 then 0 else
        (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) t.1 t.2) ^ 2) 0 := by
  classical
  have h : ∀ t : Fin r × Fin r, TendstoInProb μ
      (fun N ω => if t.1 = t.2 then 0 else
        (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) t.1 t.2) ^ 2) 0 := by
    intro t
    by_cases ht : t.1 = t.2
    · simp only [if_pos ht]
      exact TendstoInProb.const μ 0
    · simp only [if_neg ht]
      exact m.thm_rank_r_svdstack_offdiag_of_sep c β hc hβdef hR hM hsep law hG ht
  simpa using TendstoInProb.finsum h

end UnalignedModelR

/-! ### 3. The Gaussian facades -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **Item F2 under Gaussian noise.** The hypotheses are exactly those of
`thm_rank_r_svdstack_component_of_model_gaussian`: no `TableLawR` and no separation of `S`,
because the aligned model supplies both. -/
theorem thm_rank_r_svdstack_offdiag_of_model_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) {k j : Fin r} (hkj : k ≠ j) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) k j) ^ 2) 0 :=
  m.thm_rank_r_svdstack_offdiag_of_sep c β hc hβdef hR hM (m.saggSep_of_model c β hc hβdef)
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise hkj

/-- **Item F2 under Gaussian noise, the whole off-diagonal Frobenius mass.** -/
theorem thm_rank_r_svdstack_offdiag_frob_of_model_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => ∑ t : Fin r × Fin r, if t.1 = t.2 then 0 else
        (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) t.1 t.2) ^ 2) 0 :=
  m.thm_rank_r_svdstack_offdiag_frob_of_sep c β hc hβdef hR hM
    (m.saggSep_of_model c β hc hβdef) (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

end UnalignedModelR

end StackedSVD
