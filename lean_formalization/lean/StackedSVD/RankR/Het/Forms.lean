/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R2het
import StackedSVD.RMT.Het.R4het
import StackedSVD.RankR.Het.Split
import StackedSVD.RankR.RMT.Forms
import StackedSVD.Prob.GaussianMatrix

/-!
# Rank-`r` heteroscedastic resolvent forms (campaign E, stage E4)

The `r × r` resolvent-form limits of the weighted column split of `RankR/Het/Split.lean`.
Paper: `thm:rank_r_stacksvd` (`main_paper.tex:2337`, Appendix E, `sec:rank_r`); the
rank-one mirror is `thm:stacksvd_weighted` (`:463`). This file is the heteroscedastic twin of
`RankR/RMT/Forms.lean`. The rank-one counterpart in Lean is the model section of
`RMT/Het/R2het.lean` (`MultiTableModel.tendstoInProb_qform_u0Het` and its five siblings) and
the structure `ResolventLimitsHet` (`RMT/Het/R4het.lean:58`).

Objects.

1. `UtildeCol w N k`: column `k` of `Ũ = Σ^{1/2} A`, the deterministic signal column.
   `gHetCol w N ω k`: column `k` of `Σ^{1/2} E_stack V`, the weighted Gaussian column.
   `QmatHetR_col` splits column `k` of `Q` as their sum.
2. `ResolventLimitsHetR`: the edge and six `r × r` form limits on `W₀' = W0hetR` at real
   `z > b`. Diagonal entries tend to `Φ_k`, `Ψ`, `Φ_k'`, `Ψ'`; off-diagonal entries and both
   cross forms tend to `0`.
3. `ResolventLimitsHetR.cform_qcol`, `cform2_qcol`: the bilinear corollaries on the columns
   of `Q`, entrywise.
4. `tendstoInProb_dotProduct_QmatHetR_col`: the column Gram limit
   `Q_kᵀ Q_l → δ_kl (∑ w_i² θ_ik² + ∑ w_i² c_i)`.
5. `resolventLimitsHetR_of_gaussian`: the Gaussian discharge on the exactly aligned family
   (`rk = alignedRk M r`, `R_i = 1`), with the edge as a raw hypothesis. Stage E3 supplies
   the edge separately (`RankR/Het/Edge.lean`, not imported here).

Route of item 5. `exists_block_hasLaw_hetR` gives, for `d N = p N + r`, a block `B` with
`W₀' = Wsig τ d B` and the product law of `(Z_stack V, B)`. Column `k` of `Z_stack V` with
`B` has the pair law `HetR2.pairLaw`, through `measurePreserving_transpose` and the row
marginal of the transposed matrix; a unit mix `(x_k + x_l)/√2` of two columns has the same
law through `FormsR.measurePreserving_pairMix`. The general layer of `RMT/Het/R2het.lean`
(`tendstoInProb_qformC_fixed`, `qformC_gvec`, `cformC_gvec` and the three second-order
twins) gives every diagonal limit at complex `z`; polarization
(`FormsR.cformC_eq_polarizationC`) gives the off-diagonal entries from the mixed column;
`HetR2.tendstoInProb_cform_of_complex_scaled` moves each limit to the real axis.

Modeling choices.

1. The block `B` is chosen once per `N` by `Classical.choose` (`blockBR`), as the rank-one
   `MultiTableModel.blockB` does. The columns `x_k = (Z_stack V) e_k` are `xCol`.
2. The limits `Phi` and `Phi2` carry the column index `k` because the strength table
   `θ_ik` depends on it; `Psi` and `Psi2` do not.
3. The edge field is a hypothesis of `resolventLimitsHetR_of_gaussian`; the theorem does not
   import stage E3.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace UnalignedModelR

open HetR2 HetStein HetR1 MPhet R2 R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. The columns of the split -/

/-- Column `k` of `Ũ = Σ^{1/2} A`: the deterministic signal column of the weighted stack.
The rank-one mirror is `MultiTableModel.u0Het` (`RMT/Het/Split.lean:174`). -/
noncomputable def UtildeCol (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k : Fin r) : Fin (∑ i, n i N) → ℝ :=
  fun q => m.UtildeR w N q k

/-- Column `k` of `Σ^{1/2} E_stack V`: the weighted Gaussian column. The rank-one mirror is
`SigmaHalf *ᵥ eHet` (`RMT/Het/R2het.lean`, `SigmaHalf_mulVec_eHet`). -/
noncomputable def gHetCol (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (k : Fin r) : Fin (∑ i, n i N) → ℝ :=
  fun q => (m.SigmaHalfR w N * (m.stackEG N ω * m.V N)) q k

/-- Column `k` of the unscaled `Z_stack V`, the Gaussian vector `x_k = √d e_k` of the pair
law. -/
noncomputable def xCol (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) (k : Fin r) :
    Fin (∑ i, n i N) → ℝ :=
  fun q => (m.stackZG N ω * m.V N) q k

/-- Column `k` of `Q` is `ũ_k + g_k`. -/
theorem QmatHetR_col (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (k : Fin r) :
    (fun q => m.QmatHetR w N ω q k) = m.UtildeCol w N k + m.gHetCol w N ω k := by
  funext q
  simp only [QmatHetR, Matrix.add_apply, Pi.add_apply, UtildeCol, gHetCol]

/-- `g_k = gvec τ d x_k` with `τ_q = w_{blk q}`: the shape the general layer consumes. -/
theorem gHetCol_eq_gvec (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (k : Fin r) :
    m.gHetCol w N ω k = gvec (tauOf w (blkStack n N)) (d N) (m.xCol N ω k) := by
  funext q
  simp only [gHetCol, SigmaHalfR, m.stackE_eqG N ω, Matrix.smul_mul, Matrix.diagonal_mul,
    Matrix.smul_apply, smul_eq_mul, gvec, xCol, tauOf, blkStack]

/-- `Ũ` read at a row of block `i` on the exactly aligned family: `w_i θ_ik U_i(a, k)`. -/
theorem UtildeCol_apply_aligned (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (k : Fin r) (i : Fin M) (a : Fin (n i N)) :
    m.UtildeCol w N k (finSigmaFinEquiv ⟨i, a⟩)
      = w i * (m.thetaAligned i k * (m.tbl i).U N a k) := by
  have hσ : finSigmaFinEquiv.symm (finSigmaFinEquiv (⟨i, a⟩ : Σ i, Fin (n i N)))
      = (⟨i, a⟩ : Σ i, Fin (n i N)) :=
    Equiv.symm_apply_apply _ _
  simp only [UtildeCol, UtildeR, SigmaHalfR, Matrix.diagonal_mul, signalFactorG_apply', hR,
    Matrix.one_apply, thetaAligned]
  rw [hσ]
  simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, if_true]
  ring

/-- **Block bilinear sums of `Ũ`** on the exactly aligned family:
`∑_{q ∈ J_i} ũ_k(q) ũ_l(q) = δ_kl (θ_ik w_i)²`, from the orthonormal columns of `U_i`. -/
theorem sum_blockSet_UtildeCol_mul (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (i : Fin M) (k l : Fin r) :
    ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q * m.UtildeCol w N l q
      = if k = l then (m.thetaAligned i k * w i) ^ 2 else 0 := by
  rw [MultiTableModel.sum_blockSet_blkStack]
  simp only [m.UtildeCol_apply_aligned hR w N _ i]
  have hU : ∑ a : Fin (n i N), (m.tbl i).U N a k * (m.tbl i).U N a l
      = if k = l then 1 else 0 := by
    have h := congrFun (congrFun ((m.tbl i).hU N) k) l
    simpa [Matrix.mul_apply, Matrix.one_apply] using h
  have hexp : ∑ a : Fin (n i N), w i * (m.thetaAligned i k * (m.tbl i).U N a k)
      * (w i * (m.thetaAligned i l * (m.tbl i).U N a l))
      = (w i * m.thetaAligned i k) * (w i * m.thetaAligned i l)
        * ∑ a : Fin (n i N), (m.tbl i).U N a k * (m.tbl i).U N a l := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by ring
  rw [hexp, hU]
  split_ifs with hkl
  · subst hkl; ring
  · ring

/-- The block norms of `ũ_k`: `∑_{q ∈ J_i} ũ_k(q)² = (θ_ik w_i)²`. -/
theorem sum_blockSet_UtildeCol_sq (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (i : Fin M) (k : Fin r) :
    ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q ^ 2
      = (m.thetaAligned i k * w i) ^ 2 := by
  have h := m.sum_blockSet_UtildeCol_mul hR w N i k k
  simp only [if_true] at h
  rw [← h]
  exact Finset.sum_congr rfl fun q _ => sq _

/-- The block norms of the unit mix `(ũ_k + ũ_l)/√2`, `k ≠ l`. -/
theorem sum_blockSet_UtildeCol_mix_sq (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (i : Fin M) {k l : Fin r} (hkl : k ≠ l) :
    ∑ q ∈ blockSet (blkStack n N) i,
        ((Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l)) q ^ 2
      = ((m.thetaAligned i k * w i) ^ 2 + (m.thetaAligned i l * w i) ^ 2) / 2 := by
  have hkk := m.sum_blockSet_UtildeCol_mul hR w N i k k
  have hll := m.sum_blockSet_UtildeCol_mul hR w N i l l
  have hkl' := m.sum_blockSet_UtildeCol_mul hR w N i k l
  simp only [if_true, if_neg hkl] at hkk hll hkl'
  have h2 : ((Real.sqrt 2)⁻¹) ^ 2 = 1 / 2 := by
    rw [inv_pow, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2), one_div]
  have hexp : ∀ q, ((Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l)) q ^ 2
      = ((Real.sqrt 2)⁻¹) ^ 2 * (m.UtildeCol w N k q * m.UtildeCol w N k q
        + 2 * (m.UtildeCol w N k q * m.UtildeCol w N l q)
        + m.UtildeCol w N l q * m.UtildeCol w N l q) := by
    intro q
    simp only [Pi.smul_apply, Pi.add_apply, smul_eq_mul]
    ring
  simp only [hexp, ← Finset.mul_sum, Finset.sum_add_distrib, hkk, hll, hkl', h2]
  ring

/-- `ũ_k ⬝ᵥ ũ_l = δ_kl ∑_i (θ_ik w_i)²`. -/
theorem dotProduct_UtildeCol (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (k l : Fin r) :
    m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N l
      = if k = l then ∑ i, (m.thetaAligned i k * w i) ^ 2 else 0 := by
  have h : m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N l
      = ∑ i, ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q * m.UtildeCol w N l q := by
    rw [sum_blockSet_real]
    rfl
  rw [h]
  simp only [m.sum_blockSet_UtildeCol_mul hR w N _ k l]
  split_ifs <;> simp

/-- The squared norm of `ũ_k`. -/
theorem dotProduct_UtildeCol_self (m : UnalignedModelR μ M n d r (alignedRk M r))
    (hR : ∀ i, m.R i = 1) (w : Fin M → ℝ) (N : ℕ) (k : Fin r) :
    m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k = ∑ i, (m.thetaAligned i k * w i) ^ 2 := by
  have h := m.dotProduct_UtildeCol hR w N k k
  simpa using h

theorem sum_thetaW_sq_nonneg (m : UnalignedModelR μ M n d r (alignedRk M r)) (w : Fin M → ℝ)
    (k : Fin r) : 0 ≤ ∑ i, (m.thetaAligned i k * w i) ^ 2 :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

/-! ### 2. The hypothesis structure -/

/-- **The rank-`r` heteroscedastic (H1) and (H2)**, in the shape of `ResolventLimitsHet`
(`RMT/Het/R4het.lean:58`) with `r × r` entries. `W₀'` is `m.W0hetR`, the deterministic
columns `ũ_k = m.UtildeCol` play `v`, the Gaussian columns `g_k = m.gHetCol` play `g`.
Diagonal limits are `Phi k`, `Psi`, `Phi2 k`, `Psi2`; every off-diagonal and every cross
limit is `0`. On the Gaussian family the limits are `MPhet.Phihet (θ_·k) c w`,
`MPhet.Psihet c w` and the two derivatives (`resolventLimitsHetR_of_gaussian`). -/
structure ResolventLimitsHetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (b : ℝ)
    (Phi : Fin r → ℝ → ℝ) (Psi : ℝ → ℝ) (Phi2 : Fin r → ℝ → ℝ) (Psi2 : ℝ → ℝ) : Prop where
  /-- (H1), the edge bound. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1)
  /-- (H2), `ũ_kᵀ G₀' ũ_l → δ_kl Φ_k`. -/
  uu : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l))
    (if k = l then Phi k z else 0)
  /-- (H2), `g_kᵀ G₀' g_l → δ_kl Ψ`. -/
  ee : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l))
    (if k = l then Psi z else 0)
  /-- (H2), the cross form `ũ_kᵀ G₀' g_l → 0`. -/
  ue : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)) 0
  /-- (H2), `ũ_kᵀ G₀'² ũ_l → δ_kl Φ_k'`. -/
  uu2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l))
    (if k = l then Phi2 k z else 0)
  /-- (H2), `g_kᵀ G₀'² g_l → δ_kl Ψ'`. -/
  ee2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l))
    (if k = l then Psi2 z else 0)
  /-- (H2), the squared cross form `ũ_kᵀ G₀'² g_l → 0`. -/
  ue2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)) 0

/-! ### 3. The bilinear corollaries on the columns of `Q` -/

/-- **Entry `(k, l)` of `Qᵀ G₀'(z) Q`** tends to `δ_kl (Φ_k(z) + Ψ(z))`. The rank-`r`
homoscedastic mirror is `ResolventLimitsR.cform_qcol` (`RankR/RMT/Forms.lean:869`). -/
theorem ResolventLimitsHetR.cform_qcol {m : UnalignedModelR μ M n d r rk} {w : Fin M → ℝ}
    {b : ℝ} {Phi : Fin r → ℝ → ℝ} {Psi : ℝ → ℝ} {Phi2 : Fin r → ℝ → ℝ} {Psi2 : ℝ → ℝ}
    (h : m.ResolventLimitsHetR w b Phi Psi Phi2 Psi2) (k l : Fin r) {z : ℝ} (hz : b < z) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      (if k = l then Phi k z + Psi z else 0) := by
  have hfun : (fun N (ω : Ω N) => cform (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      = fun N ω =>
        cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l)
          + cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)
          + cform (m.W0hetR w N ω) z (m.UtildeCol w N l) (m.gHetCol w N ω k)
          + cform (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l) := by
    funext N ω
    rw [m.QmatHetR_col w N ω k, m.QmatHetR_col w N ω l]
    simp only [FormsR.cform_add_left, cform_add_right]
    rw [FormsR.cform_comm (m.isHermitian_W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.UtildeCol w N l)]
    ring
  have hcomb := (((h.uu k l z hz).add (h.ue k l z hz)).add (h.ue l k z hz)).add
    (h.ee k l z hz)
  rw [hfun]
  rcases eq_or_ne k l with rfl | hkl
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_true]
    ring
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_neg hkl]
    ring

/-- **Entry `(k, l)` of `Qᵀ G₀'(z)² Q`** tends to `δ_kl (Φ_k'(z) + Ψ'(z))`. -/
theorem ResolventLimitsHetR.cform2_qcol {m : UnalignedModelR μ M n d r rk} {w : Fin M → ℝ}
    {b : ℝ} {Phi : Fin r → ℝ → ℝ} {Psi : ℝ → ℝ} {Phi2 : Fin r → ℝ → ℝ} {Psi2 : ℝ → ℝ}
    (h : m.ResolventLimitsHetR w b Phi Psi Phi2 Psi2) (k l : Fin r) {z : ℝ} (hz : b < z) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      (if k = l then Phi2 k z + Psi2 z else 0) := by
  have hfun : (fun N (ω : Ω N) => cform2 (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      = fun N ω =>
        cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l)
          + cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)
          + cform2 (m.W0hetR w N ω) z (m.UtildeCol w N l) (m.gHetCol w N ω k)
          + cform2 (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l) := by
    funext N ω
    rw [m.QmatHetR_col w N ω k, m.QmatHetR_col w N ω l]
    simp only [FormsR.cform2_add_left, FormsR.cform2_add_right']
    rw [FormsR.cform2_comm (m.isHermitian_W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.UtildeCol w N l)]
    ring
  have hcomb := (((h.uu2 k l z hz).add (h.ue2 k l z hz)).add (h.ue2 l k z hz)).add
    (h.ee2 k l z hz)
  rw [hfun]
  rcases eq_or_ne k l with rfl | hkl
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_true]
    ring
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_neg hkl]
    ring

/-! ### 4. The block, the pair laws, and the regime parameters -/

section Block

variable [NeZero M] (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)

/-- The block `B` of `exists_block_hasLaw_hetR`, chosen once per `N`. The rank-one mirror is
`MultiTableModel.blockB` (`RMT/Het/R2het.lean:2004`). -/
noncomputable def blockBR (N : ℕ) : Ω N → Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ :=
  Classical.choose (m.exists_block_hasLaw_hetR hG N (hpd N))

theorem blockBR_gram (N : ℕ) (ω : Ω N) :
    m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ
      = ((d N : ℝ))⁻¹ • (m.blockBR hG hpd N ω * (m.blockBR hG hpd N ω)ᵀ) :=
  (Classical.choose_spec (m.exists_block_hasLaw_hetR hG N (hpd N))).1 ω

theorem hasLaw_ZV_blockBR (N : ℕ) :
    HasLaw (fun ω => (m.stackZG N ω * m.V N, m.blockBR hG hpd N ω))
      ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) (p N))) (μ N) :=
  (Classical.choose_spec (m.exists_block_hasLaw_hetR hG N (hpd N))).2

/-- `W₀' = Wsig τ d B` on the chosen block. The rank-one mirror is `W0het_eq_Wsig_blockB`. -/
theorem W0hetR_eq_Wsig_blockBR (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.W0hetR w N ω = Wsig (tauOf w (blkStack n N)) (d N) (m.blockBR hG hpd N ω) := by
  rw [m.W0hetR_eq_of_block w N ω _ (m.blockBR_gram hG hpd N ω)]
  rfl

/-- The pair `(x_k, B)`. -/
noncomputable def pairColR (N : ℕ) (k : Fin r) (ω : Ω N) : PairSpace (∑ i, n i N) (p N) :=
  (m.xCol N ω k, m.blockBR hG hpd N ω)

/-- The mixed pair `((x_k + s x_l)/√2, B)`. -/
noncomputable def pairMixR (N : ℕ) (k l : Fin r) (s : ℝ) (ω : Ω N) :
    PairSpace (∑ i, n i N) (p N) :=
  ((fun q => (Real.sqrt 2)⁻¹ * (m.xCol N ω k q + s * m.xCol N ω l q)), m.blockBR hG hpd N ω)

/-- **The pair law of one column.** `(x_k, B) ~ piGauss ⊗ gaussianMatrix`: transpose the
`n_tot × r` factor, read its row `k`. -/
theorem hasLaw_pairColR (N : ℕ) (k : Fin r) :
    HasLaw (m.pairColR hG hpd N k) (pairLaw (∑ i, n i N) (p N)) (μ N) := by
  have h1 : MeasurePreserving
      (Prod.map (Matrix.transpose : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ →
          Matrix (Fin r) (Fin (∑ i, n i N)) ℝ)
        (id : Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ → Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ))
      ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) (p N)))
      ((gaussianMatrix r (∑ i, n i N)).prod (gaussianMatrix (∑ i, n i N) (p N))) :=
    (measurePreserving_transpose (∑ i, n i N) r).prod (MeasurePreserving.id _)
  have h2 : MeasurePreserving
      (Prod.map (fun A : Matrix (Fin r) (Fin (∑ i, n i N)) ℝ => A k)
        (id : Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ → Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ))
      ((gaussianMatrix r (∑ i, n i N)).prod (gaussianMatrix (∑ i, n i N) (p N)))
      (pairLaw (∑ i, n i N) (p N)) :=
    (FormsR.measurePreserving_rowEval r (∑ i, n i N) k).prod (MeasurePreserving.id _)
  have h := (h2.comp h1).fun_comp_hasLaw (m.hasLaw_ZV_blockBR hG hpd N)
  exact h

/-- **The pair law of a unit mix of two columns**, `k ≠ l`, `s = ±1`: through
`FormsR.measurePreserving_pairMix` on the transposed pair. -/
theorem hasLaw_pairMixR (N : ℕ) {k l : Fin r} (hkl : k ≠ l) {s : ℝ} (hs : s = 1 ∨ s = -1) :
    HasLaw (m.pairMixR hG hpd N k l s) (pairLaw (∑ i, n i N) (p N)) (μ N) := by
  have h1 : MeasurePreserving
      (Prod.map (Matrix.transpose : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ →
          Matrix (Fin r) (Fin (∑ i, n i N)) ℝ)
        (Matrix.transpose : Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ →
          Matrix (Fin (p N)) (Fin (∑ i, n i N)) ℝ))
      ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) (p N)))
      ((gaussianMatrix r (∑ i, n i N)).prod (gaussianMatrix (p N) (∑ i, n i N))) :=
    (measurePreserving_transpose (∑ i, n i N) r).prod
      (measurePreserving_transpose (∑ i, n i N) (p N))
  have h2 := FormsR.measurePreserving_pairMix r (p N) (∑ i, n i N) hkl hs
  have h3 : MeasurePreserving
      (Prod.map (id : (Fin (∑ i, n i N) → ℝ) → (Fin (∑ i, n i N) → ℝ))
        (Matrix.transpose : Matrix (Fin (p N)) (Fin (∑ i, n i N)) ℝ →
          Matrix (Fin (∑ i, n i N)) (Fin (p N)) ℝ))
      (R2.noiseLaw (p N) (∑ i, n i N)) (pairLaw (∑ i, n i N) (p N)) :=
    (MeasurePreserving.id _).prod (measurePreserving_transpose (p N) (∑ i, n i N))
  have h := ((h3.comp h2).comp h1).fun_comp_hasLaw (m.hasLaw_ZV_blockBR hG hpd N)
  convert h using 1
  funext ω
  simp only [Function.comp_apply, Prod.map, id, Matrix.transpose_transpose,
    Matrix.transpose_apply, pairMixR, xCol]

/-- The second component of the pair is the block. -/
theorem hasLaw_blockBR (N : ℕ) :
    HasLaw (m.blockBR hG hpd N) (gaussianMatrix (∑ i, n i N) (p N)) (μ N) :=
  (MeasureTheory.measurePreserving_snd (μ := gaussianMatrix (∑ i, n i N) r)
    (ν := gaussianMatrix (∑ i, n i N) (p N))).fun_comp_hasLaw (m.hasLaw_ZV_blockBR hG hpd N)

end Block

/-! ### 5. The regime parameters of the block sequence -/

/-- `p N / d N → 1` when `d N = p N + r`. The rank-one mirror is
`MultiTableModel.tendsto_pred_div` (`RMT/Het/R2het.lean:2116`). -/
theorem tendsto_p_div {p : ℕ → ℕ} (hd : Tendsto d atTop atTop) (hpd : ∀ N, d N = p N + r) :
    Tendsto (fun N => (p N : ℝ) / d N) atTop (𝓝 1) := by
  have hdR : Tendsto (fun N => ((d N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have h1 : Tendsto (fun N => (1 : ℝ) - (r : ℝ) / (d N : ℝ)) atTop (𝓝 (1 - 0)) :=
    tendsto_const_nhds.sub (tendsto_const_nhds.div_atTop hdR)
  rw [sub_zero] at h1
  refine h1.congr' ?_
  filter_upwards [hd.eventually_gt_atTop 0] with N hN
  have hne : (d N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  have hp : (p N : ℝ) = d N - r := by
    rw [hpd N]
    push_cast
    ring
  rw [hp]
  field_simp

/-- The block proportions of the stack: `|J_i| / d N → c_i`. -/
theorem tendsto_cN_blkStackR (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (i : Fin M) :
    Tendsto (fun N => cN (fun N => blkStack n N) d N i) atTop (𝓝 (c i)) := by
  have h : (fun N => cN (fun N => blkStack n N) d N i) = fun N => (n i N : ℝ) / d N := by
    funext N
    simp only [cN, card_blockSet_blkStack]
  rw [h]
  exact (hreg i).2.2

/-! ### 6. Polarization with three limits -/

section Polarization

variable {dN : ℕ → ℕ} {z : ℂ}

/-- **The polarization step with three limits.** `a + b = √2 u`, `qformC a → La`,
`qformC b → Lb`, `qformC u → Lu` with `2 Lu = La + Lb` give `cformC a b → 0`. The one-limit
version is `FormsR.tendstoInProb_cformC_cross` (`RankR/RMT/Forms.lean:325`). -/
theorem tendstoInProb_cformC_cross3
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {La Lb Lu : ℂ} (hL : 2 * Lu = La + Lb)
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (af N ω) - La‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (bf N ω) - Lb‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (uf N ω) - Lu‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (Wf N ω) z (af N ω) (bf N ω)‖) 0 := by
  have h2C : ((Real.sqrt 2 ^ 2 : ℝ) : ℂ) = 2 := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
    norm_num
  refine TendstoInProb.of_le
    (g := fun N ω => ‖R4C.qformC (Wf N ω) z (uf N ω) - Lu‖
      + (‖R4C.qformC (Wf N ω) z (af N ω) - La‖ + ‖R4C.qformC (Wf N ω) z (bf N ω) - Lb‖) * 2⁻¹)
    (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    have hpol := FormsR.cformC_eq_polarizationC (hWf N ω) hz.ne' (af N ω) (bf N ω)
    have hkey : R4C.qformC (Wf N ω) z (af N ω + bf N ω)
        = 2 * R4C.qformC (Wf N ω) z (uf N ω) := by
      rw [hsum N ω, R2.qformC_smul, h2C]
    have heq : R4C.cformC (Wf N ω) z (af N ω) (bf N ω)
        = ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - Lu)
          - (R4C.qformC (Wf N ω) z (af N ω) - La)
          - (R4C.qformC (Wf N ω) z (bf N ω) - Lb)) / 2 := by
      rw [hpol, hkey]
      linear_combination (1 / 2 : ℂ) * hL
    rw [heq, norm_div, Complex.norm_two]
    have h4 : ‖(2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - Lu)‖
        = 2 * ‖R4C.qformC (Wf N ω) z (uf N ω) - Lu‖ := by
      rw [norm_mul, Complex.norm_two]
    have h5 := norm_sub_le ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - Lu)
        - (R4C.qformC (Wf N ω) z (af N ω) - La)) (R4C.qformC (Wf N ω) z (bf N ω) - Lb)
    have h6 := norm_sub_le ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - Lu))
      (R4C.qformC (Wf N ω) z (af N ω) - La)
    rw [h4] at h6
    linarith
  · have hlim := hqu.add ((hqa.add hqb).mul_const 2⁻¹)
    simpa using hlim

/-- The polarization step with three limits, second order. -/
theorem tendstoInProb_cform2C_cross3
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {La Lb Lu : ℂ} (hL : 2 * Lu = La + Lb)
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (af N ω) - La‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (bf N ω) - Lb‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (uf N ω) - Lu‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (Wf N ω) z (af N ω) (bf N ω)‖) 0 := by
  have h2C : ((Real.sqrt 2 ^ 2 : ℝ) : ℂ) = 2 := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
    norm_num
  refine TendstoInProb.of_le
    (g := fun N ω => ‖R4C.qform2C (Wf N ω) z (uf N ω) - Lu‖
      + (‖R4C.qform2C (Wf N ω) z (af N ω) - La‖ + ‖R4C.qform2C (Wf N ω) z (bf N ω) - Lb‖) * 2⁻¹)
    (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    have hpol := FormsR.cform2C_eq_polarizationC (hWf N ω) hz.ne' (af N ω) (bf N ω)
    have hkey : R4C.qform2C (Wf N ω) z (af N ω + bf N ω)
        = 2 * R4C.qform2C (Wf N ω) z (uf N ω) := by
      rw [hsum N ω, R2.qform2C_smul, h2C]
    have heq : R4C.cform2C (Wf N ω) z (af N ω) (bf N ω)
        = ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - Lu)
          - (R4C.qform2C (Wf N ω) z (af N ω) - La)
          - (R4C.qform2C (Wf N ω) z (bf N ω) - Lb)) / 2 := by
      rw [hpol, hkey]
      linear_combination (1 / 2 : ℂ) * hL
    rw [heq, norm_div, Complex.norm_two]
    have h4 : ‖(2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - Lu)‖
        = 2 * ‖R4C.qform2C (Wf N ω) z (uf N ω) - Lu‖ := by
      rw [norm_mul, Complex.norm_two]
    have h5 := norm_sub_le ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - Lu)
        - (R4C.qform2C (Wf N ω) z (af N ω) - La)) (R4C.qform2C (Wf N ω) z (bf N ω) - Lb)
    have h6 := norm_sub_le ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - Lu))
      (R4C.qform2C (Wf N ω) z (af N ω) - La)
    rw [h4] at h6
    linarith
  · have hlim := hqu.add ((hqa.add hqb).mul_const 2⁻¹)
    simpa using hlim

end Polarization

/-- `x + y = √2 ((√2)⁻¹ (x + y))`, the unit mix of two vectors. -/
theorem add_eq_sqrt_two_smul_mix {D : ℕ} (x y : Fin D → ℝ) :
    x + y = Real.sqrt 2 • ((Real.sqrt 2)⁻¹ • (x + y)) := by
  rw [smul_smul, mul_inv_cancel₀ (Real.sqrt_ne_zero'.mpr two_pos), one_smul]

/-- `gvec τ d x + gvec τ d y = √2 gvec τ d ((√2)⁻¹ (x + 1 y))`, the shape of `pairMixR`. -/
theorem gvec_add_eq_sqrt_two_smul {D : ℕ} (τ : Fin D → ℝ) (dd : ℕ) (x y : Fin D → ℝ) :
    gvec τ dd x + gvec τ dd y
      = Real.sqrt 2 • gvec τ dd (fun q => (Real.sqrt 2)⁻¹ * (x q + 1 * y q)) := by
  funext q
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, gvec, one_mul]
  have h2 : Real.sqrt 2 * (Real.sqrt 2)⁻¹ = 1 := mul_inv_cancel₀ (Real.sqrt_ne_zero'.mpr two_pos)
  linear_combination (-(τ q * (Real.sqrt (dd : ℕ))⁻¹ * (x q + y q))) * h2

/-! ### 7. The forms at complex `z` on the Gaussian family -/

section ComplexForms

variable [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
  (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {z : ℂ} (hz : 0 < z.im)

include hc hw hR hreg hG hpd hz in
/-- `ũ_kᵀ G₀'(z) ũ_k → Φ_k(z)` at complex `z`. The rank-one mirror is
`MultiTableModel.tendstoInProb_qformC_u0Het`. -/
theorem tendstoInProb_qformC_UtildeCol (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0hetR w N ω) z (m.UtildeCol w N k)
      - MultiTableModel.PhiC (fun i => m.thetaAligned i k) c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
    (fun N => m.UtildeCol w N k) (fun i => (m.thetaAligned i k * w i) ^ 2)
    (fun N i => m.sum_blockSet_UtildeCol_sq hR w N i k)

include hc hw hR hreg hG hpd hz in
/-- `ũ_kᵀ G₀'(z)² ũ_k → Φ_k'(z)` at complex `z`. -/
theorem tendstoInProb_qform2C_UtildeCol (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z (m.UtildeCol w N k)
      - MultiTableModel.PhiDerivC (fun i => m.thetaAligned i k) c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
    (fun N => m.UtildeCol w N k) (fun i => (m.thetaAligned i k * w i) ^ 2)
    (fun N i => m.sum_blockSet_UtildeCol_sq hR w N i k)

/-- The complex limit of the unit mix `(ũ_k + ũ_l)/√2`: the mean of `Φ_k` and `Φ_l`. -/
noncomputable def PhiMixC (m : UnalignedModelR μ M n d r (alignedRk M r)) (c w : Fin M → ℝ)
    (k l : Fin r) (z : ℂ) : ℂ :=
  ∑ i, ((((m.thetaAligned i k * w i) ^ 2 + (m.thetaAligned i l * w i) ^ 2) / 2 : ℝ) : ℂ)
    * gC w i z (sGlob c w z)

omit [NeZero M] in
theorem two_mul_PhiMixC (m : UnalignedModelR μ M n d r (alignedRk M r)) (c w : Fin M → ℝ)
    (k l : Fin r) (z : ℂ) :
    2 * m.PhiMixC c w k l z
      = MultiTableModel.PhiC (fun i => m.thetaAligned i k) c w z
        + MultiTableModel.PhiC (fun i => m.thetaAligned i l) c w z := by
  simp only [PhiMixC, MultiTableModel.PhiC, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  push_cast
  ring

include hc hw hR hreg hG hpd hz in
/-- The unit mix `(ũ_k + ũ_l)/√2`, `k ≠ l`: `qformC → (Φ_k + Φ_l)/2` at complex `z`. -/
theorem tendstoInProb_qformC_UtildeMix {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0hetR w N ω) z
      ((Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      - m.PhiMixC c w k l z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
    (fun N => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
    (fun i => ((m.thetaAligned i k * w i) ^ 2 + (m.thetaAligned i l * w i) ^ 2) / 2)
    (fun N i => m.sum_blockSet_UtildeCol_mix_sq hR w N i hkl)

/-- The second-order complex limit of the unit mix. -/
noncomputable def PhiMixDerivC (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c w : Fin M → ℝ) (k l : Fin r) (z : ℂ) : ℂ :=
  ∑ i, ((((m.thetaAligned i k * w i) ^ 2 + (m.thetaAligned i l * w i) ^ 2) / 2 : ℝ) : ℂ)
    * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹

omit [NeZero M] in
theorem two_mul_PhiMixDerivC (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c w : Fin M → ℝ) (k l : Fin r) (z : ℂ) :
    2 * m.PhiMixDerivC c w k l z
      = MultiTableModel.PhiDerivC (fun i => m.thetaAligned i k) c w z
        + MultiTableModel.PhiDerivC (fun i => m.thetaAligned i l) c w z := by
  simp only [PhiMixDerivC, MultiTableModel.PhiDerivC, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  push_cast
  ring

include hc hw hR hreg hG hpd hz in
/-- The unit mix `(ũ_k + ũ_l)/√2`, `k ≠ l`: `qform2C → (Φ_k' + Φ_l')/2` at complex `z`. -/
theorem tendstoInProb_qform2C_UtildeMix {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z
      ((Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      - m.PhiMixDerivC c w k l z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
    (fun N => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
    (fun i => ((m.thetaAligned i k * w i) ^ 2 + (m.thetaAligned i l * w i) ^ 2) / 2)
    (fun N i => m.sum_blockSet_UtildeCol_mix_sq hR w N i hkl)

include hc hw hR hreg hG hpd hz in
/-- **Off-diagonal `uu` at complex `z`**: `ũ_kᵀ G₀'(z) ũ_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cformC_UtildeCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.UtildeCol w N l)‖) 0 :=
  tendstoInProb_cformC_cross3 (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N _ => m.UtildeCol w N k) (bf := fun N _ => m.UtildeCol w N l)
    (uf := fun N _ => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
    (m.two_mul_PhiMixC c w k l z) (fun _ _ => add_eq_sqrt_two_smul_mix _ _)
    (m.tendstoInProb_qformC_UtildeCol w c hc hw hR hreg hG hpd hz k)
    (m.tendstoInProb_qformC_UtildeCol w c hc hw hR hreg hG hpd hz l)
    (m.tendstoInProb_qformC_UtildeMix w c hc hw hR hreg hG hpd hz hkl)

include hc hw hR hreg hG hpd hz in
/-- **Off-diagonal `uu2` at complex `z`**: `ũ_kᵀ G₀'(z)² ũ_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2C_UtildeCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.UtildeCol w N l)‖) 0 :=
  tendstoInProb_cform2C_cross3 (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N _ => m.UtildeCol w N k) (bf := fun N _ => m.UtildeCol w N l)
    (uf := fun N _ => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
    (m.two_mul_PhiMixDerivC c w k l z) (fun _ _ => add_eq_sqrt_two_smul_mix _ _)
    (m.tendstoInProb_qform2C_UtildeCol w c hc hw hR hreg hG hpd hz k)
    (m.tendstoInProb_qform2C_UtildeCol w c hc hw hR hreg hG hpd hz l)
    (m.tendstoInProb_qform2C_UtildeMix w c hc hw hR hreg hG hpd hz hkl)

include hc hw hreg hG hpd hz in
/-- `g_kᵀ G₀'(z) g_k → Ψ(z)` at complex `z`. The rank-one mirror is
`MultiTableModel.tendstoInProb_qformC_eHet`. -/
theorem tendstoInProb_qformC_gHetCol (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      - MultiTableModel.PsiC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_qformC_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

include hc hw hreg hG hpd hz in
/-- `g_kᵀ G₀'(z)² g_k → Ψ'(z)` at complex `z`. -/
theorem tendstoInProb_qform2C_gHetCol (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      - MultiTableModel.PsiDerivC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_qform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

include hc hw hreg hG hpd hz in
/-- The unit mix `(g_k + g_l)/√2`, `k ≠ l`: `qformC → Ψ` at complex `z`. -/
theorem tendstoInProb_qformC_gHetMix {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0hetR w N ω) z
      (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
      - MultiTableModel.PsiC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qformC_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairMixR hG hpd · k l 1)
    (fun N => m.hasLaw_pairMixR hG hpd N hkl (Or.inl rfl))

include hc hw hreg hG hpd hz in
/-- The unit mix `(g_k + g_l)/√2`, `k ≠ l`: `qform2C → Ψ'` at complex `z`. -/
theorem tendstoInProb_qform2C_gHetMix {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z
      (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
      - MultiTableModel.PsiDerivC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairMixR hG hpd · k l 1)
    (fun N => m.hasLaw_pairMixR hG hpd N hkl (Or.inl rfl))

/-- `g_k + g_l = √2 gvec τ d (x_k + x_l)/√2`, the sum in the shape of `pairMixR`. -/
theorem gHetCol_add_eq (N : ℕ) (ω : Ω N) (k l : Fin r) :
    m.gHetCol w N ω k + m.gHetCol w N ω l
      = Real.sqrt 2 • gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1 := by
  rw [m.gHetCol_eq_gvec w, m.gHetCol_eq_gvec w]
  exact gvec_add_eq_sqrt_two_smul _ _ _ _

include hc hw hreg hG hpd hz in
/-- **Off-diagonal `ee` at complex `z`**: `g_kᵀ G₀'(z) g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cformC_gHetCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)‖) 0 :=
  FormsR.tendstoInProb_cformC_cross (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N ω => m.gHetCol w N ω k) (bf := fun N ω => m.gHetCol w N ω l)
    (uf := fun N ω => gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
    (fun N ω => m.gHetCol_add_eq w hG hpd N ω k l)
    (m.tendstoInProb_qformC_gHetCol w c hc hw hreg hG hpd hz k)
    (m.tendstoInProb_qformC_gHetCol w c hc hw hreg hG hpd hz l)
    (m.tendstoInProb_qformC_gHetMix w c hc hw hreg hG hpd hz hkl)

include hc hw hreg hG hpd hz in
/-- **Off-diagonal `ee2` at complex `z`**: `g_kᵀ G₀'(z)² g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2C_gHetCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)‖) 0 :=
  FormsR.tendstoInProb_cform2C_cross (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N ω => m.gHetCol w N ω k) (bf := fun N ω => m.gHetCol w N ω l)
    (uf := fun N ω => gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
    (fun N ω => m.gHetCol_add_eq w hG hpd N ω k l)
    (m.tendstoInProb_qform2C_gHetCol w c hc hw hreg hG hpd hz k)
    (m.tendstoInProb_qform2C_gHetCol w c hc hw hreg hG hpd hz l)
    (m.tendstoInProb_qform2C_gHetMix w c hc hw hreg hG hpd hz hkl)

include hc hw hR hreg hG hpd hz in
/-- **`ue` at complex `z`**: `ũ_kᵀ G₀'(z) g_l → 0`. The rank-one mirror is
`MultiTableModel.tendstoInProb_cformC_u0Het_eHet`. -/
theorem tendstoInProb_cformC_UtildeCol_gHetCol (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.gHetCol w N ω l)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_cformC_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · l) (fun N => m.hasLaw_pairColR hG hpd N l)
    (fun N => m.UtildeCol w N k) (Ky := ∑ i, (m.thetaAligned i k * w i) ^ 2)
    (fun N => (m.dotProduct_UtildeCol_self hR w N k).le)

include hc hw hR hreg hG hpd hz in
/-- **`ue2` at complex `z`**: `ũ_kᵀ G₀'(z)² g_l → 0`. -/
theorem tendstoInProb_cform2C_UtildeCol_gHetCol (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.gHetCol w N ω l)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_cform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · l) (fun N => m.hasLaw_pairColR hG hpd N l)
    (fun N => m.UtildeCol w N k) (Ky := ∑ i, (m.thetaAligned i k * w i) ^ 2)
    (fun N => (m.dotProduct_UtildeCol_self hR w N k).le)

include hw hreg hG hpd in
/-- `‖g_k‖² → ∑ w_i² c_i` in probability. The rank-one mirror is
`MultiTableModel.tendstoInProb_SigmaHalf_eHet_norm`. -/
theorem tendstoInProb_gHetCol_norm (k : Fin r) :
    TendstoInProb μ (fun N ω => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k)
      (∑ i, w i ^ 2 * c i) := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.gHetCol_eq_gvec w]
  exact tendstoInProb_dotProduct_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hw hd
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

include hw hreg hG hpd in
/-- `‖(g_k + g_l)/√2‖² → ∑ w_i² c_i` in probability, `k ≠ l`. -/
theorem tendstoInProb_gHetMix_norm {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω =>
      gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1
        ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
      (∑ i, w i ^ 2 * c i) := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  exact tendstoInProb_dotProduct_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hw hd
    (m.tendsto_cN_blkStackR hreg) (m.pairMixR hG hpd · k l 1)
    (fun N => m.hasLaw_pairMixR hG hpd N hkl (Or.inl rfl))

end ComplexForms

/-! ### 8. Measurability and the norm events -/

section Events

/-- Entry `q` of `g_k` is measurable. The rank-one mirror is
`MultiTableModel.measurable_eHet_apply`. -/
theorem measurable_gHetCol_apply (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k : Fin r) (q : Fin (∑ i, n i N)) : Measurable fun ω => m.gHetCol w N ω k q := by
  have hE : ∀ (j : Fin (∑ i, n i N)) (a : Fin (d N)), Measurable fun ω => m.stackEG N ω j a :=
    fun j a => (measurable_pi_apply a).comp ((measurable_pi_apply j).comp (m.measurable_stackEG N))
  have h : (fun ω => m.gHetCol w N ω k q)
      = fun ω => ∑ j, m.SigmaHalfR w N q j * ∑ a, m.stackEG N ω j a * m.V N a k := rfl
  rw [h]
  exact Finset.measurable_sum _ fun j _ =>
    (Finset.measurable_sum _ fun a _ => (hE j a).mul_const _).const_mul _

/-- `g_k ⬝ᵥ g_l` is measurable. The rank-one mirror is
`MultiTableModel.measurable_SigmaHalf_eHet_dot`. -/
theorem measurable_gHetCol_dot (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k l : Fin r) : Measurable fun ω => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l :=
  Finset.measurable_sum _ fun q _ =>
    (m.measurable_gHetCol_apply w N k q).mul (m.measurable_gHetCol_apply w N l q)

variable [∀ N, IsProbabilityMeasure (μ N)] [NeZero M]
  (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
  (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)

omit [NeZero M] in
include hR in
/-- The norm event of `ũ_k` is the whole space. -/
theorem tendsto_measure_UtildeCol_norm_le (k : Fin r) :
    Tendsto (fun N => μ N {_ω : Ω N | m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k
      ≤ ∑ i, (m.thetaAligned i k * w i) ^ 2}) atTop (𝓝 1) := by
  have : ∀ N, {_ω : Ω N | m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k
      ≤ ∑ i, (m.thetaAligned i k * w i) ^ 2} = Set.univ :=
    fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_UtildeCol_self hR w N k).le
  simp only [this, measure_univ]
  exact tendsto_const_nhds

include hw hreg hG hpd in
/-- The norm event of `g_k`: `‖g_k‖² ≤ ∑ w_i² c_i + 1` with probability `→ 1`. The rank-one
mirror is `MultiTableModel.tendsto_measure_SigmaHalf_eHet_norm_le`. -/
theorem tendsto_measure_gHetCol_norm_le (k : Fin r) :
    Tendsto (fun N => μ N {ω | m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k
      ≤ ∑ i, w i ^ 2 * c i + 1}) atTop (𝓝 1) := by
  have h := m.tendstoInProb_gHetCol_norm w c hw hreg hG hpd k 1 one_pos
  refine tendsto_measure_one_of_bad (s := fun N => {ω | 1 ≤
    |m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k - ∑ i, w i ^ 2 * c i|})
    (fun N _ hω => ?_) h
  simp only [Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω ⊢
  rw [le_abs]
  left
  linarith

omit [NeZero M] in
include hc in
theorem sum_wSqC_add_one_nonneg : 0 ≤ ∑ i, w i ^ 2 * c i + 1 := by
  have : 0 ≤ ∑ i, w i ^ 2 * c i := Finset.sum_nonneg fun i _ => by
    have := (hc i).le; positivity
  linarith

end Events

/-! ### 9. The forms at real `x > b` (the fields of `ResolventLimitsHetR`) -/

section RealForms

variable [∀ N, IsProbabilityMeasure (μ N)] [NeZero M]
  (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
  (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {b : ℝ} (hb : bHet c w ≤ b)
  (hedge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1))
  {x : ℝ} (hx : b < x)

include hc hw hR hreg hG hpd hb hedge hx in
/-- **Field `uu`, diagonal, at real `x`**: `ũ_kᵀ G₀'(x) ũ_k → Φ_k(x)`. The rank-one mirror
is `MultiTableModel.tendstoInProb_qform_u0Het`. -/
theorem tendstoInProb_cform_UtildeCol_self (k : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.UtildeCol w N k)) (Phihet (fun i => m.thetaAligned i k) c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have h := tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (m.sum_thetaW_sq_nonneg w k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (MultiTableModel.PhiC (fun i => m.thetaAligned i k) c w)
    (Phihet (fun i => m.thetaAligned i k) c w x)
    (fun z hz => m.tendstoInProb_qformC_UtildeCol w c hc hw hR hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PhiC hc hw hbx)
  exact h

include hc hw hR hreg hG hpd hb hedge hx in
/-- **Field `uu2`, diagonal, at real `x`**: `ũ_kᵀ G₀'(x)² ũ_k → Φ_k'(x)`. -/
theorem tendstoInProb_cform2_UtildeCol_self (k : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.UtildeCol w N k)) (PhihetDeriv (fun i => m.thetaAligned i k) c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have h := tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (m.sum_thetaW_sq_nonneg w k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (MultiTableModel.PhiDerivC (fun i => m.thetaAligned i k) c w)
    (PhihetDeriv (fun i => m.thetaAligned i k) c w x)
    (fun z hz => m.tendstoInProb_qform2C_UtildeCol w c hc hw hR hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PhiDerivC hc hw hbx)
  exact h

include hc hw hR hreg hG hpd hedge hx in
/-- **Field `uu`, off-diagonal, at real `x`**: `ũ_kᵀ G₀'(x) ũ_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform_UtildeCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.UtildeCol w N l)) 0 :=
  tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (m.sum_thetaW_sq_nonneg w l)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cformC_UtildeCol_ne w c hc hw hR hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hR hreg hG hpd hedge hx in
/-- **Field `uu2`, off-diagonal, at real `x`**: `ũ_kᵀ G₀'(x)² ũ_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2_UtildeCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.UtildeCol w N l)) 0 :=
  tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (m.sum_thetaW_sq_nonneg w l)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cform2C_UtildeCol_ne w c hc hw hR hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `ee`, diagonal, at real `x`**: `g_kᵀ G₀'(x) g_k → Ψ(x)`. The rank-one mirror is
`MultiTableModel.tendstoInProb_qform_eHet`. -/
theorem tendstoInProb_cform_gHetCol_self (k : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω k)) (Psihet c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have h := tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (MultiTableModel.PsiC c w) (Psihet c w x)
    (fun z hz => m.tendstoInProb_qformC_gHetCol w c hc hw hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PsiC hc hw hbx)
  exact h

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `ee2`, diagonal, at real `x`**: `g_kᵀ G₀'(x)² g_k → Ψ'(x)`. -/
theorem tendstoInProb_cform2_gHetCol_self (k : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω k)) (PsihetDeriv c w x) := by
  have hbx : bHet c w < x := lt_of_le_of_lt hb hx
  have h := tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (MultiTableModel.PsiDerivC c w) (PsihetDeriv c w x)
    (fun z hz => m.tendstoInProb_qform2C_gHetCol w c hc hw hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PsiDerivC hc hw hbx)
  exact h

include hc hw hreg hG hpd hedge hx in
/-- **Field `ee`, off-diagonal, at real `x`**: `g_kᵀ G₀'(x) g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform_gHetCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cformC_gHetCol_ne w c hc hw hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hreg hG hpd hedge hx in
/-- **Field `ee2`, off-diagonal, at real `x`**: `g_kᵀ G₀'(x)² g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2_gHetCol_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cform2C_gHetCol_ne w c hc hw hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hR hreg hG hpd hedge hx in
/-- **Field `ue` at real `x`**: `ũ_kᵀ G₀'(x) g_l → 0`. The rank-one mirror is
`MultiTableModel.tendstoInProb_cform_u0Het_eHet`. -/
theorem tendstoInProb_cform_UtildeCol_gHetCol (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cformC_UtildeCol_gHetCol w c hc hw hR hreg hG hpd hz k l) hx
    (by simp)

include hc hw hR hreg hG hpd hedge hx in
/-- **Field `ue2` at real `x`**: `ũ_kᵀ G₀'(x)² g_l → 0`. -/
theorem tendstoInProb_cform2_UtildeCol_gHetCol (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_thetaW_sq_nonneg w k) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le w hR k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cform2C_UtildeCol_gHetCol w c hc hw hR hreg hG hpd hz k l)
    hx (by simp)

end RealForms

/-! ### 10. The Gaussian discharge -/

/-- **`ResolventLimitsHetR` on the exactly aligned Gaussian family.** `rk = alignedRk M r`,
`R_i = 1`, `d N = p N + r`, the edge as a raw hypothesis (stage E3 supplies it). The limits
are `Φ_k = MPhet.Phihet (θ_·k) c w`, `Ψ = MPhet.Psihet c w` and their derivatives. The
rank-one mirror is `MultiTableModel.resolventLimitsHet_of_gaussian` (`RMT/Het/R5het.lean`).
Paper: `thm:rank_r_stacksvd` (`main_paper.tex:2337`); the rank-one mirror is
`thm:stacksvd_weighted` (`:463`). -/
theorem resolventLimitsHetR_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {b : ℝ}
    (hb : MPhet.bHet c w ≤ b)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1)) :
    m.ResolventLimitsHetR w b
      (fun k => MPhet.Phihet (fun i => m.thetaAligned i k) c w) (MPhet.Psihet c w)
      (fun k => MPhet.PhihetDeriv (fun i => m.thetaAligned i k) c w)
      (MPhet.PsihetDeriv c w) := by
  refine ⟨hedge, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform_UtildeCol_self w c hc hw hR hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform_UtildeCol_ne w c hc hw hR hreg hG hpd hedge hx hkl
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform_gHetCol_self w c hc hw hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform_gHetCol_ne w c hc hw hreg hG hpd hedge hx hkl
  · intro k l x hx
    exact m.tendstoInProb_cform_UtildeCol_gHetCol w c hc hw hR hreg hG hpd hedge hx k l
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform2_UtildeCol_self w c hc hw hR hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform2_UtildeCol_ne w c hc hw hR hreg hG hpd hedge hx hkl
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform2_gHetCol_self w c hc hw hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform2_gHetCol_ne w c hc hw hreg hG hpd hedge hx hkl
  · intro k l x hx
    exact m.tendstoInProb_cform2_UtildeCol_gHetCol w c hc hw hR hreg hG hpd hedge hx k l

/-! ### 11. The column Gram limit -/

section Gram

variable [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
  (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)

include hR hreg hG hpd in
/-- **The cross Gram entry** `ũ_k ⬝ᵥ g_l → 0`, by Chebyshev on the pair law: the linear form
`ũ_kᵀ diag(τ) x / √d` has variance `≤ w²_max ‖ũ_k‖² / d`. -/
theorem tendstoInProb_UtildeCol_dotProduct_gHetCol (k l : Fin r) :
    TendstoInProb μ (fun N ω => m.UtildeCol w N k ⬝ᵥ m.gHetCol w N ω l) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  have hdR : Tendsto (fun N => ((d N : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have hid : ∀ N (y : Fin (∑ i, n i N) → ℝ),
      m.UtildeCol w N k ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) y
        = ∑ j, (m.UtildeCol w N k j * tauOf w (blkStack n N) j * (Real.sqrt (d N))⁻¹)
            * ((fun t : ℝ => t) (y j)) := by
    intro N y
    simp only [dotProduct, gvec]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hbase : TendstoInProb μ (fun N ω => |m.UtildeCol w N k ⬝ᵥ
      gvec (tauOf w (blkStack n N)) (d N) (m.pairColR hG hpd N l ω).1|) 0 := by
    refine tendstoInProb_of_integral_sq_le_pair (pN := fun N => ∑ i, n i N) (qN := p)
      (m.pairColR hG hpd · l) (fun N => m.hasLaw_pairColR hG hpd N l)
      (fun N y _ => |m.UtildeCol w N k ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) y|)
      (Filter.Eventually.of_forall fun N => ?_) (fun _ _ _ => abs_nonneg _)
      (Filter.Eventually.of_forall fun N _ => ?_)
      (K := fun N => Scalars.wSqMax w / d N * ∑ i, (m.thetaAligned i k * w i) ^ 2) ?_
      (Filter.Eventually.of_forall fun N _ => ?_)
    · refine Measurable.abs ?_
      exact Finset.measurable_sum _ fun j _ =>
        ((measurable_gvec_apply _ _ j).comp measurable_fst).const_mul _
    · simp only [sq_abs, hid]
      exact centered_id.integrable_sq_sum _
    · have := ((tendsto_const_nhds (x := Scalars.wSqMax w)).div_atTop hdR).mul_const
        (∑ i, (m.thetaAligned i k * w i) ^ 2)
      simpa using this
    · simp only [sq_abs, hid]
      rw [centered_id.integral_sq_sum, integral_sq_gauss, one_mul]
      have hsq : ((Real.sqrt (d N))⁻¹) ^ 2 = 1 / d N := by
        rw [inv_pow, Real.sq_sqrt (Nat.cast_nonneg _), one_div]
      have hterm : ∀ j, (m.UtildeCol w N k j * tauOf w (blkStack n N) j
          * (Real.sqrt (d N))⁻¹) ^ 2
          ≤ Scalars.wSqMax w / d N * m.UtildeCol w N k j ^ 2 := by
        intro j
        have h1 : (m.UtildeCol w N k j * tauOf w (blkStack n N) j * (Real.sqrt (d N))⁻¹) ^ 2
            = tauOf w (blkStack n N) j ^ 2 / d N * m.UtildeCol w N k j ^ 2 := by
          rw [mul_pow, mul_pow, hsq]
          ring
        rw [h1]
        exact mul_le_mul_of_nonneg_right
          (div_le_div_of_nonneg_right (tauOf_sq_le w _ j) (Nat.cast_nonneg _)) (sq_nonneg _)
      calc ∑ j, (m.UtildeCol w N k j * tauOf w (blkStack n N) j * (Real.sqrt (d N))⁻¹) ^ 2
          ≤ ∑ j, Scalars.wSqMax w / d N * m.UtildeCol w N k j ^ 2 :=
            Finset.sum_le_sum fun j _ => hterm j
        _ = Scalars.wSqMax w / d N * ∑ i, (m.thetaAligned i k * w i) ^ 2 := by
            rw [← Finset.mul_sum, ← m.dotProduct_UtildeCol_self hR w N k]
            congr 1
            simp only [dotProduct, sq]
  refine TendstoInProb.of_le (g := fun N ω => |m.UtildeCol w N k ⬝ᵥ
      gvec (tauOf w (blkStack n N)) (d N) (m.pairColR hG hpd N l ω).1|)
    (fun N => Filter.Eventually.of_forall fun ω => ?_) hbase
  rw [sub_zero, m.gHetCol_eq_gvec w]
  exact le_rfl

include hw hreg hG hpd in
/-- **The off-diagonal Gaussian Gram entry** `g_k ⬝ᵥ g_l → 0`, `k ≠ l`, by real polarization
through the unit mix. -/
theorem tendstoInProb_gHetCol_dotProduct_ne {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l) 0 := by
  have hu := m.tendstoInProb_gHetMix_norm w c hw hreg hG hpd hkl
  have hk := m.tendstoInProb_gHetCol_norm w c hw hreg hG hpd k
  have hl := m.tendstoInProb_gHetCol_norm w c hw hreg hG hpd l
  have hfun : (fun N (ω : Ω N) => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l)
      = fun N ω => (2 * (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1
          ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
        - m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k
        - m.gHetCol w N ω l ⬝ᵥ m.gHetCol w N ω l) * 2⁻¹ := by
    funext N ω
    have h := m.gHetCol_add_eq w hG hpd N ω k l
    have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
    have hexp : (m.gHetCol w N ω k + m.gHetCol w N ω l)
        ⬝ᵥ (m.gHetCol w N ω k + m.gHetCol w N ω l)
        = 2 * (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1
          ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1) := by
      rw [h, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc, h2]
    rw [add_dotProduct, dotProduct_add, dotProduct_add,
      dotProduct_comm (m.gHetCol w N ω l) (m.gHetCol w N ω k)] at hexp
    linear_combination (1 / 2 : ℝ) * hexp
  rw [hfun]
  have hlim := (((hu.const_mul 2).sub hk).sub hl).mul_const 2⁻¹
  refine FormsR.tendstoInProb_congr_limit ?_ hlim
  ring

include hw hR hreg hG hpd in
/-- **The column Gram limit of `Q`**: `Q_k ⬝ᵥ Q_l → δ_kl (∑ w_i² θ_ik² + ∑ w_i² c_i)`.
Paper: the normalization `Qᵀ Q → diag` of the weighted rank-`r` estimator of
`thm:rank_r_stacksvd` (`main_paper.tex:2337`); the rank-one mirror is
`thm:stacksvd_weighted` (`:463`) and `MultiTableModel.tendstoInProb_qHet_norm`
(`RMT/Het/R5het.lean`). -/
theorem tendstoInProb_dotProduct_QmatHetR_col (k l : Fin r) :
    TendstoInProb μ (fun N ω =>
        (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      (if k = l then ∑ i, w i ^ 2 * m.thetaAligned i k ^ 2 + ∑ i, w i ^ 2 * c i else 0) := by
  have hfun : (fun N (ω : Ω N) =>
        (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      = fun N ω => (if k = l then ∑ i, (m.thetaAligned i k * w i) ^ 2 else 0)
          + m.UtildeCol w N k ⬝ᵥ m.gHetCol w N ω l
          + m.UtildeCol w N l ⬝ᵥ m.gHetCol w N ω k
          + m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l := by
    funext N ω
    rw [m.QmatHetR_col w N ω k, m.QmatHetR_col w N ω l, add_dotProduct, dotProduct_add,
      dotProduct_add, m.dotProduct_UtildeCol hR w N k l,
      dotProduct_comm (m.gHetCol w N ω k) (m.UtildeCol w N l)]
    ring
  rw [hfun]
  have h1 := m.tendstoInProb_UtildeCol_dotProduct_gHetCol w c hR hreg hG hpd k l
  have h2 := m.tendstoInProb_UtildeCol_dotProduct_gHetCol w c hR hreg hG hpd l k
  rcases eq_or_ne k l with rfl | hkl
  · have h3 := m.tendstoInProb_gHetCol_norm w c hw hreg hG hpd k
    have hcomb := (((TendstoInProb.const μ (∑ i, (m.thetaAligned i k * w i) ^ 2)).add h1).add
      h2).add h3
    simp only [if_true]
    refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    rw [add_zero, add_zero]
    congr 1
    exact Finset.sum_congr rfl fun i _ => by ring
  · have h3 := m.tendstoInProb_gHetCol_dotProduct_ne w c hw hreg hG hpd hkl
    have hcomb := (((TendstoInProb.const μ (0 : ℝ)).add h1).add h2).add h3
    simp only [if_neg hkl]
    refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    ring

end Gram

end UnalignedModelR
end StackedSVD
