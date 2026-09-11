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
import StackedSVD.RankR.Het.Forms
import StackedSVD.RankR.SingleWeight.Scalars

/-!
# The matrix-valued resolvent forms of the single-weight stack (Track G, unit G1)

Mirrors `RankR/Het/Forms.lean` (the exactly aligned family, `rk = alignedRk M r`, `R_i = 1`)
at the general model `UnalignedModelR μ M n d r rk` of arbitrary alignment matrices `R_i`.
Paper: `prop:gen_rank_stacksvd_singleweight` (`main_paper.tex:2112`); the deterministic layer
is `RankR/SingleWeight/Scalars.lean`. Campaign plan: `notes/archive/trackG_plan.md`, section "G1
`Het/Forms.lean`". Date 2026-09-05.

At general `R_i` the deterministic signal column `Ũ_k` of table `i` reads
`w_i θ_ij (R_i)_{kj}` summed over `j`, so `Ũ_kᵀ Ũ_l` restricted to block `i` is
`w_i² (R_i Θ_i² R_iᵀ)_{kl} = w_i² S_i(k,l)` with `S_i := SingleWeight.sigMat`, not the
`δ_kl (θ_ik w_i)²` of the aligned case. The resolvent limit `Ũ_kᵀ G₀' Ũ_l` therefore need not
vanish off the diagonal, and the hypothesis structure `ResolventLimitsHetR` is generalized to
`ResolventLimitsSW` by making its `uu`, `uu2` fields matrix valued (`PhiM k l z`, `PhiM2 k l z`
for every pair `(k, l)`, not `if k = l then Phi k z else 0`). The Gaussian columns `g_k` never
read `R_i`, so the `ee`, `ue`, `ee2`, `ue2` fields keep their aligned shape and their proofs
are copied unchanged, retyped at general `rk`.

Objects.

1. `sum_blockSet_UtildeCol_mul_gen` and its two corollaries: the block bilinear sum of `Ũ` at
   general `R_i`, replacing the aligned `δ_kl (θ_ik w_i)²` with `w_i² S_i(k, l)`.
2. `tendstoInProb_cformC_cross3'`, `tendstoInProb_cform2C_cross3'`: the three-limit
   polarization step of `RankR/Het/Forms.lean:424,466` with the relation `2 Lu = La + Lb`
   dropped and the conclusion generalized to the value `(2 Lu - La - Lb)/2` the same algebra
   produces (needed because `Ũ_k`, `Ũ_l` no longer share the aligned coincidence).
3. `ResolventLimitsSW`, `ResolventLimitsSW.cform_qcol`, `.cform2_qcol`: the hypothesis
   structure and its column corollaries.
4. `tendstoInProb_cformC_UtildeCol_gen` and its real-axis form `tendstoInProb_cform_UtildeCol_gen`:
   the matrix-valued `uu` limit at every `(k, l)` in one lemma (the `k = l` case reduces to
   `HetR2.tendstoInProb_qformC_fixed`; `k ≠ l` goes through the unit mix and `cross3'`).
5. The `_gen` copies of the aligned `ee`/`ue` machinery at general `rk`, and
   `resolventLimitsSW_of_gaussian`: the Gaussian discharge, with the edge as a raw hypothesis
   exactly as `resolventLimitsHetR_of_gaussian` (stage E3's analog is not imported here).

Modeling choices: none beyond `RankR/Het/Forms.lean`'s own (items 1-3 of that file's header);
this file adds no new model assumption, only the general-`R_i` algebra.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace UnalignedModelR

open HetR2 HetStein HetR1 MPhet R2 R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. Block sums of `Ũ` at general `R_i` -/

/-- `Ũ` read at a row of block `i`, general `R_i`: `w_i ∑_j (U_i)_{aj} θ_ij (R_i)_{kj}`.
General analog of `UtildeCol_apply_aligned` (`RankR/Het/Forms.lean:100`). -/
theorem UtildeCol_apply_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k : Fin r) (i : Fin M) (a : Fin (n i N)) :
    m.UtildeCol w N k (finSigmaFinEquiv ⟨i, a⟩)
      = w i * ∑ j : Fin (rk i), (m.tbl i).U N a j * (m.tbl i).θ j * m.R i k j := by
  have hσ : finSigmaFinEquiv.symm (finSigmaFinEquiv (⟨i, a⟩ : Σ i, Fin (n i N)))
      = (⟨i, a⟩ : Σ i, Fin (n i N)) := Equiv.symm_apply_apply _ _
  simp only [UtildeCol, UtildeR, SigmaHalfR, Matrix.diagonal_mul, signalFactorG_apply']
  rw [hσ]

/-- Entrywise unfolding of `SingleWeight.sigMat`: `S_i(k, l) = ∑_j (R_i)_{kj} θ_ij² (R_i)_{lj}`.
Private algebraic helper for `sum_blockSet_UtildeCol_mul_gen`. -/
theorem sigMat_apply_gen {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) (k l : Fin r) :
    SingleWeight.sigMat θ R i k l = ∑ j : Fin (rk i), R i k j * θ i j ^ 2 * R i l j := by
  simp only [SingleWeight.sigMat, Matrix.mul_apply, Matrix.transpose_apply]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [Matrix.diagonal_apply, mul_ite, mul_zero, Finset.sum_ite_eq']
  simp

/-- **Block bilinear sums of `Ũ`** at general `R_i`:
`∑_{q ∈ J_i} ũ_k(q) ũ_l(q) = w_i² S_i(k, l)`, from the orthonormal columns of `U_i` and the
matrix identity `R_i Θ_i² R_iᵀ = S_i`. At `R_i = 1` this is `δ_kl (θ_ik w_i)²`, the aligned
`sum_blockSet_UtildeCol_mul` (`RankR/Het/Forms.lean:115`). -/
theorem sum_blockSet_UtildeCol_mul_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (i : Fin M) (k l : Fin r) :
    ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q * m.UtildeCol w N l q
      = w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l := by
  classical
  have hU : ∀ j j' : Fin (rk i), ∑ a : Fin (n i N), (m.tbl i).U N a j * (m.tbl i).U N a j'
      = if j = j' then 1 else 0 := by
    intro j j'
    have h := congrFun (congrFun ((m.tbl i).hU N) j) j'
    simpa [Matrix.mul_apply, Matrix.one_apply] using h
  rw [MultiTableModel.sum_blockSet_blkStack, sigMat_apply_gen]
  simp only [m.UtildeCol_apply_gen w N _ i]
  have hpt : ∀ a : Fin (n i N),
      (w i * ∑ j : Fin (rk i), (m.tbl i).U N a j * (m.tbl i).θ j * m.R i k j)
        * (w i * ∑ j' : Fin (rk i), (m.tbl i).U N a j' * (m.tbl i).θ j' * m.R i l j')
      = w i ^ 2 * ∑ j : Fin (rk i), ∑ j' : Fin (rk i),
          ((m.tbl i).θ j * m.R i k j) * ((m.tbl i).θ j' * m.R i l j')
            * ((m.tbl i).U N a j * (m.tbl i).U N a j') := by
    intro a
    calc (w i * ∑ j : Fin (rk i), (m.tbl i).U N a j * (m.tbl i).θ j * m.R i k j)
          * (w i * ∑ j' : Fin (rk i), (m.tbl i).U N a j' * (m.tbl i).θ j' * m.R i l j')
        = w i ^ 2 * ((∑ j : Fin (rk i), (m.tbl i).U N a j * (m.tbl i).θ j * m.R i k j)
            * (∑ j' : Fin (rk i), (m.tbl i).U N a j' * (m.tbl i).θ j' * m.R i l j')) := by ring
      _ = w i ^ 2 * ∑ j : Fin (rk i), ∑ j' : Fin (rk i),
            ((m.tbl i).θ j * m.R i k j) * ((m.tbl i).θ j' * m.R i l j')
              * ((m.tbl i).U N a j * (m.tbl i).U N a j') := by
          congr 1
          rw [Finset.sum_mul]
          refine Finset.sum_congr rfl fun j _ => ?_
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun j' _ => by ring
  rw [Finset.sum_congr rfl fun a (_ : a ∈ Finset.univ) => hpt a]
  rw [← Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  calc ∑ j : Fin (rk i), ∑ a : Fin (n i N), ∑ j' : Fin (rk i),
        ((m.tbl i).θ j * m.R i k j) * ((m.tbl i).θ j' * m.R i l j')
          * ((m.tbl i).U N a j * (m.tbl i).U N a j')
      = ∑ j : Fin (rk i), ∑ j' : Fin (rk i), ∑ a : Fin (n i N),
          ((m.tbl i).θ j * m.R i k j) * ((m.tbl i).θ j' * m.R i l j')
            * ((m.tbl i).U N a j * (m.tbl i).U N a j') := by
        exact Finset.sum_congr rfl fun j _ => Finset.sum_comm
    _ = ∑ j : Fin (rk i), ∑ j' : Fin (rk i),
          ((m.tbl i).θ j * m.R i k j) * ((m.tbl i).θ j' * m.R i l j')
            * (if j = j' then (1 : ℝ) else 0) := by
        refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun j' _ => ?_
        rw [← Finset.mul_sum, hU]
    _ = ∑ j : Fin (rk i), m.R i k j * (m.tbl i).θ j ^ 2 * m.R i l j := by
        refine Finset.sum_congr rfl fun j _ => ?_
        simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, if_true]
        ring

/-- The block norms of `ũ_k` at general `R_i`: `∑_{q ∈ J_i} ũ_k(q)² = w_i² S_i(k, k)`. -/
theorem sum_blockSet_UtildeCol_sq_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (i : Fin M) (k : Fin r) :
    ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q ^ 2
      = w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k := by
  have h := m.sum_blockSet_UtildeCol_mul_gen w N i k k
  rw [← h]
  exact Finset.sum_congr rfl fun q _ => sq _

/-- The block norms of the unit mix `(ũ_k + ũ_l)/√2` at general `R_i`. Unlike the aligned
`sum_blockSet_UtildeCol_mix_sq` this holds whether or not `k = l` (the hypothesis `k ≠ l` of the
mirrored plan was unused and is dropped, 2026-09-07). -/
theorem sum_blockSet_UtildeCol_mix_sq_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (i : Fin M) (k l : Fin r) :
    ∑ q ∈ blockSet (blkStack n N) i,
        ((Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l)) q ^ 2
      = w i ^ 2 * (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k
          + SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l
          + 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) / 2 := by
  have hkk := m.sum_blockSet_UtildeCol_mul_gen w N i k k
  have hll := m.sum_blockSet_UtildeCol_mul_gen w N i l l
  have hkl' := m.sum_blockSet_UtildeCol_mul_gen w N i k l
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

/-- `ũ_k ⬝ᵥ ũ_l = ∑_i w_i² S_i(k, l)` at general `R_i`. General analog of
`dotProduct_UtildeCol` (`RankR/Het/Forms.lean:169`). -/
theorem dotProduct_UtildeCol_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k l : Fin r) :
    m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N l
      = ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l := by
  have h : m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N l
      = ∑ i, ∑ q ∈ blockSet (blkStack n N) i, m.UtildeCol w N k q * m.UtildeCol w N l q := by
    rw [sum_blockSet_real]
    rfl
  rw [h]
  exact Finset.sum_congr rfl fun i _ => m.sum_blockSet_UtildeCol_mul_gen w N i k l

/-- The squared norm of `ũ_k` at general `R_i`. -/
theorem dotProduct_UtildeCol_self_gen (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k : Fin r) :
    m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k
      = ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k :=
  m.dotProduct_UtildeCol_gen w N k k

/-- The diagonal block-sum total is nonnegative: it is a squared norm. Feeds the `Kp`, `Kq`
bounds of the real-axis transfer and the `Ky` bound of the `ue` cross form. -/
theorem sum_wSqSigMat_diag_nonneg (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (k : Fin r) :
    0 ≤ ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k := by
  rw [← m.dotProduct_UtildeCol_self_gen w N k]
  exact Finset.sum_nonneg fun _ _ => mul_self_nonneg _

/-! ### 2. General boundary limits, reused for both `uu` and `ee` -/

/-- **General boundary limit.** For any real family `a`, the vertical-line limit of
`∑ i, a_i gC_i(z, sGlob(z))` as `z → x⁺` is `∑ i, a_i ghet_i(x)`. Specializes to
`MultiTableModel.tendsto_PhiC`/`tendsto_PsiC` (`RMT/Het/R2het.lean:2147,2160`) at
`a_i = (θ_i w_i)^2` resp. `c_i w_i^2`, and drives `uu`'s transfer at `a_i = w_i^2 S_i(k, l)`. -/
theorem tendsto_sum_gC_boundary {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {x : ℝ} (hx : MPhet.bHet c w < x) (a : Fin M → ℝ) :
    Tendsto (fun η : ℝ => ∑ i, ((a i : ℝ) : ℂ) * gC w i ((x : ℂ) + (η : ℂ) * Complex.I)
        (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I))) (𝓝[>] 0)
      (𝓝 ((∑ i, a i * MPhet.ghet c w i x : ℝ) : ℂ)) := by
  have hlim := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gC_sGlob hc hw hx i).const_mul (((a i : ℝ) : ℂ))
  have heq : ∑ i, ((a i : ℝ) : ℂ) * ((MPhet.ghet c w i x : ℝ) : ℂ)
      = ((∑ i, a i * MPhet.ghet c w i x : ℝ) : ℂ) := by
    simp only [Complex.ofReal_sum, Complex.ofReal_mul]
  rw [← heq]
  exact hlim

/-- The same, second order. -/
theorem tendsto_sum_gCDeriv_boundary {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) {x : ℝ} (hx : MPhet.bHet c w < x) (a : Fin M → ℝ) :
    Tendsto (fun η : ℝ => ∑ i, ((a i : ℝ) : ℂ)
        * gCDeriv w i ((x : ℂ) + (η : ℂ) * Complex.I)
          (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I))
          (zfunDerivC c w (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I)))⁻¹) (𝓝[>] 0)
      (𝓝 ((∑ i, a i * MPhet.ghetDeriv c w i x : ℝ) : ℂ)) := by
  have hlim := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
    (tendsto_gCDeriv_sGlob hc hw hx i).const_mul (((a i : ℝ) : ℂ))
  have heq : ∑ i, ((a i : ℝ) : ℂ) * ((MPhet.ghetDeriv c w i x : ℝ) : ℂ)
      = ((∑ i, a i * MPhet.ghetDeriv c w i x : ℝ) : ℂ) := by
    simp only [Complex.ofReal_sum, Complex.ofReal_mul]
  rw [← heq]
  exact hlim

/-! ### 3. Polarization with three unrelated limits -/

section Polarization

variable {dN : ℕ → ℕ} {z : ℂ}

/-- **The polarization step with three limits, no relation assumed between them.**
`a + b = √2 u`, `qformC a → La`, `qformC b → Lb`, `qformC u → Lu` give
`cformC a b → (2 Lu - La - Lb)/2`. Mirrors `tendstoInProb_cformC_cross3`
(`RankR/Het/Forms.lean:424`) with the hypothesis `2 Lu = La + Lb` dropped: the `heq` step of
that proof already produces the stated value from `hpol`, `hkey` alone, by `ring`. -/
theorem tendstoInProb_cformC_cross3'
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {La Lb Lu : ℂ}
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (af N ω) - La‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (bf N ω) - Lb‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (uf N ω) - Lu‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (Wf N ω) z (af N ω) (bf N ω)
      - (2 * Lu - La - Lb) / 2‖) 0 := by
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
    have heq : R4C.cformC (Wf N ω) z (af N ω) (bf N ω) - (2 * Lu - La - Lb) / 2
        = ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - Lu)
          - (R4C.qformC (Wf N ω) z (af N ω) - La)
          - (R4C.qformC (Wf N ω) z (bf N ω) - Lb)) / 2 := by
      rw [hpol, hkey]
      ring
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
theorem tendstoInProb_cform2C_cross3'
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {La Lb Lu : ℂ}
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (af N ω) - La‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (bf N ω) - Lb‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (uf N ω) - Lu‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (Wf N ω) z (af N ω) (bf N ω)
      - (2 * Lu - La - Lb) / 2‖) 0 := by
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
    have heq : R4C.cform2C (Wf N ω) z (af N ω) (bf N ω) - (2 * Lu - La - Lb) / 2
        = ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - Lu)
          - (R4C.qform2C (Wf N ω) z (af N ω) - La)
          - (R4C.qform2C (Wf N ω) z (bf N ω) - Lb)) / 2 := by
      rw [hpol, hkey]
      ring
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

/-! ### 4. The hypothesis structure and its column corollaries -/

/-- **The rank-`r` heteroscedastic (H1) and (H2) at general `R_i`**, in the shape of
`ResolventLimitsHetR` (`RankR/Het/Forms.lean:200`) with the two `uu` fields made matrix
valued: `PhiM k l z` is the limit of `Ũ_kᵀ G₀' Ũ_l` for every pair `(k, l)`, not only `k = l`.
The `ee`, `ue` fields are unchanged, since the Gaussian columns `g_k` never read `R_i`. -/
structure ResolventLimitsSW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (b : ℝ)
    (PhiM PhiM2 : Fin r → Fin r → ℝ → ℝ) (Psi Psi2 : ℝ → ℝ) : Prop where
  /-- (H1), the edge bound. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1)
  /-- (H2), `ũ_kᵀ G₀' ũ_l → PhiM k l`. -/
  uu : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l)) (PhiM k l z)
  /-- (H2), `g_kᵀ G₀' g_l → δ_kl Ψ`. -/
  ee : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l))
    (if k = l then Psi z else 0)
  /-- (H2), the cross form `ũ_kᵀ G₀' g_l → 0`. -/
  ue : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)) 0
  /-- (H2), `ũ_kᵀ G₀'² ũ_l → PhiM2 k l`. -/
  uu2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.UtildeCol w N l)) (PhiM2 k l z)
  /-- (H2), `g_kᵀ G₀'² g_l → δ_kl Ψ'`. -/
  ee2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.gHetCol w N ω k) (m.gHetCol w N ω l))
    (if k = l then Psi2 z else 0)
  /-- (H2), the squared cross form `ũ_kᵀ G₀'² g_l → 0`. -/
  ue2 : ∀ k l : Fin r, ∀ z, b < z → TendstoInProb μ
    (fun N ω => cform2 (m.W0hetR w N ω) z (m.UtildeCol w N k) (m.gHetCol w N ω l)) 0

/-- **Entry `(k, l)` of `Qᵀ G₀'(z) Q`** tends to `PhiM k l z + δ_kl Ψ(z)`. Mirror of
`ResolventLimitsHetR.cform_qcol` (`RankR/Het/Forms.lean:232`); no case split on `k = l` is
needed here since `h.uu` already gives the value for every pair. -/
theorem ResolventLimitsSW.cform_qcol {m : UnalignedModelR μ M n d r rk} {w : Fin M → ℝ}
    {b : ℝ} {PhiM PhiM2 : Fin r → Fin r → ℝ → ℝ} {Psi Psi2 : ℝ → ℝ}
    (h : m.ResolventLimitsSW w b PhiM PhiM2 Psi Psi2) (k l : Fin r) {z : ℝ} (hz : b < z) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      (PhiM k l z + if k = l then Psi z else 0) := by
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
  refine FormsR.tendstoInProb_congr_limit ?_ hcomb
  ring

/-- **Entry `(k, l)` of `Qᵀ G₀'(z)² Q`** tends to `PhiM2 k l z + δ_kl Ψ'(z)`. -/
theorem ResolventLimitsSW.cform2_qcol {m : UnalignedModelR μ M n d r rk} {w : Fin M → ℝ}
    {b : ℝ} {PhiM PhiM2 : Fin r → Fin r → ℝ → ℝ} {Psi Psi2 : ℝ → ℝ}
    (h : m.ResolventLimitsSW w b PhiM PhiM2 Psi Psi2) (k l : Fin r) {z : ℝ} (hz : b < z) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) z
        (fun q => m.QmatHetR w N ω q k) (fun q => m.QmatHetR w N ω q l))
      (PhiM2 k l z + if k = l then Psi2 z else 0) := by
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
  refine FormsR.tendstoInProb_congr_limit ?_ hcomb
  ring

/-! ### 5. `uu`, `uu2`: the deterministic columns at general `R_i` -/

section ComplexUU

variable [NeZero M] (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
  (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {z : ℂ} (hz : 0 < z.im)

include hc hw hreg hG hpd hz in
/-- **`uu` at complex `z`, both diagonal and off-diagonal in one lemma.** The `k = l` case is
`HetR2.tendstoInProb_qformC_fixed` (`cformC v v = qformC v` definitionally); the `k ≠ l` case
goes through the unit mix `(ũ_k + ũ_l)/√2` and `tendstoInProb_cformC_cross3'`. -/
theorem tendstoInProb_cformC_UtildeCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.UtildeCol w N k)
        (m.UtildeCol w N l)
      - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
          * gC w i z (sGlob c w z)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  rcases eq_or_ne k l with rfl | hkl
  · have h := tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N k)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i k)
    simpa only [m.W0hetR_eq_Wsig_blockBR hG hpd w, ← FormsR.qformC_eq_cformC] using h
  · have hk := tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N k)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i k)
    have hl := tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N l)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i l)
    have hmix := tendstoInProb_qformC_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      (fun i => w i ^ 2 * (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k
        + SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l
        + 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) / 2)
      (fun N i => m.sum_blockSet_UtildeCol_mix_sq_gen w N i k l)
    have hcross := tendstoInProb_cformC_cross3' (dN := fun N => ∑ i, n i N)
      (Wf := fun N ω => Wsig (tauOf w (blkStack n N)) (d N) (m.blockBR hG hpd N ω))
      (fun N ω => isHermitian_Wsig _ _ _) hz
      (af := fun N _ => m.UtildeCol w N k) (bf := fun N _ => m.UtildeCol w N l)
      (uf := fun N _ => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      (fun _ _ => add_eq_sqrt_two_smul_mix _ _) hk hl hmix
    have heq : ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
        * gC w i z (sGlob c w z)
        = (2 * (∑ i, ((w i ^ 2 * (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k
              + SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l
              + 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) / 2 : ℝ) : ℂ)
              * gC w i z (sGlob c w z))
            - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k : ℝ) : ℂ)
              * gC w i z (sGlob c w z)
            - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l : ℝ) : ℂ)
              * gC w i z (sGlob c w z)) / 2 := by
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib, Finset.sum_div]
      refine Finset.sum_congr rfl fun i _ => ?_
      push_cast
      ring
    simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, heq]
    exact hcross

include hc hw hreg hG hpd hz in
/-- **`uu2` at complex `z`**, the second-order twin of `tendstoInProb_cformC_UtildeCol_gen`. -/
theorem tendstoInProb_cform2C_UtildeCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.UtildeCol w N k)
        (m.UtildeCol w N l)
      - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
          * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  rcases eq_or_ne k l with rfl | hkl
  · have h := tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N k)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i k)
    simpa only [m.W0hetR_eq_Wsig_blockBR hG hpd w, ← FormsR.qform2C_eq_cform2C] using h
  · have hk := tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N k)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i k)
    have hl := tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => m.UtildeCol w N l)
      (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l)
      (fun N i => m.sum_blockSet_UtildeCol_sq_gen w N i l)
    have hmix := tendstoInProb_qform2C_fixed (pN := fun N => ∑ i, n i N) (qN := p)
      (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
      (m.tendsto_cN_blkStackR hreg) (m.blockBR hG hpd) (m.hasLaw_blockBR hG hpd)
      (fun N => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      (fun i => w i ^ 2 * (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k
        + SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l
        + 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) / 2)
      (fun N i => m.sum_blockSet_UtildeCol_mix_sq_gen w N i k l)
    have hcross := tendstoInProb_cform2C_cross3' (dN := fun N => ∑ i, n i N)
      (Wf := fun N ω => Wsig (tauOf w (blkStack n N)) (d N) (m.blockBR hG hpd N ω))
      (fun N ω => isHermitian_Wsig _ _ _) hz
      (af := fun N _ => m.UtildeCol w N k) (bf := fun N _ => m.UtildeCol w N l)
      (uf := fun N _ => (Real.sqrt 2)⁻¹ • (m.UtildeCol w N k + m.UtildeCol w N l))
      (fun _ _ => add_eq_sqrt_two_smul_mix _ _) hk hl hmix
    have heq : ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
        * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹
        = (2 * (∑ i, ((w i ^ 2 * (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k
              + SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l
              + 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) / 2 : ℝ) : ℂ)
              * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹)
            - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k : ℝ) : ℂ)
              * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹
            - ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l : ℝ) : ℂ)
              * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹) / 2 := by
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib, Finset.sum_div]
      refine Finset.sum_congr rfl fun i _ => ?_
      push_cast
      ring
    simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, heq]
    exact hcross

end ComplexUU

/-- **The norm event of `ũ_k` at general `R_i`**: the bound holds for every `N`, `ω`, so the
event is the whole space. General analog of `tendsto_measure_UtildeCol_norm_le`
(`RankR/Het/Forms.lean:821`) without `hR`. -/
theorem tendsto_measure_UtildeCol_norm_le_gen [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (k : Fin r) :
    Tendsto (fun N => μ N {_ω : Ω N | m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k
      ≤ ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k}) atTop (𝓝 1) := by
  have heq : ∀ N, {_ω : Ω N | m.UtildeCol w N k ⬝ᵥ m.UtildeCol w N k
      ≤ ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k} = Set.univ :=
    fun N => Set.eq_univ_of_forall fun _ => (m.dotProduct_UtildeCol_self_gen w N k).le
  simp only [heq, measure_univ]
  exact tendsto_const_nhds

section RealUU

variable [∀ N, IsProbabilityMeasure (μ N)] [NeZero M] (m : UnalignedModelR μ M n d r rk)
  (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {b : ℝ} (hb : MPhet.bHet c w ≤ b)
  (hedge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1))
  {x : ℝ} (hx : b < x)

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `uu` of `ResolventLimitsSW`, at real `x`, both diagonal and off-diagonal.** Real
transfer of `tendstoInProb_cformC_UtildeCol_gen` through
`HetR2.tendstoInProb_cform_of_complex_scaled`, with the norm bounds of `Ũ_k`, `Ũ_l` from
`sum_wSqSigMat_diag_nonneg`/`tendsto_measure_UtildeCol_norm_le_gen`. -/
theorem tendstoInProb_cform_UtildeCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.UtildeCol w N k) (m.UtildeCol w N l))
      (∑ i, w i ^ 2 * MPhet.ghet c w i x
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) := by
  have hbx : MPhet.bHet c w < x := lt_of_le_of_lt hb hx
  exact tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_wSqSigMat_diag_nonneg w 0 k) (m.sum_wSqSigMat_diag_nonneg w 0 l)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w l)
    (fun z => ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
        * gC w i z (sGlob c w z))
    (∑ i, w i ^ 2 * MPhet.ghet c w i x * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
    (fun z hz => m.tendstoInProb_cformC_UtildeCol_gen w c hc hw hreg hG hpd hz k l) hx
    (by
      have heq : (∑ i, w i ^ 2 * MPhet.ghet c w i x
            * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
          = ∑ i, (w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
              * MPhet.ghet c w i x :=
        Finset.sum_congr rfl fun i _ => by ring
      rw [heq]
      exact tendsto_sum_gC_boundary hc hw hbx
        (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l))

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `uu2` of `ResolventLimitsSW`, at real `x`.** -/
theorem tendstoInProb_cform2_UtildeCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.UtildeCol w N k)
        (m.UtildeCol w N l))
      (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i x
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l) := by
  have hbx : MPhet.bHet c w < x := lt_of_le_of_lt hb hx
  exact tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N _ => m.UtildeCol w N l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_wSqSigMat_diag_nonneg w 0 k) (m.sum_wSqSigMat_diag_nonneg w 0 l)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w k)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w l)
    (fun z => ∑ i, ((w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l : ℝ) : ℂ)
        * gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹)
    (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i x
      * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
    (fun z hz => m.tendstoInProb_cform2C_UtildeCol_gen w c hc hw hreg hG hpd hz k l) hx
    (by
      have heq : (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i x
            * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
          = ∑ i, (w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
              * MPhet.ghetDeriv c w i x :=
        Finset.sum_congr rfl fun i _ => by ring
      rw [heq]
      exact tendsto_sum_gCDeriv_boundary hc hw hbx
        (fun i => w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l))

end RealUU

/-! ### 6. `ee`, `ue`: the Gaussian columns, general `rk` copies

None of these proofs read `θ`, `R_i`, or the block-sum lemmas of section 1: the Gaussian
columns `g_k` depend only on `w`, `c`. The aligned file states them at `rk = alignedRk M r`
only because that is its section variable; the bodies below are copied unchanged. -/

section ComplexEE

variable [NeZero M] (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
  (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {z : ℂ} (hz : 0 < z.im)

include hc hw hreg hG hpd hz in
/-- `g_kᵀ G₀'(z) g_k → Ψ(z)` at complex `z`. General copy of `tendstoInProb_qformC_gHetCol`
(`RankR/Het/Forms.lean:654`). -/
theorem tendstoInProb_qformC_gHetCol_gen (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      - MultiTableModel.PsiC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_qformC_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

include hc hw hreg hG hpd hz in
/-- `g_kᵀ G₀'(z)² g_k → Ψ'(z)` at complex `z`. -/
theorem tendstoInProb_qform2C_gHetCol_gen (k : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      - MultiTableModel.PsiDerivC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_qform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

include hc hw hreg hG hpd hz in
/-- The unit mix `(g_k + g_l)/√2`, `k ≠ l`: `qformC → Ψ` at complex `z`. -/
theorem tendstoInProb_qformC_gHetMix_gen {k l : Fin r} (hkl : k ≠ l) :
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
theorem tendstoInProb_qform2C_gHetMix_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (m.W0hetR w N ω) z
      (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
      - MultiTableModel.PsiDerivC c w z‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w]
  exact tendstoInProb_qform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairMixR hG hpd · k l 1)
    (fun N => m.hasLaw_pairMixR hG hpd N hkl (Or.inl rfl))

/-- `g_k + g_l = √2 gvec τ d (x_k + x_l)/√2`, the sum in the shape of `pairMixR`. General copy
of `gHetCol_add_eq` (`RankR/Het/Forms.lean:701`). -/
theorem gHetCol_add_eq_gen (N : ℕ) (ω : Ω N) (k l : Fin r) :
    m.gHetCol w N ω k + m.gHetCol w N ω l
      = Real.sqrt 2 • gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1 := by
  rw [m.gHetCol_eq_gvec w, m.gHetCol_eq_gvec w]
  exact gvec_add_eq_sqrt_two_smul _ _ _ _

include hc hw hreg hG hpd hz in
/-- **Off-diagonal `ee` at complex `z`**: `g_kᵀ G₀'(z) g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cformC_gHetCol_ne_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)‖) 0 :=
  FormsR.tendstoInProb_cformC_cross (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N ω => m.gHetCol w N ω k) (bf := fun N ω => m.gHetCol w N ω l)
    (uf := fun N ω => gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
    (fun N ω => m.gHetCol_add_eq_gen w hG hpd N ω k l)
    (m.tendstoInProb_qformC_gHetCol_gen w c hc hw hreg hG hpd hz k)
    (m.tendstoInProb_qformC_gHetCol_gen w c hc hw hreg hG hpd hz l)
    (m.tendstoInProb_qformC_gHetMix_gen w c hc hw hreg hG hpd hz hkl)

include hc hw hreg hG hpd hz in
/-- **Off-diagonal `ee2` at complex `z`**: `g_kᵀ G₀'(z)² g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2C_gHetCol_ne_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)‖) 0 :=
  FormsR.tendstoInProb_cform2C_cross (dN := fun N => ∑ i, n i N)
    (Wf := fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω) hz
    (af := fun N ω => m.gHetCol w N ω k) (bf := fun N ω => m.gHetCol w N ω l)
    (uf := fun N ω => gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
    (fun N ω => m.gHetCol_add_eq_gen w hG hpd N ω k l)
    (m.tendstoInProb_qform2C_gHetCol_gen w c hc hw hreg hG hpd hz k)
    (m.tendstoInProb_qform2C_gHetCol_gen w c hc hw hreg hG hpd hz l)
    (m.tendstoInProb_qform2C_gHetMix_gen w c hc hw hreg hG hpd hz hkl)

include hc hw hreg hG hpd hz in
/-- **`ue` at complex `z`**: `ũ_kᵀ G₀'(z) g_l → 0`. General `Ky` bound via
`dotProduct_UtildeCol_self_gen` in place of the aligned `hR`-based bound. -/
theorem tendstoInProb_cformC_UtildeCol_gHetCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.gHetCol w N ω l)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_cformC_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · l) (fun N => m.hasLaw_pairColR hG hpd N l)
    (fun N => m.UtildeCol w N k)
    (Ky := ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
    (fun N => (m.dotProduct_UtildeCol_self_gen w N k).le)

include hc hw hreg hG hpd hz in
/-- **`ue2` at complex `z`**: `ũ_kᵀ G₀'(z)² g_l → 0`. -/
theorem tendstoInProb_cform2C_UtildeCol_gHetCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (m.W0hetR w N ω) z (m.UtildeCol w N k)
      (m.gHetCol w N ω l)‖) 0 := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.W0hetR_eq_Wsig_blockBR hG hpd w, m.gHetCol_eq_gvec w]
  exact tendstoInProb_cform2C_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hc hw hz hd (tendsto_p_div hd hpd)
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · l) (fun N => m.hasLaw_pairColR hG hpd N l)
    (fun N => m.UtildeCol w N k)
    (Ky := ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
    (fun N => (m.dotProduct_UtildeCol_self_gen w N k).le)

end ComplexEE

section NormEvents

variable [NeZero M] (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hw : ∃ i, w i ≠ 0)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)

include hw hreg hG hpd in
/-- `‖g_k‖² → ∑ w_i² c_i` in probability. General copy of `tendstoInProb_gHetCol_norm`
(`RankR/Het/Forms.lean:765`). -/
theorem tendstoInProb_gHetCol_norm_gen (k : Fin r) :
    TendstoInProb μ (fun N ω => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k)
      (∑ i, w i ^ 2 * c i) := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  simp only [m.gHetCol_eq_gvec w]
  exact tendstoInProb_dotProduct_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hw hd
    (m.tendsto_cN_blkStackR hreg) (m.pairColR hG hpd · k) (fun N => m.hasLaw_pairColR hG hpd N k)

variable [∀ N, IsProbabilityMeasure (μ N)]

include hw hreg hG hpd in
/-- The norm event of `g_k`: `‖g_k‖² ≤ ∑ w_i² c_i + 1` with probability `→ 1`. General copy
of `tendsto_measure_gHetCol_norm_le` (`RankR/Het/Forms.lean:833`). -/
theorem tendsto_measure_gHetCol_norm_le_gen (k : Fin r) :
    Tendsto (fun N => μ N {ω | m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k
      ≤ ∑ i, w i ^ 2 * c i + 1}) atTop (𝓝 1) := by
  have h := m.tendstoInProb_gHetCol_norm_gen w c hw hreg hG hpd k 1 one_pos
  refine tendsto_measure_one_of_bad (s := fun N => {ω | 1 ≤
    |m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k - ∑ i, w i ^ 2 * c i|})
    (fun N _ hω => ?_) h
  simp only [Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω ⊢
  rw [le_abs]
  left
  linarith

end NormEvents

section RealEE

variable [∀ N, IsProbabilityMeasure (μ N)] [NeZero M] (m : UnalignedModelR μ M n d r rk)
  (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {b : ℝ} (hb : MPhet.bHet c w ≤ b)
  (hedge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1))
  {x : ℝ} (hx : b < x)

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `ee`, diagonal, at real `x`**: `g_kᵀ G₀'(x) g_k → Ψ(x)`. General copy of
`tendstoInProb_cform_gHetCol_self` (`RankR/Het/Forms.lean:949`). -/
theorem tendstoInProb_cform_gHetCol_self_gen (k : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω k)) (MPhet.Psihet c w x) := by
  have hbx : MPhet.bHet c w < x := lt_of_le_of_lt hb hx
  exact tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (MultiTableModel.PsiC c w) (MPhet.Psihet c w x)
    (fun z hz => m.tendstoInProb_qformC_gHetCol_gen w c hc hw hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PsiC hc hw hbx)

include hc hw hreg hG hpd hb hedge hx in
/-- **Field `ee2`, diagonal, at real `x`**: `g_kᵀ G₀'(x)² g_k → Ψ'(x)`. -/
theorem tendstoInProb_cform2_gHetCol_self_gen (k : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω k)) (MPhet.PsihetDeriv c w x) := by
  have hbx : MPhet.bHet c w < x := lt_of_le_of_lt hb hx
  exact tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω k)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (MultiTableModel.PsiDerivC c w) (MPhet.PsihetDeriv c w x)
    (fun z hz => m.tendstoInProb_qform2C_gHetCol_gen w c hc hw hreg hG hpd hz k) hx
    (MultiTableModel.tendsto_PsiDerivC hc hw hbx)

include hc hw hreg hG hpd hedge hx in
/-- **Field `ee`, off-diagonal, at real `x`**: `g_kᵀ G₀'(x) g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform_gHetCol_ne_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cformC_gHetCol_ne_gen w c hc hw hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hreg hG hpd hedge hx in
/-- **Field `ee2`, off-diagonal, at real `x`**: `g_kᵀ G₀'(x)² g_l → 0` for `k ≠ l`. -/
theorem tendstoInProb_cform2_gHetCol_ne_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.gHetCol w N ω k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N ω => m.gHetCol w N ω k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (sum_wSqC_add_one_nonneg w c hc) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N k k)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cform2C_gHetCol_ne_gen w c hc hw hreg hG hpd hz hkl) hx
    (by simp)

include hc hw hreg hG hpd hedge hx in
/-- **Field `ue` at real `x`**: `ũ_kᵀ G₀'(x) g_l → 0`. -/
theorem tendstoInProb_cform_UtildeCol_gHetCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_wSqSigMat_diag_nonneg w 0 k) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cformC_UtildeCol_gHetCol_gen w c hc hw hreg hG hpd hz k l) hx
    (by simp)

include hc hw hreg hG hpd hedge hx in
/-- **Field `ue2` at real `x`**: `ũ_kᵀ G₀'(x)² g_l → 0`. -/
theorem tendstoInProb_cform2_UtildeCol_gHetCol_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω => cform2 (m.W0hetR w N ω) x (m.UtildeCol w N k)
      (m.gHetCol w N ω l)) 0 :=
  tendstoInProb_cform2_of_complex_scaled (dN := fun N => ∑ i, n i N)
    (fun N ω => m.W0hetR w N ω) (fun N ω => m.isHermitian_W0hetR w N ω)
    (fun N _ => m.UtildeCol w N k) (fun N ω => m.gHetCol w N ω l)
    (fun ε _ N => (m.measurableSet_lamMax_W0hetR_le w N (b + ε)).nullMeasurableSet) hedge
    (m.sum_wSqSigMat_diag_nonneg w 0 k) (sum_wSqC_add_one_nonneg w c hc)
    (fun N _ => (MultiTableModel.measurableSet_const_prop N _).nullMeasurableSet)
    (m.tendsto_measure_UtildeCol_norm_le_gen w k)
    (fun N _ => (measurableSet_le (m.measurable_gHetCol_dot w N l l)
      measurable_const).nullMeasurableSet)
    (m.tendsto_measure_gHetCol_norm_le_gen w c hw hreg hG hpd l) (fun _ => 0) 0
    (fun z hz => by
      simpa using m.tendstoInProb_cform2C_UtildeCol_gHetCol_gen w c hc hw hreg hG hpd hz k l)
    hx (by simp)

end RealEE

/-! ### 7. The Gaussian discharge -/

/-- **`ResolventLimitsSW` on the general Gaussian family.** General analog of
`resolventLimitsHetR_of_gaussian` (`RankR/Het/Forms.lean:1080`): the `uu`, `uu2` fields are
now matrix valued and hold for every `(k, l)` in one call (no case split, since
`tendstoInProb_cform_UtildeCol_gen` already covers `k = l` and `k ≠ l` uniformly); the edge is
a raw hypothesis, exactly as in the aligned theorem (stage E3's analog is not imported here). -/
theorem resolventLimitsSW_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {b : ℝ} (hb : MPhet.bHet c w ≤ b)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1)) :
    m.ResolventLimitsSW w b
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w) := by
  refine ⟨hedge, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k l x hx
    exact m.tendstoInProb_cform_UtildeCol_gen w c hc hw hreg hG hpd hb hedge hx k l
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform_gHetCol_self_gen w c hc hw hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform_gHetCol_ne_gen w c hc hw hreg hG hpd hedge hx hkl
  · intro k l x hx
    exact m.tendstoInProb_cform_UtildeCol_gHetCol_gen w c hc hw hreg hG hpd hedge hx k l
  · intro k l x hx
    exact m.tendstoInProb_cform2_UtildeCol_gen w c hc hw hreg hG hpd hb hedge hx k l
  · intro k l x hx
    rcases eq_or_ne k l with rfl | hkl
    · simp only [if_true]
      exact m.tendstoInProb_cform2_gHetCol_self_gen w c hc hw hreg hG hpd hb hedge hx k
    · simp only [if_neg hkl]
      exact m.tendstoInProb_cform2_gHetCol_ne_gen w c hc hw hreg hG hpd hedge hx hkl
  · intro k l x hx
    exact m.tendstoInProb_cform2_UtildeCol_gHetCol_gen w c hc hw hreg hG hpd hedge hx k l

/-! ### 8. The column Gram limit (follow-up G1b)

General `rk` mirror of the aligned `section Gram` (`RankR/Het/Forms.lean:1123` to `:1258`).
The cross and off-diagonal Gaussian entries are unchanged (no `hR`); the `Ũ_k ⬝ᵥ Ũ_l` term of
`tendstoInProb_dotProduct_QmatHetR_col_gen` is `dotProduct_UtildeCol_gen`'s unconditional
`∑ i, w_i² S_i(k, l)` rather than the aligned `if k = l then ... else 0`, so the final target
carries that sum outside the `if` and only the Gaussian term stays behind it. -/

section Gram

variable [NeZero M] (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hw : ∃ i, w i ≠ 0)
  (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
  {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)

include hreg hG hpd in
/-- **The cross Gram entry** `ũ_k ⬝ᵥ g_l → 0`, by Chebyshev on the pair law. General copy of
`tendstoInProb_UtildeCol_dotProduct_gHetCol` (`RankR/Het/Forms.lean:1133`): its only `hR` use
was `dotProduct_UtildeCol_self`, replaced here by `dotProduct_UtildeCol_self_gen`. -/
theorem tendstoInProb_UtildeCol_dotProduct_gHetCol_gen (k l : Fin r) :
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
      (K := fun N => Scalars.wSqMax w / d N
        * ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k) ?_
      (Filter.Eventually.of_forall fun N _ => ?_)
    · refine Measurable.abs ?_
      exact Finset.measurable_sum _ fun j _ =>
        ((measurable_gvec_apply _ _ j).comp measurable_fst).const_mul _
    · simp only [sq_abs, hid]
      exact centered_id.integrable_sq_sum _
    · have := ((tendsto_const_nhds (x := Scalars.wSqMax w)).div_atTop hdR).mul_const
        (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
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
        _ = Scalars.wSqMax w / d N
            * ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k := by
            rw [← Finset.mul_sum, ← m.dotProduct_UtildeCol_self_gen w N k]
            congr 1
            simp only [dotProduct, sq]
  refine TendstoInProb.of_le (g := fun N ω => |m.UtildeCol w N k ⬝ᵥ
      gvec (tauOf w (blkStack n N)) (d N) (m.pairColR hG hpd N l ω).1|)
    (fun N => Filter.Eventually.of_forall fun ω => ?_) hbase
  rw [sub_zero, m.gHetCol_eq_gvec w]
  exact le_rfl

include hw hreg hG hpd in
/-- `‖(g_k + g_l)/√2‖² → ∑ w_i² c_i` in probability, `k ≠ l`. General copy of
`tendstoInProb_gHetMix_norm` (`RankR/Het/Forms.lean:776`), which never reads `R_i`. -/
theorem tendstoInProb_gHetMix_norm_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω =>
      gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1
        ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
      (∑ i, w i ^ 2 * c i) := by
  have hd : Tendsto d atTop atTop := (hreg 0).2.1
  exact tendstoInProb_dotProduct_gvec (pN := fun N => ∑ i, n i N) (qN := p)
    (dN := d) (blk := fun N => blkStack n N) hw hd
    (m.tendsto_cN_blkStackR hreg) (m.pairMixR hG hpd · k l 1)
    (fun N => m.hasLaw_pairMixR hG hpd N hkl (Or.inl rfl))

include hw hreg hG hpd in
/-- **The off-diagonal Gaussian Gram entry** `g_k ⬝ᵥ g_l → 0`, `k ≠ l`, by real polarization
through the unit mix. General copy of `tendstoInProb_gHetCol_dotProduct_ne`
(`RankR/Het/Forms.lean:1193`). -/
theorem tendstoInProb_gHetCol_dotProduct_ne_gen {k l : Fin r} (hkl : k ≠ l) :
    TendstoInProb μ (fun N ω => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l) 0 := by
  have hu := m.tendstoInProb_gHetMix_norm_gen w c hw hreg hG hpd hkl
  have hk := m.tendstoInProb_gHetCol_norm_gen w c hw hreg hG hpd k
  have hl := m.tendstoInProb_gHetCol_norm_gen w c hw hreg hG hpd l
  have hfun : (fun N (ω : Ω N) => m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l)
      = fun N ω => (2 * (gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1
          ⬝ᵥ gvec (tauOf w (blkStack n N)) (d N) (m.pairMixR hG hpd N k l 1 ω).1)
        - m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω k
        - m.gHetCol w N ω l ⬝ᵥ m.gHetCol w N ω l) * 2⁻¹ := by
    funext N ω
    have h := m.gHetCol_add_eq_gen w hG hpd N ω k l
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

include hw hreg hG hpd in
/-- **The column Gram limit of `Q` at general `R_i`**:
`Q_k ⬝ᵥ Q_l → ∑ w_i² S_i(k, l) + δ_kl ∑ w_i² c_i`. The exact target `G3` consumes. Mirror of
`tendstoInProb_dotProduct_QmatHetR_col` (`RankR/Het/Forms.lean:1225`): the `Ũ_k ⬝ᵥ Ũ_l` term
is now `dotProduct_UtildeCol_gen`'s unconditional sum, so only the Gaussian term stays behind
the `if`. -/
theorem tendstoInProb_dotProduct_QmatHetR_col_gen (k l : Fin r) :
    TendstoInProb μ (fun N ω =>
        (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l
        + if k = l then ∑ i, w i ^ 2 * c i else 0) := by
  have hfun : (fun N (ω : Ω N) =>
        (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      = fun N ω => (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
          + m.UtildeCol w N k ⬝ᵥ m.gHetCol w N ω l
          + m.UtildeCol w N l ⬝ᵥ m.gHetCol w N ω k
          + m.gHetCol w N ω k ⬝ᵥ m.gHetCol w N ω l := by
    funext N ω
    rw [m.QmatHetR_col w N ω k, m.QmatHetR_col w N ω l, add_dotProduct, dotProduct_add,
      dotProduct_add, m.dotProduct_UtildeCol_gen w N k l,
      dotProduct_comm (m.gHetCol w N ω k) (m.UtildeCol w N l)]
    ring
  rw [hfun]
  have h1 := m.tendstoInProb_UtildeCol_dotProduct_gHetCol_gen w c hreg hG hpd k l
  have h2 := m.tendstoInProb_UtildeCol_dotProduct_gHetCol_gen w c hreg hG hpd l k
  rcases eq_or_ne k l with rfl | hkl
  · have h3 := m.tendstoInProb_gHetCol_norm_gen w c hw hreg hG hpd k
    have hcomb := (((TendstoInProb.const μ
        (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)).add h1).add
      h2).add h3
    simp only [if_true]
    refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    ring
  · have h3 := m.tendstoInProb_gHetCol_dotProduct_ne_gen w c hw hreg hG hpd hkl
    have hcomb := (((TendstoInProb.const μ
        (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)).add h1).add
      h2).add h3
    simp only [if_neg hkl]
    refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    ring

end Gram

end UnalignedModelR

end StackedSVD
