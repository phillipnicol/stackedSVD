/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Tie
import StackedSVD.RankR.SingleWeight.Het.OneOutCount
import StackedSVD.RankR.SingleWeight.Het.OneOutAlign
import StackedSVD.RankR.SingleWeight.Het.OneOutBulk
import StackedSVD.RankR.SingleWeight.Het.Sup

/-!
# `prop:singleweight_suboptimality` under Gaussian noise (F18b, units U3, U4f, U5)

The three weight regimes of the witness assembled: both roots detectable (`hlaw_both`,
through `prop_gen_rank_stacksvd_singleweight_gaussian`), exactly one (`hlaw_one`, through
the one-component outlier and the edge window), and the tie (`hlaw_tie`, in `Tie.lean`).
`hlaw_witness` discharges the hypothesis of `prop_singleweight_suboptimality_of_law`, and
`prop_singleweight_suboptimality_gaussian` is the unconditional statement. Plan:
`notes/archive/F18b_plan.md`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

namespace Witness

/-! ## U3: both roots detectable, untied -/

/-- U3: the law when both roots are detectable and the weights differ. -/
theorem hlaw_both (w : Fin 2 → ℝ) (h0 : 0 < w 0) (h1 : 0 < w 1) {l₀ l₁ : Fin 2}
    (hne : l₀ ≠ l₁) (hlt : w l₁ < w l₀) (hd : DetectableEx (8 / 5) 1 w l₁) :
    TendstoInProb mu (fun N ω => mdl.perfSW w N ω) (swLimitEx (8 / 5) 1 w) := by
  have hwpos : ∀ i, 0 < w i := Fin.forall_fin_two.2 ⟨h0, h1⟩
  have hd0 : DetectableEx (8 / 5) 1 w l₀ :=
    detectableEx_of_le one_pos (by norm_num) w hne (hwpos l₁) hlt.le
  have hsep : SingleWeight.EigSep (r := 2) (rk := fun _ : Fin 2 => 1)
      (fun i => (mdl.tbl i).θ) mdl.R w (fun _ => (1 : ℝ))
      (fun l => gammaEx (8 / 5) w (ordPair l₀ l₁ l))
      (fun l => EuclideanSpace.single (ordPair l₀ l₁ l) 1) :=
    eigSep_of_detectable (θ₀ := (8 / 5 : ℝ)) (c₀ := (1 : ℝ)) (by norm_num) w hwpos hlt hd0 hd
  have h := mdl.prop_gen_rank_stacksvd_singleweight_gaussian w (fun _ => 1)
    (fun _ => one_pos) mdl_regime mdl_joint _ _ hsep
  have hlim : SingleWeight.swLimit (fun i => (mdl.tbl i).θ) mdl.R w (fun _ => (1 : ℝ))
      (fun l => gammaEx (8 / 5) w (ordPair l₀ l₁ l))
      (fun l => EuclideanSpace.single (ordPair l₀ l₁ l) 1) = swLimitEx (8 / 5) 1 w :=
    (swLimitEx_eq_swLimit_ord (8 / 5) 1 w (by norm_num) hwpos hne hd0 hd).symm
  rw [hlim] at h
  exact h

/-! ## U4f: exactly one root detectable -/

/-- Both weights of the witness are nonzero. -/
private theorem hwall_wit {w : Fin 2 → ℝ} (h0 : 0 < w 0) (h1 : 0 < w 1) : ∀ i, w i ≠ 0 :=
  Fin.forall_fin_two.2 ⟨h0.ne', h1.ne'⟩

/-- The side condition `d N = p N + r` of the Gaussian chain, at `p N = N + 1`. -/
private theorem hpd_wit : ∀ N, dd N = N + 1 + 2 := fun _ => rfl

/-- `S_i = R_i Θ_i² R_iᵀ` of the witness, read through `sigMat_ex`. -/
private theorem hsig_wit (i : Fin 2) :
    SingleWeight.sigMat (fun i => (mdl.tbl i).θ) mdl.R i
      = Matrix.of fun a b => if a = i ∧ b = i then ((8 : ℝ) / 5) ^ 2 else 0 :=
  sigMat_ex (8 / 5) i

/-- The column Gram limit of the witness is diagonal: `sigMat_ex` kills the cross terms, so
the limit at `(a, b)` is `w_a² θ_0² + ∑_i w_i²` on the diagonal and `0` off it. -/
private theorem colGramLimit_wit (w : Fin 2 → ℝ) (a b : Fin 2) :
    (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (mdl.tbl i).θ) mdl.R i a b)
        + (if a = b then ∑ i, w i ^ 2 * (fun _ : Fin 2 => (1 : ℝ)) i else 0)
      = if a = b then w a ^ 2 * ((8 : ℝ) / 5) ^ 2 + (w 0 ^ 2 + w 1 ^ 2) else 0 := by
  simp only [hsig_wit, Matrix.of_apply, Fin.sum_univ_two]
  fin_cases a <;> fin_cases b <;> simp

/-- The secular positivity that `tendsto_measure_count_le_one_sw` needs: at an undetectable
root `l₁` the one-spike profile of component `l₁` has `1 + F > 0` above the edge, by
`not_assumption4_of_not_detectable` and `one_add_F_pos_of_not_assumption4`. -/
private theorem hsec_wit (w : Fin 2 → ℝ) (h0 : 0 < w 0) {l₁ : Fin 2}
    (hnd : ¬ DetectableEx (8 / 5) 1 w l₁) {z : ℝ}
    (hz : MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w < z) :
    0 < 1 + (∑ i, w i ^ 2 * MPhet.ghet (fun _ : Fin 2 => (1 : ℝ)) w i z
        * SingleWeight.sigMat (fun i => (mdl.tbl i).θ) mdl.R i l₁ l₁)
      + MPhet.Psihet (fun _ : Fin 2 => (1 : ℝ)) w z := by
  have h4 := not_assumption4_of_not_detectable (θ₀ := (8 / 5 : ℝ)) (c₀ := (1 : ℝ)) w l₁ hnd
  have hF := MPhet.one_add_F_pos_of_not_assumption4 (fun _ : Fin 2 => one_pos) ⟨0, h0.ne'⟩ h4 hz
  have hsum : (∑ i, w i ^ 2 * MPhet.ghet (fun _ : Fin 2 => (1 : ℝ)) w i z
        * SingleWeight.sigMat (fun i => (mdl.tbl i).θ) mdl.R i l₁ l₁)
      = MPhet.Phihet (fun i => if i = l₁ then (8 / 5 : ℝ) else 0)
          (fun _ : Fin 2 => (1 : ℝ)) w z := by
    unfold MPhet.Phihet
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [hsig_wit i]
    by_cases hil : i = l₁
    · subst hil
      simp
      ring
    · simp [hil, Ne.symm hil]
  rw [hsum]
  unfold MPhet.Fhet at hF
  linarith

/-- U4f.1: when root `l₁` is undetectable, the second index of the weighted stack has no
overlap with the spike directions. (Hypothesis `hnd` added 2026-09-06, decision D39: without
it the claim fails in the both-detectable regime, where `align_sw_of_gaussian` gives the
second index the positive limit `swTermEx l₁`.) -/
theorem align_bulk_sw_one (w : Fin 2 → ℝ) (h0 : 0 < w 0) (h1 : 0 < w 1) {l₁ : Fin 2}
    (hnd : ¬ DetectableEx (8 / 5) 1 w l₁) (k : Fin 2) :
    TendstoInProb mu (fun N ω => overlapIdx (mdl.stackXW w N ω) 1 (mdl.colVecG N k)) 0 := by
  classical
  have hwall : ∀ i, w i ≠ 0 := hwall_wit h0 h1
  have hw' : ∃ i, w i ≠ 0 := ⟨0, h0.ne'⟩
  have hcpos : ∀ i : Fin 2, (0 : ℝ) < (fun _ : Fin 2 => (1 : ℝ)) i := fun _ => one_pos
  have hedge := mdl.heteroEdgeR_of_gaussian_tail w (fun _ => 1) hcpos hw' mdl_regime mdl_joint
  intro ε hε
  have hN₀ : (0 : ℝ) < w 0 ^ 2 + w 1 ^ 2 := by positivity
  set N₀ : ℝ := w 0 ^ 2 + w 1 ^ 2 with hN₀def
  have hδ : (0 : ℝ) < ε * N₀ / 4 := by positivity
  -- the edge window at `δ = ε N₀ / 4`
  obtain ⟨ε₀, hε₀, hwin⟩ := mdl.tendsto_measure_normSq_specProj_edge_gt_sw w (fun _ => 1)
    hcpos hwall mdl_regime mdl_joint hpd_wit hedge k (ε * N₀ / 4) hδ
  have hwin' := hwin ε₀ hε₀ le_rfl
  -- the count event at `τ = bHet + ε₀`
  have hcnt := mdl.tendsto_measure_count_le_one_sw w (fun _ => 1) hcpos hwall mdl_regime
    mdl_joint hpd_wit hedge l₁ (z₀ := MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀)
    (by linarith) (hsec_wit w h0 hnd (by linarith))
  have hcntc : Tendsto (fun N => mu N
      {ω | ∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), 1 ≤ (q : ℕ) →
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues₀ q
          ≤ MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (mdl.measurableSet_eigenvalues₀_le_het w N 1 _).nullMeasurableSet) hcnt
  -- the column Gram event
  have hgram : ∀ q : Fin 2 × Fin 2, Tendsto (fun N => mu N
      {ω | N₀ / (2 * 2) ≤ |(fun i => mdl.QmatHetR w N ω i q.1)
            ⬝ᵥ (fun i => mdl.QmatHetR w N ω i q.2)
          - (if q.1 = q.2 then w q.1 ^ 2 * ((8 : ℝ) / 5) ^ 2 + N₀ else 0)|})
      atTop (𝓝 0) := by
    intro q
    refine FormsR.tendstoInProb_congr_limit ?_
      (mdl.tendstoInProb_dotProduct_QmatHetR_col_gen w (fun _ => 1) hw' mdl_regime mdl_joint
        hpd_wit q.1 q.2) _ (by positivity)
    rw [hN₀def]
    exact colGramLimit_wit w q.1 q.2
  -- the assembly
  refine tendsto_measure_zero_of_subset (t := fun N =>
    {ω | ∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), 1 ≤ (q : ℕ) →
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues₀ q
          ≤ MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀}ᶜ
      ∪ ({ω | ε * N₀ / 4 < ‖specProj (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
            (topEigSet (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
              (isHermitian_mul_transpose_self _) 2
              ∩ Set.Iic (MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀))
            (mdl.qColHet w N ω k)‖ ^ 2}
        ∪ ⋃ q : Fin 2 × Fin 2,
            {ω | N₀ / (2 * 2) ≤ |(fun i => mdl.QmatHetR w N ω i q.1)
                  ⬝ᵥ (fun i => mdl.QmatHetR w N ω i q.2)
                - (if q.1 = q.2 then w q.1 ^ 2 * ((8 : ℝ) / 5) ^ 2 + N₀ else 0)|})) ?_
    (tendsto_measure_zero_union hcntc
      (tendsto_measure_zero_union hwin' (tendsto_measure_zero_iUnion hgram)))
  intro N ω hω
  by_contra hbad
  simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
    not_or, not_exists, not_le, not_lt, not_not] at hbad
  obtain ⟨hcω, hwω, hGω⟩ := hbad
  have hω' : ε ≤ |overlapIdx (mdl.stackXW w N ω) 1 (mdl.colVecG N k) - 0| := hω
  rw [sub_zero, abs_of_nonneg (overlapIdx_nonneg _ _ _)] at hω'
  set X := mdl.stackXW w N ω with hXdef
  have hSh : (X * Xᵀ).IsHermitian := isHermitian_mul_transpose_self _
  -- (i) the form bound from the column Gram event
  have hclose : ∀ a b : Fin 2,
      |(fun i => mdl.QmatHetR w N ω i a) ⬝ᵥ (fun i => mdl.QmatHetR w N ω i b)
        - (if a = b then w a ^ 2 * ((8 : ℝ) / 5) ^ 2 + N₀ else 0)| ≤ N₀ / (2 * 2) :=
    fun a b => (hGω (a, b)).le
  have hQQ : ∀ y : Fin 2 → ℝ, (N₀ / 2) * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((mdl.QmatHetR w N ω)ᵀ * mdl.QmatHetR w N ω) *ᵥ y) :=
    HetBulk.form_le_of_colGram_close (by norm_num) (mdl.QmatHetR w N ω) hN₀
      (fun a => by nlinarith [sq_nonneg (w a)]) hclose
  have hrd : (2 : ℕ) ≤ dd N := by change (2 : ℕ) ≤ N + 3; omega
  have hrn : (2 : ℕ) ≤ ∑ i, nn i N := by
    rw [Fin.sum_univ_two]
    change (2 : ℕ) ≤ N + 3 + (N + 3)
    omega
  have hle : N₀ / 2 ≤ eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) 1 :=
    mdl.le_eigVal_stackGramW_of_form w N ω hQQ (by norm_num) hrn hrd
  have hμ₀ : (0 : ℝ) < N₀ / 2 := by linarith
  have hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) 1 :=
    lt_of_lt_of_le hμ₀ hle
  have hamin : (1 : ℕ) < min (∑ i, nn i N) (dd N) :=
    Het.lt_min_of_eigVal_transpose_mul_self_pos X hpos
  -- (ii) the duality
  have hdual := HetBulkDet.overlapIdx_le_of_eigVal_ge X hamin hμ₀ hle (mdl.colVecG N k)
  rw [hXdef, mdl.stackXW_mulVec_colVecG w N ω k] at hdual
  have hcol : (WithLp.toLp 2 fun q => mdl.QmatHetR w N ω q k) = mdl.qColHet w N ω k := rfl
  rw [hcol] at hdual
  -- (iii) the domination by the edge window
  have hle' : ∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), (q : ℕ) = 1 →
      hSh.eigenvalues₀ q ≤ MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀ :=
    fun q hq => hcω q (by omega)
  have hdom := normSq_specProjIdx_le_edge hSh (by norm_num : (1 : ℕ) < 2) hle'
    (mdl.qColHet w N ω k)
  -- (iv) the arithmetic
  have hwω' : ‖specProj (X * Xᵀ) (topEigSet (X * Xᵀ) hSh 2
      ∩ Set.Iic (MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w + ε₀))
      (mdl.qColHet w N ω k)‖ ^ 2 ≤ ε * N₀ / 4 := hwω
  have hfin : overlapIdx X 1 (mdl.colVecG N k) ≤ (ε * N₀ / 4) / (N₀ / 2) :=
    hdual.trans (div_le_div_of_nonneg_right (hdom.trans hwω') hμ₀.le)
  have harith : (ε * N₀ / 4) / (N₀ / 2) = ε / 2 := by
    field_simp
    ring
  rw [harith] at hfin
  linarith

/-- The count `↔` shape at a threshold just below the outlier: the tail event puts every
sorted index at or above `1` below `τ`, and the top eigenvalue sits within `δ` of `ρ`, hence
above `τ = ρ - δ`. Mirror of `Frame.lt_eigenvalues₀_iff_of_card` with the count replaced by
the two events. -/
private theorem count_one_of_good {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    (hp : 0 < p) {τ ρ δ : ℝ} (hτdef : τ = ρ - δ)
    (h1 : ∀ q : Fin (Fintype.card (Fin p)), 1 ≤ (q : ℕ) → hS.eigenvalues₀ q ≤ τ)
    (h2 : |eigVal S hS 0 - ρ| < δ) (q : Fin (Fintype.card (Fin p))) :
    τ < hS.eigenvalues₀ q ↔ (q : ℕ) < 1 := by
  have hp' : 0 < Fintype.card (Fin p) := by rw [Fintype.card_fin]; exact hp
  constructor
  · intro h
    by_contra hq
    exact absurd h (not_lt.mpr (h1 q (by omega)))
  · intro hq
    have hq0 : (q : ℕ) = 0 := by omega
    have heq : eigVal S hS 0 = hS.eigenvalues₀ q := by
      rw [eigVal_eq S hS hp']
      congr 1
      exact Fin.ext hq0.symm
    have hlow := (abs_lt.mp h2).1
    rw [heq] at hlow
    rw [hτdef]
    linarith

/-- With exactly one sorted eigenvalue above `a`, the projector at the sorted index `0` is
the half-line projector. The `r = 1` case of `specProjIdx_eq_specProj_Ioc`
(`LinAlg/SpecWindow.lean:145`) with no upper threshold: nothing sits above index `0`. -/
private theorem specProjIdx_zero_eq_specProj_Ioi {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}
    (hS : S.IsHermitian) {a : ℝ}
    (hlow : ∀ q : Fin (Fintype.card (Fin p)), a < hS.eigenvalues₀ q ↔ (q : ℕ) < 1) :
    specProjIdx S hS 0 = specProj S (Set.Ioi a) := by
  rw [specProjIdx]
  refine Frame.specProj_congr_of_iff hS fun i => ?_
  have hi : hS.eigenvalues i = hS.eigenvalues₀ ((eigIdx p).symm i) := by
    rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
  rw [hi]
  constructor
  · rintro ⟨q, hq, hqk⟩
    rw [← hqk]
    exact (hlow q).mpr (by omega)
  · intro h
    exact ⟨(eigIdx p).symm i, by have := (hlow _).mp h; omega, rfl⟩

/-- U4f.2: the law when exactly one root is detectable. -/
theorem hlaw_one (w : Fin 2 → ℝ) (h0 : 0 < w 0) (h1 : 0 < w 1) {l₀ l₁ : Fin 2}
    (hne : l₀ ≠ l₁) (hlt : w l₁ < w l₀) (hnd : ¬ DetectableEx (8 / 5) 1 w l₁) :
    TendstoInProb mu (fun N ω => mdl.perfSW w N ω) (swLimitEx (8 / 5) 1 w) := by
  classical
  have hwpos : ∀ i, 0 < w i := Fin.forall_fin_two.2 ⟨h0, h1⟩
  have hcases : (l₀ = 0 ∧ l₁ = 1) ∨ (l₀ = 1 ∧ l₁ = 0) := by
    revert hne
    fin_cases l₀ <;> fin_cases l₁ <;> simp
  have hwall : ∀ i, w i ≠ 0 := hwall_wit h0 h1
  have hw' : ∃ i, w i ≠ 0 := ⟨0, h0.ne'⟩
  have hcpos : ∀ i : Fin 2, (0 : ℝ) < (fun _ : Fin 2 => (1 : ℝ)) i := fun _ => one_pos
  have hedge := mdl.heteroEdgeR_of_gaussian_tail w (fun _ => 1) hcpos hw' mdl_regime mdl_joint
  -- 1. the detectable root and its data
  have hd0 : DetectableEx (8 / 5) 1 w l₀ :=
    detectableEx_of_le one_pos (by norm_num) w hne (hwpos l₁) hlt.le
  have hγ : Scalars.wSqMax w < gammaEx (8 / 5) w l₀ := hd0.1
  have hthresh : ∑ i, (fun _ : Fin 2 => (1 : ℝ)) i * w i ^ 4
      / (gammaEx (8 / 5) w l₀ - w i ^ 2) ^ 2 < 1 := hd0.2
  have hroot : SingleWeight.IsSecularRoot (fun i => (mdl.tbl i).θ) mdl.R w
      (gammaEx (8 / 5) w l₀) :=
    ⟨hd0.1, det_one_sub_secMat_zero (8 / 5) w (hwpos l₀) (by norm_num)⟩
  have hy : ‖(EuclideanSpace.single l₀ (1 : ℝ) : EuclideanSpace ℝ (Fin 2))‖ = 1 := by simp
  have hyeig : SingleWeight.secMat (fun i => (mdl.tbl i).θ) mdl.R w (gammaEx (8 / 5) w l₀)
        *ᵥ WithLp.ofLp (EuclideanSpace.single l₀ (1 : ℝ))
      = WithLp.ofLp (EuclideanSpace.single l₀ (1 : ℝ)) :=
    secMat_ex_mulVec_single (8 / 5) w (hwpos l₀) (by norm_num)
  have hρpos : 0 < SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀) :=
    SingleWeight.swRho_pos hcpos hw' hγ hthresh
  have hνpos : 0 < SingleWeight.swNu (fun i => (mdl.tbl i).θ) mdl.R (fun _ : Fin 2 => (1 : ℝ))
      w (gammaEx (8 / 5) w l₀) (EuclideanSpace.single l₀ 1) :=
    SingleWeight.swNu_pos hcpos hw' hγ hthresh hy hyeig
  have hbρ : MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w
      < SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀) :=
    SingleWeight.bHet_lt_swRho hcpos hw' hγ hthresh
  -- 2. the outlier count and the top eigenvalue
  have hcount : ∀ τ, MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w < τ →
      τ < SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀) →
      Tendsto (fun N => mu N {ω | ∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), 1 ≤ (q : ℕ) →
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues₀ q ≤ τ})
        atTop (𝓝 1) := fun τ hτ _ =>
    mdl.tendsto_measure_count_le_one_sw w (fun _ => 1) hcpos hwall mdl_regime mdl_joint
      hpd_wit hedge l₁ hτ (hsec_wit w h0 hnd hτ)
  have hlam := mdl.tendstoInProb_eigVal_sw_one w (fun _ => 1) hcpos hw' mdl_regime mdl_joint
    hpd_wit hedge hroot hthresh hy hyeig hcount
  -- 3. the threshold `τm` just below the outlier
  obtain ⟨δ, hδdef⟩ : ∃ t : ℝ, t = (SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w
      (gammaEx (8 / 5) w l₀) - MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w) / 2 := ⟨_, rfl⟩
  have hδ : 0 < δ := by rw [hδdef]; linarith
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ, t = SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w
      (gammaEx (8 / 5) w l₀) - δ := ⟨_, rfl⟩
  have hτmb : MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w < τm := by
    rw [hτmdef, hδdef]; linarith
  have hτmρ : τm < SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w
      (gammaEx (8 / 5) w l₀) := by rw [hτmdef]; linarith
  have habove : τm + δ ≤ SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w
      (gammaEx (8 / 5) w l₀) := by rw [hτmdef]; linarith
  -- 4. the bad family: the tail event fails, or the top eigenvalue is far from `ρ`
  have hbad : Tendsto (fun N => mu N
      ({ω : Om N | ∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), 1 ≤ (q : ℕ) →
          (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues₀ q ≤ τm}ᶜ
        ∪ {ω : Om N | δ ≤ |eigVal (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
            (isHermitian_mul_transpose_self (mdl.stackXW w N ω)) 0
            - SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀)|}))
      atTop (𝓝 0) :=
    tendsto_measure_zero_union
      (tendsto_measure_compl_zero
        (fun N => (mdl.measurableSet_eigenvalues₀_le_het w N 1 τm).nullMeasurableSet)
        (hcount τm hτmb hτmρ))
      (hlam δ hδ)
  -- 5. the count in the card form that unit U4d.2 takes
  have hcardone : Tendsto (fun N => mu N {ω : Om N | (Finset.univ.filter fun i =>
      τm < (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues i).card = 1})
      atTop (𝓝 1) := by
    refine tendsto_measure_one_of_bad (fun N ω hω => ?_) hbad
    by_contra hgood
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_ofPred_eq, not_or, not_le,
      not_not] at hgood
    refine hω ?_
    refine Frame.card_filter_of_count (isHermitian_mul_transpose_self (mdl.stackXW w N ω))
      (mdl.stack_row_pos N) ?_
    exact count_one_of_good (isHermitian_mul_transpose_self (mdl.stackXW w N ω))
      (mdl.stack_row_pos N) hτmdef hgood.1 hgood.2
  -- 6. the deterministic identity on the good event
  have hdet : ∀ (k : Fin 2) (N : ℕ) (ω : Om N),
      (∀ q : Fin (Fintype.card (Fin (∑ i, nn i N))), 1 ≤ (q : ℕ) →
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)).eigenvalues₀ q ≤ τm) →
      |eigVal (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (mdl.stackXW w N ω)) 0
        - SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀)| < δ →
      overlapIdx (mdl.stackXW w N ω) 0 (mdl.colVecG N k)
        = ‖specProj (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ) (Set.Ioi τm)
              (mdl.qColHet w N ω k)‖ ^ 2
          / eigVal (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
              (isHermitian_mul_transpose_self (mdl.stackXW w N ω)) 0 := by
    intro k N ω hg1 hg2
    have hnpos : 0 < ∑ i, nn i N := mdl.stack_row_pos N
    have hcard0 : 0 < Fintype.card (Fin (∑ i, nn i N)) := by rw [Fintype.card_fin]; exact hnpos
    have hlow := count_one_of_good (isHermitian_mul_transpose_self (mdl.stackXW w N ω))
      hnpos hτmdef hg1 hg2
    have hbpos : 0 < MPhet.bHet (fun _ : Fin 2 => (1 : ℝ)) w := MPhet.bHet_pos hcpos hw'
    have hpos' : 0 < eigVal (mdl.stackXW w N ω * (mdl.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)) 0 := by
      rw [eigVal_eq _ _ hcard0]
      have hgt := (hlow ⟨0, hcard0⟩).mpr (by norm_num)
      linarith
    have hmin : (0 : ℕ) < min (∑ i, nn i N) (dd N) :=
      lt_min hnpos (by change 0 < N + 3; omega)
    have hpos : 0 < eigVal ((mdl.stackXW w N ω)ᵀ * mdl.stackXW w N ω)
        (isHermitian_transpose_mul_self (mdl.stackXW w N ω)) 0 := by
      rw [Het.eigVal_gram_comm (mdl.stackXW w N ω) hmin]
      exact hpos'
    have hcol : (WithLp.toLp 2 fun q => mdl.QmatHetR w N ω q k) = mdl.qColHet w N ω k := rfl
    rw [Het.overlapIdx_eq_normSq_specProjIdx_div (mdl.stackXW w N ω) hpos (mdl.colVecG N k),
      mdl.stackXW_mulVec_colVecG w N ω k, hcol,
      specProjIdx_zero_eq_specProj_Ioi
        (isHermitian_mul_transpose_self (mdl.stackXW w N ω)) hlow]
  -- 7. the index-`0` limit at each spike direction
  have hq : ∀ k : Fin 2, TendstoInProb mu
      (fun N ω => overlapIdx (mdl.stackXW w N ω) 0 (mdl.colVecG N k))
      (WithLp.ofLp (EuclideanSpace.single l₀ (1 : ℝ) : EuclideanSpace ℝ (Fin 2)) k ^ 2
        / SingleWeight.swNu (fun i => (mdl.tbl i).θ) mdl.R (fun _ : Fin 2 => (1 : ℝ)) w
            (gammaEx (8 / 5) w l₀) (EuclideanSpace.single l₀ 1)
        / SingleWeight.swRho (fun _ : Fin 2 => (1 : ℝ)) w (gammaEx (8 / 5) w l₀)) := by
    intro k
    have hg := mdl.tendstoInProb_normSq_specProj_Ioi_tau_sw_one w (fun _ => 1) hcpos hw'
      mdl_regime mdl_joint hpd_wit hedge hroot hthresh hy hyeig hδ hτmb habove hcardone k
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ (hg.div hlam hρpos.ne')
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) hbad
    by_contra hgood
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_ofPred_eq, not_or, not_le,
      not_not] at hgood
    exact hω (hdet k N ω hgood.1 hgood.2)
  -- 8. the sum over the two indices and the two directions
  have hb := fun k : Fin 2 => align_bulk_sw_one w h0 h1 hnd k
  have hsum := ((hq 0).add (hq 1)).add ((hb 0).add (hb 1))
  have hfun : (fun N (ω : Om N) => mdl.perfSW w N ω)
      = fun N ω => (overlapIdx (mdl.stackXW w N ω) 0 (mdl.colVecG N 0)
            + overlapIdx (mdl.stackXW w N ω) 0 (mdl.colVecG N 1))
          + (overlapIdx (mdl.stackXW w N ω) 1 (mdl.colVecG N 0)
            + overlapIdx (mdl.stackXW w N ω) 1 (mdl.colVecG N 1)) := by
    funext N ω
    simp only [UnalignedModelR.perfSW, Fin.sum_univ_two, Fin.isValue, Fin.val_zero, Fin.val_one]
  rw [hfun]
  refine FormsR.tendstoInProb_congr_limit ?_ hsum
  -- 9. the limit arithmetic
  have hone : WithLp.ofLp (EuclideanSpace.single l₀ (1 : ℝ) : EuclideanSpace ℝ (Fin 2)) 0 ^ 2
      + WithLp.ofLp (EuclideanSpace.single l₀ (1 : ℝ) : EuclideanSpace ℝ (Fin 2)) 1 ^ 2
      = 1 := by
    fin_cases l₀ <;> simp
  have hterm : SingleWeight.swTerm (fun i => (mdl.tbl i).θ) mdl.R w (fun _ : Fin 2 => (1 : ℝ))
      (gammaEx (8 / 5) w l₀) (EuclideanSpace.single l₀ 1) = swTermEx (8 / 5) 1 w l₀ :=
    swTerm_ex (8 / 5) 1 w (by norm_num) hwpos l₀ hd0
  have hid := SingleWeight.one_div_swRho_mul_swNu (θ := fun i => (mdl.tbl i).θ) (R := mdl.R)
    hcpos hw' hγ hthresh hy hyeig
  have hlimEx : swLimitEx (8 / 5) 1 w = swTermEx (8 / 5) 1 w l₀ := by
    unfold swLimitEx
    rcases hcases with ⟨e0, e1⟩ | ⟨e0, e1⟩
    · subst e0; subst e1; rw [if_pos hd0, if_neg hnd]
    · subst e0; subst e1; rw [if_neg hnd, if_pos hd0]
  have harith : ∀ a b x y : ℝ, a + b = 1 → a / (x * y) + b / (x * y) = 1 / (y * x) := by
    intro a b x y h
    rw [← h]
    ring
  rw [hlimEx, ← hterm, ← hid, add_zero, add_zero, div_div, div_div]
  exact harith _ _ _ _ hone

/-- U5.1: the open hypothesis of `prop_singleweight_suboptimality_of_law`, discharged. -/
theorem hlaw_witness : ∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
    TendstoInProb mu (fun N ω => mdl.perfSW w N ω) (swLimitEx (8 / 5) 1 w) := by
  intro w h0 h1
  rcases lt_trichotomy (w 0) (w 1) with hlt | heq | hgt
  · by_cases hd : DetectableEx (8 / 5) 1 w 0
    · exact hlaw_both w h0 h1 (show (1 : Fin 2) ≠ 0 by decide) hlt hd
    · exact hlaw_one w h0 h1 (show (1 : Fin 2) ≠ 0 by decide) hlt hd
  · have hw : w = fun _ => w 0 := funext (Fin.forall_fin_two.2 ⟨rfl, heq.symm⟩)
    rw [hw]
    exact hlaw_tie h0
  · by_cases hd : DetectableEx (8 / 5) 1 w 1
    · exact hlaw_both w h0 h1 (show (0 : Fin 2) ≠ 1 by decide) hgt hd
    · exact hlaw_one w h0 h1 (show (0 : Fin 2) ≠ 1 by decide) hgt hd

/-- U5.2: **`prop:singleweight_suboptimality`** for Gaussian noise, unconditional. -/
theorem prop_singleweight_suboptimality_gaussian :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : UnalignedModelR μ 2 n d 2 (fun _ => 1)),
      (∀ i j, (m.tbl i).θ j = 8 / 5) ∧ (∀ i, (m.tbl i).Regime 1) ∧
      m.JointGaussianNoise ∧ m.R = Rone (RankR.Example.Rex 0) ∧
      TendstoInProb μ (fun N ω => m.perfRG N ω) (2 * betaSq (8 / 5) 1) ∧
      (∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
        TendstoInProb μ (fun N ω => m.perfSW w N ω) (swLimitEx (8 / 5) 1 w) ∧
        swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1) ∧
      2 * betaSq (8 / 5) 2 < 2 * betaSq (8 / 5) 1 :=
  prop_singleweight_suboptimality_of_law hlaw_witness

end Witness

end SingleWeight

end StackedSVD
