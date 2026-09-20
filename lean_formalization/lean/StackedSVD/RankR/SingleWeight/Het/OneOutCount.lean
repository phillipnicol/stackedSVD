/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.OneOutDet
import StackedSVD.RankR.SingleWeight.Het.Forms
import StackedSVD.RankR.SingleWeight.Het.Count
import StackedSVD.RankR.Het.Edge

/-!
# The outlier count when one of two components stays in the bulk (F18b, unit U4c)

At `r = 2`, when the secular limit of the profile of component `l₁` is positive at `z₀`
above the edge, every sorted eigenvalue of `X_W X_Wᵀ` of index at least `1` is at most `z₀`
with probability tending to one (`tendsto_measure_count_le_one_sw`). Plan:
`notes/archive/F18b_plan.md`, section 2 (U4c) and section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The two indices of `Fin 2` are distinct: the other index of `l` exists. -/
private theorem exists_ne_fin_two (l : Fin 2) : ∃ l₀ : Fin 2, l₀ ≠ l := by
  refine ⟨if l = 0 then 1 else 0, ?_⟩
  by_cases h : l = 0
  · subst h; decide
  · simp only [if_neg h]; exact Ne.symm h

/-- U4c: at `r = 2`, when the profile of component `l₁` has a positive secular limit at `z₀`,
only one eigenvalue of `X_W X_Wᵀ` stays above `z₀`. -/
theorem tendsto_measure_count_le_one_sw [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d 2 rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + 2)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w)) (l₁ : Fin 2) {z₀ : ℝ}
    (hz₀ : MPhet.bHet c w < z₀)
    (hsec : 0 < 1 + (∑ i, w i ^ 2 * MPhet.ghet c w i z₀
          * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l₁ l₁)
        + MPhet.Psihet c w z₀) :
    Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
      (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ z₀})
      atTop (𝓝 1) := by
  classical
  -- the gap `δ` between the edge and `z₀`
  have hδpos : 0 < (z₀ - MPhet.bHet c w) / 2 := by linarith
  have hδz : MPhet.bHet c w + (z₀ - MPhet.bHet c w) / 2 < z₀ := by linarith
  -- bad family 1: the edge fails
  have hB1 : Tendsto (fun N => μ N
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
        ≤ MPhet.bHet c w + (z₀ - MPhet.bHet c w) / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N _).nullMeasurableSet)
      (hedge.edge _ hδpos)
  -- bad family 2: the secular form of column `l₁` is far from its limit
  have H := m.resolventLimitsSW_of_gaussian w c hc ⟨0, hw 0⟩ hreg hG hpd le_rfl hedge.edge
  have hcf : TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) z₀
      (fun q => m.QmatHetR w N ω q l₁) (fun q => m.QmatHetR w N ω q l₁))
      ((∑ i, w i ^ 2 * MPhet.ghet c w i z₀
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l₁ l₁) + MPhet.Psihet c w z₀) :=
    FormsR.tendstoInProb_congr_limit (by simp) (H.cform_qcol l₁ l₁ hz₀)
  have hB2 : Tendsto (fun N => μ N {ω : Ω N |
      (1 + (∑ i, w i ^ 2 * MPhet.ghet c w i z₀
          * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l₁ l₁)
        + MPhet.Psihet c w z₀) / 2
      ≤ |R4.cform (m.W0hetR w N ω) z₀ (fun q => m.QmatHetR w N ω q l₁)
            (fun q => m.QmatHetR w N ω q l₁)
        - ((∑ i, w i ^ 2 * MPhet.ghet c w i z₀
              * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l₁ l₁)
          + MPhet.Psihet c w z₀)|}) atTop (𝓝 0) :=
    hcf _ (by linarith)
  -- the good event holds off the union of the two bad families
  refine tendsto_measure_one_of_bad ?_ (tendsto_measure_zero_union hB1 hB2)
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or, Set.mem_compl_iff, Set.mem_ofPred_eq, not_not] at hcon
  obtain ⟨g1, g2⟩ := hcon
  have hW := m.isHermitian_W0hetR w N ω
  have hlam : lamMax (m.W0hetR w N ω) hW < z₀ := lt_of_le_of_lt g1 hδz
  have hvm : (Matrix.vecMulVec (fun q => m.QmatHetR w N ω q l₁)
      (fun q => m.QmatHetR w N ω q l₁)).IsHermitian := by
    ext i j
    simp [Matrix.vecMulVec_apply, Matrix.conjTranspose_apply, mul_comm]
  have hA' := hW.add hvm
  have hspos : 0 < R4.secular (m.W0hetR w N ω) (fun q => m.QmatHetR w N ω q l₁) z₀ := by
    have h2 := abs_lt.mp (not_le.mp g2)
    change 0 < 1 + R4.cform (m.W0hetR w N ω) z₀ (fun q => m.QmatHetR w N ω q l₁)
      (fun q => m.QmatHetR w N ω q l₁)
    linarith [h2.1]
  have hcap := R6het.lamMax_add_vecMulVec_le_of_secular_pos (m.stack_row_pos N) hW
    (fun q => m.QmatHetR w N ω q l₁) hA' hlam hspos
  obtain ⟨l₀, hne⟩ := exists_ne_fin_two l₁
  refine hω fun k hk => ?_
  exact Frame.eigenvalues₀_le_of_two_cols (isHermitian_mul_transpose_self (m.stackXW w N ω))
    (m.gram_eq_hetR w N ω) hne hA' hcap k hk


end UnalignedModelR

end StackedSVD
