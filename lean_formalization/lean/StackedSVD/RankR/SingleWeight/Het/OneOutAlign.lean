/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.Align
import StackedSVD.RankR.SingleWeight.Het.Outliers
import StackedSVD.RankR.SingleWeight.Het.Frame
import StackedSVD.RankR.SingleWeight.Het.Count

/-!
# The single outlier of a one-component frame (F18b, units U4d.1, U4d.2)

The twins of `tendstoInProb_eigVal_sw` and `tendstoInProb_normSq_specProj_Ioi_tau_sw`
with one detectable component in place of `EigSep`: the top eigenvalue converges to
`swRho c w γ₀`, and the half-line projector of a column of `Q` has the squared limit
`(y k)² / swNu`. The outlier count enters as the hypothesis `hcount`. Plan:
`notes/archive/F18b_plan.md`, section 2 (U4d) and section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 0. Helpers for the one-component frame -/

/-- The tail event of `hcount` (U4d.1) is a finite intersection of level sets of the sorted
eigenvalues, so it is measurable. The mirror is `measurableSet_count_Ioi_het`
(`RankR/Het/Align.lean:275`); the level sets are measurable by
`UnalignedModelR.measurable_gramHet_eigenvalues₀`. -/
private theorem measurableSet_tail_le (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (τ : ℝ) :
    MeasurableSet {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
      (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ} := by
  have hset : {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}
      = ⋂ k : Fin (Fintype.card (Fin (∑ i, n i N))),
          {ω : Ω N | 1 ≤ (k : ℕ) →
            (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iInter]
  rw [hset]
  refine MeasurableSet.iInter fun k => ?_
  by_cases hk : 1 ≤ (k : ℕ)
  · have he : {ω : Ω N | 1 ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k)
            ⁻¹' Set.Iic τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic, hk, forall_const]
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N k measurableSet_Iic
  · have he : {ω : Ω N | 1 ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}
        = Set.univ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_univ, iff_true]
      intro h
      exact absurd h hk
    rw [he]
    exact MeasurableSet.univ

/-- The count event of `hcount` (U4d.2) is measurable: `Frame.card_filter_eigenvalues`
rewrites the count in the sorted index, and `card = 1` is the finite union over the single
surviving index of a finite intersection of level sets. -/
private theorem measurableSet_card_filter_eq_one (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (N : ℕ) (τ : ℝ) :
    MeasurableSet {ω : Ω N | (Finset.univ.filter fun i =>
      τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i).card = 1} := by
  classical
  have hset : {ω : Ω N | (Finset.univ.filter fun i =>
        τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i).card = 1}
      = ⋃ a : Fin (Fintype.card (Fin (∑ i, n i N))),
          ⋂ k : Fin (Fintype.card (Fin (∑ i, n i N))),
            {ω : Ω N |
              (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
                ↔ k = a)} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iUnion, Set.mem_iInter]
    rw [Frame.card_filter_eigenvalues (isHermitian_mul_transpose_self (m.stackXW w N ω))
      (fun t => τ < t), Finset.card_eq_one]
    constructor
    · rintro ⟨a, ha⟩
      refine ⟨a, fun k => ?_⟩
      have hk := Finset.ext_iff.mp ha k
      simpa using hk
    · rintro ⟨a, ha⟩
      refine ⟨a, ?_⟩
      ext k
      simpa using ha k
  rw [hset]
  refine MeasurableSet.iUnion fun a => MeasurableSet.iInter fun k => ?_
  by_cases hk : k = a
  · have he : {ω : Ω N |
          (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ↔ k = a)}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k)
            ⁻¹' Set.Ioi τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Ioi, hk, iff_true]
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N k measurableSet_Ioi
  · have he : {ω : Ω N |
          (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ↔ k = a)}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k)
            ⁻¹' Set.Iic τ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic, hk, iff_false, not_lt]
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N k measurableSet_Iic

section OneRoot

variable {m : UnalignedModelR μ M n d r rk} {w c : Fin M → ℝ} {γ₀ : ℝ}
  {y : EuclideanSpace ℝ (Fin r)}

/-- At `x = swRho c w γ₀` the first-order frame form of the single column `Q y` tends to
`-(y j)`, by the outlier equation `swFmat_mulVec_eigvec` of unit G0. The `EigSep` twin is
`tendstoInProb_cform_qcol_frame` (`RankR/SingleWeight/Het/Outliers.lean:519`). -/
private theorem tendstoInProb_cform_qcol_one
    (H : m.ResolventLimitsSW w (MPhet.bHet c w)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hγ : Scalars.wSqMax w < γ₀)
    (hthresh : ∑ i, c i * w i ^ 4 / (γ₀ - w i ^ 2) ^ 2 < 1)
    (hyeig : SingleWeight.secMat (fun i => (m.tbl i).θ) m.R w γ₀ *ᵥ WithLp.ofLp y
      = WithLp.ofLp y) (j : Fin r) :
    TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) (SingleWeight.swRho c w γ₀)
        (fun q => m.QmatHetR w N ω q j) (m.QmatHetR w N ω *ᵥ WithLp.ofLp y))
      (-(WithLp.ofLp y j)) := by
  refine FormsR.tendstoInProb_congr_limit ?_
    (m.tendstoInProb_cform_qcol_mulVec w c H j (WithLp.ofLp y)
      (SingleWeight.bHet_lt_swRho hc hw hγ hthresh))
  rw [SingleWeight.swFmat_mulVec_eigvec hc hw hγ hthresh hyeig]
  rfl

/-- At `x = swRho c w γ₀` the second-order frame form of the single column `Q y` tends to
`swNu`, by `qform_swFmat2_eq_swNu` of unit G0. The `EigSep` twin is
`tendstoInProb_cform2_qcol_frame` (`RankR/SingleWeight/Het/Outliers.lean:538`). -/
private theorem tendstoInProb_cform2_qcol_one
    (H : m.ResolventLimitsSW w (MPhet.bHet c w)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hγ : Scalars.wSqMax w < γ₀)
    (hthresh : ∑ i, c i * w i ^ 4 / (γ₀ - w i ^ 2) ^ 2 < 1)
    (hyeig : SingleWeight.secMat (fun i => (m.tbl i).θ) m.R w γ₀ *ᵥ WithLp.ofLp y
      = WithLp.ofLp y) :
    TendstoInProb μ (fun N ω => R4.cform2 (m.W0hetR w N ω) (SingleWeight.swRho c w γ₀)
        (m.QmatHetR w N ω *ᵥ WithLp.ofLp y) (m.QmatHetR w N ω *ᵥ WithLp.ofLp y))
      (SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y) :=
  FormsR.tendstoInProb_congr_limit
    (SingleWeight.qform_swFmat2_eq_swNu hc hw hγ hthresh hyeig)
    (m.tendstoInProb_cform2_qcol_mulVec w c H (WithLp.ofLp y) (WithLp.ofLp y)
      (SingleWeight.bHet_lt_swRho hc hw hγ hthresh))

/-- A unit vector of `EuclideanSpace ℝ (Fin r)` forces `1 ≤ r`. -/
private theorem one_le_of_norm_eq_one (hy : ‖y‖ = 1) : 1 ≤ r := by
  by_contra hcon
  push Not at hcon
  have hr0 : r = 0 := by omega
  subst hr0
  rw [EuclideanSpace.norm_eq] at hy
  simp at hy

end OneRoot

/-! ### 1. U4d.1, the surviving outlier -/

/-- U4d.1: the surviving outlier of a one-component frame. -/
theorem tendstoInProb_eigVal_sw_one [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ₀ : ℝ} {y : EuclideanSpace ℝ (Fin r)}
    (hroot : SingleWeight.IsSecularRoot (fun i => (m.tbl i).θ) m.R w γ₀)
    (hthresh : ∑ i, c i * w i ^ 4 / (γ₀ - w i ^ 2) ^ 2 < 1)
    (hy : ‖y‖ = 1)
    (hyeig : SingleWeight.secMat (fun i => (m.tbl i).θ) m.R w γ₀ *ᵥ WithLp.ofLp y
      = WithLp.ofLp y)
    (hcount : ∀ τ, MPhet.bHet c w < τ → τ < SingleWeight.swRho c w γ₀ →
      Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (m.stackXW w N ω)) 0)
      (SingleWeight.swRho c w γ₀) := by
  classical
  have hr1 : 1 ≤ r := one_le_of_norm_eq_one hy
  have H := m.resolventLimitsSW_of_gaussian w c hc hw hreg hG hpd le_rfl hedge.edge
  have hbρ : MPhet.bHet c w < SingleWeight.swRho c w γ₀ :=
    SingleWeight.bHet_lt_swRho hc hw hroot.1 hthresh
  have hν₀ : 0 < SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y :=
    SingleWeight.swNu_pos hc hw hroot.1 hthresh hy hyeig
  -- 1. the one-column frame and its constants
  set Zs : Matrix (Fin r) (Fin 1) ℝ := fun l _ => WithLp.ofLp y l with hZsdef
  set ρs : Fin 1 → ℝ := fun _ => SingleWeight.swRho c w γ₀ with hρsdef
  set νs : Fin 1 → ℝ := fun _ => SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y
    with hνsdef
  have hρsinj : Function.Injective ρs := fun a b _ => Subsingleton.elim a b
  have hνspos : ∀ a, 0 < νs a := fun _ => hν₀
  set Cq : ℝ := ∑ l, (m.colGramLimitSW w c l + 1) with hCqdef
  have hNl0 : ∀ l, 0 ≤ m.colGramLimitSW w c l + 1 := fun l =>
    (m.colGramLimitSW_pos w c hc hw l).le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hNl0 l
  have hCqle : ∀ l, m.colGramLimitSW w c l + 1 ≤ Cq := fun l =>
    Finset.single_le_sum (f := fun l => m.colGramLimitSW w c l + 1) (fun l' _ => hNl0 l')
      (Finset.mem_univ l)
  set Zb : ℝ := ∑ l, |WithLp.ofLp y l| with hZbdef
  have hZb0 : (0 : ℝ) ≤ Zb := Finset.sum_nonneg fun l _ => abs_nonneg _
  have hZbs : ∀ a : Fin 1, ∑ l, |Zs l a| ≤ Zb := fun _ => le_rfl
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νs := OutliersR.resCG_nonneg Cq r νs
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρs νs := OutliersR.gramC_nonneg ρs νs
  intro ε hε
  -- 2. the radius `κ` and the threshold `τ` below the outlier
  set κ : ℝ := min (ε / 2) ((SingleWeight.swRho c w γ₀ - MPhet.bHet c w) / 2) with hκdef
  have hκ0 : 0 < κ := by rw [hκdef]; exact lt_min (by linarith) (by linarith)
  have hκε : κ ≤ ε / 2 := by rw [hκdef]; exact min_le_left _ _
  have hκb : κ ≤ (SingleWeight.swRho c w γ₀ - MPhet.bHet c w) / 2 := by
    rw [hκdef]; exact min_le_right _ _
  set τ : ℝ := SingleWeight.swRho c w γ₀ - κ with hτdef
  have hτb : MPhet.bHet c w < τ := by rw [hτdef]; linarith
  have hτρ : τ < SingleWeight.swRho c w γ₀ := by rw [hτdef]; linarith
  -- 3. the accuracy `η`
  obtain ⟨η, hη0, -, hδmg, -, hδ'half, -, -⟩ :=
    SingleWeight.exists_accuracy_sw (r := r) (s := 1)
      (CR := OutliersR.resCG Cq r νs) (CG := OutliersR.gramC ρs νs) (Zb := Zb)
      (mg := κ / 2) (g := 1) (ξ := 1) hr1 hCR0 hCG0 hZb0 (by linarith) one_pos one_pos
  -- 4. the five bad families
  set ε₀ : ℝ := (τ - MPhet.bHet c w) / 2 with hε₀def
  have hε₀0 : 0 < ε₀ := by rw [hε₀def]; linarith
  have hε₀τ : MPhet.bHet c w + ε₀ < τ := by rw [hε₀def]; linarith
  have hedgeC : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + ε₀)).nullMeasurableSet)
      (H.edge ε₀ hε₀0)
  have hcountC : Tendsto (fun N => μ N ({ω : Ω N |
      ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}ᶜ))
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_tail_le w N τ).nullMeasurableSet) (hcount τ hτb hτρ)
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimitSW w c l|}) atTop (𝓝 0) :=
    fun l => m.tendstoInProb_colGramSW w c hw hreg hG hpd l 1 one_pos
  have hE1T : ∀ j : Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρs 0)
        (fun i => m.QmatHetR w N ω i j) (fun i => (m.QmatHetR w N ω * Zs) i 0)
        - -Zs j 0|}) atTop (𝓝 0) := fun j =>
    tendstoInProb_cform_qcol_one H hc hw hroot.1 hthresh hyeig j η hη0
  have hE2T : Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρs 0)
        (fun i => (m.QmatHetR w N ω * Zs) i 0) (fun i => (m.QmatHetR w N ω * Zs) i 0)
        - νs 0|}) atTop (𝓝 0) :=
    tendstoInProb_cform2_qcol_one H hc hw hroot.1 hthresh hyeig η hη0
  -- 5. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + ε₀}ᶜ
        ∪ (({ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), 1 ≤ (k : ℕ) →
              (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}ᶜ)
          ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l)
                  ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimitSW w c l|})
            ∪ ((⋃ j : Fin r, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρs 0)
                    (fun i => m.QmatHetR w N ω i j)
                    (fun i => (m.QmatHetR w N ω * Zs) i 0) - -Zs j 0|})
              ∪ {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρs 0)
                    (fun i => (m.QmatHetR w N ω * Zs) i 0)
                    (fun i => (m.QmatHetR w N ω * Zs) i 0) - νs 0|})))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hgood1, hgoodC, hgood2, hgood4, hgood5⟩ := hbad
    have hω' : ε ≤ |eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (m.stackXW w N ω)) 0
        - SingleWeight.swRho c w γ₀| := hω
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by linarith
    have hcolb : ∀ l : Fin r, ∑ i, m.QmatHetR w N ω i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      have hsq : ∑ i, m.QmatHetR w N ω i l ^ 2
          = (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) :=
        Finset.sum_congr rfl fun i _ => sq _
      rw [hsq]
      have h1 := (abs_lt.mp (hgood2 l)).2
      linarith
    have hzlt : ∀ a : Fin 1,
        lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) < ρs a := by
      intro a
      have h1 : ρs a = SingleWeight.swRho c w γ₀ := rfl
      rw [h1]
      linarith
    have hE1 : ∀ (l : Fin r) (a : Fin 1), |R4.cform (m.W0hetR w N ω) (ρs a)
        (fun i => m.QmatHetR w N ω i l) (fun i => (m.QmatHetR w N ω * Zs) i a)
        + Zs l a| ≤ η := by
      intro l a
      have ha : a = 0 := Subsingleton.elim a 0
      subst ha
      have h := (hgood4 l).le
      rwa [sub_neg_eq_add] at h
    have hE2 : ∀ a : Fin 1, |R4.cform2 (m.W0hetR w N ω) (ρs a)
        (fun i => (m.QmatHetR w N ω * Zs) i a) (fun i => (m.QmatHetR w N ω * Zs) i a)
        - νs a| ≤ η := by
      intro a
      have ha : a = 0 := Subsingleton.elim a 0
      subst ha
      exact hgood5.le
    have hres := OutliersR.residual_bound_Z (m.isHermitian_W0hetR w N ω)
      (m.gram_eq_hetR w N ω) hνspos hzlt hCq0 hη0.le hcolb hE1 0
    have hgram := OutliersR.gram_bound_Z (m.isHermitian_W0hetR w N ω) hρsinj hνspos hzlt
      hZb0 hη0.le hZbs hE1 hE2 0 0
    have hyhalf : 1 / 2 ≤ ‖OutliersR.yhatv (m.W0hetR w N ω)
        (m.QmatHetR w N ω * Zs) ρs νs 0‖ ^ 2 := by
      rw [if_pos rfl, real_inner_self_eq_norm_sq] at hgram
      have h1 := (abs_le.mp hgram).1
      linarith
    obtain ⟨i, hi⟩ := Frame.exists_eigenvalue_near_two
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) hyhalf hres
    -- the localized eigenvalue sits above `τ`, so its sorted index is `0`
    have hcard0 : (0 : ℕ) < Fintype.card (Fin (∑ i, n i N)) := by
      rw [Fintype.card_fin]; exact m.stack_row_pos N
    have hval : (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
        ((eigIdx (∑ i, n i N)).symm i)
        = (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i := by
      rw [← eigenvalues_eigIdx, Equiv.apply_symm_apply]
    have habs := abs_le.mp hi
    have hiτ : τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i := by
      have h1 : ρs 0 = SingleWeight.swRho c w γ₀ := rfl
      rw [h1] at habs
      rw [hτdef]
      linarith [habs.1]
    have hk0 : (((eigIdx (∑ i, n i N)).symm i : Fin (Fintype.card (Fin (∑ i, n i N)))) : ℕ)
        = 0 := by
      by_contra hne
      have h1 : 1 ≤ (((eigIdx (∑ i, n i N)).symm i :
          Fin (Fintype.card (Fin (∑ i, n i N)))) : ℕ) := Nat.one_le_iff_ne_zero.mpr hne
      have h2 := hgoodC _ h1
      rw [hval] at h2
      linarith
    have heigval : eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (m.stackXW w N ω)) 0
        = (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i := by
      have h4 : (eigIdx (∑ i, n i N)).symm i = ⟨0, hcard0⟩ := Fin.val_injective (by simpa using hk0)
      rw [eigVal_eq _ _ hcard0, ← hval, h4]
    rw [heigval] at hω'
    have h5 : |(isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i
        - SingleWeight.swRho c w γ₀| ≤ κ / 2 := le_trans hi hδmg
    linarith
  · exact tendsto_measure_zero_union hedgeC
      (tendsto_measure_zero_union hcountC
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T) hE2T)))

/-! ### 2. U4d.2, the half-line overlap -/

/-- U4d.2: the half-line overlap of a column of `Q` when exactly one component detaches. -/
theorem tendstoInProb_normSq_specProj_Ioi_tau_sw_one [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ₀ : ℝ} {y : EuclideanSpace ℝ (Fin r)}
    (hroot : SingleWeight.IsSecularRoot (fun i => (m.tbl i).θ) m.R w γ₀)
    (hthresh : ∑ i, c i * w i ^ 4 / (γ₀ - w i ^ 2) ^ 2 < 1)
    (hy : ‖y‖ = 1)
    (hyeig : SingleWeight.secMat (fun i => (m.tbl i).θ) m.R w γ₀ *ᵥ WithLp.ofLp y
      = WithLp.ofLp y)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : MPhet.bHet c w < τ)
    (habove : τ + mg ≤ SingleWeight.swRho c w γ₀)
    (hcount : Tendsto (fun N => μ N {ω | (Finset.univ.filter fun i =>
        τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i).card = 1})
      atTop (𝓝 1))
    (k : Fin r) :
    TendstoInProb μ
      (fun N ω => ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (m.qColHet w N ω k)‖ ^ 2)
      ((WithLp.ofLp y k) ^ 2
        / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y) := by
  classical
  have hr1 : 1 ≤ r := one_le_of_norm_eq_one hy
  have H := m.resolventLimitsSW_of_gaussian w c hc hw hreg hG hpd le_rfl hedge.edge
  have hbρ : MPhet.bHet c w < SingleWeight.swRho c w γ₀ :=
    SingleWeight.bHet_lt_swRho hc hw hroot.1 hthresh
  have hν₀ : 0 < SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y :=
    SingleWeight.swNu_pos hc hw hroot.1 hthresh hy hyeig
  -- 1. the one-column frame and its constants
  set Zs : Matrix (Fin r) (Fin 1) ℝ := fun l _ => WithLp.ofLp y l with hZsdef
  set ρs : Fin 1 → ℝ := fun _ => SingleWeight.swRho c w γ₀ with hρsdef
  set νs : Fin 1 → ℝ := fun _ => SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y
    with hνsdef
  have hρsinj : Function.Injective ρs := fun a b _ => Subsingleton.elim a b
  have hνspos : ∀ a, 0 < νs a := fun _ => hν₀
  have hρsτ : ∀ a : Fin 1, τ + mg ≤ ρs a := fun _ => habove
  set Cq : ℝ := ∑ l, (m.colGramLimitSW w c l + 1) with hCqdef
  have hNl0 : ∀ l, 0 ≤ m.colGramLimitSW w c l + 1 := fun l =>
    (m.colGramLimitSW_pos w c hc hw l).le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hNl0 l
  have hCqle : ∀ l, m.colGramLimitSW w c l + 1 ≤ Cq := fun l =>
    Finset.single_le_sum (f := fun l => m.colGramLimitSW w c l + 1) (fun l' _ => hNl0 l')
      (Finset.mem_univ l)
  set Zb : ℝ := ∑ l, |WithLp.ofLp y l| with hZbdef
  have hZb0 : (0 : ℝ) ≤ Zb := Finset.sum_nonneg fun l _ => abs_nonneg _
  have hZbs : ∀ a : Fin 1, ∑ l, |Zs l a| ≤ Zb := fun _ => le_rfl
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νs := OutliersR.resCG_nonneg Cq r νs
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρs νs := OutliersR.gramC_nonneg ρs νs
  -- 2. the scale of the column `k` and the frame scalars of the test vector
  have hNk : 0 < m.colGramLimitSW w c k := m.colGramLimitSW_pos w c hc hw k
  set sN : ℝ := Real.sqrt (m.colGramLimitSW w c k) with hsNdef
  have hsN : 0 < sN := Real.sqrt_pos.mpr hNk
  have hsN2 : sN ^ 2 = m.colGramLimitSW w c k := Real.sq_sqrt hNk.le
  have hsNi : 0 < sN⁻¹ := inv_pos.mpr hsN
  set Lb : ℝ := sN⁻¹ * Zb with hLbdef
  have hLb0 : (0 : ℝ) ≤ Lb := mul_nonneg hsNi.le hZb0
  set tv : Fin 1 → ℝ := fun _ => sN⁻¹ * -(WithLp.ofLp y k) with htvdef
  have hLbt : ∀ a, |tv a| ≤ Lb := by
    intro a
    have h1 : tv a = sN⁻¹ * -(WithLp.ofLp y k) := rfl
    rw [h1, abs_mul, abs_of_pos hsNi, abs_neg, hLbdef]
    refine mul_le_mul_of_nonneg_left ?_ hsNi.le
    rw [hZbdef]
    exact Finset.single_le_sum (f := fun l => |WithLp.ofLp y l|) (fun _ _ => abs_nonneg _)
      (Finset.mem_univ k)
  set Tv : ℝ := ∑ a : Fin 1, tv a ^ 2 / νs a with hTvdef
  have htarget : m.colGramLimitSW w c k * Tv
      = (WithLp.ofLp y k) ^ 2
        / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y := by
    have hterm : ∀ a : Fin 1, tv a ^ 2 / νs a
        = (m.colGramLimitSW w c k)⁻¹
          * ((WithLp.ofLp y k) ^ 2
            / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y) := by
      intro a
      have h1 : tv a ^ 2 = (m.colGramLimitSW w c k)⁻¹ * (WithLp.ofLp y k) ^ 2 := by
        have h2 : tv a = sN⁻¹ * -(WithLp.ofLp y k) := rfl
        rw [h2, mul_pow, neg_sq, inv_pow, hsN2]
      rw [h1]
      exact mul_div_assoc _ _ _
    rw [hTvdef, Finset.sum_congr rfl fun a (_ : a ∈ Finset.univ) => hterm a]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, one_smul]
    rw [← mul_assoc, mul_inv_cancel₀ hNk.ne', one_mul]
  -- 3. the accuracy, `ε` dependent
  intro ε hε
  have hε2 : (0 : ℝ) < ε / 2 := by linarith
  have hνsum0 : (0 : ℝ) ≤ ∑ a : Fin 1, (νs a)⁻¹ :=
    Finset.sum_nonneg fun a _ => inv_nonneg.mpr (hνspos a).le
  set B : ℝ := 5 * ((1 : ℕ) : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb)
      + OutliersR.resCG Cq r νs / mg) + (1 + 2 * Lb) * ∑ a : Fin 1, (νs a)⁻¹ with hBdef
  have hB0 : (0 : ℝ) ≤ B := by
    have h1 : (0 : ℝ) ≤ OutliersR.gramC ρs νs * (1 + Zb) + OutliersR.resCG Cq r νs / mg :=
      add_nonneg (mul_nonneg hCG0 (by linarith)) (div_nonneg hCR0 hmg.le)
    have h2 : (0 : ℝ) ≤ 5 * ((1 : ℕ) : ℝ) := by norm_num
    have h3 : (0 : ℝ) ≤ (1 + 2 * Lb) * ∑ a : Fin 1, (νs a)⁻¹ :=
      mul_nonneg (by linarith) hνsum0
    rw [hBdef]
    have h4 := mul_nonneg h2 h1
    linarith
  set K : ℝ := (m.colGramLimitSW w c k + 1) * B with hKdef
  have hK0 : (0 : ℝ) ≤ K := mul_nonneg (hNl0 k) hB0
  set η₂ : ℝ := min (m.colGramLimitSW w c k / 2) (ε / 2 / (2 * (|Tv| + 1))) with hη₂def
  have hη₂0 : 0 < η₂ := by
    rw [hη₂def]
    refine lt_min (by linarith) (div_pos hε2 ?_)
    have h1 := abs_nonneg Tv
    linarith
  have hη₂half : η₂ ≤ m.colGramLimitSW w c k / 2 := by rw [hη₂def]; exact min_le_left _ _
  have hη₂T : η₂ * |Tv| ≤ ε / 2 / 2 := by
    have h1 : η₂ ≤ ε / 2 / (2 * (|Tv| + 1)) := by rw [hη₂def]; exact min_le_right _ _
    rw [mul_comm]
    exact le_trans (mul_le_mul_of_nonneg_left h1 (abs_nonneg _))
      (OutliersR.mul_div_le_half (abs_nonneg _) hε2.le)
  obtain ⟨η, hη0, hηξ, -, -, -, -, hsmalla⟩ :=
    SingleWeight.exists_accuracy_sw (r := r) (s := 1)
      (CR := OutliersR.resCG Cq r νs) (CG := OutliersR.gramC ρs νs) (Zb := Zb)
      (mg := mg) (g := 1) (ξ := min 1 (ε / 2 / (2 * (K + 1)))) hr1 hCR0 hCG0 hZb0 hmg
      one_pos (lt_min one_pos (div_pos hε2 (by linarith)))
  have hη1 : η ≤ 1 := le_trans hηξ (min_le_left _ _)
  have hηK : K * η ≤ ε / 2 / 2 :=
    le_trans (mul_le_mul_of_nonneg_left (le_trans hηξ (min_le_right _ _)) hK0)
      (OutliersR.mul_div_le_half hK0 hε2.le)
  have hRHS : 5 * ((1 : ℕ) : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb) * η
        + OutliersR.resCG Cq r νs * η / mg)
      + (1 + 2 * Lb) * (∑ a : Fin 1, (νs a)⁻¹) * η = B * η := by
    rw [hBdef]; ring
  -- 4. the seven bad families
  set ε₀ : ℝ := (τ - MPhet.bHet c w) / 2 with hε₀def
  have hε₀0 : 0 < ε₀ := by rw [hε₀def]; linarith
  have hε₀τ : MPhet.bHet c w + ε₀ < τ := by rw [hε₀def]; linarith
  have hedgeC : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + ε₀)).nullMeasurableSet)
      (H.edge ε₀ hε₀0)
  have hcountC : Tendsto (fun N => μ N ({ω : Ω N | (Finset.univ.filter fun i =>
      τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i).card = 1}ᶜ))
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_card_filter_eq_one w N τ).nullMeasurableSet) hcount
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimitSW w c l|}) atTop (𝓝 0) :=
    fun l => m.tendstoInProb_colGramSW w c hw hreg hG hpd l 1 one_pos
  have hcol2T : Tendsto (fun N => μ N
      {ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)
        - m.colGramLimitSW w c k|}) atTop (𝓝 0) :=
    m.tendstoInProb_colGramSW w c hw hreg hG hpd k η₂ hη₂0
  have hE1T : ∀ j : Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρs 0)
        (fun i => m.QmatHetR w N ω i j) (fun i => (m.QmatHetR w N ω * Zs) i 0)
        - -Zs j 0|}) atTop (𝓝 0) := fun j =>
    tendstoInProb_cform_qcol_one H hc hw hroot.1 hthresh hyeig j η hη0
  have hE2T : Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρs 0)
        (fun i => (m.QmatHetR w N ω * Zs) i 0) (fun i => (m.QmatHetR w N ω * Zs) i 0)
        - νs 0|}) atTop (𝓝 0) :=
    tendstoInProb_cform2_qcol_one H hc hw hroot.1 hthresh hyeig η hη0
  have hE3T : Tendsto (fun N => μ N
      {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q k)
            ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
          * R4.cform (m.W0hetR w N ω) (ρs 0) (fun q => m.QmatHetR w N ω q k)
            (fun i => (m.QmatHetR w N ω * Zs) i 0)
        - tv 0|}) atTop (𝓝 0) := by
    have h2 := tendstoInProb_cform_qcol_one H hc hw hroot.1 hthresh hyeig k
    have h3 : TendstoInProb μ (fun N (ω : Ω N) =>
        (Real.sqrt ((fun q => m.QmatHetR w N ω q k)
          ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹) sN⁻¹ :=
      (m.tendstoInProb_colGramSW w c hw hreg hG hpd k).comp_continuous
        (φ := fun x => (Real.sqrt x)⁻¹) (Real.continuous_sqrt.continuousAt.inv₀ hsN.ne')
    exact (h3.mul h2) η hη0
  -- 5. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + ε₀}ᶜ
        ∪ (({ω : Ω N | (Finset.univ.filter fun i =>
              τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues i).card
                = 1}ᶜ)
          ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l)
                  ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimitSW w c l|})
            ∪ ({ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q k)
                  ⬝ᵥ (fun q => m.QmatHetR w N ω q k) - m.colGramLimitSW w c k|}
              ∪ ((⋃ j : Fin r, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρs 0)
                      (fun i => m.QmatHetR w N ω i j)
                      (fun i => (m.QmatHetR w N ω * Zs) i 0) - -Zs j 0|})
                ∪ ({ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρs 0)
                        (fun i => (m.QmatHetR w N ω * Zs) i 0)
                        (fun i => (m.QmatHetR w N ω * Zs) i 0) - νs 0|}
                  ∪ {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q k)
                          ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
                        * R4.cform (m.W0hetR w N ω) (ρs 0)
                          (fun q => m.QmatHetR w N ω q k)
                          (fun i => (m.QmatHetR w N ω * Zs) i 0)
                      - tv 0|})))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hgood1, hgoodI, hgood2, hgood3, hgood4, hgood5, hgood6⟩ := hbad
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by linarith
    have hcolb : ∀ l : Fin r, ∑ i, m.QmatHetR w N ω i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      have hsq : ∑ i, m.QmatHetR w N ω i l ^ 2
          = (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) :=
        Finset.sum_congr rfl fun i _ => sq _
      rw [hsq]
      have h1 := (abs_lt.mp (hgood2 l)).2
      linarith
    have hgpos : 0 < (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k) := by
      have h1 := (abs_lt.mp hgood3).1
      linarith
    have hgle : (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)
        ≤ m.colGramLimitSW w c k + 1 := by
      have h1 := (abs_lt.mp (hgood2 k)).2
      linarith
    have hE1s : ∀ (l : Fin r) (a : Fin 1), |R4.cform (m.W0hetR w N ω) (ρs a)
        (fun i => m.QmatHetR w N ω i l) (fun i => (m.QmatHetR w N ω * Zs) i a)
        + Zs l a| ≤ η := by
      intro l a
      have ha : a = 0 := Subsingleton.elim a 0
      subst ha
      have h := (hgood4 l).le
      rwa [sub_neg_eq_add] at h
    have hE2s : ∀ a : Fin 1, |R4.cform2 (m.W0hetR w N ω) (ρs a)
        (fun i => (m.QmatHetR w N ω * Zs) i a) (fun i => (m.QmatHetR w N ω * Zs) i a)
        - νs a| ≤ η := by
      intro a
      have ha : a = 0 := Subsingleton.elim a 0
      subst ha
      exact hgood5.le
    have hE3s : ∀ a : Fin 1, |R4.cform (m.W0hetR w N ω) (ρs a)
        ((Real.sqrt ((fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
          • fun q => m.QmatHetR w N ω q k)
        (fun i => (m.QmatHetR w N ω * Zs) i a) - tv a| ≤ η := by
      intro a
      have ha : a = 0 := Subsingleton.elim a 0
      subst ha
      rw [R4.cform_smul_left]
      exact hgood6.le
    have hdet := OutliersR.align_detZ (m.isHermitian_W0hetR w N ω)
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) (m.gram_eq_hetR w N ω) hρsinj hmg
      hlamτ hρsτ hνspos hgoodI hCq0 hcolb hZb0 hZbs (OutliersR.dot_self_normalize hgpos)
      hLbt hη0 hη1 hE1s hE2s hE3s hsmalla
    have hdet' : |‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
          (WithLp.toLp 2 ((Real.sqrt ((fun q => m.QmatHetR w N ω q k)
            ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
              • fun q => m.QmatHetR w N ω q k))‖ ^ 2 - Tv| ≤ B * η :=
      le_trans hdet (le_of_eq hRHS)
    have hclose := OutliersR.abs_mul_sub_le_of_close hgpos hgle hdet' hgood3.le
    have hω' : ε ≤ |‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2
        - (WithLp.ofLp y k) ^ 2
          / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w γ₀ y| := hω
    rw [OutliersR.normSq_specProj_eq_dot_mul _ _ hgpos, ← htarget] at hω'
    have hKe : (m.colGramLimitSW w c k + 1) * (B * η) = K * η := by rw [hKdef]; ring
    rw [hKe] at hclose
    linarith
  · exact tendsto_measure_zero_union hedgeC
      (tendsto_measure_zero_union hcountC
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
          (tendsto_measure_zero_union hcol2T
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
              (tendsto_measure_zero_union hE2T hE3T)))))

end UnalignedModelR

end StackedSVD
