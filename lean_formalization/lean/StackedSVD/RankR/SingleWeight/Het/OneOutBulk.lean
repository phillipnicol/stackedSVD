/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.OneOutDet
import StackedSVD.RankR.SingleWeight.Het.Forms
import StackedSVD.RankR.SingleWeight.Het.Outliers
import StackedSVD.RankR.Het.Bulk

/-!
# The edge window at general `R_i` (F18b, unit U4e)

The twin of `tendsto_measure_normSq_specProj_edge_gt_het` (`RankR/Het/Bulk.lean`) with no
alignment hypothesis: the top-`r` spectral window of `X_W X_Wᵀ` that sits within `ε₁` of
the edge carries no mass of a column of `Q`. The input the aligned proof takes from the
per-column profile comes here from `PhiM2 ⪰ 0` and `Psi2 → +∞`. Plan:
`notes/archive/F18b_plan.md`, section 2 (U4e), section 5 finding 2, and section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### Private helpers: the general-`R_i` twins of the aligned inputs -/

/-- Entry `(k, l)` of `∑ i, a_i S_i + t I`. -/
private theorem sum_smul_add_smul_one_apply (a : Fin M → ℝ)
    (S : Fin M → Matrix (Fin r) (Fin r) ℝ) (t : ℝ) (k l : Fin r) :
    ((∑ i, a i • S i) + t • (1 : Matrix (Fin r) (Fin r) ℝ)) k l
      = (∑ i, a i * S i k l) + if k = l then t else 0 := by
  simp [Matrix.add_apply, Matrix.sum_apply, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul,
    mul_ite]

/-- The quadratic form of `∑ i, a_i S_i + t I`; the shape of `SingleWeight.qform_swFmat2`. -/
private theorem qform_sum_smul_one (a : Fin M → ℝ) (S : Fin M → Matrix (Fin r) (Fin r) ℝ)
    (t : ℝ) (y : Fin r → ℝ) :
    y ⬝ᵥ (((∑ i, a i • S i) + t • (1 : Matrix (Fin r) (Fin r) ℝ)) *ᵥ y)
      = (∑ i, a i * (y ⬝ᵥ (S i *ᵥ y))) + t * (y ⬝ᵥ y) := by
  rw [Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_smul, smul_eq_mul, SingleWeight.qform_sum]

/-- Entry `(k, l)` of `F'(z) = ∑ i, w_i² g_i'(z) S_i + Ψ'(z) I`: exactly the limit of
`ResolventLimitsSW.cform2_qcol` at the Gaussian discharge. -/
private theorem swFmat2_apply (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (z : ℝ) (k l : Fin r) :
    SingleWeight.swFmat2 θ R c w z k l
      = (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z * SingleWeight.sigMat θ R i k l)
        + if k = l then MPhet.PsihetDeriv c w z else 0 := by
  simp only [SingleWeight.swFmat2]
  exact sum_smul_add_smul_one_apply (fun i => w i ^ 2 * MPhet.ghetDeriv c w i z)
    (SingleWeight.sigMat θ R) (MPhet.PsihetDeriv c w z) k l

/-- **A point above the edge where `Ψ'` is large.** The `θ`-free replacement of
`HetBulk.exists_z₀_FhetDeriv_gt`, from `MPhet.PsihetDeriv_tendsto_atTop`. -/
private theorem exists_z₀_PsihetDeriv_gt {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (K : ℝ) :
    ∃ z₀ : ℝ, MPhet.bHet c w < z₀ ∧ z₀ < MPhet.bHet c w + 1 ∧
      K < MPhet.PsihetDeriv c w z₀ := by
  have h : ∀ᶠ z in 𝓝[>] MPhet.bHet c w, K < MPhet.PsihetDeriv c w z :=
    (MPhet.PsihetDeriv_tendsto_atTop hc hw).eventually_gt_atTop K
  have hIoo : ∀ᶠ z in 𝓝[>] MPhet.bHet c w, z ∈ Set.Ioo (MPhet.bHet c w) (MPhet.bHet c w + 1) :=
    Ioo_mem_nhdsGT (by linarith)
  obtain ⟨z₀, hz₀K, hz₀⟩ := (h.and hIoo).exists
  exact ⟨z₀, hz₀.1, hz₀.2, hz₀K⟩

/-- **A form lower bound at `Q y` from entrywise closeness of the column forms.** The
general-`R_i` replacement of `EdgeGlueDetR.le_sum_of_close_to_diag`: the reference matrix `A`
need not be diagonal, only its own form is bounded below. -/
private theorem qform2_lower_of_close {D rr : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℝ)
    (Q : Matrix (Fin D) (Fin rr) ℝ) {A : Matrix (Fin rr) (Fin rr) ℝ} {η μ₀ : ℝ} (hη : 0 ≤ η)
    (hA : ∀ y : Fin rr → ℝ, μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ (A *ᵥ y))
    (hclose : ∀ k l, |R4.cform2 W z (fun i => Q i k) (fun i => Q i l) - A k l| ≤ η)
    (y : Fin rr → ℝ) :
    (μ₀ - η * (rr : ℝ)) * (y ⬝ᵥ y) ≤ R4.qform2 W z (Q *ᵥ y) := by
  have hG : y ⬝ᵥ ((Matrix.of fun k l =>
      R4.cform2 W z (fun i => Q i k) (fun i => Q i l)) *ᵥ y) = R4.qform2 W z (Q *ᵥ y) := by
    rw [EdgeGlueDetR.qform2_mulVec_expand]
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun x _ => by ring
  rw [← hG]
  exact HetBulkDet.form_le_of_close_to_form hη hA (fun k l => by simpa using hclose k l) y

/-- `N₀ = ∑ w_i² c_i` is a floor of the general-`R_i` column Gram limit, since the signal
part `∑ w_i² (S_i)_{ll}` is nonnegative. -/
private theorem colGramFloor_le_colGramLimitSW (m : UnalignedModelR μ M n d r rk)
    (w c : Fin M → ℝ) (l : Fin r) : ∑ i, w i ^ 2 * c i ≤ m.colGramLimitSW w c l := by
  have h := m.sum_wSqSigMat_diag_nonneg w 0 l
  unfold colGramLimitSW
  linarith

/-- **The column Gram event vanishes at general `R_i`.** The twin of
`tendsto_measure_colGram_far`: the limit matrix is `∑ i, w_i² S_i + N₀ I`, which is no longer
diagonal, and the entries come from `tendstoInProb_dotProduct_QmatHetR_col_gen`. -/
private theorem tendsto_measure_colGram_far_sw [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (w c : Fin M → ℝ) (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {η : ℝ} (hη : 0 < η) :
    Tendsto (fun N => μ N (⋃ q : Fin r × Fin r,
      {ω | η ≤ |(fun i => m.QmatHetR w N ω i q.1) ⬝ᵥ (fun i => m.QmatHetR w N ω i q.2)
        - ((∑ i, w i ^ 2 • SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i)
            + (∑ i, w i ^ 2 * c i) • (1 : Matrix (Fin r) (Fin r) ℝ)) q.1 q.2|})) atTop (𝓝 0) :=
  tendsto_measure_zero_iUnion fun q =>
    (FormsR.tendstoInProb_congr_limit
      (sum_smul_add_smul_one_apply (fun i => w i ^ 2)
        (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R) (∑ i, w i ^ 2 * c i) q.1 q.2).symm
      (m.tendstoInProb_dotProduct_QmatHetR_col_gen w c hw hreg hG hpd q.1 q.2)) η hη

/-- The eventual size condition `r ≤ min (∑ n_i) (p N)` at general `rk`. -/
private theorem eventually_r_le_min_p_sw [NeZero M] (m : UnalignedModelR μ M n d r rk)
    {c : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (c i)) {p : ℕ → ℕ}
    (hpd : ∀ N, d N = p N + r) :
    ∀ᶠ N in atTop, r ≤ min (∑ i, n i N) (p N) := by
  filter_upwards [m.eventually_r_le_min hreg] with N hN
  rwa [hpd N, Nat.add_sub_cancel] at hN

/-- The eventual size condition `r + 1 ≤ min (∑ n_i) (d N)` at general `rk`. -/
private theorem eventually_r_add_one_le_min_sw [NeZero M] (m : UnalignedModelR μ M n d r rk)
    {c : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    ∀ᶠ N in atTop, r + 1 ≤ min (∑ i, n i N) (d N) := by
  filter_upwards [(hreg ⟨0, NeZero.pos M⟩).1.eventually_ge_atTop (r + 1),
    (hreg ⟨0, NeZero.pos M⟩).2.1.eventually_ge_atTop (r + 1)] with N h1 h2
  exact le_min (h1.trans (Finset.single_le_sum (f := fun i => n i N)
    (fun i _ => Nat.zero_le _) (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))) h2

/-- `∑ (n_i N)⁻¹ → 0` at general `rk`. -/
private theorem tendsto_sum_inv_n_sw (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    Tendsto (fun N => ∑ i, ((n i N : ℝ))⁻¹) atTop (𝓝 0) := by
  have h := tendsto_finsetSum (Finset.univ : Finset (Fin M)) fun i (_ : i ∈ Finset.univ) =>
    ((tendsto_natCast_atTop_atTop (R := ℝ)).comp (hreg i).1).inv_tendsto_atTop
  simpa using h

/-! ### The edge window -/

/-- U4e: the edge window at general `R_i`, the twin of
`tendsto_measure_normSq_specProj_edge_gt_het` with no `hR` and no `alignedRk`. -/
theorem tendsto_measure_normSq_specProj_edge_gt_sw [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w)) (l : Fin r) :
    ∀ δ > 0, ∃ ε₀ > 0, ∀ ε₁, 0 < ε₁ → ε₁ ≤ ε₀ →
      Tendsto (fun N => μ N {ω |
        δ < ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (topEigSet (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
            (isHermitian_mul_transpose_self _) r ∩ Set.Iic (MPhet.bHet c w + ε₁))
          (m.qColHet w N ω l)‖ ^ 2}) atTop (𝓝 0) := by
  classical
  intro δ hδ
  have hw' : ∃ i, w i ≠ 0 := ⟨⟨0, NeZero.pos M⟩, hw _⟩
  have hlim := m.resolventLimitsSW_of_gaussian w c hc hw' hreg hG hpd le_rfl hedge.edge
  -- 1. the scalars
  have hrpos : 0 < r := l.pos
  have hrR : (0 : ℝ) < r := by exact_mod_cast hrpos
  have hyy : ∀ y : Fin r → ℝ, (0 : ℝ) ≤ y ⬝ᵥ y := fun y => by
    simp only [dotProduct]
    exact Finset.sum_nonneg fun i _ => mul_self_nonneg _
  set N₀ : ℝ := ∑ i, w i ^ 2 * c i with hN₀def
  have hN₀ : 0 < N₀ := colGramFloor_pos hc hw'
  set η : ℝ := N₀ / (2 * r) with hηdef
  have hη : 0 < η := by rw [hηdef]; positivity
  set Nl : Fin r → ℝ := fun k => m.colGramLimitSW w c k with hNldef
  have hNl : ∀ k, N₀ ≤ Nl k := fun k => m.colGramFloor_le_colGramLimitSW w c k
  set CC : Fin r → ℝ := fun k => Nl k + η with hCCdef
  have hCC0 : ∀ k, 0 ≤ CC k := fun k => by
    simp only [hCCdef]
    linarith [hNl k]
  set CQ : ℝ := CC l with hCQdef
  have hCQ : 0 < CQ := by
    simp only [hCQdef, hCCdef]
    linarith [hNl l]
  set K : ℝ := 16 * r * CQ / (δ * N₀) with hKdef
  have hK : 0 < K := by rw [hKdef]; positivity
  obtain ⟨z₀, hz₀lo, hz₀hi, hz₀K⟩ := exists_z₀_PsihetDeriv_gt hc hw' K
  -- the two reference matrices: `F'(z₀)` for the resolvent form, `∑ w_i² S_i + N₀ I` for the Gram
  have hA2form : ∀ y : Fin r → ℝ, MPhet.PsihetDeriv c w z₀ * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (SingleWeight.swFmat2 (fun i => (m.tbl i).θ) m.R c w z₀ *ᵥ y) := by
    intro y
    rw [SingleWeight.qform_swFmat2]
    have hnn : 0 ≤ ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z₀
        * (y ⬝ᵥ (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i *ᵥ y)) :=
      Finset.sum_nonneg fun i _ => mul_nonneg
        (mul_nonneg (sq_nonneg _) (MPhet.ghetDeriv_nonneg hc hw' hz₀lo i))
        (SingleWeight.qform_sigMat_nonneg _ _ i y)
    linarith
  have hAGform : ∀ y : Fin r → ℝ, N₀ * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((∑ i, w i ^ 2 • SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i)
          + (∑ i, w i ^ 2 * c i) • (1 : Matrix (Fin r) (Fin r) ℝ)) *ᵥ y) := by
    intro y
    rw [hN₀def, qform_sum_smul_one (fun i => w i ^ 2)
      (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R) (∑ i, w i ^ 2 * c i) y]
    have hnn : 0 ≤ ∑ i, w i ^ 2
        * (y ⬝ᵥ (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i *ᵥ y)) :=
      Finset.sum_nonneg fun i _ => mul_nonneg (sq_nonneg _)
        (SingleWeight.qform_sigMat_nonneg _ _ i y)
    linarith
  set ε₀ : ℝ := (z₀ - MPhet.bHet c w) / 2 with hε₀def
  have hε₀ : 0 < ε₀ := by rw [hε₀def]; linarith
  have hz₀eq : MPhet.bHet c w + 2 * ε₀ = z₀ := by rw [hε₀def]; ring
  set κ₀ : ℝ := ε₀ ^ 2 * K / (4 * r) with hκ₀def
  have hκ₀ : 0 < κ₀ := by rw [hκ₀def]; positivity
  have hbreq : K / 2 - (r : ℝ) * κ₀ / ε₀ ^ 2 = K / 4 := by
    simp only [hκ₀def]
    field_simp
    ring
  refine ⟨ε₀, hε₀, ?_⟩
  intro ε₁ hε₁ hε₁le
  -- 2. the eventual size conditions
  have hevp := m.eventually_r_le_min_p_sw hreg hpd
  have hev1 := m.eventually_r_add_one_le_min_sw hreg
  -- 3. the bad families
  set bad1 : ∀ N, Set (Ω N) := fun N =>
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + ε₀}ᶜ
    with hbad1def
  set bad2 : ∀ N, Set (Ω N) := fun N =>
    {ω | SimpleSpec (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) r}ᶜ with hbad2def
  set bad3 : ∀ N, Set (Ω N) := fun N =>
    {ω | SimpleSpec (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
      (isHermitian_mul_transpose_self _) r}ᶜ with hbad3def
  set bad4 : ∀ N, Set (Ω N) := fun N =>
    {ω | κ₀ ≤ ∑ a ∈ Finset.range r, ∑ k,
      ‖specProjIdx (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) a
        (m.qColTruncHet w N ω k (CC k))‖ ^ 2} with hbad4def
  set bad5 : ∀ N, Set (Ω N) := fun N => ⋃ q : Fin r × Fin r,
    {ω | K / (2 * r) ≤ |R4.cform2 (m.W0hetR w N ω) z₀
      (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i q.2)
      - SingleWeight.swFmat2 (fun i => (m.tbl i).θ) m.R c w z₀ q.1 q.2|} with hbad5def
  set bad6 : ∀ N, Set (Ω N) := fun N => ⋃ q : Fin r × Fin r,
    {ω | η ≤ |(fun i => m.QmatHetR w N ω i q.1) ⬝ᵥ (fun i => m.QmatHetR w N ω i q.2)
      - ((∑ i, w i ^ 2 • SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i)
          + (∑ i, w i ^ 2 * c i) • (1 : Matrix (Fin r) (Fin r) ℝ)) q.1 q.2|} with hbad6def
  set bad : ∀ N, Set (Ω N) := fun N =>
    bad1 N ∪ (bad2 N ∪ (bad3 N ∪ (bad4 N ∪ (bad5 N ∪ bad6 N)))) with hbaddef
  -- 4. each bad family vanishes
  have hv1 : Tendsto (fun N => μ N (bad1 N)) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N _).nullMeasurableSet)
      (hedge.edge ε₀ hε₀)
  have hv2 : Tendsto (fun N => μ N (bad2 N)) atTop (𝓝 0) := by
    refine EdgeGlueR.tendsto_measure_zero_of_eventually_zero ?_
    filter_upwards [hevp] with N hN
    exact ae_iff.mp (m.simpleSpec_ae_W0hetR w hw hG N (hpd N) hN)
  have hv3 : Tendsto (fun N => μ N (bad3 N)) atTop (𝓝 0) := by
    refine EdgeGlueR.tendsto_measure_zero_of_eventually_zero ?_
    filter_upwards [hev1] with N hN
    exact ae_iff.mp (m.simpleSpec_ae_stackXW_mul_transpose w hw hG N hN)
  have hv4 : Tendsto (fun N => μ N (bad4 N)) atTop (𝓝 0) := by
    have hreal : Tendsto (fun N => (r : ℝ) * (∑ k, CC k) * (∑ i, ((n i N : ℝ))⁻¹) / κ₀)
        atTop (𝓝 0) := by
      have h := ((m.tendsto_sum_inv_n_sw hreg).const_mul ((r : ℝ) * (∑ k, CC k))).div_const κ₀
      simpa using h
    have htend : Tendsto (fun N => ENNReal.ofReal
        ((r : ℝ) * (∑ k, CC k) * (∑ i, ((n i N : ℝ))⁻¹) / κ₀)) atTop (𝓝 0) := by
      have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
      rw [ENNReal.ofReal_zero] at h3
      exact h3
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend
      (Filter.Eventually.of_forall fun _ => zero_le) ?_
    filter_upwards [hevp] with N hN
    exact m.measure_kappa_ge_le_het w hw hG N (hpd N) hN hCC0 hκ₀
  have hv5 : Tendsto (fun N => μ N (bad5 N)) atTop (𝓝 0) :=
    tendsto_measure_zero_iUnion fun q =>
      (FormsR.tendstoInProb_congr_limit
        (swFmat2_apply (fun i => (m.tbl i).θ) m.R c w z₀ q.1 q.2).symm
        (hlim.cform2_qcol q.1 q.2 hz₀lo)) (K / (2 * r)) (by positivity)
  have hv6 : Tendsto (fun N => μ N (bad6 N)) atTop (𝓝 0) :=
    m.tendsto_measure_colGram_far_sw w c hw' hreg hG hpd hη
  have hvbad : Tendsto (fun N => μ N (bad N)) atTop (𝓝 0) :=
    tendsto_measure_zero_union hv1 (tendsto_measure_zero_union hv2
      (tendsto_measure_zero_union hv3 (tendsto_measure_zero_union hv4
        (tendsto_measure_zero_union hv5 hv6))))
  -- 5. the deterministic step on the intersection of the good events
  refine OutliersR.tendsto_measure_zero_of_eventually_subset (t := bad) ?_ hvbad
  filter_upwards [hevp, hev1] with N hpN hN1
  intro ω hω
  by_contra hbad
  simp only [hbaddef, hbad1def, hbad2def, hbad3def, hbad4def, hbad5def, hbad6def,
    Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
    not_or, not_exists, not_le, not_lt, not_not] at hbad
  obtain ⟨he1, he2, he3, he4, he5, he6⟩ := hbad
  have hWh := m.isHermitian_W0hetR w N ω
  have hSh : (m.stackXW w N ω * (m.stackXW w N ω)ᵀ).IsHermitian :=
    isHermitian_mul_transpose_self _
  have hrD : r ≤ ∑ i, n i N := hpN.trans (min_le_left _ _)
  have hrne : (r : ℝ) ≠ 0 := ne_of_gt hrR
  -- (i) the `μmin` bound from bad5
  have hclose5 : ∀ k k' : Fin r,
      |R4.cform2 (m.W0hetR w N ω) z₀ (fun i => m.QmatHetR w N ω i k)
          (fun i => m.QmatHetR w N ω i k')
        - SingleWeight.swFmat2 (fun i => (m.tbl i).θ) m.R c w z₀ k k'| ≤ K / (2 * r) :=
    fun k k' => (he5 (k, k')).le
  have hmin : ∀ y : Fin r → ℝ, (K / 2) * (y ⬝ᵥ y)
      ≤ R4.qform2 (m.W0hetR w N ω) z₀ (m.QmatHetR w N ω *ᵥ y) := by
    intro y
    have hstep := qform2_lower_of_close (m.W0hetR w N ω) z₀ (m.QmatHetR w N ω)
      (by positivity : (0 : ℝ) ≤ K / (2 * (r : ℝ))) hA2form hclose5 y
    have hcancel : K / (2 * (r : ℝ)) * (r : ℝ) = K / 2 := by
      field_simp
    rw [hcancel] at hstep
    have hstep2 : (K / 2) * (y ⬝ᵥ y)
        ≤ (MPhet.PsihetDeriv c w z₀ - K / 2) * (y ⬝ᵥ y) :=
      mul_le_mul_of_nonneg_right (by linarith) (hyy y)
    linarith
  -- (ii) the column norms and the κ bound from bad6 and bad4
  have hclose6 : ∀ k k' : Fin r,
      |(fun i => m.QmatHetR w N ω i k) ⬝ᵥ (fun i => m.QmatHetR w N ω i k')
        - ((∑ i, w i ^ 2 • SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i)
            + (∑ i, w i ^ 2 * c i) • (1 : Matrix (Fin r) (Fin r) ℝ)) k k'| ≤ η :=
    fun k k' => (he6 (k, k')).le
  have hcolnorm : ∀ k : Fin r, ‖m.qColHet w N ω k‖ ^ 2 ≤ CC k := by
    intro k
    rw [qColHet, EdgeDetR.norm_toLp_sq]
    have h := hclose6 k k
    rw [sum_smul_add_smul_one_apply (fun i => w i ^ 2)
      (SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R) (∑ i, w i ^ 2 * c i) k k,
      if_pos rfl] at h
    have h2 := (abs_le.mp h).2
    have hNleq : Nl k = (∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k k)
        + ∑ i, w i ^ 2 * c i := by
      simp only [hNldef, colGramLimitSW]
    simp only [hCCdef]
    linarith
  have hcolT : ∀ k : Fin r, m.qColTruncHet w N ω k (CC k) = m.qColHet w N ω k :=
    fun k => m.qColTruncHet_eq w N ω k (hcolnorm k)
  have hκ : ∑ a ∈ Finset.range r, ∑ k,
      ‖specProjIdx (m.W0hetR w N ω) hWh a
        (WithLp.toLp 2 fun i => m.QmatHetR w N ω i k)‖ ^ 2 ≤ κ₀ := by
    have hh := he4
    simp only [hcolT, qColHet] at hh
    exact hh.le
  -- (iii) the `μQ` bound from bad6
  have hμQ : (0 : ℝ) < N₀ / 2 := by positivity
  have hQQ : ∀ y : Fin r → ℝ, (N₀ / 2) * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) *ᵥ y) := by
    intro y
    have hcl : ∀ k k' : Fin r,
        |((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) k k'
          - ((∑ i, w i ^ 2 • SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i)
              + (∑ i, w i ^ 2 * c i) • (1 : Matrix (Fin r) (Fin r) ℝ)) k k'| ≤ η := by
      intro k k'
      have h := hclose6 k k'
      simp only [dotProduct] at h
      simpa only [Matrix.mul_apply, Matrix.transpose_apply] using h
    have hstep := HetBulkDet.form_le_of_close_to_form hη.le hAGform hcl y
    have hcancel : η * (r : ℝ) = N₀ / 2 := by
      rw [hηdef]
      field_simp
    rw [hcancel] at hstep
    linarith
  -- (iv) the assembly
  have hedge' : lamMax (m.W0hetR w N ω) hWh + ε₀ ≤ z₀ := by linarith
  have hτ : MPhet.bHet c w + ε₁ ≤ z₀ := by linarith
  have hbr : 0 < K / 2 - (r : ℝ) * κ₀ / ε₀ ^ 2 := by rw [hbreq]; positivity
  have hdet := HetBulkDet.normSq_specProj_edge_le_mulVec hWh hSh (m.gram_eq_hetR w N ω) he2 he3
    hrD (z₀ := z₀) (ε₀ := ε₀) (τ := MPhet.bHet c w + ε₁) (κ := κ₀) (μmin := K / 2)
    (μQ := N₀ / 2) hε₀ hedge' hτ hκ hmin hbr hμQ hQQ (Pi.single l 1)
  have hsingle : m.QmatHetR w N ω *ᵥ Pi.single l 1 = fun i => m.QmatHetR w N ω i l := by
    rw [Matrix.mulVec_single_one]
    rfl
  rw [hsingle, hbreq] at hdet
  have hcol : (WithLp.toLp 2 fun i => m.QmatHetR w N ω i l) = m.qColHet w N ω l := rfl
  rw [hcol] at hdet
  have hfinal := hdet.trans
    (HetBulk.window_arith_het hrR hCQ hN₀ hδ hKdef (hcolnorm l))
  exact absurd hω (not_lt.mpr hfinal)

end UnalignedModelR

end StackedSVD
