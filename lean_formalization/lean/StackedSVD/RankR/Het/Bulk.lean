/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.BulkDet
import StackedSVD.RankR.Het.Outliers

/-!
# Stage E8b: the heteroscedastic bulk step of `thm:rank_r_stacksvd`

Track E, stage E8b. On the exactly aligned Gaussian family with weights `w`, every column
`l` of the signal frame has overlap tending to `0` with the right singular vector of the
weighted stack `X_w` at any sorted index `a ≥ numSup` (the number of supercritical
components). This is the rank-`r` heteroscedastic twin of the Track C bulk argument
(`RankR/RMT/R6R.lean`, `align_cross_of_gaussian_subcritical_aux`) and of the edge-window
bound `RankRStack.tendsto_measure_normSq_specProj_edge_gt` (`RankR/RMT/EdgeGlueR.lean`).

Paper: the bulk half of `thm:rank_r_stacksvd` (`main_paper.tex:2337`, Appendix E,
`sec:rank_r`), `⟨v̂_a, v_l⟩² → 0` for every index `a` past the outliers. The rank-one mirror
is `thm:stacksvd_weighted` (`:463`).

## Contents

1. **Scalars.** A point `z₀ ∈ (bHet, bHet + 1)` where every column derivative
   `FhetDeriv θ_l c w z₀` exceeds a threshold `K` (`exists_z₀_FhetDeriv_gt`), and the
   arithmetic of the final bound (`window_arith_het`).
2. **The column Gram event.** `Qᵀ Q` within `η` of `diag (N_l)` at every pair, with
   probability tending to `1` (`tendsto_measure_colGram_far`), and the deterministic lower
   bound `(N₀ / 2) I ≼ Qᵀ Q` it gives (`form_le_of_colGram_close`), where
   `N₀ := ∑ w_i² c_i` bounds every `N_l` from below without any dependence on `θ`.
3. **The edge window** (`tendsto_measure_normSq_specProj_edge_gt_het`): the twin of
   `EdgeGlueR.tendsto_measure_normSq_specProj_edge_gt :1012`, on the `n`-side Gram
   `S = X_w X_wᵀ = W₀' + Q Qᵀ` (`gram_eq_hetR`) at the raw column `Q_l = X_w v_l`.
4. **The bulk theorem** (`cross_het_of_gaussian_bulk`): E7's count event dominates the
   index projector by the window projector, the column Gram event gives `λ_a(X_wᵀ X_w) ≥ N₀/2`,
   and the duality moves the bound to the `d` side.

## Deviations from Track C's shape

* The test vector of the window bound is the raw column `qColHet w N ω l` (norm² → `N_l`),
  not a unit vector; the deterministic core is E8bd's `normSq_specProj_edge_le_mulVec`
  (residual `0`), so Track C's bad7 (residual) and bad6 (form at `bulkEdge + 1`) are gone.
  The `μQ` bound comes from the column Gram event instead.
* The column norm bound `CC l := N_l + η` comes from the column Gram event too, not from
  the resolvent form at `z₀` (Track C's `hcolnorm`).
* One `z₀` serves every column: `hmin` of the core is one scalar for the whole `Q`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. Scalars -/

namespace HetBulk

/-- **A common point above the edge where every column derivative is large.** From
`MPhet.FhetDeriv_tendsto_atTop` at each of the `r` profiles and the finite intersection of
the `𝓝[>] bHet` events. -/
theorem exists_z₀_FhetDeriv_gt {M r : ℕ} {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (θ : Fin r → Fin M → ℝ) (K : ℝ) :
    ∃ z₀ : ℝ, MPhet.bHet c w < z₀ ∧ z₀ < MPhet.bHet c w + 1 ∧
      ∀ l, K < MPhet.FhetDeriv (θ l) c w z₀ := by
  have h : ∀ l : Fin r, ∀ᶠ z in 𝓝[>] MPhet.bHet c w, K < MPhet.FhetDeriv (θ l) c w z :=
    fun l => (MPhet.FhetDeriv_tendsto_atTop (θ := θ l) hc hw).eventually_gt_atTop K
  have hall := Filter.eventually_all.2 h
  have hIoo : ∀ᶠ z in 𝓝[>] MPhet.bHet c w, z ∈ Set.Ioo (MPhet.bHet c w) (MPhet.bHet c w + 1) :=
    Ioo_mem_nhdsGT (by linarith)
  obtain ⟨z₀, hz₀K, hz₀⟩ := (hall.and hIoo).exists
  exact ⟨z₀, hz₀.1, hz₀.2, hz₀K⟩

/-- The final arithmetic of the window bound: with `K = 16 r CQ / (δ N₀)`, the core bound
`Qn · 2 r / ((K/4) (N₀/2))` is at most `δ` whenever `Qn ≤ CQ`. The name says "window", not
"edge", because `R3het.edge_arith_het` (`RMT/Het/R3het.lean:480`) is a different statement:
that one squares the Sudakov-Fernique bound `(r + L + κ)² ≤ (L + sc)² + ε`. -/
theorem window_arith_het {r CQ N₀ δ K Qn : ℝ} (hr : 0 < r) (hCQ : 0 < CQ) (hN₀ : 0 < N₀)
    (hδ : 0 < δ) (hK : K = 16 * r * CQ / (δ * N₀)) (hQle : Qn ≤ CQ) :
    Qn * (2 * (r / ((K / 4) * (N₀ / 2)))) ≤ δ := by
  subst hK
  have hfac : 2 * (r / ((16 * r * CQ / (δ * N₀) / 4) * (N₀ / 2))) = δ / CQ := by
    field_simp
    ring
  rw [hfac]
  calc Qn * (δ / CQ) ≤ CQ * (δ / CQ) := by gcongr
    _ = δ := by field_simp

/-! ### 2. The column Gram event -/

/-- **`(N₀ / 2) I ≼ Qᵀ Q` from entrywise closeness to `diag (N_l)`.** With
`η := N₀ / (2 r)` and `N₀ ≤ N_l` for every `l`, `HetBulkDet.form_le_of_close_to_diag`
applies at `μ₀ := N₀ / 2`. -/
theorem form_le_of_colGram_close {D rr : ℕ} (hrr : 0 < rr) (Q : Matrix (Fin D) (Fin rr) ℝ)
    {Nl : Fin rr → ℝ} {N₀ : ℝ} (hN₀ : 0 < N₀) (hNl : ∀ k, N₀ ≤ Nl k)
    (hclose : ∀ k l, |(fun i => Q i k) ⬝ᵥ (fun i => Q i l) - (if k = l then Nl k else 0)|
      ≤ N₀ / (2 * rr)) (y : Fin rr → ℝ) :
    (N₀ / 2) * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y) := by
  have hrrR : (0 : ℝ) < rr := by exact_mod_cast hrr
  have hsymm : (Qᵀ * Q)ᵀ = Qᵀ * Q := by
    rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  have hη : 0 ≤ N₀ / (2 * rr) := by positivity
  have hL : ∀ k, N₀ / 2 + N₀ / (2 * rr) * rr ≤ Nl k := by
    intro k
    have : N₀ / (2 * rr) * rr = N₀ / 2 := by field_simp
    rw [this]
    linarith [hNl k]
  refine HetBulkDet.form_le_of_close_to_diag hsymm hη hL ?_ y
  intro k l
  have h := hclose k l
  simp only [dotProduct] at h
  simpa only [Matrix.mul_apply, Matrix.transpose_apply] using h

end HetBulk

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **The column Gram event vanishes.** The union over the `r²` pairs of the events
`η ≤ |Q_k ⬝ᵥ Q_l - δ_kl N_k|` has measure tending to `0`, from E4's
`tendstoInProb_dotProduct_QmatHetR_col` at every pair. -/
theorem tendsto_measure_colGram_far [NeZero M]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {η : ℝ} (hη : 0 < η) :
    Tendsto (fun N => μ N (⋃ q : Fin r × Fin r,
      {ω | η ≤ |(fun i => m.QmatHetR w N ω i q.1) ⬝ᵥ (fun i => m.QmatHetR w N ω i q.2)
        - (if q.1 = q.2 then m.colGramLimit w c q.1 else 0)|})) atTop (𝓝 0) :=
  tendsto_measure_zero_iUnion fun q =>
    (m.tendstoInProb_dotProduct_QmatHetR_col w c hw hR hreg hG hpd q.1 q.2) η hη

/-- The eventual size condition `r ≤ min (∑ n_i) (p N)`, from `eventually_r_le_min` and
`d N = p N + r`. -/
theorem eventually_r_le_min_p [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r))
    {c : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (c i)) {p : ℕ → ℕ}
    (hpd : ∀ N, d N = p N + r) :
    ∀ᶠ N in atTop, r ≤ min (∑ i, n i N) (p N) := by
  filter_upwards [m.eventually_r_le_min hreg] with N hN
  rwa [hpd N, Nat.add_sub_cancel] at hN

/-- The eventual size condition `r + 1 ≤ min (∑ n_i) (d N)`, from table `0`. -/
theorem eventually_r_add_one_le_min [NeZero M]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    ∀ᶠ N in atTop, r + 1 ≤ min (∑ i, n i N) (d N) := by
  filter_upwards [(hreg ⟨0, NeZero.pos M⟩).1.eventually_ge_atTop (r + 1),
    (hreg ⟨0, NeZero.pos M⟩).2.1.eventually_ge_atTop (r + 1)] with N h1 h2
  exact le_min (h1.trans (Finset.single_le_sum (f := fun i => n i N)
    (fun i _ => Nat.zero_le _) (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))) h2

/-- `∑ (n_i N)⁻¹ → 0`, from `n_i → ∞` for every table. -/
theorem tendsto_sum_inv_n (m : UnalignedModelR μ M n d r (alignedRk M r)) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    Tendsto (fun N => ∑ i, ((n i N : ℝ))⁻¹) atTop (𝓝 0) := by
  have h := tendsto_finsetSum (Finset.univ : Finset (Fin M)) fun i (_ : i ∈ Finset.univ) =>
    ((tendsto_natCast_atTop_atTop (R := ℝ)).comp (hreg i).1).inv_tendsto_atTop
  simpa using h

/-! ### 3. The edge window -/

/-- **The hetero edge window.** For every `δ > 0` there is `ε₀ > 0` such that for every
`0 < ε₁ ≤ ε₀`, the projection of the column `Q_l = X_w v_l` onto the eigenspaces of
`S = X_w X_wᵀ` in the window `topEigSet S r ∩ (-∞, bHet + ε₁]` has squared norm at most `δ`
with probability tending to `1`. The twin of
`RankRStack.tendsto_measure_normSq_specProj_edge_gt` (`RankR/RMT/EdgeGlueR.lean:1012`), with
E8bd's `normSq_specProj_edge_le_mulVec` as the deterministic core. -/
theorem tendsto_measure_normSq_specProj_edge_gt_het [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
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
  have hlim := m.resolventLimitsHetR_of_gaussian w c hc hw' hR hreg hG hpd le_rfl hedge.edge
  -- 1. the scalars
  have hrpos : 0 < r := l.pos
  have hrR : (0 : ℝ) < r := by exact_mod_cast hrpos
  set N₀ : ℝ := ∑ i, w i ^ 2 * c i with hN₀def
  have hN₀ : 0 < N₀ := colGramFloor_pos hc hw'
  set η : ℝ := N₀ / (2 * r) with hηdef
  have hη : 0 < η := by rw [hηdef]; positivity
  set Nl : Fin r → ℝ := fun k => m.colGramLimit w c k with hNldef
  have hNl : ∀ k, N₀ ≤ Nl k := fun k => m.colGramFloor_le_colGramLimit w c k
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
  have hKne : K ≠ 0 := ne_of_gt hK
  obtain ⟨z₀, hz₀lo, hz₀hi, hz₀K⟩ :=
    HetBulk.exists_z₀_FhetDeriv_gt hc hw' (fun k i => m.thetaAligned i k) K
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
  have hevp := m.eventually_r_le_min_p hreg hpd
  have hev1 := m.eventually_r_add_one_le_min hreg
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
      - (if q.1 = q.2 then MPhet.PhihetDeriv (fun i => m.thetaAligned i q.1) c w z₀
          + MPhet.PsihetDeriv c w z₀ else 0)|} with hbad5def
  set bad6 : ∀ N, Set (Ω N) := fun N => ⋃ q : Fin r × Fin r,
    {ω | η ≤ |(fun i => m.QmatHetR w N ω i q.1) ⬝ᵥ (fun i => m.QmatHetR w N ω i q.2)
      - (if q.1 = q.2 then m.colGramLimit w c q.1 else 0)|} with hbad6def
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
      have h := ((m.tendsto_sum_inv_n hreg).const_mul ((r : ℝ) * (∑ k, CC k))).div_const κ₀
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
      hlim.cform2_qcol q.1 q.2 hz₀lo (K / (2 * r)) (by positivity)
  have hv6 : Tendsto (fun N => μ N (bad6 N)) atTop (𝓝 0) :=
    m.tendsto_measure_colGram_far w c hw' hR hreg hG hpd hη
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
  -- (i) the `μmin` bound from bad5
  have hclose5 : ∀ k k' : Fin r,
      |R4.cform2 (m.W0hetR w N ω) z₀ (fun i => m.QmatHetR w N ω i k)
          (fun i => m.QmatHetR w N ω i k')
        - (if k = k' then
            ((MPhet.FhetDeriv (fun i => m.thetaAligned i k) c w z₀ / K - 1) + 1) * K else 0)|
        ≤ K / (2 * r) := by
    intro k k'
    have hkk : ((MPhet.FhetDeriv (fun i => m.thetaAligned i k) c w z₀ / K - 1) + 1) * K
        = MPhet.PhihetDeriv (fun i => m.thetaAligned i k) c w z₀
          + MPhet.PsihetDeriv c w z₀ := by
      rw [sub_add_cancel, div_mul_cancel₀ _ hKne]
      rfl
    rw [hkk]
    exact (he5 (k, k')).le
  have hlam : ∀ k : Fin r,
      0 ≤ MPhet.FhetDeriv (fun i => m.thetaAligned i k) c w z₀ / K - 1 := by
    intro k
    rw [sub_nonneg, le_div_iff₀ hK, one_mul]
    exact (hz₀K k).le
  have hmin : ∀ y : Fin r → ℝ, (K / 2) * (y ⬝ᵥ y)
      ≤ R4.qform2 (m.W0hetR w N ω) z₀ (m.QmatHetR w N ω *ᵥ y) := by
    intro y
    rw [EdgeGlueDetR.qform2_mulVec_expand]
    exact EdgeGlueDetR.le_sum_of_close_to_diag hrpos hlam hK.le hclose5 y
  -- (ii) the column norms and the κ bound from bad6 and bad4
  have hclose6 : ∀ k k' : Fin r,
      |(fun i => m.QmatHetR w N ω i k) ⬝ᵥ (fun i => m.QmatHetR w N ω i k')
        - (if k = k' then Nl k else 0)| ≤ η := fun k k' => (he6 (k, k')).le
  have hcolnorm : ∀ k : Fin r, ‖m.qColHet w N ω k‖ ^ 2 ≤ CC k := by
    intro k
    rw [qColHet, EdgeDetR.norm_toLp_sq]
    have h := hclose6 k k
    rw [if_pos rfl] at h
    have h2 := (abs_le.mp h).2
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
      ≤ y ⬝ᵥ (((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) *ᵥ y) :=
    HetBulk.form_le_of_colGram_close hrpos (m.QmatHetR w N ω) hN₀ hNl hclose6
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

/-! ### 4. The bulk theorem -/

/-- **The bulk step of `thm:rank_r_stacksvd` (`main_paper.tex:2337`), heteroscedastic rank
`r`; the rank-one mirror is `thm:stacksvd_weighted` (`:463`).** For every sorted
index `a ≥ numSup` (the number of supercritical components at the weights `w`) and every
column `l` of the signal frame, the squared overlap of the `a`-th right singular vector of the
weighted stack `X_w` with `v_l` tends to `0` in probability. Chain: E7's count event puts
every eigenvalue of `S = X_w X_wᵀ` at index `numSup` and beyond below `bHet + ε₀`; on it,
`normSq_specProjIdx_le_edge` (`RankR/RMT/R6R.lean:62`) dominates the index projector of
`Q_l = X_w v_l` by the window projector of section 3; the column Gram event gives
`λ_a(X_wᵀ X_w) ≥ N₀ / 2` (`le_eigVal_stackGramW_of_form_aligned`), and
`HetBulkDet.overlapIdx_le_of_eigVal_ge` moves the bound to the `d` side. -/
theorem cross_het_of_gaussian_bulk [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∀ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w)) {a : Fin r}
    (ha : Scalars.numSup m.thetaAligned c w ≤ a) (l : Fin r) :
    TendstoInProb μ (fun N ω => overlapIdx (m.stackXW w N ω) a (m.colVecG N l)) 0 := by
  classical
  intro ε hε
  have hw' : ∃ i, w i ≠ 0 := ⟨⟨0, NeZero.pos M⟩, hw _⟩
  have hrpos : 0 < r := l.pos
  -- the scalars
  set N₀ : ℝ := ∑ i, w i ^ 2 * c i with hN₀def
  have hN₀ : 0 < N₀ := colGramFloor_pos hc hw'
  set η : ℝ := N₀ / (2 * r) with hηdef
  have hη : 0 < η := by rw [hηdef]; positivity
  set Nl : Fin r → ℝ := fun k => m.colGramLimit w c k with hNldef
  have hNl : ∀ k, N₀ ≤ Nl k := fun k => m.colGramFloor_le_colGramLimit w c k
  set μ₀ : ℝ := N₀ / 2 with hμ₀def
  have hμ₀ : 0 < μ₀ := by rw [hμ₀def]; positivity
  set δ : ℝ := ε * μ₀ / 2 with hδdef
  have hδ : 0 < δ := by rw [hδdef]; positivity
  -- the edge window at `δ`
  obtain ⟨ε₀, hε₀, hwin⟩ :=
    m.tendsto_measure_normSq_specProj_edge_gt_het w c hc hw hR hreg hG hpd hedge l δ hδ
  have hwin' := hwin ε₀ hε₀ le_rfl
  -- the count event of E7 at `τ := bHet + ε₀`, split by `Assumption4`
  obtain ⟨t, u, e, hinl, hinr, hu⟩ := OutliersR.exists_sum_equiv_split_pred
    (fun k : Fin r => Scalars.Assumption4 (fun i => m.thetaAligned i k) c w)
  have hu' : u = Scalars.numSup m.thetaAligned c w := by
    rw [hu, Scalars.numSup]
  have hτ : MPhet.bHet c w < MPhet.bHet c w + ε₀ := by linarith
  have hcnt := m.tendsto_measure_eigenvalues₀_le_tau_het w c hc hw' hR hreg hG hpd hedge e hτ
    (fun b h4 => absurd h4 (hinl b))
  have hcntc : Tendsto (fun N => μ N
      {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
          ≤ MPhet.bHet c w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_eigenvalues₀_le_het w N u _).nullMeasurableSet) hcnt
  -- the column Gram event
  have hG1 := m.tendsto_measure_colGram_far w c hw' hR hreg hG hpd hη
  -- the assembly
  refine tendsto_measure_zero_of_subset (t := fun N =>
    {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
          ≤ MPhet.bHet c w + ε₀}ᶜ
      ∪ ({ω | δ < ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (topEigSet (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
            (isHermitian_mul_transpose_self _) r ∩ Set.Iic (MPhet.bHet c w + ε₀))
          (m.qColHet w N ω l)‖ ^ 2}
        ∪ ⋃ q : Fin r × Fin r,
          {ω | η ≤ |(fun i => m.QmatHetR w N ω i q.1) ⬝ᵥ (fun i => m.QmatHetR w N ω i q.2)
            - (if q.1 = q.2 then m.colGramLimit w c q.1 else 0)|})) ?_
    (tendsto_measure_zero_union hcntc (tendsto_measure_zero_union hwin' hG1))
  intro N ω hω
  by_contra hbad
  simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
    not_or, not_exists, not_le, not_lt, not_not] at hbad
  obtain ⟨hcω, hwω, hGω⟩ := hbad
  have hω' : ε ≤ |overlapIdx (m.stackXW w N ω) a (m.colVecG N l) - 0| := hω
  rw [sub_zero, abs_of_nonneg (sq_nonneg _)] at hω'
  set X := m.stackXW w N ω with hXdef
  have hSh : (X * Xᵀ).IsHermitian := isHermitian_mul_transpose_self _
  -- (i) `λ_a(Xᵀ X) ≥ μ₀` from the column Gram event
  have hclose6 : ∀ k k' : Fin r,
      |(fun i => m.QmatHetR w N ω i k) ⬝ᵥ (fun i => m.QmatHetR w N ω i k')
        - (if k = k' then Nl k else 0)| ≤ η := fun k k' => (hGω (k, k')).le
  have hQQ : ∀ y : Fin r → ℝ, μ₀ * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) *ᵥ y) :=
    HetBulk.form_le_of_colGram_close hrpos (m.QmatHetR w N ω) hN₀ hNl hclose6
  have hle : μ₀ ≤ eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a :=
    m.le_eigVal_stackGramW_of_form_aligned w N ω hQQ a.isLt
  have hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a := lt_of_lt_of_le hμ₀ hle
  have hamin : (a : ℕ) < min (∑ i, n i N) (d N) :=
    Het.lt_min_of_eigVal_transpose_mul_self_pos X hpos
  -- (ii) the duality
  have hdual := HetBulkDet.overlapIdx_le_of_eigVal_ge X hamin hμ₀ hle (m.colVecG N l)
  rw [hXdef, m.stackXW_mulVec_colVecG w N ω l] at hdual
  have hcol : (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l) = m.qColHet w N ω l := rfl
  rw [hcol] at hdual
  -- (iii) the domination by the window
  have hle' : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))), (q : ℕ) = a →
      hSh.eigenvalues₀ q ≤ MPhet.bHet c w + ε₀ := by
    intro q hq
    refine hcω q ?_
    rw [hq, hu']
    exact ha
  have hdom := normSq_specProjIdx_le_edge hSh a.isLt hle' (m.qColHet w N ω l)
  -- (iv) the arithmetic
  have hwω' : ‖specProj (X * Xᵀ) (topEigSet (X * Xᵀ) hSh r ∩ Set.Iic (MPhet.bHet c w + ε₀))
      (m.qColHet w N ω l)‖ ^ 2 ≤ δ := hwω
  have h1 : overlapIdx X a (m.colVecG N l) ≤ δ / μ₀ :=
    hdual.trans (div_le_div_of_nonneg_right (hdom.trans hwω') hμ₀.le)
  have h2 : δ / μ₀ = ε / 2 := by
    rw [hδdef]
    field_simp
  rw [h2] at h1
  linarith

end UnalignedModelR

end StackedSVD
