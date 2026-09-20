/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.OutliersG
import StackedSVD.RankR.RMT.EdgeR
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# Task U7, part A2a: the outlier part of the spike overlap in the mixed regime

Task A2a of the D32 rank-r campaign. The rank-1 mirror is
`RankRStack.align_of_gaussian_supercritical` (`RankR/RMT/Outliers.lean`), which needs every
spike supercritical because it runs one frame on all `r` spikes.

Here the spikes split into a subcritical part and a supercritical part. Only the
supercritical ones carry a frame vector, so the frame is indexed by an injection
`f : Fin u → Fin r` and the deterministic core is `OutliersR.align_detG`
(`RankR/RMT/OutliersG.lean`, task G2). The count of eigenvalues above the threshold
`τ = bulkEdge c + ε₁` comes from the edge theorem of task U7a
(`RankRStack.tendsto_measure_eigenvalues₀_le`, `RankR/RMT/EdgeR.lean`) through
`OutliersR.count_eq_of_forms_of_edge`.

The conclusion is about `specProj (gram) (Set.Ioi τ)`, the eigenvalues above `τ`, and not
about `specProjTop`; the edge part of the projector is task A2b (`RankR/RMT/EdgeGlueR.lean`).

Contents:

1. `OutliersR.exists_sum_equiv_split`: the split of `Fin r` into the subcritical and the
   supercritical spikes, as one `Fin t ⊕ Fin u ≃ Fin r`. This is the shape task U7a takes.
2. `RankRStack.measurable_X` and `RankRStack.measurableSet_eigenvalues₀_le`: the edge event
   of task U7a is measurable, so its complement is a null tending family.
3. `RankRStack.tendstoInProb_normSq_specProj_Ioi`: the main statement. Every spike, both
   regimes; the limit is `betaSq (√λ_j) c`, which is `0` at a subcritical spike.
4. `RankRStack.exists_margin_le_rhoSq`: the trivial `∃ ε₀` corollary that produces the
   hypothesis `hε₁'` of the main statement.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace OutliersR

/-- The split of the spikes into a subcritical block and a supercritical block, in the
`Fin t ⊕ Fin u ≃ Fin r` shape that `RankRStack.tendsto_measure_eigenvalues₀_le` takes. The
left summand collects the spikes with `λ_k ^ 2 ≤ c`, the right one those with
`c < λ_k ^ 2`. -/
theorem exists_sum_equiv_split {r : ℕ} (lam : Fin r → ℝ) (c : ℝ) :
    ∃ (t u : ℕ) (e : Fin t ⊕ Fin u ≃ Fin r),
      (∀ a : Fin t, lam (e (Sum.inl a)) ^ 2 ≤ c) ∧
        (∀ b : Fin u, c < lam (e (Sum.inr b)) ^ 2) := by
  classical
  refine ⟨Fintype.card {k : Fin r // ¬ (c < lam k ^ 2)},
    Fintype.card {k : Fin r // c < lam k ^ 2},
    (Equiv.sumCongr (Fintype.equivFin {k : Fin r // ¬ (c < lam k ^ 2)}).symm
        (Fintype.equivFin {k : Fin r // c < lam k ^ 2}).symm).trans
      ((Equiv.sumComm _ _).trans (Equiv.sumCompl fun k : Fin r => c < lam k ^ 2)),
    ?_, ?_⟩
  · intro a
    exact not_lt.mp ((Fintype.equivFin {k : Fin r // ¬ (c < lam k ^ 2)}).symm a).2
  · intro b
    exact ((Fintype.equivFin {k : Fin r // c < lam k ^ 2}).symm b).2

/-- Column norm from one diagonal resolvent form. `align_det` reads the column norms
`‖q_l‖ ^ 2 ≤ 2 ρ_l` off its diagonal accuracies `hE1 l l` at the point `ρ_l`. In the mixed
regime a subcritical column has no `ρ_l`, so `align_detG` takes the column norms as the
hypothesis `hcol` with one constant. This lemma builds that constant at one common point
`z` above the edge. -/
theorem dot_self_le_of_cform_close {p : ℕ} {W : Matrix (Fin p) (Fin p) ℝ}
    (hW : W.IsHermitian) (hpsd : ∀ a, 0 ≤ hW.eigenvalues a) {z : ℝ} (hz0 : 0 ≤ z)
    (hz : lamMax W hW < z) (y : Fin p → ℝ) {A : ℝ}
    (hA : |R4.cform W z y y - A| ≤ 1) :
    ∑ i, y i ^ 2 ≤ z * (-A + 1) := by
  have h1 : y ⬝ᵥ y ≤ z * (-R4.qform W z y) := dot_self_le_mul_neg_qform hW hpsd hz y
  have h2 : R4.qform W z y = R4.cform W z y y := (R4.cform_self W z y).symm
  have h3 : ∑ i, y i ^ 2 = y ⬝ᵥ y := by
    simp [dotProduct, sq]
  have h4 := (abs_le.mp hA).1
  rw [h3]
  refine h1.trans ?_
  rw [h2]
  exact mul_le_mul_of_nonneg_left (by linarith) hz0

end OutliersR

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 1. The outlier part of the spike overlap, both regimes

`measurable_X`, `measurable_gram_eigenvalues₀` and `measurableSet_eigenvalues₀_le` used to
stand here. `RankR/RMT/EdgeGlueR.lean` carried the same two public names, so the second
cleanup pass (2026-09-02) moved one copy of each to `RankR/RMT/EdgeR.lean`, which both files
import. -/


/-- **Task A2a, the mixed-regime outlier overlap.** For a `RankRStack` with Gaussian noise
and aspect ratio `c`, the part of the projection of the spike direction `V q_j` on the
eigenvalues of the Gram matrix above `bulkEdge c + ε₁` converges in probability to
`betaSq (√λ_j) c`. No spike has to be supercritical: a subcritical `j` gets the limit `0`,
which is `betaSq` below the threshold.

The rank-1 mirror is `align_of_gaussian_supercritical` (`RankR/RMT/Outliers.lean`), which
takes `hsup : ∀ j, c < λ_j ^ 2` and concludes about `specProjTop`. Here `hε₁'` asks only that
the threshold `bulkEdge c + ε₁` stays below every supercritical outlier `ρ_k` by the margin
`ε₁`; `exists_margin_le_rhoSq` produces it. The edge part of the projector, and so the step
from `specProj (Set.Ioi ·)` to `specProjTop`, is task A2b. -/
theorem tendstoInProb_normSq_specProj_Ioi [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) (j : Fin r)
    {ε₁ : ℝ} (hε₁ : 0 < ε₁)
    (hε₁' : ∀ k, c < s.coreEig k ^ 2 →
      bulkEdge c + 2 * ε₁ ≤ rhoSq (Real.sqrt (s.coreEig k)) c) :
    TendstoInProb μ
      (fun N ω => ‖specProj (s.gram N ω) (Set.Ioi (bulkEdge c + ε₁)) (s.spikeVec j N)‖ ^ 2)
      (betaSq (Real.sqrt (s.coreEig j)) c) := by
  classical
  obtain ⟨U, hU, hsig, hgram, h⟩ := s.resolventLimitsR_of_gaussian hG hc hdtop hns hn hp
  obtain ⟨t, u, e, hsubE, hsupE⟩ := OutliersR.exists_sum_equiv_split s.coreEig c
  set f : Fin u → Fin r := fun b => e (Sum.inr b) with hfdef
  have hfinj : Function.Injective f := fun a b hab => Sum.inr_injective (e.injective hab)
  have hsupf : ∀ k : Fin u, c < s.coreEig (f k) ^ 2 := hsupE
  have hlam0 : ∀ k, 0 ≤ s.coreEig k := s.coreEig_nonneg
  -- 1. the frame scalars, on the supercritical spikes only
  set ρv : Fin u → ℝ := fun k => rhoSq (Real.sqrt (s.coreEig (f k))) c with hρv
  set νv : Fin u → ℝ := fun k => (s.coreEig (f k) + 1) * MP.mDeriv c (ρv k) with hνv
  have hbρ : ∀ k, bulkEdge c < ρv k := fun k =>
    OutliersR.bulkEdge_lt_rho hc (hlam0 (f k)) (hsupf k)
  have hm1 : ∀ k, (s.coreEig (f k) + 1) * MP.m c (ρv k) = -1 := fun k =>
    OutliersR.mul_m_rho_eq_neg_one hc (hlam0 (f k)) (hsupf k)
  have hνpos : ∀ k, 0 < νv k := fun k => OutliersR.nu_pos hc (hlam0 (f k)) (hsupf k)
  have hτρ : ∀ k, bulkEdge c + ε₁ + ε₁ ≤ ρv k := by
    intro k
    have h1 := hε₁' (f k) (hsupE k)
    have h2 : ρv k = rhoSq (Real.sqrt (s.coreEig (f k))) c := rfl
    rw [h2]
    linarith
  -- 2. the column-norm constant, read at the common point `bulkEdge c + 1`
  have hb0 : (0 : ℝ) < bulkEdge c := MP.bulkEdge_pos hc.le
  have hbz₁ : bulkEdge c < bulkEdge c + 1 := by linarith
  have hz₁0 : (0 : ℝ) ≤ bulkEdge c + 1 := by linarith
  have hmneg : MP.m c (bulkEdge c + 1) < 0 := MP.m_neg hc hbz₁
  set Cqf : Fin r → ℝ :=
    fun l => (bulkEdge c + 1) * (-((s.coreEig l + 1) * MP.m c (bulkEdge c + 1)) + 1)
    with hCqf
  have hCqf0 : ∀ l, 0 ≤ Cqf l := by
    intro l
    have h1 : (0 : ℝ) < s.coreEig l + 1 := by linarith [hlam0 l]
    have h2 : (0 : ℝ) ≤ -((s.coreEig l + 1) * MP.m c (bulkEdge c + 1)) := by nlinarith
    have h3 : Cqf l = (bulkEdge c + 1) * (-((s.coreEig l + 1) * MP.m c (bulkEdge c + 1)) + 1) :=
      rfl
    rw [h3]
    exact mul_nonneg hz₁0 (by linarith)
  set Cq : ℝ := ∑ l, Cqf l with hCqdef
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hCqf0 l
  have hCqle : ∀ l, Cqf l ≤ Cq :=
    fun l => Finset.single_le_sum (f := Cqf) (fun l' _ => hCqf0 l') (Finset.mem_univ l)
  -- 3. the overlap target, one function for both regimes
  set Lval : ℝ :=
    Real.sqrt (s.coreEig j) * MP.m c (rhoSq (Real.sqrt (s.coreEig j)) c) with hLval
  set tv : Fin u → ℝ := fun k => if f k = j then Lval else 0 with htv
  have hLb : ∀ k, |tv k| ≤ |Lval| := by
    intro k
    have h1 : tv k = if f k = j then Lval else 0 := rfl
    rw [h1]
    rcases eq_or_ne (f k) j with hkj | hkj
    · rw [if_pos hkj]
    · rw [if_neg hkj, abs_zero]
      exact abs_nonneg Lval
  -- the scalar identity: the frame sum is `betaSq` in both regimes
  have htarget : ∑ k, tv k ^ 2 / νv k = betaSq (Real.sqrt (s.coreEig j)) c := by
    have hjcase : (∃ k₀ : Fin u, f k₀ = j) ∨ s.coreEig j ^ 2 ≤ c := by
      rcases hsym : e.symm j with a | b
      · right
        have hj : e (Sum.inl a) = j := by rw [← hsym]; exact e.apply_symm_apply j
        have hsa := hsubE a
        rwa [hj] at hsa
      · left
        refine ⟨b, ?_⟩
        have hj : e (Sum.inr b) = j := by rw [← hsym]; exact e.apply_symm_apply j
        exact hj
    rcases hjcase with ⟨k₀, hk₀⟩ | hsubj
    · have hsupj : c < s.coreEig j ^ 2 := by
        have hsk := hsupf k₀
        rwa [hk₀] at hsk
      have hcongr : ∀ k : Fin u, tv k ^ 2 / νv k = if k = k₀ then Lval ^ 2 / νv k₀ else 0 := by
        intro k
        have h1 : tv k = if f k = j then Lval else 0 := rfl
        rcases eq_or_ne k k₀ with rfl | hne
        · rw [h1, if_pos hk₀, if_pos rfl]
        · have hfk : f k ≠ j := fun hcon => hne (hfinj (hcon.trans hk₀.symm))
          rw [h1, if_neg hfk, if_neg hne]
          simp
      rw [Finset.sum_congr rfl fun k _ => hcongr k, Finset.sum_ite_eq']
      have hνk₀ : νv k₀
          = (s.coreEig j + 1) * MP.mDeriv c (rhoSq (Real.sqrt (s.coreEig j)) c) := by
        have h1 : νv k₀ = (s.coreEig (f k₀) + 1) * MP.mDeriv c (ρv k₀) := rfl
        have h2 : ρv k₀ = rhoSq (Real.sqrt (s.coreEig (f k₀))) c := rfl
        rw [h1, h2, hk₀]
      have hgoal := OutliersR.sq_div_nu_eq_betaSq hc (hlam0 j) hsupj
      simp only [Finset.mem_univ, if_true]
      rw [hLval, hνk₀]
      exact hgoal
    · have hno : ∀ k : Fin u, f k ≠ j := by
        intro k hcon
        have h1 := hsupf k
        rw [hcon] at h1
        linarith
      have hzero : ∀ k : Fin u, tv k ^ 2 / νv k = 0 := by
        intro k
        have h1 : tv k = if f k = j then Lval else 0 := rfl
        rw [h1, if_neg (hno k)]
        simp
      have hb4 : Real.sqrt (s.coreEig j) ^ 4 = s.coreEig j ^ 2 :=
        OutliersR.sqrt_pow_four (hlam0 j)
      rw [Finset.sum_congr rfl fun k _ => hzero k]
      simp only [Finset.sum_const_zero, betaSq, hb4, gt_iff_lt]
      rw [if_neg (not_lt.mpr hsubj)]
  -- 4. the unit test vector
  have hvdot : ∀ N, WithLp.ofLp (s.spikeVec j N) ⬝ᵥ WithLp.ofLp (s.spikeVec j N) = 1 := by
    intro N
    rw [← Frame.inner_eq_dot]
    have h1 := s.inner_spikeVec N j j
    rwa [if_pos rfl] at h1
  intro ε hε
  -- 5. the constants and the accuracy `η`
  have hνsum0 : (0 : ℝ) ≤ ∑ k, (νv k)⁻¹ :=
    Finset.sum_nonneg fun k _ => inv_nonneg.mpr (hνpos k).le
  have hL0 : (0 : ℝ) ≤ 1 + 2 * |Lval| := by positivity
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv := OutliersR.resCG_nonneg Cq r νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  have hX0 : (0 : ℝ) ≤ (u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / ε₁) :=
    mul_nonneg (Nat.cast_nonneg u) (add_nonneg hCG0 (div_nonneg hCR0 hε₁.le))
  have hB0 : (0 : ℝ) ≤ 5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / ε₁)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ := by
    have h1 := mul_nonneg hL0 hνsum0
    linarith
  set η : ℝ := min (min 1 (ε₁ / (OutliersR.resCG Cq r νv + 1)))
      (min (1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
          + OutliersR.resCG Cq r νv / ε₁) + 1)))
        (ε / (2 * (5 * ((u : ℝ) * (OutliersR.gramC ρv νv
            + OutliersR.resCG Cq r νv / ε₁))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    refine lt_min (lt_min one_pos (div_pos hε₁ (by linarith))) (lt_min ?_ ?_)
    · exact div_pos one_pos (by linarith)
    · exact div_pos hε (by linarith)
  have hη1 : η ≤ 1 := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_left _ _)
  have hηm : η ≤ ε₁ / (OutliersR.resCG Cq r νv + 1) := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_right _ _)
  have hηX : η ≤ 1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / ε₁) + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hηε : η ≤ ε / (2 * (5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / ε₁)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_right _ _)
  have hδm : OutliersR.resCG Cq r νv * η ≤ ε₁ :=
    OutliersR.mul_le_of_le_div_add_one hCR0 hηm hε₁.le
  have hsmall : (u : ℝ) * (OutliersR.gramC ρv νv * η
      + OutliersR.resCG Cq r νv * η / ε₁) ≤ 1 / 2 := by
    have hXe : (u : ℝ) * (OutliersR.gramC ρv νv * η + OutliersR.resCG Cq r νv * η / ε₁)
        = ((u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / ε₁)) * η := by
      ring
    rw [hXe]
    exact le_trans (mul_le_mul_of_nonneg_left hηX hX0) (OutliersR.mul_inv_le_half hX0)
  have hfinal : 5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
        + OutliersR.resCG Cq r νv * η / ε₁)
      + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η ≤ ε / 2 := by
    have hKe : 5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
          + OutliersR.resCG Cq r νv * η / ε₁)
        + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η
        = (5 * ((u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / ε₁))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹) * η := by
      ring
    rw [hKe]
    exact le_trans (mul_le_mul_of_nonneg_left hηε hB0)
      (OutliersR.mul_div_le_half hB0 hε.le)
  -- 6. the seven bad families
  have hedgeC1 : Tendsto (fun N => μ N {ω | lamMax (s.rankRW0 N ω (U N))
      (s.isHermitian_rankRW0 N ω (U N)) ≤ bulkEdge c + ε₁}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + ε₁)).nullMeasurableSet)
      (h.edge ε₁ hε₁)
  have hedgeC2 : Tendsto (fun N => μ N {ω | lamMax (s.rankRW0 N ω (U N))
      (s.isHermitian_rankRW0 N ω (U N)) ≤ bulkEdge c + 1 / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + 1 / 2)).nullMeasurableSet)
      (h.edge (1 / 2) (by norm_num))
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
        (fun i => s.qmatR N ω (U N) i l) (fun i => s.qmatR N ω (U N) i l)
        - (s.coreEig l + 1) * MP.m c (bulkEdge c + 1)|}) atTop (𝓝 0) := by
    intro l
    have h1 := s.tendstoInProb_cform_qmatR h l l hbz₁
    have h2 : TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
        (fun i => s.qmatR N ω (U N) i l) (fun i => s.qmatR N ω (U N) i l))
        ((s.coreEig l + 1) * MP.m c (bulkEdge c + 1)) :=
      FormsR.tendstoInProb_congr_limit (if_pos rfl) h1
    exact h2 1 one_pos
  have hE1T : ∀ q : Fin r × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
        (fun i => s.qmatR N ω (U N) i q.1) (fun i => s.qmatR N ω (U N) i (f q.2))
        - (if q.1 = f q.2 then -1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := s.tendstoInProb_cform_qmatR h q.1 (f q.2) (hbρ q.2)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
        (fun i => s.qmatR N ω (U N) i q.1) (fun i => s.qmatR N ω (U N) i (f q.2)))
        (if q.1 = f q.2 then -1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 (f q.2) with heq | hne'
      · rw [if_pos heq, if_pos heq]
        exact hm1 q.2
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  have hE2T : ∀ q : Fin u × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
        (fun i => s.qmatR N ω (U N) i (f q.1)) (fun i => s.qmatR N ω (U N) i (f q.2))
        - (if q.1 = q.2 then νv q.1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := s.tendstoInProb_cform2_qmatR h (f q.1) (f q.2) (hbρ q.1)
    have h2 : TendstoInProb μ (fun N ω => R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
        (fun i => s.qmatR N ω (U N) i (f q.1)) (fun i => s.qmatR N ω (U N) i (f q.2)))
        (if q.1 = q.2 then νv q.1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 q.2 with heq | hne'
      · rw [if_pos heq, if_pos (congrArg f heq), ← heq]
      · rw [if_neg hne', if_neg (fun hcon => hne' (hfinj hcon))]
    exact h2 η hη0
  have hE3T : ∀ k : Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv k)
        (WithLp.ofLp (s.spikeVec j N)) (fun i => s.qmatR N ω (U N) i (f k))
        - tv k|}) atTop (𝓝 0) := by
    intro k
    have h1 := s.tendstoInProb_cform_vmat h j (f k) (hbρ k)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) (ρv k)
        (WithLp.ofLp (s.spikeVec j N)) (fun i => s.qmatR N ω (U N) i (f k))) (tv k) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      have hts : tv k = if f k = j then Lval else 0 := rfl
      rcases eq_or_ne (f k) j with hkj | hkj
      · rw [if_pos hkj.symm, hts, if_pos hkj, hLval]
        have h3 : ρv k = rhoSq (Real.sqrt (s.coreEig (f k))) c := rfl
        rw [h3, hkj]
      · rw [if_neg (Ne.symm hkj), hts, if_neg hkj]
    exact h2 η hη0
  have hedgeC7 : Tendsto (fun N => μ N
      {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
        (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε₁}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_eigenvalues₀_le N u (bulkEdge c + ε₁)).nullMeasurableSet)
      (s.tendsto_measure_eigenvalues₀_le hc h hgram e hsubE ε₁ hε₁)
  -- 7. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
          ≤ bulkEdge c + ε₁}ᶜ
        ∪ ({ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
              ≤ bulkEdge c + 1 / 2}ᶜ
          ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
                (fun i => s.qmatR N ω (U N) i l) (fun i => s.qmatR N ω (U N) i l)
                - (s.coreEig l + 1) * MP.m c (bulkEdge c + 1)|})
            ∪ ((⋃ q : Fin r × Fin u, {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
                  (fun i => s.qmatR N ω (U N) i q.1)
                  (fun i => s.qmatR N ω (U N) i (f q.2))
                  - (if q.1 = f q.2 then -1 else 0)|})
              ∪ ((⋃ q : Fin u × Fin u, {ω | η ≤ |R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
                    (fun i => s.qmatR N ω (U N) i (f q.1))
                    (fun i => s.qmatR N ω (U N) i (f q.2))
                    - (if q.1 = q.2 then νv q.1 else 0)|})
                ∪ ((⋃ k : Fin u, {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv k)
                      (WithLp.ofLp (s.spikeVec j N))
                      (fun i => s.qmatR N ω (U N) i (f k)) - tv k|})
                  ∪ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
                      (s.isHermitian_gram N ω).eigenvalues₀ k
                        ≤ bulkEdge c + ε₁}ᶜ)))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hg1, hg2, hg3, hg4, hg5, hg6, hg7⟩ := hbad
    have hlt2 : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        < bulkEdge c + 1 := by linarith
    have hcolb : ∀ l : Fin r, ∑ i, s.qmatR N ω (U N) i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      exact OutliersR.dot_self_le_of_cform_close (s.isHermitian_rankRW0 N ω (U N))
        (fun a => s.eigenvalues_rankRW0_nonneg N ω (U N) a) hz₁0 hlt2
        (fun i => s.qmatR N ω (U N) i l) (hg3 l).le
    have hI := OutliersR.count_eq_of_forms_of_edge (Q := s.qmatR N ω (U N))
      (τ := bulkEdge c + ε₁) (mg := ε₁) (Cq := Cq) (η := η)
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω) (hgram N ω) hfinj hε₁
      hg1 hτρ hνpos hCq0 hη0 hcolb (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le)
      hsmall hδm hg7
    have hdet := OutliersR.align_detG (Q := s.qmatR N ω (U N))
      (τ := bulkEdge c + ε₁) (mg := ε₁) (Cq := Cq) (t := tv) (Lb := |Lval|) (η := η)
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω) (hgram N ω) hfinj hε₁
      hg1 hτρ hνpos hI hCq0 hcolb (hvdot N) hLb hη0 hη1
      (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le) (fun k => (hg6 k).le) hsmall
    rw [WithLp.toLp_ofLp, htarget] at hdet
    have hb := le_trans hdet hfinal
    have hω' : ε ≤ |‖specProj (s.gram N ω) (Set.Ioi (bulkEdge c + ε₁)) (s.spikeVec j N)‖ ^ 2
        - betaSq (Real.sqrt (s.coreEig j)) c| := hω
    linarith
  · exact tendsto_measure_zero_union hedgeC1
      (tendsto_measure_zero_union hedgeC2
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
              (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE3T)
                hedgeC7)))))

/-! ### 3. The margin hypothesis is not vacuous -/

/-- The hypothesis `hε₁'` of `tendstoInProb_normSq_specProj_Ioi` holds for every small
`ε₁`: every supercritical outlier `ρ_k` sits strictly above the bulk edge
(`OutliersR.bulkEdge_lt_rho`), and there are finitely many spikes. The consumer picks one
`ε₁` in `(0, ε₀]` and gets the margin at every supercritical spike at once. -/
theorem exists_margin_le_rhoSq (s : RankRStack μ ns d r) {c : ℝ} (hc : 0 < c) :
    ∃ ε₀ > 0, ∀ ε₁ : ℝ, 0 < ε₁ → ε₁ ≤ ε₀ →
      ∀ k : Fin r, c < s.coreEig k ^ 2 →
        bulkEdge c + 2 * ε₁ ≤ rhoSq (Real.sqrt (s.coreEig k)) c := by
  classical
  by_cases hne : (Finset.univ.filter fun k : Fin r => c < s.coreEig k ^ 2).Nonempty
  · refine ⟨(Finset.univ.filter fun k : Fin r => c < s.coreEig k ^ 2).inf' hne
      (fun k => (rhoSq (Real.sqrt (s.coreEig k)) c - bulkEdge c) / 2), ?_, ?_⟩
    · rw [gt_iff_lt, Finset.lt_inf'_iff]
      intro k hk
      have hk2 : c < s.coreEig k ^ 2 := (Finset.mem_filter.mp hk).2
      have hlt := OutliersR.bulkEdge_lt_rho hc (s.coreEig_nonneg k) hk2
      linarith
    · intro ε₁ hε₁ hle k hk
      have hkS : k ∈ Finset.univ.filter fun k : Fin r => c < s.coreEig k ^ 2 :=
        Finset.mem_filter.mpr ⟨Finset.mem_univ k, hk⟩
      have h1 : (Finset.univ.filter fun k : Fin r => c < s.coreEig k ^ 2).inf' hne
            (fun k => (rhoSq (Real.sqrt (s.coreEig k)) c - bulkEdge c) / 2)
          ≤ (rhoSq (Real.sqrt (s.coreEig k)) c - bulkEdge c) / 2 :=
        Finset.inf'_le _ hkS
      linarith
  · refine ⟨1, one_pos, ?_⟩
    intro ε₁ hε₁ hle k hk
    exact absurd (⟨k, Finset.mem_filter.mpr ⟨Finset.mem_univ k, hk⟩⟩ :
      (Finset.univ.filter fun k : Fin r => c < s.coreEig k ^ 2).Nonempty) hne

end RankRStack

end StackedSVD
