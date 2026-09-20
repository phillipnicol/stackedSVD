/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.EdgeTauR
import StackedSVD.RankR.RMT.AlignOutG
import StackedSVD.LinAlg.SpecWindow

/-!
# Task C3: the outlier overlap and the count at a general threshold

Task C3 of `notes/archive/rankr_plan_C.md` section 2. The rank-1 mirror is
`RankRStack.align_of_gaussian_supercritical` (`RankR/RMT/Outliers.lean`), and the fixed
threshold twin of this file is `RankRStack.tendstoInProb_normSq_specProj_Ioi`
(`RankR/RMT/AlignOutG.lean`, task A2a).

`AlignOutG.lean` runs the frame at the threshold `bulkEdge c + ε₁` and splits the spikes on
`coreEig k ^ 2 ≤ c`. Task C5 needs the same two statements at a threshold `τ` that the caller
picks, so that two thresholds trap one sorted eigenvalue. This file takes the split on the
outlier position instead: a spike is in the frame when `τ < ρ_k`, and `hsep` asks every spike
to stay off the threshold by the margin `mg`.

Contents:

1. `OutliersR.supercritical_of_lt_rhoSq`: a spike whose outlier sits above `τ > bulkEdge c` is
   supercritical. Below the threshold `rhoSq θ c = bulkEdge c` (`Defs.lean`), so the outlier
   never passes `τ` there.
2. `OutliersR.exists_sum_equiv_split_pred` and its corollary
   `OutliersR.exists_sum_equiv_split_rho`: the split of `Fin r` on `τ < ρ_k`, in the
   `Fin t ⊕ Fin u ≃ Fin r` shape that `RankRStack.tendsto_measure_eigenvalues₀_le_tau` takes,
   with the size `u` of the right summand read as a filter cardinality.
3. `RankRStack.tendstoInProb_normSq_specProj_Ioi_tau`: the overlap of the spike `j` with the
   eigenvalues above `τ`. The limit is `betaSq` when `τ < ρ_j` and `0` otherwise, because the
   spike carries a frame vector exactly in the first case.
4. `RankRStack.tendsto_measure_count_Ioi_tau`: the count. With probability tending to 1 the
   sorted eigenvalues above `τ` are the indices below the number of outliers above `τ`. The
   `↔` shape comes from `Frame.lt_eigenvalues₀_iff_of_card_eigenvalues`
   (`LinAlg/SpecWindow.lean`, task C4), which turns the count of
   `OutliersR.count_eq_of_forms_of_edge` into a statement about every sorted index.

Both proofs follow `AlignOutG.lean` line for line, with `τ` for `bulkEdge c + ε₁`, `mg` for
the margin, and `RankRStack.tendsto_measure_eigenvalues₀_le_tau`
(`RankR/RMT/EdgeTauR.lean`, task C2) for the U7a edge event. The column-norm constant `Cq` is
read at the fixed point `bulkEdge c + 1` and does not depend on `τ`.

The count is a with-probability-tending-to-one statement, not an almost sure one. Numeric
check `$SP/agents/planC/check_sandwich.py`, seed 20260902, row (vii): 13 of 20 draws at
`d = 200`, 20 of 20 at `d = 800`.

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace OutliersR

/-! ### 1. Two scalar preliminaries -/

/-- A spike whose outlier sits above a threshold `τ` above the bulk edge is supercritical.
`rhoSq θ c` is `bulkEdge c` below the threshold `c < θ ^ 4` (`Defs.lean`), so a subcritical
spike has `ρ = bulkEdge c < τ`. This is how the frame indices of task C3 inherit the
supercriticality that `OutliersR.bulkEdge_lt_rho`, `mul_m_rho_eq_neg_one` and `nu_pos`
need. -/
theorem supercritical_of_lt_rhoSq {c lam τ : ℝ} (hlam0 : 0 ≤ lam) (hτ : bulkEdge c < τ)
    (hlt : τ < rhoSq (Real.sqrt lam) c) : c < lam ^ 2 := by
  by_contra hcon
  have hne : ¬ (Real.sqrt lam ^ 4 > c) := by
    rw [gt_iff_lt, sqrt_pow_four hlam0]
    exact hcon
  unfold rhoSq at hlt
  rw [if_neg hne] at hlt
  linarith

/-- The split of `Fin r` on a decidable predicate `P`, in the `Fin t ⊕ Fin u ≃ Fin r` shape
that `RankRStack.tendsto_measure_eigenvalues₀_le_tau` takes. The left summand carries the
indices where `P` fails, the right summand the indices where `P` holds, and `u` is the filter
cardinality of `P`. `exists_sum_equiv_split_rho` is the case `P k := τ < g k`; the older
`OutliersR.exists_sum_equiv_split` (`RankR/RMT/AlignOutG.lean`) is the case
`P k := c < lam k` without the cardinality conclusion. -/
theorem exists_sum_equiv_split_pred {r : ℕ} (P : Fin r → Prop) [DecidablePred P] :
    ∃ (t u : ℕ) (e : Fin t ⊕ Fin u ≃ Fin r),
      (∀ a : Fin t, ¬ P (e (Sum.inl a))) ∧ (∀ b : Fin u, P (e (Sum.inr b))) ∧
        u = (Finset.univ.filter P).card := by
  refine ⟨Fintype.card {k : Fin r // ¬ P k}, Fintype.card {k : Fin r // P k},
    (Equiv.sumCongr (Fintype.equivFin {k : Fin r // ¬ P k}).symm
        (Fintype.equivFin {k : Fin r // P k}).symm).trans
      ((Equiv.sumComm _ _).trans (Equiv.sumCompl P)), ?_, ?_, ?_⟩
  · intro a
    exact ((Fintype.equivFin {k : Fin r // ¬ P k}).symm a).2
  · intro b
    exact ((Fintype.equivFin {k : Fin r // P k}).symm b).2
  · exact Fintype.card_subtype _

/-- The split of the spikes on the position of the outlier. The left summand collects the
spikes with `g k ≤ τ`, the right one those with `τ < g k`, and `u` is the number of the
second kind. This is `exists_sum_equiv_split_pred` at `P k := τ < g k`. -/
theorem exists_sum_equiv_split_rho {r : ℕ} (g : Fin r → ℝ) (τ : ℝ) :
    ∃ (t u : ℕ) (e : Fin t ⊕ Fin u ≃ Fin r),
      (∀ a : Fin t, g (e (Sum.inl a)) ≤ τ) ∧
        (∀ b : Fin u, τ < g (e (Sum.inr b))) ∧
        u = (Finset.univ.filter fun l : Fin r => τ < g l).card := by
  obtain ⟨t, u, e, hinl, hinr, hu⟩ := exists_sum_equiv_split_pred fun k : Fin r => τ < g k
  exact ⟨t, u, e, fun a => not_lt.mp (hinl a), hinr, hu⟩

end OutliersR

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 2. The outlier part of the spike overlap at a general threshold -/

/-- **Task C3.1, the outlier overlap at a general threshold.** For a `RankRStack` with
Gaussian noise and aspect ratio `c`, the part of the projection of the spike direction `V q_j`
on the eigenvalues of the Gram matrix above `τ` converges in probability to
`betaSq (√λ_j) c` when the outlier `ρ_j` sits above `τ`, and to `0` when it does not.

This is `RankRStack.tendstoInProb_normSq_specProj_Ioi` (`RankR/RMT/AlignOutG.lean`, task A2a)
with the free threshold `τ` in place of `bulkEdge c + ε₁`. The hypothesis `hsep` asks every
outlier to stay off `τ` by the margin `mg`; ties among the outliers cost nothing, because two
equal outliers fall on the same side. The frame runs on the spikes with `τ < ρ_k` only, so the
`if` in the limit is exactly the membership of `j` in the frame.

Task C5 uses this at two thresholds that trap one sorted eigenvalue, and reads the difference
with `normSq_specProj_Ioc` (`LinAlg/SpecWindow.lean`). -/
theorem tendstoInProb_normSq_specProj_Ioi_tau [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) (j : Fin r)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : bulkEdge c < τ)
    (hsep : ∀ k : Fin r, rhoSq (Real.sqrt (s.coreEig k)) c < τ
      ∨ τ + mg ≤ rhoSq (Real.sqrt (s.coreEig k)) c) :
    TendstoInProb μ
      (fun N ω => ‖specProj (s.gram N ω) (Set.Ioi τ) (s.spikeVec j N)‖ ^ 2)
      (if τ < rhoSq (Real.sqrt (s.coreEig j)) c then betaSq (Real.sqrt (s.coreEig j)) c
        else 0) := by
  classical
  obtain ⟨U, hU, hsig, hgram, h⟩ := s.resolventLimitsR_of_gaussian hG hc hdtop hns hn hp
  obtain ⟨t, u, e, hsubE, hsupE, hcardu⟩ :=
    OutliersR.exists_sum_equiv_split_rho (fun k => rhoSq (Real.sqrt (s.coreEig k)) c) τ
  set f : Fin u → Fin r := fun b => e (Sum.inr b) with hfdef
  have hfinj : Function.Injective f := fun a b hab => Sum.inr_injective (e.injective hab)
  have hlam0 : ∀ k, 0 ≤ s.coreEig k := s.coreEig_nonneg
  have hsupE' : ∀ k : Fin u, τ < rhoSq (Real.sqrt (s.coreEig (f k))) c := fun k => hsupE k
  have hsubE' : ∀ a : Fin t, rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c < τ := by
    intro a
    have h1 : rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c ≤ τ := hsubE a
    rcases hsep (e (Sum.inl a)) with hlt | hge
    · exact hlt
    · linarith
  have hsupf : ∀ k : Fin u, c < s.coreEig (f k) ^ 2 := fun k =>
    OutliersR.supercritical_of_lt_rhoSq (hlam0 (f k)) hτ (hsupE' k)
  -- 1. the frame scalars, on the spikes whose outlier passes `τ`
  set ρv : Fin u → ℝ := fun k => rhoSq (Real.sqrt (s.coreEig (f k))) c with hρv
  set νv : Fin u → ℝ := fun k => (s.coreEig (f k) + 1) * MP.mDeriv c (ρv k) with hνv
  have hbρ : ∀ k, bulkEdge c < ρv k := fun k =>
    OutliersR.bulkEdge_lt_rho hc (hlam0 (f k)) (hsupf k)
  have hm1 : ∀ k, (s.coreEig (f k) + 1) * MP.m c (ρv k) = -1 := fun k =>
    OutliersR.mul_m_rho_eq_neg_one hc (hlam0 (f k)) (hsupf k)
  have hνpos : ∀ k, 0 < νv k := fun k => OutliersR.nu_pos hc (hlam0 (f k)) (hsupf k)
  have hτρ : ∀ k, τ + mg ≤ ρv k := by
    intro k
    have h1 : τ < rhoSq (Real.sqrt (s.coreEig (f k))) c := hsupE' k
    have h2 : ρv k = rhoSq (Real.sqrt (s.coreEig (f k))) c := rfl
    rcases hsep (f k) with hlt | hge
    · linarith
    · rw [h2]
      exact hge
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
  -- 3. the overlap target, one function for both branches of the `if`
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
  -- the scalar identity: the frame sum is `betaSq` exactly when `j` carries a frame vector
  have htarget : ∑ k, tv k ^ 2 / νv k
      = if τ < rhoSq (Real.sqrt (s.coreEig j)) c then betaSq (Real.sqrt (s.coreEig j)) c
        else 0 := by
    by_cases hj : τ < rhoSq (Real.sqrt (s.coreEig j)) c
    · rw [if_pos hj]
      obtain ⟨k₀, hk₀⟩ : ∃ k₀ : Fin u, f k₀ = j := by
        rcases hsym : e.symm j with a | b
        · exfalso
          have hj1 : e (Sum.inl a) = j := by rw [← hsym]; exact e.apply_symm_apply j
          have h1 : rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c < τ := hsubE' a
          rw [hj1] at h1
          linarith
        · refine ⟨b, ?_⟩
          have hj1 : e (Sum.inr b) = j := by rw [← hsym]; exact e.apply_symm_apply j
          exact hj1
      have hsupj : c < s.coreEig j ^ 2 :=
        OutliersR.supercritical_of_lt_rhoSq (hlam0 j) hτ hj
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
    · rw [if_neg hj]
      have hno : ∀ k : Fin u, f k ≠ j := by
        intro k hcon
        have h1 : τ < rhoSq (Real.sqrt (s.coreEig (f k))) c := hsupE' k
        rw [hcon] at h1
        exact hj h1
      have hzero : ∀ k : Fin u, tv k ^ 2 / νv k = 0 := by
        intro k
        have h1 : tv k = if f k = j then Lval else 0 := rfl
        rw [h1, if_neg (hno k)]
        simp
      rw [Finset.sum_congr rfl fun k _ => hzero k]
      simp
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
  have hX0 : (0 : ℝ) ≤ (u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg) :=
    mul_nonneg (Nat.cast_nonneg u) (add_nonneg hCG0 (div_nonneg hCR0 hmg.le))
  have hB0 : (0 : ℝ) ≤ 5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ := by
    have h1 := mul_nonneg hL0 hνsum0
    linarith
  set η : ℝ := min (min 1 (mg / (OutliersR.resCG Cq r νv + 1)))
      (min (1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
          + OutliersR.resCG Cq r νv / mg) + 1)))
        (ε / (2 * (5 * ((u : ℝ) * (OutliersR.gramC ρv νv
            + OutliersR.resCG Cq r νv / mg))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    refine lt_min (lt_min one_pos (div_pos hmg (by linarith))) (lt_min ?_ ?_)
    · exact div_pos one_pos (by linarith)
    · exact div_pos hε (by linarith)
  have hη1 : η ≤ 1 := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_left _ _)
  have hηm : η ≤ mg / (OutliersR.resCG Cq r νv + 1) := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_right _ _)
  have hηX : η ≤ 1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg) + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hηε : η ≤ ε / (2 * (5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_right _ _)
  have hδm : OutliersR.resCG Cq r νv * η ≤ mg :=
    OutliersR.mul_le_of_le_div_add_one hCR0 hηm hmg.le
  have hsmall : (u : ℝ) * (OutliersR.gramC ρv νv * η
      + OutliersR.resCG Cq r νv * η / mg) ≤ 1 / 2 := by
    have hXe : (u : ℝ) * (OutliersR.gramC ρv νv * η + OutliersR.resCG Cq r νv * η / mg)
        = ((u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg)) * η := by
      ring
    rw [hXe]
    exact le_trans (mul_le_mul_of_nonneg_left hηX hX0) (OutliersR.mul_inv_le_half hX0)
  have hfinal : 5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
        + OutliersR.resCG Cq r νv * η / mg)
      + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η ≤ ε / 2 := by
    have hKe : 5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
          + OutliersR.resCG Cq r νv * η / mg)
        + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η
        = (5 * ((u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹) * η := by
      ring
    rw [hKe]
    exact le_trans (mul_le_mul_of_nonneg_left hηε hB0)
      (OutliersR.mul_div_le_half hB0 hε.le)
  -- 6. the seven bad families
  have hε2 : (0 : ℝ) < (τ - bulkEdge c) / 2 := by linarith
  have hedgeC1 : Tendsto (fun N => μ N {ω | lamMax (s.rankRW0 N ω (U N))
      (s.isHermitian_rankRW0 N ω (U N)) ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + (τ - bulkEdge c) / 2)).nullMeasurableSet)
      (h.edge ((τ - bulkEdge c) / 2) hε2)
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
        (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_eigenvalues₀_le N u τ).nullMeasurableSet)
      (s.tendsto_measure_eigenvalues₀_le_tau hc h hgram e hτ hsubE')
  -- 7. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
          ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ
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
                      (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ}ᶜ)))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hg1, hg2, hg3, hg4, hg5, hg6, hg7⟩ := hbad
    have hlt2 : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        < bulkEdge c + 1 := by linarith
    have hlamτ : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N)) ≤ τ := by
      linarith
    have hcolb : ∀ l : Fin r, ∑ i, s.qmatR N ω (U N) i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      exact OutliersR.dot_self_le_of_cform_close (s.isHermitian_rankRW0 N ω (U N))
        (fun a => s.eigenvalues_rankRW0_nonneg N ω (U N) a) hz₁0 hlt2
        (fun i => s.qmatR N ω (U N) i l) (hg3 l).le
    have hI := OutliersR.count_eq_of_forms_of_edge (Q := s.qmatR N ω (U N))
      (τ := τ) (mg := mg) (Cq := Cq) (η := η)
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω) (hgram N ω) hfinj hmg
      hlamτ hτρ hνpos hCq0 hη0 hcolb (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le)
      hsmall hδm hg7
    have hdet := OutliersR.align_detG (Q := s.qmatR N ω (U N))
      (τ := τ) (mg := mg) (Cq := Cq) (t := tv) (Lb := |Lval|) (η := η)
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω) (hgram N ω) hfinj hmg
      hlamτ hτρ hνpos hI hCq0 hcolb (hvdot N) hLb hη0 hη1
      (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le) (fun k => (hg6 k).le) hsmall
    rw [WithLp.toLp_ofLp, htarget] at hdet
    have hb := le_trans hdet hfinal
    have hω' : ε ≤ |‖specProj (s.gram N ω) (Set.Ioi τ) (s.spikeVec j N)‖ ^ 2
        - (if τ < rhoSq (Real.sqrt (s.coreEig j)) c then betaSq (Real.sqrt (s.coreEig j)) c
            else 0)| := hω
    linarith
  · exact tendsto_measure_zero_union hedgeC1
      (tendsto_measure_zero_union hedgeC2
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
              (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE3T)
                hedgeC7)))))

/-! ### 3. The count of the eigenvalues above the threshold -/

/-- **Task C3.2, the count at a general threshold.** With probability tending to 1 the sorted
eigenvalues of the Gram matrix above `τ` are exactly the indices below the number of outliers
above `τ`.

The `↔` shape is what the window lemmas of `LinAlg/SpecWindow.lean` (task C4) take: two
thresholds with counts `j + 1` and `j` trap the sorted eigenvalue `λ_j` in a half open window.
The upper half of the count is the edge theorem
`RankRStack.tendsto_measure_eigenvalues₀_le_tau` (`RankR/RMT/EdgeTauR.lean`, task C2), the
lower half is the frame of `OutliersR.count_eq_of_forms_of_edge`, and
`Frame.lt_eigenvalues₀_iff_of_card_eigenvalues` moves the count to every sorted index. -/
theorem tendsto_measure_count_Ioi_tau [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : bulkEdge c < τ)
    (hsep : ∀ k : Fin r, rhoSq (Real.sqrt (s.coreEig k)) c < τ
      ∨ τ + mg ≤ rhoSq (Real.sqrt (s.coreEig k)) c) :
    Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (d N))),
      (τ < (s.isHermitian_gram N ω).eigenvalues₀ k
        ↔ (k : ℕ) < (Finset.univ.filter fun l : Fin r =>
            τ < rhoSq (Real.sqrt (s.coreEig l)) c).card)}) atTop (𝓝 1) := by
  classical
  obtain ⟨U, hU, hsig, hgram, h⟩ := s.resolventLimitsR_of_gaussian hG hc hdtop hns hn hp
  obtain ⟨t, u, e, hsubE, hsupE, hcardu⟩ :=
    OutliersR.exists_sum_equiv_split_rho (fun k => rhoSq (Real.sqrt (s.coreEig k)) c) τ
  set f : Fin u → Fin r := fun b => e (Sum.inr b) with hfdef
  have hfinj : Function.Injective f := fun a b hab => Sum.inr_injective (e.injective hab)
  have hlam0 : ∀ k, 0 ≤ s.coreEig k := s.coreEig_nonneg
  have hsupE' : ∀ k : Fin u, τ < rhoSq (Real.sqrt (s.coreEig (f k))) c := fun k => hsupE k
  have hsubE' : ∀ a : Fin t, rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c < τ := by
    intro a
    have h1 : rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c ≤ τ := hsubE a
    rcases hsep (e (Sum.inl a)) with hlt | hge
    · exact hlt
    · linarith
  have hsupf : ∀ k : Fin u, c < s.coreEig (f k) ^ 2 := fun k =>
    OutliersR.supercritical_of_lt_rhoSq (hlam0 (f k)) hτ (hsupE' k)
  -- 1. the frame scalars
  set ρv : Fin u → ℝ := fun k => rhoSq (Real.sqrt (s.coreEig (f k))) c with hρv
  set νv : Fin u → ℝ := fun k => (s.coreEig (f k) + 1) * MP.mDeriv c (ρv k) with hνv
  have hbρ : ∀ k, bulkEdge c < ρv k := fun k =>
    OutliersR.bulkEdge_lt_rho hc (hlam0 (f k)) (hsupf k)
  have hm1 : ∀ k, (s.coreEig (f k) + 1) * MP.m c (ρv k) = -1 := fun k =>
    OutliersR.mul_m_rho_eq_neg_one hc (hlam0 (f k)) (hsupf k)
  have hνpos : ∀ k, 0 < νv k := fun k => OutliersR.nu_pos hc (hlam0 (f k)) (hsupf k)
  have hτρ : ∀ k, τ + mg ≤ ρv k := by
    intro k
    have h1 : τ < rhoSq (Real.sqrt (s.coreEig (f k))) c := hsupE' k
    have h2 : ρv k = rhoSq (Real.sqrt (s.coreEig (f k))) c := rfl
    rcases hsep (f k) with hlt | hge
    · linarith
    · rw [h2]
      exact hge
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
  -- 3. the accuracy `η`: the count needs the frame smallness only
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv := OutliersR.resCG_nonneg Cq r νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  have hX0 : (0 : ℝ) ≤ (u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg) :=
    mul_nonneg (Nat.cast_nonneg u) (add_nonneg hCG0 (div_nonneg hCR0 hmg.le))
  set η : ℝ := min (mg / (OutliersR.resCG Cq r νv + 1))
      (1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
        + OutliersR.resCG Cq r νv / mg) + 1))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    exact lt_min (div_pos hmg (by linarith)) (div_pos one_pos (by linarith))
  have hηm : η ≤ mg / (OutliersR.resCG Cq r νv + 1) := by
    rw [hηdef]
    exact min_le_left _ _
  have hηX : η ≤ 1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg) + 1)) := by
    rw [hηdef]
    exact min_le_right _ _
  have hδm : OutliersR.resCG Cq r νv * η ≤ mg :=
    OutliersR.mul_le_of_le_div_add_one hCR0 hηm hmg.le
  have hsmall : (u : ℝ) * (OutliersR.gramC ρv νv * η
      + OutliersR.resCG Cq r νv * η / mg) ≤ 1 / 2 := by
    have hXe : (u : ℝ) * (OutliersR.gramC ρv νv * η + OutliersR.resCG Cq r νv * η / mg)
        = ((u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg)) * η := by
      ring
    rw [hXe]
    exact le_trans (mul_le_mul_of_nonneg_left hηX hX0) (OutliersR.mul_inv_le_half hX0)
  -- 4. the six bad families
  have hε2 : (0 : ℝ) < (τ - bulkEdge c) / 2 := by linarith
  have hedgeC1 : Tendsto (fun N => μ N {ω | lamMax (s.rankRW0 N ω (U N))
      (s.isHermitian_rankRW0 N ω (U N)) ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + (τ - bulkEdge c) / 2)).nullMeasurableSet)
      (h.edge ((τ - bulkEdge c) / 2) hε2)
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
  have hedgeC7 : Tendsto (fun N => μ N
      {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
        (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_eigenvalues₀_le N u τ).nullMeasurableSet)
      (s.tendsto_measure_eigenvalues₀_le_tau hc h hgram e hτ hsubE')
  -- 5. assemble
  refine tendsto_measure_one_of_bad (s := fun N =>
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
          ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ
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
                ∪ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
                    (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ}ᶜ))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hg1, hg2, hg3, hg4, hg5, hg7⟩ := hbad
    have hlt2 : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        < bulkEdge c + 1 := by linarith
    have hlamτ : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N)) ≤ τ := by
      linarith
    have hcolb : ∀ l : Fin r, ∑ i, s.qmatR N ω (U N) i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      exact OutliersR.dot_self_le_of_cform_close (s.isHermitian_rankRW0 N ω (U N))
        (fun a => s.eigenvalues_rankRW0_nonneg N ω (U N) a) hz₁0 hlt2
        (fun i => s.qmatR N ω (U N) i l) (hg3 l).le
    have hI := OutliersR.count_eq_of_forms_of_edge (Q := s.qmatR N ω (U N))
      (τ := τ) (mg := mg) (Cq := Cq) (η := η)
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω) (hgram N ω) hfinj hmg
      hlamτ hτρ hνpos hCq0 hη0 hcolb (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le)
      hsmall hδm hg7
    apply hω
    intro k
    rw [← hcardu]
    exact Frame.lt_eigenvalues₀_iff_of_card_eigenvalues (s.isHermitian_gram N ω) hI k
  · exact tendsto_measure_zero_union hedgeC1
      (tendsto_measure_zero_union hedgeC2
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
              hedgeC7))))

end RankRStack

end StackedSVD
