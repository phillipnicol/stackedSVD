/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Markov
import StackedSVD.RMT.General.Edge.Sparse
import StackedSVD.RMT.General.Sup

/-! # Stage 3, the assembly: the sharp upper edge at a general law

The endpoint of `notes/stage3_edge.md` (route C): under `SpikedModel.GeneralNoise ν` with
`NoiseLaw ν` (mean 0, variance 1, finite fourth moment) and the regime `n/d → c > 0`,
`‖Z‖²/d ≤ (1 + √c)² + ε` with probability tending to 1
(`SpikedModel.opNorm_sq_edge_of_general`, the upper half of Bai-Yin), and its corollary
`SpikedModel.lamMax_W0_edge_of_general`, the hypothesis `hedge` of
`resolventLimits_of_general` and `singleTableLaw_of_general` (`RMT/General/Sup.lean`)
verbatim; `SpikedModel.singleTableLaw_of_moments` is the single-table law with that
hypothesis discharged. The assembly composes the Markov step (K8, `‖Ẑ‖²/d ≤ bulkEdge c + ε/2`), the
discarded part (L2, `‖R‖ ≤ ε₂ √d`) and the mean shift (T9b, `‖m_T J‖ ≤ ε₂ √d`) through
the split `Z = Ẑ + R + m_T J` (T7) at `a = 1/16`. -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace Edge

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}

/-- Squeeze: a family that contains a family of measure tending to 1 has measure tending to 1. -/
theorem tendsto_measure_one_of_subset [∀ N, IsProbabilityMeasure (μ N)] {s t : ∀ N, Set (Ω N)}
    (hst : ∀ᶠ N in atTop, s N ⊆ t N) (hs : Tendsto (fun N => μ N (s N)) atTop (𝓝 1)) :
    Tendsto (fun N => μ N (t N)) atTop (𝓝 1) :=
  tendsto_of_tendsto_of_tendsto_of_le_of_le' hs tendsto_const_nhds
    (hst.mono fun _ h => measure_mono h) (Eventually.of_forall fun _ => prob_le_one)

/-- The intersection of two families of measure tending to 1 has measure tending to 1. -/
theorem tendsto_measure_one_inter [∀ N, IsProbabilityMeasure (μ N)] {s t : ∀ N, Set (Ω N)}
    (hs' : ∀ N, NullMeasurableSet (s N) (μ N)) (ht' : ∀ N, NullMeasurableSet (t N) (μ N))
    (hs : Tendsto (fun N => μ N (s N)) atTop (𝓝 1))
    (ht : Tendsto (fun N => μ N (t N)) atTop (𝓝 1)) :
    Tendsto (fun N => μ N (s N ∩ t N)) atTop (𝓝 1) := by
  refine tendsto_measure_one_of_bad (s := fun N => (s N)ᶜ ∪ (t N)ᶜ) (fun N ω hω => ?_) ?_
  · rw [Set.compl_inter] at hω; exact hω
  · exact tendsto_measure_zero_union (tendsto_measure_compl_zero hs' hs)
      (tendsto_measure_compl_zero ht' ht)

/-- `lamMax (FᵀF) ≤ ‖F‖²`, the unscaled form of `R3.lamMax_gram_le_opNorm_sq`. -/
theorem lamMax_gram_le_opNorm_sq' {p D : ℕ} (F : Matrix (Fin p) (Fin D) ℝ) :
    lamMax (Fᵀ * F) (isHermitian_transpose_mul_self F) ≤ ‖F‖ ^ 2 := by
  have h := R3.lamMax_gram_le_opNorm_sq (r := 1) zero_le_one F
    (by simpa using isHermitian_transpose_mul_self F)
  simpa using h

/-- `‖E‖² = ‖Z‖²/d`, the scaling of `Defs.lean` (`E = d^{-1/2} Z`). -/
theorem opNorm_sq_E_eq {n d : ℕ → ℕ} (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    ‖m.E N ω‖ ^ 2 = ‖m.Z N ω‖ ^ 2 / (d N : ℝ) := by
  have hd : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.hd N
  have hE : m.E N ω = (Real.sqrt (d N))⁻¹ • m.Z N ω := rfl
  rw [hE, norm_smul, mul_pow]
  simp only [norm_inv, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _)]
  rw [inv_pow, Real.sq_sqrt hd.le]
  ring

/-- A. The canonical-law endpoint: `‖Z‖²/d ≤ bulkEdge c + ε` with probability tending to 1,
under the law `noiseMatrix ν n d`. Composes K8, L2 and T9b through the split T7. -/
theorem tendsto_measure_opNorm_edge {c : ℝ} (hc : 0 < c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    {n d : ℕ → ℕ} (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N)
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha8 : a < 1 / 8) :
    ∀ ε > 0, Tendsto (fun N => (noiseMatrix ν (n N) (d N))
      {Z | ‖Z‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  intro ε hε
  have hprob : ∀ N, IsProbabilityMeasure (noiseMatrix ν (n N) (d N)) := by
    intro N; have := hν.prob; infer_instance
  -- the scalars: `s = √(bulkEdge c + ε/2)` and the slack `ε₂`
  set b : ℝ := bulkEdge c with hbdef
  have hb0 : (0 : ℝ) ≤ b := by rw [hbdef, bulkEdge]; positivity
  set s : ℝ := Real.sqrt (b + ε / 2) with hsdef
  have hs0 : (0 : ℝ) ≤ s := Real.sqrt_nonneg _
  have hssq : s ^ 2 = b + ε / 2 := Real.sq_sqrt (by linarith)
  set e2 : ℝ := min 1 (ε / (8 * (s + 1))) with he2def
  have he2pos : (0 : ℝ) < e2 := lt_min one_pos (by positivity)
  have he2one : e2 ≤ 1 := min_le_left _ _
  have he2bd : e2 * (8 * (s + 1)) ≤ ε := by
    have h := min_le_right (1 : ℝ) (ε / (8 * (s + 1)))
    have hpos : (0 : ℝ) < 8 * (s + 1) := by linarith
    calc e2 * (8 * (s + 1)) ≤ (ε / (8 * (s + 1))) * (8 * (s + 1)) := by
          exact mul_le_mul_of_nonneg_right h hpos.le
      _ = ε := by field_simp
  have hkey : (s + 2 * e2) ^ 2 ≤ b + ε := by nlinarith [hssq, he2bd, he2pos.le, he2one, hs0]
  -- the three inputs, at their instantiations
  have hA := tendsto_measure_trunc_opNorm (c := c) hc hν hn hd hdtop hcN ha (by linarith)
    (ε := ε / 4) (by linarith)
  have hB := tendsto_measure_discard_opNorm (c := c) hν hn hd hdtop hcN ha ha8 he2pos
  have hM := meanMat_opNorm_small (c := c) hν hdtop hcN ha (by linarith) he2pos
  -- the two events are measurable
  have hmA : ∀ N, NullMeasurableSet {Z : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 / (d N : ℝ) ≤ b + 2 * (ε / 4)}
      (noiseMatrix ν (n N) (d N)) := by
    intro N
    exact (measurableSet_le
      ((((R3.measurable_opNorm (n N) (d N)).comp
        (measurable_truncMat ν (truncLevel a (d N)) (n N) (d N))).pow_const 2).div_const _)
      measurable_const).nullMeasurableSet
  have hmB : ∀ N, NullMeasurableSet {Z : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      ‖discardMat (truncLevel a (d N)) Z‖ ≤ e2 * Real.sqrt (d N)}
      (noiseMatrix ν (n N) (d N)) := by
    intro N
    exact (measurableSet_le
      ((R3.measurable_opNorm (n N) (d N)).comp
        (measurable_discardMat (truncLevel a (d N)) (n N) (d N)))
      measurable_const).nullMeasurableSet
  have hAB := tendsto_measure_one_inter hmA hmB hA hB
  -- the inclusion, valid once the mean shift is small
  refine tendsto_measure_one_of_subset ?_ hAB
  filter_upwards [hM] with N hMN Z hZ
  obtain ⟨hZ1, hZ2⟩ := hZ
  have hdR : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
  have hsd : (0 : ℝ) < Real.sqrt (d N) := Real.sqrt_pos.mpr hdR
  have hsdsq : Real.sqrt (d N) ^ 2 = (d N : ℝ) := Real.sq_sqrt hdR.le
  -- the truncated part
  have h1 : ‖truncMat ν (truncLevel a (d N)) Z‖ ≤ s * Real.sqrt (d N) := by
    have hle : ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 ≤ (s * Real.sqrt (d N)) ^ 2 := by
      have := (div_le_iff₀ hdR).mp hZ1
      calc ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 ≤ (b + 2 * (ε / 4)) * (d N : ℝ) := this
        _ = (s * Real.sqrt (d N)) ^ 2 := by rw [mul_pow, hssq, hsdsq]; ring
    nlinarith [norm_nonneg (truncMat ν (truncLevel a (d N)) Z), mul_nonneg hs0 hsd.le, hle]
  -- the three-term split
  set A : Matrix (Fin (n N)) (Fin (d N)) ℝ := truncMat ν (truncLevel a (d N)) Z with hAdef
  set B : Matrix (Fin (n N)) (Fin (d N)) ℝ := discardMat (truncLevel a (d N)) Z with hBdef
  set M : Matrix (Fin (n N)) (Fin (d N)) ℝ := meanMat ν (truncLevel a (d N)) (n N) (d N)
    with hMdef
  have hsplit : ‖Z‖ = ‖A + B + M‖ := by
    rw [hAdef, hBdef, hMdef, ← mat_split ν (truncLevel a (d N)) Z]
  have h2 : ‖Z‖ ≤ (s + 2 * e2) * Real.sqrt (d N) := by
    calc ‖Z‖ = ‖A + B + M‖ := hsplit
      _ ≤ ‖A + B‖ + ‖M‖ := norm_add_le _ _
      _ ≤ (‖A‖ + ‖B‖) + ‖M‖ := add_le_add (norm_add_le _ _) le_rfl
      _ ≤ (s * Real.sqrt (d N) + e2 * Real.sqrt (d N)) + e2 * Real.sqrt (d N) :=
            add_le_add (add_le_add h1 hZ2) hMN
      _ = (s + 2 * e2) * Real.sqrt (d N) := by ring
  -- square and divide
  have h3 : ‖Z‖ ^ 2 ≤ (b + ε) * (d N : ℝ) := by
    have hpos : (0 : ℝ) ≤ (s + 2 * e2) * Real.sqrt (d N) := by positivity
    have := pow_le_pow_left₀ (norm_nonneg Z) h2 2
    calc ‖Z‖ ^ 2 ≤ ((s + 2 * e2) * Real.sqrt (d N)) ^ 2 := this
      _ = (s + 2 * e2) ^ 2 * (d N : ℝ) := by rw [mul_pow, hsdsq]
      _ ≤ (b + ε) * (d N : ℝ) := by
          exact mul_le_mul_of_nonneg_right hkey hdR.le
  exact (div_le_iff₀ hdR).mpr h3

end Edge

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **Stage 3 target.** The sharp upper edge of the unscaled noise Gram matrix at a general
law: `‖Z‖² / d ≤ (1 + √c)² + ε` with probability tending to 1. Upper half of Bai-Yin. -/
theorem opNorm_sq_edge_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | ‖m.Z N ω‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  intro ε hε
  have h := Edge.tendsto_measure_opNorm_edge hc hν m.hn m.hd hreg.2.1 hreg.2.2
    (a := 1 / 16) (by norm_num) (by norm_num) ε hε
  have hmeas : ∀ N : ℕ, MeasurableSet {Z : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      ‖Z‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε} := fun N =>
    measurableSet_le (((R3.measurable_opNorm (n N) (d N)).pow_const 2).div_const _)
      measurable_const
  have heq : ∀ N : ℕ, μ N {ω | ‖m.Z N ω‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε}
      = (noiseMatrix ν (n N) (d N)) {Z | ‖Z‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε} :=
    fun N => (hG N).measure_eq (hmeas N)
  simpa only [heq] using h

/-- **Stage 3 corollary.** `hedge` of `RMT/General/Sup.lean` verbatim, by the contraction
`lamMax W₀ ≤ ‖E⊥‖² ≤ ‖E‖² = ‖Z‖² / d`. -/
theorem lamMax_W0_edge_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  intro ε hε
  refine Edge.tendsto_measure_one_of_subset
    (s := fun N => {ω | ‖m.Z N ω‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε})
    (Eventually.of_forall fun N ω hω => ?_)
    (m.opNorm_sq_edge_of_general hc hreg hν hG ε hε)
  have h1 : lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ ‖m.Eperp N ω‖ ^ 2 :=
    Edge.lamMax_gram_le_opNorm_sq' (m.Eperp N ω)
  have h2 : ‖m.Eperp N ω‖ ^ 2 ≤ ‖m.E N ω‖ ^ 2 :=
    pow_le_pow_left₀ (norm_nonneg _) (Edge.opNorm_Eperp_le m N ω) 2
  have h3 : ‖m.E N ω‖ ^ 2 = ‖m.Z N ω‖ ^ 2 / (d N : ℝ) := Edge.opNorm_sq_E_eq m N ω
  have h4 : ‖m.Z N ω‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε := hω
  change lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε
  linarith [h1, h2, h3.le, h3.ge, h4]

/-- **(d). `prop:single_table` at `assum:general_noise`, the edge discharged.** The Stage 1
endpoint `singleTableLaw_of_general` (`RMT/General/Sup.lean`) with its hypothesis `hedge`
supplied by `lamMax_W0_edge_of_general`. Against the paper's assumption (i.i.d. entries,
mean `0`, variance `1`, finite fourth moment: `hν`) one binder remains that the paper does
not have, the Lebesgue density `hac`, which the simplicity of the top eigenvalue
(`singleTableLaw_topSimple_of_general`) uses. Both regimes, no `hn2 : ∀ N, 2 ≤ n N`. -/
theorem singleTableLaw_of_moments [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν) :
    m.SingleTableLaw c :=
  m.singleTableLaw_of_general hc hreg hν hac hG
    (m.lamMax_W0_edge_of_general hc hreg hν hG)

end SpikedModel
end StackedSVD
