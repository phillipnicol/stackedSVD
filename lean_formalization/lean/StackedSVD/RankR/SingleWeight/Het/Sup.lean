/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.Align
import StackedSVD.RankR.SingleWeight.Het.Sub
import StackedSVD.RankR.SingleWeight.Main
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.Het.Edge

/-!
# Track G, unit G5: `SingleWeightLaw` for Gaussian noise, and the paper's facades

Unit G5 of `notes/archive/trackG_plan.md` (2026-09-05), the last unit of Track G. Mirrors
`RankR/Het/Sup.lean` (`heteroLawR_of_gaussian_aux`, `heteroLawR_of_gaussian`, and the three
`thm_rank_r_stacksvd_*_gaussian` facades) and the tail-shift lemmas of `RankR/Het/Edge.lean`
lines 300 to 400, at general `rk` (the single-weight model carries no `alignedRk`
restriction) and under `SingleWeight.EigSep` in place of a three-field hypothesis structure.

## The weights

The paper puts no condition on the weights `w_i`. Neither do the theorems below (item F18c,
2026-09-05; the first version took `hw : ∀ i, w i ≠ 0`, decision D37). `EigSep` itself
supplies a table of nonzero weight once `r ≥ 1` (`SingleWeight.EigSep.exists_ne_zero`), and
at `r = 0` every field of `SingleWeightLaw` is vacuous. The tables of weight zero drop out of
the Gram matrix; `simpleSpec_ae_stackGramW_of_exists` (`RankR/SingleWeight/Het/Sub.lean`) is
the simplicity statement on the tables that remain.

## The index range at general `rk`

`SingleWeightLaw.simpleIdx` asks for the almost sure simplicity of the sorted eigenvalue `l`
of the `d N × d N` weighted stack Gram matrix. That matrix has rank at most the total row
count of the tables of nonzero weight, `∑ i with w i ≠ 0, n i N`, so the field is false once
that count is `≤ l < d N - 1`. The aligned family avoids this because `r ≤ n_0 N` there; at
general `rk` nothing bounds `r` by `n`, so:

* the discharge for `m` itself (`singleWeightLaw_of_gaussian_aux`) takes the explicit
  hypothesis `hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N` (a coordinator decision, D38);
* the two facades take no such hypothesis: they are limits, so they run on a shift
  `m.shift k` where `r ≤ ∑ i with w i ≠ 0, n i (N + k)` holds for every `N`
  (`exists_shift_sw`), and the law is carried back by `SpikedModel.tendstoInProb_of_shift`;
* `simpleIdx_sw_of_gaussian` never needs `r ≤ d N`: it reads
  `simpleSpec_ae_stackGramW_of_exists` at `k := min r (d N)` and treats `d N ≤ l` as a
  vacuous case of `SimpleIdx`.

## Route

`SingleWeightLaw.of_shift` carries the `align` field back through
`SpikedModel.tendstoInProb_of_shift` (every object of `align` is `rfl`-equal to its shifted
twin) and takes `simpleIdx` directly, since it is a statement about one `N` at a time.
`exists_shift_sw` is the shift that makes `d N = p N + r` and `r ≤ ∑ i with w i ≠ 0, n i N`
hold together for every `N`, from `Tendsto d atTop atTop` and `Tendsto (n i₀) atTop atTop`
of one table `i₀` of nonzero weight. `singleWeightLaw_of_gaussian_aux` assembles the law
under the side condition `d N = p N + r` from `align_sw_of_gaussian` (G4) and
`simpleIdx_sw_of_gaussian`; `singleWeightLaw_of_gaussian` and
`singleWeightLaw_of_gaussian_shift` remove the side condition through the shift; the two
facades apply `prop_gen_rank_stacksvd_singleweight` and `_inner`
(`RankR/SingleWeight/Main.lean`) to the shifted law and transport the conclusion back.

Every declaration is proved in full.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. `SingleWeightLaw` transports through the tail shift -/

/-- **`SingleWeightLaw` is a tail property, except for `simpleIdx`.** The `align` field comes
from the shifted model through `SpikedModel.tendstoInProb_of_shift`; every object in it is
`rfl`-equal to its shifted twin (`shift_stackXW`, `shift_colVecG`). `simpleIdx` is a statement
about one `N` at a time, so it enters as the argument `hsimple`. Mirror: `HeteroLawR.of_shift`
(`RankR/Het/Edge.lean:365`). -/
theorem SingleWeightLaw.of_shift {m : UnalignedModelR μ M n d r rk} {w c : Fin M → ℝ}
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)} (k : ℕ)
    (law : (m.shift k).SingleWeightLaw w c γ z)
    (hsimple : ∀ (l : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N),
      SimpleIdx (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)) :
    m.SingleWeightLaw w c γ z :=
  ⟨fun l j => SpikedModel.tendstoInProb_of_shift k (law.align l j), hsimple⟩

/-- At `r = 0` there is no component, so `SingleWeightLaw` holds for every model. -/
theorem SingleWeightLaw.of_r_eq_zero (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) (hr : r = 0) :
    m.SingleWeightLaw w c γ z :=
  ⟨fun l _ => absurd l.isLt (by omega), fun l _ => absurd l.isLt (by omega)⟩

/-! ### 2. Simplicity at general `rk` -/

/-- **The `simpleIdx` field for Gaussian noise, at general `rk`.**
`simpleSpec_ae_stackGramW_of_exists` gives simplicity of every one of the top `min r (d N)`
eigenvalues; at `(l : ℕ) < d N` this covers index `l` (`simpleIdx_of_simpleSpec`); at
`d N ≤ (l : ℕ)` the goal is vacuous, since no `Fin (Fintype.card (Fin (d N)))` can equal
`(l : ℕ)`. Mirror: `simpleIdxJ_of_gaussian` (`RankR/Het/Simplicity.lean:206`), specialized
to one index with no `ellR` reindexing. -/
theorem simpleIdx_sw_of_gaussian (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (hG : m.JointGaussianNoise) (hex : ∃ i, w i ≠ 0) {N : ℕ}
    (hrn : r ≤ ∑ i with w i ≠ 0, n i N) (l : Fin r) :
    ∀ᵐ ω ∂(μ N), SimpleIdx (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ) := by
  classical
  have hae := m.simpleSpec_ae_stackGramW_of_exists w hG hex (min r (d N)) N
    (le_min ((min_le_left _ _).trans hrn) (min_le_right _ _))
  filter_upwards [hae] with ω hω
  by_cases hl : (l : ℕ) < d N
  · exact simpleIdx_of_simpleSpec hω (lt_min l.isLt hl)
  · intro q q' hq _
    exfalso
    have hq' := q.isLt
    have hcard : Fintype.card (Fin (d N)) = d N := Fintype.card_fin (d N)
    omega

/-! ### 3. The shift that supplies both side conditions -/

/-- `Tendsto d atTop atTop` and `Tendsto (n i₀) atTop atTop` of one table `i₀` of nonzero
weight give a shift after which `d (N + k) = p N + r` and `r ≤ ∑ i with w i ≠ 0, n i (N + k)`
hold together, for every `N`. Mirror: `exists_shift_add` (`RankR/Het/Edge.lean:313`), whose
second conjunct `0 < p N` is here strengthened to the row bound `simpleIdx_sw_of_gaussian`
needs at general `rk`. -/
theorem exists_shift_sw (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ}
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {w : Fin M → ℝ} {i₀ : Fin M} (hi₀ : w i₀ ≠ 0) :
    ∃ k, ∃ p : ℕ → ℕ, (∀ N, d (N + k) = p N + r) ∧
      ∀ N, r ≤ ∑ i with w i ≠ 0, n i (N + k) := by
  obtain ⟨k, hk⟩ := eventually_atTop.mp
    (((hreg i₀).2.1.eventually_ge_atTop (r + 1)).and ((hreg i₀).1.eventually_ge_atTop r))
  refine ⟨k, fun N => d (N + k) - r, fun N => ?_, fun N => ?_⟩
  · have h := (hk (N + k) (Nat.le_add_left k N)).1
    dsimp only
    omega
  · have h := (hk (N + k) (Nat.le_add_left k N)).2
    exact h.trans (Finset.single_le_sum (f := fun i => n i (N + k)) (fun i _ => Nat.zero_le _)
      (Finset.mem_filter.mpr ⟨Finset.mem_univ i₀, hi₀⟩))

/-! ### 4. The law under the side condition -/

/-- **`SingleWeightLaw` for Gaussian noise, with the side condition `d N = p N + r`.** The
`align` field is `align_sw_of_gaussian` (G4, `RankR/SingleWeight/Het/Align.lean:271`) at the
edge `heteroEdgeR_of_gaussian_tail`; the `simpleIdx` field is `simpleIdx_sw_of_gaussian`.
Both need only `hex : ∃ i, w i ≠ 0`. -/
theorem singleWeightLaw_of_gaussian_aux [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hex : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) (hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N)
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    m.SingleWeightLaw w c γ z := by
  have hedge := m.heteroEdgeR_of_gaussian_tail w c hc hex hreg hG
  exact ⟨fun l k => m.align_sw_of_gaussian w c hc hex hreg hG hpd hedge hsep l k,
    fun l N => m.simpleIdx_sw_of_gaussian w hG hex (hrn N) l⟩

/-! ### 5. The law with no side condition -/

/-- **`SingleWeightLaw` for Gaussian noise** (the target of Track G). The side condition of
`singleWeightLaw_of_gaussian_aux` is met after the tail shift `exists_shift_sw`, and
`SingleWeightLaw.of_shift` carries the law back; `simpleIdx` is supplied directly for `m`,
since it reads one `N` at a time. The weights are unconstrained: at `r = 0` the law is
vacuous, and otherwise `EigSep` supplies the table of nonzero weight the shift needs. -/
theorem singleWeightLaw_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N)
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    m.SingleWeightLaw w c γ z := by
  rcases hsep.eq_zero_or_exists with hr | hex
  · exact SingleWeightLaw.of_r_eq_zero m w c γ z hr
  obtain ⟨i₀, hi₀⟩ := hex
  obtain ⟨k, p, hpd, -⟩ := m.exists_shift_sw hreg hi₀
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact SingleWeightLaw.of_shift k
    (singleWeightLaw_of_gaussian_aux (m.shift k) w c hc ⟨i₀, hi₀⟩
      (fun i => m.shift_Regime k i (hreg i)) (m.shift_JointGaussianNoise k hG) hpd
      (fun N => hrn (N + k)) hsep)
    (fun l N => m.simpleIdx_sw_of_gaussian w hG ⟨i₀, hi₀⟩ (hrn N) l)

/-! ### 6. The law with an existential shift, for the facades -/

/-- **`SingleWeightLaw` for Gaussian noise, on a shift of the model, with no `hrn`
hypothesis.** `exists_shift_sw` already supplies the row bound for the shifted model, so this
form needs no side condition at all; the two facades below apply
`SpikedModel.tendstoInProb_of_shift` to remove the shift from their own conclusions. -/
theorem singleWeightLaw_of_gaussian_shift [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    ∃ k, (m.shift k).SingleWeightLaw w c γ z := by
  rcases hsep.eq_zero_or_exists with hr | hex
  · exact ⟨0, SingleWeightLaw.of_r_eq_zero _ w c γ z hr⟩
  obtain ⟨i₀, hi₀⟩ := hex
  obtain ⟨k, p, hpd, hrn⟩ := m.exists_shift_sw hreg hi₀
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact ⟨k, singleWeightLaw_of_gaussian_aux (m.shift k) w c hc ⟨i₀, hi₀⟩
    (fun i => m.shift_Regime k i (hreg i)) (m.shift_JointGaussianNoise k hG) hpd hrn hsep⟩

/-! ### 7. The paper's facades -/

/-- **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`), second half, for
Gaussian noise: the single-weight stacksvd performance converges in probability to the
paper's sum. Layer 1 mirror: `prop_gen_rank_stacksvd_singleweight`
(`RankR/SingleWeight/Main.lean:88`). -/
theorem prop_gen_rank_stacksvd_singleweight_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    TendstoInProb μ (fun N ω => m.perfSW w N ω)
      (SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z) := by
  obtain ⟨k, law⟩ := m.singleWeightLaw_of_gaussian_shift w c hc hreg hG hsep
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact SpikedModel.tendstoInProb_of_shift (f := fun N ω => m.perfSW w N ω) k
    (prop_gen_rank_stacksvd_singleweight (m.shift k) w c γ z hsep law)

/-- **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`) for Gaussian noise, in
the paper's own inner-product form, one pair at a time. Layer 1 mirror:
`prop_gen_rank_stacksvd_singleweight_inner` (`RankR/SingleWeight/Main.lean:120`). -/
theorem prop_gen_rank_stacksvd_singleweight_inner_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l k : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2)
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
        (WithLp.ofLp (z l) k) ^ 2) := by
  obtain ⟨k', law⟩ := m.singleWeightLaw_of_gaussian_shift w c hc hreg hG hsep
  have : ∀ N, IsProbabilityMeasure (μ (N + k')) := isProbabilityMeasure_shift k'
  exact SpikedModel.tendstoInProb_of_shift
    (f := fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2) k'
    (prop_gen_rank_stacksvd_singleweight_inner (m.shift k') w c γ z law l k)

end UnalignedModelR

end StackedSVD
