/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.GeneralMain
import StackedSVD.RankR.Frobenius
import StackedSVD.RMT.Full

/-!
# The Frobenius form and the Gaussian facade at general `r_i`

Task TBfix, the response to `notes/archive/audit_independent_trackB_2026-09-02.md` findings W1, W2,
C1 and C3. Nothing here is new mathematics: every result transcribes a rank-one result of
`RankR/Frobenius.lean` to `UnalignedModelR`, or wires an existing chain.

## Content

1. **W2, the Gaussian facade at `r_i = 1`.** `SpikedModelR.toSpiked` reads a one-spike
   `SpikedModelR` as a `SpikedModel`. The map is total, because both `SpikedModel.hθ` and
   `SpikedModelR.hθnn` ask only `0 ≤ θ` (F8, 2026-09-05). With it,
   `SpikedModelR.tableLawR_of_gaussian` discharges `TableLawR c` from Gaussian noise through
   `singleTableLaw_of_gaussian` (`RMT/Full.lean`) and `tableLawR_of_singleTableLaw`
   (`RankR/General.lean`), and `UnalignedModelR.tableLawR_of_gaussian` does it for every table
   of a model at `rk = fun _ => 1`.
2. **W1, the Frobenius form.** `vhatSvdstackG` and `vhatSvdstackGW` are the paper's
   `V̂_svdstack = Ṽᵀ Q_r Λ_r^{-1/2}` (`main_paper.tex:1982, 2016`) at general `r_i`, and
   `frobSq_vhatSvdstackG`, `frobSq_vhatSvdstackGW` identify `perfRG` and `perfRGW` with
   `‖V̂ᵀ V‖_F²`. The corollaries `..._frobenius` ask the frame only on an event of probability
   tending to one; the `..._frobenius_eig` companions take the canonical frame `topEigMat` and
   carry no frame hypothesis, because `topGap_gramG_whp` and `topGap_gramWG_whp` supply the
   gap. Mirrors: `Frobenius.lean:575-600, 700-861`.
3. **W2 (iii), the Gaussian corollaries.** `..._gaussian_one` at `rk = fun _ => 1`: the two
   Layer 1 theorems with `law` replaced by `hc`, `hreg` and `hG : m.JointGaussianNoise`.
4. **W2 (iv), the model bridge.** `UnalignedModelR.toUnaligned` reads a model at
   `rk = fun _ => 1` as an `UnalignedModel`, and `perfRG_one_eq_perfR` shows the two
   performance functionals agree along `i ↦ flat i 0`. Mirror of section 6 of
   `RankR/General.lean`, which does the same for the limits.
5. **C1 and C3, the hypothesis trims.**
   `thm_gen_rank_weight_svdstak_general_r_of_rank` drops `hrr : r ≤ rtot rk`, which follows
   from `hrankB` by `Matrix.rank_le_card_height`.
   `thm_gen_rank_weight_svdstak_general_r_paper` takes the paper's own two hypotheses,
   `β_ij > 0` and `Rank(∑_i R_i R_iᵀ) = r`, in place of `hrankB`.

Every theorem here takes `hG : UnalignedModelR.IndepNoise` where Layer 1 takes it (finding W3,
`RankR/GramR.lean`), except the `_gaussian_one` corollaries, which take
`JointGaussianNoise` and pass `JointGaussianNoise.indepNoise`.

STATUS 2026-09-02: no `sorry`, no `axiom`; see `notes/archive/agent_reports/d32_TBfix.md`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### 1. W2: a one-spike `SpikedModelR` read as a `SpikedModel` -/

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- The columns of `U N` are unit vectors. This is the model field `hU` in the shape that
`SpikedModel.hu` takes. Mirror on the left factor of `norm_col` (`RankR/GeneralMain.lean`). -/
theorem norm_ucol (t : SpikedModelR μ n d rk) (N : ℕ) (k : Fin rk) :
    ‖(WithLp.toLp 2 fun l => t.U N l k : EuclideanSpace ℝ (Fin (n N)))‖ = 1 := by
  have hUU := t.hU N
  have hdot : ∑ a, (t.U N) a k * (t.U N) a k = ((t.U N)ᵀ * t.U N) k k := by
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun a _ => by rw [Matrix.transpose_apply]
  have hsq : ‖(WithLp.toLp 2 fun l => t.U N l k : EuclideanSpace ℝ (Fin (n N)))‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct]
    change ∑ a, (t.U N) a k * (t.U N) a k = 1
    rw [hdot, hUU, Matrix.one_apply_eq]
  rw [← Real.sqrt_sq (norm_nonneg (WithLp.toLp 2 fun l => t.U N l k :
    EuclideanSpace ℝ (Fin (n N)))), hsq, Real.sqrt_one]

/-- A rank-one `SpikedModelR` read as a `SpikedModel`: the one spike, the one column of `U`
and of `V`, the same noise. The map is total in this direction, because `SpikedModel` accepts
`θ = m.θ 0 > 0`; the opposite direction does not exist, which is why
`tableLawR_of_singleTableLaw` (`RankR/General.lean`) is stated on two models. -/
noncomputable def toSpiked (m : SpikedModelR μ n d 1) : SpikedModel μ n d where
  θ := m.θ 0
  u N := WithLp.toLp 2 fun l => m.U N l 0
  v N := m.col N 0
  Z := m.Z
  hθ := m.hθnn 0
  hn := m.hn
  hd := m.hd
  hu N := m.norm_ucol N 0
  hv N := m.norm_col N 0
  hZ := m.hZ

/-- The two models carry the same data matrix: `U Θ Vᵀ = θ_1 u_1 v_1ᵀ` on `Fin 1`, and the
noise is the same function. -/
theorem toSpiked_X (m : SpikedModelR μ n d 1) (N : ℕ) (ω : Ω N) :
    m.X N ω = m.toSpiked.X N ω := by
  have hsig : m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ
      = m.toSpiked.θ • Matrix.vecMulVec (WithLp.ofLp (m.toSpiked.u N))
          (WithLp.ofLp (m.toSpiked.v N)) := by
    ext a b
    have h1 : (m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ) a b
        = m.U N a 0 * m.θ 0 * m.V N b 0 := by
      rw [Matrix.mul_apply, Fin.sum_univ_one, Matrix.mul_diagonal, Matrix.transpose_apply]
    rw [h1]
    change _ = m.θ 0 * (m.U N a 0 * m.V N b 0)
    ring
  exact congrArg (fun x => x + m.E N ω) hsig

/-- The regime is a statement about `n` and `d` only, so it transfers by `rfl`. -/
theorem toSpiked_Regime (m : SpikedModelR μ n d 1) (c : ℝ) :
    m.toSpiked.Regime c ↔ m.Regime c := Iff.rfl

/-- The two models have the same noise function, so Gaussian noise transfers by `rfl`. -/
theorem toSpiked_GaussianNoise (m : SpikedModelR μ n d 1) :
    m.toSpiked.GaussianNoise ↔ m.GaussianNoise := Iff.rfl

/-- **`prop:single_table` at `r_i = 1` under Gaussian noise.** The chain the audit asks for
(finding W2): `singleTableLaw_of_gaussian` (`RMT/Full.lean`) on the facade, then
`tableLawR_of_singleTableLaw` (`RankR/General.lean`). So at one spike per table the black box
`TableLawR` is discharged, not assumed. -/
theorem tableLawR_of_gaussian [∀ N, IsProbabilityMeasure (μ N)] (m : SpikedModelR μ n d 1)
    {c : ℝ} (hc : 0 < c) (hreg : m.Regime c) (hG : m.GaussianNoise) : m.TableLawR c :=
  tableLawR_of_singleTableLaw m m.toSpiked m.toSpiked_X (fun _ => rfl) rfl
    (SpikedModel.singleTableLaw_of_gaussian hc m.toSpiked ((m.toSpiked_Regime c).mpr hreg)
      (m.toSpiked_GaussianNoise.mpr hG))

end SpikedModelR

/-! ### 1b. C1 and C3: two facts about the rank of `B_R` -/

section BlockRank

variable {M r : ℕ} {rk : Fin M → ℕ}

/-- Finding C1: `r ≤ r̃` follows from `rank B_R = r`, because a matrix has rank at most its
number of rows (`Matrix.rank_le_card_height`). -/
theorem rank_le_card_of_rank_BBlock {β : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} (hrankB : (BBlock β R).rank = r) :
    r ≤ Fintype.card (Fin (rtot rk)) := by
  rw [← hrankB]
  exact Matrix.rank_le_card_height _

/-- `∑_i R_i R_iᵀ` is the sum of the rank-one terms of `R_stack` over the flat spike index. -/
theorem sum_mul_transpose_eq_sum_vecMulVec
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    ∑ i, R i * (R i)ᵀ
      = ∑ p : Fin (rtot rk),
          Matrix.vecMulVec (WithLp.ofLp (Rcol R p)) (WithLp.ofLp (Rcol R p)) := by
  ext k l
  rw [Matrix.sum_apply, Matrix.sum_apply]
  have hL : ∀ i : Fin M, (R i * (R i)ᵀ) k l = ∑ j : Fin (rk i), R i k j * R i l j := by
    intro i
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun j _ => by rw [Matrix.transpose_apply]
  have hR : ∀ p : Fin (rtot rk),
      Matrix.vecMulVec (WithLp.ofLp (Rcol R p)) (WithLp.ofLp (Rcol R p)) k l
        = R (blk p).1 k (blk p).2 * R (blk p).1 l (blk p).2 := fun _ => rfl
  simp only [hL, hR, blk]
  exact (sum_stack_index fun i j => R i k j * R i l j).symm

/-- Finding C3: the paper's rank condition `Rank(∑_i R_i R_iᵀ) = r` (`main_paper.tex:757`) and
the Lean form `(Rstack (Rcol R)).rank = r` are the same condition. -/
theorem rank_sum_mul_transpose (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    (∑ i, R i * (R i)ᵀ).rank = (Rstack (Rcol R)).rank := by
  rw [sum_mul_transpose_eq_sum_vecMulVec R, rank_sum_vecMulVec]

/-- Finding C3, the step the paper-literal corollary needs: `β_ij > 0` at every spike plus the
paper's rank condition give `rank B_R = r`. The converse fails, so
`thm_gen_rank_weight_svdstak_general_r` is the stronger theorem (`Weighted.lean:530`, D16). -/
theorem rank_BBlock_of_paper {β : (i : Fin M) → Fin (rk i) → ℝ}
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (hβpos : ∀ i j, 0 < β i j)
    (hrank : (∑ i, R i * (R i)ᵀ).rank = r) : (BBlock β R).rank = r := by
  rw [BBlock_eq_BR]
  refine rank_BR_of_ne_zero (Rcol R) (fun p => ne_of_gt (hβpos (blk p).1 (blk p).2)) ?_
  rw [← rank_sum_mul_transpose R]
  exact hrank

end BlockRank

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- Every table of a model with one spike per table satisfies `TableLawR` under Gaussian
noise. `gaussianNoise_of_joint_piR` (`RankR/GramR.lean`) takes the marginal of the product
law. -/
theorem tableLawR_of_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r fun _ => 1) {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∀ i, (m.tbl i).TableLawR (c i) :=
  fun i => (m.tbl i).tableLawR_of_gaussian (hc i) (hreg i)
    (gaussianNoise_of_joint_piR m.tbl hG i)

/-! ### 2. W1: `V̂_svdstack` and the Frobenius form, unweighted -/

/-- `V̂_svdstack` at general `r_i`, written as the paper writes it (`main_paper.tex:1982`):
`V̂_svdstack = Ṽᵀ Q_r Λ_r^{-1/2}` with `(Q, λ)` a top-`r` eigenframe of `Ṽ Ṽᵀ`. Mirror:
`UnalignedModel.vhatSvdstack`. -/
noncomputable def vhatSvdstackG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  (m.VtG N ω)ᵀ * Q * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹

/-- The paper's display: `V̂_svdstackᵀ V = Λ^{-1/2} Q_rᵀ Ṽ V`. -/
theorem vhatSvdstackG_transpose_mul_V (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) :
    (m.vhatSvdstackG N ω Q lam)ᵀ * m.V N
      = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) * (Qᵀ * m.VtVG N ω) := by
  simp only [vhatSvdstackG, UnalignedModelR.VtVG, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.mul_assoc]

/-- **`perfRG` is the paper's `‖V̂_svdstackᵀ V‖_F²`** (finding W1). The identity needs a top-`r`
eigenframe of `Ṽ Ṽᵀ` with nonnegative eigenvalues; `topGap_gramG_whp` and `posSemidef_gramG`
give it with probability tending to one. -/
theorem frobSq_vhatSvdstackG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    {Q : Matrix (Fin (rtot rk)) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gramG N ω) (m.isHermitian_gramG N ω) Q lam)
    (hlam : ∀ j, 0 ≤ lam j) :
    frobSq ((m.vhatSvdstackG N ω Q lam)ᵀ * m.V N) = m.perfRG N ω := by
  rw [m.vhatSvdstackG_transpose_mul_V N ω Q lam, UnalignedModelR.perfRG]
  exact frobSq_eq_trace_of_eigFrame hQ hlam (m.VtVG N ω)

/-- The Gram matrix `Ṽ Ṽᵀ` is positive semidefinite. -/
theorem posSemidef_gramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    (m.gramG N ω).PosSemidef := by
  simpa [UnalignedModelR.gramG] using Matrix.posSemidef_self_mul_conjTranspose (m.VtG N ω)

/-- The top-`r` gap of `Ṽ Ṽᵀ` holds with probability tending to one: the entries converge to
`A_{β,R}` (`gramR_general`), which has a gap by `hgap`. Mirror: `topGap_gram_whp`. At `r = 0`
the set is empty: `TopGap A hA 0` asks nothing, because no index is below `0`. The root lemma
`specTop_simple_whp_of_tendsto` still takes `0 < r`. -/
theorem topGap_gramG_whp (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hrp : r ≤ Fintype.card (Fin (rtot rk)))
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.gramG N ω) (m.isHermitian_gramG N ω) r}) atTop
      (𝓝 0) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · have hempty : ∀ N,
        {ω : Ω N | ¬ TopGap (m.gramG N ω) (m.isHermitian_gramG N ω) 0} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact fun k _ hk _ => absurd hk (Nat.not_lt_zero _)
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds
  · exact specTop_simple_whp_of_tendsto (isHermitian_ABlock β m.R)
      (fun N ω => m.isHermitian_gramG N ω)
      (fun p q => m.gramR_general c β hβdef law hG p q) hr hrp hgap

/-- **`prop:general_rank_unweighted_svdstack`** (`main_paper.tex:799`) at general `r_i`, on the
paper's own quantity `‖V̂_svdstackᵀ V‖_F²`. `(Q N ω, λ N ω)` is any selection that is a top-`r`
eigenframe of `Ṽ Ṽᵀ` with nonnegative eigenvalues on an event whose probability tends to one.
Mirror: `prop_general_rank_unweighted_svdstack_frobenius`. -/
theorem prop_general_rank_unweighted_svdstack_general_frobenius
    (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hrr : r ≤ rtot rk)
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (Q : ∀ N, Ω N → Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gramG N ω) (m.isHermitian_gramG N ω)
      (Q N ω) (lam N ω) ∧ ∀ j, 0 ≤ lam N ω j)}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackG N ω (Q N ω) (lam N ω))ᵀ * m.V N))
      (limitRG β m.R) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.prop_general_rank_unweighted_svdstack_general c β hc hβdef hrr hgap law hG)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hQ (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgood
  exact hω (m.frobSq_vhatSvdstackG N ω hgood.1 hgood.2)

/-- The same at the canonical frame `(Q_r, Λ_r) = (topEigMat, topEigVal)` of `Ṽ Ṽᵀ`, with **no
frame hypothesis**: the hypotheses are exactly those of
`prop_general_rank_unweighted_svdstack_general`. So the paper's `‖V̂_svdstackᵀ V‖_F²`
converges for a definite estimator at general `r_i`. Second cleanup pass, 2026-09-02: `0 < r`
is gone; `topGap_gramG_whp` now covers `r = 0`. -/
theorem prop_general_rank_unweighted_svdstack_general_frobenius_eig
    (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hrp : r ≤ Fintype.card (Fin (rtot rk)))
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackG N ω
        (topEigMat (m.isHermitian_gramG N ω) hrp)
        (topEigVal (m.isHermitian_gramG N ω) hrp))ᵀ * m.V N))
      (limitRG β m.R) := by
  refine m.prop_general_rank_unweighted_svdstack_general_frobenius c β hc hβdef
    (by simpa using hrp) hgap law hG _ _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topGap_gramG_whp c β hβdef hrp hgap law hG) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω ⟨isTopEigFrame_topEigMat (m.isHermitian_gramG N ω) hrp hgapN,
    fun j => topEigVal_nonneg hrp (m.posSemidef_gramG N ω) j⟩

/-! ### 3. W1: `V̂_svdstack(W)` and the Frobenius form, weighted -/

/-- `V̂_svdstack(W)` at general `r_i` (`main_paper.tex:2016`): `V̂ = Ṽ_Wᵀ Q_r Λ_r^{-1/2}` with
`(Q, λ)` a top-`r` eigenframe of `Ṽ_W Ṽ_Wᵀ`. -/
noncomputable def vhatSvdstackGW (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  (m.VtWG W N ω)ᵀ * Q * Matrix.diagonal fun j => (Real.sqrt (lam j))⁻¹

theorem vhatSvdstackGW_transpose_mul_V (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) :
    (m.vhatSvdstackGW W N ω Q lam)ᵀ * m.V N
      = Matrix.diagonal (fun j => (Real.sqrt (lam j))⁻¹) * (Qᵀ * m.VtVWG W N ω) := by
  simp only [vhatSvdstackGW, UnalignedModelR.VtVWG, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.diagonal_transpose, Matrix.mul_assoc]

/-- **`perfRGW W` is the paper's `‖V̂_svdstack(W)ᵀ V‖_F²`** (finding W1). -/
theorem frobSq_vhatSvdstackGW (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    {Q : Matrix (Fin (rtot rk)) (Fin r) ℝ} {lam : Fin r → ℝ}
    (hQ : IsTopEigFrame (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) Q lam)
    (hlam : ∀ j, 0 ≤ lam j) :
    frobSq ((m.vhatSvdstackGW W N ω Q lam)ᵀ * m.V N) = m.perfRGW W N ω := by
  rw [m.vhatSvdstackGW_transpose_mul_V W N ω Q lam, UnalignedModelR.perfRGW]
  exact frobSq_eq_trace_of_eigFrame hQ hlam (m.VtVWG W N ω)

/-- The weighted Gram matrix `Ṽ_W Ṽ_Wᵀ` is positive semidefinite. -/
theorem posSemidef_gramWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    (m.gramWG W N ω).PosSemidef := by
  simpa [UnalignedModelR.gramWG] using Matrix.posSemidef_self_mul_conjTranspose (m.VtWG W N ω)

/-- `Ṽ_W Ṽ_Wᵀ → W A_{β,R} Wᵀ` entrywise in probability. Each entry is a fixed real linear
combination of the entries of `Ṽ Ṽᵀ` (`mul_mul_transpose_apply`). The same step runs inside
`thm_gen_rank_weight_svdstak_general_r_conv`; it is named here for `topGap_gramWG_whp`.
Mirror: `UnalignedModel.gramRW` (`RankR/Frobenius.lean`). -/
theorem gramRWG (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (p q : Fin (rtot rk)) :
    TendstoInProb μ (fun N ω => m.gramWG W N ω p q) (ABlockW W β m.R p q) := by
  have hconv0 : TendstoInProbPi μ
      (fun N ω => fun t : Fin (rtot rk) × Fin (rtot rk) => m.gramG N ω t.1 t.2)
      (fun t : Fin (rtot rk) × Fin (rtot rk) => ABlock β m.R t.1 t.2) :=
    fun t => m.gramR_general c β hβdef law hG t.1 t.2
  have hcont : Continuous fun z : Fin (rtot rk) × Fin (rtot rk) → ℝ =>
      ∑ b : Fin (rtot rk), ∑ a : Fin (rtot rk), W p a * z (a, b) * W q b :=
    continuous_finsetSum _ fun b _ => continuous_finsetSum _ fun a _ =>
      (continuous_const.mul (continuous_apply _)).mul continuous_const
  have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
  have hlim : ABlockW W β m.R p q
      = ∑ b : Fin (rtot rk), ∑ a : Fin (rtot rk), W p a * ABlock β m.R a b * W q b :=
    UnalignedModel.mul_mul_transpose_apply W (ABlock β m.R) p q
  rw [hlim]
  refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
  change ∑ b : Fin (rtot rk), ∑ a : Fin (rtot rk), W p a * m.gramG N ω a b * W q b
      = m.gramWG W N ω p q
  rw [m.gramWG_eq]
  exact (UnalignedModel.mul_mul_transpose_apply W (m.gramG N ω) p q).symm

/-- The top-`r` gap of `Ṽ_W Ṽ_Wᵀ` holds with probability tending to one, from the gap of
`W A_{β,R} Wᵀ`. Weighted twin of `topGap_gramG_whp`, `r = 0` included. -/
theorem topGap_gramWG_whp (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hrp : r ≤ Fintype.card (Fin (rtot rk)))
    (hgapW : TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) r})
      atTop (𝓝 0) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · have hempty : ∀ N,
        {ω : Ω N | ¬ TopGap (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) 0} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact fun k _ hk _ => absurd hk (Nat.not_lt_zero _)
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds
  · exact specTop_simple_whp_of_tendsto (isHermitian_ABlockW W β m.R)
      (fun N ω => m.isHermitian_gramWG W N ω)
      (fun p q => m.gramRWG c β W hβdef law hG p q) hr hrp hgapW

/-- **`thm:gen_rank_weight_svdstak`** for one admissible `W` (`main_paper.tex:893`) at general
`r_i`, on the paper's own quantity. The frame is asked only on an event whose probability
tends to one. Mirror: `thm_gen_rank_weight_svdstak_frobenius`. -/
theorem thm_gen_rank_weight_svdstak_general_r_conv_frobenius
    (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hr : 0 < r) (hrr : r ≤ rtot rk)
    (hgapW : TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r)
    (hposW : 0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (Q : ∀ N, Ω N → Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gramWG W N ω)
      (m.isHermitian_gramWG W N ω) (Q N ω) (lam N ω) ∧ ∀ j, 0 ≤ lam N ω j)}) atTop (𝓝 0)) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackGW W N ω (Q N ω) (lam N ω))ᵀ * m.V N))
      (limitRGW W β m.R) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.thm_gen_rank_weight_svdstak_general_r_conv c β W hβdef hr hrr hgapW hposW law hG)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hQ (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgood
  exact hω (m.frobSq_vhatSvdstackGW W N ω hgood.1 hgood.2)

/-- The same at the canonical frame of `Ṽ_W Ṽ_Wᵀ`, with **no frame hypothesis**: the
hypotheses are those of `thm_gen_rank_weight_svdstak_general_r_conv`. -/
theorem thm_gen_rank_weight_svdstak_general_r_conv_frobenius_eig
    (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hr : 0 < r)
    (hrp : r ≤ Fintype.card (Fin (rtot rk)))
    (hgapW : TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r)
    (hposW : 0 < (isHermitian_ABlockW W β m.R).eigenvalues₀ ⟨r - 1, by omega⟩)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackGW W N ω
        (topEigMat (m.isHermitian_gramWG W N ω) hrp)
        (topEigVal (m.isHermitian_gramWG W N ω) hrp))ᵀ * m.V N))
      (limitRGW W β m.R) := by
  refine m.thm_gen_rank_weight_svdstak_general_r_conv_frobenius c β W hβdef hr
    (by simpa using hrp) hgapW hposW law hG _ _ ?_
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topGap_gramWG_whp c β W hβdef hrp hgapW law hG) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  intro hgapN
  exact hω ⟨isTopEigFrame_topEigMat (m.isHermitian_gramWG W N ω) hrp hgapN,
    fun j => topEigVal_nonneg hrp (m.posSemidef_gramWG W N ω) j⟩

/-- **`thm:gen_rank_weight_svdstak`, first display** (`main_paper.tex:893`) at general `r_i`,
on the paper's own quantity `‖V̂_svdstack(W⋆)ᵀ V‖_F²` and at the canonical frame of
`Ṽ_{W⋆} Ṽ_{W⋆}ᵀ`. The hypotheses are those of `thm_gen_rank_weight_svdstak_general_r`; the
eigengap and the positivity at `W⋆` come from `rank B_R = r`, as they do there. -/
theorem thm_gen_rank_weight_svdstak_general_r_frobenius_eig (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hr : 0 < r)
    (hrp : r ≤ Fintype.card (Fin (rtot rk))) (hrankB : (BBlock β m.R).rank = r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackGW (optWG β) N ω
        (topEigMat (m.isHermitian_gramWG (optWG β) N ω) hrp)
        (topEigVal (m.isHermitian_gramWG (optWG β) N ω) hrp))ᵀ * m.V N))
      (limitOptG β m.R hrp) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hrankB' : (BR (betaFlat β) (Rcol m.R)).rank = r := by
    rw [← BBlock_eq_BR]
    exact hrankB
  have hgapOpt : TopGap (ABlockW (optWG β) β m.R) (isHermitian_ABlockW (optWG β) β m.R) r :=
    topGap_optWR_of_rankBR (β := betaFlat β) (Rcol m.R) h0 h1 hrankB'
  have hpdOpt : (ABlockW (optWG β) β m.R).PosDef :=
    abetaRW_optWR_posDef (β := betaFlat β) (Rcol m.R) h0 h1
  have hposOpt : 0 < (isHermitian_ABlockW (optWG β) β m.R).eigenvalues₀
      ⟨r - 1, by omega⟩ := by
    rw [← eigenvalues_eigIdx]
    exact hpdOpt.eigenvalues_pos _
  have h := m.thm_gen_rank_weight_svdstak_general_r_conv_frobenius_eig c β (optWG β) hβdef hr
    hrp hgapOpt hposOpt law hG
  have hval : limitRGW (optWG β) β m.R = limitOptG β m.R hrp := by
    rw [limitRGW_eq_limitRW, limitOptG_eq_limitROpt, optWG_eq_optWR,
      limitRW_optWR (betaFlat β) (Rcol m.R) h0 h1 hrp]
  rwa [hval] at h

/-! ### 4. C1 and C3: the two hypothesis trims of `thm:gen_rank_weight_svdstak` -/

/-- **Finding C1.** `thm_gen_rank_weight_svdstak_general_r` without `hrr : r ≤ r̃`: the bound
is already inside `hrankB`, because `B_R` has `r̃` rows. -/
theorem thm_gen_rank_weight_svdstak_general_r_of_rank (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hr : 0 < r) (hrankB : (BBlock β m.R).rank = r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB)) ∧
      ∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by
            have h2 := rank_le_card_of_rank_BBlock hrankB
            omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB) :=
  m.thm_gen_rank_weight_svdstak_general_r c β hc hβdef hr
    (by simpa using rank_le_card_of_rank_BBlock hrankB) hrankB law hG

/-- **Finding C3, the paper-literal corollary.** `thm_gen_rank_weight_svdstak_general_r` with
the paper's own two hypotheses (`main_paper.tex:757, 893`): `β_ij > 0` at every spike and
`Rank(∑_i R_i R_iᵀ) = r`. They imply `rank B_R = r` and are strictly stronger, so this is the
weaker theorem of the two. -/
theorem thm_gen_rank_weight_svdstak_general_r_paper (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hr : 0 < r) (hβpos : ∀ i j, 0 < β i j)
    (hrank : (∑ i, m.R i * (m.R i)ᵀ).rank = r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R
          (rank_le_card_of_rank_BBlock (rank_BBlock_of_paper m.R hβpos hrank))) ∧
      ∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by
            have h2 := rank_le_card_of_rank_BBlock
              (rank_BBlock_of_paper (β := β) m.R hβpos hrank)
            omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R
            (rank_le_card_of_rank_BBlock (rank_BBlock_of_paper m.R hβpos hrank)) :=
  m.thm_gen_rank_weight_svdstak_general_r_of_rank c β hc hβdef hr
    (rank_BBlock_of_paper m.R hβpos hrank) law hG

/-! ### 5. W2 (iii): the Gaussian corollaries at one spike per table

At `rk = fun _ => 1` the black box `TableLawR` is discharged by
`UnalignedModelR.tableLawR_of_gaussian`, so these four statements are unconditional for
Gaussian noise in the proportional regime. -/

/-- **`prop:general_rank_unweighted_svdstack` under Gaussian noise**, at one spike per table.
No `TableLawR` hypothesis: `hreg` and `hG` produce it. -/
theorem prop_general_rank_unweighted_svdstack_general_gaussian_one
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r fun _ => 1)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin 1 → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hrr : r ≤ rtot fun _ : Fin M => (1 : ℕ))
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRG N ω) (limitRG β m.R) :=
  m.prop_general_rank_unweighted_svdstack_general c β hc hβdef hrr hgap
    (m.tableLawR_of_gaussian hc hreg hG) hG.indepNoise

/-- The same on the paper's own quantity `‖V̂_svdstackᵀ V‖_F²`, at the canonical frame. -/
theorem prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian_one
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r fun _ => 1)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin 1 → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hrp : r ≤ Fintype.card (Fin (rtot fun _ : Fin M => (1 : ℕ))))
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackG N ω
        (topEigMat (m.isHermitian_gramG N ω) hrp)
        (topEigVal (m.isHermitian_gramG N ω) hrp))ᵀ * m.V N))
      (limitRG β m.R) :=
  m.prop_general_rank_unweighted_svdstack_general_frobenius_eig c β hc hβdef hrp hgap
    (m.tableLawR_of_gaussian hc hreg hG) hG.indepNoise

/-- **`thm:gen_rank_weight_svdstak` under Gaussian noise**, at one spike per table. -/
theorem thm_gen_rank_weight_svdstak_general_r_gaussian_one
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r fun _ => 1)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin 1 → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r)
    (hrankB : (BBlock β m.R).rank = r) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB)) ∧
      ∀ W : Matrix (Fin (rtot fun _ : Fin M => (1 : ℕ))) (Fin (rtot fun _ : Fin M => (1 : ℕ))) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by
            have h2 := rank_le_card_of_rank_BBlock hrankB
            omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB) :=
  m.thm_gen_rank_weight_svdstak_general_r_of_rank c β hc hβdef hr hrankB
    (m.tableLawR_of_gaussian hc hreg hG) hG.indepNoise

/-- The same on the paper's own quantity, at `W⋆` and at the canonical frame. -/
theorem thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian_one
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r fun _ => 1)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin 1 → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r)
    (hrp : r ≤ Fintype.card (Fin (rtot fun _ : Fin M => (1 : ℕ))))
    (hrankB : (BBlock β m.R).rank = r) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackGW (optWG β) N ω
        (topEigMat (m.isHermitian_gramWG (optWG β) N ω) hrp)
        (topEigVal (m.isHermitian_gramWG (optWG β) N ω) hrp))ᵀ * m.V N))
      (limitOptG β m.R hrp) :=
  m.thm_gen_rank_weight_svdstak_general_r_frobenius_eig c β hc hβdef hr hrp hrankB
    (m.tableLawR_of_gaussian hc hreg hG) hG.indepNoise

end UnalignedModelR

/-! ### 6. W2 (iv): the model bridge at one spike per table

At `rk = fun _ => 1` an `UnalignedModelR` is an `UnalignedModel`, and the two performance
functionals agree. The flat index type `Fin (∑ i : Fin M, 1)` is not definitionally `Fin M`,
so the statements read the block objects along `i ↦ flat i 0`, exactly as section 6 of
`RankR/General.lean` does for the limits. -/

/-- `vEig` does not see the `IsHermitian` proof, so it transports along an equality of
matrices. Mirror: `eigenvalues₀_congr_mat` (`LinAlg/Eigen.lean`). -/
theorem vEig_congr_mat {p : ℕ} {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (j : ℕ) : vEig A hA j = vEig B hB j := by
  subst h; rfl

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The alignment vector of table `i` at one spike: the single column of `R_i`. -/
noncomputable def Rone (m : UnalignedModelR μ M n d r fun _ => 1) (i : Fin M) :
    EuclideanSpace ℝ (Fin r) := WithLp.toLp 2 fun k => m.R i k 0

/-- `R_i` is a unit vector, the field `hR` of `UnalignedModel`. -/
theorem norm_Rone (m : UnalignedModelR μ M n d r fun _ => 1) (i : Fin M) : ‖m.Rone i‖ = 1 := by
  have hRR := m.hR i
  have hdot : ∑ k, m.R i k 0 * m.R i k 0 = ((m.R i)ᵀ * m.R i) 0 0 := by
    rw [Matrix.mul_apply]
    exact Finset.sum_congr rfl fun k _ => by rw [Matrix.transpose_apply]
  have hsq : ‖m.Rone i‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct]
    change ∑ k, m.R i k 0 * m.R i k 0 = 1
    rw [hdot, hRR, Matrix.one_apply_eq]
  rw [← Real.sqrt_sq (norm_nonneg (m.Rone i)), hsq, Real.sqrt_one]

/-- The model bridge of finding W2 (iv): an `UnalignedModelR` with one spike per table read as
an `UnalignedModel`. Every field is the corresponding field of the facade
`SpikedModelR.toSpiked`. -/
noncomputable def toUnaligned (m : UnalignedModelR μ M n d r fun _ => 1) :
    UnalignedModel μ M n d r where
  tbl i := (m.tbl i).toSpiked
  V := m.V
  R := m.Rone
  hV := m.hV
  hR := m.norm_Rone
  hv i N := by
    have h := m.hv i N
    refine congrArg (WithLp.toLp 2) (funext fun l => ?_)
    change (m.tbl i).V N l 0 = _
    rw [h, Matrix.mul_apply]
    simp [Matrix.mulVec, dotProduct, Rone]

/-- The two models carry the same data matrix in every table. -/
theorem toUnaligned_X (m : UnalignedModelR μ M n d r fun _ => 1) (i : Fin M) (N : ℕ)
    (ω : Ω N) : (m.tbl i).X N ω = (m.toUnaligned.tbl i).X N ω :=
  (m.tbl i).toSpiked_X N ω

/-- The two models carry the same per-table Gram matrix. -/
theorem toUnaligned_tableGram (m : UnalignedModelR μ M n d r fun _ => 1) (i : Fin M) (N : ℕ)
    (ω : Ω N) : m.tableGramG i N ω = m.toUnaligned.tableGram i N ω := by
  rw [UnalignedModelR.tableGramG, UnalignedModel.tableGram, m.toUnaligned_X i N ω]

/-- The per-table estimate at the single spike index is the rank-one estimate. -/
theorem toUnaligned_vhat (m : UnalignedModelR μ M n d r fun _ => 1) (i : Fin M) (N : ℕ)
    (ω : Ω N) : m.vhatG i 0 N ω = m.toUnaligned.vhat i N ω := by
  have hgr := m.toUnaligned_tableGram i N ω
  have hv0 : vEig (m.tableGramG i N ω) (m.isHermitian_tableGramG i N ω) ((0 : Fin 1) : ℕ)
      = vMax (m.toUnaligned.tableGram i N ω) (m.toUnaligned.isHermitian_tableGram i N ω) := by
    rw [show ((0 : Fin 1) : ℕ) = 0 from rfl,
      vEig_congr_mat (m.isHermitian_tableGramG i N ω)
        (m.toUnaligned.isHermitian_tableGram i N ω) hgr, vEig_zero]
  rw [UnalignedModelR.vhatG, UnalignedModel.vhat, hv0]
  rfl

/-- `Ṽ` of the two models, along `i ↦ flat i 0`. -/
theorem VtG_one_eq_Vt (m : UnalignedModelR μ M n d r fun _ => 1) (N : ℕ) (ω : Ω N) :
    (m.VtG N ω).submatrix (fun i : Fin M => flat i 0) id = m.toUnaligned.Vt N ω := by
  ext i l
  simp only [Matrix.submatrix_apply, id_eq, UnalignedModelR.VtG, Matrix.of_apply,
    UnalignedModel.Vt, blk_flat]
  rw [m.toUnaligned_vhat i N ω]

/-- `Ṽ Ṽᵀ` of the two models, along `i ↦ flat i 0`. -/
theorem gramG_one_eq_gram (m : UnalignedModelR μ M n d r fun _ => 1) (N : ℕ) (ω : Ω N) :
    m.toUnaligned.gram N ω
      = (m.gramG N ω).submatrix (fun i : Fin M => flat i 0) fun i : Fin M => flat i 0 := by
  have h := Matrix.submatrix_mul_equiv (m.VtG N ω) (m.VtG N ω)ᵀ
    (fun i : Fin M => flat i 0) (Equiv.refl (Fin (d N))) fun i : Fin M => flat i 0
  simp only [Equiv.coe_refl] at h
  rw [UnalignedModel.gram, ← m.VtG_one_eq_Vt N ω, Matrix.transpose_submatrix, h,
    UnalignedModelR.gramG]

/-- `Ṽ V` of the two models, along `i ↦ flat i 0`. -/
theorem VtVG_one_eq_VtV (m : UnalignedModelR μ M n d r fun _ => 1) (N : ℕ) (ω : Ω N) :
    m.toUnaligned.VtV N ω = (m.VtVG N ω).submatrix (fun i : Fin M => flat i 0) id := by
  have h := Matrix.submatrix_mul_equiv (m.VtG N ω) (m.V N)
    (fun i : Fin M => flat i 0) (Equiv.refl (Fin (d N))) (id : Fin r → Fin r)
  simp only [Equiv.coe_refl, Matrix.submatrix_id_id] at h
  rw [UnalignedModel.VtV, ← m.VtG_one_eq_Vt N ω]
  change (m.VtG N ω).submatrix (fun i : Fin M => flat i 0) id * m.V N
      = (m.VtVG N ω).submatrix (fun i : Fin M => flat i 0) id
  rw [h, UnalignedModelR.VtVG]

/-- **Finding W2 (iv).** The performance of svdstack at one spike per table is the rank-one
performance: `perfRG = perfR` on the bridged model. The route is the one of
`limitRG_one_eq_limitR` (`RankR/General.lean`): the trace form does not see a simultaneous
relabeling of the flat index, and `i ↦ flat i 0` is a bijection here. -/
theorem perfRG_one_eq_perfR (m : UnalignedModelR μ M n d r fun _ => 1) (N : ℕ) (ω : Ω N) :
    m.perfRG N ω = m.toUnaligned.perfR N ω := by
  obtain ⟨e, he⟩ : ∃ e : Fin M ≃ Fin (rtot fun _ : Fin M => 1),
      ∀ i, e i = flat (rk := fun _ : Fin M => 1) i 0 := by
    refine ⟨Equiv.ofBijective (fun i : Fin M => flat (rk := fun _ : Fin M => 1) i 0) ?_,
      fun i => rfl⟩
    refine (Fintype.bijective_iff_injective_and_card _).mpr ⟨?_, by simp [rtot]⟩
    intro i i' h
    have h2 : blk (rk := fun _ : Fin M => 1) (flat i 0) = blk (flat i' 0) := congrArg blk h
    rw [blk_flat, blk_flat] at h2
    exact congrArg Sigma.fst h2
  have hcoe : (⇑e) = fun i : Fin M => flat (rk := fun _ : Fin M => 1) i 0 := funext he
  rw [UnalignedModelR.perfRG, UnalignedModel.perfR]
  refine (trace_specInvTop_congr_submatrix (m.gramG N ω) (m.isHermitian_gramG N ω)
    (m.VtVG N ω) e (m.toUnaligned.gram N ω) (m.toUnaligned.isHermitian_gram N ω)
    (m.toUnaligned.VtV N ω) ?_ ?_ r).symm
  · rw [hcoe]
    exact m.gramG_one_eq_gram N ω
  · rw [hcoe]
    exact m.VtVG_one_eq_VtV N ω

end UnalignedModelR

/-! ### 7. `Ṽ_W V → W B_{β,R}` entrywise (D4 plan item 1.2)

The weighted twin of `VtV_tendsto_general` (`RankR/GeneralMain.lean`), and the exact mirror of
`gramRWG` above. It is the second entrywise limit that the component clause of
`thm:rank_r_svdstack` consumes (`RankR/AlignedMain.lean`). -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- `Ṽ_W V → W B_{β,R}` entrywise in probability. Each entry is a fixed real linear
combination of the entries of `Ṽ V` (`VtVWG_eq`), and `VtV_tendsto_general` gives those
limits. The same step runs inline in `thm_gen_rank_weight_svdstak_general_r_conv`
(`RankR/GeneralMain.lean`); it is named here for the component clause of
`thm:rank_r_svdstack`. As `VtV_tendsto_general` does, this needs no independence across
tables, so `hG : m.IndepNoise` does not appear. Mirror: `gramRWG` above. -/
theorem VtVRWG (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (p : Fin (rtot rk)) (k : Fin r) :
    TendstoInProb μ (fun N ω => m.VtVWG W N ω p k) (BBlockW W β m.R p k) := by
  have hconv0 : TendstoInProbPi μ
      (fun N ω => fun t : Fin (rtot rk) × Fin r => m.VtVG N ω t.1 t.2)
      (fun t : Fin (rtot rk) × Fin r => BBlock β m.R t.1 t.2) :=
    fun t => m.VtV_tendsto_general c β hβdef law t.1 t.2
  have hcont : Continuous fun z : Fin (rtot rk) × Fin r → ℝ =>
      ∑ a : Fin (rtot rk), W p a * z (a, k) :=
    continuous_finsetSum _ fun a _ => continuous_const.mul (continuous_apply _)
  have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
  have hlim : BBlockW W β m.R p k = ∑ a : Fin (rtot rk), W p a * BBlock β m.R a k :=
    Matrix.mul_apply
  rw [hlim]
  refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
  change ∑ a : Fin (rtot rk), W p a * m.VtVG N ω a k = m.VtVWG W N ω p k
  rw [m.VtVWG_eq]
  exact (Matrix.mul_apply).symm

end UnalignedModelR

end StackedSVD
