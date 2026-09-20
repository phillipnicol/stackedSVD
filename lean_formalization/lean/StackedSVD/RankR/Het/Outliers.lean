/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.OutliersG
import StackedSVD.RankR.RMT.EdgeTauR
import StackedSVD.RankR.RMT.AlignTauR
import StackedSVD.RankR.Het.Forms
import StackedSVD.RankR.Het.Edge
import StackedSVD.RankR.Het.Scalars
import StackedSVD.LinAlg.SpecWindow
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# Track E, task E7: the outliers of the heteroscedastic rank-`r` Gram matrix

Paper: `thm:rank_r_stacksvd` (the weighted rank-`r` stacked SVD, `main_paper.tex:2337`),
read through `thm:stacksvd_weighted` and `eq:assumption4`. This file lands the three
Gaussian-discharged facts about the eigenvalues of `X_W X_Wᵀ` above a threshold `τ` that sits
strictly above the bulk edge `MPhet.bHet c w`:

1. **The edge at `τ`** (`UnalignedModelR.tendsto_measure_lamMax_w1_het_le_tau`,
   `UnalignedModelR.tendsto_measure_eigenvalues₀_le_tau_het`). With the columns of `Q` whose
   outlier sits below `τ` (or which are subcritical) collected in `Fin t`, the block
   `W₀' + Q_sub Q_subᵀ` has its top eigenvalue at most `τ`, with probability tending to 1.
   The Track C mirror is `RankRStack.tendsto_measure_lamMax_w1R_le_tau`
   (`RankR/RMT/EdgeTauR.lean`); the rank-one mirror is
   `MultiTableModel.align_tendstoInProb_het_subcritical` (`RMT/Het/R6het.lean`).
2. **The count** (`UnalignedModelR.tendsto_measure_count_Ioi_tau_het`). The sorted eigenvalues
   above `τ` are exactly the indices below the number of components `l` with
   `eq:assumption4` and `τ < rhoHet θ_l c w`. The Track C mirror is
   `RankRStack.tendsto_measure_count_Ioi_tau` (`RankR/RMT/AlignTauR.lean`).
3. **The half-line overlap** (`UnalignedModelR.tendstoInProb_normSq_specProj_Ioi_tau_het`).
   The squared norm of the projection of the column `Q_l` on the eigenvalues above `τ` tends
   to `1 / nuHet θ_l c w` when `l` carries an outlier above `τ`, and to `0` otherwise. The
   rank-one mirror is `MultiTableModel.align_tendstoInProb_het` (`RMT/Het/R5het.lean:680`).

## The route

The deterministic core `OutliersR.align_detG` (`RankR/RMT/OutliersG.lean`) takes a unit test
vector `v`. The heteroscedastic column `Q_l` is not a unit vector: its squared norm tends to
`colGramLimit w c l = ∑ w_i² θ_il² + ∑ w_i² c_i` (`Forms.lean`, task E4). This file applies
the core to the scaled vector `v = (√(Q_l ⬝ Q_l))⁻¹ • Q_l` on the event where the column Gram
is close to its limit, and reads the answer back through
`OutliersR.normSq_specProj_eq_dot_mul`. The frame scalars `t_k` of the core are then
`(√N_l)⁻¹ * (if f k = l then -1 else 0)`, from the `-1` identity `1 + Fhet θ (ρ) = 0`
(`MPhet.one_add_F_rhoHet`).

The count filter uses `Classical` decidability, the same instance that `Scalars.ellSup`
(`RankR/Het/Scalars.lean`) uses, so the two cardinalities agree by `rfl`.

Numeric check: `$SP/agents/E7/check_frame.py`, seed 20260902; see
`notes/archive/agent_reports/d32_E7_hetoutliers.md`.

Every declaration is proved in full.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. Deterministic helpers of the scaled route -/

namespace OutliersR

/-- The scaled error budget: if `A` is within `B` of `T` and the scale `g` is within `η₂` of
`Nl`, then `g A` is within `(Nl + 1) B + η₂ |T|` of `Nl T`. -/
theorem abs_mul_sub_le_of_close {g A T B η₂ Nl : ℝ} (hg : 0 < g) (hgle : g ≤ Nl + 1)
    (hAT : |A - T| ≤ B) (hgN : |g - Nl| ≤ η₂) :
    |g * A - Nl * T| ≤ (Nl + 1) * B + η₂ * |T| := by
  have h1 : g * A - Nl * T = g * (A - T) + (g - Nl) * T := by ring
  rw [h1]
  refine le_trans (abs_add_le _ _) ?_
  rw [abs_mul, abs_mul, abs_of_pos hg]
  have h2 : g * |A - T| ≤ (Nl + 1) * B :=
    mul_le_mul hgle hAT (abs_nonneg _) (by linarith)
  have h3 : |g - Nl| * |T| ≤ η₂ * |T| := mul_le_mul_of_nonneg_right hgN (abs_nonneg _)
  linarith

/-- The scaled vector `(√(x ⬝ x))⁻¹ • x` is a unit vector. -/
theorem dot_self_normalize {p : ℕ} {x : Fin p → ℝ} (hx : 0 < x ⬝ᵥ x) :
    ((Real.sqrt (x ⬝ᵥ x))⁻¹ • x) ⬝ᵥ ((Real.sqrt (x ⬝ᵥ x))⁻¹ • x) = 1 := by
  rw [smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc, ← mul_inv,
    Real.mul_self_sqrt hx.le, inv_mul_cancel₀ hx.ne']

/-- The projection norm of `x` is the projection norm of its scaled version, times `x ⬝ x`. -/
theorem normSq_specProj_eq_dot_mul {p : ℕ} (S : Matrix (Fin p) (Fin p) ℝ) (T : Set ℝ)
    {x : Fin p → ℝ} (hx : 0 < x ⬝ᵥ x) :
    ‖specProj S T (WithLp.toLp 2 x)‖ ^ 2
      = (x ⬝ᵥ x) * ‖specProj S T (WithLp.toLp 2 ((Real.sqrt (x ⬝ᵥ x))⁻¹ • x))‖ ^ 2 := by
  rw [WithLp.toLp_smul, map_smul, norm_smul, Real.norm_eq_abs, mul_pow, sq_abs, inv_pow,
    Real.sq_sqrt hx.le, ← mul_assoc, mul_inv_cancel₀ hx.ne', one_mul]

end OutliersR

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 2. The column Gram limit as a named scalar -/

/-- `N_l := ∑ w_i² θ_il² + ∑ w_i² c_i`, the limit of `Q_l ⬝ᵥ Q_l`
(`UnalignedModelR.tendstoInProb_dotProduct_QmatHetR_col`, task E4). The rank-one mirror is the
limit of `‖q‖²` in `MultiTableModel.tendstoInProb_qHet_norm` (`RMT/Het/R5het.lean`). -/
noncomputable def colGramLimit (m : UnalignedModelR μ M n d r (alignedRk M r))
    (w c : Fin M → ℝ) (l : Fin r) : ℝ :=
  ∑ i, w i ^ 2 * m.thetaAligned i l ^ 2 + ∑ i, w i ^ 2 * c i

/-- `N₀ := ∑ w_i² c_i`, the `θ`-free floor of every column Gram limit. -/
theorem colGramFloor_le_colGramLimit (m : UnalignedModelR μ M n d r (alignedRk M r))
    (w c : Fin M → ℝ) (l : Fin r) :
    ∑ i, w i ^ 2 * c i ≤ m.colGramLimit w c l := by
  have h1 : 0 ≤ ∑ i, w i ^ 2 * m.thetaAligned i l ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (sq_nonneg _) (sq_nonneg _)
  unfold colGramLimit
  linarith

/-- `0 < N₀` from one nonzero weight and positive `c`. -/
theorem colGramFloor_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    0 < ∑ i, w i ^ 2 * c i := by
  obtain ⟨i₀, hi₀⟩ := hw
  have hw2 : 0 < w i₀ ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi₀))
  have h2 : 0 < w i₀ ^ 2 * c i₀ := mul_pos hw2 (hc i₀)
  have h3 : w i₀ ^ 2 * c i₀ ≤ ∑ i, w i ^ 2 * c i :=
    Finset.single_le_sum (f := fun i => w i ^ 2 * c i)
      (fun i _ => mul_nonneg (sq_nonneg _) (hc i).le) (Finset.mem_univ i₀)
  linarith

/-- `0 < N_l` as soon as one weight is nonzero and every `c_i` is positive: `N_l` is `N₀`
plus a nonnegative term. -/
theorem colGramLimit_pos (m : UnalignedModelR μ M n d r (alignedRk M r))
    (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (l : Fin r) :
    0 < m.colGramLimit w c l :=
  lt_of_lt_of_le (colGramFloor_pos hc hw) (m.colGramFloor_le_colGramLimit w c l)

/-- `Q_l ⬝ᵥ Q_l → N_l` in probability, the diagonal case of
`UnalignedModelR.tendstoInProb_dotProduct_QmatHetR_col`. -/
theorem tendstoInProb_colGram [NeZero M]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) (l : Fin r) :
    TendstoInProb μ (fun N ω =>
        (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      (m.colGramLimit w c l) :=
  FormsR.tendstoInProb_congr_limit (if_pos rfl)
    (m.tendstoInProb_dotProduct_QmatHetR_col w c hw hR hreg hG hpd l l)

/-! ### 3. The edge at `τ` -/

/-- The columns of `Q` selected by `f`. The Track C mirror is `RankRStack.qsubR`
(`RankR/RMT/EdgeR.lean`). -/
noncomputable def qsubHet (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {t : ℕ} (f : Fin t → Fin r) : Matrix (Fin (∑ i, n i N)) (Fin t) ℝ :=
  (m.QmatHetR w N ω).submatrix id f

theorem qsubHet_apply (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {t : ℕ} (f : Fin t → Fin r) (q : Fin (∑ i, n i N)) (a : Fin t) :
    m.qsubHet w N ω f q a = m.QmatHetR w N ω q (f a) := rfl

/-- `W₁ = W₀' + Q_sub Q_subᵀ`, the block whose edge closes the count. The Track C mirror is
`RankRStack.w1R` (`RankR/RMT/EdgeR.lean`). -/
noncomputable def w1Het (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {t : ℕ} (f : Fin t → Fin r) :
    Matrix (Fin (∑ i, n i N)) (Fin (∑ i, n i N)) ℝ :=
  m.W0hetR w N ω + m.qsubHet w N ω f * (m.qsubHet w N ω f)ᵀ

theorem isHermitian_w1Het (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {t : ℕ} (f : Fin t → Fin r) : (m.w1Het w N ω f).IsHermitian :=
  (m.isHermitian_W0hetR w N ω).add (EdgeR.isHermitian_mul_transpose _)

/-- **The edge of the block at `τ`.** With every selected column `f a` either subcritical or
with its outlier `rhoHet` strictly below `τ`, `lamMax (W₀' + Q_sub Q_subᵀ) ≤ τ` with
probability tending to 1. The Track C mirror is `RankRStack.tendsto_measure_lamMax_w1R_le_tau`
(`RankR/RMT/EdgeTauR.lean`); the diagonal limit `1 + Fhet θ_{f a} c w τ` is positive by
`MPhet.one_add_F_pos` in the supercritical case and `MPhet.one_add_F_pos_of_not_assumption4`
otherwise. The rank-one mirror is `MultiTableModel.align_tendstoInProb_het_subcritical`
(`RMT/Het/R6het.lean`). -/
theorem tendsto_measure_lamMax_w1_het_le_tau [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {t : ℕ} {f : Fin t → Fin r} (hf : Function.Injective f)
    {τ : ℝ} (hτ : MPhet.bHet c w < τ)
    (hsub : ∀ a : Fin t, Scalars.Assumption4 (fun i => m.thetaAligned i (f a)) c w →
      MPhet.rhoHet (fun i => m.thetaAligned i (f a)) c w < τ) :
    Tendsto (fun N => μ N
        {ω | lamMax (m.w1Het w N ω f) (m.isHermitian_w1Het w N ω f) ≤ τ}) atTop (𝓝 1) := by
  classical
  have H := m.resolventLimitsHetR_of_gaussian w c hc hw hR hreg hG hpd le_rfl hedge.edge
  -- the diagonal limit is positive at every index, in both regimes
  set L : Fin t → ℝ := fun a => 1 + MPhet.Fhet (fun i => m.thetaAligned i (f a)) c w τ
    with hLdef
  have hLpos : ∀ a, 0 < L a := by
    intro a
    by_cases h4 : Scalars.Assumption4 (fun i => m.thetaAligned i (f a)) c w
    · exact MPhet.one_add_F_pos hc h4 (hsub a h4)
    · exact MPhet.one_add_F_pos_of_not_assumption4 hc hw h4 hτ
  obtain ⟨L₀, hL₀pos, hL₀le⟩ := ScalarsC.exists_pos_lower_bound L hLpos
  -- the single accuracy
  have ht1 : (0 : ℝ) < (t : ℝ) + 1 := by positivity
  set δ : ℝ := L₀ / ((t : ℝ) + 1) with hδdef
  have hδpos : 0 < δ := div_pos hL₀pos ht1
  have hδt : ∀ a, δ * (t : ℝ) ≤ L a := by
    intro a
    refine le_trans ?_ (hL₀le a)
    rw [hδdef, div_mul_eq_mul_div, div_le_iff₀ ht1]
    nlinarith [hL₀pos.le, Nat.cast_nonneg (α := ℝ) t]
  -- the edge event, at half the room between the edge and `τ`
  have hε2 : (0 : ℝ) < (τ - MPhet.bHet c w) / 2 := by linarith
  have hedgeC : Tendsto (fun N => μ N
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
        ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + (τ - MPhet.bHet c w) / 2)).nullMeasurableSet)
      (H.edge ((τ - MPhet.bHet c w) / 2) hε2)
  have hET : ∀ q : Fin t × Fin t, Tendsto (fun N => μ N
      {ω | δ ≤ |R4.cform (m.W0hetR w N ω) τ
        (fun l => m.QmatHetR w N ω l (f q.1)) (fun l => m.QmatHetR w N ω l (f q.2))
        - (if f q.1 = f q.2 then MPhet.Phihet (fun i => m.thetaAligned i (f q.1)) c w τ
            + MPhet.Psihet c w τ else 0)|}) atTop (𝓝 0) := fun q =>
    H.cform_qcol (f q.1) (f q.2) hτ δ hδpos
  have hzero : Tendsto (fun N => μ N
      ({ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
            ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (m.W0hetR w N ω) τ
            (fun l => m.QmatHetR w N ω l (f q.1)) (fun l => m.QmatHetR w N ω l (f q.2))
            - (if f q.1 = f q.2 then MPhet.Phihet (fun i => m.thetaAligned i (f q.1)) c w τ
                + MPhet.Psihet c w τ else 0)|}))) atTop (𝓝 0) :=
    tendsto_measure_zero_union hedgeC (tendsto_measure_zero_iUnion hET)
  have hincl : ∀ N, ({ω | lamMax (m.w1Het w N ω f) (m.isHermitian_w1Het w N ω f)
        ≤ τ} : Set (Ω N))ᶜ
      ⊆ {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
            ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (m.W0hetR w N ω) τ
            (fun l => m.QmatHetR w N ω l (f q.1)) (fun l => m.QmatHetR w N ω l (f q.2))
            - (if f q.1 = f q.2 then MPhet.Phihet (fun i => m.thetaAligned i (f q.1)) c w τ
                + MPhet.Psihet c w τ else 0)|}) := by
    intro N ω hω
    by_contra hbad
    have hedgeN : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
        ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2 := by
      by_contra hxx
      exact hbad (Set.mem_union_left _ hxx)
    have hclose : ∀ q : Fin t × Fin t,
        |R4.cform (m.W0hetR w N ω) τ
          (fun l => m.QmatHetR w N ω l (f q.1)) (fun l => m.QmatHetR w N ω l (f q.2))
          - (if f q.1 = f q.2 then MPhet.Phihet (fun i => m.thetaAligned i (f q.1)) c w τ
              + MPhet.Psihet c w τ else 0)| ≤ δ := by
      intro q
      by_contra hxx
      exact hbad (Set.mem_union_right _ (Set.mem_iUnion.mpr ⟨q, (not_le.mp hxx).le⟩))
    have hzlam : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) < τ := by
      linarith
    have hM : ((1 : Matrix (Fin t) (Fin t) ℝ)
        + (m.qsubHet w N ω f)ᵀ * R4.resolv (m.W0hetR w N ω) τ
          * m.qsubHet w N ω f).PosSemidef := by
      refine EdgeR.posSemidef_of_close_to_diag ?_ hδpos.le hδt ?_
      · rw [Matrix.transpose_add, Matrix.transpose_one, Matrix.transpose_mul,
          Matrix.transpose_mul, Matrix.transpose_transpose,
          R4.transpose_resolv (m.isHermitian_W0hetR w N ω), Matrix.mul_assoc]
      · intro a b
        have hentry : ((1 : Matrix (Fin t) (Fin t) ℝ)
              + (m.qsubHet w N ω f)ᵀ * R4.resolv (m.W0hetR w N ω) τ
                * m.qsubHet w N ω f) a b
            - (if a = b then L a else 0)
            = R4.cform (m.W0hetR w N ω) τ
                (fun l => m.QmatHetR w N ω l (f a)) (fun l => m.QmatHetR w N ω l (f b))
              - (if f a = f b then MPhet.Phihet (fun i => m.thetaAligned i (f a)) c w τ
                  + MPhet.Psihet c w τ else 0) := by
          have hcol : ∀ jj : Fin t, (fun k => m.qsubHet w N ω f k jj)
              = fun l => m.QmatHetR w N ω l (f jj) := fun _ => rfl
          rw [Matrix.add_apply, EdgeR.cform_eq_entry, hcol a, hcol b, Matrix.one_apply, hLdef]
          rcases eq_or_ne a b with rfl | hab
          · rw [if_pos rfl, if_pos rfl, if_pos rfl]
            simp only [MPhet.Fhet]
            ring
          · rw [if_neg hab, if_neg hab, if_neg (fun hcon => hab (hf hcon))]
            ring
        rw [hentry]
        exact hclose (a, b)
    exact hω (EdgeR.lamMax_add_le_of_posSemidef (m.stack_row_pos N)
      (m.isHermitian_W0hetR w N ω) hzlam hM (m.isHermitian_w1Het w N ω f))
  exact tendsto_measure_one_of_bad hincl hzero

/-! ### 4. Measurability of the sorted Gram eigenvalues -/

/-- The sorted eigenvalue `k` of the Gram matrix `X_W X_Wᵀ` is measurable in `ω`. The Track C
mirror is `RankRStack.measurable_gram_eigenvalues₀` (`RankR/RMT/EdgeR.lean`); the Gram
matrix here is `X_W X_Wᵀ`, which is `gramEig` of `X_Wᵀ`. -/
theorem measurable_gramHet_eigenvalues₀ (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (k : Fin (Fintype.card (Fin (∑ i, n i N)))) :
    Measurable fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k := by
  have hT : Measurable fun ω => (m.stackXW w N ω)ᵀ :=
    measurable_pi_lambda _ fun a => measurable_pi_lambda _ fun b =>
      (measurable_pi_apply a).comp ((measurable_pi_apply b).comp (m.measurable_stackXW w N))
  have h : (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k)
      = fun ω => gramEig (m.stackXW w N ω)ᵀ (k : ℕ) := by
    funext ω
    rw [gramEig_of_lt _ k.isLt]
    rfl
  rw [h]
  exact (measurable_gramEig (k : ℕ)).comp hT

/-- The edge event of the count is measurable: a finite intersection of level sets of the
sorted eigenvalues. The Track C mirror is `RankRStack.measurableSet_eigenvalues₀_le`. -/
theorem measurableSet_eigenvalues₀_le_het (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (u : ℕ) (x : ℝ) :
    MeasurableSet {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
      (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ x} := by
  classical
  have hset : {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ x}
      = ⋂ k : Fin (Fintype.card (Fin (∑ i, n i N))),
          {ω : Ω N | u ≤ (k : ℕ) →
            (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ x} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iInter]
  rw [hset]
  refine MeasurableSet.iInter fun k => ?_
  by_cases hk : u ≤ (k : ℕ)
  · have he : {ω : Ω N | u ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ x}
        = (fun ω => (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k)
            ⁻¹' Set.Iic x := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic]
      exact ⟨fun h => h hk, fun h _ => h⟩
    rw [he]
    exact m.measurable_gramHet_eigenvalues₀ w N k measurableSet_Iic
  · have he : {ω : Ω N | u ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ x}
        = Set.univ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_univ, iff_true]
      exact fun hcon => absurd hcon hk
    rw [he]
    exact MeasurableSet.univ

/-! ### 5. The count -/

/-- **The sorted eigenvalues at index `u` and above are at most `τ`.** With the `r` components
split by `e` into `Fin t` (subcritical, or outlier below `τ`) and `Fin u` others, every sorted
eigenvalue of `X_W X_Wᵀ` at index `u` or above is at most `τ`, with probability tending to 1.
The Track C mirror is `RankRStack.tendsto_measure_eigenvalues₀_le_tau`
(`RankR/RMT/EdgeTauR.lean`). -/
theorem tendsto_measure_eigenvalues₀_le_tau_het [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {t u : ℕ} (e : Fin t ⊕ Fin u ≃ Fin r) {τ : ℝ} (hτ : MPhet.bHet c w < τ)
    (hsub : ∀ a : Fin t,
      Scalars.Assumption4 (fun i => m.thetaAligned i (e (Sum.inl a))) c w →
        MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w < τ) :
    Tendsto (fun N => μ N
        {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ})
      atTop (𝓝 1) := by
  classical
  have hf : Function.Injective (fun a : Fin t => e (Sum.inl a)) :=
    e.injective.comp Sum.inl_injective
  have hmain := m.tendsto_measure_lamMax_w1_het_le_tau w c hc hw hR hreg hG hpd hedge hf hτ hsub
  have hsubset : ∀ N, {ω | lamMax (m.w1Het w N ω (fun a : Fin t => e (Sum.inl a)))
        (m.isHermitian_w1Het w N ω (fun a : Fin t => e (Sum.inl a))) ≤ τ}
      ⊆ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
          (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ} := by
    intro N ω hω k hk
    have hSeq : m.stackXW w N ω * (m.stackXW w N ω)ᵀ
        = m.w1Het w N ω (fun a : Fin t => e (Sum.inl a))
          + m.qsubHet w N ω (fun b : Fin u => e (Sum.inr b))
            * (m.qsubHet w N ω (fun b : Fin u => e (Sum.inr b)))ᵀ := by
      rw [m.gram_eq_hetR w N ω, EdgeR.mul_transpose_split (m.QmatHetR w N ω) e, w1Het, qsubHet,
        qsubHet, add_assoc]
    exact Frame.eigenvalues₀_le_of_split
      (m.isHermitian_w1Het w N ω (fun a : Fin t => e (Sum.inl a)))
      (isHermitian_mul_transpose_self _) hSeq hω k hk
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le hmain tendsto_const_nhds
    (fun N => measure_mono (hsubset N)) (fun N => prob_le_one)

open Classical in
/-- **The count of the eigenvalues above `τ`.** With probability tending to 1 the sorted
eigenvalues of `X_W X_Wᵀ` above `τ` are exactly the indices below the number of components
`l` with `eq:assumption4` and `τ < rhoHet θ_l c w`. The separation hypothesis `hsep` keeps
every supercritical outlier away from `τ` by the margin `mg`; the junk value `rhoHet = 0` of a
subcritical component never enters, because `hsep` is read only under `eq:assumption4`.
The Track C mirror is `RankRStack.tendsto_measure_count_Ioi_tau` (`RankR/RMT/AlignTauR.lean`);
the filter uses the `Classical` instance of `Scalars.ellSup`. -/
theorem tendsto_measure_count_Ioi_tau_het [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : MPhet.bHet c w < τ)
    (hsep : ∀ k : Fin r, Scalars.Assumption4 (fun i => m.thetaAligned i k) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i k) c w < τ
        ∨ τ + mg ≤ MPhet.rhoHet (fun i => m.thetaAligned i k) c w)) :
    Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
        ↔ (k : ℕ) < (Finset.univ.filter fun l : Fin r =>
            Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
              ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w).card)}) atTop (𝓝 1) := by
  have H := m.resolventLimitsHetR_of_gaussian w c hc hw hR hreg hG hpd le_rfl hedge.edge
  obtain ⟨t, u, e, hsubE, hsupE, hcardu⟩ :=
    OutliersR.exists_sum_equiv_split_pred (fun l : Fin r =>
      Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
        ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w)
  set f : Fin u → Fin r := fun b => e (Sum.inr b) with hfdef
  have hfinj : Function.Injective f := fun a b hab => Sum.inr_injective (e.injective hab)
  have hsupE' : ∀ k : Fin u, Scalars.Assumption4 (fun i => m.thetaAligned i (f k)) c w
      ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i (f k)) c w := fun k => hsupE k
  have hsubE' : ∀ a : Fin t,
      Scalars.Assumption4 (fun i => m.thetaAligned i (e (Sum.inl a))) c w →
        MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w < τ := by
    intro a h4
    have hnot : ¬ (Scalars.Assumption4 (fun i => m.thetaAligned i (e (Sum.inl a))) c w
        ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w) := hsubE a
    have h1 : MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w ≤ τ :=
      not_lt.mp fun hlt => hnot ⟨h4, hlt⟩
    rcases hsep (e (Sum.inl a)) h4 with hlt | hge
    · exact hlt
    · linarith
  -- 1. the frame scalars
  set ρv : Fin u → ℝ := fun k => MPhet.rhoHet (fun i => m.thetaAligned i (f k)) c w with hρv
  set νv : Fin u → ℝ := fun k => MPhet.nuHet (fun i => m.thetaAligned i (f k)) c w with hνv
  have hbρ : ∀ k, MPhet.bHet c w < ρv k := fun k => MPhet.bHet_lt_rhoHet hc (hsupE' k).1
  have hm1 : ∀ k, MPhet.Phihet (fun i => m.thetaAligned i (f k)) c w (ρv k)
      + MPhet.Psihet c w (ρv k) = -1 := by
    intro k
    have h1 := MPhet.one_add_F_rhoHet hc (hsupE' k).1
    simp only [MPhet.Fhet] at h1
    linarith
  have hνpos : ∀ k, 0 < νv k := fun k => MPhet.nuHet_pos hc (hsupE' k).1
  have hτρ : ∀ k, τ + mg ≤ ρv k := by
    intro k
    rcases hsep (f k) (hsupE' k).1 with hlt | hge
    · exact absurd (hsupE' k).2 (not_lt.mpr hlt.le)
    · exact hge
  -- 2. the column-norm constant, read from the column Gram limit
  set Cq : ℝ := ∑ l, (m.colGramLimit w c l + 1) with hCqdef
  have hNl0 : ∀ l, 0 ≤ m.colGramLimit w c l + 1 := fun l =>
    (m.colGramLimit_pos w c hc hw l).le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hNl0 l
  have hCqle : ∀ l, m.colGramLimit w c l + 1 ≤ Cq := fun l =>
    Finset.single_le_sum (f := fun l => m.colGramLimit w c l + 1) (fun l' _ => hNl0 l')
      (Finset.mem_univ l)
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
  -- 4. the five bad families
  have hε2 : (0 : ℝ) < (τ - MPhet.bHet c w) / 2 := by linarith
  have hedgeC1 : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + (τ - MPhet.bHet c w) / 2)).nullMeasurableSet)
      (H.edge ((τ - MPhet.bHet c w) / 2) hε2)
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimit w c l|}) atTop (𝓝 0) :=
    fun l => m.tendstoInProb_colGram w c hw hR hreg hG hpd l 1 one_pos
  have hE1T : ∀ q : Fin r × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2))
        - (if q.1 = f q.2 then -1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := H.cform_qcol q.1 (f q.2) (hbρ q.2)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2)))
        (if q.1 = f q.2 then -1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 (f q.2) with heq | hne'
      · rw [if_pos heq, if_pos heq, heq]
        exact hm1 q.2
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  have hE2T : ∀ q : Fin u × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv q.1)
        (fun i => m.QmatHetR w N ω i (f q.1)) (fun i => m.QmatHetR w N ω i (f q.2))
        - (if q.1 = q.2 then νv q.1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := H.cform2_qcol (f q.1) (f q.2) (hbρ q.1)
    have h2 : TendstoInProb μ (fun N ω => R4.cform2 (m.W0hetR w N ω) (ρv q.1)
        (fun i => m.QmatHetR w N ω i (f q.1)) (fun i => m.QmatHetR w N ω i (f q.2)))
        (if q.1 = q.2 then νv q.1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 q.2 with heq | hne'
      · rw [if_pos heq, if_pos (congrArg f heq)]
        rfl
      · rw [if_neg hne', if_neg (fun hcon => hne' (hfinj hcon))]
    exact h2 η hη0
  have hedgeC7 : Tendsto (fun N => μ N
      {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_eigenvalues₀_le_het w N u τ).nullMeasurableSet)
      (m.tendsto_measure_eigenvalues₀_le_tau_het w c hc hw hR hreg hG hpd hedge e hτ hsubE')
  -- 5. assemble
  refine tendsto_measure_one_of_bad (s := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ
        ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l)
                ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimit w c l|})
          ∪ ((⋃ q : Fin r × Fin u, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
                (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2))
                - (if q.1 = f q.2 then -1 else 0)|})
            ∪ ((⋃ q : Fin u × Fin u, {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv q.1)
                  (fun i => m.QmatHetR w N ω i (f q.1))
                  (fun i => m.QmatHetR w N ω i (f q.2))
                  - (if q.1 = q.2 then νv q.1 else 0)|})
              ∪ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
                  (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
                    ≤ τ}ᶜ)))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hg1, hg2, hg4, hg5, hg7⟩ := hbad
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by
      linarith
    have hcolb : ∀ l : Fin r, ∑ i, m.QmatHetR w N ω i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      have hsq : ∑ i, m.QmatHetR w N ω i l ^ 2
          = (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) :=
        Finset.sum_congr rfl fun i _ => sq _
      rw [hsq]
      have h1 := (abs_lt.mp (hg2 l)).2
      linarith
    have hI := OutliersR.count_eq_of_forms_of_edge (Q := m.QmatHetR w N ω)
      (τ := τ) (mg := mg) (Cq := Cq) (η := η)
      (m.isHermitian_W0hetR w N ω) (isHermitian_mul_transpose_self (m.stackXW w N ω))
      (m.gram_eq_hetR w N ω) hfinj hmg
      hlamτ hτρ hνpos hCq0 hη0 hcolb (fun l k => (hg4 (l, k)).le) (fun k l => (hg5 (k, l)).le)
      hsmall hδm hg7
    apply hω
    intro k
    rw [← hcardu]
    exact Frame.lt_eigenvalues₀_iff_of_card_eigenvalues
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) hI k
  · exact tendsto_measure_zero_union hedgeC1
      (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
            hedgeC7)))

/-! ### 6. The half-line overlap of a column of `Q` -/

open Classical in
/-- **The half-line overlap of a column of `Q`.** The squared norm of the projection of the
column `Q_l` on the eigenvalues of `X_W X_Wᵀ` above `τ` tends in probability to
`1 / nuHet θ_l c w` when `l` satisfies `eq:assumption4` with `τ < rhoHet θ_l c w`, and to `0`
otherwise (a subcritical `l`, or an outlier below `τ`). The `1 / ν` is the paper's
`1 / (ρ F'(ρ))` before the factor `ρ`, the overlap of `thm:rank_r_stacksvd` read through
`main_paper.tex:1433`; the column `Q_l` is not normalized, and its squared norm tends to
`N_l = colGramLimit w c l`, which cancels inside this proof against the normalization of the
unit vector. The consumer of the statement divides by the eigenvalue `λ_{ellSup}`, not by
`N_l`.

The route: the deterministic core `OutliersR.align_detG` is applied to the unit vector
`(√(Q_l ⬝ Q_l))⁻¹ • Q_l` on the event where `Q_l ⬝ Q_l` is close to `colGramLimit w c l`;
`OutliersR.normSq_specProj_eq_dot_mul` reads the answer back on `Q_l`. The frame scalars are
`(√N_l)⁻¹ * (if f k = l then -1 else 0)` by `MPhet.one_add_F_rhoHet`. The Track C mirror is
`RankRStack.tendstoInProb_normSq_specProj_Ioi_tau` (`RankR/RMT/AlignTauR.lean`); the rank-one
mirror is `MultiTableModel.align_tendstoInProb_het` (`RMT/Het/R5het.lean:680`). The junk value
`rhoHet = 0` of a subcritical component never enters: `hsep` and the limit read `rhoHet` only
under `eq:assumption4`. -/
theorem tendstoInProb_normSq_specProj_Ioi_tau_het [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w)) (l : Fin r)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : MPhet.bHet c w < τ)
    (hsep : ∀ k : Fin r, Scalars.Assumption4 (fun i => m.thetaAligned i k) c w →
      (MPhet.rhoHet (fun i => m.thetaAligned i k) c w < τ
        ∨ τ + mg ≤ MPhet.rhoHet (fun i => m.thetaAligned i k) c w)) :
    TendstoInProb μ
      (fun N ω => ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2)
      (if Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
          ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w
        then 1 / MPhet.nuHet (fun i => m.thetaAligned i l) c w else 0) := by
  have H := m.resolventLimitsHetR_of_gaussian w c hc hw hR hreg hG hpd le_rfl hedge.edge
  obtain ⟨t, u, e, hsubE, hsupE, hcardu⟩ :=
    OutliersR.exists_sum_equiv_split_pred (fun l' : Fin r =>
      Scalars.Assumption4 (fun i => m.thetaAligned i l') c w
        ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l') c w)
  set f : Fin u → Fin r := fun b => e (Sum.inr b) with hfdef
  have hfinj : Function.Injective f := fun a b hab => Sum.inr_injective (e.injective hab)
  have hsupE' : ∀ k : Fin u, Scalars.Assumption4 (fun i => m.thetaAligned i (f k)) c w
      ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i (f k)) c w := fun k => hsupE k
  have hsubE' : ∀ a : Fin t,
      Scalars.Assumption4 (fun i => m.thetaAligned i (e (Sum.inl a))) c w →
        MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w < τ := by
    intro a h4
    have hnot : ¬ (Scalars.Assumption4 (fun i => m.thetaAligned i (e (Sum.inl a))) c w
        ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w) := hsubE a
    have h1 : MPhet.rhoHet (fun i => m.thetaAligned i (e (Sum.inl a))) c w ≤ τ :=
      not_lt.mp fun hlt => hnot ⟨h4, hlt⟩
    rcases hsep (e (Sum.inl a)) h4 with hlt | hge
    · exact hlt
    · linarith
  -- 1. the frame scalars
  set ρv : Fin u → ℝ := fun k => MPhet.rhoHet (fun i => m.thetaAligned i (f k)) c w with hρv
  set νv : Fin u → ℝ := fun k => MPhet.nuHet (fun i => m.thetaAligned i (f k)) c w with hνv
  have hbρ : ∀ k, MPhet.bHet c w < ρv k := fun k => MPhet.bHet_lt_rhoHet hc (hsupE' k).1
  have hm1 : ∀ k, MPhet.Phihet (fun i => m.thetaAligned i (f k)) c w (ρv k)
      + MPhet.Psihet c w (ρv k) = -1 := by
    intro k
    have h1 := MPhet.one_add_F_rhoHet hc (hsupE' k).1
    simp only [MPhet.Fhet] at h1
    linarith
  have hνpos : ∀ k, 0 < νv k := fun k => MPhet.nuHet_pos hc (hsupE' k).1
  have hτρ : ∀ k, τ + mg ≤ ρv k := by
    intro k
    rcases hsep (f k) (hsupE' k).1 with hlt | hge
    · exact absurd (hsupE' k).2 (not_lt.mpr hlt.le)
    · exact hge
  -- 2. the column-norm constant and the scale of the column `l`
  set Cq : ℝ := ∑ l', (m.colGramLimit w c l' + 1) with hCqdef
  have hNl0 : ∀ l', 0 ≤ m.colGramLimit w c l' + 1 := fun l' =>
    (m.colGramLimit_pos w c hc hw l').le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l' _ => hNl0 l'
  have hCqle : ∀ l', m.colGramLimit w c l' + 1 ≤ Cq := fun l' =>
    Finset.single_le_sum (f := fun l' => m.colGramLimit w c l' + 1) (fun l'' _ => hNl0 l'')
      (Finset.mem_univ l')
  have hNl : 0 < m.colGramLimit w c l := m.colGramLimit_pos w c hc hw l
  set sN : ℝ := Real.sqrt (m.colGramLimit w c l) with hsNdef
  have hsN : 0 < sN := Real.sqrt_pos.mpr hNl
  have hsN2 : sN ^ 2 = m.colGramLimit w c l := Real.sq_sqrt hNl.le
  have hsNi : 0 < sN⁻¹ := inv_pos.mpr hsN
  -- 3. the frame scalars of the scaled column, and the scalar identity
  set tv : Fin u → ℝ := fun k => sN⁻¹ * (if f k = l then -1 else 0) with htv
  have hLb : ∀ k, |tv k| ≤ sN⁻¹ := by
    intro k
    have h1 : tv k = sN⁻¹ * (if f k = l then -1 else 0) := rfl
    rw [h1, abs_mul, abs_of_pos hsNi]
    refine mul_le_of_le_one_right hsNi.le ?_
    rcases eq_or_ne (f k) l with hkl | hkl
    · rw [if_pos hkl, abs_neg, abs_one]
    · rw [if_neg hkl, abs_zero]
      exact zero_le_one
  have htarget : m.colGramLimit w c l * ∑ k, tv k ^ 2 / νv k
      = if Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
          ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w
        then 1 / MPhet.nuHet (fun i => m.thetaAligned i l) c w else 0 := by
    by_cases hj : Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
        ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w
    · rw [if_pos hj]
      obtain ⟨k₀, hk₀⟩ : ∃ k₀ : Fin u, f k₀ = l := by
        rcases hsym : e.symm l with a | b
        · exfalso
          have hl1 : e (Sum.inl a) = l := by rw [← hsym]; exact e.apply_symm_apply l
          have h1 := hsubE a
          rw [hl1] at h1
          exact h1 hj
        · refine ⟨b, ?_⟩
          have hl1 : e (Sum.inr b) = l := by rw [← hsym]; exact e.apply_symm_apply l
          exact hl1
      have hcongr : ∀ k : Fin u, tv k ^ 2 / νv k
          = if k = k₀ then sN⁻¹ ^ 2 / νv k₀ else 0 := by
        intro k
        have h1 : tv k = sN⁻¹ * (if f k = l then -1 else 0) := rfl
        rcases eq_or_ne k k₀ with rfl | hne
        · rw [h1, if_pos hk₀, if_pos rfl]
          ring
        · have hfk : f k ≠ l := fun hcon => hne (hfinj (hcon.trans hk₀.symm))
          rw [h1, if_neg hfk, if_neg hne]
          simp
      rw [Finset.sum_congr rfl fun k _ => hcongr k, Finset.sum_ite_eq']
      have hνk₀ : νv k₀ = MPhet.nuHet (fun i => m.thetaAligned i l) c w := by
        have h1 : νv k₀ = MPhet.nuHet (fun i => m.thetaAligned i (f k₀)) c w := rfl
        rw [h1, hk₀]
      simp only [Finset.mem_univ, if_true]
      rw [hνk₀, inv_pow, hsN2, ← mul_div_assoc, mul_inv_cancel₀ hNl.ne']
    · rw [if_neg hj]
      have hno : ∀ k : Fin u, f k ≠ l := by
        intro k hcon
        have h1 := hsupE' k
        rw [hcon] at h1
        exact hj h1
      have hzero : ∀ k : Fin u, tv k ^ 2 / νv k = 0 := by
        intro k
        have h1 : tv k = sN⁻¹ * (if f k = l then -1 else 0) := rfl
        rw [h1, if_neg (hno k)]
        simp
      rw [Finset.sum_congr rfl fun k _ => hzero k]
      simp
  intro ε hε
  -- 4. the constants and the two accuracies `η`, `η₂`
  have hνsum0 : (0 : ℝ) ≤ ∑ k, (νv k)⁻¹ :=
    Finset.sum_nonneg fun k _ => inv_nonneg.mpr (hνpos k).le
  have hL0 : (0 : ℝ) ≤ 1 + 2 * sN⁻¹ := by linarith
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv := OutliersR.resCG_nonneg Cq r νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  have hX0 : (0 : ℝ) ≤ (u : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resCG Cq r νv / mg) :=
    mul_nonneg (Nat.cast_nonneg u) (add_nonneg hCG0 (div_nonneg hCR0 hmg.le))
  have hB0 : (0 : ℝ) ≤ 5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg)) + (1 + 2 * sN⁻¹) * ∑ k, (νv k)⁻¹ := by
    have h1 := mul_nonneg hL0 hνsum0
    linarith
  set K : ℝ := (m.colGramLimit w c l + 1) * (5 * ((u : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resCG Cq r νv / mg)) + (1 + 2 * sN⁻¹) * ∑ k, (νv k)⁻¹) with hKdef
  have hK0 : 0 ≤ K := mul_nonneg (hNl0 l) hB0
  have hε2 : (0 : ℝ) < ε / 2 := by linarith
  set η : ℝ := min (min 1 (mg / (OutliersR.resCG Cq r νv + 1)))
      (min (1 / (2 * ((u : ℝ) * (OutliersR.gramC ρv νv
          + OutliersR.resCG Cq r νv / mg) + 1)))
        (ε / 2 / (2 * (K + 1)))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    refine lt_min (lt_min one_pos (div_pos hmg (by linarith))) (lt_min ?_ ?_)
    · exact div_pos one_pos (by linarith)
    · exact div_pos hε2 (by linarith)
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
  have hηε : η ≤ ε / 2 / (2 * (K + 1)) := by
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
  have hfinal : (m.colGramLimit w c l + 1) * (5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
        + OutliersR.resCG Cq r νv * η / mg)
      + (1 + 2 * sN⁻¹) * (∑ k, (νv k)⁻¹) * η) ≤ ε / 2 / 2 := by
    have hKe : (m.colGramLimit w c l + 1) * (5 * (u : ℝ) * (OutliersR.gramC ρv νv * η
          + OutliersR.resCG Cq r νv * η / mg)
        + (1 + 2 * sN⁻¹) * (∑ k, (νv k)⁻¹) * η) = K * η := by
      rw [hKdef]
      ring
    rw [hKe]
    exact le_trans (mul_le_mul_of_nonneg_left hηε hK0) (OutliersR.mul_div_le_half hK0 hε2.le)
  set η₂ : ℝ := min (m.colGramLimit w c l / 2)
    (ε / 2 / (2 * (|∑ k, tv k ^ 2 / νv k| + 1))) with hη₂def
  have hη₂0 : 0 < η₂ := by
    rw [hη₂def]
    refine lt_min (by linarith) (div_pos hε2 ?_)
    have h1 := abs_nonneg (∑ k, tv k ^ 2 / νv k)
    linarith
  have hη₂half : η₂ ≤ m.colGramLimit w c l / 2 := by
    rw [hη₂def]
    exact min_le_left _ _
  have hη₂T : η₂ ≤ ε / 2 / (2 * (|∑ k, tv k ^ 2 / νv k| + 1)) := by
    rw [hη₂def]
    exact min_le_right _ _
  have hfinal2 : η₂ * |∑ k, tv k ^ 2 / νv k| ≤ ε / 2 / 2 := by
    rw [mul_comm]
    exact le_trans (mul_le_mul_of_nonneg_left hη₂T (abs_nonneg _))
      (OutliersR.mul_div_le_half (abs_nonneg _) hε2.le)
  -- 5. the seven bad families
  have hε3 : (0 : ℝ) < (τ - MPhet.bHet c w) / 2 := by linarith
  have hedgeC1 : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + (τ - MPhet.bHet c w) / 2)).nullMeasurableSet)
      (H.edge ((τ - MPhet.bHet c w) / 2) hε3)
  have hcolT : ∀ l' : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l') ⬝ᵥ (fun q => m.QmatHetR w N ω q l')
        - m.colGramLimit w c l'|}) atTop (𝓝 0) :=
    fun l' => m.tendstoInProb_colGram w c hw hR hreg hG hpd l' 1 one_pos
  have hcol2T : Tendsto (fun N => μ N
      {ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimit w c l|}) atTop (𝓝 0) :=
    m.tendstoInProb_colGram w c hw hR hreg hG hpd l η₂ hη₂0
  have hE1T : ∀ q : Fin r × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2))
        - (if q.1 = f q.2 then -1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := H.cform_qcol q.1 (f q.2) (hbρ q.2)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2)))
        (if q.1 = f q.2 then -1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 (f q.2) with heq | hne'
      · rw [if_pos heq, if_pos heq, heq]
        exact hm1 q.2
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  have hE2T : ∀ q : Fin u × Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv q.1)
        (fun i => m.QmatHetR w N ω i (f q.1)) (fun i => m.QmatHetR w N ω i (f q.2))
        - (if q.1 = q.2 then νv q.1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := H.cform2_qcol (f q.1) (f q.2) (hbρ q.1)
    have h2 : TendstoInProb μ (fun N ω => R4.cform2 (m.W0hetR w N ω) (ρv q.1)
        (fun i => m.QmatHetR w N ω i (f q.1)) (fun i => m.QmatHetR w N ω i (f q.2)))
        (if q.1 = q.2 then νv q.1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 q.2 with heq | hne'
      · rw [if_pos heq, if_pos (congrArg f heq)]
        rfl
      · rw [if_neg hne', if_neg (fun hcon => hne' (hfinj hcon))]
    exact h2 η hη0
  have hE3T : ∀ k : Fin u, Tendsto (fun N => μ N
      {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q l)
            ⬝ᵥ (fun q => m.QmatHetR w N ω q l)))⁻¹
          * R4.cform (m.W0hetR w N ω) (ρv k)
            (fun q => m.QmatHetR w N ω q l) (fun q => m.QmatHetR w N ω q (f k))
        - tv k|}) atTop (𝓝 0) := by
    intro k
    have h1 := H.cform_qcol l (f k) (hbρ k)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) (ρv k)
        (fun q => m.QmatHetR w N ω q l) (fun q => m.QmatHetR w N ω q (f k)))
        (if f k = l then -1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne (f k) l with hkl | hkl
      · rw [if_pos hkl.symm, if_pos hkl, ← hkl]
        exact hm1 k
      · rw [if_neg (Ne.symm hkl), if_neg hkl]
    have h3 : TendstoInProb μ (fun N ω => (Real.sqrt ((fun q => m.QmatHetR w N ω q l)
        ⬝ᵥ (fun q => m.QmatHetR w N ω q l)))⁻¹) sN⁻¹ :=
      (m.tendstoInProb_colGram w c hw hR hreg hG hpd l).comp_continuous
        (φ := fun x => (Real.sqrt x)⁻¹) (Real.continuous_sqrt.continuousAt.inv₀ hsN.ne')
    exact (h3.mul h2) η hη0
  have hedgeC7 : Tendsto (fun N => μ N
      {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
        (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k ≤ τ}ᶜ)
      atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_eigenvalues₀_le_het w N u τ).nullMeasurableSet)
      (m.tendsto_measure_eigenvalues₀_le_tau_het w c hc hw hR hreg hG hpd hedge e hτ hsubE')
  -- 6. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + (τ - MPhet.bHet c w) / 2}ᶜ
        ∪ ((⋃ l' : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l')
                ⬝ᵥ (fun q => m.QmatHetR w N ω q l') - m.colGramLimit w c l'|})
          ∪ ({ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q l)
                ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimit w c l|}
            ∪ ((⋃ q : Fin r × Fin u, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
                  (fun i => m.QmatHetR w N ω i q.1) (fun i => m.QmatHetR w N ω i (f q.2))
                  - (if q.1 = f q.2 then -1 else 0)|})
              ∪ ((⋃ q : Fin u × Fin u, {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv q.1)
                    (fun i => m.QmatHetR w N ω i (f q.1))
                    (fun i => m.QmatHetR w N ω i (f q.2))
                    - (if q.1 = q.2 then νv q.1 else 0)|})
                ∪ ((⋃ k : Fin u, {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q l)
                        ⬝ᵥ (fun q => m.QmatHetR w N ω q l)))⁻¹
                      * R4.cform (m.W0hetR w N ω) (ρv k)
                        (fun q => m.QmatHetR w N ω q l)
                        (fun q => m.QmatHetR w N ω q (f k)) - tv k|})
                  ∪ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))), u ≤ (k : ℕ) →
                      (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
                        ≤ τ}ᶜ)))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt, not_not] at hbad
    obtain ⟨hg1, hg2, hg3, hg4, hg5, hg6, hg7⟩ := hbad
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by
      linarith
    have hsq : ∀ l' : Fin r, ∑ i, m.QmatHetR w N ω i l' ^ 2
        = (fun q => m.QmatHetR w N ω q l') ⬝ᵥ (fun q => m.QmatHetR w N ω q l') :=
      fun l' => Finset.sum_congr rfl fun i _ => sq _
    have hcolb : ∀ l' : Fin r, ∑ i, m.QmatHetR w N ω i l' ^ 2 ≤ Cq := by
      intro l'
      refine le_trans ?_ (hCqle l')
      rw [hsq l']
      have h1 := (abs_lt.mp (hg2 l')).2
      linarith
    have hg : 0 < (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) := by
      have h1 := (abs_lt.mp hg3).1
      linarith
    have hgle : (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        ≤ m.colGramLimit w c l + 1 := by
      have h1 := (abs_lt.mp (hg2 l)).2
      linarith
    have hI := OutliersR.count_eq_of_forms_of_edge (Q := m.QmatHetR w N ω)
      (τ := τ) (mg := mg) (Cq := Cq) (η := η)
      (m.isHermitian_W0hetR w N ω) (isHermitian_mul_transpose_self (m.stackXW w N ω))
      (m.gram_eq_hetR w N ω) hfinj hmg
      hlamτ hτρ hνpos hCq0 hη0 hcolb (fun l' k => (hg4 (l', k)).le)
      (fun k l' => (hg5 (k, l')).le) hsmall hδm hg7
    have hE3' : ∀ k : Fin u, |R4.cform (m.W0hetR w N ω) (ρv k)
        ((Real.sqrt ((fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)))⁻¹
          • fun q => m.QmatHetR w N ω q l)
        (fun q => m.QmatHetR w N ω q (f k)) - tv k| ≤ η := by
      intro k
      rw [R4.cform_smul_left]
      exact (hg6 k).le
    have hdet := OutliersR.align_detG (Q := m.QmatHetR w N ω)
      (τ := τ) (mg := mg) (Cq := Cq) (t := tv) (Lb := sN⁻¹) (η := η)
      (m.isHermitian_W0hetR w N ω) (isHermitian_mul_transpose_self (m.stackXW w N ω))
      (m.gram_eq_hetR w N ω) hfinj hmg
      hlamτ hτρ hνpos hI hCq0 hcolb (OutliersR.dot_self_normalize hg) hLb hη0 hη1
      (fun l' k => (hg4 (l', k)).le) (fun k l' => (hg5 (k, l')).le) hE3' hsmall
    have hclose := OutliersR.abs_mul_sub_le_of_close hg hgle hdet hg3.le
    have hω' : ε ≤ |‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (WithLp.toLp 2 fun q => m.QmatHetR w N ω q l)‖ ^ 2
        - (if Scalars.Assumption4 (fun i => m.thetaAligned i l) c w
            ∧ τ < MPhet.rhoHet (fun i => m.thetaAligned i l) c w
          then 1 / MPhet.nuHet (fun i => m.thetaAligned i l) c w else 0)| := hω
    rw [OutliersR.normSq_specProj_eq_dot_mul _ _ hg, ← htarget] at hω'
    linarith
  · exact tendsto_measure_zero_union hedgeC1
      (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
        (tendsto_measure_zero_union hcol2T
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
              (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE3T)
                hedgeC7)))))

end UnalignedModelR

end StackedSVD
