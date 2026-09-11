/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.WeightedUpper
import StackedSVD.RankR.GeneralFrob
import StackedSVD.RankR.GeneralMain
import StackedSVD.RankR.Aligned
import StackedSVD.RMT.Full

/-!
# E2 at general `r_i`, and the aggregate clause of `thm:rank_r_svdstack`

`RankR/WeightedUpper.lean` proves the weight-free upper bound of `thm:gen_rank_weight_svdstak`
(`main_paper.tex:893`) on an `UnalignedModel`, that is at `r_i = 1`. This file is its twin on
an `UnalignedModelR`, where table `i` carries `r_i` spikes and the row index runs over the
flat block index `p ∈ [r̃]`, `r̃ = ∑_i r_i`.

Every deterministic lemma of `RankR/WeightedUpper.lean` is stated on plain matrices or on the
flat objects `BR`, `AbetaR`, `AbetaRW`, `optWR`, `limitRTrace`, `limitRWk`, so it applies here
at table count `r̃` through the `rfl` bridges of `RankR/Flatten.lean`. Only the probabilistic
lemmas are restated, because `perfRGW` reads `UnalignedModelR` and not `UnalignedModel`.

## Content

1. `rowBoundRG` and `perfRGW_le_rowBoundRG`: the weight-free bound at general `r_i`.
2. `rbFRG`, `rowBoundRG_tendsto`, `tendsto_measure_det_gramG_eq_zero`: the bound converges.
3. `perfRGWk`, `perfRGWk_tendsto`, `perfRGWk_le_perfRGW` and
   `thm_gen_rank_weight_svdstak_general_r_norank`: attainment at `W⋆` with no hypothesis on
   `rank B_R`.
4. `perfRGW_uniform_bound`, `thm_gen_rank_weight_svdstak_general_r_full` and its Gaussian
   facade at `r_i = 1`.
5. **Track D item D3**: `limitOptG_aligned` and `thm_rank_r_svdstack_aggregate`, the aggregate
   clause `‖Vᵀ V̂_svdstack‖_F² → ∑_j S_j / (S_j + 1)` of `thm:rank_r_svdstack`
   (`main_paper.tex:2401`).

The rank-one mirror of each declaration is named in its docstring.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

/-! ### 1. The weight-free bound on an `UnalignedModelR`

`rowBoundRG` is the right side of `trace_specInvTop_conj_le_inv` (`RankR/WeightedUpper.lean`)
at `G = Ṽ Ṽᵀ` and `g = Ṽ V`, with `Ṽ` the `r̃ × d` matrix of the per-spike estimates. It reads
no weight matrix, so `perfRGW_le_rowBoundRG` covers data dependent weights. -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The weight-free upper bound `tr((Ṽ V)ᵀ (Ṽ Ṽᵀ)⁻¹ (Ṽ V))` at general `r_i`
(`paper_edits.md`, item E2). Rank-one mirror: `UnalignedModel.rowBoundR`
(`RankR/WeightedUpper.lean`). -/
noncomputable def rowBoundRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVG N ω)ᵀ * (m.gramG N ω)⁻¹ * m.VtVG N ω)

/-- **The uniform bound, deterministic form, at general `r_i`.** At every `N` and every `ω`
where `Ṽ Ṽᵀ` is invertible, no weight matrix beats `rowBoundRG`. Rank-one mirror:
`UnalignedModel.perfRW_le_rowBoundR`. -/
theorem perfRGW_le_rowBoundRG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    (hdet : IsUnit (m.gramG N ω).det) :
    m.perfRGW W N ω ≤ m.rowBoundRG N ω := by
  have hS : (W * m.gramG N ω * Wᵀ).IsHermitian := by
    have h := Matrix.isHermitian_mul_mul_conjTranspose (A := m.gramG N ω) W
      (m.isHermitian_gramG N ω)
    rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have h := trace_specInvTop_conj_le_inv (m.posSemidef_gramG N ω) hdet W hS (m.VtVG N ω) r
  rw [UnalignedModelR.perfRGW, UnalignedModelR.rowBoundRG,
    specInvTop_congr_mat (m.isHermitian_gramWG W N ω) hS (m.gramWG_eq W N ω) r,
    m.VtVWG_eq W N ω]
  exact h

/-! ### 2. The bound converges

`rbFRG` is the random point of the continuous functional `rbPhiR` (`RankR/WeightedUpper.lean`)
at table count `r̃`; its limit point is `rbLimR (betaFlat β) (Rcol m.R)`. -/

/-- The random point: `Ṽ Ṽᵀ` on the left block, `Ṽ V` on the right block. Rank-one mirror:
`UnalignedModel.rbFR`. -/
noncomputable def rbFRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    (Fin (rtot rk) × Fin (rtot rk)) ⊕ (Fin (rtot rk) × Fin r) → ℝ :=
  Sum.elim (fun t : Fin (rtot rk) × Fin (rtot rk) => m.gramG N ω t.1 t.2)
    (fun t : Fin (rtot rk) × Fin r => m.VtVG N ω t.1 t.2)

theorem rbMatR_rbFRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    rbMatR (m.rbFRG N ω) = m.gramG N ω := rfl

theorem rbPhiR_rbFRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    rbPhiR (m.rbFRG N ω) = m.rowBoundRG N ω :=
  rbPhiR_eq (m.gramG N ω) (m.VtVG N ω) (m.rbFRG N ω) (fun _ _ => rfl) (fun _ _ => rfl)

/-- `gramR_general` and `VtV_tendsto_general` in one family. Rank-one mirror:
`UnalignedModel.tendstoInProbPi_rbFR`. -/
theorem tendstoInProbPi_rbFRG (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProbPi μ (fun N ω => m.rbFRG N ω) (rbLimR (betaFlat β) (Rcol m.R)) := by
  rintro (t | t)
  · exact m.gramR_general c β hβdef law hG t.1 t.2
  · exact m.VtV_tendsto_general c β hβdef law t.1 t.2

/-- The determinant of `Ṽ Ṽᵀ` converges in probability to `det A_{β,R}`. Rank-one mirror:
`UnalignedModel.det_gram_tendsto`. -/
theorem det_gramG_tendsto (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => (m.gramG N ω).det) ((ABlock β m.R).det) := by
  have hcont : ContinuousAt
      (fun z : (Fin (rtot rk) × Fin (rtot rk)) ⊕ (Fin (rtot rk) × Fin r) → ℝ =>
        (rbMatR z).det) (rbLimR (betaFlat β) (Rcol m.R)) :=
    continuous_rbMatR.matrix_det.continuousAt
  have h := TendstoInProbPi.comp_continuous hcont (m.tendstoInProbPi_rbFRG c β hβdef law hG)
  rw [rbMatR_rbLimR] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  rw [m.rbMatR_rbFRG N ω]

/-- The null determinant event vanishes, in the `IsUnit` form the squeeze consumes. Rank-one
mirror: `UnalignedModel.tendsto_measure_det_gram_eq_zero`. -/
theorem tendsto_measure_det_gramG_not_isUnit (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ IsUnit (m.gramG N ω).det}) atTop (𝓝 0) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i j, 0 ≤ β i j := fun i j => (hβ01 i j).1
  have h1 : ∀ i j, β i j < 1 := fun i j => (hβ01 i j).2
  have hpos : 0 < (ABlock β m.R).det := (ABlock_posDef m.R h0 h1).det_pos
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | (ABlock β m.R).det ≤ |(m.gramG N ω).det - (ABlock β m.R).det|}) ?_
    ((m.det_gramG_tendsto c β hβdef law hG) _ hpos)
  intro N ω hω
  have hω' : (m.gramG N ω).det = 0 := by
    by_contra hne
    exact hω (isUnit_iff_ne_zero.mpr hne)
  change (ABlock β m.R).det ≤ |(m.gramG N ω).det - (ABlock β m.R).det|
  rw [hω', zero_sub, abs_neg, abs_of_pos hpos]

/-- The null determinant event vanishes, in the `det = 0` form. -/
theorem tendsto_measure_det_gramG_eq_zero (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | (m.gramG N ω).det = 0}) atTop (𝓝 0) := by
  refine tendsto_measure_zero_of_subset (t := fun N => {ω | ¬ IsUnit (m.gramG N ω).det}) ?_
    (m.tendsto_measure_det_gramG_not_isUnit c β hc hβdef law hG)
  intro N ω hω
  have hω' : (m.gramG N ω).det = 0 := hω
  exact fun hu => (isUnit_iff_ne_zero.mp hu) hω'

/-- **The bound converges** at general `r_i` (`paper_edits.md`, item E2). No eigengap and no
rank hypothesis: `rbPhiR` is continuous at the limit point because `A_{β,R}` is positive
definite. Rank-one mirror: `UnalignedModel.rowBoundR_tendsto`. -/
theorem rowBoundRG_tendsto (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.rowBoundRG N ω) (limitRTrace (betaFlat β) (Rcol m.R)) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hpos : 0 < (AbetaR (betaFlat β) (Rcol m.R)).det := (abetaR_posDef _ h0 h1).det_pos
  have hphi : ContinuousAt (rbPhiR (M := rtot rk) (r := r))
      (rbLimR (betaFlat β) (Rcol m.R)) := by
    refine continuousAt_rbPhiR (rbLimR (betaFlat β) (Rcol m.R)) ?_
    rw [rbMatR_rbLimR]
    exact hpos.ne'
  have h := TendstoInProbPi.comp_continuous hphi (m.tendstoInProbPi_rbFRG c β hβdef law hG)
  rw [rbPhiR_rbLimR] at h
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact m.rbPhiR_rbFRG N ω

/-! ### 3. Attainment at `W⋆` with no hypothesis on `rank B_R`

The truncated performance `perfRGWk W k` uses `k` eigenvectors instead of `r`. It is monotone
in `k` and it converges whenever `W A_{β,R} Wᵀ` has a gap at the index `k`. At `W⋆` the gap at
`k = rank B_R` holds with no hypothesis (`topGap_optWR_rank`, `RankR/WeightedUpper.lean`, read
at table count `r̃`). -/

/-- The performance of weighted svdstack at general `r_i`, truncated at `k` eigenvectors
instead of `r`. At `k = r` it is `perfRGW`. Rank-one mirror: `UnalignedModel.perfRWk`. -/
noncomputable def perfRGWk (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (k : ℕ) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVWG W N ω)ᵀ *
    specInvTop (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) k * m.VtVWG W N ω)

theorem perfRGWk_r (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    m.perfRGWk W r N ω = m.perfRGW W N ω := rfl

/-- The truncated performance is monotone in the truncation index. Rank-one mirror:
`UnalignedModel.perfRWk_le_perfRW`. -/
theorem perfRGWk_le_perfRGW (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) {k : ℕ} (hk : k ≤ r) (N : ℕ) (ω : Ω N) :
    m.perfRGWk W k N ω ≤ m.perfRGW W N ω :=
  trace_specInvTop_conj_mono (m.isHermitian_gramWG W N ω) (m.posSemidef_gramWG W N ω) hk
    (m.VtVWG W N ω)

/-- `thm_gen_rank_weight_svdstak_general_r_conv` (`RankR/GeneralMain.lean`) with the truncation
index `k` free. Same proof; only the index inside `continuousAt_traceFun` changes. Rank-one
mirror: `UnalignedModel.perfRWk_tendsto`. -/
theorem perfRGWk_tendsto (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    {k : ℕ} (hk : 0 < k) (hkr : k ≤ rtot rk)
    (hgapW : TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) k)
    (hposW : 0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
      ⟨k - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGWk W k N ω)
      (limitRWk W (betaFlat β) (Rcol m.R) k) := by
  have hA : (ABlockW W β m.R).IsHermitian := isHermitian_ABlockW W β m.R
  have hkp : k ≤ Fintype.card (Fin (rtot rk)) := by simpa using hkr
  have hval : limitRWk W (betaFlat β) (Rcol m.R) k
      = Matrix.trace ((BBlockW W β m.R)ᵀ * specInvTop (ABlockW W β m.R) hA k *
        BBlockW W β m.R) := rfl
  rw [hval]
  have hconv0 : TendstoInProbPi μ
      (fun N ω => Sum.elim
        (fun t : Fin (rtot rk) × Fin (rtot rk) => m.gramG N ω t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => m.VtVG N ω t.1 t.2))
      (Sum.elim (fun t : Fin (rtot rk) × Fin (rtot rk) => ABlock β m.R t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => BBlock β m.R t.1 t.2)) := by
    rintro (t | t)
    · exact m.gramR_general c β hβdef law hG t.1 t.2
    · exact m.VtV_tendsto_general c β hβdef law t.1 t.2
  have hconv : TendstoInProbPi μ
      (fun N ω => Sum.elim
        (fun t : Fin (rtot rk) × Fin (rtot rk) => m.gramWG W N ω t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => m.VtVWG W N ω t.1 t.2))
      (Sum.elim (fun t : Fin (rtot rk) × Fin (rtot rk) => ABlockW W β m.R t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => BBlockW W β m.R t.1 t.2)) := by
    rintro (t | t)
    · obtain ⟨i, j⟩ := t
      have hcont : Continuous fun z :
          (Fin (rtot rk) × Fin (rtot rk)) ⊕ (Fin (rtot rk) × Fin r) → ℝ =>
            ∑ b : Fin (rtot rk), ∑ a : Fin (rtot rk), W i a * z (Sum.inl (a, b)) * W j b :=
        continuous_finsetSum _ fun b _ => continuous_finsetSum _ fun a _ =>
          (continuous_const.mul (continuous_apply _)).mul continuous_const
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inl] at h
      have hlim : ABlockW W β m.R i j
          = ∑ b : Fin (rtot rk), ∑ a : Fin (rtot rk), W i a * ABlock β m.R a b * W j b :=
        UnalignedModel.mul_mul_transpose_apply W (ABlock β m.R) i j
      rw [Sum.elim_inl, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inl]
      rw [m.gramWG_eq, UnalignedModel.mul_mul_transpose_apply]
    · obtain ⟨i, kk⟩ := t
      have hcont : Continuous fun z :
          (Fin (rtot rk) × Fin (rtot rk)) ⊕ (Fin (rtot rk) × Fin r) → ℝ =>
            ∑ a : Fin (rtot rk), W i a * z (Sum.inr (a, kk)) :=
        continuous_finsetSum _ fun a _ => continuous_const.mul (continuous_apply _)
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inr] at h
      have hlim : BBlockW W β m.R i kk = ∑ a : Fin (rtot rk), W i a * BBlock β m.R a kk :=
        Matrix.mul_apply
      rw [Sum.elim_inr, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inr]
      rw [m.VtVWG_eq]
      exact (Matrix.mul_apply).symm
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hk hkp hgapW hposW (BBlockW W β m.R)) hconv
  rw [traceFun_eq hA (BBlockW W β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gramWG W N ω) (m.VtVWG W N ω)

/-- **`thm:gen_rank_weight_svdstak`, attainment half, at general `r_i` and with no rank
hypothesis** (`paper_edits.md`, item E2, "Attainment"). The performance at `W⋆` is squeezed
between the truncated performance at `k = rank B_R`, whose gap is automatic, and `rowBoundRG`,
which converges with no gap at all. Rank-one mirror:
`UnalignedModel.thm_gen_rank_weight_svdstak_norank`. -/
theorem thm_gen_rank_weight_svdstak_general_r_norank (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hrr : r ≤ rtot rk)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
      (limitOptG β m.R (by simpa using hrr)) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
  have hLeq : limitOptG β m.R hrp = limitRTrace (betaFlat β) (Rcol m.R) := by
    rw [limitOptG_eq_limitROpt]
    exact (limitRTrace_eq_limitROpt (betaFlat β) (Rcol m.R) h0 h1 hrp).symm
  rw [hLeq]
  have hlow : TendstoInProb μ
      (fun N ω => m.perfRGWk (optWG β) (BR (betaFlat β) (Rcol m.R)).rank N ω)
      (limitRTrace (betaFlat β) (Rcol m.R)) := by
    rcases Nat.eq_zero_or_pos (BR (betaFlat β) (Rcol m.R)).rank with hk0 | hkpos
    · have hzero : ∀ (N : ℕ) (ω : Ω N),
          m.perfRGWk (optWG β) (BR (betaFlat β) (Rcol m.R)).rank N ω = 0 := by
        intro N ω
        rw [UnalignedModelR.perfRGWk, hk0, specInvTop_zero]
        simp
      have hlim : limitRTrace (betaFlat β) (Rcol m.R) = 0 := by
        rw [limitRTrace_eq_limitRWk_optWR (betaFlat β) (Rcol m.R) h0 h1 (k := 0)
            (le_of_eq hk0), limitRWk, specInvTop_zero]
        simp
      rw [hlim]
      exact (TendstoInProb.const μ 0).congr fun N =>
        Filter.Eventually.of_forall fun ω => (hzero N ω).symm
    · have hkr : (BR (betaFlat β) (Rcol m.R)).rank ≤ rtot rk := by
        simpa using Matrix.rank_le_height (BR (betaFlat β) (Rcol m.R))
      have h := m.perfRGWk_tendsto c β (optWG β) hβdef hkpos hkr
        (topGap_optWR_rank (Rcol m.R) h0 h1) (pos_eigenvalues₀_optWR (Rcol m.R) h0 h1 _) law hG
      have hcast : limitRWk (optWG β) (betaFlat β) (Rcol m.R)
          (BR (betaFlat β) (Rcol m.R)).rank = limitRTrace (betaFlat β) (Rcol m.R) :=
        (limitRTrace_eq_limitRWk_optWR (betaFlat β) (Rcol m.R) h0 h1 le_rfl).symm
      rwa [hcast] at h
  have hup : TendstoInProb μ (fun N ω => m.rowBoundRG N ω)
      (limitRTrace (betaFlat β) (Rcol m.R)) := m.rowBoundRG_tendsto c β hc hβdef law hG
  have hnull := m.tendsto_measure_det_gramG_not_isUnit c β hc hβdef law hG
  refine tendstoInProb_of_subset_union₃ ?_
  intro δ hδ
  refine ⟨fun N => {ω | δ ≤ |m.perfRGWk (optWG β) (BR (betaFlat β) (Rcol m.R)).rank N ω -
      limitRTrace (betaFlat β) (Rcol m.R)|},
    fun N => {ω | δ ≤ |m.rowBoundRG N ω - limitRTrace (betaFlat β) (Rcol m.R)|},
    fun N => {ω | ¬ IsUnit (m.gramG N ω).det}, ?_, hlow δ hδ, hup δ hδ, hnull⟩
  intro N ω hω
  have hωd : δ ≤ |m.perfRGW (optWG β) N ω - limitRTrace (betaFlat β) (Rcol m.R)| := hω
  simp only [Set.mem_union]
  by_cases hdet : IsUnit (m.gramG N ω).det
  · by_cases hA1 : δ ≤ |m.perfRGWk (optWG β) (BR (betaFlat β) (Rcol m.R)).rank N ω -
        limitRTrace (betaFlat β) (Rcol m.R)|
    · exact Or.inl (Or.inl hA1)
    · by_cases hA2 : δ ≤ |m.rowBoundRG N ω - limitRTrace (betaFlat β) (Rcol m.R)|
      · exact Or.inl (Or.inr hA2)
      · exfalso
        have hA1' := not_le.mp hA1
        have hA2' := not_le.mp hA2
        have hle1 : m.perfRGWk (optWG β) (BR (betaFlat β) (Rcol m.R)).rank N ω
            ≤ m.perfRGW (optWG β) N ω :=
          m.perfRGWk_le_perfRGW (optWG β) (rank_BR_le (betaFlat β) (Rcol m.R)) N ω
        have hle2 : m.perfRGW (optWG β) N ω ≤ m.rowBoundRG N ω :=
          m.perfRGW_le_rowBoundRG (optWG β) N ω hdet
        rw [abs_lt] at hA1' hA2'
        have hlt : |m.perfRGW (optWG β) N ω - limitRTrace (betaFlat β) (Rcol m.R)| < δ := by
          rw [abs_lt]
          exact ⟨by linarith [hA1'.1], by linarith [hA2'.2]⟩
        linarith
  · exact Or.inr hdet

/-! ### 4. The statements of item E2 at general `r_i`

`perfRGW_uniform_bound` is the sentence proposed for the paper: the bound
`‖V̂(W)ᵀ V‖_F² ≤ tr((Ṽ V)ᵀ (Ṽ Ṽᵀ)⁻¹ Ṽ V)` holds for every `W` at every `N`, so `W⋆` is optimal
among all weightings, data dependent ones included. No set is assumed measurable; only
monotonicity and subadditivity of `μ N` enter. -/

/-- **E2, the uniform bound, at general `r_i`.** With probability tending to one no weight
matrix beats `L⋆ + ε`. Rank-one mirror: `UnalignedModel.perfRW_uniform_bound`. -/
theorem perfRGW_uniform_bound (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hrr : r ≤ rtot rk)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
      limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W N ω}) atTop (𝓝 0) := by
  intro ε hε
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
  have hLeq : limitOptG β m.R hrp = limitRTrace (betaFlat β) (Rcol m.R) := by
    rw [limitOptG_eq_limitROpt]
    exact (limitRTrace_eq_limitROpt (betaFlat β) (Rcol m.R) h0 h1 hrp).symm
  have hrb := m.rowBoundRG_tendsto c β hc hβdef law hG
  have hnull := m.tendsto_measure_det_gramG_not_isUnit c β hc hβdef law hG
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ¬ IsUnit (m.gramG N ω).det} ∪
      {ω | ε ≤ |m.rowBoundRG N ω - limitRTrace (betaFlat β) (Rcol m.R)|}) ?_
    (tendsto_measure_zero_union hnull (hrb ε hε))
  intro N ω hω
  obtain ⟨W, hW⟩ := hω
  simp only [Set.mem_union]
  by_cases hdet : IsUnit (m.gramG N ω).det
  · right
    have hb := m.perfRGW_le_rowBoundRG W N ω hdet
    have hW' : limitRTrace (betaFlat β) (Rcol m.R) + ε ≤ m.perfRGW W N ω := by
      rw [← hLeq]
      exact hW
    change ε ≤ |m.rowBoundRG N ω - limitRTrace (betaFlat β) (Rcol m.R)|
    have hle : ε ≤ m.rowBoundRG N ω - limitRTrace (betaFlat β) (Rcol m.R) := by linarith
    exact le_trans hle (le_abs_self _)
  · left
    exact hdet

/-- The per-`W` form of the uniform bound, for one fixed weight matrix. Rank-one mirror:
`UnalignedModel.perfRW_le_opt_whp`. -/
theorem perfRGW_le_opt_whp (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hrr : r ≤ rtot rk)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W N ω}) atTop (𝓝 0) := by
  intro ε hε
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ∃ W' : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
      limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W' N ω}) ?_
    (m.perfRGW_uniform_bound c β hc hβdef hrr law hG ε hε)
  intro N ω hω
  exact ⟨W, hω⟩

/-- **`thm:gen_rank_weight_svdstak` at general `r_i`, in the form of item E2.** Three
conclusions in one declaration and **no hypothesis on `rank B_R`**:

1. the performance at `W⋆ = D^{-1/2}` tends to `L⋆`;
2. for every admissible `W` the performance tends to `limitRGW W ≤ L⋆`;
3. with probability tending to one no weight matrix at all beats `L⋆ + ε`.

Conjunct 3 replaces the paper's footnote (`main_paper.tex:888`), which excludes every `W`
without an eigengap. Rank-one mirror: `UnalignedModel.thm_gen_rank_weight_svdstak_full`. -/
theorem thm_gen_rank_weight_svdstak_general_r_full (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hr : 0 < r) (hrr : r ≤ rtot rk)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (by simpa using hrr)) ∧
      (∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (by simpa using hrr)) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W N ω}) atTop (𝓝 0) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
  refine ⟨m.thm_gen_rank_weight_svdstak_general_r_norank c β hc hβdef hrr law hG,
    fun W hgapW hposW => ⟨m.thm_gen_rank_weight_svdstak_general_r_conv c β W hβdef hr hrr
      hgapW hposW law hG, ?_⟩,
    m.perfRGW_uniform_bound c β hc hβdef hrr law hG⟩
  rw [limitRGW_eq_limitRW, limitOptG_eq_limitROpt]
  exact limitRW_le_opt W (betaFlat β) (Rcol m.R) h0 h1 hr hrp hgapW hposW

end UnalignedModelR

/-! ### 5. The Gaussian facade at one spike per table

The chain is `RankR/GeneralFrob.lean`: `SpikedModelR.norm_ucol`, `SpikedModelR.toSpiked`,
`SpikedModelR.toSpiked_X`, `SpikedModelR.tableLawR_of_gaussian` and
`UnalignedModelR.tableLawR_of_gaussian`. This file used to copy all five under primed names,
because `RankR/GeneralFrob.lean` imports `RankR/Frobenius.lean`. The two files are siblings in
the import graph, so the import is legal; the second cleanup pass (2026-09-02) added it and
deleted the copies. -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **Item E2, Layer 2 form, at one spike per table.** The three conclusions of
`thm_gen_rank_weight_svdstak_general_r_full` from the proportional regime of each table and the
joint Gaussian law, with no `TableLawR` hypothesis. Rank-one mirror:
`UnalignedModel.thm_gen_rank_weight_svdstak_full_gaussian`. -/
theorem thm_gen_rank_weight_svdstak_general_r_full_gaussian_one
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r fun _ => 1)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin 1 → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r)
    (hrr : r ≤ rtot (fun _ => 1 : Fin M → ℕ)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (by simpa using hrr)) ∧
      (∀ W : Matrix (Fin (rtot (fun _ => 1 : Fin M → ℕ)))
          (Fin (rtot (fun _ => 1 : Fin M → ℕ))) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (by simpa using hrr)) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin (rtot (fun _ => 1 : Fin M → ℕ)))
        (Fin (rtot (fun _ => 1 : Fin M → ℕ))) ℝ,
        limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W N ω}) atTop (𝓝 0) :=
  m.thm_gen_rank_weight_svdstak_general_r_full c β hc hβdef hr hrr
    (m.tableLawR_of_gaussian hc hreg hG) hG.indepNoise

end UnalignedModelR

/-! ### 6. Track D item D3: the aggregate clause of `thm:rank_r_svdstack`

`main_paper.tex:2401`. Under the rank-`r` model with every table aligned (`r_i = r`, `R_i = I`)
the svdstack limit is `∑_j S_j / (S_j + 1)` with `S_j = ∑_i β_ij² / (1 - β_ij²)`. The
deterministic half is `trace_conj_inv_AbetaR_aligned` (`RankR/Aligned.lean`); the probabilistic
half is `thm_gen_rank_weight_svdstak_general_r_norank` of section 3, which carries no rank
hypothesis and so needs no `S_j > 0`. -/

section Aligned

variable {M r : ℕ}

/-- `r̃ = M r` in the exactly aligned model. -/
theorem rtot_alignedRk : rtot (alignedRk M r) = M * r := by
  rw [rtot, Finset.sum_const, Finset.card_univ, Fintype.card_fin, smul_eq_mul]

/-- `r ≤ r̃` in the exactly aligned model, from one table. -/
theorem le_rtot_alignedRk (hM : 0 < M) : r ≤ rtot (alignedRk M r) :=
  Finset.single_le_sum (f := fun _ : Fin M => r) (fun _ _ => Nat.zero_le r)
    (Finset.mem_univ (⟨0, hM⟩ : Fin M))

/-- **`L⋆` in the exactly aligned model** (`main_paper.tex:2401`, `paper_edits.md` E4 item 1):
`L⋆ = ∑_j S_j / (S_j + 1)`. Route: `limitOptG_eq_limitROpt` (`RankR/Flatten.lean`), then
`limitRTrace_eq_limitROpt` (`RankR/WeightedUpper.lean`) backwards, then
`trace_conj_inv_AbetaR_aligned` (`RankR/Aligned.lean`). -/
theorem limitOptG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1)
    (hrr : r ≤ Fintype.card (Fin (rtot (alignedRk M r)))) :
    limitOptG (rk := alignedRk M r) β (alignedR M r) hrr = ∑ j, Sagg β j / (Sagg β j + 1) := by
  have h0' : ∀ p, 0 ≤ betaFlat (rk := alignedRk M r) β p :=
    fun p => h0 (blk p).1 (blk p).2
  have h1' : ∀ p, betaFlat (rk := alignedRk M r) β p < 1 :=
    fun p => h1 (blk p).1 (blk p).2
  rw [limitOptG_eq_limitROpt,
    ← limitRTrace_eq_limitROpt (betaFlat (rk := alignedRk M r) β) (Rcol (alignedR M r))
      h0' h1' hrr]
  exact trace_conj_inv_AbetaR_aligned β h0 h1

end Aligned

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **`thm:rank_r_svdstack`, aggregate clause** (`main_paper.tex:2401`, Track D item D3). Under
the exactly aligned rank-`r` model the weighted svdstack performance at `W⋆` converges in
probability to `∑_j S_j / (S_j + 1)`, `S_j = ∑_i β_ij² / (1 - β_ij²)`.

`perfRGW (optWG β)` is the paper's `‖Vᵀ V̂_svdstack‖_F²` whenever the top-`r` eigenframe of
`Ṽ_{W⋆} Ṽ_{W⋆}ᵀ` exists (`frobSq_vhatSvdstackGW`, `RankR/GeneralFrob.lean`).

No hypothesis on `S_j`: `thm_gen_rank_weight_svdstak_general_r_norank` needs no rank of `B_R`,
so a component with `S_j = 0` is allowed and contributes `0` to the sum. `hM : 0 < M` gives
`r ≤ r̃`; `hR` is the paper's `R_i = I_r`. -/
theorem thm_rank_r_svdstack_aggregate (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
      (∑ j, Sagg β j / (Sagg β j + 1)) := by
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i j, 0 ≤ β i j := fun i j => (hβ01 i j).1
  have h1 : ∀ i j, β i j < 1 := fun i j => (hβ01 i j).2
  have hrr : r ≤ rtot (alignedRk M r) := le_rtot_alignedRk hM
  have hRe : m.R = alignedR M r := funext hR
  have h := m.thm_gen_rank_weight_svdstak_general_r_norank c β hc hβdef hrr law hG
  rw [hRe, limitOptG_aligned β h0 h1] at h
  exact h

end UnalignedModelR

end StackedSVD
