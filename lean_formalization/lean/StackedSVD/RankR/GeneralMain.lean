/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.GramR
import StackedSVD.RankR.Flatten
import StackedSVD.RankR.Weighted
import StackedSVD.SVDStack.Deterministic

/-!
# Section 7 at general `r_i`: the five statements and the `Ṽ V` limit

STATUS 2026-09-01: all five statements are proved, 0 `sorry`.
`lem_general_rank_delocalization_general` (task T4), `VtV_tendsto_general` (T6),
`gramR_general` (T5), `prop_general_rank_unweighted_svdstack_general` (T7) and
`thm_gen_rank_weight_svdstak_general_r` (T8 of `notes/archive/rankr_plan_B.md`).

The five statements moved here verbatim from `RankR/General.lean` on 2026-09-01. Their proofs
need `RankR/GramR.lean` (the per-table delocalization twins of task B2b), and `GramR.lean`
imports `General.lean`, so the statements cannot stay where they were stated.

## Content

1. the finite sum of limits in probability comes from `TendstoInProb.finsum`
   (`Prob/TendstoInProb.lean`) since the 2026-09-02 dedupe.
2. `SpikedModelR.rk_le_d`: a table has at most `d N` spikes, so the eigenvalue index of a spike
   is in range for `vEig` and `specProjIdx`, which carry a junk value out of range.
3. `norm_vhatG`, `mem_specSpace_vhatG`, `inner_vhatG_nonneg`, `inner_vhatG_sq_le`,
   `overlapIdx_eq_inner_vhatG_sq`, `inner_vhatG_intra`, `inner_vhatG_blk_intra`: the
   interface of the signed singular vector `v̂_ij`.
   Mirrors: `norm_vhat`, `mem_topSpace_vhat`, `inner_vhat_nonneg` of `RankR/Defs.lean`, with
   `vEig` and `specProjIdx` for `vMax` and `topProj`.
4. `VtVG_eq_inner`, `gramG_eq_inner`, `inner_col_colVecG`, `norm_colVecG`: the entries of
   `Ṽ V` and of `Ṽ Ṽᵀ`, and the shared frame `V`. Mirrors: `VtV_eq_inner`, `gram_eq_inner`,
   `inner_v_colVec`, `norm_colVec`.
5. `align_innerR`, `cross_innerR`, `inner_vhatG_perp_tendsto`, `align_inner_detR`: the three
   fields of `TableLawR` in signed inner-product form, and the split of a deterministic
   direction against the signal frame of one table. Mirrors: `align_inner`,
   `inner_vhat_perp_tendsto`, `align_inner_det`.
6. `gramWG_eq`, `VtVWG_eq`: the weighted Gram matrix and the weighted overlap matrix as
   fixed linear images of the unweighted ones. Mirrors: `gramW_eq`, `VtVW_eq` of
   `RankR/Weighted.lean`.
7. The five statements, plus `thm_gen_rank_weight_svdstak_general_r_conv`, the weighted limit
   at one admissible weight matrix.

## The route of `VtV_tendsto_general`

`⟪v̂_ij, V e_k⟫` splits as the paper does (`main_paper.tex:1961`), with the **whole signal
frame** of table `i` in place of the single spike direction of the `r_i = 1` case:

```
V e_k = ∑_{l} ⟪(V R_i)_l, V e_k⟫ (V R_i)_l + perpSpan_i (V e_k).
```

`TableLawR.align j` gives the term `l = j`, `TableLawR.cross j l` kills the terms `l ≠ j`, and
`TableLawR.delocUniform j` kills the remainder, which is orthogonal to the whole span. The
coefficients are `⟪(V R_i)_l, V e_k⟫ = (R_i)_{kl}` (`inner_col_colVecG`), so the limit is
`β_ij (R_i)_{kj}`, the entry `(p, k)` of `B_R`. No independence across tables enters, so `hG`
does not appear.

The sign convention of `vhatG` is what makes the limit signed rather than squared:
`inner_vhatG_nonneg` turns `√(overlapIdx)` into the inner product itself.

`rk_le_d` and `norm_colVecG` moved out to `RankR/General.lean`;
`tableLawR_of_singleTableLaw` and `SpikedModel.tableLawR_of_singleTableLaw'` moved in
from `RankR/General.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

/-- Two eigenvectors of one Hermitian matrix at two different sorted indices are orthogonal.
`vEig` reads `Matrix.IsHermitian.eigenvectorBasis`, which is an orthonormal basis, so this
needs no simplicity hypothesis and no eigenvalue gap. -/
private theorem inner_vEig_eq_zero {p : ℕ} (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    {k l : ℕ} (hk : k < Fintype.card (Fin p)) (hl : l < Fintype.card (Fin p)) (hkl : k ≠ l) :
    ⟪vEig A hA k, vEig A hA l⟫_ℝ = 0 := by
  have hkk : vEig A hA k = hA.eigenvectorBasis (eigIdx p ⟨k, hk⟩) := by
    rw [vEig, dif_pos hk]
  have hll : vEig A hA l = hA.eigenvectorBasis (eigIdx p ⟨l, hl⟩) := by
    rw [vEig, dif_pos hl]
  rw [hkk, hll]
  refine hA.eigenvectorBasis.orthonormal.2 fun h => hkl ?_
  exact congrArg Fin.val ((eigIdx p).injective h)

/-! ### 2. A table has at most `d N` spikes -/

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- Each spike direction of a table is a unit vector. This is `inner_col` at `k = l`, in the
shape the norm hypothesis of `align_inner_detR` takes. Mirror: `SpikedModel.hv`. -/
theorem norm_col (t : SpikedModelR μ n d rk) (N : ℕ) (k : Fin rk) : ‖t.col N k‖ = 1 := by
  have h : ‖t.col N k‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, t.inner_col N k k, if_pos rfl]
  rw [← Real.sqrt_sq (norm_nonneg (t.col N k)), h, Real.sqrt_one]

end SpikedModelR

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 3. The signed singular vector `v̂_ij` -/

/-- The eigenvalue index of the spike `(i, j)` is in range for `vEig`. -/
theorem idx_lt_card (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i)) (N : ℕ) :
    (j : ℕ) < Fintype.card (Fin (d N)) := by
  simpa using lt_of_lt_of_le j.isLt ((m.tbl i).rk_le_d N)

/-- `v̂_ij` is a unit vector. Mirror: `norm_vhat`. -/
theorem norm_vhatG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i)) (N : ℕ)
    (ω : Ω N) : ‖m.vhatG i j N ω‖ = 1 := by
  rw [UnalignedModelR.vhatG]
  split_ifs with h
  · exact norm_vEig _ _ (m.idx_lt_card i j N)
  · rw [norm_neg]
    exact norm_vEig _ _ (m.idx_lt_card i j N)

/-- `v̂_ij` lies in the eigenspace of `X_iᵀ X_i` at the sorted index `j`.
Mirror: `mem_topSpace_vhat`. -/
theorem mem_specSpace_vhatG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i))
    (N : ℕ) (ω : Ω N) :
    m.vhatG i j N ω ∈ specSpace (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      (eigSetIdx (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
        (isHermitian_transpose_mul_self ((m.tbl i).X N ω)) (j : ℕ)) := by
  rw [UnalignedModelR.vhatG]
  split_ifs with h
  · exact mem_specSpace_vEig _ _ (m.idx_lt_card i j N)
  · exact Submodule.neg_mem _ (mem_specSpace_vEig _ _ (m.idx_lt_card i j N))

/-- The sign convention of `vhatG`: `⟪v̂_ij, (V R_i)_j⟫ ≥ 0`. Mirror: `inner_vhat_nonneg`. -/
theorem inner_vhatG_nonneg (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i))
    (N : ℕ) (ω : Ω N) : 0 ≤ ⟪m.vhatG i j N ω, (m.tbl i).col N j⟫_ℝ := by
  rw [UnalignedModelR.vhatG]
  split_ifs with h
  · exact h
  · rw [inner_neg_left]
    linarith [not_le.mp h]

/-- The squared inner product with `v̂_ij` never exceeds the overlap at the index `j`. The sign
of `vhatG` is squared away, so this is `overlapIdx_ge_inner_sq` of `LinAlg/SpecIdx.lean` at the
signed vector. Mirror: the two uses of `overlap_ge_inner_sq` in `RankR/Defs.lean`. -/
theorem inner_vhatG_sq_le (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i))
    (N : ℕ) (ω : Ω N) (w : EuclideanSpace ℝ (Fin (d N))) :
    ⟪m.vhatG i j N ω, w⟫_ℝ ^ 2 ≤ overlapIdx ((m.tbl i).X N ω) (j : ℕ) w := by
  rw [UnalignedModelR.vhatG]
  split_ifs with h
  · exact overlapIdx_ge_inner_sq _ _ _
  · rw [inner_neg_left, neg_sq]
    exact overlapIdx_ge_inner_sq _ _ _

/-- On the almost sure event `TableLawR.simple`, the overlap at the index `j` is exactly the
squared signed inner product with `v̂_ij`. Mirror: the use of `overlap_eq_inner_sq` in
`align_inner`. -/
theorem overlapIdx_eq_inner_vhatG_sq (m : UnalignedModelR μ M n d r rk) (i : Fin M)
    (j : Fin (rk i)) (N : ℕ) (ω : Ω N)
    (hsimple : SimpleSpec (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      (isHermitian_transpose_mul_self ((m.tbl i).X N ω)) (rk i))
    (w : EuclideanSpace ℝ (Fin (d N))) :
    overlapIdx ((m.tbl i).X N ω) (j : ℕ) w = ⟪m.vhatG i j N ω, w⟫_ℝ ^ 2 := by
  rw [overlapIdx_eq_inner_sq _ hsimple j.isLt w, UnalignedModelR.vhatG]
  split_ifs with h
  · rfl
  · rw [inner_neg_left, neg_sq]
    rfl

/-- Two spikes of **one** table give orthogonal estimates: `⟪v̂_ij, v̂_ij'⟫ = 0` for `j ≠ j'`.
This is exact, not a limit. `v̂_ij` is `± vEig` at the sorted index `j`, and `vEig` reads an
orthonormal basis, so no simplicity and no eigenvalue gap enter. It is the probabilistic side
of `ABlock_intra` (`RankR/General.lean`). -/
theorem inner_vhatG_intra (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j j' : Fin (rk i))
    (hjj : j ≠ j') (N : ℕ) (ω : Ω N) :
    ⟪m.vhatG i j N ω, m.vhatG i j' N ω⟫_ℝ = 0 := by
  have hz := inner_vEig_eq_zero (m.tableGramG i N ω) (m.isHermitian_tableGramG i N ω)
    (m.idx_lt_card i j N) (m.idx_lt_card i j' N) (fun h => hjj (Fin.val_injective h))
  simp only [UnalignedModelR.vhatG]
  split_ifs <;> simp [hz]

/-- `vhatG` at a block index that is known as a pair. The pair equality is dependent, so
`simp` cannot rewrite it inside `vhatG`; `subst` on the whole pair can. -/
private theorem vhatG_congr_blk (m : UnalignedModelR μ M n d r rk)
    {P : (i : Fin M) × Fin (rk i)} {i : Fin M} {j : Fin (rk i)} (h : P = ⟨i, j⟩) (N : ℕ)
    (ω : Ω N) : m.vhatG P.1 P.2 N ω = m.vhatG i j N ω := by
  subst h
  rfl

/-- `inner_vhatG_intra` in flat-index form: two different flat indices of one table give
orthogonal estimates. `blk_flat` moves both indices to their block pair, which is what makes
the two `vhatG` read one table. -/
theorem inner_vhatG_blk_intra (m : UnalignedModelR μ M n d r rk) {p q : Fin (rtot rk)}
    (hpq : p ≠ q) (hsame : (blk p).1 = (blk q).1) (N : ℕ) (ω : Ω N) :
    ⟪m.vhatG (blk p).1 (blk p).2 N ω, m.vhatG (blk q).1 (blk q).2 N ω⟫_ℝ = 0 := by
  obtain ⟨i, j, rfl⟩ : ∃ i : Fin M, ∃ j : Fin (rk i), flat i j = p :=
    ⟨(blk p).1, (blk p).2, Equiv.apply_symm_apply finSigmaFinEquiv p⟩
  obtain ⟨i', j', rfl⟩ : ∃ i' : Fin M, ∃ j' : Fin (rk i'), flat i' j' = q :=
    ⟨(blk q).1, (blk q).2, Equiv.apply_symm_apply finSigmaFinEquiv q⟩
  rw [m.vhatG_congr_blk (blk_flat i j) N ω, m.vhatG_congr_blk (blk_flat i' j') N ω]
  simp only [blk_flat] at hsame
  subst hsame
  exact m.inner_vhatG_intra i j j' (fun h => hpq (congrArg (flat i) h)) N ω

/-! ### 4. The entries of `Ṽ V` and the shared frame `V` -/

/-- Entry `(p, k)` of `Ṽ V` is `⟪v̂_ij, V e_k⟫` at `(i, j) = blk p`. Mirror: `VtV_eq_inner`. -/
theorem VtVG_eq_inner (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (p : Fin (rtot rk)) (k : Fin r) :
    m.VtVG N ω p k = ⟪m.vhatG (blk p).1 (blk p).2 N ω, m.colVecG N k⟫_ℝ := by
  rw [UnalignedModelR.VtVG, Matrix.mul_apply, real_inner_eq_dotProduct]
  rfl

/-- Entry `(p, q)` of `Ṽ Ṽᵀ` is `⟪v̂_ij, v̂_{i'j'}⟫` at `(i, j) = blk p`, `(i', j') = blk q`.
Mirror: `gram_eq_inner` (`RankR/Defs.lean`). -/
theorem gramG_eq_inner (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (p q : Fin (rtot rk)) :
    m.gramG N ω p q
      = ⟪m.vhatG (blk p).1 (blk p).2 N ω, m.vhatG (blk q).1 (blk q).2 N ω⟫_ℝ := by
  rw [UnalignedModelR.gramG, Matrix.mul_apply, real_inner_eq_dotProduct]
  rfl

/-- `⟪(V R_i)_j, V e_k⟫ = (R_i)_{kj}`, the paper's `(V R_i)ᵀ V = R_iᵀ`
(`main_paper.tex:1964`). Mirror: `inner_v_colVec`. -/
theorem inner_col_colVecG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i))
    (k : Fin r) (N : ℕ) : ⟪(m.tbl i).col N j, m.colVecG N k⟫_ℝ = m.R i k j := by
  have hVV := m.hV N
  rw [real_inner_eq_dotProduct]
  change ∑ l, (m.tbl i).V N l j * (m.V N) l k = m.R i k j
  rw [m.hv i N]
  calc ∑ l, (m.V N * m.R i) l j * (m.V N) l k
      = ∑ l, ∑ t, m.R i t j * ((m.V N) l t * (m.V N) l k) := by
        refine Finset.sum_congr rfl fun l _ => ?_
        rw [Matrix.mul_apply, Finset.sum_mul]
        exact Finset.sum_congr rfl fun t _ => by ring
    _ = ∑ t, m.R i t j * ((m.V N)ᵀ * m.V N) t k := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun t _ => ?_
        rw [Matrix.mul_apply, Finset.mul_sum]
        exact Finset.sum_congr rfl fun l _ => by rw [Matrix.transpose_apply]
    _ = m.R i k j := by
        rw [hVV]
        simp [Matrix.one_apply]

/-- `⟪(V R_i)_k, (V R_{i'})_l⟫ = (R_iᵀ R_{i'})_{kl}`: the shared frame `V` cancels, so the
inner products of two signal frames read the paper's `R_iᵀ R_{i'}` (`main_paper.tex:1964`).
Mirror: `inner_v_v` of `RankR/Defs.lean`, which is this statement at `r_i = r_{i'} = 1`. -/
theorem inner_col_col (m : UnalignedModelR μ M n d r rk) (i i' : Fin M) (k : Fin (rk i))
    (l : Fin (rk i')) (N : ℕ) :
    ⟪(m.tbl i).col N k, (m.tbl i').col N l⟫_ℝ = ((m.R i)ᵀ * m.R i') k l := by
  have hmat : ((m.tbl i).V N)ᵀ * (m.tbl i').V N = (m.R i)ᵀ * m.R i' := by
    rw [m.hv i N, m.hv i' N, Matrix.transpose_mul, Matrix.mul_assoc,
      ← Matrix.mul_assoc ((m.V N)ᵀ), m.hV N, Matrix.one_mul]
  rw [real_inner_eq_dotProduct]
  change ∑ p, (m.tbl i).V N p k * (m.tbl i').V N p l = _
  rw [← hmat, Matrix.mul_apply]
  exact Finset.sum_congr rfl fun p _ => by rw [Matrix.transpose_apply]

/-! ### 5. `align`, `cross` and `delocUniform` in signed inner-product form -/

/-- `⟪v̂_ij, (V R_i)_j⟫ → β_ij`, the signed form of `TableLawR.align`. The sign convention of
`vhatG` turns the square root of the overlap into the inner product itself.
Mirror: `align_inner`. -/
theorem align_innerR (m : UnalignedModelR μ M n d r rk) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).TableLawR ci) (j : Fin (rk i)) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, (m.tbl i).col N j⟫_ℝ)
      (beta ((m.tbl i).θ j) ci) := by
  have h1 : TendstoInProb μ
      (fun N ω => Real.sqrt (overlapIdx ((m.tbl i).X N ω) (j : ℕ) ((m.tbl i).col N j)))
      (Real.sqrt (betaSq ((m.tbl i).θ j) ci)) :=
    (law.align j).comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
  refine h1.congr fun N => ?_
  filter_upwards [law.simple N] with ω hω
  rw [m.overlapIdx_eq_inner_vhatG_sq i j N ω hω, Real.sqrt_sq_eq_abs,
    abs_of_nonneg (m.inner_vhatG_nonneg i j N ω)]

/-- `⟪v̂_ij, (V R_i)_l⟫ → 0` for `l ≠ j`, the signed form of `TableLawR.cross`. No sign
convention is needed: the limit is `0`. -/
theorem cross_innerR (m : UnalignedModelR μ M n d r rk) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).TableLawR ci) (j l : Fin (rk i)) (hjl : j ≠ l) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, (m.tbl i).col N l⟫_ℝ) 0 := by
  have h1 : TendstoInProb μ
      (fun N ω => Real.sqrt (overlapIdx ((m.tbl i).X N ω) (j : ℕ) ((m.tbl i).col N l))) 0 := by
    have h := (law.cross j l hjl).comp_continuous (φ := Real.sqrt)
      Real.continuous_sqrt.continuousAt
    rwa [Real.sqrt_zero] at h
  refine TendstoInProb.of_le (fun N => ?_) h1
  filter_upwards with ω
  rw [sub_zero, ← Real.sqrt_sq_eq_abs]
  exact Real.sqrt_le_sqrt (m.inner_vhatG_sq_le i j N ω _)

/-- `⟪v̂_ij, q_N⟫ → 0` for a deterministic family `q` orthogonal to the **whole** signal frame
of table `i`, of norm at most one. The normalized `q_N` lies in `orthUnitR`, so
`TableLawR.delocUniform` bounds the probability; the case `q_N = 0` gives an empty event.
Mirror: `inner_vhat_perp_tendsto`. -/
theorem inner_vhatG_perp_tendsto (m : UnalignedModelR μ M n d r rk) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).TableLawR ci) (j : Fin (rk i))
    (q : (N : ℕ) → EuclideanSpace ℝ (Fin (d N)))
    (hq : ∀ (N : ℕ) (l : Fin (rk i)), ⟪q N, (m.tbl i).col N l⟫_ℝ = 0) (hn : ∀ N, ‖q N‖ ≤ 1) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, q N⟫_ℝ) 0 := by
  intro ε hε
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (law.delocUniform j (ε ^ 2) (by positivity)) (fun _ => zero_le) (fun N => ?_)
  by_cases h0 : q N = 0
  · have hempty : {ω : Ω N | ε ≤ |⟪m.vhatG i j N ω, q N⟫_ℝ - 0|} = ∅ := by
      ext ω
      simp only [h0, inner_zero_right, sub_zero, abs_zero, Set.mem_ofPred_eq,
        Set.mem_empty_iff_false, iff_false, not_le]
      exact hε
    rw [hempty]
    simp
  · have hqpos : 0 < ‖q N‖ := norm_pos_iff.mpr h0
    have hmem : (‖q N‖⁻¹ • q N) ∈ (m.tbl i).orthUnitR N := by
      refine ⟨?_, fun l => ?_⟩
      · rw [norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
        exact inv_mul_cancel₀ hqpos.ne'
      · rw [real_inner_smul_left, hq N l, mul_zero]
    refine le_trans (measure_mono ?_) (le_iSup₂ (f := fun w (_ : w ∈ (m.tbl i).orthUnitR N) =>
      μ N {ω | ε ^ 2 ≤ overlapIdx ((m.tbl i).X N ω) (j : ℕ) w}) _ hmem)
    intro ω hω
    have h1 : ε ≤ |⟪m.vhatG i j N ω, q N⟫_ℝ| := by
      have h := hω
      rw [Set.mem_ofPred_eq, sub_zero] at h
      exact h
    have hqe : q N = ‖q N‖ • (‖q N‖⁻¹ • q N) := by
      rw [smul_smul, mul_inv_cancel₀ hqpos.ne', one_smul]
    have hval : ⟪m.vhatG i j N ω, q N⟫_ℝ
        = ‖q N‖ * ⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ := by
      conv_lhs => rw [hqe]
      rw [real_inner_smul_right]
    rw [hval, abs_mul, abs_of_nonneg (norm_nonneg _)] at h1
    have h2 : ε ≤ |⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ| := by
      nlinarith [abs_nonneg (⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ), hn N, hqpos]
    have h3 : ⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ ^ 2
        ≤ overlapIdx ((m.tbl i).X N ω) (j : ℕ) (‖q N‖⁻¹ • q N) :=
      m.inner_vhatG_sq_le i j N ω _
    have h4 : ε ^ 2 ≤ overlapIdx ((m.tbl i).X N ω) (j : ℕ) (‖q N‖⁻¹ • q N) := by
      nlinarith [sq_abs (⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ),
        abs_nonneg (⟪m.vhatG i j N ω, ‖q N‖⁻¹ • q N⟫_ℝ)]
    exact h4

/-- `⟪v̂_ij, y_N⟫ → β_ij a_j` for a deterministic family `y` of norm at most one whose inner
products with the signal frame of table `i` are the constants `a_1, …, a_{r_i}`.

The split is the paper's `y = ∑_l a_l (V R_i)_l + perpSpan_i y` (`main_paper.tex:1961`) with
the whole frame in place of the single direction of the `r_i = 1` case: `align_innerR` gives
the term `l = j`, `cross_innerR` kills the terms `l ≠ j`, and `inner_vhatG_perp_tendsto` kills
the remainder. Mirror: `align_inner_det`. -/
theorem align_inner_detR (m : UnalignedModelR μ M n d r rk) {ci : ℝ} {i : Fin M}
    (law : (m.tbl i).TableLawR ci) (j : Fin (rk i))
    (y : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))) (a : Fin (rk i) → ℝ)
    (ha : ∀ (N : ℕ) (l : Fin (rk i)), ⟪(m.tbl i).col N l, y N⟫_ℝ = a l)
    (hy : ∀ N, ‖y N‖ ≤ 1) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, y N⟫_ℝ)
      (beta ((m.tbl i).θ j) ci * a j) := by
  classical
  have hA : TendstoInProb μ
      (fun N ω => ∑ l, a l * ⟪m.vhatG i j N ω, (m.tbl i).col N l⟫_ℝ)
      (∑ l, a l * (if l = j then beta ((m.tbl i).θ j) ci else 0)) := by
    refine TendstoInProb.finsum fun l => ?_
    rcases eq_or_ne l j with rfl | hlj
    · rw [if_pos rfl]
      exact TendstoInProb.const_mul (a l) (m.align_innerR law l)
    · rw [if_neg hlj, mul_zero]
      have h := TendstoInProb.const_mul (a l) (m.cross_innerR law j l (Ne.symm hlj))
      rwa [mul_zero] at h
  have hsumval : (∑ l, a l * (if l = j then beta ((m.tbl i).θ j) ci else 0))
      = beta ((m.tbl i).θ j) ci * a j := by
    rw [Finset.sum_eq_single j]
    · rw [if_pos rfl]
      ring
    · intro l _ hl
      rw [if_neg hl, mul_zero]
    · intro h
      exact absurd (Finset.mem_univ j) h
  have hperp : TendstoInProb μ
      (fun N ω => ⟪m.vhatG i j N ω, (m.tbl i).perpSpan N (y N)⟫_ℝ) 0 := by
    refine m.inner_vhatG_perp_tendsto law j (fun N => (m.tbl i).perpSpan N (y N))
      (fun N l => ?_)
      (fun N => le_trans (norm_perpFrame_le ((m.tbl i).inner_col N) (y N)) (hy N))
    rw [real_inner_comm]
    exact inner_v_perpFrame ((m.tbl i).inner_col N) l (y N)
  have hsum := hA.add hperp
  rw [add_zero, hsumval] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  rw [(m.tbl i).perpSpan_eq N (y N)]
  simp only [inner_sub_right, inner_sum, real_inner_smul_right, ha N]
  ring

/-! ### 6. The statements -/


/-- **`lem:general_rank_delocalization`, off-diagonal half, at general `r_i`**
(`main_paper.tex:1948`): for two spikes in **different** tables,
`⟪v̂_ij, v̂_{i'j'}⟫ → β_ij β_{i'j'} ⟪(R_i)_j, (R_{i'})_{j'}⟫`. Inside one table the inner
product is exactly `0` for `j ≠ j'` and exactly `1` for `j = j'`, because the rows of `Ṽ`
coming from one table are orthonormal singular vectors; that half needs no limit and no
hypothesis, and `ABlock_intra`, `ABlock_diag` give the matching entries of `A_{β,R}`.

The route is the paper's own split
`v̂_ijᵀ v̂_{i'j'} = v̂_ijᵀ (V R_i)_j (V R_i)_jᵀ v̂_{i'j'}
+ v̂_ijᵀ (I - (V R_i)_j (V R_i)_jᵀ) v̂_{i'j'}`
(`main_paper.tex:1967`): the first term is `align_inner` of the two tables against a
deterministic direction, and the second needs `TableLawR.delocUniform` of table `i` at the
**random** direction of table `i'`, hence Fubini against the product law and `hG`. -/
theorem lem_general_rank_delocalization_general (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    {i i' : Fin M} (hii : i ≠ i') (j : Fin (rk i)) (j' : Fin (rk i')) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, m.vhatG i' j' N ω⟫_ℝ)
      (beta ((m.tbl i).θ j) (c i) * beta ((m.tbl i').θ j') (c i') *
        ((m.R i)ᵀ * m.R i') j j') := by
  classical
  obtain ⟨ν, hν, hind⟩ := hG
  have : ∀ (N : ℕ) (k : Fin M), IsProbabilityMeasure (ν N k) := hν
  -- the two Gram entries of the alignment matrices are transposes of each other
  have hRt : ((m.R i')ᵀ * m.R i) j' j = ((m.R i)ᵀ * m.R i') j j' := by
    rw [Matrix.mul_apply, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun t _ => ?_
    rw [Matrix.transpose_apply, Matrix.transpose_apply]
    ring
  -- the pointwise split against the whole signal frame of table `i`
  have hsplit : ∀ (N : ℕ) (ω : Ω N),
      ⟪m.vhatG i j N ω, m.vhatG i' j' N ω⟫_ℝ
        = (∑ l, ⟪m.vhatG i' j' N ω, (m.tbl i).col N l⟫_ℝ *
            ⟪m.vhatG i j N ω, (m.tbl i).col N l⟫_ℝ)
          + ⟪m.vhatG i j N ω, (m.tbl i).perpSpan N (m.vhatG i' j' N ω)⟫_ℝ := by
    intro N ω
    have hc : ∀ l, ⟪(m.tbl i).col N l, m.vhatG i' j' N ω⟫_ℝ
        = ⟪m.vhatG i' j' N ω, (m.tbl i).col N l⟫_ℝ := fun l => real_inner_comm _ _
    rw [(m.tbl i).perpSpan_eq N (m.vhatG i' j' N ω)]
    simp only [inner_sub_right, inner_sum, real_inner_smul_right, hc]
    ring
  -- signal part: `⟪v̂_{i'j'}, (V R_i)_l⟫` against a deterministic direction of table `i`
  have hcoef : ∀ l : Fin (rk i), TendstoInProb μ
      (fun N ω => ⟪m.vhatG i' j' N ω, (m.tbl i).col N l⟫_ℝ)
      (beta ((m.tbl i').θ j') (c i') * ((m.R i')ᵀ * m.R i) j' l) := fun l =>
    m.align_inner_detR (law i') j' (fun N => (m.tbl i).col N l)
      (fun t => ((m.R i')ᵀ * m.R i) t l) (fun N t => m.inner_col_col i' i t l N)
      (fun N => le_of_eq ((m.tbl i).norm_col N l))
  have hsig : TendstoInProb μ
      (fun N ω => ∑ l, ⟪m.vhatG i' j' N ω, (m.tbl i).col N l⟫_ℝ *
        ⟪m.vhatG i j N ω, (m.tbl i).col N l⟫_ℝ)
      (∑ l, beta ((m.tbl i').θ j') (c i') * ((m.R i')ᵀ * m.R i) j' l *
        (if l = j then beta ((m.tbl i).θ j) (c i) else 0)) := by
    refine TendstoInProb.finsum fun l => ?_
    refine (hcoef l).mul ?_
    rcases eq_or_ne l j with hlj | hlj
    · rw [if_pos hlj, hlj]
      exact m.align_innerR (law i) j
    · rw [if_neg hlj]
      exact m.cross_innerR (law i) j l hlj.symm
  have hsumval : (∑ l, beta ((m.tbl i').θ j') (c i') * ((m.R i')ᵀ * m.R i) j' l *
        (if l = j then beta ((m.tbl i).θ j) (c i) else 0))
      = beta ((m.tbl i).θ j) (c i) * beta ((m.tbl i').θ j') (c i') *
        ((m.R i)ᵀ * m.R i') j j' := by
    rw [Finset.sum_eq_single j]
    · rw [if_pos rfl, hRt]
      ring
    · intro l _ hl
      rw [if_neg hl, mul_zero]
    · intro h
      exact absurd (Finset.mem_univ j) h
  -- perpendicular part: `delocUniform` of table `i` at the random direction of table `i'`
  have hperp : TendstoInProb μ
      (fun N ω => ⟪m.vhatG i j N ω, (m.tbl i).perpSpan N (m.vhatG i' j' N ω)⟫_ℝ) 0 := by
    refine TendstoInProb.of_le
      (g := fun N ω => Real.sqrt (overlapIdx ((m.tbl i).X N ω) (j : ℕ)
        (delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ)))) (fun N => ?_) ?_
    · filter_upwards [(law i').simple N] with ω hω
      rw [sub_zero]
      have h1 : |⟪m.vhatG i j N ω, (m.tbl i).perpSpan N (m.vhatG i' j' N ω)⟫_ℝ|
          ≤ |⟪m.vhatG i j N ω,
              delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ)⟫_ℝ| :=
        abs_inner_perpFrame_le ((m.tbl i).inner_col N) _ hω j'.isLt
          (m.mem_specSpace_vhatG i' j' N ω) (m.norm_vhatG i' j' N ω) _
      have h2 : |⟪m.vhatG i j N ω,
            delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ)⟫_ℝ|
          ≤ Real.sqrt (overlapIdx ((m.tbl i).X N ω) (j : ℕ)
            (delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ))) := by
        rw [← Real.sqrt_sq_eq_abs]
        exact Real.sqrt_le_sqrt (m.inner_vhatG_sq_le i j N ω _)
      linarith
    · intro ε hε
      have hK := (law i).delocUniform j (ε ^ 2) (by positivity)
      refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hK
        (fun _ => zero_le) (fun N => ?_)
      have hset : {ω : Ω N | ε ≤ |Real.sqrt (overlapIdx ((m.tbl i).X N ω) (j : ℕ)
            (delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ))) - 0|}
          = {ω : Ω N | ε ^ 2 ≤ overlapIdx ((m.tbl i).X N ω) (j : ℕ)
            (delocDirIdx ((m.tbl i).col N) ((m.tbl i').X N ω) (j' : ℕ))} := by
        ext ω
        simp only [sub_zero, Set.mem_ofPred_eq, abs_of_nonneg (Real.sqrt_nonneg _)]
        exact Real.le_sqrt hε.le (overlapIdx_nonneg _ _ _)
      rw [hset]
      exact measure_deloc_le_of_piR m.tbl hind hii (j : ℕ) (j' : ℕ) N (by positivity)
  have hsum := hsig.add hperp
  rw [add_zero, hsumval] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  exact (hsplit N ω).symm

/-- **`lem:general_rank_delocalization`, second display, at general `r_i`**:
`Ṽ Ṽᵀ → A_{β,R}` entrywise and signed. The diagonal and the within-table entries are exact,
not limits; the cross-table entries are `lem_general_rank_delocalization_general`. -/
theorem gramR_general (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (p q : Fin (rtot rk)) :
    TendstoInProb μ (fun N ω => m.gramG N ω p q) (ABlock β m.R p q) := by
  classical
  rcases eq_or_ne p q with rfl | hpq
  · -- the diagonal: exact, because `v̂_ij` is a unit vector
    rw [ABlock_diag (β := β) m.hR p]
    refine (TendstoInProb.const μ 1).congr fun N => ?_
    filter_upwards with ω
    rw [m.gramG_eq_inner N ω p p, real_inner_self_eq_norm_sq,
      m.norm_vhatG (blk p).1 (blk p).2 N ω, one_pow]
  rcases eq_or_ne (blk p).1 (blk q).1 with hsame | hdiff
  · -- two spikes of one table: exact `0` on both sides
    rw [ABlock_intra (β := β) m.hR hpq hsame]
    refine (TendstoInProb.const μ 0).congr fun N => ?_
    filter_upwards with ω
    rw [m.gramG_eq_inner N ω p q, m.inner_vhatG_blk_intra hpq hsame N ω]
  · -- two tables: the delocalization limit
    have hA : ABlock β m.R p q
        = beta ((m.tbl (blk p).1).θ (blk p).2) (c (blk p).1) *
            beta ((m.tbl (blk q).1).θ (blk q).2) (c (blk q).1) *
            ((m.R (blk p).1)ᵀ * m.R (blk q).1) (blk p).2 (blk q).2 := by
      rw [ABlock_apply, if_neg hpq, add_zero]
      simp only [betaFlat, hβdef]
    rw [hA]
    refine (m.lem_general_rank_delocalization_general c law hG hdiff (blk p).2 (blk q).2).congr
      fun N => ?_
    filter_upwards with ω
    exact (m.gramG_eq_inner N ω p q).symm

/-- **`lem:general_rank_delocalization`, first display, at general `r_i`**:
`Ṽ V → B_R` entrywise and signed, `⟪v̂_ij, V e_k⟫ → β_ij (R_i)_{jk}`. The column `V e_k`
decomposes in the orthonormal family `(V R_i)_1, …, (V R_i)_{r_i}` plus a remainder orthogonal
to the signal span of table `i`; `TableLawR.align` gives the `j`-th coefficient,
`TableLawR.cross` kills the others and `TableLawR.delocUniform` kills the remainder. This needs
no independence across tables, so `hG` does not appear. -/
theorem VtV_tendsto_general (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (p : Fin (rtot rk)) (k : Fin r) :
    TendstoInProb μ (fun N ω => m.VtVG N ω p k) (BBlock β m.R p k) := by
  have hB : BBlock β m.R p k
      = beta ((m.tbl (blk p).1).θ (blk p).2) (c (blk p).1) * m.R (blk p).1 k (blk p).2 := by
    rw [BBlock, Matrix.of_apply, betaFlat, hβdef]
  rw [hB]
  have h := m.align_inner_detR (law (blk p).1) (blk p).2 (fun N => m.colVecG N k)
    (fun l => m.R (blk p).1 k l) (fun N l => m.inner_col_colVecG (blk p).1 l k N)
    (fun N => le_of_eq (m.norm_colVecG N k))
  refine h.congr fun N => ?_
  filter_upwards with ω
  exact (m.VtVG_eq_inner N ω p k).symm

/-- **`prop:general_rank_unweighted_svdstack`** (`main_paper.tex:799`) at general `r_i`, Layer 1
form. `hgap` is the paper's `λ_r(A_{β,R}) - λ_{r+1}(A_{β,R}) > 0` with its convention
`λ_{r̃+1} := -∞`, which `TopGap` encodes by quantifying over the indices; `hc` gives
`0 ≤ β_ij < 1`, hence `A_{β,R} ⪰ D ≻ 0` (`ABlock_posDef`). The paper's hypothesis that each
`Θ_i` has `r_i` distinct entries sits inside `TableLawR` through the model field
`SpikedModelR.hθanti`.

Route, verbatim from `prop_general_rank_unweighted_svdstack` of `RankR/Defs.lean` with `Fin M`
replaced by `Fin r̃`: `gramR_general` and `VtV_tendsto_general` give the two entrywise limits
jointly on `(Fin r̃ × Fin r̃) ⊕ (Fin r̃ × Fin r)`, `perfRG` is `traceFun r` of that family at
every `ω`, and `continuousAt_traceFun` with `TendstoInProbPi.comp_continuous` transports the
limit. No good event is needed. -/
theorem prop_general_rank_unweighted_svdstack_general (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hrr : r ≤ rtot rk)
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRG N ω) (limitRG β m.R) := by
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · -- `r = 0`: `perfRG` and `limitRG` are traces of `0 × 0` matrices, so both are `0`
    have hzero : ∀ (N : ℕ) (ω : Ω N), m.perfRG N ω = 0 := fun N ω => by
      simp [UnalignedModelR.perfRG, Matrix.trace]
    have hlim : limitRG β m.R = 0 := by simp [limitRG, Matrix.trace]
    rw [hlim]
    exact (TendstoInProb.const μ 0).congr fun N =>
      Filter.Eventually.of_forall fun ω => (hzero N ω).symm
  have hb01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have hA : (ABlock β m.R).IsHermitian := isHermitian_ABlock β m.R
  have hpd : (ABlock β m.R).PosDef :=
    ABlock_posDef m.R (fun i j => (hb01 i j).1) (fun i j => (hb01 i j).2)
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
  have hrm : r - 1 < Fintype.card (Fin (rtot rk)) := by
    simp only [Fintype.card_fin]
    omega
  have hpos : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ := by
    set t : Fin (rtot rk) := Fintype.equivOfCardEq (Fintype.card_fin _) ⟨r - 1, hrm⟩ with ht
    have hev : hA.eigenvalues t = hA.eigenvalues₀ ⟨r - 1, hrm⟩ := by
      rw [Matrix.IsHermitian.eigenvalues, ht, Equiv.symm_apply_apply]
    rw [← hev]
    exact hpd.eigenvalues_pos t
  have hconv : TendstoInProbPi μ
      (fun N ω => Sum.elim
        (fun t : Fin (rtot rk) × Fin (rtot rk) => m.gramG N ω t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => m.VtVG N ω t.1 t.2))
      (Sum.elim (fun t : Fin (rtot rk) × Fin (rtot rk) => ABlock β m.R t.1 t.2)
        (fun t : Fin (rtot rk) × Fin r => BBlock β m.R t.1 t.2)) := by
    rintro (t | t)
    · exact m.gramR_general c β hβdef law hG t.1 t.2
    · exact m.VtV_tendsto_general c β hβdef law t.1 t.2
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hr hrp hgap hpos (BBlock β m.R)) hconv
  rw [traceFun_eq hA (BBlock β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gramG N ω) (m.VtVG N ω)

/-! ### 7. The weighted estimator: the entrywise rewrites and one admissible weight matrix

`RankR/Flatten.lean` (task B0) shows that every block object of `RankR/General.lean` is the
`r_i = 1` object of `RankR/Weighted.lean` read at table count `r̃ = rtot rk`, alignment family
`Rcol R` and scalar family `betaFlat β`, and each bridge is `rfl`. So the deterministic half
of `thm:gen_rank_weight_svdstak` (Ky Fan, the eigengap at `W⋆`, the attained value `L⋆`) is
not proved again here: `limitRW_le_opt`, `limitRW_optWR`, `topGap_optWR_of_rankBR` and
`abetaRW_optWR_posDef` apply as they stand. Only the probabilistic half is new, because
`perfRGW` reads `UnalignedModelR` and not `UnalignedModel`; it is
`thm_gen_rank_weight_svdstak_general_r_conv`, the transcription of
`thm_gen_rank_weight_svdstak_general` (`RankR/Weighted.lean`) with `gramR_general` and
`VtV_tendsto_general` in place of `gramR` and `VtV_tendsto`. -/

/-- `Ṽ_W Ṽ_Wᵀ = W (Ṽ Ṽᵀ) Wᵀ`, so the weighted Gram matrix is a fixed linear image of the
unweighted one. Mirror: `UnalignedModel.gramW_eq` (`RankR/Weighted.lean`). -/
theorem gramWG_eq (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    m.gramWG W N ω = W * m.gramG N ω * Wᵀ := by
  rw [UnalignedModelR.gramWG, UnalignedModelR.VtWG, UnalignedModelR.gramG,
    Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (m.VtG N ω), ← Matrix.mul_assoc]

/-- `Ṽ_W V = W (Ṽ V)`. Mirror: `UnalignedModel.VtVW_eq` (`RankR/Weighted.lean`). -/
theorem VtVWG_eq (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    m.VtVWG W N ω = W * m.VtVG N ω := by
  rw [UnalignedModelR.VtVWG, UnalignedModelR.VtWG, UnalignedModelR.VtVG, Matrix.mul_assoc]

/-- **`thm:gen_rank_weight_svdstak`** (`main_paper.tex:893`) for one admissible weight matrix
`W`, at general `r_i`. `hgapW` is the paper's `eq:weighted_eigengap` (`main_paper.tex:2010`)
and `hposW` its `λ_r(W A_{β,R} Wᵀ) > 0`, both read on the weighted limit matrix.

Route, verbatim from `UnalignedModel.thm_gen_rank_weight_svdstak_general`
(`RankR/Weighted.lean`) with `Fin M` replaced by `Fin r̃`: `gramWG_eq` and `VtVWG_eq` write the
two random matrices as fixed real linear combinations of the entries of `Ṽ Ṽᵀ` and `Ṽ V`,
`gramR_general` and `VtV_tendsto_general` give those entrywise limits, and
`continuousAt_traceFun` with `TendstoInProbPi.comp_continuous` transports the limit twice. -/
theorem thm_gen_rank_weight_svdstak_general_r_conv (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hr : 0 < r) (hrr : r ≤ rtot rk)
    (hgapW : TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r)
    (hposW : 0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) := by
  have hA : (ABlockW W β m.R).IsHermitian := isHermitian_ABlockW W β m.R
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
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
    · obtain ⟨i, k⟩ := t
      have hcont : Continuous fun z :
          (Fin (rtot rk) × Fin (rtot rk)) ⊕ (Fin (rtot rk) × Fin r) → ℝ =>
            ∑ a : Fin (rtot rk), W i a * z (Sum.inr (a, k)) :=
        continuous_finsetSum _ fun a _ => continuous_const.mul (continuous_apply _)
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hconv0
      simp only [Sum.elim_inr] at h
      have hlim : BBlockW W β m.R i k = ∑ a : Fin (rtot rk), W i a * BBlock β m.R a k :=
        Matrix.mul_apply
      rw [Sum.elim_inr, hlim]
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      simp only [Sum.elim_inr]
      rw [m.VtVWG_eq]
      exact (Matrix.mul_apply).symm
  have hall := TendstoInProbPi.comp_continuous
    (continuousAt_traceFun hA hr hrp hgapW hposW (BBlockW W β m.R)) hconv
  rw [traceFun_eq hA (BBlockW W β m.R)] at hall
  refine hall.congr fun N => ?_
  filter_upwards with ω
  exact traceFun_eq (m.isHermitian_gramWG W N ω) (m.VtVWG W N ω)

/-- **`thm:gen_rank_weight_svdstak`** (`main_paper.tex:893`) at general `r_i`, both halves in
one declaration, as `thm_gen_rank_weight_svdstak_max` states them at `r_i = 1`.

1. With `W⋆ = D^{-1/2}` the performance of weighted svdstack tends to
   `L⋆ = r - ∑_{ℓ=1}^r λ_{r̃+1-ℓ}(A_{β,R}^{-1/2} D A_{β,R}^{-1/2})`.
2. For every admissible weight matrix `W ∈ ℝ^{r̃ × r̃}` the performance at `W` tends to
   `limitRGW W`, which is at most `L⋆`.

`W` is a full matrix, not a block-diagonal one: the paper allows every `W ∈ ℝ^{r̃ × r̃}`
(`main_paper.tex:886`) and the point of the theorem is that the optimum is still the diagonal
`W⋆`, so restricting `W` would weaken the claim. Admissibility is the paper's own eigengap
condition `eq:weighted_eigengap` (`main_paper.tex:2010`) plus the positivity of
`λ_{r-1}(W A_{β,R} Wᵀ)`, which the paper needs for `Λ_r^{-1/2}` to exist; the footnote at
`main_paper.tex:890` excludes every other `W`.

The hypothesis is `rank B_R = r`, not the paper's `β_ij > 0` plus `Rank(∑ R_i R_iᵀ) = r`:
decision D16 of `notes/FLAGGED.md` shows the paper's remark that a component with `β_ij = 0`
can be removed is false in general, and `rank B_R = r` is what the eigengap at `W⋆` needs.

Route: the flatten bridges of `RankR/Flatten.lean` read every block object as the `r_i = 1`
object of `RankR/Weighted.lean` at table count `r̃`, so the deterministic half is
`topGap_optWR_of_rankBR`, `abetaRW_optWR_posDef`, `limitRW_optWR` and `limitRW_le_opt` as they
stand. The convergence half is `thm_gen_rank_weight_svdstak_general_r_conv`, instantiated at
`W⋆` for the first conjunct and at the given `W` for the second. -/
theorem thm_gen_rank_weight_svdstak_general_r (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hr : 0 < r) (hrr : r ≤ rtot rk)
    (hrankB : (BBlock β m.R).rank = r)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (by simpa using hrr)) ∧
      ∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (by simpa using hrr) := by
  have hrp : r ≤ Fintype.card (Fin (rtot rk)) := by simpa using hrr
  have hβ01 : ∀ i j, 0 ≤ β i j ∧ β i j < 1 := fun i j => by
    rw [hβdef i j]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ p, 0 ≤ betaFlat β p := fun p => (hβ01 (blk p).1 (blk p).2).1
  have h1 : ∀ p, betaFlat β p < 1 := fun p => (hβ01 (blk p).1 (blk p).2).2
  have hrankB' : (BR (betaFlat β) (Rcol m.R)).rank = r := by
    rw [← BBlock_eq_BR]
    exact hrankB
  -- the eigengap and the positivity at `W⋆`, from `rank B_R = r`
  have hgapOpt : TopGap (ABlockW (optWG β) β m.R) (isHermitian_ABlockW (optWG β) β m.R) r :=
    topGap_optWR_of_rankBR (β := betaFlat β) (Rcol m.R) h0 h1 hrankB'
  have hpdOpt : (ABlockW (optWG β) β m.R).PosDef :=
    abetaRW_optWR_posDef (β := betaFlat β) (Rcol m.R) h0 h1
  have hposOpt : 0 < (isHermitian_ABlockW (optWG β) β m.R).eigenvalues₀
      ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ := by
    rw [← eigenvalues_eigIdx]
    exact hpdOpt.eigenvalues_pos _
  refine ⟨?_, fun W hgapW hposW => ⟨?_, ?_⟩⟩
  · -- first conjunct: the general limit at `W⋆`, then the attained value `L⋆`
    have h := m.thm_gen_rank_weight_svdstak_general_r_conv c β (optWG β) hβdef hr hrr
      hgapOpt hposOpt law hG
    have hval : limitRGW (optWG β) β m.R = limitOptG β m.R hrp := by
      rw [limitRGW_eq_limitRW, limitOptG_eq_limitROpt, optWG_eq_optWR,
        limitRW_optWR (betaFlat β) (Rcol m.R) h0 h1 hrp]
    rwa [hval] at h
  · exact m.thm_gen_rank_weight_svdstak_general_r_conv c β W hβdef hr hrr hgapW hposW law hG
  · -- Ky Fan: no admissible weight matrix beats `L⋆`
    rw [limitRGW_eq_limitRW, limitOptG_eq_limitROpt]
    exact limitRW_le_opt W (betaFlat β) (Rcol m.R) h0 h1 hr hrp hgapW hposW


end UnalignedModelR


/-- **The reduction of the black box at `r_i = 1`.** A rank-one `SpikedModel` and a
`SpikedModelR` with one spike that carry the same data matrix, the same spike direction and the
same signal strength satisfy the same law: `SingleTableLaw c` gives `TableLawR c`. So the
`r_i = 1` slice of this file needs no new hypothesis, exactly as `RankR/Defs.lean` claims.

The lemma is stated on two models tied by `hX`, `hv`, `hθ` rather than on a constructed
`SpikedModel.toRankR`, because until F8 (2026-09-05) `SpikedModelR` required `0 < θ` while
`SpikedModel` allows `θ = 0`, so no total map existed. With `hθnn` the map is total; see
`SpikedModel.toRankR` and `tableLawR_of_singleTableLaw'` below.

Field by field: `align` is `SingleTableLaw.align` after `overlapIdx_zero`; `cross` is vacuous
on `Fin 1`; `delocUniform` is `SingleTableLaw.delocUniform` after `orthUnitR = orthUnit`;
`simple` is `topSimple` after `SimpleSpec _ _ 1 ↔ TopSimple`. -/
theorem tableLawR_of_singleTableLaw {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ} (m : SpikedModelR μ n d 1) (m' : SpikedModel μ n d)
    (hX : ∀ (N : ℕ) (ω : Ω N), m.X N ω = m'.X N ω) (hv : ∀ N, m'.v N = m.col N 0)
    (hθ : m'.θ = m.θ 0) {c : ℝ} (law : m'.SingleTableLaw c) : m.TableLawR c := by
  have hval0 : ((0 : Fin 1) : ℕ) = 0 := rfl
  -- `m.orthUnitR N` and `m'.orthUnit N` are the same set once `m.col N 0 = m'.v N` (`hv`), since
  -- `Fin 1` has only the index `0`.
  have horth : ∀ N, m.orthUnitR N = m'.orthUnit N := by
    intro N
    ext w
    simp only [SpikedModelR.orthUnitR, SpikedModel.orthUnit, Set.mem_ofPred_eq]
    constructor
    · rintro ⟨hw1, hw2⟩
      exact ⟨hw1, by rw [hv N]; exact hw2 0⟩
    · rintro ⟨hw1, hw2⟩
      refine ⟨hw1, fun k => ?_⟩
      obtain rfl := Fin.eq_zero k
      rw [← hv N]
      exact hw2
  refine
    { align := ?_
      cross := ?_
      delocUniform := ?_
      simple := ?_ }
  · -- `align`: `overlapIdx_zero` collapses the index-`0` overlap to `overlap`, then transport
    -- through `hX`, `hv`, `hθ` and apply `law.align`.
    intro k
    obtain rfl := Fin.eq_zero k
    have heq : (fun N ω => overlapIdx (m.X N ω) ((0 : Fin 1) : ℕ) (m.col N 0))
        = fun N ω => overlap (m'.X N ω) (m'.v N) := by
      funext N ω
      rw [hval0, overlapIdx_zero (m.hd N), hX N ω, ← hv N]
    rw [heq, ← hθ]
    exact law.align
  · -- `cross`: `k ≠ l` is impossible on `Fin 1`.
    intro k l hkl
    exact absurd (Subsingleton.elim k l) hkl
  · -- `delocUniform`: rewrite the test-direction set by `horth`, then the body by
    -- `overlapIdx_zero` and `hX`, and apply `law.delocUniform`.
    intro k ε hε
    obtain rfl := Fin.eq_zero k
    have heq : (fun N => ⨆ w ∈ m.orthUnitR N,
          μ N {ω | ε ≤ overlapIdx (m.X N ω) ((0 : Fin 1) : ℕ) w})
        = fun N => ⨆ w ∈ m'.orthUnit N, μ N {ω | ε ≤ overlap (m'.X N ω) w} := by
      funext N
      rw [horth N]
      simp_rw [hval0, overlapIdx_zero (m.hd N), hX N]
    rw [heq]
    exact law.delocUniform ε hε
  · -- `simple`: `simpleSpec_one_of_topSimple` on the a.s. event `law.topSimple`, transported
    -- through `hX`.
    intro N
    filter_upwards [law.topSimple N] with ω hω
    rw [hX N ω]
    exact simpleSpec_one_of_topSimple _ _ hω


namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ}

/-- The bridge on the constructed model: `SingleTableLaw c` gives `TableLawR c` for
`m.toRankR`. The `r_i = 1` payoff of F8 (2026-09-05): a `SpikedModelR` at one spike needs no
law of its own; `SingleTableLaw` (proved for Gaussian noise in `RMT/Full.lean`) is enough. -/
theorem tableLawR_of_singleTableLaw' (m : SpikedModel μ n d) {c : ℝ} (law : m.SingleTableLaw c) :
    m.toRankR.TableLawR c :=
  tableLawR_of_singleTableLaw m.toRankR m (fun N ω => m.toRankR_X N ω)
    (fun N => (m.toRankR_col N).symm) rfl law

end SpikedModel

end StackedSVD
