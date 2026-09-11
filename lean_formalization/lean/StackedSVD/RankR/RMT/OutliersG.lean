/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Outliers

/-!
# Task U7, gap G2: the deterministic core with a sub-frame

Gap G2 of `notes/archive/rankr_plan_A.md` section 4. `OutliersR.align_det`
(`RankR/RMT/Outliers.lean`) is the rank-1 mirror of this file: it uses one
`Q : Matrix (Fin p) (Fin r) ℝ` for the split `S = W + Q Qᵀ` and for the frame, and it
derives the count `= r` inside the proof through
`Frame.specProjTop_frame_approx_of_split`. In the mixed case only `s ≤ r` spikes are
supercritical. The frame then has `s` vectors while the split still has `r` columns, and
the count `= s` comes from outside, from the edge theorem of task U7a
(`RankR/RMT/EdgeR.lean`, `RankRStack.tendsto_measure_eigenvalues₀_le`).

The sub-frame is indexed by an injection `f : Fin s → Fin r`. Its columns are the columns
of `Q.submatrix id f`, the same `submatrix id f` shape that `RankRStack.qsubR` uses, so one
`f` feeds both the edge theorem and this file.

Contents:

1. `resCG`, the residual constant. `align_det` reads the column norms
   `‖q_l‖² ≤ 2 ρ_l` off the diagonal forms `hE1 l l`. A subcritical column has no `ρ_l`, so
   here the column norms enter as the hypothesis `hcol` with one constant `Cq`, and
   `resCG Cq r ν = √(Cq r r) ∑ (√ν)⁻¹` replaces `resC ρ ν`. The residual identity
   `residual_eq` sums over all `r` columns of `Q`, which is why `Cq` bounds all of them.
2. `residual_bound_sub` and `gram_bound_sub`, the two frame estimates. They are public so
   that a consumer builds them once and feeds both `count_eq_of_frame_of_edge` and
   `align_detG`.
3. `align_detG`, the sub-frame core. It takes the count as the hypothesis `hI` and calls
   `Frame.specProj_frame_approx` directly, so `hrp`, `hr`, `hpsd` and `hδm` of `align_det`
   disappear. The target of the overlaps is a general `t : Fin s → ℝ` with `|t k| ≤ Lb`,
   and the conclusion is `∑ k, t k ^ 2 / ν k`. The two corollaries `align_detG_mem` and
   `align_detG_not_mem` specialize `t` to one spike and to zero.
4. `count_eq_of_frame_of_edge`, the count itself: the upper bound from the edge hypothesis
   in the sorted index, the lower bound from `Frame.le_card_filter_of_frame`.
   `count_eq_of_forms_of_edge` is the same statement with the entrywise accuracies as input.
   `eigenvalues₀_le_of_edge_le` weakens the edge hypothesis from `s` to `r`, the shape that
   `Frame.norm_sq_specProjTop_split` needs.

Ties among the `ρ_k` stay allowed: `gramC` handles `ρ k = ρ l`. No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace OutliersR

open R4

/-! ### 1. The residual constant of the sub-frame version -/

/-- Residual constant of the sub-frame version. `Cq` bounds every column norm
`∑ i, Q i l ^ 2` of the split matrix (`r` columns), while `ν` lives on the frame `Fin s`.
The rank-1 mirror is `resC`. -/
noncomputable def resCG {s : ℕ} (Cq : ℝ) (r : ℕ) (ν : Fin s → ℝ) : ℝ :=
  Real.sqrt (Cq * r * r) * ∑ k, (Real.sqrt (ν k))⁻¹

theorem resCG_nonneg {s : ℕ} (Cq : ℝ) (r : ℕ) (ν : Fin s → ℝ) : 0 ≤ resCG Cq r ν := by
  refine mul_nonneg (Real.sqrt_nonneg _) (Finset.sum_nonneg fun k _ => ?_)
  positivity

/-! ### 2. The two frame estimates -/

section Core

variable {p r s : ℕ}

/-- The columns of the sub-frame are the `f`-columns of `Q`. -/
theorem submatrix_col (Q : Matrix (Fin p) (Fin r) ℝ) (f : Fin s → Fin r) (k : Fin s) :
    (fun i => (Q.submatrix id f) i k) = fun i => Q i (f k) := rfl

/-- The residual bound of the sub-frame version (the mirror of step 2 of `align_det`).
The residual identity sums over all `r` columns of `Q`, so `hcol` bounds all of them, and
the cross entries `hE1 l k` with `l` outside the range of `f` are part of the sum. -/
theorem residual_bound_sub {W S : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {f : Fin s → Fin r} {ρ ν : Fin s → ℝ} (hν : ∀ k, 0 < ν k)
    (hz : ∀ k, lamMax W hW < ρ k) {Cq η : ℝ} (hCq : 0 ≤ Cq) (hη0 : 0 ≤ η)
    (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η) (k : Fin s) :
    ‖toOp S (yhatv W (Q.submatrix id f) ρ ν k)
        - ρ k • yhatv W (Q.submatrix id f) ρ ν k‖ ≤ resCG Cq r ν * η := by
  classical
  have hsfmul : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν k))⁻¹ = (ν k)⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (hν k).le]
  have hre : S *ᵥ xtv W (Q.submatrix id f) ρ k - ρ k • xtv W (Q.submatrix id f) ρ k
      = Q *ᵥ fun l => (if l = f k then (1 : ℝ) else 0)
          + cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k)) :=
    residual_eq hW hSeq (hz k) (f k)
  have hsv : ∀ l : Fin r, |(if l = f k then (1 : ℝ) else 0)
      + cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))| ≤ η := by
    intro l
    have hb := hE1 l k
    rcases eq_or_ne l (f k) with rfl | hlk
    · rw [if_pos rfl] at hb ⊢
      calc |(1 : ℝ) + cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f k))|
          = |cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f k)) - -1| := by
            rw [sub_neg_eq_add, add_comm]
        _ ≤ η := hb
    · rw [if_neg hlk, sub_zero] at hb
      rw [if_neg hlk, zero_add]
      exact hb
  have hsvsum : ∑ l, ((if l = f k then (1 : ℝ) else 0)
      + cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))) ^ 2 ≤ (r : ℝ) * η ^ 2 := by
    calc ∑ l, ((if l = f k then (1 : ℝ) else 0)
        + cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))) ^ 2
        ≤ ∑ _l : Fin r, η ^ 2 := by
          refine Finset.sum_le_sum fun l _ => ?_
          rw [← sq_abs]
          exact pow_le_pow_left₀ (abs_nonneg _) (hsv l) 2
      _ = (r : ℝ) * η ^ 2 := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hcolsum : ∑ l : Fin r, ∑ i, Q i l ^ 2 ≤ (r : ℝ) * Cq := by
    calc ∑ l : Fin r, ∑ i, Q i l ^ 2 ≤ ∑ _l : Fin r, Cq :=
          Finset.sum_le_sum fun l _ => hcol l
      _ = (r : ℝ) * Cq := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hcnn : (0 : ℝ) ≤ (r : ℝ) * Cq := mul_nonneg (Nat.cast_nonneg r) hCq
  have hdotres : (S *ᵥ xtv W (Q.submatrix id f) ρ k
        - ρ k • xtv W (Q.submatrix id f) ρ k)
      ⬝ᵥ (S *ᵥ xtv W (Q.submatrix id f) ρ k - ρ k • xtv W (Q.submatrix id f) ρ k)
      ≤ ((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2) := by
    rw [hre]
    have h3nn : (0 : ℝ) ≤ ∑ l, ((if l = f k then (1 : ℝ) else 0)
        + cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))) ^ 2 :=
      Finset.sum_nonneg fun l _ => sq_nonneg _
    refine le_trans (dot_mulVec_le_sum_cols Q _) ?_
    exact mul_le_mul hcolsum hsvsum h3nn hcnn
  have hofLp : WithLp.ofLp (toOp S (yhatv W (Q.submatrix id f) ρ ν k)
        - ρ k • yhatv W (Q.submatrix id f) ρ ν k)
      = (Real.sqrt (ν k))⁻¹ • (S *ᵥ xtv W (Q.submatrix id f) ρ k
        - ρ k • xtv W (Q.submatrix id f) ρ k) := by
    have h1 : WithLp.ofLp (toOp S (yhatv W (Q.submatrix id f) ρ ν k)
          - ρ k • yhatv W (Q.submatrix id f) ρ ν k)
        = S *ᵥ ((Real.sqrt (ν k))⁻¹ • xtv W (Q.submatrix id f) ρ k)
          - ρ k • ((Real.sqrt (ν k))⁻¹ • xtv W (Q.submatrix id f) ρ k) := rfl
    rw [h1, Matrix.mulVec_smul, smul_comm (ρ k) ((Real.sqrt (ν k))⁻¹), ← smul_sub]
  have hns : ‖toOp S (yhatv W (Q.submatrix id f) ρ ν k)
        - ρ k • yhatv W (Q.submatrix id f) ρ ν k‖ ^ 2
      = (ν k)⁻¹ * ((S *ᵥ xtv W (Q.submatrix id f) ρ k
          - ρ k • xtv W (Q.submatrix id f) ρ k)
        ⬝ᵥ (S *ᵥ xtv W (Q.submatrix id f) ρ k
          - ρ k • xtv W (Q.submatrix id f) ρ k)) := by
    rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot, hofLp, smul_dotProduct,
      dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc, hsfmul]
  have hCqrr : (0 : ℝ) ≤ Cq * (r : ℝ) * (r : ℝ) :=
    mul_nonneg (mul_nonneg hCq (Nat.cast_nonneg r)) (Nat.cast_nonneg r)
  have hresCGsq : resCG Cq r ν ^ 2
      = Cq * (r : ℝ) * (r : ℝ) * (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 := by
    rw [resCG, mul_pow, Real.sq_sqrt hCqrr]
  have hνle : (ν k)⁻¹ ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 := by
    have h1 := inv_sqrt_le_sum ν k
    have h2 : ((Real.sqrt (ν k))⁻¹) ^ 2 ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 :=
      pow_le_pow_left₀ (by positivity) h1 2
    rwa [inv_pow, Real.sq_sqrt (hν k).le] at h2
  have hprod_nn : (0 : ℝ) ≤ ((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2) :=
    mul_nonneg hcnn (by positivity)
  have hsq : ‖toOp S (yhatv W (Q.submatrix id f) ρ ν k)
        - ρ k • yhatv W (Q.submatrix id f) ρ ν k‖ ^ 2 ≤ (resCG Cq r ν * η) ^ 2 := by
    rw [hns]
    calc (ν k)⁻¹ * ((S *ᵥ xtv W (Q.submatrix id f) ρ k
          - ρ k • xtv W (Q.submatrix id f) ρ k)
        ⬝ᵥ (S *ᵥ xtv W (Q.submatrix id f) ρ k
          - ρ k • xtv W (Q.submatrix id f) ρ k))
        ≤ (ν k)⁻¹ * (((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2)) :=
          mul_le_mul_of_nonneg_left hdotres (inv_nonneg.mpr (hν k).le)
      _ ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 * (((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2)) :=
          mul_le_mul_of_nonneg_right hνle hprod_nn
      _ = (resCG Cq r ν * η) ^ 2 := by rw [mul_pow, hresCGsq]; ring
  exact (le_abs_self _).trans
    (abs_le_of_sq_le_sq hsq (mul_nonneg (resCG_nonneg Cq r ν) hη0))

/-- The Gram bound of the sub-frame version (the mirror of step 3 of `align_det`). The
injectivity of `f` is what turns a pair `k ≠ l` of frame indices into a pair
`f k ≠ f l` of split columns, where `hE1` gives the off-diagonal bound. -/
theorem gram_bound_sub {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    (hW : W.IsHermitian) {f : Fin s → Fin r} (hf : Function.Injective f)
    {ρ ν : Fin s → ℝ} (hν : ∀ k, 0 < ν k) (hz : ∀ k, lamMax W hW < ρ k)
    {η : ℝ} (hη0 : 0 ≤ η)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η)
    (hE2 : ∀ k l : Fin s, |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))
      - (if k = l then ν k else 0)| ≤ η) (k l : Fin s) :
    |⟪yhatv W (Q.submatrix id f) ρ ν k, yhatv W (Q.submatrix id f) ρ ν l⟫_ℝ
      - if k = l then 1 else 0| ≤ gramC ρ ν * η := by
  classical
  have hsfmul : ∀ k', (Real.sqrt (ν k'))⁻¹ * (Real.sqrt (ν k'))⁻¹ = (ν k')⁻¹ :=
    fun k' => by rw [← mul_inv, Real.mul_self_sqrt (hν k').le]
  have hinner : ∀ k' l' : Fin s,
      ⟪yhatv W (Q.submatrix id f) ρ ν k', yhatv W (Q.submatrix id f) ρ ν l'⟫_ℝ
      = (Real.sqrt (ν k'))⁻¹ * ((Real.sqrt (ν l'))⁻¹
        * (xtv W (Q.submatrix id f) ρ k' ⬝ᵥ xtv W (Q.submatrix id f) ρ l')) := by
    intro k' l'
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W (Q.submatrix id f) ρ ν k')
        = (Real.sqrt (ν k'))⁻¹ • xtv W (Q.submatrix id f) ρ k' := rfl
    have h2 : WithLp.ofLp (yhatv W (Q.submatrix id f) ρ ν l')
        = (Real.sqrt (ν l'))⁻¹ • xtv W (Q.submatrix id f) ρ l' := rfl
    rw [h1, h2, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  rcases eq_or_ne k l with rfl | hkl
  · have htie : xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ k
        = cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f k)) := by
      rw [xtv_dot_xtv hW k k]
      rfl
    have hb := hE2 k k
    rw [if_pos rfl] at hb
    have hval : ⟪yhatv W (Q.submatrix id f) ρ ν k, yhatv W (Q.submatrix id f) ρ ν k⟫_ℝ
        - (if k = k then (1 : ℝ) else 0)
        = (ν k)⁻¹ * (cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f k)) - ν k) := by
      rw [if_pos rfl, hinner k k, htie, ← mul_assoc, hsfmul k, mul_sub,
        inv_mul_cancel₀ (hν k).ne']
    rw [hval, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
    have hgt : (ν k)⁻¹ ≤ gramC ρ ν := by
      have h1 := gramTerm_le_gramC ρ ν k k
      rwa [if_pos rfl, one_mul, Real.sqrt_mul_self (hν k).le] at h1
    calc (ν k)⁻¹ * |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f k)) - ν k|
        ≤ (ν k)⁻¹ * η := mul_le_mul_of_nonneg_left hb (inv_nonneg.mpr (hν k).le)
      _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0
  · have hfkl : f k ≠ f l := fun h => hkl (hf h)
    have habs : |⟪yhatv W (Q.submatrix id f) ρ ν k, yhatv W (Q.submatrix id f) ρ ν l⟫_ℝ
          - if k = l then (1 : ℝ) else 0|
        = (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹
          * |xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l|) := by
      rw [if_neg hkl, sub_zero, hinner k l, abs_mul, abs_mul,
        abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν k))),
        abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν l)))]
    rw [habs]
    have hfac : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ = (Real.sqrt (ν k * ν l))⁻¹ := by
      rw [Real.sqrt_mul (hν k).le, mul_inv]
    rcases eq_or_ne (ρ k) (ρ l) with hρkl | hρkl
    · -- tied outliers: a `G₀²` form
      have htie : xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l
          = cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l)) := by
        rw [xtv_dot_xtv hW k l, show resolv W (ρ l) = resolv W (ρ k) by rw [hρkl]]
        rfl
      have hb := hE2 k l
      rw [if_neg hkl, sub_zero] at hb
      have hgt : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ ≤ gramC ρ ν := by
        have h1 := gramTerm_le_gramC ρ ν k l
        rw [if_pos hρkl, one_mul] at h1
        rw [hfac]
        exact h1
      calc (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹
            * |xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l|)
          ≤ (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * η) := by
            rw [htie]
            exact mul_le_mul_of_nonneg_left
              (mul_le_mul_of_nonneg_left hb (by positivity)) (by positivity)
        _ = (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * η := by ring
        _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0
    · -- distinct outliers: the resolvent identity `R4.resolv_sub_resolv`
      have hne2 : ρ l - ρ k ≠ 0 := sub_ne_zero.mpr (Ne.symm hρkl)
      have hrs := resolv_sub_resolv (W₀ := W) (z₁ := ρ k) (z₂ := ρ l)
        (isUnit_det_sub hW (hz k)) (isUnit_det_sub hW (hz l))
      have hprod : resolv W (ρ k) * resolv W (ρ l)
          = (ρ l - ρ k)⁻¹ • (resolv W (ρ l) - resolv W (ρ k)) := by
        rw [hrs, smul_smul, inv_mul_cancel₀ hne2, one_smul]
      have hval : xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l
          = (ρ l - ρ k)⁻¹ * (cform W (ρ l) (fun i => Q i (f k)) (fun i => Q i (f l))
            - cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))) := by
        rw [xtv_dot_xtv hW k l, hprod, Matrix.smul_mulVec, dotProduct_smul,
          smul_eq_mul, Matrix.sub_mulVec, dotProduct_sub]
        rfl
      have hb1 := hE1 (f k) l
      rw [if_neg hfkl, sub_zero] at hb1
      have hb2 := hE1 (f l) k
      rw [if_neg (Ne.symm hfkl), sub_zero] at hb2
      have hb2' : |cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))| ≤ η := by
        rw [FormsR.cform_comm hW]
        exact hb2
      have hdiff : |cform W (ρ l) (fun i => Q i (f k)) (fun i => Q i (f l))
          - cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))| ≤ η + η := by
        calc |cform W (ρ l) (fun i => Q i (f k)) (fun i => Q i (f l))
            - cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))|
            ≤ |cform W (ρ l) (fun i => Q i (f k)) (fun i => Q i (f l))|
              + |cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))| := abs_sub _ _
          _ ≤ η + η := add_le_add hb1 hb2'
      have hdotb : |xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l|
          ≤ 2 / |ρ k - ρ l| * η := by
        rw [hval, abs_mul, abs_inv, abs_sub_comm (ρ l) (ρ k)]
        calc |ρ k - ρ l|⁻¹ * |cform W (ρ l) (fun i => Q i (f k)) (fun i => Q i (f l))
            - cform W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))|
            ≤ |ρ k - ρ l|⁻¹ * (η + η) :=
              mul_le_mul_of_nonneg_left hdiff (by positivity)
          _ = 2 / |ρ k - ρ l| * η := by rw [div_eq_mul_inv]; ring
      have hgt : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * (2 / |ρ k - ρ l|)
          ≤ gramC ρ ν := by
        have h1 := gramTerm_le_gramC ρ ν k l
        rw [if_neg hρkl] at h1
        rw [hfac]
        calc (Real.sqrt (ν k * ν l))⁻¹ * (2 / |ρ k - ρ l|)
            = 2 / |ρ k - ρ l| * (Real.sqrt (ν k * ν l))⁻¹ := by ring
          _ ≤ gramC ρ ν := h1
      calc (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹
            * |xtv W (Q.submatrix id f) ρ k ⬝ᵥ xtv W (Q.submatrix id f) ρ l|)
          ≤ (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * (2 / |ρ k - ρ l| * η)) :=
            mul_le_mul_of_nonneg_left
              (mul_le_mul_of_nonneg_left hdotb (by positivity)) (by positivity)
        _ = (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * (2 / |ρ k - ρ l|) * η := by ring
        _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0

/-! ### 3. The sub-frame core -/

/-- **The deterministic core with a sub-frame** (gap G2 of `notes/archive/rankr_plan_A.md`). The
rank-1 mirror is `align_det`. The split `S = W + Q Qᵀ` keeps its `r` columns, the frame
takes the `s` columns `f 0, ..., f (s-1)`, and the count of eigenvalues above `τ` is the
hypothesis `hI`, in the exact shape `Frame.specProj_frame_approx` takes. `Cq` bounds every
column norm of `Q`, `t k` is the target of the `k`-th overlap and `Lb` bounds it. Ties among
the `ρ_k` are allowed. -/
theorem align_detG {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    (hW : W.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {f : Fin s → Fin r} (hf : Function.Injective f)
    {ρ ν : Fin s → ℝ} {τ mg : ℝ} (hmg : 0 < mg)
    (hτ : lamMax W hW ≤ τ) (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {Cq : ℝ} (hCq : 0 ≤ Cq) (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    {v : Fin p → ℝ} (hv : v ⬝ᵥ v = 1)
    {t : Fin s → ℝ} {Lb η : ℝ} (hLb : ∀ k, |t k| ≤ Lb) (hη0 : 0 < η) (hη1 : η ≤ 1)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η)
    (hE2 : ∀ k l : Fin s, |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))
      - (if k = l then ν k else 0)| ≤ η)
    (hE3 : ∀ k : Fin s, |cform W (ρ k) v (fun i => Q i (f k)) - t k| ≤ η)
    (hsmall : (s : ℝ) * (gramC ρ ν * η + resCG Cq r ν * η / mg) ≤ 1 / 2) :
    |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2 - ∑ k, t k ^ 2 / ν k|
      ≤ 5 * s * (gramC ρ ν * η + resCG Cq r ν * η / mg)
        + (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by
  classical
  -- 0. the frame points sit above `lamMax W`
  have hτz : ∀ k, lamMax W hW < ρ k := fun k =>
    lt_of_le_of_lt hτ (lt_of_lt_of_le (lt_add_of_pos_right τ hmg) (hρ k))
  -- 1. the two frame estimates
  have hresid : ∀ k, ‖toOp S (yhatv W (Q.submatrix id f) ρ ν k)
      - ρ k • yhatv W (Q.submatrix id f) ρ ν k‖ ≤ resCG Cq r ν * η :=
    fun k => residual_bound_sub hW hSeq hν hτz hCq hη0.le hcol hE1 k
  have hgrame : ∀ k l, |⟪yhatv W (Q.submatrix id f) ρ ν k,
        yhatv W (Q.submatrix id f) ρ ν l⟫_ℝ - if k = l then 1 else 0| ≤ gramC ρ ν * η :=
    fun k l => gram_bound_sub hW hf hν hτz hη0.le hE1 hE2 k l
  -- 2. the frame approximation, with the count taken as a hypothesis
  have hδ0 : (0 : ℝ) ≤ resCG Cq r ν * η := mul_nonneg (resCG_nonneg Cq r ν) hη0.le
  have hδ'0 : (0 : ℝ) ≤ gramC ρ ν * η := mul_nonneg (gramC_nonneg ρ ν) hη0.le
  have hframe := Frame.specProj_frame_approx hS hmg hδ0 hδ'0 hI hρ hresid hgrame hsmall
    (WithLp.toLp 2 v)
  have hx1 : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot, WithLp.ofLp_toLp, hv]
  rw [hx1, mul_one] at hframe
  -- 3. the overlaps against the target `t`
  have hover : ∀ k, ⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ
      = (Real.sqrt (ν k))⁻¹ * cform W (ρ k) v (fun i => Q i (f k)) := by
    intro k
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W (Q.submatrix id f) ρ ν k)
        = (Real.sqrt (ν k))⁻¹ • xtv W (Q.submatrix id f) ρ k := rfl
    rw [h1, WithLp.ofLp_toLp, smul_dotProduct, smul_eq_mul]
    congr 1
    rw [dotProduct_comm]
    rfl
  have hterm : ∀ k, |⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
      - t k ^ 2 / ν k| ≤ (1 + 2 * Lb) * (ν k)⁻¹ * η := by
    intro k
    have hb := hE3 k
    have h1 : ⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        = (ν k)⁻¹ * cform W (ρ k) v (fun i => Q i (f k)) ^ 2 := by
      rw [hover k, mul_pow, inv_pow, Real.sq_sqrt (hν k).le]
    have h2 : |cform W (ρ k) v (fun i => Q i (f k)) ^ 2 - t k ^ 2| ≤ (1 + 2 * Lb) * η := by
      have habs2 : |cform W (ρ k) v (fun i => Q i (f k)) ^ 2 - t k ^ 2|
          = |cform W (ρ k) v (fun i => Q i (f k)) - t k|
            * |cform W (ρ k) v (fun i => Q i (f k)) + t k| := by
        rw [← abs_mul]
        congr 1
        ring
      have h3 : |cform W (ρ k) v (fun i => Q i (f k)) + t k| ≤ η + 2 * Lb := by
        have h4 : cform W (ρ k) v (fun i => Q i (f k)) + t k
            = (cform W (ρ k) v (fun i => Q i (f k)) - t k) + 2 * t k := by ring
        rw [h4]
        calc |(cform W (ρ k) v (fun i => Q i (f k)) - t k) + 2 * t k|
            ≤ |cform W (ρ k) v (fun i => Q i (f k)) - t k| + |2 * t k| := abs_add_le _ _
          _ ≤ η + 2 * Lb := by
              rw [abs_mul, abs_two]
              exact add_le_add hb (by linarith [hLb k])
      rw [habs2]
      calc |cform W (ρ k) v (fun i => Q i (f k)) - t k|
          * |cform W (ρ k) v (fun i => Q i (f k)) + t k|
          ≤ η * (η + 2 * Lb) := mul_le_mul hb h3 (abs_nonneg _) hη0.le
        _ ≤ (1 + 2 * Lb) * η := by nlinarith [hη0.le, hη1]
    have hgoal : |⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        - t k ^ 2 / ν k|
        = (ν k)⁻¹ * |cform W (ρ k) v (fun i => Q i (f k)) ^ 2 - t k ^ 2| := by
      rw [h1, show t k ^ 2 / ν k = (ν k)⁻¹ * t k ^ 2 by rw [div_eq_mul_inv, mul_comm],
        ← mul_sub, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
    rw [hgoal]
    calc (ν k)⁻¹ * |cform W (ρ k) v (fun i => Q i (f k)) ^ 2 - t k ^ 2|
        ≤ (ν k)⁻¹ * ((1 + 2 * Lb) * η) :=
          mul_le_mul_of_nonneg_left h2 (inv_nonneg.mpr (hν k).le)
      _ = (1 + 2 * Lb) * (ν k)⁻¹ * η := by ring
  have hsum : |∑ k, ⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
      - ∑ k, t k ^ 2 / ν k| ≤ (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by
    rw [← Finset.sum_sub_distrib]
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ k, |⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - t k ^ 2 / ν k|
        ≤ ∑ k, (1 + 2 * Lb) * (ν k)⁻¹ * η := Finset.sum_le_sum fun k _ => hterm k
      _ = (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by
          rw [← Finset.sum_mul, ← Finset.mul_sum]
  -- 4. assemble
  calc |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2 - ∑ k, t k ^ 2 / ν k|
      ≤ |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2
          - ∑ k, ⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2|
        + |∑ k, ⟪yhatv W (Q.submatrix id f) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
          - ∑ k, t k ^ 2 / ν k| := abs_sub_le _ _ _
    _ ≤ 5 * s * (gramC ρ ν * η + resCG Cq r ν * η / mg)
        + (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := add_le_add hframe hsum

/-- `align_detG` at a frame index: the test vector aligns with the spike `f k₀`, so the
limit of the overlap is `L² / ν k₀`. -/
theorem align_detG_mem {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    (hW : W.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {f : Fin s → Fin r} (hf : Function.Injective f)
    {ρ ν : Fin s → ℝ} {τ mg : ℝ} (hmg : 0 < mg)
    (hτ : lamMax W hW ≤ τ) (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {Cq : ℝ} (hCq : 0 ≤ Cq) (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    {v : Fin p → ℝ} (hv : v ⬝ᵥ v = 1) (k₀ : Fin s) {L η : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η)
    (hE2 : ∀ k l : Fin s, |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))
      - (if k = l then ν k else 0)| ≤ η)
    (hE3 : ∀ k : Fin s, |cform W (ρ k) v (fun i => Q i (f k))
      - (if k = k₀ then L else 0)| ≤ η)
    (hsmall : (s : ℝ) * (gramC ρ ν * η + resCG Cq r ν * η / mg) ≤ 1 / 2) :
    |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2 - L ^ 2 / ν k₀|
      ≤ 5 * s * (gramC ρ ν * η + resCG Cq r ν * η / mg)
        + (1 + 2 * |L|) * (∑ k, (ν k)⁻¹) * η := by
  classical
  have hLb : ∀ k : Fin s, |if k = k₀ then L else 0| ≤ |L| := by
    intro k
    rcases eq_or_ne k k₀ with rfl | hne
    · rw [if_pos rfl]
    · rw [if_neg hne, abs_zero]
      exact abs_nonneg L
  have hLksum : ∑ k, (if k = k₀ then L else 0) ^ 2 / ν k = L ^ 2 / ν k₀ := by
    have hcongr : ∀ k : Fin s, (if k = k₀ then L else 0) ^ 2 / ν k
        = if k = k₀ then L ^ 2 / ν k else 0 := by
      intro k
      rcases eq_or_ne k k₀ with rfl | hne
      · rw [if_pos rfl, if_pos rfl]
      · rw [if_neg hne, if_neg hne]
        simp
    rw [Finset.sum_congr rfl fun k _ => hcongr k, Finset.sum_ite_eq']
    simp
  have h := align_detG hW hS hSeq hf hmg hτ hρ hν hI hCq hcol hv hLb hη0 hη1 hE1 hE2
    hE3 hsmall
  rwa [hLksum] at h

/-- `align_detG` away from the frame: the test vector has no supercritical component, so
every overlap target is zero and the projector norm is `O(η)`. The consumer supplies `hE3`
from the vanishing cross limit at a `j` outside the range of `f`. -/
theorem align_detG_not_mem {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    (hW : W.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {f : Fin s → Fin r} (hf : Function.Injective f)
    {ρ ν : Fin s → ℝ} {τ mg : ℝ} (hmg : 0 < mg)
    (hτ : lamMax W hW ≤ τ) (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {Cq : ℝ} (hCq : 0 ≤ Cq) (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    {v : Fin p → ℝ} (hv : v ⬝ᵥ v = 1) {η : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η)
    (hE2 : ∀ k l : Fin s, |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))
      - (if k = l then ν k else 0)| ≤ η)
    (hE3 : ∀ k : Fin s, |cform W (ρ k) v (fun i => Q i (f k))| ≤ η)
    (hsmall : (s : ℝ) * (gramC ρ ν * η + resCG Cq r ν * η / mg) ≤ 1 / 2) :
    ‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2
      ≤ 5 * s * (gramC ρ ν * η + resCG Cq r ν * η / mg) + (∑ k, (ν k)⁻¹) * η := by
  classical
  have hLb : ∀ k : Fin s, |(0 : ℝ)| ≤ 0 := fun _ => by rw [abs_zero]
  have hE3' : ∀ k : Fin s, |cform W (ρ k) v (fun i => Q i (f k)) - 0| ≤ η := by
    intro k
    rw [sub_zero]
    exact hE3 k
  have h := align_detG (t := fun _ => (0 : ℝ)) (Lb := 0) hW hS hSeq hf hmg hτ hρ hν hI
    hCq hcol hv hLb hη0 hη1 hE1 hE2 hE3' hsmall
  have hzero : ∑ k : Fin s, (0 : ℝ) ^ 2 / ν k = 0 := by simp
  rw [hzero, sub_zero] at h
  have h2 : (1 + 2 * (0 : ℝ)) * (∑ k, (ν k)⁻¹) * η = (∑ k, (ν k)⁻¹) * η := by ring
  rw [h2] at h
  exact (le_abs_self _).trans h

/-! ### 4. The count from the edge and the frame -/

/-- The count hypothesis of `align_detG`, from the edge bound of task U7a and the frame.
The upper bound is the sorted-index bookkeeping of `Frame.count_of_split_of_frame`, with
`hedge` in place of `Frame.eigenvalues₀_le_of_split`; the lower bound is
`Frame.le_card_filter_of_frame`. `hδm : δ ≤ mg` turns the frame smallness of `align_detG`
into the smallness that the lower bound needs. -/
theorem count_eq_of_frame_of_edge {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    {τ mg δ δ' : ℝ} (hmg : 0 < mg) (hδ : 0 ≤ δ) (hδ' : 0 ≤ δ')
    {y : Fin s → EuclideanSpace ℝ (Fin p)} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + mg ≤ ρ k)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : (s : ℝ) * (δ' + δ / mg) ≤ 1 / 2) (hδm : δ ≤ mg)
    (hedge : ∀ k : Fin (Fintype.card (Fin p)), s ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ) :
    (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s := by
  classical
  set b : Fin s → Fin p → ℝ := fun k i => ⟪hS.eigenvectorBasis i, y k⟫_ℝ with hb
  have hres' : ∀ k, ∑ i, (hS.eigenvalues i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2 := by
    intro k
    have h := pow_le_pow_left₀ (norm_nonneg _) (hres k) 2
    rwa [Frame.norm_sq_residual hS] at h
  have hgram' : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ' := by
    intro k l
    have h : b k ⬝ᵥ b l = ⟪y k, y l⟫_ℝ := by
      rw [Frame.inner_eq_sum_coord hS]
      rfl
    rw [h]
    exact hgram k l
  have hdm : δ / mg ≤ 1 := (div_le_one hmg).mpr hδm
  have hdm0 : 0 ≤ δ / mg := div_nonneg hδ hmg.le
  have hs0 : (0 : ℝ) ≤ s := Nat.cast_nonneg s
  have hsmall2 : (s : ℝ) * (δ' + (δ / mg) ^ 2) < 1 := by
    nlinarith [mul_nonneg hs0 (mul_nonneg hdm0 (sub_nonneg.mpr hdm))]
  have hlow : s ≤ (Finset.univ.filter fun i => τ < hS.eigenvalues i).card :=
    Frame.le_card_filter_of_frame hmg hδ' hρ hres' hgram' hsmall2
  have hsp : s ≤ p := by
    refine le_trans hlow ?_
    calc (Finset.univ.filter fun i => τ < hS.eigenvalues i).card
        ≤ (Finset.univ : Finset (Fin p)).card := Finset.card_filter_le _ _
      _ = p := by simp
  have hup : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card ≤ s := by
    rw [Frame.card_filter_eigenvalues hS (fun w => τ < w)]
    calc (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k).card
        ≤ (Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < s).card := by
          refine Finset.card_le_card fun k hk => ?_
          rw [Finset.mem_filter] at hk ⊢
          refine ⟨hk.1, ?_⟩
          by_contra hge
          exact absurd hk.2 (not_lt.mpr (hedge k (not_lt.mp hge)))
      _ = s := card_filter_lt (by simpa using hsp)
  exact le_antisymm hup hlow

/-- `count_eq_of_frame_of_edge` in the form the consumer already has: the two frame
estimates come from `residual_bound_sub` and `gram_bound_sub`, so the input is the pair of
entrywise accuracies `hE1` and `hE2` that `align_detG` also takes. -/
theorem count_eq_of_forms_of_edge {W S : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) {f : Fin s → Fin r} (hf : Function.Injective f)
    {ρ ν : Fin s → ℝ} {τ mg : ℝ} (hmg : 0 < mg) (hτ : lamMax W hW ≤ τ)
    (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    {Cq η : ℝ} (hCq : 0 ≤ Cq) (hη0 : 0 < η)
    (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => Q i (f k))
        - (if l = f k then -1 else 0)| ≤ η)
    (hE2 : ∀ k l : Fin s, |cform2 W (ρ k) (fun i => Q i (f k)) (fun i => Q i (f l))
      - (if k = l then ν k else 0)| ≤ η)
    (hsmall : (s : ℝ) * (gramC ρ ν * η + resCG Cq r ν * η / mg) ≤ 1 / 2)
    (hδm : resCG Cq r ν * η ≤ mg)
    (hedge : ∀ k : Fin (Fintype.card (Fin p)), s ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ) :
    (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s := by
  have hτz : ∀ k, lamMax W hW < ρ k := fun k =>
    lt_of_le_of_lt hτ (lt_of_lt_of_le (lt_add_of_pos_right τ hmg) (hρ k))
  exact count_eq_of_frame_of_edge hS hmg (mul_nonneg (resCG_nonneg Cq r ν) hη0.le)
    (mul_nonneg (gramC_nonneg ρ ν) hη0.le) hρ
    (fun k => residual_bound_sub hW hSeq hν hτz hCq hη0.le hcol hE1 k)
    (fun k l => gram_bound_sub hW hf hν hτz hη0.le hE1 hE2 k l) hsmall hδm hedge

/-- The edge bound at `s` weakens to the edge bound at `r` when `s ≤ r`. This is the shape
`Frame.norm_sq_specProjTop_split` and `Frame.specProjTop_eq_specProj_add_edge` need. -/
theorem eigenvalues₀_le_of_edge_le {S : Matrix (Fin p) (Fin p) ℝ} (hS : S.IsHermitian)
    {τ : ℝ} (hsr : s ≤ r)
    (hedge : ∀ k : Fin (Fintype.card (Fin p)), s ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ) :
    ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ :=
  fun k hk => hedge k (le_trans hsr hk)

end Core

end OutliersR

end StackedSVD
