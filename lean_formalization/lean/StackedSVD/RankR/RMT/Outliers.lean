/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Forms
import StackedSVD.RankR.RMT.Stack
import StackedSVD.LinAlg.Frame

/-!
# Task U4: outliers, the projector, and the supercritical discharge

Task U4 of `notes/archive/plan_subspacelaw.md` (Route P, section 1.3, and the U4 row of
section 3). The approximate eigenvectors are `x̃_k = G₀(ρ_k) Q e_k` at the deterministic
outlier `ρ_k = rhoSq (√λ_k) c`, normalized by the deterministic `√ν_k`,
`ν_k = (λ_k + 1) m'(ρ_k)`. Every spike is supercritical (`c < λ_k²`); ties among the
`ρ_k` are allowed.

1. **Scalars** (section 1). `b < ρ_k` (`MP.bulkEdge_lt_rhoSq`), `(λ_k + 1) m(ρ_k) = -1`
   (`MP.m_rhoSq`), `0 < ν_k` (`MP.mDeriv_pos`), and `(√λ_j m(ρ_j))² / ν_j = betaSq`
   (`MP.overlap_identity`).
2. **Deterministic core** (sections 2 and 3, model-free). On the event
   `lamMax W₀ ≤ τ < ρ_k` with entrywise form accuracies `η`:
   the residual `(S - ρ_k) x̃_k = Q S(ρ_k) e_k` vanishes at rate `η` (plan 1.3 item 1,
   `residual_eq`); the Gram of the normalized frame is within `O(η)` of `1`
   (item 2, ties by `cform2`, distinct `ρ` by `R4.resolv_sub_resolv`); the columns of `Q`
   have `‖q_l‖² ≤ 2 ρ_l` (`dot_self_le_mul_neg_qform`, `W₀` psd); `align_det` then feeds
   `Frame.specProjTop_frame_approx_of_split` (items 3 to 5) and compares the overlaps
   `⟪ŷ_k, v_j⟫` with `√λ_k m(ρ_k) δ_jk / √ν_k`.
3. **The discharge** (section 4). `subspaceLaw_of_gaussian_supercritical` consumes
   `UnalignedModel.resolventLimitsR_of_gaussian` (task U2) and concludes
   `m.SubspaceLaw (∑ i, cc i)`, the exact hypothesis of `prop_stacksvd_subspace`.

What remains for U7: the mixed case (a subcritical spike contributes `0` through the
edge block, plan 1.4), the D11 tail shift that removes `hn`/`hp`, and the removal of the
per-`N` side condition `r < ∑ i, n i N`. No `sorry` here: this file states and proves
only the all-supercritical case.

Numeric check before the proofs (a session script, `check_outliers.py`, not kept): seed
`2026083003`, `d = 400`, `M = 2`, `r = 2`, distinct spikes `θ = (2.0, 1.6)` and tied
spikes `θ = (1.8, 1.8)` at `c = 1.5`, all supercritical: split residual `2e-15`, count
`= r` in both, frame against projector within `0.05`, `perfStackR` against `limitStackR`
within `0.13` (distinct, weakest gap `ρ - b = 0.02`) and `0.04` (tied); the scalar
identities `L²/ν = betaSq` and `(λ+1)m(ρ) = -1` to `1e-6`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace OutliersR

open R4

/-! ### 1. Scalars of the supercritical outlier -/

section Scalars

variable {c lam : ℝ}

theorem lam_pos (hc : 0 < c) (hlam0 : 0 ≤ lam) (hsup : c < lam ^ 2) : 0 < lam := by
  rcases hlam0.lt_or_eq with h | h
  · exact h
  · exfalso
    rw [← h] at hsup
    norm_num at hsup
    linarith

theorem sqrt_pow_four (hlam0 : 0 ≤ lam) : Real.sqrt lam ^ 4 = lam ^ 2 := by
  have h : Real.sqrt lam ^ 4 = (Real.sqrt lam ^ 2) ^ 2 := by ring
  rw [h, Real.sq_sqrt hlam0]

theorem bulkEdge_lt_rho (hc : 0 < c) (hlam0 : 0 ≤ lam) (hsup : c < lam ^ 2) :
    bulkEdge c < rhoSq (Real.sqrt lam) c :=
  MP.bulkEdge_lt_rhoSq hc (Real.sqrt_pos.mpr (lam_pos hc hlam0 hsup))
    (by rw [sqrt_pow_four hlam0]; exact hsup)

/-- `(λ + 1) m(ρ) = -1` at the supercritical outlier: the secular limit vanishes. -/
theorem mul_m_rho_eq_neg_one (hc : 0 < c) (hlam0 : 0 ≤ lam) (hsup : c < lam ^ 2) :
    (lam + 1) * MP.m c (rhoSq (Real.sqrt lam) c) = -1 := by
  have hpos := lam_pos hc hlam0 hsup
  have hθ : 0 < Real.sqrt lam := Real.sqrt_pos.mpr hpos
  have h4 : c < Real.sqrt lam ^ 4 := by rw [sqrt_pow_four hlam0]; exact hsup
  rw [MP.m_rhoSq hc hθ h4, Real.sq_sqrt hlam0]
  have h1 : lam + 1 ≠ 0 := by linarith
  field_simp

/-- `ν = (λ + 1) m'(ρ) > 0`, the normalizer of the frame. -/
theorem nu_pos (hc : 0 < c) (hlam0 : 0 ≤ lam) (hsup : c < lam ^ 2) :
    0 < (lam + 1) * MP.mDeriv c (rhoSq (Real.sqrt lam) c) := by
  have hpos := lam_pos hc hlam0 hsup
  exact mul_pos (by linarith) (MP.mDeriv_pos hc (bulkEdge_lt_rho hc hlam0 hsup))

/-- `(√λ m(ρ))² / ν = betaSq (√λ) c`: `MP.overlap_identity` in the shape the frame sum
produces (plan 1.3 item 5). -/
theorem sq_div_nu_eq_betaSq (hc : 0 < c) (hlam0 : 0 ≤ lam) (hsup : c < lam ^ 2) :
    (Real.sqrt lam * MP.m c (rhoSq (Real.sqrt lam) c)) ^ 2
      / ((lam + 1) * MP.mDeriv c (rhoSq (Real.sqrt lam) c))
      = betaSq (Real.sqrt lam) c := by
  have hpos := lam_pos hc hlam0 hsup
  have hθ : 0 < Real.sqrt lam := Real.sqrt_pos.mpr hpos
  have h4 : c < Real.sqrt lam ^ 4 := by rw [sqrt_pow_four hlam0]; exact hsup
  calc (Real.sqrt lam * MP.m c (rhoSq (Real.sqrt lam) c)) ^ 2
      / ((lam + 1) * MP.mDeriv c (rhoSq (Real.sqrt lam) c))
      = Real.sqrt lam ^ 2 * MP.m c (rhoSq (Real.sqrt lam) c) ^ 2
        / ((Real.sqrt lam ^ 2 + 1) * MP.mDeriv c (rhoSq (Real.sqrt lam) c)) := by
        rw [mul_pow, Real.sq_sqrt hlam0]
    _ = betaSq (Real.sqrt lam) c := MP.overlap_identity hc hθ h4

end Scalars

/-! ### 2. The deterministic constants -/

/-- The residual constant: `‖(S - ρ_k) ŷ_k‖ ≤ resC η` on the good event. -/
noncomputable def resC {r : ℕ} (ρ ν : Fin r → ℝ) : ℝ :=
  Real.sqrt (2 * (∑ l, ρ l) * r) * ∑ k, (Real.sqrt (ν k))⁻¹

/-- The Gram constant: `|⟪ŷ_k, ŷ_l⟫ - δ_kl| ≤ gramC η` on the good event. The sum
dominates every per-pair constant, so no minimum over pairs is needed. -/
noncomputable def gramC {r : ℕ} (ρ ν : Fin r → ℝ) : ℝ :=
  ∑ k, ∑ l, (if ρ k = ρ l then 1 else 2 / |ρ k - ρ l|) * (Real.sqrt (ν k * ν l))⁻¹

theorem resC_nonneg {r : ℕ} (ρ ν : Fin r → ℝ) : 0 ≤ resC ρ ν := by
  refine mul_nonneg (Real.sqrt_nonneg _) (Finset.sum_nonneg fun k _ => ?_)
  positivity

private theorem gramTerm_nonneg {r : ℕ} (ρ ν : Fin r → ℝ) (k l : Fin r) :
    0 ≤ (if ρ k = ρ l then 1 else 2 / |ρ k - ρ l|) * (Real.sqrt (ν k * ν l))⁻¹ := by
  refine mul_nonneg ?_ (by positivity)
  split
  · norm_num
  · positivity

theorem gramC_nonneg {r : ℕ} (ρ ν : Fin r → ℝ) : 0 ≤ gramC ρ ν :=
  Finset.sum_nonneg fun k _ => Finset.sum_nonneg fun l _ => gramTerm_nonneg ρ ν k l

theorem gramTerm_le_gramC {r : ℕ} (ρ ν : Fin r → ℝ) (k l : Fin r) :
    (if ρ k = ρ l then 1 else 2 / |ρ k - ρ l|) * (Real.sqrt (ν k * ν l))⁻¹
      ≤ gramC ρ ν := by
  calc (if ρ k = ρ l then 1 else 2 / |ρ k - ρ l|) * (Real.sqrt (ν k * ν l))⁻¹
      ≤ ∑ l', (if ρ k = ρ l' then 1 else 2 / |ρ k - ρ l'|) * (Real.sqrt (ν k * ν l'))⁻¹ :=
        Finset.single_le_sum (fun l' _ => gramTerm_nonneg ρ ν k l') (Finset.mem_univ l)
    _ ≤ gramC ρ ν :=
        Finset.single_le_sum
          (fun k' _ => Finset.sum_nonneg fun l' _ => gramTerm_nonneg ρ ν k' l')
          (Finset.mem_univ k)

theorem inv_sqrt_le_sum {r : ℕ} (ν : Fin r → ℝ) (k : Fin r) :
    (Real.sqrt (ν k))⁻¹ ≤ ∑ k', (Real.sqrt (ν k'))⁻¹ :=
  Finset.single_le_sum (f := fun k' => (Real.sqrt (ν k'))⁻¹)
    (fun k' _ => by positivity) (Finset.mem_univ k)

/-! ### 3. Deterministic linear algebra -/

section Deterministic

variable {p r : ℕ} {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}

/-- For a psd `W` with `lamMax W < z`: `y ⬝ᵥ y ≤ z (-(Φ_y z))`. This replaces the raw
`‖g‖` bound of the plan: the column norms of `Q` are read off the diagonal forms. -/
theorem dot_self_le_mul_neg_qform (hW : W.IsHermitian)
    (hpsd : ∀ a, 0 ≤ hW.eigenvalues a) {z : ℝ} (hz : lamMax W hW < z) (y : Fin p → ℝ) :
    y ⬝ᵥ y ≤ z * (-qform W z y) := by
  rw [qform_eq_sum hW hz, ← dotProduct_transpose_eigU hW y]
  have hyy : ((eigU hW)ᵀ *ᵥ y) ⬝ᵥ ((eigU hW)ᵀ *ᵥ y) = ∑ a, ((eigU hW)ᵀ *ᵥ y) a ^ 2 := by
    simp [dotProduct, sq]
  have hrw : z * (-∑ a, (hW.eigenvalues a - z)⁻¹ * ((eigU hW)ᵀ *ᵥ y) a ^ 2)
      = ∑ a, z * ((z - hW.eigenvalues a)⁻¹ * ((eigU hW)ᵀ *ᵥ y) a ^ 2) := by
    rw [← Finset.sum_neg_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [show (hW.eigenvalues a - z)⁻¹ = -(z - hW.eigenvalues a)⁻¹ by
      rw [← neg_sub z (hW.eigenvalues a), inv_neg]]
    ring
  rw [hyy, hrw]
  refine Finset.sum_le_sum fun a _ => ?_
  have hla : hW.eigenvalues a < z := lt_of_le_of_lt (eigenvalues_le_lamMax hW a) hz
  have hzl : 0 < z - hW.eigenvalues a := by linarith
  have hpos : 0 < (z - hW.eigenvalues a)⁻¹ := inv_pos.mpr hzl
  have hone : (z - hW.eigenvalues a) * (z - hW.eigenvalues a)⁻¹ = 1 :=
    mul_inv_cancel₀ hzl.ne'
  nlinarith [sq_nonneg (((eigU hW)ᵀ *ᵥ y) a), hpsd a,
    mul_nonneg (mul_nonneg (hpsd a) hpos.le) (sq_nonneg (((eigU hW)ᵀ *ᵥ y) a))]

/-- Cauchy-Schwarz per row: `‖Q w‖² ≤ (∑ column norms²) ‖w‖²`. -/
theorem dot_mulVec_le_sum_cols (Q : Matrix (Fin p) (Fin r) ℝ) (w : Fin r → ℝ) :
    (Q *ᵥ w) ⬝ᵥ (Q *ᵥ w) ≤ (∑ l, ∑ i, Q i l ^ 2) * ∑ l, w l ^ 2 := by
  have hrow : ∀ i, (Q *ᵥ w) i ^ 2 ≤ (∑ l, Q i l ^ 2) * ∑ l, w l ^ 2 := by
    intro i
    have h : (∑ l, Q i l * w l) ^ 2 ≤ (∑ l, Q i l ^ 2) * ∑ l, w l ^ 2 :=
      Finset.sum_mul_sq_le_sq_mul_sq Finset.univ _ _
    simpa [Matrix.mulVec, dotProduct] using h
  have hdot : (Q *ᵥ w) ⬝ᵥ (Q *ᵥ w) = ∑ i, (Q *ᵥ w) i ^ 2 := by simp [dotProduct, sq]
  calc (Q *ᵥ w) ⬝ᵥ (Q *ᵥ w) = ∑ i, (Q *ᵥ w) i ^ 2 := hdot
    _ ≤ ∑ i, (∑ l, Q i l ^ 2) * ∑ l, w l ^ 2 := Finset.sum_le_sum fun i _ => hrow i
    _ = (∑ i, ∑ l, Q i l ^ 2) * ∑ l, w l ^ 2 := by rw [← Finset.sum_mul]
    _ = (∑ l, ∑ i, Q i l ^ 2) * ∑ l, w l ^ 2 := by rw [Finset.sum_comm]

/-- Plan 1.3 item 1: `(S - ρ) x̃_k = Q (e_k + Qᵀ G₀(ρ) Q e_k)`, the residual identity. -/
theorem residual_eq (hW : W.IsHermitian) (hSeq : S = W + Q * Qᵀ) {z : ℝ}
    (hz : lamMax W hW < z) (k : Fin r) :
    S *ᵥ (resolv W z *ᵥ fun i => Q i k) - z • (resolv W z *ᵥ fun i => Q i k)
      = Q *ᵥ fun l =>
          (if l = k then (1 : ℝ) else 0)
            + cform W z (fun i => Q i l) (fun i => Q i k) := by
  have h1 : (W - z • (1 : Matrix (Fin p) (Fin p) ℝ)) *ᵥ (resolv W z *ᵥ fun i => Q i k)
      = fun i => Q i k := by
    rw [Matrix.mulVec_mulVec, mul_resolv hW hz, Matrix.one_mulVec]
  have h2 : Qᵀ *ᵥ (resolv W z *ᵥ fun i => Q i k)
      = fun l => cform W z (fun i => Q i l) (fun i => Q i k) := by
    funext l
    simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, cform]
  have hqk : (Q *ᵥ fun l => if l = k then (1 : ℝ) else 0) = fun i => Q i k := by
    funext i
    simp [Matrix.mulVec, dotProduct, mul_ite]
  have hsplit : S *ᵥ (resolv W z *ᵥ fun i => Q i k)
        - z • (resolv W z *ᵥ fun i => Q i k)
      = (W - z • (1 : Matrix (Fin p) (Fin p) ℝ)) *ᵥ (resolv W z *ᵥ fun i => Q i k)
        + (Q * Qᵀ) *ᵥ (resolv W z *ᵥ fun i => Q i k) := by
    rw [hSeq, Matrix.add_mulVec, Matrix.sub_mulVec, Matrix.smul_mulVec,
      Matrix.one_mulVec]
    abel
  rw [hsplit, h1, ← Matrix.mulVec_mulVec, h2, ← hqk, ← Matrix.mulVec_add]
  congr 1

end Deterministic

end OutliersR

end StackedSVD

namespace StackedSVD

namespace OutliersR

open R4

/-! ### 3b. Scalar helpers for the choice of `η` -/

theorem mul_le_of_le_div_add_one {C mgv η : ℝ} (hC : 0 ≤ C) (hη : η ≤ mgv / (C + 1))
    (hmg : 0 ≤ mgv) : C * η ≤ mgv := by
  calc C * η ≤ C * (mgv / (C + 1)) := mul_le_mul_of_nonneg_left hη hC
    _ = mgv * (C / (C + 1)) := by ring
    _ ≤ mgv * 1 :=
        mul_le_mul_of_nonneg_left (by rw [div_le_one (by linarith)]; linarith) hmg
    _ = mgv := mul_one _

theorem mul_inv_le_half {X : ℝ} (hX : 0 ≤ X) : X * (1 / (2 * (X + 1))) ≤ 1 / 2 := by
  rw [mul_one_div, div_le_div_iff₀ (by linarith) (by norm_num)]
  linarith

theorem mul_div_le_half {B ε : ℝ} (hB : 0 ≤ B) (hε : 0 ≤ ε) :
    B * (ε / (2 * (B + 1))) ≤ ε / 2 := by
  have h1 : B * (ε / (2 * (B + 1))) = ε * B / (2 * (B + 1)) := by ring
  rw [h1, div_le_div_iff₀ (by linarith) (by norm_num)]
  nlinarith [mul_nonneg hε hB]

/-! ### 3c. The frame vectors and the deterministic core -/

section Core

variable {p r : ℕ}

/-- The un-normalized approximate eigenvector `x̃_k = G₀(ρ_k) Q e_k` (plan 1.3 item 1). -/
noncomputable def xtv (W : Matrix (Fin p) (Fin p) ℝ) (Q : Matrix (Fin p) (Fin r) ℝ)
    (ρ : Fin r → ℝ) (k : Fin r) : Fin p → ℝ :=
  resolv W (ρ k) *ᵥ fun i => Q i k

/-- The normalized frame vector `ŷ_k = x̃_k / √ν_k`. The normalizer is deterministic. -/
noncomputable def yhatv (W : Matrix (Fin p) (Fin p) ℝ) (Q : Matrix (Fin p) (Fin r) ℝ)
    (ρ ν : Fin r → ℝ) (k : Fin r) : EuclideanSpace ℝ (Fin p) :=
  WithLp.toLp 2 ((Real.sqrt (ν k))⁻¹ • xtv W Q ρ k)

/-- `⟨x̃_k, x̃_l⟩ = e_kᵀ Qᵀ G₀(ρ_k) G₀(ρ_l) Q e_l` (plan 1.3 item 2). -/
theorem xtv_dot_xtv {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    {ρ : Fin r → ℝ} (hW : W.IsHermitian) (k l : Fin r) :
    xtv W Q ρ k ⬝ᵥ xtv W Q ρ l
      = (fun i => Q i k) ⬝ᵥ ((resolv W (ρ k) * resolv W (ρ l)) *ᵥ fun i => Q i l) := by
  symm
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    transpose_resolv hW]
  rfl

/-- **The deterministic core of Route P** (plan 1.3 items 1 to 5). On one sample:
the split `S = W + Q Qᵀ`, a psd `W` with `lamMax W ≤ τ`, deterministic points
`ρ_k ≥ τ + mg`, and entrywise form accuracies `η` for `Qᵀ G₀ Q + 1` (`hE1`),
`Qᵀ G₀² Q - ν` (`hE2`) and `vᵀ G₀ Q - L e_jᵀ` (`hE3`). Then the top-`r` spectral
projector applied to the unit vector `v` reproduces `L² / ν_j` up to an explicit linear
function of `η`. Ties among the `ρ_k` are allowed. -/
theorem align_det (hrp : r ≤ p) (hr : 0 < r)
    {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    (hW : W.IsHermitian) (hS : S.IsHermitian) (hpsd : ∀ a, 0 ≤ hW.eigenvalues a)
    (hSeq : S = W + Q * Qᵀ) {ρ ν : Fin r → ℝ} {τ mg : ℝ} (hmg : 0 < mg)
    (hτ : lamMax W hW ≤ τ) (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    {v : Fin p → ℝ} (hv : v ⬝ᵥ v = 1) (j : Fin r) {L η : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1)
    (hE1 : ∀ jj kk : Fin r, |cform W (ρ kk) (fun i => Q i jj) (fun i => Q i kk)
      - (if jj = kk then -1 else 0)| ≤ η)
    (hE2 : ∀ kk l : Fin r, |cform2 W (ρ kk) (fun i => Q i kk) (fun i => Q i l)
      - (if kk = l then ν kk else 0)| ≤ η)
    (hE3 : ∀ k : Fin r, |cform W (ρ k) v (fun i => Q i k)
      - (if j = k then L else 0)| ≤ η)
    (hsmall : (r : ℝ) * (gramC ρ ν * η + resC ρ ν * η / mg) ≤ 1 / 2)
    (hδm : resC ρ ν * η ≤ mg) :
    |‖specProjTop S hS r (WithLp.toLp 2 v)‖ ^ 2 - L ^ 2 / ν j|
      ≤ 5 * r * (gramC ρ ν * η + resC ρ ν * η / mg)
        + (1 + 2 * |L|) * (∑ k, (ν k)⁻¹) * η := by
  classical
  -- 0. positivity bookkeeping
  have hτz : ∀ k, lamMax W hW < ρ k := fun k =>
    lt_of_le_of_lt hτ (lt_of_lt_of_le (lt_add_of_pos_right τ hmg) (hρ k))
  have hp0 : 0 < p := lt_of_lt_of_le hr hrp
  have hlm0 : 0 ≤ lamMax W hW := by
    obtain ⟨i, hi⟩ := exists_eigenvalues_eq_lamMax hW hp0
    rw [← hi]
    exact hpsd i
  have hρpos : ∀ k, 0 < ρ k := fun k =>
    lt_of_lt_of_le (by linarith : (0 : ℝ) < τ + mg) (hρ k)
  have hsρ : (0 : ℝ) ≤ ∑ l, ρ l := Finset.sum_nonneg fun l _ => (hρpos l).le
  have hsfmul : ∀ k, (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν k))⁻¹ = (ν k)⁻¹ := fun k => by
    rw [← mul_inv, Real.mul_self_sqrt (hν k).le]
  -- 1. column norms: `‖q_l‖² ≤ 2 ρ_l` from the diagonal form and `W` psd
  have hcolsum : ∀ l, ∑ i, Q i l ^ 2 ≤ 2 * ρ l := by
    intro l
    have hb := hE1 l l
    rw [if_pos rfl, sub_neg_eq_add] at hb
    have hqeq : cform W (ρ l) (fun i => Q i l) (fun i => Q i l)
        = qform W (ρ l) (fun i => Q i l) := rfl
    rw [hqeq] at hb
    have hqle : -qform W (ρ l) (fun i => Q i l) ≤ 2 := by
      have h1 := (abs_le.mp hb).1
      linarith
    have hdot := dot_self_le_mul_neg_qform hW hpsd (hτz l) (fun i => Q i l)
    have hqq : (fun i => Q i l) ⬝ᵥ (fun i => Q i l) = ∑ i, Q i l ^ 2 := by
      simp [dotProduct, sq]
    have h2 : ρ l * (-qform W (ρ l) (fun i => Q i l)) ≤ ρ l * 2 :=
      mul_le_mul_of_nonneg_left hqle (hρpos l).le
    rw [hqq] at hdot
    linarith
  -- 2. the residual bound (plan 1.3 item 1)
  have hresid : ∀ k, ‖toOp S (yhatv W Q ρ ν k) - ρ k • yhatv W Q ρ ν k‖
      ≤ resC ρ ν * η := by
    intro k
    have hre : S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k
        = Q *ᵥ fun l => (if l = k then (1 : ℝ) else 0)
            + cform W (ρ k) (fun i => Q i l) (fun i => Q i k) :=
      residual_eq hW hSeq (hτz k) k
    have hsv : ∀ l, |(if l = k then (1 : ℝ) else 0)
        + cform W (ρ k) (fun i => Q i l) (fun i => Q i k)| ≤ η := by
      intro l
      have hb := hE1 l k
      rcases eq_or_ne l k with rfl | hlk
      · rw [if_pos rfl] at hb ⊢
        calc |1 + cform W (ρ l) (fun i => Q i l) (fun i => Q i l)|
            = |cform W (ρ l) (fun i => Q i l) (fun i => Q i l) - -1| := by
              rw [sub_neg_eq_add, add_comm]
          _ ≤ η := hb
      · rw [if_neg hlk, sub_zero] at hb
        rw [if_neg hlk, zero_add]
        exact hb
    have hsvsum : ∑ l, ((if l = k then (1 : ℝ) else 0)
        + cform W (ρ k) (fun i => Q i l) (fun i => Q i k)) ^ 2 ≤ (r : ℝ) * η ^ 2 := by
      calc ∑ l, ((if l = k then (1 : ℝ) else 0)
          + cform W (ρ k) (fun i => Q i l) (fun i => Q i k)) ^ 2
          ≤ ∑ _l : Fin r, η ^ 2 := by
            refine Finset.sum_le_sum fun l _ => ?_
            rw [← sq_abs]
            exact pow_le_pow_left₀ (abs_nonneg _) (hsv l) 2
        _ = (r : ℝ) * η ^ 2 := by
            rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    have h2sum : ∑ l, ∑ i, Q i l ^ 2 ≤ 2 * ∑ l, ρ l := by
      rw [Finset.mul_sum]
      exact Finset.sum_le_sum fun l _ => hcolsum l
    have h2nn : (0 : ℝ) ≤ 2 * ∑ l, ρ l := by linarith
    have hdotres : (S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k)
        ⬝ᵥ (S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k)
        ≤ (2 * ∑ l, ρ l) * ((r : ℝ) * η ^ 2) := by
      rw [hre]
      have h3nn : (0 : ℝ) ≤ ∑ l, ((if l = k then (1 : ℝ) else 0)
          + cform W (ρ k) (fun i => Q i l) (fun i => Q i k)) ^ 2 :=
        Finset.sum_nonneg fun l _ => sq_nonneg _
      calc (Q *ᵥ fun l => (if l = k then (1 : ℝ) else 0)
            + cform W (ρ k) (fun i => Q i l) (fun i => Q i k))
          ⬝ᵥ (Q *ᵥ fun l => (if l = k then (1 : ℝ) else 0)
            + cform W (ρ k) (fun i => Q i l) (fun i => Q i k))
          ≤ (∑ l, ∑ i, Q i l ^ 2) * ∑ l, ((if l = k then (1 : ℝ) else 0)
            + cform W (ρ k) (fun i => Q i l) (fun i => Q i k)) ^ 2 :=
            dot_mulVec_le_sum_cols Q _
        _ ≤ (2 * ∑ l, ρ l) * ((r : ℝ) * η ^ 2) := mul_le_mul h2sum hsvsum h3nn h2nn
    have hofLp : WithLp.ofLp (toOp S (yhatv W Q ρ ν k) - ρ k • yhatv W Q ρ ν k)
        = (Real.sqrt (ν k))⁻¹ • (S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k) := by
      have h1 : WithLp.ofLp (toOp S (yhatv W Q ρ ν k) - ρ k • yhatv W Q ρ ν k)
          = S *ᵥ ((Real.sqrt (ν k))⁻¹ • xtv W Q ρ k)
            - ρ k • ((Real.sqrt (ν k))⁻¹ • xtv W Q ρ k) := rfl
      rw [h1, Matrix.mulVec_smul, smul_comm (ρ k) ((Real.sqrt (ν k))⁻¹), ← smul_sub]
    have hns : ‖toOp S (yhatv W Q ρ ν k) - ρ k • yhatv W Q ρ ν k‖ ^ 2
        = (ν k)⁻¹ * ((S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k)
            ⬝ᵥ (S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k)) := by
      rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot, hofLp, smul_dotProduct,
        dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc, hsfmul k]
    have hresCsq : resC ρ ν ^ 2
        = 2 * (∑ l, ρ l) * (r : ℝ) * (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 := by
      rw [resC, mul_pow, Real.sq_sqrt (by positivity)]
    have hνle : (ν k)⁻¹ ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 := by
      have h1 := inv_sqrt_le_sum ν k
      have h2 : ((Real.sqrt (ν k))⁻¹) ^ 2 ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 :=
        pow_le_pow_left₀ (by positivity) h1 2
      rwa [inv_pow, Real.sq_sqrt (hν k).le] at h2
    have hprod_nn : (0 : ℝ) ≤ (2 * ∑ l, ρ l) * ((r : ℝ) * η ^ 2) := by positivity
    have hsq : ‖toOp S (yhatv W Q ρ ν k) - ρ k • yhatv W Q ρ ν k‖ ^ 2
        ≤ (resC ρ ν * η) ^ 2 := by
      rw [hns]
      calc (ν k)⁻¹ * ((S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k)
            ⬝ᵥ (S *ᵥ xtv W Q ρ k - ρ k • xtv W Q ρ k))
          ≤ (ν k)⁻¹ * ((2 * ∑ l, ρ l) * ((r : ℝ) * η ^ 2)) :=
            mul_le_mul_of_nonneg_left hdotres (inv_nonneg.mpr (hν k).le)
        _ ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 * ((2 * ∑ l, ρ l) * ((r : ℝ) * η ^ 2)) :=
            mul_le_mul_of_nonneg_right hνle hprod_nn
        _ = (resC ρ ν * η) ^ 2 := by rw [mul_pow, hresCsq]; ring
    exact (le_abs_self _).trans
      (abs_le_of_sq_le_sq hsq (mul_nonneg (resC_nonneg ρ ν) hη0.le))
  -- 3. the Gram bound (plan 1.3 item 2)
  have hinner : ∀ k l, ⟪yhatv W Q ρ ν k, yhatv W Q ρ ν l⟫_ℝ
      = (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * (xtv W Q ρ k ⬝ᵥ xtv W Q ρ l)) := by
    intro k l
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W Q ρ ν k) = (Real.sqrt (ν k))⁻¹ • xtv W Q ρ k := rfl
    have h2 : WithLp.ofLp (yhatv W Q ρ ν l) = (Real.sqrt (ν l))⁻¹ • xtv W Q ρ l := rfl
    rw [h1, h2, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  have hgrame : ∀ k l, |⟪yhatv W Q ρ ν k, yhatv W Q ρ ν l⟫_ℝ - if k = l then 1 else 0|
      ≤ gramC ρ ν * η := by
    intro k l
    rcases eq_or_ne k l with rfl | hkl
    · have htie : xtv W Q ρ k ⬝ᵥ xtv W Q ρ k
          = cform2 W (ρ k) (fun i => Q i k) (fun i => Q i k) := by
        rw [xtv_dot_xtv hW k k]
        rfl
      have hb := hE2 k k
      rw [if_pos rfl] at hb
      have hval : ⟪yhatv W Q ρ ν k, yhatv W Q ρ ν k⟫_ℝ - (if k = k then (1 : ℝ) else 0)
          = (ν k)⁻¹ * (cform2 W (ρ k) (fun i => Q i k) (fun i => Q i k) - ν k) := by
        rw [if_pos rfl, hinner k k, htie, ← mul_assoc, hsfmul k, mul_sub,
          inv_mul_cancel₀ (hν k).ne']
      rw [hval, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
      have hgt : (ν k)⁻¹ ≤ gramC ρ ν := by
        have h1 := gramTerm_le_gramC ρ ν k k
        rwa [if_pos rfl, one_mul, Real.sqrt_mul_self (hν k).le] at h1
      calc (ν k)⁻¹ * |cform2 W (ρ k) (fun i => Q i k) (fun i => Q i k) - ν k|
          ≤ (ν k)⁻¹ * η := mul_le_mul_of_nonneg_left hb (inv_nonneg.mpr (hν k).le)
        _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0.le
    · have habs : |⟪yhatv W Q ρ ν k, yhatv W Q ρ ν l⟫_ℝ - if k = l then (1 : ℝ) else 0|
          = (Real.sqrt (ν k))⁻¹
            * ((Real.sqrt (ν l))⁻¹ * |xtv W Q ρ k ⬝ᵥ xtv W Q ρ l|) := by
        rw [if_neg hkl, sub_zero, hinner k l, abs_mul, abs_mul,
          abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν k))),
          abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν l)))]
      rw [habs]
      have hfac : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹
          = (Real.sqrt (ν k * ν l))⁻¹ := by
        rw [Real.sqrt_mul (hν k).le, mul_inv]
      rcases eq_or_ne (ρ k) (ρ l) with hρkl | hρkl
      · -- tied outliers: a `G₀²` form
        have htie : xtv W Q ρ k ⬝ᵥ xtv W Q ρ l
            = cform2 W (ρ k) (fun i => Q i k) (fun i => Q i l) := by
          rw [xtv_dot_xtv hW k l,
            show resolv W (ρ l) = resolv W (ρ k) by rw [hρkl]]
          rfl
        have hb := hE2 k l
        rw [if_neg hkl, sub_zero] at hb
        have hgt : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ ≤ gramC ρ ν := by
          have h1 := gramTerm_le_gramC ρ ν k l
          rw [if_pos hρkl, one_mul] at h1
          rw [hfac]
          exact h1
        calc (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * |xtv W Q ρ k ⬝ᵥ xtv W Q ρ l|)
            ≤ (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * η) := by
              rw [htie]
              exact mul_le_mul_of_nonneg_left
                (mul_le_mul_of_nonneg_left hb (by positivity)) (by positivity)
          _ = (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * η := by ring
          _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0.le
      · -- distinct outliers: the resolvent identity `R4.resolv_sub_resolv`
        have hne2 : ρ l - ρ k ≠ 0 := sub_ne_zero.mpr (Ne.symm hρkl)
        have hrs := resolv_sub_resolv (W₀ := W) (z₁ := ρ k) (z₂ := ρ l)
          (isUnit_det_sub hW (hτz k)) (isUnit_det_sub hW (hτz l))
        have hprod : resolv W (ρ k) * resolv W (ρ l)
            = (ρ l - ρ k)⁻¹ • (resolv W (ρ l) - resolv W (ρ k)) := by
          rw [hrs, smul_smul, inv_mul_cancel₀ hne2, one_smul]
        have hval : xtv W Q ρ k ⬝ᵥ xtv W Q ρ l
            = (ρ l - ρ k)⁻¹ * (cform W (ρ l) (fun i => Q i k) (fun i => Q i l)
              - cform W (ρ k) (fun i => Q i k) (fun i => Q i l)) := by
          rw [xtv_dot_xtv hW k l, hprod, Matrix.smul_mulVec, dotProduct_smul,
            smul_eq_mul, Matrix.sub_mulVec, dotProduct_sub]
          rfl
        have hb1 := hE1 k l
        rw [if_neg hkl, sub_zero] at hb1
        have hb2 := hE1 l k
        rw [if_neg (Ne.symm hkl), sub_zero] at hb2
        have hb2' : |cform W (ρ k) (fun i => Q i k) (fun i => Q i l)| ≤ η := by
          rw [FormsR.cform_comm hW]
          exact hb2
        have hdiff : |cform W (ρ l) (fun i => Q i k) (fun i => Q i l)
            - cform W (ρ k) (fun i => Q i k) (fun i => Q i l)| ≤ η + η := by
          calc |cform W (ρ l) (fun i => Q i k) (fun i => Q i l)
              - cform W (ρ k) (fun i => Q i k) (fun i => Q i l)|
              ≤ |cform W (ρ l) (fun i => Q i k) (fun i => Q i l)|
                + |cform W (ρ k) (fun i => Q i k) (fun i => Q i l)| := abs_sub _ _
            _ ≤ η + η := add_le_add hb1 hb2'
        have hdotb : |xtv W Q ρ k ⬝ᵥ xtv W Q ρ l| ≤ 2 / |ρ k - ρ l| * η := by
          rw [hval, abs_mul, abs_inv, abs_sub_comm (ρ l) (ρ k)]
          calc |ρ k - ρ l|⁻¹ * |cform W (ρ l) (fun i => Q i k) (fun i => Q i l)
              - cform W (ρ k) (fun i => Q i k) (fun i => Q i l)|
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
        calc (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * |xtv W Q ρ k ⬝ᵥ xtv W Q ρ l|)
            ≤ (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * (2 / |ρ k - ρ l| * η)) :=
              mul_le_mul_of_nonneg_left
                (mul_le_mul_of_nonneg_left hdotb (by positivity)) (by positivity)
          _ = (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * (2 / |ρ k - ρ l|) * η := by
              ring
          _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0.le
  -- 4. the frame approximation (plan 1.3 items 3 to 5, through `Frame`)
  have hδ0 : (0 : ℝ) ≤ resC ρ ν * η := mul_nonneg (resC_nonneg ρ ν) hη0.le
  have hδ'0 : (0 : ℝ) ≤ gramC ρ ν * η := mul_nonneg (gramC_nonneg ρ ν) hη0.le
  have hframe := Frame.specProjTop_frame_approx_of_split (δ := resC ρ ν * η)
    (δ' := gramC ρ ν * η) hW hS hSeq hrp hτ hmg hδ0 hδ'0 hρ hresid hgrame hsmall hδm
    (WithLp.toLp 2 v)
  have hx1 : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot, WithLp.ofLp_toLp, hv]
  rw [hx1, mul_one] at hframe
  -- 5. the overlaps against the target (plan 1.3 item 5)
  have hover : ∀ k, ⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ
      = (Real.sqrt (ν k))⁻¹ * cform W (ρ k) v (fun i => Q i k) := by
    intro k
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W Q ρ ν k) = (Real.sqrt (ν k))⁻¹ • xtv W Q ρ k := rfl
    rw [h1, WithLp.ofLp_toLp, smul_dotProduct, smul_eq_mul]
    congr 1
    rw [dotProduct_comm]
    rfl
  have hLk : ∀ k : Fin r, |if j = k then L else 0| ≤ |L| := by
    intro k
    rcases eq_or_ne j k with rfl | hne'
    · rw [if_pos rfl]
    · rw [if_neg hne', abs_zero]
      exact abs_nonneg L
  have hterm : ∀ k, |⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
      - (if j = k then L else 0) ^ 2 / ν k| ≤ (1 + 2 * |L|) * (ν k)⁻¹ * η := by
    intro k
    have hb := hE3 k
    have h1 : ⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        = (ν k)⁻¹ * cform W (ρ k) v (fun i => Q i k) ^ 2 := by
      rw [hover k, mul_pow, inv_pow, Real.sq_sqrt (hν k).le]
    have h2 : |cform W (ρ k) v (fun i => Q i k) ^ 2 - (if j = k then L else 0) ^ 2|
        ≤ (1 + 2 * |L|) * η := by
      have habs2 : |cform W (ρ k) v (fun i => Q i k) ^ 2 - (if j = k then L else 0) ^ 2|
          = |cform W (ρ k) v (fun i => Q i k) - (if j = k then L else 0)|
            * |cform W (ρ k) v (fun i => Q i k) + (if j = k then L else 0)| := by
        rw [← abs_mul]
        congr 1
        ring
      have h3 : |cform W (ρ k) v (fun i => Q i k) + (if j = k then L else 0)|
          ≤ η + 2 * |L| := by
        have h4 : cform W (ρ k) v (fun i => Q i k) + (if j = k then L else 0)
            = (cform W (ρ k) v (fun i => Q i k) - (if j = k then L else 0))
              + 2 * (if j = k then L else 0) := by ring
        rw [h4]
        calc |(cform W (ρ k) v (fun i => Q i k) - (if j = k then L else 0))
            + 2 * (if j = k then L else 0)|
            ≤ |cform W (ρ k) v (fun i => Q i k) - (if j = k then L else 0)|
              + |2 * (if j = k then L else 0)| := abs_add_le _ _
          _ ≤ η + 2 * |L| := by
              rw [abs_mul, abs_two]
              exact add_le_add hb (by linarith [hLk k])
      rw [habs2]
      calc |cform W (ρ k) v (fun i => Q i k) - (if j = k then L else 0)|
          * |cform W (ρ k) v (fun i => Q i k) + (if j = k then L else 0)|
          ≤ η * (η + 2 * |L|) := mul_le_mul hb h3 (abs_nonneg _) hη0.le
        _ ≤ (1 + 2 * |L|) * η := by nlinarith [abs_nonneg L]
    have hgoal : |⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        - (if j = k then L else 0) ^ 2 / ν k|
        = (ν k)⁻¹ * |cform W (ρ k) v (fun i => Q i k) ^ 2
          - (if j = k then L else 0) ^ 2| := by
      rw [h1, show (if j = k then L else 0) ^ 2 / ν k
          = (ν k)⁻¹ * (if j = k then L else 0) ^ 2 by rw [div_eq_mul_inv, mul_comm],
        ← mul_sub, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
    rw [hgoal]
    calc (ν k)⁻¹ * |cform W (ρ k) v (fun i => Q i k) ^ 2 - (if j = k then L else 0) ^ 2|
        ≤ (ν k)⁻¹ * ((1 + 2 * |L|) * η) :=
          mul_le_mul_of_nonneg_left h2 (inv_nonneg.mpr (hν k).le)
      _ = (1 + 2 * |L|) * (ν k)⁻¹ * η := by ring
  have hLksum : ∑ k, (if j = k then L else 0) ^ 2 / ν k = L ^ 2 / ν j := by
    have hcongr : ∀ k : Fin r, (if j = k then L else 0) ^ 2 / ν k
        = if j = k then L ^ 2 / ν k else 0 := by
      intro k
      rcases eq_or_ne j k with rfl | hne'
      · rw [if_pos rfl, if_pos rfl]
      · rw [if_neg hne', if_neg hne']
        simp
    rw [Finset.sum_congr rfl fun k _ => hcongr k, Finset.sum_ite_eq]
    simp
  have hsum : |∑ k, ⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - L ^ 2 / ν j|
      ≤ (1 + 2 * |L|) * (∑ k, (ν k)⁻¹) * η := by
    rw [← hLksum, ← Finset.sum_sub_distrib]
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ k, |⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        - (if j = k then L else 0) ^ 2 / ν k|
        ≤ ∑ k, (1 + 2 * |L|) * (ν k)⁻¹ * η := Finset.sum_le_sum fun k _ => hterm k
      _ = (1 + 2 * |L|) * (∑ k, (ν k)⁻¹) * η := by
          rw [← Finset.sum_mul, ← Finset.mul_sum]
  -- 6. assemble
  calc |‖specProjTop S hS r (WithLp.toLp 2 v)‖ ^ 2 - L ^ 2 / ν j|
      ≤ |‖specProjTop S hS r (WithLp.toLp 2 v)‖ ^ 2
          - ∑ k, ⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2|
        + |∑ k, ⟪yhatv W Q ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - L ^ 2 / ν j| :=
        abs_sub_le _ _ _
    _ ≤ 5 * r * (gramC ρ ν * η + resC ρ ν * η / mg)
        + (1 + 2 * |L|) * (∑ k, (ν k)⁻¹) * η := add_le_add hframe hsum

end Core

/-- Squeeze with an eventual inclusion: the side condition `r ≤ d N` holds only for
large `N`. -/
theorem tendsto_measure_zero_of_eventually_subset {Ω : ℕ → Type*}
    [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {s t : ∀ N, Set (Ω N)}
    (hst : ∀ᶠ N in atTop, s N ⊆ t N)
    (ht : Tendsto (fun N => μ N (t N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (s N)) atTop (𝓝 0) :=
  tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds ht
    (Filter.Eventually.of_forall fun _ => zero_le)
    (hst.mono fun _ h => measure_mono h)

end OutliersR

end StackedSVD

/-! ### 4. The model layer and the supercritical discharge -/

namespace StackedSVD

/-! ### 4a. The shared stack interface: the supercritical discharge

The deterministic core of section 3 is model-free already. This section runs it on
`RankRStack` (`RankR/RMT/Stack.lean`), so the discharge belongs to the interface and not to
the `r_i = 1` model. Section 4b is this theorem at `UnalignedModel.toStack`. -/

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The `vᵀ G₀ Q` limit (plan 1.3 item 5) on the shared interface: entry `k` of
`v_jᵀ G₀(z) Q` tends to `δ_jk √λ_k m(z)`, from the `vv` and `vg` fields. -/
theorem tendstoInProb_cform_vmat (s : RankRStack μ ns d r)
    {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) z
        (WithLp.ofLp (s.spikeVec j N)) (fun l => s.qmatR N ω (U N) l k))
      (if j = k then Real.sqrt (s.coreEig k) * MP.m c z else 0) := by
  have hcol : ∀ (N : ℕ) (ω : Ω N), (fun l => s.qmatR N ω (U N) l k)
      = Real.sqrt (s.coreEig k) • WithLp.ofLp (s.spikeVec k N)
        + fun l => ((s.E N ω)ᵀ * U N) l k := by
    intro N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, s.qmatR_apply N ω (U N) l k]
  have hfun : (fun N (ω : Ω N) => R4.cform (s.rankRW0 N ω (U N)) z
      (WithLp.ofLp (s.spikeVec j N)) (fun l => s.qmatR N ω (U N) l k))
      = fun N ω => Real.sqrt (s.coreEig k) * R4.cform (s.rankRW0 N ω (U N)) z
          (WithLp.ofLp (s.spikeVec j N)) (WithLp.ofLp (s.spikeVec k N))
        + R4.cform (s.rankRW0 N ω (U N)) z (WithLp.ofLp (s.spikeVec j N))
          (fun l => ((s.E N ω)ᵀ * U N) l k) := by
    funext N ω
    rw [hcol N ω, R4.cform_add_right, R4.cform_smul_right]
  rw [hfun]
  have hcomb := ((h.vv j k z hz).const_mul (Real.sqrt (s.coreEig k))).add (h.vg j k z hz)
  refine FormsR.tendstoInProb_congr_limit ?_ hcomb
  rcases eq_or_ne j k with rfl | hjk
  · simp
  · simp [hjk]

/-- **Task U4, the supercritical discharge** (plan section 3, U4 row; Route P) on the shared
interface. For a `RankRStack` with Gaussian noise, aspect ratio `c` and every spike
supercritical (`c < λ_j(C)²` for all `j`), the top-`r` eigenspace of the Gram matrix overlaps
each spike direction `V q_j` by `betaSq (√λ_j) c` in probability. That is exactly the `align`
field of `UnalignedModel.SubspaceLaw`, and of any general-`r_i` twin of it.

Ties among the `λ_j(C)` are allowed. The side condition `hn`/`hp` (that is, `r < ns N` at
every `N`) is the U2 side condition; its removal by the `TailShift` device (decision D11),
and the mixed case with subcritical spikes, are task U7. -/
theorem align_of_gaussian_supercritical [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N)
    (hsup : ∀ j, c < s.coreEig j ^ 2) :
    ∀ j : Fin r, TendstoInProb μ
      (fun N ω => ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
        (s.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (s.coreEig j)) c) := by
  classical
  obtain ⟨U, hU, hsig, hgram, h⟩ := s.resolventLimitsR_of_gaussian hG hc hdtop hns hn hp
  set ρv : Fin r → ℝ := fun k => rhoSq (Real.sqrt (s.coreEig k)) c with hρv
  set νv : Fin r → ℝ := fun k => (s.coreEig k + 1) * MP.mDeriv c (ρv k) with hνv
  have hlam0 : ∀ k, 0 ≤ s.coreEig k := s.coreEig_nonneg
  have hbρ : ∀ k, bulkEdge c < ρv k := fun k =>
    OutliersR.bulkEdge_lt_rho hc (hlam0 k) (hsup k)
  have hm1 : ∀ k, (s.coreEig k + 1) * MP.m c (ρv k) = -1 := fun k =>
    OutliersR.mul_m_rho_eq_neg_one hc (hlam0 k) (hsup k)
  have hνpos : ∀ k, 0 < νv k := fun k => OutliersR.nu_pos hc (hlam0 k) (hsup k)
  intro j
  have hne : (Finset.univ : Finset (Fin r)).Nonempty := ⟨j, Finset.mem_univ j⟩
  set ρmin : ℝ := Finset.univ.inf' hne ρv with hρmin
  have hbmin : bulkEdge c < ρmin := (Finset.lt_inf'_iff hne).mpr fun k _ => hbρ k
  set mg : ℝ := (ρmin - bulkEdge c) / 2 with hmgdef
  have hmg : 0 < mg := by rw [hmgdef]; linarith
  have hτρ : ∀ k, bulkEdge c + mg + mg ≤ ρv k := by
    intro k
    have h1 : ρmin ≤ ρv k := Finset.inf'_le _ (Finset.mem_univ k)
    rw [hmgdef]
    linarith
  set Lval : ℝ := Real.sqrt (s.coreEig j) * MP.m c (ρv j) with hLval
  have hLtarget : Lval ^ 2 / νv j = betaSq (Real.sqrt (s.coreEig j)) c :=
    OutliersR.sq_div_nu_eq_betaSq hc (hlam0 j) (hsup j)
  have hvdot : ∀ N, WithLp.ofLp (s.spikeVec j N) ⬝ᵥ WithLp.ofLp (s.spikeVec j N) = 1 := by
    intro N
    rw [← Frame.inner_eq_dot]
    have h1 := s.inner_spikeVec N j j
    rwa [if_pos rfl] at h1
  intro ε hε
  -- the constants and the accuracy `η`
  have hνsum0 : (0 : ℝ) ≤ ∑ k, (νv k)⁻¹ :=
    Finset.sum_nonneg fun k _ => inv_nonneg.mpr (hνpos k).le
  have hL0 : (0 : ℝ) ≤ 1 + 2 * |Lval| := by positivity
  have hCR0 : (0 : ℝ) ≤ OutliersR.resC ρv νv := OutliersR.resC_nonneg ρv νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  have hX0 : (0 : ℝ) ≤ (r : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resC ρv νv / mg) :=
    mul_nonneg (Nat.cast_nonneg r) (add_nonneg hCG0 (div_nonneg hCR0 hmg.le))
  have hB0 : (0 : ℝ) ≤ 5 * ((r : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resC ρv νv / mg)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ := by
    have h1 := mul_nonneg hL0 hνsum0
    linarith
  set η : ℝ := min (min 1 (mg / (OutliersR.resC ρv νv + 1)))
      (min (1 / (2 * ((r : ℝ) * (OutliersR.gramC ρv νv
          + OutliersR.resC ρv νv / mg) + 1)))
        (ε / (2 * (5 * ((r : ℝ) * (OutliersR.gramC ρv νv
            + OutliersR.resC ρv νv / mg))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    refine lt_min (lt_min one_pos (div_pos hmg (by linarith))) (lt_min ?_ ?_)
    · exact div_pos one_pos (by linarith)
    · exact div_pos hε (by linarith)
  have hη1 : η ≤ 1 := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_left _ _)
  have hηm : η ≤ mg / (OutliersR.resC ρv νv + 1) := by
    rw [hηdef]
    exact le_trans (min_le_left _ _) (min_le_right _ _)
  have hηX : η ≤ 1 / (2 * ((r : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resC ρv νv / mg) + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hηε : η ≤ ε / (2 * (5 * ((r : ℝ) * (OutliersR.gramC ρv νv
      + OutliersR.resC ρv νv / mg)) + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹ + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (min_le_right _ _)
  have hδm : OutliersR.resC ρv νv * η ≤ mg :=
    OutliersR.mul_le_of_le_div_add_one hCR0 hηm hmg.le
  have hsmall : (r : ℝ) * (OutliersR.gramC ρv νv * η + OutliersR.resC ρv νv * η / mg)
      ≤ 1 / 2 := by
    have hXe : (r : ℝ) * (OutliersR.gramC ρv νv * η + OutliersR.resC ρv νv * η / mg)
        = ((r : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resC ρv νv / mg)) * η := by
      ring
    rw [hXe]
    exact le_trans (mul_le_mul_of_nonneg_left hηX hX0) (OutliersR.mul_inv_le_half hX0)
  have hfinal : 5 * (r : ℝ) * (OutliersR.gramC ρv νv * η
        + OutliersR.resC ρv νv * η / mg)
      + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η ≤ ε / 2 := by
    have hKe : 5 * (r : ℝ) * (OutliersR.gramC ρv νv * η
          + OutliersR.resC ρv νv * η / mg)
        + (1 + 2 * |Lval|) * (∑ k, (νv k)⁻¹) * η
        = (5 * ((r : ℝ) * (OutliersR.gramC ρv νv + OutliersR.resC ρv νv / mg))
          + (1 + 2 * |Lval|) * ∑ k, (νv k)⁻¹) * η := by
      ring
    rw [hKe]
    exact le_trans (mul_le_mul_of_nonneg_left hηε hB0)
      (OutliersR.mul_div_le_half hB0 hε.le)
  -- the four bad families
  have hedgeC : Tendsto (fun N => μ N {ω | lamMax (s.rankRW0 N ω (U N))
      (s.isHermitian_rankRW0 N ω (U N)) ≤ bulkEdge c + mg}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + mg)).nullMeasurableSet)
      (h.edge mg hmg)
  have hE1T : ∀ q : Fin r × Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
        (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2)
        - (if q.1 = q.2 then -1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := s.tendstoInProb_cform_qmatR h q.1 q.2 (hbρ q.2)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
        (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2))
        (if q.1 = q.2 then -1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 q.2 with heq | hne'
      · rw [if_pos heq, if_pos heq, ← heq]
        exact hm1 q.1
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  have hE2T : ∀ q : Fin r × Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
        (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2)
        - (if q.1 = q.2 then νv q.1 else 0)|}) atTop (𝓝 0) := by
    intro q
    have h1 := s.tendstoInProb_cform2_qmatR h q.1 q.2 (hbρ q.1)
    have h2 : TendstoInProb μ (fun N ω => R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
        (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2))
        (if q.1 = q.2 then νv q.1 else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne q.1 q.2 with heq | hne'
      · rw [if_pos heq, if_pos heq, ← heq, hνv]
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  have hE3T : ∀ k : Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv k)
        (WithLp.ofLp (s.spikeVec j N)) (fun l => s.qmatR N ω (U N) l k)
        - (if j = k then Lval else 0)|}) atTop (𝓝 0) := by
    intro k
    have h1 := s.tendstoInProb_cform_vmat h j k (hbρ k)
    have h2 : TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) (ρv k)
        (WithLp.ofLp (s.spikeVec j N)) (fun l => s.qmatR N ω (U N) l k))
        (if j = k then Lval else 0) := by
      refine FormsR.tendstoInProb_congr_limit ?_ h1
      rcases eq_or_ne j k with rfl | hne'
      · rw [if_pos rfl]
      · rw [if_neg hne', if_neg hne']
    exact h2 η hη0
  -- assemble: eventual inclusion into the union of the bad families
  refine OutliersR.tendsto_measure_zero_of_eventually_subset (t := fun N =>
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
          ≤ bulkEdge c + mg}ᶜ
        ∪ ((⋃ q : Fin r × Fin r, {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv q.2)
              (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2)
              - (if q.1 = q.2 then -1 else 0)|})
          ∪ ((⋃ q : Fin r × Fin r, {ω | η ≤ |R4.cform2 (s.rankRW0 N ω (U N)) (ρv q.1)
                (fun l => s.qmatR N ω (U N) l q.1) (fun l => s.qmatR N ω (U N) l q.2)
                - (if q.1 = q.2 then νv q.1 else 0)|})
            ∪ (⋃ k : Fin r, {ω | η ≤ |R4.cform (s.rankRW0 N ω (U N)) (ρv k)
                (WithLp.ofLp (s.spikeVec j N)) (fun l => s.qmatR N ω (U N) l k)
                - (if j = k then Lval else 0)|})))) ?_ ?_
  · filter_upwards [hdtop.eventually_ge_atTop r] with N hrd
    intro ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le] at hbad
    obtain ⟨hedge, hb1, hb2, hb3⟩ := hbad
    have hdet := OutliersR.align_det (Q := s.qmatR N ω (U N)) hrd j.pos
      (s.isHermitian_rankRW0 N ω (U N)) (s.isHermitian_gram N ω)
      (fun a => s.eigenvalues_rankRW0_nonneg N ω (U N) a) (hgram N ω) hmg
      (not_lt.mp hedge) hτρ
      hνpos (hvdot N) j hη0 hη1
      (fun jj kk => (hb1 (jj, kk)).le) (fun kk l => (hb2 (kk, l)).le)
      (fun k => (hb3 k).le) hsmall hδm
    rw [WithLp.toLp_ofLp] at hdet
    have hb := le_trans hdet hfinal
    rw [hLtarget] at hb
    have hω' : ε ≤ |‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
        (s.spikeVec j N)‖ ^ 2 - betaSq (Real.sqrt (s.coreEig j)) c| := hω
    linarith
  · exact tendsto_measure_zero_union hedgeC
      (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
          (tendsto_measure_zero_iUnion hE3T)))


end RankRStack

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The block `W₀` of the split is positive semidefinite. -/
theorem eigenvalues_rankRW0_nonneg (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (a : Fin (d N)) :
    0 ≤ (m.isHermitian_rankRW0 N ω U).eigenvalues a := by
  have hpsd : (m.rankRW0 N ω U).PosSemidef := by
    change ((m.stackEperp N ω U)ᵀ * m.stackEperp N ω U).PosSemidef
    simpa using Matrix.posSemidef_conjTranspose_mul_self (m.stackEperp N ω U)
  exact hpsd.eigenvalues_nonneg a

/-- The `vᵀ G₀ Q` limit (plan 1.3 item 5): entry `k` of `v_jᵀ G₀(z) Q` tends to
`δ_jk √λ_k m(z)`, from the `vv` and `vg` fields of the interface. -/
theorem tendstoInProb_cform_vmat (m : UnalignedModel μ M n d r)
    {U : (N : ℕ) → Matrix (Fin (∑ i, n i N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => m.rankRW0 N ω (U N))
      (fun N ω => m.isHermitian_rankRW0 N ω (U N)) m.spikeVec
      (fun k N ω => fun l => ((m.stackE N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform (m.rankRW0 N ω (U N)) z
        (WithLp.ofLp (m.spikeVec j N)) (fun l => m.qmatR N ω (U N) l k))
      (if j = k then Real.sqrt (m.coreEig k) * MP.m c z else 0) := by
  have hcol : ∀ (N : ℕ) (ω : Ω N), (fun l => m.qmatR N ω (U N) l k)
      = Real.sqrt (m.coreEig k) • WithLp.ofLp (m.spikeVec k N)
        + fun l => ((m.stackE N ω)ᵀ * U N) l k := by
    intro N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, m.qmatR_apply N ω (U N) l k]
  have hfun : (fun N (ω : Ω N) => R4.cform (m.rankRW0 N ω (U N)) z
      (WithLp.ofLp (m.spikeVec j N)) (fun l => m.qmatR N ω (U N) l k))
      = fun N ω => Real.sqrt (m.coreEig k) * R4.cform (m.rankRW0 N ω (U N)) z
          (WithLp.ofLp (m.spikeVec j N)) (WithLp.ofLp (m.spikeVec k N))
        + R4.cform (m.rankRW0 N ω (U N)) z (WithLp.ofLp (m.spikeVec j N))
          (fun l => ((m.stackE N ω)ᵀ * U N) l k) := by
    funext N ω
    rw [hcol N ω, R4.cform_add_right, R4.cform_smul_right]
  rw [hfun]
  have hcomb := ((h.vv j k z hz).const_mul (Real.sqrt (m.coreEig k))).add (h.vg j k z hz)
  refine FormsR.tendstoInProb_congr_limit ?_ hcomb
  rcases eq_or_ne j k with rfl | hjk
  · simp
  · simp [hjk]

/-- **Task U4, the supercritical discharge** (plan section 3, U4 row; Route P). For an
`UnalignedModel` with joint Gaussian noise, per-table regimes `cc` with `0 < ∑ cc i`, and
every spike supercritical (`∑ cc i < λ_j(C)²` for all `j`), the rank-`r` spiked law
`SubspaceLaw` of the unit-weight stack holds at aspect ratio `∑ cc i`. This is the exact
hypothesis of `prop_stacksvd_subspace`, so the two compose with no glue.

Ties among the `λ_j(C)` are allowed. The side condition `hn`/`hp` (that is,
`r < ∑ i, n i N` at every `N`) is the U2 side condition; its removal by the `TailShift`
device (decision D11), and the mixed case with subcritical spikes, are task U7. -/
theorem subspaceLaw_of_gaussian_supercritical [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i)
    {pp : ℕ → ℕ} (hn : ∀ N, ∑ i, n i N = r + pp N) (hp : ∀ N, 0 < pp N)
    (hsup : ∀ j, ∑ i, cc i < m.coreEig j ^ 2) :
    m.SubspaceLaw (∑ i, cc i) :=
  ⟨m.toStack.align_of_gaussian_supercritical (m.gaussianNoise_toStack hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTable.stack_regime cc hreg).2.2 hn hp hsup⟩

end UnalignedModel

end StackedSVD
