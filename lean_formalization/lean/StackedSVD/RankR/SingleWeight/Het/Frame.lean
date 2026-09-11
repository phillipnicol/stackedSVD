/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.OutliersG

/-!
# Track G, unit G2: the deterministic core with frame columns `Q * Z`

Section G2 of `notes/archive/trackG_plan.md` (2026-09-05). Mirrors `OutliersR.residual_eq`
(`RankR/RMT/Outliers.lean:200`), `residual_bound_sub`, `gram_bound_sub` and `align_detG`
(`RankR/RMT/OutliersG.lean`), replacing the sub-frame selection `Q.submatrix id f` (an
injection `f : Fin s → Fin r` that picks `s` of the `r` columns of `Q`) by a general real
linear combination `Q * Z` for `Z : Matrix (Fin r) (Fin s) ℝ`. Model free: no
`UnalignedModelR`, no measure. The count lemma (part 4 of the plan's proof sketch, the `s`
eigenvalues above `τ`) is a separate unit and is not proved here.

Risk 2 of the plan: the off-diagonal `cform` limits are not `0`, unlike the aligned case.
`check_frame.py` shows `z_kᵀ F(ρ_l) z_l → -(ZᵀZ)_{kl}`, not `0` (on the check instance
`-(ZᵀZ)_{01} = 0.0026...`, matched by the two resolvent points to `1e-16`), while the
*difference* across the two resolvent points `z_kᵀ [F(ρ_l) - F(ρ_k)] z_l → 0` (`4.9e-17` on
the same instance). This is why `gram_bound_Z` carries the target `-Z l k` in `hE1`
(matching `residual_bound_Z`'s own `hE1`) and needs only a *diagonal* `cform2` accuracy
`hE2`: the tie branch of the aligned proof (`ρ k = ρ l`) is replaced by
`hρinj : Function.Injective ρ`, which makes every off-diagonal pair `k ≠ l` a distinct-`ρ`
pair, so `hE2` is never needed off the diagonal.

Contents:

1. `residual_eq_vec`: `residual_eq` with the indicator column `e_k` replaced by a general
   `y : Fin r → ℝ`, so `(S - ζ) G₀(ζ) Q y = Q (y + Qᵀ G₀(ζ) Q y)` entrywise. The `hqk` step
   of `residual_eq` is not needed: `h1`'s conclusion is already `Q *ᵥ y`, not a bare column,
   so the two summands combine directly under `Matrix.mulVec_add`.
2. Two small helpers. `mulVec_col_eq` is `(Q * Z) e_k = Q *ᵥ (Z ·, k)`, by `rfl`.
   `cform_mulVec_left` is the bilinearity of `cform` in its first argument along a column
   combination, `cform W z (Q *ᵥ u) x = ∑ j, u j * cform W z (Q_j) x`, proved through
   `Matrix.dotProduct_mulVec` and `dotProduct_comm` rather than raw `Finset.sum` algebra.
3. `residual_bound_Z`: `residual_bound_sub` with `Q.submatrix id f` replaced by `Q * Z`.
   The bounded vector loses its indicator branch entirely, since `hE1`'s target is already
   `-Z l k` and not an if-then-else, so there is no case split on `l = f k`.
4. `gram_bound_Z`: `gram_bound_sub` with the tie branch deleted (`hρinj` makes every
   off-diagonal pair a distinct-`ρ` pair) and the distinct branch's two `cform` values
   expanded by `cform_mulVec_left` against the *same* target `-(ZᵀZ)_{kl}`, so they cancel
   in the difference exactly as the two zeros do in the aligned proof, up to `2 Zb η`
   instead of `2 η`. The resulting per-pair bound is `gramC ρ ν * (1 + Zb) * η`: the
   diagonal branch keeps the aligned constant `gramC ρ ν * η` (factor `1`), the
   off-diagonal branch scales it by `Zb` (factor `Zb`), and `1 + Zb` dominates both since
   `hZb0 : 0 ≤ Zb`.
5. `align_detZ`: `align_detG` with `Q.submatrix id f` replaced by `Q * Z`, `hf` replaced by
   `hρinj`, and the Gram constant `gramC ρ ν * η` replaced by `gramC ρ ν * (1 + Zb) * η`.
   The consumer shape: `RankR/Het/Outliers.lean:922` calls `align_detG` directly (the
   general `t`/`Lb` form), not `align_detG_mem`, so only that shape is mirrored here; no
   `align_detZ_mem` or `align_detZ_col` is written.

No `sorry`.
-/

open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace OutliersR

open R4

/-! ### 1. The residual identity at a general frame vector -/

section Deterministic

variable {p r : ℕ} {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}

/-- Plan 1.3 item 1 with a general frame vector `y` in place of the indicator column `e_k`:
the mirror of `residual_eq`. `(S - ζ) G₀(ζ) Q y = Q (y + Qᵀ G₀(ζ) Q y)`, entrywise. -/
theorem residual_eq_vec (hW : W.IsHermitian) (hSeq : S = W + Q * Qᵀ) {ζ : ℝ}
    (hz : lamMax W hW < ζ) (y : Fin r → ℝ) :
    S *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y)) - ζ • (resolv W ζ *ᵥ (Q *ᵥ y))
      = Q *ᵥ fun l => y l + cform W ζ (fun i => Q i l) (Q *ᵥ y) := by
  have h1 : (W - ζ • (1 : Matrix (Fin p) (Fin p) ℝ)) *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y))
      = Q *ᵥ y := by
    rw [Matrix.mulVec_mulVec, mul_resolv hW hz, Matrix.one_mulVec]
  have h2 : Qᵀ *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y))
      = fun l => cform W ζ (fun i => Q i l) (Q *ᵥ y) := by
    funext l
    simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, cform]
  have hsplit : S *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y)) - ζ • (resolv W ζ *ᵥ (Q *ᵥ y))
      = (W - ζ • (1 : Matrix (Fin p) (Fin p) ℝ)) *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y))
        + (Q * Qᵀ) *ᵥ (resolv W ζ *ᵥ (Q *ᵥ y)) := by
    rw [hSeq, Matrix.add_mulVec, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec]
    abel
  rw [hsplit, h1, ← Matrix.mulVec_mulVec, h2, ← Matrix.mulVec_add]
  congr 1

end Deterministic

/-! ### 2. The frame with columns `Q * Z` -/

section Core

variable {p r s : ℕ}

/-- The `k`-th column of `Q * Z` is `Q` applied to the `k`-th column of `Z`. -/
private theorem mulVec_col_eq (Q : Matrix (Fin p) (Fin r) ℝ) (Z : Matrix (Fin r) (Fin s) ℝ)
    (k : Fin s) : (fun i => (Q * Z) i k) = Q *ᵥ fun l => Z l k := rfl

/-- Bilinearity of `cform` in its first argument along a column combination: this is the
`∑ j, Z j k · cform W (ρ l) (Q_j) (·)` expansion that risk 2 needs, proved through
`Matrix.dotProduct_mulVec` and `dotProduct_comm` rather than raw `Finset.sum` algebra. -/
private theorem cform_mulVec_left (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Q : Matrix (Fin p) (Fin r) ℝ) (u : Fin r → ℝ) (x : Fin p → ℝ) :
    cform W z (Q *ᵥ u) x = ∑ j, u j * cform W z (fun i => Q i j) x := by
  have h1 : cform W z (Q *ᵥ u) x = ((resolv W z *ᵥ x) ᵥ* Q) ⬝ᵥ u :=
    (dotProduct_comm (Q *ᵥ u) (resolv W z *ᵥ x)).trans (Matrix.dotProduct_mulVec _ _ _)
  have h2 : (resolv W z *ᵥ x) ᵥ* Q = fun j => cform W z (fun i => Q i j) x :=
    funext fun j => dotProduct_comm (resolv W z *ᵥ x) (fun i => Q i j)
  rw [h1, h2]
  exact Finset.sum_congr rfl fun j _ => mul_comm _ _

/-- The residual bound with frame columns `Q * Z` (the mirror of `residual_bound_sub`).
`hE1`'s target is already `-Z l k`, not an if-then-else, so there is no case split. -/
theorem residual_bound_Z {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    {Z : Matrix (Fin r) (Fin s) ℝ} (hW : W.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {ρ ν : Fin s → ℝ} (hν : ∀ k, 0 < ν k) (hz : ∀ k, lamMax W hW < ρ k)
    {Cq η : ℝ} (hCq : 0 ≤ Cq) (hη0 : 0 ≤ η) (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k) + Z l k| ≤ η) (k : Fin s) :
    ‖toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k‖ ≤ resCG Cq r ν * η := by
  classical
  have hsfmul : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν k))⁻¹ = (ν k)⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (hν k).le]
  have hre : S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k
      = Q *ᵥ fun l => Z l k + cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k) :=
    residual_eq_vec hW hSeq (hz k) (fun l => Z l k)
  have hsv : ∀ l : Fin r,
      |Z l k + cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k)| ≤ η :=
    fun l => by rw [add_comm]; exact hE1 l k
  have hsvsum : ∑ l, (Z l k + cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k)) ^ 2
      ≤ (r : ℝ) * η ^ 2 := by
    calc ∑ l, (Z l k + cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k)) ^ 2
        ≤ ∑ _l : Fin r, η ^ 2 := by
          refine Finset.sum_le_sum fun l _ => ?_
          rw [← sq_abs]
          exact pow_le_pow_left₀ (abs_nonneg _) (hsv l) 2
      _ = (r : ℝ) * η ^ 2 := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hcolsum : ∑ l : Fin r, ∑ i, Q i l ^ 2 ≤ (r : ℝ) * Cq := by
    calc ∑ l : Fin r, ∑ i, Q i l ^ 2 ≤ ∑ _l : Fin r, Cq := Finset.sum_le_sum fun l _ => hcol l
      _ = (r : ℝ) * Cq := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hcnn : (0 : ℝ) ≤ (r : ℝ) * Cq := mul_nonneg (Nat.cast_nonneg r) hCq
  have hdotres : (S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k)
        ⬝ᵥ (S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k)
      ≤ ((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2) := by
    rw [hre]
    have h3nn : (0 : ℝ) ≤ ∑ l, (Z l k + cform W (ρ k) (fun i => Q i l)
        (fun i => (Q * Z) i k)) ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
    exact le_trans (dot_mulVec_le_sum_cols Q _) (mul_le_mul hcolsum hsvsum h3nn hcnn)
  have hofLp : WithLp.ofLp (toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k)
      = (Real.sqrt (ν k))⁻¹ • (S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k) := by
    have h1 : WithLp.ofLp (toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k)
        = S *ᵥ ((Real.sqrt (ν k))⁻¹ • xtv W (Q * Z) ρ k)
          - ρ k • ((Real.sqrt (ν k))⁻¹ • xtv W (Q * Z) ρ k) := rfl
    rw [h1, Matrix.mulVec_smul, smul_comm (ρ k) ((Real.sqrt (ν k))⁻¹), ← smul_sub]
  have hns : ‖toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k‖ ^ 2
      = (ν k)⁻¹ * ((S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k)
          ⬝ᵥ (S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k)) := by
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
  have hprod_nn : (0 : ℝ) ≤ ((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2) := mul_nonneg hcnn (by positivity)
  have hsq : ‖toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k‖ ^ 2
      ≤ (resCG Cq r ν * η) ^ 2 := by
    rw [hns]
    calc (ν k)⁻¹ * ((S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k)
          ⬝ᵥ (S *ᵥ xtv W (Q * Z) ρ k - ρ k • xtv W (Q * Z) ρ k))
        ≤ (ν k)⁻¹ * (((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2)) :=
          mul_le_mul_of_nonneg_left hdotres (inv_nonneg.mpr (hν k).le)
      _ ≤ (∑ k', (Real.sqrt (ν k'))⁻¹) ^ 2 * (((r : ℝ) * Cq) * ((r : ℝ) * η ^ 2)) :=
          mul_le_mul_of_nonneg_right hνle hprod_nn
      _ = (resCG Cq r ν * η) ^ 2 := by rw [mul_pow, hresCGsq]; ring
  exact (le_abs_self _).trans (abs_le_of_sq_le_sq hsq (mul_nonneg (resCG_nonneg Cq r ν) hη0))

/-- The Gram bound with frame columns `Q * Z` (the mirror of `gram_bound_sub`). `hρinj`
replaces the tie branch: every off-diagonal pair `k ≠ l` forces `ρ k ≠ ρ l`, so `hE2` is
only ever needed on the diagonal. In the off-diagonal branch the two `cform` values at
`ρ l` and `ρ k` are each within `Zb η` of the *same* target `-(ZᵀZ)_{kl}` (`hA`, `hB`), so
they cancel in the difference up to `2 Zb η`, exactly as the two zeros do in the aligned
proof. -/
theorem gram_bound_Z {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    {Z : Matrix (Fin r) (Fin s) ℝ} (hW : W.IsHermitian) {ρ ν : Fin s → ℝ}
    (hρinj : Function.Injective ρ) (hν : ∀ k, 0 < ν k) (hz : ∀ k, lamMax W hW < ρ k)
    {Zb η : ℝ} (hZb0 : 0 ≤ Zb) (hη0 : 0 ≤ η) (hZb : ∀ k, ∑ l, |Z l k| ≤ Zb)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k) + Z l k| ≤ η)
    (hE2 : ∀ k : Fin s, |cform2 W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i k)
      - ν k| ≤ η) (k l : Fin s) :
    |⟪yhatv W (Q * Z) ρ ν k, yhatv W (Q * Z) ρ ν l⟫_ℝ - if k = l then 1 else 0|
      ≤ gramC ρ ν * (1 + Zb) * η := by
  classical
  have hsfmul : ∀ k', (Real.sqrt (ν k'))⁻¹ * (Real.sqrt (ν k'))⁻¹ = (ν k')⁻¹ :=
    fun k' => by rw [← mul_inv, Real.mul_self_sqrt (hν k').le]
  have hinner : ∀ k' l' : Fin s,
      ⟪yhatv W (Q * Z) ρ ν k', yhatv W (Q * Z) ρ ν l'⟫_ℝ
      = (Real.sqrt (ν k'))⁻¹ * ((Real.sqrt (ν l'))⁻¹
        * (xtv W (Q * Z) ρ k' ⬝ᵥ xtv W (Q * Z) ρ l')) := by
    intro k' l'
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W (Q * Z) ρ ν k')
        = (Real.sqrt (ν k'))⁻¹ • xtv W (Q * Z) ρ k' := rfl
    have h2 : WithLp.ofLp (yhatv W (Q * Z) ρ ν l')
        = (Real.sqrt (ν l'))⁻¹ • xtv W (Q * Z) ρ l' := rfl
    rw [h1, h2, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  rcases eq_or_ne k l with rfl | hkl
  · -- diagonal: `hE2` directly, exactly as the aligned tie branch
    have htie : xtv W (Q * Z) ρ k ⬝ᵥ xtv W (Q * Z) ρ k
        = cform2 W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i k) := by
      rw [xtv_dot_xtv hW k k]; rfl
    have hb := hE2 k
    have hval : ⟪yhatv W (Q * Z) ρ ν k, yhatv W (Q * Z) ρ ν k⟫_ℝ
        - (if k = k then (1 : ℝ) else 0)
        = (ν k)⁻¹ * (cform2 W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i k) - ν k) := by
      rw [if_pos rfl, hinner k k, htie, ← mul_assoc, hsfmul k, mul_sub,
        inv_mul_cancel₀ (hν k).ne']
    rw [hval, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
    have hgt : (ν k)⁻¹ ≤ gramC ρ ν := by
      have h1 := gramTerm_le_gramC ρ ν k k
      rwa [if_pos rfl, one_mul, Real.sqrt_mul_self (hν k).le] at h1
    calc (ν k)⁻¹ * |cform2 W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i k) - ν k|
        ≤ (ν k)⁻¹ * η := mul_le_mul_of_nonneg_left hb (inv_nonneg.mpr (hν k).le)
      _ ≤ gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0
      _ ≤ gramC ρ ν * (1 + Zb) * η := by
          calc gramC ρ ν * η = gramC ρ ν * η * 1 := (mul_one _).symm
            _ ≤ gramC ρ ν * η * (1 + Zb) := by
                apply mul_le_mul_of_nonneg_left _ (mul_nonneg (gramC_nonneg ρ ν) hη0)
                linarith
            _ = gramC ρ ν * (1 + Zb) * η := by ring
  · -- off-diagonal: `hρinj` forces `ρ k ≠ ρ l`
    have hρkl : ρ k ≠ ρ l := fun h => hkl (hρinj h)
    have hne2 : ρ l - ρ k ≠ 0 := sub_ne_zero.mpr (Ne.symm hρkl)
    have hrs := resolv_sub_resolv (W₀ := W) (z₁ := ρ k) (z₂ := ρ l)
      (isUnit_det_sub hW (hz k)) (isUnit_det_sub hW (hz l))
    have hprod : resolv W (ρ k) * resolv W (ρ l)
        = (ρ l - ρ k)⁻¹ • (resolv W (ρ l) - resolv W (ρ k)) := by
      rw [hrs, smul_smul, inv_mul_cancel₀ hne2, one_smul]
    have hval : xtv W (Q * Z) ρ k ⬝ᵥ xtv W (Q * Z) ρ l
        = (ρ l - ρ k)⁻¹ * (cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
          - cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)) := by
      rw [xtv_dot_xtv hW k l, hprod, Matrix.smul_mulVec, dotProduct_smul,
        smul_eq_mul, Matrix.sub_mulVec, dotProduct_sub]
      rfl
    have hA : |cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
        + ∑ j, Z j k * Z j l| ≤ Zb * η := by
      have hexp : cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
          = ∑ j, Z j k * cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l) := by
        rw [mulVec_col_eq Q Z k]
        exact cform_mulVec_left W (ρ l) Q (fun j => Z j k) (fun i => (Q * Z) i l)
      rw [hexp]
      have heq : ∑ j, Z j k * cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l)
          + ∑ j, Z j k * Z j l
          = ∑ j, Z j k * (cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l) + Z j l) := by
        rw [← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl fun j _ => by ring
      rw [heq]
      calc |∑ j, Z j k * (cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l) + Z j l)|
          ≤ ∑ j, |Z j k * (cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l) + Z j l)| :=
            Finset.abs_sum_le_sum_abs _ _
        _ = ∑ j, |Z j k| * |cform W (ρ l) (fun i => Q i j) (fun i => (Q * Z) i l) + Z j l| :=
            Finset.sum_congr rfl fun j _ => abs_mul _ _
        _ ≤ ∑ j, |Z j k| * η :=
            Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (hE1 j l) (abs_nonneg _)
        _ = (∑ j, |Z j k|) * η := by rw [← Finset.sum_mul]
        _ ≤ Zb * η := mul_le_mul_of_nonneg_right (hZb k) hη0
    have hB : |cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
        + ∑ j, Z j k * Z j l| ≤ Zb * η := by
      rw [FormsR.cform_comm hW]
      have hexp : cform W (ρ k) (fun i => (Q * Z) i l) (fun i => (Q * Z) i k)
          = ∑ j, Z j l * cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k) := by
        rw [mulVec_col_eq Q Z l]
        exact cform_mulVec_left W (ρ k) Q (fun j => Z j l) (fun i => (Q * Z) i k)
      rw [hexp]
      have heq : ∑ j, Z j l * cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k)
          + ∑ j, Z j k * Z j l
          = ∑ j, Z j l * (cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k) + Z j k) := by
        rw [← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl fun j _ => by ring
      rw [heq]
      calc |∑ j, Z j l * (cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k) + Z j k)|
          ≤ ∑ j, |Z j l * (cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k) + Z j k)| :=
            Finset.abs_sum_le_sum_abs _ _
        _ = ∑ j, |Z j l| * |cform W (ρ k) (fun i => Q i j) (fun i => (Q * Z) i k) + Z j k| :=
            Finset.sum_congr rfl fun j _ => abs_mul _ _
        _ ≤ ∑ j, |Z j l| * η :=
            Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (hE1 j k) (abs_nonneg _)
        _ = (∑ j, |Z j l|) * η := by rw [← Finset.sum_mul]
        _ ≤ Zb * η := mul_le_mul_of_nonneg_right (hZb l) hη0
    have hdiff : |cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
        - cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)| ≤ 2 * Zb * η := by
      have harg : cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
          - cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
          = (cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
              + ∑ j, Z j k * Z j l)
            - (cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
              + ∑ j, Z j k * Z j l) := by ring
      rw [harg]
      calc |(cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
              + ∑ j, Z j k * Z j l)
            - (cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
              + ∑ j, Z j k * Z j l)|
          ≤ |cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l) + ∑ j, Z j k * Z j l|
            + |cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
              + ∑ j, Z j k * Z j l| := abs_sub _ _
        _ ≤ Zb * η + Zb * η := add_le_add hA hB
        _ = 2 * Zb * η := by ring
    have hdotb : |xtv W (Q * Z) ρ k ⬝ᵥ xtv W (Q * Z) ρ l| ≤ 2 * Zb / |ρ k - ρ l| * η := by
      rw [hval, abs_mul, abs_inv, abs_sub_comm (ρ l) (ρ k)]
      calc |ρ k - ρ l|⁻¹ * |cform W (ρ l) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)
          - cform W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i l)|
          ≤ |ρ k - ρ l|⁻¹ * (2 * Zb * η) := mul_le_mul_of_nonneg_left hdiff (by positivity)
        _ = 2 * Zb / |ρ k - ρ l| * η := by rw [div_eq_mul_inv]; ring
    have hfac : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ = (Real.sqrt (ν k * ν l))⁻¹ := by
      rw [Real.sqrt_mul (hν k).le, mul_inv]
    have hgt : (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * (2 * Zb / |ρ k - ρ l|)
        ≤ Zb * gramC ρ ν := by
      have h1 := gramTerm_le_gramC ρ ν k l
      rw [if_neg hρkl] at h1
      rw [hfac]
      calc (Real.sqrt (ν k * ν l))⁻¹ * (2 * Zb / |ρ k - ρ l|)
          = Zb * (2 / |ρ k - ρ l| * (Real.sqrt (ν k * ν l))⁻¹) := by ring
        _ ≤ Zb * gramC ρ ν := mul_le_mul_of_nonneg_left h1 hZb0
    have habs : |⟪yhatv W (Q * Z) ρ ν k, yhatv W (Q * Z) ρ ν l⟫_ℝ
          - if k = l then (1 : ℝ) else 0|
        = (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹
          * |xtv W (Q * Z) ρ k ⬝ᵥ xtv W (Q * Z) ρ l|) := by
      rw [if_neg hkl, sub_zero, hinner k l, abs_mul, abs_mul,
        abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν k))),
        abs_of_nonneg (inv_nonneg.mpr (Real.sqrt_nonneg (ν l)))]
    rw [habs]
    calc (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹
          * |xtv W (Q * Z) ρ k ⬝ᵥ xtv W (Q * Z) ρ l|)
        ≤ (Real.sqrt (ν k))⁻¹ * ((Real.sqrt (ν l))⁻¹ * (2 * Zb / |ρ k - ρ l| * η)) :=
          mul_le_mul_of_nonneg_left
            (mul_le_mul_of_nonneg_left hdotb (by positivity)) (by positivity)
      _ = (Real.sqrt (ν k))⁻¹ * (Real.sqrt (ν l))⁻¹ * (2 * Zb / |ρ k - ρ l|) * η := by ring
      _ ≤ Zb * gramC ρ ν * η := mul_le_mul_of_nonneg_right hgt hη0
      _ ≤ gramC ρ ν * (1 + Zb) * η := by
          calc Zb * gramC ρ ν * η = gramC ρ ν * η * Zb := by ring
            _ ≤ gramC ρ ν * η * (1 + Zb) := by
                apply mul_le_mul_of_nonneg_left _ (mul_nonneg (gramC_nonneg ρ ν) hη0)
                linarith
            _ = gramC ρ ν * (1 + Zb) * η := by ring

/-- **The deterministic core with frame columns `Q * Z`** (mirror of `align_detG`). The
split `S = W + Q Qᵀ` keeps its `r` columns; the frame is the `s` columns of `Q * Z` for a
general `Z : Matrix (Fin r) (Fin s) ℝ` with column sums `hZb`, rather than a hard selection
`Q.submatrix id f`. `hρinj` (instead of `hf : Function.Injective f`) drives the Gram bound
through `gram_bound_Z`, whose constant `gramC ρ ν * (1 + Zb) * η` replaces `gramC ρ ν * η`
throughout. The consumer `RankR/Het/Outliers.lean:922` calls `align_detG` in this general
`t`/`Lb` shape directly, so that is the shape mirrored here. -/
theorem align_detZ {W S : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin r) ℝ}
    {Z : Matrix (Fin r) (Fin s) ℝ}
    (hW : W.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W + Q * Qᵀ)
    {ρ ν : Fin s → ℝ} (hρinj : Function.Injective ρ) {τ mg : ℝ} (hmg : 0 < mg)
    (hτ : lamMax W hW ≤ τ) (hρ : ∀ k, τ + mg ≤ ρ k) (hν : ∀ k, 0 < ν k)
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {Cq : ℝ} (hCq : 0 ≤ Cq) (hcol : ∀ l : Fin r, ∑ i, Q i l ^ 2 ≤ Cq)
    {Zb : ℝ} (hZb0 : 0 ≤ Zb) (hZb : ∀ k, ∑ l, |Z l k| ≤ Zb)
    {v : Fin p → ℝ} (hv : v ⬝ᵥ v = 1)
    {t : Fin s → ℝ} {Lb η : ℝ} (hLb : ∀ k, |t k| ≤ Lb) (hη0 : 0 < η) (hη1 : η ≤ 1)
    (hE1 : ∀ (l : Fin r) (k : Fin s),
      |cform W (ρ k) (fun i => Q i l) (fun i => (Q * Z) i k) + Z l k| ≤ η)
    (hE2 : ∀ k : Fin s, |cform2 W (ρ k) (fun i => (Q * Z) i k) (fun i => (Q * Z) i k)
      - ν k| ≤ η)
    (hE3 : ∀ k : Fin s, |cform W (ρ k) v (fun i => (Q * Z) i k) - t k| ≤ η)
    (hsmall : (s : ℝ) * (gramC ρ ν * (1 + Zb) * η + resCG Cq r ν * η / mg) ≤ 1 / 2) :
    |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2 - ∑ k, t k ^ 2 / ν k|
      ≤ 5 * s * (gramC ρ ν * (1 + Zb) * η + resCG Cq r ν * η / mg)
        + (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by
  classical
  have hτz : ∀ k, lamMax W hW < ρ k := fun k =>
    lt_of_le_of_lt hτ (lt_of_lt_of_le (lt_add_of_pos_right τ hmg) (hρ k))
  have hresid : ∀ k, ‖toOp S (yhatv W (Q * Z) ρ ν k) - ρ k • yhatv W (Q * Z) ρ ν k‖
      ≤ resCG Cq r ν * η :=
    fun k => residual_bound_Z hW hSeq hν hτz hCq hη0.le hcol hE1 k
  have hgrame : ∀ k l, |⟪yhatv W (Q * Z) ρ ν k, yhatv W (Q * Z) ρ ν l⟫_ℝ
        - if k = l then 1 else 0| ≤ gramC ρ ν * (1 + Zb) * η :=
    fun k l => gram_bound_Z hW hρinj hν hτz hZb0 hη0.le hZb hE1 hE2 k l
  have hδ0 : (0 : ℝ) ≤ resCG Cq r ν * η := mul_nonneg (resCG_nonneg Cq r ν) hη0.le
  have hδ'0 : (0 : ℝ) ≤ gramC ρ ν * (1 + Zb) * η :=
    mul_nonneg (mul_nonneg (gramC_nonneg ρ ν) (by linarith)) hη0.le
  have hframe := Frame.specProj_frame_approx hS hmg hδ0 hδ'0 hI hρ hresid hgrame hsmall
    (WithLp.toLp 2 v)
  have hx1 : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot, WithLp.ofLp_toLp, hv]
  rw [hx1, mul_one] at hframe
  have hover : ∀ k, ⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ
      = (Real.sqrt (ν k))⁻¹ * cform W (ρ k) v (fun i => (Q * Z) i k) := by
    intro k
    rw [Frame.inner_eq_dot]
    have h1 : WithLp.ofLp (yhatv W (Q * Z) ρ ν k)
        = (Real.sqrt (ν k))⁻¹ • xtv W (Q * Z) ρ k := rfl
    rw [h1, WithLp.ofLp_toLp, smul_dotProduct, smul_eq_mul]
    congr 1
    rw [dotProduct_comm]
    rfl
  have hterm : ∀ k, |⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
      - t k ^ 2 / ν k| ≤ (1 + 2 * Lb) * (ν k)⁻¹ * η := by
    intro k
    have hb := hE3 k
    have h1 : ⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
        = (ν k)⁻¹ * cform W (ρ k) v (fun i => (Q * Z) i k) ^ 2 := by
      rw [hover k, mul_pow, inv_pow, Real.sq_sqrt (hν k).le]
    have h2 : |cform W (ρ k) v (fun i => (Q * Z) i k) ^ 2 - t k ^ 2| ≤ (1 + 2 * Lb) * η := by
      have habs2 : |cform W (ρ k) v (fun i => (Q * Z) i k) ^ 2 - t k ^ 2|
          = |cform W (ρ k) v (fun i => (Q * Z) i k) - t k|
            * |cform W (ρ k) v (fun i => (Q * Z) i k) + t k| := by
        rw [← abs_mul]; congr 1; ring
      have h3 : |cform W (ρ k) v (fun i => (Q * Z) i k) + t k| ≤ η + 2 * Lb := by
        have h4 : cform W (ρ k) v (fun i => (Q * Z) i k) + t k
            = (cform W (ρ k) v (fun i => (Q * Z) i k) - t k) + 2 * t k := by ring
        rw [h4]
        calc |(cform W (ρ k) v (fun i => (Q * Z) i k) - t k) + 2 * t k|
            ≤ |cform W (ρ k) v (fun i => (Q * Z) i k) - t k| + |2 * t k| := abs_add_le _ _
          _ ≤ η + 2 * Lb := by
              rw [abs_mul, abs_two]
              exact add_le_add hb (by linarith [hLb k])
      rw [habs2]
      calc |cform W (ρ k) v (fun i => (Q * Z) i k) - t k|
          * |cform W (ρ k) v (fun i => (Q * Z) i k) + t k|
          ≤ η * (η + 2 * Lb) := mul_le_mul hb h3 (abs_nonneg _) hη0.le
        _ ≤ (1 + 2 * Lb) * η := by nlinarith [hη0.le, hη1]
    have hgoal : |⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - t k ^ 2 / ν k|
        = (ν k)⁻¹ * |cform W (ρ k) v (fun i => (Q * Z) i k) ^ 2 - t k ^ 2| := by
      rw [h1, show t k ^ 2 / ν k = (ν k)⁻¹ * t k ^ 2 by rw [div_eq_mul_inv, mul_comm],
        ← mul_sub, abs_mul, abs_of_nonneg (inv_nonneg.mpr (hν k).le)]
    rw [hgoal]
    calc (ν k)⁻¹ * |cform W (ρ k) v (fun i => (Q * Z) i k) ^ 2 - t k ^ 2|
        ≤ (ν k)⁻¹ * ((1 + 2 * Lb) * η) :=
          mul_le_mul_of_nonneg_left h2 (inv_nonneg.mpr (hν k).le)
      _ = (1 + 2 * Lb) * (ν k)⁻¹ * η := by ring
  have hsum : |∑ k, ⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - ∑ k, t k ^ 2 / ν k|
      ≤ (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by
    rw [← Finset.sum_sub_distrib]
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ k, |⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2 - t k ^ 2 / ν k|
        ≤ ∑ k, (1 + 2 * Lb) * (ν k)⁻¹ * η := Finset.sum_le_sum fun k _ => hterm k
      _ = (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := by rw [← Finset.sum_mul, ← Finset.mul_sum]
  calc |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2 - ∑ k, t k ^ 2 / ν k|
      ≤ |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 v)‖ ^ 2
          - ∑ k, ⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2|
        + |∑ k, ⟪yhatv W (Q * Z) ρ ν k, WithLp.toLp 2 v⟫_ℝ ^ 2
          - ∑ k, t k ^ 2 / ν k| := abs_sub_le _ _ _
    _ ≤ 5 * s * (gramC ρ ν * (1 + Zb) * η + resCG Cq r ν * η / mg)
        + (1 + 2 * Lb) * (∑ k, (ν k)⁻¹) * η := add_le_add hframe hsum

end Core

end OutliersR

end StackedSVD
