/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.NoiseLaw
import StackedSVD.Prob.LinFormMoments

/-! # Matrix moments at a general noise law

The two second-moment facts that the Gaussian proof of `thm:theta_est` takes from
`ThetaEst.integral_sq_dotProduct_mulVec` and `ThetaEst.integral_sq_normSqMulVec`, now for a
general noise law `ν` (`NoiseLaw ν`, `Prob/NoiseLaw.lean`): `∫ (x ⬝ᵥ Z a) ^ 2 = 1` for unit
`x` and `a`, and `∫ (D⁻¹ ‖Z a‖² - r / D)² ≤ (ν₄ + 2) r / D²` for a unit `a`, with
`ν₄ = ∫ x ^ 4 ∂ν`.

The Gaussian proof pushes `Z ↦ Z a` forward to a standard Gaussian vector and reads off the
answer by rotation invariance; that route is not available at a general `ν`. Here row `k` of
`Z` enters only through the linear form `LinForm.linForm a (Z k)`, whose first four moments
`Prob/LinFormMoments.lean` gives by an induction on the number of coordinates. The sum over
rows is then handled by `PiLaw.Centered` (`Prob/NoiseLaw.lean`), a copy of `R2.Centered` at a
general product law. The bound `(ν₄ + 2) r / D²` replaces the Gaussian exact value
`varSq · r / D²` of `ThetaEst.lean`; it comes from the fourth-moment bound on `linForm`, so it
is not claimed to be tight. -/

open MeasureTheory ProbabilityTheory
open scoped Matrix

namespace StackedSVD
namespace NoiseLaw

variable {ν : Measure ℝ}

/-- The `k`-th coordinate of `Z *ᵥ a` is the linear form of row `k`. -/
theorem mulVec_apply_eq_linForm {r D : ℕ} (Z : Matrix (Fin r) (Fin D) ℝ) (a : Fin D → ℝ)
    (k : Fin r) : (Z *ᵥ a) k = LinForm.linForm a (Z k) := by
  simp [Matrix.mulVec, dotProduct, LinForm.linForm, mul_comm]

/-! ### Helper lemmas -/

/-- `y ⬝ᵥ y` as a sum of squares. -/
private theorem dotProduct_self_eq_sum_sq {n : ℕ} (y : Fin n → ℝ) : y ⬝ᵥ y = ∑ k, y k ^ 2 :=
  Finset.sum_congr rfl fun k _ => (sq (y k)).symm

/-- The sum of squares of nonnegative terms is at most the square of the sum. -/
private theorem sum_sq_le_sq_sum {n : ℕ} (b : Fin n → ℝ) (hb : ∀ l, 0 ≤ b l) :
    ∑ l, b l ^ 2 ≤ (∑ l, b l) ^ 2 := by
  have hexpand : (∑ l, b l) ^ 2 = ∑ l, ∑ m, b l * b m := by
    rw [pow_two, Finset.sum_mul_sum]
  rw [hexpand]
  refine Finset.sum_le_sum fun l _ => ?_
  rw [pow_two]
  exact Finset.single_le_sum (fun m _ => mul_nonneg (hb l) (hb m)) (Finset.mem_univ l)

/-- The pointwise identity behind the second moment of `D⁻¹ ‖·‖² - r / D`. The model is
`ThetaEst.center_eq` (`ThetaEst.lean:397`), the Gaussian-file version at `y = Z a` directly;
here it is stated once for a general vector `y` and reused at `y k := linForm a (Z k)`. -/
private theorem center_eq (D r : ℕ) (y : Fin r → ℝ) :
    ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2
      = ((D : ℝ)⁻¹) ^ 2 * (∑ k, (1 : ℝ) * (y k ^ 2 - 1)) ^ 2 := by
  have hL : (∑ k, (1 : ℝ) * (y k ^ 2 - 1)) = (∑ k, y k ^ 2) - (r : ℝ) := by
    simp [Finset.sum_sub_distrib]
  rw [dotProduct_self_eq_sum_sq, hL, div_eq_inv_mul, ← mul_sub, mul_pow]

/-- The cross term `x ⬝ᵥ (Z *ᵥ a)`, squared, as a weighted sum of the linear form over rows. -/
private theorem dotProduct_mulVec_sq_eq {r D : ℕ} (x : Fin r → ℝ) (a : Fin D → ℝ)
    (Z : Matrix (Fin r) (Fin D) ℝ) :
    (x ⬝ᵥ (Z *ᵥ a)) ^ 2 = (∑ k, x k * LinForm.linForm a (Z k)) ^ 2 := by
  have hZ : x ⬝ᵥ (Z *ᵥ a) = ∑ k, x k * LinForm.linForm a (Z k) :=
    Finset.sum_congr rfl fun k _ => by rw [mulVec_apply_eq_linForm]
  rw [hZ]

/-- The centered squared norm `D⁻¹ ‖Z a‖² - r / D`, squared, as `center_eq` applied to the
rows of `Z` through the linear form. -/
private theorem sq_normSqMulVec_eq {r D : ℕ} (a : Fin D → ℝ) (Z : Matrix (Fin r) (Fin D) ℝ) :
    ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2
      = ((D : ℝ)⁻¹) ^ 2 * (∑ k, (1 : ℝ) * (LinForm.linForm a (Z k) ^ 2 - 1)) ^ 2 := by
  have hyy : (Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)
      = (fun k => LinForm.linForm a (Z k)) ⬝ᵥ (fun k => LinForm.linForm a (Z k)) :=
    Finset.sum_congr rfl fun k _ => by rw [mulVec_apply_eq_linForm]
  rw [hyy]
  exact center_eq D r (fun k => LinForm.linForm a (Z k))

/-! ### The two `Centered` witnesses -/

/-- The linear form of a row is centered with second moment `∑ a_l²`. -/
theorem centered_linForm (hν : NoiseLaw ν) {D : ℕ} (a : Fin D → ℝ) :
    PiLaw.Centered (Measure.pi fun _ : Fin D => ν) (LinForm.linForm a) := by
  have := hν.prob
  have h4 : Integrable (fun y => LinForm.linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν) :=
    LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 a
  refine ⟨LinForm.measurable_linForm a, ?_, ?_, ?_⟩
  · simpa using LinForm.integrable_linForm_pow_of_four h4 (by norm_num : (1 : ℕ) ≤ 4)
  · exact LinForm.integrable_linForm_pow_of_four h4 (by norm_num : (2 : ℕ) ≤ 4)
  · exact LinForm.integral_linForm hν.mean hν.var hν.mom4 a

/-- `linForm a ^ 2 - 1` is centered for a unit `a`. -/
theorem centered_linForm_sq_sub_one (hν : NoiseLaw ν) {D : ℕ} {a : Fin D → ℝ}
    (ha : a ⬝ᵥ a = 1) :
    PiLaw.Centered (Measure.pi fun _ : Fin D => ν) (fun y => LinForm.linForm a y ^ 2 - 1) := by
  have := hν.prob
  have hc := centered_linForm hν a
  have h4 : Integrable (fun y => LinForm.linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν) :=
    LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 a
  have hax : ∑ l, a l ^ 2 = 1 := (dotProduct_self_eq_sum_sq a).symm.trans ha
  refine ⟨(LinForm.measurable_linForm a).pow_const 2 |>.sub measurable_const,
    hc.sqInt.sub' (integrable_const 1), ?_, ?_⟩
  · have hexp2 : Integrable (fun y => LinForm.linForm a y ^ 4 - 2 * LinForm.linForm a y ^ 2 + 1)
        (Measure.pi fun _ : Fin D => ν) :=
      (h4.sub' (hc.sqInt.const_mul 2)).add (integrable_const 1)
    have heq : (fun y : Fin D → ℝ => (LinForm.linForm a y ^ 2 - 1) ^ 2)
        = (fun y => LinForm.linForm a y ^ 4 - 2 * LinForm.linForm a y ^ 2 + 1) := by
      funext y; ring
    rw [heq]; exact hexp2
  · rw [integral_sub hc.sqInt (integrable_const 1),
      LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 a, hax]
    simp

/-- `∫ (g² - 1)² = ∫ g⁴ - 1 ≤ ν₄ ∑ a_l⁴ + 3 (∑ a_l²)² - 1 ≤ ν₄ + 2` for a unit `a`. -/
theorem integral_sq_linForm_sq_sub_one_le (hν : NoiseLaw ν) {D : ℕ} {a : Fin D → ℝ}
    (ha : a ⬝ᵥ a = 1) :
    ∫ y, (LinForm.linForm a y ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ (∫ x, x ^ 4 ∂ν) + 2 := by
  have := hν.prob
  have hax : ∑ l, a l ^ 2 = 1 := (dotProduct_self_eq_sum_sq a).symm.trans ha
  have ha4 : ∑ l, a l ^ 4 ≤ 1 := by
    have hpt : ∀ l, a l ^ 4 = (a l ^ 2) ^ 2 := fun l => by ring
    calc ∑ l, a l ^ 4 = ∑ l, (a l ^ 2) ^ 2 := Finset.sum_congr rfl fun l _ => hpt l
      _ ≤ (∑ l, a l ^ 2) ^ 2 := sum_sq_le_sq_sum (fun l => a l ^ 2) (fun l => sq_nonneg (a l))
      _ = 1 := by rw [hax]; norm_num
  have h4 : Integrable (fun y => LinForm.linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν) :=
    LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 a
  have hsq : Integrable (fun y => LinForm.linForm a y ^ 2) (Measure.pi fun _ : Fin D => ν) :=
    LinForm.integrable_linForm_pow_of_four h4 (by norm_num : (2 : ℕ) ≤ 4)
  have heq : ∫ y, (LinForm.linForm a y ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      = (∫ y, LinForm.linForm a y ^ 4 ∂(Measure.pi fun _ : Fin D => ν)) - 1 := by
    have hpt : (fun y : Fin D → ℝ => (LinForm.linForm a y ^ 2 - 1) ^ 2)
        = (fun y => LinForm.linForm a y ^ 4 - 2 * LinForm.linForm a y ^ 2 + 1) := by
      funext y; ring
    rw [hpt, integral_add (h4.sub' (hsq.const_mul 2)) (integrable_const 1),
      integral_sub h4 (hsq.const_mul 2), integral_const_mul,
      LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 a, hax]
    simp
    ring
  have hle := LinForm.integral_linForm_pow_four_le hν.mean hν.var hν.mom4 a
  rw [hax] at hle
  have hnn4 : 0 ≤ ∫ x, x ^ 4 ∂ν := integral_nonneg fun x => by positivity
  have hCT : (∫ x, x ^ 4 ∂ν) * (∑ l, a l ^ 4) ≤ (∫ x, x ^ 4 ∂ν) * 1 :=
    mul_le_mul_of_nonneg_left ha4 hnn4
  rw [heq]
  nlinarith [hle, hCT]

/-! ### The two moment theorems -/

/-- **The cross term, second moment.** `∫ (x ⬝ᵥ Z a)² = 1` for unit `x` and unit `a`. -/
theorem integral_sq_dotProduct_mulVec (hν : NoiseLaw ν) {r D : ℕ} {x : Fin r → ℝ}
    (hx : x ⬝ᵥ x = 1) {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    ∫ Z, (x ⬝ᵥ (Z *ᵥ a)) ^ 2 ∂(noiseMatrix ν r D) = 1 := by
  have := hν.prob
  have hC := centered_linForm hν a
  have key : ∫ Z : Matrix (Fin r) (Fin D) ℝ, (∑ k, x k * LinForm.linForm a (Z k)) ^ 2
        ∂(noiseMatrix ν r D)
      = (∫ t, LinForm.linForm a t ^ 2 ∂(Measure.pi fun _ : Fin D => ν)) * ∑ k, x k ^ 2 :=
    PiLaw.Centered.integral_sq_sum hC x
  have hax : ∑ l, a l ^ 2 = 1 := (dotProduct_self_eq_sum_sq a).symm.trans ha
  have hxx : ∑ k, x k ^ 2 = 1 := (dotProduct_self_eq_sum_sq x).symm.trans hx
  simp only [dotProduct_mulVec_sq_eq]
  rw [key, LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 a, hax, hxx]
  norm_num

theorem integrable_sq_dotProduct_mulVec (hν : NoiseLaw ν) {r D : ℕ} (x : Fin r → ℝ)
    (a : Fin D → ℝ) :
    Integrable (fun Z : Matrix (Fin r) (Fin D) ℝ => (x ⬝ᵥ (Z *ᵥ a)) ^ 2) (noiseMatrix ν r D) := by
  have := hν.prob
  have hC := centered_linForm hν a
  have key : Integrable
      (fun Z : Matrix (Fin r) (Fin D) ℝ => (∑ k, x k * LinForm.linForm a (Z k)) ^ 2)
      (noiseMatrix ν r D) := PiLaw.Centered.integrable_sq_sum hC x
  simpa only [dotProduct_mulVec_sq_eq] using key

/-- **Item P, the second moment.** `∫ (D⁻¹ ‖Z a‖² - r/D)² ≤ (ν₄ + 2) r / D²` for a unit `a`.
No hypothesis `0 < D`: at `D = 0` both sides are `0`. -/
theorem integral_sq_normSqMulVec_le (hν : NoiseLaw ν) {r D : ℕ} {a : Fin D → ℝ}
    (ha : a ⬝ᵥ a = 1) :
    ∫ Z, ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2 ∂(noiseMatrix ν r D)
      ≤ ((∫ x, x ^ 4 ∂ν) + 2) * r / (D : ℝ) ^ 2 := by
  have := hν.prob
  have hcs := centered_linForm_sq_sub_one hν ha
  have key : ∫ Z : Matrix (Fin r) (Fin D) ℝ,
        (∑ k, (1 : ℝ) * (LinForm.linForm a (Z k) ^ 2 - 1)) ^ 2 ∂(noiseMatrix ν r D)
      = (∫ t, (LinForm.linForm a t ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν))
        * ∑ _k : Fin r, (1 : ℝ) ^ 2 :=
    PiLaw.Centered.integral_sq_sum hcs (fun _ => 1)
  have hone : ∑ _k : Fin r, (1 : ℝ) ^ 2 = (r : ℝ) := by simp
  have hbound := integral_sq_linForm_sq_sub_one_le hν ha
  have hrnn : (0 : ℝ) ≤ (r : ℝ) := Nat.cast_nonneg r
  have hXr : (∫ t, (LinForm.linForm a t ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)) * (r : ℝ)
      ≤ ((∫ x, x ^ 4 ∂ν) + 2) * (r : ℝ) := mul_le_mul_of_nonneg_right hbound hrnn
  have hfinal : ((D : ℝ)⁻¹) ^ 2 *
        ((∫ t, (LinForm.linForm a t ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)) * (r : ℝ))
      ≤ ((D : ℝ)⁻¹) ^ 2 * (((∫ x, x ^ 4 ∂ν) + 2) * (r : ℝ)) :=
    mul_le_mul_of_nonneg_left hXr (sq_nonneg _)
  simp only [sq_normSqMulVec_eq]
  rw [integral_const_mul, key, hone]
  calc ((D : ℝ)⁻¹) ^ 2 *
        ((∫ t, (LinForm.linForm a t ^ 2 - 1) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)) * (r : ℝ))
      ≤ ((D : ℝ)⁻¹) ^ 2 * (((∫ x, x ^ 4 ∂ν) + 2) * (r : ℝ)) := hfinal
    _ = ((∫ x, x ^ 4 ∂ν) + 2) * r / (D : ℝ) ^ 2 := by rw [inv_pow]; ring

theorem integrable_sq_normSqMulVec (hν : NoiseLaw ν) {r D : ℕ} {a : Fin D → ℝ}
    (ha : a ⬝ᵥ a = 1) :
    Integrable (fun Z : Matrix (Fin r) (Fin D) ℝ =>
      ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2) (noiseMatrix ν r D) := by
  have := hν.prob
  have hcs := centered_linForm_sq_sub_one hν ha
  have key : Integrable
      (fun Z : Matrix (Fin r) (Fin D) ℝ =>
        (∑ k, (1 : ℝ) * (LinForm.linForm a (Z k) ^ 2 - 1)) ^ 2) (noiseMatrix ν r D) :=
    PiLaw.Centered.integrable_sq_sum hcs (fun _ => 1)
  have key2 : Integrable
      (fun Z : Matrix (Fin r) (Fin D) ℝ =>
        ((D : ℝ)⁻¹) ^ 2 * (∑ k, (1 : ℝ) * (LinForm.linForm a (Z k) ^ 2 - 1)) ^ 2)
      (noiseMatrix ν r D) := key.const_mul _
  simpa only [sq_normSqMulVec_eq] using key2

end NoiseLaw
end StackedSVD
