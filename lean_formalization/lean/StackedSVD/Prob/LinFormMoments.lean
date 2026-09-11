/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# Moments of a linear form of i.i.d. coordinates

Let `ν` be a probability law on `ℝ` with mean 0, variance 1 and an integrable fourth power, and
let `y` carry the product law `Measure.pi fun _ : Fin D => ν`. This file gives the first two
moments and a bound on the fourth moment of the linear form `linForm a y = ∑ l, a l * y l`.

The proof is one induction on `D`. The step splits off the first coordinate with the measure
preserving equivalence of `MeasureTheory.measurePreserving_piFinSuccAbove`, which carries the
product law on `Fin (D + 1) → ℝ` to `ν.prod (Measure.pi fun _ : Fin D => ν)`. It then expands
`(a 0 * x + u) ^ k` and integrates each term `x ^ j * u ^ m` with `integral_prod_mul`.

The exact fourth moment is `3 * (∑ l, a l ^ 2) ^ 2 + (ν₄ - 3) * ∑ l, a l ^ 4`, with
`ν₄ = ∫ x ^ 4 ∂ν`. The form used downstream is the bound
`∫ (linForm a) ^ 4 ≤ ν₄ * ∑ l, a l ^ 4 + 3 * (∑ l, a l ^ 2) ^ 2`, which is what the induction
produces directly. It is the one ingredient of the non-Gaussian discharge of `thm:theta_est`
that the Gaussian proof took from rotation invariance.
-/

open MeasureTheory ProbabilityTheory

namespace StackedSVD
namespace LinForm

variable {ν : Measure ℝ} [IsProbabilityMeasure ν]

/-- The linear form `∑ l, a l * y l`. -/
def linForm {D : ℕ} (a : Fin D → ℝ) (y : Fin D → ℝ) : ℝ := ∑ l, a l * y l

/-- The linear form is measurable in its argument. -/
theorem measurable_linForm {D : ℕ} (a : Fin D → ℝ) : Measurable (linForm a) := by
  unfold linForm
  fun_prop

/-! ### Integrability of the low powers -/

/-- If `f ^ 4` is integrable on a finite measure, so is `f ^ k` for every `k ≤ 4`. The bound is
`|t| ^ k ≤ 1 + t ^ 4`, split at `|t| = 1`. -/
private theorem integrable_pow_of_pow_four {α : Type*} [MeasurableSpace α] {μ : Measure α}
    [IsFiniteMeasure μ] {f : α → ℝ} (hf : Measurable f)
    (h4 : Integrable (fun x => f x ^ 4) μ) {k : ℕ} (hk : k ≤ 4) :
    Integrable (fun x => f x ^ k) μ := by
  refine Integrable.mono' ((integrable_const (1 : ℝ)).add h4)
    ((hf.pow_const k).aestronglyMeasurable) ?_
  filter_upwards with x
  simp only [Pi.add_apply, Real.norm_eq_abs, abs_pow]
  have h0 : (0 : ℝ) ≤ |f x| := abs_nonneg _
  have h4' : |f x| ^ 4 = f x ^ 4 := by
    rw [← abs_pow, abs_of_nonneg (by positivity)]
  have hnn : (0 : ℝ) ≤ |f x| ^ 4 := pow_nonneg h0 4
  rcases le_total |f x| 1 with h | h
  · have hle : |f x| ^ k ≤ 1 := pow_le_one₀ h0 h
    linarith
  · have hle : |f x| ^ k ≤ |f x| ^ 4 := pow_le_pow_right₀ h hk
    linarith

/-- Every power `x ^ k` with `k ≤ 4` is integrable when `x ^ 4` is: `|x| ^ k ≤ 1 + x ^ 4`. -/
theorem integrable_pow_of_four (h4 : Integrable (fun x : ℝ => x ^ 4) ν) {k : ℕ} (hk : k ≤ 4) :
    Integrable (fun x : ℝ => x ^ k) ν :=
  integrable_pow_of_pow_four measurable_id h4 hk

/-- The same for the linear form: `(linForm a) ^ k` integrable for `k ≤ 4` once `k = 4` is. -/
theorem integrable_linForm_pow_of_four {D : ℕ} {a : Fin D → ℝ}
    (h4 : Integrable (fun y => linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν)) {k : ℕ}
    (hk : k ≤ 4) :
    Integrable (fun y => linForm a y ^ k) (Measure.pi fun _ : Fin D => ν) :=
  integrable_pow_of_pow_four (measurable_linForm a) h4 hk

/-! ### The split of the first coordinate -/

/-- Splitting off the first coordinate in an integral of a function of the linear form. -/
private theorem integral_linForm_succ {D : ℕ} (a : Fin (D + 1) → ℝ) (F : ℝ → ℝ) :
    ∫ y, F (linForm a y) ∂(Measure.pi fun _ : Fin (D + 1) => ν)
      = ∫ z : ℝ × (Fin D → ℝ), F (a 0 * z.1 + linForm (fun j => a j.succ) z.2)
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) := by
  rw [← ((measurePreserving_piFinSuccAbove (fun _ : Fin (D + 1) => ν) 0).symm).integral_comp']
  simp [linForm, Fin.sum_univ_succ]

/-- Splitting off the first coordinate in an integrability statement. -/
private theorem integrable_linForm_succ {D : ℕ} (a : Fin (D + 1) → ℝ) (F : ℝ → ℝ) :
    Integrable (fun y => F (linForm a y)) (Measure.pi fun _ : Fin (D + 1) => ν)
      ↔ Integrable (fun z : ℝ × (Fin D → ℝ) => F (a 0 * z.1 + linForm (fun j => a j.succ) z.2))
          (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
  rw [← ((measurePreserving_piFinSuccAbove (fun _ : Fin (D + 1) => ν) 0).symm).integrable_comp_emb
    (MeasurableEquiv.measurableEmbedding _)]
  simp [linForm, Fin.sum_univ_succ, Function.comp_def]

/-- The power form of `integral_linForm_succ`. -/
private theorem integral_linForm_pow_succ {D : ℕ} (a : Fin (D + 1) → ℝ) (k : ℕ) :
    ∫ y, linForm a y ^ k ∂(Measure.pi fun _ : Fin (D + 1) => ν)
      = ∫ z : ℝ × (Fin D → ℝ), (a 0 * z.1 + linForm (fun j => a j.succ) z.2) ^ k
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) :=
  integral_linForm_succ a (fun t => t ^ k)

/-- The power form of `integrable_linForm_succ`. -/
private theorem integrable_linForm_pow_succ {D : ℕ} (a : Fin (D + 1) → ℝ) (k : ℕ) :
    Integrable (fun y => linForm a y ^ k) (Measure.pi fun _ : Fin (D + 1) => ν)
      ↔ Integrable (fun z : ℝ × (Fin D → ℝ) =>
          (a 0 * z.1 + linForm (fun j => a j.succ) z.2) ^ k)
          (ν.prod (Measure.pi fun _ : Fin D => ν)) :=
  integrable_linForm_succ a (fun t => t ^ k)

/-! ### The product integrals -/

/-- One term of the expansion: a product of a power of the first coordinate and a power of the
linear form of the rest. -/
private theorem integral_term {D : ℕ} (c : ℝ) (k m : ℕ) (g : (Fin D → ℝ) → ℝ) :
    ∫ z : ℝ × (Fin D → ℝ), c * (z.1 ^ k * g z.2 ^ m)
        ∂(ν.prod (Measure.pi fun _ : Fin D => ν))
      = c * ((∫ x, x ^ k ∂ν) * ∫ y, g y ^ m ∂(Measure.pi fun _ : Fin D => ν)) := by
  rw [integral_const_mul, integral_prod_mul (fun x : ℝ => x ^ k) (fun y => g y ^ m)]

/-- Addition of five integrals. -/
private theorem integral_add₅ {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f₁ f₂ f₃ f₄ f₅ : α → ℝ} (h₁ : Integrable f₁ μ) (h₂ : Integrable f₂ μ)
    (h₃ : Integrable f₃ μ) (h₄ : Integrable f₄ μ) (h₅ : Integrable f₅ μ) :
    ∫ z, (f₁ z + f₂ z + f₃ z + f₄ z + f₅ z) ∂μ
      = (∫ z, f₁ z ∂μ) + (∫ z, f₂ z ∂μ) + (∫ z, f₃ z ∂μ) + (∫ z, f₄ z ∂μ)
        + ∫ z, f₅ z ∂μ := by
  have e₁ : Integrable (fun z => f₁ z + f₂ z) μ := h₁.add h₂
  have e₂ : Integrable (fun z => f₁ z + f₂ z + f₃ z) μ := e₁.add h₃
  have e₃ : Integrable (fun z => f₁ z + f₂ z + f₃ z + f₄ z) μ := e₂.add h₄
  rw [integral_add e₃ h₅, integral_add e₂ h₄, integral_add e₁ h₃, integral_add h₁ h₂]

/-- Addition of three integrals. -/
private theorem integral_add₃ {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f₁ f₂ f₃ : α → ℝ} (h₁ : Integrable f₁ μ) (h₂ : Integrable f₂ μ) (h₃ : Integrable f₃ μ) :
    ∫ z, (f₁ z + f₂ z + f₃ z) ∂μ
      = (∫ z, f₁ z ∂μ) + (∫ z, f₂ z ∂μ) + ∫ z, f₃ z ∂μ := by
  have e₁ : Integrable (fun z => f₁ z + f₂ z) μ := h₁.add h₂
  rw [integral_add e₁ h₃, integral_add h₁ h₂]

/-! ### The moments -/

/-- The first, second and fourth moments of the linear form under the product law, by induction
on the number `D` of coordinates. -/
theorem moments (hmean : ∫ x, x ∂ν = 0) (hvar : ∫ x, x ^ 2 ∂ν = 1)
    (h4 : Integrable (fun x : ℝ => x ^ 4) ν) :
    ∀ (D : ℕ) (a : Fin D → ℝ),
      Integrable (fun y => linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν) ∧
      ∫ y, linForm a y ∂(Measure.pi fun _ : Fin D => ν) = 0 ∧
      ∫ y, linForm a y ^ 2 ∂(Measure.pi fun _ : Fin D => ν) = ∑ l, a l ^ 2 ∧
      ∫ y, linForm a y ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
        ≤ (∫ x, x ^ 4 ∂ν) * ∑ l, a l ^ 4 + 3 * (∑ l, a l ^ 2) ^ 2 := by
  intro D
  induction D with
  | zero => intro a; refine ⟨?_, ?_, ?_, ?_⟩ <;> simp [linForm]
  | succ D ih =>
    intro a
    obtain ⟨hb4, hb1, hb2, hb4le⟩ := ih (fun j => a j.succ)
    -- moments of one coordinate
    have hx : ∀ k, k ≤ 4 → Integrable (fun x : ℝ => x ^ k) ν := fun k hk =>
      integrable_pow_of_four h4 hk
    have hg : ∀ k, k ≤ 4 → Integrable (fun y => linForm (fun j => a j.succ) y ^ k)
        (Measure.pi fun _ : Fin D => ν) := fun k hk => integrable_linForm_pow_of_four hb4 hk
    have hterm : ∀ (c : ℝ) (k m : ℕ), k ≤ 4 → m ≤ 4 →
        Integrable (fun z : ℝ × (Fin D → ℝ) =>
          c * (z.1 ^ k * linForm (fun j => a j.succ) z.2 ^ m))
          (ν.prod (Measure.pi fun _ : Fin D => ν)) :=
      fun c k m hk hm => ((hx k hk).mul_prod (hg m hm)).const_mul c
    have hx0 : ∫ x : ℝ, x ^ 0 ∂ν = 1 := by simp
    have hx1 : ∫ x : ℝ, x ^ 1 ∂ν = 0 := by simpa using hmean
    have hg0 : ∫ y, linForm (fun j => a j.succ) y ^ 0 ∂(Measure.pi fun _ : Fin D => ν) = 1 := by
      simp
    have hg1 : ∫ y, linForm (fun j => a j.succ) y ^ 1 ∂(Measure.pi fun _ : Fin D => ν) = 0 := by
      simpa using hb1
    -- the fourth power is integrable
    have hI4 : Integrable (fun y => linForm a y ^ 4)
        (Measure.pi fun _ : Fin (D + 1) => ν) := by
      rw [integrable_linForm_pow_succ a 4]
      have hfun : (fun z : ℝ × (Fin D → ℝ) =>
          (a 0 * z.1 + linForm (fun j => a j.succ) z.2) ^ 4)
          = fun z : ℝ × (Fin D → ℝ) =>
            a 0 ^ 4 * (z.1 ^ 4 * linForm (fun j => a j.succ) z.2 ^ 0)
            + 4 * a 0 ^ 3 * (z.1 ^ 3 * linForm (fun j => a j.succ) z.2 ^ 1)
            + 6 * a 0 ^ 2 * (z.1 ^ 2 * linForm (fun j => a j.succ) z.2 ^ 2)
            + 4 * a 0 * (z.1 ^ 1 * linForm (fun j => a j.succ) z.2 ^ 3)
            + 1 * (z.1 ^ 0 * linForm (fun j => a j.succ) z.2 ^ 4) := by
        funext z; ring
      rw [hfun]
      exact ((((hterm _ 4 0 (by norm_num) (by norm_num)).add
        (hterm _ 3 1 (by norm_num) (by norm_num))).add
        (hterm _ 2 2 (by norm_num) (by norm_num))).add
        (hterm _ 1 3 (by norm_num) (by norm_num))).add
        (hterm _ 0 4 (by norm_num) (by norm_num))
    refine ⟨hI4, ?_, ?_, ?_⟩
    · -- the mean
      have h1 : ∫ y, linForm a y ∂(Measure.pi fun _ : Fin (D + 1) => ν)
          = ∫ z : ℝ × (Fin D → ℝ),
              (a 0 * (z.1 ^ 1 * linForm (fun j => a j.succ) z.2 ^ 0)
                + 1 * (z.1 ^ 0 * linForm (fun j => a j.succ) z.2 ^ 1))
              ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) := by
        rw [integral_linForm_succ a (fun t => t)]
        apply integral_congr_ae
        filter_upwards with z
        ring
      rw [h1, integral_add (hterm _ 1 0 (by norm_num) (by norm_num))
        (hterm _ 0 1 (by norm_num) (by norm_num)), integral_term, integral_term,
        hx0, hx1, hg0, hg1]
      ring
    · -- the second moment
      have h2 : ∫ y, linForm a y ^ 2 ∂(Measure.pi fun _ : Fin (D + 1) => ν)
          = ∫ z : ℝ × (Fin D → ℝ),
              (a 0 ^ 2 * (z.1 ^ 2 * linForm (fun j => a j.succ) z.2 ^ 0)
                + 2 * a 0 * (z.1 ^ 1 * linForm (fun j => a j.succ) z.2 ^ 1)
                + 1 * (z.1 ^ 0 * linForm (fun j => a j.succ) z.2 ^ 2))
              ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) := by
        rw [integral_linForm_pow_succ a 2]
        apply integral_congr_ae
        filter_upwards with z
        ring
      rw [h2, integral_add₃ (hterm _ 2 0 (by norm_num) (by norm_num))
        (hterm _ 1 1 (by norm_num) (by norm_num)) (hterm _ 0 2 (by norm_num) (by norm_num)),
        integral_term, integral_term, integral_term, hx0, hx1, hvar, hg0, hg1, hb2,
        Fin.sum_univ_succ]
      ring
    · -- the fourth moment
      have h4' : ∫ y, linForm a y ^ 4 ∂(Measure.pi fun _ : Fin (D + 1) => ν)
          = ∫ z : ℝ × (Fin D → ℝ),
              (a 0 ^ 4 * (z.1 ^ 4 * linForm (fun j => a j.succ) z.2 ^ 0)
                + 4 * a 0 ^ 3 * (z.1 ^ 3 * linForm (fun j => a j.succ) z.2 ^ 1)
                + 6 * a 0 ^ 2 * (z.1 ^ 2 * linForm (fun j => a j.succ) z.2 ^ 2)
                + 4 * a 0 * (z.1 ^ 1 * linForm (fun j => a j.succ) z.2 ^ 3)
                + 1 * (z.1 ^ 0 * linForm (fun j => a j.succ) z.2 ^ 4))
              ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) := by
        rw [integral_linForm_pow_succ a 4]
        apply integral_congr_ae
        filter_upwards with z
        ring
      rw [h4', integral_add₅ (hterm _ 4 0 (by norm_num) (by norm_num))
        (hterm _ 3 1 (by norm_num) (by norm_num)) (hterm _ 2 2 (by norm_num) (by norm_num))
        (hterm _ 1 3 (by norm_num) (by norm_num)) (hterm _ 0 4 (by norm_num) (by norm_num)),
        integral_term, integral_term, integral_term, integral_term, integral_term,
        hx0, hx1, hvar, hg0, hg1, hb2, Fin.sum_univ_succ, Fin.sum_univ_succ]
      have hA : (0 : ℝ) ≤ a 0 ^ 4 := by positivity
      nlinarith [hb4le, hA]

/-- The fourth power of the linear form is integrable under the product law. -/
theorem integrable_linForm_pow_four (hmean : ∫ x, x ∂ν = 0) (hvar : ∫ x, x ^ 2 ∂ν = 1)
    (h4 : Integrable (fun x : ℝ => x ^ 4) ν) {D : ℕ} (a : Fin D → ℝ) :
    Integrable (fun y => linForm a y ^ 4) (Measure.pi fun _ : Fin D => ν) :=
  (moments hmean hvar h4 D a).1

/-- The linear form has mean zero. -/
theorem integral_linForm (hmean : ∫ x, x ∂ν = 0) (hvar : ∫ x, x ^ 2 ∂ν = 1)
    (h4 : Integrable (fun x : ℝ => x ^ 4) ν) {D : ℕ} (a : Fin D → ℝ) :
    ∫ y, linForm a y ∂(Measure.pi fun _ : Fin D => ν) = 0 :=
  (moments hmean hvar h4 D a).2.1

/-- The second moment of the linear form is `∑ l, a l ^ 2`. -/
theorem integral_linForm_sq (hmean : ∫ x, x ∂ν = 0) (hvar : ∫ x, x ^ 2 ∂ν = 1)
    (h4 : Integrable (fun x : ℝ => x ^ 4) ν) {D : ℕ} (a : Fin D → ℝ) :
    ∫ y, linForm a y ^ 2 ∂(Measure.pi fun _ : Fin D => ν) = ∑ l, a l ^ 2 :=
  (moments hmean hvar h4 D a).2.2.1

/-- The fourth moment of the linear form is at most
`(∫ x ^ 4 ∂ν) * ∑ l, a l ^ 4 + 3 * (∑ l, a l ^ 2) ^ 2`. -/
theorem integral_linForm_pow_four_le (hmean : ∫ x, x ∂ν = 0) (hvar : ∫ x, x ^ 2 ∂ν = 1)
    (h4 : Integrable (fun x : ℝ => x ^ 4) ν) {D : ℕ} (a : Fin D → ℝ) :
    ∫ y, linForm a y ^ 4 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ (∫ x, x ^ 4 ∂ν) * ∑ l, a l ^ 4 + 3 * (∑ l, a l ^ 2) ^ 2 :=
  (moments hmean hvar h4 D a).2.2.2

end LinForm
end StackedSVD
