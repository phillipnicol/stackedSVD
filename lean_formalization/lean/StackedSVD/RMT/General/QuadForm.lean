/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.NoiseLaw
import StackedSVD.Prob.LinFormMoments
import StackedSVD.Prob.NoiseMoments
import StackedSVD.RMT.R4C

/-!
# The four-moment quadratic-form bound at a general noise law

Let `ν` be a general noise law (`NoiseLaw ν`, `Prob/NoiseLaw.lean`): mean 0, variance 1, finite
fourth moment `ν₄ = ∫ x ^ 4 ∂ν`. Let `y` carry the product law `Measure.pi fun _ : Fin D => ν`.
For a real matrix `B`, the quadratic form `y ⬝ᵥ (B *ᵥ y)` has mean `B.trace` and

  `E (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2 ≤ (ν₄ + 2) * ∑ a, ∑ b, B a b ^ 2`.

This is the general-law twin of the Gaussian file's rotation-invariance argument: at a general
law there is no rotation invariance, so the bound is proved by an induction on `D`, splitting
off one coordinate at a time exactly as `LinForm.moments` (`Prob/LinFormMoments.lean`) does for
the linear form. The complex form (`bil`, a local copy of `R2.bil`, `RMT/R2.lean:75`, placed
here so the general-law files do not import `RMT/R2.lean`) follows by writing the complex
matrix as its real and imaginary parts and adding the two real bounds.

This file imports none of `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean`,
`Vendor/COLT83/` or `ThetaEst.lean`.
-/

open MeasureTheory ProbabilityTheory
open scoped Matrix

namespace StackedSVD
namespace GenRMT

variable {ν : Measure ℝ}

/-! ### The bilinear form of a complex matrix at real vectors -/

/-- The bilinear form of a complex matrix at real vectors. Same definition as `R2.bil`
(`RMT/R2.lean:75`), placed here so that the general-law files do not import `RMT/R2.lean`;
a dedup is a later item. -/
noncomputable def bil {D : ℕ} (A : Matrix (Fin D) (Fin D) ℂ) (x y : Fin D → ℝ) : ℂ :=
  R4C.cvec x ⬝ᵥ (A *ᵥ R4C.cvec y)

theorem bil_eq_sum {D : ℕ} (A : Matrix (Fin D) (Fin D) ℂ) (x y : Fin D → ℝ) :
    bil A x y = ∑ a, ∑ b, A a b * (x a : ℂ) * (y b : ℂ) := by
  simp only [bil, R4C.cvec, dotProduct, Matrix.mulVec, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- `bil` of the complex resolvent is `cformC`, by definition. -/
theorem bil_resolvC {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x y : Fin D → ℝ) :
    bil (R4C.resolvC W z) x y = R4C.cformC W z x y := rfl

theorem bil_resolvC_mul {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x y : Fin D → ℝ) :
    bil (R4C.resolvC W z * R4C.resolvC W z) x y = R4C.cform2C W z x y := rfl

theorem measurable_qform {D : ℕ} (B : Matrix (Fin D) (Fin D) ℝ) :
    Measurable (fun y : Fin D → ℝ => y ⬝ᵥ (B *ᵥ y)) := by
  simp only [dotProduct, Matrix.mulVec]
  fun_prop

theorem measurable_bil_self {D : ℕ} (A : Matrix (Fin D) (Fin D) ℂ) :
    Measurable (fun y : Fin D → ℝ => bil A y y) := by
  simp only [bil, R4C.cvec, dotProduct, Matrix.mulVec]
  fun_prop

/-! ### The one-coordinate split, matrix algebra

`f_B(y) = y ⬝ᵥ (B *ᵥ y) - B.trace` splits, for `D + 1` coordinates written as `(y 0, y ∘ succ)`,
as `B 0 0 * (y 0 ^ 2 - 1) + y 0 * linForm b (y ∘ succ) + f_{B'}(y ∘ succ)`, with
`B' = B.submatrix succ succ` and `b j = B 0 j.succ + B j.succ 0`. These two lemmas give the
undelying pointwise identities for `y ⬝ᵥ (B *ᵥ y)` and for `B.trace`. -/

/-- `y ⬝ᵥ (B *ᵥ y) = B 0 0 * y 0 ^ 2 + y 0 * linForm b (y ∘ succ)
  + (y ∘ succ) ⬝ᵥ (B' *ᵥ (y ∘ succ))`, with `b j = B 0 j.succ + B j.succ 0`
and `B' = B.submatrix succ succ`. -/
private theorem dotProduct_mulVec_succ {D : ℕ} (B : Matrix (Fin (D + 1)) (Fin (D + 1)) ℝ)
    (y : Fin (D + 1) → ℝ) :
    y ⬝ᵥ (B *ᵥ y) = B 0 0 * y 0 ^ 2
      + y 0 * LinForm.linForm (fun j => B 0 j.succ + B j.succ 0) (fun j => y j.succ)
      + (fun j => y j.succ) ⬝ᵥ (B.submatrix Fin.succ Fin.succ *ᵥ fun j => y j.succ) := by
  have h0 : (B *ᵥ y) 0 = B 0 0 * y 0 + ∑ j : Fin D, B 0 j.succ * y j.succ := by
    simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_succ]
  have hi : ∀ i : Fin D, (B *ᵥ y) i.succ
      = B i.succ 0 * y 0 + (B.submatrix Fin.succ Fin.succ *ᵥ fun j => y j.succ) i := by
    intro i
    simp only [Matrix.mulVec, dotProduct, Matrix.submatrix_apply, Fin.sum_univ_succ]
  have hLsplit : y ⬝ᵥ (B *ᵥ y)
      = y 0 * (B *ᵥ y) 0 + ∑ i : Fin D, y i.succ * (B *ᵥ y) i.succ := by
    simp only [dotProduct, Fin.sum_univ_succ]
  have hSsplit : (∑ i : Fin D, y i.succ * (B *ᵥ y) i.succ)
      = (∑ i : Fin D, y i.succ * (B i.succ 0 * y 0))
        + (fun j => y j.succ) ⬝ᵥ (B.submatrix Fin.succ Fin.succ *ᵥ fun j => y j.succ) := by
    have hstep : ∀ i : Fin D, y i.succ * (B *ᵥ y) i.succ
        = y i.succ * (B i.succ 0 * y 0)
          + y i.succ * (B.submatrix Fin.succ Fin.succ *ᵥ fun j => y j.succ) i := by
      intro i; rw [hi i, mul_add]
    rw [Finset.sum_congr rfl fun i _ => hstep i, Finset.sum_add_distrib]
    rfl
  have hcombine : y 0 * (∑ j : Fin D, B 0 j.succ * y j.succ)
      + ∑ i : Fin D, y i.succ * (B i.succ 0 * y 0)
      = y 0 * LinForm.linForm (fun j => B 0 j.succ + B j.succ 0) (fun j => y j.succ) := by
    simp only [LinForm.linForm, Finset.mul_sum]
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hy0sq : y 0 * (B 0 0 * y 0) = B 0 0 * y 0 ^ 2 := by ring
  rw [hLsplit, h0, hSsplit, mul_add]
  linarith [hcombine, hy0sq]

/-- `B.trace = B 0 0 + B'.trace`, with `B' = B.submatrix succ succ`. -/
private theorem trace_succ {D : ℕ} (B : Matrix (Fin (D + 1)) (Fin (D + 1)) ℝ) :
    B.trace = B 0 0 + (B.submatrix Fin.succ Fin.succ).trace := by
  simp only [Matrix.trace, Matrix.diag, Matrix.submatrix_apply, Fin.sum_univ_succ]

/-- The same split for the sum of squared entries:
`∑ a, ∑ b, B a b ^ 2 = B 0 0 ^ 2 + ∑ j, B 0 j.succ ^ 2 + ∑ i, B i.succ 0 ^ 2
  + ∑ i, ∑ j, B i.succ j.succ ^ 2`. -/
private theorem sum_sq_succ {D : ℕ} (B : Matrix (Fin (D + 1)) (Fin (D + 1)) ℝ) :
    (∑ a, ∑ b, B a b ^ 2) = B 0 0 ^ 2 + (∑ j : Fin D, B 0 j.succ ^ 2)
      + (∑ i : Fin D, B i.succ 0 ^ 2) + ∑ i : Fin D, ∑ j : Fin D, B i.succ j.succ ^ 2 := by
  have h0 : (∑ b : Fin (D + 1), B 0 b ^ 2) = B 0 0 ^ 2 + ∑ j : Fin D, B 0 j.succ ^ 2 := by
    rw [Fin.sum_univ_succ]
  have hi : ∀ i : Fin D, (∑ b : Fin (D + 1), B i.succ b ^ 2)
      = B i.succ 0 ^ 2 + ∑ j : Fin D, B i.succ j.succ ^ 2 := by
    intro i; rw [Fin.sum_univ_succ]
  rw [Fin.sum_univ_succ, h0, Finset.sum_congr rfl fun i _ => hi i, Finset.sum_add_distrib]
  ring

/-! ### Transfer of the integral and of integrability across the one-coordinate split

`measurePreserving_piFinSuccAbove (fun _ : Fin (D + 1) => ν) 0` carries the product law on
`Fin (D + 1) → ℝ` to `ν.prod (Measure.pi fun _ : Fin D => ν)`, splitting off coordinate `0`.
These two lemmas transfer an integral (resp. an integrability statement) of a function
`G (y 0) (y ∘ succ)` across that equivalence, for a fully generic `G`; the matrix-specific
algebra is supplied separately by `dotProduct_mulVec_succ` and `trace_succ`. The simp set that
resolves the equivalence at a pair `z` is the one `Mathlib.MeasureTheory.Integral.Pi` uses for
the same equivalence (`MeasurableEquiv.piFinSuccAbove_symm_apply`, `Fin.insertNthEquiv`,
`Fin.insertNth_zero`, `Fin.cons_zero`, `Fin.cons_succ`). -/

private theorem integral_comp_succ {D : ℕ} [SigmaFinite ν] (G : ℝ → (Fin D → ℝ) → ℝ) :
    ∫ y, G (y 0) (fun j => y j.succ) ∂(Measure.pi fun _ : Fin (D + 1) => ν)
      = ∫ z : ℝ × (Fin D → ℝ), G z.1 z.2 ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) := by
  rw [← ((measurePreserving_piFinSuccAbove (fun _ : Fin (D + 1) => ν) 0).symm).integral_comp']
  refine integral_congr_ae (Filter.Eventually.of_forall fun z => ?_)
  simp only [MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv, Fin.insertNth_zero,
    Equiv.coe_fn_mk, Fin.cons_succ, Fin.cons_zero, Fin.zero_succAbove, cast_eq]

private theorem integrable_comp_succ {D : ℕ} [SigmaFinite ν] (G : ℝ → (Fin D → ℝ) → ℝ) :
    Integrable (fun y : Fin (D + 1) → ℝ => G (y 0) (fun j => y j.succ))
        (Measure.pi fun _ : Fin (D + 1) => ν)
      ↔ Integrable (fun z : ℝ × (Fin D → ℝ) => G z.1 z.2)
          (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
  rw [← ((measurePreserving_piFinSuccAbove (fun _ : Fin (D + 1) => ν) 0).symm).integrable_comp_emb
    (MeasurableEquiv.measurableEmbedding _)]
  have hpt : ∀ z : ℝ × (Fin D → ℝ),
      (fun y : Fin (D + 1) → ℝ => G (y 0) (fun j => y j.succ))
        ((MeasurableEquiv.piFinSuccAbove (fun _ : Fin (D + 1) => ℝ) 0).symm z) = G z.1 z.2 := by
    intro z
    simp only [MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv, Fin.insertNth_zero,
      Equiv.coe_fn_mk, Fin.cons_succ, Fin.cons_zero, Fin.zero_succAbove, cast_eq]
  simp only [Function.comp_def, hpt]

/-! ### Six-term integral splitting -/

private theorem integral_add₆ {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f₁ f₂ f₃ f₄ f₅ f₆ : α → ℝ} (h₁ : Integrable f₁ μ) (h₂ : Integrable f₂ μ)
    (h₃ : Integrable f₃ μ) (h₄ : Integrable f₄ μ) (h₅ : Integrable f₅ μ) (h₆ : Integrable f₆ μ) :
    ∫ z, (f₁ z + f₂ z + f₃ z + f₄ z + f₅ z + f₆ z) ∂μ
      = (∫ z, f₁ z ∂μ) + (∫ z, f₂ z ∂μ) + (∫ z, f₃ z ∂μ) + (∫ z, f₄ z ∂μ)
        + (∫ z, f₅ z ∂μ) + ∫ z, f₆ z ∂μ := by
  have e₁ : Integrable (fun z => f₁ z + f₂ z) μ := h₁.add h₂
  have e₂ : Integrable (fun z => f₁ z + f₂ z + f₃ z) μ := e₁.add h₃
  have e₃ : Integrable (fun z => f₁ z + f₂ z + f₃ z + f₄ z) μ := e₂.add h₄
  have e₄ : Integrable (fun z => f₁ z + f₂ z + f₃ z + f₄ z + f₅ z) μ := e₃.add h₅
  rw [integral_add e₄ h₆, integral_add e₃ h₅, integral_add e₂ h₄, integral_add e₁ h₃,
    integral_add h₁ h₂]

/-- The AM-GM bound `|a * c| ≤ (a ^ 2 + c ^ 2) / 2`, used to integrate the cross term between the
linear form and the tail quadratic form. -/
private theorem abs_mul_le_half_add_sq (a c : ℝ) : |a * c| ≤ (a ^ 2 + c ^ 2) / 2 := by
  rw [abs_le]
  constructor <;> nlinarith [sq_nonneg (a + c), sq_nonneg (a - c)]

/-! ### The induction -/

/-- The combined induction: the quadratic form `y ⬝ᵥ (B *ᵥ y) - B.trace` has mean `0`, is square
integrable, and its second moment is at most `(ν₄ + 2) * ∑ a, ∑ b, B a b ^ 2`. Proved by
induction on `D`, splitting off one coordinate with `integral_comp_succ` /
`integrable_comp_succ` and `dotProduct_mulVec_succ`, exactly as `LinForm.moments`
(`Prob/LinFormMoments.lean:144`) does for the linear form. -/
private theorem qformMoments (hν : NoiseLaw ν) :
    ∀ (D : ℕ) (B : Matrix (Fin D) (Fin D) ℝ),
      Integrable (fun y => y ⬝ᵥ (B *ᵥ y) - B.trace) (Measure.pi fun _ : Fin D => ν) ∧
      ∫ y, (y ⬝ᵥ (B *ᵥ y) - B.trace) ∂(Measure.pi fun _ : Fin D => ν) = 0 ∧
      Integrable (fun y => (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2) (Measure.pi fun _ : Fin D => ν) ∧
      ∫ y, (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
        ≤ ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B a b ^ 2 := by
  have hp := hν.prob
  intro D
  induction D with
  | zero =>
    intro B
    refine ⟨?_, ?_, ?_, ?_⟩ <;>
      simp [dotProduct, Matrix.mulVec, Matrix.trace, Matrix.diag]
  | succ D ih =>
    intro B
    set B' : Matrix (Fin D) (Fin D) ℝ := B.submatrix Fin.succ Fin.succ with hB'def
    set b : Fin D → ℝ := fun j => B 0 j.succ + B j.succ 0 with hbdef
    obtain ⟨hr1, hr0, hr2, hrbound⟩ := ih B'
    -- one-variable integrability, the `x`-side
    have hx0 : Integrable (fun _ : ℝ => (1 : ℝ)) ν := integrable_const 1
    have hx1 : Integrable (fun x : ℝ => x) ν := by
      simpa using hν.integrable_pow (k := 1) (by norm_num)
    have hx2 : Integrable (fun x : ℝ => x ^ 2) ν := hν.integrable_pow (k := 2) (by norm_num)
    have hx3 : Integrable (fun x : ℝ => x ^ 3) ν := hν.integrable_pow (k := 3) (by norm_num)
    have hx4 : Integrable (fun x : ℝ => x ^ 4) ν := hν.mom4
    have hp1 : Integrable (fun x : ℝ => B 0 0 * (x ^ 2 - 1)) ν := (hx2.sub hx0).const_mul _
    have hp1sq : Integrable (fun x : ℝ => (B 0 0 * (x ^ 2 - 1)) ^ 2) ν := by
      have heq : (fun x : ℝ => (B 0 0 * (x ^ 2 - 1)) ^ 2)
          = fun x => B 0 0 ^ 2 * (x ^ 4 - 2 * x ^ 2 + 1) := by funext x; ring
      rw [heq]
      exact ((hx4.sub (hx2.const_mul 2)).add hx0).const_mul _
    have hpq1 : Integrable (fun x : ℝ => B 0 0 * (x ^ 2 - 1) * x) ν := by
      have heq : (fun x : ℝ => B 0 0 * (x ^ 2 - 1) * x) = fun x => B 0 0 * (x ^ 3 - x) := by
        funext x; ring
      rw [heq]
      exact (hx3.sub hx1).const_mul _
    -- one-variable integrability, the `u`-side
    have hb4 : Integrable (fun u => LinForm.linForm b u ^ 4) (Measure.pi fun _ : Fin D => ν) :=
      LinForm.integrable_linForm_pow_four hν.mean hν.var hν.mom4 b
    have hqu : Integrable (LinForm.linForm b) (Measure.pi fun _ : Fin D => ν) := by
      simpa using LinForm.integrable_linForm_pow_of_four hb4 (by norm_num : (1 : ℕ) ≤ 4)
    have hqusq : Integrable (fun u => LinForm.linForm b u ^ 2) (Measure.pi fun _ : Fin D => ν) :=
      LinForm.integrable_linForm_pow_of_four hb4 (by norm_num : (2 : ℕ) ≤ 4)
    have hqur1 : Integrable (fun u : Fin D → ℝ => LinForm.linForm b u * (u ⬝ᵥ (B' *ᵥ u) - B'.trace))
        (Measure.pi fun _ : Fin D => ν) := by
      have hg : Integrable (fun u : Fin D → ℝ =>
          (LinForm.linForm b u ^ 2 + (u ⬝ᵥ (B' *ᵥ u) - B'.trace) ^ 2) / 2)
          (Measure.pi fun _ : Fin D => ν) := (hqusq.add hr2).div_const 2
      have hmeas : Measurable (fun u : Fin D → ℝ =>
          LinForm.linForm b u * (u ⬝ᵥ (B' *ᵥ u) - B'.trace)) := by
        have h1 : Measurable (LinForm.linForm b) := LinForm.measurable_linForm b
        have h2 : Measurable (fun u : Fin D → ℝ => u ⬝ᵥ (B' *ᵥ u) - B'.trace) := by
          simp only [dotProduct, Matrix.mulVec]
          fun_prop
        exact h1.mul h2
      refine Integrable.mono' hg hmeas.aestronglyMeasurable
        (Filter.Eventually.of_forall fun u => ?_)
      rw [Real.norm_eq_abs]
      exact abs_mul_le_half_add_sq _ _
    -- integrability on the product measure
    have hPI : Integrable (fun z : ℝ × (Fin D → ℝ) => B 0 0 * (z.1 ^ 2 - 1))
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      simpa using hp1.mul_prod (integrable_const (1 : ℝ) (μ := Measure.pi fun _ : Fin D => ν))
    have hQI : Integrable (fun z : ℝ × (Fin D → ℝ) => z.1 * LinForm.linForm b z.2)
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := hx1.mul_prod hqu
    have hRI : Integrable (fun z : ℝ × (Fin D → ℝ) => z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      simpa using (integrable_const (1 : ℝ) (μ := ν)).mul_prod hr1
    have hPI2 : Integrable (fun z : ℝ × (Fin D → ℝ) => (B 0 0 * (z.1 ^ 2 - 1)) ^ 2)
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      simpa using hp1sq.mul_prod (integrable_const (1 : ℝ) (μ := Measure.pi fun _ : Fin D => ν))
    have hQI2 : Integrable (fun z : ℝ × (Fin D → ℝ) => (z.1 * LinForm.linForm b z.2) ^ 2)
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      have heq : (fun z : ℝ × (Fin D → ℝ) => (z.1 * LinForm.linForm b z.2) ^ 2)
          = fun z => z.1 ^ 2 * LinForm.linForm b z.2 ^ 2 := by funext z; ring
      rw [heq]
      exact hx2.mul_prod hqusq
    have hRI2 : Integrable (fun z : ℝ × (Fin D → ℝ) => (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2)
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      simpa using (integrable_const (1 : ℝ) (μ := ν)).mul_prod hr2
    have hPQI : Integrable (fun z : ℝ × (Fin D → ℝ) =>
        2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2)))
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      have heq : (fun z : ℝ × (Fin D → ℝ) =>
          2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2)))
          = fun z => 2 * (B 0 0 * (z.1 ^ 2 - 1) * z.1) * LinForm.linForm b z.2 := by
        funext z; ring
      rw [heq]
      exact (hpq1.const_mul 2).mul_prod hqu
    have hPRI : Integrable (fun z : ℝ × (Fin D → ℝ) =>
        2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)))
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      have h := (hp1.const_mul 2).mul_prod hr1
      simpa [mul_assoc] using h
    have hQRI : Integrable (fun z : ℝ × (Fin D → ℝ) =>
        2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)))
        (ν.prod (Measure.pi fun _ : Fin D => ν)) := by
      have heq : (fun z : ℝ × (Fin D → ℝ) =>
          2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)))
          = fun z => (2 * z.1) * (LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)) := by
        funext z; ring
      rw [heq]
      exact (hx1.const_mul 2).mul_prod hqur1
    -- the pointwise split
    have hsplit : ∀ y : Fin (D + 1) → ℝ,
        y ⬝ᵥ (B *ᵥ y) - B.trace
          = B 0 0 * (y 0 ^ 2 - 1) + y 0 * LinForm.linForm b (fun j => y j.succ)
            + ((fun j => y j.succ) ⬝ᵥ (B' *ᵥ fun j => y j.succ) - B'.trace) := by
      intro y
      rw [dotProduct_mulVec_succ, trace_succ]
      ring
    have hfun1 : (fun y : Fin (D + 1) → ℝ => y ⬝ᵥ (B *ᵥ y) - B.trace)
        = fun y => B 0 0 * (y 0 ^ 2 - 1) + y 0 * LinForm.linForm b (fun j => y j.succ)
            + ((fun j => y j.succ) ⬝ᵥ (B' *ᵥ fun j => y j.succ) - B'.trace) := funext hsplit
    have hfun2 : (fun y : Fin (D + 1) → ℝ => (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2)
        = fun y => (B 0 0 * (y 0 ^ 2 - 1) + y 0 * LinForm.linForm b (fun j => y j.succ)
            + ((fun j => y j.succ) ⬝ᵥ (B' *ᵥ fun j => y j.succ) - B'.trace)) ^ 2 := by
      funext y; rw [hsplit y]
    -- integrability, first moment
    have hI1 : Integrable (fun y : Fin (D + 1) → ℝ => y ⬝ᵥ (B *ᵥ y) - B.trace)
        (Measure.pi fun _ : Fin (D + 1) => ν) := by
      rw [hfun1, integrable_comp_succ (fun x u => B 0 0 * (x ^ 2 - 1)
        + x * LinForm.linForm b u + (u ⬝ᵥ (B' *ᵥ u) - B'.trace))]
      exact (hPI.add hQI).add hRI
    -- integrability, second moment
    have hI2 : Integrable (fun y : Fin (D + 1) → ℝ => (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2)
        (Measure.pi fun _ : Fin (D + 1) => ν) := by
      rw [hfun2, integrable_comp_succ (fun x u => (B 0 0 * (x ^ 2 - 1)
        + x * LinForm.linForm b u + (u ⬝ᵥ (B' *ᵥ u) - B'.trace)) ^ 2)]
      have heq : (fun z : ℝ × (Fin D → ℝ) =>
          (B 0 0 * (z.1 ^ 2 - 1) + z.1 * LinForm.linForm b z.2
            + (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)) ^ 2)
          = fun z => (B 0 0 * (z.1 ^ 2 - 1)) ^ 2 + (z.1 * LinForm.linForm b z.2) ^ 2
            + (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2
            + 2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2))
            + 2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace))
            + 2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)) := by
        funext z; ring
      rw [heq]
      exact ((((hPI2.add hQI2).add hRI2).add hPQI).add hPRI).add hQRI
    -- the mean
    have hmean0 : ∫ y : Fin (D + 1) → ℝ, y ⬝ᵥ (B *ᵥ y) - B.trace
        ∂(Measure.pi fun _ : Fin (D + 1) => ν) = 0 := by
      rw [hfun1, integral_comp_succ (fun x u => B 0 0 * (x ^ 2 - 1)
        + x * LinForm.linForm b u + (u ⬝ᵥ (B' *ᵥ u) - B'.trace))]
      have hPQint : Integrable (fun z : ℝ × (Fin D → ℝ) =>
          B 0 0 * (z.1 ^ 2 - 1) + z.1 * LinForm.linForm b z.2)
          (ν.prod (Measure.pi fun _ : Fin D => ν)) := hPI.add hQI
      rw [integral_add hPQint hRI, integral_add hPI hQI]
      have hIP : ∫ z : ℝ × (Fin D → ℝ), B 0 0 * (z.1 ^ 2 - 1)
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) => B 0 0 * (z.1 ^ 2 - 1))
              = fun z => (B 0 0 * (z.1 ^ 2 - 1)) * (1 : (Fin D → ℝ) → ℝ) z.2 from by
            funext z; simp,
          integral_prod_mul (fun x : ℝ => B 0 0 * (x ^ 2 - 1)) (1 : (Fin D → ℝ) → ℝ),
          integral_const_mul, integral_sub hx2 hx0, hν.var]
        simp
      have hIQ : ∫ z : ℝ × (Fin D → ℝ), z.1 * LinForm.linForm b z.2
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [integral_prod_mul (fun x : ℝ => x) (LinForm.linForm b), hν.mean, zero_mul]
      have hIR : ∫ z : ℝ × (Fin D → ℝ), z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) => z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)
              = fun z => (1 : ℝ) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) from by funext z; ring,
          integral_prod_mul (fun _ : ℝ => (1 : ℝ)) (fun u => u ⬝ᵥ (B' *ᵥ u) - B'.trace)]
        simp [hr0]
      rw [hIP, hIQ, hIR]
      ring
    -- the bound
    have hbound : ∫ y : Fin (D + 1) → ℝ, (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2
        ∂(Measure.pi fun _ : Fin (D + 1) => ν)
        ≤ ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B a b ^ 2 := by
      rw [hfun2, integral_comp_succ (fun x u => (B 0 0 * (x ^ 2 - 1)
        + x * LinForm.linForm b u + (u ⬝ᵥ (B' *ᵥ u) - B'.trace)) ^ 2)]
      have heq : (fun z : ℝ × (Fin D → ℝ) =>
          (B 0 0 * (z.1 ^ 2 - 1) + z.1 * LinForm.linForm b z.2
            + (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)) ^ 2)
          = fun z => (B 0 0 * (z.1 ^ 2 - 1)) ^ 2 + (z.1 * LinForm.linForm b z.2) ^ 2
            + (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2
            + 2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2))
            + 2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace))
            + 2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)) := by
        funext z; ring
      rw [heq, integral_add₆ hPI2 hQI2 hRI2 hPQI hPRI hQRI]
      have hIP2 : ∫ z : ℝ × (Fin D → ℝ), (B 0 0 * (z.1 ^ 2 - 1)) ^ 2
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = B 0 0 ^ 2 * ((∫ x, x ^ 4 ∂ν) - 1) := by
        rw [show (fun z : ℝ × (Fin D → ℝ) => (B 0 0 * (z.1 ^ 2 - 1)) ^ 2)
              = fun z => (B 0 0 * (z.1 ^ 2 - 1)) ^ 2 * (1 : (Fin D → ℝ) → ℝ) z.2 from by
            funext z; simp,
          integral_prod_mul (fun x : ℝ => (B 0 0 * (x ^ 2 - 1)) ^ 2) (1 : (Fin D → ℝ) → ℝ)]
        have hex : ∫ x : ℝ, (B 0 0 * (x ^ 2 - 1)) ^ 2 ∂ν
            = B 0 0 ^ 2 * ((∫ x, x ^ 4 ∂ν) - 2 * (∫ x, x ^ 2 ∂ν) + 1) := by
          have heq2 : (fun x : ℝ => (B 0 0 * (x ^ 2 - 1)) ^ 2)
              = fun x => B 0 0 ^ 2 * (x ^ 4 - 2 * x ^ 2 + 1) := by funext x; ring
          have hsub42 : Integrable (fun x : ℝ => x ^ 4 - 2 * x ^ 2) ν := hx4.sub (hx2.const_mul 2)
          rw [heq2, integral_const_mul, integral_add hsub42 hx0,
            integral_sub hx4 (hx2.const_mul 2), integral_const_mul]
          simp
        have h1u : ∫ _u : Fin D → ℝ, (1 : (Fin D → ℝ) → ℝ) _u
            ∂(Measure.pi fun _ : Fin D => ν) = 1 := by simp
        rw [hex, hν.var, h1u, mul_one]
        ring
      have hIQ2 : ∫ z : ℝ × (Fin D → ℝ), (z.1 * LinForm.linForm b z.2) ^ 2
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = ∑ j, b j ^ 2 := by
        have heq2 : (fun z : ℝ × (Fin D → ℝ) => (z.1 * LinForm.linForm b z.2) ^ 2)
            = fun z => z.1 ^ 2 * LinForm.linForm b z.2 ^ 2 := by funext z; ring
        rw [heq2, integral_prod_mul (fun x : ℝ => x ^ 2) (fun u => LinForm.linForm b u ^ 2),
          hν.var, LinForm.integral_linForm_sq hν.mean hν.var hν.mom4 b]
        ring
      have hIR2 : ∫ z : ℝ × (Fin D → ℝ), (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν))
          ≤ ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B' a b ^ 2 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) => (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2)
              = fun z => (1 : ℝ) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2 from by funext z; ring,
          integral_prod_mul (fun _ : ℝ => (1 : ℝ)) (fun u => (u ⬝ᵥ (B' *ᵥ u) - B'.trace) ^ 2)]
        simpa using hrbound
      have hIPQ : ∫ z : ℝ × (Fin D → ℝ),
          2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2))
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) =>
              2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.1 * LinForm.linForm b z.2)))
              = fun z => (2 * (B 0 0 * (z.1 ^ 2 - 1) * z.1)) * LinForm.linForm b z.2 from by
            funext z; ring,
          integral_prod_mul (fun x : ℝ => 2 * (B 0 0 * (x ^ 2 - 1) * x)) (LinForm.linForm b),
          LinForm.integral_linForm hν.mean hν.var hν.mom4 b, mul_zero]
      have hIPR : ∫ z : ℝ × (Fin D → ℝ),
          2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace))
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) =>
              2 * (B 0 0 * (z.1 ^ 2 - 1) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)))
              = fun z => (2 * (B 0 0 * (z.1 ^ 2 - 1))) * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) from by
            funext z; ring,
          integral_prod_mul (fun x : ℝ => 2 * (B 0 0 * (x ^ 2 - 1)))
            (fun u => u ⬝ᵥ (B' *ᵥ u) - B'.trace), hr0, mul_zero]
      have hIQR : ∫ z : ℝ × (Fin D → ℝ),
          2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace))
          ∂(ν.prod (Measure.pi fun _ : Fin D => ν)) = 0 := by
        rw [show (fun z : ℝ × (Fin D → ℝ) =>
              2 * (z.1 * LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace)))
              = fun z => (2 * z.1) * (LinForm.linForm b z.2 * (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace))
              from by funext z; ring,
          integral_prod_mul (fun x : ℝ => 2 * x)
            (fun u => LinForm.linForm b u * (u ⬝ᵥ (B' *ᵥ u) - B'.trace)),
          integral_const_mul, hν.mean, mul_zero, zero_mul]
      rw [hIP2, hIQ2, hIPQ, hIPR, hIQR]
      have hsum := sum_sq_succ B
      have hbsq_le : ∑ j : Fin D, b j ^ 2
          ≤ 2 * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2) := by
        have hpt : ∀ j : Fin D, b j ^ 2 ≤ 2 * (B 0 j.succ ^ 2 + B j.succ 0 ^ 2) := by
          intro j
          simp only [hbdef]
          nlinarith [sq_nonneg (B 0 j.succ - B j.succ 0)]
        calc ∑ j, b j ^ 2 ≤ ∑ j : Fin D, 2 * (B 0 j.succ ^ 2 + B j.succ 0 ^ 2) :=
              Finset.sum_le_sum fun j _ => hpt j
          _ = 2 * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2) := by
              rw [← Finset.mul_sum, Finset.sum_add_distrib]
      have hnu4nn : (0 : ℝ) ≤ ∫ x, x ^ 4 ∂ν := integral_pow_four_nonneg ν
      have hXnn : (0 : ℝ) ≤ ∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2 := by
        positivity
      have hstep1 : ∑ j : Fin D, b j ^ 2
          ≤ ((∫ x, x ^ 4 ∂ν) + 2)
            * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2) := by
        calc ∑ j, b j ^ 2
            ≤ 2 * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2) := hbsq_le
          _ ≤ ((∫ x, x ^ 4 ∂ν) + 2)
              * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2) :=
              mul_le_mul_of_nonneg_right (by linarith) hXnn
      have hstep2 : B 0 0 ^ 2 * ((∫ x, x ^ 4 ∂ν) - 1) ≤ ((∫ x, x ^ 4 ∂ν) + 2) * B 0 0 ^ 2 := by
        nlinarith [sq_nonneg (B 0 0)]
      have hB'sum : ∑ a, ∑ b, B' a b ^ 2
          = ∑ i : Fin D, ∑ j : Fin D, B i.succ j.succ ^ 2 := by
        simp only [hB'def, Matrix.submatrix_apply]
      calc B 0 0 ^ 2 * ((∫ x, x ^ 4 ∂ν) - 1) + ∑ j, b j ^ 2
              + ∫ z : ℝ × (Fin D → ℝ), (z.2 ⬝ᵥ (B' *ᵥ z.2) - B'.trace) ^ 2
                ∂(ν.prod (Measure.pi fun _ : Fin D => ν))
              + 0 + 0 + 0
          ≤ B 0 0 ^ 2 * ((∫ x, x ^ 4 ∂ν) - 1) + ∑ j, b j ^ 2
              + ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B' a b ^ 2 := by linarith [hIR2]
        _ ≤ ((∫ x, x ^ 4 ∂ν) + 2) * B 0 0 ^ 2
              + ((∫ x, x ^ 4 ∂ν) + 2) * (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2)
              + ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B' a b ^ 2 := by linarith [hstep1, hstep2]
        _ = ((∫ x, x ^ 4 ∂ν) + 2)
              * (B 0 0 ^ 2 + (∑ j : Fin D, B 0 j.succ ^ 2 + ∑ i : Fin D, B i.succ 0 ^ 2)
                + ∑ a, ∑ b, B' a b ^ 2) := by ring
        _ = ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B a b ^ 2 := by rw [hB'sum, hsum]; ring
    exact ⟨hI1, hmean0, hI2, hbound⟩

/-! ### The public real-form theorems -/

/-- **The real four-moment bound.** For `y` with i.i.d. coordinates of law `ν` (mean 0,
variance 1, fourth moment `ν₄`), the quadratic form `yᵀ B y` has mean `tr B` and
`E (yᵀ B y - tr B)² ≤ (ν₄ + 2) ∑_{a,b} B_ab²`. -/
theorem integral_sq_qform_sub_trace_le {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (B : Matrix (Fin D) (Fin D) ℝ) :
    Integrable (fun y => (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2) (Measure.pi fun _ : Fin D => ν) ∧
    ∫ y, (y ⬝ᵥ (B *ᵥ y) - B.trace) ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, B a b ^ 2 := by
  obtain ⟨_, _, hI2, hbound⟩ := qformMoments hν D B
  exact ⟨hI2, hbound⟩

theorem integral_qform {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ} (B : Matrix (Fin D) (Fin D) ℝ) :
    ∫ y, y ⬝ᵥ (B *ᵥ y) ∂(Measure.pi fun _ : Fin D => ν) = B.trace := by
  have hp := hν.prob
  obtain ⟨hI1, hmean0, _, _⟩ := qformMoments hν D B
  have heq : (fun y : Fin D → ℝ => y ⬝ᵥ (B *ᵥ y))
      = fun y => (y ⬝ᵥ (B *ᵥ y) - B.trace) + B.trace := by funext y; ring
  rw [heq, integral_add hI1 (integrable_const B.trace), hmean0, integral_const]
  simp

/-! ### The complex bound -/

/-- `bil A y y - A.trace` splits into its real and imaginary parts, each the real quadratic
form of the real and imaginary parts of `A`. -/
private theorem bil_self_sub_trace_eq {D : ℕ} (A : Matrix (Fin D) (Fin D) ℂ) (y : Fin D → ℝ) :
    bil A y y - A.trace
      = ((y ⬝ᵥ (A.map Complex.re *ᵥ y) - (A.map Complex.re).trace : ℝ) : ℂ)
        + ((y ⬝ᵥ (A.map Complex.im *ᵥ y) - (A.map Complex.im).trace : ℝ) : ℂ) * Complex.I := by
  have hbilsum : bil A y y = (∑ a, ∑ b, ((A a b).re : ℂ) * (y a : ℂ) * (y b : ℂ))
      + (∑ a, ∑ b, ((A a b).im : ℂ) * (y a : ℂ) * (y b : ℂ)) * Complex.I := by
    rw [bil_eq_sum, Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun b _ => ?_
    have hab : A a b = ((A a b).re : ℂ) + ((A a b).im : ℂ) * Complex.I := (Complex.re_add_im _).symm
    linear_combination (y a : ℂ) * (y b : ℂ) * hab
  have htracesum : A.trace = (∑ a, ((A a a).re : ℂ)) + (∑ a, ((A a a).im : ℂ)) * Complex.I := by
    simp only [Matrix.trace, Matrix.diag, Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun a _ => ?_
    exact (Complex.re_add_im _).symm
  have hre0 : y ⬝ᵥ (A.map Complex.re *ᵥ y) = ∑ a, ∑ b, (A a b).re * y a * y b := by
    simp only [dotProduct, Matrix.mulVec, Matrix.map_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have him0 : y ⬝ᵥ (A.map Complex.im *ᵥ y) = ∑ a, ∑ b, (A a b).im * y a * y b := by
    simp only [dotProduct, Matrix.mulVec, Matrix.map_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have htrRe : (A.map Complex.re).trace = ∑ a, (A a a).re := by
    simp only [Matrix.trace, Matrix.diag, Matrix.map_apply]
  have htrIm : (A.map Complex.im).trace = ∑ a, (A a a).im := by
    simp only [Matrix.trace, Matrix.diag, Matrix.map_apply]
  rw [hbilsum, htracesum, hre0, him0, htrRe, htrIm]
  push_cast
  ring

/-- `‖z‖ ^ 2` in terms of the real and imaginary parts. -/
private theorem normSq_eq_add_sq (z : ℂ) : ‖z‖ ^ 2 = z.re ^ 2 + z.im ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]; ring

/-- **The four-moment bound, complex form.** -/
theorem integral_normSq_bil_sub_trace_le {ν : Measure ℝ} (hν : NoiseLaw ν) {D : ℕ}
    (A : Matrix (Fin D) (Fin D) ℂ) :
    Integrable (fun y => ‖bil A y y - A.trace‖ ^ 2) (Measure.pi fun _ : Fin D => ν) ∧
    ∫ y, ‖bil A y y - A.trace‖ ^ 2 ∂(Measure.pi fun _ : Fin D => ν)
      ≤ ((∫ x, x ^ 4 ∂ν) + 2) * ∑ a, ∑ b, ‖A a b‖ ^ 2 := by
  set ReA : Matrix (Fin D) (Fin D) ℝ := A.map Complex.re with hReAdef
  set ImA : Matrix (Fin D) (Fin D) ℝ := A.map Complex.im with hImAdef
  have hfun : (fun y : Fin D → ℝ => ‖bil A y y - A.trace‖ ^ 2)
      = fun y => (y ⬝ᵥ (ReA *ᵥ y) - ReA.trace) ^ 2 + (y ⬝ᵥ (ImA *ᵥ y) - ImA.trace) ^ 2 := by
    funext y
    rw [bil_self_sub_trace_eq, normSq_eq_add_sq]
    simp [hReAdef, hImAdef]
  obtain ⟨hIre, hbre⟩ := integral_sq_qform_sub_trace_le hν ReA
  obtain ⟨hIim, hbim⟩ := integral_sq_qform_sub_trace_le hν ImA
  constructor
  · rw [hfun]; exact hIre.add hIim
  · rw [hfun, integral_add hIre hIim]
    have hentry : ∀ a b : Fin D, ‖A a b‖ ^ 2 = ReA a b ^ 2 + ImA a b ^ 2 := by
      intro a b
      simp only [hReAdef, hImAdef, Matrix.map_apply]
      exact normSq_eq_add_sq (A a b)
    have hsum : (∑ a, ∑ b, ‖A a b‖ ^ 2) = (∑ a, ∑ b, ReA a b ^ 2) + ∑ a, ∑ b, ImA a b ^ 2 := by
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun a _ => ?_
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun b _ => hentry a b
    rw [hsum, mul_add]
    linarith [hbre, hbim]

end GenRMT
end StackedSVD
