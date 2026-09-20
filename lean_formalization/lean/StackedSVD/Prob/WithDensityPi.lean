/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar
import Mathlib.MeasureTheory.Measure.WithDensity
import Mathlib.LinearAlgebra.Matrix.ToLin
import Mathlib.LinearAlgebra.Determinant

/-!
# Densities through products and linear changes of variables (L5, units U6 and U2)

Two gaps in Mathlib v4.33.0 that the Gaussian marginalization of `app:wstacksvd_mle`
(`notes/archive/L5_mle_marginal.md`) needs.

* `pi_withDensity`: a finite product of measures with densities has the product density.
  The proof is `Measure.pi_eq` plus `Measure.restrict_pi_pi` plus a Fubini for `lintegral`
  over `Fin n` coordinates, `lintegral_fin_nat_prod_eq_prod`, which mirrors Mathlib's
  `MeasureTheory.integral_fin_nat_prod_eq_prod` with `lintegral_prod_mul` in place of
  `integral_prod_mul`.
* `withDensity_map_mulVec`: on `Fin d → ℝ`, the pushforward of `volume.withDensity f` through
  an invertible matrix `A` has density `|det A|⁻¹ * f (A⁻¹ *ᵥ ·)`. The ingredients are
  `withDensity_map_equiv` (a density through a measurable equivalence) and `volume_map_mulVec`
  (`map_linearMap_addHaar_pi_eq_smul_addHaar` with `LinearMap.det_toLin'`).
-/

open MeasureTheory Set
open scoped ENNReal Matrix

namespace StackedSVD

/-- Fubini for a product of `ℝ≥0∞`-valued functions of the coordinates, `Fin n` version.
Same induction as Mathlib's `MeasureTheory.integral_fin_nat_prod_eq_prod`
(`Mathlib/MeasureTheory/Integral/Pi.lean`), with `lintegral_prod_mul` in place of
`integral_prod_mul`. -/
theorem lintegral_fin_nat_prod_eq_prod {n : ℕ} {E : Fin n → Type*}
    {mE : ∀ i, MeasurableSpace (E i)} {μ : (i : Fin n) → Measure (E i)} [∀ i, SigmaFinite (μ i)]
    (f : (i : Fin n) → E i → ℝ≥0∞) (hf : ∀ i, Measurable (f i)) :
    ∫⁻ x : (i : Fin n) → E i, ∏ i, f i (x i) ∂(Measure.pi μ) = ∏ i, ∫⁻ x, f i x ∂(μ i) := by
  induction n with
  | zero => simp
  | succ n n_ih =>
    have hprod : Measurable fun x : (i : Fin n) → E (Fin.succ i) => ∏ i, f (Fin.succ i) (x i) :=
      Finset.measurable_prod _ fun i _ => (hf i.succ).comp (measurable_pi_apply i)
    calc ∫⁻ x : (i : Fin (n + 1)) → E i, ∏ i, f i (x i) ∂(Measure.pi μ)
        = ∫⁻ x : E 0 × ((i : Fin n) → E ((0 : Fin (n + 1)).succAbove i)),
            ∏ i, f i ((MeasurableEquiv.piFinSuccAbove E 0).symm x i)
            ∂((μ 0).prod (Measure.pi fun i => μ ((0 : Fin (n + 1)).succAbove i))) := by
          rw [← ((measurePreserving_piFinSuccAbove μ 0).symm).map_eq, lintegral_map_equiv]
      _ = ∫⁻ x : E 0 × ((i : Fin n) → E (Fin.succ i)),
            f 0 x.1 * ∏ i : Fin n, f (Fin.succ i) (x.2 i)
            ∂((μ 0).prod (Measure.pi fun i => μ i.succ)) := by
          refine lintegral_congr fun x => ?_
          simp only [MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv,
            Equiv.coe_fn_mk, Fin.insertNth_zero, Fin.prod_univ_succ, Fin.cons_zero,
            Fin.cons_succ]
          -- the remaining `cast` is along `E ((0 : Fin (n + 1)).succAbove i) = E i.succ`,
          -- a reflexivity proof, so both sides are definitionally equal
          rfl
      _ = (∫⁻ x, f 0 x ∂(μ 0)) * ∏ i : Fin n, ∫⁻ x, f (Fin.succ i) x ∂(μ i.succ) := by
          rw [lintegral_prod_mul (hf 0).aemeasurable hprod.aemeasurable,
            n_ih (fun i => f i.succ) fun i => hf i.succ]
      _ = ∏ i, ∫⁻ x, f i x ∂(μ i) := by rw [Fin.prod_univ_succ]

/-- U6: a finite product of measures with densities has the product density. -/
theorem pi_withDensity {n : ℕ} {E : Fin n → Type*} [∀ i, MeasurableSpace (E i)]
    (μ : (i : Fin n) → Measure (E i)) [∀ i, SigmaFinite (μ i)]
    (f : (i : Fin n) → E i → ℝ≥0∞) (hf : ∀ i, Measurable (f i))
    [∀ i, SigmaFinite ((μ i).withDensity (f i))] :
    Measure.pi (fun i => (μ i).withDensity (f i))
      = (Measure.pi μ).withDensity (fun x => ∏ i, f i (x i)) := by
  refine Measure.pi_eq fun s hs => ?_
  rw [withDensity_apply _ (MeasurableSet.univ_pi hs), Measure.restrict_pi_pi,
    lintegral_fin_nat_prod_eq_prod f hf]
  exact Finset.prod_congr rfl fun i _ => (withDensity_apply _ (hs i)).symm

/-- U2: a density transported through a measurable equivalence. -/
theorem withDensity_map_equiv {α β : Type*} [MeasurableSpace α] [MeasurableSpace β]
    (μ : Measure α) (e : α ≃ᵐ β) {f : α → ℝ≥0∞} (hf : Measurable f) :
    (μ.withDensity f).map e = (μ.map e).withDensity (fun y => f (e.symm y)) := by
  have hfe : Measurable fun y : β => f (e.symm y) := hf.comp e.symm.measurable
  ext s hs
  rw [Measure.map_apply e.measurable hs, withDensity_apply _ (e.measurable hs),
    withDensity_apply _ hs, Measure.restrict_map e.measurable hs,
    lintegral_map hfe e.measurable]
  simp only [MeasurableEquiv.symm_apply_apply]

/-- Lebesgue measure on `Fin d → ℝ` through an invertible matrix
(`map_linearMap_addHaar_pi_eq_smul_addHaar` at `Matrix.toLin' A`, whose determinant is
`A.det` by `LinearMap.det_toLin'`). -/
theorem volume_map_mulVec {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0) :
    (volume : Measure (Fin d → ℝ)).map (fun x => A *ᵥ x)
      = ENNReal.ofReal |A.det|⁻¹ • (volume : Measure (Fin d → ℝ)) := by
  have hdet : LinearMap.det (Matrix.toLin' A) ≠ 0 := by
    rw [LinearMap.det_toLin']; exact hA
  have hcoe : ⇑(Matrix.toLin' A) = fun x => A *ᵥ x := by
    funext x; exact Matrix.toLin'_apply A x
  have h := Measure.map_linearMap_addHaar_pi_eq_smul_addHaar hdet (volume : Measure (Fin d → ℝ))
  rw [hcoe, LinearMap.det_toLin'] at h
  rw [h, abs_inv]

/-- `x ↦ A *ᵥ x` as a measurable equivalence of `Fin d → ℝ`, for `A.det ≠ 0`. -/
noncomputable def mulVecEquiv {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0) :
    (Fin d → ℝ) ≃ᵐ (Fin d → ℝ) where
  toFun x := A *ᵥ x
  invFun y := A⁻¹ *ᵥ y
  left_inv x := by
    change A⁻¹ *ᵥ (A *ᵥ x) = x
    rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr hA),
      Matrix.one_mulVec]
  right_inv y := by
    change A *ᵥ (A⁻¹ *ᵥ y) = y
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.mpr hA),
      Matrix.one_mulVec]
  measurable_toFun := (Matrix.mulVecLin A).continuous_of_finiteDimensional.measurable
  measurable_invFun := (Matrix.mulVecLin A⁻¹).continuous_of_finiteDimensional.measurable

@[simp] theorem mulVecEquiv_apply {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0)
    (x : Fin d → ℝ) : mulVecEquiv A hA x = A *ᵥ x := rfl

@[simp] theorem mulVecEquiv_symm_apply {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0)
    (y : Fin d → ℝ) : (mulVecEquiv A hA).symm y = A⁻¹ *ᵥ y := rfl

/-- U2: the pushforward of a density on `Fin d → ℝ` through an invertible matrix. -/
theorem withDensity_map_mulVec {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.det ≠ 0)
    {f : (Fin d → ℝ) → ℝ≥0∞} (hf : Measurable f) :
    ((volume : Measure (Fin d → ℝ)).withDensity f).map (fun x => A *ᵥ x)
      = (volume : Measure (Fin d → ℝ)).withDensity
          (fun y => ENNReal.ofReal |A.det|⁻¹ * f (A⁻¹ *ᵥ y)) := by
  have hmap : (volume : Measure (Fin d → ℝ)).map (mulVecEquiv A hA)
      = ENNReal.ofReal |A.det|⁻¹ • (volume : Measure (Fin d → ℝ)) := volume_map_mulVec A hA
  have hcoe : ((volume : Measure (Fin d → ℝ)).withDensity f).map (fun x => A *ᵥ x)
      = ((volume : Measure (Fin d → ℝ)).withDensity f).map (mulVecEquiv A hA) := rfl
  have hm : Measurable fun y : Fin d → ℝ => f ((mulVecEquiv A hA).symm y) :=
    hf.comp (mulVecEquiv A hA).symm.measurable
  rw [hcoe, withDensity_map_equiv _ (mulVecEquiv A hA) hf, hmap, withDensity_smul_measure,
    ← withDensity_smul _ hm]
  rfl

end StackedSVD
