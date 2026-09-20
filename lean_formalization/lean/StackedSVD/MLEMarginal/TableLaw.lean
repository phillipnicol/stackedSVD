/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLEMarginal.RowLaw

/-!
# The law of a table and of the `M` tables (L5, unit U5 and the density of the joint law)

* `reTableLaw_eq_pi`: the rows of a table are independent, `reTableLaw = Measure.pi (reRowLaw)`
  (`Measure.pi_map_pi` and `measurePreserving_arrowProdEquivProdArrow`).
* `reTableLaw_eq_withDensity`, `reJointLaw_eq_withDensity`: `pi_withDensity` twice, over the
  rows and over the tables; the joint law has Lebesgue density `reDensity`.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped ENNReal Matrix NNReal

namespace StackedSVD

/-- Helper for `reTableLaw_eq_pi`: the same map as
`MeasurableEquiv.arrowProdEquivProdArrow`, with the second factor returned at the `Matrix`
type instead of the raw function type. A named `def` (not an inline `Prod.mk` ascription)
keeps the `Matrix` measurable-space instance attached from the start, so that later lemma
applications (`Measure.map_map`) never need to unify it against the equivalence's own (defeq
but syntactically different) instance for `Fin n → (Fin d → ℝ)`. -/
def tableRowSplit (n d : ℕ) (x : Fin n → ℝ × (Fin d → ℝ)) :
    (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ :=
  (fun k => (x k).1, fun k => (x k).2)

theorem measurable_tableRowSplit (n d : ℕ) : Measurable (tableRowSplit n d) := by
  refine Measurable.prodMk ?_ ?_
  · exact measurable_pi_lambda _ fun k => measurable_fst.comp (measurable_pi_apply k)
  · exact measurable_pi_lambda _ fun k => measurable_pi_lambda _ fun l =>
      (measurable_pi_apply l).comp (measurable_snd.comp (measurable_pi_apply k))

/-- U5: the rows of a table are independent with law `reRowLaw`.

The equivalence `MeasurableEquiv.arrowProdEquivProdArrow` has codomain
`(Fin n → ℝ) × (Fin n → (Fin d → ℝ))`, using Mathlib's raw `Pi` measurable-space instance;
`reTableLaw` uses `Matrix (Fin n) (Fin d) ℝ`, whose own instance is a separate (defeq)
declaration (`instMeasurableSpaceMatrix`). Direct term ascription bridges the two (checked at
full transparency), but feeding both instances into one lemma call such as `Measure.map_map`
does not (that unification runs at a stricter transparency). So `tableRowSplit` reproves
measurability of the paired map directly at the `Matrix` type, matching `measurable_reTableMap`,
and `hbridge` (an `rfl`) swaps the equivalence for `tableRowSplit` before `Measure.map_map` is
used. -/
theorem reTableLaw_eq_pi (n d : ℕ) (θ : ℝ) (v : Fin d → ℝ) :
    reTableLaw n d θ v = Measure.pi fun _ : Fin n => reRowLaw n d θ v := by
  have hmp := measurePreserving_arrowProdEquivProdArrow ℝ (Fin d → ℝ) (Fin n)
    (fun _ : Fin n => gaussianReal 0 (n : ℝ≥0)⁻¹)
    (fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1)
  have hfeq : (fun x : Fin n → ℝ × (Fin d → ℝ) =>
        (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) (tableRowSplit n d x))
      = fun x : Fin n → ℝ × (Fin d → ℝ) => fun i : Fin n =>
          (θ * (x i).1) • v + (Real.sqrt d)⁻¹ • (x i).2 := by
    funext x
    funext i
    funext j
    change θ * ((x i).1 * v j) + (Real.sqrt d)⁻¹ * (x i).2 j
      = θ * (x i).1 * v j + (Real.sqrt d)⁻¹ * (x i).2 j
    ring
  have hpi := Measure.pi_map_pi
    (μ := fun _ : Fin n =>
      (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1))
    (f := fun _ : Fin n =>
      fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2)
    (fun _ => (measurable_reRowMap d θ v).aemeasurable)
  have step1 : reTableLaw n d θ v
      = ((Measure.pi fun _ : Fin n => gaussianReal 0 (n : ℝ≥0)⁻¹).prod
          (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) := rfl
  have step2 : ((Measure.pi fun _ : Fin n => gaussianReal 0 (n : ℝ≥0)⁻¹).prod
        (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
        (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
          θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2)
      = ((Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (MeasurableEquiv.arrowProdEquivProdArrow ℝ (Fin d → ℝ) (Fin n))).map
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) := by rw [← hmp.map_eq]
  have step3 : ((Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (MeasurableEquiv.arrowProdEquivProdArrow ℝ (Fin d → ℝ) (Fin n))).map
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2)
      = ((Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (tableRowSplit n d)).map
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) := rfl
  have step4 : ((Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (tableRowSplit n d)).map
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
            θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2)
      = (Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (fun x : Fin n → ℝ × (Fin d → ℝ) =>
            (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
                θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) (tableRowSplit n d x)) :=
    Measure.map_map (measurable_reTableMap n d θ v) (measurable_tableRowSplit n d)
  have step5 : (Measure.pi fun _ : Fin n =>
        (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
        (fun x : Fin n → ℝ × (Fin d → ℝ) =>
          (fun p : (Fin n → ℝ) × Matrix (Fin n) (Fin d) ℝ =>
              θ • vecMulVec p.1 v + (Real.sqrt d)⁻¹ • p.2) (tableRowSplit n d x))
      = (Measure.pi fun _ : Fin n =>
          (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (fun x : Fin n → ℝ × (Fin d → ℝ) => fun i : Fin n =>
            (θ * (x i).1) • v + (Real.sqrt d)⁻¹ • (x i).2) :=
    congrArg (fun f : (Fin n → ℝ × (Fin d → ℝ)) → Fin n → Fin d → ℝ =>
      (Measure.pi fun _ : Fin n =>
        (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map f)
      hfeq
  have step6 : (Measure.pi fun _ : Fin n =>
        (gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
        (fun x : Fin n → ℝ × (Fin d → ℝ) => fun i : Fin n =>
          (θ * (x i).1) • v + (Real.sqrt d)⁻¹ • (x i).2)
      = Measure.pi (fun _ : Fin n =>
          ((gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
            (fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2)) := hpi
  have step7 : Measure.pi (fun _ : Fin n =>
        ((gaussianReal 0 (n : ℝ≥0)⁻¹).prod (Measure.pi fun _ : Fin d => gaussianReal 0 1)).map
          (fun p : ℝ × (Fin d → ℝ) => (θ * p.1) • v + (Real.sqrt d)⁻¹ • p.2))
      = Measure.pi fun _ : Fin n => reRowLaw n d θ v := rfl
  exact step1.trans (step2.trans (step3.trans (step4.trans (step5.trans (step6.trans step7)))))

/-- A table has Lebesgue density `∏ₖ gaussDensity d Σ (X k)`. -/
theorem reTableLaw_eq_withDensity {n d : ℕ} (hn : 0 < n) (hd : 0 < d) (θ : ℝ)
    (v : Fin d → ℝ) :
    reTableLaw n d θ v
      = (lebesgueMatrix n d).withDensity
          (fun X => ENNReal.ofReal
            (∏ k, gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k))) := by
  have hprob : ∀ _ : Fin n, IsProbabilityMeasure
      ((volume : Measure (Fin d → ℝ)).withDensity
        (fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x))) := by
    intro _
    rw [← reRowLaw_eq_withDensity hn hd θ v]
    infer_instance
  have hmeas : ∀ _ : Fin n, Measurable
      (fun x : Fin d → ℝ =>
        ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x)) :=
    fun _ => (measurable_gaussDensity d _).ennreal_ofReal
  have key := pi_withDensity (fun _ : Fin n => (volume : Measure (Fin d → ℝ)))
    (fun _ : Fin n =>
      fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x))
    hmeas
  have hfun : (fun _ : Fin n => reRowLaw n d θ v)
      = fun _ : Fin n => (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x)) :=
    funext fun _ => reRowLaw_eq_withDensity hn hd θ v
  have hden : (fun X : Matrix (Fin n) (Fin d) ℝ => ∏ k,
        ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k)))
      = fun X => ENNReal.ofReal
          (∏ k, gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k)) := by
    funext X
    rw [ENNReal.ofReal_prod_of_nonneg (fun k _ =>
      (gaussDensity_pos (det_mleCov_pos hd (by positivity) v) (X k)).le)]
  have step1 : reTableLaw n d θ v
      = Measure.pi fun _ : Fin n => reRowLaw n d θ v := reTableLaw_eq_pi n d θ v
  have step2 : (Measure.pi fun _ : Fin n => reRowLaw n d θ v)
      = Measure.pi (fun _ : Fin n => (volume : Measure (Fin d → ℝ)).withDensity
          (fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x))) :=
    congrArg (fun μ : Fin n → Measure (Fin d → ℝ) => Measure.pi μ) hfun
  have step3 : Measure.pi (fun _ : Fin n => (volume : Measure (Fin d → ℝ)).withDensity
        (fun x => ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) x)))
      = (Measure.pi fun _ : Fin n => (volume : Measure (Fin d → ℝ))).withDensity
          (fun X => ∏ i : Fin n,
            ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X i))) := key
  have step4 : (Measure.pi fun _ : Fin n => (volume : Measure (Fin d → ℝ))).withDensity
        (fun X => ∏ i : Fin n,
          ENNReal.ofReal (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X i)))
      = (lebesgueMatrix n d).withDensity
          (fun X => ∏ k, ENNReal.ofReal
            (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k))) := rfl
  have step5 : (lebesgueMatrix n d).withDensity
        (fun X => ∏ k, ENNReal.ofReal
          (gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k)))
      = (lebesgueMatrix n d).withDensity
          (fun X => ENNReal.ofReal
            (∏ k, gaussDensity d (mleCov d (θ ^ 2 / ((n : ℝ) / d)) v) (X k))) :=
    congrArg (lebesgueMatrix n d).withDensity hden
  exact step1.trans (step2.trans (step3.trans (step4.trans step5)))

/-- The joint law of the `M` tables has Lebesgue density `reDensity`. -/
theorem reJointLaw_eq_withDensity {M : ℕ} {n : Fin M → ℕ} (hn : ∀ i, 0 < n i) {d : ℕ}
    (hd : 0 < d) (θ : Fin M → ℝ) (v : Fin d → ℝ) :
    reJointLaw n d θ v
      = (Measure.pi fun i => lebesgueMatrix (n i) d).withDensity
          (fun X => ENNReal.ofReal (reDensity n d θ v X)) := by
  have hprob : ∀ i : Fin M, IsProbabilityMeasure
      ((lebesgueMatrix (n i) d).withDensity (fun X => ENNReal.ofReal
        (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k)))) := by
    intro i
    rw [← reTableLaw_eq_withDensity (hn i) hd (θ i) v]
    infer_instance
  have hmeas : ∀ i : Fin M, Measurable (fun X : Matrix (Fin (n i)) (Fin d) ℝ =>
      ENNReal.ofReal (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k))) :=
    fun i => Measurable.ennreal_ofReal (Finset.measurable_prod _ fun k _ =>
      (measurable_gaussDensity d _).comp (measurable_pi_apply k))
  have hsf : ∀ i : Fin M, SigmaFinite (lebesgueMatrix (n i) d) := fun i => by
    unfold lebesgueMatrix
    have h1 : SigmaFinite (Measure.pi fun _ : Fin d => (volume : Measure ℝ)) :=
      Measure.pi.sigmaFinite (fun _ : Fin d => (volume : Measure ℝ))
    exact Measure.pi.sigmaFinite
      (fun _ : Fin (n i) => Measure.pi fun _ : Fin d => (volume : Measure ℝ))
  have key := pi_withDensity (fun i : Fin M => lebesgueMatrix (n i) d)
    (fun i : Fin M => fun X => ENNReal.ofReal
      (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k)))
    hmeas
  have hfun : (fun i : Fin M => reTableLaw (n i) d (θ i) v)
      = fun i : Fin M => (lebesgueMatrix (n i) d).withDensity (fun X => ENNReal.ofReal
          (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k))) :=
    funext fun i => reTableLaw_eq_withDensity (hn i) hd (θ i) v
  have hden : (fun X : (i : Fin M) → Matrix (Fin (n i)) (Fin d) ℝ =>
        ∏ i : Fin M, ENNReal.ofReal
          (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k)))
      = fun X => ENNReal.ofReal (reDensity n d θ v X) := by
    funext X
    change (∏ i : Fin M, ENNReal.ofReal
        (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k)))
      = ENNReal.ofReal (∏ i : Fin M, ∏ k,
          gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k))
    rw [ENNReal.ofReal_prod_of_nonneg (fun i _ => Finset.prod_nonneg (fun k _ =>
      (gaussDensity_pos (det_mleCov_pos hd (by positivity) v) (X i k)).le))]
  have step1 : reJointLaw n d θ v
      = Measure.pi (fun i => reTableLaw (n i) d (θ i) v) := rfl
  have step2 : Measure.pi (fun i => reTableLaw (n i) d (θ i) v)
      = Measure.pi (fun i => (lebesgueMatrix (n i) d).withDensity (fun X => ENNReal.ofReal
          (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k)))) :=
    congrArg (fun μ : (i : Fin M) → Measure (Matrix (Fin (n i)) (Fin d) ℝ) =>
      Measure.pi μ) hfun
  have step3 : Measure.pi (fun i => (lebesgueMatrix (n i) d).withDensity (fun X =>
        ENNReal.ofReal (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X k))))
      = (Measure.pi fun i => lebesgueMatrix (n i) d).withDensity
          (fun X => ∏ i : Fin M, ENNReal.ofReal
            (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k))) := key
  have step4 : (Measure.pi fun i => lebesgueMatrix (n i) d).withDensity
        (fun X => ∏ i : Fin M, ENNReal.ofReal
          (∏ k, gaussDensity d (mleCov d (θ i ^ 2 / ((n i : ℝ) / d)) v) (X i k)))
      = (Measure.pi fun i => lebesgueMatrix (n i) d).withDensity
          (fun X => ENNReal.ofReal (reDensity n d θ v X)) :=
    congrArg (Measure.pi fun i => lebesgueMatrix (n i) d).withDensity hden
  exact step1.trans (step2.trans (step3.trans step4))

end StackedSVD
