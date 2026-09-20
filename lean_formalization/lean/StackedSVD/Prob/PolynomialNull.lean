/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# Zero sets of polynomials are null, and the Gaussian matrix law is absolutely continuous

Infrastructure for item S of `notes/archive/rmt_roadmap.md` (the top eigenvalue of `Xᵀ X` is simple
almost surely). Mathlib v4.33.0 has none of the four results below.

1. `StackedSVD.volume_setOf_eval_eq_zero`: the zero set of a nonzero real multivariate
   polynomial is Lebesgue null. `StackedSVD.ae_eval_ne_zero` is the `∀ᵐ` form.
2. `StackedSVD.absolutelyContinuous_pi`: `Measure.pi` respects absolute continuity.
   Mathlib has only the binary `Measure.AbsolutelyContinuous.prod`.
3. `StackedSVD.gaussianMatrix_absolutelyContinuous`: `gaussianMatrix n d ≪ volume`.
4. `StackedSVD.ae_eval_ne_zero_gaussianMatrix`: a nonzero polynomial in the entries of the
   matrix is nonzero at almost every Gaussian matrix.

The proof of 1 is an induction on the number of variables. It peels the first variable off the
polynomial with `MvPolynomial.finSuccEquiv` and off the measure with
`MeasureTheory.volume_preserving_piFinSuccAbove`, then applies Fubini
(`MeasureTheory.Measure.prod_apply_symm`). A nonzero one variable real polynomial has a finite
root set, and a finite set is null.
-/

open MeasureTheory ProbabilityTheory Set

namespace StackedSVD

/-! ### Evaluation of a polynomial is continuous -/

/-- The evaluation map `x ↦ p(x)` of a multivariate polynomial is continuous. Mathlib has
`MvPolynomial.continuous_eval` for no such signature; this is the finite variable case that
`eval_eq'` gives at once. -/
theorem continuous_mvPolynomial_eval {σ : Type*} [Finite σ] (p : MvPolynomial σ ℝ) :
    Continuous fun x : σ → ℝ => MvPolynomial.eval x p := by
  have := Fintype.ofFinite σ
  simp only [MvPolynomial.eval_eq']
  refine continuous_finsetSum _ fun m _ => continuous_const.mul ?_
  exact continuous_finsetProd _ fun i _ => (continuous_apply i).pow _

/-- The zero set of a multivariate polynomial is measurable. -/
theorem measurableSet_setOf_mvPolynomial_eval_eq_zero {σ : Type*} [Finite σ]
    (p : MvPolynomial σ ℝ) :
    MeasurableSet {x : σ → ℝ | MvPolynomial.eval x p = 0} :=
  (continuous_mvPolynomial_eval p).measurable (measurableSet_singleton 0)

/-! ### The zero set of a nonzero polynomial is null -/

/-- The `Fin k` case of `volume_setOf_eval_eq_zero`, proved by induction on `k`. -/
private theorem volume_setOf_eval_eq_zero_fin :
    ∀ (k : ℕ) (p : MvPolynomial (Fin k) ℝ), p ≠ 0 →
      (volume : Measure (Fin k → ℝ)) {x | MvPolynomial.eval x p = 0} = 0 := by
  intro k
  induction k with
  | zero =>
      intro p hp
      have hempty : {x : Fin 0 → ℝ | MvPolynomial.eval x p = 0} = ∅ := by
        ext x
        simp only [mem_ofPred_eq, mem_empty_iff_false, iff_false]
        intro hx
        refine hp (MvPolynomial.funext fun y => ?_)
        rw [Subsingleton.elim y x, map_zero]
        exact hx
      rw [hempty, measure_empty]
  | succ k ih =>
      intro p hp
      set S : Set (Fin (k + 1) → ℝ) := {x | MvPolynomial.eval x p = 0} with hS
      have hSmeas : MeasurableSet S := measurableSet_setOf_mvPolynomial_eval_eq_zero p
      -- the polynomial as a one variable polynomial over the remaining variables
      set P : Polynomial (MvPolynomial (Fin k) ℝ) := MvPolynomial.finSuccEquiv ℝ k p with hPdef
      have hPne : P ≠ 0 := by
        intro hzero
        apply hp
        have := congrArg (MvPolynomial.finSuccEquiv ℝ k).symm hzero
        simpa [hPdef] using this
      set c : MvPolynomial (Fin k) ℝ := P.coeff P.natDegree with hcdef
      have hcne : c ≠ 0 := Polynomial.leadingCoeff_ne_zero.mpr hPne
      -- peel the first variable off the measure
      set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (k + 1) => ℝ) 0 with hedef
      have hmp : MeasurePreserving e (volume : Measure (Fin (k + 1) → ℝ)) volume :=
        volume_preserving_piFinSuccAbove (fun _ : Fin (k + 1) => ℝ) 0
      have hmp' : MeasurePreserving e.symm (volume : Measure (ℝ × (Fin k → ℝ)))
          (volume : Measure (Fin (k + 1) → ℝ)) := hmp.symm e
      have hkey : (volume : Measure (ℝ × (Fin k → ℝ))) (⇑e.symm ⁻¹' S) = volume S :=
        hmp'.measure_preimage hSmeas.nullMeasurableSet
      have hTmeas : MeasurableSet (⇑e.symm ⁻¹' S) := e.symm.measurable hSmeas
      have hcons : ∀ q : ℝ × (Fin k → ℝ), e.symm q = Fin.cons q.1 q.2 := by
        intro q
        simp only [hedef, MeasurableEquiv.piFinSuccAbove, MeasurableEquiv.symm,
          MeasurableEquiv.coe_mk, Equiv.symm_symm, Fin.insertNthEquiv_zero]
        rfl
      have hset : (⇑e.symm ⁻¹' S) =
          {q : ℝ × (Fin k → ℝ) | MvPolynomial.eval (Fin.cons q.1 q.2) p = 0} := by
        ext q
        simp only [mem_preimage, hS, mem_ofPred_eq, hcons q]
      rw [← hkey, hset, Measure.volume_eq_prod]
      rw [hset] at hTmeas
      rw [Measure.prod_apply_symm hTmeas]
      -- almost every slice is finite, hence null
      have hae : ∀ᵐ z ∂(volume : Measure (Fin k → ℝ)), MvPolynomial.eval z c ≠ 0 := by
        rw [ae_iff]
        simpa using ih c hcne
      have hzero : (fun z : Fin k → ℝ => (volume : Measure ℝ)
          ((fun y => (y, z)) ⁻¹'
            {q : ℝ × (Fin k → ℝ) | MvPolynomial.eval (Fin.cons q.1 q.2) p = 0}))
          =ᵐ[(volume : Measure (Fin k → ℝ))] 0 := by
        filter_upwards [hae] with z hz
        have hQ : Polynomial.map (MvPolynomial.eval z) P ≠ 0 := by
          intro hQ0
          apply hz
          have : (Polynomial.map (MvPolynomial.eval z) P).coeff P.natDegree
              = MvPolynomial.eval z c := by
            rw [Polynomial.coeff_map, hcdef]
          rw [hQ0] at this
          simpa using this.symm
        have hslice : ((fun y => (y, z)) ⁻¹'
            {q : ℝ × (Fin k → ℝ) | MvPolynomial.eval (Fin.cons q.1 q.2) p = 0})
            = {y : ℝ | (Polynomial.map (MvPolynomial.eval z) P).IsRoot y} := by
          ext y
          simp only [mem_preimage, mem_ofPred_eq, Polynomial.IsRoot.def]
          rw [MvPolynomial.eval_eq_eval_mv_eval' z y p, hPdef]
        simp only [Pi.zero_apply, hslice]
        exact (Polynomial.finite_setOfPred_isRoot hQ).measure_zero _
      rw [lintegral_congr_ae hzero]
      simp

/-- **The zero set of a nonzero real polynomial is Lebesgue null.** Not in Mathlib v4.33.0. -/
theorem volume_setOf_eval_eq_zero {σ : Type*} [Fintype σ] (p : MvPolynomial σ ℝ) (hp : p ≠ 0) :
    (volume : Measure (σ → ℝ)) {x | MvPolynomial.eval x p = 0} = 0 := by
  classical
  set g : Fin (Fintype.card σ) ≃ σ := (Fintype.equivFin σ).symm with hg
  set E := MeasurableEquiv.piCongrLeft (fun _ : σ => ℝ) g with hE
  have hmp : MeasurePreserving E (volume : Measure (Fin (Fintype.card σ) → ℝ))
      (volume : Measure (σ → ℝ)) :=
    volume_measurePreserving_piCongrLeft (fun _ : σ => ℝ) g
  have hF : MeasurePreserving (⇑E.symm) (volume : Measure (σ → ℝ))
      (volume : Measure (Fin (Fintype.card σ) → ℝ)) := hmp.symm E
  have hFapp : ∀ y : σ → ℝ, ⇑E.symm y = y ∘ ⇑g := by
    intro y
    funext a
    simp [hE, MeasurableEquiv.piCongrLeft]
  have hcomp : ∀ y : σ → ℝ, ((y ∘ ⇑g) ∘ ⇑g.symm) = y := by
    intro y
    funext i
    simp
  have hren : MvPolynomial.rename (⇑g.symm) p ≠ 0 := fun h =>
    hp (MvPolynomial.rename_injective _ g.symm.injective (by simpa using h))
  have hpre : ⇑E.symm ⁻¹'
      {x : Fin (Fintype.card σ) → ℝ | MvPolynomial.eval x (MvPolynomial.rename (⇑g.symm) p) = 0}
      = {y : σ → ℝ | MvPolynomial.eval y p = 0} := by
    ext y
    rw [Set.mem_preimage, Set.mem_ofPred_eq, Set.mem_ofPred_eq, hFapp y,
      MvPolynomial.eval_rename, hcomp y]
  have hmeas : MeasurableSet
      {x : Fin (Fintype.card σ) → ℝ | MvPolynomial.eval x (MvPolynomial.rename (⇑g.symm) p) = 0} :=
    measurableSet_setOf_mvPolynomial_eval_eq_zero _
  have hval := hF.measure_preimage hmeas.nullMeasurableSet
  rw [hpre] at hval
  rw [hval]
  exact volume_setOf_eval_eq_zero_fin _ _ hren

/-- **A nonzero real polynomial is nonzero at almost every point.** -/
theorem ae_eval_ne_zero {σ : Type*} [Fintype σ] (p : MvPolynomial σ ℝ) (hp : p ≠ 0) :
    ∀ᵐ x ∂(volume : Measure (σ → ℝ)), MvPolynomial.eval x p ≠ 0 := by
  rw [ae_iff]
  simpa using volume_setOf_eval_eq_zero p hp

/-! ### Absolute continuity of product measures -/

/-- Absolute continuity passes through `Measure.map`. Mathlib v4.33.0 has no such lemma. -/
theorem absolutelyContinuous_map {α β : Type*} [MeasurableSpace α] [MeasurableSpace β]
    {μ ν : Measure α} {f : α → β} (h : μ ≪ ν) (hf : Measurable f) :
    μ.map f ≪ ν.map f := by
  refine Measure.AbsolutelyContinuous.mk fun s hs hνs => ?_
  rw [Measure.map_apply hf hs] at hνs ⊢
  exact h hνs

/-- The `Fin k` case of `absolutelyContinuous_pi`, proved by induction on `k`. -/
private theorem absolutelyContinuous_pi_fin :
    ∀ (k : ℕ) {α : Fin k → Type*} [_inst : ∀ i, MeasurableSpace (α i)]
      (μ ν : ∀ i, Measure (α i)) [∀ i, SigmaFinite (μ i)] [∀ i, SigmaFinite (ν i)],
      (∀ i, μ i ≪ ν i) → Measure.pi μ ≪ Measure.pi ν := by
  intro k
  induction k with
  | zero =>
      intro α _ μ ν _ _ _
      rw [Measure.pi_of_empty μ, Measure.pi_of_empty ν]
  | succ k ih =>
      intro α _ μ ν _ _ h
      set e := MeasurableEquiv.piFinSuccAbove α 0 with hedef
      have hμ := measurePreserving_piFinSuccAbove μ 0
      have hν := measurePreserving_piFinSuccAbove ν 0
      have hprod : (μ 0).prod (Measure.pi fun j => μ ((0 : Fin (k + 1)).succAbove j))
          ≪ (ν 0).prod (Measure.pi fun j => ν ((0 : Fin (k + 1)).succAbove j)) :=
        (h 0).prod (ih _ _ fun j => h _)
      have hmap := absolutelyContinuous_map hprod e.symm.measurable
      rwa [← hμ.map_eq, ← hν.map_eq, MeasurableEquiv.map_symm_map,
        MeasurableEquiv.map_symm_map] at hmap

/-- **`Measure.pi` respects absolute continuity.** Mathlib v4.33.0 has only the binary
`MeasureTheory.Measure.AbsolutelyContinuous.prod`. -/
theorem absolutelyContinuous_pi {ι : Type*} [Fintype ι] {α : ι → Type*}
    [∀ i, MeasurableSpace (α i)] (μ ν : ∀ i, Measure (α i))
    [∀ i, SigmaFinite (μ i)] [∀ i, SigmaFinite (ν i)] (h : ∀ i, μ i ≪ ν i) :
    Measure.pi μ ≪ Measure.pi ν := by
  classical
  set g : Fin (Fintype.card ι) ≃ ι := (Fintype.equivFin ι).symm with hg
  have hμ := measurePreserving_piCongrLeft μ g
  have hν := measurePreserving_piCongrLeft ν g
  have hfin : Measure.pi (fun i => μ (g i)) ≪ Measure.pi (fun i => ν (g i)) :=
    absolutelyContinuous_pi_fin _ _ _ fun i => h (g i)
  have hmap := absolutelyContinuous_map hfin
    (MeasurableEquiv.piCongrLeft (fun i => α i) g).measurable
  rwa [hμ.map_eq, hν.map_eq] at hmap

/-! ### The Gaussian matrix law is absolutely continuous -/

/-- Lebesgue measure on matrices, the product measure over the entries. `Matrix` is a `def`, so
Mathlib v4.33.0 gives it no `MeasureSpace` instance, just as it gives it no `MeasurableSpace`
instance (see `StackedSVD.instMeasurableSpaceMatrix` in `Defs.lean`). `inferInstanceAs` keeps
this definitionally equal to the pi instance, and `MeasureTheory.volume_pi` is `rfl`, so the
sigma-algebra of this instance is the one of `instMeasurableSpaceMatrix`. -/
noncomputable instance instMeasureSpaceMatrix {m n α : Type*} [Fintype m] [Fintype n]
    [MeasureSpace α] : MeasureSpace (Matrix m n α) :=
  inferInstanceAs (MeasureSpace (m → n → α))

/-- **The canonical Gaussian matrix law is absolutely continuous with respect to Lebesgue
measure.** -/
theorem gaussianMatrix_absolutelyContinuous (n d : ℕ) :
    gaussianMatrix n d ≪ (volume : Measure (Matrix (Fin n) (Fin d) ℝ)) := by
  have hinner : ∀ _ : Fin n, (Measure.pi fun _ : Fin d => gaussianReal 0 1)
      ≪ (volume : Measure (Fin d → ℝ)) := fun _ =>
    absolutelyContinuous_pi _ _ fun _ => gaussianReal_absolutelyContinuous 0 one_ne_zero
  exact absolutelyContinuous_pi _ _ hinner

/-! ### A polynomial in the matrix entries is nonzero at almost every Gaussian matrix -/

/-- Reading a matrix as a function on the pairs of indices preserves the product measure. -/
theorem measurePreserving_uncurry {ι κ : Type*} [Fintype ι] [Fintype κ] {α : Type*}
    [MeasurableSpace α] (μ : Measure α) [SigmaFinite μ] :
    MeasurePreserving (fun (X : ι → κ → α) (q : ι × κ) => X q.1 q.2)
      (Measure.pi fun _ : ι => Measure.pi fun _ : κ => μ)
      (Measure.pi fun _ : ι × κ => μ) := by
  have hmeas : Measurable fun (X : ι → κ → α) (q : ι × κ) => X q.1 q.2 :=
    measurable_pi_lambda _ fun q => (measurable_pi_apply q.2).comp (measurable_pi_apply q.1)
  refine ⟨hmeas, ?_⟩
  refine (Measure.pi_eq fun s hs => ?_).symm
  rw [Measure.map_apply hmeas (MeasurableSet.univ_pi hs)]
  have hpre : (fun (X : ι → κ → α) (q : ι × κ) => X q.1 q.2) ⁻¹' (Set.univ.pi s)
      = Set.univ.pi fun i => Set.univ.pi fun j => s (i, j) := by
    ext X
    simp [Prod.forall]
  rw [hpre, Measure.pi_pi]
  simp only [Measure.pi_pi]
  rw [Fintype.prod_prod_type]

/-- **A nonzero polynomial in the matrix entries is nonzero at almost every Gaussian matrix.**
The polynomial is indexed by the pairs `(row, column)`; the evaluation point of the entry
`(i, j)` is `X i j`. -/
theorem ae_eval_ne_zero_gaussianMatrix {n d : ℕ} (p : MvPolynomial (Fin n × Fin d) ℝ)
    (hp : p ≠ 0) :
    ∀ᵐ X ∂(gaussianMatrix n d), MvPolynomial.eval (fun ij => X ij.1 ij.2) p ≠ 0 := by
  have hmp : MeasurePreserving (fun (X : Fin n → Fin d → ℝ) (q : Fin n × Fin d) => X q.1 q.2)
      (volume : Measure (Fin n → Fin d → ℝ))
      (volume : Measure (Fin n × Fin d → ℝ)) :=
    measurePreserving_uncurry (volume : Measure ℝ)
  have hkey := hmp.measure_preimage
    (measurableSet_setOf_mvPolynomial_eval_eq_zero p).nullMeasurableSet
  rw [volume_setOf_eval_eq_zero p hp] at hkey
  have hvol : ∀ᵐ X ∂(volume : Measure (Matrix (Fin n) (Fin d) ℝ)),
      MvPolynomial.eval (fun ij => X ij.1 ij.2) p ≠ 0 := by
    rw [ae_iff]
    simp only [ne_eq, not_not]
    exact hkey
  exact hvol.filter_mono (gaussianMatrix_absolutelyContinuous n d).ae_le

end StackedSVD
