/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-! # The general noise law

`NoiseLaw ν` is the paper's `assum:general_noise`: a law `ν` on `ℝ` with mean 0, variance 1 and
finite fourth moment, the law of one noise entry after the `√d` scaling. `noiseMatrix ν n d` is
the law of an `n × d` matrix with i.i.d. entries of law `ν`. `SpikedModel.GeneralNoise` and
`MultiTableModel.JointGeneralNoise` generalize `Defs.lean`'s `GaussianNoise` and
`JointGaussianNoise` to an arbitrary noise law `ν`; the Gaussian predicates of `Defs.lean` are
the case `ν = gaussianReal 0 1`, by `Iff.rfl`.

The `PiLaw` section below is a copy of `R2.Centered` (`RMT/R2.lean:288` to `:420`) at a general
product law `Measure.pi fun _ : Fin d => ν` over a general coordinate type `α`, rather than only
`ℝ` under `gaussianReal 0 1`, so that it applies to rows (`Fin D → ℝ`) as well as to entries.
`RMT/R2.lean` is left unchanged so that the Gaussian tree does not rebuild. -/

open MeasureTheory ProbabilityTheory
open scoped Matrix

namespace StackedSVD

/-! ### Part 1: the general noise law -/

/-- A noise law: mean 0, variance 1, finite fourth moment. The paper's `assum:general_noise`
for one entry after the `√d` scaling, with one law for every `N`, so the paper's constant `C`
is `∫ x⁴ ∂ν`. The field `mom4` makes the two identities honest: without it `∫ x ∂ν = 0` and
`∫ x² ∂ν = 1` would hold vacuously for a law without moments (Lean's convention `∫ f = 0`
for a non-integrable `f`). -/
structure NoiseLaw (ν : Measure ℝ) : Prop where
  prob : IsProbabilityMeasure ν
  mean : ∫ x, x ∂ν = 0
  var : ∫ x, x ^ 2 ∂ν = 1
  mom4 : Integrable (fun x => x ^ 4) ν

/-- The law of an `n × d` matrix with i.i.d. entries of law `ν`.
`gaussianMatrix n d = noiseMatrix (gaussianReal 0 1) n d` holds by `rfl`. -/
noncomputable def noiseMatrix (ν : Measure ℝ) (n d : ℕ) : Measure (Matrix (Fin n) (Fin d) ℝ) :=
  Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => ν

/-- `noiseMatrix` is a probability measure when `ν` is. `Matrix` is not reducible, so the
instance is found on the raw `Measure.pi` term and transported by `exact`, as in
`instIsProbabilityMeasureGaussianMatrix` (`GaussianMatrix.lean`). -/
instance instIsProbabilityMeasureNoiseMatrix (ν : Measure ℝ) [IsProbabilityMeasure ν]
    (n d : ℕ) : IsProbabilityMeasure (noiseMatrix ν n d) := by
  have h : IsProbabilityMeasure (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => ν) :=
    inferInstance
  exact h

/-- `gaussianMatrix` is `noiseMatrix` at `ν = gaussianReal 0 1`. -/
theorem gaussianMatrix_eq_noiseMatrix (n d : ℕ) :
    gaussianMatrix n d = noiseMatrix (gaussianReal 0 1) n d := rfl

section SpikedModelNoise

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- The entries of `Z N` are i.i.d. of law `ν` at every `N`; `GaussianNoise` (`Defs.lean`)
is the case `ν = gaussianReal 0 1`, by `Iff.rfl`. The moment conditions are the separate
hypothesis `NoiseLaw ν`. -/
def SpikedModel.GeneralNoise (m : SpikedModel μ n d) (ν : Measure ℝ) : Prop :=
  ∀ N, HasLaw (m.Z N) (noiseMatrix ν (n N) (d N)) (μ N)

/-- `GaussianNoise` is `GeneralNoise` at `ν = gaussianReal 0 1`. -/
theorem SpikedModel.gaussianNoise_iff (m : SpikedModel μ n d) :
    m.GaussianNoise ↔ m.GeneralNoise (gaussianReal 0 1) := Iff.rfl

end SpikedModelNoise

section MultiTableModelNoise

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- At every `N` the `M` noise matrices are independent with i.i.d. entries of law `ν`;
`JointGaussianNoise` (`Defs.lean`) is the case `ν = gaussianReal 0 1`, by `Iff.rfl`. -/
def MultiTableModel.JointGeneralNoise (m : MultiTableModel μ M n d) (ν : Measure ℝ) : Prop :=
  ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω)
    (Measure.pi fun i : Fin M => noiseMatrix ν (n i N) (d N)) (μ N)

/-- `JointGaussianNoise` is `JointGeneralNoise` at `ν = gaussianReal 0 1`. -/
theorem MultiTableModel.jointGaussianNoise_iff (m : MultiTableModel μ M n d) :
    m.JointGaussianNoise ↔ m.JointGeneralNoise (gaussianReal 0 1) := Iff.rfl

/-- The one-table law of table `i` of a joint general law, the mirror of
`MultiTableModel.gaussianNoise_of_joint` (`SVDStack/Gram.lean`): the coordinate map of a
product law is measure preserving. -/
theorem MultiTableModel.generalNoise_of_joint (m : MultiTableModel μ M n d) {ν : Measure ℝ}
    [IsProbabilityMeasure ν] (hG : m.JointGeneralNoise ν) (i : Fin M) :
    (m.tbl i).GeneralNoise ν := by
  intro N
  have : ∀ k : Fin M, IsProbabilityMeasure (noiseMatrix ν (n k N) (d N)) :=
    fun _ => inferInstance
  exact (measurePreserving_eval (fun k : Fin M => noiseMatrix ν (n k N) (d N)) i).fun_comp_hasLaw
    (hG N)

end MultiTableModelNoise

/-! ### Elementary consequences of `NoiseLaw` -/

/-- Every power `x ^ k` with `k ≤ 4` is integrable: `|x| ^ k ≤ 1 + x ^ 4`. -/
theorem NoiseLaw.integrable_pow {ν : Measure ℝ} (hν : NoiseLaw ν) {k : ℕ} (hk : k ≤ 4) :
    Integrable (fun x => x ^ k) ν := by
  have := hν.prob
  have hg : Integrable (fun x : ℝ => 1 + x ^ 4) ν := (integrable_const 1).add hν.mom4
  refine Integrable.mono' hg (measurable_id.pow_const k).aestronglyMeasurable ?_
  filter_upwards with x
  have hx4 : (0 : ℝ) ≤ x ^ 4 := by positivity
  rw [Real.norm_eq_abs, abs_pow]
  rcases le_total |x| 1 with hx | hx
  · have h1 : |x| ^ k ≤ 1 := pow_le_one₀ (abs_nonneg x) hx
    linarith
  · have h2 : |x| ^ k ≤ |x| ^ 4 := pow_le_pow_right₀ hx hk
    have h3 : |x| ^ 4 = x ^ 4 := Even.pow_abs (by norm_num) x
    linarith

/-- `∫ x ^ 4 ∂ν` is nonnegative, for every measure `ν`. -/
theorem integral_pow_four_nonneg (ν : Measure ℝ) : 0 ≤ ∫ x, x ^ 4 ∂ν :=
  integral_nonneg fun x => by positivity

/-- `1 ≤ ∫ x ^ 4 ∂ν`, from `0 ≤ ∫ (x ^ 2 - 1) ^ 2 ∂ν` and `∫ x ^ 2 ∂ν = 1`. -/
theorem NoiseLaw.one_le_integral_pow_four {ν : Measure ℝ} (hν : NoiseLaw ν) :
    1 ≤ ∫ x, x ^ 4 ∂ν := by
  have := hν.prob
  have h2 : Integrable (fun x : ℝ => x ^ 2) ν := hν.integrable_pow (by norm_num)
  have h2c : Integrable (fun x : ℝ => 2 * x ^ 2) ν := h2.const_mul 2
  have hfg : Integrable (fun x : ℝ => x ^ 4 - 2 * x ^ 2) ν := hν.mom4.sub h2c
  have hnn : 0 ≤ ∫ x, (x ^ 2 - 1) ^ 2 ∂ν := integral_nonneg fun x => sq_nonneg _
  have heq : (fun x : ℝ => (x ^ 2 - 1) ^ 2) = fun x => x ^ 4 - 2 * x ^ 2 + 1 := by
    funext x; ring
  rw [heq, integral_add hfg (integrable_const 1), integral_sub hν.mom4 h2c,
    integral_const_mul, hν.var] at hnn
  simp at hnn
  linarith

/-- The standard Gaussian is a noise law: a probability measure with mean `0`, second
moment `1` and a finite fourth moment. The proofs repeat `R2.integral_sq_gauss` and
`R2.integrable_pow_gauss` (`RMT/R2.lean`), which this file cannot import. -/
theorem noiseLaw_gaussian : NoiseLaw (gaussianReal 0 1) where
  prob := inferInstance
  mean := by simp
  var := by
    have hm : MemLp (id : ℝ → ℝ) 2 (gaussianReal 0 1) := memLp_id_gaussianReal' 2 (by simp)
    have h := variance_eq_sub (μ := gaussianReal 0 1) hm
    rw [variance_id_gaussianReal] at h
    simp only [Pi.pow_apply, id_eq] at h
    rw [integral_id_gaussianReal] at h
    simpa using h.symm
  mom4 := by
    have hm : MemLp (id : ℝ → ℝ) ((4 : ℕ) : ENNReal) (gaussianReal 0 1) :=
      memLp_id_gaussianReal' _ (by simp)
    have h := hm.integrable_norm_pow'
    rw [← integrable_norm_iff (by fun_prop)]
    simpa [norm_pow] using h

/-! ### Part 2 and 3: coordinate lemmas at a general product law

This section is a copy of `R2.Centered` (`RMT/R2.lean:288` to `:420`) at a general coordinate
type `α` and a general product law `Measure.pi fun _ : Fin d => ν`, `ν : Measure α`, in place of
`ℝ` under `gaussianReal 0 1`. -/

namespace PiLaw

section Coord

variable {α : Type*} [MeasurableSpace α] (ν : Measure α) [IsProbabilityMeasure ν] {d : ℕ}

/-- Evaluation at coordinate `j` pushes the product law `Measure.pi fun _ => ν` forward to `ν`. -/
theorem measurePreserving_coord (j : Fin d) :
    MeasurePreserving (fun x : Fin d → α => x j) (Measure.pi fun _ : Fin d => ν) ν :=
  measurePreserving_eval _ j

/-- The integral of a function of one coordinate transfers to the factor law `ν`. -/
theorem integral_comp_coord {f : α → ℝ} (hf : AEStronglyMeasurable f ν) (j : Fin d) :
    ∫ x, f (x j) ∂(Measure.pi fun _ : Fin d => ν) = ∫ t, f t ∂ν := by
  have hmp := measurePreserving_coord ν (d := d) j
  have h := integral_map (φ := fun x : Fin d → α => x j) (μ := Measure.pi fun _ : Fin d => ν)
    (measurable_pi_apply j).aemeasurable (f := f) (by rw [hmp.map_eq]; exact hf)
  rw [hmp.map_eq] at h
  exact h.symm

/-- A function of one coordinate is integrable against the product law when `f` is against `ν`. -/
theorem integrable_comp_coord {f : α → ℝ} (hf : Integrable f ν) (j : Fin d) :
    Integrable (fun x : Fin d → α => f (x j)) (Measure.pi fun _ : Fin d => ν) :=
  memLp_one_iff_integrable.mp
    ((memLp_one_iff_integrable.mpr hf).comp_measurePreserving (measurePreserving_coord ν j))

/-- Two distinct coordinates are independent under the product law. -/
theorem indepFun_coord {a b : Fin d} (hab : a ≠ b) :
    IndepFun (fun x : Fin d → α => x a) (fun x : Fin d → α => x b)
      (Measure.pi fun _ : Fin d => ν) :=
  (iIndepFun_pi (X := fun _ : Fin d => (id : α → α)) fun _ => aemeasurable_id).indepFun hab

end Coord

/-- A centered, square-integrable function of one coordinate (copy of `R2.Centered` at a
general law `ν`). -/
structure Centered {α : Type*} [MeasurableSpace α] (ν : Measure α) (h : α → ℝ) : Prop where
  meas : Measurable h
  int : Integrable h ν
  sqInt : Integrable (fun t => h t ^ 2) ν
  zero : ∫ t, h t ∂ν = 0

section CenteredLemmas

variable {α : Type*} [MeasurableSpace α] {ν : Measure α} [IsProbabilityMeasure ν] {h : α → ℝ}
  {d : ℕ}

private theorem sq_sum_expand (r u : Fin d → ℝ) :
    (∑ a, r a * u a) ^ 2 = ∑ a, ∑ b, r a * r b * (u a * u b) := by
  rw [sq, Finset.sum_mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- The product of `h` at two coordinates is integrable. -/
theorem Centered.integrable_pair (hh : Centered ν h) (a b : Fin d) :
    Integrable (fun x : Fin d → α => h (x a) * h (x b)) (Measure.pi fun _ : Fin d => ν) := by
  rcases eq_or_ne a b with rfl | hab
  · have hI := integrable_comp_coord ν hh.sqInt a
    simpa [pow_two] using hI
  · have hind := (indepFun_coord ν hab).comp hh.meas hh.meas
    simp only [Function.comp_def] at hind
    exact hind.integrable_mul (integrable_comp_coord ν hh.int a)
      (integrable_comp_coord ν hh.int b)

/-- The second moment of `h` at two coordinates: `∫ t, h t ^ 2 ∂ν` on the diagonal, `0` off it. -/
theorem Centered.integral_pair (hh : Centered ν h) (a b : Fin d) :
    ∫ x : Fin d → α, h (x a) * h (x b) ∂(Measure.pi fun _ : Fin d => ν)
      = if a = b then ∫ t, h t ^ 2 ∂ν else 0 := by
  rcases eq_or_ne a b with rfl | hab
  · have hI := integral_comp_coord ν hh.sqInt.aestronglyMeasurable a
    simpa [pow_two] using hI
  · rw [if_neg hab]
    have hind := (indepFun_coord ν hab).comp hh.meas hh.meas
    simp only [Function.comp_def] at hind
    have hmul := hind.integral_fun_mul_eq_mul_integral
      (integrable_comp_coord ν hh.int a).aestronglyMeasurable
      (integrable_comp_coord ν hh.int b).aestronglyMeasurable
    rw [hmul, integral_comp_coord ν hh.int.aestronglyMeasurable a, hh.zero, zero_mul]

/-- A weighted sum of `h` at every coordinate, squared, is integrable. -/
theorem Centered.integrable_sq_sum (hh : Centered ν h) (r : Fin d → ℝ) :
    Integrable (fun x : Fin d → α => (∑ a, r a * h (x a)) ^ 2)
      (Measure.pi fun _ : Fin d => ν) := by
  simp only [sq_sum_expand]
  refine integrable_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ => ?_
  exact ((hh.integrable_pair a b).const_mul (r a * r b))

/-- The second moment of a weighted sum of `h` at every coordinate. -/
theorem Centered.integral_sq_sum (hh : Centered ν h) (r : Fin d → ℝ) :
    ∫ x : Fin d → α, (∑ a, r a * h (x a)) ^ 2 ∂(Measure.pi fun _ : Fin d => ν)
      = (∫ t, h t ^ 2 ∂ν) * ∑ a, r a ^ 2 := by
  simp only [sq_sum_expand]
  rw [integral_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ =>
    ((hh.integrable_pair a b).const_mul (r a * r b))]
  have hb : ∀ a : Fin d,
      ∫ x : Fin d → α, ∑ b, r a * r b * (h (x a) * h (x b)) ∂(Measure.pi fun _ : Fin d => ν)
        = (∫ t, h t ^ 2 ∂ν) * r a ^ 2 := by
    intro a
    rw [integral_finsetSum _ fun b _ => ((hh.integrable_pair a b).const_mul (r a * r b))]
    have : ∀ b : Fin d,
        ∫ x : Fin d → α, r a * r b * (h (x a) * h (x b)) ∂(Measure.pi fun _ : Fin d => ν)
          = r a * r b * (if a = b then ∫ t, h t ^ 2 ∂ν else 0) := by
      intro b
      rw [integral_const_mul, hh.integral_pair a b]
    simp only [this, mul_ite, mul_zero]
    rw [Finset.sum_ite_eq Finset.univ a fun b => r a * r b * ∫ t, h t ^ 2 ∂ν]
    simp only [Finset.mem_univ, if_true]
    ring
  rw [Finset.sum_congr rfl fun a _ => hb a, ← Finset.mul_sum]

end CenteredLemmas

end PiLaw

end StackedSVD
