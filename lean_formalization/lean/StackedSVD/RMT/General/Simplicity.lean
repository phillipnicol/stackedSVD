/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Simplicity
import StackedSVD.Prob.NoiseLaw
import StackedSVD.Prob.PolynomialNull

/-!
# Item S at a general noise law: the top eigenvalue of `Xᵀ X` is simple almost surely

Unit G8 of `notes/archive/prop_single_table_general.md` section 5, Stage 1 of the non-Gaussian
extension. This file is the general-law twin of `RMT/Simplicity.lean`. Every proof is the
same argument. `gaussianMatrix n d` becomes `noiseMatrix ν n d`. The absolute continuity of
the standard Gaussian law becomes the hypothesis `hac : ν ≪ volume`. `RMT/Simplicity.lean`,
`Prob/PolynomialNull.lean` and `Prob/NoiseLaw.lean` supply every deterministic and general
lemma; this file imports them and does not edit them.

1. `noiseMatrix_absolutelyContinuous`. `noiseMatrix ν n d ≪ volume`, by
   `absolutelyContinuous_pi` applied twice.
2. `ae_eval_ne_zero_noiseMatrix`. A nonzero polynomial in the matrix entries is nonzero at
   almost every `noiseMatrix ν n d` matrix.
3. `topSimple_ae_affine_general`. Item S, canonical form, for the affine family `A + t • Z`.
4. `topSimple_ae_noiseMatrix`. Item S with no shift.
5. `SpikedModel.singleTableLaw_topSimple_of_general`. The `topSimple` field of
   `SingleTableLaw`, at a general law.

Every statement here carries `[SigmaFinite ν]`, including the last one. The frozen statement
of the plan (`notes/archive/prop_single_table_general.md`, section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35) carries no such
instance. Section 8 of this unit's report
records why the instance cannot be dropped: `topSimple_ae_affine_general` needs it to invoke
`absolutelyContinuous_pi`, and neither `hac : ν ≪ volume` nor `hG : m.GeneralNoise ν` gives a
sigma-finite `ν`.
-/

open MeasureTheory ProbabilityTheory Set
open scoped Matrix

namespace StackedSVD

/-! ### Item S, general law -/

/-- **The i.i.d. matrix law is absolutely continuous** when the entry law is. General-law twin
of `gaussianMatrix_absolutelyContinuous` (`Prob/PolynomialNull.lean:232`). -/
theorem noiseMatrix_absolutelyContinuous {ν : Measure ℝ} [SigmaFinite ν]
    (hac : ν ≪ (volume : Measure ℝ)) (n d : ℕ) :
    noiseMatrix ν n d ≪ (volume : Measure (Matrix (Fin n) (Fin d) ℝ)) := by
  have hinner : ∀ _ : Fin n, (Measure.pi fun _ : Fin d => ν)
      ≪ (volume : Measure (Fin d → ℝ)) := fun _ =>
    absolutelyContinuous_pi _ _ fun _ => hac
  exact absolutelyContinuous_pi _ _ hinner

/-- **A nonzero polynomial in the entries is nonzero at almost every i.i.d. matrix.**
General-law twin of `ae_eval_ne_zero_gaussianMatrix` (`Prob/PolynomialNull.lean:264`). -/
theorem ae_eval_ne_zero_noiseMatrix {ν : Measure ℝ} [SigmaFinite ν]
    (hac : ν ≪ (volume : Measure ℝ)) {n d : ℕ} (p : MvPolynomial (Fin n × Fin d) ℝ)
    (hp : p ≠ 0) :
    ∀ᵐ X ∂(noiseMatrix ν n d), MvPolynomial.eval (fun ij => X ij.1 ij.2) p ≠ 0 := by
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
  exact hvol.filter_mono (noiseMatrix_absolutelyContinuous hac n d).ae_le

/-- **Item S, canonical form, at a general law.** Twin of `topSimple_ae_affine`
(`RMT/Simplicity.lean:340`). -/
theorem topSimple_ae_affine_general {ν : Measure ℝ} [SigmaFinite ν]
    (hac : ν ≪ (volume : Measure ℝ)) {n d : ℕ} (hn : 0 < n) (hd : 0 < d)
    (A : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) :
    ∀ᵐ Z ∂(noiseMatrix ν n d),
      TopSimple ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) := by
  rcases le_total d n with hdn | hnd
  · filter_upwards [ae_eval_ne_zero_noiseMatrix hac (gramPolyRight A t)
      (gramPolyRight_ne_zero hd hdn A ht)] with Z hZ
    rw [eval_gramPolyRight] at hZ
    exact topSimple_of_charRes_ne_zero _ hd hZ
  · filter_upwards [ae_eval_ne_zero_noiseMatrix hac (gramPolyLeft A t)
      (gramPolyLeft_ne_zero hn hnd A ht)] with Z hZ
    rw [eval_gramPolyLeft] at hZ
    exact topSimple_of_gram_left hn hd _ (left_ne_zero_of_mul hZ) (right_ne_zero_of_mul hZ)

/-- **Item S for `noiseMatrix` with no shift.** Twin of `topSimple_ae_gaussianMatrix`
(`RMT/Simplicity.lean:355`). -/
theorem topSimple_ae_noiseMatrix {ν : Measure ℝ} [SigmaFinite ν]
    (hac : ν ≪ (volume : Measure ℝ)) (n d : ℕ) (hn : 0 < n) (hd : 0 < d) :
    ∀ᵐ Z ∂(noiseMatrix ν n d), TopSimple (Zᵀ * Z) (isHermitian_transpose_mul_self Z) := by
  have h := topSimple_ae_affine_general hac hn hd (0 : Matrix (Fin n) (Fin d) ℝ)
    (one_ne_zero (α := ℝ))
  filter_upwards [h] with Z hZ
  simpa using hZ

private theorem measurable_entry {n d : ℕ} (i : Fin n) (j : Fin d) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => Z i j := by
  change Measurable fun Z : Fin n → Fin d → ℝ => Z i j
  exact (measurable_pi_apply j).comp (measurable_pi_apply i)

private theorem measurable_affine {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A + t • Z := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  change Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A i j + t * Z i j
  exact ((measurable_entry i j).const_mul t).const_add (A i j)

namespace SpikedModel

/-- **Item S for `SingleTableLaw` at a general law.** This is the `topSimple` field. Twin of
`singleTableLaw_topSimple_of_gaussian` (`RMT/Simplicity.lean:375`).

**Deviation from the frozen statement of the plan:** this theorem carries `[SigmaFinite ν]`;
the frozen line (`notes/archive/prop_single_table_general.md` section 2, held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35) has no such
instance binder. A direct attempt at
the frozen signature (no instance) fails at the call to `topSimple_ae_affine_general`, which
itself needs `[SigmaFinite ν]` to invoke `absolutelyContinuous_pi` on the two factors of
`noiseMatrix`. Neither `hac : ν ≪ volume` nor `hG : m.GeneralNoise ν` gives a sigma-finite
`ν`: `HasLaw.ae_iff` (`Mathlib/Probability/HasLaw.lean:89`) needs no such instance, and the
`SpikedModel` structure carries no measure hypothesis on `ν`. `[SigmaFinite ν]` is the
weaker of the two candidate fixes; the report of this unit records the exact failure. -/
theorem singleTableLaw_topSimple_of_general {Ω : ℕ → Type*}
    [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}
    (m : SpikedModel μ n d) {ν : Measure ℝ} [SigmaFinite ν]
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω)) := by
  set A : Matrix (Fin (n N)) (Fin (d N)) ℝ :=
    m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) with hA
  set t : ℝ := (Real.sqrt (d N))⁻¹ with hts
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.hd N
  have ht : t ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  set p : Matrix (Fin (n N)) (Fin (d N)) ℝ → Prop := fun Z =>
    TopSimple ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) with hp
  have hpm : Measurable p := by
    rw [← measurableSet_setOfPred]
    have : {Z | p Z} = (fun Z => A + t • Z) ⁻¹'
        {Y : Matrix (Fin (n N)) (Fin (d N)) ℝ |
          TopSimple (Yᵀ * Y) (isHermitian_transpose_mul_self Y)} := rfl
    rw [this]
    exact measurable_affine A t measurableSet_topSimple
  have hae : ∀ᵐ ω ∂(μ N), p (m.Z N ω) :=
    ((hG N).ae_iff hpm).mpr (topSimple_ae_affine_general hac (m.hn N) (m.hd N) A ht)
  filter_upwards [hae] with ω hω
  exact hω

end SpikedModel

end StackedSVD
