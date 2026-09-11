/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.SimplicityR
import StackedSVD.RMT.Symmetry
import Mathlib.LinearAlgebra.Matrix.Charpoly.Basic
import Mathlib.MeasureTheory.Measure.Prod

/-!
# Task U5(b): delocalization of the top-`r` eigen-directions of a Gaussian Gram

Deliverable (b) of `notes/archive/rankr_plan_A.md` section 3 (task U5), in the **projector form**
that user decision D32 (2026-09-01) fixed: no eigenvector is ever selected, so no measurability of
an eigenvector selection is needed. The consumer is the hypothesis `hκ` of U7 part 3.

The headline, for the Gram `Zᵀ Z` of a canonical Gaussian `p × d` matrix, a sorted eigenvalue
index `k < min p d` and a unit direction `w`:

```
∫⁻ Z, ofReal (‖specProjIdx (Zᵀ * Z) _ k w‖ ^ 2) ∂(gaussianMatrix p d) ≤ ofReal (1 / d)
```

## The argument

1. **Right rotation moves the test direction** (deterministic, section 1). For `Oᵀ O = 1` the
   Gram of `Z O` is `Oᵀ (Zᵀ Z) O`. The two matrices have the same characteristic polynomial
   (`charpoly_conj`, from `Matrix.charpoly_mul_comm`), hence the same **sorted** spectrum
   (`eigenvalues₀_conj`), hence the same eigenvalue set at every index (`eigSetIdx_conj`). The
   spectral subspaces conjugate (`Symmetry.specSpace_conj`), so the projector at the index `k`
   conjugates too (`specProjIdx_conj`) and `‖specProjIdx ((Z O)ᵀ (Z O)) _ k w‖` equals
   `‖specProjIdx (Zᵀ Z) _ k (O w)‖`.
   This is the step where the index form needs more than `topProj_conj` of `Symmetry.lean`:
   `lamMax_conj` there is an `IsGreatest` argument that says nothing about index `k`.
2. **Right rotation preserves the law** (`Symmetry.measurePreserving_mul_right`), and the map
   `Z ↦ Z O` is a measurable equivalence, so the change of variables goes through
   `MeasurePreserving.lintegral_comp_emb`, which asks **nothing** of the integrand.
3. **Exchangeability** (section 4). There is no protected direction here (the signal is absent),
   so for any two unit vectors the Householder reflection with axis `w - w'` sends `w` to `w'`
   and steps 1 and 2 give equal integrals.
4. **Bessel** (section 3). On the almost sure event `SimpleSpec (Zᵀ Z) _ (k+1)` of
   `SimplicityR.simpleSpec_ae_gaussianMatrix` the projector at index `k` is the rank one
   projector on `vEig`, so the sum over any orthonormal family is at most `‖vEig‖² = 1`. The
   tie set is a null set and is handled by `filter_upwards`, never assumed away.
5. **Count** (section 5). The family is a full orthonormal basis of `ℝ^d`, which has `d`
   elements and no protected direction to drop, so the bound is `1/d`, not the `1/(d-1)` of
   `Symmetry.lintegral_overlap_le`.

## Measurability: not needed

The file proves **no** measurability of `Z ↦ ‖specProjIdx (Zᵀ Z) _ k w‖ ^ 2` (task B2a will do
that by a polynomial route). Every step above avoids it:

* `MeasurePreserving.lintegral_comp_emb` takes the integrand with no hypothesis, because the map
  is a measurable embedding;
* `le_lintegral_add` gives `∑ ∫⁻ ≤ ∫⁻ ∑` (`sum_lintegral_le`), the direction this proof needs,
  with no hypothesis (only the reverse inequality needs measurability);
* `lintegral_mono_ae` and `lintegral_const_mul'` are unconditional;
* `MeasureTheory.lintegral_prod_le` gives `∫⁻ over the product ≤ ∫⁻ ∫⁻` with no hypothesis, so
  even the conditional twin of section 6 is free of it.

## Section 6: the conditional twin

`hκ` of U7 part 3 tests the projector against the columns of `Q`, which are random and
independent of the block. Section 6 states that case on a product measure: the block carries
`gaussianMatrix p d` and the direction is any function `q` of an independent variable. The bound
is `(∫⁻ ‖q‖ ²) / d`, an integral bound and not a supremum bound, so a Gaussian direction of
unbounded length is covered; `..._of_bound` is the corollary for `‖q ω‖ ² ≤ K` on a probability
space. Both the iterated form (condition on the direction, then average) and the joint form are
given; the joint form follows from the iterated one, so the order of integration costs nothing.
`q` carries no measurability hypothesis, because none of the steps uses one.

No `sorry`, no `axiom`, no edit to any existing file.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### 1. The spectral data at one index, under an orthogonal conjugation -/

section Conj

variable {d : ℕ} {O : Matrix (Fin d) (Fin d) ℝ}

/-- Conjugation by an orthogonal matrix does not change the characteristic polynomial.
`Matrix.charpoly_mul_comm` moves the trailing `O` to the front, where it meets `Oᵀ`. -/
theorem charpoly_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ) :
    (Oᵀ * A * O).charpoly = A.charpoly := by
  rw [Matrix.charpoly_mul_comm (Oᵀ * A) O, ← Matrix.mul_assoc, mul_transpose_of_orth hO,
    Matrix.one_mul]

/-- Two Hermitian matrices with the same characteristic polynomial have the same sorted
spectrum. `Matrix.IsHermitian.sort_roots_charpoly_eq_eigenvalues₀` reads `eigenvalues₀` off the
roots of the characteristic polynomial, sorted downward; `List.ofFn` is injective. -/
theorem eigenvalues₀_congr_charpoly {A B : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A.charpoly = B.charpoly) : hA.eigenvalues₀ = hB.eigenvalues₀ := by
  rw [← List.ofFn_inj, ← hA.sort_roots_charpoly_eq_eigenvalues₀,
    ← hB.sort_roots_charpoly_eq_eigenvalues₀, h]

/-- **The sorted spectrum is invariant under conjugation by an orthogonal matrix.** This is the
index form of `Symmetry.lamMax_conj`, which only covers the index `0`. -/
theorem eigenvalues₀_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) : hB.eigenvalues₀ = hA.eigenvalues₀ :=
  eigenvalues₀_congr_charpoly hB hA (charpoly_conj hO A)

/-- The eigenvalue set at a sorted index is invariant under conjugation. -/
theorem eigSetIdx_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) (k : ℕ) :
    eigSetIdx (Oᵀ * A * O) hB k = eigSetIdx A hA k := by
  unfold eigSetIdx
  rw [eigenvalues₀_conj hO hA hB]

/-- The spectral projector at an arbitrary set of eigenvalues conjugates. Mirror of
`Symmetry.topProj_conj`, with `specSpace_conj` in place of `topSpace_conj`. -/
theorem specProj_conj (hO : Oᵀ * O = 1) (A : Matrix (Fin d) (Fin d) ℝ) (S : Set ℝ)
    (x : EuclideanSpace ℝ (Fin d)) :
    specProj (Oᵀ * A * O) S x =
      rotIso Oᵀ (transpose_orth hO) (specProj A S (rotIso O hO x)) := by
  change (specSpace (Oᵀ * A * O) S).starProjection x = _
  rw [specSpace_conj hO A S, Submodule.starProjection_map_apply, rotIso_transpose_symm hO]
  rfl

/-- Equal matrices have the same projector at an index; the two `IsHermitian` proofs are
irrelevant. Mirror of `Spectral.topProj_congr`. -/
theorem specProjIdx_congr {A B : Matrix (Fin d) (Fin d) ℝ} (h : A = B) (hA : A.IsHermitian)
    (hB : B.IsHermitian) (k : ℕ) : specProjIdx A hA k = specProjIdx B hB k := by
  subst h
  rfl

/-- **The projector at a sorted index conjugates.** -/
theorem specProjIdx_conj (hO : Oᵀ * O = 1) {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hB : (Oᵀ * A * O).IsHermitian) (k : ℕ) (x : EuclideanSpace ℝ (Fin d)) :
    specProjIdx (Oᵀ * A * O) hB k x =
      rotIso Oᵀ (transpose_orth hO) (specProjIdx A hA k (rotIso O hO x)) := by
  rw [specProjIdx, specProjIdx, eigSetIdx_conj hO hA hB k, specProj_conj hO]

end Conj

/-- **Right multiplication by an orthogonal matrix moves the test direction.** The index form
of `Symmetry.overlap_mul_right`. -/
theorem normSq_specProjIdx_mul_right {n d : ℕ} (Z : Matrix (Fin n) (Fin d) ℝ)
    {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1) (k : ℕ)
    (w : EuclideanSpace ℝ (Fin d)) :
    ‖specProjIdx ((Z * O)ᵀ * (Z * O)) (isHermitian_transpose_mul_self (Z * O)) k w‖
      = ‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (rotIso O hO w)‖ := by
  have hgram : (Z * O)ᵀ * (Z * O) = Oᵀ * (Zᵀ * Z) * O := by
    simp [Matrix.transpose_mul, Matrix.mul_assoc]
  have hA : (Zᵀ * Z).IsHermitian := isHermitian_transpose_mul_self Z
  have hB : (Oᵀ * (Zᵀ * Z) * O).IsHermitian := by
    rw [← hgram]
    exact isHermitian_transpose_mul_self _
  rw [specProjIdx_congr hgram (isHermitian_transpose_mul_self (Z * O)) hB k,
    specProjIdx_conj hO hA hB k w, LinearIsometryEquiv.norm_map]

/-! ### 2. `Z ↦ Z O` is a measurable equivalence -/

section MeasEmb

variable {d : ℕ}

/-- Right multiplication by a fixed matrix is measurable: every entry of the product is a
finite sum of coordinates times constants. -/
theorem measurable_mul_right (P : Matrix (Fin d) (Fin d) ℝ) (n : ℕ) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => Z * P := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun k => ?_
  change Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => ∑ j, Z i j * P j k
  refine Finset.measurable_sum _ fun j _ => ?_
  have hij : Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => Z i j :=
    (measurable_pi_apply j).comp (measurable_pi_apply i)
  exact hij.mul_const (P j k)

/-- Right multiplication by an orthogonal matrix, as a measurable equivalence. Its inverse is
right multiplication by the transpose. -/
def mulRightEquiv {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1) (n : ℕ) :
    Matrix (Fin n) (Fin d) ℝ ≃ᵐ Matrix (Fin n) (Fin d) ℝ where
  toEquiv :=
    { toFun := fun Z => Z * O
      invFun := fun Z => Z * Oᵀ
      left_inv := fun Z => by
        change Z * O * Oᵀ = Z
        rw [Matrix.mul_assoc, mul_transpose_of_orth hO, Matrix.mul_one]
      right_inv := fun Z => by
        change Z * Oᵀ * O = Z
        rw [Matrix.mul_assoc, hO, Matrix.mul_one] }
  measurable_toFun := measurable_mul_right O n
  measurable_invFun := measurable_mul_right Oᵀ n

/-- The change of variables of section 4 is a measurable embedding, so it needs no hypothesis
on the integrand. -/
theorem measurableEmbedding_mul_right {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ * O = 1)
    (n : ℕ) : MeasurableEmbedding (fun Z : Matrix (Fin n) (Fin d) ℝ => Z * O) :=
  (mulRightEquiv hO n).measurableEmbedding

end MeasEmb

/-! ### 3. Bessel at one index, on the simplicity event -/

/-- **Bessel.** If the eigenvalue at the sorted index `k` is simple, the projector there is rank
one, so the squared norms of its action on an orthonormal family sum to at most `1`. Index form
of `Symmetry.sum_overlap_le_one`. The simplicity hypothesis is `SimpleSpec A hA (k+1)`, which
`SimplicityR.simpleSpec_ae_gaussianMatrix` supplies almost surely in both regimes. -/
theorem sum_normSq_specProjIdx_le_one {d : ℕ} {ι : Type*} [Fintype ι]
    {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian} {k : ℕ}
    (hsimple : SimpleSpec A hA (k + 1)) (hk : k < d)
    {e : ι → EuclideanSpace ℝ (Fin d)} (he : Orthonormal ℝ e) :
    ∑ i, ‖specProjIdx A hA k (e i)‖ ^ 2 ≤ 1 := by
  have hkc : k < Fintype.card (Fin d) := by simpa using hk
  have hval : ∀ i, ‖specProjIdx A hA k (e i)‖ ^ 2 = ‖⟪e i, vEig A hA k⟫_ℝ‖ ^ 2 := by
    intro i
    rw [specProjIdx_eq_rankOne hsimple (Nat.lt_succ_self k) (e i), norm_smul,
      norm_vEig A hA hkc, mul_one, real_inner_comm]
  calc ∑ i, ‖specProjIdx A hA k (e i)‖ ^ 2 = ∑ i, ‖⟪e i, vEig A hA k⟫_ℝ‖ ^ 2 :=
        Finset.sum_congr rfl fun i _ => hval i
    _ ≤ ‖vEig A hA k‖ ^ 2 := Orthonormal.sum_inner_products_le _ he
    _ = 1 := by rw [norm_vEig A hA hkc]; norm_num

/-! ### 4. Two general facts about the lower Lebesgue integral -/

/-- Superadditivity of the lower Lebesgue integral over a finite sum. This is the direction that
needs **no** measurability; the reverse inequality does. Iterated `le_lintegral_add`. -/
theorem sum_lintegral_le {α : Type*} [MeasurableSpace α] {ν : Measure α} {ι : Type*}
    (s : Finset ι) (f : ι → α → ℝ≥0∞) :
    ∑ i ∈ s, ∫⁻ a, f i a ∂ν ≤ ∫⁻ a, ∑ i ∈ s, f i a ∂ν := by
  classical
  refine Finset.induction_on s ?_ ?_
  · simp
  · intro i t hi ih
    rw [Finset.sum_insert hi]
    refine le_trans (add_le_add (le_refl _) ih) (le_trans (le_lintegral_add _ _) (le_of_eq ?_))
    exact lintegral_congr fun a => by rw [Finset.sum_insert hi]

/-! ### 5. Exchangeability of the test direction, and the bound `1/d` -/

/-- **The law of `‖specProjIdx (Zᵀ Z) _ k w‖ ²` is the same for every unit `w`.** There is no
signal and therefore no direction to protect, so the Householder reflection needs no side
condition. Stated as equality of the two `lintegral`s. -/
theorem lintegral_normSq_specProjIdx_eq_of_unit {p d : ℕ} (k : ℕ)
    {w w' : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) (hw' : ‖w'‖ = 1) :
    ∫⁻ Z, ENNReal.ofReal (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w‖ ^ 2)
        ∂(gaussianMatrix p d)
      = ∫⁻ Z, ENNReal.ofReal (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w'‖ ^ 2)
        ∂(gaussianMatrix p d) := by
  by_cases hww : w = w'
  · rw [hww]
  have hdot : ∀ x y : EuclideanSpace ℝ (Fin d),
      WithLp.ofLp x ⬝ᵥ WithLp.ofLp y = ⟪x, y⟫_ℝ := fun x y =>
    (inner_euclidean_eq_dotProduct x y).symm
  have hww1 : WithLp.ofLp w ⬝ᵥ WithLp.ofLp w = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw]
    norm_num
  have hww'1 : WithLp.ofLp w' ⬝ᵥ WithLp.ofLp w' = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw']
    norm_num
  have hne : WithLp.ofLp w ≠ WithLp.ofLp w' := fun h => hww (by
    have := congrArg (WithLp.toLp 2) h
    simpa using this)
  set a : Fin d → ℝ := WithLp.ofLp w - WithLp.ofLp w' with ha
  set O : Matrix (Fin d) (Fin d) ℝ := householder a with hOdef
  have ha0 : a ⬝ᵥ a ≠ 0 := fun h => (sub_ne_zero.mpr hne) (dotProduct_self_eq_zero.mp h)
  have hO : Oᵀ * O = 1 := householder_orth ha0
  have hOw : rotIso O hO w = w' := by
    rw [rotIso_apply]
    exact congrArg (WithLp.toLp 2) (householder_sub_apply hww1 hww'1 hne)
  have key : ∀ Z : Matrix (Fin p) (Fin d) ℝ,
      ‖specProjIdx ((Z * O)ᵀ * (Z * O)) (isHermitian_transpose_mul_self (Z * O)) k w‖
        = ‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w'‖ := by
    intro Z
    rw [normSq_specProjIdx_mul_right Z hO k w, hOw]
  calc ∫⁻ Z, ENNReal.ofReal
          (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w‖ ^ 2)
          ∂(gaussianMatrix p d)
      = ∫⁻ Z, ENNReal.ofReal
          (‖specProjIdx ((Z * O)ᵀ * (Z * O))
            (isHermitian_transpose_mul_self (Z * O)) k w‖ ^ 2) ∂(gaussianMatrix p d) :=
        ((measurePreserving_mul_right hO p).lintegral_comp_emb
          (measurableEmbedding_mul_right hO p) _).symm
    _ = _ := by simp only [key]

/-- **U5(b), the unconditional half.** For a Gaussian `p × d` matrix, a sorted eigenvalue index
`k < min p d` and a unit direction `w`, the mean squared norm of the spectral projector at that
index, applied to `w`, is at most `1/d`.

The index bound `k < min p d` is what `SimplicityR.simpleSpec_ae_gaussianMatrix` needs; `k < d`
alone is not enough, because a `d × d` Gram built from a shorter block has a repeated zero
eigenvalue (see the deviation note in `SimplicityR.lean`). -/
theorem lintegral_normSq_specProjIdx_le {p d k : ℕ} (hk : k < min p d)
    {w : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) :
    ∫⁻ Z, ENNReal.ofReal (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w‖ ^ 2)
      ∂(gaussianMatrix p d) ≤ ENNReal.ofReal (1 / (d : ℝ)) := by
  have hp : 0 < p := by omega
  have hd : 0 < d := by omega
  set b := EuclideanSpace.basisFun (Fin d) ℝ with hb
  have hon : Orthonormal ℝ (⇑b) := b.orthonormal
  set I : ℝ≥0∞ := ∫⁻ Z, ENNReal.ofReal
      (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k w‖ ^ 2)
      ∂(gaussianMatrix p d) with hI
  have hsame : ∀ j, ∫⁻ Z, ENNReal.ofReal
      (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2)
      ∂(gaussianMatrix p d) = I := fun j =>
    (lintegral_normSq_specProjIdx_eq_of_unit k hw (hon.1 j)).symm
  -- the tie set is a null set: on its complement the projector at index `k` is rank one
  have hae : ∀ᵐ Z ∂(gaussianMatrix p d),
      ∑ j, ENNReal.ofReal
        (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2) ≤ 1 := by
    filter_upwards [simpleSpec_ae_gaussianMatrix p d hp hd (k + 1) (by omega)] with Z hZ
    have hbes : ∑ j, ‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2
        ≤ 1 := sum_normSq_specProjIdx_le_one hZ (by omega) hon
    calc ∑ j, ENNReal.ofReal
            (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2)
        = ENNReal.ofReal
            (∑ j, ‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2) :=
          (ENNReal.ofReal_sum_of_nonneg fun j _ => by positivity).symm
      _ ≤ 1 := ENNReal.ofReal_le_one.mpr hbes
  have hsum : ∑ j, ∫⁻ Z, ENNReal.ofReal
      (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2)
      ∂(gaussianMatrix p d) ≤ 1 := by
    calc ∑ j, ∫⁻ Z, ENNReal.ofReal
            (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2)
            ∂(gaussianMatrix p d)
        ≤ ∫⁻ Z, ∑ j, ENNReal.ofReal
            (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (b j)‖ ^ 2)
            ∂(gaussianMatrix p d) :=
          sum_lintegral_le _ _
      _ ≤ ∫⁻ _, (1 : ℝ≥0∞) ∂(gaussianMatrix p d) := lintegral_mono_ae hae
      _ = 1 := by simp
  simp only [hsame, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
  have hdne : (d : ℝ≥0∞) ≠ 0 := by
    simpa using hd.ne'
  have hdtop : (d : ℝ≥0∞) ≠ ⊤ := ENNReal.natCast_ne_top d
  have hdiv : I ≤ 1 / (d : ℝ≥0∞) := by
    rw [ENNReal.le_div_iff_mul_le (Or.inl hdne) (Or.inl hdtop), mul_comm]
    exact hsum
  refine hdiv.trans (le_of_eq ?_)
  have hdpos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd
  rw [ENNReal.ofReal_div_of_pos hdpos, ENNReal.ofReal_one, ENNReal.ofReal_natCast]

/-- The same bound for a direction of arbitrary length: the projector is linear, so the bound
scales by `‖x‖²`. This is the form the conditional twin of section 6 consumes. -/
theorem lintegral_normSq_specProjIdx_le_of_norm {p d k : ℕ} (hk : k < min p d)
    (x : EuclideanSpace ℝ (Fin d)) :
    ∫⁻ Z, ENNReal.ofReal (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k x‖ ^ 2)
      ∂(gaussianMatrix p d) ≤ ENNReal.ofReal (‖x‖ ^ 2 / (d : ℝ)) := by
  rcases eq_or_ne x 0 with rfl | hx
  · simp
  have hxn : ‖x‖ ≠ 0 := norm_ne_zero_iff.mpr hx
  have hun : ‖(‖x‖⁻¹ • x : EuclideanSpace ℝ (Fin d))‖ = 1 := by
    rw [norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hxn]
  have hxu : x = ‖x‖ • (‖x‖⁻¹ • x : EuclideanSpace ℝ (Fin d)) := by
    rw [smul_smul, mul_inv_cancel₀ hxn, one_smul]
  have hpt : ∀ Z : Matrix (Fin p) (Fin d) ℝ,
      ENNReal.ofReal (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k x‖ ^ 2)
        = ENNReal.ofReal (‖x‖ ^ 2) * ENNReal.ofReal
            (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k
              (‖x‖⁻¹ • x : EuclideanSpace ℝ (Fin d))‖ ^ 2) := by
    intro Z
    rw [← ENNReal.ofReal_mul (by positivity)]
    congr 1
    conv_lhs => rw [hxu]
    rw [map_smul, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]
  rw [lintegral_congr hpt, lintegral_const_mul' _ _ ENNReal.ofReal_ne_top]
  calc ENNReal.ofReal (‖x‖ ^ 2) * ∫⁻ Z, ENNReal.ofReal
          (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k
            (‖x‖⁻¹ • x : EuclideanSpace ℝ (Fin d))‖ ^ 2) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal (‖x‖ ^ 2) * ENNReal.ofReal (1 / (d : ℝ)) :=
        mul_le_mul' (le_refl _) (lintegral_normSq_specProjIdx_le hk hun)
    _ = ENNReal.ofReal (‖x‖ ^ 2 / (d : ℝ)) := by
        rw [← ENNReal.ofReal_mul (by positivity), mul_one_div]

/-! ### 6. The conditional twin: a direction independent of the block -/

section Indep

variable {p d k : ℕ} {Ω : Type*} [MeasurableSpace Ω]

/-- **The conditional twin, iterated form.** Condition on the direction, then average. For a
direction `q` on any measure space, independent of the Gaussian block, the mean of
`‖specProjIdx (Zᵀ Z) _ k (q ω)‖ ²` is at most `(∫⁻ ‖q‖ ²) / d`. The bound on `q` is an integral
bound, so a Gaussian direction of unbounded length is covered. -/
theorem lintegral_lintegral_normSq_specProjIdx_le (hk : k < min p d) (ν : Measure Ω)
    (q : Ω → EuclideanSpace ℝ (Fin d)) :
    ∫⁻ ω, (∫⁻ Z, ENNReal.ofReal
        (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (q ω)‖ ^ 2)
        ∂(gaussianMatrix p d)) ∂ν
      ≤ (∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν) * ENNReal.ofReal (1 / (d : ℝ)) := by
  have hstep : ∀ ω, (∫⁻ Z, ENNReal.ofReal
      (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (q ω)‖ ^ 2)
      ∂(gaussianMatrix p d))
      ≤ ENNReal.ofReal (‖q ω‖ ^ 2) * ENNReal.ofReal (1 / (d : ℝ)) := by
    intro ω
    refine (lintegral_normSq_specProjIdx_le_of_norm hk (q ω)).trans (le_of_eq ?_)
    rw [← ENNReal.ofReal_mul (by positivity), mul_one_div]
  calc ∫⁻ ω, (∫⁻ Z, ENNReal.ofReal
          (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (q ω)‖ ^ 2)
          ∂(gaussianMatrix p d)) ∂ν
      ≤ ∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) * ENNReal.ofReal (1 / (d : ℝ)) ∂ν :=
        lintegral_mono hstep
    _ = (∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν) * ENNReal.ofReal (1 / (d : ℝ)) :=
        lintegral_mul_const' _ _ ENNReal.ofReal_ne_top

/-- The bounded form of the conditional twin, on a probability space: `‖q ω‖ ² ≤ K` gives the
mean bound `K / d`. This is the shape `hκ` of U7 part 3 needs before Markov. -/
theorem lintegral_lintegral_normSq_specProjIdx_le_of_bound (hk : k < min p d) (ν : Measure Ω)
    [IsProbabilityMeasure ν] (q : Ω → EuclideanSpace ℝ (Fin d)) {K : ℝ} (hK0 : 0 ≤ K)
    (hK : ∀ ω, ‖q ω‖ ^ 2 ≤ K) :
    ∫⁻ ω, (∫⁻ Z, ENNReal.ofReal
        (‖specProjIdx (Zᵀ * Z) (isHermitian_transpose_mul_self Z) k (q ω)‖ ^ 2)
        ∂(gaussianMatrix p d)) ∂ν ≤ ENNReal.ofReal (K / (d : ℝ)) := by
  refine (lintegral_lintegral_normSq_specProjIdx_le hk ν q).trans ?_
  have h1 : ∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν ≤ ENNReal.ofReal K := by
    calc ∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν ≤ ∫⁻ _, ENNReal.ofReal K ∂ν :=
          lintegral_mono fun ω => ENNReal.ofReal_le_ofReal (hK ω)
      _ = ENNReal.ofReal K := by simp
  calc (∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν) * ENNReal.ofReal (1 / (d : ℝ))
      ≤ ENNReal.ofReal K * ENNReal.ofReal (1 / (d : ℝ)) := mul_le_mul' h1 (le_refl _)
    _ = ENNReal.ofReal (K / (d : ℝ)) := by
        rw [← ENNReal.ofReal_mul hK0, mul_one_div]

/-- **The conditional twin, joint form.** The same bound for the integral against the product
measure, which is the law of an independent pair. `MeasureTheory.lintegral_prod_le` needs no
measurability of the integrand, so this costs one swap of the coordinates. -/
theorem lintegral_prod_normSq_specProjIdx_le (hk : k < min p d) (ν : Measure Ω) [SFinite ν]
    (q : Ω → EuclideanSpace ℝ (Fin d)) :
    ∫⁻ zω : Matrix (Fin p) (Fin d) ℝ × Ω, ENNReal.ofReal
        (‖specProjIdx (zω.1ᵀ * zω.1) (isHermitian_transpose_mul_self zω.1) k (q zω.2)‖ ^ 2)
        ∂((gaussianMatrix p d).prod ν)
      ≤ (∫⁻ ω, ENNReal.ofReal (‖q ω‖ ^ 2) ∂ν) * ENNReal.ofReal (1 / (d : ℝ)) := by
  have hswap := lintegral_prod_swap (μ := gaussianMatrix p d) (ν := ν)
    (fun zω : Matrix (Fin p) (Fin d) ℝ × Ω => ENNReal.ofReal
      (‖specProjIdx (zω.1ᵀ * zω.1) (isHermitian_transpose_mul_self zω.1) k (q zω.2)‖ ^ 2))
  rw [← hswap]
  refine (lintegral_prod_le _).trans ?_
  exact lintegral_lintegral_normSq_specProjIdx_le hk ν q

end Indep

end StackedSVD
