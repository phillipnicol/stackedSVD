/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.PolynomialNull
import StackedSVD.LinAlg.TopProjPerturb
import StackedSVD.RMT

/-!
# Item S: the top eigenvalue of `Xᵀ X` is simple almost surely

This file discharges the `topSimple` field of `SpikedModel.SingleTableLaw` for Gaussian noise
(`notes/archive/rmt_S.md`, items S1 to S4; S5 and S6 are in `Prob/PolynomialNull.lean`).

## Route

Write `D M` (`StackedSVD.charRes`) for the resultant of the characteristic polynomial of a square
matrix `M` and its derivative. Over `ℝ` and for a monic polynomial of positive degree this
resultant is nonzero exactly when the polynomial is separable
(`charRes_ne_zero_iff_separable`), so `D M ≠ 0` says that `M` has pairwise distinct eigenvalues.

1. Pairwise distinct eigenvalues give a simple top eigenvalue
   (`topSimple_of_injective`), through `LinearMap.IsSymmetric.card_filter_eigenvalues_eq`.
2. `charRes` is a polynomial expression in the matrix entries: it commutes with a ring
   homomorphism applied entrywise (`charRes_map`). Substituting the affine family
   `Z ↦ A + t • Z` into the universal matrix over `MvPolynomial (Fin n × Fin d) ℝ`
   (`affMat`) turns `Z ↦ D (Gram (A + t • Z))` into the evaluation of one multivariate
   polynomial.
3. That polynomial is not the zero polynomial: the witness `witMat` has a diagonal Gram matrix
   with the distinct positive entries `1, 4, 9, …`.
4. `ae_eval_ne_zero_gaussianMatrix` (from `Prob/PolynomialNull.lean`) makes the zero set null.

The size split is forced. For `d - n ≥ 2` the discriminant of `charpoly (Xᵀ X)` vanishes
identically (`0` is a repeated eigenvalue), so for `n < d` the argument runs on `X Xᵀ` and
transfers with `finrank_eigenspace_transpose_mul_eq`. The transfer needs a positive top
eigenvalue, which is why the `n ≤ d` polynomial carries the extra factor `det (X Xᵀ)`.

## Rank `r` remark (no proof here)

The same three steps give pairwise distinct top `r` eigenvalues for any `r ≤ min n d`: the
witness `witMat` already has `min n d` distinct positive Gram eigenvalues, so the polynomial
`charRes` is still nonzero, and `finrank_eigenspace_eq_card` still gives
`finrank (eigenspace (toOp A) (hA.eigenvalues₀ i)) = 1` for every index `i`, not only for
`i = 0`. Only the statement of the conclusion changes (a rank `r` spectral projector in place
of `TopSimple`); no new probabilistic input is needed.

The duplicate `isHermitian_mul_transpose_self` is removed; the one copy lives in
`Defs.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. Simplicity from distinct eigenvalues -/

section Deterministic

variable {d : ℕ}

/-- The top eigenspace is the eigenspace at the largest eigenvalue. -/
private theorem topSpaceEq (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    topSpace A hA = Module.End.eigenspace (toOp A) (lamMax A hA) := by
  unfold topSpace specSpace
  simp

private theorem lamMaxEq' {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d)
    (h0 : (0 : ℕ) < Fintype.card (Fin d)) : lamMax A hA = hA.eigenvalues₀ ⟨0, h0⟩ := by
  rw [lamMax, dif_pos hd]

private theorem eig_le_lamMax {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d)
    (i : Fin (Fintype.card (Fin d))) : hA.eigenvalues₀ i ≤ lamMax A hA := by
  rw [lamMax, dif_pos hd]
  exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

private theorem eigenvalues₀_eq_eigenvalues' {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ k = hA.eigenvalues (Fintype.equivOfCardEq (Fintype.card_fin _) k) := by
  simp [Matrix.IsHermitian.eigenvalues]

/-- The dimension of an eigenspace is the number of indices carrying that eigenvalue. This is
`LinearMap.IsSymmetric.card_filter_eigenvalues_eq` in the matrix notation of `Defs.lean`. -/
theorem finrank_eigenspace_eq_card {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (t : ℝ) :
    Module.finrank ℝ (Module.End.eigenspace (toOp A) t)
      = Finset.card {i | hA.eigenvalues₀ i = t} :=
  ((symmOp hA).card_filter_eigenvalues_eq finrank_euclideanSpace t).symm

/-- Distinct eigenvalues make every eigenspace at most one dimensional. -/
theorem finrank_eigenspace_le_one {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues₀) (t : ℝ) :
    Module.finrank ℝ (Module.End.eigenspace (toOp A) t) ≤ 1 := by
  rw [finrank_eigenspace_eq_card hA t]
  refine Finset.card_le_one.mpr fun a ha b hb => ?_
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
  exact hinj (ha.trans hb.symm)

/-- An index carrying the value `t` makes the eigenspace at `t` nontrivial. -/
theorem one_le_finrank_eigenspace {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) {t : ℝ}
    (i : Fin (Fintype.card (Fin d))) (hi : hA.eigenvalues₀ i = t) :
    1 ≤ Module.finrank ℝ (Module.End.eigenspace (toOp A) t) := by
  rw [finrank_eigenspace_eq_card hA t]
  exact Finset.card_pos.mpr ⟨i, by simp [hi]⟩

/-- **S1.** Pairwise distinct eigenvalues make the top eigenvalue simple. -/
theorem topSimple_of_injective {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (hd : 0 < d)
    (hinj : Function.Injective hA.eigenvalues₀) : TopSimple A hA := by
  have h0 : (0 : ℕ) < Fintype.card (Fin d) := by simpa using hd
  unfold TopSimple
  rw [topSpaceEq]
  exact le_antisymm (finrank_eigenspace_le_one hA hinj _)
    (one_le_finrank_eigenspace hA ⟨0, h0⟩ (lamMaxEq' hA hd h0).symm)

/-- **S2.** A separable characteristic polynomial makes the eigenvalue list injective. -/
theorem injective_eigenvalues₀_of_separable {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hsep : A.charpoly.Separable) : Function.Injective hA.eigenvalues₀ := by
  have hnd : A.charpoly.roots.Nodup := Polynomial.nodup_roots hsep
  rw [hA.roots_charpoly_eq_eigenvalues₀] at hnd
  intro i j hij
  exact Multiset.inj_on_of_nodup_map hnd i (by simp) j (by simp) (by simp [hij])

end Deterministic

/-! ### 2. The resultant of the characteristic polynomial and its derivative -/

/-- `charRes M` is the resultant of `charpoly M` and its derivative, with the degrees written
out so that the definition commutes with ring homomorphisms. Over `ℝ` it is nonzero exactly
when the eigenvalues of `M` are pairwise distinct. -/
noncomputable def charRes {R : Type*} [CommRing R] {m : ℕ} (M : Matrix (Fin m) (Fin m) R) : R :=
  Polynomial.resultant M.charpoly (Polynomial.derivative M.charpoly) m (m - 1)

/-- `charRes` is a polynomial in the entries: it commutes with an entrywise ring
homomorphism. -/
theorem charRes_map {R S : Type*} [CommRing R] [CommRing S] {m : ℕ}
    (M : Matrix (Fin m) (Fin m) R) (φ : R →+* S) : φ (charRes M) = charRes (M.map φ) := by
  rw [charRes, charRes, Matrix.charpoly_map, Polynomial.derivative_map,
    Polynomial.resultant_map_map]

/-- **S3.** Over `ℝ`, `charRes M ≠ 0` says exactly that `charpoly M` is separable. -/
theorem charRes_ne_zero_iff_separable {m : ℕ} (_hm : 0 < m) (M : Matrix (Fin m) (Fin m) ℝ) :
    charRes M ≠ 0 ↔ M.charpoly.Separable := by
  have hmonic : M.charpoly.Monic := M.charpoly_monic
  have hdeg : M.charpoly.natDegree = m := by simp
  have hderiv : (Polynomial.derivative M.charpoly).natDegree = m - 1 := by
    rw [Polynomial.natDegree_derivative, hdeg]
  have hres : Polynomial.resultant M.charpoly (Polynomial.derivative M.charpoly) = charRes M := by
    rw [charRes, hdeg, hderiv]
  rw [← hres, ← isUnit_iff_ne_zero]
  exact Polynomial.isUnit_resultant_iff_isCoprime hmonic

/-- **S1 to S3 combined.** -/
theorem topSimple_of_charRes_ne_zero {d : ℕ} {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (hd : 0 < d) (h : charRes A ≠ 0) : TopSimple A hA :=
  topSimple_of_injective hA hd
    (injective_eigenvalues₀_of_separable hA ((charRes_ne_zero_iff_separable hd A).mp h))

/-- A diagonal matrix with distinct entries has a nonzero `charRes`. -/
theorem charRes_diagonal_ne_zero {m : ℕ} (hm : 0 < m) {f : Fin m → ℝ}
    (hf : Function.Injective f) : charRes (Matrix.diagonal f) ≠ 0 := by
  rw [charRes_ne_zero_iff_separable hm, Matrix.charpoly_diagonal]
  exact Polynomial.separable_prod_X_sub_C_iff.mpr hf

/-! ### 3. Transfer between `Xᵀ X` and `X Xᵀ` -/


/-- **The `n ≤ d` case.** If the `n × n` Gram matrix `X Xᵀ` has pairwise distinct nonzero
eigenvalues, then the top eigenvalue of the `d × d` Gram matrix `Xᵀ X` is simple. -/
theorem topSimple_of_gram_left {n d : ℕ} (hn : 0 < n) (hd : 0 < d)
    (X : Matrix (Fin n) (Fin d) ℝ) (hres : charRes (X * Xᵀ) ≠ 0) (hdet : (X * Xᵀ).det ≠ 0) :
    TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
  have hB : (X * Xᵀ).IsHermitian := isHermitian_mul_transpose_self X
  have hA : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  have hn0 : (0 : ℕ) < Fintype.card (Fin n) := by simpa using hn
  have hd0 : (0 : ℕ) < Fintype.card (Fin d) := by simpa using hd
  have hinjB : Function.Injective hB.eigenvalues₀ :=
    injective_eigenvalues₀_of_separable hB ((charRes_ne_zero_iff_separable hn _).mp hres)
  have hpsd : (X * Xᵀ).PosSemidef := by simpa using Matrix.posSemidef_self_mul_conjTranspose X
  have hnn : ∀ i, 0 ≤ hB.eigenvalues₀ i := by
    intro i
    rw [eigenvalues₀_eq_eigenvalues' hB i]
    exact hpsd.eigenvalues_nonneg _
  have hne : ∀ i, hB.eigenvalues₀ i ≠ 0 := by
    intro i hi
    apply hdet
    rw [hB.det_eq_prod_eigenvalues]
    refine Finset.prod_eq_zero (Finset.mem_univ
      (Fintype.equivOfCardEq (Fintype.card_fin _) i)) ?_
    rw [← eigenvalues₀_eq_eigenvalues' hB i, hi]
    simp
  have hbpos : 0 < lamMax (X * Xᵀ) hB := by
    rw [lamMaxEq' hB hn hn0]
    exact lt_of_le_of_ne (hnn _) (Ne.symm (hne _))
  have hbeig : 1 ≤ Module.finrank ℝ
      (Module.End.eigenspace (toOp (Xᵀ * X)) (lamMax (X * Xᵀ) hB)) := by
    rw [finrank_eigenspace_transpose_mul_eq X hbpos.ne']
    exact one_le_finrank_eigenspace hB ⟨0, hn0⟩ (lamMaxEq' hB hn hn0).symm
  have hble : lamMax (X * Xᵀ) hB ≤ lamMax (Xᵀ * X) hA := by
    have h3 : 0 < Finset.card {i | hA.eigenvalues₀ i = lamMax (X * Xᵀ) hB} := by
      rw [← finrank_eigenspace_eq_card hA]
      exact hbeig
    obtain ⟨i, hi⟩ := Finset.card_pos.mp h3
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hi
    rw [← hi]
    exact eig_le_lamMax hA hd i
  have hapos : 0 < lamMax (Xᵀ * X) hA := lt_of_lt_of_le hbpos hble
  unfold TopSimple
  rw [topSpaceEq]
  refine le_antisymm ?_ (one_le_finrank_eigenspace hA ⟨0, hd0⟩ (lamMaxEq' hA hd hd0).symm)
  rw [finrank_eigenspace_transpose_mul_eq X hapos.ne']
  exact finrank_eigenspace_le_one hB hinjB _

/-! ### 4. The affine family as a matrix over a polynomial ring -/

/-- The universal matrix of the affine family `Z ↦ A + t • Z`: the entry `(i, j)` is the
polynomial `A i j + t X_{(i,j)}`. -/
noncomputable def affMat {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    Matrix (Fin n) (Fin d) (MvPolynomial (Fin n × Fin d) ℝ) :=
  Matrix.of fun i j => MvPolynomial.C (A i j) + MvPolynomial.C t * MvPolynomial.X (i, j)

theorem affMat_map {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ)
    (Z : Matrix (Fin n) (Fin d) ℝ) :
    (affMat A t).map (MvPolynomial.eval fun ij : Fin n × Fin d => Z ij.1 ij.2) = A + t • Z := by
  ext i j
  simp [affMat]

/-- The polynomial used when `d ≤ n`: `charRes` of the `d × d` Gram matrix. -/
noncomputable def gramPolyRight {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    MvPolynomial (Fin n × Fin d) ℝ :=
  charRes ((affMat A t)ᵀ * affMat A t)

/-- The polynomial used when `n ≤ d`: `charRes` of the `n × n` Gram matrix times its
determinant, so that the top eigenvalue is also forced to be positive. -/
noncomputable def gramPolyLeft {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    MvPolynomial (Fin n × Fin d) ℝ :=
  charRes (affMat A t * (affMat A t)ᵀ) * (affMat A t * (affMat A t)ᵀ).det

theorem eval_gramPolyRight {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ)
    (Z : Matrix (Fin n) (Fin d) ℝ) :
    MvPolynomial.eval (fun ij : Fin n × Fin d => Z ij.1 ij.2) (gramPolyRight A t)
      = charRes ((A + t • Z)ᵀ * (A + t • Z)) := by
  rw [gramPolyRight, charRes_map, Matrix.map_mul, Matrix.transpose_map, affMat_map]

theorem eval_gramPolyLeft {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ)
    (Z : Matrix (Fin n) (Fin d) ℝ) :
    MvPolynomial.eval (fun ij : Fin n × Fin d => Z ij.1 ij.2) (gramPolyLeft A t)
      = charRes ((A + t • Z) * (A + t • Z)ᵀ) * ((A + t • Z) * (A + t • Z)ᵀ).det := by
  simp only [gramPolyLeft, map_mul, charRes_map, RingHom.map_det, RingHom.mapMatrix_apply,
    Matrix.map_mul, Matrix.transpose_map, affMat_map]

/-! ### 5. The witness matrix -/

/-- Witness with a diagonal Gram matrix: `W i j = i + 1` on the diagonal, `0` elsewhere. -/
noncomputable def witMat (n d : ℕ) : Matrix (Fin n) (Fin d) ℝ :=
  Matrix.of fun i j => if (i : ℕ) = (j : ℕ) then ((i : ℕ) : ℝ) + 1 else 0

theorem witMat_transpose (n d : ℕ) : (witMat n d)ᵀ = witMat d n := by
  ext i j
  by_cases h : (j : ℕ) = (i : ℕ)
  · simp [witMat, h]
  · have h' : (i : ℕ) ≠ (j : ℕ) := fun hh => h hh.symm
    simp [witMat, h, h']

private theorem sum_wit {n d : ℕ} (h : d ≤ n) (j k : Fin d) :
    ∑ i : Fin n, (if (i : ℕ) = (j : ℕ) then ((i : ℕ) : ℝ) + 1 else 0)
        * (if (i : ℕ) = (k : ℕ) then ((i : ℕ) : ℝ) + 1 else 0)
      = if j = k then (((j : ℕ) : ℝ) + 1) ^ 2 else 0 := by
  classical
  rw [Finset.sum_eq_single (⟨(j : ℕ), lt_of_lt_of_le j.isLt h⟩ : Fin n)]
  · by_cases hjk : j = k
    · subst hjk; simp [sq]
    · have hv : (j : ℕ) ≠ (k : ℕ) := fun hh => hjk (Fin.ext hh)
      simp [hv, hjk]
  · intro i _ hi
    have hv : (i : ℕ) ≠ (j : ℕ) := fun hh => hi (Fin.ext hh)
    simp [hv]
  · intro hcon
    exact absurd (Finset.mem_univ _) hcon

theorem witMat_gram {n d : ℕ} (h : d ≤ n) :
    (witMat n d)ᵀ * witMat n d = Matrix.diagonal fun j : Fin d => (((j : ℕ) : ℝ) + 1) ^ 2 := by
  ext j k
  rw [Matrix.mul_apply, Matrix.diagonal_apply]
  rw [← sum_wit h j k]
  exact Finset.sum_congr rfl fun i _ => by simp [witMat]

/-- The Gram eigenvalues of the witness are pairwise distinct. -/
theorem wit_inj {m : ℕ} : Function.Injective fun j : Fin m => (((j : ℕ) : ℝ) + 1) ^ 2 := by
  intro a b hab
  simp only at hab
  have hpos : (0 : ℝ) < (((a : ℕ) : ℝ) + 1) + (((b : ℕ) : ℝ) + 1) := by positivity
  have h0 : ((((a : ℕ) : ℝ) + 1) - (((b : ℕ) : ℝ) + 1))
      * ((((a : ℕ) : ℝ) + 1) + (((b : ℕ) : ℝ) + 1)) = 0 := by linear_combination hab
  rcases mul_eq_zero.mp h0 with h | h
  · have : ((a : ℕ) : ℝ) = ((b : ℕ) : ℝ) := by linarith
    exact Fin.ext (by exact_mod_cast this)
  · linarith

/-- The Gram eigenvalues of the witness are nonzero. -/
theorem wit_det_ne_zero {m : ℕ} :
    (Matrix.diagonal fun j : Fin m => (((j : ℕ) : ℝ) + 1) ^ 2).det ≠ 0 := by
  rw [Matrix.det_diagonal]
  refine Finset.prod_ne_zero_iff.mpr fun i _ => ?_
  positivity

/-! ### 6. The polynomials are not identically zero -/

private theorem affine_hits {n d : ℕ} (A W : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) :
    A + t • (t⁻¹ • (W - A)) = W := by
  rw [smul_smul, mul_inv_cancel₀ ht, one_smul]
  abel

theorem gramPolyRight_ne_zero {n d : ℕ} (hd : 0 < d) (hdn : d ≤ n)
    (A : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) : gramPolyRight A t ≠ 0 := by
  intro h0
  have hev := eval_gramPolyRight A t (t⁻¹ • (witMat n d - A))
  rw [h0, affine_hits A (witMat n d) ht, witMat_gram hdn] at hev
  simp only [map_zero] at hev
  exact charRes_diagonal_ne_zero hd wit_inj hev.symm

theorem gramPolyLeft_ne_zero {n d : ℕ} (hn : 0 < n) (hnd : n ≤ d)
    (A : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) : gramPolyLeft A t ≠ 0 := by
  intro h0
  have hev := eval_gramPolyLeft A t (t⁻¹ • (witMat n d - A))
  rw [h0, affine_hits A (witMat n d) ht] at hev
  have hW : witMat n d * (witMat n d)ᵀ
      = Matrix.diagonal fun i : Fin n => (((i : ℕ) : ℝ) + 1) ^ 2 := by
    have h1 : witMat n d * (witMat n d)ᵀ = (witMat d n)ᵀ * witMat d n := by
      rw [← witMat_transpose n d, Matrix.transpose_transpose]
    rw [h1]
    exact witMat_gram hnd
  rw [hW] at hev
  simp only [map_zero] at hev
  exact mul_ne_zero (charRes_diagonal_ne_zero hn wit_inj) wit_det_ne_zero hev.symm

/-! ### 7. Item S for the canonical Gaussian matrix measure -/

/-- **Item S, canonical form.** For every fixed shift `A` and every nonzero scale `t`, the top
eigenvalue of the Gram matrix of `A + t • Z` is simple for almost every Gaussian `Z`. -/
theorem topSimple_ae_affine {n d : ℕ} (hn : 0 < n) (hd : 0 < d)
    (A : Matrix (Fin n) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) :
    ∀ᵐ Z ∂(gaussianMatrix n d),
      TopSimple ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) := by
  rcases le_total d n with hdn | hnd
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramPolyRight A t)
      (gramPolyRight_ne_zero hd hdn A ht)] with Z hZ
    rw [eval_gramPolyRight] at hZ
    exact topSimple_of_charRes_ne_zero _ hd hZ
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramPolyLeft A t)
      (gramPolyLeft_ne_zero hn hnd A ht)] with Z hZ
    rw [eval_gramPolyLeft] at hZ
    exact topSimple_of_gram_left hn hd _ (left_ne_zero_of_mul hZ) (right_ne_zero_of_mul hZ)

/-- **Item S for `gaussianMatrix` with no shift.** -/
theorem topSimple_ae_gaussianMatrix (n d : ℕ) (hn : 0 < n) (hd : 0 < d) :
    ∀ᵐ Z ∂(gaussianMatrix n d), TopSimple (Zᵀ * Z) (isHermitian_transpose_mul_self Z) := by
  have h := topSimple_ae_affine hn hd (0 : Matrix (Fin n) (Fin d) ℝ) (one_ne_zero (α := ℝ))
  filter_upwards [h] with Z hZ
  simpa using hZ

/-! ### 8. Item S for the spiked model -/

private theorem measurable_entry {n d : ℕ} (i : Fin n) (j : Fin d) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => Z i j := by
  change Measurable fun Z : Fin n → Fin d → ℝ => Z i j
  exact (measurable_pi_apply j).comp (measurable_pi_apply i)

private theorem measurable_affine {n d : ℕ} (A : Matrix (Fin n) (Fin d) ℝ) (t : ℝ) :
    Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A + t • Z := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun j => ?_
  change Measurable fun Z : Matrix (Fin n) (Fin d) ℝ => A i j + t * Z i j
  exact ((measurable_entry i j).const_mul t).const_add (A i j)

/-- **Item S for `SingleTableLaw`.** This is the `topSimple` field, for Gaussian noise. -/
theorem singleTableLaw_topSimple_of_gaussian {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ} (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (N : ℕ) :
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
    ((hG N).ae_iff hpm).mpr (topSimple_ae_affine (m.hn N) (m.hd N) A ht)
  filter_upwards [hae] with ω hω
  exact hω

end StackedSVD
