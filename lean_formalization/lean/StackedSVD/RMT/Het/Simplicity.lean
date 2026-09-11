/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Simplicity
import StackedSVD.StackSVDWeighted

/-!
# Item H10: the top eigenvalue of the weighted stack Gram matrix is simple almost surely

Task H10 of `notes/archive/plan_heterolaw_A.md` (section 3.7 and the H10 row of section 4). This
file discharges the `topSimple` field of `MultiTableModel.HeteroLaw` for jointly Gaussian noise and
any weight vector `w` with at least one nonzero entry.

## Route

`RMT/Simplicity.lean` proves item S for one table: the resultant `charRes` of the
characteristic polynomial of the Gram matrix and its derivative is a polynomial in the noise
entries, it is not the zero polynomial (the witness `witMat` has a diagonal Gram matrix with
the distinct positive entries `1, 4, 9, ...`), and the zero set of a nonzero polynomial is
null for a law that is absolutely continuous with respect to Lebesgue measure
(`Prob/PolynomialNull.lean`).

Two changes are needed here.

1. **Row scales.** The weighted stack is `X_w = ũ₀ vᵀ + Σ^{1/2} E` with
   `Σ^{1/2} = diag(w_i I_{n_i})`, so the noise enters row `r` with the factor
   `d^{-1/2} w_{i(r)}` instead of the single scale `d^{-1/2}`. Section 1 repeats the
   construction of `RMT/Simplicity.lean` for the family `Y ↦ A + rowScale s Y`, where `s` is
   any row scale vector with no zero entry. The witness step is the only place that uses
   `s r ≠ 0`: it inverts the affine map entrywise.
2. **Zero weights.** A table with `w_i = 0` contributes zero rows, and then the discriminant
   vanishes identically (`0` is a repeated eigenvalue of `X_w X_wᵀ`, and of `X_wᵀ X_w` when
   the surviving rows are fewer than `d`). Section 3 removes those tables first:
   `stackGramW_restrict` says that the Gram matrix of the weighted stack equals the Gram
   matrix of the weighted stack of the sub-collection `S = {i | w i ≠ 0}`, which is the same
   argument as `MultiTableModel.gram_restrict` for binary weights. `hw : ∃ i, w i ≠ 0` makes
   `S` nonempty, which the sub-collection needs. The hypothesis is necessary: at `w = 0` the
   weighted stack is the zero matrix and its top eigenvalue has multiplicity `d N`.

The Gaussian law of the **unweighted** stacked noise is `MultiTableModel.stack_law` (the
product of the per-table Gaussian laws, reindexed by `finSigmaFinEquiv`), so no new absolute
continuity statement about `Measure.pi` is needed here: `Prob/PolynomialNull.lean` already
carries `absolutelyContinuous_pi` and `gaussianMatrix_absolutelyContinuous`, and `stack_law`
packages them for the stacked index.

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/Simplicity.lean` exit 0; 0 `sorry`.
The build directory held a stale `StackSVDWeighted.olean` (built 08:17, before commit 84f4c25
added `restrict_jointGaussianNoise` at 10:40), so `Scalars.lean` and `StackSVDWeighted.lean`
were rebuilt first into a private `LEAN_PATH` overlay of symbolic links. A root `lake build`
regenerates them.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 1. The affine family with row scales -/

section RowScale

variable {p q : ℕ}

/-- `rowScale s Y` multiplies row `r` of `Y` by `s r`. With `s r = d^{-1/2} w_{i(r)}` this is
the map that turns the unweighted stacked noise into the weighted one. -/
def rowScale (s : Fin p → ℝ) (Y : Matrix (Fin p) (Fin q) ℝ) : Matrix (Fin p) (Fin q) ℝ :=
  Matrix.of fun r k => s r * Y r k

@[simp]
theorem rowScale_apply (s : Fin p → ℝ) (Y : Matrix (Fin p) (Fin q) ℝ) (r : Fin p) (k : Fin q) :
    rowScale s Y r k = s r * Y r k := rfl

/-- The universal matrix of the family `Y ↦ A + rowScale s Y`: the entry `(r, k)` is the
polynomial `A r k + s r * X_{(r,k)}`. This is `StackedSVD.affMat` with a per-row scale. -/
noncomputable def affMatRow (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ) :
    Matrix (Fin p) (Fin q) (MvPolynomial (Fin p × Fin q) ℝ) :=
  Matrix.of fun r k => MvPolynomial.C (A r k) + MvPolynomial.C (s r) * MvPolynomial.X (r, k)

theorem affMatRow_map (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    (affMatRow A s).map (MvPolynomial.eval fun rk : Fin p × Fin q => Y rk.1 rk.2)
      = A + rowScale s Y := by
  ext r k
  simp [affMatRow]

/-- The polynomial used when `q ≤ p`: `charRes` of the `q × q` Gram matrix. -/
noncomputable def gramRowPolyRight (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ) :
    MvPolynomial (Fin p × Fin q) ℝ :=
  charRes ((affMatRow A s)ᵀ * affMatRow A s)

/-- The polynomial used when `p ≤ q`: `charRes` of the `p × p` Gram matrix times its
determinant, so that the top eigenvalue is also forced to be positive. -/
noncomputable def gramRowPolyLeft (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ) :
    MvPolynomial (Fin p × Fin q) ℝ :=
  charRes (affMatRow A s * (affMatRow A s)ᵀ) * (affMatRow A s * (affMatRow A s)ᵀ).det

theorem eval_gramRowPolyRight (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    MvPolynomial.eval (fun rk : Fin p × Fin q => Y rk.1 rk.2) (gramRowPolyRight A s)
      = charRes ((A + rowScale s Y)ᵀ * (A + rowScale s Y)) := by
  rw [gramRowPolyRight, charRes_map, Matrix.map_mul, Matrix.transpose_map, affMatRow_map]

theorem eval_gramRowPolyLeft (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    MvPolynomial.eval (fun rk : Fin p × Fin q => Y rk.1 rk.2) (gramRowPolyLeft A s)
      = charRes ((A + rowScale s Y) * (A + rowScale s Y)ᵀ)
        * ((A + rowScale s Y) * (A + rowScale s Y)ᵀ).det := by
  simp only [gramRowPolyLeft, map_mul, charRes_map, RingHom.map_det, RingHom.mapMatrix_apply,
    Matrix.map_mul, Matrix.transpose_map, affMatRow_map]

/-- With no zero row scale the affine family reaches every matrix. -/
theorem rowScale_hits {s : Fin p → ℝ} (hs : ∀ r, s r ≠ 0) (A W : Matrix (Fin p) (Fin q) ℝ) :
    A + rowScale s (Matrix.of fun r k => (W r k - A r k) / s r) = W := by
  ext r k
  rw [Matrix.add_apply, rowScale_apply, Matrix.of_apply, mul_div_cancel₀ _ (hs r)]
  ring

theorem gramRowPolyRight_ne_zero (hq : 0 < q) (hqp : q ≤ p) (A : Matrix (Fin p) (Fin q) ℝ)
    {s : Fin p → ℝ} (hs : ∀ r, s r ≠ 0) : gramRowPolyRight A s ≠ 0 := by
  intro h0
  have hev := eval_gramRowPolyRight A s
    (Matrix.of fun r k => ((witMat p q) r k - A r k) / s r)
  rw [h0, rowScale_hits hs A (witMat p q), witMat_gram hqp] at hev
  simp only [map_zero] at hev
  exact charRes_diagonal_ne_zero hq wit_inj hev.symm

theorem gramRowPolyLeft_ne_zero (hp : 0 < p) (hpq : p ≤ q) (A : Matrix (Fin p) (Fin q) ℝ)
    {s : Fin p → ℝ} (hs : ∀ r, s r ≠ 0) : gramRowPolyLeft A s ≠ 0 := by
  intro h0
  have hev := eval_gramRowPolyLeft A s
    (Matrix.of fun r k => ((witMat p q) r k - A r k) / s r)
  rw [h0, rowScale_hits hs A (witMat p q)] at hev
  have hW : witMat p q * (witMat p q)ᵀ
      = Matrix.diagonal fun r : Fin p => (((r : ℕ) : ℝ) + 1) ^ 2 := by
    have h1 : witMat p q * (witMat p q)ᵀ = (witMat q p)ᵀ * witMat q p := by
      rw [← witMat_transpose p q, Matrix.transpose_transpose]
    rw [h1]
    exact witMat_gram hpq
  rw [hW] at hev
  simp only [map_zero] at hev
  exact mul_ne_zero (charRes_diagonal_ne_zero hp wit_inj) wit_det_ne_zero hev.symm

/-- **Item S with row scales.** For every fixed shift `A` and every row scale vector `s` with
no zero entry, the top eigenvalue of the Gram matrix of `A + rowScale s Y` is simple for
almost every Gaussian `Y`. -/
theorem topSimple_ae_rowScale (hp : 0 < p) (hq : 0 < q) (A : Matrix (Fin p) (Fin q) ℝ)
    {s : Fin p → ℝ} (hs : ∀ r, s r ≠ 0) :
    ∀ᵐ Y ∂(gaussianMatrix p q),
      TopSimple ((A + rowScale s Y)ᵀ * (A + rowScale s Y))
        (isHermitian_transpose_mul_self (A + rowScale s Y)) := by
  rcases le_total q p with hqp | hpq
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramRowPolyRight A s)
      (gramRowPolyRight_ne_zero hq hqp A hs)] with Y hY
    rw [eval_gramRowPolyRight] at hY
    exact topSimple_of_charRes_ne_zero _ hq hY
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramRowPolyLeft A s)
      (gramRowPolyLeft_ne_zero hp hpq A hs)] with Y hY
    rw [eval_gramRowPolyLeft] at hY
    exact topSimple_of_gram_left hp hq _ (left_ne_zero_of_mul hY) (right_ne_zero_of_mul hY)

/-- The map `Y ↦ A + rowScale s Y` is measurable. -/
theorem measurable_addRowScale (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ) :
    Measurable fun Y : Matrix (Fin p) (Fin q) ℝ => A + rowScale s Y := by
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun k => ?_
  change Measurable fun Y : Matrix (Fin p) (Fin q) ℝ => A r k + s r * Y r k
  have hentry : Measurable fun Y : Matrix (Fin p) (Fin q) ℝ => Y r k := by
    change Measurable fun Y : Fin p → Fin q → ℝ => Y r k
    exact (measurable_pi_apply k).comp (measurable_pi_apply r)
  exact (hentry.const_mul (s r)).const_add (A r k)

end RowScale

/-! ### 2. The weighted stack as a row-scaled affine family -/

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

omit [NeZero M] in
/-- The weighted stacked noise is the unweighted stacked noise with row `r` scaled by
`w_{i(r)}`, the block weight of that row. -/
theorem stackZW_eq_rowScale (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.stackZW w N ω
      = rowScale (fun r => w (finSigmaFinEquiv.symm r).1) (m.stackZ N ω) := by
  ext r k
  simp [stackZW, stackZ, rowScale, Matrix.reindex_apply, Matrix.submatrix_apply]

/-- The weighted stack is the affine family `Y ↦ A + rowScale s Y` of section 1 evaluated at
the unweighted stacked noise, with `s r = d^{-1/2} w_{i(r)}`. -/
theorem stackW_X_eq_addRowScale (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    (m.stackW w).X N ω
      = ((m.stackW w).θ • Matrix.vecMulVec (WithLp.ofLp ((m.stackW w).u N))
            (WithLp.ofLp ((m.stackW w).v N)))
        + rowScale (fun r => (Real.sqrt (d N))⁻¹ * w (finSigmaFinEquiv.symm r).1)
            (m.stackZ N ω) := by
  ext r k
  have hZ : (m.stackW w).Z N ω = m.stackZW w N ω := rfl
  simp only [SpikedModel.X, SpikedModel.E, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
    hZ, m.stackZW_eq_rowScale w N ω, rowScale_apply]
  ring

/-- **The weighted stack Gram matrix has a simple top eigenvalue almost surely, when every
weight is nonzero.** Section 3 removes the zero weights. -/
theorem topSimple_ae_stackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (hG : m.JointGaussianNoise) (hw : ∀ i, w i ≠ 0) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) := by
  set A : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
    (m.stackW w).θ • Matrix.vecMulVec (WithLp.ofLp ((m.stackW w).u N))
      (WithLp.ofLp ((m.stackW w).v N)) with hA
  set s : Fin (∑ i, n i N) → ℝ :=
    fun r => (Real.sqrt (d N))⁻¹ * w (finSigmaFinEquiv.symm r).1 with hsdef
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast (m.tbl 0).hd N
  have hs : ∀ r, s r ≠ 0 := fun r =>
    mul_ne_zero (inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne') (hw _)
  set P : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ → Prop := fun Y =>
    TopSimple ((A + rowScale s Y)ᵀ * (A + rowScale s Y))
      (isHermitian_transpose_mul_self (A + rowScale s Y)) with hP
  have hPm : Measurable P := by
    rw [← measurableSet_setOfPred]
    have hpre : {Y | P Y} = (fun Y => A + rowScale s Y) ⁻¹'
        {Z : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ |
          TopSimple (Zᵀ * Z) (isHermitian_transpose_mul_self Z)} := rfl
    rw [hpre]
    exact measurable_addRowScale A s measurableSet_topSimple
  have hlaw : HasLaw (m.stackZ N) (gaussianMatrix (∑ i, n i N) (d N)) (μ N) :=
    m.stack_law hG N
  have hae : ∀ᵐ ω ∂(μ N), P (m.stackZ N ω) :=
    (hlaw.ae_iff hPm).mpr
      (topSimple_ae_rowScale (m.stack_row_pos N) ((m.tbl 0).hd N) A hs)
  filter_upwards [hae] with ω hω
  have hX : (m.stackW w).X N ω = A + rowScale s (m.stackZ N ω) :=
    m.stackW_X_eq_addRowScale w N ω
  change TopSimple (((m.stackW w).X N ω)ᵀ * (m.stackW w).X N ω)
    (isHermitian_transpose_mul_self ((m.stackW w).X N ω))
  exact (topSimple_congr_gram ((m.stackW w).X N ω) (A + rowScale s (m.stackZ N ω))
    (by rw [hX])).mpr hω

/-! ### 3. Zero weights: pass to the sub-collection of the nonzero weights -/

/-- The weighted stack of the whole collection and the weighted stack of a sub-collection that
carries every nonzero weight have the same Gram matrix: the discarded tables contribute zero
rows. This is `MultiTableModel.gram_restrict` for general weights. -/
theorem stackGramW_restrict (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (S : Finset (Fin M)) [NeZero S.card] (hS : ∀ i, i ∉ S → w i = 0) (N : ℕ) (ω : Ω N) :
    ((m.stackW w).X N ω)ᵀ * ((m.stackW w).X N ω)
      = (((m.restrict S).stackW fun j => w (S.orderEmbOfFin rfl j)).X N ω)ᵀ
        * (((m.restrict S).stackW fun j => w (S.orderEmbOfFin rfl j)).X N ω) := by
  classical
  ext k l
  rw [m.stackGramW_apply w N ω k l,
    (m.restrict S).stackGramW_apply (fun j => w (S.orderEmbOfFin rfl j)) N ω k l]
  simp only [restrict_tbl]
  rw [sum_orderEmb S
    (fun i => w i ^ 2 * ∑ j, (m.tbl i).X N ω j k * (m.tbl i).X N ω j l)]
  refine (Finset.sum_subset (Finset.subset_univ S) fun i _ hi => ?_).symm
  rw [hS i hi]
  ring

/-- **Item H10.** The `topSimple` field of `MultiTableModel.HeteroLaw` for jointly Gaussian
noise. `hw` is necessary: at `w = 0` the weighted stack is the zero matrix, and its top
eigenvalue has multiplicity `d N`. -/
theorem heteroLaw_topSimple_of_gaussian (m : MultiTableModel μ M n d) (w : Fin M → ℝ)
    (hG : m.JointGaussianNoise) (hw : ∃ i, w i ≠ 0) :
    ∀ N, ∀ᵐ ω ∂(μ N),
      TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) := by
  classical
  intro N
  set S : Finset (Fin M) := {i | w i ≠ 0} with hSdef
  obtain ⟨i₀, hi₀⟩ := hw
  have hmem : ∀ i, i ∈ S ↔ w i ≠ 0 := by
    intro i
    simp [hSdef]
  have hSne : S.Nonempty := ⟨i₀, (hmem i₀).mpr hi₀⟩
  have : NeZero S.card := ⟨(Finset.card_pos.mpr hSne).ne'⟩
  have hw' : ∀ j : Fin S.card, w (S.orderEmbOfFin rfl j) ≠ 0 := fun j =>
    (hmem _).mp (S.orderEmbOfFin_mem rfl j)
  have hout : ∀ i, i ∉ S → w i = 0 := fun i hi => by
    by_contra h
    exact hi ((hmem i).mpr h)
  have h := (m.restrict S).topSimple_ae_stackGramW (fun j => w (S.orderEmbOfFin rfl j))
    (m.restrict_jointGaussianNoise S hG) hw' N
  filter_upwards [h] with ω hω
  change TopSimple (((m.stackW w).X N ω)ᵀ * (m.stackW w).X N ω)
    (isHermitian_transpose_mul_self ((m.stackW w).X N ω))
  exact (topSimple_congr_gram ((m.stackW w).X N ω)
    (((m.restrict S).stackW fun j => w (S.orderEmbOfFin rfl j)).X N ω)
    (m.stackGramW_restrict w S hout N ω)).mpr hω

end MultiTableModel

end StackedSVD
