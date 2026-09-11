/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianInterpolation
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.MultivariateGaussian
public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.LogSumExp
public import StackedSVD.Vendor.COLT83.Mathlib.Matrix.Loewner
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianWidth
public import Mathlib.Data.Matrix.ColumnRowPartitioned
public import Mathlib.Analysis.Calculus.Deriv.MeanValue

set_option autoImplicit false

/-!
# The Sudakov–Fernique inequality

**Sudakov–Fernique inequality** (`sudakov_fernique`): if `X ~ N(0, S)` and `Y ~ N(0, T)` are
centered Gaussian vectors of `ℝ^ι` with `E[(X i - X j)²] ≤ E[(Y i - Y j)²]` for all `i, j`, i.e.
`S i i + S j j - 2 S i j ≤ T i i + T j j - 2 T i j`, then `E[max_i X i] ≤ E[max_i Y i]`.

The proof is the Gaussian interpolation method: for the soft-max `F_β = logSumExp β` and
`Z_t = √t X + √(1 - t) Y` with `X, Y` independent, `ψ(t) = E[F_β(Z_t)]` has derivative
`(1/2) E[∑ᵢⱼ (S - T)ᵢⱼ ∂ᵢⱼF_β(Z_t)]` (`hasDerivAt_integral_comp_gaussianInterp`), which is
nonpositive because `∂ᵢⱼF_β = β (δᵢⱼ pᵢ - pᵢ pⱼ)` with `p` the softmax weights
(`sum_mul_fderiv_fderiv_logSumExp_basisFun`), so that `ψ(1) ≤ ψ(0)`; letting `β → ∞` concludes.
Independent copies of `X` and `Y` are realized as the two blocks `√S (g ∘ inl)`, `√T (g ∘ inr)`
of a standard Gaussian vector `g` of `ℝ^{ι ⊕ ι}` (`blockInl`, `blockInr`).
-/

@[expose] public section

open MeasureTheory Real InnerProductSpace Matrix
open scoped RealInnerProductSpace MatrixOrder

namespace ProbabilityTheory

section monotone

variable {E E' : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [NormedAddCommGroup E']
  [InnerProductSpace ℝ E'] [FiniteDimensional ℝ E'] [MeasurableSpace E'] [BorelSpace E']
  {L₀ L₁ : E' →L[ℝ] E} {F : E → ℝ} {K : ℝ} {κ : Type*} [Fintype κ]

/-- If the Gaussian interpolation derivative is pointwise nonpositive, then
`E[F(L₁ g)] ≤ E[F(L₀ g)]`. -/
lemma integral_comp_le_of_gaussianInterp_nonpos (hF : ContDiff ℝ 2 F)
    (hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ K) (b : OrthonormalBasis κ ℝ E')
    (h : ∀ t ∈ Set.Ioo (0 : ℝ) 1, ∀ g, ∑ k, fderiv ℝ (fderiv ℝ F) (gaussianInterp L₀ L₁ t g)
      (gaussianInterp L₀ L₁ t (b k)) (gaussianInterpDeriv L₀ L₁ t (b k)) ≤ 0) :
    ∫ g, F (L₁ g) ∂stdGaussian E' ≤ ∫ g, F (L₀ g) ∂stdGaussian E' := by
  have hanti : AntitoneOn (fun t ↦ ∫ g, F (gaussianInterp L₀ L₁ t g) ∂stdGaussian E')
      (Set.Icc 0 1) := by
    refine antitoneOn_of_deriv_nonpos (convex_Icc 0 1)
      (continuousOn_integral_comp_gaussianInterp hF hK) ?_ ?_
    · rw [interior_Icc]
      exact fun t ht ↦
        (hasDerivAt_integral_comp_gaussianInterp hF hK b ht).differentiableAt.differentiableWithinAt
    · rw [interior_Icc]
      intro t ht
      rw [(hasDerivAt_integral_comp_gaussianInterp hF hK b ht).deriv]
      exact integral_nonpos fun g ↦ h t ht g
  simpa using hanti (Set.left_mem_Icc.2 zero_le_one) (Set.right_mem_Icc.2 zero_le_one) zero_le_one

end monotone

section blocks

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The linear map `ℝ^{ι ⊕ ι} → ℝ^ι`, `g ↦ A (g ∘ inl)`. -/
noncomputable def blockInl (A : Matrix ι ι ℝ) :
    EuclideanSpace ℝ (ι ⊕ ι) →L[ℝ] EuclideanSpace ℝ ι :=
  LinearMap.toContinuousLinearMap (toEuclideanLin (fromCols A 0))

/-- The linear map `ℝ^{ι ⊕ ι} → ℝ^ι`, `g ↦ B (g ∘ inr)`. -/
noncomputable def blockInr (B : Matrix ι ι ℝ) :
    EuclideanSpace ℝ (ι ⊕ ι) →L[ℝ] EuclideanSpace ℝ ι :=
  LinearMap.toContinuousLinearMap (toEuclideanLin (fromCols 0 B))

lemma blockInl_single_inl (A : Matrix ι ι ℝ) (i : ι) :
    blockInl A (EuclideanSpace.single (Sum.inl i) 1) =
      toEuclideanLin A (EuclideanSpace.single i 1) := by
  ext j
  simp [blockInl, fromCols_apply_inl]

lemma blockInl_single_inr (A : Matrix ι ι ℝ) (i : ι) :
    blockInl A (EuclideanSpace.single (Sum.inr i) 1) = 0 := by
  ext j
  simp [blockInl, fromCols_apply_inr]

lemma blockInr_single_inl (B : Matrix ι ι ℝ) (i : ι) :
    blockInr B (EuclideanSpace.single (Sum.inl i) 1) = 0 := by
  ext j
  simp [blockInr, fromCols_apply_inl]

lemma blockInr_single_inr (B : Matrix ι ι ℝ) (i : ι) :
    blockInr B (EuclideanSpace.single (Sum.inr i) 1) =
      toEuclideanLin B (EuclideanSpace.single i 1) := by
  ext j
  simp [blockInr, fromCols_apply_inr]

/-- The law of `A (g ∘ inl)` for `g ~ N(0, I)` on `ℝ^{ι ⊕ ι}` is `N(0, A Aᵀ)`. -/
lemma map_blockInl_stdGaussian (A : Matrix ι ι ℝ) :
    (stdGaussian (EuclideanSpace ℝ (ι ⊕ ι))).map (blockInl A) =
      multivariateGaussian 0 (A * Aᵀ) := by
  rw [blockInl, ← multivariateGaussian_zero_one,
    multivariateGaussian_zero_map_toEuclideanLin PosSemidef.one, Matrix.mul_one,
    transpose_fromCols, fromCols_mul_fromRows]
  simp

/-- The law of `B (g ∘ inr)` for `g ~ N(0, I)` on `ℝ^{ι ⊕ ι}` is `N(0, B Bᵀ)`. -/
lemma map_blockInr_stdGaussian (B : Matrix ι ι ℝ) :
    (stdGaussian (EuclideanSpace ℝ (ι ⊕ ι))).map (blockInr B) =
      multivariateGaussian 0 (B * Bᵀ) := by
  rw [blockInr, ← multivariateGaussian_zero_one,
    multivariateGaussian_zero_map_toEuclideanLin PosSemidef.one, Matrix.mul_one,
    transpose_fromCols, fromCols_mul_fromRows]
  simp

/-- For a bilinear form `H` on `ℝ^ι` and a matrix `A`,
`∑ᵢ H(A eᵢ, A eᵢ) = ∑ⱼₗ (A Aᵀ)ⱼₗ H(eⱼ, eₗ)`. -/
lemma sum_apply_toEuclideanLin_single_self
    (H : EuclideanSpace ℝ ι →L[ℝ] EuclideanSpace ℝ ι →L[ℝ] ℝ) (A : Matrix ι ι ℝ) :
    ∑ i, H (toEuclideanLin A (EuclideanSpace.single i 1))
        (toEuclideanLin A (EuclideanSpace.single i 1)) =
      ∑ j, ∑ l, (A * Aᵀ) j l * H (EuclideanSpace.single j 1) (EuclideanSpace.single l 1) := by
  have hcol : ∀ i, toEuclideanLin A (EuclideanSpace.single i 1) =
      ∑ j, A j i • (EuclideanSpace.single j 1 : EuclideanSpace ℝ ι) := by
    intro i
    conv_lhs => rw [← (EuclideanSpace.basisFun ι ℝ).sum_repr' (toEuclideanLin A
      (EuclideanSpace.single i 1))]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    simp [EuclideanSpace.basisFun_apply, EuclideanSpace.inner_single_left]
  simp_rw [hcol, map_sum, map_smul, _root_.sum_apply, _root_.smul_apply,
    smul_eq_mul, Finset.mul_sum, mul_apply, transpose_apply, Finset.sum_mul]
  conv_rhs => rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun l _ ↦ ?_
  conv_lhs => rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun i _ ↦ by ring

/-- The Gaussian interpolation derivative for the two blocks `A (g ∘ inl)`, `B (g ∘ inr)`:
`∑ₖ H(Zₜ eₖ, Zₜ' eₖ) = (1/2) ∑ⱼₗ (A Aᵀ)ⱼₗ Hⱼₗ - (1/2) ∑ⱼₗ (B Bᵀ)ⱼₗ Hⱼₗ`. -/
lemma sum_apply_gaussianInterp_block (H : EuclideanSpace ℝ ι →L[ℝ] EuclideanSpace ℝ ι →L[ℝ] ℝ)
    (A B : Matrix ι ι ℝ) {t : ℝ} (ht : t ∈ Set.Ioo (0 : ℝ) 1) :
    ∑ k, H (gaussianInterp (blockInr B) (blockInl A) t (EuclideanSpace.basisFun (ι ⊕ ι) ℝ k))
        (gaussianInterpDeriv (blockInr B) (blockInl A) t (EuclideanSpace.basisFun (ι ⊕ ι) ℝ k)) =
      (1 / 2) * ∑ j, ∑ l, (A * Aᵀ) j l * H (EuclideanSpace.single j 1) (EuclideanSpace.single l 1) -
        (1 / 2) * ∑ j, ∑ l, (B * Bᵀ) j l * H (EuclideanSpace.single j 1)
          (EuclideanSpace.single l 1) := by
  have ht0 : 0 < √t := Real.sqrt_pos.2 ht.1
  have ht1 : 0 < √(1 - t) := Real.sqrt_pos.2 (sub_pos.2 ht.2)
  rw [Fintype.sum_sum_type, ← sum_apply_toEuclideanLin_single_self,
    ← sum_apply_toEuclideanLin_single_self, Finset.mul_sum, Finset.mul_sum,
    ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [EuclideanSpace.basisFun_apply, gaussianInterp_apply, gaussianInterpDeriv_apply,
    blockInl_single_inl, blockInl_single_inr, blockInr_single_inl, blockInr_single_inr, smul_zero,
    add_zero, zero_add, sub_zero, zero_sub, map_smul, map_neg, _root_.smul_apply, smul_eq_mul]
  field_simp
  ring

end blocks

section softmax

variable {ι : Type*} [Fintype ι] [DecidableEq ι] [Nonempty ι] {β : ℝ}

/-- The Hessian of the coordinate soft-max `z ↦ β⁻¹ log ∑ᵢ exp (β zᵢ)` is
`β (δⱼₗ pⱼ - pⱼ pₗ)` with `p` the softmax weights. -/
lemma fderiv_fderiv_logSumExp_basisFun_single (hβ : β ≠ 0) (v : EuclideanSpace ℝ ι) (j l : ι) :
    fderiv ℝ (fderiv ℝ (logSumExp β ⇑(EuclideanSpace.basisFun ι ℝ))) v (EuclideanSpace.single j 1)
        (EuclideanSpace.single l 1) =
      β * ((if j = l then softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v j else 0) -
        softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v j *
          softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v l) := by
  rw [_root_.fderiv_fderiv_logSumExp_apply hβ]
  simp only [EuclideanSpace.basisFun_apply, EuclideanSpace.inner_single_left, conj_trivial,
    PiLp.single_apply, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  by_cases hjl : j = l
  · subst hjl
    simp
  · simp [hjl, Ne.symm hjl]

/-- The contraction of the soft-max Hessian with a symmetric matrix `S`:
`∑ⱼₗ Sⱼₗ ∂ⱼₗF = β (∑ⱼ pⱼ Sⱼⱼ - ∑ⱼₗ pⱼ pₗ Sⱼₗ)`. -/
lemma sum_mul_fderiv_fderiv_logSumExp_basisFun (hβ : β ≠ 0) (v : EuclideanSpace ℝ ι)
    (S : Matrix ι ι ℝ) :
    ∑ j, ∑ l, S j l * fderiv ℝ (fderiv ℝ (logSumExp β ⇑(EuclideanSpace.basisFun ι ℝ))) v
        (EuclideanSpace.single j 1) (EuclideanSpace.single l 1) =
      β * (∑ j, softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v j * S j j -
        ∑ j, ∑ l, softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v j *
          softmax β ⇑(EuclideanSpace.basisFun ι ℝ) v l * S j l) := by
  simp_rw [fderiv_fderiv_logSumExp_basisFun_single hβ, mul_sub, Finset.sum_sub_distrib, mul_ite,
    mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true, Finset.mul_sum]
  congr 1
  · exact Finset.sum_congr rfl fun j _ ↦ by ring
  · exact Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun l _ ↦ by ring

end softmax

section algebra

variable {ι : Type*} [Fintype ι]

/-- For a probability vector `p` and matrices `S T` with `Sᵢᵢ + Sⱼⱼ - 2 Sᵢⱼ ≤ Tᵢᵢ + Tⱼⱼ - 2 Tᵢⱼ`,
`∑ⱼ pⱼ Sⱼⱼ - ∑ⱼₗ pⱼ pₗ Sⱼₗ ≤ ∑ⱼ pⱼ Tⱼⱼ - ∑ⱼₗ pⱼ pₗ Tⱼₗ`. -/
lemma sum_mul_sub_sum_mul_le {p : ι → ℝ} (hp : ∀ i, 0 ≤ p i) (hp1 : ∑ i, p i = 1)
    {S T : Matrix ι ι ℝ} (h : ∀ i j, S i i + S j j - 2 * S i j ≤ T i i + T j j - 2 * T i j) :
    ∑ j, p j * S j j - ∑ j, ∑ l, p j * p l * S j l ≤
      ∑ j, p j * T j j - ∑ j, ∑ l, p j * p l * T j l := by
  have key : ∀ S : Matrix ι ι ℝ, 2 * (∑ j, p j * S j j - ∑ j, ∑ l, p j * p l * S j l) =
      ∑ j, ∑ l, p j * p l * (S j j + S l l - 2 * S j l) := by
    intro S
    have h1 : ∀ j, ∑ l, p j * p l * S j j = p j * S j j := fun j ↦ by
      calc ∑ l, p j * p l * S j j = (∑ l, p l) * (p j * S j j) := by
            rw [Finset.sum_mul]
            exact Finset.sum_congr rfl fun l _ ↦ by ring
        _ = p j * S j j := by rw [hp1, one_mul]
    have h2 : ∑ j, ∑ l, p j * p l * S l l = ∑ j, p j * S j j := by
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun l _ ↦ ?_
      calc ∑ j, p j * p l * S l l = (∑ j, p j) * (p l * S l l) := by
            rw [Finset.sum_mul]
            exact Finset.sum_congr rfl fun j _ ↦ by ring
        _ = p l * S l l := by rw [hp1, one_mul]
    have h3 : ∑ j, ∑ l, p j * p l * (2 * S j l) = 2 * ∑ j, ∑ l, p j * p l * S j l := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j _ ↦ ?_
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun l _ ↦ by ring
    calc 2 * (∑ j, p j * S j j - ∑ j, ∑ l, p j * p l * S j l)
        = ∑ j, ∑ l, p j * p l * S j j + ∑ j, ∑ l, p j * p l * S l l -
          ∑ j, ∑ l, p j * p l * (2 * S j l) := by
          rw [h2, h3]
          simp only [h1]
          ring
      _ = ∑ j, ∑ l, p j * p l * (S j j + S l l - 2 * S j l) := by
          simp only [mul_sub, mul_add, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  have hle := Finset.sum_le_sum fun j (_ : j ∈ Finset.univ) ↦ Finset.sum_le_sum
    fun l (_ : l ∈ Finset.univ) ↦ mul_le_mul_of_nonneg_left (h j l) (mul_nonneg (hp j) (hp l))
  linarith [key S, key T]

end algebra

section main

variable {ι : Type*} [Fintype ι] [DecidableEq ι] [Nonempty ι]

/-- The interpolation derivative of the soft-max between `N(0, S)` and `N(0, T)` is
nonpositive under the Sudakov–Fernique increment condition. -/
lemma sum_fderiv_fderiv_logSumExp_gaussianInterp_nonpos {β : ℝ} (hβ : 0 < β)
    {S T : Matrix ι ι ℝ} (hS : S.PosSemidef) (hT : T.PosSemidef)
    (h : ∀ i j, S i i + S j j - 2 * S i j ≤ T i i + T j j - 2 * T i j) {t : ℝ}
    (ht : t ∈ Set.Ioo (0 : ℝ) 1) (g : EuclideanSpace ℝ (ι ⊕ ι)) :
    ∑ k, fderiv ℝ (fderiv ℝ (logSumExp β ⇑(EuclideanSpace.basisFun ι ℝ)))
      (gaussianInterp (blockInr (CFC.sqrt T)) (blockInl (CFC.sqrt S)) t g)
      (gaussianInterp (blockInr (CFC.sqrt T)) (blockInl (CFC.sqrt S)) t
        (EuclideanSpace.basisFun (ι ⊕ ι) ℝ k))
      (gaussianInterpDeriv (blockInr (CFC.sqrt T)) (blockInl (CFC.sqrt S)) t
        (EuclideanSpace.basisFun (ι ⊕ ι) ℝ k)) ≤ 0 := by
  rw [sum_apply_gaussianInterp_block _ _ _ ht, transpose_sqrt, transpose_sqrt, hS.sqrt_mul_sqrt,
    hT.sqrt_mul_sqrt, sum_mul_fderiv_fderiv_logSumExp_basisFun hβ.ne',
    sum_mul_fderiv_fderiv_logSumExp_basisFun hβ.ne']
  have := sum_mul_sub_sum_mul_le (softmax_nonneg β ⇑(EuclideanSpace.basisFun ι ℝ)
    (gaussianInterp (blockInr (CFC.sqrt T)) (blockInl (CFC.sqrt S)) t g))
    (sum_softmax β ⇑(EuclideanSpace.basisFun ι ℝ)
      (gaussianInterp (blockInr (CFC.sqrt T)) (blockInl (CFC.sqrt S)) t g)) h
  nlinarith

/-- **Sudakov–Fernique inequality**: if `S` and `T` are positive semidefinite matrices with
`S i i + S j j - 2 S i j ≤ T i i + T j j - 2 T i j` for all `i, j` (that is,
`E[(X i - X j)²] ≤ E[(Y i - Y j)²]` for `X ~ N(0, S)` and `Y ~ N(0, T)`), then
`E[max_i X i] ≤ E[max_i Y i]`. -/
theorem sudakov_fernique {S T : Matrix ι ι ℝ} (hS : S.PosSemidef) (hT : T.PosSemidef)
    (h : ∀ i j, S i i + S j j - 2 * S i j ≤ T i i + T j j - 2 * T i j) :
    ∫ x, ⨆ i, x i ∂multivariateGaussian 0 S ≤ ∫ x, ⨆ i, x i ∂multivariateGaussian 0 T := by
  set x : ι → EuclideanSpace ℝ ι := ⇑(EuclideanSpace.basisFun ι ℝ) with hxdef
  have hx : ∀ i, ‖x i‖ ≤ (1 : NNReal) := fun i ↦ by
    simp [hxdef, EuclideanSpace.basisFun_apply, PiLp.norm_single]
  have hsup : ∀ z : EuclideanSpace ℝ ι, (⨆ i, ⟪x i, z⟫) = ⨆ i, z i := fun z ↦ by
    simp [hxdef, EuclideanSpace.basisFun_apply, EuclideanSpace.inner_single_left]
  have hint : ∀ U : Matrix ι ι ℝ,
      Integrable (fun z : EuclideanSpace ℝ ι ↦ ⨆ i, z i) (multivariateGaussian 0 U) := fun U ↦ by
    have hmeas : Measurable fun z : EuclideanSpace ℝ ι ↦ ⨆ i, z i := by
      simpa only [hsup] using measurable_ciSup_inner x
    refine IsGaussian.integrable_of_abs_le_add_mul_norm (A := 0) (B := 1)
      hmeas.aestronglyMeasurable fun z ↦ ?_
    rw [← hsup, zero_add, one_mul]
    simpa using abs_ciSup_inner_le hx z
  refine le_of_forall_pos_le_add fun ε hε ↦ ?_
  have hm : (0 : ℝ) ≤ log (Fintype.card ι) := Real.log_nonneg (by exact_mod_cast Fintype.card_pos)
  set β : ℝ := (log (Fintype.card ι) + 1) / ε with hβdef
  have hβ : 0 < β := by positivity
  set F := logSumExp β x with hFdef
  have hF : ContDiff ℝ 2 F := contDiff_logSumExp' 2 β x
  have hK : ∀ z, ‖fderiv ℝ (fderiv ℝ F) z‖ ≤ β * (1 : NNReal) ^ 2 :=
    norm_fderiv_fderiv_logSumExp_le hβ hx
  have hsand : ∀ z : EuclideanSpace ℝ ι,
      (⨆ i, z i) ≤ F z ∧ F z ≤ (⨆ i, z i) + log (Fintype.card ι) / β := fun z ↦ by
    have h1 := ciSup_inner_le_logSumExp hβ x z
    have h2 := logSumExp_le_ciSup_inner_add hβ x z
    rw [hsup] at h1 h2
    exact ⟨h1, h2⟩
  have hFint : ∀ U : Matrix ι ι ℝ, Integrable F (multivariateGaussian 0 U) := fun U ↦ by
    refine IsGaussian.integrable_of_abs_le_add_mul_norm (A := log (Fintype.card ι) / β) (B := 1)
      hF.continuous.aestronglyMeasurable fun z ↦ ?_
    have h1 := abs_ciSup_inner_le hx z
    rw [hsup] at h1
    rw [abs_le] at h1 ⊢
    simp only [NNReal.coe_one, one_mul] at h1 ⊢
    constructor <;> linarith [(hsand z).1, (hsand z).2, div_nonneg hm hβ.le]
  have hmain : ∫ z, F z ∂multivariateGaussian 0 S ≤ ∫ z, F z ∂multivariateGaussian 0 T := by
    have hSmap : multivariateGaussian 0 S =
        (stdGaussian (EuclideanSpace ℝ (ι ⊕ ι))).map (blockInl (CFC.sqrt S)) := by
      rw [map_blockInl_stdGaussian, transpose_sqrt, hS.sqrt_mul_sqrt]
    have hTmap : multivariateGaussian 0 T =
        (stdGaussian (EuclideanSpace ℝ (ι ⊕ ι))).map (blockInr (CFC.sqrt T)) := by
      rw [map_blockInr_stdGaussian, transpose_sqrt, hT.sqrt_mul_sqrt]
    rw [hSmap, hTmap, integral_map (by fun_prop) hF.continuous.aestronglyMeasurable,
      integral_map (by fun_prop) hF.continuous.aestronglyMeasurable]
    exact integral_comp_le_of_gaussianInterp_nonpos hF hK (EuclideanSpace.basisFun (ι ⊕ ι) ℝ)
      fun t ht g ↦ sum_fderiv_fderiv_logSumExp_gaussianInterp_nonpos hβ hS hT h ht g
  have hlast : log (Fintype.card ι) / β ≤ ε := by
    rw [hβdef, div_div_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith
  calc ∫ z, ⨆ i, z i ∂multivariateGaussian 0 S
      ≤ ∫ z, F z ∂multivariateGaussian 0 S := integral_mono (hint S) (hFint S) fun z ↦ (hsand z).1
    _ ≤ ∫ z, F z ∂multivariateGaussian 0 T := hmain
    _ ≤ ∫ z, ((⨆ i, z i) + log (Fintype.card ι) / β) ∂multivariateGaussian 0 T :=
        integral_mono (hFint T) ((hint T).add (integrable_const _)) fun z ↦ (hsand z).2
    _ = ∫ z, ⨆ i, z i ∂multivariateGaussian 0 T + log (Fintype.card ι) / β := by
        rw [integral_add (hint T) (integrable_const _), integral_const]
        simp
    _ ≤ ∫ z, ⨆ i, z i ∂multivariateGaussian 0 T + ε := by gcongr

end main

end ProbabilityTheory
