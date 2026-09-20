/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RankR.RMT.Stack

/-!
# Task U1: the rank-`r` splitting bridge

Task U1 of `notes/archive/plan_subspacelaw.md` (section 1.1, and the U1 row of section 3). This is
the rank-`r` twin of item R0 (`RMT/R0.lean`), stated model-free first and then instantiated
on `UnalignedModel`.

The setting is one matrix `X = A Vᵀ + t Z` with `A : Fin n × Fin r`, `V : Fin d × Fin r`
with orthonormal columns, `Z` a canonical Gaussian matrix and `t = d^{-1/2}`. With
`C = Aᵀ A = Q_C Λ Q_Cᵀ`:

1. **The factorization** (`exists_orthonormal_factor`). `A Q_C = U Λ^{1/2}` with `Uᵀ U = 1`.
   Column `j` of `U` is `A Q_C e_j / √λ_j` when `λ_j ≠ 0`. When `λ_j = 0` that column of
   `A Q_C` vanishes, and `U` gets an arbitrary orthonormal direction instead (the junk case):
   the identity still holds because both sides carry the factor `√λ_j = 0`. Nothing later
   reads such a column, so the junk direction never enters a limit. The completion is
   `Orthonormal.exists_orthonormalBasis_extension_of_card_eq`, and it needs `r ≤ n`.
2. **The split** (`gram_mul_transpose_add`, `gram_perp_eq`). With `G = Eᵀ U`,
   `E⊥ = E - U Gᵀ = (1 - U Uᵀ) E` and `Q = V Q_C Λ^{1/2} + G`, the Gram matrix is
   `Xᵀ X = W₀ + Q Qᵀ`, `W₀ = E⊥ᵀ E⊥`. Checked numerically before the proof (residual
   `1.4e-14` at `n, d, r = 23, 17, 4`, seed `2026083001`, and again with a zero spike).
3. **The laws** (`exists_frameBlock`). Complete `U` to an orthonormal basis of `ℝⁿ`, rotate
   `Z` (left invariance), and split the first `r` rows off the last `n - r`. The pair
   `(Uᵀ Z, B)` has the product law `gaussianMatrix r d ⊗ gaussianMatrix (n - r) d`, and
   `W₀ = d⁻¹ Bᵀ B`. `exists_frameBlock'` transposes the first component to the `d × r` form
   `Zᵀ U` of the plan.
4. **The aspect ratio** (`tendsto_sub_ratio`). `(n - r)/d → c` when `n/d → c`, `d → ∞`.
5. **The model** (`UnalignedModel.exists_rankR_split`). The stack `X_stack = A Vᵀ + E_stack`
   of `stackX_eq` at `A = signalFactor`, with the stacked noise law
   (`hasLaw_stackZu`) obtained from `MultiTableModel.stack_law` through a shared-spike copy
   of the tables (the noise fields are untouched, so the law transfers).

Paper: `main_paper.tex` lines 815 to 843 (`prop:stacksvd_subspace`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 1. Deterministic algebra of a rank-`r` split -/

section GramAlgebra

variable {n r c : ℕ}

/-- **The rank-`r` Gram identity.** If `F` kills the frame `U` on the left then the Gram
matrix of `U Qᵀ + F` splits into the Gram matrix of `F` and `Q Qᵀ`. This is the rank-`r`
form of `gram_vecMulVec_add` (`RMT/R0.lean`), plan section 1.1 item 2. -/
theorem gram_mul_transpose_add (U : Matrix (Fin n) (Fin r) ℝ) (Q : Matrix (Fin c) (Fin r) ℝ)
    (F : Matrix (Fin n) (Fin c) ℝ) (hU : Uᵀ * U = 1) (hF : Fᵀ * U = 0) :
    (U * Qᵀ + F)ᵀ * (U * Qᵀ + F) = Fᵀ * F + Q * Qᵀ := by
  have hF' : Uᵀ * F = 0 := by
    have := congrArg Matrix.transpose hF
    rwa [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.transpose_zero] at this
  rw [Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_transpose,
    Matrix.add_mul, Matrix.mul_add, Matrix.mul_add]
  rw [Matrix.mul_assoc Q Uᵀ (U * Qᵀ), ← Matrix.mul_assoc Uᵀ U Qᵀ, hU, Matrix.one_mul,
    Matrix.mul_assoc Q Uᵀ F, hF', Matrix.mul_zero, ← Matrix.mul_assoc Fᵀ U Qᵀ, hF,
    Matrix.zero_mul]
  abel

/-- The frame complement `Z - U Uᵀ Z` is orthogonal to the frame. Rank-`r` form of
`transpose_sub_vecMulVec_mulVec` (`RMT/R0.lean`). -/
theorem transpose_perp_mul (U : Matrix (Fin n) (Fin r) ℝ) (Z : Matrix (Fin n) (Fin c) ℝ)
    (hU : Uᵀ * U = 1) : (Z - U * (Uᵀ * Z))ᵀ * U = 0 := by
  rw [Matrix.transpose_sub, Matrix.transpose_mul, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.sub_mul, Matrix.mul_assoc, hU, Matrix.mul_one,
    sub_self]

/-- **The Gram matrix of the frame complement.** `(P_U⊥ Z)ᵀ (P_U⊥ Z) = Zᵀ Z - (Zᵀ U)(Zᵀ U)ᵀ`.
Rank-`r` form of the identity inside `SpikedModel.W0_eq_smul`. -/
theorem gram_perp_eq (U : Matrix (Fin n) (Fin r) ℝ) (Z : Matrix (Fin n) (Fin c) ℝ)
    (hU : Uᵀ * U = 1) :
    (Z - U * (Uᵀ * Z))ᵀ * (Z - U * (Uᵀ * Z)) = Zᵀ * Z - (Zᵀ * U) * (Zᵀ * U)ᵀ := by
  set F : Matrix (Fin n) (Fin c) ℝ := Z - U * (Uᵀ * Z) with hFdef
  have hF : Fᵀ * U = 0 := transpose_perp_mul U Z hU
  have hZ : U * (Zᵀ * U)ᵀ + F = Z := by
    rw [hFdef, Matrix.transpose_mul, Matrix.transpose_transpose]
    abel
  have key := gram_mul_transpose_add U (Zᵀ * U) F hU hF
  rw [hZ] at key
  rw [key]
  abel

end GramAlgebra

/-! ### 2. The factorization `A Q_C = U Λ^{1/2}` -/

section Factor

variable {n r : ℕ}

/-- The columns of `A Q_C` are pairwise orthogonal with squared norms the eigenvalues of
`C = Aᵀ A`. Plan section 1.1 item 1, first half. -/
theorem gram_mul_eigenvectorMatrix (A : Matrix (Fin n) (Fin r) ℝ)
    (Qc : Matrix (Fin r) (Fin r) ℝ) (lam : Fin r → ℝ) (hQ : Qcᵀ * Qc = 1)
    (hA : Aᵀ * A = Qc * Matrix.diagonal lam * Qcᵀ) :
    (A * Qc)ᵀ * (A * Qc) = Matrix.diagonal lam := by
  have h1 : Qcᵀ * (Qc * Matrix.diagonal lam * Qcᵀ * Qc) = Matrix.diagonal lam := by
    rw [Matrix.mul_assoc (Qc * Matrix.diagonal lam) Qcᵀ Qc, hQ, Matrix.mul_one,
      ← Matrix.mul_assoc, hQ, Matrix.one_mul]
  rw [Matrix.transpose_mul, Matrix.mul_assoc Qcᵀ Aᵀ (A * Qc), ← Matrix.mul_assoc Aᵀ A Qc, hA,
    h1]

/-- **The factorization of the signal factor**, plan section 1.1 item 1:
`A Q_C = U Λ^{1/2}` with `U` an orthonormal `r`-frame of `ℝⁿ`.

Column `j` of `U` is `A Q_C e_j / √λ_j` when `λ_j ≠ 0`. When `λ_j = 0` that column of
`A Q_C` is zero and `U` gets an arbitrary orthonormal direction from the basis extension
(the junk case); the identity survives because the right side carries the factor
`√λ_j = 0`. The side condition `r ≤ n` is the one the completion needs. -/
theorem exists_orthonormal_factor (hrn : r ≤ n) (A : Matrix (Fin n) (Fin r) ℝ)
    (Qc : Matrix (Fin r) (Fin r) ℝ) (lam : Fin r → ℝ) (hQ : Qcᵀ * Qc = 1)
    (hA : Aᵀ * A = Qc * Matrix.diagonal lam * Qcᵀ) :
    ∃ U : Matrix (Fin n) (Fin r) ℝ, Uᵀ * U = 1 ∧
      A * Qc = U * Matrix.diagonal fun j => Real.sqrt (lam j) := by
  classical
  set W : Matrix (Fin n) (Fin r) ℝ := A * Qc with hWdef
  have hWW : Wᵀ * W = Matrix.diagonal lam := gram_mul_eigenvectorMatrix A Qc lam hQ hA
  have hdot : ∀ j j' : Fin r, ∑ a : Fin n, W a j * W a j'
      = if j = j' then lam j else 0 := by
    intro j j'
    have h := congrFun (congrFun hWW j) j'
    rw [Matrix.mul_apply] at h
    simp only [Matrix.transpose_apply, Matrix.diagonal_apply] at h
    rw [h]
  have hcol : ∀ j : Fin r, ∑ a : Fin n, W a j * W a j = lam j := by
    intro j; simpa using hdot j j
  have hnonneg : ∀ j : Fin r, 0 ≤ lam j := by
    intro j
    rw [← hcol j]
    exact Finset.sum_nonneg fun a _ => mul_self_nonneg _
  -- the column of `W`, as a vector of `EuclideanSpace`
  set col : Fin r → EuclideanSpace ℝ (Fin n) :=
    fun j => WithLp.toLp 2 fun a => W a j with hcoldef
  have hinner : ∀ j j' : Fin r, ⟪col j, col j'⟫_ℝ = if j = j' then lam j else 0 := by
    intro j j'
    rw [inner_euclidean_eq_dotProduct, ← hdot j j']
    rfl
  have hnormcol : ∀ j : Fin r, ‖col j‖ = Real.sqrt (lam j) := by
    intro j
    have h1 : ‖col j‖ ^ 2 = lam j := by
      rw [← real_inner_self_eq_norm_sq, hinner j j, if_pos rfl]
    rw [← h1, Real.sqrt_sq (norm_nonneg _)]
  -- the prescribed family and the set where it is prescribed
  set emb : Fin r → Fin n := Fin.castLE hrn with hembdef
  have hembval : ∀ j : Fin r, ((emb j : Fin n) : ℕ) = (j : ℕ) := fun j => rfl
  have hembinj : ∀ j j' : Fin r, emb j = emb j' → j = j' := by
    intro j j' h
    exact Fin.ext (by rw [← hembval j, ← hembval j', h])
  set v : Fin n → EuclideanSpace ℝ (Fin n) := fun i =>
    if h : (i : ℕ) < r then (Real.sqrt (lam ⟨i, h⟩))⁻¹ • col ⟨i, h⟩ else 0 with hvdef
  set s : Set (Fin n) := {i : Fin n | ∃ h : (i : ℕ) < r, lam ⟨i, h⟩ ≠ 0} with hsdef
  have hv_pos : ∀ (i : Fin n) (h : (i : ℕ) < r),
      v i = (Real.sqrt (lam ⟨i, h⟩))⁻¹ • col ⟨i, h⟩ := by
    intro i h
    simp only [hvdef]
    rw [dif_pos h]
  have hembmem : ∀ j : Fin r, lam j ≠ 0 → emb j ∈ s := by
    intro j hj
    refine ⟨by rw [hembval j]; exact j.isLt, ?_⟩
    have : (⟨((emb j : Fin n) : ℕ), by rw [hembval j]; exact j.isLt⟩ : Fin r) = j :=
      Fin.ext (hembval j)
    rw [this]; exact hj
  have hvemb : ∀ j : Fin r, lam j ≠ 0 →
      v (emb j) = (Real.sqrt (lam j))⁻¹ • col j := by
    intro j _
    have hlt : ((emb j : Fin n) : ℕ) < r := by rw [hembval j]; exact j.isLt
    have hj : (⟨((emb j : Fin n) : ℕ), hlt⟩ : Fin r) = j := Fin.ext (hembval j)
    rw [hv_pos (emb j) hlt, hj]
  -- the prescribed family is orthonormal on `s`
  have hortho : Orthonormal ℝ (s.domRestrict v) := by
    constructor
    · rintro ⟨i, h, hne⟩
      change ‖v i‖ = 1
      rw [hv_pos i h, norm_smul, hnormcol ⟨i, h⟩]
      have hpos : 0 < lam ⟨i, h⟩ := lt_of_le_of_ne (hnonneg _) (Ne.symm hne)
      have hs : Real.sqrt (lam ⟨i, h⟩) ≠ 0 := Real.sqrt_ne_zero'.2 hpos
      rw [norm_inv, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _)]
      field_simp
    · rintro ⟨i, h, hne⟩ ⟨i', h', hne'⟩ hij
      change ⟪v i, v i'⟫_ℝ = 0
      have hne2 : (⟨i, h⟩ : Fin r) ≠ ⟨i', h'⟩ := by
        intro hcon
        have hval : (i : ℕ) = (i' : ℕ) := congrArg (Fin.val (n := r)) hcon
        exact hij (Subtype.ext (Fin.ext hval))
      rw [hv_pos i h, hv_pos i' h', real_inner_smul_left, real_inner_smul_right,
        hinner ⟨i, h⟩ ⟨i', h'⟩, if_neg hne2, mul_zero, mul_zero]
  have hcard : Module.finrank ℝ (EuclideanSpace ℝ (Fin n)) = Fintype.card (Fin n) := by
    rw [finrank_euclideanSpace_fin, Fintype.card_fin]
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq hcard
  refine ⟨Matrix.of fun a j => WithLp.ofLp (b (emb j)) a, ?_, ?_⟩
  · ext j j'
    have hb' := (orthonormal_iff_ite (𝕜 := ℝ)).1 b.orthonormal (emb j) (emb j')
    rw [inner_euclidean_eq_dotProduct] at hb'
    have hL : ((Matrix.of fun a j => WithLp.ofLp (b (emb j)) a)ᵀ *
        Matrix.of fun a j => WithLp.ofLp (b (emb j)) a) j j'
        = WithLp.ofLp (b (emb j)) ⬝ᵥ WithLp.ofLp (b (emb j')) := rfl
    rw [hL, hb', Matrix.one_apply]
    by_cases hjj : j = j'
    · rw [if_pos (by rw [hjj]), if_pos hjj]
    · rw [if_neg fun hc => hjj (hembinj j j' hc), if_neg hjj]
  · ext a j
    rw [Matrix.mul_diagonal]
    by_cases hj : lam j = 0
    · have hz : ∑ x : Fin n, W x j * W x j = 0 := by rw [hcol j, hj]
      have hall : ∀ x ∈ (Finset.univ : Finset (Fin n)), W x j * W x j = 0 :=
        (Finset.sum_eq_zero_iff_of_nonneg fun x _ => mul_self_nonneg _).1 hz
      have : W a j = 0 := by
        have := hall a (Finset.mem_univ a)
        exact mul_self_eq_zero.1 this
      rw [this, hj, Real.sqrt_zero, mul_zero]
    · have hpos : 0 < lam j := lt_of_le_of_ne (hnonneg j) (Ne.symm hj)
      have hs : Real.sqrt (lam j) ≠ 0 := Real.sqrt_ne_zero'.2 hpos
      have hbj : b (emb j) = (Real.sqrt (lam j))⁻¹ • col j := by
        rw [hb (emb j) (hembmem j hj), hvemb j hj]
      change W a j = WithLp.ofLp (b (emb j)) a * Real.sqrt (lam j)
      rw [hbj]
      change W a j = (Real.sqrt (lam j))⁻¹ * W a j * Real.sqrt (lam j)
      field_simp

end Factor

/-! ### 3. The `r`-row split and the product law -/

section RowSplit

/-- Reindexing a product of copies of one measure preserves it. This is
`measurePreserving_piCongrEquiv` (`Prob/GaussianMatrix.lean`) with the factor `ℝ` replaced by
a general factor; here the factor is one row `Fin d → ℝ` of the matrix. -/
private theorem measurePreserving_reindexPi {ι κ : Type*} [Fintype ι] [Fintype κ]
    {β : Type*} [MeasurableSpace β] (ν : Measure β) [SigmaFinite ν] (e : κ ≃ ι) :
    MeasurePreserving (fun x : ι → β => fun k => x (e k))
      (Measure.pi fun _ : ι => ν) (Measure.pi fun _ : κ => ν) := by
  have hmeas : Measurable (fun x : ι → β => fun k => x (e k)) :=
    measurable_pi_lambda _ fun k => measurable_pi_apply (e k)
  refine ⟨hmeas, ?_⟩
  refine (Measure.pi_eq fun t ht => ?_).symm
  rw [Measure.map_apply hmeas (MeasurableSet.univ_pi ht)]
  have hpre : (fun x : ι → β => fun k => x (e k)) ⁻¹' Set.univ.pi t
      = Set.univ.pi fun i => t (e.symm i) := by
    ext x
    simp only [Set.mem_preimage, Set.mem_univ_pi]
    constructor
    · intro h i
      have hi := h (e.symm i)
      rwa [Equiv.apply_symm_apply] at hi
    · intro h k
      have hk := h (e k)
      rwa [Equiv.symm_apply_apply] at hk
  rw [hpre, Measure.pi_pi]
  exact Equiv.prod_comp e.symm fun k => ν (t k)

/-- The `r`-row split on plain pi types. `Matrix` is a `def`, so instance search does not
find the product sigma-algebra on it; the `Matrix` form below is this statement, transferred
by `exact` (the device of `measurePreserving_rowSplit` in `RMT/R0.lean`). -/
private theorem measurePreserving_rowSplitAddPi (rr p d : ℕ) :
    MeasurePreserving
      (fun Z : Fin (rr + p) → Fin d → ℝ =>
        ((fun i : Fin rr => Z (Fin.castAdd p i)), fun i : Fin p => Z (Fin.natAdd rr i)))
      (Measure.pi fun _ : Fin (rr + p) => Measure.pi fun _ : Fin d => gaussianReal 0 1)
      ((Measure.pi fun _ : Fin rr => Measure.pi fun _ : Fin d => gaussianReal 0 1).prod
        (Measure.pi fun _ : Fin p => Measure.pi fun _ : Fin d => gaussianReal 0 1)) := by
  have h1 := measurePreserving_reindexPi
    (Measure.pi fun _ : Fin d => gaussianReal 0 1) (finSumFinEquiv (m := rr) (n := p))
  have h2 := measurePreserving_sumPiEquivProdPi
    (fun _ : Fin rr ⊕ Fin p => Measure.pi fun _ : Fin d => gaussianReal 0 1)
  exact h2.comp h1

/-- **Splitting the first `r` rows off a canonical Gaussian matrix gives a product law.**
The rank-`r` form of `measurePreserving_rowSplit` (`RMT/R0.lean`), through
`Fin r ⊕ Fin p ≃ Fin (r + p)`. -/
theorem measurePreserving_rowSplitAdd (rr p d : ℕ) :
    MeasurePreserving
      (fun Z : Matrix (Fin (rr + p)) (Fin d) ℝ =>
        (Matrix.of fun (i : Fin rr) (j : Fin d) => Z (Fin.castAdd p i) j,
          Matrix.of fun (i : Fin p) (j : Fin d) => Z (Fin.natAdd rr i) j))
      (gaussianMatrix (rr + p) d) ((gaussianMatrix rr d).prod (gaussianMatrix p d)) :=
  measurePreserving_rowSplitAddPi rr p d

/-- **The rank-`r` split of the noise**, plan section 1.1 item 3, model-free.

For an orthonormal `r`-frame `U` of `ℝⁿ` with `n = r + p` there is a map `Yfun` on `n × d`
matrices (the last `p` rows of `O Z`, with `O` orthogonal and first `r` rows `Uᵀ`) whose Gram
matrix is the Gram matrix of `Z` with the `U` directions removed, and the pair `(Uᵀ Z, Yfun Z)`
has the product Gaussian law. This is `exists_rowBlock_measurePreserving` (`RMT/Het/Split.lean`)
at an `r`-frame in place of one unit vector, and `SpikedModel.exists_block_hasLaw`
(`RMT/R0.lean`) with the same six steps. -/
theorem exists_frameBlock {nn p rr d : ℕ} (hn : nn = rr + p)
    (U : Matrix (Fin nn) (Fin rr) ℝ) (hU : Uᵀ * U = 1) :
    ∃ Yfun : Matrix (Fin nn) (Fin d) ℝ → Matrix (Fin p) (Fin d) ℝ,
      (∀ Z : Matrix (Fin nn) (Fin d) ℝ,
          (Yfun Z)ᵀ * Yfun Z = Zᵀ * Z - (Zᵀ * U) * (Zᵀ * U)ᵀ) ∧
        MeasurePreserving (fun Z : Matrix (Fin nn) (Fin d) ℝ => (Uᵀ * Z, Yfun Z))
          (gaussianMatrix nn d) ((gaussianMatrix rr d).prod (gaussianMatrix p d)) := by
  classical
  -- 1. an orthonormal basis of `ℝⁿ`, indexed by `Fin (rr + p)`, extending the frame.
  set ucol : Fin rr → EuclideanSpace ℝ (Fin nn) :=
    fun j => WithLp.toLp 2 fun a => U a j with hucoldef
  have hucol : ∀ j j' : Fin rr, ⟪ucol j, ucol j'⟫_ℝ = if j = j' then 1 else 0 := by
    intro j j'
    have h := congrFun (congrFun hU j) j'
    rw [Matrix.mul_apply] at h
    simp only [Matrix.transpose_apply, Matrix.one_apply] at h
    rw [inner_euclidean_eq_dotProduct, ← h]
    rfl
  set vf : Fin (rr + p) → EuclideanSpace ℝ (Fin nn) :=
    fun i => if h : (i : ℕ) < rr then ucol ⟨i, h⟩ else 0 with hvfdef
  set sset : Set (Fin (rr + p)) := {i : Fin (rr + p) | (i : ℕ) < rr} with hssetdef
  have hvf_pos : ∀ (i : Fin (rr + p)) (h : (i : ℕ) < rr), vf i = ucol ⟨i, h⟩ := by
    intro i h
    simp only [hvfdef]
    rw [dif_pos h]
  have hcastval : ∀ j : Fin rr, ((Fin.castAdd p j : Fin (rr + p)) : ℕ) = (j : ℕ) := fun _ => rfl
  have hvfcast : ∀ j : Fin rr, vf (Fin.castAdd p j) = ucol j := by
    intro j
    have hlt : ((Fin.castAdd p j : Fin (rr + p)) : ℕ) < rr := by rw [hcastval j]; exact j.isLt
    rw [hvf_pos _ hlt]
    congr 1
  have hortho : Orthonormal ℝ (sset.domRestrict vf) := by
    constructor
    · rintro ⟨i, h⟩
      change ‖vf i‖ = 1
      have hsq : ‖ucol (⟨i, h⟩ : Fin rr)‖ ^ 2 = 1 := by
        rw [← real_inner_self_eq_norm_sq, hucol ⟨i, h⟩ ⟨i, h⟩, if_pos rfl]
      rw [hvf_pos i h, ← Real.sqrt_sq (norm_nonneg (ucol (⟨i, h⟩ : Fin rr))), hsq,
        Real.sqrt_one]
    · rintro ⟨i, h⟩ ⟨i', h'⟩ hij
      change ⟪vf i, vf i'⟫_ℝ = 0
      have hne2 : (⟨i, h⟩ : Fin rr) ≠ ⟨i', h'⟩ := by
        intro hcon
        have hval : (i : ℕ) = (i' : ℕ) := congrArg (Fin.val (n := rr)) hcon
        exact hij (Subtype.ext (Fin.ext hval))
      rw [hvf_pos i h, hvf_pos i' h', hucol ⟨i, h⟩ ⟨i', h'⟩, if_neg hne2]
  have hcard : Module.finrank ℝ (EuclideanSpace ℝ (Fin nn)) = Fintype.card (Fin (rr + p)) := by
    rw [finrank_euclideanSpace_fin, Fintype.card_fin, hn]
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq hcard
  have hcastmem : ∀ i : Fin rr, (Fin.castAdd p i : Fin (rr + p)) ∈ sset := by
    intro i
    change ((Fin.castAdd p i : Fin (rr + p)) : ℕ) < rr
    rw [hcastval i]; exact i.isLt
  have hbcast : ∀ i : Fin rr, b (Fin.castAdd p i) = ucol i := by
    intro i
    rw [hb _ (hcastmem i), hvfcast i]
  -- 2. the rotation matrix: row `i` is `b i`, read through `Fin (rr + p) ≃ Fin nn`.
  set e : Fin (rr + p) → Fin nn := Fin.cast hn.symm with hedef
  set O : Matrix (Fin (rr + p)) (Fin (rr + p)) ℝ :=
    Matrix.of fun i i' => WithLp.ofLp (b i) (e i') with hOdef
  have hsum : ∀ f : Fin nn → ℝ, ∑ i' : Fin (rr + p), f (e i') = ∑ k : Fin nn, f k := by
    intro f
    exact Equiv.sum_comp (finCongr hn.symm) f
  have hOOt : O * Oᵀ = 1 := by
    ext i i'
    have hb' := (orthonormal_iff_ite (𝕜 := ℝ)).1 b.orthonormal i i'
    rw [inner_euclidean_eq_dotProduct] at hb'
    rw [Matrix.one_apply, ← hb']
    change ∑ i'' : Fin (rr + p),
        WithLp.ofLp (b i) (e i'') * WithLp.ofLp (b i') (e i'') = _
    exact hsum fun k => WithLp.ofLp (b i) k * WithLp.ofLp (b i') k
  have hO : Oᵀ * O = 1 := mul_eq_one_comm.1 hOOt
  -- 3. the composite map `Z ↦ (first r rows of O Z, last p rows of O Z)`.
  set castR : Matrix (Fin nn) (Fin d) ℝ → Matrix (Fin (rr + p)) (Fin d) ℝ :=
    fun Z => Matrix.of fun (i : Fin (rr + p)) (j : Fin d) => Z (e i) j with hcastdef
  set Y : Matrix (Fin nn) (Fin d) ℝ → Matrix (Fin (rr + p)) (Fin d) ℝ :=
    fun Z => O * castR Z with hYdef
  have hmp : MeasurePreserving (fun Z : Matrix (Fin nn) (Fin d) ℝ =>
      (Matrix.of fun (i : Fin rr) (j : Fin d) => Y Z (Fin.castAdd p i) j,
        Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j))
      (gaussianMatrix nn d) ((gaussianMatrix rr d).prod (gaussianMatrix p d)) :=
    (measurePreserving_rowSplitAdd rr p d).comp
      ((measurePreserving_mul_left hO d).comp (measurePreserving_castRow hn))
  -- 4. the rows of `O Z` in terms of `b`, and the first `r` of them as `Uᵀ Z`.
  have hrow : ∀ (Z : Matrix (Fin nn) (Fin d) ℝ) (i : Fin (rr + p)) (j : Fin d),
      Y Z i j = ∑ k : Fin nn, WithLp.ofLp (b i) k * Z k j := by
    intro Z i j
    change ∑ i' : Fin (rr + p), WithLp.ofLp (b i) (e i') * Z (e i') j = _
    exact hsum fun k => WithLp.ofLp (b i) k * Z k j
  have hrowU : ∀ Z : Matrix (Fin nn) (Fin d) ℝ,
      (Matrix.of fun (i : Fin rr) (j : Fin d) => Y Z (Fin.castAdd p i) j) = Uᵀ * Z := by
    intro Z
    ext i j
    change Y Z (Fin.castAdd p i) j = (Uᵀ * Z) i j
    rw [hrow Z (Fin.castAdd p i) j, hbcast i, Matrix.mul_apply]
    exact Finset.sum_congr rfl fun k _ => rfl
  have hent : ∀ (Z : Matrix (Fin nn) (Fin d) ℝ) (i : Fin rr) (j : Fin d),
      Y Z (Fin.castAdd p i) j = (Uᵀ * Z) i j := by
    intro Z i j
    exact congrFun (congrFun (hrowU Z) i) j
  -- 5. the Gram matrix of `O Z` splits into the first `r` rows and the block.
  have hYY : ∀ Z : Matrix (Fin nn) (Fin d) ℝ, (Y Z)ᵀ * Y Z = Zᵀ * Z := by
    intro Z
    rw [hYdef]
    simp only
    rw [Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Oᵀ O (castR Z), hO,
      Matrix.one_mul]
    exact gram_castRow hn Z
  have hgramY : ∀ Z : Matrix (Fin nn) (Fin d) ℝ,
      (Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j)ᵀ *
          (Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j)
        = Zᵀ * Z - (Zᵀ * U) * (Zᵀ * U)ᵀ := by
    intro Z
    ext j j'
    have hsplit : ∑ i : Fin (rr + p), Y Z i j * Y Z i j'
        = (∑ i : Fin rr, Y Z (Fin.castAdd p i) j * Y Z (Fin.castAdd p i) j')
          + ∑ i : Fin p, Y Z (Fin.natAdd rr i) j * Y Z (Fin.natAdd rr i) j' :=
      Fin.sum_univ_add _
    have hL : ((Y Z)ᵀ * Y Z) j j' = ∑ i : Fin (rr + p), Y Z i j * Y Z i j' := rfl
    have hR : (Zᵀ * Z) j j' = ∑ i : Fin (rr + p), Y Z i j * Y Z i j' := by
      rw [← hL, hYY Z]
    have hmat : (Zᵀ * U) * (Zᵀ * U)ᵀ = (Uᵀ * Z)ᵀ * (Uᵀ * Z) := by
      rw [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.transpose_mul,
        Matrix.transpose_transpose]
    have htop : ∑ i : Fin rr, Y Z (Fin.castAdd p i) j * Y Z (Fin.castAdd p i) j'
        = ((Zᵀ * U) * (Zᵀ * U)ᵀ) j j' := by
      rw [hmat]
      have hprod : ((Uᵀ * Z)ᵀ * (Uᵀ * Z)) j j'
          = ∑ i : Fin rr, (Uᵀ * Z) i j * (Uᵀ * Z) i j' := rfl
      rw [hprod]
      exact Finset.sum_congr rfl fun i _ => by rw [hent Z i j, hent Z i j']
    change ∑ i : Fin p, Y Z (Fin.natAdd rr i) j * Y Z (Fin.natAdd rr i) j' = (Zᵀ * Z) j j' - _
    rw [hR, hsplit, htop]
    ring
  -- 6. the block, and the two conclusions.
  refine ⟨fun Z => Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j,
    hgramY, ?_⟩
  have hfun : (fun Z : Matrix (Fin nn) (Fin d) ℝ =>
      (Matrix.of fun (i : Fin rr) (j : Fin d) => Y Z (Fin.castAdd p i) j,
        Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j))
      = fun Z : Matrix (Fin nn) (Fin d) ℝ =>
        (Uᵀ * Z, Matrix.of fun (i : Fin p) (j : Fin d) => Y Z (Fin.natAdd rr i) j) := by
    funext Z
    rw [hrowU Z]
  rwa [hfun] at hmp

/-- The same split with the first component transposed to the `d × r` form `Zᵀ U` of the
plan (section 1.1 item 2, `G = Eᵀ U`). -/
theorem exists_frameBlock' {nn p rr d : ℕ} (hn : nn = rr + p)
    (U : Matrix (Fin nn) (Fin rr) ℝ) (hU : Uᵀ * U = 1) :
    ∃ Yfun : Matrix (Fin nn) (Fin d) ℝ → Matrix (Fin p) (Fin d) ℝ,
      (∀ Z : Matrix (Fin nn) (Fin d) ℝ,
          (Yfun Z)ᵀ * Yfun Z = Zᵀ * Z - (Zᵀ * U) * (Zᵀ * U)ᵀ) ∧
        MeasurePreserving (fun Z : Matrix (Fin nn) (Fin d) ℝ => (Zᵀ * U, Yfun Z))
          (gaussianMatrix nn d) ((gaussianMatrix d rr).prod (gaussianMatrix p d)) := by
  obtain ⟨Yfun, hgram, hmp⟩ := exists_frameBlock (d := d) hn U hU
  refine ⟨Yfun, hgram, ?_⟩
  have hT : MeasurePreserving
      (Prod.map (Matrix.transpose : Matrix (Fin rr) (Fin d) ℝ → Matrix (Fin d) (Fin rr) ℝ)
        (id : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin p) (Fin d) ℝ))
      ((gaussianMatrix rr d).prod (gaussianMatrix p d))
      ((gaussianMatrix d rr).prod (gaussianMatrix p d)) :=
    (measurePreserving_transpose rr d).prod (MeasurePreserving.id _)
  have hcomp := hT.comp hmp
  have hfun : (Prod.map (Matrix.transpose : Matrix (Fin rr) (Fin d) ℝ → _) id) ∘
      (fun Z : Matrix (Fin nn) (Fin d) ℝ => (Uᵀ * Z, Yfun Z))
      = fun Z : Matrix (Fin nn) (Fin d) ℝ => (Zᵀ * U, Yfun Z) := by
    funext Z
    change ((Uᵀ * Z)ᵀ, Yfun Z) = _
    rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  rwa [hfun] at hcomp

end RowSplit

/-! ### 4. The aspect ratio of the block -/

/-- **(H0'') at rank `r`.** `(n - r)/d → c` when `n/d → c` and `d → ∞`. Rank-`r` form of
`SpikedModel.tendsto_pred_ratio` (`RMT/R0.lean`). -/
theorem tendsto_sub_ratio {nf df : ℕ → ℕ} {c : ℝ} (rr : ℕ) (hr : ∀ N, rr ≤ nf N)
    (hd : Tendsto df atTop atTop)
    (hnd : Tendsto (fun N => (nf N : ℝ) / (df N : ℝ)) atTop (𝓝 c)) :
    Tendsto (fun N => ((nf N - rr : ℕ) : ℝ) / (df N : ℝ)) atTop (𝓝 c) := by
  have hdR : Tendsto (fun N => ((df N : ℝ))) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hd
  have hinv0 : Tendsto (fun N => ((df N : ℝ))⁻¹) atTop (𝓝 0) := hdR.inv_tendsto_atTop
  have hinv : Tendsto (fun N => (rr : ℝ) / (df N : ℝ)) atTop (𝓝 0) := by
    have := hinv0.const_mul ((rr : ℝ))
    simpa [div_eq_mul_inv] using this
  have heq : ∀ N : ℕ, ((nf N - rr : ℕ) : ℝ) / (df N : ℝ)
      = (nf N : ℝ) / (df N : ℝ) - (rr : ℝ) / (df N : ℝ) := by
    intro N
    rw [Nat.cast_sub (hr N), sub_div]
  simp only [heq]
  simpa using hnd.sub hinv

/-! ### 4b. The shared stack interface

`RankRStack` (`RankR/RMT/Stack.lean`) carries the data that the chain reads: the signal
factor `A`, the shared `V`, the `N`-free core `C = AᵀA`, the noise `Zu` and its scaling.
The split below is stated on that interface, so a general-`r_i` stack instantiates it. The
`UnalignedModel` theorem of section 5 is the same statement at `UnalignedModel.toStack`. -/

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The Gram matrix of the table splits as `W₀ + Q Qᵀ` (plan 1.1 item 2). -/
theorem gram_eq_rankRW0_add (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} (hU : Uᵀ * U = 1)
    (hsig : s.signalPart N
      = U * (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j))ᵀ) :
    s.gram N ω = s.rankRW0 N ω U + s.qmatR N ω U * (s.qmatR N ω U)ᵀ := by
  have hFperp : (s.stackEperp N ω U)ᵀ * U = 0 := transpose_perp_mul U (s.E N ω) hU
  rw [gram, s.X_eq_frame_add N ω hsig,
    gram_mul_transpose_add U (s.qmatR N ω U) (s.stackEperp N ω U) hU hFperp]
  rfl

/-- On the split the block identity `W₀ = d⁻¹ Bᵀ B` holds pointwise. -/
theorem rankRW0_eq_smul_block (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} (hU : Uᵀ * U = 1)
    (hsig : s.signalPart N
      = U * (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j))ᵀ)
    {pp : ℕ} {Bm : Matrix (Fin pp) (Fin (d N)) ℝ}
    (hgram : s.gram N ω = ((d N : ℝ))⁻¹ • (Bmᵀ * Bm) + s.qmatR N ω U * (s.qmatR N ω U)ᵀ) :
    s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • (Bmᵀ * Bm) := by
  have h1 := s.gram_eq_rankRW0_add N ω hU hsig
  rw [h1] at hgram
  exact add_right_cancel hgram

/-- **The rank-`r` split of the stacked matrix**, plan section 1.1, on the shared interface.

The table is `X = A Vᵀ + E` (`hX`), so with `A Q_C = U Λ^{1/2}` the mean is
`U (V Q_C Λ^{1/2})ᵀ` and the Gram matrix splits as `W₀ + Q Qᵀ` with `W₀ = d⁻¹ Bᵀ B`. The
block `B` and the frame coefficients `Uᵀ Zu` are independent, each Gaussian.

Side condition. `ns N = r + p` is the plan's `r ≤ n` (decision D11 removes it later by the
`TailShift` device). -/
theorem exists_rankR_split (s : RankRStack μ ns d r) (hG : s.GaussianNoise) (N : ℕ)
    {pp : ℕ} (hn : ns N = r + pp) :
    ∃ (U : Matrix (Fin (ns N)) (Fin r) ℝ) (B : Ω N → Matrix (Fin pp) (Fin (d N)) ℝ),
      Uᵀ * U = 1 ∧
        s.signalPart N
          = U * (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j))ᵀ ∧
        (∀ ω, s.gram N ω = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)
          + s.qmatR N ω U * (s.qmatR N ω U)ᵀ) ∧
        HasLaw (fun ω => (Uᵀ * s.Zu N ω, B ω))
          ((gaussianMatrix r (d N)).prod (gaussianMatrix pp (d N))) (μ N) := by
  classical
  set Dm : Matrix (Fin r) (Fin r) ℝ := Matrix.diagonal fun j => Real.sqrt (s.coreEig j)
    with hDmdef
  have hrn : r ≤ ns N := by rw [hn]; exact Nat.le_add_right r pp
  have hA : (s.A N)ᵀ * s.A N
      = s.coreEigMat * Matrix.diagonal s.coreEig * (s.coreEigMat)ᵀ := by
    rw [s.hcore N, s.core_eq_conj]
  obtain ⟨U, hU, hfac⟩ := exists_orthonormal_factor hrn (s.A N) s.coreEigMat
    s.coreEig s.coreEigMat_transpose_mul_self hA
  obtain ⟨Yfun, hgramB, hmp⟩ := exists_frameBlock (d := d N) hn U hU
  -- 1. the mean of the table
  have hsig : s.signalPart N = U * (s.spikeMat N * Dm)ᵀ := by
    have h1 : s.A N = U * Dm * (s.coreEigMat)ᵀ := by
      calc s.A N
          = s.A N * (s.coreEigMat * (s.coreEigMat)ᵀ) := by
            rw [s.coreEigMat_mul_transpose_self, Matrix.mul_one]
        _ = s.A N * s.coreEigMat * (s.coreEigMat)ᵀ := by rw [Matrix.mul_assoc]
        _ = U * Dm * (s.coreEigMat)ᵀ := by rw [hfac]
    rw [signalPart, h1, spikeMat, Matrix.transpose_mul, Matrix.transpose_mul,
      hDmdef, Matrix.diagonal_transpose]
    simp only [Matrix.mul_assoc]
  -- 2. the noise, scaled and split
  have hsqd : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (by positivity)]
  have hFsmul : ∀ ω, s.E N ω - U * (Uᵀ * s.E N ω)
      = (Real.sqrt (d N))⁻¹ • (s.Zu N ω - U * (Uᵀ * s.Zu N ω)) := by
    intro ω
    rw [s.hE N ω, Matrix.mul_smul, Matrix.mul_smul, ← smul_sub]
  have hFF : ∀ ω, (s.E N ω - U * (Uᵀ * s.E N ω))ᵀ * (s.E N ω - U * (Uᵀ * s.E N ω))
      = ((d N : ℝ))⁻¹ • ((Yfun (s.Zu N ω))ᵀ * Yfun (s.Zu N ω)) := by
    intro ω
    rw [hFsmul ω, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, hsqd,
      gram_perp_eq U (s.Zu N ω) hU, ← hgramB (s.Zu N ω)]
  -- 3. the table as `U Qᵀ + E⊥`
  have hsplitX : ∀ ω, s.X N ω
      = U * (s.qmatR N ω U)ᵀ + (s.E N ω - U * (Uᵀ * s.E N ω)) := by
    intro ω
    have hT : ((s.E N ω)ᵀ * U)ᵀ = Uᵀ * s.E N ω := by
      rw [Matrix.transpose_mul, Matrix.transpose_transpose]
    rw [s.X_eq N ω, hsig, qmatR, ← hDmdef, Matrix.transpose_add, Matrix.mul_add, hT]
    abel
  refine ⟨U, fun ω => Yfun (s.Zu N ω), hU, hsig, ?_, ?_⟩
  · intro ω
    rw [gram, hsplitX ω,
      gram_mul_transpose_add U (s.qmatR N ω U) _ hU
        (transpose_perp_mul U (s.E N ω) hU), hFF ω]
  · exact hmp.fun_comp_hasLaw (hG N)

end RankRStack

/-! ### 5. The stack of `UnalignedModel` -/

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The tables of the model with every spike replaced by the spike of table `0`. Only the
noise fields are read below and those are untouched, so the stacked noise law of
`MultiTableModel` transfers. `Fin M` is empty when `M = 0`, so no shared spike has to be
produced in that case. -/
noncomputable def tblShared (m : UnalignedModel μ M n d r) (i : Fin M) :
    SpikedModel μ (n i) d :=
  { m.tbl i with
    v := fun N => (m.tbl ⟨0, i.pos⟩).v N
    hv := fun N => (m.tbl ⟨0, i.pos⟩).hv N }

/-- The shared-spike copy of the model, as a `MultiTableModel`. It has the same noise as `m`
and a different, unused, spike. -/
noncomputable def toMultiTable (m : UnalignedModel μ M n d r) : MultiTableModel μ M n d where
  tbl := m.tblShared
  hv := fun _ _ _ => rfl

/-- The stacked unscaled noise `Z_stack`, rows reindexed as `stackX` and `stackE` do. -/
noncomputable def stackZu (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (q : (i : Fin M) × Fin (n i N)) k => (m.tbl q.1).Z N ω q.2 k)

theorem stackZu_eq (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.stackZu N ω = m.toMultiTable.stackZ N ω := rfl

/-- **The law of the stacked noise.** `MultiTableModel.stack_law` on the shared-spike copy.
`NeZero M` is the side condition of `MultiTableModel.stack`, which reads table `0`. -/
theorem hasLaw_stackZu [NeZero M] (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    (N : ℕ) :
    HasLaw (m.stackZu N) (gaussianMatrix (∑ i, n i N) (d N)) (μ N) :=
  m.toMultiTable.stack_law hG N

/-- `E_stack = d^{-1/2} Z_stack`. -/
theorem stackE_eq_smul (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.stackE N ω = (Real.sqrt (d N))⁻¹ • m.stackZu N ω := by
  ext q k
  rw [stackE_apply', Matrix.smul_apply, smul_eq_mul]
  rfl

/-- The stacked unscaled noise is measurable, entry by entry from `SpikedModel.hZ`. -/
theorem measurable_stackZu (m : UnalignedModel μ M n d r) (N : ℕ) :
    Measurable (m.stackZu N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => m.stackZu N ω q l)
      = fun ω => (m.tbl (finSigmaFinEquiv.symm q).1).Z N ω
          (finSigmaFinEquiv.symm q).2 l := rfl
  rw [h]
  exact (measurable_pi_apply l).comp
    ((measurable_pi_apply _).comp ((m.tbl (finSigmaFinEquiv.symm q).1).hZ N))

/-- **The `r_i = 1` stack as `RankRStack` data.** Every field is a field or a proved
identity of the model, so each derived object of `RankRStack` is definitionally the matching
object of `UnalignedModel` and the bridge lemmas below are `rfl`. `NeZero M` supplies
`0 < d N`, which the model carries only inside a table. -/
noncomputable def toStack [NeZero M] (m : UnalignedModel μ M n d r) :
    RankRStack μ (fun N => ∑ i, n i N) d r where
  V := m.V
  A := m.signalFactor
  core := m.core
  Zu := m.stackZu
  E := m.stackE
  X := m.stackX
  hV := m.hV
  hcore := m.signalFactor_transpose_mul_self
  hE := m.stackE_eq_smul
  hX := m.stackX_eq
  hd := fun N => (m.tbl ⟨0, NeZero.pos M⟩).hd N
  hZmeas := m.measurable_stackZu

/-- Joint Gaussian noise of the tables gives Gaussian noise of the stack data. -/
theorem gaussianNoise_toStack [NeZero M] (m : UnalignedModel μ M n d r)
    (hG : m.JointGaussianNoise) : m.toStack.GaussianNoise :=
  fun N => m.hasLaw_stackZu hG N

theorem toStack_Zu [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.toStack.Zu N ω = m.stackZu N ω := rfl

theorem toStack_E [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.toStack.E N ω = m.stackE N ω := rfl

theorem toStack_signalPart [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) :
    m.toStack.signalPart N = m.signalPart N := rfl

theorem toStack_gram [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.toStack.gram N ω = m.stackGram N ω := rfl

theorem toStack_core [NeZero M] (m : UnalignedModel μ M n d r) :
    m.toStack.core = m.core := rfl

/-! #### The eigenvector matrix of the core, and the spike matrix -/

/-- `Q_C`, the matrix whose columns are the eigenvectors of the core matrix `C`. -/
noncomputable def coreEigMat (m : UnalignedModel μ M n d r) : Matrix (Fin r) (Fin r) ℝ :=
  (m.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ)

theorem coreEigMat_apply (m : UnalignedModel μ M n d r) (k j : Fin r) :
    m.coreEigMat k j = m.isHermitian_core.eigenvectorBasis j k := rfl

theorem coreEigMat_transpose_mul_self (m : UnalignedModel μ M n d r) :
    (m.coreEigMat)ᵀ * m.coreEigMat = 1 := by
  have h : star (m.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) *
      (m.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) = 1 :=
    Unitary.coe_star_mul_self _
  have hst : (m.coreEigMat)ᵀ
      = star (m.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) := by
    ext k l
    simp [coreEigMat, Matrix.star_apply]
  rw [hst]
  exact h

theorem coreEigMat_mul_transpose_self (m : UnalignedModel μ M n d r) :
    m.coreEigMat * (m.coreEigMat)ᵀ = 1 :=
  mul_eq_one_comm.1 m.coreEigMat_transpose_mul_self

/-- `C = Q_C Λ Q_Cᵀ`, the spectral theorem in the matrix form the factorization takes. -/
theorem core_eq_conj (m : UnalignedModel μ M n d r) :
    m.core = m.coreEigMat * Matrix.diagonal m.coreEig * (m.coreEigMat)ᵀ := by
  ext k l
  conv_lhs => rw [m.isHermitian_core.spectral_theorem]
  simp only [Unitary.conjStarAlgAut_apply, Matrix.mul_apply, Matrix.diagonal_apply,
    Matrix.star_apply, Matrix.transpose_apply, Matrix.IsHermitian.eigenvectorUnitary_apply,
    star_trivial, Function.comp_apply, RCLike.ofReal_real_eq_id, id_eq, mul_ite, mul_zero,
    Finset.sum_ite_eq', Finset.mem_univ, if_true, coreEigMat, coreEig]

/-- `V Q_C`, the `d × r` matrix whose columns are the spike directions `V q_j`. -/
noncomputable def spikeMat (m : UnalignedModel μ M n d r) (N : ℕ) :
    Matrix (Fin (d N)) (Fin r) ℝ := m.V N * m.coreEigMat

theorem spikeMat_apply (m : UnalignedModel μ M n d r) (N : ℕ) (k : Fin (d N)) (j : Fin r) :
    m.spikeMat N k j = WithLp.ofLp (m.spikeVec j N) k := by
  simp [spikeMat, coreEigMat, spikeVec, Matrix.mul_apply, Matrix.mulVec, dotProduct]

/-- The spike directions are orthonormal: `(V Q_C)ᵀ (V Q_C) = 1`. -/
theorem spikeMat_transpose_mul_self (m : UnalignedModel μ M n d r) (N : ℕ) :
    (m.spikeMat N)ᵀ * m.spikeMat N = 1 := by
  rw [spikeMat, Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (m.V N)ᵀ,
    m.hV N, Matrix.one_mul, m.coreEigMat_transpose_mul_self]

/-- `G = E_stackᵀ U` on the unscaled noise: `G = d^{-1/2} (Uᵀ Z_stack)ᵀ`. This is the bridge
between the split objects and the law of `exists_rankR_split`. -/
theorem stackE_transpose_mul (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    (m.stackE N ω)ᵀ * U = (Real.sqrt (d N))⁻¹ • ((Uᵀ * m.stackZu N ω)ᵀ) := by
  rw [m.stackE_eq_smul N ω, Matrix.transpose_smul, Matrix.smul_mul, Matrix.transpose_mul,
    Matrix.transpose_transpose]

/-- The `d × r` matrix `Q = V Q_C Λ^{1/2} + E_stackᵀ U` of plan section 1.1 item 2. Column
`j` is `√λ_j (V q_j) + g_j`. -/
noncomputable def qmatR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  m.spikeMat N * Matrix.diagonal (fun j => Real.sqrt (m.coreEig j)) + (m.stackE N ω)ᵀ * U

/-- Entry `(k, j)` of `Q`: column `j` is `√λ_j (V q_j) + g_j`. -/
theorem qmatR_apply (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (k : Fin (d N)) (j : Fin r) :
    m.qmatR N ω U k j
      = Real.sqrt (m.coreEig j) * WithLp.ofLp (m.spikeVec j N) k
        + ((m.stackE N ω)ᵀ * U) k j := by
  rw [qmatR, Matrix.add_apply, Matrix.mul_diagonal, m.spikeMat_apply N k j, mul_comm]

/-! #### The bridge to `RankRStack` -/

theorem toStack_coreEig [NeZero M] (m : UnalignedModel μ M n d r) :
    m.toStack.coreEig = m.coreEig := rfl

theorem toStack_coreEigMat [NeZero M] (m : UnalignedModel μ M n d r) :
    m.toStack.coreEigMat = m.coreEigMat := rfl

theorem toStack_spikeMat [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) :
    m.toStack.spikeMat N = m.spikeMat N := rfl

theorem toStack_spikeVec [NeZero M] (m : UnalignedModel μ M n d r) (j : Fin r) (N : ℕ) :
    m.toStack.spikeVec j N = m.spikeVec j N := rfl

theorem toStack_qmatR [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    m.toStack.qmatR N ω U = m.qmatR N ω U := rfl

/-! #### The split of the stack -/

/-- **The rank-`r` split of the stacked matrix**, plan section 1.1 at `A = signalFactor`.

The stack is `X_stack = A Vᵀ + E_stack` (`stackX_eq`), so with `A Q_C = U Λ^{1/2}` the mean
is `U (V Q_C Λ^{1/2})ᵀ` and the Gram matrix splits as `W₀ + Q Qᵀ` with `W₀ = d⁻¹ Bᵀ B`.
The block `B` and the frame coefficients `Uᵀ Z_stack` are independent, each Gaussian.

Side conditions. `∑ i, n i N = r + p` is the plan's `r ≤ n` (decision D11 removes it later
by the `TailShift` device). `NeZero M` comes from `MultiTableModel.stack`. -/
theorem exists_rankR_split [NeZero M] (m : UnalignedModel μ M n d r)
    (hG : m.JointGaussianNoise) (N : ℕ) {pp : ℕ} (hn : ∑ i, n i N = r + pp) :
    ∃ (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (B : Ω N → Matrix (Fin pp) (Fin (d N)) ℝ),
      Uᵀ * U = 1 ∧
        m.signalPart N
          = U * (m.spikeMat N * Matrix.diagonal fun j => Real.sqrt (m.coreEig j))ᵀ ∧
        (∀ ω, m.stackGram N ω = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)
          + m.qmatR N ω U * (m.qmatR N ω U)ᵀ) ∧
        HasLaw (fun ω => (Uᵀ * m.stackZu N ω, B ω))
          ((gaussianMatrix r (d N)).prod (gaussianMatrix pp (d N))) (μ N) :=
  m.toStack.exists_rankR_split (m.gaussianNoise_toStack hG) N hn

end UnalignedModel

end StackedSVD
