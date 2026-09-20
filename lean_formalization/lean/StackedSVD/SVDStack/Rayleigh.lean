/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Weighted

/-!
# P1: the weight-free row-space bound for weighted svdstack

Sixth module of `StackedSVD.SVDStack`. It proves the bound of `notes/archive/P1_rayleigh_bound.md`:
no weight vector, data dependent or not, beats `S/(S+1)` (`main_paper.tex:502`,
`thm:svdstack_weighted`; the paper's proof at `main_paper.tex:1274` covers only the weights
`w` for which `W A_beta W` has a simple top eigenvalue).

The mathematics. Write `Vt = m.Vt N omega` for the `M x d` matrix of the top right singular
vectors, `v` for the shared spike, `g = Vt v` and `G = Vt Vtᵀ`. For a weight vector `w` with
one nonzero entry the top eigenvalue of `Vt_wᵀ Vt_w` is positive, so its top eigenspace sits
inside the range of `Vt_wᵀ`, which sits inside the range of `Vtᵀ`. The orthogonal projector is
monotone in the subspace, so the performance is at most `rowBound`, the squared norm of the
projection of `v` on the row space of `Vt`. That number does not read `w` at all. When `G` is
invertible it equals `g ⬝ᵥ (G⁻¹ g)`, which tends in probability to
`beta ⬝ᵥ (A_beta⁻¹ beta) = S/(S+1)`.

## Content

1. `rowBound`, `rowBoundClosed`: the projector form and the closed form.
2. `svdstackPerfW_le_rowBound` (a): the almost sure bound, uniform in `w`.
3. `rowBound_eq_closed` (b): the closed form on `{det G is a unit}`.
4. `det_Abeta_pos`, `dotProduct_inv_Abeta`: the two scalar facts about `A_beta`.
5. `rowBound_tendsto` (c): `rowBound` tends in probability to `svdstackLimitOpt beta`.
6. `svdstackPerfW_uniform_bound` (d) and `svdstackPerfW_le_opt_whp`.
7. `svdstackLimitOpt_eq_of_single`: case 3 of the paper's proof (`main_paper.tex:1219`).
8. `thm_svdstack_weighted_gaussian_opt_full`: the Gaussian facade.

The three scalar lemmas `det_Abeta_pos`, `dotProduct_inv_Abeta` and
`svdstackLimitOpt_eq_of_single` sit in `namespace StackedSVD`, not in
`namespace StackedSVD.MultiTableModel`: they take no model, exactly like `abeta_gap` and
`svdstackLimitW_le_opt`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

/-! ### Helpers copied from `SVDStack/Defs.lean`

`toOp_mul_apply`, `inner_toEuclideanLin_transpose` and `sum_sq_ofLp` are private in
`SVDStack/Defs.lean` and in `SVDStack/Weighted.lean`. This file may not edit those files, so
the three proofs are copied verbatim under primed names. -/

section RayleighHelpers

variable {dd : ℕ}

/-- The operator of a matrix product is the composition of the operators. Copy of the private
`toOp_mul_apply` of `SVDStack/Defs.lean`. -/
private theorem toOp_mul_apply' {p q r : ℕ} (B : Matrix (Fin p) (Fin q) ℝ)
    (A : Matrix (Fin q) (Fin r) ℝ) (x : EuclideanSpace ℝ (Fin r)) :
    Matrix.toEuclideanLin (B * A) x = Matrix.toEuclideanLin B (Matrix.toEuclideanLin A x) := by
  apply WithLp.ofLp_injective
  change (B * A) *ᵥ WithLp.ofLp x = B *ᵥ (A *ᵥ WithLp.ofLp x)
  exact (Matrix.mulVec_mulVec _ _ _).symm

/-- The transpose is the adjoint of a real matrix operator. Copy of the private
`inner_toEuclideanLin_transpose` of `SVDStack/Defs.lean`. -/
private theorem inner_toEuclideanLin_transpose' {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (x : EuclideanSpace ℝ (Fin q)) (y : EuclideanSpace ℝ (Fin p)) :
    ⟪Matrix.toEuclideanLin A x, y⟫_ℝ = ⟪x, Matrix.toEuclideanLin Aᵀ y⟫_ℝ := by
  rw [real_inner_eq_dotProduct, real_inner_eq_dotProduct]
  change (A *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp y = WithLp.ofLp x ⬝ᵥ (Aᵀ *ᵥ WithLp.ofLp y)
  rw [Matrix.mulVec_transpose, dotProduct_comm, Matrix.dotProduct_mulVec, dotProduct_comm]

/-- `∑_l y_l² = ‖y‖²`. Copy of the private `sum_sq_ofLp` of `SVDStack/Weighted.lean`. -/
private theorem sum_sq_ofLp' {p : ℕ} (y : EuclideanSpace ℝ (Fin p)) :
    ∑ l, WithLp.ofLp y l ^ 2 = ‖y‖ ^ 2 := by
  rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, dotProduct]
  exact Finset.sum_congr rfl fun l _ => by ring

end RayleighHelpers

/-! ### Step 3 and step 4: containment and monotone projections -/

section Containment

/-- Step 4 of the note. The orthogonal projector is monotone in the subspace:
`K ≤ R` gives `‖P_K y‖ ≤ ‖P_R y‖`. Route: `P_K = P_K ∘ P_R` and `P_K` is a contraction. -/
theorem norm_starProjection_le_of_le {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] {K R : Submodule ℝ E} [K.HasOrthogonalProjection]
    [R.HasOrthogonalProjection] (h : K ≤ R) (y : E) :
    ‖K.starProjection y‖ ≤ ‖R.starProjection y‖ := by
  have h1 : K.starProjection (R.starProjection y) = K.starProjection y := by
    have h2 := Submodule.starProjection_comp_starProjection_of_le (𝕜 := ℝ) h
    exact congrArg (fun T : E →L[ℝ] E => T y) h2
  rw [← h1]
  exact Submodule.norm_starProjection_apply_le K _

/-- Rayleigh form of a Gram matrix: `⟪Bᵀ B x, x⟫ = ‖B x‖²`. -/
theorem inner_toOp_transpose_mul_self {a b : ℕ} (B : Matrix (Fin a) (Fin b) ℝ)
    (x : EuclideanSpace ℝ (Fin b)) :
    ⟪toOp (Bᵀ * B) x, x⟫_ℝ = ‖Matrix.toEuclideanLin B x‖ ^ 2 := by
  have h1 : toOp (Bᵀ * B) x = Matrix.toEuclideanLin Bᵀ (Matrix.toEuclideanLin B x) :=
    toOp_mul_apply' Bᵀ B x
  rw [h1, inner_toEuclideanLin_transpose' Bᵀ (Matrix.toEuclideanLin B x) x,
    Matrix.transpose_transpose, real_inner_self_eq_norm_sq]

/-- Step 3 of the note. For a nonzero top eigenvalue the top eigenspace of `Bᵀ B` sits inside
the range of `Bᵀ`: `x = Bᵀ (λ⁻¹ B x)`. -/
theorem topSpace_transpose_mul_self_le_range {a b : ℕ} (B : Matrix (Fin a) (Fin b) ℝ)
    (hne : lamMax (Bᵀ * B) (isHermitian_transpose_mul_self B) ≠ 0) :
    topSpace (Bᵀ * B) (isHermitian_transpose_mul_self B)
      ≤ LinearMap.range (Matrix.toEuclideanLin Bᵀ) := by
  intro x hx
  have hev : toOp (Bᵀ * B) x
      = lamMax (Bᵀ * B) (isHermitian_transpose_mul_self B) • x := toOp_of_mem_topSpace hx
  refine ⟨(lamMax (Bᵀ * B) (isHermitian_transpose_mul_self B))⁻¹ •
    Matrix.toEuclideanLin B x, ?_⟩
  rw [map_smul, ← toOp_mul_apply' Bᵀ B x, hev, smul_smul, inv_mul_cancel₀ hne, one_smul]

/-- The range of `(diag w * V)ᵀ` sits inside the range of `Vᵀ`, because
`(diag w * V)ᵀ = Vᵀ (diag w)ᵀ`. -/
theorem range_transpose_diagonal_mul_le {p q : ℕ} (w : Fin p → ℝ)
    (V : Matrix (Fin p) (Fin q) ℝ) :
    LinearMap.range (Matrix.toEuclideanLin (Matrix.diagonal w * V)ᵀ)
      ≤ LinearMap.range (Matrix.toEuclideanLin Vᵀ) := by
  rintro x ⟨c, rfl⟩
  refine ⟨Matrix.toEuclideanLin (Matrix.diagonal w)ᵀ c, ?_⟩
  rw [← toOp_mul_apply' Vᵀ (Matrix.diagonal w)ᵀ c, ← Matrix.transpose_mul]

end Containment

/-! ### The two scalar facts about `A_β` -/

section AbetaScalars

variable {M : ℕ}

/-- `A_β = β βᵀ + diag(1 - β_i²)` is positive definite when `0 ≤ β_i < 1`, so its determinant
is positive (step 7 of `notes/archive/P1_rayleigh_bound.md`). -/
theorem det_Abeta_pos (β : Fin M → ℝ) (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) :
    0 < (Abeta β).det := by
  have hpos : ∀ i, (0 : ℝ) < 1 - β i ^ 2 := fun i => by
    have h1 := (hβ i).1
    have h2 := (hβ i).2
    nlinarith
  have hd : (Matrix.diagonal fun i => 1 - β i ^ 2).PosDef :=
    Matrix.PosDef.diagonal fun i => hpos i
  have hr : (Matrix.vecMulVec β β).PosSemidef := by
    have h := Matrix.posSemidef_vecMulVec_self_star β
    have hstar : (star β : Fin M → ℝ) = β := by
      funext i
      simp
    rwa [hstar] at h
  have hpd : (Abeta β).PosDef := by
    rw [Abeta]
    exact Matrix.PosDef.posSemidef_add hr hd
  exact hpd.det_pos

/-- `βᵀ A_β⁻¹ β = S/(S+1)` (step 8 of `notes/archive/P1_rayleigh_bound.md`). The witness is
`x_i = β_i / ((S+1)(1-β_i²))`, for which `A_β x = β`. -/
theorem dotProduct_inv_Abeta (β : Fin M → ℝ) (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) :
    β ⬝ᵥ ((Abeta β)⁻¹.mulVec β) = svdstackLimitOpt β := by
  have hpos : ∀ i, (0 : ℝ) < 1 - β i ^ 2 := fun i => by
    have h1 := (hβ i).1
    have h2 := (hβ i).2
    nlinarith
  have hS0 : 0 ≤ Sval β :=
    Finset.sum_nonneg fun i _ => div_nonneg (sq_nonneg _) (hpos i).le
  have hS1 : (0 : ℝ) < Sval β + 1 := by linarith
  have hS1' : Sval β + 1 ≠ 0 := ne_of_gt hS1
  set x : Fin M → ℝ := fun i => β i / ((Sval β + 1) * (1 - β i ^ 2)) with hxdef
  have hxi : ∀ i, x i = β i / ((Sval β + 1) * (1 - β i ^ 2)) := fun i => by rw [hxdef]
  have hbx : β ⬝ᵥ x = Sval β / (Sval β + 1) := by
    have hterm : ∀ i : Fin M, β i * x i = β i ^ 2 / (1 - β i ^ 2) / (Sval β + 1) := by
      intro i
      have hb := (hpos i).ne'
      rw [hxi i]
      field_simp
    calc β ⬝ᵥ x = ∑ i, β i * x i := rfl
      _ = ∑ i, β i ^ 2 / (1 - β i ^ 2) / (Sval β + 1) :=
          Finset.sum_congr rfl fun i _ => hterm i
      _ = (∑ i, β i ^ 2 / (1 - β i ^ 2)) / (Sval β + 1) := (Finset.sum_div _ _ _).symm
      _ = Sval β / (Sval β + 1) := rfl
  have hAx : (Abeta β) *ᵥ x = β := by
    funext i
    have hrow : ((Abeta β) *ᵥ x) i = (β ⬝ᵥ x) * β i + (1 - β i ^ 2) * x i := by
      rw [Abeta, Matrix.add_mulVec, Pi.add_apply, Matrix.vecMulVec_mulVec, Pi.smul_apply,
        op_smul_eq_mul, Matrix.mulVec_diagonal]
      ring
    have hb := (hpos i).ne'
    rw [hrow, hbx, hxi i]
    field_simp
  have hdet : IsUnit (Abeta β).det := (det_Abeta_pos β hβ).ne'.isUnit
  have hinv : (Abeta β)⁻¹ *ᵥ β = x := by
    calc (Abeta β)⁻¹ *ᵥ β = (Abeta β)⁻¹ *ᵥ ((Abeta β) *ᵥ x) := by rw [hAx]
      _ = x := by
          rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
  rw [hinv, hbx]
  rfl

/-- Case 3 of the paper's proof (`main_paper.tex:1219`): when at most one table is above the
threshold, `S/(S+1) = β_k²`, so no weighting beats the best single table. -/
theorem svdstackLimitOpt_eq_of_single (β : Fin M → ℝ) (k : Fin M)
    (hβ : ∀ i, 0 ≤ β i ∧ β i < 1) (hk : ∀ i, i ≠ k → β i = 0) :
    svdstackLimitOpt β = β k ^ 2 := by
  have hpos : (0 : ℝ) < 1 - β k ^ 2 := by
    have h1 := (hβ k).1
    have h2 := (hβ k).2
    nlinarith
  have hSv : Sval β = β k ^ 2 / (1 - β k ^ 2) := by
    rw [Sval, Finset.sum_eq_single k]
    · intro i _ hi
      rw [hk i hi]
      norm_num
    · intro h
      exact absurd (Finset.mem_univ k) h
  have hb := hpos.ne'
  have hden : β k ^ 2 / (1 - β k ^ 2) + 1 ≠ 0 := by
    have hrw : β k ^ 2 / (1 - β k ^ 2) + 1 = 1 / (1 - β k ^ 2) := by
      field_simp
      ring
    rw [hrw]
    exact ne_of_gt (div_pos one_pos hpos)
  rw [svdstackLimitOpt, hSv]
  field_simp
  ring

end AbetaScalars

/-! ### The continuous map of step 6

`rbPhi` is the map `(g, G) ↦ (det G)⁻¹ ∑_{ij} g_i adj(G)_{ij} g_j`, that is `g ⬝ᵥ (G⁻¹ g)`
written so that continuity is a sum of products of polynomials. The index type carries the
`M` entries of `g` on the left and the `M²` entries of `G` on the right. -/

section Phi

variable {M : ℕ}

/-- The `M × M` block of a point of `(Fin M ⊕ Fin M × Fin M) → ℝ`. -/
def rbMat (z : (Fin M ⊕ Fin M × Fin M) → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.of fun i j => z (Sum.inr (i, j))

/-- `g ⬝ᵥ (G⁻¹ g)` as an explicit sum of polynomials divided by `det G`. -/
noncomputable def rbPhi (z : (Fin M ⊕ Fin M × Fin M) → ℝ) : ℝ :=
  ((rbMat z).det)⁻¹ *
    ∑ i, ∑ j, z (Sum.inl i) * (rbMat z).adjugate i j * z (Sum.inl j)

/-- `u ⬝ᵥ (A⁻¹ u)` in the form that `rbPhi` uses. `Matrix.inv_def` gives
`A⁻¹ = (det A)⁻¹ • adj A`. -/
theorem dotProduct_inv_eq (A : Matrix (Fin M) (Fin M) ℝ) (u : Fin M → ℝ) :
    u ⬝ᵥ (A⁻¹ *ᵥ u) = (A.det)⁻¹ * ∑ i, ∑ j, u i * A.adjugate i j * u j := by
  have hinv : ∀ i j, A⁻¹ i j = (A.det)⁻¹ * A.adjugate i j := by
    intro i j
    rw [Matrix.inv_def, Ring.inverse_eq_inv, Matrix.smul_apply, smul_eq_mul]
  have h1 : u ⬝ᵥ (A⁻¹ *ᵥ u) = ∑ i, ∑ j, u i * (A⁻¹ i j * u j) := by
    change ∑ i, u i * (A⁻¹ *ᵥ u) i = _
    refine Finset.sum_congr rfl fun i _ => ?_
    change u i * (∑ j, A⁻¹ i j * u j) = _
    rw [Finset.mul_sum]
  rw [h1, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => by rw [hinv i j]; ring

/-- `rbPhi` reads the closed form at any point that carries `A` and `u`. -/
theorem rbPhi_eq (A : Matrix (Fin M) (Fin M) ℝ) (u : Fin M → ℝ)
    (z : (Fin M ⊕ Fin M × Fin M) → ℝ) (hA : ∀ i j, z (Sum.inr (i, j)) = A i j)
    (hu : ∀ i, z (Sum.inl i) = u i) :
    rbPhi z = u ⬝ᵥ (A⁻¹ *ᵥ u) := by
  have hmat : rbMat z = A := by
    ext i j
    exact hA i j
  unfold rbPhi
  rw [hmat, dotProduct_inv_eq A u]
  congr 1
  exact Finset.sum_congr rfl fun i _ =>
    Finset.sum_congr rfl fun j _ => by rw [hu i, hu j]

/-- The `M × M` block is continuous in the point. -/
theorem continuous_rbMat :
    Continuous (fun z : (Fin M ⊕ Fin M × Fin M) → ℝ => rbMat z) :=
  continuous_matrix fun _ _ => continuous_apply _

/-- `rbPhi` is continuous at every point whose matrix block is invertible. -/
theorem continuousAt_rbPhi (u : (Fin M ⊕ Fin M × Fin M) → ℝ) (hdet : (rbMat u).det ≠ 0) :
    ContinuousAt (rbPhi (M := M)) u := by
  have hadj : Continuous (fun z : (Fin M ⊕ Fin M × Fin M) → ℝ => (rbMat z).adjugate) :=
    continuous_rbMat.matrix_adjugate
  have hnum : Continuous (fun z : (Fin M ⊕ Fin M × Fin M) → ℝ =>
      ∑ i, ∑ j, z (Sum.inl i) * (rbMat z).adjugate i j * z (Sum.inl j)) :=
    continuous_finsetSum _ fun i _ => continuous_finsetSum _ fun j _ =>
      ((continuous_apply _).mul (hadj.matrix_elem i j)).mul (continuous_apply _)
  have hdetc : Continuous (fun z : (Fin M ⊕ Fin M × Fin M) → ℝ => (rbMat z).det) :=
    continuous_rbMat.matrix_det
  change ContinuousAt (fun z => ((rbMat z).det)⁻¹ *
    ∑ i, ∑ j, z (Sum.inl i) * (rbMat z).adjugate i j * z (Sum.inl j)) u
  exact (hdetc.continuousAt.inv₀ hdet).mul hnum.continuousAt

/-- The limit point of step 6: `β` on the left, `A_β` on the right. -/
noncomputable def rbLim (β : Fin M → ℝ) : (Fin M ⊕ Fin M × Fin M) → ℝ :=
  Sum.elim β fun p => Abeta β p.1 p.2

end Phi

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-! ### The bound -/

/-- `‖P_R v‖²` with `P_R` the orthogonal projector onto the row space of `Ṽ` (the column space
of `Ṽᵀ`). It reads no weight vector (`notes/archive/P1_rayleigh_bound.md`, statement (a);
`main_paper.tex:502`). -/
noncomputable def rowBound (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) : ℝ :=
  ‖(LinearMap.range (Matrix.toEuclideanLin (m.Vt N ω)ᵀ)).starProjection ((m.tbl 0).v N)‖ ^ 2

/-- The closed form `g ⬝ᵥ (G⁻¹ g)` with `g = Ṽ v` and `G = Ṽ Ṽᵀ` (`Matrix.inv` is `0` when
`det G = 0`). The note writes the same value with a `let`; the `let` is expanded here so that
`rfl` steps see the two copies of `g`. -/
noncomputable def rowBoundClosed (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) : ℝ :=
  ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) ⬝ᵥ
    ((m.Vt N ω * (m.Vt N ω)ᵀ)⁻¹.mulVec
      ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))

omit [NeZero M] in
/-- Step 2 of the note. The Rayleigh quotient of `Ṽ_wᵀ Ṽ_w` at the unit vector `v̂_k` is at
least `w_k²`, so `λ_max(Ṽ_wᵀ Ṽ_w) ≥ w_k²`. Weighted twin of `one_le_lamMax_Abeta`. -/
theorem sq_le_lamMax_VtW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (k : Fin M) :
    w k ^ 2 ≤ lamMax ((m.VtW w N ω)ᵀ * m.VtW w N ω)
      (isHermitian_transpose_mul_self (m.VtW w N ω)) := by
  have hxn : ‖m.vhat k N ω‖ = 1 := m.norm_vhat k N ω
  have hkey : ⟪toOp ((m.VtW w N ω)ᵀ * m.VtW w N ω) (m.vhat k N ω), m.vhat k N ω⟫_ℝ
      = ‖Matrix.toEuclideanLin (m.VtW w N ω) (m.vhat k N ω)‖ ^ 2 :=
    inner_toOp_transpose_mul_self (m.VtW w N ω) (m.vhat k N ω)
  have hcomp :
      WithLp.ofLp (Matrix.toEuclideanLin (m.VtW w N ω) (m.vhat k N ω)) k = w k := by
    change (m.VtW w N ω *ᵥ WithLp.ofLp (m.vhat k N ω)) k = w k
    have h2 : (m.VtW w N ω *ᵥ WithLp.ofLp (m.vhat k N ω)) k
        = ∑ p, m.VtW w N ω k p * WithLp.ofLp (m.vhat k N ω) p := rfl
    have h3 : ∑ p, m.VtW w N ω k p * WithLp.ofLp (m.vhat k N ω) p
        = w k * ∑ p, WithLp.ofLp (m.vhat k N ω) p ^ 2 := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun p _ => ?_
      rw [m.vtW_apply w N ω k p]
      change w k * m.vhat k N ω p * m.vhat k N ω p = w k * m.vhat k N ω p ^ 2
      ring
    rw [h2, h3, sum_sq_ofLp' (m.vhat k N ω), hxn, one_pow, mul_one]
  have hlow : w k ^ 2 ≤ ‖Matrix.toEuclideanLin (m.VtW w N ω) (m.vhat k N ω)‖ ^ 2 := by
    rw [← sum_sq_ofLp' (Matrix.toEuclideanLin (m.VtW w N ω) (m.vhat k N ω)), ← hcomp]
    exact Finset.single_le_sum
      (f := fun l => WithLp.ofLp (Matrix.toEuclideanLin (m.VtW w N ω) (m.vhat k N ω)) l ^ 2)
      (fun l _ => sq_nonneg _) (Finset.mem_univ k)
  have hup := inner_toOp_self_le ((m.VtW w N ω)ᵀ * m.VtW w N ω)
    (isHermitian_transpose_mul_self (m.VtW w N ω)) (m.vhat k N ω)
  rw [hkey, hxn, one_pow, mul_one] at hup
  linarith

/-- (a) On one almost sure event, every nonzero weight vector has performance at most
`rowBound`. The event does not read `w`, so data dependent weights are covered
(`notes/archive/P1_rayleigh_bound.md`, statement (a)). -/
theorem svdstackPerfW_le_rowBound (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → m.svdstackPerfW w N ω ≤ m.rowBound N ω := by
  have hall : ∀ᵐ ω ∂(μ N), ∀ i : Fin M,
      TopSimple (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) :=
    ae_all_iff.mpr fun i => (law i).topSimple N
  filter_upwards [hall] with ω hω
  intro w hw
  obtain ⟨k, hk⟩ := hw
  -- Step 1, redone here: the identity holds for every `w` on this one event.
  have heq : m.svdstackGramW w N ω = (m.VtW w N ω)ᵀ * m.VtW w N ω := by
    ext p q
    rw [MultiTableModel.svdstackGramW, Matrix.sum_apply, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Matrix.smul_apply, m.P_eq_vecMulVec i N ω (hω i), Matrix.vecMulVec_apply,
      Matrix.transpose_apply, m.vtW_apply w N ω i p, m.vtW_apply w N ω i q]
    change w i ^ 2 * (m.vhat i N ω p * m.vhat i N ω q)
        = w i * m.vhat i N ω p * (w i * m.vhat i N ω q)
    ring
  have hH := isHermitian_transpose_mul_self (m.VtW w N ω)
  have hperf : m.svdstackPerfW w N ω
      = ‖topProj ((m.VtW w N ω)ᵀ * m.VtW w N ω) hH ((m.tbl 0).v N)‖ ^ 2 := by
    change ‖topProj (m.svdstackGramW w N ω) (m.isHermitian_svdstackGramW w N ω)
      ((m.tbl 0).v N)‖ ^ 2 = _
    rw [topProj_congr heq (m.isHermitian_svdstackGramW w N ω) hH]
  -- Step 2: the top eigenvalue is positive.
  have hwk : (0 : ℝ) < w k ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hpos : 0 < lamMax ((m.VtW w N ω)ᵀ * m.VtW w N ω) hH :=
    lt_of_lt_of_le hwk (m.sq_le_lamMax_VtW w N ω k)
  -- Step 3: containment in the row space of `Ṽ`.
  have hBdef : m.VtW w N ω = Matrix.diagonal w * m.Vt N ω := rfl
  have hle : topSpace ((m.VtW w N ω)ᵀ * m.VtW w N ω) hH
      ≤ LinearMap.range (Matrix.toEuclideanLin (m.Vt N ω)ᵀ) := by
    refine le_trans (topSpace_transpose_mul_self_le_range (m.VtW w N ω) hpos.ne') ?_
    rw [hBdef]
    exact range_transpose_diagonal_mul_le w (m.Vt N ω)
  -- Step 4: monotone projections.
  have hmono := norm_starProjection_le_of_le hle ((m.tbl 0).v N)
  have h0 : (0 : ℝ) ≤ ‖(topSpace ((m.VtW w N ω)ᵀ * m.VtW w N ω) hH).starProjection
      ((m.tbl 0).v N)‖ := norm_nonneg _
  rw [hperf]
  change ‖(topSpace ((m.VtW w N ω)ᵀ * m.VtW w N ω) hH).starProjection ((m.tbl 0).v N)‖ ^ 2
      ≤ m.rowBound N ω
  rw [MultiTableModel.rowBound]
  nlinarith [hmono, h0]

/-- (b) The closed form, when `Ṽ Ṽᵀ` is invertible (step 5 of the note, Route A: the
projection vector is identified directly, so no projector matrix is built). -/
theorem rowBound_eq_closed (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N)
    (hdet : IsUnit (m.Vt N ω * (m.Vt N ω)ᵀ).det) :
    m.rowBound N ω = m.rowBoundClosed N ω := by
  change ‖(LinearMap.range (Matrix.toEuclideanLin (m.Vt N ω)ᵀ)).starProjection
      ((m.tbl 0).v N)‖ ^ 2
    = ((m.Vt N ω) *ᵥ WithLp.ofLp ((m.tbl 0).v N)) ⬝ᵥ
      ((m.Vt N ω * (m.Vt N ω)ᵀ)⁻¹ *ᵥ ((m.Vt N ω) *ᵥ WithLp.ofLp ((m.tbl 0).v N)))
  set V := m.Vt N ω with hV
  set v := (m.tbl 0).v N with hv
  set g : Fin M → ℝ := V *ᵥ WithLp.ofLp v with hg
  set G : Matrix (Fin M) (Fin M) ℝ := V * Vᵀ with hG
  set b : Fin M → ℝ := G⁻¹ *ᵥ g with hb
  set z : EuclideanSpace ℝ (Fin (d N)) := Matrix.toEuclideanLin Vᵀ (WithLp.toLp 2 b) with hz
  have hVv : Matrix.toEuclideanLin V v = WithLp.toLp 2 g := by
    rw [hg]
    apply WithLp.ofLp_injective
    rfl
  have hVz : Matrix.toEuclideanLin V z = WithLp.toLp 2 g := by
    rw [hz]
    apply WithLp.ofLp_injective
    change V *ᵥ (Vᵀ *ᵥ b) = g
    rw [Matrix.mulVec_mulVec, ← hG, hb, Matrix.mulVec_mulVec,
      Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
  have hzero : Matrix.toEuclideanLin V (v - z) = 0 := by
    rw [map_sub, hVv, hVz, sub_self]
  have hproj : (LinearMap.range (Matrix.toEuclideanLin Vᵀ)).starProjection v = z := by
    refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero ⟨WithLp.toLp 2 b, hz.symm⟩ ?_
    rintro y ⟨cc, rfl⟩
    calc ⟪v - z, Matrix.toEuclideanLin Vᵀ cc⟫_ℝ
        = ⟪Matrix.toEuclideanLin Vᵀ cc, v - z⟫_ℝ := real_inner_comm _ _
      _ = ⟪cc, Matrix.toEuclideanLin (Vᵀ)ᵀ (v - z)⟫_ℝ :=
          inner_toEuclideanLin_transpose' Vᵀ cc (v - z)
      _ = 0 := by rw [Matrix.transpose_transpose, hzero, inner_zero_right]
  have hnorm : ‖z‖ ^ 2 = g ⬝ᵥ (G⁻¹ *ᵥ g) := by
    rw [← real_inner_self_eq_norm_sq]
    nth_rewrite 1 [hz]
    rw [inner_toEuclideanLin_transpose' Vᵀ (WithLp.toLp 2 b) z, Matrix.transpose_transpose,
      hVz, real_inner_eq_dotProduct]
    change b ⬝ᵥ g = g ⬝ᵥ (G⁻¹ *ᵥ g)
    rw [hb]
    exact dotProduct_comm _ _
  rw [hproj, hnorm]

/-! ### The limit -/

/-- The random point of step 6: `Ṽ v` on the left, `Ṽ Ṽᵀ` on the right. -/
noncomputable def rbF (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    (Fin M ⊕ Fin M × Fin M) → ℝ :=
  Sum.elim (fun i => (m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N)) i)
    fun p => (m.Vt N ω * (m.Vt N ω)ᵀ) p.1 p.2

/-- `align` and `gramEntries` in one family: the point of `rbF` tends to the point of
`rbLim β` coordinatewise. -/
theorem tendstoInProbPi_rbF (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProbPi μ (fun N ω => m.rbF N ω) (rbLim β) := by
  rintro (i | ⟨i, j⟩)
  · have h := m.align c law i
    rw [← hβdef i] at h
    exact h
  · exact m.gramEntries c β hβdef law hI i j

/-- (c) The bound converges in probability to `S/(S+1)`. No threshold hypothesis: below the
threshold every `β_i = 0`, `A_β = I` and the limit is `0` (`notes/archive/P1_rayleigh_bound.md`,
statement (c)). -/
theorem rowBound_tendsto (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.rowBound N ω) (svdstackLimitOpt β) := by
  have hβ01 : ∀ i, 0 ≤ β i ∧ β i < 1 := fun i => by
    rw [hβdef i]
    exact beta_mem_Ico (hc i)
  have hdetA : 0 < (Abeta β).det := det_Abeta_pos β hβ01
  have hmatLim : rbMat (rbLim β) = Abeta β := rfl
  have hmatF : ∀ (N : ℕ) (ω : Ω N), rbMat (m.rbF N ω) = m.Vt N ω * (m.Vt N ω)ᵀ := fun _ _ => rfl
  have hF := m.tendstoInProbPi_rbF c β hβdef law hI
  -- The closed form converges.
  have hphi : ContinuousAt (rbPhi (M := M)) (rbLim β) := by
    refine continuousAt_rbPhi (rbLim β) ?_
    rw [hmatLim]
    exact hdetA.ne'
  have hclosedLim : rbPhi (rbLim β) = svdstackLimitOpt β := by
    rw [rbPhi_eq (Abeta β) β (rbLim β) (fun i j => rfl) (fun i => rfl)]
    exact dotProduct_inv_Abeta β hβ01
  have hclosedF : ∀ (N : ℕ) (ω : Ω N), rbPhi (m.rbF N ω) = m.rowBoundClosed N ω := by
    intro N ω
    exact rbPhi_eq (m.Vt N ω * (m.Vt N ω)ᵀ)
      ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) (m.rbF N ω)
      (fun i j => rfl) (fun i => rfl)
  have hclosed : TendstoInProb μ (fun N ω => m.rowBoundClosed N ω) (svdstackLimitOpt β) := by
    have h := TendstoInProbPi.comp_continuous hphi hF
    rw [hclosedLim] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    exact hclosedF N ω
  -- The determinant converges to a positive number.
  have hdetlim : TendstoInProb μ (fun N ω => (m.Vt N ω * (m.Vt N ω)ᵀ).det) ((Abeta β).det) := by
    have hcont : ContinuousAt (fun z : (Fin M ⊕ Fin M × Fin M) → ℝ => (rbMat z).det)
        (rbLim β) := continuous_rbMat.matrix_det.continuousAt
    have h := TendstoInProbPi.comp_continuous hcont hF
    rw [hmatLim] at h
    exact h
  have hnull : Tendsto (fun N => μ N {ω | (m.Vt N ω * (m.Vt N ω)ᵀ).det = 0}) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω | (Abeta β).det ≤ |(m.Vt N ω * (m.Vt N ω)ᵀ).det - (Abeta β).det|}) ?_
      (hdetlim ((Abeta β).det) hdetA)
    intro N ω hω
    have hω' : (m.Vt N ω * (m.Vt N ω)ᵀ).det = 0 := hω
    change (Abeta β).det ≤ |(m.Vt N ω * (m.Vt N ω)ᵀ).det - (Abeta β).det|
    rw [hω', zero_sub, abs_neg, abs_of_pos hdetA]
  -- Transfer.
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hclosed
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | (m.Vt N ω * (m.Vt N ω)ᵀ).det = 0}) ?_ hnull
  intro N ω hω
  have hω' : m.rowBound N ω ≠ m.rowBoundClosed N ω := hω
  change (m.Vt N ω * (m.Vt N ω)ᵀ).det = 0
  by_contra hne
  exact hω' (m.rowBound_eq_closed N ω (isUnit_iff_ne_zero.mpr hne))

/-! ### The corollaries -/

/-- (d) Uniform corollary: with probability tending to one no nonzero weight vector, data
dependent or not, beats `S/(S+1) + ε`. The set is not assumed measurable; only monotonicity
and subadditivity of `μ N` enter (`notes/archive/P1_rayleigh_bound.md`, statement (d)). -/
theorem svdstackPerfW_uniform_bound (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
      svdstackLimitOpt β + ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0) := by
  intro ε hε
  have hrb := m.rowBound_tendsto c β hc hβdef law hI
  have hbad : ∀ N, μ N {ω | ¬ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) →
      m.svdstackPerfW w N ω ≤ m.rowBound N ω} = 0 :=
    fun N => ae_iff.mp (m.svdstackPerfW_le_rowBound c law N)
  have h1 : Tendsto (fun N => μ N {ω | ¬ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) →
      m.svdstackPerfW w N ω ≤ m.rowBound N ω}) atTop (𝓝 0) := by
    simp only [hbad]
    exact tendsto_const_nhds
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ¬ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) →
        m.svdstackPerfW w N ω ≤ m.rowBound N ω} ∪
      {ω | ε ≤ |m.rowBound N ω - svdstackLimitOpt β|}) ?_
    (tendsto_measure_zero_union h1 (hrb ε hε))
  intro N ω hω
  obtain ⟨w, hwne, hge⟩ := hω
  by_cases hall : ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → m.svdstackPerfW w N ω ≤ m.rowBound N ω
  · right
    have hb := hall w hwne
    change ε ≤ |m.rowBound N ω - svdstackLimitOpt β|
    have h2 : ε ≤ m.rowBound N ω - svdstackLimitOpt β := by linarith
    exact le_trans h2 (le_abs_self _)
  · left
    exact hall

/-- Per-`w` corollary, the form of the external audit. -/
theorem svdstackPerfW_le_opt_whp (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise)
    (w : Fin M → ℝ) (hw : ∃ k, w k ≠ 0) :
    ∀ ε > 0, Tendsto (fun N => μ N {ω | svdstackLimitOpt β + ε ≤ m.svdstackPerfW w N ω})
      atTop (𝓝 0) := by
  intro ε hε
  refine tendsto_measure_zero_of_subset
    (t := fun N => {ω | ∃ w' : Fin M → ℝ, (∃ k, w' k ≠ 0) ∧
      svdstackLimitOpt β + ε ≤ m.svdstackPerfW w' N ω}) ?_
    (m.svdstackPerfW_uniform_bound c β hc hβdef law hI ε hε)
  intro N ω hω
  exact ⟨w, hw, hω⟩

/-- Gaussian facade: the conclusion of `thm_svdstack_weighted_gaussian_opt`
(`SVDStack/Weighted.lean:1412`) with the uniform bound of P1 as a second conjunct
(`main_paper.tex:502`). -/
theorem thm_svdstack_weighted_gaussian_opt_full [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    (TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)
      ∧ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → (∀ i, 0 ≤ w i) →
          TopSimple (AbetaW w β) (isHermitian_AbetaW w β) →
          TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β)
          ∧ svdstackLimitW w β ≤ svdstackLimitOpt β)
    ∧ ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
        svdstackLimitOpt β + ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0) :=
  ⟨m.thm_svdstack_weighted_gaussian_opt c β hc hβdef hthr hreg hG,
    m.svdstackPerfW_uniform_bound c β hc hβdef
      (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
        (m.gaussianNoise_of_joint hG i)) hG.indepNoise⟩

end MultiTableModel

end StackedSVD
