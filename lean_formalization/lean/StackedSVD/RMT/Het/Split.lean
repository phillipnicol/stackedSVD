/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.StackSVDWeighted

/-!
# Item H0: the column split of the weighted stack

Task H0 of `notes/archive/plan_heterolaw_A.md` (section 2.1 and the H0 row of section 4). The
weighted stack is

```
X_w = [w_1 X_1; ...; w_M X_M] = ũ₀ vᵀ + Σ^{1/2} E,   Σ^{1/2} = diag(w_i I_{n_i}),
```

with `ũ₀ = (θ_i w_i u_i)_i` on `Fin (∑ i, n i N)` and `E` the **unweighted** stacked noise.
Item R0 splits the rows off `u`; here we split the **columns** off `v`, because a column of
`Σ^{1/2} E` has covariance `Σ/d` and the row split off `ũ₀` is not independent of the spike
unless `ũ₀` is an eigenvector of `Σ` (plan section 2.1). With

```
e = E v,   E⊥ = E - e vᵀ,   q = ũ₀ + Σ^{1/2} e,   W₀' = Σ^{1/2} E⊥ E⊥ᵀ Σ^{1/2},
```

the identity `X_w X_wᵀ = W₀' + q qᵀ` is an honest rank-one positive update of a matrix that
is independent of `q`, so the rank-one theory of `RMT/R4.lean` applies verbatim.

Two modeling choices, both deliberate:

1. `SigmaHalf` carries the **signed** weight `w_i`, not `|w_i|`. It is then exactly the
   matrix that turns the unweighted stacked noise into the weighted one
   (`stackW_E_eq`), which is what the algebra needs. The paper's `Σ^{1/2} = diag(|w_i|)`
   differs by a per-block sign that the sign symmetry of the noise law absorbs.
2. `EperpHet` keeps all `d N` columns and kills `v` (`EperpHet_mulVec_v`), rather than
   living on a `d - 1` dimensional complement. No reindexing is then needed anywhere in
   the deterministic layer; the `d - 1` block appears only inside
   `exists_block_hasLaw_het`, where the law is stated.

Paper: `main_paper.tex` lines 1409 to 1440 (`thm:stacksvd_weighted`, `eq:assumption4`).

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/Split.lean` exit 0; 0 `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### Model-free pieces of the column split -/

/-- The column analogue of the Gram computation of `R0.W0_eq_smul`: removing the `x` column
of `Z` removes the rank-one part `(Z x) (Z x)ᵀ` of `Z Zᵀ`. -/
theorem gram_sub_vecMulVec {r c : ℕ} (x : Fin c → ℝ) (Z : Matrix (Fin r) (Fin c) ℝ)
    (hx : x ⬝ᵥ x = 1) :
    (Z - Matrix.vecMulVec (Z *ᵥ x) x) * (Z - Matrix.vecMulVec (Z *ᵥ x) x)ᵀ
      = Z * Zᵀ - Matrix.vecMulVec (Z *ᵥ x) (Z *ᵥ x) := by
  have hF : (Zᵀ - Matrix.vecMulVec x (Z *ᵥ x))ᵀ *ᵥ x = 0 := by
    simpa using transpose_sub_vecMulVec_mulVec x Zᵀ hx
  have hsum : Matrix.vecMulVec x (Z *ᵥ x) + (Zᵀ - Matrix.vecMulVec x (Z *ᵥ x)) = Zᵀ := by
    abel
  have key := gram_vecMulVec_add x (Z *ᵥ x) (Zᵀ - Matrix.vecMulVec x (Z *ᵥ x)) hx hF
  rw [hsum, Matrix.transpose_transpose] at key
  have hT : (Zᵀ - Matrix.vecMulVec x (Z *ᵥ x))ᵀ = Z - Matrix.vecMulVec (Z *ᵥ x) x := by
    rw [Matrix.transpose_sub, Matrix.transpose_vecMulVec, Matrix.transpose_transpose]
  rw [← hT, Matrix.transpose_transpose, key]
  abel

/-- The argument of `R0.exists_block_hasLaw`, stated without a model: for a unit vector `x`
of `ℝ^q` with `q = p + 1` there is a map `Yfun` on `q × r` matrices (the last `p` rows of
`U Z`, with `U` orthogonal and first row `x`) whose Gram matrix is the Gram matrix of `Z`
with the `x` direction removed, and the pair `(Zᵀ x, Yfun Z)` has the product Gaussian law.
`R0.exists_block_hasLaw` runs the same six steps on a `SpikedModel`; it does not expose them
in this model-free form, and `RMT/R0.lean` is read-only here. -/
theorem exists_rowBlock_measurePreserving {q p : ℕ} (r : ℕ) (hq : q = p + 1)
    (x : EuclideanSpace ℝ (Fin q)) (hx : ‖x‖ = 1) :
    ∃ Yfun : Matrix (Fin q) (Fin r) ℝ → Matrix (Fin p) (Fin r) ℝ,
      (∀ Z : Matrix (Fin q) (Fin r) ℝ, (Yfun Z)ᵀ * Yfun Z
          = Zᵀ * Z - Matrix.vecMulVec (Zᵀ *ᵥ WithLp.ofLp x) (Zᵀ *ᵥ WithLp.ofLp x)) ∧
        MeasurePreserving
          (fun Z : Matrix (Fin q) (Fin r) ℝ => (Zᵀ *ᵥ WithLp.ofLp x, Yfun Z))
          (gaussianMatrix q r)
          ((Measure.pi fun _ : Fin r => gaussianReal 0 1).prod (gaussianMatrix p r)) := by
  classical
  -- 1. an orthonormal basis of `ℝ^q`, indexed by `Fin (p + 1)`, whose first vector is `x`.
  have hcard : Module.finrank ℝ (EuclideanSpace ℝ (Fin q)) = Fintype.card (Fin (p + 1)) := by
    rw [finrank_euclideanSpace_fin, Fintype.card_fin, hq]
  have hortho : Orthonormal ℝ
      (({0} : Set (Fin (p + 1))).domRestrict (fun _ : Fin (p + 1) => x)) := by
    refine ⟨fun _ => hx, ?_⟩
    intro i j hij
    exact absurd (Subtype.ext ((i.2 : (i : Fin (p + 1)) = 0).trans
      (j.2 : (j : Fin (p + 1)) = 0).symm)) hij
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq hcard
  have hb0 : WithLp.ofLp (b 0) = WithLp.ofLp x := by rw [hb 0 rfl]
  -- 2. the rotation matrix: row `i` is `b i`, read through `Fin (p + 1) ≃ Fin q`.
  set e : Fin (p + 1) → Fin q := Fin.cast hq.symm with hedef
  set U : Matrix (Fin (p + 1)) (Fin (p + 1)) ℝ :=
    Matrix.of fun i i' => WithLp.ofLp (b i) (e i') with hUdef
  have hsum : ∀ f : Fin q → ℝ, ∑ i' : Fin (p + 1), f (e i') = ∑ k : Fin q, f k := by
    intro f
    exact Equiv.sum_comp (finCongr hq.symm) f
  have hUUt : U * Uᵀ = 1 := by
    ext i i'
    have hb' := (orthonormal_iff_ite (𝕜 := ℝ)).1 b.orthonormal i i'
    rw [inner_euclidean_eq_dotProduct] at hb'
    rw [Matrix.one_apply, ← hb']
    change ∑ i'' : Fin (p + 1),
        WithLp.ofLp (b i) (e i'') * WithLp.ofLp (b i') (e i'') = _
    exact hsum fun k => WithLp.ofLp (b i) k * WithLp.ofLp (b i') k
  have hU : Uᵀ * U = 1 := mul_eq_one_comm.1 hUUt
  -- 3. the composite map `Z ↦ (row 0 of U Z, block of rows 1, …, p of U Z)`.
  set castR : Matrix (Fin q) (Fin r) ℝ → Matrix (Fin (p + 1)) (Fin r) ℝ :=
    fun Z => Matrix.of fun (i : Fin (p + 1)) (j : Fin r) => Z (e i) j with hcastdef
  set Y : Matrix (Fin q) (Fin r) ℝ → Matrix (Fin (p + 1)) (Fin r) ℝ :=
    fun Z => U * castR Z with hYdef
  have hmp : MeasurePreserving (fun Z : Matrix (Fin q) (Fin r) ℝ =>
      (Y Z 0, Matrix.of fun (i : Fin p) (j : Fin r) => Y Z i.succ j))
      (gaussianMatrix q r)
      ((Measure.pi fun _ : Fin r => gaussianReal 0 1).prod (gaussianMatrix p r)) :=
    (measurePreserving_rowSplit p r).comp
      ((measurePreserving_mul_left hU r).comp (measurePreserving_castRow hq))
  -- 4. the rows of `U Z` in terms of `b`.
  have hrow : ∀ (Z : Matrix (Fin q) (Fin r) ℝ) (i : Fin (p + 1)) (j : Fin r),
      Y Z i j = ∑ k : Fin q, WithLp.ofLp (b i) k * Z k j := by
    intro Z i j
    change ∑ i' : Fin (p + 1), WithLp.ofLp (b i) (e i') * Z (e i') j = _
    exact hsum fun k => WithLp.ofLp (b i) k * Z k j
  have hrow0 : ∀ Z : Matrix (Fin q) (Fin r) ℝ, Y Z 0 = Zᵀ *ᵥ WithLp.ofLp x := by
    intro Z
    funext j
    rw [hrow Z 0 j, hb0]
    change _ = ∑ k, Z k j * WithLp.ofLp x k
    exact Finset.sum_congr rfl fun k _ => mul_comm _ _
  -- 5. the Gram matrix of `U Z` splits into row `0` and the block.
  have hgramY : ∀ Z : Matrix (Fin q) (Fin r) ℝ,
      (Matrix.of fun (i : Fin p) (j : Fin r) => Y Z i.succ j)ᵀ *
          (Matrix.of fun (i : Fin p) (j : Fin r) => Y Z i.succ j)
        = Zᵀ * Z - Matrix.vecMulVec (Zᵀ *ᵥ WithLp.ofLp x) (Zᵀ *ᵥ WithLp.ofLp x) := by
    intro Z
    have hYY : (Y Z)ᵀ * Y Z = Zᵀ * Z := by
      rw [hYdef]
      simp only
      rw [Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Uᵀ U (castR Z), hU,
        Matrix.one_mul]
      exact gram_castRow hq Z
    ext j j'
    have hsplit : ∑ i : Fin (p + 1), Y Z i j * Y Z i j'
        = Y Z 0 j * Y Z 0 j' + ∑ i : Fin p, Y Z i.succ j * Y Z i.succ j' :=
      Fin.sum_univ_succ (fun i : Fin (p + 1) => Y Z i j * Y Z i j')
    have hL : ((Y Z)ᵀ * Y Z) j j' = ∑ i : Fin (p + 1), Y Z i j * Y Z i j' := rfl
    have hR : (Zᵀ * Z) j j' = ∑ i : Fin (p + 1), Y Z i j * Y Z i j' := by
      rw [← hL, hYY]
    change ∑ i : Fin p, Y Z i.succ j * Y Z i.succ j' = (Zᵀ * Z) j j' - _
    rw [hR, hsplit, hrow0 Z]
    change _ = _ + _ - (Zᵀ *ᵥ WithLp.ofLp x) j * (Zᵀ *ᵥ WithLp.ofLp x) j'
    ring
  -- 6. the block, and the two conclusions.
  refine ⟨fun Z => Matrix.of fun (i : Fin p) (j : Fin r) => Y Z i.succ j, hgramY, ?_⟩
  simpa only [hrow0] using hmp

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### The five split objects -/

/-- `Σ^{1/2} = diag(w_i I_{n_i})` on the stacked row index, through the same
`finSigmaFinEquiv` row order as `stack`. Signed, see choice 1 of the header. -/
noncomputable def SigmaHalf (_m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin (∑ i, n i N)) ℝ :=
  Matrix.diagonal fun r => w (finSigmaFinEquiv.symm r).1

/-- `ũ₀ = (θ_i w_i u_i)_i`, the stacked signal vector. Written on `stackW` so that
`X_w = ũ₀ vᵀ + E_w` is definitional. -/
noncomputable def u0Het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    Fin (∑ i, n i N) → ℝ :=
  (m.stackW w).θ • WithLp.ofLp ((m.stackW w).u N)

/-- `e = E v`, the component of the **unweighted** stacked noise along `v`. -/
noncomputable def eHet (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Fin (∑ i, n i N) → ℝ :=
  m.stack.E N ω *ᵥ WithLp.ofLp (m.stack.v N)

/-- `E⊥ = E - e vᵀ = E (I - v vᵀ)`, on all `d N` columns (choice 2 of the header). -/
noncomputable def EperpHet (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.stack.E N ω - Matrix.vecMulVec (m.eHet N ω) (WithLp.ofLp (m.stack.v N))

/-- `q = ũ₀ + Σ^{1/2} e`, the rank-one part of `X_w X_wᵀ`. -/
noncomputable def qHet (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    Fin (∑ i, n i N) → ℝ :=
  m.u0Het w N + m.SigmaHalf w N *ᵥ m.eHet N ω

/-- `W₀' = Σ^{1/2} E⊥ E⊥ᵀ Σ^{1/2}`, the noise Gram matrix on the `n` side. -/
noncomputable def W0het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (∑ i, n i N)) ℝ :=
  m.SigmaHalf w N * m.EperpHet N ω * (m.EperpHet N ω)ᵀ * m.SigmaHalf w N

/-! ### Deterministic algebra -/

omit [NeZero M] in
theorem transpose_SigmaHalf (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    (m.SigmaHalf w N)ᵀ = m.SigmaHalf w N :=
  Matrix.diagonal_transpose _

theorem dotProduct_v_self_het (m : MultiTableModel μ M n d) (N : ℕ) :
    WithLp.ofLp (m.stack.v N) ⬝ᵥ WithLp.ofLp (m.stack.v N) = 1 := by
  have h : ⟪m.stack.v N, m.stack.v N⟫_ℝ
      = WithLp.ofLp (m.stack.v N) ⬝ᵥ WithLp.ofLp (m.stack.v N) :=
    inner_euclidean_eq_dotProduct (m.stack.v N) (m.stack.v N)
  rw [real_inner_self_eq_norm_sq, m.stack.hv N] at h
  simpa using h.symm

omit [NeZero M] in
/-- The weighted stacked noise is `Σ^{1/2}` times the unweighted one. -/
theorem stackZW_eq_SigmaHalf_mul (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : m.stackZW w N ω = m.SigmaHalf w N * m.stackZ N ω := by
  ext r k
  rw [SigmaHalf, Matrix.diagonal_mul]
  simp [stackZW, stackZ, Matrix.reindex_apply, Matrix.submatrix_apply]

theorem stackW_E_eq (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.stackW w).E N ω = m.SigmaHalf w N * m.stack.E N ω := by
  change (Real.sqrt (d N))⁻¹ • m.stackZW w N ω
      = m.SigmaHalf w N * ((Real.sqrt (d N))⁻¹ • m.stackZ N ω)
  rw [stackZW_eq_SigmaHalf_mul, Matrix.mul_smul]

/-- `E = e vᵀ + E⊥`, the column analogue of `R0.E_eq_vecMulVec_add`. -/
theorem stack_E_eq_vecMulVec_add (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.stack.E N ω
      = Matrix.vecMulVec (m.eHet N ω) (WithLp.ofLp (m.stack.v N)) + m.EperpHet N ω := by
  change _ = _ + (m.stack.E N ω - Matrix.vecMulVec (m.eHet N ω) (WithLp.ofLp (m.stack.v N)))
  abel

/-- `E⊥ v = 0`, the column analogue of `R0.transpose_Eperp_mulVec_u`. -/
theorem EperpHet_mulVec_v (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.EperpHet N ω *ᵥ WithLp.ofLp (m.stack.v N) = 0 := by
  have hv : ∑ k, WithLp.ofLp (m.stack.v N) k * WithLp.ofLp (m.stack.v N) k = 1 :=
    m.dotProduct_v_self_het N
  funext r
  change ∑ k, m.EperpHet N ω r k * WithLp.ofLp (m.stack.v N) k = 0
  have key : ∀ k, m.EperpHet N ω r k * WithLp.ofLp (m.stack.v N) k
      = m.stack.E N ω r k * WithLp.ofLp (m.stack.v N) k
        - m.eHet N ω r * (WithLp.ofLp (m.stack.v N) k * WithLp.ofLp (m.stack.v N) k) := by
    intro k
    change (m.stack.E N ω r k - m.eHet N ω r * WithLp.ofLp (m.stack.v N) k)
        * WithLp.ofLp (m.stack.v N) k = _
    ring
  rw [Finset.sum_congr rfl fun k _ => key k, Finset.sum_sub_distrib, ← Finset.mul_sum, hv,
    mul_one]
  exact sub_self _

/-- `X_w = q vᵀ + Σ^{1/2} E⊥`. -/
theorem stackW_X_eq_vecMulVec_add (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    (m.stackW w).X N ω
      = Matrix.vecMulVec (m.qHet w N ω) (WithLp.ofLp (m.stack.v N))
        + m.SigmaHalf w N * m.EperpHet N ω := by
  have hE := m.stack_E_eq_vecMulVec_add N ω
  have hX : (m.stackW w).X N ω
      = Matrix.vecMulVec (m.u0Het w N) (WithLp.ofLp (m.stack.v N))
        + m.SigmaHalf w N * m.stack.E N ω := by
    rw [← m.stackW_E_eq w N ω]
    change (m.stackW w).θ • Matrix.vecMulVec (WithLp.ofLp ((m.stackW w).u N))
        (WithLp.ofLp ((m.stackW w).v N)) + (m.stackW w).E N ω
      = Matrix.vecMulVec (m.u0Het w N) (WithLp.ofLp (m.stack.v N)) + (m.stackW w).E N ω
    congr 1
    ext r k
    simp only [u0Het, Matrix.smul_apply, Matrix.vecMulVec_apply, Pi.smul_apply, smul_eq_mul]
    exact (mul_assoc _ _ _).symm
  rw [hX, hE, Matrix.mul_add, qHet, Matrix.add_vecMulVec]
  have hsm : m.SigmaHalf w N
      * Matrix.vecMulVec (m.eHet N ω) (WithLp.ofLp (m.stack.v N))
      = Matrix.vecMulVec (m.SigmaHalf w N *ᵥ m.eHet N ω) (WithLp.ofLp (m.stack.v N)) := by
    ext r k
    simp [Matrix.mul_apply, Matrix.vecMulVec_apply, Matrix.mulVec, dotProduct,
      Finset.sum_mul, mul_assoc]
  rw [hsm]
  abel

/-- **(H0) for the column split.** `X_w X_wᵀ = W₀' + q qᵀ`. Pure algebra: it uses
`‖v‖ = 1` and `E⊥ v = 0` only. -/
theorem gram_eq_het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.stackW w).X N ω * ((m.stackW w).X N ω)ᵀ
      = m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω) := by
  set V : Fin (d N) → ℝ := WithLp.ofLp (m.stack.v N) with hV
  set G : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
    m.SigmaHalf w N * m.EperpHet N ω with hG
  have hGv : G *ᵥ V = 0 := by
    rw [hG, ← Matrix.mulVec_mulVec, EperpHet_mulVec_v, Matrix.mulVec_zero]
  have hF : (Gᵀ)ᵀ *ᵥ V = 0 := by rwa [Matrix.transpose_transpose]
  have key := gram_vecMulVec_add V (m.qHet w N ω) Gᵀ (m.dotProduct_v_self_het N) hF
  rw [Matrix.transpose_add, Matrix.transpose_vecMulVec, Matrix.transpose_transpose] at key
  have hX := m.stackW_X_eq_vecMulVec_add w N ω
  rw [hX, ← hG, Matrix.transpose_add, Matrix.transpose_vecMulVec, key, W0het, hG]
  congr 1
  rw [Matrix.transpose_mul, transpose_SigmaHalf, Matrix.mul_assoc, Matrix.mul_assoc,
    Matrix.mul_assoc]

theorem isHermitian_W0het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.W0het w N ω).IsHermitian := by
  have h : m.W0het w N ω
      = (m.SigmaHalf w N * m.EperpHet N ω) * (m.SigmaHalf w N * m.EperpHet N ω)ᵀ := by
    rw [W0het, Matrix.transpose_mul, transpose_SigmaHalf, Matrix.mul_assoc, Matrix.mul_assoc]
  rw [h]
  simpa using
    Matrix.isHermitian_mul_conjTranspose_self (m.SigmaHalf w N * m.EperpHet N ω)

theorem posSemidef_W0het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.W0het w N ω).PosSemidef := by
  have h : m.W0het w N ω
      = (m.SigmaHalf w N * m.EperpHet N ω) * ((m.SigmaHalf w N * m.EperpHet N ω))ᵀ := by
    rw [W0het, Matrix.transpose_mul, transpose_SigmaHalf, Matrix.mul_assoc, Matrix.mul_assoc]
  rw [h]
  simpa using
    Matrix.posSemidef_self_mul_conjTranspose (m.SigmaHalf w N * m.EperpHet N ω)

/-! ### The block law -/

/-- `√d e = Z_stack v`: the only scaling lemma of the column split. -/
theorem sqrt_smul_eHet (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Real.sqrt (d N) • m.eHet N ω = m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N) := by
  change Real.sqrt (d N) • (((Real.sqrt (d N))⁻¹ • m.stack.Z N ω)
      *ᵥ WithLp.ofLp (m.stack.v N)) = _
  rw [Matrix.smul_mulVec, smul_smul]
  rcases eq_or_lt_of_le (Nat.zero_le (d N)) with h | h
  · exact absurd h.symm (m.stack.hd N).ne'
  · rw [mul_inv_cancel₀ (Real.sqrt_ne_zero'.mpr (by exact_mod_cast h)), one_smul]

/-- `E⊥ = d^{-1/2} (Z - (Z v) vᵀ)`, the column analogue of `R0.Eperp_eq_smul`. -/
theorem EperpHet_eq_smul (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.EperpHet N ω = (Real.sqrt (d N))⁻¹ •
      (m.stack.Z N ω - Matrix.vecMulVec (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N))
        (WithLp.ofLp (m.stack.v N))) := by
  have he : m.eHet N ω
      = (Real.sqrt (d N))⁻¹ • (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N)) := by
    change ((Real.sqrt (d N))⁻¹ • m.stack.Z N ω) *ᵥ WithLp.ofLp (m.stack.v N) = _
    rw [Matrix.smul_mulVec]
  change m.stack.E N ω - Matrix.vecMulVec (m.eHet N ω) (WithLp.ofLp (m.stack.v N)) = _
  rw [he, Matrix.smul_vecMulVec, smul_sub]
  rfl

/-- `E⊥ E⊥ᵀ` on the unscaled noise: `E⊥ E⊥ᵀ = d⁻¹ (Z Zᵀ - (Z v) (Z v)ᵀ)`. The column
analogue of `R0.W0_eq_smul`. -/
theorem EperpHet_gram_eq_smul (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.EperpHet N ω * (m.EperpHet N ω)ᵀ = ((d N : ℝ))⁻¹ •
      (m.stack.Z N ω * (m.stack.Z N ω)ᵀ -
        Matrix.vecMulVec (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N))
          (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N))) := by
  have hd : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.stack.hd N
  have hsq : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt hd.le]
  rw [m.EperpHet_eq_smul N ω, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul,
    smul_smul, hsq, gram_sub_vecMulVec _ _ (m.dotProduct_v_self_het N)]

/-- `W₀'` written on a block of the noise. Deterministic bookkeeping: it turns the block
identity of `exists_block_hasLaw_het` into the form `RMT/R4.lean` consumes. -/
theorem W0het_eq_of_block (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {p : ℕ} (B : Matrix (Fin (∑ i, n i N)) (Fin p) ℝ)
    (hB : m.EperpHet N ω * (m.EperpHet N ω)ᵀ = ((d N : ℝ))⁻¹ • (B * Bᵀ)) :
    m.W0het w N ω
      = ((d N : ℝ))⁻¹ • ((m.SigmaHalf w N * B) * (m.SigmaHalf w N * B)ᵀ) := by
  have h1 : m.W0het w N ω
      = m.SigmaHalf w N * (m.EperpHet N ω * (m.EperpHet N ω)ᵀ) * m.SigmaHalf w N := by
    simp only [W0het, Matrix.mul_assoc]
  rw [h1, hB, Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  rw [Matrix.transpose_mul, transpose_SigmaHalf]
  simp only [Matrix.mul_assoc]

/-- **(H0') for the column split.** For `d N = p + 1` there is a `(∑ n_i) × p` block `B`
with `E⊥ E⊥ᵀ = d⁻¹ B Bᵀ`, and the pair `(Z_stack v, B)` has a product law, which carries
both marginals and their independence at once.

Proof (H0 row of `notes/archive/plan_heterolaw_A.md` section 4). The column split is the row split
of the transpose. R0 states its split on a `SpikedModel`, whose split index is the **row**
index, so the argument is not a transcription; `exists_rowBlock_measurePreserving` above is
R0's six steps in model-free form and it runs on `Z_stackᵀ`.

1. `m.stack_law hG : m.stack.GaussianNoise`, so `Z_stack ~ gaussianMatrix (∑ n_i) (d N)`.
2. `measurePreserving_transpose` gives `Z_stackᵀ ~ gaussianMatrix (d N) (∑ n_i)`.
3. `exists_rowBlock_measurePreserving` at `q = d N`, `x = v`, `r = ∑ n_i` completes `v` to
   an orthonormal basis of `ℝ^{d N}`, rotates, and splits row `0` of `U Z_stackᵀ` off the
   block `B₀` of rows `1, …, p`. Row `0` is `Z_stack v`. Put `B := B₀ᵀ`.
4. `measurePreserving_transpose` again on the second factor turns
   `gaussianMatrix p (∑ n_i)` into `gaussianMatrix (∑ n_i) p`.
5. `EperpHet_gram_eq_smul` is `E⊥ E⊥ᵀ = d⁻¹ (Z_stack Z_stackᵀ - (Z_stack v) (Z_stack v)ᵀ)`,
   and the Gram half of step 3 turns the right side into `d⁻¹ B₀ᵀ B₀ = d⁻¹ B Bᵀ`.
6. `p = 0` (that is `d N = 1`) needs no special case: `B` is the empty matrix and both
   sides are `0`. -/
theorem exists_block_hasLaw_het (m : MultiTableModel μ M n d) (N : ℕ)
    (hG : m.JointGaussianNoise) {p : ℕ} (hp : d N = p + 1) :
    ∃ B : Ω N → Matrix (Fin (∑ i, n i N)) (Fin p) ℝ,
      (∀ ω, m.EperpHet N ω * (m.EperpHet N ω)ᵀ = ((d N : ℝ))⁻¹ • (B ω * (B ω)ᵀ)) ∧
        HasLaw (fun ω => (m.stack.Z N ω *ᵥ WithLp.ofLp (m.stack.v N), B ω))
          ((Measure.pi fun _ : Fin (∑ i, n i N) => gaussianReal 0 1).prod
            (gaussianMatrix (∑ i, n i N) p)) (μ N) := by
  obtain ⟨Yfun, hgram, hmp⟩ :=
    exists_rowBlock_measurePreserving (∑ i, n i N) hp (m.stack.v N) (m.stack.hv N)
  refine ⟨fun ω => (Yfun ((m.stack.Z N ω)ᵀ))ᵀ, ?_, ?_⟩
  · intro ω
    have h := hgram ((m.stack.Z N ω)ᵀ)
    rw [Matrix.transpose_transpose] at h
    rw [m.EperpHet_gram_eq_smul N ω, Matrix.transpose_transpose, h]
  · have hZt : HasLaw (fun ω => (m.stack.Z N ω)ᵀ)
        (gaussianMatrix (d N) (∑ i, n i N)) (μ N) :=
      (measurePreserving_transpose (∑ i, n i N) (d N)).fun_comp_hasLaw (m.stack_law hG N)
    have hprod : MeasurePreserving
        (Prod.map (id : (Fin (∑ i, n i N) → ℝ) → (Fin (∑ i, n i N) → ℝ))
          (Matrix.transpose :
            Matrix (Fin p) (Fin (∑ i, n i N)) ℝ → Matrix (Fin (∑ i, n i N)) (Fin p) ℝ))
        ((Measure.pi fun _ : Fin (∑ i, n i N) => gaussianReal 0 1).prod
          (gaussianMatrix p (∑ i, n i N)))
        ((Measure.pi fun _ : Fin (∑ i, n i N) => gaussianReal 0 1).prod
          (gaussianMatrix (∑ i, n i N) p)) :=
      (MeasurePreserving.id _).prod (measurePreserving_transpose p (∑ i, n i N))
    have h := (hprod.comp hmp).fun_comp_hasLaw hZt
    simpa [Function.comp_def] using h

end MultiTableModel
end StackedSVD
