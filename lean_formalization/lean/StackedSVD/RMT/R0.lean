/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.GaussianMatrix
import StackedSVD.RMT

/-!
# Item R0: the splitting bridge

See `notes/archive/rmt_R0.md`. This file supplies the two facts that
`notes/archive/SINGLE_TABLE_PLAN.md` §2 states without proof, and that no other item proves.

With `E = d^{-1/2} Z`, `X = θ u vᵀ + E`, `‖u‖ = ‖v‖ = 1`, put

```
g   = Eᵀ u,   E⊥ = E - u gᵀ,   q = θ v + g,   W₀ = E⊥ᵀ E⊥ .
```

* **(H0)** `gram_eq`: `Xᵀ X = W₀ + q qᵀ`. Deterministic, it uses `‖u‖ = 1` only.
* **(H0')** `exists_block_hasLaw`: for `n = p + 1` there is a `p × d` matrix `B` with
  `W₀ = d⁻¹ Bᵀ B`, and the pair `(Zᵀ u, B)` has the product law
  `(⨂_{j<d} N(0,1)) ⊗ gaussianMatrix p d`. One `HasLaw` on a product measure carries both
  marginals and the independence.
* **(H0'')** `tendsto_pred_ratio`: `(n - 1)/d → c` when `n/d → c` and `d → ∞`.

The route for (H0') is the one of the note: complete `u` to an orthonormal basis `b` of `ℝ^n`
with `b 0 = u`, rotate by the matrix `U` whose row `i` is `b i` (the law is invariant by
`gaussianMatrix_map_mul`), and split row `0` of `U Z` off the block of rows `1, …, p`
(`measurePreserving_rowSplit`). `p = 0`, that is `n = 1`, needs no special case: the block is
the empty matrix, `Bᵀ B = 0` and `W₀ = 0`.

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/R0.lean` exit 0, 0 `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### Deterministic algebra of a rank-one split -/

section GramAlgebra

variable {r c : ℕ}

/-- If `F` kills `x` on the left then the Gram matrix of `x yᵀ + F` splits into the Gram
matrix of `F` and the rank-one matrix `y yᵀ`. This is (H0) with no model in sight. -/
theorem gram_vecMulVec_add (x : Fin r → ℝ) (y : Fin c → ℝ) (F : Matrix (Fin r) (Fin c) ℝ)
    (hx : x ⬝ᵥ x = 1) (hF : Fᵀ *ᵥ x = 0) :
    (Matrix.vecMulVec x y + F)ᵀ * (Matrix.vecMulVec x y + F)
      = Fᵀ * F + Matrix.vecMulVec y y := by
  have hxx : ∑ k, x k * x k = 1 := hx
  have hF' : ∀ j, ∑ k, F k j * x k = 0 := by
    intro j
    have h := congrFun hF j
    simpa [Matrix.mulVec, dotProduct, Matrix.transpose_apply] using h
  ext j j'
  change ∑ k, (x k * y j + F k j) * (x k * y j' + F k j')
      = (∑ k, F k j * F k j') + y j * y j'
  have key : ∑ k, (x k * y j + F k j) * (x k * y j' + F k j')
      = ∑ k, ((y j * y j') * (x k * x k) + y j' * (F k j * x k)
        + y j * (F k j' * x k) + F k j * F k j') :=
    Finset.sum_congr rfl fun k _ => by ring
  rw [key, Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, hxx, hF' j, hF' j']
  ring

/-- The complement `M - x (Mᵀ x)ᵀ` is orthogonal to `x`. This is R0a.1. -/
theorem transpose_sub_vecMulVec_mulVec (x : Fin r → ℝ) (M : Matrix (Fin r) (Fin c) ℝ)
    (hx : x ⬝ᵥ x = 1) :
    (M - Matrix.vecMulVec x (Mᵀ *ᵥ x))ᵀ *ᵥ x = 0 := by
  have hxx : ∑ k, x k * x k = 1 := hx
  funext j
  change ∑ k, (M k j - x k * (∑ l, M l j * x l)) * x k = 0
  have key : ∑ k, (M k j - x k * (∑ l, M l j * x l)) * x k
      = ∑ k, (M k j * x k - (∑ l, M l j * x l) * (x k * x k)) :=
    Finset.sum_congr rfl fun k _ => by ring
  rw [key, Finset.sum_sub_distrib, ← Finset.mul_sum, hxx, mul_one, sub_self]

end GramAlgebra

/-! ### Splitting the first row off a canonical Gaussian matrix -/

/-- The row split on plain pi types. `Matrix` is a `def`, so instance search does not find the
product sigma-algebra on it; the `Matrix` form below is this statement, transferred by `exact`
(the same device as `measurePreserving_stackPi` in `StackSVD.lean`). -/
private theorem measurePreserving_rowSplitPi (p d : ℕ) :
    MeasurePreserving
      (fun Z : Fin (p + 1) → Fin d → ℝ => (Z 0, fun (i : Fin p) (j : Fin d) => Z i.succ j))
      (Measure.pi fun _ : Fin (p + 1) => Measure.pi fun _ : Fin d => gaussianReal 0 1)
      ((Measure.pi fun _ : Fin d => gaussianReal 0 1).prod
        (Measure.pi fun _ : Fin p => Measure.pi fun _ : Fin d => gaussianReal 0 1)) :=
  measurePreserving_piFinSuccAbove
    (fun _ : Fin (p + 1) => Measure.pi fun _ : Fin d => gaussianReal 0 1) 0

/-- **Splitting the first row off a canonical Gaussian matrix gives a product law.** -/
theorem measurePreserving_rowSplit (p d : ℕ) :
    MeasurePreserving
      (fun Z : Matrix (Fin (p + 1)) (Fin d) ℝ =>
        (Z 0, Matrix.of fun (i : Fin p) (j : Fin d) => Z i.succ j))
      (gaussianMatrix (p + 1) d)
      ((Measure.pi fun _ : Fin d => gaussianReal 0 1).prod (gaussianMatrix p d)) :=
  measurePreserving_rowSplitPi p d

/-! ### Recasting the row index along `n = p + 1` -/

/-- Row reindexing along an equality of dimensions, on plain pi types. -/
private theorem measurePreserving_castRowPi {q p d : ℕ} (h : q = p) :
    MeasurePreserving
      (fun (Z : Fin q → Fin d → ℝ) (i : Fin p) => Z (Fin.cast h.symm i))
      (Measure.pi fun _ : Fin q => Measure.pi fun _ : Fin d => gaussianReal 0 1)
      (Measure.pi fun _ : Fin p => Measure.pi fun _ : Fin d => gaussianReal 0 1) := by
  subst h
  exact MeasurePreserving.id _

/-- Row reindexing along `q = p` preserves the canonical Gaussian matrix law. -/
theorem measurePreserving_castRow {q p d : ℕ} (h : q = p) :
    MeasurePreserving
      (fun Z : Matrix (Fin q) (Fin d) ℝ =>
        Matrix.of fun (i : Fin p) (j : Fin d) => Z (Fin.cast h.symm i) j)
      (gaussianMatrix q d) (gaussianMatrix p d) :=
  measurePreserving_castRowPi h

/-- Row reindexing does not change the Gram matrix. -/
theorem gram_castRow {q p d : ℕ} (h : q = p) (Z : Matrix (Fin q) (Fin d) ℝ) :
    (Matrix.of fun (i : Fin p) (j : Fin d) => Z (Fin.cast h.symm i) j)ᵀ *
        (Matrix.of fun (i : Fin p) (j : Fin d) => Z (Fin.cast h.symm i) j)
      = Zᵀ * Z := by
  subst h
  rfl

/-- `√k ≠ 0` for a positive natural `k`. -/
private theorem sqrt_natCast_ne_zero {k : ℕ} (hk : 0 < k) : Real.sqrt (k : ℝ) ≠ 0 := by
  have hpos : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  exact Real.sqrt_ne_zero'.2 hpos

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}
variable (m : SpikedModel μ n d) (N : ℕ)

/-! ### The four split objects -/

/-- `g = Eᵀ u`, the component of the noise along `u`. -/
noncomputable def gvec (ω : Ω N) : Fin (d N) → ℝ :=
  (m.E N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)

/-- `E⊥ = E - u gᵀ`. -/
noncomputable def Eperp (ω : Ω N) : Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  m.E N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω)

/-- `W₀ = E⊥ᵀ E⊥`, the Wishart part of `Xᵀ X`. -/
noncomputable def W0 (ω : Ω N) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.Eperp N ω)ᵀ * m.Eperp N ω

/-- `q = θ v + g`, the rank-one part of `Xᵀ X`. -/
noncomputable def qvec (ω : Ω N) : Fin (d N) → ℝ :=
  m.θ • WithLp.ofLp (m.v N) + m.gvec N ω

/-! ### Small facts about the unit vector `u` -/

theorem dotProduct_u_self : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) = 1 := by
  have h : ⟪m.u N, m.u N⟫_ℝ = WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) :=
    inner_euclidean_eq_dotProduct (m.u N) (m.u N)
  rw [real_inner_self_eq_norm_sq, m.hu N] at h
  simpa using h.symm

/-! ### R0a: the deterministic identities -/

/-- `g = d^{-1/2} (Zᵀ u)`. -/
theorem gvec_eq_smul (ω : Ω N) :
    m.gvec N ω = (Real.sqrt (d N))⁻¹ • ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) := by
  change ((Real.sqrt (d N))⁻¹ • m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N) = _
  rw [Matrix.transpose_smul, Matrix.smul_mulVec]

/-- The only scaling lemma: `√d g = Zᵀ u`. -/
theorem sqrt_smul_gvec (ω : Ω N) :
    Real.sqrt (d N) • m.gvec N ω = (m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N) := by
  rw [m.gvec_eq_smul N ω, smul_smul, mul_inv_cancel₀ (sqrt_natCast_ne_zero (m.hd N)),
    one_smul]

/-- `E⊥ = d^{-1/2} (Z - u (Zᵀ u)ᵀ)`. -/
theorem Eperp_eq_smul (ω : Ω N) :
    m.Eperp N ω = (Real.sqrt (d N))⁻¹ •
      (m.Z N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N))
        ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) := by
  change m.E N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω) = _
  rw [m.gvec_eq_smul N ω, Matrix.vecMulVec_smul, smul_sub]
  rfl

/-- **R0a.1.** `E⊥ᵀ u = 0`. -/
theorem transpose_Eperp_mulVec_u (ω : Ω N) :
    (m.Eperp N ω)ᵀ *ᵥ WithLp.ofLp (m.u N) = 0 := by
  change (m.E N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω))ᵀ *ᵥ _ = _
  exact transpose_sub_vecMulVec_mulVec _ _ (m.dotProduct_u_self N)

/-- **R0a.2.** `E = u gᵀ + E⊥`. -/
theorem E_eq_vecMulVec_add (ω : Ω N) :
    m.E N ω = Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω) + m.Eperp N ω := by
  change _ = _ + (m.E N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω))
  abel

/-- **R0a.2.** `X = u qᵀ + E⊥`. -/
theorem X_eq_vecMulVec_add (ω : Ω N) :
    m.X N ω = Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.qvec N ω) + m.Eperp N ω := by
  have hq : Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.qvec N ω)
      = m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N))
        + Matrix.vecMulVec (WithLp.ofLp (m.u N)) (m.gvec N ω) := by
    change Matrix.vecMulVec (WithLp.ofLp (m.u N))
        (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
    rw [Matrix.vecMulVec_add, Matrix.vecMulVec_smul]
  rw [hq, add_assoc, ← m.E_eq_vecMulVec_add N ω]
  rfl

/-- **(H0).** The splitting identity R5 instantiates R4 with. -/
theorem gram_eq (ω : Ω N) :
    (m.X N ω)ᵀ * m.X N ω = m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω) := by
  rw [m.X_eq_vecMulVec_add N ω]
  exact gram_vecMulVec_add _ _ _ (m.dotProduct_u_self N) (m.transpose_Eperp_mulVec_u N ω)

theorem isHermitian_W0 (ω : Ω N) : (m.W0 N ω).IsHermitian :=
  isHermitian_transpose_mul_self (m.Eperp N ω)

theorem posSemidef_W0 (ω : Ω N) : (m.W0 N ω).PosSemidef := by
  change ((m.Eperp N ω)ᵀ * m.Eperp N ω).PosSemidef
  simpa using Matrix.posSemidef_conjTranspose_mul_self (m.Eperp N ω)

/-- `W₀` written on the unscaled noise: `W₀ = d⁻¹ (Zᵀ Z - w wᵀ)` with `w = Zᵀ u`. -/
theorem W0_eq_smul (ω : Ω N) :
    m.W0 N ω = ((d N : ℝ))⁻¹ •
      ((m.Z N ω)ᵀ * m.Z N ω -
        Matrix.vecMulVec ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))
          ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) := by
  set x : Fin (n N) → ℝ := WithLp.ofLp (m.u N) with hxdef
  set w : Fin (d N) → ℝ := (m.Z N ω)ᵀ *ᵥ x with hwdef
  set F : Matrix (Fin (n N)) (Fin (d N)) ℝ := m.Z N ω - Matrix.vecMulVec x w with hFdef
  have hx : x ⬝ᵥ x = 1 := m.dotProduct_u_self N
  have hF : Fᵀ *ᵥ x = 0 := transpose_sub_vecMulVec_mulVec x (m.Z N ω) hx
  have hZ : Matrix.vecMulVec x w + F = m.Z N ω := by rw [hFdef]; abel
  have hgram : (m.Z N ω)ᵀ * m.Z N ω = Fᵀ * F + Matrix.vecMulVec w w := by
    rw [← hZ]; exact gram_vecMulVec_add x w F hx hF
  have hFF : Fᵀ * F = (m.Z N ω)ᵀ * m.Z N ω - Matrix.vecMulVec w w := by
    rw [hgram]; abel
  have hsq : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (by positivity)]
  change (m.Eperp N ω)ᵀ * m.Eperp N ω = _
  rw [m.Eperp_eq_smul N ω, ← hFdef, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul,
    smul_smul, hsq, hFF]

/-! ### R0b: the laws -/

/-- **(H0').** There is a `p × d` Gaussian block `B`, independent of `g`, with
`W₀ = d⁻¹ Bᵀ B`. The `HasLaw` on a product measure carries both marginals and the
independence. -/
theorem exists_block_hasLaw (hG : m.GaussianNoise) {p : ℕ} (hp : n N = p + 1) :
    ∃ B : Ω N → Matrix (Fin p) (Fin (d N)) ℝ,
      (∀ ω, m.W0 N ω = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)) ∧
        HasLaw (fun ω => ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B ω))
          ((Measure.pi fun _ : Fin (d N) => gaussianReal 0 1).prod
            (gaussianMatrix p (d N))) (μ N) := by
  classical
  -- 1. an orthonormal basis of `ℝ^{n N}`, indexed by `Fin (p + 1)`, whose first vector is `u`.
  have hcard : Module.finrank ℝ (EuclideanSpace ℝ (Fin (n N))) = Fintype.card (Fin (p + 1)) := by
    rw [finrank_euclideanSpace_fin, Fintype.card_fin, hp]
  have hortho : Orthonormal ℝ
      (({0} : Set (Fin (p + 1))).domRestrict (fun _ : Fin (p + 1) => m.u N)) := by
    refine ⟨fun i => m.hu N, ?_⟩
    intro i j hij
    exact absurd (Subtype.ext ((i.2 : (i : Fin (p + 1)) = 0).trans
      (j.2 : (j : Fin (p + 1)) = 0).symm)) hij
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq hcard
  have hb0 : WithLp.ofLp (b 0) = WithLp.ofLp (m.u N) := by
    rw [hb 0 rfl]
  -- 2. the rotation matrix: row `i` is `b i`, read through `Fin (p + 1) ≃ Fin (n N)`.
  set e : Fin (p + 1) → Fin (n N) := Fin.cast hp.symm with hedef
  set U : Matrix (Fin (p + 1)) (Fin (p + 1)) ℝ :=
    Matrix.of fun i i' => WithLp.ofLp (b i) (e i') with hUdef
  have hsum : ∀ f : Fin (n N) → ℝ, ∑ i' : Fin (p + 1), f (e i') = ∑ k : Fin (n N), f k := by
    intro f
    exact Equiv.sum_comp (finCongr hp.symm) f
  have hUUt : U * Uᵀ = 1 := by
    ext i i'
    have hb' := (orthonormal_iff_ite (𝕜 := ℝ)).1 b.orthonormal i i'
    rw [inner_euclidean_eq_dotProduct] at hb'
    rw [Matrix.one_apply, ← hb']
    change ∑ i'' : Fin (p + 1),
        WithLp.ofLp (b i) (e i'') * WithLp.ofLp (b i') (e i'') = _
    exact hsum fun k => WithLp.ofLp (b i) k * WithLp.ofLp (b i') k
  have hU : Uᵀ * U = 1 := mul_eq_one_comm.1 hUUt
  -- 3. the composite map `Z ↦ (row 0 of U Z, block of rows 1..p of U Z)`.
  set castR : Matrix (Fin (n N)) (Fin (d N)) ℝ → Matrix (Fin (p + 1)) (Fin (d N)) ℝ :=
    fun Z => Matrix.of fun (i : Fin (p + 1)) (j : Fin (d N)) => Z (e i) j with hcastdef
  set Y : Matrix (Fin (n N)) (Fin (d N)) ℝ → Matrix (Fin (p + 1)) (Fin (d N)) ℝ :=
    fun Z => U * castR Z with hYdef
  have hmp : MeasurePreserving (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      (Y Z 0, Matrix.of fun (i : Fin p) (j : Fin (d N)) => Y Z i.succ j))
      (gaussianMatrix (n N) (d N))
      ((Measure.pi fun _ : Fin (d N) => gaussianReal 0 1).prod (gaussianMatrix p (d N))) :=
    (measurePreserving_rowSplit p (d N)).comp
      ((measurePreserving_mul_left hU (d N)).comp (measurePreserving_castRow hp))
  -- 4. the row of `U Z` in terms of `b`.
  have hrow : ∀ (Z : Matrix (Fin (n N)) (Fin (d N)) ℝ) (i : Fin (p + 1)) (j : Fin (d N)),
      Y Z i j = ∑ k : Fin (n N), WithLp.ofLp (b i) k * Z k j := by
    intro Z i j
    change ∑ i' : Fin (p + 1), WithLp.ofLp (b i) (e i') * Z (e i') j = _
    exact hsum fun k => WithLp.ofLp (b i) k * Z k j
  have hrow0 : ∀ Z : Matrix (Fin (n N)) (Fin (d N)) ℝ,
      Y Z 0 = Zᵀ *ᵥ WithLp.ofLp (m.u N) := by
    intro Z
    funext j
    rw [hrow Z 0 j, hb0]
    change _ = ∑ k, Z k j * WithLp.ofLp (m.u N) k
    exact Finset.sum_congr rfl fun k _ => mul_comm _ _
  -- 5. the Gram matrix of `U Z` splits into row `0` and the block.
  have hgramY : ∀ Z : Matrix (Fin (n N)) (Fin (d N)) ℝ,
      (Matrix.of fun (i : Fin p) (j : Fin (d N)) => Y Z i.succ j)ᵀ *
          (Matrix.of fun (i : Fin p) (j : Fin (d N)) => Y Z i.succ j)
        = Zᵀ * Z - Matrix.vecMulVec (Zᵀ *ᵥ WithLp.ofLp (m.u N))
            (Zᵀ *ᵥ WithLp.ofLp (m.u N)) := by
    intro Z
    have hYY : (Y Z)ᵀ * Y Z = Zᵀ * Z := by
      rw [hYdef]
      simp only
      rw [Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Uᵀ U (castR Z), hU,
        Matrix.one_mul]
      exact gram_castRow hp Z
    ext j j'
    have hsplit : ∑ i : Fin (p + 1), Y Z i j * Y Z i j'
        = Y Z 0 j * Y Z 0 j' + ∑ i : Fin p, Y Z i.succ j * Y Z i.succ j' :=
      Fin.sum_univ_succ (fun i : Fin (p + 1) => Y Z i j * Y Z i j')
    have hL : ((Y Z)ᵀ * Y Z) j j' = ∑ i : Fin (p + 1), Y Z i j * Y Z i j' := rfl
    have hR : (Zᵀ * Z) j j' = ∑ i : Fin (p + 1), Y Z i j * Y Z i j' := by
      rw [← hL, hYY]
    change ∑ i : Fin p, Y Z i.succ j * Y Z i.succ j' = (Zᵀ * Z) j j' - _
    rw [hR, hsplit, hrow0 Z]
    change _ = _ + _ - (Zᵀ *ᵥ WithLp.ofLp (m.u N)) j * (Zᵀ *ᵥ WithLp.ofLp (m.u N)) j'
    ring
  -- 6. the block, and the two conclusions.
  refine ⟨fun ω => Matrix.of fun (i : Fin p) (j : Fin (d N)) => Y (m.Z N ω) i.succ j, ?_, ?_⟩
  · intro ω
    rw [m.W0_eq_smul N ω, hgramY (m.Z N ω)]
  · have h := hmp.fun_comp_hasLaw (hG N)
    simpa only [hrow0] using h

/-- **(H0'), first half.** `√d g = Zᵀ u` is a standard Gaussian vector of `ℝ^d`. -/
theorem hasLaw_gvec (hG : m.GaussianNoise) :
    HasLaw (fun ω => (m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))
      (Measure.pi fun _ : Fin (d N) => gaussianReal 0 1) (μ N) := by
  obtain ⟨p, hp⟩ := Nat.exists_eq_succ_of_ne_zero (m.hn N).ne'
  obtain ⟨B, -, hlaw⟩ := m.exists_block_hasLaw N hG hp
  exact measurePreserving_fst.fun_comp_hasLaw hlaw

/-- The block alone has the canonical Gaussian matrix law. -/
theorem exists_block_hasLaw_snd (hG : m.GaussianNoise) {p : ℕ} (hp : n N = p + 1) :
    ∃ B : Ω N → Matrix (Fin p) (Fin (d N)) ℝ,
      (∀ ω, m.W0 N ω = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)) ∧
        HasLaw B (gaussianMatrix p (d N)) (μ N) := by
  obtain ⟨B, hW, hlaw⟩ := m.exists_block_hasLaw N hG hp
  exact ⟨B, hW, measurePreserving_snd.fun_comp_hasLaw hlaw⟩

/-- Independence in the form R2 conditions on. Corollary of `exists_block_hasLaw`. -/
theorem exists_block_indepFun [IsProbabilityMeasure (μ N)] (hG : m.GaussianNoise)
    {p : ℕ} (hp : n N = p + 1) :
    ∃ B : Ω N → Matrix (Fin p) (Fin (d N)) ℝ,
      (∀ ω, m.W0 N ω = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)) ∧
        IndepFun (fun ω => (m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) B (μ N) := by
  obtain ⟨B, hW, hlaw⟩ := m.exists_block_hasLaw N hG hp
  refine ⟨B, hW, ?_⟩
  exact (indepFun_iff_hasLaw_prodMk_prod (measurePreserving_fst.fun_comp_hasLaw hlaw)
    (measurePreserving_snd.fun_comp_hasLaw hlaw)).2 hlaw

/-! ### R0c: the aspect ratio of `W₀` -/

/-- **(H0'').** `W₀` has the same limiting aspect ratio as `X`. -/
theorem tendsto_pred_ratio {c : ℝ} (hreg : m.Regime c) :
    Tendsto (fun N => ((n N - 1 : ℕ) : ℝ) / (d N : ℝ)) atTop (𝓝 c) := by
  have hd : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hreg.2.1
  have hinv : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := hd.inv_tendsto_atTop
  have heq : ∀ N : ℕ, ((n N - 1 : ℕ) : ℝ) / (d N : ℝ)
      = (n N : ℝ) / (d N : ℝ) - ((d N : ℝ))⁻¹ := by
    intro N
    rw [Nat.cast_sub (m.hn N), Nat.cast_one, sub_div, one_div]
  simp only [heq]
  simpa using hreg.2.2.sub hinv

end SpikedModel

end StackedSVD

