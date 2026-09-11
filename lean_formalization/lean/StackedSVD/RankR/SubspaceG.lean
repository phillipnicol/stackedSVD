/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.General
import StackedSVD.RankR.Subspace

/-!
# `prop:stacksvd_subspace` at general `r_i`

STATUS 2026-09-01: proved, 0 `sorry` (tasks PA and PB of `notes/archive/rankr_plan_B.md` section 4).
The plan is the review note; the audit is `notes/archive/audit_rankr_plan_B_2026-09-01.md` and its
numeric section checks every identity of this file (rows P1 to P7, all PASS).

`RankR/Subspace.lean` proves `prop:stacksvd_subspace` (`main_paper.tex:834`) when every table
carries one spike. This file removes that restriction. Table `i` of an `UnalignedModelR` has
`r_i` spikes,

```
X_i = U_i Θ_i (V R_i)ᵀ + E_i,   R_i ∈ O(r, r_i),  Θ_i = diag(θ_i1, …, θ_i r_i),  U_i ∈ O(n_i, r_i),
```

so the unit-weight stack is again a rank-`r` spiked matrix (`main_paper.tex:818`)

```
X_stack = [U_1 Θ_1 R_1ᵀ Vᵀ + E_1; …; U_M Θ_M R_Mᵀ Vᵀ + E_M] = A Vᵀ + E_stack
```

with population Gram (`main_paper.tex:824`, `:826`)

```
E[X_stackᵀ] E[X_stack] = ∑_i ∑_j θ_ij² (V R_i)_j (V R_i)_jᵀ = V C Vᵀ,  C = ∑_i R_i Θ_i² R_iᵀ.
```

Every declaration of `RankR/Subspace.lean` between `Cmat` and `prop_stacksvd_subspace` has a
twin here with the suffix `G`, and the route of the proposition does not change.

## Why the flatten bridge does not apply

`RankR/Flatten.lean` reads a general-`r_i` block object as the `r_i = 1` object at table count
`r̃ = ∑_i r_i`. That bridge closes `CBlock` and `limitStackRG` at the scalar level, but it does
**not** close the stack bridge of section 2b: the `r_i` spikes of table `i` share one noise
matrix and one row block, so a general-`r_i` model is not a stack of `r̃` rank-one tables
(choice 8 of the plan). Sections 2 and 2b are therefore proved directly, with `Uᵀ U = 1` read
as an `r_i × r_i` identity instead of `‖u_i‖ = 1`.

## Content

1. `CBlock θ R = (B_θ)ᵀ B_θ`, the core matrix `C = ∑_i R_i Θ_i² R_iᵀ` at general `r_i`
   (`cblock_eq_sum` is the paper's form), and `limitStackRG`, the limit.
2. `UnalignedModelR.stackXG`, the unit-weight stack, its Gram `stackGramG`, and `coreG`.
2b. The stack bridge: `signalFactorG` (`A`), `stackEG`, `signalPartG` (`= A Vᵀ`), `stackX_eqG`,
   `signalFactorG_transpose_mul_self`, `signalGram_eq_signalPartG`, `signalGram_eqG`,
   `signalGram_eq_sum_spikeVecG`.
3. `perfStackRG`, the top-`r` subspace performance, and `SubspaceLawG`, the one new hypothesis
   structure.
4. `prop_stacksvd_subspace_general`, the Layer 1 implication.

## The performance, and the hypothesis structure

`perfStackRG` is `tr(Vᵀ P_top(X_stackᵀ X_stack, r) V) = ∑_{k<r} ‖P_top (V e_k)‖²`, the same
choice as `perfStackR` at `r_i = 1`. It equals the paper's `‖V̂_stacksvdᵀ V‖_F²` whenever no
eigenvalue of the Gram ties across the index `r` boundary; `V̂_stacksvd` is not defined in
Lean, so that identity is a doc-comment (choice 5 of the plan). It was measured to 5.3e-15 over
288 draws at `M = 2`, `r_i = (2, 1)` (row P5 of the audit).

`SubspaceLawG` has one field, as `SubspaceLaw` does. It carries no eigengap of `C` and no
distinct-spike condition, which is what the paper claims at `main_paper.tex:842`; the audit of
2026-08-30 (item 4.1) shows that a `topGap` field would be false, not almost surely false, at
any `N` with `∑_i n_i N < r ≤ d N`. The numeric check runs a cell with `spec(C) = [3.25, 3.25]`
and the bias still falls with `d` (row P7).

`prop_stacksvd_subspace_general` moved to `RankR/SubspaceMain.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 0. Helper lemmas

Private copies of five helpers of `RankR/Subspace.lean`. They are `private` there, so they
cannot be imported; the statements and the proofs are unchanged. `sum_stack_index` is public
there and is used directly. -/

/-- Orthogonality of the eigenbasis, entrywise: `Q Qᵀ = I` for the matrix `Q` whose columns
are `Matrix.IsHermitian.eigenvectorBasis`. -/
private theorem sum_eigenvectorBasis_mul {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) (k l : Fin p) :
    ∑ j : Fin p, hA.eigenvectorBasis j k * hA.eigenvectorBasis j l
      = if k = l then 1 else 0 := by
  have hU : (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) *
      star (hA.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ) = 1 :=
    Unitary.coe_mul_star_self _
  have h := congrFun (congrFun hU k) l
  simpa [Matrix.mul_apply, Matrix.star_apply, Matrix.IsHermitian.eigenvectorUnitary_apply,
    Matrix.one_apply] using h

/-- A finite sum of vectors of `EuclideanSpace` is summed coordinate by coordinate. -/
private theorem euclid_sum_apply {ι : Type*} {D : ℕ} (s : Finset ι)
    (f : ι → EuclideanSpace ℝ (Fin D)) (l : Fin D) :
    (∑ i ∈ s, f i) l = ∑ i ∈ s, f i l := by
  classical
  refine Finset.induction_on s ?_ ?_
  · simp
  · intro a t ha ih
    rw [Finset.sum_insert ha, Finset.sum_insert ha, PiLp.add_apply, ih]

/-- A real symmetric matrix is the sum of its rank-one spectral pieces, `A = ∑_j λ_j q_j q_jᵀ`.
This is `Matrix.IsHermitian.spectral_theorem` read entrywise. -/
private theorem hermitian_eq_sum_rankOne {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) :
    A = ∑ j : Fin p, hA.eigenvalues j •
      Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis j))
        (WithLp.ofLp (hA.eigenvectorBasis j)) := by
  ext k l
  rw [Matrix.sum_apply]
  conv_lhs => rw [hA.spectral_theorem]
  simp only [Unitary.conjStarAlgAut_apply, Matrix.mul_apply, Matrix.diagonal_apply,
    Matrix.star_apply, Matrix.smul_apply, Matrix.vecMulVec_apply, smul_eq_mul,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial, Function.comp_apply,
    RCLike.ofReal_real_eq_id, id_eq, mul_ite, mul_zero,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- Conjugating a rank-one matrix: `V (x yᵀ) Vᵀ = (V x)(V y)ᵀ`. -/
private theorem mul_vecMulVec_mul_transpose {D R : ℕ} (V : Matrix (Fin D) (Fin R) ℝ)
    (x y : Fin R → ℝ) :
    V * Matrix.vecMulVec x y * Vᵀ = Matrix.vecMulVec (V *ᵥ x) (V *ᵥ y) := by
  rw [Matrix.mul_vecMulVec, Matrix.vecMulVec_mul, Matrix.vecMul_transpose]

/-- One entry of `X = U Θ Vᵀ + E` at rank `rk`: the signal part is the sum over the `rk`
spikes. The rank-one twin is `SpikedModel.X` read through `Matrix.vecMulVec_apply`. -/
private theorem spikedR_X_apply {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ} {rk : ℕ} (t : SpikedModelR μ n d rk) (N : ℕ)
    (ω : Ω N) (a : Fin (n N)) (k : Fin (d N)) :
    t.X N ω a k = (∑ j, t.U N a j * t.θ j * t.V N k j) + t.E N ω a k := by
  simp only [SpikedModelR.X, Matrix.add_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Matrix.diagonal_apply, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-! ### 1. The core matrix `C` and the limit -/

/-- The core matrix `C = ∑_i R_i Θ_i² R_iᵀ` (`main_paper.tex:826`) at general `r_i`, written
as `B_θᵀ B_θ` with `B_θ = BBlock θ R` of `RankR/General.lean`. The paper's sum form is
`cblock_eq_sum`. The rank-one mirror is `Cmat` (`RankR/Subspace.lean:180`); the argument `θ`
is the family of signal strengths, not the family `β`. -/
noncomputable def CBlock {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : Matrix (Fin r) (Fin r) ℝ :=
  (BBlock θ R)ᵀ * BBlock θ R

theorem isHermitian_CBlock {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : (CBlock θ R).IsHermitian :=
  isHermitian_transpose_mul_self (BBlock θ R)

/-- One entry of the core matrix, as a double sum over the block index `(i, j)`. -/
theorem cblock_apply {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (k l : Fin r) :
    CBlock θ R k l = ∑ i, ∑ j, θ i j ^ 2 * (R i k j * R i l j) := by
  have hstep := sum_stack_index (ν := rk) (fun i j => θ i j ^ 2 * (R i k j * R i l j))
  simp only [CBlock, Matrix.mul_apply, Matrix.transpose_apply, BBlock, Matrix.of_apply,
    betaFlat, blk]
  rw [← hstep]
  exact Finset.sum_congr rfl fun p _ => by ring

/-- The paper's form of the core matrix, `C = ∑_i R_i Θ_i² R_iᵀ` (`main_paper.tex:826`).
Row `j` of `Θ_i R_iᵀ` is `θ_ij (R_i)_jᵀ`, so the `(i, j)` term is `θ_ij²` times the rank-one
matrix built from column `j` of `R_i`. -/
theorem cblock_eq_sum {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    CBlock θ R = ∑ i, ∑ j, θ i j ^ 2 •
      Matrix.vecMulVec (fun k => R i k j) (fun k => R i k j) := by
  ext k l
  rw [cblock_apply, Matrix.sum_apply]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.sum_apply]
  exact Finset.sum_congr rfl fun j _ => by simp [Matrix.vecMulVec_apply]

/-- The limit of `prop:stacksvd_subspace` (`main_paper.tex:836`) at general `r_i`: one
rank-one performance per spike of the core matrix, at signal `√(λ_j(C))` and aspect ratio
`‖c‖₁ = ∑_i c_i`. The rank-one mirror is `limitStackR`.

`c` is indexed by the **table**, not by the spike: the stack has `∑_i n_i` rows, so its aspect
ratio is `∑_i c_i`. A `c` flattened over `Fin (rtot rk)` would count table `i` `r_i` times
(choice 4 of the plan). -/
noncomputable def limitStackRG {M r : ℕ} {rk : Fin M → ℕ} (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c : Fin M → ℝ) : ℝ :=
  ∑ j : Fin r, betaSq (Real.sqrt ((isHermitian_CBlock θ R).eigenvalues j)) (∑ i, c i)

/-! ### 2. The unit-weight stack of an `UnalignedModelR` -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The unit-weight stack `[X_1; …; X_M]`, with rows reindexed from the block index
`(i : Fin M) × Fin (n i N)` by `finSigmaFinEquiv`, as `UnalignedModel.stackX` does. -/
noncomputable def stackXG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).X N ω p.2 k)

/-- Block `i` of the stack is `X_i`. -/
theorem stackXG_apply (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) (i : Fin M)
    (a : Fin (n i N)) (k : Fin (d N)) :
    m.stackXG N ω (finSigmaFinEquiv ⟨i, a⟩) k = (m.tbl i).X N ω a k := by
  simp [stackXG, Matrix.reindex_apply, Matrix.submatrix_apply]

/-- The stack read at an arbitrary row index. This is the form the sums below consume. -/
theorem stackXG_apply' (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackXG N ω q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).X N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- Gram matrix of the stack. -/
noncomputable def stackGramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.stackXG N ω)ᵀ * m.stackXG N ω

theorem isHermitian_stackGramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    (m.stackGramG N ω).IsHermitian :=
  isHermitian_transpose_mul_self (m.stackXG N ω)

/-- `X_stackᵀ X_stack = ∑_i X_iᵀ X_i` (`main_paper.tex:823`). -/
theorem stackGramG_eq_sum (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.stackGramG N ω = ∑ i, ((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω := by
  ext k l
  simp only [stackGramG, Matrix.mul_apply, Matrix.transpose_apply, Matrix.sum_apply,
    stackXG_apply']
  exact sum_stack_index (ν := fun i => n i N)
    (fun i a => (m.tbl i).X N ω a k * (m.tbl i).X N ω a l)

/-- The core matrix of the model, `C = ∑_i R_i Θ_i² R_iᵀ`. -/
noncomputable def coreG (m : UnalignedModelR μ M n d r rk) : Matrix (Fin r) (Fin r) ℝ :=
  CBlock (fun i => (m.tbl i).θ) m.R

theorem isHermitian_coreG (m : UnalignedModelR μ M n d r rk) : m.coreG.IsHermitian :=
  isHermitian_CBlock _ _

/-! ### 2b. The stack bridge at general `r_i`

The paper's justification of `prop:stacksvd_subspace` is one sentence
(`main_paper.tex:842`): the stack is a rank-`r` spiked matrix whose population Gram is
`V C Vᵀ`. This section proves that sentence with no restriction on the `r_i`.
`signalFactorG` is the paper's signal factor `A` (block row `i` is `U_i Θ_i R_iᵀ`),
`signalPartG = A Vᵀ` is the mean of the stack, and `stackEG` is the stacked noise.

The one form the paper writes that is **not** stated here is `A = W Λ^{1/2} Qᵀ` with `W`
orthonormal, for the reason given at `r_i = 1`: at a zero eigenvalue of `C` the matching
column of `A Q` vanishes. -/

/-- Noise of the stack: the `M` noise matrices placed in one column of blocks, rows reindexed
as in `stackXG`. -/
noncomputable def stackEG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).E N ω p.2 k)

theorem stackEG_apply' (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackEG N ω q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).E N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- Unscaled noise of the stack, rows reindexed as in `stackXG`. The rank-one mirror is
`MultiTableModel.stackZ` (`StackSVD.lean:100`). With `stackE_eqG` this is the `Zu` datum that
`RankRStack` (`RankR/RMT/Stack.lean`) asks of a stacked table. -/
noncomputable def stackZG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).Z N ω p.2 k)

theorem stackZG_apply' (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackZG N ω q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).Z N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- `E_stack = d^{-1/2} Z_stack`: the stack scales as one table does, because every table
carries the same `d_N`. -/
theorem stackE_eqG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.stackEG N ω = (Real.sqrt (d N))⁻¹ • m.stackZG N ω := by
  ext q k
  rw [Matrix.smul_apply, stackEG_apply', stackZG_apply', SpikedModelR.E, Matrix.smul_apply,
    smul_eq_mul]

/-- The signal factor `A ∈ ℝ^{n_stack × r}` of the stack: block row `i` is `U_i Θ_i R_iᵀ`
(`main_paper.tex:818`, the commented display `[U_i Θ_i R_iᵀ]`). Entry `((i, a), k)` is
`∑_j (U_i)_{aj} θ_ij (R_i)_{kj}`. At `r_i = 1` the sum has one term and this is
`signalFactor` (`RankR/Subspace.lean:292`). -/
noncomputable def signalFactorG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin r) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin r))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k =>
      ∑ j : Fin (rk p.1), (m.tbl p.1).U N p.2 j * (m.tbl p.1).θ j * m.R p.1 k j)

theorem signalFactorG_apply' (m : UnalignedModelR μ M n d r rk) (N : ℕ)
    (q : Fin (∑ i, n i N)) (k : Fin r) :
    m.signalFactorG N q k
      = ∑ j : Fin (rk (finSigmaFinEquiv.symm q).1),
          (m.tbl (finSigmaFinEquiv.symm q).1).U N (finSigmaFinEquiv.symm q).2 j *
            (m.tbl (finSigmaFinEquiv.symm q).1).θ j *
            m.R (finSigmaFinEquiv.symm q).1 k j := rfl

/-- The mean of the stack, `E[X_stack] = A Vᵀ`: block row `i` is `U_i Θ_i (V R_i)ᵀ`. -/
noncomputable def signalPartG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.signalFactorG N * (m.V N)ᵀ

/-- Row `(i, a)` of `A Vᵀ`, with the table index a real variable. `m.hv i N` (that is
`V_i = V R_i`) is the only model field this step reads. -/
private theorem row_mul_V (m : UnalignedModelR μ M n d r rk) (N : ℕ) (i : Fin M)
    (a : Fin (n i N)) (k : Fin (d N)) :
    (∑ l : Fin r, (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * m.R i l j) * m.V N k l)
      = ∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * (m.tbl i).V N k j := by
  have hVi : ∀ j : Fin (rk i), (m.tbl i).V N k j = ∑ l : Fin r, m.V N k l * m.R i l j := by
    intro j
    rw [m.hv i N]
    simp [Matrix.mul_apply]
  have h1 : ∀ l : Fin r,
      (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * m.R i l j) * m.V N k l
        = ∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * (m.R i l j * m.V N k l) := by
    intro l
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [Finset.sum_congr rfl fun l _ => h1 l, Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [hVi j, Finset.mul_sum]
  exact Finset.sum_congr rfl fun l _ => by ring

/-- The entries of the mean of the stack: `∑_j (U_i)_{aj} θ_ij (V R_i)_{kj}`, by
`V_i = V R_i`. -/
theorem signalPartG_apply (m : UnalignedModelR μ M n d r rk) (N : ℕ)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.signalPartG N q k
      = ∑ j : Fin (rk (finSigmaFinEquiv.symm q).1),
          (m.tbl (finSigmaFinEquiv.symm q).1).U N (finSigmaFinEquiv.symm q).2 j *
            (m.tbl (finSigmaFinEquiv.symm q).1).θ j *
            (m.tbl (finSigmaFinEquiv.symm q).1).V N k j := by
  rw [signalPartG, Matrix.mul_apply]
  simp only [signalFactorG_apply', Matrix.transpose_apply]
  exact row_mul_V m N _ _ k

/-- The stack is a rank-`r` spiked matrix: `X_stack = A Vᵀ + E_stack`
(`main_paper.tex:818`). -/
theorem stackX_eqG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.stackXG N ω = m.signalPartG N + m.stackEG N ω := by
  ext q k
  rw [Matrix.add_apply, stackXG_apply', stackEG_apply', signalPartG_apply, spikedR_X_apply]

/-- `Uᵀ U = 1` inside one table, in the shape the two Gram computations of this section
consume. At `r_i = 1` this is `‖u_i‖ = 1`; at general `r_i` it is an `r_i × r_i` identity, so
the cross terms `j ≠ j'` drop as well. -/
private theorem table_UtU (m : UnalignedModelR μ M n d r rk) (i : Fin M) (N : ℕ)
    (x y : Fin (rk i) → ℝ) :
    (∑ a : Fin (n i N),
        (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * x j) *
          (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * y j))
      = ∑ j, (m.tbl i).θ j ^ 2 * (x j * y j) := by
  have hU : ∀ j j' : Fin (rk i),
      (∑ a : Fin (n i N), (m.tbl i).U N a j * (m.tbl i).U N a j')
        = if j = j' then (1 : ℝ) else 0 := by
    intro j j'
    have h := congrFun (congrFun ((m.tbl i).hU N) j) j'
    simpa [Matrix.mul_apply, Matrix.transpose_apply, Matrix.one_apply] using h
  have expand : ∀ a : Fin (n i N),
      (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * x j) *
          (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * y j)
        = ∑ j, ∑ j', ((m.tbl i).θ j * x j) * ((m.tbl i).θ j' * y j') *
            ((m.tbl i).U N a j * (m.tbl i).U N a j') := by
    intro a
    rw [Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun j' _ => by ring
  rw [Finset.sum_congr rfl fun a _ => expand a, Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_comm]
  have hinner : ∀ j' : Fin (rk i),
      (∑ a : Fin (n i N), ((m.tbl i).θ j * x j) * ((m.tbl i).θ j' * y j') *
          ((m.tbl i).U N a j * (m.tbl i).U N a j'))
        = ((m.tbl i).θ j * x j) * ((m.tbl i).θ j' * y j') *
            (if j = j' then (1 : ℝ) else 0) := by
    intro j'
    rw [← Finset.mul_sum, hU j j']
  rw [Finset.sum_congr rfl fun j' _ => hinner j']
  simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, if_true]
  ring

/-- `AᵀA = C`: the Gram of the signal factor is the core matrix. This is where
`U_iᵀ U_i = 1` enters. -/
theorem signalFactorG_transpose_mul_self (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    (m.signalFactorG N)ᵀ * m.signalFactorG N = m.coreG := by
  ext k l
  rw [Matrix.mul_apply, coreG, cblock_apply]
  simp only [Matrix.transpose_apply, signalFactorG_apply']
  have hstep := sum_stack_index (ν := fun i => n i N)
    (fun i a => (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * m.R i k j) *
      (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * m.R i l j))
  rw [hstep]
  exact Finset.sum_congr rfl fun i _ =>
    table_UtU m i N (fun j => m.R i k j) (fun j => m.R i l j)

/-- The population Gram of the stack, `E[X_stackᵀ] E[X_stack] = ∑_i ∑_j θ_ij² v_ij v_ijᵀ`
(`main_paper.tex:824`, the matrix `S`), with `v_ij` the `j`-th column of `V R_i`. -/
noncomputable def signalGramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  ∑ i, ∑ j : Fin (rk i), (m.tbl i).θ j ^ 2 •
    Matrix.vecMulVec (WithLp.ofLp ((m.tbl i).col N j)) (WithLp.ofLp ((m.tbl i).col N j))

/-- `signalGramG` is the Gram of the mean of the stack, so it is the paper's
`E[X_stackᵀ] E[X_stack]` and not a free-standing definition. -/
theorem signalGram_eq_signalPartG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    (m.signalPartG N)ᵀ * m.signalPartG N = m.signalGramG N := by
  ext k l
  rw [Matrix.mul_apply, signalGramG, Matrix.sum_apply]
  simp only [Matrix.transpose_apply, signalPartG_apply]
  have hstep := sum_stack_index (ν := fun i => n i N)
    (fun i a => (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * (m.tbl i).V N k j) *
      (∑ j, (m.tbl i).U N a j * (m.tbl i).θ j * (m.tbl i).V N l j))
  rw [hstep]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [table_UtU m i N (fun j => (m.tbl i).V N k j) (fun j => (m.tbl i).V N l j),
    Matrix.sum_apply]
  exact Finset.sum_congr rfl fun j _ => by
    simp [Matrix.vecMulVec_apply, SpikedModelR.col]

/-- `main_paper.tex:826`: the population Gram of the stack is `V C Vᵀ`. This is the sentence
the paper's proof rests on (`main_paper.tex:842`). -/
theorem signalGram_eqG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    m.signalGramG N = m.V N * m.coreG * (m.V N)ᵀ := by
  rw [← signalGram_eq_signalPartG, signalPartG, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.mul_assoc, ← Matrix.mul_assoc ((m.signalFactorG N)ᵀ),
    signalFactorG_transpose_mul_self, ← Matrix.mul_assoc]

/-- `λ_j(C)`, the `j`-th eigenvalue of the core matrix. The index is Mathlib's unsorted one;
every statement below either sums over `j` or reads the matching eigenvector, so no ordering
is needed. -/
noncomputable def coreEigG (m : UnalignedModelR μ M n d r rk) (j : Fin r) : ℝ :=
  m.isHermitian_coreG.eigenvalues j

/-- `V q_j`, the `j`-th population right singular vector of the stack, with `q_j` the matching
eigenvector of `C`. The paper's qualifier "if the `λ_j` are distinct"
(`main_paper.tex:831`) is not needed below: `signalGram_eq_sum_spikeVecG` holds for every
spectrum. -/
noncomputable def spikeVecG (m : UnalignedModelR μ M n d r rk) (j : Fin r) (N : ℕ) :
    EuclideanSpace ℝ (Fin (d N)) :=
  WithLp.toLp 2 ((m.V N).mulVec (WithLp.ofLp (m.isHermitian_coreG.eigenvectorBasis j)))

/-- `V q_j = ∑_k (q_j)_k (V e_k)`: the spike direction is the same linear combination of the
columns of `V` that `q_j` is of the standard basis. -/
theorem spikeVecG_eq_sum (m : UnalignedModelR μ M n d r rk) (j : Fin r) (N : ℕ) :
    m.spikeVecG j N
      = ∑ k : Fin r, m.isHermitian_coreG.eigenvectorBasis j k • m.colVecG N k := by
  refine PiLp.ext fun l => ?_
  rw [euclid_sum_apply]
  simp only [spikeVecG, colVecG, PiLp.toLp_apply, PiLp.smul_apply, smul_eq_mul, Matrix.mulVec,
    dotProduct]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The last equality of `main_paper.tex:826`: `V C Vᵀ = (V Q) Λ (V Q)ᵀ`. The stack has `r`
spikes, the `j`-th of strength `√(λ_j(C))` in the direction `V q_j`. -/
theorem signalGram_eq_sum_spikeVecG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    m.signalGramG N = ∑ j : Fin r, m.coreEigG j •
      Matrix.vecMulVec (WithLp.ofLp (m.spikeVecG j N)) (WithLp.ofLp (m.spikeVecG j N)) := by
  rw [signalGram_eqG]
  conv_lhs => rw [hermitian_eq_sum_rankOne m.isHermitian_coreG]
  rw [Matrix.mul_sum, Matrix.sum_mul]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Matrix.mul_smul, Matrix.smul_mul, mul_vecMulVec_mul_transpose]
  rfl

/-- `limitStackRG` of the model reads the eigenvalues of its core matrix. -/
theorem limitStackRG_eq (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ) :
    limitStackRG (fun i => (m.tbl i).θ) m.R c
      = ∑ j : Fin r, betaSq (Real.sqrt (m.coreEigG j)) (∑ i, c i) := rfl

/-! ### 3. The performance and the hypothesis structure -/

/-- Performance of stacksvd at rank `r` and general `r_i`:
`tr(Vᵀ P_top(X_stackᵀ X_stack, r) V)`, written as `∑_{k<r} ‖P_top (V e_k)‖²`. It equals the
paper's `‖V̂_stacksvdᵀ V‖_F²` whenever no eigenvalue of the Gram ties across the index `r`
boundary. See the header. -/
noncomputable def perfStackRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) : ℝ :=
  ∑ k : Fin r,
    ‖specProjTop (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r (m.colVecG N k)‖ ^ 2

/-- The rank-`r` spiked law for the unit-weight stack at general `r_i`. It is the black box of
`prop:stacksvd_subspace`, in the style of `UnalignedModel.SubspaceLaw`
(`RankR/SubspaceMain.lean:41`). `c` is the aspect ratio of the stack, that is `‖c‖₁ = ∑_i c_i`.

`align`: the top-`r` eigenspace of `X_stackᵀ X_stack` overlaps the spike direction `V q_j` by
`betaSq (√(λ_j(C))) ‖c‖₁`. This is `prop:single_table` read once per spike of the core matrix,
and section 2b proves that the stack is the rank-`r` spiked matrix the law describes.

The structure has one field. It carries no eigengap of `C`, no distinct-spike condition and no
gap field; see the header and choice 1 of `notes/archive/rankr_plan_B.md` section 4. -/
structure SubspaceLawG (m : UnalignedModelR μ M n d r rk) (c : ℝ) : Prop where
  align : ∀ j : Fin r, TendstoInProb μ
    (fun N ω => ‖specProjTop (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r
      (m.spikeVecG j N)‖ ^ 2) (betaSq (Real.sqrt (m.coreEigG j)) c)

/-- The performance is `tr(Vᵀ P V)`, so it does not depend on which orthonormal basis of the
column span of `V` is used: the columns `V e_k` may be replaced by the spike directions
`V q_j`. This is the only step of the proposition that touches the eigenvectors of `C`, and it
needs no eigengap, because `Q` is orthogonal. The rank-one mirror is
`perfStackR_eq_sum_spikeVec` (`RankR/Subspace.lean:454`); the step runs on `Fin r` and does
not see the `r_i`. -/
theorem perfStackRG_eq_sum_spikeVec (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.perfStackRG N ω = ∑ j : Fin r,
      ‖specProjTop (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r
        (m.spikeVecG j N)‖ ^ 2 := by
  simp only [perfStackRG]
  set P := specProjTop (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r with hP
  set Q : Fin r → Fin r → ℝ := fun j k => m.isHermitian_coreG.eigenvectorBasis j k with hQ
  have hPspike : ∀ j : Fin r, P (m.spikeVecG j N) = ∑ k : Fin r, Q j k • P (m.colVecG N k) := by
    intro j
    rw [m.spikeVecG_eq_sum j N, map_sum]
    exact Finset.sum_congr rfl fun k _ => by rw [map_smul]
  have hterm : ∀ j : Fin r, ‖P (m.spikeVecG j N)‖ ^ 2
      = ∑ k : Fin r, ∑ l : Fin r,
          Q j k * (Q j l * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ) := by
    intro j
    rw [← real_inner_self_eq_norm_sq, hPspike j, sum_inner]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [real_inner_smul_left, inner_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun l _ => by rw [real_inner_smul_right]
  symm
  calc ∑ j : Fin r, ‖P (m.spikeVecG j N)‖ ^ 2
      = ∑ j : Fin r, ∑ k : Fin r, ∑ l : Fin r,
          Q j k * (Q j l * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ) :=
        Finset.sum_congr rfl fun j _ => hterm j
    _ = ∑ k : Fin r, ∑ l : Fin r, ∑ j : Fin r,
          Q j k * (Q j l * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun k _ => Finset.sum_comm
    _ = ∑ k : Fin r, ∑ l : Fin r,
          (if k = l then (1 : ℝ) else 0) * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ := by
        refine Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun l _ => ?_
        have hassoc : ∀ j : Fin r,
            Q j k * (Q j l * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ)
              = (Q j k * Q j l) * ⟪P (m.colVecG N k), P (m.colVecG N l)⟫_ℝ :=
          fun j => (mul_assoc _ _ _).symm
        rw [Finset.sum_congr rfl fun j _ => hassoc j, ← Finset.sum_mul,
          sum_eigenvectorBasis_mul]
    _ = ∑ k : Fin r, ‖P (m.colVecG N k)‖ ^ 2 := by
        refine Finset.sum_congr rfl fun k _ => ?_
        simp only [ite_mul, one_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
        exact real_inner_self_eq_norm_sq _

end UnalignedModelR

end StackedSVD
