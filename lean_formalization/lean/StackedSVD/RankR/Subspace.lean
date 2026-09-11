/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Defs
import StackedSVD.StackSVD

/-!
# `prop:stacksvd_subspace`: unweighted stacksvd in the unaligned rank-`r` model

STATUS 2026-08-30: proved, 0 `sorry`. The review note is
`notes/archive/prop_stacksvd_subspace.md`; the audit that drove the interface is
`notes/archive/audit_stacksvd_subspace_2026-08-30.md` (decision D18 of `notes/FLAGGED.md`).

`assum:unaligned` (`main_paper.tex:753`) at `r_i = 1` gives each table one spike
`v_i = V R_i` inside a shared `r`-dimensional subspace `V`. The unit-weight stack is then a
rank-`r` spiked matrix. The paper writes it as (`main_paper.tex:818`)

```
X_stack = [U_1 Θ_1 R_1ᵀ Vᵀ + E_1; …; U_M Θ_M R_Mᵀ Vᵀ + E_M] = A Vᵀ + E_stack
```

with `A` the signal factor (block row `i` is `θ_i u_i R_iᵀ`), and its population Gram is

```
E[X_stackᵀ] E[X_stack] = ∑_i θ_i² v_i v_iᵀ = V C Vᵀ = (V Q) Λ (V Q)ᵀ,   C = ∑_i θ_i² R_i R_iᵀ.
```

Section 2b proves that chain: `stackX_eq`, `signalFactor_transpose_mul_self` (`AᵀA = C`),
`signalGram_eq_signalPart`, `signalGram_eq` and `signalGram_eq_sum_spikeVec`. So the stack has
`r` spikes `√(λ_j(C))` with directions `V q_j`, and aspect ratio `‖c‖₁ = ∑_i c_i`.
`prop:stacksvd_subspace` (`main_paper.tex:834`) reads the rank-one law of `prop:single_table`
once per spike of `C`:

```
‖V̂_stacksvdᵀ V‖_F² →p ∑_{j<r} (λ_j²(C) - ‖c‖₁) / (λ_j²(C) + λ_j(C)) 1{λ_j²(C) > ‖c‖₁}
                    = ∑_{j<r} betaSq (√(λ_j(C))) ‖c‖₁.
```

The second form is `hPerf_eq_betaSq` of `RankR/Example.lean` and is `limitStackR` here.

## Content

1. `Cmat θ R = (B_R)ᵀ B_R`, the core matrix `C = ∑_i R_i Θ_i² R_iᵀ` at `r_i = 1`
   (`cmat_eq_sum` is the paper's form), and `limitStackR`, the limit above.
2. `UnalignedModel.stackX`, the unit-weight stack `[X_1; …; X_M]`, its Gram `stackGram`
   (`= ∑_i X_iᵀ X_i`), and the core matrix `core` of the model.
2b. The stack bridge: `signalFactor` (`A`), `stackE`, `signalPart` (`= A Vᵀ`, the mean of the
   stack), `stackX_eq`, `signalFactor_transpose_mul_self`, `signalGram_eq_signalPart`,
   `signalGram_eq`, `signalGram_eq_sum_spikeVec`.
3. `perfStackR`, the top-`r` subspace performance, and `SubspaceLaw`, the one new hypothesis
   structure.
4. `prop_stacksvd_subspace`, the Layer 1 implication; `limitStackR_one_eq_stackSVDLimit`, the
   paper's own reduction to `prop:stacksvd_general` (`main_paper.tex:843`);
   `limitStackR_example` and `example_perfStackR_tendsto`, the stacksvd half of the worked
   example `eq:psi_equation`.

## The performance

The paper measures `‖V̂_stacksvdᵀ V‖_F²` with `V̂_stacksvd` the top `r` right singular vectors
of `X_stack`. Those vectors are orthonormal, so their projector is the top-`r` spectral
projector of `X_stackᵀ X_stack` and

```
‖V̂_stacksvdᵀ V‖_F² = tr(Vᵀ P_top(X_stackᵀ X_stack, r) V) = ∑_{k<r} ‖P_top (V e_k)‖².
```

`perfStackR` is that last sum. No inverse appears, unlike `perfR` of `RankR/Defs.lean`, whose
rows `v̂_i` are not orthonormal. The identity holds whenever no eigenvalue of the Gram ties
across the index `r` boundary, that is under `TopGap (m.stackGram N ω) _ r`. `V̂_stacksvd`
itself is not defined in Lean, so the identity is a doc-comment and not a theorem, exactly as
for `perfR` (`notes/archive/audit_rank_r_t7_2026-08-30.md`, item 5.1). It was measured to 1.2e-14
over 160 draws (`notes/archive/prop_stacksvd_subspace.md`, necessity scan row B1).

## No distinct-spike hypothesis, and no gap field

`SubspaceLaw` is stated on the top-`r` spectral **projector**, not on one eigenvector per
index, so it carries no eigengap and no distinct-spike condition. That matches the paper,
which says at `main_paper.tex:842` that the result holds when `C` has repeated eigenvalues.
The paper's own worked example at `ψ = 0` is such a case (`spec(C) = (θ², θ²)`), and the
Monte Carlo agrees there (bias -0.0017 at `d = 300`, 24 draws). See choice 6 of the note.

The structure has one field. An earlier draft carried a second field `topGap`, the rank-`r`
twin of `SingleTableLaw.topSimple`. The audit of 2026-08-30 (item 4.1) shows that field is
**false**, not almost surely false, at any `N` with `∑_i n_i N < r ≤ d N`, which
`UnalignedModel` permits: the Gram then has rank below `r`, so `λ_{r-1} = λ_r = 0`. The law
would be unsatisfiable and the proposition vacuous on such a model. Nothing consumed the
field, so D18 drops it.

`prop_stacksvd_subspace` and the `RankR.Example` namespace moved to
`RankR/SubspaceMain.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 0. Helper lemmas

Four small facts that the file uses more than once and that no other file states. -/

/-- A sum over the stacked row index `Fin (∑ i, ν i)` is the double sum over the block index.
The two rewrites are the ones `MultiTableModel.norm_stackU` uses in `StackSVD.lean`. -/
theorem sum_stack_index {α : Type*} [AddCommMonoid α] {M : ℕ} {ν : Fin M → ℕ}
    (F : (i : Fin M) → Fin (ν i) → α) :
    ∑ q : Fin (∑ i, ν i), F (finSigmaFinEquiv.symm q).1 (finSigmaFinEquiv.symm q).2
      = ∑ i, ∑ j, F i j := by
  rw [Equiv.sum_comp finSigmaFinEquiv.symm (fun p : (i : Fin M) × Fin (ν i) => F p.1 p.2)]
  exact Fintype.sum_sigma' F

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

/-! ### 1. The core matrix `C` and the limit -/

/-- The core matrix `C = ∑_i R_i Θ_i² R_iᵀ` (`main_paper.tex:826`) at `r_i = 1`, written as
`B_Rᵀ B_R` with `B_R = diag(θ) R_stack` of `RankR/Defs.lean`. The paper's sum form is
`cmat_eq_sum`. The argument `θ` is the vector of signal strengths, not the vector `β`. -/
noncomputable def Cmat {M r : ℕ} (θ : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Matrix (Fin r) (Fin r) ℝ :=
  (BR θ R)ᵀ * BR θ R

theorem isHermitian_Cmat {M r : ℕ} (θ : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (Cmat θ R).IsHermitian :=
  isHermitian_transpose_mul_self (BR θ R)

theorem cmat_apply {M r : ℕ} (θ : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (k l : Fin r) : Cmat θ R k l = ∑ i, θ i ^ 2 * (R i k * R i l) := by
  simp only [Cmat, Matrix.mul_apply, Matrix.transpose_apply, BR, Matrix.of_apply]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- The paper's form of the core matrix, `C = ∑_i R_i Θ_i² R_iᵀ` (`main_paper.tex:826`). -/
theorem cmat_eq_sum {M r : ℕ} (θ : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Cmat θ R = ∑ i, θ i ^ 2 •
      Matrix.vecMulVec (WithLp.ofLp (R i)) (WithLp.ofLp (R i)) := by
  ext k l
  rw [cmat_apply, Matrix.sum_apply]
  exact Finset.sum_congr rfl fun i _ => by
    simp [Matrix.vecMulVec_apply]

/-- The limit of `prop:stacksvd_subspace` (`main_paper.tex:836`): one rank-one performance
per spike of the core matrix, at signal `√(λ_j(C))` and aspect ratio `‖c‖₁`. The paper's
expression `(λ_j² - ‖c‖₁)/(λ_j² + λ_j) 1{λ_j² > ‖c‖₁}` is `betaSq (√λ_j) ‖c‖₁`, because
`(√λ)⁴ = λ²` and `(√λ)² = λ` for `λ ≥ 0`. The sum runs over the unsorted index of
`Matrix.IsHermitian.eigenvalues`, so it does not depend on any ordering of the spectrum. -/
noncomputable def limitStackR {M r : ℕ} (θ : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (c : Fin M → ℝ) : ℝ :=
  ∑ j : Fin r, betaSq (Real.sqrt ((isHermitian_Cmat θ R).eigenvalues j)) (∑ i, c i)

/-! ### 2. The unit-weight stack of an `UnalignedModel` -/

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The unit-weight stack `[X_1; …; X_M]`, with rows reindexed from the block index
`(i : Fin M) × Fin (n i N)` by `finSigmaFinEquiv`, as `MultiTableModel.stackZ` does. The
model of `MultiTableModel` cannot be reused: its field `hv` forces one shared spike, while
here table `i` has the spike `V R_i`. -/
noncomputable def stackX (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).X N ω p.2 k)

/-- Block `i` of the stack is `X_i`. -/
theorem stackX_apply (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) (i : Fin M)
    (j : Fin (n i N)) (k : Fin (d N)) :
    m.stackX N ω (finSigmaFinEquiv ⟨i, j⟩) k = (m.tbl i).X N ω j k := by
  simp [stackX, Matrix.reindex_apply, Matrix.submatrix_apply]

/-- The stack read at an arbitrary row index. This is the form the sums below consume. -/
theorem stackX_apply' (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackX N ω q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).X N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- Gram matrix of the stack. -/
noncomputable def stackGram (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.stackX N ω)ᵀ * m.stackX N ω

theorem isHermitian_stackGram (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    (m.stackGram N ω).IsHermitian :=
  isHermitian_transpose_mul_self (m.stackX N ω)

/-- `X_stackᵀ X_stack = ∑_i X_iᵀ X_i` (`main_paper.tex:823`). -/
theorem stackGram_eq_sum (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.stackGram N ω = ∑ i, ((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω := by
  ext k l
  simp only [stackGram, Matrix.mul_apply, Matrix.transpose_apply, Matrix.sum_apply,
    stackX_apply']
  exact sum_stack_index (ν := fun i => n i N)
    (fun i j => (m.tbl i).X N ω j k * (m.tbl i).X N ω j l)

/-- The core matrix of the model, `C = ∑_i θ_i² R_i R_iᵀ`. -/
noncomputable def core (m : UnalignedModel μ M n d r) : Matrix (Fin r) (Fin r) ℝ :=
  Cmat (fun i => (m.tbl i).θ) m.R

theorem isHermitian_core (m : UnalignedModel μ M n d r) : m.core.IsHermitian :=
  isHermitian_Cmat _ _

/-! ### 2b. The stack bridge

The paper's own justification of `prop:stacksvd_subspace` is one sentence
(`main_paper.tex:842`): the stack is a rank-`r` spiked matrix whose population Gram is
`V C Vᵀ`. This section proves that sentence. `signalFactor` is the paper's signal factor `A`
(block row `i` is `θ_i u_i R_iᵀ`), `signalPart = A Vᵀ` is the mean of the stack, and
`stackE` is the stacked noise.

The one form the paper writes that is **not** stated here is `A = W Λ^{1/2} Qᵀ` with `W`
orthonormal. At a zero eigenvalue of `C` the matching column of `A Q` vanishes, so `W` needs
an arbitrary orthonormal completion there; the two identities that the argument uses,
`AᵀA = C` and `V C Vᵀ = (V Q) Λ (V Q)ᵀ`, are stated instead
(`signalFactor_transpose_mul_self` and `signalGram_eq_sum_spikeVec`). -/

/-- Noise of the stack: the `M` noise matrices placed in one column of blocks, rows reindexed
as in `stackX`. -/
noncomputable def stackE (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).E N ω p.2 k)

theorem stackE_apply' (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackE N ω q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).E N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- The signal factor `A ∈ ℝ^{n_stack × r}` of the stack: block row `i` is `θ_i u_i R_iᵀ`
(`main_paper.tex:818`, the commented display `[U_i Θ_i R_iᵀ]`). -/
noncomputable def signalFactor (m : UnalignedModel μ M n d r) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin r) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin r))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) a =>
      (m.tbl p.1).θ * (m.tbl p.1).u N p.2 * m.R p.1 a)

theorem signalFactor_apply' (m : UnalignedModel μ M n d r) (N : ℕ)
    (q : Fin (∑ i, n i N)) (a : Fin r) :
    m.signalFactor N q a
      = (m.tbl (finSigmaFinEquiv.symm q).1).θ *
          (m.tbl (finSigmaFinEquiv.symm q).1).u N (finSigmaFinEquiv.symm q).2 *
          m.R (finSigmaFinEquiv.symm q).1 a := rfl

/-- The mean of the stack, `E[X_stack] = A Vᵀ`: block row `i` is `θ_i u_i v_iᵀ`. -/
noncomputable def signalPart (m : UnalignedModel μ M n d r) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.signalFactor N * (m.V N)ᵀ

/-- The entries of the mean of the stack: `θ_i u_i(j) v_i(k)`, by `v_i = V R_i`. -/
theorem signalPart_apply (m : UnalignedModel μ M n d r) (N : ℕ)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.signalPart N q k
      = (m.tbl (finSigmaFinEquiv.symm q).1).θ *
          (m.tbl (finSigmaFinEquiv.symm q).1).u N (finSigmaFinEquiv.symm q).2 *
          (m.tbl (finSigmaFinEquiv.symm q).1).v N k := by
  set i := (finSigmaFinEquiv.symm q).1 with hi
  set j := (finSigmaFinEquiv.symm q).2 with hj
  have hv : (m.tbl i).v N k = ∑ a, m.V N k a * m.R i a := by
    rw [m.hv i N]
    simp [Matrix.mulVec, dotProduct]
  rw [signalPart, Matrix.mul_apply, hv, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [signalFactor_apply', Matrix.transpose_apply]
  ring

/-- The stack is a rank-`r` spiked matrix: `X_stack = A Vᵀ + E_stack` (`main_paper.tex:818`). -/
theorem stackX_eq (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.stackX N ω = m.signalPart N + m.stackE N ω := by
  ext q k
  rw [Matrix.add_apply, stackX_apply', stackE_apply', signalPart_apply, SpikedModel.X]
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.vecMulVec_apply]
  ring

/-- `AᵀA = C`: the Gram of the signal factor is the core matrix. This is where `‖u_i‖ = 1`
enters. -/
theorem signalFactor_transpose_mul_self (m : UnalignedModel μ M n d r) (N : ℕ) :
    (m.signalFactor N)ᵀ * m.signalFactor N = m.core := by
  have hu : ∀ i : Fin M, ∑ j, ((m.tbl i).u N j) ^ 2 = 1 := by
    intro i
    have hnorm := EuclideanSpace.norm_sq_eq ((m.tbl i).u N)
    rw [(m.tbl i).hu N] at hnorm
    simpa [Real.norm_eq_abs, sq_abs] using hnorm.symm
  ext a b
  rw [Matrix.mul_apply, core, cmat_apply]
  have hstep := sum_stack_index (ν := fun i => n i N)
    (fun i j => ((m.tbl i).θ * (m.tbl i).u N j * m.R i a) *
      ((m.tbl i).θ * (m.tbl i).u N j * m.R i b))
  simp only [Matrix.transpose_apply, signalFactor_apply']
  rw [hstep]
  refine Finset.sum_congr rfl fun i _ => ?_
  have hterm : ∀ j : Fin (n i N),
      ((m.tbl i).θ * (m.tbl i).u N j * m.R i a) * ((m.tbl i).θ * (m.tbl i).u N j * m.R i b)
        = ((m.tbl i).θ ^ 2 * (m.R i a * m.R i b)) * ((m.tbl i).u N j) ^ 2 := by
    intro j; ring
  rw [Finset.sum_congr rfl fun j _ => hterm j, ← Finset.mul_sum, hu i, mul_one]

/-- The population Gram of the stack, `E[X_stackᵀ] E[X_stack] = ∑_i θ_i² v_i v_iᵀ`
(`main_paper.tex:824`, the matrix `S`). -/
noncomputable def signalGram (m : UnalignedModel μ M n d r) (N : ℕ) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  ∑ i, (m.tbl i).θ ^ 2 •
    Matrix.vecMulVec (WithLp.ofLp ((m.tbl i).v N)) (WithLp.ofLp ((m.tbl i).v N))

/-- `signalGram` is the Gram of the mean of the stack, so it is the paper's
`E[X_stackᵀ] E[X_stack]` and not a free-standing definition. -/
theorem signalGram_eq_signalPart (m : UnalignedModel μ M n d r) (N : ℕ) :
    (m.signalPart N)ᵀ * m.signalPart N = m.signalGram N := by
  have hu : ∀ i : Fin M, ∑ j, ((m.tbl i).u N j) ^ 2 = 1 := by
    intro i
    have hnorm := EuclideanSpace.norm_sq_eq ((m.tbl i).u N)
    rw [(m.tbl i).hu N] at hnorm
    simpa [Real.norm_eq_abs, sq_abs] using hnorm.symm
  ext k l
  rw [Matrix.mul_apply, signalGram, Matrix.sum_apply]
  have hstep := sum_stack_index (ν := fun i => n i N)
    (fun i j => ((m.tbl i).θ * (m.tbl i).u N j * (m.tbl i).v N k) *
      ((m.tbl i).θ * (m.tbl i).u N j * (m.tbl i).v N l))
  simp only [Matrix.transpose_apply, signalPart_apply]
  rw [hstep]
  refine Finset.sum_congr rfl fun i _ => ?_
  have hterm : ∀ j : Fin (n i N),
      ((m.tbl i).θ * (m.tbl i).u N j * (m.tbl i).v N k) *
          ((m.tbl i).θ * (m.tbl i).u N j * (m.tbl i).v N l)
        = ((m.tbl i).θ ^ 2 * ((m.tbl i).v N k * (m.tbl i).v N l)) * ((m.tbl i).u N j) ^ 2 := by
    intro j; ring
  rw [Finset.sum_congr rfl fun j _ => hterm j, ← Finset.mul_sum, hu i, mul_one]
  simp [Matrix.vecMulVec_apply]

/-- `main_paper.tex:826`: the population Gram of the stack is `V C Vᵀ`. This is the sentence
the paper's proof rests on (`main_paper.tex:842`). -/
theorem signalGram_eq (m : UnalignedModel μ M n d r) (N : ℕ) :
    m.signalGram N = m.V N * m.core * (m.V N)ᵀ := by
  rw [← signalGram_eq_signalPart, signalPart, Matrix.transpose_mul, Matrix.transpose_transpose,
    Matrix.mul_assoc, ← Matrix.mul_assoc ((m.signalFactor N)ᵀ),
    signalFactor_transpose_mul_self, ← Matrix.mul_assoc]

/-- `λ_j(C)`, the `j`-th eigenvalue of the core matrix. The index is Mathlib's unsorted one;
every statement below either sums over `j` or reads the matching eigenvector, so no ordering
is needed. -/
noncomputable def coreEig (m : UnalignedModel μ M n d r) (j : Fin r) : ℝ :=
  m.isHermitian_core.eigenvalues j

/-- `V q_j`, the `j`-th population right singular vector of the stack, with `q_j` the
matching eigenvector of `C`. The paper writes at `main_paper.tex:831`: "**If the `λ_j` are
distinct**, the population (true) right singular vectors of `X_stack` are `{V q_j}` with
signal strengths `{λ_j^{1/2}}`". The qualifier is the paper's, and no statement below needs
it: `signalGram_eq_sum_spikeVec` holds for every spectrum. -/
noncomputable def spikeVec (m : UnalignedModel μ M n d r) (j : Fin r) (N : ℕ) :
    EuclideanSpace ℝ (Fin (d N)) :=
  WithLp.toLp 2 ((m.V N).mulVec (WithLp.ofLp (m.isHermitian_core.eigenvectorBasis j)))

/-- `V q_j = ∑_k (q_j)_k (V e_k)`: the spike direction is the same linear combination of the
columns of `V` that `q_j` is of the standard basis. -/
theorem spikeVec_eq_sum (m : UnalignedModel μ M n d r) (j : Fin r) (N : ℕ) :
    m.spikeVec j N
      = ∑ k : Fin r, m.isHermitian_core.eigenvectorBasis j k • m.colVec N k := by
  refine PiLp.ext fun l => ?_
  rw [euclid_sum_apply]
  simp only [spikeVec, colVec, PiLp.toLp_apply, PiLp.smul_apply, smul_eq_mul, Matrix.mulVec,
    dotProduct]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The last equality of `main_paper.tex:826`: `V C Vᵀ = (V Q) Λ (V Q)ᵀ`. The stack has `r`
spikes, the `j`-th of strength `√(λ_j(C))` in the direction `V q_j`. -/
theorem signalGram_eq_sum_spikeVec (m : UnalignedModel μ M n d r) (N : ℕ) :
    m.signalGram N = ∑ j : Fin r, m.coreEig j •
      Matrix.vecMulVec (WithLp.ofLp (m.spikeVec j N)) (WithLp.ofLp (m.spikeVec j N)) := by
  rw [signalGram_eq]
  conv_lhs => rw [hermitian_eq_sum_rankOne m.isHermitian_core]
  rw [Matrix.mul_sum, Matrix.sum_mul]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Matrix.mul_smul, Matrix.smul_mul, mul_vecMulVec_mul_transpose]
  rfl

/-- `limitStackR` of the model reads the eigenvalues of its core matrix. -/
theorem limitStackR_eq (m : UnalignedModel μ M n d r) (c : Fin M → ℝ) :
    limitStackR (fun i => (m.tbl i).θ) m.R c
      = ∑ j : Fin r, betaSq (Real.sqrt (m.coreEig j)) (∑ i, c i) := rfl

/-! ### 3. The performance and the hypothesis structure -/

/-- Performance of stacksvd at rank `r`: `tr(Vᵀ P_top(X_stackᵀ X_stack, r) V)`, written as
`∑_{k<r} ‖P_top (V e_k)‖²`. It equals the paper's `‖V̂_stacksvdᵀ V‖_F²` whenever no
eigenvalue of the Gram ties across the index `r` boundary. See the header. -/
noncomputable def perfStackR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) : ℝ :=
  ∑ k : Fin r,
    ‖specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r (m.colVec N k)‖ ^ 2

/-- The performance is `tr(Vᵀ P V)`, so it does not depend on which orthonormal basis of the
column span of `V` is used: the columns `V e_k` may be replaced by the spike directions
`V q_j`. This is the only step of the proposition that touches the eigenvectors of `C`, and
it needs no eigengap, because `Q` is orthogonal. -/
theorem perfStackR_eq_sum_spikeVec (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    m.perfStackR N ω = ∑ j : Fin r,
      ‖specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r (m.spikeVec j N)‖ ^ 2 := by
  simp only [perfStackR]
  set P := specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r with hP
  set Q : Fin r → Fin r → ℝ := fun j k => m.isHermitian_core.eigenvectorBasis j k with hQ
  have hPspike : ∀ j : Fin r, P (m.spikeVec j N) = ∑ k : Fin r, Q j k • P (m.colVec N k) := by
    intro j
    rw [m.spikeVec_eq_sum j N, map_sum]
    exact Finset.sum_congr rfl fun k _ => by rw [map_smul]
  have hterm : ∀ j : Fin r, ‖P (m.spikeVec j N)‖ ^ 2
      = ∑ k : Fin r, ∑ l : Fin r,
          Q j k * (Q j l * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ) := by
    intro j
    rw [← real_inner_self_eq_norm_sq, hPspike j, sum_inner]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [real_inner_smul_left, inner_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun l _ => by rw [real_inner_smul_right]
  symm
  calc ∑ j : Fin r, ‖P (m.spikeVec j N)‖ ^ 2
      = ∑ j : Fin r, ∑ k : Fin r, ∑ l : Fin r,
          Q j k * (Q j l * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ) :=
        Finset.sum_congr rfl fun j _ => hterm j
    _ = ∑ k : Fin r, ∑ l : Fin r, ∑ j : Fin r,
          Q j k * (Q j l * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun k _ => Finset.sum_comm
    _ = ∑ k : Fin r, ∑ l : Fin r,
          (if k = l then (1 : ℝ) else 0) * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ := by
        refine Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun l _ => ?_
        have hassoc : ∀ j : Fin r,
            Q j k * (Q j l * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ)
              = (Q j k * Q j l) * ⟪P (m.colVec N k), P (m.colVec N l)⟫_ℝ :=
          fun j => (mul_assoc _ _ _).symm
        rw [Finset.sum_congr rfl fun j _ => hassoc j, ← Finset.sum_mul,
          sum_eigenvectorBasis_mul]
    _ = ∑ k : Fin r, ‖P (m.colVec N k)‖ ^ 2 := by
        refine Finset.sum_congr rfl fun k _ => ?_
        simp only [ite_mul, one_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
        exact real_inner_self_eq_norm_sq _

/-- The rank-`r` spiked law for the unit-weight stack, on the actual random matrices. It is
the black box of `prop:stacksvd_subspace`, in the style of `SpikedModel.SingleTableLaw`
(`RMT.lean`) and of `MultiTableModel.HeteroLaw` (`StackSVDWeighted.lean`). `c` is the aspect
ratio of the stack, that is `‖c‖₁ = ∑_i c_i`.

`align`: the top-`r` eigenspace of `X_stackᵀ X_stack` overlaps the spike direction `V q_j` by
`betaSq (√(λ_j(C))) ‖c‖₁`. This is `prop:single_table` read once per spike of the core matrix,
and section 2b proves that the stack is the rank-`r` spiked matrix the law describes.

The structure carries no eigengap of `C`, no distinct-spike condition and no gap field; see
the header and choice 6 of `notes/archive/prop_stacksvd_subspace.md`. -/
structure SubspaceLaw (m : UnalignedModel μ M n d r) (c : ℝ) : Prop where
  align : ∀ j : Fin r, TendstoInProb μ
    (fun N ω => ‖specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r
      (m.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (m.coreEig j)) c)

end UnalignedModel

/-- The paper's own reduction (`main_paper.tex:843`): at `r = r_i = 1` and `R_i = 1` the core
matrix is the scalar `‖θ‖₂²`, so `prop:stacksvd_subspace` is `prop:stacksvd_general`. -/
theorem limitStackR_one_eq_stackSVDLimit {M : ℕ} (θ c : Fin M → ℝ) :
    limitStackR θ (fun _ => oneVec) c = stackSVDLimit θ c := by
  have hT : (0 : ℝ) ≤ ∑ i, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  have htr : (Cmat θ (fun _ => oneVec)).trace = ∑ i, θ i ^ 2 := by
    rw [Matrix.trace, Fin.sum_univ_one, Matrix.diag_apply, cmat_apply]
    exact Finset.sum_congr rfl fun i _ => by simp [oneVec_apply]
  have heig : (isHermitian_Cmat θ (fun _ => oneVec)).eigenvalues 0 = ∑ i, θ i ^ 2 := by
    have h := (isHermitian_Cmat θ (fun _ => oneVec)).trace_eq_sum_eigenvalues
    rw [htr, Fin.sum_univ_one] at h
    exact_mod_cast h.symm
  rw [limitStackR, Fin.sum_univ_one, heig, betaSq_sqrt hT, stackSVDLimit]

end StackedSVD

