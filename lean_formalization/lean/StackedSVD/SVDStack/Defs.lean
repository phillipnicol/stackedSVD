/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT
import StackedSVD.Spectral
import StackedSVD.Prob.TendstoInProb
import StackedSVD.LinAlg.TopProjPerturb
import StatsMLlib.LinearAlgebra.Matrix.CourantFischer

/-!
# `thm:svd_stack_general`: definitions and deterministic identities

The first module of `StackedSVD.SVDStack`. It holds the objects that every later module of
the directory needs, and the deterministic identities that need no probability.

svdstack takes the top right singular vector `v̂_i` of each table and returns the top right
singular vector of the `M × d` matrix `Ṽ` whose rows are the `v̂_iᵀ`. Its limit is
`(βᵀ v_max(A_β))² / λ_max(A_β)` with `A_β = β βᵀ + diag(1 - β_i²)`.

## Content

1. `beta`, `Abeta`, `isHermitian_Abeta`, `vMax`.
2. The spectral helpers for `lamMax`, `vMax` and `topProj` that the whole directory uses.
3. `svdstackLimit`, `svdstackPerfClosed`, `lamMax_transpose_mul_self_eq`,
   `lamMax_mul_transpose_self_eq`, and step (v) of the paper's proof,
   `svdstackPerf_eq_closed`.
4. The per-table objects `MultiTableModel.tableGram`, `vhat`, `P`, `Vt`, the svdstack objects
   `svdstackGram`, `svdstackEst`, and the performance `svdstackPerf`.

The next modules are `StackedSVD.SVDStack.Deterministic` (the facts about `A_β` and the
general limit lemmas), `StackedSVD.SVDStack.Gram` (the consequences of the single-table
laws) and `StackedSVD.SVDStack.Main` (the theorems).

`isHermitian_mul_transpose_self` moved to `Defs.lean` (F28, 2026-09-08); item 1 above
is updated to match.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

/-! ### `β` and `A_β` -/

/-- `β_i = √(β_i²)` of `prop:single_table`. -/
noncomputable def beta (θ c : ℝ) : ℝ := Real.sqrt (betaSq θ c)

/-- `A_β = β βᵀ + diag(1 - β_1², ..., 1 - β_M²)` (`eq:A_beta_main_text`). -/
noncomputable def Abeta {M : ℕ} (β : Fin M → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.vecMulVec β β + Matrix.diagonal fun i => 1 - β i ^ 2

theorem isHermitian_Abeta {M : ℕ} (β : Fin M → ℝ) : (Abeta β).IsHermitian := by
  refine Matrix.IsHermitian.add ?_ (Matrix.isHermitian_diagonal _)
  unfold Matrix.IsHermitian
  ext i j
  simp [Matrix.vecMulVec_apply, mul_comm]

/-- `A_{β,w} = W A_β W` with `W = diag(w)` (`main_paper.tex:1284`). It is the entrywise limit
in probability of the weighted Gram matrix `(W Ṽ)(W Ṽ)ᵀ`. It lives here, and not in
`SVDStack/Weighted.lean`, so that `SVDStack/Deterministic.lean` can state
`sq_le_lamMax_AbetaW` (`notes/archive/audit_weighted_2026-08-30.md`, section 3). -/
noncomputable def AbetaW {M : ℕ} (w β : Fin M → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.diagonal w * Abeta β * Matrix.diagonal w

theorem isHermitian_AbetaW {M : ℕ} (w β : Fin M → ℝ) : (AbetaW w β).IsHermitian := by
  unfold Matrix.IsHermitian AbetaW
  rw [Matrix.conjTranspose_mul, Matrix.conjTranspose_mul, Matrix.diagonal_conjTranspose,
    (isHermitian_Abeta β : (Abeta β)ᴴ = Abeta β), Matrix.mul_assoc]
  simp

/-! ### Top eigenvector of a real symmetric matrix -/

/-- A unit top eigenvector of a real symmetric matrix: the vector of
`Matrix.IsHermitian.eigenvectorBasis` at the index of `lamMax` (`eigenvalues₀` is antitone,
index `0`). The sign is arbitrary, so every statement below squares it. Junk value `0` when
`d = 0`. -/
noncomputable def vMax {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    EuclideanSpace ℝ (Fin d) :=
  if h : 0 < d then
    hA.eigenvectorBasis (Fintype.equivOfCardEq (Fintype.card_fin _) ⟨0, by simpa using h⟩)
  else 0


/-! ### Spectral helpers

These are the deterministic facts about `lamMax`, `vMax` and `topProj` that the proofs of the
directory need. The public spectral interface lives in `Spectral.lean`; the helpers here are
public only when a later module of `StackedSVD/SVDStack/` uses them.
-/

section SpectralHelpers

variable {d : ℕ}

theorem toOp_of_mem_topSpace {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian}
    {x : EuclideanSpace ℝ (Fin d)} (hx : x ∈ topSpace A hA) :
    toOp A x = lamMax A hA • x :=
  Module.End.mem_eigenspace_iff.mp (by rwa [topSpace_eq_eigenspace] at hx)

/-- A nonzero vector forces a positive dimension. -/
private theorem pos_of_ne_zero {x : EuclideanSpace ℝ (Fin d)} (hx : x ≠ 0) : 0 < d := by
  rcases Nat.eq_zero_or_pos d with h | h
  · subst h
    exact absurd (by ext i; exact i.elim0) hx
  · exact h

theorem norm_vMax (hd : 0 < d) (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    ‖vMax A hA‖ = 1 := by
  rw [vMax, dif_pos hd]
  exact hA.eigenvectorBasis.norm_eq_one _

theorem mem_topSpace_vMax (hd : 0 < d) (A : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) : vMax A hA ∈ topSpace A hA := by
  rw [topSpace_eq_eigenspace]
  refine Module.End.mem_eigenspace_iff.mpr ?_
  set j : Fin d := Fintype.equivOfCardEq (Fintype.card_fin _) ⟨0, by simpa using hd⟩ with hj
  have hev : hA.eigenvalues j = lamMax A hA := by
    rw [lamMax, dif_pos hd, Matrix.IsHermitian.eigenvalues, hj, Equiv.symm_apply_apply]
  rw [vMax, dif_pos hd, ← hj]
  apply WithLp.ofLp_injective
  rw [← hev]
  simpa using hA.mulVec_eigenvectorBasis j

/-- Rayleigh upper bound: `⟪A x, x⟫ ≤ λ_max ‖x‖²`. -/
theorem inner_toOp_self_le (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (x : EuclideanSpace ℝ (Fin d)) : ⟪toOp A x, x⟫_ℝ ≤ lamMax A hA * ‖x‖ ^ 2 := by
  rcases eq_or_ne x 0 with rfl | hx
  · simp
  have hd : 0 < d := pos_of_ne_zero hx
  have hT : (toOp A).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hA
  have hn : Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) = Fintype.card (Fin d) :=
    finrank_euclideanSpace
  have hcard : 0 < Fintype.card (Fin d) := by simpa using hd
  set i0 : Fin (Fintype.card (Fin d)) := ⟨0, hcard⟩ with hi0
  have htop : hT.trailingEigenSubspace hn i0 = ⊤ := by
    apply Submodule.eq_top_of_finrank_eq
    rw [hT.finrank_trailingEigenSubspace hn i0, hn, hi0]
    simp
  have hray := hT.rayleighQuotient_le_eigenvalues_of_mem_trailingEigenSubspace hn i0
    (x := x) (by rw [htop]; trivial) hx
  have hlam : hT.eigenvalues hn i0 = lamMax A hA := by
    rw [lamMax, dif_pos hd]
    rfl
  rw [hlam, LinearMap.rayleighQuotient, div_le_iff₀ (by positivity : (0:ℝ) < ‖x‖ ^ 2)] at hray
  simpa using hray

/-- The Rayleigh quotient of `vMax` is `lamMax`. -/
private theorem inner_toOp_vMax (hd : 0 < d) (A : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) : ⟪toOp A (vMax A hA), vMax A hA⟫_ℝ = lamMax A hA := by
  rw [toOp_of_mem_topSpace (mem_topSpace_vMax hd A hA), real_inner_smul_left,
    real_inner_self_eq_norm_sq, norm_vMax hd A hA]
  ring

/-- Under simplicity, the top projector is the rank-one projector on a unit top eigenvector. -/
theorem topProj_eq_rankOne {A : Matrix (Fin d) (Fin d) ℝ} {hA : A.IsHermitian}
    (hsimple : TopSimple A hA) {e : EuclideanSpace ℝ (Fin d)}
    (he : e ∈ topSpace A hA) (hne : ‖e‖ = 1) :
    ∀ w, topProj A hA w = ⟪e, w⟫_ℝ • e := by
  have he0 : e ≠ 0 := by
    intro h
    rw [h, norm_zero] at hne
    exact absurd hne (by norm_num)
  have hle : (ℝ ∙ e) ≤ topSpace A hA := by
    rw [Submodule.span_singleton_le_iff_mem]
    exact he
  have hrank : Module.finrank ℝ (ℝ ∙ e) = Module.finrank ℝ (topSpace A hA) := by
    rw [finrank_span_singleton he0, hsimple]
  have hspan : (ℝ ∙ e) = topSpace A hA := Submodule.eq_of_le_of_finrank_eq hle hrank
  intro w
  change (topSpace A hA).starProjection w = _
  rw [← hspan, Submodule.starProjection_unit_singleton ℝ hne]

/-- The operator of a matrix product is the composition of the operators. -/
private theorem toOp_mul_apply {p q r : ℕ} (B : Matrix (Fin p) (Fin q) ℝ)
    (A : Matrix (Fin q) (Fin r) ℝ) (x : EuclideanSpace ℝ (Fin r)) :
    Matrix.toEuclideanLin (B * A) x = Matrix.toEuclideanLin B (Matrix.toEuclideanLin A x) := by
  apply WithLp.ofLp_injective
  change (B * A) *ᵥ WithLp.ofLp x = B *ᵥ (A *ᵥ WithLp.ofLp x)
  exact (Matrix.mulVec_mulVec _ _ _).symm

/-- The transpose is the adjoint of a real matrix operator. -/
private theorem inner_toEuclideanLin_transpose {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (x : EuclideanSpace ℝ (Fin q)) (y : EuclideanSpace ℝ (Fin p)) :
    ⟪Matrix.toEuclideanLin A x, y⟫_ℝ = ⟪x, Matrix.toEuclideanLin Aᵀ y⟫_ℝ := by
  rw [real_inner_eq_dotProduct, real_inner_eq_dotProduct]
  change (A *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp y = WithLp.ofLp x ⬝ᵥ (Aᵀ *ᵥ WithLp.ofLp y)
  rw [Matrix.mulVec_transpose, dotProduct_comm, Matrix.dotProduct_mulVec, dotProduct_comm]

end SpectralHelpers

/-! ### The svdstack limit and the paper's closed form -/

/-- Limit value of `thm:svd_stack_general`: `(βᵀ v_max(A_β))² / λ_max(A_β)`.

**Warning: `vMax` reads one arbitrary top eigenvector.** When the top eigenvalue of `A_β`
repeats, the value depends on which vector `Matrix.IsHermitian.eigenvectorBasis` returns. At
`β = (0.648, 0, 0)` the matrix `A_β` is the identity, the top eigenvalue has multiplicity 3,
and the two readings are `0` and `0.419` (mechanical audit 2026-08-31, finding 5). So this
definition is the paper's limit only under `TopSimple (Abeta β)`. Every theorem that uses it
carries that assumption: `thm_svd_stack_general` gets it from `hthr` (two tables with
`0 < β_i`) through `abeta_gap`, and the weighted twins get it from `hsimple`. The rank-`r`
`limitR` of `RankR/Defs.lean` states the same caution. -/
noncomputable def svdstackLimit {M : ℕ} (β : Fin M → ℝ) : ℝ :=
  (β ⬝ᵥ WithLp.ofLp (vMax (Abeta β) (isHermitian_Abeta β))) ^ 2
    / lamMax (Abeta β) (isHermitian_Abeta β)

/-- `S = ∑_i β_i² / (1 - β_i²)` of `thm:svdstack_weighted` (`main_paper.tex:510`). The single
copy: `Scalars.lean` and `SVDStack/Weighted.lean` both read it from here
(`notes/archive/audit_weighted_2026-08-30.md`, section 4). -/
noncomputable def Sval {M : ℕ} (β : Fin M → ℝ) : ℝ := ∑ i, β i ^ 2 / (1 - β i ^ 2)

/-- The optimal weights `w_i⋆ = 1/√(1 - β_i²)` (`main_paper.tex:521`, the rewriting of
`eq:svdstack.weight` above the detectability threshold). -/
noncomputable def optW {M : ℕ} (β : Fin M → ℝ) : Fin M → ℝ :=
  fun i => 1 / Real.sqrt (1 - β i ^ 2)

/-- The limit of optimally weighted svdstack, `S/(S+1)` (`thm:svdstack_weighted`). -/
noncomputable def svdstackLimitOpt {M : ℕ} (β : Fin M → ℝ) : ℝ := Sval β / (Sval β + 1)

/-- The paper's closed form of the svdstack performance, as a deterministic function of `Ṽ`
and `v`: `(xᵀ Ṽ v)² / (xᵀ Ṽ Ṽᵀ x)` with `x = v_max(Ṽ Ṽᵀ)`. The primary definition is the
projector form `MultiTableModel.svdstackPerf`; the two agree when the top eigenvalue of
`Ṽ Ṽᵀ` is simple (choice 7 of the note). -/
noncomputable def svdstackPerfClosed {M d : ℕ} (Vt : Matrix (Fin M) (Fin d) ℝ)
    (v : EuclideanSpace ℝ (Fin d)) : ℝ :=
  (WithLp.ofLp (vMax (Vt * Vtᵀ) (isHermitian_mul_transpose_self Vt)) ⬝ᵥ
      Vt.mulVec (WithLp.ofLp v)) ^ 2 /
    (WithLp.ofLp (vMax (Vt * Vtᵀ) (isHermitian_mul_transpose_self Vt)) ⬝ᵥ
      (Vt * Vtᵀ).mulVec (WithLp.ofLp (vMax (Vt * Vtᵀ) (isHermitian_mul_transpose_self Vt))))

/-- `λ_max(Aᵀ A) = λ_max(A Aᵀ)` when the left side is positive. Rayleigh bounds both ways:
`A e₁` and `Aᵀ x` carry the top eigenvectors across, with `‖A e₁‖² = λ_max(Aᵀ A)` and
`‖Aᵀ x‖² = λ_max(A Aᵀ)`. Extracted from the proof of `svdstackPerf_eq_closed`. -/
theorem lamMax_transpose_mul_self_eq {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (hpos : 0 < lamMax (Aᵀ * A) (isHermitian_transpose_mul_self A)) :
    lamMax (Aᵀ * A) (isHermitian_transpose_mul_self A)
      = lamMax (A * Aᵀ) (isHermitian_mul_transpose_self A) := by
  have hA := isHermitian_transpose_mul_self A
  have hB := isHermitian_mul_transpose_self A
  set L := Matrix.toEuclideanLin A with hLdef
  set Lt := Matrix.toEuclideanLin Aᵀ with hLtdef
  have hAcomp : ∀ z, toOp (Aᵀ * A) z = Lt (L z) := fun z => toOp_mul_apply Aᵀ A z
  have hBcomp : ∀ y, toOp (A * Aᵀ) y = L (Lt y) := fun y => toOp_mul_apply A Aᵀ y
  have hadj : ∀ z y, ⟪L z, y⟫_ℝ = ⟪z, Lt y⟫_ℝ := by
    intro z y
    simpa [hLdef, hLtdef] using inner_toEuclideanLin_transpose A z y
  have hnormL : ∀ z, ‖L z‖ ^ 2 = ⟪toOp (Aᵀ * A) z, z⟫_ℝ := by
    intro z
    rw [hAcomp z, real_inner_comm, ← hadj z (L z), real_inner_self_eq_norm_sq]
  have hnormLt : ∀ y, ‖Lt y‖ ^ 2 = ⟪toOp (A * Aᵀ) y, y⟫_ℝ := by
    intro y
    rw [hBcomp y, hadj (Lt y) y, real_inner_self_eq_norm_sq]
  -- both dimensions are positive
  have hd : 0 < q := by
    rcases Nat.eq_zero_or_pos q with h | h
    · exact absurd hpos (by rw [lamMax, dif_neg (by omega)]; exact lt_irrefl 0)
    · exact h
  have hM : 0 < p := by
    rcases Nat.eq_zero_or_pos p with h | h
    · exfalso
      subst h
      have hzero : ∀ z : EuclideanSpace ℝ (Fin q), toOp (Aᵀ * A) z = 0 := by
        intro z
        have : Aᵀ * A = 0 := by ext k l; simp [Matrix.mul_apply]
        rw [this]
        simp
      have : lamMax (Aᵀ * A) hA = 0 := by
        rw [← inner_toOp_vMax hd (Aᵀ * A) hA, hzero, inner_zero_left]
      rw [this] at hpos
      exact lt_irrefl 0 hpos
    · exact h
  set x := vMax (A * Aᵀ) hB with hxdef
  have hxn : ‖x‖ = 1 := norm_vMax hM _ hB
  have hBx : toOp (A * Aᵀ) x = lamMax (A * Aᵀ) hB • x :=
    toOp_of_mem_topSpace (mem_topSpace_vMax hM _ hB)
  set e1 := vMax (Aᵀ * A) hA with he1def
  have he1n : ‖e1‖ = 1 := norm_vMax hd _ hA
  have hAe1 : toOp (Aᵀ * A) e1 = lamMax (Aᵀ * A) hA • e1 :=
    toOp_of_mem_topSpace (mem_topSpace_vMax hd _ hA)
  have hl1 : ‖L e1‖ ^ 2 = lamMax (Aᵀ * A) hA := by
    rw [hnormL e1, hAe1, real_inner_smul_left, real_inner_self_eq_norm_sq, he1n]
    ring
  have hl2 : ‖Lt x‖ ^ 2 = lamMax (A * Aᵀ) hB := by
    rw [hnormLt x, hBx, real_inner_smul_left, real_inner_self_eq_norm_sq, hxn]
    ring
  have step1 : lamMax (Aᵀ * A) hA ≤ lamMax (A * Aᵀ) hB := by
    have h2 : Lt (L e1) = lamMax (Aᵀ * A) hA • e1 := by rw [← hAcomp]; exact hAe1
    have h3 : ⟪toOp (A * Aᵀ) (L e1), L e1⟫_ℝ = lamMax (Aᵀ * A) hA ^ 2 := by
      rw [← hnormLt (L e1), h2, norm_smul, he1n, mul_one, Real.norm_eq_abs, sq_abs]
    have h4 := inner_toOp_self_le (A * Aᵀ) hB (L e1)
    rw [h3, hl1] at h4
    nlinarith
  have step2 : lamMax (A * Aᵀ) hB ≤ lamMax (Aᵀ * A) hA := by
    have h2 : L (Lt x) = lamMax (A * Aᵀ) hB • x := by rw [← hBcomp]; exact hBx
    have h3 : ⟪toOp (Aᵀ * A) (Lt x), Lt x⟫_ℝ = lamMax (A * Aᵀ) hB ^ 2 := by
      rw [← hnormL (Lt x), h2, norm_smul, hxn, mul_one, Real.norm_eq_abs, sq_abs]
    have h4 := inner_toOp_self_le (Aᵀ * A) hA (Lt x)
    rw [h3, hl2] at h4
    nlinarith
  exact le_antisymm step1 step2

/-- The mirror of `lamMax_transpose_mul_self_eq`: positivity on the `A Aᵀ` side. -/
theorem lamMax_mul_transpose_self_eq {p q : ℕ} (A : Matrix (Fin p) (Fin q) ℝ)
    (hpos : 0 < lamMax (A * Aᵀ) (isHermitian_mul_transpose_self A)) :
    lamMax (Aᵀ * A) (isHermitian_transpose_mul_self A)
      = lamMax (A * Aᵀ) (isHermitian_mul_transpose_self A) := by
  have h1 : lamMax (Aᵀᵀ * Aᵀ) (isHermitian_transpose_mul_self Aᵀ)
      = lamMax (A * Aᵀ) (isHermitian_mul_transpose_self A) :=
    lamMax_congr (by rw [Matrix.transpose_transpose]) _ _
  have h2 : lamMax (Aᵀ * Aᵀᵀ) (isHermitian_mul_transpose_self Aᵀ)
      = lamMax (Aᵀ * A) (isHermitian_transpose_mul_self A) :=
    lamMax_congr (by rw [Matrix.transpose_transpose]) _ _
  have hpos' : 0 < lamMax (Aᵀᵀ * Aᵀ) (isHermitian_transpose_mul_self Aᵀ) := by
    rw [h1]
    exact hpos
  rw [← h1, ← h2]
  exact (lamMax_transpose_mul_self_eq Aᵀ hpos').symm

/-- Step (v) of the paper's proof, as a deterministic identity. With `x = v_max(Ṽ Ṽᵀ)` a unit
top eigenvector of the `M × M` Gram matrix, `Ṽᵀ x / ‖Ṽᵀ x‖` is a unit top eigenvector of the
`d × d` Gram matrix `Ṽᵀ Ṽ`, and `‖Ṽᵀ x‖² = xᵀ Ṽ Ṽᵀ x`. So the projector form of the svdstack
performance equals the paper's ratio. Simplicity gives `‖P_top w‖² = ⟪v̂, w⟫²`
(`overlap_eq_inner_sq`), and `hpos` excludes `Ṽ = 0`, where the left side is `‖v‖²` and the
right side is `0/0 = 0`. This is the linking lemma that choice 7 of
`notes/archive/thm_svd_stack_general.md` asks for. -/
theorem svdstackPerf_eq_closed {M d : ℕ} (Vt : Matrix (Fin M) (Fin d) ℝ)
    (v : EuclideanSpace ℝ (Fin d))
    (hsimple : TopSimple (Vtᵀ * Vt) (isHermitian_transpose_mul_self Vt))
    (hpos : 0 < lamMax (Vtᵀ * Vt) (isHermitian_transpose_mul_self Vt)) :
    ‖topProj (Vtᵀ * Vt) (isHermitian_transpose_mul_self Vt) v‖ ^ 2 =
      svdstackPerfClosed Vt v := by
  have hA := isHermitian_transpose_mul_self Vt
  have hB := isHermitian_mul_transpose_self Vt
  set L := Matrix.toEuclideanLin Vt with hLdef
  set Lt := Matrix.toEuclideanLin Vtᵀ with hLtdef
  have hAcomp : ∀ z, toOp (Vtᵀ * Vt) z = Lt (L z) := fun z => toOp_mul_apply Vtᵀ Vt z
  have hBcomp : ∀ y, toOp (Vt * Vtᵀ) y = L (Lt y) := fun y => toOp_mul_apply Vt Vtᵀ y
  have hadj : ∀ z y, ⟪L z, y⟫_ℝ = ⟪z, Lt y⟫_ℝ := by
    intro z y
    simpa [hLdef, hLtdef] using inner_toEuclideanLin_transpose Vt z y
  have hnormLt : ∀ y, ‖Lt y‖ ^ 2 = ⟪toOp (Vt * Vtᵀ) y, y⟫_ℝ := by
    intro y
    rw [hBcomp y, hadj (Lt y) y, real_inner_self_eq_norm_sq]
  -- both dimensions are positive
  have hd : 0 < d := by
    rcases Nat.eq_zero_or_pos d with h | h
    · exact absurd hpos (by rw [lamMax, dif_neg (by omega)]; exact lt_irrefl 0)
    · exact h
  have hM : 0 < M := by
    rcases Nat.eq_zero_or_pos M with h | h
    · exfalso
      subst h
      have hzero : ∀ z : EuclideanSpace ℝ (Fin d), toOp (Vtᵀ * Vt) z = 0 := by
        intro z
        have : Vtᵀ * Vt = 0 := by ext k l; simp [Matrix.mul_apply]
        rw [this]
        simp
      have : lamMax (Vtᵀ * Vt) hA = 0 := by
        rw [← inner_toOp_vMax hd (Vtᵀ * Vt) hA, hzero, inner_zero_left]
      rw [this] at hpos
      exact lt_irrefl 0 hpos
    · exact h
  -- the two top eigenvalues agree (`lamMax_transpose_mul_self_eq`)
  set x := vMax (Vt * Vtᵀ) hB with hxdef
  have hxn : ‖x‖ = 1 := norm_vMax hM _ hB
  have hBx : toOp (Vt * Vtᵀ) x = lamMax (Vt * Vtᵀ) hB • x :=
    toOp_of_mem_topSpace (mem_topSpace_vMax hM _ hB)
  have hl2 : ‖Lt x‖ ^ 2 = lamMax (Vt * Vtᵀ) hB := by
    rw [hnormLt x, hBx, real_inner_smul_left, real_inner_self_eq_norm_sq, hxn]
    ring
  have hll : lamMax (Vtᵀ * Vt) hA = lamMax (Vt * Vtᵀ) hB :=
    lamMax_transpose_mul_self_eq Vt hpos
  have hposB : 0 < lamMax (Vt * Vtᵀ) hB := by
    rw [← hll]
    exact hpos
  -- the unit top eigenvector of `Ṽᵀ Ṽ` built from `x`
  have hrpos : 0 < ‖Lt x‖ := by
    rcases (norm_nonneg (Lt x)).lt_or_eq with h | h
    · exact h
    · exfalso; rw [← h] at hl2; simp at hl2; linarith
  set e := ‖Lt x‖⁻¹ • Lt x with hedef
  have hen : ‖e‖ = 1 := by
    rw [hedef, norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
    field_simp
  have hALtx : toOp (Vtᵀ * Vt) (Lt x) = lamMax (Vt * Vtᵀ) hB • Lt x := by
    rw [hAcomp (Lt x), ← hBcomp x, hBx, map_smul]
  have hemem : e ∈ topSpace (Vtᵀ * Vt) hA := by
    rw [topSpace_eq_eigenspace]
    refine Module.End.mem_eigenspace_iff.mpr ?_
    rw [hedef, map_smul, hALtx, hll, smul_comm]
  -- both sides equal `⟪x, Ṽ v⟫² / λ`
  have hproj := topProj_eq_rankOne hsimple hemem hen v
  have hinner : ⟪e, v⟫_ℝ = ‖Lt x‖⁻¹ * ⟪x, L v⟫_ℝ := by
    rw [hedef, real_inner_smul_left]
    congr 1
    rw [real_inner_comm, ← hadj v x, real_inner_comm]
  have hnum : WithLp.ofLp x ⬝ᵥ Vt.mulVec (WithLp.ofLp v) = ⟪x, L v⟫_ℝ := by
    rw [real_inner_eq_dotProduct]
    rfl
  have hden : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x) = lamMax (Vt * Vtᵀ) hB := by
    have : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x)
        = ⟪x, toOp (Vt * Vtᵀ) x⟫_ℝ := by
      rw [real_inner_eq_dotProduct]
      rfl
    rw [this, hBx, real_inner_smul_right, real_inner_self_eq_norm_sq, hxn]
    ring
  rw [svdstackPerfClosed, ← hxdef, hnum, hden, hproj, norm_smul, mul_pow, hen, one_pow,
    mul_one, Real.norm_eq_abs, sq_abs, hinner]
  rw [mul_pow, ← hl2]
  field_simp

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section SVDStack

/-! ### The per-table estimators -/

/-- Gram matrix `X_iᵀ X_i` of table `i`. -/
noncomputable def tableGram (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ := ((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω

theorem isHermitian_tableGram (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    (m.tableGram i N ω).IsHermitian :=
  isHermitian_transpose_mul_self ((m.tbl i).X N ω)

/-- `v̂_i = v_max(X_i)`, the top right singular vector of table `i`, with the paper's sign
convention `⟪v̂_i, v⟫ ≥ 0` (choice 2 of the note). The sign test reads the unknown `v N`, so
`vhat` is not an observable estimator: the observable one differs from it by a sign, which
`svdstackPerf` squares away, so the performance is the same. The convention is what lets
`align` and `gram` state signed limits. -/
noncomputable def vhat (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω), (m.tbl i).v N⟫_ℝ then
    vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
  else -vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)

/-- `P_i`, the top eigenprojector of `X_iᵀ X_i` as a matrix. It equals `v̂_i v̂_iᵀ` on the
almost sure event `SingleTableLaw.topSimple`. -/
noncomputable def P (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  Matrix.toEuclideanLin.symm
    (topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)).toLinearMap

/-- `Ṽ`, the `M × d` matrix whose row `i` is `v̂_iᵀ`. -/
noncomputable def Vt (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin M) (Fin (d N)) ℝ :=
  Matrix.of fun i k => m.vhat i N ω k

/-! ### The svdstack estimator -/

/-- `∑_i P_i`. It equals `Ṽᵀ Ṽ` on the almost sure event that every table has a simple top
eigenvalue, and its top eigenvector is `v̂_svdstack`. -/
noncomputable def svdstackGram (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ := ∑ i, m.P i N ω

theorem isHermitian_P (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    (m.P i N ω).IsHermitian := by
  rw [← Matrix.isSymmetric_toEuclideanLin_iff, MultiTableModel.P, LinearEquiv.apply_symm_apply]
  exact Submodule.starProjection_isSymmetric _

theorem isHermitian_svdstackGram (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    (m.svdstackGram N ω).IsHermitian := by
  change (∑ i, m.P i N ω)ᴴ = ∑ i, m.P i N ω
  rw [Matrix.conjTranspose_sum]
  exact Finset.sum_congr rfl fun i _ => m.isHermitian_P i N ω

/-- `v̂_svdstack`, a unit top eigenvector of `∑_i P_i`. The sign is arbitrary; the
performance squares it. -/
noncomputable def svdstackEst (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    EuclideanSpace ℝ (Fin (d N)) :=
  vMax (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω)

/-! ### The svdstack performance -/

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-- Performance of svdstack in projector form: `‖P_top(∑_i P_i) v‖²`, the same form as
`overlap` in `Defs.lean` (choice 7 of the note). -/
noncomputable def svdstackPerf (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) : ℝ :=
  ‖topProj (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) ((m.tbl 0).v N)‖ ^ 2

end SVDStack

end MultiTableModel

end StackedSVD
