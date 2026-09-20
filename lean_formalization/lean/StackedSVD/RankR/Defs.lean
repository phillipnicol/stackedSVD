/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Defs
import StackedSVD.SVDStack.DelocDir
import StackedSVD.LinAlg.SpecProjPerturb

/-!
# Section 7 of the paper, the `r_i = 1` slice: model, objects and statements

STATUS 2026-08-30: proved, 0 `sorry` (tasks T3 and T4 of `notes/RANK_R_PLAN.md`). The review
note is `notes/archive/rank_r_defs.md`.

`assum:unaligned` (`main_paper.tex:753`) gives every table its own right singular subspace
`V R_i` inside one shared `r`-dimensional subspace `V`. This file takes the slice `r_i = 1`,
where each table has one spike `v_i = V R_i` with `R_i` a unit vector of `ℝ^r`, so `r̃ = M`
and the per-table random matrix theory input is the already proved rank-one
`SpikedModel.SingleTableLaw`. No new hypothesis structure appears.

## Content

1. `BR`, `Dmat`, `AbetaR`: the paper's `B_R = diag(β) R_stack`, `D = diag(1 - β_ij²)` and
   `A_{β,R} = B_R B_Rᵀ + D`, with the entrywise form `abetaR_apply`.
2. `UnalignedModel`: `M` tables on one probability space per `N`, a shared `V N` with
   orthonormal columns, and per-table unit `R i` with `v_i = V R_i`. `MultiTableModel` cannot
   be reused: its field `hv` forces one shared `v`.
3. The per-table estimators `tableGram`, `vhat`, `Vt`, the two limit objects `gram`
   (`Ṽ Ṽᵀ`) and `VtV` (`Ṽ V`), and the performance `perfR`.
4. `lem_general_rank_delocalization`, `gramR`, `VtV_tendsto`,
   `prop_general_rank_unweighted_svdstack` and its Gaussian corollary.
5. Positivity: `abetaR_posSemidef_sub_Dmat`, `dmat_posDef`, `abetaR_posDef`. They give
   `A_{β,R} ⪰ D ≻ 0`, so the inverse inside `specInvTop (A_{β,R}) r` never meets a zero
   eigenvalue. Task T4 needs this for the continuity hypothesis `0 < λ_{r-1}`.
6. The rank-one reduction `abetaR_reduce`, `BR_reduce`, `limitR_one_eq_svdstackLimit`: at
   `r = 1` with every `R_i = 1` the objects of this file are the rank-one objects of
   `SVDStack/Defs.lean`. This is the paper's own consistency check (`main_paper.tex:810`).

## The performance

The paper measures `‖V̂_svdstackᵀ V‖_F²` with `V̂_svdstack` the top `r` right singular vectors
of `Ṽ`. Its own first display (`main_paper.tex:1982`) writes
`V̂_svdstackᵀ V = Λ_r^{-1/2}(Ṽ Ṽᵀ) Q_r(Ṽ Ṽᵀ)ᵀ Ṽ V`, so

```
‖V̂_svdstackᵀ V‖_F² = tr( (Ṽ V)ᵀ (specInvTop (Ṽ Ṽᵀ) r) (Ṽ V) )
```

whenever `λ_r(Ṽ Ṽᵀ) > 0`. `perfR` is that right hand side. It drops the paper's `Q`, `Λ` and
the block-diagonal alignment matrices `O^{(d)}`, which exist only because the paper carries an
eigenvector matrix instead of a projector, and it avoids Mathlib's scoped Frobenius norm
instance. The two agree to 1.8e-15 over 480 Monte Carlo draws
(`notes/archive/rank_r_defs.md`, necessity scan row M2).

## Route of the proofs (T3, T4)

The rank-`r` twin of `SVDStack/Gram.lean` and `SVDStack/Main.lean`, with `v_i = V R_i` in
place of the shared `v`:

* `align_inner_det`: `⟪v̂_i, y_N⟫ → β_i ⟪v_i, y_N⟫` for a deterministic family `y` of norm at
  most one whose inner product with `v_i` does not depend on `N`. The proof splits
  `y_N = ⟪v_i, y_N⟫ v_i + perpOf v_i y_N` and kills the second half with
  `SingleTableLaw.delocUniform` at the normalized perpendicular part.
* `VtV_tendsto` is `align_inner_det` at `y_N = ` column `k` of `V N`, where `⟪v_i, y_N⟫`
  is `(R_i)_k` by `hV`.
* `lem_general_rank_delocalization` is `MultiTableModel.lem_delocalization` with table `j`
  tested against `v_i`, not against a shared `v`: `align_inner_det` gives
  `⟪v̂_j, v_i⟫ → β_j ⟪R_j, R_i⟫`, and `measure_deloc_le` (Fubini on the product law) kills the
  perpendicular part.
* `prop_general_rank_unweighted_svdstack` needs no good event: `perfR` is
  `traceFun r` of the joint entry family `(Ṽ Ṽᵀ, Ṽ V)` at **every** `ω`, and
  `continuousAt_trace_specInvTop` is continuity at the limit point `A_{β,R}` only. The gap and
  the positivity of `λ_{r-1}(A_{β,R})` (from `abetaR_posDef`) are hypotheses of that
  continuity, not events.

`align_inner` through `prop_general_rank_unweighted_svdstack_gaussian` moved to
`RankR/Unweighted.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### 1. `B_R`, `D` and `A_{β,R}` -/

/-- `B_R = diag(β) R_stack ∈ ℝ^{M × r}`, the matrix with row `i` equal to `β_i R_iᵀ`
(`main_paper.tex:1951`). At `r_i = 1` the row index `(i, j)` of the paper is the table index
`i`, because `r̃ = ∑_i r_i = M`. -/
noncomputable def BR {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Matrix (Fin M) (Fin r) ℝ :=
  Matrix.of fun i k => β i * R i k

/-- `D = diag(1 - β_ij²)` (`main_paper.tex:2050`). -/
noncomputable def Dmat {M : ℕ} (β : Fin M → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.diagonal fun i => 1 - β i ^ 2

/-- `A_{β,R} = B_R B_Rᵀ + D` (`main_paper.tex:2053`). The paper defines the matrix entrywise
(diagonal `1`, off-diagonal `β_i β_j ⟪R_i, R_j⟫`); `abetaR_apply` and `abetaR_diag` recover
that form. At `r = 1` with every `R_i = 1` this is `Abeta β`. -/
noncomputable def AbetaR {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    Matrix (Fin M) (Fin M) ℝ :=
  BR β R * (BR β R)ᵀ + Dmat β

theorem isHermitian_AbetaR {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (AbetaR β R).IsHermitian :=
  Matrix.IsHermitian.add (isHermitian_mul_transpose_self (BR β R))
    (Matrix.isHermitian_diagonal _)

/-- The paper's entrywise definition of `A_{β,R}`. -/
theorem abetaR_apply {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    (i j : Fin M) :
    AbetaR β R i j = β i * β j * ⟪R i, R j⟫_ℝ + (if i = j then 1 - β i ^ 2 else 0) := by
  simp only [AbetaR, Matrix.add_apply, Matrix.mul_apply, Matrix.transpose_apply, BR,
    Matrix.of_apply, Dmat, Matrix.diagonal_apply]
  congr 1
  rw [real_inner_eq_dotProduct, dotProduct, Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The diagonal of `A_{β,R}` is `1`, which is the diagonal of `Ṽ Ṽᵀ` at every `N`. This is
where `‖R_i‖ = 1` is used. -/
theorem abetaR_diag {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r))
    {i : Fin M} (hR : ‖R i‖ = 1) : AbetaR β R i i = 1 := by
  rw [abetaR_apply, if_pos rfl, real_inner_self_eq_norm_sq, hR]
  ring

/-! #### Positivity of `A_{β,R}`

`A_{β,R} = B_R B_Rᵀ + D` with `B_R B_Rᵀ` positive semidefinite, so `A_{β,R} ⪰ D`, and `D` is
positive definite as soon as every `β_i < 1`, which `beta_mem_Ico` gives for `c_i > 0`. The
numeric audit measures the slack: `λ_min(A_{β,R}) - min_i (1 - β_i²) ≥ 1.6e-6` over 50000
draws (`notes/archive/audit_rank_r_2026-08-30.md`, row 3b). -/

/-- `A_{β,R} - D = B_R B_Rᵀ` is positive semidefinite, so `A_{β,R} ⪰ D` in the Loewner
order. -/
theorem abetaR_posSemidef_sub_Dmat {M r : ℕ} (β : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) : (AbetaR β R - Dmat β).PosSemidef := by
  have h : AbetaR β R - Dmat β = BR β R * (BR β R)ᵀ := by
    rw [AbetaR, add_sub_cancel_right]
  rw [h]
  simpa using Matrix.posSemidef_self_mul_conjTranspose (BR β R)

/-- `D = diag(1 - β_i²)` is positive definite when every `β_i` lies in `[0, 1)`, which
`beta_mem_Ico` gives from `0 < c_i`. -/
theorem dmat_posDef {M : ℕ} {β : Fin M → ℝ} (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) :
    (Dmat β).PosDef := by
  refine Matrix.PosDef.diagonal fun i => ?_
  nlinarith [h0 i, h1 i]

/-- A positive semidefinite matrix plus a positive definite one is positive definite. Mathlib
v4.33.0 has `Matrix.PosSemidef.add` but no mixed form; this is the same proof. -/
theorem posSemidef_add_posDef {p : ℕ} {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.PosSemidef)
    (hB : B.PosDef) : (A + B).PosDef :=
  ⟨hA.isHermitian.add hB.isHermitian, fun x hx => by
    simpa [mul_add, add_mul] using add_pos_of_nonneg_of_pos (hA.2 x) (hB.2 hx)⟩

/-- `A_{β,R}` is positive definite when every `β_i` lies in `[0, 1)`: it is `D ≻ 0` plus the
positive semidefinite `B_R B_Rᵀ`. -/
theorem abetaR_posDef {M r : ℕ} {β : Fin M → ℝ} (R : Fin M → EuclideanSpace ℝ (Fin r))
    (h0 : ∀ i, 0 ≤ β i) (h1 : ∀ i, β i < 1) : (AbetaR β R).PosDef := by
  rw [AbetaR]
  exact posSemidef_add_posDef
    (by simpa using Matrix.posSemidef_self_mul_conjTranspose (BR β R)) (dmat_posDef h0 h1)

/-- The limit of `prop:general_rank_unweighted_svdstack`, in the same trace form as `perfR`:
`tr(B_Rᵀ (specInvTop A_{β,R} r) B_R)`, which is the paper's
`‖Λ^{-1/2} Qᵀ diag(β) R_stack‖_F²`. -/
noncomputable def limitR {M r : ℕ} (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    ℝ :=
  Matrix.trace ((BR β R)ᵀ * specInvTop (AbetaR β R) (isHermitian_AbetaR β R) r * BR β R)

/-! ### 1b. The trace functional and the rank-one form of `specInvTop`

`continuousAt_trace_specInvTop` (`LinAlg/SpecProjPerturb.lean`) is stated on the disjoint union
of the two entry families, which is the shape `TendstoInProbPi.comp_continuous` consumes.
`traceFun` names that functional, `traceFun_eq` evaluates it at a pair of matrices, and
`specInvTop_one_mulVec` is the rank-one collapse that the reduction lemma of section 5 needs.
-/

section TraceFun

variable {p q : ℕ}

/-- The performance functional `z ↦ tr(Bᵀ (specInvTop A r) B)` read off the joint entry family
`z`, in the exact shape of `continuousAt_trace_specInvTop`. Public since 2026-08-30: the proof
of `thm_gen_rank_weight_svdstak_general` in `RankR/Weighted.lean` uses the same functional with
`W Ṽ` in place of `Ṽ`. -/
noncomputable def traceFun (r : ℕ)
    (z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ) : ℝ :=
  Matrix.trace ((Matrix.of fun i k => z (Sum.inr (i, k)))ᵀ *
    specInvTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) r *
    (Matrix.of fun i k => z (Sum.inr (i, k))))

/-- `traceFun` at the entry family of a symmetric `C` and an arbitrary `Y`. Public since
2026-08-30, with `traceFun`. -/
theorem traceFun_eq {r : ℕ} {C : Matrix (Fin p) (Fin p) ℝ} (hC : C.IsHermitian)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    traceFun (q := q) r (Sum.elim (fun t : Fin p × Fin p => C t.1 t.2)
        (fun t : Fin p × Fin q => Y t.1 t.2))
      = Matrix.trace (Yᵀ * specInvTop C hC r * Y) := by
  simp only [traceFun, Sum.elim_inl, Sum.elim_inr]
  rw [specInvTop_congr_mat (isHermitian_symMat _) hC (symMat_entries hC) r]
  rfl

/-- `continuousAt_trace_specInvTop` in the `traceFun` name. Public since 2026-08-30, with
`traceFun`. -/
theorem continuousAt_traceFun {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p)) (hgap : TopGap A hA r)
    (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩) (B : Matrix (Fin p) (Fin q) ℝ) :
    ContinuousAt (traceFun (p := p) (q := q) r)
      (Sum.elim (fun t : Fin p × Fin p => A t.1 t.2) (fun t : Fin p × Fin q => B t.1 t.2)) :=
  continuousAt_trace_specInvTop hA hr hrp hgap hpos B

/-- At `r = 1` and under `TopSimple`, `specInvTop A hA 1` is the rank-one matrix
`λ_max⁻¹ v_max v_maxᵀ`, written through its action on a vector. The coefficient set
`topEigSet A hA 1` is the singleton `{λ_max}`, and simplicity makes exactly one eigenvector of
`Matrix.IsHermitian.eigenvectorBasis` carry it. -/
theorem specInvTop_one_mulVec {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (hp : 0 < p)
    (hsimple : TopSimple A hA) (y : Fin p → ℝ) :
    specInvTop A hA 1 *ᵥ y
      = ((lamMax A hA)⁻¹ * (WithLp.ofLp (vMax A hA) ⬝ᵥ y)) • WithLp.ofLp (vMax A hA) := by
  classical
  have hcard : 0 < Fintype.card (Fin p) := by simpa using hp
  set i₀ : Fin p := Fintype.equivOfCardEq (Fintype.card_fin _) ⟨0, hcard⟩ with hi₀
  have hvm : vMax A hA = hA.eigenvectorBasis i₀ := by rw [vMax, dif_pos hp]
  have hlamdef : lamMax A hA = hA.eigenvalues₀ ⟨0, hcard⟩ := by rw [lamMax, dif_pos hp]
  have hev0 : hA.eigenvalues i₀ = lamMax A hA := by
    rw [hlamdef, Matrix.IsHermitian.eigenvalues, hi₀, Equiv.symm_apply_apply]
  have hsetiff : ∀ t : ℝ, t ∈ topEigSet A hA 1 ↔ t = lamMax A hA := by
    intro t
    constructor
    · rintro ⟨k, hk, rfl⟩
      have hk0 : (k : ℕ) = 0 := by omega
      rw [hlamdef]
      exact congrArg _ (Fin.val_injective hk0)
    · rintro rfl
      exact ⟨⟨0, hcard⟩, by norm_num, hlamdef.symm⟩
  have hcoef : ∀ i : Fin p, invCoef A hA 1 i
        = if hA.eigenvalues i = lamMax A hA then (lamMax A hA)⁻¹ else 0 := by
    intro i
    by_cases h : hA.eigenvalues i = lamMax A hA
    · rw [invCoef_of_mem ((hsetiff _).mpr h), if_pos h, h]
    · rw [invCoef_of_notMem (fun hc => h ((hsetiff _).mp hc)), if_neg h]
  have hmemtop : ∀ i : Fin p, hA.eigenvalues i = lamMax A hA →
      hA.eigenvectorBasis i ∈ topSpace A hA := by
    intro i hi
    rw [topSpace_eq_eigenspace A hA]
    refine Module.End.mem_eigenspace_iff.mpr ?_
    rw [← hi]
    apply WithLp.ofLp_injective
    simpa using hA.mulVec_eigenvectorBasis i
  have hnorm : ∀ i : Fin p, ‖hA.eigenvectorBasis i‖ = 1 := fun i =>
    (hA.eigenvectorBasis).orthonormal.1 i
  have hproj := topProj_eq_rankOne hsimple (hmemtop i₀ hev0) (hnorm i₀)
  have huniq : ∀ i : Fin p, hA.eigenvalues i = lamMax A hA → i = i₀ := by
    intro i hi
    by_contra hne
    have h1 : topProj A hA (hA.eigenvectorBasis i) = hA.eigenvectorBasis i :=
      Submodule.starProjection_eq_self_iff.mpr (hmemtop i hi)
    have h2 : topProj A hA (hA.eigenvectorBasis i)
        = ⟪hA.eigenvectorBasis i₀, hA.eigenvectorBasis i⟫_ℝ • hA.eigenvectorBasis i₀ :=
      hproj _
    rw [(hA.eigenvectorBasis).orthonormal.2 (Ne.symm hne), zero_smul] at h2
    have h3 := hnorm i
    rw [← h1, h2, norm_zero] at h3
    exact absurd h3 (by norm_num)
  rw [specInvTop_mulVec, Finset.sum_eq_single i₀]
  · rw [hcoef i₀, if_pos hev0, hvm]
  · intro b _ hb
    rw [hcoef b, if_neg (fun hc => hb (huniq b hc)), zero_mul, zero_smul]
  · intro h
    exact absurd (Finset.mem_univ i₀) h

end TraceFun

/-! ### 2. The unaligned model at `r_i = 1` -/

/-- `assum:unaligned` with `r_i = 1` for every table. `M` tables live on one probability space
per `N`. A shared `V N ∈ ℝ^{d_N × r}` has orthonormal columns, each table carries a fixed unit
vector `R i ∈ ℝ^r`, and the spike of table `i` is `v_i = V R_i`. Cross-table independence is
the separate predicate `JointGaussianNoise`, as for `MultiTableModel`.

`Rank(∑_i R_i R_iᵀ) = r` of the paper is **not** a field. The unweighted proposition below
does not use it (necessity scan row N3), and `thm:gen_rank_weight_svdstak` takes it as an
explicit argument (`thm_gen_rank_weight_svdstak_paper`, `RankR/Weighted.lean`). -/
structure UnalignedModel {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) (r : ℕ) where
  /-- the `M` rank-one spiked tables -/
  tbl : (i : Fin M) → SpikedModel μ (n i) d
  /-- the shared ambient subspace, as a matrix with orthonormal columns -/
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin r) ℝ
  /-- the alignment vectors `R_i ∈ O(r, 1)`, fixed in `N` -/
  R : Fin M → EuclideanSpace ℝ (Fin r)
  hV : ∀ N, (V N)ᵀ * V N = 1
  hR : ∀ i, ‖R i‖ = 1
  hv : ∀ i N, (tbl i).v N = WithLp.toLp 2 ((V N).mulVec (WithLp.ofLp (R i)))

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- Joint law of all noise matrices at level `N`: independent Gaussian tables. The same
predicate as `MultiTableModel.JointGaussianNoise`; the single-table `GaussianNoise` cannot
express independence across tables. -/
def JointGaussianNoise (m : UnalignedModel μ M n d r) : Prop :=
  ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω)
    (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N)

/-! ### 3. The per-table estimators and the svdstack objects -/

/-- Gram matrix `X_iᵀ X_i` of table `i`. -/
noncomputable def tableGram (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ := ((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω

theorem isHermitian_tableGram (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    (m.tableGram i N ω).IsHermitian :=
  isHermitian_transpose_mul_self ((m.tbl i).X N ω)

/-- `v̂_i`, the top right singular vector of table `i`, with the paper's sign convention
`⟪v̂_i, v_i⟫ ≥ 0` (`main_paper.tex:1959`). The sign test reads the unknown `v_i`, so `vhat` is
not an observable estimator; the observable one differs by a sign, which `perfR` squares
away. -/
noncomputable def vhat (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω), (m.tbl i).v N⟫_ℝ then
    vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
  else -vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)

/-- `Ṽ`, the `M × d` matrix whose row `i` is `v̂_iᵀ`. At `r_i = 1` it has `r̃ = M` rows. -/
noncomputable def Vt (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin M) (Fin (d N)) ℝ :=
  Matrix.of fun i k => m.vhat i N ω k

/-- `Ṽ Ṽᵀ`, the `M × M` Gram matrix of the per-table estimates. -/
noncomputable def gram (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin M) (Fin M) ℝ := m.Vt N ω * (m.Vt N ω)ᵀ

theorem isHermitian_gram (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    (m.gram N ω).IsHermitian :=
  isHermitian_mul_transpose_self (m.Vt N ω)

/-- `Ṽ V`, the `M × r` matrix of the overlaps `⟪v̂_i, V e_k⟫`. -/
noncomputable def VtV (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin M) (Fin r) ℝ := m.Vt N ω * m.V N

/-- Performance of svdstack: `tr( (Ṽ V)ᵀ (specInvTop (Ṽ Ṽᵀ) r) (Ṽ V) )`, which equals
`‖V̂_svdstackᵀ V‖_F²` whenever `λ_r(Ṽ Ṽᵀ) > 0`. See the header. -/
noncomputable def perfR (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtV N ω)ᵀ * specInvTop (m.gram N ω) (m.isHermitian_gram N ω) r *
    m.VtV N ω)

/-! ### 3b. Basic facts about `v̂_i`, `v_i` and the columns of `V`

The rank-`r` transcription of the corresponding section of `SVDStack/Gram.lean`. The proofs
are the same; only `v_i = V R_i` replaces the shared `v`. -/

/-- `v̂_i` is a unit vector. -/
theorem norm_vhat (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    ‖m.vhat i N ω‖ = 1 := by
  rw [UnalignedModel.vhat]
  split_ifs with h
  · exact norm_vMax ((m.tbl i).hd N) _ _
  · rw [norm_neg]
    exact norm_vMax ((m.tbl i).hd N) _ _

/-- `v̂_i` lies in the top eigenspace of `X_iᵀ X_i`. -/
theorem mem_topSpace_vhat (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    m.vhat i N ω ∈ topSpace (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      (isHermitian_transpose_mul_self ((m.tbl i).X N ω)) := by
  rw [UnalignedModel.vhat]
  split_ifs with h
  · exact mem_topSpace_vMax ((m.tbl i).hd N) _ _
  · exact Submodule.neg_mem _ (mem_topSpace_vMax ((m.tbl i).hd N) _ _)

/-- The sign convention of `vhat`: `⟪v̂_i, v_i⟫ ≥ 0`. -/
theorem inner_vhat_nonneg (m : UnalignedModel μ M n d r) (i : Fin M) (N : ℕ) (ω : Ω N) :
    0 ≤ ⟪m.vhat i N ω, (m.tbl i).v N⟫_ℝ := by
  rw [UnalignedModel.vhat]
  split_ifs with h
  · exact h
  · have hneg : ⟪-vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω),
        (m.tbl i).v N⟫_ℝ
        = -⟪vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω), (m.tbl i).v N⟫_ℝ := by
      simp
    rw [hneg]
    linarith [not_le.mp h]

/-- The spikes of two tables have the inner product of their alignment vectors:
`⟪V R_i, V R_j⟫ = ⟪R_i, R_j⟫` by `hV`. -/
theorem inner_v_v (m : UnalignedModel μ M n d r) (i j : Fin M) (N : ℕ) :
    ⟪(m.tbl i).v N, (m.tbl j).v N⟫_ℝ = ⟪m.R i, m.R j⟫_ℝ := by
  have hVV := m.hV N
  have key : ∀ a b : Fin r → ℝ,
      ((m.V N).mulVec a) ⬝ᵥ ((m.V N).mulVec b) = a ⬝ᵥ b := by
    intro a b
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec, hVV,
      Matrix.one_mulVec]
  rw [real_inner_eq_dotProduct, real_inner_eq_dotProduct, m.hv i N, m.hv j N]
  exact key _ _

/-- Column `k` of `V N`, as a vector of `EuclideanSpace`. -/
noncomputable def colVec (m : UnalignedModel μ M n d r) (N : ℕ) (k : Fin r) :
    EuclideanSpace ℝ (Fin (d N)) := WithLp.toLp 2 fun l => m.V N l k

/-- `⟪v_i, V e_k⟫ = (R_i)_k`, the paper's `(V R_i)ᵀ V = R_iᵀ` (`main_paper.tex:1964`). -/
theorem inner_v_colVec (m : UnalignedModel μ M n d r) (i : Fin M) (k : Fin r) (N : ℕ) :
    ⟪(m.tbl i).v N, m.colVec N k⟫_ℝ = m.R i k := by
  have hVV := m.hV N
  rw [real_inner_eq_dotProduct, m.hv i N]
  change ∑ l, ((m.V N).mulVec (WithLp.ofLp (m.R i))) l * (m.V N) l k = m.R i k
  calc ∑ l, ((m.V N).mulVec (WithLp.ofLp (m.R i))) l * (m.V N) l k
      = ∑ l, ∑ t, WithLp.ofLp (m.R i) t * ((m.V N) l t * (m.V N) l k) := by
        refine Finset.sum_congr rfl fun l _ => ?_
        rw [Matrix.mulVec, dotProduct, Finset.sum_mul]
        exact Finset.sum_congr rfl fun t _ => by ring
    _ = ∑ t, WithLp.ofLp (m.R i) t * ((m.V N)ᵀ * m.V N) t k := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun t _ => ?_
        rw [Matrix.mul_apply, Finset.mul_sum]
        exact Finset.sum_congr rfl fun l _ => by
          rw [Matrix.transpose_apply]
    _ = m.R i k := by
        rw [hVV]
        simp [Matrix.one_apply]

/-- The columns of `V N` are unit vectors. -/
theorem norm_colVec (m : UnalignedModel μ M n d r) (N : ℕ) (k : Fin r) :
    ‖m.colVec N k‖ = 1 := by
  have hVV := m.hV N
  have hsq : ‖m.colVec N k‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct]
    change ∑ l, (m.V N) l k * (m.V N) l k = 1
    have h : ∑ l, (m.V N) l k * (m.V N) l k = ((m.V N)ᵀ * m.V N) k k := by
      rw [Matrix.mul_apply]
      exact Finset.sum_congr rfl fun l _ => by rw [Matrix.transpose_apply]
    rw [h, hVV]
    simp
  nlinarith [norm_nonneg (m.colVec N k)]

/-- Entry `(i, k)` of `Ṽ V` is `⟪v̂_i, V e_k⟫`. -/
theorem VtV_eq_inner (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) (i : Fin M) (k : Fin r) :
    m.VtV N ω i k = ⟪m.vhat i N ω, m.colVec N k⟫_ℝ := by
  rw [UnalignedModel.VtV, Matrix.mul_apply, real_inner_eq_dotProduct]
  rfl

/-- Entry `(i, j)` of `Ṽ Ṽᵀ` is `⟪v̂_i, v̂_j⟫`. -/
theorem gram_eq_inner (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N) (i j : Fin M) :
    m.gram N ω i j = ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ := by
  rw [UnalignedModel.gram, Matrix.mul_apply, real_inner_eq_dotProduct]
  rfl

/-! ### 3c. The noise-to-data map and the Fubini step

Both facts are the model-free versions of `SVDStack/Gram.lean` at `m.tbl`: the argument reads
only the per-table data and the product law, never the shared `v` of a `MultiTableModel`.
Cleanup wave 2 replaced the 108-line transcription that stood here by these four lines. -/

/-- Table `k` as a function of its noise matrix: `X_k = θ_k u_k v_kᵀ + d^{-1/2} Z`. -/
noncomputable def dataOf (m : UnalignedModel μ M n d r) (k : Fin M) (N : ℕ)
    (Z : Matrix (Fin (n k N)) (Fin (d N)) ℝ) : Matrix (Fin (n k N)) (Fin (d N)) ℝ :=
  (m.tbl k).dataOf N Z

/-- `dataOf` at the model noise is the data matrix of that table. -/
theorem dataOf_eq (m : UnalignedModel μ M n d r) (k : Fin M) (N : ℕ) (ω : Ω N) :
    m.dataOf k N ((m.tbl k).Z N ω) = (m.tbl k).X N ω := rfl

/-- The noise-to-data map of one table is measurable. -/
theorem measurable_dataOf (m : UnalignedModel μ M n d r) (k : Fin M) (N : ℕ) :
    Measurable (m.dataOf k N) := (m.tbl k).measurable_dataOf N

/-- Marginal of the joint law: `JointGaussianNoise` gives every table its own
`SpikedModel.GaussianNoise`. -/
theorem gaussianNoise_of_joint (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    (i : Fin M) : (m.tbl i).GaussianNoise :=
  gaussianNoise_of_joint_pi m.tbl hG i

/-- The Fubini step of `lem:general_rank_delocalization`. The probability that table `i` has a
large overlap with the random direction `delocDir v_i X_j` of table `j` is at most the
supremum of the overlap over deterministic unit directions orthogonal to `v_i`, which
`SingleTableLaw.delocUniform` sends to `0`. Independence across tables enters through `hG`. -/
theorem measure_deloc_le (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) {η : ℝ} (hη : 0 < η) :
    μ N {ω | η ≤ overlap ((m.tbl i).X N ω)
        (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))}
      ≤ ⨆ w ∈ (m.tbl i).orthUnit N, μ N {ω | η ≤ overlap ((m.tbl i).X N ω) w} :=
  measure_deloc_le_of_pi m.tbl hG hij N hη


end UnalignedModel

/-! ### 5. The rank-one reduction

At `r = 1` every `R_i` is the unit vector of `ℝ^1`, so `r̃ = M = r̃` and the objects above are
the rank-one objects of `SVDStack/Defs.lean`. The paper states the same check at
`main_paper.tex:810`. The limit `limitR` needs the top eigenvalue of `A_β` to be simple: at
`β = (0.8, 0)` the matrix `A_β` is the identity, `limitR` is `0.64`, and `svdstackLimit` reads
one arbitrary eigenvector of a degenerate eigenvalue and gives `0` or `0.64`
(`notes/archive/audit_rank_r_2026-08-30.md`, row A1.4c). `thm_svd_stack_general` excludes that
through `hthr` and `abeta_gap`. -/

/-- The unit vector of `ℝ^1`: the value of every `R_i` in the rank-one reduction. -/
noncomputable def oneVec : EuclideanSpace ℝ (Fin 1) := WithLp.toLp 2 fun _ => 1

theorem oneVec_apply (k : Fin 1) : oneVec k = 1 := rfl

theorem inner_oneVec : ⟪oneVec, oneVec⟫_ℝ = 1 := by
  rw [real_inner_eq_dotProduct, dotProduct]
  simp [oneVec]

theorem norm_oneVec : ‖oneVec‖ = 1 := by
  have h := real_inner_self_eq_norm_sq oneVec
  rw [inner_oneVec] at h
  nlinarith [norm_nonneg oneVec]

/-- `B_R` at `r = 1` with every `R_i = 1` is the column vector `β`. -/
theorem BR_reduce {M : ℕ} (β : Fin M → ℝ) (i : Fin M) (k : Fin 1) :
    BR β (fun _ => oneVec) i k = β i := by
  rw [BR, Matrix.of_apply, oneVec_apply, mul_one]

/-- `A_{β,R}` at `r = 1` with every `R_i = 1` is `A_β` (`main_paper.tex:810`). -/
theorem abetaR_reduce {M : ℕ} (β : Fin M → ℝ) : AbetaR β (fun _ => oneVec) = Abeta β := by
  ext i j
  rw [abetaR_apply, inner_oneVec, mul_one]
  simp only [Abeta, Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply]

/-- Rank-one reduction of the limit: at `r = 1` with every `R_i = 1`,
`limitR β R = svdstackLimit β`. The simplicity hypothesis is not removable, see the section
header. Route: `abetaR_reduce` rewrites the matrix; `specInvTop_one_mulVec` collapses
`specInvTop (A_β) 1` to `λ_max⁻¹ v_max v_maxᵀ`; the `1 × 1` trace is then
`(βᵀ v_max)² / λ_max`, which is `svdstackLimit`. -/
theorem limitR_one_eq_svdstackLimit {M : ℕ} (β : Fin M → ℝ)
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β)) :
    limitR β (fun _ => oneVec) = svdstackLimit β := by
  have hM : 0 < M := by
    rcases Nat.eq_zero_or_pos M with h | h
    · exfalso
      subst h
      have h1 : Module.finrank ℝ (topSpace (Abeta β) (isHermitian_Abeta β))
          ≤ Module.finrank ℝ (EuclideanSpace ℝ (Fin 0)) := Submodule.finrank_le _
      rw [hsimple, finrank_euclideanSpace] at h1
      simp at h1
    · exact h
  have hS : specInvTop (AbetaR β fun _ => oneVec) (isHermitian_AbetaR β fun _ => oneVec) 1
      = specInvTop (Abeta β) (isHermitian_Abeta β) 1 :=
    specInvTop_congr_mat _ _ (abetaR_reduce β) 1
  have hmul := specInvTop_one_mulVec (isHermitian_Abeta β) hM hsimple β
  have hsm : ∀ (a : ℝ) (v : Fin M → ℝ), β ⬝ᵥ (a • v) = a * (β ⬝ᵥ v) := by
    intro a v
    simp only [dotProduct, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  have htr : Matrix.trace ((BR β fun _ => oneVec)ᵀ *
        specInvTop (Abeta β) (isHermitian_Abeta β) 1 * (BR β fun _ => oneVec))
      = β ⬝ᵥ (specInvTop (Abeta β) (isHermitian_Abeta β) 1 *ᵥ β) := by
    simp only [Matrix.trace, Matrix.diag_apply, Fin.sum_univ_one, Matrix.mul_apply,
      Matrix.transpose_apply, BR_reduce, Matrix.mulVec, dotProduct, Finset.sum_mul,
      Finset.mul_sum]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun l _ => Finset.sum_congr rfl fun k _ => by ring
  rw [limitR, hS, htr, hmul, hsm, svdstackLimit,
    dotProduct_comm (WithLp.ofLp (vMax (Abeta β) (isHermitian_Abeta β))) β, div_eq_mul_inv]
  ring

end StackedSVD
