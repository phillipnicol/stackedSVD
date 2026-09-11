/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.TopProjPerturb
import StackedSVD.Prob.TendstoInProb

/-!
# The top-`r` spectral projector: definitions and the perturbation statements

STATUS 2026-08-30: proved. No `sorry` and no `axiom`. Task T2 of `notes/RANK_R_PLAN.md`.
The review note is `notes/archive/rank_r_defs.md`; the report is
`notes/archive/agent_reports/rank_r_t2.md`.

This file is the rank-`r` twin of `LinAlg/TopProjPerturb.lean`. That file states everything
at the single index `0` with the gap `λ_1 - λ_2`. Section 7 of the paper needs the same three
facts for the projector on the top `r` eigenvalues.

## Content

1. `topEigSet`, `specTop`, `specProjTop`, `specInvTop`, `TopGap`. The projector is built from
   `Defs.lean`'s `specSpace`, which already takes a **set** of eigenvalues, so no new
   projector is defined. `topEigSet A hA r` is the set `{λ_0, ..., λ_{r-1}}` of the top `r`
   sorted eigenvalues, not an index range. That keeps `specTop` basis free when two
   eigenvalues tie across the boundary; see `notes/archive/rank_r_defs.md`, choice 5.
2. Bridges. `specSpace_eq_spectralSubspace` identifies `specSpace` with StatsMLlib's
   `spectralSubspace`, so that the Davis-Kahan bound applies. `toOp_specInvTop_eigvec` and
   `specProjTop_eigvec` diagonalize `specInvTop` and `specProjTop` in the eigenbasis, which
   gives `A · S = S · A = P` and `S · P = P · S = S` (`toCLM_comp_specInvTop`,
   `specInvTop_comp_toCLM`, `specInvTop_comp_proj`, `proj_comp_specInvTop`). That is the
   bridge item 6 of `notes/archive/audit_rank_r_2026-08-30.md` asks for, and the hard half of T2.
   These, the two coefficients `selCoef` and `invCoef` with their four values, the defining
   sum `specInvTop_eq_sum`, `specInvTop_mulVec` and `specInvTop_congr_mat` are public since
   cleanup wave 1 (2026-08-30). `RankR/Example.lean` (`specInvTop_eq_inv`) and
   `RankR/Defs.lean` (`specInvTop_one_mulVec`) re-derived them while they were private; those
   copies still stand and belong to a later wave.
3. `specProjTop_perturb`: Davis-Kahan for the top-`r` projector, on top of StatsMLlib's
   `LinearMap.IsSymmetric.davisKahan_spectralProjection_hdp`
   (`.lake/packages/StatsMLlib/StatsMLlib/LinearAlgebra/Matrix/Perturbation.lean:713`).
4. `specTop_simple_whp_of_tendsto`: the top-`r` gap holds with probability tending to one for
   a random symmetric matrix whose entries converge in probability to a matrix with a gap.
   The rank-`r` twin of `topSimple_whp_of_tendsto` of `SVDStack/Deterministic.lean`.
5. `norm_specInvTop_sub_lt` and `continuousAt_trace_specInvTop`:
   `(A, B) ↦ tr(Bᵀ (specInvTop A r) B)` is continuous at a matrix with a top-`r` gap and a
   positive `λ_{r-1}`, jointly in the entries of `A` and of `B`. This is the continuous
   mapping input of `TendstoInProbPi.comp_continuous`, in the entrywise `symMat` form of
   `norm_topProj_sq_continuousAt`. `continuousAt_trace_specInvTop_prod` is the product form
   and `abs_trace_conj_sub_le` the Lipschitz bound in the middle matrix.

## `r = p` is allowed

`specProjTop_perturb` and `specTop_simple_whp_of_tendsto` take `r ≤ Fintype.card (Fin p)`,
not `r < ...`. At `r = p` there is no eigenvalue of index `r`, `TopGap` holds for free (the
paper's convention `λ_{p+1} := -∞`), `specProjTop` is the identity (`specProjTop_eq_id`) and
both statements are trivial. The paper's own worked example `eq:psi_equation` has
`r = r̃ = M = 2`, so that case is needed (`notes/archive/audit_rank_r_2026-08-30.md`, finding 5.1).

## The constant of `specProjTop_perturb`

StatsMLlib gives `‖P_J(B) ∘ P_I(A)‖ ≤ ‖A - B‖ / δ` for `δ`-separated sets `I` and `J`. Take
`I = [λ_{r-1}(A), ∞)` and `J = (-∞, λ_r(B)]`, which are separated by `γ - δ > γ/2`, so
`‖(1 - P_top(B)) P_top(A)‖ ≤ 2 δ / γ`, and the mirror bound holds by symmetry. The identity
`Q - P = Q (1 - P) - (1 - Q) P` then gives `4 δ / γ`. The sharp constant for two projectors
of equal rank is `2 δ / γ`, so the stated bound is safe.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix Matrix.Norms.L2Operator

namespace StackedSVD

variable {p : ℕ}

/-! ### 1. The top-`r` eigenvalue set, its projector and the inverse on it -/

/-- The set of the top `r` eigenvalues of a real symmetric matrix. `eigenvalues₀` is antitone,
so the indices `0` to `r - 1` are the `r` largest values. The set, not the index range, is
what `specSpace` consumes; that makes every object below a function of `A` alone. -/
def topEigSet (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) : Set ℝ :=
  {t | ∃ k : Fin (Fintype.card (Fin p)), (k : ℕ) < r ∧ hA.eigenvalues₀ k = t}

/-- Sum of the eigenspaces of the top `r` eigenvalues. At `r = 1` it is `topSpace`. -/
noncomputable def specTop (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) :
    Submodule ℝ (EuclideanSpace ℝ (Fin p)) :=
  specSpace A (topEigSet A hA r)

/-- Orthogonal projector onto `specTop A hA r`. At `r = 1` it is `topProj`. -/
noncomputable def specProjTop (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) :
    EuclideanSpace ℝ (Fin p) →L[ℝ] EuclideanSpace ℝ (Fin p) :=
  specProj A (topEigSet A hA r)

/-- The top-`r` eigengap of `A`: every eigenvalue of index `r` or more is strictly below every
eigenvalue of index below `r`. Since `eigenvalues₀` is antitone this says
`λ_{r-1}(A) > λ_r(A)`, and it is vacuously true when `r = 0` or `r ≥ p`. The second case is
the paper's convention `λ_{r̃+1} := -∞` (`prop:general_rank_unweighted_svdstack`).

**Warning.** At `r ≥ p` no pair `(k, l)` exists, so `TopGap A hA r` holds for **every**
symmetric `A` and carries no eigengap (mechanical audit 2026-08-31, finding 4; probe P3). A
`hgap` hypothesis at `r = M` in `RankR/` therefore says nothing, and the results that take it
stay true for a different reason: at `r = M` the set `topEigSet` holds every eigenvalue,
`specInvTop` is the full inverse, and positive definiteness keeps it continuous. Do not read
`hgap` as a spectral separation unless `r < p`. -/
def TopGap (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) : Prop :=
  ∀ k l : Fin (Fintype.card (Fin p)), (k : ℕ) < r → r ≤ (l : ℕ) →
    hA.eigenvalues₀ l < hA.eigenvalues₀ k

/-- The inverse of `A` on its top-`r` invariant subspace, zero on the orthogonal complement,
as a matrix: `∑_{λ ∈ topEigSet} λ⁻¹ P_λ`. The coefficient reads the eigenvalue and not the
index, so the value does not depend on the choice of eigenvector basis. Under `TopGap A hA r`
and `0 < λ_{r-1}(A)` it is the paper's `Q_r Λ_r^{-1} Q_rᵀ`. -/
noncomputable def specInvTop (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) :
    Matrix (Fin p) (Fin p) ℝ :=
  ∑ i : Fin p,
    Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues i) •
      Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis i))
        (WithLp.ofLp (hA.eigenvectorBasis i))

/-! ### 1b. Bridge lemmas: the eigenbasis, `specSpace` and StatsMLlib's `spectralSubspace`

`specSpace A S` (a supremum of eigenspaces, `Defs.lean`) and
`LinearMap.IsSymmetric.spectralSubspace` (the span of the selected eigenvectors, StatsMLlib)
are the same submodule. StatsMLlib's Davis-Kahan bound is stated for the second one, so the
bridge is what lets `specProjTop_perturb` use it. -/

section Bridge

/-- `Matrix.IsHermitian.eigenvalues₀` is StatsMLlib's operator eigenvalue list. -/
private theorem eigenvalues₀_eq_op {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k = (symmOp hA).eigenvalues finrank_euclideanSpace k := rfl

/-- The `L2` operator norm of a matrix difference is the operator norm of the difference of the
two operators. Repeated from `TopProjPerturb.lean`, where it is private. -/
private theorem opNorm_toCLM_subP (A B : Matrix (Fin p) (Fin p) ℝ) :
    ‖(toOp A - toOp B).toContinuousLinearMap‖ = ‖A - B‖ := by
  rw [show toOp A - toOp B = toOp (A - B) from (map_sub Matrix.toEuclideanLin A B).symm]
  rfl

/-- The matrix eigenbasis of `Defs.lean` diagonalizes the operator. -/
private theorem apply_eigvecP {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) (i : Fin p) :
    toOp A (hA.eigenvectorBasis i) = hA.eigenvalues i • hA.eigenvectorBasis i := by
  apply WithLp.ofLp_injective
  simpa using hA.mulVec_eigenvectorBasis i

/-- `specSpace` (a supremum of eigenspaces) is StatsMLlib's `spectralSubspace` (a span of
eigenvectors). -/
private theorem specSpace_eq_spectralSubspace (A : Matrix (Fin p) (Fin p) ℝ)
    (hA : A.IsHermitian) (S : Set ℝ) :
    specSpace A S = (symmOp hA).spectralSubspace finrank_euclideanSpace S := by
  refine le_antisymm (iSup_le fun t => iSup_le fun ht x hx => ?_) (Submodule.span_le.mpr ?_)
  · rw [Module.End.mem_eigenspace_iff] at hx
    refine (symmOp hA).mem_spectralSubspace_of_eigenvectorBasis_repr_eq_zero
      finrank_euclideanSpace fun j hj => ?_
    rw [OrthonormalBasis.repr_apply_apply]
    refine inner_eq_zero_of_ne hA (t := t) (fun h => ?_) (apply_eigvec hA j) hx
    apply hj
    rw [show (symmOp hA).eigenvalues finrank_euclideanSpace j = t from h]
    exact ht
  · rintro y ⟨j, rfl⟩
    refine (le_iSup₂ (f := fun t (_ : t ∈ S) => Module.End.eigenspace (toOp A) t)
      ((symmOp hA).eigenvalues finrank_euclideanSpace (j : Fin (Fintype.card (Fin p))))
      j.property) ?_
    exact Module.End.mem_eigenspace_iff.mpr
      ((symmOp hA).apply_eigenvectorBasis finrank_euclideanSpace j)

/-- Two eigenvalue sets that select the same eigenvalues give the same spectral subspace. -/
private theorem spectralSubspace_congr_of_iff {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {I J : Set ℝ}
    (h : ∀ k : Fin (Fintype.card (Fin p)), hA.eigenvalues₀ k ∈ I ↔ hA.eigenvalues₀ k ∈ J) :
    (symmOp hA).spectralSubspace finrank_euclideanSpace I
      = (symmOp hA).spectralSubspace finrank_euclideanSpace J := by
  refine le_antisymm (fun x hx => ?_) (fun x hx => ?_) <;>
    refine (symmOp hA).mem_spectralSubspace_of_eigenvectorBasis_repr_eq_zero
      finrank_euclideanSpace fun j hj => ?_ <;>
    refine (symmOp hA).eigenvectorBasis_repr_eq_zero_of_mem_spectralSubspace
      finrank_euclideanSpace ?_ hx
  · exact fun hmem => hj ((h j).mp hmem)
  · exact fun hmem => hj ((h j).mpr hmem)

/-- Two eigenvalue sets that select complementary eigenvalues give orthogonal complements. -/
private theorem spectralSubspace_orthogonal_of_iff {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {I J : Set ℝ}
    (h : ∀ k : Fin (Fintype.card (Fin p)), hA.eigenvalues₀ k ∈ J ↔ hA.eigenvalues₀ k ∉ I) :
    (symmOp hA).spectralSubspace finrank_euclideanSpace J
      = ((symmOp hA).spectralSubspace finrank_euclideanSpace I)ᗮ := by
  refine le_antisymm (Submodule.span_le.mpr ?_) (fun x hx => ?_)
  · rintro y ⟨j, rfl⟩
    refine (Submodule.mem_orthogonal _ _).mpr fun u hu => ?_
    have hzero : ((symmOp hA).eigenvectorBasis finrank_euclideanSpace).repr u j = 0 :=
      (symmOp hA).eigenvectorBasis_repr_eq_zero_of_mem_spectralSubspace finrank_euclideanSpace
        ((h j).mp j.property) hu
    rw [OrthonormalBasis.repr_apply_apply] at hzero
    rw [real_inner_comm]
    exact hzero
  · refine (symmOp hA).mem_spectralSubspace_of_eigenvectorBasis_repr_eq_zero
      finrank_euclideanSpace fun j hj => ?_
    have hI : (symmOp hA).eigenvalues finrank_euclideanSpace j ∈ I := by
      by_contra hnot
      exact hj ((h j).mpr hnot)
    rw [OrthonormalBasis.repr_apply_apply]
    exact (Submodule.mem_orthogonal _ _).mp hx _
      ((symmOp hA).eigenvectorBasis_mem_spectralSubspace finrank_euclideanSpace hI)

end Bridge

/-! ### 1c. `topEigSet` and `TopGap` in terms of the sorted index -/

section IndexSets

/-- `Fin` order against an explicit index, by definitional unfolding of `Fin`'s `≤`. -/
private theorem fin_le_mk {n : ℕ} (k : Fin n) {a : ℕ} (h : a < n) (hk : (k : ℕ) ≤ a) :
    k ≤ (⟨a, h⟩ : Fin n) := hk

/-- `Fin` order against an explicit index, the other way round. -/
private theorem fin_mk_le {n : ℕ} (k : Fin n) {a : ℕ} (h : a < n) (hk : a ≤ (k : ℕ)) :
    (⟨a, h⟩ : Fin n) ≤ k := hk

/-- Under the top-`r` gap the set `topEigSet` selects exactly the indices below `r`. Public
since cleanup wave 3: `LinAlg/Eigen.lean` used to carry a copy named `mem_topEigSet_iff'`. -/
theorem mem_topEigSet_iff {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (hgap : TopGap A hA r) (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k ∈ topEigSet A hA r ↔ (k : ℕ) < r := by
  constructor
  · rintro ⟨j, hj, hjk⟩
    by_contra hk
    have hlt := hgap j k hj (by omega)
    rw [hjk] at hlt
    exact lt_irrefl _ hlt
  · intro hk
    exact ⟨k, hk, rfl⟩

/-- A strict gap between the sorted eigenvalues of index `r - 1` and `r` gives `TopGap`. -/
private theorem topGap_of_gap {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (hr : 0 < r) (hrp : r < Fintype.card (Fin p)) (hrm : r - 1 < Fintype.card (Fin p))
    (h : hA.eigenvalues₀ ⟨r, hrp⟩ < hA.eigenvalues₀ ⟨r - 1, hrm⟩) : TopGap A hA r := by
  intro k l hk hl
  have h1 : hA.eigenvalues₀ ⟨r - 1, hrm⟩ ≤ hA.eigenvalues₀ k :=
    hA.eigenvalues₀_antitone (fin_le_mk k hrm (by omega))
  have h2 : hA.eigenvalues₀ l ≤ hA.eigenvalues₀ ⟨r, hrp⟩ :=
    hA.eigenvalues₀_antitone (fin_mk_le l hrp (by omega))
  linarith

/-- The top-`r` eigenvalues are those at least `λ_{r-1}`. -/
private theorem mem_Ici_iff_mem_topEigSet {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {r : ℕ} (hr : 0 < r) (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r)
    (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k ∈ Set.Ici (hA.eigenvalues₀ ⟨r - 1, hrm⟩) ↔
      hA.eigenvalues₀ k ∈ topEigSet A hA r := by
  rw [mem_topEigSet_iff hA hgap, Set.mem_Ici]
  constructor
  · intro hk
    by_contra hlt
    exact absurd hk (not_le.mpr (hgap ⟨r - 1, hrm⟩ k (by change r - 1 < r; omega) (by omega)))
  · intro hk
    exact hA.eigenvalues₀_antitone (fin_le_mk k hrm (by omega))

/-- The eigenvalues outside the top `r` are those at most `λ_r`. -/
private theorem mem_Iic_iff_notMem_topEigSet {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {r : ℕ} (hrp : r < Fintype.card (Fin p)) (hgap : TopGap A hA r)
    (k : Fin (Fintype.card (Fin p))) :
    hA.eigenvalues₀ k ∈ Set.Iic (hA.eigenvalues₀ ⟨r, hrp⟩) ↔
      hA.eigenvalues₀ k ∉ topEigSet A hA r := by
  rw [mem_topEigSet_iff hA hgap, Set.mem_Iic, not_lt]
  constructor
  · intro hk
    by_contra hlt
    rw [not_le] at hlt
    exact absurd hk (not_le.mpr (hgap k ⟨r, hrp⟩ hlt le_rfl))
  · intro hk
    exact hA.eigenvalues₀_antitone (fin_mk_le k hrp (by omega))

/-- StatsMLlib's spectral projection is our `starProjection`, as a continuous linear map. -/
private theorem toCLM_spectralProjection {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (I : Set ℝ) :
    ((symmOp hA).spectralProjection finrank_euclideanSpace I).toContinuousLinearMap
      = ((symmOp hA).spectralSubspace finrank_euclideanSpace I).starProjection := by
  ext x
  rfl

/-- Equal submodules have equal orthogonal projectors. `Submodule.starProjection` carries an
instance argument that depends on the submodule, so `rw` cannot do this. -/
private theorem starProjection_congr {K L : Submodule ℝ (EuclideanSpace ℝ (Fin p))}
    [K.HasOrthogonalProjection] [L.HasOrthogonalProjection] (h : K = L) :
    K.starProjection = L.starProjection := by
  subst h
  rfl

/-- The projector on the orthogonal complement is `1 - P`. -/
private theorem starProjection_orthogonal_eq
    (K : Submodule ℝ (EuclideanSpace ℝ (Fin p))) :
    Kᗮ.starProjection = ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - K.starProjection := by
  ext u
  simp

/-- `specProjTop` as the `starProjection` on a spectral subspace of StatsMLlib. -/
private theorem specProjTop_eq_starProjection {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} {I : Set ℝ}
    (h : ∀ k : Fin (Fintype.card (Fin p)),
      hA.eigenvalues₀ k ∈ I ↔ hA.eigenvalues₀ k ∈ topEigSet A hA r) :
    specProjTop A hA r
      = ((symmOp hA).spectralSubspace finrank_euclideanSpace I).starProjection := by
  have hK : specSpace A (topEigSet A hA r)
      = (symmOp hA).spectralSubspace finrank_euclideanSpace I := by
    rw [specSpace_eq_spectralSubspace A hA]
    exact spectralSubspace_congr_of_iff hA fun k => (h k).symm
  exact starProjection_congr hK

/-- `1 - specProjTop` as the `starProjection` on the complementary spectral subspace. -/
private theorem one_sub_specProjTop_eq_starProjection {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} {I J : Set ℝ}
    (hI : ∀ k : Fin (Fintype.card (Fin p)),
      hA.eigenvalues₀ k ∈ I ↔ hA.eigenvalues₀ k ∈ topEigSet A hA r)
    (hJ : ∀ k : Fin (Fintype.card (Fin p)),
      hA.eigenvalues₀ k ∈ J ↔ hA.eigenvalues₀ k ∉ topEigSet A hA r) :
    ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop A hA r
      = ((symmOp hA).spectralSubspace finrank_euclideanSpace J).starProjection := by
  have hKL : (symmOp hA).spectralSubspace finrank_euclideanSpace J
      = ((symmOp hA).spectralSubspace finrank_euclideanSpace I)ᗮ :=
    spectralSubspace_orthogonal_of_iff hA fun k => by rw [hJ k, hI k]
  rw [specProjTop_eq_starProjection hA hI, starProjection_congr hKL,
    starProjection_orthogonal_eq]

/-- When `r` is at least the dimension, `TopGap` holds for every matrix: this is the paper's
convention `λ_{p+1} := -∞`. -/
theorem topGap_of_card_le {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (h : Fintype.card (Fin p) ≤ r) : TopGap A hA r := by
  intro k l _ hl
  exact absurd l.isLt (by omega)

/-- When `r` is at least the dimension, the top-`r` invariant subspace is everything. -/
private theorem specTop_eq_top {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (h : Fintype.card (Fin p) ≤ r) : specSpace A (topEigSet A hA r) = ⊤ := by
  rw [specSpace_eq_spectralSubspace A hA]
  refine eq_top_iff.mpr fun x _ => ?_
  refine (symmOp hA).mem_spectralSubspace_of_eigenvectorBasis_repr_eq_zero
    finrank_euclideanSpace fun j hj => ?_
  exact absurd (⟨j, by omega, rfl⟩ : hA.eigenvalues₀ j ∈ topEigSet A hA r) hj

/-- When `r` is at least the dimension, the top-`r` projector is the identity. -/
theorem specProjTop_eq_id {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (h : Fintype.card (Fin p) ≤ r) :
    specProjTop A hA r = ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) := by
  calc specProjTop A hA r = (⊤ : Submodule ℝ (EuclideanSpace ℝ (Fin p))).starProjection :=
        starProjection_congr (specTop_eq_top hA h)
    _ = ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) := Submodule.starProjection_top

end IndexSets

/-! ### 2. Davis-Kahan for the top-`r` projector -/

set_option linter.unusedVariables false in
/-- Davis-Kahan for the top-`r` spectral projector. If `A` has a gap `γ > 0` between the
eigenvalues of index `r - 1` and `r`, and `‖B - A‖ ≤ δ` with `2 δ < γ`, then both matrices
have a top-`r` gap and the two projectors are `4 δ / γ` apart in operator norm. The rank-`r`
twin of `topProj_perturb`.

The gap hypothesis is stated for `r < p` only. At `r = p` there is no eigenvalue of index `r`,
`TopGap` holds for free (the paper's convention `λ_{p+1} := -∞`), both projectors are the
identity and the bound is trivial; that case is the paper's own worked example `r = M = 2`.
-/
theorem specProjTop_perturb {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p))
    {γ δ : ℝ} (hγ : 0 < γ)
    (hgapA : ∀ h : r < Fintype.card (Fin p),
      γ ≤ hA.eigenvalues₀ ⟨r - 1, by omega⟩ - hA.eigenvalues₀ ⟨r, h⟩)
    (hδ : ‖B - A‖ ≤ δ) (hlt : 2 * δ < γ) :
    TopGap A hA r ∧ TopGap B hB r ∧
      ‖specProjTop B hB r - specProjTop A hA r‖ ≤ 4 * δ / γ := by
  have hδ0 : 0 ≤ δ := le_trans (norm_nonneg _) hδ
  by_cases hrlt : r < Fintype.card (Fin p)
  · have hrm : r - 1 < Fintype.card (Fin p) := by omega
    have hgapA' : γ ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩ := hgapA hrlt
    have hAB : ‖A - B‖ ≤ δ := by rwa [norm_sub_rev]
    have hAB0 : (0 : ℝ) ≤ ‖A - B‖ := norm_nonneg _
    have hw1 := abs_le.mp (abs_eigenvalues₀_sub_le A B hA hB ⟨r - 1, hrm⟩)
    have hw2 := abs_le.mp (abs_eigenvalues₀_sub_le A B hA hB ⟨r, hrlt⟩)
    have htA : TopGap A hA r := topGap_of_gap hA hr hrlt hrm (by linarith)
    have htB : TopGap B hB r := topGap_of_gap hB hr hrlt hrm (by linarith)
    refine ⟨htA, htB, ?_⟩
    have hsep : 0 < γ - δ := by linarith
    have hPA : specProjTop A hA r
        = ((symmOp hA).spectralSubspace finrank_euclideanSpace
            (Set.Ici (hA.eigenvalues₀ ⟨r - 1, hrm⟩))).starProjection :=
      specProjTop_eq_starProjection hA (mem_Ici_iff_mem_topEigSet hA hr hrm htA)
    have hPB : specProjTop B hB r
        = ((symmOp hB).spectralSubspace finrank_euclideanSpace
            (Set.Ici (hB.eigenvalues₀ ⟨r - 1, hrm⟩))).starProjection :=
      specProjTop_eq_starProjection hB (mem_Ici_iff_mem_topEigSet hB hr hrm htB)
    have hQA : ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop A hA r
        = ((symmOp hA).spectralSubspace finrank_euclideanSpace
            (Set.Iic (hA.eigenvalues₀ ⟨r, hrlt⟩))).starProjection :=
      one_sub_specProjTop_eq_starProjection hA (mem_Ici_iff_mem_topEigSet hA hr hrm htA)
        (mem_Iic_iff_notMem_topEigSet hA hrlt htA)
    have hQB : ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop B hB r
        = ((symmOp hB).spectralSubspace finrank_euclideanSpace
            (Set.Iic (hB.eigenvalues₀ ⟨r, hrlt⟩))).starProjection :=
      one_sub_specProjTop_eq_starProjection hB (mem_Ici_iff_mem_topEigSet hB hr hrm htB)
        (mem_Iic_iff_notMem_topEigSet hB hrlt htB)
    -- Davis-Kahan, the two mixed products
    have hdk1 : ‖(ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop B hB r).comp
        (specProjTop A hA r)‖ ≤ ‖A - B‖ / (γ - δ) := by
      rw [hQB, hPA]
      have hsepIJ : Set.SeparatedBy (γ - δ) (Set.Ici (hA.eigenvalues₀ ⟨r - 1, hrm⟩))
          (Set.Iic (hB.eigenvalues₀ ⟨r, hrlt⟩)) := by
        intro x hx y hy
        rw [Set.mem_Ici] at hx
        rw [Set.mem_Iic] at hy
        rw [abs_of_nonneg (by linarith)]
        linarith
      have h := (symmOp hA).davisKahan_spectralProjection_hdp (symmOp hB)
        finrank_euclideanSpace (I := Set.Ici (hA.eigenvalues₀ ⟨r - 1, hrm⟩))
        (J := Set.Iic (hB.eigenvalues₀ ⟨r, hrlt⟩)) (δ := γ - δ) Set.ordConnected_Ici hsep hsepIJ
      rwa [opNorm_toCLM_subP, toCLM_spectralProjection hB, toCLM_spectralProjection hA] at h
    have hdk2 : ‖(specProjTop B hB r).comp
        (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop A hA r)‖
          ≤ ‖A - B‖ / (γ - δ) := by
      rw [hPB, hQA]
      have hsepIJ : Set.SeparatedBy (γ - δ) (Set.Iic (hA.eigenvalues₀ ⟨r, hrlt⟩))
          (Set.Ici (hB.eigenvalues₀ ⟨r - 1, hrm⟩)) := by
        intro x hx y hy
        rw [Set.mem_Iic] at hx
        rw [Set.mem_Ici] at hy
        rw [abs_of_nonpos (by linarith), neg_sub]
        linarith
      have h := (symmOp hA).davisKahan_spectralProjection_hdp (symmOp hB)
        finrank_euclideanSpace (I := Set.Iic (hA.eigenvalues₀ ⟨r, hrlt⟩))
        (J := Set.Ici (hB.eigenvalues₀ ⟨r - 1, hrm⟩)) (δ := γ - δ) Set.ordConnected_Iic hsep
        hsepIJ
      rwa [opNorm_toCLM_subP, toCLM_spectralProjection hB, toCLM_spectralProjection hA] at h
    -- `Q - P = Q (1 - P) - (1 - Q) P`
    have key : specProjTop B hB r - specProjTop A hA r
        = (specProjTop B hB r).comp
            (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop A hA r)
          - (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop B hB r).comp
            (specProjTop A hA r) := by
      refine ContinuousLinearMap.ext fun x => ?_
      change specProjTop B hB r x - specProjTop A hA r x
        = specProjTop B hB r (x - specProjTop A hA r x)
          - (specProjTop A hA r x - specProjTop B hB r (specProjTop A hA r x))
      rw [map_sub]
      abel
    have hhalf : ‖A - B‖ / (γ - δ) ≤ 2 * δ / γ := by
      rw [div_le_div_iff₀ hsep hγ]
      nlinarith [mul_le_mul_of_nonneg_right hAB hγ.le,
        mul_nonneg hδ0 (show (0 : ℝ) ≤ γ - 2 * δ by linarith)]
    have hsum : (4 : ℝ) * δ / γ = 2 * δ / γ + 2 * δ / γ := by ring
    rw [key, hsum]
    linarith [norm_sub_le
      ((specProjTop B hB r).comp
        (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop A hA r))
      ((ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin p)) - specProjTop B hB r).comp
        (specProjTop A hA r)), hdk1, hdk2, hhalf]
  · rw [not_lt] at hrlt
    refine ⟨topGap_of_card_le hA hrlt, topGap_of_card_le hB hrlt, ?_⟩
    rw [specProjTop_eq_id hA hrlt, specProjTop_eq_id hB hrlt, sub_self, norm_zero]
    positivity

/-! ### 3. The top-`r` gap holds with probability tending to one -/

section EntrywiseLimit

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
variable {A : Matrix (Fin p) (Fin p) ℝ} {G : (N : ℕ) → Ω N → Matrix (Fin p) (Fin p) ℝ}

/-- The rank-`r` twin of `topSimple_whp_of_tendsto`. If the entries of a random symmetric
matrix converge in probability to those of `A`, and `A` has a top-`r` gap, then `G_N` has a
top-`r` gap with probability tending to one. Route: the entrywise limit gives
`‖G_N - A‖ → 0` in probability (`opNorm_sub_tendstoInProb`), then
`specProjTop_perturb` with `δ < γ / 3`. At `r = p` both sides are trivial: `TopGap` holds for
every matrix. -/
theorem specTop_simple_whp_of_tendsto (hA : A.IsHermitian)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p)) (hgap : TopGap A hA r) :
    Tendsto (fun N => μ N {ω | ¬ TopGap (G N ω) (hsymm N ω) r}) atTop (𝓝 0) := by
  by_cases hrlt : r < Fintype.card (Fin p)
  · have hrm : r - 1 < Fintype.card (Fin p) := by omega
    have hγ : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩ := by
      have h := hgap ⟨r - 1, hrm⟩ ⟨r, hrlt⟩ (show r - 1 < r by omega) (show r ≤ r from le_rfl)
      linarith
    have hop := opNorm_sub_tendstoInProb hA hsymm hconv
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
      (hop ((hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩) / 3) (by positivity))
      (fun _ => zero_le) (fun N => measure_mono fun ω hω => ?_)
    have hω' : ¬ TopGap (G N ω) (hsymm N ω) r := hω
    change (hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩) / 3
      ≤ |‖G N ω - A‖ - 0|
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    by_contra hsmall
    rw [not_le] at hsmall
    obtain ⟨-, hs, -⟩ := specProjTop_perturb hA (hsymm N ω) hr hrp hγ (fun _ => le_rfl)
      hsmall.le (by linarith)
    exact hω' hs
  · rw [not_lt] at hrlt
    have hempty : ∀ N, {ω : Ω N | ¬ TopGap (G N ω) (hsymm N ω) r} = ∅ := by
      intro N
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact topGap_of_card_le (hsymm N ω) hrlt
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds

end EntrywiseLimit


/-! ### 4a. `specInvTop` in the eigenbasis, and its bridge to `specProjTop`

`specProjTop` is built from `specSpace` and `starProjection`; `specInvTop` is an explicit sum
over `Matrix.IsHermitian.eigenvectorBasis`. Nothing in `Defs.lean` connects the two. The
lemmas of this section do: both act diagonally in the eigenbasis, with coefficients
`selCoef = 1{λ ∈ topEigSet}` and `invCoef = 1{λ ∈ topEigSet} λ⁻¹`. Once every selected
eigenvalue is nonzero this gives `specInvTop · A = A · specInvTop = specProjTop` and
`specInvTop · specProjTop = specProjTop · specInvTop = specInvTop`, the bridge that the
continuity proof rests on. -/

section SpecInv

/-- The selection indicator of the top-`r` eigenvalues, at the `i`-th eigenvector. Public
since cleanup wave 1 (2026-08-30): it names the coefficient of `specProjTop`. -/
noncomputable def selCoef (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ)
    (i : Fin p) : ℝ :=
  Set.indicator (topEigSet A hA r) (fun _ => (1 : ℝ)) (hA.eigenvalues i)

/-- The coefficient of `specInvTop` at the `i`-th eigenvector: `λ_i⁻¹` on the top `r`
eigenvalues and `0` elsewhere. Public since cleanup wave 1 (2026-08-30). -/
noncomputable def invCoef (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ)
    (i : Fin p) : ℝ :=
  Set.indicator (topEigSet A hA r) (fun t => t⁻¹) (hA.eigenvalues i)

/-- The selection indicator is `1` at a selected eigenvalue. -/
theorem selCoef_of_mem {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {r : ℕ}
    {i : Fin p} (h : hA.eigenvalues i ∈ topEigSet A hA r) : selCoef A hA r i = 1 :=
  Set.indicator_of_mem h _

/-- The selection indicator is `0` at an eigenvalue outside the top `r`. -/
theorem selCoef_of_notMem {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {r : ℕ}
    {i : Fin p} (h : hA.eigenvalues i ∉ topEigSet A hA r) : selCoef A hA r i = 0 :=
  Set.indicator_of_notMem h _

/-- The `specInvTop` coefficient is `λ_i⁻¹` at a selected eigenvalue. -/
theorem invCoef_of_mem {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {r : ℕ}
    {i : Fin p} (h : hA.eigenvalues i ∈ topEigSet A hA r) :
    invCoef A hA r i = (hA.eigenvalues i)⁻¹ :=
  Set.indicator_of_mem h _

/-- The `specInvTop` coefficient is `0` at an eigenvalue outside the top `r`. -/
theorem invCoef_of_notMem {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {r : ℕ}
    {i : Fin p} (h : hA.eigenvalues i ∉ topEigSet A hA r) : invCoef A hA r i = 0 :=
  Set.indicator_of_notMem h _

/-- `specInvTop` as its defining sum over the eigenbasis, with the coefficients named.
`RankR/Defs.lean` re-derives this as `specInvTop_eq_sum'`. -/
theorem specInvTop_eq_sum (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ) :
    specInvTop A hA r = ∑ i : Fin p, invCoef A hA r i •
      Matrix.vecMulVec (WithLp.ofLp (hA.eigenvectorBasis i))
        (WithLp.ofLp (hA.eigenvectorBasis i)) := rfl

/-- Every `eigenvalues` value is an `eigenvalues₀` value. -/
private theorem exists_eigenvalues₀_eq {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (i : Fin p) : ∃ k : Fin (Fintype.card (Fin p)), hA.eigenvalues i = hA.eigenvalues₀ k :=
  ⟨_, rfl⟩

/-- A selected eigenvalue is at least `λ_{r-1}`. -/
private theorem le_eigenvalues_of_mem {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {r : ℕ} (hr : 0 < r) (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r)
    {i : Fin p} (hmem : hA.eigenvalues i ∈ topEigSet A hA r) :
    hA.eigenvalues₀ ⟨r - 1, hrm⟩ ≤ hA.eigenvalues i := by
  obtain ⟨k, hk⟩ := exists_eigenvalues₀_eq hA i
  have hklt : (k : ℕ) < r := (mem_topEigSet_iff hA hgap k).mp (by rwa [← hk])
  rw [hk]
  exact hA.eigenvalues₀_antitone (fin_le_mk k hrm (by omega))

/-- `specInvTop` acts on a vector by the diagonal coefficients. Public since cleanup wave 1:
`RankR/Defs.lean` re-derived it for `specInvTop_one_mulVec`. -/
theorem specInvTop_mulVec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ)
    (y : Fin p → ℝ) :
    specInvTop A hA r *ᵥ y = ∑ i : Fin p,
      (invCoef A hA r i * (WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ y)) •
        WithLp.ofLp (hA.eigenvectorBasis i) := by
  rw [specInvTop_eq_sum, Matrix.sum_mulVec]
  refine Finset.sum_congr rfl fun i _ => ?_
  funext j
  simp only [Matrix.mulVec, dotProduct, Matrix.smul_apply, Matrix.vecMulVec_apply,
    Pi.smul_apply, smul_eq_mul, Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- `specInvTop` scales the `j`-th eigenvector by its coefficient. Public since cleanup
wave 1 (2026-08-30). -/
theorem toOp_specInvTop_eigvec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (r : ℕ) (j : Fin p) :
    toOp (specInvTop A hA r) (hA.eigenvectorBasis j)
      = invCoef A hA r j • hA.eigenvectorBasis j := by
  apply WithLp.ofLp_injective
  change specInvTop A hA r *ᵥ WithLp.ofLp (hA.eigenvectorBasis j)
    = WithLp.ofLp (invCoef A hA r j • hA.eigenvectorBasis j)
  rw [specInvTop_mulVec]
  have hsmul : WithLp.ofLp (invCoef A hA r j • hA.eigenvectorBasis j)
      = invCoef A hA r j • WithLp.ofLp (hA.eigenvectorBasis j) := rfl
  rw [hsmul, Finset.sum_eq_single j]
  · congr 1
    rw [← real_inner_eq_dotProduct, real_inner_self_eq_norm_sq,
      (hA.eigenvectorBasis).orthonormal.1 j]
    ring
  · intro k _ hk
    rw [← real_inner_eq_dotProduct, (hA.eigenvectorBasis).orthonormal.2 hk, mul_zero, zero_smul]
  · intro h
    exact absurd (Finset.mem_univ j) h

/-- `specProjTop` keeps the `j`-th eigenvector when its eigenvalue is selected and kills it
otherwise. Public since cleanup wave 1 (2026-08-30). -/
theorem specProjTop_eigvec (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian) (r : ℕ)
    (j : Fin p) :
    specProjTop A hA r (hA.eigenvectorBasis j) = selCoef A hA r j • hA.eigenvectorBasis j := by
  by_cases hmem : hA.eigenvalues j ∈ topEigSet A hA r
  · rw [selCoef_of_mem hmem, one_smul]
    refine Submodule.starProjection_eq_self_iff.mpr ?_
    exact Submodule.mem_iSup_of_mem (hA.eigenvalues j)
      (Submodule.mem_iSup_of_mem hmem
        (Module.End.mem_eigenspace_iff.mpr (apply_eigvecP hA j)))
  · rw [selCoef_of_notMem hmem, zero_smul]
    refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero (Submodule.zero_mem _)
      fun y hy => ?_
    rw [sub_zero]
    have hle : specSpace A (topEigSet A hA r) ≤ (ℝ ∙ hA.eigenvectorBasis j)ᗮ := by
      refine iSup_le fun t => iSup_le fun ht z hz => ?_
      rw [Module.End.mem_eigenspace_iff] at hz
      refine Submodule.mem_orthogonal_singleton_iff_inner_right.mpr ?_
      exact inner_eq_zero_of_ne hA (fun h => hmem (by rw [h]; exact ht))
        (apply_eigvecP hA j) hz
    exact Submodule.mem_orthogonal_singleton_iff_inner_right.mp (hle hy)

/-- Two continuous linear maps that agree on the eigenbasis are equal. -/
private theorem clm_ext_eigvec {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    {S T : EuclideanSpace ℝ (Fin p) →L[ℝ] EuclideanSpace ℝ (Fin p)}
    (h : ∀ j, S (hA.eigenvectorBasis j) = T (hA.eigenvectorBasis j)) : S = T := by
  refine ContinuousLinearMap.coe_injective ?_
  refine (hA.eigenvectorBasis).toBasis.ext fun j => ?_
  simpa using h j

/-- `toEuclideanCLM` is `toOp`, as a function. -/
private theorem toCLM_apply (M : Matrix (Fin p) (Fin p) ℝ) (x : EuclideanSpace ℝ (Fin p)) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) M x = toOp M x := rfl

/-- A positive lower bound on the selected eigenvalues bounds the coefficients. -/
private theorem abs_invCoef_le {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian} {r : ℕ}
    (hr : 0 < r) (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r) {L : ℝ}
    (hL : 0 < L) (hpos : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) (i : Fin p) :
    |invCoef A hA r i| ≤ L⁻¹ := by
  by_cases hmem : hA.eigenvalues i ∈ topEigSet A hA r
  · have hLle : L ≤ hA.eigenvalues i :=
      le_trans hpos (le_eigenvalues_of_mem hr hrm hgap hmem)
    rw [invCoef_of_mem hmem, abs_of_nonneg (inv_nonneg.mpr (by linarith))]
    exact inv_anti₀ hL hLle
  · rw [invCoef_of_notMem hmem, abs_zero]
    exact le_of_lt (inv_pos.mpr hL)

/-- The eigenvalue times its coefficient is the selection indicator. -/
theorem eigenvalue_mul_invCoef {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {r : ℕ} (hr : 0 < r) (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r) {L : ℝ}
    (hL : 0 < L) (hpos : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) (i : Fin p) :
    hA.eigenvalues i * invCoef A hA r i = selCoef A hA r i := by
  by_cases hmem : hA.eigenvalues i ∈ topEigSet A hA r
  · have hLle : L ≤ hA.eigenvalues i :=
      le_trans hpos (le_eigenvalues_of_mem hr hrm hgap hmem)
    rw [selCoef_of_mem hmem, invCoef_of_mem hmem]
    exact mul_inv_cancel₀ (by linarith)
  · rw [selCoef_of_notMem hmem, invCoef_of_notMem hmem, mul_zero]

/-- The selection indicator times the coefficient is the coefficient. -/
theorem selCoef_mul_invCoef {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {r : ℕ} (i : Fin p) : selCoef A hA r i * invCoef A hA r i = invCoef A hA r i := by
  by_cases hmem : hA.eigenvalues i ∈ topEigSet A hA r
  · rw [selCoef_of_mem hmem, one_mul]
  · rw [selCoef_of_notMem hmem, invCoef_of_notMem hmem, mul_zero]

end SpecInv

/-! ### 4b. The bridge identities and the operator-norm perturbation of `specInvTop`

The four identities `A · S = S · A = P` and `S · P = P · S = S` (with `S = specInvTop A r`,
`P = specProjTop A r`) are the bridge that item 6 of `notes/archive/audit_rank_r_2026-08-30.md` asks
for. They give the perturbation bound
`‖S_C - S_A‖ ≤ ‖S_C‖ ‖P_C - P_A‖ + ‖S_C‖ ‖C - A‖ ‖S_A‖ + ‖P_C - P_A‖ ‖S_A‖`,
from the identity `S_C - S_A = S_C (1 - P_A) + S_C (A - C) S_A - (1 - P_C) S_A`. -/

section SpecInvPerturb

/-- Operator-norm bound for a map that is diagonal in an orthonormal basis. -/
private theorem opNorm_le_of_diag {ι : Type*} [Fintype ι]
    (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ (Fin p)))
    (T : EuclideanSpace ℝ (Fin p) →L[ℝ] EuclideanSpace ℝ (Fin p)) {c : ι → ℝ} {K : ℝ}
    (hK : 0 ≤ K) (hT : ∀ i, T (b i) = c i • b i) (hc : ∀ i, |c i| ≤ K) : ‖T‖ ≤ K := by
  refine ContinuousLinearMap.opNorm_le_bound _ hK fun x => ?_
  have hx : x = ∑ i, ⟪b i, x⟫_ℝ • b i := (b.sum_repr' x).symm
  have hexp : T x = ∑ i, (c i * ⟪b i, x⟫_ℝ) • b i := by
    calc T x = T (∑ i, ⟪b i, x⟫_ℝ • b i) := by rw [← hx]
      _ = ∑ i, (c i * ⟪b i, x⟫_ℝ) • b i := by
          rw [map_sum]
          exact Finset.sum_congr rfl fun i _ => by
            rw [map_smul, hT i, smul_smul, mul_comm]
  have hinner : ∀ i, ⟪b i, T x⟫_ℝ = c i * ⟪b i, x⟫_ℝ := by
    intro i
    rw [hexp, inner_sum, Finset.sum_eq_single i]
    · rw [real_inner_smul_right, real_inner_self_eq_norm_sq, b.orthonormal.1 i]
      ring
    · intro k _ hk
      rw [real_inner_smul_right, b.orthonormal.2 (Ne.symm hk), mul_zero]
    · intro h
      exact absurd (Finset.mem_univ i) h
  have hnorm : ∀ z : EuclideanSpace ℝ (Fin p), ‖z‖ ^ 2 = ∑ i, ⟪b i, z⟫_ℝ ^ 2 := by
    intro z
    rw [← b.sum_sq_norm_inner_right z]
    exact Finset.sum_congr rfl fun i _ => by rw [Real.norm_eq_abs, sq_abs]
  have h1 : ‖T x‖ ^ 2 = ∑ i, (c i * ⟪b i, x⟫_ℝ) ^ 2 := by
    rw [hnorm (T x)]
    exact Finset.sum_congr rfl fun i _ => by rw [hinner i]
  have h2 : ∑ i, (c i * ⟪b i, x⟫_ℝ) ^ 2 ≤ ∑ i, K ^ 2 * ⟪b i, x⟫_ℝ ^ 2 := by
    refine Finset.sum_le_sum fun i _ => ?_
    have hcc : c i ^ 2 ≤ K ^ 2 := by nlinarith [sq_abs (c i), hc i, abs_nonneg (c i)]
    rw [mul_pow]
    exact mul_le_mul_of_nonneg_right hcc (sq_nonneg _)
  have h3 : ∑ i, K ^ 2 * ⟪b i, x⟫_ℝ ^ 2 = K ^ 2 * ‖x‖ ^ 2 := by
    rw [← Finset.mul_sum, ← hnorm x]
  have hsq : ‖T x‖ ^ 2 ≤ (K * ‖x‖) ^ 2 := by
    rw [mul_pow, h1]
    linarith
  nlinarith [norm_nonneg (T x), mul_nonneg hK (norm_nonneg x), hsq]

variable {A C : Matrix (Fin p) (Fin p) ℝ}

/-- `A · specInvTop A r = specProjTop A r`. Half of the bridge. Public since cleanup wave 1
(2026-08-30): `RankR/Example.lean` re-derived the `r ≥ p` case for `specInvTop_eq_inv`. -/
theorem toCLM_comp_specInvTop (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r)
    (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r) {L : ℝ} (hL : 0 < L)
    (hpos : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) (x : EuclideanSpace ℝ (Fin p)) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) A (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) x)
      = specProjTop A hA r x := by
  have hclm : (Matrix.toEuclideanCLM (𝕜 := ℝ) A).comp
      (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r)) = specProjTop A hA r := by
    refine clm_ext_eigvec hA fun j => ?_
    change toOp A (toOp (specInvTop A hA r) (hA.eigenvectorBasis j))
      = specProjTop A hA r (hA.eigenvectorBasis j)
    rw [toOp_specInvTop_eigvec, map_smul, apply_eigvecP, smul_smul, specProjTop_eigvec,
      ← eigenvalue_mul_invCoef hr hrm hgap hL hpos j,
      mul_comm (invCoef A hA r j) (hA.eigenvalues j)]
  exact DFunLike.congr_fun hclm x

/-- `specInvTop A r · A = specProjTop A r`. The other half of the bridge. Public since
cleanup wave 1 (2026-08-30). -/
theorem specInvTop_comp_toCLM (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r)
    (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r) {L : ℝ} (hL : 0 < L)
    (hpos : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) (x : EuclideanSpace ℝ (Fin p)) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) (Matrix.toEuclideanCLM (𝕜 := ℝ) A x)
      = specProjTop A hA r x := by
  have hclm : (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r)).comp
      (Matrix.toEuclideanCLM (𝕜 := ℝ) A) = specProjTop A hA r := by
    refine clm_ext_eigvec hA fun j => ?_
    change toOp (specInvTop A hA r) (toOp A (hA.eigenvectorBasis j))
      = specProjTop A hA r (hA.eigenvectorBasis j)
    rw [apply_eigvecP, map_smul, toOp_specInvTop_eigvec, smul_smul, specProjTop_eigvec,
      ← eigenvalue_mul_invCoef hr hrm hgap hL hpos j]
  exact DFunLike.congr_fun hclm x

/-- `specInvTop A r · specProjTop A r = specInvTop A r`. Public since cleanup wave 1
(2026-08-30). -/
theorem specInvTop_comp_proj (hA : A.IsHermitian) {r : ℕ}
    (x : EuclideanSpace ℝ (Fin p)) :
    Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) (specProjTop A hA r x)
      = Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) x := by
  have hclm : (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r)).comp (specProjTop A hA r)
      = Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) := by
    refine clm_ext_eigvec hA fun j => ?_
    change toOp (specInvTop A hA r) (specProjTop A hA r (hA.eigenvectorBasis j))
      = toOp (specInvTop A hA r) (hA.eigenvectorBasis j)
    rw [specProjTop_eigvec, map_smul, toOp_specInvTop_eigvec, smul_smul, selCoef_mul_invCoef]
  exact DFunLike.congr_fun hclm x

/-- `specProjTop A r · specInvTop A r = specInvTop A r`. Public since cleanup wave 1
(2026-08-30). -/
theorem proj_comp_specInvTop (hA : A.IsHermitian) {r : ℕ}
    (x : EuclideanSpace ℝ (Fin p)) :
    specProjTop A hA r (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) x)
      = Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) x := by
  have hclm : (specProjTop A hA r).comp (Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r))
      = Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) := by
    refine clm_ext_eigvec hA fun j => ?_
    change specProjTop A hA r (toOp (specInvTop A hA r) (hA.eigenvectorBasis j))
      = toOp (specInvTop A hA r) (hA.eigenvectorBasis j)
    rw [toOp_specInvTop_eigvec, map_smul, specProjTop_eigvec, smul_smul,
      mul_comm (invCoef A hA r j) (selCoef A hA r j), selCoef_mul_invCoef]
  exact DFunLike.congr_fun hclm x

/-- `‖specInvTop A r‖ ≤ L⁻¹` when every selected eigenvalue is at least `L > 0`. -/
private theorem norm_specInvTop_le (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r)
    (hrm : r - 1 < Fintype.card (Fin p)) (hgap : TopGap A hA r) {L : ℝ} (hL : 0 < L)
    (hpos : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) : ‖specInvTop A hA r‖ ≤ L⁻¹ := by
  rw [Matrix.cstar_norm_def]
  refine opNorm_le_of_diag (hA.eigenvectorBasis) _ (le_of_lt (inv_pos.mpr hL))
    (fun i => ?_) (fun i => abs_invCoef_le hr hrm hgap hL hpos i)
  rw [toCLM_apply]
  exact toOp_specInvTop_eigvec A hA r i

/-- The projector is idempotent. -/
private theorem specProjTop_idem (hA : A.IsHermitian) {r : ℕ} (x : EuclideanSpace ℝ (Fin p)) :
    specProjTop A hA r (specProjTop A hA r x) = specProjTop A hA r x :=
  Submodule.starProjection_eq_self_iff.mpr (Submodule.starProjection_apply_mem _ x)

/-- `‖P x‖ ≤ ‖x‖`. -/
private theorem norm_specProjTop_apply_le (hA : A.IsHermitian) {r : ℕ}
    (x : EuclideanSpace ℝ (Fin p)) : ‖specProjTop A hA r x‖ ≤ ‖x‖ :=
  Submodule.norm_starProjection_apply_le _ x

/-- The perturbation bound for `specInvTop`, from the bridge identities. -/
private theorem norm_specInvTop_sub_le (hA : A.IsHermitian) (hC : C.IsHermitian) {r : ℕ}
    (hr : 0 < r) (hrm : r - 1 < Fintype.card (Fin p)) (hgapA : TopGap A hA r)
    (hgapC : TopGap C hC r) {L : ℝ} (hL : 0 < L)
    (hposA : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩) (hposC : L ≤ hC.eigenvalues₀ ⟨r - 1, hrm⟩) :
    ‖specInvTop C hC r - specInvTop A hA r‖
      ≤ L⁻¹ * ‖specProjTop C hC r - specProjTop A hA r‖
        + L⁻¹ * ‖C - A‖ * L⁻¹
        + ‖specProjTop C hC r - specProjTop A hA r‖ * L⁻¹ := by
  have hLinv : (0 : ℝ) ≤ L⁻¹ := le_of_lt (inv_pos.mpr hL)
  have hdP : (0 : ℝ) ≤ ‖specProjTop C hC r - specProjTop A hA r‖ := norm_nonneg _
  have hSAm := norm_specInvTop_le hA hr hrm hgapA hL hposA
  have hSCm := norm_specInvTop_le hC hr hrm hgapC hL hposC
  rw [Matrix.cstar_norm_def, map_sub]
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun y => ?_
  set SA := Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop A hA r) with hSAdef
  set SC := Matrix.toEuclideanCLM (𝕜 := ℝ) (specInvTop C hC r) with hSCdef
  set PA := specProjTop A hA r with hPAdef
  set PC := specProjTop C hC r with hPCdef
  have hnSA : ‖SA‖ ≤ L⁻¹ := by rw [hSAdef, ← Matrix.cstar_norm_def]; exact hSAm
  have hnSC : ‖SC‖ ≤ L⁻¹ := by rw [hSCdef, ← Matrix.cstar_norm_def]; exact hSCm
  have hSAy : ‖SA y‖ ≤ L⁻¹ * ‖y‖ :=
    le_trans (SA.le_opNorm y) (mul_le_mul_of_nonneg_right hnSA (norm_nonneg y))
  set u1 := SC (y - PA y) with hu1
  set u2 := SC (Matrix.toEuclideanCLM (𝕜 := ℝ) A (SA y)
    - Matrix.toEuclideanCLM (𝕜 := ℝ) C (SA y)) with hu2
  set u3 := SA y - PC (SA y) with hu3
  have hsplit : SC y - SA y = u1 + u2 - u3 := by
    rw [hu1, hu2, hu3, map_sub, map_sub, toCLM_comp_specInvTop hA hr hrm hgapA hL hposA y,
      specInvTop_comp_toCLM hC hr hrm hgapC hL hposC (SA y)]
    abel
  have hp1 : ‖u1‖ ≤ L⁻¹ * ‖PC - PA‖ * ‖y‖ := by
    have hrewrite : u1 = SC (PC ((PC - PA) y)) := by
      rw [hu1, ← specInvTop_comp_proj hC (y - PA y)]
      congr 1
      change PC (y - PA y) = PC (PC y - PA y)
      rw [map_sub, map_sub, specProjTop_idem hC y]
    rw [hrewrite]
    calc ‖SC (PC ((PC - PA) y))‖ ≤ ‖SC‖ * ‖PC ((PC - PA) y)‖ := SC.le_opNorm _
      _ ≤ L⁻¹ * (‖PC - PA‖ * ‖y‖) := by
          refine mul_le_mul hnSC ?_ (norm_nonneg _) hLinv
          exact le_trans (norm_specProjTop_apply_le hC _) ((PC - PA).le_opNorm y)
      _ = L⁻¹ * ‖PC - PA‖ * ‖y‖ := by ring
  have hp2 : ‖u2‖ ≤ L⁻¹ * ‖C - A‖ * L⁻¹ * ‖y‖ := by
    have hAC : Matrix.toEuclideanCLM (𝕜 := ℝ) A (SA y)
        - Matrix.toEuclideanCLM (𝕜 := ℝ) C (SA y)
        = Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C) (SA y) := by
      rw [map_sub]
      rfl
    have hnAC : ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C)‖ = ‖C - A‖ := by
      rw [← Matrix.cstar_norm_def, norm_sub_rev]
    rw [hu2, hAC]
    calc ‖SC (Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C) (SA y))‖
        ≤ ‖SC‖ * ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C) (SA y)‖ := SC.le_opNorm _
      _ ≤ L⁻¹ * (‖C - A‖ * (L⁻¹ * ‖y‖)) := by
          refine mul_le_mul hnSC ?_ (norm_nonneg _) hLinv
          calc ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C) (SA y)‖
              ≤ ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C)‖ * ‖SA y‖ :=
                (Matrix.toEuclideanCLM (𝕜 := ℝ) (A - C)).le_opNorm _
            _ ≤ ‖C - A‖ * (L⁻¹ * ‖y‖) := by
                rw [hnAC]
                exact mul_le_mul_of_nonneg_left hSAy (norm_nonneg _)
      _ = L⁻¹ * ‖C - A‖ * L⁻¹ * ‖y‖ := by ring
  have hp3 : ‖u3‖ ≤ ‖PC - PA‖ * L⁻¹ * ‖y‖ := by
    have hz : u3 = (PA - PC) (SA y) := by
      rw [hu3]
      change SA y - PC (SA y) = PA (SA y) - PC (SA y)
      rw [proj_comp_specInvTop hA y]
    rw [hz]
    calc ‖(PA - PC) (SA y)‖ ≤ ‖PA - PC‖ * ‖SA y‖ := (PA - PC).le_opNorm _
      _ = ‖PC - PA‖ * ‖SA y‖ := by rw [norm_sub_rev]
      _ ≤ ‖PC - PA‖ * (L⁻¹ * ‖y‖) := mul_le_mul_of_nonneg_left hSAy hdP
      _ = ‖PC - PA‖ * L⁻¹ * ‖y‖ := by ring
  have htri : ‖SC y - SA y‖ ≤ ‖u1‖ + ‖u2‖ + ‖u3‖ := by
    rw [hsplit]
    calc ‖u1 + u2 - u3‖ ≤ ‖u1 + u2‖ + ‖u3‖ := norm_sub_le _ _
      _ ≤ ‖u1‖ + ‖u2‖ + ‖u3‖ := by linarith [norm_add_le u1 u2]
  change ‖SC y - SA y‖ ≤ _
  calc ‖SC y - SA y‖ ≤ ‖u1‖ + ‖u2‖ + ‖u3‖ := htri
    _ ≤ L⁻¹ * ‖PC - PA‖ * ‖y‖ + L⁻¹ * ‖C - A‖ * L⁻¹ * ‖y‖ + ‖PC - PA‖ * L⁻¹ * ‖y‖ := by
        linarith
    _ = (L⁻¹ * ‖PC - PA‖ + L⁻¹ * ‖C - A‖ * L⁻¹ + ‖PC - PA‖ * L⁻¹) * ‖y‖ := by ring

end SpecInvPerturb

/-! ### 4. Continuity of the performance functional -/

section Continuity

/-- `specInvTop` does not see the Hermitian proof, so it transports along an equality of
matrices. Public since cleanup wave 1 (2026-08-30): `RankR/Defs.lean` re-derived it as
`specInvTop_congr_mat'` and `RankR/Weighted.lean` as `specInvTop_congr_mat''`. -/
theorem specInvTop_congr_mat {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (r : ℕ) : specInvTop A hA r = specInvTop B hB r := by
  subst h
  rfl

/-- `specInvTop` is continuous in the operator norm at a matrix with a top-`r` gap and a
positive eigenvalue of index `r - 1`. This is the deterministic core of
`continuousAt_trace_specInvTop`. -/
theorem norm_specInvTop_sub_lt {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {r : ℕ}
    (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p)) (hgap : TopGap A hA r)
    (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩) {ε : ℝ} (hε : 0 < ε) :
    ∃ η > 0, ∀ (C : Matrix (Fin p) (Fin p) ℝ) (hC : C.IsHermitian), ‖C - A‖ < η →
      ‖specInvTop C hC r - specInvTop A hA r‖ < ε := by
  have hrm : r - 1 < Fintype.card (Fin p) := by omega
  have hpos' : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ := hpos
  have hL : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ / 2 := by linarith
  obtain ⟨κ, hκ0, η₀, hη₀, hkey⟩ :
      ∃ κ : ℝ, 0 ≤ κ ∧ ∃ η₀ : ℝ, 0 < η₀ ∧ ∀ (C : Matrix (Fin p) (Fin p) ℝ)
        (hC : C.IsHermitian), ‖C - A‖ < η₀ →
          TopGap C hC r ∧ ‖specProjTop C hC r - specProjTop A hA r‖ ≤ κ * ‖C - A‖ := by
    by_cases hrlt : r < Fintype.card (Fin p)
    · have hγ : 0 < hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩ := by
        have h := hgap ⟨r - 1, hrm⟩ ⟨r, hrlt⟩ (show r - 1 < r by omega) (show r ≤ r from le_rfl)
        linarith
      refine ⟨4 / (hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩), by positivity,
        (hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩) / 3, by positivity,
        fun C hC hlt => ?_⟩
      obtain ⟨-, hgC, hbound⟩ :=
        specProjTop_perturb (γ := hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩)
          (δ := ‖C - A‖) hA hC hr hrp hγ (fun _ => le_rfl) le_rfl (by linarith)
      refine ⟨hgC, ?_⟩
      calc ‖specProjTop C hC r - specProjTop A hA r‖
          ≤ 4 * ‖C - A‖ / (hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩) := hbound
        _ = 4 / (hA.eigenvalues₀ ⟨r - 1, hrm⟩ - hA.eigenvalues₀ ⟨r, hrlt⟩) * ‖C - A‖ := by ring
    · rw [not_lt] at hrlt
      refine ⟨0, le_rfl, 1, one_pos, fun C hC _ => ⟨topGap_of_card_le hC hrlt, ?_⟩⟩
      rw [specProjTop_eq_id hA hrlt, specProjTop_eq_id hC hrlt, sub_self, norm_zero, zero_mul]
  set L := hA.eigenvalues₀ ⟨r - 1, hrm⟩ / 2 with hLdef
  have hLinv : (0 : ℝ) ≤ L⁻¹ := le_of_lt (inv_pos.mpr hL)
  set M := κ * L⁻¹ + L⁻¹ * L⁻¹ + κ * L⁻¹ with hMdef
  have hM0 : 0 ≤ M := by positivity
  have hden : (0 : ℝ) < M + 1 := by linarith
  refine ⟨min (min η₀ L) (ε / (M + 1)), by positivity, fun C hC hlt => ?_⟩
  have h1 : ‖C - A‖ < η₀ := lt_of_lt_of_le hlt (le_trans (min_le_left _ _) (min_le_left _ _))
  have h2 : ‖C - A‖ < L := lt_of_lt_of_le hlt (le_trans (min_le_left _ _) (min_le_right _ _))
  have h3 : ‖C - A‖ < ε / (M + 1) := lt_of_lt_of_le hlt (min_le_right _ _)
  obtain ⟨hgC, hPbound⟩ := hkey C hC h1
  have hw := abs_le.mp (abs_eigenvalues₀_sub_le A C hA hC ⟨r - 1, hrm⟩)
  have hrev : ‖A - C‖ = ‖C - A‖ := norm_sub_rev A C
  rw [hrev] at hw
  have hposC : L ≤ hC.eigenvalues₀ ⟨r - 1, hrm⟩ := by
    rw [hLdef]
    linarith [hw.2]
  have hposA : L ≤ hA.eigenvalues₀ ⟨r - 1, hrm⟩ := by rw [hLdef]; linarith
  have hb := norm_specInvTop_sub_le hA hC hr hrm hgap hgC hL hposA hposC
  have e1 : L⁻¹ * ‖specProjTop C hC r - specProjTop A hA r‖ ≤ L⁻¹ * (κ * ‖C - A‖) :=
    mul_le_mul_of_nonneg_left hPbound hLinv
  have e2 : ‖specProjTop C hC r - specProjTop A hA r‖ * L⁻¹ ≤ κ * ‖C - A‖ * L⁻¹ :=
    mul_le_mul_of_nonneg_right hPbound hLinv
  have e3 : L⁻¹ * (κ * ‖C - A‖) + L⁻¹ * ‖C - A‖ * L⁻¹ + κ * ‖C - A‖ * L⁻¹ = M * ‖C - A‖ := by
    rw [hMdef]; ring
  have hfin : ‖specInvTop C hC r - specInvTop A hA r‖ ≤ M * ‖C - A‖ := by linarith
  have h4 : ‖C - A‖ * (M + 1) < ε := by
    calc ‖C - A‖ * (M + 1) < ε / (M + 1) * (M + 1) := mul_lt_mul_of_pos_right h3 hden
      _ = ε := by field_simp
  nlinarith [norm_nonneg (C - A), hM0, hfin, h4]

/-- The entries of a matrix are bounded by its `L2` operator norm. -/
private theorem abs_entry_le_opNorm (M : Matrix (Fin p) (Fin p) ℝ) (i j : Fin p) :
    |M i j| ≤ ‖M‖ := by
  have hx1 : ‖(EuclideanSpace.single i (1 : ℝ) : EuclideanSpace ℝ (Fin p))‖ = 1 := by simp
  have hy1 : ‖(EuclideanSpace.single j (1 : ℝ) : EuclideanSpace ℝ (Fin p))‖ = 1 := by simp
  have hval : ⟪(EuclideanSpace.single i (1 : ℝ) : EuclideanSpace ℝ (Fin p)),
      Matrix.toEuclideanCLM (𝕜 := ℝ) M (EuclideanSpace.single j (1 : ℝ))⟫_ℝ = M i j := by
    rw [real_inner_eq_dotProduct]
    simp
  calc |M i j| = |⟪(EuclideanSpace.single i (1 : ℝ) : EuclideanSpace ℝ (Fin p)),
        Matrix.toEuclideanCLM (𝕜 := ℝ) M (EuclideanSpace.single j (1 : ℝ))⟫_ℝ| := by rw [hval]
    _ ≤ ‖(EuclideanSpace.single i (1 : ℝ) : EuclideanSpace ℝ (Fin p))‖ *
          ‖Matrix.toEuclideanCLM (𝕜 := ℝ) M (EuclideanSpace.single j (1 : ℝ))‖ :=
        abs_real_inner_le_norm _ _
    _ ≤ 1 * (‖Matrix.toEuclideanCLM (𝕜 := ℝ) M‖ *
          ‖(EuclideanSpace.single j (1 : ℝ) : EuclideanSpace ℝ (Fin p))‖) := by
        rw [hx1]
        exact mul_le_mul_of_nonneg_left
          ((Matrix.toEuclideanCLM (𝕜 := ℝ) M).le_opNorm _) (by norm_num)
    _ = ‖M‖ := by rw [hy1, ← Matrix.cstar_norm_def]; ring

/-- Entrywise continuity of `specInvTop`, in the `symMat` form. -/
private theorem continuousAt_specInvTop_entry {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩) (k l : Fin p) :
    ContinuousAt (fun w : Fin p × Fin p → ℝ =>
        specInvTop (symMat w) (isHermitian_symMat w) r k l) (fun z => A z.1 z.2) := by
  refine Metric.continuousAt_iff.mpr fun ε hε => ?_
  obtain ⟨η, hη, hkey⟩ := norm_specInvTop_sub_lt hA hr hrp hgap hpos hε
  refine ⟨η / ((p : ℝ) + 1), by positivity, fun {w} hw => ?_⟩
  have hop : ‖symMat w - A‖ < η := by
    have h := norm_symMat_sub_lt hη hw
    rwa [symMat_entries hA] at h
  have hval := hkey (symMat w) (isHermitian_symMat w) hop
  rw [Real.dist_eq, specInvTop_congr_mat (isHermitian_symMat _) hA (symMat_entries hA) r]
  calc |specInvTop (symMat w) (isHermitian_symMat w) r k l - specInvTop A hA r k l|
      = |(specInvTop (symMat w) (isHermitian_symMat w) r - specInvTop A hA r) k l| := by
        rw [Matrix.sub_apply]
    _ ≤ ‖specInvTop (symMat w) (isHermitian_symMat w) r - specInvTop A hA r‖ :=
        abs_entry_le_opNorm _ k l
    _ < ε := hval

/-- `tr(Yᵀ S Y)` written out entrywise. -/
private theorem trace_conj_eq_sum {q : ℕ} (S : Matrix (Fin p) (Fin p) ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Yᵀ * S * Y) = ∑ j, ∑ l, ∑ k, Y k j * S k l * Y l j := by
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Finset.sum_mul]

/-- `ContinuousAt` of a finite sum. -/
private theorem continuousAt_sum {X ι : Type*} [TopologicalSpace X] {s : Finset ι}
    {f : ι → X → ℝ} {x : X} (h : ∀ i ∈ s, ContinuousAt (f i) x) :
    ContinuousAt (fun y => ∑ i ∈ s, f i y) x :=
  tendsto_finsetSum s h

/-- `(A, B) ↦ tr(Bᵀ (specInvTop A r) B)` is continuous at a symmetric `A` with a top-`r` gap
and a positive eigenvalue of index `r - 1`, **jointly** in the entries of `A` and of `B`.
The index type is the disjoint union of the two entry families, which is the shape
`TendstoInProbPi.comp_continuous` takes; the consumer in `RankR/Defs.lean` has a random
`Ṽ V`, so the joint form is the one it needs (`notes/archive/audit_rank_r_2026-08-30.md`, finding
5.2). Both hypotheses on `A` matter: without the gap the top-`r` invariant subspace jumps,
and without `0 < λ_{r-1}` the inverse blows up. -/
theorem continuousAt_trace_specInvTop {q : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    ContinuousAt (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
        Matrix.trace ((Matrix.of fun i k => z (Sum.inr (i, k)))ᵀ *
          specInvTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) r *
          (Matrix.of fun i k => z (Sum.inr (i, k)))))
      (Sum.elim (fun t => A t.1 t.2) (fun t => B t.1 t.2)) := by
  have hinl : Continuous fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
      fun t : Fin p × Fin p => z (Sum.inl t) :=
    continuous_pi fun t => continuous_apply (Sum.inl t)
  have hcontS : ∀ k l : Fin p,
      ContinuousAt (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
        specInvTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) r k l)
        (Sum.elim (fun t => A t.1 t.2) (fun t => B t.1 t.2)) := fun k l =>
    ContinuousAt.comp
      (g := fun w : Fin p × Fin p → ℝ => specInvTop (symMat w) (isHermitian_symMat w) r k l)
      (f := fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
        fun t : Fin p × Fin p => z (Sum.inl t))
      (continuousAt_specInvTop_entry hA hr hrp hgap hpos k l) hinl.continuousAt
  have hrw : (fun z : (Fin p × Fin p) ⊕ (Fin p × Fin q) → ℝ =>
      Matrix.trace ((Matrix.of fun i k => z (Sum.inr (i, k)))ᵀ *
        specInvTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) r *
        (Matrix.of fun i k => z (Sum.inr (i, k)))))
      = fun z => ∑ j, ∑ l, ∑ k, z (Sum.inr (k, j)) *
          specInvTop (symMat fun t => z (Sum.inl t)) (isHermitian_symMat _) r k l *
          z (Sum.inr (l, j)) := by
    funext z
    exact trace_conj_eq_sum _ _
  rw [hrw]
  refine continuousAt_sum fun j _ => continuousAt_sum fun l _ => continuousAt_sum fun k _ => ?_
  exact ((continuous_apply (Sum.inr (k, j))).continuousAt.mul (hcontS k l)).mul
    (continuous_apply (Sum.inr (l, j))).continuousAt

/-- The same continuity on the product of the two entry families. -/
theorem continuousAt_trace_specInvTop_prod {q : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} (hr : 0 < r) (hrp : r ≤ Fintype.card (Fin p))
    (hgap : TopGap A hA r) (hpos : 0 < hA.eigenvalues₀ ⟨r - 1, by omega⟩)
    (B : Matrix (Fin p) (Fin q) ℝ) :
    ContinuousAt (fun z : (Fin p × Fin p → ℝ) × (Fin p × Fin q → ℝ) =>
        Matrix.trace ((Matrix.of fun i k => z.2 (i, k))ᵀ *
          specInvTop (symMat z.1) (isHermitian_symMat z.1) r *
          (Matrix.of fun i k => z.2 (i, k))))
      (fun t => A t.1 t.2, fun t => B t.1 t.2) := by
  have hcontS : ∀ k l : Fin p,
      ContinuousAt (fun z : (Fin p × Fin p → ℝ) × (Fin p × Fin q → ℝ) =>
        specInvTop (symMat z.1) (isHermitian_symMat z.1) r k l)
        (fun t => A t.1 t.2, fun t => B t.1 t.2) := fun k l =>
    ContinuousAt.comp
      (g := fun w : Fin p × Fin p → ℝ => specInvTop (symMat w) (isHermitian_symMat w) r k l)
      (f := fun z : (Fin p × Fin p → ℝ) × (Fin p × Fin q → ℝ) => z.1)
      (continuousAt_specInvTop_entry hA hr hrp hgap hpos k l) continuous_fst.continuousAt
  have hrw : (fun z : (Fin p × Fin p → ℝ) × (Fin p × Fin q → ℝ) =>
      Matrix.trace ((Matrix.of fun i k => z.2 (i, k))ᵀ *
        specInvTop (symMat z.1) (isHermitian_symMat z.1) r *
        (Matrix.of fun i k => z.2 (i, k))))
      = fun z => ∑ j, ∑ l, ∑ k, z.2 (k, j) *
          specInvTop (symMat z.1) (isHermitian_symMat z.1) r k l * z.2 (l, j) := by
    funext z
    exact trace_conj_eq_sum _ _
  rw [hrw]
  refine continuousAt_sum fun j _ => continuousAt_sum fun l _ => continuousAt_sum fun k _ => ?_
  exact ((((continuous_apply (k, j)).comp continuous_snd).continuousAt).mul (hcontS k l)).mul
    (((continuous_apply (l, j)).comp continuous_snd).continuousAt)

/-- The Lipschitz bound in the middle matrix that task T4 also needs: the trace form of
`abs_norm_sq_topProj_sub_le`. -/
theorem abs_trace_conj_sub_le {q : ℕ} (S T : Matrix (Fin p) (Fin p) ℝ)
    (Y : Matrix (Fin p) (Fin q) ℝ) :
    |Matrix.trace (Yᵀ * S * Y) - Matrix.trace (Yᵀ * T * Y)|
      ≤ ‖S - T‖ * ∑ j, ∑ i, Y i j ^ 2 := by
  have hsub : Matrix.trace (Yᵀ * S * Y) - Matrix.trace (Yᵀ * T * Y)
      = Matrix.trace (Yᵀ * (S - T) * Y) := by
    rw [← Matrix.trace_sub]
    congr 1
    rw [Matrix.mul_sub, Matrix.sub_mul]
  rw [hsub]
  have hcol : ∀ j : Fin q, |∑ l, ∑ k, Y k j * (S - T) k l * Y l j|
      ≤ ‖S - T‖ * ∑ i, Y i j ^ 2 := by
    intro j
    set y : EuclideanSpace ℝ (Fin p) := WithLp.toLp 2 (fun i => Y i j) with hy
    have hquad : ⟪y, Matrix.toEuclideanCLM (𝕜 := ℝ) (S - T) y⟫_ℝ
        = ∑ l, ∑ k, Y k j * (S - T) k l * Y l j := by
      rw [real_inner_eq_dotProduct]
      have h1 : WithLp.ofLp (Matrix.toEuclideanCLM (𝕜 := ℝ) (S - T) y)
          = (S - T) *ᵥ WithLp.ofLp y := rfl
      have h2 : WithLp.ofLp y = fun i => Y i j := rfl
      rw [h1, h2]
      simp only [dotProduct, Matrix.mulVec, Finset.mul_sum]
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun l _ => Finset.sum_congr rfl fun k _ => by ring
    have hny : ‖y‖ ^ 2 = ∑ i, Y i j ^ 2 := EuclideanSpace.real_norm_sq_eq y
    rw [← hquad]
    calc |⟪y, Matrix.toEuclideanCLM (𝕜 := ℝ) (S - T) y⟫_ℝ|
        ≤ ‖y‖ * ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (S - T) y‖ := abs_real_inner_le_norm _ _
      _ ≤ ‖y‖ * (‖S - T‖ * ‖y‖) := by
          refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
          rw [Matrix.cstar_norm_def]
          exact (Matrix.toEuclideanCLM (𝕜 := ℝ) (S - T)).le_opNorm _
      _ = ‖S - T‖ * ‖y‖ ^ 2 := by ring
      _ = ‖S - T‖ * ∑ i, Y i j ^ 2 := by rw [hny]
  rw [trace_conj_eq_sum, Finset.mul_sum]
  calc |∑ j, ∑ l, ∑ k, Y k j * (S - T) k l * Y l j|
      ≤ ∑ j, |∑ l, ∑ k, Y k j * (S - T) k l * Y l j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ j, ‖S - T‖ * ∑ i, Y i j ^ 2 := Finset.sum_le_sum fun j _ => hcol j

end Continuity

end StackedSVD
