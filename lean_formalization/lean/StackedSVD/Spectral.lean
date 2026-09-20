/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs
import StatsMLlib.LinearAlgebra.Matrix.Perturbation

/-!
# Spectral lemmas for `overlap`, `gramLamMax` and `TopSimple`

STATUS 2026-08-29: all nine lemmas are proved. `lake env lean -j 3 StackedSVD/Spectral.lean`
exits 0 with no warning. No deferred proof and no `axiom`; every proof uses only `propext`,
`Classical.choice` and `Quot.sound`.

This file proves the nine lemmas of `notes/SERVER_TODO.md` step 5. They are the interface
that Layer 1 uses for `overlap`:

* `overlap_nonneg`, `overlap_le_norm_sq`: `overlap X w ∈ [0, ‖w‖²]` (projector contraction).
* `overlap_eq_inner_sq`, `overlap_ge_inner_sq`: the bridge to the paper's `⟨v̂, w⟩²`.
* `measurable_overlap`, `measurable_overlap₂`, `measurable_gramLamMax`,
  `measurableSet_topSimple`: the measurability facts that `TendstoInProb` and the a.s.
  statements of `SingleTableLaw` need. `measurable_overlap₂` is the joint form in `(X, w)`,
  which item D of `notes/archive/rmt_roadmap.md` needs for the Fubini step of `lem:delocalization`
  (`notes/archive/audit_scope_2026-08-29.md` section 6.7).
* `overlap_reindex`: invariance under a relabeling of rows and columns, used by stacking.

## Proof routes that the file uses

`measurable_gramLamMax`. Weyl's inequality from StatsMLlib
(`LinearMap.IsSymmetric.abs_eigenvalues_sub_le_opNorm`) bounds `|lamMax A - lamMax B|` by the
operator norm of `toEuclideanCLM (A - B)`. The map `C ↦ toEuclideanCLM C` is linear on a
finite-dimensional space, so it is continuous, and `X ↦ Xᵀ * X` is continuous. A squeeze then
gives `Continuous gramLamMax`. The product sigma-algebra on the matrix space is the Borel
sigma-algebra of the product topology (`instBorelSpaceMatrixReal`), so continuity gives
measurability. Mathlib v4.33.0 has no continuity or measurability statement for eigenvalues.

`measurable_overlap` and `measurable_overlap₂`. Route through a pointwise limit, valid because
the Gram matrix is positive semidefinite. Write `A = Xᵀ * X` and `ρ = lamMax A`.

1. In the eigenbasis of `A`, `overlap X w = ∑ i, (if λ i = ρ then ⟪ b i, w ⟫ else 0)²`
   (`normSq_topProj`) and `⟪ w, (toOp A) ^ k w ⟫ = ∑ i, (λ i) ^ k * ⟪ b i, w ⟫²`
   (`quadForm_pow`).
2. If `ρ > 0`, then `0 ≤ λ i ≤ ρ`, so `(λ i / ρ) ^ (2 k) → 1` when `λ i = ρ` and `→ 0`
   otherwise. Hence `ρ⁻¹ ^ (2 k) * ⟪ w, (toOp A) ^ (2 k) w ⟫ → overlap X w`
   (`tendsto_gram_quadForm`). Positive semidefiniteness is what rules out an eigenvalue `-ρ`.
3. If `ρ = 0`, then every `λ i` is `0 = ρ`, so `overlap X w = ‖w‖²`
   (`overlap_eq_norm_sq_of_lamMax_zero`).

Each term of the sequence is measurable, so `measurable_of_tendsto_metrizable'` applies.

`measurableSet_topSimple`. `TopSimple A` says `finrank ℝ (topSpace A) = 1`. The orthogonal
projector is an `IsProj` onto `topSpace A`, so `LinearMap.IsProj.trace` and
`LinearMap.trace_eq_sum_inner` give `∑ j, overlap X eⱼ = finrank ℝ (topSpace (Xᵀ X))` over the
standard basis (`sum_overlap_basis_eq_finrank`). The set is then the preimage of `{1}` under a
finite sum of the measurable functions of `measurable_overlap`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### Transport along an equality of matrices

`Matrix.IsHermitian` is a `Prop`, so two proofs of it for the same matrix are definitionally
equal. The five lemmas below move `lamMax`, `topSpace`, `topProj` and `TopSimple` across an
equality of matrices, for example the identity `Xᵀ X = W₀ + q qᵀ` of item R0 or the Gram
identity of the svdstack estimator. Every proof is one `subst`.

`Spectral.lean` is the home of the group because every consumer imports it. Before cleanup
wave 2 the same five statements had nine copies, under six different names, in `RMT/R5.lean`,
`RMT/R3.lean`, `RMT/Symmetry.lean`, `LinAlg/TopProjPerturb.lean`, `SVDStack/Defs.lean`,
`SVDStack/Gram.lean`, `SVDStack/Main.lean`, `SVDStack/Weighted.lean` and here. The equality
comes first in every signature. -/

section Congr

variable {d : ℕ} {A B : Matrix (Fin d) (Fin d) ℝ}

/-- The real inner product of `EuclideanSpace` is the dot product. Cleanup wave 2 moved it
here from `SVDStack/Defs.lean`, which is downstream of `LinAlg/` and so out of reach of
`LinAlg/Eigen.lean`. -/
theorem real_inner_eq_dotProduct (x y : EuclideanSpace ℝ (Fin d)) :
    ⟪x, y⟫_ℝ = WithLp.ofLp x ⬝ᵥ WithLp.ofLp y := by
  simp [EuclideanSpace.inner_eq_star_dotProduct, dotProduct_comm]

/-- `lamMax` does not see the Hermitian proof; the matrix may be rewritten under it. -/
theorem lamMax_congr (h : A = B) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    lamMax A hA = lamMax B hB := by subst h; rfl

/-- `topSpace` does not see the Hermitian proof. -/
theorem topSpace_congr (h : A = B) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    topSpace A hA = topSpace B hB := by subst h; rfl

/-- `topProj` does not see the Hermitian proof. -/
theorem topProj_congr (h : A = B) (hA : A.IsHermitian) (hB : B.IsHermitian)
    (w : EuclideanSpace ℝ (Fin d)) : topProj A hA w = topProj B hB w := by subst h; rfl

/-- `TopSimple` transports backwards along an equality of matrices. -/
theorem topSimple_congr (h : A = B) (hA : A.IsHermitian) (hB : B.IsHermitian)
    (hs : TopSimple B hB) : TopSimple A hA := by subst h; exact hs

/-- The `Iff` form of `topSimple_congr`, which `rw` and `simp` take. -/
theorem topSimple_congr_iff (h : A = B) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    TopSimple A hA ↔ TopSimple B hB := by subst h; exact Iff.rfl

end Congr

/-! ### Auxiliary spectral facts

These support the four measurability lemmas. Six of them (`specSpace_singleton`, `symmOp`,
`inner_eq_zero_of_ne`, `apply_eigvec`, `norm_sq_eq_sum_inner`, `inner_eq_sum_inner`) are
public since 2026-09-01, because `LinAlg/SpecIdx.lean` needs them. `topProj_eigvec` and
`normSq_topProj`
diagonalize the projector in the eigenbasis; `quadForm_pow` does the same for the
quadratic form of a matrix power; `tendsto_gram_quadForm` is the power limit that makes
`overlap` a pointwise limit of continuous functions. `abs_lamMax_sub_le` is Weyl's
inequality from StatsMLlib, which gives continuity of `gramLamMax`. -/

section SpectralAux

variable {d : ℕ}

theorem specSpace_singleton (A : Matrix (Fin d) (Fin d) ℝ) (t : ℝ) :
    specSpace A {t} = Module.End.eigenspace (toOp A) t := by
  unfold specSpace
  simp

theorem symmOp {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) :
    (toOp A).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hA

theorem inner_eq_zero_of_ne {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    {s t : ℝ} {x y : EuclideanSpace ℝ (Fin d)} (hst : s ≠ t)
    (hx : toOp A x = s • x) (hy : toOp A y = t • y) : ⟪x, y⟫_ℝ = 0 := by
  have h1 : ⟪toOp A x, y⟫_ℝ = ⟪x, toOp A y⟫_ℝ := symmOp hA x y
  rw [hx, hy, real_inner_smul_left, real_inner_smul_right] at h1
  have h2 : (s - t) * ⟪x, y⟫_ℝ = 0 := by linarith
  rcases mul_eq_zero.mp h2 with h | h
  · exact absurd (by linarith [sub_eq_zero.mp h] : s = t) hst
  · exact h

theorem apply_eigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (i : Fin (Fintype.card (Fin d))) :
    toOp A ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i) =
      hA.eigenvalues₀ i • (symmOp hA).eigenvectorBasis finrank_euclideanSpace i :=
  (symmOp hA).apply_eigenvectorBasis finrank_euclideanSpace i

/-- `specSpace A S` is the supremum of the eigenspaces at the points of `S`, so it contains
each one of them. -/
theorem eigenspace_le_specSpace (A : Matrix (Fin d) (Fin d) ℝ) {S : Set ℝ} {t : ℝ}
    (ht : t ∈ S) : Module.End.eigenspace (toOp A) t ≤ specSpace A S :=
  le_iSup₂ (f := fun (t : ℝ) (_ : t ∈ S) => Module.End.eigenspace (toOp A) t) t ht

open scoped Classical in
/-- The `i`-th eigenvector is fixed by the projector on `specSpace A S` when its eigenvalue is
in `S`, and killed otherwise. `topProj_eigvec` is this statement at `S = {lamMax A hA}`. -/
theorem specProj_eigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (S : Set ℝ)
    (i : Fin (Fintype.card (Fin d))) :
    specProj A S ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i) =
      if hA.eigenvalues₀ i ∈ S then
        (symmOp hA).eigenvectorBasis finrank_euclideanSpace i else 0 := by
  have hev := apply_eigvec hA i
  by_cases h : hA.eigenvalues₀ i ∈ S
  · rw [if_pos h]
    refine Submodule.starProjection_eq_self_iff.mpr ?_
    exact eigenspace_le_specSpace A h (Module.End.mem_eigenspace_iff.mpr hev)
  · rw [if_neg h]
    have horth : specSpace A S ≤
        (ℝ ∙ (symmOp hA).eigenvectorBasis finrank_euclideanSpace i)ᗮ := by
      unfold specSpace
      refine iSup₂_le fun t ht => ?_
      intro y hy
      rw [Module.End.mem_eigenspace_iff] at hy
      rw [Submodule.mem_orthogonal]
      intro u hu
      rw [Submodule.mem_span_singleton] at hu
      obtain ⟨a, rfl⟩ := hu
      rw [real_inner_smul_left,
        inner_eq_zero_of_ne hA (fun hc => h (by rw [hc]; exact ht)) hev hy, mul_zero]
    refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero (Submodule.zero_mem _) ?_
    intro y hy
    rw [sub_zero]
    have hy' := horth hy
    rw [Submodule.mem_orthogonal] at hy'
    exact hy' _ (Submodule.mem_span_singleton_self _)

/-- The `i`-th eigenvector is fixed by the top projector when its eigenvalue is the largest,
and killed otherwise. -/
private theorem topProj_eigvec {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (i : Fin (Fintype.card (Fin d))) :
    topProj A hA ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i) =
      if hA.eigenvalues₀ i = lamMax A hA then
        (symmOp hA).eigenvectorBasis finrank_euclideanSpace i else 0 := by
  have hev := apply_eigvec hA i
  by_cases h : hA.eigenvalues₀ i = lamMax A hA
  · rw [if_pos h]
    refine Submodule.starProjection_eq_self_iff.mpr ?_
    rw [specSpace_singleton, Module.End.mem_eigenspace_iff, hev, h]
  · rw [if_neg h]
    refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero (Submodule.zero_mem _) ?_
    intro y hy
    rw [specSpace_singleton, Module.End.mem_eigenspace_iff] at hy
    rw [sub_zero]
    exact inner_eq_zero_of_ne hA h hev hy

theorem norm_sq_eq_sum_inner {ι : Type*} [Fintype ι]
    (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ (Fin d))) (x : EuclideanSpace ℝ (Fin d)) :
    ‖x‖ ^ 2 = ∑ i, ⟪b i, x⟫_ℝ ^ 2 := by
  rw [← b.sum_sq_norm_inner_right x]
  exact Finset.sum_congr rfl fun i _ => by rw [Real.norm_eq_abs, sq_abs]

theorem inner_eq_sum_inner {ι : Type*} [Fintype ι]
    (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ (Fin d))) (x y : EuclideanSpace ℝ (Fin d)) :
    ⟪x, y⟫_ℝ = ∑ i, ⟪b i, x⟫_ℝ * ⟪b i, y⟫_ℝ := by
  rw [← b.sum_inner_mul_inner x y]
  exact Finset.sum_congr rfl fun i _ => by rw [real_inner_comm x (b i)]

open scoped Classical in
/-- The squared norm of the spectral projection at an arbitrary set `S` of eigenvalues, read in
the sorted eigenbasis. `normSq_topProj` is this statement at `S = {lamMax A hA}`. -/
theorem normSq_specProj (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) (S : Set ℝ)
    (w : EuclideanSpace ℝ (Fin d)) :
    ‖specProj A S w‖ ^ 2 =
      ∑ i, (if hA.eigenvalues₀ i ∈ S then
        ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ else 0) ^ 2 := by
  rw [norm_sq_eq_sum_inner ((symmOp hA).eigenvectorBasis finrank_euclideanSpace)]
  refine Finset.sum_congr rfl fun i _ => ?_
  congr 1
  have h : ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, specProj A S w⟫_ℝ
      = ⟪specProj A S ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i), w⟫_ℝ :=
    (Submodule.inner_starProjection_left_eq_right _ _ _).symm
  rw [h, specProj_eigvec hA S i]
  split <;> simp

private theorem normSq_topProj (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (w : EuclideanSpace ℝ (Fin d)) :
    ‖topProj A hA w‖ ^ 2 =
      ∑ i, (if hA.eigenvalues₀ i = lamMax A hA then
        ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ else 0) ^ 2 := by
  rw [show topProj A hA = specProj A {lamMax A hA} from rfl, normSq_specProj A hA]
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases h : hA.eigenvalues₀ i = lamMax A hA
  · rw [if_pos (Set.mem_singleton_iff.mpr h), if_pos h]
  · rw [if_neg (fun hc => h (Set.mem_singleton_iff.mp hc)), if_neg h]

private theorem toOp_pow (A : Matrix (Fin d) (Fin d) ℝ) (k : ℕ) :
    (toOp A) ^ k = toOp (A ^ k) := by
  induction k with
  | zero =>
    refine LinearMap.ext fun x => PiLp.ext fun j => ?_
    simp
  | succ k ih =>
    rw [pow_succ, pow_succ, ih]
    refine LinearMap.ext fun x => PiLp.ext fun j => ?_
    simp [Matrix.toLpLin_apply]

private theorem inner_eigvec_pow (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (k : ℕ) (w : EuclideanSpace ℝ (Fin d)) (i : Fin (Fintype.card (Fin d))) :
    ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, ((toOp A) ^ k) w⟫_ℝ
      = (hA.eigenvalues₀ i) ^ k *
        ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ := by
  induction k with
  | zero => simp
  | succ k ih =>
    have hstep : ((toOp A) ^ (k + 1)) w = toOp A (((toOp A) ^ k) w) := by
      rw [pow_succ']
      rfl
    rw [hstep]
    have hs : ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i,
          toOp A (((toOp A) ^ k) w)⟫_ℝ
        = ⟪toOp A ((symmOp hA).eigenvectorBasis finrank_euclideanSpace i),
          ((toOp A) ^ k) w⟫_ℝ :=
      (symmOp hA _ _).symm
    rw [hs, apply_eigvec hA i, real_inner_smul_left, ih]
    ring

private theorem quadForm_pow (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (k : ℕ) (w : EuclideanSpace ℝ (Fin d)) :
    ⟪w, ((toOp A) ^ k) w⟫_ℝ =
      ∑ i, (hA.eigenvalues₀ i) ^ k *
        ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ ^ 2 := by
  rw [inner_eq_sum_inner ((symmOp hA).eigenvectorBasis finrank_euclideanSpace) w
    (((toOp A) ^ k) w)]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [inner_eigvec_pow A hA k w i]
  ring

private theorem eigenvalues₀_eq_eigenvalues {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian)
    (k : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ k = hA.eigenvalues (Fintype.equivOfCardEq (Fintype.card_fin _) k) := by
  simp [Matrix.IsHermitian.eigenvalues]

private theorem eigenvalues₀_gram_nonneg {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ)
    (k : Fin (Fintype.card (Fin d))) :
    0 ≤ (isHermitian_transpose_mul_self X).eigenvalues₀ k := by
  have hpsd : (Xᵀ * X).PosSemidef := by
    simpa using Matrix.posSemidef_conjTranspose_mul_self X
  rw [eigenvalues₀_eq_eigenvalues]
  exact hpsd.eigenvalues_nonneg _

private theorem eigenvalues₀_le_lamMax (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hd : 0 < d) (i : Fin (Fintype.card (Fin d))) : hA.eigenvalues₀ i ≤ lamMax A hA := by
  rw [lamMax, dif_pos hd]
  exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

private theorem pos_of_index (i : Fin (Fintype.card (Fin d))) : 0 < d := by
  have h1 : (0 : ℕ) < Fintype.card (Fin d) := Nat.lt_of_le_of_lt (Nat.zero_le i.val) i.isLt
  simpa using h1

private theorem gramLamMax_nonneg {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) : 0 ≤ gramLamMax X := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    simp [gramLamMax, lamMax]
  · rw [gramLamMax, lamMax, dif_pos hd]
    exact eigenvalues₀_gram_nonneg X _

private theorem overlap_eq_norm_sq_of_lamMax_zero {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ)
    (w : EuclideanSpace ℝ (Fin d)) (h : gramLamMax X = 0) : overlap X w = ‖w‖ ^ 2 := by
  have hA := isHermitian_transpose_mul_self X
  change ‖topProj (Xᵀ * X) hA w‖ ^ 2 = _
  rw [normSq_topProj, norm_sq_eq_sum_inner ((symmOp hA).eigenvectorBasis finrank_euclideanSpace)]
  refine Finset.sum_congr rfl fun i _ => ?_
  have hd : 0 < d := pos_of_index i
  have hle : hA.eigenvalues₀ i ≤ lamMax (Xᵀ * X) hA := eigenvalues₀_le_lamMax _ hA hd i
  have hge : 0 ≤ hA.eigenvalues₀ i := eigenvalues₀_gram_nonneg X i
  have hzero : lamMax (Xᵀ * X) hA = 0 := h
  rw [if_pos (le_antisymm (hzero ▸ hle) (hzero ▸ hge))]

private theorem tendsto_gram_quadForm {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ)
    (w : EuclideanSpace ℝ (Fin d))
    (hpos : 0 < lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)) :
    Tendsto (fun k : ℕ =>
        ((lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X))⁻¹) ^ (2 * k) *
          ⟪w, ((toOp (Xᵀ * X)) ^ (2 * k)) w⟫_ℝ)
      atTop (𝓝 (overlap X w)) := by
  have hA := isHermitian_transpose_mul_self X
  have hd : 0 < d := by
    rcases Nat.eq_zero_or_pos d with h0 | h0
    · subst h0
      simp [lamMax] at hpos
    · exact h0
  change Tendsto _ atTop (𝓝 (‖topProj (Xᵀ * X) hA w‖ ^ 2))
  rw [normSq_topProj]
  have key : ∀ k : ℕ,
      ((lamMax (Xᵀ * X) hA)⁻¹) ^ (2 * k) * ⟪w, ((toOp (Xᵀ * X)) ^ (2 * k)) w⟫_ℝ
        = ∑ i, (hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA) ^ (2 * k) *
            ⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ ^ 2 := by
    intro k
    rw [quadForm_pow (Xᵀ * X) hA (2 * k) w, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [div_pow, inv_pow, div_eq_mul_inv]
    ring
  simp only [key]
  refine tendsto_finsetSum _ fun i _ => ?_
  by_cases hc : hA.eigenvalues₀ i = lamMax (Xᵀ * X) hA
  · rw [if_pos hc, hc, div_self (ne_of_gt hpos)]
    simp
  · rw [if_neg hc]
    have hge : 0 ≤ hA.eigenvalues₀ i := eigenvalues₀_gram_nonneg X i
    have hle : hA.eigenvalues₀ i ≤ lamMax (Xᵀ * X) hA := eigenvalues₀_le_lamMax _ hA hd i
    have hlt : hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA < 1 :=
      (div_lt_one hpos).mpr (lt_of_le_of_ne hle hc)
    have hnn : 0 ≤ hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA := div_nonneg hge hpos.le
    have hr : Tendsto (fun k : ℕ => (hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA) ^ (2 * k))
        atTop (𝓝 0) := by
      have h2 : ∀ k : ℕ, (hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA) ^ (2 * k)
          = ((hA.eigenvalues₀ i / lamMax (Xᵀ * X) hA) ^ 2) ^ k := fun k => by rw [pow_mul]
      simp only [h2]
      exact tendsto_pow_atTop_nhds_zero_of_lt_one (by positivity) (by nlinarith)
    simpa using hr.mul_const
      (⟪(symmOp hA).eigenvectorBasis finrank_euclideanSpace i, w⟫_ℝ ^ 2)

private instance instBorelSpaceMatrixReal {m k : Type*} [Countable m] [Countable k] :
    BorelSpace (Matrix m k ℝ) := inferInstanceAs (BorelSpace (m → k → ℝ))

private theorem toCLM_eq {d : ℕ} (C : Matrix (Fin d) (Fin d) ℝ) :
    LinearMap.toContinuousLinearMap (toOp C) = Matrix.toEuclideanCLM (𝕜 := ℝ) C := rfl

private theorem abs_lamMax_sub_le {d : ℕ} (A B : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) (hB : B.IsHermitian) :
    |lamMax A hA - lamMax B hB| ≤ ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - B)‖ := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    rw [lamMax, lamMax, dif_neg (by omega), dif_neg (by omega)]
    simp
  · have hcard : 0 < Fintype.card (Fin d) := by simpa using hd
    have hA' : (toOp A).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hA
    have hB' : (toOp B).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hB
    have key := hA'.abs_eigenvalues_sub_le_opNorm hB' finrank_euclideanSpace ⟨0, hcard⟩
    rw [← map_sub, toCLM_eq] at key
    rw [lamMax, lamMax, dif_pos hd, dif_pos hd]
    exact key

private theorem continuous_gramLamMax {n d : ℕ} :
    Continuous (fun X : Matrix (Fin n) (Fin d) ℝ => gramLamMax X) := by
  have hlin : Continuous (fun C : Matrix (Fin d) (Fin d) ℝ =>
      (Matrix.toEuclideanCLM (𝕜 := ℝ) C :
        EuclideanSpace ℝ (Fin d) →L[ℝ] EuclideanSpace ℝ (Fin d))) :=
    LinearMap.continuous_of_finiteDimensional
      { toFun := fun C => Matrix.toEuclideanCLM (𝕜 := ℝ) C
        map_add' := fun a b => map_add _ _ _
        map_smul' := fun c a => map_smul _ _ _ }
  have hgram : Continuous (fun X : Matrix (Fin n) (Fin d) ℝ => Xᵀ * X) :=
    (continuous_id.matrix_transpose).matrix_mul continuous_id
  rw [continuous_iff_continuousAt]
  intro X0
  have hb : ∀ X : Matrix (Fin n) (Fin d) ℝ,
      |gramLamMax X - gramLamMax X0| ≤
        ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖ :=
    fun X => abs_lamMax_sub_le (Xᵀ * X) (X0ᵀ * X0) (isHermitian_transpose_mul_self X)
      (isHermitian_transpose_mul_self X0)
  have hcont : Continuous (fun X : Matrix (Fin n) (Fin d) ℝ =>
      ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖) :=
    (hlin.comp (hgram.sub continuous_const)).norm
  have h1 : Tendsto (fun X : Matrix (Fin n) (Fin d) ℝ =>
      ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖) (𝓝 X0) (𝓝 0) := by
    have h := hcont.continuousAt (x := X0)
    simpa [ContinuousAt] using h
  have h2 : Tendsto (fun X : Matrix (Fin n) (Fin d) ℝ => |gramLamMax X - gramLamMax X0|)
      (𝓝 X0) (𝓝 0) := squeeze_zero (fun X => abs_nonneg _) hb h1
  have h3 : Tendsto (fun X : Matrix (Fin n) (Fin d) ℝ => gramLamMax X - gramLamMax X0)
      (𝓝 X0) (𝓝 0) :=
    tendsto_zero_iff_norm_tendsto_zero.2 (by simpa [Real.norm_eq_abs] using h2)
  rw [ContinuousAt, ← tendsto_sub_nhds_zero_iff]
  exact h3

private theorem continuous_gram_pow {n d : ℕ} (m : ℕ) :
    Continuous (fun X : Matrix (Fin n) (Fin d) ℝ => (Xᵀ * X) ^ m) := by
  induction m with
  | zero =>
    simp only [pow_zero]
    exact continuous_const
  | succ m ih =>
    simp only [pow_succ]
    exact ih.matrix_mul ((continuous_id.matrix_transpose).matrix_mul continuous_id)

private theorem continuous_toEuclideanCLM {d : ℕ} :
    Continuous (fun C : Matrix (Fin d) (Fin d) ℝ =>
      (Matrix.toEuclideanCLM (𝕜 := ℝ) C :
        EuclideanSpace ℝ (Fin d) →L[ℝ] EuclideanSpace ℝ (Fin d))) :=
  LinearMap.continuous_of_finiteDimensional
    { toFun := fun C => Matrix.toEuclideanCLM (𝕜 := ℝ) C
      map_add' := fun _ _ => map_add _ _ _
      map_smul' := fun _ _ => map_smul _ _ _ }

private instance instSecondCountableMatrixReal {m k : Type*} [Countable m] [Countable k] :
    SecondCountableTopology (Matrix m k ℝ) :=
  inferInstanceAs (SecondCountableTopology (m → k → ℝ))

private theorem sum_overlap_basis_eq_finrank {n d : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) :
    ∑ j, overlap X (EuclideanSpace.basisFun (Fin d) ℝ j)
      = (Module.finrank ℝ (topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) : ℝ) := by
  set K := topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) with hK
  have hproj : LinearMap.IsProj K
      (K.starProjection : EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d)) :=
    { map_mem := fun x => Submodule.coe_mem _
      map_id := fun x hx => Submodule.starProjection_eq_self_iff.mpr hx }
  rw [← hproj.trace, LinearMap.trace_eq_sum_inner _ (EuclideanSpace.basisFun (Fin d) ℝ)]
  refine Finset.sum_congr rfl fun j _ => ?_
  set v := (EuclideanSpace.basisFun (Fin d) ℝ) j with hv
  have hidem : K.starProjection (K.starProjection v) = K.starProjection v :=
    Submodule.starProjection_eq_self_iff.mpr (Submodule.coe_mem _)
  have h1 : ⟪K.starProjection v, K.starProjection v⟫_ℝ = ⟪v, K.starProjection v⟫_ℝ := by
    rw [Submodule.inner_starProjection_left_eq_right, hidem]
  change ‖K.starProjection v‖ ^ 2 = _
  rw [← real_inner_self_eq_norm_sq, h1]
  rfl

end SpectralAux

section Spectral

variable {d : ℕ}

variable {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) (w : EuclideanSpace ℝ (Fin d))

/-! ### Range of `overlap` -/

/-- `overlap` is a squared norm, so it is nonnegative. -/
theorem overlap_nonneg : 0 ≤ overlap X w := by
  unfold overlap
  positivity

/-- An orthogonal projector is a contraction, so `overlap X w ≤ ‖w‖²`. With `‖w‖ = 1` this is
the bound `|⟨v̂, w⟩|² ≤ 1` of the paper. -/
theorem overlap_le_norm_sq : overlap X w ≤ ‖w‖ ^ 2 := by
  have h : ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ≤ ‖w‖ :=
    Submodule.norm_starProjection_apply_le _ w
  have h0 : (0 : ℝ) ≤ ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ := norm_nonneg _
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2 ≤ ‖w‖ ^ 2
  nlinarith

/-- Helper: on a one-dimensional subspace with unit generator `e`, the orthogonal projection
is `u ↦ ⟪e, u⟫ • e`. -/
private theorem starProjection_of_finrank_one {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] [FiniteDimensional ℝ E] (K : Submodule ℝ E)
    (h1 : Module.finrank ℝ K = 1) {e : E} (he : e ∈ K) (hne : ‖e‖ = 1) (u : E) :
    K.starProjection u = ⟪e, u⟫_ℝ • e := by
  have he0 : e ≠ 0 := by
    intro h
    rw [h, norm_zero] at hne
    exact zero_ne_one hne
  have hspan : (ℝ ∙ e) = K :=
    Submodule.eq_of_le_of_finrank_eq ((Submodule.span_singleton_le_iff_mem _ _).2 he)
      (by rw [finrank_span_singleton he0, h1])
  refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero (K.smul_mem _ he) ?_
  intro y hy
  rw [← hspan, Submodule.mem_span_singleton] at hy
  obtain ⟨a, rfl⟩ := hy
  have hc : ⟪u, e⟫_ℝ = ⟪e, u⟫_ℝ := real_inner_comm e u
  have key : ⟪u - ⟪e, u⟫_ℝ • e, e⟫_ℝ = 0 := by
    rw [inner_sub_left, real_inner_smul_left, real_inner_self_eq_norm_sq, hne, hc]
    ring
  rw [real_inner_smul_right, key, mul_zero]

/-! ### Bridge to the paper's `⟨v̂, w⟩²` -/

/-- If the top eigenvalue of `Xᵀ X` is simple, then the projector form and the paper's
inner-product form agree: `overlap X w = ⟪e, w⟫²` for any unit top eigenvector `e`. The sign
ambiguity of `e` is absorbed by the square. -/
theorem overlap_eq_inner_sq
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {e : EuclideanSpace ℝ (Fin d)}
    (he : e ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) (hne : ‖e‖ = 1) :
    overlap X w = ⟪e, w⟫_ℝ ^ 2 := by
  have h := starProjection_of_finrank_one
    (topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) hsimple he hne w
  change ‖(topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)).starProjection w‖ ^ 2 = _
  rw [h, norm_smul, hne, mul_one, Real.norm_eq_abs, sq_abs]

/-- Without simplicity, the projector form still dominates the inner-product form. This
transfers every `overlap → 0` statement to an arbitrary selected top eigenvector. -/
theorem overlap_ge_inner_sq
    {e : EuclideanSpace ℝ (Fin d)}
    (he : e ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) (hne : ‖e‖ = 1) :
    ⟪e, w⟫_ℝ ^ 2 ≤ overlap X w := by
  set K := topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) with hK
  have h1 : ⟪e, w⟫_ℝ = ⟪e, K.starProjection w⟫_ℝ := by
    rw [← Submodule.inner_starProjection_left_eq_right,
      Submodule.starProjection_eq_self_iff.mpr he]
  have h2 : |⟪e, w⟫_ℝ| ≤ ‖K.starProjection w‖ := by
    rw [h1]
    calc |⟪e, K.starProjection w⟫_ℝ| ≤ ‖e‖ * ‖K.starProjection w‖ := abs_real_inner_le_norm _ _
      _ = ‖K.starProjection w‖ := by rw [hne, one_mul]
  change _ ≤ ‖K.starProjection w‖ ^ 2
  nlinarith [abs_nonneg (⟪e, w⟫_ℝ), sq_abs (⟪e, w⟫_ℝ)]

/-! ### Measurability -/

/-- `X ↦ overlap X w` is Borel measurable on the matrix space (product sigma-algebra,
`instMeasurableSpaceMatrix`). See the module doc-string for the route. -/
theorem measurable_overlap :
    Measurable (fun X : Matrix (Fin n) (Fin d) ℝ => overlap X w) := by
  classical
  have hset : MeasurableSet {X : Matrix (Fin n) (Fin d) ℝ |
      lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) = 0} :=
    continuous_gramLamMax.measurable (measurableSet_singleton (0 : ℝ))
  have hmeas : ∀ k : ℕ, Measurable (fun X : Matrix (Fin n) (Fin d) ℝ =>
      if lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) = 0 then ‖w‖ ^ 2
      else ((lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X))⁻¹) ^ (2 * k) *
        ⟪w, ((toOp (Xᵀ * X)) ^ (2 * k)) w⟫_ℝ) := by
    intro k
    refine Measurable.ite hset measurable_const ?_
    refine (continuous_gramLamMax.measurable.inv.pow_const _).mul ?_
    have hc : Continuous (fun X : Matrix (Fin n) (Fin d) ℝ =>
        ⟪w, (Matrix.toEuclideanCLM (𝕜 := ℝ) ((Xᵀ * X) ^ (2 * k))) w⟫_ℝ) :=
      continuous_const.inner
        ((continuous_toEuclideanCLM.comp (continuous_gram_pow (2 * k))).clm_apply
          continuous_const)
    have heq : (fun X : Matrix (Fin n) (Fin d) ℝ => ⟪w, ((toOp (Xᵀ * X)) ^ (2 * k)) w⟫_ℝ)
        = fun X : Matrix (Fin n) (Fin d) ℝ =>
          ⟪w, (Matrix.toEuclideanCLM (𝕜 := ℝ) ((Xᵀ * X) ^ (2 * k))) w⟫_ℝ := by
      funext X
      rw [toOp_pow]
      rfl
    rw [heq]
    exact hc.measurable
  refine measurable_of_tendsto_metrizable' atTop hmeas ?_
  rw [tendsto_pi_nhds]
  intro X
  by_cases h : lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) = 0
  · simp only [if_pos h]
    rw [overlap_eq_norm_sq_of_lamMax_zero X w h]
    exact tendsto_const_nhds
  · simp only [if_neg h]
    exact tendsto_gram_quadForm X w (lt_of_le_of_ne (gramLamMax_nonneg X) (Ne.symm h))

/-- Joint measurability of `overlap` in the matrix and in the test vector. `measurable_overlap`
fixes `w`, which is not enough for item D of `notes/archive/rmt_roadmap.md`: the Fubini step of
`lem:delocalization` integrates `overlap (X_i ω₁) (v̂_j ω₂)`, where the test vector is itself
random. The product sigma-algebra on the pair is `instMeasurableSpaceMatrix` times
`WithLp.measurableSpace`. Route: the same power limit as `measurable_overlap`, now applied to
the continuous map `(X, w) ↦ (ρ⁻¹ • Xᵀ X) ^ t w`. -/
theorem measurable_overlap₂ :
    Measurable (fun p : Matrix (Fin n) (Fin d) ℝ × EuclideanSpace ℝ (Fin d) =>
      overlap p.1 p.2) := by
  classical
  have hset : MeasurableSet {p : Matrix (Fin n) (Fin d) ℝ × EuclideanSpace ℝ (Fin d) |
      lamMax (p.1ᵀ * p.1) (isHermitian_transpose_mul_self p.1) = 0} :=
    (continuous_gramLamMax.measurable.comp measurable_fst) (measurableSet_singleton (0 : ℝ))
  have hmeas : ∀ k : ℕ, Measurable
      (fun p : Matrix (Fin n) (Fin d) ℝ × EuclideanSpace ℝ (Fin d) =>
        if lamMax (p.1ᵀ * p.1) (isHermitian_transpose_mul_self p.1) = 0 then ‖p.2‖ ^ 2
        else ((lamMax (p.1ᵀ * p.1) (isHermitian_transpose_mul_self p.1))⁻¹) ^ (2 * k) *
          ⟪p.2, ((toOp (p.1ᵀ * p.1)) ^ (2 * k)) p.2⟫_ℝ) := by
    intro k
    refine Measurable.ite hset (continuous_snd.norm.pow 2).measurable ?_
    refine (((continuous_gramLamMax.measurable.comp measurable_fst).inv).pow_const _).mul ?_
    have hc : Continuous (fun p : Matrix (Fin n) (Fin d) ℝ × EuclideanSpace ℝ (Fin d) =>
        ⟪p.2, (Matrix.toEuclideanCLM (𝕜 := ℝ) ((p.1ᵀ * p.1) ^ (2 * k))) p.2⟫_ℝ) :=
      continuous_snd.inner
        (((continuous_toEuclideanCLM.comp
          ((continuous_gram_pow (2 * k)).comp continuous_fst))).clm_apply continuous_snd)
    have heq : (fun p : Matrix (Fin n) (Fin d) ℝ × EuclideanSpace ℝ (Fin d) =>
          ⟪p.2, ((toOp (p.1ᵀ * p.1)) ^ (2 * k)) p.2⟫_ℝ)
        = fun p => ⟪p.2, (Matrix.toEuclideanCLM (𝕜 := ℝ) ((p.1ᵀ * p.1) ^ (2 * k))) p.2⟫_ℝ := by
      funext p
      rw [toOp_pow]
      rfl
    rw [heq]
    exact hc.measurable
  refine measurable_of_tendsto_metrizable' atTop hmeas ?_
  rw [tendsto_pi_nhds]
  rintro ⟨X, w⟩
  by_cases h : lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) = 0
  · simp only [if_pos h]
    rw [overlap_eq_norm_sq_of_lamMax_zero X w h]
    exact tendsto_const_nhds
  · simp only [if_neg h]
    exact tendsto_gram_quadForm X w (lt_of_le_of_ne (gramLamMax_nonneg X) (Ne.symm h))

/-- `X ↦ λ_max (Xᵀ X)` is Borel measurable. Route: continuity by the Rayleigh
characterization. -/
theorem measurable_gramLamMax :
    Measurable (fun X : Matrix (Fin n) (Fin d) ℝ => gramLamMax X) :=
  continuous_gramLamMax.measurable

/-- The set of matrices whose Gram matrix has a simple top eigenvalue is measurable. This is
what makes the `topSimple` field of `SingleTableLaw` an a.s. statement about a measurable set,
and milestone S then shows its complement is null. -/
theorem measurableSet_topSimple :
    MeasurableSet {X : Matrix (Fin n) (Fin d) ℝ |
      TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X)} := by
  have hset : {X : Matrix (Fin n) (Fin d) ℝ |
      TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X)}
      = (fun X => ∑ j, overlap X (EuclideanSpace.basisFun (Fin d) ℝ j)) ⁻¹' {1} := by
    ext X
    simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_singleton_iff,
      sum_overlap_basis_eq_finrank]
    exact ⟨fun h => by rw [TopSimple] at h; rw [h]; norm_num,
      fun h => by exact_mod_cast h⟩
  rw [hset]
  exact (Finset.measurable_sum _ fun j _ => measurable_overlap _) (measurableSet_singleton (1 : ℝ))

/-! ### Relabeling of rows and columns -/

/-- The operator of a reindexed matrix is conjugate to the original by the coordinate
permutation isometry. -/
private theorem toOp_reindex {d d' : ℕ} (f : Fin d ≃ Fin d') (A : Matrix (Fin d) (Fin d) ℝ)
    (x : EuclideanSpace ℝ (Fin d')) :
    toOp (Matrix.reindex f f A) x =
      LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f
        (toOp A ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm x)) := by
  apply PiLp.ext
  intro j'
  simp only [Matrix.toLpLin_apply, LinearIsometryEquiv.piLpCongrLeft_apply,
    LinearIsometryEquiv.piLpCongrLeft_symm, Equiv.piCongrLeft'_apply, Equiv.symm_symm,
    Matrix.reindex_apply, Matrix.submatrix_apply, Matrix.mulVec, dotProduct]
  exact Fintype.sum_equiv f.symm _ _ (fun j => by simp)

private theorem eigenspace_reindex {d d' : ℕ} (f : Fin d ≃ Fin d') (A : Matrix (Fin d) (Fin d) ℝ)
    (t : ℝ) :
    Module.End.eigenspace (toOp (Matrix.reindex f f A)) t =
      (Module.End.eigenspace (toOp A) t).map
        ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).toLinearEquiv :
          EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d')) := by
  rw [Submodule.map_equiv_eq_comap_symm]
  ext x
  rw [Submodule.mem_comap, Module.End.mem_eigenspace_iff, Module.End.mem_eigenspace_iff,
    toOp_reindex]
  constructor
  · intro h
    change (toOp A) ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm x)
      = t • (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm x
    have h2 := congrArg (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm h
    rwa [LinearIsometryEquiv.symm_apply_apply, map_smul] at h2
  · intro h
    have h' : (toOp A) ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm x)
        = t • (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).symm x := h
    rw [h', map_smul, LinearIsometryEquiv.apply_symm_apply]

private theorem specSpace_reindex {d d' : ℕ} (f : Fin d ≃ Fin d') (A : Matrix (Fin d) (Fin d) ℝ)
    (S : Set ℝ) :
    specSpace (Matrix.reindex f f A) S =
      (specSpace A S).map ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).toLinearEquiv :
        EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d')) := by
  unfold specSpace
  simp only [Submodule.map_iSup, eigenspace_reindex]

/-- `lamMax` is the greatest point of the real spectrum. -/
private theorem lamMax_isGreatest {d : ℕ} (hd : 0 < d) (A : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) : IsGreatest (spectrum ℝ A) (lamMax A hA) := by
  have hcard : 0 < Fintype.card (Fin d) := by simpa using hd
  have hrange : Set.range hA.eigenvalues = Set.range hA.eigenvalues₀ := by
    ext y
    constructor
    · rintro ⟨i, rfl⟩
      exact ⟨_, rfl⟩
    · rintro ⟨k, rfl⟩
      exact ⟨Fintype.equivOfCardEq (Fintype.card_fin _) k, by
        simp [Matrix.IsHermitian.eigenvalues]⟩
  constructor
  · rw [hA.spectrum_real_eq_range_eigenvalues, hrange, lamMax, dif_pos hd]
    exact ⟨_, rfl⟩
  · rintro x hx
    rw [hA.spectrum_real_eq_range_eigenvalues, hrange] at hx
    obtain ⟨k, rfl⟩ := hx
    rw [lamMax, dif_pos hd]
    exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

private theorem lamMax_reindex {d d' : ℕ} (f : Fin d ≃ Fin d') (A : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) (hB : (Matrix.reindex f f A).IsHermitian) :
    lamMax (Matrix.reindex f f A) hB = lamMax A hA := by
  have hdd : d = d' := by simpa using Fintype.card_congr f
  subst hdd
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    rw [lamMax, lamMax, dif_neg (by omega), dif_neg (by omega)]
  · have hspec : spectrum ℝ (Matrix.reindex f f A) = spectrum ℝ A := by
      have := AlgEquiv.spectrum_eq (Matrix.reindexAlgEquiv ℝ ℝ f) A
      rwa [Matrix.coe_reindexAlgEquiv] at this
    exact IsGreatest.unique (hspec ▸ lamMax_isGreatest hd _ hB) (lamMax_isGreatest hd A hA)

private theorem starProjection_map_isometry {E E' : Type*} [NormedAddCommGroup E]
    [NormedAddCommGroup E'] [InnerProductSpace ℝ E] [InnerProductSpace ℝ E']
    [FiniteDimensional ℝ E] [FiniteDimensional ℝ E']
    (L : E ≃ₗᵢ[ℝ] E') (K : Submodule ℝ E) (K' : Submodule ℝ E')
    (hK : K' = K.map (L.toLinearEquiv : E →ₗ[ℝ] E')) (x : E) :
    K'.starProjection (L x) = L (K.starProjection x) := by
  subst hK
  rw [Submodule.starProjection_map_apply L K (L x), L.symm_apply_apply]

/-- `overlap` does not change when rows and columns are relabeled, provided `w` is transported
along the column relabeling. The transport is the linear isometry equivalence
`LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f :
EuclideanSpace ℝ (Fin d) ≃ₗᵢ[ℝ] EuclideanSpace ℝ (Fin d')`,
which sends `w` to `w ∘ f.symm` (`LinearIsometryEquiv.piLpCongrLeft_apply`). Stacking uses
this with `finSigmaFinEquiv` on the rows and `Equiv.refl` on the columns. -/
theorem overlap_reindex {n' d' : ℕ} (e : Fin n ≃ Fin n') (f : Fin d ≃ Fin d') :
    overlap (Matrix.reindex e f X) (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f w) =
      overlap X w := by
  have hgram : (Matrix.reindex e f X)ᵀ * (Matrix.reindex e f X)
      = Matrix.reindex f f (Xᵀ * X) := by
    simp [Matrix.reindex_apply, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv]
  have hA : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  have hB : (Matrix.reindex f f (Xᵀ * X)).IsHermitian := by
    rw [← hgram]
    exact isHermitian_transpose_mul_self _
  have hKeq : specSpace (Matrix.reindex f f (Xᵀ * X)) {lamMax (Matrix.reindex f f (Xᵀ * X)) hB}
      = (specSpace (Xᵀ * X) {lamMax (Xᵀ * X) hA}).map
        ((LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f).toLinearEquiv :
          EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d')) := by
    rw [lamMax_reindex f (Xᵀ * X) hA hB]
    exact specSpace_reindex f (Xᵀ * X) _
  have h2 : topProj (Matrix.reindex f f (Xᵀ * X)) hB
        (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f w)
      = LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f (topProj (Xᵀ * X) hA w) :=
    starProjection_map_isometry _ _ _ hKeq w
  change ‖topProj ((Matrix.reindex e f X)ᵀ * Matrix.reindex e f X)
      (isHermitian_transpose_mul_self (Matrix.reindex e f X))
      (LinearIsometryEquiv.piLpCongrLeft 2 ℝ ℝ f w)‖ ^ 2 = _
  rw [topProj_congr hgram _ hB, h2, LinearIsometryEquiv.norm_map]
  rfl

end Spectral

/-! ### Positive scalar multiples of a matrix

`heteroLaw_const_of_singleTableLaw` (`StackSVDWeighted.lean`) needs these: a constant
weighting `w = (t, ..., t)` multiplies the stacked matrix by `t` and its Gram matrix by
`t²`. A positive factor scales every eigenvalue and leaves every eigenspace, hence the top
projector, the overlap and simplicity, alone.

`lamMax_smul` goes through `spectrum`, because `eigenvalues₀` is a sorted tuple built by
choice and Mathlib v4.33.0 has no lemma for the eigenvalues of `t • A`. `lamMax` is the
greatest point of the real spectrum (`lamMax_isGreatest`), the spectrum of `t • A` is
`t • spectrum A` for a unit `t` (`spectrum.smul_mem_smul_iff`), and the greatest point of a
set is unique.
-/

section Smul

variable {d : ℕ}

/-- A real scalar multiple of a Hermitian matrix is Hermitian. -/
theorem isHermitian_smul {A : Matrix (Fin d) (Fin d) ℝ} (hA : A.IsHermitian) (t : ℝ) :
    (t • A).IsHermitian :=
  hA.smul (IsSelfAdjoint.all t)

private theorem toOp_smul (t : ℝ) (A : Matrix (Fin d) (Fin d) ℝ) :
    toOp (t • A) = t • toOp A :=
  map_smul Matrix.toEuclideanLin t A

private theorem smul_cancel {t : ℝ} (ht : t ≠ 0) {x y : EuclideanSpace ℝ (Fin d)}
    (h : t • x = t • y) : x = y := by
  have h2 := congrArg (fun z : EuclideanSpace ℝ (Fin d) => t⁻¹ • z) h
  simpa [inv_smul_smul₀ ht] using h2

/-- The eigenspaces move with the scalar: `t • A` has eigenvalue `t * s` exactly where `A`
has eigenvalue `s`. -/
theorem eigenspace_smul {t : ℝ} (ht : t ≠ 0) (A : Matrix (Fin d) (Fin d) ℝ) (s : ℝ) :
    Module.End.eigenspace (toOp (t • A)) (t * s) = Module.End.eigenspace (toOp A) s := by
  ext x
  rw [Module.End.mem_eigenspace_iff, Module.End.mem_eigenspace_iff, toOp_smul,
    LinearMap.smul_apply, mul_smul]
  exact ⟨fun h => smul_cancel ht h, fun h => congrArg (fun z => t • z) h⟩

/-- `lamMax (t • A) = t * lamMax A` for `0 < t`. -/
theorem lamMax_smul {t : ℝ} (ht : 0 < t) (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : (t • A).IsHermitian) : lamMax (t • A) hB = t * lamMax A hA := by
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    rw [lamMax, lamMax, dif_neg (by omega), dif_neg (by omega), mul_zero]
  · have htne : t ≠ 0 := ht.ne'
    have hAeq : ((Units.mk0 t htne : ℝˣ) • A) = t • A := rfl
    have hgreat : IsGreatest (spectrum ℝ (t • A)) (t * lamMax A hA) := by
      obtain ⟨hmem, hub⟩ := lamMax_isGreatest hd A hA
      constructor
      · have h := (spectrum.smul_mem_smul_iff (a := A) (s := lamMax A hA)
          (r := Units.mk0 t htne)).mpr hmem
        have hv : ((Units.mk0 t htne : ℝˣ) • lamMax A hA : ℝ) = t * lamMax A hA := rfl
        rwa [hAeq, hv] at h
      · rintro x hx
        have hv : ((Units.mk0 t htne : ℝˣ) • (t⁻¹ * x) : ℝ) = x := by
          change t * (t⁻¹ * x) = x
          field_simp
        have hx'' : t⁻¹ * x ∈ spectrum ℝ A :=
          (spectrum.smul_mem_smul_iff (a := A) (s := t⁻¹ * x)
            (r := Units.mk0 t htne)).mp (by rw [hAeq, hv]; exact hx)
        have h2 := mul_le_mul_of_nonneg_left (hub hx'') ht.le
        rwa [← mul_assoc, mul_inv_cancel₀ htne, one_mul] at h2
    exact IsGreatest.unique (lamMax_isGreatest hd (t • A) hB) hgreat

/-- The top eigenspace does not see a positive scalar factor. -/
theorem topSpace_smul {t : ℝ} (ht : 0 < t) (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : (t • A).IsHermitian) : topSpace (t • A) hB = topSpace A hA := by
  change specSpace (t • A) {lamMax (t • A) hB} = specSpace A {lamMax A hA}
  rw [specSpace_singleton, specSpace_singleton, lamMax_smul ht A hA hB,
    eigenspace_smul ht.ne' A (lamMax A hA)]

/-- The top projector does not see a positive scalar factor. -/
theorem topProj_smul {t : ℝ} (ht : 0 < t) (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : (t • A).IsHermitian) : topProj (t • A) hB = topProj A hA := by
  change (topSpace (t • A) hB).starProjection = (topSpace A hA).starProjection
  rw [topSpace_smul ht A hA hB]

/-- Simplicity of the top eigenvalue does not see a positive scalar factor. -/
theorem topSimple_smul {t : ℝ} (ht : 0 < t) (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (hB : (t • A).IsHermitian) : TopSimple (t • A) hB ↔ TopSimple A hA := by
  change Module.finrank ℝ (topSpace (t • A) hB) = 1 ↔ Module.finrank ℝ (topSpace A hA) = 1
  rw [topSpace_smul ht A hA hB]

/-! ### The overlap sees only the Gram matrix -/

/-- `overlap` depends on `X` only through `Xᵀ X`. Two matrices with different row counts and
the same Gram matrix have the same overlap with every `w`; this is what makes the zero rows
of a discarded table invisible (`cor.2`). -/
theorem overlap_congr_gram {p q : ℕ} (X : Matrix (Fin p) (Fin d) ℝ)
    (Y : Matrix (Fin q) (Fin d) ℝ) (h : Xᵀ * X = Yᵀ * Y) (w : EuclideanSpace ℝ (Fin d)) :
    overlap X w = overlap Y w := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2
      = ‖topProj (Yᵀ * Y) (isHermitian_transpose_mul_self Y) w‖ ^ 2
  rw [topProj_congr h _ (isHermitian_transpose_mul_self Y)]

/-- `TopSimple` of the Gram matrix depends on `X` only through `Xᵀ X`. -/
theorem topSimple_congr_gram {p q : ℕ} (X : Matrix (Fin p) (Fin d) ℝ)
    (Y : Matrix (Fin q) (Fin d) ℝ) (h : Xᵀ * X = Yᵀ * Y) :
    TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X)
      ↔ TopSimple (Yᵀ * Y) (isHermitian_transpose_mul_self Y) := by
  change Module.finrank ℝ (topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)) = 1
      ↔ Module.finrank ℝ (topSpace (Yᵀ * Y) (isHermitian_transpose_mul_self Y)) = 1
  rw [topSpace_congr h (isHermitian_transpose_mul_self X) (isHermitian_transpose_mul_self Y)]

/-- The Gram matrix of `t • X` is `t² • (Xᵀ X)`. -/
theorem gram_smul {p : ℕ} (t : ℝ) (X : Matrix (Fin p) (Fin d) ℝ) :
    (t • X)ᵀ * (t • X) = (t ^ 2) • (Xᵀ * X) := by
  rw [Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, sq]

/-- `overlap (t • X) w = overlap X w` for `t ≠ 0`: the Gram matrix scales by `t² > 0`, which
neither the top eigenspace nor its projector sees. -/
theorem overlap_smul {p : ℕ} {t : ℝ} (ht : t ≠ 0) (X : Matrix (Fin p) (Fin d) ℝ)
    (w : EuclideanSpace ℝ (Fin d)) : overlap (t • X) w = overlap X w := by
  have ht2 : (0 : ℝ) < t ^ 2 := by positivity
  have hA : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  have hB : ((t ^ 2) • (Xᵀ * X)).IsHermitian := isHermitian_smul hA _
  change ‖topProj ((t • X)ᵀ * (t • X)) (isHermitian_transpose_mul_self (t • X)) w‖ ^ 2
      = ‖topProj (Xᵀ * X) hA w‖ ^ 2
  rw [topProj_congr (gram_smul t X) _ hB, topProj_smul ht2 _ hA hB]

/-- Simplicity of the top eigenvalue of the Gram matrix is invariant under `X ↦ t • X`. -/
theorem topSimple_gram_smul {p : ℕ} {t : ℝ} (ht : t ≠ 0) (X : Matrix (Fin p) (Fin d) ℝ) :
    TopSimple ((t • X)ᵀ * (t • X)) (isHermitian_transpose_mul_self (t • X))
      ↔ TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
  have ht2 : (0 : ℝ) < t ^ 2 := by positivity
  have hA : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  have hB : ((t ^ 2) • (Xᵀ * X)).IsHermitian := isHermitian_smul hA _
  have h1 : TopSimple ((t • X)ᵀ * (t • X)) (isHermitian_transpose_mul_self (t • X))
      ↔ TopSimple ((t ^ 2) • (Xᵀ * X)) hB := by
    change Module.finrank ℝ (topSpace ((t • X)ᵀ * (t • X))
        (isHermitian_transpose_mul_self (t • X))) = 1
      ↔ Module.finrank ℝ (topSpace ((t ^ 2) • (Xᵀ * X)) hB) = 1
    rw [topSpace_congr (gram_smul t X) (isHermitian_transpose_mul_self (t • X)) hB]
  rw [h1, topSimple_smul ht2 _ hA hB]

end Smul

end StackedSVD
