/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.SpecIdx

/-!
# Measurability of the overlap at one eigenvalue index

`Spectral.lean` proves `measurable_overlap` and `measurable_overlap₂` for the **top**
eigenvalue. Both proofs isolate that eigenvalue by power iteration, because
`(λ_i / λ_max) ^ (2 k)` tends to the indicator of `{λ_max}`. That route reads the top index
and does not reach an index `j > 0`. This file proves the two twins at a free sorted index
`j`, which the rank-`r_i` delocalization lemma needs.

Task B2a of `notes/archive/rankr_plan_B.md` section 2, on the route that the audit
`notes/archive/audit_rankr_plan_B_2026-09-01.md` (attack 2) selects.

## Content

| Here | Mirror in `Spectral.lean` |
|---|---|
| `gramEig`, `continuous_gramEig`, `measurable_gramEig` | `gramLamMax`, `measurable_gramLamMax` |
| `measurable_overlapIdx` | `measurable_overlap` |
| `measurable_overlapIdx₂` | `measurable_overlap₂` |

`gramEig X 0 = gramLamMax X` (`gramEig_zero`) and `overlapIdx X 0 w = overlap X w`
(`overlapIdx_zero`, `LinAlg/SpecIdx.lean`), so each name above is the `j = 0` statement of
its mirror.

## Route: the damping polynomial

Write `A = Xᵀ X`, `t = λ_j(A)` and `K = λ_max(A) + 1`, and put `q(x) = 1 - ((x - t) / K) ^ 2`.
Every eigenvalue of `A` lies in `[0, K - 1]`, so `|λ_i - t| ≤ K - 1 < K`. Therefore
`q(λ_i) ∈ (0, 1]`, and `q(λ_i) = 1` exactly when `λ_i = t`. Three steps follow.

1. `q(A) ^ m` acts diagonally in the eigenbasis of `A`, so
   `⟪w, q(A) ^ m w⟫ = ∑ i, q(λ_i) ^ m ⟪b_i, w⟫ ^ 2` (`quadForm_pow_eig`).
2. Term by term the limit is `∑ i, 1{λ_i = t} ⟪b_i, w⟫ ^ 2`, and `normSq_specProj` of
   `Spectral.lean` reads that sum as `‖specProj A {t} w‖ ^ 2`, which is `overlapIdx X j w`
   (`tendsto_dampedPow`).
3. `X ↦ q(A)` is continuous, because `X ↦ λ_j(A)` is 1-Lipschitz by Weyl at the index `j`
   (`continuous_gramEig`). So `overlapIdx` is a pointwise limit of continuous functions, and
   `measurable_of_tendsto_metrizable'` closes it, as in `measurable_overlap`.

The shift `+ 1` in `K` makes `K ≥ 1`, so this route carries no case split on a zero top
eigenvalue; `measurable_overlap` needs one. The only case split is `j < d` against `j ≥ d`,
and it is on fixed naturals, not on `X`: out of range `eigSetIdx` is empty and `overlapIdx`
is the constant `0`.

The other route of the plan, the shifted matrix `S = K ^ 2 • 1 - (A - t • 1) ^ 2`, is correct
but costs Courant-Fischer plus a kernel identity, because `i ↦ K ^ 2 - (λ_i - t) ^ 2` is not
antitone and so `eigenvalues₀_eq_of_charpoly` does not deliver the spectrum of `S` (audit,
attack 2). This file does not use it.

## Private helpers

Eight helpers of `Spectral.lean` are `private` there (`toOp_pow`, `toCLM_eq`,
`continuous_toEuclideanCLM`, `eigenvalues₀_gram_nonneg`, `eigenvalues₀_le_lamMax`,
`gramLamMax_nonneg`, and the two continuity facts they feed). They are restated here rather
than un-privated in place, so that no existing `olean` changes and no other track has to
rebuild `Spectral.lean`. Weyl at a free index is `abs_eigenvalues₀_sub_le_clm`; the public
`abs_eigenvalues₀_sub_le` of `LinAlg/TopProjPerturb.lean` states the same bound in the scoped
`L2` matrix norm, which this file does not open.
-/


open Filter Topology MeasureTheory
open scoped InnerProductSpace Matrix

namespace StackedSVD

section IdxMeasurable

variable {n p : ℕ}

/-! #### Local copies of the private helpers of `Spectral.lean` -/

private theorem toOp_one : toOp (1 : Matrix (Fin p) (Fin p) ℝ) = LinearMap.id := by
  refine LinearMap.ext fun x => PiLp.ext fun i => ?_
  simp

private theorem toOp_sub (M N : Matrix (Fin p) (Fin p) ℝ) :
    toOp (M - N) = toOp M - toOp N :=
  map_sub Matrix.toEuclideanLin M N

private theorem toOp_smul (r : ℝ) (M : Matrix (Fin p) (Fin p) ℝ) :
    toOp (r • M) = r • toOp M :=
  map_smul Matrix.toEuclideanLin r M

private theorem toOp_pow (A : Matrix (Fin p) (Fin p) ℝ) (k : ℕ) :
    (toOp A) ^ k = toOp (A ^ k) := by
  induction k with
  | zero =>
    refine LinearMap.ext fun x => PiLp.ext fun j => ?_
    simp
  | succ k ih =>
    rw [pow_succ, pow_succ, ih]
    refine LinearMap.ext fun x => PiLp.ext fun j => ?_
    simp [Matrix.toLpLin_apply]

private theorem toCLM_eq (C : Matrix (Fin p) (Fin p) ℝ) :
    LinearMap.toContinuousLinearMap (toOp C) = Matrix.toEuclideanCLM (𝕜 := ℝ) C := rfl

private theorem continuous_toEuclideanCLM :
    Continuous (fun C : Matrix (Fin p) (Fin p) ℝ =>
      (Matrix.toEuclideanCLM (𝕜 := ℝ) C :
        EuclideanSpace ℝ (Fin p) →L[ℝ] EuclideanSpace ℝ (Fin p))) :=
  LinearMap.continuous_of_finiteDimensional
    { toFun := fun C => Matrix.toEuclideanCLM (𝕜 := ℝ) C
      map_add' := fun _ _ => map_add _ _ _
      map_smul' := fun _ _ => map_smul _ _ _ }

private theorem eigenvalues₀_gram_nonneg (X : Matrix (Fin n) (Fin p) ℝ)
    (k : Fin (Fintype.card (Fin p))) :
    0 ≤ (isHermitian_transpose_mul_self X).eigenvalues₀ k := by
  have hpsd : (Xᵀ * X).PosSemidef := by
    simpa using Matrix.posSemidef_conjTranspose_mul_self X
  have h : (isHermitian_transpose_mul_self X).eigenvalues₀ k
      = (isHermitian_transpose_mul_self X).eigenvalues
        (Fintype.equivOfCardEq (Fintype.card_fin _) k) := by
    simp [Matrix.IsHermitian.eigenvalues]
  rw [h]
  exact hpsd.eigenvalues_nonneg _

private theorem eigenvalues₀_le_lamMax (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (hd : 0 < p) (i : Fin (Fintype.card (Fin p))) : hA.eigenvalues₀ i ≤ lamMax A hA := by
  rw [lamMax, dif_pos hd]
  exact hA.eigenvalues₀_antitone (by simp [Fin.le_def])

private theorem gramLamMax_nonneg (X : Matrix (Fin n) (Fin p) ℝ) : 0 ≤ gramLamMax X := by
  rcases Nat.eq_zero_or_pos p with hd | hd
  · subst hd
    simp [gramLamMax, lamMax]
  · rw [gramLamMax, lamMax, dif_pos hd]
    exact eigenvalues₀_gram_nonneg X _

/-- Weyl's inequality at a free sorted index, in the norm that `Spectral.lean` uses.
`abs_eigenvalues₀_sub_le` (`LinAlg/TopProjPerturb.lean`) states the same bound with the scoped
`L2` operator norm of a matrix; this form avoids that scoped instance. -/
private theorem abs_eigenvalues₀_sub_le_clm (A B : Matrix (Fin p) (Fin p) ℝ)
    (hA : A.IsHermitian) (hB : B.IsHermitian) (i : Fin (Fintype.card (Fin p))) :
    |hA.eigenvalues₀ i - hB.eigenvalues₀ i| ≤ ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (A - B)‖ := by
  have key := (symmOp hA).abs_eigenvalues_sub_le_opNorm (symmOp hB) finrank_euclideanSpace i
  rw [← map_sub, toCLM_eq] at key
  exact key

/-! #### The eigenvalue at a free index, as a function of the data matrix -/

/-- `λ_j(Xᵀ X)` at the sorted index `j`; junk value `0` out of range. `gramLamMax` is this
function at `j = 0` (`gramEig_zero`). Mirror: `gramLamMax` (`Defs.lean`). -/
noncomputable def gramEig (X : Matrix (Fin n) (Fin p) ℝ) (j : ℕ) : ℝ :=
  if h : j < Fintype.card (Fin p) then
    (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨j, h⟩ else 0

theorem gramEig_of_lt (X : Matrix (Fin n) (Fin p) ℝ) {j : ℕ}
    (h : j < Fintype.card (Fin p)) :
    gramEig X j = (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨j, h⟩ := dif_pos h

theorem gramEig_of_not_lt (X : Matrix (Fin n) (Fin p) ℝ) {j : ℕ}
    (h : ¬ j < Fintype.card (Fin p)) : gramEig X j = 0 := dif_neg h

/-- Index `0` is the top eigenvalue. -/
theorem gramEig_zero (X : Matrix (Fin n) (Fin p) ℝ) : gramEig X 0 = gramLamMax X := by
  by_cases hp : 0 < p
  · have hc : 0 < Fintype.card (Fin p) := by simpa using hp
    rw [gramEig_of_lt X hc, gramLamMax, lamMax, dif_pos hp]
  · have hc : ¬ 0 < Fintype.card (Fin p) := by simpa using hp
    rw [gramEig_of_not_lt X hc, gramLamMax, lamMax, dif_neg hp]

/-- `X ↦ λ_j(Xᵀ X)` is continuous, by Weyl at the index `j`. Mirror: the private
`continuous_gramLamMax` of `Spectral.lean`, which is this proof at `j = 0`. -/
theorem continuous_gramEig (j : ℕ) :
    Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => gramEig X j) := by
  by_cases hj : j < Fintype.card (Fin p)
  · have hgram : Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => Xᵀ * X) :=
      (continuous_id.matrix_transpose).matrix_mul continuous_id
    rw [continuous_iff_continuousAt]
    intro X0
    have hb : ∀ X : Matrix (Fin n) (Fin p) ℝ, |gramEig X j - gramEig X0 j| ≤
        ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖ := by
      intro X
      rw [gramEig_of_lt X hj, gramEig_of_lt X0 hj]
      exact abs_eigenvalues₀_sub_le_clm _ _ (isHermitian_transpose_mul_self X)
        (isHermitian_transpose_mul_self X0) ⟨j, hj⟩
    have hcont : Continuous (fun X : Matrix (Fin n) (Fin p) ℝ =>
        ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖) :=
      (continuous_toEuclideanCLM.comp (hgram.sub continuous_const)).norm
    have h1 : Tendsto (fun X : Matrix (Fin n) (Fin p) ℝ =>
        ‖Matrix.toEuclideanCLM (𝕜 := ℝ) (Xᵀ * X - X0ᵀ * X0)‖) (𝓝 X0) (𝓝 0) := by
      have h := hcont.continuousAt (x := X0)
      simpa [ContinuousAt] using h
    have h2 : Tendsto (fun X : Matrix (Fin n) (Fin p) ℝ => |gramEig X j - gramEig X0 j|)
        (𝓝 X0) (𝓝 0) := squeeze_zero (fun X => abs_nonneg _) hb h1
    have h3 : Tendsto (fun X : Matrix (Fin n) (Fin p) ℝ => gramEig X j - gramEig X0 j)
        (𝓝 X0) (𝓝 0) :=
      tendsto_zero_iff_norm_tendsto_zero.2 (by simpa [Real.norm_eq_abs] using h2)
    rw [ContinuousAt, ← tendsto_sub_nhds_zero_iff]
    exact h3
  · have he : (fun X : Matrix (Fin n) (Fin p) ℝ => gramEig X j) = fun _ => (0 : ℝ) := by
      funext X
      exact gramEig_of_not_lt X hj
    rw [he]
    exact continuous_const

/-- `X ↦ λ_j(Xᵀ X)` is Borel measurable. Mirror: `measurable_gramLamMax` (`Spectral.lean`). -/
theorem measurable_gramEig (j : ℕ) :
    Measurable (fun X : Matrix (Fin n) (Fin p) ℝ => gramEig X j) :=
  (continuous_gramEig j).measurable

private theorem continuous_gramLamMax :
    Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => gramLamMax X) := by
  have he : (fun X : Matrix (Fin n) (Fin p) ℝ => gramLamMax X) = fun X => gramEig X 0 := by
    funext X
    rw [gramEig_zero]
  rw [he]
  exact continuous_gramEig 0

/-! #### A power that acts diagonally in an orthonormal eigenbasis -/

private theorem apply_pow_eig {ι : Type*} [Fintype ι]
    (T : EuclideanSpace ℝ (Fin p) →ₗ[ℝ] EuclideanSpace ℝ (Fin p))
    (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ (Fin p))) (c : ι → ℝ)
    (hc : ∀ i, T (b i) = c i • b i) (m : ℕ) (w : EuclideanSpace ℝ (Fin p)) :
    (T ^ m) w = ∑ i, (c i ^ m * ⟪b i, w⟫_ℝ) • b i := by
  have hbi : ∀ (k : ℕ) (i : ι), (T ^ k) (b i) = (c i ^ k) • b i := by
    intro k
    induction k with
    | zero => intro i; simp
    | succ k ih =>
      intro i
      have hstep : (T ^ (k + 1)) (b i) = (T ^ k) (T (b i)) := by
        rw [pow_succ]
        rfl
      rw [hstep, hc i, map_smul, ih i, smul_smul, pow_succ]
      congr 1
      ring
  conv_lhs => rw [← b.sum_repr' w]
  rw [map_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [map_smul, hbi m i, smul_smul]
  congr 1
  ring

private theorem quadForm_pow_eig {ι : Type*} [Fintype ι]
    (T : EuclideanSpace ℝ (Fin p) →ₗ[ℝ] EuclideanSpace ℝ (Fin p))
    (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ (Fin p))) (c : ι → ℝ)
    (hc : ∀ i, T (b i) = c i • b i) (m : ℕ) (w : EuclideanSpace ℝ (Fin p)) :
    ⟪w, (T ^ m) w⟫_ℝ = ∑ i, c i ^ m * ⟪b i, w⟫_ℝ ^ 2 := by
  rw [apply_pow_eig T b c hc m w, inner_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [real_inner_smul_right, real_inner_comm w (b i)]
  ring

/-! #### The damped power tends to the spectral projector at one eigenvalue -/

/-- The polynomial route of task B2a. With `q(x) = 1 - ((x - t)/K)²` and `|λ_i - t| < K` for
every eigenvalue, `q(λ_i) ∈ (0, 1]` and `q(λ_i) = 1` exactly at `λ_i = t`. So the quadratic
form of `q(A)^m` tends to the squared norm of the spectral projection at `t`. Nothing here
reads the sorted order, so the statement holds at every index. -/
private theorem tendsto_dampedPow (A : Matrix (Fin p) (Fin p) ℝ) (hA : A.IsHermitian)
    (t K : ℝ) (hK : ∀ i, (hA.eigenvalues₀ i - t) ^ 2 < K ^ 2)
    (w : EuclideanSpace ℝ (Fin p)) :
    Tendsto (fun m : ℕ => ⟪w, toOp ((1 - (K ^ 2)⁻¹ •
        (A - t • (1 : Matrix (Fin p) (Fin p) ℝ)) ^ 2) ^ m) w⟫_ℝ)
      atTop (𝓝 (‖specProj A {t} w‖ ^ 2)) := by
  classical
  set b := (symmOp hA).eigenvectorBasis finrank_euclideanSpace with hbdef
  set M : Matrix (Fin p) (Fin p) ℝ :=
    1 - (K ^ 2)⁻¹ • (A - t • (1 : Matrix (Fin p) (Fin p) ℝ)) ^ 2 with hMdef
  have hone : ∀ x : EuclideanSpace ℝ (Fin p), toOp (1 : Matrix (Fin p) (Fin p) ℝ) x = x := by
    intro x
    rw [toOp_one]
    rfl
  have hAeig : ∀ i, toOp A (b i) = hA.eigenvalues₀ i • b i := by
    intro i
    rw [hbdef]
    exact apply_eigvec hA i
  have hshift : ∀ i, toOp (A - t • (1 : Matrix (Fin p) (Fin p) ℝ)) (b i)
      = (hA.eigenvalues₀ i - t) • b i := by
    intro i
    rw [toOp_sub, LinearMap.sub_apply, hAeig i, toOp_smul, LinearMap.smul_apply, hone,
      sub_smul]
  have hshift2 : ∀ i, toOp ((A - t • (1 : Matrix (Fin p) (Fin p) ℝ)) ^ 2) (b i)
      = ((hA.eigenvalues₀ i - t) ^ 2) • b i := by
    intro i
    have hstep : ((toOp (A - t • (1 : Matrix (Fin p) (Fin p) ℝ))) ^ 2) (b i)
        = toOp (A - t • (1 : Matrix (Fin p) (Fin p) ℝ))
            (toOp (A - t • (1 : Matrix (Fin p) (Fin p) ℝ)) (b i)) := by
      rw [pow_two]
      rfl
    rw [← toOp_pow, hstep, hshift i, map_smul, hshift i, smul_smul, ← pow_two]
  have hc : ∀ i, toOp M (b i) = (1 - (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2) • b i := by
    intro i
    rw [hMdef, toOp_sub, LinearMap.sub_apply, hone, toOp_smul, LinearMap.smul_apply, hshift2 i,
      smul_smul, sub_smul, one_smul,
      show ((K ^ 2)⁻¹ * (hA.eigenvalues₀ i - t) ^ 2)
        = (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2 from by ring]
  have hquad : ∀ m : ℕ, ⟪w, toOp (M ^ m) w⟫_ℝ
      = ∑ i, (1 - (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2) ^ m * ⟪b i, w⟫_ℝ ^ 2 := by
    intro m
    rw [← toOp_pow]
    exact quadForm_pow_eig (toOp M) b (fun i => 1 - (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2) hc m w
  have htarget : ‖specProj A {t} w‖ ^ 2
      = ∑ i, (if hA.eigenvalues₀ i = t then ⟪b i, w⟫_ℝ else 0) ^ 2 := by
    rw [normSq_specProj A hA {t} w, ← hbdef]
    refine Finset.sum_congr rfl fun i _ => ?_
    by_cases hi : hA.eigenvalues₀ i = t
    · rw [if_pos (Set.mem_singleton_iff.mpr hi), if_pos hi]
    · rw [if_neg (fun hcon => hi (Set.mem_singleton_iff.mp hcon)), if_neg hi]
  simp only [hquad, htarget]
  refine tendsto_finsetSum _ fun i _ => ?_
  have hKsq : (0 : ℝ) < K ^ 2 := lt_of_le_of_lt (sq_nonneg _) (hK i)
  have hlt1 : (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2 < 1 := (div_lt_one hKsq).mpr (hK i)
  have hge0 : (0 : ℝ) ≤ (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2 := by positivity
  by_cases hi : hA.eigenvalues₀ i = t
  · rw [if_pos hi, hi, sub_self]
    simp
  · rw [if_neg hi]
    have hne : (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2 ≠ 0 := by
      refine div_ne_zero ?_ (ne_of_gt hKsq)
      exact pow_ne_zero 2 (sub_ne_zero.mpr hi)
    have hpos : (0 : ℝ) < (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2 :=
      lt_of_le_of_ne hge0 (Ne.symm hne)
    have hr : Tendsto (fun m : ℕ => (1 - (hA.eigenvalues₀ i - t) ^ 2 / K ^ 2) ^ m)
        atTop (𝓝 0) :=
      tendsto_pow_atTop_nhds_zero_of_lt_one (by linarith) (by linarith)
    simpa using hr.mul_const (⟪b i, w⟫_ℝ ^ 2)

/-! #### The damping polynomial as a continuous function of the data matrix -/

/-- `q(Xᵀ X)` with `t = λ_j(Xᵀ X)` and `K = λ_max(Xᵀ X) + 1`. The shift `+ 1` makes `K` at
least `1`, so no case split on a zero top eigenvalue appears. -/
private noncomputable def polyMat (X : Matrix (Fin n) (Fin p) ℝ) (j : ℕ) :
    Matrix (Fin p) (Fin p) ℝ :=
  1 - ((gramLamMax X + 1) ^ 2)⁻¹ •
    (Xᵀ * X - gramEig X j • (1 : Matrix (Fin p) (Fin p) ℝ)) ^ 2

private theorem continuous_matrixPow {α : Type*} [TopologicalSpace α]
    {g : α → Matrix (Fin p) (Fin p) ℝ} (hg : Continuous g) (m : ℕ) :
    Continuous fun x => (g x) ^ m := by
  induction m with
  | zero => simpa using continuous_const
  | succ m ih => simpa [pow_succ] using ih.matrix_mul hg

private theorem continuous_polyMat (j : ℕ) :
    Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => polyMat X j) := by
  have hgram : Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => Xᵀ * X) :=
    (continuous_id.matrix_transpose).matrix_mul continuous_id
  have hne : ∀ X : Matrix (Fin n) (Fin p) ℝ, ((gramLamMax X + 1) ^ 2) ≠ 0 := by
    intro X
    have h := gramLamMax_nonneg X
    positivity
  have hK : Continuous (fun X : Matrix (Fin n) (Fin p) ℝ => ((gramLamMax X + 1) ^ 2)⁻¹) :=
    ((continuous_gramLamMax.add continuous_const).pow 2).inv₀ hne
  have hB : Continuous (fun X : Matrix (Fin n) (Fin p) ℝ =>
      Xᵀ * X - gramEig X j • (1 : Matrix (Fin p) (Fin p) ℝ)) :=
    hgram.sub ((continuous_gramEig j).smul continuous_const)
  exact continuous_const.sub (hK.smul (continuous_matrixPow hB 2))

/-- The quadratic form of a continuously varying matrix against a continuously varying
vector is continuous. Both measurability twins use it. -/
private theorem continuous_quadForm {α : Type*} [TopologicalSpace α]
    {f : α → Matrix (Fin p) (Fin p) ℝ} (hf : Continuous f)
    {g : α → EuclideanSpace ℝ (Fin p)} (hg : Continuous g) :
    Continuous fun x => ⟪g x, toOp (f x) (g x)⟫_ℝ := by
  have heq : (fun x => ⟪g x, toOp (f x) (g x)⟫_ℝ)
      = fun x => ⟪g x, (Matrix.toEuclideanCLM (𝕜 := ℝ) (f x)) (g x)⟫_ℝ := by
    funext x
    rfl
  rw [heq]
  exact hg.inner ((continuous_toEuclideanCLM.comp hf).clm_apply hg)

private theorem tendsto_polyMat_quadForm (X : Matrix (Fin n) (Fin p) ℝ) {j : ℕ}
    (hj : j < Fintype.card (Fin p)) (w : EuclideanSpace ℝ (Fin p)) :
    Tendsto (fun m : ℕ => ⟪w, toOp ((polyMat X j) ^ m) w⟫_ℝ) atTop
      (𝓝 (overlapIdx X j w)) := by
  have hp : 0 < p := by
    have h := Nat.lt_of_le_of_lt (Nat.zero_le j) hj
    simpa using h
  have ht : gramEig X j = (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨j, hj⟩ :=
    gramEig_of_lt X hj
  have hgoal : overlapIdx X j w = ‖specProj (Xᵀ * X) {gramEig X j} w‖ ^ 2 := by
    rw [overlapIdx, specProjIdx,
      eigSetIdx_eq_singleton (Xᵀ * X) (isHermitian_transpose_mul_self X) hj, ht]
  rw [hgoal, polyMat]
  refine tendsto_dampedPow (Xᵀ * X) (isHermitian_transpose_mul_self X) (gramEig X j)
    (gramLamMax X + 1) ?_ w
  intro i
  have h0 : 0 ≤ (isHermitian_transpose_mul_self X).eigenvalues₀ i :=
    eigenvalues₀_gram_nonneg X i
  have h1 : (isHermitian_transpose_mul_self X).eigenvalues₀ i ≤ gramLamMax X :=
    eigenvalues₀_le_lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) hp i
  have h2 : 0 ≤ gramEig X j := by
    rw [ht]
    exact eigenvalues₀_gram_nonneg X _
  have h3 : gramEig X j ≤ gramLamMax X := by
    rw [ht]
    exact eigenvalues₀_le_lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) hp _
  have f1 : (0 : ℝ) < gramLamMax X + 1
      - (isHermitian_transpose_mul_self X).eigenvalues₀ i + gramEig X j := by linarith
  have f2 : (0 : ℝ) < gramLamMax X + 1
      + (isHermitian_transpose_mul_self X).eigenvalues₀ i - gramEig X j := by linarith
  nlinarith [mul_pos f1 f2]

/-! #### The affine map, and the two measurability twins -/

/-- The affine map `Z ↦ A + t Z` on matrices is measurable. The tree used to carry two copies:
`RankRStack.measurable_affine_matrix` (`RankR/RMT/EdgeGlueR.lean`) and a private
`measurable_affine_matrix` in `RankR/RMT/DelocAffineR.lean`. Second cleanup pass,
2026-09-02. -/
theorem measurable_affine_matrix {q D : ℕ} (A : Matrix (Fin q) (Fin D) ℝ) (t : ℝ) :
    Measurable fun Z : Matrix (Fin q) (Fin D) ℝ => A + t • Z := by
  refine measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun l => ?_
  have h : (fun Z : Matrix (Fin q) (Fin D) ℝ => (A + t • Z) i l)
      = fun Z => A i l + t * Z i l := by
    funext Z
    rfl
  rw [h]
  exact measurable_const.add
    (measurable_const.mul ((measurable_pi_apply l).comp (measurable_pi_apply i)))


/-- `X ↦ overlapIdx X j w` is Borel measurable, at every fixed sorted index `j`.
Mirror: `measurable_overlap` (`Spectral.lean`), which is this statement at `j = 0`. -/
theorem measurable_overlapIdx (j : ℕ) (w : EuclideanSpace ℝ (Fin p)) :
    Measurable (fun X : Matrix (Fin n) (Fin p) ℝ => overlapIdx X j w) := by
  by_cases hj : j < Fintype.card (Fin p)
  · refine measurable_of_tendsto_metrizable' atTop
      (f := fun m X => ⟪w, toOp ((polyMat X j) ^ m) w⟫_ℝ)
      (fun m => (continuous_quadForm
        (continuous_matrixPow (continuous_polyMat j) m) continuous_const).measurable) ?_
    rw [tendsto_pi_nhds]
    intro X
    exact tendsto_polyMat_quadForm X hj w
  · have he : (fun X : Matrix (Fin n) (Fin p) ℝ => overlapIdx X j w) = fun _ => (0 : ℝ) := by
      funext X
      rw [overlapIdx, specProjIdx_apply_of_not_lt _ _ hj, norm_zero]
      norm_num
    rw [he]
    exact measurable_const

set_option maxHeartbeats 1000000 in
-- The product-space continuity step of this proof, and only this one, passes the default
-- limit; the fixed-vector twin above stays inside it.
/-- Joint measurability of `overlapIdx` in the matrix and in the test vector, at a fixed
sorted index `j`. This is the form that the Fubini step of the rank-`r` delocalization lemma
consumes, where the test direction is itself random. Mirror: `measurable_overlap₂`
(`Spectral.lean`), which is this statement at `j = 0`. -/
theorem measurable_overlapIdx₂ (j : ℕ) :
    Measurable (fun q : Matrix (Fin n) (Fin p) ℝ × EuclideanSpace ℝ (Fin p) =>
      overlapIdx q.1 j q.2) := by
  by_cases hj : j < Fintype.card (Fin p)
  · have hfst : Continuous (fun q : Matrix (Fin n) (Fin p) ℝ × EuclideanSpace ℝ (Fin p) =>
        polyMat q.1 j) := (continuous_polyMat j).comp continuous_fst
    have hm : ∀ m : ℕ, Measurable
        (fun q : Matrix (Fin n) (Fin p) ℝ × EuclideanSpace ℝ (Fin p) =>
          ⟪q.2, toOp ((polyMat q.1 j) ^ m) q.2⟫_ℝ) := fun m =>
      (continuous_quadForm (continuous_matrixPow hfst m) continuous_snd).measurable
    refine measurable_of_tendsto_metrizable' atTop hm ?_
    rw [tendsto_pi_nhds]
    rintro ⟨X, w⟩
    exact tendsto_polyMat_quadForm X hj w
  · have he : (fun q : Matrix (Fin n) (Fin p) ℝ × EuclideanSpace ℝ (Fin p) =>
        overlapIdx q.1 j q.2) = fun _ => (0 : ℝ) := by
      funext q
      rw [overlapIdx, specProjIdx_apply_of_not_lt _ _ hj, norm_zero]
      norm_num
    rw [he]
    exact measurable_const

end IdxMeasurable

end StackedSVD
