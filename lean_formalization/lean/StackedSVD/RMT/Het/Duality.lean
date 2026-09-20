/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Split
import StackedSVD.RMT.R4
import StackedSVD.RMT.R4C

/-!
# Item H2: duality between the two Gram matrices

Task H2 of `notes/archive/plan_heterolaw_A.md` (section 2.1 "Duality" and the H2 row of section 4).
The column split of `RMT/Het/Split.lean` puts the rank-one update on the `n` side, that is
on `X Xᵀ`, while `overlap` reads the top eigenvector of the `d` side `Xᵀ X`. Three lemmas
cross the gap:

* `lamMax_gram_comm`: `λ_max(Xᵀ X) = λ_max(X Xᵀ)`. Mathlib's `spectrum.nonzero_mul_comm`
  states the same fact for two square matrices of one algebra, so it does not apply to a
  rectangular `X` without padding; the proof here pushes the top eigenvector through `X`
  and uses the Rayleigh bound `R4.dotProduct_mulVec_le_lamMax` in both directions.
* `topProj_transpose_eq`: the `d`-side overlap read off an `n`-side unit top eigenvector.
* `overlap_eq_inv_lam_qform2`: `overlap X v = 1 / (λ · qᵀ G₀(λ)² q)` on the event
  `λ_max(W₀') < λ`, the identity that replaces R5's numerator bracket.

Paper: `main_paper.tex` lines 1409 to 1440. Plan check (seed 20260830): `0.60025` both ways.

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/Duality.lean` exit 0; 0 `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD
namespace Het

variable {p q : ℕ}

/-! ### The other Gram matrix -/

/-- `X Xᵀ` is Hermitian. The mirror of `isHermitian_transpose_mul_self`. -/
theorem isHermitian_mul_transpose_self (X : Matrix (Fin p) (Fin q) ℝ) :
    (X * Xᵀ).IsHermitian := by
  simpa using Matrix.isHermitian_mul_conjTranspose_self X

private theorem posSemidef_mul_transpose_self (X : Matrix (Fin p) (Fin q) ℝ) :
    (X * Xᵀ).PosSemidef := by
  simpa using Matrix.posSemidef_self_mul_conjTranspose X

private theorem posSemidef_transpose_mul_self (X : Matrix (Fin p) (Fin q) ℝ) :
    (Xᵀ * X).PosSemidef := by
  simpa using Matrix.posSemidef_conjTranspose_mul_self X

private theorem lamMax_nonneg_of_posSemidef {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) (hpsd : A.PosSemidef) (hp : 0 < p) : 0 ≤ lamMax A hA := by
  obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hp
  rw [← hj]
  exact hpsd.eigenvalues_nonneg j

/-- `(X y) ⬝ᵥ z = y ⬝ᵥ (Xᵀ z)`, the adjoint identity in `dotProduct` form. -/
private theorem dotProduct_mulVec_transpose (X : Matrix (Fin p) (Fin q) ℝ)
    (y : Fin q → ℝ) (z : Fin p → ℝ) : (X *ᵥ y) ⬝ᵥ z = y ⬝ᵥ (Xᵀ *ᵥ z) := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ => by ring

private theorem lamMax_gram_le (hp : 0 < p) (hq : 0 < q) (X : Matrix (Fin p) (Fin q) ℝ) :
    lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)
      ≤ lamMax (X * Xᵀ) (isHermitian_mul_transpose_self X) := by
  have hAh := isHermitian_transpose_mul_self X
  have hBh := isHermitian_mul_transpose_self X
  obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hAh hq
  set lam := lamMax (Xᵀ * X) hAh with hlamdef
  set y : Fin q → ℝ := ⇑(hAh.eigenvectorBasis j) with hydef
  have hmv : (Xᵀ * X) *ᵥ y = lam • y := by
    rw [← hj]; exact hAh.mulVec_eigenvectorBasis j
  have hyy : y ⬝ᵥ y = 1 := by
    have h := inner_euclidean_eq_dotProduct (hAh.eigenvectorBasis j) (hAh.eigenvectorBasis j)
    rw [real_inner_self_eq_norm_sq, hAh.eigenvectorBasis.orthonormal.1 j] at h
    simpa using h.symm
  set x : Fin p → ℝ := X *ᵥ y with hxdef
  have hxx : x ⬝ᵥ x = lam := by
    have h1 : x ⬝ᵥ x = y ⬝ᵥ ((Xᵀ * X) *ᵥ y) := by
      rw [hxdef, dotProduct_mulVec_transpose, Matrix.mulVec_mulVec]
    rw [h1, hmv, dotProduct_smul, hyy, smul_eq_mul, mul_one]
  have hxq : (X * Xᵀ) *ᵥ x = lam • x := by
    have h2 : (X * Xᵀ) *ᵥ x = X *ᵥ ((Xᵀ * X) *ᵥ y) := by
      rw [hxdef, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
    rw [h2, hmv, Matrix.mulVec_smul, ← hxdef]
  have hle := R4.dotProduct_mulVec_le_lamMax hBh x
  rw [hxq, dotProduct_smul, smul_eq_mul, hxx] at hle
  have hlam0 : 0 ≤ lam :=
    lamMax_nonneg_of_posSemidef hAh (posSemidef_transpose_mul_self X) hq
  have hB0 : 0 ≤ lamMax (X * Xᵀ) hBh :=
    lamMax_nonneg_of_posSemidef hBh (posSemidef_mul_transpose_self X) hp
  rcases eq_or_lt_of_le hlam0 with h0 | h0
  · rw [← h0]; exact hB0
  · nlinarith

/-- **The duality of plan section 2.1.** The two Gram matrices have the same top
eigenvalue. Both dimensions must be positive: `lamMax` of a `0 × 0` matrix is the junk
value `0`, and `Xᵀ X` of a `0 × q` matrix is the zero matrix, whose `lamMax` is `0` as
well, so nothing is lost. -/
theorem lamMax_gram_comm (hp : 0 < p) (hq : 0 < q) (X : Matrix (Fin p) (Fin q) ℝ) :
    gramLamMax X = lamMax (X * Xᵀ) (isHermitian_mul_transpose_self X) := by
  refine le_antisymm (lamMax_gram_le hp hq X) ?_
  have h := lamMax_gram_le hq hp Xᵀ
  have e1 : lamMax ((Xᵀ)ᵀ * Xᵀ) (isHermitian_transpose_mul_self Xᵀ)
      = lamMax (X * Xᵀ) (isHermitian_mul_transpose_self X) :=
    lamMax_congr (by rw [Matrix.transpose_transpose]) _ _
  have e2 : lamMax (Xᵀ * (Xᵀ)ᵀ) (isHermitian_mul_transpose_self Xᵀ)
      = lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) :=
    lamMax_congr (by rw [Matrix.transpose_transpose]) _ _
  rw [e1, e2] at h
  exact h

/-! ### Small helpers -/

/-- `‖x‖²` as a dot product. -/
private theorem dotProduct_self_eq_norm_sq {m : ℕ} (x : EuclideanSpace ℝ (Fin m)) :
    WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = ‖x‖ ^ 2 := by
  have h := real_inner_eq_dotProduct x x
  rw [real_inner_self_eq_norm_sq] at h
  exact h.symm

/-- A vector of unit dot product has norm one. -/
private theorem norm_eq_one_of_dotProduct {m : ℕ} (x : EuclideanSpace ℝ (Fin m))
    (h : WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = 1) : ‖x‖ = 1 := by
  have h1 : ‖x‖ ^ 2 = 1 := by rw [← dotProduct_self_eq_norm_sq, h]
  have h0 : (0 : ℝ) ≤ ‖x‖ := norm_nonneg x
  have hfac : (‖x‖ - 1) * (‖x‖ + 1) = 0 := by nlinarith [h1]
  rcases mul_eq_zero.mp hfac with h' | h'
  · linarith
  · linarith

/-- Membership in `topSpace`, in matrix language. -/
private theorem mem_topSpace_mulVec_iff {m : ℕ} (A : Matrix (Fin m) (Fin m) ℝ)
    (hA : A.IsHermitian) (x : EuclideanSpace ℝ (Fin m)) :
    x ∈ topSpace A hA ↔ A *ᵥ WithLp.ofLp x = lamMax A hA • WithLp.ofLp x := by
  rw [topSpace_eq_eigenspace, R4.mem_eigenspace_iff']

/-! ### The overlap formula on the `n` side -/

/-- **Step 2 of the duality.** The `d`-side overlap read off an `n`-side unit top
eigenvector `û` of `X Xᵀ`: `overlap X y = (ûᵀ X y)² / λ`.

Proof (H2 row of `notes/archive/plan_heterolaw_A.md` section 4). `0 < gramLamMax X` forces both
`0 < p` (else `‖û‖ = 0`) and `0 < q` (else `lamMax` is the junk value `0`).

1. `lamMax_gram_comm` gives `λ := gramLamMax X = λ_max(X Xᵀ)`.
2. `v̂ := λ^{-1/2} • Xᵀ û` is a unit vector: `‖Xᵀ û‖² = ûᵀ (X Xᵀ) û = λ ‖û‖² = λ`.
3. `v̂` is a top eigenvector of `Xᵀ X`: `(Xᵀ X) (Xᵀ û) = Xᵀ ((X Xᵀ) û) = λ Xᵀ û`.
4. `TopSimple (X Xᵀ)` transfers to `TopSimple (Xᵀ X)` without any finrank machinery: the
   top space of `X Xᵀ` is the line `ℝ û`, so `z` in the top space of `Xᵀ X` has
   `X z = c û`, hence `λ z = Xᵀ (X z) = c Xᵀ û = c √λ v̂` and `z ∈ ℝ v̂`. The two lines
   are therefore equal spans, and `overlap_eq_inner_sq` turns `overlap X y` into `⟪v̂, y⟫²`.
5. `⟪v̂, y⟫ = λ^{-1/2} (Xᵀ û) ⬝ᵥ y = λ^{-1/2} û ⬝ᵥ (X y)` by
   `dotProduct_mulVec_transpose`. Square and simplify. -/
theorem topProj_transpose_eq (X : Matrix (Fin p) (Fin q) ℝ)
    (hsimple : TopSimple (X * Xᵀ) (isHermitian_mul_transpose_self X))
    (hlam : 0 < gramLamMax X) (uh : EuclideanSpace ℝ (Fin p))
    (huh : uh ∈ topSpace (X * Xᵀ) (isHermitian_mul_transpose_self X)) (hnorm : ‖uh‖ = 1)
    (y : EuclideanSpace ℝ (Fin q)) :
    overlap X y = (WithLp.ofLp uh ⬝ᵥ (X *ᵥ WithLp.ofLp y)) ^ 2 / gramLamMax X := by
  have hA : (X * Xᵀ).IsHermitian := isHermitian_mul_transpose_self X
  have hB : (Xᵀ * X).IsHermitian := isHermitian_transpose_mul_self X
  -- both dimensions are positive
  have hp : 0 < p := by
    rcases Nat.eq_zero_or_pos p with h | h
    · subst h
      rw [EuclideanSpace.norm_eq] at hnorm
      simp at hnorm
    · exact h
  have hq0 : 0 < q := by
    rcases Nat.eq_zero_or_pos q with h | h
    · subst h
      simp [gramLamMax, lamMax] at hlam
    · exact h
  have hLA : lamMax (X * Xᵀ) hA = gramLamMax X := (lamMax_gram_comm hp hq0 X).symm
  have hLB : lamMax (Xᵀ * X) hB = gramLamMax X := rfl
  obtain ⟨s, hspos, hmul⟩ : ∃ s : ℝ, 0 < s ∧ s * s = gramLamMax X :=
    ⟨Real.sqrt (gramLamMax X), Real.sqrt_pos.mpr hlam, Real.mul_self_sqrt hlam.le⟩
  have hsne : s ≠ 0 := ne_of_gt hspos
  -- step 1: `uh` is a top eigenvector of `X Xᵀ`
  have huheig : (X * Xᵀ) *ᵥ WithLp.ofLp uh = gramLamMax X • WithLp.ofLp uh := by
    have h := (mem_topSpace_mulVec_iff (X * Xᵀ) hA uh).mp huh
    rwa [hLA] at h
  have huu : WithLp.ofLp uh ⬝ᵥ WithLp.ofLp uh = 1 := by
    rw [dotProduct_self_eq_norm_sq, hnorm, one_pow]
  -- step 2: `Xᵀ uh` has squared length `lam`
  have hxx : (Xᵀ *ᵥ WithLp.ofLp uh) ⬝ᵥ (Xᵀ *ᵥ WithLp.ofLp uh) = gramLamMax X := by
    rw [dotProduct_mulVec_transpose, Matrix.transpose_transpose, Matrix.mulVec_mulVec,
      huheig, dotProduct_smul, huu, smul_eq_mul, mul_one]
  obtain ⟨vv, hvvof⟩ : ∃ vv : EuclideanSpace ℝ (Fin q),
      WithLp.ofLp vv = s⁻¹ • (Xᵀ *ᵥ WithLp.ofLp uh) :=
    ⟨WithLp.toLp 2 (s⁻¹ • (Xᵀ *ᵥ WithLp.ofLp uh)), rfl⟩
  have hvvdot : WithLp.ofLp vv ⬝ᵥ WithLp.ofLp vv = 1 := by
    rw [hvvof, smul_dotProduct, dotProduct_smul, hxx, smul_eq_mul, smul_eq_mul, ← hmul]
    field_simp
  have hvnorm : ‖vv‖ = 1 := norm_eq_one_of_dotProduct vv hvvdot
  have hvvne : vv ≠ 0 := by
    intro h
    rw [h, norm_zero] at hvnorm
    exact zero_ne_one hvnorm
  -- step 3: `vv` is a top eigenvector of `Xᵀ X`
  have hvveig : (Xᵀ * X) *ᵥ WithLp.ofLp vv = gramLamMax X • WithLp.ofLp vv := by
    have hxeig : (Xᵀ * X) *ᵥ (Xᵀ *ᵥ WithLp.ofLp uh)
        = gramLamMax X • (Xᵀ *ᵥ WithLp.ofLp uh) := by
      rw [Matrix.mulVec_mulVec, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, huheig,
        Matrix.mulVec_smul]
    rw [hvvof, Matrix.mulVec_smul, hxeig, smul_comm]
  have hvvmem : vv ∈ topSpace (Xᵀ * X) hB := by
    rw [mem_topSpace_mulVec_iff, hLB]
    exact hvveig
  -- step 4: the top eigenvalue of `Xᵀ X` is simple as well
  have huhne : uh ≠ 0 := by
    intro h
    rw [h, norm_zero] at hnorm
    exact zero_ne_one hnorm
  have hspanU : (Submodule.span ℝ {uh} : Submodule ℝ (EuclideanSpace ℝ (Fin p)))
      = topSpace (X * Xᵀ) hA :=
    Submodule.eq_of_le_of_finrank_eq ((Submodule.span_singleton_le_iff_mem _ _).2 huh)
      (by rw [finrank_span_singleton huhne, hsimple])
  have hspanV : topSpace (Xᵀ * X) hB = Submodule.span ℝ {vv} := by
    refine le_antisymm ?_ ((Submodule.span_singleton_le_iff_mem _ _).2 hvvmem)
    intro z hz
    rw [mem_topSpace_mulVec_iff, hLB] at hz
    have hXz : (WithLp.toLp 2 (X *ᵥ WithLp.ofLp z) : EuclideanSpace ℝ (Fin p))
        ∈ topSpace (X * Xᵀ) hA := by
      rw [mem_topSpace_mulVec_iff, hLA]
      change (X * Xᵀ) *ᵥ (X *ᵥ WithLp.ofLp z) = gramLamMax X • (X *ᵥ WithLp.ofLp z)
      rw [Matrix.mulVec_mulVec, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, hz,
        Matrix.mulVec_smul]
    rw [← hspanU, Submodule.mem_span_singleton] at hXz
    obtain ⟨c, hc⟩ := hXz
    have hofc : X *ᵥ WithLp.ofLp z = c • WithLp.ofLp uh := by
      have h := congrArg WithLp.ofLp hc
      simpa using h.symm
    have hkey : gramLamMax X • WithLp.ofLp z = c • (Xᵀ *ᵥ WithLp.ofLp uh) := by
      have h1 : Xᵀ *ᵥ (X *ᵥ WithLp.ofLp z) = (Xᵀ * X) *ᵥ WithLp.ofLp z :=
        Matrix.mulVec_mulVec _ _ _
      rw [hz] at h1
      rw [← h1, hofc, Matrix.mulVec_smul]
    have hxu : Xᵀ *ᵥ WithLp.ofLp uh = s • WithLp.ofLp vv := by
      rw [hvvof, smul_smul, mul_inv_cancel₀ hsne, one_smul]
    rw [hxu, smul_smul] at hkey
    refine Submodule.mem_span_singleton.mpr ⟨(gramLamMax X)⁻¹ * (c * s), ?_⟩
    apply WithLp.ofLp_injective
    have hz2 : WithLp.ofLp z = ((gramLamMax X)⁻¹ * (c * s)) • WithLp.ofLp vv := by
      calc WithLp.ofLp z = (gramLamMax X)⁻¹ • (gramLamMax X • WithLp.ofLp z) := by
            rw [smul_smul, inv_mul_cancel₀ (ne_of_gt hlam), one_smul]
        _ = (gramLamMax X)⁻¹ • ((c * s) • WithLp.ofLp vv) := by rw [hkey]
        _ = ((gramLamMax X)⁻¹ * (c * s)) • WithLp.ofLp vv := by rw [smul_smul]
    simpa using hz2.symm
  have hsimpleB : TopSimple (Xᵀ * X) hB := by
    rw [TopSimple, hspanV, finrank_span_singleton hvvne]
  -- step 5: read the overlap off `vv`
  rw [overlap_eq_inner_sq X y hsimpleB hvvmem hvnorm]
  have hinner : ⟪vv, y⟫_ℝ = s⁻¹ * (WithLp.ofLp uh ⬝ᵥ (X *ᵥ WithLp.ofLp y)) := by
    rw [real_inner_eq_dotProduct, hvvof, smul_dotProduct, smul_eq_mul,
      dotProduct_mulVec_transpose, Matrix.transpose_transpose]
  have hinv2 : (s⁻¹) ^ 2 = (gramLamMax X)⁻¹ := by
    rw [← hmul, pow_two, ← mul_inv]
  rw [hinner, mul_pow, hinv2]
  ring

/-- **The overlap identity of plan section 2.1.** With `X Xᵀ = W + u uᵀ`, `X v = u` and
`λ_max(W) < λ = λ_max(X Xᵀ)`,

```
overlap X v = (ûᵀ u)² / λ = 1 / (λ · uᵀ G(λ)² u),   û ∝ G(λ) u,   G(z) = (W - z)⁻¹.
```

Numeric check (plan section 2.1, seed 20260830): `0.60025` both ways.

Proof (H2 row of `notes/archive/plan_heterolaw_A.md` section 4).

1. `u ≠ 0`, else `X Xᵀ = W` and the gap hypothesis is false. `0 < λ` because
   `0 < u ⬝ᵥ u = v ⬝ᵥ (Xᵀ X) v ≤ λ (v ⬝ᵥ v)` by the Rayleigh bound.
2. `λ` is an eigenvalue of `W + u uᵀ`, so `R4.secular_eq_zero_iff` gives
   `R4.secular W u λ = 0`, that is `uᵀ G(λ) u = -1`, and `R4.topSpace_eq_span` makes the
   top space of `X Xᵀ` the line through `G(λ) u`.
3. `‖G(λ) u‖² = R4.qform2 W λ u > 0`, so `û := ‖G(λ) u‖⁻¹ • G(λ) u` is a unit vector of
   that line.
4. `topProj_transpose_eq` at this `û`, with `X v = u` (hypothesis `hXv`):
   `overlap X v = (û ⬝ᵥ u)² / λ = 1 / (λ · qform2 W λ u)` by step 2.
5. `hp`, `hq` enter only through `lamMax_gram_comm` inside step 4; `hsimple` is not a
   hypothesis because `R4.topSpace_eq_span` proves it from `λ_max(W) < λ`. -/
theorem overlap_eq_inv_lam_qform2 (hp : 0 < p) (hq : 0 < q)
    (X : Matrix (Fin p) (Fin q) ℝ) (W : Matrix (Fin p) (Fin p) ℝ) (hW : W.IsHermitian)
    (u : Fin p → ℝ) (v : EuclideanSpace ℝ (Fin q))
    (hgram : X * Xᵀ = W + Matrix.vecMulVec u u)
    (hXv : X *ᵥ WithLp.ofLp v = u)
    (hlam : lamMax W hW < gramLamMax X) :
    overlap X v = 1 / (gramLamMax X * R4.qform2 W (gramLamMax X) u) := by
  have hA : (X * Xᵀ).IsHermitian := isHermitian_mul_transpose_self X
  have hAW : (W + Matrix.vecMulVec u u).IsHermitian := hgram ▸ hA
  have hLA : lamMax (X * Xᵀ) hA = gramLamMax X := (lamMax_gram_comm hp hq X).symm
  have hLW : lamMax (W + Matrix.vecMulVec u u) hAW = gramLamMax X :=
    (lamMax_congr hgram hA hAW).symm.trans hLA
  -- `u = 0` would make `X Xᵀ = W`, against the gap hypothesis
  have hu : u ≠ 0 := by
    intro h
    have hzero : Matrix.vecMulVec u u = 0 := by
      ext i j
      simp [h, Matrix.vecMulVec_apply]
    have h0 : X * Xᵀ = W := by rw [hgram, hzero, add_zero]
    have h1 := lamMax_congr h0 hA hW
    rw [hLA] at h1
    rw [← h1] at hlam
    exact lt_irrefl _ hlam
  -- the top eigenvalue is positive because `u = X v ≠ 0`
  have hgeq : gramLamMax X = lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X) := rfl
  have hpos : 0 < gramLamMax X := by
    have h1 : (X *ᵥ WithLp.ofLp v) ⬝ᵥ (X *ᵥ WithLp.ofLp v)
        = WithLp.ofLp v ⬝ᵥ ((Xᵀ * X) *ᵥ WithLp.ofLp v) := by
      rw [dotProduct_mulVec_transpose, Matrix.mulVec_mulVec]
    rw [hXv] at h1
    have h2 := R4.dotProduct_mulVec_le_lamMax (isHermitian_transpose_mul_self X)
      (WithLp.ofLp v)
    rw [← hgeq] at h2
    have h3 : 0 < u ⬝ᵥ u := by
      obtain ⟨i, hi⟩ := Function.ne_iff.mp hu
      rw [Pi.zero_apply] at hi
      simp only [dotProduct]
      exact Finset.sum_pos' (fun j _ => mul_self_nonneg _)
        ⟨i, Finset.mem_univ i, mul_self_pos.mpr hi⟩
    have h4 : 0 ≤ WithLp.ofLp v ⬝ᵥ WithLp.ofLp v := by
      simp only [dotProduct]
      exact Finset.sum_nonneg fun j _ => mul_self_nonneg _
    by_contra hcon
    push Not at hcon
    nlinarith [h1, h2, h3, h4, mul_nonneg (neg_nonneg.mpr hcon) h4]
  -- step 1: the secular equation holds at the top eigenvalue
  have hsec : R4.secular W u (gramLamMax X) = 0 := by
    rw [R4.secular_eq_zero_iff hW hlam, R4.spectrum_toOp]
    obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hAW hp
    rw [← hLW, ← hj]
    exact hAW.eigenvalues_mem_spectrum_real j
  obtain ⟨gv, hgv⟩ : ∃ t : Fin p → ℝ, R4.resolv W (gramLamMax X) *ᵥ u = t := ⟨_, rfl⟩
  have hgne : (WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p)) ≠ 0 := by
    rw [← hgv]
    exact R4.resolv_toLp_ne_zero hW hlam hu
  have hspan : topSpace (X * Xᵀ) hA
      = Submodule.span ℝ {(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))} := by
    rw [topSpace_congr hgram hA hAW, R4.topSpace_eq_span hW hu hlam hsec hAW, hgv]
  have hsimple : TopSimple (X * Xᵀ) hA := by
    rw [TopSimple, hspan, finrank_span_singleton hgne]
  -- step 2: the unit top eigenvector is `G(lam) u` normalized
  have hnpos : 0 < ‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖ := norm_pos_iff.mpr hgne
  have hnne : ‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖ ≠ 0 := ne_of_gt hnpos
  have hgvdot : gv ⬝ᵥ gv = ‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖ ^ 2 := by
    simpa using dotProduct_self_eq_norm_sq (WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))
  have hq2 : R4.qform2 W (gramLamMax X) u
      = ‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖ ^ 2 := by
    rw [← hgv]
    exact R4.dotProduct_resolv_sq_eq_norm_sq hW
  have hqf : u ⬝ᵥ gv = -1 := by
    have hs : R4.secular W u (gramLamMax X)
        = 1 + u ⬝ᵥ (R4.resolv W (gramLamMax X) *ᵥ u) := rfl
    rw [hs, hgv] at hsec
    linarith
  obtain ⟨uh, huhof⟩ : ∃ uh : EuclideanSpace ℝ (Fin p),
      WithLp.ofLp uh = ‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖⁻¹ • gv :=
    ⟨WithLp.toLp 2 (‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖⁻¹ • gv), rfl⟩
  have huhdot : WithLp.ofLp uh ⬝ᵥ WithLp.ofLp uh = 1 := by
    rw [huhof, smul_dotProduct, dotProduct_smul, hgvdot, smul_eq_mul, smul_eq_mul]
    field_simp
  have huhnorm : ‖uh‖ = 1 := norm_eq_one_of_dotProduct uh huhdot
  have huhmem : uh ∈ topSpace (X * Xᵀ) hA := by
    rw [hspan, Submodule.mem_span_singleton]
    refine ⟨‖(WithLp.toLp 2 gv : EuclideanSpace ℝ (Fin p))‖⁻¹, ?_⟩
    apply WithLp.ofLp_injective
    simpa using huhof.symm
  -- step 3: the duality, then the secular identity
  have hgz : gramLamMax X ≠ 0 := ne_of_gt hpos
  rw [topProj_transpose_eq X hsimple hpos uh huhmem huhnorm v, hXv, huhof, smul_dotProduct,
    smul_eq_mul, dotProduct_comm, hqf, hq2]
  field_simp

end Het

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-- `X_w v = q`: the column split makes the rank-one update vector the image of `v`. -/
theorem stackW_mulVec_v (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.stackW w).X N ω *ᵥ WithLp.ofLp (m.stack.v N) = m.qHet w N ω := by
  rw [m.stackW_X_eq_vecMulVec_add w N ω, Matrix.add_mulVec, R4.vecMulVec_mulVec,
    m.dotProduct_v_self_het N, one_smul, ← Matrix.mulVec_mulVec, m.EperpHet_mulVec_v N ω,
    Matrix.mulVec_zero, add_zero]

/-- The heteroscedastic overlap identity: `stackPerfW` is `1 / (λ · qᵀ G₀'(λ)² q)` on the
event `λ_max(W₀') < λ_max(X_w X_wᵀ)`. Corollary of `Het.overlap_eq_inv_lam_qform2` at the
split objects of `RMT/Het/Split.lean`. -/
theorem stackPerfW_eq_inv_lam_qform2 (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (hlam : lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω)
      < gramLamMax ((m.stackW w).X N ω)) :
    m.stackPerfW w N ω
      = 1 / (gramLamMax ((m.stackW w).X N ω)
          * R4.qform2 (m.W0het w N ω) (gramLamMax ((m.stackW w).X N ω)) (m.qHet w N ω)) :=
  Het.overlap_eq_inv_lam_qform2 (m.stack.hn N) (m.stack.hd N) ((m.stackW w).X N ω)
    (m.W0het w N ω) (m.isHermitian_W0het w N ω) (m.qHet w N ω) (m.stack.v N)
    (m.gram_eq_het w N ω) (m.stackW_mulVec_v w N ω) hlam

end MultiTableModel
end StackedSVD
