/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.SimplicityR

/-!
# Task U5(c): the affine twins of the rank-`r` simplicity lemmas

Deliverable (c) of `notes/archive/rankr_plan_A.md` section 3 (task U5). Affine twins of
`SimplicityR.injective_eigenvalues₀_ae_gaussianMatrix` (`RankR/RMT/SimplicityR.lean:86`),
`simpleSpec_ae_gaussianMatrix` (`:188`) and `topGap_ae_gaussianMatrix` (`:211`), in the pattern
of `topSimple_ae_affine` (`RMT/Simplicity.lean:339-360`): the fixed shift `A` and nonzero scale
`t` replace the specialization `A = 0, t = 1`, so the Gram matrix under the almost-sure
statement is `(A + t • Z)ᵀ * (A + t • Z)` rather than `Zᵀ * Z`.

## Route

Direct transcription of the three rank-1 proofs (`SimplicityR.lean:86-223`), with
`gramPolyRight A t` / `gramPolyLeft A t` in place of `gramPolyRight 0 1` / `gramPolyLeft 0 1`
and `gramPolyRight_ne_zero hd hdp A ht` / `gramPolyLeft_ne_zero hp hpd A ht` in place of the
specialized calls with `one_ne_zero`. Because `A` and `t` no longer collapse the affine family
to `Z` itself, the rewrite `zero_add, one_smul` that the rank-1 proofs use after
`eval_gramPolyRight` / `eval_gramPolyLeft` is not needed here: `eval_gramPolyRight A t Z` and
`eval_gramPolyLeft A t Z` already conclude about `(A + t • Z)ᵀ * (A + t • Z)` and
`(A + t • Z) * (A + t • Z)ᵀ` directly (`RMT/Simplicity.lean:237-249`).

## Hypothesis `hp`

`hp : 0 < p` is not needed by `injective_eigenvalues₀_ae_affine`: `gramPolyRight_ne_zero`
needs only `hd` and `hdp`, exactly as in the rank-1 twin
`injective_eigenvalues₀_ae_gaussianMatrix`. The plan kept the unused `hp` in both for
signature uniformity with `simpleSpec_ae_gaussianMatrix`; the linter campaign of 2026-09-07
dropped it from both (a strict generalization, no proof changed). `hp` **is** used in
`simpleSpec_ae_affine`'s `p ≤ d` branch (`gramPolyLeft_ne_zero hp hpd A ht` needs `0 < p`),
and hence in `topGap_ae_affine`, which is built from it; those two keep `hp`.

## `simpleSpec_ae_gaussianMatrix` as an `A = 0, t = 1` corollary: skipped

`SimplicityR.lean:188` already declares `simpleSpec_ae_gaussianMatrix` directly (not as a
corollary of an affine statement) in the same `StackedSVD` namespace this file also opens, so
restating it here under the same name would be a duplicate declaration, not an addition. The
one-liner check itself succeeds (`simpleSpec_ae_affine hp hd 0 (t := 1) one_ne_zero rk hrk`
reduces to the existing statement after `zero_add, one_smul`), but there is nothing new to add
under a distinct name, so this file skips it.

No `sorry`, no `axiom`, no edits to `SimplicityR.lean`, `Simplicity.lean` or any other existing
file.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-- Affine twin of `injective_eigenvalues₀_ae_gaussianMatrix` (`RankR/RMT/SimplicityR.lean:86`):
for `d ≤ p`, almost every Gaussian `Z` gives an affine shift `A + t • Z` whose `d × d` Gram has
pairwise distinct eigenvalues. The unused `hp : 0 < p` of the plan is dropped (see the header
note). -/
theorem injective_eigenvalues₀_ae_affine {p d : ℕ} (hd : 0 < d) (hdp : d ≤ p)
    (A : Matrix (Fin p) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      Function.Injective (isHermitian_transpose_mul_self (A + t • Z)).eigenvalues₀ := by
  filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramPolyRight A t)
      (gramPolyRight_ne_zero hd hdp A ht)] with Z hZ
  rw [eval_gramPolyRight] at hZ
  exact injective_eigenvalues₀_of_separable _ ((charRes_ne_zero_iff_separable hd _).mp hZ)

/-- Affine twin of `simpleSpec_ae_gaussianMatrix` (`RankR/RMT/SimplicityR.lean:188`): for
`rk ≤ min p d`, almost every Gaussian `Z` gives an affine shift `A + t • Z` whose `d × d` Gram's
top `rk` eigenvalues are pairwise distinct from each other and from every other eigenvalue.
Covers both `d ≤ p` (`injective_eigenvalues₀_ae_affine`) and `p ≤ d`
(`simpleSpec_p_of_gram_left`, monotonicity down from `rk = p`), exactly as the rank-1 proof
case-splits. -/
theorem simpleSpec_ae_affine {p d : ℕ} (hp : 0 < p) (hd : 0 < d)
    (A : Matrix (Fin p) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) (rk : ℕ) (hrk : rk ≤ min p d) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      SimpleSpec ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) rk := by
  rcases le_total d p with hdp | hpd
  · filter_upwards [injective_eigenvalues₀_ae_affine hd hdp A ht] with Z hZ
    exact simpleSpec_of_injective _ hZ rk
  · have hrkp : rk ≤ p := hrk.trans (min_le_left p d)
    filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramPolyLeft A t)
        (gramPolyLeft_ne_zero hp hpd A ht)] with Z hZ
    rw [eval_gramPolyLeft] at hZ
    have hinjB : Function.Injective (isHermitian_mul_transpose_self (A + t • Z)).eigenvalues₀ :=
      injective_eigenvalues₀_of_separable _
        ((charRes_ne_zero_iff_separable hp _).mp (left_ne_zero_of_mul hZ))
    exact simpleSpec_mono hrkp
      (simpleSpec_p_of_gram_left hpd (A + t • Z) hinjB (right_ne_zero_of_mul hZ))

/-- Affine twin of `topGap_ae_gaussianMatrix` (`RankR/RMT/SimplicityR.lean:211`). A five-line
corollary of `simpleSpec_ae_affine`, exactly as the rank-1 proof is a corollary of
`simpleSpec_ae_gaussianMatrix`: `TopGap` needs only the "index `< r` differs from every index
`≥ r`" half of `SimpleSpec`, which combined with `eigenvalues₀_antitone` gives the strict
inequality directly. -/
theorem topGap_ae_affine {p d : ℕ} (hp : 0 < p) (hd : 0 < d)
    (A : Matrix (Fin p) (Fin d) ℝ) {t : ℝ} (ht : t ≠ 0) (r : ℕ) (hr : r ≤ min p d) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      TopGap ((A + t • Z)ᵀ * (A + t • Z)) (isHermitian_transpose_mul_self (A + t • Z)) r := by
  filter_upwards [simpleSpec_ae_affine hp hd A ht r hr] with Z hZ
  intro k l hk hl
  have hval : (k : ℕ) < (l : ℕ) := lt_of_lt_of_le hk hl
  have hne : k ≠ l := by
    intro h; rw [h] at hval; exact lt_irrefl _ hval
  have hkl : k ≤ l := Fin.le_def.mpr hval.le
  exact lt_of_le_of_ne
    ((isHermitian_transpose_mul_self (A + t • Z)).eigenvalues₀_antitone hkl) (hZ k l hk hne)

end StackedSVD
