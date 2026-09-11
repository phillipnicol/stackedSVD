/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Simplicity
import StackedSVD.LinAlg.SpecIdx
import StackedSVD.LinAlg.Eigen

/-!
# Task U5(a): top-`r` simplicity of a Gaussian Gram matrix

Deliverable (a) of `notes/archive/rankr_plan_A.md` section 3 (task U5). Mirrors
`StackedSVD.RMT.Simplicity`, whose header (lines 34-41) already remarks that the same witness
and polynomial argument give distinct top `r` eigenvalues for `r ≤ min n d`, and whose section 3
comment (lines 29-32) states plainly: for `n < d` the characteristic polynomial of `Xᵀ X` has a
repeated zero eigenvalue whenever `d - n ≥ 2`.

## Deviation from the plan's literal headline

`notes/archive/rankr_plan_A.md` section 3 states the headline unconditionally in `p, d > 0`. That
statement is **false**: take `p = 1`, `d = 3`. Then `Zᵀ Z` (`3 × 3`) has rank `1` almost surely,
so `0` is an eigenvalue of multiplicity `2` and `eigenvalues₀` cannot be injective. This is
exactly the obstruction `Simplicity.lean`'s own header names (`d - n ≥ 2`). So:

1. `injective_eigenvalues₀_ae_gaussianMatrix` keeps the plan's exact conclusion but adds the
   hypothesis `hdp : d ≤ p`, the regime where the obstruction cannot occur. This is the direct,
   no-edit port the plan describes (the right-case branch of `topSimple_ae_affine`, `charRes` in
   place of `charRes` and `TopSimple`).
2. The two wrappers are stated for `rk ≤ min p d` / `r ≤ min p d`, not "every rk" as the plan's
   prose suggested for the first: at `rk = min p d` this is the general and useful statement
   (matching `Simplicity.lean`'s own "any `r ≤ min n d`" remark), and it covers **both** regimes,
   including `p ≤ d` (ambient dimension larger), which is the one `RankR/RMT/Stack.lean`'s
   `rankRW0 = E⊥ᵀ E⊥` (a `d(N) × d(N)` Gram built from an `ns(N) × d(N)` block) actually needs.

## Route

`simpleSpec_ae_gaussianMatrix` case-splits on `d ≤ p` vs `p ≤ d`, exactly as
`topSimple_ae_affine` does for `TopSimple`.

* `d ≤ p`: `injective_eigenvalues₀_ae_gaussianMatrix` gives full injectivity, which trivially
  gives `SimpleSpec` at any `rk` (`simpleSpec_of_injective`).
* `p ≤ d`: this is genuinely new deterministic linear algebra, not in `Simplicity.lean`.
  `simpleSpec_p_of_gram_left` transfers full injectivity **and** invertibility of the small
  Gram `X * Xᵀ` (`p × p`) to `SimpleSpec (Xᵀ * X) _ p`, using the rank machinery of
  `LinAlg/Eigen.lean` (`eigenvalues₀_pos_of_lt_rank`, `eigenvalues₀_eq_zero_of_rank_le`,
  already used for exactly this purpose by `RankR/Weighted.lean`'s
  `topGap_optWR_of_rankBR`) plus `finrank_eigenspace_transpose_mul_eq`
  (`LinAlg/TopProjPerturb.lean`) to match eigenspace dimensions between the two Gram matrices.
  The rank of `Xᵀ * X` is `p` (`Matrix.rank_transpose_mul_self`, `Matrix.rank_self_mul_transpose`,
  and `X * Xᵀ` invertible), which pins the top `p` sorted eigenvalues as positive and the rest
  as zero; a `Finset.image` cardinality count then upgrades this to full pairwise distinctness
  among the top `p` (not just each being positive), matching the `p` distinct eigenvalues of
  `X * Xᵀ` one for one.

`topGap_ae_gaussianMatrix` is a five-line corollary of `simpleSpec_ae_gaussianMatrix`: `TopGap`
only needs the "index `< rk` differs from every other index" half of `SimpleSpec`, combined with
`eigenvalues₀_antitone`.

No `sorry`, no `axiom`, no edits to `Simplicity.lean` or any other existing file.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. `SimpleSpec` from full injectivity, and monotonicity in `rk` -/

/-- Full injectivity of the sorted spectrum gives `SimpleSpec` at any rank. -/
theorem simpleSpec_of_injective {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues₀) (rk : ℕ) : SimpleSpec A hA rk :=
  fun _ _ _ hkl heq => hkl (hinj heq).symm

/-- `SimpleSpec` at a rank implies `SimpleSpec` at any smaller rank. -/
theorem simpleSpec_mono {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {rk rk' : ℕ} (h : rk' ≤ rk) (hs : SimpleSpec A hA rk) : SimpleSpec A hA rk' :=
  fun k l hk hkl => hs k l (lt_of_lt_of_le hk h) hkl

/-! ### 2. Full injectivity, the `d ≤ p` regime: a direct port of `topSimple_ae_affine` -/

/-- **Headline, with the added hypothesis `hdp`** (see the deviation note above): for `d ≤ p`,
almost every Gaussian `Z` has a `d × d` Gram `Zᵀ Z` with pairwise distinct eigenvalues. Direct
port of the right-case branch of `topSimple_ae_affine` (`RMT/Simplicity.lean:343-348`), with
`injective_eigenvalues₀_of_separable` in place of `topSimple_of_charRes_ne_zero`. -/
theorem injective_eigenvalues₀_ae_gaussianMatrix (p d : ℕ) (hd : 0 < d) (hdp : d ≤ p) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      Function.Injective (isHermitian_transpose_mul_self Z).eigenvalues₀ := by
  filter_upwards [ae_eval_ne_zero_gaussianMatrix
      (gramPolyRight (0 : Matrix (Fin p) (Fin d) ℝ) 1)
      (gramPolyRight_ne_zero hd hdp _ one_ne_zero)] with Z hZ
  rw [eval_gramPolyRight, zero_add, one_smul] at hZ
  exact injective_eigenvalues₀_of_separable _ ((charRes_ne_zero_iff_separable hd _).mp hZ)

/-! ### 3. The deterministic transfer for the `p ≤ d` regime -/

/-- If the small Gram `X * Xᵀ` (`p × p`, `p ≤ d`) has pairwise distinct eigenvalues and is
invertible, the large Gram `Xᵀ * X` (`d × d`) satisfies `SimpleSpec` at `rk = p`: each of its
top `p` eigenvalues differs from every other eigenvalue, including the repeated zero eigenvalue
of the null space. Pure linear algebra; the probabilistic input (`hinjB`, `hdetB`) is supplied
by `gramPolyLeft_ne_zero` at the call site. -/
theorem simpleSpec_p_of_gram_left {p d : ℕ} (hpd : p ≤ d) (X : Matrix (Fin p) (Fin d) ℝ)
    (hinjB : Function.Injective (isHermitian_mul_transpose_self X).eigenvalues₀)
    (hdetB : (X * Xᵀ).det ≠ 0) :
    SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) p := by
  classical
  set hA := isHermitian_transpose_mul_self X
  set hB := isHermitian_mul_transpose_self X
  have hcardp : Fintype.card (Fin p) = p := Fintype.card_fin p
  have hcardd : Fintype.card (Fin d) = d := Fintype.card_fin d
  -- Rank of both Gram matrices: `p`, from invertibility of the small one.
  have hunit : IsUnit (X * Xᵀ) :=
    (Matrix.isUnit_iff_isUnit_det _).mpr (isUnit_iff_ne_zero.mpr hdetB)
  have hrankB : (X * Xᵀ).rank = p := by rw [Matrix.rank_of_isUnit _ hunit, hcardp]
  have hrankA : (Xᵀ * X).rank = p :=
    (Matrix.rank_transpose_mul_self X).trans
      ((Matrix.rank_self_mul_transpose X).symm.trans hrankB)
  -- Both Gram matrices are positive semidefinite.
  have hpsdA : (Xᵀ * X).PosSemidef := by simpa using Matrix.posSemidef_conjTranspose_mul_self X
  have hnnA : ∀ k, 0 ≤ hA.eigenvalues₀ k :=
    fun k => by rw [← eigenvalues_eigIdx hA]; exact hpsdA.eigenvalues_nonneg _
  have hpsdB : (X * Xᵀ).PosSemidef := by simpa using Matrix.posSemidef_self_mul_conjTranspose X
  have hnnB : ∀ k, 0 ≤ hB.eigenvalues₀ k :=
    fun k => by rw [← eigenvalues_eigIdx hB]; exact hpsdB.eigenvalues_nonneg _
  -- The top `p` sorted eigenvalues of the large Gram are positive; the rest are zero.
  have hposA : ∀ k : Fin (Fintype.card (Fin d)), (k : ℕ) < p → 0 < hA.eigenvalues₀ k :=
    fun k hk => eigenvalues₀_pos_of_lt_rank hA hnnA (le_of_eq hrankA.symm) hk
  have hzeroA : ∀ l : Fin (Fintype.card (Fin d)), p ≤ (l : ℕ) → hA.eigenvalues₀ l = 0 :=
    fun l hl => eigenvalues₀_eq_zero_of_rank_le hA hnnA (le_of_eq hrankA) hl
  -- Every sorted eigenvalue of the small Gram is positive (it has full rank `p`).
  have hposB : ∀ i : Fin (Fintype.card (Fin p)), 0 < hB.eigenvalues₀ i :=
    fun i => eigenvalues₀_pos_of_lt_rank hB hnnB (le_of_eq (hcardp.trans hrankB.symm)) i.isLt
  -- The indices below `p` form a Finset of size `p`.
  set S1 : Finset (Fin (Fintype.card (Fin d))) :=
    Finset.univ.filter (fun j => (j : ℕ) < p) with hS1def
  have hS1card : S1.card = p := card_filter_lt (by rw [hcardd]; exact hpd)
  -- Every eigenvalue of the small Gram is attained, inside the prefix `S1`, by the large Gram:
  -- `finrank_eigenspace_transpose_mul_eq` matches eigenspace dimensions at nonzero values.
  have himg : (Finset.univ : Finset (Fin (Fintype.card (Fin p)))).image hB.eigenvalues₀
      ⊆ S1.image hA.eigenvalues₀ := by
    intro t ht
    rw [Finset.mem_image] at ht
    obtain ⟨i, -, hti⟩ := ht
    have htne : t ≠ 0 := hti ▸ (hposB i).ne'
    have hdim : 0 < Module.finrank ℝ (Module.End.eigenspace (toOp (Xᵀ * X)) t) := by
      rw [finrank_eigenspace_transpose_mul_eq X htne]
      exact lt_of_lt_of_le (one_le_finrank_eigenspace hB i hti) le_rfl
    rw [finrank_eigenspace_eq_card hA t, Finset.card_pos] at hdim
    obtain ⟨j, hj⟩ := hdim
    rw [Finset.mem_filter] at hj
    have hjlt : (j : ℕ) < p := by
      by_contra hcon
      exact htne (hj.2 ▸ hzeroA j (not_lt.mp hcon))
    rw [Finset.mem_image]
    refine ⟨j, ?_, hj.2⟩
    rw [hS1def, Finset.mem_filter]
    exact ⟨Finset.mem_univ _, hjlt⟩
  -- Cardinality: the image of `S1` under the large Gram's spectrum has exactly `p` elements,
  -- the same as `S1` itself, so the large Gram is injective on `S1`.
  have hcardB : (Finset.univ.image hB.eigenvalues₀ : Finset ℝ).card = p := by
    simp only [Finset.card_image_of_injOn hinjB.injOn, Finset.card_univ, hcardp]
  have hcardeq : (S1.image hA.eigenvalues₀).card = S1.card := by
    have h1 : p ≤ (S1.image hA.eigenvalues₀).card :=
      calc p = (Finset.univ.image hB.eigenvalues₀ : Finset ℝ).card := hcardB.symm
        _ ≤ (S1.image hA.eigenvalues₀).card := Finset.card_le_card himg
    have h2 : (S1.image hA.eigenvalues₀).card ≤ S1.card := Finset.card_image_le
    omega
  have hInjOn : Set.InjOn hA.eigenvalues₀ (S1 : Set (Fin (Fintype.card (Fin d)))) :=
    Finset.injOn_of_card_image_eq hcardeq
  -- Assemble `SimpleSpec`: two cases on whether the other index is also in the prefix.
  intro k l hk hkl heq
  by_cases hlp : (l : ℕ) < p
  · have hkS1 : k ∈ (S1 : Set (Fin (Fintype.card (Fin d)))) := by
      rw [Finset.mem_coe, hS1def, Finset.mem_filter]; exact ⟨Finset.mem_univ _, hk⟩
    have hlS1 : l ∈ (S1 : Set (Fin (Fintype.card (Fin d)))) := by
      rw [Finset.mem_coe, hS1def, Finset.mem_filter]; exact ⟨Finset.mem_univ _, hlp⟩
    exact hkl (hInjOn hlS1 hkS1 heq).symm
  · rw [hzeroA l (not_lt.mp hlp)] at heq
    exact (hposA k hk).ne' heq.symm

/-! ### 4. `SimpleSpec` and `TopGap` for `gaussianMatrix`, both regimes -/

/-- The general wrapper: for `rk ≤ min p d`, almost every Gaussian `Z` has a `d × d` Gram whose
top `rk` eigenvalues are pairwise distinct from each other and from every other eigenvalue.
Covers both `d ≤ p` (full injectivity, `simpleSpec_of_injective`) and `p ≤ d`
(`simpleSpec_p_of_gram_left`, monotonicity down from `rk = p`). -/
theorem simpleSpec_ae_gaussianMatrix (p d : ℕ) (hp : 0 < p) (hd : 0 < d) (rk : ℕ)
    (hrk : rk ≤ min p d) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      SimpleSpec (Zᵀ * Z) (isHermitian_transpose_mul_self Z) rk := by
  rcases le_total d p with hdp | hpd
  · filter_upwards [injective_eigenvalues₀_ae_gaussianMatrix p d hd hdp] with Z hZ
    exact simpleSpec_of_injective _ hZ rk
  · have hrkp : rk ≤ p := hrk.trans (min_le_left p d)
    filter_upwards [ae_eval_ne_zero_gaussianMatrix
        (gramPolyLeft (0 : Matrix (Fin p) (Fin d) ℝ) 1)
        (gramPolyLeft_ne_zero hp hpd _ one_ne_zero)] with Z hZ
    rw [eval_gramPolyLeft, zero_add, one_smul] at hZ
    have hinjB : Function.Injective (isHermitian_mul_transpose_self Z).eigenvalues₀ :=
      injective_eigenvalues₀_of_separable _
        ((charRes_ne_zero_iff_separable hp _).mp (left_ne_zero_of_mul hZ))
    exact simpleSpec_mono hrkp (simpleSpec_p_of_gram_left hpd Z hinjB (right_ne_zero_of_mul hZ))

/-- **U5(a), the `TopGap` wrapper.** For `r ≤ min p d`, almost every Gaussian `Z` has a `d × d`
Gram with a top-`r` eigengap: every eigenvalue at index `≥ r` is strictly below every eigenvalue
at index `< r`. A corollary of `simpleSpec_ae_gaussianMatrix`: `TopGap` needs only the "index
`< r` differs from every index `≥ r`" half of `SimpleSpec`, which combined with
`eigenvalues₀_antitone` gives the strict inequality directly (no separate count or gap
argument, matching the task note's suggested route). -/
theorem topGap_ae_gaussianMatrix (p d : ℕ) (hp : 0 < p) (hd : 0 < d) (r : ℕ)
    (hr : r ≤ min p d) :
    ∀ᵐ Z ∂(gaussianMatrix p d),
      TopGap (Zᵀ * Z) (isHermitian_transpose_mul_self Z) r := by
  filter_upwards [simpleSpec_ae_gaussianMatrix p d hp hd r hr] with Z hZ
  intro k l hk hl
  have hval : (k : ℕ) < (l : ℕ) := lt_of_lt_of_le hk hl
  have hne : k ≠ l := by
    intro h; rw [h] at hval; exact lt_irrefl _ hval
  have hkl : k ≤ l := Fin.le_def.mpr hval.le
  exact lt_of_le_of_ne
    ((isHermitian_transpose_mul_self Z).eigenvalues₀_antitone hkl) (hZ k l hk hne)

end StackedSVD
