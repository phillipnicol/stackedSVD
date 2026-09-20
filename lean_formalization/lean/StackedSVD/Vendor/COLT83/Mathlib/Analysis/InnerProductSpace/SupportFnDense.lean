/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Analysis.InnerProductSpace.SupportFn

set_option autoImplicit false

/-!
# The support function as a supremum over a dense sequence

For a nonempty compact set `K` of a separable inner product space there is a sequence
`q : ℕ → K` dense in `K`, and the support function `supportFn K v = sup_{x ∈ K} ⟪x, v⟫` is the
supremum of the sequence `⟪q n, v⟫` (`exists_forall_supportFn_eq_iSup`), hence the limit of the
nondecreasing partial maxima `partialSupInner q n v = max_{i ≤ n} ⟪q i, v⟫`
(`tendsto_partialSupInner`). This is the countable reduction used to pass from finite to compact
index sets in the Borell–TIS inequality (blueprint `lem:sup_countable_dense`).
-/

@[expose] public section

open Filter Topology
open scoped RealInnerProductSpace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] {K : Set E} {R : ℝ}

/-- The partial maxima `max_{i ≤ n} ⟪q i, v⟫` of the linear forms of a sequence `q`. -/
noncomputable def partialSupInner (q : ℕ → E) (n : ℕ) (v : E) : ℝ := ⨆ i : Fin (n + 1), ⟪q i, v⟫

lemma bddAbove_range_inner_seq {q : ℕ → E} (hq : ∀ n, q n ∈ K) (hR : ∀ x ∈ K, ‖x‖ ≤ R)
    (v : E) : BddAbove (Set.range fun n ↦ ⟪q n, v⟫) :=
  ⟨R * ‖v‖, by
    rintro _ ⟨n, rfl⟩
    exact (real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hR _ (hq n)) (norm_nonneg _))⟩

lemma inner_le_partialSupInner (q : ℕ → E) {n : ℕ} (i : Fin (n + 1)) (v : E) :
    ⟪q i, v⟫ ≤ partialSupInner q n v := by
  unfold partialSupInner
  exact le_ciSup (Finite.bddAbove_range fun i : Fin (n + 1) ↦ ⟪q i, v⟫) i

lemma partialSupInner_le {q : ℕ → E} {n : ℕ} {v : E} {c : ℝ} (h : ∀ i : Fin (n + 1), ⟪q i, v⟫ ≤ c) :
    partialSupInner q n v ≤ c := by
  unfold partialSupInner
  exact ciSup_le h

lemma monotone_partialSupInner (q : ℕ → E) (v : E) : Monotone fun n ↦ partialSupInner q n v :=
  monotone_nat_of_le_succ fun n ↦ partialSupInner_le fun i ↦
    (inner_le_partialSupInner q (Fin.castSucc i) v).trans_eq' (by simp)

lemma abs_partialSupInner_le {q : ℕ → E} (hq : ∀ n, q n ∈ K) (hR : ∀ x ∈ K, ‖x‖ ≤ R) (n : ℕ)
    (v : E) : |partialSupInner q n v| ≤ R * ‖v‖ := by
  have hb : ∀ i : Fin (n + 1), |⟪q i, v⟫| ≤ R * ‖v‖ := fun i ↦
    (abs_real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_right (hR _ (hq i)) (norm_nonneg _))
  refine abs_le.2 ⟨?_, partialSupInner_le fun i ↦ (le_abs_self _).trans (hb i)⟩
  calc -(R * ‖v‖) ≤ ⟪q (0 : Fin (n + 1)), v⟫ := neg_le.2 ((neg_le_abs _).trans (hb 0))
    _ ≤ partialSupInner q n v := inner_le_partialSupInner q 0 v

lemma partialSupInner_le_supportFn {q : ℕ → E} (hq : ∀ n, q n ∈ K) (hR : ∀ x ∈ K, ‖x‖ ≤ R)
    (n : ℕ) (v : E) : partialSupInner q n v ≤ supportFn K v :=
  partialSupInner_le fun i ↦ inner_le_supportFn hR (hq i) v

/-- **The support function of a compact set is a supremum over a dense sequence**
(blueprint `lem:sup_countable_dense`): for a nonempty compact `K` in a separable space there is
a sequence `q` of points of `K` with `supportFn K v = ⨆ n, ⟪q n, v⟫` for every `v`. -/
lemma exists_forall_supportFn_eq_iSup [SecondCountableTopology E] (hK : IsCompact K)
    (hne : K.Nonempty) :
    ∃ q : ℕ → E, (∀ n, q n ∈ K) ∧ ∀ v, supportFn K v = ⨆ n, ⟪q n, v⟫ := by
  obtain ⟨R, hR⟩ : ∃ R, ∀ x ∈ K, ‖x‖ ≤ R := by
    obtain ⟨r, hr⟩ := hK.isBounded.subset_closedBall (0 : E)
    exact ⟨r, fun x hx ↦ mem_closedBall_zero_iff.1 (hr hx)⟩
  have : Nonempty K := hne.to_subtype
  obtain ⟨q, hq⟩ := TopologicalSpace.exists_dense_seq K
  refine ⟨fun n ↦ (q n : E), fun n ↦ (q n).2, fun v ↦ le_antisymm ?_ ?_⟩
  · obtain ⟨x, hx, hxv⟩ := exists_supportFn_eq_inner hK hne v
    rw [hxv]
    refine le_of_forall_pos_le_add fun ε hε ↦ ?_
    obtain ⟨n, hn⟩ := Metric.denseRange_iff.1 hq ⟨x, hx⟩ (ε / (1 + ‖v‖)) (by positivity)
    rw [Subtype.dist_eq, dist_eq_norm] at hn
    have h1 : ⟪x, v⟫ - ⟪(q n : E), v⟫ ≤ ε := by
      calc ⟪x, v⟫ - ⟪(q n : E), v⟫ = ⟪x - q n, v⟫ := (inner_sub_left _ _ _).symm
        _ ≤ ‖x - q n‖ * ‖v‖ := real_inner_le_norm _ _
        _ ≤ ε / (1 + ‖v‖) * ‖v‖ := by gcongr
        _ ≤ ε := by
            rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
            nlinarith [norm_nonneg v]
    calc ⟪x, v⟫ ≤ ⟪(q n : E), v⟫ + ε := by linarith
      _ ≤ (⨆ n, ⟪(q n : E), v⟫) + ε := by
          gcongr
          exact le_ciSup (bddAbove_range_inner_seq (fun n ↦ (q n).2) hR v) n
  · exact ciSup_le fun n ↦ inner_le_supportFn hR (q n).2 v

/-- The partial maxima `max_{i ≤ n} ⟪q i, v⟫` converge to `supportFn K v` when `q` is dense in
the compact set `K` (in the sense of `exists_forall_supportFn_eq_iSup`). -/
lemma tendsto_partialSupInner {q : ℕ → E} (hq : ∀ n, q n ∈ K) (hR : ∀ x ∈ K, ‖x‖ ≤ R)
    (hsup : ∀ v, supportFn K v = ⨆ n, ⟪q n, v⟫) (v : E) :
    Tendsto (fun n ↦ partialSupInner q n v) atTop (𝓝 (supportFn K v)) := by
  have hbdd : BddAbove (Set.range fun n ↦ partialSupInner q n v) :=
    ⟨supportFn K v, by
      rintro _ ⟨n, rfl⟩
      exact partialSupInner_le_supportFn hq hR n v⟩
  have heq : (⨆ n, partialSupInner q n v) = supportFn K v := by
    refine le_antisymm (ciSup_le fun n ↦ partialSupInner_le_supportFn hq hR n v) ?_
    rw [hsup v]
    refine ciSup_le fun n ↦ ?_
    calc ⟪q n, v⟫ = ⟪q (Fin.last n : Fin (n + 1)), v⟫ := by simp
      _ ≤ partialSupInner q n v := inner_le_partialSupInner q _ v
      _ ≤ ⨆ n, partialSupInner q n v := le_ciSup hbdd n
  rw [← heq]
  exact tendsto_atTop_ciSup (monotone_partialSupInner q v) hbdd

end
