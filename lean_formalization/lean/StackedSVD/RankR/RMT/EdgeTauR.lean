/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.EdgeR

/-!
# Tasks C0 and C2: the scalars and the edge at a general threshold

Tasks C0 and C2 of `notes/archive/rankr_plan_C.md`. The audit
(`notes/archive/audit_rankr_plan_C_2026-09-02.md`, recommended order, stage 1) folds C0 into its
only consumer, so the three scalar lemmas and the two edge theorems share this file.

`RankR/RMT/EdgeR.lean` proves the edge of the block `W₁ = W₀ + Q_sub Q_subᵀ` at the fixed
threshold `bulkEdge c + ε`, with every column of `Q_sub` subcritical. Track C needs the same
bound at a threshold `τ` that a caller picks, with the columns of `Q_sub` only asked to have
their outlier `ρ` below `τ`. That is a weaker hypothesis: a supercritical spike whose outlier
sits below `τ` is allowed.

Contents:

1. **Scalars** (`ScalarsC`, task C0). `rhoSq_lt_rhoSq`: the outlier is strictly increasing in
   `θ` above the threshold `c < θ⁴`. `secular_pos_of_rhoSq_lt`: above the outlier the secular
   value `1 + (θ² + 1) m(z)` is positive, the supercritical twin of
   `MP.secular_pos_of_subcritical`. `exists_pos_lower_bound`: a finite family of positive
   reals has a positive lower bound.
2. **The edge at `τ`** (`RankRStack.tendsto_measure_lamMax_w1R_le_tau`, task C2). The proof of
   `RankRStack.tendsto_measure_lamMax_w1R_le` with `τ` in place of `bulkEdge c + ε`. One step
   changes. There the diagonal limit `L a = 1 + (coreEig (f a) + 1) m(z)` has the single lower
   bound `1 + (√c + 1) m(z) > 0`, which holds because every column is subcritical. Here `L a`
   is positive per index, by `MP.secular_pos_of_subcritical` when `coreEig (f a) ^ 2 ≤ c` and
   by `ScalarsC.secular_pos_of_rhoSq_lt` otherwise, and `exists_pos_lower_bound` turns the
   `t` bounds into the one accuracy `δ = L₀ / (t + 1)`.
3. **The count at `τ`** (`RankRStack.tendsto_measure_eigenvalues₀_le_tau`). The corollary
   through `Frame.eigenvalues₀_le_of_split`, verbatim at `τ`.

**Why `hsub` reads `ρ < τ` in both regimes.** `rhoSq θ c = bulkEdge c` when `θ⁴ ≤ c`
(`Defs.lean:49`), so at a subcritical index `hτ : bulkEdge c < τ` already gives
`rhoSq (√(coreEig (f a))) c < τ`. The hypothesis is therefore no stronger than the
`coreEig (f a) ^ 2 ≤ c` of `EdgeR.lean` on the indices that file covers.

**The edge event.** `ResolventLimitsR.edge` (`RankR/RMT/Forms.lean:501`) is in `∀ ε > 0` form,
not at a free threshold. This file instantiates it at `ε = (τ - bulkEdge c) / 2`, which is
positive and puts `bulkEdge c + ε` strictly below `τ`.

Numeric check: `$SP/agents/C02/check_scalars.py`, seed 20260902. 200000 draws per scalar
identity, 0 violations, min secular value 1.16e-4; the factorization
`ρ(θ₂) - ρ(θ₁) = (x₂ - x₁)(x₁x₂ - c)/(x₁x₂)` with `x = θ²` holds to 2.4e-14 relative.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 1. Task C0: the scalars -/

namespace ScalarsC

variable {c θ θ₁ θ₂ z : ℝ}

/-- **C0.1.** Above the threshold `c < θ⁴` the outlier `ρ²(θ)` is strictly increasing in `θ`.

With `x = θ²` the difference is `(x₂ - x₁)(x₁x₂ - c)/(x₁x₂)`, and `x₁x₂ > x₁² = θ₁⁴ > c`.
Task C5 uses this to separate the outliers of a strictly decreasing spike family. -/
theorem rhoSq_lt_rhoSq (hc : 0 < c) (hθ₁ : 0 < θ₁) (hsup : c < θ₁ ^ 4) (hlt : θ₁ < θ₂) :
    rhoSq θ₁ c < rhoSq θ₂ c := by
  have hθ₂ : 0 < θ₂ := hθ₁.trans hlt
  have hx1 : 0 < θ₁ ^ 2 := by positivity
  have hx2 : 0 < θ₂ ^ 2 := by positivity
  have hx : θ₁ ^ 2 < θ₂ ^ 2 := by nlinarith
  have hsup₂ : c < θ₂ ^ 4 := by nlinarith
  have hprod : c < θ₁ ^ 2 * θ₂ ^ 2 := by nlinarith
  rw [MP.rhoSq_eq hc hθ₁ hsup, MP.rhoSq_eq hc hθ₂ hsup₂, div_lt_div_iff₀ hx1 hx2]
  nlinarith [mul_pos (sub_pos.mpr hx) (sub_pos.mpr hprod)]

/-- **C0.2.** Above the outlier the secular value is positive: `0 < 1 + (θ² + 1) m(z)` for
`z > ρ²(θ)` and `c < θ⁴`.

This is the supercritical twin of `MP.secular_pos_of_subcritical` (`RMT/MP.lean:506`). The
value is exactly `0` at `z = ρ²(θ)` by `MP.m_rhoSq`, and `m c` is strictly increasing on
`Ici (bulkEdge c)` by `MP.m_strictMonoOn`. -/
theorem secular_pos_of_rhoSq_lt (hc : 0 < c) (hθ : 0 < θ) (hsup : c < θ ^ 4)
    (hz : rhoSq θ c < z) : 0 < 1 + (θ ^ 2 + 1) * MP.m c z := by
  have hb : bulkEdge c < rhoSq θ c := MP.bulkEdge_lt_rhoSq hc hθ hsup
  have hmono : MP.m c (rhoSq θ c) < MP.m c z :=
    MP.m_strictMonoOn hc (Set.mem_Ici.mpr hb.le) (Set.mem_Ici.mpr (hb.trans hz).le) hz
  have h1 : (0 : ℝ) < θ ^ 2 + 1 := by positivity
  have hzero : (θ ^ 2 + 1) * MP.m c (rhoSq θ c) = -1 := by
    rw [MP.m_rhoSq hc hθ hsup]
    field_simp
  nlinarith [mul_lt_mul_of_pos_left hmono h1]

/-- **C0.3.** A finite family of positive reals has a positive lower bound. The empty family
takes `L = 1`. Task C2 uses it to turn `t` per-index bounds into one accuracy. -/
theorem exists_pos_lower_bound {t : ℕ} (g : Fin t → ℝ) (hg : ∀ a, 0 < g a) :
    ∃ L : ℝ, 0 < L ∧ ∀ a, L ≤ g a := by
  classical
  rcases Nat.eq_zero_or_pos t with rfl | ht
  · exact ⟨1, one_pos, fun a => a.elim0⟩
  · have hne : (Finset.univ : Finset (Fin t)).Nonempty := ⟨⟨0, ht⟩, Finset.mem_univ _⟩
    obtain ⟨a₀, -, hmin⟩ := Finset.exists_min_image Finset.univ g hne
    exact ⟨g a₀, hg a₀, fun a => hmin a (Finset.mem_univ a)⟩

end ScalarsC

/-! ### 2. Task C2: the edge of `W₁` at a general threshold -/

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- **Gap G1 at a general threshold** (task C2 of `notes/archive/rankr_plan_C.md`). With every
column of `Q_sub` carrying an outlier `ρ²` strictly below `τ`, `lamMax (W₀ + Q_sub Q_subᵀ) ≤ τ` with
probability tending to 1.

This is `RankRStack.tendsto_measure_lamMax_w1R_le` (`RankR/RMT/EdgeR.lean:378`) with the fixed
threshold `bulkEdge c + ε` replaced by `τ`. The hypothesis `hsub` is weaker than the
subcritical one there: `rhoSq θ c = bulkEdge c` below the threshold, so `hτ` implies `hsub` at
a subcritical index, and a supercritical index is allowed whenever its outlier sits below `τ`.
The index map `f` selects the columns; it must be injective, otherwise the limit matrix is not
diagonal.

Numeric validation: `$SP/agents/C02/check_scalars.py`, seed 20260902, cell "C2 per-index
L a > 0" (200000 draws, 0 violations). -/
theorem tendsto_measure_lamMax_w1R_le_tau [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (hc : 0 < c)
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    {t : ℕ} {f : Fin t → Fin r} (hf : Function.Injective f)
    {τ : ℝ} (hτ : bulkEdge c < τ)
    (hsub : ∀ a : Fin t, rhoSq (Real.sqrt (s.coreEig (f a))) c < τ) :
    Tendsto (fun N => μ N
        {ω | lamMax (s.w1R N ω (U N) f) (s.isHermitian_w1R N ω (U N) f) ≤ τ}) atTop (𝓝 1) := by
  classical
  -- the diagonal limit is positive at every index, in both regimes
  set L : Fin t → ℝ := fun a => 1 + (s.coreEig (f a) + 1) * MP.m c τ with hLdef
  have hLpos : ∀ a, 0 < L a := by
    intro a
    have hlam0 : 0 ≤ s.coreEig (f a) := s.coreEig_nonneg (f a)
    have hsq : Real.sqrt (s.coreEig (f a)) ^ 2 = s.coreEig (f a) := Real.sq_sqrt hlam0
    have h4 : Real.sqrt (s.coreEig (f a)) ^ 4 = s.coreEig (f a) ^ 2 :=
      OutliersR.sqrt_pow_four hlam0
    rw [hLdef]
    rcases le_or_gt (s.coreEig (f a) ^ 2) c with hle | hgt
    · have hpos := MP.secular_pos_of_subcritical hc (Real.sqrt_nonneg _)
        (by rw [h4]; exact hle) hτ
      rwa [hsq] at hpos
    · have hθ : 0 < Real.sqrt (s.coreEig (f a)) :=
        Real.sqrt_pos.mpr (OutliersR.lam_pos hc hlam0 hgt)
      have hpos := ScalarsC.secular_pos_of_rhoSq_lt hc hθ (by rw [h4]; exact hgt) (hsub a)
      rwa [hsq] at hpos
  obtain ⟨L₀, hL₀pos, hL₀le⟩ := ScalarsC.exists_pos_lower_bound L hLpos
  -- the single accuracy
  have ht1 : (0 : ℝ) < (t : ℝ) + 1 := by positivity
  set δ : ℝ := L₀ / ((t : ℝ) + 1) with hδdef
  have hδpos : 0 < δ := div_pos hL₀pos ht1
  have hδt : ∀ a, δ * (t : ℝ) ≤ L a := by
    intro a
    refine le_trans ?_ (hL₀le a)
    rw [hδdef, div_mul_eq_mul_div, div_le_iff₀ ht1]
    nlinarith [hL₀pos.le, Nat.cast_nonneg (α := ℝ) t]
  -- the edge event of `ResolventLimitsR.edge`, at half the room between the edge and `τ`
  have hε2 : (0 : ℝ) < (τ - bulkEdge c) / 2 := by linarith
  have hlt2 : bulkEdge c + (τ - bulkEdge c) / 2 < τ := by linarith
  have hedgeC : Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + (τ - bulkEdge c) / 2)).nullMeasurableSet)
      (h.edge ((τ - bulkEdge c) / 2) hε2)
  have hET : ∀ q : Fin t × Fin t, Tendsto (fun N => μ N
      {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) τ
        (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
        - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
            * MP.m c τ else 0)|}) atTop (𝓝 0) := fun q =>
    s.tendstoInProb_cform_qmatR h (f q.1) (f q.2) hτ δ hδpos
  have hzero : Tendsto (fun N => μ N
      ({ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
            ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) τ
            (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
            - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
                * MP.m c τ else 0)|}))) atTop (𝓝 0) :=
    tendsto_measure_zero_union hedgeC (tendsto_measure_zero_iUnion hET)
  have hincl : ∀ N, ({ω | lamMax (s.w1R N ω (U N) f) (s.isHermitian_w1R N ω (U N) f)
        ≤ τ} : Set (Ω N))ᶜ
      ⊆ {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
            ≤ bulkEdge c + (τ - bulkEdge c) / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) τ
            (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
            - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
                * MP.m c τ else 0)|}) := by
    intro N ω hω
    by_contra hbad
    have hedge : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + (τ - bulkEdge c) / 2 := by
      by_contra hxx
      exact hbad (Set.mem_union_left _ hxx)
    have hclose : ∀ q : Fin t × Fin t,
        |R4.cform (s.rankRW0 N ω (U N)) τ
          (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
          - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
              * MP.m c τ else 0)| ≤ δ := by
      intro q
      by_contra hxx
      exact hbad (Set.mem_union_right _ (Set.mem_iUnion.mpr ⟨q, (not_le.mp hxx).le⟩))
    have hzlam : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N)) < τ := by
      linarith
    have hM : ((1 : Matrix (Fin t) (Fin t) ℝ)
        + (s.qsubR N ω (U N) f)ᵀ * R4.resolv (s.rankRW0 N ω (U N)) τ
          * s.qsubR N ω (U N) f).PosSemidef := by
      refine EdgeR.posSemidef_of_close_to_diag ?_ hδpos.le hδt ?_
      · rw [Matrix.transpose_add, Matrix.transpose_one, Matrix.transpose_mul,
          Matrix.transpose_mul, Matrix.transpose_transpose,
          R4.transpose_resolv (s.isHermitian_rankRW0 N ω (U N)), Matrix.mul_assoc]
      · intro a b
        have hentry : ((1 : Matrix (Fin t) (Fin t) ℝ)
              + (s.qsubR N ω (U N) f)ᵀ * R4.resolv (s.rankRW0 N ω (U N)) τ
                * s.qsubR N ω (U N) f) a b
            - (if a = b then L a else 0)
            = R4.cform (s.rankRW0 N ω (U N)) τ
                (fun l => s.qmatR N ω (U N) l (f a)) (fun l => s.qmatR N ω (U N) l (f b))
              - (if f a = f b then (s.coreEig (f b) + 1)
                  * MP.m c τ else 0) := by
          have hcol : ∀ jj : Fin t, (fun k => s.qsubR N ω (U N) f k jj)
              = fun l => s.qmatR N ω (U N) l (f jj) := fun _ => rfl
          rw [Matrix.add_apply, EdgeR.cform_eq_entry, hcol a, hcol b, Matrix.one_apply, hLdef]
          rcases eq_or_ne a b with rfl | hab
          · rw [if_pos rfl, if_pos rfl, if_pos rfl]
            ring
          · rw [if_neg hab, if_neg hab, if_neg (fun hcon => hab (hf hcon))]
            ring
        rw [hentry]
        exact hclose (a, b)
    exact hω (EdgeR.lamMax_add_le_of_posSemidef (s.hd N)
      (s.isHermitian_rankRW0 N ω (U N)) hzlam hM (s.isHermitian_w1R N ω (U N) f))
  exact tendsto_measure_one_of_bad hincl hzero

/-- **The count at a general threshold** (task C2 of `notes/archive/rankr_plan_C.md`). With the `r`
spikes partitioned into `Fin t` whose outliers sit below `τ` and `Fin u` others, every sorted
eigenvalue of the Gram matrix at index `u` or above is at most `τ`, with probability tending
to 1.

This is `RankRStack.tendsto_measure_eigenvalues₀_le` (`RankR/RMT/EdgeR.lean:500`) at `τ`. It
is the U7a input that task C3 consumes for the upper half of the count. -/
theorem tendsto_measure_eigenvalues₀_le_tau [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (hc : 0 < c)
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (hgram : ∀ N ω, s.gram N ω
      = s.rankRW0 N ω (U N) + s.qmatR N ω (U N) * (s.qmatR N ω (U N))ᵀ)
    {t u : ℕ} (e : Fin t ⊕ Fin u ≃ Fin r) {τ : ℝ} (hτ : bulkEdge c < τ)
    (hsub : ∀ a : Fin t, rhoSq (Real.sqrt (s.coreEig (e (Sum.inl a)))) c < τ) :
    Tendsto (fun N => μ N
        {ω | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
          (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ}) atTop (𝓝 1) := by
  classical
  have hf : Function.Injective (fun a : Fin t => e (Sum.inl a)) :=
    e.injective.comp Sum.inl_injective
  have hmain := s.tendsto_measure_lamMax_w1R_le_tau hc h hf hτ hsub
  have hsubset : ∀ N, {ω | lamMax (s.w1R N ω (U N) (fun a : Fin t => e (Sum.inl a)))
        (s.isHermitian_w1R N ω (U N) (fun a : Fin t => e (Sum.inl a))) ≤ τ}
      ⊆ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
          (s.isHermitian_gram N ω).eigenvalues₀ k ≤ τ} := by
    intro N ω hω k hk
    have hSeq : s.gram N ω
        = s.w1R N ω (U N) (fun a : Fin t => e (Sum.inl a))
          + s.qsubR N ω (U N) (fun b : Fin u => e (Sum.inr b))
            * (s.qsubR N ω (U N) (fun b : Fin u => e (Sum.inr b)))ᵀ := by
      rw [hgram N ω, EdgeR.mul_transpose_split (s.qmatR N ω (U N)) e, w1R, qsubR, qsubR,
        add_assoc]
    exact Frame.eigenvalues₀_le_of_split
      (s.isHermitian_w1R N ω (U N) (fun a : Fin t => e (Sum.inl a)))
      (s.isHermitian_gram N ω) hSeq hω k hk
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le hmain tendsto_const_nhds
    (fun N => measure_mono (hsubset N)) (fun N => prob_le_one)

end RankRStack

end StackedSVD
