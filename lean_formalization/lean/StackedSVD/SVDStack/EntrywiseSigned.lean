/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Deterministic

/-!
# `lem:entrywise_conv_eigenvec` in the paper's signed coordinate form

The paper states the lemma as: with the sign of `x̂ = v_max(Ṽ Ṽᵀ)` chosen so that
`⟨x, x̂⟩ ≥ 0`, where `x = v_max(A_β)`, every coordinate satisfies `x̂_m →_p x_m`.
`lem_entrywise_conv_eigenvec` (`SVDStack/Deterministic.lean`) proves the overlap form
`‖P_top(G_N) x‖² → ⟪x, x'⟫²`, which squares the sign away. This file derives the signed
coordinate form from it (external statement audit of 2026-09-02, finding 3, item L2).

Route. On the event `TopSimple (G N ω)`, which has probability tending to one
(`topSimple_whp_of_tendsto`, or every `ω` when `M = 1`), the overlap equals `⟪x̂_N, x⟫²`, so
`⟪x̂_N, x⟫² → 1` and `|⟪x̂_N, x⟫| → 1` in probability. With `s_N = ±1` the sign of
`⟪x̂_N, x⟫`, `‖s_N x̂_N - x‖² = 2 - 2 |⟪x̂_N, x⟫| → 0`, and each coordinate is bounded by the
norm. Nothing here is Gaussian or asymptotic in `d`: the file is model-free, like the section
`EntrywiseEigenvec` that it extends.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-- The sign that makes `s * t` nonnegative: `1` when `0 ≤ t`, `-1` otherwise. -/
noncomputable def signPos (t : ℝ) : ℝ := if 0 ≤ t then 1 else -1

theorem signPos_mul_self (t : ℝ) : signPos t * t = |t| := by
  unfold signPos
  split_ifs with h
  · rw [one_mul, abs_of_nonneg h]
  · rw [neg_one_mul, abs_of_neg (not_le.mp h)]

theorem abs_signPos (t : ℝ) : |signPos t| = 1 := by
  unfold signPos
  split_ifs <;> norm_num

section EntrywiseSigned

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {M : ℕ}
variable {A : Matrix (Fin M) (Fin M) ℝ} {G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ}

/-- `M ≥ 1` when some `M × M` matrix has a simple top eigenvalue. -/
theorem pos_of_topSimple (hA : A.IsHermitian) (hsimple : TopSimple A hA) : 0 < M := by
  rcases Nat.eq_zero_or_pos M with h0 | h0
  · exfalso
    subst h0
    have hle : Module.finrank ℝ (topSpace A hA)
        ≤ Module.finrank ℝ (EuclideanSpace ℝ (Fin 0)) := Submodule.finrank_le _
    rw [finrank_euclideanSpace, Fintype.card_fin] at hle
    rw [TopSimple] at hsimple
    omega
  · exact h0

/-- A `1 × 1` symmetric matrix has a simple top eigenvalue. -/
theorem topSimple_of_card_one (hA : A.IsHermitian) (hM : M = 1) : TopSimple A hA := by
  have hpos : 0 < M := by omega
  have hle : Module.finrank ℝ (topSpace A hA)
      ≤ Module.finrank ℝ (EuclideanSpace ℝ (Fin M)) := Submodule.finrank_le _
  rw [finrank_euclideanSpace, Fintype.card_fin] at hle
  have hne : topSpace A hA ≠ ⊥ := by
    intro hbot
    have hmem := mem_topSpace_vMax hpos A hA
    rw [hbot, Submodule.mem_bot] at hmem
    have h1 := norm_vMax hpos A hA
    rw [hmem, norm_zero] at h1
    exact zero_ne_one h1
  have hnz : Module.finrank ℝ (topSpace A hA) ≠ 0 :=
    fun h0 => hne (Submodule.finrank_eq_zero.mp h0)
  unfold TopSimple
  omega

/-- With probability tending to one the top eigenvalue of `G_N` is simple, when the top
eigenvalue of the limit `A` is simple. For `M ≥ 2` this is `topSimple_whp_of_tendsto` with the
gap `gap_pos_of_topSimple`; for `M = 1` every `ω` qualifies. -/
theorem topSimple_whp_of_topSimple (hA : A.IsHermitian) (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (A i j))
    (hsimple : TopSimple A hA) :
    Tendsto (fun N => μ N {ω | ¬ TopSimple (G N ω) (hsymm N ω)}) atTop (𝓝 0) := by
  by_cases h1 : 1 < Fintype.card (Fin M)
  · exact topSimple_whp_of_tendsto hA hsymm hconv h1 (gap_pos_of_topSimple hA h1 hsimple)
  · have hM : M = 1 := by
      have := pos_of_topSimple hA hsimple
      rw [Fintype.card_fin] at h1
      omega
    have hempty : ∀ N, {ω | ¬ TopSimple (G N ω) (hsymm N ω)} = (∅ : Set (Ω N)) := fun N =>
      Set.subset_empty_iff.mp fun ω hω => hω (topSimple_of_card_one (hsymm N ω) hM)
    simp only [hempty, measure_empty]
    exact tendsto_const_nhds

/-- `lem:entrywise_conv_eigenvec` in the paper's signed coordinate form. Let `x = v_max(A_β)`
and `x̂_N = v_max(G_N)`, and choose the sign `s_N = signPos ⟪x̂_N, x⟫`, so that
`s_N ⟪x̂_N, x⟫ = |⟪x̂_N, x⟫| ≥ 0` (the paper's "without loss of generality `⟨x, x̂⟩ ≥ 0`").
Then every coordinate of `s_N x̂_N` tends in probability to that of `x`. Same hypotheses as
`lem_entrywise_conv_eigenvec`, which supplies the overlap `‖P_top(G_N) x‖² → 1`. -/
theorem lem_entrywise_conv_eigenvec_signed (β : Fin M → ℝ)
    (G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (Abeta β i j))
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β)) (m : Fin M) :
    TendstoInProb μ
      (fun N ω => signPos ⟪vMax (G N ω) (hsymm N ω), vMax (Abeta β) (isHermitian_Abeta β)⟫_ℝ
        * vMax (G N ω) (hsymm N ω) m)
      (vMax (Abeta β) (isHermitian_Abeta β) m) := by
  set x : EuclideanSpace ℝ (Fin M) := vMax (Abeta β) (isHermitian_Abeta β) with hx
  have hM : 0 < M := pos_of_topSimple _ hsimple
  have hxn : ‖x‖ = 1 := norm_vMax hM _ _
  -- the overlap form, with `⟪x, x⟫² = 1`
  have hover := lem_entrywise_conv_eigenvec β G hsymm hconv hsimple x
  have hxx : ⟪x, x⟫_ℝ ^ 2 = 1 := by
    rw [real_inner_self_eq_norm_sq, hxn]
    norm_num
  rw [hxx] at hover
  -- `⟪x̂_N, x⟫² → 1`: the two agree on the event `TopSimple (G N ω)`
  have hwhp := topSimple_whp_of_topSimple (isHermitian_Abeta β) hsymm hconv hsimple
  have hinner2 : TendstoInProb μ (fun N ω => ⟪vMax (G N ω) (hsymm N ω), x⟫_ℝ ^ 2) 1 := by
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hover
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω | ¬ TopSimple (G N ω) (hsymm N ω)}) ?_ hwhp
    intro N ω hω hs
    exact hω (norm_topProj_sq_eq_inner_sq hs (mem_topSpace_vMax hM _ _)
      (norm_vMax hM _ _) x).symm
  -- `|⟪x̂_N, x⟫| → 1`
  have habs : TendstoInProb μ (fun N ω => |⟪vMax (G N ω) (hsymm N ω), x⟫_ℝ|) 1 := by
    have h := hinner2.comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
    rw [Real.sqrt_one] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    exact Real.sqrt_sq_eq_abs _
  -- `‖s_N x̂_N - x‖² = 2 - 2 |⟪x̂_N, x⟫| → 0`
  have hnsq : TendstoInProb μ
      (fun N ω => ‖signPos ⟪vMax (G N ω) (hsymm N ω), x⟫_ℝ • vMax (G N ω) (hsymm N ω) - x‖ ^ 2)
      0 := by
    have h := (TendstoInProb.const μ 2).sub (habs.const_mul 2)
    rw [mul_one, sub_self] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    rw [norm_sub_sq_real, norm_smul, real_inner_smul_left, signPos_mul_self, Real.norm_eq_abs,
      abs_signPos, one_mul, norm_vMax hM, hxn]
    ring
  -- `‖s_N x̂_N - x‖ → 0`
  have hn : TendstoInProb μ
      (fun N ω => ‖signPos ⟪vMax (G N ω) (hsymm N ω), x⟫_ℝ • vMax (G N ω) (hsymm N ω) - x‖)
      0 := by
    have h := hnsq.comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
    rw [Real.sqrt_zero] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    exact Real.sqrt_sq (norm_nonneg _)
  -- a coordinate is bounded by the norm
  refine TendstoInProb.of_le (fun N => ?_) hn
  filter_upwards with ω
  have hcoord := PiLp.norm_apply_le
    (signPos ⟪vMax (G N ω) (hsymm N ω), x⟫_ℝ • vMax (G N ω) (hsymm N ω) - x) m
  simpa [PiLp.sub_apply, PiLp.smul_apply, smul_eq_mul, Real.norm_eq_abs] using hcoord

end EntrywiseSigned

end StackedSVD
