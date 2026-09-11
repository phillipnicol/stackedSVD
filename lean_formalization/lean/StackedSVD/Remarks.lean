/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVDWeighted
import StackedSVD.StackSVD.Main
import StackedSVD.StackSVD.Weighted
import StackedSVD.SVDStack.Weighted

/-!
# Remark facades: the closed-form instances of the paper's two comparison remarks

`notes/archive/remark_facades.md`. The paper states two remarks with no display and no proof, and
then prints a numeric instance of each:

* `remark:stack_outperform_svd` (`main_paper.tex:566`): unweighted stacksvd can beat
  optimally weighted svdstack. Instance: `θ_i = c_i = 1` for every table. Every table sits at
  its own detection threshold, so `β = 0` and svdstack has performance `0`, while the stack of
  `M` tables clears the joint threshold and reaches `1 - 2/(M+1)`.
* `remark:svd_outperform_stack` (`main_paper.tex:588`): unweighted svdstack can beat
  binary-weighted stacksvd. Three instances: two equal supercritical tables among `M` null
  tables; `M = 3` with `c = (1, 1, c₃)` and `θ = (2, 2, c₃^{1/4})`; and `M = 2` with
  `θ = (√5, 4)`, `c = (1, 38.4)`.

Each theorem here is a parameter choice plugged into a proved closed form, so no new random
matrix theory enters. The hypotheses are those of the theorem instantiated (`hc`, `hreg`,
`hG`). The scalar work is `svdstackLimit_pair`, the twin of `Scalars.svdstackLimit_const` for
a `β` with two equal nonzero coordinates.

## What the third instance shows

At `θ₃ = c₃^{1/4}` the third table sits exactly at its own threshold, so the `cor.2`
detectable set is `{0, 1}` and binary-weighted stacksvd is the stack of tables `0` and `1`,
with limit `62/72 = 0.8611`. That is above svdstack's `6/7 = 0.8571`, so the paper's example
does not separate the two methods in the direction its heading claims; the paper's printed
formula is the *unweighted* stack of all three tables. `remark_svd_outperform_stack_three_binary`
records all three facts, and `paper_edits.md` E6 proposes `θ₃ = (c₃+1)^{1/4}` as the fix.

## Not stated here

The uniform bound over all weight vectors at `β = 0` (instance (i)(c) of the note) needs
`StackedSVD.SVDStack.Rayleigh`, which this file does not import. Only the unweighted svdstack
half of instance (i) is stated.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### The scalar limit for two equal nonzero `β` and the rest zero

`A_β = β βᵀ + diag(1 - β_i²)` at `β = β₀ (e_{i₀} + e_{i₁})` has top eigenvalue `1 + β₀²` with
eigenvector `(e_{i₀} + e_{i₁})/√2`. The five private lemmas below repeat, for a general `β`,
the pattern of `Scalars.lean` for the equicorrelated case; they are private because
`Scalars.lean` keeps its own copies and this file must not edit it.
-/

section PairScalar

variable {M : ℕ}

/-- One row of `A_β x` for a general `β`. -/
private theorem Abeta_row (β : Fin M → ℝ) (x : EuclideanSpace ℝ (Fin M)) (j : Fin M) :
    toOp (Abeta β) x j = β j * (∑ k, β k * x k) + (1 - β j ^ 2) * x j := by
  rw [Scalars.toOp_apply]
  have hterm : ∀ k : Fin M, Abeta β j k * x k
      = β j * (β k * x k) + (if j = k then (1 - β j ^ 2) * x k else 0) := by
    intro k
    have hA : Abeta β j k = β j * β k + (if j = k then 1 - β j ^ 2 else 0) := by
      simp [Abeta, Matrix.vecMulVec_apply, Matrix.diagonal_apply]
    rw [hA]
    split_ifs with h <;> ring
  rw [Finset.sum_congr rfl fun k _ => hterm k, Finset.sum_add_distrib, ← Finset.mul_sum,
    Finset.sum_ite_eq _ j (fun k => (1 - β j ^ 2) * x k), if_pos (Finset.mem_univ j)]

/-- The quadratic form of `A_β` for a general `β`. -/
private theorem quadForm_Abeta (β : Fin M → ℝ) (x : EuclideanSpace ℝ (Fin M)) :
    ⟪toOp (Abeta β) x, x⟫_ℝ = (∑ k, β k * x k) ^ 2 + ∑ i, (1 - β i ^ 2) * x i ^ 2 := by
  rw [real_inner_eq_dotProduct]
  have hd : (WithLp.ofLp (toOp (Abeta β) x)) ⬝ᵥ WithLp.ofLp x
      = ∑ i, (β i * (∑ k, β k * x k) + (1 - β i ^ 2) * x i) * x i :=
    Finset.sum_congr rfl fun i _ => by rw [← Abeta_row β x i]
  rw [hd]
  have hsplit : ∀ i : Fin M, (β i * (∑ k, β k * x k) + (1 - β i ^ 2) * x i) * x i
      = (∑ k, β k * x k) * (β i * x i) + (1 - β i ^ 2) * x i ^ 2 := fun i => by ring
  rw [Finset.sum_congr rfl fun i _ => hsplit i, Finset.sum_add_distrib, ← Finset.mul_sum]
  ring

/-- A sum over `Fin M` collapses to two terms when the summand vanishes off `{i₀, i₁}`. -/
private theorem sum_eq_pair {i₀ i₁ : Fin M} (hne : i₀ ≠ i₁) (f : Fin M → ℝ)
    (h0 : ∀ j, j ≠ i₀ → j ≠ i₁ → f j = 0) : ∑ j, f j = f i₀ + f i₁ := by
  classical
  rw [← Finset.sum_pair hne]
  refine (Finset.sum_subset (Finset.subset_univ _) ?_).symm
  intro j _ hj
  simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hj
  exact h0 j hj.1 hj.2

/-- `λ_max(A_β) = 1 + β₀²` when exactly two coordinates of `β` equal `β₀` and the rest are
zero. The quadratic form at a unit `x` is `β₀²(x_{i₀}+x_{i₁})² + 1 - β₀²(x_{i₀}²+x_{i₁}²)`, at
most `1 + β₀²` because `2ab ≤ a² + b² ≤ 1`, and equal to it at `(e_{i₀}+e_{i₁})/√2`. -/
private theorem lamMax_Abeta_pair {i₀ i₁ : Fin M} (hne : i₀ ≠ i₁) (β₀ : ℝ) (β : Fin M → ℝ)
    (hβ : ∀ i, β i = if i = i₀ ∨ i = i₁ then β₀ else 0) :
    lamMax (Abeta β) (isHermitian_Abeta β) = 1 + β₀ ^ 2 := by
  classical
  have hM : 0 < M := lt_of_le_of_lt (Nat.zero_le i₀.val) i₀.isLt
  have hb0 : β i₀ = β₀ := by rw [hβ]; simp
  have hb1 : β i₁ = β₀ := by rw [hβ]; simp
  have hbz : ∀ j, j ≠ i₀ → j ≠ i₁ → β j = 0 := by
    intro j h0 h1
    rw [hβ]
    simp [h0, h1]
  have hquad : ∀ x : EuclideanSpace ℝ (Fin M),
      ⟪toOp (Abeta β) x, x⟫_ℝ
        = (β₀ * (x i₀ + x i₁)) ^ 2 + ((∑ i, x i ^ 2) - β₀ ^ 2 * (x i₀ ^ 2 + x i₁ ^ 2)) := by
    intro x
    rw [quadForm_Abeta]
    have hS1 : ∑ k, β k * x k = β₀ * (x i₀ + x i₁) := by
      rw [sum_eq_pair hne (fun k => β k * x k)
        (fun j h0 h1 => by rw [hbz j h0 h1]; ring), hb0, hb1]
      ring
    have hS2 : ∑ i, (1 - β i ^ 2) * x i ^ 2
        = (∑ i, x i ^ 2) - β₀ ^ 2 * (x i₀ ^ 2 + x i₁ ^ 2) := by
      have h1 : ∀ i : Fin M, (1 - β i ^ 2) * x i ^ 2 = x i ^ 2 - β i ^ 2 * x i ^ 2 :=
        fun i => by ring
      rw [Finset.sum_congr rfl fun i _ => h1 i, Finset.sum_sub_distrib,
        sum_eq_pair hne (fun i => β i ^ 2 * x i ^ 2)
          (fun j h0 h1 => by rw [hbz j h0 h1]; ring), hb0, hb1]
      ring
    rw [hS1, hS2]
  refine le_antisymm (Scalars.lamMax_le_of_quadForm hM _ _ fun x hx => ?_) ?_
  · have hxsq : ∑ i, x i ^ 2 = 1 := by
      have h := Scalars.euclid_norm_sq x
      rw [hx, one_pow] at h
      exact h.symm
    have hpairle : x i₀ ^ 2 + x i₁ ^ 2 ≤ ∑ i, x i ^ 2 := by
      have h := Finset.sum_le_sum_of_subset_of_nonneg
        (Finset.subset_univ ({i₀, i₁} : Finset (Fin M))) (fun i _ _ => sq_nonneg (x i))
      rwa [Finset.sum_pair hne] at h
    rw [hquad x, hxsq]
    rw [hxsq] at hpairle
    nlinarith [mul_nonneg (sq_nonneg β₀) (sq_nonneg (x i₀ - x i₁)),
      mul_nonneg (sq_nonneg β₀) (sub_nonneg.mpr hpairle)]
  · have h2 : (0 : ℝ) < 2 := by norm_num
    have hr2 : ((Real.sqrt 2)⁻¹ : ℝ) ^ 2 = 1 / 2 := by
      rw [inv_pow, Real.sq_sqrt h2.le]
      norm_num
    set u : EuclideanSpace ℝ (Fin M) :=
      WithLp.toLp 2 (fun i : Fin M => if i = i₀ then (Real.sqrt 2)⁻¹ else
        if i = i₁ then (Real.sqrt 2)⁻¹ else 0) with hudef
    have hu0 : u i₀ = (Real.sqrt 2)⁻¹ := by simp [hudef]
    have hu1 : u i₁ = (Real.sqrt 2)⁻¹ := by simp [hudef, hne.symm]
    have huz : ∀ j, j ≠ i₀ → j ≠ i₁ → u j = 0 := by
      intro j h0 h1
      simp [hudef, h0, h1]
    have husq : ∑ i, u i ^ 2 = 1 := by
      rw [sum_eq_pair hne (fun i => u i ^ 2) (fun j h0 h1 => by rw [huz j h0 h1]; ring),
        hu0, hu1, hr2]
      norm_num
    have hunorm : ‖u‖ = 1 := by
      have h := Scalars.euclid_norm_sq u
      rw [husq] at h
      rw [← Real.sqrt_sq (norm_nonneg u), h, Real.sqrt_one]
    have hle := Scalars.le_lamMax_of_unit _ (isHermitian_Abeta β) u hunorm
    rw [hquad u, husq, hu0, hu1] at hle
    have hval : ((β₀ * ((Real.sqrt 2)⁻¹ + (Real.sqrt 2)⁻¹)) ^ 2
        + (1 - β₀ ^ 2 * ((Real.sqrt 2)⁻¹ ^ 2 + (Real.sqrt 2)⁻¹ ^ 2))) = 1 + β₀ ^ 2 := by
      have hexp : (β₀ * ((Real.sqrt 2)⁻¹ + (Real.sqrt 2)⁻¹)) ^ 2
          = 4 * β₀ ^ 2 * ((Real.sqrt 2)⁻¹ ^ 2) := by ring
      rw [hexp, hr2]
      ring
    rw [hval] at hle
    exact hle

/-- `(βᵀ v_max(A_β))² = 2β₀²` for the two-coordinate `β`. The eigenvalue equation kills every
coordinate of the top eigenvector off `{i₀, i₁}` and makes the two kept coordinates equal, so
`v_{i₀}² = 1/2` and `βᵀ v_max = 2β₀ v_{i₀}`. -/
private theorem inner_beta_vMax_pair {i₀ i₁ : Fin M} (hne : i₀ ≠ i₁) {β₀ : ℝ} (hβ₀ : β₀ ≠ 0)
    (β : Fin M → ℝ) (hβ : ∀ i, β i = if i = i₀ ∨ i = i₁ then β₀ else 0) :
    (β ⬝ᵥ WithLp.ofLp (vMax (Abeta β) (isHermitian_Abeta β))) ^ 2 = 2 * β₀ ^ 2 := by
  classical
  have hM : 0 < M := lt_of_le_of_lt (Nat.zero_le i₀.val) i₀.isLt
  have hb0 : β i₀ = β₀ := by rw [hβ]; simp
  have hb1 : β i₁ = β₀ := by rw [hβ]; simp
  have hbz : ∀ j, j ≠ i₀ → j ≠ i₁ → β j = 0 := by
    intro j h0 h1
    rw [hβ]
    simp [h0, h1]
  set v := vMax (Abeta β) (isHermitian_Abeta β) with hvdef
  have hmem : v ∈ topSpace (Abeta β) (isHermitian_Abeta β) := mem_topSpace_vMax hM _ _
  have hnv : ‖v‖ = 1 := norm_vMax hM _ _
  have heig := toOp_of_mem_topSpace hmem
  have hlam := lamMax_Abeta_pair hne β₀ β hβ
  have hrow : ∀ j : Fin M,
      β j * (∑ k, β k * v k) + (1 - β j ^ 2) * v j = (1 + β₀ ^ 2) * v j := by
    intro j
    rw [← Abeta_row β v j, heig, hlam]
    simp
  have hz : ∀ j, j ≠ i₀ → j ≠ i₁ → v j = 0 := by
    intro j h0 h1
    have h := hrow j
    rw [hbz j h0 h1] at h
    have h2 : β₀ ^ 2 * v j = 0 := by linear_combination -h
    rcases mul_eq_zero.mp h2 with h3 | h3
    · exact absurd h3 (pow_ne_zero 2 hβ₀)
    · exact h3
  have hD0 : ∑ k, β k * v k = 2 * β₀ * v i₀ := by
    have h := hrow i₀
    rw [hb0] at h
    have h2 : β₀ * ((∑ k, β k * v k) - 2 * β₀ * v i₀) = 0 := by linear_combination h
    rcases mul_eq_zero.mp h2 with h3 | h3
    · exact absurd h3 hβ₀
    · linarith
  have hD1 : ∑ k, β k * v k = 2 * β₀ * v i₁ := by
    have h := hrow i₁
    rw [hb1] at h
    have h2 : β₀ * ((∑ k, β k * v k) - 2 * β₀ * v i₁) = 0 := by linear_combination h
    rcases mul_eq_zero.mp h2 with h3 | h3
    · exact absurd h3 hβ₀
    · linarith
  have hveq : v i₀ = v i₁ := by
    have h4 : β₀ * (v i₀ - v i₁) = 0 := by linear_combination (hD1 - hD0) / 2
    rcases mul_eq_zero.mp h4 with h5 | h5
    · exact absurd h5 hβ₀
    · linarith
  have hnorm1 : ∑ j, v j ^ 2 = 1 := by
    have h := Scalars.euclid_norm_sq v
    rw [hnv, one_pow] at h
    exact h.symm
  have hpair : v i₀ ^ 2 + v i₁ ^ 2 = 1 := by
    rw [← sum_eq_pair hne (fun j => v j ^ 2) (fun j h0 h1 => by rw [hz j h0 h1]; ring)]
    exact hnorm1
  have hv0sq : v i₀ ^ 2 = 1 / 2 := by
    rw [← hveq] at hpair
    linarith
  have hdot : (β ⬝ᵥ WithLp.ofLp v) = ∑ k, β k * v k := rfl
  rw [hdot, hD0]
  rw [show (2 * β₀ * v i₀) ^ 2 = 4 * β₀ ^ 2 * v i₀ ^ 2 by ring, hv0sq]
  ring

/-- `svdstackLimit` at `β = β₀ (e_{i₀} + e_{i₁})`: the top eigenvalue of `A_β` is `1 + β₀²`
with eigenvector `(e_{i₀} + e_{i₁})/√2`, so the limit is `2β₀²/(1 + β₀²)`. Twin of
`Scalars.svdstackLimit_const` for two equal nonzero coordinates. `β₀ < 1` is not a
hypothesis: the proof does not read it, and every application takes `β₀ = beta θ c`, which
satisfies it anyway (`beta_mem_Ico`). -/
theorem svdstackLimit_pair (i₀ i₁ : Fin M) (hne : i₀ ≠ i₁) (β₀ : ℝ) (hβ₀ : 0 < β₀)
    (β : Fin M → ℝ) (hβ : ∀ i, β i = if i = i₀ ∨ i = i₁ then β₀ else 0) :
    svdstackLimit β = 2 * β₀ ^ 2 / (1 + β₀ ^ 2) := by
  unfold svdstackLimit
  rw [inner_beta_vMax_pair hne hβ₀.ne' β hβ, lamMax_Abeta_pair hne β₀ β hβ]

/-- The scalar limit of the paper's third instance: the unweighted stack of the three tables
loses all its performance as `c₃` grows. The proof squeezes between `0` and `78/√c₃`. -/
theorem remark_three_stack_tendsto_zero :
    Tendsto (fun c₃ : ℝ => (16 * Real.sqrt c₃ + 62) / (c₃ + 17 * Real.sqrt c₃ + 72)) atTop
      (nhds 0) := by
  refine squeeze_zero' ?_ ?_ (Real.tendsto_sqrt_atTop.const_div_atTop 78)
  · filter_upwards [eventually_ge_atTop (0 : ℝ)] with c hc
    have h1 : 0 ≤ Real.sqrt c := Real.sqrt_nonneg c
    positivity
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with c hc
    have hs1 : 1 ≤ Real.sqrt c := by
      rw [show (1 : ℝ) = Real.sqrt 1 by simp]
      exact Real.sqrt_le_sqrt hc
    have hsq : Real.sqrt c ^ 2 = c := Real.sq_sqrt (by linarith)
    have hspos : (0 : ℝ) < Real.sqrt c := by linarith
    rw [div_le_div_iff₀ (by nlinarith) hspos]
    nlinarith [hs1, hsq]

end PairScalar

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [∀ N, IsProbabilityMeasure (μ N)]

/-! ### `remark:stack_outperform_svd`, instance (i): `θ_i = c_i = 1` -/

/-- **`remark:stack_outperform_svd`** (`main_paper.tex:566`), stacksvd half: at `θ_i = c_i = 1`
and `M ≥ 2` the unweighted stack reaches `1 - 2/(M+1)`. The stack of `M` tables has signal
`‖θ‖₂² = M` and aspect ratio `‖c‖₁ = M`, so it clears the joint threshold `M > 1`. Instance of
`thm_simple_thm1_stacksvd_gaussian` at `θ₀ = c₀ = 1`. -/
theorem remark_stack_outperform_svd_stack [NeZero M] (m : MultiTableModel μ M n d)
    (hM : 2 ≤ M) (hθ : ∀ i, (m.tbl i).θ = 1) (hreg : ∀ i, (m.tbl i).Regime 1)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (1 - 2 / ((M : ℝ) + 1)) := by
  have h := m.thm_simple_thm1_stacksvd_gaussian 1 1 one_pos hθ hreg hG
  have hMR : (2 : ℝ) ≤ (M : ℝ) := by exact_mod_cast hM
  rw [if_pos (by nlinarith : (1 : ℝ) < (M : ℝ) * 1 ^ 4)] at h
  have heq : 1 - ((1 : ℝ) + 1 ^ 2) / ((M : ℝ) * 1 ^ 4 + 1 ^ 2) = 1 - 2 / ((M : ℝ) + 1) := by
    norm_num
  rwa [heq] at h

/-- **`remark:stack_outperform_svd`** (`main_paper.tex:566`), svdstack half: at `θ_i = c_i = 1`
every table sits at its own detection threshold `θ_i⁴ = c_i`, so `β = 0` and unweighted
svdstack has performance `0`, whatever `M` is. Instance of
`thm_svd_stack_general_zero_gaussian`. The paper says the same of *optimally weighted*
svdstack; that uniform statement needs `StackedSVD.SVDStack.Rayleigh`, which this file does
not import. -/
theorem remark_stack_outperform_svd_svdstack [NeZero M] (m : MultiTableModel μ M n d)
    (hθ : ∀ i, (m.tbl i).θ = 1) (hreg : ∀ i, (m.tbl i).Regime 1)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 :=
  m.thm_svd_stack_general_zero_gaussian (fun _ => 1) (fun _ => 0) (fun _ => one_pos)
    (fun i => by rw [hθ i, beta_eq_zero_of_not_thr (by norm_num : ¬ (1 : ℝ) < 1 ^ 4)])
    (fun _ => rfl) hreg hG

/-! ### `remark:svd_outperform_stack`, instance (ii): two equal supercritical tables -/

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:588`), stacksvd half: with two tables at
`θ₀` and `M - 2` null tables, all at aspect ratio `c₀`, the unweighted stack has performance
`0` as soon as `4θ₀⁴ ≤ M c₀`. The stack's signal is `‖θ‖₂² = 2θ₀²` and its aspect ratio is
`‖c‖₁ = M c₀`, so the guard of `stackSVDLimit` asks for `4θ₀⁴ > M c₀`. The paper's
`M > 4θ₀⁴/c₀` is the strict form; equality also gives `0`. -/
theorem remark_svd_outperform_stack_pair_stack [NeZero M] (m : MultiTableModel μ M n d)
    (i₀ i₁ : Fin M) (hne : i₀ ≠ i₁) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = if i = i₀ ∨ i = i₁ then θ₀ else 0)
    (hM : 4 * θ₀ ^ 4 ≤ (M : ℝ) * c₀)
    (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N)) 0 := by
  have hθ0 : (m.tbl i₀).θ = θ₀ := by rw [hθ]; simp
  have hθ1 : (m.tbl i₁).θ = θ₀ := by rw [hθ]; simp
  have hθz : ∀ j, j ≠ i₀ → j ≠ i₁ → (m.tbl j).θ = 0 := by
    intro j h0 h1
    rw [hθ]
    simp [h0, h1]
  have h := m.prop_stacksvd_general_gaussian (fun _ => c₀) (fun _ => hc) hreg hG
  have hsum : ∑ i, (m.tbl i).θ ^ 2 = 2 * θ₀ ^ 2 := by
    rw [sum_eq_pair hne (fun i => (m.tbl i).θ ^ 2)
      (fun j h0 h1 => by rw [hθz j h0 h1]; ring), hθ0, hθ1]
    ring
  have hcsum : ∑ _i : Fin M, c₀ = (M : ℝ) * c₀ := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hlim : stackSVDLimit (fun i => (m.tbl i).θ) (fun _ => c₀) = 0 := by
    unfold stackSVDLimit
    rw [if_neg]
    rw [hsum, hcsum, gt_iff_lt, not_lt]
    nlinarith
  rwa [hlim] at h

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:588`), svdstack half: with the same two
tables at `θ₀` above their own threshold `c₀ < θ₀⁴` and `M - 2` null tables, unweighted
svdstack reaches `2β₀²/(1 + β₀²)` with `β₀ = beta θ₀ c₀`, whatever `M` is. Instance of
`thm_svd_stack_general_gaussian`, evaluated by `svdstackLimit_pair`. -/
theorem remark_svd_outperform_stack_pair_svdstack [NeZero M] (m : MultiTableModel μ M n d)
    (i₀ i₁ : Fin M) (hne : i₀ ≠ i₁) (θ₀ c₀ : ℝ) (hc : 0 < c₀) (hthr : c₀ < θ₀ ^ 4)
    (hθ : ∀ i, (m.tbl i).θ = if i = i₀ ∨ i = i₁ then θ₀ else 0)
    (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω)
      (2 * beta θ₀ c₀ ^ 2 / (1 + beta θ₀ c₀ ^ 2)) := by
  have hβpos : 0 < beta θ₀ c₀ := beta_pos_of_thr hc hthr
  have hzero : beta 0 c₀ = 0 := beta_eq_zero_of_not_thr (by simpa using hc.le)
  have hβdef : ∀ i, (if i = i₀ ∨ i = i₁ then beta θ₀ c₀ else 0)
      = beta (m.tbl i).θ c₀ := by
    intro i
    rw [hθ i]
    by_cases h : i = i₀ ∨ i = i₁
    · rw [if_pos h, if_pos h]
    · rw [if_neg h, if_neg h, hzero]
  have hthr2 : ∃ i j : Fin M, i ≠ j
      ∧ 0 < (if i = i₀ ∨ i = i₁ then beta θ₀ c₀ else 0)
      ∧ 0 < (if j = i₀ ∨ j = i₁ then beta θ₀ c₀ else 0) := by
    refine ⟨i₀, i₁, hne, ?_, ?_⟩
    · rw [if_pos (Or.inl rfl)]; exact hβpos
    · rw [if_pos (Or.inr rfl)]; exact hβpos
  have h := m.thm_svd_stack_general_gaussian (fun _ => c₀)
    (fun i => if i = i₀ ∨ i = i₁ then beta θ₀ c₀ else 0) (fun _ => hc) hβdef hthr2 hreg hG
  rwa [svdstackLimit_pair i₀ i₁ hne (beta θ₀ c₀) hβpos _ (fun i => rfl)] at h

/-! ### `remark:svd_outperform_stack`, instance (iii): `M = 3`, `c = (1, 1, c₃)` -/

section Three

variable {n₃ : Fin 3 → ℕ → ℕ}

/-- The fourth power of `c₃^{1/4}` written with two square roots. -/
private theorem sqrt_sqrt_pow_four {c₃ : ℝ} (hc₃ : 0 ≤ c₃) :
    Real.sqrt (Real.sqrt c₃) ^ 4 = c₃ := by
  rw [show (4 : ℕ) = 2 * 2 from rfl, pow_mul, Real.sq_sqrt (Real.sqrt_nonneg c₃),
    Real.sq_sqrt hc₃]

/-- `beta 2 1 ^ 2 = 3/4`, the `β²` of the two informative tables of instance (iii). -/
private theorem beta_two_one_sq : beta 2 1 ^ 2 = 3 / 4 := by
  rw [Scalars.beta_sq, betaSq, if_pos (by norm_num : (2 : ℝ) ^ 4 > 1)]
  norm_num

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:600`), svdstack half of the paper's
`M = 3` example: `c = (1, 1, c₃)` and `θ = (2, 2, c₃^{1/4})` give `β₁² = β₂² = 3/4` and
`β₃ = 0`, so unweighted svdstack reaches `2(3/4)/(1 + 3/4) = 6/7`, whatever `c₃` is. -/
theorem remark_svd_outperform_stack_three_svdstack (m : MultiTableModel μ 3 n₃ d) (c₃ : ℝ)
    (hc₃ : 0 < c₃) (hθ : ∀ i, (m.tbl i).θ = ![2, 2, Real.sqrt (Real.sqrt c₃)] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 1, c₃] i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (6 / 7) := by
  have hq : Real.sqrt (Real.sqrt c₃) ^ 4 = c₃ := sqrt_sqrt_pow_four hc₃.le
  have hβpos : 0 < beta 2 1 := beta_pos_of_thr one_pos (by norm_num)
  have hb2 : beta (Real.sqrt (Real.sqrt c₃)) c₃ = 0 :=
    beta_eq_zero_of_not_thr (by rw [hq]; exact lt_irrefl c₃)
  have hcpos : ∀ i : Fin 3, 0 < ![(1 : ℝ), 1, c₃] i := by
    intro i
    fin_cases i <;> simp [hc₃]
  have hβdef : ∀ i : Fin 3, (if i = 0 ∨ i = 1 then beta 2 1 else 0)
      = beta (m.tbl i).θ (![(1 : ℝ), 1, c₃] i) := by
    intro i
    fin_cases i <;> simp [hθ, hb2]
  have hthr2 : ∃ i j : Fin 3, i ≠ j
      ∧ 0 < (if i = 0 ∨ i = 1 then beta 2 1 else 0)
      ∧ 0 < (if j = 0 ∨ j = 1 then beta 2 1 else 0) := by
    refine ⟨0, 1, by decide, ?_, ?_⟩
    · rw [if_pos (Or.inl rfl)]; exact hβpos
    · rw [if_pos (Or.inr rfl)]; exact hβpos
  have h := m.thm_svd_stack_general_gaussian ![(1 : ℝ), 1, c₃]
    (fun i => if i = 0 ∨ i = 1 then beta 2 1 else 0) hcpos hβdef hthr2 hreg hG
  rw [svdstackLimit_pair 0 1 (by decide) (beta 2 1) hβpos _ (fun i => rfl),
    beta_two_one_sq, show (2 : ℝ) * (3 / 4) / (1 + 3 / 4) = 6 / 7 by norm_num] at h
  exact h

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:604`), stacksvd half of the paper's
`M = 3` example: the unweighted stack of all three tables has `‖θ‖₂² = 8 + √c₃` and
`‖c‖₁ = 2 + c₃`, so its limit is the paper's `(16√c₃ + 62)/(c₃ + 17√c₃ + 72)`. That value goes
to `0` as `c₃` grows (`remark_three_stack_tendsto_zero`). -/
theorem remark_svd_outperform_stack_three_stack (m : MultiTableModel μ 3 n₃ d) (c₃ : ℝ)
    (hc₃ : 0 < c₃) (hθ : ∀ i, (m.tbl i).θ = ![2, 2, Real.sqrt (Real.sqrt c₃)] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 1, c₃] i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      ((16 * Real.sqrt c₃ + 62) / (c₃ + 17 * Real.sqrt c₃ + 72)) := by
  have hs : Real.sqrt c₃ ^ 2 = c₃ := Real.sq_sqrt hc₃.le
  have hsnn : 0 ≤ Real.sqrt c₃ := Real.sqrt_nonneg c₃
  have hcpos : ∀ i : Fin 3, 0 < ![(1 : ℝ), 1, c₃] i := by
    intro i
    fin_cases i <;> simp [hc₃]
  have h := m.prop_stacksvd_general_gaussian ![(1 : ℝ), 1, c₃] hcpos hreg hG
  have hsum : ∑ i, (m.tbl i).θ ^ 2 = 8 + Real.sqrt c₃ := by
    rw [Fin.sum_univ_three, hθ 0, hθ 1, hθ 2]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
      Matrix.cons_val_two, Matrix.tail_cons]
    rw [Real.sq_sqrt hsnn]
    ring
  have hcsum : ∑ i, ![(1 : ℝ), 1, c₃] i = 2 + c₃ := by
    rw [Fin.sum_univ_three]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
      Matrix.cons_val_two, Matrix.tail_cons]
    ring
  have hden : (0 : ℝ) < c₃ + 17 * Real.sqrt c₃ + 72 := by linarith
  have hd1 : (0 : ℝ) < (8 + Real.sqrt c₃) * ((8 + Real.sqrt c₃) + 1) := by nlinarith
  have hlim : stackSVDLimit (fun i => (m.tbl i).θ) ![(1 : ℝ), 1, c₃]
      = (16 * Real.sqrt c₃ + 62) / (c₃ + 17 * Real.sqrt c₃ + 72) := by
    unfold stackSVDLimit
    rw [hsum, hcsum, if_pos (by nlinarith : (8 + Real.sqrt c₃) ^ 2 > 2 + c₃),
      div_eq_div_iff hd1.ne' hden.ne']
    linear_combination (c₃ + Real.sqrt c₃ + 10) * hs
  rwa [hlim] at h

/-- The `cor.2` binary weighting of the paper's `M = 3` example discards table 3, because
`θ₃⁴ = c₃` is not above `c₃`. Three facts: the detectable set is `{0, 1}`; the binary-weighted
stack of tables `0` and `1` converges to `62/72`; and `62/72` is above svdstack's `6/7`. So
the paper's example does not separate the two methods in the direction its heading claims (see
the module docstring and `paper_edits.md` E6). -/
theorem remark_svd_outperform_stack_three_binary (m : MultiTableModel μ 3 n₃ d) (c₃ : ℝ)
    (hc₃ : 0 < c₃) (hθ : ∀ i, (m.tbl i).θ = ![2, 2, Real.sqrt (Real.sqrt c₃)] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 1, c₃] i)) (hG : m.JointGaussianNoise) :
    (Finset.univ.filter fun i : Fin 3 =>
        ![(1 : ℝ), 1, c₃] i < (![2, 2, Real.sqrt (Real.sqrt c₃)] i) ^ 4) = {0, 1}
    ∧ TendstoInProb μ (fun N ω =>
        m.stackPerfW (fun i => if i ∈ ({0, 1} : Finset (Fin 3)) then 1 else 0) N ω) (62 / 72)
    ∧ (6 : ℝ) / 7 < 62 / 72 := by
  have hq : Real.sqrt (Real.sqrt c₃) ^ 4 = c₃ := sqrt_sqrt_pow_four hc₃.le
  have hθ0 : (m.tbl 0).θ = 2 := by rw [hθ]; simp
  have hθ1 : (m.tbl 1).θ = 2 := by rw [hθ]; simp
  have hcpos : ∀ i : Fin 3, 0 < ![(1 : ℝ), 1, c₃] i := by
    intro i
    fin_cases i <;> simp [hc₃]
  refine ⟨?_, ?_, by norm_num⟩
  · ext i
    fin_cases i <;> simp [hq] <;> norm_num
  · have hcard : NeZero (({0, 1} : Finset (Fin 3)).card) := ⟨by decide⟩
    have h := m.stackPerfW_binary_tendsto_gaussian ({0, 1} : Finset (Fin 3))
      ![(1 : ℝ), 1, c₃] hcpos hreg hG
    have hval : Scalars.binaryStackSVDLimit ({0, 1} : Finset (Fin 3))
        (fun i => (m.tbl i).θ) ![(1 : ℝ), 1, c₃] = 62 / 72 := by
      simp only [Scalars.binaryStackSVDLimit,
        Finset.sum_pair (show (0 : Fin 3) ≠ 1 by decide), hθ0, hθ1]
      norm_num
    rwa [hval] at h

end Three

/-! ### `remark:svd_outperform_stack`, instance (iv): `M = 2`, `θ = (√5, 4)`, `c = (1, 38.4)` -/

section Two

variable {n₂ : Fin 2 → ℕ → ℕ}

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:610`), the paper's two-table example.
Both tables have `β² = 4/5`. Binary-weighted stacksvd on both tables reaches
`2008/2310 = 0.8693`, and that is the maximum over the three nonempty subsets, since each
single table gives `4/5`. Optimally weighted svdstack reaches `S/(S+1) = 8/9 = 0.8889`, which
is strictly more. The unweighted stack is the binary weighting on `Finset.univ`, so the second
conjunct is the third one at `S = Finset.univ`. -/
theorem remark_svd_outperform_stack_two (m : MultiTableModel μ 2 n₂ d)
    (hθ : ∀ i, (m.tbl i).θ = ![Real.sqrt 5, 4] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 38.4] i)) (hG : m.JointGaussianNoise) :
    (∀ i, beta ((m.tbl i).θ) (![1, 38.4] i) ^ 2 = 4 / 5)
    ∧ TendstoInProb μ (fun N ω => m.stackPerfW 1 N ω) (2008 / 2310)
    ∧ (∀ S : Finset (Fin 2), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
          (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4])
        ∧ Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4] ≤ 2008 / 2310)
    ∧ TendstoInProb μ (fun N ω =>
        m.svdstackPerfW (optW fun i => beta ((m.tbl i).θ) (![1, 38.4] i)) N ω) (8 / 9)
    ∧ (2008 : ℝ) / 2310 < 8 / 9 := by
  have h5 : Real.sqrt 5 ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  have h5' : Real.sqrt 5 ^ 4 = 25 := by
    have h : Real.sqrt 5 ^ 4 = (Real.sqrt 5 ^ 2) ^ 2 := by ring
    rw [h, h5]
    norm_num
  have hθ0 : (m.tbl 0).θ = Real.sqrt 5 := by rw [hθ]; simp
  have hθ1 : (m.tbl 1).θ = 4 := by rw [hθ]; simp
  have hcpos : ∀ i : Fin 2, 0 < ![(1 : ℝ), 38.4] i := by
    intro i
    fin_cases i <;> norm_num
  -- both tables have `β² = 4/5`
  have hb0 : beta (Real.sqrt 5) 1 ^ 2 = 4 / 5 := by
    rw [Scalars.beta_sq, betaSq, if_pos (by rw [h5']; norm_num : Real.sqrt 5 ^ 4 > (1 : ℝ)),
      h5, h5']
    norm_num
  have hb1 : beta 4 38.4 ^ 2 = 4 / 5 := by
    rw [Scalars.beta_sq, betaSq, if_pos (by norm_num : (4 : ℝ) ^ 4 > 38.4)]
    norm_num
  have hbsq : ∀ i, beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i) ^ 2 = 4 / 5 := by
    intro i
    fin_cases i
    · change beta ((m.tbl 0).θ) (![(1 : ℝ), 38.4] 0) ^ 2 = 4 / 5
      rw [hθ0]
      simpa using hb0
    · change beta ((m.tbl 1).θ) (![(1 : ℝ), 38.4] 1) ^ 2 = 4 / 5
      rw [hθ1]
      simpa using hb1
  -- the `cor.2` convergence on every nonempty subset
  have hbin : ∀ S : Finset (Fin 2), S.Nonempty →
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
        (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![(1 : ℝ), 38.4]) := by
    intro S hS
    have hcard : NeZero S.card := ⟨(Finset.card_pos.mpr hS).ne'⟩
    exact m.stackPerfW_binary_tendsto_gaussian S ![(1 : ℝ), 38.4] hcpos hreg hG
  -- the value on `Finset.univ`
  have hunivval : Scalars.binaryStackSVDLimit (Finset.univ : Finset (Fin 2))
      (fun i => (m.tbl i).θ) ![(1 : ℝ), 38.4] = 2008 / 2310 := by
    simp only [Scalars.binaryStackSVDLimit, Fin.sum_univ_two, hθ0, hθ1]
    rw [h5, if_pos (by norm_num : ((5 : ℝ) + 4 ^ 2) ^ 2 > ![(1 : ℝ), 38.4] 0
      + ![(1 : ℝ), 38.4] 1)]
    norm_num
  -- the three nonempty subsets, and the maximum
  have hSenum : ∀ S : Finset (Fin 2), S.Nonempty → S = {0} ∨ S = {1} ∨ S = {0, 1} := by
    decide
  have hle : ∀ S : Finset (Fin 2), S.Nonempty →
      Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![(1 : ℝ), 38.4]
        ≤ 2008 / 2310 := by
    intro S hS
    rcases hSenum S hS with rfl | rfl | rfl
    · simp only [Scalars.binaryStackSVDLimit, Finset.sum_singleton, hθ0]
      rw [h5]
      norm_num
    · simp only [Scalars.binaryStackSVDLimit, Finset.sum_singleton, hθ1]
      norm_num
    · simp only [Scalars.binaryStackSVDLimit,
        Finset.sum_pair (show (0 : Fin 2) ≠ 1 by decide), hθ0, hθ1]
      rw [h5]
      norm_num
  -- the unweighted stack is the binary weighting on `Finset.univ`
  have huniv : TendstoInProb μ (fun N ω => m.stackPerfW 1 N ω) (2008 / 2310) := by
    have h := hbin Finset.univ Finset.univ_nonempty
    have hw : (fun i : Fin 2 => if i ∈ (Finset.univ : Finset (Fin 2)) then (1 : ℝ) else 0)
        = 1 := by
      funext i
      simp
    rw [hw, hunivval] at h
    exact h
  -- optimally weighted svdstack
  have hopt : TendstoInProb μ (fun N ω =>
      m.svdstackPerfW (optW fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)) N ω)
      (8 / 9) := by
    have hthr : ∃ k : Fin 2, 0 < beta ((m.tbl k).θ) (![(1 : ℝ), 38.4] k) := by
      refine ⟨0, ?_⟩
      rw [hθ0]
      simp only [Matrix.cons_val_zero]
      exact beta_pos_of_thr one_pos (by rw [h5']; norm_num)
    have h := (m.thm_svdstack_weighted_gaussian_opt ![(1 : ℝ), 38.4]
      (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)) hcpos (fun i => rfl) hthr hreg hG).1
    have hSv : Sval (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)) = 8 := by
      simp only [Sval, Fin.sum_univ_two, hbsq]
      norm_num
    have hval : svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i))
        = 8 / 9 := by
      rw [svdstackLimitOpt, hSv]
      norm_num
    rwa [hval] at h
  exact ⟨hbsq, huniv, fun S hS => ⟨hbin S hS, hle S hS⟩, hopt, by norm_num⟩

end Two

end MultiTableModel

end StackedSVD
