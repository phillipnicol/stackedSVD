/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R3het
import StackedSVD.RMT.Het.R4het
import StackedSVD.RMT.Het.EdgeScalar

/-!
# Campaign E: the sharp bound for the heteroscedastic edge, and the edge itself

Tasks A, C and D of `notes/archive/heteroedge_sharp.md` (decision D28): steps E4a, E4b and E3
(task A, deterministic), then steps E5 and E6 (task C, the integral bound), then steps E8 and
E9 (task D, the model theorem `MultiTableModel.heteroEdge_of_gaussian`). The file has no
`sorry` and no `axiom`. Task A, in the first half, uses no probability and no measure.

## Mathematics

Let `S = diag(a)` on `ℝ^P`, `σ_j = a_j²` and fix `γ` with `σ_j < γ` for every `j`. The
Sudakov-Fernique comparison process of `RMT/Het/R3het.lean` is
`Y_i = ⟪S x_i, g⟫ + ‖S x_i‖ ⟪y_i, h⟫` at unit `x_i`, `y_i`. The current pointwise bound
(`integral_iSup_bilin_le_het`) splits it by the triangle inequality and pays `‖S g‖ + L ‖h‖`,
which gives `bSF`. The sharp route bounds it by a ball instead:

1. `⟪S x, g⟫ + ‖S x‖ ‖h‖ = ⟪S x, g + q₀⟫ ≤ ‖S (g + q₀)‖` at the explicit witness
   `q₀ = (‖h‖/‖S x‖) S x` of the ball `‖q‖ ≤ ‖h‖`, because `Sᵀ = S` and `‖x‖ = 1`.
2. `‖S (g + q)‖² = ∑_j σ_j (g_j + q_j)²` and, coordinatewise,
   `σ_j (g_j + q_j)² = [(σ_j - γ) q_j² + 2 σ_j g_j q_j + σ_j g_j²] + γ q_j²`. The bracket is a
   concave quadratic in `q_j` with maximum `σ_j γ/(γ - σ_j) g_j²` (step E4a), so
   `‖S (g + q)‖² ≤ γ ‖q‖² + ∑_j σ_j γ/(γ - σ_j) g_j² ≤ γ A² + ∑_j σ_j γ/(γ - σ_j) g_j²`
   whenever `‖q‖ ≤ A` (step E4b).

Step E3 puts the two together and is the pointwise replacement for the `hpt` step of
`integral_iSup_bilin_le_het` (`RMT/Het/R3het.lean`, around line 278). Only the `≤` direction
is proved; no sup identity is needed (modeling choice 2 of the plan note).

The radicand `γ ‖v‖² + ∑_j σ_j γ/(γ - σ_j) u_j²` is quadratic in the Gaussian coordinates
`(u, v)`, so task C integrates it by second moments and closes with `∫ √Z ≤ √(∫ Z)`.

STATUS 2026-08-30: `lake env lean -j 3 StackedSVD/RMT/Het/EdgeSharp.lean` exit 0, 0 `sorry`,
no `axiom`. See `notes/archive/agent_reports/edge_taskA.md`,
`notes/archive/agent_reports/edge_taskC.md` and `notes/archive/agent_reports/edge_taskD.md`.
-/

open Filter Topology
open scoped Matrix Matrix.Norms.L2Operator

namespace StackedSVD

variable {P p : ℕ}

/-! ### E4a: the scalar complete square -/

set_option linter.unusedVariables false in
/-- **E4a.** For `0 ≤ σ < γ` the concave quadratic `q ↦ (σ - γ) q² + 2 σ h q + σ h²` is at most
its maximum `σ γ/(γ - σ) h²`, attained at `q = σ h/(γ - σ)`.

The proof multiplies by `γ - σ > 0` and reads off the identity
`σ γ h² - ((σ - γ) q² + 2 σ h q + σ h²) (γ - σ) = ((γ - σ) q - σ h)²`, so `hσ : 0 ≤ σ` is not
needed. The hypothesis stays because the plan note and every call site carry it
(`σ = a_j²` there); the linter warning about it is switched off for this declaration only. -/
theorem quad_complete_square_le {σ γ : ℝ} (hσ : 0 ≤ σ) (hγ : σ < γ) (q h : ℝ) :
    (σ - γ) * q ^ 2 + 2 * σ * h * q + σ * h ^ 2 ≤ σ * γ / (γ - σ) * h ^ 2 := by
  have hd : (0 : ℝ) < γ - σ := sub_pos.mpr hγ
  have hform : σ * γ / (γ - σ) * h ^ 2 = σ * γ * h ^ 2 / (γ - σ) := by
    rw [div_mul_eq_mul_div]
  rw [hform, le_div_iff₀ hd]
  nlinarith [sq_nonneg ((γ - σ) * q - σ * h), sq_nonneg h, hσ, hd]

/-! ### E4b: the summed bound on a ball -/

/-- **E4b.** With `S = diag(a)`, `σ_j = a_j² < γ` and `‖q‖ ≤ A`:
`‖S (g + q)‖² ≤ γ A² + ∑_j σ_j γ/(γ - σ_j) g_j²`. The slack `γ (A² - ‖q‖²) ≥ 0` is added first,
then E4a is applied in every coordinate. -/
theorem weighted_ball_norm_sq_le {a : Fin P → ℝ} {γ : ℝ} (hγ0 : 0 ≤ γ)
    (ha : ∀ j, (a j) ^ 2 < γ) (g q : EuclideanSpace ℝ (Fin P)) {A : ℝ} (hq : ‖q‖ ≤ A) :
    ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (g + q)) :
        EuclideanSpace ℝ (Fin P))‖ ^ 2
      ≤ γ * A ^ 2 + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp g j) ^ 2 := by
  -- the norm, coordinate by coordinate
  have hcoord : ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (g + q)) :
        EuclideanSpace ℝ (Fin P))‖ ^ 2
      = ∑ j, (a j) ^ 2 * (WithLp.ofLp g j + WithLp.ofLp q j) ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
    refine Finset.sum_congr rfl fun j _ => ?_
    have hj : (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (g + q)) :
        EuclideanSpace ℝ (Fin P)) j = a j * (WithLp.ofLp g j + WithLp.ofLp q j) := by
      change (Matrix.diagonal a *ᵥ WithLp.ofLp (g + q)) j = _
      rw [Matrix.mulVec_diagonal, WithLp.ofLp_add]
      rfl
    rw [hj]
    ring
  -- the ball constraint, in coordinates
  have hA0 : (0 : ℝ) ≤ A := le_trans (norm_nonneg q) hq
  have hq2 : ∑ j, (WithLp.ofLp q j) ^ 2 ≤ A ^ 2 := by
    have hsq : ‖q‖ ^ 2 ≤ A ^ 2 := by nlinarith [norm_nonneg q]
    rwa [EuclideanSpace.real_norm_sq_eq] at hsq
  -- E4a in every coordinate
  have hsum : ∑ j, (a j) ^ 2 * (WithLp.ofLp g j + WithLp.ofLp q j) ^ 2
      ≤ γ * (∑ j, (WithLp.ofLp q j) ^ 2)
        + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp g j) ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_le_sum fun j _ => ?_
    have hj := quad_complete_square_le (sq_nonneg (a j)) (ha j)
      (WithLp.ofLp q j) (WithLp.ofLp g j)
    have hexp : (a j) ^ 2 * (WithLp.ofLp g j + WithLp.ofLp q j) ^ 2
        = ((a j) ^ 2 - γ) * (WithLp.ofLp q j) ^ 2
          + 2 * (a j) ^ 2 * (WithLp.ofLp g j) * (WithLp.ofLp q j)
          + (a j) ^ 2 * (WithLp.ofLp g j) ^ 2 + γ * (WithLp.ofLp q j) ^ 2 := by ring
    rw [hexp]
    linarith
  have hslack : γ * (∑ j, (WithLp.ofLp q j) ^ 2) ≤ γ * A ^ 2 :=
    mul_le_mul_of_nonneg_left hq2 hγ0
  rw [hcoord]
  linarith

/-! ### E3: the pointwise bound on the comparison process -/

/-- The radicand of E3 is nonnegative: every coefficient `σ_j γ/(γ - σ_j)` is nonnegative
because `0 ≤ σ_j < γ`. -/
theorem edge_radicand_nonneg {a : Fin P → ℝ} {γ : ℝ} (hγ0 : 0 ≤ γ) (ha : ∀ j, (a j) ^ 2 < γ)
    (t : ℝ) (ht : 0 ≤ t) (u : Fin P → ℝ) :
    0 ≤ γ * t + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (u j) ^ 2 := by
  have h1 : (0 : ℝ) ≤ γ * t := mul_nonneg hγ0 ht
  have h2 : (0 : ℝ) ≤ ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (u j) ^ 2 :=
    Finset.sum_nonneg fun j _ => by
      have hd : (0 : ℝ) < γ - (a j) ^ 2 := sub_pos.mpr (ha j)
      exact mul_nonneg (div_nonneg (mul_nonneg (sq_nonneg _) hγ0) hd.le) (sq_nonneg _)
  linarith

/-- **E3**, workhorse form with both test vectors in `EuclideanSpace`:
`⟪S x, U⟫ + ‖S x‖ ⟪y, V⟫ ≤ √(γ ‖V‖² + ∑_j σ_j γ/(γ - σ_j) U_j²)` with `S = diag(a)` and
`x`, `y` unit. Cauchy-Schwarz on the `y` block, then the explicit witness
`q₀ = (‖V‖/‖S x‖) S x` of the ball of radius `‖V‖`, then E4b. -/
theorem bilin_cmp_le_ball' {a : Fin P → ℝ} {γ : ℝ} (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ)
    {x : EuclideanSpace ℝ (Fin P)} (hx : ‖x‖ = 1)
    {y : EuclideanSpace ℝ (Fin p)} (hy : ‖y‖ = 1)
    (U : EuclideanSpace ℝ (Fin P)) (V : EuclideanSpace ℝ (Fin p)) :
    inner ℝ (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin P)) U
      + ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin P))‖
        * inner ℝ y V
      ≤ Real.sqrt (γ * ‖V‖ ^ 2
          + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp U j) ^ 2) := by
  classical
  set xs : EuclideanSpace ℝ (Fin P) :=
    WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) with hxs
  set R : ℝ :=
    γ * ‖V‖ ^ 2 + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp U j) ^ 2 with hR
  have hR0 : 0 ≤ R := by
    rw [hR]
    exact edge_radicand_nonneg hγ0.le ha _ (sq_nonneg ‖V‖) (WithLp.ofLp U)
  -- (a) Cauchy-Schwarz on the `y` block
  have hyv : inner ℝ y V ≤ ‖V‖ := by
    calc inner ℝ y V ≤ ‖y‖ * ‖V‖ := real_inner_le_norm _ _
      _ = ‖V‖ := by rw [hy, one_mul]
  have hstep : inner ℝ xs U + ‖xs‖ * inner ℝ y V ≤ inner ℝ xs U + ‖xs‖ * ‖V‖ := by
    have h := mul_le_mul_of_nonneg_left hyv (norm_nonneg xs)
    linarith
  -- `diag a` is symmetric, so `⟪S x, w⟫ = ⟪x, S w⟫`
  have hmove : ∀ w : EuclideanSpace ℝ (Fin P),
      inner ℝ xs w
        = inner ℝ x (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp w) :
            EuclideanSpace ℝ (Fin P)) := by
    intro w
    have hdot : ∀ r s : Fin P → ℝ, ∑ j, (Matrix.diagonal a *ᵥ r) j * s j
        = ∑ j, r j * (Matrix.diagonal a *ᵥ s) j := by
      intro r s
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [Matrix.mulVec_diagonal, Matrix.mulVec_diagonal]
      ring
    rw [hxs, ← R3.sum_mul_eq_inner, ← R3.sum_mul_eq_inner]
    exact hdot (WithLp.ofLp x) (WithLp.ofLp w)
  -- the ball bound at an arbitrary point of the ball
  have hball : ∀ q : EuclideanSpace ℝ (Fin P), ‖q‖ ≤ ‖V‖ →
      inner ℝ xs (U + q) ≤ Real.sqrt R := by
    intro q hqn
    have h2 : ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
        EuclideanSpace ℝ (Fin P))‖ ^ 2 ≤ R := by
      rw [hR]
      exact weighted_ball_norm_sq_le hγ0.le ha U q hqn
    have key : ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
        EuclideanSpace ℝ (Fin P))‖ ≤ Real.sqrt R := by
      calc ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
            EuclideanSpace ℝ (Fin P))‖
          = Real.sqrt (‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
              EuclideanSpace ℝ (Fin P))‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
        _ ≤ Real.sqrt R := Real.sqrt_le_sqrt h2
    calc inner ℝ xs (U + q)
        = inner ℝ x (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
            EuclideanSpace ℝ (Fin P)) := hmove (U + q)
      _ ≤ ‖x‖ * ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
            EuclideanSpace ℝ (Fin P))‖ := real_inner_le_norm _ _
      _ = ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (U + q)) :
            EuclideanSpace ℝ (Fin P))‖ := by rw [hx, one_mul]
      _ ≤ Real.sqrt R := key
  -- (b) the degenerate case `S x = 0`, and (c) the witness
  rcases eq_or_ne xs 0 with h0 | h0
  · have hzero : inner ℝ xs U + ‖xs‖ * inner ℝ y V = 0 := by
      rw [h0]; simp
    rw [hzero]
    exact Real.sqrt_nonneg R
  · have hxs0 : 0 < ‖xs‖ := norm_pos_iff.mpr h0
    have hne : ‖xs‖ ≠ 0 := ne_of_gt hxs0
    have hcancel : ‖V‖ / ‖xs‖ * ‖xs‖ = ‖V‖ := by field_simp
    set q₀ : EuclideanSpace ℝ (Fin P) := (‖V‖ / ‖xs‖) • xs with hq₀
    have hq₀n : ‖q₀‖ = ‖V‖ := by
      rw [hq₀, norm_smul, Real.norm_eq_abs,
        abs_of_nonneg (div_nonneg (norm_nonneg V) (norm_nonneg xs)), hcancel]
    have hinner : inner ℝ xs q₀ = ‖xs‖ * ‖V‖ := by
      rw [hq₀, real_inner_smul_right, real_inner_self_eq_norm_mul_norm, ← mul_assoc, hcancel]
      ring
    calc inner ℝ xs U + ‖xs‖ * inner ℝ y V
        ≤ inner ℝ xs U + ‖xs‖ * ‖V‖ := hstep
      _ = inner ℝ xs (U + q₀) := by rw [inner_add_right, hinner]
      _ ≤ Real.sqrt R := hball q₀ (le_of_eq hq₀n)

/-- **E3.** The sharp pointwise bound on the Sudakov-Fernique comparison process of
`RMT/Het/R3het.lean`, with the two Gaussian blocks given as plain coordinate functions
`u : Fin P → ℝ` and `v : Fin p → ℝ`:
`⟪S x, u⟫ + ‖S x‖ ⟪y, v⟫ ≤ √(γ ‖v‖² + ∑_j σ_j γ/(γ - σ_j) u_j²)`. -/
theorem bilin_cmp_le_ball {a : Fin P → ℝ} {γ : ℝ} (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ)
    {x : EuclideanSpace ℝ (Fin P)} (hx : ‖x‖ = 1)
    {y : EuclideanSpace ℝ (Fin p)} (hy : ‖y‖ = 1) (u : Fin P → ℝ) (v : Fin p → ℝ) :
    inner ℝ (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin P))
        (WithLp.toLp 2 u : EuclideanSpace ℝ (Fin P))
      + ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin P))‖
        * inner ℝ y (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))
      ≤ Real.sqrt (γ * ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))‖ ^ 2
          + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (u j) ^ 2) :=
  bilin_cmp_le_ball' hγ0 ha hx hy (WithLp.toLp 2 u) (WithLp.toLp 2 v)

/-- **E3, sum form.** Both inner products written as coordinate sums and `‖v‖²` written out.
This is the shape the `hpt` step of task C meets after `R3het.cmpMatHet_apply`, with
`u := fun a => w (Sum.inl a)` and `v := fun b => w (Sum.inr b)`. -/
theorem bilin_cmp_le_ball_sum {a : Fin P → ℝ} {γ : ℝ} (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ)
    {x : EuclideanSpace ℝ (Fin P)} (hx : ‖x‖ = 1)
    {y : EuclideanSpace ℝ (Fin p)} (hy : ‖y‖ = 1) (u : Fin P → ℝ) (v : Fin p → ℝ) :
    (∑ j, WithLp.ofLp (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) :
        EuclideanSpace ℝ (Fin P)) j * u j)
      + ‖(WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin P))‖
        * (∑ b, WithLp.ofLp y b * v b)
      ≤ Real.sqrt (γ * (∑ b, (v b) ^ 2)
          + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (u j) ^ 2) := by
  have h := bilin_cmp_le_ball (a := a) (γ := γ) hγ0 ha hx hy u v
  have hu : ∑ j, WithLp.ofLp (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) :
        EuclideanSpace ℝ (Fin P)) j * u j
      = inner ℝ (WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp x) :
          EuclideanSpace ℝ (Fin P)) (WithLp.toLp 2 u : EuclideanSpace ℝ (Fin P)) := by
    rw [← R3.sum_mul_eq_inner]
  have hv : ∑ b, WithLp.ofLp y b * v b
      = inner ℝ y (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p)) := by
    rw [← R3.sum_mul_eq_inner]
  have hnv : ‖(WithLp.toLp 2 v : EuclideanSpace ℝ (Fin p))‖ ^ 2 = ∑ b, (v b) ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
  rw [hu, hv, ← hnv]
  exact h

end StackedSVD

/-! ## Campaign E, task C: the sharp integral bound

Steps E5 and E6 of `notes/archive/heteroedge_sharp.md`. The Sudakov-Fernique route of
`RMT/Het/R3het.lean` is unchanged up to the pointwise bound on the comparison process; only
the last two steps differ. Where `integral_iSup_bilin_le_het` splits the comparison process by
the triangle inequality and pays `√tr(S Sᵀ) + L √p`, task C applies E3 (`bilin_cmp_le_ball_sum`)
and pays `√(γ p + ∑_j σ_j γ/(γ - σ_j))`, which is smaller for a good `γ`.

1. **E5.** The pointwise bound `⨆_i Y_i(w) ≤ √(Z(w))` with the quadratic radicand
   `Z(w) = γ ∑_b w(inr b)² + ∑_j σ_j γ/(γ - σ_j) w(inl j)²`. The bound does not depend on `i`,
   so `ciSup_le` closes it.
2. **E6.** `∫ √Z ≤ √(∫ Z)` (Jensen for the concave `√` on `[0, ∞)`), and `∫ Z = γ p + ∑_j
   σ_j γ/(γ - σ_j)` by the coordinate second moments of the standard Gaussian.
3. The net argument of `R3het.integral_opNorm_mul_le_aux` is copied with the new constant.
-/

open MeasureTheory ProbabilityTheory

namespace StackedSVD

variable {P p : ℕ}

/-! ### E6a: coordinate second moments of the standard Gaussian -/

/-- One coordinate of a standard Gaussian vector is square integrable. -/
theorem integrable_coord_sq_stdGaussian {ι : Type*} [Fintype ι] (k : ι) :
    Integrable (fun w : EuclideanSpace ℝ ι => (WithLp.ofLp w k) ^ 2)
      (stdGaussian (EuclideanSpace ℝ ι)) := by
  classical
  rw [← multivariateGaussian_zero_one (ι := ι)]
  exact (memLp_two_eval_multivariateGaussian (0 : EuclideanSpace ℝ ι)
    (1 : Matrix ι ι ℝ) k).integrable_sq

/-- `E g_k² = 1` for a standard Gaussian vector: the variance of a coordinate is the diagonal
entry `1` of the covariance, and the mean is `0`. -/
theorem integral_coord_sq_stdGaussian {ι : Type*} [Fintype ι] (k : ι) :
    ∫ w, (WithLp.ofLp w k) ^ 2 ∂(stdGaussian (EuclideanSpace ℝ ι)) = 1 := by
  classical
  have h0 : ∫ w : EuclideanSpace ℝ ι, WithLp.ofLp w k
      ∂(multivariateGaussian 0 (1 : Matrix ι ι ℝ)) = 0 := by
    simpa using integral_eval_multivariateGaussian (0 : EuclideanSpace ℝ ι) (1 : Matrix ι ι ℝ) k
  have hv := variance_eval_multivariateGaussian (μ := (0 : EuclideanSpace ℝ ι))
    (S := (1 : Matrix ι ι ℝ)) Matrix.PosSemidef.one k
  rw [variance_of_integral_eq_zero (by fun_prop) h0] at hv
  rw [← multivariateGaussian_zero_one (ι := ι)]
  simpa [Matrix.one_apply_eq] using hv

/-- The radicand of E3 is integrable against the standard Gaussian of `ℝ^P ⊕ ℝ^p`: it is a
finite sum of coordinate squares with constant coefficients. -/
theorem integrable_edge_radicand (a : Fin P → ℝ) (γ : ℝ) :
    Integrable (fun w : EuclideanSpace ℝ (Fin P ⊕ Fin p) =>
      γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
        + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2)
      (stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) := by
  refine Integrable.add ?_ ?_
  · exact (integrable_finsetSum _ fun b _ =>
      integrable_coord_sq_stdGaussian (Sum.inr b : Fin P ⊕ Fin p)).const_mul γ
  · exact integrable_finsetSum _ fun j _ =>
      (integrable_coord_sq_stdGaussian (Sum.inl j : Fin P ⊕ Fin p)).const_mul _

/-- **E6a.** `E Z = γ p + ∑_j σ_j γ/(γ - σ_j)`: every coordinate square has mean `1`. -/
theorem integral_edge_radicand (a : Fin P → ℝ) (γ : ℝ) :
    ∫ w : EuclideanSpace ℝ (Fin P ⊕ Fin p),
        (γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
          + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2)
        ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p)))
      = γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) := by
  rw [integral_add
    ((integrable_finsetSum _ fun b _ =>
      integrable_coord_sq_stdGaussian (Sum.inr b : Fin P ⊕ Fin p)).const_mul γ)
    (integrable_finsetSum _ fun j _ =>
      (integrable_coord_sq_stdGaussian (Sum.inl j : Fin P ⊕ Fin p)).const_mul _),
    integral_const_mul,
    integral_finsetSum _ (fun b _ => integrable_coord_sq_stdGaussian (Sum.inr b : Fin P ⊕ Fin p)),
    integral_finsetSum _ (fun j _ =>
      (integrable_coord_sq_stdGaussian (Sum.inl j : Fin P ⊕ Fin p)).const_mul _)]
  simp only [integral_const_mul, integral_coord_sq_stdGaussian, mul_one]
  simp

/-! ### E6b: Jensen for the square root -/

/-- The square root of a nonnegative integrable function is integrable on a finite measure,
because `√t ≤ 1 + t`. -/
theorem integrable_sqrt_comp {α : Type*} [MeasurableSpace α] {μ : Measure α}
    [IsFiniteMeasure μ] {f : α → ℝ} (hf0 : ∀ z, 0 ≤ f z) (hf : Integrable f μ) :
    Integrable (fun z => Real.sqrt (f z)) μ := by
  have hmeas : AEStronglyMeasurable (fun z => Real.sqrt (f z)) μ :=
    Real.continuous_sqrt.comp_aestronglyMeasurable hf.aestronglyMeasurable
  refine Integrable.mono' ((integrable_const (1 : ℝ)).add hf) hmeas
    (Filter.Eventually.of_forall fun z => ?_)
  simp only [Pi.add_apply]
  rw [Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _)]
  nlinarith [Real.sq_sqrt (hf0 z), Real.sqrt_nonneg (f z), sq_nonneg (Real.sqrt (f z) - 1)]

/-- **E6b.** `∫ √f ≤ √(∫ f)` for a probability measure and a nonnegative integrable `f`.
Jensen's inequality (`ConcaveOn.le_map_integral`) at the concave `√` on `[0, ∞)`. -/
theorem integral_sqrt_le_sqrt_integral {α : Type*} [MeasurableSpace α] {μ : Measure α}
    [IsProbabilityMeasure μ] {f : α → ℝ} (hf0 : ∀ z, 0 ≤ f z) (hf : Integrable f μ) :
    ∫ z, Real.sqrt (f z) ∂μ ≤ Real.sqrt (∫ z, f z ∂μ) :=
  (Real.strictConcaveOn_sqrt.concaveOn).le_map_integral
    (f := f) (μ := μ) (Real.continuous_sqrt.continuousOn) isClosed_Ici
    (Filter.Eventually.of_forall fun z => Set.mem_Ici.mpr (hf0 z)) hf
    (integrable_sqrt_comp hf0 hf)

/-! ### E5: the sharp Sudakov-Fernique bound -/

/-- **E5.** The sharp expected-supremum bound for `S = diag(a)`: for every `γ` above every
`σ_j = a_j²`,
`E ⨆_i ⟪S x_i, B y_i⟫ ≤ √(γ p + ∑_j σ_j γ/(γ - σ_j))` on unit pairs `(x_i, y_i)`.

The proof follows `R3het.integral_iSup_bilin_le_het` up to the pointwise step: the same
comparison matrix `R3het.cmpMatHet`, the same increment inequality, the same
`sudakov_fernique` application and the same change of variables to the standard Gaussian of
`ℝ^P ⊕ ℝ^p`. The pointwise step is E3 (`bilin_cmp_le_ball_sum`), whose bound does not depend
on the index, and the last step is E6b together with the second moments of E6a. -/
theorem integral_iSup_bilin_le_sharp {ι : Type*} [Finite ι] [Nonempty ι]
    (a : Fin P → ℝ) {γ : ℝ} (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ)
    (x : ι → EuclideanSpace ℝ (Fin P)) (y : ι → EuclideanSpace ℝ (Fin p))
    (hx : ∀ i, ‖x i‖ = 1) (hy : ∀ i, ‖y i‖ = 1) :
    ∫ B, ⨆ i, (Matrix.diagonal a *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i))
        ∂(gaussianMatrix P p)
      ≤ Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) := by
  classical
  have := Fintype.ofFinite ι
  set xs : ι → EuclideanSpace ℝ (Fin P) :=
    fun i => WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (x i)) with hxs
  -- the increment inequality (identical to `R3het.integral_iSup_bilin_le_het`)
  have hincr : ∀ i j,
      (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) i i + (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) j j
        - 2 * (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) i j
      ≤ (R3het.cmpMatHet xs y * (R3het.cmpMatHet xs y)ᵀ) i i
        + (R3het.cmpMatHet xs y * (R3het.cmpMatHet xs y)ᵀ) j j
        - 2 * (R3het.cmpMatHet xs y * (R3het.cmpMatHet xs y)ᵀ) i j := by
    intro i j
    simp only [R3.bilMat_gram, R3het.cmpMatHet_gram]
    rw [R3.sum_mul_self_eq_norm_sq (xs i), R3.sum_mul_self_eq_norm_sq (xs j),
      R3.sum_mul_self_eq_norm_sq (y i), R3.sum_mul_self_eq_norm_sq (y j),
      R3.sum_mul_eq_inner (xs i) (xs j), R3.sum_mul_eq_inner (y i) (y j), hy i, hy j]
    have hs : inner ℝ (xs i) (xs j) ≤ ‖xs i‖ * ‖xs j‖ := real_inner_le_norm _ _
    have ht : |inner ℝ (y i) (y j)| ≤ 1 := by
      calc |inner ℝ (y i) (y j)| ≤ ‖y i‖ * ‖y j‖ := abs_real_inner_le_norm _ _
        _ = 1 := by rw [hy i, hy j]; norm_num
    obtain ⟨-, ht2⟩ := abs_le.mp ht
    have h1 : (0 : ℝ) ≤ 1 - inner ℝ (y i) (y j) := by linarith
    nlinarith [mul_le_mul_of_nonneg_right hs h1, sq_nonneg (‖xs i‖ - ‖xs j‖)]
  -- the left side, as an integral against `N(0, S)`
  have hSmap : (stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (R3.bilMat xs y)))
      = multivariateGaussian 0 (R3.bilMat xs y * (R3.bilMat xs y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin P × Fin p),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hleft : ∫ B, ⨆ i, (Matrix.diagonal a *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i))
        ∂(gaussianMatrix P p)
      = ∫ z, ⨆ i, WithLp.ofLp z i
          ∂(multivariateGaussian 0 (R3.bilMat xs y * (R3.bilMat xs y)ᵀ)) := by
    have h1 : ∫ B, ⨆ i, (Matrix.diagonal a *ᵥ WithLp.ofLp (x i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (y i))
          ∂(gaussianMatrix P p)
        = ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3.bilMat xs y) w) i
            ∂(stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))) := by
      rw [← (R3.measurePreserving_matrixEquivP P p).integral_comp'
        (fun w => ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3.bilMat xs y) w) i)]
      exact integral_congr_ae (Filter.Eventually.of_forall fun B =>
        iSup_congr fun i => (R3.bilMat_apply xs y B i).symm)
    rw [h1, ← hSmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    rfl
  -- the right side
  have hTmap : (stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))).map
      (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (R3het.cmpMatHet xs y)))
      = multivariateGaussian 0 (R3het.cmpMatHet xs y * (R3het.cmpMatHet xs y)ᵀ) := by
    rw [← multivariateGaussian_zero_one (ι := Fin P ⊕ Fin p),
      multivariateGaussian_zero_map_toEuclideanLin Matrix.PosSemidef.one, Matrix.mul_one]
  have hright : ∫ z, ⨆ i, WithLp.ofLp z i
        ∂(multivariateGaussian 0 (R3het.cmpMatHet xs y * (R3het.cmpMatHet xs y)ᵀ))
      ≤ Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) := by
    rw [← hTmap, integral_map (by fun_prop) measurable_ciSup_apply.aestronglyMeasurable]
    -- E5, the pointwise bound; the right side does not depend on the index
    have hpt : ∀ w : EuclideanSpace ℝ (Fin P ⊕ Fin p),
        (⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3het.cmpMatHet xs y) w) i)
          ≤ Real.sqrt (γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
              + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2) := by
      intro w
      refine ciSup_le fun i => ?_
      rw [R3het.cmpMatHet_apply]
      simpa only [hxs] using bilin_cmp_le_ball_sum (a := a) (γ := γ) hγ0 ha (hx i) (hy i)
        (fun j => WithLp.ofLp w (Sum.inl j)) (fun b => WithLp.ofLp w (Sum.inr b))
    have hZ0 : ∀ w : EuclideanSpace ℝ (Fin P ⊕ Fin p),
        0 ≤ γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
          + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2 := fun w =>
      edge_radicand_nonneg hγ0.le ha _
        (Finset.sum_nonneg fun b _ => sq_nonneg _) _
    calc ∫ w, ⨆ i, WithLp.ofLp (Matrix.toEuclideanLin (R3het.cmpMatHet xs y) w) i
          ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p)))
        ≤ ∫ w, Real.sqrt (γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
              + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2)
            ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p))) := by
          exact integral_mono
            (R3.integrable_ciSup_clm
              (LinearMap.toContinuousLinearMap (Matrix.toEuclideanLin (R3het.cmpMatHet xs y))) _)
            (integrable_sqrt_comp hZ0 (integrable_edge_radicand (p := p) a γ)) hpt
      _ ≤ Real.sqrt (∫ w, (γ * (∑ b : Fin p, (WithLp.ofLp w (Sum.inr b)) ^ 2)
              + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2) * (WithLp.ofLp w (Sum.inl j)) ^ 2)
            ∂(stdGaussian (EuclideanSpace ℝ (Fin P ⊕ Fin p)))) :=
          integral_sqrt_le_sqrt_integral hZ0 (integrable_edge_radicand (p := p) a γ)
      _ = Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) := by
          rw [integral_edge_radicand (p := p) a γ]
  rw [hleft]
  exact (ProbabilityTheory.sudakov_fernique (R3.posSemidef_bilMat_gram xs y)
    (R3het.posSemidef_cmpMatHet_gram xs y) hincr).trans hright


/-! ### E5b: the net argument, with the sharp constant -/

/-- The `e`-step of the sharp bound for `diag(a) * B`. The net argument is the one of
`R3het.integral_opNorm_mul_le_aux`: an `e`-net of the two unit spheres, both signs, and
`RMT.matrixOperatorNorm_le_of_centered_bilinear_net`. Only the constant changes, because the
finite family is fed to `integral_iSup_bilin_le_sharp` instead of
`R3het.integral_iSup_bilin_le_het`. -/
theorem integral_opNorm_mul_le_sharp_aux (hP : 0 < P) (hp : 0 < p) (a : Fin P → ℝ) {γ : ℝ}
    (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ) {e : ℝ} (he0 : 0 < e) (he : 2 * e < 1) :
    ∫ B, ‖Matrix.diagonal a * B‖ ∂(gaussianMatrix P p)
      ≤ Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) / (1 - 2 * e) := by
  classical
  have hS : (Matrix.diagonal a)ᵀ = Matrix.diagonal a := Matrix.diagonal_transpose a
  obtain ⟨Nx, hNx, hsubx⟩ := R3.exists_sphere_enet P hP he0
  obtain ⟨Ny, hNy, hsuby⟩ := R3.exists_sphere_enet p hp he0
  set Nnet : RMT.CenteredMatrixBilinearNet P p e :=
    { domainNet := Ny
      codomainNet := Nx
      domain_isNet := hNy
      codomain_isNet := hNx
      domain_subset_unitBall := hsuby.trans Metric.sphere_subset_closedBall
      codomain_subset_unitBall := hsubx.trans Metric.sphere_subset_closedBall } with hNnet
  obtain ⟨vx, hvx⟩ : Nx.Nonempty := Nnet.toMatrixBilinearNet.codomainNet_nonempty_of_pos hP
  obtain ⟨vy, hvy⟩ : Ny.Nonempty := Nnet.toMatrixBilinearNet.domainNet_nonempty_of_pos hp
  have : Nonempty (↥Nx × ↥Ny × Bool) := ⟨(⟨vx, hvx⟩, ⟨vy, hvy⟩, true)⟩
  set xf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin P) :=
    fun i => if i.2.2 then (i.1 : EuclideanSpace ℝ (Fin P))
      else -(i.1 : EuclideanSpace ℝ (Fin P)) with hxf
  set yf : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin p) :=
    fun i => (i.2.1 : EuclideanSpace ℝ (Fin p)) with hyf
  have hxn : ∀ i, ‖xf i‖ = 1 := by
    intro i
    have h1 : ‖(i.1 : EuclideanSpace ℝ (Fin P))‖ = 1 :=
      RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsubx i.1.2)
    by_cases hb : i.2.2 = true <;> simp [hxf, hb, h1]
  have hyn : ∀ i, ‖yf i‖ = 1 := fun i =>
    RMT.norm_eq_one_of_mem_euclideanUnitSphere (hsuby i.2.1.2)
  set u : Matrix (Fin P) (Fin p) ℝ → ℝ :=
    fun B => ⨆ i, (Matrix.diagonal a *ᵥ WithLp.ofLp (xf i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (yf i)) with hu
  have hbdd : ∀ B, BddAbove (Set.range fun i =>
      (Matrix.diagonal a *ᵥ WithLp.ofLp (xf i)) ⬝ᵥ (B *ᵥ WithLp.ofLp (yf i))) :=
    fun B => Finite.bddAbove_range _
  have hnet : ∀ B : Matrix (Fin P) (Fin p) ℝ, ∀ q ∈ Nnet.domainNet, ∀ b ∈ Nnet.codomainNet,
      |inner ℝ (Matrix.toEuclideanLin (Matrix.diagonal a * B) q) b| ≤ u B := by
    intro B q hq b hb
    have hinner : inner ℝ (Matrix.toEuclideanLin (Matrix.diagonal a * B) q) b
        = (Matrix.diagonal a *ᵥ WithLp.ofLp b) ⬝ᵥ (B *ᵥ WithLp.ofLp q) := by
      rw [inner_euclidean_eq_dotProduct, dotProduct_comm]
      change WithLp.ofLp b ⬝ᵥ ((Matrix.diagonal a * B) *ᵥ WithLp.ofLp q) = _
      rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hS]
    have hpos := le_ciSup (hbdd B) ((⟨b, hb⟩ : ↥Nx), (⟨q, hq⟩ : ↥Ny), true)
    have hneg := le_ciSup (hbdd B) ((⟨b, hb⟩ : ↥Nx), (⟨q, hq⟩ : ↥Ny), false)
    simp only [hxf, hyf, if_true] at hpos hneg
    rw [hinner]
    refine abs_le.mpr ⟨?_, hpos⟩
    have hneg' : -((Matrix.diagonal a *ᵥ WithLp.ofLp b) ⬝ᵥ (B *ᵥ WithLp.ofLp q)) ≤ u B := by
      refine le_trans (le_of_eq ?_) hneg
      simp [Matrix.mulVec_neg, neg_dotProduct]
    linarith
  have hu0 : ∀ B, 0 ≤ u B := by
    intro B
    have := hnet B vy hvy vx hvx
    exact le_trans (abs_nonneg _) this
  have hbound : ∀ B : Matrix (Fin P) (Fin p) ℝ,
      ‖Matrix.diagonal a * B‖ ≤ u B / (1 - 2 * e) := fun B =>
    RMT.matrixOperatorNorm_le_of_centered_bilinear_net (Matrix.diagonal a * B) Nnet he0.le he
      (hu0 B) (hnet B)
  set xs : ↥Nx × ↥Ny × Bool → EuclideanSpace ℝ (Fin P) :=
    fun i => WithLp.toLp 2 (Matrix.diagonal a *ᵥ WithLp.ofLp (xf i)) with hxs
  have hint : Integrable u (gaussianMatrix P p) := by
    have hg : Integrable
        (fun w : EuclideanSpace ℝ (Fin P × Fin p) =>
          ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
            (Matrix.toEuclideanLin (R3.bilMat xs yf)) w) i)
        (stdGaussian (EuclideanSpace ℝ (Fin P × Fin p))) :=
      R3.integrable_ciSup_clm _ _
    have heq : u = (fun w : EuclideanSpace ℝ (Fin P × Fin p) =>
        ⨆ i, WithLp.ofLp (LinearMap.toContinuousLinearMap
          (Matrix.toEuclideanLin (R3.bilMat xs yf)) w) i) ∘ (R3.matrixEquivP P p) := by
      funext B
      exact (iSup_congr fun i => R3.bilMat_apply xs yf B i).symm
    rw [heq]
    exact ((R3.measurePreserving_matrixEquivP P p).integrable_comp_emb
      (R3.matrixEquivP P p).measurableEmbedding).mpr hg
  have hden : (0 : ℝ) < 1 - 2 * e := by linarith
  calc ∫ B, ‖Matrix.diagonal a * B‖ ∂(gaussianMatrix P p)
      ≤ ∫ B, u B / (1 - 2 * e) ∂(gaussianMatrix P p) :=
        integral_mono_of_nonneg (Filter.Eventually.of_forall fun B => norm_nonneg _)
          (hint.div_const _) (Filter.Eventually.of_forall hbound)
    _ = (∫ B, u B ∂(gaussianMatrix P p)) / (1 - 2 * e) := integral_div _ _
    _ ≤ Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) / (1 - 2 * e) :=
        div_le_div_of_nonneg_right
          (integral_iSup_bilin_le_sharp a hγ0 ha xf yf hxn hyn) hden.le

/-- **E5b.** The sharp expected-operator-norm bound, after `e ↓ 0`:
`E ‖diag(a) B‖ ≤ √(γ p + ∑_j σ_j γ/(γ - σ_j))` for every `γ` above every `σ_j = a_j²`. -/
theorem integral_opNorm_mul_le_sharp (hP : 0 < P) (hp : 0 < p) (a : Fin P → ℝ) {γ : ℝ}
    (hγ0 : 0 < γ) (ha : ∀ j, (a j) ^ 2 < γ) :
    ∫ B, ‖Matrix.diagonal a * B‖ ∂(gaussianMatrix P p)
      ≤ Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) := by
  set A : ℝ := Real.sqrt (γ * p + ∑ j, (a j) ^ 2 * γ / (γ - (a j) ^ 2)) with hA
  have hcont : ContinuousAt (fun e : ℝ => A / (1 - 2 * e)) 0 := by
    refine ContinuousAt.div continuousAt_const (by fun_prop) (by norm_num)
  have hf : Tendsto (fun e : ℝ => A / (1 - 2 * e)) (𝓝[>] (0 : ℝ)) (𝓝 A) := by
    have h : Tendsto (fun e : ℝ => A / (1 - 2 * e)) (𝓝[>] (0 : ℝ))
        (𝓝 (A / (1 - 2 * (0 : ℝ)))) := hcont.continuousWithinAt
    simpa using h
  refine ge_of_tendsto hf ?_
  have hsmall : ∀ᶠ e in 𝓝[>] (0 : ℝ), e ∈ Set.Iio (1 / 2 : ℝ) :=
    Filter.Eventually.filter_mono nhdsWithin_le_nhds (Iio_mem_nhds (by norm_num))
  filter_upwards [self_mem_nhdsWithin, hsmall] with e he1 he2
  have he3 : e < 1 / 2 := he2
  exact integral_opNorm_mul_le_sharp_aux hP hp a hγ0 ha he1 (by linarith)

end StackedSVD

/-! ## Campaign E, task D: the assembly and the model theorem

Steps E8 and E9 of `notes/archive/heteroedge_sharp.md`. Two declarations:

1. `tendsto_measure_lamMax_le_of_bound`, the edge statement from an abstract expectation
   bound. It is `R3het.tendsto_measure_lamMax_le_het` with the fixed constant
   `√tr(S Sᵀ) + L √p` replaced by a sequence `bd N` with `bd N ² / d N → b`. The
   concentration input (`R3het.measure_opNorm_mul_ge_le`, Lipschitz constant `L + 1`) and the
   scalar step (`R3het.edge_arith_het` at `L := 0`) are the ones of that proof. The symmetry
   hypothesis `Sᵀ = S` disappears: it served the Gordon bound only.
2. `MultiTableModel.heteroEdge_of_gaussian`, the black box `HeteroEdge` of `RMT/Het/R4het.lean`
   proved for Gaussian noise at the exact edge `MPhet.bHet c w`. The bound sequence is the
   sharp one of task C at the fixed `γ⋆ = -1/s⋆`:

   `bd N = √(γ⋆ (d N - 1) + ∑_r a_r² γ⋆/(γ⋆ - a_r²))`, `a_r = w_{i(r)}`,

   and `bd N ² / d N = γ⋆ (d N - 1)/d N + ∑_i (n_i/d) w_i² γ⋆/(γ⋆ - w_i²)` is an identity, so
   the limit is `edgeObjective c w γ⋆ = bHet c w` (`MPhet.edgeObjective_sStar_eq_bHet`) by
   `(d - 1)/d → 1` and `n_i/d → c_i`. The regrouping of `∑_r` by table is the one of
   `MultiTableModel.trace_SigmaHalf_mul_transpose`, through `finSigmaFinEquiv`.

The plumbing of the model theorem (block law, `W₀'` on a block, row count, `p + 1 = d`) is
copied from `MultiTableModel.tendsto_measure_lamMax_W0het_le_bSF` (`RMT/Het/R3het.lean`).
-/

open scoped ENNReal NNReal

namespace StackedSVD

/-! ### E8: the edge from an abstract expectation bound -/

/-- **E8.** The upper edge from any expectation bound. If `E ‖S B‖ ≤ bd N` eventually,
`bd N ≥ 0`, `bd N ² / d N → b` and `‖S N‖ ≤ L`, then for every `ε > 0`
`lamMax (d⁻¹ (S B)(S B)ᵀ) ≤ b + ε` with probability tending to `1`.

This is `R3het.tendsto_measure_lamMax_le_het` with the constant `√tr(S Sᵀ) + L √p` replaced
by `bd N`; `L` survives only as the Lipschitz constant of the concentration step.

The column count is `p N + r = d N` for any fixed `r`. The rank-one heteroscedastic split
takes `r = 1`; the rank-`r` split of `RankR/Het/Split.lean` takes the same `r` as the model.
The proof reads `hpd` only to get `0 < p N` for large `N`. -/
theorem tendsto_measure_lamMax_le_of_bound
    {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] (μ : ∀ N, Measure (Ω N))
    [∀ N, IsProbabilityMeasure (μ N)] {P p d : ℕ → ℕ}
    (S : (N : ℕ) → Matrix (Fin (P N)) (Fin (P N)) ℝ) {L b : ℝ} (hL : 0 ≤ L) (hb : 0 ≤ b)
    (hSL : ∀ N, ‖S N‖ ≤ L) (bd : ℕ → ℝ) (hbd0 : ∀ N, 0 ≤ bd N)
    (hbd : ∀ᶠ N in atTop, (∫ B', ‖S N * B'‖ ∂(gaussianMatrix (P N) (p N))) ≤ bd N)
    (hlim : Tendsto (fun N => bd N ^ 2 / (d N : ℝ)) atTop (𝓝 b))
    (B : (N : ℕ) → Ω N → Matrix (Fin (P N)) (Fin (p N)) ℝ)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (P N)) (Fin (P N)) ℝ)
    (hW₀ : ∀ N ω, W₀ N ω = ((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ))
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (hlaw : ∀ N, HasLaw (B N) (gaussianMatrix (P N) (p N)) (μ N))
    (hP : ∀ N, 0 < P N) {r : ℕ} (hd : ∀ N, 0 < d N) (hpd : ∀ N, p N + r = d N)
    (hdtop : Tendsto d atTop atTop)
    {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}) atTop (𝓝 1) := by
  classical
  set sc : ℝ := Real.sqrt b with hscdef
  have hsc0 : (0 : ℝ) ≤ sc := Real.sqrt_nonneg _
  have hscsq : sc ^ 2 = b := Real.sq_sqrt hb
  set K : ℝ := 4 * sc + 4 with hKdef
  have hK0 : (0 : ℝ) < K := by rw [hKdef]; linarith
  set κ : ℝ := min 1 (ε / K) with hkapdef
  have hκ0 : (0 : ℝ) < κ := lt_min one_pos (div_pos hε hK0)
  have hκ1 : κ ≤ 1 := min_le_left _ _
  have hκK : κ * K ≤ ε := by
    have hle : κ ≤ ε / K := min_le_right _ _
    have := mul_le_mul_of_nonneg_right hle hK0.le
    rwa [div_mul_cancel₀ _ hK0.ne'] at this
  -- the Lipschitz constant `L + 1`
  set LL : ℝ≥0 := ⟨L + 1, by positivity⟩ with hLLdef
  have hLL0 : (0 : ℝ≥0) < LL := by
    rw [← NNReal.coe_pos]; change (0 : ℝ) < L + 1; linarith
  have hLLc : (LL : ℝ) = L + 1 := rfl
  have hSLL : ∀ N, ‖S N‖ ≤ (LL : ℝ) := fun N => by rw [hLLc]; linarith [hSL N]
  -- the tail parameter and the bad event
  set t : ℕ → ℝ := fun N => κ * Real.sqrt (d N) with htdef
  have hsd : ∀ N, (0 : ℝ) < Real.sqrt (d N) := fun N =>
    Real.sqrt_pos.mpr (by exact_mod_cast hd N)
  have ht0 : ∀ N, (0 : ℝ) < t N := fun N => by
    rw [htdef]; exact mul_pos hκ0 (hsd N)
  set Bad : (N : ℕ) → Set (Matrix (Fin (P N)) (Fin (p N)) ℝ) := fun N =>
    {Z | (∫ Z', ‖S N * Z'‖ ∂(gaussianMatrix (P N) (p N))) + t N ≤ ‖S N * Z‖} with hBaddef
  have hmeasBad : ∀ N, MeasurableSet (Bad N) := fun N =>
    measurableSet_le measurable_const (R3het.measurable_opNorm_mul _)
  have hmeasure : ∀ N, μ N {ω | B N ω ∈ (Bad N)ᶜ}
      = (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ) := fun N =>
    (hlaw N).measure_eq (hmeasBad N).compl
  have hcompl : ∀ N, (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ)
      = 1 - (gaussianMatrix (P N) (p N)) (Bad N) := fun N =>
    prob_compl_eq_one_sub (hmeasBad N)
  have hd2 : ∀ᶠ N in atTop, r + 1 ≤ d N := hdtop.eventually (eventually_ge_atTop (r + 1))
  have hp0 : ∀ N, r + 1 ≤ d N → 0 < p N := fun N h => by have := hpd N; omega
  -- the tail bound
  have hsq : ∀ N, (t N) ^ 2 = κ ^ 2 * (d N : ℝ) := fun N => by
    change (κ * Real.sqrt (d N)) ^ 2 = κ ^ 2 * (d N : ℝ)
    rw [mul_pow, Real.sq_sqrt (by positivity)]
  have hbound : ∀ᶠ N in atTop, (gaussianMatrix (P N) (p N)) (Bad N)
      ≤ ENNReal.ofReal (Real.exp (-(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2))) := by
    filter_upwards [hd2] with N hN
    have h := R3het.measure_opNorm_mul_ge_le (hP N) (hp0 N hN) (S N) hLL0 (hSLL N) (ht0 N)
    rw [hLLc, hsq N] at h
    have hfin : (gaussianMatrix (P N) (p N)) (Bad N) ≠ ⊤ := measure_ne_top _ _
    rw [← ENNReal.ofReal_toReal hfin]
    exact ENNReal.ofReal_le_ofReal h
  have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hdtop
  have hexp0 : Tendsto (fun N => -(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)) atTop atBot := by
    have h2 : Tendsto (fun N => (κ ^ 2 / (2 * (L + 1) ^ 2)) * (d N : ℝ)) atTop atTop :=
      Filter.Tendsto.const_mul_atTop (by positivity) hdR
    have heq : ∀ N : ℕ, -(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)
        = -((κ ^ 2 / (2 * (L + 1) ^ 2)) * (d N : ℝ)) := fun N => by ring
    simp only [heq]
    exact tendsto_neg_atTop_atBot.comp h2
  have hzero : Tendsto (fun N => (gaussianMatrix (P N) (p N)) (Bad N)) atTop (𝓝 0) := by
    have hE : Tendsto (fun N => ENNReal.ofReal
        (Real.exp (-(κ ^ 2 * (d N : ℝ)) / (2 * (L + 1) ^ 2)))) atTop (𝓝 0) := by
      have := ENNReal.tendsto_ofReal (Real.tendsto_exp_atBot.comp hexp0)
      simpa using this
    exact tendsto_of_tendsto_of_tendsto_of_le_of_le'
      (tendsto_const_nhds (x := (0 : ℝ≥0∞)) (f := atTop)) hE
      (Eventually.of_forall fun _ => by simp) hbound
  -- the normalized bound is eventually close to `√b`
  set rn : ℕ → ℝ := fun N => bd N / Real.sqrt (d N) with hrndef
  have hrnlim : Tendsto rn atTop (𝓝 sc) := by
    refine hlim.sqrt.congr fun N => ?_
    rw [Real.sqrt_div (sq_nonneg (bd N)), Real.sqrt_sq (hbd0 N)]
  have hrsmall : ∀ᶠ N in atTop, rn N < sc + κ :=
    Filter.Tendsto.eventually_lt_const (by linarith) hrnlim
  -- the eventual inclusion
  have hev : ∀ᶠ N in atTop, (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N)
      ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} := by
    filter_upwards [hrsmall, hd2, hbd] with N hrN hN hbdN
    have hsub : {ω | B N ω ∈ (Bad N)ᶜ}
        ⊆ {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} := by
      intro ω hω
      simp only [hBaddef, Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω
      set sd : ℝ := Real.sqrt (d N) with hsddef
      have hsd0 : (0 : ℝ) < sd := hsd N
      have hr0 : (0 : ℝ) ≤ rn N := div_nonneg (hbd0 N) hsd0.le
      have hsp : rn N * sd = bd N := div_mul_cancel₀ _ hsd0.ne'
      have hM0 : (0 : ℝ) ≤ ‖S N * B N ω‖ := norm_nonneg _
      have hMle : ‖S N * B N ω‖ ≤ (rn N + κ) * sd := by
        have h1 : ‖S N * B N ω‖
            < (∫ Z', ‖S N * Z'‖ ∂(gaussianMatrix (P N) (p N))) + t N := hω
        have h3 : t N = κ * sd := rfl
        have h4 : (rn N + κ) * sd = bd N + κ * sd := by rw [← hsp]; ring
        linarith [h1, hbdN, h3, h4]
      -- the Rayleigh bound
      have hherm : (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)).IsHermitian :=
        (hW₀ N ω) ▸ (hsymm N ω)
      have hlm : lamMax (W₀ N ω) (hsymm N ω)
          = lamMax (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)) hherm :=
        lamMax_congr (hW₀ N ω) _ _
      have hray : lamMax (((d N : ℝ))⁻¹ • ((S N * B N ω) * (S N * B N ω)ᵀ)) hherm
          ≤ ((d N : ℝ))⁻¹ * ‖S N * B N ω‖ ^ 2 :=
        R3het.lamMax_gramT_le_opNorm_sq (by positivity) _ hherm
      -- the arithmetic
      have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hsdsq : sd ^ 2 = (d N : ℝ) := Real.sq_sqrt hdpos.le
      have hrk : (0 : ℝ) ≤ rn N + κ := by linarith
      have hMsq : ‖S N * B N ω‖ ^ 2 ≤ ((rn N + κ) * sd) ^ 2 := pow_le_pow_left₀ hM0 hMle 2
      have hfinal : ((d N : ℝ))⁻¹ * ‖S N * B N ω‖ ^ 2 ≤ (rn N + κ) ^ 2 := by
        rw [inv_mul_le_iff₀ hdpos]
        calc ‖S N * B N ω‖ ^ 2 ≤ ((rn N + κ) * sd) ^ 2 := hMsq
          _ = (rn N + κ) ^ 2 * (d N : ℝ) := by rw [mul_pow, hsdsq]
          _ = (d N : ℝ) * (rn N + κ) ^ 2 := by ring
      have hedge : (rn N + κ) ^ 2 ≤ b + ε := by
        have h := R3het.edge_arith_het (sc := sc) (L := 0) (κ := κ) (ε := ε) (r := rn N)
          hsc0 le_rfl hκ0 hκ1 (by linarith [hκK]) hr0 hrN.le
        simpa only [add_zero, zero_add, hscsq] using h
      simp only [Set.mem_ofPred_eq]
      rw [hlm]
      linarith [hray, hfinal, hedge]
    calc (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N)
        = (gaussianMatrix (P N) (p N)) ((Bad N)ᶜ) := (hcompl N).symm
      _ = μ N {ω | B N ω ∈ (Bad N)ᶜ} := (hmeasure N).symm
      _ ≤ μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} := measure_mono hsub
  -- squeeze
  have hlow : Tendsto (fun N => (1 : ℝ≥0∞) - (gaussianMatrix (P N) (p N)) (Bad N))
      atTop (𝓝 1) := by
    have := ENNReal.Tendsto.sub (tendsto_const_nhds (x := (1 : ℝ≥0∞)) (f := atTop))
      hzero (Or.inl (by simp))
    simpa using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds hev
    (Eventually.of_forall fun N => prob_le_one)

/-! ### E9: the model theorem -/

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

omit [NeZero M] in
/-- **E9, the black box `HeteroEdge` for Gaussian noise.** The exact heteroscedastic upper
edge: for every `ε > 0`, `lamMax W₀' ≤ bHet c w + ε` with probability tending to `1`.

The bound sequence is the sharp Sudakov-Fernique bound of task C
(`integral_opNorm_mul_le_sharp`) at the fixed `γ⋆ = -1/s⋆`, which is admissible because
`γ⋆ > max_i w_i²` (`MPhet.sq_lt_neg_inv_sStar`). Its square over `d N` is exactly
`γ⋆ (d N - 1)/d N + ∑_i (n_i/d) w_i² γ⋆/(γ⋆ - w_i²)`, so the limit is
`MPhet.edgeObjective c w γ⋆ = MPhet.bHet c w`. -/
theorem heteroEdge_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroEdge w c (MPhet.bHet c w) := by
  classical
  set γ : ℝ := -1 / MPhet.sStar c w with hγdef
  have hγ0 : (0 : ℝ) < γ := MPhet.neg_inv_sStar_pos hc hw
  have hwγ : ∀ i, w i ^ 2 < γ := fun i => MPhet.sq_lt_neg_inv_sStar hc hw i
  -- the plumbing of `tendsto_measure_lamMax_W0het_le_bSF`
  have hd : ∀ N, 0 < d N := (m.tbl 0).hd
  have hp : ∀ N, d N = (d N - 1) + 1 := fun N => (Nat.sub_add_cancel (hd N)).symm
  choose B hB hlaw using fun N => m.exists_block_hasLaw_het N hG (hp N)
  have hW : ∀ N ω, m.W0het w N ω = ((d N : ℝ))⁻¹ •
      ((m.SigmaHalf w N * B N ω) * (m.SigmaHalf w N * B N ω)ᵀ) :=
    fun N ω => m.W0het_eq_of_block w N ω (B N ω) (hB N ω)
  have hlawB : ∀ N, HasLaw (B N) (gaussianMatrix (∑ i, n i N) (d N - 1)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hlaw N)
  -- the stacked diagonal of the weights
  set a : (N : ℕ) → Fin (∑ i, n i N) → ℝ := fun N r => w (finSigmaFinEquiv.symm r).1 with hadef
  have hSig : ∀ N, m.SigmaHalf w N = Matrix.diagonal (a N) := fun N => rfl
  have haγ : ∀ (N : ℕ) (r : Fin (∑ i, n i N)), (a N r) ^ 2 < γ := fun N r => hwγ _
  -- the sharp bound sequence and its square
  set bd : ℕ → ℝ := fun N => Real.sqrt (γ * ((d N - 1 : ℕ) : ℝ)
    + ∑ r : Fin (∑ i, n i N), (a N r) ^ 2 * γ / (γ - (a N r) ^ 2)) with hbddef
  have hrad : ∀ N, (0 : ℝ) ≤ γ * ((d N - 1 : ℕ) : ℝ)
      + ∑ r : Fin (∑ i, n i N), (a N r) ^ 2 * γ / (γ - (a N r) ^ 2) := by
    intro N
    refine add_nonneg (by positivity) (Finset.sum_nonneg fun r _ => ?_)
    exact div_nonneg (mul_nonneg (sq_nonneg _) hγ0.le) (sub_pos.mpr (haγ N r)).le
  have hbd0 : ∀ N, (0 : ℝ) ≤ bd N := fun N => Real.sqrt_nonneg _
  have hbdsq : ∀ N, bd N ^ 2 = γ * ((d N - 1 : ℕ) : ℝ)
      + ∑ r : Fin (∑ i, n i N), (a N r) ^ 2 * γ / (γ - (a N r) ^ 2) := by
    intro N
    simp only [hbddef]
    exact Real.sq_sqrt (hrad N)
  -- the expectation bound of task C, once `p N = d N - 1` is positive
  have hd2 : ∀ᶠ N in atTop, 2 ≤ d N := (hreg 0).2.1.eventually (eventually_ge_atTop 2)
  have hbd : ∀ᶠ N in atTop, (∫ B', ‖m.SigmaHalf w N * B'‖
      ∂(gaussianMatrix (∑ i, n i N) (d N - 1))) ≤ bd N := by
    filter_upwards [hd2] with N hN
    have hpp : 0 < d N - 1 := by omega
    have h := integral_opNorm_mul_le_sharp (m.stack_row_pos N) hpp (a N) hγ0 (haγ N)
    rw [hSig N]
    simpa only [hbddef] using h
  -- regrouping the diagonal sum by table
  have hregroup : ∀ N, ∑ r : Fin (∑ i, n i N), (a N r) ^ 2 * γ / (γ - (a N r) ^ 2)
      = ∑ i, (n i N : ℝ) * (w i ^ 2 * γ / (γ - w i ^ 2)) := by
    intro N
    simp only [hadef]
    rw [← Equiv.sum_comp finSigmaFinEquiv]
    simp only [Equiv.symm_apply_apply]
    rw [Fintype.sum_sigma]
    simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  -- the limit of the normalized square: an exact identity, then `(d-1)/d → 1` and `n_i/d → c_i`
  have hval : γ + ∑ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2)) = MPhet.bHet c w := by
    rw [← MPhet.edgeObjective_sStar_eq_bHet hc hw, ← hγdef]
    have hterm : ∀ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2))
        = γ * (c i * w i ^ 2 / (γ - w i ^ 2)) := fun i => by ring
    simp only [hterm, MPhet.edgeObjective]
    rw [mul_add, mul_one, Finset.mul_sum]
  have hone : Tendsto (fun N => (((d N - 1 : ℕ) : ℝ)) / (d N : ℝ)) atTop (𝓝 1) := by
    have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp (hreg 0).2.1
    have h : Tendsto (fun N => 1 - ((d N : ℝ))⁻¹) atTop (𝓝 (1 - 0)) :=
      tendsto_const_nhds.sub (tendsto_inv_atTop_zero.comp hdR)
    rw [sub_zero] at h
    refine h.congr fun N => ?_
    have hdN : ((d N : ℝ)) ≠ 0 := by
      have : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      exact this.ne'
    rw [Nat.cast_sub (hd N), Nat.cast_one]
    field_simp
  have hlim : Tendsto (fun N => bd N ^ 2 / (d N : ℝ)) atTop (𝓝 (MPhet.bHet c w)) := by
    rw [← hval]
    have h1 : Tendsto (fun N => γ * ((((d N - 1 : ℕ) : ℝ)) / (d N : ℝ))) atTop (𝓝 γ) := by
      simpa using hone.const_mul γ
    have h2 : Tendsto (fun N => ∑ i, ((n i N : ℝ) / (d N : ℝ)) * (w i ^ 2 * γ / (γ - w i ^ 2)))
        atTop (𝓝 (∑ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2)))) :=
      tendsto_finsetSum _ fun i _ => (hreg i).2.2.mul_const _
    refine (h1.add h2).congr fun N => ?_
    rw [hbdsq N, hregroup N, add_div, Finset.sum_div]
    congr 1
    · ring
    · exact Finset.sum_congr rfl fun i _ => by ring
  refine ⟨fun ε hε => ?_⟩
  exact tendsto_measure_lamMax_le_of_bound μ (fun N => m.SigmaHalf w N)
    (Real.sqrt_nonneg _) (MPhet.bHet_pos hc hw).le (fun N => m.opNorm_SigmaHalf_le w N)
    bd hbd0 hbd hlim B (m.W0het w) hW (m.isHermitian_W0het w) hlawB m.stack_row_pos hd
    (fun N => Nat.sub_add_cancel (hd N)) (hreg 0).2.1 hε

end MultiTableModel
end StackedSVD
