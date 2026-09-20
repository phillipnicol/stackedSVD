/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Weighted

/-!
# `prop:binarystacksvd_inadmissable`: the svdstack clauses on the model

The proposition says that on the instance `θ_i = 1`, `c_i = 2i - 1` (0-based: `2i + 1`),
optimally weighted svdstack and unweighted svdstack are below their recovery thresholds.
`Scalars.binarystacksvd_inadmissable` proves the scalar clauses `svdstackLimitOpt β = 0` and
`svdstackLimit β = 0`; `RMT/Het/Sup.lean` proves the stacksvd clauses on the model with
Gaussian noise. This file adds the svdstack clauses on the model: under independent Gaussian
noise the performance of unweighted svdstack, and of svdstack at the paper's optimal weights
`optW β`, tends to `0` in probability (external statement audit of 2026-09-02, finding 4,
item L3).

Route. On the instance every `β_i = beta 1 (2i + 1) = 0`, because `θ_i⁴ = 1 ≤ c_i`. The
all-subcritical theorems `thm_svd_stack_general_zero_gaussian` (`SVDStack/Main.lean`) and
`thm_svdstack_weighted_zero_gaussian` (`SVDStack/Weighted.lean`) then apply; the weighted one
needs a nonzero weight, and `optW β 0 = 1`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-- On the instance every `β_i` is `0`: `θ_i⁴ = 1 ≤ c_i = 2i + 1`. -/
theorem beta_one_inadC_eq_zero (M : ℕ) (i : Fin M) : beta 1 (Scalars.inadC M i) = 0 := by
  have hi : (0:ℝ) ≤ (i.val : ℝ) := Nat.cast_nonneg _
  have hle : (1:ℝ) ^ 4 ≤ Scalars.inadC M i := by
    simp only [Scalars.inadC, one_pow]
    linarith
  unfold beta betaSq
  rw [if_neg (not_lt.mpr hle), Real.sqrt_zero]

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- `prop:binarystacksvd_inadmissable`, unweighted svdstack on the model: on the instance,
under independent Gaussian noise, the performance of svdstack tends to `0`. -/
theorem svdstackPerf_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 :=
  m.thm_svd_stack_general_zero_gaussian (Scalars.inadC M)
    (fun i => beta (m.tbl i).θ (Scalars.inadC M i)) (Scalars.inadC_pos M) (fun _ => rfl)
    (fun i => by rw [hθ i]; exact beta_one_inadC_eq_zero M i) hreg hG

/-- `prop:binarystacksvd_inadmissable`, optimally weighted svdstack on the model: on the
instance, under independent Gaussian noise, the performance of svdstack at the paper's optimal
weights `optW β` tends to `0`. On the instance `β = 0` and `optW β = 1`, so the weight vector
is nonzero, which `thm_svdstack_weighted_zero_gaussian` needs. -/
theorem svdstackPerfW_opt_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
      0 := by
  have hβ0 : ∀ i, beta (m.tbl i).θ (Scalars.inadC M i) = 0 := fun i => by
    rw [hθ i]
    exact beta_one_inadC_eq_zero M i
  refine m.thm_svdstack_weighted_zero_gaussian (Scalars.inadC M) _ _ (Scalars.inadC_pos M)
    (fun _ => rfl) hβ0 ⟨0, ?_⟩ hreg hG
  simp [optW, hβ0 0]

end MultiTableModel

end StackedSVD
