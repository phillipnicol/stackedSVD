/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Remarks
import StackedSVD.SVDStack.Rayleigh

/-!
# Two remark instances that need `SVDStack/Rayleigh.lean`

`StackedSVD/Remarks.lean` states the explicit instances of `remark:stack_outperform_svd` and
`remark:svd_outperform_stack` (`main_paper.tex:566` and `main_paper.tex:588`). Two of them do
not fit in that file, and they live here.

## Content

1. `MultiTableModel.remark_stack_outperform_svd_svdstack_uniform`, instance (i)(c). At
   `θ_i = c_i = 1` every table sits at its own detection threshold, so `β = 0`. The paper says
   that optimally weighted SVDstack then has performance `0`. `thm:svdstack_weighted` needs
   some `β_k > 0` and does not apply at `β = 0`, so the statement is P1's uniform bound
   (`svdstackPerfW_uniform_bound`, `SVDStack/Rayleigh.lean`) at `S = 0`: with probability
   tending to one no nonzero weight vector, data dependent or not, reaches `ε`. That covers
   the optimal weighting and every other one. `Remarks.lean` cannot state it because it does
   not import `SVDStack/Rayleigh.lean` (choice 1 of `notes/archive/remark_facades.md`).
2. `MultiTableModel.remark_svd_outperform_stack_two_svdstack_unweighted`, the missing
   unweighted value of instance (iv) (D33 audit, finding F5). The paper's `M = 2` example has
   `θ = (√5, 4)` and `c = (1, 38.4)`, so both tables have `β² = 4/5`.
   `remark_svd_outperform_stack_two` (`Remarks.lean`) records the optimally weighted SVDstack
   value `8/9`; the unweighted SVDstack reaches the same `8/9`, because equal `β_i` make the
   optimal weights constant. This is choice 4 of `notes/archive/remark_facades.md`, which
   marks the value as a conjunct that can be added.

`svdstackLimitOpt_zero` is the scalar step of item 1: `S = ∑_i β_i²/(1 − β_i²) = 0` at `β = 0`,
so `S/(S+1) = 0`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix

namespace StackedSVD

/-- `S/(S+1) = 0` at `β = 0`: every term of `Sval` is `0`. -/
theorem svdstackLimitOpt_zero {M : ℕ} : svdstackLimitOpt (fun _ : Fin M => (0 : ℝ)) = 0 := by
  simp [svdstackLimitOpt, Sval]

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [∀ N, IsProbabilityMeasure (μ N)]

/-! ### `remark:stack_outperform_svd`, instance (i)(c): every weighting has performance `0` -/

/-- **`remark:stack_outperform_svd`** (`main_paper.tex:566`), instance (i)(c): at
`θ_i = c_i = 1` every table sits at its own detection threshold `θ_i⁴ = c_i`, so `β = 0` and
`S = 0`. With probability tending to one no nonzero weight vector `w`, data dependent or not,
gives weighted SVDstack a performance of `ε` or more. The optimal weighting is one such `w`,
so this is the paper's "optimally weighted SVDstack falls below the detection threshold",
stated for every weighting at once. Instance of `svdstackPerfW_uniform_bound` at `β = 0`. -/
theorem remark_stack_outperform_svd_svdstack_uniform [NeZero M] (m : MultiTableModel μ M n d)
    (hθ : ∀ i, (m.tbl i).θ = 1) (hreg : ∀ i, (m.tbl i).Regime 1) (hG : m.JointGaussianNoise) :
    ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
      ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0) := by
  have hβdef : ∀ i, (fun _ : Fin M => (0 : ℝ)) i
      = beta (m.tbl i).θ ((fun _ : Fin M => (1 : ℝ)) i) := by
    intro i
    change (0 : ℝ) = beta (m.tbl i).θ 1
    rw [hθ i, beta_eq_zero_of_not_thr (by norm_num : ¬ (1 : ℝ) < 1 ^ 4)]
  have hlaw : ∀ i, (m.tbl i).SingleTableLaw ((fun _ : Fin M => (1 : ℝ)) i) := fun i =>
    SpikedModel.singleTableLaw_of_gaussian one_pos (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)
  have h := m.svdstackPerfW_uniform_bound (fun _ => 1) (fun _ => 0) (fun _ => one_pos)
    hβdef hlaw hG.indepNoise
  simp only [svdstackLimitOpt_zero, zero_add] at h
  exact h

/-! ### `remark:svd_outperform_stack`, instance (iv): the unweighted SVDstack value -/

section Two

variable {n₂ : Fin 2 → ℕ → ℕ}

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:610`), the unweighted SVDstack value of
the paper's two-table example, which `remark_svd_outperform_stack_two` of `Remarks.lean` does
not record (D33 audit, finding F5; choice 4 of `notes/archive/remark_facades.md`). With
`θ = (√5, 4)` and `c = (1, 38.4)` both tables have `β² = 4/5`, so `A_β` has equal entries and
`svdstackLimit_const` gives `2(4/5)/(1 + 4/5) = 8/9`. Optimally weighted SVDstack reaches the
same `8/9`, because equal `β_i` make the optimal weights constant. -/
theorem remark_svd_outperform_stack_two_svdstack_unweighted (m : MultiTableModel μ 2 n₂ d)
    (hθ : ∀ i, (m.tbl i).θ = ![Real.sqrt 5, 4] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 38.4] i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (8 / 9) := by
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
  have hb0 : beta (Real.sqrt 5) 1 ^ 2 = 4 / 5 := by
    rw [Scalars.beta_sq, betaSq,
      if_pos (by rw [h5']; norm_num : Real.sqrt 5 ^ 4 > (1 : ℝ)), h5, h5']
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
  have hbnn : ∀ i : Fin 2, 0 ≤ beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i) := fun i =>
    (beta_mem_Ico (hcpos i)).1
  have hβconst : (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i))
      = fun _ : Fin 2 => Real.sqrt (4 / 5) := by
    funext i
    rw [← hbsq i, Real.sqrt_sq (hbnn i)]
  have hthr : ∃ i j : Fin 2, i ≠ j
      ∧ 0 < beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)
      ∧ 0 < beta ((m.tbl j).θ) (![(1 : ℝ), 38.4] j) := by
    refine ⟨0, 1, by decide, ?_, ?_⟩
    · have h1 := hbsq 0
      have h2 := hbnn 0
      nlinarith
    · have h1 := hbsq 1
      have h2 := hbnn 1
      nlinarith
  have h := m.thm_svd_stack_general_gaussian ![(1 : ℝ), 38.4]
    (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)) hcpos (fun i => rfl) hthr hreg hG
  have hs45 : Real.sqrt (4 / 5 : ℝ) ^ 2 = 4 / 5 := Real.sq_sqrt (by norm_num)
  have hval : svdstackLimit (fun i => beta ((m.tbl i).θ) (![(1 : ℝ), 38.4] i)) = 8 / 9 := by
    rw [hβconst, Scalars.svdstackLimit_const (by norm_num) (Real.sqrt (4 / 5)), hs45]
    norm_num
  rwa [hval] at h

end Two

end MultiTableModel

end StackedSVD
