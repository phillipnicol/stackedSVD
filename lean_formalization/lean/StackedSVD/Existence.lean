/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Sat
import StackedSVD.SVDStack.Inad
import StackedSVD.RemarksUniform
import StackedSVD.StackSVD.Main

/-!
# The paper's existential claims, with the model built rather than assumed

Three statements of the paper are existential: `prop:binarystacksvd_inadmissable`
(`main_paper.tex:660`, "for any `ε ∈ (0,1)`, there exists a problem instance"), and the two
remarks `remark:stack_outperform_svd` (`main_paper.tex:566`, unweighted stacksvd *can*
outperform optimally weighted svdstack) and `remark:svd_outperform_stack`
(`main_paper.tex:588`, unweighted svdstack *can* outperform binary-weighted stacksvd). The
facades of `RMT/Het/Sup.lean`, `SVDStack/Inad.lean`, `Remarks.lean` and `RemarksUniform.lean`
take a model `m` with the instance's `θ`, regime and Gaussian noise as hypotheses. This file
supplies the model: the concrete Gaussian construction `Sat.Rank1.model` of `Sat.lean`, with
`n_i N = k_i (N+1)` and `d N = l (N+1)`, so that `c_i = k_i / l` exactly. Each theorem below
therefore concludes `∃ ... (m : MultiTableModel μ M n d), ...` with no hypothesis on a model
(external rank-one audit of 2026-09-03, findings 1 and 2).

* `MultiTableModel.prop_binarystacksvd_inadmissable_gaussian` bundles the five limits of the
  proposition on any model of the instance `θ_i = 1`, `c_i = 2i + 1`.
* `Sat.Rank1.inadModel M` is that instance as a Gaussian model (`k_i = 2i + 1`, `l = 1`), and
  `prop_binarystacksvd_inadmissable_exists` is the paper's statement: for every `ε ∈ (0,1)`
  there is a size `M = ⌈e^{-γ} e^{2/ε}⌉` and a model of `M` tables on which optimally weighted
  stacksvd tends to a value above `1 - ε` while every binary weighting of stacksvd, optimally
  weighted svdstack, unweighted stacksvd and unweighted svdstack tend to `0`.
* `remark_stack_outperform_svd_exists`: for every `M ≥ 2` there is a model of `M` tables
  (`θ_i = c_i = 1`) on which unweighted stacksvd tends to `1 - 2/(M+1) > 0` while unweighted
  svdstack tends to `0` and, with probability tending to one, no weighting of svdstack reaches
  any `ε > 0`.
* `remark_svd_outperform_stack_exists`: there is a model of two tables (`θ = (√5, 4)`,
  `c = (1, 38.4)`, so `n_1 N = 5(N+1)`, `n_2 N = 192(N+1)`, `d N = 5(N+1)`) on which every
  binary weighting of stacksvd tends to a value at most `2008/2310` while optimally weighted
  and unweighted svdstack tend to `8/9 > 2008/2310`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- `prop:binarystacksvd_inadmissable` (`main_paper.tex:660`) on a model of the instance
`θ_i = 1`, `c_i = 2i + 1` with Gaussian noise and `M ≥ e^{-γ} e^{2/ε}`: the five limits of
the proposition in one statement. In order: optimally weighted stacksvd tends to its limit,
which is above `1 - ε`; every binary weighting of stacksvd (every nonempty `S`) tends to `0`;
optimally weighted svdstack tends to `0`; unweighted stacksvd tends to `0`; unweighted
svdstack tends to `0`. The pieces are `stackPerfW_opt_inad_gaussian`,
`stackPerfW_binary_inad_gaussian` (`RMT/Het/Sup.lean`), `svdstackPerfW_opt_inad_gaussian`,
`svdstackPerf_inad_gaussian` (`SVDStack/Inad.lean`) and `prop_stacksvd_general_gaussian` with
`Scalars.inad_stackSVDLimit_eq_zero`. -/
theorem prop_binarystacksvd_inadmissable_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)] (m : MultiTableModel μ M n d)
    (hθ : ∀ i, (m.tbl i).θ = 1) {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    (TendstoInProb μ
        (fun N ω => m.stackPerfW
          (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
        (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
    (∀ S : Finset (Fin M), S.Nonempty →
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0) ∧
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
      0 ∧
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N)) 0 ∧
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 := by
  refine ⟨m.stackPerfW_opt_inad_gaussian hθ hε hε1 hM hreg hG, fun S hS => ?_,
    m.svdstackPerfW_opt_inad_gaussian hθ hreg hG, ?_, m.svdstackPerf_inad_gaussian hθ hreg hG⟩
  · have : NeZero S.card := ⟨(Finset.card_pos.mpr hS).ne'⟩
    exact m.stackPerfW_binary_inad_gaussian hθ S hreg hG
  · have hfun : (fun i => (m.tbl i).θ) = Scalars.inadTheta M := funext fun i => hθ i
    have h := m.prop_stacksvd_general_gaussian (Scalars.inadC M) (Scalars.inadC_pos M) hreg hG
    rw [hfun, Scalars.inad_stackSVDLimit_eq_zero M] at h
    exact h

end MultiTableModel

namespace Sat.Rank1

/-! ## `prop:binarystacksvd_inadmissable`: the instance as a Gaussian model -/

/-- The row-count factors of the instance: `k_i = 2i + 1`, so that with `l = 1` the aspect
ratio is `c_i = 2i + 1 = Scalars.inadC M i`. -/
abbrev inadK (M : ℕ) : Fin M → ℕ := fun i => 2 * i.val + 1

theorem inadK_pos (M : ℕ) : ∀ i, 0 < inadK M i := fun _ => Nat.succ_pos _

/-- The instance of `prop:binarystacksvd_inadmissable` as a Gaussian model: `M` tables,
`n_i N = (2i + 1)(N + 1)`, `d N = N + 1`, `θ_i = 1`. -/
noncomputable def inadModel (M : ℕ) :
    MultiTableModel (mu M (nk (inadK M)) (dl 1)) M (nk (inadK M)) (dl 1) :=
  model (inadK M) 1 (inadK_pos M) one_pos (fun _ => 1) (fun _ => zero_le_one)

theorem inadModel_θ (M : ℕ) (i : Fin M) : ((inadModel M).tbl i).θ = 1 :=
  model_θ (inadK M) 1 (inadK_pos M) one_pos (fun _ => 1) (fun _ => zero_le_one) i

theorem inadModel_regime (M : ℕ) (i : Fin M) :
    ((inadModel M).tbl i).Regime (Scalars.inadC M i) :=
  model_regime (inadK M) 1 (inadK_pos M) one_pos (fun _ => 1) (fun _ => zero_le_one)
    (Scalars.inadC M) (fun i => by simp [inadK, Scalars.inadC]) i

theorem inadModel_joint (M : ℕ) : (inadModel M).JointGaussianNoise :=
  model_joint (inadK M) 1 (inadK_pos M) one_pos (fun _ => 1) (fun _ => zero_le_one)

/-- **`prop:binarystacksvd_inadmissable`** (`main_paper.tex:660`) as the paper states it: for
every `ε ∈ (0,1)` there is a size `M = ⌈e^{-γ} e^{2/ε}⌉` and a problem instance of `M` tables
(a Gaussian model, `inadModel M`, with `θ_i = 1` and `c_i = 2i + 1`) on which optimally
weighted stacksvd tends in probability to a value above `1 - ε`, while every binary weighting
of stacksvd, optimally weighted svdstack, unweighted stacksvd and unweighted svdstack tend to
`0`. No hypothesis on a model remains: the model is built. The witness `NeZero M` (that is,
`M ≠ 0`) is part of the existential because the estimators are defined for `M ≠ 0`. -/
theorem prop_binarystacksvd_inadmissable_exists {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) :
    ∃ (M : ℕ) (_ : NeZero M) (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N))
      (μ : ∀ N, Measure (Ω N)) (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin M → ℕ → ℕ)
      (d : ℕ → ℕ) (m : MultiTableModel μ M n d),
      M = ⌈Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε)⌉₊ ∧
      (∀ i, (m.tbl i).θ = 1) ∧ (∀ i, (m.tbl i).Regime (Scalars.inadC M i)) ∧
      m.JointGaussianNoise ∧
      (TendstoInProb μ
          (fun N ω => m.stackPerfW
            (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
          (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
        1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      (∀ S : Finset (Fin M), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0) ∧
      TendstoInProb μ
        (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
        0 ∧
      TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N)) 0 ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 := by
  set M := ⌈Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε)⌉₊ with hMdef
  have hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M := Nat.le_ceil _
  have hpos : 0 < M := by
    have h1 : (0 : ℝ) < Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) :=
      mul_pos (Real.exp_pos _) (Real.exp_pos _)
    exact Nat.ceil_pos.mpr h1
  have hne : NeZero M := ⟨hpos.ne'⟩
  exact ⟨M, hne, Om M (nk (inadK M)) (dl 1), inferInstance, mu M (nk (inadK M)) (dl 1),
    inferInstance, nk (inadK M), dl 1, inadModel M, rfl, inadModel_θ M, inadModel_regime M,
    inadModel_joint M,
    (inadModel M).prop_binarystacksvd_inadmissable_gaussian (inadModel_θ M) hε hε1 hM
      (inadModel_regime M) (inadModel_joint M)⟩

/-! ## `remark:stack_outperform_svd`: `θ_i = c_i = 1` as a Gaussian model -/

/-- The instance of `remark:stack_outperform_svd` as a Gaussian model: `M` tables,
`n_i N = d N = N + 1`, `θ_i = 1`. -/
noncomputable def oneModel (M : ℕ) :
    MultiTableModel (mu M (nk fun _ : Fin M => 1) (dl 1)) M (nk fun _ => 1) (dl 1) :=
  model (fun _ => 1) 1 (fun _ => one_pos) one_pos (fun _ => 1) (fun _ => zero_le_one)

theorem oneModel_θ (M : ℕ) (i : Fin M) : ((oneModel M).tbl i).θ = 1 :=
  model_θ (fun _ => 1) 1 (fun _ => one_pos) one_pos (fun _ => 1) (fun _ => zero_le_one) i

theorem oneModel_regime (M : ℕ) (i : Fin M) : ((oneModel M).tbl i).Regime 1 :=
  model_regime (fun _ => 1) 1 (fun _ => one_pos) one_pos (fun _ => 1) (fun _ => zero_le_one)
    (fun _ => 1) (fun _ => by simp) i

theorem oneModel_joint (M : ℕ) : (oneModel M).JointGaussianNoise :=
  model_joint (fun _ => 1) 1 (fun _ => one_pos) one_pos (fun _ => 1) (fun _ => zero_le_one)

/-- **`remark:stack_outperform_svd`** (`main_paper.tex:566`) with the model built: for every
`M ≥ 2` there is a problem instance of `M` tables (a Gaussian model, `oneModel M`, with
`θ_i = c_i = 1`) on which unweighted stacksvd tends in probability to `1 - 2/(M+1) > 0`,
unweighted svdstack tends to `0`, and with probability tending to one no nonzero weighting of
svdstack, data dependent or not, reaches any `ε > 0`. So unweighted stacksvd outperforms
optimally weighted svdstack there. The pieces are `remark_stack_outperform_svd_stack`,
`remark_stack_outperform_svd_svdstack` (`Remarks.lean`) and
`remark_stack_outperform_svd_svdstack_uniform` (`RemarksUniform.lean`). The instance
argument `[NeZero M]` follows from `2 ≤ M`; it is separate because the estimators in the
statement are defined for `M ≠ 0`. -/
theorem remark_stack_outperform_svd_exists (M : ℕ) [NeZero M] (hM : 2 ≤ M) :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ)
      (m : MultiTableModel μ M n d),
      (∀ i, (m.tbl i).θ = 1) ∧ (∀ i, (m.tbl i).Regime 1) ∧ m.JointGaussianNoise ∧
      TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
        (1 - 2 / ((M : ℝ) + 1)) ∧
      0 < 1 - 2 / ((M : ℝ) + 1) ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 ∧
      (∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
        ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0)) := by
  have hMr : (2 : ℝ) ≤ M := by exact_mod_cast hM
  refine ⟨Om M (nk fun _ => 1) (dl 1), inferInstance, mu M (nk fun _ => 1) (dl 1),
    inferInstance, nk fun _ => 1, dl 1, oneModel M, oneModel_θ M, oneModel_regime M,
    oneModel_joint M,
    (oneModel M).remark_stack_outperform_svd_stack hM (oneModel_θ M) (oneModel_regime M)
      (oneModel_joint M), ?_,
    (oneModel M).remark_stack_outperform_svd_svdstack (oneModel_θ M) (oneModel_regime M)
      (oneModel_joint M),
    (oneModel M).remark_stack_outperform_svd_svdstack_uniform (oneModel_θ M)
      (oneModel_regime M) (oneModel_joint M)⟩
  have h3 : (0 : ℝ) < (M : ℝ) + 1 := by linarith
  rw [sub_pos, div_lt_one h3]
  linarith

/-! ## `remark:svd_outperform_stack`: `θ = (√5, 4)`, `c = (1, 38.4)` as a Gaussian model -/

/-- The row-count factors of the paper's two-table example: `k = (5, 192)` with `l = 5`, so
`c = (1, 192/5) = (1, 38.4)`. -/
abbrev twoK : Fin 2 → ℕ := ![5, 192]

theorem twoK_pos : ∀ i, 0 < twoK i := by
  intro i; fin_cases i <;> simp [twoK]

theorem twoTheta_nonneg : ∀ i, 0 ≤ (![Real.sqrt 5, 4] : Fin 2 → ℝ) i := by
  intro i; fin_cases i <;> simp

/-- The paper's two-table example as a Gaussian model: `n_1 N = 5(N+1)`, `n_2 N = 192(N+1)`,
`d N = 5(N+1)`, `θ = (√5, 4)`. -/
noncomputable def twoModel : MultiTableModel (mu 2 (nk twoK) (dl 5)) 2 (nk twoK) (dl 5) :=
  model twoK 5 twoK_pos (by norm_num) ![Real.sqrt 5, 4] twoTheta_nonneg

theorem twoModel_θ (i : Fin 2) : (twoModel.tbl i).θ = ![Real.sqrt 5, 4] i :=
  model_θ twoK 5 twoK_pos (by norm_num) ![Real.sqrt 5, 4] twoTheta_nonneg i

theorem twoModel_regime (i : Fin 2) : (twoModel.tbl i).Regime (![1, 38.4] i) :=
  model_regime twoK 5 twoK_pos (by norm_num) ![Real.sqrt 5, 4] twoTheta_nonneg ![1, 38.4]
    (fun i => by fin_cases i <;> norm_num [twoK]) i

theorem twoModel_joint : twoModel.JointGaussianNoise :=
  model_joint twoK 5 twoK_pos (by norm_num) ![Real.sqrt 5, 4] twoTheta_nonneg

/-- **`remark:svd_outperform_stack`** (`main_paper.tex:588`) with the model built: there is a
problem instance of two tables (a Gaussian model, `twoModel`, with `θ = (√5, 4)` and
`c = (1, 38.4)`, so `β_1² = β_2² = 4/5`) on which every binary weighting of stacksvd (every
nonempty `S`, the unweighted stack included at `S = univ`) tends in probability to a value at
most `2008/2310`, while optimally weighted svdstack and unweighted svdstack both tend to
`8/9 > 2008/2310`. So unweighted svdstack outperforms every binary weighting of stacksvd
there. The pieces are `remark_svd_outperform_stack_two` (`Remarks.lean`) and
`remark_svd_outperform_stack_two_svdstack_unweighted` (`RemarksUniform.lean`). -/
theorem remark_svd_outperform_stack_exists :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : MultiTableModel μ 2 n d),
      (∀ i, (m.tbl i).θ = ![Real.sqrt 5, 4] i) ∧ (∀ i, (m.tbl i).Regime (![1, 38.4] i)) ∧
      m.JointGaussianNoise ∧
      (∀ i, beta ((m.tbl i).θ) (![1, 38.4] i) ^ 2 = 4 / 5) ∧
      (∀ S : Finset (Fin 2), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
          (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4]) ∧
        Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4] ≤ 2008 / 2310) ∧
      TendstoInProb μ (fun N ω =>
        m.svdstackPerfW (optW fun i => beta ((m.tbl i).θ) (![1, 38.4] i)) N ω) (8 / 9) ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (8 / 9) ∧
      (2008 : ℝ) / 2310 < 8 / 9 := by
  obtain ⟨hβ, -, hbin, hopt, hlt⟩ :=
    twoModel.remark_svd_outperform_stack_two twoModel_θ twoModel_regime twoModel_joint
  exact ⟨Om 2 (nk twoK) (dl 5), inferInstance, mu 2 (nk twoK) (dl 5), inferInstance,
    nk twoK, dl 5, twoModel, twoModel_θ, twoModel_regime, twoModel_joint, hβ, hbin, hopt,
    twoModel.remark_svd_outperform_stack_two_svdstack_unweighted twoModel_θ twoModel_regime
      twoModel_joint, hlt⟩

end Sat.Rank1

end StackedSVD
