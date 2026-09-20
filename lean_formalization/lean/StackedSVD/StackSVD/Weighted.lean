/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.StackSVDWeighted
import StackedSVD.StackSVD.Main

/-!
# `thm:stacksvd_weighted` and `cor.2`: the weighted stacksvd limit

Moved out of `StackSVDWeighted.lean` (F28, 2026-09-08): `thm_stacksvd_weighted` and
its Gaussian corollaries, the constant weighting, and `cor.2` (the binary weighting)
through its Gaussian discharge. No proof changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section StackW

variable [NeZero M]

/-! ### `thm:stacksvd_weighted` -/

/-- General weights: under the heteroscedastic law the weighted stacksvd performance tends to
the paper's `L(w)`. This is `law.align`, named so that no consumer unfolds the structure. -/
theorem thm_stacksvd_weighted_general (m : MultiTableModel μ M n d) (w c : Fin M → ℝ)
    (law : m.HeteroLaw w c) :
    TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
      (Scalars.Lw (fun i => (m.tbl i).θ) c w) :=
  law.align

omit [NeZero M] in
set_option linter.unusedVariables false in
/-- `thm:stacksvd_weighted`: at `w_i⋆ = θ_i/√(θ_i²+c_i)` the weighted stacksvd performance
tends to the unique root in `(0,1)` of `∑_i θ_i⁴(1-x)/(c_i + xθ_i²) = 1`, and to `0` below
the detectability threshold `∑_i θ_i⁴/c_i > 1`. Proof: `law.align` rewritten by
`Scalars.L_optW_eq`. `Scalars.L_le_opt` says that no other weighting does better.

`hθ` records that the statement is empty when every `θ_i = 0`: there `optWstack` is the zero
vector, the weighted stack is the zero matrix and `HeteroLaw.topSimple` fails for `2 ≤ d N`
(audit item 1). It is not used by the proof.

**On the measure.** Like `prop_stacksvd_general`, this takes no `IsProbabilityMeasure`
instance and so holds for every measure family, the zero measure included, where it says
nothing (mechanical audit 2026-08-31, finding 3). The hypothesis `law : m.HeteroLaw ..` is
what a caller must produce, and every Gaussian corollary pins the measure through
`GaussianNoise`. -/
theorem thm_stacksvd_weighted [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hθ : ∃ i, (m.tbl i).θ ≠ 0)
    (law : m.HeteroLaw (Scalars.optWstack (fun i => (m.tbl i).θ) c) c) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  have h := law.align
  rwa [Scalars.L_optW_eq hc] at h

set_option linter.unusedVariables false in
/-- The paper's inner-product form of `thm:stacksvd_weighted`, for any selected unit top
eigenvector of the Gram matrix. Consumes `HeteroLaw.topSimple` through
`overlap_eq_inner_sq`. -/
theorem thm_stacksvd_weighted_inner (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hθ : ∃ i, (m.tbl i).θ ≠ 0)
    (law : m.HeteroLaw (Scalars.optWstack (fun i => (m.tbl i).θ) c) c)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈
      topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
        (m.isHermitian_stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) := by
  intro ε hε
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.thm_stacksvd_weighted c hc hθ law ε hε) (fun _ => zero_le) (fun N => ?_)
  refine measure_mono_ae ?_
  filter_upwards [law.topSimple N] with ω hω hmem'
  have heq : m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω
      = ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2 :=
    overlap_eq_inner_sq _ _ hω (hmem N ω) (hnorm N ω)
  have hfin : ε ≤ |m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω
      - Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c| := by
    rw [heq]
    exact hmem'
  exact hfin

/-- The reduction to unit weights: at `w = 1` the heteroscedastic law follows from the
`SingleTableLaw` of the stack, because `stackW 1 = stack` and `L(1) = stackSVDLimit θ c`.
This ties the new structure to the proved `prop_stacksvd_general`, and it is the sanity check
that `Lw` is normalized as the paper's `L`. -/
theorem heteroLaw_one_of_singleTableLaw (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : m.stack.SingleTableLaw (∑ i, c i)) : m.HeteroLaw 1 c := by
  have hst : m.stackW 1 = m.stack := m.stackW_one
  refine ⟨?_, ?_⟩
  · have h := m.prop_stacksvd_general c law
    rw [Scalars.Lw_one]
    simpa [stackPerfW, hst] using h
  · intro N
    filter_upwards [law.topSimple N] with ω hω
    simpa [stackGramW, hst] using hω

/-- The constant weighting: `HeteroLaw m (fun _ => t) c` follows from the `SingleTableLaw` of
the unweighted stack, for every `t ≠ 0`. -/
theorem heteroLaw_const_of_singleTableLaw (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    {t : ℝ} (ht : t ≠ 0) (law : m.stack.SingleTableLaw (∑ i, c i)) :
    m.HeteroLaw (fun _ => t) c := by
  have hX := m.stackW_const_X t
  refine ⟨?_, ?_⟩
  · have h := m.prop_stacksvd_general c law
    rw [Scalars.Lw_const ht]
    have hfun : ∀ (N : ℕ) (ω : Ω N),
        m.stackPerfW (fun _ => t) N ω = overlap (m.stack.X N ω) (m.stack.v N) := by
      intro N ω
      change overlap ((m.stackW (fun _ => t)).X N ω) ((m.tbl 0).v N) = _
      rw [hX N ω, overlap_smul ht]
      rfl
    simpa only [hfun] using h
  · intro N
    filter_upwards [law.topSimple N] with ω hω
    have h1 := topSimple_congr_gram ((m.stackW (fun _ => t)).X N ω) (t • m.stack.X N ω)
      (by rw [hX N ω])
    have h2 := topSimple_gram_smul ht (m.stack.X N ω)
    exact h1.mpr (h2.mpr hω)

/-- `cor.2`, model half: the heteroscedastic law at the binary weights `1_S` follows from the
`SingleTableLaw` of the stack of the sub-collection `S`. The limit is
`Scalars.binaryStackSVDLimit S θ c` through `Scalars.Lw_binary`. -/
theorem heteroLaw_binary_of_singleTableLaw (m : MultiTableModel μ M n d)
    (S : Finset (Fin M)) [NeZero S.card] (c : Fin M → ℝ)
    (law : (m.restrict S).stack.SingleTableLaw (∑ i ∈ S, c i)) :
    m.HeteroLaw (fun i => if i ∈ S then 1 else 0) c := by
  classical
  have hc' : ∑ j, (fun j => c (S.orderEmbOfFin rfl j)) j = ∑ i ∈ S, c i := sum_orderEmb S c
  have law' : (m.restrict S).stack.SingleTableLaw
      (∑ j, (fun j => c (S.orderEmbOfFin rfl j)) j) := by
    rw [hc']
    exact law
  have h := (m.restrict S).prop_stacksvd_general (fun j => c (S.orderEmbOfFin rfl j)) law'
  have hlim : Scalars.Lw (fun i => (m.tbl i).θ) c (fun i => if i ∈ S then 1 else 0)
      = stackSVDLimit (fun j => ((m.restrict S).tbl j).θ)
          (fun j => c (S.orderEmbOfFin rfl j)) := by
    rw [Scalars.Lw_binary,
      Scalars.binaryStackSVDLimit_eq_stackSVDLimit S (fun i => (m.tbl i).θ) c
        (S.orderIsoOfFin rfl).toEquiv]
    simp only [restrict_tbl]
    rfl
  refine ⟨?_, ?_⟩
  · rw [hlim]
    have hpf : ∀ (N : ℕ) (ω : Ω N),
        m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω
          = overlap ((m.restrict S).stack.X N ω) ((m.restrict S).stack.v N) := by
      intro N ω
      change overlap ((m.stackW (fun i => if i ∈ S then 1 else 0)).X N ω) ((m.tbl 0).v N) = _
      rw [overlap_congr_gram _ ((m.restrict S).stack.X N ω) (m.gram_restrict S N ω)]
      congr 1
      exact (m.hv (S.orderEmbOfFin rfl 0) 0 N).symm
    simpa only [hpf] using h
  · intro N
    filter_upwards [law.topSimple N] with ω hω
    exact (topSimple_congr_gram ((m.stackW (fun i => if i ∈ S then 1 else 0)).X N ω)
      ((m.restrict S).stack.X N ω) (m.gram_restrict S N ω)).mpr hω

/-- `cor.2`, model half, Layer 2: the heteroscedastic law at the binary weights `1_S` from the
proportional regime of each table and the joint Gaussian law alone. -/
theorem heteroLaw_binary_of_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (S : Finset (Fin M)) [NeZero S.card] (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroLaw (fun i => if i ∈ S then 1 else 0) c := by
  classical
  have hne : Nonempty (Fin S.card) := ⟨⟨0, Nat.pos_of_ne_zero (NeZero.ne S.card)⟩⟩
  have hcpos : 0 < ∑ j, c (S.orderEmbOfFin rfl j) :=
    Finset.sum_pos (fun j _ => hc _) Finset.univ_nonempty
  have hlaw : (m.restrict S).stack.SingleTableLaw (∑ j, c (S.orderEmbOfFin rfl j)) :=
    SpikedModel.singleTableLaw_of_gaussian hcpos (m.restrict S).stack
      ((m.restrict S).stack_regime (fun j => c (S.orderEmbOfFin rfl j)) fun j => hreg _)
      ((m.restrict S).stack_law (m.restrict_jointGaussianNoise S hG))
  rw [sum_orderEmb S c] at hlaw
  exact m.heteroLaw_binary_of_singleTableLaw S c hlaw

omit [NeZero M] in
/-- **`cor.2`** (`main_paper.tex:429`), Layer 2: on a nonempty subset `S` the binary-weighted
stacksvd performance converges to `Scalars.binaryStackSVDLimit S θ c`, with every random matrix
theory hypothesis discharged. -/
theorem stackPerfW_binary_tendsto_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (S : Finset (Fin M)) [NeZero S.card] (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
      (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) c) := by
  have h := (m.heteroLaw_binary_of_gaussian S c hc hreg hG).align
  rwa [Scalars.Lw_binary] at h

omit [NeZero M] in
/-- The subset maximum of `cor.2` (`main_paper.tex:442`), Layer 2: some nonempty subset attains
`Scalars.binaryStackSVDLimitMax`, and the binary-weighted stacksvd performance there converges
to that maximum. -/
theorem exists_binary_tendsto_max_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∃ S : Finset (Fin M), S.Nonempty ∧
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
        (Scalars.binaryStackSVDLimitMax (fun i => (m.tbl i).θ) c) := by
  obtain ⟨S, hSne, hS⟩ := Scalars.exists_nonempty_eq_max (fun i => (m.tbl i).θ) c
  have : NeZero S.card := ⟨(Finset.card_pos.mpr hSne).ne'⟩
  refine ⟨S, hSne, ?_⟩
  rw [hS]
  exact m.stackPerfW_binary_tendsto_gaussian S c hc hreg hG

omit [NeZero M] in
/-- **`thm:stacksvd_binary_optimal_svd_stack`** (`main_paper.tex:626`), Gaussian noise: with
every kept table at `c_i ≤ 1`, the binary weighting on `S`, the set of tables above their own
threshold `c_i < θ_i⁴`, converges to `Scalars.binaryStackSVDLimit S θ c`, and that limit is
at least optimally weighted svdstack. Every random matrix theory hypothesis is discharged by
`hreg` and `hG`; `hdet` only says `S` is nonempty. Components:
`stackPerfW_binary_tendsto_gaussian` for the convergence, `Scalars.svdstackOpt_le_binary` for
the inequality. -/
theorem thm_stacksvd_binary_optimal_svd_stack_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1) (hdet : ∃ i, c i < (m.tbl i).θ ^ 4)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW
        (fun i => if i ∈ Finset.univ.filter (fun i => c i < (m.tbl i).θ ^ 4) then 1 else 0)
        N ω)
      (Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
        (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
            (fun i => (m.tbl i).θ) c := by
  classical
  set S : Finset (Fin M) := Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4 with hSdef
  have hSne : S.Nonempty := by
    obtain ⟨i, hi⟩ := hdet
    exact ⟨i, hSdef ▸ Finset.mem_filter.mpr ⟨Finset.mem_univ i, hi⟩⟩
  have : NeZero S.card := ⟨(Finset.card_pos.mpr hSne).ne'⟩
  exact ⟨m.stackPerfW_binary_tendsto_gaussian S c hc hreg hG,
    Scalars.svdstackOpt_le_binary hc hc1⟩

end StackW

end MultiTableModel

end StackedSVD

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! ### `prop:binarystacksvd_inadmissable`, model half

The scalar content is in `Scalars.lean`. These two corollaries put it on the estimators of a
model whose `θ_i` are all `1`, with the two black-box laws as hypotheses: the aspect ratios
enter only through the law, as everywhere else in this file.
-/

/-- On the instance of `prop:binarystacksvd_inadmissable`, every binary weighting of stacksvd
tends to `0` in probability. `S` is any nonempty subset, so this covers the optimal binary
weighting. -/
theorem stackPerfW_binary_inad (m : MultiTableModel μ M n d)
    (hθ : ∀ i, (m.tbl i).θ = 1) (S : Finset (Fin M)) [NeZero S.card]
    (law : (m.restrict S).stack.SingleTableLaw (∑ i ∈ S, Scalars.inadC M i)) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0 := by
  have hfun : (fun i => (m.tbl i).θ) = Scalars.inadTheta M := funext fun i => hθ i
  have h := (m.heteroLaw_binary_of_singleTableLaw S (Scalars.inadC M) law).align
  rw [Scalars.Lw_binary, hfun, Scalars.inad_binaryStackSVDLimit_eq_zero S] at h
  exact h

/-- On the same instance, optimally weighted stacksvd tends to a value above `1 - ε`, as soon
as `M ≥ e^{-γ} exp(2/ε)`. -/
theorem stackPerfW_opt_inad (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M)
    (law : m.HeteroLaw (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M))
      (Scalars.inadC M)) :
    TendstoInProb μ
        (fun N ω => m.stackPerfW
          (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
        (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M) := by
  have hfun : (fun i => (m.tbl i).θ) = Scalars.inadTheta M := funext fun i => hθ i
  have hne : ∃ i, (m.tbl i).θ ≠ 0 := ⟨(0 : Fin M), by rw [hθ]; exact one_ne_zero⟩
  refine ⟨m.thm_stacksvd_weighted (Scalars.inadC M) (Scalars.inadC_pos M) hne law, ?_⟩
  rw [hfun]
  exact (Scalars.binarystacksvd_inadmissable hε hε1 M hM).1

end MultiTableModel

end StackedSVD
