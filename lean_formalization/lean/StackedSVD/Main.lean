/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Full
import StackedSVD.SVDStack.Gram
import StackedSVD.SVDStack.EntrywiseSigned
import StackedSVD.StackSVD
import StackedSVD.StackSVD.Main
import StackedSVD.SVDStack.Main
import StackedSVD.SVDStack.Simple
import StackedSVD.StackSVDWeighted
import StackedSVD.StackSVD.Weighted
import StackedSVD.RMT.Het.Sup
import StackedSVD.SVDStack.Weighted
import StackedSVD.SVDStack.Rayleigh
import StackedSVD.Secular
import StackedSVD.StrictFacades
import StackedSVD.Existence
import StackedSVD.ThetaEst
import StackedSVD.MLE
import StackedSVD.MLEConverse
import StackedSVD.MLEMarginal.Main
import StackedSVD.RankR.GramR
import StackedSVD.RankR.RMT.TableLawGaussian
import StackedSVD.RankR.GeneralMain
import StackedSVD.RankR.GeneralGaussian
import StackedSVD.RankR.SubspaceGaussian
import StackedSVD.RankR.Example
import StackedSVD.RankR.Het.Sup
import StackedSVD.RankR.SingleWeight.Het.Sup
import StackedSVD.RankR.SingleWeight.Suboptimality

/-!
# The statement file

Every result of the paper, in the form proved here for Gaussian noise, in the order of
`docs/THEOREMS.md`. Each theorem restates the signature of the theorem that proves it and
is proved by that theorem, so the checker guarantees that the statement below is exactly
what the tree proves. The proofs of the tree are elsewhere; this file has none. Run
`scripts/check_axioms.sh` to confirm that every declaration below depends only on
`propext`, `Classical.choice` and `Quot.sound` (do not add `#print axioms` here: it would
add an info line to every build).

Two restatements are proved by composition rather than by a single theorem of the tree:
`lem_delocalization` (`docs/THEOREMS.md` 3.2) and `lem_general_rank_delocalization`
(`docs/THEOREMS.md` 6.1) are both Layer 1 theorems over a `law` hypothesis; here that
hypothesis is discharged table by table by the matching `_of_gaussian` theorem, exactly as
the tree's own Gaussian facades do it elsewhere (`SVDStack/Main.lean`,
`RankR/GeneralGaussian.lean`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator MatrixOrder

namespace StackedSVD.Main

/-! ## 3.1 `prop:single_table` -/

/-- **`prop:single_table`** (`docs/THEOREMS.md` 3.1). Proved by
`SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean`). One table: the overlap
`⟨v̂, v⟩² → β²`, and the top eigenvalue of `XᵀX` converges to `θ² + 1 + c + c/θ²` above the
threshold `θ⁴ > c` and to the bulk edge `(1 + √c)²` below it. -/
theorem prop_single_table {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ} [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ}
    (hc : 0 < c) (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) :
    m.SingleTableLaw c :=
  StackedSVD.SpikedModel.singleTableLaw_of_gaussian hc m hreg hG

/-! ## 3.2 `lem:delocalization` -/

/-- **`lem:delocalization`** (`docs/THEOREMS.md` 3.2). Proved by composing
`MultiTableModel.lem_delocalization` (`SVDStack/Gram.lean`), a Layer 1 theorem, with
`SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean`) on each table and
`MultiTableModel.JointGaussianNoise.indepNoise` (`SVDStack/Gram.lean`). Per-table estimates of
two different tables, each signed so that `⟨v̂_i, v⟩ ≥ 0`: `⟨v̂_i, v̂_j⟩ → β_i β_j`. -/
theorem lem_delocalization {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j)) :=
  StackedSVD.MultiTableModel.lem_delocalization m c
    (fun k => StackedSVD.SpikedModel.singleTableLaw_of_gaussian (hc k) (m.tbl k) (hreg k)
      (StackedSVD.MultiTableModel.gaussianNoise_of_joint m hG k))
    (StackedSVD.MultiTableModel.JointGaussianNoise.indepNoise hG) hij

/-! ## 3.3 `lem:entrywise_conv_eigenvec` -/

/-- **`lem:entrywise_conv_eigenvec`** (`docs/THEOREMS.md` 3.3). Proved by
`lem_entrywise_conv_eigenvec_signed` (`SVDStack/EntrywiseSigned.lean`), the paper's signed
coordinate form: model-free, no noise law. With the sign of `x̂ = v_max(Ṽ Ṽᵀ)` chosen so that
`⟨x, x̂⟩ ≥ 0`, where `x = v_max(A_β)`, every coordinate satisfies `x̂_m → x_m` in probability. -/
theorem lem_entrywise_conv_eigenvec {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ}
    (β : Fin M → ℝ) (G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (Abeta β i j))
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β)) (m : Fin M) :
    TendstoInProb μ
      (fun N ω => signPos ⟪vMax (G N ω) (hsymm N ω), vMax (Abeta β) (isHermitian_Abeta β)⟫_ℝ
        * vMax (G N ω) (hsymm N ω) m)
      (vMax (Abeta β) (isHermitian_Abeta β) m) :=
  StackedSVD.lem_entrywise_conv_eigenvec_signed β G hsymm hconv hsimple m

/-! ## 3.4 `prop:stacksvd_general` -/

/-- **`prop:stacksvd_general`** (`docs/THEOREMS.md` 3.4). Proved by
`MultiTableModel.prop_stacksvd_general_gaussian` (`StackSVD.lean`). Unweighted stackSVD:
`⟨v̂, v⟩² → ((Σθ_i²)² − Σc_i) / ((Σθ_i²)(Σθ_i² + 1))` above the threshold `(Σθ_i²)² > Σc_i`,
`0` below. -/
theorem prop_stacksvd_general {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c) :=
  StackedSVD.MultiTableModel.prop_stacksvd_general_gaussian m c hc hreg hG

/-! ## 3.5 `thm:svd_stack_general` -/

/-- **`thm:svd_stack_general`** (`docs/THEOREMS.md` 3.5). Proved by
`MultiTableModel.thm_svd_stack_general_gaussian` (`SVDStack/Main.lean`). Unweighted SVDstack:
`⟨v̂, v⟩² → (βᵀ v_max(A_β))² / λ_max(A_β)` with `A_β = ββᵀ + diag(1 − β_i²)`, at least two
tables with `β_i > 0`. -/
theorem thm_svd_stack_general {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β) :=
  StackedSVD.MultiTableModel.thm_svd_stack_general_gaussian m c β hc hβdef hthr hreg hG

/-! ## 3.6 `thm:simple_thm1` (cor. 1) -/

/-- **`thm:simple_thm1`** (cor. 1) (`docs/THEOREMS.md` 3.6), stackSVD half. Proved by
`MultiTableModel.thm_simple_thm1_stacksvd_gaussian` (`StackSVD.lean`). Identical tables
(`θ_i = θ₀`, `c_i = c₀`): the closed form of unweighted stackSVD, both branches of the
threshold, every `M ≥ 1`. -/
theorem thm_simple_thm1_stacksvd {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (if c₀ < (M : ℝ) * θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) else 0) :=
  StackedSVD.MultiTableModel.thm_simple_thm1_stacksvd_gaussian m θ₀ c₀ hc hθ hreg hG

/-- **`thm:simple_thm1`** (cor. 1) (`docs/THEOREMS.md` 3.6), SVDstack half, every `M ≥ 1` and
both branches. Proved by `MultiTableModel.thm_simple_thm1_svdstack_gaussian_full`
(`SVDStack/Simple.lean`). Identical tables: the closed form of unweighted SVDstack. -/
theorem thm_simple_thm1_svdstack {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω)
      (if c₀ < θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) else 0) :=
  StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian_full m θ₀ c₀ hc hθ hreg hG

/-! ## 3.7 `cor.2` (binary weighting) -/

/-- **`cor.2`** (binary weighting) (`docs/THEOREMS.md` 3.7). Proved by
`MultiTableModel.stackPerfW_binary_tendsto_gaussian` (`StackSVDWeighted.lean`). Keep any
nonempty subset `S` of tables: binary-weighted stackSVD on `S` converges to the stackSVD
formula on `S`. -/
theorem cor2_binary {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (S : Finset (Fin M)) [NeZero S.card] (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
      (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) c) :=
  StackedSVD.MultiTableModel.stackPerfW_binary_tendsto_gaussian m S c hc hreg hG

/-- **`cor.2`** (binary weighting) (`docs/THEOREMS.md` 3.7), the subset maximum. Proved by
`MultiTableModel.exists_binary_tendsto_max_gaussian` (`StackSVDWeighted.lean`). The best
nonempty subset exists, and binary-weighted stackSVD there attains that maximum. -/
theorem cor2_binary_max {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∃ S : Finset (Fin M), S.Nonempty ∧
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
        (Scalars.binaryStackSVDLimitMax (fun i => (m.tbl i).θ) c) :=
  StackedSVD.MultiTableModel.exists_binary_tendsto_max_gaussian m c hc hreg hG

/-! ## 4.1 `thm:stacksvd_weighted` -/

/-- **`thm:stacksvd_weighted`** (`docs/THEOREMS.md` 4.1). Proved by
`MultiTableModel.thm_stacksvd_weighted_gaussian` (`RMT/Het/Sup.lean`). Weights
`w_i ∝ θ_i / √(θ_i² + c_i)` are optimal; the limit `γ⋆` solves
`Σ θ_i⁴ (1 − x)/(c_i + xθ_i²) = 1` when `Σ θ_i⁴/c_i > 1`, else `0`. -/
theorem thm_stacksvd_weighted {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c) :=
  StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian m c hc hθ hreg hG

/-- **`thm:stacksvd_weighted`** (`docs/THEOREMS.md` 4.1), optimality. Proved by
`MultiTableModel.thm_stacksvd_weighted_gaussian_opt` (`RMT/Het/Sup.lean`). No other nonzero
weighting does better than the optimal weights. -/
theorem thm_stacksvd_weighted_opt {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ ∀ w : Fin M → ℝ, (∃ i, w i ≠ 0) →
        TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
          (Scalars.Lw (fun i => (m.tbl i).θ) c w)
        ∧ Scalars.Lw (fun i => (m.tbl i).θ) c w
            ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c :=
  StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_opt m c hc hθ hreg hG

/-! ## 4.2 `thm:svdstack_weighted` -/

/-- **`thm:svdstack_weighted`** (`docs/THEOREMS.md` 4.2). Proved by
`MultiTableModel.thm_svdstack_weighted_gaussian` (`SVDStack/Weighted.lean`). Weights
`w_i = θ_i √((θ_i² + 1)/(θ_i² + c_i)) 1{θ_i⁴ > c_i}` give `⟨v̂, v⟩² → S/(S + 1)`,
`S = Σ β_i²/(1 − β_i²)`. -/
theorem thm_svdstack_weighted {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β) :=
  StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian m c β hc hβdef hthr hreg hG

/-- **`thm:svdstack_weighted`** (`docs/THEOREMS.md` 4.2), optimality with the uniform bound of
finding E1. Proved by `MultiTableModel.thm_svdstack_weighted_gaussian_opt_full`
(`SVDStack/Rayleigh.lean`). No nonzero weighting does better than the optimal weights, and
with probability tending to one no weight vector beats the optimum by more than `ε`. -/
theorem thm_svdstack_weighted_opt {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    (TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)
      ∧ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → (∀ i, 0 ≤ w i) →
          TopSimple (AbetaW w β) (isHermitian_AbetaW w β) →
          TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β)
          ∧ svdstackLimitW w β ≤ svdstackLimitOpt β)
    ∧ ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
        svdstackLimitOpt β + ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0) :=
  StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian_opt_full m c β hc hβdef hthr hreg hG

/-! ## 4.3 `lem:secular_equation` -/

/-- **`lem:secular_equation`**, claim 1 (`docs/THEOREMS.md` 4.3). Proved by
`Secular.det_Rmat_sub` (`Secular.lean`). The determinant identity for a diagonal-plus-rank-one
matrix `R = diag(σ) + q qᵀ`: deterministic, finite matrix algebra, no limit. -/
theorem lem_secular_equation_det {nn : ℕ} (σ q : Fin nn → ℝ) {lam : ℝ}
    (hlam : ∀ i, σ i ≠ lam) :
    (Secular.Rmat σ q - lam • (1 : Matrix (Fin nn) (Fin nn) ℝ)).det
      = (∏ i, (σ i - lam)) * Secular.secularDiag σ q lam :=
  StackedSVD.Secular.det_Rmat_sub σ q hlam

/-- **`lem:secular_equation`**, claim 3, the eigenvalue (`docs/THEOREMS.md` 4.3). Proved by
`Secular.lamMax_Rmat_eq` (`Secular.lean`). A root of the secular equation above the diagonal
entries of `σ` is the largest eigenvalue of `R` (since F30, 2026-09-07, the root hypothesis
alone: it forces `q ≠ 0` and `0 < nn`). -/
theorem lem_secular_equation_lamMax {nn : ℕ} {σ q : Fin nn → ℝ}
    {lam : ℝ} (hlam : (⨆ i, σ i) < lam) (e : Secular.secularDiag σ q lam = 0) :
    lamMax (Secular.Rmat σ q) (Secular.isHermitian_Rmat σ q) = lam :=
  StackedSVD.Secular.lamMax_Rmat_eq hlam e

/-- **`lem:secular_equation`**, claim 3, simplicity (`docs/THEOREMS.md` 4.3). Proved by
`Secular.topSimple_Rmat` (`Secular.lean`). That largest eigenvalue is simple. -/
theorem lem_secular_equation_simple {nn : ℕ} {σ q : Fin nn → ℝ}
    {lam : ℝ} (hlam : (⨆ i, σ i) < lam) (e : Secular.secularDiag σ q lam = 0) :
    TopSimple (Secular.Rmat σ q) (Secular.isHermitian_Rmat σ q) :=
  StackedSVD.Secular.topSimple_Rmat hlam e

/-! ## 5.1 `thm:stacksvd_binary_optimal_svd_stack` -/

/-- **`thm:stacksvd_binary_optimal_svd_stack`** (`docs/THEOREMS.md` 5.1). Proved by
`MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian` (`StackSVDWeighted.lean`).
Binary stackSVD (keep the detectable tables, at `c_i ≤ 1`) is at least as good as optimally
weighted SVDstack. -/
theorem thm_stacksvd_binary_optimal_svd_stack {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
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
            (fun i => (m.tbl i).θ) c :=
  StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian m c hc hc1 hdet
    hreg hG

/-- **`thm:stacksvd_binary_optimal_svd_stack`** (`docs/THEOREMS.md` 5.1), strict half. Proved
by `MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian_strict`
(`StrictFacades.lean`). With at least two tables above their own detection threshold, binary
stackSVD is strictly better than optimally weighted SVDstack. -/
theorem thm_stacksvd_binary_optimal_svd_stack_strict {Ω : ℕ → Type*}
    [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ}
    {d : ℕ → ℕ} [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1)
    (hcard : 2 ≤ (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4).card)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW
        (fun i => if i ∈ Finset.univ.filter (fun i => c i < (m.tbl i).θ ^ 4) then 1 else 0)
        N ω)
      (Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
        (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        < Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
            (fun i => (m.tbl i).θ) c :=
  StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian_strict m c hc hc1
    hcard hreg hG

/-! ## 5.2 `prop:dominance` -/

/-- **`prop:dominance`** (`docs/THEOREMS.md` 5.2). Proved by
`MultiTableModel.prop_dominance_gaussian` (`RMT/Het/Sup.lean`). Optimally weighted stackSVD is
at least as good as unweighted stackSVD and as optimally weighted SVDstack. -/
theorem prop_dominance {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c :=
  StackedSVD.MultiTableModel.prop_dominance_gaussian m c hc hθ hreg hG

/-- **`prop:dominance`** (`docs/THEOREMS.md` 5.2), both strict halves. Proved by
`MultiTableModel.prop_dominance_gaussian_strict` (`StrictFacades.lean`). Above the recovery
threshold, with at least two tables carrying signal and `θ_i²/c_i` not constant across tables,
optimally weighted stackSVD strictly dominates both comparisons. -/
theorem prop_dominance_strict {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hthr : 1 < ∑ i, (m.tbl i).θ ^ 4 / c i)
    (htwo : ∃ i j, i ≠ j ∧ (m.tbl i).θ ≠ 0 ∧ (m.tbl j).θ ≠ 0)
    (hnc : ∃ i j, (m.tbl i).θ ^ 2 * c j ≠ (m.tbl j).θ ^ 2 * c i) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c :=
  StackedSVD.MultiTableModel.prop_dominance_gaussian_strict m c hc hreg hG hthr htwo hnc

/-! ## 5.3 `prop:binarystacksvd_inadmissable` -/

/-- **`prop:binarystacksvd_inadmissable`** (`docs/THEOREMS.md` 5.3). Proved by
`Sat.Rank1.prop_binarystacksvd_inadmissable_exists` (`Existence.lean`). For every `ε ∈ (0,1)`
there is an explicit instance, built not assumed, with `M` tables where optimally weighted
stackSVD exceeds `1 − ε` while every binary weighting, optimally weighted SVDstack and both
unweighted methods have limit `0`. -/
theorem prop_binarystacksvd_inadmissable {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) :
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
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 :=
  StackedSVD.Sat.Rank1.prop_binarystacksvd_inadmissable_exists hε hε1

/-! ## 5.4 `thm:theta_est` -/

/-- **`thm:theta_est`** (`docs/THEOREMS.md` 5.4). Proved by
`MultiTableModel.thm_theta_est_gaussian` (`ThetaEst.lean`). Two tables, `θ_1⁴ > c_1`: the
estimator `θ̂_2` of `eq:theta_estimation` is consistent. -/
theorem thm_theta_est {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) (hG : m.JointGaussianNoise)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ :=
  StackedSVD.MultiTableModel.thm_theta_est_gaussian m hij hci hthr hG hregi hregj

/-! ## 5.5 `app:wstacksvd_mle`, the MLE identity and its marginalization -/

/-- **`app:wstacksvd_mle`**, the identity (`docs/THEOREMS.md` 5.5). Proved by
`MultiTableModel.mleLogLik_eq` (`MLE.lean`), a deterministic identity: the weighted stackSVD
objective is the marginal Gaussian log-likelihood of the random-effects model, at finite `N`,
up to an additive constant. -/
theorem mle_loglik_eq {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (N : ℕ) (ω : Ω N) (v : Fin (d N) → ℝ) (hv : v ⬝ᵥ v = 1) :
    m.mleLogLik c N ω v
      = m.mleConst c N ω
        + (d N : ℝ) / 2
          * (v ⬝ᵥ m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω *ᵥ v) :=
  StackedSVD.MultiTableModel.mleLogLik_eq m c hc N ω v hv

/-- **`app:wstacksvd_mle`**, the converse (`docs/THEOREMS.md` 5.5). Proved by
`MultiTableModel.mleLogLik_max_iff_mem_topSpace` (`MLEConverse.lean`), a deterministic
identity: a unit vector maximizes the log-likelihood if and only if it lies in the top
eigenspace of the weighted stack Gram matrix. -/
theorem mle_loglik_max_iff_mem_topSpace {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (N : ℕ) (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
        m.mleLogLik c N ω v ≤ m.mleLogLik c N ω (WithLp.ofLp x))
      ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
          (m.isHermitian_stackGramW _ N ω) :=
  StackedSVD.MultiTableModel.mleLogLik_max_iff_mem_topSpace m c hc N ω x hx

/-- **`app:wstacksvd_mle`**, the marginalization (`docs/THEOREMS.md` 5.5). Proved by
`MultiTableModel.thm_wstacksvd_mle_marginal` (`MLEMarginal/Main.lean`), a deterministic
identity: the row-marginal law of the random-effects model has the stated density, its
log-likelihood is the weighted stackSVD objective up to an additive constant, and its
maximizers are exactly the top eigenvectors of the weighted stack Gram matrix at the weights
`w_i = n_i N / d N`. -/
theorem thm_wstacksvd_mle_marginal {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    (m : MultiTableModel μ M n d) (N : ℕ)
    (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ,
        reJointLaw (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v
          = (Measure.pi fun i => lebesgueMatrix (n i N) (d N)).withDensity
              (fun X => ENNReal.ofReal
                (reDensity (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v X)))
    ∧ (∀ v : Fin (d N) → ℝ,
        reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v (fun i => (m.tbl i).X N ω)
          = m.mleLogLik (fun i => (n i N : ℝ) / d N) N ω v
            - ((∑ i, (n i N : ℝ)) * d N / 2) * Real.log (2 * Real.pi))
    ∧ ((∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
          reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v
              (fun i => (m.tbl i).X N ω)
            ≤ reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) (WithLp.ofLp x)
                (fun i => (m.tbl i).X N ω))
        ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ)
              (fun i => (n i N : ℝ) / d N)) N ω) (m.isHermitian_stackGramW _ N ω)) :=
  StackedSVD.MultiTableModel.thm_wstacksvd_mle_marginal m N ω x hx

/-! ## 5.6 The two remarks -/

/-- **`remark:stack_outperform_svd`** (`docs/THEOREMS.md` 5.6). Proved by
`Sat.Rank1.remark_stack_outperform_svd_exists` (`Existence.lean`). An explicit family, built
not assumed (every `M ≥ 2`, all `θ_i = c_i = 1`), where unweighted stackSVD has limit
`1 − 2/(M + 1)` and SVDstack has limit `0` under every weighting. -/
theorem remark_stack_outperform_svd (M : ℕ) [NeZero M] (hM : 2 ≤ M) :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ)
      (m : MultiTableModel μ M n d),
      (∀ i, (m.tbl i).θ = 1) ∧ (∀ i, (m.tbl i).Regime 1) ∧ m.JointGaussianNoise ∧
      TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
        (1 - 2 / ((M : ℝ) + 1)) ∧
      0 < 1 - 2 / ((M : ℝ) + 1) ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 ∧
      (∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
        ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0)) :=
  StackedSVD.Sat.Rank1.remark_stack_outperform_svd_exists M hM

/-- **`remark:svd_outperform_stack`** (`docs/THEOREMS.md` 5.6). Proved by
`Sat.Rank1.remark_svd_outperform_stack_exists` (`Existence.lean`). An explicit 2-table
instance, built not assumed (`θ = (√5, 4)`, `c = (1, 38.4)`), where unweighted SVDstack has
limit `8/9` and every binary stackSVD at most `2008/2310`. -/
theorem remark_svd_outperform_stack :
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
      (2008 : ℝ) / 2310 < 8 / 9 :=
  StackedSVD.Sat.Rank1.remark_svd_outperform_stack_exists

/-! ## 6.1 `lem:general_rank_delocalization` (Appendix D) -/

/-- **`lem:general_rank_delocalization`** (Appendix D) (`docs/THEOREMS.md` 6.1). Proved by
composing `UnalignedModelR.lem_general_rank_delocalization_general`
(`RankR/GeneralMain.lean`), a Layer 1 theorem, with `UnalignedModelR.tableLawR_of_gaussian_rk`
(`RankR/RMT/TableLawGaussian.lean`) on each table and
`UnalignedModelR.JointGaussianNoise.indepNoise` (`RankR/GramR.lean`). Cross-table entrywise
limits at general per-table rank: the off-diagonal entries of `ṼᵀV` and `ṼṼᵀ` converge to the
block objects `B_R` and `A_{β,R}`. -/
theorem lem_general_rank_delocalization {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {i i' : Fin M} (hii : i ≠ i') (j : Fin (rk i)) (j' : Fin (rk i')) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, m.vhatG i' j' N ω⟫_ℝ)
      (beta ((m.tbl i).θ j) (c i) * beta ((m.tbl i').θ j') (c i') *
        ((m.R i)ᵀ * m.R i') j j') :=
  StackedSVD.UnalignedModelR.lem_general_rank_delocalization_general m c
    (StackedSVD.UnalignedModelR.tableLawR_of_gaussian_rk m hc hreg hG)
    (StackedSVD.UnalignedModelR.JointGaussianNoise.indepNoise hG) hii j j'

/-! ## 6.2 `prop:general_rank_unweighted_svdstack` -/

/-- **`prop:general_rank_unweighted_svdstack`** (`docs/THEOREMS.md` 6.2). Proved by
`UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian`
(`RankR/GeneralGaussian.lean`). Unweighted SVDstack at general per-table rank `r_i`: the
paper's performance converges to `limitRG`, the rank-`r` analogue of `thm:svd_stack_general`.
-/
theorem prop_general_rank_unweighted_svdstack {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hrr : r ≤ rtot rk)
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRG N ω) (limitRG β m.R) :=
  StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian m c β hc
    hβdef hreg hrr hgap hG

/-! ## 6.3 `prop:stacksvd_subspace` -/

/-- **`prop:stacksvd_subspace`** (`docs/THEOREMS.md` 6.3), `r_i = 1`. Proved by
`UnalignedModel.prop_stacksvd_subspace_gaussian` (`RankR/SubspaceGaussian.lean`). StackSVD on
tables aligned inside a shared `r`-dimensional subspace: the subspace performance converges to
one rank-one performance per spike of the core matrix. -/
theorem prop_stacksvd_subspace {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω)
      (limitStackR (fun i => (m.tbl i).θ) m.R cc) :=
  StackedSVD.UnalignedModel.prop_stacksvd_subspace_gaussian m hG cc hreg hc

/-- **`prop:stacksvd_subspace`** (`docs/THEOREMS.md` 6.3), general `r_i`. Proved by
`UnalignedModelR.prop_stacksvd_subspace_general_gaussian` (`RankR/SubspaceGaussian.lean`). The
general-rank twin: the subspace performance of stackSVD converges to `limitStackRG`. -/
theorem prop_stacksvd_subspace_general {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}
    [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackRG N ω)
      (limitStackRG (fun i => (m.tbl i).θ) m.R cc) :=
  StackedSVD.UnalignedModelR.prop_stacksvd_subspace_general_gaussian m hG cc hreg hc

/-! ## 6.4 `thm:gen_rank_weight_svdstak` -/

/-- **`thm:gen_rank_weight_svdstak`** (`docs/THEOREMS.md` 6.4), uniform in the weight matrix
`W` (finding E2). Proved by
`UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian`
(`RankR/GeneralGaussian.lean`). Weighted SVDstack at general per-table rank `r_i`: the optimal
weight matrix attains the best limit `limitOptG`, every admissible `W` tends to a limit at
most `limitOptG`, and with probability tending to one no weight matrix beats it by more than
`ε`. -/
theorem thm_gen_rank_weight_svdstak {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r) (hrr : r ≤ rtot rk)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (by simpa using hrr)) ∧
      (∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by simp only [Fintype.card_fin]; omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (by simpa using hrr)) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        limitOptG β m.R (by simpa using hrr) + ε ≤ m.perfRGW W N ω}) atTop (𝓝 0) :=
  StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian m c β hc
    hβdef hreg hr hrr hG

/-! ## 6.5 The worked example of Section 7 -/

/-- The worked example of Section 7 (`docs/THEOREMS.md` 6.5), SVDstack half. Proved by
`RankR.Example.example_perfR_tendsto_gaussian` (`RankR/Example.lean`). The two-table, rank-2
example: the SVDstack rank-`r` performance has the closed form `svdstackEx` as a function of
the rotation angle `ψ`. -/
theorem example_perfR {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
    {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ} [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = RankR.Example.Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω)
      (RankR.Example.svdstackEx (beta θ c) (Real.sin ψ)) :=
  StackedSVD.RankR.Example.example_perfR_tendsto_gaussian m θ c ψ hc hθ hR hreg hG

/-- The worked example of Section 7 (`docs/THEOREMS.md` 6.5), stackSVD half. Proved by
`RankR.Example.example_perfStackR_tendsto_gaussian` (`RankR/SubspaceGaussian.lean`). The same
two-table, rank-2 example: the stackSVD rank-`r` performance has the closed form
`stacksvdEx`. -/
theorem example_perfStackR {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = RankR.Example.Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω)
      (RankR.Example.stacksvdEx θ c (Real.sin ψ)) :=
  StackedSVD.RankR.Example.example_perfStackR_tendsto_gaussian m θ c ψ hc hθ hR hreg hG

/-! ## 6.6 `thm:rank_r_svdstack` (Appendix E) -/

/-- **`thm:rank_r_svdstack`** (Appendix E) (`docs/THEOREMS.md` 6.6), aggregate clause. Proved
by `UnalignedModelR.thm_rank_r_svdstack_aggregate_gaussian` (`RankR/GeneralGaussian.lean`). At
the exactly aligned model, the aggregate Frobenius performance of optimally weighted SVDstack
sums to `Σ_j S_j/(S_j + 1)`. -/
theorem thm_rank_r_svdstack_aggregate {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
      (∑ j, Sagg β j / (Sagg β j + 1)) :=
  StackedSVD.UnalignedModelR.thm_rank_r_svdstack_aggregate_gaussian m c β hc hβdef hreg hR hM
    hG

/-- **`thm:rank_r_svdstack`** (Appendix E) (`docs/THEOREMS.md` 6.6), component clause. Proved
by `UnalignedModelR.thm_rank_r_svdstack_component_gaussian` (`RankR/GeneralGaussian.lean`). At
the exactly aligned model, with the aggregate terms `Sagg β` strictly decreasing, column `j`
of optimally weighted SVDstack converges componentwise to `S_j/(S_j + 1)`. -/
theorem thm_rank_r_svdstack_component {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hS : StrictAnti (Sagg β)) (hG : m.JointGaussianNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) :=
  StackedSVD.UnalignedModelR.thm_rank_r_svdstack_component_gaussian m c β hc hβdef hreg hR hM
    hS hG j

/-! ## 6.7 `thm:rank_r_stacksvd` (Appendix E) -/

/-- **`thm:rank_r_stacksvd`** (Appendix E) (`docs/THEOREMS.md` 6.7). Proved by
`UnalignedModelR.thm_rank_r_stacksvd_gaussian` (`RankR/Het/Sup.lean`). At the exactly aligned
model with a single weight per table, each component of weighted stackSVD converges to its
own `Scalars.gammaR` term, and the total Frobenius performance sums them. -/
theorem thm_rank_r_stacksvd {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    (∀ j : Fin r, TendstoInProb μ
        (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR m.thetaAligned c j)) ∧
      TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
        (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) :=
  StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian m c hc hR hreg hG hnz

/-! ## 6.8 `prop:gen_rank_stacksvd_singleweight` (Appendix D) -/

/-- **`prop:gen_rank_stacksvd_singleweight`** (Appendix D) (`docs/THEOREMS.md` 6.8). Proved by
`UnalignedModelR.prop_gen_rank_stacksvd_singleweight_gaussian`
(`RankR/SingleWeight/Het/Sup.lean`). StackSVD with one weight `w_i` per table and general
rotations `R_i`: the limit of `‖V̂ᵀV‖_F²` is the sum over the `r` roots `γ_ℓ` of a matrix
secular equation. -/
theorem prop_gen_rank_stacksvd_singleweight {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}
    [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    TendstoInProb μ (fun N ω => m.perfSW w N ω)
      (SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z) :=
  StackedSVD.UnalignedModelR.prop_gen_rank_stacksvd_singleweight_gaussian m w c hc hreg hG γ z
    hsep

/-- **`prop:gen_rank_stacksvd_singleweight`** (Appendix D) (`docs/THEOREMS.md` 6.8), the
paper's inner-product form. Proved by
`UnalignedModelR.prop_gen_rank_stacksvd_singleweight_inner_gaussian`
(`RankR/SingleWeight/Het/Sup.lean`). Per pair `(ℓ, k)`, the overlap `⟨v̂_ℓ, v_k⟩²` has limit
`swTerm_ℓ (z_ℓ)_k²`. -/
theorem prop_gen_rank_stacksvd_singleweight_inner {Ω : ℕ → Type*}
    [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ}
    {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ} [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l k : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2)
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
        (WithLp.ofLp (z l) k) ^ 2) :=
  StackedSVD.UnalignedModelR.prop_gen_rank_stacksvd_singleweight_inner_gaussian m w c hc hreg
    hG γ z hsep l k

/-! ## 6.9 `prop:singleweight_suboptimality` -/

/-- **`prop:singleweight_suboptimality`** (`docs/THEOREMS.md` 6.9). Proved by
`SingleWeight.Witness.prop_singleweight_suboptimality_gaussian`
(`RankR/SingleWeight/Suboptimality.lean`), for Gaussian noise with no hypothesis left. An
explicit two-table, rank-2 instance, built not assumed (`θ ≡ 8/5`, `c ≡ 1`), where unweighted
SVDstack strictly beats stackSVD restricted to one weight per table, for every choice of that
weight. -/
theorem prop_singleweight_suboptimality :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : UnalignedModelR μ 2 n d 2 (fun _ => 1)),
      (∀ i j, (m.tbl i).θ j = 8 / 5) ∧ (∀ i, (m.tbl i).Regime 1) ∧
      m.JointGaussianNoise ∧ m.R = Rone (RankR.Example.Rex 0) ∧
      TendstoInProb μ (fun N ω => m.perfRG N ω) (2 * betaSq (8 / 5) 1) ∧
      (∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
        TendstoInProb μ (fun N ω => m.perfSW w N ω) (SingleWeight.swLimitEx (8 / 5) 1 w) ∧
        SingleWeight.swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1) ∧
      2 * betaSq (8 / 5) 2 < 2 * betaSq (8 / 5) 1 :=
  StackedSVD.SingleWeight.Witness.prop_singleweight_suboptimality_gaussian

end StackedSVD.Main
