/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Weighted

/-!
# `thm:simple_thm1`, svdstack half, for every `M ≥ 1`

`notes/SERVER_TODO.md`, "After the third external review", item P3. `main_paper.tex`
(`label{thm:simple_thm1}`) states, for equal tables `θ_i = θ₀`, `c_i = c₀` and every
`M ≥ 1`:

`|⟨v̂_svdstack, v⟩|² → 1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²-(M-1)c₀)` if `θ₀ ≥ c₀^{1/4}`, else `0`.

`thm_simple_thm1_svdstack_gaussian` (`SVDStack/Main.lean:424`) proves this only for `2 ≤ M`,
because its route through `thm_svd_stack_general_gaussian` needs two distinct tables with
`0 < β_i` to pin down `TopSimple (Abeta β)` (`svdstackLimit`'s multiplicity warning,
`SVDStack/Defs.lean:200`). This file adds the `M = 1` case, following `prop:single_table`
directly: at `M = 1`, `svdstackGram = P 0` and `svdstackPerf` is the single-table overlap of
table `0`, whose limit is `betaSq θ₀ c₀` by `prop:single_table` itself, with no `Abeta`
machinery involved. The two cases are then assembled with the existing `θ₀⁴ ≤ c₀` case
(`thm_svd_stack_general_zero_gaussian`, `SVDStack/Main.lean:394`, at `β = 0`) into the
paper's if/else form.

The boundary `θ₀ = c₀^{1/4}` (the paper's `≥`) matches the Lean strict threshold `c₀ < θ₀⁴`:
at equality the closed form gives `1 - (c₀+θ₀²)/(Mc₀+θ₀²-(M-1)c₀) = 1 - (c₀+θ₀²)/(c₀+θ₀²)
= 0`, the same value as the `else` branch, so `<` versus `≤` makes no difference to the
limit.

## Content

* `topSpace_P_eq` comes from `SVDStack/Weighted.lean`, where it is public since the
  2026-09-02 dedupe; this file no longer keeps a copy.
* `svdstackPerf_eq_overlap_one`: the unweighted analogue of `SVDStack/Weighted.lean`'s private
  `svdstackPerfW_eq_overlap_one`, proved directly from `topSpace_P_eq` rather than by
  specializing the weighted lemma at `w = 1`.
* `thm_simple_thm1_svdstack_gaussian_full`: the paper's if/else statement for every `M`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-- At `M = 1` the unweighted svdstack Gram matrix `∑_i P_i` is `P_0`, so its top eigenspace
is the top eigenspace of `X_0ᵀ X_0` (`topSpace_P_eq`), and `svdstackPerf` is the single-table
overlap of table `0`. The unweighted analogue of `SVDStack/Weighted.lean`'s private
`svdstackPerfW_eq_overlap_one`. -/
private theorem svdstackPerf_eq_overlap_one (m : MultiTableModel μ M n d)
    (h1 : ¬ 1 < Fintype.card (Fin M)) (N : ℕ) (ω : Ω N) :
    m.svdstackPerf N ω = overlap ((m.tbl 0).X N ω) ((m.tbl 0).v N) := by
  have hM : M = 1 := by
    have h0 := NeZero.pos M
    simp only [Fintype.card_fin, not_lt] at h1
    omega
  subst hM
  have hEq : m.svdstackGram N ω = m.P 0 N ω := by
    rw [MultiTableModel.svdstackGram, Fin.sum_univ_one]
  have hsp : topSpace (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω)
      = topSpace (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω) := by
    rw [topSpace_congr hEq (m.isHermitian_svdstackGram N ω) (m.isHermitian_P 0 N ω)]
    exact m.topSpace_P_eq 0 N ω
  have hproj : topProj (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) ((m.tbl 0).v N)
      = topProj (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω) ((m.tbl 0).v N) := by
    change (topSpace (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω)).starProjection
        ((m.tbl 0).v N)
      = (topSpace (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω)).starProjection
        ((m.tbl 0).v N)
    rw [hsp]
  change ‖topProj (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) ((m.tbl 0).v N)‖ ^ 2 = _
  rw [hproj]
  rfl

omit [NeZero M] in
/-- **`thm:simple_thm1`**, svdstack half, for every `M ≥ 1` (`main_paper.tex`,
`label{thm:simple_thm1}`), Gaussian noise: the tables share signal `θ₀` and regime `c₀`. Above
the detectability threshold `c₀ < θ₀⁴` (the paper's `θ₀ ≥ c₀^{1/4}`, tight at the boundary,
see the module docstring) the limit is the paper's closed form
`1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²-(M-1)c₀)`; otherwise it is `0`. `2 ≤ M` reuses
`thm_simple_thm1_svdstack_gaussian`; `M = 1` reduces to `prop:single_table` for table `0`
alone (`svdstackPerf_eq_overlap_one`, `SpikedModel.singleTableLaw_of_gaussian`), matching the
closed form since `1 - (c₀+θ₀²)/(θ₀⁴+θ₀²-0) = (θ₀⁴-c₀)/(θ₀⁴+θ₀²) = betaSq θ₀ c₀`; `θ₀⁴ ≤ c₀`
reuses `thm_svd_stack_general_zero_gaussian` at `β = 0`. -/
theorem thm_simple_thm1_svdstack_gaussian_full [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω)
      (if c₀ < θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) else 0) := by
  by_cases hdet : c₀ < θ₀ ^ 4
  · rw [if_pos hdet]
    by_cases hM2 : 2 ≤ M
    · exact m.thm_simple_thm1_svdstack_gaussian θ₀ c₀ hc hdet hM2 hθ hreg hG
    · have hM1 : M = 1 := by
        have h0 := NeZero.pos M
        omega
      have h1 : ¬ 1 < Fintype.card (Fin M) := by
        rw [Fintype.card_fin]; omega
      have hfe : (fun N ω => m.svdstackPerf N ω)
          = (fun N ω => overlap ((m.tbl 0).X N ω) ((m.tbl 0).v N)) :=
        funext fun N => funext fun ω => m.svdstackPerf_eq_overlap_one h1 N ω
      have hlaw : (m.tbl 0).SingleTableLaw c₀ :=
        SpikedModel.singleTableLaw_of_gaussian hc (m.tbl 0) (hreg 0)
          (m.gaussianNoise_of_joint hG 0)
      have hconv0 : TendstoInProb μ (fun N ω => overlap ((m.tbl 0).X N ω) ((m.tbl 0).v N))
          (betaSq θ₀ c₀) := by
        have := hlaw.align
        rwa [hθ 0] at this
      have hclosed : (1 : ℝ) - (c₀ + θ₀ ^ 2)
          / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) = betaSq θ₀ c₀ := by
        have hMr : (M : ℝ) = 1 := by exact_mod_cast hM1
        have hθ4pos : (0 : ℝ) < θ₀ ^ 4 := by linarith
        have hden0 : (0 : ℝ) < θ₀ ^ 4 + θ₀ ^ 2 := by nlinarith [sq_nonneg θ₀]
        have hden : θ₀ ^ 4 + θ₀ ^ 2 ≠ 0 := hden0.ne'
        rw [betaSq, if_pos hdet, hMr,
          show (1 : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - (1 - 1) * c₀ = θ₀ ^ 4 + θ₀ ^ 2 from by ring,
          show θ₀ ^ 4 - c₀ = (θ₀ ^ 4 + θ₀ ^ 2) - (c₀ + θ₀ ^ 2) from by ring, sub_div,
          div_self hden]
      rw [hclosed, hfe]
      exact hconv0
  · rw [if_neg hdet]
    have hβ0 : beta θ₀ c₀ = 0 := by rw [beta, betaSq, if_neg hdet, Real.sqrt_zero]
    exact m.thm_svd_stack_general_zero_gaussian (fun _ => c₀) (fun _ => beta θ₀ c₀)
      (fun _ => hc) (fun i => by rw [hθ i]) (fun _ => hβ0) hreg hG

end MultiTableModel

end StackedSVD
