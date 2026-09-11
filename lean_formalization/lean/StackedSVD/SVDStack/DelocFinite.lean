/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Gram
import StackedSVD.RMT.Symmetry

/-!
# Item F1: the finite-`d` form of `lem:delocalization`

`lem_delocalization` (`SVDStack/Gram.lean`) proves the paper's limit
(`main_paper.tex:1111`): for `i ≠ j` the perpendicular part `⟪v̂_i, (I - v vᵀ) v̂_j⟫` tends to
`0` in probability. It chains two finite-`N` facts and then throws the constant away. This
file keeps the constant.

Under Gaussian noise, for `i ≠ j`, every `d_N ≥ 2` and every `ε > 0`,

`P(|⟪v̂_i, (I - v vᵀ) v̂_j⟫| ≥ ε) ≤ 1 / (ε² (d_N - 1))`.

The bound reads no `θ`, no `n_i`, no `c` and no `M`. It answers the referee question about a
rate at the one place where the paper leans on an external asymptotic theorem
(`main_paper.tex:1138`, Theorem 1 part 2 of `liu2023asymptotic`).

## Route

The composition is the one `lem_delocalization` already runs, at `η = ε²`.

1. On the almost sure event that table `j` has a simple top eigenvalue,
   `abs_inner_perpOf_le` (`SVDStack/Gram.lean`) replaces the perpendicular part by the
   measurable random direction `delocDir v_i X_j`, and `overlap_ge_inner_sq`
   (`Spectral.lean`) replaces the inner product by the overlap. So the event of the statement
   sits inside `{ε² ≤ overlap X_i (delocDir v_i X_j)}` up to a null set.
2. `measure_deloc_le` (`SVDStack/Gram.lean`) is the Fubini step: independence across tables
   bounds that measure by the supremum over deterministic unit directions orthogonal to `v_i`.
3. `measure_overlap_ge_le` (`RMT/Symmetry.lean`) bounds every term of the supremum by
   `1 / (ε² (d_N - 1))`, by right-rotation invariance and Markov.

Both inputs are proved and both are at fixed `N`. No new mathematics.

## Statements

* `MultiTableModel.measure_inner_perpOf_ge_le`, in `ℝ≥0∞`, the form `measure_overlap_ge_le`
  gives.
* `MultiTableModel.measure_inner_perpOf_ge_le_toReal`, the real-valued corollary.

No `sorry`, no `axiom`, no edit to any existing file.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- **Item F1, the finite-`d` form of `lem:delocalization`.** For `i ≠ j`, `d_N ≥ 2` and
`ε > 0`, under independent Gaussian tables,
`P(|⟪v̂_i, (I - v vᵀ) v̂_j⟫| ≥ ε) ≤ 1 / (ε² (d_N - 1))`. The bound holds at every `N`, and it
depends on no other datum of the model. `lem_delocalization` is its limit form. -/
theorem measure_inner_perpOf_ge_le (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) (hd : 2 ≤ d N) {ε : ℝ} (hε : 0 < ε) :
    μ N {ω | ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|}
      ≤ ENNReal.ofReal (1 / (ε ^ 2 * ((d N : ℝ) - 1))) := by
  have hsq : (0 : ℝ) < ε ^ 2 := by positivity
  calc μ N {ω | ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|}
      ≤ μ N {ω | ε ^ 2 ≤ overlap ((m.tbl i).X N ω)
          (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))} := by
        refine measure_mono_ae ?_
        filter_upwards [singleTableLaw_topSimple_of_gaussian (m.tbl j)
          (m.gaussianNoise_of_joint hG j) N] with ω hω hmem
        have hstart : ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ| := hmem
        have h1 : |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|
            ≤ |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ| :=
          abs_inner_perpOf_le ((m.tbl i).hv N) _ hω (m.mem_topSpace_vhat j N ω)
            (m.norm_vhat j N ω) _
        have h2 : ⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ ^ 2
            ≤ overlap ((m.tbl i).X N ω) (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)) :=
          overlap_ge_inner_sq _ _ (m.mem_topSpace_vhat i N ω) (m.norm_vhat i N ω)
        have h3 : ε ≤ |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ| :=
          hstart.trans h1
        have h4 : ε ^ 2
            ≤ |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ| ^ 2 := by
          nlinarith [hε.le, h3]
        rw [sq_abs] at h4
        exact h4.trans h2
    _ ≤ ⨆ w ∈ (m.tbl i).orthUnit N, μ N {ω | ε ^ 2 ≤ overlap ((m.tbl i).X N ω) w} :=
        m.measure_deloc_le hG.indepNoise hij N hsq
    _ ≤ ENNReal.ofReal (1 / (ε ^ 2 * ((d N : ℝ) - 1))) :=
        iSup₂_le fun w hw => (m.tbl i).measure_overlap_ge_le
          (m.gaussianNoise_of_joint hG i) N hd hsq hw

/-- The real-valued form of `measure_inner_perpOf_ge_le`. The bound `1 / (ε² (d_N - 1))` is
finite, so `toReal` is monotone here and the `ENNReal.ofReal` disappears. -/
theorem measure_inner_perpOf_ge_le_toReal (m : MultiTableModel μ M n d)
    (hG : m.JointGaussianNoise) {i j : Fin M} (hij : i ≠ j) (N : ℕ) (hd : 2 ≤ d N) {ε : ℝ}
    (hε : 0 < ε) :
    (μ N {ω | ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|}).toReal
      ≤ 1 / (ε ^ 2 * ((d N : ℝ) - 1)) := by
  have hd1 : (1 : ℝ) ≤ (d N : ℝ) - 1 := by
    have : (2 : ℝ) ≤ (d N : ℝ) := by exact_mod_cast hd
    linarith
  have hpos : (0 : ℝ) < ε ^ 2 * ((d N : ℝ) - 1) := by positivity
  have hnn : (0 : ℝ) ≤ 1 / (ε ^ 2 * ((d N : ℝ) - 1)) := by positivity
  have hmono := ENNReal.toReal_mono ENNReal.ofReal_ne_top
    (m.measure_inner_perpOf_ge_le hG hij N hd hε)
  rwa [ENNReal.toReal_ofReal hnn] at hmono

end MultiTableModel

end StackedSVD
