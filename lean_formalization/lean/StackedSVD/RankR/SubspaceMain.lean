/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/

import StackedSVD.RankR.Subspace
import StackedSVD.RankR.SubspaceG
import StackedSVD.RankR.Example

/-!
# `prop:stacksvd_subspace`: the stacksvd subspace limit at `r_i = 1` and at general `r_i`

Moved out of `RankR/Subspace.lean` and `RankR/SubspaceG.lean` (F28, 2026-09-08):
`prop_stacksvd_subspace`, the `RankR.Example` namespace (the worked example
`eq:psi_equation`, stacksvd half), and `prop_stacksvd_subspace_general`. No proof
changed.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! ### 4. `prop:stacksvd_subspace` -/

/-- `prop:stacksvd_subspace` (`main_paper.tex:834`) at `r_i = 1`, Layer 1 form: under the
rank-`r` spiked law of the stack at aspect ratio `‖c‖₁`, the stacksvd subspace performance
tends to one rank-one performance per spike of the core matrix.

Route: `perfStackR_eq_sum_spikeVec` turns the performance into the sum over the spike
directions; each term is `SubspaceLaw.align j`; a finite sum of limits in probability is the
limit of the sum (`TendstoInProb.finsum`, `Prob/TendstoInProb.lean`).

No hypothesis beyond the law is needed. In particular the eigenvalues of `C` may repeat,
which is what the paper claims at `main_paper.tex:842`. -/
theorem prop_stacksvd_subspace (m : UnalignedModel μ M n d r) (c : Fin M → ℝ)
    (law : m.SubspaceLaw (∑ i, c i)) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω)
      (limitStackR (fun i => (m.tbl i).θ) m.R c) := by
  have hfun : (fun N (ω : Ω N) => m.perfStackR N ω)
      = fun N (ω : Ω N) => ∑ j : Fin r,
          ‖specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r
            (m.spikeVec j N)‖ ^ 2 := by
    funext N ω
    exact m.perfStackR_eq_sum_spikeVec N ω
  rw [m.limitStackR_eq c, hfun]
  exact TendstoInProb.finsum (fun j => law.align j)

end UnalignedModel

end StackedSVD

/-! ### 5. The worked example `eq:psi_equation`, stacksvd half -/

namespace StackedSVD.RankR.Example

open MeasureTheory

/-- The core matrix of the paper's example is `Cex` (`main_paper.tex:2216`). -/
theorem cmat_rex (θ ψ : ℝ) : Cmat (fun _ : Fin 2 => θ) (Rex ψ) = Cex θ ψ := by
  rw [cmat_eq_sum, Cex, Finset.smul_sum]

/-- `hPerf_eq_betaSq` without the hypothesis `0 < lam`. At `lam ≤ 0` both sides are `0`:
`√lam = 0` and `betaSq 0 x = 0`, while the indicator of `hPerf` is false because `√(2c) ≥ 0`.
The endpoints `sin ψ = ±1` and `θ = 0` of the worked example need this case. -/
private theorem hPerf_eq_betaSq_all (c lam : ℝ) :
    hPerf c lam = betaSq (Real.sqrt lam) (2 * c) := by
  by_cases hlam : 0 < lam
  · exact hPerf_eq_betaSq hlam
  · have hlam' : lam ≤ 0 := not_lt.mp hlam
    have hs : Real.sqrt lam = 0 := Real.sqrt_eq_zero'.mpr hlam'
    have hb : betaSq (0 : ℝ) (2 * c) = 0 := by
      simp [betaSq]
    rw [hPerf, if_neg (not_lt.mpr (hlam'.trans (Real.sqrt_nonneg (2 * c)))), hs, hb]

/-- `eq:perf_rankr_ex_stacksvd` (`main_paper.tex:860`) as the limit of
`prop:stacksvd_subspace`: at `M = 2`, `r = 2`, `r_i = 1`, `θ_1 = θ_2 = θ`, `c_1 = c_2 = c`
and `R = Rex ψ`, the spectrum of the core matrix is `θ²(1 ± sin ψ)` and the limit is the
paper's `F(s) = h(θ²(1+s)) + h(θ²(1-s))`.

The identity needs no hypothesis: `stacksvdEx` and `limitStackR` agree at `θ = 0`, at
`c = 0` and at `sin ψ = ±1` as well, where one spike of `C` vanishes and both sides read `0`
(necessity scan rows F2 and F4 of `notes/archive/prop_stacksvd_subspace.md`, 2400 points, worst
5.6e-16).

Route: the trace and the determinant of `Cex` give the sum and the product of the two
eigenvalues, so the unsorted pair `(λ_0, λ_1)` is `(θ²(1+s), θ²(1-s))` in one of the two
orders. Both orders give the same sum. -/
theorem limitStackR_example (θ c ψ : ℝ) :
    limitStackR (fun _ : Fin 2 => θ) (Rex ψ) (fun _ => c)
      = stacksvdEx θ c (Real.sin ψ) := by
  have hpy := Real.sin_sq_add_cos_sq ψ
  have htr : (Cmat (fun _ : Fin 2 => θ) (Rex ψ)).trace
      = θ ^ 2 * (1 + Real.sin ψ) + θ ^ 2 * (1 - Real.sin ψ) := by
    rw [cmat_rex, cex_eq, Matrix.trace_smul, Matrix.trace_fin_two_of]
    simp only [smul_eq_mul]
    linear_combination θ ^ 2 * hpy
  have hdet : (Cmat (fun _ : Fin 2 => θ) (Rex ψ)).det
      = (θ ^ 2 * (1 + Real.sin ψ)) * (θ ^ 2 * (1 - Real.sin ψ)) := by
    rw [cmat_rex, cex_eq, Matrix.det_smul, Matrix.det_fin_two_of, Fintype.card_fin]
    linear_combination θ ^ 4 * hpy
  have hsum : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
      + (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 1
      = θ ^ 2 * (1 + Real.sin ψ) + θ ^ 2 * (1 - Real.sin ψ) := by
    have h := (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).trace_eq_sum_eigenvalues
    rw [htr, Fin.sum_univ_two] at h
    exact_mod_cast h.symm
  have hprod : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
      * (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 1
      = (θ ^ 2 * (1 + Real.sin ψ)) * (θ ^ 2 * (1 - Real.sin ψ)) := by
    have h := (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).det_eq_prod_eigenvalues
    rw [hdet, Fin.prod_univ_two] at h
    exact_mod_cast h.symm
  have hfac : ((isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
        - θ ^ 2 * (1 + Real.sin ψ))
      * ((isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
        - θ ^ 2 * (1 - Real.sin ψ)) = 0 := by
    linear_combination
      (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0 * hsum - hprod
  have hc2 : (∑ _i : Fin 2, c) = 2 * c := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    simp [two_mul]
  rw [limitStackR, Fin.sum_univ_two, hc2, stacksvdEx, hPerf_eq_betaSq_all,
    hPerf_eq_betaSq_all]
  rcases mul_eq_zero.mp hfac with h | h
  · have h0 : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
        = θ ^ 2 * (1 + Real.sin ψ) := by linarith
    have h1 : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 1
        = θ ^ 2 * (1 - Real.sin ψ) := by linarith
    rw [h0, h1]
  · have h0 : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 0
        = θ ^ 2 * (1 - Real.sin ψ) := by linarith
    have h1 : (isHermitian_Cmat (fun _ : Fin 2 => θ) (Rex ψ)).eigenvalues 1
        = θ ^ 2 * (1 + Real.sin ψ) := by linarith
    rw [h0, h1, add_comm]

section Convergence

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin 2 → ℕ → ℕ} {d : ℕ → ℕ}

/-- The stacksvd half of the paper's headline unaligned example, as a limit theorem: with
`eq:psi_equation` (`main_paper.tex:848`) and the rank-`r` law of the stack at aspect ratio
`2c`, the stacksvd performance converges in probability to `eq:perf_rankr_ex_stacksvd`.

The svdstack half is `example_perfR_tendsto` of `RankR/Example.lean`, and the two together
are the comparison of `main_paper.tex:864`. The proof is `prop_stacksvd_subspace` followed by
`limitStackR_example`; `∑ i : Fin 2, c = 2 * c` is `Fin.sum_univ_two`. Neither `0 < c` nor
`0 < θ` nor `ψ ∈ [0, π/2)` is needed, because `limitStackR_example` needs none of them. -/
theorem example_perfStackR_tendsto (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ) (law : m.SubspaceLaw (2 * c)) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω) (stacksvdEx θ c (Real.sin ψ)) := by
  have hc2 : (∑ i : Fin 2, (fun _ : Fin 2 => c) i) = 2 * c := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    simp [two_mul]
  have law' : m.SubspaceLaw (∑ i : Fin 2, (fun _ : Fin 2 => c) i) := by
    rw [hc2]; exact law
  have hlim := m.prop_stacksvd_subspace (fun _ => c) law'
  have hθf : (fun i => (m.tbl i).θ) = (fun _ : Fin 2 => θ) := funext hθ
  rw [hθf, hR, limitStackR_example] at hlim
  exact hlim

end Convergence

end StackedSVD.RankR.Example

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 4. `prop:stacksvd_subspace` at general `r_i` -/

/-- **`prop:stacksvd_subspace`** (`main_paper.tex:834`) at general `r_i`, Layer 1 form: under
the rank-`r` spiked law of the stack at aspect ratio `‖c‖₁`, the stacksvd subspace performance
tends to one rank-one performance per spike of the core matrix.

Route, identical to `prop_stacksvd_subspace` (`RankR/Subspace.lean:498`):
`perfStackRG_eq_sum_spikeVec` turns the performance into the sum over the spike directions;
each term is `SubspaceLawG.align j`; a finite sum of limits in probability is the limit of the
sum (`TendstoInProb.finsum`, `Prob/TendstoInProb.lean`).

No hypothesis beyond the law is needed. In particular the eigenvalues of `C` may repeat, which
is what the paper claims at `main_paper.tex:842`. -/
theorem prop_stacksvd_subspace_general (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (law : m.SubspaceLawG (∑ i, c i)) :
    TendstoInProb μ (fun N ω => m.perfStackRG N ω)
      (limitStackRG (fun i => (m.tbl i).θ) m.R c) := by
  have hfun : (fun N (ω : Ω N) => m.perfStackRG N ω)
      = fun N (ω : Ω N) => ∑ j : Fin r,
          ‖specProjTop (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r
            (m.spikeVecG j N)‖ ^ 2 := by
    funext N ω
    exact m.perfStackRG_eq_sum_spikeVec N ω
  rw [m.limitStackRG_eq c, hfun]
  exact TendstoInProb.finsum (fun j => law.align j)

end UnalignedModelR

end StackedSVD
