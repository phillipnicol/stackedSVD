/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLEMarginal.TableLaw
import StackedSVD.MLEMarginal.LogLik
import StackedSVD.MLEConverse

/-!
# `app:wstacksvd_mle`, the marginalization (L5): the headline statement

`thm_wstacksvd_mle_marginal` bundles three facts at a fixed `N` and `ω`, with
`c_i = n_i / d` as in the paper:

1. the random-effects law of the `M` tables has Lebesgue density `reDensity`
   (`reJointLaw_eq_withDensity`);
2. its logarithm at the observed matrices is `mleLogLik` up to the constant
   `(∑ᵢ n_i) d / 2 · log(2π)` (`reLogLik_eq`);
3. a unit vector maximizes that log-likelihood over unit vectors if and only if it lies in
   the top eigenspace of `stackGramW (optWstack θ c)` (`mleLogLik_max_iff_mem_topSpace`).

Together with `MLE.lean` and `MLEConverse.lean` this closes the audit item L5: the paper's
`ℓ(v)` is now derived from the model, not defined.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped ENNReal Matrix

namespace StackedSVD

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- Clause 3: the marginal MLE over unit vectors is the top eigenspace of the weighted stack
Gram matrix at the paper's weights (`optWstack θ c`, `c_i = n_i / d`). -/
theorem reLogLik_max_iff_mem_topSpace [NeZero M] (m : MultiTableModel μ M n d) (N : ℕ)
    (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
        reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v (fun i => (m.tbl i).X N ω)
          ≤ reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) (WithLp.ofLp x)
              (fun i => (m.tbl i).X N ω))
      ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ)
            (fun i => (n i N : ℝ) / d N)) N ω) (m.isHermitian_stackGramW _ N ω) := by
  have hc : ∀ i, 0 < (n i N : ℝ) / d N :=
    fun i => div_pos (Nat.cast_pos.mpr ((m.tbl i).hn N)) (Nat.cast_pos.mpr ((m.tbl i).hd N))
  refine Iff.trans ?_
    (m.mleLogLik_max_iff_mem_topSpace (fun i => (n i N : ℝ) / d N) hc N ω x hx)
  refine forall_congr' fun v => imp_congr_right fun _ => ?_
  rw [m.reLogLik_eq N ω v, m.reLogLik_eq N ω (WithLp.ofLp x)]
  exact sub_le_sub_iff_right _

/-- `app:wstacksvd_mle`, the marginalization: the random-effects law has density
`reDensity`; its log at the observed data is `mleLogLik` up to a constant; and the unit
maximizers are the top eigenspace of `stackGramW (optWstack θ c)`. -/
theorem thm_wstacksvd_mle_marginal [NeZero M] (m : MultiTableModel μ M n d) (N : ℕ)
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
  ⟨fun v => reJointLaw_eq_withDensity (fun i => (m.tbl i).hn N)
      ((m.tbl ⟨0, Nat.pos_of_ne_zero (NeZero.ne M)⟩).hd N) _ v,
    fun v => m.reLogLik_eq N ω v,
    m.reLogLik_max_iff_mem_topSpace N ω x hx⟩

end MultiTableModel

end StackedSVD
