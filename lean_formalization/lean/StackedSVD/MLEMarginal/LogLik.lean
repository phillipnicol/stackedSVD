/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.MLEMarginal.Defs

/-!
# The marginal log-likelihood is `mleLogLik` (L5, unit U7)

`reLogLik_eq`: at the observed matrices `(m.tbl i).X N ω`,
`log reDensity = mleLogLik (c_i = n_i / d) - (∑ᵢ n_i) d / 2 · log(2π)`. This is the
paper's "up to constants independent of `v`" (`main_paper.tex:1660`), made explicit. The
identity holds for every `v`; the unit norm enters only in `MLEMarginal/Main.lean`.
-/

open MeasureTheory ProbabilityTheory Matrix
open scoped Matrix

namespace StackedSVD

/-- `∑ₖ xₖᵀ S xₖ = tr(S Xᵀ X)` for the rows `xₖ` of `X` (`Matrix.trace_mul_comm`). -/
theorem sum_dotProduct_mulVec_eq_trace {n d : ℕ} (S : Matrix (Fin d) (Fin d) ℝ)
    (X : Matrix (Fin n) (Fin d) ℝ) :
    ∑ k, X k ⬝ᵥ S *ᵥ X k = (S * (Xᵀ * X)).trace := by
  have key : (X * (S * Xᵀ)).trace = ∑ k, X k ⬝ᵥ S *ᵥ X k := by
    simp [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
      dotProduct, Matrix.mulVec_apply_eq_sum]
  rw [← key, Matrix.trace_mul_comm, Matrix.mul_assoc]

/-- The logarithm of `gaussDensity` for `0 < det S`. -/
theorem log_gaussDensity {d : ℕ} {S : Matrix (Fin d) (Fin d) ℝ} (hS : 0 < S.det)
    (x : Fin d → ℝ) :
    Real.log (gaussDensity d S x)
      = -(x ⬝ᵥ S⁻¹ *ᵥ x) / 2 - ((d : ℝ) * Real.log (2 * Real.pi) + Real.log S.det) / 2 := by
  have hpi : (0 : ℝ) < (2 * Real.pi) ^ d * S.det := by positivity
  rw [gaussDensity, Real.log_div (Real.exp_ne_zero _) (Real.sqrt_ne_zero'.mpr hpi),
    Real.log_exp, Real.log_sqrt hpi.le, Real.log_mul (by positivity) hS.ne', Real.log_pow]

/-- Per-table identity behind `reLogLik_eq`: the row log-densities of one table sum to
`-(1/2)(n log det S + tr(S⁻¹ XᵀX)) - n d / 2 · log(2π)`. -/
theorem sum_log_gaussDensity {n d : ℕ} {S : Matrix (Fin d) (Fin d) ℝ} (hS : 0 < S.det)
    (X : Matrix (Fin n) (Fin d) ℝ) :
    ∑ k, Real.log (gaussDensity d S (X k))
      = -(1 / 2) * ((n : ℝ) * Real.log S.det + (S⁻¹ * (Xᵀ * X)).trace)
        - (n : ℝ) * d / 2 * Real.log (2 * Real.pi) := by
  have hlog : ∀ k, Real.log (gaussDensity d S (X k))
      = -(1 / 2) * (X k ⬝ᵥ S⁻¹ *ᵥ X k)
        + -(1 / 2) * ((d : ℝ) * Real.log (2 * Real.pi) + Real.log S.det) := fun k => by
    rw [log_gaussDensity hS (X k)]; ring
  rw [Finset.sum_congr rfl (fun k _ => hlog k), Finset.sum_add_distrib, ← Finset.mul_sum,
    Finset.sum_const, Finset.card_fin, nsmul_eq_mul, sum_dotProduct_mulVec_eq_trace]
  ring

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- U7: the marginal log-likelihood at the observed data is `mleLogLik` at `c_i = n_i / d`,
minus the constant `(∑ᵢ n_i) d / 2 · log(2π)`. -/
theorem reLogLik_eq [NeZero M] (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N)
    (v : Fin (d N) → ℝ) :
    reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v (fun i => (m.tbl i).X N ω)
      = m.mleLogLik (fun i => (n i N : ℝ) / d N) N ω v
        - ((∑ i, (n i N : ℝ)) * d N / 2) * Real.log (2 * Real.pi) := by
  have hSdet : ∀ i : Fin M,
      0 < (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v).det :=
    fun i => det_mleCov_pos ((m.tbl i).hd N) (by positivity) v
  have hgpos : ∀ (i : Fin M) (k : Fin (n i N)),
      0 < gaussDensity (d N) (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v)
        ((m.tbl i).X N ω k) :=
    fun i k => gaussDensity_pos (hSdet i) _
  have hprodpos : ∀ i : Fin M,
      0 < ∏ k, gaussDensity (d N) (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v)
        ((m.tbl i).X N ω k) :=
    fun i => Finset.prod_pos (fun k _ => hgpos i k)
  have hstep1 : reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v
      (fun i => (m.tbl i).X N ω)
      = ∑ i : Fin M, ∑ k, Real.log (gaussDensity (d N)
          (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v)
          ((m.tbl i).X N ω k)) := by
    rw [reLogLik, reDensity, Real.log_prod (fun i _ => (hprodpos i).ne')]
    exact Finset.sum_congr rfl (fun i _ => Real.log_prod (fun k _ => (hgpos i k).ne'))
  have hstep2 : ∀ i : Fin M, ∑ k, Real.log (gaussDensity (d N)
      (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v) ((m.tbl i).X N ω k))
      = -(1 / 2) * ((n i N : ℝ) * Real.log
            (mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v).det
          + ((mleCov (d N) ((m.tbl i).θ ^ 2 / ((n i N : ℝ) / d N)) v)⁻¹
              * (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)).trace)
        + (n i N : ℝ) * (-(d N : ℝ) / 2 * Real.log (2 * Real.pi)) := fun i => by
    rw [sum_log_gaussDensity (hSdet i) ((m.tbl i).X N ω)]
    ring
  rw [hstep1, Finset.sum_congr rfl (fun i _ => hstep2 i), Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.sum_mul, mleLogLik]
  ring

end MultiTableModel

end StackedSVD
