/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Regimes
import StackedSVD.RankR.SingleWeight.Het.Sub
import StackedSVD.RankR.SubspaceGaussian
import StackedSVD.RankR.SubspaceMain

/-!
# The tie of the single-weight witness (F18b, unit U2)

At `w 0 = w 1 = a` the single-weight stack is `a` times the unweighted stack, so its
performance is the unweighted subspace performance `perfStackRG` on the event that the top
two eigenvalues are simple, and `prop_stacksvd_subspace_general_gaussian` gives the limit.
`hlaw_tie` is the tie case of `hlaw`. Plan: `notes/archive/F18b_plan.md`, section 2 (U2) and
section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

section Tie

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- U2.1: the core matrix at one spike per table is the rank-one core matrix. -/
theorem CBlock_one_eq_Cmat {M r : ℕ} (θ : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    CBlock (rk := fun _ : Fin M => 1) (fun i _ => θ i) (Rone R) = Cmat θ R := by
  ext k l
  rw [cblock_apply, cmat_apply]
  exact Finset.sum_congr rfl fun i _ => by simp [Rone]

/-- `Matrix.IsHermitian.eigenvalues` does not see the proof term, so it transports along an
equality of matrices. Mirror: `eigenvalues₀_congr_mat` (`LinAlg/Eigen.lean`). -/
private theorem eigenvalues_congr_mat {p : ℕ} {A B : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) (hB : B.IsHermitian) (h : A = B) :
    hA.eigenvalues = hB.eigenvalues := by
  subst h
  rfl

/-- U2.2: the stacksvd limit at one spike per table is the rank-one one. -/
theorem limitStackRG_one_eq_limitStackR {M r : ℕ} (θ : Fin M → ℝ)
    (R : Fin M → EuclideanSpace ℝ (Fin r)) (c : Fin M → ℝ) :
    limitStackRG (rk := fun _ : Fin M => 1) (fun i _ => θ i) (Rone R) c
      = limitStackR θ R c := by
  have h : (isHermitian_CBlock (rk := fun _ : Fin M => 1) (fun i _ => θ i)
        (Rone R)).eigenvalues = (isHermitian_Cmat θ R).eigenvalues :=
    eigenvalues_congr_mat _ _ (CBlock_one_eq_Cmat θ R)
  simp only [limitStackRG, limitStackR, h]

/-- U2.3: under `SimpleSpec` the top-`r` projector splits into the `r` index projectors. -/
theorem normSq_specProjTop_eq_sum_specProjIdx {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) {r : ℕ} (hsimple : SimpleSpec A hA r) (hrp : r ≤ p)
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProjTop A hA r x‖ ^ 2 = ∑ l : Fin r, ‖specProjIdx A hA (l : ℕ) x‖ ^ 2 := by
  rw [normSq_specProjTop_eq_sum_of_simpleSpec hsimple hrp x,
    Fin.sum_univ_eq_sum_range (fun k => ‖specProjIdx A hA k x‖ ^ 2) r]

/-- U2.4: at equal weights the single-weight performance is the unweighted one. -/
theorem UnalignedModelR.perfSW_const_eq_perfStackRG (m : UnalignedModelR μ M n d r rk)
    {a : ℝ} (ha : 0 < a) (N : ℕ) (ω : Ω N) (hrd : r ≤ d N)
    (hsimple : SimpleSpec (m.stackGramW (fun _ => a) N ω)
      (m.isHermitian_stackGramW (fun _ => a) N ω) r) :
    m.perfSW (fun _ => a) N ω = m.perfStackRG N ω := by
  have ha2 : (0 : ℝ) < a ^ 2 := by positivity
  -- 1. the weighted stack is `a` times the unweighted stack
  have hX : m.stackXW (fun _ => a) N ω = a • m.stackXG N ω := by
    ext q k
    rw [m.stackXW_apply' (fun _ => a) N ω q k, Matrix.smul_apply, m.stackXG_apply' N ω q k,
      smul_eq_mul]
  -- 2. so the Gram matrices differ by the positive factor `a ^ 2`
  have hG : m.stackGramW (fun _ => a) N ω = (a ^ 2) • m.stackGramG N ω := by
    simp only [UnalignedModelR.stackGramW, UnalignedModelR.stackGramG, hX]
    rw [Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, pow_two]
  have hHt : ((a ^ 2) • m.stackGramG N ω).IsHermitian := by
    rw [← hG]
    exact m.isHermitian_stackGramW (fun _ => a) N ω
  -- 3. the factor drops out of every index projector
  have hproj : ∀ (l : ℕ) (y : EuclideanSpace ℝ (Fin (d N))),
      overlapIdx (m.stackXW (fun _ => a) N ω) l y
        = ‖specProjIdx (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) l y‖ ^ 2 := by
    intro l y
    have h1 : overlapIdx (m.stackXW (fun _ => a) N ω) l y
        = ‖specProjIdx (m.stackGramW (fun _ => a) N ω)
            (m.isHermitian_stackGramW (fun _ => a) N ω) l y‖ ^ 2 := rfl
    rw [h1, specProjIdx_congr_mat (m.isHermitian_stackGramW (fun _ => a) N ω) hHt hG l,
      EdgeGlueDetR.specProjIdx_smul (m.isHermitian_stackGramG N ω) hHt ha2 l]
  -- and out of the simplicity hypothesis
  have hinv : (a ^ 2)⁻¹ • m.stackGramW (fun _ => a) N ω = m.stackGramG N ω := by
    rw [hG, smul_smul, inv_mul_cancel₀ ha2.ne', one_smul]
  have hHinv : ((a ^ 2)⁻¹ • m.stackGramW (fun _ => a) N ω).IsHermitian := by
    rw [hinv]
    exact m.isHermitian_stackGramG N ω
  have hsimpleG : SimpleSpec (m.stackGramG N ω) (m.isHermitian_stackGramG N ω) r :=
    EdgeGlueR.simpleSpec_congr hinv hHinv (m.isHermitian_stackGramG N ω)
      (HetDeloc.simpleSpec_smul (m.isHermitian_stackGramW (fun _ => a) N ω) hHinv
        (by positivity) hsimple)
  -- 4. sum the split of U2.3 over the columns of `V`
  simp only [UnalignedModelR.perfSW, UnalignedModelR.perfStackRG]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [normSq_specProjTop_eq_sum_specProjIdx (m.isHermitian_stackGramG N ω) hsimpleG hrd]
  exact Finset.sum_congr rfl fun l _ => hproj (l : ℕ) (m.colVecG N k)

end Tie

namespace SingleWeight

namespace Witness

/-- U2.5: the law at a tie, from `prop_stacksvd_subspace_general_gaussian`. -/
theorem hlaw_tie {a : ℝ} (ha : 0 < a) :
    TendstoInProb mu (fun N ω => mdl.perfSW (fun _ => a) N ω)
      (swLimitEx (8 / 5) 1 (fun _ => a)) := by
  -- the two limits agree
  have hlim : swLimitEx (8 / 5) 1 (fun _ => a)
      = limitStackRG (fun i => (mdl.tbl i).θ) mdl.R (fun _ : Fin 2 => (1 : ℝ)) := by
    have h2 : limitStackRG (fun i => (mdl.tbl i).θ) mdl.R (fun _ : Fin 2 => (1 : ℝ))
        = limitStackR (fun _ : Fin 2 => (8 / 5 : ℝ)) (RankR.Example.Rex 0)
            (fun _ : Fin 2 => (1 : ℝ)) :=
      limitStackRG_one_eq_limitStackR (fun _ : Fin 2 => (8 / 5 : ℝ)) (RankR.Example.Rex 0)
        (fun _ => 1)
    rw [h2, RankR.Example.limitStackR_example, Real.sin_zero,
      swLimitEx_tie (by norm_num) one_pos ha (by norm_num)]
  -- the unweighted subspace law
  have hstack : TendstoInProb mu (fun N ω => mdl.perfStackRG N ω)
      (limitStackRG (fun i => (mdl.tbl i).θ) mdl.R (fun _ : Fin 2 => (1 : ℝ))) :=
    mdl.prop_stacksvd_subspace_general_gaussian mdl_joint (fun _ => 1) mdl_regime (by norm_num)
  -- the two performances agree off a null event
  have hzero : ∀ N, mu N {ω | mdl.perfSW (fun _ => a) N ω ≠ mdl.perfStackRG N ω} = 0 := by
    intro N
    have hk : (2 : ℕ) ≤ min (∑ i with (fun _ : Fin 2 => a) i ≠ 0, nn i N) (dd N) := by
      refine le_min ?_ (by change (2 : ℕ) ≤ N + 3; omega)
      refine le_trans (show (2 : ℕ) ≤ nn 0 N by change (2 : ℕ) ≤ N + 3; omega)
        (Finset.single_le_sum (f := fun i : Fin 2 => nn i N) (fun i _ => Nat.zero_le _) ?_)
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ 0, ha.ne'⟩
    have hae := UnalignedModelR.simpleSpec_ae_stackGramW_of_exists mdl (fun _ => a) mdl_joint
      ⟨0, ha.ne'⟩ 2 N hk
    refine measure_mono_null ?_ (ae_iff.mp hae)
    intro ω hω hs
    exact hω (mdl.perfSW_const_eq_perfStackRG ha N ω (by change (2 : ℕ) ≤ N + 3; omega) hs)
  rw [hlim]
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hstack
  exact tendsto_const_nhds.congr fun N => (hzero N).symm

end Witness

end SingleWeight

end StackedSVD
