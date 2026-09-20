/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Simplicity
import StackedSVD.RankR.RMT.SimplicityR
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.RankR.RMT.EdgeGlueR
import StackedSVD.RankR.StackGamma
import StackedSVD.RankR.Het.Scalars
import StackedSVD.RankR.SubspaceGStack

/-!
# Stage E5: simplicity of the weighted stack Gram matrix at every index below `r`

Step 7 of `notes/archive/rankr_TrackE_plan.md` section 2.3, the E5 row of section 3, and change 4 of
`notes/archive/audit_rankr_plan_E_2026-09-02.md`. This file discharges the `simpleIdxJ` field of
`HeteroLawR` for jointly Gaussian noise.

## Route

The file crosses two existing arguments.

1. `RMT/Het/Simplicity.lean` proves top simplicity for the row-scaled affine family
   `Y ↦ A + rowScale s Y` by the resultant of the characteristic polynomial. Every
   polynomial of that file is reused here with no change: `gramRowPolyRight`,
   `gramRowPolyLeft`, and the two `_ne_zero` witnesses.
2. `RankR/RMT/SimplicityR.lean` proves `SimpleSpec` at any rank `rk ≤ min p q`, in both
   regimes. The regime `q ≤ p` uses full injectivity of the sorted spectrum; the regime
   `p ≤ q` uses `simpleSpec_p_of_gram_left`, which handles the repeated zero eigenvalue of
   the `q - p` dimensional null space.

`simpleSpec_ae_rowScale` is the crossing: `simpleSpec_ae_affine`
(`RankR/RMT/SimplicityAffineR.lean:76`) with a per-row scale `s` in place of one scalar `t`.

Section 2 puts the weighted stack in that form. `stackXW_eq_addRowScale` writes
`X_stack(w) = rowScale (w ∘ blk) (A Vᵀ) + rowScale (d^{-1/2} w ∘ blk) Z_stack`, the rank-`r`
twin of `MultiTableModel.stackW_X_eq_addRowScale`. The law of `Z_stack` is
`UnalignedModelR.hasLaw_stackZG`, so `HasLaw.ae_iff` with the measurable predicate
`measurable_simpleSpec_rowScale` transports the Gaussian statement to `μ N`.

## Zero weights are not removed here

The rank-1 file needs a restriction step, because its `topSimple` field allows `w_i = 0` in
some tables. Every weight of the rank-`r` stack is `w_ij = θ_ij/√(θ_ij² + c_i)` with
`θ_ij > 0` and `c_i > 0`, so `wStackR_thetaAligned_ne_zero` gives `∀ i, w i ≠ 0` directly and
the sub-collection argument is not needed.

## What E6 imports

`simpleSpec_ae_rowScale` (audit change 6): the rank-one projector step of stage E6 reads it
at the block-scaled Gaussian noise.

STATUS 2026-09-02: `lean-local.sh` exit 0; 0 `sorry`; three standard axioms only.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### 1. The row-scaled affine family, at every rank below `min p q` -/

section RowScaleR

variable {p q : ℕ}

/-- **Item S at every index below `rk`, with row scales** (plan step 7). For every fixed shift
`A` and every row scale `s` with no zero entry, the Gram matrix of `A + rowScale s Y` has its
top `rk` eigenvalues pairwise distinct from every other eigenvalue, for almost every Gaussian
`Y`, whenever `rk ≤ min p q`. Rank-1 mirror: `topSimple_ae_rowScale`
(`RMT/Het/Simplicity.lean:143`); unit-scale mirror: `simpleSpec_ae_affine`
(`RankR/RMT/SimplicityAffineR.lean:76`).

The regime `p ≤ q` is the one where the `q × q` Gram matrix carries a zero eigenvalue of
multiplicity `q - p`. `simpleSpec_p_of_gram_left` still gives `SimpleSpec` at `rk = p`,
because that predicate asks only that each of the top `rk` eigenvalues differs from every
other one, and the top `p` eigenvalues are positive. -/
theorem simpleSpec_ae_rowScale (hp : 0 < p) (hq : 0 < q) (A : Matrix (Fin p) (Fin q) ℝ)
    {s : Fin p → ℝ} (hs : ∀ a, s a ≠ 0) (rk : ℕ) (hrk : rk ≤ min p q) :
    ∀ᵐ Y ∂(gaussianMatrix p q),
      SimpleSpec ((A + rowScale s Y)ᵀ * (A + rowScale s Y))
        (isHermitian_transpose_mul_self (A + rowScale s Y)) rk := by
  rcases le_total q p with hqp | hpq
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramRowPolyRight A s)
      (gramRowPolyRight_ne_zero hq hqp A hs)] with Y hY
    rw [eval_gramRowPolyRight] at hY
    exact simpleSpec_of_injective _
      (injective_eigenvalues₀_of_separable _ ((charRes_ne_zero_iff_separable hq _).mp hY)) rk
  · filter_upwards [ae_eval_ne_zero_gaussianMatrix (gramRowPolyLeft A s)
      (gramRowPolyLeft_ne_zero hp hpq A hs)] with Y hY
    rw [eval_gramRowPolyLeft] at hY
    have hinjB : Function.Injective
        (isHermitian_mul_transpose_self (A + rowScale s Y)).eigenvalues₀ :=
      injective_eigenvalues₀_of_separable _
        ((charRes_ne_zero_iff_separable hp _).mp (left_ne_zero_of_mul hY))
    exact simpleSpec_mono (hrk.trans (min_le_left _ _))
      (simpleSpec_p_of_gram_left hpq (A + rowScale s Y) hinjB (right_ne_zero_of_mul hY))

/-- Measurability of the row-scaled simplicity predicate. The row-scale twin of
`RankRStack.measurable_simpleSpec_affine` (`RankR/RMT/EdgeGlueR.lean:583`), built from the same
`EdgeGlueR.measurableSet_simpleSpec_gram` and from `measurable_addRowScale`. -/
theorem measurable_simpleSpec_rowScale (A : Matrix (Fin p) (Fin q) ℝ) (s : Fin p → ℝ)
    (rk : ℕ) :
    Measurable fun Y : Matrix (Fin p) (Fin q) ℝ =>
      SimpleSpec ((A + rowScale s Y)ᵀ * (A + rowScale s Y))
        (isHermitian_transpose_mul_self (A + rowScale s Y)) rk := by
  rw [← measurableSet_setOfPred]
  have hpre : {Y : Matrix (Fin p) (Fin q) ℝ |
      SimpleSpec ((A + rowScale s Y)ᵀ * (A + rowScale s Y))
        (isHermitian_transpose_mul_self (A + rowScale s Y)) rk}
      = (fun Y => A + rowScale s Y) ⁻¹'
        {Z : Matrix (Fin p) (Fin q) ℝ |
          SimpleSpec (Zᵀ * Z) (isHermitian_transpose_mul_self Z) rk} := rfl
  rw [hpre]
  exact measurable_addRowScale A s (EdgeGlueR.measurableSet_simpleSpec_gram rk)

end RowScaleR

/-! ### 2. From `SimpleSpec` to `SimpleIdx` -/

/-- `SimpleSpec` at `rk` gives `SimpleIdx` at every index below `rk`. `SimpleSpec`
(`LinAlg/SpecIdx.lean:88`) reads every index `k < rk`; `SimpleIdx`
(`LinAlg/SpecIdxPerturb.lean:484`) reads the single index `j`. -/
theorem simpleIdx_of_simpleSpec {p : ℕ} {A : Matrix (Fin p) (Fin p) ℝ} {hA : A.IsHermitian}
    {rk j : ℕ} (hs : SimpleSpec A hA rk) (hj : j < rk) : SimpleIdx A hA j := by
  intro k l hk hkl
  exact hs k l (by omega) hkl

/-! ### 3. The weighted stack of an `UnalignedModelR` -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The weighted stack is the affine family `Y ↦ A + rowScale s Y` at the stacked noise, with
`A = rowScale (w ∘ blk) (A Vᵀ)` and `s q = d^{-1/2} w_{blk q}`. Rank-1 mirror:
`MultiTableModel.stackW_X_eq_addRowScale` (`RMT/Het/Simplicity.lean:188`). -/
theorem stackXW_eq_addRowScale (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    m.stackXW w N ω
      = rowScale (fun q => w (finSigmaFinEquiv.symm q).1) (m.signalPartG N)
        + rowScale (fun q => (Real.sqrt (d N))⁻¹ * w (finSigmaFinEquiv.symm q).1)
            (m.stackZG N ω) := by
  ext q k
  have hX : m.stackXG N ω q k
      = m.signalPartG N q k + (Real.sqrt (d N))⁻¹ * m.stackZG N ω q k := by
    rw [m.stackX_eqG N ω, Matrix.add_apply, m.stackE_eqG N ω, Matrix.smul_apply, smul_eq_mul]
  rw [Matrix.add_apply, rowScale_apply, rowScale_apply, m.stackXW_apply' w N ω q k,
    ← m.stackXG_apply' N ω q k, hX]
  ring

/-- **The weighted stack Gram matrix has pairwise distinct top-`k` eigenvalues almost surely,
when every weight is nonzero.** Rank-1 mirror: `MultiTableModel.topSimple_ae_stackGramW`
(`RMT/Het/Simplicity.lean:203`). The two size conditions come from the model: `0 < d N` is
`SpikedModelR.hd` of table `0`, and `0 < ∑ i, n i N` is `SpikedModelR.hn` of table `0`, so
neither `0 < r` nor `0 < n 0 N` has to be assumed. -/
theorem simpleSpec_ae_stackGramW [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (hG : m.JointGaussianNoise) (hw : ∀ i, w i ≠ 0) (k : ℕ) (N : ℕ)
    (hk : k ≤ min (∑ i, n i N) (d N)) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) k := by
  classical
  have hdN : 0 < d N := (m.tbl ⟨0, NeZero.pos M⟩).hd N
  have hnN : 0 < ∑ i, n i N :=
    lt_of_lt_of_le ((m.tbl ⟨0, NeZero.pos M⟩).hn N)
      (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
        (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hdN
  set A : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
    rowScale (fun q => w (finSigmaFinEquiv.symm q).1) (m.signalPartG N) with hAdef
  set s : Fin (∑ i, n i N) → ℝ :=
    fun q => (Real.sqrt (d N))⁻¹ * w (finSigmaFinEquiv.symm q).1 with hsdef
  have hs : ∀ a, s a ≠ 0 := fun a => by
    rw [hsdef]
    exact mul_ne_zero (inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne') (hw _)
  have hae : ∀ᵐ ω ∂(μ N),
      SimpleSpec ((A + rowScale s (m.stackZG N ω))ᵀ * (A + rowScale s (m.stackZG N ω)))
        (isHermitian_transpose_mul_self (A + rowScale s (m.stackZG N ω))) k :=
    ((m.hasLaw_stackZG hG N).ae_iff (measurable_simpleSpec_rowScale A s k)).mpr
      (simpleSpec_ae_rowScale hnN hdN A hs k hk)
  filter_upwards [hae] with ω hω
  refine EdgeGlueR.simpleSpec_congr ?_ _ (m.isHermitian_stackGramW w N ω) hω
  have hXω : m.stackXW w N ω = A + rowScale s (m.stackZG N ω) := by
    rw [hAdef, hsdef]
    exact m.stackXW_eq_addRowScale w N ω
  rw [← hXω]
  rfl

/-- The stack weights `w_ij = θ_ij/√(θ_ij² + c_i)` are nonzero: `θ_ij > 0` at the weighting
component `j` by the hypothesis `hposj` (F8, 2026-09-05; inside the class only the last
component can be zero, `SpikedModelR.θ_pos_of_ne_last`) and `c_i > 0` by hypothesis. This is
`Scalars.wStackR_ne_zero` read at the model strength table `m.thetaAligned`. -/
theorem wStackR_thetaAligned_ne_zero (m : UnalignedModelR μ M n d r (alignedRk M r))
    {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) {j : Fin r} (hposj : ∀ i, 0 < (m.tbl i).θ j)
    (i : Fin M) : Scalars.wStackR m.thetaAligned c j i ≠ 0 :=
  Scalars.wStackR_ne_zero hc hposj i

/-- **The `simpleIdxJ` field of `HeteroLawR` for Gaussian noise** (audit change 4). At the
`j`-th weighting every weight is nonzero, so `simpleSpec_ae_stackGramW` holds at `k = r`, and
`ℓ_j < r` (`Scalars.ellR_lt`) puts the index the estimator reads inside that range. Paper:
the well-defined estimator `v̂_{j,stacksvd}` of `thm:rank_r_stacksvd` (`main_paper.tex:2337`);
the rank-one mirror is the simple top eigenvalue of `thm:stacksvd_weighted` (`:463`). -/
theorem simpleIdxJ_of_gaussian [NeZero M] (m : UnalignedModelR μ M n d r (alignedRk M r))
    {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hG : m.JointGaussianNoise) {j : Fin r}
    (hposj : ∀ i, 0 < (m.tbl i).θ j) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j) := by
  classical
  have hrn : r ≤ n ⟨0, NeZero.pos M⟩ N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_n N
  have hsum : r ≤ ∑ i, n i N :=
    hrn.trans (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
      (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  have hrd : r ≤ d N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_d N
  have hae := m.simpleSpec_ae_stackGramW (Scalars.wStackR m.thetaAligned c j) hG
    (fun i => m.wStackR_thetaAligned_ne_zero hc hposj i) r N (le_min hsum hrd)
  filter_upwards [hae] with ω hω
  exact simpleIdx_of_simpleSpec hω (Scalars.ellR_lt m.thetaAligned c j)

end UnalignedModelR

end StackedSVD
