/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Main
import StackedSVD.RankR.SingleWeight.Example
import StackedSVD.RankR.GeneralGaussian

/-!
# `prop:singleweight_suboptimality`: the witness model and the existence statement

Track E of `notes/archive/singleweight_plan.md` section 4.3, units E1 and E2.

## The paper

`main_paper.tex:915`, `prop:singleweight_suboptimality`: "In the general rank-`r` setting,
there exists a problem instance such that unweighted svdstack outperforms optimally weighted
stacksvd with a single weight per table." The instance is built at `main_paper.tex:2171` to
`:2194`.

## The witness

`M = 2` Gaussian tables, `r = 2`, one spike each (`rk = fun _ => 1`), with

```
R_1 = (1, 0)ᵀ,  R_2 = (0, 1)ᵀ,  θ_0 = 8/5,  c_1 = c_2 = 1,  n_i N = d N = N + 3.
```

Every constant is rational: `θ_0² = 64/25`, `θ_0⁴ = 4096/625 > 1 = c_0`, `β_0² = 39/64`,
`2 β_0² = 39/32`. The equal-weight stacksvd value is `2 betaSq (8/5) 2 = 0.999297753`, below
`2 β_0² = 1.21875`.

`R = Rone (RankR.Example.Rex 0)`, since `Rex 0 = (e_1, (sin 0, cos 0)ᵀ) = (e_1, e_2)`. So the
witness is the paper's own example of `eq:psi_equation` at `sin ψ = 0`
(`main_paper.tex:912`).

## What is stated here

1. `mdl`, the witness model. A real definition, fully proved.
2. `witness_perfRG_tendsto`: unweighted svdstack on the witness tends to `2 β_0²`. Route:
   `prop_general_rank_unweighted_svdstack_general_gaussian` (`RankR/GeneralGaussian.lean:64`),
   then `limitRG_one_eq_limitR` (`RankR/General.lean:583`), `limitR_example_beta`
   (`RankR/Example.lean:922`) and `svdstackEx_zero_s` (`RankR/Example.lean:877`). Both extra
   hypotheses of that theorem are free here: `rtot (fun _ => 1) = 2 = r`, and the top-`r`
   eigengap is vacuous at `r = r̃` by `topGap_of_card_le`.
3. `prop_singleweight_suboptimality_of_law`: the paper's existence claim, with the stacksvd
   limit on the instance taken as a hypothesis. That hypothesis is what the Gaussian discharge
   of `SingleWeightLaw` (Track G of the plan) will supply.

STATUS 2026-09-05: every definition is a real definition; both theorems are proved.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

namespace Witness

/-! ### 1. Two helpers shared by the two witnesses -/

/-- `k` orthonormal columns inside `Fin D`: the first `k` columns of the identity. Copy of
the private helper of `Sat.lean`, which this file cannot see. -/
noncomputable def frame (D k : ℕ) (h : k ≤ D) : Matrix (Fin D) (Fin k) ℝ :=
  (1 : Matrix (Fin D) (Fin D) ℝ).submatrix id (Fin.castLE h)

theorem frame_orth (D k : ℕ) (h : k ≤ D) : (frame D k h)ᵀ * frame D k h = 1 := by
  ext a b
  simp only [frame, Matrix.mul_apply, Matrix.transpose_apply, Matrix.submatrix_apply, id_eq,
    Matrix.one_apply]
  rw [Finset.sum_eq_single (Fin.castLE h a)]
  · by_cases hab : a = b
    · subst hab; simp
    · simp [hab]
  · intro c _ hc
    simp [hc]
  · intro hc
    exact absurd (Finset.mem_univ _) hc

/-- A unit vector of `ℝ^r` gives an `r × 1` matrix with orthonormal columns, which is the
field `hR` of `UnalignedModelR` at `rk = fun _ => 1`. -/
theorem Rone_orth {M r : ℕ} (R : Fin M → EuclideanSpace ℝ (Fin r)) (hR : ∀ i, ‖R i‖ = 1)
    (i : Fin M) : (Rone R i)ᵀ * Rone R i = 1 := by
  have hinner : ((Rone R i)ᵀ * Rone R i) 0 0 = ⟪R i, R i⟫_ℝ := by
    rw [real_inner_eq_dotProduct, Matrix.mul_apply]
    simp [Rone, dotProduct, Matrix.transpose_apply]
  have hnorm : ⟪R i, R i⟫_ℝ = 1 := by
    rw [real_inner_self_eq_norm_sq, hR i]; norm_num
  ext a b
  obtain rfl : a = 0 := Subsingleton.elim a 0
  obtain rfl : b = 0 := Subsingleton.elim b 0
  rw [hinner, hnorm, Matrix.one_apply_eq]

/-- `StrictAnti` is vacuous on a one-element index type. -/
theorem strictAnti_fin_one (f : Fin 1 → ℝ) : StrictAnti f := fun a b hab =>
  absurd (Subsingleton.elim a b) (ne_of_lt hab)

/-! ### 2. The shared probability space: `n_i N = d N = N + 3`, independent Gaussian noise -/

/-- `n_i N = N + 3`, the same for both tables. -/
abbrev nn : Fin 2 → ℕ → ℕ := fun _ N => N + 3

/-- `d N = N + 3`, so `c_i = 1`. -/
abbrev dd : ℕ → ℕ := fun N => N + 3

/-- The sample space at index `N`. -/
abbrev Om (N : ℕ) : Type := (i : Fin 2) → Matrix (Fin (nn i N)) (Fin (dd N)) ℝ

/-- Independent Gaussian noise for the two tables. -/
noncomputable def mu (N : ℕ) : Measure (Om N) :=
  Measure.pi fun i : Fin 2 => gaussianMatrix (nn i N) (dd N)

instance instProbMu (N : ℕ) : IsProbabilityMeasure (mu N) := by
  unfold mu; infer_instance

/-! ### 3. The witness of `prop:singleweight_suboptimality` -/

/-- `θ_0 = 8/5`, the single spike of each table. -/
noncomputable def th : Fin 1 → ℝ := fun _ => 8 / 5

/-- `R_1 = (1, 0)ᵀ`, `R_2 = (0, 1)ᵀ`, as `Rone` of the paper's `Rex 0`
(`main_paper.tex:2171`). -/
noncomputable def Rmat : (i : Fin 2) → Matrix (Fin 2) (Fin ((fun _ : Fin 2 => 1) i)) ℝ :=
  Rone (RankR.Example.Rex 0)

theorem Rmat_orth (i : Fin 2) : (Rmat i)ᵀ * Rmat i = 1 :=
  Rone_orth _ (RankR.Example.norm_Rex 0) i

/-- Table `i` of the witness: one spike of strength `8/5`, right singular direction the
`i`-th column of the shared frame. -/
noncomputable def tbl (i : Fin 2) : SpikedModelR mu (nn i) dd 1 where
  θ := th
  U := fun N => frame (nn i N) 1 (by change (1 : ℕ) ≤ N + 3; omega)
  V := fun N => frame (dd N) 2 (by change (2 : ℕ) ≤ N + 3; omega) * Rmat i
  Z := fun _ ω => ω i
  hθnn := fun _ => by norm_num [th]
  hθanti := strictAnti_fin_one th
  hn := fun N => by change (0 : ℕ) < N + 3; omega
  hd := fun N => by change (0 : ℕ) < N + 3; omega
  hU := fun N => frame_orth _ _ _
  hV := fun N => by
    rw [Matrix.transpose_mul, Matrix.mul_assoc,
      ← Matrix.mul_assoc (frame (dd N) 2 (by change (2 : ℕ) ≤ N + 3; omega))ᵀ, frame_orth,
      Matrix.one_mul, Rmat_orth]
  hZ := fun N => measurable_pi_apply i

/-- The witness model of `prop:singleweight_suboptimality`: `M = 2`, `r = 2`, one spike per
table, `R = Rone (Rex 0)`, `θ_0 = 8/5`, `c_1 = c_2 = 1`. -/
noncomputable def mdl : UnalignedModelR mu 2 nn dd 2 (fun _ => 1) where
  tbl := tbl
  V := fun N => frame (dd N) 2 (by change (2 : ℕ) ≤ N + 3; omega)
  R := Rmat
  hV := fun N => frame_orth _ _ _
  hR := Rmat_orth
  hv := fun i N => rfl

theorem mdl_theta (i : Fin 2) (j : Fin 1) : (mdl.tbl i).θ j = 8 / 5 := rfl

theorem mdl_R : mdl.R = Rone (RankR.Example.Rex 0) := rfl

theorem mdl_regime (i : Fin 2) : (mdl.tbl i).Regime 1 := by
  refine ⟨?_, ?_, ?_⟩
  · exact tendsto_atTop_mono (fun N => Nat.le_add_right N 3) tendsto_id
  · exact tendsto_atTop_mono (fun N => Nat.le_add_right N 3) tendsto_id
  · refine tendsto_const_nhds.congr fun N => ?_
    change (1 : ℝ) = ((N + 3 : ℕ) : ℝ) / ((N + 3 : ℕ) : ℝ)
    have h : ((N + 3 : ℕ) : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
    exact (div_self h).symm

theorem mdl_joint : mdl.JointGaussianNoise := fun _ => ProbabilityTheory.HasLaw.id

/-! ### 4. The unweighted svdstack half -/

/-- Unweighted svdstack on the witness tends in probability to `2 β_0² = 39/32`
(`main_paper.tex:2187`). On this instance `A_{β,R} = I_2` and `(W_opt)^{-2} = (1 - β_0²) I_2`,
so optimally weighted svdstack and unweighted svdstack agree. -/
theorem witness_perfRG_tendsto :
    TendstoInProb mu (fun N ω => mdl.perfRG N ω) (2 * betaSq (8 / 5) 1) := by
  have h := mdl.prop_general_rank_unweighted_svdstack_general_gaussian_one (fun _ => 1)
    (fun _ _ => beta (8 / 5) 1) (fun _ => one_pos) (fun _ _ => rfl) mdl_regime
    (by simp [rtot]) (topGap_of_card_le _ (by simp [rtot])) mdl_joint
  have hlim : limitRG (rk := fun _ : Fin 2 => 1) (fun _ _ => beta (8 / 5) 1) mdl.R
      = 2 * betaSq (8 / 5) 1 := by
    rw [mdl_R, limitRG_one_eq_limitR (fun _ => beta (8 / 5) 1) (RankR.Example.Rex 0),
      RankR.Example.limitR_example_beta one_pos, Real.sin_zero,
      RankR.Example.svdstackEx_zero_s, Scalars.beta_sq]
  rwa [hlim] at h

/-! ### 5. The existence proposition -/

/-- **`prop:singleweight_suboptimality`** (`main_paper.tex:915`), Layer 1 form: there is a
problem instance of the general rank-`r` setting on which unweighted svdstack outperforms
every single weighting of stacksvd.

The witness is `mdl`: two Gaussian tables, `r = 2`, one spike each of strength `θ_0 = 8/5`,
`R_1 = (1,0)ᵀ`, `R_2 = (0,1)ᵀ`, `c_1 = c_2 = 1`. Unweighted svdstack reaches
`2 β_0² = 39/32`; every weighting `w` with both weights positive reaches only
`swLimitEx (8/5) 1 w`, which is strictly below it; and the equal-weight stacksvd value
`2 betaSq (8/5) 2` is below it as well.

`hlaw` is the stacksvd limit on the instance. It is the conclusion of
`prop_gen_rank_stacksvd_singleweight` once `SingleWeightLaw` is discharged for Gaussian noise
(Track G of `notes/archive/singleweight_plan.md`), which is why it is a hypothesis here and not a
theorem.

Scope. `0 < w i` is not free: `w_i = 0` is fine for the formula but it kills the rank
condition on this instance (plan scan 3.6, modeling choice 5). -/
theorem prop_singleweight_suboptimality_of_law
    (hlaw : ∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
      TendstoInProb mu (fun N ω => mdl.perfSW w N ω) (swLimitEx (8 / 5) 1 w)) :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : UnalignedModelR μ 2 n d 2 (fun _ => 1)),
      (∀ i j, (m.tbl i).θ j = 8 / 5) ∧ (∀ i, (m.tbl i).Regime 1) ∧
      m.JointGaussianNoise ∧ m.R = Rone (RankR.Example.Rex 0) ∧
      TendstoInProb μ (fun N ω => m.perfRG N ω) (2 * betaSq (8 / 5) 1) ∧
      (∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
        TendstoInProb μ (fun N ω => m.perfSW w N ω) (swLimitEx (8 / 5) 1 w) ∧
        swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1) ∧
      2 * betaSq (8 / 5) 2 < 2 * betaSq (8 / 5) 1 := by
  refine ⟨Om, inferInstance, mu, inferInstance, nn, dd, mdl, mdl_theta, mdl_regime, mdl_joint,
    mdl_R, witness_perfRG_tendsto, fun w h0 h1 => ⟨hlaw w h0 h1, swLimitEx_lt w h0 h1⟩, ?_⟩
  norm_num [betaSq]

end Witness

end SingleWeight

end StackedSVD
