/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.Align
import StackedSVD.RankR.Het.Bulk
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.Het.Sub
import StackedSVD.RankR.StackMain

/-!
# Stage E9: `HeteroLawR` for Gaussian noise, and the corollary with no hypothesis structure

Track E, stage E9 (`notes/archive/rankr_TrackE_plan.md`, step 8). This file closes the track. It
discharges the hypothesis structure `HeteroLawR` (`RankR/StackGamma.lean:489`) on the exactly
aligned rank-`r` family with Gaussian noise, and it restates the paper's corollary
`thm:rank_r_stacksvd` (`main_paper.tex:2337`) with no hypothesis structure left.

## The paper

`main_paper.tex:2337`. On the exactly aligned rank-`r` model, rank-`r` weighted stacksvd
satisfies

```
(v_jᵀ v̂_{j,stacksvd})² →p γ_j     for j = 1, …, r,
‖Vᵀ V̂_stacksvd‖_F²     →p ∑_j γ_j.
```

The paper's separation hypothesis `θ̃_jj ≠ θ̃_jk` (`eq:stacksvd_apptildTheta`) is **not** a
hypothesis below. Inside this model class `SpikedModelR.hθanti` orders the strengths strictly
in every table, so the separation holds by `Scalars.rhoHet_lt_of_lt`
(`RankR/Het/Scalars.lean:384`). The two hypotheses that do appear are `0 < c_i` and
`R_i = 1`; the second one is the exactly aligned family (`notes/FLAGGED.md` items 16 and 17).

## The three fields

Write `w_{·j} = wStackR θ c j` for the weights of stack `j` and `ℓ_j = ellR θ c j` for the
sorted index the estimator reads.

* `align` and `crossProj` are one statement: the overlap of column `l` of the signal frame
  with the `ℓ_j`-th right singular subspace of `X_stack^(j)` tends to `γ_j` when `l = j` and
  to `0` when `l ≠ j`. That statement is `key_of_gaussian` (`RankR/Het/Sub.lean`), whose
  proof splits on `eq:assumption4` at component `j`:
  - supercritical: stage E8a (`align_cross_het_of_gaussian_sup`) at the index `ellSup`, which
    is `ℓ_j` inside the class (`Scalars.ellSup_eq_ellR`), with limit
    `Scalars.Lw θ_j c w_{·j} = γ_j` (`Scalars.Lw_wStackR_eq_gammaR`);
  - subcritical: `γ_j = 0` (`Scalars.gammaR_eq_zero_of_not_assumption4`) and every index at or
    above `numSup` is a bulk index, with `numSup ≤ ℓ_j`
    (`Scalars.numSup_le_of_not_assumption4`), so stage E8b (`cross_het_of_gaussian_bulk`)
    gives the limit `0` at both `l = j` and `l ≠ j`.
* `simpleIdxJ` is stage E5 (`simpleIdxJ_of_gaussian`, `RankR/Het/Simplicity.lean:206`), read
  through the same drop step (`simpleIdxJ_of_gaussian_of_exists`).

F8 (2026-09-05). The tables need not carry every component: `SpikedModelR.hθnn` allows
`θ_ij = 0` for the last components of a table. The endpoints assume `hnz`, that at least one
table carries each component `j`. When `hnz` fails at `j`, every weight `wStackR … j i` is
`0`, so `X_stack^(j) = 0` and the paper's own side condition `θ̃_jj ≠ θ̃_jk`
(`main_paper.tex:2338`, with `θ̃_jk = ∑_i w_ij² θ_ik²`) fails as well (for `r ≥ 2`; at
`r = 1` the side condition is empty, and the paper's `γ_1` is undefined, since its equation
reads `0 = 1`); `hnz` is therefore implied by the paper's hypothesis and adds nothing. The
proof drops the tables with `θ_ij = 0` through the sub-model of `RankR/Het/Sub.lean` (route
D of `notes/archive/F8_zero_spikes.md`).

Stages E8a and E8b both carry the side condition `d N = p N + r`, so both run on the shifted
model `m.shift k₀` of stage E3 and `HeteroLawR.of_shift` (`RankR/Het/Edge.lean:360`) carries
the two limit fields back. `simpleIdxJ` is a statement about one `N` at a time, so it is
supplied for the original model.

## Contents

1. `heteroLawR_of_gaussian_aux`: the law under the side condition `d N = p N + r`. Its two
   limit fields are `key_of_gaussian` of `RankR/Het/Sub.lean`.
2. `heteroLawR_of_gaussian`: the law with no side condition.
3. `thm_rank_r_stacksvd_gaussian` and the three single displays.

No `sorry`, no `axiom`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! ### 1. The law under the side condition of stages E8a and E8b -/

/-- **`HeteroLawR` for Gaussian noise, with the side condition `d N = p N + r`.** This is the
form stages E8a and E8b take; `heteroLawR_of_gaussian` removes the side condition by the tail
shift of stage E3.

The two limit fields come from one statement (`key_of_gaussian`, `RankR/Het/Sub.lean`) that
reads the overlap of column `l` with the `ℓ_j`-th right singular subspace of `X_stack^(j)`:
the limit is `γ_j` at `l = j` and `0` at `l ≠ j`. The case split is `eq:assumption4` at
component `j`, at the weights of component `j`. Above the threshold stage E8a applies at the
index `ellSup`, which the class identifies with `ℓ_j`; below it `γ_j = 0` and `ℓ_j` is at or
above `numSup`, so stage E8b applies at every column. `simpleIdxJ` is stage E5.

`hnz` says that at least one table carries each component (F8). The tables with `θ_ij = 0`
carry weight `0` at component `j` and are dropped by the sub-model of `RankR/Het/Sub.lean`. -/
theorem heteroLawR_of_gaussian_aux [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    m.HeteroLawR c := by
  refine ⟨fun j => ?_, fun j l hl => ?_,
    fun j N => m.simpleIdxJ_of_gaussian_of_exists hc hG (hnz j) N⟩
  · have h := m.key_of_gaussian c hc hR hreg hG hpd (hnz j) j
    rw [if_pos rfl] at h
    exact h
  · have h := m.key_of_gaussian c hc hR hreg hG hpd (hnz j) l
    rw [if_neg hl] at h
    exact h

/-! ### 2. The law -/

/-- **`HeteroLawR` for Gaussian noise** (the last stage of Track E). On the exactly aligned
rank-`r` family with jointly Gaussian noise, positive limit ratios `c_i` and `R_i = 1`, the
three fields of the hypothesis structure of `thm:rank_r_stacksvd` hold.

The side condition `d N = p N + r` of stages E8a and E8b is met after a tail shift
(`exists_shift_add`, `RankR/Het/Edge.lean:310`, from `d N → ∞`), and `HeteroLawR.of_shift`
(`RankR/Het/Edge.lean:360`) carries the two limit fields back to the original index. The
third field is a statement about one `N` at a time, so stage E5 supplies it directly for
`m`. The rank-one mirror is the pair `MultiTableModel.HeteroLaw` and its Gaussian
discharge.

Scope. `SpikedModelR.hθnn` and `hθanti` order the spikes strictly in every table and keep
them nonnegative. On this family the paper's index is `ℓ_j = j` (`main_paper.tex:2369`) at
every component that `hnz` covers. The paper (`main_paper.tex:2306`, with the
condition and the conclusion at `:2368` and `:2369`) allows an unordered `θ`, where `ℓ_j` can
differ from the outlier index; that case needs a permutation `R_i`, or a model class without
`hθanti` (`notes/FLAGGED.md` item 17 (7)). More than generality is at stake there: with the
paper's own `ℓ_j` the first display of `thm:rank_r_stacksvd` is false on an explicit
`M = 2`, `r = 2` instance, where index 0 carries `0.571 ± 0.008` and the paper's index 1
carries `0.011 ± 0.003` at `d = 1200` (`paper_edits.md` finding E7). -/
theorem heteroLawR_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    m.HeteroLawR c := by
  obtain ⟨k₀, p, hpd, -⟩ := m.exists_shift_add (hreg 0).2.1
  have hprob : ∀ N, IsProbabilityMeasure (μ (N + k₀)) := isProbabilityMeasure_shift k₀
  refine HeteroLawR.of_shift k₀ (heteroLawR_of_gaussian_aux (m.shift k₀) c hc hR
    (fun i => m.shift_Regime k₀ i (hreg i)) (m.shift_JointGaussianNoise k₀ hG) hpd hnz)
    (fun j N => m.simpleIdxJ_of_gaussian_of_exists hc hG (hnz j) N)

/-! ### 3. The corollary of the paper, with no hypothesis structure -/

/-- The first display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) for Gaussian noise,
in projector form: the squared overlap of `v_j` with the `ℓ_j`-th right singular subspace of
`X_stack^(j)` tends to `γ_j`. This form carries no simplicity side condition. Layer 1 mirror:
`thm_rank_r_stacksvd_proj` (`RankR/StackMain.lean:155`). -/
theorem thm_rank_r_stacksvd_proj_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) (j : Fin r) :
    TendstoInProb μ (fun N ω => m.stackOverlapJ c j N ω)
      (Scalars.gammaR m.thetaAligned c j) :=
  m.thm_rank_r_stacksvd_proj c (m.heteroLawR_of_gaussian c hc hR hreg hG hnz) j

/-- The first display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) for Gaussian noise,
in the paper's inner-product form: `(v_jᵀ v̂_{j,stacksvd})² →p γ_j`. Layer 1 mirror:
`thm_rank_r_stacksvd_inner` (`RankR/StackMain.lean:166`). -/
theorem thm_rank_r_stacksvd_inner_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) (j : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
      (Scalars.gammaR m.thetaAligned c j) :=
  m.thm_rank_r_stacksvd_inner c (m.heteroLawR_of_gaussian c hc hR hreg hG hnz) j

/-- The second display of **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) for Gaussian
noise: `‖Vᵀ V̂_stacksvd‖_F² →p ∑_j γ_j`. Layer 1 mirror: `thm_rank_r_stacksvd_frobenius`
(`RankR/StackMain.lean:181`). -/
theorem thm_rank_r_stacksvd_frobenius_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
      (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) :=
  m.thm_rank_r_stacksvd_frobenius c (m.heteroLawR_of_gaussian c hc hR hreg hG hnz)

/-- **`thm:rank_r_stacksvd`** (`main_paper.tex:2337`) for Gaussian noise, with no hypothesis
structure: on the exactly aligned rank-`r` family, rank-`r` weighted stacksvd satisfies both
displays of the corollary,

```
(v_jᵀ v̂_{j,stacksvd})² →p γ_j     for every j,
‖Vᵀ V̂_stacksvd‖_F²     →p ∑_j γ_j.
```

Where the paper's hypotheses sit. `assum:rank_r` is the model class
`UnalignedModelR … (alignedRk M r)` with `R_i = 1` (`hR`), and `assum:general_noise` is
`hG` together with the per-table regime `hreg` and `0 < c_i` (`hc`). The separation
`θ̃_jj ≠ θ̃_jk` of the corollary is not a hypothesis: `SpikedModelR.hθanti` orders the
strengths strictly in every table, so it holds inside the class
(`Scalars.rhoHet_lt_of_lt`). See `notes/FLAGGED.md` items 16 and 17. `hnz` asks that at least
one table carries each component `j` (F8, 2026-09-05). It is implied by the paper's own
separation: at an all-zero component every weight `w_ij` is `0`, so `θ̃_jj = θ̃_jk = 0`.

Scope. `SpikedModelR.hθnn` and `hθanti` order the spikes strictly in every table and keep
them nonnegative. On this family the paper's index is `ℓ_j = j` (`main_paper.tex:2369`) at
every component that `hnz` covers. The paper (`main_paper.tex:2306`, with the
condition and the conclusion at `:2368` and `:2369`) allows an unordered `θ`, where `ℓ_j` can
differ from the outlier index; that case is not an `UnalignedModelR` with `R_i = 1`. It needs
a permutation `R_i`, or a model class without `hθanti` (`Scalars.ellR` is stated generally so
that the statement survives the widening).
More than generality is at stake there: with the paper's own `ℓ_j` the first display of
`thm:rank_r_stacksvd` is false on an explicit `M = 2`, `r = 2` instance, where index 0
carries `0.571 ± 0.008` and the paper's index 1 carries `0.011 ± 0.003` at `d = 1200`
(`paper_edits.md` finding E7). See `notes/FLAGGED.md` item 17 (7) and
`RankR/StackGamma.lean:304` (section 3 header). -/
theorem thm_rank_r_stacksvd_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    (∀ j : Fin r, TendstoInProb μ
        (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR m.thetaAligned c j)) ∧
      TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
        (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j) :=
  m.thm_rank_r_stacksvd c (m.heteroLawR_of_gaussian c hc hR hreg hG hnz)

end UnalignedModelR

end StackedSVD
