/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.GeneralFrob
import StackedSVD.RankR.WeightedUpperG
import StackedSVD.RankR.AlignedMain
import StackedSVD.RankR.RMT.TableLawGaussian

/-!
# Section 7 under Gaussian noise, at general `r_i`

Task C10 of `notes/archive/rankr_plan_C.md`, second half. Every Layer 1 theorem of Section 7 takes
the black box `law : ∀ i, (m.tbl i).TableLawR (c i)`. Task C10 discharges that black box for
Gaussian noise at every `r_i` (`UnalignedModelR.tableLawR_of_gaussian_rk`,
`RankR/RMT/TableLawGaussian.lean`), so each theorem gains a twin that replaces `law` and `hG :
m.IndepNoise` by the proportional regime of each table and the joint Gaussian law.

This file states those twins. The pattern of every proof is one line:
`m.<name> ... (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise`.

The `_gaussian_one` block of `RankR/GeneralFrob.lean:498-580` and
`RankR/WeightedUpperG.lean:560` states the same corollaries at one spike per table
(`rk = fun _ => 1`). Those statements stay where they are: they are the `rk = 1` instances,
and they route through the rank-1 chain `singleTableLaw_of_gaussian`, not through Track C.

Contents, in the order of the originals:

1. `prop_general_rank_unweighted_svdstack_general_gaussian` (`RankR/GeneralMain.lean:603`)
2. `prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian`
   (`RankR/GeneralFrob.lean:271`)
3. `thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian`
   (`RankR/GeneralFrob.lean:419`)
4. `thm_gen_rank_weight_svdstak_general_r_of_rank_gaussian` (`RankR/GeneralFrob.lean:456`)
5. `thm_gen_rank_weight_svdstak_general_r_paper_gaussian` (`RankR/GeneralFrob.lean:478`)
6. `thm_gen_rank_weight_svdstak_general_r_full_gaussian` (`RankR/WeightedUpperG.lean:447`)
7. `thm_rank_r_svdstack_aggregate_gaussian` (`RankR/WeightedUpperG.lean:638`)
8. `thm_rank_r_svdstack_component_gaussian` and
   `thm_rank_r_svdstack_component_of_model_gaussian` (`RankR/AlignedMain.lean`, Track D
   item D4)

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology

open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. `prop:general_rank_unweighted_svdstack` -/

/-- **`prop:general_rank_unweighted_svdstack` under Gaussian noise, at general `r_i`**
(`main_paper.tex:799`). No `TableLawR` hypothesis: `hreg` and `hG` produce it through
`UnalignedModelR.tableLawR_of_gaussian_rk` (task C10). The `rk = 1` instance is
`prop_general_rank_unweighted_svdstack_general_gaussian_one`
(`RankR/GeneralFrob.lean:498`). -/
theorem prop_general_rank_unweighted_svdstack_general_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hrr : r ≤ rtot rk)
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRG N ω) (limitRG β m.R) :=
  m.prop_general_rank_unweighted_svdstack_general c β hc hβdef hrr hgap
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-- The same on the paper's own quantity `‖V̂_svdstackᵀ V‖_F²`, at the canonical top-`r`
eigenframe of `Ṽ Ṽᵀ` and under Gaussian noise, at general `r_i`; no `TableLawR` hypothesis.
The `rk = 1` instance is
`prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian_one`. -/
theorem prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk) (c : Fin M → ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hrp : r ≤ Fintype.card (Fin (rtot rk)))
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackG N ω
        (topEigMat (m.isHermitian_gramG N ω) hrp)
        (topEigVal (m.isHermitian_gramG N ω) hrp))ᵀ * m.V N))
      (limitRG β m.R) :=
  m.prop_general_rank_unweighted_svdstack_general_frobenius_eig c β hc hβdef hrp hgap
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-! ### 2. `thm:gen_rank_weight_svdstak` -/

/-- **`thm:gen_rank_weight_svdstak`, first display** (`main_paper.tex:893`) under Gaussian
noise, at general `r_i`, on the paper's own quantity `‖V̂_svdstack(W⋆)ᵀ V‖_F²` and at the
canonical frame of `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ`. The `rk = 1` instance is
`thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian_one`. -/
theorem thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r)
    (hrp : r ≤ Fintype.card (Fin (rtot rk))) (hrankB : (BBlock β m.R).rank = r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => frobSq ((m.vhatSvdstackGW (optWG β) N ω
        (topEigMat (m.isHermitian_gramWG (optWG β) N ω) hrp)
        (topEigVal (m.isHermitian_gramWG (optWG β) N ω) hrp))ᵀ * m.V N))
      (limitOptG β m.R hrp) :=
  m.thm_gen_rank_weight_svdstak_general_r_frobenius_eig c β hc hβdef hr hrp hrankB
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-- **Finding C1 under Gaussian noise, at general `r_i`**: the two conclusions of
`thm:gen_rank_weight_svdstak` with no `hrr : r ≤ r̃` (the bound sits inside `hrankB`) and no
`TableLawR` hypothesis. -/
theorem thm_gen_rank_weight_svdstak_general_r_of_rank_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r) (hrankB : (BBlock β m.R).rank = r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB)) ∧
      ∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by
            have h2 := rank_le_card_of_rank_BBlock hrankB
            omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R (rank_le_card_of_rank_BBlock hrankB) :=
  m.thm_gen_rank_weight_svdstak_general_r_of_rank c β hc hβdef hr hrankB
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-- **Finding C3, the paper-literal corollary, under Gaussian noise at general `r_i`**: the
paper's own two hypotheses (`main_paper.tex:757, 893`), `β_ij > 0` at every spike and
`Rank(∑_i R_i R_iᵀ) = r`, plus the proportional regime and the joint Gaussian law. -/
theorem thm_gen_rank_weight_svdstak_general_r_paper_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r) (hβpos : ∀ i j, 0 < β i j)
    (hrank : (∑ i, m.R i * (m.R i)ᵀ).rank = r) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
        (limitOptG β m.R
          (rank_le_card_of_rank_BBlock (rank_BBlock_of_paper m.R hβpos hrank))) ∧
      ∀ W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ,
        TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r →
        0 < (isHermitian_ABlockW W β m.R).eigenvalues₀
          ⟨r - 1, by
            have h2 := rank_le_card_of_rank_BBlock
              (rank_BBlock_of_paper (β := β) m.R hβpos hrank)
            omega⟩ →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R
            (rank_le_card_of_rank_BBlock (rank_BBlock_of_paper m.R hβpos hrank)) :=
  m.thm_gen_rank_weight_svdstak_general_r_paper c β hc hβdef hr hβpos hrank
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-- **Item E2 under Gaussian noise, at general `r_i`.** The three conclusions of
`thm_gen_rank_weight_svdstak_general_r_full` with no hypothesis on `rank B_R` and no
`TableLawR` hypothesis: the performance at `W⋆` tends to `L⋆`, every admissible `W` tends to a
limit at most `L⋆`, and with probability tending to one no weight matrix beats `L⋆ + ε`. The
`rk = 1` instance is `thm_gen_rank_weight_svdstak_general_r_full_gaussian_one`
(`RankR/WeightedUpperG.lean:560`). -/
theorem thm_gen_rank_weight_svdstak_general_r_full_gaussian
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
  m.thm_gen_rank_weight_svdstak_general_r_full c β hc hβdef hr hrr
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

end UnalignedModelR

/-! ### 3. Track D item D3: the aggregate clause -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- **`thm:rank_r_svdstack`, aggregate clause, under Gaussian noise** (`main_paper.tex:2401`,
Track D item D3). Under the exactly aligned rank-`r` model with Gaussian noise and the
proportional regime, the weighted svdstack performance at `W⋆` converges in probability to
`∑_j S_j / (S_j + 1)`, `S_j = ∑_i β_ij² / (1 - β_ij²)`. No `TableLawR` hypothesis. -/
theorem thm_rank_r_svdstack_aggregate_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
      (∑ j, Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_aggregate c β hc hβdef hR hM
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise

/-! ### 4. Track D item D4: the component clause -/

/-- **`thm:rank_r_svdstack`, component clause, under Gaussian noise** (`main_paper.tex:2405`,
Track D item D4). `v̂_j` is column `j` of `V̂_svdstack(W⋆)` at the canonical top-`r` eigenframe
of `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ`. `hS` is the paper's hypothesis (b) `S_j ≠ S_k`, in the sorted form; the
model makes it redundant, which is `thm_rank_r_svdstack_component_of_model_gaussian`. No
`TableLawR` hypothesis. -/
theorem thm_rank_r_svdstack_component_gaussian [∀ N, IsProbabilityMeasure (μ N)]
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
  m.thm_rank_r_svdstack_component_eig c β hc hβdef hR hM hS
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise j

/-- The same with the paper's hypothesis (b) dropped: inside the model the spike order
`SpikedModelR.hθanti` already separates `Sagg` wherever the component clause has content
(`UnalignedModelR.saggSep_of_model`). -/
theorem thm_rank_r_svdstack_component_of_model_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_component_of_model c β hc hβdef hR hM
    (m.tableLawR_of_gaussian_rk hc hreg hG) hG.indepNoise j

end UnalignedModelR

end StackedSVD
