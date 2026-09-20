/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.MPhet
import StackedSVD.RMT.Het.R5het
import StackedSVD.RMT.Het.R6het
import StackedSVD.RankR.StackGamma

/-!
# Stage E0: the heteroscedastic scalars of the rank-`r` weighted stack

Stage E0 of `notes/archive/rankr_TrackE_plan.md`, restated by change 1 of
`notes/archive/audit_rankr_plan_E_2026-09-02.md`. Everything here is real analysis on finite sums
over `Fin M`; no probability appears.

## The paper

Track E writes the `j`-th weighted stack `X_stack^(j)` of `eq:stacksvd_appXstack`
(`main_paper.tex:2325`) as one heteroscedastic rank-`r` spiked matrix. The `r` spikes decouple,
so component `k` under the weighting `w = w_{·j}` is the rank-1 heteroscedastic problem at the
block energy profile `(w_i² θ_ik²)_i`. Its outlier is `MPhet.rhoHet (θ_{·k}) c w` (`rho_k`), its
detectability condition is `Scalars.Assumption4 (θ_{·k}) c w` (`eq:assumption4`,
`main_paper.tex:1368`), and its overlap limit is `Scalars.Lw (θ_{·k}) c w`.

`alg:rank_r_stacksvd` line `alg_line:stacksvd_rank` (`main_paper.tex:2382`) reads the estimator
off the sorted index `ℓ_j`, the rank of `θ̃_jj` among `{θ̃_jk}`. Section 4.2 of the Track E plan
shows that at unordered `θ` this is **not** the index of the outlier of spike `j`. The index the
random matrix chain is true at is `ellSup`, the rank of `rho_j` among the outliers of the
supercritical spikes. Inside the model class (`hθnn`, `hθanti`) the two agree, which is
`ellSup_eq_ellR` below.

## Contents

### `namespace MPhet` (rank-1 facts over a strength vector `θ : Fin M → ℝ`)

1. `secular_le_secular_of_le`, `secular_lt_secular_of_lt`: the secular function of
   `lem:secular_equation` falls when the energy profile rises, above `max_i w_i²`.
2. `exists_isGammaTop_of_le`, `gammaTop_le_gammaTop`, `gammaTop_lt_gammaTop`: the top root
   `γ₁` of `eq:assumption4` is monotone in `θ²`.
3. `assumption4_mono`: `eq:assumption4` survives a rise of the energy profile.
4. `rhoHet_le_rhoHet`, `rhoHet_lt_rhoHet`: the outlier `rho` is monotone in `θ²`. This is
   fact 1 of step 4 of the Track E plan.
5. `nuHet = F'(rho)`, `nuHet_pos`, `overlap_identity_nuHet` (`1/(rho * nu) = L(w)`, the paper's
   `main_paper.tex:1433`).

### `namespace Scalars` (over a strength table `θ : Fin M → Fin r → ℝ`)

6. `numSup`, `ellSup` and their range bounds.
7. The column-`j` forms at the stack weights `wStackR`: `assumption4_wStackR_iff`,
   `Lw_wStackR_eq_gammaR`, `gammaR_eq_zero_of_not_assumption4`, `wStackR_ne_zero`.
8. The prefix property inside the model class: `assumption4_of_lt`, `rhoHet_lt_of_lt`.
9. `ellSup_eq_ellR` and `numSup_le_of_not_assumption4`, the two index theorems that stages E8
   and E9 consume.

## Already in the tree, so not restated here

The sign of `1 + F` around the outlier is item 5 of the E0 brief. Both halves exist already and
carry exactly the stated hypotheses, so this file cites them and adds nothing:

* `MPhet.one_add_F_pos` (`RMT/Het/R5het.lean:234`): `0 < 1 + Fhet θ c w z` for
  `rhoHet θ c w < z`, under `hc` and `Assumption4 θ c w`.
* `MPhet.one_add_F_neg` (`RMT/Het/R5het.lean:216`): `1 + Fhet θ c w z < 0` for
  `bHet c w < z < rhoHet θ c w`, under the same two hypotheses.
* `MPhet.one_add_F_pos_of_not_assumption4` (`RMT/Het/R6het.lean:140`) covers the subcritical
  spike on all of `(bHet, ∞)`.

`MPhet.FhetDeriv_rhoHet_pos` (`RMT/Het/R5het.lean:264`) is likewise the content of `nuHet_pos`.

No `sorry`, no new axiom.
-/

open Filter Topology

namespace StackedSVD

namespace MPhet

open Scalars

variable {M : ℕ} {θ θ' c w : Fin M → ℝ}

/-! ### 1. The secular function falls when the energy profile rises -/

/-- One summand of `f` above `max_i w_i²`: the denominator `w_j² - λ` is negative, so a larger
`θ_j²` gives a smaller (more negative) term. -/
private theorem secular_term_le_of_le (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2) {lam : ℝ}
    (hlam : wSqMax w < lam) (j : Fin M) :
    θ' j ^ 2 * w j ^ 2 / (w j ^ 2 - lam) ≤ θ j ^ 2 * w j ^ 2 / (w j ^ 2 - lam) := by
  have hd : (0:ℝ) < lam - w j ^ 2 := by
    have := le_wSqMax w j
    linarith
  have hw2 : (0:ℝ) ≤ w j ^ 2 := sq_nonneg _
  rw [show w j ^ 2 - lam = -(lam - w j ^ 2) by ring, div_neg, div_neg, neg_le_neg_iff,
    div_le_div_iff₀ hd hd]
  nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.mpr (hθ j)) hw2) hd.le]

/-- A larger energy profile keeps a table with a nonzero product `θ_i w_i`, which is what
`secular_strictMonoOn` and `existsUnique_gammaTop` read. -/
theorem exists_signal_of_le (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2) (hex : ∃ g, IsGammaTop θ w g) :
    ∃ j, θ' j * w j ≠ 0 := by
  obtain ⟨j0, hj0⟩ := exists_signal_of_root (gammaTop_of_exists hex)
  refine ⟨j0, ?_⟩
  have h1 : 0 < θ j0 ^ 2 * w j0 ^ 2 := by
    rw [show θ j0 ^ 2 * w j0 ^ 2 = (θ j0 * w j0) ^ 2 by ring]
    exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hj0))
  have h2 : 0 < θ' j0 ^ 2 * w j0 ^ 2 := by nlinarith [hθ j0, sq_nonneg (w j0)]
  intro h0
  rw [show θ' j0 ^ 2 * w j0 ^ 2 = (θ' j0 * w j0) ^ 2 by ring, h0] at h2
  norm_num at h2

/-- Above `max_i w_i²` every summand of `f(λ) = 1 + ∑_j θ_j² w_j²/(w_j² - λ)`
(`lem:secular_equation`) is nonpositive and falls when `θ_j²` rises. -/
theorem secular_le_secular_of_le (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2) {lam : ℝ}
    (hlam : wSqMax w < lam) : secular θ' w lam ≤ secular θ w lam := by
  have h := Finset.sum_le_sum
    (fun j (_ : j ∈ (Finset.univ : Finset (Fin M))) => secular_term_le_of_le hθ hlam j)
  unfold secular
  linarith

/-- The strict twin: one table with a nonzero weight and a strictly larger `θ_i²` makes the
inequality strict. -/
theorem secular_lt_secular_of_lt (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (hlt : ∃ i, w i ≠ 0 ∧ θ i ^ 2 < θ' i ^ 2) {lam : ℝ} (hlam : wSqMax w < lam) :
    secular θ' w lam < secular θ w lam := by
  obtain ⟨i0, hw0, hlt0⟩ := hlt
  have hd : (0:ℝ) < lam - w i0 ^ 2 := by
    have := le_wSqMax w i0
    linarith
  have hw2 : (0:ℝ) < w i0 ^ 2 := by positivity
  have hstrict : θ' i0 ^ 2 * w i0 ^ 2 / (w i0 ^ 2 - lam)
      < θ i0 ^ 2 * w i0 ^ 2 / (w i0 ^ 2 - lam) := by
    rw [show w i0 ^ 2 - lam = -(lam - w i0 ^ 2) by ring, div_neg, div_neg, neg_lt_neg_iff,
      div_lt_div_iff₀ hd hd]
    nlinarith [mul_pos (mul_pos (sub_pos.mpr hlt0) hw2) hd]
  have h := Finset.sum_lt_sum
    (fun j (_ : j ∈ (Finset.univ : Finset (Fin M))) => secular_term_le_of_le hθ hlam j)
    ⟨i0, Finset.mem_univ _, hstrict⟩
  unfold secular
  linarith

/-- `f ≥ 0` far above the block variances. This is the step `hfb` inside the proof of
`Scalars.existsUnique_gammaTop`, which that proof does not export. -/
theorem secular_nonneg_of_far {b : ℝ}
    (hb : ∀ j, (∑ l, θ l ^ 2 * w l ^ 2) + 1 ≤ b - w j ^ 2) : 0 ≤ secular θ w b := by
  set S := ∑ l, θ l ^ 2 * w l ^ 2 with hSdef
  have hS0 : (0:ℝ) ≤ S := by
    rw [hSdef]
    exact Finset.sum_nonneg fun l _ => by positivity
  have hS1 : (0:ℝ) < S + 1 := by linarith
  have hterm : ∀ j : Fin M,
      0 ≤ θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1) := by
    intro j
    have hle : S + 1 ≤ b - w j ^ 2 := hb j
    have h := div_le_div_of_nonneg_left (a := θ j ^ 2 * w j ^ 2) (by positivity) hS1 hle
    rw [show w j ^ 2 - b = -(b - w j ^ 2) by ring, div_neg]
    linarith
  have hpos : 0 ≤ ∑ j, (θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1)) :=
    Finset.sum_nonneg fun j _ => hterm j
  rw [Finset.sum_add_distrib, ← Finset.sum_div, ← hSdef] at hpos
  have hSS : S / (S + 1) ≤ 1 := by
    rw [div_le_one hS1]
    linarith
  unfold secular
  linarith

/-! ### 2. `γ₁` is monotone in the energy profile -/

/-- A larger energy profile still has an outlier root. `f_{θ'} ≤ 0` at `γ₁(θ)` and `f_{θ'} ≥ 0`
far out, so the intermediate value theorem gives a root above `max_i w_i²`. -/
theorem exists_isGammaTop_of_le (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (hex : ∃ g, IsGammaTop θ w g) : ∃ g, IsGammaTop θ' w g := by
  have ha : wSqMax w < gammaTop θ w := (gammaTop_of_exists hex).1
  have hga : (0:ℝ) ≤ gammaTop θ w := le_trans (wSqMax_nonneg w) ha.le
  have hS0 : (0:ℝ) ≤ ∑ l, θ' l ^ 2 * w l ^ 2 := Finset.sum_nonneg fun l _ => by positivity
  have hfa : secular θ' w (gammaTop θ w) ≤ 0 := by
    have h := secular_le_secular_of_le hθ ha
    rw [(gammaTop_of_exists hex).2] at h
    exact h
  set B := wSqMax w + (∑ l, θ' l ^ 2 * w l ^ 2) + 1 + gammaTop θ w with hBdef
  have hfB : 0 ≤ secular θ' w B := by
    refine secular_nonneg_of_far fun j => ?_
    have hj := le_wSqMax w j
    rw [hBdef]
    linarith
  have haB : gammaTop θ w ≤ B := by
    have hW := wSqMax_nonneg w
    rw [hBdef]
    linarith
  obtain ⟨g, hgmem, hgval⟩ :=
    intermediate_value_Icc haB (secular_continuousOn (θ := θ') (b := B) ha) ⟨hfa, hfB⟩
  exact ⟨g, lt_of_lt_of_le ha hgmem.1, hgval⟩

/-- `γ₁` rises with the energy profile: `f_{θ'}(γ₁(θ)) ≤ 0 = f_{θ'}(γ₁(θ'))` and `f_{θ'}` is
strictly increasing above `max_i w_i²`. -/
theorem gammaTop_le_gammaTop (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2) (hex : ∃ g, IsGammaTop θ w g) :
    gammaTop θ w ≤ gammaTop θ' w := by
  have hex' := exists_isGammaTop_of_le hθ hex
  have ha : wSqMax w < gammaTop θ w := (gammaTop_of_exists hex).1
  have ha' : wSqMax w < gammaTop θ' w := (gammaTop_of_exists hex').1
  have hle : secular θ' w (gammaTop θ w) ≤ secular θ' w (gammaTop θ' w) := by
    rw [(gammaTop_of_exists hex').2]
    have h := secular_le_secular_of_le hθ ha
    rw [(gammaTop_of_exists hex).2] at h
    exact h
  exact ((secular_strictMonoOn (exists_signal_of_le hθ hex)).le_iff_le
    (Set.mem_Ioi.2 ha) (Set.mem_Ioi.2 ha')).mp hle

/-- The strict twin of `gammaTop_le_gammaTop`. -/
theorem gammaTop_lt_gammaTop (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (hlt : ∃ i, w i ≠ 0 ∧ θ i ^ 2 < θ' i ^ 2) (hex : ∃ g, IsGammaTop θ w g) :
    gammaTop θ w < gammaTop θ' w := by
  have hex' := exists_isGammaTop_of_le hθ hex
  have ha : wSqMax w < gammaTop θ w := (gammaTop_of_exists hex).1
  have ha' : wSqMax w < gammaTop θ' w := (gammaTop_of_exists hex').1
  have hlt' : secular θ' w (gammaTop θ w) < secular θ' w (gammaTop θ' w) := by
    rw [(gammaTop_of_exists hex').2]
    have h := secular_lt_secular_of_lt hθ hlt ha
    rw [(gammaTop_of_exists hex).2] at h
    exact h
  exact ((secular_strictMonoOn (exists_signal_of_le hθ hex)).lt_iff_lt
    (Set.mem_Ioi.2 ha) (Set.mem_Ioi.2 ha')).mp hlt'

/-! ### 3. `eq:assumption4` survives a rise of the energy profile -/

/-- `eq:assumption4` is monotone in `θ²`: the sum `∑_i c_i w_i⁴/(γ₁ - w_i²)²` falls when `γ₁`
rises, and `γ₁` rises with the profile. -/
theorem assumption4_mono (hc : ∀ i, 0 < c i) (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (h4 : Assumption4 θ c w) : Assumption4 θ' c w := by
  have hex' := exists_isGammaTop_of_le hθ h4.1
  have hmono : gammaTop θ w ≤ gammaTop θ' w := gammaTop_le_gammaTop hθ h4.1
  refine ⟨hex', lt_of_le_of_lt (Finset.sum_le_sum fun i _ => ?_) h4.2⟩
  have hd : 0 < gammaTop θ w - w i ^ 2 := gammaTop_sub_pos h4.1 i
  have hd' : 0 < gammaTop θ' w - w i ^ 2 := gammaTop_sub_pos hex' i
  have hnum : (0:ℝ) ≤ c i * w i ^ 4 := by
    have := (hc i).le
    positivity
  have hsq : (gammaTop θ w - w i ^ 2) ^ 2 ≤ (gammaTop θ' w - w i ^ 2) ^ 2 := by nlinarith
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hsq hnum]

/-! ### 4. The outlier is monotone in the energy profile -/

/-- Fact 1 of step 4 of the Track E plan: the outlier `rho` rises with the energy profile.
Route: `rhoHet = zfun (-1/γ₁)`, `zfun` is strictly increasing on `[s⋆, 0)`, and
`-1/γ₁(θ) ≤ -1/γ₁(θ')` with both points on the physical branch. -/
theorem rhoHet_le_rhoHet (hc : ∀ i, 0 < c i) (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (h4 : Assumption4 θ c w) : rhoHet θ c w ≤ rhoHet θ' c w := by
  have h4' := assumption4_mono hc hθ h4
  have hw := exists_w_ne_zero_of_root h4.1
  have hb := (assumption4_iff_branch hc).mp h4
  have hb' := (assumption4_iff_branch hc).mp h4'
  have hs : -1 / gammaTop θ w ≤ -1 / gammaTop θ' w := by
    rw [neg_div, neg_div, neg_le_neg_iff]
    exact one_div_le_one_div_of_le (gammaTop_pos h4.1) (gammaTop_le_gammaTop hθ h4.1)
  rw [rhoHet_eq_zfun h4.1, rhoHet_eq_zfun h4'.1]
  exact (zfun_strictMonoOn hc hw).monotoneOn ⟨hb.1.le, hb.2⟩ ⟨hb'.1.le, hb'.2⟩ hs

/-- The strict twin of `rhoHet_le_rhoHet`. -/
theorem rhoHet_lt_rhoHet (hc : ∀ i, 0 < c i) (hθ : ∀ i, θ i ^ 2 ≤ θ' i ^ 2)
    (hlt : ∃ i, w i ≠ 0 ∧ θ i ^ 2 < θ' i ^ 2) (h4 : Assumption4 θ c w) :
    rhoHet θ c w < rhoHet θ' c w := by
  have h4' := assumption4_mono hc hθ h4
  have hw := exists_w_ne_zero_of_root h4.1
  have hb := (assumption4_iff_branch hc).mp h4
  have hb' := (assumption4_iff_branch hc).mp h4'
  have hs : -1 / gammaTop θ w < -1 / gammaTop θ' w := by
    rw [neg_div, neg_div, neg_lt_neg_iff]
    exact one_div_lt_one_div_of_lt (gammaTop_pos h4.1) (gammaTop_lt_gammaTop hθ hlt h4.1)
  rw [rhoHet_eq_zfun h4.1, rhoHet_eq_zfun h4'.1]
  exact zfun_strictMonoOn hc hw ⟨hb.1.le, hb.2⟩ ⟨hb'.1.le, hb'.2⟩ hs

/-! ### 5. `nu = F'(rho)` -/

/-- `ν = F'(ρ)`: the derivative of the limiting secular function `1 + F` at the outlier. The
overlap limit of `main_paper.tex:1433` is `1/(ρ ν)`. -/
noncomputable def nuHet (θ c w : Fin M → ℝ) : ℝ := FhetDeriv θ c w (rhoHet θ c w)

/-- `0 < ν` under `eq:assumption4`. This is `FhetDeriv_rhoHet_pos` read on `nuHet`. -/
theorem nuHet_pos (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) : 0 < nuHet θ c w :=
  FhetDeriv_rhoHet_pos hc h4

/-- The overlap identity `1/(ρ ν) = L(w)` (`main_paper.tex:1433`), read on `nuHet`. -/
theorem overlap_identity_nuHet (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    1 / (rhoHet θ c w * nuHet θ c w) = Lw θ c w :=
  overlap_identity_het hc h4

end MPhet

namespace Scalars

variable {M r : ℕ}

/-! ### 6. The two counts over the components -/

open Classical in
/-- The number of components that are supercritical at the weights `w`. -/
noncomputable def numSup (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (w : Fin M → ℝ) : ℕ :=
  (Finset.univ.filter fun l : Fin r => Assumption4 (fun i => θ i l) c w).card

open Classical in
/-- The sorted index of the outlier of component `k` at the weights `w`: the number of
supercritical components whose outlier location `rho_l` is strictly above `rho_k`. This is the
index the Track E chain is true at; `ellR` is the paper's `ℓ_j` and agrees with it inside the
model class (`ellSup_eq_ellR`). -/
noncomputable def ellSup (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (w : Fin M → ℝ) (k : Fin r) :
    ℕ :=
  (Finset.univ.filter fun l : Fin r =>
    Assumption4 (fun i => θ i l) c w ∧
      MPhet.rhoHet (fun i => θ i k) c w < MPhet.rhoHet (fun i => θ i l) c w).card

/-- At most `r` components are supercritical. -/
theorem numSup_le (θ : Fin M → Fin r → ℝ) (c w : Fin M → ℝ) : numSup θ c w ≤ r := by
  classical
  rw [numSup]
  refine le_trans (Finset.card_le_card ?_) (le_of_eq (Finset.card_fin r))
  exact Finset.filter_subset _ _

/-- The outlier index counts a subset of the supercritical components. -/
theorem ellSup_le_numSup (θ : Fin M → Fin r → ℝ) (c w : Fin M → ℝ) (k : Fin r) :
    ellSup θ c w k ≤ numSup θ c w := by
  classical
  rw [ellSup, numSup]
  refine Finset.card_le_card fun l hl => ?_
  have h := Finset.mem_filter.mp hl
  exact Finset.mem_filter.mpr ⟨h.1, h.2.1⟩

/-- `ellSup < r`: the component `k` itself never enters its own filter. Mirror: `ellR_lt`. -/
theorem ellSup_lt (θ : Fin M → Fin r → ℝ) (c w : Fin M → ℝ) (k : Fin r) :
    ellSup θ c w k < r := by
  classical
  rw [ellSup]
  refine lt_of_lt_of_le (Finset.card_lt_card ⟨Finset.filter_subset _ _, fun hsub => ?_⟩)
    (le_of_eq (Finset.card_fin r))
  have hk := Finset.mem_filter.mp (hsub (Finset.mem_univ k))
  exact absurd hk.2.2 (lt_irrefl _)

/-! ### 7. Column `j` at the stack weights -/

variable {θ : Fin M → Fin r → ℝ} {c w : Fin M → ℝ}

/-- `eq:assumption4` at the stack weights of column `j` is the paper's detectability threshold
`∑_i θ_ij⁴/c_i > 1` (`main_paper.tex:1602`). -/
theorem assumption4_wStackR_iff (hc : ∀ i, 0 < c i) (j : Fin r) :
    Assumption4 (fun i => θ i j) c (wStackR θ c j) ↔ 1 < ∑ i, θ i j ^ 4 / c i :=
  assumption4_optW_iff hc

/-- `L(w_{·j}) = γ_j` (`eq:stacksvd_gammak`): the rank-1 optimum read on column `j`. -/
theorem Lw_wStackR_eq_gammaR (hc : ∀ i, 0 < c i) (j : Fin r) :
    Lw (fun i => θ i j) c (wStackR θ c j) = gammaR θ c j :=
  L_optW_eq hc

/-- Below the threshold `γ_j = 0`, stated through `eq:assumption4` at the stack weights. -/
theorem gammaR_eq_zero_of_not_assumption4 (hc : ∀ i, 0 < c i) {j : Fin r}
    (h : ¬ Assumption4 (fun i => θ i j) c (wStackR θ c j)) : gammaR θ c j = 0 :=
  gammaR_eq_zero hc (not_lt.mp fun hlt => h ((assumption4_wStackR_iff hc j).mpr hlt))

/-- Every stack weight is nonzero when every strength at the weighting component `j` is
positive (F8). This is what supplies the witness of the strict monotonicity lemmas at
`w = w_{·j}`. -/
theorem wStackR_ne_zero (hc : ∀ i, 0 < c i) {j : Fin r} (hposj : ∀ i, 0 < θ i j) (i : Fin M) :
    wStackR θ c j i ≠ 0 := by
  rw [wStackR_apply]
  exact ne_of_gt (div_pos (hposj i)
    (Real.sqrt_pos.mpr (add_pos_of_nonneg_of_pos (sq_nonneg _) (hc i))))

/-! ### 8. The prefix property inside the model class -/

/-- Inside the model class a component to the left of a supercritical one is supercritical.
The hypotheses are those of `ellR_eq_val`: `hanti` orders the strengths the same way in every
table, so the energy profiles are ordered coordinatewise. -/
theorem assumption4_of_lt (hc : ∀ i, 0 < c i) (hnn : ∀ i k, 0 ≤ θ i k)
    (hanti : ∀ i, StrictAnti (θ i)) {j l : Fin r} (hlj : l < j)
    (h4 : Assumption4 (fun i => θ i j) c w) : Assumption4 (fun i => θ i l) c w := by
  refine MPhet.assumption4_mono hc (fun i => ?_) h4
  show θ i j ^ 2 ≤ θ i l ^ 2
  have h1 : θ i j < θ i l := hanti i hlj
  have h2 : 0 ≤ θ i j := hnn i j
  nlinarith

/-- Inside the model class the outliers are strictly ordered: a component to the left of a
supercritical one has a strictly larger outlier. `hw` holds at the stack weights by
`wStackR_ne_zero`. -/
theorem rhoHet_lt_of_lt (hc : ∀ i, 0 < c i) (hnn : ∀ i k, 0 ≤ θ i k)
    (hanti : ∀ i, StrictAnti (θ i)) (hM : 0 < M) (hw : ∀ i, w i ≠ 0) {j l : Fin r} (hlj : l < j)
    (h4 : Assumption4 (fun i => θ i j) c w) :
    MPhet.rhoHet (fun i => θ i j) c w < MPhet.rhoHet (fun i => θ i l) c w := by
  refine MPhet.rhoHet_lt_rhoHet hc (fun i => ?_) ⟨⟨0, hM⟩, hw _, ?_⟩ h4
  · show θ i j ^ 2 ≤ θ i l ^ 2
    have h1 : θ i j < θ i l := hanti i hlj
    have h2 : 0 ≤ θ i j := hnn i j
    nlinarith
  · show θ (⟨0, hM⟩ : Fin M) j ^ 2 < θ (⟨0, hM⟩ : Fin M) l ^ 2
    have h1 : θ (⟨0, hM⟩ : Fin M) j < θ (⟨0, hM⟩ : Fin M) l := hanti _ hlj
    have h2 : 0 ≤ θ (⟨0, hM⟩ : Fin M) j := hnn _ j
    nlinarith

/-! ### 9. The two index theorems -/

/-- **The outlier index is the paper's `ℓ_j` inside the model class.** At a supercritical
component `j` the filter of `ellSup` is exactly `Finset.Iio j`: a component to the left is
supercritical with a strictly larger outlier, `j` itself fails the strict inequality, and a
component to the right that is supercritical has a strictly smaller outlier. -/
theorem ellSup_eq_ellR (hc : ∀ i, 0 < c i) (hnn : ∀ i k, 0 ≤ θ i k)
    (hanti : ∀ i, StrictAnti (θ i)) (hM : 0 < M) {j : Fin r} (hposj : ∀ i, 0 < θ i j)
    (h4 : Assumption4 (fun i => θ i j) c (wStackR θ c j)) :
    ellSup θ c (wStackR θ c j) j = ellR θ c j := by
  classical
  have hw : ∀ i, wStackR θ c j i ≠ 0 := fun i => wStackR_ne_zero hc hposj i
  rw [ellSup, ellR_eq_val hc hnn hanti ⟨⟨0, hM⟩, hposj _⟩]
  refine Eq.trans (congrArg Finset.card ?_) (Fin.card_Iio j)
  ext l
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_Iio]
  constructor
  · rintro ⟨hA, hR⟩
    rcases lt_trichotomy l j with h | h | h
    · exact h
    · subst h
      exact absurd hR (lt_irrefl _)
    · exact absurd hR (not_lt.mpr (rhoHet_lt_of_lt hc hnn hanti hM hw h hA).le)
  · intro hlj
    exact ⟨assumption4_of_lt hc hnn hanti hlj h4,
      rhoHet_lt_of_lt hc hnn hanti hM hw hlj h4⟩

/-- **The subcritical branch.** When component `j` is not supercritical at its own weights,
every supercritical component sits strictly to the left of `j`, so the number of outliers is at
most `ℓ_j`. Stage E8 reads the bulk branch at every index from `numSup` up. -/
theorem numSup_le_of_not_assumption4 (hc : ∀ i, 0 < c i) (hnn : ∀ i k, 0 ≤ θ i k)
    (hanti : ∀ i, StrictAnti (θ i)) (hM : 0 < M) {j : Fin r} (hposj : ∀ i, 0 < θ i j)
    (h4 : ¬ Assumption4 (fun i => θ i j) c (wStackR θ c j)) :
    numSup θ c (wStackR θ c j) ≤ ellR θ c j := by
  classical
  rw [numSup, ellR_eq_val hc hnn hanti ⟨⟨0, hM⟩, hposj _⟩, ← Fin.card_Iio j]
  refine Finset.card_le_card fun l hl => ?_
  have hA := (Finset.mem_filter.mp hl).2
  rw [Finset.mem_Iio]
  rcases lt_trichotomy l j with h | h | h
  · exact h
  · subst h
    exact absurd hA h4
  · exact absurd (assumption4_of_lt hc hnn hanti h hA) h4

end Scalars

end StackedSVD
