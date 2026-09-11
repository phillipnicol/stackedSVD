/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVD
import StackedSVD.Scalars

/-!
# `thm:stacksvd_weighted`: optimally weighted stacked SVD

STATUS 2026-08-30: proved. `HeteroLaw` (black box 2) stays a hypothesis; everything else is
proved from it. The review note is `notes/archive/thm_stacksvd_weighted.md`, the adversarial audit
`notes/archive/audit_weighted_stacksvd_2026-08-30.md`.

## What the paper does (`main_paper.tex:1350`, `app:weightedStackSVDProof`)

The weighted stack is `X_stack(w) = [w_1 X_1; ...; w_M X_M] = ũ_0 vᵀ + Σ^{1/2} E` with
`ũ_0` the block vector `(θ_i w_i u_i)_i` and `Σ = diag(w_1²,...,w_1²,...,w_M²,...,w_M²)`,
block `i` repeated `n_i` times. The noise is therefore **block heteroscedastic**, so the
unweighted route of `StackSVD.lean` (the stack is again an i.i.d. Gaussian table) is not
available: `stackW w` has `GaussianNoise` only when every `|w_i| = 1`. The paper's input is
black box 2 of `notes/README.md` (BEJ1129 Theorem 2.3, liu2023asymptotic Theorem 2),
which is `HeteroLaw` here.

## Content

1. `Scalars.secular`, `Scalars.gammaTop`: the secular equation of `lem:secular_equation` and
   its top root `γ₁`, with existence and uniqueness (`existsUnique_gammaTop`), in the style
   of `Scalars.gW` and `Scalars.stackSVDLimitW`.
2. `Scalars.eta1`, `Scalars.Assumption4` (`eq:assumption4`) and `Scalars.Lw`, the paper's
   `L(w)`: the specialized limit of black box 2.
3. `Scalars.pOf`, the change of variables `p_j = θ_j²w_j²/(γ₁ - w_j²)` of
   `eq:var_change_wstacksvd`, with `sum_pOf` (`∑ p_j = 1`), `LwDen_eq`
   (`γ₁∑θ_j²w_j²/(γ₁-w_j²)² = 1 + ∑p_j²/θ_j²`) and `sum_c_pOf_le`, the zero-signal reduction
   of `main_paper.tex:1504`.
4. `Scalars.optWstack θ c i = θ_i/√(θ_i²+c_i)`, the paper's optimal weights, and the three
   scalar results: `L_le_opt` (the simplex optimization: no weighting beats `w⋆`),
   `L_optW_eq` (`L(w⋆)` is the root of `gW = 1`) and `assumption4_optW_iff` (the threshold
   equivalence of `main_paper.tex:1602`).
5. `MultiTableModel.stackW`, the weighted stack as a `SpikedModel`, and `stackW_X`, the
   stacking lemma. `stackW_one` says `stackW 1 = stack`.
6. `MultiTableModel.HeteroLaw`, the hypothesis structure, and the theorems
   `thm_stacksvd_weighted_general`, `thm_stacksvd_weighted`, `thm_stacksvd_weighted_inner`.
7. `heteroLaw_one_of_singleTableLaw`: at `w = 1` the new structure follows from the
   `SingleTableLaw` of the stack, which ties it to the proved `prop_stacksvd_general`.

`thm_stacksvd_weighted` and the `cor.2` Gaussian discharge moved to
`StackSVD/Weighted.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace Scalars

variable {M : ℕ}

/-! ### The secular equation (`lem:secular_equation`)

`R = ũ_0 ũ_0ᵀ + Σ` has the eigenvalue `w_i²` with multiplicity at least `n_i - 1`, and its
remaining eigenvalues are the roots of `f(λ) = 1 + ∑_j θ_j² w_j²/(w_j² - λ)`. Only the top
root `γ₁ > max_i w_i²` enters `L(w)`.
-/

/-- `max_i w_i²`, the largest block variance of `Σ`. The junk value at `M = 0` is `0`. -/
noncomputable def wSqMax (w : Fin M → ℝ) : ℝ := ⨆ i, w i ^ 2

theorem le_wSqMax (w : Fin M → ℝ) (i : Fin M) : w i ^ 2 ≤ wSqMax w :=
  le_ciSup (f := fun i => w i ^ 2) (Set.Finite.bddAbove (Set.finite_range _)) i

theorem wSqMax_nonneg (w : Fin M → ℝ) : 0 ≤ wSqMax w :=
  Real.iSup_nonneg fun _ => sq_nonneg _

/-- `max_i (t w_i)² = t² max_i w_i²`. Both sides are `0` at `t = 0` and at `M = 0`. -/
theorem wSqMax_smul (w : Fin M → ℝ) (t : ℝ) :
    wSqMax (fun i => t * w i) = t ^ 2 * wSqMax w := by
  rcases eq_or_ne t 0 with rfl | ht
  · simp [wSqMax]
  · have ht2 : (0:ℝ) < t ^ 2 := lt_of_le_of_ne (sq_nonneg t) (Ne.symm (pow_ne_zero 2 ht))
    have h1 : wSqMax (fun i => t * w i) ≤ t ^ 2 * wSqMax w := by
      refine Real.iSup_le (fun i => ?_) (mul_nonneg ht2.le (wSqMax_nonneg w))
      have h := le_wSqMax w i
      calc (t * w i) ^ 2 = t ^ 2 * w i ^ 2 := by ring
        _ ≤ t ^ 2 * wSqMax w := mul_le_mul_of_nonneg_left h ht2.le
    have h2 : wSqMax w ≤ wSqMax (fun i => t * w i) / t ^ 2 := by
      refine Real.iSup_le (fun i => ?_)
        (div_nonneg (wSqMax_nonneg _) ht2.le)
      rw [le_div_iff₀ ht2]
      have h := le_wSqMax (fun i => t * w i) i
      calc w i ^ 2 * t ^ 2 = (t * w i) ^ 2 := by ring
        _ ≤ _ := h
    have h3 : t ^ 2 * wSqMax w ≤ t ^ 2 * (wSqMax (fun i => t * w i) / t ^ 2) :=
      mul_le_mul_of_nonneg_left h2 ht2.le
    have h4 : t ^ 2 * (wSqMax (fun i => t * w i) / t ^ 2) = wSqMax (fun i => t * w i) := by
      field_simp
    linarith [h3, h4 ▸ h3]

/-- The secular function `f(λ) = 1 + ∑_j θ_j² w_j²/(w_j² - λ)` of `lem:secular_equation`. -/
noncomputable def secular (θ w : Fin M → ℝ) (lam : ℝ) : ℝ :=
  1 + ∑ j, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - lam)

/-- The summand of `f` on the positive denominator. -/
private theorem secular_term_eq (θ w : Fin M → ℝ) (j : Fin M) (x : ℝ) :
    θ j ^ 2 * w j ^ 2 / (w j ^ 2 - x) = -(θ j ^ 2 * w j ^ 2 / (x - w j ^ 2)) := by
  rw [show w j ^ 2 - x = -(x - w j ^ 2) by ring, div_neg]

/-- `g` is the top eigenvalue `γ₁` of `R`: a root of the secular equation above every block
variance. `eq:assumption4` states `γ₁ > max_i w_i²` as part of the detectability condition. -/
def IsGammaTop (θ w : Fin M → ℝ) (g : ℝ) : Prop :=
  wSqMax w < g ∧ secular θ w g = 0

/-- The condition under which the outlier root exists: a table with nonzero signal carries a
largest weight. It holds at `w = 1` when some `θ_i ≠ 0`, and at `w = optWstack θ c` when some
`θ_i ≠ 0`, since there `w_i = 0` exactly on the zero-signal tables. -/
def MaxWeightSignal (θ w : Fin M → ℝ) : Prop :=
  ∃ k, θ k * w k ≠ 0 ∧ w k ^ 2 = wSqMax w

/-- A root of the secular equation forces a table with nonzero signal and nonzero weight:
otherwise `f` is the constant `1`. -/
theorem exists_signal_of_root {θ w : Fin M → ℝ} {g : ℝ} (hg : IsGammaTop θ w g) :
    ∃ j, θ j * w j ≠ 0 := by
  by_contra hcon
  simp only [not_exists, not_not] at hcon
  have hz : ∀ j : Fin M, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - g) = 0 := by
    intro j
    rw [show θ j ^ 2 * w j ^ 2 = (θ j * w j) ^ 2 by ring, hcon j]
    simp
  have h1 : secular θ w g = 1 := by
    unfold secular
    rw [Finset.sum_congr rfl fun j _ => hz j, Finset.sum_const_zero, add_zero]
  rw [hg.2] at h1
  norm_num at h1

/-- Above `max_i w_i²` every summand of `f` increases, and the summand of a table with
`θ_j w_j ≠ 0` increases strictly, so `f` is strictly increasing there. This gives uniqueness
of `γ₁` (the paper cites the interlacing theorem of `bunch1978rank`). The hypothesis `hθw` is
needed: with `θ_j w_j = 0` for every `j` the function `f` is the constant `1`. -/
theorem secular_strictMonoOn {θ w : Fin M → ℝ} (hθw : ∃ j, θ j * w j ≠ 0) :
    StrictMonoOn (secular θ w) (Set.Ioi (wSqMax w)) := by
  obtain ⟨j0, hj0⟩ := hθw
  have hA0 : 0 < θ j0 ^ 2 * w j0 ^ 2 := by
    rw [show θ j0 ^ 2 * w j0 ^ 2 = (θ j0 * w j0) ^ 2 by ring]
    exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hj0))
  intro x hx y hy hxy
  simp only [Set.mem_Ioi] at hx hy
  have hterm : ∀ j : Fin M,
      θ j ^ 2 * w j ^ 2 / (w j ^ 2 - x) ≤ θ j ^ 2 * w j ^ 2 / (w j ^ 2 - y) := by
    intro j
    have hb := le_wSqMax w j
    have hbx : (0:ℝ) < x - w j ^ 2 := by linarith
    have h := div_le_div_of_nonneg_left (a := θ j ^ 2 * w j ^ 2) (by positivity) hbx
      (by linarith : x - w j ^ 2 ≤ y - w j ^ 2)
    rw [secular_term_eq θ w j x, secular_term_eq θ w j y]
    linarith
  have hstrict : θ j0 ^ 2 * w j0 ^ 2 / (w j0 ^ 2 - x)
      < θ j0 ^ 2 * w j0 ^ 2 / (w j0 ^ 2 - y) := by
    have hb := le_wSqMax w j0
    have hbx : (0:ℝ) < x - w j0 ^ 2 := by linarith
    have h := div_lt_div_of_pos_left hA0 hbx (by linarith : x - w j0 ^ 2 < y - w j0 ^ 2)
    rw [secular_term_eq θ w j0 x, secular_term_eq θ w j0 y]
    linarith
  have hsum := Finset.sum_lt_sum (fun j (_ : j ∈ Finset.univ) => hterm j)
    ⟨j0, Finset.mem_univ _, hstrict⟩
  unfold secular
  linarith

/-- `f` is continuous to the right of `max_i w_i²`. -/
theorem secular_continuousOn {θ w : Fin M → ℝ} {a b : ℝ} (ha : wSqMax w < a) :
    ContinuousOn (secular θ w) (Set.Icc a b) := by
  unfold secular
  refine continuousOn_const.add (continuousOn_finsetSum _ fun j _ => ContinuousOn.div
    continuousOn_const (continuousOn_const.sub continuousOn_id) fun x hx => ?_)
  have h := le_wSqMax w j
  have : w j ^ 2 < x := by
    have := hx.1
    linarith
  exact sub_ne_zero.mpr (ne_of_lt this)

/-- Existence and uniqueness of `γ₁`. `f → -∞` at the right of `max_i w_i²` (the term of the
index `k` of `MaxWeightSignal` blows up) and `f → 1` at `+∞`. -/
theorem existsUnique_gammaTop {θ w : Fin M → ℝ} (h : MaxWeightSignal θ w) :
    ∃! g : ℝ, IsGammaTop θ w g := by
  obtain ⟨k, hk, hkmax⟩ := h
  have hAk : 0 < θ k ^ 2 * w k ^ 2 := by
    rw [show θ k ^ 2 * w k ^ 2 = (θ k * w k) ^ 2 by ring]
    exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hk))
  have hS0 : 0 ≤ ∑ j, θ j ^ 2 * w j ^ 2 := Finset.sum_nonneg fun j _ => by positivity
  have hSk : θ k ^ 2 * w k ^ 2 ≤ ∑ j, θ j ^ 2 * w j ^ 2 :=
    Finset.single_le_sum (f := fun j => θ j ^ 2 * w j ^ 2)
      (fun j _ => by positivity) (Finset.mem_univ k)
  set m := wSqMax w with hm
  set S := ∑ j, θ j ^ 2 * w j ^ 2 with hSdef
  set a := m + θ k ^ 2 * w k ^ 2 with ha
  set b := m + S + 1 with hb
  have hma : m < a := by rw [ha]; linarith
  have hab : a ≤ b := by rw [ha, hb]; linarith
  -- `f(a) ≤ 0`: every term is nonpositive and the `k`-th is exactly `-1`.
  have hfa : secular θ w a ≤ 0 := by
    have hnp : ∀ j : Fin M, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - a) ≤ 0 := by
      intro j
      have hb2 := le_wSqMax w j
      rw [secular_term_eq θ w j a]
      exact neg_nonpos.mpr (div_nonneg (by positivity) (by rw [← hm] at hb2; linarith))
    have hden : w k ^ 2 - a = -(θ k ^ 2 * w k ^ 2) := by rw [ha, hkmax]; ring
    have hkval : θ k ^ 2 * w k ^ 2 / (w k ^ 2 - a) = -1 := by
      rw [hden, div_neg, div_self hAk.ne']
    have hsingle := Finset.single_le_sum
      (f := fun j => -(θ j ^ 2 * w j ^ 2 / (w j ^ 2 - a)))
      (fun j _ => neg_nonneg.mpr (hnp j)) (Finset.mem_univ k)
    rw [Finset.sum_neg_distrib] at hsingle
    unfold secular
    rw [hkval] at hsingle
    linarith
  -- `f(b) ≥ 0`: every term is at least `-θ_j²w_j²/(S+1)` and `S/(S+1) ≤ 1`.
  have hS1 : (0:ℝ) < S + 1 := by linarith
  have hfb : 0 ≤ secular θ w b := by
    have hterm : ∀ j : Fin M,
        0 ≤ θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1) := by
      intro j
      have hb2 := le_wSqMax w j
      rw [← hm] at hb2
      have hle : S + 1 ≤ b - w j ^ 2 := by rw [hb]; linarith
      have h := div_le_div_of_nonneg_left (a := θ j ^ 2 * w j ^ 2) (by positivity) hS1 hle
      rw [secular_term_eq θ w j b]
      linarith
    have hpos : 0 ≤ ∑ j, (θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1)) :=
      Finset.sum_nonneg fun j _ => hterm j
    rw [Finset.sum_add_distrib, ← Finset.sum_div, ← hSdef] at hpos
    have hSS : S / (S + 1) ≤ 1 := by rw [div_le_one hS1]; linarith
    unfold secular
    linarith
  obtain ⟨g, hgmem, hgval⟩ :=
    intermediate_value_Icc hab (secular_continuousOn (b := b) hma) ⟨hfa, hfb⟩
  have hgroot : IsGammaTop θ w g := ⟨lt_of_lt_of_le hma hgmem.1, hgval⟩
  refine ⟨g, hgroot, fun y hy => ?_⟩
  exact (secular_strictMonoOn ⟨k, hk⟩).injOn hy.1 hgroot.1 (by rw [hy.2, hgroot.2])

open Classical in
/-- `γ₁`, the top eigenvalue of `R = ũ_0 ũ_0ᵀ + Σ`: the unique root of the secular equation
above `max_i w_i²` when there is one, and `0` otherwise (then no outlier separates and the
paper's `L(w)` is `0`). Same `dite` shape as `stackSVDLimitW`. -/
noncomputable def gammaTop (θ w : Fin M → ℝ) : ℝ :=
  if h : ∃! g : ℝ, IsGammaTop θ w g then h.choose else 0

theorem gammaTop_spec {θ w : Fin M → ℝ} (h : ∃! g : ℝ, IsGammaTop θ w g) :
    IsGammaTop θ w (gammaTop θ w) := by
  unfold gammaTop
  split_ifs
  exact h.choose_spec.1

/-- Any root above `max_i w_i²` is `γ₁`. A root gives `∑_j θ_j²w_j²/(w_j² - g) = -1`, so some
`θ_j w_j ≠ 0`, and `secular_strictMonoOn` then gives uniqueness. -/
theorem gammaTop_eq {θ w : Fin M → ℝ} {g : ℝ} (hg : IsGammaTop θ w g) :
    gammaTop θ w = g := by
  have hu : ∃! g' : ℝ, IsGammaTop θ w g' := by
    refine ⟨g, hg, fun y hy => ?_⟩
    exact (secular_strictMonoOn (exists_signal_of_root hg)).injOn hy.1 hg.1 (by rw [hy.2, hg.2])
  unfold gammaTop
  split_ifs
  exact (hu.choose_spec.2 g hg).symm

/-- With no root the junk value is `0`. -/
theorem gammaTop_eq_zero {θ w : Fin M → ℝ} (h : ¬ ∃ g : ℝ, IsGammaTop θ w g) :
    gammaTop θ w = 0 := by
  unfold gammaTop
  exact dif_neg fun hu => h hu.exists

theorem gammaTop_of_exists {θ w : Fin M → ℝ} (h : ∃ g : ℝ, IsGammaTop θ w g) :
    IsGammaTop θ w (gammaTop θ w) := by
  obtain ⟨g, hg⟩ := h
  rw [gammaTop_eq hg]
  exact hg

/-- The secular equation at `γ₁`, in the form the proof of `thm:stacksvd_weighted` uses:
`∑_j p_j = 1` for `p_j = θ_j² w_j²/(γ₁ - w_j²)` (`eq:var_change_wstacksvd`). -/
theorem sum_secular_eq_neg_one {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    ∑ j, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - gammaTop θ w) = -1 := by
  obtain ⟨g, hg⟩ := h
  rw [gammaTop_eq hg]
  have h2 := hg.2
  unfold secular at h2
  linarith

/-- At unit weights `1 + ‖θ‖₂²` is a root. -/
theorem isGammaTop_one {θ : Fin M → ℝ} (hθ : ∃ i, θ i ≠ 0) :
    IsGammaTop θ 1 (1 + ∑ i, θ i ^ 2) := by
  obtain ⟨i0, hi0⟩ := hθ
  have hne : Nonempty (Fin M) := ⟨i0⟩
  have hT : 0 < ∑ i, θ i ^ 2 := by
    refine Finset.sum_pos' (fun i _ => sq_nonneg _) ⟨i0, Finset.mem_univ _, ?_⟩
    exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi0))
  constructor
  · have : wSqMax (1 : Fin M → ℝ) = 1 := by
      unfold wSqMax
      simp
    rw [this]
    linarith
  · unfold secular
    have hterm : ∀ j : Fin M,
        θ j ^ 2 * (1 : Fin M → ℝ) j ^ 2 / ((1 : Fin M → ℝ) j ^ 2 - (1 + ∑ i, θ i ^ 2))
          = -(θ j ^ 2 / ∑ i, θ i ^ 2) := by
      intro j
      have : (1 : Fin M → ℝ) j = 1 := rfl
      rw [this]
      rw [show (1:ℝ) ^ 2 - (1 + ∑ i, θ i ^ 2) = -(∑ i, θ i ^ 2) by ring]
      rw [div_neg]
      ring_nf
    rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_neg_distrib, ← Finset.sum_div,
      div_self hT.ne']
    ring

/-- At unit weights `γ₁ = 1 + ‖θ‖₂²`, the value the paper uses for `prop:stacksvd_general`. -/
theorem gammaTop_one {θ : Fin M → ℝ} (hθ : ∃ i, θ i ≠ 0) :
    gammaTop θ 1 = 1 + ∑ i, θ i ^ 2 :=
  gammaTop_eq (isGammaTop_one hθ)

/-- `f` is unchanged by `(w, λ) ↦ (t w, t² λ)`. -/
theorem secular_smul {θ w : Fin M → ℝ} {t : ℝ} (ht : t ≠ 0) (lam : ℝ) :
    secular θ (fun i => t * w i) (t ^ 2 * lam) = secular θ w lam := by
  unfold secular
  congr 1
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [show θ j ^ 2 * (t * w j) ^ 2 = t ^ 2 * (θ j ^ 2 * w j ^ 2) by ring,
    show (t * w j) ^ 2 - t ^ 2 * lam = t ^ 2 * (w j ^ 2 - lam) by ring]
  exact mul_div_mul_left _ _ (pow_ne_zero 2 ht)

theorem isGammaTop_smul {θ w : Fin M → ℝ} {t g : ℝ} (ht : t ≠ 0) :
    IsGammaTop θ (fun i => t * w i) (t ^ 2 * g) ↔ IsGammaTop θ w g := by
  have ht2 : (0:ℝ) < t ^ 2 := lt_of_le_of_ne (sq_nonneg t) (Ne.symm (pow_ne_zero 2 ht))
  unfold IsGammaTop
  rw [wSqMax_smul, secular_smul ht]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨lt_of_mul_lt_mul_left h1 ht2.le, h2⟩
  · rintro ⟨h1, h2⟩
    exact ⟨mul_lt_mul_of_pos_left h1 ht2, h2⟩

/-- `R` scales as `t²R` under `w ↦ t w`, so `γ₁` scales as `t²γ₁`. At `t = 0` both sides are
`0`: the weighted stack has no signal, so no root exists (audit item 4). -/
theorem gammaTop_smul {θ w : Fin M → ℝ} {t : ℝ} :
    gammaTop θ (fun i => t * w i) = t ^ 2 * gammaTop θ w := by
  rcases eq_or_ne t 0 with rfl | ht
  · have hno : ¬ ∃ g : ℝ, IsGammaTop θ (fun i => (0:ℝ) * w i) g := by
      rintro ⟨g, -, hs⟩
      unfold secular at hs
      rw [Finset.sum_congr rfl (fun j (_ : j ∈ Finset.univ) =>
        show θ j ^ 2 * ((0:ℝ) * w j) ^ 2 / (((0:ℝ) * w j) ^ 2 - g) = 0 by simp),
        Finset.sum_const_zero, add_zero] at hs
      norm_num at hs
    rw [gammaTop_eq_zero hno]
    ring
  · by_cases hr : ∃ g : ℝ, IsGammaTop θ w g
    · obtain ⟨g, hg⟩ := hr
      rw [gammaTop_eq hg, gammaTop_eq ((isGammaTop_smul ht).mpr hg)]
    · have hrt : ¬ ∃ g' : ℝ, IsGammaTop θ (fun i => t * w i) g' := by
        rintro ⟨g', hg'⟩
        refine hr ⟨g' / t ^ 2, (isGammaTop_smul ht).mp ?_⟩
        have he : t ^ 2 * (g' / t ^ 2) = g' := by
          field_simp
        rw [he]
        exact hg'
      rw [gammaTop_eq_zero hrt, gammaTop_eq_zero hr]
      ring

/-! ### `eq:assumption4` and the limit `L(w)` (black box 2 of `notes/README.md`) -/

/-- `η₁ = 1 - ∑_i c_i w_i⁴/(γ₁ - w_i²)²`, the numerator of `L(w)`. -/
noncomputable def eta1 (θ c w : Fin M → ℝ) : ℝ :=
  1 - ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2

/-- `γ₁ ∑_j θ_j² w_j²/(γ₁ - w_j²)²`, the denominator of `L(w)`. -/
noncomputable def LwDen (θ w : Fin M → ℝ) : ℝ :=
  gammaTop θ w * ∑ j, θ j ^ 2 * w j ^ 2 / (gammaTop θ w - w j ^ 2) ^ 2

/-- `eq:assumption4`: the outlier root exists (`γ₁ > max_i w_i²` is part of `IsGammaTop`) and
`∑_i c_i w_i⁴/(γ₁ - w_i²)² < 1`. When it fails, every singular value of the weighted stack
lies in the bulk and the overlap limit is `0`. -/
def Assumption4 (θ c w : Fin M → ℝ) : Prop :=
  (∃ g, IsGammaTop θ w g) ∧ ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 < 1

theorem eta1_pos_iff (θ c w : Fin M → ℝ) :
    0 < eta1 θ c w ↔ ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 < 1 :=
  sub_pos

open Classical in
/-- The paper's `L(w)`: the limit of `(v̂_stack(w)ᵀ v)²` for a general weighting.
`L(w) = η₁/(γ₁ ∑_j θ_j² w_j²/(γ₁ - w_j²)²)` under `eq:assumption4`, and `0` otherwise. -/
noncomputable def Lw (θ c w : Fin M → ℝ) : ℝ :=
  if Assumption4 θ c w then eta1 θ c w / LwDen θ w else 0

/-! ### The change of variables `eq:var_change_wstacksvd` -/

/-- `p_j = θ_j² w_j²/(γ₁ - w_j²)`, the simplex coordinate of `app:weightedStackSVDProof`. -/
noncomputable def pOf (θ w : Fin M → ℝ) (j : Fin M) : ℝ :=
  θ j ^ 2 * w j ^ 2 / (gammaTop θ w - w j ^ 2)

theorem gammaTop_sub_pos {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) (j : Fin M) :
    0 < gammaTop θ w - w j ^ 2 := by
  have h1 := (gammaTop_of_exists h).1
  have h2 := le_wSqMax w j
  linarith

theorem pOf_nonneg {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) (j : Fin M) :
    0 ≤ pOf θ w j :=
  div_nonneg (by positivity) (gammaTop_sub_pos h j).le

theorem pOf_eq_zero {θ w : Fin M → ℝ} {j : Fin M} (hz : θ j = 0) : pOf θ w j = 0 := by
  unfold pOf
  rw [hz]
  simp

/-- `∑_j p_j = 1`, the simplex constraint of `main_paper.tex:1517`. -/
theorem sum_pOf {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) : ∑ j, pOf θ w j = 1 := by
  have hs := sum_secular_eq_neg_one h
  have hterm : ∀ j : Fin M,
      pOf θ w j = -(θ j ^ 2 * w j ^ 2 / (w j ^ 2 - gammaTop θ w)) := by
    intro j
    rw [secular_term_eq θ w j (gammaTop θ w), neg_neg]
    rfl
  rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_neg_distrib, hs]
  ring

/-- `γ₁ ∑_j θ_j²w_j²/(γ₁-w_j²)² = 1 + ∑_j p_j²/θ_j²` (`main_paper.tex:1536`). The junk value
`0/0 = 0` covers the zero-signal tables, where `p_j = 0`. -/
theorem LwDen_eq {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    LwDen θ w = 1 + ∑ j, pOf θ w j ^ 2 / θ j ^ 2 := by
  have hterm : ∀ j : Fin M,
      gammaTop θ w * (θ j ^ 2 * w j ^ 2 / (gammaTop θ w - w j ^ 2) ^ 2)
        = pOf θ w j + pOf θ w j ^ 2 / θ j ^ 2 := by
    intro j
    have hd := gammaTop_sub_pos h j
    rcases eq_or_ne (θ j) 0 with hz | hz
    · rw [pOf_eq_zero hz, hz]
      simp
    · unfold pOf
      field_simp
      ring
  unfold LwDen
  rw [Finset.mul_sum, Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_add_distrib,
    sum_pOf h]

/-- The zero-signal reduction of `main_paper.tex:1504`. The tables with `θ_j = 0` do not
enter the secular equation, so they leave `p` and the denominator of `L(w)` unchanged and
only add the nonnegative term `c_j w_j⁴/(γ₁-w_j²)²` to `∑_i c_i w_i⁴/(γ₁-w_i²)²`. Setting
`w_j = 0` there (which `optWstack` does) therefore only increases `L(w)`. -/
theorem sum_c_pOf_le {θ c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) (h : ∃ g, IsGammaTop θ w g) :
    ∑ j, c j * pOf θ w j ^ 2 / θ j ^ 4
      ≤ ∑ j, c j * w j ^ 4 / (gammaTop θ w - w j ^ 2) ^ 2 := by
  refine Finset.sum_le_sum fun j _ => ?_
  have hd := gammaTop_sub_pos h j
  rcases eq_or_ne (θ j) 0 with hz | hz
  · rw [pOf_eq_zero hz, hz]
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, div_zero]
    exact div_nonneg (mul_nonneg (hc j) (by positivity)) (by positivity)
  · refine le_of_eq ?_
    unfold pOf
    field_simp

/-! ### Elementary facts about `L(w)` -/

theorem Lw_nonneg (θ c w : Fin M → ℝ) : 0 ≤ Lw θ c w := by
  unfold Lw
  split_ifs with h4
  · refine div_nonneg ((eta1_pos_iff θ c w).mpr h4.2).le ?_
    unfold LwDen
    refine mul_nonneg ?_ (Finset.sum_nonneg fun j _ => div_nonneg (by positivity) (by positivity))
    exact le_of_lt (lt_of_le_of_lt (wSqMax_nonneg w) (gammaTop_of_exists h4.1).1)
  · exact le_rfl

theorem one_le_LwDen {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) : 1 ≤ LwDen θ w := by
  rw [LwDen_eq h]
  have : 0 ≤ ∑ j, pOf θ w j ^ 2 / θ j ^ 2 :=
    Finset.sum_nonneg fun j _ => by positivity
  linarith

theorem Lw_le_one {θ c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) : Lw θ c w ≤ 1 := by
  unfold Lw
  split_ifs with h4
  · have hD := one_le_LwDen h4.1
    rw [div_le_one (by linarith)]
    have hE : eta1 θ c w ≤ 1 := by
      unfold eta1
      have : 0 ≤ ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 :=
        Finset.sum_nonneg fun i _ =>
          div_nonneg (mul_nonneg (hc i) (by positivity)) (by positivity)
      linarith
    linarith
  · norm_num

/-- `L(w)` is scale invariant: `R` scales as `t²R` and every ratio in `L` is unchanged. This
is what makes the paper's `w_i⋆ ∝ θ_i/√(θ_i²+c_i)` a statement about one weight vector. -/
theorem Lw_smul {θ c w : Fin M → ℝ} {t : ℝ} (ht : t ≠ 0) :
    Lw θ c (fun i => t * w i) = Lw θ c w := by
  have ht4 : (t:ℝ) ^ 4 ≠ 0 := pow_ne_zero 4 ht
  have hgs : gammaTop θ (fun i => t * w i) = t ^ 2 * gammaTop θ w := gammaTop_smul
  by_cases hr : ∃ g, IsGammaTop θ w g
  · obtain ⟨g, hg⟩ := hr
    have hrw : ∃ g' : ℝ, IsGammaTop θ w g' := ⟨g, hg⟩
    have hrt : ∃ g' : ℝ, IsGammaTop θ (fun i => t * w i) g' :=
      ⟨t ^ 2 * g, (isGammaTop_smul ht).mpr hg⟩
    have hdpos := gammaTop_sub_pos hrw
    have hden : ∀ i : Fin M,
        (gammaTop θ (fun i => t * w i) - (t * w i) ^ 2) ^ 2
          = t ^ 4 * (gammaTop θ w - w i ^ 2) ^ 2 := by
      intro i
      rw [hgs]
      ring
    have hAs : ∑ i, c i * (t * w i) ^ 4 / (gammaTop θ (fun i => t * w i) - (t * w i) ^ 2) ^ 2
        = ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 := by
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [hden i, show c i * (t * w i) ^ 4 = t ^ 4 * (c i * w i ^ 4) by ring]
      exact mul_div_mul_left _ _ ht4
    have ht2 : (t:ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 ht
    have hD : LwDen θ (fun i => t * w i) = LwDen θ w := by
      change gammaTop θ (fun i => t * w i)
            * ∑ i, θ i ^ 2 * (t * w i) ^ 2
              / (gammaTop θ (fun i => t * w i) - (t * w i) ^ 2) ^ 2
          = gammaTop θ w * ∑ i, θ i ^ 2 * w i ^ 2 / (gammaTop θ w - w i ^ 2) ^ 2
      rw [Finset.mul_sum, Finset.mul_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      have h0 : (gammaTop θ w - w i ^ 2) ≠ 0 := (hdpos i).ne'
      rw [hgs, show t ^ 2 * gammaTop θ w - (t * w i) ^ 2
        = t ^ 2 * (gammaTop θ w - w i ^ 2) by ring]
      field_simp
    have hEta : eta1 θ c (fun i => t * w i) = eta1 θ c w := by
      change (1:ℝ) - ∑ i, c i * (t * w i) ^ 4
            / (gammaTop θ (fun i => t * w i) - (t * w i) ^ 2) ^ 2
          = 1 - ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2
      rw [hAs]
    have hA4' : (∑ i, c i * (t * w i) ^ 4
          / (gammaTop θ (fun i => t * w i) - (t * w i) ^ 2) ^ 2 < 1)
        ↔ (∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 < 1) := by rw [hAs]
    have h4iff : Assumption4 θ c (fun i => t * w i) ↔ Assumption4 θ c w := by
      constructor
      · rintro ⟨-, h⟩
        exact ⟨hrw, hA4'.mp h⟩
      · rintro ⟨-, h⟩
        exact ⟨hrt, hA4'.mpr h⟩
    unfold Lw
    rw [hEta, hD]
    by_cases h4 : Assumption4 θ c w
    · rw [if_pos (h4iff.mpr h4), if_pos h4]
    · rw [if_neg fun hh => h4 (h4iff.mp hh), if_neg h4]
  · have hrt : ¬ ∃ g' : ℝ, IsGammaTop θ (fun i => t * w i) g' := by
      rintro ⟨g', hg'⟩
      refine hr ⟨g' / t ^ 2, (isGammaTop_smul ht).mp ?_⟩
      have he : t ^ 2 * (g' / t ^ 2) = g' := by field_simp
      rw [he]
      exact hg'
    unfold Lw
    rw [if_neg fun h => hrt h.1, if_neg fun h => hr h.1]

/-! ### The optimal weights and the three scalar results -/

/-- The paper's optimal weights `w_i⋆ ∝ θ_i/√(θ_i² + c_i)` (`thm:stacksvd_weighted`), at the
scale where the numerator is `θ_i`. A zero-signal table gets weight `0`, which is the paper's
first reduction (`w_j = 0` for `j ∉ J_{>0}`). -/
noncomputable def optWstack (θ c : Fin M → ℝ) : Fin M → ℝ :=
  fun i => θ i / Real.sqrt (θ i ^ 2 + c i)

/-- `L(1_M)` is the unweighted limit of `prop:stacksvd_general`. With `γ₁ = 1 + ‖θ‖₂²` the
formula collapses to `(‖θ‖₂⁴ - ‖c‖₁)/(‖θ‖₂²(‖θ‖₂² + 1))` (`main_paper.tex:1438`). The
degenerate case `θ = 0` gives `0` on both sides, so no hypothesis on `c` is needed. -/
theorem Lw_one (θ c : Fin M → ℝ) : Lw θ c 1 = stackSVDLimit θ c := by
  by_cases hθ : ∃ i, θ i ≠ 0
  · obtain ⟨i0, hi0⟩ := hθ
    have hθ' : ∃ i, θ i ≠ 0 := ⟨i0, hi0⟩
    have hT : 0 < ∑ i, θ i ^ 2 := by
      refine Finset.sum_pos' (fun i _ => sq_nonneg _) ⟨i0, Finset.mem_univ _, ?_⟩
      exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi0))
    have hg : gammaTop θ 1 = 1 + ∑ i, θ i ^ 2 := gammaTop_one hθ'
    have hroot : ∃ g : ℝ, IsGammaTop θ (1 : Fin M → ℝ) g := ⟨_, isGammaTop_one hθ'⟩
    have hone : ∀ i : Fin M, (1 : Fin M → ℝ) i = 1 := fun _ => rfl
    have hden : ∀ i : Fin M,
        gammaTop θ (1 : Fin M → ℝ) - (1 : Fin M → ℝ) i ^ 2 = ∑ l, θ l ^ 2 := by
      intro i
      rw [hg, hone i]
      ring
    have hA : ∑ i, c i * (1 : Fin M → ℝ) i ^ 4
          / (gammaTop θ (1 : Fin M → ℝ) - (1 : Fin M → ℝ) i ^ 2) ^ 2
        = (∑ i, c i) / (∑ l, θ l ^ 2) ^ 2 := by
      rw [Finset.sum_div]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [hden i, hone i]
      ring
    have hB : LwDen θ (1 : Fin M → ℝ) = (1 + ∑ l, θ l ^ 2) / ∑ l, θ l ^ 2 := by
      unfold LwDen
      have : ∀ i : Fin M, θ i ^ 2 * (1 : Fin M → ℝ) i ^ 2
          / (gammaTop θ (1 : Fin M → ℝ) - (1 : Fin M → ℝ) i ^ 2) ^ 2
          = θ i ^ 2 / (∑ l, θ l ^ 2) ^ 2 := by
        intro i
        rw [hden i, hone i]
        ring
      rw [Finset.sum_congr rfl fun i _ => this i, ← Finset.sum_div, hg]
      field_simp
    have hiff : Assumption4 θ c (1 : Fin M → ℝ) ↔ (∑ i, θ i ^ 2) ^ 2 > ∑ i, c i := by
      unfold Assumption4
      rw [hA]
      constructor
      · rintro ⟨-, h⟩
        rw [div_lt_one (by positivity)] at h
        exact h
      · intro h
        exact ⟨hroot, by rw [div_lt_one (by positivity)]; exact h⟩
    unfold Lw eta1 stackSVDLimit
    rw [hA, hB]
    by_cases hcase : (∑ i, θ i ^ 2) ^ 2 > ∑ i, c i
    · rw [if_pos (hiff.mpr hcase), if_pos hcase]
      field_simp
      ring
    · rw [if_neg fun hh => hcase (hiff.mp hh), if_neg hcase]
  · simp only [not_exists, not_not] at hθ
    have hT : ∑ i, θ i ^ 2 = 0 := Finset.sum_eq_zero fun i _ => by rw [hθ i]; ring
    have hno : ¬ ∃ g : ℝ, IsGammaTop θ (1 : Fin M → ℝ) g := by
      rintro ⟨g, -, hs⟩
      unfold secular at hs
      rw [Finset.sum_congr rfl (fun j (_ : j ∈ Finset.univ) =>
        show θ j ^ 2 * (1 : Fin M → ℝ) j ^ 2 / ((1 : Fin M → ℝ) j ^ 2 - g) = 0 by
          rw [hθ j]; simp), Finset.sum_const_zero, add_zero] at hs
      norm_num at hs
    unfold Lw stackSVDLimit
    rw [if_neg fun h => hno h.1, hT]
    simp

/-- `cor.2` at binary weights: `L(1_S)` is `binaryStackSVDLimit S θ c`, the unweighted
stacksvd limit of the sub-collection `S` (`Scalars.lean`). The scalar half of `cor.2`; the
model half is `heteroLaw_binary_of_singleTableLaw`, and, for Gaussian noise,
`heteroLaw_binary_of_gaussian`, both later in this file. -/
theorem Lw_binary (S : Finset (Fin M)) (θ c : Fin M → ℝ) :
    Lw θ c (fun i => if i ∈ S then 1 else 0) = binaryStackSVDLimit S θ c := by
  classical
  set w : Fin M → ℝ := fun i => if i ∈ S then 1 else 0 with hwdef
  have hw2 : ∀ i : Fin M, w i ^ 2 = if i ∈ S then 1 else 0 := by
    intro i
    rw [hwdef]
    by_cases hi : i ∈ S <;> simp [hi]
  have hw4 : ∀ i : Fin M, w i ^ 4 = if i ∈ S then 1 else 0 := by
    intro i
    rw [hwdef]
    by_cases hi : i ∈ S <;> simp [hi]
  by_cases hθ : ∃ i ∈ S, θ i ≠ 0
  · obtain ⟨i0, hi0S, hi0⟩ := hθ
    have hT : 0 < ∑ i ∈ S, θ i ^ 2 :=
      Finset.sum_pos' (fun i _ => sq_nonneg _) ⟨i0, hi0S,
        lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi0))⟩
    have hwmax : wSqMax w = 1 := by
      refine le_antisymm (Real.iSup_le (fun i => ?_) zero_le_one) ?_
      · rw [hw2 i]
        by_cases hi : i ∈ S <;> simp [hi]
      · have h := le_wSqMax w i0
        rw [hw2 i0, if_pos hi0S] at h
        exact h
    have hroot : IsGammaTop θ w (1 + ∑ i ∈ S, θ i ^ 2) := by
      constructor
      · rw [hwmax]
        linarith
      · change (1:ℝ) + ∑ j, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - (1 + ∑ i ∈ S, θ i ^ 2)) = 0
        have hterm : ∀ j : Fin M,
            θ j ^ 2 * w j ^ 2 / (w j ^ 2 - (1 + ∑ i ∈ S, θ i ^ 2))
              = if j ∈ S then -(θ j ^ 2 / ∑ i ∈ S, θ i ^ 2) else 0 := by
          intro j
          rw [hw2 j]
          by_cases hj : j ∈ S
          · rw [if_pos hj, if_pos hj,
              show (1:ℝ) - (1 + ∑ i ∈ S, θ i ^ 2) = -(∑ i ∈ S, θ i ^ 2) by ring, div_neg]
            ring_nf
          · rw [if_neg hj, if_neg hj]
            simp
        rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_ite_mem, Finset.univ_inter,
          Finset.sum_neg_distrib, ← Finset.sum_div, div_self hT.ne']
        ring
    have hgt : gammaTop θ w = 1 + ∑ i ∈ S, θ i ^ 2 := gammaTop_eq hroot
    have hex : ∃ g : ℝ, IsGammaTop θ w g := ⟨_, hroot⟩
    have hA : ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2
        = (∑ i ∈ S, c i) / (∑ i ∈ S, θ i ^ 2) ^ 2 := by
      have hterm : ∀ i : Fin M, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2
          = if i ∈ S then c i / (∑ l ∈ S, θ l ^ 2) ^ 2 else 0 := by
        intro i
        rw [hgt, hw2 i, hw4 i]
        by_cases hi : i ∈ S
        · rw [if_pos hi, if_pos hi,
            show (1:ℝ) + (∑ l ∈ S, θ l ^ 2) - 1 = ∑ l ∈ S, θ l ^ 2 by ring]
          ring
        · rw [if_neg hi, if_neg hi]
          simp
      rw [Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_ite_mem, Finset.univ_inter,
        ← Finset.sum_div]
    have hB : LwDen θ w = (1 + ∑ i ∈ S, θ i ^ 2) / ∑ i ∈ S, θ i ^ 2 := by
      change gammaTop θ w * ∑ j, θ j ^ 2 * w j ^ 2 / (gammaTop θ w - w j ^ 2) ^ 2 = _
      have hterm : ∀ j : Fin M, θ j ^ 2 * w j ^ 2 / (gammaTop θ w - w j ^ 2) ^ 2
          = if j ∈ S then θ j ^ 2 / (∑ l ∈ S, θ l ^ 2) ^ 2 else 0 := by
        intro j
        rw [hgt, hw2 j]
        by_cases hj : j ∈ S
        · rw [if_pos hj, if_pos hj,
            show (1:ℝ) + (∑ l ∈ S, θ l ^ 2) - 1 = ∑ l ∈ S, θ l ^ 2 by ring]
          ring
        · rw [if_neg hj, if_neg hj]
          simp
      have hTne : (∑ i ∈ S, θ i ^ 2) ≠ 0 := hT.ne'
      rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_ite_mem, Finset.univ_inter,
        ← Finset.sum_div, hgt]
      field_simp
    have hiff : Assumption4 θ c w ↔ (∑ i ∈ S, θ i ^ 2) ^ 2 > ∑ i ∈ S, c i := by
      unfold Assumption4
      rw [hA]
      constructor
      · rintro ⟨-, h⟩
        rw [div_lt_one (by positivity)] at h
        exact h
      · intro h
        exact ⟨hex, by rw [div_lt_one (by positivity)]; exact h⟩
    unfold Lw eta1 binaryStackSVDLimit
    rw [hA, hB]
    by_cases hcase : (∑ i ∈ S, θ i ^ 2) ^ 2 > ∑ i ∈ S, c i
    · rw [if_pos (hiff.mpr hcase), if_pos hcase]
      field_simp
      ring
    · rw [if_neg fun hh => hcase (hiff.mp hh), if_neg hcase]
  · simp only [not_exists, not_and, not_not] at hθ
    have hT : ∑ i ∈ S, θ i ^ 2 = 0 :=
      Finset.sum_eq_zero fun i hi => by rw [hθ i hi]; ring
    have hno : ¬ ∃ g : ℝ, IsGammaTop θ w g := by
      rintro ⟨g, -, hs⟩
      have hs' : (1:ℝ) + ∑ j, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - g) = 0 := hs
      have hterm : ∀ j : Fin M, θ j ^ 2 * w j ^ 2 / (w j ^ 2 - g) = 0 := by
        intro j
        rw [hw2 j]
        by_cases hj : j ∈ S
        · rw [if_pos hj, hθ j hj]
          simp
        · rw [if_neg hj]
          simp
      rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_const_zero, add_zero] at hs'
      norm_num at hs'
    unfold Lw binaryStackSVDLimit
    rw [if_neg fun h => hno h.1, hT]
    simp

/-- The simplex step of `app:weightedStackSVDProof`. Under `eq:assumption4` the value `L(w)`
is positive, at most `1`, and satisfies `gW θ c L(w) ≥ 1`: with `p_j` on the simplex,
`L(w) ≥ r` reads `∑_j (α_j + rβ_j)p_j² ≤ 1 - r`, and Cauchy-Schwarz bounds the left side
below by `1/∑_j 1/(α_j + rβ_j)`. -/
theorem one_le_gW_Lw {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    0 < Lw θ c w ∧ Lw θ c w ≤ 1 ∧ 1 ≤ gW θ c (Lw θ c w) := by
  classical
  obtain ⟨hroot, hA4⟩ := h4
  have h4' : Assumption4 θ c w := ⟨hroot, hA4⟩
  set L := Lw θ c w with hLdef
  have hDen := LwDen_eq hroot
  have hDpos : 0 < LwDen θ w := lt_of_lt_of_le one_pos (one_le_LwDen hroot)
  have hLeq : L = eta1 θ c w / LwDen θ w := by
    rw [hLdef]
    unfold Lw
    rw [if_pos h4']
  have hepos : 0 < eta1 θ c w := (eta1_pos_iff θ c w).mpr hA4
  have hLpos : 0 < L := by rw [hLeq]; exact div_pos hepos hDpos
  have hLle1 : L ≤ 1 := Lw_le_one fun i => (hc i).le
  refine ⟨hLpos, hLle1, ?_⟩
  -- the index set of the tables with signal
  set J : Finset (Fin M) := Finset.univ.filter (fun j => θ j ≠ 0) with hJdef
  have hmemJ : ∀ j, j ∈ J ↔ θ j ≠ 0 := by
    intro j
    rw [hJdef, Finset.mem_filter]
    simp
  have hPJ : ∑ j ∈ J, pOf θ w j = 1 := by
    rw [Finset.sum_subset (Finset.subset_univ J) ?_]
    · exact sum_pOf hroot
    · intro j _ hj
      exact pOf_eq_zero (not_not.mp fun hz => hj ((hmemJ j).mpr hz))
  have hJne : J.Nonempty := by
    by_contra hno
    rw [Finset.not_nonempty_iff_eq_empty] at hno
    rw [hno, Finset.sum_empty] at hPJ
    norm_num at hPJ
  have hSA : ∑ j ∈ J, c j * pOf θ w j ^ 2 / θ j ^ 4
      = ∑ j, c j * pOf θ w j ^ 2 / θ j ^ 4 := by
    refine Finset.sum_subset (Finset.subset_univ J) fun j _ hj => ?_
    have hz : θ j = 0 := not_not.mp fun hz => hj ((hmemJ j).mpr hz)
    rw [pOf_eq_zero hz, hz]
    simp
  have hSB : ∑ j ∈ J, pOf θ w j ^ 2 / θ j ^ 2 = ∑ j, pOf θ w j ^ 2 / θ j ^ 2 := by
    refine Finset.sum_subset (Finset.subset_univ J) fun j _ hj => ?_
    have hz : θ j = 0 := not_not.mp fun hz => hj ((hmemJ j).mpr hz)
    rw [pOf_eq_zero hz, hz]
    simp
  -- `∑_j (α_j + Lβ_j) p_j² ≤ 1 - L`
  have hkey : (∑ j ∈ J, c j * pOf θ w j ^ 2 / θ j ^ 4)
      + L * ∑ j ∈ J, pOf θ w j ^ 2 / θ j ^ 2 ≤ 1 - L := by
    have h1 : L * LwDen θ w = eta1 θ c w := by
      rw [hLeq]
      field_simp
    have h2 : eta1 θ c w ≤ 1 - ∑ j, c j * pOf θ w j ^ 2 / θ j ^ 4 := by
      unfold eta1
      linarith [sum_c_pOf_le (fun i => (hc i).le) hroot]
    rw [hDen] at h1
    rw [hSA, hSB]
    nlinarith [h1, h2]
  -- Cauchy-Schwarz on the simplex
  have hGpos : ∀ j ∈ J, 0 < (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹ := by
    intro j hj
    have hz : θ j ≠ 0 := (hmemJ j).mp hj
    have h4' : (0:ℝ) < θ j ^ 4 := by positivity
    have h2' : (0:ℝ) < θ j ^ 2 := by positivity
    exact inv_pos.mpr (add_pos (div_pos (hc j) h4') (div_pos hLpos h2'))
  have hTpos : 0 < ∑ j ∈ J, (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹ := Finset.sum_pos hGpos hJne
  have hcL : ∀ j : Fin M, c j + L * θ j ^ 2 ≠ 0 := fun j =>
    (add_pos_of_pos_of_nonneg (hc j) (mul_nonneg hLpos.le (sq_nonneg _))).ne'
  have hinv : ∀ j ∈ J, (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹ = θ j ^ 4 / (c j + L * θ j ^ 2) := by
    intro j hj
    have hz : θ j ≠ 0 := (hmemJ j).mp hj
    rw [show c j / θ j ^ 4 + L / θ j ^ 2 = (c j + L * θ j ^ 2) / θ j ^ 4 by field_simp, inv_div]
  have hCS := Finset.sq_sum_div_le_sum_sq_div J (pOf θ w)
    (g := fun j => (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹) hGpos
  have hrhs : ∑ j ∈ J, pOf θ w j ^ 2 / (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹
      = (∑ j ∈ J, c j * pOf θ w j ^ 2 / θ j ^ 4)
        + L * ∑ j ∈ J, pOf θ w j ^ 2 / θ j ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun j hj => ?_
    have hz : θ j ≠ 0 := (hmemJ j).mp hj
    have hL2 := hcL j
    rw [hinv j hj]
    field_simp
  rw [hPJ, hrhs] at hCS
  have hstep : 1 / (∑ j ∈ J, (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹) ≤ 1 - L := by
    have : (1:ℝ) ^ 2 = 1 := one_pow 2
    linarith [hCS, hkey]
  rw [div_le_iff₀ hTpos] at hstep
  have hgWJ : gW θ c L = ∑ j ∈ J, (1 - L) * (c j / θ j ^ 4 + L / θ j ^ 2)⁻¹ := by
    unfold gW
    rw [← Finset.sum_subset (Finset.subset_univ J) ?_]
    · refine Finset.sum_congr rfl fun j hj => ?_
      rw [hinv j hj]
      ring
    · intro j _ hj
      have hz : θ j = 0 := not_not.mp fun hz => hj ((hmemJ j).mpr hz)
      rw [hz]
      simp
  rw [hgWJ, ← Finset.mul_sum]
  linarith [hstep]

/-- The simplex optimization of `app:weightedStackSVDProof`: no weighting beats `w⋆`.
With `p_j = θ_j² w_j²/(γ₁ - w_j²)` on the simplex, `L(w) = (1 - ∑ α_j p_j²)/(1 + ∑ β_j p_j²)`
with `α_j = c_j/θ_j⁴` and `β_j = 1/θ_j²`, and `L(w) ≥ r` iff `∑(α_j + rβ_j)p_j² ≤ 1 - r`,
whose minimum over the simplex is `1/∑ 1/(α_j + rβ_j)` (Cauchy-Schwarz).
Task name: `L_le_opt`. -/
theorem L_le_opt {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) (w : Fin M → ℝ) :
    Lw θ c w ≤ stackSVDLimitW θ c := by
  by_cases h4 : Assumption4 θ c w
  · obtain ⟨hpos, hle1, hgw⟩ := one_le_gW_Lw hc h4
    obtain ⟨j0, hj0⟩ := exists_signal_of_root (gammaTop_of_exists h4.1)
    have hθ : ∃ i, θ i ≠ 0 := ⟨j0, (mul_ne_zero_iff.mp hj0).1⟩
    have hthr : 1 < ∑ i, θ i ^ 4 / c i := by
      have hlt := gW_strictAntiOn hc hθ (Set.left_mem_Icc.mpr zero_le_one)
        (Set.mem_Icc.mpr ⟨hpos.le, hle1⟩) hpos
      rw [gW_zero] at hlt
      linarith
    exact le_stackSVDLimitW hc hthr (Set.mem_Icc.mpr ⟨hpos.le, hle1⟩) hgw
  · unfold Lw
    rw [if_neg h4]
    exact stackSVDLimitW_nonneg θ c

/-- The optimal weights attain `stackSVDLimitW` above the threshold. With `r` the root of
`gW θ c r = 1` the change of variables of `main_paper.tex:1600` gives `γ₁ = 1/(1-r)` and
`p_j = θ_j⁴(1-r)/(c_j + rθ_j²)`, so `∑_j p_j = gW θ c r = 1`. -/
theorem assumption4_and_L_optW_of_thr {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) :
    Assumption4 θ c (optWstack θ c) ∧ Lw θ c (optWstack θ c) = stackSVDLimitW θ c := by
  obtain ⟨hrmem, hgr⟩ := stackSVDLimitW_spec (existsUnique_root hc hthr)
  obtain ⟨hr0, hr1⟩ := hrmem
  set r := stackSVDLimitW θ c with hrdef
  set w := optWstack θ c with hwdef
  have h1r : (0:ℝ) < 1 - r := by linarith
  have h1rne : (1:ℝ) - r ≠ 0 := h1r.ne'
  have hpos : ∀ i : Fin M, (0:ℝ) < θ i ^ 2 + c i := fun i => by
    have := hc i
    positivity
  have hposne : ∀ i : Fin M, θ i ^ 2 + c i ≠ 0 := fun i => (hpos i).ne'
  have hcr : ∀ i : Fin M, (0:ℝ) < c i + r * θ i ^ 2 := fun i => by
    have h := hc i
    have h2 : (0:ℝ) ≤ r * θ i ^ 2 := mul_nonneg hr0.le (sq_nonneg _)
    linarith
  have hcrne : ∀ i : Fin M, c i + r * θ i ^ 2 ≠ 0 := fun i => (hcr i).ne'
  have hw2 : ∀ i : Fin M, w i ^ 2 = θ i ^ 2 / (θ i ^ 2 + c i) := by
    intro i
    rw [hwdef]
    change (θ i / Real.sqrt (θ i ^ 2 + c i)) ^ 2 = _
    rw [div_pow, Real.sq_sqrt (hpos i).le]
  have hw4 : ∀ i : Fin M, w i ^ 4 = (θ i ^ 2 / (θ i ^ 2 + c i)) ^ 2 := fun i => by
    rw [show w i ^ 4 = (w i ^ 2) ^ 2 by ring, hw2 i]
  set g : ℝ := 1 / (1 - r) with hgdef
  have hg1 : (1:ℝ) < g := by
    rw [hgdef, lt_div_iff₀ h1r]
    linarith
  have hden : ∀ i : Fin M,
      g - w i ^ 2 = (c i + r * θ i ^ 2) / ((1 - r) * (θ i ^ 2 + c i)) := by
    intro i
    have h1 := hposne i
    rw [hw2 i, hgdef]
    field_simp
    ring
  set P : Fin M → ℝ := fun i => θ i ^ 4 * (1 - r) / (c i + r * θ i ^ 2) with hPdef
  have hPv : ∀ i : Fin M, P i = θ i ^ 4 * (1 - r) / (c i + r * θ i ^ 2) := fun _ => rfl
  have hPsum : ∑ i, P i = 1 := by
    rw [Finset.sum_congr rfl fun i _ => hPv i]
    exact hgr
  have hPnn : ∀ i : Fin M, 0 ≤ P i := fun i => by
    rw [hPv i]
    exact div_nonneg (mul_nonneg (by positivity) h1r.le) (hcr i).le
  have hroot : IsGammaTop θ w g := by
    constructor
    · have hle : wSqMax w ≤ 1 := by
        unfold wSqMax
        refine Real.iSup_le (fun i => ?_) zero_le_one
        rw [hw2 i, div_le_one (hpos i)]
        have := hc i
        linarith
      linarith
    · change (1:ℝ) + ∑ i, θ i ^ 2 * w i ^ 2 / (w i ^ 2 - g) = 0
      have hterm : ∀ i : Fin M, θ i ^ 2 * w i ^ 2 / (w i ^ 2 - g) = -P i := by
        intro i
        have h1 := hposne i
        have h2 := hcrne i
        rw [show w i ^ 2 - g = -(g - w i ^ 2) by ring, div_neg, hden i, hw2 i, hPv i]
        congr 1
        field_simp
      rw [Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_neg_distrib, hPsum]
      ring
  have hgt : gammaTop θ w = g := gammaTop_eq hroot
  have hex : ∃ g' : ℝ, IsGammaTop θ w g' := ⟨g, hroot⟩
  have hQ : ∀ i : Fin M, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2
      = P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) := by
    intro i
    have h1 := hposne i
    have h2 := hcrne i
    rw [hgt, hden i, hw4 i, hPv i]
    field_simp
  have hR : ∀ i : Fin M,
      gammaTop θ w * (θ i ^ 2 * w i ^ 2 / (gammaTop θ w - w i ^ 2) ^ 2)
        = P i * ((θ i ^ 2 + c i) / (c i + r * θ i ^ 2)) := by
    intro i
    have h1 := hposne i
    have h2 := hcrne i
    rw [hgt, hden i, hw2 i, hPv i, hgdef]
    field_simp
  set R : Fin M → ℝ := fun i => P i * ((θ i ^ 2 + c i) / (c i + r * θ i ^ 2)) with hRdef
  have hRv : ∀ i : Fin M, R i = P i * ((θ i ^ 2 + c i) / (c i + r * θ i ^ 2)) := fun _ => rfl
  have hRnn : ∀ i : Fin M, 0 ≤ R i := fun i => by
    rw [hRv i]
    exact mul_nonneg (hPnn i) (div_nonneg (hpos i).le (hcr i).le)
  have hRpos : 0 < ∑ i, R i := by
    obtain ⟨i0, hi0⟩ := exists_ne_zero_of_thr hthr
    refine Finset.sum_pos' (fun i _ => hRnn i) ⟨i0, Finset.mem_univ _, ?_⟩
    have h4 : (0:ℝ) < θ i0 ^ 4 :=
      lt_of_le_of_ne (by positivity) (Ne.symm (pow_ne_zero 4 hi0))
    have hP0 : 0 < P i0 := by
      rw [hPv i0]
      exact div_pos (mul_pos h4 h1r) (hcr i0)
    rw [hRv i0]
    exact mul_pos hP0 (div_pos (hpos i0) (hcr i0))
  have hRne : (∑ i, R i) ≠ 0 := hRpos.ne'
  have hQle : ∀ i : Fin M,
      P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) ≤ (1 - r) * P i := by
    intro i
    have hfrac : c i * (1 - r) / (c i + r * θ i ^ 2) ≤ 1 - r := by
      rw [div_le_iff₀ (hcr i)]
      nlinarith [mul_nonneg (mul_nonneg h1r.le hr0.le) (sq_nonneg (θ i))]
    calc P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) ≤ P i * (1 - r) :=
          mul_le_mul_of_nonneg_left hfrac (hPnn i)
      _ = (1 - r) * P i := mul_comm _ _
  have hQsum : ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 < 1 := by
    have h1 : ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 ≤ (1 - r) * ∑ i, P i := by
      rw [Finset.mul_sum]
      exact Finset.sum_le_sum fun i _ => by rw [hQ i]; exact hQle i
    rw [hPsum, mul_one] at h1
    linarith
  have hA4 : Assumption4 θ c w := ⟨hex, hQsum⟩
  refine ⟨hA4, ?_⟩
  have hetaeq : eta1 θ c w = 1 - ∑ i, P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) := by
    change (1:ℝ) - ∑ i, c i * w i ^ 4 / (gammaTop θ w - w i ^ 2) ^ 2 = _
    rw [Finset.sum_congr rfl fun i _ => hQ i]
  have hDeneq : LwDen θ w = ∑ i, R i := by
    change gammaTop θ w * ∑ i, θ i ^ 2 * w i ^ 2 / (gammaTop θ w - w i ^ 2) ^ 2 = _
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by rw [hR i, hRv i]
  have hpoint : ∀ i : Fin M,
      P i - P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) = r * R i := by
    intro i
    have h2 := hcrne i
    rw [hRv i, hPv i]
    field_simp
    ring
  have hsplit : (1:ℝ) - ∑ i, P i * (c i * (1 - r) / (c i + r * θ i ^ 2)) = r * ∑ i, R i := by
    have h1 : ∑ i, (P i - P i * (c i * (1 - r) / (c i + r * θ i ^ 2))) = r * ∑ i, R i := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ => hpoint i
    rw [Finset.sum_sub_distrib, hPsum] at h1
    exact h1
  unfold Lw
  rw [if_pos hA4, hetaeq, hDeneq, hsplit, mul_div_assoc, div_self hRne, mul_one]

/-- The paper's optimum: `L(w⋆)` is the unique root in `(0,1)` of `∑_i θ_i⁴(1-x)/(c_i+xθ_i²) = 1`,
which is `stackSVDLimitW`. Below the detectability threshold `∑_i θ_i⁴/c_i > 1` both sides
are `0`. Task name: `L_optW_eq`. -/
theorem L_optW_eq {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    Lw θ c (optWstack θ c) = stackSVDLimitW θ c := by
  by_cases hthr : 1 < ∑ i, θ i ^ 4 / c i
  · exact (assumption4_and_L_optW_of_thr hc hthr).2
  · have hz : stackSVDLimitW θ c = 0 := stackSVDLimitW_eq_zero hc (not_lt.mp hthr)
    refine le_antisymm ?_ ?_
    · rw [hz]
      exact hz ▸ L_le_opt hc (optWstack θ c)
    · rw [hz]
      exact Lw_nonneg θ c _

/-- The threshold equivalence of `main_paper.tex:1602`: `eq:assumption4` holds at `w⋆` exactly
when `∑_j θ_j⁴/c_j > 1`. The companion of `stackSVDLimitW_eq_zero`. -/
theorem assumption4_optW_iff {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    Assumption4 θ c (optWstack θ c) ↔ 1 < ∑ i, θ i ^ 4 / c i := by
  constructor
  · intro h4
    obtain ⟨hpos, hle1, hgw⟩ := one_le_gW_Lw hc h4
    obtain ⟨j0, hj0⟩ := exists_signal_of_root (gammaTop_of_exists h4.1)
    have hθ : ∃ i, θ i ≠ 0 := ⟨j0, (mul_ne_zero_iff.mp hj0).1⟩
    have hlt := gW_strictAntiOn hc hθ (Set.left_mem_Icc.mpr zero_le_one)
      (Set.mem_Icc.mpr ⟨hpos.le, hle1⟩) hpos
    rw [gW_zero] at hlt
    linarith
  · intro hthr
    exact (assumption4_and_L_optW_of_thr hc hthr).1

/-- `L(w)` at a constant nonzero weighting is the unweighted limit. Scale invariance
(`Lw_smul`) plus `Lw_one`. -/
theorem Lw_const {θ c : Fin M → ℝ} {t : ℝ} (ht : t ≠ 0) :
    Lw θ c (fun _ => t) = stackSVDLimit θ c := by
  have h : Lw θ c (fun i => t * (1 : Fin M → ℝ) i) = Lw θ c 1 := Lw_smul ht
  have hfun : (fun i : Fin M => t * (1 : Fin M → ℝ) i) = (fun _ : Fin M => t) := by
    funext i
    simp
  rw [hfun, Lw_one] at h
  exact h

end Scalars

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- Two spiked models with the same data fields are equal; the `Prop` fields are proof
irrelevant. -/
private theorem spikedModel_eq_of_data {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ} {A B : SpikedModel μ n d}
    (hθ : A.θ = B.θ) (hu : A.u = B.u) (hv : A.v = B.v) (hZ : A.Z = B.Z) : A = B := by
  cases A
  cases B
  simp only at hθ hu hv hZ
  subst hθ
  subst hu
  subst hv
  subst hZ
  rfl

section StackW

variable [NeZero M]

/-! ### The weighted stack as a spiked model

`[w_1 X_1; ...; w_M X_M] = θ_w ũ_w vᵀ + E_w` with `θ_w = ‖(θ_i w_i)_i‖₂`, `ũ_w` the unit
vector whose block `i` is `θ_i w_i u_i/θ_w`, and `E_w` the block vector `(w_i E_i)_i`. Every
definition below is the one of `StackSVD.lean` with `θ_i` replaced by `θ_i w_i`.
-/

/-- `∑_i θ_i² w_i²`, the squared signal strength of the weighted stack. -/
noncomputable def stackThetaSqW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) : ℝ :=
  ∑ i, ((m.tbl i).θ * w i) ^ 2

/-- `‖(θ_i w_i)_i‖₂`, the signal strength of the weighted stack. -/
noncomputable def stackThetaW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) : ℝ :=
  Real.sqrt (m.stackThetaSqW w)

/-- The weighted left singular vector in block form: block `i` is `θ_i w_i u_i/θ_w`. -/
noncomputable def stackUSigmaW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    ((i : Fin M) × Fin (n i N)) → ℝ :=
  fun p => (m.tbl p.1).θ * w p.1 * (m.tbl p.1).u N p.2 / m.stackThetaW w

/-- Left singular vector of the weighted stack. When `θ_w = 0` the stacked signal vanishes
and any unit vector serves; the first basis vector of the block index is used, as in
`stackU`. -/
noncomputable def stackUW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  if m.stackThetaSqW w = 0 then
    EuclideanSpace.single (finSigmaFinEquiv ⟨(0 : Fin M), ⟨0, (m.tbl 0).hn N⟩⟩) 1
  else WithLp.toLp 2 fun k => m.stackUSigmaW w N (finSigmaFinEquiv.symm k)

/-- Noise of the weighted stack: block `i` is `w_i Z_i`. This is the block-heteroscedastic
`Σ^{1/2} E` of the paper, written on the unscaled noise. -/
noncomputable def stackZW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => w p.1 * (m.tbl p.1).Z N ω p.2 k)

/-- The weighted left singular vector is a unit vector. Proved, because `stackW` needs it. -/
theorem norm_stackUW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    ‖m.stackUW w N‖ = 1 := by
  rw [stackUW]
  split_ifs with h
  · simp
  · have hpos : 0 < m.stackThetaSqW w :=
      lt_of_le_of_ne (Finset.sum_nonneg fun i _ => sq_nonneg _) (Ne.symm h)
    have hθ2 : m.stackThetaW w ^ 2 = m.stackThetaSqW w := Real.sq_sqrt hpos.le
    have hu : ∀ i : Fin M, ∑ j, ((m.tbl i).u N j) ^ 2 = 1 := by
      intro i
      have hnorm := EuclideanSpace.norm_sq_eq ((m.tbl i).u N)
      rw [(m.tbl i).hu N] at hnorm
      simpa [Real.norm_eq_abs, sq_abs] using hnorm.symm
    have key : ∀ i : Fin M, ∑ j, (m.stackUSigmaW w N ⟨i, j⟩) ^ 2
        = ((m.tbl i).θ * w i) ^ 2 / m.stackThetaSqW w := by
      intro i
      have hterm : ∀ j : Fin (n i N), (m.stackUSigmaW w N ⟨i, j⟩) ^ 2
          = (((m.tbl i).θ * w i) ^ 2 / m.stackThetaSqW w) * ((m.tbl i).u N j) ^ 2 := by
        intro j
        change ((m.tbl i).θ * w i * (m.tbl i).u N j / m.stackThetaW w) ^ 2 = _
        rw [div_pow, hθ2]
        ring
      rw [Finset.sum_congr rfl fun j _ => hterm j, ← Finset.mul_sum, hu i, mul_one]
    rw [EuclideanSpace.norm_eq, Real.sqrt_eq_one]
    have hpt : ∀ k : Fin (∑ i, n i N),
        ‖(WithLp.toLp 2 fun k => m.stackUSigmaW w N (finSigmaFinEquiv.symm k) :
          EuclideanSpace ℝ (Fin (∑ i, n i N))) k‖ ^ 2
          = (m.stackUSigmaW w N (finSigmaFinEquiv.symm k)) ^ 2 := by
      intro k
      simp [Real.norm_eq_abs, sq_abs]
    rw [Finset.sum_congr rfl fun k _ => hpt k,
      Equiv.sum_comp finSigmaFinEquiv.symm (fun p => (m.stackUSigmaW w N p) ^ 2),
      Fintype.sum_sigma (fun p => (m.stackUSigmaW w N p) ^ 2),
      Finset.sum_congr rfl fun i _ => key i, ← Finset.sum_div]
    exact div_self (ne_of_gt hpos)

set_option linter.unusedSectionVars false in
/-- The weighted stacked noise is measurable. Proved, because `stackW` needs it. -/
theorem measurable_stackZW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    Measurable (m.stackZW w N) := by
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun k => ?_
  exact (((measurable_pi_apply k).comp
    ((measurable_pi_apply (finSigmaFinEquiv.symm r).2).comp
      ((m.tbl (finSigmaFinEquiv.symm r).1).hZ N))).const_mul _)

/-- The weighted stack `[w_1 X_1; ...; w_M X_M]` as one rank-one spiked table, with signal
`‖(θ_i w_i)_i‖₂` and the shared right singular vector `v`. Its noise is **not** i.i.d. unless
every `|w_i| = 1`, which is why `SingleTableLaw` does not apply and `HeteroLaw` is needed. -/
noncomputable def stackW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) :
    SpikedModel μ (fun N => ∑ i, n i N) d where
  θ := m.stackThetaW w
  u := m.stackUW w
  v := (m.tbl 0).v
  Z := m.stackZW w
  hθ := Real.sqrt_nonneg _
  hn := m.stack_row_pos
  hd := (m.tbl 0).hd
  hu := m.norm_stackUW w
  hv := (m.tbl 0).hv
  hZ := m.measurable_stackZW w

@[simp]
theorem stackW_theta (m : MultiTableModel μ M n d) (w : Fin M → ℝ) :
    (m.stackW w).θ = Real.sqrt (∑ i, ((m.tbl i).θ * w i) ^ 2) := rfl

@[simp]
theorem stackW_v (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    (m.stackW w).v N = (m.tbl 0).v N := rfl

private theorem stackW_theta_eq (m : MultiTableModel μ M n d) (w : Fin M → ℝ) :
    (m.stackW w).θ = m.stackThetaW w := rfl

private theorem stackW_u_eq (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    (m.stackW w).u N = m.stackUW w N := rfl

private theorem stackW_Z_eq (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    (m.stackW w).Z N ω = m.stackZW w N ω := rfl

/-- Stacking lemma for general weights: block `i` of the weighted stack is `w_i X_i`. The
scaling `θ_i w_i u_i/θ_w` of `stackUW` is chosen for this. -/
theorem stackW_X (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (i : Fin M) (j : Fin (n i N)) (k : Fin (d N)) :
    (m.stackW w).X N ω (finSigmaFinEquiv ⟨i, j⟩) k = w i * (m.tbl i).X N ω j k := by
  have hv : (m.tbl 0).v N = (m.tbl i).v N := m.hv 0 i N
  have hZ : m.stackZW w N ω (finSigmaFinEquiv ⟨i, j⟩) k = w i * (m.tbl i).Z N ω j k := by
    simp [stackZW, Matrix.reindex_apply, Matrix.submatrix_apply]
  have hu : m.stackThetaW w * m.stackUW w N (finSigmaFinEquiv ⟨i, j⟩)
      = (m.tbl i).θ * w i * (m.tbl i).u N j := by
    rw [stackUW]
    split_ifs with hzero
    · have hθi : (m.tbl i).θ * w i = 0 := by
        have hsq := (Finset.sum_eq_zero_iff_of_nonneg
          (fun l (_ : l ∈ Finset.univ) => sq_nonneg ((m.tbl l).θ * w l))).1 hzero i
          (Finset.mem_univ i)
        exact (pow_eq_zero_iff (n := 2) (by norm_num)).1 hsq
      have hst : m.stackThetaW w = 0 := by
        rw [stackThetaW, hzero, Real.sqrt_zero]
      rw [hst, hθi]
      ring
    · have hpos : 0 < m.stackThetaSqW w :=
        lt_of_le_of_ne (Finset.sum_nonneg fun l _ => sq_nonneg _) (Ne.symm hzero)
      have hne : m.stackThetaW w ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hpos)
      rw [PiLp.toLp_apply, Equiv.symm_apply_apply, stackUSigmaW]
      field_simp
  simp only [SpikedModel.X, SpikedModel.E, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
    Matrix.vecMulVec_apply, stackW_theta_eq, stackW_u_eq, stackW_Z_eq, stackW_v]
  rw [hZ, hv, ← mul_assoc, hu]
  ring

/-- At unit weights the weighted stack is the stack of `StackSVD.lean`. -/
theorem stackW_one (m : MultiTableModel μ M n d) : m.stackW 1 = m.stack := by
  have hsq : m.stackThetaSqW 1 = m.stackThetaSq := by
    unfold stackThetaSqW stackThetaSq
    exact Finset.sum_congr rfl fun i _ => by
      show ((m.tbl i).θ * (1 : Fin M → ℝ) i) ^ 2 = (m.tbl i).θ ^ 2
      rw [show (1 : Fin M → ℝ) i = 1 from rfl, mul_one]
  have hθ : m.stackThetaW 1 = m.stackTheta := by
    unfold stackThetaW stackTheta
    rw [hsq]
  have hu : m.stackUW 1 = m.stackU := by
    funext N
    unfold stackUW stackU
    rw [hsq]
    split_ifs with h
    · rfl
    · congr 1
      funext k
      show m.stackUSigmaW 1 N _ = m.stackUSigma N _
      unfold stackUSigmaW stackUSigma
      rw [hθ, show (1 : Fin M → ℝ) _ = 1 from rfl, mul_one]
  have hZ : m.stackZW 1 = m.stackZ := by
    funext N ω
    ext rr k
    simp [stackZW, stackZ, Matrix.reindex_apply, Matrix.submatrix_apply]
  exact spikedModel_eq_of_data hθ hu rfl hZ

/-! ### The heteroscedastic law (black box 2) -/

/-- Gram matrix of the weighted stack. -/
noncomputable def stackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  ((m.stackW w).X N ω)ᵀ * (m.stackW w).X N ω

theorem isHermitian_stackGramW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.stackGramW w N ω).IsHermitian :=
  isHermitian_transpose_mul_self _

/-- Performance of weighted stacksvd: `(v̂_stack(w)ᵀ v)²` in projector form. -/
noncomputable def stackPerfW (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  overlap ((m.stackW w).X N ω) ((m.tbl 0).v N)

omit [NeZero M] in
/-- Black box 2 of `notes/README.md` (BEJ1129 Theorem 2.3, liu2023asymptotic
Theorem 2), stated on the actual random matrices of the model, in the style of
`SingleTableLaw`.

* `align`: `overlap(X_stack(w), v) → L(w)`, the specialization of `eq:weighted_norm` at
  `b = v`. Both branches of the paper are inside `Lw`: above `eq:assumption4` the closed
  form, below it `0`.
* `topSimple`: the top eigenvalue of the Gram matrix is simple a.s., so the projector form
  equals the paper's `⟨v̂, v⟩²` (`thm_stacksvd_weighted_inner`).

The paper's `b`-general half of `eq:weighted_norm` is not a field: no result below consumes
it. See the necessity scan in `notes/archive/thm_stacksvd_weighted.md`. -/
structure HeteroLaw [NeZero M] (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) : Prop where
  align : TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
    (Scalars.Lw (fun i => (m.tbl i).θ) c w)
  topSimple : ∀ N, ∀ᵐ ω ∂(μ N),
    TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω)

/-! ### The constant weighting

A constant nonzero weighting multiplies the stacked matrix by `t`, which changes neither the
overlap nor the simplicity of the top eigenvalue (`Spectral.lean`, `overlap_smul` and
`topSimple_gram_smul`), and `L(w)` is scale invariant (`Scalars.Lw_const`). So the
heteroscedastic law at `w = (t, ..., t)` follows from the `SingleTableLaw` of the unweighted
stack, exactly as at `w = 1`.
-/

/-- The weighted stack at a constant weighting is `t` times the unweighted stack. -/
theorem stackW_const_X (m : MultiTableModel μ M n d) (t : ℝ) (N : ℕ) (ω : Ω N) :
    (m.stackW (fun _ => t)).X N ω = t • m.stack.X N ω := by
  ext r k
  obtain ⟨p, rfl⟩ : ∃ p, finSigmaFinEquiv p = r :=
    ⟨finSigmaFinEquiv.symm r, Equiv.apply_symm_apply _ _⟩
  obtain ⟨i, j⟩ := p
  rw [m.stackW_X (fun _ => t) N ω i j k, Matrix.smul_apply, m.stack_X N ω i j k, smul_eq_mul]

/-! ### The binary weighting and the sub-stack (`cor.2`, model half)

`Lw_binary` is the scalar half of `cor.2`. The model half says that the binary-weighted stack
is the unweighted stack of the sub-collection `S`, up to the zero rows of the discarded
tables. The two matrices have different row counts, so they are not equal; they have the same
Gram matrix, and `overlap` and `TopSimple` see only that (`Spectral.lean`,
`overlap_congr_gram` and `topSimple_congr_gram`).
-/

/-- The tables indexed by `S`, in increasing order, as a `MultiTableModel` over
`Fin S.card`. -/
def restrict (m : MultiTableModel μ M n d) (S : Finset (Fin M)) :
    MultiTableModel μ S.card (fun j => n (S.orderEmbOfFin rfl j)) d where
  tbl j := m.tbl (S.orderEmbOfFin rfl j)
  hv _ _ N := m.hv _ _ N

set_option linter.unusedSectionVars false in
@[simp]
theorem restrict_tbl (m : MultiTableModel μ M n d) (S : Finset (Fin M)) (j : Fin S.card) :
    (m.restrict S).tbl j = m.tbl (S.orderEmbOfFin rfl j) := rfl

set_option linter.unusedSectionVars false in
/-- A sum over `Fin S.card` along the increasing enumeration is a sum over `S`. -/
theorem sum_orderEmb (S : Finset (Fin M)) (f : Fin M → ℝ) :
    ∑ j, f (S.orderEmbOfFin rfl j) = ∑ i ∈ S, f i := by
  conv_rhs => rw [← Finset.map_orderEmbOfFin_univ S rfl]
  rw [Finset.sum_map]
  rfl

/-- Entry of the Gram matrix of the weighted stack: `∑_i w_i² (X_iᵀ X_i)`. The row index of
the stack splits into the block index by `finSigmaFinEquiv`. -/
theorem stackGramW_apply (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (k l : Fin (d N)) :
    (((m.stackW w).X N ω)ᵀ * ((m.stackW w).X N ω)) k l
      = ∑ i, w i ^ 2 * ∑ j, (m.tbl i).X N ω j k * (m.tbl i).X N ω j l := by
  rw [Matrix.mul_apply, ← Equiv.sum_comp finSigmaFinEquiv
      (fun r => ((m.stackW w).X N ω)ᵀ k r * ((m.stackW w).X N ω) r l),
    Fintype.sum_sigma]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Matrix.transpose_apply, m.stackW_X w N ω i j k, m.stackW_X w N ω i j l]
  ring

/-- The binary-weighted stack and the stack of the sub-collection `S` have the same Gram
matrix: the discarded tables contribute zero rows.

`[NeZero S.card]` runs through this group of results and excludes `S = ∅`, which the paper's
maximum over `2^{[M]}` (`cor.2`, `main_paper.tex:442`) formally includes. The exclusion is
right: the empty stack is a `0 × d` matrix, it carries no data, its performance is the junk
value of an empty model, and the paper's own display gives it the value `0`, which never
attains the maximum (audit of the unaudited pieces, 2026-08-31, finding 3).
`exists_binary_tendsto_max_gaussian` reinstates `S` nonempty from
`Scalars.exists_nonempty_eq_max`, so the subset maximum is unaffected. -/
theorem gram_restrict (m : MultiTableModel μ M n d) (S : Finset (Fin M)) [NeZero S.card]
    (N : ℕ) (ω : Ω N) :
    ((m.stackW (fun i => if i ∈ S then 1 else 0)).X N ω)ᵀ *
        ((m.stackW (fun i => if i ∈ S then 1 else 0)).X N ω)
      = ((m.restrict S).stack.X N ω)ᵀ * ((m.restrict S).stack.X N ω) := by
  classical
  ext k l
  rw [m.stackGramW_apply _ N ω k l]
  have h2 := (m.restrict S).stackGramW_apply 1 N ω k l
  rw [(m.restrict S).stackW_one] at h2
  rw [h2]
  simp only [restrict_tbl, Pi.one_apply, one_pow, one_mul]
  rw [sum_orderEmb S
    (fun i => ∑ j, (m.tbl i).X N ω j k * (m.tbl i).X N ω j l)]
  have hterm : ∀ i : Fin M,
      (if i ∈ S then (1 : ℝ) else 0) ^ 2 * ∑ j, (m.tbl i).X N ω j k * (m.tbl i).X N ω j l
        = if i ∈ S then ∑ j, (m.tbl i).X N ω j k * (m.tbl i).X N ω j l else 0 := by
    intro i
    by_cases hi : i ∈ S <;> simp [hi]
  rw [Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_ite_mem, Finset.univ_inter]

/-! ### `cor.2` for Gaussian noise

The model half above consumes `(m.restrict S).stack.SingleTableLaw`, which
`SpikedModel.singleTableLaw_of_gaussian` supplies once the sub-collection is known to be
jointly Gaussian. That step is the marginalization of a product measure along the coordinates
in `S`, in two standard pieces: `measurePreserving_finsetRestrict` (checked on rectangles with
`Measure.pi_eq`; the coordinates outside `S` contribute `ν i univ = 1`) and
`measurePreserving_piCongrLeft` (reindex `S` by `Fin S.card` along the increasing
enumeration). Finding 5 of `notes/archive/audit_global_2026-08-31.md`. -/

/-- Marginalization of a product of probability measures along the coordinates of a finite
subset: `f ↦ (f i)_{i ∈ S}` carries `Measure.pi ν` to `Measure.pi (ν ∘ Subtype.val)`. -/
theorem _root_.StackedSVD.measurePreserving_finsetRestrict {ι : Type*} [Fintype ι]
    {X : ι → Type*} [∀ i, MeasurableSpace (X i)] (ν : ∀ i, Measure (X i))
    [∀ i, IsProbabilityMeasure (ν i)] (S : Finset ι) :
    MeasurePreserving (fun f : ∀ i, X i => fun j : S => f (j : ι)) (Measure.pi ν)
      (Measure.pi fun j : S => ν (j : ι)) := by
  classical
  have hmeas : Measurable (fun f : ∀ i, X i => fun j : S => f (j : ι)) :=
    measurable_pi_lambda _ fun j => measurable_pi_apply (j : ι)
  refine ⟨hmeas, (Measure.pi_eq fun s hs => ?_).symm⟩
  rw [Measure.map_apply hmeas (MeasurableSet.univ_pi hs)]
  have hpre : (fun f : ∀ i, X i => fun j : S => f (j : ι)) ⁻¹' Set.univ.pi s
      = Set.univ.pi fun i => if h : i ∈ S then s ⟨i, h⟩ else Set.univ := by
    ext f
    simp only [Set.mem_preimage, Set.mem_univ_pi]
    constructor
    · intro hf i
      by_cases h : i ∈ S
      · rw [dif_pos h]
        exact hf ⟨i, h⟩
      · rw [dif_neg h]
        exact Set.mem_univ _
    · intro hf j
      have h := hf (j : ι)
      rwa [dif_pos j.2] at h
  rw [hpre, Measure.pi_pi]
  calc ∏ i : ι, ν i (if h : i ∈ S then s ⟨i, h⟩ else Set.univ)
      = ∏ i ∈ S, ν i (if h : i ∈ S then s ⟨i, h⟩ else Set.univ) := by
        refine (Finset.prod_subset (Finset.subset_univ S) fun i _ hi => ?_).symm
        rw [dif_neg hi, measure_univ]
    _ = ∏ j : S, ν (j : ι) (if h : (j : ι) ∈ S then s ⟨(j : ι), h⟩ else Set.univ) :=
        (Finset.prod_coe_sort S _).symm
    _ = ∏ j : S, ν (j : ι) (s j) := by
        refine Finset.prod_congr rfl fun j _ => ?_
        rw [dif_pos j.2]

omit [NeZero M] in
/-- **The sub-collection is jointly Gaussian.** This is the step that `notes/FLAGGED.md` left
open, so that the inventory row `cor.2` now has a Gaussian discharge. -/
theorem restrict_jointGaussianNoise (m : MultiTableModel μ M n d) (S : Finset (Fin M))
    (hG : m.JointGaussianNoise) : (m.restrict S).JointGaussianNoise := by
  classical
  intro N
  have hres := StackedSVD.measurePreserving_finsetRestrict
    (fun i : Fin M => gaussianMatrix (n i N) (d N)) S
  have hpi := measurePreserving_piCongrLeft
    (fun j : {x // x ∈ S} => gaussianMatrix (n (j : Fin M) N) (d N))
    (S.orderIsoOfFin rfl).toEquiv
  exact ((hpi.symm (MeasurableEquiv.piCongrLeft _ _)).comp hres).fun_comp_hasLaw (hG N)

end StackW

end MultiTableModel

end StackedSVD

namespace StackedSVD.Scalars

variable {M : ℕ}

/-! ### The best binary weighting against the optimal weighting (`main_paper.tex:442`)

`binaryStackSVDLimitMax` is defined in `Scalars.lean`. The comparison with the optimally
weighted stacksvd limit needs `Lw_binary` and `L_le_opt`, which live here: a binary weighting
is one weighting among all of them, so the simplex optimization already covers it.
-/

/-- Any binary weighting is at most the optimally weighted stacksvd limit. -/
theorem binaryStackSVDLimit_le_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (S : Finset (Fin M)) :
    binaryStackSVDLimit S θ c ≤ stackSVDLimitW θ c := by
  rw [← Lw_binary S θ c]
  exact L_le_opt hc _

/-- Optimally binary-weighted stacksvd never beats optimally weighted stacksvd. With
`stackSVDLimit_le_binaryStackSVDLimitMax` this places the best binary choice between
unweighted stacksvd and the optimal weighting. -/
theorem binaryStackSVDLimitMax_le_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    binaryStackSVDLimitMax θ c ≤ stackSVDLimitW θ c :=
  binaryStackSVDLimitMax_le fun S => binaryStackSVDLimit_le_stackSVDLimitW hc S

end StackedSVD.Scalars

