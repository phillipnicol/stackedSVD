/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Example
import StackedSVD.RankR.SingleWeight.Existence
import StackedSVD.RMT.Het.R6het

/-!
# The weight regimes of the single-weight witness (F18b, unit U1)

The scalars behind `hlaw` of `prop_singleweight_suboptimality_of_law`
(`RankR/SingleWeight/Existence.lean`): which of the two secular roots of the witness
(`θ_0 = 8/5`, `c_0 = 1`, `R_1 = e_1`, `R_2 = e_2`) is detectable at a weight vector `w`,
the `EigSep` structure when both are, the bridge from `¬ DetectableEx` to
`¬ Assumption4` of the one-spike profile when one is not, and the value of `swLimitEx` at a
tie. Plan: `notes/archive/F18b_plan.md`, section 2 (U1) and section 6.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

/-- The order map of a pair of indices: `ordPair l₀ l₁ 0 = l₀`, `ordPair l₀ l₁ 1 = l₁`. -/
def ordPair (l₀ l₁ : Fin 2) : Fin 2 → Fin 2 := fun l => if l = 0 then l₀ else l₁

/-- U1.0: `wSqMax` of two weights is the square of the larger one. -/
theorem wSqMax_eq_of_le (w : Fin 2 → ℝ) {l₀ l₁ : Fin 2} (hne : l₀ ≠ l₁) (h0 : 0 ≤ w l₁)
    (hle : w l₁ ≤ w l₀) : Scalars.wSqMax w = w l₀ ^ 2 := by
  refine le_antisymm ?_ (Scalars.le_wSqMax w l₀)
  unfold Scalars.wSqMax
  refine ciSup_le fun i => ?_
  rcases eq_or_ne i l₀ with hi | hi
  · rw [hi]
  · have hi₁ : i = l₁ := by
      revert hne hi
      fin_cases i <;> fin_cases l₀ <;> fin_cases l₁ <;> decide
    rw [hi₁]
    exact pow_le_pow_left₀ h0 hle 2

/-- Helper for U1.1: at a table `l₀` with `0 < w l₀`, every index `i` with `0 ≤ w i ≤ w l₀`
contributes at most `c₀/θ₀⁴` to the threshold sum at `gammaEx θ₀ w l₀`. Both the numerator
(`w i⁴ ≤ w l₀⁴`) and the denominator (`gammaEx θ₀ w l₀ - w i² ≥ w l₀² θ₀²`) are controlled by
`w i ≤ w l₀`, and the resulting bound `c₀ w l₀⁴/(w l₀² θ₀²)² = c₀/θ₀⁴` is exact at `i = l₀`. -/
private theorem swtermEx_bound {θ₀ c₀ : ℝ} (hc : 0 < c₀) (hθ2pos : 0 < θ₀ ^ 2) (w : Fin 2 → ℝ)
    {l₀ : Fin 2} (hw0 : 0 < w l₀) (i : Fin 2) (hwi0 : 0 ≤ w i) (hwile : w i ≤ w l₀) :
    c₀ * w i ^ 4 / (gammaEx θ₀ w l₀ - w i ^ 2) ^ 2 ≤ c₀ / θ₀ ^ 4 := by
  have hθ0ne : θ₀ ≠ 0 := by
    intro h
    rw [h] at hθ2pos
    norm_num at hθ2pos
  have hsq_le : w i ^ 2 ≤ w l₀ ^ 2 := by
    nlinarith [mul_le_mul hwile hwile hwi0 (le_trans hwi0 hwile)]
  have hw4 : w i ^ 4 ≤ w l₀ ^ 4 := by
    nlinarith [mul_le_mul hsq_le hsq_le (sq_nonneg (w i)) (le_trans (sq_nonneg (w i)) hsq_le)]
  have hD : w l₀ ^ 2 * θ₀ ^ 2 ≤ gammaEx θ₀ w l₀ - w i ^ 2 := by
    unfold gammaEx; nlinarith
  have hDpos : 0 < w l₀ ^ 2 * θ₀ ^ 2 := mul_pos (by positivity) hθ2pos
  have hDsq : (w l₀ ^ 2 * θ₀ ^ 2) ^ 2 ≤ (gammaEx θ₀ w l₀ - w i ^ 2) ^ 2 := by
    nlinarith [mul_le_mul hD hD hDpos.le (le_trans hDpos.le hD)]
  have step1 : c₀ * w i ^ 4 / (gammaEx θ₀ w l₀ - w i ^ 2) ^ 2
      ≤ c₀ * w l₀ ^ 4 / (gammaEx θ₀ w l₀ - w i ^ 2) ^ 2 :=
    div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left hw4 hc.le)
      (sq_nonneg (gammaEx θ₀ w l₀ - w i ^ 2))
  have step2 : c₀ * w l₀ ^ 4 / (gammaEx θ₀ w l₀ - w i ^ 2) ^ 2
      ≤ c₀ * w l₀ ^ 4 / (w l₀ ^ 2 * θ₀ ^ 2) ^ 2 :=
    div_le_div_of_nonneg_left (by positivity) (pow_pos hDpos 2) hDsq
  have step3 : c₀ * w l₀ ^ 4 / (w l₀ ^ 2 * θ₀ ^ 2) ^ 2 = c₀ / θ₀ ^ 4 := by
    field_simp [hw0.ne', hθ0ne]
  linarith [step1, step2, step3.le]

/-- U1.1: the root of the larger weight is always detectable when `2 c₀ < θ₀⁴`. -/
theorem detectableEx_of_le {θ₀ c₀ : ℝ} (hc : 0 < c₀) (hsup : 2 * c₀ < θ₀ ^ 4)
    (w : Fin 2 → ℝ) {l₀ l₁ : Fin 2} (hne : l₀ ≠ l₁) (h1 : 0 < w l₁) (hle : w l₁ ≤ w l₀) :
    DetectableEx θ₀ c₀ w l₀ := by
  have hw0pos : 0 < w l₀ := lt_of_lt_of_le h1 hle
  have hθ4pos : 0 < θ₀ ^ 4 := by linarith
  have hθ2pos : 0 < θ₀ ^ 2 := by nlinarith [sq_nonneg θ₀, sq_nonneg (θ₀ ^ 2)]
  refine ⟨?_, ?_⟩
  · rw [wSqMax_eq_of_le w hne h1.le hle]
    unfold gammaEx
    nlinarith [mul_pos (pow_pos hw0pos 2) hθ2pos]
  · have hlt2 : c₀ / θ₀ ^ 4 + c₀ / θ₀ ^ 4 < 1 := by
      have heq : c₀ / θ₀ ^ 4 + c₀ / θ₀ ^ 4 = (2 * c₀) / θ₀ ^ 4 := by ring
      rw [heq, div_lt_one hθ4pos]
      exact hsup
    fin_cases l₀ <;> fin_cases l₁ <;> simp_all [Fin.sum_univ_two] <;>
      linarith [swtermEx_bound hc hθ2pos w hw0pos _ hw0pos.le le_rfl,
                swtermEx_bound hc hθ2pos w hw0pos _ h1.le hle, hlt2]

/-- U1.2: the roots are ordered as the weights are. -/
theorem gammaEx_lt_of_lt {θ₀ : ℝ} (hθ : 0 < θ₀) (w : Fin 2 → ℝ) {l₀ l₁ : Fin 2}
    (h0 : 0 ≤ w l₁) (hlt : w l₁ < w l₀) : gammaEx θ₀ w l₁ < gammaEx θ₀ w l₀ := by
  unfold gammaEx
  have hsq : w l₁ ^ 2 < w l₀ ^ 2 := by nlinarith
  have hpos : 0 < 1 + θ₀ ^ 2 := by positivity
  exact mul_lt_mul_of_pos_right hsq hpos

/-- `secMat` of the instance is diagonal: the `(k, k')` entry is the `k`-th secular summand
when `k = k'`, and `0` off the diagonal, since `sigMat_ex` is supported on `i = k = k'`. -/
theorem secMat_ex_apply (θ₀ γ : ℝ) (w : Fin 2 → ℝ) (k k' : Fin 2) :
    secMat (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w γ k k'
      = if k = k' then w k ^ 2 * θ₀ ^ 2 / (γ - w k ^ 2) else 0 := by
  simp only [secMat, Matrix.sum_apply, Matrix.smul_apply, sigMat_ex, Fin.sum_univ_two]
  fin_cases k <;> fin_cases k' <;> simp <;> ring

/-- At `γ = gammaEx θ₀ w j`, the `j`-th diagonal entry of `secMat` is exactly `1`, since
`γ - w j² = w j² θ₀²` cancels the summand's own denominator. -/
theorem secMat_ex_diag_eq_one (θ₀ : ℝ) (w : Fin 2 → ℝ) {j : Fin 2} (hwj : 0 < w j)
    (hθ0 : θ₀ ≠ 0) :
    w j ^ 2 * θ₀ ^ 2 / (gammaEx θ₀ w j - w j ^ 2) = 1 := by
  have hsub : gammaEx θ₀ w j - w j ^ 2 = w j ^ 2 * θ₀ ^ 2 := by unfold gammaEx; ring
  rw [hsub]
  exact div_self (by positivity)

/-- At `γ = gammaEx θ₀ w j`, `1 - secMat` is singular: its `j`-th diagonal entry vanishes. -/
theorem det_one_sub_secMat_zero (θ₀ : ℝ) (w : Fin 2 → ℝ) {j : Fin 2} (hwj : 0 < w j)
    (hθ0 : θ₀ ≠ 0) :
    (1 - secMat (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w
        (gammaEx θ₀ w j)).det = 0 := by
  have hdiag := secMat_ex_diag_eq_one θ₀ w hwj hθ0
  fin_cases j <;> simp_all [Matrix.det_fin_two, Matrix.sub_apply, secMat_ex_apply]

/-- At `γ = gammaEx θ₀ w j`, `secMat` fixes the `j`-th standard basis vector: column `j` of a
diagonal matrix with a `1` at `(j, j)` is `e_j` itself. -/
theorem secMat_ex_mulVec_single (θ₀ : ℝ) (w : Fin 2 → ℝ) {j : Fin 2} (hwj : 0 < w j)
    (hθ0 : θ₀ ≠ 0) :
    secMat (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w
        (gammaEx θ₀ w j)
      *ᵥ WithLp.ofLp (EuclideanSpace.single j (1 : ℝ))
      = WithLp.ofLp (EuclideanSpace.single j (1 : ℝ)) := by
  have hdiag := secMat_ex_diag_eq_one θ₀ w hwj hθ0
  funext k
  fin_cases k <;> fin_cases j <;> simp_all [Matrix.mulVec, dotProduct, secMat_ex_apply]

/-- U1.3: `EigSep` on the witness scalars, with the roots in decreasing order. -/
theorem eigSep_of_detectable {θ₀ c₀ : ℝ} (hθ : 0 < θ₀) (w : Fin 2 → ℝ)
    (hw : ∀ i, 0 < w i) {l₀ l₁ : Fin 2} (hlt : w l₁ < w l₀)
    (h0 : DetectableEx θ₀ c₀ w l₀) (h1 : DetectableEx θ₀ c₀ w l₁) :
    EigSep (r := 2) (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w
      (fun _ => c₀) (fun l => gammaEx θ₀ w (ordPair l₀ l₁ l))
      (fun l => EuclideanSpace.single (ordPair l₀ l₁ l) 1) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro l
    fin_cases l
    · exact ⟨h0.1, det_one_sub_secMat_zero θ₀ w (hw l₀) hθ.ne'⟩
    · exact ⟨h1.1, det_one_sub_secMat_zero θ₀ w (hw l₁) hθ.ne'⟩
  · intro a b hab
    fin_cases a <;> fin_cases b <;>
      first
      | exact absurd hab (by decide)
      | exact gammaEx_lt_of_lt hθ w (hw l₁).le hlt
  · intro l
    fin_cases l
    · exact ⟨by simp, secMat_ex_mulVec_single θ₀ w (hw l₀) hθ.ne'⟩
    · exact ⟨by simp, secMat_ex_mulVec_single θ₀ w (hw l₁) hθ.ne'⟩
  · intro l
    fin_cases l
    · exact h0.2
    · exact h1.2

/-- U1.4a: `swLimit` does not see a permutation of the components. -/
theorem swLimit_comp_equiv {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) (σ : Equiv.Perm (Fin r)) :
    swLimit θ R w c (γ ∘ σ) (z ∘ σ) = swLimit θ R w c γ z := by
  unfold swLimit
  simp only [Function.comp_apply]
  exact Equiv.sum_comp σ (fun l => swTerm θ R w c (γ l) (z l))

/-- U1.4b: the closed form of the instance is `swLimit` read in the decreasing order. -/
theorem swLimitEx_eq_swLimit_ord (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (hθ : 0 < θ₀)
    (hw : ∀ i, 0 < w i) {l₀ l₁ : Fin 2} (hne : l₀ ≠ l₁)
    (h0 : DetectableEx θ₀ c₀ w l₀) (h1 : DetectableEx θ₀ c₀ w l₁) :
    swLimitEx θ₀ c₀ w
      = swLimit (r := 2) (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀)
          (Rone (RankR.Example.Rex 0)) w (fun _ => c₀)
          (fun l => gammaEx θ₀ w (ordPair l₀ l₁ l))
          (fun l => EuclideanSpace.single (ordPair l₀ l₁ l) 1) := by
  fin_cases l₀ <;> fin_cases l₁ <;> (try exact absurd rfl hne)
  · change swLimitEx θ₀ c₀ w
        = swLimit (r := 2) (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀)
            (Rone (RankR.Example.Rex 0)) w (fun _ => c₀)
            (fun l => gammaEx θ₀ w (ordPair (0 : Fin 2) 1 l))
            (fun l => EuclideanSpace.single (ordPair (0 : Fin 2) 1 l) 1)
    have hidx : ∀ l : Fin 2, ordPair (0 : Fin 2) 1 l = l := by
      intro l; fin_cases l <;> decide
    have hγ : (fun l => gammaEx θ₀ w (ordPair (0 : Fin 2) 1 l)) = gammaEx θ₀ w := by
      funext l; rw [hidx l]
    have hz : (fun l => EuclideanSpace.single (ordPair (0 : Fin 2) 1 l) (1 : ℝ))
        = (fun l => EuclideanSpace.single l (1 : ℝ)) := by
      funext l; rw [hidx l]
    rw [hγ, hz]
    exact swLimitEx_eq_swLimit θ₀ c₀ w hθ hw h0 h1
  · change swLimitEx θ₀ c₀ w
        = swLimit (r := 2) (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀)
            (Rone (RankR.Example.Rex 0)) w (fun _ => c₀)
            (fun l => gammaEx θ₀ w (ordPair (1 : Fin 2) 0 l))
            (fun l => EuclideanSpace.single (ordPair (1 : Fin 2) 0 l) 1)
    have hidx : ∀ l : Fin 2, ordPair (1 : Fin 2) 0 l = Equiv.swap (0 : Fin 2) 1 l := by
      intro l; fin_cases l <;> decide
    have hγ : (fun l => gammaEx θ₀ w (ordPair (1 : Fin 2) 0 l))
        = (gammaEx θ₀ w) ∘ ⇑(Equiv.swap (0 : Fin 2) 1) := by
      funext l; simp only [Function.comp_apply]; rw [hidx l]
    have hz : (fun l => EuclideanSpace.single (ordPair (1 : Fin 2) 0 l) (1 : ℝ))
        = (fun l => EuclideanSpace.single l (1 : ℝ)) ∘ ⇑(Equiv.swap (0 : Fin 2) 1) := by
      funext l; simp only [Function.comp_apply]; rw [hidx l]
    rw [hγ, hz,
      swLimit_comp_equiv (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w (fun _ => c₀)
        (gammaEx θ₀ w) (fun l => EuclideanSpace.single l 1) (Equiv.swap 0 1)]
    exact swLimitEx_eq_swLimit θ₀ c₀ w hθ hw h1 h0

/-- U1.5: an undetectable root of the instance means `eq:assumption4` fails for the
one-table profile `θ^{(l₁)}`. -/
theorem not_assumption4_of_not_detectable {θ₀ c₀ : ℝ} (w : Fin 2 → ℝ) (l₁ : Fin 2)
    (hnd : ¬ DetectableEx θ₀ c₀ w l₁) :
    ¬ Scalars.Assumption4 (fun i => if i = l₁ then θ₀ else 0) (fun _ => c₀) w := by
  intro hA
  apply hnd
  obtain ⟨⟨g, hg⟩, hsum⟩ := hA
  have hgt : Scalars.wSqMax w < g := hg.1
  have hlt0 : w 0 ^ 2 < g := lt_of_le_of_lt (Scalars.le_wSqMax w 0) hgt
  have hlt1 : w 1 ^ 2 < g := lt_of_le_of_lt (Scalars.le_wSqMax w 1) hgt
  have hne0 : w 0 ^ 2 - g ≠ 0 := sub_ne_zero.mpr hlt0.ne
  have hne1 : w 1 ^ 2 - g ≠ 0 := sub_ne_zero.mpr hlt1.ne
  have h2 := hg.2
  unfold Scalars.secular at h2
  rw [Fin.sum_univ_two] at h2
  have hgeq : g = gammaEx θ₀ w l₁ := by
    unfold gammaEx
    fin_cases l₁ <;> simp_all <;> field_simp [hne0, hne1] at h2 <;> nlinarith [h2]
  have hgammaTop : Scalars.gammaTop (fun i => if i = l₁ then θ₀ else 0) w = gammaEx θ₀ w l₁ := by
    rw [Scalars.gammaTop_eq hg, hgeq]
  refine ⟨by rw [← hgeq]; exact hgt, ?_⟩
  rw [← hgammaTop]
  exact hsum

/-- U1.6: at a tie the closed form is the unweighted stacksvd value `2 β²(θ₀, 2 c₀)`. -/
theorem swLimitEx_tie {θ₀ c₀ a : ℝ} (hθ : 0 < θ₀) (hc : 0 < c₀) (ha : 0 < a)
    (hsup : 2 * c₀ < θ₀ ^ 4) :
    swLimitEx θ₀ c₀ (fun _ => a) = RankR.Example.stacksvdEx θ₀ c₀ 0 := by
  have hd0 : DetectableEx θ₀ c₀ (fun _ => a) 0 :=
    detectableEx_of_le hc hsup (fun _ => a) (show (0 : Fin 2) ≠ 1 by decide) ha (le_refl a)
  have hd1 : DetectableEx θ₀ c₀ (fun _ => a) 1 :=
    detectableEx_of_le hc hsup (fun _ => a) (show (1 : Fin 2) ≠ 0 by decide) ha (le_refl a)
  have hstack : RankR.Example.stacksvdEx θ₀ c₀ 0 = 2 * betaSq θ₀ (2 * c₀) :=
    RankR.Example.stacksvdEx_zero_s hθ
  have hbeta : betaSq θ₀ (2 * c₀) = (θ₀ ^ 4 - 2 * c₀) / (θ₀ ^ 2 * (1 + θ₀ ^ 2)) := by
    unfold betaSq
    rw [if_pos hsup]
    congr 1
    ring
  have hterm : ∀ l : Fin 2,
      swTermEx θ₀ c₀ (fun _ => a) l = (θ₀ ^ 4 - 2 * c₀) / (θ₀ ^ 2 * (1 + θ₀ ^ 2)) := by
    intro l
    have hsub : gammaEx θ₀ (fun _ : Fin 2 => a) l - (fun _ : Fin 2 => a) (other l) ^ 2
        = a ^ 2 * θ₀ ^ 2 := by
      fin_cases l <;> simp [gammaEx] <;> ring
    simp only [swTermEx, hsub]
    field_simp [ha.ne', hθ.ne']; ring
  unfold swLimitEx
  rw [if_pos hd0, if_pos hd1, hterm 0, hterm 1, hstack, hbeta]
  ring

end SingleWeight

end StackedSVD
