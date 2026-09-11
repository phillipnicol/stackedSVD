/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R5het
import StackedSVD.RMT.Het.Simplicity
import StackedSVD.RMT.R6

/-!
# Item H12: the subcritical half of `HeteroLaw.align` under `HeteroEdge`

Task H12 of `notes/archive/plan_heterolaw_A.md` (sections 2.4 and 3.7). Below `eq:assumption4` the
weighted stack has no outlier and the overlap `stackPerfW` tends to `Scalars.Lw θ c w = 0`.

Route (plan section 3.7, the `n` side of the column split, R6' step 2 only). Write
`λ := gramLamMax (X_w)`, `W₀' := m.W0het`, `q := m.qHet`, so `X_w X_wᵀ = W₀' + q qᵀ`.

1. Scalars. When `Assumption4` fails, the limiting secular function `1 + F` is positive on
   `(bHet, ∞)` (`one_add_F_pos_of_not_assumption4`): with `γ(z) = -1/s(z)`, `1 + F(z)` is
   `(γ/z) · secular θ w γ(z)`, and `γ(z) > γ₁` (branch failure) or no root exists at all.
   `F'(z) → +∞` as `z ↓ bHet` (`FhetDeriv_tendsto_atTop`, from `sPhysDeriv_tendsto_atTop`),
   and `F(z) < 0` on `(bHet, ∞)` (`Fhet_neg`).
2. Deterministic. `R6.overlap_q_le` applied to `X_wᵀ` gives, on the simple event,
   `‖P q‖² ≤ 1 / qform2 W₀' z₀ q` for every `z₀ ≥ λ`, with `P` the top projector of
   `X_w X_wᵀ`; `Het.topProj_transpose_eq` turns `stackPerfW` into `‖P q‖² / λ`. The degenerate
   case `λ = λ_max(W₀')` is inside `R6.overlap_q_le`, so no conditional symmetry is needed.
   `λ ≥ q ⬝ᵥ q` (Rayleigh at `v`) and `q ⬝ᵥ q ≥ -qform W₀' z₁ q` for `z₁ ≥ λ_max(W₀') + 1`
   keep `λ` away from `0`, so no lower edge for `W₀'` is needed and `HeteroEdge` stays
   one-sided.
3. Probability. Five bad families: the edge event, the secular form at `z₀` (positive limit,
   so `λ ≤ z₀` by `gramLamMax_le_of_secular_pos_het`), the squared form at `z₀` (limit
   `F'(z₀)`), the secular form at `z₁ = bHet + 2` (limit `1 + F(z₁) < 1`, so `q ⬝ᵥ q` is
   bounded below), and the null set where the top eigenvalue is not simple (item H10).

Numeric check (a session script, not kept): seed 20260901, `θ = (0.5, 0.4, 0.6)`,
`c = (0.8, 1.5, 2.2)`, `w = (1, 0.6, 1.4)`, `bHet = 13.6735`, `Assumption4` false: at
`d = 400` and `d = 800` (4 draws each) `λ = 13.59 ± 0.15` and `13.59 ± 0.11`, overlap
`0.008 ± 0.011` and `0.008 ± 0.008`, and every link of the chain above holds on every draw.

Paper: `main_paper.tex` line 462 (`thm:stacksvd_weighted`), line 1409 (`eq:assumption4`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### Scalars: the sign of `1 + F` and the blow-up of `F'` below threshold -/

namespace MPhet

open Scalars

variable {M : ℕ} {θ c w : Fin M → ℝ}

/-- `g_i(z) < 0` on `(bHet, ∞)`. -/
theorem ghet_neg (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hz : bHet c w < z)
    (i : Fin M) : ghet c w i z < 0 := by
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hp := one_add_mul_sPhys_pos hc hw hz i
  unfold ghet
  exact div_neg_of_neg_of_pos (by norm_num) (mul_pos hz0 hp)

/-- `Ψ(z) < 0` on `(bHet, ∞)`: every term is nonpositive and a table with `w_i ≠ 0`
contributes a negative one. -/
theorem Psihet_neg (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hz : bHet c w < z) :
    Psihet c w z < 0 := by
  obtain ⟨i₀, hi₀⟩ := hw
  have hw : ∃ i, w i ≠ 0 := ⟨i₀, hi₀⟩
  have hpos : 0 < ∑ i, -(c i * w i ^ 2 * ghet c w i z) := by
    refine Finset.sum_pos' (fun i _ => ?_) ⟨i₀, Finset.mem_univ _, ?_⟩
    · exact neg_nonneg.mpr (mul_nonpos_of_nonneg_of_nonpos
        (mul_nonneg (hc i).le (sq_nonneg _)) (ghet_neg hc hw hz i).le)
    · have h1 : 0 < c i₀ * w i₀ ^ 2 := mul_pos (hc i₀) (by positivity)
      exact neg_pos.mpr (mul_neg_of_pos_of_neg h1 (ghet_neg hc hw hz i₀))
  rw [Finset.sum_neg_distrib] at hpos
  unfold Psihet
  linarith

/-- `Φ(z) ≤ 0` on `(bHet, ∞)`. -/
theorem Phihet_nonpos (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hz : bHet c w < z) :
    Phihet θ c w z ≤ 0 := by
  unfold Phihet
  exact Finset.sum_nonpos fun i _ =>
    mul_nonpos_of_nonneg_of_nonpos (by positivity) (ghet_neg hc hw hz i).le

/-- `F(z) < 0` on `(bHet, ∞)`. -/
theorem Fhet_neg (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hz : bHet c w < z) :
    Fhet θ c w z < 0 := by
  unfold Fhet
  linarith [Phihet_nonpos (θ := θ) hc hw hz, Psihet_neg hc hw hz]

/-- Without a root of the paper's secular function above `max_i w_i²`, that function is
positive there: a nonpositive value would give a root by the intermediate value theorem
(`f ≥ 0` at `max_i w_i² + ∑ θ_j² w_j² + 1`) or by monotonicity. -/
theorem secular_pos_of_no_root {θ w : Fin M → ℝ} (hno : ¬ ∃ g, IsGammaTop θ w g) {γ : ℝ}
    (hγ : wSqMax w < γ) : 0 < secular θ w γ := by
  by_cases hθw : ∃ j, θ j * w j ≠ 0
  · by_contra hle
    rw [not_lt] at hle
    set S := ∑ j, θ j ^ 2 * w j ^ 2 with hSdef
    have hS0 : 0 ≤ S := Finset.sum_nonneg fun j _ => by positivity
    have hS1 : (0 : ℝ) < S + 1 := by linarith
    set b := wSqMax w + S + 1 with hb
    have hfb : 0 ≤ secular θ w b := by
      have hterm : ∀ j : Fin M,
          0 ≤ θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1) := by
        intro j
        have hb2 := le_wSqMax w j
        have hle' : S + 1 ≤ b - w j ^ 2 := by rw [hb]; linarith
        have h := div_le_div_of_nonneg_left (a := θ j ^ 2 * w j ^ 2) (by positivity) hS1 hle'
        have heq : θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b)
            = -(θ j ^ 2 * w j ^ 2 / (b - w j ^ 2)) := by
          rw [show w j ^ 2 - b = -(b - w j ^ 2) by ring, div_neg]
        rw [heq]
        linarith
      have hpos : 0 ≤ ∑ j, (θ j ^ 2 * w j ^ 2 / (w j ^ 2 - b) + θ j ^ 2 * w j ^ 2 / (S + 1)) :=
        Finset.sum_nonneg fun j _ => hterm j
      rw [Finset.sum_add_distrib, ← Finset.sum_div, ← hSdef] at hpos
      have hSS : S / (S + 1) ≤ 1 := by rw [div_le_one hS1]; linarith
      unfold secular
      linarith
    rcases le_or_gt γ b with hγb | hbγ
    · obtain ⟨g, hgmem, hgval⟩ :=
        intermediate_value_Icc hγb (secular_continuousOn (b := b) hγ) ⟨hle, hfb⟩
      exact hno ⟨g, lt_of_lt_of_le hγ hgmem.1, hgval⟩
    · have hbW : wSqMax w < b := by rw [hb]; linarith
      have hmono := (secular_strictMonoOn hθw).monotoneOn (Set.mem_Ioi.2 hbW)
        (Set.mem_Ioi.2 hγ) hbγ.le
      exact hno ⟨γ, hγ, le_antisymm hle (le_trans hfb hmono)⟩
  · simp only [not_exists, not_not] at hθw
    have h1 : secular θ w γ = 1 := by
      unfold secular
      rw [Finset.sum_eq_zero, add_zero]
      intro j _
      rw [show θ j ^ 2 * w j ^ 2 = (θ j * w j) ^ 2 by ring, hθw j]
      simp
    rw [h1]
    exact one_pos

/-- **The sign fact of plan section 2.4 below threshold.** When `eq:assumption4` fails,
`1 + F > 0` on `(bHet, ∞)`: either no root `γ₁` exists, or `-1/γ₁ ≤ s⋆ < s(z)` gives
`γ(z) > γ₁` and the secular function is positive above its root. -/
theorem one_add_F_pos_of_not_assumption4 (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (h4 : ¬ Assumption4 θ c w) {z : ℝ} (hz : bHet c w < z) : 0 < 1 + Fhet θ c w z := by
  rw [one_add_F_eq hc hw hz]
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hγW : wSqMax w < gamHet c w z := wSqMax_lt_gamHet hc hw hz
  have hγ : 0 < gamHet c w z := lt_of_le_of_lt (wSqMax_nonneg w) hγW
  refine mul_pos (div_pos hγ hz0) ?_
  by_cases hex : ∃ g, IsGammaTop θ w g
  · have hroot := gammaTop_of_exists hex
    have hθw := exists_signal_of_root hroot
    have hbr : ¬ (-1 / gammaTop θ w ∈ Set.Ioo (sStar c w) 0) :=
      fun h => h4 ((assumption4_iff_branch hc).mpr h)
    have hneg : -1 / gammaTop θ w < 0 := neg_inv_gammaTop_neg hex
    have hle : -1 / gammaTop θ w ≤ sStar c w := by
      by_contra h
      exact hbr ⟨not_le.mp h, hneg⟩
    have hs : -1 / gammaTop θ w < sPhys c w z := lt_of_le_of_lt hle (sStar_lt_sPhys hc hw hz)
    have hγ1 : 0 < gammaTop θ w := gammaTop_pos hex
    have hs0 : sPhys c w z < 0 := sPhys_neg hc hw hz
    have hγgt : gammaTop θ w < gamHet c w z := by
      unfold gamHet
      rw [lt_div_iff_of_neg hs0]
      rw [div_lt_iff₀ hγ1] at hs
      linarith
    have hmono := secular_strictMonoOn hθw (Set.mem_Ioi.2 hroot.1) (Set.mem_Ioi.2 hγW) hγgt
    rwa [hroot.2] at hmono
  · exact secular_pos_of_no_root hex hγW

/-- `g_i'(z) ≥ 0` on `(bHet, ∞)`. -/
theorem ghetDeriv_nonneg (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hz : bHet c w < z)
    (i : Fin M) : 0 ≤ ghetDeriv c w i z := by
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hp := one_add_mul_sPhys_pos hc hw hz i
  have hs' := sPhysDeriv_pos hc hw hz
  unfold ghetDeriv
  refine div_nonneg ?_ (sq_nonneg _)
  have : 0 ≤ z * w i ^ 2 * sPhysDeriv c w z := mul_nonneg (mul_nonneg hz0.le (sq_nonneg _)) hs'.le
  linarith

/-- **The blow-up at the edge.** `F'(z) → +∞` as `z ↓ bHet`: the table `i₀` with
`w_{i₀} ≠ 0` contributes at least `bHet c_{i₀} w_{i₀}⁴ s'(z)/(bHet + 1)²` on
`(bHet, bHet + 1)`, every other term is nonnegative, and `s'(z) → +∞`
(`sPhysDeriv_tendsto_atTop`). -/
theorem FhetDeriv_tendsto_atTop (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    Tendsto (FhetDeriv θ c w) (𝓝[>] bHet c w) atTop := by
  obtain ⟨i₀, hi₀⟩ := hw
  have hw : ∃ i, w i ≠ 0 := ⟨i₀, hi₀⟩
  have hb := bHet_pos hc hw
  have hw2 : 0 < w i₀ ^ 2 := by positivity
  have hw4 : 0 < w i₀ ^ 4 := by positivity
  set A := bHet c w * c i₀ * w i₀ ^ 4 / (bHet c w + 1) ^ 2 with hA
  have hApos : 0 < A := by
    rw [hA]
    exact div_pos (mul_pos (mul_pos hb (hc i₀)) hw4) (by positivity)
  refine tendsto_atTop_mono' _ ?_ ((sPhysDeriv_tendsto_atTop hc hw).const_mul_atTop hApos)
  filter_upwards [Ioo_mem_nhdsGT (show bHet c w < bHet c w + 1 by linarith)] with z hz
  have hz1 : bHet c w < z := hz.1
  have hz0 : 0 < z := lt_trans hb hz1
  have hs0 := sPhys_neg hc hw hz1
  have hs' := sPhysDeriv_pos hc hw hz1
  have hp := one_add_mul_sPhys_pos hc hw hz1 i₀
  have hp1 : 1 + w i₀ ^ 2 * sPhys c w z ≤ 1 := by nlinarith
  have hden : 0 < (z * (1 + w i₀ ^ 2 * sPhys c w z)) ^ 2 := pow_pos (mul_pos hz0 hp) 2
  have hden_le : (z * (1 + w i₀ ^ 2 * sPhys c w z)) ^ 2 ≤ (bHet c w + 1) ^ 2 := by
    have h1 : z * (1 + w i₀ ^ 2 * sPhys c w z) ≤ bHet c w + 1 := by
      calc z * (1 + w i₀ ^ 2 * sPhys c w z) ≤ z * 1 := mul_le_mul_of_nonneg_left hp1 hz0.le
        _ ≤ bHet c w + 1 := by linarith [hz.2]
    exact pow_le_pow_left₀ (mul_pos hz0 hp).le h1 2
  have hzs : 0 ≤ z * w i₀ ^ 2 * sPhysDeriv c w z :=
    mul_nonneg (mul_nonneg hz0.le (sq_nonneg _)) hs'.le
  have hnum : bHet c w * w i₀ ^ 2 * sPhysDeriv c w z
      ≤ 1 + w i₀ ^ 2 * sPhys c w z + z * w i₀ ^ 2 * sPhysDeriv c w z := by
    have h1 : bHet c w * w i₀ ^ 2 * sPhysDeriv c w z ≤ z * w i₀ ^ 2 * sPhysDeriv c w z :=
      mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hz1.le hw2.le) hs'.le
    linarith
  have hterm : bHet c w * w i₀ ^ 2 * sPhysDeriv c w z / (bHet c w + 1) ^ 2
      ≤ ghetDeriv c w i₀ z := by
    unfold ghetDeriv
    exact div_le_div₀ (by linarith) hnum hden hden_le
  have hPsi : c i₀ * w i₀ ^ 2 * ghetDeriv c w i₀ z ≤ PsihetDeriv c w z := by
    unfold PsihetDeriv
    exact Finset.single_le_sum (f := fun i => c i * w i ^ 2 * ghetDeriv c w i z)
      (fun i _ => mul_nonneg (mul_nonneg (hc i).le (sq_nonneg _)) (ghetDeriv_nonneg hc hw hz1 i))
      (Finset.mem_univ i₀)
  have hPhi : 0 ≤ PhihetDeriv θ c w z := by
    unfold PhihetDeriv
    exact Finset.sum_nonneg fun i _ => mul_nonneg (by positivity) (ghetDeriv_nonneg hc hw hz1 i)
  have hAeq : A * sPhysDeriv c w z
      = c i₀ * w i₀ ^ 2 * (bHet c w * w i₀ ^ 2 * sPhysDeriv c w z / (bHet c w + 1) ^ 2) := by
    rw [hA]
    field_simp
  have hcw : 0 ≤ c i₀ * w i₀ ^ 2 := mul_nonneg (hc i₀).le (sq_nonneg _)
  have h2 := mul_le_mul_of_nonneg_left hterm hcw
  change A * sPhysDeriv c w z ≤ FhetDeriv θ c w z
  unfold FhetDeriv
  linarith [hAeq, h2, hPsi, hPhi]

end MPhet

/-! ### Deterministic layer: the `n`-side overlap bound and the lower bound on `λ` -/

namespace R6het

variable {p q : ℕ}

/-- For `z ≥ λ_max(W₀) + 1`, `-Φ_y(z) = ∑ c_a²/(z - μ_a) ≤ ∑ c_a² = y ⬝ᵥ y`. -/
theorem neg_qform_le_dotProduct_self {W₀ : Matrix (Fin p) (Fin p) ℝ} (hW₀ : W₀.IsHermitian)
    {z : ℝ} (hz : lamMax W₀ hW₀ + 1 ≤ z) (y : Fin p → ℝ) : -R4.qform W₀ z y ≤ y ⬝ᵥ y := by
  have hz' : lamMax W₀ hW₀ < z := by linarith
  rw [R4.qform_eq_sum hW₀ hz', ← R4.dotProduct_transpose_eigU hW₀ y, ← Finset.sum_neg_distrib]
  simp only [dotProduct]
  refine Finset.sum_le_sum fun a _ => ?_
  have hev : hW₀.eigenvalues a ≤ lamMax W₀ hW₀ := R4.eigenvalues_le_lamMax hW₀ a
  have h1 : 1 ≤ z - hW₀.eigenvalues a := by linarith
  have h2 : (z - hW₀.eigenvalues a)⁻¹ ≤ 1 := inv_le_one_of_one_le₀ h1
  have h3 : (z - hW₀.eigenvalues a)⁻¹ = -(hW₀.eigenvalues a - z)⁻¹ := by
    rw [← inv_neg, neg_sub]
  have hc2 : 0 ≤ ((R4.eigU hW₀)ᵀ *ᵥ y) a ^ 2 := sq_nonneg _
  have h4 : (z - hW₀.eigenvalues a)⁻¹ * ((R4.eigU hW₀)ᵀ *ᵥ y) a ^ 2
      ≤ ((R4.eigU hW₀)ᵀ *ᵥ y) a ^ 2 := mul_le_of_le_one_left hc2 h2
  rw [h3, sq] at h4
  linarith

/-- Rayleigh at a unit vector: `‖X v‖² = v ⬝ᵥ (Xᵀ X v) ≤ λ_max(Xᵀ X)`. -/
theorem mulVec_dotProduct_self_le_gramLamMax (X : Matrix (Fin p) (Fin q) ℝ) {v : Fin q → ℝ}
    (hv : v ⬝ᵥ v = 1) : (X *ᵥ v) ⬝ᵥ (X *ᵥ v) ≤ gramLamMax X := by
  have h := R4.dotProduct_mulVec_le_lamMax (isHermitian_transpose_mul_self X) v
  rw [hv, mul_one, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
    Matrix.vecMul_transpose] at h
  exact h

/-- Simplicity moves from `Xᵀ X` to `X Xᵀ` when the shared top eigenvalue is positive:
`Het.topSimple_transpose_mul_self` applied to `Xᵀ`. -/
theorem topSimple_mul_transpose_self (hp : 0 < p) (hq : 0 < q) (X : Matrix (Fin p) (Fin q) ℝ)
    (hlam : 0 < gramLamMax X)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X)) :
    TopSimple (X * Xᵀ) (isHermitian_mul_transpose_self X) := by
  have hTT : Xᵀᵀ = X := Matrix.transpose_transpose X
  have hgl : gramLamMax Xᵀ = gramLamMax X := by
    rw [Het.lamMax_gram_comm hq hp Xᵀ]
    exact lamMax_congr (by rw [hTT]) _ _
  have hlam' : 0 < gramLamMax Xᵀ := by rw [hgl]; exact hlam
  have hsimple' : TopSimple (Xᵀ * Xᵀᵀ) (isHermitian_mul_transpose_self Xᵀ) :=
    topSimple_congr (by rw [hTT]) _ _ hsimple
  have h := Het.topSimple_transpose_mul_self hq hp Xᵀ hlam' hsimple'
  exact topSimple_congr (by rw [hTT]) _ _ h

/-- **The `n`-side overlap bound** (R6' step 2 for the column split). With
`X Xᵀ = W + u uᵀ`, `X v = u`, a simple top eigenvalue and `0 < λ ≤ z₀`,
`overlap X v ≤ 1 / (λ · uᵀ G(z₀)² u)`. Proof: `Het.topProj_transpose_eq` writes
`overlap X v = (û ⬝ᵥ u)² / λ` for a unit top eigenvector `û` of `X Xᵀ`, `overlap_eq_inner_sq`
identifies `(û ⬝ᵥ u)²` with `overlap Xᵀ u`, and `R6.overlap_q_le` at `Xᵀ` bounds that by
`1 / qform2 W z₀ u` (including the degenerate case `λ = λ_max(W)`). -/
theorem overlap_le_inv_lam_qform2 (hp : 0 < p) (hq : 0 < q) (X : Matrix (Fin p) (Fin q) ℝ)
    (W : Matrix (Fin p) (Fin p) ℝ) (hW : W.IsHermitian) (u : Fin p → ℝ)
    (v : EuclideanSpace ℝ (Fin q))
    (hgram : X * Xᵀ = W + Matrix.vecMulVec u u) (hXv : X *ᵥ WithLp.ofLp v = u)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    (hlam : 0 < gramLamMax X) {z₀ : ℝ} (hz₀ : gramLamMax X ≤ z₀)
    (hpos : 0 < R4.qform2 W z₀ u) :
    overlap X v ≤ 1 / (gramLamMax X * R4.qform2 W z₀ u) := by
  have hTT : Xᵀᵀ = X := Matrix.transpose_transpose X
  have hA := isHermitian_mul_transpose_self X
  have hsimpleN := topSimple_mul_transpose_self hp hq X hlam hsimple
  -- a unit top eigenvector of `X Xᵀ`
  have hne : topSpace (X * Xᵀ) hA ≠ ⊥ := by
    intro h
    have h1 : Module.finrank ℝ (topSpace (X * Xᵀ) hA) = 1 := hsimpleN
    rw [h, finrank_bot] at h1
    exact zero_ne_one h1
  obtain ⟨x, hx, hx0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hne
  set uh : EuclideanSpace ℝ (Fin p) := (‖x‖⁻¹ : ℝ) • x with huh
  have huhmem : uh ∈ topSpace (X * Xᵀ) hA := Submodule.smul_mem _ _ hx
  have huhnorm : ‖uh‖ = 1 := norm_smul_inv_norm hx0
  -- `overlap X v = (û ⬝ᵥ u)² / λ`
  have h1 := Het.topProj_transpose_eq X hsimpleN hlam uh huhmem huhnorm v
  rw [hXv] at h1
  -- `(û ⬝ᵥ u)² = overlap Xᵀ u`
  have hsimpleT : TopSimple (Xᵀᵀ * Xᵀ) (isHermitian_transpose_mul_self Xᵀ) :=
    topSimple_congr (by rw [hTT]) _ _ hsimpleN
  have huhT : uh ∈ topSpace (Xᵀᵀ * Xᵀ) (isHermitian_transpose_mul_self Xᵀ) := by
    rw [topSpace_congr (by rw [hTT] : Xᵀᵀ * Xᵀ = X * Xᵀ) _ hA]
    exact huhmem
  have h2 : overlap Xᵀ (WithLp.toLp 2 u) = (WithLp.ofLp uh ⬝ᵥ u) ^ 2 := by
    rw [overlap_eq_inner_sq Xᵀ (WithLp.toLp 2 u) hsimpleT huhT huhnorm,
      real_inner_eq_dotProduct]
  -- `R6.overlap_q_le` at `Xᵀ`
  have hgramT : Xᵀᵀ * Xᵀ = W + Matrix.vecMulVec u u := by rw [hTT]; exact hgram
  have hz₀T : gramLamMax Xᵀ ≤ z₀ := by
    have hgl : gramLamMax Xᵀ = gramLamMax X := by
      rw [Het.lamMax_gram_comm hq hp Xᵀ]
      exact lamMax_congr (by rw [hTT]) _ _
    rw [hgl]
    exact hz₀
  have h3 := R6.overlap_q_le Xᵀ hW hgramT hp hsimpleT hz₀T hpos
  rw [h2] at h3
  rw [h1]
  calc (WithLp.ofLp uh ⬝ᵥ u) ^ 2 / gramLamMax X
      ≤ 1 / R4.qform2 W z₀ u / gramLamMax X := div_le_div_of_nonneg_right h3 hlam.le
    _ = 1 / (gramLamMax X * R4.qform2 W z₀ u) := by rw [div_div, mul_comm]

end R6het

/-! ### The model: item H12 -/

namespace MultiTableModel

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-- **R5 step 9 on the `n` side.** If the secular function of `W₀' + q qᵀ` is positive at
`z₀ > λ_max(W₀')` then `λ_max(X_wᵀ X_w) ≤ z₀`. -/
theorem gramLamMax_le_of_secular_pos_het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {z₀ : ℝ} (h₀ : lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) < z₀)
    (hs : 0 < secular (m.W0het w N ω) (m.qHet w N ω) z₀) :
    gramLamMax ((m.stackW w).X N ω) ≤ z₀ := by
  by_contra hcon
  rw [not_le] at hcon
  have hW := m.isHermitian_W0het w N ω
  have hA' : (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω)).IsHermitian := by
    rw [← m.gram_eq_het w N ω]; exact Het.isHermitian_mul_transpose_self _
  have hgl : gramLamMax ((m.stackW w).X N ω)
      = lamMax (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω)) hA' := by
    rw [Het.lamMax_gram_comm (m.stack.hn N) (m.stack.hd N)]
    exact lamMax_congr (m.gram_eq_het w N ω) _ hA'
  have hlt : lamMax (m.W0het w N ω) hW < gramLamMax ((m.stackW w).X N ω) := lt_trans h₀ hcon
  obtain ⟨j, hj⟩ := exists_eigenvalues_eq_lamMax hA' (m.stack.hn N)
  have hspec : lamMax (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω)) hA'
      ∈ spectrum ℝ (toOp (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω))) := by
    rw [spectrum_toOp]
    exact hj ▸ hA'.eigenvalues_mem_spectrum_real j
  have hzero : secular (m.W0het w N ω) (m.qHet w N ω) (gramLamMax ((m.stackW w).X N ω)) = 0 := by
    rw [hgl]
    exact (secular_eq_zero_iff hW (hgl ▸ hlt)).mpr hspec
  rcases eq_or_ne (m.qHet w N ω) 0 with hq0 | hq
  · have h1 : secular (m.W0het w N ω) (m.qHet w N ω) (gramLamMax ((m.stackW w).X N ω)) = 1 := by
      change 1 + qform (m.W0het w N ω) (gramLamMax ((m.stackW w).X N ω)) (m.qHet w N ω) = 1
      rw [hq0]
      simp [qform]
    rw [h1] at hzero
    norm_num at hzero
  · have hmono := secular_strictMonoOn hW hq (Set.mem_Ioi.2 h₀) (Set.mem_Ioi.2 hlt) hcon
    linarith

/-- `q ⬝ᵥ q ≤ λ_max(X_wᵀ X_w)`: Rayleigh at the unit vector `v`, with `X_w v = q`. -/
theorem dotProduct_qHet_le_gramLamMax (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : m.qHet w N ω ⬝ᵥ m.qHet w N ω ≤ gramLamMax ((m.stackW w).X N ω) := by
  rw [← m.stackW_mulVec_v w N ω]
  exact R6het.mulVec_dotProduct_self_le_gramLamMax _ (m.dotProduct_v_self_het N)

/-- **Item H12.** Below `eq:assumption4` the weighted stack has no outlier, so the overlap
tends to `Scalars.Lw θ c w = 0`. The interface is taken at the exact edge `b = bHet`, the
only point where `F'` blows up; `RMT/Het/Sup.lean` supplies it from `HeteroEdge`. The five
bad families are in the header of this file. -/
theorem align_tendstoInProb_het_subcritical [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hG : m.JointGaussianNoise)
    (H : m.ResolventLimitsHet w (MPhet.bHet c w) (MPhet.Phihet (fun i => (m.tbl i).θ) c w)
      (MPhet.Psihet c w) (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w) (MPhet.PsihetDeriv c w))
    (h4 : ¬ Scalars.Assumption4 (fun i => (m.tbl i).θ) c w) :
    TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
      (Scalars.Lw (fun i => (m.tbl i).θ) c w) := by
  set θ : Fin M → ℝ := fun i => (m.tbl i).θ with hθ
  have hLw : Scalars.Lw θ c w = 0 := by rw [Scalars.Lw, if_neg h4]
  rw [hLw]
  intro ε hε
  set b := MPhet.bHet c w with hb
  have hb0 : 0 < b := MPhet.bHet_pos hc hw
  -- the lower bound `κ` on `λ`, from the secular form at `z₁ = b + 2`
  set z₁ := b + 2 with hz₁
  have hz₁b : b < z₁ := by linarith
  have hF1 : MPhet.Fhet θ c w z₁ < 0 := MPhet.Fhet_neg hc hw hz₁b
  set κ := -(MPhet.Fhet θ c w z₁) / 2 with hκ
  have hκpos : 0 < κ := by rw [hκ]; linarith
  -- the point `z₀` where `F'` is large
  obtain ⟨z₀, hz₀K, hz₀b⟩ : ∃ z₀, 2 / (κ * ε) < MPhet.FhetDeriv θ c w z₀ ∧ b < z₀ :=
    (((MPhet.FhetDeriv_tendsto_atTop hc hw).eventually (eventually_gt_atTop _)).and
      self_mem_nhdsWithin).exists
  set K := MPhet.FhetDeriv θ c w z₀ with hK
  have hKpos : 0 < K := lt_trans (by positivity) hz₀K
  have hF0 : 0 < 1 + MPhet.Fhet θ c w z₀ :=
    MPhet.one_add_F_pos_of_not_assumption4 hc hw h4 hz₀b
  set δ := min 1 ((z₀ - b) / 2) with hδ
  have hδpos : 0 < δ := lt_min one_pos (by linarith)
  have hδ1 : δ ≤ 1 := min_le_left _ _
  have hδz : b + δ < z₀ := by
    have := min_le_right 1 ((z₀ - b) / 2)
    linarith
  -- the five bad families
  have hB1 : Tendsto (fun N => μ N
      {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ b + δ}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0het_le w N _).nullMeasurableSet) (H.edge δ hδpos)
  have hB2 := tendstoInProb_secular_het H hz₀b ((1 + MPhet.Fhet θ c w z₀) / 2) (half_pos hF0)
  have hB3 := tendstoInProb_qform2_qHet H hz₀b (K / 2) (half_pos hKpos)
  have hB4 := tendstoInProb_secular_het H hz₁b κ hκpos
  have hB5 : Tendsto (fun N => μ N {ω : Ω N | ¬ TopSimple (m.stackGramW w N ω)
      (m.isHermitian_stackGramW w N ω)}) atTop (𝓝 0) := by
    have hz : ∀ N, μ N {ω : Ω N | ¬ TopSimple (m.stackGramW w N ω)
        (m.isHermitian_stackGramW w N ω)} = 0 := fun N =>
      ae_iff.mp (m.heteroLaw_topSimple_of_gaussian w hG hw N)
    simp only [hz]
    exact tendsto_const_nhds
  refine tendsto_measure_zero_of_subset ?_ (tendsto_measure_zero_union hB1
    (tendsto_measure_zero_union hB2 (tendsto_measure_zero_union hB3
      (tendsto_measure_zero_union hB4 hB5))))
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g1, g2, g3, g4, g5⟩ := hcon
  have hW := m.isHermitian_W0het w N ω
  have hlam0 : lamMax (m.W0het w N ω) hW ≤ b + δ := by
    by_contra hx
    exact g1 (Set.mem_compl hx)
  have hsec0 : 0 < secular (m.W0het w N ω) (m.qHet w N ω) z₀ := by
    have h : ¬ ((1 + MPhet.Fhet θ c w z₀) / 2 ≤ |secular (m.W0het w N ω) (m.qHet w N ω) z₀
        - (1 + MPhet.Phihet θ c w z₀ + MPhet.Psihet c w z₀)|) := g2
    have h2 := abs_lt.mp (not_le.mp h)
    have hF : MPhet.Fhet θ c w z₀ = MPhet.Phihet θ c w z₀ + MPhet.Psihet c w z₀ := rfl
    linarith [h2.1]
  have hQ2 : K / 2 < qform2 (m.W0het w N ω) z₀ (m.qHet w N ω) := by
    have h : ¬ (K / 2 ≤ |qform2 (m.W0het w N ω) z₀ (m.qHet w N ω)
        - (MPhet.PhihetDeriv θ c w z₀ + MPhet.PsihetDeriv c w z₀)|) := g3
    have h2 := abs_lt.mp (not_le.mp h)
    have hK' : K = MPhet.PhihetDeriv θ c w z₀ + MPhet.PsihetDeriv c w z₀ := rfl
    linarith [h2.1]
  have hsec1 : secular (m.W0het w N ω) (m.qHet w N ω) z₁ < 1 - κ := by
    have h : ¬ (κ ≤ |secular (m.W0het w N ω) (m.qHet w N ω) z₁
        - (1 + MPhet.Phihet θ c w z₁ + MPhet.Psihet c w z₁)|) := g4
    have h2 := abs_lt.mp (not_le.mp h)
    have hF : MPhet.Fhet θ c w z₁ = MPhet.Phihet θ c w z₁ + MPhet.Psihet c w z₁ := rfl
    rw [hκ]
    linarith [h2.2]
  have hsimple : TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) :=
    not_not.mp g5
  -- `λ ≤ z₀`
  have hlamz₀ : gramLamMax ((m.stackW w).X N ω) ≤ z₀ :=
    m.gramLamMax_le_of_secular_pos_het w N ω (by linarith) hsec0
  -- `κ ≤ q ⬝ᵥ q ≤ λ`
  have hqq : κ ≤ m.qHet w N ω ⬝ᵥ m.qHet w N ω := by
    have h1 := R6het.neg_qform_le_dotProduct_self hW (z := z₁) (by linarith) (m.qHet w N ω)
    have h2 : secular (m.W0het w N ω) (m.qHet w N ω) z₁
        = 1 + qform (m.W0het w N ω) z₁ (m.qHet w N ω) := rfl
    linarith
  have hkaplam : κ ≤ gramLamMax ((m.stackW w).X N ω) :=
    le_trans hqq (m.dotProduct_qHet_le_gramLamMax w N ω)
  have hlampos : 0 < gramLamMax ((m.stackW w).X N ω) := lt_of_lt_of_le hκpos hkaplam
  have hQpos : 0 < qform2 (m.W0het w N ω) z₀ (m.qHet w N ω) := lt_trans (half_pos hKpos) hQ2
  -- the bound
  have hbound : m.stackPerfW w N ω ≤ 1 / (gramLamMax ((m.stackW w).X N ω)
      * qform2 (m.W0het w N ω) z₀ (m.qHet w N ω)) :=
    R6het.overlap_le_inv_lam_qform2 (m.stack.hn N) (m.stack.hd N) _ _ hW _ _
      (m.gram_eq_het w N ω) (m.stackW_mulVec_v w N ω) hsimple hlampos hlamz₀ hQpos
  have hprod : κ * (K / 2) ≤ gramLamMax ((m.stackW w).X N ω)
      * qform2 (m.W0het w N ω) z₀ (m.qHet w N ω) :=
    mul_le_mul hkaplam hQ2.le (half_pos hKpos).le hlampos.le
  have hfin : m.stackPerfW w N ω < ε := by
    have h1 : 1 / (gramLamMax ((m.stackW w).X N ω) * qform2 (m.W0het w N ω) z₀ (m.qHet w N ω))
        ≤ 1 / (κ * (K / 2)) :=
      one_div_le_one_div_of_le (mul_pos hκpos (half_pos hKpos)) hprod
    have h2 : 1 / (κ * (K / 2)) < ε := by
      rw [div_lt_iff₀ (mul_pos hκpos (half_pos hKpos))]
      rw [div_lt_iff₀ (mul_pos hκpos hε)] at hz₀K
      linarith
    linarith
  have hnn : 0 ≤ m.stackPerfW w N ω := overlap_nonneg _ _
  have hω' : ε ≤ |m.stackPerfW w N ω - 0| := hω
  rw [sub_zero, abs_of_nonneg hnn] at hω'
  linarith

end MultiTableModel

end StackedSVD
