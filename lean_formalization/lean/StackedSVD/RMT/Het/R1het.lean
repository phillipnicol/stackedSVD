/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Stein
import StackedSVD.RMT.Het.MPhet
import StackedSVD.RMT.Het.Split
import StackedSVD.RMT.R1

/-!
# Item H6: the heteroscedastic trace law (complex Silverstein equation)

Task H6 of `notes/archive/plan_heterolaw_A.md` (sections 2.2, 3.4, the H6 row and risk 1 of
section 4). Setting of `RMT/Het/Stein.lean`: `B ~ gaussianMatrix p q`, row scale `τ`, scale
`d`, `W_σ = d⁻¹ (diag τ B)(diag τ B)ᵀ`, `G_σ(z) = (W_σ - z)⁻¹`, `s = d⁻¹ tr (W_σ' - z)⁻¹`.
Rows are grouped by a block map `blk : Fin p → Fin M` with `τ_j = w_{blk j}`.

For `Im z > 0` and along `N` with `d_N → ∞`, `q_N / d_N → 1`, `|J_i| / d_N → c_i`:

* **trace law** `s_N(z) → sGlob c w z`, the unique root in the upper half plane of the
  Silverstein equation `z = zfunC c w s = -1/s + ∑ c_i w_i² / (1 + w_i² s)`;
* **block traces** `|J_i|⁻¹ ∑_{j ∈ J_i} G_σ(z)_{jj} → -1 / (z (1 + w_i² sGlob z))`;
* **block traces of `G_σ²`** (the derivative version item R2 consumes), by Cauchy's estimate;
* **boundary values** at `z = x + iη`, `η ↓ 0`, for real `x > bHet`: `sGlob → sPhys c w x`
  and the block limits `→ ghet c w i x`, which is what item T needs.

Route (plan section 3.4). Stein (`HetStein.stein_block`) gives the mean identity
`1 + z E[g_i] = -w_i² z E[s] E[g_i] + ε_i` with `ε_i → 0` (residual plus a covariance that
the variance bounds of H5 control); the trace identity closes the system on the scalar
`E[s]`, so `zfunC (E[s]) → z`. Uniqueness and stability in `ℂ⁺` (`norm_sub_le_of_root`): the
imaginary-part identity `Im zfunC(s) = Im s (1/|s|² - ∑ c_i w_i⁴/|1 + w_i² s|²)` and one
weighted AM-GM step give `|s - σ| ≤ 2 |s| |zfunC s - z| / Im z` for any root `σ ∈ ℂ⁺`.
Existence of the root at every `z ∈ ℂ⁺` comes from the bounded sequence `E[s_N]` itself
(Bolzano-Weierstrass), so no global branch is constructed by hand; the boundary link to the
real branch `sPhys` is the local holomorphic inverse of `zfunC` at `sPhys x`.
Concentration is `measure_abs_ge_le_of_lipschitz` with the Lipschitz constants of H5.

STATUS: see `notes/archive/agent_reports/h6_r1het.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory Finset
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace HetR1

open HetStein MPhet ResolvDeriv R4

variable {M : ℕ}

/-! ### The complex Silverstein map -/

/-- `zfunC c w s = -1/s + ∑ c_i w_i² / (1 + w_i² s)` on `ℂ` (plan section 2.2). -/
noncomputable def zfunC (c w : Fin M → ℝ) (s : ℂ) : ℂ :=
  -s⁻¹ + ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹

/-- `zfunC' s = 1/s² - ∑ c_i w_i⁴ / (1 + w_i² s)²`. -/
noncomputable def zfunDerivC (c w : Fin M → ℝ) (s : ℂ) : ℂ :=
  (s ^ 2)⁻¹ - ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4 * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) ^ 2)⁻¹

/-- The block-`i` limit `g_i = -1 / (z (1 + w_i² s))` as a function of the root `s`. -/
noncomputable def gC (w : Fin M → ℝ) (i : Fin M) (z s : ℂ) : ℂ :=
  -(z * (1 + ((w i : ℝ) : ℂ) ^ 2 * s))⁻¹

/-- `A(s) = ∑ c_i w_i⁴ / |1 + w_i² s|²`, the quantity of the imaginary-part identity. -/
noncomputable def Aabs (c w : Fin M → ℝ) (s : ℂ) : ℝ :=
  ∑ i, c i * w i ^ 4 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖ ^ 2

theorem zfunC_ofReal (c w : Fin M → ℝ) (s : ℝ) : zfunC c w (s : ℂ) = (zfun c w s : ℂ) := by
  simp only [zfunC, zfun, div_eq_mul_inv, neg_mul, one_mul]
  push_cast
  rfl

theorem zfunDerivC_ofReal (c w : Fin M → ℝ) (s : ℝ) :
    zfunDerivC c w (s : ℂ) = (zfunDeriv c w s : ℂ) := by
  simp only [zfunDerivC, zfunDeriv, div_eq_mul_inv, one_mul]
  push_cast
  rfl

theorem gC_ofReal (c w : Fin M → ℝ) (i : Fin M) (x : ℝ) :
    gC w i (x : ℂ) ((sPhys c w x : ℝ) : ℂ) = (ghet c w i x : ℂ) := by
  simp only [gC, ghet, div_eq_mul_inv, neg_mul, one_mul]
  push_cast
  rfl

theorem Aabs_nonneg {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) (s : ℂ) : 0 ≤ Aabs c w s :=
  Finset.sum_nonneg fun i _ => by have := hc i; positivity

/-- **The imaginary-part identity** `Im zfunC(s) = Im s · (1/|s|² - A(s))`, with the junk
conventions `1/0 = 0` consistent on both sides. -/
theorem im_zfunC (c w : Fin M → ℝ) (s : ℂ) :
    (zfunC c w s).im = s.im * (1 / ‖s‖ ^ 2 - Aabs c w s) := by
  have hterm : ∀ i, (((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2
      * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹).im
      = -(s.im * (c i * w i ^ 4 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖ ^ 2)) := by
    intro i
    have h1 : ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 = (((c i * w i ^ 2 : ℝ)) : ℂ) := by
      push_cast; ring
    rw [h1, Complex.im_ofReal_mul, Complex.inv_im, ← Complex.normSq_eq_norm_sq]
    have h2 : (1 + ((w i : ℝ) : ℂ) ^ 2 * s).im = w i ^ 2 * s.im := by
      simp [Complex.add_im, Complex.mul_im, ← Complex.ofReal_pow]
    rw [h2]
    ring
  rw [zfunC, Complex.add_im, Complex.neg_im, Complex.inv_im, ← Complex.normSq_eq_norm_sq,
    Complex.im_sum, Finset.sum_congr rfl fun i _ => hterm i, Aabs, Finset.sum_neg_distrib,
    ← Finset.mul_sum]
  ring

/-- `Aabs s ≤ 1/|s|²` when `Im s > 0` and `Im zfunC(s) ≥ 0`. -/
theorem Aabs_le_of_im_nonneg {c w : Fin M → ℝ} {s : ℂ} (hs : 0 < s.im)
    (hz : 0 ≤ (zfunC c w s).im) : Aabs c w s ≤ 1 / ‖s‖ ^ 2 := by
  rw [im_zfunC] at hz
  have := (mul_nonneg_iff_of_pos_left hs).mp hz
  linarith

/-- `Aabs s ≤ 1/|s|² - Im zfunC(s) / |s|` when `Im s > 0` and `Im zfunC(s) ≥ 0`. -/
theorem Aabs_le_sub {c w : Fin M → ℝ} {s : ℂ} (hs : 0 < s.im) (hz : 0 ≤ (zfunC c w s).im) :
    Aabs c w s ≤ 1 / ‖s‖ ^ 2 - (zfunC c w s).im / ‖s‖ := by
  have hid := im_zfunC c w s
  have hsn : s.im ≤ ‖s‖ := le_trans (le_abs_self _) (Complex.abs_im_le_norm s)
  have hpos : 0 < ‖s‖ := lt_of_lt_of_le hs hsn
  have hA := Aabs_le_of_im_nonneg hs hz
  have h1 : (zfunC c w s).im / ‖s‖ ≤ (zfunC c w s).im / s.im :=
    div_le_div_of_nonneg_left hz hs hsn
  have h2 : (zfunC c w s).im / s.im = 1 / ‖s‖ ^ 2 - Aabs c w s := by
    rw [hid, mul_div_cancel_left₀ _ hs.ne']
  linarith

/-- `Im s > 0` forces `1 + w² s ≠ 0`. -/
theorem one_add_mul_ne_zero_of_im_pos {s : ℂ} (hs : 0 < s.im) (a : ℝ) :
    1 + ((a : ℝ) : ℂ) ^ 2 * s ≠ 0 := by
  intro h
  have him : (1 + ((a : ℝ) : ℂ) ^ 2 * s).im = a ^ 2 * s.im := by
    simp [Complex.add_im, Complex.mul_im, ← Complex.ofReal_pow]
  have hre : (1 + ((a : ℝ) : ℂ) ^ 2 * s).re = 1 + a ^ 2 * s.re := by
    simp [Complex.add_re, Complex.mul_re, ← Complex.ofReal_pow]
  rw [h] at him hre
  simp only [Complex.zero_im, zero_eq_mul, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
    pow_eq_zero_iff, Complex.zero_re] at him hre
  rcases him with ha | hs0
  · rw [ha] at hre; simp at hre
  · linarith

theorem ne_zero_of_im_pos {s : ℂ} (hs : 0 < s.im) : s ≠ 0 := fun h => by
  rw [h] at hs; simp at hs

/-- The difference of two values of `zfunC`, factored. -/
theorem zfunC_sub_zfunC (c w : Fin M → ℝ) {s σ : ℂ} (hs : s ≠ 0) (hσ : σ ≠ 0)
    (hs' : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * s ≠ 0) (hσ' : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * σ ≠ 0) :
    zfunC c w s - zfunC c w σ
      = (s - σ) * ((s * σ)⁻¹ - ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
          * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹) := by
  unfold zfunC
  have hpt : ∀ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹
      - ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ)⁻¹
      = -((s - σ) * (((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
          * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹)) := by
    intro i
    have := hs' i
    have := hσ' i
    field_simp
    ring
  have h0 : -s⁻¹ - -σ⁻¹ = (s - σ) * (s * σ)⁻¹ := by
    field_simp
    ring
  rw [mul_sub, Finset.mul_sum]
  calc -s⁻¹ + ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹
        - (-σ⁻¹ + ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ)⁻¹)
      = (-s⁻¹ - -σ⁻¹) + ∑ i, (((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2
          * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹
          - ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ)⁻¹) := by
        rw [Finset.sum_sub_distrib]; ring
    _ = _ := by
        rw [h0, Finset.sum_congr rfl fun i _ => hpt i, Finset.sum_neg_distrib]
        ring

/-- Weighted AM-GM on the cross sum: `|∑ c_i w_i⁴ / ((1 + w_i² s)(1 + w_i² σ))| ≤
(t A(s) + A(σ)/t) / 2` for every `t > 0`. -/
theorem norm_cross_sum_le {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) (s σ : ℂ) {t : ℝ}
    (ht : 0 < t) :
    ‖∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
        * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹‖
      ≤ (t * Aabs c w s + Aabs c w σ / t) / 2 := by
  refine (norm_sum_le _ _).trans ?_
  rw [Aabs, Aabs, Finset.mul_sum, Finset.sum_div, ← Finset.sum_add_distrib, Finset.sum_div]
  refine Finset.sum_le_sum fun i _ => ?_
  set X : ℝ := ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖⁻¹ with hX
  set Y : ℝ := ‖1 + ((w i : ℝ) : ℂ) ^ 2 * σ‖⁻¹ with hY
  have hcw : 0 ≤ c i * w i ^ 4 := by have := hc i; positivity
  have hnorm : ‖((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
      * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹‖
      = c i * w i ^ 4 * (X * Y) := by
    simp only [norm_mul, norm_inv, norm_pow, Complex.norm_real, Real.norm_eq_abs, mul_inv, hX, hY]
    rw [abs_of_nonneg (hc i), show |w i| ^ 4 = (|w i| ^ 2) ^ 2 by ring, sq_abs]
    ring
  rw [hnorm]
  have hXY : X * Y ≤ (t * X ^ 2 + Y ^ 2 / t) / 2 := by
    rw [le_div_iff₀ (by norm_num : (0 : ℝ) < 2), div_eq_mul_inv]
    have := sq_nonneg (t * X - Y)
    have h2 : (t * X - Y) ^ 2 * t⁻¹ ≥ 0 := by positivity
    have h3 : (t * X - Y) ^ 2 * t⁻¹ = t * X ^ 2 - 2 * X * Y + Y ^ 2 * t⁻¹ := by
      field_simp
      ring
    linarith
  have hX2 : X ^ 2 = 1 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖ ^ 2 := by rw [hX, inv_pow, one_div]
  have hY2 : Y ^ 2 = 1 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * σ‖ ^ 2 := by rw [hY, inv_pow, one_div]
  calc c i * w i ^ 4 * (X * Y) ≤ c i * w i ^ 4 * ((t * X ^ 2 + Y ^ 2 / t) / 2) :=
        mul_le_mul_of_nonneg_left hXY hcw
    _ = (t * (c i * w i ^ 4 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖ ^ 2)
          + c i * w i ^ 4 / ‖1 + ((w i : ℝ) : ℂ) ^ 2 * σ‖ ^ 2 / t) / 2 := by
        rw [hX2, hY2]; ring

/-- **Stability of the Silverstein root** (plan section 3.4). If `σ ∈ ℂ⁺` is a root at `z`
and `s` has `Im s ≥ 0`, `Im zfunC(s) > 0`, then `|s - σ| ≤ 2 |s| |zfunC(s) - z| / Im z`.
With `zfunC s = z` it is uniqueness of the root in `ℂ⁺`. -/
theorem norm_sub_le_of_root {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) {z s σ : ℂ}
    (hz : 0 < z.im) (hs : 0 ≤ s.im) (hzs : 0 < (zfunC c w s).im) (hσ : 0 < σ.im)
    (hσz : zfunC c w σ = z) :
    ‖s - σ‖ ≤ 2 * ‖s‖ * ‖zfunC c w s - z‖ / z.im := by
  have hs' : 0 < s.im := by
    rcases hs.lt_or_eq with h | h
    · exact h
    · exfalso
      have := im_zfunC c w s
      rw [← h, zero_mul] at this
      linarith
  have hs0 := ne_zero_of_im_pos hs'
  have hσ0 := ne_zero_of_im_pos hσ
  have hsn : 0 < ‖s‖ := norm_pos_iff.mpr hs0
  have hσn : 0 < ‖σ‖ := norm_pos_iff.mpr hσ0
  set D : ℂ := (s * σ)⁻¹ - ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
    * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹ with hD
  have hfac : zfunC c w s - z = (s - σ) * D := by
    rw [← hσz]
    exact zfunC_sub_zfunC c w hs0 hσ0 (fun i => one_add_mul_ne_zero_of_im_pos hs' _)
      (fun i => one_add_mul_ne_zero_of_im_pos hσ _)
  have hA1 := Aabs_le_of_im_nonneg (c := c) (w := w) hs' hzs.le
  have hA2 := Aabs_le_sub (c := c) (w := w) hσ (by rw [hσz]; exact hz.le)
  rw [hσz] at hA2
  have ht : 0 < ‖s‖ / ‖σ‖ := by positivity
  have hcross := norm_cross_sum_le (c := c) (w := w) hc s σ ht
  have hDlow : z.im / (2 * ‖s‖) ≤ ‖D‖ := by
    have h1 : ‖(s * σ)⁻¹‖ = 1 / (‖s‖ * ‖σ‖) := by rw [norm_inv, norm_mul, one_div]
    have h2 : ‖D‖ ≥ ‖(s * σ)⁻¹‖ - ‖∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
        * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) * (1 + ((w i : ℝ) : ℂ) ^ 2 * σ))⁻¹‖ :=
      norm_sub_norm_le _ _
    have h3 : (‖s‖ / ‖σ‖ * Aabs c w s + Aabs c w σ / (‖s‖ / ‖σ‖)) / 2
        ≤ 1 / (‖s‖ * ‖σ‖) - z.im / (2 * ‖s‖) := by
      have e1 : ‖s‖ / ‖σ‖ * Aabs c w s ≤ ‖s‖ / ‖σ‖ * (1 / ‖s‖ ^ 2) :=
        mul_le_mul_of_nonneg_left hA1 ht.le
      have e2 : Aabs c w σ / (‖s‖ / ‖σ‖) ≤ (1 / ‖σ‖ ^ 2 - z.im / ‖σ‖) / (‖s‖ / ‖σ‖) :=
        div_le_div_of_nonneg_right hA2 ht.le
      have e3 : ‖s‖ / ‖σ‖ * (1 / ‖s‖ ^ 2) = 1 / (‖s‖ * ‖σ‖) := by field_simp
      have e4 : (1 / ‖σ‖ ^ 2 - z.im / ‖σ‖) / (‖s‖ / ‖σ‖) = 1 / (‖s‖ * ‖σ‖) - z.im / ‖s‖ := by
        field_simp
      rw [e3] at e1
      rw [e4] at e2
      have : z.im / (2 * ‖s‖) = (z.im / ‖s‖) / 2 := by ring
      rw [this]
      linarith
    rw [h1] at h2
    linarith
  have hnorm : ‖zfunC c w s - z‖ = ‖s - σ‖ * ‖D‖ := by rw [hfac, norm_mul]
  rw [div_le_iff₀ (by positivity : (0 : ℝ) < 2 * ‖s‖)] at hDlow
  rw [le_div_iff₀ hz, hnorm]
  have := mul_le_mul_of_nonneg_left hDlow (norm_nonneg (s - σ))
  nlinarith [norm_nonneg (s - σ), norm_nonneg D]

/-! ### The closed system of the means (plan section 3.4) -/

/-- The mean system `1 + z E_i = -w_i² z S E_i + ε_i` (Stein) and `z S + κ = ∑ c_i^N (1 + z E_i)`
(trace identity), solved for `S (z - zfunC S)`: every term on the right is small. -/
theorem mul_sub_zfunC_eq (c w : Fin M → ℝ) {z S : ℂ} (hS0 : S ≠ 0)
    (hS : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * S ≠ 0) (cN : Fin M → ℝ) (κ : ℂ) (E ε : Fin M → ℂ)
    (ha : ∀ i, 1 + z * E i = -(((w i : ℝ) : ℂ) ^ 2 * z) * (S * E i) + ε i)
    (hb : z * S + κ = ∑ i, ((cN i : ℝ) : ℂ) * (1 + z * E i)) :
    S * (z - zfunC c w S)
      = (1 - κ) + ∑ i, ((cN i : ℝ) : ℂ) * ε i * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹
        + S * ∑ i, (((cN i : ℝ) : ℂ) - ((c i : ℝ) : ℂ)) * ((w i : ℝ) : ℂ) ^ 2
            * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹ := by
  have hE : ∀ i, 1 + z * E i
      = (((w i : ℝ) : ℂ) ^ 2 * S + ε i) * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹ := by
    intro i
    rw [eq_mul_inv_iff_mul_eq₀ (hS i)]
    linear_combination ha i
  have hzS : z * S = -κ + ∑ i, ((cN i : ℝ) : ℂ)
      * ((((w i : ℝ) : ℂ) ^ 2 * S + ε i) * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹) := by
    rw [Finset.sum_congr rfl fun i _ => congrArg (((cN i : ℝ) : ℂ) * ·) (hE i)] at hb
    linear_combination hb
  have hkey : ∑ i, ((cN i : ℝ) : ℂ)
        * ((((w i : ℝ) : ℂ) ^ 2 * S + ε i) * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹)
      - S * ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹
      = ∑ i, ((cN i : ℝ) : ℂ) * ε i * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹
        + S * ∑ i, (((cN i : ℝ) : ℂ) - ((c i : ℝ) : ℂ)) * ((w i : ℝ) : ℂ) ^ 2
            * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹ := by
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  have hSS : S * S⁻¹ = 1 := mul_inv_cancel₀ hS0
  unfold zfunC
  linear_combination hzS + hkey + hSS

/-- The lower bound on `|1 + w² S|` from the Stein identity and `|E| ≤ 1/η`. -/
theorem norm_one_add_mul_ge {z S E ε : ℂ} {a η : ℝ} (hz : z ≠ 0) (hη : 0 < η)
    (ha : 1 + z * E = -(((a : ℝ) : ℂ) ^ 2 * z) * (S * E) + ε) (hE : ‖E‖ ≤ 1 / η) :
    (1 - ‖ε‖) * η / ‖z‖ ≤ ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖ := by
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz
  have hid : z * E * (1 + ((a : ℝ) : ℂ) ^ 2 * S) = ε - 1 := by linear_combination ha
  have h1 : 1 - ‖ε‖ ≤ ‖z‖ * ‖E‖ * ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖ := by
    rw [← norm_mul, ← norm_mul, hid]
    have h := norm_sub_norm_le (1 : ℂ) ε
    rw [norm_one, norm_sub_rev] at h
    exact h
  have h2 : η * ‖E‖ ≤ 1 := by
    have := mul_le_mul_of_nonneg_left hE hη.le
    rwa [mul_one_div, div_self hη.ne'] at this
  rw [div_le_iff₀ hzn]
  have h3 : (1 - ‖ε‖) * η ≤ η * (‖z‖ * ‖E‖ * ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖) := by
    nlinarith
  have h4 : η * (‖z‖ * ‖E‖ * ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖)
      ≤ ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖ * ‖z‖ := by
    have : η * (‖z‖ * ‖E‖ * ‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖)
        = (η * ‖E‖) * (‖1 + ((a : ℝ) : ℂ) ^ 2 * S‖ * ‖z‖) := by ring
    rw [this]
    exact mul_le_of_le_one_left (by positivity) h2
  linarith

/-- The lower bound on `|S|` from the two mean identities. -/
theorem norm_S_ge {z S : ℂ} {κ : ℂ} {cN w : Fin M → ℝ} (hcN : ∀ i, 0 ≤ cN i) (E ε : Fin M → ℂ)
    (ha : ∀ i, 1 + z * E i = -(((w i : ℝ) : ℂ) ^ 2 * z) * (S * E i) + ε i)
    (hb : z * S + κ = ∑ i, ((cN i : ℝ) : ℂ) * (1 + z * E i)) :
    κ.re - ∑ i, cN i * ‖ε i‖ ≤ ‖z‖ * ‖S‖ * (1 + ∑ i, cN i * w i ^ 2 * ‖E i‖) := by
  have hκ : κ = -(z * S) * (1 + ∑ i, ((cN i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * E i)
      + ∑ i, ((cN i : ℝ) : ℂ) * ε i := by
    rw [Finset.sum_congr rfl fun i _ => congrArg (((cN i : ℝ) : ℂ) * ·) (ha i)] at hb
    have h1 : ∑ i, ((cN i : ℝ) : ℂ) * (-(((w i : ℝ) : ℂ) ^ 2 * z) * (S * E i) + ε i)
        = -(z * S) * ∑ i, ((cN i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * E i
          + ∑ i, ((cN i : ℝ) : ℂ) * ε i := by
      rw [Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun i _ => by ring
    rw [h1] at hb
    linear_combination hb
  have h1 : κ.re ≤ ‖κ‖ := Complex.re_le_norm κ
  have h2 : ‖κ‖ ≤ ‖z‖ * ‖S‖ * (1 + ∑ i, cN i * w i ^ 2 * ‖E i‖) + ∑ i, cN i * ‖ε i‖ := by
    rw [hκ]
    refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
    · rw [norm_mul, norm_neg, norm_mul]
      refine mul_le_mul_of_nonneg_left ?_ (by positivity)
      refine (norm_add_le _ _).trans ?_
      rw [norm_one]
      refine add_le_add le_rfl ((norm_sum_le _ _).trans (Finset.sum_le_sum fun i _ => ?_))
      rw [norm_mul, norm_mul, norm_pow, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
        Real.norm_eq_abs, abs_of_nonneg (hcN i), sq_abs]
    · refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i _ => ?_)
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (hcN i)]
  linarith

/-! ### Holomorphy of `zfunC`, the local inverse at the real branch, the global root -/

theorem hasStrictDerivAt_zfunC (c w : Fin M → ℝ) {s : ℂ} (hs : s ≠ 0)
    (hs' : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * s ≠ 0) :
    HasStrictDerivAt (zfunC c w) (zfunDerivC c w s) s := by
  have h1 : HasStrictDerivAt (fun s : ℂ => -s⁻¹) ((s ^ 2)⁻¹) s := by
    have := (hasStrictDerivAt_inv hs).neg
    rw [neg_neg] at this
    exact this
  have h2 : ∀ i, HasStrictDerivAt
      (fun s : ℂ => ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2 * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)⁻¹)
      (-(((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4 * ((1 + ((w i : ℝ) : ℂ) ^ 2 * s) ^ 2)⁻¹)) s := by
    intro i
    have ha : HasStrictDerivAt (fun s : ℂ => 1 + ((w i : ℝ) : ℂ) ^ 2 * s) (((w i : ℝ) : ℂ) ^ 2)
        s := by
      have := ((hasStrictDerivAt_id s).const_mul (((w i : ℝ) : ℂ) ^ 2)).const_add 1
      simpa using this
    have hinv := (hasStrictDerivAt_inv (hs' i)).comp s ha
    have := hinv.const_mul (((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 2)
    refine this.congr_deriv ?_
    have := hs' i
    field_simp
  have hsum := HasStrictDerivAt.fun_sum fun i (_ : i ∈ Finset.univ) => h2 i
  have := h1.add hsum
  refine (this.congr_deriv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun s => ?_)
  · rw [zfunDerivC, Finset.sum_neg_distrib, sub_eq_add_neg]
  · rfl

theorem continuousAt_zfunC (c w : Fin M → ℝ) {s : ℂ} (hs : s ≠ 0)
    (hs' : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * s ≠ 0) : ContinuousAt (zfunC c w) s :=
  (hasStrictDerivAt_zfunC c w hs hs').hasDerivAt.continuousAt

/-- `zfunC'(σ) ≠ 0` at a root `σ ∈ ℂ⁺`: `|∑ c_i w_i⁴/(1 + w_i² σ)²| ≤ A(σ) < 1/|σ|²`. -/
theorem zfunDerivC_ne_zero_of_root {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) {z σ : ℂ}
    (hz : 0 < z.im) (hσ : 0 < σ.im) (hσz : zfunC c w σ = z) : zfunDerivC c w σ ≠ 0 := by
  have hσ0 := ne_zero_of_im_pos hσ
  have hσn : 0 < ‖σ‖ := norm_pos_iff.mpr hσ0
  have hA := Aabs_le_sub (c := c) (w := w) hσ (by rw [hσz]; exact hz.le)
  rw [hσz] at hA
  have hsum : ‖∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
      * ((1 + ((w i : ℝ) : ℂ) ^ 2 * σ) ^ 2)⁻¹‖ ≤ Aabs c w σ := by
    refine (norm_sum_le _ _).trans (le_of_eq ?_)
    rw [Aabs]
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [norm_mul, norm_inv, norm_pow, Complex.norm_real, Real.norm_eq_abs]
    rw [abs_of_nonneg (hc i), show |w i| ^ 4 = (|w i| ^ 2) ^ 2 by ring, sq_abs, div_eq_mul_inv]
    ring
  intro h0
  have h1 : ‖(σ ^ 2)⁻¹‖ = 1 / ‖σ‖ ^ 2 := by rw [norm_inv, norm_pow, one_div]
  have h2 : (σ ^ 2)⁻¹ = ∑ i, ((c i : ℝ) : ℂ) * ((w i : ℝ) : ℂ) ^ 4
      * ((1 + ((w i : ℝ) : ℂ) ^ 2 * σ) ^ 2)⁻¹ := by
    rw [zfunDerivC, sub_eq_zero] at h0
    exact h0
  have h3 : z.im / ‖σ‖ > 0 := by positivity
  rw [← h1, h2] at hA
  linarith

theorem sPhys_ne_zero_C {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) : ((sPhys c w x : ℝ) : ℂ) ≠ 0 :=
  Complex.ofReal_ne_zero.mpr (sPhys_neg hc hw hx).ne

theorem one_add_mul_sPhys_ne_zero_C {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {x : ℝ} (hx : bHet c w < x) (i : Fin M) :
    1 + ((w i : ℝ) : ℂ) ^ 2 * ((sPhys c w x : ℝ) : ℂ) ≠ 0 := by
  have h := one_add_mul_sPhys_pos hc hw hx i
  have : (1 + ((w i : ℝ) : ℂ) ^ 2 * ((sPhys c w x : ℝ) : ℂ))
      = ((1 + w i ^ 2 * sPhys c w x : ℝ) : ℂ) := by push_cast; ring
  rw [this]
  exact Complex.ofReal_ne_zero.mpr h.ne'

theorem hasStrictDerivAt_zfunC_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {x : ℝ} (hx : bHet c w < x) :
    HasStrictDerivAt (zfunC c w) ((zfunDeriv c w (sPhys c w x) : ℝ) : ℂ)
      ((sPhys c w x : ℝ) : ℂ) := by
  rw [← zfunDerivC_ofReal]
  exact hasStrictDerivAt_zfunC c w (sPhys_ne_zero_C hc hw hx)
    (one_add_mul_sPhys_ne_zero_C hc hw hx)

theorem zfunDeriv_sPhys_ne_zero_C {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {x : ℝ} (hx : bHet c w < x) : ((zfunDeriv c w (sPhys c w x) : ℝ) : ℂ) ≠ 0 :=
  Complex.ofReal_ne_zero.mpr
    (zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hx) (sPhys_neg hc hw hx)).ne'

open Classical in
/-- The local holomorphic inverse of `zfunC` at the real point `sPhys c w x` (plan section
3.4): `sC c w x z` is a root of `zfunC = z` for `z` near `x`, with `sC c w x x = sPhys c w x`.
Junk `0` outside the hypotheses. -/
noncomputable def sC (c w : Fin M → ℝ) (x : ℝ) : ℂ → ℂ :=
  if h : (∀ i, 0 < c i) ∧ (∃ i, w i ≠ 0) ∧ bHet c w < x then
    HasStrictDerivAt.localInverse (zfunC c w) ((zfunDeriv c w (sPhys c w x) : ℝ) : ℂ)
      ((sPhys c w x : ℝ) : ℂ) (hasStrictDerivAt_zfunC_sPhys h.1 h.2.1 h.2.2)
      (zfunDeriv_sPhys_ne_zero_C h.1 h.2.1 h.2.2)
  else 0

section LocalInverse

variable {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ} (hx : bHet c w < x)
include hc hw hx

theorem sC_eq : sC c w x = HasStrictDerivAt.localInverse (zfunC c w)
      ((zfunDeriv c w (sPhys c w x) : ℝ) : ℂ) ((sPhys c w x : ℝ) : ℂ)
      (hasStrictDerivAt_zfunC_sPhys hc hw hx) (zfunDeriv_sPhys_ne_zero_C hc hw hx) := by
  rw [sC, dif_pos ⟨hc, hw, hx⟩]

theorem zfunC_sPhys : zfunC c w ((sPhys c w x : ℝ) : ℂ) = (x : ℂ) := by
  rw [zfunC_ofReal, zfun_sPhys hc hw hx]

theorem sC_ofReal : sC c w x (x : ℂ) = ((sPhys c w x : ℝ) : ℂ) := by
  rw [sC_eq hc hw hx, ← zfunC_sPhys hc hw hx]
  exact HasStrictFDerivAt.localInverse_apply_image _

theorem hasStrictDerivAt_sC :
    HasStrictDerivAt (sC c w x) (((zfunDeriv c w (sPhys c w x) : ℝ) : ℂ)⁻¹) (x : ℂ) := by
  rw [sC_eq hc hw hx, ← zfunC_sPhys hc hw hx]
  exact HasStrictDerivAt.to_localInverse _ _

theorem eventually_zfunC_sC : ∀ᶠ z in 𝓝 (x : ℂ), zfunC c w (sC c w x z) = z := by
  rw [sC_eq hc hw hx]
  have := HasStrictDerivAt.eventually_right_inverse (hasStrictDerivAt_zfunC_sPhys hc hw hx)
    (zfunDeriv_sPhys_ne_zero_C hc hw hx)
  rwa [zfunC_sPhys hc hw hx] at this

/-- Along the ray `x + iη`, `sC → sPhys x`. -/
theorem tendsto_sC : Tendsto (fun η : ℝ => sC c w x ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
    (𝓝 ((sPhys c w x : ℝ) : ℂ)) := by
  have hcont := (hasStrictDerivAt_sC hc hw hx).hasDerivAt.continuousAt
  rw [ContinuousAt, sC_ofReal hc hw hx] at hcont
  exact hcont.comp (MP.tendsto_ofReal_add_mul_I x)

/-- Along the ray, `Im sC(x + iη) > 0` for small `η > 0`: the derivative `1/zfun'(sPhys x)`
is a positive real. -/
theorem eventually_im_sC_pos :
    ∀ᶠ η : ℝ in 𝓝[>] 0, 0 < (sC c w x ((x : ℂ) + (η : ℂ) * Complex.I)).im := by
  set f' : ℝ := zfunDeriv c w (sPhys c w x) with hf'
  have hf'pos : 0 < f' := zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hx) (sPhys_neg hc hw hx)
  have hslope := hasDerivAt_iff_tendsto_slope_zero.mp (hasStrictDerivAt_sC hc hw hx).hasDerivAt
  have hpath : Tendsto (fun η : ℝ => (η : ℂ) * Complex.I) (𝓝[>] 0) (𝓝[≠] (0 : ℂ)) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
    · have : Tendsto (fun η : ℝ => (η : ℂ) * Complex.I) (𝓝 0) (𝓝 ((0 : ℝ) * Complex.I)) :=
        (Complex.continuous_ofReal.tendsto 0).mul_const Complex.I
      simpa using this.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with η hη
      have : (η : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (ne_of_gt hη)
      exact mul_ne_zero this Complex.I_ne_zero
  have hcomp := (hslope.comp hpath)
  have hre := (Complex.continuous_re.tendsto _).comp hcomp
  have hval : ∀ η : ℝ, 0 < η →
      ((((η : ℂ) * Complex.I)⁻¹ • (sC c w x ((x : ℂ) + (η : ℂ) * Complex.I) - sC c w x (x : ℂ))).re)
        = (sC c w x ((x : ℂ) + (η : ℂ) * Complex.I)).im / η := by
    intro η hη
    rw [sC_ofReal hc hw hx, smul_eq_mul]
    have hη' : (η : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hη.ne'
    have hinv : ((η : ℂ) * Complex.I)⁻¹ = ((η⁻¹ : ℝ) : ℂ) * (-Complex.I) := by
      rw [mul_inv, Complex.inv_I, Complex.ofReal_inv]
    rw [hinv, mul_assoc, Complex.re_ofReal_mul, Complex.mul_re]
    simp only [Complex.neg_re, Complex.neg_im, Complex.I_re, Complex.I_im, Complex.sub_im,
      Complex.ofReal_im]
    ring
  have hlim : Tendsto (fun η : ℝ => (sC c w x ((x : ℂ) + (η : ℂ) * Complex.I)).im / η)
      (𝓝[>] 0) (𝓝 ((((f' : ℝ) : ℂ)⁻¹).re)) := by
    refine hre.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with η hη
    exact hval η hη
  have hlimval : ((((f' : ℝ) : ℂ)⁻¹).re) = f'⁻¹ := by
    rw [← Complex.ofReal_inv, Complex.ofReal_re]
  rw [hlimval] at hlim
  have hev := hlim.eventually (lt_mem_nhds (inv_pos.mpr hf'pos))
  filter_upwards [hev, self_mem_nhdsWithin] with η h1 h2
  exact (div_pos_iff_of_pos_right h2).mp (lt_of_le_of_lt (by positivity) h1)

end LocalInverse

/-- `σ` is a root of the Silverstein equation at `z` in the upper half plane. -/
def IsRoot (c w : Fin M → ℝ) (z σ : ℂ) : Prop := 0 < σ.im ∧ zfunC c w σ = z

/-- Uniqueness of the root in `ℂ⁺` (`norm_sub_le_of_root` at zero residual). -/
theorem IsRoot.unique {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) {z σ σ' : ℂ} (hz : 0 < z.im)
    (h : IsRoot c w z σ) (h' : IsRoot c w z σ') : σ = σ' := by
  have := norm_sub_le_of_root hc hz h.1.le (by rw [h.2]; exact hz) h'.1 h'.2
  rw [h.2, sub_self, norm_zero, mul_zero, zero_div] at this
  exact sub_eq_zero.mp (norm_le_zero_iff.mp this)

open Classical in
/-- **The global root** `sGlob c w z`: the unique root of `zfunC = z` in `ℂ⁺` when one exists
(junk `0` otherwise). Existence at every `z ∈ ℂ⁺` is `exists_isRoot` below, from the random
matrix sequence; near the real axis it is the local inverse `sC`. -/
noncomputable def sGlob (c w : Fin M → ℝ) (z : ℂ) : ℂ :=
  if h : ∃ σ, IsRoot c w z σ then h.choose else 0

theorem isRoot_sGlob {c w : Fin M → ℝ} {z : ℂ} (h : ∃ σ, IsRoot c w z σ) :
    IsRoot c w z (sGlob c w z) := by
  rw [sGlob, dif_pos h]
  exact h.choose_spec

theorem sGlob_eq {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) {z σ : ℂ} (hz : 0 < z.im)
    (h : IsRoot c w z σ) : sGlob c w z = σ :=
  IsRoot.unique hc hz (isRoot_sGlob ⟨σ, h⟩) h

/-- **Boundary values.** For real `x > bHet`, `sGlob (x + iη) → sPhys c w x` as `η ↓ 0`
(what item T consumes as `hL`). -/
theorem tendsto_sGlob {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) :
    Tendsto (fun η : ℝ => sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0)
      (𝓝 ((sPhys c w x : ℝ) : ℂ)) := by
  refine (tendsto_sC hc hw hx).congr' ?_
  have h1 := (MP.tendsto_ofReal_add_mul_I x).eventually (eventually_zfunC_sC hc hw hx)
  filter_upwards [h1, eventually_im_sC_pos hc hw hx, self_mem_nhdsWithin] with η hη1 hη2 hη3
  refine (sGlob_eq (fun i => (hc i).le) ?_ ⟨hη2, hη1⟩).symm
  simpa using hη3

/-- Along the ray, `sGlob` is a root (so `zfunC (sGlob) = z`, `Im sGlob > 0`). -/
theorem eventually_isRoot_sGlob {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {x : ℝ} (hx : bHet c w < x) :
    ∀ᶠ η : ℝ in 𝓝[>] 0, IsRoot c w ((x : ℂ) + (η : ℂ) * Complex.I)
      (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I)) := by
  have h1 := (MP.tendsto_ofReal_add_mul_I x).eventually (eventually_zfunC_sC hc hw hx)
  filter_upwards [h1, eventually_im_sC_pos hc hw hx] with η hη1 hη2
  exact isRoot_sGlob ⟨_, hη2, hη1⟩

/-- Boundary values of the block limits: `gC w i z (sGlob z) → ghet c w i x`. -/
theorem tendsto_gC_sGlob {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) (i : Fin M) :
    Tendsto (fun η : ℝ => gC w i ((x : ℂ) + (η : ℂ) * Complex.I)
      (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I))) (𝓝[>] 0) (𝓝 ((ghet c w i x : ℝ) : ℂ)) := by
  rw [← gC_ofReal c w i x]
  have hx0 : (x : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (lt_trans (bHet_pos hc hw) hx).ne'
  have hden : (x : ℂ) * (1 + ((w i : ℝ) : ℂ) ^ 2 * ((sPhys c w x : ℝ) : ℂ)) ≠ 0 :=
    mul_ne_zero hx0 (one_add_mul_sPhys_ne_zero_C hc hw hx i)
  unfold gC
  refine Tendsto.neg (Tendsto.inv₀ ?_ hden)
  refine (MP.tendsto_ofReal_add_mul_I x).mul (tendsto_const_nhds.add (tendsto_const_nhds.mul ?_))
  exact tendsto_sGlob hc hw hx


/-! ### Blocks and the row scale -/

section Blocks

variable {p : ℕ}

/-- The rows of block `i`. -/
def blockSet (blk : Fin p → Fin M) (i : Fin M) : Finset (Fin p) :=
  Finset.univ.filter fun j => blk j = i

/-- The row scale `τ_j = w_{blk j}` (`SigmaHalf` in the model, plan section 2.1). -/
def tauOf (w : Fin M → ℝ) (blk : Fin p → Fin M) : Fin p → ℝ := fun j => w (blk j)

theorem tauOf_sq_eq (w : Fin M → ℝ) (blk : Fin p → Fin M) (i : Fin M) {j : Fin p}
    (hj : j ∈ blockSet blk i) : tauOf w blk j ^ 2 = w i ^ 2 := by
  simp only [blockSet, Finset.mem_filter, Finset.mem_univ, true_and] at hj
  simp [tauOf, hj]

theorem tauOf_sq_le (w : Fin M → ℝ) (blk : Fin p → Fin M) (j : Fin p) :
    tauOf w blk j ^ 2 ≤ Scalars.wSqMax w := by
  unfold tauOf
  exact Scalars.le_wSqMax w (blk j)

theorem sum_blockSet (blk : Fin p → Fin M) (f : Fin p → ℂ) :
    ∑ i, ∑ j ∈ blockSet blk i, f j = ∑ j, f j := by
  classical
  exact Finset.sum_fiberwise Finset.univ blk f

theorem sum_card_blockSet (blk : Fin p → Fin M) : ∑ i, (blockSet blk i).card = p := by
  classical
  have := Finset.card_eq_sum_card_fiberwise (f := blk) (s := Finset.univ) (t := Finset.univ)
    (fun _ _ => Finset.mem_univ _)
  rw [Finset.card_univ, Fintype.card_fin] at this
  exact this.symm

/-- The trace identity in block form (pointwise, plan section 2.2):
`z s + q/d = ∑_i (|J_i|/d) (1 + z g_{J_i})`. -/
theorem mul_ssig_add_eq {z : ℂ} (hz : z.im ≠ 0) {d : ℕ} (hd : 0 < d) (w : Fin M → ℝ)
    (blk : Fin p → Fin M) {q : ℕ} (B : Matrix (Fin p) (Fin q) ℝ) :
    z * ssig (tauOf w blk) d z B + ((q : ℝ) / d : ℝ)
      = ∑ i, ((((blockSet blk i).card : ℝ) / d : ℝ) : ℂ)
          * (1 + z * gsigAvg (tauOf w blk) d z (blockSet blk i) B) := by
  have hz0 : z ≠ 0 := fun h => hz (by rw [h]; simp)
  have hdC : (d : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hd.ne'
  have htr : (Gsig (tauOf w blk) d z B).trace
      = ∑ i, ((blockSet blk i).card : ℂ) * gsigAvg (tauOf w blk) d z (blockSet blk i) B := by
    rw [Matrix.trace, ← sum_blockSet blk]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [gsigAvg]
    by_cases hJ : (blockSet blk i).card = 0
    · rw [Finset.card_eq_zero.mp hJ]
      simp
    · have : ((blockSet blk i).card : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hJ
      rw [← mul_assoc, mul_inv_cancel₀ this, one_mul]
      rfl
  have hp : (p : ℂ) = ∑ i, ((blockSet blk i).card : ℂ) := by
    have := sum_card_blockSet blk
    exact_mod_cast this.symm
  rw [trace_identity hz hd, htr]
  push_cast
  have hsum : ∑ i, ((blockSet blk i).card : ℂ) / d
        * (1 + z * gsigAvg (tauOf w blk) d z (blockSet blk i) B)
      = (p : ℂ) / d + z * ((d : ℂ)⁻¹ * ∑ i, ((blockSet blk i).card : ℂ)
          * gsigAvg (tauOf w blk) d z (blockSet blk i) B) := by
    simp only [mul_add, mul_one, Finset.sum_add_distrib, ← Finset.sum_div, ← hp, Finset.mul_sum]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    ring
  rw [hsum]
  field_simp
  ring

end Blocks

/-! ### The means, their identities and bounds (fixed sizes) -/

section Means

variable {p q d : ℕ} {z : ℂ} (w : Fin M → ℝ) (blk : Fin p → Fin M)

/-- `E[s]`. -/
noncomputable def meanS (q d : ℕ) (z : ℂ) : ℂ :=
  ∫ B, ssig (tauOf w blk) d z B ∂gaussianMatrix p q

/-- `E[g_{J_i}]`. -/
noncomputable def meanG (q d : ℕ) (z : ℂ) (i : Fin M) : ℂ :=
  ∫ B, gsigAvg (tauOf w blk) d z (blockSet blk i) B ∂gaussianMatrix p q

/-- `Cov(s, g_{J_i}) = E[s g_{J_i}] - E[s] E[g_{J_i}]`. -/
noncomputable def covSG (q d : ℕ) (z : ℂ) (i : Fin M) : ℂ :=
  (∫ B, ssig (tauOf w blk) d z B * gsigAvg (tauOf w blk) d z (blockSet blk i) B
    ∂gaussianMatrix p q) - meanS w blk q d z * meanG w blk q d z i

/-- The residual of the mean identity: `ε_i = r_{J_i} - w_i² z Cov(s, g_{J_i})`. -/
noncomputable def epsN (q d : ℕ) (z : ℂ) (i : Fin M) : ℂ :=
  steinErrAvg (tauOf w blk) q d z (blockSet blk i)
    - (((w i : ℝ) : ℂ) ^ 2 * z) * covSG w blk q d z i

variable {w blk}

theorem measurable_gsigAvg (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (J : Finset (Fin p)) :
    Measurable (fun B : Matrix (Fin p) (Fin q) ℝ => gsigAvg τ d z J B) :=
  measurable_const.mul (Finset.measurable_sum _ fun j _ => measurable_gsig hz hp hd τ j)

theorem norm_gsigAvg_le (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) : ‖gsigAvg τ d z J B‖ ≤ 1 / z.im := by
  rw [gsigAvg]
  by_cases hJ : J.card = 0
  · rw [hJ]; simp only [CharP.cast_eq_zero, inv_zero, zero_mul, norm_zero, one_div, inv_nonneg]
    positivity
  have hc : (0 : ℝ) < J.card := by exact_mod_cast Nat.pos_of_ne_zero hJ
  rw [norm_mul, norm_inv, Complex.norm_natCast]
  calc ((J.card : ℝ))⁻¹ * ‖∑ j ∈ J, gsig τ d z B j‖
      ≤ ((J.card : ℝ))⁻¹ * ∑ j ∈ J, (1 / z.im) :=
        mul_le_mul_of_nonneg_left ((norm_sum_le _ _).trans
          (Finset.sum_le_sum fun j _ => norm_gsig_le hz τ d B j)) (by positivity)
    _ = 1 / z.im := by
        rw [Finset.sum_const, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hc.ne', one_mul]

theorem integrable_gsigAvg (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (J : Finset (Fin p)) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => gsigAvg τ d z J B) (gaussianMatrix p q) :=
  Integrable.mono' (integrable_const (1 / z.im))
    (measurable_gsigAvg hz hp hd τ J).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => norm_gsigAvg_le hz τ d J B)

theorem integrable_ssig (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => ssig τ d z B) (gaussianMatrix p q) :=
  Integrable.mono' (integrable_const ((q : ℝ) / ((d : ℝ) * z.im)))
    (measurable_ssig hz hp hd τ).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => norm_ssig_le hz τ d B)

theorem integrable_ssig_mul_gsigAvg (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) (τ : Fin p → ℝ)
    (J : Finset (Fin p)) :
    Integrable (fun B : Matrix (Fin p) (Fin q) ℝ => ssig τ d z B * gsigAvg τ d z J B)
      (gaussianMatrix p q) := by
  refine Integrable.mono' (integrable_const ((q : ℝ) / ((d : ℝ) * z.im) * (1 / z.im)))
    ((measurable_ssig hz hp hd τ).mul (measurable_gsigAvg hz hp hd τ J)).aestronglyMeasurable
    (Filter.Eventually.of_forall fun B => ?_)
  rw [norm_mul]
  exact mul_le_mul (norm_ssig_le hz τ d B) (norm_gsigAvg_le hz τ d J B) (norm_nonneg _)
    (by positivity)

/-- **The Stein identity for the means**: `1 + z E_i = -w_i² z E[s] E_i + ε_i`. -/
theorem mean_stein (hz : 0 < z.im) (hd : 0 < d) (i : Fin M) (hJ : (blockSet blk i).Nonempty) :
    1 + z * meanG w blk (p := p) q d z i
      = -(((w i : ℝ) : ℂ) ^ 2 * z) * (meanS w blk q d z * meanG w blk q d z i)
        + epsN w blk q d z i := by
  have hp : 0 < p := Fin.pos hJ.choose
  have h := stein_block (q := q) hz hd (tauOf w blk) hJ (σ₀ := w i ^ 2)
    (fun j hj => tauOf_sq_eq w blk i hj)
  rw [integral_add (integrable_const _) ((integrable_gsigAvg hz hp hd _ _).const_mul z),
    integral_const, integral_const_mul] at h
  simp only [probReal_univ, one_smul] at h
  push_cast at h
  rw [epsN, covSG, meanG, meanS]
  linear_combination h

/-- **The trace identity for the means**: `z E[s] + q/d = ∑ c_i^N (1 + z E_i)`. -/
theorem mean_trace (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) :
    z * meanS w blk (p := p) q d z + ((q : ℝ) / d : ℝ)
      = ∑ i, ((((blockSet blk i).card : ℝ) / d : ℝ) : ℂ) * (1 + z * meanG w blk q d z i) := by
  have hpt := fun B : Matrix (Fin p) (Fin q) ℝ => mul_ssig_add_eq hz.ne' hd w blk B
  have hint : ∫ B, (z * ssig (tauOf w blk) d z B + ((q : ℝ) / d : ℝ)) ∂gaussianMatrix p q
      = ∫ B, ∑ i, ((((blockSet blk i).card : ℝ) / d : ℝ) : ℂ)
          * (1 + z * gsigAvg (tauOf w blk) d z (blockSet blk i) B) ∂gaussianMatrix p q :=
    integral_congr_ae (Filter.Eventually.of_forall hpt)
  rw [integral_add ((integrable_ssig hz hp hd _).const_mul z) (integrable_const _),
    integral_const_mul, integral_const] at hint
  simp only [probReal_univ, one_smul] at hint
  rw [integral_finsetSum (f := fun i (B : Matrix (Fin p) (Fin q) ℝ) =>
      ((((blockSet blk i).card : ℝ) / d : ℝ) : ℂ)
        * (1 + z * gsigAvg (tauOf w blk) d z (blockSet blk i) B)) Finset.univ
      (fun i _ => ((integrable_const _).add
        ((integrable_gsigAvg hz hp hd _ _).const_mul z)).const_mul _)] at hint
  rw [meanS, hint]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [integral_const_mul, integral_add (integrable_const _)
    ((integrable_gsigAvg hz hp hd _ _).const_mul z), integral_const, integral_const_mul]
  simp only [probReal_univ, one_smul]
  rfl

theorem norm_meanG_le (hz : 0 < z.im) (i : Fin M) : ‖meanG w blk (p := p) q d z i‖ ≤ 1 / z.im := by
  have h := norm_integral_le_of_norm_le_const (μ := gaussianMatrix p q)
    (Filter.Eventually.of_forall fun B => norm_gsigAvg_le hz (tauOf w blk) d (blockSet blk i) B)
  simpa [meanG] using h

theorem norm_meanS_le (hz : 0 < z.im) :
    ‖meanS w blk (p := p) q d z‖ ≤ (q : ℝ) / ((d : ℝ) * z.im) := by
  have h := norm_integral_le_of_norm_le_const (μ := gaussianMatrix p q)
    (Filter.Eventually.of_forall fun B => norm_ssig_le hz (tauOf w blk) d B)
  simpa [meanS] using h

theorem im_ssig_nonneg (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (B : Matrix (Fin p) (Fin q) ℝ) :
    0 ≤ (ssig τ d z B).im := by
  rw [ssig, R4C.trace_resolvC (isHermitian_Wsig' τ d B) hz.ne']
  have hdc : ((d : ℕ) : ℂ)⁻¹ = (((d : ℝ)⁻¹ : ℝ) : ℂ) := by push_cast; ring
  rw [hdc, Complex.im_ofReal_mul, Complex.im_sum]
  refine mul_nonneg (by positivity) (Finset.sum_nonneg fun a _ => ?_)
  rw [Complex.inv_im]
  have : (((isHermitian_Wsig' τ d B).eigenvalues a : ℂ) - z).im = -z.im := by simp
  rw [this, neg_neg]
  exact div_nonneg hz.le (Complex.normSq_nonneg _)

theorem im_meanS_nonneg (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) :
    0 ≤ (meanS w blk (p := p) q d z).im := by
  have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM
    (integrable_ssig (q := q) hz hp hd (tauOf w blk))
  rw [meanS]
  simp only [Complex.imCLM_apply] at h
  rw [← h]
  exact integral_nonneg fun B => im_ssig_nonneg hz _ _ B

/-- The pointwise inequality `‖u‖ ≤ (t ‖u‖² + 1/t) / 2`. -/
theorem norm_le_half_mul_sq_add {u : ℂ} {t : ℝ} (ht : 0 < t) :
    ‖u‖ ≤ (t * ‖u‖ ^ 2 + 1 / t) / 2 := by
  have h := sq_nonneg (t * ‖u‖ - 1)
  have : (t * ‖u‖ - 1) ^ 2 / t = t * ‖u‖ ^ 2 - 2 * ‖u‖ + 1 / t := by
    field_simp
    ring
  have h2 : 0 ≤ (t * ‖u‖ - 1) ^ 2 / t := by positivity
  linarith

/-- **The covariance bound** `|Cov(s, g_{J_i})| ≤ 9 lipS / η`, from `Var(Re s)`, `Var(Im s)`
of H5 and `|g - E g| ≤ 2/η`. -/
theorem norm_covSG_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    {S : ℝ} (hS : ∀ j, tauOf w blk j ^ 2 ≤ S) (hS0 : 0 < S) (i : Fin M) :
    ‖covSG w blk q d z i‖ ≤ 9 * lipS z d p S / z.im := by
  rw [covSG, meanS, meanG]
  set ν := gaussianMatrix p q with hν
  set τ := tauOf w blk with hτ
  set J := blockSet blk i with hJ
  set f : Matrix (Fin p) (Fin q) ℝ → ℂ := fun B => ssig τ d z B with hf
  set g : Matrix (Fin p) (Fin q) ℝ → ℂ := fun B => gsigAvg τ d z J B with hg
  set Ef : ℂ := ∫ B, f B ∂ν with hEf
  set Eg : ℂ := ∫ B, g B ∂ν with hEg
  have hfI : Integrable f ν := integrable_ssig hz hp hd τ
  have hgI : Integrable g ν := integrable_gsigAvg hz hp hd τ J
  have hfgI : Integrable (fun B => f B * g B) ν := integrable_ssig_mul_gsigAvg hz hp hd τ J
  have hL := lipS_pos hz hp hd hS0
  -- the covariance as a centered integral
  have hI1 : Integrable (fun B => f B * g B - Ef * g B) ν := hfgI.sub (hgI.const_mul Ef)
  have hI2 : Integrable (fun B => f B * Eg - Ef * Eg) ν :=
    (hfI.mul_const Eg).sub (integrable_const _)
  have hI3 : Integrable (fun B => (f B - Ef) * (g B - Eg)) ν := by
    refine (hI1.sub hI2).congr (Filter.Eventually.of_forall fun B => ?_)
    simp only [Pi.sub_apply]
    ring
  have hcov : (∫ B, f B * g B ∂ν) - Ef * Eg = ∫ B, (f B - Ef) * (g B - Eg) ∂ν := by
    have e : ∀ B, (f B - Ef) * (g B - Eg) = (f B * g B - Ef * g B) - (f B * Eg - Ef * Eg) := by
      intro B; ring
    rw [integral_congr_ae (Filter.Eventually.of_forall e), integral_sub hI1 hI2,
      integral_sub hfgI (hgI.const_mul Ef), integral_sub (hfI.mul_const Eg) (integrable_const _),
      integral_const_mul, integral_mul_const, integral_const]
    simp only [probReal_univ, one_smul]
    ring
  -- pointwise bound
  have hgb : ∀ B, ‖g B - Eg‖ ≤ 2 / z.im := by
    intro B
    have h1 := norm_gsigAvg_le hz τ d J B
    have h2 : ‖Eg‖ ≤ 1 / z.im := by
      have h := norm_integral_le_of_norm_le_const (μ := ν)
        (Filter.Eventually.of_forall fun B => norm_gsigAvg_le hz τ d J B)
      simpa using h
    calc ‖g B - Eg‖ ≤ ‖g B‖ + ‖Eg‖ := norm_sub_le _ _
      _ ≤ 1 / z.im + 1 / z.im := add_le_add h1 h2
      _ = 2 / z.im := by ring
  have hpt : ∀ B, ‖(f B - Ef) * (g B - Eg)‖ ≤ (2 / z.im) * ‖f B - Ef‖ := by
    intro B
    rw [norm_mul, mul_comm]
    exact mul_le_mul_of_nonneg_right (hgb B) (norm_nonneg _)
  have hfc : Integrable (fun B => ‖f B - Ef‖) ν := (hfI.sub (integrable_const _)).norm
  have h1 : ‖(∫ B, f B * g B ∂ν) - Ef * Eg‖ ≤ (2 / z.im) * ∫ B, ‖f B - Ef‖ ∂ν := by
    rw [hcov]
    refine (norm_integral_le_integral_norm _).trans ?_
    rw [← integral_const_mul]
    exact integral_mono hI3.norm (hfc.const_mul _) hpt
  -- the second moment
  have hre : ∫ B, (f B).re ∂ν = Ef.re := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.reCLM hfI
    simpa using h
  have him : ∫ B, (f B).im ∂ν = Ef.im := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM hfI
    simpa using h
  have hvre := variance_re_ssig_le (q := q) hz hp hq hd τ hS hS0
  have hvim := variance_im_ssig_le (q := q) hz hp hq hd τ hS hS0
  rw [hre] at hvre
  rw [him] at hvim
  have hsq : ∀ B, ‖f B - Ef‖ ^ 2 = ((f B).re - Ef.re) ^ 2 + ((f B).im - Ef.im) ^ 2 := by
    intro B
    rw [Complex.sq_norm, Complex.normSq_apply, Complex.sub_re, Complex.sub_im]
    ring
  have hIre : Integrable (fun B => ((f B).re - Ef.re) ^ 2) ν := by
    refine Integrable.mono' (integrable_const ((2 * ((q : ℝ) / ((d : ℝ) * z.im))) ^ 2))
      (((Complex.measurable_re.comp (measurable_ssig hz hp hd τ)).sub measurable_const).pow_const
        2).aestronglyMeasurable (Filter.Eventually.of_forall fun B => ?_)
    rw [Real.norm_eq_abs, abs_pow]
    refine pow_le_pow_left₀ (abs_nonneg _) ?_ 2
    have h1 : |(f B).re| ≤ (q : ℝ) / ((d : ℝ) * z.im) :=
      (Complex.abs_re_le_norm _).trans (norm_ssig_le hz τ d B)
    have h2 : |Ef.re| ≤ (q : ℝ) / ((d : ℝ) * z.im) := by
      refine (Complex.abs_re_le_norm _).trans ?_
      have h := norm_integral_le_of_norm_le_const (μ := ν)
        (Filter.Eventually.of_forall fun B => norm_ssig_le hz τ d B)
      simpa using h
    calc |(f B).re - Ef.re| ≤ |(f B).re| + |Ef.re| := abs_sub _ _
      _ ≤ _ := by linarith
  have hIim : Integrable (fun B => ((f B).im - Ef.im) ^ 2) ν := by
    refine Integrable.mono' (integrable_const ((2 * ((q : ℝ) / ((d : ℝ) * z.im))) ^ 2))
      (((Complex.measurable_im.comp (measurable_ssig hz hp hd τ)).sub measurable_const).pow_const
        2).aestronglyMeasurable (Filter.Eventually.of_forall fun B => ?_)
    rw [Real.norm_eq_abs, abs_pow]
    refine pow_le_pow_left₀ (abs_nonneg _) ?_ 2
    have h1 : |(f B).im| ≤ (q : ℝ) / ((d : ℝ) * z.im) :=
      (Complex.abs_im_le_norm _).trans (norm_ssig_le hz τ d B)
    have h2 : |Ef.im| ≤ (q : ℝ) / ((d : ℝ) * z.im) := by
      refine (Complex.abs_im_le_norm _).trans ?_
      have h := norm_integral_le_of_norm_le_const (μ := ν)
        (Filter.Eventually.of_forall fun B => norm_ssig_le hz τ d B)
      simpa using h
    calc |(f B).im - Ef.im| ≤ |(f B).im| + |Ef.im| := abs_sub _ _
      _ ≤ _ := by linarith
  have h2m : ∫ B, ‖f B - Ef‖ ^ 2 ∂ν ≤ 8 * lipS z d p S ^ 2 := by
    rw [integral_congr_ae (Filter.Eventually.of_forall hsq), integral_add hIre hIim]
    linarith
  -- the first absolute moment through `t = 1/lipS`
  set t : ℝ := (lipS z d p S)⁻¹ with ht
  have htpos : 0 < t := inv_pos.mpr hL
  have h1m : ∫ B, ‖f B - Ef‖ ∂ν ≤ 9 * lipS z d p S / 2 := by
    have hI2 : Integrable (fun B => ‖f B - Ef‖ ^ 2) ν := by
      refine (hIre.add hIim).congr (Filter.Eventually.of_forall fun B => ?_)
      simp only [Pi.add_apply]
      rw [hsq]
    calc ∫ B, ‖f B - Ef‖ ∂ν ≤ ∫ B, (t * ‖f B - Ef‖ ^ 2 + 1 / t) / 2 ∂ν :=
          integral_mono hfc (((hI2.const_mul t).add (integrable_const _)).div_const 2)
            fun B => norm_le_half_mul_sq_add htpos
      _ = (t * ∫ B, ‖f B - Ef‖ ^ 2 ∂ν + 1 / t) / 2 := by
          rw [integral_div, integral_add (hI2.const_mul t) (integrable_const _),
            integral_const_mul, integral_const]
          simp only [probReal_univ, one_smul]
      _ ≤ (t * (8 * lipS z d p S ^ 2) + 1 / t) / 2 := by
          gcongr
      _ = 9 * lipS z d p S / 2 := by
          rw [ht]
          field_simp
          ring
  calc ‖(∫ B, f B * g B ∂ν) - Ef * Eg‖ ≤ (2 / z.im) * ∫ B, ‖f B - Ef‖ ∂ν := h1
    _ ≤ (2 / z.im) * (9 * lipS z d p S / 2) :=
        mul_le_mul_of_nonneg_left h1m (by positivity)
    _ = 9 * lipS z d p S / z.im := by ring

/-- The residual bound `‖ε_i‖ ≤ w_i²/d (1/η + ‖z‖/η²) + w_i² ‖z‖ · 9 lipS / η`. -/
theorem norm_epsN_le (hz : 0 < z.im) (hp : 0 < p) (hq : 0 < q) (hd : 0 < d)
    {S : ℝ} (hS : ∀ j, tauOf w blk j ^ 2 ≤ S) (hS0 : 0 < S) (i : Fin M)
    (hJ : (blockSet blk i).Nonempty) :
    ‖epsN w blk q d z i‖
      ≤ w i ^ 2 / d * (1 / z.im + ‖z‖ / z.im ^ 2) + w i ^ 2 * ‖z‖ * (9 * lipS z d p S / z.im) := by
  rw [epsN]
  refine (norm_sub_le _ _).trans (add_le_add ?_ ?_)
  · exact norm_steinErrAvg_le hz hd (tauOf w blk) hJ (σ₀ := w i ^ 2)
      (fun j hj => tauOf_sq_eq w blk i hj)
  · rw [norm_mul, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs]
    exact mul_le_mul_of_nonneg_left (norm_covSG_le hz hp hq hd hS hS0 i) (by positivity)

end Means


/-! ### The limit of the means along `N` (plan section 3.4) -/

section Sequence

variable {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

/-- `c_i^N = |J_i| / d_N`. -/
noncomputable def cN (blk : (N : ℕ) → Fin (pN N) → Fin M) (dN : ℕ → ℕ) (N : ℕ) (i : Fin M) : ℝ :=
  ((blockSet (blk N) i).card : ℝ) / dN N

theorem cN_nonneg (blk : (N : ℕ) → Fin (pN N) → Fin M) (dN : ℕ → ℕ) (N : ℕ) (i : Fin M) :
    0 ≤ cN blk dN N i := by unfold cN; positivity

theorem sum_cN (blk : (N : ℕ) → Fin (pN N) → Fin M) (dN : ℕ → ℕ) (N : ℕ) :
    ∑ i, cN blk dN N i = (pN N : ℝ) / dN N := by
  simp only [cN, ← Finset.sum_div]
  congr 1
  exact_mod_cast sum_card_blockSet (blk N)

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))

include hc hw hd hq hn in
/-- Eventually every size is positive and every block is nonempty. -/
theorem eventually_sizes_pos : ∀ᶠ N in atTop, 0 < dN N ∧ 0 < qN N ∧
    (∀ i, (blockSet (blk N) i).Nonempty) ∧ 0 < pN N := by
  obtain ⟨i₀, _⟩ := hw
  have h1 : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  have h2 : ∀ᶠ N in atTop, (1 / 2 : ℝ) < (qN N : ℝ) / dN N :=
    hq.eventually (lt_mem_nhds (by norm_num))
  have h3 : ∀ᶠ N in atTop, ∀ i, 0 < cN blk dN N i := by
    rw [eventually_all]
    intro i
    exact (hn i).eventually (lt_mem_nhds (hc i))
  filter_upwards [h1, h2, h3] with N hd1 hq1 hn1
  have hqpos : 0 < qN N := by
    by_contra hcon
    rw [not_lt] at hcon
    have : qN N = 0 := Nat.le_zero.mp hcon
    rw [this] at hq1
    simp at hq1
    linarith
  have hne : ∀ i, (blockSet (blk N) i).Nonempty := by
    intro i
    have := hn1 i
    rw [cN] at this
    have hcard : (0 : ℝ) < (blockSet (blk N) i).card := by
      by_contra hcon
      rw [not_lt] at hcon
      have h0 : ((blockSet (blk N) i).card : ℝ) = 0 := le_antisymm hcon (by positivity)
      rw [h0, zero_div] at this
      exact lt_irrefl _ this
    exact Finset.card_pos.mp (by exact_mod_cast hcard)
  exact ⟨hd1, hqpos, hne, Fin.pos (hne i₀).choose⟩

include hd in
theorem tendsto_lipG (S : ℝ) : Tendsto (fun N => lipG z (dN N) S) atTop (𝓝 0) := by
  have hdR : Tendsto (fun N => ((dN N : ℝ))) atTop atTop := tendsto_natCast_atTop_atTop.comp hd
  have h1 : Tendsto (fun N => 4 * S * (z.im + ‖z‖) / ((dN N : ℝ) * z.im ^ 4)) atTop (𝓝 0) := by
    have : (fun N => 4 * S * (z.im + ‖z‖) / ((dN N : ℝ) * z.im ^ 4))
        = fun N => (4 * S * (z.im + ‖z‖) / z.im ^ 4) / (dN N : ℝ) := by
      funext N; field_simp
    rw [this]
    exact Filter.Tendsto.div_atTop tendsto_const_nhds hdR
  have := (Real.continuous_sqrt.tendsto 0).comp h1
  simpa [lipG, Function.comp_def] using this

include hd hn in
omit hc hw hz hq in
theorem tendsto_lipS (S : ℝ) : Tendsto (fun N => lipS z (dN N) (pN N) S) atTop (𝓝 0) := by
  have hdR : Tendsto (fun N => ((dN N : ℝ))) atTop atTop := tendsto_natCast_atTop_atTop.comp hd
  have hp : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 (∑ i, c i)) := by
    have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) => hn i
    simpa [sum_cN] using this
  have h2 : Tendsto (fun N => (1 : ℝ) / (dN N : ℝ) ^ 2) atTop (𝓝 0) := by
    have := Filter.Tendsto.div_atTop (tendsto_const_nhds (x := (1 : ℝ)))
      ((tendsto_pow_atTop (α := ℝ) (n := 2) (by norm_num)).comp hdR)
    simpa [Function.comp_def] using this
  have h1 : Tendsto (fun N => 4 * S * (pN N : ℝ) * (z.im + ‖z‖) / ((dN N : ℝ) ^ 3 * z.im ^ 4))
      atTop (𝓝 0) := by
    have hev : ∀ᶠ N in atTop, 4 * S * (pN N : ℝ) * (z.im + ‖z‖) / ((dN N : ℝ) ^ 3 * z.im ^ 4)
        = (4 * S * (z.im + ‖z‖) / z.im ^ 4) * ((pN N : ℝ) / dN N) * (1 / (dN N : ℝ) ^ 2) := by
      filter_upwards [hd.eventually_gt_atTop 0] with N hN
      have : (dN N : ℝ) ≠ 0 := by positivity
      field_simp
    refine Tendsto.congr' (hev.mono fun N h => h.symm) ?_
    have := ((tendsto_const_nhds (x := 4 * S * (z.im + ‖z‖) / z.im ^ 4)).mul hp).mul h2
    simpa using this
  have := (Real.continuous_sqrt.tendsto 0).comp h1
  simpa [lipS, Function.comp_def] using this

include hc hw hd hq hn in
theorem tendsto_epsN (hz : 0 < z.im) (i : Fin M) :
    Tendsto (fun N => epsN w (blk N) (qN N) (dN N) z i) atTop (𝓝 0) := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  have hdR : Tendsto (fun N => ((dN N : ℝ))) atTop atTop := tendsto_natCast_atTop_atTop.comp hd
  set R : ℕ → ℝ := fun N => w i ^ 2 / dN N * (1 / z.im + ‖z‖ / z.im ^ 2)
    + w i ^ 2 * ‖z‖ * (9 * lipS z (dN N) (pN N) (Scalars.wSqMax w) / z.im) with hR
  have hRto : Tendsto R atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => w i ^ 2 / (dN N : ℝ) * (1 / z.im + ‖z‖ / z.im ^ 2))
        atTop (𝓝 0) := by
      have := (Filter.Tendsto.div_atTop (tendsto_const_nhds (x := w i ^ 2)) hdR).mul_const
        (1 / z.im + ‖z‖ / z.im ^ 2)
      simpa using this
    have h2 : Tendsto (fun N => w i ^ 2 * ‖z‖
        * (9 * lipS z (dN N) (pN N) (Scalars.wSqMax w) / z.im)) atTop (𝓝 0) := by
      have := ((tendsto_lipS hd hn (z := z) (Scalars.wSqMax w)).const_mul 9).div_const z.im
        |>.const_mul (w i ^ 2 * ‖z‖)
      simpa using this
    simpa [hR] using h1.add h2
  rw [tendsto_iff_norm_sub_tendsto_zero]
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hRto
    (Filter.Eventually.of_forall fun N => norm_nonneg _) ?_
  filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N ⟨hdN, hqN, hJ, hpN⟩
  rw [sub_zero]
  exact norm_epsN_le hz hpN hqN hdN (tauOf_sq_le w (blk N)) hS0 i (hJ i)

end Sequence


section Core

variable {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

/-- The eventual lower bound on `|E s_N|`: `1 / (2 |z| (1 + (∑ c_i w_i² + 1)/η))`. -/
noncomputable def mLow (c w : Fin M → ℝ) (z : ℂ) : ℝ :=
  1 / (2 * ‖z‖ * (1 + (∑ i, c i * w i ^ 2 + 1) / z.im))

theorem mLow_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) {z : ℂ} (hz : 0 < z.im) (hz0 : z ≠ 0) :
    0 < mLow c w z := by
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hSg : 0 ≤ ∑ i, c i * w i ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (hc i).le (sq_nonneg _)
  unfold mLow
  positivity

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))

include hc hw hz hd hq hn in
/-- **The eventual bounds and the key estimate** (plan section 3.4). Along `N`:
`|1 + w_i² S_N| ≥ η/(2|z|)`, `mLow ≤ |S_N| ≤ 2/η`, `Im S_N ≥ 0`, the Stein mean identity,
and `|zfunC(S_N) - z| ≤ R_N` with `R_N → 0`, where `S_N = E s_N`. -/
theorem eventually_bounds_meanS : ∃ R : ℕ → ℝ, Tendsto R atTop (𝓝 0) ∧ ∀ᶠ N in atTop,
    (∀ i, z.im / (2 * ‖z‖) ≤ ‖1 + ((w i : ℝ) : ℂ) ^ 2 * meanS w (blk N) (qN N) (dN N) z‖) ∧
    mLow c w z ≤ ‖meanS w (blk N) (qN N) (dN N) z‖ ∧
    ‖meanS w (blk N) (qN N) (dN N) z‖ ≤ 2 / z.im ∧
    0 ≤ (meanS w (blk N) (qN N) (dN N) z).im ∧
    ‖zfunC c w (meanS w (blk N) (qN N) (dN N) z) - z‖ ≤ R N ∧
    (∀ i, 1 + z * meanG w (blk N) (qN N) (dN N) z i
      = -(((w i : ℝ) : ℂ) ^ 2 * z) * (meanS w (blk N) (qN N) (dN N) z
          * meanG w (blk N) (qN N) (dN N) z i) + epsN w (blk N) (qN N) (dN N) z i) := by
  have hz0 : z ≠ 0 := ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hm : 0 < mLow c w z := mLow_pos hc hz hz0
  set η : ℝ := z.im with hη
  set L : ℝ := 2 * ‖z‖ / η with hL
  have hLpos : 0 < L := by positivity
  set Sig : ℝ := ∑ i, c i * w i ^ 2 + 1 with hSig
  have hSgnn : 0 ≤ ∑ i, c i * w i ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (hc i).le (sq_nonneg _)
  set R : ℕ → ℝ := fun N => (|1 - (qN N : ℝ) / dN N|
    + L * ∑ i, cN blk dN N i * ‖epsN w (blk N) (qN N) (dN N) z i‖
    + (2 / η) * L * ∑ i, |cN blk dN N i - c i| * w i ^ 2) / mLow c w z with hR
  have h2 : Tendsto (fun N => ∑ i, cN blk dN N i * ‖epsN w (blk N) (qN N) (dN N) z i‖)
      atTop (𝓝 0) := by
    have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
      (hn i).mul (tendsto_epsN hc hw hd hq hn hz i).norm
    simpa using this
  have h3 : Tendsto (fun N => ∑ i, |cN blk dN N i - c i| * w i ^ 2) atTop (𝓝 0) := by
    have := tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) =>
      ((hn i).sub_const (c i)).abs.mul_const (w i ^ 2)
    simpa using this
  have h5 : Tendsto (fun N => ∑ i, cN blk dN N i * w i ^ 2) atTop (𝓝 (∑ i, c i * w i ^ 2)) :=
    tendsto_finsetSum Finset.univ fun i (_ : i ∈ Finset.univ) => (hn i).mul_const (w i ^ 2)
  refine ⟨R, ?_, ?_⟩
  · have h1 : Tendsto (fun N => |1 - (qN N : ℝ) / dN N|) atTop (𝓝 0) := by
      have := (hq.const_sub 1).abs
      simpa using this
    have := ((h1.add (h2.const_mul L)).add (h3.const_mul ((2 / η) * L))).div_const (mLow c w z)
    simpa [hR] using this
  · have hev1 := eventually_sizes_pos hc hw hd hq hn
    have hev2 : ∀ᶠ N in atTop, ∀ i, ‖epsN w (blk N) (qN N) (dN N) z i‖ ≤ 1 / 2 := by
      rw [eventually_all]
      intro i
      have := (tendsto_epsN hc hw hd hq hn hz i).norm
      rw [norm_zero] at this
      exact (this.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 2))).mono fun N h => h.le
    have hev3 : ∀ᶠ N in atTop, (qN N : ℝ) / dN N ∈ Set.Ioo (3 / 4 : ℝ) 2 :=
      hq.eventually (Ioo_mem_nhds (by norm_num) (by norm_num))
    have hev4 : ∀ᶠ N in atTop,
        ∑ i, cN blk dN N i * ‖epsN w (blk N) (qN N) (dN N) z i‖ ≤ 1 / 4 :=
      (h2.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 4))).mono fun N h => h.le
    have hev5 : ∀ᶠ N in atTop, ∑ i, cN blk dN N i * w i ^ 2 ≤ Sig :=
      (h5.eventually (gt_mem_nhds (by rw [hSig]; linarith))).mono fun N h => h.le
    filter_upwards [hev1, hev2, hev3, hev4, hev5] with N hN hε hκ hεsum hSgN
    obtain ⟨hdN, hqN, hJ, hpN⟩ := hN
    set S := meanS w (blk N) (qN N) (dN N) z with hS
    set E : Fin M → ℂ := fun i => meanG w (blk N) (qN N) (dN N) z i with hE
    set ε : Fin M → ℂ := fun i => epsN w (blk N) (qN N) (dN N) z i with hεdef
    have ha : ∀ i, 1 + z * E i = -(((w i : ℝ) : ℂ) ^ 2 * z) * (S * E i) + ε i :=
      fun i => mean_stein hz hdN i (hJ i)
    have hb : z * S + (((qN N : ℝ) / dN N : ℝ) : ℂ)
        = ∑ i, ((cN blk dN N i : ℝ) : ℂ) * (1 + z * E i) := mean_trace hz hpN hdN
    have hEb : ∀ i, ‖E i‖ ≤ 1 / η := fun i => norm_meanG_le hz i
    -- (1) the lower bound on `|1 + w² S|`
    have h1 : ∀ i, η / (2 * ‖z‖) ≤ ‖1 + ((w i : ℝ) : ℂ) ^ 2 * S‖ := by
      intro i
      have := norm_one_add_mul_ge hz0 hz (ha i) (hEb i)
      calc η / (2 * ‖z‖) = (1 - 1 / 2) * η / ‖z‖ := by ring
        _ ≤ (1 - ‖ε i‖) * η / ‖z‖ := by
            have := hε i
            gcongr
        _ ≤ _ := this
    -- (2) the upper bound on `|S|`
    have h2 : ‖S‖ ≤ 2 / η := by
      refine (norm_meanS_le hz).trans ?_
      rw [← div_div]
      exact div_le_div_of_nonneg_right hκ.2.le hz.le
    -- (3) the lower bound on `|S|`
    have h3 : mLow c w z ≤ ‖S‖ := by
      have hge := norm_S_ge (cN_nonneg blk dN N) E ε ha hb
      rw [Complex.ofReal_re] at hge
      have hsumE : ∑ i, cN blk dN N i * w i ^ 2 * ‖E i‖ ≤ Sig / η := by
        calc ∑ i, cN blk dN N i * w i ^ 2 * ‖E i‖
            ≤ ∑ i, cN blk dN N i * w i ^ 2 * (1 / η) := by
              refine Finset.sum_le_sum fun i _ => ?_
              have := cN_nonneg blk dN N i
              exact mul_le_mul_of_nonneg_left (hEb i) (by positivity)
          _ = (∑ i, cN blk dN N i * w i ^ 2) / η := by
              rw [Finset.sum_div]
              refine Finset.sum_congr rfl fun i _ => ?_
              ring
          _ ≤ Sig / η := div_le_div_of_nonneg_right hSgN hz.le
      have hkey : 1 / 2 ≤ ‖z‖ * ‖S‖ * (1 + Sig / η) := by
        have e1 : ‖z‖ * ‖S‖ * (1 + ∑ i, cN blk dN N i * w i ^ 2 * ‖E i‖)
            ≤ ‖z‖ * ‖S‖ * (1 + Sig / η) :=
          mul_le_mul_of_nonneg_left (by linarith) (by positivity)
        linarith [hκ.1]
      rw [mLow, div_le_iff₀ (by positivity)]
      nlinarith
    -- (4) `Im S ≥ 0`
    have h4 : 0 ≤ S.im := im_meanS_nonneg hz hpN hdN
    -- (5) the key estimate
    have hS0 : S ≠ 0 := by
      intro h
      rw [h, norm_zero] at h3
      linarith
    have hSne : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * S ≠ 0 := by
      intro i h
      have := h1 i
      rw [h, norm_zero] at this
      have : 0 < η / (2 * ‖z‖) := by positivity
      linarith
    have hinv : ∀ i, ‖(1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹‖ ≤ L := by
      intro i
      rw [norm_inv, hL, ← inv_div η (2 * ‖z‖)]
      exact inv_anti₀ (by positivity) (h1 i)
    have hid := mul_sub_zfunC_eq c w hS0 hSne (cN blk dN N) _ E ε ha hb
    have hnormid : ‖S * (z - zfunC c w S)‖
        ≤ |1 - (qN N : ℝ) / dN N| + L * ∑ i, cN blk dN N i * ‖ε i‖
          + (2 / η) * L * ∑ i, |cN blk dN N i - c i| * w i ^ 2 := by
      rw [hid]
      refine (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add ?_ ?_)) ?_)
      · rw [show (1 : ℂ) - (((qN N : ℝ) / dN N : ℝ) : ℂ) = ((1 - (qN N : ℝ) / dN N : ℝ) : ℂ) by
          push_cast; ring, Complex.norm_real, Real.norm_eq_abs]
      · refine (norm_sum_le _ _).trans ?_
        rw [Finset.mul_sum]
        refine Finset.sum_le_sum fun i _ => ?_
        rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
          abs_of_nonneg (cN_nonneg blk dN N i)]
        calc cN blk dN N i * ‖ε i‖ * ‖(1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹‖
            ≤ cN blk dN N i * ‖ε i‖ * L :=
              mul_le_mul_of_nonneg_left (hinv i) (by have := cN_nonneg blk dN N i; positivity)
          _ = L * (cN blk dN N i * ‖ε i‖) := by ring
      · rw [norm_mul]
        have hsum : ‖∑ i, (((cN blk dN N i : ℝ) : ℂ) - ((c i : ℝ) : ℂ)) * ((w i : ℝ) : ℂ) ^ 2
            * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹‖ ≤ L * ∑ i, |cN blk dN N i - c i| * w i ^ 2 := by
          refine (norm_sum_le _ _).trans ?_
          rw [Finset.mul_sum]
          refine Finset.sum_le_sum fun i _ => ?_
          rw [norm_mul, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs,
            ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
          calc |cN blk dN N i - c i| * w i ^ 2 * ‖(1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹‖
              ≤ |cN blk dN N i - c i| * w i ^ 2 * L :=
                mul_le_mul_of_nonneg_left (hinv i) (by positivity)
            _ = L * (|cN blk dN N i - c i| * w i ^ 2) := by ring
        calc ‖S‖ * ‖∑ i, (((cN blk dN N i : ℝ) : ℂ) - ((c i : ℝ) : ℂ)) * ((w i : ℝ) : ℂ) ^ 2
              * (1 + ((w i : ℝ) : ℂ) ^ 2 * S)⁻¹‖
            ≤ (2 / η) * (L * ∑ i, |cN blk dN N i - c i| * w i ^ 2) :=
              mul_le_mul h2 hsum (norm_nonneg _) (by positivity)
          _ = (2 / η) * L * ∑ i, |cN blk dN N i - c i| * w i ^ 2 := by ring
    have h5 : ‖zfunC c w S - z‖ ≤ R N := by
      have hSn : 0 < ‖S‖ := norm_pos_iff.mpr hS0
      have e : ‖zfunC c w S - z‖ = ‖S * (z - zfunC c w S)‖ / ‖S‖ := by
        rw [norm_mul, norm_sub_rev]
        field_simp
      rw [e, hR]
      have hnn1 : 0 ≤ ∑ i, cN blk dN N i * ‖ε i‖ :=
        Finset.sum_nonneg fun i _ => mul_nonneg (cN_nonneg _ _ _ _) (norm_nonneg _)
      have hnn2 : 0 ≤ ∑ i, |cN blk dN N i - c i| * w i ^ 2 :=
        Finset.sum_nonneg fun i _ => by positivity
      exact div_le_div₀ (by positivity) hnormid hm h3
    exact ⟨h1, h3, h2, h4, h5, ha⟩

include hc hw hz hd hq hn in
theorem tendsto_zfunC_meanS :
    Tendsto (fun N => zfunC c w (meanS w (blk N) (qN N) (dN N) z)) atTop (𝓝 z) := by
  obtain ⟨R, hR, hev⟩ := eventually_bounds_meanS hc hw hz hd hq hn
  rw [tendsto_iff_norm_sub_tendsto_zero]
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hR
    (Filter.Eventually.of_forall fun N => norm_nonneg _) ?_
  filter_upwards [hev] with N hN
  exact hN.2.2.2.2.1

include hc hw hz hd hq hn in
/-- **Existence of the root in `ℂ⁺`**, from the bounded sequence `E s_N` (Bolzano-Weierstrass)
and the closed conditions it satisfies eventually. -/
theorem exists_isRoot : ∃ σ, IsRoot c w z σ := by
  obtain ⟨R, hR, hev⟩ := eventually_bounds_meanS hc hw hz hd hq hn
  have hz0 : z ≠ 0 := ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hm : 0 < mLow c w z := mLow_pos hc hz hz0
  have hf := tendsto_zfunC_meanS hc hw hz hd hq hn
  set S : ℕ → ℂ := fun N => meanS w (blk N) (qN N) (dN N) z with hSdef
  have hfreq : ∃ᶠ N in atTop, S N ∈ Metric.closedBall (0 : ℂ) (2 / z.im) := by
    refine (hev.mono fun N hN => ?_).frequently
    rw [Metric.mem_closedBall, dist_zero_right]
    exact hN.2.2.1
  obtain ⟨σ, -, φ, hφ, hlim⟩ :=
    tendsto_subseq_of_frequently_bounded Metric.isBounded_closedBall hfreq
  have hφtop : Tendsto φ atTop atTop := hφ.tendsto_atTop
  have hevφ : ∀ᶠ k in atTop, (∀ i, z.im / (2 * ‖z‖) ≤ ‖1 + ((w i : ℝ) : ℂ) ^ 2 * S (φ k)‖) ∧
      mLow c w z ≤ ‖S (φ k)‖ ∧ 0 ≤ (S (φ k)).im :=
    (hφtop.eventually hev).mono fun k hk => ⟨hk.1, hk.2.1, hk.2.2.2.1⟩
  have hnorm : mLow c w z ≤ ‖σ‖ :=
    ge_of_tendsto ((continuous_norm.tendsto σ).comp hlim) (hevφ.mono fun k hk => hk.2.1)
  have him : 0 ≤ σ.im :=
    ge_of_tendsto ((Complex.continuous_im.tendsto σ).comp hlim) (hevφ.mono fun k hk => hk.2.2)
  have hone : ∀ i, z.im / (2 * ‖z‖) ≤ ‖1 + ((w i : ℝ) : ℂ) ^ 2 * σ‖ := by
    intro i
    have hcont : Continuous fun s : ℂ => ‖1 + ((w i : ℝ) : ℂ) ^ 2 * s‖ := by fun_prop
    exact ge_of_tendsto ((hcont.tendsto σ).comp hlim) (hevφ.mono fun k hk => hk.1 i)
  have hσ0 : σ ≠ 0 := by
    intro h
    rw [h, norm_zero] at hnorm
    linarith
  have hσne : ∀ i, 1 + ((w i : ℝ) : ℂ) ^ 2 * σ ≠ 0 := by
    intro i h
    have := hone i
    rw [h, norm_zero] at this
    have : 0 < z.im / (2 * ‖z‖) := by positivity
    linarith
  have hroot : zfunC c w σ = z :=
    tendsto_nhds_unique ((continuousAt_zfunC c w hσ0 hσne).tendsto.comp hlim) (hf.comp hφtop)
  have hσim : 0 < σ.im := by
    rcases him.lt_or_eq with h | h
    · exact h
    · exfalso
      have := im_zfunC c w σ
      rw [hroot, ← h, zero_mul] at this
      linarith
  exact ⟨σ, hσim, hroot⟩

include hc hw hz hd hq hn in
/-- **The trace law for the means**: `E s_N → sGlob c w z`. -/
theorem tendsto_meanS :
    Tendsto (fun N => meanS w (blk N) (qN N) (dN N) z) atTop (𝓝 (sGlob c w z)) := by
  obtain ⟨R, hR, hev⟩ := eventually_bounds_meanS hc hw hz hd hq hn
  have hroot := isRoot_sGlob (exists_isRoot hc hw hz hd hq hn)
  have hf := tendsto_zfunC_meanS hc hw hz hd hq hn
  have hevim : ∀ᶠ N in atTop, 0 < (zfunC c w (meanS w (blk N) (qN N) (dN N) z)).im :=
    ((Complex.continuous_im.tendsto z).comp hf).eventually (lt_mem_nhds hz)
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have hR' : Tendsto (fun N => (4 / z.im ^ 2) * R N) atTop (𝓝 0) := by
    simpa using hR.const_mul (4 / z.im ^ 2)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hR'
    (Filter.Eventually.of_forall fun N => norm_nonneg _) ?_
  filter_upwards [hev, hevim] with N hN him
  obtain ⟨-, -, h2, h4, h5, -⟩ := hN
  have hstab := norm_sub_le_of_root (fun i => (hc i).le) hz h4 him hroot.1 hroot.2
  refine hstab.trans ?_
  rw [div_le_iff₀ hz]
  calc 2 * ‖meanS w (blk N) (qN N) (dN N) z‖ * ‖zfunC c w (meanS w (blk N) (qN N) (dN N) z) - z‖
      ≤ 2 * (2 / z.im) * R N := by gcongr
    _ = 4 / z.im ^ 2 * R N * z.im := by field_simp; ring

include hc hw hz hd hq hn in
/-- **The block trace law for the means**: `E g_{J_i} → gC w i z (sGlob c w z)`. -/
theorem tendsto_meanG (i : Fin M) :
    Tendsto (fun N => meanG w (blk N) (qN N) (dN N) z i) atTop
      (𝓝 (gC w i z (sGlob c w z))) := by
  obtain ⟨R, hR, hev⟩ := eventually_bounds_meanS hc hw hz hd hq hn
  have hroot := isRoot_sGlob (exists_isRoot hc hw hz hd hq hn)
  have hS := tendsto_meanS hc hw hz hd hq hn
  have hε := tendsto_epsN hc hw hd hq hn hz i
  have hz0 : z ≠ 0 := ne_zero_of_im_pos hz
  have hden : z * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w z) ≠ 0 :=
    mul_ne_zero hz0 (one_add_mul_ne_zero_of_im_pos hroot.1 _)
  have hlim : Tendsto (fun N => (epsN w (blk N) (qN N) (dN N) z i - 1)
      * (z * (1 + ((w i : ℝ) : ℂ) ^ 2 * meanS w (blk N) (qN N) (dN N) z))⁻¹) atTop
      (𝓝 ((0 - 1) * (z * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w z))⁻¹)) :=
    (hε.sub_const 1).mul (Tendsto.inv₀
      (tendsto_const_nhds.mul (tendsto_const_nhds.add (tendsto_const_nhds.mul hS))) hden)
  have hval : (0 - 1 : ℂ) * (z * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w z))⁻¹
      = gC w i z (sGlob c w z) := by
    rw [gC]; ring
  rw [hval] at hlim
  refine hlim.congr' ?_
  filter_upwards [hev] with N hN
  obtain ⟨h1, -, -, -, -, ha⟩ := hN
  have hne : 1 + ((w i : ℝ) : ℂ) ^ 2 * meanS w (blk N) (qN N) (dN N) z ≠ 0 := by
    intro h
    have := h1 i
    rw [h, norm_zero] at this
    have : 0 < z.im / (2 * ‖z‖) := by have := norm_pos_iff.mpr hz0; positivity
    linarith
  have hid : z * meanG w (blk N) (qN N) (dN N) z i
      * (1 + ((w i : ℝ) : ℂ) ^ 2 * meanS w (blk N) (qN N) (dN N) z)
      = epsN w (blk N) (qN N) (dN N) z i - 1 := by
    linear_combination ha i
  rw [← hid]
  field_simp

end Core


/-! ### Concentration (`measure_abs_ge_le_of_lipschitz`, choice 31 of R1) -/

section Concentration

variable {p q : ℕ}

/-- Subtracting a constant keeps the Lipschitz constant. -/
theorem lipschitzWith_sub_const {E : Type*} [PseudoMetricSpace E] {K : ℝ≥0} {F : E → ℝ}
    (hF : LipschitzWith K F) (a : ℝ) : LipschitzWith K (fun x => F x - a) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  have := hF.dist_le_mul x y
  simpa [Real.dist_eq] using this

/-- An average of `K`-Lipschitz functions is `K`-Lipschitz. -/
theorem lipschitzWith_avg {E : Type*} [PseudoMetricSpace E] {K : ℝ≥0} {ι : Type*}
    (J : Finset ι) {F : ι → E → ℝ} (hF : ∀ j ∈ J, LipschitzWith K (F j)) :
    LipschitzWith K (fun x => ((J.card : ℝ))⁻¹ * ∑ j ∈ J, F j x) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  rcases J.eq_empty_or_nonempty with hJ | hJ
  · simp only [hJ, Finset.card_empty, Nat.cast_zero, inv_zero, Finset.sum_empty, mul_zero,
      dist_self]
    positivity
  have hc : (0 : ℝ) < J.card := by exact_mod_cast Finset.card_pos.mpr hJ
  rw [Real.dist_eq, ← mul_sub, ← Finset.sum_sub_distrib, abs_mul, abs_of_pos (inv_pos.mpr hc)]
  calc ((J.card : ℝ))⁻¹ * |∑ j ∈ J, (F j x - F j y)|
      ≤ ((J.card : ℝ))⁻¹ * ∑ j ∈ J, (K : ℝ) * dist x y := by
        refine mul_le_mul_of_nonneg_left ((Finset.abs_sum_le_sum_abs _ _).trans
          (Finset.sum_le_sum fun j hj => ?_)) (by positivity)
        have := (hF j hj).dist_le_mul x y
        rwa [Real.dist_eq] at this
    _ = (K : ℝ) * dist x y := by
        rw [Finset.sum_const, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hc.ne', one_mul]

/-- The tail bound of `measure_abs_ge_le_of_lipschitz` for a model function `f = F ∘
matrixEquivE`, in `ℝ≥0∞`. -/
theorem measure_abs_sub_mean_ge_le (hpq : 0 < p * q) {L : ℝ} (hL : 0 < L)
    {F : EuclideanSpace ℝ (Fin (p * q)) → ℝ} (hF : LipschitzWith (Real.toNNReal L) F)
    {f : Matrix (Fin p) (Fin q) ℝ → ℝ} (hfF : ∀ B, f B = F (matrixEquivE p q B)) {t : ℝ}
    (ht : 0 < t) :
    gaussianMatrix p q {B | t ≤ |f B - ∫ B', f B' ∂gaussianMatrix p q|}
      ≤ ENNReal.ofReal (2 * Real.exp (-t ^ 2 / (2 * L ^ 2))) := by
  have hL' : (0 : ℝ≥0) < Real.toNNReal L := Real.toNNReal_pos.mpr hL
  have h := measure_abs_ge_le_of_lipschitz (p := p) (d := q) hpq hL' hF ht
  rw [Real.coe_toNNReal L hL.le] at h
  have heq : {B | t ≤ |f B - ∫ B', f B' ∂gaussianMatrix p q|}
      = {Z | t ≤ |F (matrixEquivE p q Z) - ∫ Z', F (matrixEquivE p q Z') ∂gaussianMatrix p q|} := by
    ext B
    simp only [Set.mem_ofPred_eq, hfF]
  rw [heq, ← ENNReal.ofReal_toReal (measure_ne_top (gaussianMatrix p q) _)]
  exact ENNReal.ofReal_le_ofReal h

/-- The complex version: real and imaginary parts each `L`-Lipschitz through `matrixEquivE`. -/
theorem measure_norm_sub_mean_ge_le (hpq : 0 < p * q) {L : ℝ} (hL : 0 < L)
    {Fre Fim : EuclideanSpace ℝ (Fin (p * q)) → ℝ} (hre : LipschitzWith (Real.toNNReal L) Fre)
    (him : LipschitzWith (Real.toNNReal L) Fim) {g : Matrix (Fin p) (Fin q) ℝ → ℂ}
    (hgI : Integrable g (gaussianMatrix p q))
    (hgre : ∀ B, (g B).re = Fre (matrixEquivE p q B))
    (hgim : ∀ B, (g B).im = Fim (matrixEquivE p q B))
    {t : ℝ} (ht : 0 < t) :
    gaussianMatrix p q {B | t ≤ ‖g B - ∫ B', g B' ∂gaussianMatrix p q‖}
      ≤ ENNReal.ofReal (4 * Real.exp (-(t / 2) ^ 2 / (2 * L ^ 2))) := by
  set m : ℂ := ∫ B', g B' ∂gaussianMatrix p q with hm
  have hmre : (∫ B', (g B').re ∂gaussianMatrix p q) = m.re := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.reCLM hgI
    simpa [hm] using h
  have hmim : (∫ B', (g B').im ∂gaussianMatrix p q) = m.im := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM hgI
    simpa [hm] using h
  have h1 := measure_abs_sub_mean_ge_le hpq hL hre hgre (half_pos ht)
  have h2 := measure_abs_sub_mean_ge_le hpq hL him hgim (half_pos ht)
  rw [hmre] at h1
  rw [hmim] at h2
  have hsub : {B | t ≤ ‖g B - m‖}
      ⊆ {B | t / 2 ≤ |(g B).re - m.re|} ∪ {B | t / 2 ≤ |(g B).im - m.im|} := by
    intro B hB
    by_contra hcon
    simp only [Set.mem_union, Set.mem_ofPred_eq, not_or, not_le] at hcon
    have hsplit := Complex.norm_le_abs_re_add_abs_im (g B - m)
    simp only [Complex.sub_re, Complex.sub_im] at hsplit
    have : t ≤ ‖g B - m‖ := hB
    linarith [hcon.1, hcon.2]
  calc gaussianMatrix p q {B | t ≤ ‖g B - m‖}
      ≤ gaussianMatrix p q ({B | t / 2 ≤ |(g B).re - m.re|} ∪ {B | t / 2 ≤ |(g B).im - m.im|}) :=
        measure_mono hsub
    _ ≤ gaussianMatrix p q {B | t / 2 ≤ |(g B).re - m.re|}
          + gaussianMatrix p q {B | t / 2 ≤ |(g B).im - m.im|} := measure_union_le _ _
    _ ≤ ENNReal.ofReal (2 * Real.exp (-(t / 2) ^ 2 / (2 * L ^ 2)))
          + ENNReal.ofReal (2 * Real.exp (-(t / 2) ^ 2 / (2 * L ^ 2))) := add_le_add h1 h2
    _ = ENNReal.ofReal (4 * Real.exp (-(t / 2) ^ 2 / (2 * L ^ 2))) := by
        rw [← ENNReal.ofReal_add (by positivity) (by positivity)]
        ring_nf

end Concentration

/-! ### From the means to convergence in probability (the R1a pattern) -/

section Transfer

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN : ℕ → ℕ}

/-- **Generic transfer.** A complex model statistic whose real and imaginary parts are
`L_N`-Lipschitz through `matrixEquivE`, with `L_N → 0` and mean `→ ℓ`, converges to `ℓ` in
probability in the norm form. -/
theorem tendstoInProb_of_lipschitz_mean (B : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (qN N)) ℝ)
    (hB : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (qN N)) (μ N))
    {g : (N : ℕ) → Matrix (Fin (pN N)) (Fin (qN N)) ℝ → ℂ} {L : ℕ → ℝ} {ℓ : ℂ}
    (hL0 : Tendsto L atTop (𝓝 0))
    (hev : ∀ᶠ N in atTop, 0 < pN N * qN N ∧ 0 < L N ∧ Measurable (g N) ∧
      Integrable (g N) (gaussianMatrix (pN N) (qN N)) ∧
      ∃ Fre Fim : EuclideanSpace ℝ (Fin (pN N * qN N)) → ℝ,
        LipschitzWith (Real.toNNReal (L N)) Fre ∧ LipschitzWith (Real.toNNReal (L N)) Fim ∧
        ∀ B, (g N B).re = Fre (matrixEquivE (pN N) (qN N) B) ∧
          (g N B).im = Fim (matrixEquivE (pN N) (qN N) B))
    (hmean : Tendsto (fun N => ∫ B, g N B ∂gaussianMatrix (pN N) (qN N)) atTop (𝓝 ℓ)) :
    TendstoInProb μ (fun N ω => ‖g N (B N ω) - ℓ‖) 0 := by
  intro ε hε
  set m : ℕ → ℂ := fun N => ∫ B, g N B ∂gaussianMatrix (pN N) (qN N) with hm
  have hev2 : ∀ᶠ N in atTop, ‖m N - ℓ‖ < ε / 2 := by
    have := (tendsto_iff_norm_sub_tendsto_zero.mp hmean).eventually
      (gt_mem_nhds (by positivity : (0 : ℝ) < ε / 2))
    exact this
  set Bd : ℕ → ℝ≥0∞ :=
    fun N => ENNReal.ofReal (4 * Real.exp (-(ε / 2 / 2) ^ 2 / (2 * L N ^ 2))) with hBd
  have hBdto : Tendsto Bd atTop (𝓝 0) := by
    have hLpos : ∀ᶠ N in atTop, 0 < L N := hev.mono fun N h => h.2.1
    have hL2 : Tendsto (fun N => 2 * L N ^ 2) atTop (𝓝[>] 0) := by
      rw [tendsto_nhdsWithin_iff]
      refine ⟨by simpa using (hL0.pow 2).const_mul 2, hLpos.mono fun N h => ?_⟩
      simp only [Set.mem_Ioi]
      positivity
    have hinv : Tendsto (fun N => (2 * L N ^ 2)⁻¹) atTop atTop :=
      tendsto_inv_nhdsGT_zero.comp hL2
    have hneg : -(ε / 2 / 2) ^ 2 < 0 := by
      have : (0 : ℝ) < (ε / 2 / 2) ^ 2 := by positivity
      linarith
    have hbot : Tendsto (fun N => -(ε / 2 / 2) ^ 2 / (2 * L N ^ 2)) atTop atBot := by
      simp only [div_eq_mul_inv]
      exact Filter.Tendsto.const_mul_atTop_of_neg hneg hinv
    have hexp := Real.tendsto_exp_atBot.comp hbot
    have h4 : Tendsto (fun N => 4 * Real.exp (-(ε / 2 / 2) ^ 2 / (2 * L N ^ 2))) atTop (𝓝 0) := by
      simpa [Function.comp_def] using hexp.const_mul 4
    have := ENNReal.tendsto_ofReal h4
    simpa [hBd] using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hBdto
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  filter_upwards [hev, hev2] with N hN hmN
  obtain ⟨hpq, hLN, hmeas, hint, Fre, Fim, hre, him, hF⟩ := hN
  have hsub : {ω | ε ≤ |‖g N (B N ω) - ℓ‖ - 0|} ⊆ {ω | ε / 2 ≤ ‖g N (B N ω) - m N‖} := by
    intro ω hω
    have h1 : ε ≤ ‖g N (B N ω) - ℓ‖ := by
      have := hω
      simp only [Set.mem_ofPred_eq, sub_zero, abs_norm] at this
      exact this
    have htri : ‖g N (B N ω) - ℓ‖ ≤ ‖g N (B N ω) - m N‖ + ‖m N - ℓ‖ := by
      have := norm_add_le (g N (B N ω) - m N) (m N - ℓ)
      simpa using this
    change ε / 2 ≤ ‖g N (B N ω) - m N‖
    linarith
  have hmeasSet : MeasurableSet {x : Matrix (Fin (pN N)) (Fin (qN N)) ℝ | ε / 2 ≤ ‖g N x - m N‖} :=
    measurableSet_le measurable_const (hmeas.sub_const _).norm
  calc μ N {ω | ε ≤ |‖g N (B N ω) - ℓ‖ - 0|}
      ≤ μ N {ω | ε / 2 ≤ ‖g N (B N ω) - m N‖} := measure_mono hsub
    _ = gaussianMatrix (pN N) (qN N) {x | ε / 2 ≤ ‖g N x - m N‖} := (hB N).measure_eq hmeasSet
    _ ≤ Bd N := measure_norm_sub_mean_ge_le hpq hLN hre him hint (fun B => (hF B).1)
        (fun B => (hF B).2) (by positivity)

end Transfer


/-! ### The trace law and the block trace law in probability (item H6, plan section 3.4) -/

section Main

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))
  (B : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (qN N)) ℝ)
  (hB : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (qN N)) (μ N))

include hc hw hz hd hq hn hB in
/-- **The heteroscedastic trace law** (plan section 2.2): `s_N(z) → sGlob c w z` in probability,
for every `z` with `Im z > 0`. -/
theorem tendstoInProb_ssig :
    TendstoInProb μ (fun N ω => ‖ssig (tauOf w (blk N)) (dN N) z (B N ω) - sGlob c w z‖) 0 := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  refine tendstoInProb_of_lipschitz_mean B hB (g := fun N B => ssig (tauOf w (blk N)) (dN N) z B)
    (L := fun N => lipS z (dN N) (pN N) (Scalars.wSqMax w)) (tendsto_lipS hd hn _) ?_
    (tendsto_meanS hc hw hz hd hq hn)
  filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
  obtain ⟨hdN, hqN, -, hpN⟩ := hN
  set τ := tauOf w (blk N) with hτ
  have hS : ∀ j, τ j ^ 2 ≤ Scalars.wSqMax w := tauOf_sq_le w (blk N)
  have hL := lipS_pos (p := pN N) hz hpN hdN hS0
  set cst : ℂ := ((qN N : ℂ) - pN N) / ((dN N : ℂ) * z) with hcst
  refine ⟨Nat.mul_pos hpN hqN, hL, measurable_ssig hz hpN hdN τ, integrable_ssig hz hpN hdN τ,
    fun x => (sfun (zeta (pN N) (dN N) z) (lam τ (qN N) x)).re - cst.re,
    fun x => (sfun (zeta (pN N) (dN N) z) (lam τ (qN N) x)).im - cst.im, ?_, ?_, fun B => ?_⟩
  · exact lipschitzWith_sub_const (lipschitzWith_sfun_lam_re hz hpN hdN τ hS hS0.le
      (Real.le_coe_toNNReal _)) _
  · exact lipschitzWith_sub_const (lipschitzWith_sfun_lam_im hz hpN hdN τ hS hS0.le
      (Real.le_coe_toNNReal _)) _
  · rw [ssig_eq_sfun hz.ne' hpN hdN τ B, Complex.sub_re, Complex.sub_im]
    exact ⟨rfl, rfl⟩

include hc hw hz hd hq hn hB in
/-- **The block trace law** (plan section 2.2): `|J_i|⁻¹ ∑_{j ∈ J_i} G_σ(z)_{jj} →
-1/(z (1 + w_i² sGlob z))` in probability, for every `z` with `Im z > 0`. -/
theorem tendstoInProb_gsigAvg (i : Fin M) :
    TendstoInProb μ (fun N ω => ‖gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (B N ω)
      - gC w i z (sGlob c w z)‖) 0 := by
  have hS0 : 0 < Scalars.wSqMax w := wSqMax_pos hw
  refine tendstoInProb_of_lipschitz_mean B hB
    (g := fun N B => gsigAvg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) B)
    (L := fun N => lipG z (dN N) (Scalars.wSqMax w)) (tendsto_lipG hd _) ?_
    (tendsto_meanG hc hw hz hd hq hn i)
  filter_upwards [eventually_sizes_pos hc hw hd hq hn] with N hN
  obtain ⟨hdN, hqN, -, hpN⟩ := hN
  set τ := tauOf w (blk N) with hτ
  set J := blockSet (blk N) i with hJ
  have hS : ∀ j, τ j ^ 2 ≤ Scalars.wSqMax w := tauOf_sq_le w (blk N)
  have hL := lipG_pos (d := dN N) hz hdN hS0
  refine ⟨Nat.mul_pos hpN hqN, hL, measurable_gsigAvg hz hpN hdN τ J,
    integrable_gsigAvg hz hpN hdN τ J,
    fun x => ((J.card : ℝ))⁻¹ * ∑ j ∈ J, (gfun τ (dN N) z j x).re,
    fun x => ((J.card : ℝ))⁻¹ * ∑ j ∈ J, (gfun τ (dN N) z j x).im, ?_, ?_, fun B => ?_⟩
  · exact lipschitzWith_avg J fun j _ =>
      lipschitzWith_gfun_re hz hpN hdN τ hS hS0.le j (Real.le_coe_toNNReal _)
  · exact lipschitzWith_avg J fun j _ =>
      lipschitzWith_gfun_im hz hpN hdN τ hS hS0.le j (Real.le_coe_toNNReal _)
  · have hcast : ((J.card : ℂ))⁻¹ = ((((J.card : ℝ))⁻¹ : ℝ) : ℂ) := by push_cast; rfl
    constructor
    · simp only [gsigAvg, hcast, Complex.re_ofReal_mul, Complex.re_sum,
        gfun_matrixEquivE hz.ne' hpN hdN]
    · simp only [gsigAvg, hcast, Complex.im_ofReal_mul, Complex.im_sum,
        gfun_matrixEquivE hz.ne' hpN hdN]

end Main

/-! ### Identification with the column split of `RMT/Het/Split.lean` -/

section Model

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-- The block map of the stacked row index: row `r` belongs to table
`(finSigmaFinEquiv.symm r).1`. -/
def blkStack (n : Fin M → ℕ → ℕ) (N : ℕ) : Fin (∑ i, n i N) → Fin M :=
  fun r => (finSigmaFinEquiv.symm r).1

omit [NeZero M] in
theorem SigmaHalf_eq_diagonal (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    m.SigmaHalf w N = Matrix.diagonal (tauOf w (blkStack n N)) := rfl

/-- **`W₀' = Wsig τ d B`** (the open point of `notes/archive/agent_reports/h5_stein.md`): on the
block `B` of `exists_block_hasLaw_het`, the noise Gram matrix of the column split is the
heteroscedastic Wishart matrix of `RMT/Het/Stein.lean` with `τ_r = w_{blk r}`. -/
theorem W0het_eq_Wsig (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {p : ℕ} (B : Matrix (Fin (∑ i, n i N)) (Fin p) ℝ)
    (hB : m.EperpHet N ω * (m.EperpHet N ω)ᵀ = ((d N : ℝ))⁻¹ • (B * Bᵀ)) :
    m.W0het w N ω = Wsig (tauOf w (blkStack n N)) (d N) B := by
  rw [m.W0het_eq_of_block w N ω B hB]
  rfl

omit [NeZero M] in
/-- The block sets of `blkStack` have the table sizes. -/
theorem card_blockSet_blkStack (N : ℕ) (i : Fin M) :
    (blockSet (blkStack n N) i).card = n i N := by
  classical
  rw [blockSet, Finset.card_filter]
  have h := Fintype.sum_equiv finSigmaFinEquiv
    (fun σ : Σ i', Fin (n i' N) => if σ.1 = i then 1 else 0)
    (fun r => if blkStack n N r = i then 1 else 0) (fun σ => by simp [blkStack])
  rw [← h, Fintype.sum_sigma]
  simp [apply_ite Finset.card]

end Model


/-! ### The derivative version: block traces of `G_σ²` (item R1b pattern, consumed by H7) -/

section Deriv

variable {p q d : ℕ} {z : ℂ}

/-- `|J|⁻¹ ∑_{j ∈ J} (G_σ²)_{jj}`, the block average of the squared resolvent. -/
noncomputable def gsig2Avg (τ : Fin p → ℝ) (d : ℕ) (z : ℂ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) : ℂ :=
  (J.card : ℂ)⁻¹ * ∑ j ∈ J, (Gsig τ d z B * Gsig τ d z B) j j

theorem norm_gsig2Avg_le (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) : ‖gsig2Avg τ d z J B‖ ≤ 1 / z.im ^ 2 := by
  rw [gsig2Avg]
  by_cases hJ : J.card = 0
  · rw [hJ]; simp only [CharP.cast_eq_zero, inv_zero, zero_mul, norm_zero, one_div, inv_nonneg]
    positivity
  have hc : (0 : ℝ) < J.card := by exact_mod_cast Nat.pos_of_ne_zero hJ
  rw [norm_mul, norm_inv, Complex.norm_natCast]
  calc ((J.card : ℝ))⁻¹ * ‖∑ j ∈ J, (Gsig τ d z B * Gsig τ d z B) j j‖
      ≤ ((J.card : ℝ))⁻¹ * ∑ j ∈ J, (1 / z.im ^ 2) :=
        mul_le_mul_of_nonneg_left ((norm_sum_le _ _).trans
          (Finset.sum_le_sum fun j _ => norm_Gsig_sq_diag_le hz τ d B j)) (by positivity)
    _ = 1 / z.im ^ 2 := by
        rw [Finset.sum_const, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hc.ne', one_mul]

/-- `ζ ↦ G(ζ)_{jj}` is holomorphic on `ℂ⁺` with derivative `(G(ζ)²)_{jj}`. -/
theorem hasDerivAt_resolvC_diag {n : ℕ} {W : Matrix (Fin n) (Fin n) ℝ} (hW : W.IsHermitian)
    (hz : 0 < z.im) (j : Fin n) :
    HasDerivAt (fun ζ => R4C.resolvC W ζ j j) ((R4C.resolvC W z * R4C.resolvC W z) j j) z := by
  have hopen : IsOpen {ζ : ℂ | 0 < ζ.im} := isOpen_lt continuous_const Complex.continuous_im
  have hentry : ∀ (f : Fin n → ℂ),
      (R4C.cmat (eigU hW) * Matrix.diagonal f * (R4C.cmat (eigU hW))ᵀ) j j
        = ∑ a, ((eigU hW j a : ℝ) : ℂ) * f a * ((eigU hW j a : ℝ) : ℂ) := by
    intro f
    rw [Matrix.mul_apply]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Matrix.mul_diagonal, Matrix.transpose_apply]
    rfl
  have hkey : HasDerivAt (fun ζ : ℂ => ∑ a, ((eigU hW j a : ℝ) : ℂ)
      * ((hW.eigenvalues a : ℂ) - ζ)⁻¹ * ((eigU hW j a : ℝ) : ℂ))
      (∑ a, ((eigU hW j a : ℝ) : ℂ) * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2
        * ((eigU hW j a : ℝ) : ℂ)) z := by
    refine HasDerivAt.fun_sum fun a _ => ?_
    have h1 : HasDerivAt (fun ζ : ℂ => (hW.eigenvalues a : ℂ) - ζ) (-1) z :=
      HasDerivAt.const_sub ((hW.eigenvalues a : ℂ)) (hasDerivAt_id z)
    have h2 := h1.inv (R4C.eigenvalue_sub_ne_zero hW hz.ne' a)
    have heq : -(-1 : ℂ) / ((hW.eigenvalues a : ℂ) - z) ^ 2
        = (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 := by
      rw [inv_pow, neg_neg, one_div]
    rw [heq] at h2
    have := (h2.const_mul ((eigU hW j a : ℝ) : ℂ)).mul_const ((eigU hW j a : ℝ) : ℂ)
    exact this
  have hval : (R4C.resolvC W z * R4C.resolvC W z) j j
      = ∑ a, ((eigU hW j a : ℝ) : ℂ) * (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2
        * ((eigU hW j a : ℝ) : ℂ) := by
    rw [R4C.resolvC_mul_resolvC_eq_conj hW hz.ne', hentry]
  rw [hval]
  refine hkey.congr_of_eventuallyEq ?_
  filter_upwards [hopen.mem_nhds hz] with ζ hζ
  show R4C.resolvC W ζ j j = _
  rw [R4C.resolvC_eq_conj hW (ne_of_gt hζ), hentry]

theorem hasDerivAt_gsigAvg (hz : 0 < z.im) (τ : Fin p → ℝ) (d : ℕ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) :
    HasDerivAt (fun ζ => gsigAvg τ d ζ J B) (gsig2Avg τ d z J B) z := by
  unfold gsigAvg gsig2Avg
  refine HasDerivAt.const_mul _ (HasDerivAt.fun_sum fun j _ => ?_)
  exact hasDerivAt_resolvC_diag (isHermitian_Wsig τ d B) hz j

end Deriv

section DerivRoot

variable {c w : Fin M → ℝ} {z : ℂ}

/-- A root in `ℂ⁺` has `|σ| ≤ 1 / Im z`, from the imaginary-part identity. -/
theorem IsRoot.norm_le {σ : ℂ} (hc : ∀ i, 0 ≤ c i) (hz : 0 < z.im) (h : IsRoot c w z σ) :
    ‖σ‖ ≤ 1 / z.im := by
  have hid := im_zfunC c w σ
  rw [h.2] at hid
  have hA := Aabs_nonneg (w := w) hc σ
  have hσn : 0 < ‖σ‖ := norm_pos_iff.mpr (ne_zero_of_im_pos h.1)
  have hsn : σ.im ≤ ‖σ‖ := le_trans (le_abs_self _) (Complex.abs_im_le_norm σ)
  have h1 : z.im ≤ σ.im * (1 / ‖σ‖ ^ 2) := by
    rw [hid]
    exact mul_le_mul_of_nonneg_left (by linarith) h.1.le
  have h2 : σ.im * (1 / ‖σ‖ ^ 2) ≤ ‖σ‖ * (1 / ‖σ‖ ^ 2) :=
    mul_le_mul_of_nonneg_right hsn (by positivity)
  have h3 : ‖σ‖ * (1 / ‖σ‖ ^ 2) = 1 / ‖σ‖ := by field_simp
  rw [le_div_iff₀ hz]
  have : z.im ≤ 1 / ‖σ‖ := by linarith
  rwa [le_div_iff₀ hσn, mul_comm] at this

/-- **Holomorphy of the global root.** At a `z ∈ ℂ⁺` where the root exists, `sGlob` is
strictly differentiable with derivative `1 / zfunC'(sGlob z)`: it agrees near `z` with the
local inverse of `zfunC` at `sGlob z`, by uniqueness of the root. -/
theorem hasStrictDerivAt_sGlob (hc : ∀ i, 0 ≤ c i) (hz : 0 < z.im)
    (hroot : IsRoot c w z (sGlob c w z)) :
    HasStrictDerivAt (sGlob c w) (zfunDerivC c w (sGlob c w z))⁻¹ z := by
  set σ := sGlob c w z with hσ
  have hσ0 := ne_zero_of_im_pos hroot.1
  have hσne := fun i => one_add_mul_ne_zero_of_im_pos hroot.1 (w i)
  have hf := hasStrictDerivAt_zfunC c w hσ0 hσne
  have hf' := zfunDerivC_ne_zero_of_root hc hz hroot.1 hroot.2
  set Linv := HasStrictDerivAt.localInverse (zfunC c w) (zfunDerivC c w σ) σ hf hf' with hL
  have hLσ : Linv z = σ := by
    rw [hL, ← hroot.2]
    exact HasStrictFDerivAt.localInverse_apply_image _
  have hder : HasStrictDerivAt Linv (zfunDerivC c w σ)⁻¹ z := by
    rw [hL, ← hroot.2]
    exact HasStrictDerivAt.to_localInverse _ _
  have hright : ∀ᶠ ζ in 𝓝 z, zfunC c w (Linv ζ) = ζ := by
    have := HasStrictDerivAt.eventually_right_inverse hf hf'
    rwa [hroot.2] at this
  have hcont : ContinuousAt Linv z := hder.hasDerivAt.continuousAt
  have him : ∀ᶠ ζ in 𝓝 z, 0 < (Linv ζ).im := by
    have := (Complex.continuous_im.continuousAt.comp hcont)
    rw [ContinuousAt, Function.comp_apply, hLσ] at this
    exact this.eventually (lt_mem_nhds hroot.1)
  have hzim : ∀ᶠ ζ in 𝓝 z, 0 < ζ.im :=
    (isOpen_lt continuous_const Complex.continuous_im).mem_nhds hz
  refine hder.congr_of_eventuallyEq ?_
  filter_upwards [hright, him, hzim] with ζ h1 h2 h3
  exact (sGlob_eq hc h3 ⟨h2, h1⟩).symm

/-- `d/dz [-1/(z (1 + w² s(z)))]` with `s' = 1/zfunC'(s)`. -/
noncomputable def gCDeriv (w : Fin M → ℝ) (i : Fin M) (z s s' : ℂ) : ℂ :=
  (1 + ((w i : ℝ) : ℂ) ^ 2 * s + z * ((w i : ℝ) : ℂ) ^ 2 * s')
    / (z * (1 + ((w i : ℝ) : ℂ) ^ 2 * s)) ^ 2

theorem gCDeriv_ofReal (c w : Fin M → ℝ) (i : Fin M) (x : ℝ) :
    gCDeriv w i (x : ℂ) ((sPhys c w x : ℝ) : ℂ) ((sPhysDeriv c w x : ℝ) : ℂ)
      = (ghetDeriv c w i x : ℂ) := by
  simp only [gCDeriv, ghetDeriv]
  push_cast
  rfl

theorem hasDerivAt_gC_sGlob (hc : ∀ i, 0 ≤ c i) (hz : 0 < z.im)
    (hroot : IsRoot c w z (sGlob c w z)) (i : Fin M) :
    HasDerivAt (fun ζ => gC w i ζ (sGlob c w ζ))
      (gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹) z := by
  have hs := (hasStrictDerivAt_sGlob hc hz hroot).hasDerivAt
  have hden : z * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w z) ≠ 0 :=
    mul_ne_zero (ne_zero_of_im_pos hz) (one_add_mul_ne_zero_of_im_pos hroot.1 _)
  have h1 : HasDerivAt (fun ζ => ζ * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w ζ))
      (1 * (1 + ((w i : ℝ) : ℂ) ^ 2 * sGlob c w z)
        + z * (((w i : ℝ) : ℂ) ^ 2 * (zfunDerivC c w (sGlob c w z))⁻¹)) z :=
    (hasDerivAt_id z).mul ((hs.const_mul _).const_add 1)
  have h2 := (h1.inv hden).neg
  refine h2.congr_deriv ?_
  rw [gCDeriv]
  field_simp

/-- Boundary values of the derivative: along the ray, `gCDeriv → ghetDeriv c w i x`. -/
theorem tendsto_gCDeriv_sGlob (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {x : ℝ}
    (hx : bHet c w < x) (i : Fin M) :
    Tendsto (fun η : ℝ => gCDeriv w i ((x : ℂ) + (η : ℂ) * Complex.I)
      (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I))
      (zfunDerivC c w (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I)))⁻¹) (𝓝[>] 0)
      (𝓝 ((ghetDeriv c w i x : ℝ) : ℂ)) := by
  rw [← gCDeriv_ofReal c w i x]
  have hx0 : (x : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (lt_trans (bHet_pos hc hw) hx).ne'
  have hden : (x : ℂ) * (1 + ((w i : ℝ) : ℂ) ^ 2 * ((sPhys c w x : ℝ) : ℂ)) ≠ 0 :=
    mul_ne_zero hx0 (one_add_mul_sPhys_ne_zero_C hc hw hx i)
  have hS := tendsto_sGlob hc hw hx
  have hpath := MP.tendsto_ofReal_add_mul_I x
  have hD : Tendsto (fun η : ℝ => (zfunDerivC c w (sGlob c w ((x : ℂ) + (η : ℂ) * Complex.I)))⁻¹)
      (𝓝[>] 0) (𝓝 (((sPhysDeriv c w x : ℝ) : ℂ))) := by
    have hcont : ContinuousAt (zfunDerivC c w) ((sPhys c w x : ℝ) : ℂ) := by
      unfold zfunDerivC
      have h0 := sPhys_ne_zero_C hc hw hx
      have h1 := one_add_mul_sPhys_ne_zero_C hc hw hx
      refine ContinuousAt.sub (ContinuousAt.inv₀ (by fun_prop) (pow_ne_zero 2 h0))
        (tendsto_finsetSum _ fun j _ => ContinuousAt.mul continuousAt_const
          (ContinuousAt.inv₀ (by fun_prop) (pow_ne_zero 2 (h1 j))))
    have h := (hcont.tendsto.comp hS).inv₀ (by
      rw [zfunDerivC_ofReal]
      exact zfunDeriv_sPhys_ne_zero_C hc hw hx)
    rw [zfunDerivC_ofReal, ← Complex.ofReal_inv] at h
    exact h
  unfold gCDeriv
  refine Tendsto.div ?_ ?_ (pow_ne_zero 2 hden)
  · exact (tendsto_const_nhds.add (tendsto_const_nhds.mul hS)).add
      ((hpath.mul tendsto_const_nhds).mul hD)
  · exact (hpath.mul (tendsto_const_nhds.add (tendsto_const_nhds.mul hS))).pow 2

end DerivRoot

section DerivMain

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN qN dN : ℕ → ℕ} {blk : (N : ℕ) → Fin (pN N) → Fin M} {c w : Fin M → ℝ} {z : ℂ}

variable (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (hz : 0 < z.im)
  (hd : Tendsto dN atTop atTop) (hq : Tendsto (fun N => (qN N : ℝ) / dN N) atTop (𝓝 1))
  (hn : ∀ i, Tendsto (fun N => cN blk dN N i) atTop (𝓝 (c i)))
  (B : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (qN N)) ℝ)
  (hB : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (qN N)) (μ N))

/-- Equicontinuity in `ζ` of the block average on the disc of radius `3η/4`: the derivative
`gsig2Avg` is bounded by `16/η²` there. -/
theorem norm_gsigAvg_sub_le {p q d : ℕ} (hz : 0 < z.im) (τ : Fin p → ℝ) (J : Finset (Fin p))
    (B : Matrix (Fin p) (Fin q) ℝ) {ζ ζ' : ℂ} (hζ : ζ ∈ Metric.closedBall z (3 * z.im / 4))
    (hζ' : ζ' ∈ Metric.closedBall z (3 * z.im / 4)) :
    ‖gsigAvg τ d ζ J B - gsigAvg τ d ζ' J B‖ ≤ (16 / z.im ^ 2) * ‖ζ - ζ'‖ := by
  set F : ℂ → ℂ := fun v => gsigAvg τ d v J B with hF
  set F' : ℂ → ℂ := fun v => gsig2Avg τ d v J B with hF'
  have him : ∀ y ∈ Metric.closedBall z (3 * z.im / 4), z.im / 4 ≤ y.im := by
    intro y hy
    have := R1.im_ge_of_mem_closedBall (z := z) hy
    linarith
  have hbound : ∀ y ∈ Metric.closedBall z (3 * z.im / 4), ‖F' y‖ ≤ 16 / z.im ^ 2 := by
    intro y hy
    have h1 : z.im / 4 ≤ y.im := him y hy
    have hypos : 0 < y.im := by linarith
    refine (norm_gsig2Avg_le hypos τ d J B).trans ?_
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [sq_nonneg (y.im - z.im / 4), sq_nonneg z.im]
  have hder : ∀ y ∈ Metric.closedBall z (3 * z.im / 4),
      HasDerivWithinAt F (F' y) (Metric.closedBall z (3 * z.im / 4)) y := by
    intro y hy
    have h1 : z.im / 4 ≤ y.im := him y hy
    have hypos : 0 < y.im := by linarith
    exact (hasDerivAt_gsigAvg hypos τ d J B).hasDerivWithinAt
  have := (convex_closedBall z (3 * z.im / 4)).norm_image_sub_le_of_norm_hasDerivWithin_le
    hder hbound hζ' hζ
  exact this

include hc hw hz hd hq hn hB in
/-- **The derivative version of the block trace law** (the R1b pattern): the block average of
`G_σ(z)²` converges to `d/dz [-1/(z (1 + w_i² sGlob z))]`. Proof: Cauchy's estimate on the
circle of radius `η/2`, a finite net of the circle, equicontinuity of the random part
(`norm_gsigAvg_sub_le`) and uniform continuity of the deterministic limit on the compact
circle, then `tendstoInProb_gsigAvg` at each net point with a union bound. -/
theorem tendstoInProb_gsig2Avg (i : Fin M) :
    TendstoInProb μ (fun N ω => ‖gsig2Avg (tauOf w (blk N)) (dN N) z (blockSet (blk N) i) (B N ω)
      - gCDeriv w i z (sGlob c w z) (zfunDerivC c w (sGlob c w z))⁻¹‖) 0 := by
  intro ε hε
  have hc' : ∀ i, 0 ≤ c i := fun i => (hc i).le
  have hroots : ∀ ζ : ℂ, 0 < ζ.im → IsRoot c w ζ (sGlob c w ζ) :=
    fun ζ hζ => isRoot_sGlob (exists_isRoot hc hw hζ hd hq hn)
  set Lf : ℂ → ℂ := fun ζ => gC w i ζ (sGlob c w ζ) with hLf
  set Lf' : ℂ → ℂ := fun ζ => gCDeriv w i ζ (sGlob c w ζ) (zfunDerivC c w (sGlob c w ζ))⁻¹
    with hLf'
  have hLfder : ∀ ζ : ℂ, 0 < ζ.im → HasDerivAt Lf (Lf' ζ) ζ :=
    fun ζ hζ => hasDerivAt_gC_sGlob hc' hζ (hroots ζ hζ) i
  set r : ℝ := z.im / 2 with hrdef
  have hrpos : (0 : ℝ) < r := by rw [hrdef]; linarith
  set K := Metric.closedBall z (3 * z.im / 4) with hK
  have himK : ∀ y ∈ K, z.im / 4 ≤ y.im := by
    intro y hy
    have := R1.im_ge_of_mem_closedBall (z := z) hy
    linarith
  have hsphereK : Metric.sphere z r ⊆ K := by
    intro y hy
    have h1 : dist y z = r := hy
    rw [hK, Metric.mem_closedBall, h1, hrdef]
    linarith
  set Csup : ℝ := ε * r / 2 with hCdef
  have hCpos : (0 : ℝ) < Csup := by rw [hCdef]; positivity
  -- uniform continuity of the limit on the circle
  have hLfcont : ContinuousOn Lf (Metric.sphere z r) := fun y hy =>
    (hLfder y (by linarith [himK y (hsphereK hy)])).continuousAt.continuousWithinAt
  obtain ⟨δ, hδpos, hδ⟩ := Metric.uniformContinuousOn_iff.mp
    ((isCompact_sphere z r).uniformContinuousOn_of_continuous hLfcont) (Csup / 4)
    (by positivity)
  set Lz : ℝ := 16 / z.im ^ 2 with hLzdef
  have hLzpos : (0 : ℝ) < Lz := by rw [hLzdef]; positivity
  set ρ : ℝ := min δ (Csup / (4 * Lz)) with hrhodef
  have hrhopos : (0 : ℝ) < ρ := by rw [hrhodef]; exact lt_min hδpos (by positivity)
  obtain ⟨b, hbsub, hbfin, hbcov⟩ :=
    (isCompact_sphere z r).elim_finite_subcover_image
      (b := Metric.sphere z r) (c := fun ζ : ℂ => Metric.ball ζ ρ)
      (fun ζ _ => Metric.isOpen_ball)
      (fun ζ hζ => Set.mem_biUnion hζ (Metric.mem_ball_self hrhopos))
  have hyim : ∀ y ∈ b, 0 < y.im := by
    intro y hy
    have := himK y (hsphereK (hbsub hy))
    linarith
  -- the block trace law at each net point
  have hpt : ∀ y ∈ hbfin.toFinset,
      Tendsto (fun N => μ N {ω | Csup / 4 ≤ ‖gsigAvg (tauOf w (blk N)) (dN N) y (blockSet (blk N) i)
        (B N ω) - Lf y‖}) atTop (𝓝 0) := by
    intro y hy
    have hyb : y ∈ b := hbfin.mem_toFinset.mp hy
    have h := tendstoInProb_gsigAvg hc hw (hyim y hyb) hd hq hn B hB i (Csup / 4) (by positivity)
    have hset : ∀ N, {ω : Ω N | Csup / 4 ≤ ‖gsigAvg (tauOf w (blk N)) (dN N) y (blockSet (blk N) i)
        (B N ω) - Lf y‖}
        = {ω | Csup / 4 ≤ |‖gsigAvg (tauOf w (blk N)) (dN N) y (blockSet (blk N) i) (B N ω)
            - gC w i y (sGlob c w y)‖ - 0|} := by
      intro N
      ext ω
      simp [hLf, abs_of_nonneg (norm_nonneg _)]
    simpa [hset] using h
  have hsum : Tendsto (fun N => ∑ y ∈ hbfin.toFinset,
      μ N {ω | Csup / 4 ≤ ‖gsigAvg (tauOf w (blk N)) (dN N) y (blockSet (blk N) i) (B N ω)
        - Lf y‖}) atTop (𝓝 0) := by
    have h := tendsto_finsetSum hbfin.toFinset hpt
    simpa using h
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hsum
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  refine Filter.Eventually.of_forall fun N => ?_
  set τ := tauOf w (blk N) with hτ
  set J := blockSet (blk N) i with hJ
  have hsubset : {ω | ε ≤ |‖gsig2Avg τ (dN N) z J (B N ω) - Lf' z‖ - 0|}
      ⊆ ⋃ y ∈ hbfin.toFinset,
          {ω | Csup / 4 ≤ ‖gsigAvg τ (dN N) y J (B N ω) - Lf y‖} := by
    intro ω hω
    by_contra hcon
    simp only [Set.mem_iUnion, Set.mem_ofPred_eq, not_exists, not_le] at hcon
    -- the sup on the circle
    have hC : ∀ ζ ∈ Metric.sphere z r, ‖gsigAvg τ (dN N) ζ J (B N ω) - Lf ζ‖ ≤ Csup := by
      intro ζ hζ
      obtain ⟨y, hyb, hyball'⟩ : ∃ y ∈ b, ζ ∈ Metric.ball y ρ := by
        have := hbcov hζ
        simpa using this
      have hlt : ‖gsigAvg τ (dN N) y J (B N ω) - Lf y‖ < Csup / 4 :=
        hcon y (hbfin.mem_toFinset.mpr hyb)
      have hdist : dist ζ y < ρ := hyball'
      have hζK : ζ ∈ K := hsphereK hζ
      have hyK : y ∈ K := hsphereK (hbsub hyb)
      have hrand : ‖gsigAvg τ (dN N) ζ J (B N ω) - gsigAvg τ (dN N) y J (B N ω)‖ ≤ Csup / 4 := by
        refine (norm_gsigAvg_sub_le hz τ J (B N ω) hζK hyK).trans ?_
        have h1 : ‖ζ - y‖ ≤ Csup / (4 * Lz) := by
          rw [← Complex.dist_eq]
          exact hdist.le.trans (min_le_right _ _)
        calc Lz * ‖ζ - y‖ ≤ Lz * (Csup / (4 * Lz)) := mul_le_mul_of_nonneg_left h1 hLzpos.le
          _ = Csup / 4 := by field_simp
      have hdet : ‖Lf ζ - Lf y‖ ≤ Csup / 4 := by
        have h := hδ ζ hζ y (hbsub hyb) (hdist.trans_le (min_le_left _ _))
        rw [Complex.dist_eq] at h
        exact h.le
      have htri : ‖gsigAvg τ (dN N) ζ J (B N ω) - Lf ζ‖
          ≤ ‖gsigAvg τ (dN N) ζ J (B N ω) - gsigAvg τ (dN N) y J (B N ω)‖
            + ‖gsigAvg τ (dN N) y J (B N ω) - Lf y‖ + ‖Lf y - Lf ζ‖ := by
        have := norm_add₃_le (a := gsigAvg τ (dN N) ζ J (B N ω) - gsigAvg τ (dN N) y J (B N ω))
          (b := gsigAvg τ (dN N) y J (B N ω) - Lf y) (c := Lf y - Lf ζ)
        simpa using this
      rw [norm_sub_rev] at hdet
      linarith
    -- Cauchy's estimate
    have hball : ∀ ζ ∈ Metric.closedBall z r, 0 < ζ.im := by
      intro ζ hζ
      have := R1.im_ge_of_mem_closedBall (z := z) hζ
      rw [hrdef] at this
      linarith
    have hDder : ∀ ζ ∈ Metric.closedBall z r,
        HasDerivAt (fun v => gsigAvg τ (dN N) v J (B N ω) - Lf v)
          (gsig2Avg τ (dN N) ζ J (B N ω) - Lf' ζ) ζ :=
      fun ζ hζ => (hasDerivAt_gsigAvg (hball ζ hζ) τ (dN N) J (B N ω)).sub (hLfder ζ (hball ζ hζ))
    have hdiffOn : DifferentiableOn ℂ (fun v => gsigAvg τ (dN N) v J (B N ω) - Lf v)
        (Metric.closedBall z r) := fun ζ hζ =>
      (hDder ζ hζ).differentiableAt.differentiableWithinAt
    have hdc : DiffContOnCl ℂ (fun v => gsigAvg τ (dN N) v J (B N ω) - Lf v) (Metric.ball z r) := by
      constructor
      · exact hdiffOn.mono Metric.ball_subset_closedBall
      · rw [closure_ball z (ne_of_gt hrpos)]
        exact hdiffOn.continuousOn
    have hderiv : deriv (fun v => gsigAvg τ (dN N) v J (B N ω) - Lf v) z
        = gsig2Avg τ (dN N) z J (B N ω) - Lf' z :=
      (hDder z (Metric.mem_closedBall_self hrpos.le)).deriv
    have hcauchy := Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hrpos hdc hC
    rw [hderiv] at hcauchy
    have hmem : ε ≤ ‖gsig2Avg τ (dN N) z J (B N ω) - Lf' z‖ := by
      have hω' : ε ≤ |‖gsig2Avg τ (dN N) z J (B N ω) - Lf' z‖ - 0| := hω
      rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
    have hval : Csup / r = ε / 2 := by
      rw [hCdef]
      field_simp
    rw [hval] at hcauchy
    linarith
  exact (measure_mono hsubset).trans (measure_biUnion_finset_le _ _)

end DerivMain

end HetR1
end StackedSVD
