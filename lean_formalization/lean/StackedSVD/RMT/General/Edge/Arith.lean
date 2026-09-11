/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Defs

/-! # Stage 3, unit K (first half): the arithmetic of the moment order

The facts about `momOrder`, `fkFactor` and `markovConst` that the Gaussian bound (unit G),
the bundled moment bound (unit M) and the Markov step (K8) read: `1 ≤ k_N` and
`8 k_N ≤ d_N` for `N` large, the Furedi-Komlos factor at most `1/2` for `N` large, and the
Markov bound `4 d ((bulkEdge c + ε)/(bulkEdge c + 2ε))^(k_N) → 0` at `C = markovConst c ε`
(`notes/stage3_edge.md`, route C; `notes/STAGE3_CAMPAIGN.md`). -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace Edge

/-- The real cast of `d_N` tends to infinity. -/
private theorem tendsto_cast_atTop {d : ℕ → ℕ} (hdtop : Tendsto d atTop atTop) :
    Tendsto (fun N => (d N : ℝ)) atTop atTop :=
  tendsto_natCast_atTop_atTop.comp hdtop

/-- `bulkEdge c = (1 + √c)²` is positive. -/
private theorem bulkEdge_pos' (c : ℝ) : 0 < bulkEdge c := by
  rw [bulkEdge]
  positivity

/-- `T² = D (D^(-a/6))¹²` at `T = D^(1/2-a)`: the form K2 reads, with `D^(-a/6) → 0`. -/
private theorem truncLevel_sq_eq {a : ℝ} {D : ℕ} (hD : 0 < D) :
    truncLevel a D ^ 2 = (D : ℝ) * ((D : ℝ) ^ (-(a / 6))) ^ 12 := by
  have hD' : (0 : ℝ) < D := by exact_mod_cast hD
  have e1 : ((D : ℝ) ^ ((1 : ℝ) / 2 - a)) ^ (2 : ℕ) = (D : ℝ) ^ ((1 : ℝ) - 2 * a) := by
    rw [← Real.rpow_natCast ((D : ℝ) ^ ((1 : ℝ) / 2 - a)) 2, ← Real.rpow_mul hD'.le]
    congr 1
    push_cast
    ring
  have e2 : ((D : ℝ) ^ (-(a / 6))) ^ (12 : ℕ) = (D : ℝ) ^ (-(2 * a)) := by
    rw [← Real.rpow_natCast ((D : ℝ) ^ (-(a / 6))) 12, ← Real.rpow_mul hD'.le]
    congr 1
    push_cast
    ring
  have e3 : (D : ℝ) * (D : ℝ) ^ (-(2 * a)) = (D : ℝ) ^ ((1 : ℝ) - 2 * a) := by
    rw [show (1 : ℝ) - 2 * a = 1 + -(2 * a) by ring, Real.rpow_add hD', Real.rpow_one]
  rw [truncLevel, e1, e2]
  exact e3.symm

/-- K0b. `1 ≤ k_N` for `N` large. -/
theorem eventually_one_le_momOrder {d : ℕ → ℕ} (hdtop : Tendsto d atTop atTop) {C : ℝ}
    (hC : 0 < C) : ∀ᶠ N in atTop, 1 ≤ momOrder C (d N) := by
  filter_upwards [hdtop.eventually_ge_atTop 2] with N hN
  have h1 : (1 : ℝ) < (d N : ℝ) := by exact_mod_cast (by omega : 1 < d N)
  have hpos : 0 < C * Real.log (d N) := mul_pos hC (Real.log_pos h1)
  have := Nat.ceil_pos.mpr hpos
  rw [momOrder]
  omega

/-- K1. `8 k_N ≤ d_N` for `N` large: `k_N = ⌈C log d⌉` grows like `log d`. -/
theorem eventually_momOrder_le {d : ℕ → ℕ} (hdtop : Tendsto d atTop atTop) {C : ℝ}
    (hC : 0 < C) : ∀ᶠ N in atTop, 8 * momOrder C (d N) ≤ d N := by
  have hlo : ∀ᶠ x : ℝ in atTop, ‖Real.log x‖ ≤ (1 / (16 * C)) * ‖id x‖ :=
    Asymptotics.IsLittleO.def Real.isLittleO_log_id_atTop (by positivity)
  filter_upwards [(tendsto_cast_atTop hdtop).eventually hlo,
    hdtop.eventually_ge_atTop 16] with N hN hd16
  have hd0 : (0 : ℝ) < (d N : ℝ) := by
    have : (16 : ℝ) ≤ (d N : ℝ) := by exact_mod_cast hd16
    linarith
  have hlogd : (0 : ℝ) ≤ Real.log (d N) := Real.log_natCast_nonneg _
  have hd16' : (16 : ℝ) ≤ (d N : ℝ) := by exact_mod_cast hd16
  have hbound : Real.log (d N) ≤ (1 / (16 * C)) * (d N : ℝ) := by
    simpa [Real.norm_eq_abs, abs_of_nonneg hlogd, abs_of_nonneg hd0.le] using hN
  have hceil : ((momOrder C (d N) : ℕ) : ℝ) < C * Real.log (d N) + 1 := by
    rw [momOrder]
    exact Nat.ceil_lt_add_one (by positivity)
  have hkey : (8 : ℝ) * ((momOrder C (d N) : ℕ) : ℝ) ≤ (d N : ℝ) := by
    have h1 : C * Real.log (d N) ≤ (d N : ℝ) / 16 := by
      have := mul_le_mul_of_nonneg_left hbound hC.le
      rw [show C * ((1 / (16 * C)) * (d N : ℝ)) = (d N : ℝ) / 16 by field_simp] at this
      exact this
    nlinarith [hceil, h1, hd16']
  exact_mod_cast hkey

/-- K2. `f_N ≤ 1/2` for `N` large at `K = 2 d^(1/2-a)` and `k = ⌈C log d⌉`:
`f_N ≍ 4 (2 k_N)^12 d^(-2a) / min(c,1)`. Every `a > 0` works. -/
theorem eventually_fkFactor_le {n d : ℕ → ℕ} {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) {C : ℝ} (hC : 0 < C) :
    ∀ᶠ N in atTop,
      fkFactor (2 * truncLevel a (d N)) (n N) (d N) (momOrder C (d N)) ≤ 1 / 2 := by
  have ha3 : (0 : ℝ) < a / 6 := by positivity
  set m : ℝ := min (c / 2) 1 with hmdef
  have hm0 : 0 < m := lt_min (by linarith) one_pos
  -- the two vanishing scalars
  have hw : Tendsto (fun N => ((d N : ℝ)) ^ (-(a / 6))) atTop (𝓝 0) :=
    (tendsto_rpow_neg_atTop ha3).comp (tendsto_cast_atTop hdtop)
  have hlw : Tendsto (fun N => Real.log (d N) * ((d N : ℝ)) ^ (-(a / 6))) atTop (𝓝 0) := by
    have h0 : Tendsto (fun x : ℝ => Real.log x / x ^ (a / 6)) atTop (𝓝 0) :=
      (isLittleO_log_rpow_atTop ha3).tendsto_div_nhds_zero
    have h1 : Tendsto (fun x : ℝ => Real.log x * x ^ (-(a / 6))) atTop (𝓝 0) := by
      refine h0.congr' ?_
      filter_upwards [eventually_gt_atTop (0 : ℝ)] with x hx
      rw [Real.rpow_neg hx.le, div_eq_mul_inv]
    exact h1.comp (tendsto_cast_atTop hdtop)
  have hsum : Tendsto (fun N => 2 * C * (Real.log (d N) * ((d N : ℝ)) ^ (-(a / 6)))
      + 2 * ((d N : ℝ)) ^ (-(a / 6))) atTop (𝓝 0) := by
    simpa using (hlw.const_mul (2 * C)).add (hw.const_mul (2 : ℝ))
  have hG : Tendsto (fun N => 4 / m * (2 * C * (Real.log (d N) * ((d N : ℝ)) ^ (-(a / 6)))
      + 2 * ((d N : ℝ)) ^ (-(a / 6))) ^ 12) atTop (𝓝 0) := by
    simpa using (hsum.pow 12).const_mul (4 / m)
  have hGle := hG.eventually_le_const (show (0 : ℝ) < 1 / 2 by norm_num)
  have hcn : ∀ᶠ N in atTop, c / 2 ≤ (n N : ℝ) / (d N : ℝ) :=
    hcN.eventually_const_le (by linarith)
  filter_upwards [hGle, hcn, hdtop.eventually_ge_atTop 1] with N hG1 hcn1 hd1
  have hd0 : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd1
  have hlogd : (0 : ℝ) ≤ Real.log (d N) := Real.log_natCast_nonneg _
  set w : ℝ := ((d N : ℝ)) ^ (-(a / 6)) with hwdef
  have hw0 : 0 < w := by rw [hwdef]; exact Real.rpow_pos_of_pos hd0 _
  -- `min(n, d)` is at least `m d`
  have hminle : m * (d N : ℝ) ≤ min (n N : ℝ) (d N : ℝ) := by
    refine le_min ?_ ?_
    · have hcd : c / 2 * (d N : ℝ) ≤ (n N : ℝ) := by
        have h := mul_le_mul_of_nonneg_right hcn1 hd0.le
        rwa [div_mul_cancel₀ _ hd0.ne'] at h
      have : m * (d N : ℝ) ≤ c / 2 * (d N : ℝ) :=
        mul_le_mul_of_nonneg_right (min_le_left _ _) hd0.le
      linarith
    · have : m * (d N : ℝ) ≤ 1 * (d N : ℝ) :=
        mul_le_mul_of_nonneg_right (min_le_right _ _) hd0.le
      linarith
  have hmd0 : (0 : ℝ) < m * (d N : ℝ) := mul_pos hm0 hd0
  have hk : ((momOrder C (d N) : ℕ) : ℝ) ≤ C * Real.log (d N) + 1 := by
    rw [momOrder]
    exact (Nat.ceil_lt_add_one (by positivity)).le
  have hknn : (0 : ℝ) ≤ ((momOrder C (d N) : ℕ) : ℝ) := Nat.cast_nonneg _
  calc fkFactor (2 * truncLevel a (d N)) (n N) (d N) (momOrder C (d N))
      = 4 * ((d N : ℝ) * w ^ 12) * (2 * ((momOrder C (d N) : ℕ) : ℝ)) ^ 12
          / min (n N : ℝ) (d N : ℝ) := by
        rw [fkFactor, mul_pow, truncLevel_sq_eq (by omega : 0 < d N), ← hwdef]
        norm_num
    _ ≤ 4 * ((d N : ℝ) * w ^ 12) * (2 * ((momOrder C (d N) : ℕ) : ℝ)) ^ 12
          / (m * (d N : ℝ)) := by
        refine div_le_div_of_nonneg_left ?_ hmd0 hminle
        positivity
    _ = 4 / m * (w * (2 * ((momOrder C (d N) : ℕ) : ℝ))) ^ 12 := by
        field_simp
    _ ≤ 4 / m * (w * (2 * (C * Real.log (d N) + 1))) ^ 12 := by
        have hle : w * (2 * ((momOrder C (d N) : ℕ) : ℝ)) ≤ w * (2 * (C * Real.log (d N) + 1)) := by
          have : (2 : ℝ) * ((momOrder C (d N) : ℕ) : ℝ) ≤ 2 * (C * Real.log (d N) + 1) := by
            linarith
          exact mul_le_mul_of_nonneg_left this hw0.le
        have hnn : (0 : ℝ) ≤ w * (2 * ((momOrder C (d N) : ℕ) : ℝ)) := by positivity
        have := pow_le_pow_left₀ hnn hle 12
        have h4m : (0 : ℝ) ≤ 4 / m := by positivity
        exact mul_le_mul_of_nonneg_left this h4m
    _ = 4 / m * (2 * C * (Real.log (d N) * w) + 2 * w) ^ 12 := by ring
    _ ≤ 1 / 2 := hG1

set_option linter.unusedVariables false in
/-- K3. `markovConst` is positive. `hc` is unused (`bulkEdge c > 0` holds at every `c`); it
stays so that the signature of the note is the one in the file. -/
theorem markovConst_pos {c ε : ℝ} (hc : 0 < c) (hε : 0 < ε) : 0 < markovConst c ε := by
  have hbe := bulkEdge_pos' c
  have hb1 : (0 : ℝ) < bulkEdge c + ε := by linarith
  have hr1 : (1 : ℝ) < (bulkEdge c + 2 * ε) / (bulkEdge c + ε) := by
    rw [one_lt_div hb1]; linarith
  rw [markovConst]
  exact div_pos (by norm_num) (Real.log_pos hr1)

set_option linter.unusedVariables false in
/-- K4. The Markov bound tends to 0 at `C = markovConst c ε`:
`4 d ((bulkEdge c + ε)/(bulkEdge c + 2ε))^(k_N) ≤ 4/d`. `hc` is unused, as in K3. -/
theorem tendsto_markov_zero {c ε : ℝ} (hc : 0 < c) (hε : 0 < ε) {d : ℕ → ℕ}
    (hdtop : Tendsto d atTop atTop) :
    Tendsto (fun N => 4 * (d N : ℝ) *
        ((bulkEdge c + ε) / (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N)))
      atTop (𝓝 0) := by
  have hbe := bulkEdge_pos' c
  have hb1 : (0 : ℝ) < bulkEdge c + ε := by linarith
  have hb2 : (0 : ℝ) < bulkEdge c + 2 * ε := by linarith
  have hq0 : (0 : ℝ) < (bulkEdge c + ε) / (bulkEdge c + 2 * ε) := div_pos hb1 hb2
  have hq1 : (bulkEdge c + ε) / (bulkEdge c + 2 * ε) ≤ 1 := by
    rw [div_le_one hb2]; linarith
  have hr1 : (1 : ℝ) < (bulkEdge c + 2 * ε) / (bulkEdge c + ε) := by
    rw [one_lt_div hb1]; linarith
  have hlogr : 0 < Real.log ((bulkEdge c + 2 * ε) / (bulkEdge c + ε)) := Real.log_pos hr1
  have hlogq : Real.log ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
      = -Real.log ((bulkEdge c + 2 * ε) / (bulkEdge c + ε)) := by
    rw [← Real.log_inv]
    congr 1
    rw [inv_div]
  have hg : Tendsto (fun N => (4 : ℝ) / (d N : ℝ)) atTop (𝓝 0) :=
    (tendsto_const_div_atTop_nhds_zero_nat (4 : ℝ)).comp hdtop
  refine squeeze_zero' (Eventually.of_forall fun N => by positivity) ?_ hg
  filter_upwards [hdtop.eventually_ge_atTop 1] with N hd1
  have hd0 : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd1
  have hlogd : (0 : ℝ) ≤ Real.log (d N) := Real.log_natCast_nonneg _
  have hk : markovConst c ε * Real.log (d N) ≤ ((momOrder (markovConst c ε) (d N) : ℕ) : ℝ) := by
    rw [momOrder]
    exact Nat.le_ceil _
  have hstep : ((bulkEdge c + ε) / (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
      ≤ (d N : ℝ) ^ (-2 : ℝ) := by
    have h1 : ((bulkEdge c + ε) / (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
        = ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
            ^ (((momOrder (markovConst c ε) (d N) : ℕ) : ℝ)) :=
      (Real.rpow_natCast _ _).symm
    have h2 : ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
          ^ (((momOrder (markovConst c ε) (d N) : ℕ) : ℝ))
        ≤ ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
            ^ (markovConst c ε * Real.log (d N)) :=
      Real.rpow_le_rpow_of_exponent_ge hq0 hq1 hk
    have h3 : ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
          ^ (markovConst c ε * Real.log (d N)) = (d N : ℝ) ^ (-2 : ℝ) := by
      rw [Real.rpow_def_of_pos hq0, Real.rpow_def_of_pos hd0, hlogq, markovConst]
      congr 1
      field_simp
    rw [h1]
    exact h2.trans_eq h3
  have hpow : (d N : ℝ) ^ (-2 : ℝ) = ((d N : ℝ) ^ (2 : ℕ))⁻¹ := by
    rw [show (-2 : ℝ) = -((2 : ℕ) : ℝ) by norm_num, Real.rpow_neg hd0.le, Real.rpow_natCast]
  have hfinal : 4 * (d N : ℝ) * ((d N : ℝ) ^ (2 : ℕ))⁻¹ = 4 / (d N : ℝ) := by
    field_simp
  calc 4 * (d N : ℝ) * ((bulkEdge c + ε) / (bulkEdge c + 2 * ε))
        ^ (momOrder (markovConst c ε) (d N))
      ≤ 4 * (d N : ℝ) * (d N : ℝ) ^ (-2 : ℝ) := by
        have : (0 : ℝ) ≤ 4 * (d N : ℝ) := by positivity
        exact mul_le_mul_of_nonneg_left hstep this
    _ = 4 / (d N : ℝ) := by rw [hpow, hfinal]

end Edge
end StackedSVD
